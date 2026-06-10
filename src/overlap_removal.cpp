// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "overlap_removal.h"

#include <algorithm>
#include <numeric>

#include "collider.h"
#include "cross_section/boolean2/predicates.h"  // IntersectSegments (step 9)
#include "disjoint_sets.h"
#include "impl.h"
#include "manifold/polygon.h"  // for Triangulate
#include "overlap_removal_internal.h"
#include "self_mesh_analysis.h"  // AnalyzeSelfMesh + WindingAt
#include "shared.h"              // for AlphaBudgetEpsilon

// Internal pipeline implementation for Manifold::RemoveSelfIntersections().
// Implements Emmett Lalish's #289 13-step sketch with a per-vert two-sided
// winding classifier (AnalyzeSelfMesh) + pair-symmetric chord enforcement
// (Phases 1, 2, 2.5, 3, 3.5) + pierce-aware cap walker + pre/post-cap pierce
// reducers + pierce/drift gate with sign-flip recovery.

namespace manifold {
namespace overlap_removal {

namespace {
// Numeric defaults shared across the pipeline. Centralized here so
// tuning + audits can change one value, not grep-and-replace.
//
// kPipelineRelTol: relative tolerance for SegmentPiercesTriInterior's
//   boundary-graze gate. 1e-12 lets eps-perturbed verts on tri edges
//   (within FP noise) NOT count as pierces, while keeping real
//   intersections detectable.
// kVolumeFloor: divide-by-zero floor when computing relative drift of
//   a tiny-volume mesh. 1e-12 = near machine precision for doubles.
// kDriftCutoff: relative volume change above which the post-pipeline
//   gate falls back to input. 50% is conservative - geometry-altering
//   pipelines (= the cap walker fills holes / drops slivers) typically
//   change volume by < 5%; > 50% means something pathological happened.
// kPostCapIters: how many iterations of PostCapPierceReducer to run.
//   8 was empirically enough for working fixtures; pierce reductions
//   plateau within 3-4 iters, the extra 4 are headroom.
// kBarycentricFloor: minimum allowed barycentric coordinate when
//   classifying a point as strictly inside a triangle. 1e-12 rejects
//   points "on the boundary" (= one barycentric coord ~= 0) without
//   bumping into FP noise. Used by BuildOnTriVertLists and
//   FindEdgeTriIntersections.
constexpr double kPipelineRelTol = 1e-12;
constexpr double kVolumeFloor = 1e-12;
constexpr double kDriftCutoff = 0.50;
constexpr int kPostCapIters = 8;
constexpr double kBarycentricFloor = 1e-12;

// Iteration caps. All loops below have an early-exit on convergence;
// these caps are tripwires that catch pathological inputs without
// allowing infinite work. DEBUG_ASSERTs check that we hit the
// convergence condition, not the cap, so a cap-hit in DEBUG builds
// surfaces as a logic error.
//
// kPairSymPhase3MaxIter / kPairSymPhase35MaxIter: monotonic
//   drop / re-key passes in pair-sym Phase 3 / 3.5. Strict upper
//   bound is the polygon count (each pass changes >=1 polygon or
//   exits). For typical mesh sizes (tens of thousands of polygons)
//   we converge in < 16 iters; cap at 64 catches a 4x slowdown
//   pathology without allowing 30k-poly meshes to grind for minutes.
// kCascadeDropMaxPass: forward-cascade drop in TriangulateAndEmit's
//   classifier post-pass. Same monotonicity argument; 16 covers all
//   working fixtures.
// kCycleWalkerGuard / kEarClipGuard: per-polygon walk guards. Max
//   work is the polygon size; 4096 covers any realistic single-tri
//   subdivision (typical: tens of verts; pathological: hundreds).
//
// (kMergeVertsMaxIter / kPierceReducerMaxIter / kDropExcessOuterMax
// / kTrimOrphansMaxRounds live in overlap_removal_internal.h as
// header constants because they're used as parameter defaults in
// the corresponding function declarations.)
constexpr int kPairSymPhase3MaxIter = 64;
constexpr int kPairSymPhase35MaxIter = 64;
constexpr int kCascadeDropMaxPass = 16;
constexpr int kCycleWalkerGuard = 4096;
constexpr int kEarClipGuard = 4096;
}  // namespace

SortedBVH BuildSortedBVH(VecView<const Box> leafBoxes) {
  SortedBVH out;
  const size_t n = leafBoxes.size();
  if (n == 0) return out;
  Box bbox;
  for (const auto& b : leafBoxes) bbox = bbox.Union(b);
  std::vector<uint32_t> rawMorton(n);
  for (size_t i = 0; i < n; ++i)
    rawMorton[i] = Collider::MortonCode(leafBoxes[i].Center(), bbox);
  out.perm.resize(n);
  std::iota(out.perm.begin(), out.perm.end(), size_t{0});
  std::stable_sort(out.perm.begin(), out.perm.end(), [&](size_t a, size_t b) {
    return rawMorton[a] < rawMorton[b];
  });
  out.boxes.resize(n);
  out.morton.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.boxes[i] = leafBoxes[out.perm[i]];
    out.morton[i] = rawMorton[out.perm[i]];
  }
  out.collider =
      Collider(VecView<const Box>(out.boxes.data(), out.boxes.size()),
               VecView<const uint32_t>(out.morton.data(), out.morton.size()));
  return out;
}

double InferEps(const Manifold& m) {
  return AlphaBudgetEpsilon(m.BoundingBox().Scale(), 1000);
}

Manifold::Impl ImplFromManifold(const Manifold& m) {
  return Manifold::Impl(m.GetMeshGL64());
}

MergeVertsResult MergeVertsEps(const Manifold& in, double eps, int maxIter) {
  if (in.IsEmpty()) return {in, 0};

  // Pull positions and triangles out via MeshGL64. Layout:
  // vertProperties = [x0,y0,z0, x1,y1,z1, ...], numProp >= 3.
  MeshGL64 mesh = in.GetMeshGL64();
  const size_t n = mesh.NumVert();
  std::vector<vec3> verts(n);
  for (size_t i = 0; i < n; ++i) {
    verts[i] = vec3(mesh.vertProperties[mesh.numProp * i + 0],
                    mesh.vertProperties[mesh.numProp * i + 1],
                    mesh.vertProperties[mesh.numProp * i + 2]);
  }

  // Per-pass: build eps/2-padded boxes, run Collider self-collision,
  // narrow-phase dist^2 < eps^2, unite via DisjointSets. Move each
  // cluster to its centroid and repeat until stable.
  int iter = 0;
  bool converged = false;
  std::vector<int> componentLabel(n);
  for (; iter < maxIter; ++iter) {
    std::vector<Box> boxes(n);
    const double halfEps = 0.5 * eps;
    for (size_t i = 0; i < n; ++i) {
      const vec3 lo(verts[i].x - halfEps, verts[i].y - halfEps,
                    verts[i].z - halfEps);
      const vec3 hi(verts[i].x + halfEps, verts[i].y + halfEps,
                    verts[i].z + halfEps);
      boxes[i] = Box(lo, hi);
    }
    SortedBVH bvh =
        BuildSortedBVH(VecView<const Box>(boxes.data(), boxes.size()));
    DisjointSets uf(static_cast<uint32_t>(n));
    int unions = 0;
    const double eps2 = eps * eps;
    auto checkPair = [&](size_t qi, size_t li) {
      if (qi >= li) return;
      const size_t va = bvh.perm[qi];
      const size_t vb = bvh.perm[li];
      const vec3 d = verts[va] - verts[vb];
      const double d2 = la::dot(d, d);
      if (d2 > eps2) return;
      uint32_t before = uf.find(static_cast<uint32_t>(va));
      uf.unite(static_cast<uint32_t>(va), static_cast<uint32_t>(vb));
      uint32_t after = uf.find(static_cast<uint32_t>(va));
      if (before != after) ++unions;
    };
    auto recorder = MakeSimpleRecorder(checkPair);
    auto qf = [&](int i) { return bvh.boxes[i]; };
    bvh.collider.Collisions<false>(recorder, qf, static_cast<int>(n),
                                   /*parallel=*/false);

    // Update positions to per-cluster centroid.
    int nComp = uf.connectedComponents(componentLabel);
    std::vector<vec3> sumByComp(nComp, vec3(0, 0, 0));
    std::vector<int> countByComp(nComp, 0);
    for (size_t i = 0; i < n; ++i) {
      sumByComp[componentLabel[i]] += verts[i];
      ++countByComp[componentLabel[i]];
    }
    bool moved = false;
    for (size_t i = 0; i < n; ++i) {
      const vec3 newPos = sumByComp[componentLabel[i]] /
                          static_cast<double>(countByComp[componentLabel[i]]);
      if (newPos.x != verts[i].x || newPos.y != verts[i].y ||
          newPos.z != verts[i].z) {
        verts[i] = newPos;
        moved = true;
      }
    }
    if (unions == 0 && !moved) {
      converged = true;
      break;
    }
  }
  DEBUG_ASSERT(converged, logicErr,
               "MergeVertsEps: hit kMergeVertsMaxIter without converging");

  // Apply merges via MeshGL64 hints. Same path manifold's sort.cpp
  // uses for eps-merging during construction, so result is consistent
  // with other Manifold-producing paths. Also overwrite each merged
  // vert's position with its cluster centroid.
  for (size_t i = 0; i < n; ++i) {
    const vec3& p = verts[i];
    mesh.vertProperties[mesh.numProp * i + 0] = p.x;
    mesh.vertProperties[mesh.numProp * i + 1] = p.y;
    mesh.vertProperties[mesh.numProp * i + 2] = p.z;
  }
  std::vector<int> compRep(componentLabel.size(), -1);
  for (size_t i = 0; i < n; ++i) {
    const int comp = componentLabel[i];
    if (compRep[comp] == -1) compRep[comp] = static_cast<int>(i);
  }
  mesh.mergeFromVert.clear();
  mesh.mergeToVert.clear();
  int mergedCount = 0;
  for (size_t i = 0; i < n; ++i) {
    const int rep = compRep[componentLabel[i]];
    if (rep != static_cast<int>(i)) {
      mesh.mergeFromVert.push_back(i);
      mesh.mergeToVert.push_back(rep);
      ++mergedCount;
    }
  }
  // Avoid round-trip when nothing was merged: GetMeshGL64 -> Manifold
  // is lossy for Subtract-derived inputs (back-side / run-transform
  // info doesn't fully survive, observed as sign-flipped volume on
  // Cray). When no merges to apply, the input is already correct.
  if (mergedCount == 0) return {in, 0};
  return {Manifold(mesh), mergedCount};
}

std::vector<Edge> EnumerateEdges(const Manifold::Impl& impl) {
  // Each manifold edge is two halfedges; the canonical (forward) one
  // is the side where startVert < endVert. Walk all halfedges and
  // collect those.
  std::vector<Edge> out;
  out.reserve(impl.halfedge_.size() / 2);
  for (size_t h = 0; h < impl.halfedge_.size(); ++h) {
    const Halfedge he = impl.halfedge_.Get(h);
    if (!he.IsForward()) continue;
    out.push_back(
        {he.startVert, he.endVert, static_cast<int>(h), he.pairedHalfedge});
  }
  return out;
}

std::vector<EdgeVertList> BuildOnEdgeVertLists(const Manifold::Impl& impl,
                                               const std::vector<Edge>& edges,
                                               double eps) {
  using la::dot;
  const size_t nE = edges.size();
  const size_t nV = impl.vertPos_.size();
  const double eps2 = eps * eps;
  std::vector<EdgeVertList> out(nE);
  if (nE == 0 || nV == 0) return out;

  // vert->neighbor adjacency for the thin-tri-apex skip in
  // GenerateChordEdges.
  std::vector<std::set<int>> adj(nV);
  for (size_t i = 0; i < impl.halfedge_.size(); ++i) {
    const int s = impl.halfedge_.Start(i);
    const int e = impl.halfedge_.End(i);
    adj[s].insert(e);
    adj[e].insert(s);
  }

  // Per-edge eps-padded AABB (BVH leaves) and per-vert eps-padded AABB
  // (queries).
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    edgeBoxes[i] = b;
  }
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }

  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(edgeBoxes.data(), edgeBoxes.size()));

  std::vector<std::vector<std::pair<double, int>>> hitsByEdge(nE);
  auto onCollision = [&](size_t vertIdxQ, size_t edgeIdxL) {
    const size_t vertIdx = vertIdxQ;
    const size_t edgeIdx = bvh.perm[edgeIdxL];
    const auto& edge = edges[edgeIdx];
    if (static_cast<int>(vertIdx) == edge.v0 ||
        static_cast<int>(vertIdx) == edge.v1)
      return;
    if (adj[vertIdx].count(edge.v0) && adj[vertIdx].count(edge.v1)) return;
    const vec3 a = impl.vertPos_[edge.v0];
    const vec3 b = impl.vertPos_[edge.v1];
    const vec3 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) return;
    const vec3 p = impl.vertPos_[vertIdx];
    const vec3 ap = p - a;
    const double t = dot(ap, ab) / abLen2;
    if (t <= 0.0 || t >= 1.0) return;
    const vec3 closest = a + ab * t;
    const vec3 d = p - closest;
    if (dot(d, d) <= eps2)
      hitsByEdge[edgeIdx].emplace_back(t, static_cast<int>(vertIdx));
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  bvh.collider.Collisions<false>(recorder, qf, static_cast<int>(nV),
                                 /*parallel=*/false);

  for (size_t e = 0; e < nE; ++e) {
    auto& hits = hitsByEdge[e];
    std::sort(hits.begin(), hits.end());
    out[e].verts.reserve(hits.size());
    out[e].ts.reserve(hits.size());
    for (const auto& [t, v] : hits) {
      out[e].verts.push_back(v);
      out[e].ts.push_back(t);
    }
  }
  return out;
}

std::vector<TriVertList> BuildOnTriVertLists(const Manifold::Impl& impl,
                                             double eps) {
  using la::cross;
  using la::dot;
  const size_t nT = impl.NumTri();
  const size_t nV = impl.vertPos_.size();
  std::vector<TriVertList> out(nT);
  if (nT == 0 || nV == 0) return out;

  std::vector<Box> triBoxes(nT);
  for (size_t t = 0; t < nT; ++t) {
    const int v0 = impl.halfedge_.Start(3 * t + 0);
    const int v1 = impl.halfedge_.Start(3 * t + 1);
    const int v2 = impl.halfedge_.Start(3 * t + 2);
    Box b(impl.vertPos_[v0], impl.vertPos_[v1]);
    b.Union(impl.vertPos_[v2]);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    triBoxes[t] = b;
  }
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));

  auto onCollision = [&](size_t vertIdx, size_t triIdxL) {
    const size_t triIdx = bvh.perm[triIdxL];
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    if (static_cast<int>(vertIdx) == t0 || static_cast<int>(vertIdx) == t1 ||
        static_cast<int>(vertIdx) == t2)
      return;
    const vec3 a = impl.vertPos_[t0];
    const vec3 b = impl.vertPos_[t1];
    const vec3 c = impl.vertPos_[t2];
    const vec3 e1 = b - a;
    const vec3 e2 = c - a;
    const vec3 n = cross(e1, e2);
    const double nMag2 = dot(n, n);
    if (nMag2 == 0) return;
    const double nMag = std::sqrt(nMag2);
    const vec3 p = impl.vertPos_[vertIdx];
    const double dist = dot(p - a, n) / nMag;
    if (std::fabs(dist) > eps) return;
    const vec3 pProj = p - dist * (n / nMag);
    const vec3 v0 = b - a, v1 = c - a, v2 = pProj - a;
    const double d00 = dot(v0, v0);
    const double d01 = dot(v0, v1);
    const double d11 = dot(v1, v1);
    const double d20 = dot(v2, v0);
    const double d21 = dot(v2, v1);
    const double denom = d00 * d11 - d01 * d01;
    if (denom == 0) return;
    const double bv = (d11 * d20 - d01 * d21) / denom;
    const double bw = (d00 * d21 - d01 * d20) / denom;
    const double bu = 1.0 - bv - bw;
    if (bu <= kBarycentricFloor || bv <= kBarycentricFloor ||
        bw <= kBarycentricFloor)
      return;
    out[triIdx].verts.push_back(static_cast<int>(vertIdx));
    out[triIdx].bary.push_back(vec3(bu, bv, bw));
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  bvh.collider.Collisions<false>(recorder, qf, static_cast<int>(nV),
                                 /*parallel=*/false);
  return out;
}

std::vector<EdgeTriIntersection> FindEdgeTriIntersections(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using la::cross;
  using la::dot;
  std::vector<EdgeTriIntersection> out;
  const size_t nE = edges.size();
  const size_t nT = impl.NumTri();
  if (nE == 0 || nT == 0) return out;
  const double eps2 = eps * eps;

  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    edgeBoxes[i] = b;
  }
  std::vector<Box> triBoxes(nT);
  for (size_t t = 0; t < nT; ++t) {
    const int v0 = impl.halfedge_.Start(3 * t + 0);
    const int v1 = impl.halfedge_.Start(3 * t + 1);
    const int v2 = impl.halfedge_.Start(3 * t + 2);
    Box b(impl.vertPos_[v0], impl.vertPos_[v1]);
    b.Union(impl.vertPos_[v2]);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    triBoxes[t] = b;
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(edgeBoxes.data(), edgeBoxes.size()));

  auto onCollision = [&](size_t triIdxQ, size_t edgeIdxL) {
    const size_t triIdx = triIdxQ;
    const size_t edgeIdx = bvh.perm[edgeIdxL];
    const Edge& edge = edges[edgeIdx];
    if (static_cast<int>(triIdx) == edge.halfedgeForward / 3 ||
        (edge.halfedgePaired >= 0 &&
         static_cast<int>(triIdx) == edge.halfedgePaired / 3))
      return;
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    const vec3 a = impl.vertPos_[edge.v0];
    const vec3 b = impl.vertPos_[edge.v1];
    const vec3 p0 = impl.vertPos_[t0];
    const vec3 p1 = impl.vertPos_[t1];
    const vec3 p2 = impl.vertPos_[t2];
    const vec3 e1 = p1 - p0;
    const vec3 e2 = p2 - p0;
    const vec3 n = cross(e1, e2);
    const vec3 d = b - a;
    const double denom = dot(n, d);
    const double nMag2 = dot(n, n);
    if (nMag2 == 0) return;
    // Scale-relative parallel-edge reject (denom has units length^3 =
    // |n| * |d|), matching the rest of the pipeline's relative tolerances.
    if (std::fabs(denom) <=
        kPipelineRelTol * std::sqrt(nMag2) * std::sqrt(dot(d, d)))
      return;
    const double s = -dot(n, a - p0) / denom;
    if (s <= 0.0 || s >= 1.0) return;
    const vec3 pos = a + s * d;
    const vec3 v0 = e1, v1 = e2, v2 = pos - p0;
    const double d00 = dot(v0, v0);
    const double d01 = dot(v0, v1);
    const double d11 = dot(v1, v1);
    const double d20 = dot(v2, v0);
    const double d21 = dot(v2, v1);
    const double bDenom = d00 * d11 - d01 * d01;
    if (bDenom == 0) return;
    const double bv = (d11 * d20 - d01 * d21) / bDenom;
    const double bw = (d00 * d21 - d01 * d20) / bDenom;
    const double bu = 1.0 - bv - bw;
    if (bu <= kBarycentricFloor || bv <= kBarycentricFloor ||
        bw <= kBarycentricFloor)
      return;

    int snapTo = -1;
    auto trySnap = [&](int v) {
      if (snapTo >= 0) return;
      const vec3 dd = pos - impl.vertPos_[v];
      if (dot(dd, dd) <= eps2) snapTo = v;
    };
    trySnap(edge.v0);
    trySnap(edge.v1);
    trySnap(t0);
    trySnap(t1);
    trySnap(t2);
    for (int v : onEdgeLists[edgeIdx].verts) trySnap(v);
    for (int v : onTriLists[triIdx].verts) trySnap(v);

    out.push_back({static_cast<int>(edgeIdx), static_cast<int>(triIdx), pos, s,
                   vec3(bu, bv, bw), snapTo});
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return triBoxes[i]; };
  bvh.collider.Collisions<false>(recorder, qf, static_cast<int>(nT),
                                 /*parallel=*/false);
  return out;
}

ChordEdges GenerateChordEdges(const Manifold::Impl& impl,
                              const std::vector<Edge>& edges,
                              const std::vector<EdgeTriIntersection>& etIsects,
                              double eps) {
  using la::dot;
  ChordEdges r;
  if (etIsects.empty()) return r;
  const double eps2 = eps * eps;
  const int baseId = static_cast<int>(impl.NumVert());

  // Resolve each etIsect to a vert id, deduping new positions
  // against each other within eps.
  r.resolvedIds.assign(etIsects.size(), -1);
  std::vector<int>& resolvedId = r.resolvedIds;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    if (x.snapTo >= 0) {
      resolvedId[i] = x.snapTo;
      continue;
    }
    int found = -1;
    for (size_t j = 0; j < r.newVertPositions.size(); ++j) {
      const vec3 d = x.position - r.newVertPositions[j];
      if (dot(d, d) <= eps2) {
        found = baseId + static_cast<int>(j);
        break;
      }
    }
    if (found >= 0) {
      resolvedId[i] = found;
    } else {
      resolvedId[i] = baseId + static_cast<int>(r.newVertPositions.size());
      r.newVertPositions.push_back(x.position);
    }
  }

  // Group resolved endpoints by tri-tri pair.
  std::map<std::pair<int, int>, std::set<int>> pairEndpoints;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    const Edge& e = edges[x.edgeIdx];
    const int triA1 = e.halfedgeForward / 3;
    const int triA2 = e.halfedgePaired >= 0 ? e.halfedgePaired / 3 : -1;
    auto add = [&](int t1, int t2) {
      auto key = (t1 < t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
      pairEndpoints[key].insert(resolvedId[i]);
    };
    add(triA1, x.triIdx);
    if (triA2 >= 0) add(triA2, x.triIdx);
  }

  // Emit chords (= 2-endpoint pairs) or interior vert records (= 1-
  // endpoint pairs); drop 0 / >= 3.
  for (auto& [key, eps_set] : pairEndpoints) {
    if (eps_set.size() == 2) {
      auto it = eps_set.begin();
      const int a = *it++;
      const int b = *it;
      const int v0 = std::min(a, b);
      const int v1 = std::max(a, b);
      r.newEdges.push_back({v0, v1, key.first, key.second});
    } else if (eps_set.size() == 1) {
      // Single shared endpoint (edge-tip touch, no through-pierce): no
      // chord edge is emitted and nothing further is recorded.
    } else {
      ++r.droppedTriTriPairsWithBadEndpointCount;
    }
  }
  return r;
}

vec3 GetPos3(int id, int baseId, const Manifold::Impl& impl,
             const std::vector<vec3>& newVertPositions) {
  if (id < baseId) {
    DEBUG_ASSERT(id >= 0 && id < static_cast<int>(impl.vertPos_.size()),
                 logicErr, "GetPos3: original vert id out of range");
    return impl.vertPos_[id];
  }
  const int j = id - baseId;
  DEBUG_ASSERT(j >= 0 && j < static_cast<int>(newVertPositions.size()),
               logicErr, "GetPos3: chord vert id beyond newVertPositions");
  if (j < 0 || j >= static_cast<int>(newVertPositions.size())) {
    return vec3(0.0, 0.0, 0.0);
  }
  return newVertPositions[j];
}

std::vector<NewEdgeWithExtras> AddInteriorVertsToNewEdges(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<PiercedNewEdge>& newEdges,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using la::dot;
  std::vector<NewEdgeWithExtras> out;
  out.reserve(newEdges.size());
  const int baseId = static_cast<int>(impl.NumVert());
  const double eps2 = eps * eps;

  auto getPos = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  for (const auto& edge : newEdges) {
    NewEdgeWithExtras nwe{edge, {}, {}};
    const vec3 a = getPos(edge.v0);
    const vec3 b = getPos(edge.v1);
    const vec3 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) {
      out.push_back(std::move(nwe));
      continue;
    }
    std::set<int> seen;
    auto check = [&](int v) {
      if (!seen.insert(v).second) return;
      if (v == edge.v0 || v == edge.v1) return;
      // Candidates come from onTriLists, which hold only original-mesh
      // verts; a chord vert id (>= baseId) here would read out of
      // bounds below (step 9 threads those separately).
      DEBUG_ASSERT(v >= 0 && v < baseId, logicErr,
                   "AddInteriorVertsToNewEdges: candidate must be an "
                   "original-mesh vert");
      const vec3 p = impl.vertPos_[v];
      const vec3 ap = p - a;
      const double t = dot(ap, ab) / abLen2;
      if (t <= 0.0 || t >= 1.0) return;
      const vec3 closest = a + ab * t;
      const vec3 d = p - closest;
      if (dot(d, d) > eps2) return;
      nwe.extraVerts.push_back(v);
      nwe.extraTs.push_back(t);
    };
    for (int v : onTriLists[edge.triA].verts) check(v);
    for (int v : onTriLists[edge.triB].verts) check(v);
    std::vector<size_t> perm(nwe.extraVerts.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](size_t i, size_t j) {
      return nwe.extraTs[i] < nwe.extraTs[j];
    });
    std::vector<int> sortedV(nwe.extraVerts.size());
    std::vector<double> sortedT(nwe.extraTs.size());
    for (size_t i = 0; i < perm.size(); ++i) {
      sortedV[i] = nwe.extraVerts[perm[i]];
      sortedT[i] = nwe.extraTs[perm[i]];
    }
    nwe.extraVerts = std::move(sortedV);
    nwe.extraTs = std::move(sortedT);
    out.push_back(std::move(nwe));
  }
  return out;
}

std::vector<std::vector<int>> GroupChordsByFace(
    const std::vector<PiercedNewEdge>& newEdges, int numTri) {
  std::vector<std::vector<int>> byFace(numTri);
  for (size_t i = 0; i < newEdges.size(); ++i) {
    const PiercedNewEdge& e = newEdges[i];
    const int chord = static_cast<int>(i);
    if (e.triA >= 0 && e.triA < numTri) byFace[e.triA].push_back(chord);
    if (e.triB >= 0 && e.triB < numTri && e.triB != e.triA) {
      byFace[e.triB].push_back(chord);
    }
  }
  return byFace;
}

std::vector<OnChordContact> FindOnChordEndpointContacts(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& chordsByFace, double tolerance,
    double eps) {
  using la::dot;
  std::vector<OnChordContact> out;
  const int baseId = static_cast<int>(impl.NumVert());
  // New-to-old snap: prior drift plus current-op error, matching
  // boolean2's newToOldThresh (tolerance + eps), not bare eps.
  const double snap = tolerance + eps;
  const double snap2 = snap * snap;
  // A chord pair shares up to two faces; record each (chord, vert)
  // contact once.
  std::set<std::pair<int, int>> seen;
  for (const std::vector<int>& faceChords : chordsByFace) {
    for (int ci : faceChords) {
      const PiercedNewEdge& c = chords[ci].edge;
      for (int di : faceChords) {
        if (di == ci) continue;
        const PiercedNewEdge& d = chords[di].edge;
        const vec3 a = GetPos3(d.v0, baseId, impl, newVertPositions);
        const vec3 b = GetPos3(d.v1, baseId, impl, newVertPositions);
        const vec3 ab = b - a;
        const double abLen2 = dot(ab, ab);
        if (abLen2 == 0) continue;
        // Endpoint-proximity zone in t-space: a contact within snap of
        // d's own endpoints is an endpoint-near-endpoint case, not an
        // interior vert (it would re-create the near-line sliver the
        // guard exists to block). snap/len >= 0.5 excludes the whole
        // chord (shorter than 2*snap).
        const double tGuard = snap / std::sqrt(abLen2);
        for (const int e : {c.v0, c.v1}) {
          if (e == d.v0 || e == d.v1) continue;
          const vec3 p = GetPos3(e, baseId, impl, newVertPositions);
          const double t = dot(p - a, ab) / abLen2;
          if (t <= tGuard || t >= 1.0 - tGuard) continue;
          const vec3 closest = a + ab * t;
          const vec3 dv = p - closest;
          if (dot(dv, dv) > snap2) continue;
          if (!seen.insert({di, e}).second) continue;
          out.push_back({di, e, t});
        }
      }
    }
  }
  // Deterministic output order: ascending chord, then ascending t.
  std::sort(out.begin(), out.end(),
            [](const OnChordContact& x, const OnChordContact& y) {
              if (x.chord != y.chord) return x.chord < y.chord;
              return x.t < y.t;
            });
  return out;
}

std::vector<ChordChordCrossing> FindChordChordCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& chordsByFace,
    VecView<const vec3> faceNormals, double eps) {
  using la::cross;
  using la::dot;
  std::vector<ChordChordCrossing> out;
  const int baseId = static_cast<int>(impl.NumVert());
  // A chord pair shares up to two faces; its crossing is recorded once
  // (lowest face wins by iteration order).
  std::set<std::pair<int, int>> seenPairs;
  for (size_t face = 0; face < chordsByFace.size(); ++face) {
    const std::vector<int>& faceChords = chordsByFace[face];
    if (faceChords.size() < 2) continue;
    DEBUG_ASSERT(face < faceNormals.size(), logicErr,
                 "FindChordChordCrossings: face normal missing");
    if (face >= faceNormals.size()) continue;
    const vec3 nRaw = faceNormals[face];
    const double nLen2 = dot(nRaw, nRaw);
    if (nLen2 == 0) continue;
    const vec3 n = nRaw / std::sqrt(nLen2);
    // Deterministic orthonormal in-plane basis - a true isometry, so
    // the kernel's 2D eps equals the pipeline's 3D eps (the axis-drop
    // projection used elsewhere is not an isometry and would contract
    // distances).
    const double ax = std::fabs(n.x);
    const double ay = std::fabs(n.y);
    const double az = std::fabs(n.z);
    const vec3 ref = (ax <= ay && ax <= az) ? vec3(1.0, 0.0, 0.0)
                     : (ay <= az)           ? vec3(0.0, 1.0, 0.0)
                                            : vec3(0.0, 0.0, 1.0);
    vec3 u = cross(n, ref);
    u = u / std::sqrt(dot(u, u));
    const vec3 v = cross(n, u);
    for (size_t i = 0; i < faceChords.size(); ++i) {
      for (size_t j = i + 1; j < faceChords.size(); ++j) {
        const int ci = faceChords[i];
        const int cj = faceChords[j];
        if (!seenPairs.insert({std::min(ci, cj), std::max(ci, cj)}).second) {
          continue;
        }
        const PiercedNewEdge& ea = chords[ci].edge;
        const PiercedNewEdge& eb = chords[cj].edge;
        // Re-project endpoints onto the face plane (through ea.v0)
        // before the 2D mapping, so out-of-plane drift from prior
        // merges is zero by construction.
        const vec3 planePt = GetPos3(ea.v0, baseId, impl, newVertPositions);
        auto to2d = [&](int id) -> vec2 {
          const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
          const vec3 inPlane = p - dot(p - planePt, n) * n;
          return vec2(dot(inPlane - planePt, u), dot(inPlane - planePt, v));
        };
        const vec2 a0 = to2d(ea.v0);
        const vec2 a1 = to2d(ea.v1);
        const vec2 b0 = to2d(eb.v0);
        const vec2 b1 = to2d(eb.v1);
        // stableEdgeId: the chord's index - deterministic (newEdges
        // order comes from GenerateChordEdges' sorted-pair map walk),
        // NOT broad-phase pair order.
        const boolean2::GraphSegment2D segA{a0, a1, ci};
        const boolean2::GraphSegment2D segB{b0, b1, cj};
        vec2 p2;
        if (!boolean2::IntersectSegments(segA, segB, eps, &p2)) continue;
        const vec2 da = a1 - a0;
        const vec2 db = b1 - b0;
        const double tA = dot(p2 - a0, da) / dot(da, da);
        const double tB = dot(p2 - b0, db) / dot(db, db);
        const vec3 pos = planePt + p2.x * u + p2.y * v;
        out.push_back({pos, ci, cj, tA, tB, static_cast<int>(face)});
      }
    }
  }
  return out;
}

Step9Threading ResolveAndThreadClusters(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions,
    const std::vector<ChordCrossing>& clusters,
    const std::vector<OnChordContact>& contacts, double tolerance, double eps) {
  using la::dot;
  const int baseId = static_cast<int>(impl.NumVert());
  const double snap = tolerance + eps;
  const double snap2 = snap * snap;
  // Canonical-id resolution: resolve-then-allocate, per cluster,
  // symmetric across ALL incident chords. Candidates are every
  // existing vert the crossing could BE - chord endpoints,
  // already-threaded on-chord verts, and the pass-0 contacts (not yet
  // threaded at this point) - within tolerance + eps. Nearest wins;
  // ties take the smallest id. Only when no candidate exists is a
  // fresh vert allocated. This is what prevents a crossing threading
  // as an endpoint id on one chord and a fresh id on another (the
  // split-identity bug).
  std::vector<ChordCrossing> crossings;
  crossings.reserve(clusters.size());
  for (const ChordCrossing& cl : clusters) {
    int best = -1;
    double bestD2 = 0.0;
    auto consider = [&](int id) {
      const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
      const vec3 dv = cl.pos - p;
      const double d2 = dot(dv, dv);
      if (d2 > snap2) return;
      if (best < 0 || d2 < bestD2 || (d2 == bestD2 && id < best)) {
        best = id;
        bestD2 = d2;
      }
    };
    for (const int side : cl.chords) {
      const PiercedNewEdge& e = chords[side].edge;
      consider(e.v0);
      consider(e.v1);
      for (const int v : chords[side].extraVerts) consider(v);
    }
    for (const OnChordContact& c : contacts) {
      if (std::find(cl.chords.begin(), cl.chords.end(), c.chord) !=
          cl.chords.end()) {
        consider(c.vertId);
      }
    }
    int id = best;
    if (id < 0) {
      id = baseId + static_cast<int>(newVertPositions.size());
      newVertPositions.push_back(cl.pos);
    }
    ChordCrossing rec = cl;
    rec.id = id;
    crossings.push_back(std::move(rec));
  }
  // Threading: per chord, the unified record list (existing extras +
  // pass-0 contacts + resolved crossings), with every t RECOMPUTED
  // from the resolved position (a snap can move the vertex by up to
  // tolerance + eps - enough to reorder a stale t-sort), the pass-0
  // endpoint-zone guard re-applied, id-dedup over the unified list,
  // then t-sort with an eps/len dedup backstop.
  std::vector<std::vector<int>> pending(chords.size());
  for (const OnChordContact& c : contacts) {
    pending[c.chord].push_back(c.vertId);
  }
  for (const ChordCrossing& cc : crossings) {
    for (const int ch : cc.chords) pending[ch].push_back(cc.id);
  }
  for (size_t ci = 0; ci < chords.size(); ++ci) {
    if (pending[ci].empty()) continue;
    NewEdgeWithExtras& nwe = chords[ci];
    std::vector<int> ids = nwe.extraVerts;
    ids.insert(ids.end(), pending[ci].begin(), pending[ci].end());
    const vec3 a = GetPos3(nwe.edge.v0, baseId, impl, newVertPositions);
    const vec3 b = GetPos3(nwe.edge.v1, baseId, impl, newVertPositions);
    const vec3 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) continue;
    const double len = std::sqrt(abLen2);
    const double tGuard = snap / len;
    const double tDedup = eps / len;
    // Pre-existing step-8 extras were admitted under step 8's weaker
    // t-guard; the step-9 endpoint-zone guard applies only to NEWLY
    // added records (contacts + crossings), so a rebuild does not
    // silently drop legitimate near-endpoint on-tri verts.
    const std::set<int> preexisting(nwe.extraVerts.begin(),
                                    nwe.extraVerts.end());
    std::vector<std::pair<double, int>> recs;
    recs.reserve(ids.size());
    std::set<int> seenIds;
    for (const int id : ids) {
      if (id == nwe.edge.v0 || id == nwe.edge.v1) continue;
      if (!seenIds.insert(id).second) continue;
      const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
      const double t = dot(p - a, ab) / abLen2;
      const bool isNew = preexisting.count(id) == 0;
      if (isNew && (t <= tGuard || t >= 1.0 - tGuard)) continue;
      recs.push_back({t, id});
    }
    std::sort(recs.begin(), recs.end());
    std::vector<int> outV;
    std::vector<double> outT;
    outV.reserve(recs.size());
    outT.reserve(recs.size());
    for (const auto& [t, id] : recs) {
      if (!outT.empty() && t - outT.back() <= tDedup) continue;
      outT.push_back(t);
      outV.push_back(id);
    }
    nwe.extraVerts = std::move(outV);
    nwe.extraTs = std::move(outT);
  }
  return {std::move(chords), std::move(newVertPositions), std::move(crossings)};
}

Step9Threading ResolveAndThreadCrossings(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<OnChordContact>& contacts, double tolerance, double eps) {
  // Singleton-cluster delegation: each raw crossing is its own
  // cluster (increment (ii) form; MergeAndPropagateCrossings supplies
  // real clusters in the full pipeline).
  std::vector<ChordCrossing> clusters;
  clusters.reserve(raw.size());
  for (const ChordChordCrossing& rc : raw) {
    clusters.push_back({rc.pos, -1, {rc.chordA, rc.chordB}, {rc.tA, rc.tB}});
  }
  return ResolveAndThreadClusters(impl, std::move(chords),
                                  std::move(newVertPositions), clusters,
                                  contacts, tolerance, eps);
}

std::vector<ChordCrossing> MergeAndPropagateCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<std::vector<int>>& chordsByFace,
    VecView<const vec3> faceNormals, double tolerance, double eps) {
  using la::dot;
  std::vector<ChordCrossing> out;
  if (raw.empty()) return out;
  const int baseId = static_cast<int>(impl.NumVert());
  const double snap = tolerance + eps;
  // New-to-new merge radius, matching boolean2's
  // kIntersectionMergeEpsFactor (covers shallow crossings to ~6 deg).
  const double mergeR = 10.0 * eps;
  const double mergeR2 = mergeR * mergeR;

  // Incident faces of a raw crossing: the hosting face plus both
  // chords' face pairs.
  const uint32_t n = static_cast<uint32_t>(raw.size());
  std::vector<std::set<int>> rcFaces(n);
  for (uint32_t i = 0; i < n; ++i) {
    rcFaces[i].insert(raw[i].face);
    for (const int ch : {raw[i].chordA, raw[i].chordB}) {
      rcFaces[i].insert(chords[ch].edge.triA);
      rcFaces[i].insert(chords[ch].edge.triB);
    }
  }

  // Face-gated union-find in sorted pair order: unite when the two
  // crossings share an incident face AND lie within 10 * eps. The
  // FACE gate (not a chord gate) is what unites a 4-chord concurrence
  // whose two crossings share no chord.
  DisjointSets uf(n);
  for (uint32_t i = 0; i < n; ++i) {
    for (uint32_t j = i + 1; j < n; ++j) {
      const vec3 d = raw[i].pos - raw[j].pos;
      if (dot(d, d) > mergeR2) continue;
      bool shareFace = false;
      for (const int f : rcFaces[i]) {
        if (rcFaces[j].count(f) > 0) {
          shareFace = true;
          break;
        }
      }
      if (!shareFace) continue;
      uf.unite(i, j);
    }
  }

  // Clusters keyed by their smallest member index, members ascending,
  // so output order and centroid summation are deterministic.
  std::map<uint32_t, uint32_t> minOfRoot;
  for (uint32_t i = 0; i < n; ++i) {
    const uint32_t r = uf.find(i);
    auto it = minOfRoot.find(r);
    if (it == minOfRoot.end()) minOfRoot[r] = i;
  }
  std::map<uint32_t, std::vector<uint32_t>> byMin;
  for (uint32_t i = 0; i < n; ++i) {
    byMin[minOfRoot[uf.find(i)]].push_back(i);
  }

  for (const auto& [minIdx, members] : byMin) {
    // Centroid in ascending member order, then re-projected onto the
    // hosting face plane (lowest hosting face) so the merged position
    // does not drift out of plane.
    vec3 centroid(0.0, 0.0, 0.0);
    int host = raw[members[0]].face;
    for (const uint32_t m : members) {
      centroid = centroid + raw[m].pos;
      host = std::min(host, raw[m].face);
    }
    centroid = centroid / static_cast<double>(members.size());
    if (host >= 0 && static_cast<size_t>(host) < faceNormals.size()) {
      const vec3 nRaw = faceNormals[host];
      const double nLen2 = dot(nRaw, nRaw);
      if (nLen2 > 0) {
        const vec3 nrm = nRaw / std::sqrt(nLen2);
        // The plane point must lie on the HOST face's plane: take it
        // from the first member hosted there. raw need not arrive in
        // face order, so members[0]'s hosting face can differ from
        // host - pairing host's normal with another face's point would
        // project the centroid onto neither plane.
        uint32_t hostMember = members[0];
        for (const uint32_t m : members) {
          if (raw[m].face == host) {
            hostMember = m;
            break;
          }
        }
        const vec3 planePt = GetPos3(chords[raw[hostMember].chordA].edge.v0,
                                     baseId, impl, newVertPositions);
        centroid = centroid - dot(centroid - planePt, nrm) * nrm;
      }
    }

    // Incident chords from the members' producing pairs...
    ChordCrossing cluster{centroid, -1, {}, {}};
    std::set<int> seenChords;
    auto addChord = [&](int ch, double t) {
      if (!seenChords.insert(ch).second) return;
      cluster.chords.push_back(ch);
      cluster.ts.push_back(t);
    };
    for (const uint32_t m : members) {
      addChord(raw[m].chordA, raw[m].tA);
      addChord(raw[m].chordB, raw[m].tB);
    }
    // ...plus eager propagation: every chord incident to any involved
    // face that the merged position lies on (point-to-segment <= eps,
    // the pass-0 endpoint-zone t-guard re-applied), so a k-fold point
    // lands on all k chords even when a pairwise crossing was missed.
    std::set<int> faces;
    for (const uint32_t m : members) {
      faces.insert(rcFaces[m].begin(), rcFaces[m].end());
    }
    const double eps2 = eps * eps;
    for (const int f : faces) {
      if (f < 0 || static_cast<size_t>(f) >= chordsByFace.size()) continue;
      for (const int ch : chordsByFace[f]) {
        if (seenChords.count(ch) > 0) continue;
        const PiercedNewEdge& e = chords[ch].edge;
        const vec3 a = GetPos3(e.v0, baseId, impl, newVertPositions);
        const vec3 b = GetPos3(e.v1, baseId, impl, newVertPositions);
        const vec3 ab = b - a;
        const double abLen2 = dot(ab, ab);
        if (abLen2 == 0) continue;
        const double tGuard = snap / std::sqrt(abLen2);
        const double t = dot(centroid - a, ab) / abLen2;
        if (t <= tGuard || t >= 1.0 - tGuard) continue;
        const vec3 closest = a + ab * t;
        const vec3 dv = centroid - closest;
        if (dot(dv, dv) > eps2) continue;
        addChord(ch, t);
      }
    }
    out.push_back(std::move(cluster));
  }
  return out;
}

void PropagateNewVertsToOnEdgeLists(
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<int>& resolvedIds, const std::vector<Edge>& edges,
    std::vector<EdgeVertList>& onEdgeLists) {
  // For each etIsect, add the resolved vert id to the on-edge list of the
  // piercing edge with parameter t = x.s. Skip if the vert is the edge's
  // endpoint or already in the list.
  std::vector<int> touched;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    auto& list = onEdgeLists[x.edgeIdx];
    const int v = resolvedIds[i];
    if (v == edges[x.edgeIdx].v0 || v == edges[x.edgeIdx].v1) continue;
    if (std::find(list.verts.begin(), list.verts.end(), v) != list.verts.end())
      continue;
    list.verts.push_back(v);
    list.ts.push_back(x.s);
    touched.push_back(x.edgeIdx);
  }
  // BuildOnEdgeVertLists left each list sorted by t, but the appends above
  // are unordered. Step 11 (BuildPerTriHalfedgeGraphs) consumes verts in
  // stored order to build consecutive sub-edges assuming monotone t, so
  // re-sort each touched list by t (carrying verts along).
  std::sort(touched.begin(), touched.end());
  touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
  for (int e : touched) {
    auto& list = onEdgeLists[e];
    std::vector<size_t> perm(list.verts.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(perm.begin(), perm.end(), [&](size_t i, size_t j) {
      return list.ts[i] < list.ts[j];
    });
    std::vector<int> sortedV(list.verts.size());
    std::vector<double> sortedT(list.ts.size());
    for (size_t k = 0; k < perm.size(); ++k) {
      sortedV[k] = list.verts[perm[k]];
      sortedT[k] = list.ts[perm[k]];
    }
    list.verts = std::move(sortedV);
    list.ts = std::move(sortedT);
  }
}

double SegmentPiercesTriInterior(vec3 a, vec3 b, vec3 v0, vec3 v1, vec3 v2,
                                 double relTol) {
  using la::cross;
  using la::dot;
  const vec3 e1 = v1 - v0;
  const vec3 e2 = v2 - v0;
  const vec3 n = cross(e1, e2);
  const double nMag = std::sqrt(dot(n, n));
  if (nMag == 0) return 0.0;
  const double scaleA = std::sqrt(dot(a - v0, a - v0));
  const double scaleB = std::sqrt(dot(b - v0, b - v0));
  const double dTol = relTol * nMag * std::max(scaleA, scaleB);
  const double dA = dot(a - v0, n);
  const double dB = dot(b - v0, n);
  if (std::fabs(dA) < dTol || std::fabs(dB) < dTol) return 0.0;
  if (dA * dB >= 0) return 0.0;
  const double t = dA / (dA - dB);
  const vec3 p = a + t * (b - a);
  const vec3 c0 = cross(v1 - v0, p - v0);
  const vec3 c1 = cross(v2 - v1, p - v1);
  const vec3 c2 = cross(v0 - v2, p - v2);
  const double d0 = dot(c0, n);
  const double d1 = dot(c1, n);
  const double d2 = dot(c2, n);
  const double bTol = relTol * nMag * nMag;
  if (!(d0 > bTol && d1 > bTol && d2 > bTol)) return 0.0;
  return std::min(std::fabs(dA), std::fabs(dB)) / nMag;
}

SelfIntersectionResult CheckSelfIntersection(const Manifold& m, double relTol) {
  SelfIntersectionResult r{};
  if (m.IsEmpty()) return r;
  MeshGL64 mesh = m.GetMeshGL64();
  const size_t nTri = mesh.NumTri();
  r.trianglesTotal = static_cast<int>(nTri);
  if (nTri < 2) return r;

  std::vector<Box> triBoxes(nTri);
  std::vector<std::array<int, 3>> triIdx(nTri);
  for (size_t t = 0; t < nTri; ++t) {
    const int i0 = mesh.triVerts[3 * t + 0];
    const int i1 = mesh.triVerts[3 * t + 1];
    const int i2 = mesh.triVerts[3 * t + 2];
    triIdx[t] = {i0, i1, i2};
    auto vp = [&](int idx) {
      return vec3(mesh.vertProperties[mesh.numProp * idx + 0],
                  mesh.vertProperties[mesh.numProp * idx + 1],
                  mesh.vertProperties[mesh.numProp * idx + 2]);
    };
    Box b(vp(i0), vp(i1));
    b.Union(vp(i2));
    triBoxes[t] = b;
  }

  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));

  auto vp = [&](int idx) {
    return vec3(mesh.vertProperties[mesh.numProp * idx + 0],
                mesh.vertProperties[mesh.numProp * idx + 1],
                mesh.vertProperties[mesh.numProp * idx + 2]);
  };

  auto checkPair = [&](size_t qi, size_t li) {
    if (qi >= li) return;
    const size_t ta = bvh.perm[qi];
    const size_t tb = bvh.perm[li];
    int shared = 0;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (triIdx[ta][a] == triIdx[tb][b]) ++shared;
    if (shared >= 2) {
      ++r.adjacentPairsSkipped;
      return;
    }
    ++r.candidatesChecked;
    const vec3 a0 = vp(triIdx[ta][0]), a1 = vp(triIdx[ta][1]),
               a2 = vp(triIdx[ta][2]);
    const vec3 b0 = vp(triIdx[tb][0]), b1 = vp(triIdx[tb][1]),
               b2 = vp(triIdx[tb][2]);
    double maxMag = 0.0;
    auto take = [&](double m) {
      if (m > maxMag) maxMag = m;
    };
    take(SegmentPiercesTriInterior(a0, a1, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(a1, a2, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(a2, a0, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(b0, b1, a0, a1, a2, relTol));
    take(SegmentPiercesTriInterior(b1, b2, a0, a1, a2, relTol));
    take(SegmentPiercesTriInterior(b2, b0, a0, a1, a2, relTol));
    if (maxMag > 0.0) {
      ++r.interiorPierces;
      if (maxMag > r.maxPierceMagnitude) r.maxPierceMagnitude = maxMag;
    }
  };

  auto recorder = MakeSimpleRecorder(checkPair);
  auto qf = [&](int i) { return bvh.boxes[i]; };
  bvh.collider.Collisions<false>(recorder, qf, static_cast<int>(nTri),
                                 /*parallel=*/false);
  return r;
}

namespace {
// Pipeline body - separated from RunOverlapRemoval so the entry
// point can wrap it in a try/catch and fall back to input on any
// internal exception (= e.g. Triangulate's CCW-check assertion
// when the polygon walker emits a degenerate sub-polygon under
// MANIFOLD_DEBUG builds).
Manifold RunOverlapRemovalImpl(const Manifold& input, double eps);
}  // namespace

Manifold RunOverlapRemoval(const Manifold& input, double eps) {
  // Outer try/catch: if any internal stage throws (= a MANIFOLD_DEBUG
  // assertion in Triangulate, the Manifold(out) constructor, or
  // std::bad_alloc), return the input unchanged. The input is already a
  // valid manifold with no more self-intersections than itself, so this
  // preserves pierce-monotonicity without running any further allocating
  // or possibly-throwing work (MergeVertsEps both asserts and can produce
  // a non-manifold result, so it must not run on the failure path).
  try {
    return RunOverlapRemovalImpl(input, eps);
  } catch (...) {
    return input;
  }
}

namespace {
Manifold RunOverlapRemovalImpl(const Manifold& input, double eps) {
  // DEMOLISHED for the faithful step-9 rewrite: steps 9-13 and the recovery
  // scaffold (pair-sym phases, cap walker, pierce reducers, gate) were removed.
  // Steps 1-8 (the arrangement front) remain and are exercised by the
  // arrangement test harness; RemoveSelfIntersections is a no-op stub until
  // the faithful arrangement + classification is rebuilt.
  (void)eps;
  return input;
}
}  // namespace
}  // namespace overlap_removal
}  // namespace manifold
