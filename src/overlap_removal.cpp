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

std::vector<PerTriHalfedgeGraph> BuildPerTriHalfedgeGraphs(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<NewEdgeWithExtras>& newEdgesWithExtras) {
  const size_t nT = impl.NumTri();
  std::vector<PerTriHalfedgeGraph> out(nT);
  for (size_t t = 0; t < nT; ++t) out[t].triId = static_cast<int>(t);

  std::vector<int> halfedgeToEdge(impl.halfedge_.size(), -1);
  for (size_t e = 0; e < edges.size(); ++e) {
    halfedgeToEdge[edges[e].halfedgeForward] = static_cast<int>(e);
    if (edges[e].halfedgePaired >= 0)
      halfedgeToEdge[edges[e].halfedgePaired] = static_cast<int>(e);
  }

  // Original edges -> one halfedge per sub-edge in T's CCW direction.
  for (size_t t = 0; t < nT; ++t) {
    for (int k : {0, 1, 2}) {
      const int h = static_cast<int>(3 * t + k);
      const Halfedge he = impl.halfedge_.Get(h);
      const int eIdx = halfedgeToEdge[h];
      if (eIdx < 0) continue;
      const Edge& edge = edges[eIdx];
      const auto& list = onEdgeLists[eIdx];
      const bool forward = (he.startVert == edge.v0);
      std::vector<int> seq;
      seq.reserve(list.verts.size() + 2);
      seq.push_back(he.startVert);
      if (forward) {
        for (int v : list.verts) seq.push_back(v);
      } else {
        for (auto it = list.verts.rbegin(); it != list.verts.rend(); ++it)
          seq.push_back(*it);
      }
      seq.push_back(he.endVert);
      for (size_t i = 0; i + 1 < seq.size(); ++i) {
        out[t].halfedges.push_back({seq[i], seq[i + 1], false});
        out[t].verts.insert(seq[i]);
        out[t].verts.insert(seq[i + 1]);
      }
    }
  }

  // New edges: each contributes BOTH directions, sub-divided by extras.
  for (const auto& nwe : newEdgesWithExtras) {
    const PiercedNewEdge& edge = nwe.edge;
    std::vector<int> seq;
    seq.reserve(nwe.extraVerts.size() + 2);
    seq.push_back(edge.v0);
    for (int v : nwe.extraVerts) seq.push_back(v);
    seq.push_back(edge.v1);
    for (int triId : {edge.triA, edge.triB}) {
      auto& g = out[triId];
      for (size_t i = 0; i + 1 < seq.size(); ++i) {
        g.halfedges.push_back({seq[i], seq[i + 1], true});
        g.halfedges.push_back({seq[i + 1], seq[i], true});
        g.verts.insert(seq[i]);
        g.verts.insert(seq[i + 1]);
      }
    }
  }

  return out;
}

void AddNextPointers(const Manifold::Impl& impl,
                     const std::vector<vec3>& newVertPositions,
                     PerTriHalfedgeGraph& g) {
  using la::dot;
  if (g.halfedges.empty()) {
    g.nextHalfedge.clear();
    return;
  }
  const int baseId = static_cast<int>(impl.NumVert());
  auto getPos3 = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Pick projection axis: drop the largest |normal| component. After
  // dropping `dropAxis`, swap the remaining axes if normal sign is
  // negative (= preserves CCW orientation in the 2D projection).
  const vec3 n = impl.faceNormal_[g.triId];
  const int dropAxis =
      (std::fabs(n.x) >= std::fabs(n.y) && std::fabs(n.x) >= std::fabs(n.z)) ? 0
      : (std::fabs(n.y) >= std::fabs(n.z))                                   ? 1
                                           : 2;
  const int axA = (dropAxis + 1) % 3;
  const int axB = (dropAxis + 2) % 3;
  const bool swap = n[dropAxis] < 0;
  auto getPos2 = [&](int id) {
    const vec3 p = getPos3(id);
    if (swap) return std::pair<double, double>(p[axB], p[axA]);
    return std::pair<double, double>(p[axA], p[axB]);
  };

  std::map<int, std::vector<int>> outgoingByVert;
  for (size_t i = 0; i < g.halfedges.size(); ++i) {
    outgoingByVert[g.halfedges[i].startVert].push_back(static_cast<int>(i));
  }

  std::vector<double> outAngle(g.halfedges.size(), 0);
  for (auto& [v, indices] : outgoingByVert) {
    auto [vx, vy] = getPos2(v);
    for (int idx : indices) {
      auto [ex, ey] = getPos2(g.halfedges[idx].endVert);
      outAngle[idx] = std::atan2(ey - vy, ex - vx);
    }
    std::sort(indices.begin(), indices.end(),
              [&](int a, int b) { return outAngle[a] < outAngle[b]; });
  }

  // For each h, next(h) = outgoing-from-h.endVert with smallest CW
  // angle from theta_rev = h.angle + pi. CW (not CCW) is correct for
  // "face on the LEFT" walk.
  g.nextHalfedge.assign(g.halfedges.size(), -1);
  constexpr double kTwoPi = 6.283185307179586;
  for (size_t i = 0; i < g.halfedges.size(); ++i) {
    const int v = g.halfedges[i].endVert;
    auto& list = outgoingByVert[v];
    if (list.empty()) continue;
    double thetaRev = outAngle[i] + 3.141592653589793;
    if (thetaRev > 3.141592653589793) thetaRev -= kTwoPi;
    int best = -1;
    double bestDelta = kTwoPi + 1.0;
    for (int idx : list) {
      if (g.halfedges[idx].endVert == g.halfedges[i].startVert) continue;
      double delta = thetaRev - outAngle[idx];
      while (delta <= 0) delta += kTwoPi;
      while (delta > kTwoPi) delta -= kTwoPi;
      if (delta < bestDelta) {
        bestDelta = delta;
        best = idx;
      }
    }
    g.nextHalfedge[i] = best;
  }
}

void AddNextPointers(const Manifold::Impl& impl,
                     const std::vector<vec3>& newVertPositions,
                     std::vector<PerTriHalfedgeGraph>& graphs) {
  for (auto& g : graphs) AddNextPointers(impl, newVertPositions, g);
}

PolygonWalkResult WalkPolygons(const PerTriHalfedgeGraph& g) {
  PolygonWalkResult r;
  if (g.halfedges.empty() || g.nextHalfedge.size() != g.halfedges.size())
    return r;
  std::vector<bool> visited(g.halfedges.size(), false);
  for (size_t start = 0; start < g.halfedges.size(); ++start) {
    if (visited[start]) continue;
    std::vector<int> polygon;
    int cur = static_cast<int>(start);
    bool stalled = false;
    const size_t maxSteps = g.halfedges.size() + 1;
    for (size_t step = 0; step < maxSteps; ++step) {
      if (cur < 0) {
        stalled = true;
        ++r.stalledHalfedges;
        break;
      }
      if (visited[cur]) break;
      visited[cur] = true;
      polygon.push_back(g.halfedges[cur].startVert);
      cur = g.nextHalfedge[cur];
    }
    // Only accept walks that returned to `start` (a closed cycle). A walk
    // that ran into a different already-visited halfedge is an open path
    // and must not be emitted as a polygon - Triangulate would treat it as
    // closed. The cap walker applies the same `cur == start` guard.
    const bool closed = (cur == static_cast<int>(start));
    if (!stalled && closed && polygon.size() >= 3) {
      r.polygons.push_back(std::move(polygon));
    } else if (!stalled && closed && polygon.size() == 2) {
      r.degeneratePolygons.push_back(std::move(polygon));
    } else if (!stalled) {
      r.stalledHalfedges += static_cast<int>(polygon.size());
    }
  }
  return r;
}

std::vector<PolygonWalkResult> WalkPolygons(
    const std::vector<PerTriHalfedgeGraph>& graphs) {
  std::vector<PolygonWalkResult> out;
  out.reserve(graphs.size());
  for (const auto& g : graphs) out.push_back(WalkPolygons(g));
  return out;
}

ChordPartnerMap BuildChordPartnerMap(
    const std::vector<PiercedNewEdge>& newEdges) {
  ChordPartnerMap r;
  for (const auto& e : newEdges) {
    int a = e.v0, b = e.v1;
    if (a > b) std::swap(a, b);
    r.partnerOf[{a, b, e.triA}] = e.triB;
    r.partnerOf[{a, b, e.triB}] = e.triA;
  }
  return r;
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

bool AnalyticalKeep(int triId, const std::vector<int>& polyVerts,
                    const ChordPartnerMap& chordPartners,
                    const std::vector<vec3>& positions3D, int baseId,
                    const Manifold::Impl& impl) {
  using la::dot;
  std::set<int> boundingPartners;
  std::set<int> chordVerts;
  const size_t n = polyVerts.size();
  for (size_t i = 0; i < n; ++i) {
    int a = polyVerts[i], b = polyVerts[(i + 1) % n];
    int aS = a, bS = b;
    if (aS > bS) std::swap(aS, bS);
    auto it = chordPartners.partnerOf.find({aS, bS, triId});
    if (it != chordPartners.partnerOf.end()) {
      boundingPartners.insert(it->second);
      chordVerts.insert(a);
      chordVerts.insert(b);
    }
  }
  if (boundingPartners.empty()) return true;

  // Conservative drop: require ALL non-chord polygon verts to be on
  // the outward side of the partner. Falls back to centroid only
  // when polygon is entirely chord-verts.
  std::vector<vec3> nonChordPts;
  for (int v : polyVerts) {
    if (chordVerts.count(v) == 0) {
      vec3 pt;
      if (v < baseId) {
        pt = impl.vertPos_[v];
      } else {
        const int j = v - baseId;
        if (j < 0 || j >= static_cast<int>(positions3D.size())) continue;
        pt = positions3D[j];
      }
      nonChordPts.push_back(pt);
    }
  }
  vec3 fallbackPt(0, 0, 0);
  if (nonChordPts.empty()) {
    int contributed = 0;
    for (int v : polyVerts) {
      if (v < baseId) {
        fallbackPt += impl.vertPos_[v];
        ++contributed;
      } else {
        const int j = v - baseId;
        if (j >= 0 && j < static_cast<int>(positions3D.size())) {
          fallbackPt += positions3D[j];
          ++contributed;
        }
      }
    }
    if (contributed > 0) {
      fallbackPt /= static_cast<double>(contributed);
      nonChordPts.push_back(fallbackPt);
    }
  }
  // Scale-relative threshold: bbox.Scale x 1e-6 ignores borderline.
  const double thresh = impl.bBox_.Scale() * 1e-6;
  for (int triB : boundingPartners) {
    const int v0 = impl.halfedge_.Start(3 * triB);
    const vec3 p_B = impl.vertPos_[v0];
    const vec3 n_B = impl.faceNormal_[triB];
    const double nMag = std::sqrt(dot(n_B, n_B));
    if (nMag == 0) continue;
    bool allOutward = true;
    for (const auto& p : nonChordPts) {
      if (dot(p - p_B, n_B) <= thresh * nMag) {
        allOutward = false;
        break;
      }
    }
    if (allOutward) return false;
  }
  return true;
}

namespace {
// sortedTriple helper used in DoCapPass forbidden check.
inline std::array<int, 3> SortedTriple(int a, int b, int c) {
  std::array<int, 3> t = {a, b, c};
  std::sort(t.begin(), t.end());
  return t;
}
}  // namespace

std::pair<int, int> DoCapPass(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples) {
  // Recompute current edge incidence + directed halfedge dirs.
  std::map<std::pair<int, int>, int> ec;
  std::map<std::pair<int, int>, int> dirCount;
  for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
    int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                static_cast<int>(out.triVerts[3 * t + 1]),
                static_cast<int>(out.triVerts[3 * t + 2])};
    for (int e : {0, 1, 2}) {
      int a = v[e], b = v[(e + 1) % 3];
      int sa = a, sb = b;
      if (sa > sb) std::swap(sa, sb);
      ++ec[{sa, sb}];
      ++dirCount[{a, b}];
    }
  }
  // Pierce-aware cap: build a BVH over current out.triVerts and
  // pierce-check each candidate fan/ear tri before emitting.
  auto getVertPos = [&](int v) -> vec3 {
    return vec3(out.vertProperties[3 * v + 0], out.vertProperties[3 * v + 1],
                out.vertProperties[3 * v + 2]);
  };
  SortedBVH capBVH;
  std::vector<std::array<int, 3>> capTris;
  auto rebuildCapBVH = [&]() {
    const size_t nT = out.triVerts.size() / 3;
    capTris.assign(nT, {0, 0, 0});
    std::vector<Box> triBoxes(nT);
    for (size_t t = 0; t < nT; ++t) {
      const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
      const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
      const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
      capTris[t] = {i0, i1, i2};
      Box b(getVertPos(i0), getVertPos(i1));
      b.Union(getVertPos(i2));
      triBoxes[t] = b;
    }
    capBVH =
        BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));
  };
  rebuildCapBVH();
  // Resync the symmetric edge-incidence counts from the current
  // out.triVerts. Called at the top of each cap cycle so cross-cycle counts
  // stay exact regardless of within-cycle fan/ear/rollback bookkeeping
  // (the ear-clip rollback and last-ear paths do not perfectly maintain ec;
  // this mirrors the rebuildCapBVH() resync).
  auto rebuildEc = [&]() {
    ec.clear();
    for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
      const int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                        static_cast<int>(out.triVerts[3 * t + 1]),
                        static_cast<int>(out.triVerts[3 * t + 2])};
      for (int e : {0, 1, 2}) {
        int sa = v[e], sb = v[(e + 1) % 3];
        if (sa > sb) std::swap(sa, sb);
        ++ec[{sa, sb}];
      }
    }
  };

  auto wouldPierce = [&](int a, int b, int c) -> bool {
    if (capTris.empty()) return false;
    const vec3 va = getVertPos(a);
    const vec3 vb = getVertPos(b);
    const vec3 vc = getVertPos(c);
    Box queryBox(va, vb);
    queryBox.Union(vc);
    bool found = false;
    auto recorderf = [&](int /*qi*/, int li) {
      if (found) return;
      const int origTri = static_cast<int>(capBVH.perm[li]);
      const int i0 = capTris[origTri][0];
      const int i1 = capTris[origTri][1];
      const int i2 = capTris[origTri][2];
      int shared = 0;
      for (int x : {a, b, c})
        for (int y : {i0, i1, i2})
          if (x == y) ++shared;
      if (shared >= 2) return;
      const vec3 vi0 = getVertPos(i0);
      const vec3 vi1 = getVertPos(i1);
      const vec3 vi2 = getVertPos(i2);
      if (SegmentPiercesTriInterior(va, vb, vi0, vi1, vi2) > 0 ||
          SegmentPiercesTriInterior(vb, vc, vi0, vi1, vi2) > 0 ||
          SegmentPiercesTriInterior(vc, va, vi0, vi1, vi2) > 0 ||
          SegmentPiercesTriInterior(vi0, vi1, va, vb, vc) > 0 ||
          SegmentPiercesTriInterior(vi1, vi2, va, vb, vc) > 0 ||
          SegmentPiercesTriInterior(vi2, vi0, va, vb, vc) > 0) {
        found = true;
      }
    };
    auto recorder = MakeSimpleRecorder(recorderf);
    auto f = [&](int) { return queryBox; };
    capBVH.collider.Collisions<false>(recorder, f, 1, /*parallel=*/false);
    return found;
  };

  std::map<std::pair<int, int>, int> k1Dir;
  for (const auto& [edge, c] : ec) {
    if (c != 1) continue;
    if (dirCount.count({edge.first, edge.second}))
      k1Dir[edge] = 1;
    else
      k1Dir[edge] = -1;
  }
  std::map<int, std::vector<int>> capN;
  for (const auto& [edge, kdir] : k1Dir) {
    if (kdir > 0)
      capN[edge.second].push_back(edge.first);
    else
      capN[edge.first].push_back(edge.second);
  }
  std::set<int> visited;
  int closed = 0, addedTris = 0;
  auto dfs = [&](int start) -> std::vector<int> {
    std::vector<int> cyc;
    cyc.push_back(start);
    std::vector<size_t> ni;
    ni.push_back(0);
    std::set<int> onP;
    onP.insert(start);
    int s = 0;
    for (; s < kCycleWalkerGuard && !cyc.empty(); ++s) {
      int cur = cyc.back();
      auto it = capN.find(cur);
      if (it == capN.end()) {
        onP.erase(cur);
        cyc.pop_back();
        ni.pop_back();
        continue;
      }
      if (ni.back() >= it->second.size()) {
        onP.erase(cur);
        cyc.pop_back();
        ni.pop_back();
        continue;
      }
      int nx = it->second[ni.back()++];
      if (nx == start && cyc.size() >= 3) return cyc;
      if (onP.count(nx)) continue;
      cyc.push_back(nx);
      ni.push_back(0);
      onP.insert(nx);
    }
    DEBUG_ASSERT(s < kCycleWalkerGuard, logicErr,
                 "cap-walker DFS hit kCycleWalkerGuard");
    return {};
  };
  auto wouldFanConflict = [&](const std::vector<int>& cyc,
                              size_t apex) -> bool {
    const size_t n = cyc.size();
    for (size_t off = 1; off + 1 < n; ++off) {
      const int vi = cyc[(apex + off) % n];
      const int vj = cyc[(apex + off + 1) % n];
      for (int newEdge : {0, 1}) {
        int va = cyc[apex];
        int vb = (newEdge == 0) ? vi : vj;
        if (va == vb) continue;
        int sa = va, sb = vb;
        if (sa > sb) std::swap(sa, sb);
        auto it = ec.find({sa, sb});
        if (it != ec.end() && it->second >= 2) return true;
      }
    }
    return false;
  };
  auto wouldFanPierce = [&](const std::vector<int>& cyc, size_t apex) -> bool {
    const size_t n = cyc.size();
    for (size_t off = 1; off + 1 < n; ++off) {
      const int vi = cyc[(apex + off) % n];
      const int vj = cyc[(apex + off + 1) % n];
      if (wouldPierce(cyc[apex], vi, vj)) return true;
    }
    return false;
  };
  auto wouldFanForbidden = [&](const std::vector<int>& cyc,
                               size_t apex) -> bool {
    if (forbiddenTriples.empty()) return false;
    const size_t n = cyc.size();
    for (size_t off = 1; off + 1 < n; ++off) {
      const int vi = cyc[(apex + off) % n];
      const int vj = cyc[(apex + off + 1) % n];
      if (forbiddenTriples.count(SortedTriple(cyc[apex], vi, vj))) return true;
    }
    return false;
  };
  for (const auto& [start, _] : capN) {
    if (visited.count(start)) continue;
    auto cyc = dfs(start);
    if (cyc.size() < 3) continue;
    // Resync ec from current geometry before this cycle's conflict checks.
    rebuildEc();
    size_t chosenApex = SIZE_MAX;
    for (size_t a = 0; a < cyc.size(); ++a) {
      if (wouldFanConflict(cyc, a)) continue;
      if (wouldFanPierce(cyc, a)) continue;
      if (wouldFanForbidden(cyc, a)) continue;
      chosenApex = a;
      break;
    }
    if (chosenApex != SIZE_MAX) {
      const size_t n = cyc.size();
      for (size_t off = 1; off + 1 < n; ++off) {
        const int va = cyc[chosenApex];
        const int vb = cyc[(chosenApex + off) % n];
        const int vc = cyc[(chosenApex + off + 1) % n];
        out.triVerts.push_back(va);
        out.triVerts.push_back(vb);
        out.triVerts.push_back(vc);
        ++addedTris;
        auto bumpEc = [&](int x, int y) {
          int sx = x, sy = y;
          if (sx > sy) std::swap(sx, sy);
          ++ec[{sx, sy}];
        };
        bumpEc(va, vb);
        bumpEc(vb, vc);
        bumpEc(va, vc);
      }
      for (int v : cyc) visited.insert(v);
      ++closed;
      rebuildCapBVH();
      continue;
    }
    // Fan failed - try greedy ear-clipping.
    std::vector<int> remaining(cyc.begin(), cyc.end());
    int earTrisAdded = 0;
    int earGuard = 0;
    while (remaining.size() >= 3 && earGuard++ < kEarClipGuard) {
      const size_t m = remaining.size();
      if (m == 3) {
        // Last ear - emit unconditionally. Applying the forbidden-triple
        // check at the final 3-cycle trades a small residual pierce for a
        // large fallback regression, so skip it here.
        out.triVerts.push_back(remaining[0]);
        out.triVerts.push_back(remaining[1]);
        out.triVerts.push_back(remaining[2]);
        ++earTrisAdded;
        ++addedTris;
        remaining.clear();
        break;
      }
      ssize_t earIdx = -1;
      for (size_t i = 0; i < m; ++i) {
        const int va = remaining[i];
        const int vb = remaining[(i + 1) % m];
        const int vc = remaining[(i + 2) % m];
        int sa = va, sc = vc;
        if (sa > sc) std::swap(sa, sc);
        auto it = ec.find({sa, sc});
        if (it != ec.end() && it->second >= 2) continue;
        if (wouldPierce(va, vb, vc)) continue;
        if (forbiddenTriples.count(SortedTriple(va, vb, vc))) continue;
        earIdx = static_cast<ssize_t>(i);
        break;
      }
      if (earIdx < 0) break;
      const int va = remaining[earIdx];
      const int vb = remaining[(earIdx + 1) % m];
      const int vc = remaining[(earIdx + 2) % m];
      out.triVerts.push_back(va);
      out.triVerts.push_back(vb);
      out.triVerts.push_back(vc);
      ++earTrisAdded;
      ++addedTris;
      int sac = va, scc = vc;
      if (sac > scc) std::swap(sac, scc);
      ++ec[{sac, scc}];
      int s1 = va, s2 = vb;
      if (s1 > s2) std::swap(s1, s2);
      ++ec[{s1, s2}];
      int s3 = vb, s4 = vc;
      if (s3 > s4) std::swap(s3, s4);
      ++ec[{s3, s4}];
      remaining.erase(remaining.begin() + (earIdx + 1) % m);
    }
    DEBUG_ASSERT(earGuard < kEarClipGuard, logicErr,
                 "ear-clip walker hit kEarClipGuard");
    if (remaining.empty() || remaining.size() < 3) {
      for (int v : cyc) visited.insert(v);
      ++closed;
      rebuildCapBVH();
    } else {
      // Couldn't fully triangulate: roll back the partial ears.
      for (int t = 0; t < earTrisAdded; ++t) {
        out.triVerts.pop_back();
        out.triVerts.pop_back();
        out.triVerts.pop_back();
        --addedTris;
      }
    }
  }
  return std::pair<int, int>{closed, addedTris};
}

namespace {
// Helper for both pierce reducers: build BVH over current
// out.triVerts and find piercing pairs. Returns (piercePairs,
// pierceCount). tris[] and perm[] are passed back so callers can resolve
// pair indices to vert triplets for forbidden-set management.
struct PierceScanResult {
  std::set<std::pair<int, int>> piercePairs;
  std::vector<int> pierceCount;
  std::vector<std::array<int, 3>> tris;
};

inline PierceScanResult ScanForPiercingPairs(const MeshGL64& out) {
  PierceScanResult r;
  const size_t nT = out.triVerts.size() / 3;
  r.pierceCount.assign(nT, 0);
  r.tris.assign(nT, {0, 0, 0});
  if (nT == 0) return r;
  auto getVPos = [&](int v) -> vec3 {
    return vec3(out.vertProperties[3 * v + 0], out.vertProperties[3 * v + 1],
                out.vertProperties[3 * v + 2]);
  };
  std::vector<Box> triBoxes(nT);
  for (size_t t = 0; t < nT; ++t) {
    const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
    const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
    const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
    r.tris[t] = {i0, i1, i2};
    Box b(getVPos(i0), getVPos(i1));
    b.Union(getVPos(i2));
    triBoxes[t] = b;
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));
  for (size_t qi = 0; qi < nT; ++qi) {
    const auto& qt = r.tris[qi];
    const vec3 qa = getVPos(qt[0]);
    const vec3 qb = getVPos(qt[1]);
    const vec3 qc = getVPos(qt[2]);
    Box queryBox(qa, qb);
    queryBox.Union(qc);
    auto recorderf = [&](int /*qiL*/, int li) {
      const int oi = static_cast<int>(bvh.perm[li]);
      if (oi <= static_cast<int>(qi)) return;
      const auto& ot = r.tris[oi];
      int shared = 0;
      for (int x : qt)
        for (int y : ot)
          if (x == y) ++shared;
      if (shared >= 2) return;
      const vec3 oa = getVPos(ot[0]);
      const vec3 ob = getVPos(ot[1]);
      const vec3 oc = getVPos(ot[2]);
      if (SegmentPiercesTriInterior(qa, qb, oa, ob, oc) > 0 ||
          SegmentPiercesTriInterior(qb, qc, oa, ob, oc) > 0 ||
          SegmentPiercesTriInterior(qc, qa, oa, ob, oc) > 0 ||
          SegmentPiercesTriInterior(oa, ob, qa, qb, qc) > 0 ||
          SegmentPiercesTriInterior(ob, oc, qa, qb, qc) > 0 ||
          SegmentPiercesTriInterior(oc, oa, qa, qb, qc) > 0) {
        r.piercePairs.insert({static_cast<int>(qi), oi});
        ++r.pierceCount[qi];
        ++r.pierceCount[oi];
      }
    };
    auto recorder = MakeSimpleRecorder(recorderf);
    auto qf = [&](int) { return queryBox; };
    bvh.collider.Collisions<false>(recorder, qf, 1, /*parallel=*/false);
  }
  return r;
}

// Common drop logic: pick which tri to drop per pair (higher
// pierce-count wins; tie -> lower idx) and produce dropMask.
inline std::vector<bool> ComputeDropMask(const PierceScanResult& scan,
                                         size_t nT) {
  std::vector<bool> dropMask(nT, false);
  for (const auto& [a, b] : scan.piercePairs) {
    if (dropMask[a] || dropMask[b]) continue;
    if (scan.pierceCount[a] > scan.pierceCount[b])
      dropMask[a] = true;
    else if (scan.pierceCount[b] > scan.pierceCount[a])
      dropMask[b] = true;
    else
      dropMask[a < b ? a : b] = true;
  }
  return dropMask;
}

// Common compaction: rebuild out.triVerts dropping marked tris.
// Returns count dropped.
inline int CompactDroppedTris(MeshGL64& out,
                              const std::vector<bool>& dropMask) {
  const size_t nT = out.triVerts.size() / 3;
  std::vector<uint64_t> newTriVerts;
  newTriVerts.reserve(out.triVerts.size());
  int dropped = 0;
  for (size_t t = 0; t < nT; ++t) {
    if (dropMask[t]) {
      ++dropped;
      continue;
    }
    newTriVerts.push_back(out.triVerts[3 * t + 0]);
    newTriVerts.push_back(out.triVerts[3 * t + 1]);
    newTriVerts.push_back(out.triVerts[3 * t + 2]);
  }
  out.triVerts = std::move(newTriVerts);
  return dropped;
}
}  // namespace

int PierceAwareReducer(MeshGL64& out, int maxIters) {
  int totalDropped = 0;
  int iter = 0;
  for (; iter < maxIters; ++iter) {
    if (out.triVerts.empty()) break;
    auto scan = ScanForPiercingPairs(out);
    if (scan.piercePairs.empty()) break;
    auto dropMask = ComputeDropMask(scan, scan.tris.size());
    int iterDropped = CompactDroppedTris(out, dropMask);
    if (iterDropped == 0) break;
    totalDropped += iterDropped;
  }
  DEBUG_ASSERT(iter < maxIters, logicErr,
               "PierceAwareReducer hit maxIters cap");
  return totalDropped;
}

int PostCapPierceReducer(MeshGL64& out,
                         std::set<std::array<int, 3>>& forbiddenTriples) {
  // No DEBUG_ASSERT on cap-hit: this loop alternates pierce-drop and
  // cap-close, which may oscillate (cap can add tris back). The cap
  // is a per-call budget, not a convergence tripwire.
  std::set<std::pair<int, int>> prevPairs;
  int totalDropped = 0;
  for (int iter = 0; iter < kPostCapIters; ++iter) {
    if (out.triVerts.empty()) break;
    auto scan = ScanForPiercingPairs(out);
    if (scan.piercePairs.empty()) break;
    const bool stuck = (iter > 0 && scan.piercePairs == prevPairs);
    prevPairs = scan.piercePairs;
    auto dropMask = ComputeDropMask(scan, scan.tris.size());
    if (stuck) {
      for (size_t t = 0; t < scan.tris.size(); ++t) {
        if (!dropMask[t]) continue;
        forbiddenTriples.insert(
            SortedTriple(scan.tris[t][0], scan.tris[t][1], scan.tris[t][2]));
      }
    }
    int iterDropped = CompactDroppedTris(out, dropMask);
    if (iterDropped == 0) break;
    totalDropped += iterDropped;
    (void)DoCapPass(out, forbiddenTriples);
  }
  return totalDropped;
}

namespace {
// Compute current edge-incidence map (sorted edge -> count) from
// out.triVerts. Used by both DropExcessHalfedgeContributors and TrimOrphans.
inline std::map<std::pair<int, int>, int> ComputeEdgeIncidence(
    const MeshGL64& out) {
  std::map<std::pair<int, int>, int> ec;
  for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
    int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                static_cast<int>(out.triVerts[3 * t + 1]),
                static_cast<int>(out.triVerts[3 * t + 2])};
    for (int e : {0, 1, 2}) {
      int a = v[e], b = v[(e + 1) % 3];
      if (a > b) std::swap(a, b);
      ++ec[{a, b}];
    }
  }
  return ec;
}
}  // namespace

EdgeReducerResult DropExcessHalfedgeContributors(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples,
    int outerMax) {
  EdgeReducerResult r;
  // Quick check: any edge with k>2?
  {
    auto ec = ComputeEdgeIncidence(out);
    int k2plus = 0;
    for (const auto& [_, c] : ec)
      if (c > 2) ++k2plus;
    if (k2plus == 0) return r;
  }
  int prevDropped = -1, prevK2Plus = -1;
  int outer = 0;
  for (; outer < outerMax; ++outer) {
    const size_t nTris = out.triVerts.size() / 3;
    if (nTris == 0) break;
    std::vector<bool> dropped(nTris, false);
    int innerDropped = 0;
    {
      // Per-iteration: drop ALL excess directional contributors
      // simultaneously. heMap[(a,b)] = list of tri indices emitting
      // halfedge (a->b). Score-based keeper pick (= prefer keeping
      // the most "load-bearing" tri).
      std::map<std::pair<int, int>, std::vector<int>> heMap;
      auto ecCur = ComputeEdgeIncidence(out);
      for (size_t t = 0; t < nTris; ++t) {
        int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                    static_cast<int>(out.triVerts[3 * t + 1]),
                    static_cast<int>(out.triVerts[3 * t + 2])};
        for (int e : {0, 1, 2}) {
          heMap[{v[e], v[(e + 1) % 3]}].push_back(static_cast<int>(t));
        }
      }
      auto dropCost = [&](int t) {
        int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                    static_cast<int>(out.triVerts[3 * t + 1]),
                    static_cast<int>(out.triVerts[3 * t + 2])};
        int cost = 0;
        for (int e : {0, 1, 2}) {
          int a = v[e], b = v[(e + 1) % 3];
          if (a > b) std::swap(a, b);
          auto it = ecCur.find({a, b});
          if (it != ecCur.end() && it->second == 2) ++cost;
        }
        return cost;
      };
      for (const auto& [he, tris] : heMap) {
        if (tris.size() <= 1) continue;
        int keepIdx = 0;
        int keepCost = dropCost(tris[0]);
        for (size_t i = 1; i < tris.size(); ++i) {
          int c = dropCost(tris[i]);
          if (c > keepCost || (c == keepCost && tris[i] < tris[keepIdx])) {
            keepIdx = static_cast<int>(i);
            keepCost = c;
          }
        }
        for (size_t i = 0; i < tris.size(); ++i) {
          if (static_cast<int>(i) == keepIdx) continue;
          if (!dropped[tris[i]]) {
            dropped[tris[i]] = true;
            ++innerDropped;
          }
        }
      }
    }
    if (innerDropped == 0) break;
    int actualDropped = CompactDroppedTris(out, dropped);
    r.totalDropped += actualDropped;
    auto [c, ct] = DoCapPass(out, forbiddenTriples);
    r.totalCapped += c;
    r.totalCapTris += ct;
    // Convergence: count current k>2. Bail if not decreasing.
    int curK2Plus = 0;
    auto ec = ComputeEdgeIncidence(out);
    for (const auto& [_, kk] : ec)
      if (kk > 2) ++curK2Plus;
    if (prevK2Plus >= 0 && curK2Plus >= prevK2Plus) break;
    prevK2Plus = curK2Plus;
    if (r.totalDropped == prevDropped) break;
    prevDropped = r.totalDropped;
  }
  DEBUG_ASSERT(outer < outerMax, logicErr,
               "DropExcessHalfedgeContributors hit outerMax cap");
  return r;
}

TriangulationResult TriangulateAndEmit(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<PolygonWalkResult>& walks,
    PolygonClassifierFn classifier) {
  using la::cross;
  using la::dot;
  TriangulationResult r{Manifold(), 0, 0, 0};
  const int baseId = static_cast<int>(impl.NumVert());
  auto getPos3 = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Build the combined vert table: original input verts + new verts.
  MeshGL64 out;
  out.numProp = 3;
  const size_t totalVerts = impl.vertPos_.size() + newVertPositions.size();
  out.vertProperties.reserve(3 * totalVerts);
  for (const auto& p : impl.vertPos_) {
    out.vertProperties.push_back(p.x);
    out.vertProperties.push_back(p.y);
    out.vertProperties.push_back(p.z);
  }
  for (const auto& p : newVertPositions) {
    out.vertProperties.push_back(p.x);
    out.vertProperties.push_back(p.y);
    out.vertProperties.push_back(p.z);
  }

  // Phase 1: tentative classification.
  std::vector<std::vector<bool>> keepFlags(walks.size());
  std::vector<bool> autoKeepFlags(walks.size(), false);
  for (size_t triId = 0; triId < walks.size(); ++triId) {
    const auto& w = walks[triId];
    const vec3 n = impl.faceNormal_[triId];
    autoKeepFlags[triId] = (w.polygons.size() == 1);
    keepFlags[triId].assign(w.polygons.size(), true);
    if (classifier && !autoKeepFlags[triId]) {
      for (size_t pi = 0; pi < w.polygons.size(); ++pi) {
        auto c = classifier(w.polygons[pi], n, static_cast<int>(triId));
        keepFlags[triId][pi] = c.keep;
      }
    }
  }
  // Phase 2: cascade-drop forward - drop auto-kept 1-poly tris whose
  // sub-edges are all on the k=1 boundary. Iterate until stable.
  int pass = 0;
  for (; pass < kCascadeDropMaxPass; ++pass) {
    std::map<std::pair<int, int>, int> edgeCount;
    for (size_t triId = 0; triId < walks.size(); ++triId) {
      for (size_t pi = 0; pi < walks[triId].polygons.size(); ++pi) {
        if (!keepFlags[triId][pi]) continue;
        const auto& poly = walks[triId].polygons[pi];
        for (size_t i = 0; i < poly.size(); ++i) {
          int a = poly[i], b = poly[(i + 1) % poly.size()];
          if (a > b) std::swap(a, b);
          ++edgeCount[{a, b}];
        }
      }
    }
    int changed = 0;
    for (size_t triId = 0; triId < walks.size(); ++triId) {
      if (!autoKeepFlags[triId]) continue;
      if (walks[triId].polygons.empty()) continue;
      if (!keepFlags[triId][0]) continue;
      const auto& poly = walks[triId].polygons[0];
      int k1 = 0, k2 = 0;
      for (size_t i = 0; i < poly.size(); ++i) {
        int a = poly[i], b = poly[(i + 1) % poly.size()];
        if (a > b) std::swap(a, b);
        auto it = edgeCount.find({a, b});
        if (it != edgeCount.end()) {
          if (it->second == 1)
            ++k1;
          else if (it->second == 2)
            ++k2;
        }
      }
      if (k1 > k2 && k1 > 0) {
        keepFlags[triId][0] = false;
        ++changed;
      }
    }
    if (!changed) break;
  }
  DEBUG_ASSERT(pass < kCascadeDropMaxPass, logicErr,
               "TriangulateAndEmit cascade-drop hit kCascadeDropMaxPass");

  // Emit loop: per tri, per kept polygon, triangulate via on-edge
  // collinear fan or manifold::Triangulate.
  for (size_t triId = 0; triId < walks.size(); ++triId) {
    const auto& w = walks[triId];
    const vec3 n = impl.faceNormal_[triId];
    const int dropAxis =
        (std::fabs(n.x) >= std::fabs(n.y) && std::fabs(n.x) >= std::fabs(n.z))
            ? 0
        : (std::fabs(n.y) >= std::fabs(n.z)) ? 1
                                             : 2;
    const int axA = (dropAxis + 1) % 3;
    const int axB = (dropAxis + 2) % 3;
    const bool swap = n[dropAxis] < 0;
    auto getPos2 = [&](int id) {
      const vec3 p = getPos3(id);
      if (swap) return vec2(p[axB], p[axA]);
      return vec2(p[axA], p[axB]);
    };

    for (size_t pi = 0; pi < w.polygons.size(); ++pi) {
      const auto& poly = w.polygons[pi];
      ++r.polygonsTriangulated;
      const bool keep = keepFlags[triId][pi];
      if (!keep) {
        ++r.polygonsDropped;
        continue;
      }
      // Drop sub-polygons collapsed to a line/point in projection (pinched
      // cycles or coincident verts): ~zero projected area. They would
      // otherwise throw inside Triangulate (caught, but noisy under
      // MANIFOLD_DEBUG) or feed the fan a degenerate apex.
      {
        double area2x = 0.0, ext = 0.0;
        const vec2 p0 = getPos2(poly[0]);
        for (size_t i = 0; i < poly.size(); ++i) {
          const vec2 a = getPos2(poly[i]);
          const vec2 b = getPos2(poly[(i + 1) % poly.size()]);
          area2x += a.x * b.y - b.x * a.y;
          ext = std::max(ext, la::length(a - p0));
        }
        if (std::fabs(area2x) <= ext * ext * 1e-12) {
          ++r.trisDroppedTooSmall;
          continue;
        }
      }
      ++r.polygonsKept;
      if (autoKeepFlags[triId] && classifier) ++r.polygonsAutoKept;
      auto emit = [&](int a, int b, int c) {
        out.triVerts.push_back(a);
        out.triVerts.push_back(b);
        out.triVerts.push_back(c);
        ++r.trianglesEmitted;
      };
      if (poly.size() == 3) {
        emit(poly[0], poly[1], poly[2]);
        continue;
      }
      // On-edge collinear apex detection: a vert V whose two perim
      // neighbors A and B are collinear with V on segment AB. If
      // present, fan from V to suppress the AB chord (which would
      // otherwise create k=4 with the opposite tri sharing AB).
      int fanIdx = -1;
      for (size_t i = 0; i < poly.size(); ++i) {
        const vec3 a = getPos3(poly[(i + poly.size() - 1) % poly.size()]);
        const vec3 v = getPos3(poly[i]);
        const vec3 b = getPos3(poly[(i + 1) % poly.size()]);
        const vec3 ab = b - a;
        const double abLen = la::length(ab);
        if (abLen <= 0) continue;
        const vec3 av = v - a;
        const double cross_mag = la::length(la::cross(ab, av));
        if (cross_mag < abLen * abLen * 1e-9) {
          const double t = la::dot(av, ab) / (abLen * abLen);
          if (t > 1e-6 && t < 1.0 - 1e-6) {
            fanIdx = static_cast<int>(i);
            break;
          }
        }
      }
      if (fanIdx >= 0) {
        const int v0 = fanIdx;
        const int N = static_cast<int>(poly.size());
        // The fan is valid only if no triangle is inverted relative to the
        // polygon's projected orientation (a polygon concave away from the
        // apex would self-overlap). The collinear apex itself is allowed
        // (zero-area tris); only strictly-opposite tris veto. If vetoed,
        // fall through to Triangulate, which handles concave polygons.
        double polyArea2 = 0.0;
        for (int i = 0; i < N; ++i) {
          const vec2 a = getPos2(poly[i]);
          const vec2 b = getPos2(poly[(i + 1) % N]);
          polyArea2 += a.x * b.y - b.x * a.y;
        }
        bool fanSimple = true;
        for (int i = 1; i + 1 < N && fanSimple; ++i) {
          const vec2 a = getPos2(poly[v0]);
          const vec2 b = getPos2(poly[(v0 + i) % N]);
          const vec2 c = getPos2(poly[(v0 + i + 1) % N]);
          const double tri2 =
              (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
          if (tri2 * polyArea2 < 0) fanSimple = false;
        }
        if (fanSimple) {
          for (int i = 1; i + 1 < N; ++i)
            emit(poly[v0], poly[(v0 + i) % N], poly[(v0 + i + 1) % N]);
          continue;
        }
      }
      SimplePolygon poly2;
      poly2.reserve(poly.size());
      for (int v : poly) poly2.push_back(getPos2(v));
      Polygons polys = {poly2};
      try {
        auto tris = Triangulate(polys);
        for (const auto& t : tris) emit(poly[t.x], poly[t.y], poly[t.z]);
      } catch (...) {
        ++r.trisDroppedTooSmall;
      }
    }
  }

  // Pre-cap pierce-aware reducer: drop classifier-output overlapping tris.
  std::set<std::array<int, 3>> forbiddenTriples;
  PierceAwareReducer(out);

  // Surface cap + downstream cleanup. Run the cap walker only when there
  // are k=1 (open boundary) cycles, gated by a quick edge-incidence scan.
  {
    int k1 = 0;
    for (const auto& [_, c] : ComputeEdgeIncidence(out)) {
      if (c == 1) ++k1;
    }
    if (k1 > 0) {
      DoCapPass(out, forbiddenTriples);
      DropExcessHalfedgeContributors(out, forbiddenTriples);
      TrimOrphans(out, forbiddenTriples);
      PostCapPierceReducer(out, forbiddenTriples);
    }
  }

  r.output = Manifold(out);
  if (r.output.Status() != Manifold::Error::NoError) {
    // Last-chance recovery: TrimOrphans (default-on) plus an extra
    // DropExcessHalfedgeContributors pass for k>2 edges that the
    // triangulator may have produced. If this still doesn't
    // construct, the caller's pierce/drift gate falls back to input.
    DropExcessHalfedgeContributors(out, forbiddenTriples);
    TrimOrphans(out, forbiddenTriples);
    r.output = Manifold(out);
  }
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
  if (input.IsEmpty()) return input;
  if (eps <= 0.0) eps = InferEps(input);

  // Step 1: eps-merge close verts.
  auto mr = MergeVertsEps(input, eps);

  // Steps 2-8: setup + pipeline structural stages.
  auto impl = ImplFromManifold(mr.manifold);
  // Defensive: if the input's inherited tolerance is wildly larger
  // than its bounding box (= polluted by an upstream Boolean op
  // involving extreme-scale geometry, e.g. cray's Subtract from a
  // [-FLT_MAX, FLT_MAX]^3 cuboid leaves tolerance=3.4e+26 on a
  // ~25400-unit-bbox result), reset to a sensible per-bbox value.
  // Without this, every geometric predicate that consults
  // impl.tolerance_ treats all positions as coincident and the
  // pipeline produces meaningless windings.
  if (impl.bBox_.IsFinite() && impl.tolerance_ > impl.bBox_.Scale()) {
    // SetEpsilon MAXes with the existing tolerance_, so reset both
    // to 0 first so the freshly-computed bBox-derived value wins.
    impl.tolerance_ = 0;
    impl.epsilon_ = 0;
    impl.SetEpsilon(-1, false);
  }
  auto edges = EnumerateEdges(impl);
  auto onEdgeLists = BuildOnEdgeVertLists(impl, edges, eps);
  auto onTriLists = BuildOnTriVertLists(impl, eps);
  auto etIsects =
      FindEdgeTriIntersections(impl, edges, onEdgeLists, onTriLists, eps);
  auto chordOutput = GenerateChordEdges(impl, edges, etIsects, eps);
  auto chordsWithExtras =
      AddInteriorVertsToNewEdges(impl, chordOutput.newVertPositions,
                                 chordOutput.newEdges, onTriLists, eps);
  PropagateNewVertsToOnEdgeLists(etIsects, chordOutput.resolvedIds, edges,
                                 onEdgeLists);

  // Steps 11p1+11p2+11p3: per-tri halfedge graphs + walk polygons.
  auto halfedgeGraphs =
      BuildPerTriHalfedgeGraphs(impl, edges, onEdgeLists, chordsWithExtras);
  AddNextPointers(impl, chordOutput.newVertPositions, halfedgeGraphs);
  auto polygonWalks = WalkPolygons(halfedgeGraphs);

  const int baseId = static_cast<int>(impl.NumVert());
  auto chordPartners = BuildChordPartnerMap(chordOutput.newEdges);

  // AnalyzeSelfMesh per-vert two-sided winding (classifier ground truth).
  Vec<std::array<int, 2>> sma_p1q2;
  sma_p1q2.reserve(etIsects.size());
  for (const auto& x : etIsects) {
    if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size())) continue;
    sma_p1q2.push_back({edges[x.edgeIdx].halfedgeForward, x.triIdx});
  }
  std::sort(sma_p1q2.begin(), sma_p1q2.end(),
            [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
              return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
            });
  auto sma = AnalyzeSelfMesh(impl, VecView<const std::array<int, 2>>(
                                       sma_p1q2.data(), sma_p1q2.size()));
  // Per-polygon two-sided winding classifier (with centroid-probe
  // fallback when per-vert windings are all zero or undecidable).
  // kProbeRayDir / kProbeOffsetCoeff / kProbeRayLengthCoeff are
  // shared with AnalyzeSelfMesh in self_mesh_analysis.h so the
  // per-vert and per-polygon probes always use the same direction
  // and scaling.
  using la::dot;
  const double meshScale = ProbeMeshScale(impl);
  const double probeEps = meshScale * kProbeOffsetCoeff;
  const double rayLen = meshScale * kProbeRayLengthCoeff;
  auto keepIfOnSurface = [&](int triId, const std::vector<int>& poly) -> bool {
    int boundaryV = 0, interiorV = 0, exteriorV = 0;
    for (int v : poly) {
      if (v < 0 || v >= static_cast<int>(sma.w_above.size())) continue;
      int wa = sma.w_above[v];
      int wb = sma.w_below[v];
      if ((wa == 0 && wb >= 1) || (wa >= 1 && wb == 0))
        ++boundaryV;
      else if (wa >= 1 && wb >= 1)
        ++interiorV;
      else if (wa <= 0 && wb <= 0)
        ++exteriorV;
    }
    if (boundaryV + interiorV + exteriorV > 0) {
      return boundaryV > 0 && boundaryV >= interiorV && boundaryV >= exteriorV;
    }
    // Per-vert tally was silent (= no verts had a clean (0,1)/(1,0)/
    // both-in/both-out signal). Fall back to centroid-probe: ray-cast
    // eps above and below the polygon centroid along the tri normal.
    vec3 c(0, 0, 0);
    for (int v : poly)
      c += GetPos3(v, baseId, impl, chordOutput.newVertPositions);
    c /= static_cast<double>(poly.size());
    const vec3 n_T = impl.faceNormal_[triId];
    const int wa = WindingAt(impl, c + n_T * probeEps, kProbeRayDir, rayLen);
    const int wb = WindingAt(impl, c - n_T * probeEps, kProbeRayDir, rayLen);
    return (wa == 0 && wb >= 1) || (wa >= 1 && wb == 0);
  };
  std::map<std::pair<int, int>, bool> precomputedKeep;
  for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
    const auto& polys = polygonWalks[triId].polygons;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
      if (polys.size() == 1) {
        precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] = true;
        continue;
      }
      precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] =
          keepIfOnSurface(static_cast<int>(triId), polys[pi]);
    }
  }
  // Degenerate (2-vert) polys: always-dropped.
  for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
    const size_t base = polygonWalks[triId].polygons.size();
    const size_t ndeg = polygonWalks[triId].degeneratePolygons.size();
    for (size_t di = 0; di < ndeg; ++di) {
      precomputedKeep[{static_cast<int>(triId), static_cast<int>(base + di)}] =
          false;
    }
  }
  // findPolyHE + isDegenerate helpers.
  auto findPolyHE = [&](int triId, int v0, int v1) -> int {
    const auto& polys = polygonWalks[triId].polygons;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
      const auto& poly = polys[pi];
      for (size_t i = 0; i < poly.size(); ++i) {
        if (poly[i] == v0 && poly[(i + 1) % poly.size()] == v1)
          return static_cast<int>(pi);
      }
    }
    const auto& deg = polygonWalks[triId].degeneratePolygons;
    for (size_t di = 0; di < deg.size(); ++di) {
      const auto& dp = deg[di];
      if (dp.size() != 2) continue;
      if ((dp[0] == v0 && dp[1] == v1) || (dp[0] == v1 && dp[1] == v0)) {
        return static_cast<int>(polys.size() + di);
      }
    }
    return -1;
  };
  auto isDegenerate = [&](int triId, int pi) -> bool {
    if (pi < 0) return false;
    const size_t base = polygonWalks[triId].polygons.size();
    return static_cast<size_t>(pi) >= base;
  };
  // Pair-sym Phase 1: per-direction chord-pair enforcement (6 branches).
  // Tie-breaks read from a snapshot of the classifier output taken before
  // the loop, so Phase 1's tie-break decisions do not depend on the order
  // chordOutput.newEdges happens to be in (no iteration sees a flag a
  // previous iteration already flipped). Writes still update the live map
  // that Phases 2-3.5 read, so the final kept set is not order-free - only
  // this phase's reads are.
  const std::map<std::pair<int, int>, bool> keepSnapshot = precomputedKeep;
  auto getSnap = [&](int triId, int pi) -> bool {
    if (pi < 0) return false;
    auto it = keepSnapshot.find({triId, pi});
    return it != keepSnapshot.end() && it->second;
  };
  for (const auto& edge : chordOutput.newEdges) {
    int piA1 = findPolyHE(edge.triA, edge.v0, edge.v1);
    int piB1 = findPolyHE(edge.triB, edge.v1, edge.v0);
    int piA2 = findPolyHE(edge.triA, edge.v1, edge.v0);
    int piB2 = findPolyHE(edge.triB, edge.v0, edge.v1);
    auto getKeep = [&](int triId, int pi) -> bool* {
      if (pi < 0) return nullptr;
      auto it = precomputedKeep.find({triId, pi});
      return it == precomputedKeep.end() ? nullptr : &it->second;
    };
    bool* kA1 = getKeep(edge.triA, piA1);
    bool* kB1 = getKeep(edge.triB, piB1);
    bool* kA2 = getKeep(edge.triA, piA2);
    bool* kB2 = getKeep(edge.triB, piB2);
    const bool snapA1 = getSnap(edge.triA, piA1);
    const bool snapB1 = getSnap(edge.triB, piB1);
    const bool snapA2 = getSnap(edge.triA, piA2);
    const bool snapB2 = getSnap(edge.triB, piB2);
    const bool A1d = isDegenerate(edge.triA, piA1);
    const bool B1d = isDegenerate(edge.triB, piB1);
    const bool A2d = isDegenerate(edge.triA, piA2);
    const bool B2d = isDegenerate(edge.triB, piB2);
    auto isContrib = [](bool* k, bool degen) { return k != nullptr && !degen; };
    const bool A1c = isContrib(kA1, A1d);
    const bool B2c = isContrib(kB2, B2d);
    const bool A2c = isContrib(kA2, A2d);
    const bool B1c = isContrib(kB1, B1d);
    const int n_d1 = (A1c ? 1 : 0) + (B2c ? 1 : 0);
    const int n_d2 = (A2c ? 1 : 0) + (B1c ? 1 : 0);
    if (n_d1 == 0 && n_d2 == 0) {
      // No contributors - nothing to enforce.
    } else if (n_d1 == 0) {
      if (A2c) *kA2 = false;
      if (B1c) *kB1 = false;
    } else if (n_d2 == 0) {
      if (A1c) *kA1 = false;
      if (B2c) *kB2 = false;
    } else if (n_d1 == 1 && n_d2 == 1) {
      if (A1c) *kA1 = true;
      if (B2c) *kB2 = true;
      if (A2c) *kA2 = true;
      if (B1c) *kB1 = true;
    } else if (n_d1 == 1 && n_d2 == 2) {
      if (A1c) *kA1 = true;
      if (B2c) *kB2 = true;
      if (A2c && B1c) {
        if (snapA2 || !snapB1) {
          *kA2 = true;
          *kB1 = false;
        } else {
          *kA2 = false;
          *kB1 = true;
        }
      }
    } else if (n_d1 == 2 && n_d2 == 1) {
      if (A2c) *kA2 = true;
      if (B1c) *kB1 = true;
      if (A1c && B2c) {
        if (snapA1 || !snapB2) {
          *kA1 = true;
          *kB2 = false;
        } else {
          *kA1 = false;
          *kB2 = true;
        }
      }
    } else {
      // Full case (n_d1=2, n_d2=2): pair-level enforcement.
      int p1Vote = (snapA1 ? 1 : 0) + (snapB1 ? 1 : 0);
      int p2Vote = (snapA2 ? 1 : 0) + (snapB2 ? 1 : 0);
      bool p1Keep, p2Keep;
      if (p1Vote == 2 && p2Vote == 0) {
        p1Keep = true;
        p2Keep = false;
      } else if (p1Vote == 0 && p2Vote == 2) {
        p1Keep = false;
        p2Keep = true;
      } else if (p1Vote >= p2Vote) {
        p1Keep = true;
        p2Keep = false;
      } else {
        p1Keep = false;
        p2Keep = true;
      }
      *kA1 = p1Keep;
      *kB1 = p1Keep;
      *kA2 = p2Keep;
      *kB2 = p2Keep;
    }
  }
  // Pair-sym Phase 2 + 2.5: mesh-edge consistency + multi-owner enforcement.
  {
    std::set<std::tuple<int, int, int>> chordHalfedges;
    for (const auto& edge : chordOutput.newEdges) {
      chordHalfedges.insert({edge.triA, edge.v0, edge.v1});
      chordHalfedges.insert({edge.triA, edge.v1, edge.v0});
      chordHalfedges.insert({edge.triB, edge.v0, edge.v1});
      chordHalfedges.insert({edge.triB, edge.v1, edge.v0});
    }
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> halfedgeMap;
    for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
      const auto& polys = polygonWalks[triId].polygons;
      for (size_t pi = 0; pi < polys.size(); ++pi) {
        const auto& poly = polys[pi];
        for (size_t i = 0; i < poly.size(); ++i) {
          int a = poly[i], b = poly[(i + 1) % poly.size()];
          halfedgeMap[{a, b}].push_back(
              {static_cast<int>(triId), static_cast<int>(pi)});
        }
      }
    }
    // Phase 2: AND-merge single-owner-per-direction non-chord pairs.
    for (const auto& [dirEdge, owners] : halfedgeMap) {
      if (owners.size() != 1) continue;
      int a = dirEdge.first, b = dirEdge.second;
      auto revIt = halfedgeMap.find({b, a});
      if (revIt == halfedgeMap.end() || revIt->second.size() != 1) continue;
      const auto& [triA, piA] = owners[0];
      const auto& [triB, piB] = revIt->second[0];
      if (triA == triB) continue;
      if (chordHalfedges.count({triA, a, b}) > 0) continue;
      auto kA_it = precomputedKeep.find({triA, piA});
      auto kB_it = precomputedKeep.find({triB, piB});
      if (kA_it == precomputedKeep.end() || kB_it == precomputedKeep.end())
        continue;
      if (kA_it->second != kB_it->second) {
        bool target = kA_it->second && kB_it->second;
        kA_it->second = target;
        kB_it->second = target;
      }
    }
    // Phase 2.5: multi-owner per-direction enforcement (non-chord).
    for (const auto& [dirEdge, owners] : halfedgeMap) {
      if (owners.size() < 2) continue;
      int a = dirEdge.first, b = dirEdge.second;
      bool isChord = false;
      for (const auto& [t, p] : owners) {
        if (chordHalfedges.count({t, a, b}) > 0) {
          isChord = true;
          break;
        }
      }
      if (isChord) continue;
      std::vector<std::pair<int, int>> kept;
      for (const auto& o : owners) {
        auto it = precomputedKeep.find(o);
        if (it == precomputedKeep.end()) continue;
        if (it->second) kept.push_back(o);
      }
      if (kept.size() <= 1) continue;
      int keepIdx = 0;
      for (size_t i = 1; i < kept.size(); ++i) {
        if (kept[i].first < kept[keepIdx].first) keepIdx = static_cast<int>(i);
      }
      for (size_t i = 0; i < kept.size(); ++i) {
        if (static_cast<int>(i) == keepIdx) continue;
        precomputedKeep[kept[i]] = false;
      }
    }
  }

  // Pair-sym Phase 3: kept-polygon manifold-pair enforcement.
  //
  // Phases 1, 2, 2.5 enforce per-chord and per-edge constraints, but
  // they do not VERIFY that the resulting kept-polygon set forms a
  // closed (manifold) surface. The triangulator faithfully turns each
  // kept polygon into sub-tris, and if the kept polygons leave any
  // halfedge unpaired, the assembled MeshGL64 is non-manifold and
  // Manifold construction rejects it (= cray, gt-7081 today).
  //
  // Iteratively drop kept polygons that contribute unpaired halfedges.
  // Drop is monotonic so the loop converges. Other fixtures already
  // have near-zero unpaired counts (offset12=0, havocglass=0,
  // gt-7863=21/2439, hull-mask=305/38101, self-intersect=293/92913)
  // and the iteration only drops a handful. cray's 69/74 unpaired-
  // to-kept ratio means most polygons get dropped, falling back to
  // input via the post-pipeline pierce/drift gate; the pipeline at
  // least no longer hands a non-manifold MeshGL64 to construction.
  {
    // Net-progress criterion: dropping a polygon with k unpaired
    // halfedges out of size total halfedges removes k unpaired and
    // adds (size - k) new unpaired (= partners of the formerly-paired
    // edges become orphans). Net change in global unpaired count is
    // k - (size - k) = 2k - size, which is < 0 (= progress) iff
    // k > size/2. Strict majority avoids the no-progress case where
    // k == size/2 (= net zero, infinite oscillation risk).
    int iter = 0, totalDropped = 0;
    while (iter < kPairSymPhase3MaxIter) {
      ++iter;
      std::map<std::pair<int, int>, int> dirCnt;
      for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
        const auto& polys = polygonWalks[triId].polygons;
        for (size_t pi = 0; pi < polys.size(); ++pi) {
          auto it = precomputedKeep.find(
              {static_cast<int>(triId), static_cast<int>(pi)});
          bool keep = (it == precomputedKeep.end()) ? true : it->second;
          if (!keep) continue;
          const auto& poly = polys[pi];
          for (size_t i = 0; i < poly.size(); ++i) {
            int a = poly[i], b = poly[(i + 1) % poly.size()];
            ++dirCnt[{a, b}];
          }
        }
      }
      int dropped = 0;
      for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
        const auto& polys = polygonWalks[triId].polygons;
        for (size_t pi = 0; pi < polys.size(); ++pi) {
          auto it = precomputedKeep.find(
              {static_cast<int>(triId), static_cast<int>(pi)});
          if (it == precomputedKeep.end() || !it->second) continue;
          const auto& poly = polys[pi];
          int unpaired = 0;
          for (size_t i = 0; i < poly.size(); ++i) {
            int a = poly[i], b = poly[(i + 1) % poly.size()];
            if (dirCnt[{b, a}] == 0) ++unpaired;
          }
          if (2 * unpaired > static_cast<int>(poly.size())) {
            it->second = false;
            ++dropped;
          }
        }
      }
      totalDropped += dropped;
      if (dropped == 0) break;
    }
    DEBUG_ASSERT(iter < kPairSymPhase3MaxIter, logicErr,
                 "pair-sym Phase 3 hit kPairSymPhase3MaxIter without "
                 "converging - drop loop should be monotonic");
    // Phase 3.5: re-key pass. Phase 3 (above) DROPS kept polygons
    // with strict-majority unpaired halfedges. The symmetric problem
    // is dropped polygons whose halfedges would COMPLETE existing
    // unpaired (b,a) entries in the kept set if they were re-keyed.
    //
    // For dropped polygon P with halfedges H, partition H by status
    // relative to the current kept set's dirCnt:
    //   - paired_to_kept: H's reverse already exists in kept set
    //     (= re-keying P closes an unpaired pair, -1 to unpaired)
    //   - novel: H is not in kept set and reverse is not in kept set
    //     (= re-keying P adds a fresh unpaired, +1 to unpaired)
    //   - duplicate: H itself already exists in kept set (= re-keying
    //     P creates a multi-owner halfedge, harmful for manifoldness)
    //
    // Re-key P when paired_to_kept > novel AND duplicate == 0. The
    // duplicate==0 guard prevents introducing k>2 edges that the
    // recovery sweep would then have to clean. Iterate to fixed
    // point. For gt-7081 this trims 1539 unpaired halfedges down
    // toward 0; for working fixtures it makes a smaller difference
    // (191 -> ~ for self-intersect).
    {
      int iter2 = 0, totalRekeyed = 0;
      while (iter2 < kPairSymPhase35MaxIter) {
        ++iter2;
        std::map<std::pair<int, int>, int> dirCnt;
        for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
          const auto& polys = polygonWalks[triId].polygons;
          for (size_t pi = 0; pi < polys.size(); ++pi) {
            auto it = precomputedKeep.find(
                {static_cast<int>(triId), static_cast<int>(pi)});
            bool keep = (it == precomputedKeep.end()) ? true : it->second;
            if (!keep) continue;
            const auto& poly = polys[pi];
            for (size_t i = 0; i < poly.size(); ++i)
              ++dirCnt[{poly[i], poly[(i + 1) % poly.size()]}];
          }
        }
        int rekeyed = 0;
        for (size_t triId = 0; triId < polygonWalks.size(); ++triId) {
          const auto& polys = polygonWalks[triId].polygons;
          for (size_t pi = 0; pi < polys.size(); ++pi) {
            auto it = precomputedKeep.find(
                {static_cast<int>(triId), static_cast<int>(pi)});
            if (it == precomputedKeep.end() || it->second) continue;
            const auto& poly = polys[pi];
            int pairedToKept = 0, novel = 0, duplicate = 0;
            for (size_t i = 0; i < poly.size(); ++i) {
              int a = poly[i], b = poly[(i + 1) % poly.size()];
              if (dirCnt.count({a, b}) > 0)
                ++duplicate;
              else if (dirCnt.count({b, a}) > 0)
                ++pairedToKept;
              else
                ++novel;
            }
            // Net change in unpaired count from re-keying P:
            //   novel - pairedToKept (new unpaired added minus existing
            //   unpaired closed). Re-key only when it strictly reduces
            //   unpaired (pairedToKept > novel) AND introduces no duplicate
            //   (multi-owner) halfedge, so we never manufacture k>2 edges
            //   the downstream reducer would then have to clean.
            if (pairedToKept > novel && duplicate == 0) {
              it->second = true;
              ++rekeyed;
            }
          }
        }
        totalRekeyed += rekeyed;
        if (rekeyed == 0) break;
      }
      DEBUG_ASSERT(iter2 < kPairSymPhase35MaxIter, logicErr,
                   "pair-sym Phase 3.5 hit kPairSymPhase35MaxIter without "
                   "converging - re-key loop should be monotonic");
    }
  }
  // Build the classifier closure that returns precomputedKeep[triId, pi]
  // (or AnalyticalKeep fallback if pi can't be located).
  PolygonClassifierFn classifier = [&](const std::vector<int>& poly,
                                       const vec3& triNormal,
                                       int triId) -> PolygonClassification {
    PolygonClassification c{false, 0, 0};
    if (poly.size() < 3) return c;
    int pi = -1;
    const auto& polys = polygonWalks[triId].polygons;
    for (size_t i = 0; i < polys.size(); ++i) {
      if (polys[i].size() == poly.size()) {
        bool match = true;
        for (size_t j = 0; j < poly.size(); ++j)
          if (polys[i][j] != poly[j]) {
            match = false;
            break;
          }
        if (match) {
          pi = static_cast<int>(i);
          break;
        }
      }
    }
    if (pi >= 0) {
      auto it = precomputedKeep.find({triId, pi});
      c.keep = (it != precomputedKeep.end()) ? it->second : true;
    } else {
      c.keep = AnalyticalKeep(triId, poly, chordPartners,
                              chordOutput.newVertPositions, baseId, impl);
    }
    (void)triNormal;
    return c;
  };
  auto step13 = TriangulateAndEmit(impl, chordOutput.newVertPositions,
                                   polygonWalks, classifier);
  Manifold out = step13.output;
  // Fallback chooser shared by every gate-reject path. Returns the merged
  // input only when it is itself a valid manifold with no more pierces
  // than the original (MergeVertsEps can both produce a non-manifold and
  // introduce pierces by eps-merging verts onto other faces); otherwise
  // the original input, which is always valid and pierce-monotonicity-safe.
  const auto pickFallback = [&]() -> Manifold {
    if (mr.manifold.Status() == Manifold::Error::NoError &&
        CheckSelfIntersection(mr.manifold, kPipelineRelTol).interiorPierces <=
            CheckSelfIntersection(input, kPipelineRelTol).interiorPierces) {
      return mr.manifold;
    }
    return input;
  };
  if (out.Status() != Manifold::Error::NoError) {
    // Pipeline output non-manifold -> fall back to the safe input.
    return pickFallback();
  }
  // Pierce-monotonicity guard: pipeline must not produce more pierces
  // than the input has. Volume-sanity: drift > kDriftCutoff or sign
  // flip also triggers fallback. This is a user-visible API contract,
  // not optional.
  //
  // Both gates anchor to the ORIGINAL input, not the merged
  // mr.manifold: MergeVertsEps can introduce pierces (by eps-merging
  // close verts onto edges of other faces) and shift volume, so
  // post-merge counts may already differ from input. The user-visible
  // regression is input vs out - that's what we gate (and the
  // sign-flip baseline below) against.
  const double inVol = input.Volume();
  const double outVol = out.Volume();
  auto siInput = CheckSelfIntersection(input, kPipelineRelTol);
  auto siOut = CheckSelfIntersection(out, kPipelineRelTol);
  const bool pierceWorse = siOut.interiorPierces > siInput.interiorPierces;
  const double drift =
      std::fabs(outVol - inVol) / std::max(std::fabs(inVol), kVolumeFloor);
  const bool signFlip = (inVol > 0) != (outVol > 0);
  const bool driftBad = (drift > kDriftCutoff) || signFlip;
  // Sign-flip recovery: when the output is sign-flipped (regardless of
  // drift), try reversing winding before falling back. Subtract-derived
  // inputs with back-side normal flips can produce inverted manifolds
  // whose mirror is the right answer.
  if (signFlip && !pierceWorse) {
    MeshGL64 flipped = step13.output.GetMeshGL64();
    for (size_t t = 0; t < flipped.triVerts.size() / 3; ++t) {
      std::swap(flipped.triVerts[3 * t + 1], flipped.triVerts[3 * t + 2]);
    }
    Manifold flipMan(flipped);
    if (flipMan.Status() == Manifold::Error::NoError) {
      const double flipVol = flipMan.Volume();
      const double flipDrift =
          std::fabs(flipVol - inVol) / std::max(std::fabs(inVol), kVolumeFloor);
      const bool flipSignOK = (inVol > 0) == (flipVol > 0);
      auto siFlip = CheckSelfIntersection(flipMan, kPipelineRelTol);
      if (flipSignOK && flipDrift <= kDriftCutoff &&
          siFlip.interiorPierces <= siInput.interiorPierces) {
        return flipMan;
      }
    }
  }
  if (pierceWorse || driftBad) {
    return pickFallback();
  }
  return out;
}
}  // namespace

EdgeReducerResult TrimOrphans(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples,
    int maxRounds) {
  // No DEBUG_ASSERT on cap-hit: this loop alternates trim and cap-
  // close, which may oscillate (cap can add tris back). The cap is
  // a per-call budget, not a convergence tripwire.
  EdgeReducerResult r;
  for (int trim = 0; trim < maxRounds; ++trim) {
    const size_t nT = out.triVerts.size() / 3;
    if (nT == 0) break;
    auto ecT = ComputeEdgeIncidence(out);
    int curK1 = 0;
    for (const auto& [_, c] : ecT)
      if (c == 1) ++curK1;
    if (curK1 == 0) break;
    std::vector<bool> dropped(nT, false);
    int trimmed = 0;
    // Two-tier: round 0 drops tris with >=2 k=1 edges; subsequent
    // rounds drop tris with >=1 k=1 edge.
    const int trimThresh = (trim == 0) ? 2 : 1;
    for (size_t t = 0; t < nT; ++t) {
      int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                  static_cast<int>(out.triVerts[3 * t + 1]),
                  static_cast<int>(out.triVerts[3 * t + 2])};
      int k1Cnt = 0;
      for (int e : {0, 1, 2}) {
        int a = v[e], b = v[(e + 1) % 3];
        if (a > b) std::swap(a, b);
        auto it = ecT.find({a, b});
        if (it != ecT.end() && it->second == 1) ++k1Cnt;
      }
      if (k1Cnt >= trimThresh) {
        dropped[t] = true;
        ++trimmed;
      }
    }
    if (trimmed == 0) break;
    int actualDropped = CompactDroppedTris(out, dropped);
    r.totalDropped += actualDropped;
    auto [c, ct] = DoCapPass(out, forbiddenTriples);
    r.totalCapped += c;
    r.totalCapTris += ct;
  }
  return r;
}

}  // namespace overlap_removal
}  // namespace manifold
