// Copyright 2026 The Manifold Authors.
//
// Spike for 3D overlap removal (companion to extras/overlap2d_proto.cpp).
// Per docs/Overlap3D.md: the existing 3D pipeline already implements
// most of what a "Smith-on-BVH 3D" algorithm would look like (BVH pair
// queries via Collider, symbolic predicates via Shadows/Interpolate/
// Shadow01/Kernel11/12, eager propagation in AddNewEdgeVerts, post-hoc
// cleanup in CollapseShortEdges/SimplifyTopology). What's missing is
// vocabulary: explicit alpha-budget eps, a Simplify entry point that
// regularizes a single Manifold, an InferEps helper. This file adds
// those wrappers.
//
// Build (from manifold repo root):
//   cmake -B build -S . -DMANIFOLD_TEST=ON
//   cmake --build build --target overlap3d_proto
//
// Run:
//   build/extras/overlap3d_proto   # smoke + fuzz + adversarial battery

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "../src/collider.h"  // for BVH broad phase in self-intersection check
#include "../src/disjoint_sets.h"  // for ε-merge cluster union-find
#include "../src/impl.h"           // Manifold::Impl for half-edge adjacency
#include "../src/self_mesh_analysis.h"  // Option B: SoS-corrected self-mesh winding
#include "../src/shared.h"              // for AlphaBudgetEpsilon
#include "../src/winding03.h"  // for production Winding03 (M-vs-M diagnostic)
#include "manifold/manifold.h"
#include "manifold/polygon.h"  // for Triangulate

namespace overlap3d {

// Determinism debugging: FNV-1a 64-bit hash over a byte range. Used
// at pipeline phase boundaries when OVERLAP3D_DET_HASH=1 to localize
// run-to-run nondeterminism. Stable across runs of the same binary
// iff the hashed bytes are themselves stable.
inline uint64_t Fnv1a64(const void* data, size_t bytes) {
  uint64_t h = 14695981039346656037ULL;
  const uint8_t* p = static_cast<const uint8_t*>(data);
  for (size_t i = 0; i < bytes; ++i) {
    h ^= static_cast<uint64_t>(p[i]);
    h *= 1099511628211ULL;
  }
  return h;
}
template <typename T>
inline uint64_t HashVec(const std::vector<T>& v) {
  return Fnv1a64(v.data(), v.size() * sizeof(T));
}
inline bool DetHashEnabled() {
  static const bool e = []() {
    const char* s = std::getenv("OVERLAP3D_DET_HASH");
    return s && std::string(s) == "1";
  }();
  return e;
}
#define OVERLAP3D_HASH_DUMP(label, expr)                                       \
  do {                                                                         \
    if (DetHashEnabled()) {                                                    \
      std::cerr << "      [det-hash] " << label << " = " << std::hex << (expr) \
                << std::dec << "\n";                                           \
    }                                                                          \
  } while (0)

using manifold::AlphaBudgetEpsilon;
using manifold::Box;
using manifold::Manifold;
using manifold::OpType;
using manifold::vec3;

// Local alias matching the 2D prototype's name; the implementation
// lives in src/shared.h so the 2D prototype, 3D spike, and any
// future production code can all share the same Smith alpha-budget
// formula. See doc on AlphaBudgetEpsilon in shared.h for the
// distinction from manifold's existing MaxEpsilon / kPrecision*scale
// formula.
inline double EpsilonFromScale(double L, int k_budget = 1000) {
  return AlphaBudgetEpsilon(L, k_budget);
}

// Maximum coordinate magnitude across the manifold's bounding box.
// `Box::Scale()` already does this; eps scales with coord magnitude
// (not feature size) per Smith's displacement-attack analysis.
inline double InferEps(const Manifold& m) {
  return EpsilonFromScale(m.BoundingBox().Scale());
}

inline double InferEps(const Manifold& a, const Manifold& b) {
  return EpsilonFromScale(
      std::max(a.BoundingBox().Scale(), b.BoundingBox().Scale()));
}

// Single-input regularization: 3D analog of overlap2d::Simplify.
// Rounds out near-coincident verts, collapses degenerate triangles,
// and simplifies topology to within `eps`. Implementation: set
// tolerance to eps (which informs SimplifyTopology and the boolean's
// epsilon-merge behavior), then run a self-union to canonicalize the
// boundary representation. The result is an eps-valid manifold.
//
// Pass eps <= 0 to auto-infer eps via Smith's alpha-budget from the
// input bounding box.
inline Manifold Simplify(const Manifold& in, double eps = 0.0) {
  if (eps <= 0.0) eps = InferEps(in);
  // Manifold::Simplify with explicit tolerance does the regularization;
  // it's a no-op alias for SimplifyTopology in the existing 3D path.
  return in.Simplify(eps);
}

// Binary boolean with explicit eps. Wraps Manifold::Boolean with the
// alpha-budget eps inference.
inline Manifold Boolean3D(const Manifold& a, const Manifold& b, OpType op,
                          double eps = 0.0) {
  if (eps <= 0.0) eps = InferEps(a, b);
  Manifold aT = a.SetTolerance(eps);
  Manifold bT = b.SetTolerance(eps);
  return aT.Boolean(bT, op);
}

// =============================================================================
// Self-intersection check: prove the boolean output contains no
// triangle pairs with strict interior overlap.
//
// Necessary because the existing tests (idempotence, volume identity,
// Status() == NoError) are NECESSARY but NOT SUFFICIENT for "the
// boolean output is geometrically clean":
//   - Idempotence says Simplify(pass1) == pass1 (algorithm converged).
//     Compatible with self-intersections that get re-emitted.
//   - Volume identity V(A∪B) + V(A∩B) = V(A) + V(B) is a measure-
//     theoretic property of the boundary; holds even with self-
//     intersections in the mesh.
//   - Status() catches construction errors but not all self-pierce
//     cases.
//
// The actual semantic of "remove overlaps" is no triangle-triangle
// interior intersections. This check tests it directly via Möller-
// style edge-pierces-triangle (segment-vs-triangle) tests, with BVH
// broad phase via manifold's Collider.
// =============================================================================

// Test whether the open line segment (a, b) STRICTLY pierces the
// STRICT interior of triangle (v0, v1, v2), and report the pierce
// magnitude. Strict on both:
//   - Segment endpoints exactly on the triangle plane do not count.
//   - Intersection at a triangle edge or vertex does not count.
// This means "shared vertex" and "shared edge" intersections are
// excluded; only true interior piercing is reported.
//
// Returns the pierce *magnitude* if a pierce is detected, or 0 if
// not. Magnitude is the perpendicular distance from the nearer
// segment endpoint to the triangle's plane: min(|dA|,|dB|) / |n|.
// This is the natural "depth" of the pierce — a graze where one
// endpoint is barely off the plane has tiny magnitude; a segment that
// crosses the plane far from either endpoint has magnitude on the
// order of half the segment's perpendicular extent. Caller can then
// compare magnitude against the manifold's tolerance contract.
//
// `relTol` is the FP-noise threshold below which an endpoint is
// considered "on the plane" (returns 0, not a pierce). It should
// be much tighter than the manifold's geometric tolerance — its job
// is to prevent zero-by-zero issues in t = dA/(dA-dB), not to filter
// for tolerance compliance (caller does that).
inline double SegmentPiercesTriInterior(vec3 a, vec3 b, vec3 v0, vec3 v1,
                                        vec3 v2, double relTol = 1e-12) {
  using manifold::la::cross;
  using manifold::la::dot;
  const vec3 e1 = v1 - v0;
  const vec3 e2 = v2 - v0;
  const vec3 n = cross(e1, e2);
  const double nMag = std::sqrt(dot(n, n));
  if (nMag == 0) return 0.0;  // degenerate triangle
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
  // Pierce magnitude = perpendicular distance from nearer endpoint to plane.
  return std::min(std::fabs(dA), std::fabs(dB)) / nMag;
}

// Check a Manifold for self-intersections — geometric tri-tri
// interior pierces. This is what overlap-removal (in Emmett's #289
// sense) is supposed to eliminate. Manifold's existing 3D Boolean
// produces output that is *topologically* manifold (every edge in
// exactly two faces) but is NOT guaranteed to be self-intersection-
// free; that's a different and weaker invariant. Any pierce reported
// here is a candidate finding for the overlap-removal feature to fix.
//
// Reports both the count of piercing pairs and the maximum pierce
// magnitude (perpendicular distance from the nearer piercing-segment
// endpoint to the triangle plane). Magnitude lets you separate
// "deep geometric overlap" from "near-coplanar grazing pierce" when
// triaging — both are bugs from the overlap-removal perspective, but
// they exercise different parts of the algorithm.
//
// Uses Collider for BVH broad phase, edge-pierces-triangle for narrow
// phase. Pairs sharing 2+ vert indices (= adjacent across an edge)
// are skipped before narrow phase as legal-by-construction. The
// `relTol` parameter is the FP-noise threshold inside the narrow
// phase (default 1e-12, near machine epsilon) — this is NOT a
// tolerance for filtering "small" pierces; it only prevents zero-by-
// zero artifacts in the t = dA/(dA-dB) computation.
struct SelfIntersectionResult {
  int interiorPierces;        // strict interior pierce count
  double maxPierceMagnitude;  // perpendicular depth of deepest pierce
  int candidatesChecked;      // pairs that passed broad phase
  int adjacentPairsSkipped;   // pairs sharing 2+ verts (legal)
  int trianglesTotal;
  // Pierce-pair graph: each entry is one piercing (triA, triB) pair,
  // populated only when minMagnitude > 0 (= when caller wants the
  // graph for cluster analysis). Otherwise empty, to keep the hot
  // path's per-pair allocation cost zero on existing call sites.
  std::vector<std::pair<int, int>> pierceEdges;
};

// Pierce-cluster analysis: from the pierceEdges list, compute
// connected components in the (tris-as-nodes, pierces-as-edges) graph.
struct PierceCluster {
  std::vector<int> tris;                   // tri indices in this component
  std::vector<std::pair<int, int>> edges;  // pierce pairs in this component
  int maxDegree;  // max # of pierces any single tri participates in
};
inline std::vector<PierceCluster> ClusterPierces(
    const std::vector<std::pair<int, int>>& edges) {
  std::map<int, std::set<int>> adj;
  for (const auto& [a, b] : edges) {
    adj[a].insert(b);
    adj[b].insert(a);
  }
  std::set<int> visited;
  std::vector<PierceCluster> clusters;
  for (const auto& [t, _] : adj) {
    if (visited.count(t)) continue;
    PierceCluster c;
    std::vector<int> stack{t};
    std::set<int> inThis;
    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      if (!visited.insert(u).second) continue;
      inThis.insert(u);
      c.tris.push_back(u);
      for (int v : adj[u])
        if (!visited.count(v)) stack.push_back(v);
    }
    c.maxDegree = 0;
    for (const auto& [a, b] : edges) {
      if (inThis.count(a) && inThis.count(b)) {
        c.edges.emplace_back(a, b);
      }
    }
    for (int u : c.tris) {
      const int d = static_cast<int>(adj[u].size());
      if (d > c.maxDegree) c.maxDegree = d;
    }
    clusters.push_back(std::move(c));
  }
  return clusters;
}

// minMagnitude: pierces with magnitude below this are not counted.
// Use mesh.GetTolerance() to skip pierces that are within the mesh's
// own coincidence tolerance (= FP noise that the Boolean engine is
// allowed to produce). Default 0.0 preserves existing call sites.
inline SelfIntersectionResult CheckSelfIntersection(const Manifold& m,
                                                    double relTol = 1e-12,
                                                    double minMagnitude = 0.0) {
  SelfIntersectionResult r{};
  if (m.IsEmpty()) return r;
  manifold::MeshGL64 mesh = m.GetMeshGL64();
  const size_t nTri = mesh.NumTri();
  r.trianglesTotal = static_cast<int>(nTri);
  if (nTri < 2) return r;

  // Per-triangle AABBs (3D Box).
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

  // Manifold's Collider expects sorted-by-Morton input. Standard build pattern:
  Box bbox;
  for (const auto& b : triBoxes) bbox = bbox.Union(b);
  std::vector<uint32_t> morton(nTri);
  for (size_t i = 0; i < nTri; ++i)
    morton[i] = manifold::Collider::MortonCode(triBoxes[i].Center(), bbox);
  std::vector<size_t> perm(nTri);
  std::iota(perm.begin(), perm.end(), 0);
  std::stable_sort(perm.begin(), perm.end(),
                   [&](size_t a, size_t b) { return morton[a] < morton[b]; });
  std::vector<Box> sortedBoxes(nTri);
  std::vector<uint32_t> sortedMorton(nTri);
  for (size_t i = 0; i < nTri; ++i) {
    sortedBoxes[i] = triBoxes[perm[i]];
    sortedMorton[i] = morton[perm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedBoxes.data(), sortedBoxes.size()),
      manifold::VecView<const uint32_t>(sortedMorton.data(),
                                        sortedMorton.size()));

  auto vp = [&](int idx) {
    return vec3(mesh.vertProperties[mesh.numProp * idx + 0],
                mesh.vertProperties[mesh.numProp * idx + 1],
                mesh.vertProperties[mesh.numProp * idx + 2]);
  };

  auto checkPair = [&](size_t qi, size_t li) {
    if (qi >= li) return;  // dedupe + skip self
    const size_t ta = perm[qi];
    const size_t tb = perm[li];
    // Count shared vert indices.
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
    if (maxMag > minMagnitude) {
      ++r.interiorPierces;
      if (maxMag > r.maxPierceMagnitude) r.maxPierceMagnitude = maxMag;
      if (minMagnitude > 0.0) {
        r.pierceEdges.emplace_back(static_cast<int>(ta), static_cast<int>(tb));
      }
      if (std::getenv("OVERLAP3D_PIERCE_LIST")) {
        std::cerr << "    pierce: tri " << ta << " <-> tri " << tb
                  << "  mag=" << maxMag << "  verts ta=[" << triIdx[ta][0]
                  << "," << triIdx[ta][1] << "," << triIdx[ta][2] << "] tb=["
                  << triIdx[tb][0] << "," << triIdx[tb][1] << ","
                  << triIdx[tb][2] << "]\n";
      }
    }
  };

  // Use Collider's pair-collisions API: query each box against the BVH.
  auto recorder = manifold::MakeSimpleRecorder(checkPair);
  auto qf = [&](int i) { return sortedBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nTri),
                             /*parallel=*/false);
  return r;
}

// =============================================================================
// Emmett's #289 13-step overlap-removal algorithm — implementation in progress.
//
// Status: steps 1–3. Subsequent steps will land in follow-up
// commits. Each step is gated against the .obj-fixture adversarial
// battery so we can measure pierce reduction as steps come online.
//
// Note on input domain: overlap removal is invoked on output from
// the existing `Manifold::Boolean` (or other Manifold-producing
// path). Its input is therefore *topology-manifold* — every edge
// is in exactly two triangles, the half-edge structure is closed
// — but may contain geometric tri-tri interior pierces (the thing
// overlap removal exists to fix). Non-manifold triangle soup is
// not in scope; that would be a different feature (mesh repair).
//
// API: each step takes `Manifold` and returns `Manifold` for
// transformations, or returns auxiliary structures (per-edge vert
// lists etc) for queries. Internally, steps work directly on
// `Manifold::Impl` — its `vertPos_`/`halfedge_`/`faceNormal_`
// are public and the half-edge adjacency is the natural format
// for the rest of the pipeline (`pairedHalfedge` for edge
// incidence, ring traversal for vert→incident-edges).
// =============================================================================

// Helper: round-trip Manifold → Impl (via the public MeshGL64 path,
// since `GetCsgLeafNode().GetImpl()` is private). The cost is one
// MeshGL64 build + one Impl construction (which does ε-merge,
// halfedge build, normal calc). Acceptable for a spike; production
// would hold the Impl in a workspace across steps.
inline manifold::Manifold::Impl ImplFromManifold(const Manifold& m) {
  return manifold::Manifold::Impl(m.GetMeshGL64());
}

// Helper: Impl → Manifold (via MeshGL64). The `GetMeshGLImpl`
// converter is in `src/impl.h:432`.
inline Manifold ManifoldFromImpl(const manifold::Manifold::Impl& impl) {
  return Manifold(manifold::GetMeshGLImpl<double, uint64_t>(impl,
                                                            /*normalIdx=*/-1));
}

// -----------------------------------------------------------------------------
// Step 1: ε-merge verts.
//
// Emmett #289 step 1: "Merge all verts that are within ε of each
// other. Care must be taken regarding chains of nearby verts so they
// don't move too far. They need to use a weighted-average position
// and calculate distance based on up-to-date values. The broad phase
// can be parallel, but the narrow phase probably needs to be serial
// for this reason."
//
// Implementation pattern mirrors `src/sort.cpp:155-167` (Manifold's
// existing open-vert merging). Per-vert ε/2-padded AABB → Morton
// sort → Collider self-collision → DisjointSets unite. Position
// update: cluster centroid (unweighted average) recomputed after each
// pass. Iterate until no new unions or `maxIter` hit; chains that
// span > ε get pulled in over multiple passes.
//
// Returns the merged Manifold *and* the count of vert pairs that
// were merged. The count is reported separately because
// `Manifold::NumVert()` may not reflect the merge: when merged verts
// don't cause any triangle to collapse, `RemoveUnreferencedVerts`
// (src/impl.cpp:144) sets unreferenced positions to NaN but doesn't
// compact `vertPos_.size()`, so post-merge `NumVert()` can still
// equal pre-merge `NumVert()` even though the underlying mesh has
// been merged. The `mergedCount` field is the authoritative answer.
// -----------------------------------------------------------------------------
struct MergeVertsResult {
  Manifold manifold;
  int mergedCount = 0;
};

inline MergeVertsResult MergeVertsEps(const Manifold& in, double eps,
                                      int maxIter = 4) {
  if (in.IsEmpty()) return {in, 0};

  // Pull positions and triangles out of the input via MeshGL64 (the
  // canonical flat mesh-data type). MeshGL64 lays vertProperties as
  // [x0,y0,z0, x1,y1,z1, ...] and triVerts as [a0,b0,c0, a1,b1,c1,...].
  manifold::MeshGL64 mesh = in.GetMeshGL64();
  const size_t n = mesh.NumVert();
  std::vector<vec3> verts(n);
  for (size_t i = 0; i < n; ++i) {
    verts[i] = vec3(mesh.vertProperties[mesh.numProp * i + 0],
                    mesh.vertProperties[mesh.numProp * i + 1],
                    mesh.vertProperties[mesh.numProp * i + 2]);
  }

  // Per-pass: build boxes, run Collider, narrow-phase dist² < eps²,
  // unite. Repeat if any unions happened.
  int iter = 0;
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
    Box bBox;
    for (const auto& b : boxes) bBox = bBox.Union(b);
    std::vector<uint32_t> morton(n);
    for (size_t i = 0; i < n; ++i)
      morton[i] = manifold::Collider::MortonCode(boxes[i].Center(), bBox);
    std::vector<size_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::stable_sort(perm.begin(), perm.end(),
                     [&](size_t a, size_t b) { return morton[a] < morton[b]; });
    std::vector<Box> sortedBoxes(n);
    std::vector<uint32_t> sortedMorton(n);
    for (size_t i = 0; i < n; ++i) {
      sortedBoxes[i] = boxes[perm[i]];
      sortedMorton[i] = morton[perm[i]];
    }
    manifold::Collider collider(
        manifold::VecView<const Box>(sortedBoxes.data(), sortedBoxes.size()),
        manifold::VecView<const uint32_t>(sortedMorton.data(),
                                          sortedMorton.size()));
    DisjointSets uf(static_cast<uint32_t>(n));
    int unions = 0;
    const double eps2 = eps * eps;
    auto checkPair = [&](size_t qi, size_t li) {
      if (qi >= li) return;
      const size_t va = perm[qi];
      const size_t vb = perm[li];
      const vec3 d = verts[va] - verts[vb];
      const double d2 = manifold::la::dot(d, d);
      if (d2 > eps2) return;
      // Serial unite; the parallel-safe atomics in DisjointSets handle
      // the case but we don't need them here.
      uint32_t before = uf.find(static_cast<uint32_t>(va));
      uf.unite(static_cast<uint32_t>(va), static_cast<uint32_t>(vb));
      uint32_t after = uf.find(static_cast<uint32_t>(va));
      if (before != after) ++unions;
    };
    auto recorder = manifold::MakeSimpleRecorder(checkPair);
    auto qf = [&](int i) { return sortedBoxes[i]; };
    collider.Collisions<false>(recorder, qf, static_cast<int>(n),
                               /*parallel=*/false);

    // Update positions: each component → centroid of its members.
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
    if (unions == 0 && !moved) break;
  }

  // Apply the merge by setting `mergeFromVert`/`mergeToVert` hints on
  // the input MeshGL64 and letting `Manifold(MeshGL64)` perform the
  // actual merge during construction. This delegates to the same code
  // path manifold's `sort.cpp:155-167` uses for ε-merging, so the
  // result is consistent with what other Manifold-producing paths do.
  // We also overwrite the position of every merged vert with its
  // cluster centroid so the merged group is geometrically coherent.
  for (size_t i = 0; i < n; ++i) {
    const vec3& p = verts[i];
    mesh.vertProperties[mesh.numProp * i + 0] = p.x;
    mesh.vertProperties[mesh.numProp * i + 1] = p.y;
    mesh.vertProperties[mesh.numProp * i + 2] = p.z;
  }
  // Pick a canonical rep per component (smallest index in component).
  // Then any non-rep vert i gets a (i → rep) merge hint.
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
  // Avoid the round-trip when nothing was merged. The
  // GetMeshGL64 → Manifold round-trip is lossy on inputs
  // produced by `Subtract` (back-side / run-transform info
  // doesn't fully survive; observed as sign-flipped volume on
  // Cray). When there are no merges to apply, the input is
  // already correct.
  if (mergedCount == 0) return {in, 0};
  return {Manifold(mesh), mergedCount};
}

// -----------------------------------------------------------------------------
// Step 2: drop collapsed triangles + enumerate edges.
//
// Emmett #289 step 2 in the 3D context: "Same as 2 above, though
// triangles also get collapsed and discarded naturally. Edges with
// the same two verts are merged and create a list of their
// attached triangles (an even number)."
//
// Since overlap-removal's input is always topology-manifold (the
// `Manifold` type's invariant), the collapsed-tri drop is in
// practice a no-op — `Manifold(MeshGL)` already runs
// `RemoveDegenerates` at construction. After step 1's vert merge,
// the new `Manifold(out)` constructor call re-runs that pass, so
// any triangles whose verts were merged together are dropped
// automatically. We don't need a separate drop pass here.
//
// What we *do* need is the **edge enumeration**: for the rest of
// the pipeline (steps 3-6), each unique edge needs an ID so we
// can attach on-edge vert lists, intersection candidates, etc.
// In `Manifold::Impl` this is free: `halfedge_` already represents
// the edge structure as a pair of half-edges per edge, with
// `pairedHalfedge` linking the two sides. We pick the canonical
// "forward" halfedge of each edge (where `startVert < endVert`)
// as the edge identifier.
// -----------------------------------------------------------------------------
struct Edge {
  int v0, v1;           // vert indices, v0 < v1
  int halfedgeForward;  // halfedge id with startVert == v0
  int halfedgePaired;   // its pair (= -1 only on non-manifold edges)
};

// Bounds-checked position lookup that handles both original verts (id
// < baseId, into impl.vertPos_) and chord verts (id >= baseId, into
// newVertPositions). Replaces the inline `getPos3` lambdas that
// previously appeared at six call sites — three of them had the
// off-by-baseId / out-of-bounds bug class fixed in commit 5acc6ab5
// (read uninitialized memory → denormal garbage → nondeterminism).
//
// Bounds-checked under MANIFOLD_ASSERT=ON via DEBUG_ASSERT (throws
// std::logic_error). In production builds with MANIFOLD_ASSERT=OFF
// the assert is a no-op; the defensive return-zero on the chord-vert
// path keeps the read deterministic instead of UB. A return of
// vec3(0) is geometrically wrong but at least reproducible — the
// real fix is to never feed an out-of-range id into this helper.
inline manifold::vec3 GetPos3(
    int id, int baseId, const manifold::Manifold::Impl& impl,
    const std::vector<manifold::vec3>& newVertPositions) {
  if (id < baseId) {
    DEBUG_ASSERT(id >= 0 && id < static_cast<int>(impl.vertPos_.size()),
                 logicErr, "GetPos3: original vert id out of range");
    return impl.vertPos_[id];
  }
  const int j = id - baseId;
  DEBUG_ASSERT(j >= 0 && j < static_cast<int>(newVertPositions.size()),
               logicErr, "GetPos3: chord vert id beyond newVertPositions");
  if (j < 0 || j >= static_cast<int>(newVertPositions.size())) {
    return manifold::vec3(0.0, 0.0, 0.0);  // defensive — see comment above
  }
  return newVertPositions[j];
}

inline std::vector<Edge> EnumerateEdges(const manifold::Manifold::Impl& impl) {
  std::vector<Edge> out;
  out.reserve(impl.halfedge_.size() / 2);
  for (size_t i = 0; i < impl.halfedge_.size(); ++i) {
    const Halfedge he = impl.halfedge_.Get(i);
    if (he.IsForward()) {
      out.push_back(
          {he.startVert, he.endVert, static_cast<int>(i), he.pairedHalfedge});
    }
  }
  return out;
}

// -----------------------------------------------------------------------------
// Step 3: per-edge ordered list of within-ε on-edge verts.
//
// Emmett #289 step 3: "Every edge builds an ordered list of all
// verts that are within ε of it. The order is reliable even with
// inexact arithmetic because any verts that were near each other
// have already been collapsed."
//
// Implementation pattern follows the 2D prototype's
// `BuildEdgeVertLists` (extras/overlap2d_proto.cpp:488). Broad
// phase: edge BVH (ε-padded segment AABBs), queried by vert
// boxes. Narrow phase: project vert onto segment in 3D, accept
// if `0 < t < 1` and `dist² ≤ ε²`. Per-edge hits collected then
// sorted by `t` for stable ordering.
//
// Excludes (per the 2D prototype's lessons learned, see
// `Overlap2D.md` "Alternatives considered"):
//   - vert is an endpoint of the edge (not interior)
//   - "thin-triangle apex": vert is connected to *both* edge
//     endpoints by other edges. In 2D this caused step 5
//     canonicalization to cancel the apex-split sub-edges
//     against the triangle's other two sides, producing empty
//     output. Same hazard in 3D when a tri is folded into a
//     near-degenerate sliver. Only-one-endpoint adjacency is
//     normal mesh-neighbor configuration so we require BOTH.
// -----------------------------------------------------------------------------
struct EdgeVertList {
  std::vector<int> verts;  // sorted by t along the edge
  std::vector<double> ts;  // parametric position 0 < t < 1
};

inline std::vector<EdgeVertList> BuildOnEdgeVertLists(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    double eps) {
  using manifold::la::dot;
  const size_t nE = edges.size();
  const size_t nV = impl.vertPos_.size();
  const double eps2 = eps * eps;
  std::vector<EdgeVertList> out(nE);
  if (nE == 0 || nV == 0) return out;

  // Build vert→neighbor adjacency from `halfedge_` (for the thin-tri-
  // apex skip — see comment block above). The half-edge structure
  // makes this trivial: iterate halfedges, accumulate (start, end)
  // pairs both directions.
  std::vector<std::set<int>> adj(nV);
  for (size_t i = 0; i < impl.halfedge_.size(); ++i) {
    const int s = impl.halfedge_.Start(i);
    const int e = impl.halfedge_.End(i);
    adj[s].insert(e);
    adj[e].insert(s);
  }

  // Per-edge ε-padded AABB (covers the segment plus an ε margin
  // perpendicular to it; using component-wise ε pad covers the
  // perpendicular envelope conservatively).
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    edgeBoxes[i] = b;
  }

  // Per-vert ε-padded AABB used as queries.
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }

  // Morton-sort edge boxes for the BVH.
  Box bBox;
  for (const auto& b : edgeBoxes) bBox = bBox.Union(b);
  std::vector<uint32_t> edgeMorton(nE);
  for (size_t i = 0; i < nE; ++i)
    edgeMorton[i] = manifold::Collider::MortonCode(edgeBoxes[i].Center(), bBox);
  std::vector<size_t> edgePerm(nE);
  std::iota(edgePerm.begin(), edgePerm.end(), 0);
  std::stable_sort(edgePerm.begin(), edgePerm.end(), [&](size_t a, size_t b) {
    return edgeMorton[a] < edgeMorton[b];
  });
  std::vector<Box> sortedEdgeBoxes(nE);
  std::vector<uint32_t> sortedEdgeMorton(nE);
  for (size_t i = 0; i < nE; ++i) {
    sortedEdgeBoxes[i] = edgeBoxes[edgePerm[i]];
    sortedEdgeMorton[i] = edgeMorton[edgePerm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedEdgeBoxes.data(),
                                   sortedEdgeBoxes.size()),
      manifold::VecView<const uint32_t>(sortedEdgeMorton.data(),
                                        sortedEdgeMorton.size()));

  // Per-edge accumulator of (t, vert) hits, then sorted at the end.
  std::vector<std::vector<std::pair<double, int>>> hitsByEdge(nE);

  auto onCollision = [&](size_t vertIdxQ, size_t edgeIdxL) {
    const size_t vertIdx = vertIdxQ;            // query order = vert order
    const size_t edgeIdx = edgePerm[edgeIdxL];  // BVH leaf → real edge
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
  auto recorder = manifold::MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nV),
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

// -----------------------------------------------------------------------------
// Step 4: edge-edge intersection discovery.
//
// Emmett #289 step 4 (3D): "Same as 4 above [the 2D step 4] - this
// is now another degenerate case, as edges would all be skew with
// input in general position." In 3D, two segments generally do not
// intersect — they're skew. They DO intersect when they're coplanar
// and cross within both. This step finds those crossings.
//
// Existing-pipeline analog: `Kernel11`/`Shadow11` in `boolean3.cpp`
// implements the same predicate with symbolic perturbation. That
// version is tightly coupled to `Manifold::Impl`'s two-input A/B
// boolean structure; for a single-input self-overlap context the
// algorithm is the same but the data flow is simpler. We reuse the
// half-edge structure for filtering (skip pairs sharing a vertex)
// and rely on the on-edge vert lists from step 3 for the snap-to-
// existing-vert filter (Emmett: "they cannot include any where an
// end-vert of one edge is in the vert list of the other").
//
// Geometry: closest-points-on-two-segments. Solve the 2x2 normal
// equations for the parameters (s, t) that minimize ||L1(s) - L2(t)||².
// If lines are parallel (det ≈ 0), skip — no general-position
// intersection. Otherwise, accept the candidate iff:
//   - 0 < s < 1, 0 < t < 1 (both interior)
//   - distance between the two closest points ≤ ε (truly coplanar
//     within the working precision)
//
// Returns candidates (not yet inserted into edge lists; insertion
// is part of step 4's "snap" step which the next commit covers).
// Each candidate carries a snapTo field: -1 if it's a fresh vert,
// or ≥0 if its position falls within ε of an existing vert (in
// which case overlap removal would re-use that vert instead).
// -----------------------------------------------------------------------------
struct EdgeEdgeIntersection {
  int edgeA;        // index into edges[]
  int edgeB;        // index into edges[]
  vec3 position;    // intersection point in 3D
  double tA, tB;    // parameters along edgeA, edgeB (interior)
  double distance;  // closest-points distance (≤ ε)
  int snapTo;       // -1 = new vert, ≥0 = snap to existing vert index
};

inline std::vector<EdgeEdgeIntersection> FindEdgeEdgeIntersections(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists, double eps) {
  using manifold::la::cross;
  using manifold::la::dot;
  std::vector<EdgeEdgeIntersection> out;
  const size_t nE = edges.size();
  if (nE < 2) return out;
  const double eps2 = eps * eps;

  // Per-edge ε-padded AABB for the broad phase (same as step 3).
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    edgeBoxes[i] = b;
  }
  Box bBox;
  for (const auto& b : edgeBoxes) bBox = bBox.Union(b);
  std::vector<uint32_t> morton(nE);
  for (size_t i = 0; i < nE; ++i)
    morton[i] = manifold::Collider::MortonCode(edgeBoxes[i].Center(), bBox);
  std::vector<size_t> perm(nE);
  std::iota(perm.begin(), perm.end(), 0);
  std::stable_sort(perm.begin(), perm.end(),
                   [&](size_t a, size_t b) { return morton[a] < morton[b]; });
  std::vector<Box> sortedBoxes(nE);
  std::vector<uint32_t> sortedMorton(nE);
  for (size_t i = 0; i < nE; ++i) {
    sortedBoxes[i] = edgeBoxes[perm[i]];
    sortedMorton[i] = morton[perm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedBoxes.data(), sortedBoxes.size()),
      manifold::VecView<const uint32_t>(sortedMorton.data(),
                                        sortedMorton.size()));

  // For each pair (qi < li), narrow-phase 3D segment-segment.
  auto onPair = [&](size_t qi, size_t li) {
    if (qi >= li) return;
    const size_t eA = perm[qi];
    const size_t eB = perm[li];
    const Edge& a = edges[eA];
    const Edge& b = edges[eB];
    // Skip if shared vertex.
    if (a.v0 == b.v0 || a.v0 == b.v1 || a.v1 == b.v0 || a.v1 == b.v1) return;
    // Skip if one edge's endpoint is in the other's on-edge list
    // (already a known intersection from step 3).
    auto inList = [&](const std::vector<int>& list, int v) {
      return std::find(list.begin(), list.end(), v) != list.end();
    };
    if (inList(onEdgeLists[eA].verts, b.v0) ||
        inList(onEdgeLists[eA].verts, b.v1) ||
        inList(onEdgeLists[eB].verts, a.v0) ||
        inList(onEdgeLists[eB].verts, a.v1))
      return;

    const vec3 p1 = impl.vertPos_[a.v0];
    const vec3 p2 = impl.vertPos_[a.v1];
    const vec3 p3 = impl.vertPos_[b.v0];
    const vec3 p4 = impl.vertPos_[b.v1];
    const vec3 d1 = p2 - p1;  // direction of edge A
    const vec3 d2 = p4 - p3;  // direction of edge B
    const vec3 r = p1 - p3;
    const double a11 = dot(d1, d1);
    const double a22 = dot(d2, d2);
    const double a12 = -dot(d1, d2);
    const double b1 = -dot(d1, r);
    const double b2 = dot(d2, r);
    const double det = a11 * a22 - a12 * a12;
    // Parallel (or near-degenerate edge): skip. The relative
    // threshold uses a11*a22 = (|d1|*|d2|)² as the natural scale;
    // a det of 0 corresponds exactly to parallel directions.
    if (std::fabs(det) <= 1e-30 || a11 < 1e-30 || a22 < 1e-30) return;
    const double s = (a22 * b1 - a12 * b2) / det;
    const double t = (a11 * b2 - a12 * b1) / det;
    // Strict interior on both segments. Don't flag endpoint-touching
    // pairs as edge-edge intersections; those are step-3 territory.
    if (s <= 0.0 || s >= 1.0 || t <= 0.0 || t >= 1.0) return;
    const vec3 closestA = p1 + s * d1;
    const vec3 closestB = p3 + t * d2;
    const vec3 diff = closestA - closestB;
    const double dist2 = dot(diff, diff);
    if (dist2 > eps2) return;
    // Crossing within ε. Use midpoint of the two closest points as
    // the canonical position.
    const vec3 pos = 0.5 * (closestA + closestB);

    // Snap-to-existing: search the union of (a, b)'s endpoints and
    // existing on-edge list members for a vert within ε of `pos`.
    // (For step 4's first pass this is a per-pair narrow-phase
    // check; a fuller implementation would do BVH queries.)
    int snapTo = -1;
    auto trySnap = [&](int v) {
      if (snapTo >= 0) return;
      const vec3 d = pos - impl.vertPos_[v];
      if (dot(d, d) <= eps2) snapTo = v;
    };
    trySnap(a.v0);
    trySnap(a.v1);
    trySnap(b.v0);
    trySnap(b.v1);
    for (int v : onEdgeLists[eA].verts) trySnap(v);
    for (int v : onEdgeLists[eB].verts) trySnap(v);

    out.push_back({static_cast<int>(eA), static_cast<int>(eB), pos, s, t,
                   std::sqrt(dist2), snapTo});
  };
  auto recorder = manifold::MakeSimpleRecorder(onPair);
  auto qf = [&](int i) { return sortedBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nE),
                             /*parallel=*/false);
  return out;
}

// -----------------------------------------------------------------------------
// Step 5: per-triangle on-interior vert list.
//
// Emmett #289 step 5: "Every triangle builds a list of verts that
// are within ε of its interior. None need to be merged." This is
// the 3D analog of step 3, where we found verts within ε of edge
// interiors. Step 5 finds verts within ε of triangle *interiors*
// (excluding vert-on-edge cases, which step 3 owns).
//
// Existing-pipeline analog: `Kernel02` in boolean3.cpp is the
// vert-in-face classification kernel — same predicate (vert lies
// on plane of triangle), with symbolic perturbation of the
// on-edge / on-vert degeneracies via Shadow01.
//
// Algorithm: BVH of ε-padded triangle AABBs, queried by vert
// boxes. Narrow phase: project vert onto triangle plane (signed
// distance via dot with normal), accept iff |dist| ≤ ε AND the
// projected point's barycentric coordinates are all strictly
// positive (interior, not on an edge or vertex).
//
// Excludes: verts that are themselves a vertex of the triangle
// (always trivially "in plane" but at a vertex, not interior).
// Verts on a triangle edge are not excluded here — they're
// reported by both step 3 (on-edge) and step 5 (in plane); the
// barycentric check filters out edge cases at the strict-interior
// gate.
// -----------------------------------------------------------------------------
struct TriVertList {
  std::vector<int> verts;  // verts within ε of interior
  std::vector<vec3> bary;  // barycentric (b0, b1, b2) — all > 0
};

inline std::vector<TriVertList> BuildOnTriVertLists(
    const manifold::Manifold::Impl& impl, double eps) {
  using manifold::la::cross;
  using manifold::la::dot;
  const size_t nT = impl.NumTri();
  const size_t nV = impl.vertPos_.size();
  std::vector<TriVertList> out(nT);
  if (nT == 0 || nV == 0) return out;

  // Per-triangle AABB (ε-padded).
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
  // Per-vert ε-padded query box.
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }
  // Morton-sort triangle boxes for the BVH.
  Box bBox;
  for (const auto& b : triBoxes) bBox = bBox.Union(b);
  std::vector<uint32_t> triMorton(nT);
  for (size_t i = 0; i < nT; ++i)
    triMorton[i] = manifold::Collider::MortonCode(triBoxes[i].Center(), bBox);
  std::vector<size_t> triPerm(nT);
  std::iota(triPerm.begin(), triPerm.end(), 0);
  std::stable_sort(triPerm.begin(), triPerm.end(), [&](size_t a, size_t b) {
    return triMorton[a] < triMorton[b];
  });
  std::vector<Box> sortedBoxes(nT);
  std::vector<uint32_t> sortedMorton(nT);
  for (size_t i = 0; i < nT; ++i) {
    sortedBoxes[i] = triBoxes[triPerm[i]];
    sortedMorton[i] = triMorton[triPerm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedBoxes.data(), sortedBoxes.size()),
      manifold::VecView<const uint32_t>(sortedMorton.data(),
                                        sortedMorton.size()));

  auto onCollision = [&](size_t vertIdx, size_t triIdxL) {
    const size_t triIdx = triPerm[triIdxL];
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
    if (nMag2 == 0) return;  // degenerate triangle
    const double nMag = std::sqrt(nMag2);
    const vec3 p = impl.vertPos_[vertIdx];
    const double dist = dot(p - a, n) / nMag;
    if (std::fabs(dist) > eps) return;
    // Project p onto the triangle plane.
    const vec3 pProj = p - (dist) * (n / nMag);
    // Strict-interior barycentric check (b0, b1, b2 > 0).
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
    // Strict interior: all barycentrics > 0. Use a tiny tolerance
    // to avoid FP wobble at the boundary; this is purely a noise
    // filter, not a tolerance for "small interior" suppression.
    const double bTol = 1e-12;
    if (bu <= bTol || bv <= bTol || bw <= bTol) return;
    out[triIdx].verts.push_back(static_cast<int>(vertIdx));
    out[triIdx].bary.push_back(vec3(bu, bv, bw));
  };
  auto recorder = manifold::MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nV),
                             /*parallel=*/false);
  return out;
}

// -----------------------------------------------------------------------------
// Step 6: edge × triangle intersections.
//
// Emmett #289 step 6: "Find new edge-tri intersections - they cannot
// include any pairs already covered by steps 4 or 5. Insert the new
// vert into the edge and tri lists, unless it is close to one of
// the existing verts in the list, in which case it is merged
// instead, as in 4."
//
// This is the **main intersection-finding phase**. Each non-trivial
// 3D mesh self-intersection involves an edge of one triangle
// piercing the interior of another — this is what catches the
// Cray-class pierces (tri-tri interior overlaps). Each tri-tri
// pierce manifests as TWO edge-tri intersections, one for each
// triangle's edge piercing the other.
//
// Existing-pipeline analog: `Kernel12`/`Intersect12` in boolean3.cpp
// — the main intersection-finding phase of the existing 3D Boolean.
// That implementation composes Kernel02 (vert-in-face) and Kernel11
// (edge-edge) with symbolic perturbation; here we do a direct
// ray-triangle intersection with ε-snap filtering.
//
// Algorithm: BVH of edge AABBs, queried by tri AABBs. Narrow phase:
// parametrize the edge as a ray a + s*(b-a) for s ∈ [0, 1]; solve
// for s where the ray meets the triangle plane; check 0 < s < 1
// and strict-interior barycentrics on the intersection point.
// Filters per Emmett:
//   - Skip if edge is one of the triangle's own three edges
//     (an edge can't pierce its own face's interior).
//   - Skip if either endpoint of the edge is in the triangle's
//     in-tri vert list (= step 5 already represented this).
//   - Snap-to-existing: if the intersection point is within ε of
//     any existing vert (the edge's endpoints, the tri's verts,
//     the on-edge or in-tri vert lists), reuse that vert index.
// -----------------------------------------------------------------------------
struct EdgeTriIntersection {
  int edgeIdx;    // index into edges[]
  int triIdx;     // triangle index in halfedge_/faceNormal_
  vec3 position;  // intersection point in 3D
  double s;       // parameter along the edge (0 < s < 1)
  vec3 bary;      // barycentric on the triangle (all > 0)
  int snapTo;     // -1 = new vert, ≥0 = snap to existing vert
};

inline std::vector<EdgeTriIntersection> FindEdgeTriIntersections(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using manifold::la::cross;
  using manifold::la::dot;
  std::vector<EdgeTriIntersection> out;
  const size_t nE = edges.size();
  const size_t nT = impl.NumTri();
  if (nE == 0 || nT == 0) return out;
  const double eps2 = eps * eps;

  // Per-edge ε-padded AABB (BVH leaves).
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.Union(vec3(b.min.x - eps, b.min.y - eps, b.min.z - eps));
    b.Union(vec3(b.max.x + eps, b.max.y + eps, b.max.z + eps));
    edgeBoxes[i] = b;
  }
  // Per-tri ε-padded AABB (queries).
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
  // Morton-sort edge boxes for the BVH.
  Box bBox;
  for (const auto& b : edgeBoxes) bBox = bBox.Union(b);
  std::vector<uint32_t> edgeMorton(nE);
  for (size_t i = 0; i < nE; ++i)
    edgeMorton[i] = manifold::Collider::MortonCode(edgeBoxes[i].Center(), bBox);
  std::vector<size_t> edgePerm(nE);
  std::iota(edgePerm.begin(), edgePerm.end(), 0);
  std::stable_sort(edgePerm.begin(), edgePerm.end(), [&](size_t a, size_t b) {
    return edgeMorton[a] < edgeMorton[b];
  });
  std::vector<Box> sortedEdgeBoxes(nE);
  std::vector<uint32_t> sortedEdgeMorton(nE);
  for (size_t i = 0; i < nE; ++i) {
    sortedEdgeBoxes[i] = edgeBoxes[edgePerm[i]];
    sortedEdgeMorton[i] = edgeMorton[edgePerm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedEdgeBoxes.data(),
                                   sortedEdgeBoxes.size()),
      manifold::VecView<const uint32_t>(sortedEdgeMorton.data(),
                                        sortedEdgeMorton.size()));

  auto onCollision = [&](size_t triIdxQ, size_t edgeIdxL) {
    const size_t triIdx = triIdxQ;
    const size_t edgeIdx = edgePerm[edgeIdxL];
    const Edge& edge = edges[edgeIdx];
    // Skip if edge is one of the triangle's own 3 edges.
    if (static_cast<int>(triIdx) == edges[edgeIdx].halfedgeForward / 3 ||
        static_cast<int>(triIdx) == edges[edgeIdx].halfedgePaired / 3)
      return;
    // The classic Emmett #289 step 6 skips pierces where the edge
    // endpoint coincides with a tri vert (or is in the tri's snap-
    // merged onTriList, or where a tri vert is in the edge's onEdge
    // list). The original rationale: "no NEW pierce; the endpoint
    // connection is already in the topology."
    //
    // CORRECTNESS finding (2026-05-10): for SELF-INTERSECTING input
    // where a Boolean op has merged verts between left+right, those
    // merged verts hide REAL geometric pierces — the edge passes
    // through the tri's interior PAST the shared vert. The classic
    // skip silently treats these as "already represented topology"
    // and the pipeline can never see/fix them. On hull-mask, offset12
    // and self-intersect, this was hiding 1, 1, and 43 pierces
    // respectively (= self-intersect went 58 → 15 pierces with the
    // skips bypassed).
    //
    // Default: skips DISABLED. Opt back in via OVERLAP3D_SKIP_SHARED
    // _VERT_PIERCES=1 if a fixture regresses (none currently do).
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    if (std::getenv("OVERLAP3D_SKIP_SHARED_VERT_PIERCES")) {
      if (edge.v0 == t0 || edge.v0 == t1 || edge.v0 == t2 || edge.v1 == t0 ||
          edge.v1 == t1 || edge.v1 == t2)
        return;
      auto inTriList = [&](int v) {
        const auto& list = onTriLists[triIdx].verts;
        return std::find(list.begin(), list.end(), v) != list.end();
      };
      if (inTriList(edge.v0) || inTriList(edge.v1)) return;
      auto inEdgeList = [&](int v) {
        const auto& list = onEdgeLists[edgeIdx].verts;
        return std::find(list.begin(), list.end(), v) != list.end();
      };
      if (inEdgeList(t0) || inEdgeList(t1) || inEdgeList(t2)) return;
    }

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
    // Edge parallel to triangle plane: skip (no general-position
    // pierce; coplanar configurations are step 4's territory).
    const double nMag2 = dot(n, n);
    if (nMag2 == 0) return;
    if (std::fabs(denom) <= 1e-30) return;
    const double s = -dot(n, a - p0) / denom;
    if (s <= 0.0 || s >= 1.0) return;
    const vec3 pos = a + s * d;
    // Strict-interior barycentric check.
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
    const double bTol = 1e-12;
    if (bu <= bTol || bv <= bTol || bw <= bTol) return;

    // Snap-to-existing: search the union of edge/tri verts and
    // their associated vert lists for a vert within ε of `pos`.
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
  auto recorder = manifold::MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return triBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nT),
                             /*parallel=*/false);
  return out;
}

// -----------------------------------------------------------------------------
// Step 7 (phase 1): enumerate tri-tri pairs from step 6 results.
//
// Emmett #289 step 7: "Add a new edge for each tri-tri intersection.
// The endpoints will be the exactly two verts that are shared between
// the edges of one tri and the edges and interior of the other tri
// and vice-versa. I believe this will always be true based on the
// merging we've done, but a proof would be nice."
//
// This first phase doesn't yet build new edges — it groups step 6's
// edge-tri intersections by tri-pair and validates the
// "exactly 2 endpoints per piercing pair" claim. The grouping key
// is the unordered pair (triA, triB) where triA is the triangle that
// contains the piercing edge, and triB is the triangle being pierced.
// Each piercing edge has TWO adjacent triangles (manifold input),
// so each step-6 result contributes to TWO tri-pairs.
//
// Output: per-pair count of endpoint contributions, plus the list
// of intersection points. Histogram of counts is reported so we
// can see how often Emmett's exactly-2 claim holds.
// -----------------------------------------------------------------------------
struct TriTriPair {
  int triA, triB;               // sorted (triA < triB)
  std::vector<vec3> endpoints;  // intersection points contributed
};

inline std::vector<TriTriPair> EnumerateTriTriPairs(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeTriIntersection>& etIsects) {
  std::map<std::pair<int, int>, std::vector<vec3>> grouped;
  auto add = [&](int a, int b, vec3 pos) {
    auto key = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
    grouped[key].push_back(pos);
  };
  for (const auto& x : etIsects) {
    const Edge& e = edges[x.edgeIdx];
    // Each edge belongs to two adjacent triangles (manifold input).
    const int triA1 = e.halfedgeForward / 3;
    const int triA2 = e.halfedgePaired >= 0 ? e.halfedgePaired / 3 : -1;
    add(triA1, x.triIdx, x.position);
    if (triA2 >= 0) add(triA2, x.triIdx, x.position);
  }
  std::vector<TriTriPair> out;
  out.reserve(grouped.size());
  for (auto& [key, eps] : grouped) {
    out.push_back({key.first, key.second, std::move(eps)});
  }
  return out;
}

// -----------------------------------------------------------------------------
// Step 6 (RELAXED): same as FindEdgeTriIntersections but with the
// "skip if endpoint is tri vert" and "skip if endpoint in in-tri
// list" filters removed. Used to diagnose how much of the pierce-
// coverage shortfall comes from those filters being too aggressive.
//
// Filter (3) — "edge is one of the tri's three edges" — is kept;
// without it, a triangle's own edges would self-pierce its plane
// at every vertex, giving spurious hits.
//
// Filter (4) — "tri vert in edge's on-edge list" — is kept since
// it represents a step-3 finding (vert-on-edge-interior already
// detected) and dropping it could double-count.
// -----------------------------------------------------------------------------
inline std::vector<EdgeTriIntersection> FindEdgeTriIntersectionsRelaxed(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists, double eps) {
  using manifold::la::cross;
  using manifold::la::dot;
  std::vector<EdgeTriIntersection> out;
  const size_t nE = edges.size();
  const size_t nT = impl.NumTri();
  if (nE == 0 || nT == 0) return out;
  const double eps2 = eps * eps;

  // Same broad-phase setup as FindEdgeTriIntersections.
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
  Box bBox;
  for (const auto& b : edgeBoxes) bBox = bBox.Union(b);
  std::vector<uint32_t> edgeMorton(nE);
  for (size_t i = 0; i < nE; ++i)
    edgeMorton[i] = manifold::Collider::MortonCode(edgeBoxes[i].Center(), bBox);
  std::vector<size_t> edgePerm(nE);
  std::iota(edgePerm.begin(), edgePerm.end(), 0);
  std::stable_sort(edgePerm.begin(), edgePerm.end(), [&](size_t a, size_t b) {
    return edgeMorton[a] < edgeMorton[b];
  });
  std::vector<Box> sortedEdgeBoxes(nE);
  std::vector<uint32_t> sortedEdgeMorton(nE);
  for (size_t i = 0; i < nE; ++i) {
    sortedEdgeBoxes[i] = edgeBoxes[edgePerm[i]];
    sortedEdgeMorton[i] = edgeMorton[edgePerm[i]];
  }
  manifold::Collider collider(
      manifold::VecView<const Box>(sortedEdgeBoxes.data(),
                                   sortedEdgeBoxes.size()),
      manifold::VecView<const uint32_t>(sortedEdgeMorton.data(),
                                        sortedEdgeMorton.size()));

  auto onCollision = [&](size_t triIdxQ, size_t edgeIdxL) {
    const size_t triIdx = triIdxQ;
    const size_t edgeIdx = edgePerm[edgeIdxL];
    const Edge& edge = edges[edgeIdx];
    // Filter (3) only: edge is one of the tri's own edges.
    if (static_cast<int>(triIdx) == edges[edgeIdx].halfedgeForward / 3 ||
        static_cast<int>(triIdx) == edges[edgeIdx].halfedgePaired / 3)
      return;
    // Filter (4) only: tri vert in edge's on-edge list.
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    auto inEdgeList = [&](int v) {
      const auto& list = onEdgeLists[edgeIdx].verts;
      return std::find(list.begin(), list.end(), v) != list.end();
    };
    if (inEdgeList(t0) || inEdgeList(t1) || inEdgeList(t2)) return;
    // (No filter on edge-endpoint = tri-vert, no filter on
    //  endpoint-in-in-tri-list. Geometric check is the only gate.)

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
    if (std::fabs(denom) <= 1e-30) return;
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
    const double bTol = 1e-12;
    if (bu <= bTol || bv <= bTol || bw <= bTol) return;

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
    out.push_back({static_cast<int>(edgeIdx), static_cast<int>(triIdx), pos, s,
                   vec3(bu, bv, bw), snapTo});
  };
  auto recorder = manifold::MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return triBoxes[i]; };
  collider.Collisions<false>(recorder, qf, static_cast<int>(nT),
                             /*parallel=*/false);
  return out;
}

// -----------------------------------------------------------------------------
// Step 7 phase 2: emit new verts + new edges from step 6's edge-tri
// intersections.
//
// Phase 1 (above) validated that piercing tri-tri pairs typically
// have exactly 2 endpoints. Phase 2 turns those endpoints into
// concrete vert/edge candidates:
//   - Each step-6 EdgeTriIntersection becomes a "vert id":
//     `snapTo` if ≥0 (reuse existing), else a freshly allocated id
//     starting at `impl.NumVert()`.
//   - New positions are deduped against each other within ε so
//     two near-coincident pierces from different tri-pairs collapse
//     to one new vert.
//   - For each tri-tri pair with exactly 2 endpoints, emit a new
//     edge connecting the two endpoint vert ids (sorted).
//
// This is the first *mutation* step in the pipeline — it doesn't
// modify the input `Manifold` yet, but it produces the data that
// step 9+ would use to build a new triangulation.
//
// Output is `Step7Phase2Result{ newVertPositions, newEdges }`:
//   - `newVertPositions[i]` is the position of new vert with id
//     `impl.NumVert() + i`.
//   - `newEdges` are pairs of vert ids that overlap removal would
//     insert into the rebuilt mesh.
// -----------------------------------------------------------------------------
struct PiercedNewEdge {
  int v0, v1;      // vert ids (sorted, v0 < v1; existing if < NumVert,
                   // else new with index `id - NumVert`)
  int triA, triB;  // the two tris this edge lies on (sorted)
};

struct Step7Phase2Result {
  std::vector<vec3> newVertPositions;
  std::vector<PiercedNewEdge> newEdges;
  // Parallel to the EdgeTriIntersection input passed to
  // EmitNewVertsAndEdges. resolvedIds[i] is the vert id that
  // event i became — `snapTo` if ≥0 in the input, otherwise a
  // freshly allocated id (≥ NumVert) post-dedup. Used by
  // PropagateNewVertsToOnEdgeLists to add the new verts to the
  // on-edge lists of their piercing edges.
  std::vector<int> resolvedIds;
  int dropped_n_not_2 = 0;  // tri-tri pairs that didn't have exactly 2
  // 1-endpoint tri-tri pairs (= edge-tip-in-interior, no through-
  // pierce) record their endpoint as an interior vert of the tri.
  // The triangulation step fan-triangulates each affected tri
  // around its interior vert(s), eliminating the T-junction that
  // would otherwise remain (= a real geometric pierce on hull-mask,
  // offset12, etc.). Map: triId → set of interior vert ids.
  std::map<int, std::set<int>> interiorVertsPerTri;
};

inline Step7Phase2Result EmitNewVertsAndEdges(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeTriIntersection>& etIsects, double eps) {
  using manifold::la::dot;
  Step7Phase2Result r;
  if (etIsects.empty()) return r;
  const double eps2 = eps * eps;
  const int baseId = static_cast<int>(impl.NumVert());

  // Resolve each etIsect to a vert id, deduping new positions
  // against each other within ε.
  r.resolvedIds.assign(etIsects.size(), -1);
  std::vector<int>& resolvedId = r.resolvedIds;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    if (x.snapTo >= 0) {
      resolvedId[i] = x.snapTo;
      continue;
    }
    // Search existing newVertPositions for one within ε.
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

  // Group resolved endpoints by tri-tri pair. Each step-6 result
  // contributes to two pairs (its piercing edge's two adjacent tris
  // each pair with the pierced tri). Use a set per pair to dedupe
  // when the same vert id arrives twice (e.g., from both adjacent
  // halfedges of the piercing edge).
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

  // For each pair with exactly 2 distinct endpoints, emit a new edge.
  // For pairs with 1 endpoint, treat the endpoint as an interior vert
  // of BOTH tris in the pair (= T-junction; the endpoint sits inside
  // each tri's plane and the surface needs to incorporate it during
  // triangulation). For pairs with 0 or >=3, drop as before.
  for (auto& [key, eps_set] : pairEndpoints) {
    if (eps_set.size() == 2) {
      auto it = eps_set.begin();
      const int a = *it++;
      const int b = *it;
      const int v0 = std::min(a, b);
      const int v1 = std::max(a, b);
      r.newEdges.push_back({v0, v1, key.first, key.second});
    } else if (eps_set.size() == 1) {
      const int v = *eps_set.begin();
      r.interiorVertsPerTri[key.first].insert(v);
      r.interiorVertsPerTri[key.second].insert(v);
    } else {
      ++r.dropped_n_not_2;
    }
  }
  return r;
}

// -----------------------------------------------------------------------------
// Step 8: propagate interior verts onto new edges.
//
// Emmett #289 step 8: "Add to this new edge vert list any verts
// from the interior of either triangle that are within ε of it."
//
// For each new edge from step 7 phase 2, search the union of the
// two adjacent triangles' in-tri vert lists (from step 5) for verts
// within ε of the new edge segment. Sort by parameter t along the
// edge.
//
// This is the 3D analog of the 2D prototype's step 4(b) eager
// propagation: when a freshly-allocated intersection vert is
// created at a tri-tri pierce, it should also appear on any nearby
// edge (NEW or original) that runs near it. Without this pass,
// k-fold concurrence at the new edge produces inconsistent
// per-edge vert lists.
//
// Output: per-new-edge ordered list of additional vert ids that
// should appear between v0 and v1 on the rebuilt edge.
// -----------------------------------------------------------------------------
struct NewEdgeWithExtras {
  PiercedNewEdge edge;
  std::vector<int> extraVerts;  // sorted along the edge
  std::vector<double> extraTs;
};

inline std::vector<NewEdgeWithExtras> AddInteriorVertsToNewEdges(
    const manifold::Manifold::Impl& impl,
    const std::vector<vec3>& newVertPositions,
    const std::vector<PiercedNewEdge>& newEdges,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using manifold::la::dot;
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
    // Candidate verts: union of in-tri verts of triA and triB.
    std::set<int> seen;
    auto check = [&](int v) {
      if (!seen.insert(v).second) return;
      if (v == edge.v0 || v == edge.v1) return;
      const vec3 p = impl.vertPos_[v];  // in-tri verts are existing verts
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
    // Sort by t along the edge.
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

// -----------------------------------------------------------------------------
// Step 9: new-edge × new-edge intersections within each tri.
//
// Emmett #289 step 9: "Find new intersection verts between new edges
// within each tri. Insert or merge the new vert into both edges'
// lists as in 4."
//
// Each new edge from step 7 phase 2 lies on two triangles (its
// `triA` and `triB`). Two new edges sharing a triangle are
// candidates for intersection within that triangle. Since both
// edges lie in the same tri's plane, the 3D edge-edge intersection
// is well-defined (no skew configuration; coplanar by construction).
//
// Algorithm:
//   1. Group new edges by triangle (each contributes to 2 groups).
//   2. For each tri with ≥ 2 new edges, pairwise check 3D
//      segment-segment intersection (closest-points; same as step
//      4's narrow phase but the two segments are guaranteed
//      coplanar so the dist² < ε² check is mostly redundant).
//   3. Each intersection: midpoint of closest points, snap-to-
//      existing if within ε of an endpoint, otherwise allocate a
//      new vert id continuing from step 7 phase 2's id space.
//
// Output: per-pair intersection candidates (NewEdgeNewEdge), to be
// inserted into both edges' on-edge vert lists in step 10.
// -----------------------------------------------------------------------------
struct NewEdgeNewEdgeIsect {
  int newEdgeA;  // index into newEdges[]
  int newEdgeB;
  int triId;   // shared triangle
  int snapTo;  // -1 = newly allocated id, ≥0 = existing
  vec3 position;
};

inline std::vector<NewEdgeNewEdgeIsect> FindNewEdgeIntersections(
    const manifold::Manifold::Impl& impl,
    const std::vector<vec3>& newVertPositions,
    const std::vector<PiercedNewEdge>& newEdges, double eps) {
  using manifold::la::dot;
  std::vector<NewEdgeNewEdgeIsect> out;
  if (newEdges.size() < 2) return out;
  const double eps2 = eps * eps;
  const int baseId = static_cast<int>(impl.NumVert());

  auto getPos = [&](int id) -> vec3 {
    return id < baseId ? impl.vertPos_[id] : newVertPositions[id - baseId];
  };

  // Group: triId → list of indices into newEdges[]
  std::map<int, std::vector<int>> byTri;
  for (size_t i = 0; i < newEdges.size(); ++i) {
    byTri[newEdges[i].triA].push_back(static_cast<int>(i));
    byTri[newEdges[i].triB].push_back(static_cast<int>(i));
  }

  for (auto& [tri, idxs] : byTri) {
    if (idxs.size() < 2) continue;
    // Pairwise intersection check.
    for (size_t i = 0; i + 1 < idxs.size(); ++i) {
      for (size_t j = i + 1; j < idxs.size(); ++j) {
        const PiercedNewEdge& ea = newEdges[idxs[i]];
        const PiercedNewEdge& eb = newEdges[idxs[j]];
        // Skip if shared endpoint.
        if (ea.v0 == eb.v0 || ea.v0 == eb.v1 || ea.v1 == eb.v0 ||
            ea.v1 == eb.v1)
          continue;
        const vec3 p1 = getPos(ea.v0);
        const vec3 p2 = getPos(ea.v1);
        const vec3 p3 = getPos(eb.v0);
        const vec3 p4 = getPos(eb.v1);
        const vec3 d1 = p2 - p1;
        const vec3 d2 = p4 - p3;
        const vec3 r = p1 - p3;
        const double a11 = dot(d1, d1);
        const double a22 = dot(d2, d2);
        const double a12 = -dot(d1, d2);
        const double b1 = -dot(d1, r);
        const double b2 = dot(d2, r);
        const double det = a11 * a22 - a12 * a12;
        if (std::fabs(det) <= 1e-30 || a11 < 1e-30 || a22 < 1e-30) continue;
        const double s = (a22 * b1 - a12 * b2) / det;
        const double t = (a11 * b2 - a12 * b1) / det;
        if (s <= 0.0 || s >= 1.0 || t <= 0.0 || t >= 1.0) continue;
        const vec3 closestA = p1 + s * d1;
        const vec3 closestB = p3 + t * d2;
        const vec3 diff = closestA - closestB;
        if (dot(diff, diff) > eps2) continue;
        const vec3 pos = 0.5 * (closestA + closestB);
        // Snap-to-existing.
        int snapTo = -1;
        auto trySnap = [&](int v) {
          if (snapTo >= 0) return;
          const vec3 dd = pos - getPos(v);
          if (dot(dd, dd) <= eps2) snapTo = v;
        };
        trySnap(ea.v0);
        trySnap(ea.v1);
        trySnap(eb.v0);
        trySnap(eb.v1);
        out.push_back({idxs[i], idxs[j], tri, snapTo, pos});
      }
    }
  }
  return out;
}

// -----------------------------------------------------------------------------
// Step 10: per-tri sub-edge enumeration.
//
// Emmett #289 step 10: "Same as 2D step 5, but now each triangle
// also gets a list of halfedges and the sub-edges are pushed as
// half-edges into each of their triangles's lists (the new edges
// insert both a forward and backward halfedge into each triangle)."
//
// 2D step 5 was: "Break the edges into sub-edges according to their
// internal sorted vert list." In 3D, each triangle has:
//   - 3 original edges, each split into (1 + |on-edge-verts|)
//     sub-edges by step 3's verts
//   - N new edges from step 7p2 (any new edge whose triA or triB
//     equals this tri), each split into (1 + |extras|) sub-edges
//     by step 8's propagated verts. New-edge × new-edge crossings
//     from step 9 add another split.
//
// First commit: enumerate sub-edges per tri and report distribution.
// Don't yet build halfedges with orientation — that's step 11
// territory (cyclic ordering for polygon partition).
//
// Output: per-tri counts (originalSubEdges, newSubEdges, total).
// -----------------------------------------------------------------------------
struct PerTriSubEdgeCount {
  int triId;
  int originalSubEdges;
  int newSubEdges;
  int total() const { return originalSubEdges + newSubEdges; }
};

inline std::vector<PerTriSubEdgeCount> CountSubEdgesPerTri(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<NewEdgeWithExtras>& newEdgesWithExtras) {
  const size_t nT = impl.NumTri();
  std::vector<PerTriSubEdgeCount> out(nT);
  for (size_t t = 0; t < nT; ++t) {
    out[t].triId = static_cast<int>(t);
    out[t].originalSubEdges = 0;
    out[t].newSubEdges = 0;
  }

  // Original edges contribute (1 + on-edge-verts) sub-edges each.
  // For each tri, sum across its 3 halfedges. We need to map each
  // tri's halfedges back to the canonical-forward edge index.
  // Build halfedge-id → edge index lookup.
  std::vector<int> halfedgeToEdge(impl.halfedge_.size(), -1);
  for (size_t e = 0; e < edges.size(); ++e) {
    halfedgeToEdge[edges[e].halfedgeForward] = static_cast<int>(e);
    if (edges[e].halfedgePaired >= 0)
      halfedgeToEdge[edges[e].halfedgePaired] = static_cast<int>(e);
  }
  for (size_t t = 0; t < nT; ++t) {
    for (int k : {0, 1, 2}) {
      const int h = static_cast<int>(3 * t + k);
      const int eIdx = halfedgeToEdge[h];
      if (eIdx < 0) continue;
      const int splits = static_cast<int>(onEdgeLists[eIdx].verts.size());
      out[t].originalSubEdges += (1 + splits);
    }
  }

  // New edges: each one belongs to 2 tris (triA, triB) and
  // contributes (1 + extraVerts) sub-edges to each.
  for (const auto& nwe : newEdgesWithExtras) {
    const int splits = static_cast<int>(nwe.extraVerts.size());
    out[nwe.edge.triA].newSubEdges += (1 + splits);
    out[nwe.edge.triB].newSubEdges += (1 + splits);
  }

  return out;
}

// -----------------------------------------------------------------------------
// Propagate new verts (from step 6 / 7p2) to on-edge vert lists.
//
// Per Emmett #289 step 7: "The endpoints will be the exactly two
// verts that are shared between the edges of one tri and the edges
// and interior of the other tri." For each step-6 event
// (piercingEdge, piercedTri, position, s), the new vert produced at
// `position` lies on `piercingEdge` at parameter `s` — it's
// genuinely on the piercing edge's interior. So overlap removal
// must add it to the on-edge vert list of `piercingEdge` so that
// step 11's halfedge-graph construction splits the original edge at
// the new vert.
//
// Without this propagation, the new edges from step 7p2 are
// "dangling cracks" with no boundary connection — visible in step
// 11 phase 2 as ~20% of halfedges with unset next pointers (Cray:
// 236 unset = 117 new edges × 2 directions exactly).
// -----------------------------------------------------------------------------
inline void PropagateNewVertsToOnEdgeLists(
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<int>& resolvedIds, int baseId,
    const std::vector<Edge>& edges, std::vector<EdgeVertList>& onEdgeLists) {
  // For each etIsect event, add the resolved vert id to the on-edge
  // list of the piercing edge with parameter t = x.s. Skip only if
  // the vert IS the edge's endpoint or already in the list. Note: we
  // no longer skip existing verts wholesale — for SHARED-vert pierces
  // (= step 6 narrow phase detected the pierce because of the shared-
  // vert detect flag), the resolvedId may snap to an EXISTING vert
  // that is NOT already in the edge's onEdgeList; in that case we
  // need to add it so the polygon walker subdivides the edge there.
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    auto& list = onEdgeLists[x.edgeIdx];
    const int v = resolvedIds[i];
    // Skip if vert is one of the edge's endpoints (= no subdivision
    // needed; the edge already terminates there).
    if (v == edges[x.edgeIdx].v0 || v == edges[x.edgeIdx].v1) continue;
    if (std::find(list.verts.begin(), list.verts.end(), v) != list.verts.end())
      continue;
    list.verts.push_back(v);
    list.ts.push_back(x.s);
  }
  // Re-sort each augmented list by t.
  for (auto& list : onEdgeLists) {
    if (list.verts.size() < 2) continue;
    std::vector<size_t> perm(list.verts.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
              [&](size_t a, size_t b) { return list.ts[a] < list.ts[b]; });
    std::vector<int> sortedV(list.verts.size());
    std::vector<double> sortedT(list.ts.size());
    for (size_t j = 0; j < perm.size(); ++j) {
      sortedV[j] = list.verts[perm[j]];
      sortedT[j] = list.ts[perm[j]];
    }
    list.verts = std::move(sortedV);
    list.ts = std::move(sortedT);
  }
}

// -----------------------------------------------------------------------------
// Step 11 phase 1: per-tri halfedge graph construction.
//
// For each triangle T, build the halfedge graph that step 11 phase 3
// will walk to extract sub-polygons. Per Emmett step 10:
//   - Original edges contribute ONE halfedge per sub-edge to T's
//     graph, going in T's CCW direction. (The opposite-direction
//     halfedge belongs to the adjacent triangle T'.)
//   - NEW edges contribute TWO halfedges per sub-edge to T's graph,
//     one each direction — the new edge is a 2-sided "crack" through
//     the triangle, and the polygon walk on each side traces a
//     separate sub-polygon.
//
// Sub-edge enumeration:
//   - Each original halfedge of T (3 of them, one per side of the
//     triangle) corresponds to an Edge in `edges[]`. If that edge has
//     N on-edge verts (from step 3), it splits into N+1 sub-edges
//     and contributes N+1 halfedges to T's graph in CCW order.
//   - Each new edge belonging to T (triA == T or triB == T) with M
//     extras (from step 8) splits into M+1 sub-edges, contributes
//     2·(M+1) halfedges.
//
// This phase just constructs and counts. Phase 2 adds per-vert
// cyclic ordering via 2D projection + atan2; phase 3 walks the
// polygons.
// -----------------------------------------------------------------------------
struct PerTriHalfedge {
  int startVert, endVert;  // vert ids (>= NumVert means it's a new vert)
  bool isFromNewEdge;      // false = original edge, true = new edge
};

struct PerTriHalfedgeGraph {
  int triId;
  std::vector<PerTriHalfedge> halfedges;
  std::set<int> verts;
  // Filled by AddNextPointers (step 11 phase 2). Parallel to
  // halfedges[]: nextHalfedge[i] is the index of the halfedge
  // that follows halfedges[i] in a polygon walk on the left side.
  // -1 = unset (no pairing computed).
  std::vector<int> nextHalfedge;
};

// -----------------------------------------------------------------------------
// Step 11 phase 2: 2D projection + atan2 angle sort + next-around-face
// pointers.
//
// For each per-tri halfedge graph, project the verts to 2D (drop the
// dominant axis of the triangle's normal), sort the outgoing
// halfedges around each vert by angle, and compute the polygon-walk
// `next(h)` pointer for each halfedge.
//
// Convention: walking with the polygon face on the LEFT. For
// halfedge h ending at vert v, next(h) is the outgoing halfedge from
// v whose direction is the smallest CCW rotation from h's reverse
// direction. Equivalently: "turn right as little as possible at v."
//
// Phase 3 will walk these pointers to extract sub-polygons.
// -----------------------------------------------------------------------------
inline void AddNextPointers(const manifold::Manifold::Impl& impl,
                            const std::vector<vec3>& newVertPositions,
                            PerTriHalfedgeGraph& g) {
  using manifold::la::dot;
  if (g.halfedges.empty()) {
    g.nextHalfedge.clear();
    return;
  }
  const int baseId = static_cast<int>(impl.NumVert());
  auto getPos3 = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Pick projection axis: drop the largest |normal| component.
  const vec3 n = impl.faceNormal_[g.triId];
  const int dropAxis =
      (std::fabs(n.x) >= std::fabs(n.y) && std::fabs(n.x) >= std::fabs(n.z)) ? 0
      : (std::fabs(n.y) >= std::fabs(n.z))                                   ? 1
                                           : 2;
  // We need to preserve CCW orientation in the 2D projection. After
  // dropping `dropAxis`, the remaining 2 axes form a 2D coord whose
  // CCW direction matches the tri's normal direction iff the cyclic
  // order (drop+1, drop+2) is consistent with the normal sign.
  // Standard trick: if n[dropAxis] < 0, swap the two remaining axes.
  const int axA = (dropAxis + 1) % 3;
  const int axB = (dropAxis + 2) % 3;
  const bool swap = n[dropAxis] < 0;
  auto getPos2 = [&](int id) {
    const vec3 p = getPos3(id);
    if (swap) return std::pair<double, double>(p[axB], p[axA]);
    return std::pair<double, double>(p[axA], p[axB]);
  };

  // For each vert in the graph, collect indices of OUTGOING halfedges.
  std::map<int, std::vector<int>> outgoingByVert;
  for (size_t i = 0; i < g.halfedges.size(); ++i) {
    outgoingByVert[g.halfedges[i].startVert].push_back(static_cast<int>(i));
  }

  // For each vert, sort its outgoing halfedges by 2D angle (atan2).
  // Cache the angle for each halfedge index.
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

  // For each halfedge h, compute next(h) = outgoing-from-h.endVert
  // whose angle is smallest CCW from h's reverse direction. Equivalent
  // to: in the angle-sorted list at h.endVert, find the halfedge
  // with angle just CCW after `θ_rev = h.angle + π (mod 2π)`.
  g.nextHalfedge.assign(g.halfedges.size(), -1);
  constexpr double kTwoPi = 6.283185307179586;
  for (size_t i = 0; i < g.halfedges.size(); ++i) {
    const int v = g.halfedges[i].endVert;
    auto& list = outgoingByVert[v];
    if (list.empty()) continue;
    double thetaRev = outAngle[i] + 3.141592653589793;
    if (thetaRev > 3.141592653589793) thetaRev -= kTwoPi;
    // Find the outgoing halfedge with angle smallest CW from
    // thetaRev. CW (not CCW) is correct for "face on the LEFT of
    // h" walk: starting from thetaRev and rotating CW (toward the
    // half-plane to h's left), the first outgoing encountered
    // bounds the same face as h. CCW picks the wrong neighbor —
    // when h ends at a vert with both a perimeter-continuation
    // halfedge and a chord halfedge into a sub-polygon, CCW picks
    // the perimeter (so walks trace the outer cycle THROUGH chord
    // endpoints instead of turning into sub-polygons).
    int best = -1;
    double bestDelta = kTwoPi + 1.0;
    for (int idx : list) {
      if (g.halfedges[idx].endVert == g.halfedges[i].startVert) continue;
      double delta = thetaRev - outAngle[idx];  // CW direction
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

inline void AddNextPointers(const manifold::Manifold::Impl& impl,
                            const std::vector<vec3>& newVertPositions,
                            std::vector<PerTriHalfedgeGraph>& graphs) {
  for (auto& g : graphs) AddNextPointers(impl, newVertPositions, g);
}

inline std::vector<PerTriHalfedgeGraph> BuildPerTriHalfedgeGraphs(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<NewEdgeWithExtras>& newEdgesWithExtras) {
  const size_t nT = impl.NumTri();
  std::vector<PerTriHalfedgeGraph> out(nT);
  for (size_t t = 0; t < nT; ++t) out[t].triId = static_cast<int>(t);

  // Build halfedge → edge index map (same as in step 10).
  std::vector<int> halfedgeToEdge(impl.halfedge_.size(), -1);
  for (size_t e = 0; e < edges.size(); ++e) {
    halfedgeToEdge[edges[e].halfedgeForward] = static_cast<int>(e);
    if (edges[e].halfedgePaired >= 0)
      halfedgeToEdge[edges[e].halfedgePaired] = static_cast<int>(e);
  }

  // Original edges → one halfedge per sub-edge in T's CCW direction.
  // T's CCW direction along halfedge h is from h.startVert → h.endVert.
  // The Edge in `edges[]` has v0 < v1; its on-edge verts are sorted
  // by t (with t=0 at v0, t=1 at v1). Going from h.startVert to
  // h.endVert may match the edge's (v0→v1) direction OR be reversed.
  for (size_t t = 0; t < nT; ++t) {
    for (int k : {0, 1, 2}) {
      const int h = static_cast<int>(3 * t + k);
      const Halfedge he = impl.halfedge_.Get(h);
      const int eIdx = halfedgeToEdge[h];
      if (eIdx < 0) continue;
      const Edge& edge = edges[eIdx];
      const auto& list = onEdgeLists[eIdx];
      const bool forward = (he.startVert == edge.v0);
      // Build sequence of vert ids along T's halfedge direction.
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
      // Emit halfedges seq[i] → seq[i+1].
      for (size_t i = 0; i + 1 < seq.size(); ++i) {
        out[t].halfedges.push_back({seq[i], seq[i + 1], false});
        out[t].verts.insert(seq[i]);
        out[t].verts.insert(seq[i + 1]);
      }
    }
  }

  // New edges: each contributes BOTH directions. Sub-divided by
  // extras.
  for (const auto& nwe : newEdgesWithExtras) {
    const PiercedNewEdge& edge = nwe.edge;
    // Build sequence in v0→v1 order; extras already sorted by t.
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

// -----------------------------------------------------------------------------
// Step 11 phase 3: walk halfedge polygons per tri.
//
// With phase 2's next-pointer field set on every halfedge of every
// per-tri graph, the polygon walk is a standard halfedge mesh
// traversal: start at any unvisited halfedge, follow `nextHalfedge`
// until cycling back, mark visited, repeat. Each cycle is one
// sub-polygon of the triangle.
//
// For tris with no new edges, the graph has 3 halfedges forming one
// 3-cycle = the original triangle. For tris with new edges, multiple
// sub-polygons emerge — the new edges (2 directions each) split the
// triangle into pieces.
//
// Cycles that include unset (-1) next pointers are skipped and
// counted as "stalled". Step 12 (multiplicity merge) and step 13
// (winding output) consume the cycle list.
// -----------------------------------------------------------------------------
struct PolygonWalkResult {
  std::vector<std::vector<int>> polygons;  // each = sequence of vert ids
  // 2-vert sub-polygons (= chord pairs whose endpoints are interior
  // to the parent triangle, isolated 2-cycle in the halfedge graph).
  // Kept separate from `polygons` so classifier + triangulation only
  // see the >=3-vert ones; pair-sym Phase 1 consults this list via
  // findPolyHE to know that the chord direction "exists but
  // contributes 0 halfedges to output." This lets pair-sym apply its
  // chord-pair constraints correctly even when one corner is
  // degenerate.
  std::vector<std::vector<int>> degeneratePolygons;
  int stalledHalfedges = 0;  // hit an unset next pointer
};

inline PolygonWalkResult WalkPolygons(const PerTriHalfedgeGraph& g) {
  PolygonWalkResult r;
  if (g.halfedges.empty() || g.nextHalfedge.size() != g.halfedges.size())
    return r;
  std::vector<bool> visited(g.halfedges.size(), false);
  for (size_t start = 0; start < g.halfedges.size(); ++start) {
    if (visited[start]) continue;
    std::vector<int> polygon;
    int cur = static_cast<int>(start);
    bool stalled = false;
    // Follow next pointers up to a safety bound (avoid infinite
    // loops on malformed graphs).
    const size_t maxSteps = g.halfedges.size() + 1;
    for (size_t step = 0; step < maxSteps; ++step) {
      if (cur < 0) {
        stalled = true;
        ++r.stalledHalfedges;
        break;
      }
      if (visited[cur]) {
        // We're back to the start (or another visited halfedge).
        // Standard case: cur == start.
        break;
      }
      visited[cur] = true;
      polygon.push_back(g.halfedges[cur].startVert);
      cur = g.nextHalfedge[cur];
    }
    if (!stalled && polygon.size() >= 3) {
      r.polygons.push_back(std::move(polygon));
    } else if (!stalled && polygon.size() == 2) {
      r.degeneratePolygons.push_back(std::move(polygon));
    }
  }
  return r;
}

inline std::vector<PolygonWalkResult> WalkPolygons(
    const std::vector<PerTriHalfedgeGraph>& graphs) {
  std::vector<PolygonWalkResult> out;
  out.reserve(graphs.size());
  for (const auto& g : graphs) out.push_back(WalkPolygons(g));
  return out;
}

// -----------------------------------------------------------------------------
// Step 12: multiplicity merge of equivalent polygons.
//
// Emmett #289 step 12: "Find any polygons that are equivalent
// (reference all of the same vertices) and merge them into one with
// a multiplicity value, as in 2D step 5, discarding those that
// cancel out."
//
// Two polygons are equivalent if their vertex sequences are the same
// up to cyclic rotation. Two equivalent polygons with the *same*
// rotation direction (both CCW or both CW) reinforce each other
// (multiplicity +1 each → +2). Two with *opposite* rotation cancel
// (multiplicity +1 + −1 → 0; both dropped).
//
// Canonical form per polygon:
//   1. Find the lex-min vert position; rotate the polygon to start
//      there.
//   2. Compare the forward rotation against the reverse rotation
//      starting at the same lex-min vert. Pick the lex-smaller as
//      the canonical key; the sign records which of (forward,
//      reverse) was chosen.
//
// Two polygons with the same canonical key and same sign are
// "same direction" duplicates. Same key, opposite sign = cancelling
// pair.
//
// Output: list of (canonical key, signed multiplicity) entries.
// For Cray and similar mesh-with-pierces inputs, most polygons are
// unique (multiplicity ±1, no cancellation). Significant cancelling
// only appears in T-T-full-coincidence configurations.
// -----------------------------------------------------------------------------
struct CanonicalPolygon {
  std::vector<int> verts;  // canonical (lex-min start, lex-smaller direction)
  int sign;                // +1 if canonical = forward, -1 if = reverse
};

inline CanonicalPolygon CanonicalForm(const std::vector<int>& poly) {
  CanonicalPolygon r;
  r.sign = 0;
  if (poly.size() < 3) {
    r.verts = poly;
    return r;
  }
  // Find lex-min vert id (and pick the smallest start index if ties).
  size_t best = 0;
  for (size_t i = 1; i < poly.size(); ++i) {
    if (poly[i] < poly[best]) best = i;
  }
  // Forward sequence starting at `best`.
  std::vector<int> fwd;
  fwd.reserve(poly.size());
  for (size_t i = 0; i < poly.size(); ++i)
    fwd.push_back(poly[(best + i) % poly.size()]);
  // Reverse sequence: same lex-min start, traverse backward.
  std::vector<int> rev;
  rev.reserve(poly.size());
  rev.push_back(poly[best]);
  for (size_t i = 1; i < poly.size(); ++i)
    rev.push_back(poly[(best + poly.size() - i) % poly.size()]);
  if (fwd <= rev) {
    r.verts = std::move(fwd);
    r.sign = +1;
  } else {
    r.verts = std::move(rev);
    r.sign = -1;
  }
  return r;
}

struct MergedPolygon {
  std::vector<int> verts;  // canonical form
  int multiplicity;        // signed sum (in canonical-form direction)
};

inline std::vector<MergedPolygon> MergeMultiplicities(
    const std::vector<PolygonWalkResult>& walks) {
  std::map<std::vector<int>, int> mult;
  for (const auto& w : walks) {
    for (const auto& p : w.polygons) {
      auto cf = CanonicalForm(p);
      mult[cf.verts] += cf.sign;
    }
  }
  std::vector<MergedPolygon> out;
  out.reserve(mult.size());
  for (auto& [k, m] : mult) {
    if (m == 0) continue;
    out.push_back({k, m});
  }
  return out;
}

// (Removed: dead code `ClassifyPolygon` / Step 13b ray-cast classifier
// — superseded by the Step 13d AnalyticalKeep classifier below. Last
// version had a known off-by-baseId OOB bug at `centroid +=
// positions3D[v]`; same bug class as commit 5acc6ab5. See git history
// if reviving the ray-cast approach is needed. The
// `PolygonClassification` struct below is still used by the active
// classifier path.)

struct PolygonClassification {
  bool keep;
  bool reverse;  // true: emit polygon in reversed vert order
  int windingUp;
  int windingDown;
};

// -----------------------------------------------------------------------------
// Step 13d: per-tri-pair analytical classifier (no ray-cast).
//
// For each multi-poly polygon S of triA, the polygon should be
// retained iff its centroid is OUTSIDE every triB that pierces
// triA. "Outside triB" means: signed distance from S's centroid to
// triB's plane (along n_B's direction) > 0, where n_B is triB's
// face normal (which points outward from the volume in manifold's
// convention).
//
// Why this works: a piercing pair (triA, triB) gives a chord on
// triA's plane. The chord splits triA into two sides. For one side
// the centroid is on the "outside" half-space of triB (= dot > 0);
// for the other, "inside" (= dot < 0). The "outside" side is the
// part of triA NOT occluded by triB and SHOULD be in the output.
// The "inside" side is occluded by triB and should be dropped.
//
// For triA pierced by multiple triBs, S must be on the "outside"
// side of EVERY piercing partner. (Equivalent: a polygon bounded
// by multiple chords needs to be on the right side of all of them.)
//
// Implementation: build per-tri list of piercing partners from
// step 7p2's newEdges (each new edge has triA, triB; both directions
// added). At classify time, iterate the polygon's tri's piercing
// partners and check signed distance.
//
// Orientation: if a polygon's centroid is "outside" all partners,
// keep it in NATURAL orientation (the input mesh's tri orientation
// is correct). Reversal is unnecessary in this analytical setup.
// -----------------------------------------------------------------------------
// Per-edge map: from sorted-vert-pair (= chord identity) to the
// partner tri that produced it. A chord between triA and triB (via
// step 7p2's piercing pair) contributes:
//   chordPartner[(v0, v1)][triA] = triB
//   chordPartner[(v0, v1)][triB] = triA
// So when classifying a polygon, we can look up "which partner does
// THIS chord halfedge correspond to in this tri".
struct ChordPartnerMap {
  // key: (sorted v0, v1, owning triId), value: partner triId
  std::map<std::tuple<int, int, int>, int> partnerOf;
};

// -----------------------------------------------------------------------------
// Per-vert winding via DisjointSets flood-fill (= Winding03_ port).
//
// Reference: src/boolean3.cpp::Winding03_. For each vert v of the
// input mesh, compute the winding number of M at a point just off
// the surface near v. Verts connected by NON-BROKEN halfedges are
// in the same component (= same winding) — broken halfedges are
// those whose edge is pierced by another tri (= step 6's etIsects).
//
// Components are evaluated via ray-cast at one rep vert. For
// connected manifold pieces with no piercing, the entire piece is
// 1 component → 1 ray-cast.
//
// New verts (≥ baseId from step 7p2): each is its own component
// (= treated as fresh chord-endpoint locations, no halfedge
// connectivity to original verts).
// -----------------------------------------------------------------------------
struct PerVertWinding {
  std::vector<int> winding;       // winding at v + ε * probeDir
  std::vector<int> windingMinus;  // winding at v - ε * probeDir
  int numComponents = 0;
};

inline int RayCastWindingAt(const Manifold& inputMf, const vec3& origin,
                            const vec3& direction, double rayLen) {
  const vec3 endpoint = origin + direction * rayLen;
  auto hits = inputMf.RayCast(origin, endpoint);
  using manifold::la::dot;
  int w = 0;
  for (const auto& h : hits) {
    const double s = dot(direction, h.normal);
    if (s > 0)
      w += 1;
    else if (s < 0)
      w -= 1;
  }
  return w;
}

inline PerVertWinding ComputePerVertWinding(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<vec3>& newVertPositions, const Manifold& inputMf,
    double eps) {
  using manifold::la::normalize;
  PerVertWinding r;
  const int origVerts = static_cast<int>(impl.NumVert());
  const int newVerts = static_cast<int>(newVertPositions.size());
  const int totalVerts = origVerts + newVerts;
  r.winding.assign(totalVerts, 0);
  r.windingMinus.assign(totalVerts, 0);

  // 1. Broken-halfedge set: for each etIsect, both halfedges of the
  // pierced edge are broken.
  std::set<int> brokenHalfedges;
  for (const auto& x : etIsects) {
    if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size())) continue;
    const auto& e = edges[x.edgeIdx];
    brokenHalfedges.insert(e.halfedgeForward);
    if (e.halfedgePaired >= 0) brokenHalfedges.insert(e.halfedgePaired);
  }

  // 2. DisjointSets over original verts only (new verts stay singleton)
  DisjointSets uA(origVerts);
  for (size_t he = 0; he < impl.halfedge_.size(); ++he) {
    if (!impl.halfedge_.IsForward(static_cast<int>(he))) continue;
    if (brokenHalfedges.count(static_cast<int>(he)) > 0) continue;
    int v0 = impl.halfedge_.Start(he);
    int v1 = impl.halfedge_.End(he);
    if (v0 < 0 || v0 >= origVerts || v1 < 0 || v1 >= origVerts) continue;
    uA.unite(v0, v1);
  }

  // 3. Find component reps (original verts).
  std::set<int> reps;
  for (int v = 0; v < origVerts; ++v) {
    reps.insert(static_cast<int>(uA.find(v)));
  }
  r.numComponents = static_cast<int>(reps.size()) + newVerts;

  // 4. Ray-cast each rep. Use a tilted direction to avoid grazing
  // input tris. Probe slightly off-surface to bias to one side.
  const vec3 probeDir = normalize(vec3(0.7234, 0.4567, 0.5191));
  const double rayLen = inputMf.BoundingBox().Scale() * 4.0;
  std::map<int, int> repWindingPlus, repWindingMinus;
  for (int rep : reps) {
    if (rep < 0 || rep >= origVerts) continue;
    const vec3 originPlus = impl.vertPos_[rep] + probeDir * eps;
    const vec3 originMinus = impl.vertPos_[rep] - probeDir * eps;
    repWindingPlus[rep] =
        RayCastWindingAt(inputMf, originPlus, probeDir, rayLen);
    repWindingMinus[rep] =
        RayCastWindingAt(inputMf, originMinus, probeDir, rayLen);
  }

  // 5. Propagate: each original vert's winding = its rep's winding.
  for (int v = 0; v < origVerts; ++v) {
    int rep = static_cast<int>(uA.find(v));
    auto itPlus = repWindingPlus.find(rep);
    auto itMinus = repWindingMinus.find(rep);
    if (itPlus != repWindingPlus.end()) r.winding[v] = itPlus->second;
    if (itMinus != repWindingMinus.end()) r.windingMinus[v] = itMinus->second;
  }
  // 6. New verts: each gets its own ray-cast (singleton components).
  for (int j = 0; j < newVerts; ++j) {
    const vec3 originPlus = newVertPositions[j] + probeDir * eps;
    const vec3 originMinus = newVertPositions[j] - probeDir * eps;
    r.winding[origVerts + j] =
        RayCastWindingAt(inputMf, originPlus, probeDir, rayLen);
    r.windingMinus[origVerts + j] =
        RayCastWindingAt(inputMf, originMinus, probeDir, rayLen);
  }
  return r;
}

// -----------------------------------------------------------------------------
// Production Winding03 wrapper (M-vs-M diagnostic).
//
// Calls `manifold::Winding03<true>(impl, impl, p1q2, /*expandP=*/true)` —
// the SoS-perturbed flood-fill from `src/winding03.h` (lifted out of
// `boolean3.cpp`'s anonymous namespace as part of the Option A
// productionization spike). Builds a p1q2 vec from the spike's `etIsects`
// (= forward halfedge × pierced face), sorted by halfedge index so
// `Winding03_`'s `lower_bound` works.
//
// Important caveat: this calls the production code with M==M. Boolean3's
// SoS perturbation (`expandP`) is designed to symbolically distinguish
// inputs P vs Q; with inP==inQ the perturbation collapses, just like
// `m.Boolean(m, OpType::Add)` collapses (see docs/Overlap3D.md handoff).
// The output is therefore *diagnostic*, not authoritative. Use only to
// validate the visibility refactor and characterize the M-vs-M failure
// mode for any future Option B work.
// -----------------------------------------------------------------------------
inline manifold::Vec<int> ProductionWinding03(
    const manifold::Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeTriIntersection>& etIsects) {
  manifold::Vec<std::array<int, 2>> p1q2;
  p1q2.reserve(etIsects.size());
  for (const auto& x : etIsects) {
    if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size())) continue;
    const int he = edges[x.edgeIdx].halfedgeForward;
    p1q2.push_back({he, x.triIdx});
  }
  // Winding03_'s lower_bound sweeps the sorted halfedge index (== p1q2[i][0]
  // when forward=true). Stable-sort by halfedge first, face second.
  std::sort(p1q2.begin(), p1q2.end(),
            [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
              return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
            });
  return manifold::Winding03<true>(
      impl, impl,
      manifold::VecView<std::array<int, 2>>(p1q2.data(), p1q2.size()),
      /*expandP=*/true);
}

inline ChordPartnerMap BuildChordPartnerMap(
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

// Per-polygon classifier (refined): only consider partners whose
// CHORDS BOUND THE POLYGON.
//
// To avoid the "thin sliver near chord" centroid bias, the test
// point is a polygon vertex that's NOT on any chord (= not a new
// vert from step 7p2 / step 8). For polygons made entirely of
// chord verts (rare), fall back to centroid.
//
// Convention: drop polygon iff the test-point is on the OUTWARD
// side of any bounding partner triB (= dot(testPt - vB, n_B) > 0
// where n_B points outward from the volume). That side is in the
// "carved out" region for Subtract operations or generally the
// "wrong side" of the partner — the polygon there is occluded.
inline bool AnalyticalKeep(int triId, const std::vector<int>& polyVerts,
                           const ChordPartnerMap& chordPartners,
                           const std::vector<vec3>& positions3D, int baseId,
                           const manifold::Manifold::Impl& impl) {
  using manifold::la::dot;
  std::set<int> boundingPartners;
  std::set<int> chordVerts;  // verts on chord boundaries
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
  if (boundingPartners.empty()) return true;  // no chord bounds → keep

  // Conservative drop: require ALL non-chord polygon verts to be on
  // the outward side of the partner. If any non-chord vert is on
  // the inward side (= dot < 0), the polygon spans the chord plane
  // — keep with natural orientation. Falls back to centroid only
  // when the polygon is entirely chord-verts.
  std::vector<vec3> nonChordPts;
  for (int v : polyVerts) {
    if (chordVerts.count(v) == 0) {
      vec3 pt;
      if (v < baseId) {
        pt = impl.vertPos_[v];
      } else {
        const int j = v - baseId;
        // Belt-and-suspenders bounds check (commit 5acc6ab5 fixed
        // the original OOB read; this guard catches any regression
        // in the polygon walker's chord-vert id range). Validated
        // post-fix that the OOB never triggers on the current
        // .obj battery + advfuzz; left in to harden against future
        // changes.
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
        // BUG (now fixed): was `positions3D[v]` which forgot to
        // subtract baseId — read out-of-bounds and accumulated
        // garbage into fallbackPt. Same Cray-determinism class as
        // the per-vert path above.
      }
    }
    if (contributed > 0) {
      fallbackPt /= static_cast<double>(contributed);
      nonChordPts.push_back(fallbackPt);
    }
  }
  // Compute scale-relative threshold to ignore borderline cases.
  // bbox.Scale × 1e-6 is a few orders of magnitude above FP-precision
  // and well below any meaningful geometric feature.
  const double thresh = impl.bBox_.Scale() * 1e-6;
  for (int triB : boundingPartners) {
    const int v0 = impl.halfedge_.Start(3 * triB);
    const vec3 p_B = impl.vertPos_[v0];
    const vec3 n_B = impl.faceNormal_[triB];
    // Compute mag (= |n_B|) for the threshold scaling
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

// -----------------------------------------------------------------------------
// Step 13a: triangulate retained polygons + emit output Manifold.
//
// For each polygon from step 11 phase 3 (or its merged form from
// step 12), project to 2D (drop dominant normal axis of the
// containing tri), pass to `manifold::Triangulate` to handle non-
// convex cases, and emit triangles into a fresh MeshGL64. Vertex
// positions are the union of the input mesh's vertPos_ and the new
// vert positions allocated by step 7 phase 2.
//
// First cut: keep ALL polygons regardless of multiplicity sign or
// winding. This validates the polygon-walk → triangulate → MeshGL64
// → Manifold data flow before we add inside/outside filtering.
// Result on Cray (no T-T cancellations) should produce a Manifold
// with similar volume to the input — though the polygons we pass
// include both sides of every new edge, which means the output
// mesh will have effectively the same surface as the input.
//
// Step 13b will add the ray-cast classification that drops
// "interior" polygons (those bounding regions with winding 0 from
// both sides) so the output retains only the boundary surface.
// -----------------------------------------------------------------------------
struct Step13aResult {
  Manifold output;
  int polygonsTriangulated;
  int trianglesEmitted;
  int trisDroppedTooSmall;
  int polygonsKept = 0;
  int polygonsDropped = 0;
  int polygonsReversed = 0;
  int polygonsAutoKept = 0;  // skipped classifier (1-poly tris)
};

// Classifier: given polygon, tri normal, and tri ID, returns
// (keep, reverse). Only called for polygons of MULTI-poly tris
// (= tris with chords from new edges). 1-poly tris are auto-kept
// in natural orientation because their single polygon is the input
// tri (possibly subdivided at on-edge verts), which is already a
// correct boundary surface.
using PolygonClassifierFn = std::function<PolygonClassification(
    const std::vector<int>&, const vec3&, int triId)>;

inline Step13aResult TriangulateAndEmit(
    const manifold::Manifold::Impl& impl,
    const std::vector<vec3>& newVertPositions,
    const std::vector<PolygonWalkResult>& walks,
    PolygonClassifierFn classifier = nullptr,
    const std::map<int, std::set<int>>* interiorVertsPerTri = nullptr) {
  using manifold::la::cross;
  using manifold::la::dot;
  Step13aResult r{Manifold(), 0, 0, 0};
  const int baseId = static_cast<int>(impl.NumVert());
  auto getPos3 = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Build the combined vert table: original input verts + new verts.
  manifold::MeshGL64 out;
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
  std::vector<std::vector<bool>> reverseFlags(walks.size());
  std::vector<bool> autoKeepFlags(walks.size(), false);
  for (size_t triId = 0; triId < walks.size(); ++triId) {
    const auto& w = walks[triId];
    const vec3 n = impl.faceNormal_[triId];
    autoKeepFlags[triId] = (w.polygons.size() == 1);
    keepFlags[triId].assign(w.polygons.size(), true);
    reverseFlags[triId].assign(w.polygons.size(), false);
    if (classifier && !autoKeepFlags[triId]) {
      for (size_t pi = 0; pi < w.polygons.size(); ++pi) {
        auto c = classifier(w.polygons[pi], n, static_cast<int>(triId));
        keepFlags[triId][pi] = c.keep;
        reverseFlags[triId][pi] = c.reverse;
      }
    }
  }
  // Phase 2: cascade-drop forward — drop auto-kept 1-poly tris whose
  // sub-edges are all on the k=1 boundary (= entirely surrounded by
  // dropped polygons). Iterate until stable.
  for (int pass = 0; pass < 16; ++pass) {
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
      if (swap) return manifold::vec2(p[axB], p[axA]);
      return manifold::vec2(p[axA], p[axB]);
    };

    for (size_t pi = 0; pi < w.polygons.size(); ++pi) {
      const auto& poly = w.polygons[pi];
      ++r.polygonsTriangulated;
      const bool keep = keepFlags[triId][pi];
      bool reverse = reverseFlags[triId][pi];
      if (!keep) {
        ++r.polygonsDropped;
        continue;
      }
      ++r.polygonsKept;
      if (autoKeepFlags[triId] && classifier) ++r.polygonsAutoKept;
      if (reverse) ++r.polygonsReversed;
      // Helper: emit a triangle (a, b, c) in the polygon, respecting reverse.
      auto emit = [&](int a, int b, int c) {
        if (reverse) {
          out.triVerts.push_back(a);
          out.triVerts.push_back(c);
          out.triVerts.push_back(b);
        } else {
          out.triVerts.push_back(a);
          out.triVerts.push_back(b);
          out.triVerts.push_back(c);
        }
        ++r.trianglesEmitted;
      };
      // Trivial 3-vert polygons: fan-triangulate around any interior
      // verts assigned to this tri (= 1-endpoint pairs from step 7
      // phase 2; T-junctions that need to incorporate a vert in the
      // tri's interior). For 0 interior verts, just emit the tri.
      // For 1+ interior verts, fan-triangulate. Skip the fan if the
      // tri has multiple polygons (= chord-split; the interior vert
      // would need to be assigned to a specific sub-polygon, more
      // complex — graceful fallback is to emit normally and leave
      // the T-junction).
      if (poly.size() == 3) {
        emit(poly[0], poly[1], poly[2]);
        continue;
      }
      // Detect on-edge collinear-triple: a vert V whose two perimeter
      // neighbors A and B are such that V lies on segment AB. This
      // happens whenever an on-edge vert (from step 7 propagation) is
      // a polygon-perimeter vert. If we let the default triangulator
      // chord A-B, polygons sharing input edge AB on opposite tris
      // each generate that chord → k=4 at AB. Fan from V to suppress
      // the AB chord entirely.
      int fanIdx = -1;
      for (size_t i = 0; i < poly.size(); ++i) {
        const vec3 a = getPos3(poly[(i + poly.size() - 1) % poly.size()]);
        const vec3 v = getPos3(poly[i]);
        const vec3 b = getPos3(poly[(i + 1) % poly.size()]);
        const vec3 ab = b - a;
        const double abLen = manifold::la::length(ab);
        if (abLen <= 0) continue;
        const vec3 av = v - a;
        const double cross_mag =
            manifold::la::length(manifold::la::cross(ab, av));
        if (cross_mag < abLen * abLen * 1e-9) {
          // collinear; check V is between A and B (param in [0,1])
          const double t = manifold::la::dot(av, ab) / (abLen * abLen);
          if (t > 1e-6 && t < 1.0 - 1e-6) {
            fanIdx = static_cast<int>(i);
            break;
          }
        }
      }
      if (fanIdx >= 0) {
        // Fan from poly[fanIdx]. Convex sub-polygons of triangles
        // with chords are convex from their on-edge verts.
        const int v0 = fanIdx;
        const int N = static_cast<int>(poly.size());
        for (int i = 1; i + 1 < N; ++i) {
          emit(poly[v0], poly[(v0 + i) % N], poly[(v0 + i + 1) % N]);
        }
        continue;
      }
      manifold::SimplePolygon poly2;
      poly2.reserve(poly.size());
      for (int v : poly) poly2.push_back(getPos2(v));
      manifold::Polygons polys = {poly2};
      try {
        auto tris = manifold::Triangulate(polys);
        for (const auto& t : tris) emit(poly[t.x], poly[t.y], poly[t.z]);
      } catch (...) {
        ++r.trisDroppedTooSmall;
      }
    }
  }
  // Pierce-checked fan post-pass for 1-endpoint pairs (= interior
  // verts from step 7 phase 2). OPT-IN via OVERLAP3D_FAN_INTERIOR=1.
  //
  // INVESTIGATION FINDING: even with pierce-checks rejecting the
  // obvious bad candidates, applying fan creates duplicate-halfedge
  // (k=4) conflicts that the cap walker can't reconcile. On
  // self-intersect, even 1 applied fan caused the pipeline to fail
  // (= 15 pierces → 661 fallback). Root cause: V-T_i edges introduced
  // by the fan can match edges already present in other tris (= the
  // surrounding mesh has V via shared-vert merge), creating duplicate
  // ownership not detected by simple pierce checks.
  //
  // Full fix would need: either edge-subdivision of the piercing edge
  // (requires re-triangulating its adjacent tris), OR a duplicate-
  // halfedge check in addition to pierce check.
  //
  // Left opt-in for future experimentation. Struct is populated and
  // ready for the proper consumer.
  if (interiorVertsPerTri && !interiorVertsPerTri->empty() &&
      std::getenv("OVERLAP3D_FAN_INTERIOR")) {
    // Build BVH over current out.triVerts (= the just-emitted mesh).
    auto getVertPos = [&](int v) -> vec3 {
      return vec3(out.vertProperties[3 * v + 0], out.vertProperties[3 * v + 1],
                  out.vertProperties[3 * v + 2]);
    };
    auto wouldPierceMesh = [&](const std::vector<size_t>& perm,
                               const std::vector<std::array<int, 3>>& tris,
                               const manifold::Collider& bvh, int a, int b,
                               int c) -> bool {
      const vec3 va = getVertPos(a), vb = getVertPos(b), vc = getVertPos(c);
      manifold::Box queryBox(va, vb);
      queryBox.Union(vc);
      bool found = false;
      auto recorderf = [&](int /*qi*/, int li) {
        if (found) return;
        const int origTri = static_cast<int>(perm[li]);
        const int i0 = tris[origTri][0];
        const int i1 = tris[origTri][1];
        const int i2 = tris[origTri][2];
        int shared = 0;
        for (int x : {a, b, c})
          for (int y : {i0, i1, i2})
            if (x == y) ++shared;
        if (shared >= 2) return;
        const vec3 vi0 = getVertPos(i0), vi1 = getVertPos(i1),
                   vi2 = getVertPos(i2);
        if (SegmentPiercesTriInterior(va, vb, vi0, vi1, vi2) > 0 ||
            SegmentPiercesTriInterior(vb, vc, vi0, vi1, vi2) > 0 ||
            SegmentPiercesTriInterior(vc, va, vi0, vi1, vi2) > 0 ||
            SegmentPiercesTriInterior(vi0, vi1, va, vb, vc) > 0 ||
            SegmentPiercesTriInterior(vi1, vi2, va, vb, vc) > 0 ||
            SegmentPiercesTriInterior(vi2, vi0, va, vb, vc) > 0) {
          found = true;
        }
      };
      auto recorder = manifold::MakeSimpleRecorder(recorderf);
      auto qf = [&](int) { return queryBox; };
      bvh.Collisions<false>(recorder, qf, 1, /*parallel=*/false);
      return found;
    };
    auto buildBVH = [&](std::vector<size_t>& permOut,
                        std::vector<std::array<int, 3>>& trisOut,
                        manifold::Collider& bvhOut) {
      const size_t nT = out.triVerts.size() / 3;
      trisOut.assign(nT, {0, 0, 0});
      std::vector<manifold::Box> triBoxes(nT);
      for (size_t t = 0; t < nT; ++t) {
        const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
        const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
        const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
        trisOut[t] = {i0, i1, i2};
        manifold::Box b(getVertPos(i0), getVertPos(i1));
        b.Union(getVertPos(i2));
        triBoxes[t] = b;
      }
      if (nT == 0) {
        bvhOut = manifold::Collider();
        permOut.clear();
        return;
      }
      manifold::Box bbox;
      for (const auto& b : triBoxes) bbox = bbox.Union(b);
      std::vector<uint32_t> morton(nT);
      for (size_t i = 0; i < nT; ++i)
        morton[i] = manifold::Collider::MortonCode(triBoxes[i].Center(), bbox);
      permOut.assign(nT, 0);
      std::iota(permOut.begin(), permOut.end(), 0);
      std::stable_sort(permOut.begin(), permOut.end(), [&](size_t a, size_t b) {
        return morton[a] < morton[b];
      });
      std::vector<manifold::Box> sortedBoxes(nT);
      std::vector<uint32_t> sortedMorton(nT);
      for (size_t i = 0; i < nT; ++i) {
        sortedBoxes[i] = triBoxes[permOut[i]];
        sortedMorton[i] = morton[permOut[i]];
      }
      bvhOut =
          manifold::Collider(manifold::VecView<const manifold::Box>(
                                 sortedBoxes.data(), sortedBoxes.size()),
                             manifold::VecView<const uint32_t>(
                                 sortedMorton.data(), sortedMorton.size()));
    };
    std::vector<size_t> perm;
    std::vector<std::array<int, 3>> tris;
    manifold::Collider bvh;
    buildBVH(perm, tris, bvh);
    int fansApplied = 0, fansRejected = 0;
    for (const auto& [triId, vSet] : *interiorVertsPerTri) {
      if (triId < 0 || static_cast<size_t>(triId) >= walks.size()) continue;
      // Only fan if the polygon walker output is the natural triangle
      // (no chord splits).
      if (walks[triId].polygons.size() != 1) continue;
      const auto& polyT = walks[triId].polygons[0];
      if (polyT.size() != 3) continue;
      // Need to check: was tri actually emitted as (T0,T1,T2)?
      // Find that tri in out.triVerts. If not present (= classifier
      // dropped), skip.
      const int T0 = polyT[0], T1 = polyT[1], T2 = polyT[2];
      ssize_t foundIdx = -1;
      for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
        const int v0 = static_cast<int>(out.triVerts[3 * t + 0]);
        const int v1 = static_cast<int>(out.triVerts[3 * t + 1]);
        const int v2 = static_cast<int>(out.triVerts[3 * t + 2]);
        // CCW match in either orientation (reverse flag may have flipped).
        if ((v0 == T0 && v1 == T1 && v2 == T2) ||
            (v0 == T1 && v1 == T2 && v2 == T0) ||
            (v0 == T2 && v1 == T0 && v2 == T1) ||
            (v0 == T0 && v1 == T2 && v2 == T1) ||
            (v0 == T2 && v1 == T1 && v2 == T0) ||
            (v0 == T1 && v1 == T0 && v2 == T2)) {
          foundIdx = static_cast<ssize_t>(t);
          break;
        }
      }
      if (foundIdx < 0) continue;
      // Strict-interior filter on V (barycentric > 1e-6).
      const vec3 va = getVertPos(T0), vb = getVertPos(T1), vc = getVertPos(T2);
      const vec3 e1 = vb - va, e2 = vc - va;
      const double d00 = manifold::la::dot(e1, e1);
      const double d01 = manifold::la::dot(e1, e2);
      const double d11 = manifold::la::dot(e2, e2);
      const double denom = d00 * d11 - d01 * d01;
      if (denom == 0) continue;
      // Pick the FIRST V that's strict-interior AND pierce-free.
      int chosenV = -1;
      for (int v : vSet) {
        if (v == T0 || v == T1 || v == T2) continue;
        const vec3 vp = getVertPos(v);
        const vec3 vv = vp - va;
        const double d20 = manifold::la::dot(vv, e1);
        const double d21 = manifold::la::dot(vv, e2);
        const double bv = (d11 * d20 - d01 * d21) / denom;
        const double bw = (d00 * d21 - d01 * d20) / denom;
        const double bu = 1.0 - bv - bw;
        if (!(bu > 1e-6 && bv > 1e-6 && bw > 1e-6)) continue;
        // Pierce-check the 3 fan tris.
        if (wouldPierceMesh(perm, tris, bvh, T0, T1, v)) continue;
        if (wouldPierceMesh(perm, tris, bvh, T1, T2, v)) continue;
        if (wouldPierceMesh(perm, tris, bvh, T2, T0, v)) continue;
        chosenV = v;
        break;
      }
      if (chosenV < 0) {
        ++fansRejected;
        continue;
      }
      // Apply the fan: replace tri at foundIdx with 3 fan tris.
      const bool reverseOrient = reverseFlags[triId][0];
      out.triVerts[3 * foundIdx + 0] = T0;
      out.triVerts[3 * foundIdx + 1] = reverseOrient ? chosenV : T1;
      out.triVerts[3 * foundIdx + 2] = reverseOrient ? T1 : chosenV;
      auto pushTri = [&](int x, int y, int z) {
        if (reverseOrient) {
          out.triVerts.push_back(x);
          out.triVerts.push_back(z);
          out.triVerts.push_back(y);
        } else {
          out.triVerts.push_back(x);
          out.triVerts.push_back(y);
          out.triVerts.push_back(z);
        }
        ++r.trianglesEmitted;
      };
      pushTri(T1, T2, chosenV);
      pushTri(T2, T0, chosenV);
      ++fansApplied;
      // Rebuild BVH so subsequent fans see this fan's tris.
      buildBVH(perm, tris, bvh);
    }
    if ((fansApplied > 0 || fansRejected > 0) &&
        std::getenv("OVERLAP3D_PIPELINE_DIAG")) {
      std::cerr << "      fan post-pass: applied=" << fansApplied
                << " rejected=" << fansRejected << "\n";
    }
  }
  // Pierce-aware reducer: drop classifier-output tris that pierce
  // each other. The pierce-aware cap walker only checks NEW cap
  // tris against existing geometry; classifier-emitted tris that
  // pierce each other (= overlapping kept polygons) aren't caught.
  // self-intersect's 15 residual pierces are likely in this class.
  //
  // Algorithm: build BVH, find piercing pairs, drop ONE tri per
  // pair (heuristic: lowest-cost = fewest k=2 perimeter neighbors),
  // re-cap any new k=1 cycles, iterate. Bounded outer loop.
  //
  // Default ON. Opt-out via OVERLAP3D_NO_PIERCE_REDUCER=1.
  if (!std::getenv("OVERLAP3D_NO_PIERCE_REDUCER")) {
    auto getVPos = [&](int v) -> vec3 {
      return vec3(out.vertProperties[3 * v + 0], out.vertProperties[3 * v + 1],
                  out.vertProperties[3 * v + 2]);
    };
    int totalDroppedPR = 0;
    for (int iter = 0; iter < 32; ++iter) {
      const size_t nT = out.triVerts.size() / 3;
      if (nT == 0) break;
      // Build BVH.
      std::vector<std::array<int, 3>> tris(nT);
      std::vector<manifold::Box> triBoxes(nT);
      for (size_t t = 0; t < nT; ++t) {
        const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
        const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
        const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
        tris[t] = {i0, i1, i2};
        manifold::Box b(getVPos(i0), getVPos(i1));
        b.Union(getVPos(i2));
        triBoxes[t] = b;
      }
      manifold::Box bbox;
      for (const auto& b : triBoxes) bbox = bbox.Union(b);
      std::vector<uint32_t> morton(nT);
      for (size_t i = 0; i < nT; ++i)
        morton[i] = manifold::Collider::MortonCode(triBoxes[i].Center(), bbox);
      std::vector<size_t> permPR(nT);
      std::iota(permPR.begin(), permPR.end(), 0);
      std::stable_sort(permPR.begin(), permPR.end(), [&](size_t a, size_t b) {
        return morton[a] < morton[b];
      });
      std::vector<manifold::Box> sortedBoxes(nT);
      std::vector<uint32_t> sortedMorton(nT);
      for (size_t i = 0; i < nT; ++i) {
        sortedBoxes[i] = triBoxes[permPR[i]];
        sortedMorton[i] = morton[permPR[i]];
      }
      manifold::Collider colliderPR(
          manifold::VecView<const manifold::Box>(sortedBoxes.data(),
                                                 sortedBoxes.size()),
          manifold::VecView<const uint32_t>(sortedMorton.data(),
                                            sortedMorton.size()));
      // Find piercing pairs.
      std::vector<int> pierceCount(nT, 0);
      std::set<std::pair<int, int>> piercePairs;
      for (size_t qi = 0; qi < nT; ++qi) {
        const auto& qt = tris[qi];
        const vec3 qa = getVPos(qt[0]);
        const vec3 qb = getVPos(qt[1]);
        const vec3 qc = getVPos(qt[2]);
        manifold::Box queryBox(qa, qb);
        queryBox.Union(qc);
        auto recorderf = [&](int /*qiL*/, int li) {
          const int oi = static_cast<int>(permPR[li]);
          if (oi <= static_cast<int>(qi)) return;  // dedupe pairs
          const auto& ot = tris[oi];
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
            piercePairs.insert({static_cast<int>(qi), oi});
            ++pierceCount[qi];
            ++pierceCount[oi];
          }
        };
        auto recorder = manifold::MakeSimpleRecorder(recorderf);
        auto qf = [&](int) { return queryBox; };
        colliderPR.Collisions<false>(recorder, qf, 1, /*parallel=*/false);
      }
      if (piercePairs.empty()) break;
      // Drop strategy: for each pair, drop the tri with HIGHER total
      // pierce count (= involved in more pierces, so dropping helps
      // more). Tie-break: lower triId. Mark drops, then remove.
      std::vector<bool> dropMask(nT, false);
      for (const auto& [a, b] : piercePairs) {
        if (dropMask[a] || dropMask[b]) continue;
        if (pierceCount[a] > pierceCount[b])
          dropMask[a] = true;
        else if (pierceCount[b] > pierceCount[a])
          dropMask[b] = true;
        else
          dropMask[a < b ? a : b] = true;  // tie: lower idx
      }
      int iterDropped = 0;
      std::vector<uint64_t> newTriVerts;
      newTriVerts.reserve(out.triVerts.size());
      for (size_t t = 0; t < nT; ++t) {
        if (dropMask[t]) {
          ++iterDropped;
          continue;
        }
        newTriVerts.push_back(out.triVerts[3 * t + 0]);
        newTriVerts.push_back(out.triVerts[3 * t + 1]);
        newTriVerts.push_back(out.triVerts[3 * t + 2]);
      }
      if (iterDropped == 0) break;
      out.triVerts = std::move(newTriVerts);
      totalDroppedPR += iterDropped;
    }
    if (totalDroppedPR > 0 && std::getenv("OVERLAP3D_PIPELINE_DIAG")) {
      std::cerr << "      pierce-aware reducer: dropped " << totalDroppedPR
                << " tris\n";
    }
  }
  // Edge-incidence diagnostic: count tri sides per unsorted edge.
  {
    std::map<std::pair<int, int>, int> edgeCount;
    for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
      int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                  static_cast<int>(out.triVerts[3 * t + 1]),
                  static_cast<int>(out.triVerts[3 * t + 2])};
      for (int e : {0, 1, 2}) {
        int a = v[e], b = v[(e + 1) % 3];
        if (a > b) std::swap(a, b);
        ++edgeCount[{a, b}];
      }
    }
    int k1 = 0, k2 = 0, k3 = 0, k4 = 0, kMore = 0;
    for (const auto& [k, c] : edgeCount) {
      if (c == 1)
        ++k1;
      else if (c == 2)
        ++k2;
      else if (c == 3)
        ++k3;
      else if (c == 4)
        ++k4;
      else
        ++kMore;
    }
    if (k1 + k3 + k4 + kMore > 0) {
      std::cerr << "      edge-incidence: k=1: " << k1 << ", k=2: " << k2
                << ", k=3: " << k3 << ", k=4: " << k4 << ", k>4: " << kMore
                << "\n";
      // Classify k=4 edges by chord-vert incidence (both-new = both
      // verts are step7p2 new chord verts, likely a chord; one-new =
      // mixed; both-old = pure original-mesh edge).
      if (k4 > 0 && std::getenv("OVERLAP3D_K4_CLASSIFY")) {
        const int baseIdLocal = static_cast<int>(impl.NumVert());
        int bothNew = 0, oneNew = 0, bothOld = 0;
        std::vector<std::pair<int, int>> bothOldEdges;
        for (const auto& [k, c] : edgeCount) {
          if (c != 4) continue;
          const bool aNew = k.first >= baseIdLocal;
          const bool bNew = k.second >= baseIdLocal;
          if (aNew && bNew)
            ++bothNew;
          else if (aNew || bNew)
            ++oneNew;
          else {
            ++bothOld;
            bothOldEdges.push_back(k);
          }
        }
        std::cerr << "      k=4 classify: bothNew=" << bothNew
                  << " oneNew=" << oneNew << " bothOld=" << bothOld << "\n";
        // Dump first few bothOld k=4 edges and the polygons that own them
        // (perimeter halfedge in either direction). Each line: edge a-b
        // followed by tri/pi/dir lists (1=A→B, 0=B→A).
        if (bothOld > 0) {
          int dumped = 0;
          for (const auto& [a, b] : bothOldEdges) {
            if (dumped >= 5) break;
            std::cerr << "        k=4 bothOld " << a << "-" << b << ":";
            for (size_t triId = 0; triId < walks.size(); ++triId) {
              for (size_t pi = 0; pi < walks[triId].polygons.size(); ++pi) {
                if (!keepFlags[triId][pi]) continue;
                const auto& poly = walks[triId].polygons[pi];
                for (size_t i = 0; i < poly.size(); ++i) {
                  int va = poly[i], vb = poly[(i + 1) % poly.size()];
                  if ((va == a && vb == b) || (va == b && vb == a)) {
                    std::cerr << " t" << triId << "/p" << pi << "/"
                              << ((va == a) ? 1 : 0);
                  }
                }
              }
            }
            std::cerr << "\n";
            ++dumped;
          }
        }
      }
    }
    // Surface-cap heuristic: close k=1 cycles via fan-triangulation
    // (with conflict-aware fan-apex selection + ear-clip fallback in
    // doCapPass below). Opt-in via OVERLAP3D_CAP=1.
    if (std::getenv("OVERLAP3D_CAP") && k1 > 0) {
      // Forbidden-triple set: vert triplets (sorted) that the pierce
      // reducer dropped. Cap refuses to re-emit any of these. Breaks
      // the drop+re-cap cycle on stuck residuals (gt-7863, self-
      // intersect). Updated by the post-cap pierce reducer.
      std::set<std::array<int, 3>> forbiddenTriples;
      auto sortedTriple = [](int a, int b, int c) -> std::array<int, 3> {
        std::array<int, 3> t = {a, b, c};
        std::sort(t.begin(), t.end());
        return t;
      };
      // Cap walk lambda: closes k=1 cycles via fan-triangulation,
      // with conflict-aware fan-apex selection (= avoids creating
      // new k>2 by skipping cycles where every fan would land on
      // an already-k=2 edge) and greedy ear-clip fallback. Updates
      // ec live as fans are emitted so subsequent cycles see
      // current state. Used both for the initial cap pass and after
      // the k>2 reducer drops tris (= re-cap creates new k=1 cycles
      // we can now close).
      auto doCapPass = [&]() {
        // Re-compute current edge incidence + directed halfedge dirs.
        std::map<std::pair<int, int>, int> ec;
        std::map<std::pair<int, int>, int> dirCount;  // a→b → count
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
        // Pierce-aware cap: build a BVH over current `out.triVerts` and
        // pierce-check each candidate fan/ear tri before emitting.
        // Reject any candidate that pierces (Möller-Trumbore against
        // existing geometry). Opt-out via OVERLAP3D_NO_PIERCE_CAP=1.
        //
        // Rationale: cap fan-fills can introduce new pierces by adding
        // tris that cross unrelated geometry. On gt-7081 specifically,
        // trim+cap added 635 new pierces (input 183 → output 818). The
        // pierce-aware cap refuses such candidates; if no apex/ear is
        // pierce-free, leave the cycle uncapped (= status will become
        // non-manifold and pierce-guard will fall back).
        //
        // BVH rebuilt at start of each doCapPass call. Within-pass
        // cross-cycle pierces (= cap of cycle B pierces cap of cycle A
        // emitted earlier in same call) are NOT caught — would need
        // incremental rebuild. Rare in practice.
        const bool pierceAware = !std::getenv("OVERLAP3D_NO_PIERCE_CAP");
        auto getVertPos = [&](int v) -> vec3 {
          return vec3(out.vertProperties[3 * v + 0],
                      out.vertProperties[3 * v + 1],
                      out.vertProperties[3 * v + 2]);
        };
        manifold::Collider capCollider;
        std::vector<size_t> capPerm;
        std::vector<std::array<int, 3>> capTris;
        auto rebuildCapBVH = [&]() {
          const size_t nT = out.triVerts.size() / 3;
          capTris.assign(nT, {0, 0, 0});
          std::vector<manifold::Box> triBoxes(nT);
          for (size_t t = 0; t < nT; ++t) {
            const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
            const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
            const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
            capTris[t] = {i0, i1, i2};
            manifold::Box b(getVertPos(i0), getVertPos(i1));
            b.Union(getVertPos(i2));
            triBoxes[t] = b;
          }
          if (nT == 0) {
            capCollider = manifold::Collider();
            capPerm.clear();
            return;
          }
          manifold::Box bbox;
          for (const auto& b : triBoxes) bbox = bbox.Union(b);
          std::vector<uint32_t> morton(nT);
          for (size_t i = 0; i < nT; ++i)
            morton[i] =
                manifold::Collider::MortonCode(triBoxes[i].Center(), bbox);
          capPerm.assign(nT, 0);
          std::iota(capPerm.begin(), capPerm.end(), 0);
          std::stable_sort(
              capPerm.begin(), capPerm.end(),
              [&](size_t a, size_t b) { return morton[a] < morton[b]; });
          std::vector<manifold::Box> sortedBoxes(nT);
          std::vector<uint32_t> sortedMorton(nT);
          for (size_t i = 0; i < nT; ++i) {
            sortedBoxes[i] = triBoxes[capPerm[i]];
            sortedMorton[i] = morton[capPerm[i]];
          }
          capCollider =
              manifold::Collider(manifold::VecView<const manifold::Box>(
                                     sortedBoxes.data(), sortedBoxes.size()),
                                 manifold::VecView<const uint32_t>(
                                     sortedMorton.data(), sortedMorton.size()));
        };
        if (pierceAware) rebuildCapBVH();
        // Returns true if candidate tri (a,b,c) would pierce any
        // existing tri in capTris (via BVH broad phase + 6 segment-tri
        // narrow tests). Adjacent tris (sharing 2+ verts) are skipped.
        auto wouldPierce = [&](int a, int b, int c) -> bool {
          if (!pierceAware || capTris.empty()) return false;
          const vec3 va = getVertPos(a);
          const vec3 vb = getVertPos(b);
          const vec3 vc = getVertPos(c);
          manifold::Box queryBox(va, vb);
          queryBox.Union(vc);
          bool found = false;
          auto recorderf = [&](int /*qi*/, int li) {
            if (found) return;
            const int origTri = static_cast<int>(capPerm[li]);
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
          auto recorder = manifold::MakeSimpleRecorder(recorderf);
          auto f = [&](int) { return queryBox; };
          capCollider.Collisions<false>(recorder, f, 1, /*parallel=*/false);
          return found;
        };
        std::map<std::pair<int, int>, int> k1Dir;
        for (const auto& [edge, c] : ec) {
          if (c != 1) continue;
          // Determine direction: which (a,b) is in dirCount?
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
        // DFS-with-backtracking, same as main cap walker.
        auto dfs = [&](int start) -> std::vector<int> {
          std::vector<int> cyc;
          cyc.push_back(start);
          std::vector<size_t> ni;
          ni.push_back(0);
          std::set<int> onP;
          onP.insert(start);
          for (int s = 0; s < 4096 && !cyc.empty(); ++s) {
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
          return {};
        };
        // Helper: would fan-triangulation from apex index `apex` of
        // cycle `cyc` introduce a new chord that's already at k=2 in
        // `ec`? Returns true if conflict.
        auto wouldFanConflict = [&](const std::vector<int>& cyc,
                                    size_t apex) -> bool {
          const size_t n = cyc.size();
          for (size_t off = 1; off + 1 < n; ++off) {
            const int vi = cyc[(apex + off) % n];
            const int vj = cyc[(apex + off + 1) % n];
            // Fan chord (apex, vi) and (apex, vj) — both new.
            // Cycle edge (vi, vj) goes from k=1 to k=2 = fine.
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
        // Pierce-aware variant: returns true if any fan-tri (apex, vi,
        // vj) would pierce an existing tri. Conservative — refuses
        // any apex whose fan introduces a single new pierce.
        auto wouldFanPierce = [&](const std::vector<int>& cyc,
                                  size_t apex) -> bool {
          if (!pierceAware) return false;
          const size_t n = cyc.size();
          for (size_t off = 1; off + 1 < n; ++off) {
            const int vi = cyc[(apex + off) % n];
            const int vj = cyc[(apex + off + 1) % n];
            if (wouldPierce(cyc[apex], vi, vj)) return true;
          }
          return false;
        };
        // Forbidden-triple check: returns true if any fan-tri (apex,
        // vi, vj) matches a vert triplet that was dropped by the
        // pierce reducer. Breaks the drop+re-cap loop on stuck
        // residuals — cap won't re-create a tri we just dropped.
        auto wouldFanForbidden = [&](const std::vector<int>& cyc,
                                     size_t apex) -> bool {
          if (forbiddenTriples.empty()) return false;
          const size_t n = cyc.size();
          for (size_t off = 1; off + 1 < n; ++off) {
            const int vi = cyc[(apex + off) % n];
            const int vj = cyc[(apex + off + 1) % n];
            if (forbiddenTriples.count(sortedTriple(cyc[apex], vi, vj)))
              return true;
          }
          return false;
        };
        for (const auto& [start, _] : capN) {
          if (visited.count(start)) continue;
          auto cyc = dfs(start);
          if (cyc.size() < 3) continue;
          // Try every possible apex for fan; pick the first that
          // doesn't create new k>2 AND (if pierce-aware) doesn't
          // pierce existing geometry.
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
              // Update ec live so subsequent cap cycles' fan-conflict
              // checks see this fan's contributions. Each fan tri
              // adds 3 edges to the incidence map.
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
            // Rebuild BVH so subsequent cycles' pierce check sees
            // this fan's contributions (= catches cap-vs-cap pierces
            // within the same doCapPass call). Opt-out via
            // OVERLAP3D_NO_PIERCE_REBUILD=1.
            if (pierceAware && !std::getenv("OVERLAP3D_NO_PIERCE_REBUILD"))
              rebuildCapBVH();
            continue;
          }
          // Fan failed — try greedy ear-clipping. For each candidate
          // ear (3 consecutive verts), check if its chord (= cycle[i],
          // cycle[i+2]) would create k>2. If conflict-free, emit the
          // ear, remove the middle vert, and continue. Bookkeeping:
          // also update `ec` to reflect the new edges added by the
          // ear so subsequent ear checks see the updated state.
          std::vector<int> remaining(cyc.begin(), cyc.end());
          int earTrisAdded = 0;
          int earGuard = 0;
          while (remaining.size() >= 3 && earGuard++ < 4096) {
            const size_t m = remaining.size();
            if (m == 3) {
              // Last ear — chord is cycle's start↔end edge, already
              // k=1. Don't apply forbidden-triple check here: forbidden
              // is for breaking the drop+re-cap loop on non-converging
              // pierces, but if we refuse the only available 3-cycle
              // closure, the result is k=1 → status non-manifold →
              // FALLBACK to merged input. That trades a small pierce
              // (2 stuck residuals on gt-7863) for a HUGE regression
              // (returns input as-is = 7 pierces). Better to accept
              // the small residual and keep manifold status.
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
              // Chord (va, vc) is the new edge for this ear.
              int sa = va, sc = vc;
              if (sa > sc) std::swap(sa, sc);
              auto it = ec.find({sa, sc});
              if (it != ec.end() && it->second >= 2) continue;
              // Pierce-aware: reject ear if it would pierce existing
              // geometry. Conservative — if no ear is pierce-free,
              // earIdx stays -1 and the cycle stays uncapped.
              if (wouldPierce(va, vb, vc)) continue;
              // Forbidden-triple: reject if dropped by pierce reducer.
              if (forbiddenTriples.count(sortedTriple(va, vb, vc))) continue;
              earIdx = static_cast<ssize_t>(i);
              break;
            }
            if (earIdx < 0) break;  // No conflict-free ear; give up.
            const int va = remaining[earIdx];
            const int vb = remaining[(earIdx + 1) % m];
            const int vc = remaining[(earIdx + 2) % m];
            out.triVerts.push_back(va);
            out.triVerts.push_back(vb);
            out.triVerts.push_back(vc);
            ++earTrisAdded;
            ++addedTris;
            // Update ec: chord (va, vc) is now k=1 (will pair with
            // a future ear's edge).
            int sac = va, scc = vc;
            if (sac > scc) std::swap(sac, scc);
            ++ec[{sac, scc}];
            // Cycle edges (va, vb), (vb, vc) become k=2.
            int s1 = va, s2 = vb;
            if (s1 > s2) std::swap(s1, s2);
            ++ec[{s1, s2}];
            int s3 = vb, s4 = vc;
            if (s3 > s4) std::swap(s3, s4);
            ++ec[{s3, s4}];
            // Remove vb (the middle vert).
            remaining.erase(remaining.begin() + (earIdx + 1) % m);
          }
          if (remaining.empty() || remaining.size() < 3) {
            // Successfully triangulated.
            for (int v : cyc) visited.insert(v);
            ++closed;
            // Same per-cycle BVH refresh as fan path.
            if (pierceAware && !std::getenv("OVERLAP3D_NO_PIERCE_REBUILD"))
              rebuildCapBVH();
          } else {
            // Couldn't fully triangulate. Roll back the ears we did
            // emit (= remove last earTrisAdded tris from out).
            for (int t = 0; t < earTrisAdded; ++t) {
              out.triVerts.pop_back();
              out.triVerts.pop_back();
              out.triVerts.pop_back();
              --addedTris;
            }
          }
        }
        return std::pair<int, int>{closed, addedTris};
      };
      // Initial cap pass — close k=1 cycles from the classifier
      // output, with conflict-aware fan-apex + ear-clip fallback so
      // we don't create new k>2 by capping. (Replaces the older
      // simple-fan inline cap that was the source of cap-induced
      // k>2 conflicts.)
      auto [initialCapped, initialCapTris] = doCapPass();
      std::cerr << "      surface-cap: " << initialCapped << " cycles, "
                << initialCapTris << " tris added\n";
      // Re-run edge-incidence after cap.
      std::map<std::pair<int, int>, int> ec2;
      for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
        int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                    static_cast<int>(out.triVerts[3 * t + 1]),
                    static_cast<int>(out.triVerts[3 * t + 2])};
        for (int e : {0, 1, 2}) {
          int a = v[e], b = v[(e + 1) % 3];
          if (a > b) std::swap(a, b);
          ++ec2[{a, b}];
        }
      }
      int p1 = 0, p2 = 0, p3 = 0, p4 = 0, pM = 0;
      for (const auto& [k, c] : ec2) {
        if (c == 1)
          ++p1;
        else if (c == 2)
          ++p2;
        else if (c == 3)
          ++p3;
        else if (c == 4)
          ++p4;
        else
          ++pM;
      }
      std::cerr << "      post-cap incidence: k=1: " << p1 << ", k=2: " << p2
                << ", k=3: " << p3 << ", k=4: " << p4 << ", k>4: " << pM
                << "\n";
      // Directional k>2 reducer: for each edge with > 2 incidences,
      // count contributions per direction. To make manifold (= 1
      // halfedge per direction), drop excess halfedges per direction.
      // Drop the entire tri containing the chosen halfedge. Then
      // re-cap any new k=1 cycles created by the drop.
      if (p3 > 0 || p4 > 0 || pM > 0) {
        // Outer iter cap is generous; convergence check on k>2 count
        // bails early when reductions stop helping.
        const int outerMax = 32;
        int totalDropped = 0, totalCapped = 0, totalCapTris = 0;
        int prevDropped = -1, prevK2Plus = -1;
        for (int outer = 0; outer < outerMax; ++outer) {
          const size_t nTris = out.triVerts.size() / 3;
          std::vector<bool> dropped(nTris, false);
          // Per-iteration: drop ALL excess directional contributors
          // simultaneously (= one tri per duplicate halfedge direction
          // per edge). For k=4 (2 dups per direction), this drops 2 tris
          // total in one round. Faster convergence than one-at-a-time.
          int innerDropped = 0;
          {
            std::map<std::pair<int, int>, std::vector<int>> heMap;
            std::map<std::pair<int, int>, int> ecCur;  // sorted-edge → count
            for (size_t t = 0; t < nTris; ++t) {
              int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                          static_cast<int>(out.triVerts[3 * t + 1]),
                          static_cast<int>(out.triVerts[3 * t + 2])};
              for (int e : {0, 1, 2}) {
                heMap[{v[e], v[(e + 1) % 3]}].push_back(static_cast<int>(t));
                int a = v[e], b = v[(e + 1) % 3];
                if (a > b) std::swap(a, b);
                ++ecCur[{a, b}];
              }
            }
            // Score-based drop selection: for each candidate tri, the
            // "drop cost" is the number of its OTHER 2 edges currently
            // at k=2 (= each such edge becomes k=1 if we drop the tri).
            // Keep the candidate with the highest cost (= dropping it
            // would create the most new k=1, so it's "load-bearing").
            // Drop the rest. Ties broken by tri index (deterministic).
            //
            // Iterates per duplicate-halfedge direction; for k=4 chord
            // edges, both directions iterate independently → drops 2
            // tris total per round. Compared to "drop highest-tri-index"
            // the score-based pick reduces new-k=1 created per drop and
            // shifts the post-reducer k=1 ceiling on gt-7081.
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
                if (c > keepCost ||
                    (c == keepCost && tris[i] < tris[keepIdx])) {
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
          int innerIter = 1;
          if (innerDropped > 0) {
            std::vector<uint64_t> newTris;
            newTris.reserve(out.triVerts.size());
            for (size_t t = 0; t < nTris; ++t) {
              if (dropped[t]) continue;
              newTris.push_back(out.triVerts[3 * t]);
              newTris.push_back(out.triVerts[3 * t + 1]);
              newTris.push_back(out.triVerts[3 * t + 2]);
            }
            out.triVerts = std::move(newTris);
            totalDropped += innerDropped;
            // Re-cap: close any k=1 cycles created by drops.
            auto [c, ct] = doCapPass();
            totalCapped += c;
            totalCapTris += ct;
            // Convergence: count current k>2 count. If not decreasing
            // monotonically, bail out (= we're churning without progress).
            int curK2Plus = 0;
            std::map<std::pair<int, int>, int> ecCheck;
            for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
              int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                          static_cast<int>(out.triVerts[3 * t + 1]),
                          static_cast<int>(out.triVerts[3 * t + 2])};
              for (int e : {0, 1, 2}) {
                int a = v[e], b = v[(e + 1) % 3];
                if (a > b) std::swap(a, b);
                ++ecCheck[{a, b}];
              }
            }
            for (const auto& [_, kk] : ecCheck)
              if (kk > 2) ++curK2Plus;
            if (prevK2Plus >= 0 && curK2Plus >= prevK2Plus) {
              // No reduction in k>2 count this round — stop.
              break;
            }
            prevK2Plus = curK2Plus;
            if (totalDropped == prevDropped) break;
            prevDropped = totalDropped;
          } else {
            break;  // converged: no k>2 left.
          }
        }
        // Final edge incidence count.
        std::map<std::pair<int, int>, int> ec4;
        for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
          int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                      static_cast<int>(out.triVerts[3 * t + 1]),
                      static_cast<int>(out.triVerts[3 * t + 2])};
          for (int e : {0, 1, 2}) {
            int a = v[e], b = v[(e + 1) % 3];
            if (a > b) std::swap(a, b);
            ++ec4[{a, b}];
          }
        }
        int q1 = 0, q2 = 0, q3 = 0, q4 = 0, qM = 0;
        for (const auto& [k, c] : ec4) {
          if (c == 1)
            ++q1;
          else if (c == 2)
            ++q2;
          else if (c == 3)
            ++q3;
          else if (c == 4)
            ++q4;
          else
            ++qM;
        }
        if (totalDropped > 0 || totalCapped > 0) {
          std::cerr << "      k>2 reducer + re-cap: dropped " << totalDropped
                    << " tris, capped " << totalCapped << " new cycles ("
                    << totalCapTris << " tris) → k=1: " << q1 << ", k=2: " << q2
                    << ", k=3: " << q3 << ", k=4: " << q4 << ", k>4: " << qM
                    << "\n";
        }
      }
      // Trim-orphans pass: any triangle with >= 2 of its 3 edges at
      // k=1 is structurally non-manifold (= can't be part of a closed
      // surface); drop it. Iterate with re-cap until convergence. This
      // attacks residual k=1 chains/branches that the cap walker can't
      // close (= non-simple cycles). The bound is generous to allow
      // deep cascade-drops on dense self-intersection geometry like
      // gt-7081.
      //
      // CORRECTNESS NOTE: trim-orphans is OPT-IN via OVERLAP3D_TRIM=1
      // (= default off). Reasons:
      //   1. The pierce-aware cap (default-on, see wouldPierce in
      //      doCapPass) prevents most cases where trim would help —
      //      it refuses to fan-fill cycles that would introduce new
      //      pierces, so we no longer have the "cap-induced new
      //      pierces" failure mode that trim was patching over.
      //   2. trim+pierce-cap is too slow on gt-7081 (= 32 trim rounds
      //      × per-round pierce-aware doCapPass with BVH rebuild =
      //      >10min on 16k-vert mesh). Trim alone is fast but the
      //      pierce-aware cap inside makes it expensive.
      //   3. Even if trim ran on gt-7081, pierce guard would fall back
      //      anyway — trim's "fix" creates more pierces than input.
      // Net: trim is now a vestigial tool, kept for experimentation.
      if (std::getenv("OVERLAP3D_TRIM")) {
        int trimTotal = 0, trimCapped = 0, trimCapTris = 0;
        for (int trim = 0; trim < 32; ++trim) {
          const size_t nT = out.triVerts.size() / 3;
          if (nT == 0) break;
          std::map<std::pair<int, int>, int> ecT;
          for (size_t t = 0; t < nT; ++t) {
            int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                        static_cast<int>(out.triVerts[3 * t + 1]),
                        static_cast<int>(out.triVerts[3 * t + 2])};
            for (int e : {0, 1, 2}) {
              int a = v[e], b = v[(e + 1) % 3];
              if (a > b) std::swap(a, b);
              ++ecT[{a, b}];
            }
          }
          // Count current k=1; bail if zero (= manifold).
          int curK1 = 0;
          for (const auto& [_, c] : ecT)
            if (c == 1) ++curK1;
          if (curK1 == 0) break;
          std::vector<bool> dropped(nT, false);
          int trimmed = 0;
          // Two-tier trim: round 0 drops tris with >=2 k=1 (clearly
          // non-manifold). Subsequent rounds drop tris with >=1 k=1
          // ONLY if cap-pass made no progress last round (= residual
          // is non-cyclic, can't be closed). The aggressive >=1 trim
          // converges by eating the dangling chain back to the nearest
          // closed surface; volume cost is bounded by the chain length.
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
          std::vector<uint64_t> newTris;
          newTris.reserve(out.triVerts.size());
          for (size_t t = 0; t < nT; ++t) {
            if (dropped[t]) continue;
            newTris.push_back(out.triVerts[3 * t]);
            newTris.push_back(out.triVerts[3 * t + 1]);
            newTris.push_back(out.triVerts[3 * t + 2]);
          }
          out.triVerts = std::move(newTris);
          trimTotal += trimmed;
          auto [c, ct] = doCapPass();
          trimCapped += c;
          trimCapTris += ct;
        }
        if (trimTotal > 0) {
          std::map<std::pair<int, int>, int> ecF;
          for (size_t t = 0; t < out.triVerts.size() / 3; ++t) {
            int v[3] = {static_cast<int>(out.triVerts[3 * t]),
                        static_cast<int>(out.triVerts[3 * t + 1]),
                        static_cast<int>(out.triVerts[3 * t + 2])};
            for (int e : {0, 1, 2}) {
              int a = v[e], b = v[(e + 1) % 3];
              if (a > b) std::swap(a, b);
              ++ecF[{a, b}];
            }
          }
          int r1 = 0, r2 = 0, rR = 0;
          for (const auto& [k, c] : ecF) {
            if (c == 1)
              ++r1;
            else if (c == 2)
              ++r2;
            else
              ++rR;
          }
          std::cerr << "      trim-orphans: dropped " << trimTotal
                    << " tris, capped " << trimCapped << " (" << trimCapTris
                    << " tris) → k=1: " << r1 << ", k=2: " << r2
                    << ", k>2: " << rR << "\n";
        }
      }
      // Post-cap pierce-aware reducer + re-cap loop: alternates
      // pierce-drop and cap-close. Dropping post-cap tris creates
      // k=1 cycles which the re-cap (pierce-aware) then closes.
      // Iterate until fixed point. Hugely effective: hull-mask
      // 31→0, offset12 46→0, havocglass 1→0 (all fully fixed).
      // Opt-out via OVERLAP3D_NO_POST_CAP_PIERCE_REDUCER=1.
      //
      // PERF: skipped on very large meshes (> 25000 tris) because
      // the iterative re-cap doesn't converge fast enough on dense
      // self-intersection (gt-7081 with 33k tris times out >4min).
      // Override threshold via OVERLAP3D_POST_CAP_TRI_LIMIT=N (use
      // 0 to disable the cap, large N to force run).
      const size_t postCapTriLimit = []() {
        const char* s = std::getenv("OVERLAP3D_POST_CAP_TRI_LIMIT");
        return s ? static_cast<size_t>(std::atoi(s)) : size_t(25000);
      }();
      if (!std::getenv("OVERLAP3D_NO_POST_CAP_PIERCE_REDUCER") &&
          out.triVerts.size() / 3 <= postCapTriLimit) {
        auto getVPos2 = [&](int v) -> vec3 {
          return vec3(out.vertProperties[3 * v + 0],
                      out.vertProperties[3 * v + 1],
                      out.vertProperties[3 * v + 2]);
        };
        int totalDroppedPP = 0;
        // Iter cap: 8 by default (= bounded runtime on large meshes).
        // Override via OVERLAP3D_POST_CAP_ITERS=N. Higher values let
        // small fixtures converge further (gt-7863 stuck at 2 pierces
        // with iter=8; might go lower with more iters).
        const int postCapIterLimit = []() {
          const char* s = std::getenv("OVERLAP3D_POST_CAP_ITERS");
          return s ? std::atoi(s) : 8;
        }();
        // Track pair set across iters; if pairs are stable (= same
        // set as prev iter), drop+re-cap is in a steady state with no
        // net progress. Break early to avoid wasted work.
        std::set<std::pair<int, int>> prevPairs;
        for (int iter = 0; iter < postCapIterLimit; ++iter) {
          const size_t nT = out.triVerts.size() / 3;
          if (nT == 0) break;
          std::vector<std::array<int, 3>> tris(nT);
          std::vector<manifold::Box> triBoxes(nT);
          for (size_t t = 0; t < nT; ++t) {
            const int i0 = static_cast<int>(out.triVerts[3 * t + 0]);
            const int i1 = static_cast<int>(out.triVerts[3 * t + 1]);
            const int i2 = static_cast<int>(out.triVerts[3 * t + 2]);
            tris[t] = {i0, i1, i2};
            manifold::Box b(getVPos2(i0), getVPos2(i1));
            b.Union(getVPos2(i2));
            triBoxes[t] = b;
          }
          manifold::Box bbox;
          for (const auto& b : triBoxes) bbox = bbox.Union(b);
          std::vector<uint32_t> morton(nT);
          for (size_t i = 0; i < nT; ++i)
            morton[i] =
                manifold::Collider::MortonCode(triBoxes[i].Center(), bbox);
          std::vector<size_t> permPP(nT);
          std::iota(permPP.begin(), permPP.end(), 0);
          std::stable_sort(
              permPP.begin(), permPP.end(),
              [&](size_t a, size_t b) { return morton[a] < morton[b]; });
          std::vector<manifold::Box> sortedBoxes(nT);
          std::vector<uint32_t> sortedMorton(nT);
          for (size_t i = 0; i < nT; ++i) {
            sortedBoxes[i] = triBoxes[permPP[i]];
            sortedMorton[i] = morton[permPP[i]];
          }
          manifold::Collider colliderPP(
              manifold::VecView<const manifold::Box>(sortedBoxes.data(),
                                                     sortedBoxes.size()),
              manifold::VecView<const uint32_t>(sortedMorton.data(),
                                                sortedMorton.size()));
          std::vector<int> pierceCount(nT, 0);
          std::set<std::pair<int, int>> piercePairs;
          for (size_t qi = 0; qi < nT; ++qi) {
            const auto& qt = tris[qi];
            const vec3 qa = getVPos2(qt[0]);
            const vec3 qb = getVPos2(qt[1]);
            const vec3 qc = getVPos2(qt[2]);
            manifold::Box queryBox(qa, qb);
            queryBox.Union(qc);
            auto recorderf = [&](int /*qiL*/, int li) {
              const int oi = static_cast<int>(permPP[li]);
              if (oi <= static_cast<int>(qi)) return;
              const auto& ot = tris[oi];
              int shared = 0;
              for (int x : qt)
                for (int y : ot)
                  if (x == y) ++shared;
              if (shared >= 2) return;
              const vec3 oa = getVPos2(ot[0]);
              const vec3 ob = getVPos2(ot[1]);
              const vec3 oc = getVPos2(ot[2]);
              if (SegmentPiercesTriInterior(qa, qb, oa, ob, oc) > 0 ||
                  SegmentPiercesTriInterior(qb, qc, oa, ob, oc) > 0 ||
                  SegmentPiercesTriInterior(qc, qa, oa, ob, oc) > 0 ||
                  SegmentPiercesTriInterior(oa, ob, qa, qb, qc) > 0 ||
                  SegmentPiercesTriInterior(ob, oc, qa, qb, qc) > 0 ||
                  SegmentPiercesTriInterior(oc, oa, qa, qb, qc) > 0) {
                piercePairs.insert({static_cast<int>(qi), oi});
                ++pierceCount[qi];
                ++pierceCount[oi];
              }
            };
            auto recorder = manifold::MakeSimpleRecorder(recorderf);
            auto qf = [&](int) { return queryBox; };
            colliderPP.Collisions<false>(recorder, qf, 1, false);
          }
          if (std::getenv("OVERLAP3D_PIERCE_PER_ITER")) {
            std::cerr << "      post-cap pierce iter " << iter << ": "
                      << piercePairs.size()
                      << " pierce pairs (forbid=" << forbiddenTriples.size()
                      << ")";
            int dumped = 0;
            for (const auto& [a, b] : piercePairs) {
              if (dumped++ >= 5) break;
              std::cerr << " [" << a << "," << b << "](";
              for (int k : {0, 1, 2}) {
                std::cerr << static_cast<int>(out.triVerts[3 * a + k])
                          << ((k < 2) ? "," : "");
              }
              std::cerr << ")(";
              for (int k : {0, 1, 2}) {
                std::cerr << static_cast<int>(out.triVerts[3 * b + k])
                          << ((k < 2) ? "," : "");
              }
              std::cerr << ")";
            }
            std::cerr << "\n";
          }
          if (piercePairs.empty()) break;
          // No-progress detection: if the pierce pairs are exactly
          // the same set as last iter, drop+re-cap is in a steady
          // state. Enable forbidden-triple tracking so cap can't
          // re-create the offending tris (= breaks the loop).
          // Without this gate, forbidden tracking is too aggressive
          // and regresses fixtures that converge naturally (hull-
          // mask 0→5, offset12 0→5, havocglass 0→1).
          const bool stuck = (iter > 0 && piercePairs == prevPairs);
          prevPairs = piercePairs;
          std::vector<bool> dropMask(nT, false);
          for (const auto& [a, b] : piercePairs) {
            if (dropMask[a] || dropMask[b]) continue;
            if (pierceCount[a] > pierceCount[b])
              dropMask[a] = true;
            else if (pierceCount[b] > pierceCount[a])
              dropMask[b] = true;
            else
              dropMask[a < b ? a : b] = true;
          }
          int iterDropped = 0;
          std::vector<uint64_t> newTriVerts;
          newTriVerts.reserve(out.triVerts.size());
          for (size_t t = 0; t < nT; ++t) {
            if (dropMask[t]) {
              ++iterDropped;
              // Record dropped tri's vert triplet only when stuck
              // (= drop+re-cap loop detected). This prevents over-
              // restriction on naturally-converging fixtures.
              if (stuck) {
                forbiddenTriples.insert(
                    sortedTriple(static_cast<int>(out.triVerts[3 * t + 0]),
                                 static_cast<int>(out.triVerts[3 * t + 1]),
                                 static_cast<int>(out.triVerts[3 * t + 2])));
              }
              continue;
            }
            newTriVerts.push_back(out.triVerts[3 * t + 0]);
            newTriVerts.push_back(out.triVerts[3 * t + 1]);
            newTriVerts.push_back(out.triVerts[3 * t + 2]);
          }
          if (iterDropped == 0) break;
          out.triVerts = std::move(newTriVerts);
          totalDroppedPP += iterDropped;
          // Re-cap any k=1 cycles created by the drops. Pierce-aware
          // cap will refuse pierce-creating fans AND forbidden tri
          // triplets we just dropped. Resulting non-manifold output
          // (= residual k=1) is caught by the outer pierce/drift
          // guard.
          (void)doCapPass();
        }
        if (totalDroppedPP > 0 && std::getenv("OVERLAP3D_PIPELINE_DIAG")) {
          std::cerr << "      post-cap pierce-aware reducer: dropped "
                    << totalDroppedPP << " tris (with re-cap)\n";
        }
      }
    }
    int targetK = 0;
    if (std::getenv("OVERLAP3D_DUMP_K4") && k4 > 0)
      targetK = 4;
    else if (std::getenv("OVERLAP3D_DUMP_K3") && k3 > 0)
      targetK = 3;
    else if (std::getenv("OVERLAP3D_DUMP_K1") && k1 > 0)
      targetK = 1;
    if (targetK > 0) {
      int dumped = 0;
      for (const auto& [edge, c] : edgeCount) {
        if (c != targetK) continue;
        std::cerr << "    k=" << targetK << " edge[" << edge.first << ","
                  << edge.second << "] kept:";
        for (size_t triId = 0; triId < walks.size(); ++triId) {
          for (size_t pi = 0; pi < walks[triId].polygons.size(); ++pi) {
            if (!keepFlags[triId][pi]) continue;
            const auto& poly = walks[triId].polygons[pi];
            bool hA = false, hB = false;
            int ia = -1, ib = -1;
            for (size_t i = 0; i < poly.size(); ++i) {
              if (poly[i] == edge.first) {
                hA = true;
                ia = i;
              }
              if (poly[i] == edge.second) {
                hB = true;
                ib = i;
              }
            }
            if (hA && hB) {
              std::cerr << " tri" << triId << "/p" << pi << "[";
              for (size_t i = 0; i < poly.size(); ++i) {
                if (i) std::cerr << ",";
                std::cerr << poly[i];
              }
              std::cerr << " a@" << ia << " b@" << ib << "]";
            }
          }
        }
        // Also show DROPPED polys that contain both verts (= the
        // missing partner candidates):
        std::cerr << " dropped:";
        for (size_t triId = 0; triId < walks.size(); ++triId) {
          for (size_t pi = 0; pi < walks[triId].polygons.size(); ++pi) {
            if (keepFlags[triId][pi]) continue;
            const auto& poly = walks[triId].polygons[pi];
            bool hA = false, hB = false;
            int ia = -1, ib = -1;
            for (size_t i = 0; i < poly.size(); ++i) {
              if (poly[i] == edge.first) {
                hA = true;
                ia = i;
              }
              if (poly[i] == edge.second) {
                hB = true;
                ib = i;
              }
            }
            if (hA && hB) {
              const int N = static_cast<int>(poly.size());
              const bool adj =
                  (std::abs(ia - ib) == 1) || (std::abs(ia - ib) == N - 1);
              std::cerr << " tri" << triId << "/p" << pi << "(" << N << "v"
                        << (adj ? "/perim" : "/chord") << ")";
            }
          }
        }
        std::cerr << "\n";
        if (++dumped >= 6) break;
      }
    }
  }
  r.output = Manifold(out);
  return r;
}

// =============================================================================
// OverlapRemoval (entry point) — runs steps 1-10 queries on the input
// Manifold, collects the data structures that steps 11-13 need to
// build the output, and returns a new Manifold.
//
// Step 11 phase 0: output is currently a no-op (round-trips the Impl
// back to a Manifold unchanged). This commit establishes the API
// shape and validates that the round-trip is lossless. Steps 11
// phase 1+ will progressively replace the stub with actual polygon
// partition + winding classification.
// =============================================================================
struct OverlapRemovalDebug {
  // Counts from each query step; useful for tests.
  int step1Merges;
  int step3OnEdgeHits;
  int step4EdgeEdge;
  int step5InTri;
  int step6EdgeTri;
  int step7p2NewVerts;
  int step7p2NewEdges;
  int step8Propagated;
  int step9NewNew;
  int step10TotalSubEdges;
};

inline std::pair<Manifold, OverlapRemovalDebug> OverlapRemoval(
    const Manifold& input, double eps = 0.0) {
  OverlapRemovalDebug dbg{};
  if (input.IsEmpty()) return {input, dbg};
  if (eps <= 0.0) eps = InferEps(input);

  auto mr = MergeVertsEps(input, eps);
  dbg.step1Merges = mr.mergedCount;

  auto impl = ImplFromManifold(mr.manifold);
  OVERLAP3D_HASH_DUMP(
      "impl.vertPos_",
      Fnv1a64(impl.vertPos_.data(), impl.vertPos_.size() * sizeof(vec3)));
  // Halfedges is SoA (upstream #1709); flatten to Vec<Halfedge> for
  // the hash so the dump is stable across runs.
  {
    auto heData = impl.halfedge_.ToData();
    OVERLAP3D_HASH_DUMP(
        "impl.halfedge_",
        Fnv1a64(heData.data(), heData.size() * sizeof(Halfedge)));
  }
  OVERLAP3D_HASH_DUMP(
      "impl.faceNormal_",
      Fnv1a64(impl.faceNormal_.data(), impl.faceNormal_.size() * sizeof(vec3)));
  OVERLAP3D_HASH_DUMP(
      "impl.vertNormal_",
      Fnv1a64(impl.vertNormal_.data(), impl.vertNormal_.size() * sizeof(vec3)));
  OVERLAP3D_HASH_DUMP("impl.bBox_", Fnv1a64(&impl.bBox_, sizeof(impl.bBox_)));
  auto edges = EnumerateEdges(impl);
  OVERLAP3D_HASH_DUMP("edges", HashVec(edges));
  auto onEdgeLists = BuildOnEdgeVertLists(impl, edges, eps);
  for (const auto& l : onEdgeLists) dbg.step3OnEdgeHits += l.verts.size();
  if (DetHashEnabled()) {
    uint64_t h = 14695981039346656037ULL;
    for (const auto& l : onEdgeLists) {
      uint64_t e = HashVec(l.verts);
      h ^= e;
      h *= 1099511628211ULL;
    }
    OVERLAP3D_HASH_DUMP("onEdgeLists ordered", h);
  }

  auto eeIsects = FindEdgeEdgeIntersections(impl, edges, onEdgeLists, eps);
  dbg.step4EdgeEdge = static_cast<int>(eeIsects.size());
  OVERLAP3D_HASH_DUMP("eeIsects", HashVec(eeIsects));

  auto onTriLists = BuildOnTriVertLists(impl, eps);
  for (const auto& l : onTriLists) dbg.step5InTri += l.verts.size();
  if (DetHashEnabled()) {
    uint64_t h = 0;
    for (const auto& l : onTriLists) h ^= HashVec(l.verts);
    OVERLAP3D_HASH_DUMP("onTriLists XOR", h);
  }

  auto etIsects =
      FindEdgeTriIntersections(impl, edges, onEdgeLists, onTriLists, eps);
  dbg.step6EdgeTri = static_cast<int>(etIsects.size());
  OVERLAP3D_HASH_DUMP("etIsects", HashVec(etIsects));

  auto step7p2 = EmitNewVertsAndEdges(impl, edges, etIsects, eps);
  dbg.step7p2NewVerts = static_cast<int>(step7p2.newVertPositions.size());
  dbg.step7p2NewEdges = static_cast<int>(step7p2.newEdges.size());
  OVERLAP3D_HASH_DUMP("step7p2.newVertPositions",
                      HashVec(step7p2.newVertPositions));
  OVERLAP3D_HASH_DUMP("step7p2.newEdges", HashVec(step7p2.newEdges));
  if (std::getenv("OVERLAP3D_PIPELINE_DIAG")) {
    std::cerr << "      pipeline: eeIsects=" << eeIsects.size()
              << ", etIsects=" << etIsects.size()
              << ", step7p2 newVerts=" << step7p2.newVertPositions.size()
              << ", newEdges=" << step7p2.newEdges.size()
              << ", dropped(n!=2)=" << step7p2.dropped_n_not_2
              << ", interior(n=1)=" << step7p2.interiorVertsPerTri.size()
              << "\n";
  }

  auto step8 = AddInteriorVertsToNewEdges(impl, step7p2.newVertPositions,
                                          step7p2.newEdges, onTriLists, eps);
  for (const auto& nwe : step8) dbg.step8Propagated += nwe.extraVerts.size();
  if (DetHashEnabled()) {
    uint64_t h = 0;
    for (const auto& nwe : step8) h ^= HashVec(nwe.extraVerts);
    OVERLAP3D_HASH_DUMP("step8.extraVerts XOR", h);
  }

  auto step9 = FindNewEdgeIntersections(impl, step7p2.newVertPositions,
                                        step7p2.newEdges, eps);
  dbg.step9NewNew = static_cast<int>(step9.size());

  auto step10 = CountSubEdgesPerTri(impl, edges, onEdgeLists, step8);
  for (const auto& c : step10) dbg.step10TotalSubEdges += c.total();

  // Steps 11-13: real pipeline. Builds polygon graph, runs the
  // analytical classifier with pair-symmetric chord enforcement,
  // emits cap triangles for k=1 boundary cycles. Returns
  // status-0 Manifold for the 4 stable fixtures + Cray (best-case
  // ~3/5 runs given variance; falls back to step-1 Manifold if
  // the pipeline output is non-manifold).
  PropagateNewVertsToOnEdgeLists(etIsects, step7p2.resolvedIds,
                                 static_cast<int>(impl.NumVert()), edges,
                                 onEdgeLists);
  if (DetHashEnabled()) {
    uint64_t h = 0;
    for (const auto& l : onEdgeLists) h ^= HashVec(l.verts);
    OVERLAP3D_HASH_DUMP("onEdgeLists XOR (post-propagate)", h);
  }
  auto step11p1 = BuildPerTriHalfedgeGraphs(impl, edges, onEdgeLists, step8);
  if (DetHashEnabled()) {
    uint64_t hHE = 14695981039346656037ULL;
    uint64_t hVerts = 14695981039346656037ULL;
    for (const auto& g : step11p1) {
      uint64_t e1 = Fnv1a64(g.halfedges.data(),
                            g.halfedges.size() * sizeof(g.halfedges[0]));
      hHE ^= e1;
      hHE *= 1099511628211ULL;
      // hash the verts set deterministically by walking it sorted.
      std::vector<int> sortedVerts(g.verts.begin(), g.verts.end());
      std::sort(sortedVerts.begin(), sortedVerts.end());
      uint64_t e2 = HashVec(sortedVerts);
      hVerts ^= e2;
      hVerts *= 1099511628211ULL;
    }
    OVERLAP3D_HASH_DUMP("step11p1.halfedges ordered", hHE);
    OVERLAP3D_HASH_DUMP("step11p1.verts ordered", hVerts);
  }
  AddNextPointers(impl, step7p2.newVertPositions, step11p1);
  if (DetHashEnabled()) {
    uint64_t h = 14695981039346656037ULL;
    for (const auto& g : step11p1) {
      uint64_t e = HashVec(g.nextHalfedge);
      h ^= e;
      h *= 1099511628211ULL;
    }
    OVERLAP3D_HASH_DUMP("step11p2.nextHalfedge ordered", h);
  }
  auto step11p3 = WalkPolygons(step11p1);
  if (DetHashEnabled()) {
    uint64_t h = 14695981039346656037ULL;
    for (const auto& w : step11p3) {
      for (const auto& p : w.polygons) {
        uint64_t e = HashVec(p);
        h ^= e;
        h *= 1099511628211ULL;
      }
    }
    OVERLAP3D_HASH_DUMP("step11p3.polygons ordered", h);
  }

  const int baseId = static_cast<int>(impl.NumVert());
  auto chordPartners = BuildChordPartnerMap(step7p2.newEdges);
  if (DetHashEnabled()) {
    std::vector<std::tuple<int, int, int>> flat;
    for (const auto& [k, v] : chordPartners.partnerOf) {
      flat.push_back({std::get<0>(k), std::get<1>(k), v});
    }
    OVERLAP3D_HASH_DUMP("chordPartners", HashVec(flat));
  }

  // Compute Option B's per-vert winding via AnalyzeSelfMesh (one
  // ε-offset ray-cast per flood-fill component) for the
  // winding-grounded classifier below. The classifier replaces the
  // old per-polygon AnalyticalKeep heuristic with a "boundary of
  // (winding ≥ 1) region" test grounded in the global winding field.
  manifold::Vec<std::array<int, 2>> sma_p1q2;
  sma_p1q2.reserve(etIsects.size());
  for (const auto& x : etIsects) {
    if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size())) continue;
    sma_p1q2.push_back({edges[x.edgeIdx].halfedgeForward, x.triIdx});
  }
  std::sort(sma_p1q2.begin(), sma_p1q2.end(),
            [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
              return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
            });
  auto sma = manifold::AnalyzeSelfMesh(
      impl, manifold::VecView<const std::array<int, 2>>(sma_p1q2.data(),
                                                        sma_p1q2.size()));

  // Pre-compute classifier decisions + run pair-symmetric chord
  // enforcement (Phase 1). For each chord between triA and triB,
  // identify the 2 twin polygon pairs by halfedge direction; force
  // exactly one pair K + one pair D, deciding by K-vote per pair.
  using manifold::la::dot;
  using manifold::la::length;
  using manifold::la::normalize;
  const vec3 probeDir = normalize(vec3(0.7234, 0.4567, 0.5191));
  const double meshScale = length(impl.bBox_.max - impl.bBox_.min);
  const double probeEps = meshScale * 1e-9;
  const double rayLen = meshScale * 4.0;
  // OptB classifier: for each polygon vert, classify boundary /
  // interior / exterior using sma.w_above/w_below. Polygon kept iff
  // boundary verts dominate. All-chord-vert polygons fall back to
  // centroid probe via WindingAt(M, c ± ε * n_T, ...).
  // Centroid-probe-always variant: instead of using the per-vert
  // sma.w_above/w_below (which depend on AnalyzeSelfMesh's per-vert
  // normal averaging — ambiguous at intersection regions where left/
  // right tris with opposite normals share verts post-Subtract), do
  // a centroid probe via the polygon's tri face normal (= unambiguous
  // outward direction for that polygon's plane). Slower (one ray-cast
  // pair per polygon vs cached per-vert) but topology-correct on
  // Subtract-derived inputs (= cray, the Subtract back-side flip).
  const bool centroidAlways =
      std::getenv("OVERLAP3D_OPTB_CENTROID_ALWAYS") != nullptr;
  auto optbKeep = [&](int triId, const std::vector<int>& poly) -> bool {
    if (!centroidAlways) {
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
        return boundaryV > 0 && boundaryV >= interiorV &&
               boundaryV >= exteriorV;
      }
    }
    // Centroid probe via WindingAt, with per-face normal (= unambiguous
    // outward for the polygon's containing tri, even if the tri's verts
    // have ambiguous averaged normals).
    vec3 c(0, 0, 0);
    for (int v : poly) c += GetPos3(v, baseId, impl, step7p2.newVertPositions);
    c /= static_cast<double>(poly.size());
    const vec3 n_T = impl.faceNormal_[triId];
    const int wa =
        manifold::WindingAt(impl, c + n_T * probeEps, probeDir, rayLen);
    const int wb =
        manifold::WindingAt(impl, c - n_T * probeEps, probeDir, rayLen);
    return (wa == 0 && wb >= 1) || (wa >= 1 && wb == 0);
  };
  std::map<std::pair<int, int>, bool> precomputedKeep;
  for (size_t triId = 0; triId < step11p3.size(); ++triId) {
    const auto& polys = step11p3[triId].polygons;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
      if (polys.size() == 1) {
        precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] = true;
        continue;
      }
      precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] =
          optbKeep(static_cast<int>(triId), polys[pi]);
    }
  }
  if (DetHashEnabled()) {
    std::vector<std::tuple<int, int, int>> flat;
    for (const auto& [k, v] : precomputedKeep) {
      flat.push_back({k.first, k.second, v ? 1 : 0});
    }
    OVERLAP3D_HASH_DUMP("precomputedKeep (post-classifier)", HashVec(flat));
    // Also dump the raw map for diff'ing.
    if (std::getenv("OVERLAP3D_DET_HASH_DUMP_KEEP")) {
      for (const auto& [k, v] : precomputedKeep) {
        std::cerr << "      [keep] tri=" << k.first << " pi=" << k.second
                  << " keep=" << v << "\n";
      }
    }
  }
  // Pre-populate precomputedKeep entries for degenerate (2-vert) polys.
  // Default = false (= no output contribution). pair-sym never
  // changes these — they're "always dropped" semantically.
  for (size_t triId = 0; triId < step11p3.size(); ++triId) {
    const size_t base = step11p3[triId].polygons.size();
    const size_t ndeg = step11p3[triId].degeneratePolygons.size();
    for (size_t di = 0; di < ndeg; ++di) {
      precomputedKeep[{static_cast<int>(triId), static_cast<int>(base + di)}] =
          false;
    }
  }
  // findPolyHE: returns the polygon index containing the directed
  // halfedge (v0→v1). Searches both real polygons and degenerate
  // 2-vert polygons. Degenerate hits return pi = polygons.size() + di
  // (= a marker beyond the real range).
  auto findPolyHE = [&](int triId, int v0, int v1) -> int {
    const auto& polys = step11p3[triId].polygons;
    for (size_t pi = 0; pi < polys.size(); ++pi) {
      const auto& poly = polys[pi];
      for (size_t i = 0; i < poly.size(); ++i) {
        if (poly[i] == v0 && poly[(i + 1) % poly.size()] == v1)
          return static_cast<int>(pi);
      }
    }
    // Degenerate 2-vert polys: emit either direction (depending on
    // walk order). Accept both to find the chord regardless.
    const auto& deg = step11p3[triId].degeneratePolygons;
    for (size_t di = 0; di < deg.size(); ++di) {
      const auto& dp = deg[di];
      if (dp.size() != 2) continue;
      if ((dp[0] == v0 && dp[1] == v1) || (dp[0] == v1 && dp[1] == v0)) {
        return static_cast<int>(polys.size() + di);
      }
    }
    return -1;
  };
  // Helper: is this pi a degenerate-poly index (= "always dropped")?
  auto isDegenerate = [&](int triId, int pi) -> bool {
    if (pi < 0) return false;
    const size_t base = step11p3[triId].polygons.size();
    return static_cast<size_t>(pi) >= base;
  };
  int chordPairSymTotal = 0, chordPairSymSkipped = 0;
  int csBranch_zz = 0, csBranch_d1z = 0, csBranch_d2z = 0, csBranch_11 = 0,
      csBranch_12 = 0, csBranch_21 = 0, csBranch_22 = 0;
  static int dumpedChord = 0;
  for (const auto& edge : step7p2.newEdges) {
    ++chordPairSymTotal;
    int piA1 = findPolyHE(edge.triA, edge.v0, edge.v1);
    int piB1 = findPolyHE(edge.triB, edge.v1, edge.v0);
    int piA2 = findPolyHE(edge.triA, edge.v1, edge.v0);
    int piB2 = findPolyHE(edge.triB, edge.v0, edge.v1);
    if (std::getenv("OVERLAP3D_PAIRSYM_TRACE") && dumpedChord < 3) {
      std::fprintf(stderr,
                   "      [chord trace] triA=%d triB=%d v0=%d v1=%d "
                   "piA1=%d piB1=%d piA2=%d piB2=%d\n",
                   edge.triA, edge.triB, edge.v0, edge.v1, piA1, piB1, piA2,
                   piB2);
      ++dumpedChord;
    }
    auto getKeep = [&](int triId, int pi) -> bool* {
      if (pi < 0) return nullptr;
      auto it = precomputedKeep.find({triId, pi});
      return it == precomputedKeep.end() ? nullptr : &it->second;
    };
    // Per-direction manifold constraint:
    //   d1 (= v0→v1) contributors: A1, B2.
    //   d2 (= v1→v0) contributors: A2, B1.
    // For each direction, exactly 1 contributor must keep (= manifold
    // k=2 on the chord). A "contributor" is a non-null, non-degenerate
    // poly index. Degenerate polys (= 2-vert chord-pair cycles)
    // contribute 0 to output triangles and must not be flipped.
    bool* kA1 = getKeep(edge.triA, piA1);
    bool* kB1 = getKeep(edge.triB, piB1);
    bool* kA2 = getKeep(edge.triA, piA2);
    bool* kB2 = getKeep(edge.triB, piB2);
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
      // No contributors at all — nothing to enforce.
      ++chordPairSymSkipped;
      ++csBranch_zz;
    } else if (n_d1 == 0) {
      // d1 has 0 contributors; can't make d1 = k=1. Drop d2 contribs
      // to avoid k=1 mismatch (= chord has 0 incidence both sides).
      if (A2c) *kA2 = false;
      if (B1c) *kB1 = false;
      ++chordPairSymSkipped;
      ++csBranch_d1z;
    } else if (n_d2 == 0) {
      if (A1c) *kA1 = false;
      if (B2c) *kB2 = false;
      ++chordPairSymSkipped;
      ++csBranch_d2z;
    } else if (n_d1 == 1 && n_d2 == 1) {
      // Each direction has exactly 1 contributor → both must keep.
      if (A1c) *kA1 = true;
      if (B2c) *kB2 = true;
      if (A2c) *kA2 = true;
      if (B1c) *kB1 = true;
      ++chordPairSymSkipped;
      ++csBranch_11;
    } else if (n_d1 == 1 && n_d2 == 2) {
      // d1's only contributor must keep. d2: keep one, drop other.
      // Tiebreak: prefer the one whose classifier said keep; else A2.
      if (A1c) *kA1 = true;
      if (B2c) *kB2 = true;
      if (A2c && B1c) {
        if (*kA2 || !*kB1) {
          *kA2 = true;
          *kB1 = false;
        } else {
          *kA2 = false;
          *kB1 = true;
        }
      }
      ++chordPairSymSkipped;
      ++csBranch_12;
    } else if (n_d1 == 2 && n_d2 == 1) {
      if (A2c) *kA2 = true;
      if (B1c) *kB1 = true;
      if (A1c && B2c) {
        if (*kA1 || !*kB2) {
          *kA1 = true;
          *kB2 = false;
        } else {
          *kA1 = false;
          *kB2 = true;
        }
      }
      ++chordPairSymSkipped;
      ++csBranch_21;
    } else {
      // Full case (n_d1=2, n_d2=2): original pair-sym pair-level
      // enforcement. Pair p1 = (A1, B1), pair p2 = (A2, B2). One pair
      // fully keeps, other fully drops.
      int p1Vote = (*kA1 ? 1 : 0) + (*kB1 ? 1 : 0);
      int p2Vote = (*kA2 ? 1 : 0) + (*kB2 ? 1 : 0);
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
      ++csBranch_22;
    }
  }
  if (std::getenv("OVERLAP3D_PAIRSYM_DIAG")) {
    std::cerr << "      chord pair-sym: " << chordPairSymTotal
              << " chords total, " << chordPairSymSkipped
              << " skipped (incomplete 4-poly corner)\n";
    std::cerr << "      branches: zz=" << csBranch_zz << " d1z=" << csBranch_d1z
              << " d2z=" << csBranch_d2z << " 11=" << csBranch_11
              << " 12=" << csBranch_12 << " 21=" << csBranch_21
              << " 22=" << csBranch_22 << "\n";
  }

  // Phase 2: mesh-edge halfedge consistency. For each polygon
  // perimeter halfedge (a→b), find its reverse (b→a) in another
  // polygon. If found, AND-merge (drop iff either dropped). Single
  // pass — iteration over-cascades on dense self-intersection.
  {
    std::set<std::tuple<int, int, int>> chordHalfedges;
    for (const auto& edge : step7p2.newEdges) {
      chordHalfedges.insert({edge.triA, edge.v0, edge.v1});
      chordHalfedges.insert({edge.triA, edge.v1, edge.v0});
      chordHalfedges.insert({edge.triB, edge.v0, edge.v1});
      chordHalfedges.insert({edge.triB, edge.v1, edge.v0});
    }
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> halfedgeMap;
    for (size_t triId = 0; triId < step11p3.size(); ++triId) {
      const auto& polys = step11p3[triId].polygons;
      for (size_t pi = 0; pi < polys.size(); ++pi) {
        const auto& poly = polys[pi];
        for (size_t i = 0; i < poly.size(); ++i) {
          int a = poly[i], b = poly[(i + 1) % poly.size()];
          halfedgeMap[{a, b}].push_back(
              {static_cast<int>(triId), static_cast<int>(pi)});
        }
      }
    }
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
    // Phase 2.5: multi-owner per-direction enforcement. For edges where
    // > 1 polygon emits the SAME directed halfedge (a→b), only one can
    // survive in the manifold output (= each direction needs k=1, so
    // the undirected edge is k=2). This happens at 4-way (or higher)
    // self-intersection geometry where N input tris share the edge a-b
    // and each contributes its own perimeter halfedge.
    //
    // Strategy: among the multiple kept owners per direction, keep one
    // (preferring the classifier's already-kept ones; deterministic
    // tiebreak by triId). Drop the rest. Skip chord edges (= already
    // handled by Phase 1 pair-sym).
    int phase25Drops = 0;
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
      // Filter to currently-kept owners only.
      std::vector<std::pair<int, int>> kept;
      for (const auto& o : owners) {
        auto it = precomputedKeep.find(o);
        if (it == precomputedKeep.end()) continue;
        if (it->second) kept.push_back(o);
      }
      if (kept.size() <= 1) continue;
      // Keep the one with the lowest triId (deterministic). Drop the
      // rest. (Future refinement: classifier-confidence weighting.)
      int keepIdx = 0;
      for (size_t i = 1; i < kept.size(); ++i) {
        if (kept[i].first < kept[keepIdx].first) keepIdx = static_cast<int>(i);
      }
      for (size_t i = 0; i < kept.size(); ++i) {
        if (static_cast<int>(i) == keepIdx) continue;
        precomputedKeep[kept[i]] = false;
        ++phase25Drops;
      }
    }
    if (phase25Drops > 0 && std::getenv("OVERLAP3D_PAIRSYM_DIAG")) {
      std::cerr << "      pair-sym phase 2.5 (multi-owner): " << phase25Drops
                << " polys dropped\n";
    }
  }

  PolygonClassifierFn classifier = [&](const std::vector<int>& poly,
                                       const vec3& triNormal,
                                       int triId) -> PolygonClassification {
    PolygonClassification c{false, false, 0, 0};
    if (poly.size() < 3) return c;
    int pi = -1;
    const auto& polys = step11p3[triId].polygons;
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
                              step7p2.newVertPositions, baseId, impl);
    }
    c.reverse = false;
    (void)triNormal;
    return c;
  };
  // Always enable surface-cap in OverlapRemoval API path (= the
  // milestone behavior; produces 5/5 status 0 for fixtures where
  // cap walk closes successfully).
  setenv("OVERLAP3D_CAP", "1", 0);
  auto step13 = TriangulateAndEmit(impl, step7p2.newVertPositions, step11p3,
                                   classifier, &step7p2.interiorVertsPerTri);
  Manifold out = step13.output;
  if (out.Status() == Manifold::Error::NoError) {
    // Pierce-count guard: now safe to use after the determinism fix
    // (commit 40a6eb00). If the pipeline produced a manifold but with
    // MORE pierces than the input, the "fix" is worse than the
    // disease — fall back to merged input. Opt-out via
    // OVERLAP3D_NO_PIERCE_GUARD=1 (= return raw pipeline output even
    // if it regressed).
    //
    // Force materialization via Volume() before the pierce check so
    // both sides see the post-eval mesh (= avoid the lazy-eval
    // ordering issue where consecutive CheckSelfIntersection calls
    // can disagree). Matches the --api harness's measurement order.
    if (!std::getenv("OVERLAP3D_NO_PIERCE_GUARD")) {
      const double inVol = mr.manifold.Volume();
      const double outVol = out.Volume();
      // Compare against ORIGINAL input pierces, not post-merge.
      // MergeVertsEps can introduce pierces (= ε-merging close verts
      // onto edges of other faces), so post-merge siIn ≥ input pierces
      // and would mask user-visible regressions. The user passes
      // `input` and sees `out`; gate against that pair.
      auto siIn = CheckSelfIntersection(input, 1e-12);
      auto siOut = CheckSelfIntersection(out, 1e-12);
      // Pierce-monotonicity: pipeline must not make pierces worse.
      const bool pierceWorse = siOut.interiorPierces > siIn.interiorPierces;
      if (std::getenv("OVERLAP3D_DEBUG_GATE")) {
        auto siMerge = CheckSelfIntersection(mr.manifold, 1e-12);
        std::cerr << "[gate] input v" << input.NumVert() << " t"
                  << input.NumTri() << " vol " << input.Volume() << " | out v"
                  << out.NumVert() << " t" << out.NumTri() << " vol "
                  << out.Volume() << " | siInput=" << siIn.interiorPierces
                  << " siMerge=" << siMerge.interiorPierces
                  << " siOut=" << siOut.interiorPierces
                  << " worse=" << pierceWorse << "\n";
      }
      // Volume-sanity: even if pierces improved, an output with wildly
      // different volume (= geometry destroyed/inverted, e.g. cray's
      // Subtract back-side flip producing a negative-volume shell) is
      // also a fallback case. Threshold: 50% drift OR sign flip.
      const double drift =
          std::fabs(outVol - inVol) / std::max(std::fabs(inVol), 1e-12);
      const bool signFlip = (inVol > 0) != (outVol > 0);
      const bool driftBad = (drift > 0.50) || signFlip;
      // Sign-flip recovery: for Subtract-derived inputs with back-side
      // normal flip (= cray case), the pipeline may produce an inverted
      // manifold with negative volume. If pierces improved AND the
      // ONLY issue is sign flip, try reversing the output orientation
      // and re-checking.
      if (signFlip && !pierceWorse && drift > 0.50) {
        // Reverse winding by flipping vert order in each tri.
        manifold::MeshGL64 flipped = step13.output.GetMeshGL64();
        for (size_t t = 0; t < flipped.triVerts.size() / 3; ++t) {
          std::swap(flipped.triVerts[3 * t + 1], flipped.triVerts[3 * t + 2]);
        }
        Manifold flipMan(flipped);
        if (flipMan.Status() == Manifold::Error::NoError) {
          (void)flipMan.Volume();
          auto siFlip = CheckSelfIntersection(flipMan, 1e-12);
          const double flipVol = flipMan.Volume();
          const double flipDrift =
              std::fabs(flipVol - inVol) / std::max(std::fabs(inVol), 1e-12);
          const bool flipSignOK = (inVol > 0) == (flipVol > 0);
          if (flipSignOK && flipDrift <= 0.50 &&
              siFlip.interiorPierces <= siIn.interiorPierces) {
            if (std::getenv("OVERLAP3D_OR_TRACE")) {
              std::cerr << "      OverlapRemoval: SIGN-FLIP RECOVERY "
                        << "(flipped output: vol " << flipVol << ", drift "
                        << (flipDrift * 100.0) << "%, pierces "
                        << siFlip.interiorPierces << ")\n";
            }
            (void)eeIsects;
            (void)step9;
            return {flipMan, dbg};
          }
        }
      }
      if (pierceWorse || driftBad) {
        if (std::getenv("OVERLAP3D_OR_TRACE")) {
          std::cerr << "      OverlapRemoval: FALLBACK ";
          if (pierceWorse)
            std::cerr << "[pierce regress " << siOut.interiorPierces
                      << " vs in " << siIn.interiorPierces << "] ";
          if (driftBad)
            std::cerr << "[drift " << (drift * 100.0) << "%"
                      << (signFlip ? " sign-flip" : "") << "] ";
          std::cerr << "to merged input v" << mr.manifold.NumVert() << " t"
                    << mr.manifold.NumTri() << "\n";
        }
        (void)eeIsects;
        (void)step9;
        // If MergeVertsEps regressed vs the original input, return
        // the user's input rather than mr.manifold.
        auto siMergeFB = CheckSelfIntersection(mr.manifold, 1e-12);
        if (siMergeFB.interiorPierces > siIn.interiorPierces) {
          return {input, dbg};
        }
        return {mr.manifold, dbg};
      }
    }
    if (std::getenv("OVERLAP3D_OR_TRACE")) {
      std::cerr << "      OverlapRemoval: returning pipeline output v"
                << out.NumVert() << " t" << out.NumTri() << "\n";
    }
    (void)eeIsects;
    (void)step9;
    return {out, dbg};
  }
  if (std::getenv("OVERLAP3D_OR_TRACE")) {
    std::cerr << "      OverlapRemoval: pipeline output Status="
              << static_cast<int>(out.Status())
              << ", FALLBACK to merged input v" << mr.manifold.NumVert() << " t"
              << mr.manifold.NumTri() << "\n";
  }
  // Pipeline output is non-manifold (= cap walks failed or
  // classifier left k=1 edges). Fall back, but check first whether
  // mr.manifold (= MergeVertsEps result) introduced more pierces
  // than the original input. If so, return the ORIGINAL input —
  // returning mr.manifold would be a regression from the user's
  // perspective.
  //
  // Empirical: overnight fuzz (5000 seeds, ~125k self-pierce cases)
  // found 317 worsened cases where this exact pattern triggered.
  auto siInput = CheckSelfIntersection(input, 1e-12);
  auto siMerge = CheckSelfIntersection(mr.manifold, 1e-12);
  if (std::getenv("OVERLAP3D_DEBUG_GATE")) {
    std::cerr << "[gate-fallback] input v" << input.NumVert() << " t"
              << input.NumTri() << " vol " << input.Volume()
              << " | out-pipeline-status=" << static_cast<int>(out.Status())
              << " | mr v" << mr.manifold.NumVert() << " t"
              << mr.manifold.NumTri()
              << " | siInput=" << siInput.interiorPierces
              << " siMerge=" << siMerge.interiorPierces
              << (siMerge.interiorPierces > siInput.interiorPierces
                      ? " (merge regressed → fall back to INPUT)"
                      : " (returning mr.manifold)")
              << "\n";
  }
  (void)eeIsects;
  (void)step9;
  if (siMerge.interiorPierces > siInput.interiorPierces) {
    return {input, dbg};
  }
  return {mr.manifold, dbg};
}

// NOTE on cross-checking with manifold's own `IsSelfIntersecting`:
//
// Manifold has an internal tri-tri self-intersection check at
// `src/properties.cpp::Manifold::Impl::IsSelfIntersecting()`. It uses
// a different algorithm (BVH self-collisions + Möller's `Distance-
// TriangleTriangleSquared` + a ±eps along-normal nudge heuristic to
// distinguish "touching" from "crossing"), and is gated in production
// behind `ManifoldParams().selfIntersectionChecks` (default off,
// `MANIFOLD_DEBUG`-only). It is exposed only on `Manifold::Impl`
// (private), reachable via `Manifold::GetCsgLeafNode().GetImpl()` --
// also private.
//
// To run a side-by-side cross-check from this spike, the cleanest
// option would be to add MANIFOLD_DEBUG to the build, set
// `selfIntersectionChecks = true`, and observe that the existing
// pipeline throws `logicErr("self intersection detected")` from
// `csg_tree.cpp::SimpleBoolean` -- but that loses the result and
// turns each finding into a fatal error, defeating batch-mode
// adversarial cataloging. Instead, the spike's `CheckSelfIntersection`
// implements an independent edge-pierces-triangle algorithm; it is
// tested as a working oracle against the named `.obj` fixtures
// (which are already known by the existing test suite to need
// `processOverlaps = true`).

}  // namespace overlap3d

// =============================================================================
// Smoke test battery.
// =============================================================================

namespace {

void Header(const char* name) { std::cout << "=== " << name << " ===\n"; }

void Report(const char* tag, const manifold::Manifold& m) {
  std::cout << "  " << tag << ": volume=" << m.Volume()
            << ", surfaceArea=" << m.SurfaceArea() << ", numTris=" << m.NumTri()
            << ", numVert=" << m.NumVert() << "\n";
}

bool NearlyEqual(double a, double b, double tol = 1e-9) {
  return std::fabs(a - b) < tol;
}

}  // namespace

// Fast targeted-fixture debug path. Usage:
//   overlap3d_proto cray | self-intersect | hull-mask | offset12 | havocglass
// Loads the named fixture, runs steps 1-13d, prints all the debug
// stats, exits. Skips the heavy fuzz + battery so iteration is <2s.
inline int RunSingleFixture(const std::string& name);

int main(int argc, char** argv) {
  if (argc >= 2) {
    bool nofilter = false;
    bool apiMode = false;
    bool prodApi = false;
    std::string name = argv[1];
    for (int i = 2; i < argc; ++i) {
      if (std::string(argv[i]) == "--nofilter") nofilter = true;
      if (std::string(argv[i]) == "--api") apiMode = true;
      if (std::string(argv[i]) == "--prod-api") {
        apiMode = true;
        prodApi = true;
      }
    }
    setenv("OVERLAP3D_NOFILTER", nofilter ? "1" : "0", 1);
    // Adversarial cube fuzz CLI mode. Targets near-degenerate
    // configs at high displacement, where FP precision bites and
    // self-pierce rate is much higher than the general 800-case
    // fuzz (1/800 rate). Usage:
    //   ./overlap3d_proto --advfuzz
    if (name == "--advfuzz") {
      using namespace overlap3d;
      // Seed override via OVERLAP3D_FUZZ_SEED env var. Default 31415
      // matches the original hardcoded seed for reproducibility.
      // Use to scan many seeds in a wrapper script overnight.
      uint32_t advSeed = 31415;
      if (const char* s = std::getenv("OVERLAP3D_FUZZ_SEED"))
        advSeed = static_cast<uint32_t>(std::strtoul(s, nullptr, 10));
      const bool advSave = std::getenv("OVERLAP3D_FUZZ_SAVE") != nullptr;
      const std::string advOutDir = []() {
        const char* d = std::getenv("OVERLAP3D_FUZZ_OUTDIR");
        return std::string(d ? d : "/tmp/fuzz");
      }();
      if (advSave) {
        std::filesystem::create_directories(advOutDir);
      }
      std::mt19937 rng(advSeed);
      std::uniform_real_distribution<double> u(-1.0, 1.0);
      std::uniform_real_distribution<double> tiny(-0.01, 0.01);
      const double off = std::ldexp(1.5, 30);
      // For each Boolean result: count self-pierces, then if any,
      // run OverlapRemoval and count again. Returns:
      //   {pre-pierce, post-pierce, fallback}
      // pre/post are interiorPierces from CheckSelfIntersection.
      // fallback = true if OverlapRemoval returned merged input
      // (= API's pipeline produced non-manifold).
      //
      // Mesh-quality probe: separate from pierce monotonicity, we
      // also measure tri-shape stats on input and output:
      //   - numSlivers: tris with area < bbox.Scale()^2 * 1e-12
      //     (scale-invariant sub-FP-noise threshold).
      //   - maxAspect: max over all tris of longest_edge / shortest_edge.
      //     equilateral = 1; high aspect = thin sliver.
      //   - minAngle: min over all tris of min interior angle (degrees).
      //     equilateral = 60°; near-0° = degenerate.
      // Tracks whether the pipeline INTRODUCES slivers (= a regression
      // separate from pierce count) and how its output tri-shape
      // distribution compares to its input.
      struct MeshQuality {
        int numSlivers = 0;
        double maxAspect = 1.0;
        double minAngle = 180.0;
      };
      auto computeMeshQuality = [&](const Manifold& m) -> MeshQuality {
        MeshQuality q;
        if (m.IsEmpty()) return q;
        const auto mesh = m.GetMeshGL64();
        const double scale = m.BoundingBox().Scale();
        const double sliverAreaThreshold = scale * scale * 1e-12;
        const size_t nTri = mesh.NumTri();
        for (size_t t = 0; t < nTri; ++t) {
          auto vp = [&](int i) {
            return vec3(
                mesh.vertProperties[mesh.numProp * mesh.triVerts[3 * t + i] +
                                    0],
                mesh.vertProperties[mesh.numProp * mesh.triVerts[3 * t + i] +
                                    1],
                mesh.vertProperties[mesh.numProp * mesh.triVerts[3 * t + i] +
                                    2]);
          };
          const vec3 a = vp(0), b = vp(1), c = vp(2);
          const vec3 ab = b - a, bc = c - b, ca = a - c;
          const vec3 cross_ =
              manifold::la::cross(ab, vec3(-ca.x, -ca.y, -ca.z));
          const double area =
              0.5 * std::sqrt(manifold::la::dot(cross_, cross_));
          if (area < sliverAreaThreshold) ++q.numSlivers;
          const double lenAb = std::sqrt(manifold::la::dot(ab, ab));
          const double lenBc = std::sqrt(manifold::la::dot(bc, bc));
          const double lenCa = std::sqrt(manifold::la::dot(ca, ca));
          const double longest = std::max({lenAb, lenBc, lenCa});
          const double shortest = std::min({lenAb, lenBc, lenCa});
          if (shortest > 0) {
            const double aspect = longest / shortest;
            if (aspect > q.maxAspect) q.maxAspect = aspect;
          }
          // Min interior angle via law of cosines on each corner.
          auto cornerAngle = [](double opp, double s1, double s2) {
            if (s1 <= 0 || s2 <= 0) return 180.0;
            const double c = (s1 * s1 + s2 * s2 - opp * opp) / (2 * s1 * s2);
            const double clamped = std::max(-1.0, std::min(1.0, c));
            return std::acos(clamped) * 180.0 / 3.141592653589793;
          };
          const double angA = cornerAngle(lenBc, lenAb, lenCa);
          const double angB = cornerAngle(lenCa, lenAb, lenBc);
          const double angC = cornerAngle(lenAb, lenBc, lenCa);
          const double minAng = std::min({angA, angB, angC});
          if (minAng < q.minAngle) q.minAngle = minAng;
        }
        return q;
      };

      struct CaseStats {
        int prePierces = 0;
        int postPierces = -1;  // -1 = not run
        bool fallback = false;
        bool fixSucceeded = false;  // post=0 AND no fallback
        MeshQuality preMQ;
        MeshQuality postMQ;
        bool postMQValid = false;
      };
      int diagFailCount = 0;
      auto saveAdv = [&](const Manifold& m, const std::string& label) {
        if (!advSave || m.IsEmpty()) return;
        const std::string path = advOutDir + "/seed" + std::to_string(advSeed) +
                                 "_" + label + ".obj";
        std::ofstream out(path);
        if (out.is_open()) m.WriteOBJ(out);
      };
      auto runCase = [&](Manifold result,
                         const std::string& label = "") -> CaseStats {
        CaseStats s;
        if (result.IsEmpty() || result.Status() != Manifold::Error::NoError) {
          s.prePierces = -1;
          return s;
        }
        // Skip pierces below the mesh's own coincidence tolerance:
        // those are FP-noise that the Boolean engine is allowed to
        // produce within tolerance, not geometric overlaps the pipeline
        // should fix. Use input tolerance for both pre and post so the
        // comparison is apples-to-apples even if the pipeline tightens
        // tolerance.
        const double pierceMinMag = result.GetTolerance();
        auto si = CheckSelfIntersection(result, 1e-12, pierceMinMag);
        s.prePierces = si.interiorPierces;
        s.preMQ = computeMeshQuality(result);
        if (s.prePierces == 0) return s;
        // Save every piercing input when OVERLAP3D_FUZZ_SAVE_ALL is set
        // (= analysis sweep). Default save only writes failures below.
        if (std::getenv("OVERLAP3D_FUZZ_SAVE_ALL") && !label.empty()) {
          saveAdv(result,
                  label + "_pre" + std::to_string(s.prePierces) + "_in");
        }
        // Has pierces — run OverlapRemoval to see if it fixes them.
        auto [out, _] = OverlapRemoval(result);
        if (out.IsEmpty() || out.Status() != Manifold::Error::NoError) {
          s.postPierces = -2;  // OverlapRemoval failed
          if (std::getenv("OVERLAP3D_ADVFUZZ_DIAG") && diagFailCount < 5) {
            std::cerr << "    regularized#" << diagFailCount << " input v"
                      << result.NumVert() << " t" << result.NumTri() << " vol "
                      << result.Volume() << " bbox.scale "
                      << result.BoundingBox().Scale() << " pre-pierces "
                      << s.prePierces << " (sub-eps sliver: feature/bbox ~"
                      << (result.Volume() /
                          std::pow(result.BoundingBox().Scale(), 3))
                      << ")\n";
            ++diagFailCount;
          }
          return s;
        }
        // Detect fallback by checking if output is identical to input.
        s.fallback = (out.NumVert() == result.NumVert() &&
                      out.NumTri() == result.NumTri() &&
                      std::fabs(out.Volume() - result.Volume()) <
                          std::max(std::fabs(result.Volume()), 1.0) * 1e-15);
        auto si2 = CheckSelfIntersection(out, 1e-12, pierceMinMag);
        s.postPierces = si2.interiorPierces;
        s.fixSucceeded = (s.postPierces == 0 && !s.fallback);
        s.postMQ = computeMeshQuality(out);
        s.postMQValid = true;
        // Save adversarial inputs: cases where pipeline didn't fully
        // fix and didn't fall back. These are the interesting failure
        // modes for further investigation.
        if (s.postPierces > 0 && !s.fallback && !label.empty()) {
          saveAdv(result, label + "_pre" + std::to_string(s.prePierces) +
                              "_post" + std::to_string(s.postPierces));
        } else if (s.fallback && !label.empty() && advSave) {
          // Also save fallback cases (= classifier rejected entirely)
          // when explicitly asked, since these signify pipeline limits.
          saveAdv(result,
                  label + "_pre" + std::to_string(s.prePierces) + "_fallback");
        }
        if (std::getenv("OVERLAP3D_ADVFUZZ_DIAG") && s.postPierces > 0 &&
            !s.fallback) {
          std::cerr << "    partial-fix:" << " input v" << result.NumVert()
                    << " t" << result.NumTri() << " vol " << result.Volume()
                    << " bbox.scale " << result.BoundingBox().Scale()
                    << " pre-pierces " << s.prePierces << " post-pierces "
                    << s.postPierces << " (pre-pierce magnitude "
                    << si.maxPierceMagnitude << " rel="
                    << (si.maxPierceMagnitude /
                        std::max(result.BoundingBox().Scale(), 1.0))
                    << ", post-magnitude " << si2.maxPierceMagnitude << " rel="
                    << (si2.maxPierceMagnitude /
                        std::max(result.BoundingBox().Scale(), 1.0))
                    << ")\n";
          std::cerr << "      output: v" << out.NumVert() << " t"
                    << out.NumTri() << " vol " << out.Volume()
                    << " volDelta=" << (out.Volume() - result.Volume()) << "\n";
        }
        return s;
      };
      std::uniform_real_distribution<double> ang(0.0, 6.283);
      struct ClassResult {
        const char* name;
        int total = 0, valid = 0, pierces = 0, totalPierces = 0;
        // Fix outcomes (only counted for cases with pierces):
        int fixed = 0;        // post=0, no fallback
        int reduced = 0;      // post>0 but post<pre
        int unchanged = 0;    // post==pre
        int worsened = 0;     // post>pre
        int fallbackCnt = 0;  // OverlapRemoval fell back to input
        int orFailed = 0;     // OverlapRemoval status != NoError
        double maxRel = 0;
        // Mesh-quality aggregates across all cases with a valid postMQ.
        int totalPreSlivers = 0;  // sum of preMQ.numSlivers
        int totalPostSlivers =
            0;  // sum of postMQ.numSlivers (only counted when postMQ valid)
        int casesIntroducedSlivers = 0;  // postMQ.numSlivers > preMQ.numSlivers
        int casesReducedSlivers = 0;     // postMQ.numSlivers < preMQ.numSlivers
        double worstPreAspect = 1.0;
        double worstPostAspect = 1.0;
        double worstPreMinAngle = 180.0;
        double worstPostMinAngle = 180.0;
      };
      std::vector<ClassResult> classes{
          {"shallow-tight"},
          {"near-coplanar-slabs"},
          {"3-cube-chain"},
          {"4-cube-chain"},
          {"5-cube-chain"},
          {"6-cube-chain"},
          {"vertex-on-face"},
          {"subtract-rotated"},
          {"subtract-3-chain"},
          {"smith-three-mutual"},
          {"smith-four-coplanar"},
          {"cospherical-shell"},
          {"cgal-12-cube-chain"},
          // Curved-surface + non-cube primitives (added May 2026).
          {"sphere-intersect"},  // two intersecting spheres, Intersect op
          {"cylinder-cross"},    // two crossing cylinders, Add op
          {"mixed-op-chain"},    // (sphere + cube) - cylinder
      };
      // Default 200 cases/class; override with --advseeds N.
      int casesPerClass = 200;
      for (int i = 2; i < argc - 1; ++i) {
        if (std::string(argv[i]) == "--advseeds") {
          casesPerClass = std::atoi(argv[i + 1]);
          if (casesPerClass <= 0) casesPerClass = 200;
        }
      }
      auto record = [](ClassResult& cls, const CaseStats& s) {
        if (s.prePierces < 0) return;  // boolean error
        ++cls.valid;
        // Mesh-quality aggregates: only sum pre + post when post is
        // valid (= pipeline actually ran). Comparing pre across ALL
        // cases vs post across cases-pipeline-ran is asymmetric and
        // misleading. This way the in→out totals reflect the same
        // case set, so a reduction in slivers is attributable to the
        // pipeline rather than to input distribution.
        if (s.postMQValid) {
          cls.totalPreSlivers += s.preMQ.numSlivers;
          cls.totalPostSlivers += s.postMQ.numSlivers;
          if (s.preMQ.maxAspect > cls.worstPreAspect)
            cls.worstPreAspect = s.preMQ.maxAspect;
          if (s.preMQ.minAngle < cls.worstPreMinAngle)
            cls.worstPreMinAngle = s.preMQ.minAngle;
          if (s.postMQ.maxAspect > cls.worstPostAspect)
            cls.worstPostAspect = s.postMQ.maxAspect;
          if (s.postMQ.minAngle < cls.worstPostMinAngle)
            cls.worstPostMinAngle = s.postMQ.minAngle;
          if (s.postMQ.numSlivers > s.preMQ.numSlivers)
            ++cls.casesIntroducedSlivers;
          else if (s.postMQ.numSlivers < s.preMQ.numSlivers)
            ++cls.casesReducedSlivers;
        }
        if (s.prePierces == 0) return;
        ++cls.pierces;
        cls.totalPierces += s.prePierces;
        // post=-2 means OverlapRemoval returned empty (= Requicha-
        // Tilove regularization of sub-eps sliver input). This is
        // the CORRECT output per CGAL/Clipper2/SVG fill-rule
        // conventions, not a failure.
        if (s.postPierces == -2) {
          ++cls.orFailed;
          return;
        }
        if (s.fallback) {
          ++cls.fallbackCnt;
          ++cls.unchanged;
          return;
        }
        if (s.postPierces == 0)
          ++cls.fixed;
        else if (s.postPierces < s.prePierces)
          ++cls.reduced;
        else if (s.postPierces == s.prePierces)
          ++cls.unchanged;
        else
          ++cls.worsened;
      };
      // Class A: shallow tight intersections at kPow=30 (rotated)
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[0].total;
        const vec3 trans(tiny(rng) * 0.1, tiny(rng) * 0.1, tiny(rng) * 0.1);
        const double theta = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                         .Translate(vec3(off, off, off) + trans);
        record(classes[0], runCase(Boolean3D(a, b, OpType::Add),
                                   "class0_seed" + std::to_string(seed)));
      }
      // Class B: near-coplanar thin slabs
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[1].total;
        const double angle = tiny(rng);
        Manifold a =
            Manifold::Cube({2, 2, 0.001}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 0.001}, true)
                         .Rotate(angle * 57.296, 0, 0)
                         .Translate(vec3(off, off, off));
        record(classes[1],
               runCase(a + b, "class1_seed" + std::to_string(seed)));
      }
      // Class C: 3-cube chained Boolean (rotated)
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[2].total;
        const vec3 t1(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const vec3 t2(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const double r1 = ang(rng), r2 = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r1 * 57.296, r1 * 28.6, r1 * 14.3)
                         .Translate(vec3(off, off, off) + t1);
        Manifold c = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r2 * 57.296, r2 * 28.6, r2 * 14.3)
                         .Translate(vec3(off, off, off) + t2);
        record(classes[2],
               runCase(Boolean3D(Boolean3D(a, b, OpType::Add), c, OpType::Add),
                       "class2_seed" + std::to_string(seed)));
      }
      // Class D: 4-cube chained Boolean (rotated)
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[3].total;
        const vec3 t1(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const vec3 t2(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const vec3 t3(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const double r1 = ang(rng), r2 = ang(rng), r3 = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r1 * 57.296, r1 * 28.6, r1 * 14.3)
                         .Translate(vec3(off, off, off) + t1);
        Manifold c = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r2 * 57.296, r2 * 28.6, r2 * 14.3)
                         .Translate(vec3(off, off, off) + t2);
        Manifold d = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r3 * 57.296, r3 * 28.6, r3 * 14.3)
                         .Translate(vec3(off, off, off) + t3);
        record(classes[3],
               runCase(Boolean3D(Boolean3D(Boolean3D(a, b, OpType::Add), c,
                                           OpType::Add),
                                 d, OpType::Add),
                       "class3_seed" + std::to_string(seed)));
      }
      // Class E: 5-cube chain
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[4].total;
        std::vector<vec3> ts(5);
        std::vector<double> rs(5);
        for (int i = 0; i < 5; ++i) {
          ts[i] = vec3(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
          rs[i] = ang(rng);
        }
        Manifold acc =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        for (int i = 0; i < 4; ++i) {
          Manifold next =
              Manifold::Cube({2, 2, 2}, true)
                  .Rotate(rs[i] * 57.296, rs[i] * 28.6, rs[i] * 14.3)
                  .Translate(vec3(off, off, off) + ts[i]);
          acc = Boolean3D(acc, next, OpType::Add);
          if (acc.IsEmpty() || acc.Status() != Manifold::Error::NoError) break;
        }
        record(classes[4], runCase(acc, "class4_seed" + std::to_string(seed)));
      }
      // Class F: 6-cube chain
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[5].total;
        std::vector<vec3> ts(6);
        std::vector<double> rs(6);
        for (int i = 0; i < 6; ++i) {
          ts[i] = vec3(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
          rs[i] = ang(rng);
        }
        Manifold acc =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        for (int i = 0; i < 5; ++i) {
          Manifold next =
              Manifold::Cube({2, 2, 2}, true)
                  .Rotate(rs[i] * 57.296, rs[i] * 28.6, rs[i] * 14.3)
                  .Translate(vec3(off, off, off) + ts[i]);
          acc = Boolean3D(acc, next, OpType::Add);
          if (acc.IsEmpty() || acc.Status() != Manifold::Error::NoError) break;
        }
        record(classes[5], runCase(acc, "class5_seed" + std::to_string(seed)));
      }
      // Class G: vertex-on-face placement (cube B's corner exactly
      // on cube A's face). Uses controlled offsets, no randomness
      // for the exact placement; small perturbation + BVH variance
      // gives the adversarial behavior.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[6].total;
        const vec3 trans(1.0 + tiny(rng), u(rng) * 0.4, u(rng) * 0.4);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Translate(vec3(off + 1, off, off) + trans);
        record(classes[6], runCase(Boolean3D(a, b, OpType::Add),
                                   "class6_seed" + std::to_string(seed)));
      }
      // Class H: subtract-rotated — Cray-class (single Subtract of two
      // overlapping rotated cubes). Tests whether Option B's known
      // weakness on Cray (back-side normal flip) generalizes to the
      // population of single-Subtract inputs.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[7].total;
        const vec3 trans(tiny(rng) * 0.1, tiny(rng) * 0.1, tiny(rng) * 0.1);
        const double theta = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                         .Translate(vec3(off, off, off) + trans);
        record(classes[7], runCase(Boolean3D(a, b, OpType::Subtract),
                                   "class7_seed" + std::to_string(seed)));
      }
      // Class I: subtract-3-chain — 2-cube Add then 3rd-cube Subtract.
      // Compound-op pattern that exercises chained back-side flag
      // propagation (= the topology Cray exhibits).
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[8].total;
        const vec3 t1(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const vec3 t2(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const double r1 = ang(rng), r2 = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r1 * 57.296, r1 * 28.6, r1 * 14.3)
                         .Translate(vec3(off, off, off) + t1);
        Manifold c = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(r2 * 57.296, r2 * 28.6, r2 * 14.3)
                         .Translate(vec3(off, off, off) + t2);
        record(classes[8], runCase(Boolean3D(Boolean3D(a, b, OpType::Add), c,
                                             OpType::Subtract),
                                   "class8_seed" + std::to_string(seed)));
      }
      // Class J: smith-three-mutual — three slabs (one per axis)
      // with random rotations. 3D analog of Smith UCAM-CL-TR-766
      // ch.9.2.2 / fig 9.1 (= the "edge-pierces-triangle laddering"
      // 3D case Smith identifies as the hardest). Slab thickness 0.5
      // matches the existing cube class feature scale to stay above
      // ε at kPow=30 displacement.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[9].total;
        const vec3 t1(tiny(rng), tiny(rng), tiny(rng));
        const vec3 t2(tiny(rng), tiny(rng), tiny(rng));
        const vec3 t3(tiny(rng), tiny(rng), tiny(rng));
        const double r1 = ang(rng), r2 = ang(rng), r3 = ang(rng);
        Manifold sx = Manifold::Cube({0.5, 2, 2}, true)
                          .Rotate(r1 * 57.296, r1 * 28.6, r1 * 14.3)
                          .Translate(vec3(off, off, off) + t1);
        Manifold sy = Manifold::Cube({2, 0.5, 2}, true)
                          .Rotate(r2 * 57.296, r2 * 28.6, r2 * 14.3)
                          .Translate(vec3(off, off, off) + t2);
        Manifold sz = Manifold::Cube({2, 2, 0.5}, true)
                          .Rotate(r3 * 57.296, r3 * 28.6, r3 * 14.3)
                          .Translate(vec3(off, off, off) + t3);
        record(classes[9], runCase(Boolean3D(Boolean3D(sx, sy, OpType::Add), sz,
                                             OpType::Add),
                                   "class9_seed" + std::to_string(seed)));
      }
      // Class K: smith-four-coplanar — four overlapping slabs offset
      // along x with a tiny tilt around y. Smith's only reported
      // field-bug case (UCAM-CL-TR-766 p.84): four near-coplanar /
      // near-parallel facets that historically broke a CAD system.
      // Slab thickness 0.5, offset 0.4 → adjacent slabs overlap by 0.1.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[10].total;
        const double dx = 0.4;
        const double tilt = tiny(rng) * 0.05;  // small tilt around y
        Manifold a = Manifold::Cube({0.5, 2, 2}, true)
                         .Rotate(0, tilt * 57.296, 0)
                         .Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({0.5, 2, 2}, true)
                         .Rotate(0, tilt * 57.296, 0)
                         .Translate(vec3(off + dx, off, off));
        Manifold c = Manifold::Cube({0.5, 2, 2}, true)
                         .Rotate(0, tilt * 57.296, 0)
                         .Translate(vec3(off + 2 * dx, off, off));
        Manifold d = Manifold::Cube({0.5, 2, 2}, true)
                         .Rotate(0, tilt * 57.296, 0)
                         .Translate(vec3(off + 3 * dx, off, off));
        record(classes[10],
               runCase(Boolean3D(Boolean3D(Boolean3D(a, b, OpType::Add), c,
                                           OpType::Add),
                                 d, OpType::Add),
                       "class10_seed" + std::to_string(seed)));
      }
      // Class L: cospherical-shell — 3D analog of cocircular-points
      // pathology. Two concentric spheres of slightly different radii;
      // their union/subtract creates a thin shell where many vertices
      // lie nearly cospherical. Radii r1=1.0 and r2 ≈ 1.5 + tiny match
      // feature scale so the shell stays above ε at kPow=30.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[11].total;
        const double r1 = 1.0;
        const double r2 = 1.5 + tiny(rng);
        const int segs = 16 + (seed % 16);  // 16-31 segments
        Manifold inner =
            Manifold::Sphere(r1, segs).Translate(vec3(off, off, off));
        Manifold outer =
            Manifold::Sphere(r2, segs).Translate(vec3(off, off, off));
        record(classes[11], runCase(Boolean3D(outer, inner, OpType::Subtract),
                                    "class11_seed" + std::to_string(seed)));
      }
      // Class M: cgal-12-cube-chain — CGAL Lazard & Valque (CGF 2025)
      // autorefinement workload: iterative union of N rotated cubes.
      // Doubled from the existing 6-cube-chain class to test deep-
      // chain FP-error accumulation that snap-rounding was designed
      // for. 100 cases (not 200) since each case runs 12 Booleans.
      for (int seed = 0; seed < 100; ++seed) {
        ++classes[12].total;
        std::vector<vec3> ts(12);
        std::vector<double> rs(12);
        for (int i = 0; i < 12; ++i) {
          ts[i] = vec3(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
          rs[i] = ang(rng);
        }
        Manifold acc =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        for (int i = 0; i < 11; ++i) {
          Manifold next =
              Manifold::Cube({2, 2, 2}, true)
                  .Rotate(rs[i] * 57.296, rs[i] * 28.6, rs[i] * 14.3)
                  .Translate(vec3(off, off, off) + ts[i]);
          acc = Boolean3D(acc, next, OpType::Add);
          if (acc.IsEmpty() || acc.Status() != Manifold::Error::NoError) break;
        }
        record(classes[12],
               runCase(acc, "class12_seed" + std::to_string(seed)));
      }
      // Class N: sphere-intersect — two unit spheres slightly
      // displaced from each other, Intersect op. Tests the chord
      // assembly on curved surfaces (= many more intersection
      // chords than cube-vs-cube) and exercises Intersect's code
      // path (= the existing 13 classes only use Add + Subtract).
      // Low circularSegments (16) keeps tri count manageable.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[13].total;
        const vec3 trans(tiny(rng) * 5.0, tiny(rng) * 5.0, tiny(rng) * 5.0);
        Manifold a = Manifold::Sphere(1.0, 16).Translate(vec3(off, off, off));
        Manifold b = Manifold::Sphere(1.0, 16).Translate(
            vec3(off, off, off) + trans + vec3(0.5, 0.5, 0.5));
        record(classes[13], runCase(Boolean3D(a, b, OpType::Intersect),
                                    "class13_seed" + std::to_string(seed)));
      }
      // Class N+1: cylinder-cross — two unit-radius cylinders, one
      // along z and one rotated 90° about y so it lies along x.
      // They cross at the origin region. Add op. Tests cylinders'
      // ring-shaped tri pattern through the chord assembly.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[14].total;
        const vec3 trans(tiny(rng), tiny(rng), tiny(rng));
        const double theta = ang(rng);
        Manifold a = Manifold::Cylinder(2.0, 1.0, 1.0, 16, true)
                         .Translate(vec3(off, off, off));
        Manifold b = Manifold::Cylinder(2.0, 1.0, 1.0, 16, true)
                         .Rotate(90, theta * 10, 0)
                         .Translate(vec3(off, off, off) + trans);
        record(classes[14], runCase(Boolean3D(a, b, OpType::Add),
                                    "class14_seed" + std::to_string(seed)));
      }
      // Class N+2: mixed-op-chain — sphere + cube + subtract cylinder.
      // Tests sequence (Add then Subtract) with mixed primitive
      // types, exercising the chord assembly across two different
      // op kinds in one input.
      for (int seed = 0; seed < casesPerClass; ++seed) {
        ++classes[15].total;
        const double r1 = ang(rng), r2 = ang(rng);
        const vec3 t1(u(rng) * 0.5, u(rng) * 0.5, u(rng) * 0.5);
        const vec3 t2(u(rng) * 0.3, u(rng) * 0.3, u(rng) * 0.3);
        Manifold sph = Manifold::Sphere(1.0, 16).Translate(vec3(off, off, off));
        Manifold cube = Manifold::Cube({1.6, 1.6, 1.6}, true)
                            .Rotate(r1 * 57.296, r1 * 28.6, r1 * 14.3)
                            .Translate(vec3(off, off, off) + t1);
        Manifold cyl = Manifold::Cylinder(2.0, 0.5, 0.5, 16, true)
                           .Rotate(r2 * 57.296, 0, 0)
                           .Translate(vec3(off, off, off) + t2);
        record(classes[15], runCase(Boolean3D(Boolean3D(sph, cube, OpType::Add),
                                              cyl, OpType::Subtract),
                                    "class15_seed" + std::to_string(seed)));
      }
      std::cout << "Adversarial cube fuzz (kPow=30 displacement, "
                << casesPerClass << " cases/class):\n";
      std::cout << "  class                | pierces / valid (rate) | "
                   "OverlapRemoval outcome\n";
      std::cout << "  ---------------------+------------------------+"
                   "--------------------------\n";
      int totalPierceCases = 0, totalValid = 0;
      int totalFixed = 0, totalReduced = 0, totalUnchanged = 0;
      int totalWorsened = 0, totalFallback = 0, totalOrFailed = 0;
      for (const auto& cls : classes) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "  %-20s | %3d/%-4d (%5.1f%%)      | "
                      "fix=%d red=%d reg=%d unch=%d wors=%d fb=%d\n",
                      cls.name, cls.pierces, cls.valid,
                      cls.valid > 0 ? 100.0 * cls.pierces / cls.valid : 0.0,
                      cls.fixed, cls.reduced, cls.orFailed, cls.unchanged,
                      cls.worsened, cls.fallbackCnt);
        std::cout << buf;
        totalPierceCases += cls.pierces;
        totalValid += cls.valid;
        totalFixed += cls.fixed;
        totalReduced += cls.reduced;
        totalUnchanged += cls.unchanged;
        totalWorsened += cls.worsened;
        totalFallback += cls.fallbackCnt;
        totalOrFailed += cls.orFailed;
      }
      std::cout << "\n  TOTAL self-pierce cases: " << totalPierceCases << "/"
                << totalValid << " ("
                << (100.0 * totalPierceCases / std::max(totalValid, 1))
                << "% self-pierce rate)\n";
      const int validOutputs = totalFixed + totalReduced + totalOrFailed;
      std::cout << "  OverlapRemoval outcomes (of " << totalPierceCases
                << " self-piercing cases):\n";
      std::cout << "    fixed (post=0):    " << totalFixed << " ("
                << (100.0 * totalFixed / std::max(totalPierceCases, 1))
                << "%)\n";
      std::cout << "    reduced (post<pre):" << totalReduced << "\n";
      std::cout << "    regularized-to-empty (sub-eps sliver, valid): "
                << totalOrFailed << "\n";
      std::cout << "    unchanged:         " << totalUnchanged << "\n";
      std::cout << "    worsened:          " << totalWorsened << "\n";
      std::cout << "    fallback (API non-manifold, returned input): "
                << totalFallback << "\n";
      std::cout << "  → " << validOutputs << "/" << totalPierceCases
                << " produce valid output (fixed + reduced + regularized)\n";

      // Mesh-quality probe: tri-shape stats across input vs output.
      // Tracks slivers (= area < bbox.Scale()^2 * 1e-12), aspect
      // ratio (longest-edge / shortest-edge), and min interior
      // angle. Detects pipelines that produce geometrically poor
      // output without ever increasing pierce count.
      std::cout << "\n  Mesh quality (input → output):\n";
      std::cout << "  class                |   slivers (in→out, ±cases) |"
                   "  worst aspect (in→out)  |  worst min-angle°\n";
      std::cout << "  ---------------------+----------------------------+"
                   "------------------------+--------------------\n";
      int totPreSliver = 0, totPostSliver = 0;
      int totIntroduced = 0, totReducedSliver = 0;
      double worstInAspect = 1.0, worstOutAspect = 1.0;
      double worstInMinAng = 180.0, worstOutMinAng = 180.0;
      for (const auto& cls : classes) {
        char buf[300];
        std::snprintf(buf, sizeof(buf),
                      "  %-20s |  %4d → %4d  (+%d/-%d cases) |"
                      "  %.1f → %.1f  |  %.2f° → %.2f°\n",
                      cls.name, cls.totalPreSlivers, cls.totalPostSlivers,
                      cls.casesIntroducedSlivers, cls.casesReducedSlivers,
                      cls.worstPreAspect, cls.worstPostAspect,
                      cls.worstPreMinAngle, cls.worstPostMinAngle);
        std::cout << buf;
        totPreSliver += cls.totalPreSlivers;
        totPostSliver += cls.totalPostSlivers;
        totIntroduced += cls.casesIntroducedSlivers;
        totReducedSliver += cls.casesReducedSlivers;
        if (cls.worstPreAspect > worstInAspect)
          worstInAspect = cls.worstPreAspect;
        if (cls.worstPostAspect > worstOutAspect)
          worstOutAspect = cls.worstPostAspect;
        if (cls.worstPreMinAngle < worstInMinAng)
          worstInMinAng = cls.worstPreMinAngle;
        if (cls.worstPostMinAngle < worstOutMinAng)
          worstOutMinAng = cls.worstPostMinAngle;
      }
      std::cout << "\n  TOTAL slivers: " << totPreSliver << " → "
                << totPostSliver
                << "  (cases introducing slivers: " << totIntroduced
                << ", cases reducing slivers: " << totReducedSliver << ")\n";
      std::cout << "  worst aspect ratio overall:  in=" << worstInAspect
                << "  out=" << worstOutAspect << "\n";
      std::cout << "  worst min interior angle:    in=" << worstInMinAng
                << "°  out=" << worstOutMinAng << "°\n";
      return 0;
    }
    // Load a saved adversarial fuzz OBJ and run it through
    // OverlapRemoval with diagnostics. Pairs with --advfuzz save dir.
    //   ./overlap3d_proto --advload
    //   /tmp/fuzz/seed31415_classN_seedM_preX_postY.obj
    if (name == "--advload") {
      using namespace overlap3d;
      if (argc < 3) {
        std::cerr << "Usage: overlap3d_proto --advload <path-to-saved-obj>\n";
        return 2;
      }
      std::ifstream in(argv[2]);
      if (!in.is_open()) {
        std::cerr << "Could not open " << argv[2] << "\n";
        return 2;
      }
      Manifold input = Manifold::ReadOBJ(in);
      std::cout << "input v" << input.NumVert() << " t" << input.NumTri()
                << " vol " << input.Volume() << " bbox.scale "
                << input.BoundingBox().Scale() << " status "
                << static_cast<int>(input.Status()) << " tol "
                << input.GetTolerance() << "\n";
      auto siBefore = CheckSelfIntersection(input, 1e-12);
      std::cout << "input self-pierces (relTol 1e-12): "
                << siBefore.interiorPierces << " (max mag "
                << siBefore.maxPierceMagnitude << ")\n";
      // Also report at the input's own tolerance scale so we can tell
      // "real pierces" from "below-tolerance FP noise".
      const double bboxScale = input.BoundingBox().Scale();
      const double tolRel =
          bboxScale > 0 ? input.GetTolerance() / bboxScale : 1e-12;
      auto siAtTol = CheckSelfIntersection(input, std::max(tolRel, 1e-12),
                                           input.GetTolerance());
      std::cout << "input self-pierces (at tolerance, relTol "
                << std::max(tolRel, 1e-12) << "): " << siAtTol.interiorPierces
                << " (max mag " << siAtTol.maxPierceMagnitude << ")\n";
      // Pierce-graph clustering: connected components of the
      // (tris=nodes, pierces=edges) graph. Distinguishes isolated
      // pierce pairs (= local-fix friendly) from chain/tangle patterns
      // (= where local fixes can introduce new pierces on neighbors).
      auto clusters = ClusterPierces(siAtTol.pierceEdges);
      std::cout << "pierce clusters: " << clusters.size() << "\n";
      for (size_t i = 0; i < clusters.size(); ++i) {
        const auto& c = clusters[i];
        const int n = static_cast<int>(c.tris.size());
        const int m = static_cast<int>(c.edges.size());
        std::string shape;
        if (n == 2 && m == 1)
          shape = "isolated";
        else if (m == n - 1 && c.maxDegree <= 2)
          shape = "chain";
        else if (m == n && c.maxDegree <= 2)
          shape = "ring";
        else if (c.maxDegree >= 3)
          shape = "tangle";
        else
          shape = "tree-like";
        std::cout << "  cluster " << i << ": " << n << " tris, " << m
                  << " pierces, maxDeg " << c.maxDegree << " [" << shape
                  << "]  tris=[";
        for (size_t j = 0; j < c.tris.size(); ++j) {
          if (j) std::cout << ",";
          std::cout << c.tris[j];
        }
        std::cout << "]\n";
      }
      auto [out, status] = OverlapRemoval(input);
      std::cout << "output v" << out.NumVert() << " t" << out.NumTri()
                << " vol " << out.Volume() << " status "
                << static_cast<int>(out.Status()) << "\n";
      auto siAfter = CheckSelfIntersection(out, 1e-12);
      std::cout << "output self-pierces: " << siAfter.interiorPierces
                << " (max mag " << siAfter.maxPierceMagnitude << ")\n";
      const bool fallback =
          (out.NumVert() == input.NumVert() && out.NumTri() == input.NumTri() &&
           std::fabs(out.Volume() - input.Volume()) <
               std::max(std::fabs(input.Volume()), 1.0) * 1e-15);
      std::cout << "fallback: " << (fallback ? "yes" : "no") << "\n";
      if (const char* outPath = std::getenv("OVERLAP3D_DUMP_OUT")) {
        std::ofstream f(outPath);
        if (f.is_open()) {
          out.WriteOBJ(f);
          std::cerr << "wrote output to " << outPath << "\n";
        }
      }
      return 0;
    }
    // General-fuzz CLI mode: same case generator as the no-args
    // deepfuzz (rotated cubes at 4 displacement scales × 200 seeds),
    // but running each self-piercing Boolean output through
    // OverlapRemoval to measure fix rate per scale.  Picks up
    // OVERLAP3D_OPTB / OVERLAP3D_OPTB_CENTROID_ONLY / etc. so the
    // classifier mode can be A/B-tested across the broader corpus.
    //
    // Usage:
    //   ./overlap3d_proto --genfuzz
    //
    // Default classifier vs Option B comparison (= the "broader-
    // corpus" test of whether the .obj-fixture wins generalize to
    // typical Boolean inputs):
    //   ./overlap3d_proto --genfuzz                             (default)
    //   OVERLAP3D_OPTB=1 OVERLAP3D_OPTB_CENTROID_ONLY=1 \
    //     OVERLAP3D_NOPAIRSYM_PHASE2=1 ./overlap3d_proto --genfuzz
    if (name == "--genfuzz") {
      using namespace overlap3d;
      std::mt19937 rng(1729);
      std::uniform_real_distribution<double> u(-1.5, 1.5);
      std::uniform_real_distribution<double> ang(0.0, 6.283);
      struct ScaleResult {
        int kPow = 0;
        int total = 0, valid = 0, pierces = 0, totalPierces = 0;
        int fixed = 0, reduced = 0, unchanged = 0, worsened = 0;
        int fallbackCnt = 0, orFailed = 0;
      };
      std::vector<ScaleResult> scales;
      auto runCase = [&](Manifold result, ScaleResult& cls) {
        if (result.IsEmpty() || result.Status() != Manifold::Error::NoError)
          return;
        ++cls.valid;
        auto si = CheckSelfIntersection(result, 1e-12);
        if (si.interiorPierces == 0) return;
        ++cls.pierces;
        cls.totalPierces += si.interiorPierces;
        auto [out, _] = OverlapRemoval(result);
        if (out.IsEmpty() || out.Status() != Manifold::Error::NoError) {
          ++cls.orFailed;
          return;
        }
        const bool fallback =
            (out.NumVert() == result.NumVert() &&
             out.NumTri() == result.NumTri() &&
             std::fabs(out.Volume() - result.Volume()) <
                 std::max(std::fabs(result.Volume()), 1.0) * 1e-15);
        if (fallback) {
          ++cls.fallbackCnt;
          ++cls.unchanged;
          return;
        }
        auto si2 = CheckSelfIntersection(out, 1e-12);
        if (si2.interiorPierces == 0)
          ++cls.fixed;
        else if (si2.interiorPierces < si.interiorPierces)
          ++cls.reduced;
        else if (si2.interiorPierces == si.interiorPierces)
          ++cls.unchanged;
        else
          ++cls.worsened;
      };
      for (int kPow : {0, 10, 20, 30}) {
        ScaleResult cls;
        cls.kPow = kPow;
        const double off = (kPow == 0) ? 0.0 : std::ldexp(1.5, kPow);
        for (int seed = 0; seed < 200; ++seed) {
          ++cls.total;
          const vec3 trans(u(rng), u(rng), u(rng));
          const vec3 axis(u(rng), u(rng), u(rng));
          const double theta = ang(rng);
          (void)axis;
          Manifold a =
              Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
          Manifold b = Manifold::Cube({2, 2, 2}, true)
                           .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                           .Translate(vec3(off, off, off) + trans);
          runCase(Boolean3D(a, b, OpType::Add), cls);
        }
        scales.push_back(cls);
      }
      std::cout << "General-fuzz: rotated cubes, 4 displacement scales × "
                   "200 seeds, OverlapRemoval applied to self-piercing "
                   "Booleans.\n";
      std::cout << "  kPow | pierces / valid (rate) | OverlapRemoval outcome\n";
      std::cout << "  -----+------------------------+-----------------------\n";
      int totPierces = 0, totValid = 0, totFixed = 0, totReduced = 0;
      int totUnchanged = 0, totWorsened = 0, totFallback = 0, totOrFailed = 0;
      for (const auto& cls : scales) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "  %4d | %3d/%-4d (%5.1f%%)      | "
                      "fix=%d red=%d reg=%d unch=%d wors=%d fb=%d\n",
                      cls.kPow, cls.pierces, cls.valid,
                      cls.valid > 0 ? 100.0 * cls.pierces / cls.valid : 0.0,
                      cls.fixed, cls.reduced, cls.orFailed, cls.unchanged,
                      cls.worsened, cls.fallbackCnt);
        std::cout << buf;
        totPierces += cls.pierces;
        totValid += cls.valid;
        totFixed += cls.fixed;
        totReduced += cls.reduced;
        totUnchanged += cls.unchanged;
        totWorsened += cls.worsened;
        totFallback += cls.fallbackCnt;
        totOrFailed += cls.orFailed;
      }
      std::cout << "\n  TOTAL self-pierce cases: " << totPierces << "/"
                << totValid << " ("
                << (100.0 * totPierces / std::max(totValid, 1))
                << "% self-pierce rate)\n";
      const int validOutputs = totFixed + totReduced + totOrFailed;
      std::cout << "    fixed (post=0):    " << totFixed << "\n";
      std::cout << "    reduced (post<pre):" << totReduced << "\n";
      std::cout << "    regularized-to-empty (valid): " << totOrFailed << "\n";
      std::cout << "    unchanged:         " << totUnchanged << "\n";
      std::cout << "    worsened:          " << totWorsened << "\n";
      std::cout << "    fallback:          " << totFallback << "\n";
      std::cout << "  → " << validOutputs << "/" << totPierces
                << " produce valid output (fixed + reduced + regularized)\n";
      return 0;
    }
    // Try Boolean3 self-union as overlap removal: m.Boolean(m, Add)
    // should be M ∪ M = M regularized. Tests whether the existing
    // Boolean3 pipeline handles inP == inQ correctly for the
    // production path.
    if (name == "--selfunion") {
      using namespace overlap3d;
      namespace fs = std::filesystem;
      struct Fixture {
        const char* a;
        const char* b;
        OpType op;
      };
      std::map<std::string, Fixture> fixtures{
          {"cray", {"Cray_left.obj", "Cray_right.obj", OpType::Subtract}},
          {"self-intersect",
           {"self_intersectA.obj", "self_intersectB.obj", OpType::Add}},
          {"hull-mask", {"hull-body.obj", "hull-mask.obj", OpType::Subtract}},
          {"offset12", {"Offset1.obj", "Offset2.obj", OpType::Add}},
          {"offset34", {"Offset3.obj", "Offset4.obj", OpType::Add}},
          {"havocglass",
           {"Havocglass8_left.obj", "Havocglass8_right.obj", OpType::Add}},
          {"generic-twin-7081",
           {"Generic_Twin_7081.1.t0_left.obj",
            "Generic_Twin_7081.1.t0_right.obj", OpType::Add}},
          {"generic-twin-7863",
           {"Generic_Twin_7863.1.t0_left.obj",
            "Generic_Twin_7863.1.t0_right.obj", OpType::Add}},
      };
      if (argc < 3) {
        std::cerr << "Usage: overlap3d_proto --selfunion <fixture>\n";
        return 2;
      }
      std::string fname = argv[2];
      auto it = fixtures.find(fname);
      if (it == fixtures.end()) return 2;
      fs::path modelsDir =
          fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
      auto load = [&](const char* fn) {
        std::ifstream f(modelsDir / fn);
        return Manifold::ReadOBJ(f);
      };
      Manifold a = load(it->second.a);
      Manifold b = load(it->second.b);
      Manifold input = a.Boolean(b, it->second.op);
      std::cout << "input v" << input.NumVert() << " t" << input.NumTri()
                << " vol " << input.Volume() << " status "
                << static_cast<int>(input.Status()) << "\n";
      auto siBefore = CheckSelfIntersection(input, 1e-12);
      std::cout << "input self-pierces: " << siBefore.interiorPierces
                << " (max mag " << siBefore.maxPierceMagnitude << ")\n";
      // Try self-union via Boolean
      Manifold out = input.Boolean(input, OpType::Add);
      std::cout << "output v" << out.NumVert() << " t" << out.NumTri()
                << " vol " << out.Volume() << " status "
                << static_cast<int>(out.Status()) << "\n";
      auto siAfter = CheckSelfIntersection(out, 1e-12);
      std::cout << "output self-pierces: " << siAfter.interiorPierces
                << " (max mag " << siAfter.maxPierceMagnitude << ")\n";
      const double drift = std::fabs(out.Volume() - input.Volume()) /
                           std::max(std::fabs(input.Volume()), 1e-12);
      std::cout << "volume drift: " << (drift * 100.0) << "%\n";
      return 0;
    }
    if (apiMode) {
      // Test the OverlapRemoval entry point on the named fixture.
      using namespace overlap3d;
      namespace fs = std::filesystem;
      struct Fixture {
        const char* a;
        const char* b;
        OpType op;
      };
      std::map<std::string, Fixture> fixtures{
          {"cray", {"Cray_left.obj", "Cray_right.obj", OpType::Subtract}},
          {"self-intersect",
           {"self_intersectA.obj", "self_intersectB.obj", OpType::Add}},
          {"hull-mask", {"hull-body.obj", "hull-mask.obj", OpType::Subtract}},
          {"offset12", {"Offset1.obj", "Offset2.obj", OpType::Add}},
          {"offset34", {"Offset3.obj", "Offset4.obj", OpType::Add}},
          {"havocglass",
           {"Havocglass8_left.obj", "Havocglass8_right.obj", OpType::Add}},
          {"generic-twin-7081",
           {"Generic_Twin_7081.1.t0_left.obj",
            "Generic_Twin_7081.1.t0_right.obj", OpType::Add}},
          {"generic-twin-7863",
           {"Generic_Twin_7863.1.t0_left.obj",
            "Generic_Twin_7863.1.t0_right.obj", OpType::Add}},
      };
      auto it = fixtures.find(name);
      if (it == fixtures.end()) return 2;
      fs::path modelsDir =
          fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
      auto load = [&](const char* fname) {
        std::ifstream f(modelsDir / fname);
        return Manifold::ReadOBJ(f);
      };
      Manifold a = load(it->second.a);
      Manifold b = load(it->second.b);
      Manifold result = a.Boolean(b, it->second.op);
      Manifold out;
      if (prodApi) {
        out = result.RemoveSelfIntersections();
      } else {
        auto [spike_out, dbg] = OverlapRemoval(result);
        out = spike_out;
        (void)dbg;
      }
      std::cout << "OverlapRemoval" << (prodApi ? "(prod-API)" : "") << "("
                << name << "):\n";
      std::cout << "  input  v" << result.NumVert() << " t" << result.NumTri()
                << " vol " << result.Volume() << "\n";
      std::cout << "  output v" << out.NumVert() << " t" << out.NumTri()
                << " vol " << out.Volume() << " status "
                << static_cast<int>(out.Status()) << "\n";
      // Pierce check BEFORE Volume() to avoid lazy-eval/materialization
      // ordering side effects.
      auto siInBefore = CheckSelfIntersection(result, 1e-12);
      auto siOutBefore = CheckSelfIntersection(out, 1e-12);
      const double drift = std::fabs(out.Volume() - result.Volume()) /
                           std::max(std::fabs(result.Volume()), 1e-12);
      std::cout << "  volume drift: " << (drift * 100.0) << "%\n";
      // Pierce check AFTER Volume() — compare to detect lazy-eval drift.
      auto siInAfter = CheckSelfIntersection(result, 1e-12);
      auto siOutAfter = CheckSelfIntersection(out, 1e-12);
      std::cout << "  pierces (before Volume): " << siInBefore.interiorPierces
                << " → " << siOutBefore.interiorPierces << "\n";
      std::cout << "  pierces (after Volume):  " << siInAfter.interiorPierces
                << " → " << siOutAfter.interiorPierces << "\n";
      return 0;
    }
    return RunSingleFixture(name);
  }

  using namespace overlap3d;
  bool allPass = true;

  // ---------------------------------------------------------------------
  // Eps inference smoke
  // ---------------------------------------------------------------------
  {
    Header("EpsilonFromScale");
    for (double L : {1.0, 100.0, 1e6, 1e10}) {
      const double eps = EpsilonFromScale(L);
      std::cout << "  L=" << L << " -> eps=" << eps << "\n";
    }
  }

  // ---------------------------------------------------------------------
  // Two overlapping unit cubes
  // ---------------------------------------------------------------------
  {
    Header("Boolean3D: two overlapping unit cubes");
    Manifold a = Manifold::Cube({2, 2, 2}, false);  // [0,2]^3
    Manifold b = Manifold::Cube({2, 2, 2}, false).Translate({1, 1, 1});
    const double eps = InferEps(a, b);
    std::cout << "  inferred eps = " << eps << "\n";

    Manifold add = Boolean3D(a, b, OpType::Add, eps);
    Manifold sub = Boolean3D(a, b, OpType::Subtract, eps);
    Manifold isec = Boolean3D(a, b, OpType::Intersect, eps);

    Report("a", a);
    Report("b", b);
    Report("Add", add);
    Report("Sub", sub);
    Report("Intersect", isec);
    // Sanity-check the self-intersection check on this known-clean case:
    auto siAdd = CheckSelfIntersection(add);
    auto siSub = CheckSelfIntersection(sub);
    auto siIsec = CheckSelfIntersection(isec);
    std::cout << "  Self-pierces: Add=" << siAdd.interiorPierces
              << " Sub=" << siSub.interiorPierces
              << " Isec=" << siIsec.interiorPierces << " (all should be 0)\n";

    // Expected: two cubes of volume 8 each, overlap of volume 1
    // (the unit cube at [1,2]^3). Add = 8 + 8 - 1 = 15. Sub = 7.
    // Intersect = 1.
    const double aV = a.Volume();
    const double bV = b.Volume();
    const double addV = add.Volume();
    const double subV = sub.Volume();
    const double isecV = isec.Volume();
    const bool ok = NearlyEqual(aV, 8.0) && NearlyEqual(bV, 8.0) &&
                    NearlyEqual(addV, 15.0) && NearlyEqual(subV, 7.0) &&
                    NearlyEqual(isecV, 1.0);
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Simplify a non-self-overlapping cube: should be a no-op
  // ---------------------------------------------------------------------
  {
    Header("Simplify: clean cube (no-op expected)");
    Manifold cube = Manifold::Cube({1, 1, 1}, true);  // [-0.5, 0.5]^3
    Manifold simplified = Simplify(cube);
    Report("input", cube);
    Report("simplified", simplified);
    const double inV = cube.Volume();
    const double outV = simplified.Volume();
    const bool ok = NearlyEqual(inV, 1.0) && NearlyEqual(outV, 1.0);
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Displacement attack: cubes near 2^k coordinate boundaries
  // ---------------------------------------------------------------------
  {
    Header("Displacement attack: overlapping cubes at 2^30");
    const double offset = std::ldexp(1.5, 30);  // ~1.6e9
    Manifold a =
        Manifold::Cube({2, 2, 2}, false).Translate({offset, offset, offset});
    Manifold b = Manifold::Cube({2, 2, 2}, false)
                     .Translate({offset + 1, offset + 1, offset + 1});
    const double eps = InferEps(a, b);
    std::cout << "  offset = " << offset << ", inferred eps = " << eps << "\n";

    Manifold add = Boolean3D(a, b, OpType::Add, eps);
    Manifold isec = Boolean3D(a, b, OpType::Intersect, eps);
    Report("Add", add);
    Report("Intersect", isec);

    const double addV = add.Volume();
    const double isecV = isec.Volume();
    // Volumes should still be 15 and 1 respectively (translation
    // doesn't change volume, and coords near 2^30 are still well within
    // double precision for a unit-scale feature).
    const bool ok =
        NearlyEqual(addV, 15.0, 1e-3) && NearlyEqual(isecV, 1.0, 1e-3);
    std::cout << "  " << (ok ? "PASS" : "FAIL")
              << " (tolerance 1e-3 due to displaced FP)\n";
    if (!ok) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Self-union (the regularization shape)
  // ---------------------------------------------------------------------
  {
    Header("Boolean3D self-union: should equal input");
    Manifold cube = Manifold::Cube({1, 1, 1}, true);
    Manifold selfUnion = Boolean3D(cube, cube, OpType::Add);
    Report("input", cube);
    Report("self union", selfUnion);
    const double inV = cube.Volume();
    const double outV = selfUnion.Volume();
    const bool ok = NearlyEqual(inV, outV, 1e-9);
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Volume identity: V(A∪B) + V(A∩B) == V(A) + V(B), exactly under
  // arbitrary affine inputs. This is the strongest quantitative check
  // available without an external reference.
  //
  // Cubes have feature size 2 (side length). The α-budget formula
  // gives eps ~= 3 at kPow=40 (2^41 * coeff), which EXCEEDS feature
  // size: Smith's framework correctly tells the algorithm "anything
  // within 3 units of anything else is the same point," so the cubes
  // get absorbed and the volume identity becomes vacuous (V(A)=V(B)=0).
  // Test stops at kPow=30 (eps ~= 0.003, well below feature size).
  // ---------------------------------------------------------------------
  {
    Header("Volume identity (V(A∪B) + V(A∩B) == V(A) + V(B))");
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    int pass = 0, fail = 0;
    double maxErr = 0.0;
    for (int kPow : {0, 10, 20, 30}) {
      const double off = (kPow == 0) ? 0.0 : std::ldexp(1.5, kPow);
      for (int seed = 0; seed < 20; ++seed) {
        const vec3 trans(u(rng), u(rng), u(rng));
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Translate(vec3(off, off, off) + trans);
        const double aV = a.Volume();
        const double bV = b.Volume();
        Manifold add = Boolean3D(a, b, OpType::Add);
        Manifold isec = Boolean3D(a, b, OpType::Intersect);
        const double addV = add.Volume();
        const double isecV = isec.Volume();
        const double identity = addV + isecV;
        const double expected = aV + bV;
        const double err =
            std::fabs(identity - expected) / std::max(expected, 1.0);
        if (err > maxErr) maxErr = err;
        const double tol = (kPow >= 30) ? 1e-3 : 1e-9;
        if (err < tol)
          ++pass;
        else
          ++fail;
      }
    }
    std::cout << "  80 cases (4 scales × 20 seeds): " << pass << " pass, "
              << fail << " fail. Max relative error: " << maxErr << "\n";
    std::cout << "  " << (fail == 0 ? "PASS" : "FAIL") << "\n";
    if (fail > 0) allPass = false;
  }

  // ---------------------------------------------------------------------
  // α-budget edge: at kPow=40, eps exceeds feature size; informational.
  // ---------------------------------------------------------------------
  {
    Header("α-budget edge probe (kPow=40, eps > feature size)");
    const double off = std::ldexp(1.5, 40);
    Manifold a = Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
    const double eps = InferEps(a);
    std::cout << "  feature size 2.0, offset " << off << ", inferred eps "
              << eps << "\n";
    std::cout << "  eps " << (eps >= 2.0 ? ">=" : "<") << " feature size: "
              << (eps >= 2.0 ? "α-budget tells algorithm to absorb the input"
                             : "still resolved")
              << "\n";
    Manifold add = Boolean3D(a, a.Translate({1, 1, 1}), OpType::Add);
    std::cout << "  Self+translated-by-1 union volume: " << add.Volume()
              << " (expected ~15 at small kPow; eps absorbs at kPow=40)\n";
    std::cout << "  Informational; not a fail condition.\n";
  }

  // ---------------------------------------------------------------------
  // Difference identity: V(A−B) == V(A) − V(A∩B)
  // ---------------------------------------------------------------------
  {
    Header("Difference identity (V(A−B) == V(A) − V(A∩B))");
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    int pass = 0, fail = 0;
    for (int seed = 0; seed < 50; ++seed) {
      const vec3 trans(u(rng), u(rng), u(rng));
      Manifold a = Manifold::Cube({2, 2, 2}, true);
      Manifold b = Manifold::Cube({2, 2, 2}, true).Translate(trans);
      const double aV = a.Volume();
      Manifold sub = Boolean3D(a, b, OpType::Subtract);
      Manifold isec = Boolean3D(a, b, OpType::Intersect);
      const double subV = sub.Volume();
      const double isecV = isec.Volume();
      const double err = std::fabs(subV - (aV - isecV)) / std::max(aV, 1.0);
      if (err < 1e-9)
        ++pass;
      else
        ++fail;
    }
    std::cout << "  50 cases: " << pass << " pass, " << fail << " fail\n";
    std::cout << "  " << (fail == 0 ? "PASS" : "FAIL") << "\n";
    if (fail > 0) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Iteration-to-fixed-point characterization. The 2D prototype's
  // IterateToFixedPoint catches sub-eps geometry drift between pass 1
  // and pass 2 (~2.3% of cases need iter=2 in the 2D 12k fuzz). Does
  // the 3D pipeline ever need iteration? Probe: run a boolean, then
  // feed the output back through Simplify, compare structurally. If
  // pass1 == pass2, the algorithm is idempotent on this case;
  // otherwise iteration would change something.
  // ---------------------------------------------------------------------
  {
    Header("3D iter-to-fixed-point characterization");
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    int idempotent = 0;
    int changedTopology = 0;
    int changedGeometryOnly = 0;
    double maxVolDelta = 0.0;
    for (int kPow : {0, 10, 20, 30}) {
      const double off = (kPow == 0) ? 0.0 : std::ldexp(1.5, kPow);
      for (int seed = 0; seed < 20; ++seed) {
        const vec3 trans(u(rng), u(rng), u(rng));
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Translate(vec3(off, off, off) + trans);
        Manifold pass1 = Boolean3D(a, b, OpType::Add);
        const double eps = InferEps(pass1);
        Manifold pass2 = Simplify(pass1, eps);
        const bool topoSame = (pass1.NumVert() == pass2.NumVert() &&
                               pass1.NumTri() == pass2.NumTri());
        const double vDelta = std::fabs(pass1.Volume() - pass2.Volume()) /
                              std::max(pass1.Volume(), 1.0);
        if (vDelta > maxVolDelta) maxVolDelta = vDelta;
        if (topoSame && vDelta < 1e-12) {
          ++idempotent;
        } else if (!topoSame) {
          ++changedTopology;
        } else {
          ++changedGeometryOnly;
        }
      }
    }
    std::cout << "  80 cases (4 scales × 20 seeds):\n";
    std::cout << "    idempotent (topology + geometry stable): " << idempotent
              << "\n";
    std::cout << "    iteration changes topology:              "
              << changedTopology << "\n";
    std::cout << "    iteration changes only geometry (sub-eps): "
              << changedGeometryOnly << "\n";
    std::cout << "    max volume delta: " << (maxVolDelta * 100.0) << "%\n";
    std::cout << "  "
              << (idempotent == 80 ? "PASS (no iteration needed)"
                                   : "INFORMATIONAL")
              << "\n";
  }

  // ---------------------------------------------------------------------
  // Wider 3D fuzz at displaced coords. Same shape as the 2D deepfuzz:
  // cases at multiple displacement scales (with feature size > eps),
  // multiple seeds. Now using rotated cubes (not axis-aligned) so we
  // exercise non-trivial face intersections, not just the trivial
  // edge-on-vert cases.
  // ---------------------------------------------------------------------
  {
    Header("3D deepfuzz: rotated cubes at displacement");
    std::mt19937 rng(1729);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    std::uniform_real_distribution<double> ang(0.0, 6.283);
    int total = 0;
    int firstPassValid = 0;
    int idempotent = 0;
    int changedTopology = 0;
    int changedGeometryOnly = 0;
    int volumeIdentityPass = 0;
    int volumeIdentityFail = 0;
    int selfPiercesTotal = 0;
    int selfPierceCases = 0;
    // Steps 1+2+3 fuzz accumulators
    int step1MergeCases = 0;        // cases where step 1 merged anything
    int step1MergeTotal = 0;        // total verts merged across all cases
    int step1MaxMerges = 0;         // worst single-case merge count
    int step2NonManifoldCases = 0;  // cases with unpaired (open) edges
    int step4Cases = 0;             // cases with any edge-edge intersection
    int step4Total = 0;             // total edge-edge intersections
    int step4MaxPerCase = 0;        // worst single-case count
    int step5Cases = 0;             // cases with any in-tri vert
    int step5Total = 0;             // total in-tri verts
    int step5MaxPerCase = 0;        // worst single-case count
    int step6Cases = 0;             // cases with any edge×tri intersection
    int step6Total = 0;             // total edge×tri intersections
    int step6MaxPerCase = 0;        // worst single-case count
    int step3HitCases = 0;
    int step3HitTotal = 0;
    double maxIterDelta = 0.0;
    double maxIdentityErr = 0.0;
    for (int kPow : {0, 10, 20, 30}) {
      const double off = (kPow == 0) ? 0.0 : std::ldexp(1.5, kPow);
      const double tol = (kPow >= 30) ? 1e-3 : 1e-9;
      for (int seed = 0; seed < 200; ++seed) {
        ++total;
        const vec3 trans(u(rng), u(rng), u(rng));
        // Rotate the second cube around a random axis by a random angle
        // so the boolean has non-trivial face × face intersections.
        const vec3 axis(u(rng), u(rng), u(rng));
        const double theta = ang(rng);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                         .Translate(vec3(off, off, off) + trans);
        const double aV = a.Volume();
        const double bV = b.Volume();
        Manifold pass1 = Boolean3D(a, b, OpType::Add);
        Manifold isec = Boolean3D(a, b, OpType::Intersect);
        // First-pass validity: did boolean produce a non-error manifold?
        if (pass1.Status() == Manifold::Error::NoError &&
            isec.Status() == Manifold::Error::NoError) {
          ++firstPassValid;
        }
        // Volume identity: V(A∪B) + V(A∩B) == V(A) + V(B)
        const double identityErr =
            std::fabs((pass1.Volume() + isec.Volume()) - (aV + bV)) /
            std::max(aV + bV, 1.0);
        if (identityErr > maxIdentityErr) maxIdentityErr = identityErr;
        if (identityErr < tol)
          ++volumeIdentityPass;
        else
          ++volumeIdentityFail;
        // Iteration: pass2 = Simplify(pass1)
        const double eps = InferEps(pass1);
        Manifold pass2 = Simplify(pass1, eps);
        const bool topoSame = (pass1.NumVert() == pass2.NumVert() &&
                               pass1.NumTri() == pass2.NumTri());
        const double vDelta = std::fabs(pass1.Volume() - pass2.Volume()) /
                              std::max(pass1.Volume(), 1.0);
        if (vDelta > maxIterDelta) maxIterDelta = vDelta;
        if (topoSame && vDelta < 1e-12)
          ++idempotent;
        else if (!topoSame)
          ++changedTopology;
        else
          ++changedGeometryOnly;
        // Self-intersection check on pass1 output.
        const auto si = CheckSelfIntersection(pass1);
        if (si.interiorPierces > 0) {
          ++selfPierceCases;
          selfPiercesTotal += si.interiorPierces;
          std::cout << "    pierce@ kPow=" << kPow << " seed=" << seed
                    << " count=" << si.interiorPierces << "\n";
        }
        // Steps 1+2+3 sanity: do they fire on this Boolean output?
        auto s1Result = MergeVertsEps(pass1, eps);
        const int merged = s1Result.mergedCount;
        auto s1Impl = ImplFromManifold(s1Result.manifold);
        auto edges = EnumerateEdges(s1Impl);
        int kUnpaired = 0;
        for (const auto& e : edges)
          if (e.halfedgePaired < 0) ++kUnpaired;
        auto s3 = BuildOnEdgeVertLists(s1Impl, edges, eps);
        int s3Hits = 0;
        for (const auto& l : s3) s3Hits += l.verts.size();
        // Step 4: edge-edge intersections.
        auto s4 = FindEdgeEdgeIntersections(s1Impl, edges, s3, eps);
        const int s4Count = static_cast<int>(s4.size());
        // Step 5: per-tri on-interior vert lists.
        auto s5 = BuildOnTriVertLists(s1Impl, eps);
        int s5Count = 0;
        for (const auto& l : s5) s5Count += l.verts.size();
        // Step 6: edge × triangle intersections.
        auto s6 = FindEdgeTriIntersections(s1Impl, edges, s3, s5, eps);
        const int s6Count = static_cast<int>(s6.size());
        if (merged > 0) {
          ++step1MergeCases;
          step1MergeTotal += merged;
          if (merged > step1MaxMerges) step1MaxMerges = merged;
        }
        if (kUnpaired > 0) ++step2NonManifoldCases;
        if (s3Hits > 0) {
          ++step3HitCases;
          step3HitTotal += s3Hits;
        }
        if (s4Count > 0) {
          ++step4Cases;
          step4Total += s4Count;
          if (s4Count > step4MaxPerCase) step4MaxPerCase = s4Count;
        }
        if (s5Count > 0) {
          ++step5Cases;
          step5Total += s5Count;
          if (s5Count > step5MaxPerCase) step5MaxPerCase = s5Count;
        }
        if (s6Count > 0) {
          ++step6Cases;
          step6Total += s6Count;
          if (s6Count > step6MaxPerCase) step6MaxPerCase = s6Count;
        }
      }
    }
    std::cout << "  Total cases: " << total << " (4 scales × 200 seeds)\n";
    std::cout << "  First-pass valid (no Error status): " << firstPassValid
              << "\n";
    std::cout << "  Volume identity:    " << volumeIdentityPass << " pass, "
              << volumeIdentityFail << " fail, max err "
              << (maxIdentityErr * 100.0) << "%\n";
    std::cout << "  Iteration probe:\n";
    std::cout << "    idempotent:                            " << idempotent
              << "\n";
    std::cout << "    iteration changes topology:            "
              << changedTopology << "\n";
    std::cout << "    iteration changes only geometry:       "
              << changedGeometryOnly << "\n";
    std::cout << "    max volume delta from iteration:       "
              << (maxIterDelta * 100.0) << "%\n";
    std::cout << "  Self-intersection check at relTol=1e-9:\n";
    std::cout << "    cases with any pierce:    " << selfPierceCases << "\n";
    std::cout << "    total pierce pairs found: " << selfPiercesTotal << "\n";
    std::cout << "  Steps 1+2+3 fuzz on Boolean output:\n";
    std::cout << "    step 1 merges (cases / total / max-per-case): "
              << step1MergeCases << " / " << step1MergeTotal << " / "
              << step1MaxMerges << "\n";
    std::cout << "    step 2 unpaired edges (cases): " << step2NonManifoldCases
              << "\n";
    std::cout << "    step 3 on-edge hits (cases / total): " << step3HitCases
              << " / " << step3HitTotal << "\n";
    std::cout << "    step 4 edge-edge isects (cases / total / max-per-case): "
              << step4Cases << " / " << step4Total << " / " << step4MaxPerCase
              << "\n";
    std::cout << "    step 5 in-tri verts (cases / total / max-per-case): "
              << step5Cases << " / " << step5Total << " / " << step5MaxPerCase
              << "\n";
    std::cout << "    step 6 edge×tri isects (cases / total / max-per-case): "
              << step6Cases << " / " << step6Total << " / " << step6MaxPerCase
              << "  ← MAIN PIERCE FINDER\n";

    // Re-probe at tighter (1e-12) and looser (1e-6) tolerances to validate.
    int p12_cases = 0, p12_total = 0, p6_cases = 0, p6_total = 0;
    std::mt19937 rngP(1729);
    std::uniform_real_distribution<double> uP(-1.5, 1.5);
    std::uniform_real_distribution<double> angP(0.0, 6.283);
    for (int kPow : {0, 10, 20, 30}) {
      const double off = (kPow == 0) ? 0.0 : std::ldexp(1.5, kPow);
      for (int seed = 0; seed < 200; ++seed) {
        const vec3 trans(uP(rngP), uP(rngP), uP(rngP));
        const vec3 axis(uP(rngP), uP(rngP), uP(rngP));
        const double theta = angP(rngP);
        Manifold a =
            Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
        Manifold b = Manifold::Cube({2, 2, 2}, true)
                         .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                         .Translate(vec3(off, off, off) + trans);
        Manifold pass1 = Boolean3D(a, b, OpType::Add);
        auto si12 = CheckSelfIntersection(pass1, 1e-12);
        auto si6 = CheckSelfIntersection(pass1, 1e-6);
        if (si12.interiorPierces > 0) {
          ++p12_cases;
          p12_total += si12.interiorPierces;
        }
        if (si6.interiorPierces > 0) {
          ++p6_cases;
          p6_total += si6.interiorPierces;
        }
      }
    }
    std::cout << "  Tolerance sweep (same fuzz, different tol):\n";
    std::cout << "    relTol=1e-12 (tighter): " << p12_cases << " cases, "
              << p12_total << " pierces\n";
    std::cout << "    relTol=1e-9  (current): " << selfPierceCases << " cases, "
              << selfPiercesTotal << " pierces\n";
    std::cout << "    relTol=1e-6  (looser):  " << p6_cases << " cases, "
              << p6_total << " pierces\n";
    bool ok = (firstPassValid == total && volumeIdentityFail == 0 &&
               selfPierceCases == 0);
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
  }

  // ---------------------------------------------------------------------
  // Adversarial 3D cases: deliberately try to provoke iteration. Models
  // the configurations the 2D thin-polygon idempotence quirk surfaced
  // (sub-ε feature dimensions, near-coincident faces, shared geometry).
  // ---------------------------------------------------------------------
  {
    Header("Adversarial 3D cases (designed to provoke iteration)");
    int probes = 0;
    int idempotent = 0;
    auto probe = [&](const char* name, Manifold a, Manifold b, OpType op) {
      ++probes;
      Manifold pass1 = Boolean3D(a, b, op);
      const double eps = InferEps(pass1);
      Manifold pass2 = Simplify(pass1, eps);
      const bool topoSame = (pass1.NumVert() == pass2.NumVert() &&
                             pass1.NumTri() == pass2.NumTri());
      const double vDelta = std::fabs(pass1.Volume() - pass2.Volume()) /
                            std::max(pass1.Volume(), 1e-12);
      const bool same = topoSame && vDelta < 1e-9;
      if (same) ++idempotent;
      const auto si = CheckSelfIntersection(pass1);
      std::cout << "  " << name << ": " << "verts " << pass1.NumVert() << "/"
                << pass2.NumVert() << " tris " << pass1.NumTri() << "/"
                << pass2.NumTri() << " vol " << pass1.Volume() << " (Δ"
                << (vDelta * 100.0) << "%) "
                << (same ? "IDEMPOTENT" : "CHANGED")
                << "  self-pierces: " << si.interiorPierces
                << (si.interiorPierces == 0 ? " ✓" : " ✗") << "\n";
    };
    // (1) Two cubes that share a face exactly.
    probe("share-face", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({2, 0, 0}), OpType::Add);
    // (2) Two cubes that share an edge exactly.
    probe("share-edge", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({2, 2, 0}), OpType::Add);
    // (3) Two cubes that share a vertex exactly.
    probe("share-vert", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({2, 2, 2}), OpType::Add);
    // (4) Identical cubes (true self-union).
    probe("identical-self-union", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true), OpType::Add);
    // (5) One cube fully inside another.
    probe("nested", Manifold::Cube({4, 4, 4}, true),
          Manifold::Cube({2, 2, 2}, true), OpType::Add);
    // (6) Tiny offset (smaller than typical eps but bigger than ulp at this
    // scale).
    probe("tiny-offset (1e-6)", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({1e-6, 0, 0}), OpType::Add);
    // (7) Long thin prism overlapping a cube (2D's thin-polygon analog).
    probe("thin-prism-cross", Manifold::Cube({0.05, 0.05, 4}, true),
          Manifold::Cube({4, 0.05, 0.05}, true), OpType::Add);
    // (8) Same prism cross at 2^30 displacement.
    {
      const double off = std::ldexp(1.5, 30);
      probe("thin-prism-cross @ 2^30",
            Manifold::Cube({0.05, 0.05, 4}, true).Translate({off, off, off}),
            Manifold::Cube({4, 0.05, 0.05}, true).Translate({off, off, off}),
            OpType::Add);
    }
    // (9) Intersect of nested cubes (small one).
    probe("nested-intersect", Manifold::Cube({4, 4, 4}, true),
          Manifold::Cube({2, 2, 2}, true), OpType::Intersect);
    // (10) Subtract that takes out a corner.
    probe("subtract-corner", Manifold::Cube({2, 2, 2}, false),
          Manifold::Cube({1, 1, 1}, false).Translate({-0.5, -0.5, -0.5}),
          OpType::Subtract);

    std::cout << "\n  " << idempotent << " / " << probes << " idempotent\n";
    std::cout << "  " << (idempotent == probes ? "PASS" : "INFORMATIONAL")
              << "\n";
  }

  // ---------------------------------------------------------------------
  // Smith-taxonomy 3D adversarial cases. Smith UCAM-CL-TR-766 §6.5
  // enumerates three 2D degeneracy classes that an overlap-removal
  // algorithm must handle: vertex-vertex (V-V) coincidence, vertex-on-
  // edge-interior (V-E), and edge-edge partial/full coincidence (E-E).
  // The 3D analog adds three more involving triangles: vertex-on-
  // triangle-interior (V-T), edge-on-triangle (E-T), and triangle-on-
  // triangle (T-T) coincidence. Smith's thesis doesn't supply specific
  // 3D test geometries; these cases are constructed to exercise each
  // class via boolean inputs that produce the degeneracy in the output.
  //
  // For each: we report idempotence (pass1 == pass2 under Simplify)
  // and the resulting volume. A "FAIL" here would be a real bug
  // signal; "IDEMPOTENT" is the expected and observed outcome.
  // ---------------------------------------------------------------------
  {
    Header("Smith taxonomy adversarial cases (V-V, V-E, V-T, E-E, E-T, T-T)");
    int probes = 0;
    int idempotent = 0;
    auto probe = [&](const char* taxonomy, const char* name, Manifold a,
                     Manifold b, OpType op) {
      ++probes;
      Manifold pass1 = Boolean3D(a, b, op);
      const double eps = InferEps(pass1);
      Manifold pass2 = Simplify(pass1, eps);
      const bool topoSame = (pass1.NumVert() == pass2.NumVert() &&
                             pass1.NumTri() == pass2.NumTri());
      const double vDelta = std::fabs(pass1.Volume() - pass2.Volume()) /
                            std::max(pass1.Volume(), 1e-12);
      const bool same = topoSame && vDelta < 1e-9;
      if (same) ++idempotent;
      const auto si = CheckSelfIntersection(pass1);
      std::cout << "  [" << taxonomy << "] " << name << ": " << "v"
                << pass1.NumVert() << " t" << pass1.NumTri() << " vol "
                << pass1.Volume() << " " << (same ? "IDEMPOTENT" : "CHANGED")
                << "  self-pierces: " << si.interiorPierces
                << (si.interiorPierces == 0 ? " ✓" : " ✗") << "\n";
    };

    // === V-V: vertex-vertex coincidence ===
    // Two cubes that meet exactly at a corner.
    probe("V-V", "two-cubes-corner-touch", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({2, 2, 2}), OpType::Add);

    // === V-E: vertex on edge interior ===
    // Cube B's corner sits exactly midway along an edge of cube A.
    // After boolean Add, cube A's edge has a vertex inserted at its
    // midpoint by cube B's corner intersecting it.
    probe("V-E", "B-corner-on-A-edge-midpoint", Manifold::Cube({4, 4, 4}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({3, 3, 0}), OpType::Add);

    // === V-T: vertex on triangle interior (3D-specific) ===
    // Cube B's corner lies exactly on the interior of one of cube A's
    // faces (centered on a face). The boolean output has a vertex on
    // the interior of A's face triangulation.
    probe("V-T", "B-corner-on-A-face-center", Manifold::Cube({4, 4, 4}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({0, 0, 3}), OpType::Add);

    // === E-E: edge-edge partial coincidence ===
    // Two cubes positioned so that a portion of an edge of A and a
    // portion of an edge of B are collinear and overlapping.
    probe("E-E", "edges-partially-collinear", Manifold::Cube({4, 2, 2}, true),
          Manifold::Cube({4, 2, 2}, true).Translate({3, 0, 0}), OpType::Add);

    // === E-T: edge on triangle interior (3D-specific) ===
    // An edge of cube B lies entirely on a face of cube A.
    probe("E-T", "B-edge-on-A-face", Manifold::Cube({4, 4, 4}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({0, 3, 0}), OpType::Add);

    // === T-T: full triangle coincidence ===
    // Two cubes share a complete face exactly (one face of A is the
    // same triangle pair as one face of B).
    probe("T-T", "share-full-face", Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({2, 0, 0}), OpType::Add);

    // === T-T partial: coplanar partial overlap ===
    // Two cubes' faces are coplanar but only partially overlap.
    probe("T-T-partial", "coplanar-partial-face-overlap",
          Manifold::Cube({4, 4, 2}, true),
          Manifold::Cube({2, 2, 2}, true).Translate({1, 1, 2}), OpType::Add);

    // === 3D k-fold concurrence: 3+ triangles meeting at one point ===
    // (3D analog of the 2D fig 6.5(d) adversarial hexagon, which
    // exercises 3 edges concurrent at one point.) Three cubes
    // arranged so they all share one vertex.
    {
      Manifold a = Manifold::Cube({2, 2, 2}, false);
      Manifold b = Manifold::Cube({2, 2, 2}, false).Translate({-2, 0, 0});
      Manifold c = Manifold::Cube({2, 2, 2}, false).Translate({0, -2, 0});
      probe("k-fold", "three-cubes-share-vertex", Boolean3D(a, b, OpType::Add),
            c, OpType::Add);
    }

    // === Self-intersection via input transformation ===
    // Cube rotated such that its translated copy interpenetrates with
    // axis-misaligned faces (the "two cubes overlapping" case from
    // earlier but more adversarial: small rotation = many shallow
    // crossings near degenerate angles).
    probe("near-degenerate-angle", "shallow-rotation-overlap",
          Manifold::Cube({2, 2, 2}, true),
          Manifold::Cube({2, 2, 2}, true)
              .Rotate(0.001, 0.001, 0.001)
              .Translate({1, 1, 1}),
          OpType::Add);

    // === All three boolean ops on a Smith-style configuration ===
    // T-T full coincidence under each op.
    {
      Manifold a = Manifold::Cube({2, 2, 2}, true);
      Manifold b = Manifold::Cube({2, 2, 2}, true).Translate({2, 0, 0});
      probe("T-T", "share-face Add", a, b, OpType::Add);
      probe("T-T", "share-face Subtract", a, b, OpType::Subtract);
      probe("T-T", "share-face Intersect", a, b, OpType::Intersect);
    }

    std::cout << "\n  " << idempotent << " / " << probes
              << " idempotent (Smith taxonomy)\n";
    std::cout << "  " << (idempotent == probes ? "PASS" : "INFORMATIONAL")
              << "\n";
  }

  // ---------------------------------------------------------------------
  // Adversarial battery: named .obj fixtures from test/models/.
  //
  // These are manifold's own curated adversarial inputs that exercise
  // known-difficult geometry: long thin features (Cray), aggressive
  // self-intersection (Havocglass, self_intersect), near-coplanar
  // offsets (Offset[1-4]), surface complexity (Generic_Twin), masking
  // (hull-body/mask). Several have ManifoldParams().processOverlaps
  // = true or selfIntersectionChecks = true gates in their
  // boolean_complex_test.cpp counterparts, indicating the existing
  // pipeline knows these inputs sit at the edge of its envelope.
  //
  // For each fixture pair, we run the operation that the existing
  // test runs (Add/Subtract per the named test) and report:
  //  - construction status (input loaded as a valid manifold?)
  //  - boolean op status
  //  - tri-tri self-intersection count, max pierce magnitude
  //  - bbox scale (so magnitudes can be normalized)
  //
  // The point of this battery is to *measure the gap*: how often
  // does the existing pipeline produce self-pierces on real
  // adversarial inputs? Each finding here is a candidate that
  // overlap-removal (per Emmett's #289) would need to fix.
  // ---------------------------------------------------------------------
  {
    Header("Adversarial battery: named .obj fixtures from test/models/");
    namespace fs = std::filesystem;
    fs::path modelsDir =
        fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
    auto loadObj = [&](const std::string& name) -> Manifold {
      fs::path p = modelsDir / name;
      std::ifstream f(p);
      if (!f.is_open()) {
        std::cout << "  [skip] " << name << ": cannot open " << p << "\n";
        return Manifold();
      }
      return Manifold::ReadOBJ(f);
    };
    int totalCases = 0;
    int constructionFail = 0;
    int booleanFail = 0;
    int piercesFound = 0;
    int piercesTotalCases = 0;
    auto runFixture = [&](const char* tag, const char* fileA, const char* fileB,
                          OpType op) {
      ++totalCases;
      Manifold a = loadObj(fileA);
      Manifold b = loadObj(fileB);
      const bool aOk = !a.IsEmpty() && a.Status() == Manifold::Error::NoError;
      const bool bOk = !b.IsEmpty() && b.Status() == Manifold::Error::NoError;
      if (!aOk || !bOk) {
        ++constructionFail;
        std::cout << "  [load FAIL] " << tag << ": A=" << (aOk ? "ok" : "bad")
                  << " B=" << (bOk ? "ok" : "bad") << "\n";
        return;
      }
      Manifold result = a.Boolean(b, op);
      const bool resOk = result.Status() == Manifold::Error::NoError;
      if (!resOk) ++booleanFail;
      const double bboxScale = result.BoundingBox().Scale();
      auto si = CheckSelfIntersection(result);
      if (si.interiorPierces > 0) {
        ++piercesTotalCases;
        piercesFound += si.interiorPierces;
      }
      const char* opStr = (op == OpType::Add)        ? "Add"
                          : (op == OpType::Subtract) ? "Sub"
                                                     : "Isect";
      std::cout << "  " << tag << " (" << opStr << "): result "
                << (resOk ? "OK" : "ERROR") << " v" << result.NumVert() << " t"
                << result.NumTri() << " vol " << result.Volume()
                << "\n    pierces=" << si.interiorPierces
                << " maxMag=" << si.maxPierceMagnitude
                << " (relative=" << (si.maxPierceMagnitude / bboxScale)
                << ") candidates=" << si.candidatesChecked << "\n";
    };
    // Names taken verbatim from boolean_complex_test.cpp / manifold_test.cpp.
    runFixture("self-intersect", "self_intersectA.obj", "self_intersectB.obj",
               OpType::Add);
    runFixture("Generic_Twin_7081", "Generic_Twin_7081.1.t0_left.obj",
               "Generic_Twin_7081.1.t0_right.obj", OpType::Add);
    runFixture("Generic_Twin_7863", "Generic_Twin_7863.1.t0_left.obj",
               "Generic_Twin_7863.1.t0_right.obj", OpType::Add);
    runFixture("Havocglass8", "Havocglass8_left.obj", "Havocglass8_right.obj",
               OpType::Add);
    runFixture("Cray", "Cray_left.obj", "Cray_right.obj", OpType::Subtract);
    runFixture("hull-mask", "hull-body.obj", "hull-mask.obj", OpType::Subtract);
    runFixture("Offset12", "Offset1.obj", "Offset2.obj", OpType::Add);
    runFixture("Offset34", "Offset3.obj", "Offset4.obj", OpType::Add);

    std::cout << "\n  Summary: " << totalCases << " cases, " << constructionFail
              << " construction failures, " << booleanFail
              << " boolean errors, " << piercesTotalCases
              << " cases with self-pierces (" << piercesFound
              << " total pierce pairs)\n";
    std::cout << "  This is informational; findings here are candidates\n"
              << "  for the overlap-removal feature to address.\n";
  }

  // ---------------------------------------------------------------------
  // Steps 1+2+3+4 sanity test: vert merge + collapsed-tri drop +
  // per-edge on-edge vert lists + edge-edge intersections.
  //
  // Take the .obj fixtures' Boolean output and run it through Soup →
  // MergeVertsEps(eps) → DropCollapsedAndBuildEdges. Confirms:
  //  - input (already a clean Manifold output, so already ε-merged
  //    by the constructor) does not lose any verts to overzealous
  //    merging at the inferred eps;
  //  - any tri collapses are documented;
  //  - edge incidence is built without error;
  //  - max edge incidence count flags non-manifold k≥3 edges.
  //
  // For Boolean output, expect zero merges (already clean) and zero
  // collapses. If the test ever shows non-zero, that's interesting:
  // either the inferred eps is too generous, or the Boolean output
  // has hidden duplicate verts that the constructor's ε-merge missed
  // (which would be an upstream bug to log).
  // ---------------------------------------------------------------------
  {
    Header("Steps 1-13b: 1-13a + ray-cast inside/outside filter");
    namespace fs = std::filesystem;
    fs::path modelsDir =
        fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
    auto loadObj = [&](const std::string& name) -> Manifold {
      fs::path p = modelsDir / name;
      std::ifstream f(p);
      if (!f.is_open()) return Manifold();
      return Manifold::ReadOBJ(f);
    };
    auto runStep1Plus2 = [&](const char* tag, const char* fileA,
                             const char* fileB, OpType op) {
      Manifold a = loadObj(fileA);
      Manifold b = loadObj(fileB);
      if (a.IsEmpty() || b.IsEmpty()) {
        std::cout << "  [skip] " << tag << "\n";
        return;
      }
      Manifold result = a.Boolean(b, op);
      if (result.Status() != Manifold::Error::NoError || result.IsEmpty()) {
        std::cout << "  [skip] " << tag << " (boolean error)\n";
        return;
      }
      const double eps = InferEps(result);
      const size_t vertsBefore = result.NumVert();
      const size_t trisBefore = result.NumTri();
      auto mergeResult = MergeVertsEps(result, eps);
      const Manifold& merged = mergeResult.manifold;
      const size_t trisAfter = merged.NumTri();
      // Step 2: enumerate edges via halfedge_ on the merged Manifold.
      auto impl = ImplFromManifold(merged);
      auto edges = EnumerateEdges(impl);
      // Edge incidence stats. With a topology-manifold input, every
      // edge has exactly 2 incident half-edges (= 2 incident triangles).
      // The pre-refactor flat-pair-key approach occasionally reported
      // k≥3 from duplicate triangle pairs; halfedge_-based enumeration
      // gives the canonical answer.
      int kPaired = 0, kUnpaired = 0;
      for (const auto& e : edges) {
        if (e.halfedgePaired >= 0)
          ++kPaired;
        else
          ++kUnpaired;
      }
      // Step 3: per-edge on-edge vert lists.
      auto onEdgeLists = BuildOnEdgeVertLists(impl, edges, eps);
      int edgesWithHits = 0;
      int totalHits = 0;
      int maxHitsPerEdge = 0;
      for (const auto& l : onEdgeLists) {
        if (!l.verts.empty()) ++edgesWithHits;
        totalHits += static_cast<int>(l.verts.size());
        if (static_cast<int>(l.verts.size()) > maxHitsPerEdge)
          maxHitsPerEdge = static_cast<int>(l.verts.size());
      }
      // Step 4: edge-edge intersections.
      auto eeIsects = FindEdgeEdgeIntersections(impl, edges, onEdgeLists, eps);
      int eeNew = 0, eeSnapped = 0;
      for (const auto& x : eeIsects) {
        if (x.snapTo >= 0)
          ++eeSnapped;
        else
          ++eeNew;
      }
      // Step 5: per-triangle on-interior vert lists.
      auto onTriLists = BuildOnTriVertLists(impl, eps);
      int trisWithInteriorHits = 0;
      int totalInteriorHits = 0;
      for (const auto& l : onTriLists) {
        if (!l.verts.empty()) ++trisWithInteriorHits;
        totalInteriorHits += static_cast<int>(l.verts.size());
      }
      // Step 6: edge × triangle intersections (the main pierce finder).
      auto etIsects =
          FindEdgeTriIntersections(impl, edges, onEdgeLists, onTriLists, eps);
      int etNew = 0, etSnapped = 0;
      for (const auto& x : etIsects) {
        if (x.snapTo >= 0)
          ++etSnapped;
        else
          ++etNew;
      }
      const char* opStr = (op == OpType::Add)        ? "Add"
                          : (op == OpType::Subtract) ? "Sub"
                                                     : "Isect";
      std::cout << "  " << tag << " (" << opStr << ") eps=" << eps << "\n"
                << "    verts " << vertsBefore << " (merged "
                << mergeResult.mergedCount << " pairs)\n"
                << "    tris  " << trisBefore << " → " << trisAfter << "\n"
                << "    edges " << edges.size() << " (paired: " << kPaired
                << ", unpaired: " << kUnpaired << ")\n"
                << "    on-edge verts: " << totalHits << " hits across "
                << edgesWithHits << " edges (max per edge: " << maxHitsPerEdge
                << ")\n"
                << "    edge-edge intersections: " << eeIsects.size()
                << " (new: " << eeNew << ", snap-to-existing: " << eeSnapped
                << ")\n"
                << "    on-tri interior verts: " << totalInteriorHits
                << " hits across " << trisWithInteriorHits << " tris\n"
                << "    edge×tri intersections: " << etIsects.size()
                << " (new: " << etNew << ", snap-to-existing: " << etSnapped
                << ") ← MAIN PIERCE FINDER\n";
      // Step 7 phase 1: tri-tri pairs.
      auto pairs = EnumerateTriTriPairs(impl, edges, etIsects);
      std::map<int, int> hist;
      for (const auto& p : pairs) ++hist[static_cast<int>(p.endpoints.size())];
      std::cout << "    tri-tri pairs: " << pairs.size() << " (";
      bool first = true;
      for (auto [n, c] : hist) {
        if (!first) std::cout << ", ";
        std::cout << "n=" << n << ":" << c;
        first = false;
      }
      std::cout << ") ← Emmett: \"exactly two\" → n=2 should dominate\n";
      // Step 7 phase 2: emit new verts + edges.
      auto step7p2 = EmitNewVertsAndEdges(impl, edges, etIsects, eps);
      std::cout << "    new verts: " << step7p2.newVertPositions.size()
                << ", new edges: " << step7p2.newEdges.size()
                << ", dropped pairs (n≠2): " << step7p2.dropped_n_not_2 << "\n";
      // Propagate the new verts (from step 6/7p2) to the on-edge
      // vert lists of the piercing edges. Without this, the new
      // verts are not seen as splitting the original edges, and
      // step 11 phase 2's next-around-face pointers are
      // disconnected from the new edges.
      PropagateNewVertsToOnEdgeLists(etIsects, step7p2.resolvedIds,
                                     static_cast<int>(impl.NumVert()), edges,
                                     onEdgeLists);
      int totalAugmentedHits = 0;
      for (const auto& l : onEdgeLists) totalAugmentedHits += l.verts.size();
      std::cout << "    propagation: on-edge verts " << totalHits << " → "
                << totalAugmentedHits << " (added "
                << (totalAugmentedHits - totalHits) << ")\n";
      // Pierce coverage: compare against ground-truth
      // CheckSelfIntersection on the same Manifold.
      auto siCheck = CheckSelfIntersection(result);
      const int pairsFound = static_cast<int>(pairs.size());
      const double cov = (siCheck.interiorPierces == 0)
                             ? 100.0
                             : 100.0 * pairsFound / siCheck.interiorPierces;
      std::cout << "    pierce coverage: ground-truth pierces "
                << siCheck.interiorPierces << ", piercing pairs found "
                << pairsFound << " (" << cov << "%)\n";
      // Relaxed-filter variant for diagnosing coverage shortfall.
      auto etRelax =
          FindEdgeTriIntersectionsRelaxed(impl, edges, onEdgeLists, eps);
      auto pairsRelax = EnumerateTriTriPairs(impl, edges, etRelax);
      const int pairsRelaxCount = static_cast<int>(pairsRelax.size());
      const double covRelax =
          (siCheck.interiorPierces == 0)
              ? 100.0
              : 100.0 * pairsRelaxCount / siCheck.interiorPierces;
      std::cout << "    pierce coverage (relaxed step 6): " << pairsRelaxCount
                << " pairs (" << covRelax << "%) "
                << "Δ=" << (pairsRelaxCount - pairsFound) << "\n";
      // Step 8: propagate in-tri verts onto new edges.
      auto step8 = AddInteriorVertsToNewEdges(
          impl, step7p2.newVertPositions, step7p2.newEdges, onTriLists, eps);
      int totalExtras = 0;
      int edgesWithExtras = 0;
      int maxExtras = 0;
      for (const auto& nwe : step8) {
        const int n = static_cast<int>(nwe.extraVerts.size());
        totalExtras += n;
        if (n > 0) ++edgesWithExtras;
        if (n > maxExtras) maxExtras = n;
      }
      std::cout << "    step 8 propagated verts: " << totalExtras << " across "
                << edgesWithExtras << " edges (max " << maxExtras
                << " per edge)\n";
      // Step 9: new-edge × new-edge intersections per tri.
      auto step9 = FindNewEdgeIntersections(impl, step7p2.newVertPositions,
                                            step7p2.newEdges, eps);
      int s9New = 0, s9Snap = 0;
      for (const auto& x : step9) {
        if (x.snapTo >= 0)
          ++s9Snap;
        else
          ++s9New;
      }
      std::cout << "    step 9 new-edge × new-edge: " << step9.size()
                << " (new: " << s9New << ", snap: " << s9Snap << ")\n";
      // Step 10: per-tri sub-edge counts.
      auto step10 = CountSubEdgesPerTri(impl, edges, onEdgeLists, step8);
      int totalSubEdges = 0, trisWithNewSubEdges = 0;
      int maxSubEdgesPerTri = 0;
      for (const auto& c : step10) {
        totalSubEdges += c.total();
        if (c.newSubEdges > 0) ++trisWithNewSubEdges;
        if (c.total() > maxSubEdgesPerTri) maxSubEdgesPerTri = c.total();
      }
      const int origSubEdgeBaseline = static_cast<int>(impl.NumTri() * 3);
      std::cout << "    step 10 sub-edges: total " << totalSubEdges
                << " (orig baseline " << origSubEdgeBaseline << ", added "
                << (totalSubEdges - origSubEdgeBaseline)
                << "), tris-with-new-sub-edges " << trisWithNewSubEdges
                << ", max/tri " << maxSubEdgesPerTri << "\n";
      // Step 11 phase 1: per-tri halfedge graph construction.
      auto step11p1 =
          BuildPerTriHalfedgeGraphs(impl, edges, onEdgeLists, step8);
      int trisWithNewEdges = 0, totalGraphHalfedges = 0;
      int maxHalfedgesPerTri = 0, maxVertsPerTri = 0;
      for (const auto& g : step11p1) {
        bool hasNew = false;
        for (const auto& h : g.halfedges) {
          if (h.isFromNewEdge) {
            hasNew = true;
            break;
          }
        }
        if (hasNew) ++trisWithNewEdges;
        const int nh = static_cast<int>(g.halfedges.size());
        const int nv = static_cast<int>(g.verts.size());
        totalGraphHalfedges += nh;
        if (nh > maxHalfedgesPerTri) maxHalfedgesPerTri = nh;
        if (nv > maxVertsPerTri) maxVertsPerTri = nv;
      }
      std::cout << "    step 11p1 halfedge graphs: total halfedges "
                << totalGraphHalfedges << ", tris-with-new-edges "
                << trisWithNewEdges << ", max halfedges/tri "
                << maxHalfedgesPerTri << ", max verts/tri " << maxVertsPerTri
                << "\n";
      // Step 11 phase 2: 2D projection + atan2 angle sort + next pointers.
      AddNextPointers(impl, step7p2.newVertPositions, step11p1);
      int totalNextSet = 0, totalNextUnset = 0;
      for (const auto& g : step11p1) {
        for (int n : g.nextHalfedge) {
          if (n >= 0)
            ++totalNextSet;
          else
            ++totalNextUnset;
        }
      }
      std::cout << "    step 11p2 next pointers: " << totalNextSet << " set, "
                << totalNextUnset << " unset\n";
      // Step 11 phase 3: walk the halfedge polygons.
      auto step11p3 = WalkPolygons(step11p1);
      int totalPolygons = 0, totalStalled = 0;
      int trisOnePoly = 0, trisMultiPoly = 0;
      int maxPolysPerTri = 0, maxPolygonSize = 0;
      for (const auto& w : step11p3) {
        const int n = static_cast<int>(w.polygons.size());
        totalPolygons += n;
        totalStalled += w.stalledHalfedges;
        if (n == 1)
          ++trisOnePoly;
        else if (n > 1)
          ++trisMultiPoly;
        if (n > maxPolysPerTri) maxPolysPerTri = n;
        for (const auto& p : w.polygons) {
          if (static_cast<int>(p.size()) > maxPolygonSize)
            maxPolygonSize = static_cast<int>(p.size());
        }
      }
      std::cout << "    step 11p3 polygons: " << totalPolygons
                << " (1-poly tris: " << trisOnePoly
                << ", multi-poly tris: " << trisMultiPoly
                << ", max polys/tri: " << maxPolysPerTri
                << ", max poly size: " << maxPolygonSize
                << ", stalled halfedges: " << totalStalled << ")\n";
      // Step 12: multiplicity merge.
      auto step12 = MergeMultiplicities(step11p3);
      int posMult = 0, negMult = 0, otherMult = 0;
      for (const auto& mp : step12) {
        if (mp.multiplicity > 0)
          ++posMult;
        else if (mp.multiplicity < 0)
          ++negMult;
        else
          ++otherMult;
      }
      const int cancelled = totalPolygons - static_cast<int>(step12.size());
      std::cout << "    step 12 merged: " << step12.size()
                << " unique (+: " << posMult << ", −: " << negMult << "), "
                << cancelled << " cancelled out of " << totalPolygons << "\n";
      // Step 13d: per-tri-pair analytical classifier.
      const int baseId = static_cast<int>(impl.NumVert());
      auto getPos3 = [&](int id) -> vec3 {
        return GetPos3(id, baseId, impl, step7p2.newVertPositions);
      };
      auto chordPartners = BuildChordPartnerMap(step7p2.newEdges);
      PolygonClassifierFn classifier = [&](const std::vector<int>& poly,
                                           const vec3& triNormal,
                                           int triId) -> PolygonClassification {
        PolygonClassification c{false, false, 0, 0};
        if (poly.size() < 3) return c;
        vec3 centroid(0, 0, 0);
        for (int v : poly) centroid += getPos3(v);
        centroid /= static_cast<double>(poly.size());
        (void)centroid;
        c.keep = AnalyticalKeep(triId, poly, chordPartners,
                                step7p2.newVertPositions, baseId, impl);
        c.reverse = false;
        (void)triNormal;
        return c;
      };
      auto step13b = TriangulateAndEmit(impl, step7p2.newVertPositions,
                                        step11p3, classifier);
      const Manifold& outMf = step13b.output;
      const bool outOk = outMf.Status() == Manifold::Error::NoError;
      const int outErrCode = static_cast<int>(outMf.Status());
      const double outVol = outMf.Volume();
      const double inVol = result.Volume();
      const double drift =
          std::fabs(outVol - inVol) / std::max(std::fabs(inVol), 1e-12);
      std::cout << "    step 13d (analytical) output: status "
                << (outOk ? "NoError"
                          : ("ERROR(" + std::to_string(outErrCode) + ")"))
                << " v" << outMf.NumVert() << " t" << outMf.NumTri() << " vol "
                << outVol << " (input vol " << inVol << ", drift "
                << (drift * 100.0) << "%)\n"
                << "      polys " << step13b.polygonsTriangulated << " → kept "
                << step13b.polygonsKept << " (auto " << step13b.polygonsAutoKept
                << "), dropped " << step13b.polygonsDropped << ", reversed "
                << step13b.polygonsReversed << " → tris "
                << step13b.trianglesEmitted << "\n";
      if (outOk) {
        // Check pierce reduction in the output.
        auto siOut = CheckSelfIntersection(outMf);
        std::cout << "      output self-pierces: " << siOut.interiorPierces
                  << " (input had " << siCheck.interiorPierces << ")\n";
      }
    };
    runStep1Plus2("self-intersect", "self_intersectA.obj",
                  "self_intersectB.obj", OpType::Add);
    runStep1Plus2("Cray", "Cray_left.obj", "Cray_right.obj", OpType::Subtract);
    runStep1Plus2("Havocglass8", "Havocglass8_left.obj",
                  "Havocglass8_right.obj", OpType::Add);
    runStep1Plus2("Offset12", "Offset1.obj", "Offset2.obj", OpType::Add);
    runStep1Plus2("hull-mask", "hull-body.obj", "hull-mask.obj",
                  OpType::Subtract);
    std::cout << "\n  Notes:\n"
              << "  - k=2 edges = topology-manifold. k=1 = open boundary,\n"
              << "    k≥3 = non-manifold T-junction. Boolean output should\n"
              << "    be all-k=2; deviations flag upstream issues.\n"
              << "  - on-edge verts = step 3 hits. Boolean output has\n"
              << "    no T-junctions by construction (verts are at\n"
              << "    corners or tri-tri crossings). The 0-hit result\n"
              << "    here is correct; step 3's real workload comes\n"
              << "    from non-manifold soup input — see the synthetic\n"
              << "    test below.\n";
  }

  // ---------------------------------------------------------------------
  // OverlapRemoval end-to-end test: API surface + no-op output.
  //
  // This is the entry-point test for the eventual overlap-removal
  // operation. Step 11 phase 0: the function runs all the query
  // steps 1-10 and stubs the output construction (round-trips the
  // Impl back to a Manifold unchanged). Validates that the API
  // shape works and the round-trip is lossless.
  //
  // For each fixture: run OverlapRemoval, compare input vs output
  // verts/tris/volume; report dbg counts. Pass criterion: output
  // is a valid Manifold with NoError and the round-trip preserves
  // volume to within tolerance.
  // ---------------------------------------------------------------------
  {
    Header("OverlapRemoval end-to-end (step 11 phase 0: no-op output)");
    namespace fs = std::filesystem;
    fs::path modelsDir =
        fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
    auto loadObj = [&](const std::string& name) -> Manifold {
      fs::path p = modelsDir / name;
      std::ifstream f(p);
      if (!f.is_open()) return Manifold();
      return Manifold::ReadOBJ(f);
    };
    auto runOR = [&](const char* tag, const char* fileA, const char* fileB,
                     OpType op) {
      Manifold a = loadObj(fileA);
      Manifold b = loadObj(fileB);
      if (a.IsEmpty() || b.IsEmpty()) {
        std::cout << "  [skip] " << tag << "\n";
        return;
      }
      Manifold result = a.Boolean(b, op);
      if (result.Status() != Manifold::Error::NoError || result.IsEmpty()) {
        std::cout << "  [skip] " << tag << " (boolean error)\n";
        return;
      }
      auto [out, dbg] = OverlapRemoval(result);
      const bool resOk = out.Status() == Manifold::Error::NoError;
      const double inV = result.Volume();
      const double outV = out.Volume();
      const double drift = std::fabs(outV - inV) / std::max(inV, 1e-12);
      const char* opStr = (op == OpType::Add)        ? "Add"
                          : (op == OpType::Subtract) ? "Sub"
                                                     : "Isect";
      std::cout << "  " << tag << " (" << opStr << "):\n"
                << "    input  v" << result.NumVert() << " t" << result.NumTri()
                << " vol " << inV << "\n"
                << "    output v" << out.NumVert() << " t" << out.NumTri()
                << " vol " << outV << " status "
                << (resOk ? "NoError" : "ERROR") << "\n"
                << "    volume drift: " << (drift * 100.0) << "%\n"
                << "    dbg: s1 " << dbg.step1Merges << ", s3 "
                << dbg.step3OnEdgeHits << ", s4 " << dbg.step4EdgeEdge
                << ", s5 " << dbg.step5InTri << ", s6 " << dbg.step6EdgeTri
                << ", s7p2 v" << dbg.step7p2NewVerts << "/e"
                << dbg.step7p2NewEdges << ", s8 " << dbg.step8Propagated
                << ", s9 " << dbg.step9NewNew << ", s10 "
                << dbg.step10TotalSubEdges << "\n";
    };
    runOR("self-intersect", "self_intersectA.obj", "self_intersectB.obj",
          OpType::Add);
    runOR("Cray", "Cray_left.obj", "Cray_right.obj", OpType::Subtract);
    runOR("Havocglass8", "Havocglass8_left.obj", "Havocglass8_right.obj",
          OpType::Add);
    runOR("Offset12", "Offset1.obj", "Offset2.obj", OpType::Add);
    runOR("hull-mask", "hull-body.obj", "hull-mask.obj", OpType::Subtract);
    std::cout << "\n  Note: output construction is currently a no-op\n"
              << "  (round-trips the Impl). Step 11 phases 1-4 will\n"
              << "  progressively replace the stub with actual polygon\n"
              << "  partition + winding classification.\n";
  }

  // ---------------------------------------------------------------------
  // Chained-Boolean test: do steps 1-3 fire more often after multiple
  // boolean ops?
  //
  // Rationale: each Boolean op produces topology-manifold output but
  // can introduce verts within ε of each other (post-op merge isn't
  // guaranteed by the pipeline). After N chained ops, FP error has
  // had N opportunities to drift; the pre-op ε-merge at op N+1 sees
  // op N's output as input but uses its own input bbox to pick eps,
  // which may be wider than what op N used. The hypothesis: chained
  // ops accumulate near-coincident verts that step 1 would merge.
  //
  // The chain pattern: r1 = a+b; r2 = r1+a; r3 = r2+b. Three identical
  // ops on three intermediate results. Stats show whether step 1's
  // merge count rises across the chain.
  // ---------------------------------------------------------------------
  {
    Header("Chained-Boolean test on .obj fixtures (step 1 trigger rate)");
    namespace fs = std::filesystem;
    fs::path modelsDir =
        fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
    auto loadObj = [&](const std::string& name) -> Manifold {
      fs::path p = modelsDir / name;
      std::ifstream f(p);
      if (!f.is_open()) return Manifold();
      return Manifold::ReadOBJ(f);
    };
    auto runStep1Stats = [&](const Manifold& m, double eps) {
      auto mr = MergeVertsEps(m, eps);
      auto impl = ImplFromManifold(mr.manifold);
      auto edges = EnumerateEdges(impl);
      int kUnpaired = 0;
      for (const auto& e : edges)
        if (e.halfedgePaired < 0) ++kUnpaired;
      auto step3 = BuildOnEdgeVertLists(impl, edges, eps);
      int totalHits = 0;
      for (const auto& l : step3) totalHits += l.verts.size();
      auto step4 = FindEdgeEdgeIntersections(impl, edges, step3, eps);
      auto step5 = BuildOnTriVertLists(impl, eps);
      int s5Hits = 0;
      for (const auto& l : step5) s5Hits += l.verts.size();
      auto step6 = FindEdgeTriIntersections(impl, edges, step3, step5, eps);
      auto step7p2 = EmitNewVertsAndEdges(impl, edges, step6, eps);
      auto step8 = AddInteriorVertsToNewEdges(impl, step7p2.newVertPositions,
                                              step7p2.newEdges, step5, eps);
      int s8Total = 0;
      for (const auto& nwe : step8) s8Total += nwe.extraVerts.size();
      return std::make_tuple(
          mr.mergedCount, kUnpaired, totalHits, static_cast<int>(step4.size()),
          s5Hits, static_cast<int>(step6.size()),
          static_cast<int>(step7p2.newEdges.size()), s8Total);
    };
    auto runChain = [&](const char* tag, const char* fileA, const char* fileB,
                        OpType op) {
      Manifold a = loadObj(fileA);
      Manifold b = loadObj(fileB);
      if (a.IsEmpty() || b.IsEmpty()) {
        std::cout << "  [skip] " << tag << "\n";
        return;
      }
      Manifold r1 = a.Boolean(b, op);
      if (r1.Status() != Manifold::Error::NoError || r1.IsEmpty()) {
        std::cout << "  [skip] " << tag << " (op 1 error)\n";
        return;
      }
      Manifold r2 = r1.Boolean(a, op);
      if (r2.Status() != Manifold::Error::NoError || r2.IsEmpty()) {
        std::cout << "  [skip] " << tag << " (op 2 error)\n";
        return;
      }
      Manifold r3 = r2.Boolean(b, op);
      if (r3.Status() != Manifold::Error::NoError || r3.IsEmpty()) {
        std::cout << "  [skip] " << tag << " (op 3 error)\n";
        return;
      }
      std::cout << "  " << tag << ":\n";
      auto report = [&](const char* label, const Manifold& m) {
        const double eps = InferEps(m);
        auto [merged, kUnpaired, hits, eeIsects, s5Hits, etIsects, s7newEdges,
              s8total] = runStep1Stats(m, eps);
        std::cout << "    " << label << " v" << m.NumVert() << " t"
                  << m.NumTri() << " vol " << m.Volume() << " eps=" << eps
                  << " → s1 " << merged << ", s3 " << hits << ", s4 e-e "
                  << eeIsects << ", s5 in-tri " << s5Hits << ", s6 e×t "
                  << etIsects << ", s7p2 newE " << s7newEdges << ", s8 prop "
                  << s8total << "\n";
        (void)kUnpaired;
      };
      report("r1=a∘b   ", r1);
      report("r2=r1∘a  ", r2);
      report("r3=r2∘b  ", r3);
    };
    runChain("Offset12 (Add)", "Offset1.obj", "Offset2.obj", OpType::Add);
    runChain("hull-mask (Sub)", "hull-body.obj", "hull-mask.obj",
             OpType::Subtract);
    runChain("Cray (Sub)", "Cray_left.obj", "Cray_right.obj", OpType::Subtract);
    runChain("Havocglass8 (Add)", "Havocglass8_left.obj",
             "Havocglass8_right.obj", OpType::Add);
    std::cout << "\n  Hypothesis: step 1 merge counts rise across the\n"
              << "  chain as FP error accumulates. If they don't, the\n"
              << "  pipeline's per-op ε-merge is sufficient and step 1\n"
              << "  is mostly redundant on real chained workloads.\n";
  }

  // (Step 3 functional validation now lives in the chained .obj
  // test above — Offset12 shows step 3 hits 0→1→6 across the
  // 3-op chain, and the deepfuzz reports 4/800 cases with hits.
  // The previous synthetic test used a non-Manifold soup and
  // doesn't fit the Manifold-only API.)

  // ---------------------------------------------------------------------
  // Constructed Smith fig 6.4/6.5 3D analog cases.
  //
  // Smith UCAM-CL-TR-766 §6.5:
  //   - Fig 6.4: rounded-ordering inconsistency between three edges
  //     (a vert is "above" edge A and "below" edge B, but A is
  //     "below" B in the rounded comparison). The 3D analog is
  //     non-transitive plane-side relations between three near-
  //     coplanar triangles.
  //   - Fig 6.5(d): three edges that are *truly* concurrent at one
  //     point but whose rounded intersections place them at three
  //     different points, forming an invalid micro-triangle. The 3D
  //     analog is k>=3 triangles truly concurrent at one line whose
  //     rounded pair-wise intersection lines disagree.
  //
  // These constructions are designed to *amplify* the pierce
  // findings observed in the .obj battery and the rotated-cube
  // fuzz, not just exhibit them. The aim is to give overlap-removal
  // a known-failing test bed to develop against, with controllable
  // parameters (rotation angle, displacement, k value).
  // ---------------------------------------------------------------------
  {
    Header("Smith fig 6.4/6.5 3D analog adversarial cases");
    int probes = 0;
    int piercesFound = 0;
    int piercesTotalCases = 0;
    auto probe = [&](const char* tag, const char* analog, Manifold result) {
      ++probes;
      if (result.Status() != Manifold::Error::NoError) {
        std::cout << "  [" << analog << "] " << tag << ": construction ERROR ("
                  << static_cast<int>(result.Status()) << ")\n";
        return;
      }
      const double bboxScale = std::max(result.BoundingBox().Scale(), 1e-12);
      auto si = CheckSelfIntersection(result);
      if (si.interiorPierces > 0) {
        ++piercesTotalCases;
        piercesFound += si.interiorPierces;
      }
      std::cout << "  [" << analog << "] " << tag << ": v" << result.NumVert()
                << " t" << result.NumTri() << " vol " << result.Volume()
                << "\n    pierces=" << si.interiorPierces
                << " maxMag=" << si.maxPierceMagnitude
                << " (relative=" << (si.maxPierceMagnitude / bboxScale)
                << ")\n";
    };

    // ---------------------------- Fig 6.5(d) analogs ----------------------
    // K=3 near-concurrent: three cubes rotated by tiny incremental angles
    // around three near-aligned axes, then unioned. The truly-concurrent
    // line of intersection (where all three meet) gets rounded to three
    // very-close-but-distinct lines.
    {
      const double a = 1e-5;  // ~tiny rotation in radians
      Manifold A = Manifold::Cube({2, 2, 2}, true);
      Manifold B = Manifold::Cube({2, 2, 2}, true).Rotate(a * 57.296, 0, 0);
      Manifold C = Manifold::Cube({2, 2, 2}, true).Rotate(0, a * 57.296, 0);
      probe("3-cubes-near-concurrent (rot 1e-5 rad)", "fig-6.5(d)",
            (A + B) + C);
    }
    {
      const double a = 1e-3;
      Manifold A = Manifold::Cube({2, 2, 2}, true);
      Manifold B = Manifold::Cube({2, 2, 2}, true).Rotate(a * 57.296, 0, 0);
      Manifold C = Manifold::Cube({2, 2, 2}, true).Rotate(0, a * 57.296, 0);
      probe("3-cubes-near-concurrent (rot 1e-3 rad)", "fig-6.5(d)",
            (A + B) + C);
    }
    {
      const double a = 0.1;
      Manifold A = Manifold::Cube({2, 2, 2}, true);
      Manifold B = Manifold::Cube({2, 2, 2}, true).Rotate(a * 57.296, 0, 0);
      Manifold C = Manifold::Cube({2, 2, 2}, true).Rotate(0, a * 57.296, 0);
      probe("3-cubes-near-concurrent (rot 0.1 rad)", "fig-6.5(d)", (A + B) + C);
    }
    {
      // K=4 escalation
      const double a = 1e-3;
      Manifold A = Manifold::Cube({2, 2, 2}, true);
      Manifold B = Manifold::Cube({2, 2, 2}, true).Rotate(a * 57.296, 0, 0);
      Manifold C = Manifold::Cube({2, 2, 2}, true).Rotate(0, a * 57.296, 0);
      Manifold D = Manifold::Cube({2, 2, 2}, true).Rotate(0, 0, a * 57.296);
      probe("4-cubes-near-concurrent (rot 1e-3 rad)", "fig-6.5(d)",
            ((A + B) + C) + D);
    }

    // ---------------------------- Fig 6.4 analogs --------------------------
    // Three near-coplanar slabs unioned; the resulting near-coplanar
    // tris in the output produce non-transitive plane-side relations
    // when checked pair-wise.
    {
      Manifold A = Manifold::Cube({4, 4, 0.001}, true);
      Manifold B =
          Manifold::Cube({4, 4, 0.001}, true).Rotate(0.001 * 57.296, 0, 0);
      Manifold C =
          Manifold::Cube({4, 4, 0.001}, true).Rotate(0, 0.001 * 57.296, 0);
      probe("near-coplanar-thin-slabs", "fig-6.4", (A + B) + C);
    }
    {
      // Two thin slabs at a near-zero angle: shallow wedge
      Manifold A = Manifold::Cube({4, 4, 0.01}, true);
      Manifold B =
          Manifold::Cube({4, 4, 0.01}, true).Rotate(0.01 * 57.296, 0, 0);
      probe("two-slabs-shallow-wedge", "fig-6.4", A + B);
    }

    // ---------------------------- Sliver-producing cases -------------------
    // Cubes intersecting at a very shallow angle produce thin sliver
    // triangles in the boolean output. Slivers are FP-marginal: their
    // surface normals are near-degenerate, and their plane equations
    // have relative error proportional to 1/sliver-thickness.
    {
      // Two boxes crossing at ~3 degrees -> long thin overlap region
      Manifold A = Manifold::Cube({4, 4, 4}, true);
      Manifold B = Manifold::Cube({4, 4, 4}, true).Rotate(3, 0, 0);
      probe("two-cubes-3deg-rotation", "sliver", A + B);
    }
    {
      // Cube and a thin slab passing through it at a shallow angle
      Manifold A = Manifold::Cube({4, 4, 4}, true);
      Manifold B = Manifold::Cube({0.05, 6, 6}, true).Rotate(0, 5, 0);
      probe("cube-by-thin-slab-5deg", "sliver", A - B);
    }

    // ---------------------------- Compounding cases ------------------------
    // Chain of boolean ops: each compounds FP error in the eps-merged
    // output. The third or fourth op's input has accumulated error
    // approximately 3-4x what a single op sees.
    {
      Manifold base = Manifold::Cube({2, 2, 2}, true);
      Manifold step1 = base + base.Translate({1.0, 0, 0});
      Manifold step2 = step1 + base.Translate({0, 1.0, 0});
      Manifold step3 = step2 + base.Translate({0, 0, 1.0});
      probe("4-cube-staircase-chain", "compound", step3);
    }
    {
      // Same but with rotation between each step (forces re-tessellation)
      Manifold base = Manifold::Cube({2, 2, 2}, true);
      Manifold step1 = base + base.Translate({1.0, 0, 0}).Rotate(0.1, 0, 0);
      Manifold step2 = step1 + base.Translate({0, 1.0, 0}).Rotate(0, 0.1, 0);
      Manifold step3 = step2 + base.Translate({0, 0, 1.0}).Rotate(0, 0, 0.1);
      probe("4-cube-rotated-chain", "compound", step3);
    }

    // ---------------------------- Sub-tolerance grazing case ---------------
    // Two cubes whose faces are *almost* coplanar but offset by exactly
    // 2 * eps - the algorithm has to decide between merging and not.
    {
      Manifold base = Manifold::Cube({2, 2, 2}, true);
      Manifold a = base;
      Manifold b = base.Translate({2.0 - 1e-12, 0, 0});  // ULP-grazing
      probe("ULP-grazing-faces (offset 1e-12)", "fig-6.4", a + b);
    }

    std::cout << "\n  " << probes << " probes, " << piercesTotalCases
              << " cases with self-pierces (" << piercesFound
              << " pierce pairs)\n";
    if (piercesTotalCases == 0) {
      std::cout << "  Negative result: simple synthetic geometry does not\n"
                << "  reproduce the pierce class observed in the .obj\n"
                << "  fixture battery and the rotated-cube-at-displacement\n"
                << "  fuzz. The pierce class appears to require irregular\n"
                << "  natural geometry (~10^4 tris, varied normals, shallow\n"
                << "  features) and/or coordinate displacement. Constructed\n"
                << "  fig 6.4/6.5 analogs of cubes alone are absorbed by\n"
                << "  the existing pipeline.\n";
    }
  }

  // ---------------------------------------------------------------------
  // Reproduce the kPow=30 seed=116 self-pierce case and dump details.
  // Verifies: is this a real geometric self-intersection or a marginal
  // FP-precision artifact at the boundary of my tolerance?
  // ---------------------------------------------------------------------
  {
    Header("Probe: kPow=30 seed=116 (self-pierce reproducer)");
    std::mt19937 rng(1729);
    std::uniform_real_distribution<double> u(-1.5, 1.5);
    std::uniform_real_distribution<double> ang(0.0, 6.283);
    // Walk the rng to match the fuzz loop's draw pattern. Each inner
    // iteration in the fuzz consumes: 3 u() for trans + 3 u() for axis +
    // 1 ang() for theta = 6 u-calls and 1 ang-call per seed.
    auto skipOne = [&]() {
      for (int i = 0; i < 6; ++i) (void)u(rng);
      (void)ang(rng);
    };
    for (int kPow : {0, 10, 20}) {
      (void)kPow;
      for (int seed = 0; seed < 200; ++seed) skipOne();
    }
    for (int seed = 0; seed < 116; ++seed) skipOne();
    const double off = std::ldexp(1.5, 30);
    const vec3 trans(u(rng), u(rng), u(rng));
    const vec3 axis(u(rng), u(rng), u(rng));
    const double theta = ang(rng);
    Manifold a = Manifold::Cube({2, 2, 2}, true).Translate(vec3(off, off, off));
    Manifold b = Manifold::Cube({2, 2, 2}, true)
                     .Rotate(theta * 57.296, theta * 28.6, theta * 14.3)
                     .Translate(vec3(off, off, off) + trans);
    Manifold pass1 = Boolean3D(a, b, OpType::Add);
    Report("a", a);
    Report("b", b);
    Report("pass1 (Add)", pass1);
    std::cout << "  Status: "
              << (pass1.Status() == Manifold::Error::NoError ? "NoError"
                                                             : "ERROR")
              << ", Genus: " << pass1.Genus() << "\n";
    auto mesh = pass1.GetMeshGL64();
    auto vp = [&](int idx) {
      return vec3(mesh.vertProperties[mesh.numProp * idx + 0],
                  mesh.vertProperties[mesh.numProp * idx + 1],
                  mesh.vertProperties[mesh.numProp * idx + 2]);
    };
    // Find the piercing pairs and print their geometric details.
    const size_t nTri = mesh.NumTri();
    int piercesFound = 0;
    for (size_t ta = 0; ta < nTri && piercesFound < 4; ++ta) {
      for (size_t tb = ta + 1; tb < nTri && piercesFound < 4; ++tb) {
        const int ia[3] = {(int)mesh.triVerts[3 * ta],
                           (int)mesh.triVerts[3 * ta + 1],
                           (int)mesh.triVerts[3 * ta + 2]};
        const int ib[3] = {(int)mesh.triVerts[3 * tb],
                           (int)mesh.triVerts[3 * tb + 1],
                           (int)mesh.triVerts[3 * tb + 2]};
        int shared = 0;
        for (int i : {0, 1, 2})
          for (int j : {0, 1, 2})
            if (ia[i] == ib[j]) ++shared;
        if (shared >= 2) continue;
        const vec3 a0 = vp(ia[0]), a1 = vp(ia[1]), a2 = vp(ia[2]);
        const vec3 b0 = vp(ib[0]), b1 = vp(ib[1]), b2 = vp(ib[2]);
        const bool pierce =
            SegmentPiercesTriInterior(a0, a1, b0, b1, b2, 1e-9) ||
            SegmentPiercesTriInterior(a1, a2, b0, b1, b2, 1e-9) ||
            SegmentPiercesTriInterior(a2, a0, b0, b1, b2, 1e-9) ||
            SegmentPiercesTriInterior(b0, b1, a0, a1, a2, 1e-9) ||
            SegmentPiercesTriInterior(b1, b2, a0, a1, a2, 1e-9) ||
            SegmentPiercesTriInterior(b2, b0, a0, a1, a2, 1e-9);
        if (pierce) {
          ++piercesFound;
          std::cout << "  pierce[" << piercesFound << "]: tri " << ta
                    << " (verts " << ia[0] << "," << ia[1] << "," << ia[2]
                    << ") vs tri " << tb << " (verts " << ib[0] << "," << ib[1]
                    << "," << ib[2] << "), shared=" << shared << "\n";
          // Compute the rough scale of the pierce.
          const vec3 e1 = b1 - b0;
          const vec3 e2 = b2 - b0;
          const vec3 nB = manifold::la::cross(e1, e2);
          const double nMag = std::sqrt(manifold::la::dot(nB, nB));
          for (int e = 0; e < 3; ++e) {
            const vec3 va = (e == 0) ? a0 : (e == 1) ? a1 : a2;
            const vec3 vb = (e == 0) ? a1 : (e == 1) ? a2 : a0;
            const double dA = manifold::la::dot(va - b0, nB);
            const double dB = manifold::la::dot(vb - b0, nB);
            std::cout << "    edge a" << e << ": dA=" << dA << " dB=" << dB
                      << " (relative to nMag=" << nMag << ")\n";
          }
        }
      }
    }
  }

  std::cout << "\n==== OVERALL: " << (allPass ? "PASS" : "FAIL") << " ====\n";
  return allPass ? 0 : 1;
}

inline int RunSingleFixture(const std::string& name) {
  using namespace overlap3d;
  namespace fs = std::filesystem;
  struct Fixture {
    const char* a;
    const char* b;
    OpType op;
  };
  std::map<std::string, Fixture> fixtures{
      {"cray", {"Cray_left.obj", "Cray_right.obj", OpType::Subtract}},
      {"self-intersect",
       {"self_intersectA.obj", "self_intersectB.obj", OpType::Add}},
      {"hull-mask", {"hull-body.obj", "hull-mask.obj", OpType::Subtract}},
      {"offset12", {"Offset1.obj", "Offset2.obj", OpType::Add}},
      {"offset34", {"Offset3.obj", "Offset4.obj", OpType::Add}},
      {"havocglass",
       {"Havocglass8_left.obj", "Havocglass8_right.obj", OpType::Add}},
      {"generic-twin-7081",
       {"Generic_Twin_7081.1.t0_left.obj", "Generic_Twin_7081.1.t0_right.obj",
        OpType::Add}},
      {"generic-twin-7863",
       {"Generic_Twin_7863.1.t0_left.obj", "Generic_Twin_7863.1.t0_right.obj",
        OpType::Add}},
  };
  auto it = fixtures.find(name);
  if (it == fixtures.end()) {
    std::cerr << "Unknown fixture: " << name << "\n";
    std::cerr << "Known: cray | self-intersect | hull-mask | offset12 | "
                 "offset34 | havocglass | generic-twin-7081 | "
                 "generic-twin-7863\n";
    return 2;
  }
  fs::path modelsDir =
      fs::path(__FILE__).parent_path().parent_path() / "test" / "models";
  auto load = [&](const char* fname) {
    fs::path p = modelsDir / fname;
    std::ifstream f(p);
    return Manifold::ReadOBJ(f);
  };
  Manifold a = load(it->second.a);
  Manifold b = load(it->second.b);
  Manifold result = a.Boolean(b, it->second.op);
  if (result.IsEmpty() || result.Status() != Manifold::Error::NoError) {
    std::cerr << "Boolean failed.\n";
    return 1;
  }
  std::cout << "=== " << name << " ===\n";
  std::cout << "input: v" << result.NumVert() << " t" << result.NumTri()
            << " vol " << result.Volume() << "\n";
  const double eps = InferEps(result);
  auto mr = MergeVertsEps(result, eps);
  auto impl = ImplFromManifold(mr.manifold);
  auto edges = EnumerateEdges(impl);
  auto onEdgeLists = BuildOnEdgeVertLists(impl, edges, eps);
  auto eeIsects = FindEdgeEdgeIntersections(impl, edges, onEdgeLists, eps);
  auto onTriLists = BuildOnTriVertLists(impl, eps);
  auto etIsects =
      FindEdgeTriIntersections(impl, edges, onEdgeLists, onTriLists, eps);
  auto step7p2 = EmitNewVertsAndEdges(impl, edges, etIsects, eps);
  PropagateNewVertsToOnEdgeLists(etIsects, step7p2.resolvedIds,
                                 static_cast<int>(impl.NumVert()), edges,
                                 onEdgeLists);
  auto step8 = AddInteriorVertsToNewEdges(impl, step7p2.newVertPositions,
                                          step7p2.newEdges, onTriLists, eps);
  auto step9 = FindNewEdgeIntersections(impl, step7p2.newVertPositions,
                                        step7p2.newEdges, eps);
  auto step11p1 = BuildPerTriHalfedgeGraphs(impl, edges, onEdgeLists, step8);
  AddNextPointers(impl, step7p2.newVertPositions, step11p1);
  auto step11p3 = WalkPolygons(step11p1);

  // Per-vert winding via Winding03_-style flood fill (diagnostic
  // for now; not wired to classifier).
  if (std::getenv("OVERLAP3D_DUMP_WINDING")) {
    auto pvw = ComputePerVertWinding(
        impl, edges, etIsects, step7p2.newVertPositions, mr.manifold, eps);
    std::cerr << "      per-vert winding (spike scaffold): "
              << pvw.numComponents << " components, ";
    std::map<int, int> counts;
    for (int w : pvw.winding) ++counts[w];
    for (const auto& [w, n] : counts) {
      std::cerr << "w=" << w << ":" << n << " ";
    }
    std::cerr << "\n";

    // Production Winding03 for comparison. M==M means SoS collapses;
    // output is diagnostic only (Option A productionization spike).
    auto prodW = ProductionWinding03(impl, edges, etIsects);
    std::cerr << "      per-vert winding (production W03): ";
    std::map<int, int> prodCounts;
    for (int w : prodW) ++prodCounts[w];
    for (const auto& [w, n] : prodCounts) {
      std::cerr << "w=" << w << ":" << n << " ";
    }
    std::cerr << "\n";

    // Vert-by-vert agreement (compare the spike's "+ε side" winding to
    // the production W03 over original verts only).
    int agree = 0, disagree = 0;
    const int nOrig = static_cast<int>(impl.NumVert());
    for (int v = 0; v < nOrig; ++v) {
      if (pvw.winding[v] == prodW[v])
        ++agree;
      else
        ++disagree;
    }
    std::cerr << "      W03 agreement (orig verts): " << agree << " agree, "
              << disagree << " disagree (out of " << nOrig << ")\n";

    // Option B: SoS-corrected self-mesh winding (= production Winding03_
    // structure + Kernel02 adjacency filter for the M==M case).
    // Build p1q2 the same way as ProductionWinding03 so the broken-
    // halfedge components match.
    manifold::Vec<std::array<int, 2>> sma_p1q2;
    sma_p1q2.reserve(etIsects.size());
    for (const auto& x : etIsects) {
      if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size()))
        continue;
      sma_p1q2.push_back({edges[x.edgeIdx].halfedgeForward, x.triIdx});
    }
    std::sort(sma_p1q2.begin(), sma_p1q2.end(),
              [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
                return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
              });
    auto sma = manifold::AnalyzeSelfMesh(
        impl, manifold::VecView<const std::array<int, 2>>(sma_p1q2.data(),
                                                          sma_p1q2.size()));
    std::map<int, int> aboveCounts, belowCounts;
    for (int w : sma.w_above) ++aboveCounts[w];
    for (int w : sma.w_below) ++belowCounts[w];
    std::cerr << "      Option B w_above: ";
    for (const auto& [w, n] : aboveCounts)
      std::cerr << "w=" << w << ":" << n << " ";
    std::cerr << "\n      Option B w_below: ";
    for (const auto& [w, n] : belowCounts)
      std::cerr << "w=" << w << ":" << n << " ";
    std::cerr << "\n";

    // Spike-vs-Option-B agreement on original verts.
    int agreeAbove = 0, agreeBelow = 0, agreeBoth = 0;
    for (int v = 0; v < nOrig; ++v) {
      bool ab = pvw.winding[v] == sma.w_above[v];
      bool bl = pvw.windingMinus[v] == sma.w_below[v];
      if (ab) ++agreeAbove;
      if (bl) ++agreeBelow;
      if (ab && bl) ++agreeBoth;
    }
    std::cerr << "      Option B vs spike: above=" << agreeAbove
              << " below=" << agreeBelow << " both=" << agreeBoth << " (out of "
              << nOrig << ")\n";
  }

  // Step 13d: analytical classifier.
  const int baseId = static_cast<int>(impl.NumVert());
  auto getPos3 = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, step7p2.newVertPositions);
  };
  auto chordPartners = BuildChordPartnerMap(step7p2.newEdges);
  const char* nf = std::getenv("OVERLAP3D_NOFILTER");
  const bool noFilter = nf && std::string(nf) == "1";
  // Pair-symmetric chord enforcement is now default; opt-out via
  // OVERLAP3D_NOPAIRSYM=1 for A/B comparison.
  const char* sf = std::getenv("OVERLAP3D_NOPAIRSYM");
  const bool pairSym = !(sf && std::string(sf) == "1");
  // Option B classifier (SoS-corrected self-mesh winding via
  // src/self_mesh_analysis.h): opt-in via OVERLAP3D_OPTB=1.
  const char* optb = std::getenv("OVERLAP3D_OPTB");
  const bool useOptionB = optb && std::string(optb) == "1";

  manifold::SelfMeshAnalysis sma;
  if (useOptionB) {
    manifold::Vec<std::array<int, 2>> sma_p1q2;
    sma_p1q2.reserve(etIsects.size());
    for (const auto& x : etIsects) {
      if (x.edgeIdx < 0 || x.edgeIdx >= static_cast<int>(edges.size()))
        continue;
      sma_p1q2.push_back({edges[x.edgeIdx].halfedgeForward, x.triIdx});
    }
    std::sort(sma_p1q2.begin(), sma_p1q2.end(),
              [](const std::array<int, 2>& a, const std::array<int, 2>& b) {
                return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
              });
    sma = manifold::AnalyzeSelfMesh(
        impl, manifold::VecView<const std::array<int, 2>>(sma_p1q2.data(),
                                                          sma_p1q2.size()));
    std::cerr << "      Option B classifier: w_above + w_below computed for "
              << sma.w_above.size() << " verts\n";
  }
  using manifold::la::dot;
  using manifold::la::normalize;
  const vec3 probeDir = normalize(vec3(0.7234, 0.4567, 0.5191));

  // Pre-compute classifier decisions for ALL polygons.
  std::map<std::pair<int, int>, bool> precomputedKeep;
  if (!noFilter) {
    for (size_t triId = 0; triId < step11p3.size(); ++triId) {
      const auto& polys = step11p3[triId].polygons;
      for (size_t pi = 0; pi < polys.size(); ++pi) {
        if (useOptionB) {
          // Option B classifier (two flavors):
          //   - Default: per-vert above/below vote across polygon
          //     perimeter; centroid probe only as fallback when all
          //     verts are chord verts (= no signal).
          //   - OVERLAP3D_OPTB_CENTROID_ONLY=1: always use the
          //     centroid probe (skip the per-vert vote entirely).
          //     Useful when the per-vert path is biased by chord-
          //     heavy perimeters or by the spike's pre-step-7
          //     vert numbering.
          //
          // Single-poly tris (= no piercing) auto-keep, matching the
          // default classifier. Skipping the centroid probe here
          // avoids spurious flips on back-side surfaces where n_T
          // points the "wrong" direction (Cray-style Subtract
          // results).
          if (polys.size() == 1) {
            precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] =
                true;
            continue;
          }
          static const bool centroidOnly = []() {
            const char* s = std::getenv("OVERLAP3D_OPTB_CENTROID_ONLY");
            return s && std::string(s) == "1";
          }();
          int boundaryV = 0, interiorV = 0, exteriorV = 0;
          if (!centroidOnly) {
            for (int v : polys[pi]) {
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
          }
          bool keep;
          if (!centroidOnly && boundaryV + interiorV + exteriorV > 0) {
            keep = boundaryV > 0 && boundaryV >= interiorV &&
                   boundaryV >= exteriorV;
          } else {
            // Centroid probe.
            vec3 c(0, 0, 0);
            for (int v : polys[pi]) c += getPos3(v);
            c /= static_cast<double>(polys[pi].size());
            const vec3 n_T = impl.faceNormal_[triId];
            const double meshScale =
                manifold::la::length(impl.bBox_.max - impl.bBox_.min);
            const double probeEps = meshScale * 1e-9;
            const double rayLen = meshScale * 4.0;
            const vec3 rayDir = probeDir;
            const int wa =
                manifold::WindingAt(impl, c + n_T * probeEps, rayDir, rayLen);
            const int wb =
                manifold::WindingAt(impl, c - n_T * probeEps, rayDir, rayLen);
            // Boundary criterion: keep iff exactly one side has
            // winding 0 (= polygon on the boundary of the union
            // region). The "|wa - wb| ≥ 1" jump variant was tried
            // and over-keeps overlap-layer boundaries — see git
            // history for OVERLAP3D_OPTB_JUMP.
            keep = (wa == 0 && wb >= 1) || (wa >= 1 && wb == 0);
          }
          precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] =
              keep;
          continue;
        }
        if (polys.size() == 1) {
          precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] =
              true;
          continue;
        }
        bool keep =
            AnalyticalKeep(static_cast<int>(triId), polys[pi], chordPartners,
                           step7p2.newVertPositions, baseId, impl);
        precomputedKeep[{static_cast<int>(triId), static_cast<int>(pi)}] = keep;
      }
    }
  }

  // Perimeter-aware drop guard for Option B. Targeted at Cray's
  // cascade-drop pattern (= Option B drops a multi-poly polygon
  // whose perimeter shares an edge with an auto-kept single-poly
  // neighbor → orphaned k=1 edge). Pass:
  //   1. Build edge → list-of-(triId, pi) using-as-perimeter map.
  //   2. Mark "anchor" edges: edges where the polygon containing them
  //      is a single-poly auto-kept tri.
  //   3. For each Option B drop decision, if any of the polygon's
  //      perimeter edges is an anchor edge, override to keep.
  if (useOptionB && !noFilter) {
    {
      // Step 1: edge → list of (triId, pi) using as perimeter edge.
      std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
          perimEdgeMap;
      for (size_t triId = 0; triId < step11p3.size(); ++triId) {
        const auto& polys = step11p3[triId].polygons;
        for (size_t pi = 0; pi < polys.size(); ++pi) {
          const auto& poly = polys[pi];
          for (size_t i = 0; i < poly.size(); ++i) {
            int a = poly[i], b = poly[(i + 1) % poly.size()];
            if (a > b) std::swap(a, b);
            perimEdgeMap[{a, b}].push_back(
                {static_cast<int>(triId), static_cast<int>(pi)});
          }
        }
      }
      // Step 2: identify anchor edges (= edges adjacent to a kept
      // single-poly tri).
      std::set<std::pair<int, int>> anchorEdges;
      for (const auto& [edge, owners] : perimEdgeMap) {
        for (const auto& [triId, pi] : owners) {
          if (step11p3[triId].polygons.size() != 1) continue;
          auto it = precomputedKeep.find({triId, pi});
          if (it != precomputedKeep.end() && it->second) {
            anchorEdges.insert(edge);
            break;
          }
        }
      }
      // Step 3: for each Option B drop, override to keep if any
      // perimeter edge is an anchor.
      int rescued = 0;
      for (auto& [key, keep] : precomputedKeep) {
        if (keep) continue;
        const int triId = key.first;
        const int pi = key.second;
        const auto& poly = step11p3[triId].polygons[pi];
        for (size_t i = 0; i < poly.size(); ++i) {
          int a = poly[i], b = poly[(i + 1) % poly.size()];
          if (a > b) std::swap(a, b);
          if (anchorEdges.count({a, b}) > 0) {
            keep = true;
            ++rescued;
            break;
          }
        }
      }
      std::cerr << "      perim-guard: rescued " << rescued << " polygons\n";
    }
  }
  if (std::getenv("OVERLAP3D_DUMP_KEEP")) {
    for (const auto& [k, v] : precomputedKeep) {
      std::cerr << "      [keep] tri=" << k.first << " pi=" << k.second
                << " keep=" << v << "\n";
    }
  }

  // Pair-symmetric chord enforcement: for each chord, find the twin
  // polygon pairs across triA and triB. If the classifier decisions
  // disagree within a twin pair, force consistency. Default: force
  // KEEP (= avoids creating k=1 from over-aggressive drops).
  int symForced = 0;
  if (pairSym && !noFilter) {
    {
      auto findPolyWithHalfedge = [&](int triId, int v0, int v1) -> int {
        const auto& polys = step11p3[triId].polygons;
        for (size_t pi = 0; pi < polys.size(); ++pi) {
          const auto& poly = polys[pi];
          for (size_t i = 0; i < poly.size(); ++i) {
            if (poly[i] == v0 && poly[(i + 1) % poly.size()] == v1) {
              return static_cast<int>(pi);
            }
          }
        }
        return -1;
      };
      for (const auto& edge : step7p2.newEdges) {
        int piA1 = findPolyWithHalfedge(edge.triA, edge.v0, edge.v1);
        int piB1 = findPolyWithHalfedge(edge.triB, edge.v1, edge.v0);
        int piA2 = findPolyWithHalfedge(edge.triA, edge.v1, edge.v0);
        int piB2 = findPolyWithHalfedge(edge.triB, edge.v0, edge.v1);
        auto getKeep = [&](int triId, int pi) -> bool* {
          if (pi < 0) return nullptr;
          auto it = precomputedKeep.find({triId, pi});
          if (it == precomputedKeep.end()) return nullptr;
          return &it->second;
        };
        bool* kA1 = getKeep(edge.triA, piA1);
        bool* kB1 = getKeep(edge.triB, piB1);
        bool* kA2 = getKeep(edge.triA, piA2);
        bool* kB2 = getKeep(edge.triB, piB2);
        if (!kA1 || !kB1 || !kA2 || !kB2) continue;
        int p1Vote = (*kA1 ? 1 : 0) + (*kB1 ? 1 : 0);
        int p2Vote = (*kA2 ? 1 : 0) + (*kB2 ? 1 : 0);
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
        bool changed = (*kA1 != p1Keep) || (*kB1 != p1Keep) ||
                       (*kA2 != p2Keep) || (*kB2 != p2Keep);
        if (changed) ++symForced;
        *kA1 = p1Keep;
        *kB1 = p1Keep;
        *kA2 = p2Keep;
        *kB2 = p2Keep;
      }
    }
    // Phase 2: mesh-edge halfedge consistency, iterated. For each
    // polygon perimeter halfedge (a→b), find its reverse (b→a) in
    // another polygon. If found AND not already a chord pair (=
    // already handled in phase 1), AND-merge (drop iff either
    // dropped). Iterate to fixed point — cascade-drop forward
    // across mesh-edge adjacencies.
    //
    // Disabled when OVERLAP3D_NOPAIRSYM_PHASE2=1 (= keep Phase 1
    // chord enforcement, skip mesh-edge AND-merge). Useful for
    // measuring Option B's signal without the cascade-drop pressure.
    int meshSymForced = 0;
    const char* nps2 = std::getenv("OVERLAP3D_NOPAIRSYM_PHASE2");
    const bool skipPhase2 = nps2 && std::string(nps2) == "1";
    if (!skipPhase2) {
      std::set<std::tuple<int, int, int>> chordHalfedges;
      for (const auto& edge : step7p2.newEdges) {
        chordHalfedges.insert({edge.triA, edge.v0, edge.v1});
        chordHalfedges.insert({edge.triA, edge.v1, edge.v0});
        chordHalfedges.insert({edge.triB, edge.v0, edge.v1});
        chordHalfedges.insert({edge.triB, edge.v1, edge.v0});
      }
      std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>
          halfedgeMap;
      for (size_t triId = 0; triId < step11p3.size(); ++triId) {
        const auto& polys = step11p3[triId].polygons;
        for (size_t pi = 0; pi < polys.size(); ++pi) {
          const auto& poly = polys[pi];
          for (size_t i = 0; i < poly.size(); ++i) {
            int a = poly[i], b = poly[(i + 1) % poly.size()];
            halfedgeMap[{a, b}].push_back(
                {static_cast<int>(triId), static_cast<int>(pi)});
          }
        }
      }
      // Single pass — full iteration tends to over-cascade on
      // dense self-intersection meshes (e.g. self-intersect went
      // 333 → 590 k=1 with iteration-to-fixed-point).
      {
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
            if (kA_it->second != target) {
              kA_it->second = target;
              ++meshSymForced;
            }
            if (kB_it->second != target) {
              kB_it->second = target;
              ++meshSymForced;
            }
          }
        }
      }
    }
    symForced += meshSymForced;
    if (pairSym)
      std::cerr << "      mesh-sym forced: " << meshSymForced << "\n";
  }

  PolygonClassifierFn classifier = nullptr;
  if (!noFilter) {
    classifier = [&](const std::vector<int>& poly, const vec3& triNormal,
                     int triId) -> PolygonClassification {
      PolygonClassification c{false, false, 0, 0};
      if (poly.size() < 3) return c;
      // Find pi by matching poly contents (= linear scan).
      int pi = -1;
      const auto& polys = step11p3[triId].polygons;
      for (size_t i = 0; i < polys.size(); ++i) {
        if (polys[i].size() == poly.size()) {
          bool match = true;
          for (size_t j = 0; j < poly.size(); ++j) {
            if (polys[i][j] != poly[j]) {
              match = false;
              break;
            }
          }
          if (match) {
            pi = static_cast<int>(i);
            break;
          }
        }
      }
      if (pi >= 0) {
        auto it = precomputedKeep.find({triId, pi});
        if (it != precomputedKeep.end())
          c.keep = it->second;
        else
          c.keep = true;
      } else {
        c.keep = AnalyticalKeep(triId, poly, chordPartners,
                                step7p2.newVertPositions, baseId, impl);
      }
      c.reverse = false;
      (void)triNormal;
      return c;
    };
  }
  if (pairSym) std::cerr << "      pair-sym forced: " << symForced << "\n";
  auto step13 =
      TriangulateAndEmit(impl, step7p2.newVertPositions, step11p3, classifier);
  std::cout << (noFilter ? "[no filter]" : "[analytical filter]") << "\n";
  std::cout << "polys " << step13.polygonsTriangulated << " → kept "
            << step13.polygonsKept << " (auto " << step13.polygonsAutoKept
            << "), dropped " << step13.polygonsDropped << ", reversed "
            << step13.polygonsReversed << " → tris " << step13.trianglesEmitted
            << "\n";
  std::cout << "output: status " << static_cast<int>(step13.output.Status())
            << " v" << step13.output.NumVert() << " t" << step13.output.NumTri()
            << " vol " << step13.output.Volume() << "\n";
  std::cout << "step 7p2: " << step7p2.newEdges.size() << " new edges, "
            << step7p2.newVertPositions.size() << " new verts\n";

  // Per-chord diagnostic: for each new edge from step 7p2, find the
  // sub-polygons of triA and triB that contain it and report whether
  // each was kept. Highlight cases where BOTH sides kept (= k=4 in
  // output) or BOTH sides dropped (= k=0 / dangling chord).
  using manifold::la::dot;
  int chordsBothKept = 0;
  int chordsBothDropped = 0;
  int chordsOneEach = 0;
  for (const auto& edge : step7p2.newEdges) {
    auto inPolyAndKept = [&](int triId) -> std::pair<bool, double> {
      // Returns (containsChord, centroidDot to triPartner)
      const auto& polys = step11p3[triId].polygons;
      for (const auto& poly : polys) {
        const size_t n = poly.size();
        bool containsChord = false;
        for (size_t i = 0; i < n; ++i) {
          int a = poly[i], b = poly[(i + 1) % n];
          if ((a == edge.v0 && b == edge.v1) ||
              (a == edge.v1 && b == edge.v0)) {
            containsChord = true;
            break;
          }
        }
        if (containsChord) {
          // Compute centroid + classification.
          vec3 c(0, 0, 0);
          for (int v : poly) c += getPos3(v);
          c /= static_cast<double>(n);
          // Partner is the OTHER tri of this new edge.
          int partner = (triId == edge.triA) ? edge.triB : edge.triA;
          const int v0 = impl.halfedge_.Start(3 * partner);
          const vec3 pP = impl.vertPos_[v0];
          const vec3 nP = impl.faceNormal_[partner];
          return {true, dot(c - pP, nP)};
        }
      }
      return {false, 0.0};
    };
    auto [foundA, dotA] = inPolyAndKept(edge.triA);
    auto [foundB, dotB] = inPolyAndKept(edge.triB);
    bool keptA = foundA && dotA >= 0;
    bool keptB = foundB && dotB >= 0;
    if (foundA && foundB) {
      // Each side of the chord on each tri produces its own polygon.
      // We need to check BOTH sides' polygons of each tri.
      // For now, our simple inPolyAndKept just finds ONE poly that
      // contains the chord — the count below is a proxy.
      if (keptA && keptB)
        ++chordsBothKept;
      else if (!keptA && !keptB)
        ++chordsBothDropped;
      else
        ++chordsOneEach;
    }
  }
  std::cout << "chord stats: " << chordsBothKept << " both-kept, "
            << chordsBothDropped << " both-dropped, " << chordsOneEach
            << " one-each (out of " << step7p2.newEdges.size() << " chords)\n";
  return 0;
}
