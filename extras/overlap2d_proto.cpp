// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Prototype for 2D overlap removal (issue #289).
//
// Single-translation-unit prototype following Julian Smith,
// "Towards robust inexact geometric computation", UCAM-CL-TR-766
// (2009), Chapter 7, with the BVH adaptation sketched by elalish in
// github.com/elalish/manifold/issues/289.
//
// References (cited inline as "Smith §X" / "Smith fig X" below):
//   [Smith09] Julian M. Smith, "Towards robust inexact geometric
//             computation", PhD dissertation, University of Cambridge,
//             Technical Report UCAM-CL-TR-766, 2009.
//             https://www.cl.cam.ac.uk/techreports/UCAM-CL-TR-766.pdf
//
//             - Chapter 6: figs 6.4 and 6.5 (degenerate-vertex test
//               cases used by `Smith fig 6.5(d)` in the test battery).
//             - Chapter 7: overall algorithm.
//             - §7.4 onwards: robust Bentley-Ottmann.
//             - §7.7 / fig 7.16: iterate-to-fixed-point bound (≤2
//               iterations under symbolic positions).
//             - Table 7.3: PolySet2 (canonical sub-edges with signed
//               multiplicity; see `Canonicalize` and `OutEdge` here).
//             - Chapter 8: floating-point error analysis;
//               u = 2⁻⁵³ for IEEE 754 doubles, per-intersection
//               bound α = √153·u·L (used by `EpsilonFromScale`).
//
//   [Issue289] https://github.com/elalish/manifold/issues/289
//             elalish's 6-step BVH-adapted sketch (no sweep line),
//             which this prototype implements with Smith's symbolic
//             perturbation as the FP-robust core.
//
// Build (from the manifold repo root):
//   g++ -std=c++17 -O2 -I include -I src -DMANIFOLD_PAR=-1 \
//     -ffp-contract=off -fexcess-precision=standard \
//     extras/overlap2d_proto.cpp -o overlap2d_proto
// Run:
//   ./overlap2d_proto              # full test battery
//   ./overlap2d_proto diagnose 0   # diagnostic dump for one case
//   ./overlap2d_proto deepfuzz 100 # broader randomized verification
//
// Algorithm shape (matches manifold internals where possible):
//   - Spatial queries via manifold's `Collider` BVH (Morton-sorted leaves).
//   - Edge-edge intersection via a trim-and-`Interpolate` symbolic kernel
//     atop the `Shadows` orientation primitive from src/shared.h.
//   - Vertex equality via `manifold::DisjointSets` lock-free union-find.
//   - Step 6 face traversal via DCEL (the same structure manifold's 3D
//     `Manifold::Impl::halfedge_` uses), winding-rayed per face from a
//     point on the LEFT side of any boundary half-edge.
//   - Iterate to fixed point per Smith §7.7 (default maxIter=2, his bound).
//
// Single-threaded. Header-only manifold dependency (no link step, no TBB).
//
// Deferred for graduation to manifold's mainline build:
//   - Mechanical std::vector → manifold::Vec rename (Vec is the project's
//     CPU/GPU-portable vector). Drop-in API-compatible for our usage; the
//     conversion is ~135 declarations and is best done in the same patch
//     that wires the prototype into the build system, so the type churn
//     and the build-graph churn land together.
//   - ZoneScoped Tracy markers at phase boundaries.
//   - ExecutionContext::Impl* ctx threading for parallelism dispatch.
//   - Internal namespace (overlap2d::detail) hiding everything except
//     the public Simplify / Boolean2D entry points.

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include "../src/collider.h"
#include "../src/disjoint_sets.h"
#include "../src/iters.h"
#include "../src/parallel.h"
#include "../src/shared.h"
#include "../src/utils.h"
#include "manifold/common.h"

namespace overlap2d {

using manifold::Box;
using manifold::CCW;
using manifold::Collider;
using manifold::MakeSimpleRecorder;
using manifold::vec2;
using manifold::vec3;
using manifold::VecView;
using manifold::la::dot;
using manifold::la::length;

// Centered-shoelace signed area of a closed polygon loop. Same FP trick
// as the per-face area computation in step 6: subtract a reference vert
// before multiplying so products stay at edge-length scale instead of
// blowing up to O(L^2) at displaced coordinates. Total telescopes to the
// same answer as the raw shoelace because Sigma(b - a) around any closed
// loop is zero. Used by area-preservation regression tests and the
// area-drift tracking in DeepFuzz.
inline double SignedArea(const manifold::SimplePolygon& loop) {
  if (loop.size() < 3) return 0.0;
  const auto& r = loop[0];
  double sum = 0.0;
  for (size_t i = 0; i < loop.size(); ++i) {
    const auto& a = loop[i];
    const auto& b = loop[(i + 1) % loop.size()];
    const double ax = a.x - r.x, ay = a.y - r.y;
    const double bx = b.x - r.x, by = b.y - r.y;
    sum += ax * by - bx * ay;
  }
  return 0.5 * sum;
}

inline double TotalSignedArea(const manifold::Polygons& polys) {
  double total = 0.0;
  for (const auto& loop : polys) total += SignedArea(loop);
  return total;
}

// Per chapter 8: u = 2^-53 for double-precision IEEE 754.
constexpr double kU = 1.110223024625156540423631668e-16;
// Smith's per-intersection bound: alpha = sqrt(153) * u * L; sqrt(153)
// ~= 12.37.
constexpr double kAlphaCoeff = 12.37;

// Edge with signed multiplicity. v0/v1 index into the vertex vector.
// Used both as algorithm input and as the oriented sub-edge type in the
// step-6 output (alias `OutEdge` below).
struct EdgeM {
  int v0;
  int v1;
  int mult;  // +1 default; -1 for reversed contribution; etc.
};
using OutEdge = EdgeM;

// Signed area of an OverlapResult-style edge set, computed directly from
// the OutEdge list without round-tripping through OutEdgesToPolygons.
// Each directed edge contributes its centered shoelace term; the sum
// telescopes to the same answer as walking the boundary in order
// (because Sigma over a closed face's edges of (b - a) is zero).
// Avoids the round-trip drop-on-degenerate behavior (the 0.17%
// mismatch documented elsewhere) so it's the right primitive for the
// area-drift regression test in DeepFuzz.
inline double SignedAreaFromOutEdges(const std::vector<vec2>& verts,
                                     const std::vector<OutEdge>& edges) {
  if (edges.empty() || verts.empty()) return 0.0;
  const vec2 r = verts[edges[0].v0];
  double sum = 0.0;
  for (const auto& oe : edges) {
    const vec2& a = verts[oe.v0];
    const vec2& b = verts[oe.v1];
    const double ax = a.x - r.x, ay = a.y - r.y;
    const double bx = b.x - r.x, by = b.y - r.y;
    sum += oe.mult * (ax * by - bx * ay);
  }
  return 0.5 * sum;
}

// Choose epsilon for the operation. L = bounding box half-extent rounded
// up to power of 2. k_budget is the user's expected upper bound on how
// many times any one edge may be adjusted (default 1000 ~= 10^-12 L).
inline double EpsilonFromScale(double L, int k_budget = 1000) {
  // Round L up to power of 2 (Smith's analysis assumes this).
  if (L <= 0) return 0;
  int expBits;
  std::frexp(L, &expBits);
  const double L_pow2 = std::ldexp(1.0, expBits);
  return (k_budget + 1) * kAlphaCoeff * kU * L_pow2;
}

// Cross-detection via manifold's `CCW` (orientation predicate with
// tolerance). Two segments AB and CD strictly cross iff CCW(A, B, C) and
// CCW(A, B, D) have opposite signs AND CCW(C, D, A) and CCW(C, D, B)
// have opposite signs. CCW returns 0 within `tol` of collinear; one
// non-zero zero means a single endpoint lies on the other segment
// ("touch but don't cross"), which we reject. All four zero means the
// segments are mutually collinear; that case is resolved by the
// Edelsbrunner-Mücke SoS perturbation in `IntersectSegments` below;
// permutation parity of the four ranks decides whether the perturbed
// segments cross. Position computation in `IntersectSegments` is the
// trim-and-`Interpolate` kernel; FP noise in the resulting position is
// absorbed by the eps-radius merge in step 1 of the next iteration.
inline int CCWPerturbed(vec2 a, vec2 b, vec2 c, int rA, int rB, int rC,
                       double tol) {
  const int s = CCW(a, b, c, tol);
  if (s != 0) return s;
  // Sort (rA, rB, rC), counting transpositions. Even parity → +1.
  int sw = 0;
  if (rA > rB) { std::swap(rA, rB); ++sw; }
  if (rB > rC) { std::swap(rB, rC); ++sw; }
  if (rA > rB) { std::swap(rA, rB); ++sw; }
  return (sw & 1) ? -1 : 1;
}

// =============================================================================
// 2D edge-edge symbolic intersection (BVH-friendly).
//
// The classical Kernel11 from boolean3.cpp can't be used in BVH-pair-query
// context: it requires "one endpoint inside, one outside" the other
// segment's projection, which sweep-line guarantees but BVH pair queries
// don't (most pairs have one segment fully nested in the other's
// projection).
//
// This kernel works for any pair by **trimming both segments to their
// projection-axis overlap before applying Intersect**. Steps:
//   1. CCW + SoS for cross-or-not (existing predicate, untouched).
//   2. Pick the axis (x or y) where BOTH segments have non-zero spread,
//      preferring the larger min-spread for stability.
//   3. Sort each segment's endpoints L→R along that axis.
//   4. Compute axis-overlap interval [overlapL, overlapR] = intersection
//      of a's and b's axis spans.
//   5. Use `Interpolate` (from shared.h) to evaluate each segment's
//      orthogonal coord at overlapL and overlapR. This produces four
//      (axis, ortho) points all spanning the same axis interval, which
//      is the precondition Intersect's closed-form expects.
//   6. Apply boolean3 Intersect's closed-form (smaller |dy| endpoint
//      picked for FP stability) to compute the intersection position.
//
// Trimming makes both segments span the same axis interval by
// construction, so the Kernel11 "inside/outside endpoint" precondition is
// satisfied even for the nested-axis cases that arise from BVH pairs.
// =============================================================================
inline double Coord(vec2 p, int axis) { return axis == 0 ? p.x : p.y; }

inline bool IntersectSegments(vec2 a0, vec2 a1, vec2 b0, vec2 b1, int rA0,
                              int rA1, int rB0, int rB1, double eps,
                              vec2* out) {
  // Step 1: SoS-aware cross-detection (same as the non-symbolic path).
  const int s1 = CCW(a0, a1, b0, eps);
  const int s2 = CCW(a0, a1, b1, eps);
  const int s3 = CCW(b0, b1, a0, eps);
  const int s4 = CCW(b0, b1, a1, eps);
  const int zeros = (s1 == 0) + (s2 == 0) + (s3 == 0) + (s4 == 0);
  if (zeros > 0 && zeros < 4) return false;
  if (zeros == 4) {
    const int p1 = CCWPerturbed(a0, a1, b0, rA0, rA1, rB0, eps);
    const int p2 = CCWPerturbed(a0, a1, b1, rA0, rA1, rB1, eps);
    if (p1 == p2) return false;
    const int p3 = CCWPerturbed(b0, b1, a0, rB0, rB1, rA0, eps);
    const int p4 = CCWPerturbed(b0, b1, a1, rB0, rB1, rA1, eps);
    if (p3 == p4) return false;
  } else {
    if (s1 == s2 || s3 == s4) return false;
  }

  // Step 2: pick the axis where BOTH segments have non-zero spread, with
  // the larger spread of the two preferring stability (smaller |dy| works
  // better in Intersect when the trimmed segments are well-separated).
  // The min-spread per axis is what matters: a vertical segment has zero
  // x-spread, so x is unusable; we'd pick y. Without this check (just
  // bbox spread), the Smith hexagon's vertical CE segment causes
  // degenerate axis-overlap (overlapL == overlapR) and the kernel falsely
  // reports no intersection.
  const double aSpreadX = std::fabs(a1.x - a0.x);
  const double aSpreadY = std::fabs(a1.y - a0.y);
  const double bSpreadX = std::fabs(b1.x - b0.x);
  const double bSpreadY = std::fabs(b1.y - b0.y);
  const double xUsable = std::min(aSpreadX, bSpreadX);
  const double yUsable = std::min(aSpreadY, bSpreadY);

  // Special case: both segments are axis-aligned to opposite axes (one
  // horizontal, one vertical, exactly), so neither axis has both
  // segments contributing spread. Trim-and-Interpolate would degenerate
  // (zero-width overlap interval). Compute the intersection directly:
  // it's the cross of the vertical segment's constant x and the
  // horizontal segment's constant y. Common in real CAD/SVG inputs
  // (axis-aligned rectangles overlapping each other).
  if (xUsable == 0 && yUsable == 0) {
    const bool aHoriz = aSpreadX > 0 && aSpreadY == 0;
    const bool aVert = aSpreadY > 0 && aSpreadX == 0;
    const bool bHoriz = bSpreadX > 0 && bSpreadY == 0;
    const bool bVert = bSpreadY > 0 && bSpreadX == 0;
    if ((aHoriz && bVert) || (aVert && bHoriz)) {
      const double ix = aVert ? a0.x : b0.x;
      const double iy = aHoriz ? a0.y : b0.y;
      *out = vec2(ix, iy);
      return std::isfinite(out->x) && std::isfinite(out->y);
    }
    return false;  // both points, or other degenerate config
  }

  const int axis = xUsable >= yUsable ? 0 : 1;

  // Step 3: sort each segment along the chosen axis.
  vec2 aL = a0, aR = a1, bL = b0, bR = b1;
  if (Coord(aR, axis) < Coord(aL, axis)) std::swap(aL, aR);
  if (Coord(bR, axis) < Coord(bL, axis)) std::swap(bL, bR);

  // Step 4: axis-overlap interval. CCW already confirmed crossing, so the
  // overlap is non-empty (modulo FP noise; clamp on hairline).
  const double overlapL = std::max(Coord(aL, axis), Coord(bL, axis));
  const double overlapR = std::min(Coord(aR, axis), Coord(bR, axis));
  if (overlapR <= overlapL) return false;

  // Step 5: trim each segment to the overlap. Interpolate(aL, aR, x) wants
  // a vec3 with the projection axis as its x-component, so we permute when
  // axis == 1 (y becomes x in the call frame, x becomes the orthogonal
  // coord). The returned vec2's first component is the orthogonal coord
  // at the requested projection value.
  auto embed = [&](vec2 p) {
    return axis == 0 ? vec3(p.x, p.y, 0.0) : vec3(p.y, p.x, 0.0);
  };
  const vec3 aL3 = embed(aL), aR3 = embed(aR);
  const vec3 bL3 = embed(bL), bR3 = embed(bR);
  const double aOL = manifold::Interpolate(aL3, aR3, overlapL).x;
  const double aOR = manifold::Interpolate(aL3, aR3, overlapR).x;
  const double bOL = manifold::Interpolate(bL3, bR3, overlapL).x;
  const double bOR = manifold::Interpolate(bL3, bR3, overlapR).x;

  // Step 6: Intersect closed-form. dyL/dyR are the ortho gaps at the two
  // overlap boundaries; pick whichever has smaller |dy| as the lambda
  // basepoint for FP stability (port from boolean3.cpp:36-54).
  const double dyL = bOL - aOL;
  const double dyR = bOR - aOR;
  const bool useL = std::fabs(dyL) < std::fabs(dyR);
  const double dProj = overlapR - overlapL;
  double lambda = (useL ? dyL : dyR) / (dyL - dyR);
  if (!std::isfinite(lambda)) return false;
  const double outProj = lambda * dProj + (useL ? overlapL : overlapR);
  const double aDy = aOR - aOL;
  const double bDy = bOR - bOL;
  const bool useA = std::fabs(aDy) < std::fabs(bDy);
  const double outOrtho = lambda * (useA ? aDy : bDy) +
                          (useL ? (useA ? aOL : bOL) : (useA ? aOR : bOR));
  *out = axis == 0 ? vec2(outProj, outOrtho) : vec2(outOrtho, outProj);
  return std::isfinite(out->x) && std::isfinite(out->y);
}


// =============================================================================
// Tier 1: BVH wrappers.
//
// Mechanical adapters around manifold's `Collider` (src/collider.h). The
// collider operates on 3D AABBs; we embed 2D inputs at z=0. These helpers
// are the only place where manifold::Box / manifold::Collider show up; the
// algorithm steps below use them as a black-box "give me overlapping pairs"
// service.
// =============================================================================

// 2D point as eps-padded 3D AABB (z=0 plane).
inline Box BoxOf2DPoint(vec2 p, double eps) {
  vec3 p3(p.x, p.y, 0);
  vec3 pad(eps, eps, 0);
  return Box(p3 - pad, p3 + pad);
}

// 2D segment as eps-padded 3D AABB (z=0 plane).
inline Box BoxOf2DEdge(vec2 p0, vec2 p1, double eps) {
  Box b(vec3(p0.x, p0.y, 0), vec3(p1.x, p1.y, 0));
  vec3 pad(eps, eps, 0);
  return Box(b.min - pad, b.max + pad);
}

// BVH built from a flat list of leaf boxes. The Collider expects leaves in
// Morton-sorted order, so we sort internally and remember the permutation
// so callers can keep using their original indices throughout.
struct BVH {
  Collider collider;
  std::vector<int> leafToOrig;  // sortedLeafIdx -> caller's input index
};

inline BVH BVHBuildFromBoxes(const std::vector<Box>& boxes) {
  const int n = static_cast<int>(boxes.size());
  BVH out;
  out.leafToOrig.resize(n);
  for (int i = 0; i < n; ++i) out.leafToOrig[i] = i;
  if (n == 0) return out;
  Box bbox;
  for (const auto& b : boxes) bbox = bbox.Union(b);
  std::vector<uint32_t> morton(n);
  for (int i = 0; i < n; ++i)
    morton[i] = Collider::MortonCode(boxes[i].Center(), bbox);
  std::stable_sort(out.leafToOrig.begin(), out.leafToOrig.end(),
                   [&](int a, int b) { return morton[a] < morton[b]; });
  std::vector<Box> sortedBB(n);
  std::vector<uint32_t> sortedMorton(n);
  for (int i = 0; i < n; ++i) {
    sortedBB[i] = boxes[out.leafToOrig[i]];
    sortedMorton[i] = morton[out.leafToOrig[i]];
  }
  out.collider =
      Collider(VecView<const Box>(sortedBB.data(), sortedBB.size()),
               VecView<const uint32_t>(sortedMorton.data(), sortedMorton.size()));
  return out;
}

// Collide every query box in `queries` against the BVH and call
// `f(queryIdx, origLeafIdx)` once per overlap. `origLeafIdx` is in the
// caller's input space (we undo the Morton permutation here). Always
// sequential; the prototype runs with MANIFOLD_PAR=-1 anyway.
template <typename F>
inline void CollidePairs(const BVH& bvh, const std::vector<Box>& queries,
                         F&& f) {
  if (bvh.leafToOrig.empty() || queries.empty()) return;
  auto adapter = [&](int qi, int leafIdx) { f(qi, bvh.leafToOrig[leafIdx]); };
  auto recorder = MakeSimpleRecorder(adapter);
  auto qf = [&](int i) { return queries[i]; };
  bvh.collider.Collisions<false>(recorder, qf,
                                 static_cast<int>(queries.size()),
                                 /*parallel=*/false);
}

// =============================================================================
// Per-phase wall-clock accumulator. Atomic so concurrent RemoveOverlaps2D
// calls (e.g. from a multi-threaded corpus driver, future) can share.
// Currently RemoveOverlaps2D itself is single-threaded; the accumulator
// captures cumulative time across many cases for the `time` CLI mode.
// Defined here so step 3/4 helpers can reference it.
// =============================================================================
struct PhaseAcc {
  std::atomic<int64_t> mergeNs{0};
  std::atomic<int64_t> remapNs{0};
  std::atomic<int64_t> buildListsNs{0};
  std::atomic<int64_t> findIxNs{0};
  std::atomic<int64_t> restructNs{0};
  std::atomic<int64_t> canonNs{0};
  std::atomic<int64_t> filterDcelNs{0};
  std::atomic<int64_t> totalNs{0};
  std::atomic<int64_t> cases{0};
  // Fine-grain sub-phase timers: BVH build, broad phase (collide-pairs),
  // narrow phase (per-pair work), eager propagation. Steps 3 and 4 share
  // edgeBoxes + bvh, so build/build is counted once at the call site.
  std::atomic<int64_t> bvhBuildNs{0};
  std::atomic<int64_t> step3BroadNs{0};
  std::atomic<int64_t> step3NarrowNs{0};
  std::atomic<int64_t> step4BroadNs{0};
  std::atomic<int64_t> step4NarrowNs{0};
  std::atomic<int64_t> step4PropNs{0};
  std::atomic<int64_t> step1BvhBuildNs{0};
  std::atomic<int64_t> step1CollideNs{0};
  std::atomic<int64_t> step1RestNs{0};
  void Reset() {
    mergeNs = 0; remapNs = 0; buildListsNs = 0; findIxNs = 0;
    restructNs = 0; canonNs = 0; filterDcelNs = 0; totalNs = 0; cases = 0;
    bvhBuildNs = 0; step3BroadNs = 0; step3NarrowNs = 0;
    step4BroadNs = 0; step4NarrowNs = 0; step4PropNs = 0;
    step1BvhBuildNs = 0; step1CollideNs = 0; step1RestNs = 0;
  }
};
inline PhaseAcc& GlobalPhases() { static PhaseAcc p; return p; }

namespace timing_detail {
using Clock = std::chrono::steady_clock;
inline int64_t Ns(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
}
}  // namespace timing_detail

// =============================================================================
// Step 1: vertex merge.
// Returns: remap[oldIdx] = newIdx, and merged vert positions.
// =============================================================================
struct VertexMerge {
  std::vector<int> remap;
  std::vector<vec2> verts;
};

VertexMerge MergeVerts(const std::vector<vec2>& in, double eps) {
  using Clock = std::chrono::steady_clock;
  auto Ns = [](Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
  };
  auto t0 = Clock::now();
  const int n = static_cast<int>(in.size());
  DisjointSets uf(n);
  // BVH broad phase: each vert's eps-padded box queried against all eps-
  // padded vert boxes. A pair (i, j) overlapping in box space gets a precise
  // distance check; pairs within eps unite. We collect candidates first,
  // sort by (min(i,j), max(i,j)), then unite. Sorting matches the order of
  // the original O(n^2) loop so DisjointSets's union-by-rank tie-break
  // produces the same cluster roots and the centroid map iterates the same
  // way (which controls assigned new-vertex indices).
  const double eps2 = eps * eps;
  // Broad-phase: collect candidate (i, j) pairs whose padded boxes overlap.
  // For small n, brute-force O(n²) is faster than the BVH build + Morton
  // sort + tree traversal — the constant factors for those routines
  // outweigh n² at this scale. The threshold (~32) was tuned empirically
  // on JTS-shape workloads. The pairs are sorted lex-ascending in both
  // paths so the unite step below is deterministic.
  std::vector<std::pair<int, int>> pairs;
  if (n < 32) {
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        // AABB overlap test (the same the BVH would do): both axes must
        // overlap on [v[i] - eps, v[i] + eps] ∩ [v[j] - eps, v[j] + eps].
        const vec2 d = in[i] - in[j];
        if (std::fabs(d.x) <= 2 * eps && std::fabs(d.y) <= 2 * eps)
          pairs.emplace_back(i, j);
      }
    }
  } else {
    // Sort-by-x sweep replaces the BVH self-collide for vertex-vertex
    // proximity. The BVH build (Morton sort + tree construction +
    // leaf-permutation table) costs ~20µs even when there are no
    // candidate pairs, and the typical CrossSection input is
    // *already* canonical with no close-vert pairs to find. For 1D
    // proximity (we only need pairs within 2·eps in BOTH axes),
    // sorting by x and scanning forward until x-distance exceeds 2·eps
    // is O(n log n) sort + O(n) average sweep — typically 5-10x
    // faster than BVH at the n's seen in JTS. Pairs are sorted
    // lex-ascending at the end so the unite step below stays
    // deterministic.
    auto tBvhStart = Clock::now();
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b) { return in[a].x < in[b].x; });
    auto tBvhEnd = Clock::now();
    GlobalPhases().step1BvhBuildNs.fetch_add(Ns(tBvhStart, tBvhEnd),
                                              std::memory_order_relaxed);
    auto tCollideStart = Clock::now();
    const double thresh = 2 * eps;
    for (int i = 0; i < n; ++i) {
      const int ai = idx[i];
      const double ax = in[ai].x;
      const double ay = in[ai].y;
      for (int j = i + 1; j < n; ++j) {
        const int bi = idx[j];
        const double dx = in[bi].x - ax;
        if (dx > thresh) break;
        if (std::fabs(in[bi].y - ay) > thresh) continue;
        if (ai < bi) pairs.emplace_back(ai, bi);
        else pairs.emplace_back(bi, ai);
      }
    }
    std::sort(pairs.begin(), pairs.end());
    auto tCollideEnd = Clock::now();
    GlobalPhases().step1CollideNs.fetch_add(Ns(tCollideStart, tCollideEnd),
                                             std::memory_order_relaxed);
  }
  // Fast path: no candidates → no merges, identity remap.
  if (pairs.empty()) {
    auto tRest = Clock::now();
    GlobalPhases().step1RestNs.fetch_add(Ns(t0, tRest) - 0,
                                          std::memory_order_relaxed);
    std::vector<int> remap(n);
    std::iota(remap.begin(), remap.end(), 0);
    return {std::move(remap), in};
  }
  // Parallelize the geometric distance gate (read-only on `in`); unite
  // serially in sorted pair order so cluster roots are deterministic
  // regardless of thread scheduling. See the matching pattern in step 4b.
  std::vector<uint8_t> doUnite(pairs.size(), 0);
  manifold::for_each_n(
      manifold::autoPolicy(pairs.size()), manifold::countAt(0),
      pairs.size(), [&](size_t k) {
        const auto [i, j] = pairs[k];
        vec2 d = in[i] - in[j];
        if (dot(d, d) <= eps2) doUnite[k] = 1;
      });
  bool anyMerge = false;
  for (size_t k = 0; k < pairs.size(); ++k) {
    if (doUnite[k]) {
      uf.unite(pairs[k].first, pairs[k].second);
      anyMerge = true;
    }
  }
  if (!anyMerge) {
    std::vector<int> remap(n);
    std::iota(remap.begin(), remap.end(), 0);
    return {std::move(remap), in};
  }
  // Compute centroid per cluster. Replace the previous std::map<int, ...>
  // (RB-tree, O(log n) per op) with a vector indexed by root id (O(1) per
  // op). Sparse — only root ids are populated — but root ids fit in [0, n)
  // so the address space is dense enough; the constant-factor win on the
  // 2*n accesses (sum-by-root + remap lookup) more than offsets the
  // memory.
  std::vector<vec2> sumPos(n, vec2{0, 0});
  std::vector<int> sumCnt(n, 0);
  for (int i = 0; i < n; ++i) {
    int r = uf.find(i);
    sumPos[r] = sumPos[r] + in[i];
    sumCnt[r] += 1;
  }
  // Assign new indices in ascending root-id order so output ordering is
  // deterministic and matches what the old std::map iteration produced.
  std::vector<int> rootToNew(n, -1);
  std::vector<vec2> verts;
  verts.reserve(n);
  for (int r = 0; r < n; ++r) {
    if (sumCnt[r] == 0) continue;
    rootToNew[r] = static_cast<int>(verts.size());
    verts.push_back(sumPos[r] * (1.0 / sumCnt[r]));
  }
  std::vector<int> remap(n);
  for (int i = 0; i < n; ++i) remap[i] = rootToNew[uf.find(i)];
  return {std::move(remap), std::move(verts)};
}

// vertEdges[v] (filled by step 4) and adj[v] (filled by step 3) both
// hold a small set of int ids per vertex. Almost always 2-4 elements;
// occasionally larger at concurrent intersection points. A sorted
// std::vector<int> beats a std::set<int> by 5-10x on per-op cost for
// sets this small (no node allocation, no tree rebalancing, contiguous
// memory). Helpers keep the "set" semantics: idempotent insert, fast
// contains, ordered iteration.
inline bool VESetContains(const std::vector<int>& vec, int x) {
  return std::binary_search(vec.begin(), vec.end(), x);
}
inline void VESetInsert(std::vector<int>* vec, int x) {
  auto it = std::lower_bound(vec->begin(), vec->end(), x);
  if (it == vec->end() || *it != x) vec->insert(it, x);
}

// =============================================================================
// Step 2: collapse edges whose endpoints map to the same vertex.
// =============================================================================
std::vector<EdgeM> RemapAndCollapse(const std::vector<EdgeM>& edges,
                                    const std::vector<int>& remap) {
  std::vector<EdgeM> out;
  out.reserve(edges.size());
  for (const auto& e : edges) {
    int a = remap[e.v0];
    int b = remap[e.v1];
    if (a != b) out.push_back({a, b, e.mult});
  }
  return out;
}

// =============================================================================
// Step 3: per-edge ordered list of vertices within eps of the edge interior.
// Returns vertList[edgeIdx] = sorted list of vert indices along the edge.
// =============================================================================
// `edgeBoxes` and `bvh` are the eps-padded segment AABBs and the BVH
// built over them; they are passed in from the caller so step 3 and
// step 4 can share a single build (the edges array doesn't change
// between them, so the boxes don't either).
//
// (Tried two alternative broad phases here, both regressed:
//   1. Uniform grid: long edges span many cells → excessive
//      candidate counts even after per-vert dedup.
//   2. Interior-only BVH (boxes shrunk by eps along the edge to
//      exclude endpoints): the eps perpendicular padding cancels the
//      along-shrink in the AABB, so endpoints stay inside the box and
//      no callbacks are filtered. AABB can't cleanly separate
//      endpoint from "valid vert on edge near endpoint" without
//      OBB-style geometry.)
std::vector<std::vector<int>> BuildEdgeVertLists(
    const std::vector<EdgeM>& edges, const std::vector<vec2>& verts,
    double eps, const std::vector<Box>& edgeBoxes, const BVH& bvh) {
  const int nE = static_cast<int>(edges.size());
  const int nV = static_cast<int>(verts.size());
  const double eps2 = eps * eps;
  std::vector<std::vector<int>> lists(nE);
  // BVH broad phase: edges as eps-padded segment AABBs, queried by vert
  // points (eps-padded boxes). Each candidate (edge, vert) pair runs the
  // exact projection test below. Per-edge `hits` are sorted by parameter
  // at the end so the result is independent of broad-phase visit order.
  // (edgeBoxes and bvh are now passed in from the caller — see signature.)
  std::vector<Box> vertBoxes(nV);
  for (int v = 0; v < nV; ++v) vertBoxes[v] = BoxOf2DPoint(verts[v], eps);
  // Build vert→neighbors adjacency from the input edges. Used below to
  // detect "thin triangle apex" cases (see comment in CollidePairs).
  // Sorted vector here too — a polygon edge graph has degree ~2-3 per
  // vert, so std::set is pure overhead.
  std::vector<std::vector<int>> adj(nV);
  for (const auto& e : edges) {
    VESetInsert(&adj[e.v0], e.v1);
    VESetInsert(&adj[e.v1], e.v0);
  }
  // Collect (edge, vert) candidate pairs first; then process per edge.
  // (An earlier two-phase attempt that ran the gate in parallel
  // regressed wall-clock by ~25% on JTS — the extra allocation for the
  // candidate buffer and the per-hit Keep flag outweighed the parallel
  // gain because most cases have <1e4 candidates and stay serial under
  // autoPolicy. The serial form is faster across the corpus.)
  // Single flat (edge, t, vert) buffer instead of vector-of-vectors.
  // For typical input most edges have zero hits, so the per-edge vector
  // construction (10170 cases × ~200 edges = ~2M empty-vector
  // constructions) was a real allocator-churn cost. Sorting one flat
  // vector and slicing per-edge groups is faster.
  struct Hit { int e; double t; int v; };
  std::vector<Hit> flatHits;
  // Per-(v, e) narrow-phase test. Captures the same logic whether
  // candidates come from BVH or brute-force broad phase.
  auto narrow = [&](int v, int e) {
    if (v == edges[e].v0 || v == edges[e].v1) return;
    // Thin-triangle-apex skip: when V is connected to BOTH edge endpoints
    // by other edges, V is the apex of a triangle (V, e.v0, e.v1) whose
    // base is this edge. With non-tiny eps (large displacement), the apex
    // can fall within eps of its base; without the skip, step 5
    // canonicalization cancels the apex-split sub-edges against the
    // triangle's other two sides, producing empty output. Only-one-
    // endpoint adjacency is normal polygon-neighbor configuration, so
    // we require BOTH to be conservative.
    if (VESetContains(adj[v], edges[e].v0) &&
        VESetContains(adj[v], edges[e].v1))
      return;
    const vec2 a = verts[edges[e].v0];
    const vec2 b = verts[edges[e].v1];
    const vec2 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) return;
    const vec2 ap = verts[v] - a;
    const double t = dot(ap, ab) / abLen2;
    if (t <= 0 || t >= 1) return;
    const vec2 closest = a + ab * t;
    const vec2 d = verts[v] - closest;
    if (dot(d, d) <= eps2) flatHits.push_back({e, t, v});
  };
  if (bvh.leafToOrig.empty()) {
    // Brute-force broad phase: O(V·E). Faster than building+querying
    // a BVH for small E (caller decided E < threshold).
    for (int v = 0; v < nV; ++v) {
      for (int e = 0; e < nE; ++e) {
        if (vertBoxes[v].DoesOverlap(edgeBoxes[e])) narrow(v, e);
      }
    }
  } else {
    CollidePairs(bvh, vertBoxes, [&](int v, int e) { narrow(v, e); });
  }
  // Sort flat hits by (edge, t) and emit per-edge subsequences.
  std::sort(flatHits.begin(), flatHits.end(), [](const Hit& a, const Hit& b) {
    if (a.e != b.e) return a.e < b.e;
    return a.t < b.t;
  });
  for (size_t i = 0; i < flatHits.size();) {
    const int e = flatHits[i].e;
    size_t j = i;
    while (j < flatHits.size() && flatHits[j].e == e) ++j;
    lists[e].reserve(j - i);
    for (size_t k = i; k < j; ++k) lists[e].push_back(flatHits[k].v);
    i = j;
  }
  return lists;
}

// =============================================================================
// Step 4: edge-edge intersection discovery; insert new verts (or snap to
// existing nearby verts) into per-edge lists.
// =============================================================================
// Step 4 broad phase only: find candidate (i, j) edge pairs whose AABBs
// overlap. Output is sorted lex-ascending. Extracted as a separate
// function so it can run concurrently with step 3 (BuildEdgeVertLists)
// via tbb::parallel_invoke.
//
// For BVH path, uses Collider's parallel-capable Collisions API with a
// thread-local accumulator (tbb::combinable). Each thread emits pairs
// into its own vector; combine_each merges them, then a single sort
// re-establishes lex order. This is faster than the previous serial
// CollidePairs adapter when the candidate count justifies the parallel
// overhead.
#if (MANIFOLD_PAR == 1)
namespace step4_recorder {
struct PairsRecorder {
  using Local = std::vector<std::pair<int, int>>;
  const std::vector<int>& leafToOrig;
  tbb::combinable<Local> tls;
  inline void record(int queryIdx, int leafIdx, Local& local) const {
    const int li = leafToOrig[leafIdx];
    if (queryIdx < li) local.emplace_back(queryIdx, li);
  }
  Local& local() { return tls.local(); }
};
}  // namespace step4_recorder
#endif

inline std::vector<std::pair<int, int>> CollectStep4Pairs(
    const std::vector<EdgeM>& edges, const std::vector<Box>& edgeBoxes,
    const BVH& bvh) {
  const int nE = static_cast<int>(edges.size());
  std::vector<std::pair<int, int>> pairs;
  if (bvh.leafToOrig.empty()) {
    for (int i = 0; i < nE; ++i) {
      for (int j = i + 1; j < nE; ++j) {
        if (edgeBoxes[i].DoesOverlap(edgeBoxes[j])) pairs.emplace_back(i, j);
      }
    }
    return pairs;
  }
#if (MANIFOLD_PAR == 1)
  // Same threshold story as parallel_invoke: only enable Collider's
  // parallel mode when the work justifies the per-thread accumulator +
  // combine_each overhead. Below the threshold, use the serial adapter.
  if (nE >= 256) {
    step4_recorder::PairsRecorder rec{bvh.leafToOrig, {}};
    auto qf = [&](int i) { return edgeBoxes[i]; };
    bvh.collider.Collisions<false>(rec, qf, nE, /*parallel=*/true);
    rec.tls.combine_each([&](const auto& localPairs) {
      pairs.insert(pairs.end(), localPairs.begin(), localPairs.end());
    });
    std::sort(pairs.begin(), pairs.end());
  } else {
    CollidePairs(bvh, edgeBoxes, [&](int qi, int li) {
      if (qi >= li) return;
      pairs.emplace_back(qi, li);
    });
    std::sort(pairs.begin(), pairs.end());
  }
#else
  CollidePairs(bvh, edgeBoxes, [&](int qi, int li) {
    if (qi >= li) return;
    pairs.emplace_back(qi, li);
  });
  std::sort(pairs.begin(), pairs.end());
#endif
  return pairs;
}

void FindAndInsertIntersections(const std::vector<EdgeM>& edges,
                                std::vector<vec2>* verts,
                                std::vector<std::vector<int>>* lists,
                                std::vector<std::vector<int>>* vertEdges,
                                double eps,
                                const std::vector<Box>& edgeBoxes,
                                const BVH& bvh,
                                const std::vector<std::pair<int, int>>& pairs) {
  using Clock = std::chrono::steady_clock;
  auto Ns = [](Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count();
  };
  auto tNarrowStart = Clock::now();
  const int nE = static_cast<int>(edges.size());
  (void)nE;
  const double eps2 = eps * eps;
  vertEdges->resize(verts->size());
  const int origNumVerts = static_cast<int>(verts->size());
  // Broad phase (BVH self-collide → `pairs`) is now caller-supplied so it
  // can run concurrently with step 3 via tbb::parallel_invoke. The
  // narrow-phase loop below is order-dependent (snap-to-existing reads
  // `lists[*]` mutated by earlier pairs), so it stays serial.
  // (An earlier two-phase split — parallel IntersectSegments compute,
  // serial bookkeeping — measured as a wash on JTS: the per-call
  // allocation of intersection-point buffers cancelled the parallel
  // gain because most cases have small pair counts that stay serial
  // under autoPolicy. Kept the simple sequential form.)
  for (auto [i, j] : pairs) {
    // Skip if shares an endpoint.
    if (edges[i].v0 == edges[j].v0 || edges[i].v0 == edges[j].v1 ||
        edges[i].v1 == edges[j].v0 || edges[i].v1 == edges[j].v1)
      continue;
    vec2 a0 = (*verts)[edges[i].v0];
    vec2 a1 = (*verts)[edges[i].v1];
    vec2 b0 = (*verts)[edges[j].v0];
    vec2 b1 = (*verts)[edges[j].v1];
    vec2 p;
    if (!IntersectSegments(a0, a1, b0, b1, edges[i].v0, edges[i].v1,
                           edges[j].v0, edges[j].v1, eps, &p))
      continue;
    // Snap: is p within eps of any existing vert? Search the union of
    // (i,j)'s endpoints and existing list members of i and j.
    auto nearVert = [&](int candidate) -> bool {
      vec2 d = p - (*verts)[candidate];
      return dot(d, d) <= eps2;
    };
    int snapTo = -1;
    for (int v : {edges[i].v0, edges[i].v1, edges[j].v0, edges[j].v1}) {
      if (nearVert(v)) {
        snapTo = v;
        break;
      }
    }
    if (snapTo < 0) {
      for (int v : (*lists)[i]) {
        if (nearVert(v)) {
          snapTo = v;
          break;
        }
      }
    }
    if (snapTo < 0) {
      for (int v : (*lists)[j]) {
        if (nearVert(v)) {
          snapTo = v;
          break;
        }
      }
    }
    int vNew;
    if (snapTo >= 0) {
      vNew = snapTo;
    } else {
      vNew = static_cast<int>(verts->size());
      verts->push_back(p);
      vertEdges->emplace_back();
    }
    // Record edge incidence: this vert is now known to lie on edges i and j.
    VESetInsert(&(*vertEdges)[vNew], i);
    VESetInsert(&(*vertEdges)[vNew], j);
    // Insert into both edges' lists, sorted by parametric position.
    auto insertSorted = [&](int eIdx) {
      if (vNew == edges[eIdx].v0 || vNew == edges[eIdx].v1) return;
      vec2 a = (*verts)[edges[eIdx].v0];
      vec2 b = (*verts)[edges[eIdx].v1];
      vec2 ab = b - a;
      double abLen2 = dot(ab, ab);
      if (abLen2 == 0) return;
      double tNew = dot(p - a, ab) / abLen2;
      auto& lst = (*lists)[eIdx];
      auto pos = std::lower_bound(
          lst.begin(), lst.end(), tNew, [&](int v, double t) {
            double tv = dot((*verts)[v] - a, ab) / abLen2;
            return tv < t;
          });
      if (pos == lst.end() || *pos != vNew) lst.insert(pos, vNew);
    };
    insertSorted(i);
    insertSorted(j);
  }

  // Eager propagation pass (2D analog of boolean_result.cpp::AddNewEdgeVerts).
  // The pair-by-pair loop above only adds each new intersection vert to the
  // two edges in its source pair. But k>=3 input edges can be concurrent at
  // one true point: e.g., edges A, B, C all crossing at P generates three
  // pair-intersections (A,B), (A,C), (B,C); the snap-to-existing logic
  // catches duplicates within eps when they share an edge in the iteration
  // order, but misses the "no shared edge" case (A,B vs E,F both at P).
  // Without propagation, edge C ends up with no on-edge vert at P even though
  // P geometrically lies on C; step 5 then emits C as one un-split sub-edge
  // and the resulting DCEL has the wrong incidence at P.
  //
  // For each freshly-allocated intersection vert, query the edge BVH for all
  // edges within eps and add the vert to their on-edge lists. The query box
  // is eps-padded; the per-candidate exact projection-distance gate matches
  // step 3's. Step 4b's structural+geometric merge stays as a fallback for
  // FP noise that exceeds eps between near-duplicate intersections (which
  // this propagation does not address).
  // Fast-path: if no new intersection verts were created in the main
  // loop, the propagation pass has nothing to do.
  auto tNarrowEnd = Clock::now();
  GlobalPhases().step4NarrowNs.fetch_add(Ns(tNarrowStart, tNarrowEnd),
                                          std::memory_order_relaxed);
  if ((int)verts->size() == origNumVerts) return;
  auto tPropStart = Clock::now();
  std::vector<Box> queryBoxes;
  std::vector<int> queryVerts;
  queryBoxes.reserve((int)verts->size() - origNumVerts);
  queryVerts.reserve((int)verts->size() - origNumVerts);
  for (int v = origNumVerts; v < (int)verts->size(); ++v) {
    queryBoxes.push_back(BoxOf2DPoint((*verts)[v], eps));
    queryVerts.push_back(v);
  }
  // Per-(qi, eIdx) propagation step. Same logic for BVH and brute-force
  // broad phases.
  auto propagateNarrow = [&](int qi, int eIdx) {
    const int v = queryVerts[qi];
    if (v == edges[eIdx].v0 || v == edges[eIdx].v1) return;
    if (VESetContains((*vertEdges)[v], eIdx)) return;  // already incident
    const vec2 a = (*verts)[edges[eIdx].v0];
    const vec2 b = (*verts)[edges[eIdx].v1];
    const vec2 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) return;
    const vec2 p = (*verts)[v];
    const double t = dot(p - a, ab) / abLen2;
    if (t <= 0 || t >= 1) return;
    const vec2 closest = a + ab * t;
    const vec2 d = p - closest;
    if (dot(d, d) > eps2) return;
    auto& lst = (*lists)[eIdx];
    auto pos = std::lower_bound(
        lst.begin(), lst.end(), t, [&](int vv, double tQ) {
          double tv = dot((*verts)[vv] - a, ab) / abLen2;
          return tv < tQ;
        });
    if (pos == lst.end() || *pos != v) lst.insert(pos, v);
    VESetInsert(&(*vertEdges)[v], eIdx);
  };
  if (bvh.leafToOrig.empty()) {
    for (size_t qi = 0; qi < queryBoxes.size(); ++qi) {
      for (int e = 0; e < nE; ++e) {
        if (queryBoxes[qi].DoesOverlap(edgeBoxes[e]))
          propagateNarrow((int)qi, e);
      }
    }
  } else {
    CollidePairs(bvh, queryBoxes,
                 [&](int qi, int eIdx) { propagateNarrow(qi, eIdx); });
  }
  auto tPropEnd = Clock::now();
  GlobalPhases().step4PropNs.fetch_add(Ns(tPropStart, tPropEnd),
                                        std::memory_order_relaxed);
}

// =============================================================================
// Step 5: break edges into sub-edges; merge duplicates with multiplicity sum.
// Smith's PolySet2 (Table 7.3): map<lex-ordered segment, signed multiplicity>.
// =============================================================================
struct CanonEdge {
  int vMin, vMax;
  int mult;
};

// Canonical sub-edges, stored as a sorted-by-(vMin,vMax) flat vector so
// the multiple downstream sequential scans (step 6's halfedge build,
// FastEdge build, and final filter) hit contiguous memory instead of
// walking a std::map RB-tree three times. Build-side: append unsorted
// during `Canonicalize`, then `Finalize` sorts and merges duplicates by
// summing signed mults (matches the older map-insert-with-fold behavior).
struct CanonicalSubEdges {
  std::vector<CanonEdge> edges;  // sorted by (vMin, vMax) after Finalize

  inline void Add(int v0, int v1, int mult) {
    if (v0 == v1) return;
    int vMin = std::min(v0, v1);
    int vMax = std::max(v0, v1);
    int signedMult = (v0 < v1) ? mult : -mult;
    edges.push_back({vMin, vMax, signedMult});
  }

  // Sort by (vMin, vMax), merge consecutive duplicates by summing mults,
  // drop any whose summed mult is zero.
  void Finalize() {
    std::sort(edges.begin(), edges.end(),
              [](const CanonEdge& a, const CanonEdge& b) {
                if (a.vMin != b.vMin) return a.vMin < b.vMin;
                return a.vMax < b.vMax;
              });
    size_t w = 0;
    for (size_t r = 0; r < edges.size();) {
      size_t k = r;
      int sumMult = 0;
      while (k < edges.size() && edges[k].vMin == edges[r].vMin &&
             edges[k].vMax == edges[r].vMax) {
        sumMult += edges[k].mult;
        ++k;
      }
      if (sumMult != 0) {
        edges[w] = {edges[r].vMin, edges[r].vMax, sumMult};
        ++w;
      }
      r = k;
    }
    edges.resize(w);
  }
};

CanonicalSubEdges Canonicalize(const std::vector<EdgeM>& edges,
                               const std::vector<std::vector<int>>& lists) {
  CanonicalSubEdges out;
  // Pre-reserve. Each input edge contributes (1 + lists[e].size()) sub-edges.
  size_t total = edges.size();
  for (const auto& l : lists) total += l.size();
  out.edges.reserve(total);
  for (size_t e = 0; e < edges.size(); ++e) {
    int prev = edges[e].v0;
    for (int v : lists[e]) {
      out.Add(prev, v, edges[e].mult);
      prev = v;
    }
    out.Add(prev, edges[e].v1, edges[e].mult);
  }
  out.Finalize();
  return out;
}

// =============================================================================
// Step 6: winding-number filter via per-edge midpoint ray-cast.
// =============================================================================

// Cast a horizontal ray to +x from `origin`. Count signed crossings of
// directed sub-edges. For each sub-edge (vMin, vMax) with signed
// multiplicity m, the contribution is +m if the edge goes upward (in
// canonical (vMin, vMax) form, "upward" = vMin's y < vMax's y when
// crossed from below), -m if downward.
//
// Used per-face by `FilterByWindingDCEL`: the origin is offset
// perpendicularly LEFT of a boundary half-edge, which puts it strictly
// inside the face and away from any vertex, so the ray never hits a
// vertex exactly under non-adversarial inputs.
int CastWindingRay(vec2 origin, const CanonicalSubEdges& canon,
                   const std::vector<vec2>& verts) {
  int winding = 0;
  for (const auto& edge : canon.edges) {
    const int mult = edge.mult;
    vec2 p0 = verts[edge.vMin];
    vec2 p1 = verts[edge.vMax];
    // Order so p0.y <= p1.y for crossing test.
    bool upward = p0.y < p1.y;
    if (!upward) std::swap(p0, p1);
    // Strictly half-open in y to avoid double-counting at vertices.
    if (origin.y < p0.y || origin.y >= p1.y) continue;
    // Compute x of the segment at origin.y.
    double t = (origin.y - p0.y) / (p1.y - p0.y);
    double xCross = p0.x + t * (p1.x - p0.x);
    if (xCross < origin.x) continue;
    // Crossing direction: original direction was upward iff key.first <
    // key.second and positions matched -- already encoded in `mult`'s sign. We
    // need the signed contribution to the winding number when crossing
    // left-to-right. For a positive-multiplicity edge oriented (vMin -> vMax)
    // in canonical form, an upward crossing (with origin to the left of the
    // edge) increments the winding number on the right side by mult. Since we
    // cast +x, we are computing winding on the LEFT of the ray, which is the
    // side we're at the origin. This is +mult per upward crossing, -mult per
    // downward crossing.
    winding += upward ? mult : -mult;
  }
  return winding;
}

// Pre-flattened canonical sub-edges for the per-face ray-cast loop. Step
// 6 calls CastWindingRay once per face, each call iterating canon.map
// (an RB-tree). For F faces and E canonical edges that's F·E tree
// traversals. Hoisting to a flat array of (yMin, yMax, x0, dx_dy,
// signedMult) once before the parallel for_each turns the inner loop
// into a tight cache-friendly scan and removes per-call indirection
// through canon.map / verts[].
//
// signedMult encodes both direction and magnitude: positive when the
// canonical (vMin → vMax) form is upward (p0.y < p1.y), negative when
// downward. The ray-cast then just sums signedMult over crossed edges.
struct FastEdge {
  double yMin, yMax;       // canonical orientation: yMin < yMax
  double x0;               // x at y = yMin
  double dxdy;             // (x1 - x0) / (yMax - yMin)
  int signedMult;          // mult if upward in canonical form, else -mult
};

inline std::vector<FastEdge> BuildFastEdges(const CanonicalSubEdges& canon,
                                            const std::vector<vec2>& verts) {
  std::vector<FastEdge> out;
  out.reserve(canon.edges.size());
  for (const auto& edge : canon.edges) {
    const int mult = edge.mult;
    vec2 p0 = verts[edge.vMin];
    vec2 p1 = verts[edge.vMax];
    const bool upward = p0.y < p1.y;
    if (!upward) std::swap(p0, p1);
    if (p0.y == p1.y) continue;  // horizontal edges contribute nothing
    FastEdge e;
    e.yMin = p0.y;
    e.yMax = p1.y;
    e.x0 = p0.x;
    e.dxdy = (p1.x - p0.x) / (p1.y - p0.y);
    e.signedMult = upward ? mult : -mult;
    out.push_back(e);
  }
  return out;
}

inline int CastWindingRayFast(vec2 origin,
                              const std::vector<FastEdge>& edges) {
  int winding = 0;
  for (const auto& e : edges) {
    if (origin.y < e.yMin || origin.y >= e.yMax) continue;
    const double xCross = e.x0 + e.dxdy * (origin.y - e.yMin);
    if (xCross < origin.x) continue;
    winding += e.signedMult;
  }
  return winding;
}

// =============================================================================
// Step 6: planar face traversal (DCEL).
//
// Builds a DCEL (doubly-connected edge list, the same structure
// manifold's 3D mesh `Manifold::Impl::halfedge_` uses) from the
// canonical sub-edges, walks face cycles to identify each planar face,
// and assigns ONE winding number per face by ray-casting from a point
// in each face's interior. An edge is kept iff its left-face and
// right-face windings disagree on the >0 predicate. Per-vertex
// consistency is structural: all edges incident to a vertex see the
// same face windings on each side, so no per-edge ray-cast can disagree
// with its neighbors.
//
// Complexity: O(E log E) for the per-vertex angular sort + O(F · E)
// for the per-face winding ray-casts (F = face count). On the deepfuzz
// 14000 cases, this is 0 first-pass topology FAILs vs 270 for the
// per-edge ray-cast approach this replaces.
// =============================================================================

namespace dcel_internal {
struct HalfEdge {
  int twin;    // index of the twin half-edge in halfedges[]
  int next;    // next half-edge along the same face's CCW boundary
  int origin;  // vertex this half-edge starts at
  int face;    // face id (-1 until assigned)
  int mult;    // signed multiplicity in this direction
};
}  // namespace dcel_internal

// Internal predicate over a face's winding number, deciding whether
// the face is "inside" the result region. Three predicates are needed:
//   - WindAdd:       w > 0     (default; Smith's wind > 0 union)
//   - WindIntersect: w > 1     (Boolean2D::Intersect, both inputs cover)
//   - WindEvenOdd:   w & 1     (used internally by Xor)
// An edge is retained iff its left and right faces disagree on the
// predicate. Kept as an internal mechanism rather than a public
// FillRule enum: CrossSection's existing API doesn't expose a fill
// rule (Simplify is (eps), Boolean is (other, OpType)), so a public
// FillRule on this prototype would be API expansion beyond the
// landing target. Boolean2D / Xor / Simplify pick the right
// predicate internally based on the operation requested.
using WindPredicate = std::function<bool(int)>;

inline WindPredicate WindAdd() { return [](int w) { return w > 0; }; }
inline WindPredicate WindIntersect() { return [](int w) { return w > 1; }; }
inline WindPredicate WindEvenOdd() {
  return [](int w) { return (w & 1) != 0; };
}

std::vector<OutEdge> FilterByWindingDCEL(
    const CanonicalSubEdges& canon, const std::vector<vec2>& verts,
    bool debug = false, const WindPredicate& isInside = WindAdd()) {
  using dcel_internal::HalfEdge;
  if (debug) {
    std::cerr << "[FilterByWindingDCEL] canon.edges.size()="
              << canon.edges.size() << " verts.size()=" << verts.size()
              << "\n";
  }
  // 1. Build half-edges. Each canonical (vMin, vMax) with mult m becomes:
  //    - hA: vMin → vMax, mult = m
  //    - hB: vMax → vMin, mult = -m
  //    Twins are paired (hA.twin = hB, hB.twin = hA).
  std::vector<HalfEdge> halfedges;
  halfedges.reserve(2 * canon.edges.size());
  for (const auto& edge : canon.edges) {
    int hA = static_cast<int>(halfedges.size());
    halfedges.push_back({hA + 1, -1, edge.vMin, -1, edge.mult});
    halfedges.push_back({hA, -1, edge.vMax, -1, -edge.mult});
  }
  if (halfedges.empty()) return {};

  // 2. Group half-edges by origin vertex; sort each group by direction angle
  //    (CCW). This is the "rotational order" needed to compute next pointers.
  //    Vert ids are dense [0, verts.size()), so a vector-of-vector indexed
  //    by vert id beats std::map on cache locality and lookup cost.
  std::vector<std::vector<int>> outgoing(verts.size());
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    outgoing[halfedges[i].origin].push_back(i);
  }
  // Per-vertex angular sort: each vertex's outgoing list is independent
  // (read-only on halfedges/verts; writes only its own slot). atan2 is
  // expensive and the inner sort is O(d log d) per vertex of degree d;
  // for big arrangements the total is the second-largest serial cost
  // inside step 6. Output is deterministic because the sort is pure.
  manifold::for_each(
      manifold::autoPolicy(outgoing.size()), manifold::countAt(0),
      manifold::countAt(outgoing.size()), [&](size_t v) {
        auto& hes = outgoing[v];
        if (hes.size() < 2) return;
        const vec2 vp = verts[v];
        // atan2-free angular comparator: split the plane into two half-
        // planes (bucket 0 = upper + +x axis, bucket 1 = lower + -x
        // axis); within a bucket, compare by sign of the cross product.
        // Sorts CCW from +x. Same monotone order as atan2 but no
        // transcendental in the per-comparison hot path.
        auto bucket = [](const vec2& d) {
          return (d.y > 0 || (d.y == 0 && d.x > 0)) ? 0 : 1;
        };
        std::sort(hes.begin(), hes.end(), [&](int a, int b) {
          const vec2 dA = verts[halfedges[halfedges[a].twin].origin] - vp;
          const vec2 dB = verts[halfedges[halfedges[b].twin].origin] - vp;
          const int bA = bucket(dA), bB = bucket(dB);
          if (bA != bB) return bA < bB;
          return dA.x * dB.y - dA.y * dB.x > 0;
        });
      });

  // 3. Compute next pointers. For half-edge h arriving at vertex v
  //    (= h.twin.origin), h.next must be the outgoing edge that makes
  //    the SMALLEST LEFT TURN from h's incoming direction. h.incoming
  //    direction is opposite of h.twin's outgoing direction; the smallest
  //    CCW rotation from h.incoming visits half-edges starting from
  //    "h.incoming + small CCW" and finds the first entry. In the sorted
  //    CCW outgoing list, this corresponds to **one step CW** from
  //    h.twin (with wraparound).
  //
  //    Equivalent: starting at angle (h.twin + π) = h.incoming, sweep CCW
  //    by ε, look up the first sorted entry. That entry is at sorted-list
  //    position one before h.twin (= it - 1, with wraparound).
  //
  //    Using "it+1" instead of "it-1" picks the half-edge ALMOST A FULL
  //    REVOLUTION CCW from h.twin = a RIGHT turn at v. For degree-2
  //    vertices (chains, simple polygon corners), N=2 symmetry makes
  //    "it+1" and "it-1" equivalent and both work. For degree-≥3
  //    vertices (intersection points after step 4), they differ and
  //    only "it-1" produces correctly-oriented face cycles.
  // Each halfedge's .next is determined independently by reading the
  // (now-sorted) outgoing list at its destination vertex. Writes are to
  // independent slots; reads are read-only. Tried caching position-of-
  // self in an int-per-halfedge array — for the typical degree-2
  // polygon-corner case the std::find on a 2-element vector is faster
  // than building the cache, so kept the find.
  manifold::for_each(
      manifold::autoPolicy(halfedges.size()), manifold::countAt(0),
      manifold::countAt(static_cast<int>(halfedges.size())), [&](int i) {
        const int twinIdx = halfedges[i].twin;
        const int destV = halfedges[twinIdx].origin;
        auto& sorted = outgoing[destV];
        auto it = std::find(sorted.begin(), sorted.end(), twinIdx);
        if (it == sorted.end()) return;
        auto prevIt = (it == sorted.begin()) ? (sorted.end() - 1) : (it - 1);
        halfedges[i].next = *prevIt;
      });

  // 4. Walk face cycles, assign face IDs. Each unmarked half-edge starts a
  //    new face; follow `next` chain back to the start.
  int nFaces = 0;
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    if (halfedges[i].face != -1) continue;
    int h = i;
    int safety = 0;
    do {
      if (halfedges[h].next == -1 || safety++ > (int)halfedges.size()) {
        // Malformed cycle; bail rather than infinite-loop.
        break;
      }
      halfedges[h].face = nFaces;
      h = halfedges[h].next;
    } while (h != i);
    ++nFaces;
  }

  // 5. Compute signed area per face. Outer face has the most negative area.
  //    CRITICAL: at displaced coords (e.g. 1.5e9), the raw shoelace
  //    `a.x * b.y - b.x * a.y` has each product on order 2.25e18 with
  //    ULP ~2000, so summation precision swamps a typical face area of
  //    O(1). The sign becomes random and outer-face detection breaks.
  //    Fix: center each face's coordinates relative to its first vertex
  //    before summing. With centered coords O(edge length), products
  //    are O((edge length)²) and the sum is precise.
  // First half-edge encountered per face. Used both as the centering
  // reference for the shoelace area (below) and as the starting half-edge
  // for the per-face winding ray-cast (step 6).
  std::vector<int> faceStartHE(nFaces, -1);
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    if (halfedges[i].face >= 0 && faceStartHE[halfedges[i].face] == -1)
      faceStartHE[halfedges[i].face] = i;
  }
  std::vector<double> faceArea(nFaces, 0.0);
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    if (halfedges[i].face < 0) continue;
    const int faceRefHE = faceStartHE[halfedges[i].face];
    if (faceRefHE < 0) continue;
    const vec2 ref = verts[halfedges[faceRefHE].origin];
    const vec2 a = verts[halfedges[i].origin] - ref;
    const vec2 b = verts[halfedges[halfedges[i].twin].origin] - ref;
    faceArea[halfedges[i].face] += (a.x * b.y - b.x * a.y) * 0.5;
  }
  int outerFace = 0;
  for (int f = 1; f < nFaces; ++f) {
    if (faceArea[f] < faceArea[outerFace]) outerFace = f;
  }
  if (debug) {
    std::cerr << "DCEL: " << halfedges.size() << " halfedges, " << nFaces
              << " faces\n";
    int negAreaCount = 0;
    for (int f = 0; f < nFaces; ++f) {
      std::cerr << "  face " << f << " area=" << faceArea[f]
                << (f == outerFace ? "  <-- outer" : "") << "\n";
      if (faceArea[f] < 0 && f != outerFace) ++negAreaCount;
    }
    if (negAreaCount > 0) {
      std::cerr << "  WARNING: " << negAreaCount
                << " bounded face(s) have negative signed area; cycle "
                   "convention may be inverted\n";
    }
    // Group half-edges by face, count mults.
    std::map<int, std::map<int, int>> faceMults;
    for (int i = 0; i < (int)halfedges.size(); ++i) {
      faceMults[halfedges[i].face][halfedges[i].mult]++;
    }
    for (auto& [f, m] : faceMults) {
      std::cerr << "  face " << f << " mults:";
      for (auto& [mu, c] : m) std::cerr << " " << mu << "x" << c;
      std::cerr << "\n";
    }
  }

  // 6. Compute winding per face by propagation through twin pointers.
  //
  //    Earlier this did F independent ray-casts (one per face), which
  //    is O(F·E) work. Replace with: ray-cast a single seed face per
  //    connected component of the DCEL, then BFS-propagate to all
  //    other faces in the component. Stepping from face_A to face_B
  //    across a halfedge h (h.face=A, h.twin.face=B) crosses h's
  //    canonical-edge contribution, so faceWind[B] = faceWind[A] -
  //    h.mult (we're stepping from LEFT of h, where +mult contributes,
  //    to RIGHT, where 0 does). BFS reaches every face in the
  //    component exactly once, so total work is O(E + F) instead of
  //    O(F·E).
  //
  //    Multi-component arrangements (real for self-intersecting input
  //    whose result has multiple disjoint regions): when BFS finishes,
  //    any unvisited face is in another component; ray-cast its
  //    interior to seed and BFS again. Most cases are single-component
  //    so this fallback rarely fires.
  //
  //    We keep the FastEdge hoist (still needed for the seed ray-casts)
  //    but the F-fold dependency on it is gone — propagation is a
  //    constant per face.
  const auto fastEdges = BuildFastEdges(canon, verts);
  std::vector<int> faceWind(nFaces, 0);
  std::vector<uint8_t> wAssigned(nFaces, 0);
  auto seedRayCast = [&](int f) {
    int h = faceStartHE[f];
    if (h < 0) {
      faceWind[f] = 0;
      wAssigned[f] = 1;
      return;
    }
    const vec2 a = verts[halfedges[h].origin];
    const vec2 b = verts[halfedges[halfedges[h].twin].origin];
    const vec2 mid = (a + b) * 0.5;
    const vec2 d = b - a;
    const double len = length(d);
    if (len == 0) {
      faceWind[f] = 0;
    } else {
      const vec2 perp(-d.y / len, d.x / len);
      const vec2 pInF = mid + perp * (len * 1e-3);
      faceWind[f] = CastWindingRayFast(pInF, fastEdges);
    }
    wAssigned[f] = 1;
  };
  // Outer face is the convention-fixed seed: by topology it sits
  // outside the arrangement, so faceWind = 0. The shoelace area test
  // identified it above.
  //
  // Don't ray-cast — at large coord magnitudes (e.g. JTS GIS data with
  // coords ~4.5e6 and polygon extents only ~80), the perpendicular step
  // (`len * 1e-3`) can be small enough relative to coordinate ULPs that
  // the FP-rounded `pInF` doesn't reliably land in the *outer* half-
  // plane of the boundary edge — the ray-cast then returns a non-zero
  // winding and every other face's BFS-propagated winding inherits the
  // bias. By topology the outer face is unbounded, so winding=0 is the
  // correct convention regardless of input orientation.
  faceWind[outerFace] = 0;
  wAssigned[outerFace] = 1;
  std::vector<int> bfsQ;
  bfsQ.reserve(nFaces);
  auto propagateFrom = [&](int seed) {
    bfsQ.clear();
    bfsQ.push_back(seed);
    size_t head = 0;
    while (head < bfsQ.size()) {
      const int f = bfsQ[head++];
      const int h0 = faceStartHE[f];
      if (h0 < 0) continue;
      int hh = h0;
      int safety = 0;
      do {
        const int twinH = halfedges[hh].twin;
        const int adj = halfedges[twinH].face;
        if (adj >= 0 && !wAssigned[adj]) {
          // Stepping LEFT of hh (= f) → RIGHT of hh (= adj):
          // winding loses the +mult contribution that the LEFT side saw.
          faceWind[adj] = faceWind[f] - halfedges[hh].mult;
          wAssigned[adj] = 1;
          bfsQ.push_back(adj);
        }
        hh = halfedges[hh].next;
        if (hh < 0 || ++safety > (int)halfedges.size()) break;
      } while (hh != h0);
    }
  };
  propagateFrom(outerFace);
  // Pick up any disconnected components with their own seed ray-cast.
  for (int f = 0; f < nFaces; ++f) {
    if (!wAssigned[f]) {
      seedRayCast(f);
      propagateFrom(f);
    }
  }

  if (debug) {
    std::cerr << "  face windings:";
    for (int f = 0; f < nFaces; ++f) {
      std::cerr << " f" << f << "=" << faceWind[f];
    }
    std::cerr << "\n";
  }

  // 7. Filter canonical sub-edges by left/right face windings. The first
  //    half-edge of each pair (the (vMin → vMax) direction) is at index
  //    2*i; its twin (vMax → vMin) is at 2*i + 1.
  std::vector<OutEdge> out;
  out.reserve(canon.edges.size());
  int hi = 0;
  for (const auto& edge : canon.edges) {
    const int hA = hi;
    const int hB = hi + 1;
    hi += 2;
    const int leftFace = halfedges[hA].face;
    const int rightFace = halfedges[hB].face;
    if (leftFace < 0 || rightFace < 0) continue;
    const int wL = faceWind[leftFace];
    const int wR = faceWind[rightFace];
    const bool leftIn = isInside(wL);
    const bool rightIn = isInside(wR);
    if (leftIn == rightIn) continue;
    if (leftIn) {
      out.push_back({edge.vMin, edge.vMax, 1});
    } else {
      out.push_back({edge.vMax, edge.vMin, 1});
    }
  }
  return out;
}

// =============================================================================
// End-to-end driver.
// =============================================================================
struct OverlapResult {
  std::vector<vec2> verts;
  std::vector<OutEdge> edges;
  std::vector<int> inputRemap;  // input vert idx -> output vert idx
  int numMergedVerts;           // count of verts before step-4 intersections
};

// Per-phase wall-clock accumulator. Atomic so concurrent RemoveOverlaps2D
// calls (e.g. from a multi-threaded corpus driver, future) can share.
// Currently RemoveOverlaps2D itself is single-threaded; the accumulator
// captures cumulative time across many cases for the `time` CLI mode.
// Defined earlier in the file — see the PhaseAcc declaration above the
// step 1 helpers.

OverlapResult RemoveOverlaps2D(const std::vector<vec2>& vertsIn,
                               const std::vector<EdgeM>& edgesIn, double eps,
                               bool debug = false,
                               const WindPredicate& isInside = WindAdd()) {
  using timing_detail::Clock;
  using timing_detail::Ns;
  auto& P = GlobalPhases();
  const auto tStart = Clock::now();
  // Step 1: vertex merge.
  auto t0 = Clock::now();
  auto merge = MergeVerts(vertsIn, eps);
  auto t1 = Clock::now();
  P.mergeNs.fetch_add(Ns(t0, t1), std::memory_order_relaxed);
  const int numMerged = static_cast<int>(merge.verts.size());
  // Step 2: collapse edges.
  auto edges = RemapAndCollapse(edgesIn, merge.remap);
  auto t2 = Clock::now();
  P.remapNs.fetch_add(Ns(t1, t2), std::memory_order_relaxed);
  // Build the edge BVH once for steps 3 and 4. For very small edge
  // counts (< 32), skip the BVH build — Morton sort + tree
  // construction overhead beats brute-force O(n²) on tiny inputs. The
  // step 3/4 functions detect an empty bvh.leafToOrig and fall back to
  // their brute-force broad phase.
  std::vector<Box> edgeBoxes(edges.size());
  for (size_t e = 0; e < edges.size(); ++e) {
    edgeBoxes[e] =
        BoxOf2DEdge(merge.verts[edges[e].v0], merge.verts[edges[e].v1], eps);
  }
  auto tBvhStart = Clock::now();
  BVH bvh;
  if (edges.size() >= 32) bvh = BVHBuildFromBoxes(edgeBoxes);
  auto tBvhEnd = Clock::now();
  P.bvhBuildNs.fetch_add(Ns(tBvhStart, tBvhEnd), std::memory_order_relaxed);
  // Steps 3 (BuildEdgeVertLists) and 4-broad (CollectStep4Pairs) are
  // independent: step 3 writes `lists`, step 4-broad writes `pairs`,
  // both read-only on (edges, verts, edgeBoxes, bvh). Run them in
  // parallel via tbb::parallel_invoke when the work is large enough to
  // amortize the parallel-task overhead — TBB pays a few µs to spin up
  // a parallel region, which can exceed the gain on tiny cases. The
  // E≥256 threshold was tuned empirically on a contended cloud
  // machine (where blind parallelism regressed wall-clock vs serial).
  std::vector<std::vector<int>> lists;
  std::vector<std::pair<int, int>> step4Pairs;
  // Sub-phase timing: time step 3 (full) and step 4 broad separately
  // even when they run concurrently — measures the total CPU time, not
  // wall-clock, so the sum of the two roughly equals the work done.
  auto tStep3Start = Clock::now();
#if MANIFOLD_PAR == 1
  if (edges.size() >= 256) {
    tbb::parallel_invoke(
        [&] {
          lists = BuildEdgeVertLists(edges, merge.verts, eps, edgeBoxes, bvh);
        },
        [&] { step4Pairs = CollectStep4Pairs(edges, edgeBoxes, bvh); });
  } else {
    lists = BuildEdgeVertLists(edges, merge.verts, eps, edgeBoxes, bvh);
    step4Pairs = CollectStep4Pairs(edges, edgeBoxes, bvh);
  }
#else
  lists = BuildEdgeVertLists(edges, merge.verts, eps, edgeBoxes, bvh);
  step4Pairs = CollectStep4Pairs(edges, edgeBoxes, bvh);
#endif
  // Charge to step3BroadNs as a stand-in for the combined step 3 work
  // (broad + narrow are interleaved inside BuildEdgeVertLists' callback;
  // splitting them would require a structural change). step 4 broad
  // phase is included in this same wall-clock segment because of
  // parallel_invoke; the breakdown is approximate.
  auto tStep3End = Clock::now();
  P.step3BroadNs.fetch_add(Ns(tStep3Start, tStep3End),
                            std::memory_order_relaxed);
  auto t3 = Clock::now();
  P.buildListsNs.fetch_add(Ns(t2, t3), std::memory_order_relaxed);
  // Step 4 narrow phase: order-dependent serial pass.
  std::vector<std::vector<int>> vertEdges;
  FindAndInsertIntersections(edges, &merge.verts, &lists, &vertEdges, eps,
                             edgeBoxes, bvh, step4Pairs);
  auto t4 = Clock::now();
  P.findIxNs.fetch_add(Ns(t3, t4), std::memory_order_relaxed);

  // Step 4b: structural re-merge of intersection verts. Step 4 inserts each
  // intersection at the time its parent edge pair is processed; if pairs
  // (A, B) and (A, C) both produce intersections at the same true point P*
  // (i.e. three edges meet there) they may land FP-close but not snap to
  // each other because neither saw the other yet at insertion time.
  //
  // Two intersections that should be the same true point share at least one
  // common edge in their incidence list (any two of {AB, AC, BC} share an
  // edge). Two intersections from disjoint edge sets cannot be the same
  // true point, regardless of geometric distance. So we union-find verts
  // that share an edge AND fall within eps; this avoids the angle-dependent
  // threshold of a pure-geometric merge (a fixed factor like 1.5*eps fails
  // for shallow crossings; large factors over-merge legitimately-distinct
  // intersections from unrelated edge pairs).
  //
  // Caveat: 4+ edges concurrent at one point produces some intersection
  // pairs (e.g. AB vs CD) that share no edge. Such configurations are
  // adversarial and not covered here; iterate-to-fixed-point catches some
  // but not all cases.
  {
    DisjointSets uf(static_cast<int>(merge.verts.size()));
    // The geometric upper bound for "same true point" is eps/sin(theta)
    // where theta is the crossing angle. For shallow crossings this can
    // be large; we use a generous 10*eps cutoff which covers theta down
    // to ~6 degrees. The structural gate prevents over-merging
    // legitimately-distinct intersections (e.g. edge A crosses B at one
    // point and C at a different point along A: vAB and vAC share edge A
    // but are at different true points and shouldn't merge unless they
    // ALSO geometrically coincide). A sweep across the displacement fuzz
    // showed 10*eps gives the best iteration count (1:448 2:2) without
    // over-merging; tightening below 3*eps causes single-pass failures,
    // loosening to 100*eps causes new over-merge failures.
    //
    // Note: when union-find creates a multi-vert cluster, the centroid
    // computed below is offset from the original positions by up to ~eps;
    // that's intentional (we WANT the merged point to land at the average)
    // and is the source of the residual iter=2 cases. Smith's bound
    // proves convergence in ≤2 iterations under his α-budget framework.
    const double mergeThresh = 10.0 * eps;
    const double mergeThresh2 = mergeThresh * mergeThresh;
    // BVH broad phase over intersection verts only (those with non-empty
    // edge incidence). Each box is padded by `mergeThresh` so any pair
    // within geometric distance ≤ mergeThresh has overlapping boxes.
    // Candidate pairs then run the structural gate (shared-edge
    // incidence) and the exact distance gate.
    std::vector<int> ixVerts;
    std::vector<Box> ixBoxes;
    ixVerts.reserve(merge.verts.size());
    ixBoxes.reserve(merge.verts.size());
    for (size_t i = 0; i < merge.verts.size(); ++i) {
      if (i < vertEdges.size() && !vertEdges[i].empty()) {
        ixVerts.push_back(static_cast<int>(i));
        ixBoxes.push_back(BoxOf2DPoint(merge.verts[i], mergeThresh));
      }
    }
    BVH bvh = BVHBuildFromBoxes(ixBoxes);
    std::vector<std::pair<int, int>> pairs;
    CollidePairs(bvh, ixBoxes, [&](int qi, int li) {
      if (qi >= li) return;
      pairs.emplace_back(ixVerts[qi], ixVerts[li]);
    });
    std::sort(pairs.begin(), pairs.end());
    // Parallelize the gate checks (read-only on vertEdges/merge.verts,
    // independent per pair); collect pass/fail into a flag vector. Then
    // unite() serially in the sorted pair order so the choice of cluster
    // root (rank-based merge) is deterministic regardless of thread
    // scheduling. The serial-unite cost is trivial vs the gate cost (set
    // intersection + sqrt-free distance test).
    std::vector<uint8_t> doUnite(pairs.size(), 0);
    manifold::for_each_n(
        manifold::autoPolicy(pairs.size()), manifold::countAt(0),
        pairs.size(), [&](size_t i) {
          const auto [a, b] = pairs[i];
          bool shared = false;
          for (int e : vertEdges[a]) {
            if (VESetContains(vertEdges[b], e)) { shared = true; break; }
          }
          if (!shared) return;
          vec2 d = merge.verts[b] - merge.verts[a];
          if (dot(d, d) > mergeThresh2) return;
          doUnite[i] = 1;
        });
    for (size_t i = 0; i < pairs.size(); ++i) {
      if (doUnite[i]) uf.unite(pairs[i].first, pairs[i].second);
    }
    // Build remap from union-find clusters; cluster position is centroid.
    // Vector-by-root-id replaces std::map (same fix as MergeVerts).
    const int nV = static_cast<int>(merge.verts.size());
    std::vector<vec2> sumPos(nV, vec2{0, 0});
    std::vector<int> sumCnt(nV, 0);
    int distinctClusters = 0;
    for (int i = 0; i < nV; ++i) {
      int r = uf.find(i);
      if (sumCnt[r] == 0) ++distinctClusters;
      sumPos[r] = sumPos[r] + merge.verts[i];
      sumCnt[r] += 1;
    }
    if (distinctClusters < nV) {
      std::vector<int> rootToNew(nV, -1);
      std::vector<vec2> newVerts;
      newVerts.reserve(distinctClusters);
      for (int r = 0; r < nV; ++r) {
        if (sumCnt[r] == 0) continue;
        rootToNew[r] = static_cast<int>(newVerts.size());
        newVerts.push_back(sumPos[r] * (1.0 / sumCnt[r]));
      }
      std::vector<int> remap(nV);
      for (int i = 0; i < nV; ++i) remap[i] = rootToNew[uf.find(i)];
      // Apply remap to edges + lists + composed input remap.
      for (auto& e : edges) {
        e.v0 = remap[e.v0];
        e.v1 = remap[e.v1];
      }
      for (auto& list : lists) {
        for (auto& v : list) v = remap[v];
        list.erase(std::unique(list.begin(), list.end()), list.end());
      }
      for (auto& r : merge.remap) r = remap[r];
      merge.verts = std::move(newVerts);
    }
  }

  auto t4b = Clock::now();
  P.restructNs.fetch_add(Ns(t4, t4b), std::memory_order_relaxed);
  // Step 5: sub-edge canonicalization.
  auto canon = Canonicalize(edges, lists);
  auto t5 = Clock::now();
  P.canonNs.fetch_add(Ns(t4b, t5), std::memory_order_relaxed);
  // Step 6: DCEL face-traversal winding filter.
  auto out = FilterByWindingDCEL(canon, merge.verts, debug, isInside);
  auto t6 = Clock::now();
  P.filterDcelNs.fetch_add(Ns(t5, t6), std::memory_order_relaxed);
  P.totalNs.fetch_add(Ns(tStart, t6), std::memory_order_relaxed);
  P.cases.fetch_add(1, std::memory_order_relaxed);
  return {std::move(merge.verts), std::move(out), std::move(merge.remap),
          numMerged};
}

// =============================================================================
// Polygons API wrapper. The public-facing shape matches manifold's
// `CrossSection` / `Polygons` typedef (`std::vector<std::vector<vec2>>`).
// Each `SimplePolygon` is a closed loop of vertices: CCW outer, CW holes.
// Internally we flatten to a vert/edge list, run the overlap-removal
// pipeline, then walk the directed output edges back into closed loops.
// =============================================================================

// Flatten manifold::Polygons into the lower-level (verts, edges) input.
// Each loop becomes a sequence of edges (v0→v1, v1→v2, …, v_{n-1}→v_0)
// with mult=+1 each. Smith's wind = ±1 convention then assigns the right
// sign for CCW outer (interior wind=+1) vs CW hole (hole-interior wind=0,
// surrounding-polygon-interior wind=+1).
inline std::pair<std::vector<vec2>, std::vector<EdgeM>>
PolygonsToInput(const manifold::Polygons& polys) {
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  for (const auto& loop : polys) {
    // Drop degenerate loops (fewer than 3 verts can't form a closed
    // simple polygon under Smith's wind > 0 convention). Silent skip
    // matches `Polygons` consumer expectations elsewhere in manifold;
    // callers that want hard validation should pre-check.
    if (loop.size() < 3) continue;
    const int base = static_cast<int>(verts.size());
    const int n = static_cast<int>(loop.size());
    for (const auto& v : loop) verts.push_back(v);
    for (int i = 0; i < n; ++i) {
      edges.push_back({base + i, base + ((i + 1) % n), 1});
    }
  }
  return {std::move(verts), std::move(edges)};
}

// Walk the directed sub-edges of an OverlapResult into closed loops.
// Output is **regularized** in the Requicha-Tilove (1978) sense: zero-
// area features (lens-shaped 2-vert loops where two oriented sub-edges
// trace the same line segment in opposite directions, or 1-vert
// degenerate loops) are dropped because they can't be represented in
// `manifold::Polygons` (= `vector<vector<vec2>>`, a list of CCW outer
// + CW hole loops with no way to encode a 1D feature without an
// enclosing face). Matches CGAL `Polygon_set_2`, Clipper2, and SVG
// fill-rule conventions; consumers that need non-regularized output
// require a richer type (CGAL `Arrangement_2`, Clipper2 `PolyTree64`).
//
// At a vertex of degree ≥4 (e.g., an X-cross between two triangles in a
// figure-8 boundary), arbitrarily picking "any unvisited outgoing edge"
// would jump between distinct loops. Same DCEL convention as step 6:
// the next outgoing edge that continues the same loop is the one
// **immediately CW from the incoming half-edge's reverse direction** in
// the vertex's CCW-sorted angular order, i.e., "smallest left turn"
// from the incoming direction.
inline manifold::Polygons OutEdgesToPolygons(
    const std::vector<vec2>& verts, const std::vector<OutEdge>& edges) {
  const int nE = static_cast<int>(edges.size());
  // Per-vertex outgoing edges, sorted CCW by direction angle. Vert ids
  // are dense in [0, verts.size()), so a vector-of-vector indexed by
  // vert id beats std::map on cache locality and lookup cost.
  std::vector<std::vector<int>> outgoing(verts.size());
  for (int i = 0; i < nE; ++i) outgoing[edges[i].v0].push_back(i);
  // atan2-free angular comparator (see step 6 sort site for derivation).
  auto bucket = [](const vec2& d) {
    return (d.y > 0 || (d.y == 0 && d.x > 0)) ? 0 : 1;
  };
  for (size_t v = 0; v < outgoing.size(); ++v) {
    auto& lst = outgoing[v];
    if (lst.size() < 2) continue;
    const vec2 vp = verts[v];
    std::sort(lst.begin(), lst.end(), [&](int a, int b) {
      const vec2 da = verts[edges[a].v1] - vp;
      const vec2 db = verts[edges[b].v1] - vp;
      const int bA = bucket(da), bB = bucket(db);
      if (bA != bB) return bA < bB;
      return da.x * db.y - da.y * db.x > 0;
    });
  }
  std::vector<bool> visited(nE, false);
  manifold::Polygons polys;
  for (int start = 0; start < nE; ++start) {
    if (visited[start]) continue;
    manifold::SimplePolygon loop;
    int cur = start;
    while (cur >= 0 && !visited[cur]) {
      visited[cur] = true;
      loop.push_back(verts[edges[cur].v0]);
      const int destV = edges[cur].v1;
      if (destV < 0 || destV >= (int)outgoing.size() ||
          outgoing[destV].empty()) {
        cur = -1;
        break;
      }
      // To continue the same loop at destV, pick the outgoing edge
      // that's "smallest left turn" from cur's incoming direction.
      // The incoming direction is `verts[edges[cur].v0] → verts[destV]`,
      // i.e., angle = atan2(verts[destV].y - verts[edges[cur].v0].y,
      //                     verts[destV].x - verts[edges[cur].v0].x).
      // The next outgoing edge whose direction is just past this
      // (going CCW) is the one we want. Linear scan with the atan2
      // delta-from-rev calculation — kept as-is because the walk-step
      // runs few times per case (avg ~50 steps per RemoveOverlaps2D
      // output) and atan2 here is not a measurable hotspot vs the
      // angular-sort comparator above (which executes n log n times
      // per sort, called once per output vertex).
      const vec2 vp = verts[destV];
      const vec2 inDir = vp - verts[edges[cur].v0];
      const double inAngle = std::atan2(inDir.y, inDir.x);
      double rev = inAngle + M_PI;
      if (rev > M_PI) rev -= 2 * M_PI;
      const auto& lst = outgoing[destV];
      int next = -1;
      double bestDelta = std::numeric_limits<double>::infinity();
      for (int e : lst) {
        if (visited[e]) continue;
        const vec2 d = verts[edges[e].v1] - vp;
        double ang = std::atan2(d.y, d.x);
        double delta = ang - rev;
        if (delta <= 0) delta += 2 * M_PI;
        if (delta < bestDelta) {
          bestDelta = delta;
          next = e;
        }
      }
      cur = next;
    }
    if (loop.size() >= 3) {
      polys.push_back(std::move(loop));
    } else {
      // Regularization drop: the loop is a zero-area degenerate (1-vert
      // self-loop or 2-vert lens). With straight-line-segment edges,
      // both cases enclose zero area; drop matches CGAL/Clipper2/SVG
      // convention. The assert exists to flag if a future change ever
      // produces a positive-area sub-3-vert loop, which would be an
      // upstream bug.
      assert(std::fabs(SignedArea(loop)) < 1e-12 &&
             "regularized-drop loop should have zero area");
    }
  }
  return polys;
}

// Single-input regularization. Matches `CrossSection::Simplify(eps)`
// (the eventual landing target) exactly: one input, one eps, returns
// the canonical wind > 0 boundary. No public fill-rule parameter,
// since CrossSection's existing API has none.
inline manifold::Polygons Simplify(const manifold::Polygons& in, double eps) {
  auto [verts, edges] = PolygonsToInput(in);
  if (verts.empty()) return {};
  auto r = RemoveOverlaps2D(verts, edges, eps);
  return OutEdgesToPolygons(r.verts, r.edges);
}

// Use manifold's existing public `OpType` from <manifold/common.h> rather
// than defining a parallel local enum. `OpType::{Add, Subtract, Intersect}`
// is what manifold's 3D `Manifold::Boolean` and 2D `CrossSection::Boolean`
// both already accept; aligning here means callers can pass the same
// enum across 3D and 2D codepaths.
using manifold::OpType;

// Infer eps from a polygon set's bounding-box half-extent via Smith's
// α-budget formula. Returns 0 for empty input. Used by Boolean2D when
// the caller passes eps <= 0.
inline double InferEps(const manifold::Polygons& a,
                       const manifold::Polygons& b) {
  double xMin = std::numeric_limits<double>::infinity();
  double yMin = xMin, xMax = -xMin, yMax = -xMin;
  auto bound = [&](const manifold::Polygons& polys) {
    for (const auto& loop : polys) {
      for (const auto& v : loop) {
        xMin = std::min(xMin, v.x);
        yMin = std::min(yMin, v.y);
        xMax = std::max(xMax, v.x);
        yMax = std::max(yMax, v.y);
      }
    }
  };
  bound(a);
  bound(b);
  if (!std::isfinite(xMin)) return 0.0;
  const double L = std::max(std::max(std::fabs(xMin), std::fabs(xMax)),
                            std::max(std::fabs(yMin), std::fabs(yMax)));
  return EpsilonFromScale(L);
}

// Binary boolean. Combines A and B into a single edge set with B's
// multiplicity flipped for Subtract, then applies the same overlap-
// removal pipeline with an op-specific winding predicate at step 6:
//
//   - Add(A, B):       w > 0       (any input covers the face)
//   - Subtract(A, B):  w > 0       (A covers but B's flipped mult cancels)
//   - Intersect(A, B): w > 1       (BOTH inputs cover the face)
//
// Pass `eps <= 0` to auto-infer eps from the combined bounding box.
inline manifold::Polygons Boolean2D(const manifold::Polygons& a,
                                    const manifold::Polygons& b, OpType op,
                                    double eps = 0.0) {
  if (eps <= 0.0) eps = InferEps(a, b);
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  auto append = [&](const manifold::Polygons& polys, int mult) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i) {
        edges.push_back({base + i, base + ((i + 1) % n), mult});
      }
    }
  };
  append(a, 1);
  append(b, op == OpType::Subtract ? -1 : 1);
  if (verts.empty()) return {};
  WindPredicate pred = (op == OpType::Intersect) ? WindIntersect() : WindAdd();
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, pred);
  return OutEdgesToPolygons(r.verts, r.edges);
}

// Symmetric difference (XOR): the region covered by A or B but not both.
// `manifold::OpType` only has three values (Add/Subtract/Intersect), so
// XOR doesn't fit through Boolean2D. Implementation: combine A and B
// with mult=+1 each, then apply EvenOdd predicate at step 6. Regions
// covered once (only A or only B) have w=1 (odd, kept); regions covered
// twice (both) have w=2 (even, dropped); regions covered zero times
// have w=0 (even, dropped). Equivalent to Clipper2's `Xor` clip type
// and the wiring target for `CrossSection::operator^`.
inline manifold::Polygons Xor(const manifold::Polygons& a,
                              const manifold::Polygons& b, double eps = 0.0) {
  if (eps <= 0.0) eps = InferEps(a, b);
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  auto append = [&](const manifold::Polygons& polys) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i) {
        edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
    }
  };
  append(a);
  append(b);
  if (verts.empty()) return {};
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, WindEvenOdd());
  return OutEdgesToPolygons(r.verts, r.edges);
}

// =============================================================================
// Iterate-to-fixed-point.
//
// Smith §7.7 / figure 7.16: residual eps-scale edge intersections from
// rounded arithmetic in the first pass are resolved by re-applying the
// algorithm. We iterate up to maxIter additional times beyond the initial
// pass. Termination:
//   - Converged: new fingerprint == previous fingerprint (no change).
//   - Cycle: new fingerprint == any earlier (non-immediate) fingerprint.
//     Pick the iteration with the lex-smallest fingerprint as canonical
//     (deterministic choice among the cycle's equivalent outputs).
//   - MaxedOut: hit maxIter without convergence/cycle; return current.
//
// Fingerprint quantizes vert positions to multiples of eps/100 so that
// cross-iteration vert renumbering doesn't break comparison.
// =============================================================================
// Quantum-parameterized fingerprint. `Fingerprint` (default) uses eps/100;
// `CoarseFingerprint` (eps quantum) captures only topology, i.e., two
// runs differing by sub-eps position drift will have identical coarse
// fingerprints, which is the right test for "did pass 2 actually change
// the topology, or just shift positions by a few ULPs?"
inline std::string FingerprintAt(const OverlapResult& r, double quantum) {
  auto q = [quantum](double x) {
    return static_cast<int64_t>(std::round(x / quantum));
  };
  std::vector<std::tuple<int64_t, int64_t, int64_t, int64_t, int>> subs;
  subs.reserve(r.edges.size());
  for (const auto& oe : r.edges) {
    vec2 p0 = r.verts[oe.v0];
    vec2 p1 = r.verts[oe.v1];
    auto ka = std::make_pair(q(p0.x), q(p0.y));
    auto kb = std::make_pair(q(p1.x), q(p1.y));
    int mult = oe.mult;
    if (kb < ka) {
      std::swap(ka, kb);
      mult = -mult;
    }
    subs.emplace_back(ka.first, ka.second, kb.first, kb.second, mult);
  }
  std::sort(subs.begin(), subs.end());
  std::ostringstream oss;
  for (const auto& [ax, ay, bx, by, m] : subs) {
    oss << ax << "," << ay << "," << bx << "," << by << ":" << m << ";";
  }
  return oss.str();
}

inline std::string CoarseFingerprint(const OverlapResult& r, double eps) {
  return FingerprintAt(r, eps);
}

// Fine-grained fingerprint used for idempotence detection. The eps/100
// factor lets two iterations at the same coarse-eps tolerance still
// produce distinguishable fingerprints if their snap decisions landed
// at different sub-eps positions; CoarseFingerprint above lumps any
// such pair together. Empirical: anything finer than ~eps/100 starts
// catching pure FP-noise differences.
inline std::string Fingerprint(const OverlapResult& r, double eps) {
  return FingerprintAt(r, eps * 0.01);
}

enum class IterStatus {
  Converged,
  Cycled,
  MaxedOut,
};

// Smith §7.7 / fig 7.16 proves convergence in ≤2 iterations under his
// α-budget framework when intersection positions are tracked symbolically.
// With FP-rounded positions out of `IntersectSegments`, ~1-3% of cases
// need iter=2 for fingerprint convergence; topology is correct from pass
// 1. The production default is 2 (Smith's bound). Tests and the deepfuzz
// pass higher maxIter explicitly to measure tail behavior.
OverlapResult IterateToFixedPoint(const std::vector<vec2>& vIn,
                                  const std::vector<EdgeM>& eIn, double eps,
                                  int maxIter = 2, int* outIters = nullptr,
                                  IterStatus* outStatus = nullptr) {
  std::vector<OverlapResult> history;
  std::vector<std::string> fps;
  history.push_back(RemoveOverlaps2D(vIn, eIn, eps));
  fps.push_back(Fingerprint(history.back(), eps));
  // composedRemap[orig_input_vert] = current iteration's vert idx. Updated
  // each iteration so callers can validate the final result against the
  // original input. Without this, only first-pass `inputRemap` is meaningful.
  std::vector<int> composedRemap = history.back().inputRemap;
  for (int iter = 1; iter <= maxIter; ++iter) {
    std::vector<EdgeM> nextEdges;
    nextEdges.reserve(history.back().edges.size());
    for (const auto& oe : history.back().edges)
      nextEdges.push_back({oe.v0, oe.v1, oe.mult});
    auto next = RemoveOverlaps2D(history.back().verts, nextEdges, eps);
    // Compose: orig→prev_iter via composedRemap, then prev_iter→next via
    // next.inputRemap.
    for (auto& v : composedRemap) v = next.inputRemap[v];
    auto nextFp = Fingerprint(next, eps);
    if (nextFp == fps.back()) {
      if (outIters) *outIters = iter;
      if (outStatus) *outStatus = IterStatus::Converged;
      next.inputRemap = std::move(composedRemap);
      return next;
    }
    // Cycle detection: same fingerprint seen earlier (not just last).
    for (size_t k = 0; k + 1 < fps.size(); ++k) {
      if (fps[k] == nextFp) {
        if (outIters) *outIters = iter;
        if (outStatus) *outStatus = IterStatus::Cycled;
        // Lex-smallest fingerprint wins. The remap is correct for the
        // returned iteration only if we picked `next`; for older history
        // entries we can't reconstruct the composed remap (we threw away
        // intermediate composedRemap states). The caller passing
        // numMergedVerts and edges remains consistent for `next`.
        size_t minIdx = 0;
        for (size_t j = 1; j < fps.size(); ++j) {
          if (fps[j] < fps[minIdx]) minIdx = j;
        }
        if (nextFp < fps[minIdx]) {
          next.inputRemap = std::move(composedRemap);
          return next;
        }
        return std::move(history[minIdx]);
      }
    }
    history.push_back(std::move(next));
    fps.push_back(std::move(nextFp));
  }
  if (outIters) *outIters = maxIter;
  if (outStatus) *outStatus = IterStatus::MaxedOut;
  history.back().inputRemap = std::move(composedRemap);
  return std::move(history.back());
}

// =============================================================================
// Property checks.
// =============================================================================

// Compute per-vertex (out_mult - in_mult) balance for an edge set.
template <typename Edge>
std::map<int, int> ComputeBalance(const std::vector<Edge>& edges) {
  std::map<int, int> b;
  for (const auto& e : edges) {
    b[e.v0] += e.mult;
    b[e.v1] -= e.mult;
  }
  return b;
}

// Topological invariant: for vertices that exist in BOTH input and output
// (i.e. original input verts surviving the merge), the output balance must
// equal the input balance. New intersection-created verts must have balance
// zero. This generalizes "sum in == sum out" to handle open inputs.
//
// Note: empty output trivially passes; closed input edges contribute
// equal in/out at every vertex, expected balance is 0 everywhere, and an
// empty edge list also has 0 balance everywhere. Callers that care about
// "did we drop everything?" must additionally check `r.edges` non-empty.
//
// numMergedVerts = number of verts emitted by step 1's merge (the first N
// in r.verts are either original or merged-from-original; intersection
// verts come after).
bool CheckTopologicalValidity(const OverlapResult& r,
                              const std::vector<EdgeM>& inputEdges,
                              const std::vector<int>& inputRemap,
                              int numMergedVerts) {
  // Expected balance: remap input edges through `inputRemap`, compute
  // per-output-vertex balance contribution from the input.
  std::vector<EdgeM> remapped;
  for (const auto& e : inputEdges) {
    int a = inputRemap[e.v0];
    int b = inputRemap[e.v1];
    if (a != b) remapped.push_back({a, b, e.mult});
  }
  auto expected = ComputeBalance(remapped);
  auto actual = ComputeBalance(r.edges);
  bool ok = true;
  for (int v = 0; v < static_cast<int>(r.verts.size()); ++v) {
    int exp = expected.count(v) ? expected[v] : 0;
    int act = actual.count(v) ? actual[v] : 0;
    // Original-input verts (idx < numMergedVerts): output balance must
    // match input. Intersection verts (idx >= numMergedVerts): must be 0.
    int target = (v < numMergedVerts) ? exp : 0;
    if (act != target) {
      std::cerr << "[FAIL] vertex " << v << " (orig=" << (v < numMergedVerts)
                << ") expected balance " << target << ", got " << act
                << std::endl;
      ok = false;
    }
  }
  return ok;
}

// =============================================================================
// Test patterns from Smith §7.7.
// =============================================================================

// Random topological polygon: N points on unit circle in random angular
// order. Produces a self-intersecting closed curve.
std::pair<std::vector<vec2>, std::vector<EdgeM>> RandomTopologicalPolygon(
    int n, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u01(0, 1);
  std::vector<vec2> verts(n);
  for (int i = 0; i < n; ++i) {
    double theta = 2 * M_PI * u01(rng);
    verts[i] = {std::cos(theta), std::sin(theta)};
  }
  // Random permutation for connectivity.
  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::shuffle(order.begin(), order.end(), rng);
  std::vector<EdgeM> edges;
  edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    edges.push_back({order[i], order[(i + 1) % n], 1});
  }
  return {std::move(verts), std::move(edges)};
}

// Polygonal star: n equispaced points on unit circle, connected by skipping
// k vertices. Highly self-intersecting and symmetric.
std::pair<std::vector<vec2>, std::vector<EdgeM>> PolygonalStar(int n,
                                                               int skip) {
  std::vector<vec2> verts(n);
  for (int i = 0; i < n; ++i) {
    double theta = 2 * M_PI * i / n;
    verts[i] = {std::cos(theta), std::sin(theta)};
  }
  std::vector<EdgeM> edges;
  edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    edges.push_back({i, (i + skip) % n, 1});
  }
  return {std::move(verts), std::move(edges)};
}

// Apply a translation that brings the bounding box near a power of 2.
// Smith §7.7's "displacement attack" makes representable grid spacing
// comparable to feature sizes.
void Displace(std::vector<vec2>* verts, double offset) {
  for (auto& v : *verts) {
    v.x += offset;
    v.y += offset;
  }
}

// L-shaped concave polygon (CCW). 6 verts, no self-intersection.
std::pair<std::vector<vec2>, std::vector<EdgeM>> LShape() {
  std::vector<vec2> v = {
      {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2},
  };
  std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1},
                          {3, 4, 1}, {4, 5, 1}, {5, 0, 1}};
  return {std::move(v), std::move(e)};
}

// Square (CCW) with a square hole (CW). Two disjoint loops.
std::pair<std::vector<vec2>, std::vector<EdgeM>> SquareWithHole() {
  std::vector<vec2> v = {
      // Outer CCW square
      {0, 0},
      {4, 0},
      {4, 4},
      {0, 4},
      // Inner CW square (hole)
      {1, 1},
      {1, 3},
      {3, 3},
      {3, 1},
  };
  std::vector<EdgeM> e = {
      {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // outer CCW
      {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},  // inner CW
  };
  return {std::move(v), std::move(e)};
}

// Two disjoint CCW squares.
std::pair<std::vector<vec2>, std::vector<EdgeM>> TwoSquares() {
  std::vector<vec2> v = {
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {3, 0}, {3, 1}, {2, 1},
  };
  std::vector<EdgeM> e = {
      {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},
      {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},
  };
  return {std::move(v), std::move(e)};
}

// =============================================================================
// Main: run a battery of tests and report.
// =============================================================================
struct TestCase {
  std::string name;
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  double eps;
};

// Run RemoveOverlaps2D and verify topological invariant. Returns the result
// for further property checks.
OverlapResult RunCase(const TestCase& tc, bool* allPass) {
  std::cout << "=== " << tc.name << " ===" << std::endl;
  std::cout << "  Input: " << tc.verts.size() << " verts, " << tc.edges.size()
            << " edges, eps=" << tc.eps << std::endl;
  auto r = RemoveOverlaps2D(tc.verts, tc.edges, tc.eps);
  std::cout << "  Output: " << r.verts.size() << " verts, " << r.edges.size()
            << " edges" << std::endl;
  bool topo =
      CheckTopologicalValidity(r, tc.edges, r.inputRemap, r.numMergedVerts);
  std::cout << "  Topological validity: " << (topo ? "PASS" : "FAIL")
            << std::endl;
  if (!topo) *allPass = false;
  return r;
}

// Idempotence: feeding the output back through RemoveOverlaps2D should
// produce the same result. Uses Fingerprint (position-based, eps/100
// quantized) for cross-iteration vert-renumbering tolerance, matching
// the displacement fuzz's check.
bool CheckIdempotence(const OverlapResult& first, double eps) {
  std::vector<EdgeM> asEdgeM;
  for (const auto& e : first.edges) {
    asEdgeM.push_back({e.v0, e.v1, e.mult});
  }
  auto second = RemoveOverlaps2D(first.verts, asEdgeM, eps);
  bool ok = (Fingerprint(first, eps) == Fingerprint(second, eps));
  std::cout << "  Idempotence: " << (ok ? "PASS" : "FAIL") << std::endl;
  return ok;
}

}  // namespace overlap2d

// Diagnostic: run ONE failing case and dump intermediate state.
// Invoked via `./overlap2d_proto diagnose <seed> [kPow] [n]`. Default
// kPow=30, n=50, but any DeepFuzz parameter combo can be targeted.
// IMPORTANT: the seed→RNG mapping must match DeepFuzz exactly.
namespace overlap2d {
void Diagnose(uint64_t seed, int kPow = 30, int n = 50) {
  const double offset = std::ldexp(1.5, kPow);
  const double eps = EpsilonFromScale(offset);
  std::cerr << "=== Diagnose: kPow=" << kPow << " n=" << n << " seed=" << seed
            << " offset=" << offset << " eps=" << eps << " ===\n";

  // Match the seed mapping used by both main()'s displacement fuzz and
  // DeepFuzz: `seed + 1000 * kPow`.
  auto [vIn, eIn] = RandomTopologicalPolygon(n, seed + 1000ull * kPow);
  Displace(&vIn, offset);
  std::cerr << "Input: " << vIn.size() << " verts, " << eIn.size()
            << " edges\n";

  // Run pipeline manually to inspect intermediate state.
  auto merge = MergeVerts(vIn, eps);
  std::cerr << "After step 1 (merge): " << merge.verts.size() << " verts (was "
            << vIn.size() << ")\n";

  auto edges = RemapAndCollapse(eIn, merge.remap);
  std::cerr << "After step 2 (collapse): " << edges.size() << " edges (was "
            << eIn.size() << ")\n";

  std::vector<Box> diagEdgeBoxes(edges.size());
  for (size_t i = 0; i < edges.size(); ++i)
    diagEdgeBoxes[i] = BoxOf2DEdge(merge.verts[edges[i].v0],
                                   merge.verts[edges[i].v1], eps);
  BVH diagBvh = BVHBuildFromBoxes(diagEdgeBoxes);
  auto lists =
      BuildEdgeVertLists(edges, merge.verts, eps, diagEdgeBoxes, diagBvh);
  size_t totalListSize = 0;
  for (const auto& l : lists) totalListSize += l.size();
  std::cerr << "After step 3 (lists): " << totalListSize
            << " total list entries across " << edges.size() << " edges\n";

  const int beforeStep4 = static_cast<int>(merge.verts.size());
  std::vector<std::vector<int>> vertEdges;
  auto diagPairs = CollectStep4Pairs(edges, diagEdgeBoxes, diagBvh);
  FindAndInsertIntersections(edges, &merge.verts, &lists, &vertEdges, eps,
                             diagEdgeBoxes, diagBvh, diagPairs);
  const int afterStep4 = static_cast<int>(merge.verts.size());
  std::cerr << "After step 4 (intersections): " << merge.verts.size()
            << " verts (added " << (afterStep4 - beforeStep4) << ")\n";

  // Step 4b (structural merge) is omitted here: production
  // RemoveOverlaps2D uses union-find over verts that share a parent edge
  // and fall within 10*eps. The diagnostic shows the post-step-4
  // (pre-structural-merge) state so any residual near-coincident
  // intersection clusters surface in the dump below.
  std::cerr << "After step 4b (re-merge): skipped (production uses "
               "structural merge; see RemoveOverlaps2D)\n";

  auto canon = Canonicalize(edges, lists);
  std::cerr << "After step 5 (canonicalize): " << canon.edges.size()
            << " sub-edges (after multiplicity collapse)\n";

  // Vertex balance from canonicalized sub-edges.
  std::map<int, int> canonBalance;
  for (const auto& edge : canon.edges) {
    canonBalance[edge.vMin] += edge.mult;
    canonBalance[edge.vMax] -= edge.mult;
  }

  // Expected balance from remapped input.
  std::vector<EdgeM> remappedInput;
  for (const auto& e : eIn) {
    int a = merge.remap[e.v0];
    int b = merge.remap[e.v1];
    if (a != b) remappedInput.push_back({a, b, e.mult});
  }
  auto expBalance = ComputeBalance(remappedInput);

  // Find imbalanced vertices.
  const int numMergedVerts = beforeStep4;
  std::cerr << "\n=== Per-vertex balance check ===\n";
  int badCount = 0;
  for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
    int target =
        (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
    int actual = canonBalance.count(v) ? canonBalance[v] : 0;
    if (actual != target) {
      ++badCount;
      if (badCount <= 8) {
        std::cerr << "  v" << v << " (orig=" << (v < numMergedVerts)
                  << ") pos=(" << merge.verts[v].x << "," << merge.verts[v].y
                  << ") expected=" << target << " actual=" << actual << "\n";
        // Dump sub-edges touching this vertex.
        std::cerr << "    sub-edges touching v" << v << ":\n";
        int touched = 0;
        for (const auto& edge : canon.edges) {
          if (edge.vMin == v || edge.vMax == v) {
            std::cerr << "      (" << edge.vMin << "," << edge.vMax
                      << ") mult=" << edge.mult << "\n";
            if (++touched >= 12) {
              std::cerr << "      ... (more truncated)\n";
              break;
            }
          }
        }
      }
    }
  }
  std::cerr << "Total imbalanced vertices at step 5: " << badCount << "\n";

  // Now also run step 6 (winding filter) and check balance again.
  auto out = FilterByWindingDCEL(canon, merge.verts);
  std::cerr << "\nAfter step 6 (winding filter): " << out.size()
            << " output edges (from " << canon.edges.size()
            << " canonical sub-edges)\n";
  std::map<int, int> outBalance;
  for (const auto& e : out) {
    outBalance[e.v0] += e.mult;
    outBalance[e.v1] -= e.mult;
  }
  int outBadCount = 0;
  for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
    int target =
        (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
    int actual = outBalance.count(v) ? outBalance[v] : 0;
    if (actual != target) {
      ++outBadCount;
      if (outBadCount <= 8) {
        std::cerr << "  v" << v << " (orig=" << (v < numMergedVerts)
                  << ") pos=(" << merge.verts[v].x << "," << merge.verts[v].y
                  << ") expected=" << target << " actual=" << actual << "\n";
        std::cerr << "    output edges touching v" << v << ":\n";
        int touched = 0;
        for (const auto& e : out) {
          if (e.v0 == v || e.v1 == v) {
            std::cerr << "      (" << e.v0 << "->" << e.v1
                      << ") mult=" << e.mult << "\n";
            if (++touched >= 12) {
              std::cerr << "      ... (more truncated)\n";
              break;
            }
          }
        }
      }
    }
  }
  std::cerr << "Total imbalanced vertices at step 6: " << outBadCount << "\n";

  // Deep dive on first imbalanced vertex: dump winding decisions for every
  // canonical sub-edge touching it, including the wL/wR values.
  if (outBadCount > 0) {
    int targetV = -1;
    for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
      int target =
          (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
      int actual = outBalance.count(v) ? outBalance[v] : 0;
      if (actual != target) {
        targetV = v;
        break;
      }
    }
    std::cerr << "\n=== Deep dive on v" << targetV << " ===\n";
    // Also dump nearby vertices and check distances for §7.6-style grouping.
    std::cerr << "  v" << targetV << " neighbors (within 10*eps):\n";
    const double thresh2 = (10 * eps) * (10 * eps);
    for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
      if (v == targetV) continue;
      vec2 d = merge.verts[v] - merge.verts[targetV];
      double dist2 = dot(d, d);
      if (dist2 <= thresh2) {
        std::cerr << "    v" << v << " dist=" << std::sqrt(dist2) << " ("
                  << (std::sqrt(dist2) / eps) << "*eps)\n";
      }
    }
    // Per-edge perpendicular-offset ray-cast for diagnostic visualization
    // only (the actual step 6 uses DCEL face traversal).
    const double ofs = eps * 1e-1;
    for (const auto& edge : canon.edges) {
      if (edge.vMin != targetV && edge.vMax != targetV) continue;
      vec2 p0 = merge.verts[edge.vMin];
      vec2 p1 = merge.verts[edge.vMax];
      vec2 mid = (p0 + p1) * 0.5;
      vec2 d = p1 - p0;
      double len = length(d);
      vec2 perp = {-d.y / len, d.x / len};
      vec2 leftPt = mid + perp * ofs;
      vec2 rightPt = mid + perp * -ofs;
      int wL = CastWindingRay(leftPt, canon, merge.verts);
      int wR = CastWindingRay(rightPt, canon, merge.verts);
      std::cerr << "  sub-edge (" << edge.vMin << "," << edge.vMax
                << ") mult=" << edge.mult << " | mid=(" << mid.x - offset
                << "," << mid.y - offset << ") perp=(" << perp.x << ","
                << perp.y << ") | wL=" << wL << " wR=" << wR
                << (((wL > 0) != (wR > 0)) ? " => KEPT" : " => DROPPED")
                << "\n";
    }
  }

  // Cluster proximity: how many distinct verts are within eps of each other?
  int closeCount = 0;
  const double eps2 = eps * eps;
  for (int i = 0; i < static_cast<int>(merge.verts.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(merge.verts.size()); ++j) {
      vec2 d = merge.verts[i] - merge.verts[j];
      if (dot(d, d) <= eps2) ++closeCount;
    }
  }
  std::cerr << "After step 4b: " << closeCount
            << " pairs of distinct verts within eps of each other\n";

  // Now run iterate-to-fixed-point and report on the iterated result.
  std::cerr << "\n=== IterateToFixedPoint ===\n";
  int iters = 0;
  IterStatus iterStatus = IterStatus::Converged;
  auto rIter = IterateToFixedPoint(vIn, eIn, eps, /*maxIter=*/8, &iters,
                                   &iterStatus);
  const char* statusStr =
      iterStatus == IterStatus::Converged
          ? "Converged"
          : (iterStatus == IterStatus::Cycled ? "Cycled" : "MaxedOut");
  std::cerr << "Iters: " << iters << " status: " << statusStr << "\n";
  std::cerr << "Iterated result: " << rIter.verts.size() << " verts, "
            << rIter.edges.size() << " edges, numMergedVerts="
            << rIter.numMergedVerts << "\n";

  // Check topology validity on the iterated result with composed remap.
  bool iterValid = CheckTopologicalValidity(rIter, eIn, rIter.inputRemap,
                                            rIter.numMergedVerts);
  std::cerr << "Iterated topology valid: " << (iterValid ? "YES" : "NO")
            << "\n";

  // Sanity-check pass 1 the way DeepFuzz does it (the production
  // RemoveOverlaps2D pipeline, NOT the legacy diagnostic above).
  std::cerr << "\n=== Pass 1 (RemoveOverlaps2D, the production path) ===\n";
  auto pass1 = RemoveOverlaps2D(vIn, eIn, eps);
  bool pass1Valid = CheckTopologicalValidity(pass1, eIn, pass1.inputRemap,
                                             pass1.numMergedVerts);
  std::cerr << "Pass 1 verts=" << pass1.verts.size()
            << " edges=" << pass1.edges.size()
            << " numMergedVerts=" << pass1.numMergedVerts << " valid="
            << (pass1Valid ? "YES" : "NO") << "\n";

  // Find imbalanced verts and dump their canonical sub-edge environment.
  if (!pass1Valid) {
    std::cerr << "\n=== Imbalanced vert analysis ===\n";
    // Compute expected balance from input edges remapped through
    // pass1.inputRemap (same as CheckTopologicalValidity).
    std::vector<EdgeM> remapped;
    for (const auto& ie : eIn) {
      int a = pass1.inputRemap[ie.v0], b = pass1.inputRemap[ie.v1];
      if (a != b) remapped.push_back({a, b, ie.mult});
    }
    auto exp = ComputeBalance(remapped);
    std::map<int, int> act;
    for (const auto& oe : pass1.edges) {
      act[oe.v0] += oe.mult;
      act[oe.v1] -= oe.mult;
    }
    for (int v = 0; v < (int)pass1.verts.size(); ++v) {
      int target = (v < pass1.numMergedVerts)
                       ? (exp.count(v) ? exp[v] : 0)
                       : 0;
      int actual = act.count(v) ? act[v] : 0;
      if (actual != target) {
        const vec2 p = pass1.verts[v];
        std::cerr << "  v" << v << " (orig=" << (v < pass1.numMergedVerts)
                  << ") pos=(" << (p.x - offset) << "," << (p.y - offset)
                  << ") expected=" << target << " got=" << actual << "\n";
        // Edges in pass1.edges touching this vert.
        std::cerr << "    output edges touching:\n";
        int tch = 0;
        for (const auto& oe : pass1.edges) {
          if (oe.v0 == v || oe.v1 == v) {
            std::cerr << "      " << oe.v0 << "→" << oe.v1
                      << " mult=" << oe.mult << "\n";
            if (++tch >= 6) break;
          }
        }
      }
    }
    // For each imbalanced vert, find canonical sub-edges incident to it,
    // and check how step 6 voted on each.
    std::cerr << "\n=== Step 6 verdicts for imbalanced verts ===\n";
    // Re-run the production pipeline to canonical sub-edges (structural
    // step 4b included; same as RemoveOverlaps2D up to step 5).
    // We want canon as the production pipeline produces it.
    // Note: pass1.edges is post-step-6; we need post-step-5. So re-run.
    auto m3 = MergeVerts(vIn, eps);
    auto e3 = RemapAndCollapse(eIn, m3.remap);
    std::vector<Box> e3Boxes(e3.size());
    for (size_t i = 0; i < e3.size(); ++i)
      e3Boxes[i] =
          BoxOf2DEdge(m3.verts[e3[i].v0], m3.verts[e3[i].v1], eps);
    BVH e3Bvh = BVHBuildFromBoxes(e3Boxes);
    auto l3 = BuildEdgeVertLists(e3, m3.verts, eps, e3Boxes, e3Bvh);
    std::vector<std::vector<int>> ve3;
    auto e3Pairs = CollectStep4Pairs(e3, e3Boxes, e3Bvh);
    FindAndInsertIntersections(e3, &m3.verts, &l3, &ve3, eps, e3Boxes, e3Bvh,
                               e3Pairs);
    // structural step 4b (copy of the production code)
    {
      DisjointSets uf(static_cast<int>(m3.verts.size()));
      const double mergeThresh2 = (10.0 * eps) * (10.0 * eps);
      for (size_t a = 0; a < m3.verts.size(); ++a) {
        if (a >= ve3.size() || ve3[a].empty()) continue;
        for (size_t b = a + 1; b < m3.verts.size(); ++b) {
          if (b >= ve3.size() || ve3[b].empty()) continue;
          bool shared = false;
          for (int e : ve3[a]) if (VESetContains(ve3[b], e)) { shared = true; break; }
          if (!shared) continue;
          vec2 d = m3.verts[b] - m3.verts[a];
          if (dot(d, d) > mergeThresh2) continue;
          uf.unite((int)a, (int)b);
        }
      }
      std::map<int, std::pair<vec2, int>> sums;
      for (size_t i = 0; i < m3.verts.size(); ++i) {
        int r = uf.find((int)i);
        auto& s = sums[r];
        s.first = s.first + m3.verts[i];
        s.second += 1;
      }
      if (sums.size() < m3.verts.size()) {
        std::map<int, int> rootToNew;
        std::vector<vec2> nv;
        for (const auto& [root, sc] : sums) {
          rootToNew[root] = (int)nv.size();
          nv.push_back(sc.first * (1.0 / sc.second));
        }
        std::vector<int> rmp(m3.verts.size());
        for (size_t i = 0; i < m3.verts.size(); ++i)
          rmp[i] = rootToNew[uf.find((int)i)];
        for (auto& ed : e3) { ed.v0 = rmp[ed.v0]; ed.v1 = rmp[ed.v1]; }
        for (auto& list : l3) {
          for (auto& v : list) v = rmp[v];
          list.erase(std::unique(list.begin(), list.end()), list.end());
        }
        m3.verts = std::move(nv);
      }
    }
    auto canon = Canonicalize(e3, l3);
    std::cerr << "  canonical sub-edges (post step 5): "
              << canon.edges.size() << "\n";

    // For each imbalanced vert, dump every canonical edge touching it and
    // step 6's verdict.
    std::set<int> imbalanced{101, 110, 112, 195};
    const double ofs = eps * 1e-1;
    for (int target : imbalanced) {
      if (target >= (int)m3.verts.size()) continue;
      vec2 vp = m3.verts[target];
      std::cerr << "  v" << target << " pos=(" << (vp.x - offset) << ","
                << (vp.y - offset) << ")\n";
      int sub_count = 0, kept = 0, dropped = 0;
      for (const auto& cedge : canon.edges) {
        if (cedge.vMin != target && cedge.vMax != target) continue;
        ++sub_count;
        const int other = (cedge.vMin == target) ? cedge.vMax : cedge.vMin;
        const vec2 op = m3.verts[other];
        vec2 p0 = vp, p1 = op;
        if (cedge.vMin != target) std::swap(p0, p1);
        vec2 mid = (p0 + p1) * 0.5;
        vec2 d = p1 - p0;
        const double len = length(d);
        if (len == 0) continue;
        vec2 perp = {-d.y / len, d.x / len};
        vec2 lpt = mid + perp * ofs, rpt = mid - perp * ofs;
        int wL = CastWindingRay(lpt, canon, m3.verts);
        int wR = CastWindingRay(rpt, canon, m3.verts);
        bool keep = (wL > 0) != (wR > 0);
        if (keep) ++kept; else ++dropped;
        std::cerr << "    (" << cedge.vMin << "↔" << cedge.vMax
                  << ") mult=" << cedge.mult << " wL=" << wL << " wR=" << wR
                  << (keep ? " KEEP" : " DROP") << " other.pos=("
                  << (op.x - offset) << "," << (op.y - offset) << ")\n";
        if (sub_count >= 8) {
          std::cerr << "    ...\n";
          break;
        }
      }
      std::cerr << "    " << kept << " kept, " << dropped
                << " dropped at this vert\n";
    }
  }
  // Also report what the test is expecting per-vertex.
  std::cerr << "Composed remap (original input vert → final vert idx):\n";
  for (size_t i = 0; i < rIter.inputRemap.size(); ++i) {
    if (i < 10 || i == rIter.inputRemap.size() - 1)
      std::cerr << "  v" << i << " → " << rIter.inputRemap[i] << "\n";
  }
}
}  // namespace overlap2d

// =============================================================================
// Polygon corpus runner.
//
// Reads the manifold project's existing curated test corpus at
// test/polygons/polygon_corpus.txt - 100 named adversarial polygons
// (Eberly, Sliver, Comb, KissingZigzag, BarelyValid, Tricky, Sponge4a,
// CondensedMatter64, etc.) used in the project's triangulation tests.
// Each entry has form (matches test/polygon_test.cpp's
// RegisterPolygonTestsFile):
//
//   <name> <expected_tri_count> <eps> <num_loops>
//   <loop1_vert_count>
//   <x> <y>     (loop1_vert_count times)
//   <loop2_vert_count>   (only if multi-loop)
//   <x> <y>     (loop2_vert_count times)
//   ...
//
// eps = -1.0 means "use a sensible default for this scale". For the
// overlap-removal tests we infer one from the bounding box.
//
// Runs each polygon through Simplify, checks topology
// validity, and verifies signed area is preserved within a tolerance
// (since these are mostly already-valid polygons - any large area drift
// indicates the algorithm modified something it shouldn't have, or
// failed to clean up a self-intersection it should have). The
// expected_tri_count is informational for now (it's a triangulation
// metric, not directly an overlap-removal one).
// =============================================================================
namespace overlap2d {

struct CorpusEntry {
  std::string name;
  manifold::Polygons polys;
  double eps;          // -1 if not specified in file (caller infers)
  int expected_tris;   // triangulation triangle count (informational)
};

inline std::vector<CorpusEntry> LoadCorpus(const std::string& path) {
  std::vector<CorpusEntry> entries;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return entries;
  }
  while (in.good()) {
    CorpusEntry entry;
    int num_loops = 0;
    if (!(in >> entry.name >> entry.expected_tris >> entry.eps >> num_loops))
      break;
    for (int i = 0; i < num_loops; ++i) {
      int loop_size = 0;
      if (!(in >> loop_size)) break;
      manifold::SimplePolygon loop;
      loop.reserve(loop_size);
      for (int j = 0; j < loop_size; ++j) {
        double x, y;
        if (!(in >> x >> y)) break;
        loop.push_back({x, y});
      }
      entry.polys.push_back(std::move(loop));
    }
    entries.push_back(std::move(entry));
  }
  return entries;
}

inline void RunCorpus(const std::string& path) {
  auto entries = LoadCorpus(path);
  std::cout << "=== Polygon corpus: " << entries.size() << " entries from "
            << path << " ===\n";
  int total = 0;
  int topo_ok = 0;
  int topo_fail = 0;
  int area_preserved = 0;     // |drift| < 1%
  int area_drifted = 0;       // 1% <= |drift| < 10%
  int area_collapsed = 0;     // |drift| >= 10%
  int area_input_zero = 0;    // input had zero net area
  std::vector<std::tuple<std::string, double>> driftList;
  for (const auto& entry : entries) {
    ++total;
    double eps = entry.eps > 0 ? entry.eps : InferEps(entry.polys, {});
    if (eps <= 0) eps = EpsilonFromScale(1.0);
    auto out = Simplify(entry.polys, eps);
    // Topology: convert input polys to (verts, edges) form, run pipeline,
    // run the project's CheckTopologicalValidity.
    auto [v, e] = PolygonsToInput(entry.polys);
    auto r = RemoveOverlaps2D(v, e, eps);
    bool topoOk = CheckTopologicalValidity(r, e, r.inputRemap,
                                           r.numMergedVerts);
    if (topoOk) ++topo_ok;
    else ++topo_fail;
    const double inArea = std::fabs(TotalSignedArea(entry.polys));
    const double outArea = std::fabs(TotalSignedArea(out));
    if (inArea < 1e-12) {
      ++area_input_zero;
    } else {
      const double drift = std::fabs(outArea - inArea) / inArea;
      if (drift >= 0.10) {
        ++area_collapsed;
        driftList.emplace_back(entry.name, drift);
      } else if (drift >= 0.01) {
        ++area_drifted;
        driftList.emplace_back(entry.name, drift);
      } else {
        ++area_preserved;
      }
    }
  }
  std::cout << "  Total: " << total << "\n";
  std::cout << "  Topology valid:   " << topo_ok << " / " << total << "\n";
  std::cout << "  Topology INVALID: " << topo_fail << " / " << total << "\n";
  std::cout << "  Area preserved (drift < 1%):       " << area_preserved
            << "\n";
  std::cout << "  Area drifted   (1% <= drift < 10%): " << area_drifted
            << "\n";
  std::cout << "  Area collapsed (drift >= 10%):       " << area_collapsed
            << "\n";
  std::cout << "  Input had zero net area: " << area_input_zero
            << " (degenerate or pure-CW polygons; not an oracle case)\n";
  if (!driftList.empty()) {
    std::sort(driftList.begin(), driftList.end(),
              [](const auto& a, const auto& b) {
                return std::get<1>(a) > std::get<1>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 20):\n";
    for (size_t i = 0; i < driftList.size() && i < 20; ++i) {
      auto& [name, dr] = driftList[i];
      std::cout << "    " << name << "  drift=" << (dr * 100.0) << "%\n";
    }
  }
}

}  // namespace overlap2d

// =============================================================================
// Clipper2 Polygons.txt corpus runner.
//
// Independent third-party oracle from Clipper2's test suite (195 numbered
// cases at build/_deps/clipper2-src/Tests/Polygons.txt). Each case has:
//
//   CAPTION: <n>. [<note>]
//   CLIPTYPE: UNION | INTERSECTION | DIFFERENCE
//   FILLRULE: NONZERO | EVENODD | NEGATIVE | POSITIVE
//   SOL_AREA: <expected absolute area, or -1 if not validated>
//   SOL_COUNT: <expected loop count, or -1>
//   SUBJECTS
//   x,y, x,y, ...        (one polygon per line; integer coords)
//   ...
//   [CLIPS
//   x,y, x,y, ...]
//   <blank line>
//
// Dispatch model: Clipper2's CLIPTYPE × FILLRULE semantics is
//
//   result = boolean_op(fill_rule(subjects), fill_rule(clips))
//
// i.e. each input is independently canonicalized to a CCW arrangement
// under the declared fill rule, then the boolean op combines them. The
// prototype's public Simplify/Boolean2D/Xor APIs assume CCW-canonical
// input (predicates use Smith's w>0 union convention), so they don't
// match Clipper2's semantics directly when inputs are self-intersecting
// or in screen-coord (y-down → CW in y-up math) orientation.
//
// To fairly compare, we run each case through a pre-fill stage:
//
//   - FILLRULE: NONZERO → RemoveOverlaps2D with predicate w != 0
//   - FILLRULE: EVENODD → RemoveOverlaps2D with predicate w & 1
//   - FILLRULE: NEGATIVE/POSITIVE → unsupported (6 cases skipped).
//
// After pre-fill, both sides are CCW-canonical, and we combine via:
//
//   - UNION        → Boolean2D::Add        (w>0 on combined CCW)
//   - INTERSECTION → Boolean2D::Intersect  (w>1 on combined CCW)
//   - DIFFERENCE   → Boolean2D::Subtract   (B mult flipped, then w>0)
//
// The eps inferred from the bounding box (Smith α-budget at L = max
// |coord|) lands at ~1e-12 for these int-coordinate inputs — well below
// the unit grid the corpus uses, so iteration shouldn't be triggered
// by snap-induced changes.
// =============================================================================
namespace overlap2d {

struct Clipper2Case {
  int n = 0;
  std::string caption;
  std::string cliptype;
  std::string fillrule;
  double solArea = -1.0;
  int solCount = -1;
  manifold::Polygons subjects;
  manifold::Polygons clips;
};

// Parse one polygon line "x,y, x,y, x,y" into a SimplePolygon. Coordinates
// are comma-separated; pairs are flat in the list. Tolerates trailing
// whitespace/commas.
inline manifold::SimplePolygon ParseClipper2PolyLine(const std::string& line) {
  manifold::SimplePolygon poly;
  std::vector<double> nums;
  std::string buf;
  for (char c : line) {
    if (c == ',' || std::isspace(static_cast<unsigned char>(c))) {
      if (!buf.empty()) {
        try { nums.push_back(std::stod(buf)); } catch (...) {}
        buf.clear();
      }
    } else {
      buf.push_back(c);
    }
  }
  if (!buf.empty()) {
    try { nums.push_back(std::stod(buf)); } catch (...) {}
  }
  for (size_t i = 0; i + 1 < nums.size(); i += 2) {
    poly.push_back({nums[i], nums[i + 1]});
  }
  return poly;
}

inline std::vector<Clipper2Case> LoadClipper2Cases(const std::string& path) {
  std::vector<Clipper2Case> cases;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return cases;
  }
  Clipper2Case cur;
  bool inCase = false;
  enum class Section { None, Subjects, Clips };
  Section sec = Section::None;
  std::string line;
  auto flush = [&]() {
    if (inCase && !cur.caption.empty()) cases.push_back(std::move(cur));
    cur = Clipper2Case{};
    inCase = false;
    sec = Section::None;
  };
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
      line.pop_back();
    if (line.empty()) {
      flush();
      continue;
    }
    auto starts = [&](const char* k) {
      size_t n = std::strlen(k);
      return line.size() >= n && line.compare(0, n, k) == 0;
    };
    if (starts("CAPTION:")) {
      flush();
      inCase = true;
      cur.caption = line.substr(8);
      // Strip leading whitespace and trailing dot.
      size_t p = 0;
      while (p < cur.caption.size() && std::isspace(static_cast<unsigned char>(cur.caption[p]))) ++p;
      cur.caption = cur.caption.substr(p);
      cur.n = std::atoi(cur.caption.c_str());
      continue;
    }
    if (!inCase) continue;
    if (starts("CLIPTYPE:")) { cur.cliptype = line.substr(9); auto p = cur.cliptype.find_first_not_of(" \t"); cur.cliptype = (p == std::string::npos) ? "" : cur.cliptype.substr(p); continue; }
    if (starts("FILLRULE:")) { cur.fillrule = line.substr(9); auto p = cur.fillrule.find_first_not_of(" \t"); cur.fillrule = (p == std::string::npos) ? "" : cur.fillrule.substr(p); continue; }
    if (starts("SOL_AREA:")) { cur.solArea = std::stod(line.substr(9)); continue; }
    if (starts("SOL_COUNT:")) { cur.solCount = std::atoi(line.c_str() + 10); continue; }
    if (line == "SUBJECTS" || line == "SUBJECTS_OPEN") { sec = Section::Subjects; continue; }
    if (line == "CLIPS" || line == "CLIPS_OPEN") { sec = Section::Clips; continue; }
    auto poly = ParseClipper2PolyLine(line);
    if (poly.size() >= 3) {
      if (sec == Section::Clips) cur.clips.push_back(std::move(poly));
      else cur.subjects.push_back(std::move(poly));
    }
  }
  flush();
  return cases;
}

enum class Clipper2Reach { Supported, Unsupported };

inline Clipper2Reach ClassifyCase(const Clipper2Case& c) {
  if (c.fillrule == "NEGATIVE" || c.fillrule == "POSITIVE")
    return Clipper2Reach::Unsupported;
  if (c.fillrule != "NONZERO" && c.fillrule != "EVENODD")
    return Clipper2Reach::Unsupported;
  if (c.cliptype != "UNION" && c.cliptype != "INTERSECTION" &&
      c.cliptype != "DIFFERENCE")
    return Clipper2Reach::Unsupported;
  return Clipper2Reach::Supported;
}

// Pre-fill polys under the declared fill rule. Output is a CCW-canonical
// arrangement (no self-overlap). Reuses RemoveOverlaps2D with a custom
// predicate so we don't perturb the public API surface.
inline manifold::Polygons FillUnderRule(const manifold::Polygons& polys,
                                        const std::string& rule, double eps) {
  if (polys.empty()) return {};
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  for (const auto& loop : polys) {
    if (loop.size() < 3) continue;
    const int base = static_cast<int>(verts.size());
    const int n = static_cast<int>(loop.size());
    for (const auto& v : loop) verts.push_back(v);
    for (int i = 0; i < n; ++i)
      edges.push_back({base + i, base + ((i + 1) % n), 1});
  }
  if (verts.empty()) return {};
  WindPredicate pred = (rule == "EVENODD")
                           ? WindEvenOdd()
                           : WindPredicate([](int w) { return w != 0; });
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, pred);
  return OutEdgesToPolygons(r.verts, r.edges);
}

// Run one case under the universal pre-fill model. Returns the output
// polygons; sets `reach` to indicate whether the case was dispatched.
inline manifold::Polygons RunClipper2Case(const Clipper2Case& c,
                                          Clipper2Reach* reach,
                                          double* epsOut = nullptr) {
  *reach = ClassifyCase(c);
  if (*reach == Clipper2Reach::Unsupported) return {};
  const double eps = InferEps(c.subjects, c.clips);
  if (epsOut) *epsOut = eps;
  auto a = FillUnderRule(c.subjects, c.fillrule, eps);
  auto b = FillUnderRule(c.clips, c.fillrule, eps);
  if (c.cliptype == "UNION") {
    if (b.empty()) return a;
    if (a.empty()) return b;
    return Boolean2D(a, b, OpType::Add, eps);
  }
  if (c.cliptype == "INTERSECTION") {
    if (a.empty() || b.empty()) return {};
    return Boolean2D(a, b, OpType::Intersect, eps);
  }
  if (c.cliptype == "DIFFERENCE") {
    if (a.empty()) return {};
    if (b.empty()) return a;
    return Boolean2D(a, b, OpType::Subtract, eps);
  }
  return {};
}

inline void RunClipper2Corpus(const std::string& path) {
  auto cases = LoadClipper2Cases(path);
  std::cout << "=== Clipper2 corpus: " << cases.size() << " cases from "
            << path << " ===\n\n";
  // Per-case TSV (case_n, op, fillrule, structure, sol_area, output_area)
  // is written to /tmp/clipper2_proto_areas.tsv so the per-case area can
  // be diffed against Clipper2's own current output. SOL_AREA in the
  // corpus is older Clipper2 stored output and may differ from what
  // current Clipper2 produces for the same input.
  std::ofstream tsv("/tmp/clipper2_proto_areas.tsv");
  tsv << "n\top\tfillrule\tstructure\tsol_area\tproto_area\n";
  int total = 0;
  int run = 0, unsupported = 0;
  int topoOk = 0, topoFail = 0;
  int areaPreserved = 0;     // < 1% of SOL_AREA
  int areaDrifted = 0;       // 1% .. 10%
  int areaCollapsed = 0;     // > 10% (or zero output where SOL_AREA > 0)
  int solSkipped = 0;        // SOL_AREA <= 0 (no oracle)
  int outZeroOk = 0;         // SOL_AREA == 0 and we returned ~0 area
  std::vector<std::tuple<int, std::string, double, double, double>> drifters;
  // (case_n, descriptor, drift_pct, outArea, solArea)
  for (const auto& c : cases) {
    ++total;
    Clipper2Reach reach;
    double eps = 0.0;
    manifold::Polygons out;
    try {
      out = RunClipper2Case(c, &reach, &eps);
    } catch (const std::exception& ex) {
      std::cerr << "  case " << c.n << " (" << c.cliptype << "+" << c.fillrule
                << "): EXCEPTION " << ex.what() << "\n";
      ++topoFail;
      continue;
    }
    std::string desc = c.cliptype + "+" + c.fillrule +
                       (c.clips.empty() ? "/sub-only" : "/sub+clip");
    if (reach == Clipper2Reach::Unsupported) {
      ++unsupported;
      continue;
    }
    ++run;
    {
      const double outArea = std::fabs(TotalSignedArea(out));
      tsv << c.n << '\t' << c.cliptype << '\t' << c.fillrule << '\t'
          << (c.clips.empty() ? "subonly" : "subclip") << '\t'
          << c.solArea << '\t' << outArea << '\n';
    }

    // Topology check on the output polygons (treat as a fresh input).
    auto [v, e] = PolygonsToInput(out);
    if (!v.empty()) {
      auto r2 = RemoveOverlaps2D(v, e, eps);
      bool topo = CheckTopologicalValidity(r2, e, r2.inputRemap,
                                           r2.numMergedVerts);
      // For a valid output, re-running RemoveOverlaps2D should be a no-op
      // (idempotent under sub-eps drift). Topology failure here flags
      // either an invalid output or a corner case the pipeline mishandles
      // when fed back its own result.
      if (topo) ++topoOk;
      else ++topoFail;
    } else {
      // Empty output: only "valid" when SOL_AREA == 0 or unsupported.
      ++topoOk;
    }

    const double outArea = std::fabs(TotalSignedArea(out));
    if (c.solArea < 0) {
      ++solSkipped;
    } else if (c.solArea == 0.0) {
      if (outArea < 1.0) ++outZeroOk;
      else { ++areaCollapsed; drifters.emplace_back(c.n, desc, outArea, outArea, c.solArea); }
    } else {
      const double drift = std::fabs(outArea - c.solArea) / c.solArea;
      if (drift < 0.01) ++areaPreserved;
      else if (drift < 0.10) {
        ++areaDrifted;
        drifters.emplace_back(c.n, desc, drift, outArea, c.solArea);
      } else {
        ++areaCollapsed;
        drifters.emplace_back(c.n, desc, drift, outArea, c.solArea);
      }
    }
  }
  std::cout << "  Total cases:     " << total << "\n";
  std::cout << "  Run:             " << run << "\n";
  std::cout << "  Unsupported:     " << unsupported
            << " (NEGATIVE/POSITIVE fill rule)\n";
  std::cout << "\n  Of " << run << " cases run:\n";
  std::cout << "    Topology valid:   " << topoOk << " / " << run << "\n";
  std::cout << "    Topology INVALID: " << topoFail << " / " << run << "\n";
  std::cout << "    Area within  1% of SOL_AREA: " << areaPreserved << "\n";
  std::cout << "    Area within  1-10% of SOL_AREA: " << areaDrifted << "\n";
  std::cout << "    Area > 10% from SOL_AREA:     " << areaCollapsed << "\n";
  std::cout << "    SOL_AREA == 0, returned ~0:    " << outZeroOk << "\n";
  std::cout << "    SOL_AREA < 0 (no oracle):      " << solSkipped << "\n";
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 30):\n";
    for (size_t i = 0; i < drifters.size() && i < 30; ++i) {
      auto& [n, desc, drift, outA, solA] = drifters[i];
      std::cout << "    #" << n << "  " << desc
                << "  drift=" << (drift * 100.0) << "%"
                << "  out=" << outA << "  sol=" << solA << "\n";
    }
  }
}

}  // namespace overlap2d

// =============================================================================
// mfogel/polygon-clipping end-to-end fixtures.
//
// Independent third-party corpus harvested from real bug reports against
// the JS polygon-clipping library. Each fixture lives in its own dir
// under build/_deps/mfogel-polygon-clipping-src/test/end-to-end/<name>/:
//
//   args.geojson   — GeoJSON FeatureCollection of N input operands.
//                    Each Feature is a Polygon or MultiPolygon.
//   union.geojson, intersection.geojson, difference.geojson,
//   xor.geojson    — expected outputs for each op (at most one per file).
//   all.geojson    — expected output that's the same for every op
//                    (typical for single-input self-cleanup cases).
//   broken-issue-* — args only; no expected output. Skip these.
//
// GeoJSON convention (RFC 7946): exterior rings CCW, interior (hole)
// rings CW. NONZERO fill rule under that orientation isolates outer−hole.
// The runner pre-fills each feature under NONZERO, then chains the
// boolean op left-to-right across features.
//
// The corpus is small (~120 dirs) but its naming is itself a degeneracy
// taxonomy: triple-coincident-segments, almost-colinear-segments,
// infinitely-thin-polygon, self-intersects-but-doesnt-cross, ...
// This is exactly the layer Clipper2's int-coord corpus doesn't cover.
// =============================================================================
namespace overlap2d {

// ---------------------------------------------------------------------------
// Minimal JSON value + parser. Hand-rolled to avoid a third-party
// dependency for what is otherwise a single-translation-unit prototype.
// Supports only what GeoJSON needs: object, array, string, number, bool,
// null. No escape parsing in strings beyond \" and \\ — none of the
// fixture files use richer escapes.
// ---------------------------------------------------------------------------
struct JsonValue {
  enum Kind { Null, Bool, Number, String, Array, Object };
  Kind kind = Null;
  bool b = false;
  double num = 0.0;
  std::string str;
  std::vector<JsonValue> arr;
  std::map<std::string, JsonValue> obj;
};

struct JsonParser {
  const std::string& s;
  size_t i = 0;
  explicit JsonParser(const std::string& src) : s(src) {}
  void skipWs() {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  }
  bool match(char c) {
    skipWs();
    if (i < s.size() && s[i] == c) { ++i; return true; }
    return false;
  }
  void expect(char c) {
    if (!match(c)) throw std::runtime_error(
        std::string("expected '") + c + "' at " + std::to_string(i));
  }
  JsonValue parseValue() {
    skipWs();
    if (i >= s.size()) throw std::runtime_error("unexpected eof");
    char c = s[i];
    if (c == '{') return parseObject();
    if (c == '[') return parseArray();
    if (c == '"') return parseString();
    if (c == 't' || c == 'f') return parseBool();
    if (c == 'n') return parseNull();
    return parseNumber();
  }
  JsonValue parseObject() {
    JsonValue v; v.kind = JsonValue::Object;
    expect('{');
    skipWs();
    if (match('}')) return v;
    for (;;) {
      JsonValue key = parseString();
      expect(':');
      v.obj[key.str] = parseValue();
      skipWs();
      if (match(',')) continue;
      expect('}');
      break;
    }
    return v;
  }
  JsonValue parseArray() {
    JsonValue v; v.kind = JsonValue::Array;
    expect('[');
    skipWs();
    if (match(']')) return v;
    for (;;) {
      v.arr.push_back(parseValue());
      skipWs();
      if (match(',')) continue;
      expect(']');
      break;
    }
    return v;
  }
  JsonValue parseString() {
    JsonValue v; v.kind = JsonValue::String;
    expect('"');
    while (i < s.size() && s[i] != '"') {
      if (s[i] == '\\' && i + 1 < s.size()) {
        char n = s[i + 1];
        if (n == '"' || n == '\\') { v.str.push_back(n); i += 2; continue; }
        v.str.push_back(s[i]); ++i;
      } else {
        v.str.push_back(s[i]); ++i;
      }
    }
    expect('"');
    return v;
  }
  JsonValue parseBool() {
    JsonValue v; v.kind = JsonValue::Bool;
    if (s.compare(i, 4, "true") == 0) { v.b = true; i += 4; }
    else if (s.compare(i, 5, "false") == 0) { v.b = false; i += 5; }
    else throw std::runtime_error("bad bool");
    return v;
  }
  JsonValue parseNull() {
    JsonValue v; v.kind = JsonValue::Null;
    if (s.compare(i, 4, "null") == 0) { i += 4; return v; }
    throw std::runtime_error("bad null");
  }
  JsonValue parseNumber() {
    JsonValue v; v.kind = JsonValue::Number;
    size_t start = i;
    if (s[i] == '-' || s[i] == '+') ++i;
    while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) ||
                            s[i] == '.' || s[i] == 'e' || s[i] == 'E' ||
                            s[i] == '-' || s[i] == '+'))
      ++i;
    v.num = std::stod(s.substr(start, i - start));
    return v;
  }
};

inline JsonValue ParseJsonFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("open: " + path);
  std::ostringstream oss; oss << in.rdbuf();
  std::string text = oss.str();
  JsonParser p(text);
  return p.parseValue();
}

// mfogel's polygon-clipping library normalizes ring winding by position:
// the first ring of each polygon is outer (CCW), subsequent rings are
// holes (CW). It does this regardless of the input's actual winding —
// the corpus has explicit "backward-ring-winding-order" test cases that
// rely on this auto-normalization. Strict GeoJSON conformance (RFC 7946)
// agrees with this convention but doesn't mandate enforcement; mfogel
// enforces. We do the same here so NONZERO pre-fill on the combined
// rings reproduces outer−holes correctly.
inline void NormalizeMfogelWinding(manifold::Polygons& rings) {
  for (size_t i = 0; i < rings.size(); ++i) {
    const double a = SignedArea(rings[i]);
    const bool wantCCW = (i == 0);
    if ((wantCCW && a < 0) || (!wantCCW && a > 0))
      std::reverse(rings[i].begin(), rings[i].end());
  }
}

// Extract a list of polygons from one Feature's geometry.
// Polygon: coordinates = [ring1, ring2, ...] → one polygon with multiple
//          rings (exterior + holes).
// MultiPolygon: coordinates = [[ring1, ring2], [ring1], ...] → one
//          polygon per outer entry.
// Returns vector<Polygons> (each entry = one polygon, possibly with holes).
inline std::vector<manifold::Polygons> ExtractFeatureGeometry(
    const JsonValue& feature) {
  std::vector<manifold::Polygons> polys;
  auto gIt = feature.obj.find("geometry");
  if (gIt == feature.obj.end() || gIt->second.kind != JsonValue::Object)
    return polys;
  const auto& geom = gIt->second;
  auto tIt = geom.obj.find("type");
  auto cIt = geom.obj.find("coordinates");
  if (tIt == geom.obj.end() || cIt == geom.obj.end()) return polys;
  const std::string& gtype = tIt->second.str;
  auto parseRing = [](const JsonValue& ring) {
    manifold::SimplePolygon out;
    for (const auto& pt : ring.arr) {
      if (pt.arr.size() >= 2)
        out.push_back({pt.arr[0].num, pt.arr[1].num});
    }
    // Drop the closing duplicate vertex if present (RFC 7946 requires it).
    if (out.size() >= 2 && out.front().x == out.back().x &&
        out.front().y == out.back().y)
      out.pop_back();
    return out;
  };
  // mfogel drops the entire polygon if its outer (first) ring is
  // degenerate (< 3 verts after closure-removal). The corpus has the
  // explicit `rings-with-no-area` case for this. Subsequent (hole)
  // rings can be silently dropped if degenerate; the polygon survives.
  auto collect = [&](const JsonValue& ringList,
                     manifold::Polygons& dest) -> bool {
    bool first = true;
    for (const auto& ring : ringList.arr) {
      auto r = parseRing(ring);
      if (first) {
        first = false;
        if (r.size() < 3) return false;
        dest.push_back(std::move(r));
      } else {
        if (r.size() >= 3) dest.push_back(std::move(r));
      }
    }
    return !dest.empty();
  };
  if (gtype == "Polygon") {
    manifold::Polygons one;
    if (collect(cIt->second, one)) {
      NormalizeMfogelWinding(one);
      polys.push_back(std::move(one));
    }
  } else if (gtype == "MultiPolygon") {
    for (const auto& poly : cIt->second.arr) {
      manifold::Polygons one;
      if (collect(poly, one)) {
        NormalizeMfogelWinding(one);
        polys.push_back(std::move(one));
      }
    }
  }
  return polys;
}

// Flatten a list-of-Polygons (one per outer polygon in a MultiPolygon)
// into a single Polygons (concat all rings). NONZERO pre-fill on the
// concatenation reproduces the multipolygon's covered region under
// GeoJSON's CCW-outer / CW-hole convention.
inline manifold::Polygons FlattenMulti(
    const std::vector<manifold::Polygons>& multi) {
  manifold::Polygons flat;
  for (const auto& p : multi)
    for (const auto& r : p) flat.push_back(r);
  return flat;
}

// Top-level extractor: returns one polygon-group per Feature in the file.
// Handles either a single Feature or a FeatureCollection at the root.
inline std::vector<std::vector<manifold::Polygons>> ExtractFeatures(
    const std::string& path) {
  JsonValue root = ParseJsonFile(path);
  std::vector<std::vector<manifold::Polygons>> out;
  if (root.kind != JsonValue::Object) return out;
  auto tIt = root.obj.find("type");
  if (tIt == root.obj.end()) return out;
  if (tIt->second.str == "FeatureCollection") {
    auto fIt = root.obj.find("features");
    if (fIt == root.obj.end()) return out;
    for (const auto& f : fIt->second.arr)
      out.push_back(ExtractFeatureGeometry(f));
  } else if (tIt->second.str == "Feature") {
    out.push_back(ExtractFeatureGeometry(root));
  }
  return out;
}

struct MfogelCase {
  std::string name;
  std::vector<manifold::Polygons> features;  // one entry per arg feature
  std::map<std::string, manifold::Polygons> expected;  // op -> expected (flat)
};

inline std::vector<MfogelCase> LoadMfogelCorpus(const std::string& dir) {
  std::vector<MfogelCase> cases;
  // Walk subdirs without <filesystem> dependency: ls via popen.
  std::string cmd = "ls -1 " + dir;
  FILE* fp = popen(cmd.c_str(), "r");
  if (!fp) return cases;
  char buf[1024];
  std::vector<std::string> names;
  while (std::fgets(buf, sizeof(buf), fp)) {
    std::string n = buf;
    if (!n.empty() && n.back() == '\n') n.pop_back();
    if (!n.empty()) names.push_back(n);
  }
  pclose(fp);
  for (const auto& name : names) {
    const std::string caseDir = dir + "/" + name;
    const std::string args = caseDir + "/args.geojson";
    std::ifstream chk(args);
    if (!chk) continue;
    MfogelCase c;
    c.name = name;
    try {
      auto feats = ExtractFeatures(args);
      for (auto& f : feats) {
        // Collapse each feature's polygon list into one rings-list. We
        // pre-fill under NONZERO so the orientation-encoded holes work.
        c.features.push_back(FlattenMulti(f));
      }
    } catch (const std::exception& e) {
      std::cerr << "  parse failed for " << name << ": " << e.what() << "\n";
      continue;
    }
    for (const char* op :
         {"union", "intersection", "difference", "xor", "all"}) {
      const std::string p = caseDir + "/" + op + ".geojson";
      std::ifstream chk2(p);
      if (!chk2) continue;
      try {
        auto feats = ExtractFeatures(p);
        manifold::Polygons exp;
        for (auto& f : feats) {
          auto flat = FlattenMulti(f);
          for (auto& r : flat) exp.push_back(std::move(r));
        }
        c.expected[op] = std::move(exp);
      } catch (...) {}
    }
    cases.push_back(std::move(c));
  }
  return cases;
}

// Run one (case, op) pair through the prototype as a single-pass
// combined-input pipeline.
//
// Earlier this chained `Boolean2D(acc, filled[i], op)` left-to-right.
// The chain breaks on N>=3-feature cases like `issue-44` because the
// intermediate `acc` accumulates T-junction vertices (extra colinear
// verts on shared edges) that the next feature doesn't have, and the
// pipeline mishandles those collinear runs at the next merge.
//
// Single-pass dispatch sidesteps the chain by feeding all features into
// one RemoveOverlaps2D call with appropriate multiplicities + predicate:
//
//   union          → mult=+1 each, predicate w != 0
//   xor            → mult=+1 each, predicate (w & 1)
//   intersection   → mult=+1 each, predicate w >= N (covered by every input)
//   difference     → A mult=+1, others mult=-1, predicate w > 0
//
// Each feature is pre-filled under NONZERO first so its rings are
// CCW-canonical before being merged into the combined input.
inline manifold::Polygons RunMfogelOp(const MfogelCase& c,
                                      const std::string& op, double eps) {
  if (c.features.empty()) return {};
  std::vector<manifold::Polygons> filled;
  for (const auto& f : c.features) {
    if (f.empty()) { filled.push_back({}); continue; }
    filled.push_back(FillUnderRule(f, "NONZERO", eps));
  }
  if (op == "all") return filled[0];

  // Combine all features' edges into one input list with appropriate
  // per-feature multiplicities.
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  size_t nNonEmpty = 0;
  auto append = [&](const manifold::Polygons& polys, int mult) {
    if (polys.empty()) return;
    ++nNonEmpty;
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i)
        edges.push_back({base + i, base + ((i + 1) % n), mult});
    }
  };
  if (op == "difference") {
    append(filled[0], +1);
    for (size_t i = 1; i < filled.size(); ++i) append(filled[i], -1);
  } else {
    for (auto& f : filled) append(f, +1);
  }
  if (verts.empty()) return {};
  WindPredicate pred;
  if (op == "union") {
    pred = [](int w) { return w != 0; };
  } else if (op == "xor") {
    pred = WindEvenOdd();
  } else if (op == "intersection") {
    const int N = static_cast<int>(filled.size());
    pred = [N](int w) { return w >= N; };
  } else if (op == "difference") {
    pred = WindAdd();  // w > 0
  } else {
    return {};
  }
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, pred);
  return OutEdgesToPolygons(r.verts, r.edges);
}

inline double InferEpsForFeatures(
    const std::vector<manifold::Polygons>& features) {
  manifold::Polygons all;
  for (const auto& f : features)
    for (const auto& r : f) all.push_back(r);
  return InferEps(all, {});
}

inline void RunMfogelCorpus(const std::string& dir) {
  auto cases = LoadMfogelCorpus(dir);
  std::cout << "=== mfogel/polygon-clipping corpus: " << cases.size()
            << " dirs from " << dir << " ===\n\n";
  int casesWithExpected = 0;
  int totalOps = 0;
  int topoOk = 0, topoFail = 0;
  int areaPreserved = 0, areaDrifted = 0, areaCollapsed = 0;
  int oracleEmptyOk = 0;          // expected = 0 area, we returned ~0
  int oracleEmptyFail = 0;        // expected = 0, we returned non-zero
  std::vector<std::tuple<std::string, std::string, double, double, double>>
      drifters;
  for (const auto& c : cases) {
    if (c.expected.empty()) continue;
    ++casesWithExpected;
    const double eps = InferEpsForFeatures(c.features);
    for (const auto& [op, expected] : c.expected) {
      ++totalOps;
      manifold::Polygons out;
      try {
        out = RunMfogelOp(c, op, eps);
      } catch (const std::exception& e) {
        std::cerr << "  EXC " << c.name << "/" << op << ": " << e.what()
                  << "\n";
        ++topoFail;
        continue;
      }
      // Topology check by re-feeding output into the pipeline.
      auto [v, e] = PolygonsToInput(out);
      if (!v.empty()) {
        auto r2 = RemoveOverlaps2D(v, e, eps);
        bool topo = CheckTopologicalValidity(r2, e, r2.inputRemap,
                                             r2.numMergedVerts);
        if (topo) ++topoOk;
        else ++topoFail;
      } else {
        ++topoOk;
      }
      const double outA = std::fabs(TotalSignedArea(out));
      const double expA = std::fabs(TotalSignedArea(expected));
      if (expA < 1e-9) {
        if (outA < 1e-9) ++oracleEmptyOk;
        else { ++oracleEmptyFail; drifters.emplace_back(c.name, op, outA, outA, 0.0); }
      } else {
        const double drift = std::fabs(outA - expA) / expA;
        if (drift < 0.01) ++areaPreserved;
        else if (drift < 0.10) {
          ++areaDrifted;
          drifters.emplace_back(c.name, op, drift, outA, expA);
        } else {
          ++areaCollapsed;
          drifters.emplace_back(c.name, op, drift, outA, expA);
        }
      }
    }
  }
  std::cout << "  Cases with expected-output files: " << casesWithExpected
            << " / " << cases.size() << "\n";
  std::cout << "  (Skipped: cases with only args.geojson — known broken in"
               " mfogel itself, e.g. broken-issue-*.)\n";
  std::cout << "  Total (case, op) pairs run:       " << totalOps << "\n\n";
  std::cout << "  Topology valid:   " << topoOk << " / " << totalOps << "\n";
  std::cout << "  Topology INVALID: " << topoFail << " / " << totalOps << "\n";
  std::cout << "  Area within  1%:  " << areaPreserved << "\n";
  std::cout << "  Area 1-10%:        " << areaDrifted << "\n";
  std::cout << "  Area >= 10%:       " << areaCollapsed << "\n";
  std::cout << "  Expected empty, returned ~0:  " << oracleEmptyOk << "\n";
  std::cout << "  Expected empty, returned >0: " << oracleEmptyFail << "\n";
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 30):\n";
    for (size_t i = 0; i < drifters.size() && i < 30; ++i) {
      auto& [name, op, drift, outA, expA] = drifters[i];
      std::cout << "    " << name << " / " << op
                << "  drift=" << (drift * 100.0) << "%"
                << "  out=" << outA << "  exp=" << expA << "\n";
    }
  }
}

}  // namespace overlap2d

// =============================================================================
// JTS robust/overlay corpus runner.
//
// JTS Topology Suite (LocationTech) ships an XML+WKB test format with
// cases harvested from real GIS bug reports (GEOS, PostGIS, QGIS, JTS
// itself, Shapely, OSGeo). The dominant op is `overlayAreaTest`, which
// checks the invariant
//
//     |A∪B| + |A∩B| = |A| + |B|
//
// This is a self-checking oracle — no expected-output WKB needed —
// and it's the strongest single cross-check available: a bug in either
// union or intersection that creates or loses area shows up immediately.
// Real-bug-derived inputs from production GIS systems mean the invariant
// is the right shape of test for them.
//
// XML+WKB parsing lives outside this prototype in test/polygons/
// jts_to_corpus.py, which produces test/polygons/jts_overlay_corpus.txt.
// That keeps the WKB/XML/EWKB/WKT plumbing out of the algorithm code
// and makes the corpus inspectable as plain text. The same converter
// handles other JTS-format-compatible corpora (GEOS testxml, etc.) if
// added later.
//
// Format consumed here:
//
//   CASE <n> <op> <source>[#<idx>] [<expected_value>]
//   A <num_rings>
//   <num_pts> x y x y ...    (one ring per line; coords space-separated)
//   B <num_rings>            (omit if no B)
//   <num_pts> ...
//   <blank line>
//
// Supported ops: `overlayAreaTest` (invariant) and `unionArea` (oracle).
// =============================================================================
namespace overlap2d {

struct JtsCase {
  int n = 0;
  std::string op;
  std::string source;
  double expected = 0.0;
  bool hasExpected = false;
  manifold::Polygons a;  // flat ring list (orientation preserved)
  manifold::Polygons b;
};

inline std::vector<JtsCase> LoadJtsCorpus(const std::string& path) {
  std::vector<JtsCase> cases;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return cases;
  }
  std::string line;
  JtsCase cur;
  bool inCase = false;
  // When we see "A <n>" or "B <n>" we then expect <n> ring lines.
  int ringsRemaining = 0;
  bool readingA = false;
  auto flush = [&]() {
    if (inCase && cur.n > 0) cases.push_back(std::move(cur));
    cur = JtsCase{};
    inCase = false;
    ringsRemaining = 0;
    readingA = false;
  };
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
      line.pop_back();
    if (line.empty() || line[0] == '#') {
      if (line.empty()) flush();
      continue;
    }
    std::istringstream ss(line);
    if (ringsRemaining > 0) {
      // Ring line: <num_pts> x y x y ...
      int npts;
      ss >> npts;
      manifold::SimplePolygon ring;
      ring.reserve(npts);
      for (int i = 0; i < npts; ++i) {
        double x, y;
        if (!(ss >> x >> y)) break;
        ring.push_back({x, y});
      }
      if (ring.size() >= 3) {
        if (readingA) cur.a.push_back(std::move(ring));
        else cur.b.push_back(std::move(ring));
      }
      --ringsRemaining;
      continue;
    }
    std::string token;
    ss >> token;
    if (token == "CASE") {
      flush();
      inCase = true;
      ss >> cur.n >> cur.op >> cur.source;
      double v;
      if (ss >> v) { cur.expected = v; cur.hasExpected = true; }
    } else if (token == "A") {
      ss >> ringsRemaining;
      readingA = true;
    } else if (token == "B") {
      ss >> ringsRemaining;
      readingA = false;
    }
    // Other lines silently ignored (forward-compat with format extensions).
  }
  flush();
  return cases;
}

inline void RunJtsCorpus(const std::string& path) {
  auto cases = LoadJtsCorpus(path);
  std::cout << "=== JTS overlay corpus: " << cases.size()
            << " cases from " << path << " ===\n\n";
  size_t areaTests = 0, unionAreaTests = 0, skipped = 0;
  size_t topoOk = 0, topoFail = 0;
  size_t invariantHeld = 0;
  size_t invariantDriftSmall = 0;  // < 1%
  size_t invariantDriftLarge = 0;  // >= 1%
  size_t unionAreaOk = 0, unionAreaDrift = 0;
  std::vector<std::tuple<std::string, double, double, double, double>>
      drifters;  // (source, drift_pct, A∪B, A∩B, |A|+|B|)
  for (const auto& c : cases) {
    if (c.op != "overlayAreaTest" && c.op != "overlayareatest" &&
        c.op != "unionArea") {
      ++skipped;
      continue;
    }
    try {
      const double eps = InferEps(c.a, c.b);
      auto fa = FillUnderRule(c.a, "NONZERO", eps);
      if (c.op == "unionArea") {
        ++unionAreaTests;
        auto [v, e] = PolygonsToInput(fa);
        if (!v.empty()) {
          auto r2 = RemoveOverlaps2D(v, e, eps);
          if (CheckTopologicalValidity(r2, e, r2.inputRemap,
                                       r2.numMergedVerts)) ++topoOk;
          else ++topoFail;
        } else ++topoOk;
        const double a = std::fabs(TotalSignedArea(fa));
        if (c.hasExpected && c.expected > 0) {
          double drift = std::fabs(a - c.expected) / c.expected;
          if (drift < 0.01) ++unionAreaOk; else ++unionAreaDrift;
        } else {
          ++unionAreaOk;
        }
        continue;
      }
      ++areaTests;
      auto fb = FillUnderRule(c.b, "NONZERO", eps);
      const double Aa = std::fabs(TotalSignedArea(fa));
      const double Ab = std::fabs(TotalSignedArea(fb));
      manifold::Polygons unionAB, interAB;
      if (fa.empty()) unionAB = fb;
      else if (fb.empty()) unionAB = fa;
      else unionAB = Boolean2D(fa, fb, OpType::Add, eps);
      if (!fa.empty() && !fb.empty())
        interAB = Boolean2D(fa, fb, OpType::Intersect, eps);
      const double Au = std::fabs(TotalSignedArea(unionAB));
      const double Ai = std::fabs(TotalSignedArea(interAB));
      auto [v, e] = PolygonsToInput(unionAB);
      if (!v.empty()) {
        auto r2 = RemoveOverlaps2D(v, e, eps);
        if (CheckTopologicalValidity(r2, e, r2.inputRemap,
                                     r2.numMergedVerts)) ++topoOk;
        else ++topoFail;
      } else ++topoOk;
      const double lhs = Au + Ai;
      const double rhs = Aa + Ab;
      const double scale = std::max(rhs, 1e-12);
      const double drift = std::fabs(lhs - rhs) / scale;
      if (drift < 1e-9) ++invariantHeld;
      else if (drift < 0.01) ++invariantDriftSmall;
      else {
        ++invariantDriftLarge;
        drifters.emplace_back(c.source, drift, Au, Ai, rhs);
      }
    } catch (const std::exception& ex) {
      std::cerr << "  EXC " << c.source << ": " << ex.what() << "\n";
      ++topoFail;
    }
  }
  std::cout << "  overlayAreaTest run:   " << areaTests << "\n";
  std::cout << "  unionArea run:         " << unionAreaTests << "\n";
  std::cout << "  Skipped (other ops):   " << skipped << "\n\n";
  std::cout << "  Topology valid:   " << topoOk << "\n";
  std::cout << "  Topology INVALID: " << topoFail << "\n\n";
  std::cout << "  Invariant holds (drift < 1e-9):  " << invariantHeld
            << " / " << areaTests << "\n";
  std::cout << "  Invariant near-holds (< 1%):     " << invariantDriftSmall
            << " / " << areaTests << "\n";
  std::cout << "  Invariant violated (>= 1%):      " << invariantDriftLarge
            << " / " << areaTests << "\n";
  if (unionAreaTests > 0) {
    std::cout << "  unionArea matches expected:   " << unionAreaOk << " / "
              << unionAreaTests << "\n";
    std::cout << "  unionArea drifts:             " << unionAreaDrift << " / "
              << unionAreaTests << "\n";
  }
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<1>(a) > std::get<1>(b);
              });
    std::cout << "\n  Invariant violators (worst first, up to 20):\n";
    for (size_t i = 0; i < drifters.size() && i < 20; ++i) {
      auto& [src, drift, Au, Ai, rhs] = drifters[i];
      std::cout << "    " << src
                << "  drift=" << (drift * 100.0) << "%"
                << "  A∪B=" << Au << "  A∩B=" << Ai << "  |A|+|B|=" << rhs
                << "\n";
    }
  }
}

}  // namespace overlap2d

// =============================================================================
// Head-to-head benchmark vs Clipper2 (the library overlap2d is intended to
// eventually replace as `CrossSection`'s boolean/Simplify backend). Built
// only when -DOVERLAP2D_WITH_CLIPPER2 is defined and Clipper2 is linked.
//
// Compile/link:
//   g++ ... -DOVERLAP2D_WITH_CLIPPER2 -DMANIFOLD_PAR=1 \
//     -I build/_deps/clipper2-src/CPP/Clipper2Lib/include \
//     extras/overlap2d_proto.cpp \
//     build/_deps/clipper2-build/libClipper2.a -ltbb \
//     -o overlap2d_proto
//
// What's timed: the boolean call itself, on already-loaded inputs. We
// exclude geometry conversion (manifold::Polygons → Clipper2 PathsD)
// because that's a one-shot cost both sides would amortize out under
// any real-world wiring.
//
// Workloads: each of the three corpora we already validate against,
// run a UNION on every (sub, clip) pair (or every subject-only case).
// Clipper2's NONZERO is the cross-comparable fill rule; cases with
// EVENODD or NEGATIVE/POSITIVE are skipped here so both libraries do
// the same work. (Correctness comparison happens in the corpus runs;
// this benchmark is purely about wall-clock at parity.)
//
// Reference: `Clipper64` (native int64 API), not `ClipperD` (real-coord
// wrapper). Per elalish/manifold#1683 the planned post-precision-bug
// `CrossSection` wiring is `Paths64`-based, so `Clipper64` is the right
// baseline; `ClipperD`'s per-op rescaling would flatter our number
// against an end state nobody plans to ship.
// =============================================================================
#ifdef OVERLAP2D_WITH_CLIPPER2
#include "clipper2/clipper.h"
namespace overlap2d {

inline Clipper2Lib::PathsD ToPathsD(const manifold::Polygons& p) {
  Clipper2Lib::PathsD out;
  out.reserve(p.size());
  for (const auto& loop : p) {
    Clipper2Lib::PathD path;
    path.reserve(loop.size());
    for (const auto& v : loop) path.emplace_back(v.x, v.y);
    out.push_back(std::move(path));
  }
  return out;
}

// Scale doubles to int64 with a fixed factor and convert to Paths64 (the
// native Clipper2 API; ClipperD wraps this with rescaling at every op).
// Per elalish/manifold#1683 the planned post-Clipper2-precision-bug
// CrossSection wiring is Paths64-based, so the head-to-head we want
// reports against this faster, full-int-precision path. Scale chosen
// so a 1e6-coordinate input still has ~1 unit precision in int64.
constexpr double kPaths64Scale = 1e9;
inline Clipper2Lib::Paths64 ToPaths64(const manifold::Polygons& p) {
  Clipper2Lib::Paths64 out;
  out.reserve(p.size());
  for (const auto& loop : p) {
    Clipper2Lib::Path64 path;
    path.reserve(loop.size());
    for (const auto& v : loop) {
      path.emplace_back(static_cast<int64_t>(v.x * kPaths64Scale),
                        static_cast<int64_t>(v.y * kPaths64Scale));
    }
    out.push_back(std::move(path));
  }
  return out;
}

struct BenchTotals {
  int64_t oursNs = 0;
  int64_t clipperNs = 0;
  int cases = 0;
  int oursDrops = 0;     // we returned empty
  int clipperDrops = 0;  // clipper2 returned empty
};

inline void TimeOneUnion(const manifold::Polygons& a,
                         const manifold::Polygons& b, BenchTotals* totals) {
  using Clock = std::chrono::steady_clock;
  const double eps = InferEps(a, b);
  // Ours: single-pass dispatch (mfogel-runner shape). Combine subject
  // and clip edges with mult=+1 each, run RemoveOverlaps2D once with
  // a NONZERO `w != 0` predicate — matches Clipper2's UNION+NONZERO
  // semantics directly without paying for two extra FillUnderRule
  // pipeline runs that the multi-pass Boolean2D path would incur.
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  auto append = [&](const manifold::Polygons& polys) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i)
        edges.push_back({base + i, base + ((i + 1) % n), 1});
    }
  };
  manifold::Polygons ours;
  auto t0 = Clock::now();
  append(a);
  append(b);
  if (!verts.empty()) {
    auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                              [](int w) { return w != 0; });
    ours = OutEdgesToPolygons(r.verts, r.edges);
  }
  auto t1 = Clock::now();
  totals->oursNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        t1 - t0).count();
  if (ours.empty()) ++totals->oursDrops;

  // Clipper2 Clipper64 (native int64 API; what CrossSection is targeting
  // post-#1683 to avoid PathsD's precision cap and rescale overhead).
  // Conversion to Paths64 is excluded from the timed region; we only
  // time Execute, the same as for overlap2d.
  auto a2 = ToPaths64(a);
  auto b2 = ToPaths64(b);
  Clipper2Lib::Clipper64 clip;
  clip.AddSubject(a2);
  if (!b.empty()) clip.AddClip(b2);
  Clipper2Lib::Paths64 sol;
  auto t2 = Clock::now();
  clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero,
               sol);
  auto t3 = Clock::now();
  totals->clipperNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           t3 - t2).count();
  if (sol.empty()) ++totals->clipperDrops;
  ++totals->cases;
}

inline void ReportBench(const std::string& label, const BenchTotals& t) {
  auto fmt = [](int64_t ns) {
    std::ostringstream o; o.setf(std::ios::fixed); o.precision(3);
    o << ns * 1e-6 << " ms";
    return o.str();
  };
  std::cout << "  " << label << ":\n";
  std::cout << "    cases:        " << t.cases << "\n";
  std::cout << "    overlap2d:    " << fmt(t.oursNs)
            << " (" << (t.cases ? t.oursNs / t.cases / 1000.0 : 0)
            << " µs/case avg)";
  if (t.oursDrops) std::cout << "  drops=" << t.oursDrops;
  std::cout << "\n";
  std::cout << "    Clipper2 :    " << fmt(t.clipperNs)
            << " (" << (t.cases ? t.clipperNs / t.cases / 1000.0 : 0)
            << " µs/case avg)";
  if (t.clipperDrops) std::cout << "  drops=" << t.clipperDrops;
  std::cout << "\n";
  if (t.clipperNs > 0) {
    const double ratio = static_cast<double>(t.oursNs) / t.clipperNs;
    std::cout << "    ratio:        ";
    std::cout.setf(std::ios::fixed); std::cout.precision(2);
    std::cout << ratio << "x  ("
              << (ratio < 1 ? "overlap2d faster" : "Clipper2 faster")
              << ")\n";
    std::cout.unsetf(std::ios::fixed);
  }
}

inline void RunVsClipper2_Clipper2Corpus(const std::string& path) {
  auto cases = LoadClipper2Cases(path);
  BenchTotals t;
  for (const auto& c : cases) {
    if (c.fillrule != "NONZERO" || c.cliptype != "UNION") continue;
    TimeOneUnion(c.subjects, c.clips, &t);
  }
  ReportBench("Clipper2 corpus / UNION+NONZERO", t);
}

inline void RunVsClipper2_JtsCorpus(const std::string& path) {
  auto cases = LoadJtsCorpus(path);
  BenchTotals t;
  for (const auto& c : cases) {
    // overlayAreaTest cases run union(A, B); unionArea is single-input.
    if (c.op == "overlayAreaTest" || c.op == "overlayareatest") {
      TimeOneUnion(c.a, c.b, &t);
    } else if (c.op == "unionArea") {
      TimeOneUnion(c.a, {}, &t);
    }
  }
  ReportBench("JTS overlay / UNION", t);
}

inline void RunVsClipper2_MfogelCorpus(const std::string& dir) {
  auto cases = LoadMfogelCorpus(dir);
  BenchTotals t;
  for (const auto& c : cases) {
    // Run a UNION on the first two features (the most common 2-input
    // shape; subject-only single-feature cases are skipped here since
    // that's a self-clean and Clipper2's fill-rule handling differs).
    if (c.features.size() >= 2 && !c.features[0].empty() &&
        !c.features[1].empty()) {
      TimeOneUnion(c.features[0], c.features[1], &t);
    }
  }
  ReportBench("mfogel end-to-end / UNION", t);
}

inline void RunVsClipper2(const std::string& which) {
  std::cout << "=== overlap2d vs Clipper2 head-to-head (UNION+NONZERO) ===\n\n";
  if (which == "all" || which == "clipper2") {
    RunVsClipper2_Clipper2Corpus(
        "build/_deps/clipper2-src/Tests/Polygons.txt");
  }
  if (which == "all" || which == "jts") {
    RunVsClipper2_JtsCorpus("test/polygons/jts_overlay_corpus.txt");
  }
  if (which == "all" || which == "mfogel") {
    RunVsClipper2_MfogelCorpus(
        "build/_deps/mfogel-polygon-clipping-src/test/end-to-end");
  }
}

// Diagnostic: walk the JTS corpus and print every case where overlap2d's
// output diverges from Clipper2's (one drops while the other doesn't, or
// area differs by more than a small tolerance). Used to triage drop cases
// without polluting the timed benchmark path.
inline void RunJtsDropDiag(const std::string& path) {
  auto cases = LoadJtsCorpus(path);
  std::cout << "=== JTS drop / divergence diagnostic ===\n";
  int oursDrops = 0, clipperDrops = 0, bothEmpty = 0, bothNonEmpty = 0;
  int ourEmptyTheyNot = 0, theirEmptyOursNot = 0;
  for (const auto& c : cases) {
    bool isUnion = (c.op == "overlayAreaTest" || c.op == "overlayareatest");
    bool isUnionArea = (c.op == "unionArea");
    if (!isUnion && !isUnionArea) continue;
    const manifold::Polygons& A = c.a;
    const manifold::Polygons& B = isUnion ? c.b : manifold::Polygons{};
    const double eps = InferEps(A, B);
    // Build inputs for overlap2d.
    std::vector<vec2> verts;
    std::vector<EdgeM> edges;
    auto append = [&](const manifold::Polygons& polys) {
      for (const auto& loop : polys) {
        if (loop.size() < 3) continue;
        const int base = static_cast<int>(verts.size());
        const int n = static_cast<int>(loop.size());
        for (const auto& v : loop) verts.push_back(v);
        for (int i = 0; i < n; ++i)
          edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
    };
    append(A); append(B);
    manifold::Polygons ours;
    if (!verts.empty()) {
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                                [](int w) { return w != 0; });
      ours = OutEdgesToPolygons(r.verts, r.edges);
    }
    // Clipper2 reference.
    auto a2 = ToPaths64(A);
    auto b2 = ToPaths64(B);
    Clipper2Lib::Clipper64 clip;
    clip.AddSubject(a2);
    if (!B.empty()) clip.AddClip(b2);
    Clipper2Lib::Paths64 sol;
    clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero,
                 sol);
    const bool oursEmpty = ours.empty();
    const bool theirsEmpty = sol.empty();
    if (oursEmpty) ++oursDrops;
    if (theirsEmpty) ++clipperDrops;
    if (oursEmpty && theirsEmpty) ++bothEmpty;
    if (!oursEmpty && !theirsEmpty) ++bothNonEmpty;
    if (oursEmpty && !theirsEmpty) {
      ++ourEmptyTheyNot;
      // Clipper2 area for context.
      double clipArea = 0;
      for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      // Total input vert/edge count for size context.
      int totalVerts = 0, totalRings = 0;
      for (const auto& p : A) { totalVerts += p.size(); ++totalRings; }
      for (const auto& p : B) { totalVerts += p.size(); ++totalRings; }
      std::cout << "OURS-DROP case=" << c.n << " op=" << c.op
                << " src=" << c.source
                << " A.rings=" << A.size() << " B.rings=" << B.size()
                << " verts=" << totalVerts
                << " eps=" << eps
                << " clipperArea=" << clipArea << "\n";
    }
    if (!oursEmpty && theirsEmpty) {
      ++theirEmptyOursNot;
      std::cout << "THEY-DROP case=" << c.n << " op=" << c.op
                << " src=" << c.source << "\n";
    }
  }
  std::cout << "\nSummary:\n";
  std::cout << "  total cases:           " << cases.size() << "\n";
  std::cout << "  both nonempty:         " << bothNonEmpty << "\n";
  std::cout << "  both empty:            " << bothEmpty << "\n";
  std::cout << "  ours drops, theirs ok: " << ourEmptyTheyNot << "\n";
  std::cout << "  theirs drops, ours ok: " << theirEmptyOursNot << "\n";
  std::cout << "  total ours drops:      " << oursDrops << "\n";
  std::cout << "  total theirs drops:    " << clipperDrops << "\n";
}

}  // namespace overlap2d
#endif  // OVERLAP2D_WITH_CLIPPER2

// =============================================================================
// Deep fuzz: wider seed sweep + iter-vs-topology verification.
//
// The standard displacement fuzz uses 30 seeds per (kPow, n) cell. This
// mode runs `seedsPerCell` seeds per cell across a wider parameter grid
// and, for every case the iterate-to-fixed-point reports as iter ≥ 2,
// compares pass-1 vs pass-2 outputs at:
//
//   - eps/100 quantum (the standard fingerprint, what triggers iter ≥ 2)
//   - eps quantum (topology-only; captures sub-edges modulo sub-eps drift)
//
// If pass-1 and pass-2 differ at eps/100 but match at eps for every iter≥2
// case, the iteration is shifting positions only, not topology. That's the
// cut-corner answer for whether iter=2 is doing real work.
// =============================================================================
namespace overlap2d {
void DeepFuzz(int seedsPerCell) {
  const std::vector<int> kPows = {10, 20, 30, 35, 40, 45, 49};
  const std::vector<int> sizes = {8, 20, 50, 100};
  std::cout << "=== DeepFuzz: " << seedsPerCell << " seeds × " << sizes.size()
            << " sizes × " << kPows.size() << " kPow = "
            << (seedsPerCell * sizes.size() * kPows.size()) << " cases ===\n";

  int total = 0;
  int firstPassFail = 0;
  int convergedPassValid = 0;     // converged result: topo valid
  int convergedPassInvalid = 0;   // converged result: topo invalid
  // Polygons-API regularization counter. Each fuzz case is run through
  // both the lower-level pipeline (counts every retained sub-edge) and
  // the public Simplify (regularizes via OutEdgesToPolygons, dropping
  // zero-area lens loops). When the counts disagree, the difference is
  // exactly the regularization drop, which is correct CGAL/Clipper2/SVG
  // behavior. We track the rate as a regularization-frequency stat
  // rather than as a "mismatch" (it isn't a bug), and as a tripwire
  // for any future drift in regularization behavior.
  int regUnchanged = 0;       // post-regularization edge count == lower-level
  int regDropped = 0;         // dropped one or more zero-area lens loops
  std::map<int, int> iterDist;
  int iterGE2 = 0;
  int iterGE2_topo_match = 0;     // pass1 == pass2 at eps quantum
  int iterGE2_fine_match = 0;     // pass1 == pass2 at eps/100 quantum
  int iterGE2_topo_diff = 0;
  int iterGE2_pass1_invalid = 0;  // iter≥2 cases where pass 1 was wrong
  int iterGE2_pass1_valid = 0;    // iter≥2 cases where pass 1 was already valid
  int iter_repaired = 0;          // pass-1 invalid → converged valid
  int iter_degraded = 0;          // pass-1 valid → converged INVALID
  int pass1_invalid_unfixed = 0;  // pass-1 invalid AND converged invalid
  // Geometric-emptiness check: counts cases where pass 2 outputs zero edges
  // even though pass 1 output non-zero. Indicates pass-1 → pass-2 winding-
  // sign disagreement (output orientation not Smith-convention-compatible).
  int pass2_collapsed = 0;
  int pass1_nonempty = 0;
  // Area-drift tracking: a quantitative oracle on top of the topology
  // check. Pass 1 produces an output area; iterate-to-fixed-point
  // produces another. If they differ by > 1% we flag it as a regression
  // (per the design doc's quantitative oracle recommendation). Catches
  // thin-polygon-style collapses where iteration empties a non-empty
  // pass-1 output, which topology balance alone misses.
  int area_pass1_nonzero = 0;       // pass-1 had a measurable area
  int area_drift_over_1pct = 0;     // > 1% drift between pass 1 and converged
  int area_drift_over_10pct = 0;    // > 10% drift (typically a collapse)
  double area_drift_max = 0.0;      // largest fractional drift seen
  std::vector<std::tuple<int, int, uint64_t, double>> areaDriftList;
  std::vector<std::tuple<int, int, uint64_t>> topoMismatch;
  std::vector<std::tuple<int, int, uint64_t, int>> convergedInvalidList;
  std::vector<std::tuple<int, int, uint64_t, int>> degradedList;
  // (kPow, n, seed, iters) for both
  std::vector<std::tuple<int, int, uint64_t>> geomCollapseList;

  for (int kPow : kPows) {
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    for (int n : sizes) {
      for (int seed = 0; seed < seedsPerCell; ++seed) {
        ++total;
        // Seed mapping matches the standard displacement fuzz in main()
        // so any failing (kPow, n, seed) here can be passed verbatim to
        // `./overlap2d_proto diagnose <seed> <kPow> <n>`.
        auto [v, e] = RandomTopologicalPolygon(
            n, static_cast<uint64_t>(seed) + 1000ull * kPow);
        Displace(&v, offset);

        auto pass1 = RemoveOverlaps2D(v, e, eps);
        const bool pass1Valid = CheckTopologicalValidity(
            pass1, e, pass1.inputRemap, pass1.numMergedVerts);
        if (!pass1Valid) ++firstPassFail;

        // Polygons round-trip: convert (v, e) (which is a single closed
        // cycle from RandomTopologicalPolygon) into manifold::Polygons,
        // run the public Simplify API, and check the
        // boundary-edge count matches the lower-level RemoveOverlaps2D
        // output. Catches regressions in PolygonsToInput's drop-degenerate
        // logic and OutEdgesToPolygons's loop-walking.
        {
          // Walk the edge cycle (v0 -> v1 chain) into a single SimplePolygon.
          std::vector<int> nextV(v.size(), -1);
          for (const auto& edge : e) nextV[edge.v0] = edge.v1;
          manifold::SimplePolygon loop;
          int cur = e.empty() ? -1 : e[0].v0;
          for (size_t step = 0; step < v.size() && cur >= 0; ++step) {
            loop.push_back(v[cur]);
            cur = nextV[cur];
            if (cur == e[0].v0) break;
          }
          manifold::Polygons polys = {std::move(loop)};
          auto out = Simplify(polys, eps);
          // Compare lower-level edge count to post-regularization count.
          // A drop indicates OutEdgesToPolygons regularized one or more
          // zero-area lens loops out of the result, which is correct
          // CGAL/Clipper2/SVG behavior; track the rate as a stat, not
          // a bug.
          size_t outEdges = 0;
          for (const auto& l : out) outEdges += l.size();
          if (outEdges == pass1.edges.size())
            ++regUnchanged;
          else
            ++regDropped;
        }

        // Run iterate-to-fixed-point to classify.
        int iters = 0;
        IterStatus status = IterStatus::Converged;
        auto rIter =
            IterateToFixedPoint(v, e, eps, /*maxIter=*/8, &iters, &status);
        ++iterDist[iters];
        // Note: rIter.inputRemap is the LAST iteration's remap, not the
        // composed original->final. The fuzz inputs are random topological
        // polygons where step 1 produces an identity remap (no input verts
        // merge), so this is OK here.
        const bool convergedValid = CheckTopologicalValidity(
            rIter, e, rIter.inputRemap, rIter.numMergedVerts);
        if (convergedValid) {
          ++convergedPassValid;
          if (!pass1Valid) ++iter_repaired;
        } else {
          ++convergedPassInvalid;
          convergedInvalidList.emplace_back(kPow, n,
                                            static_cast<uint64_t>(seed),
                                            iters);
          if (pass1Valid) {
            ++iter_degraded;
            degradedList.emplace_back(kPow, n,
                                      static_cast<uint64_t>(seed), iters);
          } else {
            ++pass1_invalid_unfixed;
          }
        }

        // Area-drift check: compare pass 1's signed area against the
        // converged result's signed area. If pass 1 was meaningful
        // (non-empty, non-zero area) and the iteration changed it by
        // more than 1%, flag as a regression candidate. Special case:
        // pass-1-nonempty -> converged-empty is the thin-polygon
        // collapse mode (drift = 100%, area_drift_over_10pct picks it up).
        const double pass1Area =
            SignedAreaFromOutEdges(pass1.verts, pass1.edges);
        const double convArea =
            SignedAreaFromOutEdges(rIter.verts, rIter.edges);
        if (std::fabs(pass1Area) > 0) {
          ++area_pass1_nonzero;
          const double drift =
              std::fabs(convArea - pass1Area) / std::fabs(pass1Area);
          if (drift > area_drift_max) area_drift_max = drift;
          if (drift > 0.01) {
            ++area_drift_over_1pct;
            areaDriftList.emplace_back(kPow, n,
                                       static_cast<uint64_t>(seed), drift);
          }
          if (drift > 0.10) ++area_drift_over_10pct;
        }

        if (iters >= 2) {
          ++iterGE2;
          if (pass1Valid)
            ++iterGE2_pass1_valid;
          else
            ++iterGE2_pass1_invalid;
          // Hypothesis check: pass1 vs pass2 differ at eps/100 but match
          // at eps quantum (positions drift, topology doesn't)?
          std::vector<EdgeM> p2edges;
          for (const auto& oe : pass1.edges)
            p2edges.push_back({oe.v0, oe.v1, oe.mult});
          auto pass2 = RemoveOverlaps2D(pass1.verts, p2edges, eps);
          if (!pass1.edges.empty()) ++pass1_nonempty;
          if (!pass1.edges.empty() && pass2.edges.empty()) {
            ++pass2_collapsed;
            geomCollapseList.emplace_back(kPow, n,
                                          static_cast<uint64_t>(seed));
          }
          const std::string fp1_fine = FingerprintAt(pass1, eps * 0.01);
          const std::string fp2_fine = FingerprintAt(pass2, eps * 0.01);
          const std::string fp1_topo = CoarseFingerprint(pass1, eps);
          const std::string fp2_topo = CoarseFingerprint(pass2, eps);
          if (fp1_fine == fp2_fine) ++iterGE2_fine_match;
          if (fp1_topo == fp2_topo)
            ++iterGE2_topo_match;
          else {
            ++iterGE2_topo_diff;
            topoMismatch.emplace_back(kPow, n,
                                      static_cast<uint64_t>(seed));
          }
        }
      }
    }
  }

  std::cout << "  Total cases: " << total << "\n";
  std::cout << "  First-pass topology FAIL: " << firstPassFail << "\n";
  std::cout << "  Converged result topology valid:   " << convergedPassValid
            << "\n";
  std::cout << "  Converged result topology INVALID: " << convergedPassInvalid
            << "\n";
  std::cout << "    breakdown:\n";
  std::cout << "      iteration REPAIRED (invalid→valid):   " << iter_repaired
            << "\n";
  std::cout << "      iteration DEGRADED (valid→invalid):   " << iter_degraded
            << "\n";
  std::cout << "      pass-1 invalid + iteration didn't fix: "
            << pass1_invalid_unfixed << "\n";
  std::cout << "  Iter distribution: ";
  for (auto [k, c] : iterDist) std::cout << k << ":" << c << " ";
  std::cout << "\n";
  std::cout << "  iter≥2 cases: " << iterGE2 << "\n";
  std::cout << "    of which pass 1 was already topo-valid: "
            << iterGE2_pass1_valid << "\n";
  std::cout << "    of which pass 1 was topo-INVALID (iteration repaired): "
            << iterGE2_pass1_invalid << "\n";
  std::cout << "    where pass1 == pass2 at eps quantum (TOPO match): "
            << iterGE2_topo_match << "\n";
  std::cout << "    where pass1 == pass2 at eps/100 quantum (FINE match): "
            << iterGE2_fine_match << "\n";
  std::cout << "    where TOPO differs (real algorithm change): "
            << iterGE2_topo_diff << "\n";
  std::cout << "  Geometric collapse: pass1 nonempty=" << pass1_nonempty
            << " of which pass2=empty: " << pass2_collapsed
            << " (output orientation bug)\n";
  std::cout << "  Polygons regularization (zero-area lens drops): "
            << regDropped << " of " << total
            << " cases regularized away one or more zero-area loops "
               "(correct CGAL/Clipper2/SVG behavior; not a bug)\n";
  std::cout << "  Area drift (pass-1 vs converged, |delta|/|pass1|):\n";
  std::cout << "    pass-1 with measurable area: " << area_pass1_nonzero
            << "\n";
  std::cout << "    drift > 1%:  " << area_drift_over_1pct << "\n";
  std::cout << "    drift > 10%: " << area_drift_over_10pct
            << "  (typically thin-polygon collapse)\n";
  std::cout << "    max drift seen: " << (area_drift_max * 100.0) << "%\n";
  // For idempotence diagnosis: dump first 5 iter=2 cases (pass 1 valid).
  if (!topoMismatch.empty()) {
    std::cout << "\n  First 5 iter≥2 cases (for idempotence probe):\n";
    for (size_t i = 0; i < topoMismatch.size() && i < 5; ++i) {
      auto [kp, nn, sd] = topoMismatch[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd << "\n";
    }
  }
  // Geometric collapse cases: pass 1 nonempty, pass 2 empty. These are
  // the genuine bugs where iteration is changing the algorithm output.
  // Track them by re-running and checking. (Slow: requires re-running
  // pass 2 for each case. Done only for first 30 cases.)
  if (!geomCollapseList.empty()) {
    std::cout << "\n  Geometric COLLAPSE cases (pass 2 empty when pass 1 "
                 "nonempty):\n";
    for (size_t i = 0; i < geomCollapseList.size() && i < 20; ++i) {
      auto [kp, nn, sd] = geomCollapseList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd << "\n";
    }
  }
  if (!degradedList.empty()) {
    std::cout << "\n  DEGRADED cases (pass 1 valid, iteration broke it; first 20):\n";
    for (size_t i = 0; i < degradedList.size() && i < 20; ++i) {
      auto [kp, nn, sd, it] = degradedList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << " iters=" << it << "\n";
    }
  }
  if (!areaDriftList.empty()) {
    std::cout << "\n  Area-drift cases > 1% (first 10):\n";
    // Sort by drift magnitude so the worst cases surface first.
    std::sort(areaDriftList.begin(), areaDriftList.end(),
              [](const auto& a, const auto& b) {
                return std::get<3>(a) > std::get<3>(b);
              });
    for (size_t i = 0; i < areaDriftList.size() && i < 10; ++i) {
      auto [kp, nn, sd, dr] = areaDriftList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << "  drift=" << (dr * 100.0) << "%\n";
    }
  }
  if (!convergedInvalidList.empty()) {
    std::cout
        << "\n  CONVERGED-INVALID cases (iteration didn't fix; first 30):\n";
    for (size_t i = 0; i < convergedInvalidList.size() && i < 30; ++i) {
      auto [kp, nn, sd, it] = convergedInvalidList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << " iters=" << it << "\n";
    }
  }
}
}  // namespace overlap2d

int main(int argc, char** argv) {
  using namespace overlap2d;
  if (argc > 1 && std::string(argv[1]) == "diagnose") {
    uint64_t seed = (argc > 2) ? std::stoull(argv[2]) : 0;
    int kPow = (argc > 3) ? std::atoi(argv[3]) : 30;
    int n = (argc > 4) ? std::atoi(argv[4]) : 50;
    Diagnose(seed, kPow, n);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "deepfuzz") {
    int seedsPerCell = (argc > 2) ? std::atoi(argv[2]) : 200;
    DeepFuzz(seedsPerCell);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "corpus") {
    const std::string path = (argc > 2)
                                  ? argv[2]
                                  : "test/polygons/polygon_corpus.txt";
    RunCorpus(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "clipper2corpus") {
    std::string path =
        "build/_deps/clipper2-src/Tests/Polygons.txt";
    if (argc > 2) path = argv[2];
    RunClipper2Corpus(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "mfogelcorpus") {
    std::string dir =
        "build/_deps/mfogel-polygon-clipping-src/test/end-to-end";
    if (argc > 2) dir = argv[2];
    RunMfogelCorpus(dir);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtscorpus") {
    std::string path = "test/polygons/jts_overlay_corpus.txt";
    if (argc > 2) path = argv[2];
    RunJtsCorpus(path);
    return 0;
  }
#ifdef OVERLAP2D_WITH_CLIPPER2
  if (argc > 1 && std::string(argv[1]) == "vsclipper2") {
    std::string which = (argc > 2) ? argv[2] : "all";
    RunVsClipper2(which);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtsdrops") {
    std::string path = (argc > 2) ? argv[2] : "test/polygons/jts_overlay_corpus.txt";
    RunJtsDropDiag(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtsdiag") {
    if (argc < 3) {
      std::cerr << "Usage: jtsdiag <case_n>\n";
      return 1;
    }
    int targetCase = std::atoi(argv[2]);
    auto cases = LoadJtsCorpus("test/polygons/jts_overlay_corpus.txt");
    for (const auto& c : cases) {
      if (c.n != targetCase) continue;
      std::cout << "=== JTS case " << c.n << " ===\n";
      std::cout << "  op=" << c.op << " src=" << c.source << "\n";
      std::cout << "  A.rings=" << c.a.size() << " B.rings=" << c.b.size() << "\n";
      double xMin = 1e300, yMin = 1e300, xMax = -1e300, yMax = -1e300;
      int totalVerts = 0;
      auto bound = [&](const manifold::Polygons& polys) {
        for (const auto& loop : polys) {
          totalVerts += loop.size();
          for (const auto& v : loop) {
            xMin = std::min(xMin, v.x); yMin = std::min(yMin, v.y);
            xMax = std::max(xMax, v.x); yMax = std::max(yMax, v.y);
          }
        }
      };
      bound(c.a); bound(c.b);
      const double eps = InferEps(c.a, c.b);
      std::cout << "  bbox: x=[" << xMin << ".." << xMax << "] y=[" << yMin << ".." << yMax << "]\n";
      std::cout << "  total verts: " << totalVerts << " eps: " << eps << "\n";
      // Run our pipeline with debug=true.
      const manifold::Polygons& A = c.a;
      const manifold::Polygons& B = (c.op == "overlayAreaTest" || c.op == "overlayareatest") ? c.b : manifold::Polygons{};
      std::vector<vec2> verts;
      std::vector<EdgeM> edges;
      auto append = [&](const manifold::Polygons& polys) {
        for (const auto& loop : polys) {
          if (loop.size() < 3) continue;
          const int base = static_cast<int>(verts.size());
          const int n = static_cast<int>(loop.size());
          for (const auto& v : loop) verts.push_back(v);
          for (int i = 0; i < n; ++i)
            edges.push_back({base + i, base + ((i + 1) % n), 1});
        }
      };
      append(A); append(B);
      std::cout << "  Pipeline input: " << verts.size() << " verts, " << edges.size() << " edges\n";
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/true,
                                [](int w) { return w != 0; });
      auto out = OutEdgesToPolygons(r.verts, r.edges);
      std::cout << "  Pipeline output: " << r.verts.size() << " verts, " << r.edges.size() << " edges, " << out.size() << " polygons\n";
      // Clipper2 reference for comparison.
      auto a2 = ToPaths64(A);
      auto b2 = ToPaths64(B);
      Clipper2Lib::Clipper64 clip;
      clip.AddSubject(a2);
      if (!B.empty()) clip.AddClip(b2);
      Clipper2Lib::Paths64 sol;
      clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero, sol);
      double clipArea = 0; for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      std::cout << "  Clipper2 output: " << sol.size() << " paths, area=" << clipArea << "\n";
      return 0;
    }
    std::cerr << "Case " << targetCase << " not found\n";
    return 1;
  }
#endif
  if (argc > 1 && std::string(argv[1]) == "time") {
    // Per-phase wall-clock timing on a representative workload. Default
    // workload is DeepFuzz with seedsPerCell=200 (5,600 cases × ~3
    // RemoveOverlaps2D calls each ≈ 17k pipeline runs); selectable by a
    // second positional arg matching another mode name. Examples:
    //
    //   ./overlap2d_proto time                 # deepfuzz 200
    //   ./overlap2d_proto time deepfuzz 1000   # 28k cases
    //   ./overlap2d_proto time jtscorpus       # JTS overlay corpus
    //   ./overlap2d_proto time clipper2corpus
    //   ./overlap2d_proto time mfogelcorpus
    //   ./overlap2d_proto time corpus          # polygon_corpus.txt
    //
    // Output: total wall-clock + per-phase ns + percentage breakdown.
    GlobalPhases().Reset();
    const auto wallStart = std::chrono::steady_clock::now();
    std::string sub = (argc > 2) ? argv[2] : "deepfuzz";
    // Silence per-mode stdout so timing summary is uncluttered.
    std::ostringstream sink;
    std::streambuf* oldCout = std::cout.rdbuf(sink.rdbuf());
    if (sub == "deepfuzz") {
      int seeds = (argc > 3) ? std::atoi(argv[3]) : 200;
      DeepFuzz(seeds);
    } else if (sub == "corpus") {
      const std::string p = (argc > 3) ? argv[3] : "test/polygons/polygon_corpus.txt";
      RunCorpus(p);
    } else if (sub == "clipper2corpus") {
      const std::string p = (argc > 3) ? argv[3]
        : std::string("build/_deps/clipper2-src/Tests/Polygons.txt");
      RunClipper2Corpus(p);
    } else if (sub == "mfogelcorpus") {
      const std::string p = (argc > 3) ? argv[3]
        : std::string("build/_deps/mfogel-polygon-clipping-src/test/end-to-end");
      RunMfogelCorpus(p);
    } else if (sub == "jtscorpus") {
      const std::string p = (argc > 3) ? argv[3] : "test/polygons/jts_overlay_corpus.txt";
      RunJtsCorpus(p);
    } else {
      std::cout.rdbuf(oldCout);
      std::cerr << "unknown subcommand: " << sub << "\n";
      return 2;
    }
    std::cout.rdbuf(oldCout);
    const auto wallEnd = std::chrono::steady_clock::now();
    const int64_t wallNs = std::chrono::duration_cast<std::chrono::nanoseconds>(
        wallEnd - wallStart).count();
    const auto& P = GlobalPhases();
    const int64_t total = P.totalNs.load();
    const int64_t cases = P.cases.load();
    auto pct = [&](int64_t ns) {
      return total > 0 ? (ns * 100.0 / total) : 0.0;
    };
    auto fmt = [](int64_t ns) {
      double ms = ns * 1e-6;
      std::ostringstream o;
      o.setf(std::ios::fixed);
      o.precision(3);
      o << ms << " ms";
      return o.str();
    };
    std::cout << "=== timing: " << sub << " ===\n";
    std::cout << "  Wall (driver+pipeline): " << fmt(wallNs) << "\n";
    std::cout << "  Pipeline only (sum):    " << fmt(total) << "\n";
    std::cout << "  RemoveOverlaps2D calls: " << cases << "\n";
    if (cases > 0) {
      std::cout << "  Avg per call:           " << fmt(total / cases) << "\n";
    }
    std::cout << "\n  Per-phase breakdown (% of pipeline total):\n";
    auto row = [&](const char* name, int64_t ns) {
      std::cout << "    " << std::left << std::setw(28) << name
                << std::right << std::setw(12) << fmt(ns)
                << "  " << std::setw(6);
      std::cout.setf(std::ios::fixed);
      std::cout.precision(2);
      std::cout << pct(ns) << "%\n";
      std::cout.unsetf(std::ios::fixed);
    };
    row("step 1  MergeVerts",          P.mergeNs.load());
    row("step 2  RemapAndCollapse",    P.remapNs.load());
    row("step 3  BuildEdgeVertLists",  P.buildListsNs.load());
    row("step 4  FindAndInsertIxs",    P.findIxNs.load());
    row("step 4b structural re-merge", P.restructNs.load());
    row("step 5  Canonicalize",        P.canonNs.load());
    row("step 6  FilterByWindingDCEL", P.filterDcelNs.load());
    std::cout << "\n  Sub-phase breakdown (% of pipeline total):\n";
    row("  step 1 BVH build",           P.step1BvhBuildNs.load());
    row("  step 1 BVH self-collide",    P.step1CollideNs.load());
    row("  step 3+4 BVH build",         P.bvhBuildNs.load());
    row("  step 3+4-broad concurrent",  P.step3BroadNs.load());
    row("  step 4 narrow",              P.step4NarrowNs.load());
    row("  step 4 propagation",         P.step4PropNs.load());
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "pentagon") {
    // Sanity: feed a clean pentagon (displaced if requested), see what DCEL
    // does. With kPow arg, displaces to that scale to mimic re-feed conditions.
    int kPow = (argc > 2) ? std::atoi(argv[2]) : 0;
    const double offset = (kPow > 0) ? std::ldexp(1.5, kPow) : 0.0;
    const double scale = (kPow > 0) ? offset : 1.0;
    std::vector<vec2> v = {
        {std::cos(0.0), std::sin(0.0)},
        {std::cos(2 * M_PI / 5), std::sin(2 * M_PI / 5)},
        {std::cos(4 * M_PI / 5), std::sin(4 * M_PI / 5)},
        {std::cos(6 * M_PI / 5), std::sin(6 * M_PI / 5)},
        {std::cos(8 * M_PI / 5), std::sin(8 * M_PI / 5)},
    };
    if (kPow > 0)
      for (auto& p : v) {
        p.x += offset;
        p.y += offset;
      }
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1},
                            {3, 4, 1}, {4, 0, 1}};
    auto r = RemoveOverlaps2D(v, e, EpsilonFromScale(scale));
    std::cerr << "Pentagon (kPow=" << kPow << "): " << r.verts.size()
              << " verts, " << r.edges.size() << " edges\n";
    for (const auto& oe : r.edges) {
      std::cerr << "  " << oe.v0 << "→" << oe.v1 << " mult=" << oe.mult
                << "\n";
    }
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "idempotence") {
    // ./overlap2d_proto idempotence <seed> [kPow] [n]
    // Compares pass 1 vs pass 2 sub-edge sets to identify what's drifting.
    uint64_t seed = (argc > 2) ? std::stoull(argv[2]) : 0;
    int kPow = (argc > 3) ? std::atoi(argv[3]) : 10;
    int n = (argc > 4) ? std::atoi(argv[4]) : 8;
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    auto [vIn, eIn] = RandomTopologicalPolygon(n, seed + 1000ull * kPow);
    Displace(&vIn, offset);
    std::cerr << "=== idempotence: kPow=" << kPow << " n=" << n
              << " seed=" << seed << " eps=" << eps << " ===\n";
    std::cerr << "--- pass 1 ---\n";
    auto pass1 = RemoveOverlaps2D(vIn, eIn, eps, /*debug=*/true);
    std::vector<EdgeM> p2in;
    for (const auto& oe : pass1.edges)
      p2in.push_back({oe.v0, oe.v1, oe.mult});
    std::cerr << "--- pass 2 (input: " << pass1.verts.size() << " verts, "
              << p2in.size() << " edges) ---\n";
    // Manually run pass 2's pipeline to see where edges disappear.
    auto p2_mrg = MergeVerts(pass1.verts, eps);
    std::cerr << "  step 1 mergeverts: " << pass1.verts.size() << "→"
              << p2_mrg.verts.size() << "\n";
    auto p2_e = RemapAndCollapse(p2in, p2_mrg.remap);
    std::cerr << "  step 2 remap: " << p2in.size() << "→" << p2_e.size()
              << " edges\n";
    std::vector<Box> p2_eBoxes(p2_e.size());
    for (size_t i = 0; i < p2_e.size(); ++i)
      p2_eBoxes[i] = BoxOf2DEdge(p2_mrg.verts[p2_e[i].v0],
                                 p2_mrg.verts[p2_e[i].v1], eps);
    BVH p2_bvh = BVHBuildFromBoxes(p2_eBoxes);
    auto p2_l = BuildEdgeVertLists(p2_e, p2_mrg.verts, eps, p2_eBoxes, p2_bvh);
    int p2_totalList = 0;
    for (auto& l : p2_l) p2_totalList += l.size();
    std::cerr << "  step 3 lists: " << p2_totalList << " total entries\n";
    std::vector<std::vector<int>> p2_ve;
    auto p2_pairs = CollectStep4Pairs(p2_e, p2_eBoxes, p2_bvh);
    FindAndInsertIntersections(p2_e, &p2_mrg.verts, &p2_l, &p2_ve, eps,
                               p2_eBoxes, p2_bvh, p2_pairs);
    std::cerr << "  step 4 intersections: " << p2_mrg.verts.size()
              << " verts after\n";
    auto p2_canon = Canonicalize(p2_e, p2_l);
    std::cerr << "  step 5 canon: " << p2_canon.edges.size()
              << " sub-edges\n";

    auto pass2 = RemoveOverlaps2D(pass1.verts, p2in, eps, /*debug=*/true);
    std::cerr << "pass1: " << pass1.verts.size() << " verts, "
              << pass1.edges.size() << " edges\n";
    std::cerr << "pass2: " << pass2.verts.size() << " verts, "
              << pass2.edges.size() << " edges\n";
    // Dump pass 1 edges as (qx, qy)→(qx, qy) at eps/100 quantum.
    auto qedges = [&](const OverlapResult& r) {
      const double q = eps * 0.01;
      auto qq = [q](double x) { return (int64_t)std::round(x / q); };
      std::set<std::tuple<int64_t, int64_t, int64_t, int64_t, int>> s;
      for (const auto& oe : r.edges) {
        const vec2 p0 = r.verts[oe.v0], p1 = r.verts[oe.v1];
        auto a = std::make_pair(qq(p0.x), qq(p0.y));
        auto b = std::make_pair(qq(p1.x), qq(p1.y));
        int m = oe.mult;
        if (b < a) {
          std::swap(a, b);
          m = -m;
        }
        s.emplace(a.first, a.second, b.first, b.second, m);
      }
      return s;
    };
    auto e1 = qedges(pass1), e2 = qedges(pass2);
    std::set<std::tuple<int64_t, int64_t, int64_t, int64_t, int>> only1, only2;
    std::set_difference(e1.begin(), e1.end(), e2.begin(), e2.end(),
                        std::inserter(only1, only1.end()));
    std::set_difference(e2.begin(), e2.end(), e1.begin(), e1.end(),
                        std::inserter(only2, only2.end()));
    std::cerr << "Edges in pass1 but not pass2: " << only1.size() << "\n";
    for (const auto& [ax, ay, bx, by, m] : only1) {
      std::cerr << "  (" << ax << "," << ay << ")→(" << bx << "," << by
                << ") mult=" << m << "\n";
    }
    std::cerr << "Edges in pass2 but not pass1: " << only2.size() << "\n";
    for (const auto& [ax, ay, bx, by, m] : only2) {
      std::cerr << "  (" << ax << "," << ay << ")→(" << bx << "," << by
                << ") mult=" << m << "\n";
    }
    return 0;
  }
  bool allPass = true;

  // ---- Simple sanity ----
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
    auto r = RunCase({"CCW square (no overlap)", v, e, EpsilonFromScale(1.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Boolean2D entry point: Add / Subtract / Intersect ----
  // Two CCW unit squares overlapping by a 0.5×1 strip in the middle.
  {
    manifold::Polygons a = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    manifold::Polygons b = {{{0.5, 0}, {1.5, 0}, {1.5, 1}, {0.5, 1}}};
    auto add = Boolean2D(a, b, OpType::Add);
    auto sub = Boolean2D(a, b, OpType::Subtract);
    auto isec = Boolean2D(a, b, OpType::Intersect);
    // Expected boundary edge counts (each square is 4 edges; expected
    // post-overlap: union 4 edges (rectangle 0..1.5 by 0..1 with 4 sides
    // but each side has a midpoint vert from the cross-cut so actually
    // 6 boundary edges); difference 4 edges (square 0..0.5 by 0..1 with
    // 4 sides), intersection 4 edges (square 0.5..1 by 0..1 with 4
    // sides). Allow some flex from intersection-vert insertion.
    auto edgeCount = [](const manifold::Polygons& p) {
      size_t s = 0;
      for (const auto& l : p) s += l.size();
      return s;
    };
    std::cout << "=== Boolean2D: two overlapping squares ===\n";
    std::cout << "  Add:       " << add.size() << " loop(s), "
              << edgeCount(add) << " edges\n";
    std::cout << "  Subtract:  " << sub.size() << " loop(s), "
              << edgeCount(sub) << " edges\n";
    std::cout << "  Intersect: " << isec.size() << " loop(s), "
              << edgeCount(isec) << " edges\n";
    bool ok = !add.empty() && !sub.empty() && !isec.empty() &&
              edgeCount(add) >= 4 && edgeCount(sub) >= 4 &&
              edgeCount(isec) >= 4;
    std::cout << "  Smoke: " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  // Two axis-aligned CCW squares overlapping at an L-corner. Their
  // boundaries cross at TWO perpendicular axis-aligned edge pairs
  // (A's right edge × B's bottom edge at (2,1); A's top edge × B's
  // left edge at (1,2)). The Boolean2D smoke test above offsets the
  // squares horizontally so their crossing edges are *collinear*
  // (both top edges at y=1), not perpendicular - it doesn't exercise
  // axis-aligned perpendicular intersection at all. This test does.
  {
    std::vector<vec2> v = {{0, 0}, {2, 0}, {2, 2}, {0, 2},   // square A CCW
                           {1, 1}, {3, 1}, {3, 3}, {1, 3}};  // square B CCW
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // A
        {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},  // B
    };
    auto r = RunCase(
        {"Two axis-aligned squares overlapping (L-shape union)", v, e,
         EpsilonFromScale(3.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(3.0))) allPass = false;
    // Expected union: single L-shaped loop with 8 boundary verts and
    // total area 4 + 4 − 1 = 7. Step 4 must detect the two
    // perpendicular axis-aligned intersections.
    manifold::Polygons in = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                             {{1, 1}, {3, 1}, {3, 3}, {1, 3}}};
    auto out = Simplify(in, EpsilonFromScale(3.0));
    size_t totalEdges = 0;
    for (const auto& l : out) totalEdges += l.size();
    bool shapeOk = (out.size() == 1 && totalEdges == 8);
    std::cout << "  Polygons union: " << out.size() << " loop(s), "
              << totalEdges << " edges, "
              << (shapeOk ? "PASS" : "FAIL") << "\n";
    if (!shapeOk) allPass = false;
    std::cout << std::endl;
  }

  // Two CCW unit squares overlapping diagonally, run through all three
  // Boolean2D ops with quantitative area checks. Smith fig 6.2 in
  // spirit: a single configuration that exercises every winding rule.
  // Critically, Intersect drives the `w > 1` predicate (the central
  // 1x1 overlap region has winding 2 because both inputs cover it).
  // The existing Boolean2D smoke test (two squares offset *horizontally*)
  // produces a wind=2 strip but its boundaries are collinear, not
  // perpendicular, and step-6 Intersect was untested against
  // axis-aligned perpendicular crossings until this case.
  {
    manifold::Polygons a = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}}};   // CCW
    manifold::Polygons b = {{{1, 1}, {3, 1}, {3, 3}, {1, 3}}};   // CCW
    const double eps = EpsilonFromScale(3.0);
    auto add = Boolean2D(a, b, OpType::Add, eps);
    auto sub = Boolean2D(a, b, OpType::Subtract, eps);
    auto isec = Boolean2D(a, b, OpType::Intersect, eps);
    auto xorOp = Xor(a, b, eps);
    const double areaAdd = TotalSignedArea(add);
    const double areaSub = TotalSignedArea(sub);
    const double areaIsec = TotalSignedArea(isec);
    const double areaXor = TotalSignedArea(xorOp);
    // Expected: A area 4, B area 4, overlap 1.
    // Add: 4 + 4 - 1 = 7.  Subtract: 4 - 1 = 3.  Intersect: 1.
    // Xor (symmetric difference): (A | B) - (A & B) = 7 - 1 = 6.
    auto near = [](double a, double b) { return std::fabs(a - b) < 1e-9; };
    bool ok = near(areaAdd, 7.0) && near(areaSub, 3.0) &&
              near(areaIsec, 1.0) && near(areaXor, 6.0);
    std::cout << "=== Boolean2D area regression: diagonal squares ===\n";
    std::cout << "  Add area:       " << areaAdd << " (expect 7)\n";
    std::cout << "  Subtract area:  " << areaSub << " (expect 3)\n";
    std::cout << "  Intersect area: " << areaIsec << " (expect 1)\n";
    std::cout << "  Xor area:       " << areaXor << " (expect 6)\n";
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = LShape();
    auto r = RunCase({"L-shape (concave, CCW)", v, e, EpsilonFromScale(2.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(2.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = SquareWithHole();
    auto r = RunCase(
        {"Square with hole (outer CCW, inner CW)", v, e, EpsilonFromScale(4.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(4.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = TwoSquares();
    auto r = RunCase({"Two disjoint CCW squares", v, e, EpsilonFromScale(3.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(3.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Polygons-API smoke tests ----
  // Verify the public Simplify wrapper round-trips simple
  // closed loops through the manifold::Polygons type.
  {
    std::cout << "=== Polygons API: square ===" << std::endl;
    manifold::Polygons in = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    auto out = Simplify(in, EpsilonFromScale(1.0));
    bool ok = (out.size() == 1 && out[0].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    if (!out.empty()) std::cout << ", " << out[0].size() << " verts in loop 0";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: square with hole ===" << std::endl;
    manifold::Polygons in = {
        {{0, 0}, {4, 0}, {4, 4}, {0, 4}},          // outer CCW
        {{1, 1}, {1, 3}, {3, 3}, {3, 1}},          // inner CW (hole)
    };
    auto out = Simplify(in, EpsilonFromScale(4.0));
    // Expect two loops (outer + hole) with 4 verts each.
    bool ok = (out.size() == 2 && out[0].size() == 4 && out[1].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: bowtie keeps the wind>0 lobe ==="
              << std::endl;
    // Self-intersecting bowtie input: under Smith's `wind > 0`
    // convention, only the CCW lobe is in the boolean union; the
    // other lobe has wind = −1. So output = 1 triangle.
    manifold::Polygons in = {{{0, 0}, {2, 2}, {2, 0}, {0, 2}}};
    auto out = Simplify(in, EpsilonFromScale(2.0));
    bool ok = (out.size() == 1 && out[0].size() == 3);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: two disjoint squares ===" << std::endl;
    manifold::Polygons in = {
        {{0, 0}, {1, 0}, {1, 1}, {0, 1}},          // square A (CCW)
        {{2, 0}, {3, 0}, {3, 1}, {2, 1}},          // square B (CCW)
    };
    auto out = Simplify(in, EpsilonFromScale(3.0));
    bool ok = (out.size() == 2 && out[0].size() == 4 && out[1].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  // ---- Non-closed / non-manifold input tests ----
  // The prototype was designed for closed polygons. These tests probe
  // what happens at the edges of that assumption: open polylines,
  // shared verts/edges between polygons, T-junctions, etc.
  //
  // For *open* inputs (polylines, T-junctions) the per-vertex balance
  // check `CheckTopologicalValidity` doesn't apply: it compares input
  // edge contributions per vertex against output edge contributions,
  // and for open inputs the input has end-vertex balance ±1 while the
  // algorithm's correct output is empty (balance 0). So we accept
  // those as observational and only fail closed-input cases.
  enum class NMExpect { ClosedTopo, OpenObservational };
  auto runNonManifold = [&allPass](
                            const char* name, const std::vector<vec2>& v,
                            const std::vector<EdgeM>& e, double eps,
                            const char* expectation, NMExpect mode) {
    std::cout << "=== Non-manifold: " << name << " ===" << std::endl;
    std::cout << "  Input: " << v.size() << " verts, " << e.size() << " edges, "
              << "expectation: " << expectation << std::endl;
    auto r = RemoveOverlaps2D(v, e, eps);
    std::cout << "  Output: " << r.verts.size() << " verts, " << r.edges.size()
              << " edges" << std::endl;
    bool topoOk =
        CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
    if (mode == NMExpect::ClosedTopo) {
      std::cout << "  Topological validity: " << (topoOk ? "PASS" : "FAIL")
                << std::endl;
      if (!topoOk) allPass = false;
    } else {
      std::cout << "  Topological validity: " << (topoOk ? "PASS" : "FAIL")
                << " (observational; open input, balance check N/A)"
                << std::endl;
    }
    std::cout << std::endl;
    return r;
  };

  // (1) Open polyline (zigzag): no closure, no enclosed region.
  // Expected: empty output (wind > 0 is empty everywhere).
  {
    std::vector<vec2> v = {{0, 0}, {1, 1}, {2, 0}, {3, 1}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
    runNonManifold("open polyline (Z-zigzag)", v, e, EpsilonFromScale(3.0),
                   "empty (no enclosed region)",
                   NMExpect::OpenObservational);
  }

  // (2) Two polygons touching at a single vertex (kissing corner).
  // Squares share v=(1,1); each is otherwise CCW and disjoint.
  // Expected: 2 disjoint loops (8 output edges total) + the shared
  // vertex appears in both.
  {
    std::vector<vec2> v = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1},  // square A
        {2, 1}, {2, 2}, {1, 2},          // square B (shares (1,1))
    };
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // square A CCW
        {2, 4, 1}, {4, 5, 1}, {5, 6, 1}, {6, 2, 1},  // square B CCW (via (1,1))
    };
    runNonManifold("two polys touching at a vertex", v, e,
                   EpsilonFromScale(2.0), "2 disjoint loops",
                   NMExpect::ClosedTopo);
  }

  // (3) Two polygons sharing an edge (non-manifold along the shared
  // segment, two CCW squares glued along (1,0)→(1,1)). The shared
  // edge has canonical mult = +2 in the combined input. Expected: 1
  // outer loop = the rectangle [0,2]×[0,1] (the two squares fused).
  {
    std::vector<vec2> v = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1},  // square A
        {2, 0}, {2, 1},                  // square B's other corners
    };
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // square A CCW
        {1, 4, 1}, {4, 5, 1}, {5, 2, 1}, {2, 1, 1},  // square B CCW
    };
    runNonManifold("two polys sharing an edge", v, e, EpsilonFromScale(2.0),
                   "1 fused outer loop", NMExpect::ClosedTopo);
  }

  // (4) T-junction: vertex with degree 3, no closure. Three edges
  // dangling from a common point. No enclosed region.
  // Expected: empty output.
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}};
    runNonManifold("T-junction (degree-3 dangling)", v, e,
                   EpsilonFromScale(1.0), "empty",
                   NMExpect::OpenObservational);
  }

  // (5) Open polyline crossing itself (figure-8 missing the closing
  // edge). Topology-wise: 5 verts, 4 edges, one self-crossing.
  // Expected: empty output (no closed boundary).
  {
    std::vector<vec2> v = {{0, 0}, {2, 2}, {2, 0}, {0, 2}, {0, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
    // Note: although the input has 5 verts forming an *open* zigzag,
    // verts 0 and 4 share the same coordinate (0,0). Step 1's eps-merge
    // collapses them, turning the input into a *closed* bowtie. So
    // this case actually exercises the closed bowtie path, not an
    // open-polyline path. The topology check applies.
    runNonManifold("polyline w/ duplicate-position endpoints (=> closed bowtie)",
                   v, e, EpsilonFromScale(2.0), "1 triangle (= bowtie lobe)",
                   NMExpect::ClosedTopo);
  }

  // (6) Degenerate triangle: three collinear verts (zero-area input).
  // Expected: empty (no enclosed region).
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {2, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    runNonManifold("collinear triangle (zero area)", v, e,
                   EpsilonFromScale(2.0), "empty", NMExpect::ClosedTopo);
  }

  // (7) Duplicate-vertex polygon: same coordinate listed twice. Step 1
  // should merge them; step 2 should drop the resulting self-loop edge.
  // Expected: 1 valid loop (the remaining triangle).
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {0, 1}, {0, 1}};  // v2 == v3
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
    runNonManifold("polygon with duplicate vertex", v, e, EpsilonFromScale(1.0),
                   "1 triangle (after merge)", NMExpect::ClosedTopo);
  }

  // (8) Adversarial 4+ concurrent edges: 4 line segments all passing
  // through origin at different orientations. Step 4 produces
  // C(4, 2) = 6 pairwise intersections that should all snap to one
  // point in step 4b. Some pairs share no edge in their incidence sets
  // (e.g. seg AB×CD vs seg EF×GH share no input edge), so step 4b's
  // structural gate cannot merge them directly. The iterate-to-
  // fixed-point pass picks up any leftover near-duplicates. Open
  // input → empty output (no enclosed region).
  {
    std::vector<vec2> v = {{-1, 0},     {1, 0},      {0, -1},     {0, 1},
                           {-0.7, -0.7}, {0.7, 0.7},  {-0.7, 0.7}, {0.7, -0.7}};
    std::vector<EdgeM> e = {{0, 1, 1}, {2, 3, 1}, {4, 5, 1}, {6, 7, 1}};
    runNonManifold("4 segments concurrent at origin (asterisk)", v, e,
                   EpsilonFromScale(1.0), "empty (open input)",
                   NMExpect::OpenObservational);
  }

  // (9) Near-equal-coords at displaced scale: two triangle verts
  // separated by a few ULPs at 2^49 displacement. Step 1's eps-radius
  // merge should fold them, collapsing the triangle to empty output
  // (degenerate after merge).
  {
    const double offset = std::ldexp(1.5, 49);
    const double ulp = std::ldexp(1.0, 49 - 52);  // ~1 ULP at offset
    std::vector<vec2> v = {{offset + 1.0, offset},
                           {offset + 1.0 + ulp, offset},
                           {offset, offset + 1.0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    runNonManifold("near-equal coords at 2^49 displacement", v, e,
                   EpsilonFromScale(offset), "empty (degenerate after merge)",
                   NMExpect::ClosedTopo);
  }

  // ---- Self-intersecting cases ----
  {
    std::vector<vec2> v = {{0, 0}, {2, 0}, {0, 2}, {2, 2}};
    std::vector<EdgeM> e = {{0, 3, 1}, {3, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    auto r = RunCase(
        {"Bowtie (one interior intersection)", v, e, EpsilonFromScale(2.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(2.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = RandomTopologicalPolygon(8, 42);
    auto r = RunCase(
        {"Random topological polygon (n=8)", v, e, EpsilonFromScale(1.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = PolygonalStar(7, 3);
    auto r = RunCase({"Heptagram (n=7, skip=3)", v, e, EpsilonFromScale(1.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = PolygonalStar(7, 3);
    const double offset = std::ldexp(1.5, 49);
    Displace(&v, offset);
    RunCase({"Heptagram displaced by 1.5*2^49", v, e, EpsilonFromScale(offset)},
            &allPass);
    std::cout << std::endl;
  }

  // ---- Smith §7.2 figure 6.5(d) variants ----
  // Smith catalogs two inconsistent-edge sets that exhibit the cyclic-
  // ordering failure: {AF, BG, CE} and {AG, BF, CE}. Both embedded in
  // closed self-intersecting hexagons.
  const double u = std::ldexp(1.0, -53);
  std::vector<vec2> smithVerts = {
      {0.5 + 9 * u, 0.5 + 23 * u},   // A=0
      {0.5 + 23 * u, 0.5 + 25 * u},  // B=1
      {9.0, 0.5},                    // C=2
      {9.0, 32.0},                   // E=3
      {32.0, 32.0 - 32 * u},         // F=4
      {32.0, 32.0},                  // G=5
  };
  {
    // A→F→C→E→B→G→A (uses AF, BG, with CE inside)
    std::vector<EdgeM> e = {{0, 4, 1}, {4, 2, 1}, {2, 3, 1},
                            {3, 1, 1}, {1, 5, 1}, {5, 0, 1}};
    auto r = RunCase({"Smith fig 6.5(d) hexagon {AF, BG, CE}", smithVerts, e,
                      EpsilonFromScale(32.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(32.0))) allPass = false;
    std::cout << std::endl;
  }
  {
    // A→G→C→E→B→F→A (uses AG, BF, with CE inside)
    std::vector<EdgeM> e = {{0, 5, 1}, {5, 2, 1}, {2, 3, 1},
                            {3, 1, 1}, {1, 4, 1}, {4, 0, 1}};
    auto r = RunCase({"Smith fig 6.5(d) hexagon {AG, BF, CE}", smithVerts, e,
                      EpsilonFromScale(32.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(32.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Fuzz: random topological polygons at multiple seeds and sizes ----
  std::cout << "=== Fuzz: random topological polygons ===" << std::endl;
  int fuzzPass = 0, fuzzFail = 0;
  for (int n : {5, 8, 12, 20, 50, 100, 200}) {
    for (uint64_t seed = 0; seed < 50; ++seed) {
      auto [v, e] = RandomTopologicalPolygon(n, seed);
      auto r = RemoveOverlaps2D(v, e, EpsilonFromScale(1.0));
      bool ok = CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
      if (ok) {
        ++fuzzPass;
      } else {
        ++fuzzFail;
        std::cout << "  FAIL: n=" << n << " seed=" << seed << std::endl;
      }
    }
  }
  std::cout << "  Topological invariant: " << fuzzPass << " pass, " << fuzzFail
            << " fail (7 sizes × 50 seeds = 350 cases)" << std::endl;
  if (fuzzFail > 0) allPass = false;

  // ---- Fuzz: displaced random polygons (Smith §7.7 displacement attack) ----
  std::cout << "=== Fuzz: random polygons displaced near 2^k boundaries ==="
            << std::endl;
  int dispPass = 0, dispFail = 0;
  int singlePassIdem = 0, iterConverged = 0, iterCycled = 0, iterMaxedOut = 0;
  std::map<int, int> iterDistribution;  // iter count -> count
  for (int kPow : {10, 20, 30, 40, 49}) {
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    for (int n : {8, 20, 50}) {
      for (uint64_t seed = 0; seed < 30; ++seed) {
        auto [v, e] = RandomTopologicalPolygon(n, seed + 1000 * kPow);
        Displace(&v, offset);
        // Single-pass: topological validity only.
        auto r = RemoveOverlaps2D(v, e, eps);
        bool ok =
            CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
        if (ok) {
          ++dispPass;
        } else {
          ++dispFail;
          std::cout << "  FAIL: kPow=" << kPow << " n=" << n << " seed=" << seed
                    << std::endl;
        }
        // Single-pass idempotence: re-run and compare fingerprint.
        std::vector<EdgeM> rEdges;
        for (const auto& oe : r.edges)
          rEdges.push_back({oe.v0, oe.v1, oe.mult});
        auto r2 = RemoveOverlaps2D(r.verts, rEdges, eps);
        if (Fingerprint(r, eps) == Fingerprint(r2, eps)) ++singlePassIdem;
        // Iterate-to-fixed-point.
        int iters = 0;
        IterStatus status = IterStatus::Converged;
        auto rIter =
            IterateToFixedPoint(v, e, eps, /*maxIter=*/8, &iters, &status);
        ++iterDistribution[iters];
        switch (status) {
          case IterStatus::Converged:
            ++iterConverged;
            break;
          case IterStatus::Cycled:
            ++iterCycled;
            break;
          case IterStatus::MaxedOut:
            ++iterMaxedOut;
            break;
        }
        // Topological validity on the iterated result. NOTE: rIter.inputRemap
        // is the LAST iteration's input->output remap, not the composed
        // original->final. Works for the random-topology displaced inputs
        // here because step 1 produces an identity remap when no input verts
        // merge (true for these inputs); a production version of
        // IterateToFixedPoint should compose remaps across iterations.
        if (!CheckTopologicalValidity(rIter, e, rIter.inputRemap,
                                      rIter.numMergedVerts)) {
          std::cout << "  ITER FAIL: kPow=" << kPow << " n=" << n
                    << " seed=" << seed << std::endl;
        }
      }
    }
  }
  std::cout << "  Single-pass idempotence: " << singlePassIdem << "/450"
            << std::endl;
  std::cout << "  Iterate-to-fixed-point: " << iterConverged << " converged, "
            << iterCycled << " cycled, " << iterMaxedOut << " hit max-iter (8)"
            << std::endl;
  std::cout << "  Iter distribution: ";
  for (const auto& [k, v] : iterDistribution) {
    std::cout << k << ":" << v << " ";
  }
  std::cout << std::endl;
  std::cout << "  Displacement invariant: " << dispPass << " pass, " << dispFail
            << " fail (5 scales × 3 sizes × 30 seeds = 450 cases)" << std::endl;
  if (dispFail > 0) allPass = false;
  std::cout << std::endl;

  std::cout << "==== OVERALL: " << (allPass ? "PASS" : "FAIL")
            << " ====" << std::endl;
  return allPass ? 0 : 1;
}
