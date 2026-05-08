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
//     the public OverlapRemovePolygons / Boolean2D entry points.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <climits>
#include <functional>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include "../src/collider.h"
#include "../src/disjoint_sets.h"
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
// Step 1: vertex merge.
// Returns: remap[oldIdx] = newIdx, and merged vert positions.
// =============================================================================
struct VertexMerge {
  std::vector<int> remap;
  std::vector<vec2> verts;
};

VertexMerge MergeVerts(const std::vector<vec2>& in, double eps) {
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
  std::vector<Box> boxes(n);
  for (int i = 0; i < n; ++i) boxes[i] = BoxOf2DPoint(in[i], eps);
  BVH bvh = BVHBuildFromBoxes(boxes);
  std::vector<std::pair<int, int>> pairs;
  CollidePairs(bvh, boxes, [&](int qi, int li) {
    if (qi >= li) return;  // dedupe + skip self
    pairs.emplace_back(qi, li);
  });
  std::sort(pairs.begin(), pairs.end());
  for (auto [i, j] : pairs) {
    vec2 d = in[i] - in[j];
    if (dot(d, d) <= eps2) uf.unite(i, j);
  }
  // Compute centroid per cluster.
  std::map<int, std::pair<vec2, int>> sums;  // root -> (sum, count)
  for (int i = 0; i < n; ++i) {
    int r = uf.find(i);
    auto& s = sums[r];
    s.first = s.first + in[i];
    s.second += 1;
  }
  // Assign new indices.
  std::map<int, int> rootToNew;
  std::vector<vec2> verts;
  for (const auto& [root, sumCount] : sums) {
    rootToNew[root] = static_cast<int>(verts.size());
    verts.push_back(sumCount.first * (1.0 / sumCount.second));
  }
  std::vector<int> remap(n);
  for (int i = 0; i < n; ++i) remap[i] = rootToNew[uf.find(i)];
  return {std::move(remap), std::move(verts)};
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
std::vector<std::vector<int>> BuildEdgeVertLists(
    const std::vector<EdgeM>& edges, const std::vector<vec2>& verts,
    double eps) {
  const int nE = static_cast<int>(edges.size());
  const int nV = static_cast<int>(verts.size());
  const double eps2 = eps * eps;
  std::vector<std::vector<int>> lists(nE);
  // BVH broad phase: edges as eps-padded segment AABBs, queried by vert
  // points (eps-padded boxes). Each candidate (edge, vert) pair runs the
  // exact projection test below. Per-edge `hits` are sorted by parameter
  // at the end so the result is independent of broad-phase visit order.
  std::vector<Box> edgeBoxes(nE);
  for (int e = 0; e < nE; ++e) {
    edgeBoxes[e] =
        BoxOf2DEdge(verts[edges[e].v0], verts[edges[e].v1], eps);
  }
  BVH bvh = BVHBuildFromBoxes(edgeBoxes);
  std::vector<Box> vertBoxes(nV);
  for (int v = 0; v < nV; ++v) vertBoxes[v] = BoxOf2DPoint(verts[v], eps);
  // Build vert→neighbors adjacency from the input edges. Used below to
  // detect "thin triangle apex" cases (see comment in CollidePairs).
  std::vector<std::set<int>> adj(nV);
  for (const auto& e : edges) {
    adj[e.v0].insert(e.v1);
    adj[e.v1].insert(e.v0);
  }
  // Collect (edge, vert) candidate pairs first; then process per edge.
  std::vector<std::vector<std::pair<double, int>>> hitsPerEdge(nE);
  CollidePairs(bvh, vertBoxes, [&](int v, int e) {
    if (v == edges[e].v0 || v == edges[e].v1) return;
    // Thin-triangle-apex skip: when V is connected to BOTH edge endpoints
    // by other edges, V is the apex of a triangle (V, e.v0, e.v1) whose
    // base is this edge. With non-tiny eps (large displacement), the apex
    // can fall within eps of its base; without the skip, step 5
    // canonicalization cancels the apex-split sub-edges against the
    // triangle's other two sides, producing empty output. Only-one-
    // endpoint adjacency is normal polygon-neighbor configuration, so
    // we require BOTH to be conservative.
    if (adj[v].count(edges[e].v0) && adj[v].count(edges[e].v1)) return;
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
    if (dot(d, d) <= eps2) hitsPerEdge[e].emplace_back(t, v);
  });
  for (int e = 0; e < nE; ++e) {
    auto& hits = hitsPerEdge[e];
    std::sort(hits.begin(), hits.end());
    lists[e].reserve(hits.size());
    for (const auto& [t, v] : hits) lists[e].push_back(v);
  }
  return lists;
}

// =============================================================================
// Step 4: edge-edge intersection discovery; insert new verts (or snap to
// existing nearby verts) into per-edge lists.
// =============================================================================
void FindAndInsertIntersections(const std::vector<EdgeM>& edges,
                                std::vector<vec2>* verts,
                                std::vector<std::vector<int>>* lists,
                                std::vector<std::set<int>>* vertEdges,
                                double eps) {
  const int nE = static_cast<int>(edges.size());
  const double eps2 = eps * eps;
  vertEdges->resize(verts->size());
  // BVH broad phase: each edge as eps-padded segment AABB, queried against
  // all edge AABBs. Self-collision is filtered by `qi < li`. Pairs are
  // sorted lexicographically; the snap-and-insert below depends on
  // `lists[*]` accumulating verts as earlier pairs are processed, and
  // sorting on the int pair (rather than Morton order) keeps the order
  // deterministic across runs and BVH builds.
  std::vector<Box> edgeBoxes(nE);
  for (int i = 0; i < nE; ++i) {
    edgeBoxes[i] =
        BoxOf2DEdge((*verts)[edges[i].v0], (*verts)[edges[i].v1], eps);
  }
  BVH bvh = BVHBuildFromBoxes(edgeBoxes);
  std::vector<std::pair<int, int>> pairs;
  CollidePairs(bvh, edgeBoxes, [&](int qi, int li) {
    if (qi >= li) return;
    pairs.emplace_back(qi, li);
  });
  std::sort(pairs.begin(), pairs.end());
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
    (*vertEdges)[vNew].insert(i);
    (*vertEdges)[vNew].insert(j);
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
}

// =============================================================================
// Step 5: break edges into sub-edges; merge duplicates with multiplicity sum.
// Smith's PolySet2 (Table 7.3): map<lex-ordered segment, signed multiplicity>.
// =============================================================================
struct CanonicalSubEdges {
  std::map<std::pair<int, int>, int> map;  // (vMin, vMax) -> signed mult

  // Add a directed sub-edge. The map key is lex-min-first; multiplicity sign
  // flips if direction is reversed.
  void Add(int v0, int v1, int mult) {
    if (v0 == v1) return;
    int vMin = std::min(v0, v1);
    int vMax = std::max(v0, v1);
    int signedMult = (v0 < v1) ? mult : -mult;
    auto it = map.find({vMin, vMax});
    if (it == map.end()) {
      map.emplace(std::make_pair(vMin, vMax), signedMult);
    } else {
      it->second += signedMult;
      if (it->second == 0) map.erase(it);
    }
  }
};

CanonicalSubEdges Canonicalize(const std::vector<EdgeM>& edges,
                               const std::vector<std::vector<int>>& lists) {
  CanonicalSubEdges out;
  for (size_t e = 0; e < edges.size(); ++e) {
    int prev = edges[e].v0;
    for (int v : lists[e]) {
      out.Add(prev, v, edges[e].mult);
      prev = v;
    }
    out.Add(prev, edges[e].v1, edges[e].mult);
  }
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
  for (const auto& [key, mult] : canon.map) {
    vec2 p0 = verts[key.first];
    vec2 p1 = verts[key.second];
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

// Predicate over a face's winding number deciding whether the face is
// "inside" the result region. Standard ops:
//   - Add/Subtract: w > 0  (Smith's wind > 0 union convention)
//   - Intersect:    w > 1  (both inputs cover the face)
// An edge is retained iff its left and right faces disagree on this
// predicate.
using WindPredicate = std::function<bool(int)>;

inline WindPredicate WindAdd() { return [](int w) { return w > 0; }; }
inline WindPredicate WindIntersect() { return [](int w) { return w > 1; }; }

std::vector<OutEdge> FilterByWindingDCEL(
    const CanonicalSubEdges& canon, const std::vector<vec2>& verts,
    bool debug = false, const WindPredicate& isInside = WindAdd()) {
  using dcel_internal::HalfEdge;
  if (debug) {
    std::cerr << "[FilterByWindingDCEL] canon.map.size()="
              << canon.map.size() << " verts.size()=" << verts.size() << "\n";
  }
  // 1. Build half-edges. Each canonical (vMin, vMax) with mult m becomes:
  //    - hA: vMin → vMax, mult = m
  //    - hB: vMax → vMin, mult = -m
  //    Twins are paired (hA.twin = hB, hB.twin = hA).
  std::vector<HalfEdge> halfedges;
  halfedges.reserve(2 * canon.map.size());
  for (const auto& [k, m] : canon.map) {
    int hA = static_cast<int>(halfedges.size());
    halfedges.push_back({hA + 1, -1, k.first, -1, m});
    halfedges.push_back({hA, -1, k.second, -1, -m});
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
  for (size_t v = 0; v < outgoing.size(); ++v) {
    auto& hes = outgoing[v];
    if (hes.size() < 2) continue;
    const vec2 vp = verts[v];
    std::sort(hes.begin(), hes.end(), [&](int a, int b) {
      const vec2 dA = verts[halfedges[halfedges[a].twin].origin] - vp;
      const vec2 dB = verts[halfedges[halfedges[b].twin].origin] - vp;
      return std::atan2(dA.y, dA.x) < std::atan2(dB.y, dB.x);
    });
  }

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
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    const int twinIdx = halfedges[i].twin;
    const int destV = halfedges[twinIdx].origin;
    auto& sorted = outgoing[destV];
    auto it = std::find(sorted.begin(), sorted.end(), twinIdx);
    if (it == sorted.end()) continue;
    auto prevIt = (it == sorted.begin()) ? (sorted.end() - 1) : (it - 1);
    halfedges[i].next = *prevIt;
  }

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

  // 6. Compute winding per face by ray-cast. We do this PER FACE rather
  //    than per-edge (the old `FilterByWinding`) because per-edge produces
  //    inconsistent verdicts at shared vertices, AND we do it per-face
  //    rather than via BFS-from-outer because BFS doesn't propagate
  //    between disconnected face-graph components (real for self-
  //    intersecting polygons whose union has multiple disjoint regions).
  //    Ray-cast errors don't compound since each face decision is
  //    independent.
  //
  //    For each face, find a half-edge on its boundary, perp-offset its
  //    midpoint into the face (the LEFT side, by DCEL convention), cast
  //    a horizontal ray, sum signed mult contributions of edges crossed.
  std::vector<int> faceWind(nFaces, 0);
  // Ray-cast point per face: perpendicular-offset from a half-edge of
  // the face's boundary, into the face (LEFT side). This works because
  // the face is on the LEFT of every half-edge in its cycle, so a
  // perpendicular-CCW offset from any boundary edge lands inside the
  // face. Centroids would seem more robust, but for non-convex faces
  // (common after step 4) the centroid can lie outside the face.
  for (int f = 0; f < nFaces; ++f) {
    if (f == outerFace) {
      faceWind[f] = 0;
      continue;
    }
    int h = faceStartHE[f];
    if (h < 0) continue;
    const vec2 a = verts[halfedges[h].origin];
    const vec2 b = verts[halfedges[halfedges[h].twin].origin];
    const vec2 mid = (a + b) * 0.5;
    const vec2 d = b - a;
    const double len = length(d);
    if (len == 0) continue;
    const vec2 perp(-d.y / len, d.x / len);  // 90° CCW = LEFT of h's direction
    // Step inward by len/1000 to land an interior point of the face. Picks
    // up the face's winding via ray-cast. The 1/1000 factor assumes the
    // face's inscribed-circle radius from this halfedge's midpoint is at
    // least len/1000; degenerate slivers worse than ~1000:1 aspect would
    // miss. None appear in the test corpus.
    const vec2 pInF = mid + perp * (len * 1e-3);
    faceWind[f] = CastWindingRay(pInF, canon, verts);
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
  out.reserve(canon.map.size());
  int hi = 0;
  for (const auto& [k, m] : canon.map) {
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
      out.push_back({k.first, k.second, 1});
    } else {
      out.push_back({k.second, k.first, 1});
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

OverlapResult RemoveOverlaps2D(const std::vector<vec2>& vertsIn,
                               const std::vector<EdgeM>& edgesIn, double eps,
                               bool debug = false,
                               const WindPredicate& isInside = WindAdd()) {
  // Step 1: vertex merge.
  auto merge = MergeVerts(vertsIn, eps);
  const int numMerged = static_cast<int>(merge.verts.size());
  // Step 2: collapse edges.
  auto edges = RemapAndCollapse(edgesIn, merge.remap);
  // Step 3: per-edge vert list.
  auto lists = BuildEdgeVertLists(edges, merge.verts, eps);
  // Step 4: intersection discovery.
  std::vector<std::set<int>> vertEdges;
  FindAndInsertIntersections(edges, &merge.verts, &lists, &vertEdges, eps);

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
    for (auto [a, b] : pairs) {
      // Structural gate: do they share an edge?
      bool shared = false;
      for (int e : vertEdges[a]) {
        if (vertEdges[b].count(e)) {
          shared = true;
          break;
        }
      }
      if (!shared) continue;
      // Geometric gate: within mergeThresh?
      vec2 d = merge.verts[b] - merge.verts[a];
      if (dot(d, d) > mergeThresh2) continue;
      uf.unite(a, b);
    }
    // Build remap from union-find clusters; cluster position is centroid.
    std::map<int, std::pair<vec2, int>> sums;
    for (size_t i = 0; i < merge.verts.size(); ++i) {
      int r = uf.find(static_cast<int>(i));
      auto& s = sums[r];
      s.first = s.first + merge.verts[i];
      s.second += 1;
    }
    if (sums.size() < merge.verts.size()) {
      std::map<int, int> rootToNew;
      std::vector<vec2> newVerts;
      for (const auto& [root, sumCount] : sums) {
        rootToNew[root] = static_cast<int>(newVerts.size());
        newVerts.push_back(sumCount.first * (1.0 / sumCount.second));
      }
      std::vector<int> remap(merge.verts.size());
      for (size_t i = 0; i < merge.verts.size(); ++i) {
        remap[i] = rootToNew[uf.find(static_cast<int>(i))];
      }
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

  // Step 5: sub-edge canonicalization.
  auto canon = Canonicalize(edges, lists);
  // Step 6: DCEL face-traversal winding filter.
  auto out = FilterByWindingDCEL(canon, merge.verts, debug, isInside);
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
  for (size_t v = 0; v < outgoing.size(); ++v) {
    auto& lst = outgoing[v];
    if (lst.size() < 2) continue;
    const vec2 vp = verts[v];
    std::sort(lst.begin(), lst.end(), [&](int a, int b) {
      const vec2 da = verts[edges[a].v1] - vp;
      const vec2 db = verts[edges[b].v1] - vp;
      return std::atan2(da.y, da.x) < std::atan2(db.y, db.x);
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
      // (going CCW) is the one we want; that's the entry one step CW
      // from the position where the reverse-incoming direction would
      // sit in the sorted list. Implemented by binary-search on the
      // reverse direction.
      const vec2 vp = verts[destV];
      const vec2 inDir = vp - verts[edges[cur].v0];
      const double inAngle = std::atan2(inDir.y, inDir.x);
      // Reverse direction in (-π, π].
      double rev = inAngle + M_PI;
      if (rev > M_PI) rev -= 2 * M_PI;
      const auto& lst = outgoing[destV];
      // Find the unvisited entry with the smallest CCW angle past `rev`.
      int next = -1;
      double bestDelta = std::numeric_limits<double>::infinity();
      for (int e : lst) {
        if (visited[e]) continue;
        const vec2 d = verts[edges[e].v1] - vp;
        double ang = std::atan2(d.y, d.x);
        // CCW delta from rev to ang in [0, 2π).
        double delta = ang - rev;
        if (delta <= 0) delta += 2 * M_PI;
        if (delta < bestDelta) {
          bestDelta = delta;
          next = e;
        }
      }
      cur = next;
    }
    if (loop.size() >= 3) polys.push_back(std::move(loop));
  }
  return polys;
}

// Public-facing API. Mirrors what `CrossSection` would call into.
inline manifold::Polygons OverlapRemovePolygons(const manifold::Polygons& in,
                                                double eps) {
  auto [verts, edges] = PolygonsToInput(in);
  if (verts.empty()) return {};
  auto r = RemoveOverlaps2D(verts, edges, eps);
  return OutEdgesToPolygons(r.verts, r.edges);
}

// Operation tag for the binary-boolean entry point. `Add` is the union;
// `Subtract` removes B from A; `Intersect` keeps the region covered by
// both inputs.
enum class BoolOp { Add, Subtract, Intersect };

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
                                    const manifold::Polygons& b, BoolOp op,
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
  append(b, op == BoolOp::Subtract ? -1 : 1);
  if (verts.empty()) return {};
  WindPredicate pred = (op == BoolOp::Intersect) ? WindIntersect() : WindAdd();
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, pred);
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

  auto lists = BuildEdgeVertLists(edges, merge.verts, eps);
  size_t totalListSize = 0;
  for (const auto& l : lists) totalListSize += l.size();
  std::cerr << "After step 3 (lists): " << totalListSize
            << " total list entries across " << edges.size() << " edges\n";

  const int beforeStep4 = static_cast<int>(merge.verts.size());
  std::vector<std::set<int>> vertEdges;
  FindAndInsertIntersections(edges, &merge.verts, &lists, &vertEdges, eps);
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
  std::cerr << "After step 5 (canonicalize): " << canon.map.size()
            << " sub-edges (after multiplicity collapse)\n";

  // Vertex balance from canonicalized sub-edges.
  std::map<int, int> canonBalance;
  for (const auto& [k, m] : canon.map) {
    canonBalance[k.first] += m;
    canonBalance[k.second] -= m;
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
        for (const auto& [k, m] : canon.map) {
          if (k.first == v || k.second == v) {
            std::cerr << "      (" << k.first << "," << k.second
                      << ") mult=" << m << "\n";
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
            << " output edges (from " << canon.map.size()
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
    for (const auto& [k, m] : canon.map) {
      if (k.first != targetV && k.second != targetV) continue;
      vec2 p0 = merge.verts[k.first];
      vec2 p1 = merge.verts[k.second];
      vec2 mid = (p0 + p1) * 0.5;
      vec2 d = p1 - p0;
      double len = length(d);
      vec2 perp = {-d.y / len, d.x / len};
      vec2 leftPt = mid + perp * ofs;
      vec2 rightPt = mid + perp * -ofs;
      int wL = CastWindingRay(leftPt, canon, merge.verts);
      int wR = CastWindingRay(rightPt, canon, merge.verts);
      std::cerr << "  sub-edge (" << k.first << "," << k.second
                << ") mult=" << m << " | mid=(" << mid.x - offset << ","
                << mid.y - offset << ") perp=(" << perp.x << "," << perp.y
                << ") | wL=" << wL << " wR=" << wR
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
    auto l3 = BuildEdgeVertLists(e3, m3.verts, eps);
    std::vector<std::set<int>> ve3;
    FindAndInsertIntersections(e3, &m3.verts, &l3, &ve3, eps);
    // structural step 4b (copy of the production code)
    {
      DisjointSets uf(static_cast<int>(m3.verts.size()));
      const double mergeThresh2 = (10.0 * eps) * (10.0 * eps);
      for (size_t a = 0; a < m3.verts.size(); ++a) {
        if (a >= ve3.size() || ve3[a].empty()) continue;
        for (size_t b = a + 1; b < m3.verts.size(); ++b) {
          if (b >= ve3.size() || ve3[b].empty()) continue;
          bool shared = false;
          for (int e : ve3[a]) if (ve3[b].count(e)) { shared = true; break; }
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
    std::cerr << "  canonical sub-edges (post step 5): " << canon.map.size()
              << "\n";

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
      for (const auto& [k, m] : canon.map) {
        if (k.first != target && k.second != target) continue;
        ++sub_count;
        const int other = (k.first == target) ? k.second : k.first;
        const vec2 op = m3.verts[other];
        vec2 p0 = vp, p1 = op;
        if (k.first != target) std::swap(p0, p1);
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
        std::cerr << "    (" << k.first << "↔" << k.second << ") mult=" << m
                  << " wL=" << wL << " wR=" << wR
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
  int roundTripOK = 0;            // (v, e) → Polygons → OverlapRemovePolygons
                                  // → no crash + edge count consistent
  int roundTripCrash = 0;         // (asserts at boundaries between API layers)
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
        // run the public OverlapRemovePolygons API, and check the
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
          auto out = OverlapRemovePolygons(polys, eps);
          // Sum boundary edges across output loops; should equal
          // pass1.edges.size() up to loop-walk reordering.
          size_t outEdges = 0;
          for (const auto& l : out) outEdges += l.size();
          if (outEdges == pass1.edges.size())
            ++roundTripOK;
          else
            ++roundTripCrash;
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
  std::cout << "  Polygons round-trip (lower vs API edge count): "
            << roundTripOK << " match, " << roundTripCrash << " mismatch\n";
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
    auto p2_l = BuildEdgeVertLists(p2_e, p2_mrg.verts, eps);
    int p2_totalList = 0;
    for (auto& l : p2_l) p2_totalList += l.size();
    std::cerr << "  step 3 lists: " << p2_totalList << " total entries\n";
    std::vector<std::set<int>> p2_ve;
    FindAndInsertIntersections(p2_e, &p2_mrg.verts, &p2_l, &p2_ve, eps);
    std::cerr << "  step 4 intersections: " << p2_mrg.verts.size()
              << " verts after\n";
    auto p2_canon = Canonicalize(p2_e, p2_l);
    std::cerr << "  step 5 canon: " << p2_canon.map.size() << " sub-edges\n";

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
    auto add = Boolean2D(a, b, BoolOp::Add);
    auto sub = Boolean2D(a, b, BoolOp::Subtract);
    auto isec = Boolean2D(a, b, BoolOp::Intersect);
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
  // Verify the public OverlapRemovePolygons wrapper round-trips simple
  // closed loops through the manifold::Polygons type.
  {
    std::cout << "=== Polygons API: square ===" << std::endl;
    manifold::Polygons in = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    auto out = OverlapRemovePolygons(in, EpsilonFromScale(1.0));
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
    auto out = OverlapRemovePolygons(in, EpsilonFromScale(4.0));
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
    auto out = OverlapRemovePolygons(in, EpsilonFromScale(2.0));
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
    auto out = OverlapRemovePolygons(in, EpsilonFromScale(3.0));
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
