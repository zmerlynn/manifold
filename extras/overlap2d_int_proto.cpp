// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Integer prototype for 2D overlap removal (issue #289).
//
// Sibling to extras/overlap2d_proto.cpp (the FP prototype). Same 6-step
// architecture (BVH broad-phase / DCEL face arrangement / winding-number
// filter), but the geometric kernel is exact int128 instead of
// FP+symbolic-perturbation. Goal: validate the hypothesis that an
// integer kernel closes most of the perf gap to Clipper64 without
// changing the high-level shape.
//
// Coordinate convention: int64_t coordinates. Caller scales doubles to
// int (factor SCALE configurable; default 1e9 — gives sub-nanometer
// precision for UTM-million inputs and fits comfortably in int64).
//
// Predicate strategy:
//   - Orient (sign of cross product) — exact via __int128. Drives all
//     topological decisions.
//   - Intersection point — double-FP computation, rounded to int grid.
//     Exact-int parametric coordinates would need >128 bits in the
//     general case (numT * Δx with each term ~62 bits = 186 bits).
//     Snap-rounding to int means our final answer is "the closest int
//     point" anyway, which double precision delivers to ~1 ULP at the
//     int grid. Topology stays correct because orient is exact.
//
// Build (single-TU prototype, no MANIFOLD_PAR dependency for the spike):
//   g++ -std=c++17 -O2 -I include -I src -DMANIFOLD_PAR=-1 \
//     -ffp-contract=off -fexcess-precision=standard \
//     extras/overlap2d_int_proto.cpp -o overlap2d_int_proto
//
// Status: spike. Only `Simplify` (single-input regularization) is
// implemented. Multi-input booleans, parallelism, and the corpus
// runners are deferred until the perf hypothesis is validated.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "../src/collider.h"
#include "manifold/common.h"

namespace overlap2d_int {

using i64 = int64_t;
using i128 = __int128;
using manifold::Box;
using manifold::Collider;
using manifold::MakeSimpleRecorder;
using manifold::Polygons;
using manifold::SimplePolygon;
using manifold::vec3;

struct Vec2i {
  i64 x, y;
  bool operator==(Vec2i o) const { return x == o.x && y == o.y; }
  bool operator<(Vec2i o) const {
    return x != o.x ? x < o.x : y < o.y;
  }
};

// ============================================================================
// Exact predicates (i128).
// ============================================================================

// Sign of (b - a) × (c - a). Returns +1 (CCW), -1 (CW), 0 (collinear).
// Inputs must satisfy max(|coord|) < 2^62 for the i128 product to not
// overflow (|Δ| ≤ 2^63, product ≤ 2^126 fits in i128).
inline int Orient(Vec2i a, Vec2i b, Vec2i c) {
  const i128 cross =
      (i128)(b.x - a.x) * (i128)(c.y - a.y) -
      (i128)(b.y - a.y) * (i128)(c.x - a.x);
  return (cross > 0) - (cross < 0);
}

// Do segments AB and CD properly cross (interiors meet at a single point)?
// Strict crossing only — endpoint coincidence returns false.
inline bool SegmentsProperlyCross(Vec2i a, Vec2i b, Vec2i c, Vec2i d) {
  const int o1 = Orient(a, b, c);
  const int o2 = Orient(a, b, d);
  const int o3 = Orient(c, d, a);
  const int o4 = Orient(c, d, b);
  // Both pairs straddle the other segment's line, with no collinearity.
  return o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4;
}

// Is point p strictly between (not equal to) endpoints a and b on the line ab?
// Caller must have already established collinearity (Orient(a, b, p) == 0).
inline bool BetweenCollinear(Vec2i a, Vec2i b, Vec2i p) {
  if (p == a || p == b) return false;
  // Use the longer axis to avoid degenerate axis projection.
  if (a.x != b.x) {
    return (a.x < p.x && p.x < b.x) || (b.x < p.x && p.x < a.x);
  }
  return (a.y < p.y && p.y < b.y) || (b.y < p.y && p.y < a.y);
}

// Squared distance between two int verts, as i128 to avoid overflow on
// large-magnitude coords. Used for the eps-merge gate (which becomes
// "snap to int grid" in the integer formulation: eps² = 0).
inline i128 Dist2(Vec2i a, Vec2i b) {
  const i128 dx = (i128)(a.x - b.x);
  const i128 dy = (i128)(a.y - b.y);
  return dx * dx + dy * dy;
}

// Compute intersection point of AB and CD as integers (rounded to grid).
// Returns true if a finite intersection exists. Caller should check
// SegmentsProperlyCross first to know the result is meaningful.
//
// Uses double for the parametric divide; orient was already exact so
// topology is correct, and the intersection POINT is being snap-rounded
// to int anyway.
inline bool IntersectPointSnapped(Vec2i a, Vec2i b, Vec2i c, Vec2i d,
                                  Vec2i* out) {
  const double denom =
      (double)(b.x - a.x) * (double)(d.y - c.y) -
      (double)(b.y - a.y) * (double)(d.x - c.x);
  if (denom == 0) return false;
  const double t =
      ((double)(c.x - a.x) * (double)(d.y - c.y) -
       (double)(c.y - a.y) * (double)(d.x - c.x)) /
      denom;
  const double px = (double)a.x + t * (double)(b.x - a.x);
  const double py = (double)a.y + t * (double)(b.y - a.y);
  out->x = (i64)std::llround(px);
  out->y = (i64)std::llround(py);
  return true;
}

// ============================================================================
// BVH broad-phase (reuses manifold::Collider, which expects FP boxes).
// Int verts are widened to double for the box; broad-phase doesn't need
// exact arithmetic, only the narrow-phase predicates do (and those stay
// in i128). Padding by 1 unit on each side handles snap-to-grid boundary
// cases without missing any candidate.
// ============================================================================

inline Box BoxOf2DPointI(Vec2i p) {
  vec3 lo((double)p.x - 1, (double)p.y - 1, 0);
  vec3 hi((double)p.x + 1, (double)p.y + 1, 0);
  return Box(lo, hi);
}

inline Box BoxOf2DEdgeI(Vec2i p0, Vec2i p1) {
  Box b(vec3((double)p0.x, (double)p0.y, 0),
        vec3((double)p1.x, (double)p1.y, 0));
  vec3 pad(1, 1, 0);
  return Box(b.min - pad, b.max + pad);
}

struct BVH {
  Collider collider;
  std::vector<int> leafToOrig;  // sortedLeafIdx -> caller's input index
};

inline BVH BVHBuildFromBoxes(const std::vector<Box>& boxes) {
  const int n = (int)boxes.size();
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
  out.collider = Collider(
      manifold::VecView<const Box>(sortedBB.data(), sortedBB.size()),
      manifold::VecView<const uint32_t>(sortedMorton.data(),
                                         sortedMorton.size()));
  return out;
}

template <typename F>
inline void CollidePairs(const BVH& bvh, const std::vector<Box>& queries,
                         F&& f) {
  if (bvh.leafToOrig.empty() || queries.empty()) return;
  auto adapter = [&](int qi, int leafIdx) { f(qi, bvh.leafToOrig[leafIdx]); };
  auto recorder = MakeSimpleRecorder(adapter);
  auto qf = [&](int i) { return queries[i]; };
  bvh.collider.Collisions<false>(recorder, qf, (int)queries.size(),
                                 /*parallel=*/false);
}

// ============================================================================
// Step 0: convert manifold::Polygons (FP) ↔ integer representation.
// Scale factor SCALE multiplies input doubles before snapping to int.
// ============================================================================

constexpr i64 kDefaultScale = 1'000'000'000;  // 1e9: nanometer precision @ unit-1.

struct EdgeM {
  int v0, v1;
  int mult;
};

inline std::pair<std::vector<Vec2i>, std::vector<EdgeM>> PolygonsToInput(
    const Polygons& in, i64 scale) {
  std::vector<Vec2i> verts;
  std::vector<EdgeM> edges;
  for (const auto& loop : in) {
    if (loop.size() < 3) continue;
    const int base = (int)verts.size();
    for (const auto& v : loop) {
      verts.push_back({(i64)std::llround(v.x * scale),
                       (i64)std::llround(v.y * scale)});
    }
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i) {
      edges.push_back({base + i, base + ((i + 1) % n), 1});
    }
  }
  return {verts, edges};
}

// ============================================================================
// Step 1: vertex merge. With integer coords and a snap-to-grid model, "merge"
// is "coalesce verts at identical int positions" — eps == 0 effectively.
// O(n log n) via sort + sweep on the (x, y) key.
// ============================================================================

struct VertMergeResult {
  std::vector<int> remap;   // input idx → output idx
  std::vector<Vec2i> verts; // deduped output positions
};

inline VertMergeResult MergeVerts(const std::vector<Vec2i>& in) {
  const int n = (int)in.size();
  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](int a, int b) { return in[a] < in[b]; });
  VertMergeResult r;
  r.remap.assign(n, -1);
  for (int k = 0; k < n;) {
    int j = k;
    while (j < n && in[order[j]] == in[order[k]]) ++j;
    const int newIdx = (int)r.verts.size();
    r.verts.push_back(in[order[k]]);
    for (int m = k; m < j; ++m) r.remap[order[m]] = newIdx;
    k = j;
  }
  return r;
}

// ============================================================================
// Step 2: collapse edges whose endpoints map to the same merged vert.
// ============================================================================

inline std::vector<EdgeM> RemapAndCollapse(const std::vector<EdgeM>& in,
                                            const std::vector<int>& remap) {
  std::vector<EdgeM> out;
  out.reserve(in.size());
  for (const auto& e : in) {
    const int a = remap[e.v0];
    const int b = remap[e.v1];
    if (a != b) out.push_back({a, b, e.mult});
  }
  return out;
}

// ============================================================================
// Step 3: per-edge ordered list of verts that lie on the edge interior.
// O(n_edges · n_verts) brute-force for the spike — replace with BVH later.
// ============================================================================

inline std::vector<std::vector<int>> BuildEdgeVertLists(
    const std::vector<EdgeM>& edges, const std::vector<Vec2i>& verts) {
  const int nE = (int)edges.size();
  const int nV = (int)verts.size();
  std::vector<std::vector<int>> lists(nE);

  // BVH on edges, query each vert against it. Saves O(E·V) → O((E+V) log E).
  std::vector<Box> edgeBoxes(nE);
  for (int e = 0; e < nE; ++e)
    edgeBoxes[e] = BoxOf2DEdgeI(verts[edges[e].v0], verts[edges[e].v1]);
  BVH bvh = BVHBuildFromBoxes(edgeBoxes);
  std::vector<Box> vertBoxes(nV);
  for (int v = 0; v < nV; ++v) vertBoxes[v] = BoxOf2DPointI(verts[v]);

  std::vector<std::vector<std::pair<i128, int>>> hitsPerEdge(nE);
  CollidePairs(bvh, vertBoxes, [&](int v, int e) {
    if (v == edges[e].v0 || v == edges[e].v1) return;
    const Vec2i a = verts[edges[e].v0];
    const Vec2i b = verts[edges[e].v1];
    const Vec2i p = verts[v];
    if (Orient(a, b, p) != 0) return;
    if (!BetweenCollinear(a, b, p)) return;
    const i128 key = (i128)(p.x - a.x) * (b.x - a.x) +
                     (i128)(p.y - a.y) * (b.y - a.y);
    hitsPerEdge[e].emplace_back(key, v);
  });
  for (int e = 0; e < nE; ++e) {
    auto& hits = hitsPerEdge[e];
    std::sort(hits.begin(), hits.end());
    lists[e].reserve(hits.size());
    for (const auto& [_, v] : hits) lists[e].push_back(v);
  }
  return lists;
}

// ============================================================================
// Step 4: edge-edge intersection discovery. For each crossing pair, compute
// the int-snapped intersection, dedupe against existing verts (eps = 0 here:
// equality is exact), and insert into both edges' on-edge lists.
// O(n²) brute-force; BVH replacement later.
// ============================================================================

inline void FindAndInsertIntersections(
    const std::vector<EdgeM>& edges, std::vector<Vec2i>* verts,
    std::vector<std::vector<int>>* lists) {
  const int nE = (int)edges.size();
  // Map from int position → vert index for snap-dedup.
  std::map<Vec2i, int> posToIdx;
  for (int i = 0; i < (int)verts->size(); ++i) posToIdx[(*verts)[i]] = i;

  // BVH on edges, self-collide. Saves O(E²) → O((E + k) log E).
  std::vector<Box> edgeBoxes(nE);
  for (int e = 0; e < nE; ++e)
    edgeBoxes[e] = BoxOf2DEdgeI((*verts)[edges[e].v0], (*verts)[edges[e].v1]);
  BVH bvh = BVHBuildFromBoxes(edgeBoxes);
  std::vector<std::pair<int, int>> pairs;
  CollidePairs(bvh, edgeBoxes, [&](int qi, int li) {
    if (qi >= li) return;
    pairs.emplace_back(qi, li);
  });
  std::sort(pairs.begin(), pairs.end());

  for (auto [i, j] : pairs) {
    const Vec2i a = (*verts)[edges[i].v0];
    const Vec2i b = (*verts)[edges[i].v1];
    const Vec2i c = (*verts)[edges[j].v0];
    const Vec2i d = (*verts)[edges[j].v1];
    if (edges[i].v0 == edges[j].v0 || edges[i].v0 == edges[j].v1 ||
        edges[i].v1 == edges[j].v0 || edges[i].v1 == edges[j].v1)
      continue;
    if (!SegmentsProperlyCross(a, b, c, d)) continue;
    Vec2i p;
    if (!IntersectPointSnapped(a, b, c, d, &p)) continue;
    int vNew;
    auto it = posToIdx.find(p);
    if (it != posToIdx.end()) {
      vNew = it->second;
    } else {
      vNew = (int)verts->size();
      verts->push_back(p);
      posToIdx[p] = vNew;
    }
    auto insertSorted = [&](int eIdx) {
      if (vNew == edges[eIdx].v0 || vNew == edges[eIdx].v1) return;
      const Vec2i ea = (*verts)[edges[eIdx].v0];
      const Vec2i eb = (*verts)[edges[eIdx].v1];
      const i128 keyP =
          (i128)(p.x - ea.x) * (eb.x - ea.x) +
          (i128)(p.y - ea.y) * (eb.y - ea.y);
      auto& lst = (*lists)[eIdx];
      auto pos = std::lower_bound(
          lst.begin(), lst.end(), keyP, [&](int v, i128 k) {
            const Vec2i vv = (*verts)[v];
            const i128 kv =
                (i128)(vv.x - ea.x) * (eb.x - ea.x) +
                (i128)(vv.y - ea.y) * (eb.y - ea.y);
            return kv < k;
          });
      if (pos == lst.end() || *pos != vNew) lst.insert(pos, vNew);
    };
    insertSorted(i);
    insertSorted(j);
  }
}

// ============================================================================
// Step 5: canonicalize sub-edges. Each input edge is split into runs of
// consecutive (a → list[0] → list[1] → … → b) sub-edges. Duplicate sub-edges
// (same canonical (vMin, vMax)) accumulate signed multiplicities.
// ============================================================================

struct CanonicalSubEdges {
  std::map<std::pair<int, int>, int> map;  // (vMin, vMax) → signed mult
};

inline CanonicalSubEdges Canonicalize(
    const std::vector<EdgeM>& edges,
    const std::vector<std::vector<int>>& lists) {
  CanonicalSubEdges canon;
  auto add = [&](int v0, int v1, int mult) {
    if (v0 == v1) return;
    int a = std::min(v0, v1), b = std::max(v0, v1);
    int signedMult = (v0 < v1) ? mult : -mult;
    canon.map[{a, b}] += signedMult;
  };
  for (int e = 0; e < (int)edges.size(); ++e) {
    const auto& list = lists[e];
    int prev = edges[e].v0;
    for (int v : list) {
      add(prev, v, edges[e].mult);
      prev = v;
    }
    add(prev, edges[e].v1, edges[e].mult);
  }
  // Drop zero-mult entries.
  for (auto it = canon.map.begin(); it != canon.map.end();) {
    if (it->second == 0) it = canon.map.erase(it);
    else ++it;
  }
  return canon;
}

// ============================================================================
// Step 6: DCEL face traversal + winding-number filter.
//
// For each canonical sub-edge with signed mult m, emit two half-edges
// (vMin → vMax, mult=m) and (vMax → vMin, mult=-m). For each vertex,
// sort outgoing half-edges by angle. Walk face cycles via "smallest left
// turn" next-pointer rule. For each face except the outer, ray-cast
// horizontally from a point inside to compute winding.
// ============================================================================

struct OutEdge {
  int v0, v1, mult;
};

namespace dcel {
struct HalfEdge {
  int twin, next, origin, face, mult;
};
}  // namespace dcel

// Cast +x ray from origin in int-FP arithmetic. For each canonical edge
// (vMin → vMax) with signed mult m: if the ray crosses it, contribute
// ±m depending on direction. Edges with vMin.y == origin.y are handled
// by half-open [yMin, yMax) convention.
inline int CastWindingRay(Vec2i origin, const CanonicalSubEdges& canon,
                          const std::vector<Vec2i>& verts) {
  int winding = 0;
  for (const auto& [key, mult] : canon.map) {
    Vec2i p0 = verts[key.first];
    Vec2i p1 = verts[key.second];
    bool upward = p0.y < p1.y;
    if (!upward) std::swap(p0, p1);
    if (origin.y < p0.y || origin.y >= p1.y) continue;
    // Does the ray (origin.x, origin.y) → (+inf, origin.y) cross edge?
    // Check exactly: at y=origin.y, x_cross = p0.x + (p1.x - p0.x) *
    // (origin.y - p0.y) / (p1.y - p0.y). Need x_cross > origin.x.
    // Multiply through by (p1.y - p0.y) > 0:
    //   (origin.x - p0.x) * (p1.y - p0.y) < (p1.x - p0.x) * (origin.y - p0.y)
    const i128 lhs = (i128)(origin.x - p0.x) * (i128)(p1.y - p0.y);
    const i128 rhs = (i128)(p1.x - p0.x) * (i128)(origin.y - p0.y);
    if (lhs < rhs) winding += upward ? mult : -mult;
  }
  return winding;
}

inline std::vector<OutEdge> FilterByWindingDCEL(
    const CanonicalSubEdges& canon, const std::vector<Vec2i>& verts,
    const std::function<bool(int)>& isInside) {
  using dcel::HalfEdge;
  std::vector<HalfEdge> halfedges;
  halfedges.reserve(2 * canon.map.size());
  for (const auto& [k, m] : canon.map) {
    const int hA = (int)halfedges.size();
    halfedges.push_back({hA + 1, -1, k.first, -1, m});
    halfedges.push_back({hA, -1, k.second, -1, -m});
  }
  if (halfedges.empty()) return {};

  // Group + angular-sort outgoing half-edges per vertex.
  std::vector<std::vector<int>> outgoing(verts.size());
  for (int i = 0; i < (int)halfedges.size(); ++i)
    outgoing[halfedges[i].origin].push_back(i);
  for (size_t v = 0; v < outgoing.size(); ++v) {
    auto& hes = outgoing[v];
    if (hes.size() < 2) continue;
    const Vec2i vp = verts[v];
    std::sort(hes.begin(), hes.end(), [&](int a, int b) {
      const Vec2i dA = {verts[halfedges[halfedges[a].twin].origin].x - vp.x,
                        verts[halfedges[halfedges[a].twin].origin].y - vp.y};
      const Vec2i dB = {verts[halfedges[halfedges[b].twin].origin].x - vp.x,
                        verts[halfedges[halfedges[b].twin].origin].y - vp.y};
      // Lexicographic angle bucket: half-plane (y > 0, or y == 0 && x > 0) =
      // upper. Within a half-plane, sort by orient against a reference axis.
      auto bucket = [](const Vec2i& d) {
        if (d.y > 0 || (d.y == 0 && d.x > 0)) return 0;
        return 1;
      };
      const int bA = bucket(dA), bB = bucket(dB);
      if (bA != bB) return bA < bB;
      // Same half-plane: dA before dB iff cross(dA, dB) > 0 (dA is CW of dB).
      // Wait — we want CCW order, so dA before dB iff dA is "CCW-earlier".
      // cross(dA, dB) > 0 means dB is CCW of dA, so dA comes first.
      const i128 cross = (i128)dA.x * dB.y - (i128)dA.y * dB.x;
      return cross > 0;
    });
  }

  // Compute next pointers.
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    const int twinIdx = halfedges[i].twin;
    const int destV = halfedges[twinIdx].origin;
    auto& sorted = outgoing[destV];
    auto it = std::find(sorted.begin(), sorted.end(), twinIdx);
    if (it == sorted.end()) continue;
    auto prevIt = (it == sorted.begin()) ? (sorted.end() - 1) : (it - 1);
    halfedges[i].next = *prevIt;
  }

  // Walk faces.
  int nFaces = 0;
  std::vector<int> faceStartHE;
  for (int i = 0; i < (int)halfedges.size(); ++i) {
    if (halfedges[i].face != -1) continue;
    int h = i;
    int safety = 0;
    do {
      if (halfedges[h].next == -1 || safety++ > (int)halfedges.size()) break;
      halfedges[h].face = nFaces;
      h = halfedges[h].next;
    } while (h != i);
    faceStartHE.push_back(i);
    ++nFaces;
  }

  // Identify the outer face: the one with most-negative signed area.
  // Signed area via shoelace on the face's vertex cycle.
  int outerFace = -1;
  i128 worstArea = 0;
  for (int f = 0; f < nFaces; ++f) {
    int h = faceStartHE[f];
    if (h < 0) continue;
    const Vec2i ref = verts[halfedges[h].origin];
    i128 area2 = 0;
    int hh = h;
    int safety = 0;
    do {
      const Vec2i p0 = verts[halfedges[hh].origin];
      const Vec2i p1 = verts[halfedges[halfedges[hh].twin].origin];
      area2 += (i128)(p0.x - ref.x) * (p1.y - ref.y) -
               (i128)(p1.x - ref.x) * (p0.y - ref.y);
      hh = halfedges[hh].next;
      if (hh < 0 || safety++ > (int)halfedges.size()) break;
    } while (hh != h);
    if (area2 < worstArea) { worstArea = area2; outerFace = f; }
  }

  // Per-face winding via ray-cast from a point inside.
  std::vector<int> faceWind(nFaces, 0);
  for (int f = 0; f < nFaces; ++f) {
    if (f == outerFace) continue;
    const int h = faceStartHE[f];
    if (h < 0) continue;
    const Vec2i a = verts[halfedges[h].origin];
    const Vec2i b = verts[halfedges[halfedges[h].twin].origin];
    // Point inside the face: midpoint of edge h, perp-LEFT-offset by 1.
    // Perp-LEFT of (b - a) is (-(b.y-a.y), (b.x-a.x)). For int coords
    // we step by sign(perp) units (1 in the larger axis); face must be
    // at least 1-unit wide in that direction or the offset misses it.
    // For typical inputs at scale 1e9, all faces have far more than
    // 1-unit width.
    const i64 dx = b.x - a.x;
    const i64 dy = b.y - a.y;
    const i64 mx = (a.x + b.x) / 2;
    const i64 my = (a.y + b.y) / 2;
    // Perp-LEFT: (-dy, dx). Step by sign of each component.
    i64 px, py;
    if (std::abs(dy) > std::abs(dx)) {
      // Edge is more vertical — step in x-direction (negative dy → +x or -x).
      px = mx + (dy > 0 ? -1 : 1);
      py = my;
    } else {
      px = mx;
      py = my + (dx > 0 ? 1 : -1);
    }
    faceWind[f] = CastWindingRay({px, py}, canon, verts);
  }

  // Filter: keep edges where left/right faces disagree on the predicate.
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
    if (leftIn) out.push_back({k.first, k.second, 1});
    else out.push_back({k.second, k.first, 1});
  }
  return out;
}

// ============================================================================
// Driver.
// ============================================================================

inline std::pair<std::vector<Vec2i>, std::vector<OutEdge>> RemoveOverlaps2D(
    const std::vector<Vec2i>& vertsIn, const std::vector<EdgeM>& edgesIn,
    const std::function<bool(int)>& isInside) {
  auto m = MergeVerts(vertsIn);
  auto e = RemapAndCollapse(edgesIn, m.remap);
  auto lists = BuildEdgeVertLists(e, m.verts);
  FindAndInsertIntersections(e, &m.verts, &lists);
  auto canon = Canonicalize(e, lists);
  auto out = FilterByWindingDCEL(canon, m.verts, isInside);
  return {std::move(m.verts), std::move(out)};
}

// ============================================================================
// Output edges → polygons (DCEL face cycle walk).
// ============================================================================

inline Polygons OutEdgesToPolygons(const std::vector<Vec2i>& verts,
                                    const std::vector<OutEdge>& edges,
                                    i64 scale) {
  const int nE = (int)edges.size();
  std::vector<std::vector<int>> outgoing(verts.size());
  for (int i = 0; i < nE; ++i) outgoing[edges[i].v0].push_back(i);
  for (size_t v = 0; v < outgoing.size(); ++v) {
    auto& lst = outgoing[v];
    if (lst.size() < 2) continue;
    const Vec2i vp = verts[v];
    std::sort(lst.begin(), lst.end(), [&](int a, int b) {
      const Vec2i dA = {verts[edges[a].v1].x - vp.x,
                        verts[edges[a].v1].y - vp.y};
      const Vec2i dB = {verts[edges[b].v1].x - vp.x,
                        verts[edges[b].v1].y - vp.y};
      auto bucket = [](const Vec2i& d) {
        if (d.y > 0 || (d.y == 0 && d.x > 0)) return 0;
        return 1;
      };
      const int bA = bucket(dA), bB = bucket(dB);
      if (bA != bB) return bA < bB;
      const i128 cross = (i128)dA.x * dB.y - (i128)dA.y * dB.x;
      return cross > 0;
    });
  }
  std::vector<bool> visited(nE, false);
  Polygons polys;
  const double inv = 1.0 / (double)scale;
  for (int start = 0; start < nE; ++start) {
    if (visited[start]) continue;
    SimplePolygon loop;
    int cur = start;
    while (cur >= 0 && !visited[cur]) {
      visited[cur] = true;
      const Vec2i p = verts[edges[cur].v0];
      loop.push_back({(double)p.x * inv, (double)p.y * inv});
      const int destV = edges[cur].v1;
      if (destV < 0 || destV >= (int)outgoing.size() ||
          outgoing[destV].empty()) {
        cur = -1;
        break;
      }
      const Vec2i vp = verts[destV];
      const Vec2i inDir = {vp.x - verts[edges[cur].v0].x,
                           vp.y - verts[edges[cur].v0].y};
      const auto& lst = outgoing[destV];
      // Find unvisited entry with smallest CCW delta from reverse-incoming.
      int next = -1;
      for (int e : lst) {
        if (visited[e]) continue;
        // Pick the one that's "smallest left turn" from inDir. Use orient
        // sign relative to inDir; we need the CCW-nearest.
        // Simpler: walk lst, for each candidate compute angle delta;
        // pick min positive. With int coords, compare via cross product
        // signs only — exact.
        // For the spike, just take any unvisited edge — works for simple
        // arrangements. (Production needs proper face walking.)
        next = e;
        break;
      }
      cur = next;
    }
    if (loop.size() >= 3) polys.push_back(std::move(loop));
  }
  return polys;
}

// ============================================================================
// Public API.
// ============================================================================

inline Polygons Simplify(const Polygons& in, i64 scale = kDefaultScale) {
  auto [verts, edges] = PolygonsToInput(in, scale);
  if (verts.empty()) return {};
  auto pred = [](int w) { return w != 0; };  // NONZERO fill
  auto [outVerts, outEdges] =
      RemoveOverlaps2D(verts, edges, pred);
  return OutEdgesToPolygons(outVerts, outEdges, scale);
}

inline double SignedArea(const SimplePolygon& p) {
  double a = 0;
  const int n = (int)p.size();
  for (int i = 0; i < n; ++i) {
    const auto& p0 = p[i];
    const auto& p1 = p[(i + 1) % n];
    a += p0.x * p1.y - p1.x * p0.y;
  }
  return a * 0.5;
}

inline double TotalSignedArea(const Polygons& polys) {
  double a = 0;
  for (const auto& loop : polys) a += SignedArea(loop);
  return a;
}

}  // namespace overlap2d_int

// ============================================================================
// Bigger workload: synthesize a self-intersecting polygon with N vertices
// arranged on a circle but visited in a star-skip order, producing many
// crossings. Useful for both correctness checking (does the int kernel
// stay correct as N grows?) and timing.
// ============================================================================

namespace overlap2d_int {

inline manifold::Polygons MakeStar(int n, int skip, double radius = 100.0) {
  manifold::SimplePolygon loop;
  for (int i = 0; i < n; ++i) {
    int j = (i * skip) % n;
    if (n > 1) {
      // ensure we visit n distinct verts when gcd(skip, n) == 1
    }
    const double ang = 2 * M_PI * j / n;
    loop.push_back({radius * std::cos(ang), radius * std::sin(ang)});
  }
  return {loop};
}

}  // namespace overlap2d_int

// ============================================================================
// Main: smoke tests + a head-to-head timing vs Clipper64 (when built with
// -DOVERLAP2D_INT_WITH_CLIPPER2).
// ============================================================================
#include <chrono>

#ifdef OVERLAP2D_INT_WITH_CLIPPER2
#include "clipper2/clipper.h"
namespace overlap2d_int {

inline Clipper2Lib::Paths64 ToPaths64(const manifold::Polygons& p,
                                      i64 scale) {
  Clipper2Lib::Paths64 out;
  out.reserve(p.size());
  for (const auto& loop : p) {
    Clipper2Lib::Path64 path;
    path.reserve(loop.size());
    for (const auto& v : loop) {
      path.emplace_back((i64)std::llround(v.x * scale),
                        (i64)std::llround(v.y * scale));
    }
    out.push_back(std::move(path));
  }
  return out;
}

inline manifold::Polygons FromPaths64(const Clipper2Lib::Paths64& p,
                                       i64 scale) {
  manifold::Polygons out;
  out.reserve(p.size());
  const double inv = 1.0 / (double)scale;
  for (const auto& path : p) {
    manifold::SimplePolygon loop;
    loop.reserve(path.size());
    for (const auto& pt : path) loop.push_back({pt.x * inv, pt.y * inv});
    out.push_back(std::move(loop));
  }
  return out;
}

inline manifold::Polygons ClipperSimplify(const manifold::Polygons& in,
                                          i64 scale) {
  Clipper2Lib::Clipper64 c;
  c.AddSubject(ToPaths64(in, scale));
  Clipper2Lib::Paths64 sol;
  c.Execute(Clipper2Lib::ClipType::Union,
            Clipper2Lib::FillRule::NonZero, sol);
  return FromPaths64(sol, scale);
}

}  // namespace overlap2d_int
#endif

namespace {

void RunCase(const char* name, const manifold::Polygons& in) {
  auto t0 = std::chrono::steady_clock::now();
  auto out = overlap2d_int::Simplify(in);
  auto t1 = std::chrono::steady_clock::now();
  const double inA = std::fabs(overlap2d_int::TotalSignedArea(in));
  const double outA = std::fabs(overlap2d_int::TotalSignedArea(out));
  const double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
  std::printf("=== %s ===\n", name);
  std::printf("  in:  %zu loop(s), |area|=%.4f\n", in.size(), inA);
  std::printf("  out: %zu loop(s), |area|=%.4f  (%.1f µs)\n", out.size(), outA,
              us);
}

}  // namespace

int main() {
  using manifold::Polygons;
  using manifold::SimplePolygon;

  // 1. Bowtie
  RunCase("bowtie", Polygons{
                        SimplePolygon{{0, 0}, {10, 10}, {10, 0}, {0, 10}}});

  // 2. Two overlapping squares (UNION+NONZERO equivalent on combined input)
  RunCase("two_overlapping_squares",
          Polygons{
              SimplePolygon{{0, 0}, {4, 0}, {4, 4}, {0, 4}},
              SimplePolygon{{2, 2}, {6, 2}, {6, 6}, {2, 6}}});

  // 3. Three nested triangles (mfogel issue-44)
  RunCase("issue_44_three_nested_triangles",
          Polygons{
              SimplePolygon{{0, 0}, {3, -2}, {5, 0}},
              SimplePolygon{{0, 0}, {3, -1}, {4, 0}},
              SimplePolygon{{0, 0}, {3, -3}, {5, 0}}});

  // 4. Pinched rectangle
  RunCase("pinched_rectangle",
          Polygons{
              SimplePolygon{{0, 0}, {4, 0}, {4, 4}, {2, 2}, {0, 4}}});

  // 5. Clean non-overlapping square (sanity)
  RunCase("clean_square",
          Polygons{
              SimplePolygon{{0, 0}, {10, 0}, {10, 10}, {0, 10}}});

  // 6. Heavy self-crossing star polygon — exercises the kernel at scale.
  //    {n/skip} = visit n verts on a circle in (skip)-step order.
  for (int n : {12, 24, 50}) {
    int skip = (n / 2) - 1;  // skip-(n/2-1) gives many crossings
    if (skip < 2) skip = 3;
    auto star = overlap2d_int::MakeStar(n, skip);
    char name[64];
    std::snprintf(name, sizeof(name), "star_n%d_skip%d", n, skip);
    RunCase(name, star);
  }

#ifdef OVERLAP2D_INT_WITH_CLIPPER2
  // Head-to-head vs Clipper64 on a larger workload. Repeated K times so the
  // wall-clock signal isn't dominated by setup.
  std::printf("\n=== head-to-head vs Clipper64 ===\n");
  struct Bench {
    const char* name;
    manifold::Polygons input;
    int reps;
  };
  std::vector<Bench> benches = {
      {"two_overlapping_squares",
       Polygons{SimplePolygon{{0, 0}, {4, 0}, {4, 4}, {0, 4}},
                SimplePolygon{{2, 2}, {6, 2}, {6, 6}, {2, 6}}},
       10000},
      {"issue_44_three_nested_triangles",
       Polygons{SimplePolygon{{0, 0}, {3, -2}, {5, 0}},
                SimplePolygon{{0, 0}, {3, -1}, {4, 0}},
                SimplePolygon{{0, 0}, {3, -3}, {5, 0}}},
       10000},
      {"star_n12", overlap2d_int::MakeStar(12, 5), 5000},
      {"star_n24", overlap2d_int::MakeStar(24, 11), 1000},
      {"star_n50", overlap2d_int::MakeStar(50, 23), 200},
  };
  for (const auto& b : benches) {
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();
    for (int i = 0; i < b.reps; ++i) {
      auto _ = overlap2d_int::Simplify(b.input);
      (void)_;
    }
    auto t1 = Clock::now();
    for (int i = 0; i < b.reps; ++i) {
      auto _ = overlap2d_int::ClipperSimplify(b.input,
                                              overlap2d_int::kDefaultScale);
      (void)_;
    }
    auto t2 = Clock::now();
    const double oursMs =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double c2Ms =
        std::chrono::duration<double, std::milli>(t2 - t1).count();
    const double oursUs = (oursMs * 1000.0) / b.reps;
    const double c2Us = (c2Ms * 1000.0) / b.reps;
    std::printf("  %-40s ours=%.2f us  c2=%.2f us  ratio=%.2fx %s\n", b.name,
                oursUs, c2Us, oursUs / c2Us,
                oursUs < c2Us ? "(ours faster)" : "(c2 faster)");
  }
#endif

  return 0;
}
