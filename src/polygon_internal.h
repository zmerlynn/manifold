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

#pragma once

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "manifold/polygon.h"
#include "shared.h"

namespace manifold {

struct HalfedgeTriangulation {
  std::vector<Halfedge> halfedges;
  size_t contourEnd = 0;
  double epsilon = -1;

  void AddContours(const PolygonsIdx& polys) {
    size_t numContourEdges = 0;
    for (const SimplePolygonIdx& poly : polys) numContourEdges += poly.size();
    halfedges.reserve(halfedges.size() + numContourEdges);
    edge2halfedge.reserve(edge2halfedge.size() + numContourEdges);
    for (const SimplePolygonIdx& poly : polys) {
      for (size_t i = 0; i < poly.size(); ++i) {
        const int start = poly[i].idx;
        const int end = poly[i + 1 < poly.size() ? i + 1 : 0].idx;
        // Store the exterior contour halfedge, opposite the filled contour.
        AddHalfedge(end, start);
      }
    }
    contourEnd = halfedges.size();
  }

  void ReserveTriangles(size_t numTri) {
    halfedges.reserve(contourEnd + 3 * numTri);
    edge2halfedge.reserve(edge2halfedge.size() + numTri);
  }

  void AddTriangle(int first, int second, int third) {
    AddHalfedge(first, second);
    AddHalfedge(second, third);
    AddHalfedge(third, first);
  }

  size_t NumTri() const { return (halfedges.size() - contourEnd) / 3; }

  std::vector<ivec3> Triangles() const {
    std::vector<ivec3> triangles;
    triangles.reserve(NumTri());
    for (size_t edge = contourEnd; edge < halfedges.size(); edge += 3) {
      triangles.push_back({halfedges[edge].startVert,
                           halfedges[edge + 1].startVert,
                           halfedges[edge + 2].startVert});
    }
    return triangles;
  }

  void Finalize() {
#ifdef MANIFOLD_DEBUG
    DEBUG_ASSERT(edge2halfedge.empty(), topologyErr,
                 "triangulation has unpaired halfedges");
    for (size_t i = 0; i < halfedges.size(); ++i) {
      const int pair = halfedges[i].pairedHalfedge;
      DEBUG_ASSERT(pair >= 0 && pair < static_cast<int>(halfedges.size()),
                   topologyErr, "invalid paired halfedge");
      DEBUG_ASSERT(halfedges[pair].pairedHalfedge == static_cast<int>(i),
                   topologyErr, "halfedge pair is not reciprocal");
      DEBUG_ASSERT(halfedges[i].startVert == halfedges[pair].endVert &&
                       halfedges[i].endVert == halfedges[pair].startVert,
                   topologyErr, "halfedge pair endpoints do not match");
    }
#endif
    edge2halfedge.clear();
    edge2halfedge.rehash(0);
  }

 private:
  std::unordered_map<uint64_t, std::vector<int>> edge2halfedge;

  static uint64_t EdgeKey(int start, int end) {
    return (uint64_t{static_cast<uint32_t>(start)} << 32) |
           static_cast<uint32_t>(end);
  }

  void AddHalfedge(int start, int end) {
    const int halfedge = halfedges.size();
    Halfedge data = {start, end, -1, -1};
    auto reverse = edge2halfedge.find(EdgeKey(end, start));
    if (reverse != edge2halfedge.end() && !reverse->second.empty()) {
      data.pairedHalfedge = reverse->second.back();
      halfedges[data.pairedHalfedge].pairedHalfedge = halfedge;
      reverse->second.pop_back();
      if (reverse->second.empty()) edge2halfedge.erase(reverse);
    } else {
      edge2halfedge[EdgeKey(start, end)].push_back(halfedge);
    }
    halfedges.push_back(data);
  }
};

HalfedgeTriangulation TriangulateIdxHalfedges(const PolygonsIdx& polys,
                                              double epsilon = -1,
                                              bool allowConvex = true);

// ---------------------------------------------------------------------------
// EXACT (opt-in) triangulation mode - manifold::exacttri (polygon.cpp).
//
// The exact counterpart of the tolerance triangulator above.  Operates on 3D
// positions projected by DOMINANT-AXIS DROP (pure coordinate selection -
// exact), with every orientation decision made through the one blessed exact
// drop-frame orient2d form below.  Opt-in via a separate entry point: the
// default tolerance path (TriangulateIdx / EarClip) is byte-for-byte
// unaffected.  See the module comment in polygon.cpp for WHY the exact mode
// exists (it is a boundary-pairing contract, not a triangle-quality knob).
// ---------------------------------------------------------------------------

// The one blessed exact drop-frame orient2d FORM (filter-first, exact
// escalation) - defined in overlap3.cpp beside the exact kernel (one
// predicate, one implementation).
int ExactOrient2DDrop(const vec3& p, const vec3& q, const vec3& r, int axis);

namespace exacttri {

// 2D projection by dominant-axis DROP - pure coordinate selection (exact),
// matching ExactOrient2DDrop's convention.
inline vec2 Drop2(const vec3& v, int axis) {
  if (axis == 0) return {v.y, v.z};
  if (axis == 1) return {v.x, v.z};
  return {v.x, v.y};
}

// Exact sign of the projected orientation (a,b,c) - the landed exact-on-
// doubles orient2d (filter-first).  NEGATED: ExactOrient2DDrop's raw
// determinant (orient3d against a +axis lift) is NEGATIVE for a CCW triple
// in the dropped frame (its production callers are straddle-only,
// sign-agnostic); this module needs the standard CCW-positive convention
// (verified: the un-negated form traced every group's OUTER contour as the
// positive loop and failed every ear test).
inline int O2(const vec3& a, const vec3& b, const vec3& c, int axis) {
  return -ExactOrient2DDrop(a, b, c, axis);
}

// v strictly interior to the open segment (a,b) in the projected frame:
// exactly collinear and strictly between in the wider coordinate.
inline bool OnOpenSeg2(const vec3& a, const vec3& b, const vec3& v, int axis) {
  if (O2(a, b, v, axis) != 0) return false;
  const vec2 a2 = Drop2(a, axis), b2 = Drop2(b, axis), v2 = Drop2(v, axis);
  if (std::abs(b2.x - a2.x) >= std::abs(b2.y - a2.y))
    return (a2.x < v2.x) != (b2.x < v2.x);
  return (a2.y < v2.y) != (b2.y < v2.y);
}

// Proper crossing of open segments (p,q) x (a,b) in the projected frame.
inline bool ProperCross2(const vec3& p, const vec3& q, const vec3& a,
                         const vec3& b, int axis) {
  const int d1 = O2(p, q, a, axis), d2 = O2(p, q, b, axis);
  const int d3 = O2(a, b, p, axis), d4 = O2(a, b, q, axis);
  return d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0 && d1 != d2 && d3 != d4;
}

// Signed shoelace sum (2x area) of a projected loop, accumulated about the
// loop's OWN first vertex.  The shoelace is translation-invariant in exact
// arithmetic; translating collapses the roundoff floor from ulp(|coord|^2)
// (raw terms ~coord^2 cancel catastrophically) to ~ulp(span^2).  On raw
// coordinates a micro cell far from the origin has |true s| at or below the
// term noise and its SIGN is garbage - a real CCW cell then misroutes into
// the hole-ring/dust arms and its half-edges vanish unpaired (measured: the
// openscad double-vertex-fan tip triangle, |s|=5.4e-14 against term ulp
// 5.7e-14, silently dropped -> the corner stub family).
inline double LoopShoelace(const std::vector<int>& L,
                           const std::vector<vec3>& pos3, int axis) {
  if (L.empty()) return 0.0;
  const vec2 o = Drop2(pos3[L[0]], axis);
  double s = 0.0;
  for (size_t k = 0; k < L.size(); ++k) {
    const vec2 p1 = Drop2(pos3[L[k]], axis) - o;
    const vec2 p2 = Drop2(pos3[L[(k + 1) % L.size()]], axis) - o;
    s += p1.x * p2.y - p2.x * p1.y;
  }
  return s;
}

// EARCLIP-first exact triangulation of a weakly-simple CCW polygon (collinear
// runs, pinch-repeated vertices, keyhole-duplicated bridges all allowed).
// loop = vertex indices into pos3; emits index triples into out.  eps is the
// weld radius (dust adjudications).  Returns false when the loop cannot be
// covered (the caller adjudicates dust vs macro failure).
bool Triangulate(const std::vector<int>& loop, const std::vector<vec3>& pos3,
                 int axis, double eps, std::vector<ivec3>& out);

}  // namespace exacttri

}  // namespace manifold
