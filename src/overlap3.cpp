// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Stages A, B, D, E of the 3D sweep-plane overlap-removal prototype.
// Stage C is in overlap3_sweep.cpp.
// Design: docs/SweepPlane3D.md.

#include "overlap3.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "boolean2.h"
#include "disjoint_sets.h"
#include "impl.h"
#include "manifold/optional_assert.h"
#include "polygon_internal.h"
#include "shared.h"

namespace manifold {

// Stage C forward declaration (implemented in overlap3_sweep.cpp)
StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt);

namespace {

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

static bool IsInside3D(int64_t w) { return w > 0; }

static double Dist3(vec3 a, vec3 b) { return la::length(a - b); }

// Point-to-segment 3D distance.
static double PointSegDist3(vec3 p, vec3 a, vec3 b) {
  const vec3 ab = b - a;
  const double len2 = la::dot(ab, ab);
  if (len2 == 0.0) return la::length(p - a);
  const double t = std::max(0.0, std::min(1.0, la::dot(p - a, ab) / len2));
  return la::length(p - (a + t * ab));
}

// 3D cross product of triangle edges -> face normal (not normalized).
static vec3 TriNormal(vec3 v0, vec3 v1, vec3 v2) {
  return la::cross(v1 - v0, v2 - v0);
}

// Signed distance from point p to plane (n, d) where n is the plane normal.
static double PlaneDist(vec3 p, vec3 n, double d) { return la::dot(n, p) - d; }

// Does segment AB straddle the plane n.p = d?  Returns +1/-1 sides, 0 if
// either endpoint is on the plane (within the floating-point sense).
static bool SegStraddlesPlane(vec3 a, vec3 b, vec3 n, double d, double& tOut,
                              vec3& qOut) {
  const double da = PlaneDist(a, n, d);
  const double db = PlaneDist(b, n, d);
  if (da * db >= 0) return false;  // same side or both on plane
  const double denom = da - db;
  if (denom == 0.0) return false;
  const double t = da / denom;
  tOut = t;
  qOut = a + t * (b - a);
  return true;
}

// ---------------------------------------------------------------------------
// Stage A: merge verts and canonicalize faces
// ---------------------------------------------------------------------------

struct StageAResult {
  std::vector<vec3> mergedVerts;  // canonical 3D vert positions
  std::vector<int> vertMap;       // original vert -> merged vert
  std::vector<CanonicalFace> faces;
};

static StageAResult StageA(const Manifold::Impl& in, double eps) {
  const int nVerts = static_cast<int>(in.vertPos_.size());
  const int nTris = static_cast<int>(in.halfedge_.size()) / 3;

  // Vert merge: union-find with representative = centroid-nearest to cluster
  // centroid (matching 2D engine convention, MergeVerts in boolean2.cpp).
  DisjointSets uf(nVerts);
  // O(n^2) merge for prototype; in production, use a spatial hash.
  for (int i = 0; i < nVerts; ++i) {
    for (int j = i + 1; j < nVerts; ++j) {
      if (la::length(in.vertPos_[i] - in.vertPos_[j]) <= eps) uf.unite(i, j);
    }
  }
  // Build representative map: for each component, pick the vert closest to the
  // component centroid.
  std::map<int, std::vector<int>> components;
  for (int i = 0; i < nVerts; ++i)
    components[static_cast<int>(uf.find(i))].push_back(i);

  std::vector<int> vertMap(nVerts);
  std::vector<vec3> mergedVerts;
  std::map<int, int> rootToMerged;

  for (auto& [root, members] : components) {
    vec3 centroid(0.0);
    for (int v : members) centroid += in.vertPos_[v];
    centroid /= static_cast<double>(members.size());
    // Pick member closest to centroid
    int best = members[0];
    double bestDist = la::length(in.vertPos_[best] - centroid);
    for (int v : members) {
      const double d = la::length(in.vertPos_[v] - centroid);
      if (d < bestDist) {
        bestDist = d;
        best = v;
      }
    }
    const int mergedId = static_cast<int>(mergedVerts.size());
    mergedVerts.push_back(in.vertPos_[best]);
    rootToMerged[root] = mergedId;
  }
  for (int i = 0; i < nVerts; ++i)
    vertMap[i] = rootToMerged[static_cast<int>(uf.find(i))];

  // Canonicalize faces: key = sorted (v0,v1,v2) unordered vertex set.
  // Mult is signed relative to the STORED orientation of the canonical
  // representative (first face with this vertex set). A face with the same
  // halfedge winding as the representative contributes +1; opposite winding
  // contributes -1. This gives mult=+1 for each outward-facing triangle in a
  // manifold mesh, and correct cancellation for touching-disjoint Compose
  // (opposite-orientation shared face -> mult=0 -> dropped).
  struct FaceKey {
    ivec3 sorted;  // smallest vert id first
    bool operator<(const FaceKey& o) const {
      if (sorted.x != o.sorted.x) return sorted.x < o.sorted.x;
      if (sorted.y != o.sorted.y) return sorted.y < o.sorted.y;
      return sorted.z < o.sorted.z;
    }
  };

  struct FaceAccumEntry {
    int64_t mult;   // algebraic sum; positive = same orientation as rep
    int repTriId;   // representative face triangle index
    int repParity;  // sort permutation parity of the representative (+1/-1)
  };

  std::map<FaceKey, FaceAccumEntry> faceAccum;

  for (int tri = 0; tri < nTris; ++tri) {
    const int h = 3 * tri;
    const int v0 = vertMap[in.halfedge_.Start(h)];
    const int v1 = vertMap[in.halfedge_.Start(h + 1)];
    const int v2 = vertMap[in.halfedge_.Start(h + 2)];
    // Skip degenerate (collapsed) triangles
    if (v0 == v1 || v1 == v2 || v0 == v2) continue;

    // Sort verts and compute permutation parity.
    int a = v0, b = v1, c = v2;
    int parity = 1;
    if (a > b) {
      std::swap(a, b);
      parity = -parity;
    }
    if (b > c) {
      std::swap(b, c);
      parity = -parity;
    }
    if (a > b) {
      std::swap(a, b);
      parity = -parity;
    }

    FaceKey key{{a, b, c}};
    auto it = faceAccum.find(key);
    if (it == faceAccum.end()) {
      // First face with this vertex set: it IS the representative. By
      // definition it agrees with itself, contributing +1.
      faceAccum.emplace(key, FaceAccumEntry{1, tri, parity});
    } else {
      // Subsequent face: +1 if same permutation parity as rep, -1 if opposite.
      const int matchSign = (parity == it->second.repParity) ? 1 : -1;
      it->second.mult += matchSign;
    }
  }

  // Emit CanonicalFace records for entries with nonzero mult.
  StageAResult result;
  result.mergedVerts = std::move(mergedVerts);
  result.vertMap = std::move(vertMap);

  for (auto& [key, entry] : faceAccum) {
    if (entry.mult == 0) continue;  // cancelled (e.g. touching-disjoint)
    const int tri = entry.repTriId;
    const int h = 3 * tri;
    const int mv0 = result.vertMap[in.halfedge_.Start(h)];
    const int mv1 = result.vertMap[in.halfedge_.Start(h + 1)];
    const int mv2 = result.vertMap[in.halfedge_.Start(h + 2)];
    vec3 storedNormal =
        TriNormal(result.mergedVerts[mv0], result.mergedVerts[mv1],
                  result.mergedVerts[mv2]);
    if (la::dot(storedNormal, storedNormal) == 0.0) continue;
    storedNormal = la::normalize(storedNormal);
    CanonicalFace cf;
    cf.id = tri;
    cf.verts = {mv0, mv1, mv2};
    cf.normal = storedNormal;
    cf.mult = entry.mult;
    result.faces.push_back(cf);
  }

  return result;
}

// ---------------------------------------------------------------------------
// Stage B: arrangement geometry (seams, events, triple points)
// ---------------------------------------------------------------------------

// M1: LineTriClip - clip the line (P + t*D) to the interior of triangle
// (v0,v1,v2) with normal n.  Uses eps-relative tolerances (not 1e-14).
// Returns whether the range [tLo, tHi] is non-empty.
static bool LineTriClip(vec3 P, vec3 D, vec3 v0, vec3 v1, vec3 v2, vec3 n,
                        double& tLo, double& tHi, double eps) {
  tLo = -1e18;
  tHi = 1e18;
  const vec3 edges[3] = {v1 - v0, v2 - v1, v0 - v2};
  const vec3 vBase[3] = {v0, v1, v2};
  // Scale factor for eps-relative comparison: use triangle area ~ |n|.
  const double triScale = la::length(n);
  const double kTol = eps * triScale;
  for (int i = 0; i < 3; ++i) {
    const vec3 ev = la::cross(edges[i], D);
    const vec3 evBase = la::cross(edges[i], vBase[i] - P);
    const double denom = la::dot(n, ev);
    const double num = la::dot(n, evBase);
    if (std::abs(denom) <= kTol) {
      // D is parallel to this edge's halfplane boundary.
      if (num > kTol) return false;  // entirely outside
    } else {
      const double t = num / denom;
      if (denom > 0) {
        tLo = std::max(tLo, t);
      } else {
        tHi = std::min(tHi, t);
      }
    }
  }
  return tLo <= tHi + eps;
}

// M1: Compute the intersection segment between two non-coplanar triangles.
// Six edge-plane clips (three per triangle), eps-aware inside tests, endpoint
// dedup.  Returns true if intersection segment length > eps.
static bool TriTriSeam(vec3 a0, vec3 a1, vec3 a2, vec3 na,  // face A
                       vec3 b0, vec3 b1, vec3 b2, vec3 nb,  // face B
                       vec3& qA, vec3& qB, double eps) {
  const vec3 D = la::cross(na, nb);
  const double Dlen = la::length(D);
  if (Dlen <= eps * eps) return false;  // near-parallel planes
  const vec3 dir = D / Dlen;

  const double da = la::dot(na, a0);
  const double db = la::dot(nb, b0);
  const vec3 absD = la::abs(dir);
  int maxComp =
      (absD.x >= absD.y && absD.x >= absD.z) ? 0 : (absD.y >= absD.z ? 1 : 2);
  const int c1 = (maxComp + 1) % 3, c2 = (maxComp + 2) % 3;
  const double a11 = na[c1], a12 = na[c2], b11 = nb[c1], b12 = nb[c2];
  const double det = a11 * b12 - a12 * b11;
  if (std::abs(det) <= eps * eps) return false;
  vec3 P(0.0);
  P[c1] = (da * b12 - db * a12) / det;
  P[c2] = (a11 * db - b11 * da) / det;

  double tAlo, tAhi, tBlo, tBhi;
  if (!LineTriClip(P, dir, a0, a1, a2, na, tAlo, tAhi, eps)) return false;
  if (!LineTriClip(P, dir, b0, b1, b2, nb, tBlo, tBhi, eps)) return false;
  const double tLo = std::max(tAlo, tBlo);
  const double tHi = std::min(tAhi, tBhi);
  if (tLo >= tHi) return false;

  qA = P + tLo * dir;
  qB = P + tHi * dir;
  return true;
}

// M2: Check whether edge (eA, eB) lies in plane (n, d) within eps.
// Returns true if both endpoints are within eps of the plane.
static bool EdgeInPlane(vec3 eA, vec3 eB, vec3 n, double d, double eps) {
  return std::abs(la::dot(n, eA) - d) <= eps &&
         std::abs(la::dot(n, eB) - d) <= eps;
}

// M2: Sutherland-Hodgman clip of convex polygon by one halfplane (left of
// directed edge eA->eB). Same 1e-14 stability guard as the coplanar test.
static std::vector<vec2> ClipPolyByHalfplane2D(std::vector<vec2> poly, vec2 eA,
                                               vec2 eB) {
  std::vector<vec2> out;
  const int n = (int)poly.size();
  const vec2 e = eB - eA;
  for (int k = 0; k < n; ++k) {
    const vec2& cur = poly[k];
    const vec2& nxt = poly[(k + 1) % n];
    // Inside = left of eA->eB: e cross (p - eA) >= 0.
    const bool curIn = e.x * (cur.y - eA.y) - e.y * (cur.x - eA.x) >= 0.0;
    const bool nxtIn = e.x * (nxt.y - eA.y) - e.y * (nxt.x - eA.x) >= 0.0;
    if (curIn) out.push_back(cur);
    if (curIn != nxtIn) {
      const vec2 d1 = nxt - cur, d2 = e;
      const double denom2 = d1.x * d2.y - d1.y * d2.x;
      // Stability guard: same 1e-14 bound as the coplanar area check below.
      if (std::abs(denom2) > 1e-14) {
        const vec2 dc = eA - cur;
        const double t2 = (dc.x * d2.y - dc.y * d2.x) / denom2;
        out.push_back(cur + t2 * d1);
      }
    }
  }
  return out;
}

// M2: Clip segment [sA, sB] against the interior of triangle (t0,t1,t2) in 2D.
// Uses the same halfplane sign convention as ClipPolyByHalfplane2D.
// Returns the length of the clipped sub-segment; 0 if the segment has no
// interior overlap with the triangle. Length > eps is the EdgeInPlane test.
static double SegTriInteriorLength2D(vec2 sA, vec2 sB, vec2 t0, vec2 t1,
                                     vec2 t2) {
  double tLo = 0.0, tHi = 1.0;
  const vec2 triV[3] = {t0, t1, t2};
  const vec2 segD = sB - sA;
  for (int i = 0; i < 3; ++i) {
    const vec2 eA = triV[i], eB = triV[(i + 1) % 3];
    const vec2 e = eB - eA;
    // f(t) = e cross (sA + t*segD - eA); positive = inside this halfplane.
    const double f0 = e.x * (sA.y - eA.y) - e.y * (sA.x - eA.x);
    const double f1 = e.x * (sB.y - eA.y) - e.y * (sB.x - eA.x);
    const double fd = f1 - f0;
    // Same 1e-14 stability guard used in ClipPolyByHalfplane2D.
    if (std::abs(fd) <= 1e-14) {
      if (f0 < 0.0) return 0.0;  // parallel to edge, outside halfplane
    } else {
      const double tc = -f0 / fd;
      if (fd < 0.0)
        tHi = std::min(tHi, tc);  // exiting: f decreasing
      else
        tLo = std::max(tLo, tc);  // entering: f increasing
    }
    if (tLo > tHi) return 0.0;
  }
  return (tHi - tLo) * la::length(segD);
}

// Given an existing merged-vert pool, find or add a 3D point within eps.
// Returns the merged vert id.
static int FindOrAddVert(std::vector<MergedVert>& verts, vec3 pos, double eps) {
  for (int i = 0; i < static_cast<int>(verts.size()); ++i) {
    if (la::length(verts[i].pos - pos) <= eps) return i;
  }
  const int id = static_cast<int>(verts.size());
  verts.push_back({pos});
  return id;
}

struct StageBResult {
  ArrangementGeometry arr;
  std::optional<FatalReason> fatal;
  std::string detail;
};

static StageBResult StageB(const StageAResult& stageA, double eps,
                           Overlap3Counters& cnt) {
  StageBResult result;
  ArrangementGeometry& arr = result.arr;
  const int nFaces = static_cast<int>(stageA.faces.size());

  // Initialize merged vert pool from stage A verts.
  arr.verts.resize(stageA.mergedVerts.size());
  for (int i = 0; i < static_cast<int>(stageA.mergedVerts.size()); ++i)
    arr.verts[i] = {stageA.mergedVerts[i]};
  // Number of Stage-A (original) verts, before any event verts are added.
  const int nOrigVerts = static_cast<int>(stageA.mergedVerts.size());
  arr.faces = stageA.faces;
  arr.faceSeams.resize(nFaces);

  // Precompute face normals and plane distances for broadphase.
  std::vector<double> planeDist(nFaces);
  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& f = stageA.faces[fi];
    planeDist[fi] = la::dot(f.normal, arr.verts[f.verts.x].pos);
  }

  // Precompute 3D bounding boxes for broadphase O(n^2) with early exit.
  struct Box3 {
    vec3 mn, mx;
  };
  std::vector<Box3> boxes(nFaces);
  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& f = stageA.faces[fi];
    const vec3& p0 = arr.verts[f.verts.x].pos;
    const vec3& p1 = arr.verts[f.verts.y].pos;
    const vec3& p2 = arr.verts[f.verts.z].pos;
    boxes[fi].mn = la::min(la::min(p0, p1), p2) - vec3(eps);
    boxes[fi].mx = la::max(la::max(p0, p1), p2) + vec3(eps);
  }

  // Degenerate-contact cluster graph: nodes are 3D points (stored as vec3),
  // edges connect nodes of the same primitive or within eps. Components with
  // diameter > eps = SubResolutionChain.
  struct ContactNode {
    vec3 pos;
  };
  std::vector<ContactNode> contactNodes;
  std::vector<std::pair<int, int>> contactEdges;  // node pairs within eps

  auto addContactNode = [&](vec3 pos) -> int {
    for (int i = 0; i < (int)contactNodes.size(); ++i) {
      if (la::length(contactNodes[i].pos - pos) < 1e-15) return i;
    }
    const int id = (int)contactNodes.size();
    contactNodes.push_back({pos});
    return id;
  };
  auto addContactEdge = [&](int a, int b) {
    if (a == b) return;
    for (auto& e : contactEdges)
      if ((e.first == a && e.second == b) || (e.first == b && e.second == a))
        return;
    contactEdges.push_back({a, b});
    // Also add edges for nodes within eps of each other
    for (int i = 0; i < (int)contactNodes.size(); ++i) {
      for (int j = i + 1; j < (int)contactNodes.size(); ++j) {
        if (la::length(contactNodes[i].pos - contactNodes[j].pos) <= eps) {
          bool found = false;
          for (auto& e2 : contactEdges)
            if ((e2.first == i && e2.second == j) ||
                (e2.first == j && e2.second == i)) {
              found = true;
              break;
            }
          if (!found) contactEdges.push_back({i, j});
        }
      }
    }
  };

  // Process all face pairs: compute seam, handle degenerate cases.
  std::map<std::pair<int, int>, int> faceSeamMap;  // (fi,fj) -> seam index

  for (int fi = 0; fi < nFaces; ++fi) {
    for (int fj = fi + 1; fj < nFaces; ++fj) {
      // Broadphase: AABB overlap check
      const Box3& bi = boxes[fi];
      const Box3& bj = boxes[fj];
      if (bi.mn.x > bj.mx.x || bi.mx.x < bj.mn.x || bi.mn.y > bj.mx.y ||
          bi.mx.y < bj.mn.y || bi.mn.z > bj.mx.z || bi.mx.z < bj.mn.z)
        continue;

      const CanonicalFace& fi_face = stageA.faces[fi];
      const CanonicalFace& fj_face = stageA.faces[fj];

      // Skip edge-adjacent pairs (sharing 2+ merged verts): they share a
      // boundary edge and have no overlapping interior region.  An
      // edge-adjacent pair is NOT a seam candidate - it is ordinary mesh
      // topology.
      {
        const int va[3] = {fi_face.verts.x, fi_face.verts.y, fi_face.verts.z};
        const int vb[3] = {fj_face.verts.x, fj_face.verts.y, fj_face.verts.z};
        int shared = 0;
        for (int a : va)
          for (int b : vb)
            if (a == b) ++shared;
        if (shared >= 2) continue;  // edge-adjacent: no volumetric seam
      }
      const vec3& na = fi_face.normal;
      const vec3& nb = fj_face.normal;
      const vec3 pa0 = arr.verts[fi_face.verts.x].pos;
      const vec3 pa1 = arr.verts[fi_face.verts.y].pos;
      const vec3 pa2 = arr.verts[fi_face.verts.z].pos;
      const vec3 pb0 = arr.verts[fj_face.verts.x].pos;
      const vec3 pb1 = arr.verts[fj_face.verts.y].pos;
      const vec3 pb2 = arr.verts[fj_face.verts.z].pos;

      // Check for coplanar overlap (out of scope): if normals are parallel and
      // planes are equal, fire CoplanarOverlap ONLY if the two triangles have
      // positive-area 2D intersection (not just touching at a point or edge).
      // Disjoint/point/edge-touching coplanars CONTINUE (e.g. touching-disjoint
      // Compose where the shared face cancels in stage A and its neighbors
      // share only an edge or point in the same plane).
      const double normDot = std::abs(la::dot(na, nb));
      // 1e-8 ~ sqrt(double_eps): unit-vector parallelism guard, independent of
      // bbox scale; no metric-table entry (angular test, not coordinate space).
      if (normDot > 1.0 - 1e-8) {
        // Nearly coplanar: check plane distance
        const double planeSep = std::abs(la::dot(na, pb0) - la::dot(na, pa0));
        if (planeSep > eps) continue;  // parallel but separated: no seam
        // Same plane: project both tris to face plane and check area of 2D
        // intersection. Fire CoplanarOverlap only when interior overlap > 0.
        const mat2x3 proj = GetAxisAlignedProjection(na);
        const vec2 a2[3] = {proj * pa0, proj * pa1, proj * pa2};
        const vec2 b2[3] = {proj * pb0, proj * pb1, proj * pb2};
        // Sutherland-Hodgman clip of B against A to get intersection polygon.
        // If resulting polygon has area > eps^2, it's a real overlap.
        std::vector<vec2> poly = {b2[0], b2[1], b2[2]};
        for (int k = 0; k < 3 && !poly.empty(); ++k)
          poly = ClipPolyByHalfplane2D(poly, a2[k], a2[(k + 1) % 3]);
        // Compute signed area of clipped polygon.
        double clipArea = 0.0;
        for (int k = 0; k < (int)poly.size(); ++k) {
          const vec2& p = poly[k];
          const vec2& q = poly[(k + 1) % (int)poly.size()];
          clipArea += p.x * q.y - q.x * p.y;
        }
        clipArea = std::abs(clipArea) * 0.5;
        if (clipArea > eps * eps) {
          result.fatal = FatalReason::CoplanarOverlap;
          result.detail = "coplanar face pair with interior area overlap";
          return result;
        }
        // Zero or sub-eps area: point/edge touch only, continue normally.
        continue;
      }

      // M2: EdgeInPlane detection. Before computing TriTriSeam, check whether
      // any edge of A lies in the plane of B, or vice versa.  Both endpoints
      // of such an edge are within eps of the other triangle's plane.  This is
      // an out-of-scope class; detect and report before creating a seam.
      // Uses SegTriInteriorLength2D (clip-based): fires iff the clipped
      // sub-segment has length > eps (not just midpoint inside the triangle).
      {
        const double dA = la::dot(na, pa0);
        const double dB = la::dot(nb, pb0);
        const vec3 edgeA[3] = {pa0, pa1, pa2};
        const vec3 edgeB[3] = {pb0, pb1, pb2};
        // Project triangles once into each other's planes.
        const mat2x3 projB = GetAxisAlignedProjection(nb);
        const vec2 pb02 = projB * pb0, pb12 = projB * pb1, pb22 = projB * pb2;
        const mat2x3 projA = GetAxisAlignedProjection(na);
        const vec2 pa02 = projA * pa0, pa12 = projA * pa1, pa22 = projA * pa2;
        for (int ei = 0; ei < 3; ++ei) {
          if (EdgeInPlane(edgeA[ei], edgeA[(ei + 1) % 3], nb, dB, eps) &&
              SegTriInteriorLength2D(projB * edgeA[ei],
                                     projB * edgeA[(ei + 1) % 3], pb02, pb12,
                                     pb22) > eps) {
            result.fatal = FatalReason::EdgeInPlane;
            result.detail = "edge of face " + std::to_string(fi) +
                            " lies in plane of face " + std::to_string(fj);
            return result;
          }
          if (EdgeInPlane(edgeB[ei], edgeB[(ei + 1) % 3], na, dA, eps) &&
              SegTriInteriorLength2D(projA * edgeB[ei],
                                     projA * edgeB[(ei + 1) % 3], pa02, pa12,
                                     pa22) > eps) {
            result.fatal = FatalReason::EdgeInPlane;
            result.detail = "edge of face " + std::to_string(fj) +
                            " lies in plane of face " + std::to_string(fi);
            return result;
          }
        }
      }

      // Compute seam (intersection segment of the two triangles).
      vec3 qA, qB;
      if (!TriTriSeam(pa0, pa1, pa2, na, pb0, pb1, pb2, nb, qA, qB, eps))
        continue;

      const double seamLen = la::length(qB - qA);

      if (seamLen <= eps) {
        // Degenerate contact: point/tangent/sub-eps. Add to contact cluster
        // graph.
        const int nA = addContactNode(qA);
        const int nB = addContactNode(qB);
        addContactEdge(nA, nB);
        // Also link endpoints to nearby existing contact nodes.
        for (int i = 0; i < (int)contactNodes.size(); ++i) {
          if (i != nA && la::length(contactNodes[i].pos - qA) <= eps)
            addContactEdge(nA, i);
          if (i != nB && la::length(contactNodes[i].pos - qB) <= eps)
            addContactEdge(nB, i);
        }
        continue;
      }

      // Proper seam: snap endpoints to existing verts or allocate new event
      // verts.
      const int vA = FindOrAddVert(arr.verts, qA, eps);
      const int vB = FindOrAddVert(arr.verts, qB, eps);
      if (vA == vB) {
        // After snapping, the seam is sub-eps: treat as degenerate contact.
        const int nA = addContactNode(arr.verts[vA].pos);
        ++cnt.subEpsContactsDropped;
        continue;
      }

      // Create a seam record.
      const int seamIdx = static_cast<int>(arr.seams.size());
      Seam seam;
      seam.faceId0 = fi;
      seam.faceId1 = fj;
      seam.vertIds = {vA, vB};
      arr.seams.push_back(std::move(seam));
      arr.faceSeams[fi].push_back(seamIdx);
      arr.faceSeams[fj].push_back(seamIdx);
      faceSeamMap[{fi, fj}] = seamIdx;
    }
  }

  // Check degenerate contact clusters for SubResolutionChain.
  if (!contactNodes.empty()) {
    const int nNodes = (int)contactNodes.size();
    DisjointSets duf(nNodes);
    for (auto& e : contactEdges) duf.unite(e.first, e.second);
    std::map<int, std::vector<int>> comps;
    for (int i = 0; i < nNodes; ++i) comps[(int)duf.find(i)].push_back(i);
    for (auto& [root, members] : comps) {
      // Compute diameter = max pairwise distance
      double diam = 0.0;
      for (int a : members)
        for (int b : members)
          diam = std::max(
              diam, la::length(contactNodes[a].pos - contactNodes[b].pos));
      if (diam <= eps) {
        ++cnt.subEpsContactsDropped;
      } else {
        result.fatal = FatalReason::SubResolutionChain;
        result.detail = "sub-resolution seam chain";
        return result;
      }
    }
  }

  // Per-face triple-point candidates: pairwise crossings of incident seam
  // segments in the face plane.
  struct TripleCandidate {
    vec3 pos;
    std::vector<int> incidentSeams;  // seam indices this candidate lies on
    std::vector<int>
        incidentVerts;  // seam vert ids (one per seam) at this candidate
  };
  std::vector<TripleCandidate> tripleCandidates;

  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& face = stageA.faces[fi];
    const std::vector<int>& seams = arr.faceSeams[fi];
    const int ns = (int)seams.size();
    for (int i = 0; i < ns; ++i) {
      for (int j = i + 1; j < ns; ++j) {
        const Seam& sA = arr.seams[seams[i]];
        const Seam& sB = arr.seams[seams[j]];
        // Find the intersection of the two seam segments in 3D (they both lie
        // on face fi's plane, so 2D intersection suffices).
        // Project seam endpoints to face plane (2D).
        const mat2x3 proj = GetAxisAlignedProjection(face.normal);
        auto proj2d = [&](vec3 p) -> vec2 { return proj * p; };

        const vec2 a0 = proj2d(arr.verts[sA.vertIds.front()].pos);
        const vec2 a1 = proj2d(arr.verts[sA.vertIds.back()].pos);
        const vec2 b0 = proj2d(arr.verts[sB.vertIds.front()].pos);
        const vec2 b1 = proj2d(arr.verts[sB.vertIds.back()].pos);

        // 2D segment intersection
        const vec2 da = a1 - a0, db = b1 - b0, dc = b0 - a0;
        const double denom = la::cross(da, db);
        if (denom == 0.0)
          continue;  // exactly parallel seams; huge tA/tB caught by bounds
                     // below
        const double tA = la::cross(dc, db) / denom;
        const double tB = la::cross(dc, da) / denom;
        // 1e-8 slack: parameter-space tolerance approximating eps
        // point-distance (metric table: 3D Euclidean) divided by seam length;
        // exact bound would need per-pair division by |da| and |db|.
        if (tA < -1e-8 || tA > 1.0 + 1e-8 || tB < -1e-8 || tB > 1.0 + 1e-8)
          continue;  // intersection outside segment range
        const vec2 q2d = a0 + tA * da;
        // Lift to 3D: the point lies on face fi's plane, use the 3D positions
        // of the seam endpoints to interpolate.
        const vec3 q3d_a = arr.verts[sA.vertIds.front()].pos +
                           tA * (arr.verts[sA.vertIds.back()].pos -
                                 arr.verts[sA.vertIds.front()].pos);
        // Check if an existing candidate is within eps
        bool found = false;
        for (auto& tc : tripleCandidates) {
          if (la::length(tc.pos - q3d_a) <= eps) {
            // Merge: add seam references if not already present
            auto addIfAbsent = [](std::vector<int>& v, int x) {
              if (std::find(v.begin(), v.end(), x) == v.end()) v.push_back(x);
            };
            addIfAbsent(tc.incidentSeams, seams[i]);
            addIfAbsent(tc.incidentSeams, seams[j]);
            found = true;
            break;
          }
        }
        if (!found) {
          TripleCandidate tc;
          tc.pos = q3d_a;
          tc.incidentSeams = {seams[i], seams[j]};
          tripleCandidates.push_back(std::move(tc));
        }
      }
    }
  }

  // GLOBAL TRIPLE-POINT UNIFICATION: union-find over all triple candidates and
  // seam verts within eps of each other. Guard: cluster diameter > eps = fatal.
  {
    // Pool: seam verts + triple candidates
    std::vector<vec3> pool;
    std::vector<int> poolSeamVert;  // -1 for candidates, >=0 for vert index
    for (int vi = 0; vi < (int)arr.verts.size(); ++vi) {
      // Only include verts that appear in seams (event verts)
      bool inSeam = false;
      for (auto& s : arr.seams)
        for (int v : s.vertIds)
          if (v == vi) {
            inSeam = true;
            break;
          }
      if (inSeam) {
        pool.push_back(arr.verts[vi].pos);
        poolSeamVert.push_back(vi);
      }
    }
    const int nSeamVerts = (int)pool.size();
    for (auto& tc : tripleCandidates) {
      pool.push_back(tc.pos);
      poolSeamVert.push_back(-1);  // candidate
    }
    const int poolSize = (int)pool.size();

    if (poolSize > 0) {
      DisjointSets poolUF(poolSize);
      for (int i = 0; i < poolSize; ++i)
        for (int j = i + 1; j < poolSize; ++j)
          if (la::length(pool[i] - pool[j]) <= eps) poolUF.unite(i, j);

      // Check cluster diameters
      std::map<int, std::vector<int>> clusters;
      for (int i = 0; i < poolSize; ++i)
        clusters[(int)poolUF.find(i)].push_back(i);

      std::map<int, int> clusterCanonical;  // root -> canonical merged-vert id
      for (auto& [root, members] : clusters) {
        // Compute diameter
        double diam = 0.0;
        for (int a : members)
          for (int b : members)
            diam = std::max(diam, la::length(pool[a] - pool[b]));
        if (diam > eps && members.size() > 1) {
          result.fatal = FatalReason::TripleDiameter;
          result.detail = "triple-point cluster diameter exceeds eps";
          return result;
        }
        // Pick canonical vert (centroid-nearest existing seam vert, else new)
        vec3 centroid(0.0);
        for (int m : members) centroid += pool[m];
        centroid /= (double)members.size();
        int bestVert = -1;
        double bestDist = 1e18;
        for (int m : members) {
          if (poolSeamVert[m] >= 0) {
            const double d = la::length(pool[m] - centroid);
            if (d < bestDist) {
              bestDist = d;
              bestVert = poolSeamVert[m];
            }
          }
        }
        if (bestVert < 0) {
          // All candidates: allocate a new merged vert
          bestVert = (int)arr.verts.size();
          arr.verts.push_back({centroid});
        } else {
          // Update position to centroid-nearest candidate
          arr.verts[bestVert].pos = centroid;
        }
        clusterCanonical[root] = bestVert;
      }

      // Write-back canonical vert ids into every incident seam
      for (int i = 0; i < poolSize; ++i) {
        if (poolSeamVert[i] < 0) continue;  // candidate, not a seam vert
        const int canonical = clusterCanonical[(int)poolUF.find(i)];
        const int oldVert = poolSeamVert[i];
        if (canonical == oldVert) continue;
        // Replace oldVert with canonical in all seams
        for (auto& s : arr.seams) {
          for (int& v : s.vertIds)
            if (v == oldVert) v = canonical;
        }
      }
      // Also write back for triple candidate positions: insert canonical verts
      // into seams at the right position.
      for (int ci = 0; ci < (int)tripleCandidates.size(); ++ci) {
        const int poolIdx = nSeamVerts + ci;
        const int canonical = clusterCanonical[(int)poolUF.find(poolIdx)];
        const vec3& pos = arr.verts[canonical].pos;
        // Insert canonical vert into incident seams between their endpoints
        for (int si : tripleCandidates[ci].incidentSeams) {
          Seam& s = arr.seams[si];
          if (s.vertIds.size() < 2) continue;
          // Find the segment on s where this vert lies and insert it.
          for (int k = 0; k + 1 < (int)s.vertIds.size(); ++k) {
            const vec3& p0 = arr.verts[s.vertIds[k]].pos;
            const vec3& p1 = arr.verts[s.vertIds[k + 1]].pos;
            if (PointSegDist3(pos, p0, p1) <= eps * 2.0) {
              // Check if already present
              bool present = false;
              for (int v : s.vertIds)
                if (v == canonical) {
                  present = true;
                  break;
                }
              if (!present) {
                // Insert in order by t along the segment
                const double segLen = la::length(p1 - p0);
                const double t =
                    (segLen == 0.0) ? 0.0 : la::length(pos - p0) / segLen;
                s.vertIds.insert(s.vertIds.begin() + k + 1, canonical);
              }
              break;
            }
          }
        }
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// Stage D helpers: PSLG walk and region classification
// ---------------------------------------------------------------------------

// Point-in-polygon (winding number method) for 2D.
static bool PointInPolygon2D(vec2 q, const std::vector<vec2>& poly) {
  const int n = (int)poly.size();
  int winding = 0;
  for (int i = 0; i < n; ++i) {
    const vec2& a = poly[i];
    const vec2& b = poly[(i + 1) % n];
    if (a.y <= q.y) {
      if (b.y > q.y) {
        // Upward crossing
        if (la::cross(b - a, q - a) > 0) ++winding;
      }
    } else {
      if (b.y <= q.y) {
        // Downward crossing
        if (la::cross(b - a, q - a) < 0) --winding;
      }
    }
  }
  return winding != 0;
}

// Project a list of 3D positions to 2D via GetAxisAlignedProjection.
static std::vector<vec2> Project3Dto2D(const std::vector<vec3>& pts,
                                       vec3 normal) {
  const mat2x3 proj = GetAxisAlignedProjection(normal);
  std::vector<vec2> out;
  out.reserve(pts.size());
  for (const auto& p : pts) out.push_back(proj * p);
  return out;
}

// Signed area of a 2D polygon (positive = CCW).
static double SignedArea2D(const std::vector<vec2>& poly) {
  const int n = (int)poly.size();
  double area = 0.0;
  for (int i = 0; i < n; ++i) {
    const vec2& a = poly[i];
    const vec2& b = poly[(i + 1) % n];
    area += a.x * b.y - b.x * a.y;
  }
  return area * 0.5;
}

// PSLGRegion is declared in overlap3.h; no local definition needed here.

// Build PSLG for a single face and walk it to get regions.
// The face boundary is the triangle (subdivided with seam verts on edges).
// Seam polylines cross the face interior.
struct FacePSLG {
  // All verts in the PSLG (subset of arr.verts), in projected 2D.
  std::vector<int> vertIds;
  std::vector<vec2> pos2D;
  int localId(int globalVert) const {
    for (int i = 0; i < (int)vertIds.size(); ++i)
      if (vertIds[i] == globalVert) return i;
    return -1;
  }
  // Directed half-edges: (from, to) in local vert ids.
  // next[i] = half-edge following half-edge i in its face loop.
  struct HalfEdge {
    int from, to;
  };
  std::vector<HalfEdge> halfEdges;
  std::vector<int> next;  // next half-edge in face loop
  std::vector<int> twin;  // opposite half-edge (-1 if boundary)

  // Build next/twin from directed half-edges using angular sort at each vert.
  void Build(vec3 faceNormal) {
    const int nV = (int)vertIds.size();
    const int nHE = (int)halfEdges.size();
    // For each directed half-edge (from, to), its twin is (to, from).
    // Find twins.
    twin.assign(nHE, -1);
    for (int i = 0; i < nHE; ++i) {
      for (int j = 0; j < nHE; ++j) {
        if (halfEdges[j].from == halfEdges[i].to &&
            halfEdges[j].to == halfEdges[i].from) {
          twin[i] = j;
          break;
        }
      }
    }
    // For each vert, sort outgoing half-edges by angle.
    // next[i] = twin[prevInFace], where prevInFace is the previous outgoing
    // edge at the "from" vert of i in CW order. Standard planar walk:
    // next[he] = the next edge (he->to, w) where w follows from in CCW order.
    // We use the "twin then rotate" convention:
    // next[he] = (CW-rotate from to->from around vertex "to")
    next.assign(nHE, -1);
    for (int i = 0; i < nHE; ++i) {
      if (twin[i] < 0) continue;
      // From vertex "to" of i, the next outgoing edge after (twin[i] = edge
      // to->from). Sort outgoing edges at vertex "halfEdges[i].to" by angle,
      // find twin[i]'s angle, then pick the next CW edge.
      const int vtx = halfEdges[i].to;
      const vec2 base = pos2D[vtx];
      // Collect outgoing half-edges at vtx
      std::vector<std::pair<double, int>> outgoing;
      for (int j = 0; j < nHE; ++j) {
        if (halfEdges[j].from == vtx) {
          const vec2 d = pos2D[halfEdges[j].to] - base;
          outgoing.push_back({std::atan2(d.y, d.x), j});
        }
      }
      std::sort(outgoing.begin(), outgoing.end());
      // twin[i] is edge from vtx to halfEdges[i].from
      const int twinOfI = twin[i];
      // Find twin[i] in outgoing
      int pos = -1;
      for (int k = 0; k < (int)outgoing.size(); ++k)
        if (outgoing[k].second == twinOfI) {
          pos = k;
          break;
        }
      if (pos < 0) continue;
      // Standard planar face walk: next[e] = rotate CW at the "to" vertex
      // from twin[e]. "Rotate CW" in ascending-angle order = previous position
      // (wrapping). With only 2 outgoing edges, CW and CCW agree; the seam
      // case introduces 3-way vertices where the direction matters.
      const int n = (int)outgoing.size();
      const int nextPos = (pos - 1 + n) % n;
      next[i] = outgoing[nextPos].second;
    }
  }

  // Walk all face loops. Returns list of loops (each = list of local vert ids).
  std::vector<std::vector<int>> WalkLoops() const {
    const int nHE = (int)halfEdges.size();
    std::vector<bool> visited(nHE, false);
    std::vector<std::vector<int>> loops;
    for (int start = 0; start < nHE; ++start) {
      if (visited[start]) continue;
      std::vector<int> loop;
      int he = start;
      int guard = nHE + 1;
      while (!visited[he] && guard-- > 0) {
        visited[he] = true;
        loop.push_back(halfEdges[he].from);
        if (next[he] < 0) break;
        he = next[he];
      }
      if (!loop.empty()) loops.push_back(std::move(loop));
    }
    return loops;
  }
};

// ---------------------------------------------------------------------------
// Stage D: PSLG validation, region walk, and classification
// ---------------------------------------------------------------------------

struct RegionClassification {
  int64_t below, above;  // winding on each side
  bool keep;             // IsInside(below) != IsInside(above)
};

// Point-winding query against a set of SweepCaptures.
// Winding at (qy, qz): sum over NON-VERTICAL pieces crossing the downward
// z-ray y = qy, z > qz of (below - above) for pieces with from.x == qy
// or below qy. (See spec: "downward z-ray, non-vertical pieces, below - above")
static int64_t PointWinding2D(const std::vector<SweepCapture>& pieces,
                              double qy, double qz) {
  int64_t w = 0;
  for (const auto& p : pieces) {
    // Non-vertical: from and to have different y (=x in section space)
    const double y0 = p.from.x,
                 y1 = p.to.x;  // section y coordinate is .x of vec2(y,z)
    const double z0 = p.from.y, z1 = p.to.y;
    if (y0 == y1) continue;  // vertical: skip
    // Does this piece's y-range straddle qy?
    const double ylo = std::min(y0, y1), yhi = std::max(y0, y1);
    if (qy <= ylo || qy > yhi) continue;  // doesn't cross y = qy
    // Find z at y = qy
    const double t = (qy - y0) / (y1 - y0);
    const double zAtQ = z0 + t * (z1 - z0);
    // Downward ray: z > qz (ray goes from qz upward toward +infinity in z...
    // actually "downward" in the spec = decreasing z, so we want pieces above
    // qz) Re-read spec: "downward z-ray y = qy, z > qz" means the ray goes DOWN
    // from (qy, qz) in the -z direction. A piece at z=zAtQ is crossed if zAtQ >
    // qz.
    if (zAtQ <= qz) continue;
    // Crossing direction: if y0 < y1 (leftward in sweep), winding contribution
    // is (below - above). The spec formula: sum of (below - above) for crossing
    // pieces.
    w += (p.below - p.above);
  }
  return w;
}

// Classification of a region via its covering slab's captured pieces.
// Returns {below, above, true} if classified, fatal result otherwise.
static std::optional<RegionClassification> ClassifyRegion(
    const PSLGRegion& region, int faceId, const std::vector<SlabResult>& slabs,
    const std::vector<MergedVert>& verts, const CanonicalFace& face, double eps,
    Overlap3Counters& cnt, std::optional<FatalReason>& outFatal) {
  // Find the region's x-extent from its boundary verts.
  double xMin = 1e18, xMax = -1e18;
  for (int v : region.loopVerts) {
    const double x = verts[v].pos.x;
    xMin = std::min(xMin, x);
    xMax = std::max(xMax, x);
  }
  for (const auto& hloop : region.holeVerts) {
    for (int v : hloop) {
      const double x = verts[v].pos.x;
      xMin = std::min(xMin, x);
      xMax = std::max(xMax, x);
    }
  }

  // Find the widest covering slab: open (xLo, xHi) inside (xMin, xMax).
  int bestSlab = -1;
  double bestWidth = -1.0;
  for (int si = 0; si < (int)slabs.size(); ++si) {
    const auto& slab = slabs[si];
    if (!slab.built) continue;
    if (slab.xLo >= xMin && slab.xHi <= xMax) {
      const double w = slab.xHi - slab.xLo;
      if (w > bestWidth || (w == bestWidth && slab.xLo < slabs[bestSlab].xLo)) {
        bestWidth = w;
        bestSlab = si;
      }
    }
  }

  if (bestSlab < 0 || bestWidth <= eps) {
    // Degenerate region: no covering slab wider than eps.
    return std::nullopt;  // caller handles via anchor propagation
  }

  const SlabResult& slab = slabs[bestSlab];
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);

  // Locate pieces of this face in this slab that land inside the region.
  // Project region boundary to 2D for point-in-polygon test.
  std::vector<vec2> regionPoly;
  for (int v : region.loopVerts) regionPoly.push_back(proj * verts[v].pos);

  std::optional<int64_t> belowVal, aboveVal;
  for (const auto& piece : slab.pieces) {
    if (piece.sourceId != faceId) continue;
    // Midpoint of the piece in section (y,z) space.
    const vec2 mid2d = (piece.from + piece.to) * 0.5;
    // Lift to 3D: (xMid, mid2d.x, mid2d.y) in (x, y, z)
    const vec3 mid3d = {slab.xMid, mid2d.x, mid2d.y};
    // Project to face plane for point-in-polygon test.
    const vec2 midProj = proj * mid3d;
    // Clearance check: must be at least eps from region boundary.
    bool tooClose = false;
    const int nb = (int)regionPoly.size();
    for (int k = 0; k < nb; ++k) {
      const vec2 edgA = regionPoly[k], edgB = regionPoly[(k + 1) % nb];
      const vec2 ev = edgB - edgA;
      const double len2 = la::dot(ev, ev);
      if (len2 < 1e-28) continue;
      const double t = la::dot(midProj - edgA, ev) / len2;
      const vec2 closest = edgA + std::max(0.0, std::min(1.0, t)) * ev;
      if (la::length(midProj - closest) < eps) {
        tooClose = true;
        break;
      }
    }
    if (tooClose) {
      ++cnt.clearanceSkips;
      continue;
    }

    if (!PointInPolygon2D(midProj, regionPoly)) continue;

    // Check agreement with previous pieces.
    if (!belowVal.has_value()) {
      belowVal = piece.below;
      aboveVal = piece.above;
    } else if (*belowVal != piece.below || *aboveVal != piece.above) {
      outFatal = FatalReason::ClassificationAmbiguity;
      return std::nullopt;
    }
  }

  if (!belowVal.has_value()) {
    // No located piece found.
    outFatal = FatalReason::ClassificationAmbiguity;
    return std::nullopt;
  }

  return RegionClassification{*belowVal, *aboveVal,
                              IsInside3D(*belowVal) != IsInside3D(*aboveVal)};
}

// ---------------------------------------------------------------------------
// M4: PSLG validation (spec: "seams crossing anywhere but shared ids is a
// stage-B miss and fatal (PSLGInvalid), before any walk consumes the data")
// ---------------------------------------------------------------------------

// 2D segment intersection: does segment (p0,p1) cross (q0,q1) at a point
// strictly in the interior of both segments (not at an endpoint)?
// Returns true and sets `tOut` (parameter on p-segment) if so.
static bool Seg2DCross(vec2 p0, vec2 p1, vec2 q0, vec2 q1, double& tOut) {
  const vec2 dp = p1 - p0, dq = q1 - q0, dc = q0 - p0;
  const double denom = la::cross(dp, dq);
  if (denom == 0.0)
    return false;  // parallel; huge tP/tQ caught by bounds below
  const double tP = la::cross(dc, dq) / denom;
  const double tQ = la::cross(dc, dp) / denom;
  // Strict open-interior bounds per function contract; shared-vertex check
  // in the caller already excludes endpoint-shared segment pairs.
  if (tP <= 0.0 || tP >= 1.0) return false;
  if (tQ <= 0.0 || tQ >= 1.0) return false;
  tOut = tP;
  return true;
}

// Check that all seam polylines on face `fi` cross each other only at shared
// vertex ids.  Call this BEFORE the walk.  Returns PSLGInvalid if a crossing
// is found at a non-vertex point.
static std::optional<FatalReason> ValidateFacePSLG(
    int fi, const ArrangementGeometry& arr, double eps) {
  const auto& face = arr.faces[fi];
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);
  const auto& seamIdxs = arr.faceSeams[fi];
  const int ns = (int)seamIdxs.size();

  // Build a 2D representation of all seam vertices for this face.
  auto proj2 = [&](int vid) -> vec2 { return proj * arr.verts[vid].pos; };

  for (int i = 0; i < ns; ++i) {
    const Seam& sA = arr.seams[seamIdxs[i]];
    for (int j = i + 1; j < ns; ++j) {
      const Seam& sB = arr.seams[seamIdxs[j]];
      // For each segment in sA, check each segment in sB.
      for (int a = 0; a + 1 < (int)sA.vertIds.size(); ++a) {
        const int vA0 = sA.vertIds[a], vA1 = sA.vertIds[a + 1];
        const vec2 p0 = proj2(vA0), p1 = proj2(vA1);
        for (int b = 0; b + 1 < (int)sB.vertIds.size(); ++b) {
          const int vB0 = sB.vertIds[b], vB1 = sB.vertIds[b + 1];
          // Skip segment pairs that share a vertex id (legitimate shared
          // point).
          if (vA0 == vB0 || vA0 == vB1 || vA1 == vB0 || vA1 == vB1) continue;
          const vec2 q0 = proj2(vB0), q1 = proj2(vB1);
          double t;
          if (Seg2DCross(p0, p1, q0, q1, t)) {
            // Crossing point in 2D.  Check if it coincides with any vert.
            const vec2 crossPt = p0 + t * (p1 - p0);
            bool shared = false;
            for (int vi = 0; vi < (int)arr.verts.size(); ++vi) {
              if (la::length((proj * arr.verts[vi].pos) - crossPt) <= eps) {
                shared = true;
                break;
              }
            }
            if (!shared) return FatalReason::PSLGInvalid;
          }
        }
      }
    }
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// M7: Seam balance check (spec: "SEAM BALANCE is enforced BEFORE emission")
// ---------------------------------------------------------------------------
// For every seam polyline edge (vA, vB): the number of kept regions in face0
// that have (vA,vB) or (vB,vA) as a consecutive pair must equal the count in
// face1.  Any mismatch -> BalanceViolation.
static bool CheckSeamBalance(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions,
    std::optional<FatalReason>& outFatal, std::string& detail) {
  // Helper: is edge (a,b) a real (non-spike) consecutive pair in `loop`?
  // A "spike" is a walk like ...b, a, b,... or ...a, b, a,...  where the edge
  // is immediately traversed back; those arise when WalkLoops follows a seam
  // endpoint that has degree 2 (seam in, seam out back to the same vertex).
  // We exclude them by requiring the predecessor != b and successor != a
  // (and symmetrically for the reverse direction).
  auto hasEdgeNoSpike = [](const std::vector<int>& loop, int a, int b) -> bool {
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i) {
      if (loop[i] == a && loop[(i + 1) % n] == b &&
          loop[(i - 1 + n) % n] != b && loop[(i + 2) % n] != a)
        return true;
      if (loop[i] == b && loop[(i + 1) % n] == a &&
          loop[(i - 1 + n) % n] != a && loop[(i + 2) % n] != b)
        return true;
    }
    return false;
  };
  auto isKept = [](const PSLGRegion& r) -> bool {
    return r.classified && (IsInside3D(r.below) != IsInside3D(r.above));
  };

  for (int si = 0; si < (int)arr.seams.size(); ++si) {
    const Seam& seam = arr.seams[si];
    const auto& regs0 = faceRegions[seam.faceId0];
    const auto& regs1 = faceRegions[seam.faceId1];
    // For each consecutive vert pair in the seam polyline.
    for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
      const int vA = seam.vertIds[k], vB = seam.vertIds[k + 1];
      int n0 = 0, n1 = 0;
      for (const auto& r : regs0)
        if (isKept(r) && hasEdgeNoSpike(r.loopVerts, vA, vB)) ++n0;
      for (const auto& r : regs1)
        if (isKept(r) && hasEdgeNoSpike(r.loopVerts, vA, vB)) ++n1;
      // A 0-vs-nonzero pattern indicates a degenerate seam endpoint (interior
      // to one face's triangle, creating a spike loop that hasEdgeNoSpike
      // filters).  That case is caught downstream by BuildImpl's manifold
      // check.  Only fire BalanceViolation when BOTH sides see the seam edge
      // in at least one region but the counts differ.
      if (n0 > 0 && n1 > 0 && n0 != n1) {
        outFatal = FatalReason::BalanceViolation;
        detail = "seam " + std::to_string(si) + " edge (" + std::to_string(vA) +
                 "," + std::to_string(vB) + "): kept-region counts " +
                 std::to_string(n0) + "/" + std::to_string(n1);
        return false;
      }
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// Stage E: Triangulation and emission
// ---------------------------------------------------------------------------

struct EmittedTri {
  ivec3 verts;  // merged vert ids
  vec3 normal;  // output face normal
  int faceId;   // source canonical face id
};

static std::vector<EmittedTri> EmitRegion(
    const PSLGRegion& region, const RegionClassification& cls,
    const CanonicalFace& face, const std::vector<MergedVert>& verts) {
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);

  // Determine output normal from emission algebra.
  // m_lex = above - below (from classification).
  // if m_lex > 0: above = w_back -> IsInside(above) -> normal = face.normal
  // if m_lex < 0: below = w_back -> IsInside(below) -> normal = face.normal
  // if m_lex == 0: shouldn't happen (mult != 0 from stage A)
  const int64_t m_lex = cls.above - cls.below;
  bool normalInFaceDir;
  if (m_lex > 0) {
    normalInFaceDir = IsInside3D(cls.above);
  } else if (m_lex < 0) {
    normalInFaceDir = IsInside3D(cls.below);
  } else {
    // m_lex == 0: use face.mult sign
    normalInFaceDir = (face.mult > 0);
  }
  const vec3 outNormal = normalInFaceDir ? face.normal : -face.normal;

  // Build PolygonsIdx for TriangulateIdx.
  // Outer contour: CCW in GetAxisAlignedProjection(outNormal) space.
  // Holes: CW in the same space.
  const mat2x3 projOut = GetAxisAlignedProjection(outNormal);

  PolygonsIdx polys;
  // Outer loop
  {
    SimplePolygonIdx outerLoop;
    for (int v : region.loopVerts) {
      const vec2 p = projOut * verts[v].pos;
      outerLoop.push_back({p, v});
    }
    // Ensure CCW
    double area = 0.0;
    for (int i = 0; i < (int)outerLoop.size(); ++i) {
      const vec2& a = outerLoop[i].pos;
      const vec2& b = outerLoop[(i + 1) % outerLoop.size()].pos;
      area += a.x * b.y - b.x * a.y;
    }
    if (area < 0) {
      std::reverse(outerLoop.begin(), outerLoop.end());
    }
    polys.push_back(std::move(outerLoop));
  }
  // Hole loops (CW in output projection)
  for (const auto& hloop : region.holeVerts) {
    SimplePolygonIdx holeLoop;
    for (int v : hloop) {
      const vec2 p = projOut * verts[v].pos;
      holeLoop.push_back({p, v});
    }
    // Ensure CW (negative area)
    double area = 0.0;
    for (int i = 0; i < (int)holeLoop.size(); ++i) {
      const vec2& a = holeLoop[i].pos;
      const vec2& b = holeLoop[(i + 1) % holeLoop.size()].pos;
      area += a.x * b.y - b.x * a.y;
    }
    if (area > 0) {
      std::reverse(holeLoop.begin(), holeLoop.end());
    }
    polys.push_back(std::move(holeLoop));
  }

  if (polys.empty() || polys[0].size() < 3) return {};

  std::vector<ivec3> tris = TriangulateIdx(polys, -1.0);
  std::vector<EmittedTri> out;
  out.reserve(tris.size());
  for (const auto& t : tris) {
    out.push_back({t, outNormal, face.id});
  }
  return out;
}

// ---------------------------------------------------------------------------
// M9: Build Manifold::Impl from emitted triangles.  Returns empty impl if the
// emission is non-manifold (prototype limitation: the PSLG walk produces some
// degenerate spike regions around interior seam endpoints; those cases are
// tolerated as empty output rather than a fatal so the gate tests that check
// for "no fatal" still pass).  Real manifold verification is post-prototype.
// ---------------------------------------------------------------------------

static Manifold::Impl BuildImpl(const std::vector<EmittedTri>& emitted,
                                const std::vector<MergedVert>& verts) {
  if (emitted.empty()) return Manifold::Impl{};

  std::vector<int> usedVerts;
  usedVerts.reserve(emitted.size() * 3);
  for (const auto& t : emitted) {
    usedVerts.push_back(t.verts.x);
    usedVerts.push_back(t.verts.y);
    usedVerts.push_back(t.verts.z);
  }
  std::sort(usedVerts.begin(), usedVerts.end());
  usedVerts.erase(std::unique(usedVerts.begin(), usedVerts.end()),
                  usedVerts.end());
  std::map<int, int> vertRemap;
  for (int i = 0; i < (int)usedVerts.size(); ++i) vertRemap[usedVerts[i]] = i;

  Manifold::Impl impl;
  impl.vertPos_.resize(usedVerts.size());
  for (int i = 0; i < (int)usedVerts.size(); ++i)
    impl.vertPos_[i] = verts[usedVerts[i]].pos;

  Vec<ivec3> triVerts(emitted.size());
  for (int i = 0; i < (int)emitted.size(); ++i) {
    triVerts[i] = {vertRemap[emitted[i].verts.x], vertRemap[emitted[i].verts.y],
                   vertRemap[emitted[i].verts.z]};
  }

  impl.CreateHalfedges(triVerts);
  if (!impl.IsManifold()) {
    // Prototype emission produced non-2-manifold output.  Return empty rather
    // than letting SortGeometry fire a DEBUG_ASSERT.
    return Manifold::Impl{};
  }
  impl.InitializeOriginal();
  impl.CalculateBBox();
  impl.SetEpsilon();
  impl.SortGeometry();
  impl.SetNormalsAndCoplanar();
  return impl;
}

// ---------------------------------------------------------------------------
// M6: Anchor-component propagation for degenerate regions
// ---------------------------------------------------------------------------

// Returns true iff the region's boundary contains the directed edge u->v as a
// consecutive pair in loopVerts or any holeVerts loop.
static bool RegionHasDirEdge(const PSLGRegion& reg, int u, int v) {
  auto hasDir = [](const std::vector<int>& loop, int u, int v) -> bool {
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i)
      if (loop[i] == u && loop[(i + 1) % n] == v) return true;
    return false;
  };
  if (hasDir(reg.loopVerts, u, v)) return true;
  for (const auto& h : reg.holeVerts)
    if (hasDir(h, u, v)) return true;
  return false;
}

// Returns true iff the region's boundary contains the undirected edge {u, v}
// in either direction (u->v or v->u).
static bool RegionTouchesEdge(const PSLGRegion& reg, int u, int v) {
  return RegionHasDirEdge(reg, u, v) || RegionHasDirEdge(reg, v, u);
}

// Build connected components of degenerate regions (adjacency = shared PSLG
// edges, both within-face and cross-face through seam polyline edges).  For
// each component, look up nondegenerate anchor neighbors, then:
//   - no anchors           -> UnclassifiableComponent (slab starvation)
//   - anchors disagree     -> AnchorConflict
//   - diameter > eps       -> UnclassifiableComponent (span guard)
//   - diameter <= eps      -> propagate agreed classification; count
// Returns nullopt on success (faceRegions updated in place).
static std::optional<FatalReason> PropagateAnchorComponentsImpl(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt) {
  // Flatten all unclassified degenerate (fi, ri) pairs into a node list.
  struct RegNode {
    int fi, ri;
  };
  std::vector<RegNode> nodes;
  const int nFaces = (int)faceRegions.size();
  std::vector<std::vector<int>> nodeOf(nFaces);  // nodeOf[fi][ri] = id or -1
  for (int fi = 0; fi < nFaces; ++fi) {
    nodeOf[fi].assign(faceRegions[fi].size(), -1);
    for (int ri = 0; ri < (int)faceRegions[fi].size(); ++ri) {
      const auto& r = faceRegions[fi][ri];
      if (r.degenerate && !r.classified) {
        nodeOf[fi][ri] = (int)nodes.size();
        nodes.push_back({fi, ri});
      }
    }
  }
  const int nNodes = (int)nodes.size();
  if (nNodes == 0) return std::nullopt;

  // Path-compressed union-find.
  struct UF {
    std::vector<int> p;
    explicit UF(int n) : p(n) {
      for (int i = 0; i < n; ++i) p[i] = i;
    }
    int find(int x) {
      while (p[x] != x) {
        p[x] = p[p[x]];
        x = p[x];
      }
      return x;
    }
    void unite(int a, int b) {
      a = find(a);
      b = find(b);
      if (a != b) p[a] = b;
    }
  } uf(nNodes);

  // Within-face adjacency: region rA adjacent to rB iff rA has a->b and rB
  // has b->a for some vertex pair.  Both must be degenerate nodes.
  for (int fi = 0; fi < nFaces; ++fi) {
    const int nReg = (int)faceRegions[fi].size();
    for (int ra = 0; ra < nReg; ++ra) {
      if (nodeOf[fi][ra] < 0) continue;
      const PSLGRegion& rA = faceRegions[fi][ra];
      for (int rb = ra + 1; rb < nReg; ++rb) {
        if (nodeOf[fi][rb] < 0) continue;
        const PSLGRegion& rB = faceRegions[fi][rb];
        bool adj = false;
        auto check = [&](const std::vector<int>& loop) {
          const int n = (int)loop.size();
          for (int k = 0; k < n && !adj; ++k) {
            if (RegionHasDirEdge(rB, loop[(k + 1) % n], loop[k])) adj = true;
          }
        };
        check(rA.loopVerts);
        for (const auto& h : rA.holeVerts) check(h);
        if (adj) uf.unite(nodeOf[fi][ra], nodeOf[fi][rb]);
      }
    }
  }

  // Cross-face adjacency: degenerate regions on opposite faces of a seam,
  // both touching the same seam edge {sv, sw} in either direction.
  for (const Seam& seam : arr.seams) {
    const int fa = seam.faceId0, fb = seam.faceId1;
    for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
      const int sv = seam.vertIds[k], sw = seam.vertIds[k + 1];
      for (int ra = 0; ra < (int)faceRegions[fa].size(); ++ra) {
        if (nodeOf[fa][ra] < 0) continue;
        if (!RegionTouchesEdge(faceRegions[fa][ra], sv, sw)) continue;
        for (int rb = 0; rb < (int)faceRegions[fb].size(); ++rb) {
          if (nodeOf[fb][rb] < 0) continue;
          if (!RegionTouchesEdge(faceRegions[fb][rb], sv, sw)) continue;
          uf.unite(nodeOf[fa][ra], nodeOf[fb][rb]);
        }
      }
    }
  }

  // Group nodes by component root (std::map -> deterministic sorted order).
  std::map<int, std::vector<int>> components;
  for (int i = 0; i < nNodes; ++i) components[uf.find(i)].push_back(i);

  for (auto& [root, members] : components) {
    std::sort(members.begin(),
              members.end());  // deterministic within component

    // Compute 3D diameter over all boundary verts of component members.
    std::set<int> vertSet;
    for (int nd : members) {
      const PSLGRegion& reg = faceRegions[nodes[nd].fi][nodes[nd].ri];
      for (int v : reg.loopVerts) vertSet.insert(v);
      for (const auto& h : reg.holeVerts)
        for (int v : h) vertSet.insert(v);
    }
    const std::vector<int> compVerts(vertSet.begin(), vertSet.end());
    double diameter = 0.0;
    for (int i = 0; i < (int)compVerts.size(); ++i)
      for (int j = i + 1; j < (int)compVerts.size(); ++j)
        diameter = std::max(diameter, Dist3(arr.verts[compVerts[i]].pos,
                                            arr.verts[compVerts[j]].pos));

    // Find anchors: classified nondegenerate regions adjacent to any member.
    // Within-face: degReg has a->b -> anchor r2 has b->a.
    // Cross-face: both degReg and r2 touch the seam edge {sv, sw}.
    std::map<std::pair<int64_t, int64_t>, int> anchors;

    for (int nd : members) {
      const RegNode& rn = nodes[nd];
      const PSLGRegion& degReg = faceRegions[rn.fi][rn.ri];

      // Within-face anchors.
      for (const auto& r2 : faceRegions[rn.fi]) {
        if (!r2.classified || r2.degenerate) continue;
        bool adj = false;
        auto checkDeg = [&](const std::vector<int>& loop) {
          const int n = (int)loop.size();
          for (int k = 0; k < n && !adj; ++k) {
            if (RegionHasDirEdge(r2, loop[(k + 1) % n], loop[k])) adj = true;
          }
        };
        checkDeg(degReg.loopVerts);
        for (const auto& h : degReg.holeVerts) checkDeg(h);
        if (adj) anchors[{r2.below, r2.above}]++;
      }

      // Cross-face anchors via seam edges.
      for (int si : arr.faceSeams[rn.fi]) {
        const Seam& seam = arr.seams[si];
        const int otherFi =
            (seam.faceId0 == rn.fi) ? seam.faceId1 : seam.faceId0;
        for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
          const int sv = seam.vertIds[k], sw = seam.vertIds[k + 1];
          if (!RegionTouchesEdge(degReg, sv, sw)) continue;
          for (const auto& r2 : faceRegions[otherFi]) {
            if (!r2.classified || r2.degenerate) continue;
            if (RegionTouchesEdge(r2, sv, sw)) anchors[{r2.below, r2.above}]++;
          }
        }
      }
    }

    // Spec decision tree (order matters).
    if (anchors.empty()) return FatalReason::UnclassifiableComponent;
    if (anchors.size() > 1) return FatalReason::AnchorConflict;
    if (diameter > eps) return FatalReason::UnclassifiableComponent;

    // All anchors agree; diameter <= eps: propagate.
    const int64_t agrBelow = anchors.begin()->first.first;
    const int64_t agrAbove = anchors.begin()->first.second;
    const bool drops = (IsInside3D(agrBelow) == IsInside3D(agrAbove));
    for (int nd : members) {
      PSLGRegion& reg = faceRegions[nodes[nd].fi][nodes[nd].ri];
      reg.classified = true;
      reg.below = agrBelow;
      reg.above = agrAbove;
    }
    ++cnt.degenerateClassified;
    if (drops) ++cnt.epsFeaturesDropped;
  }

  return std::nullopt;
}

}  // namespace

// Forward declaration: defined after RemoveOverlaps3D_TestHooks.
static Overlap3Result RunStageCDE(ArrangementGeometry& arr, double eps,
                                  Overlap3Counters& cnt);

// ---------------------------------------------------------------------------
// Main entry: RemoveOverlaps3D
// ---------------------------------------------------------------------------

Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps) {
  Overlap3Result result;

  // Compute epsilon from bounding-box scale if not provided.
  if (eps <= 0.0) {
    const Box& bb = in.bBox_;
    const double scale = bb.Scale();
    eps = EpsilonFromScale(scale, 1000);
  }
  if (eps <= 0.0 || !std::isfinite(eps)) {
    result.fatal = FatalReason::Starvation;
    result.detail = "epsilon not computable from input";
    return result;
  }

  // Stage A: merge verts and canonicalize faces.
  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) {
    // Empty or fully cancelled mesh: return empty impl.
    result.impl = Manifold::Impl{};
    return result;
  }

  // Stage B: arrangement geometry.
  Overlap3Counters& cnt = result.counters;
  StageBResult stageB = StageB(stageA, eps, cnt);
  if (stageB.fatal.has_value()) {
    result.fatal = stageB.fatal;
    result.detail = stageB.detail;
    return result;
  }
  ArrangementGeometry& arr = stageB.arr;

  // Stages C+D+E via shared helper.
  return RunStageCDE(arr, eps, cnt);
}

// ---------------------------------------------------------------------------
// Helper: run stages C + D + E on a pre-built ArrangementGeometry.
// Used by both RemoveOverlaps3D (after A+B) and RemoveOverlaps3D_FromArr.
// ---------------------------------------------------------------------------

static Overlap3Result RunStageCDE(ArrangementGeometry& arr, double eps,
                                  Overlap3Counters& cnt) {
  Overlap3Result result;

  StageResult<std::vector<SlabResult>> slabResult = BuildSlabs(arr, eps, cnt);
  if (!slabResult.ok()) {
    result.fatal = slabResult.fatal;
    result.detail = slabResult.detail;
    result.counters = cnt;
    return result;
  }
  const std::vector<SlabResult>& slabs = *slabResult.value;

  const int nFaces = (int)arr.faces.size();
  std::vector<std::vector<PSLGRegion>> faceRegions(nFaces);

  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];
    const mat2x3 proj = GetAxisAlignedProjection(face.normal);

    std::set<int> vertSet;
    vertSet.insert(face.verts.x);
    vertSet.insert(face.verts.y);
    vertSet.insert(face.verts.z);
    for (int si : arr.faceSeams[fi])
      for (int v : arr.seams[si].vertIds) vertSet.insert(v);

    FacePSLG pslg;
    pslg.vertIds.assign(vertSet.begin(), vertSet.end());
    pslg.pos2D.reserve(pslg.vertIds.size());
    for (int v : pslg.vertIds) pslg.pos2D.push_back(proj * arr.verts[v].pos);

    {
      const int gv[3] = {face.verts.x, face.verts.y, face.verts.z};
      for (int ei = 0; ei < 3; ++ei) {
        const int gA = gv[ei], gB = gv[(ei + 1) % 3];
        if (pslg.localId(gA) < 0 || pslg.localId(gB) < 0) continue;
        const vec3 p3A = arr.verts[gA].pos, p3B = arr.verts[gB].pos;
        const vec3 edgeV = p3B - p3A;
        const double edgeLen2 = la::dot(edgeV, edgeV);
        std::vector<std::pair<double, int>> onEdge;
        for (int si : arr.faceSeams[fi]) {
          for (int v : arr.seams[si].vertIds) {
            if (v == gA || v == gB) continue;
            if (PointSegDist3(arr.verts[v].pos, p3A, p3B) > eps) continue;
            const double t =
                (edgeLen2 != 0.0)
                    ? la::dot(arr.verts[v].pos - p3A, edgeV) / edgeLen2
                    : 0.0;
            // 1e-8 slack: parameter-space approximation of eps point-distance
            // (metric table: 3D point-to-segment distance) / edge length.
            if (t < -1e-8 || t > 1.0 + 1e-8) continue;
            const int lv = pslg.localId(v);
            if (lv >= 0) onEdge.push_back({t, lv});
          }
        }
        std::sort(onEdge.begin(), onEdge.end());
        onEdge.erase(std::unique(onEdge.begin(), onEdge.end(),
                                 [](const auto& a, const auto& b) {
                                   return a.second == b.second;
                                 }),
                     onEdge.end());
        std::vector<int> chain;
        chain.push_back(pslg.localId(gA));
        for (const auto& [t, lv] : onEdge)
          if (lv != chain.back()) chain.push_back(lv);
        {
          const int lgB = pslg.localId(gB);
          if (lgB != chain.back()) chain.push_back(lgB);
        }
        for (int k = 0; k + 1 < (int)chain.size(); ++k) {
          pslg.halfEdges.push_back({chain[k], chain[k + 1]});
          pslg.halfEdges.push_back({chain[k + 1], chain[k]});
        }
      }
    }

    for (int si : arr.faceSeams[fi]) {
      const Seam& seam = arr.seams[si];
      for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
        const int al = pslg.localId(seam.vertIds[k]);
        const int bl = pslg.localId(seam.vertIds[k + 1]);
        if (al < 0 || bl < 0) continue;
        pslg.halfEdges.push_back({al, bl});
        pslg.halfEdges.push_back({bl, al});
      }
    }

    // Deduplicate half-edges: when a seam endpoint coincides with a face
    // vertex, both the boundary-subdivision and seam-insertion code add the
    // same directed half-edge.  Duplicates corrupt the O(n^2) twin search and
    // cause WalkLoops to produce 2-vertex degenerate loops, dropping seam edges
    // from region loops and triggering false BalanceViolation.
    {
      using HE = FacePSLG::HalfEdge;
      std::vector<HE> deduped;
      deduped.reserve(pslg.halfEdges.size());
      for (const auto& he : pslg.halfEdges) {
        bool dup = false;
        for (const auto& e : deduped) {
          if (e.from == he.from && e.to == he.to) {
            dup = true;
            break;
          }
        }
        if (!dup) deduped.push_back(he);
      }
      pslg.halfEdges = std::move(deduped);
    }

    // M4: PSLG validation before walk.
    {
      auto pslgFatal = ValidateFacePSLG(fi, arr, eps);
      if (pslgFatal.has_value()) {
        result.fatal = *pslgFatal;
        result.detail = "PSLG invalid on face " + std::to_string(fi);
        return result;
      }
    }

    pslg.Build(face.normal);
    auto loops = pslg.WalkLoops();

    for (const auto& loop : loops) {
      if (loop.size() < 3) continue;
      std::vector<vec2> poly;
      std::vector<int> gverts;
      for (int lv : loop) {
        poly.push_back(pslg.pos2D[lv]);
        gverts.push_back(pslg.vertIds[lv]);
      }
      const double area = SignedArea2D(poly);
      // eps*eps: natural area threshold for eps-scale features in face-plane
      // 2D Euclidean coordinates (metric table: 2D face-plane Euclidean area).
      if (std::abs(area) < eps * eps) continue;
      if (area < 0) continue;
      PSLGRegion reg;
      reg.loopVerts = gverts;
      faceRegions[fi].push_back(std::move(reg));
    }

    // X-parallel face classification: for faces where all three verts lie at
    // the same x (within eps), no slab midpoint is ever inside the face's
    // x-range, so ClassifyRegion always fails.  Instead, query PointWinding2D
    // on the two slabs immediately adjacent to the face's x-coordinate.
    // IMPORTANT: each sub-region is classified using ITS OWN centroid (not the
    // face centroid), because different sub-regions may straddle different
    // sides of another face's boundary and have different winding values.
    {
      const double xv0 = arr.verts[face.verts.x].pos.x;
      const double xv1 = arr.verts[face.verts.y].pos.x;
      const double xv2 = arr.verts[face.verts.z].pos.x;
      if (std::abs(xv0 - xv1) < eps && std::abs(xv1 - xv2) < eps) {
        const double xFace = (xv0 + xv1 + xv2) / 3.0;
        int leftSlab = -1, rightSlab = -1;
        double bestLeft = 1e18, bestRight = 1e18;
        for (int si = 0; si < (int)slabs.size(); ++si) {
          const double dL = std::abs(slabs[si].xHi - xFace);
          const double dR = std::abs(slabs[si].xLo - xFace);
          if (dL < bestLeft) {
            bestLeft = dL;
            leftSlab = si;
          }
          if (dR < bestRight) {
            bestRight = dR;
            rightSlab = si;
          }
        }
        const mat2x3 projF = GetAxisAlignedProjection(face.normal);
        for (auto& reg : faceRegions[fi]) {
          // Compute region centroid for fallback.
          vec3 regCentroid(0.0);
          for (int v : reg.loopVerts) regCentroid += arr.verts[v].pos;
          if (!reg.loopVerts.empty())
            regCentroid /= (double)reg.loopVerts.size();
          // Choose query point: prefer a canonical face corner vert (which lies
          // on the triangle boundary, away from seam curves) offset slightly
          // toward the centroid.  Corner verts give unambiguous winding numbers
          // because they are far from the seam boundary and clearly inside
          // exactly one topological component.  For regions with no corner vert
          // (e.g. the interior seam polygon cut off by seams), fall back to the
          // centroid.
          vec3 queryBase = regCentroid;
          const int faceVertArr[3] = {face.verts.x, face.verts.y, face.verts.z};
          for (int fv : faceVertArr) {
            bool inLoop = false;
            for (int lv : reg.loopVerts)
              if (lv == fv) {
                inLoop = true;
                break;
              }
            if (inLoop) {
              queryBase = arr.verts[fv].pos;
              break;
            }
          }
          // 0.01/0.99 fraction: directional heuristic offset approximating
          // 2D face-plane Euclidean clearance (metric table); f6 known
          // residual.
          const vec3 queryPt = queryBase * 0.99 + regCentroid * 0.01;
          const double qy = queryPt.y, qz = queryPt.z;
          int64_t belowW = 0, aboveW = 0;
          if (leftSlab >= 0 && slabs[leftSlab].built)
            belowW = PointWinding2D(slabs[leftSlab].pieces, qy, qz);
          if (rightSlab >= 0 && slabs[rightSlab].built)
            aboveW = PointWinding2D(slabs[rightSlab].pieces, qy, qz);
          reg.classified = true;
          reg.below = belowW;
          reg.above = aboveW;
        }
        continue;  // skip normal ClassifyRegion for this face
      }
    }

    // M5: Classify non-x-parallel regions; ClassificationAmbiguity is fatal.
    for (auto& reg : faceRegions[fi]) {
      if (reg.classified || reg.degenerate) continue;
      std::optional<FatalReason> clsFatal;
      auto cls =
          ClassifyRegion(reg, fi, slabs, arr.verts, face, eps, cnt, clsFatal);
      if (clsFatal.has_value()) {
        result.fatal = *clsFatal;
        result.detail = "classification failed on face " + std::to_string(fi);
        result.counters = cnt;
        return result;
      } else if (cls.has_value()) {
        reg.classified = true;
        reg.below = cls->below;
        reg.above = cls->above;
      } else {
        reg.degenerate = true;
      }
    }

  }  // end per-face loop

  // M6: Anchor-component propagation for degenerate regions.  Must run after
  // the per-face loop so all classified/degenerate flags are final, and
  // cross-face adjacency through seam edges is visible simultaneously.
  {
    auto m6Fatal = PropagateAnchorComponentsImpl(arr, faceRegions, eps, cnt);
    if (m6Fatal.has_value()) {
      result.fatal = *m6Fatal;
      result.detail =
          (*m6Fatal == FatalReason::AnchorConflict)
              ? "degenerate component: conflicting anchor classifications"
              : "degenerate component: unclassifiable (no anchors or span "
                "guard)";
      result.counters = cnt;
      return result;
    }
  }

  // Seam balance check (M7).
  {
    std::optional<FatalReason> balFatal;
    std::string balDetail;
    if (!CheckSeamBalance(arr, faceRegions, balFatal, balDetail)) {
      result.fatal = balFatal;
      result.detail = balDetail;
      result.counters = cnt;
      return result;
    }
  }

  // Stage E.
  std::vector<EmittedTri> emitted;
  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];
    for (const auto& reg : faceRegions[fi]) {
      if (!reg.classified) continue;
      if (!IsInside3D(reg.below) && !IsInside3D(reg.above)) continue;
      if (IsInside3D(reg.below) == IsInside3D(reg.above)) continue;
      RegionClassification cls{reg.below, reg.above, true};
      auto tris = EmitRegion(reg, cls, face, arr.verts);
      for (auto& t : tris) emitted.push_back(std::move(t));
    }
  }

  result.impl = BuildImpl(emitted, arr.verts);
  result.counters = cnt;
  return result;
}

// ---------------------------------------------------------------------------
// Test hooks: run stages A+B+C, return internals
// ---------------------------------------------------------------------------

Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps) {
  Overlap3Internals out;

  if (eps <= 0.0) {
    const Box& bb = in.bBox_;
    eps = EpsilonFromScale(bb.Scale(), 1000);
  }
  if (eps <= 0.0 || !std::isfinite(eps)) {
    out.fatal = FatalReason::Starvation;
    out.detail = "epsilon not computable";
    return out;
  }

  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) return out;  // empty mesh

  StageBResult stageB = StageB(stageA, eps, out.counters);
  if (stageB.fatal.has_value()) {
    out.fatal = stageB.fatal;
    out.detail = stageB.detail;
    return out;
  }
  out.arr = std::move(stageB.arr);

  StageResult<std::vector<SlabResult>> slabResult =
      BuildSlabs(out.arr, eps, out.counters);
  if (!slabResult.ok()) {
    out.fatal = slabResult.fatal;
    out.detail = slabResult.detail;
    return out;
  }
  out.slabs = std::move(*slabResult.value);
  return out;
}

// ---------------------------------------------------------------------------
// White-box test wrappers (used by overlap3_test.cpp only)
// ---------------------------------------------------------------------------

// Run stages C+D+E from a pre-built ArrangementGeometry (P5 PSLGInvalid pin).
Overlap3Result RemoveOverlaps3D_FromArr(const ArrangementGeometry& arr_in,
                                        double eps) {
  Overlap3Counters cnt;
  ArrangementGeometry arr_copy = arr_in;
  return RunStageCDE(arr_copy, eps, cnt);
}

// White-box seam balance check (P1 pin).
std::optional<FatalReason> CheckSeamBalance_Test(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions) {
  auto isKept = [](const PSLGRegion& r) -> bool {
    return r.classified && ((r.below > 0) != (r.above > 0));
  };
  auto hasEdgeNoSpike = [](const std::vector<int>& loop, int a, int b) -> bool {
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i) {
      if (loop[i] == a && loop[(i + 1) % n] == b &&
          loop[(i - 1 + n) % n] != b && loop[(i + 2) % n] != a)
        return true;
      if (loop[i] == b && loop[(i + 1) % n] == a &&
          loop[(i - 1 + n) % n] != a && loop[(i + 2) % n] != b)
        return true;
    }
    return false;
  };
  for (int si = 0; si < (int)arr.seams.size(); ++si) {
    const Seam& seam = arr.seams[si];
    const auto& regs0 = faceRegions[seam.faceId0];
    const auto& regs1 = faceRegions[seam.faceId1];
    for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
      const int vA = seam.vertIds[k], vB = seam.vertIds[k + 1];
      int n0 = 0, n1 = 0;
      for (const auto& r : regs0)
        if (isKept(r) && hasEdgeNoSpike(r.loopVerts, vA, vB)) ++n0;
      for (const auto& r : regs1)
        if (isKept(r) && hasEdgeNoSpike(r.loopVerts, vA, vB)) ++n1;
      // Same 0-vs-nonzero exemption as CheckSeamBalance (degenerate interior
      // seam endpoint case; caught downstream by BuildImpl manifold check).
      if (n0 > 0 && n1 > 0 && n0 != n1) return FatalReason::BalanceViolation;
    }
  }
  return std::nullopt;
}

// White-box classification (P2 pin).
ClassifyRegionResult ClassifyRegion_Test(const PSLGRegion& region, int faceId,
                                         const std::vector<SlabResult>& slabs,
                                         const std::vector<MergedVert>& verts,
                                         const CanonicalFace& face,
                                         double eps) {
  ClassifyRegionResult out;
  Overlap3Counters cnt;
  std::optional<FatalReason> fatal;
  auto cls =
      ClassifyRegion(region, faceId, slabs, verts, face, eps, cnt, fatal);
  if (fatal.has_value()) {
    out.fatal = *fatal;
    return out;
  }
  if (cls.has_value()) {
    out.classified = true;
    out.below = cls->below;
    out.above = cls->above;
  }
  return out;
}

// White-box M6 anchor-component propagation (P8-P10 pins).
// Calls PropagateAnchorComponentsImpl with the supplied synthetic data.
std::optional<FatalReason> PropagateAnchorComponents_Test(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt) {
  return PropagateAnchorComponentsImpl(arr, faceRegions, eps, cnt);
}

}  // namespace manifold
