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

// Stages A, B', D', E' of the 3D sweep-native emission prototype.
// Stage C' is in overlap3_sweep.cpp.
// Design: docs/SweepEmit3D.md (three crucible rounds).

#include "overlap3.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
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

// Stage C' forward declaration.
StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt);

namespace {

// ---------------------------------------------------------------------------
// Stage B' helpers: coplanar clip, edge-in-plane detection, seam computation.
// ---------------------------------------------------------------------------

// Sutherland-Hodgman clip of a convex polygon against one half-plane.
static std::vector<vec2> ClipPolyByHalfplane(std::vector<vec2> poly, vec2 eA,
                                             vec2 eB) {
  std::vector<vec2> out;
  const int n = (int)poly.size();
  const vec2 e = eB - eA;
  for (int k = 0; k < n; ++k) {
    const vec2 &cur = poly[k], &nxt = poly[(k + 1) % n];
    const bool ci = e.x * (cur.y - eA.y) - e.y * (cur.x - eA.x) >= 0.0;
    const bool ni = e.x * (nxt.y - eA.y) - e.y * (nxt.x - eA.x) >= 0.0;
    if (ci) out.push_back(cur);
    if (ci != ni) {
      const vec2 d1 = nxt - cur, d2 = e;
      const double den = d1.x * d2.y - d1.y * d2.x;
      if (den != 0.0) {
        const vec2 dc = eA - cur;
        out.push_back(cur + (dc.x * d2.y - dc.y * d2.x) / den * d1);
      }
    }
  }
  return out;
}

// Length of the intersection of segment [sA,sB] with the interior of triangle
// (t0,t1,t2) in 2D. Used for EdgeInPlane detection.
static double SegTriInteriorLen2D(vec2 sA, vec2 sB, vec2 t0, vec2 t1, vec2 t2) {
  double tLo = 0.0, tHi = 1.0;
  const vec2 triV[3] = {t0, t1, t2}, segD = sB - sA;
  for (int i = 0; i < 3; ++i) {
    const vec2 eA = triV[i], eB = triV[(i + 1) % 3], e = eB - eA;
    const double f0 = e.x * (sA.y - eA.y) - e.y * (sA.x - eA.x);
    const double f1 = e.x * (sB.y - eA.y) - e.y * (sB.x - eA.x);
    const double fd = f1 - f0;
    if (fd == 0.0) {
      if (f0 < 0.0) return 0.0;
    } else {
      const double tc = -f0 / fd;
      if (fd < 0.0)
        tHi = std::min(tHi, tc);
      else
        tLo = std::max(tLo, tc);
    }
    if (tLo > tHi) return 0.0;
  }
  return (tHi - tLo) * la::length(segD);
}

// Clip line P + t*D against face triangle; return true when the overlap
// interval is non-empty.  tLo/tHi receive the parameter range.
static bool LineTriClip(vec3 P, vec3 D, vec3 v0, vec3 v1, vec3 v2, vec3 n,
                        double& tLo, double& tHi, double eps) {
  tLo = -std::numeric_limits<double>::infinity();
  tHi = std::numeric_limits<double>::infinity();
  const vec3 edges[3] = {v1 - v0, v2 - v1, v0 - v2};
  const vec3 vBase[3] = {v0, v1, v2};
  const double kTol = eps * la::length(n);
  for (int i = 0; i < 3; ++i) {
    const vec3 ev = la::cross(edges[i], D),
               evB = la::cross(edges[i], vBase[i] - P);
    const double den = la::dot(n, ev), num = la::dot(n, evB);
    if (std::abs(den) <= kTol) {
      if (num > kTol) return false;
    } else {
      const double t = num / den;
      if (den > 0)
        tLo = std::max(tLo, t);
      else
        tHi = std::min(tHi, t);
    }
  }
  return tLo <= tHi + eps;
}

// Compute the seam segment [qA, qB] as the intersection of face triangles
// (a0,a1,a2) with normal na and (b0,b1,b2) with normal nb.
// Returns false if no interior seam exists (non-intersecting or degenerate).
static bool TriTriSeam(vec3 a0, vec3 a1, vec3 a2, vec3 na, vec3 b0, vec3 b1,
                       vec3 b2, vec3 nb, vec3& qA, vec3& qB, double eps) {
  const vec3 D = la::cross(na, nb);
  const double Dlen = la::length(D);
  if (Dlen == 0.0) return false;
  const vec3 dir = D / Dlen;
  const double da = la::dot(na, a0), db = la::dot(nb, b0);
  const vec3 absD = la::abs(dir);
  const int mc = (absD.x >= absD.y && absD.x >= absD.z) ? 0
                 : (absD.y >= absD.z)                   ? 1
                                                        : 2;
  const int c1 = (mc + 1) % 3, c2 = (mc + 2) % 3;
  const double a11 = na[c1], a12 = na[c2], b11 = nb[c1], b12 = nb[c2];
  const double det = a11 * b12 - a12 * b11;
  if (det == 0.0) return false;
  vec3 P(0.0);
  P[c1] = (da * b12 - db * a12) / det;
  P[c2] = (a11 * db - b11 * da) / det;
  double tAlo, tAhi, tBlo, tBhi;
  if (!LineTriClip(P, dir, a0, a1, a2, na, tAlo, tAhi, eps)) return false;
  if (!LineTriClip(P, dir, b0, b1, b2, nb, tBlo, tBhi, eps)) return false;
  const double tLo = std::max(tAlo, tBlo), tHi = std::min(tAhi, tBhi);
  if (tLo >= tHi) return false;
  qA = P + tLo * dir;
  qB = P + tHi * dir;
  return true;
}

// Compute the 3D crossing point of two line segments [A0,A1] and [B0,B1]
// assumed to be coplanar (both lie in the shared face plane).
// Returns true and sets t (param on A) and x (the crossing x-coordinate) when:
//   - the lines are non-parallel (|cross(dA,dB)| > 0)
//   - the crossing is strictly interior to both seams (t, s in (eps_t,
//   1-eps_t))
// The x is the only output that matters per spec B' (over-inclusion harmless).
static bool SeamSeamCrossX(vec3 A0, vec3 A1, vec3 B0, vec3 B1, double eps,
                           double& xOut) {
  const vec3 dA = A1 - A0, dB = B1 - B0, dC = B0 - A0;
  const vec3 cAB = la::cross(dA, dB);
  const double cABlen2 = la::dot(cAB, cAB);
  if (cABlen2 < eps * eps * eps * eps) return false;  // parallel or degenerate
  // t on A: t * |cAB|^2 = dot(cross(dC, dB), cAB)
  const double t = la::dot(la::cross(dC, dB), cAB) / cABlen2;
  // s on B: s * |cAB|^2 = dot(cross(dC, dA), cAB)
  const double s = la::dot(la::cross(dC, dA), cAB) / cABlen2;
  // Require strictly interior to both seams (not at endpoints).
  const double lenA = la::length(dA), lenB = la::length(dB);
  if (lenA < eps || lenB < eps) return false;
  const double tEps = eps / lenA, sEps = eps / lenB;
  if (t <= tEps || t >= 1.0 - tEps) return false;
  if (s <= sEps || s >= 1.0 - sEps) return false;
  xOut = A0.x + t * dA.x;
  return true;
}

// Find or add a vert within eps of pos; return its index.
static int FindOrAddVert(std::vector<MergedVert>& verts, vec3 pos, double eps) {
  for (int i = 0; i < (int)verts.size(); ++i)
    if (la::length(verts[i].pos - pos) <= eps) return i;
  const int id = (int)verts.size();
  verts.push_back({pos});
  return id;
}

// ---------------------------------------------------------------------------
// Stage A: vert merge, multiplicity accumulation.
// ---------------------------------------------------------------------------

struct StageAResult {
  std::vector<vec3> mergedVerts;
  std::vector<int> vertMap;
  std::vector<CanonicalFace> faces;
};

static StageAResult StageA(const Manifold::Impl& in, double eps) {
  const int nVerts = (int)in.vertPos_.size();
  const int nTris = (int)in.halfedge_.size() / 3;

  DisjointSets uf(nVerts);
  for (int i = 0; i < nVerts; ++i)
    for (int j = i + 1; j < nVerts; ++j)
      if (la::length(in.vertPos_[i] - in.vertPos_[j]) <= eps) uf.unite(i, j);

  std::map<int, std::vector<int>> comps;
  for (int i = 0; i < nVerts; ++i) comps[(int)uf.find(i)].push_back(i);

  std::vector<int> vertMap(nVerts);
  std::vector<vec3> mergedVerts;
  std::map<int, int> rootToMerged;
  for (auto& [root, members] : comps) {
    vec3 c(0.0);
    for (int v : members) c += in.vertPos_[v];
    c /= (double)members.size();
    int best = members[0];
    double bestD = la::length(in.vertPos_[best] - c);
    for (int v : members) {
      double d = la::length(in.vertPos_[v] - c);
      if (d < bestD) {
        bestD = d;
        best = v;
      }
    }
    rootToMerged[root] = (int)mergedVerts.size();
    mergedVerts.push_back(in.vertPos_[best]);
  }
  for (int i = 0; i < nVerts; ++i) vertMap[i] = rootToMerged[(int)uf.find(i)];

  struct FaceKey {
    ivec3 s;
    bool operator<(const FaceKey& o) const {
      if (s.x != o.s.x) return s.x < o.s.x;
      if (s.y != o.s.y) return s.y < o.s.y;
      return s.z < o.s.z;
    }
  };
  struct FaceEntry {
    int64_t mult;
    int repTri;
    int repPar;
  };
  std::map<FaceKey, FaceEntry> faceMap;

  for (int tri = 0; tri < nTris; ++tri) {
    const int h = 3 * tri;
    int v0 = vertMap[in.halfedge_.Start(h)],
        v1 = vertMap[in.halfedge_.Start(h + 1)],
        v2 = vertMap[in.halfedge_.Start(h + 2)];
    if (v0 == v1 || v1 == v2 || v0 == v2) continue;
    int a = v0, b = v1, c = v2, par = 1;
    if (a > b) {
      std::swap(a, b);
      par = -par;
    }
    if (b > c) {
      std::swap(b, c);
      par = -par;
    }
    if (a > b) {
      std::swap(a, b);
      par = -par;
    }
    FaceKey key{{a, b, c}};
    auto it = faceMap.find(key);
    if (it == faceMap.end())
      faceMap.emplace(key, FaceEntry{1, tri, par});
    else
      it->second.mult += (par == it->second.repPar) ? 1 : -1;
  }

  StageAResult res;
  res.mergedVerts = std::move(mergedVerts);
  res.vertMap = std::move(vertMap);
  for (auto& [key, e] : faceMap) {
    if (e.mult == 0) continue;
    const int h = 3 * e.repTri;
    const int mv0 = res.vertMap[in.halfedge_.Start(h)];
    const int mv1 = res.vertMap[in.halfedge_.Start(h + 1)];
    const int mv2 = res.vertMap[in.halfedge_.Start(h + 2)];
    const vec3 n = la::cross(res.mergedVerts[mv1] - res.mergedVerts[mv0],
                             res.mergedVerts[mv2] - res.mergedVerts[mv0]);
    if (la::dot(n, n) == 0.0) continue;
    CanonicalFace cf;
    cf.id = e.repTri;
    cf.verts = {mv0, mv1, mv2};
    cf.normal = la::normalize(n);
    cf.mult = e.mult;
    res.faces.push_back(cf);
  }
  return res;
}

// ---------------------------------------------------------------------------
// Stage B': CoplanarOverlap + EdgeInPlane detection, seam computation.
// No triple-point unification, no edge subdivision, no faceSeams.
// ---------------------------------------------------------------------------

struct StageBResult {
  ArrangementGeometry arr;
  std::optional<FatalReason> fatal;
  std::string detail;
};

static StageBResult StageBPrime(const StageAResult& stageA, double eps,
                                Overlap3Counters& cnt) {
  StageBResult res;
  ArrangementGeometry& arr = res.arr;
  const int nFaces = (int)stageA.faces.size();

  arr.verts.resize(stageA.mergedVerts.size());
  for (int i = 0; i < (int)stageA.mergedVerts.size(); ++i)
    arr.verts[i] = {stageA.mergedVerts[i]};
  arr.faces = stageA.faces;

  struct Box3 {
    vec3 mn, mx;
  };
  std::vector<Box3> boxes(nFaces);
  std::vector<double> planeDist(nFaces);
  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& f = stageA.faces[fi];
    planeDist[fi] = la::dot(f.normal, arr.verts[f.verts.x].pos);
    const vec3 p0 = arr.verts[f.verts.x].pos, p1 = arr.verts[f.verts.y].pos,
               p2 = arr.verts[f.verts.z].pos;
    boxes[fi] = {la::min(la::min(p0, p1), p2) - vec3(eps),
                 la::max(la::max(p0, p1), p2) + vec3(eps)};
  }

  for (int fi = 0; fi < nFaces; ++fi) {
    for (int fj = fi + 1; fj < nFaces; ++fj) {
      const Box3 &bi = boxes[fi], &bj = boxes[fj];
      if (bi.mn.x > bj.mx.x || bi.mx.x < bj.mn.x || bi.mn.y > bj.mx.y ||
          bi.mx.y < bj.mn.y || bi.mn.z > bj.mx.z || bi.mx.z < bj.mn.z)
        continue;

      const CanonicalFace &FA = stageA.faces[fi], &FB = stageA.faces[fj];
      {
        int shared = 0;
        const int va[3] = {FA.verts.x, FA.verts.y, FA.verts.z};
        const int vb[3] = {FB.verts.x, FB.verts.y, FB.verts.z};
        for (int a : va)
          for (int b : vb)
            if (a == b) ++shared;
        if (shared >= 2) continue;
      }

      const vec3 &na = FA.normal, &nb = FB.normal;
      const vec3 pa0 = arr.verts[FA.verts.x].pos,
                 pa1 = arr.verts[FA.verts.y].pos,
                 pa2 = arr.verts[FA.verts.z].pos;
      const vec3 pb0 = arr.verts[FB.verts.x].pos,
                 pb1 = arr.verts[FB.verts.y].pos,
                 pb2 = arr.verts[FB.verts.z].pos;

      // Coplanar detection: called for both directions (B in A's plane,
      // A in B's plane). Returns true if the 2D clipped area > eps^2.
      auto coplanarInteriorOverlap = [&](vec3 n, const vec3 pA[3],
                                         const vec3 pB[3]) -> bool {
        const mat2x3 proj = GetAxisAlignedProjection(n);
        const vec2 a2[3] = {proj * pA[0], proj * pA[1], proj * pA[2]};
        const vec2 b2[3] = {proj * pB[0], proj * pB[1], proj * pB[2]};
        std::vector<vec2> poly = {b2[0], b2[1], b2[2]};
        for (int k = 0; k < 3 && !poly.empty(); ++k)
          poly = ClipPolyByHalfplane(poly, a2[k], a2[(k + 1) % 3]);
        double area = 0.0;
        for (int k = 0; k < (int)poly.size(); ++k) {
          const vec2 &p = poly[k], &q = poly[(k + 1) % (int)poly.size()];
          area += p.x * q.y - q.x * p.y;
        }
        return std::abs(area) * 0.5 > eps * eps;
      };
      const vec3 pA[3] = {pa0, pa1, pa2}, pB[3] = {pb0, pb1, pb2};

      // B in A's plane?
      if (std::abs(la::dot(na, pb0) - planeDist[fi]) <= eps &&
          std::abs(la::dot(na, pb1) - planeDist[fi]) <= eps &&
          std::abs(la::dot(na, pb2) - planeDist[fi]) <= eps) {
        if (coplanarInteriorOverlap(na, pA, pB)) {
          res.fatal = FatalReason::CoplanarOverlap;
          res.detail = "coplanar face pair with interior overlap";
          return res;
        }
        continue;
      }
      // A in B's plane?
      const double dB = la::dot(nb, pb0);
      if (std::abs(la::dot(nb, pa0) - dB) <= eps &&
          std::abs(la::dot(nb, pa1) - dB) <= eps &&
          std::abs(la::dot(nb, pa2) - dB) <= eps) {
        if (coplanarInteriorOverlap(nb, pB, pA)) {
          res.fatal = FatalReason::CoplanarOverlap;
          res.detail = "coplanar face pair with interior overlap";
          return res;
        }
        continue;
      }

      const vec3 Dcross = la::cross(na, nb);
      if (la::length(Dcross) == 0.0) continue;

      // EdgeInPlane detection.
      {
        const double dA = planeDist[fi], dB = la::dot(nb, pb0);
        const vec3 eA[3] = {pa0, pa1, pa2}, eB[3] = {pb0, pb1, pb2};
        const mat2x3 projB = GetAxisAlignedProjection(nb),
                     projA = GetAxisAlignedProjection(na);
        const vec2 pb02 = projB * pb0, pb12 = projB * pb1, pb22 = projB * pb2;
        const vec2 pa02 = projA * pa0, pa12 = projA * pa1, pa22 = projA * pa2;
        for (int ei = 0; ei < 3; ++ei) {
          const int ej = (ei + 1) % 3;
          if (std::abs(la::dot(nb, eA[ei]) - dB) <= eps &&
              std::abs(la::dot(nb, eA[ej]) - dB) <= eps &&
              SegTriInteriorLen2D(projB * eA[ei], projB * eA[ej], pb02, pb12,
                                  pb22) > eps) {
            res.fatal = FatalReason::EdgeInPlane;
            res.detail = "edge of face " + std::to_string(fi) +
                         " in plane of face " + std::to_string(fj);
            return res;
          }
          if (std::abs(la::dot(na, eB[ei]) - dA) <= eps &&
              std::abs(la::dot(na, eB[ej]) - dA) <= eps &&
              SegTriInteriorLen2D(projA * eB[ei], projA * eB[ej], pa02, pa12,
                                  pa22) > eps) {
            res.fatal = FatalReason::EdgeInPlane;
            res.detail = "edge of face " + std::to_string(fj) +
                         " in plane of face " + std::to_string(fi);
            return res;
          }
        }
      }

      // Seam computation.
      vec3 qA, qB;
      if (!TriTriSeam(pa0, pa1, pa2, na, pb0, pb1, pb2, nb, qA, qB, eps))
        continue;
      const double seamLen = la::length(qB - qA);

      if (seamLen <= eps) {
        // Degenerate contact: add endpoints to arr.verts for use as
        // x-criticals. SubEpsFeature guard in stage C' will catch any
        // macro-scale hazard.
        FindOrAddVert(arr.verts, qA, eps);
        FindOrAddVert(arr.verts, qB, eps);
        ++cnt.subEpsContactsDropped;
        continue;
      }

      const int vA = FindOrAddVert(arr.verts, qA, eps);
      const int vB = FindOrAddVert(arr.verts, qB, eps);
      if (vA == vB) {
        // Collapsed after snap: treat as degenerate.
        ++cnt.subEpsContactsDropped;
        continue;
      }

      arr.seams.push_back({fi, fj, vA, vB});
    }
  }

  // M1: seam-seam triple criticals (spec B'). For each pair of seams sharing
  // a face, compute the crossing of their in-plane projections. The crossing x
  // is the only quantity consumed (over-inclusion is harmless).
  const int nSeams = (int)arr.seams.size();
  for (int si = 0; si < nSeams; ++si) {
    for (int sj = si + 1; sj < nSeams; ++sj) {
      const Seam& SA = arr.seams[si];
      const Seam& SB = arr.seams[sj];
      // Check for shared face.
      if (SA.faceId0 != SB.faceId0 && SA.faceId0 != SB.faceId1 &&
          SA.faceId1 != SB.faceId0 && SA.faceId1 != SB.faceId1)
        continue;
      const vec3 A0 = arr.verts[SA.vertId0].pos;
      const vec3 A1 = arr.verts[SA.vertId1].pos;
      const vec3 B0 = arr.verts[SB.vertId0].pos;
      const vec3 B1 = arr.verts[SB.vertId1].pos;
      double xCross;
      if (!SeamSeamCrossX(A0, A1, B0, B1, eps, xCross)) continue;
      // Add a vertex carrying the crossing x into the critical set.
      // Compute the 3D crossing from the parameter on A (recovered from x).
      const double dAx = A1.x - A0.x;
      const vec3 crossPt =
          (std::abs(dAx) > eps)
              ? A0 + ((xCross - A0.x) / dAx) * (A1 - A0)
              : vec3{xCross, (A0.y + A1.y) * 0.5, (A0.z + A1.z) * 0.5};
      FindOrAddVert(arr.verts, crossPt, eps);
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// Stage D': strip emission.
// For each built slab, for each retained piece, emit a quad (2 triangles)
// by track-extending the (y,z) endpoints to xLo and xHi.
// ---------------------------------------------------------------------------

// Evaluate a piece endpoint in (y,z) at x=xTarget using the face's FaceTrack.
// t = parameterization of pt on segment [ft.p0, ft.p1].
static vec2 ExtendPt(vec2 pt, const FaceTrack& ft, double xTarget) {
  const vec2 segD = ft.p1 - ft.p0;
  const double segLen2 = la::dot(segD, segD);
  double t = 0.5;
  if (segLen2 > 0.0) {
    t = la::dot(pt - ft.p0, segD) / segLen2;
    t = std::max(0.0, std::min(1.0, t));
  }
  const vec2 yz0 = InterpolateSafe(ft.va0, ft.vb0, xTarget);
  const vec2 yz1 = InterpolateSafe(ft.va1, ft.vb1, xTarget);
  return (1.0 - t) * yz0 + t * yz1;
}

// Like ExtendPt but uses the seam track for class-ii endpoints (spec D'/E').
// A class-ii endpoint is interior to the face's section segment (not at t~0
// or t~1) and lies on a seam.  The seam track gives the correct extension;
// the face-edge track gives the wrong yz because it follows a diagonal edge
// instead of the seam line.
static vec2 ExtendPtWithSeams(vec2 pt, const FaceTrack& ft,
                              const std::vector<SeamTrackEntry>& seamTracks,
                              double xTarget, double eps) {
  // Check if pt is strictly interior to the face segment (class-ii indicator).
  const vec2 segD = ft.p1 - ft.p0;
  const double segLen2 = la::dot(segD, segD);
  if (segLen2 > eps * eps) {
    const double atP0 = la::length(pt - ft.p0);
    const double atP1 = la::length(pt - ft.p1);
    if (atP0 > eps && atP1 > eps) {
      // Interior point: search for matching seam track (class-ii).
      for (const auto& st : seamTracks) {
        if (la::length(pt - st.yzMid) <= eps) {
          return InterpolateSafe(st.vA, st.vB, xTarget);
        }
      }
      // No seam track matches: forced-through weld endpoint.
      // Spec D' (R2-fold): extend CONSTANT in (y,z) across the slab.
      return pt;
    }
  }
  // Class-i endpoint: extend via the face edge track.
  return ExtendPt(pt, ft, xTarget);
}

struct OutTri3D {
  vec3 v[3];
};

// Find interior split t-values in (0,1) where cap arrangement vertices lie on
// the strip segment [p0c, p1c] (within eps).  Returned values are sorted.
static std::vector<double> FindCapSplits(const std::vector<vec2>& capVerts,
                                         vec2 p0c, vec2 p1c, double eps) {
  std::vector<double> result;
  const vec2 segD = p1c - p0c;
  const double segLen2 = la::dot(segD, segD);
  if (segLen2 <= eps * eps) return result;
  const double invLen2 = 1.0 / segLen2;
  const double tEps = eps / std::sqrt(segLen2);
  for (const auto& v : capVerts) {
    const double t = la::dot(v - p0c, segD) * invLen2;
    if (t <= tEps || t >= 1.0 - tEps) continue;
    const vec2 proj = p0c + t * segD;
    if (la::length(v - proj) > eps) continue;
    result.push_back(t);
  }
  std::sort(result.begin(), result.end());
  auto last = std::unique(result.begin(), result.end(),
                          [](double a, double b) { return b - a < 1e-10; });
  result.erase(last, result.end());
  return result;
}

// Zipper-triangulate a strip between two polylines parameterized from t=0 to 1.
// s_params: sorted xLo polyline split params (must include 0 and 1).
// t_params: sorted xHi polyline split params (must include 0 and 1).
// Each side is subdivided independently by the cap arrangement at that
// critical x (spec R1-fold: no spurious interior vertices on the other side).
static void ZipperEmit(double xLo, double xHi, vec2 yz_a0, vec2 yz_a1,
                       vec2 yz_b0, vec2 yz_b1,
                       const std::vector<double>& s_params,
                       const std::vector<double>& t_params,
                       std::vector<OutTri3D>& out) {
  auto ptA = [&](double s) -> vec3 {
    const vec2 yz = (1.0 - s) * yz_a0 + s * yz_a1;
    return {xLo, yz.x, yz.y};
  };
  auto ptB = [&](double t) -> vec3 {
    const vec2 yz = (1.0 - t) * yz_b0 + t * yz_b1;
    return {xHi, yz.x, yz.y};
  };
  auto emitTri = [&](vec3 v0, vec3 v1, vec3 v2) {
    if (v0 != v1 && v1 != v2 && v0 != v2) out.push_back({v0, v1, v2});
  };

  const int m = (int)s_params.size() - 1;
  const int n = (int)t_params.size() - 1;
  int i = 0, j = 0;
  while (i < m || j < n) {
    const double s_next = (i < m) ? s_params[i + 1] : 2.0;
    const double t_next = (j < n) ? t_params[j + 1] : 2.0;
    const vec3 Ai = ptA(s_params[i]);
    const vec3 Bj = ptB(t_params[j]);
    if (t_next < s_next - 1e-10) {
      emitTri(Ai, ptB(t_next), Bj);  // advance right
      ++j;
    } else if (s_next < t_next - 1e-10) {
      emitTri(Ai, ptA(s_next), Bj);  // advance left
      ++i;
    } else {
      emitTri(Ai, ptA(s_next), ptB(t_next));  // advance both (quad tri 1)
      emitTri(Ai, ptB(t_next), Bj);           // advance both (quad tri 2)
      ++i;
      ++j;
    }
  }
}

// Stage D': emit strips with cap-vertex subdivision at each c-side.
// capVertsMap[xCrit] = all yz vertices from the cap arrangement at xCrit.
// Per spec R1-fold: each c-side takes its polyline from the cap arrangement
// at that critical x; the two sides are subdivided independently.
static void EmitStrips(const std::vector<SlabResult>& slabs,
                       const std::map<double, std::vector<vec2>>& capVertsMap,
                       std::vector<OutTri3D>& out, double eps) {
  const std::vector<vec2> emptyVerts;
  // eps-based map lookup: seam endpoints can land at a slightly different
  // x than the original vertex, making the slab boundary differ from the
  // capVertsMap key by a sub-eps amount.
  auto capVertsLookup = [&](double x) -> const std::vector<vec2>& {
    auto it = capVertsMap.lower_bound(x - eps);
    if (it != capVertsMap.end() && std::abs(it->first - x) <= eps)
      return it->second;
    return emptyVerts;
  };
  for (const auto& slab : slabs) {
    if (!slab.built) continue;
    const double xLo = slab.xLo, xHi = slab.xHi;
    const std::vector<vec2>& capLo = capVertsLookup(xLo);
    const std::vector<vec2>& capHi = capVertsLookup(xHi);

    for (const auto& piece : slab.pieces) {
      const FaceTrack* ft = nullptr;
      for (const auto& tk : slab.faceTracks) {
        if (tk.faceId == piece.sourceId) {
          ft = &tk;
          break;
        }
      }
      if (!ft) continue;

      // Piece is emission-oriented (interior-on-left) per spec engine contract.
      const vec2 p0 = piece.from;
      const vec2 p1 = piece.to;

      const vec2 yz_a0 = ExtendPtWithSeams(p0, *ft, slab.seamTracks, xLo, eps);
      const vec2 yz_a1 = ExtendPtWithSeams(p1, *ft, slab.seamTracks, xLo, eps);
      const vec2 yz_b0 = ExtendPtWithSeams(p0, *ft, slab.seamTracks, xHi, eps);
      const vec2 yz_b1 = ExtendPtWithSeams(p1, *ft, slab.seamTracks, xHi, eps);

      // Independent splits at each c-side from the cap at that critical x.
      std::vector<double> s_params = {0.0};
      for (double t : FindCapSplits(capLo, yz_a0, yz_a1, eps))
        s_params.push_back(t);
      s_params.push_back(1.0);

      std::vector<double> t_params = {0.0};
      for (double t : FindCapSplits(capHi, yz_b0, yz_b1, eps))
        t_params.push_back(t);
      t_params.push_back(1.0);

      ZipperEmit(xLo, xHi, yz_a0, yz_a1, yz_b0, yz_b1, s_params, t_params, out);
    }
  }
}

// ---------------------------------------------------------------------------
// Stage E': cap emission.
// At each critical x=c, run ONE combined 2D arrangement over the extended
// pieces of both adjacent slabs (spec "ONE cap arrangement").  This avoids
// the double-arrangement failure: if L and R are assembled separately then
// combined in Boolean2D, the two independent vertex-snap passes can produce
// slightly different coordinates for shared boundary points, causing
// OutEdgesToPolygons to see open walks.
// ---------------------------------------------------------------------------

// Triangulate a Polygons and append cap triangles to out, emitting at x=xCap.
// If flipWinding is true, reverse each triangle (for -x normal caps).
// idx is global across all loops (flat concatenation of allVerts).
static void TriangulateCap(const Polygons& polys, double xCap, bool flipWinding,
                           double eps, std::vector<OutTri3D>& out) {
  if (polys.empty()) return;

  PolygonsIdx pidx;
  std::vector<vec2> allVerts;
  {
    int gIdx = 0;
    for (const auto& loop : polys) {
      SimplePolygonIdx spi;
      spi.reserve(loop.size());
      for (const auto& v : loop) {
        allVerts.push_back(v);
        spi.push_back({v, gIdx++});
      }
      pidx.push_back(std::move(spi));
    }
  }

  const auto tris = TriangulateIdx(pidx, eps);

  for (const auto& t : tris) {
    if (t.x < 0 || t.x >= (int)allVerts.size()) continue;
    if (t.y < 0 || t.y >= (int)allVerts.size()) continue;
    if (t.z < 0 || t.z >= (int)allVerts.size()) continue;
    const vec3 v0 = {xCap, allVerts[t.x].x, allVerts[t.x].y};
    const vec3 v1 = {xCap, allVerts[t.y].x, allVerts[t.y].y};
    const vec3 v2 = {xCap, allVerts[t.z].x, allVerts[t.z].y};
    if (!flipWinding)
      out.push_back({v0, v1, v2});
    else
      out.push_back({v0, v2, v1});
  }
}

// Build one raw-vertex list for both L and R slab pieces extended to xCap.
// Returns (rawVerts, lFrom, lTo, rFrom, rTo) where lFrom[i]/lTo[i] are
// indices into rawVerts for the i-th L piece, and similarly for R.
// NO vertex pre-dedup: RemoveOverlaps2D's MergeVerts handles it.
// Both cap_plus and cap_minus arrangements below receive the SAME rawVerts
// in the SAME ORDER, so MergeVerts produces the SAME snap for both.
struct CapEdgeSet {
  std::vector<vec2> rawVerts;
  std::vector<int> lFrom, lTo;
  std::vector<int> rFrom, rTo;
};

static CapEdgeSet BuildCapEdgeSet(const SlabResult* leftSlab,
                                  const SlabResult* rightSlab, double xCap,
                                  double eps) {
  CapEdgeSet ces;
  auto addPieces = [&](const SlabResult& slab, std::vector<int>& froms,
                       std::vector<int>& tos) {
    for (const auto& piece : slab.pieces) {
      const FaceTrack* ft = nullptr;
      for (const auto& tk : slab.faceTracks)
        if (tk.faceId == piece.sourceId) {
          ft = &tk;
          break;
        }
      if (!ft) continue;
      // Piece is emission-oriented (interior-on-left) per spec engine contract.
      const vec2 eFrom =
          ExtendPtWithSeams(piece.from, *ft, slab.seamTracks, xCap, eps);
      const vec2 eTo =
          ExtendPtWithSeams(piece.to, *ft, slab.seamTracks, xCap, eps);
      froms.push_back((int)ces.rawVerts.size());
      ces.rawVerts.push_back(eFrom);
      tos.push_back((int)ces.rawVerts.size());
      ces.rawVerts.push_back(eTo);
    }
  };
  // L verts appended first, R verts second - consistent order for both
  // cap_plus and cap_minus arrangements below.
  if (leftSlab && leftSlab->built) addPieces(*leftSlab, ces.lFrom, ces.lTo);
  if (rightSlab && rightSlab->built) addPieces(*rightSlab, ces.rFrom, ces.rTo);
  return ces;
}

// Run TWO arrangements over the cap at xCap using the SAME raw-vert set built
// by BuildCapEdgeSet. Passing the same rawVerts to both ensures MergeVerts
// yields the same snap, so cap_plus and cap_minus share a consistent vertex
// geometry (spec M4 "one source of truth per critical").
// capVerts accumulates all yz vertices from both arrangements for strip
// subdivision at this critical (spec E' R1-fold).
static void ComputeCap(const SlabResult* leftSlab, const SlabResult* rightSlab,
                       double xCap, double eps, std::vector<OutTri3D>& out,
                       std::vector<vec2>& capVerts) {
  const CapEdgeSet ces = BuildCapEdgeSet(leftSlab, rightSlab, xCap, eps);
  if (ces.lFrom.empty() && ces.rFrom.empty()) return;

  auto collectVerts = [&](const std::vector<vec2>& rv) {
    for (const auto& v : rv) {
      bool found = false;
      for (const auto& cv : capVerts)
        if (la::length(v - cv) <= eps) {
          found = true;
          break;
        }
      if (!found) capVerts.push_back(v);
    }
  };

  auto buildEdges = [&](int lSign, int rSign) -> std::vector<EdgeM> {
    std::vector<EdgeM> edges;
    for (int i = 0; i < (int)ces.lFrom.size(); ++i)
      edges.push_back({ces.lFrom[i], ces.lTo[i], lSign});
    for (int i = 0; i < (int)ces.rFrom.size(); ++i)
      edges.push_back({ces.rFrom[i], ces.rTo[i], rSign});
    return edges;
  };

  // cap_plus = L - R: L(+1), R(-1).  cap_minus = R - L: R(+1), L(-1).
  // Both use ces.rawVerts so MergeVerts snaps identically.
  {
    const std::vector<EdgeM> cpEdges = buildEdges(+1, -1);
    if (!cpEdges.empty()) {
      const OverlapResult r = RemoveOverlaps2D(ces.rawVerts, cpEdges, eps);
      collectVerts(r.verts);
      if (!r.edges.empty()) {
        const Polygons cp = OutEdgesToPolygons(r.verts, r.edges);
        TriangulateCap(cp, xCap, /*flipWinding=*/false, eps, out);
      }
    }
  }
  {
    const std::vector<EdgeM> cmEdges = buildEdges(-1, +1);
    if (!cmEdges.empty()) {
      const OverlapResult r = RemoveOverlaps2D(ces.rawVerts, cmEdges, eps);
      collectVerts(r.verts);
      if (!r.edges.empty()) {
        const Polygons cm = OutEdgesToPolygons(r.verts, r.edges);
        TriangulateCap(cm, xCap, /*flipWinding=*/true, eps, out);
      }
    }
  }
}

// Returns a map from each critical x to the yz vertex set of its cap
// arrangement, for use by EmitStrips (spec R1-fold).
static std::map<double, std::vector<vec2>> EmitCaps(
    const std::vector<SlabResult>& slabs, const std::vector<double>& crits,
    double eps, std::vector<OutTri3D>& out) {
  std::map<double, std::vector<vec2>> capVertsMap;
  const int nSlabs = (int)slabs.size();
  for (double c : crits) {
    int leftIdx = -1, rightIdx = -1;
    for (int si = 0; si < nSlabs; ++si) {
      if (!slabs[si].built) continue;
      if (slabs[si].xHi <= c + eps * 0.1) leftIdx = si;
    }
    for (int si = 0; si < nSlabs; ++si) {
      if (!slabs[si].built) continue;
      if (slabs[si].xLo >= c - eps * 0.1) {
        rightIdx = si;
        break;
      }
    }
    const SlabResult* leftSlab = (leftIdx >= 0) ? &slabs[leftIdx] : nullptr;
    const SlabResult* rightSlab = (rightIdx >= 0) ? &slabs[rightIdx] : nullptr;
    ComputeCap(leftSlab, rightSlab, c, eps, out, capVertsMap[c]);
  }
  return capVertsMap;
}

// ---------------------------------------------------------------------------
// Assembly: collect emitted triangles, dedup verts, build Impl.
// ---------------------------------------------------------------------------

static StageResult<Manifold::Impl> BuildImpl(const std::vector<OutTri3D>& tris,
                                             double eps) {
  if (tris.empty()) return StageResult<Manifold::Impl>::Ok(Manifold::Impl{});

  // Collect verts with eps-dedup.
  std::vector<vec3> verts;
  auto getVertIdx = [&](vec3 p) -> int {
    for (int i = 0; i < (int)verts.size(); ++i)
      if (la::length(verts[i] - p) <= eps) return i;
    int id = (int)verts.size();
    verts.push_back(p);
    return id;
  };

  // Filter degenerate and duplicate triangles.
  // Degenerate: strip quads whose corners collapse within eps produce v0==v1
  // etc., which would crash CreateHalfedges.
  // Duplicate: sub-eps critical pairs can emit caps at x=c1 and x=c2 where
  // |c2-c1|<eps; both cap triangulations collapse to the same vertex indices
  // after eps-dedup.  Emitting each triangle twice is non-manifold.  Spec E'
  // says "their differences are empty" - dropping the second copy is correct.
  Vec<ivec3> tv;
  tv.reserve(tris.size());
  std::set<std::tuple<int, int, int>> seenTris;
  for (const auto& tri : tris) {
    const int v0 = getVertIdx(tri.v[0]);
    const int v1 = getVertIdx(tri.v[1]);
    const int v2 = getVertIdx(tri.v[2]);
    if (v0 == v1 || v1 == v2 || v0 == v2) continue;
    // Canonical key: rotate so smallest vertex is first, preserving
    // orientation.
    int a = v0, b = v1, c = v2;
    if (b < a && b < c) {
      a = v1;
      b = v2;
      c = v0;
    } else if (c < a && c < b) {
      a = v2;
      b = v0;
      c = v1;
    }
    if (!seenTris.insert({a, b, c}).second) continue;  // exact duplicate
    tv.push_back({v0, v1, v2});
  }

  Manifold::Impl impl;
  impl.vertPos_.resize(verts.size());
  for (int i = 0; i < (int)verts.size(); ++i) impl.vertPos_[i] = verts[i];

  impl.CreateHalfedges(tv);
  if (!impl.IsManifold()) {
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::NonManifoldEmission,
        "emitted triangulation not 2-manifold");
  }
  impl.InitializeOriginal();
  impl.CalculateBBox();
  impl.SetEpsilon();
  impl.SortGeometry();
  impl.SetNormalsAndCoplanar();
  return StageResult<Manifold::Impl>::Ok(std::move(impl));
}

}  // namespace

// ---------------------------------------------------------------------------
// Pipeline: stages C' + D' + E'.
// ---------------------------------------------------------------------------

static Overlap3Result RunCDEPrime(ArrangementGeometry& arr, double eps,
                                  Overlap3Counters& cnt) {
  Overlap3Result result;

  auto slabRes = BuildSlabs(arr, eps, cnt);
  if (!slabRes.ok()) {
    result.fatal = slabRes.fatal;
    result.detail = slabRes.detail;
    result.counters = cnt;
    return result;
  }
  const std::vector<SlabResult>& slabs = *slabRes.value;

  // Collect x-criticals (all vert x-values) for cap generation.
  // Use eps-based dedup: seam endpoints can land slightly off a mesh vertex
  // x-coordinate, creating sub-eps duplicate criticals that would invoke
  // ComputeCap twice with the same adjacent slabs and emit double caps.
  // Spec E' (R3-fold): every critical gets its own cap, runs are NOT merged.
  // Only exact duplicates are collapsed (same IEEE bit pattern after sort).
  std::vector<double> crits;
  crits.reserve(arr.verts.size());
  for (const auto& v : arr.verts) crits.push_back(v.pos.x);
  std::sort(crits.begin(), crits.end());
  crits.erase(std::unique(crits.begin(), crits.end()), crits.end());

  // E': emit caps at each critical, collecting arrangement vertices.
  // Caps run first so their vertex set can subdivide adjacent strip edges
  // (spec R1-fold: one source of truth per critical).
  std::vector<OutTri3D> emitted;
  const auto capVertsMap = EmitCaps(slabs, crits, eps, emitted);

  // D': emit strips, subdividing each c-side by the cap arrangement vertices.
  EmitStrips(slabs, capVertsMap, emitted, eps);

  // Assembly.
  auto buildRes = BuildImpl(emitted, eps);
  if (!buildRes.ok()) {
    result.fatal = buildRes.fatal;
    result.detail = buildRes.detail;
    result.counters = cnt;
    return result;
  }
  result.impl = std::move(buildRes.value);
  result.counters = cnt;
  return result;
}

// ---------------------------------------------------------------------------
// Public entry points.
// ---------------------------------------------------------------------------

Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps) {
  Overlap3Result result;
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    result.fatal = FatalReason::SubEpsInput;
    result.detail = "epsilon not computable";
    return result;
  }
  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) {
    result.impl = Manifold::Impl{};
    return result;
  }
  Overlap3Counters& cnt = result.counters;
  StageBResult stageB = StageBPrime(stageA, eps, cnt);
  if (stageB.fatal.has_value()) {
    result.fatal = stageB.fatal;
    result.detail = stageB.detail;
    return result;
  }
  return RunCDEPrime(stageB.arr, eps, cnt);
}

Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps) {
  Overlap3Internals out;
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    out.fatal = FatalReason::SubEpsInput;
    out.detail = "epsilon not computable";
    return out;
  }
  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) return out;
  StageBResult stageB = StageBPrime(stageA, eps, out.counters);
  if (stageB.fatal.has_value()) {
    out.fatal = stageB.fatal;
    out.detail = stageB.detail;
    return out;
  }
  out.arr = std::move(stageB.arr);
  auto slabRes = BuildSlabs(out.arr, eps, out.counters);
  if (!slabRes.ok()) {
    out.fatal = slabRes.fatal;
    out.detail = slabRes.detail;
    return out;
  }
  out.slabs = std::move(*slabRes.value);
  return out;
}

}  // namespace manifold
