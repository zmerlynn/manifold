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
// Returns false only when the triangles are disjoint (or parallel-plane);
// a point/tangent contact admitted by the eps clip model returns true with
// qA == qB, so the caller's degenerate-contact arm records its x-critical.
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
  // dir is unit length, so the parameter interval is in length units and
  // eps applies directly (matching LineTriClip's own slop).
  if (tLo > tHi + eps) return false;  // genuinely disjoint
  if (tLo >= tHi) {
    // Point/tangent contact: collapse to the interval midpoint.
    qA = qB = P + (0.5 * (tLo + tHi)) * dir;
    return true;
  }
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
  const double lenA = la::length(dA), lenB = la::length(dB);
  if (lenA < eps || lenB < eps) return false;
  const vec3 cAB = la::cross(dA, dB);
  const double cABlen2 = la::dot(cAB, cAB);
  // Near-parallel gate, dimensionally correct: |cross| has units len^2, so
  // compare against eps * (lenA + lenB).  Endpoint noise on near-parallel
  // seams otherwise yields a pseudo-crossing at a meaningless x.
  const double parTol = eps * (lenA + lenB);
  if (cABlen2 <= parTol * parTol) return false;
  // t on A: t * |cAB|^2 = dot(cross(dC, dB), cAB)
  const double t = la::dot(la::cross(dC, dB), cAB) / cABlen2;
  // s on B: s * |cAB|^2 = dot(cross(dC, dA), cAB)
  const double s = la::dot(la::cross(dC, dA), cAB) / cABlen2;
  // Require strictly interior to both seams (not at endpoints).
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
        // Faces sharing a merged edge are exempt from pair analysis.  Their
        // planes meet along the shared edge, so seam computation would return
        // that edge (its verts are already criticals) and the edge-in-plane
        // detector would flag the shared edge itself.  This exemption also
        // covers the coplanar shared-edge cases, in both directions:
        // opposite-diagonal triangulations of two solids touching on a
        // common plane are LEGAL (pinned by TouchingDisjoint) and locally
        // indistinguishable from a folded flap here - the hazardous
        // same-winding case fails closed downstream instead, where its
        // coincident section edges carry two source ids with nonzero net
        // multiplicity (EngineIdConflict); the legal touching case cancels
        // (net zero) before any conflict is recorded.
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
      // A in B's plane). Returns true if the clipped region is a genuine
      // interior overlap: area > perimeter * eps, the house dimensionally-
      // correct threshold (an overlap wider than eps somewhere).  Clip
      // roundoff on adjacent coplanar triangles produces slivers of area
      // ~ scale * machine-eps, far below it; sub-eps-wide contact is legal
      // touching under the eps contract.
      auto coplanarInteriorOverlap = [&](vec3 n, const vec3 pA[3],
                                         const vec3 pB[3]) -> bool {
        const mat2x3 proj = GetAxisAlignedProjection(n);
        const vec2 a2[3] = {proj * pA[0], proj * pA[1], proj * pA[2]};
        const vec2 b2[3] = {proj * pB[0], proj * pB[1], proj * pB[2]};
        std::vector<vec2> poly = {b2[0], b2[1], b2[2]};
        for (int k = 0; k < 3 && !poly.empty(); ++k)
          poly = ClipPolyByHalfplane(poly, a2[k], a2[(k + 1) % 3]);
        if (poly.size() < 3) return false;
        double area = 0.0, perim = 0.0;
        for (int k = 0; k < (int)poly.size(); ++k) {
          const vec2 &p = poly[k], &q = poly[(k + 1) % (int)poly.size()];
          area += p.x * q.y - q.x * p.y;
          perim += la::length(q - p);
        }
        return std::abs(area) * 0.5 > perim * eps;
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
        // Degenerate contact: record the endpoint x's as criticals (no vert
        // identity needed). SubEpsFeature guard in stage C' will catch any
        // macro-scale hazard.
        arr.criticalXs.push_back(qA.x);
        arr.criticalXs.push_back(qB.x);
        ++cnt.subEpsContactsDropped;
        continue;
      }

      const int vA = FindOrAddVert(arr.verts, qA, eps);
      const int vB = FindOrAddVert(arr.verts, qB, eps);
      if (vA == vB) {
        // Collapsed after snap: same degenerate-contact handling as above -
        // the endpoint x's still enter the critical set.
        arr.criticalXs.push_back(qA.x);
        arr.criticalXs.push_back(qB.x);
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
      arr.criticalXs.push_back(xCross);
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// Stage D': strip emission.
// Track extension evaluates a piece endpoint's (y,z) at a boundary critical;
// the extended limits feed the cap arrangements (stage E'), whose output
// chains the strips then zip (spec [R2-fold]: strips consume the cap
// subdivision, they do not re-derive it).
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
//
// The class test is an eps-match over the slab's own tracks because the
// engine capture carries bare (y,z) endpoints - piece endpoints have no
// provenance.  Adjudicated in spec [R2-fold]: threading per-endpoint
// provenance through the engine's split/merge machinery is heavier than
// matching against the slab's bounded seam-track set, and the classes are
// mutually exclusive at scale > eps (a seam crossing interior to a face
// segment is > eps from its endpoints, else the block rule collapsed it).
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

// The third consumer of the cap arrangement (spec E' [R2-fold]): its output
// verts subdivide the adjacent strip edges.  `arrVerts` is the arrangement's
// FULL vert set (merged inputs + every collected-arrangement vertex, per the
// engine's negEdges contract).  Returns the verts lying strictly interior to
// the chord [p0, p1] (within eps - the tolerance-model on-edge test, same
// posture as the engine's own incidence rule), sorted by projection parameter
// (exact ties by lex order).  The returned POSITIONS become the strip edge's
// chain verts, bitwise equal to the cap triangulation corners they pair with.
static std::vector<vec2> ChainSplitVerts(const std::vector<vec2>& arrVerts,
                                         vec2 p0, vec2 p1, double eps) {
  const vec2 segD = p1 - p0;
  const double segLen2 = la::dot(segD, segD);
  if (segLen2 <= eps * eps) return {};
  const double invLen2 = 1.0 / segLen2;
  const double tEps = eps / std::sqrt(segLen2);
  std::vector<std::pair<double, vec2>> hits;
  for (const auto& v : arrVerts) {
    const double t = la::dot(v - p0, segD) * invLen2;
    if (t <= tEps || t >= 1.0 - tEps) continue;
    const vec2 proj = p0 + t * segD;
    if (la::length(v - proj) > eps) continue;
    hits.push_back({t, v});
  }
  std::sort(
      hits.begin(), hits.end(),
      [](const std::pair<double, vec2>& a, const std::pair<double, vec2>& b) {
        if (a.first != b.first) return a.first < b.first;
        if (a.second.x != b.second.x) return a.second.x < b.second.x;
        return a.second.y < b.second.y;
      });
  std::vector<vec2> out;
  out.reserve(hits.size());
  for (const auto& h : hits) out.push_back(h.second);
  return out;
}

// Zipper-triangulate a strip between two chains (position polylines): a at
// x=xLo, b at x=xHi, each ordered along the strip edge.  A single-vert chain
// is a piece that vanished at that critical; the strip degenerates to a fan
// from the point.  Chain verts are emitted bitwise, so strip corners coincide
// exactly with cap triangulation corners.  The merge advances by projection
// parameter on each chain's own chord; exact ties advance the xLo side first
// (the other diagonal of the same quad).
static void ZipperEmit(double xLo, double xHi, const std::vector<vec2>& a,
                       const std::vector<vec2>& b, std::vector<OutTri3D>& out) {
  const int m = (int)a.size(), n = (int)b.size();
  auto params = [](const std::vector<vec2>& c) {
    std::vector<double> t(c.size(), 0.0);
    const vec2 d = c.back() - c.front();
    const double len2 = la::dot(d, d);
    if (c.size() >= 2 && len2 > 0.0) {
      for (int i = 1; i + 1 < (int)c.size(); ++i)
        t[i] = la::dot(c[i] - c.front(), d) / len2;
      t.back() = 1.0;
    }
    return t;
  };
  auto emitTri = [&](vec3 v0, vec3 v1, vec3 v2) {
    if (v0 != v1 && v1 != v2 && v0 != v2) out.push_back({v0, v1, v2});
  };
  const std::vector<double> ta = params(a), tb = params(b);
  int i = 0, j = 0;
  while (i + 1 < m || j + 1 < n) {
    const bool advanceA = i + 1 < m && (j + 1 >= n || ta[i + 1] <= tb[j + 1]);
    if (advanceA) {
      emitTri({xLo, a[i].x, a[i].y}, {xLo, a[i + 1].x, a[i + 1].y},
              {xHi, b[j].x, b[j].y});
      ++i;
    } else {
      emitTri({xLo, a[i].x, a[i].y}, {xHi, b[j + 1].x, b[j + 1].y},
              {xHi, b[j].x, b[j].y});
      ++j;
    }
  }
}

// Per built slab: the strip-edge chains at each boundary critical, indexed
// like slab.pieces.  Written by the cap at the corresponding critical.
struct StripChains {
  std::vector<std::vector<vec2>> lo, hi;
};

// Stage D': emit strips.  Each strip zips its two c-side chains, built by the
// cap arrangements at its slab's boundary-pair canonical criticals (spec E'
// [R2-fold] third consumer).  No geometry is computed here: the chains carry
// the cap arrangements' vert positions bitwise.  Every built slab has both
// sides bound (its own xHi is always its pair's canonical; its lo is bound at
// the preceding gap's canonical), with one chain per piece - attribution
// failures fail closed in EmitCaps before this runs.
static std::optional<std::pair<FatalReason, std::string>> EmitStrips(
    const std::vector<SlabResult>& slabs,
    const std::vector<StripChains>& chains, std::vector<OutTri3D>& out) {
  for (int si = 0; si < (int)slabs.size(); ++si) {
    if (!slabs[si].built) continue;
    const StripChains& ch = chains[si];
    // One non-empty chain per retained piece on each side is the EmitCaps
    // contract (every built slab's sides bind at their pair-canonical
    // criticals); violation means emission would drop or mispair strips -
    // fail closed.  Both sides are checked against the piece count itself,
    // not just each other, so equal truncation cannot slip through.
    const size_t nPieces = slabs[si].pieces.size();
    if (ch.lo.size() != nPieces || ch.hi.size() != nPieces) {
      DEBUG_ASSERT(false, logicErr,
                   "strip chain count disagrees with piece count");
      return std::make_pair(
          FatalReason::NonManifoldEmission,
          std::string("strip chain count disagrees with piece count"));
    }
    for (size_t k = 0; k < ch.lo.size(); ++k) {
      if (ch.lo[k].empty() || ch.hi[k].empty()) {
        DEBUG_ASSERT(false, logicErr, "empty strip chain");
        return std::make_pair(FatalReason::NonManifoldEmission,
                              std::string("empty strip chain"));
      }
      ZipperEmit(slabs[si].xLo, slabs[si].xHi, ch.lo[k], ch.hi[k], out);
    }
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Stage E': cap emission.
// At each critical x=c, ONE 2D arrangement runs over the extended pieces of
// both adjacent slabs, and everything downstream consumes it (spec [R2-fold]
// "one arrangement per critical, three consumers"): the positive measure's
// boundary is cap_plus (L - R), the negated measure's is cap_minus (R - L),
// and the arrangement's retained verts subdivide both adjacent slabs'
// strip-edge chains.  One vertex snap, one crossing discovery - cap corners
// and strip chain verts are the same positions bitwise.
// ---------------------------------------------------------------------------

// Triangulate a Polygons and append cap triangles to out, emitting at x=xCap.
// If flipWinding is true, reverse each triangle (for -x normal caps).
// idx is global across all loops (flat concatenation of allVerts).
// Returns false if the triangulator produced an invalid index (fail closed).
static bool TriangulateCap(const Polygons& polys, double xCap, bool flipWinding,
                           double eps, std::vector<OutTri3D>& out) {
  if (polys.empty()) return true;

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

  const int nAll = (int)allVerts.size();
  for (const auto& t : tris) {
    const bool inRange = t.x >= 0 && t.x < nAll && t.y >= 0 && t.y < nAll &&
                         t.z >= 0 && t.z < nAll;
    DEBUG_ASSERT(inRange, logicErr,
                 "cap triangulation returned out-of-range index");
    if (!inRange) return false;  // fail closed: no silent cap loss
    const vec3 v0 = {xCap, allVerts[t.x].x, allVerts[t.x].y};
    const vec3 v1 = {xCap, allVerts[t.y].x, allVerts[t.y].y};
    const vec3 v2 = {xCap, allVerts[t.z].x, allVerts[t.z].y};
    if (!flipWinding)
      out.push_back({v0, v1, v2});
    else
      out.push_back({v0, v2, v1});
  }
  return true;
}

// Extended cap input (spec [R2-fold] step (1)): both slabs extend their
// retained pieces to x=c.  Endpoint slots align 1:1 with slab.pieces.  A
// piece whose source face has no track is a broken attribution; ok goes
// false and the caller fails closed (BuildSlabs made id conflicts fatal, so
// this is unreachable in a consistent pipeline - asserted AND checked).
// NO vertex pre-dedup: the cap arrangement's MergeVerts owns snapping.
struct CapEdgeSet {
  std::vector<vec2> rawVerts;
  std::vector<std::pair<int, int>> lPieces, rPieces;  // rawVerts index pairs
  bool ok = true;
};

static CapEdgeSet BuildCapEdgeSet(const SlabResult* leftSlab,
                                  const SlabResult* rightSlab, double xCap,
                                  double eps) {
  CapEdgeSet ces;
  auto addPieces = [&](const SlabResult& slab,
                       std::vector<std::pair<int, int>>& slots) {
    std::map<int, const FaceTrack*> trackOf;
    for (const auto& tk : slab.faceTracks) trackOf[tk.faceId] = &tk;
    slots.reserve(slab.pieces.size());
    for (const auto& piece : slab.pieces) {
      const auto it = trackOf.find(piece.sourceId);
      DEBUG_ASSERT(it != trackOf.end(), logicErr,
                   "retained piece has no face track");
      if (it == trackOf.end()) {
        ces.ok = false;
        return;
      }
      // Piece is emission-oriented (interior-on-left) per spec engine contract.
      const vec2 eFrom = ExtendPtWithSeams(piece.from, *it->second,
                                           slab.seamTracks, xCap, eps);
      const vec2 eTo =
          ExtendPtWithSeams(piece.to, *it->second, slab.seamTracks, xCap, eps);
      const int i0 = (int)ces.rawVerts.size();
      ces.rawVerts.push_back(eFrom);
      ces.rawVerts.push_back(eTo);
      slots.push_back({i0, i0 + 1});
    }
  };
  // L verts appended first, then R: one deterministic input order for the ONE
  // arrangement.
  if (leftSlab && leftSlab->built) addPieces(*leftSlab, ces.lPieces);
  if (ces.ok && rightSlab && rightSlab->built)
    addPieces(*rightSlab, ces.rPieces);
  return ces;
}

// The cap at x=xCap: ONE arrangement, three consumers (spec [R2-fold]).
// L pieces enter at +1, R pieces at -1; the positive measure's boundary is
// cap_plus (L - R), the negated measure's is cap_minus (R - L), and the
// retained verts subdivide the adjacent strip-edge chains.  `leftChains` /
// `rightChains`, when bound, receive one chain per slab piece (positions from
// r.verts bitwise); a single-vert chain is a piece that vanished at this
// critical.
// Returns nullopt on success, or the fatal (reason, detail) to propagate.
static std::optional<std::pair<FatalReason, std::string>> ComputeCap(
    const SlabResult* leftSlab, const SlabResult* rightSlab, double xCap,
    double eps, std::vector<OutTri3D>& out,
    std::vector<std::vector<vec2>>* leftChains,
    std::vector<std::vector<vec2>>* rightChains, Overlap3Counters& cnt) {
  const CapEdgeSet ces = BuildCapEdgeSet(leftSlab, rightSlab, xCap, eps);
  if (!ces.ok)
    return std::make_pair(FatalReason::EngineIdConflict,
                          std::string("retained piece has no face track"));
  // Pre-size bound chain outputs so every piece has a slot.
  if (leftChains) leftChains->assign(ces.lPieces.size(), {});
  if (rightChains) rightChains->assign(ces.rPieces.size(), {});
  if (ces.rawVerts.empty()) return std::nullopt;

  std::vector<EdgeM> edges;
  edges.reserve(ces.lPieces.size() + ces.rPieces.size());
  for (const auto& s : ces.lPieces) edges.push_back({s.first, s.second, +1});
  for (const auto& s : ces.rPieces) edges.push_back({s.first, s.second, -1});

  std::vector<OutEdge> negEdges;
  ++cnt.capArrangements;
  const OverlapResult r =
      RemoveOverlaps2D(ces.rawVerts, edges, eps, /*debug=*/false, WindRule::Add,
                       /*trace=*/nullptr, &negEdges);
  bool loopsClosed = true;
  bool trisOk = true;
  if (!r.edges.empty()) {
    const Polygons cp = OutEdgesToPolygons(r.verts, r.edges, &loopsClosed);
    trisOk = TriangulateCap(cp, xCap, /*flipWinding=*/false, eps, out);
  }
  if (loopsClosed && trisOk && !negEdges.empty()) {
    const Polygons cm = OutEdgesToPolygons(r.verts, negEdges, &loopsClosed);
    trisOk = TriangulateCap(cm, xCap, /*flipWinding=*/true, eps, out);
  }
  if (!loopsClosed)
    return std::make_pair(FatalReason::NonManifoldEmission,
                          std::string("cap boundary walk failed to close"));
  if (!trisOk)
    return std::make_pair(
        FatalReason::NonManifoldEmission,
        std::string("cap triangulation returned invalid indices"));

  if (!leftChains && !rightChains) return std::nullopt;
  // The subdivision authority for strip edges is the arrangement's FULL vert
  // set (r.verts: merged input verts plus every collected-arrangement vertex,
  // per the engine's negEdges contract).  It cannot be narrowed to
  // retained-edge verts - where L and R coincide, their edges annihilate and
  // retention is empty, yet the strips on both sides still need the SAME
  // subdivision to pair across the critical; the merge/incidence machinery is
  // exactly what guarantees both sides see the same verts.
  auto buildChains = [&](const std::vector<std::pair<int, int>>& slots,
                         std::vector<std::vector<vec2>>& chains) {
    for (size_t k = 0; k < slots.size(); ++k) {
      const int m0 = r.inputVert2Merged[slots[k].first];
      const int m1 = r.inputVert2Merged[slots[k].second];
      std::vector<vec2>& chain = chains[k];
      chain.push_back(r.verts[m0]);
      if (m0 != m1) {
        for (const vec2& v :
             ChainSplitVerts(r.verts, r.verts[m0], r.verts[m1], eps))
          chain.push_back(v);
        chain.push_back(r.verts[m1]);
      }
    }
  };
  if (leftChains) buildChains(ces.lPieces, *leftChains);
  if (rightChains) buildChains(ces.rPieces, *rightChains);
  return std::nullopt;
}

// Stage E' driver: one cap per critical (spec [R3-fold]: runs are never
// merged; only IEEE-exact duplicate criticals collapse).  Adjacent built
// slabs are found by INDEX: slab bounds are these exact critical values by
// construction, so no tolerance enters the lookup.
//
// Chains bind per adjacent-built-slab PAIR, at the pair's canonical critical:
// the first critical of the gap between the slabs (index li+1, which is the
// critical itself under direct adjacency - the generic case).  Both sides of
// the pair consume that ONE arrangement, so the strip-strip weld across a
// sub-eps run pairs exactly: the two arrangements of a sub-eps critical pair
// can disagree macroscopically about the subdivision of a shared (cancelled)
// edge - geometry changes discontinuously at a critical by definition - so
// one arrangement must own the seam.  The run's other caps are emitted from
// their own arrangements but are sub-eps slivers that collapse in the
// assembly weld.
static std::optional<std::pair<FatalReason, std::string>> EmitCaps(
    const std::vector<SlabResult>& slabs, const std::vector<double>& crits,
    double eps, std::vector<OutTri3D>& out, std::vector<StripChains>& chains,
    Overlap3Counters& cnt) {
  const int nSlabs = (int)slabs.size();
  for (int ci = 0; ci < (int)crits.size(); ++ci) {
    int li = ci - 1;
    while (li >= 0 && !slabs[li].built) --li;
    int ri = ci;
    while (ri < nSlabs && !slabs[ri].built) ++ri;
    const SlabResult* left = li >= 0 ? &slabs[li] : nullptr;
    const SlabResult* right = ri < nSlabs ? &slabs[ri] : nullptr;
    const bool canonical = ci == li + 1;
    if (auto fatal =
            ComputeCap(left, right, crits[ci], eps, out,
                       (canonical && left) ? &chains[li].hi : nullptr,
                       (canonical && right) ? &chains[ri].lo : nullptr, cnt))
      return fatal;
  }
  return std::nullopt;
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
  // Duplicate: per-critical caps (spec [R3-fold]: runs are never merged) mean
  // a sub-eps critical pair computes the SAME macro difference twice - at
  // x=c1 and x=c2 with |c2-c1| <= eps - and both triangulations collapse to
  // identical vertex triples after the weld.  This exact-duplicate drop is
  // the [R3] sentence "their differences are empty" realized empirically: it
  // suppresses only post-weld identical triangles, never distinct caps (the
  // run-merge mistake this replaced suppressed by slab-pair identity).
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

  // The criticals are the slab boundaries: the same sorted-unique vert
  // x-values BuildSlabs derived the slabs from, recovered exactly (no
  // recompute, no tolerance).
  std::vector<double> crits;
  if (!slabs.empty()) {
    crits.reserve(slabs.size() + 1);
    for (const auto& s : slabs) crits.push_back(s.xLo);
    crits.push_back(slabs.back().xHi);
  }

  // E' then D': caps run first; each cap's ONE arrangement is the source of
  // truth at its critical - cap_plus, cap_minus, and both adjacent slabs'
  // strip chains all consume it (spec [R2-fold]).
  std::vector<OutTri3D> emitted;
  std::vector<StripChains> chains(slabs.size());
  if (auto capFatal = EmitCaps(slabs, crits, eps, emitted, chains, cnt)) {
    result.fatal = capFatal->first;
    result.detail = std::move(capFatal->second);
    result.counters = cnt;
    return result;
  }
  if (auto stripFatal = EmitStrips(slabs, chains, emitted)) {
    result.fatal = stripFatal->first;
    result.detail = std::move(stripFatal->second);
    result.counters = cnt;
    return result;
  }

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
