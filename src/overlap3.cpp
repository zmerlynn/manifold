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
// Design: docs/SweepPlane3D.md (crucible round 3).

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

// Stage C forward declaration.
StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt);

namespace {

// ---------------------------------------------------------------------------
// Shared utilities
// ---------------------------------------------------------------------------

static bool IsInside3D(int64_t w) { return w > 0; }

static double Dist3(vec3 a, vec3 b) { return la::length(a - b); }

static double PointSegDist3(vec3 p, vec3 a, vec3 b) {
  const vec3 ab = b - a;
  const double len2 = la::dot(ab, ab);
  if (len2 == 0.0) return la::length(p - a);
  const double t = std::max(0.0, std::min(1.0, la::dot(p - a, ab) / len2));
  return la::length(p - (a + t * ab));
}

static double SignedArea2D(const std::vector<vec2>& p) {
  const int n = (int)p.size();
  double a = 0.0;
  for (int i = 0; i < n; ++i) {
    const vec2 &u = p[i], &v = p[(i + 1) % n];
    a += u.x * v.y - v.x * u.y;
  }
  return a * 0.5;
}

static bool PointInPoly2D(vec2 q, const std::vector<vec2>& p) {
  const int n = (int)p.size();
  int w = 0;
  for (int i = 0; i < n; ++i) {
    const vec2 &a = p[i], &b = p[(i + 1) % n];
    if (a.y <= q.y) {
      if (b.y > q.y && la::cross(b - a, q - a) > 0) ++w;
    } else {
      if (b.y <= q.y && la::cross(b - a, q - a) < 0) --w;
    }
  }
  return w != 0;
}

static std::vector<int> CollapseSpikes(std::vector<int> loop) {
  bool changed = true;
  while (changed && (int)loop.size() >= 3) {
    changed = false;
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i) {
      if (loop[(i - 1 + n) % n] == loop[(i + 1) % n] &&
          loop[(i - 1 + n) % n] != loop[i]) {
        loop.erase(loop.begin() + i);
        changed = true;
        break;
      }
    }
    if ((int)loop.size() < 3) break;
    const int m = (int)loop.size();
    bool hasDup = false;
    for (int i = 0; i < m; ++i)
      if (loop[i] == loop[(i + 1) % m]) {
        hasDup = true;
        break;
      }
    if (hasDup) {
      std::vector<int> nxt;
      nxt.reserve(m);
      for (int i = 0; i < m; ++i)
        if (loop[i] != loop[(i + 1) % m]) nxt.push_back(loop[i]);
      loop = nxt;
      changed = true;
    }
  }
  return loop;
}

// Point-winding query (spec "AXIS-PARALLEL FACES").
static int64_t PointWinding2D(const std::vector<SweepCapture>& pieces,
                              double qy, double qz) {
  int64_t w = 0;
  for (const auto& pc : pieces) {
    const double y0 = pc.from.x, y1 = pc.to.x;
    const double z0 = pc.from.y, z1 = pc.to.y;
    if (y0 == y1) continue;
    const double ylo = std::min(y0, y1), yhi = std::max(y0, y1);
    if (qy <= ylo || qy > yhi) continue;
    const double zAtQ = z0 + (qy - y0) / (y1 - y0) * (z1 - z0);
    if (zAtQ <= qz) continue;
    w += (pc.below - pc.above);
  }
  return w;
}

// Minimum 2D distance from point to polygon boundary.
static double DistToBoundary2D(vec2 p, const std::vector<vec2>& poly) {
  const int n = (int)poly.size();
  double minD = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    const vec2 a = poly[i], b = poly[(i + 1) % n];
    const vec2 ab = b - a;
    const double len2 = la::dot(ab, ab);
    const vec2 closest =
        (len2 == 0.0)
            ? a
            : a + std::max(0.0, std::min(1.0, la::dot(p - a, ab) / len2)) * ab;
    minD = std::min(minD, la::length(p - closest));
  }
  return minD;
}

// ---------------------------------------------------------------------------
// Stage A: merge verts, canonicalize faces
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
// Stage B helpers
// ---------------------------------------------------------------------------

static bool LineTriClip(vec3 P, vec3 D, vec3 v0, vec3 v1, vec3 v2, vec3 n,
                        double& tLo, double& tHi, double eps) {
  tLo = -1e18;
  tHi = 1e18;
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

// Exact Dlen==0.0 check: coplanar was already handled by planeSep<=eps before
// call.
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
  const double tLo = std::max(tAlo, tBlo), tHi2 = std::min(tAhi, tBhi);
  if (tLo >= tHi2) return false;
  qA = P + tLo * dir;
  qB = P + tHi2 * dir;
  return true;
}

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

static int FindOrAddVert(std::vector<MergedVert>& verts, vec3 pos, double eps) {
  for (int i = 0; i < (int)verts.size(); ++i)
    if (la::length(verts[i].pos - pos) <= eps) return i;
  const int id = (int)verts.size();
  verts.push_back({pos});
  return id;
}

// ---------------------------------------------------------------------------
// Stage B: seam computation, contact graph, triple unification, subdEdges
// ---------------------------------------------------------------------------

struct StageBResult {
  ArrangementGeometry arr;
  std::optional<FatalReason> fatal;
  std::string detail;
};

static StageBResult StageB(const StageAResult& stageA, double eps,
                           Overlap3Counters& cnt) {
  StageBResult res;
  ArrangementGeometry& arr = res.arr;
  const int nFaces = (int)stageA.faces.size();

  arr.verts.resize(stageA.mergedVerts.size());
  for (int i = 0; i < (int)stageA.mergedVerts.size(); ++i)
    arr.verts[i] = {stageA.mergedVerts[i]};
  arr.faces = stageA.faces;
  arr.faceSeams.resize(nFaces);

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

  // Contact-graph: raw 3D points from degenerate intersections.
  // The union-find below handles clustering. No dedup here.
  std::vector<vec3> contactPts;

  // Triple-point candidates with seam parameters for deterministic write-back.
  struct TripleCand {
    vec3 pos;
    int seamA, seamB;
    double tA, tB;
  };
  std::vector<TripleCand> tripleCands;

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

      // Coplanar detection: ALL THREE vertices of FB within eps of FA's plane.
      // Checking only one vertex (pb0) is insufficient - it fires for
      // EdgeInPlane configurations (where 1-2 verts of fj are on fi's plane)
      // and for near-parallel faces where individual vertices happen to be
      // close.
      const double planeSep0 = std::abs(la::dot(na, pb0) - planeDist[fi]);
      const double planeSep1 = std::abs(la::dot(na, pb1) - planeDist[fi]);
      const double planeSep2 = std::abs(la::dot(na, pb2) - planeDist[fi]);
      if (planeSep0 <= eps && planeSep1 <= eps && planeSep2 <= eps) {
        const mat2x3 proj = GetAxisAlignedProjection(na);
        const vec2 a2[3] = {proj * pa0, proj * pa1, proj * pa2};
        const vec2 b2[3] = {proj * pb0, proj * pb1, proj * pb2};
        std::vector<vec2> poly = {b2[0], b2[1], b2[2]};
        for (int k = 0; k < 3 && !poly.empty(); ++k)
          poly = ClipPolyByHalfplane(poly, a2[k], a2[(k + 1) % 3]);
        double clipArea = 0.0;
        for (int k = 0; k < (int)poly.size(); ++k) {
          const vec2 &p = poly[k], &q = poly[(k + 1) % (int)poly.size()];
          clipArea += p.x * q.y - q.x * p.y;
        }
        if (std::abs(clipArea) * 0.5 > eps * eps) {
          res.fatal = FatalReason::CoplanarOverlap;
          res.detail = "coplanar face pair with interior overlap";
          return res;
        }
        continue;
      }
      // Also check if all three vertices of FA are within eps of FB's plane
      // (symmetric coplanar check - handles the reversed orientation case).
      {
        const double dB = la::dot(nb, pb0);
        const double sepA0 = std::abs(la::dot(nb, pa0) - dB);
        const double sepA1 = std::abs(la::dot(nb, pa1) - dB);
        const double sepA2 = std::abs(la::dot(nb, pa2) - dB);
        if (sepA0 <= eps && sepA1 <= eps && sepA2 <= eps) {
          const mat2x3 proj = GetAxisAlignedProjection(nb);
          const vec2 a2[3] = {proj * pa0, proj * pa1, proj * pa2};
          const vec2 b2[3] = {proj * pb0, proj * pb1, proj * pb2};
          std::vector<vec2> poly = {a2[0], a2[1], a2[2]};
          for (int k = 0; k < 3 && !poly.empty(); ++k)
            poly = ClipPolyByHalfplane(poly, b2[k], b2[(k + 1) % 3]);
          double clipArea = 0.0;
          for (int k = 0; k < (int)poly.size(); ++k) {
            const vec2 &p = poly[k], &q = poly[(k + 1) % (int)poly.size()];
            clipArea += p.x * q.y - q.x * p.y;
          }
          if (std::abs(clipArea) * 0.5 > eps * eps) {
            res.fatal = FatalReason::CoplanarOverlap;
            res.detail = "coplanar face pair with interior overlap";
            return res;
          }
          continue;
        }
      }

      const vec3 Dcross = la::cross(na, nb);
      if (la::length(Dcross) == 0.0) continue;  // exact parallel, separated

      // EdgeInPlane detection (spec Stage B).
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

      vec3 qA, qB;
      if (!TriTriSeam(pa0, pa1, pa2, na, pb0, pb1, pb2, nb, qA, qB, eps))
        continue;
      const double seamLen = la::length(qB - qA);

      if (seamLen <= eps) {
        // Degenerate: add BOTH original computed endpoints (spec: "keep
        // original endpoints in the diameter graph").
        contactPts.push_back(qA);
        contactPts.push_back(qB);
        continue;
      }

      const int vA = FindOrAddVert(arr.verts, qA, eps);
      const int vB = FindOrAddVert(arr.verts, qB, eps);
      if (vA == vB) {
        // Collapsed after snapping: add original endpoints to contact graph.
        contactPts.push_back(qA);
        contactPts.push_back(qB);
        continue;
      }

      const int seamIdx = (int)arr.seams.size();
      arr.seams.push_back({fi, fj, {vA, vB}});
      arr.faceSeams[fi].push_back(seamIdx);
      arr.faceSeams[fj].push_back(seamIdx);
    }
  }

  // Contact-graph SubResolutionChain check.
  if (!contactPts.empty()) {
    const int nPts = (int)contactPts.size();
    DisjointSets duf(nPts);
    for (int i = 0; i < nPts; ++i)
      for (int j = i + 1; j < nPts; ++j)
        if (la::length(contactPts[i] - contactPts[j]) <= eps) duf.unite(i, j);
    std::map<int, std::vector<int>> comps;
    for (int i = 0; i < nPts; ++i) comps[(int)duf.find(i)].push_back(i);
    for (auto& [root, members] : comps) {
      double diam = 0.0;
      for (int a : members)
        for (int b : members)
          diam = std::max(diam, la::length(contactPts[a] - contactPts[b]));
      if (diam <= eps)
        ++cnt.subEpsContactsDropped;
      else {
        res.fatal = FatalReason::SubResolutionChain;
        res.detail = "sub-resolution seam chain";
        return res;
      }
    }
  }

  // Triple-point candidates: pairwise seam crossings in each face plane.
  // Store ALL without pre-merging (spec global unification).
  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = stageA.faces[fi];
    const mat2x3 proj = GetAxisAlignedProjection(face.normal);
    const auto& fseams = arr.faceSeams[fi];
    const int ns = (int)fseams.size();
    for (int i = 0; i < ns; ++i) {
      for (int j = i + 1; j < ns; ++j) {
        const Seam &sA = arr.seams[fseams[i]], &sB = arr.seams[fseams[j]];
        const vec3 a3s = arr.verts[sA.vertIds.front()].pos,
                   a3e = arr.verts[sA.vertIds.back()].pos;
        const vec3 b3s = arr.verts[sB.vertIds.front()].pos,
                   b3e = arr.verts[sB.vertIds.back()].pos;
        const vec2 a0 = proj * a3s, a1 = proj * a3e, b0 = proj * b3s,
                   b1 = proj * b3e;
        const vec2 da = a1 - a0, db = b1 - b0, dc = b0 - a0;
        const double den = la::cross(da, db);
        if (den == 0.0) continue;
        const double tA = la::cross(dc, db) / den, tB = la::cross(dc, da) / den;
        const vec3 q3d = a3s + tA * (a3e - a3s);
        // Validate with 3D Euclidean (metric table) instead of parameter slack.
        if (PointSegDist3(q3d, a3s, a3e) > eps) continue;
        if (PointSegDist3(q3d, b3s, b3e) > eps) continue;
        tripleCands.push_back({q3d, fseams[i], fseams[j], tA, tB});
      }
    }
  }

  // Global triple-point unification (spec Stage B).
  {
    std::vector<vec3> pool;
    std::vector<int> poolVert;  // >=0: seam vert index, -1: candidate
    std::set<int> seamVertSet;
    for (auto& s : arr.seams)
      for (int v : s.vertIds) seamVertSet.insert(v);
    for (int v : seamVertSet) {
      pool.push_back(arr.verts[v].pos);
      poolVert.push_back(v);
    }
    const int nSeamPts = (int)pool.size();
    for (auto& tc : tripleCands) {
      pool.push_back(tc.pos);
      poolVert.push_back(-1);
    }
    const int poolSz = (int)pool.size();

    if (poolSz > 0) {
      DisjointSets poolUF(poolSz);
      for (int i = 0; i < poolSz; ++i)
        for (int j = i + 1; j < poolSz; ++j)
          if (la::length(pool[i] - pool[j]) <= eps) poolUF.unite(i, j);

      std::map<int, std::vector<int>> clusters;
      for (int i = 0; i < poolSz; ++i)
        clusters[(int)poolUF.find(i)].push_back(i);

      std::map<int, int> clusterCanon;
      for (auto& [root, members] : clusters) {
        double diam = 0.0;
        for (int a : members)
          for (int b : members)
            diam = std::max(diam, la::length(pool[a] - pool[b]));
        if (diam > eps && members.size() > 1) {
          res.fatal = FatalReason::TripleDiameter;
          res.detail = "triple-point cluster diameter exceeds eps";
          return res;
        }
        // Centroid-nearest existing seam vert (spec convention).
        vec3 centroid(0.0);
        for (int m : members) centroid += pool[m];
        centroid /= (double)members.size();
        int bestVert = -1;
        double bestD = 1e30;
        for (int m : members) {
          if (poolVert[m] >= 0) {
            const double d = la::length(pool[m] - centroid);
            if (d < bestD) {
              bestD = d;
              bestVert = poolVert[m];
            }
          }
        }
        if (bestVert < 0) {
          bestVert = (int)arr.verts.size();
          arr.verts.push_back({centroid});
        }
        // Do NOT overwrite the representative's position (spec "pick and
        // retain").
        clusterCanon[root] = bestVert;
      }

      // Write-back: remap seam verts to canonical ids.
      for (int i = 0; i < nSeamPts; ++i) {
        const int canon = clusterCanon[(int)poolUF.find(i)];
        const int oldV = poolVert[i];
        if (canon == oldV) continue;
        for (auto& s : arr.seams)
          for (int& v : s.vertIds)
            if (v == oldV) v = canon;
      }

      // Insert triple candidates into incident seams by stored parameter.
      for (int ci = 0; ci < (int)tripleCands.size(); ++ci) {
        const int canon = clusterCanon[(int)poolUF.find(nSeamPts + ci)];
        const TripleCand& tc = tripleCands[ci];
        auto insertAt = [&](int sIdx, double t) {
          Seam& s = arr.seams[sIdx];
          for (int v : s.vertIds)
            if (v == canon) return;  // already present
          if (s.vertIds.size() < 2) return;
          const vec3& start = arr.verts[s.vertIds.front()].pos;
          const vec3& end = arr.verts[s.vertIds.back()].pos;
          const double segLen = la::length(end - start);
          int insertPos = 1;
          for (int k = 0; k + 1 < (int)s.vertIds.size(); ++k) {
            const double tk =
                (segLen > 0)
                    ? la::length(arr.verts[s.vertIds[k]].pos - start) / segLen
                    : 0.0;
            if (t > tk) insertPos = k + 1;
          }
          s.vertIds.insert(s.vertIds.begin() + insertPos, canon);
        };
        insertAt(tc.seamA, tc.tA);
        insertAt(tc.seamB, tc.tB);
      }
    }
  }

  // Populate arr.subdEdges: per-face-edge subdivision with event/triple verts.
  // Point-to-segment distance (metric table) for on-edge incidence.
  {
    const int nV = (int)arr.verts.size();
    for (int fi = 0; fi < nFaces; ++fi) {
      const CanonicalFace& face = arr.faces[fi];
      const int gv[3] = {face.verts.x, face.verts.y, face.verts.z};
      for (int ei = 0; ei < 3; ++ei) {
        const int gA = gv[ei], gB = gv[(ei + 1) % 3];
        const vec3 pA = arr.verts[gA].pos, pB = arr.verts[gB].pos;
        const vec3 ev = pB - pA;
        const double len2 = la::dot(ev, ev);
        std::vector<std::pair<double, int>> onEdge;
        for (int v = 0; v < nV; ++v) {
          if (v == gA || v == gB) continue;
          if (PointSegDist3(arr.verts[v].pos, pA, pB) > eps) continue;
          const double t =
              (len2 > 0) ? la::dot(arr.verts[v].pos - pA, ev) / len2 : 0.5;
          if (t <= 0.0 || t >= 1.0) continue;
          onEdge.push_back({t, v});
        }
        std::sort(onEdge.begin(), onEdge.end());
        onEdge.erase(std::unique(onEdge.begin(), onEdge.end(),
                                 [](const auto& a, const auto& b) {
                                   return a.second == b.second;
                                 }),
                     onEdge.end());
        SubdividedEdge se;
        se.faceId = fi;
        se.edgeIdx = ei;
        se.vertIds.push_back(gA);
        for (auto& [t, v] : onEdge) se.vertIds.push_back(v);
        se.vertIds.push_back(gB);
        arr.subdEdges.push_back(std::move(se));
      }
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// Stage D helpers: PSLG walk, classification, balance, propagation
// ---------------------------------------------------------------------------

// Planar half-edge graph (standard DCEL face-walk).
struct PlanarGraph {
  struct HE {
    int from, to;
  };
  std::vector<HE> he;
  std::vector<vec2> pos;
  std::vector<int> next_;

  void Build() {
    const int n = (int)he.size();
    std::vector<int> twin(n, -1);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (he[j].from == he[i].to && he[j].to == he[i].from) {
          twin[i] = j;
          break;
        }
    next_.assign(n, -1);
    for (int i = 0; i < n; ++i) {
      if (twin[i] < 0) continue;
      const int v = he[i].to;
      const vec2& base = pos[v];
      std::vector<std::pair<double, int>> out;
      for (int j = 0; j < n; ++j)
        if (he[j].from == v)
          out.push_back(
              {std::atan2((pos[he[j].to] - base).y, (pos[he[j].to] - base).x),
               j});
      std::sort(out.begin(), out.end());
      int idx = -1;
      for (int k = 0; k < (int)out.size(); ++k)
        if (out[k].second == twin[i]) {
          idx = k;
          break;
        }
      if (idx < 0) continue;
      const int m = (int)out.size();
      next_[i] = out[(idx - 1 + m) % m].second;
    }
  }

  std::vector<std::vector<int>> WalkLoops() const {
    const int n = (int)he.size();
    std::vector<bool> vis(n, false);
    std::vector<std::vector<int>> loops;
    for (int s = 0; s < n; ++s) {
      if (vis[s]) continue;
      std::vector<int> loop;
      int e = s, guard = n + 1;
      while (!vis[e] && guard-- > 0) {
        vis[e] = true;
        loop.push_back(he[e].from);
        if (next_[e] < 0) break;
        e = next_[e];
      }
      if ((int)loop.size() >= 3) loops.push_back(std::move(loop));
    }
    return loops;
  }
};

static bool Seg2DCross(vec2 p0, vec2 p1, vec2 q0, vec2 q1) {
  const vec2 dp = p1 - p0, dq = q1 - q0, dc = q0 - p0;
  const double den = la::cross(dp, dq);
  if (den == 0.0) return false;
  const double tP = la::cross(dc, dq) / den, tQ = la::cross(dc, dp) / den;
  return tP > 0.0 && tP < 1.0 && tQ > 0.0 && tQ < 1.0;
}

static std::optional<FatalReason> ValidateFacePSLG(
    int fi, const ArrangementGeometry& arr) {
  const CanonicalFace& face = arr.faces[fi];
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);
  auto p2 = [&](int vid) -> vec2 { return proj * arr.verts[vid].pos; };
  const auto& sIdxs = arr.faceSeams[fi];
  const int ns = (int)sIdxs.size();
  for (int i = 0; i < ns; ++i) {
    const Seam& sA = arr.seams[sIdxs[i]];
    for (int j = i + 1; j < ns; ++j) {
      const Seam& sB = arr.seams[sIdxs[j]];
      for (int a = 0; a + 1 < (int)sA.vertIds.size(); ++a) {
        const int vA0 = sA.vertIds[a], vA1 = sA.vertIds[a + 1];
        for (int b = 0; b + 1 < (int)sB.vertIds.size(); ++b) {
          const int vB0 = sB.vertIds[b], vB1 = sB.vertIds[b + 1];
          if (vA0 == vB0 || vA0 == vB1 || vA1 == vB0 || vA1 == vB1) continue;
          if (Seg2DCross(p2(vA0), p2(vA1), p2(vB0), p2(vB1)))
            return FatalReason::PSLGInvalid;
        }
      }
    }
  }
  return std::nullopt;
}

// Build face PSLG and walk to get PSLGRegion list.
// Uses arr.subdEdges for boundary subdivision (spec "On-edge incidences").
static std::vector<PSLGRegion> BuildFacePSLG(int fi,
                                             const ArrangementGeometry& arr,
                                             double eps) {
  const CanonicalFace& face = arr.faces[fi];
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);

  // Local vert index and 2D positions.
  std::vector<int> gVerts;     // local id -> global vert id
  std::map<int, int> localOf;  // global vert id -> local id
  auto ensureVert = [&](int gv) -> int {
    auto it = localOf.find(gv);
    if (it != localOf.end()) return it->second;
    const int id = (int)gVerts.size();
    localOf[gv] = id;
    gVerts.push_back(gv);
    return id;
  };

  // Canonical half-edge accumulator (no duplicates by construction).
  std::set<std::pair<int, int>> heSet;
  std::vector<std::pair<int, int>> heList;
  auto addHE = [&](int gFrom, int gTo) {
    const int lf = ensureVert(gFrom), lt = ensureVert(gTo);
    if (lf == lt) return;
    if (heSet.insert({lf, lt}).second) heList.push_back({lf, lt});
    if (heSet.insert({lt, lf}).second) heList.push_back({lt, lf});
  };

  // Boundary edges from subdEdges.
  for (const auto& se : arr.subdEdges)
    if (se.faceId == fi)
      for (int k = 0; k + 1 < (int)se.vertIds.size(); ++k)
        addHE(se.vertIds[k], se.vertIds[k + 1]);

  // Interior seam edges.
  for (int sIdx : arr.faceSeams[fi]) {
    const Seam& s = arr.seams[sIdx];
    for (int k = 0; k + 1 < (int)s.vertIds.size(); ++k)
      addHE(s.vertIds[k], s.vertIds[k + 1]);
  }

  if (heList.empty()) return {};

  // Build 2D positions now that all verts are known.
  std::vector<vec2> localPos(gVerts.size());
  for (int i = 0; i < (int)gVerts.size(); ++i)
    localPos[i] = proj * arr.verts[gVerts[i]].pos;

  PlanarGraph graph;
  graph.pos = localPos;
  graph.he.reserve(heList.size());
  for (auto [f, t] : heList) graph.he.push_back({f, t});
  graph.Build();
  auto rawLoops = graph.WalkLoops();

  // Separate CCW and CW loops; build CCW (sub-regions) and CW (candidates).
  struct LD {
    std::vector<int> gverts;
    std::vector<vec2> poly;
  };
  std::vector<LD> ccwLoops;
  std::vector<LD> cwLoops;

  for (const auto& raw : rawLoops) {
    std::vector<int> loop = CollapseSpikes(raw);
    if ((int)loop.size() < 3) continue;
    std::vector<int> gv;
    std::vector<vec2> poly;
    for (int lv : loop) {
      gv.push_back(gVerts[lv]);
      poly.push_back(localPos[lv]);
    }
    if ((int)gv.size() < 3) continue;
    const double area = SignedArea2D(poly);
    if (std::abs(area) < eps * eps) continue;  // sub-eps-area: drop
    if (area > 0)
      ccwLoops.push_back({gv, poly});
    else
      cwLoops.push_back({gv, poly});
  }

  // Drop the outer face CW loop: it is the CW traversal of the face's outer
  // boundary and has the largest absolute area of any CW loop. The outer face
  // wraps the face's exterior; its vertex centroid coincidentally lies inside
  // the CCW region (shared vertices), so centroid-based assignment would
  // incorrectly treat it as a hole. Removing it by max-area identification is
  // safe: genuine hole CW loops are sub-regions of a CCW face and always have
  // smaller absolute area than the outer face.
  if (!cwLoops.empty()) {
    int outerIdx = 0;
    double outerA = std::abs(SignedArea2D(cwLoops[0].poly));
    for (int i = 1; i < (int)cwLoops.size(); ++i) {
      const double a = std::abs(SignedArea2D(cwLoops[i].poly));
      if (a > outerA) {
        outerA = a;
        outerIdx = i;
      }
    }
    cwLoops.erase(cwLoops.begin() + outerIdx);
  }

  // Assign remaining CW loops (genuine holes) to the smallest containing CCW
  // region.
  std::vector<std::vector<std::vector<int>>> regHoles(ccwLoops.size());
  for (const auto& cw : cwLoops) {
    // Use edge-midpoint probe offset to LEFT (exterior of CW boundary) to find
    // the surrounding CCW region, avoiding the vertex-centroid confusion.
    // Compute centroid of the CW loop verts as a fallback sample.
    vec2 sample(0.0);
    for (int gv : cw.gverts) sample += proj * arr.verts[gv].pos;
    sample /= (double)cw.gverts.size();
    int bestIdx = -1;
    double bestArea = std::numeric_limits<double>::infinity();
    for (int ri = 0; ri < (int)ccwLoops.size(); ++ri) {
      if (PointInPoly2D(sample, ccwLoops[ri].poly)) {
        const double a = SignedArea2D(ccwLoops[ri].poly);
        if (a < bestArea) {
          bestArea = a;
          bestIdx = ri;
        }
      }
    }
    if (bestIdx >= 0) regHoles[bestIdx].push_back(cw.gverts);
    // else: not inside any CCW loop -> exterior complement, drop.
  }

  std::vector<PSLGRegion> regions;
  regions.reserve(ccwLoops.size());
  for (int ri = 0; ri < (int)ccwLoops.size(); ++ri) {
    PSLGRegion reg;
    reg.loopVerts = std::move(ccwLoops[ri].gverts);
    reg.holeVerts = std::move(regHoles[ri]);
    regions.push_back(std::move(reg));
  }
  return regions;
}

// Region classification via the widest covering slab's captured pieces.
static std::optional<std::pair<int64_t, int64_t>> ClassifyRegion(
    const PSLGRegion& region, int faceId, const std::vector<SlabResult>& slabs,
    const std::vector<MergedVert>& verts, const CanonicalFace& face, double eps,
    Overlap3Counters& cnt, std::optional<FatalReason>& outFatal) {
  double xMin = 1e18, xMax = -1e18;
  for (int v : region.loopVerts) {
    xMin = std::min(xMin, verts[v].pos.x);
    xMax = std::max(xMax, verts[v].pos.x);
  }
  for (const auto& h : region.holeVerts)
    for (int v : h) {
      xMin = std::min(xMin, verts[v].pos.x);
      xMax = std::max(xMax, verts[v].pos.x);
    }

  int bestSlab = -1;
  double bestW = -1.0;
  for (int si = 0; si < (int)slabs.size(); ++si) {
    const auto& sl = slabs[si];
    if (!sl.built) continue;
    if (sl.xLo >= xMin && sl.xHi <= xMax) {
      const double w = sl.xHi - sl.xLo;
      if (w > bestW ||
          (w == bestW && bestSlab >= 0 && sl.xLo < slabs[bestSlab].xLo)) {
        bestW = w;
        bestSlab = si;
      }
    }
  }
  if (bestSlab < 0 || bestW <= eps) return std::nullopt;

  const SlabResult& sl = slabs[bestSlab];
  const mat2x3 proj = GetAxisAlignedProjection(face.normal);
  std::vector<vec2> rPoly;
  for (int v : region.loopVerts) rPoly.push_back(proj * verts[v].pos);
  std::vector<std::vector<vec2>> hPolys;
  for (const auto& h : region.holeVerts) {
    std::vector<vec2> hp;
    for (int v : h) hp.push_back(proj * verts[v].pos);
    hPolys.push_back(std::move(hp));
  }

  std::optional<int64_t> belowVal, aboveVal;
  for (const auto& piece : sl.pieces) {
    if (piece.sourceId != faceId) continue;
    const vec2 mid2d = (piece.from + piece.to) * 0.5;
    const vec3 mid3d = {sl.xMid, mid2d.x, mid2d.y};
    const vec2 midProj = proj * mid3d;
    // Clearance: 2D face-plane Euclidean (metric table).
    if (DistToBoundary2D(midProj, rPoly) < eps) {
      ++cnt.clearanceSkips;
      continue;
    }
    bool tooClose = false;
    for (const auto& hp : hPolys)
      if (DistToBoundary2D(midProj, hp) < eps) {
        tooClose = true;
        break;
      }
    if (tooClose) {
      ++cnt.clearanceSkips;
      continue;
    }
    if (!PointInPoly2D(midProj, rPoly)) continue;
    {
      bool inHole = false;
      for (const auto& hp : hPolys)
        if (PointInPoly2D(midProj, hp)) {
          inHole = true;
          break;
        }
      if (inHole) continue;
    }
    if (!belowVal.has_value()) {
      belowVal = piece.below;
      aboveVal = piece.above;
    } else if (*belowVal != piece.below || *aboveVal != piece.above) {
      outFatal = FatalReason::ClassificationAmbiguity;
      return std::nullopt;
    }
  }
  if (!belowVal.has_value()) {
    outFatal = FatalReason::ClassificationAmbiguity;
    return std::nullopt;
  }
  return std::make_pair(*belowVal, *aboveVal);
}

// Axis-parallel face classification (spec "AXIS-PARALLEL FACES").
// "below" = winding in the nearest built slab to the left of xFace (from
// PointWinding2D). "above" = below + delta, where delta is the sum of winding
// contributions from ALL axis-parallel faces at xFace for which the probe is
// inside their (y,z) polygon. Each such face contributes -fi.mult *
// sign(fi.normal.x) (outward face normal.x>0 decreases winding as you exit;
// inward face normal.x<0 increases winding as you enter). Summing ALL faces
// handles the touching case (faces cancel: delta=0) and avoids querying the
// right slab's xMid which may be too far from xFace for faces near another
// solid's boundary.
static std::optional<std::pair<int64_t, int64_t>> ClassifyAxisParallel(
    const PSLGRegion& region, const CanonicalFace& face,
    const std::vector<SlabResult>& slabs, const std::vector<MergedVert>& verts,
    const std::vector<CanonicalFace>& allFaces, double xFace, double eps) {
  // Nearest built slab to the left: rightmost built slab with xHi <= xFace.
  // Seam vertices can produce sub-eps (unbuilt) slabs right next to xFace;
  // skip those.
  int leftSlab = -1;
  for (int si = 0; si < (int)slabs.size(); ++si) {
    if (!slabs[si].built) continue;
    if (slabs[si].xHi <= xFace) leftSlab = si;
  }
  if (leftSlab < 0) return std::nullopt;

  // Probe: centroid of largest triangle in region triangulation.
  // PointWinding2D operates in (y,z) section coordinates (piece.from.x = y,
  // piece.from.y = z). GetAxisAlignedProjection for negative-x normals flips y,
  // placing the probe outside the section and returning winding 0.
  // Fix: always use (y,z) from vertex positions for the winding query.
  auto secYZ = [&](int v) -> vec2 { return {verts[v].pos.y, verts[v].pos.z}; };

  // Build polygon in (y,z). Ensure CCW so TriangulateIdx doesn't throw.
  PolygonsIdx polys;
  {
    SimplePolygonIdx outer;
    for (int v : region.loopVerts) outer.push_back({secYZ(v), v});
    if ((int)outer.size() < 3) return std::nullopt;
    {
      std::vector<vec2> tmp;
      for (auto& pv : outer) tmp.push_back(pv.pos);
      if (SignedArea2D(tmp) < 0) std::reverse(outer.begin(), outer.end());
    }
    polys.push_back(std::move(outer));
  }
  for (const auto& h : region.holeVerts) {
    SimplePolygonIdx hole;
    for (int v : h) hole.push_back({secYZ(v), v});
    {
      std::vector<vec2> tmp;
      for (auto& pv : hole) tmp.push_back(pv.pos);
      if (SignedArea2D(tmp) > 0) std::reverse(hole.begin(), hole.end());
    }
    polys.push_back(std::move(hole));
  }
  const auto tris = TriangulateIdx(polys, -1.0);
  if (tris.empty()) return std::nullopt;

  double maxArea = -1.0;
  vec2 probe(0.0);
  for (const auto& t : tris) {
    const vec2 a = secYZ(t.x), b = secYZ(t.y), c = secYZ(t.z);
    const double area = std::abs(SignedArea2D({a, b, c}));
    if (area > maxArea) {
      maxArea = area;
      probe = (a + b + c) / 3.0;
    }
  }

  // Clearance from region boundary in (y,z) section space.
  std::vector<vec2> rPoly;
  for (int v : region.loopVerts) rPoly.push_back(secYZ(v));
  {
    const double clr = DistToBoundary2D(probe, rPoly);
    if (clr < eps) return std::nullopt;
  }

  const int64_t belowW =
      PointWinding2D(slabs[leftSlab].pieces, probe.x, probe.y);
  // Sum winding contributions from all axis-parallel faces at xFace.
  // Each contributes -fi.mult*sign(fi.normal.x) when the probe is inside its
  // (y,z) polygon.
  int64_t delta = 0;
  for (const auto& fi : allFaces) {
    const double x0 = verts[fi.verts.x].pos.x, x1 = verts[fi.verts.y].pos.x,
                 x2 = verts[fi.verts.z].pos.x;
    if (std::abs(x0 - xFace) > eps || std::abs(x1 - xFace) > eps ||
        std::abs(x2 - xFace) > eps)
      continue;
    const vec2 fv0 = {verts[fi.verts.x].pos.y, verts[fi.verts.x].pos.z};
    const vec2 fv1 = {verts[fi.verts.y].pos.y, verts[fi.verts.y].pos.z};
    const vec2 fv2 = {verts[fi.verts.z].pos.y, verts[fi.verts.z].pos.z};
    if (PointInPoly2D(probe, {fv0, fv1, fv2}))
      delta += -fi.mult * (fi.normal.x > 0 ? 1 : -1);
  }
  const int64_t aboveW = belowW + delta;
  return std::make_pair(belowW, aboveW);
}

// M7: seam balance check.
// 0-vs-nonzero guard retained: per-seam two-face counting cannot observe
// third-face pairing at triple intersections; the output manifold gate owns
// that class (spec crucible implementation records).
static bool CheckSeamBalance(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions,
    std::optional<FatalReason>& outFatal, std::string& detail) {
  auto hasEdgeNoSpike = [](const std::vector<int>& loop, int a, int b) -> bool {
    const int n = (int)loop.size();
    for (int i = 0; i < n; ++i) {
      const int cur = loop[i], nxt = loop[(i + 1) % n],
                prev = loop[(i - 1 + n) % n], nn = loop[(i + 2) % n];
      if (cur == a && nxt == b && prev != b && nn != a) return true;
      if (cur == b && nxt == a && prev != a && nn != b) return true;
    }
    return false;
  };
  auto countEdge = [&](const PSLGRegion& r, int a, int b) -> int {
    int c = hasEdgeNoSpike(r.loopVerts, a, b) ? 1 : 0;
    for (const auto& h : r.holeVerts) c += hasEdgeNoSpike(h, a, b) ? 1 : 0;
    return c;
  };
  auto isKept = [](const PSLGRegion& r) -> bool {
    return r.classified && (IsInside3D(r.below) != IsInside3D(r.above));
  };
  for (int si = 0; si < (int)arr.seams.size(); ++si) {
    const Seam& s = arr.seams[si];
    for (int k = 0; k + 1 < (int)s.vertIds.size(); ++k) {
      const int vA = s.vertIds[k], vB = s.vertIds[k + 1];
      int n0 = 0, n1 = 0;
      for (const auto& r : faceRegions[s.faceId0])
        if (isKept(r)) n0 += countEdge(r, vA, vB);
      for (const auto& r : faceRegions[s.faceId1])
        if (isKept(r)) n1 += countEdge(r, vA, vB);
      if (n0 > 0 && n1 > 0 && n0 != n1) {
        outFatal = FatalReason::BalanceViolation;
        detail = "seam " + std::to_string(si) + " edge (" + std::to_string(vA) +
                 "," + std::to_string(vB) + "): " + std::to_string(n0) + "/" +
                 std::to_string(n1);
        return false;
      }
    }
  }
  return true;
}

// ---------------------------------------------------------------------------
// M6: Anchor-component propagation
// ---------------------------------------------------------------------------

static bool RegionHasDirEdge(const PSLGRegion& r, int u, int v) {
  auto has = [](const std::vector<int>& lp, int u, int v) {
    const int n = (int)lp.size();
    for (int i = 0; i < n; ++i)
      if (lp[i] == u && lp[(i + 1) % n] == v) return true;
    return false;
  };
  if (has(r.loopVerts, u, v)) return true;
  for (const auto& h : r.holeVerts)
    if (has(h, u, v)) return true;
  return false;
}
static bool RegionTouchesEdge(const PSLGRegion& r, int u, int v) {
  return RegionHasDirEdge(r, u, v) || RegionHasDirEdge(r, v, u);
}

static std::optional<FatalReason> PropagateAnchorComponentsImpl(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt) {
  struct RN {
    int fi, ri;
  };
  std::vector<RN> nodes;
  const int nF = (int)faceRegions.size();
  std::vector<std::vector<int>> nodeOf(nF);
  for (int fi = 0; fi < nF; ++fi) {
    nodeOf[fi].assign(faceRegions[fi].size(), -1);
    for (int ri = 0; ri < (int)faceRegions[fi].size(); ++ri)
      if (faceRegions[fi][ri].degenerate && !faceRegions[fi][ri].classified)
        nodeOf[fi][ri] = (int)nodes.size(), nodes.push_back({fi, ri});
  }
  const int nN = (int)nodes.size();
  if (nN == 0) return std::nullopt;

  std::vector<int> ufP(nN);
  for (int i = 0; i < nN; ++i) ufP[i] = i;
  std::function<int(int)> ufFind = [&](int x) -> int {
    return ufP[x] == x ? x : ufP[x] = ufFind(ufP[x]);
  };
  auto ufUnite = [&](int a, int b) {
    a = ufFind(a);
    b = ufFind(b);
    if (a != b) ufP[a] = b;
  };

  for (int fi = 0; fi < nF; ++fi) {
    const int nr = (int)faceRegions[fi].size();
    for (int ra = 0; ra < nr; ++ra) {
      if (nodeOf[fi][ra] < 0) continue;
      const PSLGRegion& rA = faceRegions[fi][ra];
      for (int rb = ra + 1; rb < nr; ++rb) {
        if (nodeOf[fi][rb] < 0) continue;
        const PSLGRegion& rB = faceRegions[fi][rb];
        bool adj = false;
        auto chk = [&](const std::vector<int>& lp) {
          const int n = (int)lp.size();
          for (int k = 0; k < n && !adj; ++k)
            if (RegionHasDirEdge(rB, lp[(k + 1) % n], lp[k])) adj = true;
        };
        chk(rA.loopVerts);
        for (const auto& h : rA.holeVerts) chk(h);
        if (adj) ufUnite(nodeOf[fi][ra], nodeOf[fi][rb]);
      }
    }
  }
  for (const Seam& seam : arr.seams) {
    const int fa = seam.faceId0, fb = seam.faceId1;
    for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
      const int sv = seam.vertIds[k], sw = seam.vertIds[k + 1];
      for (int ra = 0; ra < (int)faceRegions[fa].size(); ++ra) {
        if (nodeOf[fa][ra] < 0 ||
            !RegionTouchesEdge(faceRegions[fa][ra], sv, sw))
          continue;
        for (int rb = 0; rb < (int)faceRegions[fb].size(); ++rb) {
          if (nodeOf[fb][rb] < 0 ||
              !RegionTouchesEdge(faceRegions[fb][rb], sv, sw))
            continue;
          ufUnite(nodeOf[fa][ra], nodeOf[fb][rb]);
        }
      }
    }
  }

  std::map<int, std::vector<int>> comps;
  for (int i = 0; i < nN; ++i) comps[ufFind(i)].push_back(i);

  for (auto& [root, members] : comps) {
    std::sort(members.begin(), members.end());
    std::set<int> vset;
    for (int nd : members) {
      const PSLGRegion& rg = faceRegions[nodes[nd].fi][nodes[nd].ri];
      for (int v : rg.loopVerts) vset.insert(v);
      for (const auto& h : rg.holeVerts)
        for (int v : h) vset.insert(v);
    }
    const std::vector<int> cverts(vset.begin(), vset.end());
    double diam = 0.0;
    for (int i = 0; i < (int)cverts.size(); ++i)
      for (int j = i + 1; j < (int)cverts.size(); ++j)
        diam = std::max(
            diam, Dist3(arr.verts[cverts[i]].pos, arr.verts[cverts[j]].pos));

    std::map<std::pair<int64_t, int64_t>, int> anchors;
    for (int nd : members) {
      const RN& rn = nodes[nd];
      const PSLGRegion& deg = faceRegions[rn.fi][rn.ri];
      for (const auto& r2 : faceRegions[rn.fi]) {
        if (!r2.classified || r2.degenerate) continue;
        bool adj = false;
        auto chk = [&](const std::vector<int>& lp) {
          const int n = (int)lp.size();
          for (int k = 0; k < n && !adj; ++k)
            if (RegionHasDirEdge(r2, lp[(k + 1) % n], lp[k])) adj = true;
        };
        chk(deg.loopVerts);
        for (const auto& h : deg.holeVerts) chk(h);
        if (adj) anchors[{r2.below, r2.above}]++;
      }
      for (int si : arr.faceSeams[rn.fi]) {
        const Seam& seam = arr.seams[si];
        const int ofi = (seam.faceId0 == rn.fi) ? seam.faceId1 : seam.faceId0;
        for (int k = 0; k + 1 < (int)seam.vertIds.size(); ++k) {
          const int sv = seam.vertIds[k], sw = seam.vertIds[k + 1];
          if (!RegionTouchesEdge(deg, sv, sw)) continue;
          for (const auto& r2 : faceRegions[ofi])
            if (r2.classified && !r2.degenerate &&
                RegionTouchesEdge(r2, sv, sw))
              anchors[{r2.below, r2.above}]++;
        }
      }
    }

    // Diameter-first ordering (spec M6 crucible):
    // - No anchors + large -> Unclassifiable (P9)
    // - No anchors + sub-eps -> drop as isolated eps feature
    // - Large + conflicting anchors -> AnchorConflict (P8)
    // - Large + single anchor -> Unclassifiable (span guard: can't trust anchor
    // alone)
    // - Sub-eps (any anchor count) -> classify from first anchor; conflicts
    // allowed
    //   because the region is truly tiny and won't affect emission materially
    //   (P10).
    if (anchors.empty()) {
      if (diam > eps) return FatalReason::UnclassifiableComponent;
      // Sub-eps with no anchors: isolated eps artifact, drop.
      for (int nd : members) {
        PSLGRegion& reg = faceRegions[nodes[nd].fi][nodes[nd].ri];
        reg.classified = true;
        reg.below = 0;
        reg.above = 0;
      }
      ++cnt.epsFeaturesDropped;
      continue;
    }
    if (diam > eps) {
      if (anchors.size() > 1) return FatalReason::AnchorConflict;
      return FatalReason::UnclassifiableComponent;
    }
    // Sub-eps with one or more anchors: take first anchor (map ordered by key).
    const int64_t agrB = anchors.begin()->first.first,
                  agrA = anchors.begin()->first.second;
    for (int nd : members) {
      PSLGRegion& reg = faceRegions[nodes[nd].fi][nodes[nd].ri];
      reg.classified = true;
      reg.below = agrB;
      reg.above = agrA;
    }
    ++cnt.degenerateClassified;
    if (IsInside3D(agrB) == IsInside3D(agrA)) ++cnt.epsFeaturesDropped;
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Stage E: emission and Impl construction
// ---------------------------------------------------------------------------

struct EmittedTri {
  ivec3 verts;
  vec3 normal;
};

static std::vector<EmittedTri> EmitRegion(
    const PSLGRegion& region, int64_t below, int64_t above,
    const CanonicalFace& face, const std::vector<MergedVert>& verts) {
  const int64_t m_lex = above - below;
  // m_lex==0 is impossible: the region is kept because
  // IsInside(below)!=IsInside(above).
  DEBUG_ASSERT(m_lex != 0, logicErr,
               "EmitRegion: m_lex==0 (spec: impossible for kept region)");
  const bool normInFaceDir =
      (m_lex > 0) ? IsInside3D(above) : IsInside3D(below);
  const vec3 outNormal = normInFaceDir ? face.normal : -face.normal;
  const mat2x3 projOut = GetAxisAlignedProjection(outNormal);

  PolygonsIdx polys;
  {
    SimplePolygonIdx outer;
    for (int v : region.loopVerts) outer.push_back({projOut * verts[v].pos, v});
    double area = 0.0;
    for (int i = 0; i < (int)outer.size(); ++i) {
      const vec2 &a = outer[i].pos, &b = outer[(i + 1) % outer.size()].pos;
      area += a.x * b.y - b.x * a.y;
    }
    if (area < 0) std::reverse(outer.begin(), outer.end());
    polys.push_back(std::move(outer));
  }
  for (const auto& hloop : region.holeVerts) {
    SimplePolygonIdx hole;
    for (int v : hloop) hole.push_back({projOut * verts[v].pos, v});
    double area = 0.0;
    for (int i = 0; i < (int)hole.size(); ++i) {
      const vec2 &a = hole[i].pos, &b = hole[(i + 1) % hole.size()].pos;
      area += a.x * b.y - b.x * a.y;
    }
    if (area > 0) std::reverse(hole.begin(), hole.end());
    polys.push_back(std::move(hole));
  }
  if (polys.empty() || polys[0].size() < 3) return {};
  const auto tris = TriangulateIdx(polys, -1.0);
  std::vector<EmittedTri> out;
  out.reserve(tris.size());
  for (const auto& t : tris) out.push_back({t, outNormal});
  return out;
}

static StageResult<Manifold::Impl> BuildImpl(
    const std::vector<EmittedTri>& emitted,
    const std::vector<MergedVert>& verts) {
  if (emitted.empty()) return StageResult<Manifold::Impl>::Ok(Manifold::Impl{});
  std::vector<int> usedV;
  usedV.reserve(emitted.size() * 3);
  for (const auto& t : emitted) {
    usedV.push_back(t.verts.x);
    usedV.push_back(t.verts.y);
    usedV.push_back(t.verts.z);
  }
  std::sort(usedV.begin(), usedV.end());
  usedV.erase(std::unique(usedV.begin(), usedV.end()), usedV.end());
  std::map<int, int> remap;
  for (int i = 0; i < (int)usedV.size(); ++i) remap[usedV[i]] = i;

  Manifold::Impl impl;
  impl.vertPos_.resize(usedV.size());
  for (int i = 0; i < (int)usedV.size(); ++i)
    impl.vertPos_[i] = verts[usedV[i]].pos;

  Vec<ivec3> tv(emitted.size());
  for (int i = 0; i < (int)emitted.size(); ++i)
    tv[i] = {remap[emitted[i].verts.x], remap[emitted[i].verts.y],
             remap[emitted[i].verts.z]};

  impl.CreateHalfedges(tv);
  if (!impl.IsManifold())
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::NonManifoldEmission,
        "emitted triangulation not 2-manifold");
  impl.InitializeOriginal();
  impl.CalculateBBox();
  impl.SetEpsilon();
  impl.SortGeometry();
  impl.SetNormalsAndCoplanar();
  return StageResult<Manifold::Impl>::Ok(std::move(impl));
}

}  // namespace

// ---------------------------------------------------------------------------
// Pipeline: stages C + D + E
// ---------------------------------------------------------------------------

static Overlap3Result RunStageCDE(ArrangementGeometry& arr, double eps,
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
  const int nFaces = (int)arr.faces.size();
  std::vector<std::vector<PSLGRegion>> faceRegions(nFaces);

  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];

    if (auto pf = ValidateFacePSLG(fi, arr)) {
      result.fatal = *pf;
      result.detail = "PSLG invalid on face " + std::to_string(fi);
      result.counters = cnt;
      return result;
    }

    faceRegions[fi] = BuildFacePSLG(fi, arr, eps);

    const double xv0 = arr.verts[face.verts.x].pos.x;
    const double xv1 = arr.verts[face.verts.y].pos.x;
    const double xv2 = arr.verts[face.verts.z].pos.x;
    const bool axisParallel =
        (std::abs(xv0 - xv1) < eps && std::abs(xv1 - xv2) < eps);

    if (axisParallel) {
      for (auto& reg : faceRegions[fi]) {
        auto cls = ClassifyAxisParallel(reg, face, slabs, arr.verts, arr.faces,
                                        xv0, eps);
        if (cls.has_value()) {
          reg.classified = true;
          reg.below = cls->first;
          reg.above = cls->second;
        } else
          reg.degenerate = true;
      }
      continue;
    }

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
      }
      if (cls.has_value()) {
        reg.classified = true;
        reg.below = cls->first;
        reg.above = cls->second;
      } else
        reg.degenerate = true;
    }
  }

  if (auto m6f = PropagateAnchorComponentsImpl(arr, faceRegions, eps, cnt)) {
    result.fatal = *m6f;
    result.detail = (*m6f == FatalReason::AnchorConflict)
                        ? "degenerate component: conflicting anchors"
                        : "degenerate component: unclassifiable";
    result.counters = cnt;
    return result;
  }

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

  std::vector<EmittedTri> emitted;
  for (int fi = 0; fi < nFaces; ++fi) {
    const CanonicalFace& face = arr.faces[fi];
    for (const auto& reg : faceRegions[fi]) {
      if (!reg.classified) continue;
      if (!IsInside3D(reg.below) && !IsInside3D(reg.above)) continue;
      if (IsInside3D(reg.below) == IsInside3D(reg.above)) continue;
      auto tris = EmitRegion(reg, reg.below, reg.above, face, arr.verts);
      for (auto& t : tris) emitted.push_back(std::move(t));
    }
  }

  auto buildRes = BuildImpl(emitted, arr.verts);
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
// Public entry points
// ---------------------------------------------------------------------------

Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps) {
  Overlap3Result result;
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    result.fatal = FatalReason::Starvation;
    result.detail = "epsilon not computable";
    return result;
  }
  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) {
    result.impl = Manifold::Impl{};
    return result;
  }
  Overlap3Counters& cnt = result.counters;
  StageBResult stageB = StageB(stageA, eps, cnt);
  if (stageB.fatal.has_value()) {
    result.fatal = stageB.fatal;
    result.detail = stageB.detail;
    return result;
  }
  return RunStageCDE(stageB.arr, eps, cnt);
}

Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps) {
  Overlap3Internals out;
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    out.fatal = FatalReason::Starvation;
    out.detail = "epsilon not computable";
    return out;
  }
  const StageAResult stageA = StageA(in, eps);
  if (stageA.faces.empty()) return out;
  StageBResult stageB = StageB(stageA, eps, out.counters);
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

// ---------------------------------------------------------------------------
// White-box test wrappers
// ---------------------------------------------------------------------------

Overlap3Result RemoveOverlaps3D_FromArr(const ArrangementGeometry& arr_in,
                                        double eps) {
  Overlap3Counters cnt;
  ArrangementGeometry arr_copy = arr_in;
  return RunStageCDE(arr_copy, eps, cnt);
}

std::optional<FatalReason> CheckSeamBalance_Test(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions) {
  std::optional<FatalReason> f;
  std::string d;
  CheckSeamBalance(arr, faceRegions, f, d);
  return f;
}

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
    out.below = cls->first;
    out.above = cls->second;
  }
  return out;
}

std::optional<FatalReason> PropagateAnchorComponents_Test(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt) {
  return PropagateAnchorComponentsImpl(arr, faceRegions, eps, cnt);
}

}  // namespace manifold
