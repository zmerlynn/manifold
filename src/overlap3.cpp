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

// The 3D sweep-native emission prototype: canonicalize -> seams -> slabs ->
// caps -> strips -> assembly.  Slab construction lives in overlap3_sweep.cpp.
// Design: docs/SweepEmit3D.md (three crucible rounds).

#include "overlap3.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "boolean2.h"
#include "disjoint_sets.h"
#include "impl.h"
#include "manifold/optional_assert.h"
#include "polygon_internal.h"
#include "shared.h"

namespace manifold {

// Slabs stage (overlap3_sweep.cpp) forward declaration.
StageResult<std::vector<SlabResult>> BuildSlabs(const ArrangementGeometry& arr,
                                                double eps,
                                                Overlap3Counters& cnt);

namespace {

// ---------------------------------------------------------------------------
// Seams-stage helpers: seam segments and in-plane crossing computation.
// ---------------------------------------------------------------------------

// Clip line P + t*D against face triangle; returns the parameter interval,
// or nullopt when the overlap is empty.
std::optional<std::pair<double, double>> LineTriClip(vec3 P, vec3 D, vec3 v0,
                                                     vec3 v1, vec3 v2, vec3 n,
                                                     double eps) {
  double tLo = -std::numeric_limits<double>::infinity();
  double tHi = std::numeric_limits<double>::infinity();
  const vec3 edges[3] = {v1 - v0, v2 - v1, v0 - v2};
  const vec3 vBase[3] = {v0, v1, v2};
  const double kTol = eps * la::length(n);
  for (const int i : {0, 1, 2}) {
    const vec3 ev = la::cross(edges[i], D),
               evB = la::cross(edges[i], vBase[i] - P);
    const double den = la::dot(n, ev), num = la::dot(n, evB);
    if (std::abs(den) <= kTol) {
      if (num > kTol) return std::nullopt;
    } else {
      const double t = num / den;
      if (den > 0)
        tLo = std::max(tLo, t);
      else
        tHi = std::min(tHi, t);
    }
  }
  if (tLo > tHi + eps) return std::nullopt;
  return std::make_pair(tLo, tHi);
}

// Compute the seam segment [first, second] as the intersection of face
// triangles (a0,a1,a2) with normal na and (b0,b1,b2) with normal nb.
// nullopt only when the triangles are disjoint (or parallel-plane); a
// point/tangent contact admitted by the eps clip model returns an equal pair,
// so the caller's degenerate-contact arm records its x-critical.
std::optional<std::pair<vec3, vec3>> TriTriSeam(vec3 a0, vec3 a1, vec3 a2,
                                                vec3 na, vec3 b0, vec3 b1,
                                                vec3 b2, vec3 nb, double eps) {
  const vec3 D = la::cross(na, nb);
  const double Dlen = la::length(D);
  if (Dlen == 0.0) return std::nullopt;
  const vec3 dir = D / Dlen;
  const double da = la::dot(na, a0), db = la::dot(nb, b0);
  const vec3 absD = la::abs(dir);
  const int mc = (absD.x >= absD.y && absD.x >= absD.z) ? 0
                 : (absD.y >= absD.z)                   ? 1
                                                        : 2;
  const int c1 = (mc + 1) % 3, c2 = (mc + 2) % 3;
  const double a11 = na[c1], a12 = na[c2], b11 = nb[c1], b12 = nb[c2];
  const double det = a11 * b12 - a12 * b11;
  if (det == 0.0) return std::nullopt;
  vec3 P(0.0);
  P[c1] = (da * b12 - db * a12) / det;
  P[c2] = (a11 * db - b11 * da) / det;
  const auto clipA = LineTriClip(P, dir, a0, a1, a2, na, eps);
  if (!clipA) return std::nullopt;
  const auto clipB = LineTriClip(P, dir, b0, b1, b2, nb, eps);
  if (!clipB) return std::nullopt;
  const double tLo = std::max(clipA->first, clipB->first);
  const double tHi = std::min(clipA->second, clipB->second);
  // dir is unit length, so the parameter interval is in length units and
  // eps applies directly (matching LineTriClip's own slop).
  if (tLo > tHi + eps) return std::nullopt;  // genuinely disjoint
  if (tLo >= tHi) {
    // Point/tangent contact: collapse to the interval midpoint.
    const vec3 q = P + (0.5 * (tLo + tHi)) * dir;
    return std::make_pair(q, q);
  }
  return std::make_pair(P + tLo * dir, P + tHi * dir);
}

// The crossing x-coordinate of two coplanar line segments [A0,A1] and
// [B0,B1] (both in a shared face plane), or nullopt unless the lines are
// non-parallel and the crossing is strictly interior to both seams.  The x
// is the only output per spec SEAMS (over-inclusion harmless).
std::optional<double> SeamSeamCrossX(vec3 A0, vec3 A1, vec3 B0, vec3 B1,
                                     double eps) {
  const vec3 dA = A1 - A0, dB = B1 - B0, dC = B0 - A0;
  const double lenA = la::length(dA), lenB = la::length(dB);
  if (lenA < eps || lenB < eps) return std::nullopt;
  const vec3 cAB = la::cross(dA, dB);
  const double cABlen2 = la::length2(cAB);
  // Near-parallel gate, dimensionally correct: |cross| has units len^2, so
  // compare against eps * (lenA + lenB).  Endpoint noise on near-parallel
  // seams otherwise yields a pseudo-crossing at a meaningless x.
  const double parTol = eps * (lenA + lenB);
  if (cABlen2 <= parTol * parTol) return std::nullopt;
  // t on A: t * |cAB|^2 = dot(cross(dC, dB), cAB)
  const double t = la::dot(la::cross(dC, dB), cAB) / cABlen2;
  // s on B: s * |cAB|^2 = dot(cross(dC, dA), cAB)
  const double s = la::dot(la::cross(dC, dA), cAB) / cABlen2;
  // Require strictly interior to both seams (not at endpoints).
  const double tEps = eps / lenA, sEps = eps / lenB;
  if (t <= tEps || t >= 1.0 - tEps) return std::nullopt;
  if (s <= sEps || s >= 1.0 - sEps) return std::nullopt;
  return A0.x + t * dA.x;
}

// Find or add a vert within eps of pos; return its index.
int FindOrAddVert(std::vector<MergedVert>& verts, vec3 pos, double eps) {
  for (int i = 0; i < static_cast<int>(verts.size()); ++i)
    if (la::length(verts[i].pos - pos) <= eps) return i;
  const int id = static_cast<int>(verts.size());
  verts.push_back({pos});
  return id;
}

// ---------------------------------------------------------------------------
// Canonicalize stage: vert merge, face multiplicity accumulation.
// ---------------------------------------------------------------------------

struct CanonicalGeometry {
  std::vector<vec3> mergedVerts;
  std::vector<int> origVert2Merged;
  std::vector<CanonicalFace> faces;
};

CanonicalGeometry Canonicalize(const Manifold::Impl& in, double eps) {
  const int nVerts = static_cast<int>(in.vertPos_.size());
  const int nTris = static_cast<int>(in.halfedge_.size()) / 3;

  DisjointSets uf(nVerts);
  for (int i = 0; i < nVerts; ++i)
    for (int j = i + 1; j < nVerts; ++j)
      if (la::length(in.vertPos_[i] - in.vertPos_[j]) <= eps) uf.unite(i, j);

  std::map<int, std::vector<int>> comps;
  for (int i = 0; i < nVerts; ++i)
    comps[static_cast<int>(uf.find(i))].push_back(i);

  std::vector<int> origVert2Merged(nVerts);
  std::vector<vec3> mergedVerts;
  std::map<int, int> rootToMerged;
  for (auto& [root, members] : comps) {
    vec3 c(0.0);
    for (int v : members) c += in.vertPos_[v];
    c /= static_cast<double>(members.size());
    int best = members[0];
    double bestD = la::length(in.vertPos_[best] - c);
    for (int v : members) {
      double d = la::length(in.vertPos_[v] - c);
      if (d < bestD) {
        bestD = d;
        best = v;
      }
    }
    rootToMerged[root] = static_cast<int>(mergedVerts.size());
    mergedVerts.push_back(in.vertPos_[best]);
  }
  for (int i = 0; i < nVerts; ++i)
    origVert2Merged[i] = rootToMerged[static_cast<int>(uf.find(i))];

  // Sorted vert triple; std::array's lexicographic order keys the map.
  using FaceKey = std::array<int, 3>;
  struct FaceEntry {
    int64_t mult;
    int repTri;
    int repPar;
  };
  std::map<FaceKey, FaceEntry> faceMap;

  for (int tri = 0; tri < nTris; ++tri) {
    const int h = 3 * tri;
    int v0 = origVert2Merged[in.halfedge_.Start(h)],
        v1 = origVert2Merged[in.halfedge_.Start(h + 1)],
        v2 = origVert2Merged[in.halfedge_.Start(h + 2)];
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
    const FaceKey key{a, b, c};
    auto it = faceMap.find(key);
    if (it == faceMap.end())
      faceMap.emplace(key, FaceEntry{1, tri, par});
    else
      it->second.mult += (par == it->second.repPar) ? 1 : -1;
  }

  CanonicalGeometry res;
  res.mergedVerts = std::move(mergedVerts);
  res.origVert2Merged = std::move(origVert2Merged);
  for (auto& [key, e] : faceMap) {
    if (e.mult == 0) continue;
    const int h = 3 * e.repTri;
    const int mv0 = res.origVert2Merged[in.halfedge_.Start(h)];
    const int mv1 = res.origVert2Merged[in.halfedge_.Start(h + 1)];
    const int mv2 = res.origVert2Merged[in.halfedge_.Start(h + 2)];
    const vec3 n = la::cross(res.mergedVerts[mv1] - res.mergedVerts[mv0],
                             res.mergedVerts[mv2] - res.mergedVerts[mv0]);
    if (la::length2(n) == 0.0) continue;
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
// Seams stage: coplanar grouping, seam segments, vertex-free extra criticals
// (seam-seam and in-plane skeleton crossings).  No triple-point unification,
// no edge subdivision.
// ---------------------------------------------------------------------------

struct SeamsResult {
  ArrangementGeometry arr;
  std::optional<FatalReason> fatal;
  std::string detail;
};

SeamsResult FindSeams(const CanonicalGeometry& canon, double eps,
                      Overlap3Counters& cnt) {
  SeamsResult res;
  ArrangementGeometry& arr = res.arr;
  const int nFaces = static_cast<int>(canon.faces.size());

  arr.verts.resize(canon.mergedVerts.size());
  for (int i = 0; i < static_cast<int>(canon.mergedVerts.size()); ++i)
    arr.verts[i] = {canon.mergedVerts[i]};
  arr.faces = canon.faces;

  std::vector<Box> boxes(nFaces);
  std::vector<double> planeDist(nFaces);
  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& f = canon.faces[fi];
    planeDist[fi] = la::dot(f.normal, arr.verts[f.verts.x].pos);
    const vec3 p0 = arr.verts[f.verts.x].pos, p1 = arr.verts[f.verts.y].pos,
               p2 = arr.verts[f.verts.z].pos;
    boxes[fi] = Box(la::min(la::min(p0, p1), p2) - vec3(eps),
                    la::max(la::max(p0, p1), p2) + vec3(eps));
  }

  // Coplanar grouping pre-pass (spec COPLANAR mechanism 1): union faces
  // whose planes coincide within eps - the symmetric OR of the two direction
  // tests - INCLUDING shared-edge pairs (a folded flap's pair must join
  // its group before any exemption skips it).  Grouping is by PLANE, not by
  // overlap: non-overlapping same-plane faces grouping is benign (their
  // sections never coincide).  Components of two or more faces become
  // groups; their in-plane content resolves through group-id seeding and
  // the caps/strips machinery instead of failing closed.
  {
    DisjointSets faceUf(nFaces);
    auto inPlaneOf = [&](int host, const CanonicalFace& Q) {
      const vec3& n = canon.faces[host].normal;
      const double d = planeDist[host];
      return std::abs(la::dot(n, arr.verts[Q.verts.x].pos) - d) <= eps &&
             std::abs(la::dot(n, arr.verts[Q.verts.y].pos) - d) <= eps &&
             std::abs(la::dot(n, arr.verts[Q.verts.z].pos) - d) <= eps;
    };
    for (int fi = 0; fi < nFaces; ++fi)
      for (int fj = fi + 1; fj < nFaces; ++fj)
        if (inPlaneOf(fi, canon.faces[fj]) || inPlaneOf(fj, canon.faces[fi]))
          faceUf.unite(fi, fj);
    arr.face2Group.assign(nFaces, -1);
    std::map<int, int> rootCount;
    for (int fi = 0; fi < nFaces; ++fi)
      ++rootCount[static_cast<int>(faceUf.find(fi))];
    std::map<int, int> root2Group;
    for (const auto& [root, count] : rootCount)
      if (count >= 2) {
        const int g = arr.numGroups++;
        root2Group[root] = g;
      }
    for (int fi = 0; fi < nFaces; ++fi) {
      const auto it = root2Group.find(static_cast<int>(faceUf.find(fi)));
      if (it != root2Group.end()) arr.face2Group[fi] = it->second;
    }
  }

  for (int fi = 0; fi < nFaces; ++fi) {
    for (int fj = fi + 1; fj < nFaces; ++fj) {
      if (!boxes[fi].DoesOverlap(boxes[fj])) continue;

      // Same plane group: no seam exists (identical planes), and the pair's
      // in-plane interaction is the group machinery's job (spec COPLANAR).
      if (arr.face2Group[fi] >= 0 && arr.face2Group[fi] == arr.face2Group[fj])
        continue;

      const CanonicalFace &FA = canon.faces[fi], &FB = canon.faces[fj];
      {
        // Faces sharing a merged edge produce no useful seam: their planes
        // meet along the shared edge, whose verts are already criticals.
        // Coplanar shared-edge pairs (touching diagonal triangulations,
        // folded flaps) were unioned by the grouping pre-pass above and
        // resolve through the group machinery (spec COPLANAR
        // re-adjudication: a same-winding flap is ordinary +2 content).
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

      // Parallel distinct planes: no seam.  (Same-plane pairs were grouped
      // above; a pair reaching here with parallel planes is separated by
      // more than eps - ordinary disjoint-parallel geometry.)
      const vec3 Dcross = la::cross(na, nb);
      if (la::length(Dcross) == 0.0) continue;

      // Seam computation.
      const auto seam = TriTriSeam(pa0, pa1, pa2, na, pb0, pb1, pb2, nb, eps);
      if (!seam) continue;
      const vec3 qA = seam->first, qB = seam->second;
      const double seamLen = la::length(qB - qA);

      if (seamLen <= eps) {
        // Degenerate contact: record the endpoint x's as criticals (no vert
        // identity needed). The slabs stage's SubEpsFeature guard will catch
        // any macro-scale hazard.
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

  // M1: seam-seam triple criticals (spec SEAMS). For each pair of seams sharing
  // a face, compute the crossing of their in-plane projections. The crossing x
  // is the only quantity consumed (over-inclusion is harmless).
  const int nSeams = static_cast<int>(arr.seams.size());
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
      const auto xCross = SeamSeamCrossX(A0, A1, B0, B1, eps);
      if (xCross) arr.criticalXs.push_back(*xCross);
    }
  }

  // In-plane skeleton criticals (spec COPLANAR mechanism 4).  Each group's
  // skeleton = its members' edges (all three per face - internal shared
  // edges included, over-inclusion is harmless) plus the in-plane seams of
  // transversal faces with members (such a seam lies in the shared plane by
  // construction).  Pairwise skeleton crossings are section-combinatorics
  // events; SeamSeamCrossX's contract (coplanar 3D segments) computes them.
  if (arr.numGroups > 0) {
    std::vector<std::vector<std::pair<vec3, vec3>>> skeleton(arr.numGroups);
    std::set<std::tuple<int, int, int>> seenEdges;  // (group, loVert, hiVert)
    for (int fi = 0; fi < nFaces; ++fi) {
      const int g = arr.face2Group[fi];
      if (g < 0) continue;
      const int vs[3] = {canon.faces[fi].verts.x, canon.faces[fi].verts.y,
                         canon.faces[fi].verts.z};
      for (const int k : {0, 1, 2}) {
        int a = vs[k], b = vs[(k + 1) % 3];
        if (a > b) std::swap(a, b);
        if (!seenEdges.insert({g, a, b}).second) continue;
        skeleton[g].push_back({arr.verts[a].pos, arr.verts[b].pos});
      }
    }
    for (const auto& seam : arr.seams) {
      const std::pair<vec3, vec3> seg = {arr.verts[seam.vertId0].pos,
                                         arr.verts[seam.vertId1].pos};
      const int g0 = arr.face2Group[seam.faceId0];
      const int g1 = arr.face2Group[seam.faceId1];
      if (g0 >= 0) skeleton[g0].push_back(seg);
      if (g1 >= 0 && g1 != g0) skeleton[g1].push_back(seg);
    }
    for (const auto& segs : skeleton) {
      const int n = static_cast<int>(segs.size());
      for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
          const auto xCross =
              SeamSeamCrossX(segs[i].first, segs[i].second, segs[j].first,
                             segs[j].second, eps);
          if (xCross) arr.criticalXs.push_back(*xCross);
        }
      }
    }
  }

  return res;
}

// ---------------------------------------------------------------------------
// Strips stage: strip emission.
// Track extension evaluates a piece endpoint's (y,z) at a boundary critical;
// the extended limits feed the cap arrangements (caps stage), whose output
// chains the strips then zip (spec [R2-fold]: strips consume the cap
// subdivision, they do not re-derive it).
// ---------------------------------------------------------------------------

// Per-slab extension resolver (spec COPLANAR mechanism 3): each section
// vertex resolves ONCE to a 3D track (or a weld), and every incident piece
// endpoint extends through that same resolution - same-side closure holds by
// construction.  Candidate classes in priority order:
//   (i)   face cutting-edge tracks, matched at their segment ends;
//   (ii)  seam tracks (class-ii crossings);
//   (iii) in-plane edges of coplanar-group members crossing this slab;
//   (iv)  weld - constant (y,z), the forced-through fallback; its deviation
//         is cap-covered (spec STRIPS [R2-fold]).
// Ties resolve to the lowest class, then track order (face id / seam order /
// edge order) - deterministic.  A track evaluated at a critical where it
// terminates lands bitwise on the 3D vertex (InterpolateSafe at the
// endpoint's own x): the spec's 3D-identity preference falls out of class
// priority (terminating tracks beat welds).
class SlabResolver {
 public:
  SlabResolver() = default;
  SlabResolver(const SlabResult& slab, const ArrangementGeometry& arr,
               double eps)
      : eps_(eps) {
    if (!slab.built) return;
    // Class-iii candidates: edges of grouped faces crossing this slab's
    // section plane (all three edges per member; over-inclusion harmless).
    if (arr.numGroups > 0) {
      std::set<std::pair<int, int>> seen;
      for (int fi = 0; fi < static_cast<int>(arr.faces.size()); ++fi) {
        if (arr.face2Group[fi] < 0) continue;
        const int vs[3] = {arr.faces[fi].verts.x, arr.faces[fi].verts.y,
                           arr.faces[fi].verts.z};
        for (const int k : {0, 1, 2}) {
          int a = vs[k], b = vs[(k + 1) % 3];
          if (a > b) std::swap(a, b);
          if (!seen.insert({a, b}).second) continue;
          const vec3 pa = arr.verts[a].pos, pb = arr.verts[b].pos;
          if (pa.x == pb.x) continue;  // no crossing trajectory
          if (std::min(pa.x, pb.x) > slab.xMid ||
              std::max(pa.x, pb.x) < slab.xMid)
            continue;
          planar_.push_back({pa, pb});
        }
      }
    }
    // Resolve every piece endpoint once.
    for (const SweepCapture& piece : slab.pieces) {
      for (const vec2 pt : {piece.from, piece.to}) {
        const auto key = std::make_pair(pt.x, pt.y);
        if (tracks_.count(key)) continue;
        tracks_.emplace(key, Resolve(pt, slab));
      }
    }
  }

  // Extend a piece endpoint (exact section-vert coordinates) to xTarget.
  vec2 Extend(vec2 pt, double xTarget) const {
    const auto it = tracks_.find({pt.x, pt.y});
    DEBUG_ASSERT(it != tracks_.end(), logicErr,
                 "unresolved section vertex in SlabResolver");
    if (it == tracks_.end() || it->second.weld) return pt;
    return InterpolateSafe(it->second.a, it->second.b, xTarget);
  }

 private:
  struct Track {
    bool weld;
    vec3 a, b;
  };

  Track Resolve(vec2 pt, const SlabResult& slab) const {
    for (const FaceTrack& ft : slab.faceTracks) {
      if (la::length(pt - ft.p0) <= eps_) return {false, ft.va0, ft.vb0};
      if (la::length(pt - ft.p1) <= eps_) return {false, ft.va1, ft.vb1};
    }
    for (const SeamTrackEntry& st : slab.seamTracks) {
      if (la::length(pt - st.yzMid) <= eps_) return {false, st.vA, st.vB};
    }
    for (const auto& e : planar_) {
      const vec2 q = InterpolateSafe(e.first, e.second, slab.xMid);
      if (la::length(pt - q) <= eps_) return {false, e.first, e.second};
    }
    return {true, vec3(0.0), vec3(0.0)};
  }

  double eps_ = 0.0;
  std::vector<std::pair<vec3, vec3>> planar_;
  std::map<std::pair<double, double>, Track> tracks_;
};

struct OutTri3D {
  vec3 v[3];
};

// The third consumer of the cap arrangement (spec CAPS [R2-fold]): its output
// verts subdivide the adjacent strip edges.  `arrVerts` is the arrangement's
// FULL vert set (merged inputs + every collected-arrangement vertex, per the
// engine's negEdges contract).  Returns the verts lying strictly interior to
// the chord [p0, p1] (within eps - the tolerance-model on-edge test, same
// posture as the engine's own incidence rule), sorted by projection parameter
// (exact ties by lex order).  The returned POSITIONS become the strip edge's
// chain verts, bitwise equal to the cap triangulation corners they pair with.
std::vector<vec2> ChainSplitVerts(const std::vector<vec2>& arrVerts, vec2 p0,
                                  vec2 p1, double eps) {
  const vec2 segD = p1 - p0;
  const double segLen2 = la::length2(segD);
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
void ZipperEmit(std::vector<OutTri3D>& out, double xLo, double xHi,
                const std::vector<vec2>& a, const std::vector<vec2>& b) {
  const int m = static_cast<int>(a.size()), n = static_cast<int>(b.size());
  auto params = [](const std::vector<vec2>& c) {
    std::vector<double> t(c.size(), 0.0);
    const vec2 d = c.back() - c.front();
    const double len2 = la::length2(d);
    if (c.size() >= 2 && len2 > 0.0) {
      for (int i = 1; i + 1 < static_cast<int>(c.size()); ++i)
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

// Strips stage: each strip zips its two c-side chains, built by the cap
// arrangements at its slab's boundary-pair canonical criticals (spec CAPS
// [R2-fold] third consumer).  No geometry is computed here: the chains carry
// the cap arrangements' vert positions bitwise.  Every built slab has both
// sides bound (its own xHi is always its pair's canonical; its lo is bound at
// the preceding gap's canonical), with one chain per piece - attribution
// failures fail closed in EmitCaps before this runs.
std::optional<std::pair<FatalReason, std::string>> EmitStrips(
    std::vector<OutTri3D>& out, const std::vector<SlabResult>& slabs,
    const std::vector<StripChains>& chains) {
  for (int si = 0; si < static_cast<int>(slabs.size()); ++si) {
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
      ZipperEmit(out, slabs[si].xLo, slabs[si].xHi, ch.lo[k], ch.hi[k]);
    }
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Caps stage: cap emission.
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
bool TriangulateCap(std::vector<OutTri3D>& out, const Polygons& polys,
                    double xCap, bool flipWinding, double eps) {
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

  const int nAll = static_cast<int>(allVerts.size());
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
// retained pieces to x=c through their SlabResolvers (spec COPLANAR
// mechanism 3), so every piece endpoint at one section vertex extends
// identically.  Endpoint slots align 1:1 with slab.pieces.
// NO vertex pre-merge: the cap arrangement's MergeVerts owns snapping.
struct CapEdgeSet {
  std::vector<vec2> rawVerts;
  std::vector<std::pair<int, int>> lPiece2RawVerts,
      rPiece2RawVerts;  // rawVerts index pairs
};

CapEdgeSet BuildCapEdgeSet(const SlabResult* leftSlab,
                           const SlabResolver* leftResolver,
                           const SlabResult* rightSlab,
                           const SlabResolver* rightResolver, double xCap) {
  CapEdgeSet ces;
  auto addPieces = [&](const SlabResult& slab, const SlabResolver& resolver,
                       std::vector<std::pair<int, int>>& slots) {
    slots.reserve(slab.pieces.size());
    for (const auto& piece : slab.pieces) {
      // Piece is emission-oriented (interior-on-left) per spec engine contract.
      const vec2 eFrom = resolver.Extend(piece.from, xCap);
      const vec2 eTo = resolver.Extend(piece.to, xCap);
      const int i0 = static_cast<int>(ces.rawVerts.size());
      ces.rawVerts.push_back(eFrom);
      ces.rawVerts.push_back(eTo);
      slots.push_back({i0, i0 + 1});
    }
  };
  // L verts appended first, then R: one deterministic input order for the ONE
  // arrangement.
  if (leftSlab && leftSlab->built)
    addPieces(*leftSlab, *leftResolver, ces.lPiece2RawVerts);
  if (rightSlab && rightSlab->built)
    addPieces(*rightSlab, *rightResolver, ces.rPiece2RawVerts);
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
std::optional<std::pair<FatalReason, std::string>> ComputeCap(
    std::vector<OutTri3D>& out, std::vector<std::vector<vec2>>* leftChains,
    std::vector<std::vector<vec2>>* rightChains, Overlap3Counters& cnt,
    const SlabResult* leftSlab, const SlabResolver* leftResolver,
    const SlabResult* rightSlab, const SlabResolver* rightResolver, double xCap,
    double eps) {
  const CapEdgeSet ces =
      BuildCapEdgeSet(leftSlab, leftResolver, rightSlab, rightResolver, xCap);
  // Pre-size bound chain outputs so every piece has a slot.
  if (leftChains) leftChains->assign(ces.lPiece2RawVerts.size(), {});
  if (rightChains) rightChains->assign(ces.rPiece2RawVerts.size(), {});
  if (ces.rawVerts.empty()) return std::nullopt;

  std::vector<EdgeM> edges;
  edges.reserve(ces.lPiece2RawVerts.size() + ces.rPiece2RawVerts.size());
  for (const auto& s : ces.lPiece2RawVerts)
    edges.push_back({s.first, s.second, +1});
  for (const auto& s : ces.rPiece2RawVerts)
    edges.push_back({s.first, s.second, -1});

  std::vector<OutEdge> negEdges;
  ++cnt.capArrangements;
  const OverlapResult r =
      RemoveOverlaps2D(ces.rawVerts, edges, eps, /*debug=*/false, WindRule::Add,
                       /*trace=*/nullptr, &negEdges);
  bool loopsClosed = true;
  bool trisOk = true;
  if (!r.edges.empty()) {
    const Polygons cp = OutEdgesToPolygons(r.verts, r.edges, &loopsClosed);
    trisOk = TriangulateCap(out, cp, xCap, /*flipWinding=*/false, eps);
  }
  if (loopsClosed && trisOk && !negEdges.empty()) {
    const Polygons cm = OutEdgesToPolygons(r.verts, negEdges, &loopsClosed);
    trisOk = TriangulateCap(out, cm, xCap, /*flipWinding=*/true, eps);
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
  if (leftChains) buildChains(ces.lPiece2RawVerts, *leftChains);
  if (rightChains) buildChains(ces.rPiece2RawVerts, *rightChains);
  return std::nullopt;
}

// Caps-stage driver: one cap per critical (spec [R3-fold]: runs are never
// merged; only IEEE-exact duplicate criticals collapse).  Adjacent built
// slabs are found by INDEX: slab bounds are these exact critical values by
// construction, so no slack enters the lookup.
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
std::optional<std::pair<FatalReason, std::string>> EmitCaps(
    std::vector<OutTri3D>& out, std::vector<StripChains>& chains,
    Overlap3Counters& cnt, const std::vector<SlabResult>& slabs,
    const std::vector<SlabResolver>& resolvers,
    const std::vector<double>& crits, double eps) {
  const int nSlabs = static_cast<int>(slabs.size());
  for (int ci = 0; ci < static_cast<int>(crits.size()); ++ci) {
    int li = ci - 1;
    while (li >= 0 && !slabs[li].built) --li;
    int ri = ci;
    while (ri < nSlabs && !slabs[ri].built) ++ri;
    const SlabResult* left = li >= 0 ? &slabs[li] : nullptr;
    const SlabResult* right = ri < nSlabs ? &slabs[ri] : nullptr;
    const bool canonical = ci == li + 1;
    if (auto fatal =
            ComputeCap(out, (canonical && left) ? &chains[li].hi : nullptr,
                       (canonical && right) ? &chains[ri].lo : nullptr, cnt,
                       left, li >= 0 ? &resolvers[li] : nullptr, right,
                       ri < nSlabs ? &resolvers[ri] : nullptr, crits[ci], eps))
      return fatal;
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Assembly: weld verts, drop exact-duplicate triangles, build Impl.
// ---------------------------------------------------------------------------

StageResult<Manifold::Impl> BuildImpl(const std::vector<OutTri3D>& tris,
                                      double eps) {
  if (tris.empty()) return StageResult<Manifold::Impl>::Ok(Manifold::Impl{});

  // Collect verts with an eps-weld.
  std::vector<vec3> verts;
  auto getVertIdx = [&](vec3 p) -> int {
    for (int i = 0; i < static_cast<int>(verts.size()); ++i)
      if (la::length(verts[i] - p) <= eps) return i;
    int id = static_cast<int>(verts.size());
    verts.push_back(p);
    return id;
  };

  // Filter degenerate and duplicate triangles.
  // Degenerate: strip quads whose corners collapse within eps produce v0==v1
  // etc., which would crash CreateHalfedges.
  // Duplicate: per-critical caps (spec [R3-fold]: runs are never merged) mean
  // a sub-eps critical pair computes the SAME macro difference twice - at
  // x=c1 and x=c2 with |c2-c1| <= eps - and both triangulations collapse to
  // identical vertex triples after the weld.  This drop suppresses only
  // post-weld identical triangles, never distinct caps (the [R3] sentence
  // "their differences are empty").
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
  for (int i = 0; i < static_cast<int>(verts.size()); ++i)
    impl.vertPos_[i] = verts[i];

  impl.CreateHalfedges(tv);
  // The full 2-manifold gate (edge pairing AND vertex links): edge-on-face
  // touching contact welds into a genuinely non-manifold union - the honest
  // outcome is this named fatal, not a downstream assertion.
  if (!impl.IsManifold() || !impl.Is2Manifold()) {
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

// ---------------------------------------------------------------------------
// Pipeline: slabs -> caps -> strips -> assembly.
// ---------------------------------------------------------------------------

Overlap3Result SweepEmit(ArrangementGeometry& arr, double eps,
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
  // recompute, no slack).
  std::vector<double> crits;
  if (!slabs.empty()) {
    crits.reserve(slabs.size() + 1);
    for (const auto& s : slabs) crits.push_back(s.xLo);
    crits.push_back(slabs.back().xHi);
  }

  // Caps before strips: each cap's ONE arrangement is the source of
  // truth at its critical - cap_plus, cap_minus, and both adjacent slabs'
  // strip chains all consume it (spec [R2-fold]).
  // Per-slab extension resolvers (spec COPLANAR mechanism 3), built once;
  // caps on both sides of a slab consume the same resolution table.
  std::vector<SlabResolver> resolvers;
  resolvers.reserve(slabs.size());
  for (const SlabResult& slab : slabs) resolvers.emplace_back(slab, arr, eps);

  std::vector<OutTri3D> emitted;
  std::vector<StripChains> chains(slabs.size());
  if (auto capFatal =
          EmitCaps(emitted, chains, cnt, slabs, resolvers, crits, eps)) {
    result.fatal = capFatal->first;
    result.detail = std::move(capFatal->second);
    result.counters = cnt;
    return result;
  }
  if (auto stripFatal = EmitStrips(emitted, slabs, chains)) {
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

}  // namespace

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
  const CanonicalGeometry canon = Canonicalize(in, eps);
  if (canon.faces.empty()) {
    result.impl = Manifold::Impl{};
    return result;
  }
  Overlap3Counters& cnt = result.counters;
  SeamsResult seams = FindSeams(canon, eps, cnt);
  if (seams.fatal.has_value()) {
    result.fatal = seams.fatal;
    result.detail = seams.detail;
    return result;
  }
  return SweepEmit(seams.arr, eps, cnt);
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
  const CanonicalGeometry canon = Canonicalize(in, eps);
  if (canon.faces.empty()) return out;
  SeamsResult seams = FindSeams(canon, eps, out.counters);
  if (seams.fatal.has_value()) {
    out.fatal = seams.fatal;
    out.detail = seams.detail;
    return out;
  }
  out.arr = std::move(seams.arr);
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
