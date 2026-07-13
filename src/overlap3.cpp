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
#include <unordered_map>
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

// A fail-closed outcome for the emission stages (caps/strips): a reason plus
// detail, or nullopt on success.  These stages produce no product - only the
// possibility of failing closed - so they return MaybeFatal rather than a
// StageResult<T>.
struct Fatal {
  FatalReason reason;
  std::string detail;
};
using MaybeFatal = std::optional<Fatal>;

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

// Deduped 3D edges of coplanar-group members, tagged by group index - the
// shared source for the skeleton criticals and the resolver's class-iii
// candidates (one concept, one collection).
struct GroupedEdge {
  int group;
  vec3 a, b;
};

std::vector<GroupedEdge> CollectGroupedMemberEdges(
    const ArrangementGeometry& arr) {
  std::vector<GroupedEdge> out;
  if (arr.numGroups == 0) return out;
  std::set<std::tuple<int, int, int>> seen;  // (group, loVert, hiVert)
  for (int fi = 0; fi < static_cast<int>(arr.faces.size()); ++fi) {
    const int g = arr.face2Group[fi];
    if (g < 0) continue;
    const int vs[3] = {arr.faces[fi].verts.x, arr.faces[fi].verts.y,
                       arr.faces[fi].verts.z};
    for (const int k : {0, 1, 2}) {
      int a = vs[k], b = vs[(k + 1) % 3];
      if (a > b) std::swap(a, b);
      if (!seen.insert({g, a, b}).second) continue;
      out.push_back({g, arr.verts[a], arr.verts[b]});
    }
  }
  return out;
}

// Pairwise in-plane crossing x's of coplanar segments, appended as criticals.
void AppendPairwiseCrossXs(const std::vector<std::pair<vec3, vec3>>& segs,
                           double eps, std::vector<double>& criticalXs) {
  const int n = static_cast<int>(segs.size());
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const auto xCross = SeamSeamCrossX(segs[i].first, segs[i].second,
                                         segs[j].first, segs[j].second, eps);
      if (xCross) criticalXs.push_back(*xCross);
    }
  }
}

// Find or add a vert within eps of pos; return its index.
int FindOrAddVert(std::vector<vec3>& verts, vec3 pos, double eps) {
  for (int i = 0; i < static_cast<int>(verts.size()); ++i)
    if (la::length(verts[i] - pos) <= eps) return i;
  const int id = static_cast<int>(verts.size());
  verts.push_back(pos);
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
    arr.verts[i] = canon.mergedVerts[i];
  arr.faces = canon.faces;

  std::vector<Box> boxes(nFaces);
  std::vector<double> planeDist(nFaces);
  for (int fi = 0; fi < nFaces; ++fi) {
    const auto& f = canon.faces[fi];
    planeDist[fi] = la::dot(f.normal, arr.verts[f.verts.x]);
    const vec3 p0 = arr.verts[f.verts.x], p1 = arr.verts[f.verts.y],
               p2 = arr.verts[f.verts.z];
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
      return std::abs(la::dot(n, arr.verts[Q.verts.x]) - d) <= eps &&
             std::abs(la::dot(n, arr.verts[Q.verts.y]) - d) <= eps &&
             std::abs(la::dot(n, arr.verts[Q.verts.z]) - d) <= eps;
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
      const vec3 pa0 = arr.verts[FA.verts.x], pa1 = arr.verts[FA.verts.y],
                 pa2 = arr.verts[FA.verts.z];
      const vec3 pb0 = arr.verts[FB.verts.x], pb1 = arr.verts[FB.verts.y],
                 pb2 = arr.verts[FB.verts.z];

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
      const vec3 A0 = arr.verts[SA.vertId0];
      const vec3 A1 = arr.verts[SA.vertId1];
      const vec3 B0 = arr.verts[SB.vertId0];
      const vec3 B1 = arr.verts[SB.vertId1];
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
    for (const GroupedEdge& e : CollectGroupedMemberEdges(arr))
      skeleton[e.group].push_back({e.a, e.b});
    for (const auto& seam : arr.seams) {
      const std::pair<vec3, vec3> seg = {arr.verts[seam.vertId0],
                                         arr.verts[seam.vertId1]};
      const int g0 = arr.face2Group[seam.faceId0];
      const int g1 = arr.face2Group[seam.faceId1];
      if (g0 >= 0) skeleton[g0].push_back(seg);
      if (g1 >= 0 && g1 != g0) skeleton[g1].push_back(seg);
    }
    for (const auto& segs : skeleton)
      AppendPairwiseCrossXs(segs, eps, arr.criticalXs);
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
  // groupedEdges is the arr-constant CollectGroupedMemberEdges result, hoisted
  // out of the per-slab loop by the caller (it does not depend on the slab).
  SlabResolver(const SlabResult& slab,
               const std::vector<GroupedEdge>& groupedEdges, double eps)
      : eps_(eps) {
    if (!slab.built) return;
    // Class-iii candidates: edges of grouped faces crossing this slab's
    // section plane (all three edges per member; over-inclusion harmless).
    for (const GroupedEdge& e : groupedEdges) {
      if (e.a.x == e.b.x) continue;  // no crossing trajectory
      if (std::min(e.a.x, e.b.x) > slab.xMid ||
          std::max(e.a.x, e.b.x) < slab.xMid)
        continue;
      planar_.push_back({e.a, e.b});
    }
    // Resolve piece endpoints per EPS-CLUSTER, not per exact vertex: the
    // winding pass constructs T-junction verts exactly and never merges
    // them, so a section can carry sub-eps-adjacent twins.  Resolved
    // independently they can match DIFFERENT tracks, and divergent tracks
    // amplify a sub-eps section gap into a many-eps gap at the cap
    // (corpus-diagnosed).  One ball, one track, every member extends
    // identically - the same one-entity rule as the vert merges.
    std::vector<vec2> pts;
    for (const SweepCapture& piece : slab.pieces) {
      for (const vec2 pt : {piece.from, piece.to}) {
        const auto key = std::make_pair(pt.x, pt.y);
        if (tracks_.count(key)) continue;
        tracks_.emplace(key, Track{true, vec3(0.0), vec3(0.0)});
        pts.push_back(pt);
      }
    }
    const int n = static_cast<int>(pts.size());
    DisjointSets clusterUf(n);
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        if (la::length(pts[i] - pts[j]) <= eps) clusterUf.unite(i, j);
    std::map<int, Track> root2Track;
    for (int i = 0; i < n; ++i) {
      const int root = static_cast<int>(clusterUf.find(i));
      auto it = root2Track.find(root);
      if (it == root2Track.end())
        it = root2Track.emplace(root, Resolve(pts[root], slab)).first;
      tracks_[{pts[i].x, pts[i].y}] = it->second;
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
// loX/hiX are the cap PLANE x's at which each side's chains were bound (the
// pair-canonical critical) - not necessarily the slab's own xLo/xHi.  A slab
// following a run of unbuilt slabs has its lo side bound at the run's canonical
// critical (the preceding built slab's xHi), so its lo strip edge spans the
// skipped run back to that plane (spec CHAIN-PLANE RULE): the cap and both
// adjacent strip boundaries then share one exact plane and closure is
// constructional, not weld-dependent.  Under the strict-FP slab gate the
// spanned run is only ulp-scale, but that exact placement is STILL load-bearing
// (spec STRICT-FP SLAB BUILDING): the money fixture Coplanar_PerpFacesSubEps
// resolves only because the post-run strip corner is placed bitwise at the cap
// plane - emitting at the slab bound leaves the sheet splitter an unpaired fan
// even at an ulp gap.  Set by EmitCaps at bind time; a built slab whose sides
// were not bound trips the count check in EmitStrips before these are read.
struct StripChains {
  std::vector<std::vector<vec2>> lo, hi;
  double loX = 0.0, hiX = 0.0;
};

// Strips stage: each strip zips its two c-side chains, built by the cap
// arrangements at its slab's boundary-pair canonical criticals (spec CAPS
// [R2-fold] third consumer).  No geometry is computed here: the chains carry
// the cap arrangements' vert positions bitwise.  Every built slab has both
// sides bound (its own xHi is always its pair's canonical; its lo is bound at
// the preceding gap's canonical), with one chain per piece - attribution
// failures fail closed in EmitCaps before this runs.
MaybeFatal EmitStrips(std::vector<OutTri3D>& out,
                      const std::vector<SlabResult>& slabs,
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
      return Fatal{FatalReason::NonManifoldEmission,
                   "strip chain count disagrees with piece count"};
    }
    for (size_t k = 0; k < ch.lo.size(); ++k) {
      if (ch.lo[k].empty() || ch.hi[k].empty()) {
        DEBUG_ASSERT(false, logicErr, "empty strip chain");
        return Fatal{FatalReason::NonManifoldEmission, "empty strip chain"};
      }
      // Emit at the chains' CAP PLANE x's (spec CHAIN-PLANE RULE), not the
      // slab's own boundary: a post-gap slab's lo edge spans back across the
      // skipped run to the pair-canonical critical so its corners are bitwise
      // the cap's.  hiX is always the slab's xHi (a slab is the left member of
      // its own pair); only loX can differ, when an unbuilt run precedes -
      // which under the strict-FP slab gate is ulp-scale but still load-bearing
      // for closure (spec STRICT-FP SLAB BUILDING).
      ZipperEmit(out, ch.loX, ch.hiX, ch.lo[k], ch.hi[k]);
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
// `wideRun` is true when a run of unbuilt slabs wider than eps in total
// separates the two built slabs (spec CHAIN-PLANE RULE guard): the chain-plane
// rule then spans that run by linear interpolation, which is exact only for a
// COINCIDENT transition (cap empty).  A non-empty cap here is a genuine macro
// geometry change collapsed into a sub-eps-per-slab interval - fail closed as
// SubEpsFeature.  Under the strict-FP slab gate (spec STRICT-FP SLAB BUILDING)
// an unbuilt slab is ulp-wide, so a run exceeds eps only for a pathological
// pile of ~1000 consecutive adjacent-double criticals; wideRun is therefore
// unreachable on realistic input and this guard is retained as the chain-plane
// rule's fidelity backstop (the alternative would be a silent oracle-wrong
// resolve there), not a boundary any fixture reaches.
// Returns nullopt on success, or the fatal to propagate.
MaybeFatal ComputeCap(std::vector<OutTri3D>& out,
                      std::vector<std::vector<vec2>>* leftChains,
                      std::vector<std::vector<vec2>>* rightChains,
                      Overlap3Counters& cnt, const SlabResult* leftSlab,
                      const SlabResolver* leftResolver,
                      const SlabResult* rightSlab,
                      const SlabResolver* rightResolver, double xCap,
                      bool wideRun, double eps) {
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

  // Cap inputs are twice-constructed (seam endpoints are TriTriSeam
  // constructions, extensions interpolate along tracks), so two derivations
  // of one junction can land several eps apart - above the engine's merge
  // ball at eps, which leaves micro-edges and sliver fans in the output
  // (corpus-falsified on real geometry).  The cap arrangement therefore runs
  // at the construction-noise radius: alpha ~ eps per kernel, depth two,
  // both sides -> 4 eps, doubled for headroom.  Chains use the same radius
  // so strip subdivision matches the cap identity model.
  const double capEps = 8 * eps;
  std::vector<OutEdge> negEdges;
  std::vector<std::vector<vec2>> edgeSubdiv;
  ++cnt.capArrangements;
  const OverlapResult r = RemoveOverlaps2D(
      ces.rawVerts, edges, capEps, /*debug=*/false, WindRule::Add,
      /*trace=*/nullptr, &negEdges,
      (leftChains || rightChains) ? &edgeSubdiv : nullptr);

  // CHAIN-PLANE RULE guard (spec).  Over a run wider than eps the chain-plane
  // rule spans the gap by linear interpolation of the flanking sections; that
  // is exact only when they COINCIDE (region(L) == region(R), so both cap
  // measures are empty - the corpus dead-zone-c case, where interpolating one
  // loop reproduces it).  A NON-EMPTY cap means the flanks differ macro-
  // scopically: a real geometry change squeezed into a sub-eps-per-slab
  // interval wider than eps that no linear span can carry.  Fail closed by
  // name - the single-face SubEpsFeature guard (BuildSlabs) misses this because
  // the change is distributed across straddling faces.  Sub-capEps noise
  // already annihilated in the arrangement, so only genuine macro differences
  // reach here.
  if (wideRun && (!r.edges.empty() || !negEdges.empty()))
    return Fatal{FatalReason::SubEpsFeature,
                 "macro cap content over a skipped run wider than eps"};

  bool loopsClosed = true;
  bool trisOk = true;
  // Emit one signed cap side: walk its retained boundary into loops and
  // triangulate at x=xCap.  cap_plus (L - R) faces +x; cap_minus (R - L) faces
  // -x, hence the flipped winding.  Runs only while the prior side stayed
  // closed and valid; a failure propagates through loopsClosed/trisOk.
  auto emitCapSide = [&](const std::vector<OutEdge>& capEdges, bool flip) {
    if (capEdges.empty() || !loopsClosed || !trisOk) return;
    const Polygons cap = OutEdgesToPolygons(r.verts, capEdges, &loopsClosed);
    trisOk = TriangulateCap(out, cap, xCap, flip, eps);
  };
  emitCapSide(r.edges, /*flip=*/false);
  emitCapSide(negEdges, /*flip=*/true);
  if (!loopsClosed)
    return Fatal{FatalReason::NonManifoldEmission,
                 "cap boundary walk failed to close"};
  if (!trisOk)
    return Fatal{FatalReason::NonManifoldEmission,
                 "cap triangulation returned invalid indices"};

  if (!leftChains && !rightChains) return std::nullopt;
  // Strip chains are PROVENANCE-EXACT (spec PROVENANCE CHAINS): each piece's
  // chain is the engine's own subdivision of that input edge (edgeSubdiv,
  // index-aligned with `edges` = L pieces then R pieces).  Coincident L/R
  // pieces (shared merged endpoints) get identical interior sequences by
  // construction, so strips pair across the critical exactly.
  const size_t nL = ces.lPiece2RawVerts.size();
  if (leftChains)
    for (size_t k = 0; k < nL; ++k) (*leftChains)[k] = edgeSubdiv[k];
  if (rightChains)
    for (size_t k = 0; k < ces.rPiece2RawVerts.size(); ++k)
      (*rightChains)[k] = edgeSubdiv[nL + k];
  return std::nullopt;
}

// Caps-stage driver: PAIR-CANONICAL emission.  Each adjacent-built-slab pair
// emits exactly ONE cap, at the pair's canonical critical ci == li+1 - the
// first critical of the gap between the slabs (the critical itself under
// direct adjacency).  Adjacent built slabs are found by INDEX: slab bounds are
// these exact critical values by construction, so no slack enters the lookup.
// Both sides of the pair consume that one arrangement, so strips weld exactly
// across the critical even where a sub-eps run's two arrangements disagree
// about a shared cancelled edge (geometry is discontinuous at a critical).
MaybeFatal EmitCaps(std::vector<OutTri3D>& out,
                    std::vector<StripChains>& chains, Overlap3Counters& cnt,
                    const std::vector<SlabResult>& slabs,
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
    // Non-canonical criticals (interior to an unbuilt run) emit NOTHING:
    // their cap would be the same L-R difference the run's canonical cap
    // already computed, re-derived at an x within eps of it - and the two
    // arrangements coincide only up to snap noise, which track slopes can
    // amplify past eps (corpus-falsified: the doubled-cap emission-closure
    // defect).  One cap per adjacent-built-slab pair, at the critical where
    // the chains bind.
    if (ci != li + 1) continue;
    // Record the cap PLANE x for the strip chains this cap binds (spec
    // CHAIN-PLANE RULE): the left slab's hi edge and the right slab's lo edge
    // both live at crits[ci], so their strips emit there and share the cap's
    // plane exactly - even when a skipped sub-eps run separates the slabs.
    if (left) chains[li].hiX = crits[ci];
    if (right) chains[ri].loX = crits[ci];
    // A WIDE skipped run separates two BUILT slabs when the gap between the
    // left slab's xHi (= crits[li+1], the cap plane) and the right slab's xLo
    // (= crits[ri]) exceeds eps (spec CHAIN-PLANE RULE guard).  Under the
    // strict-FP slab gate every unbuilt run is ulp-scale, so this is
    // unreachable on realistic input (see ComputeCap) and evaluates false; it
    // is kept as the chain-plane rule's fidelity backstop.  Exterior caps (one
    // flank null) are the genuine geometry end - never guarded.
    const bool wideRun = left && right && crits[ri] - crits[li + 1] > eps;
    if (auto fatal = ComputeCap(out, left ? &chains[li].hi : nullptr,
                                right ? &chains[ri].lo : nullptr, cnt, left,
                                li >= 0 ? &resolvers[li] : nullptr, right,
                                ri < nSlabs ? &resolvers[ri] : nullptr,
                                crits[ci], wideRun, eps))
      return fatal;
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// Assembly: weld verts, drop exact-duplicate triangles, build Impl.
// ---------------------------------------------------------------------------

// Split geometrically-welded touching sheets back into topologically
// separate manifolds - the epsilon-valid posture Boolean3 itself emits for
// touching solids (coincident geometry, separate topology).  The eps-weld
// fuses surfaces that touch on measure-zero sets (edge-on-face, edge-edge,
// self-touch) into complexes with 2k-halfedge fan edges and non-disk vertex
// links.  Resolution in two steps:
//   1. RADIAL PAIRING per fan edge: sort incident faces by angle around the
//      edge axis.  A forward halfedge (lo->hi) has its material just BELOW
//      its angle, a backward one just ABOVE (from the outward-normal
//      convention), so material wedges alternate with empty ones and each
//      backward halfedge pairs with the NEXT forward one CCW.
//   2. VERTEX SPLIT: union corners around each vert through PAIRED
//      halfedges only; each connected component becomes its own vert copy.
// Returns false (the caller fails closed) when the fan cannot be paired:
// unbalanced counts, a sliver third-vert on the edge line, radial ties
// (tangent sheets - either pairing is a coin flip and the wrong one keeps
// volume while garbling topology), or a non-alternating pattern
// (overlapping material - upstream resolution failed).
bool SplitTouchingSheets(std::vector<vec3>& verts, Vec<ivec3>& tv) {
  const int nTri = static_cast<int>(tv.size());
  auto heFrom = [&](int h) { return tv[h / 3][h % 3]; };
  auto heTo = [&](int h) { return tv[h / 3][(h % 3 + 1) % 3]; };

  std::map<std::pair<int, int>, std::vector<int>> edge2He;
  for (int h = 0; h < 3 * nTri; ++h) {
    const int a = heFrom(h), b = heTo(h);
    edge2He[{std::min(a, b), std::max(a, b)}].push_back(h);
  }

  std::vector<int> pairedHe(3 * nTri, -1);
  for (const auto& [edge, hes] : edge2He) {
    std::vector<int> fwd, bwd;
    for (const int h : hes) (heFrom(h) == edge.first ? fwd : bwd).push_back(h);
    if (fwd.size() != bwd.size()) return false;
    if (fwd.size() == 1) {
      pairedHe[fwd[0]] = bwd[0];
      pairedHe[bwd[0]] = fwd[0];
      continue;
    }

    const vec3 pa = verts[edge.first], pb = verts[edge.second];
    const vec3 ax = la::normalize(pb - pa);
    struct RingEntry {
      double angle;
      int he;
      bool fwd;
    };
    std::vector<RingEntry> ring;
    ring.reserve(hes.size());
    vec3 u(0.0), v(0.0);
    for (const int h : hes) {
      const int c = tv[h / 3][(h % 3 + 2) % 3];
      vec3 d = verts[c] - pa;
      d -= la::dot(d, ax) * ax;
      const double len = la::length(d);
      if (len == 0.0) return false;  // sliver: third vert on the edge line
      d /= len;
      if (ring.empty()) {
        u = d;
        v = la::cross(ax, u);
      }
      double angle = std::atan2(la::dot(d, v), la::dot(d, u));
      if (angle < 0.0) angle += kTwoPi;
      ring.push_back({angle, h, heFrom(h) == edge.first});
    }
    std::sort(ring.begin(), ring.end(),
              [](const RingEntry& a, const RingEntry& b) {
                return a.angle < b.angle;
              });
    const int k = static_cast<int>(ring.size());
    // Ties (including the wraparound pair) and alternation.  Genuine
    // touching contacts have macro dihedral separation; kAngleTie guards
    // the tangent-sheet coin flip.
    constexpr double kAngleTie = 1e-9;
    for (int i = 0; i < k; ++i) {
      const RingEntry& cur = ring[i];
      const RingEntry& nxt = ring[(i + 1) % k];
      const double gap =
          (i + 1 < k) ? nxt.angle - cur.angle : nxt.angle + kTwoPi - cur.angle;
      if (gap < kAngleTie) return false;
      if (cur.fwd == nxt.fwd) return false;  // material overlap
    }
    for (int i = 0; i < k; ++i) {
      if (ring[i].fwd) continue;
      const RingEntry& partner = ring[(i + 1) % k];  // next CCW is forward
      pairedHe[ring[i].he] = partner.he;
      pairedHe[partner.he] = ring[i].he;
    }
  }

  // Vertex split by paired-fan connectivity.
  DisjointSets cornerUf(3 * nTri);
  for (int h = 0; h < 3 * nTri; ++h) {
    const int p = pairedHe[h];
    if (p < 0) return false;
    if (p < h) continue;
    const int t1 = h / 3, i1 = h % 3, t2 = p / 3, i2 = p % 3;
    cornerUf.unite(3 * t1 + i1, 3 * t2 + (i2 + 1) % 3);  // at heFrom(h)
    cornerUf.unite(3 * t1 + (i1 + 1) % 3, 3 * t2 + i2);  // at heTo(h)
  }
  std::map<std::pair<int, int>, int> vertComp2Out;
  std::vector<bool> vertKept(verts.size(), false);
  for (int corner = 0; corner < 3 * nTri; ++corner) {
    const int t = corner / 3, k = corner % 3;
    const int vsrc = tv[t][k];
    const int root = static_cast<int>(cornerUf.find(corner));
    auto it = vertComp2Out.find({vsrc, root});
    if (it == vertComp2Out.end()) {
      int out;
      if (!vertKept[vsrc]) {
        vertKept[vsrc] = true;
        out = vsrc;
      } else {
        out = static_cast<int>(verts.size());
        verts.push_back(verts[vsrc]);
      }
      it = vertComp2Out.emplace(std::make_pair(vsrc, root), out).first;
    }
    tv[t][k] = it->second;
  }
  return true;
}

// Integer cell key for the assembly weld's uniform hash grid (cell = eps).
struct GridCell {
  int64_t x, y, z;
  bool operator==(const GridCell& o) const {
    return x == o.x && y == o.y && z == o.z;
  }
};
struct GridCellHash {
  size_t operator()(const GridCell& c) const {
    size_t h = 1469598103934665603ull;  // FNV-1a mix of the three coords
    for (const int64_t v : {c.x, c.y, c.z}) {
      h = (h ^ static_cast<size_t>(v)) * 1099511628211ull;
    }
    return h;
  }
};

StageResult<Manifold::Impl> BuildImpl(const std::vector<OutTri3D>& tris,
                                      double eps) {
  if (tris.empty()) return StageResult<Manifold::Impl>::Ok(Manifold::Impl{});

  // Collect verts with an eps-weld.  A uniform hash grid (cell = eps) over the
  // growing vert list makes each query O(1) amortized instead of O(n); it
  // returns the MINIMUM-INDEX vert within eps, identical to the linear
  // first-match scan it replaces - a point within eps of the query lies in the
  // query's cell or one of its 26 neighbours, so the 3x3x3 sweep sees every
  // candidate and the min-index tie-break reproduces first-match exactly.
  std::vector<vec3> verts;
  std::unordered_map<GridCell, std::vector<int>, GridCellHash> grid;
  auto cellOf = [&](vec3 p) -> GridCell {
    return {static_cast<int64_t>(std::floor(p.x / eps)),
            static_cast<int64_t>(std::floor(p.y / eps)),
            static_cast<int64_t>(std::floor(p.z / eps))};
  };
  auto getVertIdx = [&](vec3 p) -> int {
    const GridCell c = cellOf(p);
    int best = -1;
    for (const int64_t dx : {-1, 0, 1})
      for (const int64_t dy : {-1, 0, 1})
        for (const int64_t dz : {-1, 0, 1}) {
          const auto it = grid.find({c.x + dx, c.y + dy, c.z + dz});
          if (it == grid.end()) continue;
          for (const int i : it->second)
            if ((best < 0 || i < best) && la::length(verts[i] - p) <= eps)
              best = i;
        }
    if (best >= 0) return best;
    const int id = static_cast<int>(verts.size());
    verts.push_back(p);
    grid[c].push_back(id);
    return id;
  };

  // Filter degenerate and duplicate triangles.
  // Degenerate: strip quads whose corners collapse within eps produce v0==v1
  // etc., which would crash CreateHalfedges.
  // Duplicate: after the weld two emitted triangles can share the same vertex
  // triple; drop exact duplicates (keeping one) - a repeated face would break
  // the 2-manifold topology.
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

  // Touching sheets separate BEFORE the topology is built (spec COPLANAR
  // implementation close: touching contacts).
  if (!SplitTouchingSheets(verts, tv)) {
    return StageResult<Manifold::Impl>::Fatal(FatalReason::NonManifoldEmission,
                                              "unresolvable sheet contact");
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
  const std::vector<GroupedEdge> groupedEdges = CollectGroupedMemberEdges(arr);
  for (const SlabResult& slab : slabs)
    resolvers.emplace_back(slab, groupedEdges, eps);

  std::vector<OutTri3D> emitted;
  std::vector<StripChains> chains(slabs.size());
  if (auto capFatal =
          EmitCaps(emitted, chains, cnt, slabs, resolvers, crits, eps)) {
    result.fatal = capFatal->reason;
    result.detail = std::move(capFatal->detail);
    result.counters = cnt;
    return result;
  }
  if (auto stripFatal = EmitStrips(emitted, slabs, chains)) {
    result.fatal = stripFatal->reason;
    result.detail = std::move(stripFatal->detail);
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

// The shared prefix of both entry points: resolve eps (in/out), canonicalize,
// and find seams.  A SubEpsInput fatal rides SeamsResult.fatal; trivially-empty
// input (canonicalize dropped every face) returns a fatal-free result whose
// arr.faces is empty, which the callers emit as an empty manifold.
SeamsResult PrepareArrangement(const Manifold::Impl& in, double& eps,
                               Overlap3Counters& cnt) {
  SeamsResult res;
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    res.fatal = FatalReason::SubEpsInput;
    res.detail = "epsilon not computable";
    return res;
  }
  const CanonicalGeometry canon = Canonicalize(in, eps);
  if (canon.faces.empty()) return res;  // arr.faces stays empty
  return FindSeams(canon, eps, cnt);
}

}  // namespace

// ---------------------------------------------------------------------------
// Public entry points.
// ---------------------------------------------------------------------------

Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps) {
  Overlap3Result result;
  SeamsResult seams = PrepareArrangement(in, eps, result.counters);
  if (seams.fatal.has_value()) {
    result.fatal = seams.fatal;
    result.detail = seams.detail;
    return result;
  }
  if (seams.arr.faces.empty()) {
    result.impl = Manifold::Impl{};
    return result;
  }
  return SweepEmit(seams.arr, eps, result.counters);
}

Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps) {
  Overlap3Internals out;
  SeamsResult seams = PrepareArrangement(in, eps, out.counters);
  if (seams.fatal.has_value()) {
    out.fatal = seams.fatal;
    out.detail = seams.detail;
    return out;
  }
  if (seams.arr.faces.empty()) return out;
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

// ---------------------------------------------------------------------------
// Regularization operator (docs/Regularize3D.md) - Stage-1 GATE + DISPATCH.
// ---------------------------------------------------------------------------

namespace {

// Forward decl (defined below): per-face exactly-coplanar overlap cluster id
// (or -1).  Used by DecomposeComponents to merge coplanar-overlapping
// components and by GateComponent to route a coplanar-overlap component to B.
std::vector<int> DetectCoplanarClusters(const Manifold::Impl& in);

// Split `in` into connected components by halfedge connectivity - the Decompose
// primitive (constructors.cpp:455) mirrored at the Impl level so the operator
// never round-trips through the CSG layer.  Each returned component is a
// finished Impl (bbox/normals/collider) whose epsilon_ is pinned to the
// resolved global machine scale, so the per-component gate (IsSelfIntersecting
// reads collider_, faceNormal_, epsilon_) runs directly.  A single connected
// input returns exactly one component that IS a copy of `in`, so a clean
// input of one component passes through unchanged.
//
// COPLANAR-OVERLAP MERGE (docs/Regularize3D.md coplanar axis): two connectivity
// components that overlap in a shared plane (a buried plug with a coincident
// cap, coplanar stacked boxes) are a genuine SELF-OVERLAP defect the coplanar
// fold resolves - but decompose-by-connectivity would split them so the
// per-component gate never sees the overlap (it is CROSS-component).  When
// there is more than one component, the exactly-coplanar overlap clusters are
// detected on the whole input and the components sharing a cluster are UNITED,
// so the overlap becomes internal to one super-component (the fold's resolve
// target). Scoped to genuine 2D-AREA coplanar overlap (DetectCoplanarClusters'
// overlap2D witness): touching contacts (edge/vertex, zero area) and disjoint
// objects never cluster, so the non-fusion posture holds for everything but a
// real coplanar self-overlap.
std::vector<Manifold::Impl> DecomposeComponents(const Manifold::Impl& in,
                                                double eps) {
  std::vector<Manifold::Impl> out;
  const int numVert = static_cast<int>(in.NumVert());
  if (numVert == 0 || in.halfedge_.size() == 0) return out;

  DisjointSets uf(numVert);
  for (size_t e = 0; e < in.halfedge_.size(); ++e)
    if (in.halfedge_.IsForward(e))
      uf.unite(in.halfedge_.Start(static_cast<int>(e)),
               in.halfedge_.End(static_cast<int>(e)));
  std::vector<int> vertLabel;
  int numComponents = uf.connectedComponents(vertLabel);

  // Coplanar-overlap merge (see the header comment): only when the connectivity
  // split produced more than one component - a single component is already the
  // finest unit and running the O(nTri^2) cluster scan on it would be pure
  // cost. Uniting one vertex of each cluster member folds the components that
  // share a coplanar cluster together; a cluster entirely within one component
  // is a no-op (its verts are already connected).
  if (numComponents > 1) {
    const std::vector<int> face2cluster = DetectCoplanarClusters(in);
    std::vector<int> clusterRep;  // first vert seen per cluster id
    for (int f = 0; f < static_cast<int>(face2cluster.size()); ++f) {
      const int c = face2cluster[f];
      if (c < 0) continue;
      const int v = in.halfedge_.Start(3 * f);
      if (c >= static_cast<int>(clusterRep.size()))
        clusterRep.resize(c + 1, -1);
      if (clusterRep[c] < 0)
        clusterRep[c] = v;
      else
        uf.unite(clusterRep[c], v);
    }
    numComponents = uf.connectedComponents(vertLabel);
  }

  if (numComponents == 1) {
    // The whole input is one component; copy it through unchanged so a clean
    // single-component input is bitwise pass-through.  Its finished state
    // (collider/normals from construction) drives the gate; only epsilon_ is
    // pinned to the resolved machine scale.
    out.push_back(in);
    out.back().epsilon_ = eps;
    return out;
  }

  const int numTri = static_cast<int>(in.NumTri());
  for (int c = 0; c < numComponents; ++c) {
    // Compact this component's verts; vertNew2Old feeds ReindexVerts.
    Vec<int> vertNew2Old;
    for (int v = 0; v < numVert; ++v)
      if (vertLabel[v] == c) vertNew2Old.push_back(v);
    if (vertNew2Old.empty()) continue;

    // Faces whose first vert carries this label; halfedge connectivity
    // guarantees all three verts of a face share it.
    Vec<int> faceNew2Old;
    for (int f = 0; f < numTri; ++f)
      if (vertLabel[in.halfedge_.Start(3 * f)] == c) faceNew2Old.push_back(f);
    if (faceNew2Old.empty()) continue;

    Manifold::Impl comp;
    comp.vertPos_.resize(vertNew2Old.size());
    for (size_t i = 0; i < vertNew2Old.size(); ++i)
      comp.vertPos_[i] = in.vertPos_[vertNew2Old[i]];
    comp.GatherFaces(in, faceNew2Old);  // halfedges with OLD vert ids
    comp.ReindexVerts(vertNew2Old, in.NumVert());  // remap to the compacted ids
    // Finish so the gate has bbox/collider/normals; pin epsilon_ to the global
    // machine scale for a consistent 2*eps relaxation across components.
    comp.CalculateBBox();
    comp.SetEpsilon();
    comp.SortGeometry();
    comp.SetNormalsAndCoplanar();
    comp.epsilon_ = eps;
    out.push_back(std::move(comp));
  }
  return out;
}

enum class GateVerdict { Clean, Dirty, Invalid };

// The per-component gate (docs/Regularize3D.md step 2): valid AND
// non-self-intersecting.  Components of a valid oriented 2-manifold are
// themselves valid, so Invalid is defensive (never expected on a valid input).
//
// R2 CLEAN-BIAS CAVEAT (docs/Regularize3D.md R2): IsSelfIntersecting is
// systematically CLEAN-biased - its 2*eps shares-vertex relaxation SUPPRESSES
// near-miss detection (it returns non-intersecting when an eps normal nudge
// separates the pair), it does not flag near-misses.  So a genuine crossing
// whose two verts sit in the thin band just outside the eps weld can pass the
// gate, EARLY-EXIT as clean, and carry an unregularized self-overlap through
// silently (the R2(i) = R1 blind spot).  This gate does not close that narrow
// hole; it is a recorded open, not a claim of completeness.
GateVerdict GateComponent(const Manifold::Impl& comp) {
  if (!comp.IsManifold() || !comp.Is2Manifold()) return GateVerdict::Invalid;
  if (comp.IsSelfIntersecting()) return GateVerdict::Dirty;
  // A component that PASSES the self-intersection test can still carry a pure
  // COPLANAR overlap (coincident / overlapping coplanar faces):
  // IsSelfIntersecting does NOT flag coplanar coincidence (docs/Regularize3D.md
  // R2(i)), so such a component would early-exit "clean" yet is NOT the
  // {w_S>=1} boundary (a buried plug's coincident cap doubles the cover).
  // DetectCoplanarClusters (bbox-overlap prefilter, then exact orient3d
  // coplanarity + a 2D-area overlap witness) routes it to B, where the
  // exact-coplanar fold consumes the overlap.  A clean solid with no coplanar
  // overlap detects nothing and stays Clean (bitwise pass-through); the cost is
  // the prefiltered scan, proportional to the input.
  for (int c : DetectCoplanarClusters(comp))
    if (c >= 0) return GateVerdict::Dirty;
  return GateVerdict::Clean;
}

// ---------------------------------------------------------------------------
// Candidate B mechanism (docs/Regularize3D.md "B's mechanism"), ported from the
// FRAGMENT-VALIDATED reference (v5b fragment drivers + v5b-r3/r4 notebooks):
// operand-agnostic ENUMERATION (level-0 pierce predicates through a static
// Shewchuk filter) + coupled integer-delta WINDING.  Every crossing DECISION is
// a level-0 orient3d on INPUT coordinates.  These are the substrate of B; the
// cell-complex + halfedge {w_S>=1} boundary EMISSION (THE BUILD) sits on top
// and is the doc's named largest-unbuilt-piece.
// ---------------------------------------------------------------------------

// orient3d on double INPUT coords through the Shewchuk STATIC error-bound
// filter (o3derrboundA = (7 + 56u)u, u = 2^-53).  Returns the CERTIFIED sign
// (+/-1) when |det| exceeds the permanent-scaled error bound; returns 0 when
// the filtered sign is UNCERTAIN (sub-bound or exact-zero).  The exact-Fraction
// fallback that the fragment used for the 0 case is NOT ported (net-new exact
// kernel); B treats an uncertain deciding predicate as a fail-closed boundary
// rather than guess a sign.  On the corpus single-shell self-intersectors this
// filter certifies every enumeration predicate (v5b-r4: siA/siB 100% certified,
// zero fallback).
inline int Orient3DFilterSign(const vec3& a, const vec3& b, const vec3& c,
                              const vec3& d) {
  const vec3 ad = a - d, bd = b - d, cd = c - d;
  const double bc = bd.y * cd.z, cb = cd.y * bd.z;
  const double ca = cd.y * ad.z, ac = ad.y * cd.z;
  const double ab = ad.y * bd.z, ba = bd.y * ad.z;
  const double det = ad.x * (bc - cb) + bd.x * (ca - ac) + cd.x * (ab - ba);
  const double perm = (std::abs(bc) + std::abs(cb)) * std::abs(ad.x) +
                      (std::abs(ca) + std::abs(ac)) * std::abs(bd.x) +
                      (std::abs(ab) + std::abs(ba)) * std::abs(cd.x);
  constexpr double u = 0x1p-53;
  const double errb = (7.0 + 56.0 * u) * u * perm;
  if (errb > 0.0 && std::abs(det) > errb) return det > 0.0 ? 1 : -1;
  return 0;  // uncertain -> exact fallback (unbuilt); caller fails closed
}

// Does edge (u,v) pierce the INTERIOR of triangle (a,b,c)?  Level-0, decomposed
// into the D1 straddle (orient3d of the plane vs each endpoint) and the D3
// edge-edge z-order (orient3d of the edge vs each triangle edge) - the v5b-r4
// P4 centerpiece.  Returns 1 (genuine pierce), 0 (no pierce), or -1 (a deciding
// predicate was filter-uncertain / an exact-zero boundary: the SoS axis).
int EdgePiercesTri(const vec3& u, const vec3& v, const vec3& a, const vec3& b,
                   const vec3& c) {
  const int su = Orient3DFilterSign(a, b, c, u);
  const int sv = Orient3DFilterSign(a, b, c, v);
  if (su == 0 || sv == 0) return -1;  // endpoint on/near the plane
  if (su == sv) return 0;             // both same side -> no straddle
  const int o1 = Orient3DFilterSign(u, v, a, b);
  const int o2 = Orient3DFilterSign(u, v, b, c);
  const int o3 = Orient3DFilterSign(u, v, c, a);
  if (o1 == 0 || o2 == 0 || o3 == 0) return -1;
  return (o1 == o2 && o2 == o3) ? 1 : 0;
}

struct BEnumeration {
  int seamCount = 0;  // genuine non-adjacent self-crossings
  int boundaryTouchPairs =
      0;  // pairs hitting an exact-zero/uncertain predicate
};

// Enumerate the component's genuine self-crossing arrangement: bbox broadphase
// + shared-vertex (self-adjacency) skip + level-0 pierce test.  Reproduces the
// fragment's seam set (v5b-r4: siA/siB = 338 seams, 0 boundary-touch).  Brute
// bbox broadphase is O(F^2) (17k tris ~0.2s); the collider_ broadphase is the
// perf path (unbuilt here - correctness-first).
BEnumeration EnumerateSelfCrossings(const Manifold::Impl& in) {
  BEnumeration out;
  const int nTri = static_cast<int>(in.NumTri());
  std::vector<std::array<vec3, 3>> tri(nTri);
  std::vector<std::array<int, 3>> vid(nTri);
  std::vector<vec3> lo(nTri), hi(nTri);
  for (int t = 0; t < nTri; ++t) {
    for (int k = 0; k < 3; ++k) {
      vid[t][k] = in.halfedge_.Start(3 * t + k);
      tri[t][k] = in.vertPos_[vid[t][k]];
    }
    lo[t] = la::min(la::min(tri[t][0], tri[t][1]), tri[t][2]);
    hi[t] = la::max(la::max(tri[t][0], tri[t][1]), tri[t][2]);
  }
  auto bboxOverlap = [&](int i, int j) {
    return !(hi[i].x < lo[j].x || hi[j].x < lo[i].x || hi[i].y < lo[j].y ||
             hi[j].y < lo[i].y || hi[i].z < lo[j].z || hi[j].z < lo[i].z);
  };
  auto sharesVert = [&](int i, int j) {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (vid[i][a] == vid[j][b]) return true;
    return false;
  };
  for (int i = 0; i < nTri; ++i) {
    for (int j = i + 1; j < nTri; ++j) {
      if (!bboxOverlap(i, j)) continue;
      if (sharesVert(i, j)) continue;  // self-adjacency skip (S4a)
      const auto& A = tri[i];
      const auto& B = tri[j];
      bool genuine = false, boundary = false;
      for (int e = 0; e < 3 && !genuine; ++e) {
        const int r = EdgePiercesTri(A[e], A[(e + 1) % 3], B[0], B[1], B[2]);
        if (r == 1) genuine = true;
        if (r == -1) boundary = true;
      }
      for (int e = 0; e < 3 && !genuine; ++e) {
        const int r = EdgePiercesTri(B[e], B[(e + 1) % 3], A[0], A[1], A[2]);
        if (r == 1) genuine = true;
        if (r == -1) boundary = true;
      }
      if (genuine)
        ++out.seamCount;
      else if (boundary)
        ++out.boundaryTouchPairs;
    }
  }
  return out;
}

// EXACT-COPLANAR FOLD, cluster detection (docs/Regularize3D.md coplanar axis).
// Union non-self-adjacent, bbox-overlapping faces that are EXACTLY coplanar -
// every vertex of each lies on the other's plane, decided by the level-0
// orient3d filter (all six cross-checks certified 0).  The filter returns 0
// only when the vertex is within ~1 ULP (relative) of the plane, i.e. the
// coplanarity gap is far below eps (the machine weld radius), so projecting the
// cluster onto one plane is eps-valid.  The NEAR-coplanar thin band (gap above
// the filter's error bound but below eps) has a NONZERO filter sign, is NOT
// clustered here, and stays the transversal / fail-closed residue the mission
// leaves open.  Returns a per-face cluster id, or -1 for a face in no
// multi-face coplanar cluster (the ordinary transversal path).
std::vector<int> DetectCoplanarClusters(const Manifold::Impl& in) {
  const int nTri = static_cast<int>(in.NumTri());
  std::vector<std::array<vec3, 3>> tri(nTri);
  std::vector<std::array<int, 3>> vid(nTri);
  std::vector<vec3> lo(nTri), hi(nTri);
  for (int t = 0; t < nTri; ++t) {
    for (int k = 0; k < 3; ++k) {
      vid[t][k] = in.halfedge_.Start(3 * t + k);
      tri[t][k] = in.vertPos_[vid[t][k]];
    }
    lo[t] = la::min(la::min(tri[t][0], tri[t][1]), tri[t][2]);
    hi[t] = la::max(la::max(tri[t][0], tri[t][1]), tri[t][2]);
  }
  auto bboxOverlap = [&](int i, int j) {
    return !(hi[i].x < lo[j].x || hi[j].x < lo[i].x || hi[i].y < lo[j].y ||
             hi[j].y < lo[i].y || hi[i].z < lo[j].z || hi[j].z < lo[i].z);
  };
  auto sharesVert = [&](int i, int j) {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (vid[i][a] == vid[j][b]) return true;
    return false;
  };
  auto coplanar = [&](int i, int j) {
    for (int k = 0; k < 3; ++k)
      if (Orient3DFilterSign(tri[i][0], tri[i][1], tri[i][2], tri[j][k]) != 0)
        return false;
    for (int k = 0; k < 3; ++k)
      if (Orient3DFilterSign(tri[j][0], tri[j][1], tri[j][2], tri[i][k]) != 0)
        return false;
    return true;
  };
  // Do the two coplanar triangles share positive 2D area?  Only genuinely
  // OVERLAPPING coplanar faces need folding; the coplanar tiles of one flat
  // face (an annulus, a subdivided facet) merely abut and must NOT cluster - a
  // vertex strictly inside the other, or a properly-crossing edge pair, is the
  // area-overlap witness (triangles are convex, so this is exhaustive).
  auto overlap2D = [&](int i, int j) {
    const vec3 e1 = la::normalize(tri[i][1] - tri[i][0]);
    const vec3 nrm = la::cross(tri[i][1] - tri[i][0], tri[i][2] - tri[i][0]);
    const double nl = la::length(nrm);
    if (!(nl > 0.0)) return false;
    const vec3 e2 = la::cross(nrm / nl, e1);
    const vec3 o = tri[i][0];
    auto pr = [&](const vec3& P) {
      return vec2(la::dot(P - o, e1), la::dot(P - o, e2));
    };
    vec2 A[3] = {pr(tri[i][0]), pr(tri[i][1]), pr(tri[i][2])};
    vec2 B[3] = {pr(tri[j][0]), pr(tri[j][1]), pr(tri[j][2])};
    auto cr = [](const vec2& u, const vec2& v) {
      return u.x * v.y - u.y * v.x;
    };
    auto strictIn = [&](const vec2& p, const vec2* T) {
      const double d0 = cr(T[1] - T[0], p - T[0]);
      const double d1 = cr(T[2] - T[1], p - T[1]);
      const double d2 = cr(T[0] - T[2], p - T[2]);
      const bool neg = d0 < 0 || d1 < 0 || d2 < 0;
      const bool pos = d0 > 0 || d1 > 0 || d2 > 0;
      return !(neg && pos) && d0 != 0 && d1 != 0 && d2 != 0;
    };
    for (int k = 0; k < 3; ++k)
      if (strictIn(A[k], B) || strictIn(B[k], A)) return true;
    auto proper = [&](const vec2& p1, const vec2& p2, const vec2& p3,
                      const vec2& p4) {
      const double d1 = cr(p2 - p1, p3 - p1), d2 = cr(p2 - p1, p4 - p1);
      const double d3 = cr(p4 - p3, p1 - p3), d4 = cr(p4 - p3, p2 - p3);
      return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 &&
             d2 != 0 && d3 != 0 && d4 != 0;
    };
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (proper(A[a], A[(a + 1) % 3], B[b], B[(b + 1) % 3])) return true;
    return false;
  };
  DisjointSets uf(nTri);
  bool any = false;
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (!bboxOverlap(i, j) || sharesVert(i, j)) continue;
      if (coplanar(i, j) && overlap2D(i, j)) {
        uf.unite(i, j);
        any = true;
      }
    }
  std::vector<int> face2cluster(nTri, -1);
  if (!any) return face2cluster;
  std::map<int, int> rootCount;
  for (int f = 0; f < nTri; ++f) ++rootCount[static_cast<int>(uf.find(f))];
  std::map<int, int> rootId;
  int nc = 0;
  for (int f = 0; f < nTri; ++f) {
    const int r = static_cast<int>(uf.find(f));
    if (rootCount[r] > 1) {
      auto it = rootId.find(r);
      if (it == rootId.end()) it = rootId.emplace(r, nc++).first;
      face2cluster[f] = it->second;
    }
  }
  return face2cluster;
}

// Coupled soup winding w_S(p): the signed count of oriented-face crossings on
// the ray p->seed, every crossing decided by level-0 orient3d through the
// static filter (the Winding03 discipline, boolean3.cpp:388).  The delta per
// crossed face is sign(dot(seed-p, n_f)) on the input normal - a +/-1 integer,
// FP-safe by construction.  Returns nullopt if any deciding predicate was
// filter-uncertain (grazed a vertex/edge): the caller re-seeds or fails closed.
std::optional<int> WindingAt(const Manifold::Impl& in, const vec3& p,
                             const vec3& seed) {
  int w = 0;
  const int nTri = static_cast<int>(in.NumTri());
  for (int t = 0; t < nTri; ++t) {
    const vec3 a = in.vertPos_[in.halfedge_.Start(3 * t)];
    const vec3 b = in.vertPos_[in.halfedge_.Start(3 * t + 1)];
    const vec3 c = in.vertPos_[in.halfedge_.Start(3 * t + 2)];
    const int da = Orient3DFilterSign(a, b, c, p);
    const int db = Orient3DFilterSign(a, b, c, seed);
    if (da == 0 || db == 0) return std::nullopt;
    if (da == db) continue;  // p and seed on the same side of the plane
    const int o1 = Orient3DFilterSign(p, seed, a, b);
    const int o2 = Orient3DFilterSign(p, seed, b, c);
    const int o3 = Orient3DFilterSign(p, seed, c, a);
    if (o1 == 0 || o2 == 0 || o3 == 0) return std::nullopt;
    if (o1 == o2 && o2 == o3) {
      const vec3 n = la::cross(b - a, c - a);
      w += (la::dot(seed - p, n) > 0.0) ? 1 : -1;
    }
  }
  return w;
}

// Robust soup winding: try the coupled ray winding from a few unrelated seeds
// and take the first that grazes no vertex/edge (the winding is single-valued
// off-surface, so any certified seed is authoritative).  nullopt only if EVERY
// seed hit a filter-uncertain deciding predicate (SoS / near-degenerate).
std::optional<int> RobustWinding(const Manifold::Impl& in, const vec3& p,
                                 const std::vector<vec3>& seeds) {
  for (const vec3& s : seeds) {
    const std::optional<int> w = WindingAt(in, p, s);
    if (w) return w;
  }
  return std::nullopt;
}

// ---------------------------------------------------------------------------
// THE BUILD: {w_S>=1} halfedge-boundary EMISSION by per-face 2D arrangement.
// The steer's reuse (ComputeCap already runs this exact shape for the sweep's
// caps): each crossed face projects its triangle boundary + its seam segments
// into the face plane; RemoveOverlaps2D (the Smith sweep) arranges them and the
// Add winding rule retains the sub-region on the {w_S>=1} boundary; the
// retained loops triangulate and emit at canonical 3D positions.
//
// Retention rule (uniform, derived): a sub-face of an oriented mult-1 face f is
// on d{w_S>=1} iff w_S on the +n_f side == 0, kept with the ORIGINAL
// orientation
// - because w_below = w_above + 1 always, so the solid is always on the -n_f
// side at a retained face (no flips).  Encoded as a 2D winding: triangle
// boundary CCW mult +1, each seam mult -1 with its higher-winding
// (inside-the-other-lump) side on the LEFT, so Add (net > 0) keeps {1 - G(q) >
// 0} = {G(q) == 0}.  The seam sign is DERIVED from the crossing face's normal
// (not a hand global sign); the output is re-gated + volume-checked, so a wrong
// sign fails closed, never wrong.
//
// ONCE-ONLY construction (doc R1): each seam endpoint is a pierce point keyed
// by (undirected mesh edge, pierced triangle), built ONCE and shared by both
// faces of the seam AND by the two seams that chain at it; it lies on the two
// faces' plane-intersection line, so its 2D projection is exact in either face
// and the cross-face weld is bit-identical.  The winding half stays
// exact-by-integer.
// ---------------------------------------------------------------------------

// One seam segment as seen from a specific face: its two canonical 3D endpoints
// (on this face's plane) and the OUTWARD normal of the crossing face (the sign
// source that orients the 2D edge's winding contribution).  p*Interior flags an
// endpoint that lies in this face's INTERIOR (the crossing face's edge pierced
// this triangle) rather than on this face's boundary edge - the marker the spur
// prune reads (an interior degree-1 end dangles and does not bound {G==0}).
struct BuildSeam {
  vec3 p0, p1;
  vec3 nOther;
  bool p0Interior, p1Interior;
};

using PierceKey =
    std::tuple<int, int, int>;  // (min edge vert, max, pierced tri)

// The intersection point of segment (u,v) with the plane of triangle (a,b,c).
// Precondition: the segment straddles the plane (EdgePiercesTri == 1), so the
// denominator is nonzero.
vec3 SegPlanePoint(const vec3& u, const vec3& v, const vec3& a, const vec3& b,
                   const vec3& c) {
  const vec3 n = la::cross(b - a, c - a);
  const double du = la::dot(u - a, n);
  const double dv = la::dot(v - a, n);
  return u + (du / (du - dv)) * (v - u);
}

struct BuildArrangement {
  std::vector<std::array<vec3, 3>> tri;
  std::vector<std::array<int, 3>> vid;
  std::vector<vec3> faceN;  // la::cross(b-a,c-a), unnormalized outward
  std::vector<std::vector<BuildSeam>> faceSeams;
  std::vector<char> seamed;
  bool boundaryTouch = false;
  bool ok = true;  // false = a structural anomaly (fail closed)
};

// Enumerate the self-crossing arrangement AND record each seam's canonical 3D
// segment per incident face (the geometry EnumerateSelfCrossings only counted).
// `face2cluster` (from DetectCoplanarClusters) marks exactly-coplanar face
// groups: same-cluster pairs are SKIPPED here (the in-plane fold resolves them,
// not the transversal seam machinery), so their exact-zero pierce ties do not
// raise A.boundaryTouch.  A cross-cluster / non-coplanar exact-zero tie still
// sets boundaryTouch (the SoS residue).
BuildArrangement RecordSeams(const Manifold::Impl& in,
                             const std::vector<int>& face2cluster) {
  BuildArrangement A;
  const int nTri = static_cast<int>(in.NumTri());
  A.tri.resize(nTri);
  A.vid.resize(nTri);
  A.faceN.resize(nTri);
  A.faceSeams.resize(nTri);
  A.seamed.assign(nTri, 0);
  std::vector<vec3> lo(nTri), hi(nTri);
  for (int t = 0; t < nTri; ++t) {
    for (int k = 0; k < 3; ++k) {
      A.vid[t][k] = in.halfedge_.Start(3 * t + k);
      A.tri[t][k] = in.vertPos_[A.vid[t][k]];
    }
    A.faceN[t] =
        la::cross(A.tri[t][1] - A.tri[t][0], A.tri[t][2] - A.tri[t][0]);
    lo[t] = la::min(la::min(A.tri[t][0], A.tri[t][1]), A.tri[t][2]);
    hi[t] = la::max(la::max(A.tri[t][0], A.tri[t][1]), A.tri[t][2]);
  }
  // Once-only pierce cache (doc R1): the UNDIRECTED (min,max) edge key makes a
  // shared mesh edge traversed in opposite order by two adjacent faces resolve
  // to ONE bit-identical pierce point, so chaining seams meet exactly and the
  // cross-face weld cannot manufacture a twin.  This is a CORRECTNESS-by-
  // construction backstop for the near-parallel tail (R1), not a runtime check:
  // it is verified NON-load-bearing on the shipped corpus (s2-verify) AND on a
  // broad general-position sphere family (reg3d-s3: bit-identical emitted
  // volume under a directed-key mutation, the eps weld absorbing the sub-ULP
  // divergence), so there is no DEBUG_ASSERT to demote it to - it stays as the
  // R1 backstop.
  std::map<PierceKey, vec3> cache;
  auto pierce = [&](int edgeV0, int edgeV1, int piercedTri) -> vec3 {
    const PierceKey key{std::min(edgeV0, edgeV1), std::max(edgeV0, edgeV1),
                        piercedTri};
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    const vec3 p = SegPlanePoint(in.vertPos_[edgeV0], in.vertPos_[edgeV1],
                                 A.tri[piercedTri][0], A.tri[piercedTri][1],
                                 A.tri[piercedTri][2]);
    cache.emplace(key, p);
    return p;
  };
  auto bboxOverlap = [&](int i, int j) {
    return !(hi[i].x < lo[j].x || hi[j].x < lo[i].x || hi[i].y < lo[j].y ||
             hi[j].y < lo[i].y || hi[i].z < lo[j].z || hi[j].z < lo[i].z);
  };
  auto sharesVert = [&](int i, int j) {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (A.vid[i][a] == A.vid[j][b]) return true;
    return false;
  };
  // Vertices used by each coplanar cluster (all lie on that cluster's plane):
  // an edge touching a folded plane at one of ITS OWN cluster vertices is a
  // riser vertex, not a transversal vertex-on-face SoS tie.
  int nClusters = 0;
  for (int c : face2cluster) nClusters = std::max(nClusters, c + 1);
  std::vector<std::set<int>> clusterVerts(nClusters);
  for (int f = 0; f < nTri; ++f)
    if (face2cluster[f] >= 0)
      for (int k = 0; k < 3; ++k)
        clusterVerts[face2cluster[f]].insert(A.vid[f][k]);
  for (int i = 0; i < nTri; ++i) {
    for (int j = i + 1; j < nTri; ++j) {
      if (!bboxOverlap(i, j)) continue;
      if (sharesVert(i, j)) continue;
      // Same exactly-coplanar cluster: the in-plane fold owns this pair; its
      // exact-zero pierce ties are not a seam and not an SoS boundary-touch.
      if (face2cluster[i] >= 0 && face2cluster[i] == face2cluster[j]) continue;
      const auto& T0 = A.tri[i];
      const auto& T1 = A.tri[j];
      // Collect the up-to-two seam endpoints: i's edges piercing tri j (keyed
      // to plane j) and j's edges piercing tri i (keyed to plane i).  Every
      // endpoint lands on the plane-i/\plane-j intersection line, so it is
      // exact in both faces' bases.
      std::array<vec3, 4> pts;
      std::array<int, 4> ptTri;  // piercedTri per endpoint (interiority source)
      int nPts = 0;
      bool boundary = false;
      // An edge of `owner` that only TOUCHES `tgt`'s plane (does not cross it)
      // is a boundary contact, not a transversal crossing, so it adds no
      // winding jump.  Two flavors are benign:
      //  - EDGE-IN-PLANE (both endpoints on the plane): in a valid 2-manifold
      //  the
      //    edge lies outside tgt's triangle or on a shared (skipped) boundary,
      //    so benign when tgt is unfolded; when tgt is a folded cluster face it
      //    is benign only if the edge is that cluster's own boundary (opposite
      //    face a cluster member) - a wall rising off the fold, already a
      //    constraint.
      //  - CLUSTER-VERTEX-ON-PLANE (one endpoint on the plane, and that
      //  endpoint
      //    is one of tgt's cluster's OWN vertices): a riser vertex of the fold,
      //    not the transversal vertex-on-face SoS tie.
      // Anything else (a non-cluster vertex on a face, or an edge-edge
      // crossing) stays the single-global-SoS residue.
      auto benignInPlane = [&](int owner, int e, int tgt) {
        const int e1 = (e + 1) % 3;
        const int su = Orient3DFilterSign(A.tri[tgt][0], A.tri[tgt][1],
                                          A.tri[tgt][2], A.tri[owner][e]);
        const int sv = Orient3DFilterSign(A.tri[tgt][0], A.tri[tgt][1],
                                          A.tri[tgt][2], A.tri[owner][e1]);
        if (su == 0 && sv == 0) {  // edge-in-plane
          if (face2cluster[tgt] < 0) return true;
          const int pairFace =
              static_cast<int>(in.halfedge_.Pair(3 * owner + e)) / 3;
          return face2cluster[pairFace] == face2cluster[tgt];
        }
        if ((su == 0) != (sv == 0)) {  // vertex-on-plane
          const int onV = (su == 0) ? A.vid[owner][e] : A.vid[owner][e1];
          // A cluster vertex on its own folded plane is a fold-owned riser.
          if (face2cluster[tgt] >= 0 &&
              clusterVerts[face2cluster[tgt]].count(onV) > 0)
            return true;
          // A vertex on tgt's plane but strictly OUTSIDE tgt's triangle does
          // not touch tgt's face - benign (a valid-manifold corner grazing an
          // adjacent face's plane).  Only a vertex inside / on tgt's triangle
          // is the vertex-on-face SoS tie.
          const vec3& p =
              A.vid[owner][e] == onV ? A.tri[owner][e] : A.tri[owner][e1];
          const vec3 n = la::cross(A.tri[tgt][1] - A.tri[tgt][0],
                                   A.tri[tgt][2] - A.tri[tgt][0]);
          const double area2 = la::length2(n);
          if (!(area2 > 0.0)) return false;
          const double margin = area2 * 1e-9;
          const double s0 = la::dot(
              n, la::cross(A.tri[tgt][1] - A.tri[tgt][0], p - A.tri[tgt][0]));
          const double s1 = la::dot(
              n, la::cross(A.tri[tgt][2] - A.tri[tgt][1], p - A.tri[tgt][1]));
          const double s2 = la::dot(
              n, la::cross(A.tri[tgt][0] - A.tri[tgt][2], p - A.tri[tgt][2]));
          return s0 < -margin || s1 < -margin || s2 < -margin;
        }
        return false;  // edge-edge crossing tie: real
      };
      // An edge whose BOTH endpoints are vertices of one coplanar cluster lies
      // in that cluster's folded plane; its exact-zero grazes are coplanar
      // in-plane incidences the fold owns (cluster-face edges and the cap edges
      // of walls rising off the fold).  A genuine transversal cross of a
      // cluster face is a PIERCE (r==1) recorded as a seam and caught as
      // entanglement; only a graze OUTSIDE every cluster plane is the SoS
      // residue.
      auto edgeInClusterPlane = [&](int owner, int e) {
        const int a = A.vid[owner][e], b = A.vid[owner][(e + 1) % 3];
        for (int c = 0; c < nClusters; ++c)
          if (clusterVerts[c].count(a) && clusterVerts[c].count(b)) return true;
        return false;
      };
      for (int e = 0; e < 3; ++e) {
        const int r =
            EdgePiercesTri(T0[e], T0[(e + 1) % 3], T1[0], T1[1], T1[2]);
        if (r == 1 && nPts < 4) {
          ptTri[nPts] = j;  // pierces tri j -> on i's edge, interior to j
          pts[nPts++] = pierce(A.vid[i][e], A.vid[i][(e + 1) % 3], j);
        } else if (r == -1 && !edgeInClusterPlane(i, e) &&
                   !benignInPlane(i, e, j))
          boundary = true;
      }
      for (int e = 0; e < 3; ++e) {
        const int r =
            EdgePiercesTri(T1[e], T1[(e + 1) % 3], T0[0], T0[1], T0[2]);
        if (r == 1 && nPts < 4) {
          ptTri[nPts] = i;  // pierces tri i -> on j's edge, interior to i
          pts[nPts++] = pierce(A.vid[j][e], A.vid[j][(e + 1) % 3], i);
        } else if (r == -1 && !edgeInClusterPlane(j, e) &&
                   !benignInPlane(j, e, i))
          boundary = true;
      }
      if (boundary) {
        // A COPLANAR exact-zero tie is never a transversal crossing: the two
        // faces share a plane, so this is either a folded cluster pair (skipped
        // above) or a benign non-overlapping coplanar contact.  Only a
        // NON-coplanar exact-zero tie (a vertex-on-face / edge-on-edge
        // incidence between transversal faces) is the single-global-SoS
        // residue.
        bool coplanar = true;
        for (int k = 0; k < 3 && coplanar; ++k)
          if (Orient3DFilterSign(T0[0], T0[1], T0[2], T1[k]) != 0)
            coplanar = false;
        for (int k = 0; k < 3 && coplanar; ++k)
          if (Orient3DFilterSign(T1[0], T1[1], T1[2], T0[k]) != 0)
            coplanar = false;
        if (!coplanar) A.boundaryTouch = true;
        continue;
      }
      if (nPts == 0) continue;  // no genuine crossing
      if (nPts != 2) {
        // A genuine seam has exactly two endpoints; anything else is a
        // degenerate incidence the level-0 filter did not flag - fail closed.
        A.ok = false;
        continue;
      }
      // Interior-to-face flag: an endpoint is interior to face i iff it was
      // built by piercing tri i (ptTri==i), else it lies on face i's edge.
      A.faceSeams[i].push_back(
          {pts[0], pts[1], A.faceN[j], ptTri[0] == i, ptTri[1] == i});
      A.faceSeams[j].push_back(
          {pts[0], pts[1], A.faceN[i], ptTri[0] == j, ptTri[1] == j});
      A.seamed[i] = 1;
      A.seamed[j] = 1;
    }
  }
  return A;
}

// Extract the bounded CELLS of a planar subdivision: undirected `uedges` over
// `pts` (no interior crossings on this corpus).  Dangling spurs (degree-1
// chains) bound no cell and are pruned first (so triangulation sees no
// zero-area spike).  Returns one CCW vertex-index loop per bounded cell (the
// unbounded outer face is dropped by its negative signed area), or false on a
// malformed walk.  Standard halfedge face traversal: the outgoing half-edges at
// each vertex are angularly ordered and next(u->v) is the outgoing edge at v
// immediately CLOCKWISE from v->u, which keeps the cell interior on the left.
// `holes`, when non-null, receives the CW (negative-area) boundary walks: the
// unbounded outer face AND any interior hole loops (a face nested inside
// another, e.g. one coplanar triangle contained in another).  The seamed-face
// path passes nullptr (its cells are simply connected).  The coplanar fold uses
// them to triangulate multiply-connected cells so a contained region is not
// double-covered.
bool ExtractCells(const std::vector<vec2>& pts,
                  const std::vector<std::pair<int, int>>& uedges,
                  std::vector<std::vector<int>>& cells,
                  std::vector<std::vector<int>>* holes = nullptr) {
  const int n = static_cast<int>(pts.size());
  std::vector<std::set<int>> nbr(n);
  for (const auto& e : uedges) {
    if (e.first == e.second) continue;
    nbr[e.first].insert(e.second);
    nbr[e.second].insert(e.first);
  }
  for (bool changed = true; changed;) {  // prune degree-1 spurs
    changed = false;
    for (int v = 0; v < n; ++v)
      if (nbr[v].size() == 1) {
        const int w = *nbr[v].begin();
        nbr[v].clear();
        nbr[w].erase(v);
        changed = true;
      }
  }
  std::vector<std::vector<int>> order(n);  // CCW-sorted neighbors
  std::vector<std::unordered_map<int, int>> at(n);
  for (int v = 0; v < n; ++v) {
    order[v].assign(nbr[v].begin(), nbr[v].end());
    std::sort(order[v].begin(), order[v].end(), [&](int p, int q) {
      return std::atan2(pts[p].y - pts[v].y, pts[p].x - pts[v].x) <
             std::atan2(pts[q].y - pts[v].y, pts[q].x - pts[v].x);
    });
    for (int k = 0; k < static_cast<int>(order[v].size()); ++k)
      at[v][order[v][k]] = k;
  }
  std::set<std::pair<int, int>> visited;
  for (int s0 = 0; s0 < n; ++s0)
    for (const int s1 : order[s0]) {
      if (visited.count({s0, s1})) continue;
      std::vector<int> loop;
      int u = s0, v = s1;
      bool bad = false;
      do {
        visited.insert({u, v});
        loop.push_back(u);
        const auto it = at[v].find(u);
        if (it == at[v].end()) {
          bad = true;
          break;
        }
        const int deg = static_cast<int>(order[v].size());
        const int w = order[v][(it->second - 1 + deg) % deg];
        u = v;
        v = w;
      } while (!(u == s0 && v == s1) &&
               loop.size() <= static_cast<size_t>(2 * uedges.size() + 4));
      if (bad || !(u == s0 && v == s1)) return false;
      double area = 0.0;
      for (size_t k = 0; k < loop.size(); ++k) {
        const vec2 p = pts[loop[k]], q = pts[loop[(k + 1) % loop.size()]];
        area += p.x * q.y - q.x * p.y;
      }
      if (area > 0.0)
        cells.push_back(std::move(loop));
      else if (holes && loop.size() >= 3)
        holes->push_back(std::move(loop));
    }
  return true;
}

// Emit the retained sub-faces of one SEAMED face into `out` (3D triangles at
// canonical positions), or set ok=false to fail closed.  Reuse RemoveOverlaps2D
// as the ARRANGEMENT primitive (robust crossing-split + eps-merge via
// edgeSubdiv), extract the cells, and classify EACH by the real 3D coupled
// winding: a cell is retained iff w_S on the +n_f side == 0 (the {w_S>=1}
// boundary criterion, uniform and flip-free for oriented mult-1 faces).
void EmitSeamedFace(std::vector<OutTri3D>& out, const BuildArrangement& A,
                    int f, const Manifold::Impl& in,
                    const std::vector<vec3>& seeds, double eps, bool& ok) {
  const vec3 a = A.tri[f][0], b = A.tri[f][1], c = A.tri[f][2];
  const double nLen = la::length(A.faceN[f]);
  if (!(nLen > 0.0)) {
    ok = false;
    return;
  }
  const vec3 nHat = A.faceN[f] / nLen;
  const vec3 e1raw = b - a;
  const double e1Len = la::length(e1raw);
  if (!(e1Len > 0.0)) {
    ok = false;
    return;
  }
  const vec3 e1 = e1raw / e1Len;
  const vec3 e2 = la::cross(nHat, e1);  // e1 x e2 == nHat: 2D-CCW -> +nHat
  auto proj = [&](const vec3& P) {
    return vec2(la::dot(P - a, e1), la::dot(P - a, e2));
  };

  // Vertices, deduped by canonical 3D bit pattern so shared endpoints (chain
  // junctions, corners) collapse to one input vertex with one projection.
  std::vector<vec2> verts2;
  std::vector<vec3> canon3;
  std::map<std::tuple<double, double, double>, int> vidx;
  auto getV = [&](const vec3& P) {
    const std::tuple<double, double, double> key{P.x, P.y, P.z};
    auto it = vidx.find(key);
    if (it != vidx.end()) return it->second;
    const int id = static_cast<int>(verts2.size());
    verts2.push_back(proj(P));
    canon3.push_back(P);
    vidx.emplace(key, id);
    return id;
  };
  const int iA = getV(a), iB = getV(b), iC = getV(c);
  const vec2 pa = verts2[iA], pb = verts2[iB], pc = verts2[iC];
  if (!(0.5 * ((pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x)) >
        0.0)) {
    ok = false;  // left-handed basis or degenerate projection
    return;
  }
  std::vector<EdgeM> edges = {{iA, iB, 1}, {iB, iC, 1}, {iC, iA, 1}};
  for (const BuildSeam& s : A.faceSeams[f]) {
    const int v0 = getV(s.p0), v1 = getV(s.p1);
    if (v0 != v1) edges.push_back({v0, v1, 1});
  }

  // ARRANGEMENT via RemoveOverlaps2D: edgeSubdiv gives the exact per-input-edge
  // subdivision (triangle edges split at the on-edge seam endpoints; seams
  // split at any crossing).  The winding rule is irrelevant here (we take only
  // the subdivision), so any pred works.  eps is the component weld radius; the
  // canonical seam points are exact-shared, so no extra construction headroom
  // is needed and a tight eps avoids merging genuinely-distinct sub-eps
  // features.
  std::vector<std::vector<vec2>> sub;
  RemoveOverlaps2D(verts2, edges, eps, /*debug=*/false, WindRule::Add,
                   /*trace=*/nullptr, /*edgesNeg=*/nullptr, &sub);

  // Reconstruct the planar subdivision, mapping each arrangement position back
  // to its input vertex bit-exactly.  A position with no input preimage is a
  // NEW crossing vertex (a >2-sheet triple point, unbuilt) -> fail closed.
  std::map<std::tuple<double, double>, int> pos2in;
  for (int k = 0; k < static_cast<int>(verts2.size()); ++k)
    pos2in.emplace(std::tuple<double, double>{verts2[k].x, verts2[k].y}, k);
  std::vector<std::pair<int, int>> uedges;
  for (const auto& poly : sub) {
    for (size_t k = 0; k + 1 < poly.size(); ++k) {
      const auto i0 = pos2in.find({poly[k].x, poly[k].y});
      const auto i1 = pos2in.find({poly[k + 1].x, poly[k + 1].y});
      if (i0 == pos2in.end() || i1 == pos2in.end()) {
        ok = false;
        return;
      }
      uedges.push_back({i0->second, i1->second});
    }
  }

  std::vector<std::vector<int>> cells;
  if (!ExtractCells(verts2, uedges, cells)) {
    ok = false;
    return;
  }
  for (const std::vector<int>& cell : cells) {
    if (cell.size() < 3) continue;
    // Triangulate the cell; the largest triangle's centroid is a guaranteed
    // interior classify point.
    PolygonsIdx pidx(1);
    for (int idx : cell) pidx[0].push_back({verts2[idx], idx});
    const std::vector<ivec3> tris = TriangulateIdx(pidx, eps);
    if (tris.empty()) continue;
    int best = 0;
    double bestArea = -1.0;
    for (int t = 0; t < static_cast<int>(tris.size()); ++t) {
      const vec2 p0 = verts2[tris[t].x], p1 = verts2[tris[t].y],
                 p2 = verts2[tris[t].z];
      const double ar = std::abs((p1.x - p0.x) * (p2.y - p0.y) -
                                 (p1.y - p0.y) * (p2.x - p0.x));
      if (ar > bestArea) {
        bestArea = ar;
        best = t;
      }
    }
    const vec2 cen2 =
        (verts2[tris[best].x] + verts2[tris[best].y] + verts2[tris[best].z]) /
        3.0;
    const vec3 cen3 = a + cen2.x * e1 + cen2.y * e2;
    const std::optional<int> g = RobustWinding(in, cen3 + eps * nHat, seeds);
    if (!g) {  // filter-uncertain deciding predicate (SoS axis): fail closed
      ok = false;
      return;
    }
    // Witness theorem, general form: for an oriented mult-1 face w_below =
    // w_above + 1 universally (the +/-1 crossing delta, seed-independent), so a
    // cell is on d{w_S>=1} iff EXACTLY ONE side has w>=1, which for mult-1 is
    // exactly w_above == 0 (w_above<=0 && w_above+1>=1).  A negative w_above
    // means w_below = w_above+1 <= 0: BOTH sides exterior to {w_S>=1}, so the
    // cell is dropped, NOT failed closed (openscad/subtraction absorption).
    // The solid, when retained, is always on the -n_f side -> original
    // orientation, flip-free.
    if (*g != 0)
      continue;  // w_above != 0: buried (>=1) or exterior (<=0): drop
    // Retained: emit at canonical 3D, 2D-CCW -> +nHat = original orientation.
    for (const ivec3& t : tris)
      out.push_back({canon3[t.x], canon3[t.y], canon3[t.z]});
  }
}

// 2D point-in-triangle (inclusive), orientation-agnostic: true iff p is on the
// same side (or on) all three directed edges under either winding.
bool PointInTri2D(const vec2& p, const vec2& a, const vec2& b, const vec2& c) {
  auto cr = [](const vec2& u, const vec2& v) { return u.x * v.y - u.y * v.x; };
  const double d1 = cr(b - a, p - a), d2 = cr(c - b, p - b),
               d3 = cr(a - c, p - c);
  const bool neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
  const bool pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
  return !(neg && pos);
}

// EXACT-COPLANAR IN-PLANE FOLD (docs/Regularize3D.md coplanar axis).  Each
// cluster of exactly-coplanar faces (DetectCoplanarClusters) is overlaid in its
// shared plane: RemoveOverlaps2D arranges the members' triangle boundaries (an
// arrangement primitive - coplanar edge crossings ARE new vertices, all exact
// in-plane reconstructions, unlike a transversal triple point), ExtractCells
// yields the non-overlapping sub-faces, and each sub-face carries an integer
// MULT m = the net signed in-plane cover (anti-oriented content cancels, same-
// oriented sums).  Retention generalizes B's mult-1 rule: with w_below =
// w_above + m (the 3D winding jump across the plane equals the coincident
// cover, so it is SELF-CHECKED against the real coupled winding on both sides),
// a sub-face is on d{w_S>=1} iff EXACTLY ONE side is inside {w>=1}; the solid
// side fixes the emitted orientation.  m==0 (pure cancellation) drops.  A
// cluster face that is ALSO transversally seamed is the coplanar/transversal
// ENTANGLEMENT (a seam line would straddle a cell); that stays fail-closed with
// its own named reason, distinct from the near-coplanar residue.  Any
// degenerate projection, malformed cell walk, filter-uncertain winding, or
// self-check mismatch fails closed (never a silent wrong resolve).
void FoldCoplanarClusters(std::vector<OutTri3D>& out, const Manifold::Impl& in,
                          const BuildArrangement& A,
                          const std::vector<int>& face2cluster,
                          const std::vector<vec3>& seeds, double eps,
                          bool& ok) {
  int nc = 0;
  for (int c : face2cluster) nc = std::max(nc, c + 1);
  if (nc == 0) return;
  std::vector<std::vector<int>> clusters(nc);
  for (int f = 0; f < static_cast<int>(face2cluster.size()); ++f)
    if (face2cluster[f] >= 0) clusters[face2cluster[f]].push_back(f);

  for (const std::vector<int>& faces : clusters) {
    // ENTANGLEMENT: a cluster face pierced transversally by a non-coplanar face
    // (a seam that would split a fold cell). Fail closed with a distinct
    // reason.
    for (int f : faces)
      if (A.seamed[f]) {
        ok = false;
        return;
      }
    const int f0 = faces[0];
    const double nLen = la::length(A.faceN[f0]);
    if (!(nLen > 0.0)) {
      ok = false;
      return;
    }
    const vec3 nHat = A.faceN[f0] / nLen;
    const vec3 a0 = A.tri[f0][0];
    const vec3 e1raw = A.tri[f0][1] - a0;
    const double e1Len = la::length(e1raw);
    if (!(e1Len > 0.0)) {
      ok = false;
      return;
    }
    const vec3 e1 = e1raw / e1Len;
    const vec3 e2 = la::cross(nHat, e1);
    auto proj = [&](const vec3& P) {
      return vec2(la::dot(P - a0, e1), la::dot(P - a0, e2));
    };

    // Input verts (dedup by canonical 3D bit pattern) + triangle-boundary
    // edges; each member triangle carries its signed orientation vs nHat.
    std::vector<vec2> verts2;
    std::vector<vec3> canon3;
    std::map<std::tuple<double, double, double>, int> vidx;
    auto getV = [&](const vec3& P) {
      const std::tuple<double, double, double> key{P.x, P.y, P.z};
      auto it = vidx.find(key);
      if (it != vidx.end()) return it->second;
      const int id = static_cast<int>(verts2.size());
      verts2.push_back(proj(P));
      canon3.push_back(P);
      vidx.emplace(key, id);
      return id;
    };
    struct FaceTri {
      vec2 p0, p1, p2;
      int s;
    };
    std::vector<FaceTri> ftris;
    // Directed triangle-boundary edges; an edge shared by two triangles of the
    // SAME member facet (a quad's diagonal) appears in both directions and is
    // NOT a real arrangement boundary - cancel it, else it would split the fold
    // caps (crossing another face's outline) without splitting that face's
    // walls, manufacturing a T-junction.
    std::map<std::pair<int, int>, int> dir;  // (min,max) -> net signed count
    for (int f : faces) {
      const int i0 = getV(A.tri[f][0]), i1 = getV(A.tri[f][1]),
                i2 = getV(A.tri[f][2]);
      const int tv[3] = {i0, i1, i2};
      for (int k = 0; k < 3; ++k) {
        int a = tv[k], b = tv[(k + 1) % 3];
        dir[{std::min(a, b), std::max(a, b)}] += (a < b) ? 1 : -1;
      }
      const int s = (la::dot(A.faceN[f], nHat) > 0.0) ? 1 : -1;
      ftris.push_back({verts2[i0], verts2[i1], verts2[i2], s});
    }
    std::vector<EdgeM> segEdges;
    for (const auto& [e, net] : dir)
      if (net != 0)
        segEdges.push_back({e.first, e.second, 1});  // boundary only

    // ARRANGEMENT via RemoveOverlaps2D (the same primitive the seamed path
    // reuses at EmitSeamedFace): the member triangles' boundary edges overlap
    // in-plane, and edgeSubdiv splits each input edge at BOTH proper crossings
    // AND T-junctions / collinear incidences (a hand-rolled pairwise-crossing
    // arrangement would miss the T-junctions, splitting one edge without
    // splitting the edge that ends on it - RO2DProbe, reg3d-s4-verify Audit
    // 2a). The winding rule is irrelevant here (we consume only the
    // subdivision). Every arrangement vertex lies in this exact plane, so a NEW
    // crossing position's 3D image is a0 + x*e1 + y*e2 (an input vert keeps its
    // canonical 3D via getP; the eps-box match folds a merged endpoint back
    // onto its input vert).
    std::vector<vec2> pts = verts2;
    std::vector<vec3> pts3 = canon3;
    auto getP = [&](const vec2& q) {
      for (int k = 0; k < static_cast<int>(pts.size()); ++k)
        if (std::abs(pts[k].x - q.x) <= eps && std::abs(pts[k].y - q.y) <= eps)
          return k;
      const int id = static_cast<int>(pts.size());
      pts.push_back(q);
      pts3.push_back(a0 + q.x * e1 + q.y * e2);
      return id;
    };
    std::vector<std::vector<vec2>> edgeSubdiv;
    RemoveOverlaps2D(verts2, segEdges, eps, /*debug=*/false, WindRule::Add,
                     /*trace=*/nullptr, /*edgesNeg=*/nullptr, &edgeSubdiv);
    std::vector<std::pair<int, int>> uedges;
    for (const std::vector<vec2>& poly : edgeSubdiv)
      for (size_t k = 0; k + 1 < poly.size(); ++k) {
        const int i0 = getP(poly[k]), i1 = getP(poly[k + 1]);
        if (i0 != i1) uedges.push_back({i0, i1});
      }

    std::vector<std::vector<int>> cells, holeLoops;
    if (!ExtractCells(pts, uedges, cells, &holeLoops)) {
      ok = false;
      return;
    }
    // Signed 2D area (CCW>0) helper, and point-strictly-inside-loop test, to
    // attach each interior hole loop to the smallest cell that contains it (a
    // coplanar face fully inside another becomes a hole, not a separate cell).
    auto loopArea = [&](const std::vector<int>& L) {
      double ar = 0.0;
      for (size_t k = 0; k < L.size(); ++k) {
        const vec2 p = pts[L[k]], q = pts[L[(k + 1) % L.size()]];
        ar += p.x * q.y - q.x * p.y;
      }
      return 0.5 * ar;
    };
    auto inLoop = [&](const vec2& p, const std::vector<int>& L) {
      bool in = false;
      for (size_t k = 0, j = L.size() - 1; k < L.size(); j = k++) {
        const vec2 a = pts[L[k]], b = pts[L[j]];
        if (((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
          in = !in;
      }
      return in;
    };
    std::vector<double> cellArea(cells.size());
    for (size_t c = 0; c < cells.size(); ++c) cellArea[c] = loopArea(cells[c]);
    std::vector<std::vector<int>> cellHoles(cells.size());  // hole loop indices
    for (int h = 0; h < static_cast<int>(holeLoops.size()); ++h) {
      // A CW loop is the same edge cycle as the CCW cell on its other side, so
      // attach it to the smallest cell that STRICTLY contains it (larger area)
      // - never to the identical cell across the same edges.  Use the loop
      // centroid for a robust interior probe.
      vec2 hc(0.0, 0.0);
      for (int idx : holeLoops[h]) hc += pts[idx];
      hc /= static_cast<double>(holeLoops[h].size());
      const double hArea = std::abs(loopArea(holeLoops[h]));
      int best = -1;
      double bestA = 0.0;
      for (size_t c = 0; c < cells.size(); ++c)
        if (cellArea[c] > hArea * (1.0 + 1e-9) && inLoop(hc, cells[c]) &&
            (best < 0 || cellArea[c] < bestA)) {
          best = static_cast<int>(c);
          bestA = cellArea[c];
        }
      if (best >= 0) cellHoles[best].push_back(h);  // else the unbounded face
    }
    for (size_t ci = 0; ci < cells.size(); ++ci) {
      const std::vector<int>& cell = cells[ci];
      if (cell.size() < 3) continue;
      PolygonsIdx pgon(1 + cellHoles[ci].size());
      for (int idx : cell) pgon[0].push_back({pts[idx], idx});
      for (size_t hh = 0; hh < cellHoles[ci].size(); ++hh)
        for (int idx : holeLoops[cellHoles[ci][hh]])
          pgon[hh + 1].push_back({pts[idx], idx});
      // A malformed cell (degenerate sliver, non-simple loop) makes Triangulate
      // throw; that is a decline, not a crash.
      std::vector<ivec3> tris;
      try {
        tris = TriangulateIdx(pgon, eps);
      } catch (...) {
        ok = false;
        return;
      }
      if (tris.empty()) continue;
      int best = 0;
      double bestArea = -1.0;
      for (int t = 0; t < static_cast<int>(tris.size()); ++t) {
        const vec2 p0 = pts[tris[t].x], p1 = pts[tris[t].y],
                   p2 = pts[tris[t].z];
        const double ar = std::abs((p1.x - p0.x) * (p2.y - p0.y) -
                                   (p1.y - p0.y) * (p2.x - p0.x));
        if (ar > bestArea) {
          bestArea = ar;
          best = t;
        }
      }
      // Classify at the largest sub-triangle's INCENTER (strictly interior, as
      // far from every edge as possible), then jitter by a fraction of the
      // incircle radius in a GENERIC direction so the point does not land on a
      // member wall's plane (e.g. an axis-aligned overlap boundary), where the
      // winding orient3d would be near-degenerate and fail closed spuriously.
      const vec2 va = pts[tris[best].x], vb = pts[tris[best].y],
                 vc = pts[tris[best].z];
      const double la0 = la::length(vc - vb), lb = la::length(va - vc),
                   lc = la::length(vb - va);
      const double lsum = la0 + lb + lc;
      vec2 cen2 = lsum > 0.0 ? (la0 * va + lb * vb + lc * vc) / lsum
                             : (va + vb + vc) / 3.0;
      const double inrad = lsum > 0.0 ? bestArea / lsum : 0.0;  // 2*area/perim
      cen2 += 0.25 * inrad * la::normalize(vec2(0.4359, 0.9000));
      const vec3 cen3 = a0 + cen2.x * e1 + cen2.y * e2;
      int m = 0;
      for (const FaceTri& ft : ftris)
        if (PointInTri2D(cen2, ft.p0, ft.p1, ft.p2)) m += ft.s;
      if (m == 0) continue;  // net cover cancels: emit nothing
      const std::optional<int> wa = RobustWinding(in, cen3 + eps * nHat, seeds);
      const std::optional<int> wb = RobustWinding(in, cen3 - eps * nHat, seeds);
      if (!wa || !wb) {  // filter-uncertain classify (SoS): fail closed
        ok = false;
        return;
      }
      if (*wb - *wa != m) {  // 3D winding jump must equal the in-plane cover
        ok = false;
        return;
      }
      const bool aboveIn = *wa >= 1, belowIn = *wb >= 1;
      if (aboveIn == belowIn)
        continue;  // both sides same class: not a boundary
      // Retained. Solid on the {w>=1} side fixes orientation: solid on the
      // -nHat side (belowIn) keeps the CCW 2D winding (+nHat); solid on +nHat
      // reverses.
      for (const ivec3& t : tris) {
        if (belowIn)
          out.push_back({pts3[t.x], pts3[t.y], pts3[t.z]});
        else
          out.push_back({pts3[t.x], pts3[t.z], pts3[t.y]});
      }
    }
  }
}

// Classify + emit the CLEAN (un-seamed, un-folded) faces.  g_above is constant
// across any mesh edge that carries no seam transition, and clean-clean edges
// never do, so clean faces partition into patches of uniform coverage; one
// winding probe per patch decides keep-whole vs drop by the SAME witness rule
// as the seamed path (keep iff w_above == 0).  A NEGATIVE w_above is exterior
// on both sides and DROPS (the axis-1 subtraction absorption), not a
// fail-closed; only a filter-uncertain probe (the SoS axis) fails closed.
bool EmitCleanFaces(std::vector<OutTri3D>& out, const Manifold::Impl& in,
                    const BuildArrangement& A,
                    const std::vector<int>& face2cluster,
                    const std::vector<vec3>& seeds, double eps) {
  const int nTri = static_cast<int>(in.NumTri());
  // A face is "clean" only if it is neither seamed nor part of a coplanar
  // cluster (the fold owns cluster faces); cluster faces bound the flood so a
  // patch never crosses into folded content.
  auto isClean = [&](int t) { return !A.seamed[t] && face2cluster[t] < 0; };
  DisjointSets patches(nTri);
  for (int h = 0; h < static_cast<int>(in.halfedge_.size()); ++h) {
    const int t = h / 3, u = in.halfedge_.Pair(h) / 3;
    if (isClean(t) && isClean(u)) patches.unite(t, u);
  }
  std::unordered_map<int, int> patchKeep;  // root -> 0 drop, 1 keep
  for (int t = 0; t < nTri; ++t) {
    if (!isClean(t)) continue;
    const int root = static_cast<int>(patches.find(t));
    auto it = patchKeep.find(root);
    int keep;
    if (it == patchKeep.end()) {
      const vec3 cen = (A.tri[t][0] + A.tri[t][1] + A.tri[t][2]) / 3.0;
      const double nLen = la::length(A.faceN[t]);
      if (!(nLen > 0.0)) return false;
      const vec3 nHat = A.faceN[t] / nLen;
      const std::optional<int> g = RobustWinding(in, cen + eps * nHat, seeds);
      if (!g)
        return false;  // filter-uncertain deciding predicate (SoS): closed
      // Same witness rule as the seamed path: keep iff w_above == 0.  A
      // negative w_above (subtraction / openscad-class) is exterior on both
      // sides -> drop (keep=0), NOT a fail-closed.
      keep = (*g == 0) ? 1 : 0;
      patchKeep.emplace(root, keep);
    } else {
      keep = it->second;
    }
    if (keep == 1) out.push_back({A.tri[t][0], A.tri[t][1], A.tri[t][2]});
  }
  return true;
}

// THE BUILD driver: fold exactly-coplanar clusters in-plane, emit seamed
// sub-faces + clean faces, assemble + weld.  Returns the regularized Impl, or a
// fatal.
StageResult<Manifold::Impl> RunCandidateBBuild(
    const Manifold::Impl& in, const BuildArrangement& A,
    const std::vector<int>& face2cluster, double eps) {
  // Winding seeds: a few far points in unrelated directions off the bbox.
  const vec3 c = in.bBox_.Center();
  const double L = in.bBox_.Scale() + 1.0;
  const std::vector<vec3> seeds = {
      c + L * vec3(3.13, 5.71, 1.37),   c + L * vec3(-2.71, 1.41, 4.19),
      c + L * vec3(1.73, -3.31, -2.23), c + L * vec3(-4.27, -1.19, 2.83),
      c + L * vec3(2.39, -4.61, 3.07),  c + L * vec3(-1.51, 3.89, -4.43)};

  std::vector<OutTri3D> emitted;
  bool ok = true;
  // Exact-coplanar clusters first (transversal seams on their faces are the
  // entanglement decline).
  FoldCoplanarClusters(emitted, in, A, face2cluster, seeds, eps, ok);
  if (!ok)
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "candidate B: exact-coplanar fold declined (coplanar/transversal "
        "entanglement, degenerate projection, or filter-uncertain classify) - "
        "fail-closed");
  const int nTri = static_cast<int>(in.NumTri());
  for (int f = 0; f < nTri && ok; ++f)
    if (A.seamed[f]) EmitSeamedFace(emitted, A, f, in, seeds, eps, ok);
  if (!ok)
    // B declined to build this face's arrangement exactly: a >2-sheet triple
    // point, a coplanar/degenerate projection, a malformed cell walk, or a
    // filter-uncertain classify probe (the SoS axis).  Negative winding is NOT
    // a decline - the witness rule absorbs it (w_above==0 retain).  Fail closed
    // - never emit geometry B could not verify.
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "candidate B: seam sub-face arrangement not exactly resolvable "
        "(triple point / degenerate / filter-uncertain) - fail-closed");
  if (!EmitCleanFaces(emitted, in, A, face2cluster, seeds, eps))
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "candidate B: clean-face winding probe was filter-uncertain (SoS) - "
        "fail-closed");
  return BuildImpl(emitted, eps);
}

// Candidate B (docs/Regularize3D.md "B's mechanism") - the dirty-core resolver.
// The validated MECHANISM (enumeration + coupled winding) is ported; THE BUILD
// (the {w_S>=1} halfedge boundary emission) reuses RemoveOverlaps2D per crossed
// face (RunCandidateBBuild above).  B enumerates + records the seam geometry,
// runs the build, and re-gates; anything it cannot resolve exactly (an
// exact-zero pierce tie = single-global SoS, a >2-sheet triple point, a
// coplanar seam, a negative-winding patch) FAILS CLOSED with a named reason -
// never a silent wrong result.  The caller re-gates the output once more
// (IsSelfIntersecting).
StageResult<Manifold::Impl> RunCandidateB(const Manifold::Impl& dirty,
                                          double eps) {
  // Exactly-coplanar face clusters are resolved by the in-plane fold; the
  // transversal seam enumeration skips their pairs so their exact-zero coplanar
  // ties do not raise boundaryTouch.
  const std::vector<int> face2cluster = DetectCoplanarClusters(dirty);
  const BuildArrangement A = RecordSeams(dirty, face2cluster);
  if (A.boundaryTouch) {
    // A deciding pierce predicate hit an exact-zero / filter-uncertain boundary
    // that the coplanar fold does NOT consume: a NON-coplanar vertex-on-face /
    // edge-in-face incidence (the residual single-global SoS tie family,
    // PokedCube-class), or the near-coplanar thin band.  Guessing a sign would
    // risk an oracle-wrong resolve.  Fail closed.
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "candidate B: non-coplanar exact-zero tie; single-global SoS "
        "(PokedCube-class) unbuilt - fail-closed");
  }
  if (!A.ok) {
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "candidate B: a self-crossing pair had a non-2-endpoint seam "
        "(degenerate incidence) - fail-closed");
  }
  return RunCandidateBBuild(dirty, A, face2cluster, eps);
}

// Compose the surviving components back into one Impl by CONCATENATION - no
// cross-component weld, no fusion (docs/Regularize3D.md step 6).  A single
// component is returned as-is (bitwise pass-through for a clean
// single-component input); multiple components are concatenated through
// MeshGL64, whose halfedge pairing is per-component (distinct-position verts
// across components never merge, so touching contacts stay separate).
Manifold::Impl ComposeComponents(std::vector<Manifold::Impl>& parts) {
  if (parts.empty()) return Manifold::Impl{};
  if (parts.size() == 1) return std::move(parts[0]);

  MeshGL64 combined;
  combined.numProp = 3;
  for (const Manifold::Impl& p : parts) {
    const MeshGL64 mg = GetMeshGLImpl<double, uint64_t>(p, -1);
    const uint64_t base = combined.NumVert();
    for (size_t i = 0; i < mg.vertProperties.size(); ++i)
      combined.vertProperties.push_back(mg.vertProperties[i]);
    for (size_t i = 0; i < mg.triVerts.size(); ++i)
      combined.triVerts.push_back(mg.triVerts[i] + base);
  }
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(combined);
}

}  // namespace

RegularizeResult RegularizeImpl(const Manifold::Impl& in, double eps) {
  RegularizeResult result;

  // Empty input -> empty output (matches RemoveOverlaps3D's trivially-empty
  // path).
  if (in.NumTri() == 0) {
    result.impl = Manifold::Impl{};
    return result;
  }

  // Resolve the machine-scale weld radius once; every component's gate uses it.
  if (eps <= 0.0) eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  if (eps <= 0.0 || !std::isfinite(eps)) {
    result.fatal = FatalReason::SubEpsInput;
    result.detail = "epsilon not computable";
    return result;
  }

  // 1. DECOMPOSE by connectivity.
  std::vector<Manifold::Impl> components = DecomposeComponents(in, eps);
  result.counters.components = static_cast<int>(components.size());

  // 2-5. Gate + dispatch every component.  We gate ALL components (the gate is
  // cheap) so the white-box dispatch counters are complete regardless of the
  // decompose order; the first fail-closed is the reported fatal, and no
  // partial output is composed once any component fails.
  std::optional<FatalReason> firstFatal;
  std::string firstDetail;
  std::vector<Manifold::Impl> outComponents;

  for (Manifold::Impl& comp : components) {
    const GateVerdict verdict = GateComponent(comp);
    if (verdict == GateVerdict::Invalid) {
      // Defensive: a component of a valid input is valid; a non-manifold one is
      // neither early-exitable nor a case B resolves.  Fail closed.
      ++result.counters.failClosed;
      if (!firstFatal) {
        firstFatal = FatalReason::NonManifoldEmission;
        firstDetail = "input component is not 2-manifold";
      }
      continue;
    }
    if (verdict == GateVerdict::Clean) {
      // 3. EARLY-EXIT: already the boundary of a simple solid.
      ++result.counters.clean;
      outComponents.push_back(std::move(comp));
      continue;
    }

    // 4. DIRTY -> candidate B.
    ++result.counters.dirty;
    StageResult<Manifold::Impl> bRes = RunCandidateB(comp, eps);
    if (!bRes.ok()) {
      ++result.counters.failClosed;
      if (!firstFatal) {
        firstFatal = bRes.fatal;
        firstDetail = std::move(bRes.detail);
      }
      continue;
    }
    // 5. RE-GATE B's output once (same gate as the input; B's coords are
    // double-rounded).  A clean pass composes in; a failure is the honest
    // fail-closed, never a silent wrong result.  This is a PRODUCTION
    // fail-closed backstop for the R1/R2 weld-fold blind spot (BuildImpl
    // already gates non-manifold emission; the re-gate's non-redundant job is
    // catching a MANIFOLD-but-self-intersecting output = a weld-manufactured
    // fold).  It is verified UNREACHED on constructible general-position
    // fixtures (reg3d-s3: 0/180 sphere variants produce re-gate-catchable
    // output - every bad case is caught earlier by BuildImpl's manifold gate),
    // i.e. it fires only in the unbuilt weld-fold regime.  It is deliberately
    // NOT demoted to a DEBUG_ASSERT: it must fail closed in RELEASE, not
    // compile out and admit wrong geometry.
    Manifold::Impl bImpl = std::move(*bRes.value);
    bImpl.epsilon_ = eps;
    if (GateComponent(bImpl) != GateVerdict::Clean) {
      ++result.counters.failClosed;
      if (!firstFatal) {
        firstFatal = FatalReason::NonManifoldEmission;
        firstDetail = "candidate B output failed the re-gate";
      }
      continue;
    }
    ++result.counters.regularized;
    outComponents.push_back(std::move(bImpl));
  }

  if (firstFatal) {
    // Fail-closed: a recorded reason, no partial output.
    result.fatal = firstFatal;
    result.detail = std::move(firstDetail);
    return result;
  }

  // 6. COMPOSE BACK by concatenation (no fusion).
  result.impl = ComposeComponents(outComponents);
  return result;
}

// Test hook: exercise B's ported mechanism directly (enumeration + coupled
// winding) so it can be graded against the fragment's recorded numbers.
CandidateBProbe RegularizeB_Probe(const Manifold::Impl& dirty,
                                  const std::vector<vec3>& probes,
                                  const vec3& seed) {
  CandidateBProbe out;
  const BEnumeration enu = EnumerateSelfCrossings(dirty);
  out.seamCount = enu.seamCount;
  out.boundaryTouchPairs = enu.boundaryTouchPairs;
  for (int c : DetectCoplanarClusters(dirty))
    if (c >= 0) ++out.coplanarClusterFaces;
  out.probeWinding.reserve(probes.size());
  for (const vec3& p : probes) {
    const std::optional<int> w = WindingAt(dirty, p, seed);
    out.probeWinding.push_back(w.has_value() ? *w : kWindingUncertain);
  }
  return out;
}

RegularizeResult RegularizeDirtyDirect(const Manifold::Impl& soup, double eps) {
  RegularizeResult result;
  if (eps <= 0.0) eps = EpsilonFromScale(soup.bBox_.Scale(), 1000);
  result.counters.components = 1;
  result.counters.dirty = 1;
  StageResult<Manifold::Impl> bRes = RunCandidateB(soup, eps);
  if (!bRes.ok()) {
    ++result.counters.failClosed;
    result.fatal = bRes.fatal;
    result.detail = std::move(bRes.detail);
    return result;
  }
  Manifold::Impl bImpl = std::move(*bRes.value);
  bImpl.epsilon_ = eps;
  if (GateComponent(bImpl) != GateVerdict::Clean) {
    ++result.counters.failClosed;
    result.fatal = FatalReason::NonManifoldEmission;
    result.detail = "candidate B output failed the re-gate";
    return result;
  }
  ++result.counters.regularized;
  result.impl = std::move(bImpl);
  return result;
}

}  // namespace manifold
