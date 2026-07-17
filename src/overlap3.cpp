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

// RemoveOverlaps3D: the regularization operator (valid oriented soup ->
// boundary of {w_S >= 1}).  Phases: decompose -> gate -> planarize ->
// arrange -> wind -> emit -> re-gate.  Design: docs/Regularize3D.md.

#include "overlap3.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
#include "polygon_internal.h"
#include "shared.h"

namespace manifold {

namespace {

struct OutTri3D {
  vec3 v[3];
};

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
// MEASUREMENT ONLY (f4-b1): the once-only triple positions, so the emission-
// wall census can classify each open-boundary edge by triple-point incidence
// (design-b's decisive census on a BUILT arrangement).  Populated by
// EnumerateTriplePoints under F4B_DUMP; empty in production.
static std::set<std::tuple<double, double, double>> gF4BTriplePts;

// `sheetN` (OPTIONAL, per-triangle, parallel to tv): pre-weld ORIENTED sheet
// normals - the plane-group provenance of each emitted triangle (the group's
// representative normal, signed by the emitted winding).  When present, the
// NEAR-TANGENT RADIAL BRANCH is enabled:
//   1. fan ring angles come from the provenance ray cross(n, heDir) instead
//      of the third-vertex chord (the eps-weld bends sub-eps-thin triangles,
//      so their chord reads up to radians wrong - measured 2.4 rad at a
//      micro corner - while sheets 0.004 rad apart must order correctly);
//   2. doubled directed edges surviving the vertex split (a true X-contact
//      whose two material wedges legitimately reconnect around BOTH
//      endpoints) are repaired by SUBDIVIDING one copy's triangle pair at a
//      point on the edge - a pointwise-identical surface in representable
//      form.
// Identity/provenance-over-distance: the pairing decisions use pre-weld
// information the weld destroys, never wider tolerances.  A same-group
// forward/backward pair does NOT imply a continuing sheet (oracle-refuted:
// a crease pairing each side with a transversal wall also presents as
// 1fwd+1bwd of one group), so no identity-based pre-pairing exists - the
// ring alternation on honest angles is the pairing rule.  When `sheetN` is
// absent (the per-face production path), behavior is BYTE-IDENTICAL to the
// chord-angle form.
bool SplitTouchingSheets(std::vector<vec3>& verts, Vec<ivec3>& tv,
                         const std::vector<vec3>* sheetN) {
  const int nTri = static_cast<int>(tv.size());
  auto heFrom = [&](int h) { return tv[h / 3][h % 3]; };
  auto heTo = [&](int h) { return tv[h / 3][(h % 3 + 1) % 3]; };

  std::map<std::pair<int, int>, std::vector<int>> edge2He;
  for (int h = 0; h < 3 * nTri; ++h) {
    const int a = heFrom(h), b = heTo(h);
    edge2He[{std::min(a, b), std::max(a, b)}].push_back(h);
  }

  // MEASUREMENT ONLY (f4-b1): F4B_DUMP counts the emission fan anomalies at
  // this wall (open boundary = fwd!=bwd, slivers, radial ties, material
  // overlaps) instead of failing closed on the first, prints one census line,
  // then still fails closed.  Env-gated: production behavior (the plain `return
  // false`) is byte-unchanged when it is unset.
  static const bool kF4BDump = std::getenv("F4B_DUMP") != nullptr;
  int cOpen = 0, cSliver = 0, cTie = 0, cOverlap = 0;
  int cOpenAtTriple = 0, cOpenOneTriple = 0, cOpenNoTriple = 0;
  int cOpenFan1 = 0, cOpenFan2 = 0, cOpenFan3 = 0, cOpenFan4 = 0,
      cOpenFanBig = 0;
  double tieMinGap = 1e300;
  auto atTriple = [&](int v) {
    const vec3& p = verts[v];
    return gF4BTriplePts.count({p.x, p.y, p.z}) > 0;
  };
  std::vector<int> pairedHe(3 * nTri, -1);
  for (const auto& [edge, hes] : edge2He) {
    bool badFan = false;
    std::vector<int> fwd, bwd;
    for (const int h : hes) (heFrom(h) == edge.first ? fwd : bwd).push_back(h);
    if (fwd.size() != bwd.size()) {
      if (!kF4BDump) return false;
      ++cOpen;
      const int nt =
          (atTriple(edge.first) ? 1 : 0) + (atTriple(edge.second) ? 1 : 0);
      if (nt == 2)
        ++cOpenAtTriple;
      else if (nt == 1)
        ++cOpenOneTriple;
      else
        ++cOpenNoTriple;
      const size_t fan = hes.size();
      if (fan == 1)
        ++cOpenFan1;
      else if (fan == 2)
        ++cOpenFan2;
      else if (fan == 3)
        ++cOpenFan3;
      else if (fan == 4)
        ++cOpenFan4;
      else
        ++cOpenFanBig;
      if (std::getenv("F4B_OPENDUMP") != nullptr) {
        std::fprintf(stderr, "F4B_OPENEDGE fan=%zu fwd=%zu bwd=%zu faces=[",
                     hes.size(), fwd.size(), bwd.size());
        for (const int h : hes) std::fprintf(stderr, "%d ", h / 3);
        std::fprintf(stderr, "] p=(%.6g,%.6g,%.6g)\n", verts[edge.first].x,
                     verts[edge.first].y, verts[edge.first].z);
      }
      continue;
    }
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
      // Ring direction: the third-vertex CHORD, in BOTH modes.  A
      // provenance-ray variant (cross(sheetN, heDir)) was built and
      // MEASURED-REFUTED: at the micro X-contact the fan structure is
      // RADIUS-DEPENDENT (the near-tangent wedge sheets pass within
      // ~1e-11 of the edge, not through it, so their crossing order at
      // r=1e-8 is the REVERSE of their r->0 ray order) - no single per-
      // sheet angle exists and the ray ordering broke alternation on a fan
      // the chord ordering pairs oracle-TRUE.  The chord pairing was
      // verified against exact sector windings on every multi-sheet fan of
      // the carrier corpus.
      const int c = tv[h / 3][(h % 3 + 2) % 3];
      vec3 d = verts[c] - pa;
      d -= la::dot(d, ax) * ax;
      const double len = la::length(d);
      if (len == 0.0) {  // sliver: third vert on the edge line
        if (!kF4BDump) return false;
        ++cSliver;
        badFan = true;
        break;
      }
      d /= len;
      if (ring.empty()) {
        u = d;
        v = la::cross(ax, u);
      }
      double angle = std::atan2(la::dot(d, v), la::dot(d, u));
      if (angle < 0.0) angle += kTwoPi;
      ring.push_back({angle, h, heFrom(h) == edge.first});
    }
    if (badFan) continue;  // dump mode: sliver counted, skip this fan
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
      if (gap < kAngleTie) {
        if (!kF4BDump) return false;
        ++cTie;
        tieMinGap = std::min(tieMinGap, gap);
        badFan = true;
        break;
      }
      if (cur.fwd == nxt.fwd) {  // material overlap
        if (!kF4BDump) return false;
        ++cOverlap;
        std::fprintf(
            stderr, "F4B_OVERLAP at (%.9g,%.9g,%.9g)->(%.9g,%.9g,%.9g) ring:",
            verts[edge.first].x, verts[edge.first].y, verts[edge.first].z,
            verts[edge.second].x, verts[edge.second].y, verts[edge.second].z);
        for (const RingEntry& r : ring) {
          std::fprintf(stderr, " %.6f%s", r.angle, r.fwd ? "f" : "b");
          {
            const vec3& n = (*sheetN)[r.he / 3];
            std::fprintf(stderr, "[n=%.3g,%.3g,%.3g]", n.x, n.y, n.z);
          }
        }
        std::fprintf(stderr, "\n");
        badFan = true;
        break;
      }
    }
    if (badFan) continue;  // dump mode: tie/overlap counted, skip this fan
    for (int i = 0; i < k; ++i) {
      if (ring[i].fwd) continue;
      const RingEntry& partner = ring[(i + 1) % k];  // next CCW is forward
      pairedHe[ring[i].he] = partner.he;
      pairedHe[partner.he] = ring[i].he;
    }
  }
  if (kF4BDump) {
    std::fprintf(stderr,
                 "F4B_CENSUS openEdges=%d slivers=%d ties=%d overlaps=%d\n",
                 cOpen, cSliver, cTie, cOverlap);
    std::fprintf(
        stderr,
        "F4B_OPEN atTriple=%d oneTriple=%d noTriple=%d | fan1=%d fan2=%d "
        "fan3=%d fan4=%d fanBig=%d | tieMinGap=%g\n",
        cOpenAtTriple, cOpenOneTriple, cOpenNoTriple, cOpenFan1, cOpenFan2,
        cOpenFan3, cOpenFan4, cOpenFanBig, cTie ? tieMinGap : 0.0);
    if (cOpen + cSliver + cTie + cOverlap > 0) return false;
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
  // DOUBLED-EDGE SUBDIVISION (radial branch only): a true X-contact edge
  // whose two material wedges legitimately RECONNECT around both endpoints
  // survives the vertex split as two undirected copies of one vertex pair -
  // a surface the halfedge representation cannot carry (Is2Manifold forbids
  // duplicate directed edges).  Subdivide every copy beyond the first at a
  // distinct interior point of the (shared) segment: the surface is
  // pointwise unchanged, the topology becomes representable, and the
  // pairing already proved each copy a coherent two-triangle sheet pair.
  // (Post-flip: the engine's provenance is the ONLY mode - the identity-
  // free chord-angle-only caller died with the per-face world.)
  {
    std::map<std::pair<int, int>, std::vector<int>> und;  // undirected -> hes
    for (int h = 0; h < 3 * nTri; ++h) {
      const int a = tv[h / 3][h % 3], b = tv[h / 3][(h % 3 + 1) % 3];
      und[{std::min(a, b), std::max(a, b)}].push_back(h);
    }
    for (const auto& kv : und) {
      if (kv.second.size() <= 2) continue;
      // group the halfedges into their paired two-triangle sheets; keep the
      // first pair on the original edge, subdivide each further pair at a
      // distinct parameter (1/2, 1/3, ...).
      std::set<int> seen;
      int extra = 0;
      for (const int h : kv.second) {
        if (seen.count(h)) continue;
        const int p = pairedHe[h];
        seen.insert(h);
        seen.insert(p);
        // a triangle already subdivided under another doubled edge no longer
        // carries this halfedge's original corners - leave it (the manifold
        // gate stays the fail-closed backstop for the unhandled residue)
        auto still = [&](int hh) {
          const int a = tv[hh / 3][hh % 3], b = tv[hh / 3][(hh % 3 + 1) % 3];
          return std::minmax(a, b) ==
                 std::minmax(kv.first.first, kv.first.second);
        };
        if (!still(h) || !still(p)) continue;
        if (extra++ == 0) continue;  // first copy keeps the edge
        const double tSplit = 1.0 / static_cast<double>(extra);
        const int u = kv.first.first, w = kv.first.second;
        const int m = static_cast<int>(verts.size());
        verts.push_back(verts[u] + tSplit * (verts[w] - verts[u]));
        for (const int hh : {h, p}) {
          const int t = hh / 3, k = hh % 3;
          const ivec3 tri = tv[t];
          // replace edge (tri[k], tri[k+1]) by (tri[k], m) in place and
          // append the (m, tri[k+1]) half
          ivec3 t2 = tri;
          t2[k] = m;
          tv[t][(k + 1) % 3] = m;
          tv.push_back(t2);
        }
      }
    }
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

// `sheetN` (OPTIONAL, parallel to `tris`): pre-weld ORIENTED sheet normals,
// carried through the weld/filter into SplitTouchingSheets' fan pairing (the
// near-tangent radial branch).  Absent = byte-identical to the identity-free
// pipeline.
StageResult<Manifold::Impl> BuildImpl(const std::vector<OutTri3D>& tris,
                                      double eps,
                                      const std::vector<vec3>* sheetN) {
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
  std::vector<vec3> tvN;  // sheet normal per KEPT triangle (when supplied)
  std::set<std::tuple<int, int, int>> seenTris;
  for (size_t i = 0; i < tris.size(); ++i) {
    const OutTri3D& tri = tris[i];
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
    tvN.push_back((*sheetN)[i]);
  }

  // Touching sheets separate BEFORE the topology is built (spec COPLANAR
  // implementation close: touching contacts).
  if (!SplitTouchingSheets(verts, tv, &tvN)) {
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
    if (std::getenv("E1_DUMP") != nullptr) {
      std::fprintf(stderr, "E1 GATE IsManifold=%d Is2Manifold=%d\n",
                   impl.IsManifold() ? 1 : 0, impl.Is2Manifold() ? 1 : 0);
      std::map<std::pair<int, int>, int> dcnt;
      for (const ivec3& t : tv)
        for (int k = 0; k < 3; ++k) ++dcnt[{t[k], t[(k + 1) % 3]}];
      for (const auto& kv : dcnt)
        if (kv.second > 1) {
          const vec3& p = verts[kv.first.first];
          const vec3& q = verts[kv.first.second];
          std::fprintf(stderr,
                       "E1 DUPEDGE x%d v%d(%.9g,%.9g,%.9g) -> "
                       "v%d(%.9g,%.9g,%.9g)\n",
                       kv.second, kv.first.first, p.x, p.y, p.z,
                       kv.first.second, q.x, q.y, q.z);
        }
    }
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

// Forward decl (defined below): per-face exactly-coplanar overlap cluster id
// (or -1).  Used by GateComponent to route a WITHIN-component coplanar-overlap
// component to the resolver (the R2(i) blind spot the self-intersection test
// misses).
std::vector<int> DetectCoplanarClusters(const Manifold::Impl& in);

// Forward decl (defined below): does `in` carry ANY within-component coplanar
// overlap?  Exactly DetectCoplanarClusters' pass-1 seed existence: pass 2 only
// EXTENDS pass-1 clusters (it unites only when a side is already clustered), so
// a face gets a cluster id iff pass 1 seeded one.  The gate reads this bool
// with an early-out instead of computing the full labeling (which the resolver
// still does on the dirty path).
bool HasCoplanarOverlap(const Manifold::Impl& in);

// Split `in` into connected components by halfedge connectivity - the Decompose
// primitive (constructors.cpp:455) mirrored at the Impl level so the operator
// never round-trips through the CSG layer.  Each returned component is a
// finished Impl (bbox/normals/collider) whose epsilon_ is pinned to the
// resolved global machine scale, so the per-component gate (IsSelfIntersecting
// reads collider_, faceNormal_, epsilon_) runs directly.  A single connected
// input returns exactly one component that IS a copy of `in`, so a clean
// input of one component passes through unchanged.
//
// NON-FUSION POSTURE (docs/Regularize3D.md contract): the connectivity split is
// the unit of scope, PERIOD.  Cross-component interaction - overlapping,
// touching, coplanar OR transversal - is NEVER regularized here; two distinct
// components that happen to coplanar-overlap (a buried plug, coplanar stacked
// boxes) are each regularized on their own and composed back by concatenation,
// never united.  A global {w_S>=1} read WOULD fuse them, but that is the
// Boolean's job (already done upstream in any operation chain); the per-
// component read is the deliberate choice, matching the touching-contact
// posture.  The coplanar fold reaches a defect only when it is INTERNAL to one
// connected component (a doubled wall / folded flap stitched into the shell).
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
  const int numComponents = uf.connectedComponents(vertLabel);

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
  // WITHIN-COMPONENT coplanar overlap.  A single connected component that
  // passes the self-intersection test can still carry a pure COPLANAR overlap
  // between its OWN faces (a doubled internal wall / folded-flat flap whose two
  // sheets coincide): IsSelfIntersecting does NOT flag coplanar coincidence
  // (docs/Regularize3D.md R2(i)), so such a component would early-exit "clean"
  // yet is NOT the {w_S>=1} boundary.  DetectCoplanarClusters (bbox-overlap
  // prefilter, then exact orient3d coplanarity + a 2D-area overlap witness),
  // run on THIS component only, routes it to the resolver where the
  // exact-coplanar fold consumes the overlap.  This is WITHIN-component by
  // construction: `comp` is one connectivity component, so a CROSS-component
  // coplanar overlap (two distinct components that happen to coincide) is
  // invisible here BY DESIGN - the non-fusion posture (fusion is the Boolean's
  // job).  A clean solid with no internal coplanar overlap detects nothing and
  // stays Clean (bitwise pass-through); the cost is the prefiltered scan,
  // proportional to the input.
  if (HasCoplanarOverlap(comp)) return GateVerdict::Dirty;
  return GateVerdict::Clean;
}

// ---------------------------------------------------------------------------
// The resolver mechanism (docs/Regularize3D.md "the resolver's mechanism"):
// operand-agnostic ENUMERATION (level-0 pierce predicates through a static
// Shewchuk filter) + coupled integer-delta WINDING.  Every crossing DECISION is
// a level-0 orient3d on INPUT coordinates.  The cell-complex + halfedge
// {w_S>=1} boundary EMISSION (THE BUILD) sits on top.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// SINGLE GLOBAL TIE-BREAK CONVENTION (SoS), docs/Regularize3D.md stage 6.
// ONE exact integer implementation, threaded through every enumeration
// predicate so a ray or edge grazing a shared boundary resolves identically for
// every probe.  Each coord is decomposed to (mantissa * 2^exp) via frexp; the
// perturbed 4x4 orient3d cascade (the e^0 exact orient3d plus the symbolic
// higher-e terms) evaluates on one integer path - products of <= 3 mantissas
// accumulated in an adaptive-width two's-complement bigint.  TOTAL for every
// finite-double input: no window-fail, no refusal (0 means an exact geometric
// tie, never "uncertain").  The derivation - the e^K expansion, the
// accumulator-width bound, and the local-rank reduction that lets a per-
// predicate LOCAL rank stand in for the global one - is in docs/Regularize3D.md
// stage 6.  Only invoked on the RARE filter-uncertain (0) fallback, so the
// certified fast path is untouched.
//
// RELUCTANT ACCEPTANCE (owner contract): this integer exact kernel is net-new
// surface the design had declined ("no exact kernel in the tree").  Accepted
// NARROWLY as the stage-6 filter-0 fallback - one implementation, single call
// site discipline, off every certified path.  QUEUED FOR REVISIT
// (docs/Regularize3D.md open list): whether the arrangement can be structured
// to avoid an exact orient3d kernel at all remains open.
namespace sos {
constexpr int kPerm[24][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
    {0, 3, 2, 1}, {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0},
    {1, 3, 0, 2}, {1, 3, 2, 0}, {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3},
    {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0}, {3, 0, 1, 2}, {3, 0, 2, 1},
    {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}};
inline int Parity(const int s[4]) {
  int p = 1;
  for (int i = 0; i < 4; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (s[i] > s[j]) p = -p;
  return p;
}
// A determinant term: sign * |product of <= NMag*(64/53) mantissas| * 2^e.  The
// magnitude limb width NMag is a compile-time parameter of the ONE homogeneous
// form (homog-design): NMag==3 holds a <=3-factor product (2^159, the
// historical orient3d) and NMag==8 holds a <=9-factor product (2^477, the
// degree-9 constructed-point orient2d).  The degree-3 instantiation is
// bit-for-bit the historical width-3 term.
template <int NMag>
struct TermN {
  uint64_t mag[NMag];
  long e;
  int sign;
};
using Term = TermN<3>;  // the historical width-3 term (orient3d / SoS)
// Adaptive two's-complement accumulator width.  Each coord is m*2^E with |m| <
// 2^53 and E = ex - 53, ex the frexp exponent in [-1073, 1024], so E in
// [-1126, 971].  For the degree-3 form a three-factor esum in [-3378, 2913]
// gives a <= 6291-bit spread and 112 limbs was a provable bound.  The degree-9
// constructed-point orient2d (homog-design) needs more: 9*53 = 477 magnitude
// bits and a worst-case exponent spread 9*(971-(-1126)) = 18873 bits, so
// (477 + 18873)/64 + 3 ~= 305 limbs -> 320 keeps the "total by construction,
// never window-fail" guarantee for BOTH forms.  Measured active width on
// mesh-scale data (incl. the near-parallel wedge family): 12 limbs.  The active
// limb count is computed per call from the actual spread; the fixed storage
// just guarantees totality (a ~2.5KB stack array; zero runtime change on mesh
// data).
constexpr int kAccumLimbs = 320;
// |product of the cnt (nonzero) mantissas| -> mag[NMag]; returns the product
// sign.  NMag==3 is the historical <=3-factor orient3d specialization
// (unchanged, bit-for-bit); NMag>3 is the general schoolbook multiply used by
// the degree-9 constructed-point orient2d (<= 9 factors -> 8 limbs).
template <int NMag>
inline int MulMagN(const int64_t* pm, int cnt, uint64_t mag[NMag]) {
  if constexpr (NMag == 3) {
    mag[0] = 1;
    mag[1] = 0;
    mag[2] = 0;
    int sg = 1;
    uint64_t a[3];
    for (int i = 0; i < cnt; ++i) {
      a[i] = (uint64_t)(pm[i] < 0 ? -pm[i] : pm[i]);
      if (pm[i] < 0) sg = -sg;
    }
    if (cnt == 1) {
      mag[0] = a[0];
    } else if (cnt == 2) {
      const unsigned __int128 p = (unsigned __int128)a[0] * a[1];
      mag[0] = (uint64_t)p;
      mag[1] = (uint64_t)(p >> 64);
    } else if (cnt == 3) {
      const unsigned __int128 p12 = (unsigned __int128)a[0] * a[1];
      const uint64_t lo = (uint64_t)p12, hi = (uint64_t)(p12 >> 64);
      const unsigned __int128 pLo = (unsigned __int128)lo * a[2];
      const unsigned __int128 pHi = (unsigned __int128)hi * a[2] + (pLo >> 64);
      mag[0] = (uint64_t)pLo;
      mag[1] = (uint64_t)pHi;
      mag[2] = (uint64_t)(pHi >> 64);
    }
    return sg;
  } else {
    for (int L = 0; L < NMag; ++L) mag[L] = 0;
    mag[0] = 1;
    int sg = 1;
    for (int i = 0; i < cnt; ++i) {
      const uint64_t f = (uint64_t)(pm[i] < 0 ? -pm[i] : pm[i]);
      if (pm[i] < 0) sg = -sg;
      unsigned __int128 carry = 0;
      for (int L = 0; L < NMag; ++L) {
        const unsigned __int128 p = (unsigned __int128)mag[L] * f + carry;
        mag[L] = (uint64_t)p;
        carry = p >> 64;
      }
    }
    return sg;
  }
}
// Historical name kept for the width-3 callers (orient3d / SoS): a forward to
// the specialization above (instruction-identical after inlining).
inline int MulMag(const int64_t* pm, int cnt, uint64_t mag[3]) {
  return MulMagN<3>(pm, cnt, mag);
}
// acc[0..nLimbs) += sign * (mag << shift), two's complement.  NMag = term
// magnitude width (3 = orient3d, 8 = degree-9 orient2d).
template <int NMag>
inline void AddShiftedMagN(uint64_t* acc, int nLimbs, const uint64_t mag[NMag],
                           int shift, int sign) {
  const int limbShift = shift / 64, bitShift = shift % 64;
  uint64_t tmp[kAccumLimbs] = {0};
  for (int i = 0; i < NMag; ++i) {
    const int dst = i + limbShift;
    if (dst >= 0 && dst < nLimbs) {
      tmp[dst] |= mag[i] << bitShift;
      if (bitShift && dst + 1 < nLimbs)
        tmp[dst + 1] |= mag[i] >> (64 - bitShift);
    }
  }
  if (sign < 0) {
    unsigned __int128 carry = 1;
    for (int i = 0; i < nLimbs; ++i) {
      const unsigned __int128 s = (unsigned __int128)(~tmp[i]) + carry;
      tmp[i] = (uint64_t)s;
      carry = s >> 64;
    }
  }
  unsigned __int128 carry = 0;
  for (int i = 0; i < nLimbs; ++i) {
    const unsigned __int128 s = (unsigned __int128)acc[i] + tmp[i] + carry;
    acc[i] = (uint64_t)s;
    carry = s >> 64;
  }
}
// Sign of the exact integer sum of the terms.  Adaptive width, total (no
// refusal): the true sum fits the two's-complement window by construction. NMag
// = term magnitude width; the magnitude-bit budget is NMag*64 (192 for the
// historical width-3 form, 512 for the degree-9 form).
template <int NMag>
inline int SumSignN(const TermN<NMag>* t, int nT) {
  if (nT == 0) return 0;
  long emin = t[0].e, emax = t[0].e;
  for (int i = 1; i < nT; ++i) {
    emin = std::min(emin, t[i].e);
    emax = std::max(emax, t[i].e);
  }
  int nLimbs = (int)((NMag * 64 + (emax - emin)) / 64) + 3;
  if (nLimbs > kAccumLimbs) nLimbs = kAccumLimbs;  // never hit (proven bound)
  uint64_t acc[kAccumLimbs] = {0};
  for (int i = 0; i < nT; ++i)
    AddShiftedMagN<NMag>(acc, nLimbs, t[i].mag, (int)(t[i].e - emin),
                         t[i].sign);
  if (acc[nLimbs - 1] >> 63) return -1;
  for (int i = 0; i < nLimbs; ++i)
    if (acc[i]) return 1;
  return 0;
}
// Historical name kept for the width-3 callers (orient3d / SoS).
inline int SumSign(const Term* t, int nT) { return SumSignN<3>(t, nT); }
inline void Decompose(const double pts[4][3], int64_t M[4][3], int E[4][3]) {
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 3; ++c) {
      const double d = pts[r][c];
      if (d == 0.0) {
        M[r][c] = 0;
        E[r][c] = 0;
        continue;
      }
      int ex;
      const double f = std::frexp(d, &ex);
      M[r][c] = (int64_t)std::ldexp(f, 53);
      E[r][c] = ex - 53;
    }
}
// Homogeneous companion of Decompose: split each of the 4x4 homogeneous
// coordinates d into (mantissa M, exponent E) with d == M*2^E, |M| < 2^53.
inline void DecomposeH(const double pts[4][4], int64_t M[4][4], int E[4][4]) {
  for (int r = 0; r < 4; ++r)
    for (int c = 0; c < 4; ++c) {
      const double d = pts[r][c];
      if (d == 0.0) {
        M[r][c] = 0;
        E[r][c] = 0;
        continue;
      }
      int ex;
      const double f = std::frexp(d, &ex);
      M[r][c] = (int64_t)std::ldexp(f, 53);
      E[r][c] = ex - 53;
    }
}
// The ONE blessed exact predicate FORM (homog-design): the sign of the
// homogeneous orientation determinant of four points (X,Y,Z,W), corrected by
// sign(prod W_i) - see the TRIPWIRE at Orient3DExactSign.  TrivialW=true is the
// INPUT-POINT orient3d instantiation: each point has W==1 (an input double
// point is the intersection of its three trivial axis planes, denominator 1),
// so the W column is the literal ones column (skipped exactly as the historical
// ExactOrient3D) and sign(prod W)=+1 - the loop is then BIT-FOR-BIT the
// historical orient3d (24 permutations, <=3-factor width-3 accumulation, the 0
// iff the four points are exactly coplanar).  TrivialW=false keeps the W column
// (4-factor terms) and multiplies in the weight-product sign; it is not
// instantiated in production (the two production instantiations are this
// TrivialW=true orient3d and the degree-9 HomogOrient2DExact).  TOTAL: no
// window-fail; an exact 0 is a genuine geometric tie.
template <bool TrivialW>
inline int HomogOrient3DSign(const double pts[4][4]) {
  int64_t M[4][4];
  int E[4][4];
  DecomposeH(pts, M, E);
  constexpr int NMag = TrivialW ? 3 : 4;
  TermN<NMag> t[24];
  int nT = 0;
  for (const auto& s : kPerm) {
    int64_t pm[4];
    long pe = 0;
    int cnt = 0;
    bool zero = false;
    for (int r = 0; r < 4; ++r) {
      const int c = s[r];
      if constexpr (TrivialW) {
        if (c == 3) continue;  // the W column is the literal ones column
      }
      if (M[r][c] == 0) {
        zero = true;
        break;
      }
      pm[cnt++] = M[r][c];
      pe += E[r][c];
    }
    if (zero) continue;
    t[nT].sign = Parity(s) * MulMagN<NMag>(pm, cnt, t[nT].mag);
    t[nT].e = pe;
    ++nT;
  }
  const int detSign = SumSignN<NMag>(t, nT);
  if constexpr (TrivialW) {
    return detSign;  // prod W_i == 1, sign +1
  } else {
    int w = 1;
    for (int r = 0; r < 4; ++r)
      if (pts[r][3] < 0.0) w = -w;
    return detSign * w;
  }
}
// Symbolically-perturbed orient3d sign (never 0 for distinct idx).  Enumerates
// all 24*8 monomials, groups by the e-exponent K, and returns the sign of the
// lowest-K nonzero coefficient.  When the e^0 (K==0) coefficient - the exact
// orient3d - is nonzero this returns that exact sign; otherwise the
// perturbation decides.  idx = the four points' global vertex indices
// (distinct).
inline int SoSOrient3D(const double pts[4][3], const int idx[4]) {
  int64_t M[4][3];
  int E[4][3];
  Decompose(pts, M, E);
  int order[4] = {0, 1, 2, 3};
  std::sort(order, order + 4, [&](int a, int b) { return idx[a] < idx[b]; });
  int rankOf[4];
  for (int r = 0; r < 4; ++r) rankOf[order[r]] = r;
  auto keyOf = [&](int r, int c) { return 1 << (rankOf[r] * 3 + (2 - c)); };
  struct Mono {
    int K;
    Term t;
  };
  Mono mono[192];
  int nM = 0;
  for (const auto& s : kPerm) {
    const int par = Parity(s);
    int rr[3], cc[3], nf = 0;
    for (int r = 0; r < 4; ++r)
      if (s[r] != 3) {  // the three non-ones factors
        rr[nf] = r;
        cc[nf] = s[r];
        ++nf;
      }
    for (int mask = 0; mask < 8; ++mask) {  // subset perturbed vs real
      int K = 0;
      int64_t pm[3];
      long pe = 0;
      int cnt = 0;
      bool zero = false;
      for (int f = 0; f < nf; ++f) {
        if (mask & (1 << f)) {
          K += keyOf(rr[f], cc[f]);  // perturbed factor
        } else {
          const int64_t m = M[rr[f]][cc[f]];
          if (m == 0) {
            zero = true;
            break;
          }
          pm[cnt++] = m;
          pe += E[rr[f]][cc[f]];
        }
      }
      if (zero) continue;
      mono[nM].K = K;
      mono[nM].t.sign = par * MulMag(pm, cnt, mono[nM].t.mag);
      mono[nM].t.e = pe;
      ++nM;
    }
  }
  std::sort(mono, mono + nM,
            [](const Mono& a, const Mono& b) { return a.K < b.K; });
  int i = 0;
  while (i < nM) {
    Term grp[192];
    int ng = 0;
    const int k0 = mono[i].K;
    while (i < nM && mono[i].K == k0) grp[ng++] = mono[i++].t;
    const int sg = SumSign(grp, ng);
    if (sg != 0) return sg;
  }
  return 0;  // unreachable for distinct points; caller treats 0 as fail-closed
}

// ===========================================================================
// THE GENERAL (degree-9) INSTANTIATION of the ONE homogeneous form (homog-
// design): exact orient2d of three IN-FACE crossing points, each the Cramer
// intersection of a plane triple {F,g,h}.  A crossing point's homogeneous
// coords (X,Y,Z,W) are 3x3 determinants of the plane coefficients (degree 3);
// the 3x3 orient2d determinant of three of them is degree 9, summed on the SAME
// adaptive accumulator (SumSignN<8>) and corrected by sign(W0 W1 W2). Validated
// vs exact rationals over 1.3M random + near-parallel-wedge plane triples: ZERO
// disagreements, ZERO wrong filter certs (harness h2).  An INPUT vertex is the
// trivial triple (W==1), so mixed orient2d (a constructed crossing vs input
// triangle vertices) is the same form at lower degree.
//
// A degree-<=9 monomial: |product of nf mantissas| * sign * 2^e.
struct Mono2D {
  uint64_t mant[9];  // abs mantissas, each < 2^53
  int nf;            // number of factors (<= 9)
  long e;            // sum of exponents
  int sign;          // +-1
};
using Poly = std::vector<Mono2D>;  // a sum of monomials
inline Poly PConst(double d) {
  if (d == 0.0) return {};
  int ex;
  const double f = std::frexp(d, &ex);
  const int64_t M = (int64_t)std::ldexp(f, 53);
  Mono2D m;
  m.nf = 1;
  m.mant[0] = (uint64_t)(M < 0 ? -M : M);
  m.e = ex - 53;
  m.sign = M < 0 ? -1 : 1;
  return {m};
}
inline Poly PMul(const Poly& a, const Poly& b) {
  Poly r;
  r.reserve(a.size() * b.size());
  for (const auto& x : a)
    for (const auto& y : b) {
      Mono2D m;
      m.nf = x.nf + y.nf;  // <= 9 by construction
      for (int i = 0; i < x.nf; ++i) m.mant[i] = x.mant[i];
      for (int i = 0; i < y.nf; ++i) m.mant[x.nf + i] = y.mant[i];
      m.e = x.e + y.e;
      m.sign = x.sign * y.sign;
      r.push_back(m);
    }
  return r;
}
inline Poly PAdd(Poly a, const Poly& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}
inline Poly PNeg(Poly a) {
  for (auto& m : a) m.sign = -m.sign;
  return a;
}
inline Poly PSub(const Poly& a, const Poly& b) { return PAdd(a, PNeg(b)); }
// Sign of the exact integer sum of a Poly's monomials, on the ONE accumulator.
inline int PolySumSign(const Poly& p) {
  if (p.empty()) return 0;
  std::vector<TermN<8>> t;
  t.reserve(p.size());
  for (const auto& m : p) {
    int64_t pm[9];
    for (int i = 0; i < m.nf; ++i) pm[i] = (int64_t)m.mant[i];  // < 2^53
    TermN<8> term;
    term.sign = m.sign * MulMagN<8>(pm, m.nf, term.mag);
    term.e = m.e;
    t.push_back(term);
  }
  return SumSignN<8>(t.data(), (int)t.size());
}
// A homogeneous point as Polys (X,Y,Z,W).  A plane triple point (cramer) is
// degree 3; an input vertex (Trivial) is degree 1 with W==1.
struct HPoint {
  Poly X, Y, Z, W;
};
inline HPoint TrivialHPoint(const vec3& v) {
  return {PConst(v.x), PConst(v.y), PConst(v.z), PConst(1.0)};
}
// The crossing point F /\ g /\ h (matches Intersect3Planes: c12=cross(g,h),
// c20=cross(h,F), c01=cross(F,g); W=F.c12; X=dF*c12+dg*c20+dh*c01).
inline HPoint CramerHPoint(const vec3& nF, double dF, const vec3& ng, double dg,
                           const vec3& nh, double dh) {
  const Poly Fx = PConst(nF.x), Fy = PConst(nF.y), Fz = PConst(nF.z);
  const Poly Gx = PConst(ng.x), Gy = PConst(ng.y), Gz = PConst(ng.z);
  const Poly Hx = PConst(nh.x), Hy = PConst(nh.y), Hz = PConst(nh.z);
  const Poly Fd = PConst(dF), Gd = PConst(dg), Hd = PConst(dh);
  auto cross = [](const Poly& ax, const Poly& ay, const Poly& az,
                  const Poly& bx, const Poly& by, const Poly& bz, Poly& cx,
                  Poly& cy, Poly& cz) {
    cx = PSub(PMul(ay, bz), PMul(az, by));
    cy = PSub(PMul(az, bx), PMul(ax, bz));
    cz = PSub(PMul(ax, by), PMul(ay, bx));
  };
  Poly c12x, c12y, c12z, c20x, c20y, c20z, c01x, c01y, c01z;
  cross(Gx, Gy, Gz, Hx, Hy, Hz, c12x, c12y, c12z);
  cross(Hx, Hy, Hz, Fx, Fy, Fz, c20x, c20y, c20z);
  cross(Fx, Fy, Fz, Gx, Gy, Gz, c01x, c01y, c01z);
  HPoint p;
  p.W = PAdd(PAdd(PMul(Fx, c12x), PMul(Fy, c12y)), PMul(Fz, c12z));
  p.X = PAdd(PAdd(PMul(Fd, c12x), PMul(Gd, c20x)), PMul(Hd, c01x));
  p.Y = PAdd(PAdd(PMul(Fd, c12y), PMul(Gd, c20y)), PMul(Hd, c01y));
  p.Z = PAdd(PAdd(PMul(Fd, c12z), PMul(Gd, c20z)), PMul(Hd, c01z));
  return p;
}
inline const Poly& KeepPoly(const HPoint& p, int axis, int which) {
  const int a0 = (axis == 0) ? 1 : 0;
  const int a1 = (axis == 2) ? 1 : 2;
  const int a = which == 0 ? a0 : a1;
  return a == 0 ? p.X : (a == 1 ? p.Y : p.Z);
}
// EXACT orient2d(P0,P1,P2) in face F's plane, dropping `axis` (the face
// normal's dominant axis).  = sign(det[[A0,B0,W0],...]) * sign(W0 W1 W2).
// Returns 0 iff the three points are exactly collinear in the face (a genuine
// coincidence - concurrency / aliased crossing - routed by the caller, never
// perturbed).
inline int HomogOrient2DExact(const HPoint& P0, const HPoint& P1,
                              const HPoint& P2, int axis) {
  const Poly& A0 = KeepPoly(P0, axis, 0);
  const Poly& B0 = KeepPoly(P0, axis, 1);
  const Poly& A1 = KeepPoly(P1, axis, 0);
  const Poly& B1 = KeepPoly(P1, axis, 1);
  const Poly& A2 = KeepPoly(P2, axis, 0);
  const Poly& B2 = KeepPoly(P2, axis, 1);
  const Poly det = PAdd(PSub(PMul(A0, PSub(PMul(B1, P2.W), PMul(B2, P1.W))),
                             PMul(A1, PSub(PMul(B0, P2.W), PMul(B2, P0.W)))),
                        PMul(A2, PSub(PMul(B0, P1.W), PMul(B1, P0.W))));
  const int sDet = PolySumSign(det);
  const int sW = PolySumSign(P0.W) * PolySumSign(P1.W) * PolySumSign(P2.W);
  return sDet * sW;
}

// --- construction-aware filter (EBD running forward-error bound) ------------
// The naive final-determinant permanent is UNSOUND: in the near-parallel wedge
// regime the constructed coords are ~eps by O(1) cancellation, so the permanent
// collapses to ~eps^3 while the error floor stays ~u (ratio blows up ~1e15).
// The certified filter propagates the running forward error THROUGH the
// construction; certified iff the det AND all three W's are sign-certified.
struct EBD {
  double v, err;  // |true - v| <= err
};
constexpr double kEbdU = 0x1p-53;
inline EBD EIn(double d) { return {d, 0.0}; }
inline EBD EAdd(EBD a, EBD b) {
  const double v = a.v + b.v;
  return {v, a.err + b.err + kEbdU * std::abs(v)};
}
inline EBD ESub(EBD a, EBD b) {
  const double v = a.v - b.v;
  return {v, a.err + b.err + kEbdU * std::abs(v)};
}
inline EBD EMul(EBD a, EBD b) {
  const double v = a.v * b.v;
  return {v, std::abs(a.v) * b.err + std::abs(b.v) * a.err + a.err * b.err +
                 kEbdU * std::abs(v)};
}
struct EHPoint {
  EBD X, Y, Z, W;
};
inline EHPoint ETrivialHPoint(const vec3& v) {
  return {EIn(v.x), EIn(v.y), EIn(v.z), EIn(1.0)};
}
inline EHPoint ECramerHPoint(const vec3& nF, double dF, const vec3& ng,
                             double dg, const vec3& nh, double dh) {
  auto cross = [](EBD ax, EBD ay, EBD az, EBD bx, EBD by, EBD bz, EBD& cx,
                  EBD& cy, EBD& cz) {
    cx = ESub(EMul(ay, bz), EMul(az, by));
    cy = ESub(EMul(az, bx), EMul(ax, bz));
    cz = ESub(EMul(ax, by), EMul(ay, bx));
  };
  const EBD Fx = EIn(nF.x), Fy = EIn(nF.y), Fz = EIn(nF.z);
  const EBD Gx = EIn(ng.x), Gy = EIn(ng.y), Gz = EIn(ng.z);
  const EBD Hx = EIn(nh.x), Hy = EIn(nh.y), Hz = EIn(nh.z);
  EBD c12x, c12y, c12z, c20x, c20y, c20z, c01x, c01y, c01z;
  cross(Gx, Gy, Gz, Hx, Hy, Hz, c12x, c12y, c12z);
  cross(Hx, Hy, Hz, Fx, Fy, Fz, c20x, c20y, c20z);
  cross(Fx, Fy, Fz, Gx, Gy, Gz, c01x, c01y, c01z);
  EHPoint p;
  p.W = EAdd(EAdd(EMul(Fx, c12x), EMul(Fy, c12y)), EMul(Fz, c12z));
  p.X =
      EAdd(EAdd(EMul(EIn(dF), c12x), EMul(EIn(dg), c20x)), EMul(EIn(dh), c01x));
  p.Y =
      EAdd(EAdd(EMul(EIn(dF), c12y), EMul(EIn(dg), c20y)), EMul(EIn(dh), c01y));
  p.Z =
      EAdd(EAdd(EMul(EIn(dF), c12z), EMul(EIn(dg), c20z)), EMul(EIn(dh), c01z));
  return p;
}
inline const EBD& EKeep(const EHPoint& p, int axis, int which) {
  const int a0 = (axis == 0) ? 1 : 0;
  const int a1 = (axis == 2) ? 1 : 2;
  const int a = which == 0 ? a0 : a1;
  return a == 0 ? p.X : (a == 1 ? p.Y : p.Z);
}
// Filter verdict: +-1 if certified, 0 if uncertain (escalate to exact).
inline int HomogOrient2DFilter(const EHPoint& P0, const EHPoint& P1,
                               const EHPoint& P2, int axis) {
  auto certSign = [](EBD e) -> int {
    if (std::abs(e.v) > e.err) return e.v > 0 ? 1 : -1;
    return 0;
  };
  const int sW0 = certSign(P0.W), sW1 = certSign(P1.W), sW2 = certSign(P2.W);
  const EBD& A0 = EKeep(P0, axis, 0);
  const EBD& B0 = EKeep(P0, axis, 1);
  const EBD& A1 = EKeep(P1, axis, 0);
  const EBD& B1 = EKeep(P1, axis, 1);
  const EBD& A2 = EKeep(P2, axis, 0);
  const EBD& B2 = EKeep(P2, axis, 1);
  const EBD det = EAdd(ESub(EMul(A0, ESub(EMul(B1, P2.W), EMul(B2, P1.W))),
                            EMul(A1, ESub(EMul(B0, P2.W), EMul(B2, P0.W)))),
                       EMul(A2, ESub(EMul(B0, P1.W), EMul(B1, P0.W))));
  const int sDet = certSign(det);
  if (sDet == 0 || sW0 == 0 || sW1 == 0 || sW2 == 0) return 0;
  return sDet * sW0 * sW1 * sW2;
}

// ===========================================================================
// INPUT-EXACT compositional evaluator (depthk lane, ported).  The CAVEAT the
// theory lanes flagged: CramerHPoint above takes ROUNDED faceN doubles, so it
// is exact-relative-to-rounded-planes.  This evaluator composes each plane's
// (n = cross, d = dot) EXACTLY from the input vertex doubles, so every
// construction and decision is exact w.r.t. the true input coordinates (the
// degree-20-in-inputs instantiation of the ONE homogeneous form).  Carried as
// exact-integer sign-magnitude bigints (value = sign*mag*2^e), no monomial
// expansion (depthk: d2-vs-Fractions 0 disagreements, EBD depth-general/sound).
// Construction only here (position + the offline-validated decisions); the
// production decision predicate reuses the same form.
struct Big {
  int sign = 0;
  long e = 0;
  std::vector<uint64_t> m;  // little-endian, no high zero limbs
};
inline void BigTrim(std::vector<uint64_t>& m) {
  while (!m.empty() && m.back() == 0) m.pop_back();
}
inline Big BigFromME(int64_t mant, long exp) {
  Big b;
  if (mant == 0) return b;
  b.sign = mant < 0 ? -1 : 1;
  b.e = exp;
  b.m.push_back((uint64_t)(mant < 0 ? -(unsigned long long)mant : mant));
  return b;
}
inline Big BigFromDouble(double d) {
  if (d == 0.0) return Big{};
  int ex;
  const double f = std::frexp(d, &ex);
  return BigFromME((int64_t)std::ldexp(f, 53), ex - 53);
}
inline int BigCmpMag(const std::vector<uint64_t>& a,
                     const std::vector<uint64_t>& b) {
  if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
  for (size_t i = a.size(); i-- > 0;)
    if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
  return 0;
}
inline std::vector<uint64_t> BigAddMag(const std::vector<uint64_t>& a,
                                       const std::vector<uint64_t>& b) {
  std::vector<uint64_t> r;
  const size_t n = std::max(a.size(), b.size());
  r.resize(n + 1, 0);
  unsigned __int128 carry = 0;
  for (size_t i = 0; i < n; ++i) {
    unsigned __int128 s = carry;
    if (i < a.size()) s += a[i];
    if (i < b.size()) s += b[i];
    r[i] = (uint64_t)s;
    carry = s >> 64;
  }
  r[n] = (uint64_t)carry;
  BigTrim(r);
  return r;
}
inline std::vector<uint64_t> BigSubMag(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b) {  // a>=b
  std::vector<uint64_t> r(a.size(), 0);
  __int128 borrow = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    __int128 s = (__int128)a[i] - borrow - (i < b.size() ? b[i] : 0);
    if (s < 0) {
      s += ((__int128)1 << 64);
      borrow = 1;
    } else
      borrow = 0;
    r[i] = (uint64_t)s;
  }
  BigTrim(r);
  return r;
}
inline std::vector<uint64_t> BigShl(const std::vector<uint64_t>& m, long k) {
  if (m.empty() || k == 0) return m;
  const long limbShift = k / 64, bitShift = k % 64;
  std::vector<uint64_t> r(m.size() + limbShift + 1, 0);
  for (size_t i = 0; i < m.size(); ++i) {
    r[i + limbShift] |= m[i] << bitShift;
    if (bitShift) r[i + limbShift + 1] |= m[i] >> (64 - bitShift);
  }
  BigTrim(r);
  return r;
}
inline std::vector<uint64_t> BigMulMag(const std::vector<uint64_t>& a,
                                       const std::vector<uint64_t>& b) {
  if (a.empty() || b.empty()) return {};
  std::vector<uint64_t> r(a.size() + b.size(), 0);
  for (size_t i = 0; i < a.size(); ++i) {
    unsigned __int128 carry = 0;
    for (size_t j = 0; j < b.size(); ++j) {
      unsigned __int128 s = (unsigned __int128)a[i] * b[j] + r[i + j] + carry;
      r[i + j] = (uint64_t)s;
      carry = s >> 64;
    }
    r[i + b.size()] += (uint64_t)carry;
  }
  BigTrim(r);
  return r;
}
inline Big BigMul(const Big& a, const Big& b) {
  Big r;
  if (a.sign == 0 || b.sign == 0) return r;
  r.sign = a.sign * b.sign;
  r.e = a.e + b.e;
  r.m = BigMulMag(a.m, b.m);
  if (r.m.empty()) r.sign = 0;
  return r;
}
inline Big BigAddSub(const Big& a, const Big& b, bool sub) {
  const int bsign = sub ? -b.sign : b.sign;
  if (a.sign == 0) {
    Big r = b;
    r.sign = bsign;
    return r;
  }
  if (bsign == 0) return a;
  const long emin = std::min(a.e, b.e);
  std::vector<uint64_t> ma = BigShl(a.m, a.e - emin);
  std::vector<uint64_t> mb = BigShl(b.m, b.e - emin);
  Big r;
  r.e = emin;
  if (a.sign == bsign) {
    r.m = BigAddMag(ma, mb);
    r.sign = a.sign;
  } else {
    const int c = BigCmpMag(ma, mb);
    if (c == 0) return Big{};
    if (c > 0) {
      r.m = BigSubMag(ma, mb);
      r.sign = a.sign;
    } else {
      r.m = BigSubMag(mb, ma);
      r.sign = bsign;
    }
  }
  if (r.m.empty()) r.sign = 0;
  return r;
}
inline Big BigAdd(const Big& a, const Big& b) { return BigAddSub(a, b, false); }
inline Big BigSub(const Big& a, const Big& b) { return BigAddSub(a, b, true); }
inline int BigSign(const Big& a) { return a.sign; }
// Leading-bits conversion to long double (value = sign*mag*2^e); exact enough
// for a position quotient (result rounds to a double).
inline long double BigToLD(const Big& b) {
  if (b.sign == 0) return 0.0L;
  const int n = (int)b.m.size();
  long double v = (long double)b.m[n - 1];
  long lowExp = b.e + 64L * (n - 1);
  if (n >= 2) {
    v = v * 18446744073709551616.0L + (long double)b.m[n - 2];  // *2^64
    lowExp = b.e + 64L * (n - 2);
  }
  return (b.sign < 0 ? -1.0L : 1.0L) * std::ldexp(v, (int)lowExp);
}
// A homogeneous point as 4 bigints.
struct BigHPoint {
  Big X, Y, Z, W;
};
inline BigHPoint TrivialBigHPoint(const vec3& v) {
  return {BigFromDouble(v.x), BigFromDouble(v.y), BigFromDouble(v.z),
          BigFromME(1, 0)};
}
// INPUT-EXACT plane triple point: each plane's n=cross(v1-v0,v2-v0),
// d=dot(n,v0) composed exactly as Big from the three input vertices, then the
// Cramer point (matches Intersect3Planes: c12=cross(g,h), W=nF.c12,
// X=dF*c12+dG*c20+dH*c01).  Degree 6/7 in inputs (untruncated).
inline BigHPoint CramerBigHPointIX(const vec3 tf[3], const vec3 tg[3],
                                   const vec3 th[3]) {
  auto planeBig = [](const vec3 t[3], Big& nx, Big& ny, Big& nz, Big& d) {
    const Big ax = BigSub(BigFromDouble(t[1].x), BigFromDouble(t[0].x));
    const Big ay = BigSub(BigFromDouble(t[1].y), BigFromDouble(t[0].y));
    const Big az = BigSub(BigFromDouble(t[1].z), BigFromDouble(t[0].z));
    const Big bx = BigSub(BigFromDouble(t[2].x), BigFromDouble(t[0].x));
    const Big by = BigSub(BigFromDouble(t[2].y), BigFromDouble(t[0].y));
    const Big bz = BigSub(BigFromDouble(t[2].z), BigFromDouble(t[0].z));
    nx = BigSub(BigMul(ay, bz), BigMul(az, by));
    ny = BigSub(BigMul(az, bx), BigMul(ax, bz));
    nz = BigSub(BigMul(ax, by), BigMul(ay, bx));
    d = BigAdd(BigAdd(BigMul(nx, BigFromDouble(t[0].x)),
                      BigMul(ny, BigFromDouble(t[0].y))),
               BigMul(nz, BigFromDouble(t[0].z)));
  };
  Big Fx, Fy, Fz, Fd, Gx, Gy, Gz, Gd, Hx, Hy, Hz, Hd;
  planeBig(tf, Fx, Fy, Fz, Fd);
  planeBig(tg, Gx, Gy, Gz, Gd);
  planeBig(th, Hx, Hy, Hz, Hd);
  auto cross = [](const Big& ax, const Big& ay, const Big& az, const Big& bx,
                  const Big& by, const Big& bz, Big& cx, Big& cy, Big& cz) {
    cx = BigSub(BigMul(ay, bz), BigMul(az, by));
    cy = BigSub(BigMul(az, bx), BigMul(ax, bz));
    cz = BigSub(BigMul(ax, by), BigMul(ay, bx));
  };
  Big c12x, c12y, c12z, c20x, c20y, c20z, c01x, c01y, c01z;
  cross(Gx, Gy, Gz, Hx, Hy, Hz, c12x, c12y, c12z);
  cross(Hx, Hy, Hz, Fx, Fy, Fz, c20x, c20y, c20z);
  cross(Fx, Fy, Fz, Gx, Gy, Gz, c01x, c01y, c01z);
  BigHPoint p;
  p.W = BigAdd(BigAdd(BigMul(Fx, c12x), BigMul(Fy, c12y)), BigMul(Fz, c12z));
  p.X = BigAdd(BigAdd(BigMul(Fd, c12x), BigMul(Gd, c20x)), BigMul(Hd, c01x));
  p.Y = BigAdd(BigAdd(BigMul(Fd, c12y), BigMul(Gd, c20y)), BigMul(Hd, c01y));
  p.Z = BigAdd(BigAdd(BigMul(Fd, c12z), BigMul(Gd, c20z)), BigMul(Hd, c01z));
  return p;
}
// INPUT-EXACT point-pair /\ plane point (P3): plane g composed exactly from its
// rep face's three input vertices, then W = ng.(w-u), coord = u*W +
// (dg-ng.u)*(w-u) over the two input edge vertices u,w.  Degree 4/3 in inputs
// (< the deg-6/7 triple), the maxdepth P3 instantiation.  Same accumulator
// (adaptive Big).
inline BigHPoint SegPlaneBigHPoint(const vec3& u, const vec3& w,
                                   const vec3 tg[3]) {
  const Big ax = BigSub(BigFromDouble(tg[1].x), BigFromDouble(tg[0].x));
  const Big ay = BigSub(BigFromDouble(tg[1].y), BigFromDouble(tg[0].y));
  const Big az = BigSub(BigFromDouble(tg[1].z), BigFromDouble(tg[0].z));
  const Big bx = BigSub(BigFromDouble(tg[2].x), BigFromDouble(tg[0].x));
  const Big by = BigSub(BigFromDouble(tg[2].y), BigFromDouble(tg[0].y));
  const Big bz = BigSub(BigFromDouble(tg[2].z), BigFromDouble(tg[0].z));
  const Big ngx = BigSub(BigMul(ay, bz), BigMul(az, by));
  const Big ngy = BigSub(BigMul(az, bx), BigMul(ax, bz));
  const Big ngz = BigSub(BigMul(ax, by), BigMul(ay, bx));
  const Big dg = BigAdd(BigAdd(BigMul(ngx, BigFromDouble(tg[0].x)),
                               BigMul(ngy, BigFromDouble(tg[0].y))),
                        BigMul(ngz, BigFromDouble(tg[0].z)));
  const Big ux = BigFromDouble(u.x), uy = BigFromDouble(u.y),
            uz = BigFromDouble(u.z);
  const Big ex = BigSub(BigFromDouble(w.x), ux),
            ey = BigSub(BigFromDouble(w.y), uy),
            ez = BigSub(BigFromDouble(w.z), uz);  // w - u
  const Big W =
      BigAdd(BigAdd(BigMul(ngx, ex), BigMul(ngy, ey)), BigMul(ngz, ez));
  const Big ngu =
      BigAdd(BigAdd(BigMul(ngx, ux), BigMul(ngy, uy)), BigMul(ngz, uz));
  const Big num = BigSub(dg, ngu);
  BigHPoint p;
  p.W = W;
  p.X = BigAdd(BigMul(ux, W), BigMul(num, ex));
  p.Y = BigAdd(BigMul(uy, W), BigMul(num, ey));
  p.Z = BigAdd(BigMul(uz, W), BigMul(num, ez));
  return p;
}
// Homogeneous X/W,Y/W,Z/W -> committed double position (input-exact position).
inline vec3 BigHPointToPos(const BigHPoint& p) {
  const long double w = BigToLD(p.W);
  return {(double)(BigToLD(p.X) / w), (double)(BigToLD(p.Y) / w),
          (double)(BigToLD(p.Z) / w)};
}
// axis-dropped 2D homogeneous coords (A,B,W) for the in-face orient2d.
struct BigP2 {
  Big A, B, W;
};
inline BigP2 BigExtract2D(const BigHPoint& p, int axis) {
  const int a0 = (axis == 0) ? 1 : 0;
  const int a1 = (axis == 2) ? 1 : 2;
  auto coord = [&](int a) { return a == 0 ? p.X : (a == 1 ? p.Y : p.Z); };
  return {coord(a0), coord(a1), p.W};
}
// INPUT-EXACT orient2d(P0,P1,P2) = sign(det[[A,B,W]...]) * sign(prod W).  The
// degree-20-in-inputs instantiation (operands degree 6/7); 0 iff exactly
// collinear in the face.
inline int BigOrient2D(const BigP2& p0, const BigP2& p1, const BigP2& p2) {
  const Big det = BigAdd(
      BigSub(BigMul(p0.A, BigSub(BigMul(p1.B, p2.W), BigMul(p2.B, p1.W))),
             BigMul(p1.A, BigSub(BigMul(p0.B, p2.W), BigMul(p2.B, p0.W)))),
      BigMul(p2.A, BigSub(BigMul(p0.B, p1.W), BigMul(p1.B, p0.W))));
  const int sW = BigSign(p0.W) * BigSign(p1.W) * BigSign(p2.W);
  return BigSign(det) * sW;
}
}  // namespace sos

// orient3d on double INPUT coords through the Shewchuk STATIC error-bound
// filter (o3derrboundA = (7 + 56u)u, u = 2^-53 - the bound constant is per
// Shewchuk's error analysis; NO kernel of his is in the build, only the
// constant).  Returns the CERTIFIED sign (+/-1) when |det| exceeds the
// permanent-scaled error bound; returns 0 when the filtered sign is UNCERTAIN
// (sub-bound or exact-zero).  The 0 case is no longer a fail-closed boundary:
// it routes to the micro exact tie-test (Orient3DExactSign, filter-0-only) and
// - on a genuine exact zero - the single-global SoS (Orient3DSoS), per the
// stage-6 owner contract.  On the corpus single-shell self-intersectors this
// filter certifies every enumeration predicate (v5b-r4: siA/siB 100% certified,
// zero fallback), so the exact kernel below is never touched on the certified
// fast path.
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
  return 0;  // uncertain -> Orient3DExactSign, then SoS (Orient3DSoS)
}

// The micro exact tie-test: the EXACT orient3d sign, 0 iff the four points are
// exactly coplanar.  TOTAL for every finite-double input.
// TRIPWIRE (owner contract, docs/Regularize3D.md open list): this is the ONE
// blessed exact predicate FORM - the sign of a homogeneous orientation
// determinant corrected by sign(prod W_i), summed on the ONE adaptive-width
// integer accumulator (sos::SumSignN).  It has exactly TWO instantiations:
//  (1) INPUT-POINT ORIENT3D (degree 3): four points with W==1 (an input double
//      point == the intersection of its three trivial axis planes, denominator
//      1).  The W column is the literal ones column and the weight-product sign
//      is +1, so this instantiation (sos::HomogOrient3DSign<true>) is
//      BIT-IDENTICAL to the historical orient3d (same 24-permutation,
//      <=3-factor width-3 accumulation) - every existing caller unchanged: the
//      EdgePiercesTriSoS edge-in-plane guard; the winding-crossing escalation
//      (WindCrossTri, shared by the O(nTri) walk, the winding broadphase, the
//      once-per-component seed-sign precompute, and the RecordSeams
//      phantom-seam guard's strict-interior pierce test); the junction
//      registry's input-vertex-on-edge arm (InputVertexStrictlyOnEdge, via
//      ExactOrient2DDrop); ExactSegProperCross/ExactOrient2DDrop (exact
//      in-plane orient2d on the ROUNDED seam endpoint doubles, drop the
//      dominant normal axis), now the E1_OFF/E1_MEASURE seam-crossing levers
//      only; the tie cascade Orient3DSoS (whose SoS K==0 group IS this exact
//      sign); and the test probe.
//  (2) CONSTRUCTED-POINT ORIENT2D (degree 9): three in-face crossing points,
//      each the Cramer intersection of a plane triple {F,g,h} (homogeneous
//      X,Y,W are 3x3 determinants of the plane coefficients).  This is the SAME
//      form at degree 9 - sos::HomogOrient2DFilter (the construction-aware EBD
//      running forward-error filter; the naive final-determinant permanent is
//      UNSOUND, it collapses in the near-parallel wedge regime, so the
//      certified bound propagates the error THROUGH the construction)
//      filter-first, then on a filter-0 the exact sos::HomogOrient2DExact fires
//      (the degree-9 monomial determinant summed on sos::SumSignN<8>).  An
//      exact zero is a GENUINE coincidence (concurrent triple points / aliased
//      crossing) routed to the level-0 incidence path (nomerge), never a
//      perturbation - the new site needs no SoS (SoS stays input-point-scoped).
//      Its production caller is ExactSeamsCross (via HPointStrictlyInTri): the
//      EXACT segment x segment straddle deciding the seam-crossing
//      enumeration's existence decision - the constructed crossing X={f,g,h}
//      strictly interior to all three seam triangles f, g, h (X within both
//      seams' symbolic extent), which carries the seam SEGMENT extent, not just
//      the line.
// A THIRD instantiation of a NEW DEGREE (or any new constructed-point form) is
// an OWNER DECISION - it widens the accumulator's proven totality bound and
// needs a new per-degree filter constant; never add one silently.  Vendoring
// rule intact: if an exact primitive OUTSIDE this one form is ever needed,
// VENDOR Shewchuk's public-domain predicates.c - do NOT rebuild expansion
// arithmetic piecemeal.  The certified fast path never touches the exact
// kernel.
inline int Orient3DExactSign(const vec3& a, const vec3& b, const vec3& c,
                             const vec3& d) {
  // The w==1 instantiation of the ONE form: an input point is the intersection
  // of its three trivial axis planes (denominator 1), so W==1 and this is
  // bit-identical to the historical orient3d (harness h1/h1b: 0 diff / 1.5e7).
  const double pts[4][4] = {{a.x, a.y, a.z, 1.0},
                            {b.x, b.y, b.z, 1.0},
                            {c.x, c.y, c.z, 1.0},
                            {d.x, d.y, d.z, 1.0}};
  return sos::HomogOrient3DSign<true>(pts);
}

// The complete orient3d decision (docs/Regularize3D.md stage 6): the certified
// filter sign on the fast path; else the single-global SoS decides.  NEVER 0.
// `i*` are the four points' global vertex indices.  No pre-SoS exact shortcut:
// SoSOrient3D's e^0 (K==0) monomial group IS the exact orient3d (the same 24
// real terms, summed by the same accumulator), so when the four points are NOT
// exactly coplanar SoS returns that exact sign from its lowest-K group - a
// separate ExactSign call ahead of it would return the identical value and is
// provably redundant (a bitwise no-op, verified).  A genuine exact zero (K==0
// group sums to zero) falls through to the perturbation, as it must.
inline int Orient3DSoS(const vec3& a, const vec3& b, const vec3& c,
                       const vec3& d, int ia, int ib, int ic, int id) {
  const int s = Orient3DFilterSign(a, b, c, d);
  if (s != 0) return s;
  const double pts[4][3] = {
      {a.x, a.y, a.z}, {b.x, b.y, b.z}, {c.x, c.y, c.z}, {d.x, d.y, d.z}};
  const int idx[4] = {ia, ib, ic, id};
  return sos::SoSOrient3D(pts, idx);
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

// The SoS completion of EdgePiercesTri for the residue the FILTER refuses (-1):
// a GENUINE NON-COPLANAR transversal exact-zero tie (vertex-on-face /
// edge-on-edge), which the single-global SoS now DECIDES to a definite pierce
// (1) or non-pierce (0) - never refuses.  `i*` are the vertices' global
// indices.  The caller restricts this to the transversal residue (the coplanar
// / cluster-riser / benign ties are the fold's, handled before this is
// reached); the edge-in-plane guard below is a defensive second gate so a
// coplanar incidence can never manufacture a phantom seam.
int EdgePiercesTriSoS(const vec3& u, const vec3& v, const vec3& a,
                      const vec3& b, const vec3& c, int iu, int iv, int ia,
                      int ib, int ic) {
  const int fu = Orient3DFilterSign(a, b, c, u);
  const int fv = Orient3DFilterSign(a, b, c, v);
  if (fu == 0 && fv == 0 && Orient3DExactSign(a, b, c, u) == 0 &&
      Orient3DExactSign(a, b, c, v) == 0)
    return 0;  // edge exactly in the tri plane: the fold's, never a seam
  const int su = Orient3DSoS(a, b, c, u, ia, ib, ic, iu);
  const int sv = Orient3DSoS(a, b, c, v, ia, ib, ic, iv);
  if (su == sv) return 0;  // both same (perturbed) side -> no straddle
  const int o1 = Orient3DSoS(u, v, a, b, iu, iv, ia, ib);
  const int o2 = Orient3DSoS(u, v, b, c, iu, iv, ib, ic);
  const int o3 = Orient3DSoS(u, v, c, a, iu, iv, ic, ia);
  return (o1 == o2 && o2 == o3) ? 1 : 0;
}

// Per-triangle geometry cache + broadphase prefilters shared by the O(F^2)
// self-crossing scans (enumeration, exact-coplanar detection, near-coplanar
// snap): each triangle's three vertex positions and global ids plus its AABB,
// with the bbox-overlap and shared-vertex (self-adjacency) skips every scan
// applies before any predicate.  Selection/compare only (la::min/max, integer
// equality), so it is bit-for-bit the inline builds it replaces.  RecordSeams
// keeps its own copy because it stores tri/vid/faceN into the BuildArrangement.
struct TriSoup {
  std::vector<std::array<vec3, 3>> tri;
  std::vector<std::array<int, 3>> vid;
  std::vector<Box> box;
  explicit TriSoup(const Manifold::Impl& in) {
    const int nTri = static_cast<int>(in.NumTri());
    tri.resize(nTri);
    vid.resize(nTri);
    box.resize(nTri);
    for (int t = 0; t < nTri; ++t) {
      for (int k = 0; k < 3; ++k) {
        vid[t][k] = in.halfedge_.Start(3 * t + k);
        tri[t][k] = in.vertPos_[vid[t][k]];
      }
      box[t].min = la::min(la::min(tri[t][0], tri[t][1]), tri[t][2]);
      box[t].max = la::max(la::max(tri[t][0], tri[t][1]), tri[t][2]);
    }
  }
  bool BBoxOverlap(int i, int j) const { return box[i].DoesOverlap(box[j]); }
  bool SharesVert(int i, int j) const {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (vid[i][a] == vid[j][b]) return true;
    return false;
  }
};

// Are two triangles EXACTLY coplanar at level-0?  Each of the six cross checks
// (every vertex of one against the other's plane) is a certified orient3d
// filter sign; a 0 means the vertex is within ~1 ULP of the plane, i.e. the
// coplanarity gap is far below eps.  Returns a SIGN-derived bool (no
// constructed geometry), so it is the single implementation of the
// "filter-coplanar pair" test the exact-fold detector, the seam recorder, and
// the near-coplanar snap all share.
bool FacesFilterCoplanar(const std::array<vec3, 3>& Ti,
                         const std::array<vec3, 3>& Tj) {
  for (int k = 0; k < 3; ++k)
    if (Orient3DFilterSign(Ti[0], Ti[1], Ti[2], Tj[k]) != 0) return false;
  for (int k = 0; k < 3; ++k)
    if (Orient3DFilterSign(Tj[0], Tj[1], Tj[2], Ti[k]) != 0) return false;
  return true;
}

// Do two coplanar (or near-coplanar) triangles share positive 2D area,
// projected into Ti's plane?  Only genuinely OVERLAPPING faces need folding;
// the coplanar tiles of one flat face (an annulus, a subdivided facet) merely
// abut and must NOT cluster - a vertex strictly inside the other, or a
// properly-crossing edge pair, is the area-overlap witness (triangles are
// convex, so this is exhaustive).  Shared by the exact and near-coplanar
// cluster detectors.
bool TrianglesOverlap2D(const std::array<vec3, 3>& Ti,
                        const std::array<vec3, 3>& Tj) {
  const vec3 e1 = la::normalize(Ti[1] - Ti[0]);
  const vec3 nrm = la::cross(Ti[1] - Ti[0], Ti[2] - Ti[0]);
  const double nl = la::length(nrm);
  if (!(nl > 0.0)) return false;
  const vec3 e2 = la::cross(nrm / nl, e1);
  const vec3 o = Ti[0];
  auto pr = [&](const vec3& P) {
    return vec2(la::dot(P - o, e1), la::dot(P - o, e2));
  };
  vec2 A[3] = {pr(Ti[0]), pr(Ti[1]), pr(Ti[2])};
  vec2 B[3] = {pr(Tj[0]), pr(Tj[1]), pr(Tj[2])};
  auto strictIn = [&](const vec2& p, const vec2* T) {
    const double d0 = la::cross(T[1] - T[0], p - T[0]);
    const double d1 = la::cross(T[2] - T[1], p - T[1]);
    const double d2 = la::cross(T[0] - T[2], p - T[2]);
    const bool neg = d0 < 0 || d1 < 0 || d2 < 0;
    const bool pos = d0 > 0 || d1 > 0 || d2 > 0;
    return !(neg && pos) && d0 != 0 && d1 != 0 && d2 != 0;
  };
  for (int k = 0; k < 3; ++k)
    if (strictIn(A[k], B) || strictIn(B[k], A)) return true;
  auto proper = [&](const vec2& p1, const vec2& p2, const vec2& p3,
                    const vec2& p4) {
    const double d1 = la::cross(p2 - p1, p3 - p1),
                 d2 = la::cross(p2 - p1, p4 - p1);
    const double d3 = la::cross(p4 - p3, p1 - p3),
                 d4 = la::cross(p4 - p3, p2 - p3);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 &&
           d2 != 0 && d3 != 0 && d4 != 0;
  };
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      if (proper(A[a], A[(a + 1) % 3], B[b], B[(b + 1) % 3])) return true;
  return false;
}

// EXACT-COPLANAR FOLD, cluster detection (docs/Regularize3D.md coplanar axis).
// Union non-self-adjacent, bbox-overlapping faces that are EXACTLY coplanar -
// every vertex of each lies on the other's plane, decided by the level-0
// orient3d filter (all six cross-checks certified 0).  The filter returns 0
// only when the vertex is within ~1 ULP (relative) of the plane, i.e. the
// coplanarity gap is far below eps (the machine weld radius), so projecting the
// cluster onto one plane is eps-valid.  The NEAR-coplanar thin band (gap above
// the filter's error bound but below eps) has a NONZERO filter sign and is NOT
// clustered here - it is PLANARIZED upstream by SnapNearCoplanarClusters (stage
// 5) so that by the time this exact detector runs its clusters are exactly
// coplanar again.  Returns a per-face cluster id, or -1 for a face in no
// multi-face coplanar cluster (the ordinary transversal path).
std::vector<int> DetectCoplanarClusters(const Manifold::Impl& in) {
  const int nTri = static_cast<int>(in.NumTri());
  const TriSoup soup(in);
  const auto& tri = soup.tri;
  DisjointSets uf(nTri);
  bool any = false;
  // PASS 1 - SEED clusters from DISTINCT-PATCH overlaps: coplanar faces that
  // overlap in area and share NO vertex.  These are unambiguously two separate
  // coplanar sheets covering common ground (a buried plug's cap, the twin
  // composition's doubled flat face).  TrianglesOverlap2D (exact
  // strict-interior / proper-crossing) already excludes the mere ABUTMENT of
  // one flat face's tiles.  A SHARES-VERTEX overlap is deferred to pass 2.
  std::vector<std::pair<int, int>> sharedCornerOverlaps;
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (!soup.BBoxOverlap(i, j)) continue;
      if (!FacesFilterCoplanar(tri[i], tri[j]) ||
          !TrianglesOverlap2D(tri[i], tri[j]))
        continue;
      if (soup.SharesVert(i, j)) {
        sharedCornerOverlaps.emplace_back(i, j);
        continue;
      }
      uf.unite(i, j);
      any = true;
    }
  // PASS 2 - EXTEND a seeded cluster through SHARED-CORNER overlaps.  A
  // coplanar self-overlap group's faces meet at shared corners/edges too (the
  // GT7863 flat face's tiles that touch a distinct-patch overlap at a vertex);
  // leaving them out made the fold's cluster INCOMPLETE, so its in-plane cover
  // disagreed with the 3D winding (the self-check fired) and the partial
  // re-triangulation T-junctioned against the un-clustered coplanar neighbours
  // (an OPEN-BOUNDARY emission fail).  But a shares-vertex overlap with NO
  // distinct-patch seed is a LOCAL FOLD-BACK (an everted-spike face pair,
  // PokedCube) that the per-face winding rule owns - do NOT seed a cluster from
  // it.  So only join a shared-corner pair when one side already belongs to a
  // seeded cluster; iterate to a fixpoint so a chain extends fully.  (This
  // mirrors the shares-vertex arrangement-completeness fix reg3d-wjump made in
  // RecordSeams; reg3d-7863c1.)
  if (any) {
    std::vector<int> setSize(nTri, 0);
    for (int f = 0; f < nTri; ++f) ++setSize[static_cast<int>(uf.find(f))];
    auto inCluster = [&](int f) {
      return setSize[static_cast<int>(uf.find(f))] > 1;
    };
    bool changed = true;
    while (changed) {
      changed = false;
      for (const auto& [i, j] : sharedCornerOverlaps) {
        if (uf.same(i, j)) continue;
        if (!inCluster(i) && !inCluster(j)) continue;
        const int ri = static_cast<int>(uf.find(i)),
                  rj = static_cast<int>(uf.find(j));
        const int merged = setSize[ri] + setSize[rj];
        uf.unite(i, j);
        setSize[static_cast<int>(uf.find(i))] = merged;
        changed = true;
      }
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

// The gate's cheap existence probe: DetectCoplanarClusters' pass-1 seed test
// with an early-out - the FIRST distinct-patch (non-self-adjacent) coplanar
// overlap makes the component dirty, without building the DisjointSets, the
// shared-corner extension, or the id renumber.  Answer-identical to
// "any face id >= 0 from DetectCoplanarClusters" (proven above: pass 2 only
// extends pass-1 seeds, so a face is clustered iff pass 1 seeded).
bool HasCoplanarOverlap(const Manifold::Impl& in) {
  const int nTri = static_cast<int>(in.NumTri());
  const TriSoup soup(in);
  const auto& tri = soup.tri;
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (!soup.BBoxOverlap(i, j)) continue;
      if (soup.SharesVert(i, j))
        continue;  // distinct-patch seeds only (pass 1)
      if (FacesFilterCoplanar(tri[i], tri[j]) &&
          TrianglesOverlap2D(tri[i], tri[j]))
        return true;
    }
  return false;
}

// Coupled soup winding w_S(p): the signed count of oriented-face crossings on
// the ray p->seed, every crossing decided by level-0 orient3d through the
// static filter (the Winding03 discipline, boolean3.cpp:388).  The delta per
// crossed face is sign(dot(seed-p, n_f)) on the input normal - a +/-1 integer,
// FP-safe by construction.  A filter-uncertain deciding predicate ESCALATES to
// the exact tie-test (below), which decides the near-tangent graze the filter
// refuses (reg3d-c2bx: GT7081's near-coplanar shallow-dihedral faces).  Returns
// nullopt only when the escalation finds a GENUINE exact-zero tie (the probe
// grazes a vertex/edge/plane exactly): there the caller re-seeds or fails
// closed, since the soup winding is single-valued only OFF the surface.
// One triangle's oriented crossing contribution for the winding ray p->seed,
// factored so the O(nTri) walk and the winding broadphase (below) share ONE
// exact predicate chain (the "one predicate / one implementation" discipline).
// TRIPWIRE caller (winding probe, docs/Regularize3D.md open list): a filter-0
// plane-side tie is NOT a fail-closed boundary - it ESCALATES to the exact
// tie-test (filter-first: exact fires only on filter-0).  A near-tangent
// shallow-dihedral face grazes the static filter's uncertainty band while the
// exact kernel decides the constructed probe DECIDABLY off the plane
// (reg3d-c2b: GT7081's minGap is a few eps, exactly ONE such face per shell).
// The probe is a constructed double the exact kernel reads verbatim, so no SoS
// / vertex index is needed.  Returns false on a GENUINE exact-zero tie (p/seed
// lie ON a face plane, or the segment grazes an edge/vertex exactly) - the soup
// winding is single-valued only OFF the surface, so the caller re-seeds / fails
// closed; on true, `delta` is the signed crossing (0 none, +/-1).
// dbCached == kWindDbLive: compute the seed's plane-side sign live (the walk).
// Otherwise dbCached IS that sign (already filter-then-exact escalated) - the
// per-seed precompute (PrecomputeSeedSign): for a FIXED seed the plane-side
// sign is pure per triangle, so it is computed ONCE instead of once per
// clean-face query.  The looked-up value is bit-identical to the live one (same
// predicate chain), so the winding is unchanged; a stored 0 reproduces the live
// genuine tie.
inline constexpr int kWindDbLive = 2;  // not a valid orient3d sign (-1/0/+1)
inline bool WindCrossTri(const vec3& a, const vec3& b, const vec3& c,
                         const vec3& p, const vec3& seed, int& delta,
                         int dbCached = kWindDbLive) {
  delta = 0;
  int da = Orient3DFilterSign(a, b, c, p);
  int db =
      (dbCached == kWindDbLive) ? Orient3DFilterSign(a, b, c, seed) : dbCached;
  if (da == 0) da = Orient3DExactSign(a, b, c, p);
  if (db == 0 && dbCached == kWindDbLive) db = Orient3DExactSign(a, b, c, seed);
  if (da == 0 || db == 0) return false;
  if (da == db) return true;  // p and seed on the same side of the plane
  int o1 = Orient3DFilterSign(p, seed, a, b);
  int o2 = Orient3DFilterSign(p, seed, b, c);
  int o3 = Orient3DFilterSign(p, seed, c, a);
  if (o1 == 0) o1 = Orient3DExactSign(p, seed, a, b);
  if (o2 == 0) o2 = Orient3DExactSign(p, seed, b, c);
  if (o3 == 0) o3 = Orient3DExactSign(p, seed, c, a);
  if (o1 == 0 || o2 == 0 || o3 == 0) return false;
  if (o1 == o2 && o2 == o3) {
    const vec3 n = la::cross(b - a, c - a);
    delta = (la::dot(seed - p, n) > 0.0) ? 1 : -1;
  }
  return true;
}

// ---------------------------------------------------------------------------
// WINDING BROADPHASE (perf: docs/Regularize3D.md winding-query axis).  The
// clean-face winding walks EVERY triangle per query (O(nTri * cleanFaces), the
// resolve hot loop on single-component shells).  A per-component Morton
// triangle collider shrinks each query to the triangles whose AABB overlaps the
// winding SEGMENT's AABB - a PROVEN EXACT SUPERSET of the segment's crossings,
// so the result is unchanged:
//   If the segment p->seed genuinely crosses triangle T (WindCrossTri delta!=0)
//   the crossing point x lies in BOTH the segment box Box(p,seed) and T's box
//   (x is a convex combination of the segment endpoints, and of T's vertices).
//   Both boxes are the componentwise min/max of the SAME doubles the exact
//   predicates read (min/max of doubles is exact), so x witnesses the closed-
//   interval Box::DoesOverlap and T is a candidate.  No box inflation, no eps.
// The only behavior the broadphase drops is the walk's OVER-conservative
// re-seed when p or seed lies exactly on a DISTANT triangle's plane the segment
// never reaches (that triangle contributes no crossing); the winding is
// single-valued off-surface so the value from this seed already equals what the
// walk returns from the next seed - strictly more permissive, never wrong
// (validated: windDiverge=0 and crossing-superset violations=0 across the
// corpus).
struct TriWindBVH {
  Collider collider;
  Vec<int> leaf2tri;  // leaf index -> original triangle index
};

TriWindBVH BuildTriWindBVH(const std::vector<std::array<vec3, 3>>& tri,
                           const Box& bBox) {
  TriWindBVH out;
  const int nTri = static_cast<int>(tri.size());
  if (nTri == 0) return out;
  std::vector<Box> box(nTri);
  std::vector<uint32_t> morton(nTri);
  for (int t = 0; t < nTri; ++t) {
    Box bx;
    bx.Union(tri[t][0]);
    bx.Union(tri[t][1]);
    bx.Union(tri[t][2]);
    box[t] = bx;
    morton[t] =
        Collider::MortonCode((tri[t][0] + tri[t][1] + tri[t][2]) / 3.0, bBox);
  }
  // LBVH construction needs Morton-sorted leaves; the sort only orders leaves
  // (a tree-quality heuristic), never the true leaf boxes, so it cannot drop a
  // candidate.
  std::vector<int> perm(nTri);
  for (int i = 0; i < nTri; ++i) perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(),
                   [&](int x, int y) { return morton[x] < morton[y]; });
  Vec<Box> sbox(nTri);
  Vec<uint32_t> smort(nTri);
  out.leaf2tri.resize(nTri);
  for (int i = 0; i < nTri; ++i) {
    sbox[i] = box[perm[i]];
    smort[i] = morton[perm[i]];
    out.leaf2tri[i] = perm[i];
  }
  out.collider = Collider(sbox, smort);
  return out;
}

// Triangles whose AABB overlaps the winding segment p->seed (the crossing
// superset).  Single-query sequential descent (safe to call concurrently on a
// const collider).
void WindCandidates(const TriWindBVH& bvh, const vec3& p, const vec3& seed,
                    std::vector<int>& out) {
  out.clear();
  const Box qbox(p, seed);
  auto rec = [&](int, int leaf) { out.push_back(bvh.leaf2tri[leaf]); };
  auto recorder = MakeSimpleRecorder(rec);
  auto f = [&](int) { return qbox; };
  bvh.collider.Collisions<false>(recorder, f, 1, /*parallel=*/false, nullptr);
}

// Winding over the candidate triangle list (SAME exact per-triangle chain as
// the walk; only the iteration set shrinks).
std::optional<int> WindingAtCands(const std::vector<std::array<vec3, 3>>& tri,
                                  const vec3& p, const vec3& seed,
                                  const std::vector<int>& cands,
                                  const signed char* seedSign) {
  int w = 0;
  for (int t : cands) {
    int delta;
    const int db = seedSign ? static_cast<int>(seedSign[t]) : kWindDbLive;
    if (!WindCrossTri(tri[t][0], tri[t][1], tri[t][2], p, seed, delta, db))
      return std::nullopt;
    w += delta;
  }
  return w;
}

std::optional<int> WindingAt(const Manifold::Impl& in, const vec3& p,
                             const vec3& seed) {
  int w = 0;
  const int nTri = static_cast<int>(in.NumTri());
  for (int t = 0; t < nTri; ++t) {
    const vec3 a = in.vertPos_[in.halfedge_.Start(3 * t)];
    const vec3 b = in.vertPos_[in.halfedge_.Start(3 * t + 1)];
    const vec3 c = in.vertPos_[in.halfedge_.Start(3 * t + 2)];
    int delta;
    if (!WindCrossTri(a, b, c, p, seed, delta)) return std::nullopt;
    w += delta;
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
// Each crossed face projects its triangle boundary + its seam segments into the
// face plane; RemoveOverlaps2D (the Smith sweep) arranges them and the
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
// on this face's plane.  (Retention is decided by the real 3D coupled winding,
// so the crossing face's normal and endpoint-interiority markers a seam-sign
// scheme once needed are no longer carried.)
struct BuildSeam {
  vec3 p0, p1;
  // The PARTNER face whose plane cut this seam (B1 triple-point naming): a seam
  // on face f from the pair (f,other) lies on plane(f) INT plane(other).  Two
  // seams on f with distinct partner planes that cross in f's interior name a
  // 3-face triple point {f, other_1, other_2}.
  int other = -1;
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
  // B1 once-only triple points (EnumerateTriplePoints).  seamTriples[f][k] =
  // one (2D on-seam crossing position, canonical 3D triple position) per triple
  // point on seam k of face f, used to PRE-SPLIT that seam at the shared
  // vertex.  The 2D position is the exact intersection of the two crossing
  // seams in f's frame (strictly on-segment, so the split never folds back);
  // the 3D position is the ONCE-ONLY point every incident face welds onto.
  // Empty for every seam when the component has no 3-face triple point (the
  // whole corpus off openscad), so EmitSeamedFace is a byte-identical no-op.
  std::vector<std::vector<std::vector<std::pair<vec2, vec3>>>> seamTriples;
  // GLOBAL JUNCTION REGISTRY (f4-junction): every once-only arrangement vertex
  // of the component - all seam endpoints AND all triple points - deduped to a
  // canonical 3D position.  A T-junction opens the emission fan when such a
  // vertex sits strictly interior to a NEIGHBOUR / PARTNER / THIRD face's
  // emitted edge without a shared split (openscad's dominant open residue).
  // Every emit path (seamed, clean, fold) pre-splits each of its emitted edges
  // at every registry vertex strictly interior to it (level-0 on-segment test
  // in the face frame, keyed by the once-only 3D bits), so all incident faces
  // split at the IDENTICAL point and the fans close.  Empty on any component
  // whose seams never terminate interior to another edge (the whole corpus off
  // openscad), so all three paths stay byte-identical there.
  std::vector<vec3> junctions;
  // planeId[f] = face2cluster[f]>=0 ? nTri+cluster : f - the plane a face lies
  // on, coplanar clusters collapsed to one id.  Populated once in
  // ResolveComponent; the seam / wedge / provenance passes read it (hoisted
  // from three identical local recomputations).  A seam of face f from partner
  // g lies on planes {planeId[f], planeId[g]}.
  std::vector<int> planeId;
  // EX2 exact-overlay flag (exact2d lane): faces carrying a GAP-FREE wedge
  // chain (a diverging near-parallel wedge whose exact interior crossing is >
  // the coordinated-collapse envelope from every seam endpoint - the openscad
  // residue the double-precision RemoveOverlaps2D eps-merge cannot emit).  A
  // flagged face's per-face arrangement is built by the EXACT PSLG overlay
  // (EmitSeamedFace) instead of RemoveOverlaps2D.  Empty / all-zero on every
  // resolving carrier (zero gap-free chains - measured), so the overlay is a
  // byte-clean no-op off openscad.
  std::vector<char> wedgeFace;
  bool ok = true;  // false = a structural anomaly (fail closed)
};

// INPUT-EXACT basis lever (this lane): route every wedge-overlay POSITION and
// DECISION through the input-exact compositional evaluator (planes composed
// exactly from input vertex doubles) instead of the rounded faceN.  The
// rounded-plane basis misplaces near-tangent triple points by up to ~1400x eps,
// scrambling the order of clusters separated by only ~22x eps; input-exact
// places and orders them correctly (offline: exact-vs-Fractions 0
// disagreements, every openscad cluster genuinely distinct + above the merge
// radius).  Byte- clean off the wedge population (provOf/planeTri empty there
// -> never reached).
inline bool IXEnabled() {
  static const bool on = std::getenv("IX_OFF") == nullptr;  // default ON
  return on;
}
// Input-exact plane-triple point of three rep faces (input vertices).
inline sos::BigHPoint IXTripleHPoint(const BuildArrangement& A, int rf, int rg,
                                     int rh) {
  return sos::CramerBigHPointIX(A.tri[rf].data(), A.tri[rg].data(),
                                A.tri[rh].data());
}
// Input-exact in-face orient2d(P0,P1,P2), dropping the face axis.
inline int IXOrient2D(const sos::BigHPoint& p0, const sos::BigHPoint& p1,
                      const sos::BigHPoint& p2, int axis) {
  return sos::BigOrient2D(sos::BigExtract2D(p0, axis),
                          sos::BigExtract2D(p1, axis),
                          sos::BigExtract2D(p2, axis));
}

// 3D point-vs-triangle classification with an area-scaled margin - the shared
// core of RecordSeams' vertex-on-face graze test and the cap-plane seam
// endpoint test.  P is assumed on the triangle's plane; n = the raw
// (unnormalized) normal, area2 = |n|^2, margin = area2 * 1e-9.  Returns
// kDegenerate if area2 == 0, else kInside if P is within margin of all three
// directed edges, else kOutside.
enum class TriSide { kDegenerate, kInside, kOutside };
inline TriSide ClassifyPointInTri3D(const vec3& P, const vec3& t0,
                                    const vec3& t1, const vec3& t2) {
  const vec3 n = la::cross(t1 - t0, t2 - t0);
  const double area2 = la::length2(n);
  if (!(area2 > 0.0)) return TriSide::kDegenerate;
  const double margin = area2 * 1e-9;
  const double s0 = la::dot(n, la::cross(t1 - t0, P - t0));
  const double s1 = la::dot(n, la::cross(t2 - t1, P - t1));
  const double s2 = la::dot(n, la::cross(t0 - t2, P - t2));
  return (s0 >= -margin && s1 >= -margin && s2 >= -margin) ? TriSide::kInside
                                                           : TriSide::kOutside;
}

// Enumerate the self-crossing arrangement AND record each seam's canonical 3D
// segment per incident face.  `face2cluster` (from DetectCoplanarClusters)
// marks exactly-coplanar face groups: same-cluster pairs are SKIPPED here (the
// in-plane fold resolves them, not the transversal seam machinery).  The
// optional seamCountOut/boundaryTouchOut fold the test probe's level-0
// self-crossing counts onto this scan (nullptr in production).
BuildArrangement RecordSeams(const Manifold::Impl& in,
                             const std::vector<int>& face2cluster, double eps,
                             int* seamCountOut = nullptr,
                             int* boundaryTouchOut = nullptr) {
  BuildArrangement A;
  const int nTri = static_cast<int>(in.NumTri());
  A.tri.resize(nTri);
  A.vid.resize(nTri);
  A.faceN.resize(nTri);
  A.faceSeams.resize(nTri);
  A.seamed.assign(nTri, 0);
  std::vector<Box> box(nTri);
  for (int t = 0; t < nTri; ++t) {
    for (int k = 0; k < 3; ++k) {
      A.vid[t][k] = in.halfedge_.Start(3 * t + k);
      A.tri[t][k] = in.vertPos_[A.vid[t][k]];
    }
    A.faceN[t] =
        la::cross(A.tri[t][1] - A.tri[t][0], A.tri[t][2] - A.tri[t][0]);
    box[t].min = la::min(la::min(A.tri[t][0], A.tri[t][1]), A.tri[t][2]);
    box[t].max = la::max(la::max(A.tri[t][0], A.tri[t][1]), A.tri[t][2]);
  }
  // Once-only pierce cache (doc R1): the UNDIRECTED (min,max) edge key makes a
  // shared mesh edge traversed in opposite order by two adjacent faces resolve
  // to ONE bit-identical pierce point, so chaining seams meet exactly and the
  // cross-face weld cannot manufacture a twin.  A correctness-by-construction
  // backstop for the near-parallel tail, not a runtime check (measured
  // non-load-bearing on the corpus, but retained; docs/Regularize3D.md R1).
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
  auto bboxOverlap = [&](int i, int j) { return box[i].DoesOverlap(box[j]); };
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
      // Probe counters (test hook): the level-0 self-crossing classification
      // EnumerateComponent_Probe reports, folded onto THIS scan so no separate
      // enumeration pass is needed.  Non-adjacent (shares-vertex skipped, S4a)
      // bbox-overlap pairs are graded by the raw level-0 pierce test: a genuine
      // pierce (r==1) is a seam, else a filter-uncertain touch (r==-1) is a
      // boundary-touch pair.  nullptr (production) skips this entirely, so the
      // recorded arrangement is unchanged.
      if (seamCountOut && boundaryTouchOut) {
        // Index-keyed shares-vertex skip, PROBE-ONLY (inlined at its one use so
        // it does not read as a shared production helper).  The production
        // recovery and skip below key on POSITION coincidence; the probe's
        // welded synthetic meshes carry no unwelded duplicates, so index
        // equality is the right adjacency test here.
        bool sharesVertIdx = false;
        for (int a = 0; a < 3 && !sharesVertIdx; ++a)
          for (int b = 0; b < 3; ++b)
            if (A.vid[i][a] == A.vid[j][b]) sharesVertIdx = true;
        if (!sharesVertIdx) {
          const auto& Ti = A.tri[i];
          const auto& Tj = A.tri[j];
          bool genuine = false, boundary = false;
          for (int e = 0; e < 3 && !genuine; ++e) {
            const int r =
                EdgePiercesTri(Ti[e], Ti[(e + 1) % 3], Tj[0], Tj[1], Tj[2]);
            if (r == 1) genuine = true;
            if (r == -1) boundary = true;
          }
          for (int e = 0; e < 3 && !genuine; ++e) {
            const int r =
                EdgePiercesTri(Tj[e], Tj[(e + 1) % 3], Ti[0], Ti[1], Ti[2]);
            if (r == 1) genuine = true;
            if (r == -1) boundary = true;
          }
          if (genuine)
            ++*seamCountOut;
          else if (boundary)
            ++*boundaryTouchOut;
        }
      }
      // SHARES-VERTEX GENUINE-CROSSING RECOVERY (reg3d-wjump
      // decision-completion 1), keyed on POSITION coincidence (reg3d-oscad
      // REOPEN).  The shared-vertex broadphase skip assumes shared-vertex =>
      // self-adjacent => no transversal crossing, which is UNSOUND for a
      // self-intersecting soup where two faces sharing a vertex fold back and
      // cross elsewhere (the everted-corner / near-triple-point arrangement
      // incompleteness: PokedCube's spike, reg3d-c1a/s7b's "unmatched seam
      // pierce points").  An imported soup additionally carries UNWELDED
      // DUPLICATE vertices (identical position, distinct global id); two faces
      // meeting only at such a coincident corner are self-adjacent-at-a-point
      // exactly like an index-shared pair, but an index-keyed skip MISSES them
      // (measuring shares-vertex by INDEX reads no shared vertex and wrongly
      // rules the family out - it IS the shares-vertex family, keyed on
      // POSITION).  Key the skip/recovery on position coincidence (index
      // equality implies it).  Recover the pair ONLY when an edge with NEITHER
      // endpoint at a coincident-vertex position genuinely pierces the other
      // triangle's interior (a real OFF-VERTEX transversal crossing); the
      // coincident-vertex / shared-edge touches the skip correctly drops have
      // no such off-vertex pierce and stay skipped.  The SoS convention
      // (EdgePiercesTriSoS) decides the shared-vertex edge incidences when the
      // seam is assembled below (EXACT-INCIDENT TIES).  (The probe classify
      // above keys sharesVert by INDEX: its welded synthetic meshes carry no
      // unwelded duplicates, and it never runs on the production path.)
      bool recoveredSV = false;
      vec3 svPos(0.0);
      int nShared = 0;
      bool anyCoinc = false;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          if (A.tri[i][a] == A.tri[j][b]) anyCoinc = true;
      if (anyCoinc) {
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            if (A.tri[i][a] == A.tri[j][b]) {
              ++nShared;
              svPos = A.tri[i][a];
            }
        auto offVertexPierces = [&](int owner, int tgt) {
          const auto& To = A.tri[owner];
          const auto& T = A.tri[tgt];
          for (int e = 0; e < 3; ++e) {
            const vec3& u = To[e];
            const vec3& w = To[(e + 1) % 3];
            bool inc = false;
            for (int k = 0; k < 3; ++k)
              if (T[k] == u || T[k] == w)
                inc = true;  // edge at a coincident vtx
            if (inc)
              continue;  // edge incident to a coincident vertex: trivial touch
            int r = EdgePiercesTri(u, w, T[0], T[1], T[2]);
            if (r == -1)
              r = EdgePiercesTriSoS(u, w, T[0], T[1], T[2], A.vid[owner][e],
                                    A.vid[owner][(e + 1) % 3], A.vid[tgt][0],
                                    A.vid[tgt][1], A.vid[tgt][2]);
            if (r == 1) return true;
          }
          return false;
        };
        if (!offVertexPierces(i, j) && !offVertexPierces(j, i)) continue;
        recoveredSV = true;
      }
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
      int nPts = 0;
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
          // Benign iff p lies strictly OUTSIDE tgt's triangle; a degenerate
          // tgt is treated as a real tie (not benign), matching the old guard.
          return ClassifyPointInTri3D(p, A.tri[tgt][0], A.tri[tgt][1],
                                      A.tri[tgt][2]) == TriSide::kOutside;
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
      // Pair coplanarity (filter): a coplanar pair's exact-zero ties belong to
      // the in-plane FOLD, never the SoS transversal path (the s3/s4adj sliver
      // rail).  Only a NON-coplanar transversal residue is SoS-decided.
      const bool pairCoplanar = FacesFilterCoplanar(T0, T1);
      // Record the crossing points of `owner`'s edges through `tgt`.  A FILTER-
      // certified pierce (r==1) records directly.  A filter-refused (-1)
      // GENUINE transversal exact-zero tie - not a cluster riser, not a benign
      // graze, and the pair is non-coplanar - is DECIDED by the single-global
      // SoS (EdgePiercesTriSoS): the stage-6 completion of the vertex-on-face /
      // edge-on-edge tie family.  Coplanar / cluster-riser / benign ties record
      // nothing (the fold or the valid-manifold structure owns them).
      auto recordEdge = [&](int owner, int e, int tgt) {
        const vec3& u = A.tri[owner][e];
        const vec3& w = A.tri[owner][(e + 1) % 3];
        const auto& T = A.tri[tgt];
        const int r = EdgePiercesTri(u, w, T[0], T[1], T[2]);
        bool hit = (r == 1);
        if (r == -1 && !pairCoplanar && !edgeInClusterPlane(owner, e) &&
            !benignInPlane(owner, e, tgt)) {
          // Genuine non-coplanar transversal exact-zero residue: the single-
          // global SoS (stage 6) DECIDES it - pierce (1) or no-pierce (0),
          // never a refusal (EdgePiercesTriSoS is TOTAL: every leg is
          // Orient3DSoS, which never returns 0).
          const int sp =
              EdgePiercesTriSoS(u, w, T[0], T[1], T[2], A.vid[owner][e],
                                A.vid[owner][(e + 1) % 3], A.vid[tgt][0],
                                A.vid[tgt][1], A.vid[tgt][2]);
          if (sp == 1) hit = true;
        }
        // CAP-PLANE SEAM ENDPOINT (coplanar/transversal junction completion):
        // when BOTH faces of a wall-wall pair are non-cluster (neither is a
        // folded cap) but `owner`'s edge LIES in a coplanar cap cluster plane
        // and PROPERLY CROSSES `tgt`'s plane (both endpoints strictly off it,
        // on opposite sides) at a point on/inside `tgt`'s triangle, that
        // crossing is a genuine transversal seam endpoint sitting ON the cap
        // fold's in-plane arrangement (an overlap-corner vertex the fold also
        // emits). edgeInClusterPlane's blanket suppression would drop it,
        // truncating the seam (nPts==1) at the junction; recording it gives the
        // endpoint the fold arrangement's identity so the seamed wall and the
        // folded cap weld shut at the reentrant corner.  The symmetric
        // double-pierce (both walls' cap edges meet here) is deduped below; any
        // mis-/over-recovery only ever makes nPts!=2 -> fail closed, never a
        // wrong resolve.
        if (!hit && r == -1 && !pairCoplanar && face2cluster[i] < 0 &&
            face2cluster[j] < 0 && edgeInClusterPlane(owner, e)) {
          const int su = Orient3DFilterSign(T[0], T[1], T[2], u);
          const int sv = Orient3DFilterSign(T[0], T[1], T[2], w);
          if (su != 0 && sv != 0 && su != sv) {
            const vec3 P = SegPlanePoint(u, w, T[0], T[1], T[2]);
            if (ClassifyPointInTri3D(P, T[0], T[1], T[2]) == TriSide::kInside)
              hit = true;
          }
        }
        if (hit && nPts < 4)
          pts[nPts++] = pierce(A.vid[owner][e], A.vid[owner][(e + 1) % 3], tgt);
      };
      for (int e = 0; e < 3; ++e) recordEdge(i, e, j);
      for (int e = 0; e < 3; ++e) recordEdge(j, e, i);
      // Dedup endpoints within the weld radius: a cap-plane reentrant junction
      // is pierced by BOTH walls' cap edges (the symmetric corner incidence)
      // and so is recorded twice; the emission weld would merge them anyway.
      // Genuine distinct seam endpoints are far more than eps apart, so an
      // ordinary seam is untouched (bitwise).
      const int nPtsPre = nPts;
      for (int a = 0; a + 1 < nPts; ++a)
        for (int b = nPts - 1; b > a; --b)
          if (la::length(pts[a] - pts[b]) <= eps) {
            pts[b] = pts[nPts - 1];
            --nPts;
          }
      if (nPts == 0) continue;  // no genuine crossing
      // Recovered shares-vertex crossing: two triangles sharing exactly one
      // vertex V and crossing transversally intersect along [V, P], where P is
      // the recorded off-vertex crossing and V is the shared corner (on both
      // planes, on both boundaries).  The edges incident to V only TOUCH at V
      // (not a transversal pierce), so V is not collected; supply it as the
      // second seam endpoint (on both faces' boundary).
      if (nPts == 1 && recoveredSV && nShared == 1 &&
          la::length(pts[0] - svPos) > eps) {
        pts[1] = svPos;
        nPts = 2;
      }
      if (nPts != 2) {
        // PHANTOM-SEAM GUARD (reg3d-oscad REOPEN).  A lone / absent endpoint is
        // a genuine seam only where the pair CROSSES TRANSVERSALLY.  Exact
        // rational reconstruction of the imported-soup residue showed the
        // dominant nPts==1 truncations are MEASURE-ZERO contacts (coincident-
        // vertex duplicates, vertex-on-edge / vertex-on-face T-junctions,
        // collinear edge-on-edge overlaps) whose exact tri-tri intersection is
        // a POINT or a boundary segment - the filter/SoS recorded a PHANTOM
        // endpoint but a valid arrangement has NO seam there (proven exactly:
        // not one has a clean off-plane edge piercing the other's strict
        // interior).
        //
        // SUB-EPS SEAM COLLAPSE (witness theorem): if two or more endpoints
        // were collected but the dedup merged them under the weld radius to a
        // single point (nPtsPre>=2 -> nPts==1), the seam degenerates to a point
        // below eps and its endpoints weld - the cancel/collapse case.  It
        // contributes NO split; skip it.  Exact (no re-derivation): only a pair
        // that actually recorded a second endpoint within eps collapses, so a
        // genuine >eps seam cannot be dropped here.
        if (nPts == 1 && nPtsPre >= 2) continue;
        // MEASURE-ZERO CONTACT: a seam is real only where an edge CLEANLY
        // pierces the other's STRICT interior - both endpoints strictly off the
        // plane on OPPOSITE sides AND the crossing strictly inside the triangle
        // (all three edge-edge orientations one nonzero sign).  This is the
        // EXACT, NON-PERTURBING completion of EdgePiercesTri: a filter-REFUSED
        // sign is completed by the EXACT predicate (Orient3DExactSign, filter-
        // first), and ANY exact-zero - an endpoint exactly ON the plane
        // (vertex-on-face / edge T-junction) or a crossing exactly ON the
        // triangle boundary (collinear edge-on-edge graze) - is a MEASURE-ZERO
        // contact, NOT a clean pierce, so that edge is skipped.  It must NOT be
        // handed to the SoS convention (EdgePiercesTriSoS): SoS answers "which
        // way under perturbation", which manufactures a PHANTOM pierce out of
        // exactly these measure-zero contacts (the openscad soup is dense with
        // them; the oscad-reopen exact reconstruction proved NONE has a clean
        // strict-interior pierce).  A genuine near-tangent crossing the filter
        // cannot certify IS caught (exact-completed to a strict pierce) and, at
        // nPts!=2, fails closed - never a silent drop.
        // A strict-interior pierce IS a DECIDED WindCrossTri with a nonzero
        // crossing delta, so this rides the ONE blessed WindCrossTri escalation
        // chain instead of re-rolling the predicate.  For owner edge (u,w) vs
        // triangle Tt, WindCrossTri(Tt0,Tt1,Tt2, u, w, delta) runs the SAME
        // plane-side (da/db) and edge-edge (o1/o2/o3) orient chain in the SAME
        // order with the SAME filter-then-exact completion.  It returns false
        // on any exact-zero (endpoint on-plane / crossing on the tri boundary =
        // a measure-zero contact) and on same-side (no straddle); delta!=0
        // exactly on the strict-interior crossing (da!=db forces dot(w-u,n)!=0,
        // so the sign is +/-1).  The extra cross/dot only sets delta's sign,
        // which the delta!=0 test collapses - bitwise-identical verdict, no new
        // exact call.
        auto cleanPierce = [&](int owner, int tgt) {
          const auto& To = A.tri[owner];
          const auto& Tt = A.tri[tgt];
          for (int e = 0; e < 3; ++e) {
            int delta;
            if (WindCrossTri(Tt[0], Tt[1], Tt[2], To[e], To[(e + 1) % 3],
                             delta) &&
                delta != 0)
              return true;  // strict-interior pierce
          }
          return false;
        };
        if (!cleanPierce(i, j) && !cleanPierce(j, i))
          continue;  // certified (filter/exact) no clean pierce: measure-zero
        // A genuine transversal crossing whose two eps-separated endpoints did
        // not both record is an honestly-open degenerate incidence - fail
        // closed.
        A.ok = false;
        continue;
      }
      A.faceSeams[i].push_back({pts[0], pts[1], j});
      A.faceSeams[j].push_back({pts[0], pts[1], i});
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
        area += la::cross(p, q);
      }
      if (area > 0.0)
        cells.push_back(std::move(loop));
      else if (holes && loop.size() >= 3)
        holes->push_back(std::move(loop));
    }
  return true;
}

// Dominant axis of a face normal (index of the largest-magnitude component).
inline int DominantAxis(const vec3& n) {
  const double ax = std::abs(n.x), ay = std::abs(n.y), az = std::abs(n.z);
  return (ax >= ay && ax >= az) ? 0 : (ay >= az ? 1 : 2);
}

// In-plane 2D flattening + canonical-vertex table shared by the per-face
// arrangement drivers (EmitSeamedFace, FoldCoplanarClusters, the clean-face
// junction split).  Flattening is an EXACT axis-drop: the 2D projection SELECTS
// the two non-dominant coordinates (pure coordinate selection - no arithmetic,
// no rounding), so every orientation/merge predicate inside RemoveOverlaps2D
// runs on exact input/construction doubles at level 0.  The dropped axis is the
// dominant component of the face normal, so |nHat[axis]| is the largest
// component of a unit vector (>= 1/sqrt(3)) and the plane is never edge-on to
// the projection.  The surviving pair is (axis+1, axis+2) cyclically so 2D-CCW
// maps to +axis, SWAPPED when the dropped normal component is negative so 2D-
// CCW still maps to +nHat (the parity fix).  `lift` inverts the drop for a
// 2D-BORN point by solving the plane equation for the dropped coordinate.  Owns
// the dedup-by-3D-bit-pattern insertion (shared endpoints - chain junctions,
// corners - collapse to one input vertex with one projection).
struct AxisDropFrame {
  int axis = 2;         // dropped (dominant normal) axis
  bool flip = false;    // dropped normal component < 0: swap the 2D pair
  vec3 nHat;            // unit face normal (the winding-probe offset direction)
  double planeD = 0.0;  // nHat . (an in-plane point), for the lift
  std::vector<vec2> verts2;
  std::vector<vec3> canon3;
  std::map<std::tuple<double, double, double>, int> vidx;
  // EXACT axis-drop: select the (axis+1, axis+2) coordinate pair, swapped when
  // the dropped normal component is negative so 2D-CCW maps to +nHat.  Affine
  // (a linear projection), so on-line feet and param sorts stay consistent.
  vec2 proj(const vec3& P) const {
    const double c1 = P[(axis + 1) % 3], c2 = P[(axis + 2) % 3];
    return flip ? vec2(c2, c1) : vec2(c1, c2);
  }
  // Lift a 2D-BORN point (a RemoveOverlaps2D crossing / an interior probe with
  // no 3D preimage) back onto the plane: undo the pair swap, then solve
  // nHat . P = planeD for the dropped (dominant, so nonzero) coordinate.
  vec3 lift(const vec2& q) const {
    const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
    const double c1 = flip ? q.y : q.x, c2 = flip ? q.x : q.y;
    vec3 P;
    P[a1] = c1;
    P[a2] = c2;
    P[axis] = (planeD - nHat[a1] * c1 - nHat[a2] * c2) / nHat[axis];
    return P;
  }
  int add(const vec3& P) {
    const std::tuple<double, double, double> key{P.x, P.y, P.z};
    auto it = vidx.find(key);
    if (it != vidx.end()) return it->second;
    const int id = static_cast<int>(verts2.size());
    verts2.push_back(proj(P));
    canon3.push_back(P);
    vidx.emplace(key, id);
    return id;
  }
  // Insert at an EXPLICIT 2D position (e.g. an exact on-seam crossing) but
  // dedup by the 3D canonical bits, so every reference to one once-only triple
  // point collapses to a single face vertex even though its projection is off
  // the constructed point by rounding (B1 pre-split).
  int addAt(const vec2& p2, const vec3& canon) {
    const std::tuple<double, double, double> key{canon.x, canon.y, canon.z};
    auto it = vidx.find(key);
    if (it != vidx.end()) return it->second;
    const int id = static_cast<int>(verts2.size());
    verts2.push_back(p2);
    canon3.push_back(canon);
    vidx.emplace(key, id);
    return id;
  }
};
// Build the axis-drop frame for a face with unnormalized outward normal faceN
// and an in-plane point planePt (a triangle vertex).  Returns false (caller
// fails closed) on a degenerate normal.  The dominant-axis pick keeps the
// projection non-degenerate BY CONSTRUCTION (|nHat[axis]| >= 1/sqrt(3)), so an
// edge-on plane is impossible here - no branch, the invariant is asserted.
bool BuildAxisDropFrame(const vec3& faceN, const vec3& planePt,
                        AxisDropFrame& pf) {
  const double nLen = la::length(faceN);
  if (!(nLen > 0.0)) return false;
  pf.nHat = faceN / nLen;
  pf.axis = DominantAxis(faceN);
  pf.flip = pf.nHat[pf.axis] < 0.0;
  pf.planeD = la::dot(pf.nHat, planePt);
  DEBUG_ASSERT(std::abs(pf.nHat[pf.axis]) > 0.0, logicErr,
               "axis-drop: dominant normal component is zero");
  return true;
}

// ---------------------------------------------------------------------------
// B1: ONCE-ONLY 3-FACE TRIPLE POINTS (arrangement-vertex-first, f4-design-b/c).
//
// A transversal 3-face triple point T = plane(f) INT plane(g) INT plane(h) is a
// 0-cell where two seams of face f (from partners g and h) cross in f's
// interior.  The per-face resolver builds it THREE times (once per incident
// face's frame), so the three double-rounded images and three independent
// radial subdivisions disagree and the emission fans open (the 121-open-edge
// wall, oscad-f4).  B1 enumerates each triple ONCE globally, constructs ONE
// canonical double-precision position keyed by the sorted plane triple, and
// threads it into every incident face's overlay by PRE-SPLITTING the seams at
// the shared vertex - so all three faces reference the identical 3D point and
// their sub-faces weld bit-identically (the arrangement is complete at the
// 0-cell).  Tripwire-free: the position is a CONSTRUCTION (like SegPlanePoint),
// never a decision; the crossing test and the winding classify stay the
// existing level-0 predicates.
// ---------------------------------------------------------------------------

// The double solution of the three planes n_i . (x - a_i) = 0, i in {0,1,2},
// via Cramer over the (unnormalized) input normals.  Symmetric in the three
// planes, so a sorted key feeds identical bits for all three incident faces.
// Returns false on a (near-)degenerate triple (planes not independent): the
// caller then declines to pre-split and the crossing falls to the pos2in
// fail-closed backstop - never a garbage constructed vertex.
bool Intersect3Planes(const vec3& n0, const vec3& a0, const vec3& n1,
                      const vec3& a1, const vec3& n2, const vec3& a2,
                      vec3& out) {
  const vec3 c12 = la::cross(n1, n2), c20 = la::cross(n2, n0),
             c01 = la::cross(n0, n1);
  const double det = la::dot(n0, c12);
  if (!(std::abs(det) > 0.0)) return false;
  const double d0 = la::dot(n0, a0), d1 = la::dot(n1, a1), d2 = la::dot(n2, a2);
  out = (d0 * c12 + d1 * c20 + d2 * c01) / det;
  return std::isfinite(out.x) && std::isfinite(out.y) && std::isfinite(out.z);
}

// 2D intersection point of lines through [a,b] and [c,d] (precondition: they
// properly cross, so the denominator is nonzero).  The exact on-seam crossing:
// B1 splits the seam here (strictly interior, no fold-back) and keys the vertex
// by the once-only 3D triple point.
vec2 SegLineIntersect2D(const vec2& a, const vec2& b, const vec2& c,
                        const vec2& d) {
  const vec2 r = b - a, s = d - c;
  const double rxs = la::cross(r, s);
  const double t = la::cross(c - a, s) / rxs;
  return a + t * r;
}

// EXACT in-plane orientation sign of three coplanar 3D points: drop the
// dominant axis of the (unnormalized) face normal and take the exact 2D orient
// of the surviving axis pair as the padded orient3d (embed at z=0, lift the
// first point in +z) through the ONE blessed exact predicate FORM
// (Orient3DExactSign), FILTER-FIRST - the SAME exact helper the input-vertex
// T-junction arm (InputVertexStrictlyOnEdge) now uses for its collinearity
// test.  Dropping the dominant normal axis keeps the projection non-degenerate.
// NO new predicate FORM, NO exact-on-constructed (dyadic input coords only).
// The crossing test below compares only relative signs, so handedness is moot.
// EXPORTED (declared in polygon_internal.h, hence the namespace close/reopen):
// the shared exact triangulation module (manifold::exacttri, polygon.cpp)
// makes every orientation decision through this same predicate form - one
// predicate, one implementation.
}  // namespace

int ExactOrient2DDrop(const vec3& p, const vec3& q, const vec3& r, int axis) {
  auto proj = [&](const vec3& v) -> vec3 {
    if (axis == 0) return {v.y, v.z, 0.0};
    if (axis == 1) return {v.x, v.z, 0.0};
    return {v.x, v.y, 0.0};
  };
  const vec3 pp = proj(p), qp = proj(q), rp = proj(r);
  const vec3 lift{pp.x, pp.y, 1.0};
  const int s = Orient3DFilterSign(pp, qp, rp, lift);
  return s != 0 ? s : Orient3DExactSign(pp, qp, rp, lift);
}

namespace {

// EXACT strict proper crossing of two coplanar 3D segments [p0,p1] and [q0,q1]
// on a face whose normal's dominant axis is `axis`: each segment's endpoints
// strictly straddle the other's supporting line (four ExactOrient2DDrop signs),
// excluding collinear / endpoint-touching contacts.  Replaced the double
// SegProperCross2D (level-0 la::cross on the former rounded orthonormal
// per-face projection) in the triple-point enumeration, which OVER-detected
// near-tangent phantom crossings (dbl=1/ex=0) the exact sign refutes -
// byte-identical on every resolving carrier (double and exact agree there),
// removing only phantom triples on the near-tangent openscad residue.
inline bool ExactSegProperCross(const vec3& p0, const vec3& p1, const vec3& q0,
                                const vec3& q1, int axis) {
  const int o1 = ExactOrient2DDrop(p0, p1, q0, axis);
  const int o2 = ExactOrient2DDrop(p0, p1, q1, axis);
  const int o3 = ExactOrient2DDrop(q0, q1, p0, axis);
  const int o4 = ExactOrient2DDrop(q0, q1, p1, axis);
  return o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && o1 != o2 && o3 != o4;
}

// Is the constructed crossing point X (given as its filter EHPoint eT and, when
// escalated, its exact HPoint via exactT()) strictly INTERIOR to triangle
// (t0,t1,t2) whose normal's dominant axis is `axis`?  Three orient2d of X
// against the directed edges (the general degree-9 instantiation of the ONE
// homogeneous form, filter-first: sos::HomogOrient2DFilter escalating to the
// exact sos::HomogOrient2DExact only on a filter-0).  Strictly interior iff the
// three orientations share ONE nonzero sign (the common sign(prod W) flip
// cancels across the three, so handedness is moot); a zero is ON an edge (not
// strictly interior); all-zero is W==0 (X at infinity - parallel planes).
template <class ExactT>
inline bool HPointStrictlyInTri(const sos::EHPoint& eT, ExactT&& exactT,
                                const vec3& t0, const vec3& t1, const vec3& t2,
                                int axis) {
  auto orient = [&](const vec3& a, const vec3& b) -> int {
    const sos::EHPoint eA = sos::ETrivialHPoint(a), eB = sos::ETrivialHPoint(b);
    const int fs = sos::HomogOrient2DFilter(eA, eB, eT, axis);
    if (fs != 0) return fs;
    const sos::HPoint hA = sos::TrivialHPoint(a), hB = sos::TrivialHPoint(b);
    return sos::HomogOrient2DExact(hA, hB, exactT(), axis);
  };
  const int o0 = orient(t0, t1);
  const int o1 = orient(t1, t2);
  const int o2 = orient(t2, t0);
  if (o0 == 0 || o1 == 0 || o2 == 0) return false;
  return o0 == o1 && o1 == o2;
}

// E1 SYMBOLIC EXTENT (homog-design + e1-plumb): do the two seam SEGMENTS on
// face f - seam_g = tri(f) INT tri(g) and seam_h = tri(f) INT tri(h), sharing
// carrier plane f - cross in f's interior EXACTLY?  Their supporting lines
// (f INT g and f INT h) meet at the constructed triple point X = {f,g,h}
// (Cramer over the three plane coefficients).  Each seam's extent is carried
// SYMBOLICALLY by its bounding triangles: X lies within seam_g iff X in tri(f)
// AND tri(g), within seam_h iff X in tri(f) AND tri(h).  So the two SEGMENTS
// cross iff X is strictly interior to ALL THREE triangles f, g, h - the exact
// segment x segment question, decided by clipping X to each triangle via the
// LANDED degree-9 instantiation (HomogOrient2D) against that triangle's edges.
// This makes the straddle EXACT: unlike ExactSegProperCross on the ROUNDED
// seam-endpoint doubles (which lost near-tangent crossings), and unlike a pure
// line-crossing existence test on the constructed point (X in tri(f) only,
// which IGNORES the segment extent and OVER-detects ~30x - the seam LINES cross
// in f far more often than the finite overlap SEGMENTS do), the tri(g)/tri(h)
// clip IS the exact extent. A boundary landing (an orient 0) is X coincident
// with a seam ENDPOINT = an ALIAS with an existing junction, declined here and
// owned by the junction registry (level-0 incidence, nomerge).  NO new
// predicate FORM: only the two landed instantiations (input-point orient3d
// inside the filter, degree-9 constructed-point orient2d).  rF/rG/rH are the
// canonical coplanar reps for the PLANE bits (coplanar-consistent construction
// of X); f/g/h are the ACTUAL seam triangles whose extents clip it.  `clip`
// selects the extent: 3 (DEFAULT) = X in tri(f) AND tri(g) AND tri(h) = the
// exact segment x segment straddle; 1 = X in tri(f) ALONE = the seam LINES
// cross ignoring the segment extent (the ~30x over-detect - the E1_NOEXTENT
// mutation and the E1_MEASURE census).
// The canonical once-only triple position of a plane triple (rep faces
// r0,r1,r2).  INPUT-EXACT under IX (the Cramer point composed exactly from the
// nine input vertices, rounded to double), else Intersect3Planes on rounded
// faceN.  Used everywhere a triple position is committed so all keys (provOf,
// seamTriples, canon3) share ONE consistent position.  Byte-clean off the wedge
// population (IX-gated, and triples are openscad-only).
inline bool CanonTriplePos(const BuildArrangement& A, int r0, int r1, int r2,
                           vec3& pos) {
  if (IXEnabled()) {
    const sos::BigHPoint bp = IXTripleHPoint(A, r0, r1, r2);
    if (sos::BigSign(bp.W) == 0) return false;  // parallel: degenerate triple
    const vec3 p = sos::BigHPointToPos(bp);
    if (!(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z)))
      return false;
    pos = p;
    return true;
  }
  return Intersect3Planes(A.faceN[r0], A.tri[r0][0], A.faceN[r1], A.tri[r1][0],
                          A.faceN[r2], A.tri[r2][0], pos);
}

bool ExactSeamsCross(const BuildArrangement& A, int f, int g, int h, int rF,
                     int rG, int rH, int clip = 3) {
  const vec3 nF = A.faceN[rF];
  const double dF = la::dot(nF, A.tri[rF][0]);
  const vec3 nG = A.faceN[rG];
  const double dG = la::dot(nG, A.tri[rG][0]);
  const vec3 nH = A.faceN[rH];
  const double dH = la::dot(nH, A.tri[rH][0]);
  const sos::EHPoint eT = sos::ECramerHPoint(nF, dF, nG, dG, nH, dH);
  // The exact Cramer point (Poly) is built at most once per pair and only when
  // a filter escalates (a filter-0), so the certified fast path never touches
  // it.
  bool haveExact = false;
  sos::HPoint hT;
  auto exactT = [&]() -> const sos::HPoint& {
    if (!haveExact) {
      hT = sos::CramerHPoint(nF, dF, nG, dG, nH, dH);
      haveExact = true;
    }
    return hT;
  };
  if (!HPointStrictlyInTri(eT, exactT, A.tri[f][0], A.tri[f][1], A.tri[f][2],
                           DominantAxis(A.faceN[f])))
    return false;
  if (clip < 3) return true;  // tri(f) only (over-detect lever / census)
  return HPointStrictlyInTri(eT, exactT, A.tri[g][0], A.tri[g][1], A.tri[g][2],
                             DominantAxis(A.faceN[g])) &&
         HPointStrictlyInTri(eT, exactT, A.tri[h][0], A.tri[h][1], A.tri[h][2],
                             DominantAxis(A.faceN[h]));
}

// Enumerate the component's 3-face triple points ONCE and record, per seam, the
// canonical 3D positions that split it.  Populates A.seamTriples (parallel to
// A.faceSeams).  A no-op (all-empty) on any component whose seams never cross
// in a face interior (the whole corpus off openscad), so EmitSeamedFace stays
// byte-identical there.
void EnumerateTriplePoints(BuildArrangement& A, double eps) {
  const int nTri = static_cast<int>(A.tri.size());
  A.seamTriples.assign(nTri, {});
  for (int f = 0; f < nTri; ++f)
    A.seamTriples[f].assign(A.faceSeams[f].size(), {});
  gF4BTriplePts.clear();  // measurement only (per component)

  // planeId collapses exactly-coplanar faces onto one plane id, so a triple
  // reached through different coplanar representatives is keyed - and thus
  // constructed - ONCE (openscad is coplanar-dominated).  Coplanar face ids and
  // cluster ids are kept disjoint by the nTri offset.
  const std::vector<int>& planeId = A.planeId;  // hoisted (ResolveComponent)
  // Canonical plane (normal + a point) per plane id = the lowest-index face
  // carrying it, so all three incident faces read identical plane bits.
  std::map<int, int> rep;
  for (int f = 0; f < nTri; ++f) {
    auto it = rep.find(planeId[f]);
    if (it == rep.end()) rep.emplace(planeId[f], f);
  }

  // MUTATION (measurement only): per-face back-projection of the crossing
  // instead of the once-only canonical point.  The three incident faces then
  // get three double-rounded images -> the fans reopen (proves the once-only
  // keying is load-bearing; f4-design-b P4, f4-design-a baseline).
  const bool perFace = std::getenv("F4B_PERFACE") != nullptr;

  // E1 SYMBOLIC EXTENT (homog-design + e1-plumb): the seam-seam crossing
  // EXISTENCE test is the EXACT segment x segment straddle - the constructed
  // triple point X = {f,g,h} strictly interior to ALL THREE seam triangles f,
  // g, h (ExactSeamsCross), each seam's extent carried symbolically by its
  // bounding triangles.  This is decided on the symbolic construction, not the
  // eps- truncated seam-endpoint doubles, and the tri(g)/tri(h) extent clip
  // stops the ~30x line-crossing over-detection (the seam LINES cross in f far
  // more often than the finite overlap SEGMENTS do).
  //
  // DECISIVE MEASUREMENT (e1-plumb): the PURE exact straddle (E1_PURE, X in all
  // three triangles with NO rounded-straddle gate) collapses the over-detect to
  // the true crossing set, but wiring it as the sole existence test is NOT
  // byte-clean and BREAKS a resolving carrier: on GT7081's 0.002deg
  // near-parallel twins it admits genuine near-tangent crossings ~1.4e-7 from a
  // near-coincident corner (endpoints diverged ABOVE eps ~4.5e-8 - the
  // near-parallel-plane wall), and on openscad it splits the near-tangent thin
  // cells the rounded straddle left unsplit; the DOUBLE-PRECISION downstream
  // (RemoveOverlaps2D seam sub-face / SplitTouchingSheets) cannot represent
  // those sub-eps crossings, so GT7081 regresses to a 4-open-edge fail-closed
  // and openscad worsens (21 -> 35 opens / b3 in 12 faces).  This DECISIVELY
  // confirms the terminal is the exact-rational 2D arrangement (it blocks even
  // a RESOLVING carrier, not just openscad), not the enumeration predicate.
  //
  // DEFAULT (representability-safe): keep the rounded straddle
  // (ExactSegProperCross) as a REPRESENTABILITY GATE - split only where the
  // finite rounded segments actually cross, i.e. where the double-precision
  // arrangement can place the crossing - and REFINE that set with the exact
  // extent (drop the rounded false positives whose exact X falls outside seam
  // g's or h's triangle).  Byte-clean on the whole resolving corpus (no rounded
  // seam pair crosses there - zero triple points, so ExactSeamsCross is never
  // reached), and a strict refinement of the pre-plumb X-in-f default (can only
  // drop, never add).  Mutation levers: E1_OFF reverts to the pre-plumb
  // rounded-straddle-and-X-in-f default; E1_PURE is the ungated exact straddle
  // (the decisive negative - breaks GT7081); E1_MEASURE censuses the crossing
  // set.
  static const bool kE1Measure = std::getenv("E1_MEASURE") != nullptr;
  static const bool kE1Pure = std::getenv("E1_PURE") != nullptr;
  static const bool kE1Off = std::getenv("E1_OFF") != nullptr;
  int e1New = 0, e1Old = 0, e1NewOnly = 0, e1OldOnly = 0, e1InFOnly = 0;
  // EX2 CENSUS (exact-2d overlay lane, measurement only): per seamed face, the
  // gap-free wedge crossings = the exact segment straddle (ExactSeamsCross
  // clip=3) that the DEFAULT rounded-straddle gate DECLINES, at dEnd > 14 eps
  // from every seam endpoint (gap-free, not a gap-bounded near-corner / alias).
  static const bool kEx2Census = std::getenv("EX2_CENSUS") != nullptr;
  std::map<int, std::vector<double>> ex2GapFree;  // face -> dEnd(eps) list

  std::map<std::array<int, 3>, vec3> tripleTab;  // sorted plane triple -> pos
  for (int f = 0; f < nTri; ++f) {
    if (!A.seamed[f]) continue;
    const int ns = static_cast<int>(A.faceSeams[f].size());
    if (ns < 2) continue;
    AxisDropFrame pf;
    if (!BuildAxisDropFrame(A.faceN[f], A.tri[f][0], pf)) continue;
    std::vector<std::array<vec2, 2>> seg(ns);
    for (int k = 0; k < ns; ++k) {
      seg[k][0] = pf.proj(A.faceSeams[f][k].p0);
      seg[k][1] = pf.proj(A.faceSeams[f][k].p1);
    }
    // Exact in-plane orientation drops the dominant normal axis (the frame's
    // dropped axis, computed once per face).
    const int axis = pf.axis;
    for (int k1 = 0; k1 < ns; ++k1) {
      const int pg = planeId[A.faceSeams[f][k1].other];
      if (pg == planeId[f]) continue;
      for (int k2 = k1 + 1; k2 < ns; ++k2) {
        const int ph = planeId[A.faceSeams[f][k2].other];
        if (ph == planeId[f] || ph == pg) continue;  // collinear / degenerate
        const int gFace = A.faceSeams[f][k1].other;
        const int hFace = A.faceSeams[f][k2].other;
        const int rF = rep[planeId[f]], rG = rep[pg], rH = rep[ph];
        // The rounded-endpoint straddle = the REPRESENTABILITY GATE (the finite
        // rounded segments actually cross, so the double-precision arrangement
        // can place the crossing).  Exact-on-rounded via ExactOrient2DDrop,
        // cheap.
        const bool straddle =
            kE1Pure ? true
                    : ExactSegProperCross(
                          A.faceSeams[f][k1].p0, A.faceSeams[f][k1].p1,
                          A.faceSeams[f][k2].p0, A.faceSeams[f][k2].p1, axis);
        if (kE1Measure) {
          // Census the exact-extent (segment) set, the pre-plumb X-in-f set,
          // and the over-detect line-crossing set (X-in-f alone).  Lazy,
          // measure-only.
          const bool inFGH = ExactSeamsCross(A, f, gFace, hFace, rF, rG, rH);
          const bool inF =
              ExactSeamsCross(A, f, gFace, hFace, rF, rG, rH, /*clip=*/1);
          const bool str = ExactSegProperCross(
              A.faceSeams[f][k1].p0, A.faceSeams[f][k1].p1,
              A.faceSeams[f][k2].p0, A.faceSeams[f][k2].p1, axis);
          if (str && inFGH) ++e1New;  // the shipped default set
          if (str && inF) ++e1Old;    // the pre-plumb default set
          if (inFGH) ++e1NewOnly;     // the pure exact-segment set (ungated)
          if (str && inF && !inFGH)
            ++e1OldOnly;                   // rounded phantoms the clip drops
          if (inF && !inFGH) ++e1InFOnly;  // extent-clip removals
        }
        if (kEx2Census && !straddle &&
            ExactSeamsCross(A, f, gFace, hFace, rF, rG, rH)) {
          // A DECLINED exact segment crossing (rounded straddle false): the
          // candidate gap-free wedge population.  Measure dEnd to the 4 seam
          // endpoints on the once-only 3D triple point.
          std::array<int, 3> ck = {planeId[f], pg, ph};
          std::sort(ck.begin(), ck.end());
          vec3 cpos;
          if (Intersect3Planes(A.faceN[rep[ck[0]]], A.tri[rep[ck[0]]][0],
                               A.faceN[rep[ck[1]]], A.tri[rep[ck[1]]][0],
                               A.faceN[rep[ck[2]]], A.tri[rep[ck[2]]][0],
                               cpos)) {
            const double dEnd =
                std::min(std::min(la::length(cpos - A.faceSeams[f][k1].p0),
                                  la::length(cpos - A.faceSeams[f][k1].p1)),
                         std::min(la::length(cpos - A.faceSeams[f][k2].p0),
                                  la::length(cpos - A.faceSeams[f][k2].p1)));
            ex2GapFree[f].push_back(dEnd / eps);
          }
        }
        // DEFAULT: representability-gated exact extent (rounded straddle AND
        // the exact segment straddle).  E1_PURE drops the gate (straddle==true)
        // for the decisive negative; the clip (X in tri(g)/tri(h)) is the exact
        // extent.
        const int clip = kE1Off ? 1 : 3;
        if (!(straddle &&
              ExactSeamsCross(A, f, gFace, hFace, rF, rG, rH, clip)))
          continue;
        // The on-seam split position: the 2D crossing of the two ROUNDED seam
        // segments (SegLineIntersect2D), which lies ON both 2D seam lines so
        // the pre-split chain has no kink (the split coordinate stays collinear
        // with each seam's endpoints - projecting the exact 3D point instead
        // kinks the rounded chain and RemoveOverlaps2D then re-crosses it).  It
        // only orders the pre-split; the emitted vertex is keyed by the
        // once-only 3D triple point pos, not this 2D position.  The
        // exact-extent test guarantees the crossing is strictly interior to
        // both segments, so x is interior (no extrapolation).
        const vec2 x =
            SegLineIntersect2D(seg[k1][0], seg[k1][1], seg[k2][0], seg[k2][1]);
        if (!(std::isfinite(x.x) && std::isfinite(x.y))) continue;
        vec3 pos;
        bool cached = false;
        std::array<int, 3> key = {planeId[f], pg, ph};
        if (perFace) {
          pos = pf.lift(x);  // this face's own image of the crossing
          if (!(std::isfinite(pos.x) && std::isfinite(pos.y) &&
                std::isfinite(pos.z)))
            continue;
        } else {
          std::sort(key.begin(), key.end());
          auto it = tripleTab.find(key);
          if (it != tripleTab.end()) {
            pos = it->second;
            cached = true;
          } else {
            const int r0 = rep[key[0]], r1 = rep[key[1]], r2 = rep[key[2]];
            if (!CanonTriplePos(A, r0, r1, r2, pos))
              continue;  // degenerate triple: leave to the pos2in backstop
          }
        }
        // ALIASING / SUB-EPS COLLAPSE (nomerge witness theorem): a proper
        // crossing must be strictly interior to both SEGMENTS, whose extents
        // end at the ROUNDED seam endpoints.  If X lands within the weld radius
        // of a seam ENDPOINT (of k1 or k2 - both on the same seam line as X,
        // sharing the {f,g}/{f,h} carrier), the two seams MEET at that shared
        // junction within eps rather than crossing in the interior; registering
        // a split there manufactures a sub-eps-degenerate seam the emission
        // weld would collapse anyway (the same within-construction eps collapse
        // RecordSeams does on its endpoints).  Decline: the endpoint is already
        // carried as a junction. Under the shipped straddle-gated default this
        // is a near-no-op (a rounded proper crossing is interior); it is
        // load-bearing under E1_PURE, where it collapses the ungated
        // near-endpoint crossings (e.g. SelfIntersectB's 4.5e-8-long seam whose
        // exact X lands 2e-15 from an endpoint).  Checked before the tripleTab
        // emplace so a declined pair registers no triple.
        if (!kE1Off && (la::length(pos - A.faceSeams[f][k1].p0) <= eps ||
                        la::length(pos - A.faceSeams[f][k1].p1) <= eps ||
                        la::length(pos - A.faceSeams[f][k2].p0) <= eps ||
                        la::length(pos - A.faceSeams[f][k2].p1) <= eps))
          continue;
        if (!perFace && !cached) tripleTab.emplace(key, pos);
        A.seamTriples[f][k1].push_back({x, pos});
        A.seamTriples[f][k2].push_back({x, pos});
        if (std::getenv("F4B_DUMP") != nullptr)
          gF4BTriplePts.insert({pos.x, pos.y, pos.z});
      }
    }
  }
  if (std::getenv("F4B_DUMP") != nullptr) {  // measurement only
    int inc = 0;
    for (const auto& ff : A.seamTriples)
      for (const auto& sk : ff) inc += static_cast<int>(sk.size());
    std::fprintf(stderr, "F4B_TRIPLES distinct=%d incidences=%d perFace=%d\n",
                 static_cast<int>(tripleTab.size()), inc, perFace ? 1 : 0);
  }
  if (kE1Measure)
    std::fprintf(
        stderr,
        "E1_CENSUS default=%d prePlumb=%d pureExact=%d phantomsDropped=%d "
        "extentClipRemovals=%d\n",
        e1New, e1Old, e1NewOnly, e1OldOnly, e1InFOnly);
  if (kEx2Census) {
    int nGapFreeFaces = 0, nGapFreeCross = 0, nBoundedDeclined = 0;
    for (const auto& fv : ex2GapFree) {
      int gf = 0;
      for (double d : fv.second)
        if (d > 14.0)
          ++gf;
        else
          ++nBoundedDeclined;
      if (gf > 0) {
        ++nGapFreeFaces;
        nGapFreeCross += gf;
        std::fprintf(stderr,
                     "EX2_FACE f=%d gapFreeCross=%d dEnds(eps)=", fv.first, gf);
        for (double d : fv.second) std::fprintf(stderr, "%.1f ", d);
        std::fprintf(stderr, "\n");
      }
    }
    std::fprintf(
        stderr,
        "EX2_CENSUS gapFreeFaces=%d gapFreeCross=%d boundedDeclined=%d\n",
        nGapFreeFaces, nGapFreeCross, nBoundedDeclined);
  }
}

// ---------------------------------------------------------------------------
// EX2 WEDGE SPLITS (exact2d overlay lane): the coordinated per-triple pass that
// registers the GAP-FREE wedge crossings the DEFAULT rounded-straddle gate
// declines, and FLAGS the incident faces for the exact PSLG overlay.
//
// EnumerateTriplePoints (default) registers only crossings whose ROUNDED seam
// segments straddle - the crossings the double-precision RemoveOverlaps2D can
// place.  openscad's diverging near-parallel wedges cross EXACTLY
// (ExactSeamsCross clip=3) at an interior point the rounded segments do NOT
// straddle (the wedge is sub-eps thin there), so the default declines them and
// RemoveOverlaps2D's eps-merge / sweep mangles the wide part (the 21-open-edge
// residue).  This pass finds every exact interior crossing (ungated
// ExactSeamsCross), keys it by the sorted plane triple, and CLASSIFIES it ONCE
// (coordinated, not per-face) by the MIN distance from its once-only 3D triple
// point (Intersect3Planes, the same bits the default constructs) to any
// incident seam endpoint:
//   dEnd  >  kWedgeEnv eps  => GAP-FREE (a diverging wedge; pivot's scale-free
//     spectrum): register the split on ALL incident faces + flag them.
//   dEnd  <= kWedgeEnv eps  => GAP-BOUNDED (a near-coincident corner / alias -
//     GT7081's near-corner E1_PURE crossings, the 330 cross-provenance welds):
//     DECLINE - collapse to the existing endpoint junction (the coordinated
//     collapse; no split, no flag).  This is what keeps the RESOLVING carriers
//     byte-clean: GT7081's crossings are all <= kWedgeEnv eps -> zero flagged
//     faces (canary); the whole corpus off openscad has zero gap-free
//     crossings.
// The envelope sits in the measured spectral GAP (pivot/hybrid): GT7081's
// bounded crossings <= 14 eps, openscad's wedges >= 26.8 eps.
// ---------------------------------------------------------------------------
constexpr double kWedgeEnv = 14.0;  // coordinated-collapse envelope, eps units
void EnumerateWedgeSplits(BuildArrangement& A, double eps) {
  const int nTri = static_cast<int>(A.tri.size());
  A.wedgeFace.assign(nTri, 0);
  if (std::getenv("EX2_OFF") != nullptr) return;  // mutation lever: overlay off
  const std::vector<int>& planeId = A.planeId;    // hoisted (ResolveComponent)
  std::map<int, int> rep;
  for (int f = 0; f < nTri; ++f)
    rep.emplace(planeId[f], f);  // lowest-index rep
  struct Cand {
    vec3 pos;
    double minDEnd = std::numeric_limits<double>::max();
    bool havePos = false;
    std::vector<std::array<int, 3>> inc;  // (face, seam k1, seam k2)
    std::vector<vec2> x2d;        // matching on-seam 2D crossing per inc
    std::vector<double> dEndInc;  // per-incidence dEnd (this face's ends)
  };
  std::map<std::array<int, 3>, Cand> tab;
  for (int f = 0; f < nTri; ++f) {
    if (!A.seamed[f]) continue;
    const int ns = static_cast<int>(A.faceSeams[f].size());
    if (ns < 2) continue;
    AxisDropFrame pf;
    if (!BuildAxisDropFrame(A.faceN[f], A.tri[f][0], pf)) continue;
    for (int k1 = 0; k1 < ns; ++k1) {
      const int pg = planeId[A.faceSeams[f][k1].other];
      if (pg == planeId[f]) continue;
      for (int k2 = k1 + 1; k2 < ns; ++k2) {
        const int ph = planeId[A.faceSeams[f][k2].other];
        if (ph == planeId[f] || ph == pg) continue;
        const int gFace = A.faceSeams[f][k1].other;
        const int hFace = A.faceSeams[f][k2].other;
        const int rF = rep[planeId[f]], rG = rep[pg], rH = rep[ph];
        if (!ExactSeamsCross(A, f, gFace, hFace, rF, rG, rH)) continue;
        std::array<int, 3> key = {planeId[f], pg, ph};
        std::sort(key.begin(), key.end());
        Cand& c = tab[key];
        if (!c.havePos) {
          const int r0 = rep[key[0]], r1 = rep[key[1]], r2 = rep[key[2]];
          if (!CanonTriplePos(A, r0, r1, r2, c.pos))
            continue;  // degenerate triple: leave to the default backstop
          c.havePos = true;
        }
        const double dEnd =
            std::min(std::min(la::length(c.pos - A.faceSeams[f][k1].p0),
                              la::length(c.pos - A.faceSeams[f][k1].p1)),
                     std::min(la::length(c.pos - A.faceSeams[f][k2].p0),
                              la::length(c.pos - A.faceSeams[f][k2].p1)));
        c.minDEnd = std::min(c.minDEnd, dEnd);
        // The wedge crossing's 2D position for the pre-split = the rounded
        // seam- line intersection.  NOTE (the sharpened terminal): the DECLINED
        // gap-free crossings have the rounded segments NOT straddling, so NO 2D
        // position is interior to BOTH rounded segments - this line
        // intersection extrapolates beyond a rounded endpoint
        // (SegLineIntersect2D) and pf.proj(pos) lands just off-segment too.
        // Either way the exact overlay's ExtractCells (which walks the ROUNDED
        // verts2) folds back, even though the exact PSLG completeness check (on
        // the exact canon3 triple point) passes.  The crossing is genuine and
        // the arrangement is combinatorially complete; the double-precision
        // cell walk cannot realize it.  Only symbolic (exact) seam endpoints in
        // the cell walk close it (the priced completion).
        const vec2 x = SegLineIntersect2D(
            pf.proj(A.faceSeams[f][k1].p0), pf.proj(A.faceSeams[f][k1].p1),
            pf.proj(A.faceSeams[f][k2].p0), pf.proj(A.faceSeams[f][k2].p1));
        if (!(std::isfinite(x.x) && std::isfinite(x.y))) continue;
        c.inc.push_back({f, k1, k2});
        c.x2d.push_back(x);
        c.dEndInc.push_back(dEnd);
      }
    }
  }
  int nFlagFaces = 0, nSplitsAdded = 0;
  std::set<int> flagged;
  for (auto& kv : tab) {
    Cand& c = kv.second;
    if (!c.havePos || c.inc.empty()) continue;
    if (c.minDEnd <= kWedgeEnv * eps)
      continue;  // gap-bounded: coordinated collapse
    for (size_t i = 0; i < c.inc.size(); ++i) {
      const int f = c.inc[i][0], k1 = c.inc[i][1], k2 = c.inc[i][2];
      auto addSplit = [&](int k) -> bool {
        for (const auto& xp : A.seamTriples[f][k])
          if (la::length(xp.second - c.pos) <= eps) return false;  // already
        A.seamTriples[f][k].push_back({c.x2d[i], c.pos});
        ++nSplitsAdded;
        return true;
      };
      // A face routes through the exact overlay ONLY if it gets a genuinely NEW
      // split - i.e. it carries a DECLINED gap-free crossing RemoveOverlaps2D
      // cannot place.  An incident face whose view of the same triple was
      // already default-registered (rounded straddle true) keeps
      // RemoveOverlaps2D (which handles that wide crossing) and still splits at
      // the SAME once-only 3D point, so the shared seam welds across the mixed
      // routing.
      const bool a1 = addSplit(k1);
      const bool a2 = addSplit(k2);
      if (a1 || a2) {
        A.wedgeFace[f] = 1;
        flagged.insert(f);
      }
    }
  }
  // CLUSTER COVERAGE (coverage lane): flag every face incident to a NEAR-
  // COINCIDENT triple cluster (>=2 genuinely-distinct triple points, DIFFERENT
  // plane triples, within R eps).  Under input-exact positions the crowded
  // openscad cluster's crossings are all gap-BOUNDED (dEnd <= 14 eps) so the
  // gap-free path above declines them - yet they ARE the 21-open residue.  This
  // routes the whole crowded cluster through the exact overlay (the provenance-
  // derived coverage the input-exact basis enables), not the 8 hand-named gap-
  // free faces; the completeness pass below then registers each flagged face's
  // interior crossings.  Measurement lever EX2_CLUSTER=<R in eps> (default
  // OFF).
  if (const char* cs = std::getenv("EX2_CLUSTER")) {
    const double R = (std::atof(cs) > 0.0 ? std::atof(cs) : 1000.0) * eps;
    std::vector<std::array<int, 3>> keys;
    std::vector<vec3> cpos;
    for (auto& kv : tab)
      if (kv.second.havePos && !kv.second.inc.empty()) {
        keys.push_back(kv.first);
        cpos.push_back(kv.second.pos);
      }
    std::vector<char> crowded(keys.size(), 0);
    for (size_t i = 0; i < keys.size(); ++i)
      for (size_t j = i + 1; j < keys.size(); ++j)
        if (keys[i] != keys[j] && la::length(cpos[i] - cpos[j]) <= R)
          crowded[i] = crowded[j] = 1;
    size_t ci = 0;
    for (auto& kv : tab) {
      if (!(kv.second.havePos && !kv.second.inc.empty())) continue;
      if (crowded[ci])
        for (const auto& in3 : kv.second.inc) {
          A.wedgeFace[in3[0]] = 1;
          flagged.insert(in3[0]);
        }
      ++ci;
    }
  }
  // COMPLETENESS ON FLAGGED FACES (symwalk): a flagged wedge face's exact
  // overlay must include EVERY genuine exact crossing INTERIOR TO IT, not only
  // the gap-free one that flagged it.  The wedge crossings are
  // gap-free/scale-free (pivot HALF-1): a crossing can be a sub-eps ALIAS on
  // one incident face (GLOBAL minDEnd <= eps, collapsed) yet a genuine
  // near-corner crossing (14.77 eps) on face f.  Classify PER-FACE: register
  // the split on flagged face f whenever it is > eps from f's OWN seam
  // endpoints (genuine here), so f's arrangement is exact-complete; the alias
  // faces collapse it via the eps snap-round.  An unflagged face (the whole
  // resolving corpus; GT7081) is never touched, so the E1_PURE canary holds.
  for (auto& kv : tab) {
    Cand& c = kv.second;
    if (!c.havePos) continue;
    for (size_t i = 0; i < c.inc.size(); ++i) {
      const int f = c.inc[i][0], k1 = c.inc[i][1], k2 = c.inc[i][2];
      if (!A.wedgeFace[f]) continue;      // flagged faces only (canary safety)
      if (c.dEndInc[i] <= eps) continue;  // alias on this face: snap-collapse
      auto addSplit = [&](int k) {
        for (const auto& xp : A.seamTriples[f][k])
          if (la::length(xp.second - c.pos) <= eps) return;  // already there
        A.seamTriples[f][k].push_back({c.x2d[i], c.pos});
        ++nSplitsAdded;
      };
      addSplit(k1);
      addSplit(k2);
    }
  }
  if (std::getenv("EX2_ALL") != nullptr)  // diagnostic: overlay on ALL seamed
    for (int f = 0; f < nTri; ++f)
      if (A.seamed[f]) A.wedgeFace[f] = 1;
  nFlagFaces = static_cast<int>(flagged.size());
  if (std::getenv("EX2_CENSUS") != nullptr) {
    std::fprintf(stderr,
                 "EX2_WEDGE flaggedFaces=%d splitsAdded=%d triples=%d\n",
                 nFlagFaces, nSplitsAdded, static_cast<int>(tab.size()));
    for (int fl : flagged) {
      const vec3 c = (A.tri[fl][0] + A.tri[fl][1] + A.tri[fl][2]) * (1.0 / 3.0);
      std::fprintf(stderr, "EX2_FLAGPOS f=%d centroid=(%.5g,%.5g,%.5g)\n", fl,
                   c.x, c.y, c.z);
    }
  }
}

// ---------------------------------------------------------------------------
// GLOBAL JUNCTION REGISTRY (f4-junction): the once-only completion of the
// per-face arrangement at NON-proper-crossing junctions.  EnumerateTriplePoints
// only welds PROPER seam X-crossings; the dominant openscad open residue is
// T-junctions where a seam ENDPOINT (or a triple that TERMINATES a seam on the
// third face) lands strictly interior to a neighbour / partner / third face's
// emitted edge, which that face fails to split -> an unbalanced fan.  The
// registry gathers EVERY once-only arrangement vertex (all seam endpoints + all
// triple points), deduped to a canonical 3D position, so all three emit paths
// (seamed, clean, fold) can split each emitted edge at every registry vertex
// strictly interior to it - the SAME point on every incident face.
// ---------------------------------------------------------------------------

// A1 INPUT-VERTEX-ON-EDGE T-JUNCTIONS (f4-r5): the seam/triple registry omits
// plain input vertices, but a self-overlap can land an INPUT vertex STRICTLY
// INTERIOR to a FOREIGN triangle's edge.  That foreign face emits the edge
// unsplit while the vertex's own emission fan terminates there -> an unbalanced
// (open) fan.  Registering the vertex makes the foreign face split its edge at
// the SAME once-only point, closing the fan.  The decision is EXACT and level-0
// (input mesh doubles only, never a constructed point): (a) collinear - the 3D
// cross (V-a)x(b-a)==0, tested as the three axis-projected 2D orients, each the
// padded orient3d through the ONE blessed exact predicate FORM, filter-first;
// and (b) strictly between the endpoints - an exact strict compare on the
// widest axis, well-defined once collinear.  The eps registry extension (all
// verts, eps on-segment) splits near-touching corners and is MEASURED
// rail-breaking (f4-r4); the exact predicate registers ONLY genuine
// T-junctions, so the arm is a no-op wherever no input vertex is EXACTLY on a
// foreign edge (the whole corpus off openscad).

// V exactly collinear with [a,b] AND strictly interior (V != a, V != b).
bool InputVertexStrictlyOnEdge(const vec3& V, const vec3& a, const vec3& b) {
  // Collinear iff all three axis-drop 2D orients vanish - the same exact FORM
  // as the seam-crossing test's ExactOrient2DDrop (the drop-y case is coord-
  // order transposed, a sign flip that is invariant under the == 0 test).
  if (ExactOrient2DDrop(V, a, b, 2) != 0) return false;  // drop z (x,y)
  if (ExactOrient2DDrop(V, a, b, 0) != 0) return false;  // drop x (y,z)
  if (ExactOrient2DDrop(V, a, b, 1) != 0) return false;  // drop y (x,z)
  // Collinear: strict between-ness on the widest axis (a != b -> nonzero span),
  // an exact strict double compare (endpoints excluded).
  const vec3 d = b - a;
  const double adx = std::abs(d.x), ady = std::abs(d.y), adz = std::abs(d.z);
  double av, bv, vv;
  if (adx >= ady && adx >= adz) {
    av = a.x;
    bv = b.x;
    vv = V.x;
  } else if (ady >= adz) {
    av = a.y;
    bv = b.y;
    vv = V.y;
  } else {
    av = a.z;
    bv = b.z;
    vv = V.z;
  }
  return (av < vv && vv < bv) || (bv < vv && vv < av);
}

// The genuine input-vertex-on-edge T-junction positions of the component: for
// each distinct input vertex, a triangle-AABB broadphase (the winding BVH)
// gates the exact on-edge confirm; the vertex's OWN incident triangles are
// skipped (adjacency), so only foreign self-overlap T-junctions register.
std::vector<vec3> CollectInputVertexTJunctions(const BuildArrangement& A,
                                               double eps) {
  std::vector<vec3> out;
  const int nTri = static_cast<int>(A.tri.size());
  if (nTri == 0) return out;
  Box bBox;
  for (int t = 0; t < nTri; ++t)
    for (int k = 0; k < 3; ++k) bBox.Union(A.tri[t][k]);
  const TriWindBVH bvh = BuildTriWindBVH(A.tri, bBox);
  std::unordered_map<int, vec3> vpos;  // global vid -> position
  for (int t = 0; t < nTri; ++t)
    for (int k = 0; k < 3; ++k) vpos[A.vid[t][k]] = A.tri[t][k];
  const vec3 pad(eps, eps, eps);
  std::vector<int> cand;
  for (const auto& kv : vpos) {
    const int vid = kv.first;
    const vec3& V = kv.second;
    cand.clear();
    auto rec = [&](int, int leaf) { cand.push_back(bvh.leaf2tri[leaf]); };
    auto recorder = MakeSimpleRecorder(rec);
    const Box qbox(V - pad, V + pad);
    auto qf = [&](int) { return qbox; };
    bvh.collider.Collisions<false>(recorder, qf, 1, /*parallel=*/false,
                                   nullptr);
    bool hit = false;
    for (int t : cand) {
      if (A.vid[t][0] == vid || A.vid[t][1] == vid || A.vid[t][2] == vid)
        continue;  // adjacency: skip the vertex's own incident triangles
      for (int e = 0; e < 3; ++e) {
        const vec3& a = A.tri[t][e];
        const vec3& b = A.tri[t][(e + 1) % 3];
        // Broadphase gate: V within eps of the segment interior (cheap double
        // test) before the exact on-edge confirm.
        const vec3 d3 = b - a;
        const double len2 = la::dot(d3, d3);
        if (!(len2 > 0.0)) continue;
        const vec3 w = V - a;
        const double tp = la::dot(w, d3) / len2;
        if (!(tp > 0.0 && tp < 1.0)) continue;
        if (la::length(w - tp * d3) > eps) continue;
        if (InputVertexStrictlyOnEdge(V, a, b)) {
          hit = true;
          break;
        }
      }
      if (hit) break;
    }
    if (hit) out.push_back(V);
  }
  return out;
}

// Build A.junctions: all seam endpoints + all triple points + the exact
// input-vertex-on-edge T-junctions, deduped so no two registry vertices are
// within eps (a within-eps merge to the sorted-lowest representative -
// deterministic, order-independent; the seam endpoints and triples are already
// once-only constructions, so equal geometry gives equal or within-eps bits).
// Empty of interior landings on any complete arrangement (the whole corpus off
// openscad), so consumption stays byte-identical there.
void BuildJunctionRegistry(BuildArrangement& A, double eps) {
  // MUTATION lever (measurement only): F4J_NOREG leaves the registry empty, so
  // every emit path reverts to its pre-junction single-edge push (proves the
  // registry is load-bearing - the openscad open residue jumps back up).
  if (std::getenv("F4J_NOREG") != nullptr) return;
  const int nTri = static_cast<int>(A.tri.size());
  std::vector<vec3> raw;
  for (int f = 0; f < nTri; ++f)
    for (const auto& s : A.faceSeams[f]) {
      raw.push_back(s.p0);
      raw.push_back(s.p1);
    }
  for (const auto& ff : A.seamTriples)
    for (const auto& sk : ff)
      for (const auto& xp : sk) raw.push_back(xp.second);
  // A1 arm (f4-r5): the exact input-vertex-on-edge T-junctions.  F4R_NOVJUNC
  // leaves them out (mutation lever: the openscad A1 fans reopen).
  if (std::getenv("F4R_NOVJUNC") == nullptr)
    for (const vec3& v : CollectInputVertexTJunctions(A, eps)) raw.push_back(v);
  std::sort(raw.begin(), raw.end(), [](const vec3& p, const vec3& q) {
    return std::tie(p.x, p.y, p.z) < std::tie(q.x, q.y, q.z);
  });
  A.junctions.clear();
  for (const vec3& p : raw) {
    bool dup = false;
    for (int k = static_cast<int>(A.junctions.size()) - 1; k >= 0; --k) {
      // sorted by x first: once an accepted vertex is more than eps below in x,
      // no earlier one can be within eps.
      if (A.junctions[k].x < p.x - eps) break;
      if (la::length(A.junctions[k] - p) <= eps) {
        dup = true;
        break;
      }
    }
    if (!dup) A.junctions.push_back(p);
  }
}

// Winding seeds: a few far points in unrelated directions off the bbox (the
// coupled winding is single-valued off-surface, so any certified seed is
// authoritative; several give the robust-winding graze fallbacks).
std::vector<vec3> WindingSeeds(const Box& bBox) {
  const vec3 c = bBox.Center();
  const double L = bBox.Scale() + 1.0;
  return {c + L * vec3(3.13, 5.71, 1.37),   c + L * vec3(-2.71, 1.41, 4.19),
          c + L * vec3(1.73, -3.31, -2.23), c + L * vec3(-4.27, -1.19, 2.83),
          c + L * vec3(2.39, -4.61, 3.07),  c + L * vec3(-1.51, 3.89, -4.43)};
}

// ========================= E1 COORDINATED ENGINE ===========================
// (e1engine notebook.)  Fallback resolver for the component class the per-face
// emission FAILS CLOSED on (openscad's near-tangent multiply-wound cluster).
// It is reached ONLY after EmitComponentBoundary returns a fatal, so it is
// byte-clean on the whole resolving corpus BY CONSTRUCTION (never executed),
// and it cannot weaken fail-closed (on its own failure the ORIGINAL fatal is
// returned unchanged; its success passes the same BuildImpl gate + component
// re-gate as any resolve).
//
// The recipe (validated offline in exact rationals, e1engine notebook: the
// coordinated emission closes the openscad component at 0 open edges,
// GWN-checked against the exact winding oracle):
//   1. PER-PLANE (not per-face) exact 2D arrangements: faces grouped by exact
//      GEOMETRIC coplanarity (anti-oriented coplanar faces share a group with
//      sign -1), so coincident sheets are classified ONCE with a net covering
//      jump - the per-face two-sided-consistency gap dissolves structurally.
//   2. Segments per group: member triangle edges + the recorded seams
//      (A.faceSeams).  Cross-plane subdivision T-consistency is BY SHARED
//      DOUBLES: seam endpoints are stored once per pair (RecordSeams), triple
//      crossings once per sorted plane triple (CanonTriplePos, input-exact),
//      junction splits from the shared registry (A.junctions) - every group
//      derives the identical split point bits, so coincident sub-edges weld.
//   3. Crossing EXISTENCE is INPUT-EXACT and symmetric (IXOrient2D strict
//      point-in-triangle of the triple against all four bounding triangles),
//      so any two groups agree on every split by construction.  No filter in
//      front (a rounded-basis certificate can disagree with input-exact truth
//      - inexact-basis); the engine runs once, on one failing component.
//   4. Cells by rotation-system trace with an EXACT angular comparator
//      (ExactOrient2DDrop on the shared doubles - no atan2, the exact2d
//      fold-back lesson).  Pinched loops are kept WHOLE (a cell's islands ring
//      it as pinch-connected holes; splitting them off and dropping the
//      negative lobes drops the island boundary = the offline v3 bug);
//      disconnected negative loops are island hole-rings, attached to their
//      containing cell by parity + keyhole bridge.
//   5. Per-cell classify: covering jump (member orientation signs at an
//      interior point) + winding probes at an ADAPTIVE normal offset (below
//      half the distance to every nearby foreign plane - the sound
//      "infinitesimal"; the fixed eps*nHat probe overshoots near-tangent
//      sheets).  COMPLETENESS CERTIFICATE: the probed winding delta must equal
//      the combinatorial jump, else fail closed (zero-oracle-wrong).
//   6. Boundary cells ({w>=1} transition) triangulated as polygons-WITH-HOLES
//      by the shared EXACT triangulator (manifold::exacttri, polygon.cpp;
//      exact diagonal-split - earclip is not hole-safe), oriented per-triangle
//      toward the exterior; assembled by the ordinary BuildImpl eps-weld (the
//      offline run proved the once-rounded coordinated soup position-welds
//      CLOSED - identity plumbing through the weld is not needed: shared
//      doubles make every coincident vertex byte-equal).
namespace e1 {

// The exact drop-frame helpers + the exact earclip/diagonal-split triangulator
// RELOCATED to the shared triangulation module (manifold::exacttri,
// polygon.cpp + polygon_internal.h) - the s2-tri unification.  The resolver's
// e1:: spelling stays valid at every call site via these using-declarations;
// the module comment in polygon.cpp records why the engine's cell
// triangulation must be the EXACT mode (dust-cell orientation is noise at the
// weld scale; only the exact deterministic diagonalization keeps shared-edge
// emission anti-correlated across adjacent cells).
using exacttri::Drop2;
using exacttri::LoopShoelace;
using exacttri::O2;
using exacttri::OnOpenSeg2;
using exacttri::ProperCross2;

using K3 = std::tuple<double, double, double>;
inline K3 KeyOf(const vec3& p) { return {p.x, p.y, p.z}; }

}  // namespace e1

// engine seam pair-record (tri-tri intersection segment; `other` = partner)
struct E1Seam {
  vec3 p0, p1;
  int other;
};

// TWO-PHASE ENGINE (registry-first): buildPhase writes every split-producing
// event into the ONE per-line registry (no cells, no classification); the
// consume phase reads the registry ONLY (zero local reconciliation) and runs
// the walk/classify/emit machinery.  The wrapper iterates build to a global
// fixpoint, then consumes once.  seamCache: the ungated seam enumeration
// depends only on A, so the driver computes it once and every fixpoint
// round reuses it (it was re-enumerated per round - measured waste).
StageResult<Manifold::Impl> EmitCoordinatedBoundaryImpl(
    const Manifold::Impl& in, const BuildArrangement& A, double eps,
    std::map<std::tuple<int, int, int>, std::map<e1::K3, vec3>>& lineReg,
    std::vector<std::vector<E1Seam>>& seamCache, bool buildPhase) {
  using e1::K3;
  using e1::KeyOf;
  const int nTri = static_cast<int>(A.tri.size());
  const std::vector<vec3> seeds = WindingSeeds(in.bBox_);
  auto fail = [](const char* msg) {
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved, msg);
  };
  static const bool kDump = std::getenv("E1_DUMP") != nullptr;

  // E1_TIME: coarse per-phase wall clocks (stage-2 hot-spot attribution)
  static const bool kFlTime = std::getenv("E1_TIME") != nullptr;
  auto flNow = []() { return std::chrono::steady_clock::now(); };
  auto flMs = [](std::chrono::steady_clock::time_point a,
                 std::chrono::steady_clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  const auto flT0 = flNow();
  double tSeam = 0, tGroups = 0, tWalk = 0, tClassify = 0, tSolve = 0;
  // ---- 1. geometric plane groups (exact coplanarity union-find) ----
  std::vector<int> uf(nTri);
  for (int f = 0; f < nTri; ++f) uf[f] = f;
  std::function<int(int)> find = [&](int x) {
    while (uf[x] != x) x = uf[x] = uf[uf[x]];
    return x;
  };
  auto coplanarExact = [&](int i, int j) -> bool {
    for (int k = 0; k < 3; ++k) {
      const vec3& p = A.tri[j][k];
      // filter-first; a filter 0 (uncertain or true zero) escalates exact
      int s = Orient3DFilterSign(A.tri[i][0], A.tri[i][1], A.tri[i][2], p);
      if (s == 0)
        s = Orient3DExactSign(A.tri[i][0], A.tri[i][1], A.tri[i][2], p);
      if (s != 0) return false;
    }
    return true;
  };
  // SNAPPED-CLUSTER UNIONS (stage 3): faces of one near-coplanar cluster are
  // ONE semantic plane - the planarize snap's contract, carried by A.planeId.
  // The snap's projected coordinates are coplanar only to ~1e-13, so the
  // exact-coplanarity test below cannot re-derive the cluster; without this
  // union the engine emitted each cluster face as its own group and the
  // mutual footprints rang unpaired (measured: the near-coplanar folds'
  // 4+3 fan1 opens at SplitTouchingSheets).
  if (!A.planeId.empty()) {
    std::map<int, int> firstOfPlane;
    for (int f = 0; f < nTri; ++f) {
      const auto it = firstOfPlane.find(A.planeId[f]);
      if (it == firstOfPlane.end())
        firstOfPlane.emplace(A.planeId[f], f);
      else
        uf[find(f)] = find(it->second);
    }
  }
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (find(i) == find(j)) continue;
      const vec3 cr = la::cross(A.faceN[i], A.faceN[j]);
      const double nn = la::length(A.faceN[i]) * la::length(A.faceN[j]);
      if (la::length(cr) > 1e-9 * nn) continue;  // clearly non-parallel
      if (coplanarExact(i, j)) uf[find(i)] = find(j);
    }
  std::map<int, int> root2g;
  std::vector<int> gid(nTri, -1), rep;
  for (int f = 0; f < nTri; ++f) {
    const int r = find(f);
    auto it = root2g.find(r);
    if (it == root2g.end()) {
      it = root2g.emplace(r, static_cast<int>(rep.size())).first;
      rep.push_back(f);  // lowest face index = canonical rep
    }
    gid[f] = it->second;
  }
  const int nG = static_cast<int>(rep.size());
  std::vector<std::vector<int>> members(nG);
  std::vector<int> fsgn(nTri, 1);
  for (int f = 0; f < nTri; ++f) {
    members[gid[f]].push_back(f);
    fsgn[f] = la::dot(A.faceN[f], A.faceN[rep[gid[f]]]) >= 0.0 ? 1 : -1;
  }
  if (kDump)
    std::fprintf(stderr, "E1 groups=%d faces=%d junctions=%d\n", nG, nTri,
                 static_cast<int>(A.junctions.size()));

  // per-face float bboxes (probe-offset scan + crossing pre-filter)
  std::vector<Box> fbox(nTri);
  for (int f = 0; f < nTri; ++f) {
    Box b;
    for (int k = 0; k < 3; ++k) b.Union(A.tri[f][k]);
    fbox[f] = b;
  }
  // input-vertex identity set (canonV's never-snap rule), HOISTED: it
  // depends only on A.tri, and the former per-group rebuild dominated the
  // whole engine wall on the large self-intersector (measured 91%: 17k
  // groups x 51k ordered-set inserts per round)
  std::set<e1::K3> inputVerts;
  for (int f2 = 0; f2 < nTri; ++f2)
    for (int k = 0; k < 3; ++k) inputVerts.insert(KeyOf(A.tri[f2][k]));
  // the shared-collider broadphase over the component's triangles (the
  // proven exact box-overlap superset - the winding BVH reused for the
  // stack-face and near-sheet scans; candidates sorted for determinism)
  std::vector<std::array<vec3, 3>> flTriArr(nTri);
  for (int f = 0; f < nTri; ++f) flTriArr[f] = A.tri[f];
  const TriWindBVH flBvh = BuildTriWindBVH(flTriArr, in.bBox_);
  std::vector<int> flCands;
  auto flBoxQuery = [&](const vec3& lo, const vec3& hi) {
    WindCandidates(flBvh, lo, hi, flCands);
    std::sort(flCands.begin(), flCands.end());
  };

  // ---- 1b. ENGINE SEAMS: UNGATED exact tri-tri intersection segments ----
  // RecordSeams' production seams are representability-GATED (F11 sub-eps
  // collapse, phantom guard) - correct for the per-face path, but the
  // coordinated engine needs the FULL crossing structure: the tiny-dihedral
  // wedge chords the gates collapse are exactly the in-plane splits whose
  // absence the completeness certificate caught (a plane 7e-15 away at the
  // probe yet transversal at scale).  Enumerate exactly and ungated:
  //  - vertex side vs a group plane: filter-first exact orient3d against the
  //    group REP triangle (plane-consistent for every pair), cached per
  //    (vertex id, group);
  //  - a strictly-straddling edge contributes its pierce point, constructed
  //    ONCE per (edge, plane) key via the landed P3 SegPlaneBigHPoint and
  //    rounded once - byte-identical across all pairs and groups (the shared-
  //    double identity that makes cross-group sub-edges weld);
  //  - an exactly-on-plane vertex is itself a candidate (the measure-zero
  //    contact families: seam-through-vertex, edge-in-plane);
  //  - a candidate joins the seam iff INCLUSIVELY inside the OTHER triangle
  //    (input-exact IXOrient2D signs); 2 survivors = the seam segment, >2
  //    (degenerate contacts) = the extremes along the plane-pair direction.
  using ESeam = E1Seam;
  const auto flTSeam0 = flNow();
  const double scaleTop = in.bBox_.Scale();
  std::vector<vec3> cand;  // hoisted pair-candidate scratch (alloc churn)
  const bool seamCached = !seamCache.empty();
  if (!seamCached) seamCache.resize(nTri);
  std::vector<std::vector<ESeam>>& eSeams = seamCache;
  // hoisted memos (the completion pass constructs edge x seam crossings as
  // PIERCE identities - the exact once-only construction on BOTH lines)
  // packed-key hash memos (profiled: the ordered-map lookups per pair were
  // a top-of-profile cost on GT7081's dense near-tangent pair set)
  std::unordered_map<int64_t, int> sideCache;  // (vid<<32|gid) -> sign
  std::unordered_map<int64_t, vec3> pierceCache;
  auto sideKey = [](int vid2, int q) {
    return (static_cast<int64_t>(vid2) << 32) | static_cast<uint32_t>(q);
  };
  auto pierceKey = [](int v0, int v1, int q) {
    return (static_cast<int64_t>(v0) << 42) | (static_cast<int64_t>(v1) << 21) |
           q;
  };
  auto pierceGet = [&](int vlo, int vhi, const vec3& u, const vec3& w,
                       int q) -> std::optional<vec3> {
    const int64_t key = pierceKey(vlo, vhi, q);
    const auto it = pierceCache.find(key);
    if (it != pierceCache.end()) return it->second;
    const sos::BigHPoint X = sos::SegPlaneBigHPoint(u, w, A.tri[rep[q]].data());
    if (sos::BigSign(X.W) == 0) return std::nullopt;
    const vec3 P = sos::BigHPointToPos(X);
    if (!(std::isfinite(P.x) && std::isfinite(P.y) && std::isfinite(P.z)))
      return std::nullopt;
    pierceCache.emplace(key, P);
    return P;
  };
  if (!seamCached) {
    auto sideOf = [&](int f, int k, int q) -> int {
      const int64_t key = sideKey(A.vid[f][k], q);
      const auto it = sideCache.find(key);
      if (it != sideCache.end()) return it->second;
      const int r = rep[q];
      int s = Orient3DFilterSign(A.tri[r][0], A.tri[r][1], A.tri[r][2],
                                 A.tri[f][k]);
      if (s == 0)
        s = Orient3DExactSign(A.tri[r][0], A.tri[r][1], A.tri[r][2],
                              A.tri[f][k]);
      sideCache.emplace(key, s);
      return s;
    };
    auto bigInTriIncl = [&](const sos::BigHPoint& X, int t) -> bool {
      if (sos::BigSign(X.W) == 0) return false;
      const int ax = DominantAxis(A.faceN[t]);
      int pos = 0, neg = 0;
      for (int e = 0; e < 3; ++e) {
        const int o =
            IXOrient2D(sos::TrivialBigHPoint(A.tri[t][e]),
                       sos::TrivialBigHPoint(A.tri[t][(e + 1) % 3]), X, ax);
        if (o > 0) ++pos;
        if (o < 0) ++neg;
      }
      return !(pos && neg);
    };
    for (int f = 0; f < nTri; ++f) {
      for (int f2 = f + 1; f2 < nTri; ++f2) {
        if (gid[f2] == gid[f]) continue;
        const Box &ba = fbox[f], &bb = fbox[f2];
        if (ba.min.x > bb.max.x + eps || bb.min.x > ba.max.x + eps ||
            ba.min.y > bb.max.y + eps || bb.min.y > ba.max.y + eps ||
            ba.min.z > bb.max.z + eps || bb.min.z > ba.max.z + eps)
          continue;
        // EARLY PAIR REJECT (sound): a genuine tri-tri crossing needs a
        // plane straddle; if every vertex of one triangle sits STRICTLY one
        // side of the other's rep plane with a margin nine-plus orders above
        // the double dot's rounding, the pair has no seam.  (Profiled: the
        // per-pair memo/alloc constant on the dense GT7081 pair set was the
        // enumeration wall; most pairs are box-close but plane-separated.)
        auto farOneSide = [&](int fa, int fb) -> bool {
          const int q = gid[fb];
          const vec3& nq = A.faceN[rep[q]];
          const vec3& pq = A.tri[rep[q]][0];
          const double M = 1e-6 * (1.0 + scaleTop) * la::length(nq);
          int abv = 0, blw = 0;
          for (int k = 0; k < 3; ++k) {
            const double d = la::dot(nq, A.tri[fa][k] - pq);
            if (d > M)
              ++abv;
            else if (d < -M)
              ++blw;
          }
          return abv == 3 || blw == 3;
        };
        if (farOneSide(f, f2) || farOneSide(f2, f)) continue;
        cand.clear();
        auto addCand = [&](const vec3& p) {
          for (const auto& c : cand)
            if (KeyOf(c) == KeyOf(p)) return;
          cand.push_back(p);
        };
        auto collect = [&](int fa, int fb) {  // fa's edges vs plane(gid[fb])
          const int q = gid[fb];
          for (int e = 0; e < 3; ++e) {
            const int e2 = (e + 1) % 3;
            const int s0 = sideOf(fa, e, q), s1 = sideOf(fa, e2, q);
            if (s0 == 0 &&
                bigInTriIncl(sos::TrivialBigHPoint(A.tri[fa][e]), fb))
              addCand(A.tri[fa][e]);
            if (s1 == 0 &&
                bigInTriIncl(sos::TrivialBigHPoint(A.tri[fa][e2]), fb))
              addCand(A.tri[fa][e2]);
            if (s0 != 0 && s1 != 0 && s0 != s1) {
              int v0 = A.vid[fa][e], v1 = A.vid[fa][e2];
              vec3 u = A.tri[fa][e], w = A.tri[fa][e2];
              if (v0 > v1) {
                std::swap(v0, v1);
                std::swap(u, w);
              }
              const int64_t key = pierceKey(v0, v1, q);
              auto it = pierceCache.find(key);
              if (it == pierceCache.end()) {
                const sos::BigHPoint X =
                    sos::SegPlaneBigHPoint(u, w, A.tri[rep[q]].data());
                if (sos::BigSign(X.W) == 0) continue;  // parallel: no pierce
                it = pierceCache.emplace(key, sos::BigHPointToPos(X)).first;
              }
              const vec3& P = it->second;
              if (!(std::isfinite(P.x) && std::isfinite(P.y) &&
                    std::isfinite(P.z)))
                continue;
              // COARSE INCLUSION PRE-CHECK on the once-rounded P: the three
              // dropped-frame edge dets in doubles.  P is within ~2^-52
              // relative of the exact pierce, so each det's error is below
              // ~1e-15*scale^2; the margin sits NINE-plus orders above it -
              // a verdict outside the margin band cannot flip.  Only the
              // margin band pays the Big construction + exact inclusion
              // (the per-pair Big tests were the GT7081 seam wall: 177s).
              const int axF = DominantAxis(A.faceN[fb]);
              const double mIncl =
                  1e-6 * (1.0 + in.bBox_.Scale()) * (1.0 + in.bBox_.Scale());
              int cPos = 0, cNeg = 0;
              double dMin = std::numeric_limits<double>::infinity();
              for (int ee = 0; ee < 3; ++ee) {
                const vec2 a2 = e1::Drop2(A.tri[fb][ee], axF);
                const vec2 b2 = e1::Drop2(A.tri[fb][(ee + 1) % 3], axF);
                const vec2 p2 = e1::Drop2(P, axF);
                const double det = la::cross(b2 - a2, p2 - a2);
                dMin = std::min(dMin, std::abs(det));
                if (det > 0.0) ++cPos;
                if (det < 0.0) ++cNeg;
              }
              if (dMin > mIncl) {
                if (!(cPos && cNeg)) addCand(P);
                continue;
              }
              const sos::BigHPoint X =
                  sos::SegPlaneBigHPoint(u, w, A.tri[rep[q]].data());
              if (bigInTriIncl(X, fb)) addCand(P);
            }
          }
        };
        collect(f, f2);
        collect(f2, f);
        if (cand.size() < 2) continue;  // point contact / no crossing
        vec3 p0 = cand[0], p1 = cand[1];
        if (cand.size() > 2) {  // degenerate contact: extremes along the line
          const vec3 dir =
              la::cross(A.faceN[rep[gid[f]]], A.faceN[rep[gid[f2]]]);
          double lo = la::dot(cand[0], dir), hi = lo;
          for (const vec3& c : cand) {
            const double t = la::dot(c, dir);
            if (t < lo) {
              lo = t;
              p0 = c;
            }
            if (t > hi) {
              hi = t;
              p1 = c;
            }
          }
          if (KeyOf(p0) == KeyOf(p1)) continue;
        }
        eSeams[f].push_back({p0, p1, f2});
        eSeams[f2].push_back({p0, p1, f});
      }
    }
    if (const char* ef = std::getenv("E1_ESEAMS")) {
      const int ff = std::atoi(ef);
      for (const ESeam& s : eSeams[ff])
        std::fprintf(stderr,
                     "E1 ESEAM f=%d other=%d gid=%d p0=(%.10g,%.10g,%.10g) "
                     "p1=(%.10g,%.10g,%.10g)\n",
                     ff, s.other, gid[s.other], s.p0.x, s.p0.y, s.p0.z, s.p1.x,
                     s.p1.y, s.p1.z);
    }
    if (kDump) {
      int nseam = 0;
      for (const auto& v : eSeams) nseam += static_cast<int>(v.size());
      std::fprintf(stderr, "E1 engine seams=%d (pair-records)\n", nseam);
    }
  }

  tSeam += flMs(flTSeam0, flNow());
  // canonical engine triple positions, once per sorted gid triple
  std::map<std::array<int, 3>, std::pair<bool, vec3>> tripleCache;
  auto triplePos = [&](int g0, int g1, int g2, vec3& pos) -> bool {
    std::array<int, 3> key = {g0, g1, g2};
    std::sort(key.begin(), key.end());
    auto it = tripleCache.find(key);
    if (it == tripleCache.end()) {
      vec3 p;
      const bool ok =
          CanonTriplePos(A, rep[key[0]], rep[key[1]], rep[key[2]], p);
      it = tripleCache.emplace(key, std::make_pair(ok, p)).first;
    }
    pos = it->second.second;
    return it->second.first;
  };
  // input-exact strict interior of triple {g0,g1,g2} in triangle t
  auto tripleInTri = [&](int g0, int g1, int g2, int t) -> bool {
    const sos::BigHPoint X = IXTripleHPoint(A, rep[g0], rep[g1], rep[g2]);
    if (sos::BigSign(X.W) == 0) return false;
    const int axis = DominantAxis(A.faceN[t]);
    int o[3];
    for (int e = 0; e < 3; ++e) {
      o[e] = IXOrient2D(sos::TrivialBigHPoint(A.tri[t][e]),
                        sos::TrivialBigHPoint(A.tri[t][(e + 1) % 3]), X, axis);
      if (o[e] == 0) return false;
    }
    return o[0] == o[1] && o[1] == o[2];
  };

  std::vector<OutTri3D> out;
  std::vector<int> outG;   // per-tri emitting group (diagnostics)
  std::vector<vec3> outN;  // per-tri ORIENTED sheet normal (radial branch)
  int suspectTotal = 0;    // snap-grid identity audit (owner invariant 2)
  // MEMOIZED CONSTRUCTION (owner directive; the campaign's memoize-values
  // principle = the stage-2 shape): every constructed point is committed
  // ONCE under its canonical identity; every later path LOOKS IT UP -
  // bit-equality by construction.  Identities: seam endpoints (edge vids,
  // plane gid) -> pierceCache; triples (sorted gid triple) -> tripleCache;
  // pair crossings (unordered segment-identity pair) -> crossMemo; feet
  // (vertex bits, segment identity) -> footMemo.  A segment's identity is
  // its endpoint bit-pair (endpoints are themselves committed identities).
  using SegKey = std::pair<e1::K3, e1::K3>;
  auto segKeyOf = [](const vec3& p0, const vec3& p1) -> SegKey {
    const e1::K3 a = e1::KeyOf(p0), b = e1::KeyOf(p1);
    return a < b ? SegKey{a, b} : SegKey{b, a};
  };
  std::map<std::pair<e1::K3, SegKey>, vec3> footMemo;
  // registry provenance (E1_DUMP diagnostics): position bits -> producer tag
  std::map<e1::K3, const char*> regProv;
  const char* curProducer = "?";
  int dustTri = 0, triFail = 0, spliceFail = 0;
  int nCells = 0, nNeg = 0, nJump = 0, nBoundary = 0, nOwned = 0, nDustCell = 0;
  int nPancake = 0;
  const double scale = in.bBox_.Scale();

  // ==== FLOOD WINDING FIELD (flip arc stage 1) ==============================
  // The component-global integer winding field over the arrangement cells.
  // Field value E(c) = the EXACT winding of the epsilon-layer immediately on
  // the +nHat side of cell c's group plane.  The eps-layer is the right field
  // variable because its lateral deltas are purely COMBINATORIAL: crossing a
  // sub-edge crosses exactly the chord-partner sheets ON that sub-edge (a
  // sheet without a chord there cannot separate the two eps-layers - the
  // mesh-edge cleanliness argument, exact as eps -> 0), with sign
  // -sgn(dot(t_out, n_partner)).  The per-cell probe value converts
  // VERTICALLY: wA = E - sPosMid (sPosMid = signed count of covering foreign
  // sheets crossed by the cenP vertical between the plane and the probe at
  // 2T; the gap-finder certifies that set at cenP).  Anchors: (a) the HULL
  // SEED - at the lex-max input vertex v* the point v* + (d, d^2, d^3) is
  // exterior (w=0, no ray cast: the lex-max of a compact polyhedron is a
  // vertex) and lies on the lexsign side of every plane through v*, so a
  // group with exactly one node cell touching v* anchors combinatorially;
  // (b) residual probes (one certified probe anchors each subgraph).  BFS
  // propagates over RELIABILITY-FILTERED edges; a conflict clears the field
  // and the emission probes every cell (the pre-flood per-cell semantics).
  // FLOOD-PRIMARY (session 3): the field values every cell; the exact probe
  // remains as (a) the residual subgraph anchor, (b) the per-cell
  // specialist where edge transport is unreliable (the near-tangent ladder
  // class), (c) the E1_FLOODDIFF shadow validator.
  struct FlNode {
    int g = -1;
    int jump = 0, ownJump = 0;
    int sPosMid = 0;       // covering sheets with stackWin < tPos < 2T at cenP
    int sZero = 0;         // covering foreign sheets with tPos == 0 at cenP
    bool hasBand = false;  // any covering sheet inside the sub-resolution
                           // band (|tPos| <= stackWin): its side is not
                           // double-resolvable (the engine's stack class)
    bool pancake = false;
    bool emit = false;  // carries emission payload (jump!=0 or pancake)
    // probeState: 0 unprobed, 1 certified (probeWA/WB valid), 2 grazed,
    // 3 probed but certificate failed
    signed char probeState = 0;
    // anchor-grade: certified at the BASE offset with no retry/acceptance -
    // ladder-retried readings are trusted for the cell's own emission but
    // not as the global field base (they sit in delicate corridors)
    bool anchorOK = false;
    bool certified = false;  // == (probeState == 1), kept for the census
    int probeWA = 0, probeWB = 0;
    vec3 cenP;        // projected probe center (residual probing)
    double off = 0;   // probe offset 2T
    double extW = 0;  // largest-subtri altitude (the dust-width rule)
    vec3 Nrep;
    std::vector<std::pair<double, int>> stackSheets;  // pancake layer cert
    // covering foreign sheets at cenP: (face, state, sgn(dot(n,nHat)))
    // with state +1 ABOVE (tPos > stackWin), 0 MID (in-band), -1 BELOW -
    // sorted by face; the field level E' sits just above the MID band, so
    // only ABOVE-state changes are lateral crossings (band sides are noise)
    std::vector<std::array<int, 3>> nearSign;
    std::vector<ivec3> tris;  // emission payload
    std::map<int, vec3> pos;  // emission payload
  };
  std::vector<FlNode> flNodes;
  // (n0, n1, d, tag): E(n1) = E(n0) + d; tag = contributing-arm bitmask
  // (1 straddle-chord, 2 touch-chord, 4 ridge, 8 parity, 32 handoff,
  // 64 handoff-leg-correction) - the differential's arm-level census
  std::vector<std::array<int, 4>> flEdges;
  // handoff records: (input vid pair, sub-edge position bits) ->
  // (node, eOff) with the invariant E(node) + eOff equal on both sides
  std::map<std::pair<std::pair<int, int>, SegKey>,
           std::vector<std::pair<int, int>>>
      flHand;
  std::vector<std::pair<int, int>> flSeeds;  // (node, anchored E value)
  auto flLexsign = [](const vec3& n) -> int {
    if (n.x != 0.0) return n.x > 0.0 ? 1 : -1;
    if (n.y != 0.0) return n.y > 0.0 ? 1 : -1;
    if (n.z != 0.0) return n.z > 0.0 ? 1 : -1;
    return 0;
  };
  static const bool kFloodDiff = std::getenv("E1_FLOODDIFF") != nullptr;
  // v* = lex-max input vertex (position identity; ties collapse to one K3)
  K3 flStarKey{0, 0, 0};
  {
    bool have = false;
    for (size_t v = 0; v < in.vertPos_.size(); ++v) {
      const K3 k = KeyOf(in.vertPos_[v]);
      if (!have || flStarKey < k) {
        flStarKey = k;
        have = true;
      }
    }
  }

  for (int g = 0; g < nG; ++g) {
    const auto flTG0 = flNow();
    int gProbeFail = 0, gProbeCert = 0;  // flood graze census
    // flood graph collection (this group's arrangement):
    // sub-edge (lo,hi vertex ids) -> contributing segment indices
    std::map<std::pair<int, int>, std::vector<int>> flEdgeSegs;
    // sub-edge -> (node, traversed lo->hi?) per adjacent node cell
    std::map<std::pair<int, int>, std::vector<std::pair<int, bool>>>
        flEdgeCells;
    std::vector<int> flStarNodes;  // node cells touching v* in this group
    int flVstarLocal = -1;         // v*'s walk-graph id in this group
    // ---- 2. the group's segment set (3D endpoint pairs, shared doubles) ----
    struct Seg {
      vec3 p0, p1;
      int planeQ;  // partner group (-1 = member triangle edge)
      int fOwn, fOther;
      int vidLo = -1, vidHi = -1;  // mesh-edge line identity (planeQ < 0)
    };
    std::vector<Seg> segs;
    for (const int f : members[g]) {
      for (int e = 0; e < 3; ++e) {
        const int v0 = A.vid[f][e], v1 = A.vid[f][(e + 1) % 3];
        // vids POSITION-MATCHED to p0/p1 (pierce identities need the
        // canonical vid order with matching positions); sorted at use sites
        segs.push_back({A.tri[f][e], A.tri[f][(e + 1) % 3], -1, f, -1, v0, v1});
      }
      for (const ESeam& s : eSeams[f])
        segs.push_back({s.p0, s.p1, gid[s.other], f, s.other});
    }
    const vec3 Nrep = A.faceN[rep[g]];
    const int axis = DominantAxis(Nrep);
    const vec3 nHat = Nrep / la::length(Nrep);
    // STACK (sub-double-resolution near-coincident sheets): a foreign face
    // whose plane sits closer than any double probe can separate is treated
    // as part of ONE effective sheet with this group: its EDGES join the
    // arrangement (as shadow segments - the same input-vert doubles on both
    // groups, so the cross-group boundary welds), its covering sign joins the
    // NET jump, and only the LOWEST gid of the covering stack emits.  (The
    // exact-rational offline engine probed between such sheets; doubles
    // cannot - the certificate measured planes 1.7e-14 apart, below ULP.)
    const double stackWin =
        std::max(eps / 100.0, 128.0 * std::numeric_limits<double>::epsilon() *
                                  (1.0 + scale));
    Box gbox;
    for (const int f : members[g]) gbox.Union(fbox[f]);
    std::vector<int> stackFaces;
    for (int f2 = 0; f2 < nTri; ++f2) {
      if (gid[f2] == g) continue;
      const Box& b = fbox[f2];
      if (b.min.x > gbox.max.x + eps || b.max.x < gbox.min.x - eps ||
          b.min.y > gbox.max.y + eps || b.max.y < gbox.min.y - eps ||
          b.min.z > gbox.max.z + eps || b.max.z < gbox.min.z - eps)
        continue;
      bool near = true;
      for (int k = 0; k < 3 && near; ++k)
        near = std::abs(la::dot(Nrep, A.tri[f2][k] - A.tri[rep[g]][0])) /
                   la::length(Nrep) <=
               stackWin;
      if (!near) continue;
      stackFaces.push_back(f2);
      for (int e = 0; e < 3; ++e) {
        const int v0 = A.vid[f2][e], v1 = A.vid[f2][(e + 1) % 3];
        segs.push_back(
            {A.tri[f2][e], A.tri[f2][(e + 1) % 3], -1, f2, -1, v0, v1});
      }
    }
    const int nS = static_cast<int>(segs.size());

    // ---- BUILD PHASE ONLY: the split producers (crossing enumeration,
    // T-junction pool, foot exchange, planarity completion) write into the
    // registry via addSplit; the consume phase reads the registry only.
    // ---- 3. splits: engine triples (input-exact, symmetric) + registry ----
    std::vector<std::vector<std::pair<double, vec3>>> splits(nS);
    // CANONICAL VERTEX GRID (rounding-scale): every split point is snapped
    // onto any existing vertex within ~64 ULP before entering the
    // arrangement.  Derivation noise (retry interpolations, triple-vs-drawn
    // constructions) otherwise accumulates ULP-VARIANT chains of the same
    // geometric vertex with degenerate sliver cells between them; the weld
    // drops those slivers and their pairing halfedges (measured: four
    // 5e-16-apart variants of one vertex around an open edge).  The radius is
    // far below representable structure - this collapses noise, not geometry.
    // OWNER INVARIANTS (snap-grid review): (1) the snap applies to
    // CONSTRUCTED vertices only - INPUT vertices are excluded (identity =
    // the given bits / the F11 position rule; snapping inputs would move
    // real data); (2) the grid must never merge two vertices of DIFFERENT
    // canonical identity - audited cheaply below: snaps beyond the ~few-ULP
    // derivation-noise band are counted and reported under E1_DUMP (any such
    // event indicates an identity-keying bug, not a tolerance issue).  The
    // 64-ULP radius is a margin INSIDE the measured safe band (noise ~few
    // ULP; smallest real structure 0.586 eps and ~1e-8 pairs, orders above);
    // an eps-scale snap would merge real structure and is banned.
    const double rho =
        64.0 * std::numeric_limits<double>::epsilon() * (1.0 + scale);
    const double rhoNoise =
        8.0 * std::numeric_limits<double>::epsilon() * (1.0 + scale);
    int suspectSnaps = 0;
    std::map<std::tuple<long long, long long, long long>, std::vector<vec3>>
        canonGrid;
    auto gridInsert = [&](const vec3& v) {
      canonGrid[{static_cast<long long>(std::floor(v.x / rho)),
                 static_cast<long long>(std::floor(v.y / rho)),
                 static_cast<long long>(std::floor(v.z / rho))}]
          .push_back(v);
    };
    auto canonV = [&](const vec3& v) -> vec3 {
      if (inputVerts.count(KeyOf(v))) return v;  // inputs are never snapped
      const long long cx = static_cast<long long>(std::floor(v.x / rho));
      const long long cy = static_cast<long long>(std::floor(v.y / rho));
      const long long cz = static_cast<long long>(std::floor(v.z / rho));
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz) {
            const auto it = canonGrid.find({cx + dx, cy + dy, cz + dz});
            if (it == canonGrid.end()) continue;
            for (const vec3& w : it->second) {
              const double d = la::length(w - v);
              if (d <= rho) {
                if (d > rhoNoise) ++suspectSnaps;  // identity audit trail
                return w;
              }
            }
          }
      gridInsert(v);
      return v;
    };
    // SEED + ENDPOINT ADOPTION, two passes.  Pass 1: input endpoints are the
    // unconditional anchors (identity = the given bits, never snapped).
    // Pass 2: constructed endpoints ADOPT an existing anchor within the
    // rounding band or become one - and the record is REWRITTEN to the
    // adopted bits.  Seeding without rewriting left a junction represented
    // TWICE: a pierce constructed within a few ULP of an input vertex kept
    // its own rounding as a segment endpoint, so one line's chain ended at
    // the input bits while the crossing line's chain ran through the pierce
    // bits with NO adjacency between them - the rotation walk then sews the
    // two fans into one self-overlapping macro cell whose boundary is
    // oracle-false (measured: the (-17.9,0.6,-204) junction, twins 3.3e-15
    // apart, the whole z=-204 over/under-emission family).  Input-first
    // ordering makes the adopted identity deterministic and input-preferring.
    for (const Seg& s : segs) {
      if (inputVerts.count(KeyOf(s.p0))) gridInsert(s.p0);
      if (inputVerts.count(KeyOf(s.p1))) gridInsert(s.p1);
    }
    for (Seg& s : segs) {
      if (!inputVerts.count(KeyOf(s.p0))) s.p0 = canonV(s.p0);
      if (!inputVerts.count(KeyOf(s.p1))) s.p1 = canonV(s.p1);
    }
    // PER-LINE REGISTRY key: seams by plane-gid pair; mesh edges by vid pair
    // (one geometric line = one split list, shared by every segment record
    // and every group on it)
    auto lineKeyOf = [&](int si) -> std::tuple<int, int, int> {
      if (segs[si].planeQ >= 0) {
        const auto pq = std::minmax(g, segs[si].planeQ);
        return {1, pq.first, pq.second};
      }
      return {0, std::min(segs[si].vidLo, segs[si].vidHi),
              std::max(segs[si].vidLo, segs[si].vidHi)};
    };
    // ENDPOINT REGISTRATION (build): every segment record's endpoints are
    // committed identities (pierces / input verts); registering them on the
    // line lets every OTHER record of the line split at record boundaries -
    // otherwise a longer record spans past a shorter one's end and the
    // handoff stub opens (the measured residue class).
    if (buildPhase)
      for (int i = 0; i < nS; ++i) {
        auto& ent = lineReg[lineKeyOf(i)];
        ent.emplace(KeyOf(segs[i].p0), segs[i].p0);
        ent.emplace(KeyOf(segs[i].p1), segs[i].p1);
      }
    auto addSplit = [&](int si, const vec3& Vraw) -> bool {
      const vec3 V = canonV(Vraw);
      const Seg& s = segs[si];
      const vec3 d3 = s.p1 - s.p0;
      const double len2 = la::dot(d3, d3);
      if (!(len2 > 0.0)) return false;
      const double t = la::dot(V - s.p0, d3) / len2;
      // interiority at ROUNDING scale, not eps: rejecting a split within eps
      // of an endpoint leaves a REAL crossing in the drawn graph (non-planar
      // walk - measured in the near-duplicate zigzag zones where crossings
      // crowd the endpoints); splitting there merely creates a sub-eps
      // sliver sub-edge that welds away.
      const double tlo = 64.0 * std::numeric_limits<double>::epsilon() *
                         (1.0 + scale) / std::sqrt(len2);
      if (!(t > tlo && t < 1.0 - tlo)) return false;  // strictly interior
      for (const auto& pr : splits[si])
        if (KeyOf(pr.second) == KeyOf(V)) return false;
      splits[si].push_back({t, V});
      if (buildPhase) lineReg[lineKeyOf(si)].emplace(KeyOf(V), V);
      if (buildPhase && kDump) regProv.emplace(KeyOf(V), curProducer);
      return true;
    };
    // PER-LINE REGISTRY READ (both phases): every line's accumulated split
    // union - seams by plane pair, mesh edges by vid pair - lands on every
    // segment record of that line, across groups.  ON-LINE GUARD: one
    // plane-pair key can span near-tangent lens records microns apart; an
    // entry only lands where it lies on THIS segment's line.
    for (int i = 0; i < nS; ++i) {
      const auto it = lineReg.find(lineKeyOf(i));
      if (it == lineReg.end()) continue;
      const Seg& s = segs[i];
      const vec3 d3 = s.p1 - s.p0;
      const double len2 = la::dot(d3, d3);
      if (!(len2 > 0.0)) continue;
      for (const auto& kv : it->second) {
        const vec3& V = kv.second;
        const double t = la::dot(V - s.p0, d3) / len2;
        if (!(t > 0.0 && t < 1.0)) continue;
        if (la::length(V - s.p0 - t * d3) > 2.0 * rho) continue;
        addSplit(i, V);
      }
    }
    if (buildPhase) {
      for (int i = 0; i < nS; ++i) {
        if (segs[i].planeQ < 0) continue;
        for (int j = i + 1; j < nS; ++j) {
          if (segs[j].planeQ < 0 || segs[j].planeQ == segs[i].planeQ) continue;
          // quick reject: 2D bboxes of the two seams disjoint
          const vec2 a0 = e1::Drop2(segs[i].p0, axis),
                     a1 = e1::Drop2(segs[i].p1, axis),
                     b0 = e1::Drop2(segs[j].p0, axis),
                     b1 = e1::Drop2(segs[j].p1, axis);
          if (std::max(a0.x, a1.x) < std::min(b0.x, b1.x) - eps ||
              std::max(b0.x, b1.x) < std::min(a0.x, a1.x) - eps ||
              std::max(a0.y, a1.y) < std::min(b0.y, b1.y) - eps ||
              std::max(b0.y, b1.y) < std::min(a0.y, a1.y) - eps)
            continue;
          // exact symmetric existence: X = {g, Qi, Qj} strictly interior to all
          // four bounding triangles (both seams' extents)
          if (!tripleInTri(g, segs[i].planeQ, segs[j].planeQ, segs[i].fOwn))
            continue;
          if (!tripleInTri(g, segs[i].planeQ, segs[j].planeQ, segs[i].fOther))
            continue;
          if (segs[j].fOwn != segs[i].fOwn &&
              !tripleInTri(g, segs[i].planeQ, segs[j].planeQ, segs[j].fOwn))
            continue;
          if (!tripleInTri(g, segs[i].planeQ, segs[j].planeQ, segs[j].fOther))
            continue;
          vec3 X;
          if (!triplePos(g, segs[i].planeQ, segs[j].planeQ, X)) continue;
          addSplit(i, X);
          addSplit(j, X);
        }
      }
      // T-junction splits (shared canonical doubles): the group's OWN segment
      // endpoints (engine seam ends + member corners - the offline pool) plus
      // the production registry; the same on-line rule everywhere, consistent
      // across groups because every group reads identical (vertex, segment)
      // bits.
      // pool = the group's OWN segment endpoints (exact once-only
      // constructions, ULP-accurate - the offline pool).  NOT A.junctions: its
      // eps-deduped representatives sit up to eps off the exact lines, and
      // inserting them with an eps window BENDS the chains by eps - the bent
      // sub-chains then properly cross and fragment the walk (measured: a split
      // cascade). OUTER FIXPOINT: the exchange and the completion feed each
      // other - a completion-added crossing (ill-conditioned on near-collinear
      // pairs: the same line's overlapping segments each get their own noisy
      // crossing position) must be EXCHANGED onto every collinear twin, and
      // exchanged points can expose new crossings.  (Measured: identical-
      // endpoint twin segments carrying different completion splits 2.4e-5
      // apart - the T-junction/lens class.)
      for (int outer = 0; outer < 4; ++outer) {
        size_t nsplit0 = 0;
        for (int i = 0; i < nS; ++i) nsplit0 += splits[i].size();
        // T-junction pool pass, iterated to FIXPOINT with ALL SPLIT POINTS in
        // the pool: near-collinear overlapping chains (the same line reached
        // via different member pairs, ULP apart) must carry IDENTICAL
        // subdivisions - endpoint-only exchange left one chain split where its
        // twin spanned whole (measured: intra-group T-junctions with the
        // on-vertex 9e-16 off the unsplit edge).  Split values are canonical
        // (snap grid), so the exchange converges. 8*rho: the canonical snap
        // grid moves split points up to rho off their segment lines, so on-line
        // tests must budget the snap displacement
        const double tolLine = 2.0 * rho;
        for (int round = 0; round < 4; ++round) {
          std::vector<vec3> pool;
          pool.reserve(2 * segs.size());
          for (const Seg& s : segs) {
            pool.push_back(s.p0);
            pool.push_back(s.p1);
          }
          for (int i = 0; i < nS; ++i)
            for (const auto& pr : splits[i]) pool.push_back(pr.second);
          bool added = false;
          for (int i = 0; i < nS; ++i) {
            const Seg& s = segs[i];
            const vec3 d3 = s.p1 - s.p0;
            const double len2 = la::dot(d3, d3);
            if (!(len2 > 0.0)) continue;
            const double len = std::sqrt(len2);
            for (const vec3& V : pool) {
              const vec3 w = V - s.p0;
              const double t = la::dot(w, d3) / len2;
              if (!(t > eps / len && t < 1.0 - eps / len)) continue;
              if (la::length(w - t * d3) > tolLine) continue;
              added |= addSplit(i, V);
            }
          }
          if (!added) break;
        }

        // PLANARITY COMPLETION (all remaining segment-pair crossings): the
        // exact seam-x-seam enumeration and the endpoint pool cover the
        // canonical crossings, but the drawn (rounded) graph must be PLANAR for
        // the face walk - overlapping coplanar members (the fold structure)
        // cross member edges and same-line seams in ways the passes above miss
        // (measured: properly-crossing sub-edges -> bowtie walks ->
        // untriangulable cells). Detect every remaining proper crossing exactly
        // on the shared rounded endpoints and split both segments; the split
        // point uses the canonical triple when both carriers are seams of
        // distinct planes, else the in-segment interpolation (identity across
        // groups holds within weld tolerance via the shared-endpoint
        // constructions).
        for (int i = 0; i < nS; ++i) {
          const vec2 a0 = e1::Drop2(segs[i].p0, axis),
                     a1 = e1::Drop2(segs[i].p1, axis);
          for (int j = i + 1; j < nS; ++j) {
            const vec2 b0 = e1::Drop2(segs[j].p0, axis),
                       b1 = e1::Drop2(segs[j].p1, axis);
            if (std::max(a0.x, a1.x) < std::min(b0.x, b1.x) - eps ||
                std::max(b0.x, b1.x) < std::min(a0.x, a1.x) - eps ||
                std::max(a0.y, a1.y) < std::min(b0.y, b1.y) - eps ||
                std::max(b0.y, b1.y) < std::min(a0.y, a1.y) - eps)
              continue;
            if (!e1::ProperCross2(segs[i].p0, segs[i].p1, segs[j].p0,
                                  segs[j].p1, axis))
              continue;
            vec3 X;
            bool have = false;
            if (segs[i].planeQ >= 0 && segs[j].planeQ >= 0 &&
                segs[i].planeQ != segs[j].planeQ)
              have = triplePos(g, segs[i].planeQ, segs[j].planeQ, X);
            if (!have) {
              // EXACT-CONSTRUCTION crossings (the interpolated-memo
              // refutation is the design constraint): an interpolated point
              // lies exactly on ONE line only, so wherever a committed
              // identity exists, use it -
              //   edge x seam  -> the PIERCE of the edge through the seam's
              //                   partner plane (once-only, on BOTH lines)
              //   edge x edge / seam-twin pairs -> canonical interpolation
              //                   (deterministic carrier order, bit-equal
              //                   across groups)
              const bool iEdge = segs[i].planeQ < 0, jEdge = segs[j].planeQ < 0;
              bool built = false;
              if (iEdge != jEdge) {
                const int se = iEdge ? i : j;  // the edge carrier
                const int ss = iEdge ? j : i;  // the seam carrier
                const int q = segs[ss].planeQ;
                const bool fwd = segs[se].vidLo <= segs[se].vidHi;
                const int vlo = fwd ? segs[se].vidLo : segs[se].vidHi;
                const int vhi = fwd ? segs[se].vidHi : segs[se].vidLo;
                const vec3& u = fwd ? segs[se].p0 : segs[se].p1;
                const vec3& w = fwd ? segs[se].p1 : segs[se].p0;
                if (const auto P = pierceGet(vlo, vhi, u, w, q)) {
                  X = *P;
                  built = true;
                }
              }
              if (!built) {
                // canonical carrier: the smaller segment identity
                const int ci = segKeyOf(segs[i].p0, segs[i].p1) <=
                                       segKeyOf(segs[j].p0, segs[j].p1)
                                   ? i
                                   : j;
                const int cj = ci == i ? j : i;
                const vec2 c0 = e1::Drop2(segs[ci].p0, axis),
                           c1 = e1::Drop2(segs[ci].p1, axis);
                const vec2 e0 = e1::Drop2(segs[cj].p0, axis),
                           e1v = e1::Drop2(segs[cj].p1, axis);
                const double dax = c1.x - c0.x, day = c1.y - c0.y;
                const double dbx = e1v.x - e0.x, dby = e1v.y - e0.y;
                const double den = dax * dby - day * dbx;
                if (den == 0.0) continue;
                const double t =
                    ((e0.x - c0.x) * dby - (e0.y - c0.y) * dbx) / den;
                X = segs[ci].p0 + t * (segs[ci].p1 - segs[ci].p0);
              }
            }
            addSplit(i, X);
            addSplit(j, X);
          }
        }

        size_t nsplit1 = 0;
        for (int i = 0; i < nS; ++i) nsplit1 += splits[i].size();
        if (nsplit1 == nsplit0) break;
      }
      suspectTotal += suspectSnaps;
      tGroups += flMs(flTG0, flNow());
      continue;  // build phase: no cells, no classification
    }
    // NOTE: no sub-edge-level completion pass is needed: with the exact
    // endpoint pool and the rounding-scale on-line tolerance, chain bends are
    // ~ULP, so residual sub-edge crossings are ULP-scale bowties the
    // triangulation's weld-dust remainder acceptance absorbs (a full
    // sub-edge crossing fixpoint was measured to CASCADE: each round's
    // interpolated splits create new bent sub-edges - and cost minutes).

    if (const char* sd = std::getenv("E1_SEGDUMP")) {
      int pg = -1;
      double pu = 0, pv = 0, pr = 0;
      if (std::sscanf(sd, "%d,%lf,%lf,%lf", &pg, &pu, &pv, &pr) == 4 &&
          pg == g) {
        const vec2 P{pu, pv};
        const double pr2r = pr;
        for (int i = 0; i < nS; ++i) {
          const vec2 a = e1::Drop2(segs[i].p0, axis),
                     b = e1::Drop2(segs[i].p1, axis);
          // distance from P to segment ab
          const vec2 d2 = b - a;
          const double L2 = la::dot(d2, d2);
          double t = L2 > 0 ? la::dot(P - a, d2) / L2 : 0.0;
          t = std::max(0.0, std::min(1.0, t));
          if (la::length(P - (a + t * d2)) > pr) continue;
          std::fprintf(stderr,
                       "E1 SEG g=%d i=%d planeQ=%d fOwn=%d fOther=%d "
                       "a=(%.9g,%.9g) b=(%.9g,%.9g) splits=%d\n",
                       g, i, segs[i].planeQ, segs[i].fOwn, segs[i].fOther, a.x,
                       a.y, b.x, b.y, static_cast<int>(splits[i].size()));
          for (const auto& pr : splits[i]) {
            const vec2 q = e1::Drop2(pr.second, axis);
            if (la::length(q - P) <= 8.0 * pr2r)
              std::fprintf(stderr, "      split t=%.6f (%.12g,%.12g)\n",
                           pr.first, q.x, q.y);
          }
        }
      }
    }
    tGroups += flMs(flTG0, flNow());
    const auto flTW0 = flNow();
    // ---- 4. 2D graph (verts keyed by 3D bits) + exact rotation walk ----
    // SUB-RHO JUNCTION CLUSTERING: distinct committed identities inside the
    // rounding band (unmergeable anchors - e.g. an input twin pair 3.3e-15
    // apart in this corpus) become ONE walk-graph node.  Their chain systems
    // can properly CROSS within the band (measured: 4.7e-13 from the twin
    // endpoints - below the split-interiority floor, so the crossing is
    // unrepresentable as a graph vertex); a graph that keeps the copies
    // separate is then NON-PLANAR and the rotation walk sews the fans into
    // self-overlapping macro cells whose classification paints oracle-false
    // boundary (measured: the (-17.9,0.6,-204) input twin pair, the whole
    // z=-204 over/under-emission family).  Collapsing the band collapses the
    // crossing INTO the junction, exactly where the exact arrangement's
    // structure lands once rounded.  The REGISTRY identities stay distinct
    // (this is emission topology, not an identity merge); the node's
    // representative is a committed identity (first-seen anchor, inputs
    // first via the endpoint-adoption pass), and the assembly weld (radius
    // eps, two orders above rho) identifies the pair in the output
    // regardless.
    std::map<K3, int> vidOf;
    std::vector<vec3> pos3;
    std::map<std::tuple<long long, long long, long long>, std::vector<int>>
        vidGrid;
    auto vidCellOf = [&](const vec3& p) {
      return std::make_tuple(static_cast<long long>(std::floor(p.x / rho)),
                             static_cast<long long>(std::floor(p.y / rho)),
                             static_cast<long long>(std::floor(p.z / rho)));
    };
    auto vid = [&](const vec3& p) -> int {
      auto it = vidOf.find(KeyOf(p));
      if (it != vidOf.end()) return it->second;
      const auto [cx, cy, cz] = vidCellOf(p);
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz) {
            const auto git = vidGrid.find({cx + dx, cy + dy, cz + dz});
            if (git == vidGrid.end()) continue;
            for (const int w : git->second)
              if (la::length(pos3[w] - p) <= rho) {
                vidOf.emplace(KeyOf(p), w);  // alias into the cluster
                return w;
              }
          }
      const int id = static_cast<int>(pos3.size());
      vidOf.emplace(KeyOf(p), id);
      pos3.push_back(p);
      vidGrid[{cx, cy, cz}].push_back(id);
      return id;
    };
    std::vector<std::set<int>> adj;
    auto link = [&](int a, int b) {
      if (a == b) return;
      const int mx = std::max(a, b);
      if (static_cast<int>(adj.size()) <= mx) adj.resize(mx + 1);
      adj[a].insert(b);
      adj[b].insert(a);
    };
    auto flEdgeKey = [](int a, int b) {
      return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
    };
    for (int i = 0; i < nS; ++i) {
      std::sort(splits[i].begin(), splits[i].end(),
                [](const auto& x, const auto& y) { return x.first < y.first; });
      int prev = vid(segs[i].p0);
      for (const auto& pr : splits[i]) {
        const int v = vid(pr.second);
        link(prev, v);
        if (prev != v) flEdgeSegs[flEdgeKey(prev, v)].push_back(i);
        prev = v;
      }
      const int last = vid(segs[i].p1);
      link(prev, last);
      if (prev != last) flEdgeSegs[flEdgeKey(prev, last)].push_back(i);
    }
    adj.resize(pos3.size());
    {  // v*'s walk-graph id (cluster rep), if v* is a vertex of this group
      const auto vit = vidOf.find(flStarKey);
      flVstarLocal = vit != vidOf.end() ? vit->second : -1;
    }
    // exact CCW angular order around each vertex (half-plane + orient sign)
    auto angLess = [&](int v, int a, int b) -> bool {
      const vec2 pv = e1::Drop2(pos3[v], axis);
      const vec2 pa = e1::Drop2(pos3[a], axis), pb = e1::Drop2(pos3[b], axis);
      const double ax = pa.x - pv.x, ay = pa.y - pv.y;
      const double bx = pb.x - pv.x, by = pb.y - pv.y;
      const int ha = (ay > 0 || (ay == 0 && ax > 0)) ? 0 : 1;
      const int hb = (by > 0 || (by == 0 && bx > 0)) ? 0 : 1;
      if (ha != hb) return ha < hb;
      return e1::O2(pos3[v], pos3[a], pos3[b], axis) > 0;
    };
    std::vector<std::map<int, int>> cwprev(pos3.size());
    for (size_t v = 0; v < pos3.size(); ++v) {
      std::vector<int> nb(adj[v].begin(), adj[v].end());
      std::sort(nb.begin(), nb.end(), [&](int a, int b) {
        return angLess(static_cast<int>(v), a, b);
      });
      for (size_t k = 0; k < nb.size(); ++k)
        cwprev[v][nb[k]] = nb[(k + nb.size() - 1) % nb.size()];
    }
    // WALK TRACER (E1_WALKAT="g,u1,v1,u2,v2"): follow the rotation walk
    // from the half-edge whose endpoints are nearest the two 2D points -
    // the face boundary containing that half-edge, step by step.
    if (const char* wa = std::getenv("E1_WALKAT")) {
      int pg = -1;
      double u1, v1, u2, v2;
      if (std::sscanf(wa, "%d,%lf,%lf,%lf,%lf", &pg, &u1, &v1, &u2, &v2) == 5 &&
          pg == g) {
        auto nearest = [&](double uu, double vv) -> int {
          int best = -1;
          double bd = 1e300;
          for (size_t k = 0; k < pos3.size(); ++k) {
            const vec2 p = e1::Drop2(pos3[k], axis);
            const double d = la::length(p - vec2{uu, vv});
            if (d < bd) {
              bd = d;
              best = static_cast<int>(k);
            }
          }
          return best;
        };
        const int va = nearest(u1, v1), vb = nearest(u2, v2);
        std::fprintf(stderr, "E1 WALK start %d->%d\n", va, vb);
        int ca = va, cb = vb;
        for (int s = 0; s < 60; ++s) {
          const vec2 p = e1::Drop2(pos3[ca], axis);
          std::fprintf(stderr, "  step %d: v%d (%.12g,%.12g) deg=%d\n", s, ca,
                       p.x, p.y, static_cast<int>(adj[ca].size()));
          const auto it = cwprev[cb].find(ca);
          if (it == cwprev[cb].end()) {
            std::fprintf(stderr, "  walk broke\n");
            break;
          }
          ca = cb;
          cb = it->second;
          if (ca == va && cb == vb) {
            std::fprintf(stderr, "  walk closed after %d steps\n", s + 1);
            break;
          }
        }
      }
    }
    std::set<std::pair<int, int>> used;
    std::vector<std::vector<int>> cells, negloops;
    for (size_t a0 = 0; a0 < pos3.size(); ++a0)
      for (const int b0 : adj[a0]) {
        if (used.count({static_cast<int>(a0), b0})) continue;
        std::vector<int> loop;
        int ca = static_cast<int>(a0), cb = b0;
        bool ok = true;
        for (int guard = 0; guard < 4 * static_cast<int>(pos3.size()) + 16;
             ++guard) {
          used.insert({ca, cb});
          loop.push_back(ca);
          const auto it = cwprev[cb].find(ca);
          if (it == cwprev[cb].end()) {
            ok = false;
            break;
          }
          ca = cb;
          cb = it->second;
          if (ca == static_cast<int>(a0) && cb == b0) break;
          if (guard == 4 * static_cast<int>(pos3.size()) + 15) ok = false;
        }
        if (!ok || loop.size() < 3) continue;
        // BRIDGE + SPUR EXCISION: an undirected edge the walk traverses in
        // BOTH directions separates nothing (the dangling-chord family - a
        // seam with the same face on both sides); the walk is then two lobes
        // joined by a zero-width corridor, which no polygon triangulation can
        // cover soundly (measured: bridge-spanning ears/diagonals emit
        // geometry OUTSIDE the face).  Split the walk at every doubled edge
        // into its lobes (dropping the bridge), then excise 2-step spur tips.
        std::vector<std::vector<int>> work{loop};
        while (!work.empty()) {
          std::vector<int> L = std::move(work.back());
          work.pop_back();
          if (L.size() < 3) continue;
          const int m = static_cast<int>(L.size());
          int k1 = -1, k2 = -1;
          std::map<std::pair<int, int>, int> firstAt;
          for (int k = 0; k < m && k1 < 0; ++k) {
            const int u = L[k], v = L[(k + 1) % m];
            const auto it = firstAt.find({v, u});
            if (it != firstAt.end()) {
              k1 = it->second;
              k2 = k;
            } else {
              firstAt.emplace(std::make_pair(u, v), k);
            }
          }
          if (k1 >= 0) {  // doubled edge e_k1=(u,v), e_k2=(v,u): split lobes
            std::vector<int> l1, l2;
            for (int k = k1 + 1; k < k2; ++k) l1.push_back(L[k]);
            for (int k = (k2 + 1) % m; k != k1; k = (k + 1) % m) {
              l2.push_back(L[k]);
              if (static_cast<int>(l2.size()) > m) break;  // guard
            }
            work.push_back(std::move(l1));
            work.push_back(std::move(l2));
            continue;
          }
          // PINCH DECOMPOSITION: a vertex visited twice (no doubled edge)
          // joins lobes at a point; split there.  Both-positive lobes are
          // separate cells; a negative lobe is a hole ring and is ROUTED
          // through containment + keyhole attachment below (the offline v3
          // island lesson: never dropped).
          // ... but only when both lobes have the SAME area sign (bowtie /
          // interleaved-rotation artifacts).  An OPPOSITE-sign pinch is a
          // hole ring touching its cell at a vertex: keep the walk WHOLE -
          // it is already the keyhole form the triangulators handle, and
          // decomposing forces containment + splice round-trips that create
          // multi-keyholes with shared ring vertices (measured stall).
          {
            auto lobeSign = [&](const std::vector<int>& lb) -> int {
              const double s = e1::LoopShoelace(lb, pos3, axis);
              return s > 0 ? 1 : (s < 0 ? -1 : 0);
            };
            std::map<int, int> seen;
            bool split = false;
            for (int k = 0; k < m && !split; ++k) {
              const auto it = seen.find(L[k]);
              if (it == seen.end()) {
                seen.emplace(L[k], k);
                continue;
              }
              const int pi = it->second, pk = k;
              std::vector<int> l1(L.begin() + pi, L.begin() + pk);
              std::vector<int> l2(L.begin(), L.begin() + pi);
              l2.insert(l2.end(), L.begin() + pk, L.end());
              if (l1.size() < 3 || l2.size() < 3) {
                work.push_back(std::move(l1));
                work.push_back(std::move(l2));
                split = true;
                break;
              }
              const int s1 = lobeSign(l1), s2 = lobeSign(l2);
              if (s1 == s2 || s1 == 0 || s2 == 0) {
                work.push_back(std::move(l1));
                work.push_back(std::move(l2));
                split = true;
              }
            }
            if (split) continue;
          }
          // spur tips
          bool changed = true;
          while (changed && L.size() >= 3) {
            changed = false;
            const int mm = static_cast<int>(L.size());
            for (int k = 0; k < mm; ++k)
              if (L[(k + mm - 1) % mm] == L[(k + 1) % mm]) {
                const int hi = std::max((k + 1) % mm, k);
                const int lo = std::min((k + 1) % mm, k);
                L.erase(L.begin() + hi);
                L.erase(L.begin() + lo);
                changed = true;
                break;
              }
          }
          if (L.size() < 3) continue;
          const double s = e1::LoopShoelace(L, pos3, axis);
          if (s > 0)
            cells.push_back(std::move(L));
          else if (s < 0)
            negloops.push_back(std::move(L));
        }
      }

    // ---- 5. disconnected island hole-rings: containment + keyhole ----
    auto loopArea = [&](const std::vector<int>& loop) -> double {
      return e1::LoopShoelace(loop, pos3, axis);
    };
    auto inLoop = [&](const vec2& p, const std::vector<int>& loop) -> bool {
      static const double kSlope[] = {1.0 / 7919.0, 3.0 / 104729.0,
                                      -5.0 / 1299709.0, 7.0 / 15485863.0};
      for (const double r : kSlope) {
        int cnt = 0;
        bool ok = true;
        for (size_t k = 0; k < loop.size() && ok; ++k) {
          const vec2 A2 = e1::Drop2(pos3[loop[k]], axis);
          const vec2 B2 = e1::Drop2(pos3[loop[(k + 1) % loop.size()]], axis);
          const double dx = B2.x - A2.x, dy = B2.y - A2.y;
          const double det = dy - r * dx;
          if (det == 0.0) continue;
          const double u = (-(r) * (p.x - A2.x) + (p.y - A2.y)) / det;
          const double t = (dx * (p.y - A2.y) - dy * (p.x - A2.x)) / det;
          if (u == 0.0 || u == 1.0 || t == 0.0)
            ok = false;
          else if (u > 0 && u < 1 && t > 0)
            ++cnt;
        }
        if (ok) return (cnt % 2) == 1;
      }
      return false;
    };
    if (!negloops.empty()) {
      std::map<int, std::vector<std::vector<int>>> holeof;
      for (const auto& nl : negloops) {
        // try several ring vertices: the parity test degenerates when the
        // probe vertex sits ON the containing cell's boundary (the pinch
        // vertex of a decomposed island ring) - a single-vertex test then
        // orphans the ring and its edges vanish (measured: micro-corner
        // island opens)
        int best = -1;
        double bestA = 0.0;
        for (size_t pv = 0; pv < nl.size() && pv < 4 && best < 0; ++pv) {
          const vec2 p = e1::Drop2(pos3[nl[pv]], axis);
          for (size_t ci = 0; ci < cells.size(); ++ci)
            if (inLoop(p, cells[ci])) {
              const double a = std::abs(loopArea(cells[ci]));
              if (best < 0 || a < bestA) {
                best = static_cast<int>(ci);
                bestA = a;
              }
            }
        }
        // DUST HOLE RING: a ring whose width is below the weld radius is a
        // sub-representable near-duplicate zigzag (measured: macro-long,
        // ~1e-12-wide rings in the near-tangent fold overlap); splicing it
        // strangles the triangulation on sub-ULP structure, and at the weld
        // it vanishes anyway - drop the ring, keep the cell solid.
        {
          const double s = e1::LoopShoelace(nl, pos3, axis);
          double ext = 0.0;
          vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
          for (size_t k = 0; k < nl.size(); ++k) {
            const vec2 p1 = e1::Drop2(pos3[nl[k]], axis);
            lo = la::min(lo, p1);
            hi = la::max(hi, p1);
          }
          ext = std::max(hi.x - lo.x, hi.y - lo.y);
          // extent below the weld radius is dust outright (dust-dot rings'
          // shoelace is rounding noise, the width ratio garbage)
          if (ext <= 0.99 * eps || std::abs(0.5 * s) / ext <= 0.99 * eps)
            continue;
        }
        if (best >= 0) holeof[best].push_back(nl);
        // best<0: the component's outer contour - contained in nothing.
      }
      for (auto& kv : holeof) {
        std::vector<int>& merged = cells[kv.first];
        std::vector<std::vector<int>>& pend = kv.second;
        while (!pend.empty()) {
          bool spliced = false;
          for (size_t hi = 0; hi < pend.size() && !spliced; ++hi) {
            const std::vector<int>& hole = pend[hi];
            for (size_t i = 0; i < merged.size() && !spliced; ++i) {
              const vec3& b = pos3[merged[i]];
              for (size_t j = 0; j < hole.size() && !spliced; ++j) {
                const vec3& d = pos3[hole[j]];
                if (KeyOf(b) == KeyOf(d)) continue;
                bool ok = true;
                auto checkRing = [&](const std::vector<int>& ring) {
                  const size_t m = ring.size();
                  for (size_t k = 0; k < m && ok; ++k) {
                    const vec3& P = pos3[ring[k]];
                    const vec3& Q = pos3[ring[(k + 1) % m]];
                    if (KeyOf(P) == KeyOf(b) || KeyOf(P) == KeyOf(d) ||
                        KeyOf(Q) == KeyOf(b) || KeyOf(Q) == KeyOf(d))
                      continue;
                    if (e1::ProperCross2(b, d, P, Q, axis)) ok = false;
                  }
                  for (size_t k = 0; k < m && ok; ++k) {
                    const vec3& v = pos3[ring[k]];
                    if (KeyOf(v) == KeyOf(b) || KeyOf(v) == KeyOf(d)) continue;
                    if (e1::OnOpenSeg2(b, d, v, axis)) ok = false;
                  }
                };
                checkRing(merged);
                for (const auto& h2 : pend)
                  if (ok) checkRing(h2);
                if (!ok) continue;
                std::vector<int> nm(merged.begin(), merged.begin() + i + 1);
                for (size_t k = 0; k <= hole.size(); ++k)
                  nm.push_back(hole[(j + k) % hole.size()]);
                nm.insert(nm.end(), merged.begin() + i, merged.end());
                merged = nm;
                pend.erase(pend.begin() + hi);
                spliced = true;
              }
            }
          }
          if (!spliced) {
            ++spliceFail;
            break;
          }
        }
      }
    }

    tWalk += flMs(flTW0, flNow());
    const auto flTC0 = flNow();
    // ---- 6. classify + emit ----
    nCells += static_cast<int>(cells.size());
    nNeg += static_cast<int>(negloops.size());
    // window dump (E1_CELLWIN="g,u1,v1,u2,v2"): every cell whose bbox
    // intersects the window, with its vertex list - collection-vs-classify
    // disposition tracing.
    static const char* kCellWin = std::getenv("E1_CELLWIN");
    if (kCellWin) {
      int pg = -1;
      double wu1, wv1, wu2, wv2;
      if (std::sscanf(kCellWin, "%d,%lf,%lf,%lf,%lf", &pg, &wu1, &wv1, &wu2,
                      &wv2) == 5 &&
          pg == g) {
        auto dumpLoop = [&](const std::vector<int>& L, const char* tag) {
          vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
          for (const int v : L) {
            const vec2 p = e1::Drop2(pos3[v], axis);
            lo = la::min(lo, p);
            hi = la::max(hi, p);
          }
          if (hi.x < wu1 || lo.x > wu2 || hi.y < wv1 || lo.y > wv2) return;
          std::fprintf(stderr, "E1 CELLWIN %s n=%d:", tag,
                       static_cast<int>(L.size()));
          for (size_t k = 0; k < L.size() && k < 16; ++k)
            std::fprintf(stderr, " v%d", L[k]);
          std::fprintf(stderr, "\n");
          for (size_t k = 0; k < L.size() && k < 16; ++k) {
            const vec2 p = e1::Drop2(pos3[L[k]], axis);
            std::fprintf(stderr, "    v%d (%.17g,%.17g)\n", L[k], p.x, p.y);
          }
        };
        for (const auto& L : cells) dumpLoop(L, "cell");
        for (const auto& L : negloops) dumpLoop(L, "neg");
      }
    }
    for (const std::vector<int>& loop : cells) {
      std::vector<ivec3> tris;
      if (!exacttri::Triangulate(loop, pos3, axis, eps, tris)) {
        // ALL-OR-NOTHING: a partial covering emits unpaired interior edges
        // (the offline lesson) - discard, then adjudicate by WIDTH: a loop
        // whose area/extent is below the weld scale is a rounded-degenerate
        // sliver (non-simple at double precision) that welds away; anything
        // wider is an honest triangulation failure - DEMOTED (owner
        // adjudication, s2-tri): DEBUG assert; in release SKIP the cell and
        // let the emission's re-gate refuse the resulting unpaired boundary
        // (fail-closed, with the skip count breadcrumbed into the fatal
        // detail).  NO fallback triangulator - a macro failure here is an
        // upstream invariant violation and masking it is banned.
        tris.clear();
        const double s = e1::LoopShoelace(loop, pos3, axis);
        double ext = 0.0;
        vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
        for (size_t k = 0; k < loop.size(); ++k) {
          const vec2 p1 = e1::Drop2(pos3[loop[k]], axis);
          lo = la::min(lo, p1);
          hi = la::max(hi, p1);
        }
        ext = std::max(hi.x - lo.x, hi.y - lo.y);
        const double width = ext > 0.0 ? std::abs(0.5 * s) / ext : 0.0;
        // an extent below the weld radius is dust outright (the shoelace of
        // a dust-dot loop is rounding noise and the width ratio is garbage)
        if (ext <= 0.99 * eps || width <= 0.99 * eps) {
          ++nDustCell;  // rounded-degenerate sliver cell: weld dust
          continue;
        }
        DEBUG_ASSERT(false, geometryErr,
                     "regularize3d: macro arrangement cell failed exact "
                     "triangulation");
        ++triFail;  // breadcrumb: rides into the re-gate's fatal detail
        if (kDump && triFail <= 8) {
          const vec2 c0 = e1::Drop2(pos3[loop[0]], axis);
          std::fprintf(stderr,
                       "E1 TRIFAIL g=%d n=%d width=%.3g ext=%.3g s=%.3g "
                       "at2d=(%.9g,%.9g) p0=(%.9g,%.9g,%.9g)\n",
                       g, static_cast<int>(loop.size()), width, ext, s, c0.x,
                       c0.y, pos3[loop[0]].x, pos3[loop[0]].y, pos3[loop[0]].z);
          if (triFail <= 3)
            for (size_t k = 0; k < loop.size(); ++k) {
              const vec2 p = e1::Drop2(pos3[loop[k]], axis);
              std::fprintf(stderr, "    v%zu id=%d (%.17g,%.17g)\n", k, loop[k],
                           p.x, p.y);
            }
        }
      }
      if (tris.empty()) continue;  // dust cell (collinear at double precision)
      // interior point: largest sub-triangle's centroid
      int best = 0;
      double bestA = -1.0;
      for (size_t k = 0; k < tris.size(); ++k) {
        const vec2 p0 = e1::Drop2(pos3[tris[k].x], axis),
                   p1 = e1::Drop2(pos3[tris[k].y], axis),
                   p2 = e1::Drop2(pos3[tris[k].z], axis);
        const double a2 = std::abs(la::cross(p1 - p0, p2 - p0));
        if (a2 > bestA) {
          bestA = a2;
          best = static_cast<int>(k);
        }
      }
      const vec3 cen3 =
          (pos3[tris[best].x] + pos3[tris[best].y] + pos3[tris[best].z]) / 3.0;
      // PROJECT the probe center onto the group plane FIRST: junction-split
      // vertices sit up to eps OFF it (the eps-snapped T-junction family), so
      // the raw sub-tri centroid can be further from the member sheet than
      // the probe offset - both probes landing one-sided (measured).  The
      // projection also moves the LATERAL position (nHat has in-plane
      // components), so coverage MUST be evaluated at the SAME cenP the
      // probes use (measured: an edge-adjacent cell classified covering at
      // cen3 while the probe ran just outside the member at cenP).
      const double d0 =
          la::dot(Nrep, cen3 - A.tri[rep[g]][0]) / la::length(Nrep);
      const vec3 cenP = cen3 - d0 * nHat;
      // NET covering jump: members + covering STACK sheets (owner rule: the
      // lowest covering gid emits the stack's net transition).  Containment
      // projects along the TESTED face's OWN dominant axis (a group-axis
      // projection is meaningless for a transversal face whose plane merely
      // passes near the probe - measured: a steep wall z-projected "covering"
      // a cap probe its sheet never touches).
      auto covers = [&](int f) -> bool {
        const int ax = DominantAxis(A.faceN[f]);
        const int o0 = e1::O2(A.tri[f][0], A.tri[f][1], cenP, ax);
        const int o1 = e1::O2(A.tri[f][1], A.tri[f][2], cenP, ax);
        const int o2 = e1::O2(A.tri[f][2], A.tri[f][0], cenP, ax);
        const bool neg = o0 < 0 || o1 < 0 || o2 < 0;
        const bool pos = o0 > 0 || o1 > 0 || o2 > 0;
        return !(neg && pos);
      };
      // targeted classify dump: E1_PROBEAT="g,u,v" (group + 2D point)
      static const char* kProbeAt = std::getenv("E1_PROBEAT");
      bool probeHit = false;
      if (kProbeAt) {
        int pg = -1;
        double pu = 0, pv = 0;
        if (std::sscanf(kProbeAt, "%d,%lf,%lf", &pg, &pu, &pv) == 3 &&
            pg == g) {
          // containment by crossing parity (bbox was too coarse)
          int cnt = 0;
          const int m = static_cast<int>(loop.size());
          for (int k = 0; k < m; ++k) {
            const vec2 A2 = e1::Drop2(pos3[loop[k]], axis);
            const vec2 B2 = e1::Drop2(pos3[loop[(k + 1) % m]], axis);
            const double dx = B2.x - A2.x, dy = B2.y - A2.y;
            const double r = 1.0 / 7919.0;
            const double det = dy - r * dx;
            if (det == 0.0) continue;
            const double uu = (-(r) * (pu - A2.x) + (pv - A2.y)) / det;
            const double tt = (dx * (pv - A2.y) - dy * (pu - A2.x)) / det;
            if (uu > 0 && uu < 1 && tt > 0) ++cnt;
          }
          probeHit = (cnt % 2) == 1;
        }
      }
      int jump = 0;
      bool anyOwn = false;
      for (const int f : members[g])
        if (covers(f)) {
          jump += fsgn[f];
          anyOwn = true;
        }
      if (probeHit)
        std::fprintf(stderr,
                     "E1 PROBEAT g=%d n=%d anyOwn=%d jump0=%d cen=(%.9g,%.9g,"
                     "%.9g)\n",
                     g, static_cast<int>(loop.size()), anyOwn ? 1 : 0, jump,
                     cenP.x, cenP.y, cenP.z);
      if (!anyOwn) continue;  // no member sheet here (shadow-only region)
      ++nJump;
      // PER-CELL STACK + probe offset by GAP-FINDING over the nearby foreign
      // plane distances: grow the stack threshold T until an 8x gap opens,
      // probe at 2T (above the whole stack, below a quarter of everything
      // else).  Stack sheets covering the point join the NET jump; the lowest
      // covering gid owns the cell.
      std::vector<std::pair<double, int>> nearD;
      const double m = 1e-6 * (1.0 + scale);
      flBoxQuery(cenP - vec3(m), cenP + vec3(m));
      for (const int f2 : flCands) {
        if (gid[f2] == g) continue;
        // ALONG-nHat crossing distance: the probe segment runs along nHat,
        // so the relevant quantity is where the segment meets the plane, not
        // the perpendicular plane distance (a steep transversal wall has a
        // tiny perpendicular distance near its seam line yet its crossing
        // sits far outside the probe window - measured).  A plane parallel
        // to the probe direction is never crossed: skip.
        const double denom = std::abs(la::dot(A.faceN[f2], nHat));
        if (!(denom > 0.0)) continue;
        const double dist =
            std::abs(la::dot(A.faceN[f2], cenP - A.tri[f2][0])) / denom;
        if (dist < eps * 1e5) nearD.push_back({dist, f2});
      }
      std::sort(nearD.begin(), nearD.end());
      double T = stackWin;
      double dNext = std::numeric_limits<double>::infinity();
      for (const auto& dn : nearD) {
        if (dn.first <= T) continue;
        // HARD CAP at the weld radius: a stack may only absorb sheets that
        // weld together anyway (unrepresentably close).  Un-capped 8x
        // chaining absorbed ladders of REPRESENTABLE distinct sheets into
        // one net emission (measured: systematic opens across dozens of
        // groups wherever the sheet-distance ladder had no 8x gap).
        if (dn.first <= 8.0 * T && dn.first <= 0.99 * eps) {
          T = dn.first;
        } else {
          dNext = dn.first;
          break;
        }
      }
      // The probe must sit in the VERIFIED gap: above the whole stack, below
      // the FIRST non-absorbed sheet.  A sheet in (T, 2T) exists when the
      // 8x-chain hits the 0.99*eps hard cap (a representable sheet just
      // above the weld radius); probing at a blind 2T then crosses it
      // without it being in the jump - the certificate refused honestly but
      // wrongly (measured: GT7081 shell2, jump=1 probed 0|0 at a 6mm cell).
      const double off = 2.0 * T < dNext ? 2.0 * T : 0.5 * (T + dNext);
      const int ownJump = jump;  // own-group net, before the stack merge
      bool owned = false;
      bool stackMerged = false;
      int ownerFace = -1;
      // merged stack sheets as (signed along-nHat position, crossing sign):
      // the LAYER structure between the sheets decides whether a
      // net-cancelled stack is a material slab or a void one (sheet order is
      // constant across a cell - cells are split at the mutual seams).
      std::vector<std::pair<double, int>> stackSheets;
      for (const auto& dn : nearD) {
        if (dn.first > T) break;
        if (!covers(dn.second)) continue;
        if (gid[dn.second] < g) {
          owned = true;
          ownerFace = dn.second;
          break;
        }
        const int s = la::dot(A.faceN[dn.second], nHat) >= 0.0 ? 1 : -1;
        const double tPos =
            -la::dot(A.faceN[dn.second], cenP - A.tri[dn.second][0]) /
            la::dot(A.faceN[dn.second], nHat);
        jump += s;
        stackSheets.push_back({tPos, s});
        stackMerged = true;
      }
      if (owned) {
        if (probeHit)
          std::fprintf(stderr, "E1 PROBEAT g=%d OWNED by f=%d gid=%d\n", g,
                       ownerFace, gid[ownerFace]);
        ++nOwned;
        continue;  // a lower covering group owns this stack cell
      }
      // PANCAKE ARM: a net-cancelled STACK (own sheet + sub-weld-close
      // foreign sheets summing to zero) over a cell of REPRESENTABLE width
      // is a material/void slab thinner than the weld radius but wider than
      // it: its bounding sheets collapse to one membrane at assembly, yet
      // its SIDE CHAINS (up to a few eps apart) do not weld - dropping the
      // cell leaves the neighbours' sheets ringing an unpairable hole
      // (measured: the two near-tangent corner strips, oracle transects
      // 0|1|0 across a 2.5e-11 slab).  The stack OWNER emits the collapsed
      // membrane ONCE, oriented as its OWN sheet so it continues the
      // own-plane neighbours; the certificate for an invisible slab is
      // probe agreement (wA == wB).  A same-plane cancellation
      // (ownJump == 0, no stack) stays dropped: exactly-coplanar sheets
      // cancel pointwise and the neighbouring fans pair without a cap.
      // A sub-weld-WIDE pancake also stays dropped: its side chains weld
      // together and the hole closes in the assembly.
      const bool pancake = (jump == 0 && stackMerged && ownJump != 0);
      // ---- flood node (flip arc stage 1): every classified cell joins the
      // component-global field graph, including net-cancelled conduits ----
      const int flNode = static_cast<int>(flNodes.size());
      {
        int sPosMid = 0, sZero = 0;
        bool hasBand = false;
        std::vector<std::array<int, 3>> nearSign;
        for (const auto& dn : nearD) {
          if (!covers(dn.second)) continue;
          const double den = la::dot(A.faceN[dn.second], nHat);
          const double tp =
              -la::dot(A.faceN[dn.second], cenP - A.tri[dn.second][0]) / den;
          const int s = den >= 0.0 ? 1 : -1;
          const int state = tp > stackWin ? 1 : (tp < -stackWin ? -1 : 0);
          if (state == 0) {
            hasBand = true;
            if (tp == 0.0) sZero += s;
          }
          // E' (the field) sits just above the sub-resolution band: the
          // anchor conversion crosses exactly the robust ABOVE sheets below
          // the probe (band sides are double-noise and never counted).
          if (state > 0 && dn.first < off) sPosMid += s;
          nearSign.push_back({dn.second, state, s});
        }
        std::sort(nearSign.begin(), nearSign.end());
        FlNode fn;
        fn.nearSign = std::move(nearSign);
        fn.hasBand = hasBand;
        fn.g = g;
        fn.jump = jump;
        fn.ownJump = ownJump;
        fn.sPosMid = sPosMid;
        fn.sZero = sZero;
        fn.pancake = pancake;
        fn.Nrep = Nrep;
        fn.stackSheets = stackSheets;
        flNodes.push_back(std::move(fn));
        bool touch = false;
        const int m = static_cast<int>(loop.size());
        for (int k = 0; k < m; ++k) {
          if (loop[k] == flVstarLocal) touch = true;
          const int a = loop[k], b = loop[(k + 1) % m];
          if (a == b) continue;
          flEdgeCells[a < b ? std::make_pair(a, b) : std::make_pair(b, a)]
              .push_back({flNode, a < b});
        }
        if (touch && flVstarLocal >= 0) flStarNodes.push_back(flNode);
      }
      // FLOOD-PRIMARY (flip arc session 3): the classify loop no longer
      // probes; every emission-eligible cell stores its payload and the
      // field solve values it (residual probes anchor subgraphs; a BFS
      // conflict falls back to probe-everything = the old semantics).
      bool emitEligible = !(jump == 0 && !pancake);
      if (emitEligible && pancake) {
        // Dust rule for pancakes is EXTENT-outright (the dust-dot rule),
        // NOT altitude: a strip 0.96 eps wide does NOT weld closed when its
        // side-chain vertices are staggered along the strip (measured: the
        // corner strip, chains 0.96 eps apart laterally, nearest vertices
        // 3.4 eps apart) - the vertex weld never pairs the chains and only
        // the cap can.  Only a cell welding to a single POINT is skippable.
        vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
        for (const int lv : loop) {
          const vec2 p = e1::Drop2(pos3[lv], axis);
          lo = la::min(lo, p);
          hi = la::max(hi, p);
        }
        if (std::max(hi.x - lo.x, hi.y - lo.y) <= 0.99 * eps) {
          ++nDustCell;
          emitEligible = false;  // dust dot: welds to a point
        } else {
          ++nPancake;
        }
      }
      {
        FlNode& fn = flNodes[flNode];
        fn.cenP = cenP;
        fn.off = off;
        if (emitEligible) {
          fn.emit = true;
          fn.tris = tris;
          for (const ivec3& t : tris) {
            fn.pos.emplace(t.x, pos3[t.x]);
            fn.pos.emplace(t.y, pos3[t.y]);
            fn.pos.emplace(t.z, pos3[t.z]);
          }
          double eMax = 0.0;
          {
            const vec2 q0 = e1::Drop2(pos3[tris[best].x], axis);
            const vec2 q1 = e1::Drop2(pos3[tris[best].y], axis);
            const vec2 q2 = e1::Drop2(pos3[tris[best].z], axis);
            eMax = std::max({la::length(q1 - q0), la::length(q2 - q1),
                             la::length(q0 - q2)});
          }
          fn.extW = eMax > 0.0 ? bestA / eMax : 0.0;
        }
      }
      if (!kFloodDiff) continue;
      // FULL DIFFERENTIAL (E1_FLOODDIFF): probe + certify every cell as a
      // shadow validator; the values are NOT anchors (the solve anchors from
      // residual probes + seeds, so the differential grades flood-primary).
      const std::optional<int> wA = RobustWinding(in, cenP + off * nHat, seeds);
      const std::optional<int> wB = RobustWinding(in, cenP - off * nHat, seeds);
      if (!wA || !wB) {
        flNodes[flNode].probeState = 2;
        continue;
      }
      const bool certified = pancake ? (*wA == *wB) : (*wB - *wA == jump);
      {
        FlNode& fn = flNodes[flNode];
        fn.probeState = certified ? 1 : 3;
        fn.certified = certified;
        fn.anchorOK = certified;  // base-offset probe: anchor-grade
        fn.probeWA = *wA;
        fn.probeWB = *wB;
      }
      if (probeHit)
        std::fprintf(stderr,
                     "E1 PROBEAT g=%d jump=%d wA=%d wB=%d off=%.3g cert=%d\n",
                     g, jump, *wA, *wB, off, certified ? 1 : 0);
      {
        static const char* kFE = std::getenv("E1_FLOODEDGE");
        if (kFE && std::atoi(kFE) == g) {
          std::fprintf(stderr,
                       "E1 FLNODE n=%d g=%d wA=%d wB=%d jump=%d ownJump=%d "
                       "sPos=%d cenP=(%.6g,%.6g,%.6g) nHat=(%.3g,%.3g,%.3g)\n",
                       flNode, g, *wA, *wB, jump, ownJump,
                       flNodes[flNode].sPosMid, cenP.x, cenP.y, cenP.z, nHat.x,
                       nHat.y, nHat.z);
          for (const auto& e2 : flNodes[flNode].nearSign)
            std::fprintf(stderr, "    near f=%d gid=%d tSign=%d q=%d\n", e2[0],
                         gid[e2[0]], e2[1], e2[2]);
        }
      }
    }
    tClassify += flMs(flTC0, flNow());
    // ---- flood graph (flip arc stage 1): this group's edges + handoffs ----
    // // ---- flood graph (flip arc stage 1): this group's edges + handoffs
    // ----
    {
      // Drop2 has no parity swap, so cell-CCW-in-frame is about +nHat only up
      // to the frame sign: (y,z)/( x,y) are right-handed about +x/+z, (x,z) is
      // LEFT-handed about +y; and the dominant component of Nrep signs the
      // projection direction.  t_out (from the cell traversing lo->hi toward
      // its neighbor) = sigmaF * cross(d3, Nrep).
      const int sigmaF = (axis == 1 ? -1 : 1) * (Nrep[axis] > 0.0 ? 1 : -1);
      // Exact plane side of a point vs the group's rep triangle (a,b,c):
      // Orient3DFilterSign's det is (a-d).((b-d)x(c-d)) = -dot(n_abc, d-a),
      // so the +Nrep side is a NEGATIVE sign.  Filter-first, exact on 0 -
      // the blessed WindCrossTri escalation pattern.
      auto flAbovePlane = [&](const vec3& x) -> int {  // +1 above, -1 below
        const std::array<vec3, 3>& rt = A.tri[rep[g]];
        int s = Orient3DFilterSign(rt[0], rt[1], rt[2], x);
        if (s == 0) s = Orient3DExactSign(rt[0], rt[1], rt[2], x);
        return -s;
      };
      static const char* kFloodEdge = std::getenv("E1_FLOODEDGE");
      const bool flDump = kFloodEdge && std::atoi(kFloodEdge) == g;
      for (const auto& ec : flEdgeCells) {
        const auto& cl = ec.second;
        const auto es = flEdgeSegs.find(ec.first);
        if (flDump) {
          std::fprintf(stderr, "E1 FLEDGE g=%d (%d,%d) cells=[", g,
                       ec.first.first, ec.first.second);
          for (const auto& pr : cl)
            std::fprintf(stderr, " n%d%s", pr.first, pr.second ? "+" : "-");
          std::fprintf(stderr, " ] segs=%d\n",
                       es == flEdgeSegs.end()
                           ? -1
                           : static_cast<int>(es->second.size()));
          if (es != flEdgeSegs.end())
            for (const int si : es->second)
              std::fprintf(
                  stderr, "    seg%d planeQ=%d fOwn=%d gidOwn=%d vid=(%d,%d)\n",
                  si, segs[si].planeQ, segs[si].fOwn, gid[segs[si].fOwn],
                  segs[si].vidLo, segs[si].vidHi);
        }
        if (es == flEdgeSegs.end()) continue;
        const vec3 d3 = pos3[ec.first.second] - pos3[ec.first.first];
        const vec3 tOut = static_cast<double>(sigmaF) * la::cross(d3, Nrep);
        // Scan the segs riding this sub-edge: chord partners cross the plane
        // HERE (each is one eps-layer crossing, sign -sgn(dot(t_out, n)));
        // member mesh edges whose manifold twin RISES off-plane (+Nrep side,
        // exact sign) are RIDGES - the rising twin cuts the eps-layer at the
        // edge exactly like a chord.  Its crossing sign is combinatorial:
        // the twin halfedge runs anti-parallel to the member's, and for a
        // CCW twin rising with interior direction m, n_twin = cross(d_twin,
        // m), giving delta(c_lo->hi -> other) = sigmaF * sgn(along).
        int delta = 0, tag = 0;
        bool bad = false;
        std::set<int> seenF;  // distinct chord partner faces
        std::set<int> seenG;  // their GROUPS (parity-arm exclusion)
        std::set<std::pair<int, int>> vps;  // distinct member vid pairs
        std::set<int> ridgeTf;              // twins counted via member segs
        int handSi = -1, handTf = -1;       // single-pair handoff candidate
        bool foreignChord = false;          // chord partner != the twin
        // PASS A: member mesh-edge segs -> manifold-twin ridges.  A twin
        // rising to the +nHat side cuts the eps-layer at the edge (sign
        // combinatorial: the twin halfedge runs anti-parallel to the
        // member's, so delta = sigmaF * sgn(along) - see the derivation in
        // the flood comment).  A dipping twin leaves the layer intact.
        for (const int si : es->second) {
          if (segs[si].planeQ >= 0 || segs[si].vidLo < 0 ||
              gid[segs[si].fOwn] != g)
            continue;
          const int vlo = std::min(segs[si].vidLo, segs[si].vidHi);
          const int vhi = std::max(segs[si].vidLo, segs[si].vidHi);
          if (!vps.insert({vlo, vhi}).second) continue;  // twin's own record
          const int fOwn = segs[si].fOwn;
          int he = -1;  // recover the halfedge carrying this edge
          for (int e = 0; e < 3; ++e)
            if (A.vid[fOwn][e] == segs[si].vidLo &&
                A.vid[fOwn][(e + 1) % 3] == segs[si].vidHi)
              he = 3 * fOwn + e;
          if (he < 0) continue;
          const int tf = in.halfedge_.Pair(he) / 3;
          if (tf < 0 || gid[tf] == g) continue;  // coplanar neighbor: delta 0
          int tv = -1;  // twin's third vertex (not on the shared edge)
          for (int k = 0; k < 3; ++k)
            if (A.vid[tf][k] != segs[si].vidLo &&
                A.vid[tf][k] != segs[si].vidHi)
              tv = k;
          const double along = la::dot(d3, segs[si].p1 - segs[si].p0);
          if (tv < 0 || !(along != 0.0)) {
            bad = true;
            break;
          }
          const int s3 = flAbovePlane(A.tri[tf][tv]);
          if (s3 == 0) {
            bad = true;  // twin's third vertex exactly on-plane: degenerate
            break;
          }
          if (flDump)
            std::fprintf(stderr,
                         "    ridge fOwn=%d tf=%d gidT=%d s3=%d al=%s\n", fOwn,
                         tf, gid[tf], s3, along > 0 ? "+" : "-");
          // band check: a twin rising but staying inside the sub-resolution
          // band never crosses the E' level (h3 in doubles; s3 exact)
          const double h3 =
              std::abs(la::dot(Nrep, A.tri[tf][tv] - A.tri[rep[g]][0])) /
              la::length(Nrep);
          if (s3 > 0 && h3 > stackWin) {
            delta += sigmaF * (along > 0.0 ? 1 : -1);
            tag |= 4;
          }
          ridgeTf.insert(tf);
          seenG.insert(gid[tf]);
          handSi = si;
          handTf = tf;
        }
        // PASS B: chord partners.  The ungated seam enumeration records
        // touching contacts (edge-attached partners) as chords too, so the
        // partner's EXACT vertex plane-side census decides the topology:
        // straddle = proper transversal crossing; touch-from-above = an
        // attached riser (T-junction class, not the halfedge twin); touch-
        // from-below leaves the eps-layer intact.
        if (!bad)
          for (const int si : es->second) {
            if (segs[si].planeQ < 0) continue;
            const int f2 = segs[si].fOther;
            if (ridgeTf.count(f2)) continue;  // already counted as the twin
            if (!seenF.insert(f2).second) continue;
            seenG.insert(gid[f2]);
            int nAb = 0, nBe = 0, vAb = -1;
            std::array<int, 3> side;
            for (int k = 0; k < 3; ++k) {
              side[k] = flAbovePlane(A.tri[f2][k]);
              if (side[k] > 0) {
                ++nAb;
                vAb = k;
              } else if (side[k] < 0) {
                ++nBe;
              }
            }
            if (flDump)
              std::fprintf(stderr, "    chord fOther=%d gidO=%d ab=%d be=%d\n",
                           f2, gid[f2], nAb, nBe);
            foreignChord = true;     // a non-twin partner reaches this sub-edge
            if (nAb == 0) continue;  // touches/dips below: layer intact
            // band check (as the ridge arm): the highest above-vertex must
            // clear the band for the sheet to cross the E' level
            double hMax = 0.0;
            for (int k = 0; k < 3; ++k)
              hMax = std::max(
                  hMax,
                  std::abs(la::dot(Nrep, A.tri[f2][k] - A.tri[rep[g]][0])) /
                      la::length(Nrep));
            if (hMax <= stackWin) continue;  // in-band sheet: E' uncut
            if (nBe == 0) {
              // TOUCH-FROM-ABOVE: attached riser along its on-plane edge
              // (coincident-position attachment, not the halfedge twin).
              int e0 = -1;
              for (int k = 0; k < 3; ++k)
                if (side[k] == 0 && side[(k + 1) % 3] == 0) e0 = k;
              if (e0 < 0) continue;  // vertex-touch only: measure-zero
              const double along =
                  la::dot(d3, A.tri[f2][(e0 + 1) % 3] - A.tri[f2][e0]);
              if (!(along != 0.0)) {
                bad = true;
                break;
              }
              delta += -sigmaF * (along > 0.0 ? 1 : -1);
              tag |= 2;
              continue;
            }
            // STRADDLE: robust two-branch crossing sign.  For steep partners
            // -sgn(dot(tOut, n)) is solid; for near-tangent partners that
            // dot is cancellation noise, but q = sgn(dot(Nrep, n)) is solid
            // and the rising side r comes from the EXACT in-frame O2 of the
            // above-vertex vs the sub-edge (delta = q * r).
            const double tv2 = la::dot(tOut, A.faceN[f2]);
            const double qv = la::dot(Nrep, A.faceN[f2]);
            if (std::abs(tv2) * la::length(Nrep) >=
                std::abs(qv) * la::length(tOut)) {
              if (!(tv2 != 0.0)) {
                bad = true;
                break;
              }
              delta += tv2 > 0.0 ? -1 : 1;
              tag |= 1;
            } else {
              const int o = e1::O2(pos3[ec.first.first], pos3[ec.first.second],
                                   A.tri[f2][vAb], axis);
              if (o == 0 || !(qv != 0.0)) {
                bad = true;
                break;
              }
              delta += (qv > 0.0 ? 1 : -1) * (o < 0 ? 1 : -1);
              tag |= 1;
            }
          }
        if (flDump)
          std::fprintf(stderr, "    => delta=%d bad=%d foreign=%d vps=%d\n",
                       delta, bad ? 1 : 0, foreignChord ? 1 : 0,
                       static_cast<int>(vps.size()));
        if (bad) continue;
        // INTRA-GROUP EDGE: exactly two distinct node cells, opposite
        // traversal directions.
        if (cl.size() == 2 && cl[0].first != cl[1].first &&
            cl[0].second != cl[1].second) {
          const int nA = cl[0].second ? cl[0].first : cl[1].first;
          const int nB = cl[0].second ? cl[1].first : cl[0].first;
          // STACK-CROSSING PARITY ARM: a covering foreign SHEET whose signed
          // tPos FLIPS between the two cenPs crossed the group plane between
          // them with no representable chord (the sub-eps stack class).  The
          // eps-layer path crosses it an odd number of times (parity =
          // endpoint sign flip, path-independent).  Sheets are identified by
          // their plane GROUP (one sheet may cover the two cenPs with
          // different member triangles - measured on the near-coplanar
          // twins), aggregating q = sgn(dot(n, nHat)) into above/below sums
          // per group; under coverage conservation (same net total on both
          // sides) the contribution is aboveB - aboveA (= q*r per flipped
          // sheet, r = +1 iff above on nB's side).  Faces already counted as
          // chords/ridges on this sub-edge are excluded (a represented
          // crossing lies on every shared sub-edge of the pair).  A group
          // with an exact-zero tPos or non-conserved coverage contributes 0
          // (ends laterally; the differential owns the residue).
          {
            // {aboveA, restA, aboveB, restB} - GLOBAL aggregate (a warped
            // foreign quad's triangles land in different plane groups yet
            // form ONE crossing sheet; the global sums still conserve).
            // Only ABOVE-state (tPos > stackWin, robust) sheets sit above
            // the E' level; MID (band) and BELOW merge - a band sheet's
            // side is double-noise and its transit is not a crossing of E'
            // (the sub-ULP openscad class, measured).
            std::array<int, 4> a{0, 0, 0, 0};
            auto addTo = [&](const FlNode& fn2, int ia0, int ib0) {
              for (const auto& e2 : fn2.nearSign) {
                if (seenG.count(gid[e2[0]])) continue;
                a[e2[1] > 0 ? ia0 : ib0] += e2[2];
              }
            };
            addTo(flNodes[nA], 0, 1);
            addTo(flNodes[nB], 2, 3);
            if (a[0] + a[1] == a[2] + a[3] && a[2] != a[0]) {
              delta += a[2] - a[0];
              tag |= 8;
              if (flDump)
                std::fprintf(stderr, "    parity dAbove=%d (nA=%d nB=%d)\n",
                             a[2] - a[0], nA, nB);
            }
          }
          flEdges.push_back({nA, nB, delta, tag});
        }
        // CROSS-GROUP HANDOFF RECORD: a clean manifold mesh-edge sub-edge -
        // exactly ONE member vid pair, no THIRD sheet reaching the edge (the
        // halfedge twin's own touching-contact chord is benign; any other
        // partner is not) - whose twin face lives in another group.  Both
        // groups record E(node) + eOff; equal at the sub-edge by the
        // mesh-edge dihedral rule (eps-layer form).  Multi-pair (geometric
        // T-junction) sub-edges are skipped: the halfedge twin need not be
        // fan-adjacent there.
        if (foreignChord || vps.size() != 1 || handSi < 0) continue;
        const int fOwn = segs[handSi].fOwn;
        const double along = la::dot(d3, segs[handSi].p1 - segs[handSi].p0);
        // f's interior-side cell traverses lo->hi iff the sub-edge runs along
        // f's own winding direction XNOR f's winding is CCW-in-frame.
        const bool wantFwd = (along > 0.0) == (fsgn[fOwn] * sigmaF > 0);
        int node = -1;
        for (const auto& pr : cl)
          if (pr.second == wantFwd) {
            node = node < 0 ? pr.first : -2;  // ambiguous: drop
            if (node == -2) break;
          }
        if (node < 0) continue;
        const FlNode& fn = flNodes[node];
        // Band contamination: a covering sheet inside the sub-resolution
        // band makes the edge-level transport double-ambiguous (its side is
        // noise); the record is dropped and connectivity routes elsewhere.
        if (flNodes[node].hasBand) continue;
        // COMPOSITE-BAND GUARD: group members other than fOwn covering
        // either side of the sub-edge put the whole coplanar band into the
        // edge fan (the doubled-cap class) and the two-sheet dihedral rule
        // does not apply.  In the exactly-coplanar case the WALL side is
        // refused by its foreignChord (the band's touching-contact chord
        // rides the same sub-edge), but a sub-eps TILT moves that chord off
        // the edge and evades it (measured: the near-coplanar fold's 4 bad
        // handoffs, d=0 vs true +-1).  Skip conservatively; connectivity
        // falls to the residual probes.
        {
          int farNode = -1;
          for (const auto& pr : cl)
            if (pr.first != node && pr.second != wantFwd) farNode = pr.first;
          if (fn.ownJump != fsgn[fOwn] ||
              (farNode >= 0 && flNodes[farNode].ownJump != 0))
            continue;
        }
        // HANDOFF PARITY CORRECTION: the dihedral rule equates the eps-layer
        // values AT the edge, but the node's E is cenP-anchored; a covering
        // foreign sheet that flips plane-side between cenP and the sub-edge
        // (the stack-crossing class) shifts the edge value by q per flip.
        // Same aggregation + conservation guard as the intra parity arm,
        // with "B" = the sub-edge midpoint (sheets ending laterally break
        // conservation and skip - the differential owns that residue).  The
        // below-layer path crosses the same transiting sheets, so the
        // correction applies uniformly for both fsgn cases.
        int corr = 0;
        if (!fn.nearSign.empty()) {
          const vec3 mid = 0.5 * (pos3[ec.first.first] + pos3[ec.first.second]);
          int aCen = 0, totCen = 0, aMid = 0, totMid = 0;
          bool degen = false;
          for (const auto& e2 : fn.nearSign) {
            const int f2 = e2[0];
            // the pair's own twin (and any ridge twin attached AT this
            // sub-edge) is the dihedral rule's sector structure, not a
            // transiting foreign sheet - excluded from the leg correction
            if (f2 == handTf || ridgeTf.count(f2)) continue;
            totCen += e2[2];
            if (e2[1] > 0) aCen += e2[2];
            const int ax2 = DominantAxis(A.faceN[f2]);
            const int o0 = e1::O2(A.tri[f2][0], A.tri[f2][1], mid, ax2);
            const int o1 = e1::O2(A.tri[f2][1], A.tri[f2][2], mid, ax2);
            const int o2 = e1::O2(A.tri[f2][2], A.tri[f2][0], mid, ax2);
            const bool neg = o0 < 0 || o1 < 0 || o2 < 0;
            const bool pos = o0 > 0 || o1 > 0 || o2 > 0;
            if (neg && pos) continue;  // not covering the midpoint
            const double den = la::dot(A.faceN[f2], nHat);
            if (!(den != 0.0)) {
              degen = true;
              break;
            }
            const double tp = -la::dot(A.faceN[f2], mid - A.tri[f2][0]) / den;
            if (std::abs(tp) <= stackWin) {
              degen = true;  // transiting the band at the edge: ambiguous
              break;
            }
            totMid += e2[2];
            if (tp > stackWin) aMid += e2[2];
          }
          if (!degen && totCen == totMid) corr = aMid - aCen;
        }
        const int eOff = (fsgn[fOwn] > 0 ? 0 : fn.ownJump + fn.sZero) + corr;
        const int vlo = std::min(segs[handSi].vidLo, segs[handSi].vidHi);
        const int vhi = std::max(segs[handSi].vidLo, segs[handSi].vidHi);
        flHand[{{vlo, vhi},
                segKeyOf(pos3[ec.first.first], pos3[ec.first.second])}]
            .push_back({node, eOff});
      }
      // THE HULL SEED: exactly one node cell touching v* in this group
      // anchors combinatorially (see the field comment above).  GUARD: a
      // covering foreign sheet at the corner cell (nearSign non-empty) can
      // cross the group plane INSIDE the cell without a chord (the sub-eps
      // stack class - measured on the GT7863 twins), invalidating the
      // corner-to-cenP transport of the lexsign value; the seed is skipped
      // there (certified anchors + handoffs carry those carriers).
      if (flStarNodes.size() == 1) {
        const FlNode& fn = flNodes[flStarNodes[0]];
        const int sg = flLexsign(Nrep);
        if (sg != 0 && fn.nearSign.empty())
          flSeeds.push_back(
              {flStarNodes[0], sg > 0 ? 0 : -(fn.ownJump + fn.sZero)});
      }
    }
    if (kDump && gProbeFail > 0)
      std::fprintf(stderr, "E1 PROBEGRP g=%d fail=%d cert=%d cells=%d\n", g,
                   gProbeFail, gProbeCert, static_cast<int>(cells.size()));
    suspectTotal += suspectSnaps;
  }
  // ---- FLOOD SOLVE + EMISSION (flood-primary, flip arc session 3) --------
  // Reliability-filtered edges; anchors = residual probes (one certified
  // probe per unreached subgraph, in node order) then the guarded hull
  // seeds; a BFS conflict clears the field so the emission below probes
  // EVERY cell (the pre-flood per-cell semantics, fail-closed intact).
  // E1_FLOODDIFF probes every cell in the classify loop as a shadow
  // validator and prints the arm census; it never affects anchoring.
  const auto flTS0 = flNow();
  if (!buildPhase) {
    static const bool kFlHandDump = std::getenv("E1_FLOODHAND") != nullptr;
    int flHandEdges = 0, flHandOrphan = 0;
    for (const auto& kv : flHand) {
      if (kv.second.size() != 2) {
        if (kv.second.size() > 2) ++flHandOrphan;
        continue;
      }
      // E(n0) + eOff0 == E(n1) + eOff1  =>  E(n1) = E(n0) + (eOff0 - eOff1)
      flEdges.push_back({kv.second[0].first, kv.second[1].first,
                         kv.second[0].second - kv.second[1].second, 32});
      ++flHandEdges;
      if (kFlHandDump)
        std::fprintf(stderr,
                     "E1 FLHAND vid=(%d,%d) n%d(g%d,eOff%d) <-> n%d(g%d,"
                     "eOff%d)\n",
                     kv.first.first.first, kv.first.first.second,
                     kv.second[0].first, flNodes[kv.second[0].first].g,
                     kv.second[0].second, kv.second[1].first,
                     flNodes[kv.second[1].first].g, kv.second[1].second);
    }
    if (kFlHandDump)
      for (const auto& s : flSeeds)
        std::fprintf(stderr, "E1 FLSEED n=%d g=%d E=%d\n", s.first,
                     flNodes[s.first].g, s.second);
    const int nN = static_cast<int>(flNodes.size());
    // EDGE RELIABILITY (the near-tangent-ladder adjudication, measured on
    // openscad): (a) edges between one cell pair must AGREE - the ladder's
    // per-sub-edge chord counts are each incomplete in a different way and
    // disagree, violating path-independence; (b) an edge where BOTH a chord
    // and the parity arm fired double-counts sheets whose chords ride other
    // sub-edges of the same corridor (measured 335/335 wrong).  Unreliable
    // edges are dropped; their cells fall to residual probes.
    std::vector<bool> flOk(flEdges.size(), true);
    int flPairDrop = 0, flTagDrop = 0;
    {
      std::map<std::pair<int, int>, std::vector<size_t>> byPair;
      for (size_t i = 0; i < flEdges.size(); ++i) {
        const auto& e = flEdges[i];
        if ((e[3] & 8) && (e[3] & 3)) {
          flOk[i] = false;
          ++flTagDrop;
          continue;
        }
        byPair[e[0] < e[1] ? std::make_pair(e[0], e[1])
                           : std::make_pair(e[1], e[0])]
            .push_back(i);
      }
      for (const auto& kv : byPair) {
        bool same = true;
        int d0 = 0;
        for (size_t k = 0; k < kv.second.size(); ++k) {
          const auto& e = flEdges[kv.second[k]];
          const int d = e[0] == kv.first.first ? e[2] : -e[2];
          if (k == 0)
            d0 = d;
          else if (d != d0)
            same = false;
        }
        if (!same) {
          for (const size_t i : kv.second) flOk[i] = false;
          ++flPairDrop;
        }
      }
    }
    std::vector<std::vector<std::pair<int, int>>> nadj(nN);
    for (size_t i = 0; i < flEdges.size(); ++i) {
      if (!flOk[i]) continue;
      const auto& e = flEdges[i];
      nadj[e[0]].push_back({e[1], e[2]});
      nadj[e[1]].push_back({e[0], -e[2]});
    }
    // ARM-LEVEL DIFFERENTIAL CENSUS: every edge whose BOTH endpoints carry a
    // certified probe is directly checkable (certE(n1)-certE(n0) vs d); the
    // per-tag failure counts attribute wrong deltas to their producing arm.
    if (kFloodDiff && kDump) {
      std::map<int, std::pair<int, int>> tagCensus;  // tag -> (checked, bad)
      int badEdges = 0;
      for (const auto& e : flEdges) {
        const FlNode& a = flNodes[e[0]];
        const FlNode& b = flNodes[e[1]];
        if (!a.certified || !b.certified) continue;
        const int ca = a.probeWA + a.sPosMid, cb = b.probeWA + b.sPosMid;
        auto& tc = tagCensus[e[3]];
        ++tc.first;
        if (cb - ca != e[2]) {
          ++tc.second;
          ++badEdges;
          if (badEdges <= 12)
            std::fprintf(stderr,
                         "E1 FLOOD BADEDGE n%d(g%d)->n%d(g%d) d=%d true=%d "
                         "tag=%d\n",
                         e[0], a.g, e[1], b.g, e[2], cb - ca, e[3]);
        }
      }
      for (const auto& kv : tagCensus)
        std::fprintf(stderr, "E1 FLOOD ARMCENSUS tag=%d checked=%d bad=%d\n",
                     kv.first, kv.second.first, kv.second.second);
    }
    constexpr int kFlUnset = std::numeric_limits<int>::min();
    std::vector<int> E(nN, kFlUnset);
    std::vector<int> q;
    q.reserve(nN);
    int flMismatch = 0;
    auto flAnchor = [&](int n, int v) {
      if (E[n] == kFlUnset) {
        E[n] = v;
        q.push_back(n);
      } else if (E[n] != v) {
        ++flMismatch;
        if (kDump && flMismatch <= 8)
          std::fprintf(stderr,
                       "E1 FLOOD MISMATCH anchor n=%d g=%d E=%d vs %d\n", n,
                       flNodes[n].g, E[n], v);
      }
    };
    size_t qh = 0;
    auto flDrain = [&]() {
      for (; qh < q.size(); ++qh) {
        const int n = q[qh];
        for (const auto& [m, d] : nadj[n]) {
          const int v = E[n] + d;
          if (E[m] == kFlUnset) {
            E[m] = v;
            q.push_back(m);
          } else if (E[m] != v) {
            ++flMismatch;
            if (kDump && flMismatch <= 8)
              std::fprintf(stderr,
                           "E1 FLOOD MISMATCH n=%d g=%d E=%d vs %d (from n=%d "
                           "g=%d cert=%d)\n",
                           m, flNodes[m].g, E[m], v, n, flNodes[n].g,
                           flNodes[n].certified ? 1 : 0);
          }
        }
      }
    };
    int flResProbe = 0, flResAnchor = 0;
    auto flProbe = [&](int pn) {
      FlNode& fn = flNodes[pn];
      if (fn.probeState != 0) return;
      ++flResProbe;
      const vec3 nH = fn.Nrep / la::length(fn.Nrep);
      // OFFSET RETRY LADDER: the probe offset is a heuristic (the gap-finder
      // can only see cenP-covering sheets); the CERTIFICATE is the authority
      // (probed delta == combinatorial jump).  A cert failure at the base
      // offset retries at 2x/4x and accepts the first certifying reading
      // (measured: GT7081 shell2's twin corridor - a boundary cell whose
      // base probes both sat inside a sub-band gap read 0|0, while 2x reads
      // the true 0|1 and certifies).
      int consistent = 0, firstA = 0, firstB = 0;
      for (const double mul : {1.0, 2.0, 4.0, 8.0, 16.0}) {
        // stay under the weld radius: beyond it the probe leaves the cell's
        // vertical neighborhood entirely
        if (mul > 1.0 && mul * fn.off > 0.5 * eps) break;
        const std::optional<int> pa =
            RobustWinding(in, fn.cenP + mul * fn.off * nH, seeds);
        const std::optional<int> pb =
            RobustWinding(in, fn.cenP - mul * fn.off * nH, seeds);
        if (!pa || !pb) {
          if (fn.probeState == 0) fn.probeState = 2;
          continue;
        }
        fn.probeWA = *pa;
        fn.probeWB = *pb;
        if (consistent == 0) {
          firstA = *pa;
          firstB = *pb;
          consistent = 1;
        } else if (consistent > 0 && *pa == firstA && *pb == firstB) {
          ++consistent;
        } else {
          consistent = -99;  // readings diverge across the ladder
        }
        const bool cert = fn.pancake ? (*pa == *pb) : (*pb - *pa == fn.jump);
        fn.probeState = cert ? 1 : 3;
        if (cert) {
          fn.anchorOK = mul == 1.0;
          break;
        }
      }
      // LADDER-CONSISTENT NON-BOUNDARY acceptance: identical readings at
      // every offset with wA == wB mean the vertical crosses nothing - the
      // combinatorial jump was an illusory inclusive-coverage graze (the
      // near-collinear corridor slivers, measured on GT7081 shell2).  The
      // cell drops as non-boundary; a WRONG drop cannot be silent (the
      // unpaired neighbours surface at SplitTouchingSheets - fail-closed).
      if (fn.probeState == 3 && consistent >= 2 && fn.probeWA == fn.probeWB &&
          !fn.pancake)
        fn.probeState = 1;
      fn.certified = fn.probeState == 1;
    };
    // RESIDUAL-PROBE ANCHORS: in node order, the first certifiable cell of
    // each unreached subgraph anchors it (typically ONE probe per connected
    // component - the perf collapse; ladder-class cells whose edges were
    // dropped each probe individually, exactly the pre-flood cost there).
    for (int n = 0; n < nN; ++n) {
      if (E[n] != kFlUnset) continue;
      flProbe(n);
      if (flNodes[n].probeState == 1 && flNodes[n].anchorOK) {
        ++flResAnchor;
        flAnchor(n, flNodes[n].probeWA + flNodes[n].sPosMid);
        flDrain();
      }
    }
    // THE HULL SEED, guarded, for subgraphs where every probe grazed
    // (the seed's corner-to-cenP transport is refuted on twin corners -
    // GT7863 measured - so probes take precedence everywhere they certify).
    for (const auto& s : flSeeds)
      if (E[s.first] == kFlUnset) flAnchor(s.first, s.second);
    flDrain();
    // BFS conflict: clear the field - the emission below probes EVERY cell
    // (the pre-flood per-cell semantics; fail-closed intact).
    if (flMismatch > 0) std::fill(E.begin(), E.end(), kFlUnset);
    // FLOOD-PRIMARY DIFFERENTIAL (E1_FLOODDIFF): field vs the shadow probes
    int flDiffBad = 0;
    if (kFloodDiff) {
      for (int n = 0; n < nN; ++n)
        if (flNodes[n].certified && E[n] != kFlUnset &&
            E[n] != flNodes[n].probeWA + flNodes[n].sPosMid) {
          ++flDiffBad;
          if (kDump && flDiffBad <= 8)
            std::fprintf(stderr, "E1 FLOOD DIFF n=%d g=%d E=%d cert=%d\n", n,
                         flNodes[n].g, E[n],
                         flNodes[n].probeWA + flNodes[n].sPosMid);
        }
    }
    // ---- EMISSION off the field (per-cell probes where unvalued) ----
    for (int n = 0; n < nN; ++n) {
      FlNode& fn = flNodes[n];
      if (!fn.emit) continue;
      int wA, wB;
      if (E[n] != kFlUnset) {
        wA = E[n] - fn.sPosMid;
        wB = wA + fn.jump;
      } else {
        flProbe(n);
        if (fn.probeState == 2) {
          if (kDump)
            std::fprintf(stderr, "E1 FLOOD GRAZE n=%d g=%d\n", n, fn.g);
          return fail("e1: winding probe filter-uncertain - fail-closed");
        }
        if (fn.probeState == 3) {
          // certificate failed: a sub-weld-width sliver welds away; a
          // macro cell is an honest completeness failure
          if (fn.extW <= 0.99 * eps) continue;
          if (kDump) {
            std::fprintf(stderr,
                         "E1 FAIL cert n=%d g=%d jump=%d own=%d wA=%d wB=%d "
                         "ext=%.3g off=%.3g cenP=(%.9g,%.9g,%.9g)\n",
                         n, fn.g, fn.jump, fn.ownJump, fn.probeWA, fn.probeWB,
                         fn.extW, fn.off, fn.cenP.x, fn.cenP.y, fn.cenP.z);
            for (const auto& e2 : fn.nearSign)
              std::fprintf(stderr, "    near f=%d tSign=%d q=%d\n", e2[0],
                           e2[1], e2[2]);
            for (const auto& pv : fn.pos)
              std::fprintf(stderr, "    pos id=%d (%.9g,%.9g,%.9g)\n", pv.first,
                           pv.second.x, pv.second.y, pv.second.z);
            const vec3 nH = fn.Nrep / la::length(fn.Nrep);
            for (double mul : {0.25, 0.5, 1.0, 2.0, 8.0}) {
              const std::optional<int> qa =
                  RobustWinding(in, fn.cenP + mul * fn.off * nH, seeds);
              const std::optional<int> qb =
                  RobustWinding(in, fn.cenP - mul * fn.off * nH, seeds);
              std::fprintf(stderr, "    ladder %.2f: wA=%d wB=%d\n", mul,
                           qa ? *qa : -99, qb ? *qb : -99);
            }
          }
          return fail(
              "e1: coordinated-arrangement completeness certificate failed "
              "(winding delta != covering jump) - fail-closed");
        }
        wA = fn.probeWA;
        wB = fn.probeWB;
      }
      {
        const bool aIn = wA >= 1, bIn = wB >= 1;
        if (!fn.pancake && aIn == bIn) continue;  // not a {w>=1} boundary
        if (fn.pancake) {
          // the layer certificate, verbatim from the probed path
          if (aIn) continue;
          std::vector<std::pair<double, int>> layers = fn.stackSheets;
          layers.push_back({0.0, fn.ownJump});
          std::sort(layers.begin(), layers.end(),
                    [](const std::pair<double, int>& x,
                       const std::pair<double, int>& y) {
                      return x.first > y.first;
                    });
          int w = wA;
          bool material = false;
          for (size_t k = 0; k < layers.size();) {
            size_t j = k;
            int dsum = 0;
            while (j < layers.size() && layers[j].first == layers[k].first) {
              dsum += layers[j].second;
              ++j;
            }
            w += dsum;
            if (j < layers.size() && w >= 1) material = true;
            k = j;
          }
          if (!material) continue;
        }
        const int orient =
            fn.pancake ? (fn.ownJump > 0 ? 1 : -1) : (bIn ? 1 : -1);
        ++nBoundary;
        for (const ivec3& t : fn.tris) {
          const vec3 p0 = fn.pos.at(t.x), p1 = fn.pos.at(t.y),
                     p2 = fn.pos.at(t.z);
          const vec3 nr = la::cross(p1 - p0, p2 - p0);
          const double sd = la::dot(nr, fn.Nrep);
          if ((sd > 0.0 || (sd == 0.0 && orient > 0)) == (orient > 0))
            out.push_back({{p0, p1, p2}});
          else
            out.push_back({{p0, p2, p1}});
          outG.push_back(fn.g);
          outN.push_back(static_cast<double>(orient) * fn.Nrep);
        }
      }
    }
    if (kDump)
      std::fprintf(stderr,
                   "E1 FLOOD nodes=%d edges=%d hand=%d orphan=%d seeds=%d "
                   "reached=%d mismatch=%d pairDrop=%d tagDrop=%d "
                   "resProbe=%d resAnchor=%d diffBad=%d\n",
                   nN, static_cast<int>(flEdges.size()), flHandEdges,
                   flHandOrphan, static_cast<int>(flSeeds.size()),
                   static_cast<int>(q.size()), flMismatch, flPairDrop,
                   flTagDrop, flResProbe, flResAnchor, flDiffBad);
  }
  tSolve = flMs(flTS0, flNow());
  if (kFlTime)
    std::fprintf(stderr,
                 "E1 TIME build=%d total=%.0f seam=%.0f groups=%.0f "
                 "walk=%.0f classify=%.0f solve=%.0f\n",
                 buildPhase ? 1 : 0, flMs(flT0, flNow()), tSeam, tGroups, tWalk,
                 tClassify, tSolve);
  if (kDump)
    std::fprintf(stderr,
                 "E1 emitted=%d dust=%d triFail=%d spliceFail=%d cells=%d "
                 "neg=%d jump=%d owned=%d boundary=%d dustCell=%d "
                 "pancake=%d suspectSnaps=%d\n",
                 static_cast<int>(out.size()), dustTri, triFail, spliceFail,
                 nCells, nNeg, nJump, nOwned, nBoundary, nDustCell, nPancake,
                 suspectTotal);
  if (buildPhase && kDump && std::getenv("E1_REGDUMP")) {
    for (const auto& kv : lineReg)
      for (const auto& e2 : kv.second) {
        const vec3& P = e2.second;
        if (std::abs(P.x + 18.0) < 1e-9 && std::abs(P.y - 2.36572) < 5e-5 &&
            P.z > -204.4 && P.z < -203.9) {
          const auto ip = regProv.find(e2.first);
          std::fprintf(stderr,
                       "E1 REG line=(%d,%d,%d) p=(%.10g,%.10g,%.10g) %s\n",
                       std::get<0>(kv.first), std::get<1>(kv.first),
                       std::get<2>(kv.first), P.x, P.y, P.z,
                       ip != regProv.end() ? ip->second : "?");
        }
      }
  }
  if (buildPhase)
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved, "e1: build phase");
  // Triangulation failures are DEMOTED (skip-and-continue at the cell; the
  // re-gate below refuses the resulting unpaired boundary) - only the
  // hole-splice failures keep their own census fatal.
  if (spliceFail > 0) return fail("e1: hole-splice incomplete - fail-closed");
  if (const char* sf = std::getenv("E1_SOUPFILE")) {  // offline diagnostics
    if (FILE* fp = std::fopen(sf, "w")) {
      for (size_t k = 0; k < out.size(); ++k) {
        const OutTri3D& t = out[k];
        std::fprintf(fp, "%d %la %la %la %la %la %la %la %la %la\n", outG[k],
                     t.v[0].x, t.v[0].y, t.v[0].z, t.v[1].x, t.v[1].y, t.v[1].z,
                     t.v[2].x, t.v[2].y, t.v[2].z);
      }
      std::fclose(fp);
    }
  }
  // Sheet provenance rides through the weld: the radial branch's fan pairing
  // is provenance-driven where the chord angle is weld-bent noise.
  StageResult<Manifold::Impl> built = BuildImpl(out, eps, &outN);
  // Demoted-skip breadcrumb: when skipped cells leave the boundary unpaired,
  // the re-gate's fatal should name the true first cause.
  if (!built.ok() && triFail > 0)
    built.detail += " [" + std::to_string(triFail) +
                    " arrangement cell(s) skipped after exact-triangulation "
                    "failure]";
  return built;
}

// Two-pass driver: pass 1 collects the failing walks' observed self-crossing
// points (residual drawn-graph non-planarity at the representable-thin
// scale); pass 2 injects them as group-global splits.  Converges because the
// injected points lie ON the crossing sub-edges (both incident cells split
// identically); a second failure is an honest fail-closed.
StageResult<Manifold::Impl> EmitCoordinatedBoundary(const Manifold::Impl& in,
                                                    const BuildArrangement& A,
                                                    double eps) {
  // BUILD to a global registry fixpoint, then CONSUME once.
  std::map<std::tuple<int, int, int>, std::map<e1::K3, vec3>> lineReg;
  std::vector<std::vector<E1Seam>> seamCache;  // enumerated once, reused
  size_t prev = static_cast<size_t>(-1);
  for (int round = 0; round < 8; ++round) {
    EmitCoordinatedBoundaryImpl(in, A, eps, lineReg, seamCache, true);
    size_t sz = 0;
    for (const auto& kv : lineReg) sz += kv.second.size();
    if (sz == prev) break;  // registry stable
    prev = sz;
  }
  return EmitCoordinatedBoundaryImpl(in, A, eps, lineReg, seamCache, false);
}

// NEAR-COPLANAR WIDEN + GLOBAL-PLANARITY GUARD (docs/Regularize3D.md stage-5;
// reg3d-nearcoplanar-research candidate (a)).  The exact-coplanar fold only
// admits pairs whose six cross orient3d filter signs are all 0 (coplanarity gap
// below ~1 ULP).  Faces within eps of coplanar but ABOVE that bound are
// DECIDABLE non-coplanar yet their arrangement cells are sub-eps thin:
// enumerated transversally they double-round to slivers (unresolvable sheet
// contact) - the thin band this pass closes.
//
// This is an INPUT-SIDE PLANARIZATION, run BEFORE the resolver's
// enumeration/winding/emit, so the resolver RE-DERIVES the whole arrangement
// from the snapped input (the thin cell ceases to exist).  It is NOT an
// emission-time snap (those fight decisions the arrangement already made,
// ExactArrangement3D variant-E kill); it perturbs the INPUT by <= eps inside
// the standing epsilon-valid contract, coordinated by construction.
//
//  1. WIDEN: union bbox-overlapping, non-self-adjacent faces that OVERLAP in 2D
//     and whose max cross vertex-plane distance is < eps (the near band).  A
//     cluster is EXACT (skip - the existing fold owns it, bitwise) unless some
//     admitted pair is filter-non-coplanar (a genuine near-band pair).
//  2. GLOBAL-PLANARITY GUARD (the curvature safety net + snap target): fit ONE
//     plane to a near cluster (centroid + area-weighted normal); if any member
//     vertex deviates > eps, FAIL CLOSED (a curved near-coplanar chain,
//     distinct from the SoS residue) - a strictly-narrower refusal than today's
//     blanket.
//  3. SNAP each near cluster's verts onto its fitted plane (<= eps move). After
//     the snap the cluster is EXACTLY coplanar (to ~1 ULP), so the landed exact
//     fold in EmitComponentBoundary handles it verbatim and the m ==
//     winding-jump self-check holds by the exact argument.
//
// Returns: {value} = the snapped copy when a near cluster was snapped;
//          {} (no value, no fatal) when there is no near-band cluster (the
//          caller runs on the ORIGINAL input, bitwise); Fatal on a guard
//          failure or a vertex shared by two near clusters (an inconsistent
//          snap).
StageResult<Manifold::Impl> SnapNearCoplanarClusters(const Manifold::Impl& in,
                                                     double eps) {
  const int nTri = static_cast<int>(in.NumTri());
  const TriSoup soup(in);
  const auto& tri = soup.tri;
  const auto& vid = soup.vid;
  // Unit normal of face f (or false if degenerate).
  auto unitN = [&](int f, vec3& n) {
    const vec3 raw = la::cross(tri[f][1] - tri[f][0], tri[f][2] - tri[f][0]);
    const double l = la::length(raw);
    if (!(l > 0.0)) return false;
    n = raw / l;
    return true;
  };
  // Local coplanarity gap of the pair: the smaller of the two directional maxes
  // (each triangle's verts to the OTHER's plane).  The min direction is the
  // smaller face's verts against the larger plane = the true gap over the
  // shared footprint; the max direction extrapolates one plane across the
  // other's full extent (diameter-amplified, a red herring per the research
  // memo).  A macro transversal crossing has BOTH directions large, so min
  // still rejects it. Returns +inf if either plane is degenerate.
  auto pairGap = [&](int i, int j) {
    vec3 ni, nj;
    if (!unitN(i, ni) || !unitN(j, nj)) {
      return std::numeric_limits<double>::infinity();
    }
    double di = 0.0, dj = 0.0;
    for (int k = 0; k < 3; ++k) {
      di = std::max(di, std::abs(la::dot(ni, tri[j][k] - tri[i][0])));
      dj = std::max(dj, std::abs(la::dot(nj, tri[i][k] - tri[j][0])));
    }
    return std::min(di, dj);
  };
  DisjointSets uf(nTri);
  std::vector<char> nearFace(nTri, 0);  // touched a genuine near-band pair
  bool anyNear = false;
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (!soup.BBoxOverlap(i, j) || soup.SharesVert(i, j)) continue;
      if (pairGap(i, j) >= eps) continue;
      if (!TrianglesOverlap2D(tri[i], tri[j])) continue;
      uf.unite(i, j);
      if (!FacesFilterCoplanar(tri[i],
                               tri[j])) {  // near band, not exact-coplanar
        nearFace[i] = nearFace[j] = 1;
        anyNear = true;
      }
    }
  // No near-band cluster: the caller runs on the ORIGINAL input, bitwise.  The
  // exact path is thus wholly unperturbed by this pass.
  if (!anyNear) return StageResult<Manifold::Impl>{};

  // Group faces by root; a cluster is NEAR (to be snapped) iff it holds a
  // near-band face, EXACT (left bitwise for the existing fold) otherwise.
  std::map<int, std::vector<int>> members;
  for (int f = 0; f < nTri; ++f)
    members[static_cast<int>(uf.find(f))].push_back(f);

  // Snap into a copy; reads stay against the pristine `in`.  A vertex in two
  // NEAR clusters is an inconsistent snap (fail closed, never silently pick
  // one).
  Manifold::Impl work = in;
  std::vector<int> vertCluster(in.NumVert(), -1);
  int clusterId = 0;
  for (auto& [root, faces] : members) {
    if (faces.size() < 2) continue;
    bool clusterNear = false;
    for (int f : faces) clusterNear = clusterNear || nearFace[f];
    if (!clusterNear) continue;  // exact cluster: leave bitwise

    // Fitted plane: area-weighted normal (aligned to face 0 so anti-oriented
    // members do not cancel) through the member-vertex centroid.
    std::set<int> verts;
    for (int f : faces)
      for (int k = 0; k < 3; ++k) verts.insert(vid[f][k]);
    // faces[0] is non-degenerate by construction: a degenerate face has
    // pairGap == +inf (unitN false), so it never unites and stays a singleton
    // (skipped at faces.size() < 2), never entering a size>=2 cluster.
    vec3 ref;
    unitN(faces[0], ref);
    vec3 nSum(0.0, 0.0, 0.0);
    for (int f : faces) {
      const vec3 raw = la::cross(tri[f][1] - tri[f][0], tri[f][2] - tri[f][0]);
      nSum += (la::dot(raw, ref) < 0.0) ? -raw : raw;
    }
    // nl > 0 always: every raw is sign-aligned to the unit ref, and the
    // faces[0] term contributes dot(raw_0, ref) = |raw_0| > 0, so
    // dot(nSum, ref) >= |raw_0| > 0 and thus |nSum| > 0.
    const double nl = la::length(nSum);
    const vec3 N = nSum / nl;
    vec3 cen(0.0, 0.0, 0.0);
    for (int v : verts) cen += in.vertPos_[v];
    cen /= static_cast<double>(verts.size());

    // GLOBAL-PLANARITY GUARD: max member deviation from the fitted plane must
    // be within eps, else this is a curved near-coplanar chain the fold cannot
    // carry (the anti-chain-reaction net).  Fail closed, distinctly named.
    for (int v : verts)
      if (std::abs(la::dot(N, in.vertPos_[v] - cen)) > eps) {
        return StageResult<Manifold::Impl>::Fatal(
            FatalReason::DirtyComponentUnresolved,
            "resolver: near-coplanar cluster fails the global-planarity "
            "guard (curved chain, max deviation > eps) - fail-closed");
      }
    // SNAP onto the fitted plane; flag an inconsistent multi-cluster vertex.
    for (int v : verts) {
      if (vertCluster[v] >= 0 && vertCluster[v] != clusterId) {
        return StageResult<Manifold::Impl>::Fatal(
            FatalReason::DirtyComponentUnresolved,
            "resolver: a vertex lies in two near-coplanar clusters "
            "(inconsistent snap) - fail-closed");
      }
      vertCluster[v] = clusterId;
      work.vertPos_[v] = in.vertPos_[v] - la::dot(N, in.vertPos_[v] - cen) * N;
    }
    ++clusterId;
  }

  work.CalculateBBox();
  work.epsilon_ = in.epsilon_;
  return StageResult<Manifold::Impl>::Ok(std::move(work));
}

// The resolver (docs/Regularize3D.md "the resolver's mechanism") - the
// dirty-core resolver. The validated MECHANISM (enumeration + coupled winding)
// is ported; THE BUILD (the {w_S>=1} halfedge boundary emission) reuses
// RemoveOverlaps2D per crossed face (EmitComponentBoundary above).  the
// resolver enumerates + records the seam geometry, runs the build, and
// re-gates; anything it cannot resolve exactly (an exact-zero pierce tie =
// single-global SoS, a >2-sheet triple point, a coplanar seam, a
// negative-winding patch) FAILS CLOSED with a named reason - never a silent
// wrong result.  The caller re-gates the output once more (IsSelfIntersecting).
StageResult<Manifold::Impl> ResolveComponent(const Manifold::Impl& dirty,
                                             double eps) {
  // NEAR-COPLANAR PRE-PASS (docs/Regularize3D.md stage-5): planarize any
  // within-eps near-coplanar overlap cluster onto its fitted plane so the
  // resolver re-derives the arrangement from an exactly-coplanar input.  No
  // near cluster
  // -> `dirty` is used bitwise (the exact path is unperturbed); a curved chain
  // (global-planarity guard failure) fails closed here, distinctly named.
  StageResult<Manifold::Impl> snapped = SnapNearCoplanarClusters(dirty, eps);
  if (snapped.fatal) return snapped;
  const Manifold::Impl& in = snapped.value ? *snapped.value : dirty;

  // Exactly-coplanar face clusters are resolved by the in-plane fold; the
  // transversal seam enumeration skips their pairs.  The single-global SoS
  // decides every non-coplanar transversal exact-zero tie (pierce / no-pierce),
  // so there is no "SoS refused" fail-closed slot left here - a residue instead
  // surfaces as a non-2-endpoint seam (A.ok) or downstream at emission.
  const std::vector<int> face2cluster = DetectCoplanarClusters(in);
  BuildArrangement A = RecordSeams(in, face2cluster, eps);
  if (!A.ok) {
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "resolver: a self-crossing pair had a non-2-endpoint seam "
        "(degenerate incidence) - fail-closed");
  }
  // Hoist planeId (coplanar clusters collapsed to one id) once, so the seam /
  // wedge / provenance / emit passes share it instead of recomputing locally.
  {
    const int nTri = static_cast<int>(A.tri.size());
    A.planeId.assign(nTri, 0);
    for (int f = 0; f < nTri; ++f)
      A.planeId[f] = face2cluster[f] >= 0 ? nTri + face2cluster[f] : f;
  }
  // B1: enumerate the once-only 3-face triple points before per-face emission
  // (no-op off openscad; the whole point on the triple-point-dense soup).
  EnumerateTriplePoints(A, eps);
  // EX2 (exact2d): register the GAP-FREE wedge crossings the default gate
  // declines and flag their incident faces for the exact PSLG overlay.  Runs
  // AFTER the default enumeration (adds to A.seamTriples) and BEFORE the
  // registry (so the registry picks up the wedge split points).  No-op off
  // openscad (zero gap-free crossings on every resolving carrier).
  EnumerateWedgeSplits(A, eps);
  // f4-junction: gather the once-only junction registry (seam endpoints +
  // triples) so every emit path splits its edges at the non-proper-crossing
  // junctions the triple enumeration misses (no-op off openscad).
  BuildJunctionRegistry(A, eps);
  // THE FLIP (owner-ordered, flip arc stage 4): the coordinated engine is THE
  // emission for every dirty component - the per-face path, its E1_ENGINE
  // fallback gate and the E1_FLIP preview lever are deleted.  Fail-closed is
  // preserved: the engine's own certificates + the re-gate refuse with named
  // reasons; there is no second code path.
  return EmitCoordinatedBoundary(in, A, eps);
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

// Test hook (overlap3.h): expose the micro exact tie-test so the property pin
// can grade it directly.
int Orient3DExactSignProbe(const vec3& a, const vec3& b, const vec3& c,
                           const vec3& d) {
  return Orient3DExactSign(a, b, c, d);
}

// EXACT RE-GATE ARM (filter-first).  GateComponent's heuristics
// (IsSelfIntersecting's PhysX distance test, HasCoplanarOverlap's witness)
// are the FILTER; when they flag resolver output, THIS is the arbiter.
// Contract: a contact whose interpenetration is below the weld radius is
// COINCIDENT/TOUCHING - a valid double-precision rendering of geometry the
// representation cannot separate (measured on the graded openscad resolve:
// max crossing depth 0.23 eps, coplanar lens widths <= 4e-7 eps); anything
// at or beyond eps is a genuine violation and refuses.  TOPOLOGICAL
// decisions (does an edge strictly cross a triangle; is a coplanar pair
// overlapping) use the exact kernel (Orient3DFilterSign ->
// Orient3DExactSign, ExactOrient2DDrop - landed forms only); only the eps-
// scale MAGNITUDES (depth, lens width) are double arithmetic, with four-plus
// orders of margin to the threshold.  Adjacency needs no special casing:
// shared vertices/edges yield orient zeros and strictness rejects them.
// Returns true iff every contact is within-weld (the component may pass).
static bool RegateContactsWithinWeld(const Manifold::Impl& m, double eps) {
  const size_t nTri = m.halfedge_.size() / 3;
  std::vector<std::array<vec3, 3>> tri(nTri);
  for (size_t t = 0; t < nTri; ++t)
    for (int k = 0; k < 3; ++k)
      tri[t][k] = m.vertPos_[m.halfedge_.Start(3 * t + k)];
  auto o3 = [](const vec3& a, const vec3& b, const vec3& c,
               const vec3& d) -> int {
    const int s = Orient3DFilterSign(a, b, c, d);
    return s != 0 ? s : Orient3DExactSign(a, b, c, d);
  };
  // uniform grid over the triangle boxes
  const double H = 0.5;
  std::map<std::tuple<long long, long long, long long>, std::vector<int>> grid;
  std::vector<Box> boxes(nTri);
  for (size_t t = 0; t < nTri; ++t) {
    Box b;
    for (int k = 0; k < 3; ++k) b.Union(tri[t][k]);
    boxes[t] = b;
    for (long long x = static_cast<long long>(std::floor(b.min.x / H));
         x <= static_cast<long long>(std::floor(b.max.x / H)); ++x)
      for (long long y = static_cast<long long>(std::floor(b.min.y / H));
           y <= static_cast<long long>(std::floor(b.max.y / H)); ++y)
        for (long long z = static_cast<long long>(std::floor(b.min.z / H));
             z <= static_cast<long long>(std::floor(b.max.z / H)); ++z)
          grid[{x, y, z}].push_back(static_cast<int>(t));
  }
  // strict segment-triangle crossing + penetration depth of the edge tip
  auto crossingBeyondWeld = [&](const std::array<vec3, 3>& A,
                                const std::array<vec3, 3>& B) -> bool {
    for (int k = 0; k < 3; ++k) {
      const vec3 &p = A[k], &q = A[(k + 1) % 3];
      const int s1 = o3(B[0], B[1], B[2], p);
      const int s2 = o3(B[0], B[1], B[2], q);
      if (s1 == 0 || s2 == 0 || s1 == s2) continue;
      const int t0 = o3(p, q, B[0], B[1]);
      const int t1 = o3(p, q, B[1], B[2]);
      const int t2 = o3(p, q, B[2], B[0]);
      if (t0 == 0 || t0 != t1 || t1 != t2) continue;
      // strict crossing: how far does the offending edge extend past B's
      // plane?  (a lower bound of the surfaces' interpenetration)
      const vec3 nB = la::cross(B[1] - B[0], B[2] - B[0]);
      const double nl = la::length(nB);
      if (!(nl > 0.0)) continue;
      const double depth = std::min(std::abs(la::dot(nB, p - B[0])),
                                    std::abs(la::dot(nB, q - B[0]))) /
                           nl;
      if (depth >= eps) return true;
    }
    return false;
  };
  // exact-coplanar 2D overlap + lens width
  auto coplanarBeyondWeld = [&](const std::array<vec3, 3>& A,
                                const std::array<vec3, 3>& B) -> bool {
    for (int k = 0; k < 3; ++k)
      if (o3(A[0], A[1], A[2], B[k]) != 0) return false;
    const vec3 nA = la::cross(A[1] - A[0], A[2] - A[0]);
    const int ax = DominantAxis(nA);
    auto d2 = [&](const vec3& v) -> vec2 {
      if (ax == 0) return {v.y, v.z};
      if (ax == 1) return {v.x, v.z};
      return {v.x, v.y};
    };
    auto o2 = [&](const vec3& a, const vec3& b, const vec3& c) -> int {
      return ExactOrient2DDrop(a, b, c, ax);
    };
    auto strictInside = [&](const vec3& p, const std::array<vec3, 3>& T) {
      const int a = o2(T[0], T[1], p), b = o2(T[1], T[2], p),
                c = o2(T[2], T[0], p);
      return (a > 0 && b > 0 && c > 0) || (a < 0 && b < 0 && c < 0);
    };
    auto properX = [&](const vec3& a, const vec3& b, const vec3& c,
                       const vec3& d) {
      const int d1 = o2(a, b, c), da = o2(a, b, d), d3 = o2(c, d, a),
                d4 = o2(c, d, b);
      return d1 != 0 && da != 0 && d3 != 0 && d4 != 0 && d1 != da && d3 != d4;
    };
    bool overlap = false;
    for (int k = 0; k < 3 && !overlap; ++k)
      overlap = strictInside(A[k], B) || strictInside(B[k], A);
    for (int k = 0; k < 3 && !overlap; ++k)
      for (int j = 0; j < 3 && !overlap; ++j)
        overlap = properX(A[k], A[(k + 1) % 3], B[j], B[(j + 1) % 3]);
    if (!overlap) return false;
    // lens width via Sutherland-Hodgman clip (doubles about a local origin;
    // error orders below the eps threshold)
    const vec2 o = d2(A[0]);
    std::vector<vec2> P, Q;
    for (int k = 0; k < 3; ++k) P.push_back(d2(A[k]) - o);
    for (int k = 0; k < 3; ++k) Q.push_back(d2(B[k]) - o);
    auto ccw = [](std::vector<vec2>& R) {
      double s = 0;
      for (size_t k = 0; k < R.size(); ++k) {
        const vec2 &u = R[k], &v = R[(k + 1) % R.size()];
        s += u.x * v.y - v.x * u.y;
      }
      if (s < 0) std::reverse(R.begin(), R.end());
    };
    ccw(P);
    ccw(Q);
    std::vector<vec2> out = P;
    for (size_t k = 0; k < Q.size() && !out.empty(); ++k) {
      const vec2 &Aq = Q[k], &Bq = Q[(k + 1) % Q.size()];
      std::vector<vec2> cur = std::move(out);
      out.clear();
      for (size_t mI = 0; mI < cur.size(); ++mI) {
        const vec2 &C = cur[mI], &D = cur[(mI + 1) % cur.size()];
        const double dc =
            (Bq.x - Aq.x) * (C.y - Aq.y) - (Bq.y - Aq.y) * (C.x - Aq.x);
        const double dd =
            (Bq.x - Aq.x) * (D.y - Aq.y) - (Bq.y - Aq.y) * (D.x - Aq.x);
        if (dc >= 0) out.push_back(C);
        if ((dc > 0 && dd < 0) || (dc < 0 && dd > 0)) {
          const double tt = dc / (dc - dd);
          out.push_back({C.x + tt * (D.x - C.x), C.y + tt * (D.y - C.y)});
        }
      }
    }
    if (out.size() < 3) return false;
    double s = 0;
    for (size_t k = 0; k < out.size(); ++k) {
      const vec2 &u = out[k], &v = out[(k + 1) % out.size()];
      s += u.x * v.y - v.x * u.y;
    }
    const double area = std::abs(0.5 * s);
    if (!(area > 0.0)) return false;
    vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
    for (int k = 0; k < 3; ++k) {
      lo = la::min(la::min(lo, d2(A[k]) - o), d2(B[k]) - o);
      hi = la::max(la::max(hi, d2(A[k]) - o), d2(B[k]) - o);
    }
    const double ext = std::max(hi.x - lo.x, hi.y - lo.y);
    const double width = ext > 0.0 ? area / ext : 0.0;
    return width >= eps;
  };
  for (const auto& kv : grid) {
    const std::vector<int>& v = kv.second;
    for (size_t a = 0; a < v.size(); ++a)
      for (size_t b = a + 1; b < v.size(); ++b) {
        const int i = v[a], j = v[b];
        if (boxes[i].min.x > boxes[j].max.x + eps ||
            boxes[j].min.x > boxes[i].max.x + eps ||
            boxes[i].min.y > boxes[j].max.y + eps ||
            boxes[j].min.y > boxes[i].max.y + eps ||
            boxes[i].min.z > boxes[j].max.z + eps ||
            boxes[j].min.z > boxes[i].max.z + eps)
          continue;
        if (crossingBeyondWeld(tri[i], tri[j]) ||
            crossingBeyondWeld(tri[j], tri[i]) ||
            coplanarBeyondWeld(tri[i], tri[j])) {
          if (std::getenv("E1_DUMP") != nullptr)
            std::fprintf(stderr, "E1 EXACTARM violation t%d x t%d\n", i, j);
          return false;
        }
      }
  }
  return true;
}

RegularizeResult RemoveOverlaps3D(const Manifold::Impl& in, double eps) {
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

  // Components are independent by contract (no cross-component weld), so
  // gate+resolve is order-free.  A single two-pass path covers every component
  // count: resolve each component in parallel into a slot (manifold::for_each_n
  // / autoPolicy), then reduce in index order.  This keeps the counters,
  // firstFatal, and output order bitwise-identical to a sequential loop; at one
  // component autoPolicy is Seq, so it degenerates to exactly that loop and the
  // reduce over a single slot reproduces the old single-component result.
  // (meshID is assigned once, sequentially, in ComposeComponents, so the
  // atomic-counter order is output-invariant.)
  {
    enum { KClean, KReg, KFail };
    struct CompOut {
      int kind = KFail;
      bool dirty = false;
      std::optional<FatalReason> fatal;
      std::string detail;
      Manifold::Impl impl;
    };
    std::vector<CompOut> co(components.size());
    for_each_n(autoPolicy(components.size(), 1), countAt(0), components.size(),
               [&](int i) {
                 Manifold::Impl& comp = components[i];
                 const GateVerdict verdict = GateComponent(comp);
                 // Reachable defensive guard on the direct-Impl path only: the
                 // public Impl(MeshGL64) ctor sanitizes a non-manifold soup to
                 // empty (early-returns before here), but an internally-built
                 // Impl can present a non-manifold component. Fail closed.
                 if (verdict == GateVerdict::Invalid) {
                   co[i].fatal = FatalReason::NonManifoldEmission;
                   co[i].detail = "input component is not 2-manifold";
                   return;
                 }
                 if (verdict == GateVerdict::Clean) {
                   co[i].kind = KClean;
                   co[i].impl = std::move(comp);
                   return;
                 }
                 co[i].dirty = true;
                 StageResult<Manifold::Impl> bRes = ResolveComponent(comp, eps);
                 if (!bRes.ok()) {
                   co[i].fatal = bRes.fatal;
                   co[i].detail = std::move(bRes.detail);
                   return;
                 }
                 Manifold::Impl bImpl = std::move(*bRes.value);
                 bImpl.epsilon_ = eps;
                 // RE-GATE the resolver's output once (same gate as the input;
                 // the resolver's coords are double-rounded).  A clean pass
                 // composes in; a failure is the honest fail-closed, never a
                 // silent wrong result.  This is a PRODUCTION fail-closed
                 // backstop for the R1/R2 weld-fold blind spot (BuildImpl
                 // already gates non-manifold emission; the re-gate's
                 // non-redundant job is catching a MANIFOLD-but-self-
                 // intersecting output = a weld-manufactured fold).  It is
                 // verified UNREACHED on constructible general-position
                 // fixtures (reg3d-s3: 0/180 sphere variants produce
                 // re-gate-catchable output - every bad case is caught earlier
                 // by BuildImpl's manifold gate), i.e. it fires only in the
                 // unbuilt weld-fold regime.  It is deliberately NOT demoted to
                 // a DEBUG_ASSERT: it must fail closed in RELEASE, not compile
                 // out and admit wrong geometry.
                 // FILTER-FIRST re-gate: the heuristic gate flags, the EXACT
                 // arm arbitrates (owner-blessed).  Manifoldness must hold
                 // regardless; a flagged component passes only when every
                 // contact is within the weld radius (sub-eps coincidence -
                 // the correct rendering of unseparable geometry), and any
                 // genuine crossing/overlap at or beyond eps refuses exactly
                 // as before.  Placement is SCOPED to the resolver's re-gate:
                 // the shared IsSelfIntersecting keeps its heuristic
                 // sensitivity for input ROUTING (a false-positive there
                 // routes to the resolver, the safe direction) and for its
                 // external test assertions.
                 const GateVerdict rgv = GateComponent(bImpl);
                 const bool regateOk = rgv == GateVerdict::Clean ||
                                       (rgv == GateVerdict::Dirty &&
                                        RegateContactsWithinWeld(bImpl, eps));
                 if (!regateOk) {
                   if (std::getenv("E1_DUMP") != nullptr)
                     std::fprintf(
                         stderr, "E1 REGATE manifold=%d selfx=%d coplanar=%d\n",
                         (bImpl.IsManifold() && bImpl.Is2Manifold()) ? 1 : 0,
                         bImpl.IsSelfIntersecting() ? 1 : 0,
                         HasCoplanarOverlap(bImpl) ? 1 : 0);
                   if (const char* rf = std::getenv("E1_REGATEDUMP")) {
                     if (FILE* fp = std::fopen(rf, "w")) {
                       for (size_t t = 0; t < bImpl.halfedge_.size() / 3; ++t) {
                         for (int k = 0; k < 3; ++k) {
                           const vec3& p =
                               bImpl.vertPos_[bImpl.halfedge_.Start(3 * t + k)];
                           std::fprintf(fp, "%la %la %la ", p.x, p.y, p.z);
                         }
                         std::fprintf(fp, "\n");
                       }
                       std::fclose(fp);
                     }
                   }
                   co[i].fatal = FatalReason::NonManifoldEmission;
                   co[i].detail = "resolver output failed the re-gate";
                   return;
                 }
                 co[i].kind = KReg;
                 co[i].impl = std::move(bImpl);
               });
    for (auto& c : co) {
      if (c.dirty) ++result.counters.dirty;
      if (c.kind == KClean) {
        ++result.counters.clean;
        outComponents.push_back(std::move(c.impl));
      } else if (c.kind == KReg) {
        ++result.counters.regularized;
        outComponents.push_back(std::move(c.impl));
      } else {
        ++result.counters.failClosed;
        if (!firstFatal) {
          firstFatal = c.fatal;
          firstDetail = std::move(c.detail);
        }
      }
    }
  }

  if (firstFatal) {
    // Fail-closed: a recorded reason, no partial output.
    result.fatal = firstFatal;
    result.detail = std::move(firstDetail);
    return result;
  }

  // 6. COMPOSE BACK by concatenation (no fusion).
  result.impl = ComposeComponents(outComponents);
  if (const char* of = std::getenv("E1_OUTFILE")) {  // offline oracle grading
    if (FILE* fp = std::fopen(of, "w")) {
      const Manifold::Impl& m = *result.impl;
      for (size_t t = 0; t < m.halfedge_.size() / 3; ++t) {
        std::fprintf(fp, "%la %la %la %la %la %la %la %la %la\n",
                     m.vertPos_[m.halfedge_.Start(3 * t)].x,
                     m.vertPos_[m.halfedge_.Start(3 * t)].y,
                     m.vertPos_[m.halfedge_.Start(3 * t)].z,
                     m.vertPos_[m.halfedge_.Start(3 * t + 1)].x,
                     m.vertPos_[m.halfedge_.Start(3 * t + 1)].y,
                     m.vertPos_[m.halfedge_.Start(3 * t + 1)].z,
                     m.vertPos_[m.halfedge_.Start(3 * t + 2)].x,
                     m.vertPos_[m.halfedge_.Start(3 * t + 2)].y,
                     m.vertPos_[m.halfedge_.Start(3 * t + 2)].z);
      }
      std::fclose(fp);
    }
  }
  return result;
}

// Test hook: exercise the resolver's ported mechanism directly (enumeration +
// coupled winding) so it can be graded against the fragment's recorded numbers.
ComponentEnumProbe EnumerateComponent_Probe(const Manifold::Impl& dirty,
                                            const std::vector<vec3>& probes,
                                            const vec3& seed) {
  ComponentEnumProbe out;
  // Fold the level-0 self-crossing counts onto the real seam scan (no separate
  // enumeration pass): DetectCoplanarClusters ONCE, then RecordSeams with the
  // probe counters.  The counters grade the raw pierce test independently of
  // the cluster skip, so they reproduce the standalone enumerator's numbers.
  const std::vector<int> face2cluster = DetectCoplanarClusters(dirty);
  const double eps = EpsilonFromScale(dirty.bBox_.Scale(), 1000);
  RecordSeams(dirty, face2cluster, eps, &out.seamCount,
              &out.boundaryTouchPairs);
  for (int c : face2cluster)
    if (c >= 0) ++out.coplanarClusterFaces;
  out.probeWinding.reserve(probes.size());
  for (const vec3& p : probes) {
    const std::optional<int> w = WindingAt(dirty, p, seed);
    out.probeWinding.push_back(w.has_value() ? *w : kWindingUncertain);
  }
  return out;
}

RegularizeResult ResolveComponentDirect(const Manifold::Impl& soup,
                                        double eps) {
  RegularizeResult result;
  if (eps <= 0.0) eps = EpsilonFromScale(soup.bBox_.Scale(), 1000);
  result.counters.components = 1;
  result.counters.dirty = 1;
  StageResult<Manifold::Impl> bRes = ResolveComponent(soup, eps);
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
    result.detail = "resolver output failed the re-gate";
    return result;
  }
  ++result.counters.regularized;
  result.impl = std::move(bImpl);
  return result;
}

}  // namespace manifold
