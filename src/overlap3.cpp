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

bool SplitTouchingSheets(std::vector<vec3>& verts, Vec<ivec3>& tv) {
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
// POINT-PAIR /\ PLANE point (P3, the coplanar-cap seam endpoint): where the
// input EDGE SEGMENT (u,w) - a point-pair line, two input vertices - pierces
// plane g. W = ng.(w-u); coord = u*W + (dg - ng.u)*(w-u) so coord/W = u +
// t(w-u), t = (dg-ng.u)/(ng.(w-u)).  Arises ONLY on an exactly-coplanar
// cluster's internal edge (across-neighbor coplanar => the edge is not a plane
// pair), where a transversal seam f/\g ends and the {f,g,neighbor} triple
// degenerates; the pierce plane is the seam partner NOT on the owner's cluster
// plane.  Degree <=3 in the plane-coeff basis / <=4 in inputs (STRICTLY under
// the deg-9/20 triple), SAME homogeneous orient2d form, SAME accumulator (no
// monomial exceeds 9 plane factors when a P3 row mixes with triple rows).
// maxdepth's P3 instantiation.
inline HPoint SegPlaneHPoint(const vec3& u, const vec3& w, const vec3& ng,
                             double dg) {
  const Poly ux = PConst(u.x), uy = PConst(u.y), uz = PConst(u.z);
  const Poly ex = PSub(PConst(w.x), ux), ey = PSub(PConst(w.y), uy),
             ez = PSub(PConst(w.z), uz);  // w - u
  const Poly Gx = PConst(ng.x), Gy = PConst(ng.y), Gz = PConst(ng.z);
  const Poly W =
      PAdd(PAdd(PMul(Gx, ex), PMul(Gy, ey)), PMul(Gz, ez));  // ng.(w-u)
  const Poly ngu =
      PAdd(PAdd(PMul(Gx, ux), PMul(Gy, uy)), PMul(Gz, uz));  // ng.u
  const Poly num = PSub(PConst(dg), ngu);                    // dg - ng.u
  HPoint p;
  p.W = W;
  p.X = PAdd(PMul(ux, W), PMul(num, ex));
  p.Y = PAdd(PMul(uy, W), PMul(num, ey));
  p.Z = PAdd(PMul(uz, W), PMul(num, ez));
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
// Filter (running-EBD) companion of SegPlaneHPoint.
inline EHPoint ESegPlaneHPoint(const vec3& u, const vec3& w, const vec3& ng,
                               double dg) {
  const EBD ux = EIn(u.x), uy = EIn(u.y), uz = EIn(u.z);
  const EBD ex = ESub(EIn(w.x), ux), ey = ESub(EIn(w.y), uy),
            ez = ESub(EIn(w.z), uz);  // w - u
  const EBD Gx = EIn(ng.x), Gy = EIn(ng.y), Gz = EIn(ng.z);
  const EBD W = EAdd(EAdd(EMul(Gx, ex), EMul(Gy, ey)), EMul(Gz, ez));
  const EBD ngu = EAdd(EAdd(EMul(Gx, ux), EMul(Gy, uy)), EMul(Gz, uz));
  const EBD num = ESub(EIn(dg), ngu);
  EHPoint p;
  p.W = W;
  p.X = EAdd(EMul(ux, W), EMul(num, ex));
  p.Y = EAdd(EMul(uy, W), EMul(num, ey));
  p.Z = EAdd(EMul(uz, W), EMul(num, ez));
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

// The seed's plane-side sign per triangle, precomputed once (all clean-face
// probes share seeds[0]).  filter-then-exact, so the value the winding looks up
// is bit-identical to the live path.
std::vector<signed char> PrecomputeSeedSign(
    const std::vector<std::array<vec3, 3>>& tri, const vec3& seed) {
  const int nTri = static_cast<int>(tri.size());
  std::vector<signed char> sign(nTri);
  for (int t = 0; t < nTri; ++t) {
    int db = Orient3DFilterSign(tri[t][0], tri[t][1], tri[t][2], seed);
    if (db == 0) db = Orient3DExactSign(tri[t][0], tri[t][1], tri[t][2], seed);
    sign[t] = static_cast<signed char>(db);
  }
  return sign;
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

// RobustWinding through the winding broadphase (clean-face path).  `cands` is a
// caller-owned scratch buffer reused across queries.  Bit-identical winding
// VALUE to the walk (crossing-superset proof above), hence identical keep/drop.
std::optional<int> RobustWindingBVH(const std::vector<std::array<vec3, 3>>& tri,
                                    const TriWindBVH& bvh, const vec3& p,
                                    const std::vector<vec3>& seeds,
                                    std::vector<int>& cands,
                                    const signed char* seedSign0) {
  for (size_t i = 0; i < seeds.size(); ++i) {
    WindCandidates(bvh, p, seeds[i], cands);
    // The precomputed table is for seeds[0]; the rare seeds[1..] retries go
    // live.
    const signed char* sign = (i == 0) ? seedSign0 : nullptr;
    if (const std::optional<int> w =
            WindingAtCands(tri, p, seeds[i], cands, sign))
      return w;
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
  // SYMBOLIC PROVENANCE (symwalk lane): the exact HPoint of every arrangement
  // vertex a WEDGE face's exact-overlay cell walk needs, keyed by the canonical
  // 3D double bits (the same bits pf.add/addAt deduplicate on).  A vertex is
  // either an INPUT point (trivial, W==1) or a 3-plane CRAMER intersection
  // (triangle corner + input-vertex T-junction = trivial; interior triple point
  // + wedge split = cramer{f,g,h}; seam endpoint = cramer{f,g,edge-neighbor}).
  // Populated ONLY when a wedge face is flagged (BuildSymbolicProvenance), so
  // it is an empty no-op on the whole resolving corpus.  The symbolic cell walk
  // (ExtractCellsSymbolic) orders vertices by the LANDED HomogOrient2D on these
  // HPoints, realizing the exact arrangement regardless of rounded straddle.
  struct VProv {
    bool trivial = true;
    bool segPlane = false;  // P3: point-pair (su,sw) /\ plane (rep face rg)
    vec3 v;                 // trivial input point (W==1)
    vec3 nF, nG, nH;        // cramer: three plane normals (unnormalized)
    double dF = 0.0, dG = 0.0, dH = 0.0;  // cramer: three plane offsets
    int pf = -1, pg = -1, ph = -1;  // cramer: the three plane ids (rep faces)
    int rf = -1, rg = -1,
        rh = -1;  // cramer: the three rep FACE indices (for
                  // the INPUT-EXACT Cramer over input verts);
                  // segPlane reuses rg = the pierce plane rep.
    vec3 su, sw;  // segPlane: the two input edge vertices
  };
  std::map<std::tuple<double, double, double>, VProv> provOf;
  // Plane id -> (unnormalized normal, a point on the plane), for constructing
  // the exact crossing point of two symbolic seam sub-edges during the wedge
  // arrangement completion (symwalk).  Populated with provOf (wedge faces
  // only).
  std::map<int, std::pair<vec3, vec3>> planeTab;
  // Plane id -> rep FACE index (the canonical face whose three input vertices
  // define the plane), for the INPUT-EXACT basis: every overlay position and
  // decision is composed exactly from these input doubles, not the rounded
  // faceN.  Populated alongside planeTab (wedge faces only).
  std::map<int, int> planeTri;
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
// Input-exact HPoint of a provenance vertex (trivial input point W==1, or the
// plane triple of its three rep faces).
inline sos::BigHPoint IXProvHPoint(const BuildArrangement& A,
                                   const BuildArrangement::VProv& pv) {
  if (pv.trivial) return sos::TrivialBigHPoint(pv.v);
  if (pv.segPlane)
    return sos::SegPlaneBigHPoint(pv.su, pv.sw, A.tri[pv.rg].data());
  return IXTripleHPoint(A, pv.rf, pv.rg, pv.rh);
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

// SYMBOLIC CELL WALK (symwalk lane): ExtractCells with the ROUNDED-double atan2
// radial sort + double area sign REPLACED by the LANDED exact HomogOrient2D on
// per-vertex symbolic HPoints (eh = filter EHPoint, hp = exact HPoint).  For
// the wedge faces the rounded seam segments do NOT straddle at their exact
// interior crossing (that is why the default gate declined them), so the double
// atan2 order folds back and emits sub-eps slivers that do not pair; the exact
// angular order realizes the true arrangement regardless.  The comparator uses
// ONLY HomogOrient2D(center, a, b) on REAL arrangement points (an arbitrary
// neighbor as the reference ray, collinear ties broken by a non-collinear
// witness) - no point at infinity, no coordinate-difference sign, NO new
// predicate form. `axis` = the face-normal dominant axis; `signMul` (+/-1, from
// the CCW triangle corners) aligns HomogOrient2D's unflipped keep-pair frame to
// verts2's proj frame so the walk's handedness matches the rest of the
// pipeline.  The single unbounded outer face (uniquely most-negative double
// area, |outer| ~ model^2 >> any sub-eps sliver) is dropped; every other loop
// is a bounded cell in walked (combinatorially CCW) order.  Returns false on a
// malformed walk (fail closed).
bool ExtractCellsSymbolic(const std::vector<vec2>& pts,
                          const std::vector<sos::EHPoint>& eh,
                          const std::vector<sos::HPoint>& hp,
                          const std::vector<sos::BigHPoint>& bhp, bool ix,
                          int axis, int signMul,
                          const std::vector<std::pair<int, int>>& uedges,
                          std::vector<std::vector<int>>& cells) {
  const int n = static_cast<int>(pts.size());
  auto orient = [&](int i, int a, int b) -> int {
    int s;
    if (ix)
      s = IXOrient2D(bhp[i], bhp[a], bhp[b], axis);
    else {
      s = sos::HomogOrient2DFilter(eh[i], eh[a], eh[b], axis);
      if (s == 0) s = sos::HomogOrient2DExact(hp[i], hp[a], hp[b], axis);
    }
    return signMul * s;  // sign in verts2's CCW frame
  };
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
  std::vector<std::vector<int>> order(n);  // exact-CCW-sorted neighbors
  std::vector<std::unordered_map<int, int>> at(n);
  for (int v = 0; v < n; ++v) {
    order[v].assign(nbr[v].begin(), nbr[v].end());
    if (order[v].empty()) continue;
    const int ref = order[v][0];  // reference ray (angle 0); any neighbor works
    int witness = -1;             // a neighbor NOT collinear with ref through v
    for (int w : order[v])
      if (w != ref && orient(v, ref, w) != 0) {
        witness = w;
        break;
      }
    auto opposite = [&](int w) -> bool {  // w collinear with ref, opposite ray
      if (witness < 0) return true;       // 1D star: the other direction
      const int sw = orient(v, witness, w), sr = orient(v, witness, ref);
      return sw != 0 && sr != 0 && sw != sr;
    };
    // Ascending CCW angle from ref in [0, 2pi): group 0 = ref (angle 0), 1 =
    // CCW half (0,pi), 2 = opposite ref (pi), 3 = CW half (pi,2pi); within a
    // half a total order by orient sign (span < pi).  Only the CYCLIC order is
    // used, so the reference/branch-cut choice is immaterial.
    auto gkey = [&](int w) -> int {
      if (w == ref) return 0;
      const int s = orient(v, ref, w);
      if (s > 0) return 1;
      if (s < 0) return 3;
      return opposite(w) ? 2 : 0;
    };
    std::sort(order[v].begin(), order[v].end(), [&](int w1, int w2) {
      if (w1 == w2) return false;
      const int g1 = gkey(w1), g2 = gkey(w2);
      if (g1 != g2) return g1 < g2;
      if (g1 == 0 || g1 == 2) return false;  // single-angle groups: equal
      return orient(v, w1, w2) > 0;
    });
    for (int k = 0; k < static_cast<int>(order[v].size()); ++k)
      at[v][order[v][k]] = k;
  }
  std::set<std::pair<int, int>> visited;
  std::vector<std::vector<int>> loops;
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
      loops.push_back(std::move(loop));
    }
  int outer = -1;
  double outerArea = 0.0;
  for (int i = 0; i < static_cast<int>(loops.size()); ++i) {
    double area = 0.0;
    const std::vector<int>& L = loops[i];
    for (size_t k = 0; k < L.size(); ++k)
      area += la::cross(pts[L[k]], pts[L[(k + 1) % L.size()]]);
    if (outer < 0 || area < outerArea) {
      outer = i;
      outerArea = area;
    }
  }
  for (int i = 0; i < static_cast<int>(loops.size()); ++i)
    if (i != outer && loops[i].size() >= 3)
      cells.push_back(std::move(loops[i]));
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
inline int ExactOrient2DDrop(const vec3& p, const vec3& q, const vec3& r,
                             int axis) {
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

// SYMBOLIC PROVENANCE (symwalk lane): populate A.provOf with the exact HPoint
// of every arrangement vertex the wedge faces' symbolic cell walk will need,
// keyed by the canonical 3D double bits pf.add/addAt deduplicate on.  Runs ONLY
// when a wedge face is flagged (EnumerateWedgeSplits), so it is an empty no-op
// on the whole resolving corpus.  Three sources, registered so the
// most-specific wins by emplace (first-wins) order:
//  (a) INPUT points (trivial, W==1): every triangle corner - covers the
//  triangle
//      boundary vertices AND any seam endpoint / T-junction that IS an input
//      vertex (shares-vertex recovery, input-vertex-on-edge).
//  (b) INTERIOR TRIPLE POINTS (cramer {f,g,h}): re-enumerated exactly as
//      EnumerateTriplePoints / EnumerateWedgeSplits construct them (sorted
//      plane-triple reps -> Intersect3Planes), so the registered bits are
//      bit-identical to the seamTriples canon the emit path deduplicates on.
//  (c) SEAM ENDPOINTS: a seam on plane(f) INT plane(g) exits triangle f or g
//      across an edge e.  TRANSVERSAL neighbor => that edge's line is
//      plane(owner) INT plane(pair-across-e), so the endpoint is EXACTLY
//      plane(f) INT plane(g) INT plane(neighbor) = cramer {f,g,neighbor}.
//      COPLANAR-CLUSTER-INTERNAL neighbor (the cap-provenance case) => the
//      across-neighbor is coplanar with the owner, so the edge is a POINT-PAIR
//      line (two input verts), not a plane pair; the endpoint is the P3 point
//      (u,w) INT plane(seam-partner) = SegPlaneHPoint (formerly left
//      unregistered
//      -> the face failed closed).  The crossed edge is found geometrically.
// Junction-registry splits dedup to an actual raw member (a seam endpoint /
// triple / input vertex), so they are already covered.
void BuildSymbolicProvenance(BuildArrangement& A, const Manifold::Impl& in,
                             double eps) {
  bool anyWedge = false;
  for (char c : A.wedgeFace)
    if (c) {
      anyWedge = true;
      break;
    }
  if (!anyWedge) return;  // byte-clean: no provenance off the wedge population
  const int nTri = static_cast<int>(A.tri.size());
  using Key = std::tuple<double, double, double>;
  auto key3 = [](const vec3& p) { return Key{p.x, p.y, p.z}; };
  auto planeOf = [&](int f, vec3& n, double& d) {
    n = A.faceN[f];
    d = la::dot(n, A.tri[f][0]);
  };
  // (a) input points -> trivial.
  for (int f = 0; f < nTri; ++f)
    for (int k = 0; k < 3; ++k) {
      BuildArrangement::VProv pv;
      pv.trivial = true;
      pv.v = A.tri[f][k];
      A.provOf.emplace(key3(A.tri[f][k]), pv);
    }
  const std::vector<int>& planeId = A.planeId;  // hoisted (ResolveComponent)
  std::map<int, int> rep;
  for (int f = 0; f < nTri; ++f) rep.emplace(planeId[f], f);
  for (const auto& [id, r] : rep) {
    A.planeTab.emplace(id, std::make_pair(A.faceN[r], A.tri[r][0]));
    A.planeTri.emplace(id, r);  // input-exact: plane id -> rep face
  }
  auto regCramer = [&](const Key& k, int r0, int r1, int r2) {
    BuildArrangement::VProv pv;
    pv.trivial = false;
    planeOf(r0, pv.nF, pv.dF);
    planeOf(r1, pv.nG, pv.dG);
    planeOf(r2, pv.nH, pv.dH);
    pv.pf = planeId[r0];
    pv.pg = planeId[r1];
    pv.ph = planeId[r2];
    pv.rf = r0;
    pv.rg = r1;
    pv.rh = r2;
    A.provOf.emplace(k, pv);  // trivial (a) already there wins
  };
  // P3: a coplanar-cluster-internal seam endpoint = point-pair line (u,w)
  // piercing the rep face rPierce of the seam-partner plane (the cap-provenance
  // identity).
  auto regSegPlane = [&](const Key& k, const vec3& u, const vec3& w,
                         int rPierce) {
    BuildArrangement::VProv pv;
    pv.trivial = false;
    pv.segPlane = true;
    pv.su = u;
    pv.sw = w;
    pv.rg = rPierce;
    A.provOf.emplace(k, pv);  // trivial (a) already there wins
  };
  // (b) interior triple points -> cramer (sorted plane-triple reps).
  for (int f = 0; f < nTri; ++f) {
    if (!A.seamed[f]) continue;
    const int ns = static_cast<int>(A.faceSeams[f].size());
    for (int k1 = 0; k1 < ns; ++k1) {
      const int pg = planeId[A.faceSeams[f][k1].other];
      if (pg == planeId[f]) continue;
      for (int k2 = k1 + 1; k2 < ns; ++k2) {
        const int ph = planeId[A.faceSeams[f][k2].other];
        if (ph == planeId[f] || ph == pg) continue;
        const int gF = A.faceSeams[f][k1].other, hF = A.faceSeams[f][k2].other;
        const int rF = rep[planeId[f]], rG = rep[pg], rH = rep[ph];
        if (!ExactSeamsCross(A, f, gF, hF, rF, rG, rH)) continue;
        std::array<int, 3> key = {planeId[f], pg, ph};
        std::sort(key.begin(), key.end());
        const int r0 = rep[key[0]], r1 = rep[key[1]], r2 = rep[key[2]];
        vec3 pos;
        if (!CanonTriplePos(A, r0, r1, r2, pos)) continue;
        regCramer(key3(pos), r0, r1, r2);
      }
    }
  }
  // (c) seam endpoints -> cramer {f, g, edge-neighbor} when the crossed edge is
  // a plane pair (transversal neighbor), else P3 point-pair /\ pierce-plane
  // (the coplanar-cluster-internal edge: the across-neighbor is coplanar so the
  // edge is NOT a plane pair - the cap-provenance case that formerly stayed
  // unregistered).
  auto onEdge = [&](const vec3& P, const vec3& u, const vec3& w) -> bool {
    const vec3 d = w - u;
    const double len2 = la::dot(d, d);
    if (!(len2 > 0.0)) return false;
    const double t = la::dot(P - u, d) / len2;
    if (t < -1e-6 || t > 1.0 + 1e-6) return false;
    return la::length(P - (u + t * d)) <= 8.0 * eps;  // on the edge line
  };
  for (int f = 0; f < nTri; ++f) {
    for (const BuildSeam& s : A.faceSeams[f]) {
      const int g = s.other;
      if (g < 0 || g >= nTri) continue;
      const vec3 ends[2] = {s.p0, s.p1};
      for (const vec3& P : ends) {
        if (A.provOf.count(key3(P))) continue;  // trivial / triple already
        bool done = false;
        for (int owner : {f, g}) {
          for (int e = 0; e < 3 && !done; ++e) {
            const vec3& u = A.tri[owner][e];
            const vec3& w = A.tri[owner][(e + 1) % 3];
            if (!onEdge(P, u, w)) continue;
            const int nb =
                static_cast<int>(in.halfedge_.Pair(3 * owner + e)) / 3;
            if (nb < 0 || nb >= nTri) continue;
            const int pf = planeId[f], pgId = planeId[g], pn = planeId[nb];
            bool reg = false;
            if (pn != pf && pn != pgId) {  // transversal neighbor: plane triple
              const int rF = rep[pf], rG = rep[pgId], rN = rep[pn];
              vec3 chk;
              if (Intersect3Planes(A.faceN[rF], A.tri[rF][0], A.faceN[rG],
                                   A.tri[rG][0], A.faceN[rN], A.tri[rN][0],
                                   chk)) {
                regCramer(key3(P), rF, rG, rN);
                reg = true;
              }
            }
            if (!reg) {  // coplanar-cluster-internal edge (or non-independent
                         // triple): P3 = (u,w) /\ seam-partner-plane.  The
                         // pierce plane is the seam plane NOT on the owner's
                         // cluster plane; P lies on it (seam endpoint) and on
                         // line(u,w).
              const int ownerPlane = planeId[owner];
              const int piercePlane = (ownerPlane == pf) ? pgId : pf;
              regSegPlane(key3(P), u, w, rep[piercePlane]);
            }
            done = true;
          }
          if (done) break;
        }
      }
    }
  }
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

// The registry vertices strictly interior to segment [P0,P1] in frame pf, each
// returned as {on-line 2D foot, canonical 3D vertex} ordered along the segment.
// The 2D position is the perpendicular FOOT on the segment line (NOT proj(V),
// which rounds off-line and would make RemoveOverlaps2D manufacture a fold-back
// crossing - the B1 lesson), keyed by the once-only 3D bits so every incident
// face welds onto the identical vertex.  Strict-interior in PARAMETER (never
// within eps of an endpoint -> no degenerate sub-edge) and on-line within eps.
std::vector<std::pair<vec2, vec3>> JunctionSplitsOnSegment(
    const vec2& q0, const vec2& q1, const vec3& P0, const vec3& P1,
    const std::vector<vec3>& junctions, double eps,
    const std::vector<vec3>* skipNear = nullptr) {
  // q0/q1 are the caller's 2D projection of P0/P1 (any affine in-plane frame -
  // the on-line foot below is frame-agnostic), so this helper is independent of
  // which face-flattening the driver uses.
  std::vector<std::pair<double, std::pair<vec2, vec3>>> ord;
  const vec3 d3 = P1 - P0;
  const double len2 = la::dot(d3, d3);
  if (!(len2 > 0.0)) return {};
  const double len = std::sqrt(len2);
  const double tlo = eps / len, thi = 1.0 - eps / len;
  for (const vec3& V : junctions) {
    // The split DECISION is a pure 3D on-segment test against the segment's own
    // endpoints, so the two incident faces (which see the IDENTICAL 3D
    // endpoints of a shared seam / mesh edge) split at the SAME junction set -
    // no per-frame disagreement near an endpoint (which manufactures new
    // T-junctions).  An over-inclusive collinear hit is then HARMLESS: both
    // faces take it, so the sub-edges still pair.
    const vec3 w = V - P0;
    const double t = la::dot(w, d3) / len2;
    if (!(t > tlo && t < thi)) continue;         // strictly interior in param
    if (la::length(w - t * d3) > eps) continue;  // on the segment line
    // Skip a junction the face already introduces itself: a vertex coincident
    // with one of the face's OWN seam endpoints is split into this edge by the
    // face's own RemoveOverlaps2D, so pre-splitting it would only re-tessellate
    // an already-complete arrangement (a spurious non-no-op on carriers like
    // EntangledBars whose seam endpoints land interior to their own wall
    // edges).  Only FOREIGN junctions - vertices no seam of this face ends at -
    // are the T-junctions the neighbour fails to split.
    if (skipNear) {
      bool own = false;
      for (const vec3& s : *skipNear)
        if (la::length(V - s) <= eps) {
          own = true;
          break;
        }
      if (own) continue;
    }
    // 2D position = the on-line foot in THIS face's frame (proj is affine, so
    // proj(P0 + t*d3) = q0 + t*(q1-q0)), keyed by the once-only 3D bits.  On-
    // line by construction -> no RemoveOverlaps2D fold-back (the B1 lesson).
    const vec2 foot = q0 + t * (q1 - q0);
    ord.push_back({t, {foot, V}});
  }
  std::sort(ord.begin(), ord.end(),
            [](const auto& x, const auto& y) { return x.first < y.first; });
  std::vector<std::pair<vec2, vec3>> out;
  out.reserve(ord.size());
  for (auto& e : ord) out.push_back(e.second);
  return out;
}

// Index of the largest-area triangle of `tris` over `pts` (its
// centroid/incenter is the most robust interior classify point); |2*area|
// returned in area2.
int LargestSubTri(const std::vector<ivec3>& tris, const std::vector<vec2>& pts,
                  double& area2) {
  int best = 0;
  area2 = -1.0;
  for (int t = 0; t < static_cast<int>(tris.size()); ++t) {
    const vec2 p0 = pts[tris[t].x], p1 = pts[tris[t].y], p2 = pts[tris[t].z];
    const double ar = std::abs(la::cross(p1 - p0, p2 - p0));
    if (ar > area2) {
      area2 = ar;
      best = t;
    }
  }
  return best;
}

// The both-sides {w_S>=1} boundary witness - the single retention rule shared
// by the seamed, fold, and clean emit paths.  A sub-face is on d{w_S>=1} iff
// EXACTLY ONE side is inside {w_S>=1}; the solid side fixes orientation.
// wAbove / wBelow are the +nHat / -nHat winding probes.  Returns 0 (drop: both
// sides same class), 1 (emit ORIGINAL orientation: solid on the -nHat side), or
// 2 (emit REVERSED: solid on +nHat) - matching EmitCleanFaces' status codes.
inline int BothSidesRetain(int wAbove, int wBelow) {
  const bool aboveIn = wAbove >= 1, belowIn = wBelow >= 1;
  if (aboveIn == belowIn) return 0;
  return belowIn ? 1 : 2;
}

// MEASUREMENT ONLY (f4-b1): per-branch fail-closed census for EmitSeamedFace,
// gated by F4B_DUMP, printed by EmitComponentBoundary.  Never touched in
// production (kF4BDump false).
struct F4BSeamCensus {
  int b3faces = 0, b3verts = 0, b4faces = 0, b5faces = 0, okfaces = 0;
};
static F4BSeamCensus gF4BSeam;

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
  AxisDropFrame pf;
  if (!BuildAxisDropFrame(A.faceN[f], a, pf)) {
    ok = false;
    return;
  }
  const vec3& nHat = pf.nHat;
  std::vector<vec2>& verts2 = pf.verts2;
  std::vector<vec3>& canon3 = pf.canon3;
  auto getV = [&](const vec3& P) { return pf.add(P); };
  const int iA = getV(a), iB = getV(b), iC = getV(c);
  const vec2 pa = verts2[iA], pb = verts2[iB], pc = verts2[iC];
  if (!(0.5 * la::cross(pb - pa, pc - pa) > 0.0)) {
    ok = false;  // parity flip did not yield CCW, or degenerate projection
    return;
  }
  std::vector<EdgeM> edges;
  static const std::vector<std::pair<vec2, vec3>> kNoSplits;
  // Chain vFrom -> (interior splits ordered along [P0,P1]) -> vTo.  `pre` are
  // the exact-crossing TRIPLE splits already computed for this segment (B1's
  // seam splits: the on-seam 2D crossing, on-line for BOTH crossing seams, so
  // all three incident faces weld onto the identical 3D vertex the pos2in map
  // would otherwise refuse).  The registry then adds every OTHER junction
  // strictly interior to the segment at its on-line FOOT, so the
  // non-proper-crossing T-junctions (seam-endpoint-on-edge / seam-endpoint-on-
  // seam) that open the fan get the SAME shared split on every incident face.
  // With no triple and no interior junction this is byte-identical to the plain
  // single-edge push (the whole corpus off openscad).
  // The face's own seam endpoints: on the DEFAULT path those are split into the
  // triangle edges by the face's own RemoveOverlaps2D (eps-incidence), so the
  // registry must NOT re-inject them (no-op rail).  On the EXACT-OVERLAY path
  // (wedge) there is no RemoveOverlaps2D, so the registry MUST split every edge
  // at every own seam endpoint interior to it (skipNear = nullptr) - otherwise
  // a seam endpoint landing on a triangle edge stays unsplit and the two
  // segments "cross" at it (the residual-crossing fail-closed).
  const bool wedge = f < static_cast<int>(A.wedgeFace.size()) && A.wedgeFace[f];
  std::vector<vec3> ownEnds;
  for (const auto& s : A.faceSeams[f]) {
    ownEnds.push_back(s.p0);
    ownEnds.push_back(s.p1);
  }
  const std::vector<vec3>* skipEnds = wedge ? nullptr : &ownEnds;
  // GRAZE ADJUDICATION (exact-overlay branch only): the junction registry
  // injects every registry vertex within eps of a face-f carrier - injection-
  // by-eps.  On a WEDGE face that snaps FOREIGN vertices (a triple point / seam
  // endpoint of ANOTHER face) that are only eps-NEAR plane f, not exactly on
  // it, onto face f's edges - the near-tangent grazes that form the phantom
  // cross- seam fragments halting the exact arrangement completion
  // (CARRIER_FAIL). Adjudicate by EXACTNESS instead: inject a junction ONLY if
  // it is a genuine vertex of face f's in-plane PSLG - a trivial input point /
  // registry T-junction (exactly on a face-f edge by the exact on-edge arm), or
  // a constructed triple that carries planeId[f] (exactly on plane f by
  // construction, a real interior/edge junction).  A constructed junction whose
  // triple EXCLUDES planeId[f] is a one-sided graze off plane f: DROP it here
  // (its geometry lives in the foreign face).  Any GENUINE seam-seam crossing
  // dropped this way is re-derived by the arrangement completion below with the
  // landed degree-9 form (self-healing), so no real junction is lost; the
  // skew-near-crossing grazes never re-appear (they do not cross in 3D).
  // Combinatorial (degree-0) plane membership - no new predicate, no nested
  // construction.  DEFAULT path (non-wedge) is untouched (byte-clean).
  const int pfIdAdj = A.planeId[f];
  const bool grazeAdj = std::getenv("GRAZE_ADJ_OFF") == nullptr;
  // STITCH (two-sided consistency): the graze adjudication drops a registry
  // junction that is eps-NEAR but not EXACTLY on this carrier (its provenance
  // triple excludes the carrier's two planes - a foreign seam endpoint of
  // another face rounded near this segment).  Applied ONLY to the wedge path it
  // makes the wedge face drop grazes its DEFAULT neighbour still injects (the
  // boundary subdivision diverges -> the 32-open residue).  onCarrier is
  // SYMMETRIC in the two incident faces (both need {planeId[f], carrierQ} = the
  // same plane pair), so applying it to the DEFAULT path too makes both sides
  // make the IDENTICAL keep/drop decision - the shared boundary subdivides
  // identically and the fans pair.  Byte-clean off the wedge population: provOf
  // is EMPTY there (BuildSymbolicProvenance returns early with no wedge face),
  // so onCarrier keeps every junction and the default injection is unchanged.
  const bool stitchDefault = std::getenv("STITCH_DEF_OFF") == nullptr;
  // A registry junction V belongs on face f's carrier segment (carrier line =
  // plane f INT plane carrierQ) iff it lies EXACTLY on that carrier.  A trivial
  // input point / registry T-junction is exact-on-edge by the exact on-edge arm
  // (kept).  A constructed junction lies on the carrier iff its plane triple
  // carries BOTH planeId[f] AND carrierQ; a triple that carries planeId[f] but
  // NOT carrierQ is a genuine on-plane interior junction that is only eps-NEAR
  // this carrier (v13 = an interior triple {440,447,467} 0.13 eps off edge2's
  // {467,nb} line) - injecting it forms the phantom cross-seam chord; and a
  // triple carrying NEITHER is a foreign off-plane graze.  Both are dropped.
  auto onCarrier = [&](const vec3& V, int carrierQ) -> bool {
    const auto it = A.provOf.find({V.x, V.y, V.z});
    if (it == A.provOf.end() || it->second.trivial) return true;
    const int p0 = it->second.pf, p1 = it->second.pg, p2 = it->second.ph;
    const bool hasF = p0 == pfIdAdj || p1 == pfIdAdj || p2 == pfIdAdj;
    const bool hasQ = p0 == carrierQ || p1 == carrierQ || p2 == carrierQ;
    return hasF && hasQ;
  };
  auto pushChain = [&](int vFrom, int vTo, const vec3& P0, const vec3& P1,
                       const std::vector<std::pair<vec2, vec3>>& pre,
                       int carrierQ) {
    std::vector<std::pair<vec2, vec3>> splits = pre;
    for (const auto& js : JunctionSplitsOnSegment(
             pf.proj(P0), pf.proj(P1), P0, P1, A.junctions, eps, skipEnds)) {
      if ((wedge || stitchDefault) && grazeAdj &&
          !onCarrier(js.second, carrierQ)) {
        if (std::getenv("GRAZE_DROP") != nullptr)
          std::fprintf(stderr, "GRAZE_DROP f=%d junction=(%.17g,%.17g,%.17g)\n",
                       f, js.second.x, js.second.y, js.second.z);
        continue;  // not on this carrier: drop (graze / off-carrier interior)
      }
      bool dup = false;
      for (const auto& s : splits)
        if (la::length(s.second - js.second) <= eps) {
          dup = true;
          break;
        }
      if (!dup) splits.push_back(js);
    }
    if (splits.empty()) {
      if (vFrom != vTo) edges.push_back({vFrom, vTo, 1});
      return;
    }
    const vec2 q0 = pf.proj(P0), dir = pf.proj(P1) - q0;
    const double len2 = la::dot(dir, dir);
    auto param = [&](const vec2& p) {
      return len2 > 0.0 ? la::dot(p - q0, dir) / len2 : 0.0;
    };
    std::sort(
        splits.begin(), splits.end(),
        [&](const std::pair<vec2, vec3>& x, const std::pair<vec2, vec3>& y) {
          return param(x.first) < param(y.first);
        });
    int prev = vFrom;
    for (const auto& e : splits) {
      const int v = pf.addAt(e.first, e.second);
      if (v != prev) edges.push_back({prev, v, 1});
      prev = v;
    }
    if (vTo != prev) edges.push_back({prev, vTo, 1});
  };
  // Carrier partner plane per triangle edge e = the neighbour across it (the
  // triangle edge line is plane f INT plane(nb)); -1 (no wedge carrier) off the
  // adjudication.  Edge e runs corner e -> corner (e+1)%3.
  auto edgeCarrierQ = [&](int e) -> int {
    const int nb = static_cast<int>(in.halfedge_.Pair(3 * f + e)) / 3;
    return (nb >= 0 && nb < static_cast<int>(A.planeId.size())) ? A.planeId[nb]
                                                                : -1;
  };
  pushChain(iA, iB, a, b, kNoSplits, edgeCarrierQ(0));
  pushChain(iB, iC, b, c, kNoSplits, edgeCarrierQ(1));
  pushChain(iC, iA, c, a, kNoSplits, edgeCarrierQ(2));
  for (int k = 0; k < static_cast<int>(A.faceSeams[f].size()); ++k) {
    const BuildSeam& s = A.faceSeams[f][k];
    const int v0 = getV(s.p0), v1 = getV(s.p1);
    // EnumerateTriplePoints (always run by ResolveComponent before any seamed
    // emit) sizes seamTriples to nTri, and a seamed face has index < nTri, so
    // this is never empty here.
    const std::vector<std::pair<vec2, vec3>>& pre = A.seamTriples[f][k];
    const int carrierQ =
        (s.other >= 0 && s.other < static_cast<int>(A.planeId.size()))
            ? A.planeId[s.other]
            : -1;
    pushChain(v0, v1, s.p0, s.p1, pre, carrierQ);
  }

  static const bool kF4BDump = std::getenv("F4B_DUMP") != nullptr;
  std::vector<std::pair<int, int>> uedges;
  std::vector<std::vector<int>> cells;
  bool cellsBuilt = false;
  if (wedge) {
    // EXACT PSLG OVERLAY (exact2d, the gap-free wedge closure): the pre-split
    // `edges` ARE the arrangement 1-skeleton - every interior crossing is
    // pre-split at a once-only triple point (default + the EX2 wedge splits)
    // and every vertex-on-edge incidence at a registry junction, so no two
    // non-adjacent sub-edges properly cross in the face interior.  Build uedges
    // DIRECTLY, SKIPPING RemoveOverlaps2D's MergeVerts eps-merge (which
    // collapses the sub-eps wedge split vertex into the seam endpoint, DROPPING
    // the split - the 21-open-edge residue) and its sweep (which re-crosses the
    // rounded sub-segments -> the face-467 no-preimage phantom).  VERIFY
    // completeness EXACTLY first (the landed exact orient2d on the canonical 3D
    // points, drop the face axis): any residual proper crossing means the
    // enumeration missed a triple point -> fail closed (never emit unverified
    // geometry). WELD near-coincident vertices (the cross-provenance ALIASES:
    // seam endpoints built by different pierce paths that round to
    // bit-different but geometrically identical points - the 330 the default
    // RemoveOverlaps2D MergeVerts collapses; nomerge/sliver-emit).  Union-find
    // within eps to the lowest-index rep; the gap-free wedge splits are >
    // kWedgeEnv eps from everything, so the weld never merges a genuine split.
    // This is the coordinated collapse of the GAP-BOUNDED alias population, the
    // exact-overlay analogue of MergeVerts WITHOUT the sweep that re-crosses
    // the sub-segments.
    const int nv = static_cast<int>(verts2.size());
    std::vector<int> rep(nv);
    for (int k = 0; k < nv; ++k) rep[k] = k;
    auto find = [&](int x) {
      while (rep[x] != x) x = rep[x] = rep[rep[x]];
      return x;
    };
    for (int i = 0; i < nv; ++i)
      for (int j = i + 1; j < nv; ++j)
        if (la::length(verts2[i] - verts2[j]) <= eps) {
          const int ri = find(i), rj = find(j);
          if (ri != rj) rep[std::max(ri, rj)] = std::min(ri, rj);
        }
    const bool symOff = std::getenv("SYMWALK_OFF") != nullptr;
    // SYMBOLIC HPOINTS (symwalk): the exact point of every welded-rep vertex
    // (A.provOf), on which BOTH the completeness check and the cell walk
    // decide. A wedge vertex with no symbolic provenance fails the face closed
    // (never orders / verifies on a non-symbolic point).
    std::vector<sos::EHPoint> eh(nv);
    std::vector<sos::HPoint> hp(nv);
    std::vector<sos::BigHPoint> bhp(nv);  // INPUT-EXACT parallel HPoints
    std::vector<char> haveH(nv, 0);
    const bool ix = IXEnabled();
    auto buildH = [&](int k) -> bool {
      if (haveH[k]) return true;
      const auto it = A.provOf.find({canon3[k].x, canon3[k].y, canon3[k].z});
      if (it == A.provOf.end()) return false;
      const BuildArrangement::VProv& pv = it->second;
      if (pv.trivial) {
        eh[k] = sos::ETrivialHPoint(pv.v);
        hp[k] = sos::TrivialHPoint(pv.v);
        if (ix) bhp[k] = sos::TrivialBigHPoint(pv.v);
      } else if (pv.segPlane) {
        const vec3 ng = A.faceN[pv.rg];
        const double dg = la::dot(ng, A.tri[pv.rg][0]);
        eh[k] = sos::ESegPlaneHPoint(pv.su, pv.sw, ng, dg);
        hp[k] = sos::SegPlaneHPoint(pv.su, pv.sw, ng, dg);
        if (ix)
          bhp[k] = sos::SegPlaneBigHPoint(pv.su, pv.sw, A.tri[pv.rg].data());
      } else {
        eh[k] = sos::ECramerHPoint(pv.nF, pv.dF, pv.nG, pv.dG, pv.nH, pv.dH);
        hp[k] = sos::CramerHPoint(pv.nF, pv.dF, pv.nG, pv.dG, pv.nH, pv.dH);
        if (ix) bhp[k] = IXProvHPoint(A, pv);
      }
      haveH[k] = 1;
      return true;
    };
    // Input-exact decision when IX on (the degree-20-in-inputs orient2d, exact
    // - no rounded-basis filter, which would certify a rounded-basis sign that
    // can DISAGREE with the input-exact truth in the near-tangent wedge); else
    // the landed rounded-plane filter/exact.
    auto orientH = [&](int i, int a, int b) -> int {
      if (ix) return IXOrient2D(bhp[i], bhp[a], bhp[b], pf.axis);
      int s = sos::HomogOrient2DFilter(eh[i], eh[a], eh[b], pf.axis);
      if (s == 0) s = sos::HomogOrient2DExact(hp[i], hp[a], hp[b], pf.axis);
      return s;
    };
    // Exact proper crossing of two symbolic segments (the landed degree-9
    // straddle): four HomogOrient2D signs.  signMul cancels in the relative
    // o1!=o2 / o3!=o4 tests, so raw signs suffice.
    auto properCrossH = [&](int a0, int a1, int b0, int b1) -> bool {
      const int o1 = orientH(a0, a1, b0), o2 = orientH(a0, a1, b1);
      const int o3 = orientH(b0, b1, a0), o4 = orientH(b0, b1, a1);
      return o1 && o2 && o3 && o4 && o1 != o2 && o3 != o4;
    };
    if (!symOff) {
      bool allH = true;
      for (const EdgeM& e : edges)
        if (!buildH(find(e.v0)) || !buildH(find(e.v1))) {
          allH = false;
          break;
        }
      if (allH && (!buildH(iA) || !buildH(iB) || !buildH(iC))) allH = false;
      if (!allH) {  // a wedge arrangement vertex with no symbolic point
        if (kF4BDump) ++gF4BSeam.b3faces;
        if (std::getenv("EX2_SKIPB3") != nullptr) return;
        ok = false;
        return;
      }
      if (std::getenv("SYM_VALIDATE") !=
          nullptr) {  // provenance sanity (debug)
        double worst = 0.0;
        int worstK = -1;
        for (int k = 0; k < nv; ++k) {
          if (!haveH[k]) continue;
          const auto it =
              A.provOf.find({canon3[k].x, canon3[k].y, canon3[k].z});
          if (it == A.provOf.end() || it->second.trivial || it->second.segPlane)
            continue;
          const BuildArrangement::VProv& pv = it->second;
          auto rel = [&](const vec3& n, double d) {
            const double s = std::abs(la::dot(n, canon3[k]) - d);
            const double sc = la::length(n) * (1.0 + la::length(canon3[k]));
            return sc > 0.0 ? s / sc : s;
          };
          const double r = std::max(
              {rel(pv.nF, pv.dF), rel(pv.nG, pv.dG), rel(pv.nH, pv.dH)});
          if (r > worst) {
            worst = r;
            worstK = k;
          }
        }
        std::fprintf(stderr, "SYM_VALIDATE f=%d worstRelResid=%.3e atV=%d\n", f,
                     worst, worstK);
      }
    }
    // Symbolic plane triple of every welded rep vertex (the foreign-ness key
    // for the graze adjudication AND the completion loop's carriers()).  Grows
    // with verts2 as the completion appends crossing vertices.  A vertex is
    // FOREIGN to this face iff it is constructed (a plane triple) and its
    // triple EXCLUDES planeId[f] - eps-near plane f but not exactly on it (a
    // one-sided graze).
    std::vector<std::array<int, 3>> vertPlanes(verts2.size(), {-1, -1, -1});
    for (int k = 0; k < static_cast<int>(verts2.size()); ++k) {
      const auto it = A.provOf.find({canon3[k].x, canon3[k].y, canon3[k].z});
      if (it != A.provOf.end() && !it->second.trivial)
        vertPlanes[k] = {it->second.pf, it->second.pg, it->second.ph};
    }
    const int pfIdFor = A.planeId[f];
    auto foreignVert = [&](int V) -> bool {
      if (V < 0 || V >= static_cast<int>(vertPlanes.size())) return false;
      const auto& t = vertPlanes[V];
      return t[0] >= 0 && t[0] != pfIdFor && t[1] != pfIdFor && t[2] != pfIdFor;
    };
    // Build the arrangement 1-skeleton `uedges`.  Default (symbolic):
    // SNAP-ROUND each pre-split edge at every OTHER rep vertex within eps of
    // its interior - the exact-overlay analogue of RemoveOverlaps2D's
    // MergeVerts eps-merge, scoped to the alias / T-junction collapse.  It
    // absorbs the sub-eps ALIAS crossings (two edges crossing within eps of a
    // shared-but-unsplit vertex - the nomerge collapse the rounded
    // ExactSegProperCross hid but the exact symbolic check flags) to that
    // vertex, so no residual sub-eps proper crossing survives.  A genuine wedge
    // split is > kWedgeEnv eps from every vertex, so it is never snapped.
    // SYMWALK_OFF keeps the raw edges (its downstream double ExtractCells does
    // the eps-merge itself).
    if (symOff) {
      for (const EdgeM& e : edges) {
        const int a = find(e.v0), b = find(e.v1);
        if (a != b) uedges.push_back({a, b});
      }
    } else {
      // STITCH (mechanism b): snap V onto edge E ONLY if V lies on E's exact
      // CARRIER line - V's plane triple carries BOTH of E's two carrier planes
      // (the planes common to E's two constructed endpoints).  Without this a
      // genuine interior triple point of face f's OTHER seams that rounds
      // within eps of a SHARED seam's line near the crowded convergence gets
      // snapped ONTO the shared seam, wiring an interior-only vertex onto the
      // boundary - which the partner face (not on that vertex's third plane)
      // cannot carry, so the fans do not pair (the wedge-wedge divergence).
      // Degree-0 plane membership; keeps the aliases (same carrier) and the
      // on-carrier T-junctions, drops only the off-carrier near-tangent leak.
      // A trivial V (no triple) keeps the eps test; a trivial edge endpoint (no
      // carrier) keeps the eps test (fallback).  MEASURED NEGATIVE (opt-in
      // lever STITCH_SNAP_ON, off by default): unmasks f=467's genuine
      // skew-near- crossing (the wide double-phantom, a CARRIER_FAIL the snap
      // was papering over) -> b3faces=1, opens 56; and does not close the
      // crowded near-tangent cluster (the wedge-wedge divergence survives).
      const bool stitchSnap = std::getenv("STITCH_SNAP_ON") != nullptr;
      auto edgeCarrier2 = [&](int a, int b, int& p, int& q) -> bool {
        const auto& pa = vertPlanes[a];
        const auto& pb = vertPlanes[b];
        if (pa[0] < 0 || pb[0] < 0) return false;
        int cc[3], n = 0;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            if (pa[i] >= 0 && pa[i] == pb[j]) {
              bool seen = false;
              for (int t = 0; t < n; ++t)
                if (cc[t] == pa[i]) seen = true;
              if (!seen && n < 3) cc[n++] = pa[i];
            }
        if (n != 2) return false;
        p = cc[0];
        q = cc[1];
        return true;
      };
      std::set<int> repVerts;
      for (const EdgeM& e : edges) {
        repVerts.insert(find(e.v0));
        repVerts.insert(find(e.v1));
      }
      for (const EdgeM& e : edges) {
        const int a = find(e.v0), b = find(e.v1);
        if (a == b) continue;
        const vec3 P0 = canon3[a], P1 = canon3[b], d = P1 - P0;
        const double len2 = la::dot(d, d);
        int cp = -1, cq = -1;
        const bool haveCarrier = stitchSnap && edgeCarrier2(a, b, cp, cq);
        std::vector<std::pair<double, int>> mids;  // (param, interior vertex)
        if (len2 > 0.0)
          for (int V : repVerts) {
            if (V == a || V == b) continue;
            if (grazeAdj && foreignVert(V))
              continue;  // never snap a foreign graze onto a face-f carrier
            if (haveCarrier) {
              const auto& tv = vertPlanes[V];
              if (tv[0] >=
                  0) {  // constructed V: require it on E's carrier line
                const bool hasP = tv[0] == cp || tv[1] == cp || tv[2] == cp;
                const bool hasQ = tv[0] == cq || tv[1] == cq || tv[2] == cq;
                if (!hasP || !hasQ) continue;  // off-carrier near-tangent leak
              }
            }
            const double t = la::dot(canon3[V] - P0, d) / len2;
            if (t <= 1e-12 || t >= 1.0 - 1e-12) continue;  // strictly interior
            if (la::length(canon3[V] - (P0 + t * d)) > eps)
              continue;  // on edge
            mids.push_back({t, V});
          }
        std::sort(mids.begin(), mids.end());
        int prev = a;
        for (const auto& m : mids) {
          if (m.second != prev) uedges.push_back({prev, m.second});
          prev = m.second;
        }
        if (b != prev) uedges.push_back({prev, b});
      }
    }
    // ARRANGEMENT COMPLETION (symwalk): the pre-split PSLG (pushChain + snap)
    // can still carry proper crossings of its symbolic seam sub-edges - the
    // gap-free scale-free wedge cluster (pivot HALF-1) whose triple points
    // ExactSeamsCross's triangle-clip declined.  Find them exactly
    // (properCrossH), construct the once-only triple point cramer{f,g,h} of the
    // two crossing seams' distinct partner planes, and split EVERY seam at
    // EVERY crossing on it - keyed by exact PLANE MEMBERSHIP (cramer{f,g,h}
    // lies on both seams {f,g} and {f,h}), so the near-parallel split is
    // decided integrally, not by an unreliable double perpendicular distance.
    // Collect-then-split per pass (O(E^2), NOT O(splits*E^2)); one pass
    // completes straight segments, iterate to a fixpoint.
    if (!symOff) {
      // GRAZE STAGE-1 CENSUS (measurement only): every EPS-INJECTED element on
      // this wedge face.  The snap-round (above) inserts a rep vertex V as an
      // interior split of an edge whenever V is within eps of the edge in the
      // ROUNDED 2D frame - injection-by-eps.  Adjudicate each injection by
      // EXACTNESS instead: a constructed V is a genuine vertex of THIS face's
      // in-plane PSLG iff planeId[f] is one of V's three carrier planes (V is
      // then exactly on plane f, by construction - degree-0, no predicate); a
      // trivial input V is on plane f iff it is a corner of f or an exact
      // registry T-junction on a face-f edge.  A constructed V whose triple
      // EXCLUDES planeId[f] is FOREIGN - eps-near but not exactly on plane f, a
      // one-sided graze.  Reported here with: dPlaneF (exact perpendicular
      // offset of V's 3D point to plane f, eps units - a genuine graze is > 0;
      // an exact 0 would be a 4-plane concurrency that IS on plane f and must
      // be KEPT), the near edge + rounded on-edge distance, and projCollinear
      // (the landed degree-9 HomogOrient2DExact on the edge endpoints + V: does
      // V PROJECT onto the edge line exactly, the 2D-projection incidence the
      // snap-round used).  The verdict FOREIGN(off-plane) vs GENUINE(on-plane)
      // is the drop decision the adjudication wires.
      if (std::getenv("GRAZE_CENSUS") != nullptr) {
        const int pfId = A.planeId[f];
        const vec3 nf = A.faceN[f];
        const double nfl = la::length(nf);
        const double d0f = la::dot(nf, A.tri[f][0]);
        std::map<int, int> deg;  // rep vertex -> #incident edges
        for (const EdgeM& e : edges) {
          const int a = find(e.v0), b = find(e.v1);
          if (a == b) continue;
          ++deg[a];
          ++deg[b];
        }
        int nV = 0, nForeign = 0, nGenuine = 0, nTriv = 0, nConcurrent = 0;
        for (const auto& kv : deg) {
          const int V = kv.first;
          ++nV;
          const auto& tv = vertPlanes[V];
          const bool triv = tv[0] < 0;
          const bool onF =
              !triv && (tv[0] == pfId || tv[1] == pfId || tv[2] == pfId);
          if (triv) {
            ++nTriv;
            continue;
          }
          if (onF) {
            ++nGenuine;
            continue;
          }
          const double dPlaneF =
              std::abs(la::dot(nf, canon3[V]) - d0f) / nfl / eps;
          if (dPlaneF <= 1.0) ++nConcurrent;  // near-4-concurrency: flag/KEEP
          ++nForeign;
          // LANDED degree-9 orient2d of V against each triangle EDGE (drop the
          // face axis): does V PROJECT exactly onto a face-f carrier line?  A
          // genuine on-plane 4-concurrency (V=15) is exactly collinear on its
          // edge (== 0); a true off-plane graze projects OFF the line (!= 0).
          // This is the decisive Stage-1 separation: can the two landed
          // instantiations tell keep from drop, with NO nested construction?
          buildH(iA);
          buildH(iB);
          buildH(iC);
          const int oAB =
              sos::HomogOrient2DExact(hp[iA], hp[iB], hp[V], pf.axis);
          const int oBC =
              sos::HomogOrient2DExact(hp[iB], hp[iC], hp[V], pf.axis);
          const int oCA =
              sos::HomogOrient2DExact(hp[iC], hp[iA], hp[V], pf.axis);
          const int nZero = (oAB == 0) + (oBC == 0) + (oCA == 0);
          // Is there a GENUINE provenance alias within eps (a provOf entry
          // whose triple carries planeId[f])?  If so, this vertex IS a genuine
          // face-f point (seam endpoint / triple) that the weld represented
          // with a FOREIGN rep - a provenance-aliasing miss, not a graze.
          int aliasGen = -1;
          double aliasD = 1e300;
          for (const auto& pv : A.provOf) {
            if (pv.second.trivial) continue;
            if (pv.second.pf != pfId && pv.second.pg != pfId &&
                pv.second.ph != pfId)
              continue;
            const vec3 q{std::get<0>(pv.first), std::get<1>(pv.first),
                         std::get<2>(pv.first)};
            const double dd = la::length(q - canon3[V]);
            if (dd <= eps && dd < aliasD) {
              aliasD = dd;
              aliasGen = pv.second.pf == pfId   ? pv.second.pf
                         : pv.second.pg == pfId ? pv.second.pg
                                                : pv.second.ph;
            }
          }
          std::fprintf(stderr,
                       "GRAZE_FOREIGN f=%d V=%d planes=[%d,%d,%d] deg=%d"
                       " dPlaneF=%.3feps o2dEdges=[%d,%d,%d] onEdgeLine=%d"
                       " aliasGenuine=%d aliasD=%.3feps\n",
                       f, V, tv[0], tv[1], tv[2], kv.second, dPlaneF, oAB, oBC,
                       oCA, nZero, aliasGen,
                       aliasGen >= 0 ? aliasD / eps : -1.0);
        }
        std::fprintf(stderr,
                     "GRAZE_CENSUS f=%d PSLGverts=%d trivial=%d genuine=%d"
                     " foreign=%d nearConcurrent=%d\n",
                     f, nV, nTriv, nGenuine, nForeign, nConcurrent);
      }
      // EXJUNCT STAGE-1 CENSUS (measurement only): for every CROSS-FACE
      // arrangement vertex (a foreign triple whose planes exclude planeId[f]),
      // measure the two candidate exact-junction re-derivations onto face f's
      // carrier structure: (1) dLine = the perpendicular distance to the
      // nearest face-f carrier LINE (the projection floor - the minimum move to
      // put it on any carrier); (2) the best PLANE-TRIPLE re-derivation move:
      // min over the three pairs {a,b} of its triple, of |V -
      // cramer{a,b,planeId[f]}| (the exact seam-pair-pierces-face-plane point
      // the task prescribes).  A move
      // <= eps is close-able (a genuine near-coincidence); a move > eps is a
      // FINDING (a near-tangent grazing / skew near-crossing with no exact
      // plane-triple within eps).
      if (std::getenv("EXJ_CENSUS") != nullptr) {
        const int pf_id = A.planeId[f];
        auto lineDist = [&](const vec3& P, const vec3& u, const vec3& w) {
          const vec3 d = w - u;
          const double l2 = la::dot(d, d);
          if (!(l2 > 0.0)) return la::length(P - u);
          const double t = la::dot(P - u, d) / l2;
          return la::length(P - (u + t * d));
        };
        for (int k = 0; k < static_cast<int>(verts2.size()); ++k) {
          if (vertPlanes[k][0] < 0) continue;  // trivial / no triple
          if (vertPlanes[k][0] == pf_id || vertPlanes[k][1] == pf_id ||
              vertPlanes[k][2] == pf_id)
            continue;  // face-f's own triple (on-plane already)
          const vec3 P = canon3[k];
          double dLine = 1e300;
          for (int e = 0; e < 3; ++e)
            dLine = std::min(dLine,
                             lineDist(P, A.tri[f][e], A.tri[f][(e + 1) % 3]));
          for (int s = 0; s < (int)A.faceSeams[f].size(); ++s)
            dLine = std::min(
                dLine, lineDist(P, A.faceSeams[f][s].p0, A.faceSeams[f][s].p1));
          // model A: {seam-pair, planeId[f]} - seam pierces the face plane.
          double bestMove = 1e300;
          const auto Pf = A.planeTab.find(pf_id);
          if (Pf != A.planeTab.end())
            for (int i = 0; i < 3; ++i)
              for (int j = i + 1; j < 3; ++j) {
                const auto Pa = A.planeTab.find(vertPlanes[k][i]);
                const auto Pb = A.planeTab.find(vertPlanes[k][j]);
                if (Pa == A.planeTab.end() || Pb == A.planeTab.end()) continue;
                vec3 X;
                if (!Intersect3Planes(Pa->second.first, Pa->second.second,
                                      Pb->second.first, Pb->second.second,
                                      Pf->second.first, Pf->second.second, X))
                  continue;
                bestMove = std::min(bestMove, la::length(P - X));
              }
          // model B: nearest face-f carrier LINE {p,q} INT one foreign plane -
          // where a foreign plane crosses face f's own edge/seam line.  Gather
          // the {p,q} of every carrier line (triangle edges {f,nb}, seams
          // {f,other}), pick the line nearest V, re-derive {p,q,X} for X in V's
          // triple, min move.
          double lineMove = 1e300;
          auto tryLine = [&](int p, int q, const vec3& u, const vec3& w) {
            if (lineDist(P, u, w) > 8.0 * eps) return;
            const auto Pp = A.planeTab.find(p), Pq = A.planeTab.find(q);
            if (Pp == A.planeTab.end() || Pq == A.planeTab.end()) return;
            for (int t = 0; t < 3; ++t) {
              const auto Px = A.planeTab.find(vertPlanes[k][t]);
              if (Px == A.planeTab.end()) continue;
              vec3 X;
              if (!Intersect3Planes(Pp->second.first, Pp->second.second,
                                    Pq->second.first, Pq->second.second,
                                    Px->second.first, Px->second.second, X))
                continue;
              lineMove = std::min(lineMove, la::length(P - X));
            }
          };
          for (int e = 0; e < 3; ++e) {
            const int nb = static_cast<int>(in.halfedge_.Pair(3 * f + e)) / 3;
            if (nb >= 0 && nb < static_cast<int>(A.planeId.size()))
              tryLine(pf_id, A.planeId[nb], A.tri[f][e], A.tri[f][(e + 1) % 3]);
          }
          for (int s = 0; s < (int)A.faceSeams[f].size(); ++s) {
            const int ot = A.faceSeams[f][s].other;
            if (ot >= 0 && ot < static_cast<int>(A.planeId.size()))
              tryLine(pf_id, A.planeId[ot], A.faceSeams[f][s].p0,
                      A.faceSeams[f][s].p1);
          }
          std::fprintf(stderr,
                       "EXJ_CENSUS f=%d v=%d planes=[%d,%d,%d] dLine=%.3feps "
                       "tripleMove=%.3feps lineMove=%.3feps\n",
                       f, k, vertPlanes[k][0], vertPlanes[k][1],
                       vertPlanes[k][2], dLine / eps,
                       bestMove > 1e299 ? -1.0 : bestMove / eps,
                       lineMove > 1e299 ? -1.0 : lineMove / eps);
        }
      }
      // GRAZE ADJUDICATION (arrangement-vertex level, exact-overlay only): the
      // pushChain filter drops foreign junction SPLITS, but the eps-weld can
      // still make a FOREIGN triple point the min-index REP of an arrangement
      // vertex (its bits win the union-find over a coincident junction whose
      // own bits passed the filter), leaving a residual foreign vertex -
      // constructed, its triple excluding planeId[f], no genuine within-eps
      // alias (measured aliasGenuine<0) - as a low-degree SPUR whose phantom
      // edge crosses a real carrier (the f=467 CARRIER_FAIL that survived the
      // split filter).  Such a vertex is off plane f (a one-sided graze), NOT
      // part of this face's in-plane PSLG: drop every uedge incident to it.  A
      // deg<=1 spur simply drops (ExtractCellsSymbolic would prune it anyway -
      // here we drop it BEFORE the completion/completeness check so its phantom
      // crossing never fires); a deg-2 interior foreign vertex splices its two
      // GENUINE neighbours so the carrier stays whole.  A deg>=3 foreign hub is
      // not expected once the split filter has run - fail closed (never emit an
      // ambiguous splice).  Combinatorial (degree-0) plane membership on the
      // welded rep - no new predicate, no nested construction.
      if (grazeAdj) {
        const int pfIdV = A.planeId[f];
        auto isForeign = [&](int V) -> bool {
          if (V < 0 || V >= static_cast<int>(vertPlanes.size())) return false;
          const auto& t = vertPlanes[V];
          return t[0] >= 0 && t[0] != pfIdV && t[1] != pfIdV && t[2] != pfIdV;
        };
        std::map<int, std::vector<int>> adj;
        bool anyForeign = false;
        for (const auto& e : uedges) {
          adj[e.first].push_back(e.second);
          adj[e.second].push_back(e.first);
          if (isForeign(e.first) || isForeign(e.second)) anyForeign = true;
        }
        if (anyForeign) {
          std::vector<std::pair<int, int>> nu;
          nu.reserve(uedges.size());
          for (const auto& e : uedges)  // drop every foreign-incident edge
            if (!isForeign(e.first) && !isForeign(e.second)) nu.push_back(e);
          bool hub = false;
          for (const auto& kv : adj) {
            if (!isForeign(kv.first)) continue;
            if (kv.second.size() == 2) {  // splice the carrier across it
              const int a = kv.second[0], b = kv.second[1];
              if (a != b && !isForeign(a) && !isForeign(b))
                nu.push_back({a, b});
            } else if (kv.second.size() >= 3) {
              hub = true;
            }
          }
          if (std::getenv("GRAZE_REMOVE") != nullptr) {
            int nf2 = 0, maxdeg = 0;
            for (const auto& kv : adj)
              if (isForeign(kv.first)) {
                ++nf2;
                maxdeg = std::max(maxdeg, (int)kv.second.size());
              }
            std::fprintf(stderr,
                         "GRAZE_REMOVE f=%d foreignVerts=%d maxForeignDeg=%d"
                         " hub=%d uedges %d->%d\n",
                         f, nf2, maxdeg, hub ? 1 : 0,
                         static_cast<int>(uedges.size()),
                         static_cast<int>(nu.size()));
          }
          if (hub) {  // a foreign hub survived the split filter: fail closed
            if (kF4BDump) ++gF4BSeam.b3faces;
            if (std::getenv("EX2_SKIPB3") != nullptr) return;
            ok = false;
            return;
          }
          uedges.swap(nu);
        }
      }
      auto carriers = [&](int u, int v, int& p, int& q) -> bool {
        const auto &pu = vertPlanes[u], &pv = vertPlanes[v];
        if (pu[0] < 0 || pv[0] < 0) return false;
        int cc[3], n = 0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            if (pu[a] == pv[b]) {
              bool seen = false;
              for (int t = 0; t < n; ++t)
                if (cc[t] == pu[a]) seen = true;
              if (!seen && n < 3) cc[n++] = pu[a];
            }
        if (n != 2) return false;
        p = cc[0];
        q = cc[1];
        return true;
      };
      auto hasPlanes = [&](int V, int p, int q) -> bool {
        const auto& pv = vertPlanes[V];
        bool a = false, b = false;
        for (int t = 0; t < 3; ++t) {
          if (pv[t] == p) a = true;
          if (pv[t] == q) b = true;
        }
        return a && b;
      };
      bool completed = true;
      for (int pass = 0; pass < 6 && completed; ++pass) {
        std::vector<int> crossVerts;  // newly-added crossing vertices this pass
        const int nu0 = static_cast<int>(uedges.size());
        for (int i = 0; i < nu0 && completed; ++i)
          for (int j = i + 1; j < nu0 && completed; ++j) {
            const int a0 = uedges[i].first, a1 = uedges[i].second;
            const int b0 = uedges[j].first, b1 = uedges[j].second;
            if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) continue;
            if (a0 == a1 || b0 == b1) continue;
            if (!properCrossH(a0, a1, b0, b1)) continue;
            int ca[2], cb[2];
            if (!carriers(a0, a1, ca[0], ca[1]) ||
                !carriers(b0, b1, cb[0], cb[1])) {
              if (std::getenv("SYM_CDIAG") != nullptr)
                std::fprintf(stderr,
                             "SYM_CDIAG f=%d CARRIER_FAIL e0=(%d[%d,%d,%d],"
                             "%d[%d,%d,%d]) e1=(%d[%d,%d,%d],%d[%d,%d,%d])\n",
                             f, a0, vertPlanes[a0][0], vertPlanes[a0][1],
                             vertPlanes[a0][2], a1, vertPlanes[a1][0],
                             vertPlanes[a1][1], vertPlanes[a1][2], b0,
                             vertPlanes[b0][0], vertPlanes[b0][1],
                             vertPlanes[b0][2], b1, vertPlanes[b1][0],
                             vertPlanes[b1][1], vertPlanes[b1][2]);
              if (std::getenv("SYM_GEOM") != nullptr) {
                const vec3 nf = A.faceN[f];
                const double nfl = la::length(nf);
                const double d0f = la::dot(nf, A.tri[f][0]);
                auto lineDist = [&](const vec3& P, const vec3& u,
                                    const vec3& w) {
                  const vec3 d = w - u;
                  const double l2 = la::dot(d, d);
                  if (!(l2 > 0.0)) return la::length(P - u);
                  const double t = la::dot(P - u, d) / l2;
                  return la::length(P - (u + t * d));
                };
                for (int vv : {a0, a1, b0, b1}) {
                  const vec3 P = canon3[vv];
                  const double dP = std::abs(la::dot(nf, P) - d0f) / nfl / eps;
                  // nearest face-f carrier line (edges + seams)
                  double best = 1e300;
                  int bk = -2;  // -1..-3 edges (by e), else seam k
                  for (int e = 0; e < 3; ++e) {
                    const double dd =
                        lineDist(P, A.tri[f][e], A.tri[f][(e + 1) % 3]);
                    if (dd < best) {
                      best = dd;
                      bk = -1 - e;
                    }
                  }
                  for (int k = 0; k < (int)A.faceSeams[f].size(); ++k) {
                    const double dd =
                        lineDist(P, A.faceSeams[f][k].p0, A.faceSeams[f][k].p1);
                    if (dd < best) {
                      best = dd;
                      bk = k;
                    }
                  }
                  // param along the near line + full coords
                  const vec3 lu =
                      bk < 0 ? A.tri[f][-bk - 1] : A.faceSeams[f][bk].p0;
                  const vec3 lw = bk < 0 ? A.tri[f][(-bk - 1 + 1) % 3]
                                         : A.faceSeams[f][bk].p1;
                  const vec3 ld = lw - lu;
                  const double tt = la::dot(P - lu, ld) / la::dot(ld, ld);
                  std::fprintf(
                      stderr,
                      "SYM_GEOM f=%d v=%d dPlane=%.2feps nearLine=%s%d "
                      "dLine=%.2feps t=%.6f P=(%.17g,%.17g,%.17g)\n",
                      f, vv, dP, bk < 0 ? "edge" : "seam",
                      bk < 0 ? -bk - 1 : bk, best / eps, tt, P.x, P.y, P.z);
                }
                // mutual 3D distances among the 4 endpoints
                const int vs[4] = {a0, a1, b0, b1};
                for (int x = 0; x < 4; ++x)
                  for (int y = x + 1; y < 4; ++y)
                    std::fprintf(
                        stderr, "SYM_GEOM f=%d dist v%d-v%d=%.3feps\n", f,
                        vs[x], vs[y],
                        la::length(canon3[vs[x]] - canon3[vs[y]]) / eps);
              }
              completed = false;  // a non-seam edge crossing: cannot complete
              break;
            }
            int s = -1;
            for (int x = 0; x < 2; ++x)
              for (int y = 0; y < 2; ++y)
                if (ca[x] == cb[y]) s = ca[x];
            const int g = (ca[0] == s) ? ca[1] : ca[0];
            const int h = (cb[0] == s) ? cb[1] : cb[0];
            const auto Ps = A.planeTab.find(s), Pg = A.planeTab.find(g),
                       Ph = A.planeTab.find(h);
            if (s < 0 || g == h || g == s || h == s || Ps == A.planeTab.end() ||
                Pg == A.planeTab.end() || Ph == A.planeTab.end()) {
              completed = false;
              break;
            }
            vec3 X;
            if (!Intersect3Planes(Ps->second.first, Ps->second.second,
                                  Pg->second.first, Pg->second.second,
                                  Ph->second.first, Ph->second.second, X)) {
              completed = false;
              break;
            }
            // INPUT-EXACT completion-crossing position (the once-only {s,g,h}
            // triple composed exactly from the rep triangles): correctly
            // ordered vs the rounded X (off by up to ~1400x eps in the wedge).
            sos::BigHPoint bX;
            bool haveBX = false;
            if (ix) {
              const auto Ts = A.planeTri.find(s), Tg = A.planeTri.find(g),
                         Th = A.planeTri.find(h);
              if (Ts != A.planeTri.end() && Tg != A.planeTri.end() &&
                  Th != A.planeTri.end()) {
                bX = IXTripleHPoint(A, Ts->second, Tg->second, Th->second);
                const vec3 Xix = sos::BigHPointToPos(bX);
                if (std::isfinite(Xix.x) && std::isfinite(Xix.y) &&
                    std::isfinite(Xix.z)) {
                  X = Xix;
                  haveBX = true;
                }
              }
            }
            bool exists =
                false;  // dedup: already present (this pass's split uses it)
            for (int k = 0; k < static_cast<int>(canon3.size()); ++k)
              if (la::length(canon3[k] - X) <= eps) {
                exists = true;
                break;
              }
            if (exists) continue;
            const vec3 nS = Ps->second.first, nG = Pg->second.first,
                       nH = Ph->second.first;
            const double dS = la::dot(nS, Ps->second.second),
                         dG = la::dot(nG, Pg->second.second),
                         dH = la::dot(nH, Ph->second.second);
            crossVerts.push_back(static_cast<int>(verts2.size()));
            verts2.push_back(pf.proj(X));
            canon3.push_back(X);
            eh.push_back(sos::ECramerHPoint(nS, dS, nG, dG, nH, dH));
            hp.push_back(sos::CramerHPoint(nS, dS, nG, dG, nH, dH));
            bhp.push_back(haveBX ? bX : sos::BigHPoint{});
            haveH.push_back(1);
            vertPlanes.push_back({s, g, h});
          }
        if (!completed) break;
        if (crossVerts.empty()) break;  // no crossing left: complete
        // Re-split every seam sub-edge at every crossing vertex on its carrier
        // line (exact plane membership), ordered by param.
        std::vector<std::pair<int, int>> nu;
        nu.reserve(uedges.size() + crossVerts.size() * 2);
        for (const auto& e : uedges) {
          const int a = e.first, b = e.second;
          int p, q;
          if (!carriers(a, b, p, q)) {
            nu.push_back(e);
            continue;
          }
          const vec3 P0 = canon3[a], P1 = canon3[b], d = P1 - P0;
          const double len2 = la::dot(d, d);
          std::vector<std::pair<double, int>> mids;
          if (len2 > 0.0)
            for (int V : crossVerts) {
              if (V == a || V == b || !hasPlanes(V, p, q)) continue;
              const double t = la::dot(canon3[V] - P0, d) / len2;
              if (t <= 1e-12 || t >= 1.0 - 1e-12) continue;
              mids.push_back({t, V});
            }
          std::sort(mids.begin(), mids.end());
          int prev = a;
          for (const auto& m : mids) {
            if (m.second != prev) nu.push_back({prev, m.second});
            prev = m.second;
          }
          if (b != prev) nu.push_back({prev, b});
        }
        uedges.swap(nu);
      }
      if (std::getenv("SYM_CDIAG") != nullptr)
        std::fprintf(stderr, "SYM_CDIAG f=%d nEdges=%d completed=%d\n", f,
                     static_cast<int>(uedges.size()), completed ? 1 : 0);
      if (!completed) {  // could not build a crossing-free arrangement
        if (kF4BDump) ++gF4BSeam.b3faces;
        if (std::getenv("EX2_SKIPB3") != nullptr) return;
        ok = false;
        return;
      }
    }
    // COMPLETENESS: no two non-adjacent 1-skeleton edges properly cross in the
    // face interior (a residual crossing = an unbuilt vertex -> fail closed).
    // EXACT on the symbolic HPoints (default): the f=467 phantom (ROUNDED seams
    // cross where the EXACT seams do NOT) verifies CLEAN, so the symbolic walk
    // simply builds the crossing-free arrangement.  SYMWALK_OFF uses the
    // rounded-canon3 ExactSegProperCross (which false-positives that phantom).
    const int nu = static_cast<int>(uedges.size());
    for (int i = 0; i < nu; ++i)
      for (int j = i + 1; j < nu; ++j) {
        const int a0 = uedges[i].first, a1 = uedges[i].second;
        const int b0 = uedges[j].first, b1 = uedges[j].second;
        if (a0 == b0 || a0 == b1 || a1 == b0 || a1 == b1) continue;  // adjacent
        if (a0 == a1 || b0 == b1) continue;  // degenerate
        const bool cross =
            symOff ? ExactSegProperCross(canon3[a0], canon3[a1], canon3[b0],
                                         canon3[b1], pf.axis)
                   : properCrossH(a0, a1, b0, b1);
        if (cross) {
          if (std::getenv("EX2_DIAG") != nullptr) {
            const vec2 xc = SegLineIntersect2D(verts2[a0], verts2[a1],
                                               verts2[b0], verts2[b1]);
            double dmin = 1e300;
            for (const vec2& q : verts2)
              dmin = std::min(dmin, la::length(q - xc));
            auto kind = [&](int k) -> char {
              const auto it =
                  A.provOf.find({canon3[k].x, canon3[k].y, canon3[k].z});
              return (it == A.provOf.end()) ? '?'
                     : it->second.trivial   ? 'T'
                                            : 'C';
            };
            std::fprintf(stderr,
                         "EX2_DIAG f=%d residualCross dNearestVert=%.3feps "
                         "e0=(%d%c,%d%c) e1=(%d%c,%d%c) len0=%.2feps "
                         "len1=%.2feps\n",
                         f, dmin / eps, a0, kind(a0), a1, kind(a1), b0,
                         kind(b0), b1, kind(b1),
                         la::length(canon3[a0] - canon3[a1]) / eps,
                         la::length(canon3[b0] - canon3[b1]) / eps);
          }
          if (kF4BDump) ++gF4BSeam.b3faces;
          // Measurement lever: emit nothing but do NOT fail the component, so
          // the open-edge count of the RESOLVED faces past this one is
          // measurable.
          if (std::getenv("EX2_SKIPB3") != nullptr) return;
          ok = false;
          return;
        }
      }
    // SYMBOLIC CELL WALK: order the cell walk by the exact HomogOrient2D on the
    // per-vertex HPoints (ExtractCellsSymbolic), realizing the exact
    // arrangement regardless of the rounded straddle - unlike ExtractCells'
    // rounded-double atan2, which folds back on the near-parallel wedges (the
    // 56-open residue). SYMWALK_OFF keeps the double ExtractCells (the mutation
    // anchor).
    if (!symOff) {
      // signMul aligns HomogOrient2D's keep-pair frame to verts2's CCW proj
      // frame using the triangle corners (verts2-CCW by the earlier assert).
      const int sc =
          ix ? IXOrient2D(bhp[iA], bhp[iB], bhp[iC], pf.axis)
             : sos::HomogOrient2DExact(hp[iA], hp[iB], hp[iC], pf.axis);
      const int signMul = sc >= 0 ? 1 : -1;
      if (!ExtractCellsSymbolic(verts2, eh, hp, bhp, ix, pf.axis, signMul,
                                uedges, cells)) {
        if (kF4BDump) ++gF4BSeam.b4faces;
        if (std::getenv("EX2_SKIPB3") != nullptr) return;
        ok = false;
        return;
      }
      cellsBuilt = true;
    }
  } else {
    // DEFAULT (RemoveOverlaps2D): edgeSubdiv gives the exact per-input-edge
    // subdivision (triangle edges split at the on-edge seam endpoints; seams
    // split at any crossing).  The winding rule is irrelevant here (we take
    // only the subdivision), so any pred works.  eps is the component weld
    // radius; the canonical seam points are exact-shared, so no extra
    // construction headroom is needed and a tight eps avoids merging
    // genuinely-distinct sub-eps features.
    std::vector<std::vector<vec2>> sub;
    RemoveOverlaps2D(verts2, edges, eps, /*debug=*/false, WindRule::Add,
                     /*trace=*/nullptr, /*edgesNeg=*/nullptr, &sub);

    // Reconstruct the planar subdivision, mapping each arrangement position
    // back to its input vertex bit-exactly.  A position with no input preimage
    // is a NEW crossing vertex (a >2-sheet triple point, unbuilt) -> fail
    // closed.
    std::map<std::tuple<double, double>, int> pos2in;
    for (int k = 0; k < static_cast<int>(verts2.size()); ++k)
      pos2in.emplace(std::tuple<double, double>{verts2[k].x, verts2[k].y}, k);
    for (const auto& poly : sub) {
      for (size_t k = 0; k + 1 < poly.size(); ++k) {
        const auto i0 = pos2in.find({poly[k].x, poly[k].y});
        const auto i1 = pos2in.find({poly[k + 1].x, poly[k + 1].y});
        if (i0 == pos2in.end() || i1 == pos2in.end()) {
          if (kF4BDump) {  // census: count every residual constructed vertex
            std::set<std::tuple<double, double>> miss;
            for (const auto& p : sub)
              for (const vec2& q : p)
                if (pos2in.find({q.x, q.y}) == pos2in.end())
                  miss.insert({q.x, q.y});
            ++gF4BSeam.b3faces;
            gF4BSeam.b3verts += static_cast<int>(miss.size());
          }
          ok = false;
          return;
        }
        uedges.push_back({i0->second, i1->second});
      }
    }
  }

  // STITCH CENSUS (measurement only): the boundary subdivision this face places
  // on each of its shared boundary segments (3 triangle edges + each seam), as
  // the FINAL uedges realize it.  Emitted per-face; grouped OFFLINE by the
  // shared segment (sorted endpoint bits) to diff the two incident faces' split
  // lists - the two-sided-consistency divergence.  inReg = the split is within
  // eps of a global registry junction (a canonical shared vertex both sides
  // consume); inReg=0 = exact-only structure (a wedge completion crossing /
  // weld rep no neighbour has).  Env-gated, byte-clean off STITCH_CENSUS.
  if (std::getenv("STITCH_CENSUS") != nullptr) {
    std::set<int> used;
    for (const auto& e : uedges) {
      used.insert(e.first);
      used.insert(e.second);
    }
    auto regNear = [&](const vec3& P) -> int {
      for (const vec3& j : A.junctions)
        if (la::length(j - P) <= eps) return 1;
      return 0;
    };
    auto onLine = [&](int k, const vec3& P0, const vec3& d, double l2,
                      double& t) -> bool {
      if (l2 <= 0.0) return false;
      t = la::dot(canon3[k] - P0, d) / l2;
      if (t <= 1e-9 || t >= 1.0 - 1e-9) return false;
      return la::length(canon3[k] - (P0 + t * d)) <= eps;
    };
    auto dumpSeg = [&](const char* kind, int sub, const vec3& P0,
                       const vec3& P1) {
      const vec3 d = P1 - P0;
      const double l2 = la::dot(d, d);
      // The ACTUAL subdivision: interior vertices that are an endpoint of a
      // uedge whose OTHER endpoint is ALSO on this segment line - i.e. a vertex
      // the cell walk connects ALONG this seam (not a foreign vertex merely
      // eps-near the line but wired into a different seam's chain).
      std::set<int> onSeamEdge;
      for (const auto& e : uedges) {
        double ta, tb;
        const bool a = onLine(e.first, P0, d, l2, ta) ||
                       la::length(canon3[e.first] - P0) <= eps ||
                       la::length(canon3[e.first] - P1) <= eps;
        const bool b = onLine(e.second, P0, d, l2, tb) ||
                       la::length(canon3[e.second] - P0) <= eps ||
                       la::length(canon3[e.second] - P1) <= eps;
        if (a && b) {
          if (onLine(e.first, P0, d, l2, ta)) onSeamEdge.insert(e.first);
          if (onLine(e.second, P0, d, l2, tb)) onSeamEdge.insert(e.second);
        }
      }
      std::vector<std::pair<double, int>> on;
      for (int k : onSeamEdge) {
        double t;
        if (onLine(k, P0, d, l2, t)) on.push_back({t, k});
      }
      std::sort(on.begin(), on.end());
      std::fprintf(stderr,
                   "STITCH_SEG f=%d wedge=%d %s%d nsplit=%d "
                   "e0=(%.17g,%.17g,%.17g) e1=(%.17g,%.17g,%.17g)\n",
                   f, wedge ? 1 : 0, kind, sub, static_cast<int>(on.size()),
                   P0.x, P0.y, P0.z, P1.x, P1.y, P1.z);
      for (const auto& m : on) {
        const vec3& P = canon3[m.second];
        int tp0 = -1, tp1 = -1, tp2 = -1, triv = 1;
        const auto it = A.provOf.find({P.x, P.y, P.z});
        if (it != A.provOf.end() && !it->second.trivial) {
          triv = 0;
          tp0 = it->second.pf;
          tp1 = it->second.pg;
          tp2 = it->second.ph;
        }
        std::fprintf(stderr,
                     "STITCH_SPLIT f=%d %s%d t=%.9f inReg=%d triv=%d "
                     "planes=[%d,%d,%d] P=(%.17g,%.17g,%.17g)\n",
                     f, kind, sub, m.first, regNear(P), triv, tp0, tp1, tp2,
                     P.x, P.y, P.z);
      }
    };
    dumpSeg("edge", 0, a, b);
    dumpSeg("edge", 1, b, c);
    dumpSeg("edge", 2, c, a);
    for (int k = 0; k < static_cast<int>(A.faceSeams[f].size()); ++k)
      dumpSeg("seam", k, A.faceSeams[f][k].p0, A.faceSeams[f][k].p1);
  }

  if (!cellsBuilt && !ExtractCells(verts2, uedges, cells)) {
    if (kF4BDump) ++gF4BSeam.b4faces;
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
    double area2;
    const int best = LargestSubTri(tris, verts2, area2);
    // Interior classify point = the DIRECT 3D centroid of the largest sub-
    // triangle's canonical points (projection-independent - the sub-triangle
    // verts are all shared 3D constructions), not a back-projection of the 2D
    // centroid.  Strictly interior to the cell for the winding probe.
    const vec3 cen3 =
        (canon3[tris[best].x] + canon3[tris[best].y] + canon3[tris[best].z]) /
        3.0;
    // BOTH-SIDES RETENTION (f4-junction, port of the coplanar fold's rule):
    // probe w_S on BOTH sides of the sub-cell and retain iff EXACTLY ONE side
    // is inside {w_S>=1}; the solid side fixes the emitted orientation.  This
    // is the general d{w_S>=1} boundary criterion, not the mult-1
    // specialization.  For an oriented mult-1 face w_below = w_above + 1 (the
    // +/-1 crossing delta, seed-independent), so aboveIn!=belowIn holds EXACTLY
    // at w_above==0 with the solid on the -n side (belowIn) -> original
    // orientation: BITWISE-IDENTICAL to the former w_above==0 rule on every
    // jump==1 cell (the whole corpus off openscad).  Where the true jump is NOT
    // 1 (a coplanar coincidence a transversal seam sub-cell can carry -
    // openscad's tangent/overlap fans), the one-sided rule mis-orients or
    // wrongly drops the sub-face (the F4B census's radial ties + material
    // overlaps); the both-sides read decides them correctly.  A negative
    // w_above (both sides exterior) drops, NOT fail-closed (subtraction
    // absorption).
    const std::optional<int> g = RobustWinding(in, cen3 + eps * nHat, seeds);
    const std::optional<int> gb = RobustWinding(in, cen3 - eps * nHat, seeds);
    if (!g || !gb) {  // filter-uncertain deciding predicate (SoS axis): fail
      if (kF4BDump) ++gF4BSeam.b5faces;  // closed
      ok = false;
      return;
    }
    const int keep = BothSidesRetain(*g, *gb);
    if (wedge && std::getenv("EX2_CELLDUMP") != nullptr) {
      double aC = 0.0;
      for (size_t k = 0; k < cell.size(); ++k)
        aC += la::cross(verts2[cell[k]], verts2[cell[(k + 1) % cell.size()]]);
      std::fprintf(
          stderr, "EX2_CELL f=%d n=%d area=%.3geps2 w+=%d w-=%d keep=%d\n", f,
          static_cast<int>(cell.size()), 0.5 * aC / (eps * eps), *g, *gb, keep);
    }
    if (keep == 0) continue;  // both sides same class: not a boundary
    // Retained: keep==1 (solid on the -nHat side) keeps the CCW 2D winding
    // (+nHat = original orientation); keep==2 (solid on +nHat) reverses.
    for (const ivec3& t : tris) {
      if (keep == 1)
        out.push_back({canon3[t.x], canon3[t.y], canon3[t.z]});
      else
        out.push_back({canon3[t.x], canon3[t.z], canon3[t.y]});
    }
  }
  if (kF4BDump) ++gF4BSeam.okfaces;
}

// 2D point-in-triangle (inclusive), orientation-agnostic: true iff p is on the
// same side (or on) all three directed edges under either winding.
bool PointInTri2D(const vec2& p, const vec2& a, const vec2& b, const vec2& c) {
  const double d1 = la::cross(b - a, p - a), d2 = la::cross(c - b, p - b),
               d3 = la::cross(a - c, p - c);
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
// oriented sums).  Retention generalizes the resolver's mult-1 rule: with
// w_below = w_above + m (the 3D winding jump across the plane equals the
// coincident cover, so it is SELF-CHECKED against the real coupled winding on
// both sides), a sub-face is on d{w_S>=1} iff EXACTLY ONE side is inside
// {w>=1}; the solid side fixes the emitted orientation.  m==0 (pure
// cancellation) drops.  A cluster face that is ALSO transversally seamed is the
// coplanar/transversal ENTANGLEMENT (a seam line would straddle a cell); that
// stays fail-closed with its own named reason, distinct from the near-coplanar
// residue.  Any degenerate projection, malformed cell walk, filter-uncertain
// winding, or self-check mismatch fails closed (never a silent wrong resolve).
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
    const vec3 a0 = A.tri[f0][0];
    // The cluster plane -> the dominant axis of the fitted (f0) normal; a0
    // fixes the plane offset for lift.
    AxisDropFrame pf;
    if (!BuildAxisDropFrame(A.faceN[f0], a0, pf)) {
      ok = false;
      return;
    }
    const vec3& nHat = pf.nHat;
    // Input verts (dedup by canonical 3D bit pattern) + triangle-boundary
    // edges; each member triangle carries its signed orientation vs nHat.
    std::vector<vec2>& verts2 = pf.verts2;
    std::vector<vec3>& canon3 = pf.canon3;
    auto getV = [&](const vec3& P) { return pf.add(P); };
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
    // Boundary segments (net!=0), used ONLY by the retired FOLD_MEMBER_GUARD
    // mutation lever below (the member-member-crossing skip that is UNSOUND at
    // a near-degenerate cluster; see the canonical-per-line comment below).
    std::vector<std::pair<vec3, vec3>> bsegs;
    for (const auto& [e, net] : dir)
      if (net != 0) bsegs.push_back({canon3[e.first], canon3[e.second]});
    auto onBoundarySeg = [&](const vec3& V, const vec3& S0, const vec3& S1) {
      const vec3 d = S1 - S0;
      const double l2 = la::dot(d, d);
      if (!(l2 > 0.0)) return false;
      const double margin = eps / std::sqrt(l2);
      const double t = la::dot(V - S0, d) / l2;
      if (t < -margin || t > 1.0 + margin) return false;  // inclusive endpoints
      return la::length((V - S0) - t * d) <= eps;
    };
    std::vector<EdgeM> segEdges;
    // CANONICAL PER-LINE SUBDIVISION (perline): split this fold-boundary edge
    // at EVERY registry junction strictly interior to it, so the fold's caps
    // and the incident walls subdivide the shared boundary line at the
    // IDENTICAL once-only registry vertex and their emission fans pair.
    //
    // The prior f4-junction guard skipped a junction lying on >=2 member
    // boundary segments ("a member-member crossing RemoveOverlaps2D owns it").
    // That is UNSOUND at a near-degenerate coplanar cluster: RO2D re-derives
    // the member-member crossing per-face in ROUNDED doubles, landing it
    // sub-eps off the CANONICAL registry junction the incident wall split at,
    // so the cap sub-face spanned the shared line unsplit and its fan opened
    // against the wall's (openscad's cap-wall T-junction opens).  Splitting at
    // the canonical junction unconditionally is once-only-correct: the registry
    // vertex is the shared 3D construction both the cap and the wall reference,
    // so presplitting there never disagrees with a neighbour.  A cap-cap
    // crossing RO2D would also split at is now pre-split at the SAME canonical
    // point instead of RO2D's rounded recomputation (an equivalence-preserving
    // re-triangulation on the rotated/irrational-junction fold carriers,
    // byte-identical on the axis-aligned ones).  FOLD_MEMBER_GUARD restores the
    // old skip (mutation lever: the openscad cap-wall opens reappear).
    static const bool kMemberGuard =
        std::getenv("FOLD_MEMBER_GUARD") != nullptr;
    for (const auto& [e, net] : dir) {
      if (net == 0) continue;  // interior (cancelled) diagonal: not a boundary
      const std::vector<std::pair<vec2, vec3>> splits = JunctionSplitsOnSegment(
          pf.proj(canon3[e.first]), pf.proj(canon3[e.second]), canon3[e.first],
          canon3[e.second], A.junctions, eps);
      int prev = e.first;
      for (const auto& sp : splits) {
        if (kMemberGuard) {
          int onCount = 0;
          for (const auto& bs : bsegs)
            if (onBoundarySeg(sp.second, bs.first, bs.second)) ++onCount;
          if (onCount >= 2) continue;  // member-member crossing: RO2D owns it
        }
        const int v = pf.addAt(sp.first, sp.second);
        if (v != prev) segEdges.push_back({prev, v, 1});
        prev = v;
      }
      if (e.second != prev) segEdges.push_back({prev, e.second, 1});
    }

    // ARRANGEMENT via RemoveOverlaps2D (the same primitive the seamed path
    // reuses at EmitSeamedFace): the member triangles' boundary edges overlap
    // in-plane, and edgeSubdiv splits each input edge at BOTH proper crossings
    // AND T-junctions / collinear incidences (a hand-rolled pairwise-crossing
    // arrangement would miss the T-junctions, splitting one edge without
    // splitting the edge that ends on it - RO2DProbe, reg3d-s4-verify Audit
    // 2a). The winding rule is irrelevant here (we consume only the
    // subdivision). Every arrangement vertex lies in this exact plane, so a NEW
    // crossing position's 3D image is the axis-drop lift of its 2D coords (an
    // input vert keeps its canonical 3D via getP; the eps-box match folds a
    // merged endpoint back onto its input vert).
    std::vector<vec2> pts = verts2;
    std::vector<vec3> pts3 = canon3;
    auto getP = [&](const vec2& q) {
      for (int k = 0; k < static_cast<int>(pts.size()); ++k)
        if (std::abs(pts[k].x - q.x) <= eps && std::abs(pts[k].y - q.y) <= eps)
          return k;
      const int id = static_cast<int>(pts.size());
      pts.push_back(q);
      pts3.push_back(pf.lift(q));
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
        ar += la::cross(p, q);
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
      double bestArea;
      const int best = LargestSubTri(tris, pts, bestArea);
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
      // cen2 is a 2D-born interior point (no 3D preimage): lift it onto the
      // cluster plane.  Still strictly interior to the sub-triangle (the drop
      // is an affine image), so a valid winding probe.
      const vec3 cen3 = pf.lift(cen2);
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
      const int keep = BothSidesRetain(*wa, *wb);
      if (keep == 0) continue;  // both sides same class: not a boundary
      // Retained. keep==1 (solid on the -nHat side) keeps the CCW 2D winding
      // (+nHat); keep==2 (solid on +nHat) reverses.
      for (const ivec3& t : tris) {
        if (keep == 1)
          out.push_back({pts3[t.x], pts3[t.y], pts3[t.z]});
        else
          out.push_back({pts3[t.x], pts3[t.z], pts3[t.y]});
      }
    }
  }
}

// Classify + emit the CLEAN (un-seamed, un-folded) faces.  Retention rule (SAME
// witness rule as the seamed path): BOTH-SIDES retention - probe w_S on each
// side of the face and keep iff EXACTLY ONE side is inside {w_S>=1}; the solid
// side fixes orientation (solid on -n -> original, on +n -> reversed).  For an
// oriented mult-1 clean face w_below = w_above + 1, so this reduces EXACTLY to
// the former w_above==0 / original-orientation rule (bit-identical off
// openscad); it additionally drops a clean tile that is coplanar-coincident
// with an unclustered seamed face (w_above==w_below==0, exterior on both
// sides), the one-sided-rule over-emission the whole-face read cannot see.  A
// NEGATIVE w_above (both sides exterior) still DROPS (the axis-1 subtraction
// absorption), not a fail-closed.
//
// PER-FACE, not per-patch: coverage is constant across an uncrossed clean-clean
// edge only when the arrangement is complete, so a per-patch representative
// flood is silently wrong on a folded soup (an everted corner can put a
// w_above==0 boundary face and a w_above==-1 exterior face in ONE patch, so a
// w<0 representative would drop genuine boundary faces).  Probing each face by
// its OWN winding is uniformly sound - a wrong retain is impossible because
// every emitted face's winding is directly measured.  When the centroid probe
// GRAZES every seed (axis-aligned integer geometry), re-probe at other interior
// points of the SAME triangle: the winding cell just above a clean triangle is
// constant, so this dodges the graze without any patch-uniformity assumption; a
// face that grazes at EVERY interior point fails closed.  Cost: O(nTri) winding
// queries worst case (the winding-query perf axis; a BVH is the later pass).
// The everted-corner rationale is in docs/Regularize3D.md.
bool EmitCleanFaces(std::vector<OutTri3D>& out, const Manifold::Impl& in,
                    const BuildArrangement& A,
                    const std::vector<int>& face2cluster,
                    const std::vector<vec3>& seeds, double eps) {
  const int nTri = static_cast<int>(in.NumTri());
  // Build the winding broadphase once for this component (a proven exact
  // crossing-superset; the clean-face winding hot loop).
  const TriWindBVH bvh = BuildTriWindBVH(A.tri, in.bBox_);
  const std::vector<signed char> seedSign0 =
      PrecomputeSeedSign(A.tri, seeds[0]);
  // A face is "clean" only if it is neither seamed nor part of a coplanar
  // cluster (the fold owns cluster faces).
  auto isClean = [&](int t) { return !A.seamed[t] && face2cluster[t] < 0; };
  // Interior barycentric samples, centroid first (so a non-grazing face is
  // bitwise-identical to the old single-centroid probe); the rest are the
  // graze-dodge fallbacks on the SAME (constant-winding) cell above the
  // triangle.
  static constexpr double kBary[][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3},
                                        {0.6, 0.2, 0.2},
                                        {0.2, 0.6, 0.2},
                                        {0.2, 0.2, 0.6},
                                        {0.5, 0.3, 0.2},
                                        {0.2, 0.5, 0.3}};
  // The winding of a clean face is read-only on `in` and independent per face,
  // so the classification is order-free.  Two-pass keeps the emit ORDER (and
  // thus the output) bitwise-identical regardless of thread count: classify in
  // parallel (manifold::for_each_n / autoPolicy - sequential in a series
  // build), then append in ascending face index.
  //   status: -2 not-clean, -1 fail-closed (degenerate/SoS), 0 drop,
  //           1 emit original orientation, 2 emit reversed orientation.
  // BOTH-SIDES RETENTION (f4-r5, port of EmitSeamedFace's rule to clean faces):
  // probe w_S on BOTH sides of the face and retain iff EXACTLY ONE side is
  // inside {w_S>=1}; the solid side fixes the emitted orientation.  For an
  // oriented mult-1 clean face w_below = w_above + 1, so aboveIn!=belowIn holds
  // EXACTLY at w_above==0 with the solid on the -nHat side (belowIn) ->
  // original orientation: BITWISE-IDENTICAL to the former w_above==0 rule on
  // the whole corpus off openscad.  Under a coplanar coincidence (a clean tile
  // oppositely coincident with a seamed face the shares-vertex skip left
  // unclustered) w_above==w_below==0, and the one-sided rule over-emits it; the
  // both-sides read drops it.  F4R_ONESIDE reverts to the one-sided rule
  // (mutation lever).
  static const bool kOneSide = std::getenv("F4R_ONESIDE") != nullptr;
  auto classify = [&](int t) -> int {
    if (!isClean(t)) return -2;
    const double nLen = la::length(A.faceN[t]);
    if (!(nLen > 0.0)) return -1;
    const vec3 nHat = A.faceN[t] / nLen;
    std::vector<int> cands;  // per-task candidate scratch
    std::optional<int> g, gb;
    for (const auto& w : kBary) {
      const vec3 p =
          w[0] * A.tri[t][0] + w[1] * A.tri[t][1] + w[2] * A.tri[t][2];
      g = RobustWindingBVH(A.tri, bvh, p + eps * nHat, seeds, cands,
                           seedSign0.data());
      if (kOneSide) {
        if (g) break;  // any interior sample measures the (constant) cell
        continue;
      }
      gb = RobustWindingBVH(A.tri, bvh, p - eps * nHat, seeds, cands,
                            seedSign0.data());
      if (g && gb) break;  // both interior samples measure the constant cells
    }
    if (kOneSide) {
      if (!g) return -1;  // grazes at every interior point (SoS): fail closed
      return (*g == 0) ? 1 : 0;
    }
    if (!g || !gb) return -1;         // grazes everywhere (SoS): fail closed
    return BothSidesRetain(*g, *gb);  // 0 drop / 1 original / 2 reversed
  };
  std::vector<signed char> status(nTri, -2);
  for_each_n(autoPolicy(nTri, 256), countAt(0), static_cast<size_t>(nTri),
             [&](int t) { status[t] = static_cast<signed char>(classify(t)); });
  for (int t = 0; t < nTri; ++t) {
    if (status[t] == -1) return false;  // fail closed (same as the walk)
    if (status[t] != 1 && status[t] != 2) continue;
    // The solid side fixes orientation: status 1 (solid on -nHat) keeps the
    // face's original CCW (+nHat) winding; status 2 (solid on +nHat) reverses.
    // On the whole corpus off openscad every retained clean face is status 1
    // (mult-1 -> belowIn), so the reversed branch is a defensive completion,
    // never a byte-identity divergence.
    const bool rev = (status[t] == 2);
    // f4-junction: split the retained clean face's boundary at any registry
    // junction strictly interior to one of its edges (a neighbour's seam
    // endpoint / triple lands here), so the shared mesh edge is split on BOTH
    // sides and the emission fan closes.  A clean face is a single winding cell
    // (uniform coverage), so every sub-triangle carries its OWN retained
    // orientation (2D-CCW -> +nHat = original).  No interior junction -> the
    // single-triangle emit, byte-identical (the whole corpus off openscad).
    AxisDropFrame pf;
    if (A.junctions.empty() ||
        !BuildAxisDropFrame(A.faceN[t], A.tri[t][0], pf)) {
      if (rev)
        out.push_back({A.tri[t][0], A.tri[t][2], A.tri[t][1]});
      else
        out.push_back({A.tri[t][0], A.tri[t][1], A.tri[t][2]});
      continue;
    }
    const vec3 P[3] = {A.tri[t][0], A.tri[t][1], A.tri[t][2]};
    const int corner[3] = {pf.add(P[0]), pf.add(P[1]), pf.add(P[2])};
    std::vector<int> loop;
    bool anySplit = false;
    for (int e = 0; e < 3; ++e) {
      loop.push_back(corner[e]);
      for (const auto& js :
           JunctionSplitsOnSegment(pf.proj(P[e]), pf.proj(P[(e + 1) % 3]), P[e],
                                   P[(e + 1) % 3], A.junctions, eps)) {
        const int v = pf.addAt(js.first, js.second);
        if (v != loop.back()) {
          loop.push_back(v);
          anySplit = true;
        }
      }
    }
    if (!anySplit) {
      if (rev)
        out.push_back({A.tri[t][0], A.tri[t][2], A.tri[t][1]});
      else
        out.push_back({A.tri[t][0], A.tri[t][1], A.tri[t][2]});
      continue;
    }
    PolygonsIdx pidx(1);
    for (int idx : loop) pidx[0].push_back({pf.verts2[idx], idx});
    std::vector<ivec3> ctris;
    try {
      ctris = TriangulateIdx(pidx, eps);
    } catch (...) {
      return false;  // malformed split polygon: fail closed
    }
    for (const ivec3& tr : ctris) {
      if (rev)
        out.push_back({pf.canon3[tr.x], pf.canon3[tr.z], pf.canon3[tr.y]});
      else
        out.push_back({pf.canon3[tr.x], pf.canon3[tr.y], pf.canon3[tr.z]});
    }
  }
  return true;
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
//      (exact diagonal-split; earclip is not hole-safe), oriented per-triangle
//      toward the exterior; assembled by the ordinary BuildImpl eps-weld (the
//      offline run proved the once-rounded coordinated soup position-welds
//      CLOSED - identity plumbing through the weld is not needed: shared
//      doubles make every coincident vertex byte-equal).
namespace e1 {

// weld radius for the dust adjudications inside the triangulators (set by
// the engine before each run; single-threaded fallback path)
inline double kDustEps = 0.0;

using K3 = std::tuple<double, double, double>;
inline K3 KeyOf(const vec3& p) { return {p.x, p.y, p.z}; }

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
// sign-agnostic); the engine needs the standard CCW-positive convention
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

// Exact diagonal-split triangulation of a weakly-simple CCW polygon (collinear
// runs, pinch-repeated vertices, keyhole-duplicated bridges all allowed).
// poly = vertex indices into pos3; emits index triples.  Returns false when no
// valid diagonal exists (caller counts + fails closed via the census).
inline bool DiagSplitFwd(const std::vector<int>& poly,
                         const std::vector<vec3>& pos3, int axis,
                         std::vector<ivec3>& out, int depth);

// EARCLIP-first triangulation of a weakly-simple CCW polygon (the offline-
// validated combination): clip strictly-convex ears whose closed triangle
// contains no other polygon vertex (inclusive blocking - duplicates block
// conservatively), fall back to the exact diagonal splitter on a stall.
// The slit/lollipop walks (a chord traversed on both sides with structure at
// its end) stall a pure diagonal search but always expose ears elsewhere.
inline bool Triangulate(const std::vector<int>& loop,
                        const std::vector<vec3>& pos3, int axis, double eps,
                        std::vector<ivec3>& out) {
  std::vector<int> idx = loop;  // ids into pos3; slots may repeat positions
  while (idx.size() > 3) {
    const int m = static_cast<int>(idx.size());
    bool found = false;
    for (int k = 0; k < m && !found; ++k) {
      const int a = idx[(k + m - 1) % m], b = idx[k], c = idx[(k + 1) % m];
      if (O2(pos3[a], pos3[b], pos3[c], axis) <= 0) continue;
      bool ok = true;
      for (int j = 0; j < m && ok; ++j) {
        if (j == k || j == (k + m - 1) % m || j == (k + 1) % m) continue;
        const vec3& v = pos3[idx[j]];
        // twin instance of a corner is not a blocker - compare in the DRAWN
        // (2D) frame: distinct 3D identities can project to the identical 2D
        // point (they differ only along the dropped axis - measured: such a
        // co-projected pair blocked every adjacent ear)
        const vec2 v2 = Drop2(v, axis);
        if (v2 == Drop2(pos3[a], axis) || v2 == Drop2(pos3[b], axis) ||
            v2 == Drop2(pos3[c], axis))
          continue;
        const int d1 = O2(pos3[a], pos3[b], v, axis);
        const int d2 = O2(pos3[b], pos3[c], v, axis);
        const int d3 = O2(pos3[c], pos3[a], v, axis);
        if (d1 >= 0 && d2 >= 0 && d3 >= 0) ok = false;
      }
      if (!ok) continue;
      out.push_back({a, b, c});
      idx.erase(idx.begin() + k);
      found = true;
    }
    if (!found) {
      // stalled remainder: exact diagonal split; if THAT stalls (an eps-scale
      // bowtie from bent sub-chains - the rounded arrangement's residue),
      // accept the partial covering iff the remainder is WELD-DUST (width
      // below the weld radius: it welds away; the clipped macro ears stand).
      std::vector<ivec3> rem;
      if (DiagSplitFwd(idx, pos3, axis, rem, 0)) {
        out.insert(out.end(), rem.begin(), rem.end());
        return true;
      }
      double s = 0.0, ext = 0.0;
      vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
      const int mm = static_cast<int>(idx.size());
      for (int k = 0; k < mm; ++k) {
        const vec2 p1 = Drop2(pos3[idx[k]], axis);
        const vec2 p2 = Drop2(pos3[idx[(k + 1) % mm]], axis);
        s += p1.x * p2.y - p2.x * p1.y;
        lo = la::min(lo, p1);
        hi = la::max(hi, p1);
      }
      ext = std::max(hi.x - lo.x, hi.y - lo.y);
      return ext <= 0.0 || std::abs(0.5 * s) / ext <= 0.99 * eps;
    }
  }
  if (idx.size() == 3) {
    if (O2(pos3[idx[0]], pos3[idx[1]], pos3[idx[2]], axis) >= 0)
      out.push_back({idx[0], idx[1], idx[2]});
    return true;
  }
  return true;
}

inline bool DiagSplit(const std::vector<int>& poly,
                      const std::vector<vec3>& pos3, int axis,
                      std::vector<ivec3>& out, int depth = 0) {
  const int n = static_cast<int>(poly.size());
  if (n < 3 || depth > 4 * n + 64) return n < 3;
  if (n == 3) {
    // emit COLLINEAR leaves too: a rounded-degenerate cap triangle carries
    // REAL boundary sub-edges - (a,m),(m,b) stitch a subdivided neighbour to
    // an unsubdivided one through the weld (dropping them opened the fans:
    // the emitting cell classified boundary yet its edge had no partner).
    // A NEGATIVE leaf is a rounded fold-back: when dust-THIN its sub-edges
    // are equally real boundary (emit, loop order); a macro-negative leaf is
    // a genuine triangulation error (fail).
    const int s3 = O2(pos3[poly[0]], pos3[poly[1]], pos3[poly[2]], axis);
    if (s3 >= 0) {
      out.push_back({poly[0], poly[1], poly[2]});
      return true;
    }
    const vec2 q0 = Drop2(pos3[poly[0]], axis), q1 = Drop2(pos3[poly[1]], axis),
               q2 = Drop2(pos3[poly[2]], axis);
    const double a2 = std::abs(la::cross(q1 - q0, q2 - q0));
    const double L = std::max(
        {la::length(q1 - q0), la::length(q2 - q1), la::length(q0 - q2)});
    if (L > 0.0 && a2 / L <= kDustEps) {
      out.push_back({poly[0], poly[1], poly[2]});
      return true;
    }
    return false;
  }
  for (int i = 0; i < n; ++i) {
    const vec3& a = pos3[poly[(i + n - 1) % n]];
    const vec3& b = pos3[poly[i]];
    const vec3& c = pos3[poly[(i + 1) % n]];
    for (int jo = 2; jo <= n - 2; ++jo) {
      const int j = (i + jo) % n;
      const vec3& d = pos3[poly[j]];
      if (Drop2(d, axis) == Drop2(b, axis)) continue;  // 2D twin (bridge)
      // in-cone at i (collinear prev/next treated as convex half-plane)
      if (O2(a, b, c, axis) >= 0) {
        if (!(O2(b, c, d, axis) > 0 && O2(b, d, a, axis) > 0)) continue;
      } else {
        if (O2(b, c, d, axis) <= 0 && O2(b, d, a, axis) <= 0) continue;
      }
      bool ok = true;
      for (int k = 0; k < n && ok; ++k) {
        const int k2 = (k + 1) % n;
        if (k == i || k2 == i || k == j || k2 == j) continue;
        if (ProperCross2(b, d, pos3[poly[k]], pos3[poly[k2]], axis)) ok = false;
      }
      for (int k = 0; k < n && ok; ++k) {
        if (k == i || k == j) continue;
        const vec3& v = pos3[poly[k]];
        if (Drop2(v, axis) == Drop2(b, axis) ||
            Drop2(v, axis) == Drop2(d, axis))
          continue;
        if (OnOpenSeg2(b, d, v, axis)) ok = false;
      }
      if (!ok) continue;
      std::vector<int> p1, p2;
      for (int k = i;; k = (k + 1) % n) {
        p1.push_back(poly[k]);
        if (k == j) break;
      }
      for (int k = j;; k = (k + 1) % n) {
        p2.push_back(poly[k]);
        if (k == i) break;
      }
      return DiagSplit(p1, pos3, axis, out, depth + 1) &&
             DiagSplit(p2, pos3, axis, out, depth + 1);
    }
  }
  return false;
}

inline bool DiagSplitFwd(const std::vector<int>& poly,
                         const std::vector<vec3>& pos3, int axis,
                         std::vector<ivec3>& out, int depth) {
  return DiagSplit(poly, pos3, axis, out, depth);
}

}  // namespace e1

StageResult<Manifold::Impl> EmitCoordinatedBoundaryImpl(
    const Manifold::Impl& in, const BuildArrangement& A, double eps,
    const std::vector<std::pair<int, vec3>>* extraJ,
    std::vector<std::pair<int, vec3>>* collectX,
    std::vector<vec3>* collectT = nullptr) {
  using e1::K3;
  using e1::KeyOf;
  const int nTri = static_cast<int>(A.tri.size());
  const std::vector<vec3> seeds = WindingSeeds(in.bBox_);
  auto fail = [](const char* msg) {
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved, msg);
  };
  static const bool kDump = std::getenv("E1_DUMP") != nullptr;
  e1::kDustEps = eps;

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
  struct ESeam {
    vec3 p0, p1;
    int other;
  };
  std::vector<std::vector<ESeam>> eSeams(nTri);
  {
    std::map<std::pair<int, int>, int> sideCache;  // (vert id, gid) -> sign
    std::map<std::tuple<int, int, int>, vec3> pierceCache;
    auto sideOf = [&](int f, int k, int q) -> int {
      const auto key = std::make_pair(A.vid[f][k], q);
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
        std::vector<vec3> cand;
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
              const auto key = std::make_tuple(v0, v1, q);
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
    if (kDump) {
      int nseam = 0;
      for (const auto& v : eSeams) nseam += static_cast<int>(v.size());
      std::fprintf(stderr, "E1 engine seams=%d (pair-records)\n", nseam);
    }
  }

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
  std::vector<int> outG;  // per-tri emitting group (diagnostics)
  int dustTri = 0, triFail = 0, spliceFail = 0;
  int nCells = 0, nNeg = 0, nJump = 0, nBoundary = 0, nOwned = 0, nDustCell = 0;
  const double scale = in.bBox_.Scale();

  for (int g = 0; g < nG; ++g) {
    // ---- 2. the group's segment set (3D endpoint pairs, shared doubles) ----
    struct Seg {
      vec3 p0, p1;
      int planeQ;  // partner group (-1 = member triangle edge)
      int fOwn, fOther;
    };
    std::vector<Seg> segs;
    for (const int f : members[g]) {
      for (int e = 0; e < 3; ++e)
        segs.push_back({A.tri[f][e], A.tri[f][(e + 1) % 3], -1, f, -1});
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
      for (int e = 0; e < 3; ++e)
        segs.push_back({A.tri[f2][e], A.tri[f2][(e + 1) % 3], -1, f2, -1});
    }
    const int nS = static_cast<int>(segs.size());

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
    const double rho =
        64.0 * std::numeric_limits<double>::epsilon() * (1.0 + scale);
    std::map<std::tuple<long long, long long, long long>, std::vector<vec3>>
        canonGrid;
    auto canonV = [&](const vec3& v) -> vec3 {
      const long long cx = static_cast<long long>(std::floor(v.x / rho));
      const long long cy = static_cast<long long>(std::floor(v.y / rho));
      const long long cz = static_cast<long long>(std::floor(v.z / rho));
      for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
          for (long long dz = -1; dz <= 1; ++dz) {
            const auto it = canonGrid.find({cx + dx, cy + dy, cz + dz});
            if (it == canonGrid.end()) continue;
            for (const vec3& w : it->second)
              if (la::length(w - v) <= rho) return w;
          }
      canonGrid[{cx, cy, cz}].push_back(v);
      return v;
    };
    for (const Seg& s : segs) {  // seed: endpoints are the canonical anchors
      canonV(s.p0);
      canonV(s.p1);
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
      return true;
    };
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
    // pool = the group's OWN segment endpoints (exact once-only constructions,
    // ULP-accurate - the offline pool).  NOT A.junctions: its eps-deduped
    // representatives sit up to eps off the exact lines, and inserting them
    // with an eps window BENDS the chains by eps - the bent sub-chains then
    // properly cross and fragment the walk (measured: a split cascade).
    // OUTER FIXPOINT: the exchange and the completion feed each other - a
    // completion-added crossing (ill-conditioned on near-collinear pairs:
    // the same line's overlapping segments each get their own noisy
    // crossing position) must be EXCHANGED onto every collinear twin, and
    // exchanged points can expose new crossings.  (Measured: identical-
    // endpoint twin segments carrying different completion splits 2.4e-5
    // apart - the T-junction/lens class.)
    for (int outer = 0; outer < 4; ++outer) {
      size_t nsplit0 = 0;
      for (int i = 0; i < nS; ++i) nsplit0 += splits[i].size();
      // T-junction pool pass, iterated to FIXPOINT with ALL SPLIT POINTS in
      // the pool: near-collinear overlapping chains (the same line reached via
      // different member pairs, ULP apart) must carry IDENTICAL subdivisions -
      // endpoint-only exchange left one chain split where its twin spanned
      // whole (measured: intra-group T-junctions with the on-vertex 9e-16 off
      // the unsplit edge).  Split values are canonical (snap grid), so the
      // exchange converges.
      // 8*rho: the canonical snap grid moves split points up to rho off their
      // segment lines, so on-line tests must budget the snap displacement
      const double tolLine = 8.0 * rho;
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
      // RETRY SPLITS (second pass): the first pass's failing walks' observed
      // self-crossing points, injected group-globally so every incident cell
      // subdivides identically.
      if (extraJ) {
        int landed = 0;
        for (int i = 0; i < nS; ++i) {
          const Seg& s = segs[i];
          const vec3 d3 = s.p1 - s.p0;
          const double len2 = la::dot(d3, d3);
          if (!(len2 > 0.0)) continue;
          const double tolLine2 = 8.0 * rho;
          for (const auto& gv : *extraJ) {
            // group-scoped (loop-derived points splash); g==-1 = GLOBAL entries
            // (canonical triple constructions: exact, safe everywhere)
            if (gv.first != g && gv.first != -1) continue;
            const vec3& V = gv.second;
            const vec3 w = V - s.p0;
            const double t = la::dot(w, d3) / len2;
            if (!(t > 0.0 && t < 1.0)) continue;
            if (la::length(w - t * d3) > tolLine2) continue;
            if (addSplit(i, V)) ++landed;
          }
        }
        if (kDump && landed)
          std::fprintf(stderr, "E1 retry g=%d landed=%d/%d\n", g, landed,
                       static_cast<int>(extraJ->size()));
      }

      // PLANARITY COMPLETION (all remaining segment-pair crossings): the exact
      // seam-x-seam enumeration and the endpoint pool cover the canonical
      // crossings, but the drawn (rounded) graph must be PLANAR for the face
      // walk - overlapping coplanar members (the fold structure) cross member
      // edges and same-line seams in ways the passes above miss (measured:
      // properly-crossing sub-edges -> bowtie walks -> untriangulable cells).
      // Detect every remaining proper crossing exactly on the shared rounded
      // endpoints and split both segments; the split point uses the canonical
      // triple when both carriers are seams of distinct planes, else the
      // in-segment interpolation (identity across groups holds within weld
      // tolerance via the shared-endpoint constructions).
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
          if (!e1::ProperCross2(segs[i].p0, segs[i].p1, segs[j].p0, segs[j].p1,
                                axis))
            continue;
          vec3 X;
          bool have = false;
          if (segs[i].planeQ >= 0 && segs[j].planeQ >= 0 &&
              segs[i].planeQ != segs[j].planeQ)
            have = triplePos(g, segs[i].planeQ, segs[j].planeQ, X);
          if (!have) {
            // in-plane 2D crossing, interpolated along segment i's 3D span
            const double dax = a1.x - a0.x, day = a1.y - a0.y;
            const double dbx = b1.x - b0.x, dby = b1.y - b0.y;
            const double den = dax * dby - day * dbx;
            if (den == 0.0) continue;
            const double t = ((b0.x - a0.x) * dby - (b0.y - a0.y) * dbx) / den;
            X = segs[i].p0 + t * (segs[i].p1 - segs[i].p0);
          }
          addSplit(i, X);
          addSplit(j, X);
        }
      }

      size_t nsplit1 = 0;
      for (int i = 0; i < nS; ++i) nsplit1 += splits[i].size();
      if (nsplit1 == nsplit0) break;
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
    // ---- 4. 2D graph (verts keyed by 3D bits) + exact rotation walk ----
    std::map<K3, int> vidOf;
    std::vector<vec3> pos3;
    auto vid = [&](const vec3& p) -> int {
      auto it = vidOf.find(KeyOf(p));
      if (it == vidOf.end()) {
        it = vidOf.emplace(KeyOf(p), static_cast<int>(pos3.size())).first;
        pos3.push_back(p);
      }
      return it->second;
    };
    std::vector<std::set<int>> adj;
    auto link = [&](int a, int b) {
      if (a == b) return;
      const int mx = std::max(a, b);
      if (static_cast<int>(adj.size()) <= mx) adj.resize(mx + 1);
      adj[a].insert(b);
      adj[b].insert(a);
    };
    for (int i = 0; i < nS; ++i) {
      std::sort(splits[i].begin(), splits[i].end(),
                [](const auto& x, const auto& y) { return x.first < y.first; });
      int prev = vid(segs[i].p0);
      for (const auto& pr : splits[i]) {
        const int v = vid(pr.second);
        link(prev, v);
        prev = v;
      }
      link(prev, vid(segs[i].p1));
    }
    adj.resize(pos3.size());
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
          {
            std::map<int, int> seen;
            int pi = -1, pk = -1;
            for (int k = 0; k < m && pi < 0; ++k) {
              const auto it = seen.find(L[k]);
              if (it != seen.end()) {
                pi = it->second;
                pk = k;
              } else {
                seen.emplace(L[k], k);
              }
            }
            if (pi >= 0) {
              std::vector<int> l1(L.begin() + pi, L.begin() + pk);
              std::vector<int> l2(L.begin(), L.begin() + pi);
              l2.insert(l2.end(), L.begin() + pk, L.end());
              work.push_back(std::move(l1));
              work.push_back(std::move(l2));
              continue;
            }
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
          double s = 0.0;
          for (size_t k = 0; k < L.size(); ++k) {
            const vec2 p1 = e1::Drop2(pos3[L[k]], axis);
            const vec2 p2 = e1::Drop2(pos3[L[(k + 1) % L.size()]], axis);
            s += p1.x * p2.y - p2.x * p1.y;
          }
          if (s > 0)
            cells.push_back(std::move(L));
          else if (s < 0)
            negloops.push_back(std::move(L));
        }
      }

    // ---- 5. disconnected island hole-rings: containment + keyhole ----
    auto loopArea = [&](const std::vector<int>& loop) -> double {
      double s = 0.0;
      for (size_t k = 0; k < loop.size(); ++k) {
        const vec2 p1 = e1::Drop2(pos3[loop[k]], axis);
        const vec2 p2 = e1::Drop2(pos3[loop[(k + 1) % loop.size()]], axis);
        s += p1.x * p2.y - p2.x * p1.y;
      }
      return s;
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
        const vec2 p = e1::Drop2(pos3[nl[0]], axis);
        int best = -1;
        double bestA = 0.0;
        for (size_t ci = 0; ci < cells.size(); ++ci)
          if (inLoop(p, cells[ci])) {
            const double a = std::abs(loopArea(cells[ci]));
            if (best < 0 || a < bestA) {
              best = static_cast<int>(ci);
              bestA = a;
            }
          }
        // DUST HOLE RING: a ring whose width is below the weld radius is a
        // sub-representable near-duplicate zigzag (measured: macro-long,
        // ~1e-12-wide rings in the near-tangent fold overlap); splicing it
        // strangles the triangulation on sub-ULP structure, and at the weld
        // it vanishes anyway - drop the ring, keep the cell solid.
        {
          double s = 0.0, ext = 0.0;
          vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
          for (size_t k = 0; k < nl.size(); ++k) {
            const vec2 p1 = e1::Drop2(pos3[nl[k]], axis);
            const vec2 p2 = e1::Drop2(pos3[nl[(k + 1) % nl.size()]], axis);
            s += p1.x * p2.y - p2.x * p1.y;
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

    // ---- 6. classify + emit ----
    nCells += static_cast<int>(cells.size());
    nNeg += static_cast<int>(negloops.size());
    for (const std::vector<int>& loop : cells) {
      std::vector<ivec3> tris;
      if (!e1::Triangulate(loop, pos3, axis, eps, tris)) {
        // ALL-OR-NOTHING: a partial covering emits unpaired interior edges
        // (the offline lesson) - discard, then adjudicate by WIDTH: a loop
        // whose area/extent is below the weld scale is a rounded-degenerate
        // sliver (non-simple at double precision) that welds away; anything
        // wider is an honest triangulation failure (fail-closed below).
        tris.clear();
        double s = 0.0, ext = 0.0;
        vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
        for (size_t k = 0; k < loop.size(); ++k) {
          const vec2 p1 = e1::Drop2(pos3[loop[k]], axis);
          const vec2 p2 = e1::Drop2(pos3[loop[(k + 1) % loop.size()]], axis);
          s += p1.x * p2.y - p2.x * p1.y;
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
        if (collectX) {
          // planarize-detect: record the failing walk's own proper crossings
          const int mm = static_cast<int>(loop.size());
          for (int u = 0; u < mm; ++u) {
            const vec3 &ua = pos3[loop[u]], &ub = pos3[loop[(u + 1) % mm]];
            for (int v = u + 1; v < mm; ++v) {
              const vec3 &va = pos3[loop[v]], &vb = pos3[loop[(v + 1) % mm]];
              if (KeyOf(ua) == KeyOf(va) || KeyOf(ua) == KeyOf(vb) ||
                  KeyOf(ub) == KeyOf(va) || KeyOf(ub) == KeyOf(vb))
                continue;
              if (!e1::ProperCross2(ua, ub, va, vb, axis)) continue;
              const vec2 a0 = e1::Drop2(ua, axis), a1 = e1::Drop2(ub, axis);
              const vec2 b0 = e1::Drop2(va, axis), b1 = e1::Drop2(vb, axis);
              const double dax = a1.x - a0.x, day = a1.y - a0.y;
              const double dbx = b1.x - b0.x, dby = b1.y - b0.y;
              const double den = dax * dby - day * dbx;
              if (den == 0.0) continue;
              const double t =
                  ((b0.x - a0.x) * dby - (b0.y - a0.y) * dbx) / den;
              collectX->push_back({g, ua + t * (ub - ua)});
            }
          }
        }
        ++triFail;
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
      for (int f2 = 0; f2 < nTri; ++f2) {
        if (gid[f2] == g) continue;
        const Box& b = fbox[f2];
        const double m = 1e-6 * (1.0 + scale);
        if (cenP.x < b.min.x - m || cenP.x > b.max.x + m ||
            cenP.y < b.min.y - m || cenP.y > b.max.y + m ||
            cenP.z < b.min.z - m || cenP.z > b.max.z + m)
          continue;
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
      for (const auto& dn : nearD) {
        if (dn.first <= T) continue;
        // HARD CAP at the weld radius: a stack may only absorb sheets that
        // weld together anyway (unrepresentably close).  Un-capped 8x
        // chaining absorbed ladders of REPRESENTABLE distinct sheets into
        // one net emission (measured: systematic opens across dozens of
        // groups wherever the sheet-distance ladder had no 8x gap).
        if (dn.first <= 8.0 * T && dn.first <= 0.99 * eps)
          T = dn.first;
        else
          break;
      }
      const double off = 2.0 * T;
      bool owned = false;
      for (const auto& dn : nearD) {
        if (dn.first > T) break;
        if (!covers(dn.second)) continue;
        if (gid[dn.second] < g) {
          owned = true;
          break;
        }
        jump += la::dot(A.faceN[dn.second], nHat) >= 0.0 ? 1 : -1;
      }
      if (owned) {
        ++nOwned;
        continue;  // a lower covering group owns this stack cell
      }
      if (jump == 0) continue;  // net-cancelled: not a sheet
      const std::optional<int> wA = RobustWinding(in, cenP + off * nHat, seeds);
      const std::optional<int> wB = RobustWinding(in, cenP - off * nHat, seeds);
      if (!wA || !wB) {
        if (kDump)
          std::fprintf(stderr,
                       "E1 FAIL probe g=%d cen=(%.9g,%.9g,%.9g) off=%.3g\n", g,
                       cen3.x, cen3.y, cen3.z, off);
        return fail("e1: winding probe filter-uncertain - fail-closed");
      }
      const bool certified = (*wB - *wA == jump);
      // COMPLETENESS CERTIFICATE: probed delta == combinatorial covering jump
      if (!certified) {
        // a DUST cell - WIDTH below the weld scale - cannot be probed (its
        // interior hugs seam-line structure closer than the probe window, and
        // containment of off-plane points is unreliable within that scale)
        // and welds away regardless: its long sides are closer than eps, so
        // the assembly weld collapses it and the neighbours carry the
        // boundary.  Width = altitude of the largest sub-triangle.
        double eMax = 0.0;
        {
          const vec2 q0 = e1::Drop2(pos3[tris[best].x], axis);
          const vec2 q1 = e1::Drop2(pos3[tris[best].y], axis);
          const vec2 q2 = e1::Drop2(pos3[tris[best].z], axis);
          eMax = std::max(
              {la::length(q1 - q0), la::length(q2 - q1), la::length(q0 - q2)});
        }
        const double ext = eMax > 0.0 ? bestA / eMax : 0.0;  // = width
        if (ext <= 0.99 * eps) continue;  // sub-weld-width sliver: dust
        if (kDump) {
          std::fprintf(stderr,
                       "E1 FAIL cert g=%d jump=%d wA=%d wB=%d n=%d off=%.3g "
                       "ext=%.3g cen=(%.9g,%.9g,%.9g) d0=%.3g\n",
                       g, jump, *wA, *wB, static_cast<int>(loop.size()), off,
                       ext, cen3.x, cen3.y, cen3.z, d0);
          for (const auto& dn : nearD) {
            if (dn.first > T) break;
            const int f2 = dn.second;
            std::fprintf(stderr,
                         "  E1 stack f=%d gid=%d dist=%.3g covers=%d sgn=%d\n",
                         f2, gid[f2], dn.first, covers(f2) ? 1 : 0,
                         la::dot(A.faceN[f2], nHat) >= 0.0 ? 1 : -1);
            std::fprintf(stderr,
                         "    tri |N|=%.3g v0=(%.9g,%.9g,%.9g) "
                         "v1=(%.9g,%.9g,%.9g) v2=(%.9g,%.9g,%.9g)\n",
                         la::length(A.faceN[f2]), A.tri[f2][0].x,
                         A.tri[f2][0].y, A.tri[f2][0].z, A.tri[f2][1].x,
                         A.tri[f2][1].y, A.tri[f2][1].z, A.tri[f2][2].x,
                         A.tri[f2][2].y, A.tri[f2][2].z);
            const int oA0 = e1::O2(A.tri[f2][0], A.tri[f2][1], cenP, axis);
            const int oA1 = e1::O2(A.tri[f2][1], A.tri[f2][2], cenP, axis);
            const int oA2 = e1::O2(A.tri[f2][2], A.tri[f2][0], cenP, axis);
            std::fprintf(stderr, "    o=%d,%d,%d\n", oA0, oA1, oA2);
          }
          for (const int f : members[g]) {
            const int o0 = e1::O2(A.tri[f][0], A.tri[f][1], cenP, axis);
            const int o1 = e1::O2(A.tri[f][1], A.tri[f][2], cenP, axis);
            const int o2 = e1::O2(A.tri[f][2], A.tri[f][0], cenP, axis);
            const bool neg = o0 < 0 || o1 < 0 || o2 < 0;
            const bool pos = o0 > 0 || o1 > 0 || o2 > 0;
            if (!(neg && pos) || (!neg == !pos))
              std::fprintf(stderr,
                           "  E1 member f=%d sgn=%d o=%d,%d,%d covers=%d\n", f,
                           fsgn[f], o0, o1, o2, !(neg && pos) ? 1 : 0);
          }
          for (double mul = 1.0; mul <= 1000.0; mul *= 10.0) {
            const std::optional<int> wa =
                RobustWinding(in, cenP + mul * off * nHat, seeds);
            const std::optional<int> wb =
                RobustWinding(in, cenP - mul * off * nHat, seeds);
            std::fprintf(stderr, "  E1 ladder off=%.3g wA=%d wB=%d\n",
                         mul * off, wa ? *wa : -99, wb ? *wb : -99);
          }
        }
        return fail(
            "e1: coordinated-arrangement completeness certificate failed "
            "(winding delta != covering jump) - fail-closed");
      }
      if (probeHit)
        std::fprintf(stderr,
                     "E1 PROBEAT g=%d jump=%d wA=%d wB=%d off=%.3g cert=%d\n",
                     g, jump, *wA, *wB, off, certified ? 1 : 0);
      const bool aIn = *wA >= 1, bIn = *wB >= 1;
      if (aIn == bIn) continue;         // not a {w>=1} boundary here
      const int orient = bIn ? 1 : -1;  // +1: solid below, outward = +nHat
      ++nBoundary;
      for (const ivec3& t : tris) {
        const vec3 nr = la::cross(pos3[t.y] - pos3[t.x], pos3[t.z] - pos3[t.x]);
        const double sd = la::dot(nr, Nrep);
        // sd == 0 (rounded-degenerate): keep it - the cap triangle's edges
        // are real boundary; loop order is the cell's CCW, so orient decides
        if ((sd > 0.0 || (sd == 0.0 && orient > 0)) == (orient > 0))
          out.push_back({{pos3[t.x], pos3[t.y], pos3[t.z]}});
        else
          out.push_back({{pos3[t.x], pos3[t.z], pos3[t.y]}});
        outG.push_back(g);
      }
    }
  }
  if (kDump)
    std::fprintf(stderr,
                 "E1 emitted=%d dust=%d triFail=%d spliceFail=%d cells=%d "
                 "neg=%d jump=%d owned=%d boundary=%d dustCell=%d\n",
                 static_cast<int>(out.size()), dustTri, triFail, spliceFail,
                 nCells, nNeg, nJump, nOwned, nBoundary, nDustCell);
  if (collectT)
    for (const auto& kv : tripleCache)
      if (kv.second.first) collectT->push_back(kv.second.second);
  if (triFail > 0 || spliceFail > 0)
    return fail("e1: cell triangulation/hole-splice incomplete - fail-closed");
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
  return BuildImpl(out, eps);
}

// Two-pass driver: pass 1 collects the failing walks' observed self-crossing
// points (residual drawn-graph non-planarity at the representable-thin
// scale); pass 2 injects them as group-global splits.  Converges because the
// injected points lie ON the crossing sub-edges (both incident cells split
// identically); a second failure is an honest fail-closed.
StageResult<Manifold::Impl> EmitCoordinatedBoundary(const Manifold::Impl& in,
                                                    const BuildArrangement& A,
                                                    double eps) {
  std::vector<std::pair<int, vec3>> crossX;
  std::vector<vec3> triples;
  StageResult<Manifold::Impl> r =
      EmitCoordinatedBoundaryImpl(in, A, eps, nullptr, &crossX, &triples);
  if (!r.fatal) return r;
  // GLOBAL TRIPLE INJECTION: a committed triple splits the two seams of the
  // group that DISCOVERED it; the partner groups' copies of shared lines
  // must split at the identical canonical position (the 4-tri existence test
  // is member-pair-specific and legitimately asymmetric across groups).
  // Triples are exact canonical constructions - global injection with the
  // rounding-scale on-line tolerance cannot bend chains.
  for (const vec3& t : triples) crossX.push_back({-1, t});
  for (int pass = 0; pass < 6 && r.fatal && !crossX.empty(); ++pass) {
    const size_t before = crossX.size();
    r = EmitCoordinatedBoundaryImpl(in, A, eps, &crossX, &crossX);
    if (crossX.size() == before) break;  // no new crossings: converged/stuck
  }
  return r;
}

// THE BUILD driver: fold exactly-coplanar clusters in-plane, emit seamed
// sub-faces + clean faces, assemble + weld.  Returns the regularized Impl, or a
// fatal.
StageResult<Manifold::Impl> EmitComponentBoundary(
    const Manifold::Impl& in, const BuildArrangement& A,
    const std::vector<int>& face2cluster, double eps) {
  const std::vector<vec3> seeds = WindingSeeds(in.bBox_);

  std::vector<OutTri3D> emitted;
  bool ok = true;
  // Exact-coplanar clusters first (transversal seams on their faces are the
  // entanglement decline).
  FoldCoplanarClusters(emitted, in, A, face2cluster, seeds, eps, ok);
  if (!ok)
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "resolver: exact-coplanar fold declined (coplanar/transversal "
        "entanglement, degenerate projection, or filter-uncertain classify) - "
        "fail-closed");
  const int nTri = static_cast<int>(in.NumTri());
  static const bool kF4BDump = std::getenv("F4B_DUMP") != nullptr;
  if (kF4BDump) {  // census: run EVERY seamed face, count the fail-closed
                   // branch
    gF4BSeam = F4BSeamCensus{};
    bool anyFail = false;
    for (int f = 0; f < nTri; ++f) {
      if (!A.seamed[f]) continue;
      bool fok = true;
      EmitSeamedFace(emitted, A, f, in, seeds, eps, fok);
      if (!fok) anyFail = true;
    }
    std::fprintf(stderr,
                 "F4B_SEAM b3faces=%d b3verts=%d b4faces=%d b5faces=%d "
                 "okfaces=%d\n",
                 gF4BSeam.b3faces, gF4BSeam.b3verts, gF4BSeam.b4faces,
                 gF4BSeam.b5faces, gF4BSeam.okfaces);
    ok = !anyFail;
  } else {
    for (int f = 0; f < nTri && ok; ++f)
      if (A.seamed[f]) EmitSeamedFace(emitted, A, f, in, seeds, eps, ok);
  }
  if (!ok)
    // The resolver declined to build this face's arrangement exactly: a
    // >2-sheet triple point, a coplanar/degenerate projection, a malformed cell
    // walk, or a filter-uncertain classify probe (the SoS axis).  Negative
    // winding is NOT a decline - the witness rule absorbs it (w_above==0
    // retain).  Fail closed
    // - never emit geometry the resolver could not verify.
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "resolver: seam sub-face arrangement not exactly resolvable "
        "(triple point / degenerate / filter-uncertain) - fail-closed");
  if (!EmitCleanFaces(emitted, in, A, face2cluster, seeds, eps))
    return StageResult<Manifold::Impl>::Fatal(
        FatalReason::DirtyComponentUnresolved,
        "resolver: clean-face winding probe was filter-uncertain (SoS) - "
        "fail-closed");
  return BuildImpl(emitted, eps);
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
  // symwalk: build the exact HPoint provenance of every wedge-face arrangement
  // vertex (input / triple / seam endpoint), so the wedge overlay's cell walk
  // orders by HomogOrient2D on symbolic points.  No-op unless a wedge face is
  // flagged (empty on the whole resolving corpus).
  BuildSymbolicProvenance(A, in, eps);
  // f4-junction: gather the once-only junction registry (seam endpoints +
  // triples) so every emit path splits its edges at the non-proper-crossing
  // junctions the triple enumeration misses (no-op off openscad).
  BuildJunctionRegistry(A, eps);
  StageResult<Manifold::Impl> r =
      EmitComponentBoundary(in, A, face2cluster, eps);
  // E1 COORDINATED ENGINE (e1engine): the fallback resolver for the component
  // class the per-face emission fails closed on.  Reached ONLY on a fatal, so
  // byte-clean on every resolving carrier by construction; on its own failure
  // the ORIGINAL fatal is preserved (fail-closed never weakened).  Dev-gated:
  // E1_ENGINE=1 enables (default OFF until the openscad closure flips the pin).
  if (r.fatal && std::getenv("E1_ENGINE") != nullptr) {
    StageResult<Manifold::Impl> e1r = EmitCoordinatedBoundary(in, A, eps);
    if (!e1r.fatal) return e1r;
  }
  return r;
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
                 if (GateComponent(bImpl) != GateVerdict::Clean) {
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

CleanFaceProbe ClassifyCleanFaces_Probe(const Manifold::Impl& soup) {
  CleanFaceProbe out;
  double eps = EpsilonFromScale(soup.bBox_.Scale(), 1000);
  if (!(eps > 0.0) || !std::isfinite(eps)) return out;
  // Mirror ResolveComponent's prefix so the classify inputs match production
  // exactly (snap near-coplanar, detect exact-coplanar clusters, record seams).
  StageResult<Manifold::Impl> snapped = SnapNearCoplanarClusters(soup, eps);
  if (snapped.fatal) return out;
  const Manifold::Impl& in = snapped.value ? *snapped.value : soup;
  const std::vector<int> face2cluster = DetectCoplanarClusters(in);
  const BuildArrangement A = RecordSeams(in, face2cluster, eps);
  const std::vector<vec3> seeds = WindingSeeds(in.bBox_);
  // Run the REAL (compiled) EmitCleanFaces and record which clean faces it kept
  // by matching their exact vertex triple (clean faces emit verbatim, unsplit),
  // so this hook reflects whichever classification is compiled - the mutation
  // (revert to the flood) flips `kept` for the mislabeled boundary faces.
  std::vector<OutTri3D> emitted;
  EmitCleanFaces(emitted, in, A, face2cluster, seeds, eps);
  using TriKey = std::tuple<double, double, double, double, double, double,
                            double, double, double>;
  auto keyOf = [](const vec3& a, const vec3& b, const vec3& d) {
    return TriKey{a.x, a.y, a.z, b.x, b.y, b.z, d.x, d.y, d.z};
  };
  std::set<TriKey> keptSet;
  for (const OutTri3D& t : emitted)
    keptSet.insert(keyOf(t.v[0], t.v[1], t.v[2]));
  const int nTri = static_cast<int>(in.NumTri());
  for (int t = 0; t < nTri; ++t) {
    if (A.seamed[t] || face2cluster[t] >= 0) continue;  // clean only
    out.faceIdx.push_back(t);
    const vec3 cen = (A.tri[t][0] + A.tri[t][1] + A.tri[t][2]) / 3.0;
    const double nLen = la::length(A.faceN[t]);
    std::optional<int> g;
    if (nLen > 0.0)
      g = RobustWinding(in, cen + eps * (A.faceN[t] / nLen), seeds);
    out.ownWinding.push_back(g ? *g : kWindingUncertain);
    out.kept.push_back(
        keptSet.count(keyOf(A.tri[t][0], A.tri[t][1], A.tri[t][2])) ? 1 : 0);
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
