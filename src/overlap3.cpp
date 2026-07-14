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

// Forward decl (defined below): per-face exactly-coplanar overlap cluster id
// (or -1).  Used by GateComponent to route a WITHIN-component coplanar-overlap
// component to B (the R2(i) blind spot the self-intersection test misses).
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
  // run on THIS component only, routes it to B where the exact-coplanar fold
  // consumes the overlap.  This is WITHIN-component by construction: `comp` is
  // one connectivity component, so a CROSS-component coplanar overlap (two
  // distinct components that happen to coincide) is invisible here BY DESIGN -
  // the non-fusion posture (fusion is the Boolean's job).  A clean solid with
  // no internal coplanar overlap detects nothing and stays Clean (bitwise
  // pass-through); the cost is the prefiltered scan, proportional to the input.
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

// ---------------------------------------------------------------------------
// SINGLE GLOBAL TIE-BREAK CONVENTION (SoS), docs/Regularize3D.md stage 6.
// ONE exact integer implementation.  The perturbed 4x4 orient3d determinant is
// sum_K coeff_K * e^K; the SoS sign is the sign of the LOWEST-K nonzero
// coefficient.  Perturb matrix entry (row r, coord c) by e^(2^(rank(r)*3 +
// (2-c))), rank = order of the four points' GLOBAL vertex indices; the twelve
// keys are DISTINCT powers of two, so K is a bitmask of the perturbed entries
// and every coefficient coeff_K is a signed sum of products of <= 3 of the
// twelve coordinate MANTISSAS (a minor of degree <= 3 over the same windowed
// coords).  So the whole cascade - the e^0 exact orient3d AND the higher-e
// perturbation terms - evaluates on ONE integer path: decompose each coord to
// (mantissa*2^exp) via frexp, form each term as a 192-bit |product of <= 3
// mantissas| at a 2-exponent, align a coefficient's terms to their min
// exponent, and accumulate the sign in a two's-complement bigint.  The
// accumulator width is ADAPTIVE (limbs computed from the term-exponent spread)
// and sized to the worst-case finite-double spread, so there is NO window-fail
// refusal: the sign is TOTAL for every finite-double input (0 means an exact
// geometric tie, never "uncertain").  Local-rank reduction: the map (local rank
// r)->(global rank G(r)) is strictly monotonic on the bit positions r*3+(2-c),
// so the leading monomial - hence the sign - is invariant to using local ranks
// 0..3 in place of the true global ranks; a per-predicate local computation
// decides consistently with ONE global perturbation.  This is the ONLY
// tie-break convention, threaded through every enumeration predicate so a ray
// or edge grazing a shared boundary resolves identically for every probe.
// Validated (reg3d-s6r notebook): BITWISE-identical to the prior
// Shewchuk-expansion cascade on the mesh domain (same-scale + mixed-magnitude,
// e^0 and SoS, zero mismatch), and exact-correct vs an arbitrary-precision
// oracle on wild inputs where the expansion overflowed double; the fixed-window
// predecessor refused a zero-straddling adversary this path decides.  Only
// invoked on the RARE filter-uncertain (0) fallback, so the certified fast path
// (and siA/siB, fully certified) is untouched.
//
// RELUCTANT ACCEPTANCE (owner contract): this integer exact kernel is net-new
// surface the design had declined ("no exact kernel in the tree",
// docs/Regularize3D.md B mechanism).  It is accepted NARROWLY as the stage-6
// filter-0 fallback - one implementation, single call site discipline, off
// every certified path.  QUEUED FOR REVISIT (docs/Regularize3D.md open list):
// whether the arrangement can be structured to avoid needing an exact orient3d
// kernel at all remains an open question; this kernel is the current,
// reluctantly-accepted answer, not a settled one.
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
// A determinant term: sign * |product of <= 3 mantissas| * 2^e.  |product| <
// 2^159 (three 53-bit mantissas) fits in three 64-bit limbs.
struct Term {
  uint64_t mag[3];
  long e;
  int sign;
};
// Adaptive two's-complement accumulator width.  Each coord is m*2^E with |m| <
// 2^53 and E = ex - 53, ex the frexp exponent in [-1073, 1024], so E in
// [-1126, 971] and a three-factor exponent esum in [-3378, 2913]: the term
// spread is <= 6291 bits, the accumulator <= 159 + 6291 + carry bits, so 112
// limbs is a PROVABLE upper bound (never approached on same-scale mesh data,
// where the spread is a handful of bits).  The active limb count is computed
// per call from the actual spread; the fixed storage just guarantees totality.
constexpr int kAccumLimbs = 112;
// |product of the cnt (nonzero) mantissas| -> mag[3]; returns the product sign.
inline int MulMag(const int64_t* pm, int cnt, uint64_t mag[3]) {
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
}
// acc[0..nLimbs) += sign * (mag << shift), two's complement.
inline void AddShiftedMag(uint64_t* acc, int nLimbs, const uint64_t mag[3],
                          int shift, int sign) {
  const int limbShift = shift / 64, bitShift = shift % 64;
  uint64_t tmp[kAccumLimbs] = {0};
  for (int i = 0; i < 3; ++i) {
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
// refusal): the true sum fits the two's-complement window by construction.
inline int SumSign(const Term* t, int nT) {
  if (nT == 0) return 0;
  long emin = t[0].e, emax = t[0].e;
  for (int i = 1; i < nT; ++i) {
    emin = std::min(emin, t[i].e);
    emax = std::max(emax, t[i].e);
  }
  int nLimbs = (int)((192 + (emax - emin)) / 64) + 3;
  if (nLimbs > kAccumLimbs) nLimbs = kAccumLimbs;  // never hit (proven bound)
  uint64_t acc[kAccumLimbs] = {0};
  for (int i = 0; i < nT; ++i)
    AddShiftedMag(acc, nLimbs, t[i].mag, (int)(t[i].e - emin), t[i].sign);
  if (acc[nLimbs - 1] >> 63) return -1;
  for (int i = 0; i < nLimbs; ++i)
    if (acc[i]) return 1;
  return 0;
}
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
// Exact orient3d sign (the e^0 coefficient), 0 iff the four points are exactly
// coplanar.  TOTAL: no window-fail; the 0 is a genuine geometric tie.
inline int ExactOrient3D(const double pts[4][3]) {
  int64_t M[4][3];
  int E[4][3];
  Decompose(pts, M, E);
  Term t[24];
  int nT = 0;
  for (const auto& s : kPerm) {
    int64_t pm[3];
    long pe = 0;
    int cnt = 0;
    bool zero = false;
    for (int r = 0; r < 4; ++r) {
      const int c = s[r];
      if (c == 3) continue;  // the ones column
      if (M[r][c] == 0) {
        zero = true;
        break;
      }
      pm[cnt++] = M[r][c];
      pe += E[r][c];
    }
    if (zero) continue;
    t[nT].sign = Parity(s) * MulMag(pm, cnt, t[nT].mag);
    t[nT].e = pe;
    ++nT;
  }
  return SumSign(t, nT);
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

// The micro exact tie-test (owner contract): the EXACT orient3d sign, 0 iff the
// four points are exactly coplanar.  ONE integer path (sos::ExactOrient3D),
// TOTAL for every finite-double input - no window-fail, no expansion fallback.
// RELUCTANT ACCEPTANCE: this exact kernel is net-new surface the design had
// declined ("no exact kernel in the tree", docs/Regularize3D.md B mechanism);
// the owner accepted it NARROWLY for the stage-6 tie residue - its ONLY call
// sites are behind a filter 0 (Orient3DSoS, the EdgePiercesTriSoS edge-in-plane
// guard, the test probe), never the certified fast path.
// QUEUED FOR REVISIT / TRIPWIRE (docs/Regularize3D.md open list): this integer
// path is legitimate ONLY as ONE predicate at ONE call site (~100 lines,
// exhaustively testable).  If a SECOND exact predicate or a SECOND call site is
// ever needed, VENDOR Shewchuk's public-domain predicates.c instead of growing
// this - do NOT rebuild expansion arithmetic piecemeal.
inline int Orient3DExactSign(const vec3& a, const vec3& b, const vec3& c,
                             const vec3& d) {
  const double pts[4][3] = {
      {a.x, a.y, a.z}, {b.x, b.y, b.z}, {c.x, c.y, c.z}, {d.x, d.y, d.z}};
  return sos::ExactOrient3D(pts);
}

// The complete orient3d decision (docs/Regularize3D.md stage 6): the certified
// filter sign on the fast path; else the micro exact tie-test decides the
// filter-uncertain-but-nonzero band exactly; else (a genuine exact zero) the
// single-global SoS breaks the tie.  NEVER 0.  `i*` are the four points'
// global vertex indices.
inline int Orient3DSoS(const vec3& a, const vec3& b, const vec3& c,
                       const vec3& d, int ia, int ib, int ic, int id) {
  const int s = Orient3DFilterSign(a, b, c, d);
  if (s != 0) return s;
  const int ex = Orient3DExactSign(a, b, c, d);
  if (ex != 0) return ex;
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
  std::vector<vec3> lo, hi;
  explicit TriSoup(const Manifold::Impl& in) {
    const int nTri = static_cast<int>(in.NumTri());
    tri.resize(nTri);
    vid.resize(nTri);
    lo.resize(nTri);
    hi.resize(nTri);
    for (int t = 0; t < nTri; ++t) {
      for (int k = 0; k < 3; ++k) {
        vid[t][k] = in.halfedge_.Start(3 * t + k);
        tri[t][k] = in.vertPos_[vid[t][k]];
      }
      lo[t] = la::min(la::min(tri[t][0], tri[t][1]), tri[t][2]);
      hi[t] = la::max(la::max(tri[t][0], tri[t][1]), tri[t][2]);
    }
  }
  bool BBoxOverlap(int i, int j) const {
    return !(hi[i].x < lo[j].x || hi[j].x < lo[i].x || hi[i].y < lo[j].y ||
             hi[j].y < lo[i].y || hi[i].z < lo[j].z || hi[j].z < lo[i].z);
  }
  bool SharesVert(int i, int j) const {
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (vid[i][a] == vid[j][b]) return true;
    return false;
  }
};

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
  const TriSoup soup(in);
  const auto& tri = soup.tri;
  for (int i = 0; i < nTri; ++i) {
    for (int j = i + 1; j < nTri; ++j) {
      if (!soup.BBoxOverlap(i, j)) continue;
      if (soup.SharesVert(i, j)) continue;  // self-adjacency skip (S4a)
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
  auto cr = [](const vec2& u, const vec2& v) { return u.x * v.y - u.y * v.x; };
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
  for (int i = 0; i < nTri; ++i)
    for (int j = i + 1; j < nTri; ++j) {
      if (!soup.BBoxOverlap(i, j) || soup.SharesVert(i, j)) continue;
      if (FacesFilterCoplanar(tri[i], tri[j]) &&
          TrianglesOverlap2D(tri[i], tri[j])) {
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
                             const std::vector<int>& face2cluster, double eps) {
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
          // global SoS (stage 6) DECIDES it (pierce / no-pierce).  A NON-
          // decision is only reachable with the convention disabled, and is the
          // documented fail-closed slot (A.boundaryTouch, RunCandidateB).
          const int sp =
              EdgePiercesTriSoS(u, w, T[0], T[1], T[2], A.vid[owner][e],
                                A.vid[owner][(e + 1) % 3], A.vid[tgt][0],
                                A.vid[tgt][1], A.vid[tgt][2]);
          if (sp == 1)
            hit = true;
          else if (sp != 0)
            A.boundaryTouch = true;
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
            const vec3 n = la::cross(T[1] - T[0], T[2] - T[0]);
            const double area2 = la::length2(n);
            if (area2 > 0.0) {
              const double margin = area2 * 1e-9;
              const double s0 = la::dot(n, la::cross(T[1] - T[0], P - T[0]));
              const double s1 = la::dot(n, la::cross(T[2] - T[1], P - T[1]));
              const double s2 = la::dot(n, la::cross(T[0] - T[2], P - T[2]));
              if (s0 >= -margin && s1 >= -margin && s2 >= -margin) hit = true;
            }
          }
        }
        if (hit && nPts < 4) {
          ptTri[nPts] = tgt;  // pierces tri tgt -> interior to tgt
          pts[nPts++] = pierce(A.vid[owner][e], A.vid[owner][(e + 1) % 3], tgt);
        }
      };
      for (int e = 0; e < 3; ++e) recordEdge(i, e, j);
      for (int e = 0; e < 3; ++e) recordEdge(j, e, i);
      // Dedup endpoints within the weld radius: a cap-plane reentrant junction
      // is pierced by BOTH walls' cap edges (the symmetric corner incidence)
      // and so is recorded twice; the emission weld would merge them anyway.
      // Genuine distinct seam endpoints are far more than eps apart, so an
      // ordinary seam is untouched (bitwise).
      for (int a = 0; a + 1 < nPts; ++a)
        for (int b = nPts - 1; b > a; --b)
          if (la::length(pts[a] - pts[b]) <= eps) {
            pts[b] = pts[nPts - 1];
            ptTri[b] = ptTri[nPts - 1];
            --nPts;
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

// Classify + emit the CLEAN (un-seamed, un-folded) faces.  Retention rule (SAME
// witness rule as the seamed path): keep iff the face's +n winding w_above == 0
// (the witness theorem makes every single face m=+1, so the solid is on the -n
// side of a retained face -> original orientation, flip-free).  A NEGATIVE
// w_above is exterior on both sides and DROPS (the axis-1 subtraction
// absorption), not a fail-closed.
//
// PER-FACE, not per-patch (reg3d-s7b/reg3d-arr).  Clean faces partition into
// clean-clean-connected patches; the coverage is CONSTANT across an uncrossed
// clean-clean edge UNLESS the arrangement is incomplete (a shares-vertex-skip
// crossing, the everted-corner defect).  The earlier flood decided a whole
// patch by ONE representative probe, ASSUMING uniformity - which holds on mild
// self-intersectors (siA/siB) but BREAKS on a folded soup: PokedCube's everted
// corner puts a w_above==0 boundary face and a w_above==-1 exterior face in ONE
// clean patch, so a w<0 representative wrongly DROPPED the boundary faces (a
// latent SILENT-WRONGNESS class - the mirror config would EMIT unverified
// faces), leaving the open fan the carrier fails on.
//
// So: probe each face by its OWN winding (this RETAINS the everted-corner
// boundary faces the flood dropped).  The coverage just above a CLEAN
// (uncrossed, non-near-coplanar post stage-5) triangle is CONSTANT across its
// interior - no face is crossed moving the query point at height eps over the
// triangle - so when the centroid probe GRAZES every seed (the
// component-local-seed axis, hit by axis-aligned integer geometry:
// BridgedCaps/TJunction/BarsCrossZ) we re-probe at OTHER interior points of the
// SAME triangle, which sample the SAME winding cell.  This dodges the graze
// WITHOUT any patch-uniformity assumption (fully sound: a wrong retain is
// impossible - every emitted face's winding was directly measured).  Only a
// face that grazes at EVERY interior point fails closed (never emit unverified
// geometry).  Cost: O(nTri) winding queries worst case (the doc's named
// winding-query perf axis; BVH is the later pass).
bool EmitCleanFaces(std::vector<OutTri3D>& out, const Manifold::Impl& in,
                    const BuildArrangement& A,
                    const std::vector<int>& face2cluster,
                    const std::vector<vec3>& seeds, double eps) {
  const int nTri = static_cast<int>(in.NumTri());
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
  for (int t = 0; t < nTri; ++t) {
    if (!isClean(t)) continue;
    const double nLen = la::length(A.faceN[t]);
    if (!(nLen > 0.0)) return false;
    const vec3 nHat = A.faceN[t] / nLen;
    std::optional<int> g;
    for (const auto& w : kBary) {
      const vec3 p =
          w[0] * A.tri[t][0] + w[1] * A.tri[t][1] + w[2] * A.tri[t][2];
      g = RobustWinding(in, p + eps * nHat, seeds);
      if (g) break;  // any interior sample measures the (constant) cell winding
    }
    if (!g) return false;  // grazes at every interior point (SoS): fail closed
    if (*g == 0) out.push_back({A.tri[t][0], A.tri[t][1], A.tri[t][2]});
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

// NEAR-COPLANAR WIDEN + GLOBAL-PLANARITY GUARD (docs/Regularize3D.md stage-5;
// reg3d-nearcoplanar-research candidate (a)).  The exact-coplanar fold only
// admits pairs whose six cross orient3d filter signs are all 0 (coplanarity gap
// below ~1 ULP).  Faces within eps of coplanar but ABOVE that bound are
// DECIDABLE non-coplanar yet their arrangement cells are sub-eps thin:
// enumerated transversally they double-round to slivers (unresolvable sheet
// contact) - the thin band this pass closes.
//
// This is an INPUT-SIDE PLANARIZATION, run BEFORE B's enumeration/winding/emit,
// so B RE-DERIVES the whole arrangement from the snapped input (the thin cell
// ceases to exist).  It is NOT an emission-time snap (those fight decisions the
// arrangement already made, ExactArrangement3D variant-E kill); it perturbs the
// INPUT by <= eps inside the standing epsilon-valid contract, coordinated by
// construction.
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
//     fold in RunCandidateBBuild handles it verbatim and the m == winding-jump
//     self-check holds by the exact argument.
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
    vec3 ref;
    if (!unitN(faces[0], ref)) {
      return StageResult<Manifold::Impl>::Fatal(
          FatalReason::DirtyComponentUnresolved,
          "candidate B: near-coplanar cluster has a degenerate face - "
          "fail-closed");
    }
    vec3 nSum(0.0, 0.0, 0.0);
    for (int f : faces) {
      const vec3 raw = la::cross(tri[f][1] - tri[f][0], tri[f][2] - tri[f][0]);
      nSum += (la::dot(raw, ref) < 0.0) ? -raw : raw;
    }
    const double nl = la::length(nSum);
    if (!(nl > 0.0)) {
      return StageResult<Manifold::Impl>::Fatal(
          FatalReason::DirtyComponentUnresolved,
          "candidate B: near-coplanar cluster normal is degenerate - "
          "fail-closed");
    }
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
            "candidate B: near-coplanar cluster fails the global-planarity "
            "guard (curved chain, max deviation > eps) - fail-closed");
      }
    // SNAP onto the fitted plane; flag an inconsistent multi-cluster vertex.
    for (int v : verts) {
      if (vertCluster[v] >= 0 && vertCluster[v] != clusterId) {
        return StageResult<Manifold::Impl>::Fatal(
            FatalReason::DirtyComponentUnresolved,
            "candidate B: a vertex lies in two near-coplanar clusters "
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
  // NEAR-COPLANAR PRE-PASS (docs/Regularize3D.md stage-5): planarize any
  // within-eps near-coplanar overlap cluster onto its fitted plane so B
  // re-derives the arrangement from an exactly-coplanar input.  No near cluster
  // -> `dirty` is used bitwise (the exact path is unperturbed); a curved chain
  // (global-planarity guard failure) fails closed here, distinctly named.
  StageResult<Manifold::Impl> snapped = SnapNearCoplanarClusters(dirty, eps);
  if (snapped.fatal) return snapped;
  const Manifold::Impl& in = snapped.value ? *snapped.value : dirty;

  // Exactly-coplanar face clusters are resolved by the in-plane fold; the
  // transversal seam enumeration skips their pairs so their exact-zero coplanar
  // ties do not raise boundaryTouch.
  const std::vector<int> face2cluster = DetectCoplanarClusters(in);
  const BuildArrangement A = RecordSeams(in, face2cluster, eps);
  if (A.boundaryTouch) {
    // A deciding pierce predicate hit an exact-zero / filter-uncertain boundary
    // that the coplanar fold does NOT consume: a NON-coplanar vertex-on-face /
    // edge-in-face incidence (the residual single-global SoS tie family,
    // PokedCube-class), or a near-coplanar tie the stage-5 pre-pass did not
    // fold (no 2D overlap, or a guard-refused curved chain).  Guessing a sign
    // would risk an oracle-wrong resolve.  Fail closed.
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
  return RunCandidateBBuild(in, A, face2cluster, eps);
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

CleanFaceProbe RegularizeCleanFaces_Probe(const Manifold::Impl& soup) {
  CleanFaceProbe out;
  double eps = EpsilonFromScale(soup.bBox_.Scale(), 1000);
  if (!(eps > 0.0) || !std::isfinite(eps)) return out;
  // Mirror RunCandidateB's prefix so the classify inputs match production
  // exactly (snap near-coplanar, detect exact-coplanar clusters, record seams).
  StageResult<Manifold::Impl> snapped = SnapNearCoplanarClusters(soup, eps);
  if (snapped.fatal) return out;
  const Manifold::Impl& in = snapped.value ? *snapped.value : soup;
  const std::vector<int> face2cluster = DetectCoplanarClusters(in);
  const BuildArrangement A = RecordSeams(in, face2cluster, eps);
  // Same winding seeds as RunCandidateBBuild.
  const vec3 c = in.bBox_.Center();
  const double L = in.bBox_.Scale() + 1.0;
  const std::vector<vec3> seeds = {
      c + L * vec3(3.13, 5.71, 1.37),   c + L * vec3(-2.71, 1.41, 4.19),
      c + L * vec3(1.73, -3.31, -2.23), c + L * vec3(-4.27, -1.19, 2.83),
      c + L * vec3(2.39, -4.61, 3.07),  c + L * vec3(-1.51, 3.89, -4.43)};
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
