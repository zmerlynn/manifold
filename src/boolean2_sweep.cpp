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
//
// Sweep-line arrangement and winding for the Boolean2 pipeline, after J.R.
// Smith, "Robustness in geometric algorithms" (2009), ch. 6-8. SweepWinding
// runs two lexicographic sweeps: an arrangement pass that discovers crossings
// and emits a true arrangement, then a winding pass that emits the retained
// boundary. Smith's 7.6.2 block rule drives both: at each event p it splits
// every status edge bracketing p through p, collapsing a dense near-concurrence
// to one shared vertex and forcing any sub-eps residual crossing onto it.
//
// Engine extension (3D overlap-removal): PolySet2's value carries {m, srcId}
// so that each piece in the winding pass can be attributed to its source face.
// MergeVerticals1D performs per-interval active-contributor tracking (7.6.1
// footnote 9 rewrite): each emitted vertical interval carries the id of its
// unique contributor; multi-source intervals record a conflict (srcId=-1, count
// incremented). SweepWinding accepts an optional capture out-channel and an
// optional conflict-count out-pointer; default-null makes existing 2D callers
// unchanged.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "boolean2.h"
#include "manifold/optional_assert.h"
#include "shared.h"

namespace manifold {

namespace {

// Lexicographic (x then y) order on exact vec2 coordinates.
struct LexLess {
  bool operator()(const vec2& a, const vec2& b) const {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  }
};
constexpr LexLess kLexLess{};

// Table 7.1 (p.107): y-value of a non-vertical segment at abscissa xP, via the
// shared interpolation kernel (shared.h; lambda measured from the nearer
// endpoint, 8.5).
double YAtX(const vec2& l, const vec2& r, double xP) {
  return Interpolate(vec3(l, 0.0), vec3(r, 0.0), xP).x;
}

// Table 7.2 (p.108): the crossing point on the common x-range, via the shared
// intersection kernel (shared.h). Contract: A=(xL,yAL)->(xR,yAR),
// B=(xL,yBL)->(xR,yBR) with yAL <= yBL and yBR <= yAR, at least one strict; the
// result lies in both bounding boxes (8.5 eqs 8.18/8.19).
vec2 Interpolate2(double xL, double yAL, double yBL, double xR, double yAR,
                  double yBR) {
  const vec4 q = Intersect(vec3(xL, yAL, 0.0), vec3(xR, yAR, 0.0),
                           vec3(xL, yBL, 0.0), vec3(xR, yBR, 0.0));
  return {q.x, q.y};
}

// 6.5 vertex-on-edge-interior test, exact equality only.
bool OnInterior(const vec2& v, const vec2& a, const vec2& b) {
  const vec2& lo = kLexLess(a, b) ? a : b;
  const vec2& hi = kLexLess(a, b) ? b : a;
  if (v == lo || v == hi) return false;
  if (lo.x == hi.x) {  // vertical: xP==xL and yL < yP < yR
    return v.x == lo.x && lo.y < v.y && v.y < hi.y;
  }
  if (v.x < lo.x || v.x > hi.x) return false;
  return YAtX(lo, hi, v.x) == v.y;
}

enum class Side : uint8_t { UNDER, OVER, ON, ENDS };

// Positive fill (Add: w > 0), intersection (Intersect: w > 1). Subtract feeds
// the second operand with negative multiplicity, so it uses the Add rule.
bool IsInside(WindRule rule, int64_t w) {
  switch (rule) {
    case WindRule::Add:
      return w > 0;
    case WindRule::Intersect:
      return w > 1;
  }
  return false;
}

// Multiplicity + source-id pair. The PolySet2 mapped value (engine extension):
// m is the signed net multiplicity in lex-normalized direction; srcId is the
// face id from the 3D caller (0 = unattributed, -1 = conflicted).
struct PolyVal {
  int64_t m = 0;
  int32_t srcId = 0;
};

// Lex-normalized directed-edge key with PolyVal: the PolySet2 of Smith
// 7.1.2 / 6.5 (reversal negates m, srcId unchanged; exact-coincident edges
// sum m and merge srcIds per the conflict rule; zero-mult and zero-length
// annihilate). Keyed by (lo, hi) in lexicographic point order.
struct PairLexLess {
  bool operator()(const std::pair<vec2, vec2>& a,
                  const std::pair<vec2, vec2>& b) const {
    if (kLexLess(a.first, b.first)) return true;
    if (kLexLess(b.first, a.first)) return false;
    return kLexLess(a.second, b.second);
  }
};
using PolySet2 = std::map<std::pair<vec2, vec2>, PolyVal, PairLexLess>;

// Merge two srcIds under the conflict rule: same id or one is already -1 ->
// propagate (no new conflict); two distinct valid ids -> conflict (-1), return
// true to signal a new conflict event.
static bool MergeSrcId(int32_t& existing, int32_t incoming) {
  if (existing == incoming) return false;
  if (existing == -1 || incoming == -1) {
    existing = -1;
    return false;  // conflict already recorded upstream
  }
  existing = -1;
  return true;  // new conflict
}

// Add the directed edge a->b with PolyVal pv into ps. Normalizes to lex order
// (negating m on reversal, srcId unchanged). On an exact-coincident key: sums
// m (erasing if zero); updates srcId via MergeSrcId; returns true if a new
// conflict was recorded.
bool PolySetAdd(PolySet2& ps, vec2 a, vec2 b, PolyVal pv) {
  if (a == b || pv.m == 0) return false;
  if (kLexLess(b, a)) {
    std::swap(a, b);
    pv.m = -pv.m;
  }
  auto it = ps.find({a, b});
  if (it == ps.end()) {
    ps.emplace(std::make_pair(a, b), pv);
    return false;
  }
  const int64_t newM = it->second.m + pv.m;
  if (newM == 0) {
    ps.erase(it);
    return false;
  }
  const bool conflict = MergeSrcId(it->second.srcId, pv.srcId);
  it->second.m = newM;
  return conflict;
}

// 7.6.1 footnote 9 (rewritten for per-interval contributor tracking): resolve
// each x-group of coincident/overlapping vertical edges into signed-coverage
// segments between consecutive breakpoints. Each emitted interval carries the
// srcId of its unique active contributor; if multiple distinct srcIds are
// active in an interval, srcId = -1 and *conflicts is incremented (if
// non-null). Adjacent same-plane triangles making multi-source vertical groups
// is ordinary in the 3D use-case, not a corner case. All input breakpoints are
// preserved.
void MergeVerticals1D(PolySet2& ps, int* conflicts = nullptr) {
  using SegEntry = std::tuple<double, double, int64_t, int32_t>;
  std::map<double, std::vector<SegEntry>> groups;
  std::vector<std::pair<vec2, vec2>> toErase;
  for (const auto& kv : ps) {
    if (kv.first.first.x == kv.first.second.x) {
      // PolySet2 lex-normalizes so kv.first.first.y <= kv.first.second.y
      groups[kv.first.first.x].emplace_back(kv.first.first.y, kv.first.second.y,
                                            kv.second.m, kv.second.srcId);
      toErase.push_back(kv.first);
    }
  }
  for (const auto& k : toErase) ps.erase(k);

  for (const auto& g : groups) {
    const double x = g.first;
    // Build delta events: at each y-breakpoint record (delta_m, srcId)
    std::map<double, std::vector<std::pair<int64_t, int32_t>>> events;
    for (const auto& seg : g.second) {
      double yLo, yHi;
      int64_t m;
      int32_t srcId;
      std::tie(yLo, yHi, m, srcId) = seg;
      events[yLo].emplace_back(m, srcId);
      events[yHi].emplace_back(-m, srcId);
    }
    // Active contributors: srcId -> net coverage contribution
    std::map<int32_t, int64_t> activeSrcs;
    int64_t totalCover = 0;
    double prevY = 0.0;
    bool have = false;
    for (const auto& ev : events) {
      if (have && totalCover != 0) {
        // Determine srcId for this (prevY, ev.first) interval.
        int32_t resolvedId = 0;
        bool have_id = false;
        bool multi = false;
        for (const auto& as : activeSrcs) {
          if (as.second != 0) {
            if (!have_id) {
              resolvedId = as.first;
              have_id = true;
            } else {
              multi = true;
              break;
            }
          }
        }
        if (multi) {
          resolvedId = -1;
          if (conflicts) ++(*conflicts);
        } else {
          DEBUG_ASSERT(
              have_id, logicErr,
              "MergeVerticals1D: nonzero totalCover but no active src");
        }
        PolySetAdd(ps, {x, prevY}, {x, ev.first}, {totalCover, resolvedId});
      }
      // Apply the delta events at ev.first
      for (const auto& e : ev.second) {
        activeSrcs[e.second] += e.first;
        if (activeSrcs[e.second] == 0) activeSrcs.erase(e.second);
        totalCover += e.first;
      }
      prevY = ev.first;
      have = true;
    }
  }
}

// A live sweep edge. Roles (7.6.1): l = processed end, r = pending end. `m` is
// the winding increment crossing l->r; LexMultiplicity normalizes it to lex
// direction. `srcId` is the face attribution carried through splits.
struct SweepEdge {
  vec2 l, r;
  int64_t m;
  int32_t srcId;
  uint64_t seq;
};

// One lexicographic sweep. Arrangement mode discovers crossings (TestPair) and
// emits every finalized piece with its signed multiplicity. Winding mode runs
// over the finished arrangement without crossing discovery and emits a piece
// iff the fill (IsInside) differs across it, oriented interior-on-left; the
// 7.6.2 forced-through block rule still splits bracketed edges at the event.
enum class SweepMode { Arrangement, Winding };

class SweepPass {
 public:
  SweepPass(WindRule rule, SweepMode mode,
            std::vector<SweepCapture>* capture = nullptr)
      : rule_(rule), mode_(mode), capture_(capture) {}

  void Seed(vec2 a, vec2 b, PolyVal pv) { PendingAdd(a, b, pv); }

  // `mergeVerticalOutput` runs the footnote-9 vertical resolve on the result.
  void Run(bool mergeVerticalOutput = true) {
    while (!events_.empty()) {
      const vec2 p = *events_.begin();
      events_.erase(events_.begin());
      ProcessEvent(p);
    }
    DEBUG_ASSERT(status_.empty() && pending_.empty(), logicErr,
                 "Boolean2 sweep left live edges after draining");
    if (mergeVerticalOutput) MergeVerticals1D(out_, &conflictCount_);
  }

  PolySet2& Out() { return out_; }
  int ConflictCount() const { return conflictCount_; }

 private:
  int64_t LexMultiplicity(const SweepEdge& e) const {
    return kLexLess(e.l, e.r) ? e.m : -e.m;
  }

  void PendingAdd(vec2 a, vec2 b, PolyVal pv) {
    if (a == b || pv.m == 0) return;
    if (kLexLess(b, a)) {
      std::swap(a, b);
      pv.m = -pv.m;
    }
    auto& inner = pending_[a];
    auto it = inner.find(b);
    if (it == inner.end()) {
      inner.emplace(b, pv);
    } else {
      const int64_t newM = it->second.m + pv.m;
      if (newM == 0) {
        inner.erase(it);
      } else {
        if (MergeSrcId(it->second.srcId, pv.srcId)) ++conflictCount_;
        it->second.m = newM;
      }
    }
    if (inner.empty()) pending_.erase(a);
    events_.insert(a);
    events_.insert(b);
  }

  Side Classify(const SweepEdge& e, const vec2& p) const {
    if (e.r == p) return Side::ENDS;
    if (e.l.x == e.r.x) {
      if (p.x != e.l.x) return Side::ON;
      if (e.l == p) return (e.r.y > p.y) ? Side::OVER : Side::UNDER;
      const double ylo = std::min(e.l.y, e.r.y), yhi = std::max(e.l.y, e.r.y);
      if (yhi <= p.y) return Side::UNDER;
      if (ylo >= p.y) return Side::OVER;
      return Side::ON;
    }
    const vec2& lo = kLexLess(e.l, e.r) ? e.l : e.r;
    const vec2& hi = kLexLess(e.l, e.r) ? e.r : e.l;
    if (p.x < lo.x || p.x > hi.x) return p.x < lo.x ? Side::OVER : Side::ON;
    const double ys = YAtX(lo, hi, p.x);
    if (ys < p.y) return Side::UNDER;
    if (ys > p.y) return Side::OVER;
    return Side::ON;
  }

  static int GradientRank(const SweepEdge& e) {
    if (e.r.x == e.l.x) return e.r.y > e.l.y ? 1 : -1;
    return 0;
  }
  static bool GradientLess(const SweepEdge& a, const SweepEdge& b) {
    const int ra = GradientRank(a), rb = GradientRank(b);
    if (ra != rb) return ra < rb;
    if (ra != 0) return a.seq < b.seq;
    const double c = la::cross(a.r - a.l, b.r - b.l);
    if (c != 0) return c > 0;
    return a.seq < b.seq;
  }

  // 7.4.3 (+collect override): in arrangement mode the piece goes to `out_`
  // with its role PolyVal. In winding mode emit iff the fill (IsInside) differs
  // across it, oriented interior-on-left; capture into `capture_` if non-null.
  void EmitBoundary(const vec2& from, const vec2& to, int64_t m, int64_t below,
                    int64_t above, int32_t srcId) {
    if (from == to) return;
    if (mode_ == SweepMode::Arrangement) {
      if (PolySetAdd(out_, from, to, {m, srcId})) ++conflictCount_;
      return;
    }
    const bool insB = IsInside(rule_, below);
    const bool insA = IsInside(rule_, above);
    if (insB == insA) return;
    // Retained-piece out-channel: only pieces crossing the fill boundary are
    // captured. Store EMISSION-ORIENTED (interior-on-left of from->to).
    // insA: above side is inside, edge goes in the original sweep direction.
    // insB: below side is inside, edge goes in the reversed direction.
    if (capture_) {
      const vec2 emitFrom = insA ? from : to;
      const vec2 emitTo = insA ? to : from;
      capture_->push_back({emitFrom, emitTo, srcId});
    }
    // above inside -> lex-forward; below inside -> lex-backward.
    if (PolySetAdd(out_, from, to,
                   {insA ? (kLexLess(from, to) ? 1 : -1)
                         : (kLexLess(from, to) ? -1 : 1),
                    srcId}))
      ++conflictCount_;
  }

  // 7.4.1: split status edge idx at q; shorten in place, requeue the remainder.
  void SplitAt(size_t idx, const vec2& q) {
    SweepEdge e = status_[idx];
    if (q == e.r) return;
    if (q == e.l) {
      status_.erase(status_.begin() + idx);
      PendingAdd(q, e.r, {e.m, e.srcId});
      return;
    }
    status_[idx].r = q;
    PendingAdd(q, e.r, {e.m, e.srcId});
    events_.insert(q);
  }

  // 7.4.2 preliminary + 7.1.1/8.5 crossing test on a neighbour pair.
  void TestPair(size_t i, size_t j) {
    if (j >= status_.size() || i >= j) return;
    SweepEdge a = status_[i], b = status_[j];
    if (a.l == b.l || a.l == b.r || a.r == b.l || a.r == b.r) return;
    if (OnInterior(b.r, a.l, a.r)) return SplitAt(i, b.r);
    if (OnInterior(a.r, b.l, b.r)) return SplitAt(j, a.r);
    const double axlo = std::min(a.l.x, a.r.x), axhi = std::max(a.l.x, a.r.x);
    const double bxlo = std::min(b.l.x, b.r.x), bxhi = std::max(b.l.x, b.r.x);
    const double xL = std::max(axlo, bxlo), xR = std::min(axhi, bxhi);
    if (xL > xR) return;
    double yAL, yAR, yBL, yBR;
    if (axlo == axhi) {
      yAL = std::min(a.l.y, a.r.y);
      yAR = std::max(a.l.y, a.r.y);
    } else {
      const vec2& lo = kLexLess(a.l, a.r) ? a.l : a.r;
      const vec2& hi = kLexLess(a.l, a.r) ? a.r : a.l;
      yAL = YAtX(lo, hi, xL);
      yAR = YAtX(lo, hi, xR);
    }
    if (bxlo == bxhi) {
      yBL = std::min(b.l.y, b.r.y);
      yBR = std::max(b.l.y, b.r.y);
    } else {
      const vec2& lo = kLexLess(b.l, b.r) ? b.l : b.r;
      const vec2& hi = kLexLess(b.l, b.r) ? b.r : b.l;
      yBL = YAtX(lo, hi, xL);
      yBR = YAtX(lo, hi, xR);
    }
    const bool aLower =
        (yAL <= yBL) && (yBR <= yAR) && (yAL < yBL || yBR < yAR);
    const bool bLower =
        (yBL <= yAL) && (yAR <= yBR) && (yBL < yAL || yAR < yBR);
    if (!aLower && !bLower) return;
    const vec2 q = aLower ? Interpolate2(xL, yAL, yBL, xR, yAR, yBR)
                          : Interpolate2(xL, yBL, yAL, xR, yBR, yAR);
    events_.insert(q);
    SplitAt(j, q);  // j first: index safety
    SplitAt(i, q);
  }

  // 7.6.2 block rule: bracket the contiguous status block straddling p (from
  // the last strictly-under edge to the first strictly-over edge), split each
  // bracketed edge through p, and re-enter the remainders sorted by gradient.
  // In arrangement mode, test the newly-adjacent pairs for crossings.
  void ProcessEvent(const vec2& p) {
    const size_t n = status_.size();
    std::vector<Side> cls(n);
    for (size_t i = 0; i < n; ++i) cls[i] = Classify(status_[i], p);
    size_t lo = 0;
    while (lo < n && cls[lo] == Side::UNDER) lo++;
    size_t hi = n;
    while (hi > lo && cls[hi - 1] == Side::OVER) hi--;
    std::vector<SweepEdge> reinsert;
    int64_t w = 0;
    for (size_t i = 0; i < lo; ++i) w += LexMultiplicity(status_[i]);
    for (size_t i = lo; i < hi; ++i) {
      SweepEdge& e = status_[i];
      const int64_t lm = LexMultiplicity(e);
      const int64_t below = w, above = w + lm;
      w = above;
      if (cls[i] == Side::ENDS) {
        EmitBoundary(e.l, e.r, e.m, below, above, e.srcId);
      } else {  // forced through p (7.6.2): commit [e.l, p], re-enter [p, e.r]
        if (e.l != p) EmitBoundary(e.l, p, e.m, below, above, e.srcId);
        reinsert.push_back({p, e.r, e.m, e.srcId, seqCounter_++});
      }
    }
    const bool removedAny = hi > lo;
    status_.erase(status_.begin() + lo, status_.begin() + hi);
    auto pit = pending_.find(p);
    if (pit != pending_.end()) {
      for (const auto& kv : pit->second)
        reinsert.push_back(
            {p, kv.first, kv.second.m, kv.second.srcId, seqCounter_++});
      pending_.erase(pit);
    }
    std::stable_sort(reinsert.begin(), reinsert.end(), GradientLess);
    status_.insert(status_.begin() + lo, reinsert.begin(), reinsert.end());
    const size_t k = reinsert.size();
    // 8.4: test the pairs that just became adjacent. The arrangement is static
    // in winding mode, so no adjacency tests and no constructed points.
    if (mode_ == SweepMode::Arrangement) {
      if (k > 0) {
        TestPair(lo + k - 1, lo + k);
        if (lo > 0) TestPair(lo - 1, lo);
      } else if (removedAny && lo > 0) {
        TestPair(lo - 1, lo);
      }
    }
  }

  WindRule rule_;
  SweepMode mode_;
  std::vector<SweepCapture>* capture_;
  std::set<vec2, LexLess> events_;
  std::map<vec2, std::map<vec2, PolyVal, LexLess>, LexLess> pending_;
  std::vector<SweepEdge> status_;
  PolySet2 out_;
  uint64_t seqCounter_ = 0;
  int conflictCount_ = 0;
};

// Arrangement pass: the collect sweep discovers crossings (TestPair), resolves
// near-concurrences with the block rule, merges coincident verticals, and emits
// every finalized piece. The result is a true arrangement - every crossing is a
// shared vertex, no two pieces cross.
PolySet2 CollectArrangement(PolySet2 arr, WindRule rule,
                            int* conflicts = nullptr) {
  SweepPass collect(rule, SweepMode::Arrangement);
  for (const auto& kv : arr)
    collect.Seed(kv.first.first, kv.first.second, kv.second);
  collect.Run();
  if (conflicts) *conflicts += collect.ConflictCount();
  return std::move(collect.Out());
}

// Arrangement then winding over the resulting true arrangement. The optional
// `capture` receives one SweepCapture per retained boundary piece in the
// winding pass. `conflicts` is incremented by the total number of attribution
// conflicts. When `negOut` is non-null, a SECOND winding measure runs over the
// same collected arrangement with all multiplicities negated and *negOut
// receives its retained boundary - the other signed difference from the one
// arrangement. The winding pass constructs no points, so both measures share
// the arrangement's geometry exactly. Its conflict events are sign-independent
// duplicates of the first measure's and are not re-counted.
PolySet2 CollectThenMeasure(PolySet2 arr, WindRule rule,
                            std::vector<SweepCapture>* capture = nullptr,
                            int* conflicts = nullptr,
                            PolySet2* negOut = nullptr) {
  PolySet2 clean = CollectArrangement(std::move(arr), rule, conflicts);
  SweepPass measure(rule, SweepMode::Winding, capture);
  for (const auto& kv : clean)
    measure.Seed(kv.first.first, kv.first.second, kv.second);
  measure.Run(/*mergeVerticalOutput=*/false);  // preserve T-junction verts
  if (conflicts) *conflicts += measure.ConflictCount();
  if (negOut) {
    SweepPass measureNeg(rule, SweepMode::Winding, nullptr);
    for (const auto& kv : clean)
      measureNeg.Seed(kv.first.first, kv.first.second,
                      {-kv.second.m, kv.second.srcId});
    measureNeg.Run(/*mergeVerticalOutput=*/false);
    *negOut = std::move(measureNeg.Out());
  }
  return std::move(measure.Out());
}

}  // namespace

// Arranges `edges` (referencing `verts`) and returns the retained boundary as
// directed OutEdges under `rule`. `verts` is in/out: the constructed crossing
// vertices the returned edges reference are appended to it.
// If `capture` is non-null, appends one SweepCapture per retained winding-pass
// piece (see SweepCapture definition in boolean2.h).
// If `conflictCount` is non-null, *conflictCount receives the total number of
// source-id attribution conflicts.
// If `negEdges` is non-null, *negEdges receives the retained boundary of the
// negated measure over the same arrangement, referencing the same `verts`
// (see boolean2.h). 2D callers pass nullptr for all three.
std::vector<OutEdge> SweepWinding(const std::vector<EdgeM>& edges,
                                  std::vector<vec2>& verts, WindRule rule,
                                  std::vector<SweepCapture>* capture,
                                  int* conflictCount,
                                  std::vector<OutEdge>* negEdges) {
  std::map<std::pair<double, double>, int> vertId;
  for (int v = 0; v < static_cast<int>(verts.size()); ++v)
    vertId.emplace(std::make_pair(verts[v].x, verts[v].y), v);
  auto getId = [&](const vec2& p) -> int {
    const auto key = std::make_pair(p.x, p.y);
    auto it = vertId.find(key);
    if (it != vertId.end()) return it->second;
    const int id = static_cast<int>(verts.size());
    verts.push_back(p);
    vertId.emplace(key, id);
    return id;
  };

  // Seed the sweep from the input edges as a lex-normalized signed-multiplicity
  // PolySet2 (coincident edges sum, opposite ones cancel). The sweep discovers
  // the exact crossings and applies the 7.6.2 block rule, collapsing dense
  // near-concurrences to a single shared vertex; eps-scale vertex-on-edge
  // incidences are pre-split by the caller.
  PolySet2 arr;
  int conflicts = 0;
  for (const auto& e : edges) {
    if (e.v0 == e.v1) continue;
    if (PolySetAdd(arr, verts[e.v0], verts[e.v1], {e.mult, e.srcId}))
      ++conflicts;
  }
  MergeVerticals1D(arr, &conflicts);

  PolySet2 negOut;
  const PolySet2 out = CollectThenMeasure(
      std::move(arr), rule, capture, &conflicts, negEdges ? &negOut : nullptr);
  if (conflictCount) *conflictCount = conflicts;

  // Materialize retained boundaries as directed OutEdges via the shared getId,
  // so both measures reference one vert list. A key stores the lex-min/lex-max
  // endpoints with signed multiplicity: positive runs lo->hi.
  auto materialize = [&](const PolySet2& ps, std::vector<OutEdge>& result) {
    for (const auto& kv : ps) {
      const int64_t m = kv.second.m;
      if (m == 0) continue;
      const int loId = getId(kv.first.first);
      const int hiId = getId(kv.first.second);
      const int from = m > 0 ? loId : hiId;
      const int to = m > 0 ? hiId : loId;
      for (int64_t c = 0; c < std::abs(m); ++c)
        result.push_back({from, to, 1, kv.second.srcId});
    }
  };
  std::vector<OutEdge> result;
  materialize(out, result);
  if (negEdges) materialize(negOut, *negEdges);
  return result;
}

}  // namespace manifold
