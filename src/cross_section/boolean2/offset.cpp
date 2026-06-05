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
// Polygon offset backing `CrossSection::Offset`.

#include "offset.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "boolean2.h"
#include "predicates.h"

namespace manifold {
namespace boolean2 {

namespace {

// Outward normal of a directed edge (right-perpendicular, unit length).
// For a CCW polygon, this points away from the interior.
vec2 OutwardNormal(vec2 edge) {
  const double len = std::sqrt(edge.x * edge.x + edge.y * edge.y);
  // Keep this exact: public CrossSection inputs are regularized before
  // Offset, and a local eps threshold here can erase valid tiny-edge corners.
  if (len == 0) return vec2(0, 0);
  return vec2(edge.y / len, -edge.x / len);
}

vec2 RotateDegrees(vec2 v, double angle) {
  const double c = cosd(angle);
  const double s = sind(angle);
  return vec2(v.x * c - v.y * s, v.x * s + v.y * c);
}

bool StraightTurn(vec2 ePrev, vec2 eNext) {
  const double prevLen2 = dot(ePrev, ePrev);
  const double nextLen2 = dot(eNext, eNext);
  const double maxLen2 = std::max(prevLen2, nextLen2);
  if (maxLen2 == 0) return true;
  const double eps = EpsilonFromScale(std::sqrt(maxLen2));
  const double cross = la::cross(ePrev, eNext);
  return 4.0 * cross * cross <= maxLen2 * eps * eps;
}

// Number of chords in a full circle at radius `r` such that each chord's
// perpendicular sagitta error stays <= arcTol. The public Offset() path derives
// arcTol from a requested segment count via the same deterministic
// math::cos(kPi/n) sagitta formula `withinTol` uses below, so the round-trip is
// exact; use the analytic inverse as a search bound, then evaluate the formula
// directly to recover the minimal integer count despite roundoff.
int FullCircleChordCount(double r, double arcTol) {
  if (!std::isfinite(r) || r <= 0) return 1;
  if (!std::isfinite(arcTol) || arcTol <= 0) {
    return Quality::GetCircularSegments(r);
  }
  if (arcTol >= 2.0 * r) return 1;
  auto withinTol = [&](int n) {
    return (1.0 - math::cos(kPi / n)) * r <= arcTol;
  };
  const double cosHalfStep = std::clamp(1.0 - arcTol / r, -1.0, 1.0);
  const double halfStep = std::acos(cosHalfStep);
  if (!std::isfinite(halfStep) || halfStep <= 0) {
    return Quality::GetCircularSegments(r);
  }
  const double estimate = std::ceil(kPi / halfStep);
  if (!std::isfinite(estimate) || estimate > std::numeric_limits<int>::max()) {
    return std::numeric_limits<int>::max();
  }
  int hi = std::max(2, static_cast<int>(estimate));
  while (!withinTol(hi)) {
    if (hi > std::numeric_limits<int>::max() / 2) {
      return std::numeric_limits<int>::max();
    }
    hi *= 2;
  }
  int lo = 1;
  while (lo + 1 < hi) {
    const int mid = lo + (hi - lo) / 2;
    if (withinTol(mid)) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return hi;
}

bool BeforeSweepTarget(vec2 dir, vec2 target, double rotSign) {
  const double cross = la::cross(dir, target);
  return rotSign > 0 ? cross > 0 : cross < 0;
}

// Append a round (arc) join between the offset endpoints `endPrev`
// (start of arc) and `startNext` (end of arc), centered at `V` with
// radius |delta|. Emits deterministic full-circle chord steps until the
// next step would reach or pass startNext; caller appends startNext.
void AppendRoundJoin(SimplePolygon& out, vec2 V, vec2 nPrev, vec2 nNext,
                     double delta, double arcTol) {
  // The arc sweeps the short side between nPrev and nNext on a circle
  // of radius |delta| around V. Caller invokes only on the convex
  // side, where this short side is < pi. The step direction is sign(delta)
  // so inset sweeps the mirror direction.
  const double absDelta = std::fabs(delta);
  const int fullCircleCount = FullCircleChordCount(absDelta, arcTol);
  const double step = 360.0 / fullCircleCount;
  const double rotSign = (delta >= 0) ? 1.0 : -1.0;
  for (int i = 1; i < fullCircleCount; ++i) {
    const vec2 dir = RotateDegrees(nPrev, rotSign * i * step);
    if (!BeforeSweepTarget(dir, nNext, rotSign)) break;
    out.push_back(V + delta * dir);
  }
}

// Square join: two extra vertices forming a chord tangent to a circle
// of radius |delta| around V at the bisector, capping the corner with a
// flat. Matches Clipper2's `DoSquare`. For a corner with half-angle
// alpha (angle between bisector and either normal), the chord half-
// length is |delta| * tan(alpha/2) = |delta| * (1 - cos alpha) /
// sin alpha; using the more FP-stable form sin / (1 + cos).
void AppendSquareJoin(SimplePolygon& out, vec2 V, vec2 nPrev, vec2 nNext,
                      double delta) {
  vec2 bis(nPrev.x + nNext.x, nPrev.y + nNext.y);
  const double bisLen = std::sqrt(bis.x * bis.x + bis.y * bis.y);
  if (bisLen == 0) return;  // 180-deg reversal; nothing reasonable to emit
  bis.x /= bisLen;
  bis.y /= bisLen;
  const vec2 tang(-bis.y, bis.x);
  const double cosHalf = bis.x * nPrev.x + bis.y * nPrev.y;
  if (cosHalf <= 0) return;  // reflex; convex caller should not invoke this
  const double sinHalf = std::sqrt(std::max(0.0, 1.0 - cosHalf * cosHalf));
  const double half = std::fabs(delta) * sinHalf / (1.0 + cosHalf);
  // Signed delta: the cap follows the offset side (outward for delta > 0,
  // inward for inset). Using std::fabs here would mirror the cap to the wrong
  // side for inset, which the Positive union cannot recover.
  const vec2 mid = V + delta * bis;
  out.push_back(mid - half * tang);
  out.push_back(mid + half * tang);
}

// Intersect the two offset edges (prev edge's offset line and next edge's
// offset line) at the vertex V. Returns the miter point. Falls back to V
// + delta * average-normal if lines are parallel.
vec2 MiterPoint(vec2 V, vec2 nPrev, vec2 nNext, double delta) {
  // The two offset lines pass through V + delta*nPrev (perp to ePrev)
  // and V + delta*nNext (perp to eNext). Their intersection lies along
  // the bisector at distance delta / cos(half-angle).
  const double dotN = nPrev.x * nNext.x + nPrev.y * nNext.y;
  const double denom = 1.0 + dotN;
  if (denom <= 0) {
    // Opposite normals make the miter unbounded. The caller's miter-limit
    // check handles near-opposite normals before this point.
    return V + delta * nPrev;
  }
  return V +
         delta * vec2((nPrev.x + nNext.x) / denom, (nPrev.y + nNext.y) / denom);
}

double ValidMiterLimit(double miterLimit) {
  return std::isfinite(miterLimit) && miterLimit >= 2.0 ? miterLimit : 2.0;
}

// Offset a single input contour. Positive `delta` inflates the solid
// region (outer CCW rings expand outward; inner CW holes shrink as
// their boundary moves into the hole). Negative `delta` does the
// reverse. Orientation handling falls out automatically from
// OutwardNormal's right-of-edge-direction convention:
//   - CCW outer: right-of-edge = outward of solid.
//   - CW hole:   right-of-edge = into the hole = outward of solid.
// So `V + delta * OutwardNormal` always moves the boundary outward of
// the solid for delta > 0, regardless of ring orientation, and the
// convex/concave decision depends only on `cross * sign(delta)`.
SimplePolygon OffsetContour(const SimplePolygon& contour, double delta,
                            OffsetJoinType jt, double miterLimit,
                            double arcTol) {
  const int n = static_cast<int>(contour.size());
  if (n < 3 || delta == 0) return contour;
  const double deltaSign = (delta >= 0) ? 1.0 : -1.0;
  miterLimit = ValidMiterLimit(miterLimit);

  SimplePolygon out;
  out.reserve(static_cast<size_t>(n) * 2);
  for (int i = 0; i < n; ++i) {
    const vec2 V = contour[i];
    const vec2 P = contour[(i + n - 1) % n];
    const vec2 N = contour[(i + 1) % n];
    const vec2 ePrev = vec2(V.x - P.x, V.y - P.y);
    const vec2 eNext = vec2(N.x - V.x, N.y - V.y);
    const vec2 nPrev = OutwardNormal(ePrev);
    const vec2 nNext = OutwardNormal(eNext);
    if (nPrev == vec2(0, 0) || nNext == vec2(0, 0)) continue;
    const vec2 endPrev = V + delta * nPrev;
    const vec2 startNext = V + delta * nNext;
    // Solid's convex/concave at this vertex: cross > 0 is a convex
    // outer corner of the solid regardless of ring orientation (CCW
    // outer left-turn = convex outer; CW hole left-turn = hole's
    // concave indent = solid's convex bulge into the hole). Negative
    // delta flips the offset role (a solid-convex corner becomes a
    // shrinking-corner that needs a miter), hence the `* deltaSign`.
    const double cross = la::cross(ePrev, eNext);
    const double convex = cross * deltaSign;
    // Use the same scale-derived collinearity shape as CCW(): the tolerance is
    // a length from the larger adjacent edge, not an absolute cross-product.
    if (StraightTurn(ePrev, eNext)) {
      // Zero cross is either a straight continuation (dot >= 0: nPrev == nNext,
      // so endPrev == startNext - one point suffices) or an antiparallel
      // reversal (dot < 0: nPrev == -nNext, so the offset endpoints are
      // distinct and on opposite sides). Emit both for the reversal to bevel
      // across the spike instead of collapsing it to a single point.
      out.push_back(endPrev);
      if (dot(ePrev, eNext) < 0) out.push_back(startNext);
      continue;
    }
    if (convex < 0) {
      // Concave joins intentionally create a negative region that the final
      // Positive union removes, matching Clipper2's offset cleanup.
      out.push_back(endPrev);
      out.push_back(V);
      out.push_back(startNext);
      continue;
    }
    // Convex corner: apply join.
    out.push_back(endPrev);
    switch (jt) {
      case OffsetJoinType::Round:
        AppendRoundJoin(out, V, nPrev, nNext, delta, arcTol);
        break;
      case OffsetJoinType::Miter: {
        // miterLen / |delta| = 1 / cos(half_angle) =
        // sqrt(2 / (1 + dot(nPrev, nNext))). The limit
        // miterLen <= miterLimit * |delta| rearranges to
        // dot >= 2 / miterLimit^2 - 1. Comparing on the dot product
        // directly avoids materialising the miter point in the
        // clamped case (and avoids FP-fragility of MiterPoint when
        // the denominator (1 + dot) approaches zero - exactly the
        // case where we'd clamp anyway).
        const double dotN = nPrev.x * nNext.x + nPrev.y * nNext.y;
        const double miterCosThresh = 2.0 / (miterLimit * miterLimit) - 1.0;
        // Equality is allowed by the miter limit. `dotN` comes from rounded
        // unit normals, so use only the baseline unit-scale predicate epsilon
        // to avoid squaring a corner that is exactly on the limit.
        const double miterTieTol = EpsilonFromScale(1.0, /*k_budget=*/0);
        // Near-opposite normals make MiterPoint unbounded (it scales as
        // 1/(1 + dotN)), and the miterLimit gate stops bounding it once
        // 2/miterLimit^2 underflows miterTieTol (i.e. for very large
        // miterLimit). Independently square such degenerate-sharp corners so
        // the emitted coordinate stays bounded by sqrt(2 / kMinMiterDenom) *
        // |delta|; 2e-12 caps it at ~1e6 * |delta|, far beyond any meaningful
        // miter, so honest finite limits are unaffected.
        constexpr double kMinMiterDenom = 2e-12;
        if (dotN + miterTieTol < miterCosThresh ||
            1.0 + dotN < kMinMiterDenom) {
          AppendSquareJoin(out, V, nPrev, nNext, delta);
        } else {
          out.push_back(MiterPoint(V, nPrev, nNext, delta));
        }
        break;
      }
      case OffsetJoinType::Square:
        AppendSquareJoin(out, V, nPrev, nNext, delta);
        break;
      case OffsetJoinType::Bevel:
        break;
    }
    out.push_back(startNext);
  }
  return out;
}

Polygons RemoveCollinear(Polygons polys, double eps) {
  const double eps2 = eps * eps;
  for (auto& loop : polys) {
    if (loop.size() < 3) continue;
    SimplePolygon kept;
    kept.reserve(loop.size());
    const int n = static_cast<int>(loop.size());
    for (int i = 0; i < n; ++i) {
      const vec2 P = kept.empty() ? loop[(i + n - 1) % n] : kept.back();
      const vec2 V = loop[i];
      const vec2 N = loop[(i + 1) % n];
      const vec2 ePrev = vec2(V.x - P.x, V.y - P.y);
      const vec2 eNext = vec2(N.x - V.x, N.y - V.y);
      if (dot(ePrev, ePrev) < eps2) continue;  // zero-length back-edge
      if (dot(eNext, eNext) < eps2) continue;  // zero-length forward edge
      // Perpendicular squared distance from V to line PN; if less than
      // eps^2, V is essentially collinear with P-N and can be dropped.
      const vec2 pn = vec2(N.x - P.x, N.y - P.y);
      const double pnLen2 = dot(pn, pn);
      if (pnLen2 > 0) {
        const double cross = la::cross(ePrev, eNext);
        if (cross * cross < eps2 * pnLen2) continue;
      }
      kept.push_back(V);
    }
    // Wrap-around check: the first kept vertex may be collinear with
    // last-kept and kept[1] now that earlier collapses settled.
    while (kept.size() >= 3) {
      const vec2 P = kept.back();
      const vec2 V = kept.front();
      const vec2 N = kept[1];
      const vec2 ePrev = vec2(V.x - P.x, V.y - P.y);
      const vec2 eNext = vec2(N.x - V.x, N.y - V.y);
      const vec2 pn = vec2(N.x - P.x, N.y - P.y);
      const double pnLen2 = dot(pn, pn);
      const double cross = la::cross(ePrev, eNext);
      if (pnLen2 > 0 && cross * cross < eps2 * pnLen2) {
        kept.erase(kept.begin());
      } else {
        break;
      }
    }
    loop = std::move(kept);
  }
  // Drop sub-3-vertex degenerate rings that the collinear pass collapsed.
  polys.erase(
      std::remove_if(polys.begin(), polys.end(),
                     [](const SimplePolygon& l) { return l.size() < 3; }),
      polys.end());
  return polys;
}

}  // namespace

// Public Offset API. `delta` is the offset distance (positive = inflate,
// negative = inset). `miterLimit` is relative to |delta|. `arcTol` is
// the maximum perpendicular chord-error for Round joins; same semantics
// as Clipper2's `arc_tolerance`.
//
// Each input contour produces one offset ring, then all rings are regularized
// with Positive/Add filling. A final pass strips collinear vertices (matching
// Clipper2's `InflatePaths` finishing behaviour so callers see the same
// NumVert).
Polygons Offset(const Polygons& in, double delta, OffsetJoinType jt,
                double miterLimit, double arcTol, double tolerance) {
  if (delta == 0 || in.empty()) return in;
  // Reject NaN/Inf delta and input coordinates.
  if (!std::isfinite(delta)) return {};
  for (const auto& ring : in) {
    for (const auto& v : ring) {
      if (!std::isfinite(v.x) || !std::isfinite(v.y)) return {};
    }
  }
  Polygons offsetRings;
  offsetRings.reserve(in.size());
  for (const auto& ring : in) {
    auto off = OffsetContour(ring, delta, jt, miterLimit, arcTol);
    if (off.size() >= 3) offsetRings.push_back(std::move(off));
  }
  if (offsetRings.empty()) return {};
  const double eps = InferEps(offsetRings, {});
  // Resolve self-intersecting offset rings (e.g. when delta exceeds a
  // thin feature's half-width and the offset pinches itself into multiple
  // loops). CrossSection storage is normal-oriented before reaching Offset,
  // so Positive/Add cleanup keeps the filled side for both dilation and inset.
  Polygons unioned = Simplify(offsetRings, eps, tolerance);
  return RemoveCollinear(std::move(unioned), eps);
}

}  // namespace boolean2
}  // namespace manifold
