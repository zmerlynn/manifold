# Boolean output inherits polluted `tolerance_` from FLT_MAX-scale operands

> **Filed upstream as elalish/manifold#1714, converted to discussion
> #1715 and closed as INTENTIONAL (May 2026).**
>
> Maintainer (Emmett Lalish) response:
> > This is intentional behavior. The idea is to bound floating-point
> > error conservatively — if you actually make a new vertex by
> > intersecting a really huge face with a really tiny edge, in
> > general its meaningful precision will be the lower of the two,
> > even if all of those intersections happen to cluster into a
> > small bounding box. There are potentially fancier things we
> > could do here, but I'd need a pretty clear user need to go down
> > that road. I think in general it's better to warn people that
> > their tolerance is getting wide and hopefully they can transform
> > their inputs appropriately.
>
> Framing: `tolerance_` is "worst-case FP error from the arithmetic
> that produced this vert," NOT "geometric precision at the vert's
> location." Intersecting FLT_MAX-scale geometry with normal-scale
> geometry produces verts whose worst-case FP error IS bounded by
> the FLT_MAX scale, regardless of where they happen to land.
> Conservatively reporting that as `tolerance_` is the design
> intent. Cray is degenerate at the FP level and out of scope for
> clean handling.
>
> Implications: the defensive clamp in `RunOverlapRemovalImpl`
> (commit `4678ad72`) is a **permanent workaround** for our
> pipeline, not a transitional port of a future upstream fix. We're
> claiming more precision than upstream's contract supports because
> our per-vertex predicates need a meaningful tolerance to function,
> and we know cray-shaped inputs are inherently pathological for
> RemoveSelfIntersections regardless of what the FP-error bound
> would technically be. Keep the clamp.
>
> Companion: the sliver bug (`docs/boolean-sliver-bug.md` / #1706)
> is also a Boolean-engine output property, but is being fixed
> upstream via a planned QEM decimator. Both bugs are independent.

## Summary

When `Manifold::Boolean` produces an output whose operands had wildly
different scales, the result's `tolerance_` field is inherited from the
larger-scale operand rather than re-inferred from the actual output
bbox. Downstream consumers that consult `impl.tolerance_` then treat
all positions in the (correctly small-scale) output as coincident.

This is a 2-year-old latent bug that has likely never been hit hard
because most downstream consumers don't depend on `tolerance_` for
high-impact decisions. Our work on `Manifold::RemoveSelfIntersections`
exposed it because the pipeline runs many per-vertex geometric
predicates that the polluted value can break.

## Reproduction

Test models are checked in at `test/models/Cray_{left,right}.obj`
(originally added in PR #714, see History below).

```cpp
#include "manifold/manifold.h"
#include <fstream>
using namespace manifold;

std::ifstream cl("test/models/Cray_left.obj");
std::ifstream cr("test/models/Cray_right.obj");
Manifold m1 = Manifold::ReadOBJ(cl);
Manifold m2 = Manifold::ReadOBJ(cr);
Manifold cray = m1 - m2;
// cray.GetTolerance() == 3.40282e+26
// cray.BoundingBox() spans only [-24534, 25400]^3
```

Observed values:

| | bbox.x range | volume | epsilon | tolerance |
|---|---|---|---|---|
| Cray_left | [-25400, 25400] | 1.96e+12 | 2.54e-08 | 2.54e-08 |
| Cray_right | [-3.40e+38, 3.40e+38] | 1.58e+116 | 3.40e+26 | 3.40e+26 |
| Cray_left - Cray_right | [-24534, 25400] | 4.91e+10 | **3.40e+26** | **3.40e+26** |

The result has a sensible ~25400-unit bbox and reasonable volume but
carries forward Cray_right's `tolerance_` of `3.40e+26`. Anything
within `3.4e+26` of anything else is "the same point" by tolerance,
which is far larger than the entire output mesh.

`Cray_right.obj` is 8 vertices and 12 triangles (a single cube), with
extents at `+/- FLT_MAX` in x and y and `[-FLT_MAX, -100]` in z. It is
used as a near-infinite chopping half-space to clip Cray_left to
`z > -100`.

## Why it survived

1. **PR #714 (Jan 2024) fixed the immediate symptom.** Issue #712 was
   a `Status()` crash from NaN propagation when the Boolean engine
   processed FLT_MAX-scale geometry. The PR added `isfinite(...)`
   guards in `boolean3.cpp`'s `Interpolate` and `Intersect` helpers
   so `0 * inf` couldn't propagate into output vertex positions. The
   crash was fixed; nobody asked whether the resulting `tolerance_`
   was meaningful.

2. **Most downstream consumers don't depend on `tolerance_`.** The
   field is informational for many code paths. Visualization,
   serialization, and most queries don't use it as a decision
   threshold.

3. **`RemoveSelfIntersections` is the first heavy consumer.** Our
   pipeline runs per-vertex epsilon-offset winding probes
   (`AnalyzeSelfMesh`) and many edge-face / edge-edge predicates.
   The polluted tolerance broke the per-vertex classifier signal
   for cray (= histogram shows no clean `(0,1)` boundary verts,
   only `(-4,-4)`, `(0,0)`, `(1,3)`).

4. **No follow-up PR has touched the area.** `git log --grep
   tolerance` between 2024-01 and now shows nothing relevant.

## Affected code paths

The `tolerance_` field is set in two places:

- `src/impl.h:308` (`Impl::Impl(MeshGLP)`): `tolerance_ = meshGL.tolerance`
  copies whatever the input MeshGL64 carries. This is correct for
  user-supplied meshes but wrong for Boolean outputs that get
  serialized through MeshGL.
- `src/impl.cpp:617` (`Impl::SetEpsilon`): `tolerance_ = std::max(
  tolerance_, minTol)`. The `std::max` is the load-bearing line, it
  prevents `SetEpsilon` from ever SHRINKING `tolerance_`. Designed to
  preserve user-set lower bounds, but it makes recovery from a
  polluted value impossible without manually zeroing first.

The Boolean engine's output construction lives in
`src/boolean_result.cpp` and `src/csg_tree.cpp` (CsgOpNode finalize).
Look for where the result `Impl` gets its `tolerance_` assigned.

## Why the propagation rule is intentional (per maintainer)

`Manifold::Boolean`'s output `tolerance_` semantic is "worst-case
FP error bound across the arithmetic that produced this output."
For an output vertex computed by intersecting a face from a
FLT_MAX-scale operand with an edge from a normal-scale operand,
the FP arithmetic involved scales with the larger operand —
catastrophic cancellation when subtracting near-equal FLT_MAX
values can lose 20+ decimal digits regardless of where the result
lands geometrically. Reporting `tolerance_ = max(operand
tolerances)` conservatively captures this.

A fix that re-infers from the output bbox would claim more
precision than the arithmetic actually has, even when the
geometric output looks well-conditioned. That trade-off doesn't
match the contract — `tolerance_` is the consumer-visible
guarantee that "any two points within this distance might be the
same point at FP precision."

The maintainer's suggested user-facing path: when your input
tolerances are getting wide, transform inputs to a normalized
scale before the Boolean.

## Permanent workaround in `RemoveSelfIntersections`

Commit `4678ad72` (`overlap_removal: defensive tolerance reset for
upstream-polluted inputs`) adds a clamp at the entry of
`RunOverlapRemovalImpl`:

```cpp
if (impl.bBox_.IsFinite() &&
    impl.tolerance_ > impl.bBox_.Scale()) {
  impl.tolerance_ = 0;
  impl.epsilon_ = 0;
  impl.SetEpsilon(-1, false);
}
```

This catches cray (= the only fixture with this pathology) and
re-infers a sensible value. It is intentionally narrow (only fires
when tolerance exceeds bbox scale, a clear sign that the input came
from a much-larger-scale Boolean operand). The `tolerance_ = 0;
epsilon_ = 0;` reset is needed because `SetEpsilon`'s `std::max`
would otherwise preserve the inherited value.

The clamp is technically lying to the rest of the pipeline about
the FP-error bound (= upstream's `tolerance_` is conservatively
correct, ours is "what makes our per-vert predicates work"). For
inputs where the pathology applies (= cray), the FP-error bound
already says "any answer is fine" — so our claim is no worse than
the conservative one for our purposes.

## Diagnostic helpers added in this branch

`src/overlap_removal.cpp` has env-gated diagnostics that are useful
for verifying the fix:

- `OVERLAP3D_DEBUG_CLASSIFIER=1` prints the per-vertex `(w_above,
  w_below)` histogram. For a healthy input most verts should be
  `(0,1)` (clean boundary). cray currently shows `(-4,-4):9
  (0,0):68 (1,3):46` (= no boundary signal). After the upstream fix
  the histogram should change and may show `(0,1)` verts.
- `OVERLAP3D_DEBUG_FALLBACK=1` prints exception messages and
  post-pipeline non-manifold status reasons.

Run via:
```
OVERLAP3D_DEBUG_CLASSIFIER=1 build/extras/overlap3d_proto cray --prod-api
```

## Verification (current state)

Our defensive clamp is the end state, not a workaround pending
upstream fix. Verification is that:

1. `Manifold::ReadOBJ("Cray_right.obj").GetTolerance()` is `~3.4e+26`
   (unchanged, upstream-correct).
2. `(Manifold::ReadOBJ(cl) - Manifold::ReadOBJ(cr)).GetTolerance()`
   is `~3.4e+26` (unchanged, upstream-intentional).
3. Inside our pipeline, the clamp fires for cray and we use a
   bbox-derived tolerance instead. Verifiable via
   `OVERLAP3D_DEBUG_CLASSIFIER` — the `(wa, wb)` histogram for cray
   still shows no clean `(0,1)` verts because the underlying input
   is geometrically degenerate (= the clamp lets the predicate run
   but doesn't manufacture signal that isn't there).
4. Cray remains a permanent passthrough fixture (= the gate falls
   back to input). This is the expected outcome.

## History

- **PR #714 (Jan 2024)**: "Fixed float overflow." Closed issue #712.
  Author: Emmett Lalish. Diff was in `boolean3.cpp`'s `Interpolate`
  and `Intersect`. Did not touch `tolerance_` handling.
- **PR #1529 (Feb 2026)**: "make MeshIO private." Touched the test
  models (= the format conversion that put `Cray_*.obj` files in
  the current location). Did not address tolerance.
- **commit 4678ad72 (May 2026, this branch)**: Defensive clamp in
  `RemoveSelfIntersections`. Treats the symptom; awaits root-cause
  fix.
- **Issue #1714 → Discussion #1715 (May 2026)**: filed against
  elalish/manifold with the explicit-vs-inferred framing. Maintainer
  closed as INTENTIONAL — propagating max-of-operand tolerances is
  the conservative FP-error bound the contract promises, and a fix
  that re-infers from output bbox would claim more precision than
  the arithmetic supports. Our defensive clamp stays as a
  permanent workaround for the pipeline.

## Out of scope

- The "warn people their tolerance is getting wide and they can
  transform inputs appropriately" path the maintainer suggested is
  a user-facing concern, not a library code change. If we wanted
  to surface a warning to RemoveSelfIntersections callers when
  input tolerance is suspiciously wide, that would be additive on
  top of the current defensive clamp.
