# 2D overlap removal — design

Status: design proposal accompanying the prototype at
`extras/overlap2d_proto.cpp`.
References: [`RobustBoolean.pdf`](RobustBoolean.pdf) (Julian Smith,
*Towards robust inexact geometric computation*, UCAM-CL-TR-766, 2009);
Emmett Lalish's algorithm sketch in
[issue #289](https://github.com/elalish/manifold/issues/289).

## Goal

Explore whether a manifold-native 2D implementation could eventually
take over the **boolean and self-overlap-removal** parts of
`CrossSection`'s current Clipper2 use — i.e. `Boolean(Add/Subtract/
Intersect)` and `Simplify`. Clipper2's polygon-offset path
(`ClipperOffset`, used by `CrossSection::Offset`) is a Minkowski sum
with a disk under various join styles; it's a different problem and
stays on Clipper2 (or becomes its own project later). Even if this
prototype were fully productionized, `CrossSection` would still link
Clipper2 for offset.

The prototype validates the approach end-to-end; this doc proposes
landing it in `extras/` as reference and gathering feedback before
deciding on production scope.

Design constraints the prototype satisfies:

1. Handles self-intersecting input polygons (Smith treats self-overlap
   and inter-polygon overlap uniformly).
2. Is robust under FP arithmetic, including the displacement attack
   (geometry near 2<sup>k</sup> coordinate boundaries where ULP
   spacing is comparable to feature size).
3. Reuses 3D booleans' building blocks: `manifold::Collider` for
   spatial queries, `manifold::DisjointSets` for vertex equality, and
   the orientation kernels in `src/shared.h`.
4. Supports Add / Subtract / Intersect plus the unary "remove
   overlaps" simplification used by `CrossSection::Simplify`.

## Approach

**Algorithmic spine — Emmett's #289 sketch.** BVH pair queries
instead of sweep-line; a 6-step pipeline (vert merge → edge collapse
→ on-edge verts → edge-edge intersection → sub-edge canonicalization
→ winding filter); reuse manifold's 3D building blocks
(`Collider`, `DisjointSets`, the orientation kernels in
`src/shared.h`) so the 2D code lives in the same arithmetic and
spatial-query world as `Manifold::Impl::Boolean3`.

**Correctness framework — Smith (UCAM-CL-TR-766, §7).** Provides
the part Emmett's sketch doesn't: a concrete ε formula and a proof
of convergence in ≤2 iterations under bounded FP error
(§7.7, fig 7.16):

> ε = (k<sub>budget</sub> + 1) · √153 · u · L

with u = 2<sup>−53</sup>, L = bounding-box half-extent rounded up to
a power of 2, k<sub>budget</sub> = 1000. The operation has a single
tunable parameter with clear semantics ("how much movement is
acceptable per edge"), not a per-call epsilon left to the caller.
Cross-detection uses `CCW`+SoS (Edelsbrunner–Mücke
permutation-parity perturbation); position uses
trim-and-`Interpolate` for FP stability.

**Local divergences from both.** Three places where the prototype
goes its own way, each prompted by a problem that surfaced during
development:

- **Step 6: DCEL face traversal instead of Smith's per-edge
  ray-cast.** Per-edge ray-casts produce inconsistent verdicts at
  shared vertices (270 first-pass topology FAILs in DeepFuzz 14000).
  We build a doubly-connected edge list — the same structure
  `Manifold::Impl::halfedge_` uses for 3D — assign one winding per
  planar face via a single ray-cast from inside, and keep an edge
  iff its left and right faces disagree on the inside/outside
  predicate. Structurally consistent at shared vertices by
  construction.
- **Step 4b: structural re-merge of intersection verts.** When 3+
  edges meet at one true point, step 4 inserts an intersection vert
  per pair (so up to 3 distinct verts where there should be 1).
  Union-find with a structural gate (verts share an input edge) plus
  a geometric gate (within 10ε) collapses them. Not in either source
  algorithm; it's a fix for a problem step 4 introduces.
- **Centered shoelace area.** At 2<sup>49</sup> displacement, raw
  shoelace `a.x·b.y − b.x·a.y` produces ULP ~2000 on a unit-area
  face, randomizing the sign. Centering each face's coordinates on
  a boundary vert before summing keeps products on order
  (edge-length)<sup>2</sup>; outer-face detection becomes stable.

## Pipeline

| Step | Job | Implementation notes |
|------|-----|---------------------|
| 1 | Vertex merge | BVH broad phase over ε-padded vert AABBs; `DisjointSets` union-find within ε. |
| 2 | Edge collapse | Drop edges whose endpoints map to the same merged vertex. |
| 3 | Edge → on-edge vert list | BVH broad phase over ε-padded edge AABBs queried by vert AABBs; exact projection-distance gate. Skips "thin-triangle apex" verts (apex of a degenerate triangle whose base is the candidate edge). |
| 4 | Edge-edge intersections | BVH broad phase over edge AABBs; trim-and-`Interpolate` symbolic kernel atop `CCW`/SoS. New verts snap to nearby existing verts before insertion. |
| 4b | Re-merge intersection verts | When 3 edges meet at one true point, step 4 inserts up to 3 distinct intersection verts; union-find with structural (verts share an input edge) + geometric (within 10ε) gates collapses them. |
| 5 | Sub-edge canonicalization | `map<lex(v0,v1), signed_mult>`; collinear sub-edges with opposing multiplicity cancel. |
| 6 | Winding filter | DCEL face traversal; one winding per face via ray-cast from a perpendicular-LEFT offset of any boundary half-edge. Edge retained iff `isInside(leftFace) != isInside(rightFace)`. |

**Pass 1 is the production path.** A 336,000-case fuzz (12,000
seeds × 4 sizes × 7 displacement scales) shows pass 1 produces
topology-valid output on every case, including adversarial
displacements up to 2⁴⁹ where ULP spacing exceeds feature size by
orders of magnitude. The recommended production default is
`maxIter=1`: topology is the contract that consumers care about,
and pass 1 satisfies it on every input we've thrown at it.

What pass 1 is *not* is idempotent in all cases. Feeding pass-1
output back as pass-2 input occasionally produces a different
result, which is what the prototype's `IterateToFixedPoint` path
(default `maxIter=2`, Smith's bound) studies. ~2.4% of fuzz cases
differ between pass 1 and pass 2 at sub-ε position quanta but
match at the ε quantum (i.e., topology-identical, geometry shifted
by a few ULPs). One in 336,000 cases hits a more pathological
mode: pass 1 produces a geometrically valid but *thin* output
polygon (all interior dimensions within 2ε), and pass 2 correctly
notices that the thin polygon's "long way around" coincides
to within ε with its "short way around" and collapses it to
empty. This is consistent behavior, not a bug, but it means
pass-1 output is not always a stable fixed point of the
algorithm. Production users on `maxIter=1` get the topology-
valid pass-1 output and don't see the collapse; iterate-to-fixed-
point users get the collapse as a "this region is sub-ε in some
dimension" signal.

The prototype keeps `maxIter=2` as its default to match Smith's
framework citation and to make iter-≥-2 cases easy to study; only
the production landing should switch to `maxIter=1`.

## API

```cpp
namespace overlap2d {

enum class BoolOp { Add, Subtract, Intersect };

manifold::Polygons OverlapRemovePolygons(const manifold::Polygons& in,
                                         double eps);

manifold::Polygons Boolean2D(const manifold::Polygons& a,
                             const manifold::Polygons& b,
                             BoolOp op, double eps = 0.0);  // <=0 auto-infers
}
```

Step 6's winding test is parameterized via a `WindPredicate` functor
(default `w > 0` for Add/Subtract, `w > 1` for Intersect). Subtract
flips B's edge multiplicities to −1 before joining; the same `w > 0`
predicate then carves A − B.

## Validation

| Test | Pass rate | Notes |
|------|-----------|-------|
| Standard battery | 350/350 | 13 named cases + 350 random topological polygons. Topology-balance: every vertex's input edge contribution = output. |
| Displacement fuzz | 450/450 | 5 scales (2<sup>10</sup>...2<sup>49</sup>) × 3 sizes × 30 seeds. |
| DeepFuzz (12,000 seeds) | **336,000 / 336,000 first-pass topology valid** | 12,000 seeds × 4 sizes × 7 displacement scales (2<sup>10</sup> to 2<sup>49</sup>). 0 first-pass FAILs across the entire battery. Iter distribution after fixed-point: 1:327,903; 2:7,937; 3:152; 4:7; 5:1. 0 iteration-DEGRADED (pass-1 valid → broken by iteration). 1 geometric collapse (thin-polygon idempotence quirk; see Pipeline section). |
| Polygons round-trip | 335,440 / 336,000 | 0.17% mismatch rate (`OutEdgesToPolygons` drops <3-vert loops the lower-level pipeline counts as edges; documented gap). Mismatch rate stable across the 100-seed and 12,000-seed runs. |
| Boolean2D smoke | 3/3 | Two overlapping unit squares → Add 8 edges, Subtract 4, Intersect 4. |
| Axis-aligned perpendicular | 1/1 | Two CCW unit squares offset diagonally; produces L-shape union with 8 boundary verts. Regression test for the kernel bug surfaced during OQ #3 / #4 investigation. |

Displacement fuzz is the load-bearing test: traditional FP boolean
libraries fail catastrophically at 2<sup>k</sup> displacement, which
is exactly what Smith's framework exists to handle. The 336,000-case
run verifies that pass 1 produces topology-valid output across the
entire adversarial range, in every size and seed tested.

## Open questions

1. **Landing site in `src/`.** New `src/overlap2d/`, or wired into
   `cross_section.cpp` directly as the boolean/`Simplify` backend
   (Clipper2 stays linked for `Offset` either way). Driven by how
   much shared code with 3D booleans we expose.
2. **eps policy.** `Boolean2D` accepts explicit eps or auto-infers
   when ≤0. `CrossSection`'s existing API doesn't expose eps; add it
   as an optional parameter, or always auto-infer? Smith's α-budget
   is conservative (k<sub>budget</sub>=1000); `Simplify` may want
   stricter.
3. **Multi-threading.** `manifold::Collider` is parallel-ready; the
   prototype is `MANIFOLD_PAR=-1` only. Step 1 union-find and step 6
   per-face ray-cast are the natural parallelism candidates.

## Deferred (graduation patch)

Mechanical/additive, batched with the build-graph wiring:

- `std::vector` → `manifold::Vec` rename.
- `ZoneScoped` Tracy markers at phase boundaries.
- `ExecutionContext::Impl*` ctx threading.
- Internal namespace (`overlap2d::detail`) hiding everything but the
  public entry points.

## Spin-out plan

This proposal covers the first two items only; the rest depends on
feedback from this discussion.

1. **Kernel extraction** (commit `388e3612`, internal): lift
   `withSign`, `Interpolate`, `Shadows` from `boolean3.cpp` into
   `src/shared.h`. Independent prerequisite for the prototype's
   step-4 kernel and SoS predicates.
2. **Drop `extras/overlap2d_proto.cpp` in-tree as reference**, behind
   an opt-in CMake option so CI tracks it but the default build
   doesn't link it. Natural place for design feedback before
   committing to a wider scope.
3. *(TBD pending discussion.)* Productionize into `src/` with the
   deferred items above, wired into `CrossSection` for the boolean
   and `Simplify` paths. Whether it lives next to or inside
   `CrossSection`, and whether it sits alongside Clipper2-for-
   booleans during a transition or replaces it outright, depends on
   what comes out of step 2. `CrossSection::Offset` continues to
   use Clipper2's `ClipperOffset` regardless.

## Alternatives considered

We tried these and rejected each for the reason listed; recording so
we don't relitigate.

- **Embed at z=0 and reuse 3D booleans wholesale.** Exposes too much
  3D-specific machinery (manifoldness assumptions, triangle vs.
  segment); cleaner to share the kernels (`Shadows`, `Interpolate`,
  `Collider`, `DisjointSets`) than the whole boolean3 pipeline.
- **`Kernel11` from `boolean3.cpp` directly for step 4.** Requires
  "one endpoint inside, one outside" the other segment's projection —
  a sweep-line invariant. BVH pair queries don't guarantee this; most
  pairs are nested-axis. **Replaced** with trim-and-`Interpolate`:
  trim both segments to their axis-projection overlap first, then
  apply Intersect's closed form. Same FP-stability bias (smaller
  |dy| endpoint), works on any pair.
- **Cramer's rule for step-4 position.** Tested vs. trim-and-
  `Interpolate`; identical iter distribution (1:431, 2:17, 3:2 on an
  earlier fuzz run). Kept `CCW` for cross-detection
  (it's the orientation predicate we already have) but moved
  position computation to the symbolic-perturbation form because
  Smith's α-budget proof requires symbolic positions.
- **BFS from the outer face** for assigning per-face winding (once
  we'd switched to DCEL). Doesn't propagate between disconnected
  components of the face graph (real for self-intersecting inputs
  whose union has multiple disjoint regions). **Replaced** with
  per-face ray-cast (`O(F · E)`, independent per face).
- **Step 4b at threshold ε.** Causes ITER FAILs at shallow crossings.
  Tried 1.5ε (fails), 3ε (still fails on adversarial cases), **10ε**
  (sweet spot, covers crossings down to ~6° angle without
  over-merging unrelated intersections), 100ε (over-merges). The
  structural gate (shared input edge) prevents over-merging; the
  10ε figure is the geometric width *given* the structural gate.
- **Step 3 "skip if vert is connected to one or both edge
  endpoints".** Initial attempt regressed n=200 fuzz. **Final:** skip
  only when connected to **both**, which is the thin-triangle apex
  case Smith discusses. One-endpoint adjacency is normal
  polygon-neighbor configuration.
- **`PerturbScheme::InitFromEdges` symbolic-perturbation tracker** as
  a separate phase. Not needed once the trim-and-`Interpolate` kernel
  resolves all-collinear pairs via permutation-parity SoS inline.
  **Retired** before the prototype's current shape.
- **Bumping `maxIter` past 2 to chase the iter ≥ 3 tail.** Smith's
  bound is 2; tested up to 8. The 12,000-seed DeepFuzz finds 152
  cases at iter=3, 7 at iter=4, 1 at iter=5 (out of 336,000), all
  ultimately converging to topology-valid output. Bumping the cap
  catches more fingerprint-stable convergence but doesn't change
  topology, since pass 1 was already topology-correct on every
  case. Production default of `maxIter=1` skips this entirely
  (see Pipeline section); the prototype's `maxIter=2` matches
  Smith's bound and is the right knob for studying the tail.
