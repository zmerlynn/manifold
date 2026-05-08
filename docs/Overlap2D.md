# 2D overlap removal: design

Status: design proposal accompanying the prototype at
`extras/overlap2d_proto.cpp`.

References:

- [`RobustBoolean.pdf`](RobustBoolean.pdf): Julian Smith, *Towards
  robust inexact geometric computation*, UCAM-CL-TR-766, 2009. The
  α-budget framework, symbolic-perturbation predicates, and ≤2-iter
  convergence proof come from §7.
- [Emmett Lalish's algorithm sketch in issue #289](https://github.com/elalish/manifold/issues/289#issuecomment-2111069955).
  The 6-step BVH-on-Smith pipeline and the choice not to use a
  sweep-line.
- Edelsbrunner & Mücke, *Simulation of Simplicity*, ACM TOG 9(1), 1990.
  The permutation-parity SoS perturbation used by `CCWPerturbed`.
- Shewchuk, *Adaptive Precision Floating-Point Arithmetic and Fast
  Robust Geometric Predicates*, Discrete & Computational Geometry
  18(3), 1997. The standard reference for `orient2d`-style
  predicates; the prototype's `CCW`+SoS approach is in this lineage.
- de Berg, Cheong, van Kreveld & Overmars, *Computational Geometry:
  Algorithms and Applications* (3rd ed.), Springer 2008. Chapter 2
  introduces the DCEL specifically for boolean operations on
  subdivisions; chapter 8 covers the arrangement-and-classify
  template the prototype follows.
- Margalit & Knott, *An algorithm for computing the union, intersection
  or difference of two polygons*, Computers & Graphics 13(2), 1989.
  Pioneers the winding-rule classification (NonZero / EvenOdd) that
  Clipper2 directly inherits.

The shape is the canonical "build a planar arrangement, classify
faces by fill rule, reconstruct the boundary" template (de Berg ch.8),
implemented with the same DCEL data structure CGAL's `Arrangement_2`
uses.

## Goal

Explore whether a manifold-native 2D implementation could eventually
take over the **boolean and self-overlap-removal** parts of
`CrossSection`'s current Clipper2 use, specifically
`Boolean(Add/Subtract/Intersect)` and `Simplify`.

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

**Algorithmic spine: Emmett's #289 sketch.** BVH pair queries
instead of sweep-line; a 6-step pipeline (vert merge → edge collapse
→ on-edge verts → edge-edge intersection → sub-edge canonicalization
→ winding filter); reuse manifold's 3D building blocks
(`Collider`, `DisjointSets`, the orientation kernels in
`src/shared.h`) so the 2D code lives in the same arithmetic and
spatial-query world as `Manifold::Impl::Boolean3`.

**Correctness framework: Smith (UCAM-CL-TR-766, §7).** Provides
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
  We build a doubly-connected edge list (the same structure
  `Manifold::Impl::halfedge_` uses for 3D), assign one winding per
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

**Pass 1 is topology-correct on every fuzz case.** A 336,000-case
fuzz (12,000 seeds × 4 sizes × 7 displacement scales) shows pass 1
produces topology-valid output on every input, including adversarial
displacements up to 2⁴⁹ where ULP spacing exceeds feature size by
orders of magnitude. Topology is the hard contract; pass 1 satisfies
it across the entire battery.

What pass 1 is *not* is idempotent in all cases. Feeding pass-1
output back as pass-2 input occasionally produces a different
result, which is what the prototype's `IterateToFixedPoint` path
(default `maxIter=2`, Smith's bound) studies. ~2.4% of fuzz cases
differ between pass 1 and pass 2 at sub-ε position quanta but
match at the ε quantum (topology-identical, geometry shifted by a
few ULPs). One in 336,000 cases hits a more pathological mode:
pass 1 produces a topologically-valid but geometrically *thin*
output polygon (all interior dimensions within 2ε), and pass 2
correctly notices that the thin polygon's "long way around"
coincides to within ε with its "short way around" and collapses
it to empty. Consistent behavior under the algorithm's own rules,
but means pass-1 output is not always a stable fixed point.

**The right `maxIter` default for production is operation-dependent,
not a single number.** Manifold's 3D analog of `Simplify` was
explicitly switched from one-pass to iterate-to-fixed-point
([#675](https://github.com/elalish/manifold/issues/675), closed by
[#894](https://github.com/elalish/manifold/pull/894)) because users
filed bugs ([#1262](https://github.com/elalish/manifold/issues/1262),
[#1239](https://github.com/elalish/manifold/issues/1239)) about
topology-valid-but-user-hostile outputs. The thin-polygon collapse
above is exactly that bug class. So the defensible split:

- `Boolean(Add/Subtract/Intersect)`: `maxIter=1` is reasonable. Inputs
  are presumed clean; iteration mostly refines positions.
- `Simplify`: iterate-to-fixed-point. Inputs are dirty by definition;
  users are explicitly asking for cleanup, and the 3D-side precedent
  is direct.

The prototype keeps `maxIter=2` as its default for tail-behavior
study; the production landing should pick per-op as above. The
landing should also include an **idempotency regression test** in
addition to topology validity, since the thin-polygon collapse and
its kin slip past topology-only checks.

## API

```cpp
namespace overlap2d {

// Standard CSG fill rules: SVG/PostScript/Clipper2/CGAL convention.
enum class FillRule { Positive, Negative, EvenOdd, NonZero };

// Single-input regularization. Same name as `CrossSection::Simplify`
// (the eventual landing target).
manifold::Polygons Simplify(const manifold::Polygons& in, double eps,
                            FillRule fill = FillRule::Positive);

// Binary boolean. Uses manifold's existing OpType enum from
// `include/manifold/common.h` (the same one Manifold::Boolean and
// CrossSection::Boolean accept), so callers can pass the same value
// through 3D and 2D codepaths.
manifold::Polygons Boolean2D(const manifold::Polygons& a,
                             const manifold::Polygons& b,
                             manifold::OpType op,
                             double eps = 0.0);  // <=0 auto-infers

// Symmetric difference. Standard 4th binary boolean alongside
// Add/Subtract/Intersect. `manifold::OpType` doesn't include Xor, so
// it's a separate entry point.
manifold::Polygons Xor(const manifold::Polygons& a,
                       const manifold::Polygons& b, double eps = 0.0);
}
```

Step 6's classification is internally a `WindPredicate` functor that
maps a face's winding number to inside/outside. The public `FillRule`
enum is the user-facing wrapper. Subtract flips B's edge
multiplicities to −1 before joining, so the same `Positive` rule
carves A − B. Intersect uses an internal "w > 1" predicate (not a
public fill rule, since it depends on Boolean2D maintaining mult=+1
per input). Xor combines inputs with mult=+1 each and applies
`EvenOdd`.

## Validation

| Test | Pass rate | Notes |
|------|-----------|-------|
| Standard battery | 350/350 | 13 named cases + 350 random topological polygons. Topology-balance: every vertex's input edge contribution = output. |
| Displacement fuzz | 450/450 | 5 scales (2<sup>10</sup>...2<sup>49</sup>) × 3 sizes × 30 seeds. |
| DeepFuzz (12,000 seeds) | **336,000 / 336,000 first-pass topology valid** | 12,000 seeds × 4 sizes × 7 displacement scales (2<sup>10</sup> to 2<sup>49</sup>). 0 first-pass FAILs across the entire battery. Iter distribution after fixed-point: 1:327,903; 2:7,937; 3:152; 4:7; 5:1. 0 iteration-DEGRADED (pass-1 valid → broken by iteration). 1 geometric collapse (thin-polygon idempotence quirk; see Pipeline section). |
| Polygons round-trip | 335,440 / 336,000 | 0.17% mismatch rate (`OutEdgesToPolygons` drops <3-vert loops the lower-level pipeline counts as edges; documented gap). Mismatch rate stable across the 100-seed and 12,000-seed runs. |
| Boolean2D area regression | 3/3 | Two CCW unit squares offset diagonally; Add area=7, Subtract=3, Intersect=1 (exact). Drives the wind > 1 predicate against perpendicular axis-aligned crossings. |
| Axis-aligned perpendicular | 1/1 | Two CCW unit squares offset diagonally; produces L-shape union with 8 boundary verts. Regression test for the kernel bug surfaced during OQ #3 / #4 investigation. |
| `polygon_corpus.txt` (100 entries) | 100/100 topology valid | Manifold's existing curated corpus, run via `./overlap2d_proto corpus`. Area-preservation breakdown: 87/100 < 1% drift, 1/100 mid-drift (`Precision4` at 2.6%), 7/100 collapsed by name (`Degenerate4`, `DegenerateRing`, `NearlyLinear`, `Looping1/2`, `Precision`, `Sweep` - cases the algorithm correctly identifies as eps-thin). 5/100 had zero-area input. The `Precision4` 2.6% drift is correct behavior: the input is a sliver polygon with an eps-thin tail, and the algorithm trims sub-eps features as eps says it should. |
| DeepFuzz area drift (100-seed) | 8 cases > 1% (0.29%) | Quantitative oracle on top of the topology check. All 8 cases at kPow=35 (where ULP precision is weakest); max drift 43.9%. Topology balance alone showed 1 geometric-collapse case at the same fuzz size; area oracle revealed 8x more cases of meaningful iteration drift, all in the thin-polygon class. |

Displacement fuzz is the load-bearing test: traditional FP boolean
libraries fail catastrophically at 2<sup>k</sup> displacement, which
is exactly what Smith's framework exists to handle. The 336,000-case
run verifies that pass 1 produces topology-valid output across the
entire adversarial range, in every size and seed tested.

## Open questions

1. **Landing site in `src/`.** New `src/overlap2d/`, or wired into
   `cross_section.cpp` directly as the boolean/`Simplify` backend.
   Driven by how much shared code with 3D booleans we expose.
2. **eps policy.** `Boolean2D` accepts explicit eps or auto-infers
   when ≤0. `CrossSection`'s existing API doesn't expose eps; add it
   as an optional parameter, or always auto-infer? Smith's α-budget
   is conservative (k<sub>budget</sub>=1000); `Simplify` may want
   stricter.
3. **Multi-threading.** `manifold::Collider` is parallel-ready; the
   prototype is `MANIFOLD_PAR=-1` only. Step 1 union-find and step 6
   per-face ray-cast are the natural parallelism candidates.
4. **Additional quantitative oracles.** The prototype now has area
   preservation and area-drift tracking on top of the topology check
   (see Validation), and consumes `polygon_corpus.txt` directly. The
   remaining easy add is Clipper2's `Polygons.txt` (195 numbered
   cases with `SOL_AREA` references), which gives an independent
   third-party oracle on a wider input distribution than
   `polygon_corpus.txt`. Worth doing for production landing; not
   strictly required for the design discussion.

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
   what comes out of step 2.

## Alternatives considered

We tried these and rejected each for the reason listed; recording so
we don't relitigate.

- **Embed at z=0 and reuse 3D booleans wholesale.** Exposes too much
  3D-specific machinery (manifoldness assumptions, triangle vs.
  segment); cleaner to share the kernels (`Shadows`, `Interpolate`,
  `Collider`, `DisjointSets`) than the whole boolean3 pipeline.
- **`Kernel11` from `boolean3.cpp` directly for step 4.** Requires
  "one endpoint inside, one outside" the other segment's projection,
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
