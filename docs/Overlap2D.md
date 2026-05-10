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

*Why `Collider` (3D BVH at z=0) and not `tree2d` (manifold's existing
2D kd-tree used by the polygon triangulator)?* Tree2d's API is
points-only (`BuildTwoDTree(VecView<PolyVert>)` +
`QueryTwoDTree(rect, F)`); 3 of the 5 spatial-query sites in this
pipeline operate on points (steps 1, 4b, plus the eager-propagation
sub-phase of step 4) and would migrate cleanly. The other 2 sites
(step 4 main edge-edge broad phase, step 3 edge AABB queries) need
an edge-AABB tree, which only `Collider` provides. Mixing both
tree2d (for point sites) and `Collider` (for edge sites) would mean
two spatial data structures and two query patterns in the same
algorithm, and tree2d's serial kd-tree would lose the parallel
broad-phase Emmett's #289 sketch explicitly cited as a reason for
choosing BVH over sweep-line. Sticking to `Collider` everywhere is
the lower-complexity, parallel-ready choice; if the production
landing wants a different split (e.g., tree2d for point sites with
a custom edge-AABB tree elsewhere), that's reasonable but distinct
from the prototype's "use the existing parallel-ready BVH" decision.

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
- **Step 4 has an eager-propagation phase + step 4b is a residual
  cleanup.** When k≥3 input edges are concurrent at one true point,
  pair-by-pair processing (one BVH-pair-query at a time) produces
  per-pair intersection verts that should geometrically be one. Step
  4 addresses this in two phases: (a) the standard pair-by-pair
  insertion with snap-to-existing on the current pair's edges, and
  (b) eager propagation of each freshly-allocated intersection vert
  to *all* other edges within ε via a BVH range query. Step 4b is a
  union-find pass on the residue: near-duplicate intersection verts
  whose FP-rounded positions disagree by more than ε (typical at
  displaced coords with shallow crossings) get merged via a
  structural gate (verts share an input edge) plus a geometric gate
  (within 10ε). The two-phase step 4 + cleanup-pass step 4b mirror
  manifold's 3D pipeline exactly: step 4's propagation is the 2D
  analog of `boolean_result.cpp::AddNewEdgeVerts` (which adds an
  edge×face intersection vert to three lists, not just the source
  pair); step 4b is the 2D analog of `edge_op.cpp::CollapseShortEdges`
  inside `SimplifyTopology`.
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
| 4 | Edge-edge intersections | BVH broad phase over edge AABBs; trim-and-`Interpolate` symbolic kernel atop `CCW`/SoS. Two sub-phases: (a) pair-by-pair intersection insertion with snap-to-existing on the current pair's edges; (b) eager propagation: for each freshly-allocated intersection vert, query the edge BVH for *all* other edges within ε and add the vert to their on-edge vert lists. The propagation pass mirrors `boolean_result.cpp::AddNewEdgeVerts` from manifold's 3D boolean (where an edge × face intersection is added to three lists, not just the source pair). |
| 4b | Residual cleanup of near-duplicate intersection verts | When pair-by-pair processing produces near-duplicate intersection verts whose FP-rounded positions disagree by *more than ε* (typical at displaced coords with shallow crossings), step 4's eager propagation can't snap them. Union-find with structural (verts share an input edge) + geometric (within 10ε) gates collapses them. Mirrors `edge_op.cpp::CollapseShortEdges` in manifold's 3D `SimplifyTopology` cleanup. |
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
(default `maxIter=2`, Smith's bound) studies. ~2.3% of fuzz cases
differ between pass 1 and pass 2 at sub-ε position quanta but match
at the ε quantum (topology-identical, geometry shifted by a few
ULPs). A small fraction hit a more pathological mode: pass 1
produces a topologically-valid but geometrically *thin* output
polygon (all interior dimensions within 2ε), and pass 2 correctly
notices that the thin polygon's "long way around" coincides to
within ε with its "short way around" and collapses it to empty.
Consistent behavior under the algorithm's own rules, but means
pass-1 output is not always a stable fixed point. This is a *step
3* artifact (edge → on-edge vert detection): step 4's eager
propagation does not address it because the relevant edges have no
new intersection verts to propagate.

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

// Single-input regularization. Matches CrossSection::Simplify(eps)
// exactly: one input, one eps, returns the canonical wind > 0 boundary.
manifold::Polygons Simplify(const manifold::Polygons& in, double eps);

// Binary boolean. Uses manifold's existing OpType enum from
// `include/manifold/common.h` (the same one Manifold::Boolean and
// CrossSection::Boolean accept), so callers can pass the same value
// through 3D and 2D codepaths.
manifold::Polygons Boolean2D(const manifold::Polygons& a,
                             const manifold::Polygons& b,
                             manifold::OpType op,
                             double eps = 0.0);  // <=0 auto-infers

// Symmetric difference. Standard 4th binary boolean alongside
// Add/Subtract/Intersect; wiring target for `CrossSection::operator^`
// (which today delegates directly to Clipper2's Xor).
manifold::Polygons Xor(const manifold::Polygons& a,
                       const manifold::Polygons& b, double eps = 0.0);
}
```

The public surface mirrors `CrossSection`'s existing API exactly:
`Simplify(eps)` matches `CrossSection::Simplify(eps)`, `Boolean2D` is
the wiring target for `CrossSection::Boolean`, and `Xor` is the
wiring target for `CrossSection::operator^`. No fill-rule parameter
is exposed publicly, since `CrossSection`'s existing API has none.

Step 6's classification is internally a `WindPredicate` functor with
three named instances: `WindAdd` (w > 0, the default for `Simplify`
and Add/Subtract), `WindIntersect` (w > 1, used by Boolean2D's
Intersect path once both inputs contribute mult=+1), and `WindEvenOdd`
(w & 1, used by Xor). Subtract flips B's edge multiplicities to −1
before joining, so the same `WindAdd` predicate carves A − B. The
predicates stay internal because they depend on per-input
multiplicity conventions that the public entry points maintain.

## Validation

| Test | Pass rate | Notes |
|------|-----------|-------|
| Standard battery | 350/350 | 13 named cases + 350 random topological polygons. Topology-balance: every vertex's input edge contribution = output. Iter dist: 1:449 2:1. |
| Displacement fuzz | 450/450 | 5 scales (2<sup>10</sup>...2<sup>49</sup>) × 3 sizes × 30 seeds. |
| DeepFuzz (12,000 seeds, post-propagation) | **336,000 / 336,000 first-pass topology valid** | Iter dist: 1:328,495 2:7,359 3:139 4:6 5:1. 0 iteration-DEGRADED (pass-1 valid → broken by iteration). 1 geometric collapse (the thin-polygon idempotence quirk at `kPow=35 n=20 seed=8993` documented in Pipeline; survives propagation because it's a step-3 artifact, not k-fold concurrence). Compared to pre-propagation: **+592 cases promoted to iter=1**, fewer iter-2/3/4 cases at every level. |
| Polygons regularization (12,000 seeds) | 556 / 336,000 (0.17%) | The Polygons API regularizes output (Requicha-Tilove sense): zero-area lens loops (two oriented sub-edges between the same vert pair, tracing the same straight-line path in opposite directions) are dropped because `manifold::Polygons` can't encode 1D features without an enclosing face. Matches CGAL `Polygon_set_2`, Clipper2, and SVG fill-rule convention. Rate stable across the pre- and post-propagation 12k runs (560 → 556). A defensive assert verifies dropped loops always have zero area. |
| Boolean2D area regression | 4/4 | Two CCW unit squares offset diagonally; Add area=7, Subtract=3, Intersect=1, **Xor=6** (all exact). Drives the wind > 1 predicate against perpendicular axis-aligned crossings, plus EvenOdd via Xor. |
| Axis-aligned perpendicular | 1/1 | Two CCW unit squares offset diagonally; produces L-shape union with 8 boundary verts. Regression test for the kernel bug surfaced during OQ #3 / #4 investigation. |
| `polygon_corpus.txt` (100 entries) | 100/100 topology valid | Manifold's existing curated corpus, run via `./overlap2d_proto corpus`. Area breakdown: 87/100 < 1% drift, 1/100 mid-drift (`Precision4` at 2.6%, correct behavior on a sliver-with-eps-thin-tail), 7/100 collapsed by name (`Degenerate4`, `DegenerateRing`, `NearlyLinear`, `Looping1/2`, `Precision`, `Sweep`), 5/100 had zero-area input. |
| Clipper2 `Polygons.txt` corpus (195 entries) | 189/189 topology valid; 173/189 within 1% of current Clipper2 output | Independent third-party oracle, run via `./overlap2d_proto clipper2corpus`. 6 cases skipped (NEGATIVE/POSITIVE fill rule, no predicate analog). 173/189 (91.5%) within 1% of *current* Clipper2's output, 10 in 1-10%, 6 ≥10% — the ≥10% set are exclusively tiny intersections (areas 3-70 units) where 1-13 absolute units of FP error gives a large relative drift. Two findings emerged: (1) the corpus's stored `SOL_AREA` differs from current Clipper2 by ≥1% on 12 cases — it's stale, not a fresh oracle. (2) The `Boolean2D` API as written assumes CCW-canonical, single-orientation input (predicates use Smith's `w > 0` union convention). To match Clipper2's `boolean(fill_rule(subj), fill_rule(clip))` semantics on arbitrary inputs, callers must pre-fill each side with the declared fill rule. The corpus runner does this transparently; consumers wiring overlap2d as a `CrossSection` backend will hit the same need. |
| mfogel/polygon-clipping end-to-end corpus (134 case-op pairs across 87 named cases) | 134/134 topology valid; 112/134 within 1% of expected | Real-bug-derived JS corpus (run via `./overlap2d_proto mfogelcorpus`). 84% match within 1% after pre-fill under NONZERO with mfogel's first-ring-outer/rest-holes auto-normalization. The 2 cases ≥10% drift are both `coincident-with-invalid-segments` — mfogel-specific quirk on intentionally-malformed input. **Found a real framework bug**: chained `Boolean2D(Boolean2D(A, B), C)` mishandles intermediate T-junction vertices when the next input also has a coincident edge (`issue-44`). Fixed by switching the corpus runner's N-feature dispatch to a single-pass `RemoveOverlaps2D` call with combined inputs and a per-op predicate (`w != 0` for union, `w >= N` for N-input intersection, `w > 0` after mult-flipping for N-input difference, `w & 1` for xor). This is also more semantically faithful to N-input ops than left-to-right chaining. |
| Performance head-to-head vs Clipper2 (UNION+NONZERO) | overlap2d **~12x slower than Clipper2** on JTS overlay (the most realistic workload), ~14-50x on smaller corpora | Run via `./overlap2d_proto vsclipper2 [all\|clipper2\|jts\|mfogel]` (build with `-DOVERLAP2D_WITH_CLIPPER2` + link `libClipper2.a`). Per-corpus wall-clock for the boolean call only (excludes geometry conversion). On JTS overlay (2042 cases): overlap2d 1.45s parallel / 3.17s serial, Clipper2 0.13s, ratio 12x parallel / 19x serial — parallelism closes about a third of the gap. Clipper2 corpus (21 cases): overlap2d 0.085s, Clipper2 0.0015s, ratio ~50x; mfogel (62): ratio ~14x. Caveats: (1) Clipper2's `ClipperD` runs at default precision=2 (sub-0.01-unit truncation, integer arithmetic internally) — strictly less work than overlap2d's full-FP pipeline. (2) Clipper2 returned empty on 16 JTS cases that overlap2d preserved as non-empty (the `stmlf` snap-collapse — Clipper2's int rounding throws geometry away that overlap2d keeps; correctness going the other way). (3) Three concrete optimization avenues remain unexplored: arena allocation across cases (per-call malloc churn dominates small cases), an eps short-circuit when input is already CCW-canonical (skip the FillUnderRule pass entirely), and parallelizing step 4's intersection discovery (currently serial because of snap-as-you-go logic). Production landing should treat 5-10x as the realistic perf target after standard hardening, not the current spike-maturity number. |
| JTS robust/overlay corpus (2042 cases) | 2042/2042 topology valid; 1999/2034 invariant satisfied within 1% | Real-world GIS bug corpus from JTS Topology Suite (cases harvested from production tickets in GEOS, PostGIS, QGIS, JTS itself, Shapely). Run via `./overlap2d_proto jtscorpus` after generating the corpus locally with `python3 test/polygons/jts_to_corpus.py <jts_overlay_dir> test/polygons/jts_overlay_corpus.txt` (the 11 MB generated file is `.gitignore`d to keep the repo lean; the converter is checked in). The dominant op (~2034 cases) is `overlayAreaTest`, an invariant check `\|A∪B\| + \|A∩B\| = \|A\| + \|B\|` — a self-checking oracle; a bug in either union or intersection that creates or loses area shows up immediately. 1238 cases hold the invariant to within 1e-9, 761 within 1%, 35 violate by ≥1%. The 35 violators cluster in two datasets: `TestOverlay-stmlf.xml` (Bavarian land parcels at UTM-million coordinate scale, marked "all handled by snapping" in JTS — exposes the prototype's bbox-α-budget eps being too aggressive at large displacement, OQ #2) and `TestOverlay-isochrone.xml` (very-small-area intersections, ~1e-5 units, where FP noise dominates). XML+WKB+WKT parsing lives in `test/polygons/jts_to_corpus.py` (~150 LOC); the C++ runner just consumes a flat text format, keeping the algorithm code unencumbered by data-format plumbing. The same converter pattern absorbs other XML+WKB corpora (e.g. GEOS testxml) without C++ changes. |
| DeepFuzz area drift (12,000 seeds, post-propagation) | 406 cases > 1% (0.23% of 177,895 measurable), 145 > 10% | Quantitative oracle on top of topology. Of the cases with measurable input area, 0.23% drift > 1% between pass 1 and converged. Almost all cases at kPow=35 (where ULP precision is weakest); a few at kPow=30. Max drift 100% (the seed=8993 thin-polygon collapse). The 100-seed sample showed 6 such cases; the 12k run reveals a longer thin-polygon-class tail invisible at smaller scale. |

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
3. **Multi-threading.** Resolved. The prototype builds with either
   `MANIFOLD_PAR=1` (TBB) or `MANIFOLD_PAR=-1` (serial) and parallelizes
   the data-parallel hot loops via `manifold::for_each(autoPolicy(N), ...)`
   plus `DisjointSets`'s lock-free atomic `unite()` /  `find()`.
   Parallelized phases:

   - **Step 6 face-winding ray-cast** (per-face independent; was 82% of
     pipeline time, now 33-50%).
   - **Step 6 per-vertex angular sort** (per-vertex independent).
   - **Step 6 next-pointer computation** (per-halfedge independent).
   - **Step 1 vertex-merge** geometric gate (parallel filter, serial
     `unite()` in sorted pair order for deterministic cluster roots).
   - **Step 4b structural re-merge** geometric+structural gate (same
     parallel-filter / serial-unite pattern as step 1).

   `autoPolicy()`'s `kSeqThreshold = 1e4` keeps small inputs serial so
   there's no overhead penalty on per-call latency. Determinism is
   preserved on all 4 corpora (Clipper2/mfogel/JTS/polygon_corpus.txt)
   bit-identically.

   Step 4 (`FindAndInsertIntersections`) was not parallelized: the
   snap-to-existing logic mutates `lists[*]` and reads it on every
   subsequent pair, making the algorithm inherently sequential without
   a structural rewrite (e.g., a parallel intersection-discovery phase
   followed by a separate snap-merge pass).

   An orthogonal serial optimization landed in the same patch: the
   per-face ray-cast was switched from iterating `canon.map` (a std::map
   RB-tree) to a pre-flattened `FastEdge` array (`yMin`, `yMax`, `x0`,
   `dxdy`, `signedMult`) built once before the per-face loop. This is
   what most of the wall-clock improvement is — the parallel `for_each`
   on top adds an additional ~1.3x on big workloads.

   **Speedups on the standard workloads** (8-core machine, all corpora
   produce identical output to the serial baseline):

   | Workload | Original serial | Serial + FastEdge | Parallel + FastEdge |
   |---|---|---|---|
   | DeepFuzz 200 (22.7k cases, mostly small) | 71.2 s | 25.1 s (2.84x) | 24.5 s (2.91x) |
   | JTS overlay (2042 cases, mixed sizes) | 36.95 s | 4.86 s (7.6x) | **3.6 s (10.3x)** |

   Most JTS cases have N < 1e4 verts and stay serial under
   `autoPolicy()`; the wall-clock win comes from the few large cases
   (e.g. `TestOverlay-osmwater.xml`, `stmlf.xml`) where step 6's F faces
   exceed the threshold and parallelize ~4-5x. DeepFuzz cases are all
   small enough to stay serial — its speedup is purely from FastEdge.

   Per-phase wall-clock instrumentation lives behind a `time` CLI mode
   (`./overlap2d_proto time deepfuzz | jtscorpus | clipper2corpus | ...`)
   that prints accumulated nanoseconds per pipeline phase plus a
   percentage breakdown.
4. **Additional quantitative oracles.** The prototype consumes
   `polygon_corpus.txt` (100 cases) and Clipper2's `Polygons.txt`
   (195 cases; 189 reachable) — see Validation. The Clipper2 run
   surfaced one design item worth flagging: the `Boolean2D` API
   currently assumes CCW-canonical input (predicates `w > 0` /
   `w > 1`), so reproducing Clipper2's `boolean(fill_rule(subj),
   fill_rule(clip))` semantics on arbitrary fill-rule input requires
   each side to be pre-filled under its declared fill rule first.
   The corpus runner does this with a `FillUnderRule` helper; whether
   `CrossSection`'s wiring should expose a fill-rule parameter (or
   always auto-canonicalize internally) is the production-landing
   question.

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
  bound is 2; tested up to 8. The 12,000-seed DeepFuzz (post-
  propagation) finds 139 cases at iter=3, 6 at iter=4, 1 at iter=5
  (out of 336,000), all ultimately converging to topology-valid
  output. Bumping the cap catches more fingerprint-stable convergence
  but doesn't change topology, since pass 1 was already topology-
  correct on every case. Production default of `maxIter=1` skips this
  entirely (see Pipeline section); the prototype's `maxIter=2` matches
  Smith's bound and is the right knob for studying the tail.
