# gt-7081 root cause: 615 sliver triangles from twin-overlap Boolean

## Summary

`Generic_Twin_7081.1.t0_left + .t0_right` is a "twin" Boolean Add
where the two operands share a small overlap region. Boolean Add
produces 615 sliver triangles (area < 1e-9, on a mesh with bbox
~25400 units) along the near-coincident face boundaries in the
overlap region. These slivers overwhelm the polygon walker and
pair-symmetric chord enforcement in `RemoveSelfIntersections`,
producing 1500 unpaired halfedges that the cap walker can't close.
Pipeline falls back to input.

This is fundamentally different from cray (which has a polluted
upstream tolerance from a `FLT_MAX`-scale Subtract operand). gt-7081
is a real-shape pathology, not a metadata bug.

## Shape characteristics

```
m1 (left.obj):  v=9874  t=19704 vol=101591  bbox=[8250, 14723, 1845]
m2 (right.obj): v=5830  t=11656 vol= 43004  bbox=[4382,  8305,  845]
m1 ^ m2 (intersection): vol=984  (= 1% of m1, status NoError, t=2934)
m1 + m2 (add):  v=16635 t=33230 vol=143611
```

m1 and m2 are "twin" shapes whose bboxes overlap and whose volumes
intersect by ~1%. The intersection region (vol 984) has near-
coincident faces between m1 and m2 — the Boolean Add chord runs
along these faces and produces sliver triangles.

## Sliver distribution

Tri area histogram (log10 buckets):

```
bucket  | count  | comment
--------|--------|--------
10^-13  |   21   | severe sliver (area ≈ 1e-13, edge ≈ 3e-7)
10^-12  |  181   |
10^-11  |  258   |
10^-10  |  111   |
10^-9   |   44   |   <-- threshold for "sliver" in our analysis
10^-8   |   18   |
10^-7   |   42   |
10^-6   |   51   |
10^-5   |   66   |
10^-4   |  114   |
10^-3   |  710   |
10^-2   |  214   |
10^-1   |  222   |
10^0    | 31051  | <-- bulk of the mesh, well-formed unit-area tris
10^1    |   52   |
10^2    |   45   |
10^3    |   30   | largest tri area ≈ 4914
```

Total slivers (area < 1e-9): **615**. Ratio largest/smallest area:
**1.24e+17** (= 17 orders of magnitude). For comparison:

| fixture | tri count | slivers (<1e-9) | min area |
|---|---|---|---|
| self-intersect | 33542 | **29** | 9.4e-20 |
| hull-mask | 13974 | **4** | 6.4e-10 |
| **gt-7081** | 33230 | **615** | 4.0e-14 |
| gt-7863 | 940 | **3** | 5.6e-13 |

gt-7081 has 21x more slivers than the next-worst fixture
(self-intersect) at the same tri-count scale.

## Why slivers break the pipeline

The pipeline in `src/overlap_removal.cpp` runs the Emmett #289 13-step
process: enumerate edges, find tri-tri intersections, build per-tri
halfedge graphs, walk polygons, classify per polygon (two-sided
winding), enforce pair-symmetric chord constraints, triangulate, cap.

Slivers sabotage this in several ways:

1. **Polygon walker.** A sliver tri's edges enter the per-tri
   halfedge graph. If the chord intersects the sliver, sub-polygons
   are emitted that are themselves degenerate. Downstream classifier
   and pair-sym phases process these as if they were real polygons.

2. **Pair-sym chord enforcement.** Phase 1 enforces "exactly one
   direction kept per chord," picking arbitrarily when ties occur.
   Slivers create many chord halfedges with ambiguous pair status
   (= the polygon walker emitted both directions in degenerate
   sub-polygons). Phase 1 picks arbitrarily; the pick is often
   wrong for slivers.

3. **Phase 2 AND-cascade.** Single-owner non-chord pairs that
   disagree on keep state get dropped together. For gt-7081 this
   cascades to ~3375 polygon drops (mostly from Phase 1's drops
   propagating through halfedge pairings).

4. **Phase 3 strict-majority equilibrium.** After Phase 1/2 drops,
   Phase 3 drops kept polygons with strict-majority unpaired
   halfedges. Equilibrium leaves 1500 unpaired residuals — the
   sliver-induced drops created so much disagreement that no local
   move can reach a manifold subset.

5. **Cap walker exhaustion.** The recovery sweep at the end of
   `TriangulateAndEmit` tries to close k=1 cycles. For 1500 unpaired
   it can't converge.

## Why our 5-minute runtime

Each phase is roughly O(N) or O(N log N) in polygon count, but the
sliver-induced complexity multiplies the constant factor:

- Polygon walker emits ~34000 polygons (vs ~33000 input tris) but
  with very different size distribution — many tiny polygons
  triggered by slivers.
- Phase 3 iterates 10 times to reach equilibrium (vs 6-7 for
  working fixtures of similar size).
- Phase 3.5 re-key pass adds another ~4 iterations.
- Recovery trim + cap iterates 8 more times.
- Each iteration scans all 34000 polygons.

Working fixtures of similar tri count (self-intersect: 33542 tris)
finish in ~2 seconds. gt-7081 takes ~5 minutes. The 150x slowdown is
the sliver-induced iteration count, not asymptotic complexity.

## What was tried this session that did NOT help

- **Phase 3 strict-majority drop**: reduced unpaired 2051 → 1539
  (initial).
- **Phase 3.5 re-key dropped polygons**: reduced 1539 → 1468.
- **Loose Phase 3 (drop any polygon with unpaired halfedges)**:
  REGRESSED self-intersect (661 → 0 → 11). Math: dropping a
  triangle with k=1 unpaired adds (size − 2k) = 1 NEW unpaired.
  Triangles can't be dropped at the strict-majority equilibrium
  without making things worse.
- **Re-keying with relaxed duplicate constraint**: small additional
  benefit (1499 → 1468 unpaired).
- **`Manifold::Simplify(tol)` with various tolerances**:
  Simplify(0.001) collapses short edges and creates MORE slivers
  (615 → 737). Simplify(1.0) destroys ~24% of volume. Simplify
  doesn't target zero-area tris with non-short edges (= the
  long-thin sliver case that dominates gt-7081).
- **Drop-slivers-only pre-process**: dropping the 615 area-<1e-9
  triangles produces a NotManifold MeshGL64 (= leaves 615 holes).
  The cap walker can't close that many simultaneously.

## What might actually solve gt-7081

Three plausible approaches, increasing in scope:

1. **Sliver-collapse pre-process**: for each sliver tri, identify
   the two near-collinear vertex pairs and merge them (= edge
   collapse, NOT just drop). Preserves topology. Implementation:
   call `Manifold::Impl::CollapseEdge` for the shortest edge of
   each sliver. Existing internal API in `src/edge_op.cpp:521`.
   Need to thread sliver detection into pipeline entry and call
   `SimplifyTopology` with a custom predicate that targets area,
   not edge length.

2. **Make Boolean op produce fewer slivers**: upstream fix in
   `src/boolean3.cpp` / `src/boolean_result.cpp`. When chord
   intersection produces a sub-edge with length below FP noise,
   collapse it into the existing endpoint instead of emitting a
   degenerate triangle. This is the right architectural fix but
   requires deep familiarity with the Boolean engine's invariants.

3. **Global SAT/ILP for kept-polygon manifold-pair selection**:
   model the problem as integer programming, solve to find
   minimal-cost kept subset. Exponential worst-case but tractable
   for ~30k polygons with a good solver. Would replace Phase 1, 2,
   2.5, 3 with a single optimization pass. Significant
   implementation cost.

## Comparison with the spike

The spike's default-mode pipeline ALSO fails on gt-7081 (status 2
NotManifold output, 2051 unpaired k=1 edges). Our prod pipeline
adds Phase 3 + 3.5 which reduces unpaired to 1500. So this is not
a porting bug — it's a fundamental limitation of the heuristic
family. Both pipelines fall back to input.

## History

- **PR #742 (Feb 2024, pca006132)**: "Fix sparseindices." Closed
  issue #741. The fix addressed a `SparseIndices` integer overflow
  triggered by the boolean engine on this fixture. The crash was
  fixed; the resulting mesh's 615 slivers were not addressed — they
  are within `tolerance_` and considered valid output by Manifold.
- **commit 7f79cb59 (May 2026, this branch)**: Phase 3 added.
  Reduced unpaired residual but doesn't reach manifold.
- **commit 2f8b7335 (May 2026, this branch)**: Phase 3.5 added.
  Marginal improvement.

## Recommendation

For shipping, mark gt-7081 as known-pathological alongside cray.
The fixture remains useful for regression testing (= the pipeline
must not crash, must fall back cleanly). Pierce-monotonicity is
preserved (output pierces ≤ input pierces, by virtue of the gate).

For deep work, **approach 2 above (Boolean op slivers)** is likely
the highest leverage. The 615 slivers should not exist in a clean
Boolean output. Fixing them in the engine benefits any downstream
consumer, not just RemoveSelfIntersections. Approach 1 (sliver-
collapse pre-process) is a defensive measure that could go in
RemoveSelfIntersections but doesn't address the root cause.
