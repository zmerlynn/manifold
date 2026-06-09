# Boolean engine emits sliver triangles at near-coincident face boundaries

> **Filed upstream as elalish/manifold#1706 (May 2026).**
>
> Maintainer (Emmett Lalish) responded that **slivers are expected
> output of Smith's Boolean** due to symbolic perturbation —
> attempting to snap or remove them inside the Boolean kernels would
> risk breaking the manifold-ness guarantee that SoS provides. The
> intended cleanup mechanism has always been
> `Manifold::Impl::SimplifyTopology` (called from
> `Manifold::Simplify(tolerance)`), but it's currently too
> conservative (recently tightened to avoid stacking FP errors that
> were collapsing legitimate high-res geometry like spheres into
> garbage), so it under-cleans on inputs like gt-7081.
>
> **Agreed-on fix path: replace `SimplifyTopology` with a proper
> QEM-based decimator** following Hoppe's
> [New Quadric Metric for Simplifying Meshes With Appearance
> Attributes](https://hhoppe.com/newqem.pdf), structured as a
> semi-parallel recursive algorithm:
>
>   1. Compute the QEM metric for unprocessed edges in parallel
>      (initially: all edges).
>   2. Sort by metric, keep only those below the desired threshold.
>   3. Run `CollapseEdge` (with vert update — see below) serially,
>      smallest-metric first, marking endpoints visited; skip any
>      edge whose endpoints are already visited.
>   4. Reset visited and recurse with the skipped queue.
>
> Key implementation note from Emmett: the existing `CollapseEdge`
> doesn't move the surviving vert to a new position, which is why
> `SimplifyTopology`'s repeated passes were stacking FP error. A
> QEM-aware collapse moves the survivor to the position that
> minimizes `pᵀ(Qa + Qb)p`, keeping cumulative geometric error
> bounded across many collapses.
>
> Maintainer is planning to do this work himself; deferred from this
> branch.

## Symptom

When `Manifold::Boolean` operates on inputs whose surfaces nearly
coincide in a region (= "twin" geometry, common in CSG kit-bashing
and offset-Boolean compositions), the result mesh contains very
small ("sliver") triangles along the chord between the near-
coincident faces. These slivers have area below FP noise (typically
`10^-9` to `10^-14` square units on a normal-scale mesh) but are
topologically wired into the mesh, so they can't simply be dropped
without leaving non-manifold holes.

For most downstream consumers the slivers are tolerated (= still a
valid Manifold by `tolerance_`), but they sabotage any algorithm
that traverses the mesh via per-vertex geometric probes or
per-polygon enumeration. The current case study is
`Manifold::RemoveSelfIntersections` (see
`docs/gt7081-shape-analysis.md`), where 615 slivers in
`Generic_Twin_7081 (m1+m2)` overwhelm the polygon walker and
pair-symmetric chord enforcement.

## Reproduction

```cpp
#include "manifold/manifold.h"
#include <fstream>
using namespace manifold;

std::ifstream a("test/models/Generic_Twin_7081.1.t0_left.obj");
std::ifstream b("test/models/Generic_Twin_7081.1.t0_right.obj");
Manifold m1 = Manifold::ReadOBJ(a);
Manifold m2 = Manifold::ReadOBJ(b);
Manifold result = m1 + m2;
// result.NumTri() == 33230
// 615 of these triangles have area < 1e-9 (mesh bbox is ~25400 units)
// 21 of those have area < 1e-13 (= edge length ~3e-7, near FP noise)
```

Tri area histogram (log10 buckets):

```
10^-13:   21    <-- worst slivers
10^-12:  181
10^-11:  258
10^-10:  111
10^-9:    44
10^-8:    18
...
10^0:  31051   <-- bulk of the mesh
```

## Why slivers appear

The Boolean engine's chord assembly path:

1. `src/boolean3.cpp:Intersect12_<true, forward>(inP, inQ)`: finds
   each edge of P that pierces a face of Q. For each piercing, the
   `Kernel12` SoS kernel computes a 3D intersection vertex `v12`
   (line `boolean3.cpp:184`). Direction sign goes into `i12`,
   intersection point into `v12`.
2. `src/boolean_result.cpp:AddNewEdgeVerts` (line 213) distributes
   each new vertex into three `EdgePos` bins:
   - The original edge of P (= it gets subdivided at v12).
   - The two new face-face chord edges (one per face of P
     adjacent to the source edge).
3. `src/boolean_result.cpp:AppendPartialEdges` (line 305) and
   `AppendNewEdges` (~line 385) sort each bin by `edgePos` (= the
   parametric coordinate along the edge) and emit sub-edges
   between consecutive entries.
4. The triangulator (`src/face_op.cpp` and friends) emits triangles
   using these sub-edges as boundaries.

The pathology arises in step 3: when an edge of P is **nearly
tangent** to face Q (= the input surfaces nearly coincide along the
chord), `Kernel12` produces a `v12` that lies very close to one of
the existing endpoints of edge P. After sorting by `edgePos`, this
new vertex appears between the existing vertex and its true
neighbor, with a `edgePos` delta below FP noise. The triangulator
emits a triangle that includes this near-zero sub-edge as one of
its sides → sliver triangle (= long-thin or flat).

This is **by design**: SoS guarantees that distinct piercings
produce distinct vertices even when they're geometrically near-
identical, so manifoldness is preserved through the kernel cascade.
The slivers are the visible cost of that guarantee. Cleaning them
post-hoc (= the QEM decimator path above) is the right layer.

Slivers are particularly common when:

- Two inputs share a face that's been shifted/rotated by FP-noise
  amounts (= "almost aligned" CSG operands).
- An offset-Boolean composition where the offset distance is much
  smaller than the input feature scale.
- Twin geometry (= two near-identical shapes added or intersected).

## Why existing simplification doesn't catch them

Per Emmett's response on #1706:

> `SimplifyTopology` was designed precisely for this, but it's not
> doing a great job currently. It used to be better but could stack
> FP-size errors causing e.g. a high-res sphere to collapse into
> garbage. I fixed that, but the net result is we're now too
> conservative about collapsing edges.

Specifically, the existing `CollapseEdge` doesn't move the surviving
vert to a new position — it keeps one of the endpoints. That means
each collapse contributes geometric error bounded only by the
endpoint distance, and chained collapses can drift the surface
arbitrarily far from the original. The recent tightening prevents
that drift but at the cost of leaving real slivers in place.

We confirmed empirically:
`Manifold::Simplify(0.001)` on gt-7081 produces 737 slivers
(more than the 615 in raw output) because the edge-length-keyed
collapses generate fresh sliver tris. Larger tolerances destroy
legitimate mesh structure (e.g., `Simplify(1.0)` drops vol from
143611 to 109868). No existing tolerance value cleans gt-7081
without geometric distortion.

## Pre-existing defense in `RemoveSelfIntersections`

None today. The downstream `Manifold::RemoveSelfIntersections`
pipeline is a heuristic family that local-equilibrium-converges
on small sliver counts (= self-intersect with 29 slivers, hull-mask
with 4) but cannot handle gt-7081's 615. There is no sliver-aware
pre-process in `src/overlap_removal.cpp` — pierce-aware Phase 3
helps but doesn't reach manifoldness on gt-7081, so the pierce/
drift gate falls back to input.

If the upstream QEM decimator is delayed, a defensive sliver-area
pre-pass at the entry of `RunOverlapRemovalImpl` could be added,
but it would be a workaround for the symptom rather than the cause,
and it would re-implement what the upstream decimator is meant to
provide.

## Verification once the QEM decimator lands

1. `Manifold::Simplify(tolerance)` on `Generic_Twin_7081 (m1+m2)`
   should eliminate the 615 slivers without measurably changing
   the bulk geometry (= total volume drift < 1%, max vert
   displacement bounded by `tolerance`).
2. `Cray_left - Cray_right` is a separate pathology
   (= polluted `tolerance_`, see `docs/tolerance-pollution-bug.md`)
   and should be unaffected by the QEM work.
3. `Manifold::RemoveSelfIntersections()` on a Simplify'd
   gt-7081 should produce a non-passthrough result. It may not
   eliminate all 183 piercings (= the underlying input may have
   real self-intersections beyond the sliver region), but the
   pipeline should at least produce a NoError manifold whose
   pierce count is ≤ input.
4. Existing
   `test/boolean_complex_test.cpp:GenericTwinBooleanTest7081`
   continues to pass (= the Boolean op itself doesn't crash, and
   `result.GetMeshGL()` works).

## History

- **Issue #741 / PR #742 (Feb 2024, pca006132)**: "SparseIndices
  too large" core-dump fix. The fix made the Boolean engine survive
  this input without crashing; the sliver output was not addressed
  because crashes were the immediate symptom.
- **commit 7f79cb59 (May 2026, this branch)**: Phase 3 added to
  RemoveSelfIntersections to handle small sliver counts. Reduces
  but does not eliminate the gt-7081 unpaired residual.
- **commit 8f56cd45 (May 2026, this branch)**:
  `docs/gt7081-shape-analysis.md` — companion analysis identifying
  the slivers as the root cause.
- **Issue #1706 (May 2026)**: filed against elalish/manifold;
  Emmett confirmed slivers are by design and outlined the QEM
  decimator path as the intended fix.

## Out of scope

- The "real" geometric question (= are gt-7081's two operands
  actually meant to be perfectly-coincident in the overlap region?)
  is between the user and the input data. The Boolean engine
  cannot tell intent from FP-noise-sized differences. The fix is
  purely about cleaner output, not about second-guessing intent.
- This bug and the tolerance-pollution bug
  (`docs/tolerance-pollution-bug.md`) are independent: tolerance
  pollution affects `Cray_*` (= `FLT_MAX`-scale operand pollutes
  output `tolerance_`); slivers affect `Generic_Twin_7081` (=
  near-coincident operands produce sliver tris). Both manifest in
  the Boolean engine's output, but the fixes live in different
  code paths.
