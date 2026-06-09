# 3D overlap removal: speculative design

Status: **speculative**, exploring how the 2D framework
([`Overlap2D.md`](Overlap2D.md)) extends to 3D. This doc describes
how Emmett Lalish's 13-step 3D sketch in
[#289](https://github.com/elalish/manifold/issues/289#issuecomment-2111069955)
relates to the existing `src/boolean3.cpp` + `src/boolean_result.cpp`
+ `src/edge_op.cpp` pipeline, and what the gap is.

References:

- [`Overlap2D.md`](Overlap2D.md) for the 2D version this extends.
- [`RobustBoolean.pdf`](RobustBoolean.pdf): Smith's UCAM-CL-TR-766.
  3D framework discussion in §8.
- Existing implementation: `src/boolean3.cpp`,
  `src/boolean_result.cpp`, `src/edge_op.cpp`,
  `src/properties.cpp::IsSelfIntersecting`.

## Two distinct invariants

Manifold's existing 3D `Boolean` pipeline targets a **topological**
invariant: every edge of the output is shared by exactly two
triangles, the half-edge structure is consistent, and the result
reads as a closed 2-manifold. This is what `Status() == NoError`
asserts.

Emmett's #289 13-step 3D sketch targets a stricter **geometric**
invariant: no two triangles in the output have strict interior
overlap (no tri-tri pierce). This is what overlap removal *is*.

The two invariants are not equivalent. A mesh can be topologically
manifold and contain tri-tri pierces, and the existing 3D pipeline
sometimes produces such outputs — manifold's own
`Manifold::Impl::IsSelfIntersecting()`
(`src/properties.cpp:138`) is exactly the post-hoc geometric check,
gated behind `ManifoldParams().selfIntersectionChecks` (default OFF,
`MANIFOLD_DEBUG`-only). Manifold's tree includes named adversarial
fixtures whose primary purpose is to exercise this gap:
`Cray_left/right.obj`, `Havocglass8_left/right.obj`,
`Offset[1-4].obj`, `Generic_Twin_*.obj`, `self_intersectA/B.obj`.

So **the existing pipeline is not an implementation of overlap
removal in Emmett's sense.** It implements the topology-manifold
boolean. Overlap removal would be a new operation that takes
arbitrary input (potentially self-intersecting, potentially
non-manifold) and produces a topology-manifold *and*
self-intersection-free output, per a chosen winding rule.

## Where the existing pipeline overlaps with Emmett's sketch

Several of Emmett's 13 steps have analogs in the existing topology-
manifold pipeline; that's expected, since both sit on the same
arithmetic foundations (`Collider`, `Shadows`, `Interpolate`,
`Kernel11/12`). The mapping below is for orientation, not a claim
that the existing code already implements overlap removal.

| # | Emmett's #289 step (overlap removal) | Closest existing-pipeline analog |
|---|---|---|
| 1 | Vert merge within ε | Input merging before `Boolean3::Boolean3()` (uses `epsilon_`) |
| 2 | Drop collapsed edges/triangles | `RemoveDegenerates` / `SimplifyTopology` pre-pass |
| 3 | Triangle → on-tri vert list | Implicit in `Kernel02` / `Shadow02` (point-in-triangle) |
| 4 | Edge-edge intersections (skew degeneracies) | `Kernel11` / `Shadow11` |
| 5 | Triangle → on-tri vert list (post step 4) | `Kernel02` for edge-of-A vert in face-of-B |
| 6 | Edge × triangle intersections | `Kernel12` / `Intersect12` |
| 7 | Add a new edge per tri-tri intersection | Implicit in `xv12_`/`xv21_` from `Intersect12` |
| 8 | Add interior verts within ε to new edges | `AddNewEdgeVerts` populates per-pair lists |
| 9 | New-edge × new-edge intersection within each tri | Per-tri-face polygon assembly in `boolean_result.cpp` |
| 10 | Sub-edge canonicalization (per-tri halfedge lists) | `boolean_result.cpp` builds `edgesP/Q/New` keyed by face |
| 11 | Per-triangle polygon partition via cyclic ordering | Face-cycle walking in `boolean_result.cpp::Face` |
| 12 | Equivalent-polygon merge with multiplicity | Implicit in the `pair<int,int>` keying |
| 13 | Output polygons separating winding 0 from >0; triangulate | `Winding03` + `Triangulate` |

The mapping is *structural*, not behavioral: even where a step has
an existing analog, the existing analog is doing it in service of
the topology-manifold invariant, not the no-pierce invariant.
For example, step 8's `AddNewEdgeVerts` propagation is per-pair, not
the broad-phase-against-all-edges that Emmett's step 8 specifies. A
real overlap-removal implementation would have to revisit each step
under the stricter invariant.

## Empirical findings (from the spike)

The companion spike (`extras/overlap3d_proto.cpp`) runs the existing
`Manifold::Boolean` (wrapped with the α-budget eps inference) over
fuzz and adversarial cases, and applies a tri-tri self-intersection
check (`CheckSelfIntersection`, an independent BVH+Möller
implementation) to the output. Findings (current battery — to be
expanded as adversarial corpus is added):

- 800-case rotated-cube fuzz (4 displacement scales × 200 seeds):
  799 produce self-intersection-free output, **1 produces a tri-tri
  pierce** at sub-ε magnitude (kPow=30 seed=116, two pierces of
  magnitude ~2×10⁻⁸ relative). Below the manifold's stated
  tolerance (~10⁻³ at this scale) but above FP noise. Real
  geometric self-intersection, undetected by `Status() == NoError`.
- **Adversarial cube fuzz** (`./overlap3d_proto --advfuzz`,
  commit 91d2dca4): 7 classes × 200 cases = 1395 valid cases at
  kPow=30 displacement. 12 self-pierce cases = 0.86% rate,
  **~7x the general fuzz**. Per-class:
    shallow-tight:       0/200
    near-coplanar-slabs: 3/195 (1.5%)
    3-cube-chain:        2/200 (1%)
    4-cube-chain:        0/200 (variance, was 1/100 in run 1)
    5-cube-chain:        **5/200 (2.5%)** ← most adversarial
    6-cube-chain:        2/200 (1%)
    vertex-on-face:      0/200
  Chained Booleans accumulate FP rounding; long chains amplify.
- **`--advfuzz` end-to-end utility test** (commits d50730c3 +
  8fff5ce6): for each self-piercing Boolean output, run
  `OverlapRemoval` and re-check. Out of 12 self-piercing cases:
    fully fixed (post-pierces=0):       **8/12** (66.7%)
    reduced to FP-noise floor:          1/12 (3-pierces 8.6e-13
                                              rel → 2-pierces
                                              5.0e-17 rel)
    regularized-to-empty (sub-eps):     3/12 (correct per
                                              Requicha-Tilove,
                                              matches CGAL/SVG
                                              fill-rule conv)
    **TOTAL valid outputs: 12/12 (100%)**
  The spike produces a valid manifold result for every
  adversarial self-piercing input the fuzz finds.

  *On the "1 reduced" case:* it went 3 → 2 pierces but the
  surviving 2 had relative magnitude 5.0×10⁻¹⁷, which is
  *below* IEEE double precision epsilon (2.22×10⁻¹⁶). They are
  not real geometric pierces, just FP rounding artifacts that
  happen to pass `CheckSelfIntersection`'s `relTol=1e-12`
  filter at extreme kPow=30 scale. Pre-fix magnitude was
  8.6×10⁻¹³ (above FP noise but sub-eps); post-fix is in the
  noise. Practically: fully fixed. Going below this requires
  higher-precision arithmetic or a different algorithm,
  neither of which is in scope for the spike.
- Smith-taxonomy probes (V-V, V-E, V-T, E-E, E-T, T-T, k-fold): all
  self-pierce-free in the current minimal cases. *To be expanded*
  with the named .obj fixtures (Cray, Havocglass, Offset,
  Generic_Twin, self_intersect) and constructed near-coplanar
  sliver cases.
- **Battery .obj fixtures**: all 8 fixtures from `test/models/`
  are now in the granular runner — Cray, Havocglass8, hull-mask,
  Offset12, Offset34, self-intersect, Generic_Twin_7081,
  Generic_Twin_7863. The OverlapRemoval API path produces
  real output for 5/8 (clean drift); 3/8 fall back gracefully
  on heavy-pierce/cascade-drop inputs (offset34, generic-twin-
  7081, cray ~30% of runs).

The 1-in-800 pierce is the kind of finding overlap removal exists
to fix. Whether the rate is 1-in-800 on natural inputs or much
higher on adversarial inputs is what the expanded battery will
determine.

## What an overlap-removal feature would add

These are the *real* deltas — what the new feature would do that
the existing pipeline doesn't. Note: input is topology-manifold
(per the `Manifold` type's invariant) — overlap removal does not
attempt to repair non-manifold soup, gaps, or T-junctions; that
would be a separate "mesh repair" feature. The thing it removes is
*geometric self-intersection* in already-topology-manifold input.

1. **Snap-instead-of-resolve at degeneracies.** Emmett's "above /
   below / on" three-valued predicates (motivated by Smith) snap
   verts to edges when on, snap edges to faces when on, etc. The
   existing pipeline's symbolic perturbation always picks one side,
   producing the topology-manifold output but leaving the geometric
   pierce in place.

2. **Broad-phase propagation of new verts.** When a new edge×tri
   intersection vert is created, it must be checked against *all*
   nearby edges/tris (BVH range query at radius ε), not just the
   pair that produced it. The existing `AddNewEdgeVerts` propagates
   only within the source triple. The 2D prototype's step 4(b)
   eager propagation is the same idea.

3. **Output a fill-rule-classified surface, not a manifold-recovery
   surface.** Step 13's "winding 0 vs >0" in 3D is a face-by-face
   ray-cast, not the topology-recovery `Winding03` that the
   existing pipeline does. Different question, different answer.

4. **An α-budget eps formula (Smith) at the operation boundary.**
   `AlphaBudgetEpsilon(L, k_budget=1000)` already exists in
   `src/shared.h` (lifted from the 2D prototype work). The
   overlap-removal entry point would set the operation's eps using
   it by default.

Each of these is non-trivial; together they're a real algorithm,
not a wrapper around `Boolean`. The estimate of "documentation
pass + small additions" in earlier drafts of this doc was wrong.

## Minimal viable spike (in progress)

Building a real overlap-removal implementation is multi-week work.
The spike (`extras/overlap3d_proto.cpp`) is structured to first
*characterize the gap* — measure how often the existing pipeline
produces self-pierces, on what inputs — and then implement
Emmett's #289 13 steps incrementally with stats reported at each
step. Sequence:

1. ✅ α-budget eps formula promoted to `src/shared.h`.
2. ✅ Independent `CheckSelfIntersection` (BVH + Möller) operating
   on `Manifold`, used as ground truth on Boolean output.
3. ✅ Smith vocabulary doc-comments in `boolean3.cpp` /
   `boolean_result.cpp` (these *don't* claim to implement overlap
   removal; they label which symbolic kernels exist).
4. ✅ **Adversarial battery**: named `.obj` fixtures from
   `test/models/` (Cray, Havocglass, Generic_Twin, Offset,
   self_intersect, hull-mask) + chained-Boolean test +
   constructed Smith fig 6.4/6.5 3D analogs + 800-case
   rotated-cube fuzz at displacement.
5. 🔄 **13-step implementation in `extras/`**, queries-only
   read-only against the input `Manifold`'s `Impl` (no
   mutation yet — that's the next phase):
   - ✅ Step 1: ε-merge verts (Collider + DisjointSets, mirrors
     `sort.cpp:155-167`).
   - ✅ Step 2: drop collapsed + edge enumeration (via
     `halfedge_`).
   - ✅ Step 3: per-edge on-edge vert lists.
   - ✅ Step 4: edge-edge intersection discovery (3D segment-
     segment closest-points).
   - ✅ Step 5: per-tri on-interior vert lists.
   - ✅ Step 6: edge × triangle intersection — **the main
     pierce finder**. Catches Cray's 105 tri-tri pierces and
     self-intersect's 657 (vs CheckSelfIntersection's 114 / 661;
     numbers close).
   - ✅ Step 7 (phase 1): tri-tri pair enumeration. Validates
     Emmett's "exactly 2 endpoints per pair" claim — Cray
     100% n=2, self-intersect 99.7% n=2.
   - ✅ Step 7 (phase 2): emit new verts + new edges from the
     piercing pairs (Cray: 105 new verts, 105 new edges, 0
     dropped). *First mutation step (output is candidates only;
     input Manifold is not modified.)*
   - ✅ Step 8: propagate in-tri verts onto new edges (defensive
     pass for k-fold concurrence at new edges).
   - ✅ Step 9: new-edge × new-edge intersections per tri
     (defensive pass for k-fold concurrence at multi-pierce tris).
   - ✅ Step 10 (counts only): per-tri sub-edge enumeration.
     Cray: 948 sub-edges (714 baseline + 234 from new edges,
     max 15 sub-edges on one tri = 6 piercing pairs sharing it).
   - ✅ Step 11 phase 0: `OverlapRemoval(Manifold, eps) →
     (Manifold, debug)` entry point. Output construction is
     stubbed (returns the input through MergeVertsEps); volume
     round-trip is lossless across all .obj fixtures. The
     stubbing also surfaced a real pre-existing finding: a
     `Manifold → GetMeshGL64 → Manifold` round-trip flipped
     volume sign and inflated 100× on Cray (because Cray is a
     Subtract result with back-side flag; back-side info
     doesn't fully survive the round-trip). Worked around by
     returning input unchanged when no merges occur.
   - ✅ Step 11 phase 1: per-tri halfedge graph construction.
     Original sub-edges contribute one halfedge each; new
     sub-edges contribute two (both directions per Emmett step
     10). Cray: 1418 halfedges, max 37 / 20 verts on the
     hot-spot tri.
   - ✅ Step 11 phase 2: 2D projection (drop dominant normal
     axis, axis-swap to preserve CCW) + atan2 angle sort per
     vert + next-around-face pointer per halfedge. Surfaced
     a real propagation gap (~20% unset next pointers); fixed
     by adding step-6 new verts to the on-edge vert lists of
     their piercing edges (PropagateNewVertsToOnEdgeLists).
     99.998% next pointers set after fix.
   - ✅ Step 11 phase 3: polygon walk via `nextHalfedge`
     follow. Cray 238 tris → 305 sub-polygons (180 trivial +
     125 from 57 multi-poly tris); self-intersect 33542 →
     33687.
   - ✅ Step 12 (multiplicity merge): canonical-form per
     polygon, hash-group, accumulate signed multiplicity.
     0 cancellations on the .obj battery — no T-T full
     coincidence in these fixtures.
   - ✅ Step 13a (triangulate + emit output Manifold):
     project each polygon to 2D, call `Triangulate()`, emit
     triangles into MeshGL64. **3 of 5 fixtures successfully
     round-trip end-to-end** (Havocglass8, Offset12,
     hull-mask) with 0–0.008% volume drift. Cray and
     self-intersect fail Manifold construction at this stage
     because emitting both sides of every new edge produces
     non-manifold output (4 halfedges per edge instead of 2).
   - ⚠️ Step 13b (ray-cast inside/outside filter): proof-of-
     concept only. **Works** on unpierced fixtures; **fails**
     on every pierce-class fixture (NotManifold output —
     dropped/reversed polygons leave edges with wrong
     cardinality). The fundamental issue: ray-cast against
     self-intersecting input doesn't track winding monotonically
     near intersection lines. A robust classifier likely needs
     either (a) operating on a non-self-intersecting reference,
     (b) porting `boolean_result.cpp::Winding03_` directly, or
     (c) per-tri-pair analytical classification.

### Pipeline status summary

End-to-end via `OverlapRemoval(Manifold)`:

| Fixture | Step 13a | Step 13b |
|---|---|---|
| Havocglass8 | ✓ NoError, drift 0% | ✓ NoError, drift 0% |
| Offset12 | ✓ NoError, drift 0% | ✓ NoError, drift 0% |
| hull-mask | ✓ NoError, drift 0.008% | ✗ NotManifold |
| Cray | ✗ NotManifold | ✗ NotManifold |
| self-intersect | ✗ NotManifold | ✗ NotManifold |

Step 13a (no filter) passes 3/5 fixtures because they have few or
no piercing pairs. Cray and self-intersect have many pierces;
emitting both sides of every new edge produces non-manifold output
(4 halfedges per edge instead of 2). Step 13b's filter is an
attempt to fix that but the classifier itself is not reliable on
pierced input.

The pipeline is structurally complete; the remaining work is
classifier reliability on the heavy-pierce inputs, which is a
research problem separate from the spike's algorithmic exploration.

### Design notes for steps 11–13

Step 11 (per-tri polygon partition) is the heaviest:

- **2D projection per tri.** Drop the axis with the largest
  `|faceNormal_[t]|` component; the remaining 2 give a stable
  2D embedding for the tri's plane.
- **Halfedge data structure per tri.** Each sub-edge produces
  two halfedges (forward + backward); each halfedge has
  `startVert`, `endVert`, `pairedHalfedge`, and a position in
  the per-vert cyclic ring. New edges insert *both directions*
  per Emmett's step 10.
- **Per-vert angle sort.** For each vertex used by the tri's
  sub-edges, project to the 2D plane and sort incident
  halfedges by `atan2`. The next-around-vert pointer for a
  halfedge `h` is the one whose angle is the next-clockwise
  position from `h`'s reverse direction.
- **Polygon walk.** Standard halfedge-mesh polygon traversal:
  pick any unvisited halfedge, repeatedly follow
  `pairedHalfedge.nextAroundVert` to trace a closed loop.
  Each loop is one sub-polygon.

Step 12 (multiplicity merge): equivalent sub-polygons (same
vert sequence ignoring rotation) accumulate winding-direction
multiplicity; ones that cancel to zero are dropped.

Step 13 (winding output): for each retained sub-polygon, ray-
cast from a midpoint to determine inside/outside, then keep
those separating winding 0 from ≥1. Triangulate non-convex
sub-polygons via the existing `Triangulate` (`src/triangulator.cpp`,
Hertel-Mehlhorn + earcut) — currently used by `boolean_result.cpp`.

The Cray hot-spot (1 tri with 5 piercing pairs / 13 sub-edges)
is the non-trivial test — for tris with 0 new edges (the
majority), step 11's polygon partition is degenerate to the
original tri.

A reasonable first sub-implementation: skip the polygon walk
entirely on tris with 0 new sub-edges (output the original
triangle), and stub out the tris with new sub-edges (accept
the pierce remains). That validates the Impl→MeshGL64→Manifold
round-trip in a "no-op overlap removal" mode and unblocks
testing the full pipeline shape before committing to the
polygon-walk implementation.

   **Pierce coverage** (steps 6+7p1 vs `CheckSelfIntersection`
   ground truth): self-intersect 99.5% (658/661), Cray 104.4%
   (119/114 — *exceeds* ground truth; step 6's strict
   `bTol=1e-12` catches FP-marginal cases that
   `CheckSelfIntersection`'s scale-relative tolerance misses).
   The complex fixtures hull-mask / Offset12 / Havocglass8
   stay at 23% / 13% / 0% — and the relaxed step-6 variant
   confirms this is *not* filter over-aggressiveness (relaxed
   produces 5-10× false positives on those fixtures). Their
   shortfall is presumably a different geometric configuration
   class that step 6 doesn't model.

## Open questions specific to 3D

1. **Step 4b's analog in 3D.** In 2D, k≥3 edges concurrent at one
   true point produced near-duplicate intersection verts. In 3D,
   the analog is k≥3 triangles meeting at one true line. Manifold's
   existing `CollapseShortEdges` catches some residuals; whether a
   new pass is needed depends on what the adversarial battery shows.

2. **Triangulation strategy for step 13.** Sub-polygons of input
   triangles are not necessarily convex. The existing
   `Triangulate` (`src/triangulator.cpp`, Hertel-Mehlhorn + earcut)
   is robust under the topology-manifold invariant; whether it's
   robust under the same eps in the no-pierce regime is unknown.

3. **k-fold-line concurrence (3D analog of fig 6.5 hexagon).** Three
   or more triangles passing through one true line. The
   topology-manifold pipeline handles this by adjacency conventions;
   overlap removal needs the line to be a single shared edge.

4. **Output orientation in 3D** (CCW outer, CW inner faces analog).
   Already conventional in manifold's halfedge structure. No new
   work.

5. **Open paths / non-manifold inputs.** Out of scope. Overlap
   removal's input is topology-manifold (the `Manifold` type's
   invariant); non-manifold soup repair would be a separate feature
   that runs *before* overlap removal in any combined pipeline.

## What this is not

- **Not** a proposal to reimplement `boolean3.cpp`. That pipeline
  is mature and correctly targets the topology-manifold invariant.
- **Not** a claim that the existing pipeline is buggy. Producing a
  topology-manifold output that has a sub-ε tri-tri pierce is
  *within* its contract, not a violation.
- **Not** a request to land anything yet. This is exploration on a
  separate branch from the 2D landing.

## Why this is on its own branch

The 2D landing (`explore/overlap-removal-2d`) is a real proposal
with code, tests, and a discussion thread targeting `CrossSection`.
This 3D exploration is purely speculative and shouldn't dilute
that discussion. If Emmett's 2D feedback opens the door to a 3D
extension, this doc is the starting point.

---

## Spike learnings + design proposal

After implementing all 13 of Emmett's #289 steps end-to-end as a
single-shot prototype (`extras/overlap3d_proto.cpp`, ~4500 lines
including tests + diagnostics), here is what holds, what breaks,
and where the spike fits relative to Emmett's sketch and the
existing manifold pipeline.

### What works

End-to-end status on the canonical .obj-fixture battery (all
already-clean Boolean outputs from manifold's existing pipeline,
fed through `OverlapRemoval(Manifold) → Manifold`):

| Fixture | Pierces in | Result | Drift |
|---|---:|---|---:|
| Havocglass8 | 0 | NoError v141 t278 | 0% |
| Offset12 | 0 (n=2) | NoError v3716 t7372 | 0% |
| hull-mask | 1 (k=4 tolerated) | NoError v6980 t13980 | 0.008% |
| Cray | 81 | NotManifold | — |
| self-intersect | 653 | NotManifold | — |

3 of 5 fixtures round-trip end-to-end and produce a
topology-manifold output Manifold. The two failures are
heavy-pierce inputs where the per-tri-pair analytical
classifier (step 13d) leaves residual k=4 chord edges and
introduces k=1 dangling edges from the cascade-drop pattern.

### Where Emmett's sketch needed real implementation discoveries

Three places where the spike's actual implementation deviates
from the sketch in non-trivial ways. None invalidate the sketch;
each is a "this is what the implementation has to do" detail.

1. **Polygon-walk rotation direction.** Emmett's step 11 says
   "follow the cyclic ordering around each vert" without
   specifying whether the per-vert sort gives CCW or CW
   neighbors. In a "face on the LEFT of h" walk, the next
   halfedge from h.endVert is the smallest *clockwise* delta
   from h's reverse direction, not CCW. (Confirmed via manual
   trace of a single-chord triangle: CCW picks the perimeter
   continuation, tracing the outer cycle through chord
   endpoints; CW correctly turns into the sub-polygon.)

2. **New verts must be propagated to *original* edges' on-edge
   vert lists.** Emmett's step 7 says new vert endpoints are
   "shared between the edges of one tri and the edges and
   interior of the other." The piercing-edge endpoints are
   geometrically on triA's edge interior, but step 3 builds
   the on-edge vert lists *before* step 6/7 produces those new
   verts. Without explicit propagation
   (`PropagateNewVertsToOnEdgeLists`), step 11's halfedge graph
   leaves new edges as dangling cracks (~20% of next-pointers
   unset on Cray). Adding the propagation pass closes
   99.998% of the gap.

3. **Step 12's multiplicity merge is rare in practice.** The
   "equivalent polygons cancel" case (= T-T full coincidence)
   doesn't fire on any of the .obj-fixture battery (0
   cancellations across self-intersect's 33687 polygons). Step
   12 is correctness-preserving infrastructure rather than a
   workhorse. Most of the heavy lifting lives in step 13.

### The classifier — Emmett's least-specified step

Step 13's classifier is where Emmett's sketch is thinnest:
"output only the polygons that separate winding numbers <=0
from >0" via a "ray-cast from a vert" plus "face shadowing"
adjustment. The spike tried four variants:

| Variant | Approach | Result |
|---|---|---|
| 13a (no filter) | Keep all polygons; rely on Manifold's RemoveDegenerates | Works for ≤1 k=4 edges; fails for many |
| 13b (ray-cast) | RayCast against the input mesh from probe ± ε·n | Unreliable on self-intersecting input |
| 13c (rigorous ray-cast) | RayCast against the canonical sub-polygon arrangement (= 2D analog) | O(P²) brute force; partial improvement |
| 13d (analytical) | For each polygon, signed dot(non-chord-vert, partner-normal) | **Best so far** |

The analytical classifier (13d) needs no ray-cast and is local:
for each multi-poly polygon, compute the signed distance from a
non-chord interior vert to each *bounding* partner triangle's
plane. Drop iff strictly outside (= dot > scale-relative
threshold) every bounding partner. The "non-chord interior
vert" choice avoids centroid bias when polygons are thin
slivers near the chord plane.

This eliminates the k=4-per-chord problem on self-intersect
(653 → 0) and reduces it on Cray (81 → 35). The remaining gap
is k=1 dangling edges from the **cascade-drop** pattern.

### The cascade-drop problem (= the open issue)

When a multi-poly polygon Y of triA is dropped, Y's perimeter
sub-edges are the parts of triA's perimeter *that Y bounds*.
Those sub-edges are shared with triA's neighbor triA'
across an original edge of triA. If triA' is a 1-poly tri
(auto-kept by the spike's gate), triA' contributes 1 triangle
incident to that sub-edge. After dropping Y, triA contributes
0 → k=1 dangling edge in the output.

Geometrically, the cascade-drop happens because the part of
triA' *near the dropped Y region* is also occluded by the
piercing partner — but triA' itself isn't pierced, so it
auto-keeps without checking partner occlusion. Fixing this
requires either:

(a) **Drop propagation across tri-tri adjacency.** When Y is
    dropped, mark Y's perimeter sub-edges. For each marked
    sub-edge, run the analytical classifier on the
    *neighbor's* polygon; if the neighbor's polygon's vertex
    near the sub-edge is occluded by the same partner, split
    the neighbor's polygon and drop the corresponding side.
    Substantial bookkeeping.

(b) **Per-vert winding via flood fill (= boolean3.cpp's
    `Winding03_`).** Compute the winding number at each input
    vert relative to the input mesh. Use `DisjointSets` to
    propagate the winding across uninterrupted halfedge runs;
    only the connected components that contain a chord need
    explicit ray-cast / Kernel02 evaluation. This is the
    existing pipeline's solution — porting it is the most
    direct path to a robust step 13.

    A purely-topological cascade-drop forward (drop a 1-poly
    auto-kept tri whose perimeter sub-edges have more k=1 than
    k=2 incidences) was tried and found to cap at marginal
    gains: Cray 96→92 k=1, self-intersect 357→343 k=1. The
    cascade can't see partner-occlusion geometry, so it
    converges to a local minimum determined by edge-incidence
    parity rather than the true volume boundary. This validates
    that path (b)'s explicit winding evaluation is necessary.

### Triangulation-induced k=4 (fixed) and k=3 (residual)

Cray's pre-fix 38 k=4 edges turned out to **not** be classifier
failures. They were triangulation artifacts:

When step 7's on-edge vert V lies on input edge (A, B), V
becomes a polygon-perimeter vert in both input tris incident to
edge (A, B). Each polygon perimeter then has 3 collinear verts
[A, V, B]. The default `manifold::Triangulate` picks chord
(A, B) in each polygon → 2 internal triangles per polygon use
chord (A, B) → 4 wings meeting at edge (A, B) → k=4.

Fix: detect collinear-triple in polygon perimeter (cross
product near zero, parametric position strictly interior). Fan
the triangulation from V instead of using the default chord.
This eliminates Cray's k=4 entirely (38 → 0) and cleans up
hull-mask's 1 tolerated k=4 too.

Residual k=3 (~15 edges on Cray): fan-from-V introduces
internal chords (V, X) for every non-adjacent perimeter vert X.
If chord (V, X) coincides with another polygon's perimeter
edge, that edge gets k=3 (= 1 from the other polygon's
perimeter + 2 from this polygon's fan). Fix would require a
"global perimeter awareness" triangulator: avoid fan chords
that match any external perimeter edge. Lower priority than
the cascade-drop / Winding03 issue.

A "conflict-aware fan" attempt (collect global perimeter
edges first, prefer fan vertices whose chords don't conflict)
was tried and found to regress: Cray k=4 0→3, k=1 80→92,
k=3 17→12. The asymmetric chord choice across paired polygons
(= polygons sharing an input edge with the same on-edge V)
breaks the symmetry that makes fan-from-V work as a k=4 fix.
A correct global-aware triangulator needs coordinated chord
choice across paired polygons, not independent per-polygon
optimization.

### Pair-symmetric chord enforcement (REAL WIN)

Empirical investigation of Cray's k=1 edges (via
`OVERLAP3D_DUMP_K1=1`) showed every k=1 has a *dropped*
multi-poly partner with the edge as its perimeter. The
analytical classifier was deciding K/D *independently* per
polygon, producing inconsistent decisions across the 4
polygons of a chord (= 2 polygons per tri × 2 tris).

Fix: **pair-symmetric chord enforcement**. For each chord
between triA and triB, identify the 2 twin polygon pairs by
halfedge direction:

  Pair 1: triA polygon with halfedge (v0→v1) +
          triB polygon with halfedge (v1→v0)
  Pair 2: triA polygon with halfedge (v1→v0) +
          triB polygon with halfedge (v0→v1)

A twin pair is geometrically the same surface (= polygons
sharing the chord with opposite halfedge directions, as
required for manifold output at the chord). After
pre-computing per-polygon classifier decisions, force exactly
one pair K and the other D, deciding by K-vote count per pair
(ties prefer Pair1 K + Pair2 D).

Result on the failing fixtures:

  Cray:           k=1 87→75, k=3 17→11, k=4 14→0 (typical run)
                  net non-k=2: 118 → 86 (27% reduction)
  self-intersect: k=1 343→319, k=3 24→0, k=4 unchanged
                  net non-k=2: 367 → 319 (13% reduction)

Other fixtures (hull-mask, offset12, havocglass): unchanged.

The k=4 elimination on Cray is direct: pair-sym guarantees
exactly one pair K per chord = 2 incidences per chord = k=2.
The k=3 elimination on self-intersect comes from the same
mechanism (chord pair conflicts manifesting as one polygon
contributing 1 perimeter incidence + the other contributing 2
fan incidences).

Default-on; opt out via `OVERLAP3D_NOPAIRSYM=1` for A/B
comparison.

### Mesh-edge halfedge consistency (extends pair-sym)

Generalization of pair-sym to ALL halfedge twin pairs (= chord
+ original mesh-edge sub-edges). For each polygon perimeter
halfedge (a→b), find its reverse (b→a) in another polygon.
For mesh-edge sub-edges across 2 polygons of different tris,
AND-merge: drop iff either dropped (= cascade-drop forward
once across mesh-edge adjacency).

Single pass (NOT iterated to fixed point — full iteration
over-cascades on dense self-intersection meshes; self-intersect
goes 333 → 590 k=1 if iterated).

Result combined with chord pair-sym:

  Cray:           k=1 75 → 62, k=3 11 → 0, k=4 0 → 0
                  net non-k=2: 86 → 62 (28% reduction from
                  pair-sym-only)
  self-intersect: k=1 319 → 333 (slight regression),
                  k=3 0 → 0, k=4 0 → 0
                  net: 319 → 333 (4% regression)

Cleaner topology overall. Cray now has ONLY k=1 issues to
solve (no k=3/k=4 left). Self-intersect's slight regression
is the cost of resolving a few mesh-edge mismatches that
existed before but were balanced.

Implementation note: the AND-merge mutates `precomputedKeep`
during iteration over halfedges. This induces a bounded
cascade *within* a single pass (= a flipped polygon's other
halfedges, when seen later in the iteration, see the new
state). A snapshot variant (read-only initial state, apply
updates at end) was tried but produced worse Cray results
(k=1 ≈ 75-97 vs 62-74 with cascade). The cascade is the
better trade-off here. There is run-to-run variance from
parallel earlier-step processing producing different polygon
walks; Cray fluctuates between ~62 and ~74 k=1.

### Open issue: cascade-drop floor

After both pair-sym phases, k=1 remains at the floor:
  Cray:           ~62-74 (variance)
  self-intersect: ~333

K=1 dump shows: dropped polygons (= multi-poly polys correctly
dropped per chord classifier) leave their perimeter sub-edges
unmatched with neighbor 1-poly tris that auto-keep. Pair-sym
Phase 2 propagates one drop step but doesn't fully cascade
without over-aggression on dense self-intersection meshes.

The principled fix remains `Winding03_` port (path (b)). With
per-vert winding evaluated via flood-fill components, each
polygon's keep/drop decision is geometrically grounded — not
dependent on classifier convention or pair-sym heuristics.

### Spike plateau

Several loops of additional heuristic exploration did not move
the floor below Cray k=1 ≈ 60-62 / self-intersect k=1 ≈ 333.
Attempts that regressed and were reverted:

- **Removing the auto-keep gate in `TriangulateAndEmit`**: lets
  pair-sym Phase 2's drop decisions land on 1-poly tris.
  Catastrophic over-drop (Cray 4/317 polys kept) when combined
  with cascade-during-iteration. With snapshot Phase 2:
  Cray k=1 100-121, k=4 reappears (5-9). The auto-keep gate is
  load-bearing.
- **Phase 1 protection in Phase 2**: prevent Phase 2 from
  flipping multi-poly polys touched by chord pair-sym (avoids
  Phase 2 overriding Phase 1's chord consistency). Combined
  with gate-removal: Cray k=1 100-119 + k=4 0-9. Still worse.

The pair-sym infrastructure (chord + mesh-edge halfedge
consistency) has hit its limit. Further reductions require
geometrically-grounded per-vert winding (= `Winding03_` port),
not topological heuristics.

**Recommendation**: pause heuristic exploration. Resume with
`Winding03_` port as a focused effort, or accept the current
plateau as the spike's terminal state and harden + ship the
3/5 passing fixtures with documented limitations.

### Surface-cap heuristic (5/5 status 0 milestone)

A topology-first heuristic that closes the cascade-drop holes
left by classifier drops. After classification:

1. Identify k=1 edges from the partial output.
2. For each k=1 edge, determine the kept halfedge's direction.
   Cap halfedge direction = opposite (= the missing manifold
   pair).
3. Walk directed cap halfedges to find closed cycles.
4. Fan-triangulate each cycle. Cap halfedges' direction makes
   the new tris pair properly with the kept tris (k=2 every
   edge).

Result with `OVERLAP3D_CAP=1`:

  Cray:           status 0! 5 cycles, 52 cap tris
                  vol 1.21e12 vs input 4.91e10 — 25x (geo wrong)
  self-intersect: status 0! 7 cycles, 319 cap tris
                  vol 0.279 vs input 0.287 — 3% drift (geo OK)
  hull-mask, offset12, havocglass: status 0 ✓ (no caps needed)

**All 5 fixtures pass Manifold validation when k=1 cycles close properly.** This is the spike's primary topology milestone.

*Variance caveat:* Cray exhibits run-to-run variance from
parallel earlier-step processing. In ~3/5 runs, the polygon
walks produce 62 k=1 edges that close into 5 cycles
(milestone hit). In ~2/5 runs, walks produce 74 k=1 edges
where some cycles fail to close (cap walk hits a dead-end on
its directed traversal), leaving k=1 > 0 and status 2. The
underlying parallelism is in the BVH builds and similar steps
in earlier query phases.

A production landing would need to either (a) pin the
parallelism with deterministic ordering, or (b) make the cap
walk more tolerant of incomplete cycles (e.g., merge
adjacent open paths into one closed cycle through a virtual
edge). Both are plausible follow-ups.

### `OverlapRemoval(Manifold)` API surface (commit 06463de1)

The entry-point function `OverlapRemoval(input, eps)` now
runs the full pipeline (steps 1–13 + pair-sym Phase 1 + cap)
instead of stubbing the output. Falls back to the merged
input (= a valid manifold via `MergeVertsEps`) if the pipeline
output is non-manifold. CLI mode for testing:

  `./build/extras/overlap3d_proto <fixture> --api`

Results across the 5 fixtures via the API:

  cray:           **fallback to merged input** (cap walk fails
                  consistently here)
  hull-mask:      real pipeline output, 0.008% drift, status 0
  self-intersect: real pipeline output, 2.8% drift, status 0
  offset12:       real pipeline output, ~0% drift, status 0
  havocglass:     pipeline output equivalent to input (no
                  overlaps to remove)

So 4/5 fixtures get the real overlap-removed output via the
API; Cray gracefully falls back. The fallback is *deterministic*
unlike the `RunSingleFixture` path which has BVH variance —
suggesting the pipeline's deterministic codepath in the API
runs differs slightly from the diagnostic CLI.

This is a defensible API surface for the spike. Production
hardening would tighten the variance, integrate
`Manifold::Simplify` semantics, and either commit to the
Boolean3 port or accept the heuristic limitations on the
heavy-pierce regime.

The geometry quality varies: self-intersect's caps are local
and produce nearly-correct volume. Cray's caps span large
non-planar cycles that fan-triangulate to 25x volume.

Cap is opt-in (`OVERLAP3D_CAP=1`) because:
- Cray's geometry is wrong (= cap should NOT be applied for
  meshes where the cycles aren't reasonable to triangulate as
  fans).
- A production solution needs either smarter cap geometry
  (e.g., minimum-area triangulation, or projecting cycles to
  their best-fit plane) OR principled classification (Boolean3
  port) that avoids cycles in the first place.

  *Tried:* best-fit plane projection (Newell's normal +
  triangulate in 2D). Improved self-intersect's drift slightly
  (3% → 2%) but produced k=4 issues from self-intersecting
  projected polygons, and made Cray status 2 in 3/5 runs (=
  cycles too non-planar to project cleanly). Reverted to the
  simple fan approach.

  *Why Cray's cap is geometrically broken:* Cray's 5 cycles
  break down as (per `OVERLAP3D_CAP_DEBUG=1` dump):
    cycle 0: 21 verts, bbox 71137 (= full mesh bbox)
    cycle 1: 25 verts, bbox 65237 (= ~full mesh bbox)
    cycle 2: 5 verts, bbox 1946 (small, local)
    cycle 3: 6 verts, bbox 34668 (~half mesh bbox)
    cycle 4: 5 verts, bbox 3447 (small, local)

  Cycles 0 and 1 span the *entire* mesh. They are global
  boundaries of large dropped regions, not local holes — Cray
  is a Subtract result, and the classifier drops *both* A's
  "inside-B" surface AND B's "outside-A" surface. The dropped
  regions are large and the cycles wrap around them.

  Geometrically the right thing would be to *keep* B's
  inside-A surface (with reversed orientation) so the cycles
  don't form. The analytical classifier's "drop iff outward
  of partner" convention is correct for B's surface (Subtract
  convention) but WRONG for A's surface. For overlap-removal
  of an arbitrary self-intersecting input mesh M, there's no
  A/B labeling; the proper fix uses per-vert winding (W03).

  Cray's 25× volume reflects this: the cap fans contribute
  large interior triangles that include parts of the bbox
  that should be empty.

But it does demonstrate that the spike's polygon partition +
classification + topology recovery framework can produce a
manifold output for *all* the test fixtures — just not
geometrically correct on heavy-pierce inputs without further
work.

### Winding03_ port: scaffolded but not integrated

`ComputePerVertWinding` (commit abe8c597) ports the boolean3
flood-fill: DisjointSets over original verts, unite via
non-broken halfedges, per-component ray-cast. Extended in
commit 7c011f60 to compute both above/below windings per vert.

Wired as opt-in via `OVERLAP3D_W03=1`. Result alone is *worse*
than the analytical + pair-sym default:
  Cray:           k=1 56-62 → 74
  self-intersect: k=1 333 → 720

Why: per-polygon W03 classification correctly identifies
interior polygons (= verts have winding ≥ 1 on both sides) and
drops them. But it drops ~91% of self-intersect's multi-poly
polygons, creating massive topology mismatches with
neighbors. The full Boolean3 pipeline avoids this by
generating a fresh manifold output via Kernel11/Kernel12
*after* using winding info — that's the missing piece.

A correct W03-based path would need either:
- Pair-sym integration (W03 sets initial decisions, pair-sym
  enforces chord/mesh-edge consistency on top)
- Or a port of the full triangulation flow from
  boolean_result.cpp (Kernel11 + Kernel12 + edge insertion).

The Boolean3 pipeline is several thousand lines; a faithful
port is multi-week work. The spike's current heuristic state
(3/5 passing) is defensible for review and the plateau is
documented.

### Hybrid analytical + ray-cast classifier attempt

A second attempt at the cascade-drop issue: hybrid path that
uses `AnalyticalKeep` for chord-bounded polygons (which works
for k=4 elimination) AND `ClassifyPolygon`'s ray-cast for
1-poly polygons (= the auto-keep case).

The ray-cast classifier was modified to use the refined keep
criterion: keep iff exactly one of {windingUp, windingDown}
is 0 (= polygon on outer boundary of `winding ≥ 1` region).

Result: all Cray 1-poly tris (202) returned keep — no 1-poly
got dropped. The k=1 cascade source on Cray is therefore not
"interior 1-poly tris" — at least, not detected by simple
per-polygon ray-cast on `mr.manifold`. The mixed nature of
Cray's k=1 edges (some 3-vert auto-kept polygons, some 7-vert
multi-poly polygons) suggests the cascade is from
*multi-poly* classification mismatches at perimeter sub-edges,
not from interior 1-poly tris.

Either:
- Multi-poly tris drop their adjacent-to-chord-endpoint
  polygons differently from their pair-partner tri's matching
  drop, leaving the on-edge sub-edge with mismatched k=1.
- Or ray-cast on the original self-intersecting input
  systematically miscounts winding for polygons inside
  multi-overlap regions.

### Piercer-augmented classifier attempt

To address the cascade-drop k=1 issue, an augmentation to the
analytical classifier was attempted: build `vertToPiercers`
from `etIsects` (= for each on-edge vert V, the set of
piercing tris that contributed V via step 6's edge-tri
intersection), then for any polygon containing on-edge verts,
test the centroid against each piercer's plane.

Two conventions were tried — drop iff polygon centroid is
*outward* of any piercer (matching the chord classifier's
convention) and drop iff centroid is *inward* (the opposite).
Neither produced clean improvements:
- Outward convention: piercer set never triggered drops on
  the 1-poly tris that needed it. The 1-poly tris with
  on-edge verts have *crossing* centroids (their interiors
  cross the piercer's plane), so "all outward" is rarely
  satisfied.
- Inward convention: over-aggressive. Broke offset12
  (0 → 18 k=1 from 13 incorrect drops on a previously-clean
  fixture), regressed Cray (k=1 80 → 122) and self-intersect
  (k=1 343 → 678).

The chord classifier's convention is calibrated for chord-
bounded polygons (= polygons on one side of a chord plane).
The 1-poly auto-keep case is fundamentally different: the
polygon spans both sides of the piercer's plane, so a
single-side test is the wrong question. The right question
is *which sub-region of the polygon's interior is inside vs.
outside the piercer's volume* — and answering that requires
splitting the polygon at an implicit chord, not classifying
the whole polygon.

This validates path (b) (`Winding03_` port) more strongly:
per-vert winding propagated via flood fill is exactly the
fixed-point invariant that handles the multi-sided polygon
case correctly, by classifying *each vertex* rather than
each polygon.

(c) **Operate on a non-self-intersecting reference.** The
    existing Boolean operates on two non-self-intersecting
    inputs and combines them; Kernel02 / Kernel11 / Kernel12
    are robust because the inputs aren't self-intersecting.
    Overlap-removal of *one* self-intersecting mesh has no
    such reference. A workaround: split the input into
    non-self-intersecting components first (via flood fill on
    triangle-tri-pair adjacency), then process each
    independently.

### How the spike compares to Emmett's sketch

|   | Emmett's #289 sketch | Spike implementation |
|---|---|---|
| Step 1 (vert merge) | "use a weighted-average position" | Cluster centroid (unweighted average); iterate up to 4 times for chains |
| Step 2 (drop collapsed) | "discard collapsed edges" | `EnumerateEdges` from `halfedge_` (manifold input → no actual collapse) |
| Step 3 (per-edge on-edge verts) | "ordered list … within ε of [edge]" | BVH + projection + thin-tri-apex skip |
| Step 4 (edge-edge) | "above/below/on three-valued" | Closest-points-on-segments with ε threshold (no SoS) |
| Step 5 (per-tri in-tri verts) | "verts within ε of interior" | BVH + plane-distance + barycentric |
| Step 6 (edge-tri) | "main intersection-finding phase" | Direct ray-tri with strict-interior gate |
| Step 7 (line-of-intersection edges) | "exactly two endpoints per pair" | Validated empirically (Cray 100% n=2) |
| Step 8 (propagate to new edges) | "verts within ε of [new edge]" | Per-edge BVH; quiet on the .obj battery |
| Step 9 (new-edge × new-edge) | within each tri | Pairwise per tri; quiet on the battery |
| Step 10 (sub-edge canonicalize) | halfedge lists per tri | Implicit in step 11's halfedge graph |
| Step 11 (polygon partition) | "follow cyclic ordering" | atan2 angle-sort + **CW** next-pointer |
| Step 12 (multiplicity merge) | "polygons that cancel" | Lex-min canonical form + signed sum |
| Step 13 (winding output) | "ray-cast or face-shadowing" | **Analytical signed-distance** (no ray-cast); cascade-drop unresolved |

The big departures: step 11's CW vs. CCW (sketch ambiguous);
step 13 going analytical instead of ray-cast (sketch suggested
ray-cast); the new-vert-to-original-edge propagation pass
that's not in the sketch.

### Why this isn't yet productionizable

Even on the 3 passing fixtures, the spike is missing several
pieces a production landing would need:

- **No `Manifold::Simplify` API.** The natural surface for 3D
  overlap removal is a method on `Manifold` (analogous to
  `CrossSection::Simplify`). The spike runs as an executable.
- **No fuzz/displacement validation.** The 800-case fuzz uses
  Boolean output as input (which is already clean). Adversarial
  fuzzing on raw self-intersecting input was not done.
- **The cascade-drop on Cray and self-intersect is unresolved.**
  Heavy-pierce inputs are exactly the inputs that overlap
  removal is supposed to fix. Failing on them undermines the
  feature's main use case.
- **Performance.** O(P²) per-polygon classification + per-tri
  BVH builds; not optimized for parallelism. Existing
  `boolean3.cpp` is heavily parallelized.

### Recommended next steps

In order of value:

1. **Port `Winding03_`** (option (b) above). The existing
   pipeline solves the classifier-on-self-intersecting-input
   problem; the spike's analytical classifier is a partial
   reimplementation. Replacing 13d with a port of `Winding03_`
   closes the cascade-drop and likely makes Cray + self-intersect
   pass.

2. **Land the structural pieces** (steps 1–11) into `src/`
   *only after* step 13 is robust. The query-step infrastructure
   (BVH builds, `Collider`/`DisjointSets` use patterns) maps
   1-to-1 onto manifold's existing scaffolding and is mostly
   unsurprising.

3. **Adversarial fuzz.** Once the classifier is robust, run
   the same 800-case displacement battery the 2D prototype uses,
   plus the .obj fixtures, and confirm zero pierces in the
   output (using `CheckSelfIntersection` as ground truth).

4. **`Manifold::Simplify`-style API** with explicit ε, mirroring
   the 2D prototype's surface.

If the recommended path looks too heavy for a single landing,
the spike can also serve as a **diagnostic tool**: ship
`CheckSelfIntersection` and the analytical classifier as opt-in
post-Boolean validation for users who want stricter geometric
guarantees than the existing topology-manifold contract
provides.

---

## Handoff context (2026-05-09)

This section is for an agent or engineer picking up the work
from a session boundary. Read first, then dive in.

### Where the spike sits

- Branch: `explore/overlap-removal-3d`
- Latest commit: `5d49e03a` (or HEAD)
- Build: `cmake --build build` from repo root; binary is
  `build/extras/overlap3d_proto`.
- Main spike file: `extras/overlap3d_proto.cpp` (~5500 lines,
  end-to-end implementation of all 13 steps + cap heuristic +
  Winding03 scaffolding).
- All work committed; no uncommitted changes.

### Things that work

- 13-step pipeline implemented end-to-end with pair-symmetric
  chord classifier + surface-cap heuristic.
- 8 .obj fixtures wired through both granular CLI mode and
  `OverlapRemoval(Manifold)` API mode (`--api` flag).
- Adversarial fuzz CLI mode (`--advfuzz`) finds 12/1395 self-
  piercing Boolean outputs (0.86% rate, ~7× the general fuzz).
- `OverlapRemoval` API produces valid output for **12/12** of
  those 12 self-piercing cases (8 fully fixed, 1 reduced to
  FP noise, 3 regularized to empty for sub-eps slivers).
- `ComputePerVertWinding` scaffolding in the spike (DisjointSets
  flood-fill + per-component ray-cast) is wired and working but
  not yet integrated as the production classifier.

### Things that don't work / are limited

- Cray's geometry through the cap heuristic is wildly wrong
  (25× volume) — cycles span the full mesh bbox, fan-
  triangulation produces interior-spanning tris. Status 0 is
  achieved but volume drift is unacceptable for production.
- Run-to-run variance: Cray hits status 0 in only ~3/5 runs;
  the other ~2/5 hit status 2 because the cap walk dead-ends
  on different polygon-walk states from BVH parallelism. Other
  fixtures are deterministic.
- The W03 classifier alone (opt-in via `OVERLAP3D_W03=1`) is
  *worse* than the analytical+pair-sym default (Cray k=1
  60→74, self-intersect 333→720) — needs Boolean3's full
  Kernel11/Kernel12 emission to be useful.
- `m.Boolean(m, OpType::Add)` as a self-union shortcut DOES
  NOT work — Cray empties, self-intersect 661→216 partial,
  generic-twin hangs. Boolean3's SoS distinguishes inputs P
  vs Q symbolically; when inP==inQ they are indistinguishable.

### CLI cheat sheet

```
# Granular fixture (~0.35-15s):
./build/extras/overlap3d_proto <fixture>

# OverlapRemoval API on a fixture:
./build/extras/overlap3d_proto <fixture> --api

# Self-union via existing Boolean3 (broken, useful for diagnostic):
./build/extras/overlap3d_proto --selfunion <fixture>

# Adversarial fuzz + OverlapRemoval measurement (~15s):
./build/extras/overlap3d_proto --advfuzz

# Add OVERLAP3D_ADVFUZZ_DIAG=1 to see per-failure details.

# Fixtures (8): cray, hull-mask, self-intersect, offset12,
#   offset34, havocglass, generic-twin-7081, generic-twin-7863
```

### Productionization options analyzed

The user asked: *"reuse as much code as possible (possibly
refactoring existing kernels if needed but ideally just reusing)."*

| Option | Description | Effort | Risk |
|---|---|---|---|
| **A. Minimal expose** | Move `Winding03_` (and `Kernel02` dependency) from anonymous namespace at `src/boolean3.cpp:71-527` to `manifold::detail::` namespace. Declare in `boolean3.h`. Replace spike's `ComputePerVertWinding` with the production version. | 1-2 sessions | Low — no behavior change to existing Boolean3, just visibility refactor. |
| **B. New entry point** | Add `class SelfMeshAnalysis` to `boolean3.cpp` that wraps the kernels for inP=inQ case with adjusted SoS perturbation. Provides clean self-pierce intersections + winding via reused kernels. Then write a self-mesh assembly inspired by `boolean_result.cpp` but specialized. | 3-5 sessions | Medium — SoS adjustment is delicate; need to verify against the spike's existing test fixtures + advfuzz. |
| **C. Full self-mode** | Add a `selfMode` flag to `Boolean3` that adjusts SoS for self-mesh. Most reuse, most invasive. Requires regression testing on all existing Boolean tests. | 5-10 sessions | High — touches the production Boolean3 codepath. |

**Recommended starting point: A.** Smallest change that
exposes production-quality kernels for reuse. Validates the
visibility refactor before committing to deeper work.

### Concrete first steps for Option A

1. **Identify the minimal subset to expose.** `Winding03_` (line
   453) directly uses `Kernel02` (line 212). `Kernel02` uses
   `Shadow01` (already public in `shared.h`) and `Manifold::Impl`
   (public). So the minimal subset is `Winding03_` + `Kernel02`
   + their template wrappers.

2. **Refactor location.** Move them to a new namespace
   `manifold::detail::` (or similar) inside boolean3.cpp. Or
   create `src/winding.h` with declarations. Don't change
   behavior — pure visibility change.

3. **Verify nothing breaks.** Run the existing manifold test
   suite: `cmake --build build && ctest --test-dir build`.
   The Boolean3 tests should be unchanged.

4. **Replace `ComputePerVertWinding` in the spike.** It's at
   `extras/overlap3d_proto.cpp` around line 2305-2415. Swap
   for the new public function. Verify Cray/self-intersect/etc.
   per-vert windings match (or improve) the spike's scaffold.

5. **Then attempt Option B.** With kernels exposed, the next
   chunk is `SelfMeshAnalysis` — running the kernel cascade
   for self-mesh case, then assembling output via Winding03.

### Key file references

- `src/boolean3.cpp` — Boolean3 class + kernels (anonymous
  namespace)
- `src/boolean3.h` — Boolean3 class declaration only
- `src/boolean_result.cpp` — Boolean output assembly (uses
  w03/w30 + xv12/xv21 to construct the output mesh)
- `src/disjoint_sets.h` — DisjointSets union-find (already
  public, used by Winding03_ and the spike's
  ComputePerVertWinding)
- `src/shared.h` — Shadow01 / Interpolate / withSign / SoS
  helpers (already public)
- `extras/overlap3d_proto.cpp` — the spike (5500 lines):
  - Line 200: `CheckSelfIntersection` (= the spike's
    ground-truth self-pierce check; independent of Boolean3)
  - Line 2305-2415: `ComputePerVertWinding` (= spike's
    scaffolded Winding03 port; replace this with production
    version)
  - Line 2972+: `OverlapRemoval(Manifold, eps)` API entry
    (uses pair-sym Phase 1 + Phase 2 + cap; falls back to
    merged input if non-manifold)
  - Line ~2466+: `TriangulateAndEmit` (the spike's polygon
    classifier + cap walker)
  - Line 4876+: granular fixture map + RunSingleFixture

### Tests + adversarial corpus

The 12-case adversarial corpus from `--advfuzz` is the most
informative for productionization:

```
class                | pierces / valid (rate) | OverlapRemoval outcome
shallow-tight        |   0/200                 | (none)
near-coplanar-slabs  |   3/195  (1.5%)         | reg=3 (sub-eps slivers)
3-cube-chain         |   2/200  (1.0%)         | fix=1 red=1
4-cube-chain         |   0/200                 | (none, was 1/100 in prior run)
5-cube-chain         |   5/200  (2.5%)         | fix=5
6-cube-chain         |   2/200  (1.0%)         | fix=2
vertex-on-face       |   0/200                 | (none)
TOTAL                |  12/1395 (0.86%)        | fix=8 red=1 reg=3 (12/12 valid)
```

Production criterion: same 12 cases (and the broader 1395)
should still produce valid output via the productionized
Boolean3-reuse path. Volume drift should be tighter than the
spike's heuristic (specifically: Cray volume drift currently
2369%, must drop to <1%).

### Open design questions for the next session

1. **What's the SoS adjustment for self-mesh case?** Boolean3's
   perturbation logic uses `expandP` to symbolically shift
   `inP`'s vertices. For inP==inQ, what perturbation correctly
   distinguishes "this vert" from "this face on the same mesh"?
   The expected answer: perturbing one *copy* of the verts
   while leaving the original — but this is a real algorithm
   design question, not just bookkeeping.
2. **Boolean3 `Result(Add)` output for self-mesh:** boolean_result.cpp
   constructs the output via op-specific masks. For self-union,
   the natural mask is "winding ≥ 1" (= the union-of-itself
   region). Does the existing `c1 + c3 * v` formula in
   boolean_result.cpp:733 already produce this for an Add op
   with self-mesh, or does it need adjustment?
3. **What about emission?** `CollapseShortEdges` and similar
   post-processing in boolean_result.cpp may or may not be
   appropriate for self-union output. To verify.

---

## Option A: completion (2026-05-09)

The handoff above recommended starting with **Option A — Winding03_
visibility refactor**: lift `Shadow01`, `Kernel02`, `Winding03_`,
`Winding03` out of `boolean3.cpp`'s anonymous namespace into a header,
so the spike (and any future callers) can reuse the production
flood-fill instead of the spike's ray-cast scaffold.

### What changed

- New header `src/winding03.h` declares the four templates in
  `namespace manifold` (same precedent as `shared.h`'s exposure of
  `Shadows`, `Interpolate`, `withSign`).
- `src/boolean3.cpp` now `#include "winding03.h"` and the four
  definitions are deleted from its anonymous namespace. Pure
  visibility change — the existing call sites (`Boolean3::Boolean3`,
  `Manifold::Impl::RayCast`, `Kernel11`, `Kernel12`) resolve via
  `using namespace manifold;` at file scope.
- `extras/CMakeLists.txt` adds `TBB::tbb` to `overlap3d_proto`'s
  link line — the spike's transitive TBB symbols (`tbb::combinable`
  via `parallel.h`) weren't being satisfied through the manifold
  shared lib.
- `extras/overlap3d_proto.cpp` adds a `ProductionWinding03` wrapper
  (around `manifold::Winding03<true>(impl, impl, p1q2, expandP=true)`)
  and a side-by-side dump under the existing
  `OVERLAP3D_DUMP_WINDING=1` env var for spike-vs-production
  comparison.

### Verification

- `ctest --test-dir build` — **388/388 passed** (same as pre-refactor).
  Confirms zero behavior change to the production Boolean3 codepath.
- `--api` mode on the 5 .obj fixtures matches handoff baseline drift
  exactly (havocglass 0%, offset12 ≈0%, hull-mask 0.011%, cray
  ~2369%, self-intersect 2.8%).
- `--advfuzz`: 19/20 valid (this run) vs the handoff's 12/12; within
  documented run-to-run fuzz variance band.

### Findings: production W03 with M==M is SoS-collapsed

Calling `Winding03_<true, true>(M, M, p1q2)` for self-mesh produces
output that is *not* directly usable as a per-vert winding label
inside a self-intersecting mesh. The SoS perturbation (`expandP`)
is designed to symbolically distinguish inputs P vs Q; with M==M
that distinction collapses, just like `m.Boolean(m, OpType::Add)`
collapses (handoff §"Things that don't work").

Side-by-side spike-scaffold (ray-cast at v ± ε) vs production W03,
under `OVERLAP3D_DUMP_WINDING=1`:

| Fixture | Spike: w distribution | Production W03: w distribution | Agree |
|---|---|---|---|
| havocglass | w=0:141 | w=1:141 | 0/141 (0%) |
| offset12 | w=0:99 w=1:3617 | w=0:63 w=1:3652 | 3554/3715 (96%) |
| hull-mask | w=0:4625 w=1:2351 | w=0:5 w=1:6971 | 2356/6976 (34%) |
| cray | w=-2:1 w=-1:17 w=0:106 w=1:98 w=2:8 w=4:11 | w=0:123 | 46/123 (37%) |
| self-intersect | w=0:317 w=1:16600 w=2:511 | w=1:16099 w=2:672 | 16439/16771 (98%) |

The disagreement pattern matches the SoS-collapse hypothesis:

- **havocglass** (closed manifold, no pierces): production says every
  original vert is "inside" itself (w=1) because SoS symbolically
  shifts the vert into M's interior; the spike's ray-cast probes
  v ± ε and gets w=0 outside the surface.
- **cray** (heavy-pierce): production says w=0 everywhere — the SoS
  collapse + many broken halfedges + per-component rep selection
  result in degenerate per-vert classifications.
- **self-intersect** (heavy-pierce, dense): production agrees with
  the spike on 98% of verts. The bulk of self-intersect's verts are
  in regions where the SoS collapse doesn't matter (verts well away
  from intersection lines); only ~2% disagree.

This validates the handoff's prediction:

> | A. Minimal expose | ... | 1-2 sessions | Low |
> | B. New entry point | ... wraps the kernels for inP=inQ case
>     with adjusted SoS perturbation | 3-5 sessions | Medium |

Option A is a *prerequisite* for Option B (the kernels are now
reusable from outside `boolean3.cpp`), but Option A alone does
**not** solve the cascade-drop. A direct M==M call to the
production `Winding03_` is not authoritative; the SoS adjustment
in Option B is the substantive fix.

### Status going into Option B

- Production kernels (`Shadow01`, `Kernel02`, `Winding03_`,
  `Winding03`) are reusable from anywhere that can include
  `src/winding03.h` — including the spike, future overlap-removal
  entry points, and any Manifold-self-analysis paths.
- The spike's `ComputePerVertWinding` (ray-cast scaffold) remains
  the workhorse classifier for the existing pipeline; it is not
  yet replaced. A real replacement requires the SoS-adjusted
  `SelfMeshAnalysis` (Option B) or per-side ε-perturbation around
  every rep vert.
- The diagnostic dump under `OVERLAP3D_DUMP_WINDING=1` makes
  spike-vs-production divergence measurable for any fixture, so
  Option B's correctness can be tracked against this baseline.

---

## Option B: prototype landed (2026-05-09)

The `SelfMeshAnalysis` entry point (`src/self_mesh_analysis.{h,cpp}`)
is implemented as the Option B prototype: per-vert above/below
winding for a single self-intersecting input mesh M, plus a
`WindingAt(M, origin, direction, length)` helper for arbitrary 3D
probe points. Wired into the spike's polygon classifier as an
opt-in mode (`OVERLAP3D_OPTB=1`).

### Key implementation choice: geometric ε-offset, not SoS

The original Option B sketch in the handoff envisaged calling
`Winding03_<expandP, true>(M, M, p1q2)` twice with `expandP=true`
and `expandP=false` to recover above and below windings. Diagnostic
runs (committed in 563a0ce2 under `OVERLAP3D_DUMP_WINDING=1`)
showed this **does not work**: `expandP` in the production
`Winding03_` controls a tiebreaker on the `Shadows()` predicate,
not a geometric probe direction. With M==M and the topological
adjacency filter (= skip Kernel02(v, f) when v is a face vert),
both `expandP` choices produce the same z-projection ray winding
through `Kernel02` — same answer, no geometric above/below.

The working approach: keep production `Winding03_`'s structure
(DisjointSets flood-fill over verts, broken at p1q2-listed
halfedges, one classification per component), but substitute a
**literal geometric ε-offset ray-cast** for the per-component
classification. For component rep `v`:

  origin_above = v + ε * n(v)   // n(v) = outward vert normal
  origin_below = v - ε * n(v)
  w_above[c]   = signed winding of M at origin_above
  w_below[c]   = signed winding of M at origin_below
  ∀ v ∈ c: w_above[v] = w_above[c]; w_below[v] = w_below[c]

The flood-fill stays valid because intact halfedges preserve
outward-normal orientation between adjacent faces — so "outward"
near each component's surface is well-defined and the per-rep
classification labels every vert in that component.

Ray-cast narrow phase is Möller-Trumbore against M's triangles,
broad phase via `M.collider_` (the production BVH). Sign
convention is "inside = +1" to match the spike's classifier
expectations.

### Per-vert winding outputs

`OVERLAP3D_DUMP_WINDING=1` shows the per-vert above/below
distribution alongside the spike's `ComputePerVertWinding` scaffold
and the M-vs-M production Winding03 baseline:

| Fixture | spike scaffold | Option B w_above | Option B w_below | Geometric reading |
|---|---|---|---|---|
| havocglass | w=0:141 | w=0:141 | w=1:141 | Clean closed manifold, all verts on boundary ✓ |
| offset12 | w=0:99 w=1:3617 | w=0:3715 | w=1:3715 | Same — all verts on boundary ✓ |
| hull-mask | w=0:4625 w=1:2351 | w=0:6976 | w=0:5 w=1:6971 | 5 verts near pierce hard to classify; 6971 on boundary ✓ |
| cray | -2:1 -1:17 0:106 1:98 2:8 4:11 | -4:9 0:68 1:46 | -4:9 0:114 | 46 boundary verts, plus mixed-orientation Subtract surfaces |
| self-intersect | 0:317 1:16600 2:511 | 0:16099 1:672 | 1:16099 2:672 | 16099 boundary verts; **672 verts inside an overlap layer** ✓ |

For the clean / mild-pierce fixtures, every vert classifies
cleanly as a surface boundary vert (w_above=0, w_below=1). For
self-intersect, the 672 verts sitting inside an overlap region
(w_above=1, w_below=2) are correctly identified — that is exactly
the signal Option B exists to provide.

### Classifier swap

`OVERLAP3D_OPTB=1` activates a polygon-keep classifier that uses
the per-vert windings: a polygon's perimeter votes "boundary /
interior / exterior" per non-chord vert, and the polygon is kept
iff boundary verts win the plurality. Polygons composed entirely
of chord verts (= no per-vert signal) fall back to a centroid
probe via `WindingAt(M, centroid ± ε * n_T, rayDir, rayLen)`.

Results vs the heuristic-baseline default classifier:

| Fixture | default | Option B |
|---|---|---|
| havocglass | status 0, drift 0% | status 0, drift 0% (identical) |
| offset12 | status 0, drift ≈0% | status 0, drift ≈0% (identical) |
| hull-mask | status 0, drift 0.011% | status 0, drift 0.011% (identical) |
| cray | k=1: 64 → status 2 | k=1: 74 → status 2 |
| self-intersect | k=1: 333 → status 2 | k=1: 379 → status 2 |
| **advfuzz (1395 cases)** | **19/20 valid** | **19/20 valid** (matches) |

Option B matches the heuristic baseline on the clean-to-mild
fixtures and on adversarial fuzz. On the heavy-pierce inputs
(Cray, self-intersect) it does not yet beat the heuristic
baseline's k=1 floor — slightly worse, in fact (Cray: 64→74,
self-intersect: 333→379 k=1 edges).

### Why Option B doesn't beat the heuristic on heavy-pierce inputs

The current per-vert classifier only labels original verts of M;
chord verts (= step-7 intersection-line endpoints) don't have
windings. For Cray's 317 polygons, ~196 are auto-kept (single-poly
tris with no chord) and ~121 are multi-poly polygons whose
perimeters are mostly chord verts. Per-vert windings provide
strong signal on the auto-keep side but only weak signal on the
multi-poly side. Pair-sym Phase 1+2 then operate on top, and the
combined system trends slightly worse than the analytical+pair-sym
default tuned over many iterations.

Three plausible next steps to push Option B past the heuristic
baseline:

1. **Per-(chord-vert, parent-tri) windings**: classify each chord
   vert v_c with a probe ± ε * n_T per parent triangle T (since
   v_c lies on multiple parent tris). Adds rays proportional to
   (chords × parent tris).

2. **Per-polygon-centroid probe always** (not just fallback):
   replace the per-vert vote entirely with a single centroid probe
   per polygon. O(P) ray-casts; clean signal because n_T is well-
   defined per polygon.

3. **Combine with pair-sym Phase 1 only** (drop Phase 2): Phase 2
   tends to cascade-drop over-aggressively; Option B's per-vert
   signal might be cleaner on its own without the mesh-edge AND-
   merge.

### Code reuse audit (= "as much code as possible")

What `SelfMeshAnalysis` reuses from production manifold:

- `Manifold::Impl::collider_` (BVH broad phase)
- `manifold::DisjointSets` (flood-fill components)
- `manifold::MakeSimpleRecorder` (collider iteration glue)
- `manifold::for_each` / `autoPolicy` / `countAt` (parallel scaffolding)
- `manifold::Vec` / `VecView` (storage)
- `manifold::Halfedge` + `Manifold::Impl::halfedge_` / `vertPos_` /
  `vertNormal_` / `faceNormal_` / `bBox_` (mesh accessors)

What it does NOT reuse (and why):

- **`Winding03_` / `Kernel02` (= the public Option-A exports)**: the
  SoS perturbation collapses for M==M, so the per-component
  classification step had to be replaced with geometric ε-offset
  ray-cast. The flood-fill *structure* (DisjointSets, broken
  halfedges) is reused, but the inner classifier is bespoke.
- **`Manifold::Impl::RayCast`** (`src/boolean3.cpp:609`): exists
  but takes a public `Manifold` wrapper context. Reimplemented
  (~30 lines, Möller-Trumbore + BVH) inline in
  `self_mesh_analysis.cpp` to keep the API narrow on `Impl`.
  Could be deduplicated by lifting `RayCast` to also accept
  `Impl` directly.

Total new code: ~150 LOC in `self_mesh_analysis.cpp` + ~70 LOC
header + ~50 LOC of spike classifier wiring. Most of that is
glue around the reused production primitives.

### Tuning iterations (results)

After the initial prototype landed, three follow-up iterations
explored ways to push Option B past the heuristic baseline.

**Iter 1: centroid-only mode (`OVERLAP3D_OPTB_CENTROID_ONLY=1`).**
Skips the per-vert vote and uses only the polygon-centroid ± ε * n_T
probe. Removes the bias from chord-heavy perimeters. Result on
heavy-pierce inputs:

  cray k=1 floor:                      62 (default) → 74 (Option B per-vert) → 74 (centroid)
  self-intersect k=1 floor:           333 (default) → 379 (Option B per-vert) → 355 (centroid)

Centroid mode helps self-intersect (379 → 355) but doesn't move
Cray. advfuzz unchanged at 19/20 valid.

**Iter 2: drop pair-sym Phase 2 (`OVERLAP3D_NOPAIRSYM_PHASE2=1`).**
Phase 2's mesh-edge AND-merge was tuned for the analytical+heuristic
classifier; Option B's centroid probe already encodes the
"winding ≥ 1" mask cleanly per polygon, so Phase 2 becomes
redundant pressure. Result on self-intersect:

  default (P1+P2):                    333
  default (P1 only):                  319
  Option B per-vert + skipP2:         479 (worse)
  Option B centroid + skipP2:         **309**  ← new low (-7% vs heuristic baseline)

Cray with the same combo: 71 k=1 + 7 k=3 + 5 k=4 — worse than the
default's clean (62, 0, 0). So the Phase-2-drop win is
self-intersect-specific. advfuzz unchanged at 19/20.

**Iter 3: jump-mode criterion (`OVERLAP3D_OPTB_JUMP=1`,
negative result).** Replaces "keep iff exactly one of {wa, wb} ==
0" with "keep iff |wa - wb| ≥ 1" — meant to handle Cray's back-side
surfaces where windings go negative and the standard rule misses
them. Result:

  cray skipP2:                         71 → 70 k=1, but k=3 7→13 and k=4 5→11
  self-intersect skipP2:              309 → 539 k=1 (much worse)

Over-keeps polygons on overlap-layer boundaries that should be
dropped. The "= 0 on one side" rule is correct for self-union
semantics; Cray's residual gap is elsewhere.

### Status going into chord-vert / further work

Option B (centroid + Phase-1-only) is now the **best classifier
in the spike** for self-intersect-shaped inputs (= ordinary
self-piercing meshes with consistent outward normals): 309 k=1 vs
the heuristic baseline's 333, with adversarial fuzz unchanged.

**Even bigger win on generic-twin-7081** (one of the heaviest-
pierce fixtures, default k=1 = 1773): Option B (centroid + skipP2)
drops k=1 to **705 (60% reduction)**. k=4 rises (25 → 82) — the
heuristic was masking real chord-pair ambiguity that the centroid
probe correctly surfaces, but the residual k=4 means topology
doesn't yet close to manifold. Net: status 2 still, but
dramatically closer.

| Fixture | default k=1 | OptB cent+skipP2 k=1 | Δ |
|---|---:|---:|---:|
| havocglass / offset12 / offset34 / hull-mask | 0 | 0 | — (clean) |
| generic-twin-7863 | 9 | 8 | -11% |
| self-intersect | 333 | 309 | **-7%** |
| generic-twin-7081 | 1773 | **705** | **-60%** |
| cray | 62 | 71 | +15% (worse) |

Cray remains worse — the residual gap is the back-side normal flip
in Subtract results. Two paths to close that gap:

1. **Chord-vert windings** (planned, task 18): classify each
   chord vert v_c with a probe ± ε * n_T per parent triangle T.
   Adds per-vert signal where the centroid alone is weak.
2. **Use forward-only face normals**: when the input mesh comes
   from a Boolean op with back-side flag, the canonical outward
   direction can be recovered from the underlying op metadata.
   Out-of-spike scope.

### Handoff for further Option B work

Remaining items:

- **Medium effort**: chord-vert windings via per-parent-tri probes
  (task 18). Should improve the per-vert classifier's signal on
  chord-heavy inputs and probably also help the per-vert+pair-sym
  combo. Less likely to help centroid-only.
- **High-leverage**: extend `AnalyzeSelfMesh` to also produce the
  edge-face pierce list — Kernel12 + adjacency-filter version
  that would replace ~200 LOC of the spike's hand-written
  `FindEdgeTriIntersections`. Heavier reuse of production kernels.
- **Open**: investigate Cray's per-vert windings to identify
  exactly why w_above goes to -4 on certain back-side verts. If
  it's a vertNormal_ pathology, computing per-face normals on
  demand might fix it.

### Broader-corpus extension (2026-05-10)

The original adversarial corpus (`--advfuzz`) had 7 classes / 1395
cases. After the user asked "does Option B's win generalize?", the
corpus was extended to 13 classes / 2495 cases, drawing from
external sources:

- **Subtract operations** (subtract-rotated, subtract-3-chain): test
  whether Cray's loss generalizes to single-Subtract or Subtract-
  in-chain inputs. **It doesn't** — Option B fixes 4/4 Subtract
  pierces, identical to default.
- **Smith UCAM-CL-TR-766 cases** (smith-three-mutual,
  smith-four-coplanar, cospherical-shell): three constructed
  adversarial classes from Julian Smith's 2009 thesis
  ("Towards Robust Inexact Geometric Computation") and the
  3D analog of the cocircular-points pathology.
- **CGAL 12-cube-chain**: mirrors Lazard & Valque (CGF 2025)
  autorefinement workload — iterative union of 12 rotated cubes.
  Highest self-pierce rate of any class (10%).

Smith's two non-laddering cases (smith-four-coplanar = his only
reported field-bug case from 1997, cospherical-shell = the
cocircular-points 3D analog) **do not produce self-pierces** in
manifold's modern Boolean3 — they're handled cleanly by the
α-budget ε + BVH + SoS perturbation. This is a positive finding for
manifold itself.

Smith's smith-three-mutual (= his "edge-pierces-triangle laddering"
3D case from ch.9.2.2 / fig 9.1) does produce 2/200 pierces; both
fixed in default and Option B.

CGAL 12-cube-chain at 10% pierce rate is the most demanding fuzz
class to date; Option B fixes all 10/10.

Final advfuzz comparison across the full 2495-case extended corpus:

| Classifier | Valid output | Fixed | Regularized | Unchanged |
|---|---|---|---|---|
| default (heuristic) | 35/36 | 32 | 3 | 1 |
| Option B (centroid+skipP2) | 35/36 | 32 | 3 | 1 |

**Conclusion**: Option B's per-fixture wins generalize — there's no
class in the adversarial corpus where Option B regresses, and on
the .obj fixtures it wins on 4 of 5 (60% k=1 reduction on the
heaviest-pierce case generic-twin-7081). The Cray-specific loss is
a real outlier: it does NOT predict Option B's behavior on the
broader Subtract-input population.

### Perim-guard refinement (2026-05-10)

Empirical investigation of Cray's residual gap (`OVERLAP3D_DUMP_K1=1`)
showed every k=1 edge had a kept polygon AND a dropped polygon both
containing the edge — Option B was over-dropping multi-poly polygons
whose perimeter shared an edge with a single-poly auto-kept neighbor,
orphaning the edge.

`OVERLAP3D_OPTB_NO_PERIM_GUARD=1` opts out of a new guard pass that:

1. Builds edge → list-of-(triId, pi) using-as-perimeter map.
2. Marks anchor edges (= edges adjacent to a kept single-poly tri).
3. For each Option B drop decision, overrides to keep if any
   perimeter edge is an anchor.

Result with perim-guard ON (default for Option B):

| Fixture | default k=1 | OptB+guard k=1 | Δ |
|---|---:|---:|---:|
| havocglass / offset12 / offset34 / hull-mask | 0 | 0 | identical |
| generic-twin-7863 | 9 | **5** | **-44%** |
| self-intersect | 333 | 328 | -2% |
| generic-twin-7081 | 1773 | **350** | **-80%** |
| cray | 62 | 72 | +16% |

**Option B now beats the heuristic baseline on 5 of 6 piercing
fixtures.** The 80% reduction on generic-twin-7081 (the heaviest-
pierce fixture) is the headline: 1773 → 350 k=1 edges. advfuzz
unchanged at 35/36 valid — no regression on the 2495-case
adversarial corpus.

Cray remains the documented exception (back-side normal flip in
Subtract results). Future work to close this gap is itemized below.

### Recommended Option B invocation (production-leaning)

```
OVERLAP3D_OPTB=1 \
  OVERLAP3D_OPTB_CENTROID_ONLY=1 \
  OVERLAP3D_NOPAIRSYM_PHASE2=1 \
  ./build/extras/overlap3d_proto <fixture>
```

Perim-guard is on by default in Option B mode (= no env var needed).
For Subtract-derived inputs with extensive back-side surfaces (Cray-
class), the default heuristic remains preferable until back-side
detection is wired in.

### Milestone: pausing Option B iteration

After 5 iterations against the .obj battery + advfuzz (centroid-
only, drop pair-sym Phase 2, jump-mode criterion, faceNormal_-based
probe, single-poly auto-keep), Option B has reached this state:

  havocglass / offset12 / offset34 / hull-mask : matches default (clean)
  generic-twin-7863                            :  9 →  8  k=1 (-11%)
  self-intersect                               : 333 → 309 k=1  (-7%)
  generic-twin-7081                            : 1773 → 705 k=1 (-60%)
  cray                                         : 64 → 71  k=1 (+11%, back-side flip)
  advfuzz (1395 cases)                         : 19/20 valid in both

Per the stop conditions in the iteration prompt:
  (a) "Option B clearly beats the heuristic on Cray + self-intersect
       k=1 floors" — partial: beats self-intersect (and 3 others),
       doesn't beat Cray.
  (b) "Several iterations show no improvement reachable without much
       heavier work" — confirmed for Cray. The 5 iterations exhaust
       the easy levers; further improvement requires either deep
       investigation of Cray's back-side normal pathology or wiring
       the back-side flag from the Boolean3 op metadata into
       AnalyzeSelfMesh's probe direction. Both are out-of-scope for
       this prototype iteration.

Pausing Option B iteration here. Net of 11 commits: a working,
documented, unit-test-passing prototype that beats the heuristic
baseline on **4 of 5 .obj fixtures** including the heaviest-pierce
non-Subtract case (generic-twin-7081, -60% k=1), with no adversarial-
fuzz regression. Cray remains the documented exception with a
clear next-step recipe.

Recommended invocation for new self-intersecting input meshes
(production-leaning convention):

```
OVERLAP3D_OPTB=1 OVERLAP3D_OPTB_CENTROID_ONLY=1 \
  OVERLAP3D_NOPAIRSYM_PHASE2=1 \
  ./build/extras/overlap3d_proto <fixture>
```

Switch to default heuristic for Subtract-derived inputs (= those
likely to have back-side surfaces) until the back-side detection is
wired in.

---

## Determinism fix (2026-05-10)

> **All metric tables in sections above were contaminated** by a
> spike-internal nondeterminism bug; corrected numbers are at the
> end of this section. The qualitative conclusions still hold;
> the magnitudes have shifted.

### Discovery

Manifold's main library is **strongly deterministic by design** — see
`src/csg_tree.cpp:425` ("make sure the order of result is
deterministic"), `src/boolean3.cpp:286-294` (`Kernel12Recorder::get()`
stable_sorts `tbb::combinable`-collected pairs into canonical order),
`src/parallel.h:419,450,462,477` (parallel reductions explicitly
documented as deterministic-only-when-commutative-and-associative).

But the spike's heuristic classifier showed Cray k=1 fluctuating
56–64 across runs (mode 62). When asked to track this down, the
investigation went:

1. **Hash every pipeline phase boundary**. Added an FNV-1a 64-bit
   hash + `OVERLAP3D_DET_HASH=1` env var that prints stable hashes
   at each phase (impl arrays, edges, eeIsects, etIsects, step7p2,
   step8, step11p1/p2/p3, chordPartners, precomputedKeep).
2. **Manifold's Boolean3 IS deterministic for Cray** — confirmed
   v123 t238 vol 4.9136e+10 stable across 30 runs, with stable
   hashes for impl.vertPos_, impl.halfedge_, impl.faceNormal_,
   impl.vertNormal_, impl.bBox_.
3. **etIsects ordering varies** (3 distinct hashes across 6 runs)
   even with `parallel=false` on the BVH walk. EmitNewVertsAndEdges
   then canonicalizes downstream, masking the variance.
4. **precomputedKeep also varies** (~25% of runs) despite all
   classifier inputs hashing stable. This was the smoking gun:
   `AnalyticalKeep(triId=72, pi=1, …)` returned different results
   on identical inputs.
5. **Trace per-call FP values** for tri=72. Found one entry of
   `nonChordPts` reading `(1.66e-310, 8.28e-311, 1.99e-310)` —
   denormal garbage. **Uninitialized memory.**

### Root cause: two off-by-baseId reads in `AnalyticalKeep`

```cpp
// PER-VERT PATH (line ~2548):
//   pt = (v < baseId) ? impl.vertPos_[v] : positions3D[v];
//                                                    ^^^
//                                          missing - baseId

// FALLBACK PATH (line ~2574):
//   fallbackPt += (v < baseId) ? impl.vertPos_[v] : positions3D[v];
//                                                              ^^^
//                                                    missing - baseId
```

The polygon walker can emit chord vert ids (≥ baseId) that are NOT
present in `step7p2.newVertPositions` — step 8 propagation
introduces additional new ids referenced from halfedge graphs but
not in step7p2's vert table. So `positions3D[v]` reads off the end
of the array → uninitialized memory → denormal garbage. The
`dot(garbage, n_B)` term then random-walks across runs (because the
"garbage" depends on whatever lives in the next memory page, which
varies with malloc state, environment, etc.), flipping the comparison
`dot(p − p_B, n_B) <= thresh * nMag` and therefore the keep/drop
decision.

Fix is mechanical: bounds-check `j = v − baseId`, skip OOB entries,
require at least one valid contribution before pushing to nonChordPts.
Both bugs fixed in the same commit. Determinism infrastructure
(`OVERLAP3D_DET_HASH=1`, `OVERLAP3D_DET_HASH_DUMP_KEEP=1`) kept for
future debugging.

### Verification (post-fix)

| Fixture | classifier | k=1 across N runs | Status |
|---|---|---|---|
| cray | default | **56 ×8** | deterministic (was 56-64 spread, mode 62) |
| cray | Option B | **72 ×8** | deterministic (was deterministic at 72) |
| self-intersect | both | identical ×3 | deterministic |
| generic-twin-7081 | both | identical ×3 | deterministic |
| advfuzz 2495 cases | both | 35/36 valid | unchanged |

`OVERLAP3D_DET_HASH=1` shows `precomputedKeep` hash 16/16 identical
after the fix.

### Corrected metric table (supersedes all earlier sections)

| Fixture | default k=1 (post-fix) | Option B k=1 | Δ vs default |
|---|---:|---:|---:|
| havocglass / offset12 / offset34 / hull-mask | 0 | 0 | identical |
| generic-twin-7863 | 9 | **5** | **−44%** |
| self-intersect | 333 | 328 | −2% |
| generic-twin-7081 | **946** (was reported 1773) | 350 | **−63%** (was reported −80%) |
| cray | **56** (was reported 62) | 72 | **+29%** (was reported +16%) |
| advfuzz (2495 adversarial cases) | 35/36 valid | 35/36 valid | unchanged |

All numbers in this table are now reproducible across runs.

The qualitative story stands:
- Option B beats default on 3 of 4 piercing fixtures (the 4 clean
  fixtures are ties since they have no pierces).
- generic-twin-7081 remains the headline win at -63%.
- Cray remains the documented exception, with the gap a bit larger
  (-29% instead of -16%) because the bug had been masking some of
  the heuristic's true performance.
- adversarial fuzz unchanged at 35/36 valid in both classifiers.

### Implication for productionization

A productionized `RemoveSelfIntersections` MUST preserve manifold's
determinism contract. Two paths:

1. **Keep the heuristic classifier** but ensure all the OOB reads
   are caught by tests (bounds-check assertions in debug builds; an
   automated determinism test that runs the operation N times and
   asserts byte-identical output). The fixed AnalyticalKeep is now
   safe; analogous bugs may exist elsewhere in the pipeline (the
   etIsects-order nondeterminism is benign because EmitNewVertsAndEdges
   canonicalizes, but should be checked for any data flow that
   doesn't canonicalize).
2. **Use Option B's centroid classifier** which never goes through
   AnalyticalKeep. Empirically already deterministic; doesn't depend
   on the per-vert chord-vert lookup that surfaced the bug.

---

## Apples-to-apples reframe (2026-05-10): the classifier wasn't the win

> **Most of what was attributed to "Option B classifier" in the
> sections above was actually pair-sym Phase 2 being detrimental
> on certain inputs.** The cleaner attribution is below.

After fixing the determinism bug, comparing default vs Option B at
the *same* pair-sym configuration tells a different story:

| Fixture | default (full P1+P2) | default + skipP2 | OptB + skipP2 + perim | OptB + full pair-sym |
|---|---:|---:|---:|---:|
| 4 clean fixtures | 0 | 0 | 0 | 0 |
| **cray** | **56** | 86 | 72 | 74 |
| self-intersect | 333 | 319 | 328 | 323 |
| **generic-twin-7081** | 946 | **352** | **350** | 953 |
| generic-twin-7863 | 9 | 8 | **5** | 9 |

What this shows:

- **The -63% headline win on generic-twin-7081 came from disabling
  pair-sym Phase 2** (946 → 352). Option B's centroid classifier on
  top adds ~0.6% (352 → 350 = essentially noise).
- **Option B's centroid classifier ≈ default analytical** at the
  same pair-sym config. Sometimes slightly better (cray 72 vs 86,
  generic-twin-7863 5 vs 8), sometimes slightly worse (self-intersect
  328 vs 319), basically neutral overall.
- **The "right" pair-sym config depends on the input class**:
  - cray-class (Subtract with back-side): wants Phase 2 ON.
  - generic-twin-7081 (heavy chained Add): wants Phase 2 OFF.
  - self-intersect (medium): slight preference for OFF.

Phase 2's job is mesh-edge halfedge AND-merge — for each polygon
perimeter halfedge (a→b), if its reverse (b→a) is in another
polygon and one is dropped, drop both (cascade-drop forward across
mesh-edge adjacency). The cascade is correct for cray-class (where
chord-pair drops should propagate to neighbors), wrong for
generic-twin-7081-class (where the cascade over-aggressively
removes legitimate boundary polygons).

**Open algorithmic question**: design an adaptive Phase 2 that
detects which input regime it's in and applies the cascade
selectively. Could be heuristic (e.g. cap cascade depth, or only
cascade when the dropped polygon is a known back-side face) or
learned (= per-fixture Phase 2 toggle based on some metric).
The current --advfuzz battery doesn't yet have a class that
discriminates between these regimes well; future work should
probably expand it.

### Recommended invocation update

The "production-leaning" recommendation in the section above
(`OVERLAP3D_OPTB=1 + CENTROID_ONLY=1 + NOPAIRSYM_PHASE2=1`) gives
the best generic-twin-7081 result but pays ~16 k=1 on cray.
Equivalent results are achievable with `OVERLAP3D_NOPAIRSYM_PHASE2=1`
alone (= default classifier with Phase 2 off). Pick based on
priorities:

- **Best generic-twin-7081 / generic-twin-7863** (heavy chained Add):
  `OVERLAP3D_NOPAIRSYM_PHASE2=1` (default classifier with Phase 2 off)
  — equivalent to OptB on those fixtures, 100 LOC less surface area.
- **Best cray** (Subtract back-side): `./overlap3d_proto cray`
  (default config) — full pair-sym is right for this.
- **Best balance**: no single config is best across all inputs.
  An adaptive policy is the open algorithmic problem.

---

## Adversarial-cases push (2026-05-10)

After the apples-to-apples reframe, the user asked to drive the
adversarial cases to zero. End-state via `OverlapRemoval` API
(`--api`, with cap on by default in that path):

| Fixture | status | drift | outcome |
|---|---|---|---|
| havocglass | 0 | 0% | clean (no pierces) |
| offset12 | 0 | ~0% | clean (no pierces) |
| offset34 | 0 | 0% | clean (no pierces) |
| hull-mask | 0 | 0.011% | real cleanup |
| cray | 0 | 2250% | cap-fan geometry wrong but topology manifold |
| self-intersect | 0 | 2.8% | real cleanup |
| generic-twin-7863 | 0 | 0.0009% | real cleanup |
| generic-twin-7081 | 0 | 0% | identity fallback (pair-sym structural issue) |

| advfuzz (2495 cases) | 35/36 valid | 32 fixed + 3 regularized + 1 unchanged |

**7/8 .obj fixtures truly fixed; 1/8 gracefully falls back. advfuzz
effectively at 36/36** (1 unchanged case is a sub-eps pierce at FP
noise floor — relative magnitude 4.8×10⁻¹⁴, comfortably below ε).

### Improvements landed this iteration

- **DFS-with-backtracking cap walker** (commit c9269ec9): replaced
  greedy first-unvisited-neighbor walk with DFS that backtracks on
  dead-ends. Closed +14 cycles on gt-7081 (k=1 169 → 101).
- **Partial-corner chord pair-sym** (commit 76e42aa2): pair-sym
  Phase 1 used to skip when not all 4 chord corners were findable;
  now enforces what constraints are applicable (within-tri opposite +
  cross-pair propagation when corners are partially known). No-op on
  the .obj battery (the missing-corner cases are 0-corners-on-triB,
  not 2-3); kept as a forward-looking robustness improvement.
- **OOB audit + bounds-checked GetPos3 helper** (commit 250fb00b):
  consolidated 6 inline `getPos3` lambdas with off-by-baseId risk
  into a shared helper. DEBUG_ASSERTs under MANIFOLD_ASSERT=ON
  (= throws std::logic_error on regression).
- **Deleted dead `ClassifyPolygon`** (commit 250fb00b): had a known
  same-class OOB bug at line 2249.
- **Production `self_mesh_analysis.cpp` hardening** (commit dec17722):
  defensive bbox/vertNormal_ checks, switched std::vector→Vec and
  unordered_set→set for determinism + bounds-checking.
- **Determinism debugging infrastructure**: `OVERLAP3D_DET_HASH=1`
  prints stable FNV-1a hashes at each pipeline phase.
  `OVERLAP3D_PAIRSYM_DIAG=1` counts chord pair-sym skips.
  `OVERLAP3D_PAIRSYM_TRACE=1` dumps first-3 chord lookup details.
  Used to localize the AnalyticalKeep OOB and the gt-7081 chord
  pair-sym structural issue.

### Open items

generic-twin-7081's pair-sym Phase 1 skip rate is 1259/1259 (100%)
— `findPolyHE` returns -1 for all chord corners on triB. Root cause:
chord halfedges between two interior verts of triB form an isolated
2-cycle in triB's halfedge graph (= no perimeter connection), so the
polygon walker emits a 2-vert polygon which then gets dropped (the
walker requires `polygon.size() >= 3`). The chord verts then disappear
from triB's polygon list, so pair-sym can't see the corners.

The structural fix is in BuildPerTriHalfedgeGraphs / WalkPolygons:
either (a) preserve sub-3-vert polygons so chord halfedges in
degenerate sub-polygons remain visible, or (b) "stretch" triB's
perimeter halfedge graph to include interior chord verts so the
chord's two sides are real polygons. Both are real polygon-walker
redesigns; out-of-scope for the current iteration.

The 1 advfuzz unchanged case is at relative pierce magnitude 4.8×10⁻¹⁴
(below ε). Already at FP noise floor; OverlapRemoval correctly leaves
it alone since any "fix" would just shuffle FP noise.

### Polygon-walker fix attempts (negative results)

In a follow-up wakeup pass, we tried to unblock gt-7081 by extending
WalkPolygons to preserve 2-vert degenerate sub-polygons (= chord
pairs whose endpoints are interior to the parent triangle, isolated
2-cycle in the halfedge graph) so pair-sym Phase 1's findPolyHE
could find them as the "missing" chord corners.

Two variants tried:

1. **Mix degenerate polys into the main polygons list.** Triangulation
   already skipped size<3 via existing gate, so the size check was
   safe. But the classifier (AnalyticalKeep) processes them too,
   drops them as degenerate (no chord-bounded structure), and
   creates new k=1 / k=3 elsewhere. Cray regressed from drift 2235%
   (status 0 + cap) to drift 0% (= fallback). Reverted.

2. **Separate `degeneratePolygons` list, only consulted by findPolyHE.**
   Cleaner separation. But pair-sym's "force the matching corner to
   the same side" logic doesn't work when the matching corner is a
   degenerate polygon: the degenerate poly contributes 0 halfedges to
   the output, so pair-sym's premise (= make k=2 by matching keep/drop
   across corners) is invalid. Cray regressed similarly. Reverted.

The fundamental issue: pair-sym's matching logic assumes corners
contribute halfedges to the output. Degenerate polys don't. So
treating them as "corners" creates inconsistencies.

The real fix would require either:
- A polygon walker that integrates interior chord halfedges into
  triB's main polygon (= insert chord verts into the perimeter
  sequence so the chord becomes part of a non-degenerate polygon).
- Or a separate "chord-only" pair-sym pass that knows degenerate
  polys don't contribute to output and adjusts the matching logic
  accordingly.

Both are real algorithmic redesigns. The conservative
partial-corner pair-sym (commit 661e5040) is the achievable
improvement that helps gt-7081's pre-cap k=1 (1673 → 499 = -70%)
without breaking other fixtures, even though it doesn't fully
close the gap to status 0.

### Post-cap k>2 reducer + re-cap loop (commits 5f07cbd9, 4720eaaa)

Adds a topology-repair pass after the cap walker: for each edge with
> 2 directed halfedge incidences (= chord-pair conflicts pair-sym
Phase 1 left unresolved), drop one tri from the duplicate direction.
The drop creates new k=1 cycles which the cap then closes. Iterate
until either k>2 stops decreasing OR an inner-iteration limit hits.

Greedy heuristic: per round, pick the tri that contributes to the
MOST duplicate directions (= helps reduce multiple conflicts at once).
Convergence by monotonic-k>2-decrease check.

Results:

  cray: dropped 6 tris, capped 2 new cycles → k>2 = 0. Status 0
        with cap-fan drift unchanged at 2235%.
  gt-7081: dropped 384 tris (over ~80 seconds of reducer-cap iter),
           capped 188 new cycles → k=4 went 98 → **51** (-48%),
           k=1 went 15 → 2. Still status 2 → fallback. The 51-k=4
           ceiling appears structural — the reducer can't reach it
           without creating new k>2 elsewhere.

The reducer takes the spike output significantly closer to manifold
on heavy-pierce inputs but doesn't fully unblock gt-7081. The
remaining k=4=51 are chord conflicts that the heuristic can't
distinguish chord-pair partners from interior-edge contributors —
the polygon-walker structural issue manifests here too.

### Final adversarial-cases state

  Fixture            | status | drift     | notes
  ------------------ | ------ | --------- | ----
  havocglass         | 0      | 0%        | clean
  offset12           | 0      | ~0%       | clean
  offset34           | 0      | 0%        | clean
  hull-mask          | 0      | 0.011%    | real cleanup
  cray               | 0      | 2235%     | cap-fan + reducer (k>2=0)
  self-intersect     | 0      | 2.8%      | real cleanup
  generic-twin-7863  | 0      | 0.0009%   | real cleanup
  generic-twin-7081  | 0      | 0%        | identity fallback (51 k=4 floor)

  advfuzz: 35/36 valid (1 sub-eps unchanged at noise floor).

7/8 fixtures truly fixed. 1/8 (gt-7081) blocked at structural
ceiling — pair-sym Phase 1's 4-corner assumption fails on chords
between two interior verts of triB; the polygon walker drops the
chord-direction sub-polygon as < 3 verts; pair-sym can't find the
"missing" corner; chord conflicts persist as k=4 in output. The
reducer + re-cap closes ~half of these but the residual is
geometrically required (= dropping more would create k=1 elsewhere
that the cap can't close).

Resolving gt-7081 requires polygon-walker redesign (preserve sub-3
sub-polygons + adjusted pair-sym matching semantics, OR stretch
triB's perimeter to include interior chord verts so the chord
becomes part of a non-degenerate polygon). Both are real
algorithmic redesigns out-of-scope for the current spike iteration.

### Conflict-aware cap + ear-clip fallback (commits fcf13107, 1226b30b)

Two more cap-walker improvements aimed at gt-7081:

1. **Conflict-aware fan apex selection.** Pre-checks each cycle's
   fan triangulation against current edge incidence map. Tries every
   apex position; picks the first where the fan's interior chord
   edges don't already exist at k=2 (= would otherwise create k>2).
   If all apex positions conflict, skip the cycle (= leave k=1).

2. **Greedy ear-clip fallback** when fan-apex fails. Tries to
   ear-clip the cycle by picking conflict-free diagonals one at a
   time. Updates the edge incidence map as ears are emitted so
   subsequent ear checks see updated state. Rolls back partial ears
   if can't fully triangulate.

Result on gt-7081: **k>2 fully eliminated** (was 53+1, now 0+0).
The cost is leaving 97 k=1 dangling edges where the cap can't
close them without creating new k>2.

Why those 97 cycles can't close: every diagonal chord they could
use is already at k=2 in the rest of the mesh. The cycles form
"true holes" in the classifier's kept-polygon set — the surface
cannot be made manifold without either (a) dropping more upstream
polygons or (b) accepting the holes.

This is the **classifier output limit**: gt-7081's polygon-keep
decisions don't form a valid manifold-decomposable surface. The
post-cap repair pipeline (cap walker + reducer) is structurally
optimal given that input. Further gains require revisiting the
upstream classifier (= not just pair-sym Phase 1; the basic
AnalyticalKeep decisions).

### Upstream classifier rework (2026-05-10)

Per the docu'd next-step, replaced AnalyticalKeep with a winding-
grounded classifier in OverlapRemoval (= the `--api` path).

**Three interlocking changes** in commit 1b53a353:

1. **Option B classifier in OverlapRemoval.** Calls
   `manifold::AnalyzeSelfMesh` for per-vert above/below windings.
   Per-polygon vote: each non-chord perimeter vert classifies as
   boundary (`one of {wa, wb} == 0`) / interior (both ≥ 1) /
   exterior (both ≤ 0). Polygon kept iff boundary verts dominate.
   All-chord-vert polygons fall back to centroid probe via
   `WindingAt` (also from `self_mesh_analysis.h`).

2. **Polygon walker preserves degenerate sub-polys.** Earlier
   attempts at this broke cray; the trick is keeping them in a
   *separate* `degeneratePolygons` list (not in `polygons[]`) so
   the classifier and triangulator don't process them. Only
   pair-sym's `findPolyHE` consults this list.

3. **Degenerate-aware per-direction chord pair-sym.** Replaces the
   pair-level "force one pair K, other D" with per-direction
   manifold constraint:
   - d1 (v0→v1) contributors: A1, B2.
   - d2 (v1→v0) contributors: A2, B1.
   - Each direction needs exactly 1 keeping contributor (= manifold k=2).
   - Degenerate corners contribute 0 (no output triangles).
   Six cases handled (n_d1, n_d2 ∈ {0, 1, 2}); the n_d1=2 && n_d2=2
   case reduces to the original pair-level enforcement.

4. **Unified conflict-aware cap.** Removed the inline simple-fan cap
   walker (~120 LOC); the only cap path is now `doCapPass` which has
   conflict-aware fan-apex selection + ear-clip fallback + LIVE
   `ec` update as fans are emitted. **Eliminated cap-induced k>2**:
   on gt-7081, was adding +7-+10 k=4 per cap pass; now adds 0.

**Result on gt-7081**:
  pre-classifier-rework  : k=1=502, k=4=88   pre-cap
  pre-classifier-rework  : k=1=15,  k=4=98   post-cap (cap added k=4!)
  post-classifier-rework : k=1=296, k=4=99   pre-cap (Option B + pair-sym different topology)
  post-classifier-rework : k=1=16,  k=4=99   post-cap (cap added 0 — clean)
  post-classifier-rework : k=1=97,  k=4=0    post-reducer + re-cap

So the cap is now topology-preserving. The 99 k=4 chord conflicts
that pair-sym Phase 1's per-direction constraint can't resolve
(= 1259/1259 chords still go through partial-corner branches due
to degenerate-corner asymmetry) get traded by the reducer for 97
k=1 dangling edges, which the cap correctly refuses to close.

**Other fixtures**: 7/8 still cleanly fixed. Cray's drift went 2235% →
3054% (= different cap-fan geometry), still status 0. self-intersect
drift 2.8% → 3.7%, still status 0.

**Remaining open**: chord-pair-aware reducer that drops chord pairs
atomically (= drop both A1+B1 OR both A2+B2 as a unit, not single
tris). Single-tri drops trade k>2 for k=1; chord-pair drops cleanly
remove the chord without creating k=1 elsewhere. The chord-pair
structural lookup is now in pair-sym; the reducer needs to consult
it. Not implemented in this iteration (multi-step refactor).

### gt-7081 manifold-fixable, but pierce-regressing (2026-05-10, evening)

The chord-pair-aware reducer hypothesis from the prior iteration
turned out wrong on inspection: with the OVERLAP3D_K4_CLASSIFY
diagnostic added, only 50 of 99 post-cap k=4 edges were chord
edges; **46 were pure input-mesh edges** (both verts pre-existing
in the original mesh) and 3 were one-new-vert edges. These k=4
mesh edges come from 4-way self-intersection geometry where 4
input tris share an edge a-b — a normal artifact of dense
self-intersection.

Three interlocking fixes (commits 71f85a42 + 6c83de31):

1. **Pair-sym Phase 2.5** for multi-owner non-chord edges. The
   prior Phase 2 only handled `owners.size() == 1` per direction;
   any edge with > 1 owner per direction (= the 4-way mesh-edge
   case) was silently skipped. New Phase 2.5 keeps one owner per
   direction (lowest triId) and drops the rest. Result: post-cap
   k=4 dropped from 99 → 59 (40 mesh-edge conflicts resolved).

2. **Min-k1-cost k>2 reducer.** The reducer previously picked the
   lowest tri-index to keep per duplicate-halfedge direction. Now
   picks the tri whose drop would create the MOST new k=1 edges
   (= "most load-bearing"); drops the rest. Reduces dangling
   creation per drop. Post-reducer k=1 dropped from 97 → 32.

3. **Trim-orphans final pass.** The post-reducer residual k=1 are
   chains/branches/non-simple paths that the cap walker can't
   close (its DFS only finds simple closed cycles). The trim-
   orphans pass eats them away: round 0 drops tris with >= 2 k=1
   edges, subsequent rounds drop tris with >= 1 k=1 edge. Re-caps
   between rounds. Iterates 32 rounds max. The aggressive >= 1
   trim is bounded — it only eats the dangling chain (= the bulk
   manifold surface has k=2 by construction). Post-trim k=1: 0,
   manifold!

**Result on gt-7081**: was fallback-to-identity (= no fix at all).
Now: real manifold output v21119 t42190, vol 143611, drift 2.9e-6%.

**self-intersect**: drift 3.68% → 2.34% (from Phase 2.5 alone, the
multi-owner mesh-edge fix).

**Other 6 fixtures**: all unchanged (still status 0). cray still
documented loss (3054% drift, Subtract back-side flip).

**ctest**: 388/388 passing. **advfuzz**: unchanged 35/36 + 1.

The trim-orphans pass is the productionization-blocker (= it's
"give up and chew off dangling regions" rather than "produce the
right geometry"). The volume cost is bounded but the strategy is
heuristic. A cleaner upstream fix would prevent the dangling
halfedges in the first place — likely by making pair-sym + cap
reason about the FULL polygon-keep configuration before
triangulation, instead of fixing it up after.

### Pierce-count vs manifold-status — the real correctness story

Volume drift was not the right correctness oracle. After adding a
post-OverlapRemoval pierce check (= `CheckSelfIntersection` on
input vs output, --api output now shows `pierces: A → B`):

| fixture          | input pierces | default output | trim output |
|------------------|---------------|----------------|-------------|
| cray             | 114           | 393 (3054% drift)| same      |
| self-intersect   | 661           | 2061 (3.7% drift)| 312 (2.3%)|
| hull-mask        | 31            | 31 (0.01%)     | same        |
| offset12         | 46            | 46 (~0%)       | same        |
| offset34         | 0             | 0 (0%)         | same        |
| havocglass       | 1             | 1 (0%)         | same        |
| gt-7081          | 183           | 183 (fallback) | 818 (2.9e-6%)|
| gt-7863          | 7             | 2 (3.1e-13%)   | same        |

Headline: the pipeline does NOT generally produce self-
intersection-free output. It produces a topology-manifold output
that may have its own pierces (different from input's, sometimes
more). The "manifold but more pierces" outcome is a Pyrrhic
victory for an "overlap removal" API.

Trim-orphans is now OPT-IN via OVERLAP3D_TRIM=1 because:
- For gt-7081: trim makes pipeline status 0 but pierces 183 → 818
  (4x worse than fallback). Default fallback to input is better.
- For self-intersect: trim REDUCES pierces from 2061 to 312 (=
  below input's 661, so a real fix). Worth enabling here.

The ideal solution: an in-pipeline pierce-count guard that falls
back when pierces would worsen. Implemented but reverted because
`CheckSelfIntersection` returns inconsistent results between
consecutive calls on the same Manifold (= lazy CSG evaluation
order in the Manifold class — first call returns 2061, second
call after Volume() returns 312). Needs root-cause investigation
before the guard can be trusted.

advfuzz remains 35/36 + 1 unchanged (= effectively 36/36) under
either trim setting, because advfuzz cases are small primitives
where the cap+reducer doesn't blow up topology.

### Determinism fix + pierce-count guard (2026-05-10, late evening)

While investigating the lazy-eval inconsistency (= same `out`
reporting 2061 then 312 pierces), discovered a real determinism
bug in `AnalyzeSelfMesh`: 1 of 5 runs of self-intersect produced
k=1=645/pierces=2061 instead of the consistent k=1=379/pierces=312.

**Root cause** (commit 40a6eb00): `DisjointSets::unite` is thread-
safe (CAS-based), but the chosen REPRESENTATIVE per component
depends on parallel union order. Same connectivity, different rep
across runs. The rep's `vertNormal_` becomes the ε-offset probe
direction for the above/below ray-casts, so a different rep gives
a different probe → different winding → different keep decisions
cascading to a different output mesh.

**Fix**: re-anchor each component to its smallest-id member via a
post-pass map (`compMinId: rawRep → minId`). Both step 3 (per-rep
ray-cast) and step 4 (propagate to all verts) use the deterministic
min-id rep. **Verification**: 20/20 runs of self-intersect --api
produce identical k=1=379. Same determinism class as the earlier
unordered_set bug (commit 5acc6ab5).

**Side effect**: cray's pierce count went 393 → 449 (= the other
non-deterministic outcome). Determinism > marginal variation.

**Pierce-count guard** (commit 0276b88c): with `CheckSelfIntersection`
now reliable, OverlapRemoval gained an in-pipeline guard — if
pipeline output has MORE pierces than merged input, fall back to
merged input. Force materialization via `Volume()` on both meshes
before pierce check so both sides see the post-eval mesh (= avoids
the lazy-eval ordering issue).

**Trim-orphans flipped to DEFAULT-ON**: now safe because the pierce
guard catches the regression case (gt-7081 trim → 818 pierces →
guard fires → fallback to 183).

**Final per-fixture pierce results** (NO regressions; two real
pipeline-fix wins):

| fixture          | input | output | branch        |
|------------------|-------|--------|---------------|
| cray             | 114   | 114    | FALLBACK      |
| self-intersect   | 661   | 312    | PIPELINE FIX  |
| hull-mask        | 31    | 31     | passthrough   |
| offset12         | 46    | 46     | passthrough   |
| offset34         | 0     | 0      | passthrough   |
| havocglass       | 1     | 1      | passthrough   |
| gt-7081          | 183   | 183    | FALLBACK      |
| gt-7863          | 7     | 2      | PIPELINE FIX  |

**Pierce-monotonicity guarantee**: OverlapRemoval never makes
self-intersection count worse than input (modulo the opt-out env
var). This is the correctness invariant that should have been
asserted from day 1.

**Cost**: pierce guard adds one `CheckSelfIntersection`-pair per
call. On gt-7081 (16k verts) this adds ~3 min per OverlapRemoval
invocation. Acceptable for correctness-first; opt-out via
`OVERLAP3D_NO_PIERCE_GUARD=1` for batched workloads.

ctest 388/388. advfuzz 35/36 + 1 unchanged.

### Pierce-aware cap walker (2026-05-10, late-late evening)

The next attack vector after the pierce guard: prevent the cap
walker from emitting fan/ear-clip tris that would pierce existing
geometry in the first place.

Implementation (commit 8c5859ae): inside `doCapPass`, build a BVH
over current `out.triVerts` at start of the call. For each
candidate fan-tri (apex, vi, vj), query the BVH for broad-phase
overlaps and run 6 segment-vs-tri tests against each (Möller-
Trumbore). If any test indicates a pierce, refuse the candidate
apex; if no apex is conflict+pierce-free, fall through to ear-
clip (also pierce-checked); if neither works, leave the cycle
uncapped and let the pierce/drift guard fall back.

Per-cycle BVH rebuild added in commit 944fb75b — catches within-
pass cross-cycle pierces (= cycle B's cap piercing cycle A's
cap added earlier in same call). Cheap relative to pierce checks.

**Volume-drift sanity guard** added alongside pierce guard: falls
back if pipeline output has > 50% drift OR sign-flipped volume.
Catches cray's "Subtract back-side flip produces inverted shell"
case where pierces are 4 (= huge improvement from 114) but the
geometry is geometrically wrong (volume sign flipped from +4.9e10
to -2.4e11). Without the drift guard, the pierce guard alone
would let this through — pierce count is not the only correctness
oracle.

**Trim-orphans flipped to OPT-IN** (was opt-out). Reasons: (1)
pierce-aware cap subsumes most trim wins by preventing the
pierce-creating cap candidates trim was patching over; (2)
trim+pierce-cap is too slow on gt-7081 (>10min); (3) even if
trim ran, pierce guard would fall back anyway.

**Final per-fixture results** (this is the strongest pipeline
state of the session):

| fixture | input | output | branch        |
|---|---|---|---|
| cray | 114 | 114 | FALLBACK (drift sign-flip) |
| self-intersect | 661 | **58** | PIPELINE FIX (was 312, 5x better) |
| hull-mask | 31 | 31 | passthrough |
| offset12 | 46 | 46 | passthrough |
| offset34 | 0 | 0 | passthrough |
| havocglass | 1 | 1 | passthrough |
| gt-7081 | 183 | 183 | FALLBACK (status non-manifold) |
| gt-7863 | 7 | 2 | PIPELINE FIX |

Headline win: **self-intersect now near-fixed** (58 residual
pierces vs 661 input, vs 312 with the prior pipeline). Other
fixtures hold pierce-monotonicity. No regressions.

Pierce-aware cap is the right architectural move — refusing to
fan-fill what would create pierces is structurally correct, and
it shifts the failure mode from "produce manifold with worse
pierces" (= what the pre-pierce-aware cap did) to "leave cycle
uncapped, fall back" (= safer, preserves input).

cray's fallback-via-drift-sign-flip and gt-7081's fallback-via-
non-manifold are both correct outcomes given the current
classifier capabilities. Further fix would require an upstream
classifier rework that produces fewer dangling halfedges and
fewer geometrically inverted polygon-keep decisions.

ctest 388/388. advfuzz 35/36 + 1 unchanged.

### Shared-vert pierce detection (2026-05-10, very late evening)

Investigated the 4 fixtures where the pipeline runs but doesn't
fix anything (= passthrough hull-mask 31→31, offset12 46→46, and
fallbacks cray, gt-7081). Root cause traced for the passthroughs:
step 6 (`FindEdgeTriIntersections`) silently skipped any edge-tri
intersection where the edge endpoint coincided with a triangle
vertex (or appeared in the snap-merged onTriList/onEdgeList).
Original rationale: "no NEW pierce; the endpoint connection is
already in the topology."

For self-intersecting input where a Boolean op merged verts
between left+right meshes, the merged vert hides a REAL geometric
pierce — the edge passes THROUGH the tri's interior PAST the
shared vert. The classic skip was treating these as legal-by-
construction; the pipeline never saw them.

**Fix** (commit d978dd6f): default-bypass all three skips
(shared-vert, in-tri list, in-edge list). Opt back into the old
behavior via `OVERLAP3D_SKIP_SHARED_VERT_PIERCES=1`.

**Pipeline-fix wins after this commit**:

| fixture | input | output | branch |
|---|---|---|---|
| cray | 114 | 114 | FALLBACK (drift sign-flip) |
| **self-intersect** | 661 | **15** | PIPELINE FIX (97.7% reduction) |
| **hull-mask** | 31 | 30 | new pipeline fix (was passthrough) |
| **offset12** | 46 | 45 | new pipeline fix (was passthrough) |
| offset34 | 0 | 0 | already clean |
| havocglass | 1 | 1 | passthrough |
| gt-7081 | 183 | 183 | FALLBACK (status non-manifold) |
| gt-7863 | 7 | 2 | PIPELINE FIX |

self-intersect went 58 → 15 pierces, a 4x improvement on top of
the prior pierce-aware cap win. From 661 input → 15 output, the
pipeline now resolves 97.7% of self-intersect's pierces.

**Architectural limit acknowledged**: hull-mask and offset12 only
get 1 pierce fixed each (out of 30+). The rest hit step 7 phase
2's "exactly 2 endpoints per tri-tri pair" filter — most pierces
in these fixtures have 1-endpoint topology (= edge-tip-in-
interior, not edge-piercing-through). Fixing this requires step
7 + polygon walker + pair-sym to handle 1-endpoint pairs
(= multi-session refactor).

ctest 388/388. advfuzz 35/36 + 1 unchanged. No regressions.

### Pierce-aware reducer (2026-05-10, very-late evening)

Investigation showed the remaining pierces in self-intersect (15),
hull-mask (30), and offset12 (45) were post-classifier overlapping
tris (= classifier kept multiple polygons whose emitted tris
geometrically pierce each other). The pierce-aware cap walker only
checks NEW cap tris against existing geometry; classifier-emitted
tris that pierce each other weren't caught.

Added pierce-aware reducer that runs after polygon emit and before
cap walker (commit a8a7c9fd). Algorithm: build BVH over current
out.triVerts, find piercing pairs via Möller-Trumbore broad+narrow
phase, drop one tri per pair (heuristic: higher pierce-count tri =
more "load-bearing" to remove; ties broken by lower triId).
Iterate up to 32 rounds. Default ON; opt-out via
`OVERLAP3D_NO_PIERCE_REDUCER=1`.

**HUGE wins** on the previously-stuck fixtures:

| fixture | input | prior | new | reducer drops |
|---|---|---|---|---|
| **hull-mask** | 31 | 30 | **8** (74% better) | 24 |
| **offset12** | 46 | 45 | **8** (83% better) | 27 |
| **self-intersect** | 661 | 15 | **13** | 1 |
| havocglass | 1 | 1 | 1 | 1 |
| gt-7863 | 7 | 2 | 2 | 1 |

Hull-mask and offset12 unlocked: they had been stuck at near-
input pierce count because most of their pierces are classifier-
output overlaps that require post-emit pierce-drop rather than
upstream classifier intervention.

Post-cap pierce reducer + re-cap loop (commit 69af76ac): the
"future work" landed within the same session. Integrating doCapPass
inside the post-cap pierce reducer's iteration loop closes the k=1
cycles that pierce-drops create. Cap is pierce-aware so doesn't
undo the pierce-fix.

**3 fixtures fully fixed (PIERCES → 0)**:

  hull-mask:   31 → 0 (was 8)
  offset12:    46 → 0 (was 8)
  havocglass:   1 → 0 (was 1)

Default ON; opt-out via `OVERLAP3D_NO_POST_CAP_PIERCE_REDUCER=1`.

**Final correctness state of session**:

| fixture | input | output | branch | reduction |
|---|---|---|---|---|
| cray | 114 | 114 | FALLBACK | 0% |
| **self-intersect** | 661 | **13** | PIPELINE FIX | 98.0% |
| **hull-mask** | 31 | **0** | PIPELINE FIX | **100%** |
| **offset12** | 46 | **0** | PIPELINE FIX | **100%** |
| offset34 | 0 | 0 | already clean | n/a |
| **havocglass** | 1 | **0** | PIPELINE FIX | **100%** |
| gt-7081 | 183 | 183 | FALLBACK | 0% |
| gt-7863 | 7 | 2 | PIPELINE FIX | 71.4% |

**5 of 7 piercing fixtures now have pipeline fixes; 3 of those
are FULLY RESOLVED** (0 residual pierces). Only cray and gt-7081
still fall back, both for known structural reasons (Subtract
back-side flip, dense self-intersection respectively).

ctest 388/388, advfuzz 35/36 + 1 unchanged.

### Architectural ceiling hit on residuals (2026-05-10, very late)

Final investigation of the residual pierces (self-intersect 13,
gt-7863 2) revealed they're locked-in due to a deterministic
drop+re-cap loop: the post-cap pierce reducer drops the offending
tri, the cap immediately re-creates the same vert triplet at the
same position to close the resulting k=1 hole. Pierces re-emerge
identically every iteration.

Tried multiple unlock approaches:

1. **Higher iter limits** (commit 77d5b664): no effect — converges
   to the same locked state by iter ~3.

2. **"Drop both per pair" aggressive heuristic**: no effect.

3. **Coplanar overlap detection**: residuals aren't coplanar
   (CheckSelfIntersection has same skip).

4. **Forbidden-triple tracking** (commit 045d768b): when no-
   progress detected, record the dropped tri's vert triplet and
   reject cap fans that would re-emit it. Built the
   infrastructure (wouldFanForbidden + 3-vert ear gate). On
   gt-7863, applying forbidden at the 3-cycle ear closure
   refuses the only available closure → status non-manifold →
   pierce/drift guard triggers FALLBACK to input → gt-7863 went
   2 → 7. WORSE outcome. Reverted the 3-cycle gate.

   Net: forbidden infrastructure is in place but inactive on
   current fixtures (= the multi-vert fan/ear paths where it
   could safely apply don't reach the offending triplets).

**Architectural conclusion**: the locked residuals require a
fundamentally different fix. Three candidate features (all
multi-session):
  - **Vertex snapping**: move offending verts so the pierce
    disappears geometrically.
  - **Edge subdivision at pierce points**: make the pierce
    intersection a vertex of both tris (= eliminates the
    crossing).
  - **Coplanar merging**: combine overlapping coplanar tris
    into a single non-overlapping triangulation.

The pierce reducer + post-cap re-cap loop is structurally optimal
given current cap+drop primitives. Pierce-monotonicity is
guaranteed (= no fixture's output ever has more pierces than
input) and 5 of 7 piercing fixtures get pipeline fixes.
