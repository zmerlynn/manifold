# RemoveOverlaps3D as a regularization operator (current design)

STATUS: DESIGN, post-owner-decision. This is the small, current statement of what
RemoveOverlaps3D should be, distilled from the campaign in docs/ExactArrangement3D.md
(kept as the historical record - read this doc for the design, that doc for how it was
reached). Candidate B below is validated at FRAGMENT scale (exact-rational probes with
a double mirror), not landed production C++: `RemoveOverlaps3D`
(src/overlap3.cpp:1293) today still runs the v3 sweep pipeline (SweepEmit / BuildSlabs
/ EmitCaps / EmitStrips). This doc specifies the operator that replaces the sweep as
the resolver; the sweep stays as reference and fixture harness until B lands.

Conventions: magnitudes only (exact figures cite the lab notebooks), ASCII, `eps`
means the machine-scale weld radius (kPrecision * bBox.Scale(), utils.h), `tol` means
a user tolerance. w_S is the signed winding of the oriented input face soup.


## Contract and semantics

RemoveOverlaps3D is the REGULARIZATION OPERATOR the library never had: the map from
any valid oriented face soup to the boundary of the solid {p : w_S(p) >= 1}. It takes
a `Manifold::Impl` that may be self-overlapping (a valid oriented 2-manifold whose
interior is multiply covered) and returns the boundary of the region covered at least
once - the boundary-of-a-simple-solid reading downstream consumers assume.

What it is NOT. It is not the Boolean. It never FUSES separate objects: two disjoint
components that happen to touch are each regularized on their own and composed back
by concatenation, never unioned. Fusion is the Boolean's job, already done upstream in
any operation chain; what arrives at RemoveOverlaps3D is one or a few epsilon-valid
manifolds whose only unresolved defect is INTERNAL self-overlap. (Whole-soup {w_S>=1}
would fuse; per-component {w_S>=1} does not - the deliberate choice, matching the
touching-contacts posture.)

The one scoped exception is a genuine COPLANAR self-overlap that connectivity would
split. A buried plug (a box whose coincident cap doubles the cover) is two connectivity
components with a CROSS-component coplanar overlap, and connecting them into one
manifold instead forces a transversal entanglement (measured) - so a per-component gate
can never see the defect. The decompose step therefore UNITES the components that share
an exactly-coplanar 2D-area overlap cluster (DetectCoplanarClusters' overlap2D witness),
routing the united super-component to the coplanar fold. This is scoped to a real
coplanar self-overlap: touching contacts (edge/vertex, zero area) and disjoint objects
never cluster, so the non-fusion posture holds for everything else. (Owner-review note:
this widens the per-component boundary for the coplanar-overlap class only; it is the
change that makes the fold reachable in production instead of test-hook-only.)

Precision contract. Output COORDINATES may be eps-noisy (constructed intersection
points round to double); the output TOPOLOGY is exact - decided from input data, not
from rounded positions. Epsilon-validity is the standing library contract and remains
the output guarantee; strict self-intersection-freedom is re-checked once (below) since
rounding can create eps-scale self-crossings.


## The pipeline

1. DECOMPOSE by connectivity, then MERGE coplanar-overlapping components. Split the
   input into connected components (the existing Decompose primitive), and when there
   is more than one, unite the components that share an exactly-coplanar 2D-area
   overlap cluster (see "What it is NOT"). The (possibly merged) component is the unit
   of work and the unit of scope.

2. PER-COMPONENT GATE. Each component is tested: valid (`IsManifold` &&
   `Is2Manifold`, overlap3.cpp:1189), non-self-intersecting
   (`Manifold::Impl::IsSelfIntersecting`, properties.cpp:138 - a Morton/AABB broadphase
   over the collider plus a triangle-triangle distance test, shares-vertex skip,
   2*eps relaxation), AND free of coplanar overlap (`DetectCoplanarClusters`, which
   `IsSelfIntersecting` does not flag - R2(i)). A component with a coplanar overlap
   routes DIRTY even when self-intersection-free, so the exact-coplanar fold reaches
   it.

3. EARLY-EXIT clean components. A component that passes the gate is already the
   boundary of a simple solid; it is copied through untouched.

4. RUN CANDIDATE B per DIRTY component. A component that fails the self-intersection
   test is regularized by B (next section): compute its local arrangement + per-cell
   w_S, emit the oriented {w_S>=1} boundary.

5. RE-GATE B's output once. B's emitted geometry is double-rounded, so it gets the
   SAME gate as the input (validity + IsSelfIntersecting). A clean pass composes in; a
   failure is the honest fail-closed outcome (a recorded FatalReason, overlap3.h:49),
   never a silent wrong result.

6. COMPOSE BACK. Concatenate the early-exit components and the regularized components.
   No cross-component weld, no fusion.

This is candidate A reduced to its GATE. The owner dropped candidate A's Boolean fold
(the pair corpus decomposes into separate valid solids only because the fixtures were
built by composition; a production chain has already done any wanted union upstream).
What survives of A is decompose + test + early-exit + route-to-B + compose.


## B's mechanism (the dirty-core resolver)

B computes the arrangement and winding of ONE self-overlapping component in DOUBLES,
kernel-free, reusing boolean3's discipline rather than introducing an exact kernel.

- ENUMERATION. The collider's broadphase returns candidate face pairs; keep genuine
  non-adjacent crossings, SKIP self-adjacent (shared-vertex) pairs. Every crossing
  DECISION (face-pair straddle, edge-plane crossing, edge-edge z-order) is a level-0
  orient2d/orient3d on INPUT coordinates through a static Shewchuk error filter; the
  exact-check fallback fires only when the filtered sign is uncertain. Constructed
  intersection POSITIONS (output geometry) are built once each, in double.

- WINDING. Per-cell w_S by signed INTEGER face-crossing deltas: crossing an oriented
  face f along a path changes w by sign(dot(t, n_f)) on the input normal, coupled
  across the arrangement (unite non-crossing edges, one seed ray per connected
  cell-component, flood - the Winding03 discipline, boolean3.cpp:388). The deltas are
  +-1 integers, never near-zero, so the winding half is FP-safe by construction.
  The cost the localizer's dissolution ("the component is the unit of scope") hid:
  each seed cast is an O(ntri) winding query scoped to the component, so ONE large
  component pays the whole-component cost - measured at seconds-per-query scale on the
  largest corpus component (GT7081), with no adaptive orient3d kernel in the tree
  today, so this is net-new surface area (.claude/lane-reports/
  v5-verify-probe-1783913337.md).

- ASSEMBLY. Build the arrangement cell complex; radial-order the halfedges around each
  arrangement edge; classify each cell by its integer w_S; emit the oriented boundary
  of {w_S >= 1}. The radial order is decided by an input-predicate determinant. On the
  corpus every arrangement edge is exactly two distinct planes (a two-sheet edge, the
  PairUp shape), so the radial sort is a trivial 2-element sort; the >2-sheet branch
  carries the specified determinant rule but is never forced.

Why B replaces PairUp rather than inheriting it: PairUp pairs an edge's crossings by
1D start-end alternation, which is genuinely violated on the doubly-covered (w_S=2)
stratum. B does not pair - it thresholds the coupled integer w_S, which is indifferent
to crossing order, so the double-cover and the negative-winding regime fall out of the
same read.

The one global symbolic-perturbation convention (the `Shadows(p,q,dir)` pattern,
shared.h, lifted to orient2d/orient3d = 0) governs EVERY predicate family, so a ray or
edge grazing a shared boundary resolves identically for every probe. Each intersection
/ triple point is constructed ONCE via a single canonical derivation and referenced
everywhere; two derivations of one point diverge without bound in the near-parallel
tail, so once-only is a required design element, not advice.

What B reuses from v3, honestly: Canonicalize (overlap3.cpp:218 - input quantization +
the signed-multiplicity fold; the winding delta sign is exactly CanonicalFace.mult
times the crossed face orientation), the collider broadphase (also the gate's), the
gate itself, the Shadows tie-break pattern, and the Winding03 coupling. B does NO plane
sweep: it does not touch BuildSlabs / SlabResolver / track-extension / EmitCaps /
EmitStrips / FindSeams. The sweep's dense-critical ArrangementBudget refusal is a
sweep artifact (the O(seam^2) x-crossing density is a sweep-PROJECTION effect, not
radial structure) and does not transfer to B.


## What is proven (evidence by pointer - one line each, not re-narrated)

- The correct topology is a pure function of the oriented soup; no ambiguity witness
  exists: .claude/lane-reports/a0-verify-witness-1783911000.md (V4ImplPlan.md, A0).
- Every eps-quantized identity channel is degenerate on the hard carriers, so the
  barrier is PRECISION, not information: V4ImplPlan.md stage A0 + a0-verify-prefold.
- The lost material is real: exact w_S = 2 at the Havoc lump, reproduced bit-for-bit by
  two independent algorithms: .claude/lane-reports/v5-verify-probe-1783913337.md.
- B's deciding predicates are FP-safe exhaustively (static filter proven sound over
  ~1.6B predicates, zero certified flips): .claude/lane-reports/v5b-s12-1783918690.md.
- Every combinatorial decision restructures to a level-0 input predicate (constructed-
  point precision dissolved for decisions): same lane (S2).
- Zero genuine multi-sheet edges on the corpus (the triple-point mechanism gap
  dissolves): .claude/lane-reports/v5b-s34-1783917823.md (S3).
- The coupled winding half paper-executes against Havoc's exact ground truth: same
  lane (S4).
- The w=2 alternation residual is absorbed by the threshold read (measured non-
  alternating SSEE): .claude/lane-reports/v5b-r3-1783919852.md (R3-i).
- The first built fragment runs end-to-end on Havoc AND a true single-shell self-
  intersector (siA), matching exact ground truth: same lane (R3-ii).
- The filtered enumeration is kernel-free end-to-end (the edge-edge z-order 100%
  certified, fallback fires only on genuine exact ties): v5b-r4-1783923225.md.
- Candidate A resolves the whole decomposable corpus (topology); the fold is dropped:
  .claude/lane-reports/v5b-fold-1783918998.md + ExactArrangement3D.md kernel section.
- Zero silent wrong-resolves on the union corpus (fail-closed is the only failure
  mode observed): ExactArrangement3D.md "Zero silent wrong-resolves".


## Open list (honest)

- THE BUILD. The cell complex + halfedge-boundary extraction is the largest unbuilt
  piece; the fragment did point classification of recorded cells, not the boundary
  build. B is fragment-validated, not landed.
- SINGLE GLOBAL SoS, unexercised. The cross-operand vertex-on-face tie family that
  GT7863 (coplanar exact-zeros) and openscad (coincident verts) need is specified but
  UNEXERCISED at fragment scale; the fragment hit only same-operand incidences.
- NEGATIVE WINDING / subtraction, untested. openscad's soup winding reaches -1; the
  {w_S>=1} threshold read should absorb it, but no subtraction carrier has exercised
  it.
- COMPONENT-LOCAL SEED policy. One ray per connected cell-component with a direction
  policy that avoids the far-seed near-grazing measured on siA; micro cost, unbuilt.
- ONCE-ONLY under a genuine triple point. The corpus forces zero triple points, so the
  once-only construction rule is unstressed; a synthetic co-axial-3-face carrier is the
  attack that would exercise it.


## Risks (attack surface for the adversarial round)

R1. THE SELECTIVE WELD MAY HAVE BEEN DISSOLVED WRONGLY. B replaces a whole component
    and composes back, so the patch-to-exterior stitch is gone - but B's own
    constructed intersection coordinates, rounded to double, can land sub-eps-distinct.
    If the ordinary assembly weld (the global eps grid) merges two exact-distinct B
    verts it re-manufactures the twin the arrangement just resolved; if it merges a
    B-interior vert onto a boundary vert it corrupts B's topology. The re-gate does NOT
    backstop this failure: when the weld merges two exact-distinct B verts, the
    triangles straddling the merge now share that vertex position, so
    IsSelfIntersecting's shares-vertex skip drops the pair - a self-FOLD manufactured by
    the merge is invisible to it. Only the validity half (IsManifold/Is2Manifold)
    catches a weld outcome, and only the non-manifold one (a pinch/tear), never a fold.
    Nor is a constant-radius bounded weld the escape: the wall-A weld-bump probe killed
    it at every multiplier (too small re-manufactures the twin, too big collapses
    slivers to holes). R1's honest rebuttal is B's own STRUCTURAL defense, not the gate:
    the once-only construction rule (each intersection built once, referenced
    everywhere) never makes two rounded images of one point, and radial assembly (not
    projection, not self-location) has no O(seam^2) projected twins. The residual that
    survives - two GENUINELY DISTINCT arrangement points rounding within eps and being
    merged by the uniform weld - is open, neither cheaply weldable nor reliably
    re-gated.

R2. THE GATE'S SELF-COLLISION TEST MAY BE MORE EXPENSIVE / DIFFERENT THAN CLAIMED.
    IsSelfIntersecting is a broadphase + tri-tri distance over every component, run
    once per input component AND again on every B output. At GT7081 scale (tens of
    thousands of tris) "cheap" needs measuring. Worse, the gate is systematically
    CLEAN-biased: the 2*eps shares-vertex relaxation SUPPRESSES near-miss detection (it
    returns non-intersecting when an eps normal nudge separates the pair) - it does not
    flag near-misses - and that bias does not agree with B's arrangement notion. The
    silent-miss direction is GATE-CLEAN but B-DIRTY: a genuine crossing whose two verts
    sit in a near-degenerate band just outside the eps weld (distinct enough not to
    merge, close enough to trip the 2*eps skip) passes the gate, EARLY-EXITS as clean,
    and carries an unregularized self-overlap through silently - worse than fail-closed,
    though narrow (it needs a near-degenerate config in that thin band). The reverse
    (gate-dirty, B finds nothing) is near-empty: B's shared-vertex skip set is a subset
    of the gate's distance skip set, so nothing the gate robustly flags is skipped by B;
    the only case is a legitimate touching contact the gate flags and B re-emits clean -
    wasteful, not wrong. R2(i) (gate-clean/B-dirty) and R1 are the SAME blind spot: the
    weld-merge fold in R1 manufactures exactly the shared-vertex config the gate skips,
    so the re-gate that would "catch" R1 IS the clean-biased gate R2 flags.

R3. PER-COMPONENT SEMANTICS VS THE CONNECTIVITY BOUNDARY. Decompose splits by halfedge
    connectivity. Two solids touching at a shared edge or vertex are ONE connected
    component, so B would fuse their overlap where whole-soup semantics is what the
    caller wanted - or vice versa, a caller wanting no fusion is served only if the
    Boolean already separated them upstream. Connectivity may be the wrong
    decomposition boundary at exactly the touching-contact cases the operator claims to
    respect.

R4. REPO-TERM MISUSE. This doc leans on Canonicalize / IsSelfIntersecting / the
    collider / Shadows / Winding03 as reusable primitives; if any is described with the
    wrong contract (e.g. treating IsSelfIntersecting's epsilon relaxation as exact, or
    calling the signed-multiplicity fold an operand label it is not) the reuse story is
    thinner than stated. eps-vs-tol and halfedge-not-DCEL are the standing traps.

R5. THE DESIGN READS AS MORE SETTLED THAN THE CODE. Everything past the gate is
    validated at fragment scale with exact-rational probes, not built. The unbuilt
    axes (SoS, subtraction, the halfedge boundary build, the seed policy) are named
    opens, and a reader should not mistake "fragment-validated" for "landed".
