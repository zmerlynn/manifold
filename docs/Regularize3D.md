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

What it is NOT. It is not the Boolean. It never FUSES separate objects, and the contract
is UNIFORM: cross-component interaction - overlapping, touching, coplanar OR transversal -
is NEVER this operator's job. Two connectivity components are each regularized on their
own and composed back by concatenation, never unioned, whatever their spatial relation.
Fusion is the Boolean's job, already done upstream in any operation chain; what arrives at
RemoveOverlaps3D is one or a few epsilon-valid manifolds whose only unresolved defect is
INTERNAL (within-component) self-overlap.

Cross-component adjudication (owner decision). A global {w_S >= 1} read over the whole
soup WOULD fuse two components that face-touch or coplanar-overlap into their union - that
is a real semantic choice, and it is deliberately NOT taken here: the per-component read
does not fuse. A caller that wants fusion runs the Boolean. Concretely, a buried plug (a
box whose coincident cap doubles the cover of a second box) is TWO connectivity components
with a cross-component coplanar overlap; this operator passes both through unchanged. The
coplanar fold reaches a coplanar self-overlap only when it is INTERNAL to one connected
component (a doubled internal wall / folded-flat flap stitched into the shell), which the
per-component gate below detects.

Precision contract. Output COORDINATES may be eps-noisy (constructed intersection
points round to double); the output TOPOLOGY is exact - decided from input data, not
from rounded positions. Epsilon-validity is the standing library contract and remains
the output guarantee; strict self-intersection-freedom is re-checked once (below) since
rounding can create eps-scale self-crossings.


## The pipeline

1. DECOMPOSE by connectivity. Split the input into connected components (the existing
   Decompose primitive). The connectivity component is the unit of work and the unit of
   scope - PERIOD; there is no cross-component merge (non-fusion, above).

2. PER-COMPONENT GATE. Each component is tested: valid (`IsManifold` &&
   `Is2Manifold`, overlap3.cpp:1189), non-self-intersecting
   (`Manifold::Impl::IsSelfIntersecting`, properties.cpp:138 - a Morton/AABB broadphase
   over the collider plus a triangle-triangle distance test, shares-vertex skip,
   2*eps relaxation), AND free of WITHIN-component coplanar overlap
   (`DetectCoplanarClusters`, run on this component only, which `IsSelfIntersecting` does
   not flag - R2(i)). A component whose OWN faces coplanar-overlap routes DIRTY even when
   self-intersection-free, so the exact-coplanar fold reaches it; a CROSS-component
   coplanar overlap is invisible to this gate by design (non-fusion).

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
kernel-free on the certified path, reusing boolean3's discipline.  AMENDMENT (stage 6,
owner contract): a micro exact tie-test (Orient3DExactSign) plus the single-global SoS
cascade were RELUCTANTLY ACCEPTED as net-new exact-kernel surface, confined strictly to
the filter-0 fallback - the certified fast path never touches them.  Both are ONE
integer implementation (the whole cascade evaluates as signed sums of products of the
coordinate mantissas over an adaptive-width two's-complement accumulator); the earlier
Shewchuk-expansion substrate and the fixed int256 window were removed in favour of it.
The filter itself uses only Shewchuk's error-bound CONSTANT (per his analysis), no
kernel of his in the build.

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


## Standing adjudication (how refusals are treated)

A fail-closed on an unbuilt axis is NOT a defect to paper over - it is a COMPLETENESS
COUNTEREXAMPLE queued as correctness work. The rule: every solution must be a DECISION-
COMPLETION (extend the exact decision procedure to cover the case), never a repair pass on
rounded geometry, and never a lossy fallback. The queued completions:

- Stage 5 - NEAR-COPLANAR widen + guard: LANDED (reg3d-s5).  The thin band whose
  coplanarity gap is above the filter's error bound but below eps is neither clustered by
  the exact fold nor safely transversal; enumerated transversally its sub-eps-thin cell
  double-rounds to a sliver (the "unresolvable sheet contact" emission).  It is now closed
  by INPUT-SIDE PLANARIZATION in front of candidate B (research memo candidate (a)): a
  pre-pass clusters bbox-overlapping, non-adjacent, 2D-overlapping faces whose LOCAL
  coplanarity gap (the smaller of the two directional vertex-plane maxes - the larger is
  diameter-amplified, a red herring) is below eps; a GLOBAL-PLANARITY GUARD fits one plane
  per cluster (centroid + area-weighted normal) and FAILS CLOSED, distinctly named, when any
  member vertex deviates beyond eps (the anti-chain-reaction net that catches a curved
  near-tangent tessellation); admitted clusters SNAP onto the fitted plane (a <= eps input
  perturbation inside the standing epsilon-valid contract).  B then RE-DERIVES the whole
  arrangement / coupled winding / cell-classify / emit from the snapped, now
  exactly-coplanar input, so the landed exact fold handles it verbatim and the m ==
  winding-jump self-check holds by the exact argument.  This is coordinated-by-construction,
  NOT an emission-time snap (those fight decisions the arrangement already made,
  ExactArrangement3D variant-E kill) - it perturbs the INPUT before enumeration.  The exact
  path is bitwise-unperturbed (a cluster with no near-band pair, and any input with no near
  cluster, is left untouched).  MEASURED: a same-oriented mult-2 buried plug and an
  anti-oriented cancellation with sub-eps-tilted caps resolve oracle-true (GWN membership +
  volume band + tol-invariance; mutation-verified: disabling the snap reverts to the sliver
  fail-closed, disabling the guard makes a curved chain WRONGLY fold to one plane).  A pure
  near-coplanar overlap (deviation < 2eps) is invisible to IsSelfIntersecting (its 2*eps
  normal-nudge separates two near-coplanar faces), so - like the exact fold - the widen is
  reached when B already runs on a dirty component; the corpus's real near-coplanar geometry
  (the hull's body/mask facets) is CROSS-component and stays pass-through under non-fusion.
- Stage 6 - SINGLE-CONVENTION SoS: LANDED (reg3d-s6, reworked reg3d-s6r). One global
  symbolic-perturbation convention (the Shadows pattern lifted to orient3d = 0:
  Edelsbrunner-Mucke perturbation keyed by global vertex index; local-rank reduction
  proven sign-invariant, so per-predicate evaluation decides consistently with ONE global
  perturbation) now DECIDES every genuine non-coplanar transversal exact-zero tie
  (vertex-on-face, edge-on-edge) instead of refusing.  The whole cascade evaluates on ONE
  integer path: each e-coefficient is a signed sum of products of the coordinate
  mantissas (a minor of degree <= 3 over the same windowed coords), accumulated in an
  ADAPTIVE-WIDTH two's-complement integer sized to the worst-case finite-double exponent
  spread.  So the exact orient3d sign is TOTAL: there is no window-fail refusal branch,
  and an exact-zero result means a genuine geometric TIE everywhere, never "uncertain" -
  the ambiguity gate is gone BY CONSTRUCTION (an assertability goal met structurally, not
  by a runtime assert).  The earlier Shewchuk-expansion substrate and the fixed int256
  window were deleted for this one implementation (bitwise-stable on the mesh domain,
  and strictly more total on wide-magnitude inputs where the expansion overflowed
  double).  The coplanar family stays the fold's - coplanar pairs are never SoS-perturbed
  (perturbing them manufactures the sub-eps sliver cells of the thin-cell axis).  Level-0
  discipline holds: the tie is detected exactly (filter, then the exact tie-test
  Orient3DExactSign, RELUCTANTLY-ACCEPTED net-new exact-kernel surface, owner contract,
  confined to filter-0 call sites), then broken by the convention; once-only
  constructions unchanged.  MEASURED: the genus-handle bridged carrier
  (BridgedCaps) resolves oracle-true through the real entry (GWN membership + volume +
  tol-invariance; mutation-verified: disabling the convention reverts to the old
  fail-closed, flipping its direction stays oracle-true).  The former SoS-gate refusals
  NARROW to three precisely-named residues, each still a hard fail-closed with no output:
  (i) coplanar-DOMINATED soups (GT7863, PokedCube) pass the tie gate and fail at EMISSION
  (touching-sheet / sub-eps slivers - the stage-7 thin-cell axis); (ii) COPLANAR/TRANSVERSAL
  ENTANGLEMENT - a transversal wall-wall seam whose endpoint lands ON a coplanar cap cluster
  plane was truncated (nPts==1) by the cap-cluster suppression and failed closed; this was
  RE-DIAGNOSED as NOT a >2-sheet triple point (entangled bars has zero book-of-pages lines,
  every arrangement edge is exactly two walls; the coincident caps are the fold's, not a
  radial sheet - reg3d-radial) and is now CLOSED for that 2-sheet cap-plane case by the
  SEAM-ENDPOINT JUNCTION COMPLETION (reg3d-ent): RecordSeams gives the dropped endpoint the
  fold's in-plane arrangement identity (records the cap-plane crossing, dedups the symmetric
  double-pierce), so the seamed wall drops its interior span and the folded cap + seamed
  walls weld shut at the reentrant corner.  MEASURED oracle-true: EntangledBars (axis-aligned,
  exact union vol) AND a rotated variant (vol == the library boolean union, closure via the
  weld on irrational junctions) - GWN membership + one solid + tol-invariance; mutation-
  verified (disabling the recovery OR the dedup reverts to the truncation fail-closed).  No
  new predicate / no radial rule.  openscad NARROWS strictly (RecordSeams truncations 496 ->
  344) but stays fail-closed - its residue is the near-coplanar-sliver (GT7863-class) and the
  cap-INTERIOR endpoint (a seam piercing a cap face interior, which would need the pierce
  injected as new fold input), both separate axes; (iii) axis-aligned integer geometry can
  graze every winding probe
  seed - NARROWED (reg3d-arr): the CLEAN-face classify now re-probes other interior
  points of the same uncrossed triangle (a constant winding cell above it), so BarsCrossZ
  resolves oracle-true; only the SEAMED path's per-cell centroid probe can still graze
  (the residual seed-policy open).
- Stage 7 - THIN-CELL representability at rounding. A cell thinner than eps has no
  representable double boundary; decide its retention exactly from input data, do not emit a
  sub-eps sliver.

Within-component note. An internal coplanar overlap (two coincident sheets, no shared
edges) is one connected 2-manifold ONLY via a genus handle (a solid bridge, or a hole
through the coincident caps): sharing the coincidence boundary is a non-manifold pinch or a
transversal crossing, and a buried plug's three winding levels (0/1/2) cannot be bridged
consistently (a handle would weld a w0|w1 sheet to a w1|w2 sheet). Only a stacked pair (both
interiors w=1) is bridgeable, and its junction is on the stage-6 SoS axis - now landed:
the bridged carrier is the family's oracle-true resolve pin.


## Parallelism (order-freeness, measured)

The operator is order-free: the output is bit-identical regardless of component / face /
query order (the compose is concatenation, the winding is a coupled integer read, the
predicates are level-0 on input coords). So it is embarrassingly parallel - per-component,
per-face-pair enumeration, and per-query winding are all independent. A BVH-accelerated
winding query plus threading over components/queries is a named LATER perf pass, not a
correctness change (the sequential result is the spec the parallel one must match bitwise).


## Open list (honest)

- THE BUILD. The cell complex + halfedge-boundary extraction is the largest unbuilt
  piece; the fragment did point classification of recorded cells, not the boundary
  build. B is fragment-validated, not landed.
- SINGLE GLOBAL SoS: LANDED (stage 6 above) - the bridge-junction family resolves
  oracle-true.  What remains open is the NARROWED residue behind it: the stage-7
  thin-cell emission (coplanar-dominated soups) and the component-local seed policy
  (both named, hard fail-closed).  The former "entangled bars / >2-sheet triple
  point" residue was RE-DIAGNOSED as a 2-sheet cap-plane seam TRUNCATION (not a
  triple point) and CLOSED by the seam-endpoint junction completion (reg3d-ent,
  standing adjudication (ii)); a genuine >2-sheet radial junction is not forced by
  any corpus carrier (reg3d-radial verified the reduction but found zero
  book-of-pages lines), so the radial rule stays unbuilt (kernel-tripwire-gated).
- EXACT-KERNEL SURFACE, QUEUED FOR REVISIT (owner contract, reluctant acceptance).
  Whether the arrangement can be structured to avoid needing an exact orient3d kernel
  at all is still open; the current integer Orient3DExactSign / SoS cascade is the
  reluctantly-accepted answer, not a settled one.  TRIPWIRE: this integer path is
  legitimate ONLY as ONE predicate at ONE call site (~100 lines, exhaustively testable).
  If a SECOND exact predicate or a SECOND call site is ever needed, VENDOR Shewchuk's
  public-domain predicates.c instead of growing this - do NOT rebuild expansion
  arithmetic piecemeal.
- NEGATIVE WINDING / subtraction, untested. openscad's soup winding reaches -1; the
  {w_S>=1} threshold read should absorb it, but no subtraction carrier has exercised
  it.
- COMPONENT-LOCAL SEED policy. One ray per connected cell-component with a direction
  policy that avoids the far-seed near-grazing measured on siA; micro cost.  PARTIALLY
  CLOSED (reg3d-arr): the CLEAN-face classify dodges the graze by re-sampling the
  constant winding cell above an uncrossed triangle at several interior points (sound,
  no uniformity assumption) - closes BarsCrossZ.  The SEAMED path's per-cell centroid
  probe is the remaining graze site.
- CLEAN-FACE classification is PER-FACE (reg3d-s7b/arr), not a per-patch representative
  flood.  The flood assumed uniform coverage per clean-clean-connected patch; a severely
  folded soup (PokedCube's everted corner) breaks that (a shares-vertex-skip crossing
  puts differently-wound faces in one patch), so the representative wrongly dropped
  genuine boundary faces = a latent silent-wrongness class.  Per-face directly measures
  each face's own winding.  Bitwise-identical to the flood on uniform patches (siA/siB).
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
    connectivity, and that split IS the scope - uniformly (the owner adjudicated against a
    coplanar-fusion exception). The residual tension is only that connectivity, not
    geometry, draws the line: two solids sharing an edge/vertex are ONE component (so their
    overlap is within-scope) while two coincident-but-disjoint solids are TWO (out of
    scope, pass-through). The operator does not second-guess this - a caller wanting a
    different grouping (fusion) runs the Boolean first. What connectivity cannot deliver is
    a within-component coplanar RESOLVE: the connection that makes such a defect one
    component is on the stage-6 SoS axis, so it fails closed (never wrong) until stage 6.

R4. REPO-TERM MISUSE. This doc leans on Canonicalize / IsSelfIntersecting / the
    collider / Shadows / Winding03 as reusable primitives; if any is described with the
    wrong contract (e.g. treating IsSelfIntersecting's epsilon relaxation as exact, or
    calling the signed-multiplicity fold an operand label it is not) the reuse story is
    thinner than stated. eps-vs-tol and halfedge-not-DCEL are the standing traps.

R5. THE DESIGN READS AS MORE SETTLED THAN THE CODE. Everything past the gate is
    validated at fragment scale with exact-rational probes, not built. The unbuilt
    axes (SoS, subtraction, the halfedge boundary build, the seed policy) are named
    opens, and a reader should not mistake "fragment-validated" for "landed".

KNOWN OUTPUT ARTIFACT (ent-verify, disclosure): on rotated/irrational junctions
the emission weld can materialize a single inert zero-area triangle (paired
edges, topology sound, oracle-exact volume; IsSelfIntersecting-blind).  The
clean fix is the once-only fold-vertex identity hardening (queued, simplicity
pass); the R1 weld blind spot remains the underlying open.
