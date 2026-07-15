# RemoveOverlaps3D as a regularization operator (current design)

STATUS: SPEC (landed). This is the statement of what RemoveOverlaps3D IS, distilled
from the campaign in docs/ExactArrangement3D.md (kept as the historical record - read
this doc for the design, that doc for how it was reached). `RemoveOverlaps3D`
(src/overlap3.cpp) IS this operator: the v3 sweep pipeline it replaced (SweepEmit /
BuildSlabs / EmitCaps / EmitStrips) has been DELETED from the branch - its history
lives in git and on the explore/sweep-plane-3d-v3/-v4 branches. The dirty-component
resolver (enumeration + coupled winding + the {w_S>=1} halfedge-boundary emission) is
built and runs on the corpus; the honest open axes (a >2-sheet radial junction, the
thin-cell emission, negative-winding subtraction, the component-local seed policy)
fail closed with a named reason and are listed in the open list below.

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
   `Is2Manifold`, overlap3.cpp), non-self-intersecting
   (`Manifold::Impl::IsSelfIntersecting`, properties.cpp - a Morton/AABB broadphase
   over the collider plus a triangle-triangle distance test, shares-vertex skip,
   2*eps relaxation), AND free of WITHIN-component coplanar overlap
   (`DetectCoplanarClusters`, run on this component only, which `IsSelfIntersecting` does
   not flag - R2(i)). A component whose OWN faces coplanar-overlap routes DIRTY even when
   self-intersection-free, so the exact-coplanar fold reaches it; a CROSS-component
   coplanar overlap is invisible to this gate by design (non-fusion).

3. EARLY-EXIT clean components. A component that passes the gate is already the
   boundary of a simple solid; it is copied through untouched.

4. RUN THE RESOLVER per DIRTY component. A component that fails the self-intersection
   test is regularized by the resolver (next section): compute its local arrangement + per-cell
   w_S, emit the oriented {w_S>=1} boundary.

5. RE-GATE the resolver's output once. The resolver's emitted geometry is double-rounded, so it gets the
   SAME gate as the input (validity + IsSelfIntersecting). A clean pass composes in; a
   failure is the honest fail-closed outcome (a recorded FatalReason, overlap3.h),
   never a silent wrong result.

6. COMPOSE BACK. Concatenate the early-exit components and the regularized components.
   No cross-component weld, no fusion.

This is candidate A reduced to its GATE (candidate A was an earlier design with a
Boolean fold; the owner dropped the fold - the pair corpus decomposes into separate
valid solids only because the fixtures were built by composition, and a production
chain has already done any wanted union upstream). What survives is decompose + test +
early-exit + route-to-the-resolver + compose.


## Phase map (code)

The operator's phases and where they live in src/overlap3.cpp (the sole impl):

- DECOMPOSE  -> DecomposeComponents (connectivity split; the unit of scope).
- GATE       -> GateComponent (validity + IsSelfIntersecting + DetectCoplanarClusters).
- DISPATCH   -> RemoveOverlaps3D (public entry: gate every component, early-exit clean,
                route dirty, re-gate, compose).
- RESOLVE    -> ResolveComponent (the dirty-core resolver, next section), which runs:
    - PLANARIZE -> SnapNearCoplanarClusters (near-coplanar widen + global-planarity guard).
    - ENUMERATE -> RecordSeams (+ EnumerateSelfCrossings / EdgePiercesTri, sos::* kernel).
    - WIND      -> WindingAt / RobustWinding (coupled integer-delta w_S).
    - EMIT      -> EmitComponentBoundary: FoldCoplanarClusters, EmitSeamedFace,
                   EmitCleanFaces, then the assembly (BuildImpl + SplitTouchingSheets).
- RE-GATE    -> a second GateComponent on the resolver's double-rounded output.
- COMPOSE    -> ComposeComponents (concatenation, no fusion).


## The resolver's mechanism (the dirty-core resolver)

The resolver computes the arrangement and winding of ONE self-overlapping component in DOUBLES,
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
  cell-component, flood - the Winding03 discipline, boolean3.cpp). The deltas are
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

Why the resolver replaces PairUp rather than inheriting it: PairUp pairs an edge's crossings by
1D start-end alternation, which is genuinely violated on the doubly-covered (w_S=2)
stratum. The resolver does not pair - it thresholds the coupled integer w_S, which is indifferent
to crossing order, so the double-cover and the negative-winding regime fall out of the
same read.

The one global symbolic-perturbation convention (the `Shadows(p,q,dir)` pattern,
shared.h, lifted to orient2d/orient3d = 0) governs EVERY predicate family, so a ray or
edge grazing a shared boundary resolves identically for every probe. Each intersection
/ triple point is constructed ONCE via a single canonical derivation and referenced
everywhere; two derivations of one point diverge without bound in the near-parallel
tail, so once-only is a required design element, not advice.

What the resolver reuses, honestly: the collider broadphase (also the gate's), the
gate itself (IsSelfIntersecting + the coplanar-cluster detector), the Shadows
tie-break pattern, the Winding03 coupling, and the assembly (BuildImpl +
SplitTouchingSheets) that welds emitted triangles into a 2-manifold and re-separates
touching sheets.  It does its own input read (a TriSoup over the component's
halfedges), NOT the deleted sweep's Canonicalize quantization.  It does NO plane
sweep - there is no BuildSlabs / SlabResolver / EmitCaps / EmitStrips / FindSeams
(all deleted).  The v3 sweep's dense-critical ArrangementBudget refusal was a
sweep-PROJECTION artifact (an O(seam^2) x-crossing density, not radial structure) and
has no analog here.


## What is proven (evidence by pointer - one line each, not re-narrated)

- The correct topology is a pure function of the oriented soup; no ambiguity witness
  exists: .claude/lane-reports/a0-verify-witness-1783911000.md (V4ImplPlan.md, A0).
- Every eps-quantized identity channel is degenerate on the hard carriers, so the
  barrier is PRECISION, not information: V4ImplPlan.md stage A0 + a0-verify-prefold.
- The lost material is real: exact w_S = 2 at the Havoc lump, reproduced bit-for-bit by
  two independent algorithms: .claude/lane-reports/v5-verify-probe-1783913337.md.
- the resolver's deciding predicates are FP-safe exhaustively (static filter proven sound over
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
  by INPUT-SIDE PLANARIZATION in front of the resolver (research memo candidate (a)): a
  pre-pass clusters bbox-overlapping, non-adjacent, 2D-overlapping faces whose LOCAL
  coplanarity gap (the smaller of the two directional vertex-plane maxes - the larger is
  diameter-amplified, a red herring) is below eps; a GLOBAL-PLANARITY GUARD fits one plane
  per cluster (centroid + area-weighted normal) and FAILS CLOSED, distinctly named, when any
  member vertex deviates beyond eps (the anti-chain-reaction net that catches a curved
  near-tangent tessellation); admitted clusters SNAP onto the fitted plane (a <= eps input
  perturbation inside the standing epsilon-valid contract).  the resolver then RE-DERIVES the whole
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
  reached when the resolver already runs on a dirty component; the corpus's real near-coplanar geometry
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
  (touching-sheet / sub-eps slivers - the stage-7 thin-cell axis) [SUPERSEDED: GT7863 and
  PokedCube now RESOLVE oracle-true - the emission "wall" was ARRANGEMENT INCOMPLETENESS
  from unsound shares-vertex skips, closed by reg3d-wjump (RecordSeams off-vertex-pierce
  recovery) + reg3d-7863c1 (DetectCoplanarClusters two-pass); see Gap ledger closure plan
  pt 3]; (ii) COPLANAR/TRANSVERSAL
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
  new predicate / no radial rule.  openscad NARROWED strictly with the recovery and
  reg3d-7863c1's re-clustering (magnitudes only), and its F11 seam-truncation wall is now CLOSED
  [reg3d-oscad REOPEN, LANDED: exact rational reconstruction of every truncated seam's TRUE
  second endpoint refuted the earlier "ZERO shares-vertex / no fold-owned completion /
  plane-based-representation PROOF wall" reading as a MEASUREMENT GAP (it searched only for a
  dropped interior pierce and measured shares-vertex by INDEX).  The dominant truncations are
  MEASURE-ZERO contacts - coincident-position DUPLICATE vertices (a shares-vertex family the
  index-keyed skip misses), vertex-on-edge / vertex-on-face T-junctions, and collinear
  edge-on-edge overlaps - whose exact tri-tri intersection is a POINT or a boundary segment, so
  a valid arrangement has NO seam there (proven exactly: none has a clean off-plane edge piercing
  the other's strict interior).  Three RecordSeams decision-completions close it: the
  shares-vertex skip/recovery keyed on POSITION coincidence (recovering the coincident-vertex
  transversal crossings, the shared vertex the second endpoint), a phantom-seam guard (no clean
  off-plane interior pierce -> measure-zero -> skip), and a sub-eps collapse (endpoints dedup'd
  under eps -> witness-theorem collapse).  openscad now fails one wall DEEPER at F4, the seam
  sub-face / winding CLASSIFY probe = the SAME near-tangent winding-probe phenomenon as GT7081];
  still fail-closed, never a silent wrong resolve.  (iii)
  axis-aligned integer geometry can
  graze every winding probe
  seed - NARROWED (reg3d-arr): the CLEAN-face classify now re-probes other interior
  points of the same uncrossed triangle (a constant winding cell above it), so BarsCrossZ
  resolves oracle-true; only the SEAMED path's per-cell centroid probe can still graze
  (the residual seed-policy open).  A SEPARATE winding-probe residue - a filter-PRECISION
  graze where a near-tangent shallow-dihedral face is filter-uncertain but exact-decidably
  off-plane (GT7081's 0.002deg geometry) - is now CLOSED (reg3d-c2bx): WindingAt's filter-0
  terms ESCALATE to Orient3DExactSign (the blessed winding-probe caller, filter-first), deciding the
  graze exactly.  This is distinct from the seed-policy graze (moving the probe cannot escape
  a face near-tangent over its whole extent; deciding it exactly can).  MEASURED on GT7081:
  both dirty shells now clear the winding classify with zero genuine ties, revealing the
  Cluster-1 emission wall underneath (below).
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
per-face-pair enumeration, and per-query winding are all independent.

LANDED (winding-query perf axis): the clean-face winding walked every triangle per query
(the resolve hot loop on single-component shells). It now queries a per-component Morton
triangle collider with each winding ray's SEGMENT AABB - a PROVEN EXACT crossing-superset
(a genuine crossing point lies in both the segment box and the triangle box, and both are
the componentwise min/max of the same doubles the exact predicates read, so the closed-
interval overlap test cannot exclude it; no box inflation, no eps). The seed plane-side
sign is precomputed once per component (pure per triangle for a fixed seed). Both are
bit-identical to the walk on the corpus (validated: winding-value divergence and
crossing-superset violations both zero). LANDED (order-freeness): the clean-face classify
and the per-component gate+resolve run through manifold::for_each_n / autoPolicy as
two-pass (classify in parallel, reduce in ascending index), so the output is bitwise-
identical across thread counts and equal to the sequential build - the determinism the
order-freeness promises, now a standing rail.


## Open list (honest)

- THE BUILD. LANDED. The per-face 2D-arrangement + halfedge {w_S>=1} boundary
  emission (EmitComponentBoundary: FoldCoplanarClusters, EmitSeamedFace,
  EmitCleanFaces, then BuildImpl) is built and runs the corpus (siA/siB resolve
  oracle-true through the real entry). What is NOT built is the >2-sheet radial
  junction branch (no corpus carrier forces it - kernel-tripwire-gated) and the
  thin-cell emission; both fail closed, named, below.
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
  reluctantly-accepted answer, not a settled one.  TRIPWIRE, RELAXED by the owner to
  ONE PREDICATE, ONE IMPLEMENTATION (reg3d-c2bx): the integer path is legitimate as
  ONE exact predicate FORM (this ~100-line adaptive-integer orient3d), and ADDITIONAL
  CALLERS are within the blessed pattern (zero new arithmetic) as long as each stays
  FILTER-FIRST (exact fires only behind a filter 0).  FIRE CENSUS (orient-land
  necessity attack): the exact cascade earns its keep in two roles.  (i) STRUCTURAL-TIE
  DETECTION - the bulk of the tie-cascade fallback: genuine exact-zeros, overwhelmingly
  STRUCTURAL (two boxes sharing an axis coordinate, or a repeated point), where the e^0
  sum is provably zero; these route to the single-global SoS, which breaks them by
  perturbation.  (ii) NEAR-TANGENT SIGN DECISIONS - a filter-uncertain but
  exact-decidably-NONZERO e^0, concentrated almost entirely on the largest corpus
  component (GT7081) and its ~0.002deg twin faces; this is the winding-probe graze the
  exact sign resolves.  CALLER INVENTORY of the exact e^0 predicate (Orient3DExactSign):
  the EdgePiercesTriSoS edge-in-plane guard and the WindCrossTri winding-crossing escalation
  (both the WindingAt winding-probe AND the RecordSeams phantom-seam guard ride it -
  cleanPierce is no longer a distinct exact-call site) - the production callers (plus the test
  probe), each FILTER-FIRST (exact fires only behind a filter 0).  The tie cascade Orient3DSoS is NO LONGER
  among them (orient-land item 1): its SoS e^0 (K==0) monomial group IS the exact
  orient3d - the same real terms over the same accumulator - so a pre-SoS ExactSign
  shortcut returned the identical sign and was provably redundant; it was dropped, a
  bitwise no-op on the corpus carriers.  The >2-sheet radial rule remains a BANKED
  (unbuilt) prospective caller.
  Why these exact calls stay, not filter-only (orient-land item 2, an explicit soundness
  adjudication): the EG-guard exact call is a NEVER-EXERCISED-ON-CORPUS SOUNDNESS
  BACKSTOP, deliberately KEPT.  On every corpus carrier the guard's filter-0 endpoints
  are also exact-0 (the exact check always confirms the filter's "both near plane"), so a
  FILTER-ONLY guard would be corpus-identical.  But off-corpus it is a lossy guess: a
  filter-0 edge whose two endpoints straddle the plane exactly (a near-tangent transversal
  pierce - both endpoints within filter-eps of the plane but on OPPOSITE exact sides)
  genuinely pierces, and the exact call is what routes it to the TOTAL SoS instead of
  suppressing it.  A filter-only guard would return "no seam" for that edge WITHOUT
  routing to SoS - a definite geometric answer on uncertainty, i.e. a silent guess, not a
  fail-closed refusal.  Zero-oracle-wrong is corpus-graded; fail-closed is global; a guard
  that guesses off-corpus is banned.  So the exact call stays.  The phantom-seam guard's
  cleanPierce is the same backstop kind (reg3d-phantom-close): it now RIDES the WindCrossTri
  filter-then-exact chain, completing each filter-refused straddle/interior sign with
  Orient3DExactSign behind a filter 0 - never the SoS convention,
  whose perturbation would manufacture a PHANTOM strict-interior pierce out of a
  measure-zero contact (endpoint on-plane or crossing on the triangle boundary) - so on
  corpus it never flips a verdict (measured: exact fires but zero pierce deltas, bitwise
  no-op) yet off-corpus it catches a near-tangent genuine crossing the filter cannot
  certify and fails closed instead of silently dropping it.
  PRICING the WindingAt caller (the 1 -> 0 the attack priced): it is load-bearing for
  GT7081's resolve ALONE - removing it regresses GT7081 to a fail-closed
  DirtyComponentUnresolved, and is a no-op everywhere else.  Its fires are
  da = Orient3DFilterSign(a, b, c, p), which takes NO seed argument, so they are
  SEED-INVARIANT: a margin-max ray or a probe/path reroute changes only the seed-side
  terms and cannot clear a da graze (refuted by construction, not by measurement).  The
  only genuine route to zero is a coupled-integer WIND re-architecture that removes
  ray-casting entirely, or accepting the loss of GT7081's oracle-true resolve - a research
  axis, not a tweak.
  What stays tripwired is a SECOND predicate FORM: if one is ever needed, VENDOR
  Shewchuk's public-domain predicates.c instead of growing this - do NOT rebuild
  expansion arithmetic piecemeal.
- NEGATIVE WINDING / subtraction, untested. openscad's soup winding reaches -1; the
  {w_S>=1} threshold read should absorb it, but no subtraction carrier has exercised
  it.
- COMPONENT-LOCAL SEED policy. One ray per connected cell-component with a direction
  policy that avoids the far-seed near-grazing measured on siA; micro cost.  PARTIALLY
  CLOSED (reg3d-arr): the CLEAN-face classify dodges the graze by re-sampling the
  constant winding cell above an uncrossed triangle at several interior points (sound,
  no uniformity assumption) - closes BarsCrossZ.  A DISTINCT graze - the winding-probe
  FILTER-PRECISION residue (a near-tangent face filter-uncertain but exact-decidably
  off-plane) - is closed by the WindingAt exact escalation (reg3d-c2bx), which decides
  the graze rather than moving the probe.  What remains open is the pure SEED-POSITION
  graze on the SEAMED per-cell centroid probe (a genuine exact-zero tie moving the
  probe could escape but the escalation fails closed on): no corpus carrier forces it
  (GT7081's grazes were all filter-precision, not genuine ties).
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

R1. THE SELECTIVE WELD MAY HAVE BEEN DISSOLVED WRONGLY. The resolver replaces a whole
    component and composes back, so the patch-to-exterior stitch is gone - but the
    resolver's own constructed intersection coordinates, rounded to double, can land
    sub-eps-distinct.  If the ordinary assembly weld (the global eps grid) merges two
    exact-distinct resolver verts it re-manufactures the twin the arrangement just
    resolved; if it merges a resolver-interior vert onto a boundary vert it corrupts
    the resolver's topology. The re-gate does NOT
    backstop this failure: when the weld merges two exact-distinct resolver verts, the
    triangles straddling the merge now share that vertex position, so
    IsSelfIntersecting's shares-vertex skip drops the pair - a self-FOLD manufactured by
    the merge is invisible to it. Only the validity half (IsManifold/Is2Manifold)
    catches a weld outcome, and only the non-manifold one (a pinch/tear), never a fold.
    Nor is a constant-radius bounded weld the escape: the wall-A weld-bump probe killed
    it at every multiplier (too small re-manufactures the twin, too big collapses
    slivers to holes). R1's honest rebuttal is the resolver's own STRUCTURAL defense, not the gate:
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
    flag near-misses - and that bias does not agree with the resolver's arrangement notion. The
    silent-miss direction is GATE-CLEAN but RESOLVER-DIRTY: a genuine crossing whose two verts
    sit in a near-degenerate band just outside the eps weld (distinct enough not to
    merge, close enough to trip the 2*eps skip) passes the gate, EARLY-EXITS as clean,
    and carries an unregularized self-overlap through silently - worse than fail-closed,
    though narrow (it needs a near-degenerate config in that thin band). The reverse
    (gate-dirty, the resolver finds nothing) is near-empty: the resolver's shared-vertex skip set is a subset
    of the gate's distance skip set, so nothing the gate robustly flags is skipped by the resolver;
    the only case is a legitimate touching contact the gate flags and the resolver re-emits clean -
    wasteful, not wrong. R2(i) (gate-clean/resolver-dirty) and R1 are the SAME blind spot: the
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

R5. WHAT IS LANDED VS OPEN. The gate, dispatch, coplanar fold, seam recording,
    coupled winding, the {w_S>=1} boundary emission, and the assembly are built and
    run the corpus. The named-open axes (a >2-sheet radial junction, the thin-cell
    emission, negative-winding subtraction, the component-local seed policy) each
    FAIL CLOSED with a distinct reason - never a silent wrong resolve - so a reader
    should read the open list as "fails closed here", not "unimplemented and unsafe".

KNOWN OUTPUT ARTIFACT (ent-verify, disclosure): on rotated/irrational junctions
the emission weld can materialize a single inert zero-area triangle (paired
edges, topology sound, oracle-exact volume; IsSelfIntersecting-blind).  The
clean fix is the once-only fold-vertex identity hardening (queued, simplicity
pass); the R1 weld blind spot remains the underlying open.


## Gap ledger (the authoritative census)

The complete enumeration of what still fails closed, folded to root mechanism.
Magnitudes only here; the exact per-carrier fatal strings, counters, and the
instrumented measurement live in .claude/lane-reports/reg3d-census-*.md (the
notebook owns the numbers).  Every entry below is either queued behind a closure
crucible or carries a proof-shaped adjudication; "research-grade" is never a bare
terminal state.

Census (each fatal-string SITE in overlap3.cpp -> its axis; F-ids are the
notebook's).  The sites, by kind, are: two EMISSION strings (BuildImpl: "unresolvable
sheet contact" from SplitTouchingSheets, "emitted triangulation not 2-manifold");
one COPLANAR-FOLD decline; the SEAM-SUBFACE "not exactly resolvable"; the
non-2-endpoint-seam "degenerate incidence"; the clean-face "filter-uncertain"
graze; four NEAR-COPLANAR planarize refusals (degenerate face / degenerate normal
/ global-planarity guard / inconsistent snap); the non-coplanar exact-zero SoS
tie; the "epsilon not computable" input guard; the "input not 2-manifold"
defensive gate; and the "resolver output failed the re-gate" R1/R2 backstop.

Live fail-closed carriers (measured at HEAD, ZZZ_CENSUS instrumented then
reverted):

  carrier          | cluster            | fatal (measured)              | reach
  -----------------|--------------------|-------------------------------|--------
  PokedCube        | (RESOLVES wjump)   | - (vol 0.25, GWN-oracle)      | constructed
  GT7863 pair      | (RESOLVES 7863c1)  | - (both dirty comps, GWN)     | corpus
  openscad soup    | 1 EMISSION (f4-junc| unresolvable sheet contact    | corpus
  GT7081 pair      | (RESOLVES wjump)   | - (both shells, vol-preserved)| corpus
  NearCoplanarChain| 4 PLANARIZE-GUARD  | global-planarity guard        | constructed
  BridgedCaps      | (RESOLVES)         | -                             | resolves

  reg3d-wjump LANDED (the winding-jump / degenerate-incidence EMISSION
  completion): the SHARES-VERTEX GENUINE-CROSSING RECOVERY (RecordSeams) closes
  the Cluster-1 emission wall for its self-intersecting carriers.  The shared-
  vertex broadphase skip was UNSOUND for self-intersecting soups (two faces
  sharing a corner that fold back and cross OFF the corner had their seam dropped,
  leaving an incomplete arrangement whose emission opened - the unbalanced-fan
  "8-edge hole").  Recovering those crossings (only genuine off-vertex transversal
  pierces, SoS-decided; the seam is [V, offVertexP], V the shared corner) COMPLETES
  the arrangement so the existing per-face witness rule emits the everted
  double-sheet correctly (decision-completion 2 is SUBSUMED: no new emission rule -
  s7b was right that the +1-only rule is already general on a COMPLETE
  arrangement).  RESOLVES (oracle-graded, zero-oracle-wrong): PokedCube (vol=0.25
  == GWN MC oracle, tol-invariant, non-self-int); GT7081 both dirty shells
  (volume-preserved, tol-invariant); GT7863's SELF-INTERSECTING component (the
  8-edge hole, vol 855, tol-invariant).  All resolving fixtures preserved
  (siA/siB stay in the MC volume band; the retention rule is byte-untouched).
  Mutation-verified (disable the recovery -> the resolves revert to "unresolvable
  sheet contact").  NO new predicate FORM (EdgePiercesTriSoS is the existing
  blessed SoS caller); no kernel escalation - REFUTING the prior lanes'
  "research-grade / plane-based-rep / kernel-tripwire" conclusion for these
  carriers.

  RECLASSIFICATION (reg3d-c2b -> reg3d-c2bx, measured): GT7081's original F4
  string read "seam sub-face arrangement not exactly resolvable"; reg3d-c2b's
  anatomy REFUTED the Cluster-2 reading (the 2D seam sub-face arrangement resolves
  EXACTLY - zero triple points) and re-scoped it to the WINDING-PROBE
  filter-precision residue (both dirty shells fail at the winding CLASSIFY probe
  grazing a near-coplanar shallow-dihedral face the FILTER cannot decide but the
  exact kernel decides NONZERO).  reg3d-c2bx then CLOSED that residue with the
  winding-probe escalation (WindingAt filter-0 -> Orient3DExactSign): both shells
  now decide the classify exactly (measured zero genuine ties) and fail at the
  DEEPER pre-existing wall - "unresolvable sheet contact" (SplitTouchingSheets ->
  NonManifoldEmission), the CLUSTER-1 emission representability wall IDENTICAL to
  GT7863's near-coplanar sliver (GT7081's 0.002deg near-tangent geometry reduces
  the {w>=1} boundary to touching sheets with no representable double-manifold).
  So GT7081 RECLASSIFIES AGAIN to Cluster 1a and folds into the C-1a crucible with
  GT7863.  Its winding-probe residue is closed; its terminal wall is emission.

  reg3d-f4-b1 LANDED (the ONCE-ONLY 3-face triple-point construction): oscad-f4
  measured that openscad's F4 "seam sub-face not exactly resolvable" wall is
  ENTIRELY the 3-face TRIPLE-POINT case (the winding probe fires ZERO times -
  REFUTING the "winding-probe / GT7081-class" reading logged for openscad in
  Cluster 2a/3 above).  Its large dirty component is triple-point dense, and the
  per-face resolver builds each triple point three times (once per incident
  face's frame), so the pos2in input-preimage map refuses the constructed
  crossings.  EnumerateTriplePoints now enumerates each triple ONCE, constructs
  ONE canonical double point keyed by the sorted plane triple (Intersect3Planes,
  a construction sibling of SegPlanePoint, NOT a decision), and threads it into
  all three incident faces by PRE-SPLITTING the seams at the exact on-seam
  crossing keyed to the shared point - so the arrangement completes at the
  0-cells and the triple-point refusal is DISSOLVED.  Tripwire-free: no exact
  arithmetic on constructed points, the crossing test and the winding classify
  stay the existing level-0 predicates, and the retention rule is BYTE-UNCHANGED.
  Bitwise no-op off openscad (zero triple points on every other corpus carrier).
  openscad then fails one wall DEEPER at the EMISSION wall ("unresolvable sheet
  contact", SplitTouchingSheets -> NonManifoldEmission), joining GT7081/GT7863 in
  Cluster 1.  MEASURED (f4-b1 census): the once-only welding is LOAD-BEARING (the
  per-face-reconstruct MUTATION reopens the triple-incident open edges) but
  INSUFFICIENT - it closes only the triple-incident holes.  The DOMINANT
  open-boundary residue is the DROPPED-BOUNDARY class AWAY from the triple points
  (the everted / high-cover strata; most open edges are single dangling
  halfedges with neither endpoint at a triple), plus a few EXACT-coincident
  radial-tangent ties (gap exactly zero = genuine tangent sheets, not double-
  rounding) and material overlaps.  This CONFIRMS design-a's radial-arrangement
  residue and design-c's stated FALSIFIER (open edges PERSIST on a once-only-
  consistent arrangement with the emission rule unchanged), and REFUTES design-b's
  "once-only closes the near-triple majority" prediction (its near/far split was
  confounded, as its own notebook warned) AND design-a's winding-jump re-emission
  rule (measured to move the count the WRONG way).  Closing the residue is the
  exact-radial substrate (the tripwire), RECORDED not built.  (An earlier aside
  here claimed the pre-B1 wall was MANIFOLD_PAR-dependent; the verification lane
  could not reproduce that - pre-B1 is F4 under both settings in from-scratch
  builds.  With B1 the terminal is the emission wall under both settings.)

  reg3d-f4-junction LANDED (the emission residue NARROWED in-form, byte-for-byte
  no-op off openscad).  (1) BOTH-SIDES RETENTION: EmitSeamedFace now probes w_S
  on BOTH sides of each sub-cell and retains iff exactly one is inside {w_S>=1}
  (the coplanar fold's general rule, the m==1 case of which is byte-identical to
  the former w_above==0 rule), closing the radial ties + material overlaps - a
  one-sided-rule artifact where the seam sub-cell's true jump != 1.  (2) GLOBAL
  JUNCTION REGISTRY: the per-face arrangement wove PROPER seam crossings (B1) but
  not the NON-proper-crossing junctions - a seam endpoint (or a triple that
  TERMINATES a seam on the third face) landing strictly interior to a neighbour /
  partner / third face's emitted edge without a shared split, which opens the
  fan.  BuildJunctionRegistry dedups every once-only arrangement vertex (seam
  endpoints + triples) to a canonical 3D point, and every emit path (seamed,
  clean, fold) pre-splits its edges at each registry vertex strictly interior to
  it - the split DECISION a pure 3D on-segment test on the SHARED endpoints (so
  both incident faces decide identically; a per-frame 2D test disagreed near
  endpoints and manufactured fresh T-junctions), the 2D foot only the on-line
  addAt position (no RemoveOverlaps2D fold-back).  No new predicate FORM, no
  exact-on-constructed-point (level-0 3D on-segment on existing once-only
  constructions), so the design lanes' "exact-radial substrate needed for the
  T-junction residue" reading is REFUTED for the dominant class: the registry
  closes ~two thirds of the open edges in-form, mutation-verified load-bearing
  (registry off reopens them), FNV-identical on all resolving carriers.  The
  SURVIVING residue is near-degenerate and stays fail-closed: near-tangent
  >2-SHEET RADIAL junctions (the unbuilt radial branch), genuinely-distinct
  near-coincident vertices from near-PARALLEL planes, a coplanar/transversal
  ENTANGLEMENT line (the fold's cap-cap crossing and a seamed wall's endpoint
  differ by more than eps), and collinear clean-clean overlaps - the
  plane-based-representation escalation the campaign gates behind the kernel
  tripwire.  openscad stays at the emission wall, now much narrower.

MEASURED CORRECTIONS to earlier anatomy (were logged at a pre-entanglement,
pre-stage-6 HEAD): (1) PokedCube AND GT7863 fail at the SplitTouchingSheets
"unresolvable sheet contact" EMISSION string, which fires BEFORE the 2-manifold
gate - not at "not 2-manifold" as reg3d-arr's open-halfedge dump implied.  (2)
BridgedCaps RESOLVES (stage-6 SoS closed it); its in-test comment predates that
landing and is stale.

Cluster map (root mechanism, not which string fires):

- CLUSTER 1 - EMISSION REPRESENTABILITY (stage-7 thin-cell / touching sheet).
  The resolved {w_S>=1} boundary reduces to sub-eps / touching sheets with no
  representable double-manifold; BuildImpl declines.  UNIFIED ANATOMY (reg3d-c1a):
  every carrier fails SplitTouchingSheets at an UNBALANCED fan = an OPEN BOUNDARY
  (dropped faces), NOT a touching/doubled sheet - measured zero surviving weld
  twins and zero non-manifold edges on all of them; the wall needs MORE emitted
  faces, so C-5's twin-identity unification is a measured no-op here (see closure
  plan 5).  Two sub-shapes:
  1a EXACT-INCIDENT / everted-emission (GT7863 AND GT7081 post-reg3d-c2bx;
  root label CORRECTED by the bake-off, reg3d-1a/1b: NOT a near-coplanar input
  sliver - the "887-2259 ULP sliver" is CROSS-mesh pass-through; the dirty-component
  wall is exact-incident edge-on-edge seam endpoints <1 ULP + a double-sheet /
  everted-region winding drop, the reg3d-wjump completion; options a/b both refuted
  with empty patches).  CLOSED for its named carriers: reg3d-wjump resolved
  GT7863's 8-edge hole + GT7081 + PokedCube; reg3d-7863c1 resolved GT7863's
  coplanar comp#1 - both were ARRANGEMENT-INCOMPLETENESS from an unsound
  shares-vertex skip, not an emission-representability wall.  1b negative-winding
  double sheet
  (PokedCube: an everted corner makes a genuine w_S=-1 region, so the w=-1|w=1
  junction is a double sheet the mult-1 per-face emission opens - s7b refuted the
  simple orientation flip).  Production-reachable.
- CLUSTER 2 - SEAM SUB-FACE EXACT RESOLVABILITY.  RecordSeams / EmitSeamedFace
  cannot build a face's 2D arrangement exactly at a degenerate incidence.  2a
  openscad's residue was RECLASSIFIED OUT of this cluster (reg3d-oscad REOPEN, LANDED): the
  earlier per-pair reading (reg3d-c2a) treated its nPts==1 truncations as a cap-seaming
  entanglement / degenerate-endpoint wall, but exact rational reconstruction of every truncated
  seam's true second endpoint showed the dominant residue is MEASURE-ZERO contacts
  (coincident-position DUPLICATE vertices, vertex-on-edge / vertex-on-face T-junctions,
  collinear edge-on-edge overlaps) - a POINT or boundary-segment intersection with NO
  transversal crossing, a PHANTOM seam the filter/SoS manufactured (proven exactly: none has a
  clean off-plane edge piercing the other's strict interior).  Three RecordSeams
  decision-completions - the shares-vertex skip/recovery keyed on POSITION coincidence, a
  phantom-seam guard, a sub-eps collapse - close the F11 seam-truncation wall; openscad then
  fails at the WINDING CLASSIFY probe (F4), moving to the winding-probe axis (Cluster 3) below,
  the SAME near-tangent phenomenon as GT7081.  2b
  GT7081's shells were RECLASSIFIED OUT of this cluster (reg3d-c2b): the anatomy
  shows their 2D seam sub-face arrangement resolves EXACTLY (preimage-strict
  reconstruct, zero triple points); the fatal is the WINDING CLASSIFY probe
  grazing a near-coplanar face on the filter, so they move to the WINDING-PROBE
  axis (Cluster 3 / O4 seed policy) below.  Production-reachable.
- CLUSTER 3 - SoS / EXACT-TIE RESIDUE + the unbuilt >2-sheet radial branch.  The
  >2-sheet radial branch is not corpus-forced (reg3d-radial: zero book-of-pages
  lines; every corpus nPts!=2 is a 2-sheet truncation, not a triple point); its
  reduction is proven sound but fires zero times.  The WINDING-PROBE
  FILTER-PRECISION residue (reg3d-c2b, GT7081) that once lived here is CLOSED
  (reg3d-c2bx): the winding CLASSIFY probe (WindingAt, both the seamed per-cell and
  clean-face paths) previously decided ray crossings on the FILTER only
  (Orient3DFilterSign) and failed closed on a filter-0; it now ESCALATES each
  filter-0 term to Orient3DExactSign (the blessed winding-probe caller, filter-first).
  GT7081's near-tangent 0.002deg geometry put the constructed probe (cell centroid
  + eps*n) within the filter's uncertainty band of a face whose exact sign is
  NONZERO (decidably off-plane); the escalation decides it exactly (measured zero
  genuine ties).  The O4 component-local seed reprobe was MEASURED insufficient
  (reg3d-c2b: it cleared multi-triangle seamed cells but not a single-triangle
  clean face near-coplanar over its whole extent) AND is now redundant (the
  escalation clears the graze exactly, so it was adjudicated OUT, reg3d-c2bx).
  Closing the winding-probe residue did NOT resolve GT7081 end-to-end: it revealed
  the pre-existing Cluster-1a emission wall underneath (both shells now fail at
  "unresolvable sheet contact"), so GT7081 moved to Cluster 1a.  The
  coupled-integer-flood re-architecture of the WIND phase (which would remove the
  ray-cast winding altogether) remains a named research axis, no longer needed for
  this residue.  openscad's large component moved PAST this axis: reg3d-f4-b1
  (once-only triple points) then reg3d-f4-junction (both-sides retention + the
  global junction registry) carried it through the seam sub-face winding classify
  and the T-junction arrangement completion to the EMISSION wall (Cluster 1a),
  where its near-degenerate residue - near-tangent >2-sheet radial junctions plus
  the coplanar/transversal entanglement - joins GT7081/GT7863 behind the
  plane-based-representation tripwire.  Fail-closed, no wrong resolve.
- CLUSTER 4 - NEAR-COPLANAR PLANARIZE GUARD / SNAP refusals.  A curved near-band
  is genuinely not one plane; the guard fails closed rather than fold to a wrong
  plane.  These are CORRECT, decision-complete refusals (the exact procedure DOES
  cover them), not unbuilt axes.  Not corpus-forced (the corpus's real
  near-coplanar geometry is cross-component pass-through).
- CLUSTER 5 - WELD / GATE BLIND SPOTS (R1 / R2 / the zero-area artifact).  The
  uniform emission weld can merge two genuinely-distinct arrangement points within
  eps (a fold the shares-vertex-skip gate is blind to) or leave a 1-ULP-distinct
  vertex unmerged (the disclosed zero-area triangle).  No corpus carrier
  demonstrated; live on the synthetic rotated pin.  reg3d-c1a REFUTED the
  hypothesis that the once-only fold-vertex identity (C-5) could double as the
  emission-wall fix (measured no-op; the emission fans carry no surviving twins)
  AND showed C-5 "proper" is not a bounded DC but the plane-based-rep escalation
  (the reentrant-corner point has three distinct cross-phase constructions;
  unifying them by construction needs one canonical symbolic vertex).
- CLUSTER 6 - DEGENERATE-INPUT / DEFENSIVE GUARDS ("epsilon not computable",
  "input not 2-manifold").  Terminal-correct contract edges, not gaps.

Closure plan (ordered by production-reachability x carriers-unlocked x closure
shape; DC = decision-completion, PROOF = proof-shaped terminal adjudication,
RESEARCH = memo with a required proof sketch, TRIPWIRE = kernel-vendor decision):

1. CRUCIBLE C-2a (RESEARCH - reclassified from DC; reg3d-c2a).  The census scoped
   this as a bounded cap-INTERIOR pierce injection; the per-pair measurement + a
   constructed minimal fixture REFUTE that.  A cap-interior seam endpoint arises
   ONLY from a wall piercing a cap face, which SEAMS the cap = the coplanar/
   transversal entanglement the fold declines by design; the F11 seam truncation is
   only a symptom, with the F3 fold-decline underneath (proven: the same fixture
   with the plug RESTING on the cap - no pierce - resolves; piercing it fails
   closed; openscad's failing component was once read as likewise cap-seam-entangled -
   SUPERSEDED, see TERMINAL below).  Closing THIS constructed fixture needs the fold to
   arrange around the transversal seam and classify each split sub-cell by the real 3D
   coupled winding = the Cluster-1 coordinated-emission wall, not a bounded completion.
   openscad's F11 wall, by contrast, WAS closed by a RecordSeams-only fix (reg3d-oscad
   REOPEN, TERMINAL below): its truncations were not this fixture's cap-seam-entanglement
   but measure-zero phantom contacts.  Fail-closed, pinned
   (Regularize_CapSeamingEntanglement_CapInteriorPierce_FailClosed).  The
   cap-seaming-entanglement PROOF-PAIR (this constructed fixture: plug PIERCING a
   cap fails closed, plug RESTING on the cap resolves) stands as the class exemplar.
   TERMINAL SUPERSEDED (reg3d-oscad REOPEN, LANDED): the earlier
   "plane-based-representation PROOF wall / ZERO shares-vertex / no fold-owned
   completion" reading was a MEASUREMENT GAP - it searched only for a dropped interior
   pierce and measured shares-vertex by INDEX, but the unwelded-DUPLICATE family IS the
   shares-vertex family keyed on POSITION.  Exact rational reconstruction of every
   truncated seam's true second endpoint showed the dominant residue is MEASURE-ZERO
   contacts (coincident-position duplicate vertices, vertex-on-edge / vertex-on-face
   T-junctions, collinear edge-on-edge overlaps) - a POINT or boundary-segment
   intersection, a phantom seam with no transversal crossing (proven exactly: none has a
   clean off-plane edge piercing the other's strict interior).  Three RecordSeams
   decision-completions close the F11 wall (no new predicate FORM): position-keyed
   shares-vertex skip/recovery, phantom-seam guard, sub-eps collapse.  openscad then fails
   one wall DEEPER at F4, the seam sub-face WINDING CLASSIFY probe = the near-tangent
   winding-probe class, the SAME exact-kernel-load phenomenon as GT7081 (Cluster 3).
   Still fail-closed.  (The capvert recovery WIDENING remains separately refuted - a
   benign wall-corner-on-cap-plane touch is not a seam endpoint, so recording touches
   over-recovers; the ent proper-cross gate stays load-bearing for soundness.)
2. CRUCIBLE C-2b (LANDED escalation; residue -> C-1a; reg3d-c2b -> reg3d-c2bx).
   reg3d-c2b's anatomy REFUTED the Cluster-2 reading (GT7081's 2D seam sub-face
   arrangement resolves EXACTLY; both shells fail at the WINDING CLASSIFY probe
   grazing a near-coplanar face on the filter, exact-decidable NONZERO).  The owner
   RELAXED the tripwire to ONE PREDICATE, ONE IMPLEMENTATION (same
   Orient3DExactSign, additional CALLERS blessed, zero new arithmetic), so
   reg3d-c2bx LANDED closure (a): WindingAt's filter-0 terms ESCALATE to
   Orient3DExactSign (the winding-probe caller, filter-first).  MEASURED: both dirty shells
   clear the winding classify exactly (zero genuine ties, so no wrong resolve), but
   the escalation reveals a DEEPER pre-existing wall - both fail at emission
   ("unresolvable sheet contact"), the Cluster-1a near-coplanar-sliver wall
   IDENTICAL to GT7863.  So GT7081 does NOT resolve end-to-end; it RECLASSIFIES to
   Cluster 1a and its terminal closure is C-1a (below).  The pin flips
   FatalReason DirtyComponentUnresolved -> NonManifoldEmission (narrower, mutation-
   verified: disabling the escalation reverts it).  The O4 seed-policy reprobe was
   ADJUDICATED OUT (the escalation makes it redundant; measured, reg3d-c2bx).  The
   coupled-integer-flood WIND re-architecture (b) is no longer needed for this
   residue.  SURFACED earlier (owner triage, out of scope): the capvert
   recovery-WIDENING recovers openscad's fold-owned vertex-on-plane touches but F3
   is UNAFFECTED, so it does NOT unlock openscad alone - an emission-side follow-up.
   (openscad is UNCHANGED by the escalation: it fails UPSTREAM at RecordSeams
   "degenerate incidence"; the winding probe fires zero times, reg3d-c2bx.)
3. CRUCIBLE C-1a (PROOF; reg3d-c1a UNIFIED EMISSION ANATOMY).  The consolidated
   emission wall ("unresolvable sheet contact", SplitTouchingSheets ->
   NonManifoldEmission) was instrumented on ALL its carriers (PokedCube, GT7863,
   GT7081 both shells; cap-seaming fixture).  MEASURED (ZZZ_C1A, reverted
   byte-clean): every emission carrier fails at the SAME branch - an UNBALANCED
   fan (fwd/bwd count mismatch) = an OPEN BOUNDARY / hole, with ZERO surviving
   near-duplicate vertices in any band (eps..16*eps) near the failing edge or
   anywhere in the emitted vert set, and zero non-manifold edges (no doubled
   sheets).  PokedCube's 8-triangle emission reconstructs into three closed
   boundary loops (a 7-gon through the everted spike + two ISOLATED islands whose
   connecting neighbors were dropped).  GT7863 = an 8-edge hole in a near-coplanar
   x-sliver (x span ~0.76 at x~-31165); GT7081 = the SAME wall at a 0.002deg
   shallower dihedral (112 + 6 open edges across two shells, winding classify now
   exact per reg3d-c2bx, so only the near-tangent sliver emission remains).
   PROOF (genuinely unpairable from the information present): each fan is an open
   boundary because faces on d{w>=1} were DROPPED by the per-face emission at a
   degenerate corner (PokedCube's w=-1|w=1 double-sheet everted spike) or a
   near-coplanar / near-tangent site (GT7863/GT7081 slivers) - the dropped faces
   are unrecoverable from the rounded-vec3 per-face arrangement + winding.  The
   coordinated fixes converge: PokedCube needs the winding-jump-of-2 double-sheet
   re-emission + an EXACT arrangement at the degenerate corner (triple points);
   GT7863/GT7081 need a stage-5 near-coplanar fold WIDENING above eps under a
   per-sub-face thinness bound (which risks the global-planarity-guard /
   curved-chain oracle-wrongness).  Both fix shapes are the PLANE-BASED
   REPRESENTATION escalation (exact symbolic vertices -> exact triple points +
   at-rounding thin-cell/sliver decisions), which crosses the exact-kernel tripwire
   (vendor Shewchuk).  RECORDED (build nothing on it).  The cap-seaming fixture is
   NOT on this wall - its failing component dies UPSTREAM at RecordSeams F11/F3
   (its non-piercing control emits clean); the C-2a research adjudication stands.
   BAKE-OFF CORRECTION (reg3d-1a / reg3d-1b-optB, both REFUTED-WITH-EVIDENCE, empty
   patches): the "1a near-coplanar sliver" ROOT LABEL is corrected by
   dirty-component measurement.  Option (a) [stage-5 near-coplanar WIDENING] has NO
   admissible target - zero near-coplanar 2D-AREA overlap below eps (or 100*eps) in
   EITHER carrier's dirty components; a 1000x widen mutation still leaves GT7863
   fail-closed.  Option (b) [thinness-aware WELD] has NO target - ZERO
   face-collapsing weld merges corpus-wide; the defect is ~500-1000x the weld radius.
   The documented "887-2259 ULP near-coplanar sliver" (planerep probe) is CROSS-MESH
   (LxR): the resolver decomposes L/R into separate components (dirty comps carry
   ~0-1 self-crossings, not the ~266 LxR interpenetration), so that wall is
   CROSS-COMPONENT PASS-THROUGH under the standing non-fusion contract - it never
   reaches per-component emission (folding it = breaking non-fusion, out of scope).
   The ACTUAL dirty-component wall is: (i) EXACT-INCIDENT edge-on-edge seam endpoints
   (<1 ULP = genuine form-independent ties, SoS jurisdiction) + winding classify;
   (ii) a DOUBLE-SHEET / everted-region emission (GT7863: 2 faces dropped at
   w_above=-1, GWN-confirmed everted; GT7081: 18 dropped faces) whose closure needs
   the winding-jump / degenerate-incidence EMISSION completion (this crucible,
   reg3d-wjump), NOT input-snap (a) and NOT weld (b).
   CLOSED (reg3d-wjump LANDED).  The root was NOT a near-coplanar sliver at all: it
   was ARRANGEMENT INCOMPLETENESS from the SHARES-VERTEX broadphase skip.  The
   "8-edge hole" fans opened because seams between faces sharing a corner that
   ALSO cross transversally off the corner (the everted-spike / near-triple-point
   geometry) were DROPPED by `sharesVert => continue` - unsound for a
   self-intersecting soup.  Recovering only the genuine OFF-VERTEX transversal
   crossings (SoS-decided pierce; seam [V, offVertexP], V the shared corner
   supplied as the second endpoint when the seam assembly yields nPts=1) COMPLETES
   the arrangement, and the EXISTING per-face witness rule then emits the everted
   double sheet correctly (the everted w=-1 region drops, the {w>=1} boundary is
   kept - decision-completion 2 is SUBSUMED, no new emission rule).  RESOLVES
   oracle-graded: PokedCube (vol=0.25 == GWN oracle), GT7081 (both shells,
   volume-preserved), GT7863's self-intersecting component (vol 855).  No new
   predicate FORM, no kernel escalation - so the earlier plane-based-rep /
   Shewchuk-tripwire conclusion is REFUTED for these carriers.
   GT7863 comp#1 CLOSED (reg3d-7863c1).  The prior "honest wall" reading was
   WRONG about the mechanism (the fold-escalation was the wrong lever - it was
   BANKED and stays UNBUILT, never landed).  comp#1 is NOT a clean doubled sheet:
   it is the twin composition's FLAT FACE triangulated with overlapping tiles
   (including a near-collinear sliver).  DetectCoplanarClusters was SKIPPING
   shares-vertex pairs, so the flat face's tiles that touch the overlap AT A
   CORNER were left OUT of the fold cluster.  The partial cluster's in-plane
   cover then disagreed with the 3D winding (the self-check fired, m=2 vs jump=1)
   and forcing emission T-junctioned against the un-clustered coplanar neighbours
   (an OPEN-BOUNDARY emission fail - the Cluster-1a anatomy).  FIX (decision-
   completion, arrangement-completeness - the exact shape reg3d-wjump used on the
   RecordSeams shares-vertex skip): a two-pass cluster detect - SEED from
   distinct-patch overlaps (share no vertex), then EXTEND a seeded cluster
   through shared-corner overlaps, but do NOT seed a cluster from a
   shared-corner-only overlap (that is a local FOLD-BACK, e.g. an everted-spike
   face pair, which the per-face winding rule owns - PokedCube stays bitwise).
   Completing the cluster makes the in-plane cover match the winding (m=1==jump),
   the self-check PASSES with NO escalation, and the existing mult+winding rule
   emits one sheet.  comp#1 RESOLVES oracle-true (vol ~2830 == signed volume ==
   GWN, tol-invariant, non-self-int); both GT7863 dirty components now resolve so
   the whole compose regularizes.  Mutation-verified (revert the skip narrowing
   -> comp#1 reverts to the coplanar-fold decline).
4. CRUCIBLE C-1b (CLOSED by reg3d-wjump; was folded into C-1a).  PokedCube's
   everted corner IS a w=-1|w=1 double-sheet AND a degenerate near-triple-point -
   but the near-triple-point was the shares-vertex-skipped crossings, not a
   >2-sheet exact triple point: recovering them completes the arrangement and the
   witness rule emits it (vol=0.25, GWN-oracle-confirmed).  s7b's clean-patch
   per-face fix (already landed) + this recovery together close it.
5. CRUCIBLE C-5 (REFUTED as the emission-wall key; residual = the plane-based-rep
   escalation).  reg3d-c1a MEASURED C-5's predicted effect on every emission fan
   and found it NULL: C-5 unifies within-eps twin identities, but it only changes
   the SplitTouchingSheets input where a twin SURVIVES the BuildImpl eps-weld
   (>eps apart) - and there are ZERO such twins on any carrier (the fans need MORE
   faces, not merged twins).  So the census hypothesis that construction-keyed
   identity could make the fans pairable is refuted.  A SOURCE read (RecordSeams
   recovery + dedup) further shows C-5 "proper"
   (bit-identical BY CONSTRUCTION rather than the current positional dedup) is NOT
   the bounded DC the census assumed: the symmetric double-pierce is two different
   SegPlanePoint constructions of one reentrant corner, and the fold overlay builds
   it a third way - unifying them requires ONE canonical construction shared across
   phases (RecordSeams runs BEFORE FoldCoplanarClusters), i.e. the plane-based
   representation.  The A4 zero-area artifact (synthetic rotated pin, inert, no
   corpus carrier) is the only thing C-5 would cure, and it too routes to the same
   escalation.  The honest R1 residual (two genuinely-distinct points rounding
   within eps) stays open.
6. ADJUDICATION - Cluster 3 (TRIPWIRE).  The >2-sheet radial rule is proven sound
   and fires zero times on the corpus; building it adds a second exact call site,
   so it stays unbuilt behind the vendor-Shewchuk decision (owner's).  Fail-closed
   posture stands.

Terminal-correct, no crucible (PROOF): Cluster 4 (a curved chain is not a single
planar fold - snapping it is the oracle-wrong outcome the guard prevents), Cluster
6 (degenerate-input guards), and Cluster 3's fail-closed posture until vendored.
