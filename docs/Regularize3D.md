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
  ONE PREDICATE FORM, ONE IMPLEMENTATION (reg3d-c2bx; generalized in homog-design):
  the integer path is legitimate as ONE exact predicate FORM - the sign of a
  homogeneous orientation determinant corrected by sign(prod W_i), summed on the ONE
  adaptive-width integer accumulator (sos::SumSignN).  It has exactly TWO
  instantiations: (1) INPUT-POINT orient3d (degree 3, W==1: an input double point is
  the intersection of its three trivial axis planes, denominator 1, so the W column
  is the literal ones column and the weight-product sign is +1 - the loop is
  BIT-IDENTICAL to the historical adaptive-integer orient3d, every existing caller
  unchanged; differential harness 0 diff over >1e7 inputs incl forced ties), and
  (2) CONSTRUCTED-POINT orient2d (degree 9: three in-face seam crossing points, each
  the Cramer intersection of a plane triple {F,g,h}; the SAME form at higher degree,
  filter-first via a degree-9 construction-aware static bound, on the widened
  accumulator - the naive final-determinant permanent is UNSOUND in the near-parallel
  wedge regime and must be replaced by the all-abs companion Pdet).  ADDITIONAL
  CALLERS of either instantiation are within the blessed pattern (zero new arithmetic)
  as long as each stays FILTER-FIRST (exact fires only behind a filter 0).  FIRE CENSUS (orient-land
  necessity attack): the exact cascade earns its keep in two roles.  (i) STRUCTURAL-TIE
  DETECTION - the bulk of the tie-cascade fallback: genuine exact-zeros, overwhelmingly
  STRUCTURAL (two boxes sharing an axis coordinate, or a repeated point), where the e^0
  sum is provably zero; these route to the single-global SoS, which breaks them by
  perturbation.  (ii) NEAR-TANGENT SIGN DECISIONS - a filter-uncertain but
  exact-decidably-NONZERO e^0, concentrated almost entirely on the largest corpus
  component (GT7081) and its ~0.002deg twin faces; this is the winding-probe graze the
  exact sign resolves.  CALLER INVENTORY of the exact e^0 predicate (Orient3DExactSign):
  the EdgePiercesTriSoS edge-in-plane guard, the WindCrossTri winding-crossing escalation
  (both the WindingAt winding-probe AND the RecordSeams phantom-seam guard ride it -
  cleanPierce is no longer a distinct exact-call site), the junction registry's
  input-vertex-on-edge arm (exact collinearity/between-ness on input doubles), and
  the triple-point seam-crossing test (ExactSegProperCross - exact in-plane
  orient2d, drop the dominant normal axis, refuting the near-tangent phantom
  crossings the former rounded-projection double crossing test it replaced
  over-detected in EnumerateTriplePoints; byte-identical on the resolving corpus, removes only
  phantom triples on the near-tangent openscad residue) - the production callers (plus the test
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
  What stays tripwired is a THIRD instantiation of a NEW DEGREE (or any new
  constructed-point form beyond the two above): an OWNER DECISION - it widens the
  accumulator's proven totality bound (kAccumLimbs, now 320 to keep the degree-9 form
  total by construction) and needs a new per-degree filter constant; never add one
  silently.  And if an exact primitive OUTSIDE this one homogeneous form is ever
  needed, VENDOR Shewchuk's public-domain predicates.c instead of growing this - do
  NOT rebuild expansion arithmetic piecemeal.
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

Census (each fatal-string SITE in overlap3.cpp -> its axis).  The sites, by
kind, are: two EMISSION strings (BuildImpl: "unresolvable sheet contact" from
SplitTouchingSheets, "emitted triangulation not 2-manifold"); one COPLANAR-FOLD
decline; the SEAM-SUBFACE "not exactly resolvable"; the non-2-endpoint-seam
"degenerate incidence"; the clean-face "filter-uncertain" graze; two
NEAR-COPLANAR planarize refusals (global-planarity guard / inconsistent snap);
the "epsilon not computable" input guard; the "input not 2-manifold" reachable
defensive guard (fires only on the direct-Impl path); and the "resolver output
failed the re-gate" R1/R2 backstop.  (There is no "SoS refused" fatal slot: the
single-global SoS is TOTAL by construction - it decides every non-coplanar
transversal exact-zero tie - so a residue surfaces instead as the non-2-endpoint
seam or downstream at emission.)

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

The census BY KIND (current state).  The lane reports own the history and the
numbers - .claude/lane-reports/reg3d-*, f4-*, attack22-*, census2-* - so this
ledger states, per gap: what forces it, the counterexample it refutes, and (for
the one terminal wall) the spec and price to cross it.  No history-narration
here; the notebooks carry the landings.

REACHABLY-TRIPPED ARMS (a live carrier forces each):

- E1 EMISSION - "unresolvable sheet contact" (SplitTouchingSheets open
  boundary).  openscad's large dirty component (corpus), the sole
  PRODUCTION-reachable completeness terminal.  COUNTEREXAMPLE: the resolver
  cannot emit a representable double-precision {w_S>=1} boundary for it.  Its
  arrangement COMPLETES (once-only triple points + both-sides retention + the
  global junction registry carry it past the seam sub-face winding classify and
  the T-junction completion), but the emitted boundary reduces to an
  open-boundary fan whose missing sub-faces are bounded by NEAR-TANGENT seam
  crossings - clusters of genuinely-distinct triple points spread from ~1e-4 down
  to ~2*eps that RemoveOverlaps2D's eps-merge collapses to one cell spanning two
  winding regions.  A finer per-sub-triangle read is measured WORSE (it unpairs
  real neighbours).  Zero-oracle-wrong holds: no geometry is emitted, and every
  probed winding-only completion face was winding-validated and rejected, so the
  residue is a representability wall, not a solvable pairing.  Refuted closure
  families, each with evidence: stage-5 near-coplanar widening (no admissible
  sub-eps target), thinness-aware weld (no face-collapsing merge), winding-only
  boundary completion (over-fills non-uniform regions the winding rejects),
  branch-node radial order as the blocker (the order is level-0 FREE from
  input-plane normals + the emitted fwd bit), and the exact axis-drop projection
  (relocates the fan imbalance, does not resolve).
  SPEC (the only path across the wall): an EXACT per-face near-tangent 2D
  arrangement over CONSTRUCTED intersection points.  The single tripwire crossing
  is orient2d on a constructed near-tangent point (not on input coords).
  STATUS (homog-design): the exact predicate is now LANDED and VALIDATED - the
  degree-9 constructed-point orient2d (sos::HomogOrient2DExact / the
  construction-aware HomogOrient2DFilter, instantiation (2) of the ONE homogeneous
  form, on the widened SumSignN<8> accumulator; validated ZERO disagreements vs an
  exact-rational oracle over 1.3M random + near-parallel-wedge plane triples).  It
  is wired at the seam-crossing enumeration (ExactTripleStrictlyInFace) to REFINE
  the existence decision - the exact constructed crossing must be strictly
  interior to the face - which refutes the rounded-straddle phantoms byte-cleanly
  on the whole resolving corpus.  MEASURED (this landing): the exact predicate
  alone does NOT close the openscad residue.  Replacing the rounded-endpoint
  segment straddle OUTRIGHT with the constructed-crossing test OVER-detects (it
  ignores the seam SEGMENT extent - the seam LINES cross inside the face at ~30x
  the rate the finite overlap SEGMENTS do), manufacturing spurious triples that
  break the resolving carriers.  So the residue's remaining net-new surface is NOT
  the crossing predicate (landed) but the seam ENDPOINTS carried as SYMBOLIC
  plane-triples {F, edge-face-a, edge-face-b} (their once-only exact extent), so
  the straddle itself becomes exact - the ~600-1300 LOC exact-arrangement plumbing
  f4-r6/nomerge priced separately.  The radial order and the winding half remain
  FREE.
  (GT7863, GT7081, and PokedCube reached this wall historically but now RESOLVE
  oracle-graded: their residue was ARRANGEMENT-INCOMPLETENESS from an unsound
  shares-vertex broadphase skip, not an emission-representability wall.  openscad
  is the sole surviving corpus carrier at E1.)

- A1 SEAM - "non-2-endpoint seam / degenerate incidence" (RecordSeams).
  CapSeamingEntanglement's cap-interior pierce (constructed).  A wall piercing a
  cap face SEAMS the cap - the coplanar/transversal entanglement the fold
  declines by design.  COUNTEREXAMPLE-PAIR: the same fixture with the plug
  RESTING on the cap (no pierce) RESOLVES; piercing it fails closed.  Closing the
  constructed pierce is the same coordinated-emission wall as E1, not a bounded
  completion.  Decision-correct.

- S3 PLANARIZE - "global-planarity guard, curved chain"
  (SnapNearCoplanarClusters).  NearCoplanarChain (constructed).  A curved
  near-band is genuinely not one plane; the guard fails closed rather than fold
  to a wrong plane.  A CORRECT, decision-complete refusal, not an unbuilt axis.
  Not corpus-forced (the corpus's real near-coplanar geometry is cross-component
  pass-through).

- D1 INPUT - "epsilon not computable" (RemoveOverlaps3D).  Zero-scale input (a
  cube collapsed to a point -> bbox scale 0 -> eps 0).  Fires before decompose.
  A degenerate-input contract edge.

OFF-CORPUS / DEFENSIVE BACKSTOPS (no corpus or constructed carrier forces them;
each fails CLOSED on its trigger - correctness by code-read, never a wrong
resolve):

- R1 coplanar-fold decline - fires only if RecordSeams passes but
  FoldCoplanarClusters declines (degenerate projection / filter-uncertain
  in-plane classify / coplanar-transversal entanglement).  The one entanglement
  carrier trips A1 first, masking R1.
- R2 seam sub-face not exactly resolvable - the historical openscad F4,
  dissolved by the once-only triple-point construction.  Remaining trigger = an
  UNBUILT >2-sheet radial junction; a co-axial 3-sheet carrier is the named
  synthetic.  The reduction is proven sound but fires zero times on the corpus.
- R3 clean-face winding probe filter-uncertain (SoS) - fires only on a genuine
  exact-zero seed-position graze over a clean face's whole extent (the
  component-local seed-policy open).  No corpus carrier: the winding probe now
  escalates filter-0 to Orient3DExactSign, so corpus grazes decide exactly.
- S4 vertex in two near-coplanar clusters (inconsistent snap) - reachable in
  principle (the shares-vertex skip lets two near clusters share a vertex), but
  needs a closed soup with two near-coplanar internal overlaps meeting at a
  shared vertex, not cheaply realized.
- E2 emitted triangulation not 2-manifold - BuildImpl's second emission backstop
  after SplitTouchingSheets passes.  Every corpus emission failure is E1 (open
  boundary), never E2.
- D3 resolver output failed the re-gate - a resolve that passes BuildImpl's
  manifold gate but is self-intersecting (a weld-manufactured fold).  Verified
  unreached on constructible general-position fixtures; a RELEASE fail-closed,
  deliberately not a DEBUG_ASSERT.

DEAD-BY-CONSTRUCTION (the two near-coplanar refusals DELETED in the simp4 fold):
the "degenerate face" and "degenerate normal" branches in
SnapNearCoplanarClusters were unreachable - a degenerate face has pairGap ==
+inf and never unites, so a size>=2 cluster's first face is always
non-degenerate and the sign-aligned member-normal sum is always nonzero.  See
the in-code invariant comments.

LATENT LIABILITIES (measured-dangerous, oracle-masked today; census2).  L1 - the
sharpest, the only one that had been influencing resolving-carrier output - is now
RESOLVED (below); L2-L4 remain:

- L1 PLANEFRAME PROJECTION - RESOLVED (exact axis-drop face-flattening).  The
  resolver formerly flattened each face to 2D through PlaneFrame's orthonormal
  e1/e2 basis, whose dot-product projection carried rounding; on openscad's large
  dirty component that rounding made crossing decisions and a few eps-merges that
  DIFFERED from an exact axis-drop (the convicted overlay liability).  It now
  flattens by EXACT axis-drop: the 2D projection SELECTS the two non-dominant
  coordinates of the dominant-normal-axis drop (pure coordinate selection, no
  arithmetic -> no rounding), with a parity swap so 2D-CCW still maps to +nHat
  and a plane-solve LIFT (solve nHat.P = planeD for the dropped dominant
  coordinate, |nHat[axis]| >= 1/sqrt(3) so never near-zero) for the only 2D-BORN
  points - a fold new-crossing vertex and the fold interior classify point; every
  other emitted vertex keeps its exact canonical 3D.  PlaneFrame is DELETED: the
  whole resolver's overlay predicates now run on exact input/construction doubles
  at level 0.
  RE-BLESS (equivalence, not regression): oracle-graded per carrier against the
  pre-swap output - identical outcome class, triangle/vertex counts, enclosed
  volume, GWN membership, and re-gate.  The change re-triangulates INSIDE a face
  (the shared canonical vertices are unchanged, so a projection-rounded crossing
  only chose a different tessellation of the SAME enclosed solid), so several
  resolving carriers are byte-DIFFERENT and others stay byte-identical while all
  enclose the identical {w>=1} solid.  Bitwise-deterministic across thread counts.
  openscad still fails closed - the swap makes the INPUT predicates exact but does
  not resolve the projection-invariant near-tangent residue.
- L2 WELD-FOLD BLIND SPOT (the R1/R2 root).  The uniform emission weld can merge
  two genuinely-distinct arrangement points that round within eps (a self-fold
  IsSelfIntersecting's shares-vertex skip is blind to) or leave a 1-ULP-distinct
  vertex unmerged (the disclosed inert zero-area triangle on rotated/irrational
  junctions).  The once-only construction rule avoids two rounded IMAGES of one
  point; D3 is the release backstop for a weld-manufactured NON-manifold outcome
  (it does not catch a manifold self-fold).  No corpus carrier; live on the
  synthetic rotated pin.
- L3 GATE CLEAN-BIAS.  IsSelfIntersecting's 2*eps shares-vertex relaxation can
  suppress a near-miss crossing in a thin band just outside the eps weld, so a
  near-degenerate self-overlap could pass the gate as clean (worse than
  fail-closed).  No corpus carrier; the gate is deliberately clean-biased.

Cluster map (root mechanism, by kind):

- CLUSTER 1 - EMISSION REPRESENTABILITY (the E1 wall above).  The resolved
  {w_S>=1} boundary reduces to an open-boundary fan of near-tangent extent that
  no double-precision per-face arrangement + winding recovers; BuildImpl
  declines.  Production-reachable (openscad); the terminal behind the kernel
  tripwire.
- CLUSTER 2 - SEAM SUB-FACE EXACT RESOLVABILITY (A1 / R2).  RecordSeams /
  EmitSeamedFace at a degenerate incidence (the coplanar/transversal
  entanglement) or an unbuilt >2-sheet radial junction.
- CLUSTER 3 - SoS / EXACT-TIE + the unbuilt >2-sheet radial branch (R2 / R3).
  The winding-probe filter-precision residue is CLOSED (the classify escalates
  each filter-0 term to Orient3DExactSign); the radial branch is proven sound and
  fires zero times, gated behind the vendor decision.
- CLUSTER 4 - NEAR-COPLANAR PLANARIZE GUARD / SNAP (S3 / S4).  Curved-chain and
  inconsistent-snap refusals - correct, decision-complete, not unbuilt axes.
- CLUSTER 5 - WELD / GATE BLIND SPOTS (L2 / L3, R1).  The weld-fold and
  clean-bias blind spots above; no corpus carrier, live on the synthetic pin.
- CLUSTER 6 - DEGENERATE-INPUT / DEFENSIVE GUARDS (D1 / D2 / E2 / D3).
  Terminal-correct contract edges, not gaps.

Closure = the plane-based-representation escalation (exact symbolic vertices ->
the exact constructed-crossing 2D arrangement of E1's spec), owner-gated behind
the kernel tripwire (vendor Shewchuk, or the one {f,g,h} predicate FORM above).
Everything not on that wall is decision-correct (Cluster 4), a defensive
contract edge (Cluster 6), or an off-corpus backstop.  No silent wrong resolve
at any arm.
