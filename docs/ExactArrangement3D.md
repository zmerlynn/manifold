# Exact-Arithmetic Arrangement for the Hard Carriers (design)

HISTORICAL RECORD. The CURRENT design is docs/Regularize3D.md (RemoveOverlaps3D as a
regularization operator: decompose by connectivity, per-component gate, run candidate
B per dirty component, re-gate, compose - never fuse). This document is the campaign
that reached it, kept for the evidence and the kill tables; read Regularize3D.md for
the design as it stands.

Status: DESIGN, load-bearing probe validated by two independent algorithms,
round-1 adversarial review folded (rev 2), kernel-avoidance adjudication folded
(rev 3), round-2 kernel-avoidance crucible (adversarial verification + S1-S4)
folded (rev 4). Candidate B is now specified, precision-cleared, and mechanism
paper-verified; its residual has narrowed to one 1D probe. The open work is
NAMED, not hand-waved: a per-region localizer, a selective weld, and
clean-boundary growth measurement (all stage-0, measurable without the exact
kernel). It is grounded on three converged results (do not re-derive them here):

  1. The hard carriers' correct topology is a PURE FUNCTION of the input soup:
     the oracle union is {p : w_S(p) >= 1}, w_S the signed winding of the
     oriented input faces. No ambiguity witness exists
     (.claude/lane-reports/a0-verify-witness-...).
  2. Every eps-quantized channel is DEGENERATE on those carriers (V4ImplPlan.md
     stage-A0 RECORDED RESULT + verification: incidence constant, fold bijective,
     positions non-monotone).
  3. Therefore THE BARRIER IS PRECISION OF DECISIONS, NOT INFORMATION.

"Exact arithmetic" here means EXACT-VALUED DECISIONS, not infinite-precision
storage: predicates on INPUT data (doubles are exact rationals; orient3d and
winding are evaluable exactly by adaptive precision or rationals), with
constructions either avoided (decisions re-derived from input predicates) or
represented exactly where unavoidable. Output COORDINATES may stay double /
eps-noisy; the TOPOLOGY must be exact. eps remains the honest weld radius for
constructed output geometry.

This doc is the exact-decision counterpart to docs/MaintainedEmission3D.md (which
converged on variant iii = the RSI-#3 arrangement completion at EPS precision,
and named exactness as the upgrade path). It re-adjudicates the one FP-premised
RSI-#3 kill under exactness, and adjudicates where exact-valued decisions enter
the pipeline.

Round-1 review lanes folded here (exact numbers live in the notebooks; this doc
carries magnitudes and cites):
  - probe : .claude/lane-reports/v5-verify-probe-1783913337.md (probe reproduction,
            per-carrier pivot, integer shells, query cost)
  - c     : .claude/lane-reports/v5-verify-c-1783911848.md     (localizer, selective
            weld, downstream care, per-carrier scope, zero silent wrong-resolve)
  - bvc   : .claude/lane-reports/v5-verify-bvc-1783911600.md   (variant E kill,
            B/C merge, reopen phrasing, staged-sketch stage-0)

Round-2 kernel-avoidance crucible lanes folded (rev 4):
  - verify : .claude/lane-reports/v5b-verify-1783917600.md (adversarial audit -
             GT7081 bit-exact retraction, fold-order determinism, GT7863 exact-zeros)
  - s12    : .claude/lane-reports/v5b-s12-1783918690.md    (S1 exhaustive input
             predicates, S2 constructed-point consistency + once-only rule)
  - s34    : .claude/lane-reports/v5b-s34-1783917823.md    (S3 radial-assembly gap,
             S4 self-adapted enumeration + Havoc paper-execution)

STANDING DOC CONVENTION. An equality claim between two computed floating-point
quantities (e.g. "fold == a+b") means BIT-PATTERN identity, verified by comparing
IEEE-754 bit patterns - NEVER a fixed-precision decimal print agreeing. A
fixed-precision print ("%.10g") hid a false "bit-exact" claim THREE times on this
branch; equal decimal digits are not equal bits. Magnitudes carry "~"; exact
figures cite the notebooks named above.

The KERNEL-AVOIDANCE ADJUDICATION section below (rev 3) re-adjudicates this
design's central claim - that an exact kernel is REQUIRED - against pure and deep
REUSE of the in-tree Boolean. It is grounded in
.claude/lane-reports/v5b-design-1783914836.md and, for rev 4, in the round-2
crucible lanes (verify / s12 / s34 above). Read it first: it refines the barrier,
resolves most carriers with zero new mechanism, and relocates the kernel
requirement to a much smaller residual than rev 2 claimed. The ROUND-2 RESULTS
section (rev 4) then executes the S1-S4 staged plan and upgrades candidate B from
"plausible, live" to "specified, precision-cleared, mechanism paper-verified."


## KERNEL-AVOIDANCE ADJUDICATION (rev 3/4: can REUSE reach the same outcomes?)

The owner's question: reach v5's outcomes by REUSING existing paradigms and
machinery - do in 3D what 2D did (STRUCTURE, not a new precision mechanism) -
and where reuse cannot, draw the boundary honestly. This section supersedes the
rev-2 "an exact kernel is required" framing IN PLACE: the kernel is not required
for the winding half and is not required for the common-case arrangement; it
survives, if at all, only as a rare filtered fallback for one predicate family.

### The refined barrier (the seed, verified)

THE ORACLE IS ALREADY A DOUBLE-PRECISION PROGRAM. boolean3 (src/boolean3.cpp,
boolean_result.cpp) resolves every two-operand union in the corpus correctly
TODAY, in doubles, with NO exact arithmetic. Its consistency comes from two
places, both structural, neither precision:
  - Shadows(p,q,dir) (shared.h:121): one global symbolic-perturbation convention
    that breaks EXACT FP ties by a sign. Not higher precision - a shared tie-break.
  - Winding03 (boolean3.cpp:387-459): winding is NOT probed per vertex. Edges
    that do not cross the other surface UNITE their endpoints (DisjointSets); one
    representative per connected component is ray-cast; the value FLOODS the
    component. Two verts in a component cannot disagree - they are identical by
    construction. Decisions are COUPLED through mesh connectivity.

So exact-valued decisions are SUFFICIENT (rev 2's probe) but demonstrably NOT
NECESSARY on these inputs. The barrier statement sharpens:

  rev 2:  the barrier is PRECISION OF DECISIONS.
  rev 3:  the barrier is precision GIVEN INDEPENDENT decisions. Structural
          COUPLING is the FP-native alternative - and it is exactly the 2D law
          (order-free local decisions do not compose; the maintained/coupled
          structure is load-bearing). The 2D fix was STRUCTURE, not precision.

The reopen already reached the exact-arithmetic version of this (one convention
makes all probes agree). Rev 3 observes the FP version is the same shape with a
DIFFERENT consistency mechanism: agreement by COUPLING (one seed cell + integer
winding deltas propagated across shared faces), not agreement by exactness.

### Candidate A - decompose + Boolean fold (PURE reuse, zero new mechanism)

If the input soup decomposes by connectivity into components that are each a
valid, non-self-intersecting manifold, fold-union them through the EXISTING
Boolean = the oracle itself, as a front-end path inside RemoveOverlaps3D.

WHY THE GATE IS LOAD-BEARING (owner review): self-intersection is NOT
invalidity - a self-intersecting mesh is a valid manifold, and the Boolean
accepts it and emits another valid manifold.  What the Boolean does not do is
REGULARIZE: multi-winding regions pass through uninterpreted (a doubly-covered
lump stays doubly covered), so the output is wrong only relative to the
boundary-of-a-simple-solid reading downstream consumers assume - an invisible
mismatch, which is the original motivation of this whole campaign.  (Measured
fine print: near self-overlap zones PairUp's alternation precondition is
violated - round 3 - so assembly there can additionally mis-pair.)  The
per-component gate is a check that exists nowhere in the library today, and
the system it feeds is best read as the REGULARIZATION OPERATOR the library
never had: the map from any valid oriented soup to the boundary of {w>=1},
which is exactly candidate B's specification.

OUTPUT-CONTRACT CAVEAT (owner review): the Boolean's output is epsilon-valid -
guaranteed manifold, NOT guaranteed self-intersection-free (rounding can create
eps-scale self-crossings; accumulation across chained operations is the original
motivation for this whole campaign).  So A's fold output gets the SAME gate as
its inputs (validity + self-intersection check), and a fold output that fails
routes back through candidate B as a single dirty component - A and B form a
small fixed-point loop, not a one-way pipe.  The corpus fold outputs below were
verified valid-manifold at oracle volumes; their self-intersection audit is an
UNMEASURED cell (production gate requirement).  Epsilon-valid output remains the
library's standing contract, so the un-repaired fold is acceptable where callers
already accept it; the B pass is the strict mode.

OWNER DECISION (supersedes the fold): candidate A's Boolean fold is DROPPED.
The production usage model is RemoveOverlaps at the END of operation chains -
any wanted cross-component union has already happened upstream (the chain's
operations ARE Booleans), so what arrives is one or a few epsilon-valid
manifolds whose only unresolved defect is internal self-overlap (the pair
corpus decomposes only because the fixtures were built by composition - BR4's
fixture-artifact warning, confirmed).  A reduces to its gate: DECOMPOSE by
connectivity, TEST each component (valid + non-self-intersecting), EARLY-EXIT
clean components, route self-intersecting ones to B per component, compose
back.  This is also a SEMANTIC choice: per-component B regularizes each object
and never fuses separate overlapping components (whole-soup {w>=1} would);
fusion remains the Boolean's job - matching the touching-contacts posture.
Consequences: the corpus PAIR gates tested union semantics and become B stress
fixtures (deliberate one-soup treatment, test-only); the production shape is
the self_intersect class.  Residue kept: B's double-rounded output passes the
same gate once (rounding can self-cross at eps scale).

MEASURED (probeA/probeA2/probeGT, scratchpad; the fold is instant everywhere):

  carrier   decomposes to        candidate-A verdict
  -------   ------------------   ------------------------------------------
  Havoc     2 valid operands     RESOLVES (winding-exact vs {w_S>=1} and a+b)
  GT7863    4 valid components   RESOLVES (INTRINSIC vertex-on-face - FP path
                                  handles it; not a precision wall)
  hull      6 valid components   RESOLVES (the SWEEP pipeline fails closed here)
  Cray      operands valid       RESOLVES (fold == a+b; extreme-scale coords)
  GT7081    13 valid components  RESOLVES (fold == a+b to ~1e-10 rel, NOT
                                  bit-exact - rev-4 retraction; both print
                                  identically at 10 digits but the bits differ
                                  (verify lane); containment agrees; the SWEEP
                                  pipeline hangs/budgets here)
  Offset1   many valid shells    RESOLVES (disjoint; the SWEEP pipeline TIMES OUT)
  Offset2   many valid shells    RESOLVES (disjoint)
  Offset3/4 one clean shell      no overlap to resolve
  openscad  2 comps, w in [-1,3] INSUFFICIENT - negative-multiplicity regime
                                  (R5): a plain Add-fold overshoots; needs
                                  per-component orientation (Add vs Subtract)
  siA/siB   one self-intersecting SINGLE SHELL - INAPPLICABLE (nothing to fold);
            shell, w in [0,2]     the genuine residual

READING. Candidate A resolves EVERY corpus carrier that decomposes into
non-self-intersecting components - which is all of them except the single
self-intersecting shells (siA/siB) and the negative-multiplicity case (openscad,
fixable by orientation-aware folding). It even covers three carriers the sweep
pipeline fails-closed or times-out on. The "intrinsic" degeneracies (GT7863's
vertex exactly on a face) do NOT break the FP+Shadows Boolean.

THE HONEST COST. The pairs and multi-shell singles decompose because they were
BUILT by composing/offsetting valid parts. A genuinely self-intersecting single
shell does not decompose. So candidate A is a PRODUCT win (a real front-end that
delegates a large slice of the corpus to the machinery this campaign set out to
replace, with zero new code beyond the dispatch) and a research NON-ANSWER for
true single-shell self-overlap. That residual is candidate B's target.

FOLD DETERMINISM (rev-4 verify-lane correction). The Add-fold's output TOPOLOGY
(manifoldness, component structure) is order-STABLE, but its output VOLUME is
order-DEPENDENT at the floating-point level: reordering / shuffling the decomposed
components before the fold shifts the reported volume from a couple of ULP (Havoc)
through ~1e-5-scale (openscad), because a different fold ASSOCIATION rounds
differently. "fold == a+b" is therefore an order-contingent FP statement, not a
determinism guarantee; the topological result is what A actually promises. State
this honestly - rev 3 presented A as clean/deterministic with no such note.

Kill-table / relocation: A introduces no mechanism, re-attempts no wall-A kill,
and cannot silently wrong-resolve (it runs the existing validated Boolean or does
not apply). Its relocation risk is nil - it either resolves fully or defers to B.

### Candidate B - Boolean-paradigm dirty-core resolver (the DEEP reuse)

For the residual (siA/siB, and any flagged region), compute the local arrangement
+ winding using boolean3's OWN discipline in DOUBLES: Shadows-convention
predicates, intersections via its Kernel shapes, winding propagated through
connectivity (COUPLED) rather than by independent probes. The localizer and
selective-weld preconditions (rev 2) apply unchanged. The crux is the self case
(boolean3 assumes two distinct operands P, Q). Adjudicated by source-level
enumeration + empirical probe (probeB.py, exact rationals vs double):

WHAT TRANSFERS unchanged (geometry-only, and FP-safe by measurement):
  - Shadow01 / Kernel02 / Kernel11 / Kernel12 and the Shadows tie-break. They
    read vertex positions and normals; they are operand-agnostic.
  - The Winding03 coupling: unite across non-crossing edges, one ray-cast per
    component, flood fill. The propagated deltas are +-1 INTEGERS (the crossed
    face's orientation on a robust dominant axis), never a near-zero FP quantity.

WHAT NEEDS SELF-ADAPTATION (bounded, not a wall):
  - The intersection enumeration (collider) returns an edge's OWN incident faces
    when P == Q; it must skip self-adjacent (vertex-sharing) faces and keep only
    genuine non-adjacent crossings.
  - Winding-of-own-vertex is a 100%-boundary evaluation (every soup vertex lies
    on the soup). The Shadows expandP perturbation already resolves this per
    query; boolean3 does the SPARSE cross-operand version today, and GT7863's
    vertex-exactly-on-a-face (resolved winding-exact by candidate A) is the
    existence proof that the convention handles the on-surface incidence.

WHAT BREAKS (the keystone - a MECHANISM gap, not precision):
  - Assembly. PairUp (boolean_result.cpp:285-301) pairs an edge's crossings by
    1D position along the edge (start-end alternation). Its own comment: if the
    order is not a clean alternation "this algorithm becomes a heuristic." That
    holds for two valid operands. A general self-arrangement needs RADIAL
    ORDERING where more than two sheets meet at an edge (triple points). Generic
    2-sheet transversal self-intersection reduces to the two-operand shape and
    PairUp may serve; the >2-sheet junctions do not, and are exactly the June
    memo's standing serial-seam-matcher / triple-point kill and this doc's own
    "the radial halfedge-ordering work REMAINS." This is new mechanism, not reuse.
    [rev-4: ROUND-2 RESULTS/S3 measured this on the corpus and found ZERO genuine
    >2-sheet junctions on any failing carrier - the keystone DISSOLVES on the
    corpus (the "dense triple points" premise was a sweep-projection artifact, not
    radial multi-sheet structure). The residual narrows to whether PairUp's 1D
    start-end alternation survives the w=2 double-covered stratum (round-3 probe).]

THE EMPIRICAL CRUX - is PRECISION the barrier for the residual? Measured NO:
  - Winding SEED anchor: double +z-ray winding equals the recorded exact w_S at
    every Havoc crux point (lump, centroids, outside). The coupling's absolute
    anchor is FP-safe.
  - Coupling CONTENT: on self_intersectA (the genuine residual), double winding
    equals exact w_S at hundreds of points spanning interior AND near-surface
    (perturbed to a tiny fraction of eps off-plane), no disagreement, the
    double-covered w=2 stratum included. The winding decision is FP-safe where
    it matters.
  - Input orient3d MARGINS in the dirty cores (Havoc, self_intersectA, and
    GT7081, the worst near-tangency carrier): across many thousands of SAMPLED
    near-incidence predicates the minimum relative margin sat around 1e-10, i.e.
    orders of magnitude above the double orient3d error bound, with no sign flip
    in the sample. The near-tangency is real but the sampled INPUT predicates
    stayed FP-recoverable. [rev-4 correction, see ROUND-2 RESULTS/S1: this sample
    was NOT exhaustive and scoped to three carriers. Exhaustive enumeration DOES
    contain raw sub-bound cases and sign flips (GT7081), and GT7863's exact-zeros
    were never sampled - so "no exact-zero appeared" was scope-limited, not a
    corpus fact. The claim survives only once restricted to the DECIDING
    predicates, which S1 certifies; and two safety MECHANISMS live under
    "FP-safe" - safe-by-MARGIN for near-tangency, safe-by-CONVENTION (Shadows SoS)
    for the exact-zeros. See the per-carrier safety table below.]

So the rev-2 "barrier is precision" does NOT hold for the winding half or the
input-predicate half of the residual. The precision concern that survives is the
rev-2 Stage-3 LEVEL-2 predicate on CONSTRUCTED intersection points (bit-growth) -
which my probes did NOT exercise (they touched input predicates only), and which
is SHARED with the exact design (both must construct intersection points). The
open question is whether Shadows on constructed points stays consistent for the
self case, as it does for two operands.

B VERDICT (rev 3): NEEDS-CHANGE, live - iterate. Kernel avoidance is established
for the winding classification (coupled integer deltas + a robust seed need no
exact arithmetic) and plausible for the common-case arrangement. The two open
items are the assembly's radial ordering at multi-sheet edges (mechanism) and the
constructed-point predicates (precision, unmeasured).

B VERDICT (rev 4): UPGRADED to SPECIFIED, PRECISION-CLEARED, MECHANISM
PAPER-VERIFIED. Round-2 executed S1-S4 (see ROUND-2 RESULTS below). Precision is
cleared for the deciding predicates (S1: exhaustive sweep, min deciding margin
thousands of times the error bound, zero certified flips); the constructed-point
precision concern is dissolved for DECISIONS (S2: every combinatorial decision
restructures to a level-0 input predicate); the radial-order mechanism does not
block on the corpus (S3: zero genuine multi-sheet edges on any failing carrier);
and the enumeration + coupled winding paper-execute against Havoc's bit-for-bit
exact ground truth (S4). The residual has narrowed from the rev-3 keystone to ONE
1D question - PairUp's start-end alternation on the w=2 double-covered stratum -
plus the first empirical B-resolver fragment. Both are round-3 targets.

B VERDICT (rev 6): EMPIRICALLY VALIDATED AT FRAGMENT SCALE, KERNEL-FREE-NESS NOW
DEMONSTRATED END-TO-END. Round-3 executed both targets (see ROUND-3 RESULTS); the
alternation residual is CLOSED - the w=2 stratum is genuinely non-alternating
(measured SSEE) but the coupled winding-delta formulation absorbs it by
thresholding w_S rather than pairing, so B replaces PairUp instead of inheriting
it. The first built fragment runs end-to-end on both Havoc's composed core AND
self_intersectA (a true single-shell self-intersector), reproducing exact ground
truth. Round-3's one scope caveat - the fragment's ENUMERATION half ran in exact
rationals while only the WINDING half was filtered doubles - is CLOSED by round 4
(see ROUND-4 RESULTS): the enumeration decisions (face-pair straddle, edge-plane
crossing, edge-edge z-order) are rebuilt as LEVEL-0 orient3d on input through the
static filter, reproduce round-3's exact seam set on both cores, and are certified
to a fraction of a percent of exact-fallback (all genuine exact-ties on the
coplanar-contact carrier, zero precision sub-bound, zero sign flips; the z-order
centerpiece is fully certified everywhere). So the kernel is avoided end-to-end on
the reachable corpus. The residual is now an ENGINEERING build-out (localizer +
selective weld + cell-complex + halfedge boundary + component-local seed) plus two
specified-but-untested axes (single-global SoS for exact-zero/coplanar carriers,
subtraction/negative winding) - no mechanism gap remains.

### Candidate C - minimal old-paradigm kernels (boundary fallback, no depth)

If a worst-case sub-errbound predicate is confirmed (see the named probe in
Risks), the cheapest acceptable kernel is an ADAPTIVE orient3d: an FP filter with
an exact fallback that fires only when the filtered sign is uncertain - which, on
the measured samples, is ~never. The winding coupling needs NO kernel at all.
This RELOCATES and SHRINKS the rev-2 kernel requirement: not a general exact
arrangement kernel, but a rarely-triggered filtered fallback for the one
constructed-point predicate family. Plain big-rational evaluation in the flagged
region remains the correctness backstop / offline oracle (the B scope of rev 2).

### Per-carrier coverage map (which reuse path resolves each carrier)

  carrier    candidate A (fold)      residual path if A does not fully cover
  -------    --------------------    ---------------------------------------
  Havoc      RESOLVES                -
  GT7863     RESOLVES                -
  hull       RESOLVES                -
  Cray       RESOLVES                -
  GT7081     RESOLVES                -
  Offset1/2  RESOLVES                -
  Offset3/4  no overlap              -
  openscad   INSUFFICIENT (R5)       A + per-component orientation (Add/Sub)
  siA/siB    INAPPLICABLE            B (dirty-core resolver) - the research case

The reuse front (A) covers the entire pair/multi-shell corpus. The research
residual (B) is the single self-intersecting shell - and there, precision is
measurably NOT the barrier; mechanism (self-enumeration + radial assembly) is.

### ROUND-2 RESULTS (S1-S4 executed, rev 4)

Round-2 built the S1-S4 probes rev 3 named. All measurement is exact-rational /
certified-forward-error analysis in the crucible lanes (s12 = S1+S2, s34 = S3+S4,
verify = adversarial audit); the doc carries magnitudes, the notebooks the exact
figures. Bottom line: precision is REMOVED as B's barrier and the triple-point
mechanism gap does not materialize on the corpus. B's residual is one 1D question.

S1 - DECIDING PREDICATES ARE FP-SAFE (exhaustive, not sampled). B's resolver does
  not evaluate a monolithic orient3d; it runs a ladder of 1D Shadows comparisons
  plus interpolated operands (Shadow01 / Kernel02 / Kernel11), which S1 decomposed
  by construction level. The static Shewchuk filter (orient3d permanent-scaled
  error bound; sign CERTIFIED when |det|/errbound > 1) was proven SOUND
  exhaustively: across the full sweep (~1.6B orient2d/orient3d predicates, every
  carrier) ZERO sign flips occurred at a certified ratio, and an mpmath oracle
  agreed with exact Fractions on every genuine flip. So candidate C's adaptive
  orient3d (FP filter + exact fallback when the ratio is sub-bound) is a CORRECT
  design; its fallback fires exactly on the sub-bound + exact-zero set.

  The load-bearing distinction is DECIDING vs NON-DECIDING. boolean3 consults the
  edge-edge z-order predicate ONLY when the shadows already overlap (a genuine
  xy-crossing) and the vertex-vs-face sidedness ONLY for a vertex projecting inside
  the face. Restricting the exhaustive sweep to that DECIDING subset:
    - The GLOBAL minimum deciding margin/errbound is ~3000x (GT7081's genuine
      edge-edge crossings; its vertex-vs-face family ~1e7 x). GT7081 was the feared
      worst case and its DECIDING predicates are fully certified.
    - GT7081's RAW sweep does contain sub-bound cases and sign flips (~175k raw
      edge-edge flips) - exactly the worst case BR1 feared the sample hid - but
      EVERY one is on a NON-crossing near-coplanar edge pair or an out-of-face
      vertex: a predicate INSTANCE the resolver never evaluates as a decision. The
      near-tangent seam's ambiguity lives in geometry the resolver does not decide.
    - The one ambiguity in the whole deciding set is openscad's ~2 edge-edge
      sub-bound instances (near-coplanar genuine crossings). The doubles were in
      fact CORRECT there; the filter simply cannot certify - which sizes candidate
      C's fallback at MICRO scale (a rarely-fired exact check, not a kernel), and
      it never fired a wrong certified sign anywhere.

  MODEL CAVEAT (stated honestly). The deciding/non-deciding boundary is defined by
  boolean3's own evaluation rule (shadows-overlap gate + vertex-inside gate), which
  S1 validated as itself FP-safe where it gates (the crossing-existence orient2d is
  certified or resolves to an exact SoS tie). It is NOT defined by running the
  resolver on the residual. If a future resolver evaluates a predicate the current
  rule skips, that instance re-enters scope. The claim is "the predicates boolean3's
  rule decides on are FP-safe," grounded on the rule (itself validated FP-safe), not
  on a resolver run.

  PER-CARRIER SAFETY - two mechanisms under "FP-safe" (rev-4 correction). Rev 3
  conflated near-tangency margins with exact-zero ties. They are DIFFERENT safety
  mechanisms and both must be named:

    carrier    deciding-predicate safety       mechanism
    -------    ----------------------------    ----------------------------------
    Havoc      margins >> errbound             safe-by-MARGIN
    siA/siB    margins >> errbound (pristine)  safe-by-MARGIN
    GT7081     deciding margins ~3000x+        safe-by-MARGIN (raw flips all
                                               non-deciding)
    openscad   ~2 deciding sub-bound +         safe-by-MARGIN + C micro-fallback,
               coincident exact-zeros          and safe-by-CONVENTION (SoS)
    GT7863     coplanar exact-zeros (hundreds) safe-by-CONVENTION (Shadows SoS) -
               + a score of vertex-on-face     NOT by margin

  GT7863 is the sharp case (verify lane, exact hunt): hundreds of exact-coplanar
  orient3d zeros and a score of genuine vertex-on-face incidences. The round-1
  sweep NEVER sampled GT7863, and its centroid-grid method is STRUCTURALLY BLIND to
  coplanar exact-zeros (a grid never lands on the measure-zero coplanar set). Those
  zeros are resolved not by more bits but by the Shadows one-convention perturbation
  (an exact FP tie -> a fixed sign) - the same SoS pattern candidate A already rides
  to resolve GT7863 to a valid manifold. Safe-by-convention, not safe-by-margin;
  stated so the rev-3 "comfortable margins" framing is not misread corpus-wide.

S2 - CONSTRUCTED-POINT PRECISION DISSOLVED FOR DECISIONS. BR2 (the rev-2 level-2
  bit-growth concern on constructed intersection points) was the unmeasured axis.
  S2 measured it and reframed it:
    - 100% of B's combinatorial / topology DECISIONS restructure to a LEVEL-0
      predicate on INPUT coordinates. In particular the edge-edge z-order that
      boolean3 evaluates via interpolation is PROVEN identical to a sign on an
      input-only determinant: sign(z1 - z2) = -sign(det[Q-P, S-R, R-P]) * sign(den),
      all operands input vertices (identity checked exact on hundreds of real
      crossings, zero mismatches). boolean3's interpolation is an EVALUATION
      artifact, replaceable by an input-only orient2d/orient3d, so no topological
      decision needs a constructed operand - BR2 dissolves for decisions.
    - Constructed operands remain ONLY for (i) OUTPUT coordinates (eps-noisy,
      already permitted behind the validity gate) and (ii) intersection /
      triple-point POSITIONS. Neither is a decision.
    - Measured cost of NOT restructuring: computing the z-order on constructed
      operands compounds error (~7x margin loss on siA) and grows bits (degree-3
      -> ~500-bit), strictly worse than the level-0 form. Confirms "keep decisions
      as INPUT predicates."
    - ONCE-ONLY CONSTRUCTION RULE is now a REQUIRED design element, not advice. Two
      independent derivations of the SAME triple point (three-plane solve vs two
      line-meets-plane paths) diverge up to MILLIONS of eps in the ill-conditioned
      near-parallel tail (measured). So each intersection / triple point must be
      constructed ONCE via a single canonical derivation, stored, and referenced
      everywhere; never compare two constructions. This mirrors the winding coupling
      (compute once, flood-fill) and eliminates the constructed-vs-constructed
      comparison STRUCTURALLY.

S3 - NO GENUINE MULTI-SHEET EDGE ON THE CORPUS (the rev-3 keystone dissolves). The
  rev-3 keystone was PairUp's radial ordering at >2-sheet (triple-point) edges -
  the June serial-seam-matcher territory. S3 enumerated every face-face
  intersection LINE on Havoc and the two genuine residual shells (siA/siB), grouped
  by exact supporting line, and counted DISTINCT PLANES per line (the exact
  broadphase reproduced the pipeline's seam counts bit-for-bit, so the enumeration
  is complete, not sampled):
    - ZERO genuine multi-sheet edges on any carrier. EVERY arrangement edge is a
      TWO-DISTINCT-PLANE edge = PairUp's two-operand shape (4 half-sheets). Not one
      edge carries more than two distinct planes.
    - The lines with >2 FACES are coplanar TESSELLATION of one sheet (adjacent
      triangles of the same plane), folded into one sheet by the coplanar grouping
      BEFORE any radial decision - a decided operation (gap=0 means "same sheet,
      group it"), never an undecidable radial tie.
    - Distinct-sheet dihedral margins run ~1e-1 (siA/siB, healthily transversal
      self-crossings) down to ~1e-7 (Havoc's right operand grazing the left surface,
      the ~50-eps edge-edge contact) - all far above the double error bound. None
      sub-errbound on any carrier.

  THE "DENSE TRIPLE POINTS" PREMISE IS REFUTED, and the correction matters for the
  wall-A / June narrative: the near-coplanar density observed on siA is a
  SWEEP-PROJECTION artifact (O(seam^2) pairwise x-crossings piling into a narrow
  x-band -> ArrangementBudget), NOT radial near-tangency. Candidate B does no plane
  sweep, so that density re-appears only as many DECIDABLE 2-sheet edges, not as
  dense >2-sheet junctions. And the wall-A "bitwise-zero angular tie" was between
  EMISSION-MANUFACTURED twin cap images (one 3D junction emitted twice by two cap
  planes ~1.7 eps apart) - a twin of ONE sheet made by the sweep machinery, NOT two
  input sheets. B assembles from input intersection edges and never self-locates via
  cap-plane emission, so that tie is not in B's world (exact enumeration confirms:
  zero exactly-coplanar inter-sheet pairs, zero sub-1e-6 near-tangent self-crossings
  in siA/siB).

S4 - SOUP-VS-SOUP SPEC + PAPER-EXECUTION. S4 specified the self version of
  boolean3's pipeline on a flagged dirty core and paper-executed its winding half:
    - ENUMERATION: boolean3's AABB collider is operand-agnostic; the only
      self-adaptation is SKIP shared-vertex (self-adjacent) face pairs, keep genuine
      non-adjacent crossings. Verified: a broadphase with that one filter reproduces
      the pipeline seam counts on all three carriers.
    - WINDING: signed INTEGER face-crossing deltas - crossing an oriented face f
      along path direction t changes w by sign(dot(t, n_f)) on an input normal, plus
      one seed probe per connected cell-component. Path-independence is a COCYCLE
      invariant: for a closed oriented 2-cycle soup the signed crossing count is a
      cocycle, so deltas sum to zero around any loop and the net delta between cells
      is path-independent (this is the faithfulness gate the probe lanes assert).
    - ASSEMBLY: radial pairing sorts incident half-sheets by the input-predicate
      determinant sign(dot(e, cross(r_i, r_j))) on input face normals, then pairs by
      winding transition. S3 makes this a trivial 2-element sort at every real edge
      (no triple points) = exactly PairUp's two-operand shape corpus-wide; the
      >2-sheet branch carries the SPECIFIED determinant rule for any input that
      forces it.
    - SEED: one double-precision ray per connected cell-component (boolean3's
      Winding03 discipline unchanged); ray degeneracy is boundary-only, so cell
      seeds almost never invoke symbolic perturbation.
    PAPER-EXECUTION on Havoc's dirty core (exact + double side by side): the coupled
    deltas reproduce the recorded exact ground truth - mission_lump w_S = 2 (the lost
    lump), left_centroid = 1 (the retained control) - with all four distinct seed
    paths agreeing (cocycle holds) and double == exact at every crossing on every
    path (FP-safe). The winding half survives the paper execution.

    JUNE-KILL DIFFERENTIATION (recorded). The June serial-seam-matcher was killed
    because its triple-point rule was UNSPECIFIED (a proximity/order heuristic) and
    because independent per-cell raycasts could disagree across a shared edge. The
    S4 spec differs on both axes: the junction rule IS an input-predicate determinant
    (not a proximity guess), and the winding is COUPLED (one seed per component +
    integer deltas over ONCE-ONLY-shared constructions), the structural cure for the
    independent-raycast kill. And S3 shows the corpus never even forces the >2-sheet
    branch. Dependency stated: this is only as consistent as the once-only
    construction guarantee (S2's rule).

### Round-3 plan for B (the next round builds/probes these)

S1-S4 removed precision and the triple-point mechanism as barriers. What remains:
  (i)  THE ALTERNATION PROBE (the narrowed residual). PairUp assumes a clean 1D
       start-end ALTERNATION of an edge's crossings. With w_S reaching 2 (the
       double-covered stratum the probe confirmed) the retention multiplicity is
       > 1; AppendPartialEdges already emits abs(inclusion) copies, so multiplicity
       is representable, but whether the emitted start/end sequence stays a valid
       alternation under w=2 in the SELF case is a per-edge 1D property S3 did not
       test. This is B's one remaining open. Cheap, no kernel.
  (ii) THE FIRST EMPIRICAL FRAGMENT. Prototype the B resolver (self-adjacency-
       filtered collider + coupled integer-delta winding + once-only construction +
       the static filter with the micro exact-check fallback) on Havocglass8's dirty
       core, and check its cell classification against the recorded exact ground
       truth. This is the first BUILT fragment (S1-S4 were analytical); it converts
       the paper-execution into running code on the smallest real core.
  Exit: either a doubles-only dirty-core resolver matching exact w_S on the Havoc
  core end-to-end (kernel fully avoided on the reachable corpus), or a precise BREAK
  at the alternation probe naming the one 1D junction class doubles cannot serve.

REQUIRED ADDITIONS for a B build (the design elements S1-S4 pinned):
  - a SINGLE GLOBAL SoS convention, reusing the Shadows(p,q,dir) pattern lifted to
    orient2d/orient3d = 0, governing EVERY predicate family (this is R4; GT7863's
    and openscad's exact-zeros depend on it);
  - ONCE-ONLY shared constructions for every intersection / triple point (S2);
  - the static Shewchuk error filter on orient2d/orient3d (S1, proven sound);
  - the micro exact-check fallback for the ~2 openscad deciding sub-bound instances
    (candidate C at micro scale - a rarely-fired exact check, not a kernel).

### ROUND-3 RESULTS: alternation + empirical fragment (rev 5)

Round-3 executed the two items the plan named (probe lane v5b-r3;
scratchpad/alternation_probe.py + fragment_havoc.py + siA_fragment.py, all
exact-rational with a double mirror). Both resolve in B's favor. Bottom line: the
alternation residual is CLOSED (the coupled winding-delta formulation replaces
PairUp rather than inheriting it), and the first BUILT fragment runs end-to-end -
matching exact ground truth with the kernel avoided - on BOTH a composed dirty
core (Havoc) and a genuine single-shell self-intersector (siA). B moves from
paper-verified to EMPIRICALLY VALIDATED AT FRAGMENT SCALE.

R3-i - THE ALTERNATION PROBE: alternation FAILS on the w=2 stratum, and the
  winding-delta formulation ABSORBS it. PairUp (boolean_result.cpp:285-301) tags
  each retained crossing start/end by the LOCAL SIGN of its winding inclusion
  (AddNewEdgeVerts sets direction = inclusion<0 and pushes abs(inclusion) copies)
  and pairs by 1D position, assuming a start-end-start-end order. Cast generic
  transversals through Havoc's w=2 lump (mission_lump, right_centroid) and read
  the exact per-cell soup winding: generic lines read 0 -> 1 -> 2 -> 1 -> 0
  (a small fraction of directions read the longer 0 -> 1 -> 2 -> 1 -> 2 -> 1 -> 0
  profile - verification round), whose
  local-sign token order is SSEE (or SSESEE) - NOT the assumed alternation. PairUp's stated
  precondition is literally violated on the double-covered stratum, exactly as
  BR-a feared. But:
    - The coupled winding-delta assembly (S4c) does not pair: it classifies each
      cell by the integer w_S and emits the boundary of the solid {w_S>=1}. On the
      double-cover the interior 1|2 and 2|1 sheets are w>=1 on both sides, so they
      are NOT boundaries and are DISCARDED; only the outer 0|1 and 1|0 sheets are
      emitted. This threshold read is indifferent to token order - no alternation
      is invoked. (Verified against the spec text: "retain the sector whose
      w_S >= 1.") A naive PairUp self-reuse instead keeps the interior sheets as
      overlapping doubled edges - two spurious sheets per w=2 line - the wrong
      output for a clean self-removal boundary.
    - WHY the two cases differ: two-operand PairUp works because each operand's
      edges carry inclusion w.r.t. the OTHER operand in {0,1} (clean SESE), and
      w=2 arises only on composition where the inclusion filter has already
      dropped the interior. The SELF case has no second operand - the natural
      field is the TOTAL soup winding, which reaches 2 in ONE field - so
      local-sign pairing over-emits and B must threshold instead.
    - On the union corpus the winding stays non-negative (measured: min 0 on every
      line), so PairUp's abs(inclusion)-copy machinery still yields a
      manifold-but-doubled result (a ballot/Dyck generalization of alternation,
      not a crash). The assumption's HARD break is NEGATIVE winding (subtraction;
      openscad's soup winding reaches -1, BR-d), where the ballot condition fails
      and local-sign pairing mislabels. The threshold read is indifferent there
      too. Double == exact on every cell (winding is the validated integer field).
  ANSWER: B does not inherit PairUp; it REPLACES the 1D position-pairing with a
  winding threshold. The "1D alternation" open is a property of the retired
  mechanism, not of B. BR-a is closed.

R3-ii - THE FIRST EMPIRICAL FRAGMENT: the S4 resolver built and run end-to-end.
  All six spec elements wired with their asserts ACTIVE.  SCOPE CAVEAT
  (verification round; RESOLVED in ROUND-4 RESULTS below): the
  level-0-plus-static-filter discipline was demonstrated for the WINDING half
  only; the fragment's ENUMERATION half ran in exact rational arithmetic
  throughout (more exact predicates than the filtered winding side, including the
  constructed-point overlap decisions S2 proved restructurable to level-0 but
  which this fragment computed exactly instead).  B's correctness is demonstrated
  end-to-end; its kernel-free-ness was demonstrated for winding and rested on S2's
  proven identities - not yet exercised - for enumeration.  Round 4 rebuilt the
  enumeration on those level-0 identities through the filter and closed the gap.
    On HAVOC's dirty core:
      - ENUMERATION (self-adjacency skip) reproduces S3's self-crossing seam count
        bit-for-bit, all cross-operand; every recorded exact w_S (the double-
        covered lump = 2, the retained control = 1, exterior = 0) is reproduced.
      - COUPLED WINDING is cocycle-unanimous across several independent exterior
        seed rays (path-independence violations: zero) and the double seed ray
        equals exact (FP-safety violations: zero).
      - RADIAL: every arrangement line is exactly two distinct planes; the
        >2-sheet branch fires ZERO times. ONCE-ONLY construction: zero double-
        constructions (no triple points to stress).
      - STATIC FILTER: EVERY level-0 WINDING orient3d certified (ratio>1); the
        exact fallback NEVER fired; min certified ratio ~5e2; zero
        certified-wrong signs. The winding half is safe-by-margin; the
        enumeration half's filtered form is the round-4 demonstration (see the
        scope caveat above).
    On self_intersectA (a TRUE single-shell self-intersector, the genuine residual):
      - ENUMERATION reproduces S3 bit-for-bit: thousands of arrangement lines, ALL
        exactly two distinct planes; genuine multi-sheet edges = ZERO (the
        >2-sheet radial branch never fires on the single shell either).
      - The coupled winding classifies the FIRST w=2 SELF-OVERLAP stratum on a
        genuine single shell (profile 0 -> 1 -> 0 -> 1 -> 2 -> 1 -> 0), cocycle-
        unanimous and double == exact.
      - STATIC FILTER: the winding rays touch near-grazing configurations - a
        MICRO fraction (order 1e-4) fell sub-bound and one double sign flipped -
        ALL caught by the filter and routed to exact; the classification is
        correct and the double-only winding also matched (the grazing predicate is
        non-deciding for the integer winding). So candidate C's micro exact-check
        is genuinely EXERCISED and SOUND on siA's winding rays (not dead code as on
        Havoc), never a wrong certified sign (S1's exhaustive soundness holds).
        Part of the grazing is a probe artifact of a far seed; a component-local
        seed (Winding03 discipline) shrinks it - a sizing input, not a wall.
  RESULT: the fragment converts S4's paper-execution into running code and every
  predicted number holds on both carriers. The kernel is avoided on the reachable
  corpus; the exact fallback is either dead (Havoc) or micro-scale and sound (siA).

B VERDICT (rev 5): EMPIRICALLY VALIDATED AT FRAGMENT SCALE. Enumeration, once-only
construction, coupled integer-delta winding + live cocycle, double seed ray +
FP-safety, the radial 2-sheet rule, and the static filter + micro exact fallback
all RUN and match exact ground truth on a composed dirty core AND a genuine
single-shell self-intersector. The residual has moved off MECHANISM entirely
(triple points dissolved in S3; alternation absorbed in R3-i) onto an ENGINEERING
build-out plus two specified-but-untested axes. No new wall.

NEXT-STEP ADJUDICATION for the owner (what remains between the fragment and a
production resolver, sized honestly):
  DEMONSTRATED in round 4 (was unbuilt at rev 5):
    - FILTERED ENUMERATION. The level-0-restructured + static-filtered enumeration
      (S2's proven identities, the edge-edge z-order the centerpiece) is now BUILT
      and RUN end-to-end on both cores, reproducing round-3's exact seam set
      seam-for-seam and certified to a fraction-of-a-percent exact-fallback (all
      genuine exact-ties, zero precision sub-bound, zero flips). See ROUND-4.
  UNBUILT - the v5 preconditions (rev-2), the bulk of the resolver:
    - LOCALIZER. The fragment is O(ntri) per winding query (whole soup); production
      must bound the query to the flagged dirty submesh. The B-scales-with-input
      vs C-bounds-to-submesh cost adjudication is now concrete (cost is real).
    - SELECTIVE WELD. Identify the flagged dirty region + weld the clean boundary.
    - CELL COMPLEX + HALFEDGE {w_S>=1} BOUNDARY EMISSION. The fragment did point
      classification of the recorded cells; production must build the arrangement
      halfedge structure and emit the retained-solid boundary. S3/S4 spec it and the
      fragment confirmed it is 2-sheet-only everywhere, but the halfedge build + boundary
      extraction is the largest unbuilt piece.
    - COMPONENT-LOCAL SEED POLICY (Winding03's unite-non-crossing-edges +
      one-ray-per-component), with a direction policy that avoids the far-seed
      near-grazing the fragment measured on siA. Micro cost, but must be built.
  SPECIFIED but UNTESTED at fragment scale (named opens, not new walls):
    - SINGLE GLOBAL SoS (R4). The round-3 winding rays hit no exact-zeros;
      round-4's filtered ENUMERATION hit only SAME-operand ones (final
      verification corrected the attribution: all are intra-mesh vertex-on-plane
      incidences on NON-CROSSING pairs, LL/RR, zero cross-operand), handled by
      exact fallback (sign 0). The cross-operand vertex-on-face family that
      GT7863/openscad need therefore remains fully UNEXERCISED at fragment
      scale. A single global SoS is the production
      requirement to decide those boundary edges deterministically; GT7863's
      coplanar zeros and openscad's coincident verts REQUIRE it too. Now partly
      EXERCISED (enumeration), still untested as a unified convention.
    - SUBTRACTION / NEGATIVE WINDING (BR-d). The fragment is all-union (w>=0); the
      alternation probe pinpoints negative winding as PairUp's genuine break and
      the point B's threshold read must absorb - untested on a subtraction carrier.
    - TRIPLE-POINT ONCE-ONLY (BR-b/BR-c). Zero triple points on the corpus, so
      S2c's two-path divergence is unstressed; a synthetic co-axial-3-face carrier
      is the attack that would exercise it.

### ROUND-4 RESULTS: filtered enumeration - the scope caveat resolved (rev 6)

Round 4 (probe lane v5b-r4; scratchpad/v5b_r4_havoc.py + v5b_r4_siA.py +
v5b_r4_perm.py) closed the R3-ii scope caveat: the ENUMERATION half, which ran in
exact rationals in round 3, is rebuilt on the S2 LEVEL-0 restructuring through the
S1 static filter (exact-Fraction fallback armed and counted per class), rerun
end-to-end on BOTH cores, and graded seam-for-seam against round-3's exact
enumeration and cell-for-cell against recorded exact ground truth.

THE RESTRUCTURING. Every enumeration decision is expressed as orient3d on INPUT
coordinates, decomposed into the three decision classes the plan named:
  - D1 FACE-PAIR OVERLAP (plane straddle): orient3d(face, vertex) - which side of
    one triangle's plane each of the other's vertices lies on.
  - D2 EDGE-PLANE CROSSING EXISTENCE: a read of the D1 signs (opposite signs =>
    the edge crosses the plane); no new predicate.
  - D3 SEGMENT-vs-TRIANGLE / EDGE-EDGE Z-ORDER (the P4 centerpiece): a genuine
    transversal crossing is "some edge of one triangle pierces the interior of the
    other", decided WITHOUT constructing the plane-plane line by orient3d(edge,
    triangle-edge) x3 (the edge-edge determinant class, det[u-v, ...] on four input
    vertices). This is the same orient3d shape B already uses+trusts for the winding
    ray crossing, and is provably equivalent to the exact interval-overlap the
    round-3 enumeration used (an overlap endpoint is always an interior pierce of
    one triangle by an edge of the other). It sidesteps the Guigue-Devillers
    canonical permutation while remaining a level-0-orient3d theorem.
The filtered enumeration's DECISIONS construct nothing (level-0 on input); only the
output segment POSITIONS are constructed, in double, once each (once-only,
zero double-constructions) - decisions level-0, geometry double, exactly the design
split S2 prescribed.

THE GRADE (both cores). Three genuine-crossing sets - the round-3 exact
interval-overlap oracle, the exact pierce formulation, and the level-0 FILTERED
pierce - coincide PAIR-FOR-PAIR on both Havoc and self_intersectA (formulation-
equivalence, filter-fidelity, and end-to-end all hold; the round-3 seam sets are
reproduced exactly, all cross-operand on Havoc, all clean 2-sheet self-crossings on
siA; the >2-sheet radial branch fires zero times). The cell/winding classification
still matches recorded exact ground truth (Havoc's four recorded cells; siA's w=2
self-overlap stratum, cocycle-unanimous, double == exact). Input-permutation
invariance re-checked and stronger than round 3: under triangle-order shuffle AND
within-triangle vertex rotation+reflection (which flips every face normal and
reorders the orient3d arguments) the geometric seam set and the filter's zero-flip
soundness are unchanged.

THE FALLBACK RATE (measured honestly, the caveat's alternative resolution). Across
the enumeration deciding predicates of both cores, the static filter CERTIFIES
about 99.8 percent. The centerpiece edge-edge z-order (D3) is 100 percent CERTIFIED
on BOTH cores - zero fallback, with large margins (min ratio ~3e2 on Havoc, ~7e6 on
siA) - so the kernel is demonstrably avoided for the z-order decision. The only
enumeration fallback is a fraction of a percent of Havoc's D1 straddle predicates,
and every one of them is an EXACT-ZERO (a level-0 predicate evaluating to exactly
0; final verification corrected the attribution: all are SAME-operand intra-mesh
vertex-on-plane incidences on non-crossing pairs - LL/RR, none cross-operand):
ZERO precision sub-bound, ZERO sign flips. siA's enumeration (an organic single shell, no
coplanar-operand degeneracy) is 100 percent certified, zero fallback, the exact
kernel dead code. No enumeration decision was found that cannot be expressed
level-0 - the FINDING class is empty.

READING. The exact-zeros are handled by exact fallback (returning sign 0) plus a
cleanly-piercing sibling edge, which keeps the filtered set exactly equal to the
oracle; they are the point where a SINGLE GLOBAL SoS (R4, in the specified-untested
list) is the production requirement to decide the boundary edges deterministically
rather than leaning on a sibling. That is a named, already-listed open, not a new
wall and not a level-0-expressibility failure. RESULT: B's kernel-free-ness is now
DEMONSTRATED END-TO-END on the reachable corpus - enumeration and winding both -
with the exact fallback either dead (the z-order everywhere; all of siA) or firing
only on genuine exact ties (Havoc's same-operand incidences), never on precision
and never producing a certified-wrong sign (final verification added a
certified-sign audit over every certified predicate: zero wrong; and filter
tightening up to 1e12x moves no classification). The remaining production work is
unchanged and
re-sized in the NEXT-STEP list above: the FILTERED ENUMERATION line moves from
UNBUILT to DEMONSTRATED; the bulk (localizer + selective weld + cell-complex + halfedge boundary +
boundary emission + component-local seed) and the two specified axes (a unified
global SoS - now partly exercised on the enumeration - and subtraction/negative
winding) remain.

### Risks refreshed to B's current opens (rev 4)

Round-2 closed the rev-3 attack surface: BR1 (input margins) is now exhaustive and
resolved to the deciding set; BR2 (constructed-point precision) is dissolved for
decisions; BR3 (dense triple points) is refuted on the corpus. B's current opens:

  BR-a. ALTERNATION UNDER MULTIPLICITY (the narrowed residual) - RESOLVED by
       round-3 probe (i), see ROUND-3 RESULTS/R3-i. PairUp's clean 1D start-end
       alternation is genuinely VIOLATED on the w=2 stratum (measured SSEE), but B
       does not pair - it thresholds the coupled w_S and emits the {w_S>=1}
       boundary, which is indifferent to token order. So B REPLACES PairUp rather
       than needing a 1D fix; alternation is a property of the retired mechanism.
       The residue folds into BR-d: the ballot condition (and thus any naive
       PairUp reuse) breaks only at NEGATIVE winding (subtraction), which the
       threshold read also absorbs but is untested on a subtraction carrier.
  BR-b. CORPUS-BOUNDED, NOT GENERAL. S3's "zero genuine multi-sheet edges" is a
       measurement of the three residual carriers, not a proof that no
       self-intersecting shell can carry a genuine triple point (a shell with three
       faces exactly co-axial would). The specified radial rule covers the >2-sheet
       branch; it is just never exercised on THIS corpus. Attack: a synthetic shell
       that forces the branch.
  BR-c. ONCE-ONLY IS LOAD-BEARING. The whole spec's consistency depends on the
       once-only-construction guarantee (S2). A build that constructs a triple point
       twice reintroduces the millions-of-eps divergence. It is a required design
       element, so the risk is an IMPLEMENTATION discipline to enforce, not an open
       question - but it must be enforced, not assumed.
  BR-d. SUBTRACTION / NEGATIVE-MULT (R5, still open). Every measurement is
       union-composed (all +1 multiplicity). openscad's soup winding reaches -1
       (verified - a local sheet reversal, not an inward shell; both decomposed
       components carry positive signed volume), so the {w_S>=1} rule and the winding
       strata under subtraction are UNMEASURED for both A and B. Attack: a
       subtraction carrier with an inward shell.
  BR-e. RAY / SoS CONVENTION UNIFORMITY (R4 residue). The winding validation used a
       shared +z projection axis on both the double and exact sides, so a shared
       axis-aligned degeneracy would be invisible; the anatomy lane's multi-direction
       unanimity partially covers this. And the single-global-SoS requirement is a
       correctness constraint: a convention consistent within orient3d but different
       for the winding ray reintroduces the cross-boundary disagreement kill.

  Candidate A's fixture-artifact / non-decomposability-dispatch risk (rev-3 BR4) is
  unchanged and lives in the candidate-A section; it is not a B open.


## The load-bearing probe (run first; it decides whether the arc is live)

Question: is the material the escape kill LOST present in the exact input, or is
the input genuinely ambiguous there? If w_S = 0 at the lost lump, the whole arc
is dead. Exact rational winding on Havocglass8's composed soup S = faces(left) ++
faces(right):

  - Faithfulness gates PASS. Both operands are closed oriented 2-manifolds
    (every directed edge paired) with positive signed volume, so orientation is
    preserved and OBJ-direct winding is faithful to the pipeline's soup (winding
    is a topological invariant of the oriented surface; ComposeImpl/GetMeshGL64
    merges nothing here - no intra-operand near-duplicate verts).

  - At the escape kill's lost-lump coordinate: w_left = 1, w_right = 1, w_S = 2
    (inside both operands). MATERIAL. The oracle keeps it.
  - Correctly-retained control region: w_S = 1. MATERIAL, kept.
  - Outside control: w_S = 0.
  - The right operand pokes a large wedge OUT of the left: most of its vertices
    carry right-only material (w_left = 0, w_right = 1, w_S = 1). The lost lump
    is a slice of this wedge.

INDEPENDENTLY CONFIRMED (round-1, probe lane). Two genuinely different algorithms
reproduced every value bit-for-bit against the original segment-cast driver:
exact solid-angle summation (Van Oosterom-Strackee, high-precision), and exact
rational semi-infinite ray casting over nine independent generic directions
(unanimity across directions is an internal ray-invariance proof). Margins are
astronomically clean (solid-angle residuals from the nearest integer are ~1e-100;
the rays are unanimous). INTERIORITY is exact: these are not boundary-adjacent
points - the lump sits millions of eps deep inside both surfaces, the control
centroids hundreds of millions of eps deep. "Strictly inside" is an exact, robust
statement here, not a sampling artifact.

VERDICT: the theory holds emphatically. The lost region and the correctly-
retained regions carry the SAME KIND of w_S value (both >= 1, both kept). eps
cannot tell them apart from quantized emitted positions; exact w_S treats them
identically. The lost lump is a pure TOPOLOGY-SELECTION error at eps, which exact
w_S >= 1 decides trivially. The decision A0 proved unreachable at eps is made
correctly from input data alone. THE ARC IS LIVE.

INTEGER SHELLS - caveat (probe lane, ATTACK 3). Sampling exact w_S up a column
through the near-coincident sheets stratifies the region into clean integer
shells, and the thin right-only sliver survives - the architectural conclusion is
correct. But the JUSTIFICATION is exact evaluation, NOT the sampling. The shell
sequence is location-dependent (0 -> 1 -> 0, 2 -> 1 -> 0, and 0 -> 1 -> 2 all
occur along different columns) with thicknesses from tens of eps (the edge-edge
surface minimum) to over ten thousand eps. Critically, at the mission-lump column
a shell hundreds of eps thick was SKIPPED by both a 400-sample and a 700-sample
scan (both reported a clean 0 -> 2 jump); only bisecting the two sheet crossings
separately resolved it. A thin shell can hide between samples. So the manifold-by-
construction claim rests on the winding being a clean integer OFF-surface
(verified exactly, margins ~1e-80..1e-100), not on any stratification sample.

QUERY COST - measured (probe lane, ATTACK 4). Exact w_S is O(ntri) PER QUERY. The
Python/high-precision baseline runs in tens of ms per query at Havoc scale (176
tris), which projects to seconds per query at GT7081 scale (~31k tris). Those
absolute numbers are Python-specific; a C++ adaptive-orient3d kernel (FP filter,
exact only when the sign is uncertain) is the production path and would be orders
of magnitude faster in the common case. The O(ntri) full-soup scaling is REAL and
makes the B-scales-with-input / C-bounds-to-submesh adjudication concrete. Ray
degeneracy (a grazing ray needing symbolic perturbation) is BOUNDARY-ONLY: zero
occurrences in thousands of random interior queries, but 100% on points placed
exactly on a face. Because the design classifies CELLS by interior points, the
symbolic-perturbation machinery is almost never invoked at runtime - the SoS load
is light (R4's correctness concern still stands; it is just not a hotspot).


## The pivot, re-scoped per carrier (was: emission-manufactured, not intrinsic)

The rev-1 pivot claimed the hard carriers' near-degeneracy is a property of the
eps sweep-emission MACHINERY, not of the input, and offered that as the general
structural argument for an exact arrangement. Round-1 direct measurement of the
minimum INPUT feature separation per carrier REFUTES the universal claim. The
pivot is HAVOC-ONLY:

  carrier    input degeneracy                                  verdict
  -------    ----------------------------------------------    ------------
  Havoc      no sub-eps input feature; surface min ~50 eps     MANUFACTURED
             (edge-edge), so the ~1.7-eps emitted twin is
             genuinely made by emission
  GT7863     a vertex exactly ON a face (zero separation,      INTRINSIC
             exact coplanar + inside)
  Offset1    exactly-coincident vertices (dx = dy = dz = 0)    INTRINSIC
  openscad   exactly-coincident vertex pairs                   INTRINSIC
  GT7081     over a million near-tangent seam contacts +       INTRINSIC
             hundreds of near-coplanar overlapping faces       (surface
             (surface near-tangency)                            near-tangency)

For HAVOC the pivot holds, but the rev-1 cleanliness numbers were OVERSTATED:
  - "parallel sheets ~7e-4 apart (tens of thousands of eps)" was the VERTEX-vertex
    gap. The true SURFACE-surface minimum is ~50 eps (edge-edge) - hundreds of
    times closer. Emission amplifies ~50 eps into the ~1.7-eps twin (a modest
    factor), which is plausible; the "tens of thousands of eps" headline was a
    vertex artifact.
  - "right verts >= 23 eps from the nearest left face plane" was a weak proxy:
    plane distance, not face distance (a vert can sit tens of eps from a far
    triangle's infinite plane while being tens of thousands of eps from any actual
    face). The 23-eps figure is not evidence of near-tangency and is dropped.

For Havoc, exact arithmetic still sidesteps the amplification: it resolves the
~50-eps input offset with certainty, never manufactures the ~1.7-eps twin, and so
never faces the sub-eps collapse decision. Both sheets are kept distinct, the
right-only sliver (the lump) survives, the {w_S >= 1} boundary is the single outer
sheet. The eps-collapse of the manufactured twin is what deletes it.

ARCHITECTURAL CONSEQUENCE (the narrowing STRENGTHENS C). The four intrinsic
carriers carry exact-coincident vertices, a vertex exactly on a face, or genuine
sub-eps surface contacts AT THE INPUT. That is precisely what a local exact
arrangement handles and what no placement shortcut ever could - you cannot "place"
your way out of two vertices that are exactly equal, or a vertex that is exactly
on a face. So the pivot's Havoc-only rationale does NOT weaken the arc; it removes
the emission-manufactured hope as a GENERAL argument (which is what kills variant
E below) and leaves the local exact arrangement as the only mechanism that covers
both regimes - manufactured (Havoc) and intrinsic (the other four).


## The reopen (the one sanctioned re-adjudication)

The June RSI-#3 memo killed option 3 - "complete the arrangement + classify each
sub-face independently by a 3D interior-point winding probe = manifold-by-
construction" - as FALSE: independent FP raycasts from different origins have no
cross-origin symbolic-perturbation coupling, so adjacent sub-faces across a
shared edge can disagree, leaving an unpaired halfedge. The kill is CORRECT under
FP.

Under EXACT evaluation with ONE global symbolic-perturbation convention the
premise is removed. All probes evaluate the SAME exact function w_S:
  - Each cell's classification is a well-defined integer; two adjacent cells
    cannot disagree about w_S. Exactness provides CONSISTENT CLASSIFICATION.
  - Around every arrangement edge the winding jumps sum to zero exactly
    (conservation of the signed multiplicities of the input faces meeting there),
    so the number of boundary sectors around each edge is EVEN. Exactness
    provides PARITY.

What exactness does NOT provide (rev-1 overclaim, corrected). The rev-1 text said
the halfedges "pair by construction... the cross-edge coupling the halfedge structure provided
for free, now provided by exactness." That oversells it. The June memo's durable
result is that manifold output needs a GLOBAL halfedge structure with RADIAL ORDERING of the
face-sides around each arrangement edge. Exactness gives consistent classification
and parity; it does NOT give the radial ORDER itself. You still must sort the
face-sides cyclically around each edge to PAIR the in/out transitions and to
DETECT genuine non-manifold edges (more than two boundary sectors). The radial-
ordering halfedge work REMAINS; exactness makes its predicates tie-free, and that
tie-freedom is exactly what removes the June FP-disagreement kill.

REOPEN VERDICT: the FP-disagreement kill does NOT survive exact winding with a
global convention. It does not resurrect a shortcut. "Manifold by construction"
still requires the COMPLETE EXACT ARRANGEMENT (all intersection curves
constructed, all faces split, radial order sorted). Corrected option-3 therefore
collapses into a COMPLETE exact arrangement - global (variant B) or local (variant
C) depending on whether a clean local boundary exists (see the B/C merge below).
The kill's mechanism dies; its deeper wisdom - you cannot skip the arrangement -
survives.

SCOPE. Only the FP-raycast-disagreement kill is reopened. The memo's other kills
stand and are NOT FP-disagreement kills: reusing boolean2's per-plane 2D winding
(the per-plane 2D winding is not the 3D solid winding - wrong quantity, not an FP
tie), and the serial seam-assembly matcher (unspecified triple-point rule). Every
measured kill in the campaign (the seven wall-A kills, the escape kill's lump
invariance, A0, the witness theorem, the naked snap) stands. [rev-4 refinement,
ROUND-2 RESULTS/S3-S4: the serial-seam-matcher kill stands as a kill of the
UNSPECIFIED heuristic, but the corpus residual carries ZERO genuine triple points -
its observed "dense" near-coplanar structure is a SWEEP-PROJECTION artifact, not
radial multi-sheet junctions - and the S4 spec replaces the matcher's unspecified
rule with an input-predicate determinant. So the kill's mechanism (an unprincipled
junction guess) dies under a specified rule the corpus never even forces.]


## In-tree machinery inventory (what an exact design can and cannot reuse)

Grepped src/: there is NO exact-predicate kernel. No orient3d/orient2d, no
adaptive/Shewchuk expansions, no interval arithmetic, no rationals, no int128.
The robustness surface is:
  - shared.h Shadows(p,q,dir): a 1D FP predicate that breaks EXACT ties by a
    `dir` sign - manifold's entire "one global perturbation convention", used
    pervasively in boolean3.cpp's sweep.
  - eps = kPrecision * bBox.Scale() welding; FP shoelace areas; FP AABB collider.

So an exact-decision architecture must INTRODUCE adaptive orient3d + exact winding
(new surface area, not a reuse). [SUPERSEDED at TEXT level by rev 3/4 - read this
sentence as the rev-2 framing it replaces. The kernel-avoidance adjudication and
ROUND-2 RESULTS show the WINDING half needs NO kernel at all (coupled integer
deltas + a robust seed), and every arrangement DECISION restructures to a level-0
input predicate the static filter certifies. An adaptive orient3d survives, if at
all, only as a MICRO exact-check fallback for the ~2 openscad deciding sub-bound
instances - not as introduced exact-winding surface area.] It CAN and MUST reuse
the Shadows PATTERN - one global symbolic-perturbation convention - which is
exactly what the reopen requires so all probes agree at shared boundaries. This is
a design constraint, not a free lunch. (Round-1 grep confirmed: no adaptive/exact
kernel in tree.)


## Variant space and verdicts

Adjudicated against: correctness composition (the relocation test - does an eps
decision downstream undo the exact one?), cost (should scale with DEGENERACY, not
input size), integration (what of v3 survives), and the reopen.

  (A) EXACT DECISION CORE in the existing eps pipeline. The eps pipeline runs as
      today; the named decisions (fan pairing, survive-vs-collapse, cancelled-edge
      agreement) are re-posed as exact input predicates consumed at the existing
      sites.
      VERDICT: KILL - relocates. The load-bearing survive-vs-collapse call is made
      at EMISSION (EmitCaps weld on constructed twin positions), DOWNSTREAM of the
      arr.vert sites where an exact answer would be injected. The eps weld undoes
      the exact answer (merge = lose lump / keep = reopen fan). To reach the
      deciding site you must exactify emission, which is the exact arrangement.
      Subsumed.

  (B/C) EXACT ARRANGEMENT WITH A SCOPE DIAL. One mechanism - exact intersection
      points, complete arrangement, radial halfedge order, exact per-cell winding
      classification, extract the {w_S >= 1} boundary, round outputs to double
      behind a validity gate - run at one of two SCOPES. These are NOT rival
      architectures (rev-1 presented them as two variants; that overstated the
      distinction). The scope dial:
        - C = flagged-region PRODUCTION mode. The eps pipeline runs; where it
          fails closed, the exact arrangement runs on ONLY the faces incident to
          the flagged region and welds to the eps exterior at a clean boundary.
          Cost scales with DEGENERACY. This is the escape kill's shape with the
          information supplied: the escape kill failed because it was positions-
          only (it guessed collapse-or-not on emitted positions, an invariant
          topology-selection error the positions could not arbitrate); exact w_S
          arbitrates it (probe: the lump is material, the sheets stratify).
        - B = whole-mesh / offline-oracle mode. The same mechanism over the whole
          soup. Cost scales with INPUT size and arrangement complexity, not with
          degeneracy. Role: the correctness backstop (run offline to validate C),
          and the unbounded LIMIT of C when no clean boundary exists.
      VERDICT: ADVANCE - front-runner is the C scope, with the B scope as oracle
      and fallback. GT7081's anatomy is delocalized (over 90% of critical-x gaps
      are sub-eps, one contiguous run of ~34k near-coincident planes packed into a
      sub-eps band), so it lands at the WHOLE-MESH end of the dial: C = B there,
      no speedup. The relocation test passes for TOPOLOGY at both scopes; the only
      eps decision is the output-coordinate snap (gated by Is2Manifold / self-
      intersection). Integration: the C scope is MAX v3 reuse (bulk runs eps); the
      B scope dies most of v3. See the C preconditions and staged sketch.

  (D) EXACT 2D SECTIONS. Per-slab sections computed as exact functions of input
      for flagged slabs.
      VERDICT: KILL - insufficient. Exactifies the in-plane arrangement but leaves
      the CROSS-PLANE emission (junction self-location across cap planes) in eps.
      The twin is a cross-plane phenomenon - two cap planes at distinct criticals
      emit one junction twice; exact-per-slab does not unify the cap images. The
      twins SURVIVE exactness-per-slab. It does not reach the deciding decision.

  (E) EXACT VERTEX PLACEMENT (no arrangement). PROPOSED-AND-KILLED (round-1, bvc
      lane). Proposal: at flagged junctions only, compute the emitted vertex
      position from INPUT data exactly, round once, place both cap images there -
      no local arrangement, no patch - so the twins never diverge because both
      derive from one exact point.
      VERDICT: KILL. Decisive argument: exact PLACEMENT of a vertex is not exact
      CLASSIFICATION of a cell. The lump's survival is a per-cell WINDING fact
      (w_S = 2, independently re-derived; column 0 -> 2 -> 0), proven INVARIANT
      across representative choice by the escape kill - an input-exact point is
      merely a THIRD representative, covered by that invariance. So on the images E
      can touch, making the two images bitwise-equal is a MERGE = the killed weld-
      bump / probe-3 collapse, which loses the lump; leaving them distinct reopens
      the fan (E built no arrangement to pair them). And on Havoc's ACTUAL blockers
      - weld images that name no canonical arr.vert (bare 2D track crossings with
      no junction identity) and far-track interpolations - placement has nothing to
      anchor to; E cannot even place them (this is R2's already-measured never-
      bound wall). It also does nothing for GT7863's bracketed-macro dead zone,
      which is a cell/winding question. E is the naked snap with cleaner
      coordinates on the images it can touch and the never-bound wall on the images
      that block. Placement is not classification. E does NOT demote C - it
      sharpens why C's local arrangement + per-cell winding is necessary.

RECOMMENDATION: build the C scope (flagged-region production: bounded exact
arrangement + clean-boundary weld). The B scope is the offline oracle that
validates it and the fallback for any region C cannot bound with a clean boundary
(percolation). C and B are the same object at two scopes; C is B made to pay only
for degeneracy.


## C preconditions (round-1 review: what must exist before build)

Round-1 (c lane) broke C on three of four attacks. The mechanism survives (exact
w_S arbitrates the lump; exact-arrangement cost is face-count-bounded corpus-wide;
no silent wrong-resolve found), but three preconditions must be built and one
scope claim corrected before C is buildable.

1. LOCALIZER (Stage 2 has no seed today). C's flag-and-grow needs to know WHERE it
   failed and which INPUT faces are incident. Neither exists. Today's flags are
   GLOBAL: NonManifoldEmission is `!IsManifold() || !Is2Manifold()` on the whole
   mesh (overlap3.cpp:1189) with no offending edge/vertex; ArrangementBudget is a
   cumulative piece counter (overlap3_sweep.cpp:246) with no region;
   SplitTouchingSheets returns a bare bool. And emitted OutTri3D is `vec3 v[3]`
   (overlap3.cpp:603) with NO provenance to input/canonical faces. So there is
   neither a SEED (which region failed) nor a MAP (emitted geometry -> input
   faces). The rev-1 line "the flagging machinery already exists as FatalReason
   returns" conflated knowing-it-failed with knowing-WHERE. Building the localizer
   = instrument emission provenance + make each gate report the offending halfedge
   + walk it back to the incident input faces. Named design work, stage-0 below.

2. SELECTIVE WELD (Stage 4 is unsound as sketched). There is ONE assembly weld:
   BuildImpl getVertIdx (overlap3.cpp:1119-1141), a uniform hash grid returning
   the min-index prior vert within Euclidean eps, through which EVERY emitted
   triangle passes - exterior AND any C patch. A faithful toy re-impl shows a
   single uniform radius CANNOT serve an exact patch:
     - two exact-distinct patch verts ~0.4 eps apart are MERGED - the arrangement
       resolved them, the weld re-manufactures the twin;
     - a patch-INTERIOR constructed vert ~0.3 eps from a shared boundary vert is
       MERGED - the interior collapses onto the boundary, exact topology corrupted.
   The two requirements collide: one radius must simultaneously MERGE patch-
   boundary onto the exterior vert (to stitch) and NOT merge sub-eps-distinct
   patch-interior verts (to preserve topology). The corpus DOES emit sub-eps-
   distinct arrangement verts (GT7081's sub-eps critical structure is direct
   evidence), so the collision is real, not hypothetical. The patch needs weld
   EXEMPTION with boundary-conformance constraints: the patch boundary adopts the
   eps exterior's vert positions (so the stitch conforms even where the exterior
   emitted constructed, eps-noisy positions along an input edge); interior verts
   are weld-exempt. A bounded local weld distinct from the global eps grid.

3. DOWNSTREAM CARE (Stage 4 walk-order). The assembly walk is: getVertIdx eps-weld
   -> degenerate/exact-dup-tri drop -> SplitTouchingSheets -> CreateHalfedges ->
   IsManifold/Is2Manifold gate -> SortGeometry/SetNormals/SetEpsilon. Two entries
   need care:
     - getVertIdx eps-weld: NEEDS-BYPASS (precondition 2). The dup-tri drop runs on
       WELDED indices, so it is coupled to the bypass.
     - SplitTouchingSheets: NEEDS-CARE. Its radial pairing returns false ->
       NonManifoldEmission on radial ties, zero-length slivers, or non-alternation.
       A legitimate exact patch with GENUINE tangent sheets (R6 touching contacts)
       is REJECTED, not composed - exact-correct patches can fail closed here.
     - The IsManifold/Is2Manifold gate is patch-safe but GLOBAL: one bad weld
       anywhere aborts the WHOLE result (the fail-closed floor is preserved, but a
       single bad weld loses the model, not just the patch).
     - SortGeometry/SetNormals/SetEpsilon are pure reorder/recompute, no re-weld -
       they genuinely compose.

4. PER-CARRIER SCOPE (correct the "localized" claim). Two locality proxies per
   carrier (dirty core = largest near-degenerate face cluster):

     carrier    fatal reason        core faces   core span (frac of diagonal)
     -------    ----------------    ----------   ----------------------------
     Havoc      NonManifoldEmit     ~8.5%        ~0.96   (model-spanning)
     GT7863     NonManifoldEmit     <1%          ~0.12
     Offset1    NonManifoldEmit     <1%          ~0.03
     openscad   NonManifoldEmit     ~2.6%        ~0.41   (model-spanning)
     GT7081     ArrangementBudget   <1%          ~0.56   (>90% sub-eps crit gaps)

   READING: COST (face count) IS bounded corpus-wide - every dirty core is under
   9% of faces, so C's exact-arrangement region does not swallow the model by face
   count. The "scales with degeneracy" cost story SURVIVES on face count. BUT the
   core is SPATIALLY MODEL-SPANNING on the single-feature carriers (Havoc ~0.96 of
   the bbox diagonal, GT7081 ~0.56, openscad ~0.41). The rev-1 "near-coincidences
   are localized, a clean boundary exists" conflated "one contact region" with
   "spatially small"; they differ. The near-degenerate contact is a model-spanning
   SHEET (few tessellated faces, huge reach), so the clean-boundary weld seam
   threads across the whole model - R2 weld exposure scales with MODEL size, not
   degeneracy. And GT7081 is delocalized at the ARRANGEMENT level (the sub-eps
   critical run), so C grows to most of the arrangement = B; today it fails
   ArrangementBudget, so C would run full B to upgrade it, at B's cost.


## Zero silent wrong-resolves (recorded finding)

Round-1 (c lane, ATTACK 4) spot-checked C's foundational premise - that the eps
pipeline never resolves WRONG silently, it only fails closed. For three RESOLVING
corpus cases (Offset2, Offset3, Cray), exact rational w_S(input soup) was compared
against the eps pipeline's output containment at sampled interior/exterior points.
Everywhere measured, they AGREE - zero disagreements. The zero-oracle-wrong record
holds on this spot check; refusal (fail-closed) is the only observed failure mode.
This supports C's premise where tested and means C only ever UPGRADES a fail-closed
region to correct, never corrects a silent error (there were none to find). SCOPE:
all three are union / self-overlap (all-positive multiplicity); the subtraction /
negative-multiplicity regime (R5) is untested here.


## Surviving architecture (C): staged sketch (rev 2)

Rebuilt with the v4-plan discipline: a stage-0 go/no-go measurable WITHOUT the
exact kernel (the second cheap decision point - the probe was the first), an
objectively-checkable gate and stub-proof pin per stage, and an explicit
stage -> v3-machinery kept/deleted map so "MAX v3 reuse" is a checked claim, not
an adjective.

STAGE 0 - GO/NO-GO WITHOUT THE KERNEL (build this before Stage 1). C's own load-
  bearing assumption is R1 (clean-boundary growth is bounded on the corpus), and
  R1 is measurable with NO exact arithmetic: run the eps pipeline, build the
  localizer (precondition 1), grow each flagged region along guard-clean input
  edges, and measure the size distribution. Also produce the selective-weld design
  note (precondition 2). This is the deciding experiment the plan must not defer.
  GATE: bounded clean-boundary growth -> build C (Stage 1). Unbounded
  (percolation) -> only the whole-mesh B scope survives, re-decide before spending
  on the kernel.
  INITIAL ANSWER (c lane, folded): face-count bounded (< 9% everywhere) = GO on
  cost; but spatially MODEL-SPANNING on Havoc/openscad and DELOCALIZED on GT7081
  (C = B there). So the gate is a PARTIAL go: build C for the face-count-bounded,
  spatially-compact carriers; dial GT7081 to the whole-mesh B scope explicitly.
  The localizer and the selective weld are the two pieces of design work this stage
  must finish before any kernel spend.
  v3 map: reuses the FatalReason returns at BuildSlabs/EmitCaps/EmitStrips as the
  failure SEED, plus the NEW provenance channel; the eps pipeline is untouched.

STAGE 1 - EXACT PREDICATE KERNEL (the reusable primitive; cost scales with
  degeneracy, not size).
  - Adaptive orient3d (FP filter first; exact expansion only when the FP
    determinant is not sign-certain). Cheap in the common non-degenerate case.
  - Exact w_S(p): segment-crossing count on INPUT vertices (orient3d signs), one
    fixed global symbolic-perturbation convention (the Shadows pattern lifted to
    3D orient) so a ray grazing a shared edge/vertex resolves identically for
    every probe. This is the load-bearing consistency guarantee of the reopen.
  PIN (stub-proof): differential sign-agreement vs a throwaway rational driver on
  the probe's query points, bit-for-bit. The probe lanes already pass a partial
  instance (two independent algorithms agree bit-for-bit).
  v3 map: NEW surface area alongside the bulk path; DELETES nothing (the tree has
  no exact kernel to reuse).

STAGE 2 - FLAG + BOUND. The eps pipeline runs unchanged. On a guard fire, the
  localizer (Stage 0) supplies the seed and the incident INPUT faces; GROW the
  region until its boundary is guard-clean: every boundary edge is an INPUT edge
  the eps pipeline resolved without a guard. The region is a submesh whose
  interface with the exterior is shared input geometry (bitwise equal on both
  sides).
  GATE: the grown boundary is entirely guard-clean input edges (checkable).
  v3 map: reuses the localizer's seed + the FatalReason returns; adds the grower.

STAGE 3 - LOCAL EXACT ARRANGEMENT. Compute the exact arrangement of the submesh:
  rational intersection points among its faces, cell decomposition, radial-order
  halfedge order around each internal edge (tie-free predicates - the reopen), classify each
  half-face by exact w_S (Stage 1), extract the oriented {w_S >= 1} boundary. The
  topology is exact; constructed intersection coordinates round to double only at
  emission, behind a LOCAL validity gate (2-manifold, no self-intersection within
  the patch).
  GATE: local 2-manifold + no self-intersection (checkable).
  v3 map: NEW.
  COST (bvc lane, folded): level-1 predicates on INPUT are a few machine words
  (a real Havoc coord is an exact m/2^k with m ~53-bit; orient3d exact fallback
  ~165-bit). Building the arrangement reaches level-2 predicates on CONSTRUCTED
  points (order intersection points along an edge; point-in-face) whose operands
  are ~165-bit rationals, so a degree-3 op yields ~500-bit numerators; feeding
  constructed coordinates into further predicates cascades (~1500-bit at level-3).
  Tractable (GMP handles 500-bit) but NAME it: keep decisions as INPUT predicates
  wherever possible so bit-growth stays bounded.

STAGE 4 - WELD. Stitch the local patch to the eps exterior along the shared
  input-edge boundary using the SELECTIVE weld (precondition 2): patch boundary
  adopts the exterior's vert positions; interior verts are weld-exempt (a bounded
  local weld, NOT the global eps grid). Handle the downstream-care entries
  (precondition 3): SplitTouchingSheets must not reject legitimate tangent patches.
  GATE: boundary vertices are INPUT vertices, shared BITWISE (checkable); no
  patch-interior vert merged.
  v3 map: reuses the assembly weld WITH a bounded local exemption, and
  SplitTouchingSheets WITH a needs-care patch for genuine tangents.

STAGE 5 - VALIDITY GATE + FALLBACK. If the welded result has unpaired halfedges or
  fails validity, fall back to the B scope on that connected component, or fail
  closed. The honest fail-closed contract is the floor; the exact resolver only
  ever UPGRADES a fail-closed region to correct.
  v3 map: reuses the global IsManifold/Is2Manifold gate as the floor.

Second-probe status: the local mechanism's CORE is validated at toy scale on the
real cluster (exact winding stratifies the sheet region cleanly; the lump
survives; the boundary is manifold-by-construction). The UNVALIDATED load-bearing
assumption is Stage 0's clean-boundary growth at corpus scale (partially answered:
face-count GO, spatial-spanning caveat; see Risks).


## Risks (the attack surface for the adversarial round)

R1. CLEAN-BOUNDARY GROWTH (the load-bearing assumption). C's cost and correctness
    both hinge on growing every flagged region to a guard-clean input-edge
    boundary. Round-1 (c lane) gave the initial measurement, now Stage 0's job:
    face-count bounded corpus-wide (< 9%), but SPATIALLY model-spanning on
    Havoc/openscad and DELOCALIZED on GT7081 (C = B). So the risk is REAL and
    partly realized: on single-feature carriers the boundary threads across the
    whole model, and on GT7081 no bounded neighborhood is clean. Open work: finish
    Stage 0's growth measurement with the localizer built.

R2. WELD RELOCATION. A boundary input edge may be crossed by a third face near-
    degenerately, so the eps exterior subdivides it at an eps-placed vertex the
    exact patch places exactly -> an unpaired halfedge at the weld. Round-1
    sharpened this: because the dirty core is spatially model-spanning, the weld
    seam is long, so R2 exposure scales with MODEL size, not degeneracy. The
    selective weld (precondition 2) is the mechanism-level mitigation (boundary
    conformance to exterior positions); the residual attack is a boundary edge with
    a hidden near-degenerate crossing.

R3. OUTPUT SNAP-ROUNDING. Rounding constructed intersection coordinates to double
    can create NEW near-coincidences at output (two exact-distinct patch vertices
    round together), silently breaking 2-manifoldness below the validity gate's
    resolution. Iterated snap-rounding may be required. ADDED (bvc lane): to bound
    bit-growth (Stage 3), the cascade may need INTERMEDIATE snapping, and
    intermediate snapping can change TOPOLOGY, not just output coordinates - a
    stronger failure than the output-only case. Attack: a patch whose exact
    vertices are closer than eps after the arrangement, so no snap is both valid
    and faithful.

R4. ONE CONVENTION ACROSS ALL PREDICATES (the subtle reopen residue). Exactness
    cures the FP-disagreement kill ONLY IF a SINGLE symbolic-perturbation
    convention governs EVERY predicate family (orient3d, the winding segment-
    crossing, the arrangement's intersection ordering). A convention consistent
    within orient3d but different for the winding ray reintroduces cross-boundary
    disagreement - the kill in disguise. Round-1: the runtime LOAD is light
    (ray degeneracy is boundary-only, ~0% on interior cell points), so this is a
    correctness constraint, not a hotspot. Attack: construct a mixed-convention gap.

R5. SUBTRACTION / NEGATIVE-MULT REGIME. Every measurement is union-composed (all
    multiplicities +1), including the zero-silent-wrong-resolve spot check.
    Subtraction inputs carry negative multiplicities; w_S can go negative; the
    {w_S >= 1} rule still holds but the strata and the winding jumps differ.
    Unmeasured. Attack: a subtraction carrier with an inward shell.

R6. COINCIDENT / COPLANAR FACES (measure-zero overlaps). Exact w_S is undefined ON
    a coincident face (measure zero); the fill needs a tie-break there. Round-1
    made this concrete, not hypothetical: GT7863 has a vertex exactly on a face,
    Offset1 and openscad have exactly-coincident verts, and the intrinsic carriers
    carry these AT THE INPUT. The tie-break must match the library's (the
    coincident-wall cleaning the witness lane's Case B relied on). Attack: coplanar
    coincident sheets where w_S jumps by 2 across a doubled wall, or an operand wall
    buried inside another.

R7. KERNEL COST + CERTIFICATION. Adaptive orient3d must be certified correct
    (Shewchuk-class); the exact arrangement's rational constructions reach ~500-bit
    at level-2 even when bounded (Stage 3 cost). Attack: a flagged region that is
    bounded by face count but internally DENSE - GT7081's sub-eps critical run
    (~34k near-coincident planes in a thin band) is the concrete instance, where
    the local arrangement is combinatorially large despite a small face count.
    [rev-4, ROUND-2 RESULTS/S3: the "dense" structure is a SIZE / budget axis, not
    a radial-mechanism axis. The near-coplanar density is a SWEEP-PROJECTION
    artifact (O(seam^2) x-crossings in a narrow band); the underlying arrangement
    edges are still 2-sheet-decidable. Candidate B does no plane sweep, so it does
    not inherit that O(seam^2) budget explosion - the density re-appears only as
    many decidable 2-sheet edges. The kernel-COST concern (bit-growth,
    certification) stands; the "combinatorially large radial junction" reading does
    not.]


## What is settled (do not re-litigate)

The probe results (now confirmed bit-for-bit by two independent algorithms) and
every measured kill across the campaign; A0; the witness theorem; the escape
kill's lump invariance. Variant E is KILLED (placement is not classification). B
and C are MERGED into one mechanism with a scope dial. The zero-silent-wrong-
resolve spot check holds on the union corpus. The v3/v4/v5 branches stay as they
are; zero-oracle-wrong; correct-or-honest. The only object reopened is the
FP-premised RSI-#3 option-3 kill, which does not survive exact winding with one
global convention - and which, corrected, is the exact arrangement this doc designs
C to pay for only where the eps pipeline fails closed.

Settled at rev 4 (round-2 crucible, do not re-litigate): candidate A resolves the
whole decomposable corpus (topology), with the fold's VOLUME order-dependent but its
topology stable; candidate B's deciding predicates are FP-safe exhaustively (static
filter proven sound over ~1.6B predicates, zero certified flips; the raw GT7081
flips are all non-deciding); B's decisions restructure 100% to level-0 input
predicates; the once-only construction rule is required; the corpus carries zero
genuine multi-sheet edges (the "dense triple points" premise refuted, a
sweep-projection artifact); B's winding half paper-executes against Havoc's exact
ground truth. B's status is SPECIFIED, PRECISION-CLEARED, MECHANISM PAPER-VERIFIED.
Its one open residual is the w=2 alternation probe (round-3). The GT7081 "bit-exact"
claim is RETRACTED (rel ~1e-10, not bit-exact).
