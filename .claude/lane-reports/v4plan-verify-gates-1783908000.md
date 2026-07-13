# v4 impl crucible round-1 ADVERSARIAL LANE: plan conformance + gate integrity

Branch explore/sweep-plane-3d-v4 @ 2d404d57. Artifact under attack:
docs/V4ImplPlan.md. Design of record: docs/MaintainedEmission3D.md. Stage-0
notebook: .claude/lane-reports/v4impl-stage0-1783901754.md.

Method: reasoning lane + targeted code reads (src/overlap3.cpp, overlap3.h,
overlap3_sweep.cpp, test/overlap3_test.cpp). rsync copy + fresh cmake configured
under /tmp/v4plan-gates (Release, MANIFOLD_TEST=ON) for optional fence confirm;
findings are code-read based, build not load-bearing.

VERDICT: NEED-CHANGE. The plan is unusually honest and the gates-first frame is
real, but four load-bearing softenings let an implementing agent declare stages
DONE while resolving nothing on the real corpus and while the identity table -
the whole point - is not proven load-bearing. All are fixable with tighter gate
wording; none is a fatal design break.

---

## ATTACK 1 - THE DONE-GATE SOFTENING (central adjudication)

### (a) What the design ACTUALLY claims - necessity vs sufficiency

The escape kill (probe 3) and convergent verdict assert NECESSITY ONLY, never
sufficiency. Quotes from MaintainedEmission3D.md:

- L471-472: "The oracle-TRUE topology is UNREACHABLE from emitted positions;
  reaching it NEEDS the upstream arr.vert identities that emission discards."
  ("needs" = necessary condition on the identities, WITHOUT them unreachable.)
- L582: "only the global identity collect CAN reach the oracle-true topology."
  ("only X can reach" = X is the unique candidate, a necessity/capability claim,
  NOT "X DOES reach" = sufficiency.)
- L596-598 (convergent verdict): "the oracle-true topology NEEDS the upstream
  identity emission throws away. Variant (iii) is the one architecture that
  CARRIES that identity." (carries the necessary info; does not claim reaches.)

So the design proves: WITHOUT the upstream identities, oracle-true is
unreachable. It does NOT prove: WITH them, oracle-true IS reachable. The plan
(L33-40) represents this correctly and explicitly: "It did NOT settle that the
identity-carrying emission reaches ORACLE-TRUE on the hardest carriers ... That
claim is the RESEARCH FRONTIER ... it is UNPROVEN on Havoc."

### (b) Reachable-in-principle? NO - only necessity is proven.

So carrier resolution cannot be pre-declared a proven stage-C/D target. But it
equally cannot be pre-declared unreachable. It is a genuinely open question.

### (c) Genuinely unproven -> plan MUST decide it EARLY. IT DOES NOT.

This is the finding. The sufficiency question (does feeding the stage-A identity
table into emission MOVE Havoc's ~336 oracle error?) is first touched only at
stage C's "corpus re-adjudicated" - AFTER the full A/B/C coordinated-re-emission
architecture is built. There is no cheap early throwaway probe (the emission-side
analog of the stage-0 skeleton) scheduled to de-risk sufficiency before the
expensive build. The stage-0 skeleton proved the table is bounded and
order-invariant UPSTREAM; it did NOT test whether consuming it CLOSES anything.
The escape kill operated positions-only DOWNSTREAM. The gap between "upstream
identity exists" and "consuming it reaches oracle-true" is UNBRIDGED by any
probe, and the plan schedules the bridge for the END.

Worse: the one gate that RESEMBLES a sufficiency test - pin 4b
(Pin_IsolatedTwin_ClosesOracleTrue, stage C) - tests a case the design ALREADY
proved closes WITHOUT the identity table. Probe 3's POSITIVE result
(MaintainedEmission3D L483-485): "An ISOLATED case-(1) twin IS boundedly CLOSABLE
by the local re-mesh: collapse + per-plane planar cap re-triangulation closes it
to a clean 2-manifold, volume-neutral." That is positions-only (killed candidate
6). So pin 4b's mutation ("disable cap-plane unification -> NonManifold") proves
only that cap-plane unification closes the twin - NOT that the identity table is
load-bearing. An implementation carrying NO identity passes pin 4b.

### VERDICT on attack 1

HONEST epistemics on the necessity/sufficiency distinction (correctly lifted from
the design), but SCOPE-SANDBAGGING IN EFFECT via two gaps: (1) branch (c)
violated - no early sufficiency experiment; the project's central unproven
premise is only tested after building A-D; (2) the sufficiency-looking gate (pin
4b) is satisfiable by identity-free machinery, so passing it certifies nothing
about the frontier.

CORRECTED GATE WORDING (add to plan before stage B/C):
- EARLY SUFFICIENCY PROBE (env-gated, throwaway, stage-0-skeleton style): bind ONE
  entangled carrier's emission to the stage-A table's canonical junction positions
  and record whether the oracle-volume error MOVES from the escape kill's ~336
  baseline. If it moves materially toward zero -> carrier resolution becomes a
  stage-C red-first TARGET. If it does NOT move (identity present, still
  oracle-wrong by ~336) -> that is THE project-defining discovery, surfaced early;
  re-scope honestly rather than building A-D to learn it.
- STRENGTHEN pin 4b: add mutation "blank the identity table -> isolated twin
  REOPENS" so the pin proves the table is load-bearing for closure, not just
  cap-plane unification; and pin the twin fixture's EXACT coordinates in the plan
  so the code author cannot detect-and-dodge (attack 2).
- PROGRESS RATCHET: name at least one currently-fail-closed entangled carrier as a
  must-MOVE target with its recorded oracle error, so "zero-oracle-wrong" cannot
  be satisfied by fail-closing everything (attack 2, pin 4c).

---

## ATTACK 2 - PIN-BY-PIN STUB ATTACK

Stub-passable pins: 4 of 5 (pins 1, 2, 4a, 4b), plus a structural loophole
shared by pins 4c/4d/5. Only pin 3 is nearly sound and even it has a hole.

PIN 1 (OrderInvariant, C1) - STUB-PASSABLE. The pin asserts the sorted junction
set is bitwise-identical under input permutation on the fixtures. But "LexMin over
v3 arr.verts" (exactly what the skeleton used) passes this WITHOUT being C1-proper
("pos is a PURE FUNCTION of the member INPUT coordinates," invariant #1). On the
small fixtures arr.verts come out bitwise under permutation, so LexMin-over-
constructed-coords passes; the pin never runs on GT7081 (plan admits too heavy),
which is exactly where arr.verts may be order-dependent (R-ORDER-UNMEASURED). The
pin tests the OUTCOME (bitwise on small fixtures), not the MECHANISM (from-input
derivation). Tighten: assert the canonical is recomputable from INPUT-vertex ids,
independent of the arr build - not merely stable on the two fixtures.

PIN 2 (OnceOnly, C2 counter) - STUB-PASSABLE. The per-junction construction
counter == junctionCount catches "moved the only canonicalization out"; it does
NOT catch "ADDED a redundant re-derivation in the consumer" - which violates C2
("canonicalized exactly once") but leaves the construction counter untouched. The
spec says "assert a per-consumer path does NOT re-derive it" but gives no
mechanism and only the move-out mutation. Tighten: a GLOBAL once-only counter on
the representative-selection PRIMITIVE, asserted == junctionCount across the whole
pipeline run (construction AND emission).

PIN 3 (Classifier_Fence) - PARTIALLY STUB-PASSABLE. Best-specified pin (two
directions, two-sided mutation) but the mutation for the LOAD-BEARING conjunct is
vacuous on the specified fixture. The fixture's distinct pair is "tens of eps
apart in x" - excluded by the dx threshold ALONE. So "drop maxTr>=maxDx" does NOT
merge it (still dx-excluded); only "widen dxThresh" reds. A dx-threshold-only stub
(no transverse test) passes pin 3. Yet section 1.2 calls the transverse conjunct
"load-bearing." Tighten: the fixture must ALSO carry a sub-eps-x, transverse-
SEPARATED pair (maxDx small, maxTr < maxDx) that an eps-ball merges but the
conjunct excludes, so "drop maxTr>=maxDx" genuinely reds.

PIN 4a (CapImage_CanonicalTransverse) - STUB-PASSABLE. Asserts the junction's two
cap images have bit-identical (y,z). A stub that snaps ALL verts to a shared grid
makes both images equal WITHOUT reading the table - passes. Mutation ("restore
self-location -> diverge") reds a self-locator but not a global-snap stub.
Tighten: assert image.yz == identityTable.junction(v).pos.yz BITWISE (white-box
table read), + mutation "blank the table -> pin cannot find the canonical."

PIN 4b (IsolatedTwin_ClosesOracleTrue) - STUB-PASSABLE two ways. (i) Fixture is
"hand-built," unspecified coordinates, author-controlled alongside the code -> a
detect-and-dodge special case (match bbox/vert count, emit hardcoded answer)
passes. (ii) Deeper: the isolated twin closes positions-only per probe 3's
positive result, so the pin's mutation certifies cap-plane unification, not the
identity table. Tighten: pin exact fixture coords in the plan; add mutation "blank
the identity table -> twin reopens." (See attack 1.)

PINS 4c / 4d / 5 - THE FAIL-CLOSED LOOPHOLE (structural, shared).
- Pin 4c accepts, for EVERY pair (Havoc, GT7863, GT7081, hull, Cray), a resolve
  matching the oracle OR a named-guard fail-closed. There is NO corpus pair that
  MUST resolve (hull/Gate4c is currently SKIPPED). So an implementation that
  regresses ALL corpus pairs to fail-closed passes pin 4c AND pin 5. The current
  corpus already fail-closes (CorpusPairGate: Havoc/7863 -> NonManifoldEmission,
  7081 -> ArrangementBudget; test comments confirm). "Zero-oracle-wrong" is
  therefore satisfiable by DOING NOTHING to the hard cases - the exact vacuous
  green the plan claims to guard against. Needs a progress ratchet (attack 1).
- Pin 5 is a REGRESSION/wrongness guard, not a progress gate - honest, but it must
  be PAIRED with a must-move target or the DONE criteria are vacuous-green-safe.
  Also: the CURRENT CorpusSingleGate (test L1619-1637) checks single-mesh resolves
  on VALIDITY ONLY (ours.Status()==NoError) - NO oracle, NO tol-invariance. So pin
  5 / pin 4d is NOT enforced by today's harness; it REQUIRES rewriting
  CorpusSingleGate. A build agent claiming "pin 5 green" without that rewrite is
  vacuous - orchestrator must verify the rewrite landed.
- Pin 4d/pin 5 TENSION: pin 5 says "any carrier that cannot be checked stays
  FAIL-CLOSED," but pin 4d lets SINGLES resolve on tol-invariance, which
  R-THINNING-ORACLE itself admits "does NOT prove it is geometrically right."
  Tol-invariance is the escape hatch that lets an oracle-unverifiable single-mesh
  resolve count as a pass; a deterministic-wrong stub passes it. Corrected: a
  single-mesh resolve with no a+b oracle stays FAIL-CLOSED (pin 5's own rule), or
  a tol-invariant resolve is labeled PROVISIONAL and never counted as a win.

---

## ATTACK 3 - V3-FATE COMPLETENESS

Walked src/overlap3.cpp + overlap3_sweep.cpp mechanism by mechanism against the
plan's section 1.4 fate table + stage-D per-guard list. Grep confirms plan
mentions (SplitTouchingSheets 0, sheet 0, face2Group/GroupedEdge/class-iii/loX/hiX
0, matching 0).

ACCOUNTED (explicit fate): SlabResolver self-location (table read at C, pin-dead
at D); provenance chains / edgeSubdiv (survive, compose - section 1.3);
ArrangementBudget (stays, domain shrinks); SubEpsFeature single-face coverage
guard (KEEP at D); chain-plane WIDE-RUN backstop (re-adjudicate at D); assembly
weld hash grid (no-op for junction images); Is2Manifold gate (KEPT);
EngineIdConflict/CoplanarOverlap/EdgeInPlane (retired); strict-FP slab gate
(survives as BuildSlabs kernel, tension named as R-ONCE-VS-STRICTFP).

UNACCOUNTED (no fate - agent will improvise): 2 primary + 2 loose.

1. SplitTouchingSheets (overlap3.cpp:985) - ZERO plan mentions. The design names
   it as the mechanism the wall breaks: MaintainedEmission3D L51-52, the twin is
   "a twin cap sheet plus a micro-edge, which the SHEET SPLITTER fails to pair."
   The plan's fate table jumps from "assembly weld no-op" straight to "Is2Manifold
   KEPT," skipping the radial-pairing + vertex-split step that sits BETWEEN them
   and is the wall's proximate failure. It ALSO handles genuine touching contacts
   (Touch_EdgeEdge_Cubes, Touch_VertexOnly_Cubes fences). Its fate under identity
   emission (dead for junctions? still runs for touching contacts?) is undefined.
   An agent deleting "structurally dead machinery" at stage D improvises here.
   STRONGEST gap.

2. The section-vert -> junction MATCHING/BIND mechanism - "matching" 0 mentions;
   "bind" is ASSERTED 6x, never mechanized. Section 1.4 / R-ONCE-VS-STRICTFP say
   section verts "must BIND into the global identity (a junction-imaging vert takes
   the junction's canonical)" - but HOW a section vert is identified AS junction
   V's image is exactly the matching the design flagged as a near-degeneracy
   decision that RELOCATES (variant ii, MaintainedEmission L216-222, L353-357).
   The whole architecture rests on this bind; the plan states it as given. This is
   the single most load-bearing UNDER-SPECIFIED mechanism - where the wall could
   re-enter unseen. R-ONCE-VS-STRICTFP names the risk but the plan supplies no
   mechanism, so stage B/C's builder improvises the matching.

3. Coplanar-group family (face2Group / GroupedEdge / CollectGroupedMemberEdges /
   SlabResolver class-iii planar / in-plane skeleton criticals, overlap3.cpp:345-
   395, 478-485, 528-536, 591-594) - "coplanar" appears only as the RETIRED
   FatalReason (CoplanarOverlap, L244), NOT the live grouping mechanism, which
   fences ~10 Coplanar_*_Oracle tests. Folded implicitly under "BuildSlabs
   survives," but whether a coplanar-group junction enters the identity table, and
   how class-iii planar resolution interacts with the canonical read, is
   unaddressed. LOOSE.

4. chain-plane loX/hiX exact cap-plane placement (StripChains, overlap3.cpp:662-
   665, 697-704, 943-944) - distinct from the wide-run backstop that stage D
   re-adjudicates. This "constructional closure" placement (money fixture
   Coplanar_PerpFacesSubEps depends on it) may die under name-joined strips, but
   stage D's guard list names only the wide-run rule. LOOSE.

Also loose: pair-canonical emission driver (EmitCaps `ci != li+1`) is "unified"
in stage C but its post-modification form is not pinned.

Count: 2 primary unaccounted (SplitTouchingSheets; the match/bind mechanism) + 2
loose (coplanar family; chain-plane loX/hiX placement).

---

## ATTACK 4 - STAGE-GATE EXECUTABILITY

Stage A DONE - EXECUTABLE. "three pins green (mutation-verified) + suite green +
emission bitwise-unchanged (diff of emitted geometry vs HEAD empty)." Objective
and checkable (needs a geometry-dump hook that does not yet exist - assumed
buildable, not pinned). Weakness is pin vacuity (attack 2), not checkability.

Stage B DONE - PARTLY JUDGMENT. "emitted-geometry diff CONFINED TO THE TRANSVERSE
BIND (no new resolves yet)." "Confined to the transverse bind" requires judgment
about which diff region counts - an implementing agent can self-certify. Make
executable: diff must be empty EXCEPT at the enumerated junction-imaging verts (a
named set), bitwise. "No new resolves" is checkable.

Stage C DONE - VACUOUS-GREEN-COMPATIBLE. "isolated twin oracle-true + corpus
re-adjudicated with zero oracle-wrong (SOME MAY newly resolve ... OTHERS STAY
recorded contracts)." The parenthetical is permissive, not a criterion: zero new
resolves satisfies it. Combined with pin 4b closing the synthetic twin
positions-only (attack 1), stage C is declarable DONE by machinery that closes a
synthetic twin AND fail-closes the entire real corpus - the v3 vacuous-green
pattern reproduced. Needs the progress ratchet.

Stage D DONE - MOSTLY EXECUTABLE, two soft spots. "guards deleted and pinned-dead"
+ "perf within bound" are checkable. "thinning landed+oracle-fenced OR deferred"
- the OR makes it trivially satisfiable by deferring. "ArrangementBudget domain
RE-ADJUDICATED" - "re-adjudicated" is a judgment word with no named check.

Cross-cutting: every "mutation-verified" DONE depends on the orchestrator
INDEPENDENTLY re-running the mutation (the plan mandates this in section 3, but it
lives outside the gate text). Executable only if enforced by a distinct verifier,
not the build agent.

---

## SUMMARY

VERDICT: NEED-CHANGE.
- Attack 1: softening is honest on necessity/sufficiency but sandbags in effect -
  no early sufficiency experiment (branch c violated) and pin 4b tests an
  identity-free-closable case. Add early sufficiency probe + load-bearing pin 4b
  mutation + progress ratchet.
- Attack 2: 4/5 pins stub-passable (1,2,4a,4b) + pin 3 partial + the
  pin-4c/4d/5 fail-closed loophole (zero-oracle-wrong satisfiable by fail-closing
  everything; CorpusSingleGate not yet oracle/tol-checked).
- Attack 3: 2 primary unaccounted mechanisms (SplitTouchingSheets; the
  section-vert->junction match/bind) + 2 loose (coplanar family; chain-plane
  loX/hiX).
- Attack 4: stage B "confined to transverse bind" and stage C "some may resolve"
  are judgment/permissive; stage D has an OR-escape and a "re-adjudicated"
  judgment word. Stage A is clean.

The design is sound and the plan is honest; the gates just need teeth: a
progress ratchet so zero-oracle-wrong is not fail-close-satisfiable, load-bearing
mutations so the pins prove the identity table does work, an early sufficiency
probe so the frontier is discovered cheaply, and explicit fates for
SplitTouchingSheets and the match/bind mechanism.
