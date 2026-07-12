# STRICT-FP SLAB BUILDING - lab notebook

Branch explore/sweep-plane-3d-v3 @ 89f52a6f. Solo, sole editor. Build build/ (Release).
Suite green at start (baseline claim): Overlap3 = 167 pass + 1 skip (Gate4c) per mission.

## Mission

Adopt STRICT-FP slab building: lower the BuildSlabs slab-width gate
(overlap3_sweep.cpp:146, `if (slab.xHi - slab.xLo <= eps)`) to strict FP ordering
(build every slab whose midpoint is strictly between its bounds in double
arithmetic), and DISSOLVE the skipped-run machinery the eps gate necessitated.

Three prior arcs unblocked this:
- EngineIdConflict demotion (wall-B): the near-coplanar attribution fatal that
  the predecessor probe hit is now a benign counter. srcId unconsumed in 3D.
- ArrangementBudget guard: max(4M, 64*nFaces) cumulative retained pieces; bounds
  the heavy strict-FP runs (was NOT present at the probe's b6443958).
- sub-eps probe: proved the mechanics (bitwise-identical resolves on clean
  geometry: Cray/Offset2/3/4).

The prize: content in formerly-skipped sub-eps runs gets PROCESSED (sectioned at
each critical) instead of refused; the chain-plane rule / wide-run guard /
coverage guard / non-canonical-skip machinery deletes IF structurally dead.

## Key facts established from grounding read (done)

- Gate is at overlap3_sweep.cpp:146. xMid computed line 144 = (xLo+xHi)*0.5.
  Strict-FP gate: build iff `xLo < xMid && xMid < xHi` (probe recipe).
- Machinery that exists FOR skipped runs (to re-adjudicate in Phase B):
  1. EmitCaps non-canonical skip `if (ci != li + 1) continue;` (overlap3.cpp:927)
  2. chain-plane loX/hiX (StripChains, EmitCaps bind :932-933, ZipperEmit :698)
  3. wide-run guard (ComputeCap wideRun param, SubEpsFeature fatal :860; wideRun
     computed EmitCaps :942 = crits[ri]-crits[li+1] > eps)
  4. coverage guard (SubEpsFeature single-face loop, overlap3_sweep.cpp:246-274)
- Under strict-FP an unbuilt slab is ~1-ulp wide (xMid rounds to a bound). A
  "skipped run" is a maximal sequence of criticals each 1 ulp apart -> run width
  is (#gaps)*1ulp, which can only exceed eps (~1000 ulps) with ~1000 consecutive
  1-ulp criticals. So wideRun `crits[ri]-crits[li+1] > eps` is PREDICTED
  unreachable; ci==li+1 is NOT always-true (1-ulp gaps persist) so the
  non-canonical skip is NOT dead; loX/hiX displacement shrinks to ~1 ulp (weld-
  valid) but is not structurally zero. TO VERIFY EMPIRICALLY, not assume.
- Money fixture: Coplanar_PerpFacesSubEpsApart_Recorded (test:1502). Run of
  slabs[1..2] ~1.00004 eps, flanks a-square vs b-square (MACRO diff). Under
  eps gate the content collapses into the run (resolve-wrong without the guard,
  fail-closed with it). Under strict-FP those two ~0.5eps slabs BUILD -> content
  sectioned -> does it resolve ORACLE-CORRECT? The arc's money question.
- Cray/Offset2/3/4 are NOT Overlap3 tests (they live in BooleanComplex, the a+b
  path). Measured via the driver. Bitwise targets: Cray 1.5760799710533451e+116,
  Offset2 209026.71902485067, Offset3 10803.472688795122, Offset4 15240.582916271358.
- Budget guard is IN the tree (kPieceBudget max(4M,64*nFaces)). So heavy strict-FP
  runs (self_A/B growing 5-6x built slabs) may trip ArrangementBudget FAST instead
  of the probe's 434s emission blowup. C3 interplay = measure + adjudicate.

## Plan (phases per mission brief)

A. RE-BASELINE + measure. Env-gated OV3_STRICTFP probe in BuildSlabs (marked TEMP
   PROBE), commit checkpoint. Adapt predecessor driver (/tmp/subeps_driver.cpp,
   fix deleted EngineIdConflict, env->OV3_STRICTFP). Per-case OFF vs ON table:
   outcome class, oracle vol (a+b), wall time, retained pieces vs budget. Full
   Overlap3 suite ON census. Heavy cases (self_A/B, Offset1, openscad, Havoc,
   hull) first run under (ulimit -v 6000000; timeout 900); GT7081 under
   (ulimit -v 4000000; timeout 900). Output to scratchpad files, bounded reads.
B. ADJUDICATE adoption shape + STRUCTURAL DEATH. Verify each of the 4 mechanisms:
   dead (delete + re-adjudicate pins + spec supersede) or live (keep, do not force).
C. IMPLEMENT. C1 red-first money fixture (OracleCompare before flipping contract);
   guard re-labels honest; mutation-verify surviving pins. C2 deletion commits
   SEPARATE from gate change. C3 budget interplay adjudication. C4 spec amend
   (chain-plane / [R3] / WALL-B / corpus + new STRICT-FP section).
D. PERF BOUND. No single test > ~120s; filtered Overlap3+Boolean2 within ~2x
   baseline. If exceeded: BuildSlabs face-range sweep (sort by x-interval). No
   eps-gate-as-fast-path fallback. Fidelity gate ABSOLUTE: zero oracle-wrong.
E. WRAP. Verdict, class map, suite counts, machinery delta, perf table, risks.

Probe commits reverted before landing. Landing history clean.

## Log

### Step 0: probe wired + driver adapted (DONE)

Env-gated OV3_STRICTFP in BuildSlabs (static getenv, strict-FP gate
`xLo < xMid < xHi`), TEMP PROBE marked, committed as probe checkpoint. Driver =
predecessor /tmp/subeps_driver.cpp adapted (env->OV3_STRICTFP, EngineIdConflict
deleted, ArrangementBudget added, + PerpFaces/InPlaneSkel synthetic composed
pairs + interior-run RUNSTATS). Builds clean, links libmanifold.so.

### Step 1: per-case OFF vs ON evidence (DONE)

built = built slabs; pieces = retained section pieces; budget = 4M floor.

| case        | OFF class          | OFF vol/detail   | ON class            | ON vol/detail  | built OFF->ON | pieces ON | t OFF->ON |
|-------------|--------------------|------------------|---------------------|----------------|---------------|-----------|-----------|
| Cray        | RESOLVE            | 1.576079971e116  | RESOLVE             | SAME (bitwise) | 1->1          | 8         | 0.0/0.0   |
| Offset2     | RESOLVE            | 209026.719       | RESOLVE             | SAME (bitwise) | 583->1283     | 40074     | 0.26/0.31 |
| Offset3     | RESOLVE            | 10803.47269      | RESOLVE             | SAME (bitwise) | 15->32        | 518       | ~0        |
| Offset4     | RESOLVE            | 15240.58292      | RESOLVE             | SAME (bitwise) | 21->53        | 1155      | ~0        |
| PerpFaces   | **SubEpsFeature**  | fail-closed      | **RESOLVE**         | vol=2=oracle   | 2->3          | 30        | ~0        |
| InPlaneSkel | RESOLVE            | 1.364501847      | RESOLVE             | SAME correct   | 12->22        | 270       | ~0        |
| Gate4a_W8   | RESOLVE (manifold) | -                | RESOLVE (manifold)  | manifold       | -             | -         | (126ms)   |
| Havoc       | NonManifoldEmission| sheet contact    | NonManifoldEmission | SAME           | 153->409      | 11319     | 0.06/0.07 |
| GT7863      | **SubEpsFeature**  | wide-run guard   | **NonManifoldEmission** | sheet fan  | 245->741      | 87451     | 0.19/0.99 |
| Offset1     | NonManifoldEmission| sheet fan        | NonManifoldEmission | SAME           | 1598->6046    | 497041    | 1.45/4.04 |
| openscad    | **SubEpsFeature**  | wide-run guard   | **NonManifoldEmission** | sheet fan  | 3492->7766    | 1165051   | 4.45/17.7 |
| hull        | **SubEpsFeature**  | wide-run guard   | **ArrangementBudget** | >4M dense    | 5014->budget  | >4M       | 9.6/27.6  |
| self_A      | NonManifoldEmission| sheet fan        | **ArrangementBudget** | >4M dense    | 8932->budget  | >4M       | 36.3/13.0 |
| self_B      | NonManifoldEmission| sheet fan        | **ArrangementBudget** | >4M dense    | 8926->budget  | >4M       | 29.3/12.0 |
| GT7081      | ArrangementBudget  | >4M dense        | ArrangementBudget   | SAME          | budget        | >4M       | 22.7/29.6 |

FIDELITY: zero oracle-wrong resolves in either mode. Every resolve oracle-correct
(clean 4 bitwise; PerpFaces=2; InPlaneSkel/Gate4a correct). The eps gate was NOT
load-bearing for correctness anywhere. Money fixture confirmed: PerpFaces
SubEpsFeature -> oracle-correct RESOLVE = the arc's prize (content in the former
sub-eps run now sectioned faithfully instead of collapsed-and-refused).

### Step 2: full Overlap3 ON census (DONE, minus self_A/B/GT7081 run separately)

51 pass, 2 fail:
- Gate4c_HullMask: hull -> ArrangementBudget (fatal=3), not in its accept list
  {SubEpsFeature-chainplane, NonManifoldEmission} -> honest RE-LABEL needed.
- Pin_ChainPlaneRule_WideRunResolves: white-box ASSERT_GT(widestInteriorRun,eps)
  FAILS (1.11e-16 vs 5.5e-12) - the fixture's wide run VANISHES under strict-FP;
  the pin retires WITH the machinery it exercises (predicted by the probe).
Everything else GREEN incl. Coplanar_PerpFacesSubEpsApart (now resolves),
Gate4a_Wedges8, Pin_InPlaneSkeletonCriticals, GT7863/Offset1/openscad/Havoc
(all accepted named guards). Boolean2 untouched (2D seeds srcId 0).

### Step 3: interior unbuilt-run structure under strict-FP (DONE) - the Phase B pivot

RUNSTATS (interior runs flanked by built both sides):
| case | nRuns | nNonCanonCrit | maxRunW/eps |
| Cray | 0 | 0 | 0 | Offset2 | 510 | 591 | 0.0002 | Offset3 | 14 | 24 | 0.0002 |
| Offset4 | 21 | 36 | 0.0003 | PerpFaces | 1 | 0 | 0.0000 | InPlaneSkel | 8 | 3 | 0.0001 |
| Havoc | 59 | 62 | 0.0005 | GT7863 | 60 | 89 | 0.0007 | Offset1 | 1480 | 2354 | 0.0002 |
| openscad | 581 | 447 | 0.0002 |

FINDINGS (structural death adjudication input):
1. Skipped runs are NOT impossible - they persist at ULP scale (a few adjacent-
   double criticals; maxRunW/eps <= 0.0007 across every case). FP-floor retention,
   not gate removal (1-ulp slabs can't be built: xMid would coincide with a
   critical plane -> degenerate section).
2. ci==li+1 is NOT always true - multi-slab runs exist (nNonCanonCrit>0), so the
   NON-CANONICAL SKIP is exercised (Offset2: 591 non-canonical criticals). Deleting
   it re-introduces doubled caps -> ALIVE, keep.
3. The WIDE-RUN GUARD (fires iff run > eps) is UNREACHABLE in practice: widest
   interior run is 0.0007 eps. Not provably impossible (needs ~1000 consecutive
   adjacent-double criticals) but empirically dead on all real+synthetic geometry.
4. The CHAIN-PLANE loX/hiX placement now bridges only a few-ulp displacement
   (<< eps weld tolerance) - no longer LOAD-BEARING (the weld would close it), but
   still active. Deletion candidate PENDING the Step-4 experiment.

Probe checkpoint commit; heavy cases wrapped (ulimit -v; timeout). Budget guard
(in tree) bounded every heavy run - no OOM, no 434s blowup.

### Step 4: PHASE B structural-death experiment (DONE) - the probe's dissolution hope is REFUTED

Made the gate unconditional (strict-FP), rebuilt. Full Overlap3 suite = 54 pass,
2 fail (Gate4c budget re-label; Pin_ChainPlaneRule white-box). Gate change alone
is sound; every MUST-RESOLVE test still resolves (only Recorded tests could mask
a degrade, checked below).

Then DELETED the two run-spanning mechanisms (chain-plane loX/hiX placement +
wide-run guard) to test the probe's "they degenerate to no-ops" prediction:
- Suite stayed 54 pass / 2 fail (no NEW test failures) - BUT:
- The MONEY FIXTURE REGRESSED SILENTLY: PerpFaces RESOLVE(vol=2) -> FATAL
  NonManifoldEmission. The Recorded contract accepts the fatal, so the suite
  masked it; the driver caught it.

DIAGNOSIS. Removing a fatal (the wide-run guard) CANNOT create a fatal, so the
culprit is provably the CHAIN-PLANE PLACEMENT (ZipperEmit at slab bounds vs at
the cap plane). Under strict-FP, PerpFaces has a single 1-ULP interior run; the
chain-plane rule makes the post-run strip's lo corner BITWISE the cap corner,
and that exact coincidence is LOAD-BEARING for closure - the 1-ulp displacement
at slab bounds does NOT reliably weld (the sheet splitter leaves an unpaired
fan -> "unresolvable sheet contact"). The probe ASSUMED loX/hiX == slab bounds
under strict-FP (a no-op); RUNSTATS shows they differ by the ulp-scale run width,
and PerpFaces proves that difference matters. "Verify rather than assume" earned
its keep here.

VERDICT on structural death (evidence-driven, corrects the probe):
- CHAIN-PLANE RULE (loX/hiX placement): NOT dead - LOAD-BEARING at ulp scale
  (PerpFaces). KEEP. Runs shrink from eps-scale to ulp-scale but do not vanish,
  and the rule's exact-coincidence closure is still required.
- WIDE-RUN GUARD: empirically UNREACHABLE (widest run 0.0007 eps; a >eps run needs
  ~1000 consecutive adjacent-double criticals). Deleting it is a no-op on all real
  input. BUT with the chain-plane rule KEPT (it spans runs by linear interp), the
  guard is the ONLY thing converting a pathological wide-run macro change from a
  silent oracle-WRONG resolve into a fail-closed. Fidelity gate is ABSOLUTE ->
  KEEP as the rule's coupled safety net (honest comment: unreachable-in-practice).
- NON-CANONICAL SKIP (ci != li+1): ALIVE (multi-slab ulp runs; deletion -> doubled
  caps). KEEP.
- COVERAGE GUARD (single-face SubEpsFeature): dormant (no corpus case trips it ON)
  but not structurally dead (a near-x-perpendicular macro face landing wholly in a
  ulp run has no built-slab coverage). KEEP.

So NO run machinery is structurally dead. The arc's CORE prize - content in
formerly-skipped runs gets PROCESSED (PerpFaces resolves oracle-true; GT7863/
openscad process into emission) - is delivered by the GATE CHANGE ALONE. The
probe's hoped-for machinery dissolution does not happen; the residual ulp-scale
runs still need the rule + guard. Restored overlap3.cpp to all-machinery-present.

### Step 5: PHASE B adoption shape (DECIDED)

- ADOPTION: FP-floor RETENTION (keep a minimal `xLo < xMid < xHi` gate) - forced,
  since 1-ulp slabs cannot be built (xMid would coincide with a critical plane).
  Gate removal is not a thing.
- BUDGET INTERPLAY (C3): hull/self_A/self_B flip NonManifoldEmission/SubEpsFeature
  -> ArrangementBudget (they genuinely retain >4M pieces once every ulp slab in
  the dense near-coplanar bands builds). DO NOT recalibrate up: relaxing the
  budget would let them run to the 434s emission blowup and STILL fail
  (NonManifold) - strictly worse. Accept the re-label as an equally-honest,
  FASTER refusal (self_A 36->13s, self_B 29->12s). 4M floor stays.
- CLASS MAP (old -> new):
  Cray/Offset2/3/4: RESOLVE -> RESOLVE (bitwise).
  PerpFaces: SubEpsFeature(fail) -> RESOLVE oracle-true. [STRENGTHEN - money]
  InPlaneSkel/Gate4a: RESOLVE -> RESOLVE (probe's EngineIdConflict breakage gone
    via wall-B demotion).
  GT7863/openscad: SubEpsFeature(wide-run guard) -> NonManifoldEmission (wall-A
    fan; formerly-refused content now processed).
  hull: SubEpsFeature(wide-run guard) -> ArrangementBudget.
  self_A/self_B: NonManifoldEmission -> ArrangementBudget.
  Offset1/Havoc: NonManifoldEmission -> NonManifoldEmission (same).
  GT7081: ArrangementBudget -> ArrangementBudget (same).
- TEST/PIN plan: (a) promote Coplanar_PerpFacesSubEpsApart to MUST-RESOLVE (red-
  first: it fatals under the eps gate, resolves under strict-FP; doubles as the
  chain-plane rule's mutation pin). (b) RETIRE Pin_ChainPlaneRule (RingedBox's
  eps-scale ring gaps all BUILD under strict-FP -> no interior run to exercise;
  the rule is now pinned by PerpFaces at ulp scale). (c) Gate4c: accept
  ArrangementBudget. (d) corpus comments: self_A/B residual now ArrangementBudget.

### Step 6: PHASE C implement (DONE)

Gate change unconditional (overlap3_sweep.cpp strict-FP gate, cstdlib probe
include removed). overlap3.cpp: ALL machinery KEPT (verified: git diff is
comments-only); comments on the chain-plane rule / wide-run guard / StripChains /
EmitStrips updated to the ulp-scale + fidelity-backstop framing. Tests:
- Coplanar_PerpFacesSubEpsApart_Resolves: promoted to MUST-RESOLVE. RED-FIRST
  verified (eps gate: fatal=1 wide-run guard). MUTATION verified (revert ZipperEmit
  to slab bounds: fatal=2 unresolvable sheet contact) - confirms the chain-plane
  rule is load-bearing at ulp scale.
- Pin_ChainPlaneRule_WideRunResolves + RingedBox helper: DELETED.
- Gate4c: skip clause -> ArrangementBudget (pinned to detail).
- Corpus comments (single-gate header, Offset1/openscad/self_A/self_B, pair gate).
Spec (C4): new STRICT-FP SLAB BUILDING section (decision/evidence/what-dissolves-
nothing/residuals/tests) + superseding notes on CHAIN-PLANE RULE and Corpus
fail-closed attribution.

Full Overlap3 suite: 55 tests = 54 PASS + 1 SKIP (Gate4c, arrangement budget).
Every heavy case bounded by the 4M budget (no OOM). C3 budget: NOT recalibrated
(relaxing would run self_A/B to the emission blowup and still fail - worse);
re-label accepted as equally-honest + faster.

### Step 7: PHASE D perf (DONE) - within bound, no optimization needed

Overlap3 + Boolean2 filtered suite: 90 PASS + 1 SKIP, ZERO failures, 103.5s.
Slowest single test GT7081 29.7s (< 120s cap). Boolean2 fully green + unchanged
(2D seeds srcId 0). Per-case vs baseline:
| case | base | strictFP | ratio |
| self_A | 33.6 | 12.5 | 0.37x | self_B | 29.3 | 11.8 | 0.40x |
| GT7081 | 22.7 | 29.7 | 1.31x | hull/Gate4c | ~9 | 27.4 | 3.0x |
| openscad | ~4.5 | 16.3 | 3.6x | Offset1 | 1.5 | 3.9 | 2.6x |
Corpus subtotal ~101s -> ~102s (~1.02x): the self_A/B early-budget-trip speedups
offset the hull/openscad slab-count slowdowns. Suite total within 2x, no single
test over 120s -> the pre-profiled BuildSlabs face-range sweep is NOT needed
(recorded as an available future optimization if a denser case pushes hull/openscad
past the cap). Fidelity gate: zero oracle-wrong resolves anywhere, both modes.

## VERDICT (Phase E)

STRICT-FP SLAB BUILDING ADOPTED. The BuildSlabs gate lowers from eps-width to the
FP floor (xLo < xMid < xHi). The PRIZE is delivered: content in formerly-skipped
sub-eps runs is now SECTIONED, not refused - the money fixture PerpFaces flips
SubEpsFeature-refusal -> oracle-true RESOLVE (red-first + mutation verified), and
GT7863/openscad process their in-run content into emission. Zero oracle-wrong
resolves anywhere; clean pairs bitwise. The eps gate was NOT load-bearing for
correctness - only for containing near-coplanar density (now handled by the
wall-B demotion + budget).

MACHINERY DELTA: ZERO deletions, ZERO additions. This CORRECTS the sub-eps probe's
central prediction. Unbuilt runs do not vanish - they persist at ulp scale, and
the empirical test refutes the "chain-plane rule degenerates to a no-op" claim:
deleting the placement reds PerpFaces (1-ulp displacement does not reliably weld).
Chain-plane rule KEPT (load-bearing at ulp scale); wide-run guard KEPT (unreachable
but the rule's fidelity backstop); non-canonical skip + coverage guard KEPT (alive).

CLASS MAP (eps -> strict-FP): Cray/Offset2/3/4 resolve->resolve (bitwise);
PerpFaces SubEpsFeature->resolve; InPlaneSkel/Gate4a resolve->resolve;
GT7863/openscad SubEpsFeature->NonManifoldEmission; hull SubEpsFeature->
ArrangementBudget; self_A/self_B NonManifoldEmission->ArrangementBudget;
Offset1/Havoc NonManifoldEmission (same); GT7081 ArrangementBudget (same).

SUITE: Overlap3 55 = 54 pass + 1 skip (Gate4c on budget). Overlap3+Boolean2 = 90
pass + 1 skip, 0 fail, 103.5s (~1.02x baseline, max single test 29.7s). Boolean2
unchanged. overlap3.h enum untouched -> no other suite affected.

COMMIT RANGE: 4115351e (TEMP PROBE, marked; env-gating reverted in the landing)
.. HEAD. Landing = 94eb88df (gate+tests+comments), 0075c473 (spec), + this
notebook. Branch explore/sweep-plane-3d-v3, local only.

OPEN RISKS (plain):
- The wide-run guard is now UNTESTED (no fixture reaches it under strict-FP) and
  unreachable on realistic input. It is retained on a fidelity argument, not a
  test. A future reviewer could reasonably challenge keeping unreachable code;
  the counter is that deleting it, with the interpolating chain-plane rule kept,
  opens a pathological silent-oracle-wrong path.
- hull/self_A/self_B ArrangementBudget re-label: honest but it means these no
  longer reach their emission residual (masked behind the budget). If the wall-A
  arrangement arc later reduces density, they would surface NonManifoldEmission.
- The chain-plane placement's ulp-scale load-bearingness suggests a latent weld/
  closure fragility (a 1-ulp gap that should weld but does not pair into a sheet).
  Not investigated here - the rule papers over it correctly. A future arc could
  root-cause that and THEN the rule might genuinely delete.
- Perf headroom is thin on hull/openscad (3x, ~16-27s). A denser future corpus
  case could push past 120s; the BuildSlabs face-range sweep is the ready fix.
