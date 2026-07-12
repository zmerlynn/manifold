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
