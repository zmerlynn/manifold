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
