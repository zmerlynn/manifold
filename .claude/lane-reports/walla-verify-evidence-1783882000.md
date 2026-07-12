# WALL-A tier-3 landing - adversarial verification lane

Artifact under audit: commit b68c6ea4 (doc-only) on explore/sweep-plane-3d-v3.
The WALL-A memo (docs/SweepEmit3D.md "WALL-A: the last correctness wall") +
notebook .claude/lane-reports/walla-arc-1783880487.md. Generation lane ran
env-gated probes, recorded evidence, reverted all code. This lane rebuilds the
probes in an isolated copy and tries to break the numbers.

Workspace: rsync copy at /tmp/walla-verify-evidence, fresh cmake -B vbuild
Release, -j4. Canonical repo READ-ONLY (this notebook is the only write).
HEAD b68c6ea4, git diff src/ EMPTY at HEAD (landing is doc-only, confirmed).

## Probe reconstruction (all env-gated, in the copy only; reverted at end)

The notebook records the probes by MECHANISM + LOCATION, not exact diffs. I
reconstructed each from the description:
- OV3_THIN=tol : BuildSlabs (overlap3_sweep.cpp) - thin ONLY arr.criticalXs
  (crossings) to representatives >= tol*eps apart, never vert x's. Notebook says
  "in FindSeams"; I applied it on the criticalXs copy in BuildSlabs, which is
  the consumer - the resulting crits set (verts.x U thinned criticalXs, sorted,
  deduped) is IDENTICAL either way, so the resolve outcome is faithful.
- OV3_FANDBG / OV3_FANDUMP : SplitTouchingSheets (overlap3.cpp) - read-only
  whole-surface pass classifying each edge (simple / boundary-hole / UNBAL /
  SLIVER / TIE / NONALT), min tie gap, cap-plane test, nearest-twin distance.
- OV3_WELDMUL=mult : BuildImpl getVertIdx (overlap3.cpp) - scale the weld radius
  AND the grid cell (else the 3x3x3 sweep misses enlarged-radius candidates -
  an unstated detail the notebook does not mention).

REPRODUCIBILITY AUDIT verdict: the probes ARE reconstructible from the notebook
descriptions and the headline numbers reproduce (below). Two unstated details a
reconstructor must recover: (i) TIE and NONALT are counted as INDEPENDENT arms
(an edge can be both) - only then do the histogram counts match; (ii) the weld
bump must scale the grid cell too. Neither is in the notebook. The probes are
described unambiguously enough to reconstruct, but NOT pinned to exact diffs.

## Claim 5 (cheap anchor): suite green at HEAD - SURVIVE
Overlap3+Boolean2+CrossSection: 166 passed, 1 skipped (Gate4c_HullMask). Exactly
as claimed. Landing is doc-only; git diff src/ empty.

## Claim 3 (fan anatomy) - SURVIVE (decisive)
fan_driver + OV3_FANDBG, whole-surface histogram:
- Havoc:  nTri=9738 simple=14566 boundHole=0 UNBAL=0 SLIVER=0 TIE=18 NONALT=13
- GT7863: nTri=69320 UNBAL=74 SLIVER=0 TIE=270 NONALT=198 boundHole=0
Both match the notebook table EXACTLY once TIE/NONALT are counted as independent
arms. boundaryHoles = 0 on both -> the dead-zone-c macro-hole claim (strict-FP
cured them) CONFIRMED.
TWIN GEOMETRY (OV3_FANDUMP, Havoc): NONALT edge(4,5) is cap-plane (dx=0); its
twin cap images sit at nearestTwin/eps a=1.6984, b=1.6815 - matching the
notebook's v4-v509=1.698, v5-v510=1.682 to the digit. Twins are ~1.7 eps apart,
just over the weld eps. Headline residual claim reproduced.

## Claim 4b (kAngleTie kill) - SURVIVE, with a scope nuance
minTieGap/kTie = 0 on both cases. RAW SCOPE (whole surface):
- Havoc:  tieGaps=23, exactZero=21 (91%), maxSubThreshGap=3.99e-13 = 0.0004 kTie
- GT7863: tieGaps=393, exactZero=370 (94%), maxSubThreshGap=5.18e-10 = 0.518 kTie
The DOMINANT mode is exact bitwise zero (91-94% of tie-gaps). The remainder are
nonzero but sub-threshold, up to 0.52 kTie on GT7863. The memo's "EXACTLY zero"
is the modal/minimum truth, not literally every tie; but the kill argument holds
- the exact-zero majority need kTie=0 (guard off) to un-fire, and coincident
sheets cannot be paired regardless of threshold.

## Claim 4a (weld-radius bump kill) - SURVIVE (outcome + mechanism)
Havoc with OV3_WELDMUL sweep, all fatal=2 (NonManifoldEmission), never resolves:
  mult=1.0 nTri=9738 TIE=18 NONALT=13
  mult=1.8 nTri=9520 TIE=8  NONALT=7   (218 triangles dropped)
  mult=2.5 nTri=9356 TIE=8  NONALT=6   (382 dropped)
  mult=4.0 nTri=9128 TIE=6  NONALT=5   (610 dropped)
Mechanism confirmed: merging the twins collapses cap triangles to degenerate,
they drop (nTri falls by 610 at 4x = the "vanishing pieces"), and the fan shrinks
but never closes (TIE/NONALT stay > 0). Bounded output weld does not resolve; it
destroys geometry. Matches the notebook's "does NOT resolve at ANY multiplier".

## Claim 2a (resolve preservation) - MIXED: core SURVIVE, "bitwise" NEED-CHANGE
verify_bitwise (%.17g + raw uint64 bit compare):
- Cray oracle: vol=1.5760799710533451e+116 bits=580fffffa1b5fe14, ORACLE-BITWISE
  at every tol AND bitwise-identical vs unthinned baseline. Matches the task's
  quoted number exactly. SURVIVE.
- Offset3: BITWISE-IDENTICAL across all tols. SURVIVE.
- Offset2: tol=off vol=209026.71902479883 -> thinned 209026.71902479368
  (bits ...12e0 -> ...122f), and nVert 16580 -> 16579. DIFFERS by ULPs.
- Offset4: tol=off 15240.582916271353 -> thinned ...271351 (1 ULP). DIFFERS.
FINDING: the memo/notebook claim "Offset2/3/4 BITWISE-IDENTICAL volume with and
without thinning" is FALSE for Offset2 and Offset4 - they shift by ULPs (Offset2
also drops a vertex). Only Offset3 (+ the Cray oracle) is bitwise. This is the
SAME fixed-precision-print artifact the STRICT-FP verification round already
corrected once ("only Cray survives an exact comparison") - repeated in the
wall-A memo. The SUBSTANTIVE claim (thinning is resolve-preserving: valid, genus
and component count preserved, volume oracle-correct to ULPs, tol-stable) HOLDS;
the word "bitwise" for the Offset trio is overstated. NEED-CHANGE (doc precision).
- Synthetic gates under OV3_THIN=1: 47 Overlap3 non-corpus gates PASS. SURVIVE.

## Claim 2b (self_A/B flip budget -> valid resolve) - SURVIVE (self_A)
self_intersectA: tol=off FATAL=3 (ArrangementBudget, "retained section content
exceeds budget"); tol=0.5/1/2/4 ALL RESOLVE valid=1 vol=0.14266215546791941
bits=3fc242c0e60e484d genus=1 comps=1, BITWISE tol-invariant across all four
tols (nVert drifts 2097027..2090345 but the volume bits are identical). Matches
notebook's 0.142662155468. Is2Manifold-equivalent: Status()==NoError (valid=1),
single component (comps=1), genus=1 - a clean closed manifold.
self_intersectB: tol=off FATAL=3 (budget); tol=0.5/1/2/4 ALL RESOLVE valid=1
vol=0.14266443142431906 bits=3fc242d3fda28b98 genus=1 comps=1, bitwise
tol-invariant. Matches notebook's 0.142664431424 (notebook listed 0.5/1/2; my
run adds tol=4, also identical).
A-vs-B CROSS-CHECK: |volA - volB| = 2.28e-6 (rel 1.6e-5) - matches the
notebook's "differ 2e-6". Two independent self-intersection meshes of the same
object resolving to matching volumes: corroboration reproduced.

## Claim 3 addendum (epistemics-lane tiebreaker): twin separation BY AXIS
Requested raw numbers, Havocglass8 headline twin pairs (independent axis-split
FANDUMP reproduction, eps units):
  v4 -> v509: dist/eps=1.6984  |dx|=0.0107  |dy|=0.4145  |dz|=1.6470
  v5 -> v510: dist/eps=1.6815  |dx|=0.0107  |dy|=0.4104  |dz|=1.6306
  other twin sites: v8->v13 |dx|=0.0000 transverse=8.07; v1452->v1451
  |dx|=0.0000 transverse=23.65; v428->v454 |dx|=0.2667 transverse=1.0488.
The EMITTED twin separation is TRANSVERSE-DOMINATED (dz ~1.65 eps, dx ~0.01
eps) - agreeing with the epistemics lane's reading of the arc's own dump, and
NOT with a reading of the memo's part-1 sentence "more than eps apart in x
while sub-eps in the transverse plane" as describing these twins. (That
sentence describes the Canonicalize-stage arr.verts pair, which I did not
instrument; but the memo's "dominated by the steep coordinate" phrasing next to
it invites the x-reading, which my measurement contradicts for the emitted
twins.) Raw numbers reported; interpretation is the epistemics lane's call.

## Claim 2c (hull budget -> fan) - behavior SURVIVES, "slower" DOES NOT
fan_driver on hull-body/hull-mask, two runs each (ulimit 6GB):
  tol=off : fatal=3 ArrangementBudget, 26.69s / 27.46s wall, 2.35GB peak
  tol=1   : fatal=2 NonManifoldEmission "unresolvable sheet contact",
            26.41s / 25.86s wall, 2.28GB peak
The BEHAVIOR change reproduces exactly: thinning lets hull section and land on
the cap-plane sheet fan (still fail-closed) instead of the retained-piece
budget. But the claimed SLOWDOWN does not reproduce: wall-time PARITY (thinned
marginally faster, and slightly less memory). Both paths are dominated by the
shared arrangement/seams phase; the thinned run's smaller crit set makes
BuildSlabs cheaper, offsetting the emission stages it unlocks. The memo's
not-landing reason #2 ("converts hull from a fast budget refusal into a slower
emission-stage refusal", notebook: "fast-budget -> 29s emission") is NOT
supported on this machine - the refusal is different, not slower. The OTHER
not-landing reasons (no oracle for the singles, owner-review gate) stand
untouched; but the perf argument should be dropped or re-measured by the owner.

## Claim 2 periphery: GT7081 thinned still budget - SURVIVE
GT7081 pair with OV3_THIN=1, budget intact, ulimit 4GB: fatal=3
ArrangementBudget, 24.35s, 1.35GB peak. Crossing thinning does not flip GT7081
(endpoint density remains), exactly as the notebook records.
## Anatomy reproduction (explosion driver, from the arc's own recorded driver)
All three light rows of the notebook table reproduce EXACTLY:
  Havoc:   critXs=1431 (8.13x nFaces) maxBundle=148  maxSpread/eps=1.540
  GT7863:  critXs=4476 (9.82x)        maxBundle=1369 maxSpread/eps=1.633
  Offset1: critXs=36098 (8.18x)       maxBundle=282  maxSpread/eps=1.635
Crossing fraction 87-94% of criticals - "crossings dominate" confirmed.

## Probe-neutrality check
Full suite re-run against the PATCHED library with all probe env vars unset:
166 pass + 1 skip, identical to HEAD. The reconstructed probes are env-gated,
zero default-path behavior change - the same property the arc claimed of its
instrumentation. Copy's src/ then reverted (git checkout, diff empty).

## VERDICT: NEED-CHANGE (two doc corrections; every mechanism claim SURVIVES)

Per claim:
1. REPRODUCIBILITY: SURVIVE with a caveat. Probes reconstructible from the
   notebook's mechanism descriptions; headline numbers reproduce. Not exact
   diffs; two unstated details (TIE/NONALT independent-arm counting, weld-bump
   grid-cell scaling) had to be recovered by matching outputs. Below the "exact
   diff" bar the lane brief asks about, but above the unreconstructible bar -
   I reconstructed everything without reading any generation-lane code.
2. THINNING: (a) resolve preservation - Cray ORACLE-BITWISE at every tol
   (1.5760799710533451e+116, bits 580fffffa1b5fe14); 47 synthetic gates green
   under OV3_THIN=1; BUT "Offset2/3/4 bitwise-identical" is FALSE for Offset2
   and Offset4 (ULP shifts, Offset2 drops a vert) - only Offset3 is bitwise.
   NEED-CHANGE: the memo repeats the exact fixed-precision-print artifact the
   strict-FP verification round already corrected; change "BITWISE-IDENTICAL
   volume" to "identical to ULPs" for the Offset trio.
   (b) self_A/B flips - SURVIVE, bitwise tol-invariant across ALL tols
   (A=0.14266215546791941, B=0.14266443142431906, |A-B|=2.28e-6), valid
   single-component genus-1 manifolds.
   (c) hull - behavior flip SURVIVES (budget -> sheet-fan refusal), but the
   claimed SLOWDOWN does not reproduce: wall parity ~27s vs ~26s over two runs
   each. NEED-CHANGE: drop or re-measure the "slower emission-stage refusal"
   justification; the honest statement is "a different refusal at the same
   cost (on this machine), still fail-closed".
3. FAN ANATOMY: SURVIVE, decisive. Histograms match exactly; boundaryHoles=0
   (dead-zone-c gone under strict-FP); twin pair at 1.6984/1.6815 eps matches
   the notebook to 4 digits. AXIS SPLIT (epistemics tiebreaker): the emitted
   twin separation is TRANSVERSE-dominated (|dz|~1.65 eps, |dx|~0.01 eps) -
   the memo's axis phrasing needs the epistemics lane's correction.
4. KILLS: SURVIVE. (a) weld bump: no resolve at 1.8/2.5/4x; mechanism confirmed
   (610 triangles dropped at 4x, fan shrinks but never closes = degenerate
   collapse). (b) kAngleTie: 91-94% of tie gaps are EXACT bitwise zero
   (21/23 Havoc, 370/393 GT7863; residue max 0.52 kTie, still far below any
   sane threshold) - no threshold value fixes coincident sheets. SCOPE: kill
   evidence verified on Havoc+GT7863 (2 of 4 sheet-fan carriers) - a SAMPLE,
   not the population; Offset1/openscad not re-probed here (the arc's own
   histogram covered all four; my reproduction covers two).
5. SUITE: SURVIVE. 166 pass + 1 skip (Gate4c_HullMask) at HEAD in a fresh
   build; landing is doc-only (git diff src/ empty at b68c6ea4).

The tier-3 memo's SCIENCE is sound - anatomy, kills, thinning mechanism, and
the irreducible-coupling conclusion all reproduce. The NEED-CHANGE is confined
to two overstated evidence sentences (Offset-trio "bitwise", hull "slower") and
the axis phrasing; none of the three affects the verdict structure (killed
probes stay killed, validated-safe stays validated-safe, not-landing still has
standing reasons without the hull-perf one).

## Cleanup
Probe patches reverted in the copy (git checkout src/); canonical repo
untouched except this notebook. Probe drivers remain in the session scratchpad;
verify_bitwise.cpp (this lane's bitwise driver) copied alongside them.

