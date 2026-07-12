# STRICT-FP slab-gate adoption - ADVERSARIAL FIDELITY VERIFY lane

Auditor lane. Canonical repo READ-ONLY; work in /tmp/strictfp-verify-fidelity
(rsync + FRESH `cmake -S . -B vbuild -DCMAKE_BUILD_TYPE=Release`, touch-test
verified the build tracks the copied tree, not a stale cache). Artifact under
audit: commits 89f52a6f..ef10175e, landing 94eb88df (gate+tests) + 0075c473
(spec). VERDICT = NEED-CHANGE (one empirical overclaim; fidelity gate itself
holds).

## Method
- BUILD1 = HEAD (strict-FP gate, unconditional). Full Overlap3 + Boolean2 +
  CrossSection suites, per-test timing. Fresh driver (vdriver.cpp) linking the
  copy lib: per-case class + full-precision (%.17g + %a hex) resolved vol +
  a+b oracle vol + eps*SA bound + wall time.
- BUILD2 = eps-gate restored in my copy (`slab.xHi - slab.xLo <= eps`). Driver
  clean-4 + PerpFaces red-first.
- BUILD3 = strict-FP gate + chain-plane mutation (EmitStrips ZipperEmit loX/hiX
  -> slabs[si].xLo/xHi). PerpFaces mutation-pin.
- All heavy runs (hull/self_A/self_B/GT7081) under (ulimit -v 6000000; timeout
  900). Mutations reverted; copy source diff-clean vs HEAD at end.

## Claim 1 - suite green + counts: SURVIVE
- Overlap3: 54 PASS + 1 SKIP (Gate4c_HullMask). Matches the lane report.
- Boolean2 + CrossSection: 111 PASS, 0 fail (36 + 75).
- Count delta 56->55 tests accounted test-by-test (git diff 89f52a6f..HEAD):
  ONE rename (PerpFacesSubEpsApart_Recorded -> _Resolves) + ONE deletion
  (Pin_ChainPlaneRule_WideRunResolves), and the RingedBox helper it solely used.
  No other TEST added/removed. The -1 net = the sanctioned pin-retirement.
  PerpFaces contract Recorded(accept-fatal-or-resolve) -> MUST-RESOLVE is a
  strict narrowing. No fence violation.

## Claim 2 - money claim + fidelity: MOSTLY SURVIVE, ONE OVERCLAIM (NEED-CHANGE)
OracleCompare (test:292) checks THREE things, not volume only: volume within
eps*max(SA) bound, Genus() equality, AND a 17^3 winding-number grid (<5
disagreements). Strong oracle. Coplanar_PerpFacesSubEpsApart_Resolves passes at
HEAD (4 ms).

RED-FIRST (eps-gate build): PerpFaces FATALs with fatal=1 SubEpsFeature "macro
cap content over a skipped run wider than eps" (the wide-run guard). Confirms
the money claim: eps gate refuses, strict-FP resolves oracle-true.

"zero oracle-wrong resolves anywhere": HOLDS. Every resolving fixture is
oracle-correct in BOTH modes; every corpus resolve/refusal is a named guard;
the whole Overlap3 suite (which runs OracleCompare on all MUST-RESOLVE tests) is
green. Clean-pair oracle bounds at HEAD:
  Cray    voldiff 0.000e+00  (bitwise == oracle)
  Offset2 voldiff 4.28e-08   bound 1.08e-02
  Offset3 voldiff 5.46e-12   bound 2.92e-07
  Offset4 voldiff 3.64e-12   bound 3.57e-07
All oracle-correct with 5+ orders of margin.

*** OVERCLAIM (NEED-CHANGE): "Cray/Offset2/3/4 BITWISE" is FALSE for 3 of 4. ***
The mission's target constants are the EPS-GATE outputs (verified: eps-gate
Offset2=209026.71902485067, Offset3=10803.472688795122, Offset4=15240.582916271358
- exact matches; both eps-gate runs deterministic). The strict-FP (HEAD) outputs
DIFFER by a few ulps and are deterministic too:
  case     strict-FP hex            eps-gate hex             bitwise?
  Cray     0x1.fffffa1b5fe14p+385   0x1.fffffa1b5fe14p+385   YES
  Offset2  0x1.98415c09012e p+17    0x1.98415c09019d5p+17    NO (~ulp e-10 rel)
  Offset3  0x1.519bc8111021b p+13   0x1.519bc8111021e p+13   NO (3 ulp)
  Offset4  0x1.dc44a9d0018e2 p+13   0x1.dc44a9d0018e5 p+13   NO (3 ulp)
Building more slabs perturbs the arrangement/emission rounding; ONLY Cray is
bitwise-identical between modes. Both modes stay oracle-correct, so this is NOT a
fidelity break - it is a wording error. Root cause: the arc's driver printed
%.10g, at which all three round to a shared display value; "bitwise" was inferred
from a 10-sig-fig match. The claim appears in THREE places, all should read
"oracle-correct / identical to display precision; only Cray bitwise":
  - spec docs/SweepEmit3D.md:1166 "resolve BITWISE-identically"
  - commit 94eb88df msg "Every resolve stays oracle-correct (Cray/Offset2/3/4 bitwise)"
  - notebook strictfp-arc Step 1 "ON SAME (bitwise)" for Offset2/3/4

## Claim 3 - pin retirement re-anchors on PerpFaces: SURVIVE
Chain-plane mutation (ZipperEmit at slab bounds) reds PerpFaces with fatal=2
NonManifoldEmission "unresolvable sheet contact" - exactly the deleted
Pin_ChainPlaneRule's mutation duty, now carried by PerpFaces. The chain-plane
rule is NOT unpinned. (fatal=2 not 3: the EngineIdConflict enum retirement
shifted NonManifoldEmission from index 3 to 2; same guard.) Gate4a_Wedges8 and
Pin_InPlaneSkeletonCriticals both GREEN at HEAD (the demotion cured the probe-era
breakage). RingedBox is unreferenced elsewhere - its deletion is sanctioned.

## Claim 4 - class map: SURVIVE (all observed classes match)
Driver at HEAD, every class as documented:
  GT7863 NonManifoldEmission | openscad NonManifoldEmission
  hull ArrangementBudget | self_A ArrangementBudget | self_B ArrangementBudget
  GT7081 ArrangementBudget (unchanged) | Havoc NonManifoldEmission | Offset1 NonManifoldEmission
Acceptance-list notes (contract judgment owned by other lanes, flagged not
adjudicated): CorpusSingleGate accepts ANY fatal or valid manifold - the class
comments are documentation, NOT enforced (a future class flip would not red the
test); I confirmed the actual classes match the comments empirically.
CorpusPairGate retains SubEpsFeature in its OR-list though no pair fixture fires
it under strict-FP (reassigned in-comment to the still-live coverage guard).
Gate4c retains a NonManifoldEmission skip clause that does not currently fire
(framed as a forward boundary); the OLD replaced class (SubEpsFeature chain-plane)
was correctly removed. None is a stale widening of the replaced class.

## Claim 5 - perf decomposition: SURVIVE, thin GT7081 headroom flag
No single test > 120s. Driver wall times at HEAD (shared box, runs high):
  Cray 0.00 Offset2 0.31 Offset3/4 0.01 GT7863 0.86 Havoc 0.07 Offset1 3.97
  openscad 16.49 hull 30.89 self_A 12.91 self_B 12.38 GT7081 35.10
The ~1.02x corpus parity is a NET, not resolve-path parity - split:
  - genuinely-unchanged-ANSWER resolves (2-2.5x more slabs, still <0.5s):
    Cray/Offset2/3/4 (oracle-correct; see claim 2 they are NOT bitwise).
  - REFUSE-EARLIER speedups (budget trip during slab build vs former full-section
    NonManifoldEmission run): self_A 33.6->12.9, self_B 29.3->12.4 (~-38s).
  - process-MORE slowdowns (former wide-run-guard refusals now run deeper):
    openscad 4.5->16.5, hull 9->30.9, Offset1 1.5->4.0, GT7863 (~+38s).
  - unchanged: GT7081 (ArrangementBudget both modes).
The self_A/B speedups OFFSET the hull/openscad/Offset1 slowdowns -> the ~1.02x is
arithmetic coincidence, correctly framed by the notebook ("early-budget-trip
speedups offset the slab-count slowdowns") and NOT overstated as resolve parity
by the spec (its PERF line mentions only the slowdown). FLAG: GT7081 in the full
suite measured 105.4s on this contended box (vs notebook 29.7s / driver 35s) -
under the 120s cap but thin, environment-dependent headroom; on a quiet box it is
~30-35s. Not a claim contradiction, a robustness note.

## Claim 6 - docs/style: SURVIVE
STRICT-FP spec section (1144-1220) ASCII-clean, magnitudes-only (2-2.5x, ~1000
criticals, <0.001 eps, "hundreds of non-canonical criticals" - no stale exact
counts; my Offset2 built 583->1283 = 2.2x matches "2-2.5x"). The wide-run-guard
"unreachable on realistic input" phrasing (spec:1192, comment overlap3.cpp:814)
is hedged ("on realistic input", "No fixture reaches it now") - NOT stated as
absolute proof; the ulp-run reachability adjudication is another lane's. Only
doc-accuracy defect is the "bitwise" overclaim above.

## VERDICT: NEED-CHANGE
Fidelity gate (zero oracle-wrong) SURVIVES; red-first, pin re-anchor, class map,
counts, perf-cap, doc-ASCII all SURVIVE. The single substantive defect: the
"Cray/Offset2/3/4 resolve BITWISE" claim (spec + commit + arc notebook) is
empirically FALSE for Offset2/3/4 (they shift a few ulps; only Cray is bitwise) -
correct the wording to "oracle-correct within eps-bound, only Cray bitwise".
No fence violation; no test-acceptance regression.
