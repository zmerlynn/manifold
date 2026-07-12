# WALL-B verification lane (fidelity/build) - lab notebook

Adversarial audit of the wall-B landing on explore/sweep-plane-3d-v3.
Artifact: commit range de01ee24..2224ead1, landing f0dab0fb ("retire EngineIdConflict").
Auditing the generation lane's notebook (.claude/lane-reports/wallb-arc-1783861900.md).
Default skeptical; try to BREAK the green claims.

Workspace: /tmp/wallb-verify-fidelity (isolated rsync copy, FRESH build vbuild/).
Canonical repo is READ-ONLY. This notebook lives in the canonical repo (file write only,
no git commit; orchestrator commits).

## Plan

- CLAIM 1 Suite green: build fresh, run Overlap3 (expect 54 pass + 1 skip Gate4c on
  SubEpsFeature chain-plane guard), Boolean2 unchanged, CrossSection 2D if present.
- CLAIM 2 Mutation re-verify: re-arm the fatal at BuildSlabs conflict site -> Gate4c must
  change (skip/fail on conflict-shaped msg); coplanar oracle pins must STAY green.
  If nothing changes, pins are decorative -> BREAK. Revert after.
- CLAIM 3 Resolving-volume fidelity: Cray pair + Offset2/3/4 bitwise per corpus tests.
- CLAIM 4 Hull real-content: instrument chain-plane guard, confirm it fires on REAL macro
  cap content (nonempty minus-side cap edges), not noise. If noise -> BREAK.
- CLAIM 5 Test-contract narrowing audit: git diff de01..2224 test/overlap3_test.cpp;
  every change must NARROW (drop EngineIdConflict) or be a sanctioned change. Scrutinize
  Gate4c skip: is SubEpsFeature skip gated on the guard's detail message or would ANY
  SubEpsFeature silently skip? NEED-CHANGE if loose.
- CLAIM 6 Docs/style: [WALL-B] section - no raw counts, ASCII only, claims match measured.

## Log

### Setup (DONE)
- rsync copy at /tmp/wallb-verify-fidelity, HEAD=2224ead1. FRESH build vbuild/ (Release,
  -ffp-contract=off -fexcess-precision=standard). Confirmed isolated: CMakeCache
  manifold_SOURCE_DIR=/tmp/wallb-verify-fidelity; touch src/overlap3_sweep.cpp ->
  recompiled that TU from the isolated path + relinked. Build clean -j4.

### CLAIM 1 suite green (SURVIVE)
- about-to-try: run Overlap3.* + Boolean2.* + CrossSection.* in fresh build.
- observed: Overlap3 = 54 PASS + 1 SKIP. The skip is Gate4c_HullMask_MustResolve, detail
  "Gate4c hull: chain-plane wide-run guard: macro cap content over a skipped run wider
  than eps". Boolean2 36/36. CrossSection 75/75. Zero failures.
- conclusion: matches the lane's claim exactly. SURVIVE.

### CLAIM 2 mutation re-verify (SURVIVE)
- about-to-try: re-arm the retired fatal at the BuildSlabs conflict site (return
  Fatal{SubEpsInput,"MUTATION re-arm conflict site"} on conflictCount>0 - SubEpsInput is
  NOT in Gate4c's skip set, so the test must FAIL not skip). Rebuild, run Gate4c + the two
  coplanar pins + Gate4d/e/f.
- observed: Gate4c FAILED "Gate4c hull MUST RESOLVE but got fatal=0 MUTATION re-arm
  conflict site" -> the hull DOES reach the conflict site; demotion is load-bearing.
  Coplanar_SameOriented_Oracle OK, Pin_TouchingDisjoint OK (both assert !fatal; they pin
  the flap through coplanar grouping, never conflict). Gate4d/e/f OK (they hit
  SubEpsFeature via a different path, not the conflict site).
- conclusion: pins are NOT decorative; re-arming changes Gate4c and leaves the pins green.
  Mutation reverted, workspace matches HEAD. SURVIVE.

### CLAIM 3 resolving-volume fidelity (SURVIVE)
- about-to-try: instrument (in-tree corpus tests don't print volumes) - add VL_ tests
  mirroring CorpusPairGate/CorpusSingleGate, print %.17g volume of impl Property and
  reconstructed Manifold.
- observed (bitwise, both metrics agree):
  Cray    implVol=oursVol=1.5760799710533451e+116  (== expected)
  Offset2 implVol=oursVol=209026.71902485067       (== expected)
  Offset3 implVol=oursVol=10803.472688795122       (== expected)
  Offset4 implVol=oursVol=15240.582916271358       (== expected)
- conclusion: all four resolve (no fatal) at the recorded volumes, bit-for-bit. (These
  cases never conflict, so demotion cannot move them; confirmed empirically.) VL_ tests
  reverted. SURVIVE.

### CLAIM 4 hull real-content (SURVIVE)
- about-to-try: instrument the chain-plane guard (overlap3.cpp:860) - print run width,
  eps, capEps, plus/minus edge counts, rawVerts, bbox + max edge length of surviving
  content. Run the hull (Gate4c).
- observed: exactly ONE wideRun in the whole hull run, and it fires the guard:
  [WB-RUN]  xCap=-99.6843 runWidth=8.008954e-10 eps=3.519281e-10 ratio=2.28
  [WB-GUARD] xCap=-99.6843 eps=3.519e-10 capEps=2.815e-09 plus=0 minus=14 rawVerts=396
             bboxYZ=[-144.35,-135.15]x[142.53,150.4] maxEdgeLen=4.4925 (=1.6e9 x capEps)
  Reproduces the lane's [WALLB-CAP] plusEdges=0 minusEdges=14 rawVerts=396 exactly. The
  run is genuinely wider than eps (2.28x); the surviving minus-side cap content is MACRO
  (edges up to ~4.5 units, ~1.6e9 x capEps), not sub-8*eps noise.
- conclusion: guard fires on REAL macro cap content over a real wide run. Adjudication
  sound. Instrumentation reverted (all sources match HEAD). SURVIVE.

### CLAIM 5 test-contract narrowing (NEED-CHANGE, scoped)
- git diff de01..2224 -- test/overlap3_test.cpp: five OR-list edits all DROP
  EngineIdConflict only (Gate4d/e/f, RazorBand, PerpFaces) -> strictly NARROWING (fewer
  accepted fatals = STRICTER). Gate4c skip swap Engine->SubEpsFeature is the sanctioned
  change. No assertion deleted, no acceptance widened; no net TEST added/removed in range
  (TEMP_GT7081_Probe added in f0bb7357 + removed in f0dab0fb net-cancels).
- FLAG: the Gate4c SubEpsFeature skip is gated on the BARE ENUM, not the detail message.
  SubEpsFeature has TWO fire sites (BuildSlabs single-face guard "face in unbuilt x-range
  ..."; ComputeCap chain-plane guard "macro cap content over a skipped run wider than
  eps"). The retired EngineIdConflict had ONE fire site, so the old skip was
  site-specific; the new one is NOT. A future regression that moves the hull to the
  single-face SubEpsFeature guard would SILENTLY SKIP (masked), and the test does not pin
  the landing's own specific "lands on the chain-plane guard" claim.
- severity: LOW (Gate4c is aspirational-skip; the detail IS printed; claim-2 proved the
  "fail on unexpected fatal" arm works). But the fix is one line and strictly better:
  gate the skip on result.detail containing "macro cap content over a skipped run wider
  than eps". Recommend tightening -> NEED-CHANGE for this sub-point only.

### CLAIM 6 docs/style (SURVIVE)
- WALL-B section (docs/SweepEmit3D.md) + the entire added diff: ASCII-clean (no unicode
  dashes/arrows/quotes/ellipsis). No stale-prone raw counts - digit tokens are dates,
  sentinel -1, multiplicity +2, spec refs [D2-4]/mechanism-5, algorithm constant 8*eps,
  fixture name 7081. Describes by kind + magnitude, per the no-counts convention.
- claims cross-checked against measurements: chain-plane landing (claim 4), Gate4c skip
  narrowing (claim 5), suite green (claim 1), Boolean2 unchanged (claim 1 + boolean2
  sources untouched in range). srcId-no-3D-reader claim independently grep-confirmed:
  overlap3_sweep.cpp:170-171 is the sole WRITE (group-id seeding); overlap3.cpp has ZERO
  srcId/sourceId occurrences. SURVIVE.

### VERDICT: NEED-CHANGE (scoped to Gate4c skip tightness)
Claims 1,2,3,4,6 SURVIVE empirically. Claim 5: the narrowings are sound and stricter, but
the Gate4c SubEpsFeature skip is gated on the bare enum (2 fire sites) and does not pin
the landing's specific chain-plane claim; recommend a one-line detail-substring tightening.
No correctness BREAK: the suite is green, resolving volumes are bit-identical, the guard
fires on real macro content, and the demotion is mutation-verified load-bearing.
Workspace left at /tmp/wallb-verify-fidelity, all sources reverted to HEAD, clean binary.
