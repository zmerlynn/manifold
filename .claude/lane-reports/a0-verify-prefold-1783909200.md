# A0 verification lane: PRE-CANONICALIZE (pre-fold) separability attack

Mission: attack the stage-A0 NOT-SEPARABLE verdict. A0 measured the signed-mult
channel POST-Canonicalize and found it constant +1. But Canonicalize
(overlap3.cpp:218-309) FOLDS coincident faces by signed mult and DISCARDS the
cancelled pairs (line 293 `if (e.mult == 0) continue;`). PRE-fold the composed input
holds those faces individually. PRIMARY ATTACK: does any PRE-CANONICALIZE channel
(cancelled-coincident-content marker, pre-fold mult spectrum, coincident-face counts,
identity of what cancelled) separate the must-survive lump cluster from the
must-collapse controls, across the junction population, both carriers, fenced both
ways? If yes -> A0 PREMATURE (report SEPARABLE-PRE-FOLD). If pre-fold is also
degenerate -> A0 HARDENED.

Workspace: /tmp/a0-verify-prefold (rsync copy, canonical READ-ONLY except this
notebook). Fresh vbuild. Output CSVs under scratchpad.

## Step 0: grounding (COMPLETE)

Structural facts confirmed from source, not prose:
- Canonicalize (overlap3.cpp:218): merges verts by eps-ball uf; folds tris into
  faceMap keyed by sorted merged-vert triple; `mult += (par==repPar)?1:-1`; then at
  line 293 DROPS `mult==0` faces and at 300 drops zero-area normals. arr.verts ==
  canon.mergedVerts (FindSeams copies all merged verts, line 329-331), so a merged
  vert whose faces ALL cancel survives as a ZERO-INCIDENCE arr.vert.
- THE HOLE: A0 (v4-a0 notebook) found Havoc's 8 right-operand lump verts map EXACTLY
  (distEps=0) to ZERO-INCIDENCE arr.verts (seamDeg=0, nFaces=0), and concluded "no
  incidence -> no information; lump loss is at emission below the table's reach." But
  nFaces==0 post-fold does NOT mean empty pre-fold: every input vert belongs to tris,
  so a zero-incidence arr.vert MUST have had pre-fold content that all cancelled or
  degenerated. The post-fold histogram (all +1) cannot distinguish "no coincidence
  ever" from "coincidence that cancelled to 0 (dropped)". That distinction is the
  candidate separator A0 never measured.
- A0's mult histogram: Havoc 1:176, GT7863 1:456, all +1 (one +2 in openscad). A0
  itself flags: "windings cancel to zero (dropped at Canonicalize) OR stay +1" - it
  ACKNOWLEDGES cancellation may have happened and treats the result as informationless.
- Fixtures: CorpusPairGate loads left+right .obj, ComposeImpl (union-composed),
  eps=EpsilonFromScale(bbox,1000). Havoc/GT7863 are UNION-composed (A0 caveat: a
  subtraction carrier could carry -1; out of A0 scope - noted, will keep in verdict).

PLAN: tap Canonicalize's faceMap BEFORE the mult==0 drop (record nSame/nOpp per key =
pre-fold fold spectrum, including cancelled faces + their normals). Expose via global
+ test hook. Probe: junction grouping (eps-ball uf @ 8*eps over arr.verts, A0 recipe),
per-vert pre-fold aggregate (cancelledFaceCount, coincidentFaceCount, foldCount
spectrum, cancelled-face normals), cross-checked against post-fold arr incidence.
Locate control (dx~0.445 twin) + lump verts (right-operand x~12940-12942). Compare +
full-population dump (secondary attack 2). Env-gated, reverted after.

## Step 1: probe built + STRUCTURAL CORRECTION

Src tap: overlap3.cpp Canonicalize FaceEntry gains nSame/nOpp; after the fold loop,
BEFORE the mult==0 drop, capture ALL faceMap entries (key, nSame, nOpp, mult, rep-tri
normal) + degenerate-at-merge tris + nInputTris into a global PreFoldCapture (env-gated
via SetPreFoldCaptureEnabled). Probe: test/overlap3_test.cpp A0PreFold_Probe
(OV3_PREFOLD=havoc|gt7863|poscancel|poscoin). Build clean; only survivors==arr.faces
asserted (tap faithfulness).

STRUCTURAL CORRECTION found immediately: arr.verts != canon.mergedVerts. arr.verts =
[merged verts 0..nMergedV-1] + [seam CROSSING-POINT verts appended by FindSeams
FindOrAddVert, overlap3.cpp:481]. Havoc nMergedV=92, arrVerts=212 (120 seam verts);
GT7863 nMergedV=236, arrVerts=596 (360 seam verts). So A0's "nFaces/seamDeg" is
SEAM-based incidence (vert -> incident seams -> faceId0/1), NOT triangle membership.
Reproduced A0's seam-based channel faithfully; the control (dx~0.445 twin) is a pair of
near-coincident SEAM crossings (preTotal=0 = no triangle membership, seamDeg=2), while
the lump verts are INPUT triangle verts (preTotal>0, seamDeg=0). Different levels - as
A0 concluded.

## Step 2: MEASURED - the pre-fold channel is EMPTY on both carriers

TAP FAITHFULNESS (airtight): inputTris == preFoldFaces == survivors == arrFaces
EXACTLY. Havoc 176==176==176==176; GT7863 456==456==456==456. The Canonicalize fold
DROPS NOTHING on either carrier - it is a pure BIJECTION (each input tri -> one unique
canonical face).

  carrier | inputTris arrFaces | cancelledFaces coincidentFaces(fold>=2) maxFold degenTris
  havoc   |    176      176    |      0              0                     1        0
  gt7863  |    456      456    |      0              0                     1        0

POPULATION-WIDE (secondary attack 2 - FULL scan, not a sample):
- Havoc nV=212 (92 merged): vertsWithCancel=0 vertsWithCoin=0 vertsWithDegen=0. All 92
  zero-post-incidence verts: 0 cancel, 0 degen.
- GT7863 nV=596 (236 merged): vertsWithCancel=0 vertsWithCoin=0 vertsWithDegen=0. All
  210 zero-post-incidence verts: 0 cancel, 0 degen.
- ZERO markers anywhere. No outlier the small A0 sample missed. Not one merged vert on
  either carrier bears cancelled or coincident pre-fold content.

LUMP vs CONTROL (primary attack, the direct comparison):
- Havoc: 8 right-operand (lump) verts all MERGED, preTotal 4-5, preCancel=0 preCoin=0
  preDegen=0. Control (cid200, dx=0.445eps, 2 seam verts): preTotal=0, preCancel=0
  preCoin=0. IDENTICAL on every pre-fold channel (both all-zero cancellation).
- GT7863: 126 right-sited verts, preCancel=0 preCoin=0 preDegen=0 throughout. Control
  (cid254, dx=0.000eps coincident seam pair): preCancel=0. Identical.

POSITIVE CONTROLS (prove the detector is not just always-zero):
- poscancel = Cube composed with its winding-INVERTED self: inputTris=24, preFold=12,
  survivors=0, cancelledFaces=12, coincidentFaces=12, maxFold=2. ALL 12 faces cancel to
  mult 0 and are dropped - the detector sees exactly the "cancelled coincident content
  the fold discards" the attack posits.
- poscoin = two IDENTICAL cubes (same winding): survivors=12, cancelledFaces=0,
  coincidentFaces=12, maxFold=2. Coincidence WITHOUT cancellation, correctly separated.
The instrument fires on cancellation AND on coincidence; it reads zero on the carriers
because there IS none.

INTERPRETATION: the carriers are UNION-composed; their two operands OVERLAP but share
NO coincident faces (surfaces cross transversally, not tangentially-coplanar). So no two
input tris merge to one face key -> no fold -> no cancellation, no +2, no mult spectrum.
The attack's premise ("Canonicalize discarded cancelled coincident content") is
empirically FALSE for these carriers: nothing was discarded (inputTris==arrFaces). The
Havoc ~1.7-eps twin is a GEOMETRIC near-coincidence >eps apart that never merges into a
cancellable pair - EXACTLY the "coincident content never reaches Canonicalize as
cancellable pairs" hardening outcome. To see the twin pre-fold you would need a
positions/normal near-duplicate-FACE channel = A0's already-measured-and-killed
face-normal-spread / the emission cap-image clustering (the killed classifier). No NEW
pre-fold separator exists.

VERDICT (primary attack): A0 HARDENED. The pre-fold signed-mult / coincidence / mult-
spectrum / degenerate-at-merge channels are ALL degenerate (identically zero) on both
carriers, across the full junction population, fenced by calibrated positive controls.
A0's post-fold "constant +1" is not an artifact of the fold hiding cancellation - the
fold has nothing to hide (bijection). CAVEAT: this is measured on the UNION-composed A0
targets; a SUBTRACTION-composed carrier would carry mult -1 / cancellation (poscancel
demonstrates the channel is real) and is out of A0 scope - so the verdict is
"NOT-SEPARABLE with the channels present on the union carriers", not "no pre-fold
channel can ever exist". (This caveat is in the A0 notebook L132-133 but NOT in the
plan - see secondary attack 4.)

## Step 3: secondary attacks

ATTACK 2 (population coverage) - CLEARED, and STRENGTHENS A0. A0's step-1/2 separability
compared the control against a handful of lump candidates. This lane scanned the FULL
merged-vert population on both carriers (Havoc 92 merged / 212 arr; GT7863 236 merged /
596 arr) for the pre-fold markers: vertsWithCancel/Coin/Degen = 0/0/0 everywhere, incl.
all 92+210 zero-post-incidence verts. No outlier the small sample missed. The pre-fold
degeneracy is UNIVERSAL, not sample-limited.

ATTACK 3 (pin-4b vacuity) - SOUND for the current fixture family; SCOPED, not universal.
The vacuity rests on the dichotomy {isolated twin closes positions-only (table not
load-bearing) vs entangled carrier is table-load-bearing but does not close}. My probe
CONFIRMS the premise for the current fixtures: the isolated twin (A0 measured perptwin
all +1) and both carriers carry NO cancellation channel, so blanking a constant-+1 table
changes nothing - mutation (2) is genuinely vacuous AS SPECIFIED. BUT the vacuity is not
a theorem for ANY constructible fixture: poscancel proves a subtraction/coincident-face
composition carries a REAL mult-0 cancellation channel. A fixture that (i) uses that
channel as the must-survive marker AND (ii) closes oracle-true would red mutation (2) -
the non-vacuous pin. Whether such a fixture is constructible depends on the UNBUILT
stage C (A0 did not build the identity-informed closure), so it stays hypothetical. The
plan's pin-4b language "VACUOUS as specified" is appropriately hedged; but its
justification clause "no fixture is both table-load-bearing AND closing" reads as
universal and should be scoped to the union/positions-only family. CONNECTION to primary
attack: no pre-fold separator on the union carriers => no non-vacuous pin from THEM; the
only route to a non-vacuous pin is a deliberately-constructed subtraction/coincident
fixture, which is exactly where the cancellation channel lives.

ATTACK 4 (plan-edit fidelity) - LARGELY FAITHFUL, two gaps.
- The RECORDED RESULT / ratchet edits correctly attribute NOT-SEPARABLE to the measured
  channels (mult constant, other channels control==lump, jaccard falsified by control,
  lump on zero-incidence verts collapsing at emission). Histogram numbers match. The
  "positions-only spread classifier the escape kill already had" attribution is correct.
  No overreach on the core separability claim - if anything my pre-fold scan makes it
  STRONGER (A0 measured post-fold; pre-fold is also degenerate).
- GAP A (the one the mission flagged): the UNION-COMPOSED scope caveat is in the A0
  notebook (L132-133) but is NOT in docs/V4ImplPlan.md - not in the RECORDED RESULT
  block, not in the ratchet A0 VERDICT ("NOT-SEPARABLE for both A0 targets. THE MUST-MOVE
  SET IS EMPTY" is stated full-stop). poscancel makes this load-bearing: the channel IS
  real for non-union compositions. The plan should carry "NOT-SEPARABLE with the channels
  present on these union-composed carriers; subtraction-composed carriers out of scope."
- GAP B (minor): RECORDED RESULT L361-362 "union/self-overlap windings cancel to zero
  (dropped at Canonicalize) or stay +1" implies cancellation-drop occurred on the
  carriers. MEASURED: nothing cancelled (inputTris==arrFaces, fold is a bijection); the
  true reason for all-+1 is "no coincident faces to fold", not "cancelled and dropped".
  The "or stay +1" hedge keeps it from being wrong, but it conflates the hypothetical
  subtraction case with the measured union case.
- C/D RE-SCOPE value statement ("density wins + match-consistency + machinery
  simplification") MATCHES the round-1 density evidence (R-DENSITY/R-TRIPLE SURVIVEs:
  table bounded at vertex scale, thinning shrinks criticals ~12x/~7x, zero triple
  points) and correctly claims NO carrier resolution. No overreach.

RESIDUAL (noted, out of the commissioned pre-fold scope): the composed input Impl retains
per-triangle runOriginalID (operand label), which Canonicalize discards (reads only
vertPos_/halfedge_). That is a distinct channel from the fold. FACT 1 (closure lane)
already adjudicated it non-disambiguating for union content; my data corroborates -
Havoc's lump verts are ALL right-operand, and GT7863 has 126 right-operand-sited verts,
so "operand==right" cannot separate a lump from the many safe-to-collapse right verts.
Not re-measured; flagged as already-adjudicated.

## WRAP / fence
- VERDICT: A0 HARDENED. Pre-fold channel EMPTY on both union carriers (fold is a
  bijection; inputTris==arrFaces; cancelledFaces=0 coincidentFaces=0 degenTris=0
  population-wide), calibrated by positive controls (poscancel: 12 cancelled; poscoin:
  12 coincident). No pre-fold separator. A0's NOT-SEPARABLE stands and strengthens.
- Secondary: attack2 CLEARED (full-pop scan, no outlier). attack3 vacuity SOUND but
  SCOPED (poscancel = the non-vacuous route via a constructed subtraction fixture).
  attack4 two fidelity gaps (missing union-scope caveat in plan; "cancel to zero
  dropped" imprecise for the carriers).
- Evidence CSVs: scratchpad prefold_{havoc,gt7863}_{pop,verts}.csv.
- FENCE: canonical repo untouched (only this notebook written there). All probe code
  lives in the /tmp/a0-verify-prefold COPY (src tap + test probe), reverted below.
  No git commit.
