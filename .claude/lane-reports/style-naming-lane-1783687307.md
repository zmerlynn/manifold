# Style/naming lane - stage-rename coherence review

Tree: /home/zml/src/manifold-overlap @ d9df4c85 (explore/sweep-plane-3d-v3). READ-ONLY except this file.

Rename under review: A -> CANONICALIZE (Canonicalize/CanonicalGeometry), B' -> SEAMS
(FindSeams/SeamsResult), C' -> SLABS (BuildSlabs), D' -> STRIPS (EmitStrips),
E' -> CAPS (EmitCaps/ComputeCap), driver RunCDEPrime -> SweepEmit.
Doc: SweepEmit3D.md architecture renamed w/ legend; crucible-record/M4/RISKS sections
intentionally keep letters (do NOT flag there).

## Plan

1. Read skill Naming and Vocabulary section (.claude/skills/manifold-style-review/SKILL.md).
2. RESIDUALS: exhaustive grep for stage letters in src/overlap3.h, src/overlap3.cpp,
   src/overlap3_sweep.cpp, test/overlap3_test.cpp + doc normative sections.
   Patterns: stage A/B'/C'/D'/E', A+B', A-C, B/C, C/D, "spec B'", Stage-B, CDE, etc.
3. COHERENCE: read each renamed function body; judge name vs behavior
   (Canonicalize, FindSeams incl. fatal detectors + extra criticals, BuildSlabs,
   EmitStrips, EmitCaps, ComputeCap, SweepEmit driver). Newcomer mental-model test:
   canonicalize -> seams -> slabs -> caps -> strips -> assembly.
4. NAME RULES: srcType2DstType maps (hunt direction-less cross-stage maps), precise
   verbs (merge vs dedup), spelled-out geometry, bool predicates, eps-vs-tol.
5. TYPE/FIELD: CanonicalGeometry vs ArrangementGeometry vs SeamsResult family;
   StageResult<T> comment; any struct/field implying old letters.

Findings format: file:line, incoherence, smallest fix. No churn-for-churn renames.

## Log

### Step 1: skill Naming and Vocabulary - DONE
Read SKILL.md:138-155. Key rules loaded: srcType2DstType for cross-stage remaps
(direction-less only flagged when crossing stages; short loop locals exempt);
no rename-only churn; spelled-out geometric terms; plain bool predicates; past
tense consumed state; merge=geometric vs dedup=topological; pithy test names.
Also Comments rules (contracts expected on large fns, accuracy gate) - relevant
to "should FindSeams' contract carry the fatal-detector behavior".

### Step 2: RESIDUALS grep - about to run
Why: renames leave stragglers in comments/messages. Patterns: [Ss]tage[- ]?[A-E],
[A-E]' prime forms, CDE, spec [A-E], stages [A-E], A+B, B/C, C/D, phase [A-E],
RunCDEPrime, EPrime/CPrime/DPrime identifiers. Single bare letters too noisy;
target contexts. Files: the 4 code files + full doc (then classify doc hits by
section: normative vs settled-historical).

### Step 2 results: residual greps
Pass 1 (stage letters, primes, CDE, spec X, A+B, B/C): code hits at
test/overlap3_test.cpp:618 ("B':"), :945 ("stageA.faces"), src/overlap3.cpp:575
("Spec D'"), :1022 ("E' then D'"). The "B in A's plane" / "tet B's face" hits are
SOLID A/B (two input solids), not stages - excluded.
Pass 2 (identifiers, prime forms): test/overlap3_test.cpp:607 "StageBPrime" in a
comment - the function is now FindSeams; no identifier residuals in src (all new
names present: Canonicalize/CanonicalGeometry, FindSeams/SeamsResult, BuildSlabs,
EmitStrips, EmitCaps/ComputeCap, SweepEmit, ArrangementGeometry, StageResult).
Doc hits: intro lines 3-11 reference predecessor doc SweepPlane3D.md's stages
(D/E, A-C) - cross-doc reference, old doc still uses letters; legend at 22-23;
"(formerly X)" tags at 26/29/35/52/61/89 = the legend approach; line 182
"stage-D failure" is in "What dies" (normative!) - need to read context: it
names the OLD SweepPlane3D stage D that died, probably legitimate.
Next: read contexts of the 5 code hits to confirm staleness + smallest fixes.

### Step 2 conclusions: 5 confirmed code residuals (all comments, no identifiers)
R1 test/overlap3_test.cpp:607 "arr produced by StageBPrime" - fn is FindSeams.
   Fix: s/StageBPrime/FindSeams/.
R2 test/overlap3_test.cpp:617-618 "(spec B': triple criticals...)" - doc section
   is SEAMS now; nearby lines 615/628 already say "seams stage". Fix: s/spec B'/spec SEAMS/.
R3 test/overlap3_test.cpp:945 "-> stageA.faces empty ->" - no stageA identifier
   exists; the canonicalize output local is `canon` (CanonicalGeometry.faces).
   Fix: s/stageA.faces/canon.faces/ (line 942 comment already says "canonicalize stage").
   Also lines 942-944 have a bad rewrap artifact ("merges the\nduplicate positions") -
   trivial, from the rename edit.
R4 src/overlap3.cpp:575 "Spec D' (R2-fold)" - line 587 shows post-rename convention
   "spec CAPS [R2-fold]". The rule is the strip constant-extension rule. Fix:
   s/Spec D'/Spec STRIPS/.
R5 src/overlap3.cpp:1022 "E' then D': caps run first" - banner at 996 already uses
   new vocab (slabs -> caps -> strips -> assembly). Fix: s/E' then D'/CAPS then
   STRIPS/ (or drop the prefix; "caps run first" carries it).
Doc line 182 "stage-D failure" still to classify (What dies section). Intro lines
3-11 reference SweepPlane3D.md's own letters - cross-doc, coherent, not flagging.

### Step 3: COHERENCE - about to read overlap3.h fully then each renamed fn body
Why: judge name-vs-behavior for Canonicalize/FindSeams/BuildSlabs/EmitStrips/
EmitCaps/ComputeCap/SweepEmit + newcomer mental-model test + contract comments.

### Step 3 progress: overlap3.h + overlap3.cpp stage bodies read
Header (all 213 lines): clean of letters. Banners use stage names
("Canonicalize-stage output", "Slabs-stage types", "Seams-stage types",
"Arrangement geometry (canonicalize + seams output; slabs, caps, and strips
input)"). StageResult:67 comment "Per-stage result: either a product or a fatal
reason" - letter-agnostic, fine as briefed. FaceTrack "(strips stage)",
SeamTrackEntry "(spec STRIPS/CAPS)", ArrangementGeometry "spec SEAMS" - new
vocab used consistently. Header type order puts Slabs-stage types before
Seams-stage types (dependency-driven: MergedVert/Seam sit next to
ArrangementGeometry) - ordering not naming, no flag.
- Canonicalize (202-308): UF eps-merge of verts + face parity/multiplicity
  accumulation, drops mult==0 + degenerate. Name honest; banner carries detail.
  origVert2Merged: direction-ful, good. "merge" used for eps-scale = correct verb.
- FindSeams (310-519): banner 310-314 explicitly carries the extra duties
  ("CoplanarOverlap + EdgeInPlane detection, seam segments, vertex-free extra
  criticals"). Fatal channel visible in SeamsResult type. Judgment: name honest
  WITH the banner; detectors are pair-analysis byproducts of the same face-pair
  loop that finds seams - splitting the name would misdescribe the fusion.
  "spec SEAMS" at 496. OK.
- EmitStrips (669-705): zips chains, no geometry computed - "emit" honest;
  contract present. EmitCaps (887-922) driver vs ComputeCap (809-885) one cap:
  considered flagging ComputeCap vs EmitCap - rejected: EmitCap/EmitCaps would
  be a one-char confusable pair, and ComputeCap does arrangement+chains too.
- SweepEmit (999-1051): banner "Pipeline: slabs -> caps -> strips -> assembly"
  = exact contract of what it runs; RemoveOverlaps3D reads Canonicalize ->
  FindSeams -> SweepEmit. Mental-model composition correct. Driver name covers
  sweep(slabs)+emission(caps/strips) = old C'+D'+E' scope, honest.
- NEW FINDING (verb precision, skill 152-153): BuildImpl comments call the
  eps-scale GEOMETRIC vert merge "dedup": src/overlap3.cpp:925 "dedup verts",
  :932 "Collect verts with eps-dedup" (body: la::length <= eps merge). House
  vocab elsewhere says weld/merge (901 "assembly weld", 948 "after the weld").
  Fix: s/dedup verts/merge verts/ (925), s/eps-dedup/eps-merge/ (932). The TRI
  dedup at 942-954 (seenTris exact set) is correctly "duplicate"/topological.
- CapEdgeSet:766 "NO vertex pre-dedup ... MergeVerts owns snapping" - "pre-dedup"
  refers to exact-duplicate removal it deliberately does NOT do; borderline but
  the sentence contrasts with (geometric) snapping correctly. Not flagged.

### Step 4: about to read overlap3_sweep.cpp (BuildSlabs) + overlap3.cpp:1-195
Why: BuildSlabs coherence + hunt direction-less cross-stage maps in helpers.

### Step 4 results: BuildSlabs + helpers + test vocabulary
- overlap3_sweep.cpp read fully: signpost 15-17 "Slabs stage ... Design:
  SweepEmit3D.md, SLABS"; banner 93-95 carries duties (x-criticals, sections,
  engine calls, face tracks). BuildSlabs name honest (builds SlabResult vec;
  fatal guards documented). vertIdx map (135) = slab-local exact-key position
  map, loop-local -> exempt per skill. trackOf (overlap3.cpp:779) local too.
- overlap3.cpp:1-195: signpost line 15-16 states the FULL new pipeline order
  "canonicalize -> seams -> slabs -> caps -> strips -> assembly" - newcomer
  mental-model anchor present at top of the main file. Helpers under
  "Seams-stage helpers" banner; geometric names spelled out (TriTriSeam,
  SeamSeamCrossX, ClipPolyByHalfplane, SegTriInteriorLen2D); tolerances derived
  from eps named kTol/parTol - eps-vs-tol vocab OK.
- ZipperEmit/ChainSplitVerts coherent; comments use spec CAPS/STRIPS refs.
- 'stage' word audit across all 4 files: every hit uses new stage nouns
  (canonicalize/seams/slabs/caps/strips stage) - no orphan letters beyond the
  5 residuals already logged.
- Test names: Gate1-5, Pin_M1/P4/P6/P7/P12/S3a/S3b, EmissionAlgebra - gate/pin/
  milestone tags, not stage letters. Pithy per skill. No flag.

### Step 5: about to read doc normative sections (1-233) + line 182 context
Why: classify doc letter hits (intro 3-11, legend 22-23, formerly-tags, 182
"stage-D failure" inside "What dies" which IS normative) and check the legend
maps every letter used.

### Step 5 results: doc normative sections read (1-233) + settled sections scanned
- Intro 3-14: letters there reference the PREDECESSOR doc ("stage-D/E
  architecture of docs/SweepPlane3D.md", "stages A-C's proofed mechanics",
  "stage D had no 2D analog") - cross-doc, attributed, coherent. Not residuals.
- Legend 18-23 maps letters; "(formerly X)" tags at 25/29/52/61/89 are the
  legend approach - intentional per briefing.
- Line 182 "the whole stage-D failure taxonomy (ClassificationAmbiguity, ...)"
  in What dies (normative): names the PREDECESSOR's dead taxonomy by the
  predecessor's term; the enumerated reasons disambiguate. Not flagged (the
  dead thing has no new name; qualifying it would be churn).
- Eps posture 210 says "at CANONICALIZE" (was "at stage A" pre-rename) - good.
- Output and tests: clean.
- FINDING (doc): legend sentence 20-23 claims "The crucible records later in
  this doc predate the names and keep their original stage letters" - FALSE
  empirically. grep of lines 234-452 (RISKS + all crucible/M4/post-M4): zero
  stage letters, zero prime marks, zero old uppercase names (MERGE/SECTIONS/
  CRITICALS AND SEAMS: 0 hits); they use R1-R4/L1-L3/M4 ids, and 5 hits of the
  NEW code names (EmitStrips, ArrangementGeometry...). Verified pre-rename too
  (c96212a5: letters only ever lived in intro + architecture headings + eps
  posture). The letters the legend actually serves are the intro's
  SweepPlane3D.md references. Smallest fix: retarget one sentence, e.g.
  "The predecessor doc (SweepPlane3D.md) and older records use the original
  stage letters: A = CANONICALIZE, ...". Accuracy gate: inaccurate > absent.

### CONCURRENT DRIFT observed mid-review
src/overlap3.cpp modified under me during the review (git diff vs d9df4c85):
la::length2/lerp adoption, FaceKey -> std::array, lPieces/rPieces ->
lPiece2RawVerts/rPiece2RawVerts, and MY VERB FINDING PRE-EMPTED: "dedup verts"
-> "weld verts", "eps-dedup" -> "eps-weld", "pre-dedup" -> "pre-merge",
"tolerance" -> "slack" x2. Re-verified all letter residuals in the CURRENT
tree: all 5 still present (src lines shifted 575->569, 1022->1016; test lines
unchanged). Doc untouched. test:1085 "BuildImpl deduplicates the resulting
identical triangles" checked - exact-duplicate topological drop, correct verb,
no flag.

### Step 6: residual-hunt completeness + minor observations closed out
- Uppercase stage-name grep across code: only new names (spec SEAMS/STRIPS/
  CAPS, SLABS) - no MERGE/SECTIONS/CRITICALS.
- Cross-stage maps: origVert2Merged, inputVert2Merged, rootToMerged direction-
  ful; vertIdx/trackOf/faceMap/seenTris loop-local (exempt). CanonicalFace.id
  (h:39): direction-less back-ref to orig Impl AND unconsumed (write-only:
  cf.id=e.repTri, no reader) - noted as minor; rename-to-repTri only if kept,
  deletion call belongs to a dead-code/simplify lane, not naming.
- Bools: built/ok/canonical/loopsClosed/trisOk/hasCoverage - plain predicates,
  good. eps-vs-tol: kTol/parTol derived tolerances, fine.
- Considered + NOT flagged (churn or wrong-lane): ArrangementGeometry name
  overloads "arrangement" vs the 2D cap arrangements (and the design pointedly
  builds NO 3D arrangement) - but header comment "canonicalize + seams output;
  slabs, caps, strips input" disambiguates, `arr` pervasive, rename = big churn;
  SeamsResult could structurally be StageResult<ArrangementGeometry> (simplify
  lane, not naming); ComputeCap vs EmitCap (EmitCap/EmitCaps one-char
  confusable - current split is right); header type-block order (dependency-
  driven); doc presents STRIPS before CAPS vs execution caps-first (narrative
  order; the normative dataflow paragraph and code banners state the true
  order); test P6 comment bad rewrap at 942-944 (formatting artifact, noted
  inside finding R3).
- Newcomer mental-model test PASSES: overlap3.cpp:15-16 signpost states
  "canonicalize -> seams -> slabs -> caps -> strips -> assembly"; SweepEmit
  banner "slabs -> caps -> strips -> assembly"; RemoveOverlaps3D reads
  Canonicalize -> FindSeams -> SweepEmit. FindSeams' extra duties (fatal
  detectors, extra criticals) carried by its banner 310-314 + SeamsResult.fatal
  in the signature - name honest with contract, detectors are byproducts of
  the same face-pair loop.

## Verdict

FINDINGS: 7 at d9df4c85 (6 still open in working tree - the verb cluster was
fixed under me mid-review):
1. test/overlap3_test.cpp:607 StageBPrime -> FindSeams
2. test/overlap3_test.cpp:617-618 "spec B'" -> "spec SEAMS"
3. test/overlap3_test.cpp:945 "stageA.faces" -> "canon.faces"
4. src/overlap3.cpp:569 (575@d9df4c85) "Spec D'" -> "Spec STRIPS"
5. src/overlap3.cpp:1016 (1022@d9df4c85) "E' then D':" -> "CAPS then STRIPS:"
   (or delete the prefix - "caps run first" already says it)
6. [fixed concurrently] dedup-for-geometric-merge verbs (925/932/766@d9df4c85)
7. docs/SweepEmit3D.md:20-23 legend justification clause is factually wrong
   (crucible records contain no letters); retarget at SweepPlane3D.md/history

NEED-CHANGE - but scheme-level SURVIVE: the five names + driver name are
honest, the mental model composes correctly, vocabulary is consistently
converted everywhere except 5 straggler comment lines (3 of 5 in the test
file, which the rename pass under-swept) and one false legend claim in a
normative section. All fixes are one-line.
