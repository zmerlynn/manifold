# Style/simplicity crucible - main-agent generation notebook

Branch explore/sweep-plane-3d-v3 from c96212a5 (arc-closed M4 state).
GOAL: behavior-preserving style+naming+reuse pass per
.claude/skills/manifold-style-review; user directives: (1) house
style, (2) dedicated "let linalg do its job"/reuse lane, (3) NAME
THE STAGES (A/B'/C'/D'/E' are meaningless).
DONE GATE: suite 580+1 green, zero behavior change, lanes converge.

## Plan (generation worklist, before edits)

STAGE NAMES (code + doc + tests, primes die):
- A -> CANONICALIZE (vert merge + face multiplicity):
  StageA/StageAResult -> Canonicalize/CanonicalGeometry-ish.
- B' -> SEAMS (face-pair detectors + seam segments + criticalXs):
  StageBPrime/StageBResult -> FindSeams/SeamsResult-ish.
- C' -> SLABS (BuildSlabs keeps its name; comments renamed).
- D' -> STRIPS (EmitStrips keeps).
- E' -> CAPS (EmitCaps keeps).
- RunCDEPrime -> SweepEmit (matches the doc title).
- docs/SweepEmit3D.md: architecture sections renamed with a legend;
  historical crucible records keep their letters verbatim (they are
  records), one mapping line at the top.
- Test messages "stages A-C" updated.

STYLE (skill items I can see from here):
- for (const int i : {0,1,2}) for fixed tri loops.
- C-style (int) casts -> match type or static_cast (mechanical).
- vertMap (crosses stages) -> origVert2Merged per src2Dst rule.
- static inside anonymous namespace: drop (redundant; match
  boolean2's plain functions).
- Comment hygiene: contracts present (mostly), hunt stale ones.

LINALG/REUSE (my own list; the lane hunts more):
- Box3 hand-rolled in FindSeams -> manifold::Box + DoesOverlap.
- Hand-rolled 2D cross products (ClipPolyByHalfplane,
  SegTriInteriorLen2D, coplanarInteriorOverlap area loop) ->
  la::cross(vec2, vec2).

SETTLED (grounding for lanes): the spec adjudications + audit-round
rebuttals; fail-closed arms are spec-demanded; fence rule (no test
weakening); behavior preservation trumps style (skill's own rule).

## Step 1: generation

DONE, committed d9df4c85, full suite 580+1 green. Delivered:
- Stage rename end-to-end: Canonicalize/CanonicalGeometry,
  FindSeams/SeamsResult, SweepEmit (was RunCDEPrime); banners +
  spec refs (SEAMS/SLABS/STRIPS/CAPS); test messages; doc
  architecture headers renamed with code-entry-point legend,
  records keep letters.
- vertMap -> origVert2Merged (11 sites).
- Box3 struct deleted -> manifold::Box + DoesOverlap (6-way
  hand-rolled test died).
- Hand-rolled 2D crosses -> la::cross at 4 sites (checked
  expression-identical: same term order, bitwise-preserving).
- for (const int i : {0,1,2}) at 4 fixed-tri loops.
- (int) -> static_cast<int> (25 sites, regex + arithmetic-site
  hand-check: cast binds to the chain before / and -, same
  semantics).
- Anonymous-namespace statics dropped (boolean3.cpp precedent:
  zero statics); SweepEmit moved inside the anon namespace.

## Step 2: lanes (round 1)

- Style judge (Codex, /tmp/codex-branch-review/review): skill-file
  as authority, verify pass completeness + hunt misses, cap 3-5.
- Reuse hunter (Codex, review-reuse COPY - relaunched after I
  caught my own one-editor-per-tree violation: the first launch
  told it to edit the same tree the style lane reads; killed
  within seconds, cloned the tree, retargeted, relaunched).
- Naming coherence (Sonnet, read-only canonical): residual letter
  refs, name honesty (does FindSeams' name cover its detector
  duties?), skill naming rules, type-family coherence.
Awaiting verdicts.

## Step 3: round-1 verdicts + validation + fold

STYLE: NEED-CHANGE/5. REUSE: NEED-CHANGE/6 (lane built+tested each
claim in its own copy). NAMING: NEED-CHANGE/6, scheme SURVIVES.

Validation before folding (Phase 4):
- Gate4e/4f "accept success" claim: VERIFIED by reading - the
  assertions are guard-shape (wrong fatal = red) with success
  documented-legal; the MustFail NAME lies. Rename-only fix, no
  assertion touched (fence rule intact).
- linalg equivalences: read linalg.h definitions - lerp = a*(1-c) +
  b*c (component products commute -> bitwise vs (1-t)*a + t*b),
  cross(T, vec2) = {-a*b.y, a*b.x} (x1.0 exact), length2 = dot(a,a),
  .yz() swizzle exists. All five reuse items bitwise ✓.
- ADJUDICATED AGAINST the lanes: (a) style-3's SweepWinding
  reorder - its trailing defaulted optional outputs ARE the
  skill-correct defaultable shape; reordering would break the
  default-arg pattern. Emission-cluster internals reordered
  instead (out, in/out, const-inputs-last). (b) reuse-1 Box2->Rect:
  upstream-landed engine internals, out of this pass's scope;
  recorded as possible upstream follow-up. cnt stays trailing
  everywhere (consistency with BuildSlabs/SweepWinding).

FOLDED: R2 swizzle, R3 cross(1.0,yz) (+comment), R4 length2 x6,
R5 lerp, R6 FaceKey->std::array; S1 cast sweep ((double) +
22 test-file casts), S2 statics (MergeSrcId + test anon-ns),
S3 emission-cluster out-params-first, S4
lPieces->lPiece2RawVerts/rPiece2RawVerts, S5 M1 banner rewrite
(was describing the pre-criticalXs mechanism!) + spec-B' ref +
tolerance->slack + Gate4e/4f -> _FailClosed; naming stragglers
(canon.faces, Spec STRIPS, caps-before-strips comment) + doc
legend accuracy fix (the records contain NO letters - the legend
claimed they keep them; retargeted at SweepPlane3D.md).

Suite after full fold: 580 + 1 skip. Python part-B abort lesson
(again): assert-before-write means a failed match loses the whole
batch - reread the drifted site, re-run.

## Step 4: round 2 (convergence, final per cap)

REUSE: SURVIVE, FINDINGS 0 - the dedicated lane converged; its own
build+filter run as evidence. STYLE: NEED-CHANGE/3.

Validation + disposition:
- F1 (boolean2.cpp reopened anon namespaces at 254/500/630):
  VERIFIED upstream-landed (this branch's boolean2.cpp diff is 10
  lines - the allClosed/edgesNeg threading only). Deferred as an
  upstream follow-up, same class as Box2->Rect. Not re-litigated.
- F2 (seam-cluster out-refs): FOLDED as return values, the skill's
  preferred form - LineTriClip -> optional<interval>, TriTriSeam ->
  optional<segment pair> (equal pair = point contact, contract
  preserved), SeamSeamCrossX -> optional<double>;
  ComputeSectionSegment outputs-first AND its nullable track
  pointer was dead generality (single call site always passes it)
  -> mandatory ref, guard deleted.
- F3 (review-history prose in comments): FOLDED bounded - the
  ExtendPtWithSeams justification narrative and the BuildImpl
  "run-merge mistake this replaced" clause trimmed; invariants and
  spec [fold] pointers KEPT (they are contract pointers, not
  narrative); "(mutation-verified)" tag dropped, the stubbing
  property sentence kept (it is the pin's meaningfulness contract).

Suite after round-2 fold: 580 + 1 skip. CLOSED AT CAP (2 rounds).

## Endpoint

Three commits: d9df4c85 (generation: stage names + house style +
linalg), 5d45de7b (round-1 fold), + the round-2 close commit.
Lane trajectory: style 5 -> 3 -> folded/deferred; reuse 6 ->
SURVIVE/0; naming 6 (scheme SURVIVE) -> folded. Residuals, both
recorded upstream-deferrals: boolean2.cpp namespace structure,
Box2 vs Rect. Process notes: caught my own one-editor-per-tree
violation at launch (reuse lane retargeted to a clone within
seconds); two python assert-before-write batch aborts from
formatter drift (reread the site, re-run).
