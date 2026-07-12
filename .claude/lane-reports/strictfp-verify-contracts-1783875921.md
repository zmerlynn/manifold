# STRICT-FP SLAB GATE - adversarial verification (contracts / spec / simplicity)

Lane: judgment layer. Artifact = strict-FP slab-gate adoption, branch
explore/sweep-plane-3d-v3, commits 89f52a6f..ef10175e (landing 94eb88df gate+tests,
0075c473 spec). Other lanes own empirical fidelity + the ulp-run reachability analysis.

GROUNDING READ (done): full diff 89f52a6f..HEAD; strictfp-arc notebook;
docs/SweepEmit3D.md STRICT-FP + chain-plane/[R3]/WALL-B/corpus sections;
chain-plane-1783819552 (width-guard original adjudication); wallb-arc (hull
SubEpsFeature real-content verification; budget calibration); subeps-slab-probe
(the corrected dissolution prediction). Did NOT re-run heavy fixtures - my
strongest finding is an INTERNAL CONTRADICTION provable from the artifact itself;
empirical fidelity is another lane's. All changed files verified ASCII-clean.

VERDICT: NEED-CHANGE. The core adoption is honest and the machinery-delta-zero
finding (correcting the probe) is well-supported and prominently recorded. One
clear defect: the budget rationale COMMENT is stale, contradicted by the arc's own
class map. Plus soft coherence gaps: WALL-B and [R3] spec sections left with
present-tense stale statements (no supersession note, unlike chain-plane/corpus),
and the spec perf paragraph is thin/understated relative to the notebook.

---

## ATTACK 1 - REFUSAL-QUALITY REGRESSIONS

### 1a hull: SubEpsFeature-on-verified-real-content -> ArrangementBudget
Information preserved, weld loose. The wallb-arc VERIFIED the hull's SubEpsFeature
fired on REAL macro cap content (14 nonempty minus-side cap edges over a wide run,
after 8*eps noise annihilation), and the spec [WALL-B] section (lines 1104-1106)
preserves that as historical fact. The STRICT-FP section (line 1175) explains the
re-label: "dense near-coplanar bands retain past the ceiling once every ulp slab
builds." So BOTH halves exist (finding preserved; re-label explained).

BUT the causal thread is implicit. ArrangementBudget is a RESOURCE diagnosis that
does not name the geometry; the spec never explicitly says "the same WALL-B-verified
real macro cap content is now SECTIONED rather than refused, and the section of that
content is what trips the ceiling." A reader must infer that "dense near-coplanar
bands" == the verified real content. Information at the refusal SITE is lost (budget
says "too dense," not "macro cap transition here"), but recoverable from the spec.
RECOMMEND: the STRICT-FP hull line weld the two - "the WALL-B-verified real macro
cap content is now processed into a section too dense for the ceiling." Minor.

### 1b GT7863/openscad: SubEpsFeature (specific) -> NonManifoldEmission (generic)
Progress, adequately pointed. This is genuinely more faithful: the in-run content
is now SECTIONED and processed to a real emission failure (wall-A sheet fan) rather
than pre-refused. The thin-band -> processed -> wall-A story IS preserved in the
test comments (Corpus_OpenscadNonmanifold: "the former sub-eps runs build, so the
section runs to completion and lands on the same wall-A sheet fan... not the old
wide-run guard") and the STRICT-FP class map.

Honest trade to note: the RUNTIME detail string regresses from content-specific
("macro cap content over a skipped run wider than eps") to generic ("unresolvable
sheet contact"). A debugger hitting a live NonManifoldEmission on GT7863 learns the
in-run origin only from the test/spec, not the message. Same generic string
Offset1/self already emit, so not a NEW runtime regression - acceptable, framed as
progress. No change required.

## ATTACK 2 - GATE4C TRAJECTORY (third identity)

SURVIVE. Skip discipline held. The primary skip is PINNED to the detail substring
"retained section content exceeds budget" (test:531), which byte-matches the
ArrangementBudget detail in overlap3_sweep.cpp:246 - the fidelity lane's
pin-to-detail demand survived the rewrite. The comment tells the truth about WHY
(too dense to section; real fix is arrangement robustness) and what regresses
("Anything else is a regression").

One pre-existing looseness (NOT this arc's): the SECONDARY accept
`if (result.fatal == NonManifoldEmission) GTEST_SKIP` is unpinned (accepts any
NonManifold detail). It was present verbatim before this arc and is forward-looking
("the other honest boundary a future arrangement fix could move it to"). The current
landing (ArrangementBudget) is the pinned one. Not a NEED-CHANGE for this arc; worth
a future tighten if the hull ever actually lands on NonManifold.

## ATTACK 3 - C3 NON-RECALIBRATION  *** STRONGEST FINDING ***

The re-label decision is honest; the CODE COMMENT justifying the floor is STALE.

Decision honesty: self_A/B NonManifoldEmission -> ArrangementBudget is defensible.
The counter-argument holds: relaxing the budget would run them to the emission
blowup and STILL fail NonManifold, strictly worse and slower. Spec + notebook both
state this. Fine.

Information loss (real but bounded): self_A/B used to demonstrate FULL-PIPELINE
TERMINATION (they completed BuildSlabs, ran the whole section, hit NonManifold at
emission). They now refuse in BuildSlabs, so that end-to-end-termination datum is
masked behind the ceiling. Also, spec line 1201 lumps them under "fail closed
DEEPER" - imprecise: for the self-intersection PAIR the budget trips DURING slab
build, an EARLIER pipeline stage than the emission where NonManifold fired before.
"Deeper" is right for GT7863/openscad (they reach emission now), wrong-signed for
self_A/B (they fail earlier). Minor wording.

STALE BUDGET RATIONALE = NEED-CHANGE. overlap3_sweep.cpp:131-134 (UNTOUCHED by this
arc - diff confirms) still reads:
  "FLOOR 4M: clears the legitimate corpus maximum (~1.8M pieces on the
   self-intersection pair, which COMPLETE this stage) with margin ...
   A 1M floor false-trips the self-intersection pair - do not lower it."
The arc's OWN class map says the self-intersection pair now retains >4M and TRIPS
this budget under strict-FP. So:
- "clears the legitimate corpus maximum (~1.8M ... which COMPLETE this stage)" is
  now FALSE - the pair no longer completes; it trips.
- "A 1M floor false-trips the self-intersection pair" - the 4M floor now trips them
  too, but as an INTENDED true trip, not a false trip. The whole false-trip framing
  is obsolete for the pair.
This is an internal contradiction between a load-bearing calibration comment and the
landed behavior; the mission predicted it exactly. The spec [WALL-B] line 1131 even
points readers here ("see kPieceBudget for the calibration"), propagating the stale
story. SUGGESTED REWORDING: under strict-FP the self-intersection pair's dense
near-coplanar bands all section, so it retains >4M and legitimately TRIPS this
budget (an honest refusal, not a false trip). Re-justify the 4M floor against the
current completing maximum (the clean/medium cases that still complete, e.g. the
Offset pairs), and keep GT7081's ~30M uncapped need as the upper landmark.

## ATTACK 4 - SPEC / NOTEBOOK COHERENCE

### 4a perf 2-12x probe vs ~1.02x landing
Notebook explains it; spec does not, and understates. The notebook Step 7 gives the
mechanism: "self_A/B early-budget-trip speedups offset the hull/openscad slab-count
slowdowns" - i.e. the ~1.02x corpus subtotal is an ARTIFACT of the budget truncating
the heavy cases (self_A 33.6->12.5s, self_B 29.3->11.8s, FASTER by refusing sooner),
not a sign the work didn't grow. The probe's 434s/12x was measured WITHOUT the budget.

The SPEC perf paragraph (lines 1208-1210) does NOT carry this - it just says
"building 2-2.5x more slabs (clean cases) up to ~1.5x more (dense) ... measured
within the suite's perf budget." Two problems: (i) no budget-truncation caveat, so a
spec-only reader who saw the probe's 12x cannot reconcile; (ii) "~1.5x more (dense)"
UNDERSTATES - the completing dense cases measure 2.2-3.8x more slabs (Offset1 3.8x,
GT7863 3x, Havoc 2.7x, openscad 2.2x per the notebook table); the only cases near
~1.5x are the budget-truncated ones (they never build the full explosion). Minor
NEED-CHANGE: spec perf paragraph should note the budget truncates the heavy cases
(that is why the subtotal stays flat) and correct/clarify the "~1.5x" figure.

### 4b premise correction (dissolution prediction)
SURVIVE. The spec records it squarely: STRICT-FP "WHAT DISSOLVES. Nothing, and this
CORRECTS the sub-eps probe's prediction that the run machinery would delete" (line
1183), then itemizes chain-plane-rule/wide-run-guard/non-canonical-skip/coverage as
all kept-and-why. The probe's dissolution plan lived in its recommendation
("WHAT WOULD DISSOLVE if adopted: the CHAIN-PLANE RULE ... wideRun guard ... coverage
guard ... non-canonical skip"); the living doc (spec) now carries the refutation
prominently. The probe NOTEBOOK is append-only history (correctly not rewritten).
Caveat outside this arc's scope: the memory pointer that references the probe should
note the refutation so a memory-first reader doesn't re-derive the dissolution plan.

### 4c ASCII / magnitudes / overclaims
ASCII: clean (grep -P '[^\x00-\x7F]' empty across spec/code/tests/notebook).
Magnitudes: "< 0.001 eps" consistent with notebook's measured <=0.0007 eps. Good.
Overclaim - MINOR: "unreachable" is hedged ("on realistic input", "No fixture reaches
it now", mechanism "~1000 consecutive adjacent-double criticals") but not explicitly
OBSERVED-NOT-PROVEN. The notebook is more careful ("empirically UNREACHABLE ... Not
provably impossible"). Since the ulp lane is actively testing reachability, the spec
+ ComputeCap comment should mirror the notebook's hedge (observed across the measured
corpus; not proven impossible) rather than the near-assertive "unreachable on
realistic input." Wording tighten, not a blocker.

## ATTACK 5 - SIMPLICITY VERDICT

SURVIVE with a coherence caveat. Net machinery delta is genuinely zero (one gate
predicate changed + comments; tests net -125 lines by deleting RingedBox + its pin).
The two candidate smells are both HONESTLY DISCLOSED:
- Retained-but-unreachable wide-run guard: kept as the chain-plane rule's fidelity
  backstop, no fixture reaches it. The notebook lists this as an OPEN RISK in plain
  terms ("retained on a fidelity argument, not a test"). The guard + rule are a
  coupled pair (rule interpolates across runs; guard catches when interpolation would
  be oracle-wrong); keeping both is coherent, not sediment.
- Retired-pin reanchoring: Pin_ChainPlaneRule (eps-scale RingedBox wide run, which
  vanishes under strict-FP) DELETED; the rule now pinned by PerpFaces at ulp scale,
  mutation-verified (revert ZipperEmit to slab bounds reds it). Legitimate move of a
  pin to the fixture that can still exercise the mechanism - a simplification.

COHERENCE CAVEAT (soft NEED-CHANGE): the spec was amended UNEVENLY. Chain-plane and
Corpus sections got forward-pointing supersession NOTEs; WALL-B and [R3] did NOT,
and both now contain present-tense statements the arc falsified:
- WALL-B (lines 1102-1110): hull "lands on the CHAIN-PLANE wide-run guard
  (SubEpsFeature)" and Gate4c "pinned to the chain-plane guard's detail (any other
  SubEpsFeature must fail)" - but the current hull lands on ArrangementBudget and the
  current Gate4c pins to the budget detail (test:531). Stale without a note.
- [R3] DEGENERATE SLABS (line 152): "width <= eps: no section is built" - literally
  false under strict-FP (only ulp-wide slabs skip). A top-to-bottom reader hits this
  ~1000 lines before the STRICT-FP correction.
A fresh reviewer CAN reconstruct the current slab policy (strict-FP floor + ulp-run
regime + budget) from the layered sections, but only by trusting that the LAST
section wins and mentally superseding the earlier present-tense claims. Under the
project's append-with-supersession-note convention this is acceptable IF the notes
are applied consistently; here two sections the mission named were left un-noted.
RECOMMEND: add one-line supersession notes to WALL-B (hull residual + Gate4c pin ->
ArrangementBudget) and [R3] DEGENERATE SLABS (gate is now the FP floor), mirroring
the chain-plane/corpus treatment.

---

## TALLY

NEED-CHANGE (must):
- A3: budget rationale comment overlap3_sweep.cpp:131-134 stale (self_A/B "COMPLETE
  this stage @ ~1.8M" contradicted by the arc's own >4M-trips class map). Reword.

NEED-CHANGE (soft / coherence):
- A5 + A1a: add supersession notes to WALL-B (hull -> ArrangementBudget; Gate4c pin)
  and [R3] DEGENERATE SLABS (gate = FP floor now); weld STRICT-FP hull line to the
  WALL-B real-content finding.
- A4a: spec perf paragraph - add budget-truncation caveat, correct "~1.5x (dense)"
  (measured 2.2-3.8x on completing dense cases).
- A4c: soften "unreachable on realistic input" to observed-not-proven (ulp lane is
  testing reachability).
- A3: "fail closed DEEPER" mis-signed for self_A/B (they fail EARLIER in stage).

SURVIVE clean: A1b (honest progress, story preserved in tests), A2 (skip pinned to
detail, discipline survived), A4b (probe dissolution refutation recorded in spec),
A5 core (machinery-delta-zero, one design, smells disclosed).

No BREAK: nothing dishonest or unsound; the corrections are wording/comment
freshness, not a design defect.
