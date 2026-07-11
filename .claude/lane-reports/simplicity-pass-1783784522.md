# Simplicity pass - overlap3 sweep-native emission (v3)

Branch: explore/sweep-plane-3d-v3 @ 83a80a48
Role: GENERATOR for adversarial simplicity crucible. Sole editor.
Bar (owner): "simplicity trumps... if there's a lot of special case garbage we'll rewrite it."

## Plan

1. Read authorities: SKILL.md (house style), SweepEmit3D.md (spec), git log
   bd2a330e..83a80a48 (accretion window, seven arcs). DONE.
2. Sweep every in-scope file end to end, cataloguing candidates with file:line:
   - src/overlap3.h, src/overlap3.cpp, src/overlap3_sweep.cpp
   - branch-added portions of boolean2.{h,cpp}, boolean2_sweep.cpp (diff vs 8bffb521)
   - test/overlap3_test.cpp (structure/comments only)
3. Classify each candidate: behavior-preserving (implement) vs behavior-changing
   (record as finding for the judge). Priority = review friction first.
4. Implement the safe ones in logical commits; format + build + gtest_filter each.
5. Full suite once before final commit. Commit notebook with final commit.

Behavior preservation is ABSOLUTE: 601 total green, prefer bitwise-identical,
no test assertion/fixture changes.

## Known suspects (from prompt; verify don't assume)
- ComputeCap signature + optional<pair<FatalReason,string>> return repeated
  across EmitCaps/EmitStrips/ComputeCap - result struct or alias?
- Orphans from deleted mechanisms (ChainSplitVerts died in provenance arc);
  dangling helpers/includes/comments re dead machinery.
- SlabResolver candidate cascade vs provenance channel - duplicated concepts.
- Corpus test helpers CorpusPairGate + perf-arc singles - near-duplicate merge.
- Comment staleness across seven arcs' seams.
- capEps=8, kAngleTie, adjudicated constants - each carries justification comment?

## Catalogue (file:line -> candidate)

### IMPLEMENT (clearly behavior-preserving, review-friction first)

A. overlap3.cpp EmitStrips/ComputeCap/EmitCaps all return
   `std::optional<std::pair<FatalReason,std::string>>` (652,783,869); SweepEmit
   reads `->first`/`->second` (1184-5,1190-1). Name it: `struct Fatal{reason,
   detail}` + `using MaybeFatal = std::optional<Fatal>`. Suspect #1 confirmed.
   Style skill prefers small structs over pairs; kills .first/.second friction.
B. overlap3.cpp:842-844 ComputeCap comment trails history: "This replaces the
   geometric on-chord projection (ChainSplitVerts), which diverged... corpus-
   falsified." ChainSplitVerts is deleted; docs own the history. Trim to the
   standing contract. Comment rule: delete vestigial narrative naming dead code.
C. overlap3.cpp RemoveOverlaps3D (1215) and RemoveOverlaps3D_TestHooks (1238)
   duplicate a ~9-line prefix: eps resolve+validate, Canonicalize, empty-faces
   early-out, FindSeams, seams.fatal check. Fold into one helper. BONUS: revives
   the DEAD SeamsResult.fatal/detail (never written by FindSeams - confirmed via
   grep; the two `seams.fatal.has_value()` checks are currently dead) by routing
   SubEpsInput through it. Traced empty-faces equiv: FindSeams copies
   canon.faces->arr.faces, so arr.faces.empty() iff canon.faces.empty().
D. test/overlap3_test.cpp:389 "sentinels-only slabs" is STALE (sentinel slabs
   removed; BuildSlabs comment at overlap3_sweep.cpp:111 says "No sentinel
   slabs"). Trim to mechanism-neutral wording. In-scope (comment, not assertion).

### RECORD as findings (behavior-changing or judgment calls; NOT implemented)

E. overlap3.h:153 MergedVert is a single-field `{vec3 pos;}` wrapper; every use
   is `.pos`. Flatten to std::vector<vec3>? HIGH churn across .cpp + test hooks'
   `.pos` accesses; style skill discourages rename-only churn and it is in the
   test-hook surface. Finding, not implemented.
F. Four fatal-carrying idioms coexist: StageResult<T>, SeamsResult,
   optional<pair> (fixed by A), Overlap3Result. SeamsResult could be
   StageResult<ArrangementGeometry>. Deferred: FindSeams always populates arr
   even on (never-fired) fatal; unifying is invasive. C partially addresses by
   giving SeamsResult.fatal a live use.
G. boolean2.cpp:483 AssembleChain sits in a branch-added REOPENED anonymous
   namespace. File already reopens anon ns 5x (upstream, 255/545/675); style
   crucible explicitly DEFERRED "upstream boolean2.cpp's reopened anonymous
   namespaces" as upstream follow-up. Local convention is the tiebreaker -> leave.
H. test CorpusPairGate vs CorpusSingleGate (1531,1584): share OBJ-open
   boilerplate but differ semantically (pair has a+b oracle + NonManifoldEmission-
   only guard; single has NO oracle + any-guard). Already factored as two
   helpers; the residual dup is 4 lines of model-open (also in Gate4c). Low
   value, test-only. Leave (merging would parameterize away real differences).
I. ComputeCap cap_plus/cap_minus blocks (819-836) near-duplicate
   OutEdgesToPolygons+TriangulateCap. A 2-instance emit-lambda is possible and
   bitwise-safe (fatal path discards `out`), but marginal. Deferred unless time.

### Examined, NO ACTION (spec-demanded / already deduped)
- [R2-fold]/[R3-fold] etc. code comments = spec anchors to SweepEmit3D.md
  adjudications; prompt says spec-demanded arms are NOT garbage. Keep.
- capEps=8*eps (811) and kAngleTie=1e-9 (977) each carry justification. Keep.
- CollectGroupedMemberEdges already hoisted (perf arc); two call sites are two
  stages (FindSeams skeleton, SweepEmit resolvers) that cannot share via arr.
- SlabResolver cascade vs provenance channel: different concerns (resolver =
  endpoint track; edgeSubdiv = interior subdivision). No dup.
- CollectThenMeasure 7-param signature: each param spec-demanded (M4 negOut,
  provenance classSubdiv, capture, conflicts, cleanOut), all null-defaulted,
  2D-caller-unchanged. Keep.
- corpus-falsified/diagnosed parentheticals at 534,807,886 justify numeric
  constants/rules (why-breaks-if-removed) - style skill protects these. Keep.

## Log

- Sweep complete. Baseline green: filter Overlap3.*:Boolean2* = 89 pass + 1
  Gate4c skip (corpus singles slow ~37s each -> exclude during iteration).
- Plan: implement A, B, C, D as logical commits; full suite before final.

### Implemented (each built + filtered-green before commit)

- A+B (1c28d74f): struct Fatal{reason,detail} + MaybeFatal alias replacing
  optional<pair<FatalReason,string>> in EmitCaps/ComputeCap/EmitStrips;
  ->first/->second at the SweepEmit sites become ->reason/->detail. Trimmed the
  ChainSplitVerts history tail in the ComputeCap provenance comment.
- C (39546c76): PrepareArrangement hoists the 9-line canon+seams prefix out of
  both entry points; revives the previously-dead SeamsResult.fatal by routing
  SubEpsInput through it. Traced empty-faces equivalence and counters threading
  before/after - identical.
- I (8554dfc2): ComputeCap's cap_plus/cap_minus emits folded into one
  emitCapSide lambda. Bitwise-identical (same walks/triangulations/order;
  verified the fatal-path `out` discard makes the entry-guard reordering moot).
- D (ef56d8dc): dropped stale "sentinels-only slabs" test comment.

- Did NOT implement E/F/G/H (recorded above): MergedVert flatten (churn into
  test hooks), full fatal-idiom unification (invasive), AssembleChain namespace
  (deferred upstream per style-crucible close), corpus helper merge (real
  semantic asymmetry). Left for the judge round.

- Filtered gate (Overlap3.*:Boolean2*, slow corpus singles excluded): 85 pass +
  1 Gate4c skip after every commit. FULL SUITE (601): 600 pass + 1 Gate4c skip -
  identical to baseline. Behavior preserved.

### Close
Five commits, all overlap3.cpp/test readability + dedup, zero behavior change.
No "special case garbage" removed because the accreted special cases (SubEps
guards, canonical-cap rule, capEps floor, coplanar mechanisms, sheet splitter)
are all spec-demanded and carry justifications - flagged none as garbage. The
review-friction wins were: the pair->struct fatal return, the duplicated
entry-point prefix (with a dead field revived), and the two-cap-side dedup.
Findings E/F/G/H/(SeamsResult unification) handed to the judge round.
