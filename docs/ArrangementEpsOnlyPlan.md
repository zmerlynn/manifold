# Plan: tolerance out of the arrangement (the #1757 posture)

Directive: assume upstream #1757 ("Make the intersection-vertex merge
eps-only") is more correct, and refactor RSI to match. #1757's two 2D
findings, both empirically demolished there: tolerance-scale
new-to-old snapping lets a PROPAGATED bound (unbounded, input-
dependent) move arrangement geometry and fold it; 10x-eps new-to-new
fusion welds genuinely distinct crossings into a non-manifold pinch.
Its fix: the arrangement is eps-only (new-to-new at eps, new-to-old
at 2 eps, clusters with an old vert pin to the old vert), and
tolerance-scale decimation is Simplify's job.

RSI already half-holds this posture. Steps 6/7 snap at bare eps with
comments naming the exact hazard ("a tolerance-scale radius on a
tolerance-inflated input would drag pierce events onto far verts and
deform the arrangement"); every new-to-old resolution pins to the
existing id (originals are never centroid-dragged - #1757's anchor
rule is already RSI's); conditioned radii widen only where the
conditioning is computable at the source. The violations are four
allocation-time snap sites using `tolerance + eps`, three
triangulation-eps seeds mixing `tolerance_` in, the rebuild
construction sweep consuming the seeded claim (the latter two
classes found by review rounds 1 and 2; the draft claimed four
sites total), one 10x-eps fusion radius, and the docs/comments
citing old-boolean2 semantics (`newToOldThresh`,
`kIntersectionMergeEpsFactor`) as precedent - semantics #1757
abolishes, leaving the citations dangling.

## Complete tolerance inventory (review round 1 demanded it; verified)

Every read of `tolerance`/`tolerance_` in the pipeline, classified:

- ARRANGEMENT-AFFECTING, Phase 1 removes: the four allocation snaps
  (6.5 corner-snap base :1552; step-9 contacts :1904; step-9 cluster
  resolution + threading guards :2020; step-9 propagation guard
  :2169) AND three triangulation-eps seeds review round 1 caught the
  draft missing - the emit retry-ladder seed (:622), the step-10
  island TriangulateIdx seed (:2629), and the seed-cast target
  triangulation (:3644/:3649). The triangulator's eps controls
  short-edge clipping, ear validity, and hole admission, so a
  tolerance-scale seed lets a propagated claim reshape emitted
  triangulations and seed-target choices - same posture violation,
  different mechanism. The adjacent comments justify only the eps
  term (pipeline jitter); the tolerance_ term was unmotivated. These
  become `max(impl.epsilon_, eps)`.
- CONSTRUCTION-SWEEP, Phase 1 reorders (review round 2 caught the
  draft classifying this as claim-only): BuildImplFromTris seeds
  `tolerance_` BEFORE its construction sweep, and core's
  SetNormalsAndCoplanar / RemoveDegenerates CONSUME tolerance_ - so
  both RSI rebuilds (the MergeVertsEps rebuild at :857 and the final
  emit at :667) reshape geometry at claim scale. The fix: the sweep
  runs with the epsilon-floored default (tolerance untouched until
  after RemoveUnreferencedVerts), then
  `tolerance_ = max(tolerance_, toleranceSeed)` - the claim is
  still exported, epsilon floor preserved, but it no longer licenses
  construction-time geometry changes. This deliberately diverges
  from the Manifold(MeshGL) ctor chain ON THE SWEEP-SCALE ONLY, with
  the #1757 rationale: a claim DESCRIBES movement already applied,
  it must not LICENSE more. PROBE D (recorded): the reorder is
  fully green, 553/553.
- OUTPUT-CLAIM, kept: the output tolerance floor (:669) and the
  toleranceSeed values both callers pass (the input's claim; the
  emit formula). The driver's `tolerance` local (:451) survives only
  for :669, moved next to it.

## Probe results (recorded before this plan was written)

- PROBE A - `tolerance + eps` -> `2 * eps` at all four sites (step 6.5
  corner-snap base; step-9 contacts snap; step-9 propagation t-guard;
  step-9 cluster resolution + threading t-guard): suite 551/553.
  Both failures are white-box pins whose FIXTURES use tolerance as
  the amplifier (the t-recompute reorder pin sets tolerance 0.03 to
  invert a t-sort; the 6.5 snapped-endpoint-subdivides pin sets
  tolerance 0.05 to reach a far apex). The rules they pin are
  radius-independent. Every feature-level fixture green.
- PROBE B - additionally `mergeR` 10x -> 1x eps: suite 547/553. The
  four additional reds are step-9 merge-mechanics unit pins whose
  fixtures use 10x-scale member spacing, including the collapse pin
  that documents itself as breaking below 9x. Feature level still
  fully green.
- PROBE C (run after review round 1 found the triangulation seeds) -
  `tolerance_` dropped from all three seeds: suite 553/553, fully
  green. Nothing rides on tolerance-scale triangulation seeding.
- Conclusion the probes license: nothing at feature level rides on
  tolerance-scale arrangement snaps or the 10x fusion. The recovery
  scaffold (gates, cap walker, reducers) guards FP-degeneracy
  classes, not self-inflicted tolerance damage - no gate becomes
  removable on this evidence, and the plan claims no such
  simplification.

## Design

Phase 1 (the core - tolerance leaves the arrangement):

- Principle, per site (review round 1 demanded the candidate-class
  analysis): the base budget is `2 * eps` (#1757's: the combined eps
  reach of two represented points), and the CONDITIONED radius - the
  computable, source-gated uncertainty the old `tolerance + eps`
  covered only by accident when tolerance happened to be large - is
  the sole legitimate widener. Site by site:
  - 6.5 corner snap (:1552): candidates are the pair's six ORIGINAL
    corners; the crossing's own conditioning is in hand. Stays
    `max(2 * eps, condOnly)` as drafted.
  - Step-9 contacts (:1904): the tested endpoint is NEW-capable (6.5
    allocations are chord endpoints) and carries a conditioned
    radius in `trace.newVertSnapR`. The pass gains a `snapR` view
    (net signature size unchanged - it loses `tolerance`); per
    contact the radius is `max(2 * eps, snapR(endpoint))`.
  - Step-9 cluster resolution (:2020): candidates (endpoints,
    extras, contacts) are NEW-capable, `newVertSnapR` is already in
    hand (the function resizes it). The radius splits by candidate
    class (review round 2 caught the draft's max-everywhere rule
    re-introducing a wide new-to-new weld): ORIGINAL candidates
    resolve at `max(2 * eps, cl.snapR)` - original positions are
    exact, the cluster's uncertainty disc is the only spread, the
    9.5 original-snap convention. NEW candidates resolve at
    `max(2 * eps, min(cl.snapR, snapR(candidate)))` - the 6.5
    new-to-new SOURCE-GATING convention (:1568-1578): both claims
    must be conditioned before the radius widens, else an
    ill-conditioned cluster absorbs a well-conditioned DISTINCT new
    vert - the exact weld class 6.5's min-gate exists to prevent.
  - Step-9 propagation guard (:2169): endpoint-zone exclusion for
    the cluster being propagated; its aggregated `cl.snapR` is
    in-function. `max(2 * eps, cl.snapR)`.
  - Threading guard + t-recompute (:2103): applies to newly added
    records whose per-record radii are in `newVertSnapR`; per record
    `max(2 * eps, snapR(record))`, mirroring the contacts pass it
    backstops. (Reviewer may judge flat `2 * eps` sufficient here -
    the failure mode is a fail-safe rim - but the per-record form
    costs nothing: the vector is already indexed in the loop.)
- The three triangulation-eps seeds drop `tolerance_`:
  `max(impl.epsilon_, eps)` (probe C: fully green).
- Signature simplification, mirroring #1757's own: the `tolerance`
  parameter drops from CoplanarTraceChords,
  FindOnChordEndpointContacts, MergeAndPropagateCrossings,
  ResolveAndThreadClusters, and ResolveAndThreadCrossings
  (FindOnChordEndpointContacts swaps it for the `snapR` view). The
  driver's `tolerance = max(impl.tolerance_, eps)` local survives
  ONLY in the output-tolerance claim (the output floor must still
  honor the input's claim) - moved next to that use, commented as
  output-claim-only.
- Step 9.5 is OUT OF SCOPE and unchanged, stated with its reason:
  its new-new unite radius is frame-to-frame FP spread of ONE
  computed point (not distinct-crossing fusion), its new-onto-
  original band is load-bearing for step 6's contract ("events in
  the (eps, 10 eps] band are unified onto originals by step 9.5
  anyway"), and both radii are bounded pipeline-error scale - the
  #1757 indictment is specifically the UNBOUNDED propagated
  tolerance. Touching 9.5 re-arms verification of the whole
  unification pass for no posture gain.
- Steps 1/6/7, the kernel radii, the conditioned-radius machinery,
  and the output formula are already conformant and untouched.

Phase 2 (the fusion radius - same pass, separable decision):

- `mergeR` in MergeAndPropagateCrossings: 10x eps -> eps. True
  k-fold concurrences still cluster (their raw crossings land
  sub-eps apart).
- What Phase 2 actually changes, stated precisely (review round 1
  resolved the draft's convergence question AGAINST the draft's
  distinctness claim): step 9.5's new-new unite is flat 10x eps -
  conditioned radii widen only its new-onto-ORIGINAL snap - so two
  crossings 9 eps apart that step 9 now keeps distinct are STILL
  united by 9.5. Final identity-distinctness is NOT a Phase-2
  outcome and the plan does not claim it. The real differences:
  step-9 fusion emits a manufactured CENTROID re-projected onto the
  host face (the collapse pin documents the off-chord cost) and
  eagerly propagates it onto every chord of the involved faces at
  10x scale; 9.5 unification picks a REAL member vert (an original,
  else the smallest new id) and only remaps existing consumers. So
  the eps fusion means: no manufactured positions, no 10x-scale
  eager propagation to chords that did not already consume a member
  - a 9-eps-spread pair was never a true concurrence and no longer
  gets a concurrence's treatment - while identity still fuses
  downstream at a real member's position. Stated as the full trade
  (round 3): off-chord consumption in the 9.5 band is NOT
  eliminated - after the remap, the chord whose own crossing lost
  the representative election consumes a point off itself by up to
  the FULL pair spacing, versus the old centroid's half-spacing off
  BOTH chords. The win is categorical (a real computed point, no
  manufactured one; no eager propagation), not metric. That is
  #1757's "don't let the merge manufacture geometry" posture,
  scoped to what step 9 controls.
- If review judges that benefit too thin to carry six test
  re-fixtures, dropping Phase 2 leaves Phase 1 intact - they share
  no code. The 10x mergeR then keeps an OWN-rationale comment (an
  error-budget fusion trade this pipeline accepts, measured by its
  pin) instead of the dead boolean2 citation.

Phase 0 (separable commit, recommended): port #1757 verbatim into
this branch's split-layout boolean2 (`src/cross_section/boolean2/`:
driver.cpp thresholds + anchor rule + the constant's deletion from
intersections.h, the tolerance parameter removal from driver.h /
boolean2.h / boolean2.cpp, and the cross_section.cpp Boolean/
BatchBoolean callers). The branch's 2D backend currently carries both
pre-#1757 annihilation bugs - live exposure for CrossSection users of
this branch (the demo agent). Clearly labeled as a port that drops
out when the stack rebases past #1757's upstream landing. The
ApplyFillRule(contours, tolerance_) constructor call is outside
#1757's diff and stays as-is.

## Simplifications gained (and not)

- Five RSI signatures drop a parameter; the threading prose about
  snaps moving verts "up to tolerance + eps" tightens to the 2-eps /
  conditioned scale; the 6.5 "kept SEPARATE from the snap base"
  guard comment dissolves into one line (the conditioned radius is
  the only widener left).
- The docs' eps-policy taxonomy loses its `tolerance + eps` category
  entirely: eps (identity), 2 eps (allocation new-to-old),
  conditioned (computable widening), 10x (9.5 frame-spread
  unification). One category fewer to reason about at review time.
- Tolerance-inflated inputs - including RSI's OWN resolved outputs,
  whose claims run 10x-1000x base eps after the retry ladder - stop
  having their re-pass arrangements built with claim-scale snap
  radii. A second RemoveSelfIntersections pass now arranges at the
  same radii as the first. Previously an undocumented hazard; the
  idempotence posture strengthens for free.
- NOT gained, stated honestly: no gate, reducer, or recovery pass
  becomes removable - the probes show the suite identical at feature
  level, so that machinery was never compensating for tolerance
  damage on covered inputs. The open question whether the 10x fusion
  manufactures members of the PINCHED class (stack-rank #4) on real
  geometry needs a feature-level pinch fixture, which does not exist
  yet; that probe belongs to the #4 arc, run under BOTH radii
  postures.

## Tests

- TDD anchor (red first, committed DISABLED_ per the house seed-queue
  convention): arrangement tolerance-independence at feature level -
  the same piercing geometry processed with a low and a high seeded
  `tolerance_` must produce identical RemoveSelfIntersections
  geometry. Red today (the high seed inflates the snap radii, the
  triangulation seeds, AND the rebuild construction sweeps); green
  after Phase 1, which is defensible only because the inventory
  above is complete. Seeding path review-verified: MeshGL64's
  tolerance field reaches Impl::tolerance_, survives SetEpsilon's
  floor and the MergeVertsEps rebuild, and lands in the driver's
  max(impl.tolerance_, eps). FIXTURE DISCIPLINE (review round 2):
  core's OWN ctor chain also consumes tolerance, upstream of RSI -
  the anchor asserts as its PREMISE that the two
  differently-seeded inputs are geometry-identical AFTER
  construction, pre-RSI (a fixture whose ctor sweeps are
  tolerance-neutral no-ops), so the pin isolates RSI's sensitivity
  and self-reports a bad fixture instead of misattributing a ctor
  difference. COMPARATOR (round 3): the two runs legitimately
  differ in EXPORTED tolerance - the output claim is seeded by the
  input's - so the anchor compares topology and positions with a
  comparator that EXCLUDES MeshGL64::tolerance
  (ExpectMeshGL64GeometryIdentical checks the tolerance field and
  would spuriously fail) and asserts the two output claims
  separately. This is THE pin of the plan's posture, the analogue
  of #1757's "tolerance no longer affects the arrangement at all".
- Re-fixture, not delete, the two probe-A reds: the t-recompute
  reorder pin re-amplifies via a short host chord (a 2-eps-scale snap
  on a short chord still inverts t-order - the rule stays
  load-bearing); the 6.5 subdivides pin re-fixtures with the apex
  inside the CONDITIONED radius (the widener that survives) instead
  of inside tolerance.
- Phase 2 flips: the collapse pin INVERTS - two distinct crossings
  9 eps apart now stay distinct at step 9 (the direct port of
  #1757's lesson), with a comment noting 9.5's separate
  unification role; the face-gate, cluster-max-radius, and
  reprojection pins re-fixture at sub-eps member spacing (every
  mechanism they pin survives - only the spacing amplifier changes).
- Phase 0 carries no new tests here: the port's behavior pins live
  upstream with #1757 (and its test PR); this branch's CrossSection
  suite must stay green over the port, which the upstream PR already
  validated on master's layout.
- MECHANICAL SWEEP, compile-surfaced (rounds 1 and 3): every direct
  seam caller of the five signatures losing the tolerance parameter
  updates mechanically. The caller GROUPS: the step-9 white-box
  block (contacts / threading / merge / conditioning tests) and the
  step-6.5 white-box block (trace / corner-snap / subdivision
  tests) - the compiler surfaces the full list, the plan does not
  claim to enumerate it - and their tolerance+eps PROSE updates
  with them (a comment describing a radius the call no longer takes
  is a lie).
- All RemoveSelfIntersections feature fixtures, the seam tests, and
  the eps-retry ladder pins must stay green untouched throughout.

## Docs

OverlapRemoval.md: the eps-policy section rewrites - the
`tolerance + eps` bullet is deleted, the 2-eps allocation budget
documented with #1757 cited as the posture source instead of
old-boolean2 as precedent; the 10x bullet narrows to 9.5's
frame-spread rationale (its "matching kIntersectionMergeEpsFactor"
citation goes); the step-9 sections update. Code comments at the
four snap sites + the merge + the three triangulation seeds +
internal.h's MergeAndPropagateCrossings block rewrite the same way,
and the step-7 comment's "NOT the tolerance + eps allocation radius"
contrast (:1040) updates - the radius it contrasts against no longer
exists. docs/InteriorIslandPlan.md's two descriptions of the
TriangulateIdx/emit seed formula (:65-66, :114-116) update with the
seeds they describe. The known-limitations section is re-checked for
tolerance-snap mentions. manifold.h needs nothing (no public
contract changes - output tolerance semantics are unchanged).

## Out of scope / stop conditions

No step 9.5 changes, no gate changes, no new conditioned-radius
MACHINERY (Phase 1 only plumbs the existing per-vert radii to sites
that should have used them), no kernel/predicate changes, no public
API changes. STOP if implementing Phase 1 flips any feature-level
fixture beyond the recorded unit re-fixtures (the probes say it
will not - a divergence is new information), or if the
tolerance-independence anchor stays red after Phase 1 (that would
mean a NINTH tolerance consumer; the inventory above - four snaps,
three triangulation seeds, one construction sweep - is claimed
complete after rounds 1 and 2 each found consumers the draft
missed), or if Phase 0 flips any CrossSection test (upstream
validated green; a flip here means the split-layout port diverged
from the diff).
