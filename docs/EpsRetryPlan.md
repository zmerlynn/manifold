# Plan: eps-retry ladder for gate fallbacks (deliberately small)

Stack-rank #3: tangent-degenerate inputs that gate at the inferred
eps but resolve at an inflated one - proven empirically before this
plan was written (TDD anchor
`DISABLED_RemoveSelfIntersectionsHullResolvesAtOriginScale`,
committed red): the origin-scale hull falls back at 1x-40x the
inferred eps and RESOLVES at 100x (strict pierce reduction, all three
hulls kept, volume preserved) - the same outcome the identical
geometry already gets at 1e4 scale, where the bbox-derived eps is
naturally that much larger.

The simplicity orders from the last pass apply: no new seam
parameters, no reason enums, no per-gate plumbing. The whole feature
is a short loop in ONE function.

## Design

In `RemoveOverlaps` (the wrapper - not the pipeline body, not the
member):

- Applies ONLY when the caller passed `eps <= 0` (the inferred-eps
  path, i.e. the public member). Explicit-eps callers keep exact
  single-attempt semantics - every existing seam test's behavior is
  frozen by construction.
- First attempt passes the caller's eps through UNCHANGED (the impl
  infers the base for eps <= 0 internally, as today). If it returns
  a value (success OR a Cancelled-status Impl), return it -
  cancellation never retries.
- On nullopt: compute `CheckSelfIntersection(input)` ONCE, here on
  the rare fallback path. Zero STRICT-INTERIOR pierces makes the
  nullopt retry-ineligible - return it. Not "clean": island/hazard
  gates can fail closed without strict pierces, and those nullopts
  are equally non-retried. Safety and cost agree here: a retry
  acceptance must strictly reduce a positive count, so a zero-pierce
  retry could never be accepted anyway - the guard skips attempts
  that cannot change the outcome, and wider-eps welding of
  pierce-free geometry is impossible under the same rule.
- Retry-safety framing, stated precisely: a retry is the FULL
  pipeline at a coarser eps with every gate and the monotonicity
  check active - the same posture as the far-from-origin success,
  where the coarser eps arises from scale instead of a factor. No
  gate is bypassed; a class that gates for correctness (island,
  hazard, fold) either gates again or genuinely conforms at the
  coarser eps. The claim is NOT "all piercing nullopts resolve
  safely" - it is "retrying cannot produce anything a direct coarser
  run could not". IMPLEMENTATION PROBE: run the known piercing
  fallback fixtures (the ovoid; the pinched-class geometry at
  feature level if constructible; the hull) under 10x/100x and
  record per fixture whether it resolves or re-gates - surprises are
  stop conditions, not footnotes.
- Otherwise retry at 10x, then 100x the inferred eps (the probed
  ladder: the known class needs 100x; 10x first so inputs that
  resolve earlier keep the tighter tolerance). A RETRY-rung success
  is accepted only if it STRICTLY reduces the interior pierce count
  (one CheckSelfIntersection on the candidate, rare path) - the base
  attempt keeps the standing equal-pierce acceptance, but a wider-eps
  rebuild that fails to reduce pierces is churn, not progress, and
  review simulation showed accepting it breaks idempotence (an
  equal-pierce second pass with wider tolerance, then component
  drift). A Cancelled-status Impl from a retry rung is returned
  immediately, exempt from the acceptance check - the valued-return
  rule applies to every attempt, not just the first. Rejected
  candidates are treated as nullopt and the ladder continues. All
  attempts nullopt -> nullopt (the bit-identical fallback contract
  is the member's, unchanged).
- Mechanics (signature-free; revised in round 5 - the hoist is OUT):
  the draft hoisted InferEps into the wrapper, and review showed the
  hoist must then replicate the impl's exact entry order (cancel
  poll BEFORE the empty-input nullopt, per CancelBeforeRunEmptyInput,
  and the empty check before InferEps - an empty Impl's default Box
  is infinities, which AlphaBudgetEpsilon feeds to frexp). Rather
  than duplicate that sequence in two places, the wrapper does not
  hoist at all: the first attempt is the caller's eps verbatim, so
  the impl's live entry order and its eps<=0 inference stay the
  single copy. Retry rungs launch only after "nullopt + input
  pierces > 0", which implies a non-empty input, and the wrapper
  computes base = InferEps(input) THERE - the same deterministic
  function the impl just used, recomputed on the rare path rather
  than threaded out (no seam surface). One comment at that call
  site pins the coupling: the rungs scale the inference the first
  attempt used. The seam header's eps<=0 doc updates to describe
  the ladder semantics.
- Cancellation between attempts, precedence-preserving: the wrapper
  polls IsCancelled ONLY immediately before launching a retry rung -
  new work - never after a clean-class nullopt (pierces == 0) and
  never after the final rung's nullopt. Completed work keeps winning
  over a racing cancel exactly as the empty-chord arm documents; a
  cancel observed where new work would begin returns the
  Cancelled-status Impl.
- Debug-throw handling per attempt: the current MANIFOLD_DEBUG
  try/catch wraps the single impl call; the ladder factors it into a
  per-attempt helper so a debug-only assertion throw is treated as
  THAT RUNG's nullopt and the ladder continues, instead of aborting
  the whole ladder.
- Roughly twenty lines plus a `kEpsRetryFactors` constant.

## Cost, stated honestly

A clean input takes zero extra ATTEMPTS but one extra pierce check:
the impl computes the input pierce count internally and the wrapper
cannot see it (no reason channel by design), so on a nullopt it
recomputes once to make the retry decision. Threading the count out
would need seam surface; rejected for simplicity. Successful first
attempts pay nothing extra.

## What the contract already handles

- Tolerance honesty is automatic: the output tolerance formula floors
  at 10x the WORKING eps of the attempt that succeeded, so a 100x
  retry claims a 1000x-base floor. One doc sentence states that
  retries widen the claim.
- Pierce-monotonicity, fail-closed gates, ctx polls: each attempt is
  a full ordinary pipeline run; nothing inside changes.
- Idempotence: the strict-reduction rule restores it for the probed
  class (review re-simulation: second pass base nullopt, the 100x
  candidate's 7->7 rejected, input returned unchanged). The claim is
  scoped to that class plus one stated assumption: a BASE-attempt
  success may still accept equal pierces (the standing single-attempt
  contract, unchanged by this plan) - the general idempotence pin
  remains the existing test, flipped to two-pass semantics.
- Determinism: a fixed two-entry ladder.

## Tests

- The TDD anchor un-DISABLEs and becomes THE surviving public hull
  pin; `RemoveSelfIntersectionsHullMaskFixture` (the bit-identical
  fallback pin over the same fixture) is DELETED, its comment history
  folded into the anchor's - one fixture, one pin, no duplicate.
- Clean-input-never-retried, stated honestly: the existing
  CleanInputUnchanged / BooleanResult bit-identical passthroughs
  must stay green untouched, but they pin the OUTPUT, not the
  absence of retries - an implementation that uselessly retried a
  clean fixture and still fell back would pass them. Never-retried
  is a code-shape invariant (pierces == 0 returns before any rung),
  carried by review like the inferred-path first attempt below; no
  instrumentation pin is added for it. The eps-chain weld test's LE
  bound discriminates over-inflation on the EXPLICIT-eps path only
  (it calls the seam directly); an inflated FIRST attempt on the
  inferred path has no direct pin - the ladder's code shape (base
  attempt first, factors applied only on retry) plus review carries
  it, stated here rather than overclaimed.
- The ovoid dense-sliver fixture: run it under the ladder during
  implementation; keep its fallback pin if it still falls back at
  100x, flip it with justification if it resolves. Either outcome is
  recorded, not assumed (the ovoid does NOT inherit fold coverage -
  which gate fires for it is unproven and stays that way).
  Feature-level fold coverage, stated honestly: there is no cheap
  beyond-ladder fold fixture - uniform scaling moves the clusters
  and the eps base together, so the resolving factor is
  scale-invariant (review caught the draft's 0.01x arm as exactly
  this error; FarFromOrigin works only because TRANSLATION grows the
  bbox without growing the clusters). The public fold-gate-identity
  pin is therefore LOST with the hull flip and the fold gate's
  discrimination narrows to its unit arms - accepted and recorded.
  The public FALLBACK-OUTCOME pin (some piercing input falls back
  bit-identically through the member under the full ladder) is
  carried by whichever existing fixture still gates when probed at
  implementation time (the ovoid is the candidate); if none does,
  that is reported, not papered over.
- Idempotent and Deterministic currently use the origin hull as a
  bit-identical FALLBACK fixed point; under the ladder it becomes a
  RESOLUTION. Their origin-hull arms flip to success-path semantics:
  geometry-only determinism (the far-scale arms' existing pattern)
  and the explicit two-pass expectation (first pass resolves, second
  pass returns its input unchanged under strict reduction).
- Explicit-eps single-attempt freezing needs no new pin: every seam
  test IS that pin.

## Docs

manifold.h's RemoveSelfIntersections comment gains one sentence (the
inferred-eps path retries at a small fixed ladder of wider epsilons
before falling back; resolved outputs report correspondingly wider
tolerance). OverlapRemoval.md: the driver section describes the
ladder; Known limitations item 1 narrows (the hull class now resolves
at origin scale via retry; the residual is the conditioned-band
pierce remainder, already documented); AND the validation coverage
map's hull line (which today says the fixture pins folded-shell
bit-identical fallback) plus the known-limitation fallback and
idempotence wording all update to the resolve semantics - the full
stale-prose sweep, not just the two named sections. The sweep
includes TEST comments: the clean-passthrough, Idempotent,
Deterministic, and FarFromOrigin blocks (the last describes far
scale as the success-path complement of origin-scale fallback)
still describe origin-hull folded-shell fallback semantics in
prose and must flip with their assertions; the BooleanResult
comment's "HullMask and SelfIntersect below cover the
actually-piercing paths" pointer retargets to the hull resolution
anchor plus the ovoid fallback pin once HullMaskFixture is deleted.

## Out of scope / stop conditions

No per-gate retry policies, no adaptive factor search, no exposing
the ladder publicly, no changes inside RemoveOverlapsImpl at all
(its entry order and eps<=0 inference stay live and untouched).
STOP if the ladder interacts with cancellation or tolerance in any
way the existing formula does not already express. A draft stop
condition about losing the fold gate's last feature-level fallback
pin is RESOLVED, not pending: flipping the hull pin does lose it
(the bent-open-fold arms are unit-level, not feature-level), and
that loss is the accepted, recorded outcome of the Tests section -
the ovoid survives only as probed public fallback-OUTCOME coverage,
not fold coverage.
