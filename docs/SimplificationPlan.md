# Plan: simplification pass over the overlap-removal feature

The interior-island arc landed at roughly six times the weight of its
core idea, almost all of it validation accreted across nine
adversarial review rounds - each round ADDED a defense, and nothing
asked whether the total was proportionate to house norms. This pass
removes the disproportionate weight across the WHOLE feature while
preserving every behavioral contract except TWO named, argued
changes: the PSLG output-validation deletion (a user-approved trust
decision) and the edge-count cap deletion that rides with it (the cap
bounded only the deleted quadratic checks; large valid holed regions
move from fail-closed to decomposed).

## The proportionality argument (what the deletion actually is)

Two claims, kept separate (the first draft conflated them):
- TRUST PRECEDENT (true): manifold consumes `TriangulateIdx` output
  unvalidated everywhere else - the RSI emit in this same file
  (`Triangulate`, which wraps `TriangulateIdx`)
  (which runs AFTER the gates and is not backstopped by them either;
  only empty/status/volume/pierce checks follow it) and Extrude,
  which feeds the output straight into CreateHalfedges with no
  fail-closed wrapper at all (CrossSection-backed extrusion is the
  softer comparison: its contours are Boolean2-regularized first; raw
  Extrude(Polygons) is the direct precedent). The island call site's
  INPUT is proven before the call (clean simple loops, strict
  containment, sub-resolution holes gated) - guarantees at least as
  strong as any other call site's.
- NETS DO NOT COVER THIS (also true): a PSLG-violating triangulation
  that passes the kept triad would NOT be caught downstream -
  BuildCellComplex and the emit-topology check key by exact endpoint
  pairs (geometric crossings are invisible to them) and the final
  pierce gate returns zero for coplanar contacts, which in-face
  overlaps are.

So the deletion is a TRUST decision, stated plainly and APPROVED by
the user on exactly that framing: on proven-valid input we extend
`TriangulateIdx` the same unchecked trust the rest of the library
extends it on far weaker input guarantees. The residual
risk is a triangulator defect on valid input - a class the codebase
accepts universally. The kept triad still catches the REALISTIC
failure modes verified in review: a misclassified (non-hole) contour
fails positive-area or boundary-direction coverage; an empty result
gates; a filled 3-vert hole is rejected by the hole-key check.

## Cuts

1. ISLAND OUTPUT-VALIDATION DIET (the bulk). Keep, verbatim: the
   per-triangle checks (three distinct ids, strictly positive area,
   not canonical-equal to a hole contour), the undirected boundary
   coverage (boundary edges exactly once in the required direction,
   no reverse occurrences; interior edges paired opposite), and the
   signed area-preservation check - this triad is cheap and is
   exactly what catches a hole the triangulator misclassified.
   DELETE: the full PSLG suite over the triangulation output (strict
   crossings, vertex-on-nonincident-edge, collinear overlap,
   duplicate-edge analysis, and all of its adjacent-pair machinery)
   and the edge-count cap. The cap deletion is a SECOND, named
   behavior change: today a large valid holed region fail-closes at
   the cap; afterwards it decomposes. The cap existed solely to bound
   the quadratic PSLG checks - the surviving triad is linear, so no
   replacement cap is warranted, and the success class deliberately
   widens for large holed regions (gates beneath unchanged). The INPUT-side
   checks are untouched: `isLoopSimple2D` (with its adjacent-pair
   handling - it proves OUR walk's input, not the triangulator's
   output), strict PIP with border rejection, the sub-resolution
   fast-fail, the hazard flag.
2. SHARED 2D PREDICATES - exactly TWO live duplicates: pointOnSeg
   (isLoopSimple2D) and pointOnSeg2D (strictPIP), byte-equivalent
   exact predicates, consolidate to one file-local helper. The PSLG
   suite's copies delete with it. Signed-area consolidation only
   where the code is identical p2(id) loop/cycle form - the
   edge-difference triangle-area sites differ and stay. The
   adjacency-sensitive segment logic inside isLoopSimple2D is
   loop-specific and is NOT merged with anything.
3. DEAD CODE AFTER CUT 1: anything reachable only from the deleted
   suite (helpers, constants, comments referencing it), found by
   compiler + grep, not by guess.
4. TESTS: review establishes the deletable-pin list is EMPTY - no
   test pins the PSLG output suite directly; every island test pins
   surviving behavior (resolve, gates, nested, hazard, public API)
   and stays byte-untouched in its assertions. ONE addition: a
   sub-resolution island pin (an island loop whose doubled area sits
   below the fast-fail threshold gates rather than decomposes) - the
   only validation arm that is deterministically forceable from the
   input side and currently unpinned. Otherwise test changes are
   fixture sharing only: where island fixtures hand-roll the same
   impl-building boilerplate, extract a small shared builder - never
   pin merging or assertion weakening.
5. DOCS: OverlapRemoval.md's validation passage shrinks to the triad
   + the trusted-triangulator statement (the nets-coverage claim is
   dead - see above); InteriorIslandPlan.md gains a
   short addendum recording this diet and why (the plan doc is the
   history; do not rewrite its review record). The PUBLIC comment in
   include/manifold/manifold.h is deliberately unchanged: its named
   fail-closed examples (pinched / coplanar-hazard island loops)
   remain accurate, and the deleted edge-count cap was never part of
   the public contract.

## Out of scope

The mature, TP-hardened core (steps 1-9.5, gates, emit, driver,
member): no restructuring, no perf work, no behavior changes. The
recorded post-landing items (D3 diagnostic family, D5 file split, V1
write-shape) stay recorded - this pass is weight removal, not
architecture. No eps/contract changes anywhere.

## Invariants (the review's checklist)

- Full suite green; the only test deletions are subject-deleted pins,
  each named in the implementation report.
- Fail-closed posture unchanged: every surviving failure arm still
  sets `interiorIslandVerts` and returns; the driver gate untouched.
- The triad's arms (misclassified contour, empty result, degenerate
  or hole-equal triangle) are REVIEWED-NOT-PINNED residuals, stated
  as such: they guard output the trusted triangulator cannot be made
  to produce on valid input, so no deterministic pin exists - the
  review record (island-impl rounds, simplify-plan rounds) is their
  evidence. The implementation must not break their code paths
  (compile coverage + the review diff is the check).
- The resolve pins (seam, public-API, nested, two-island) and all
  EXISTING gating pins (pinched, branchy, vert-sharing, zero-area,
  hazard) are byte-untouched; the sub-resolution pin is the NEW test
  cut 4 adds, not an existing one.
- Determinism untouched (no container or iteration-order changes in
  decision paths).

## Stop conditions

A "dead" helper turns out to have a live caller outside the deleted
suite; a test deletion would remove the LAST pin of any surviving
behavior; any cut that would change an eps, a gate order, or an
emitted result on the existing suite.
