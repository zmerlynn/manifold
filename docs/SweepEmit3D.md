# Sweep-native emission for 3D overlap removal (design, v3)

> SUPERSEDED (history). RemoveOverlaps3D is now the per-component regularization
> operator specified in docs/Regularize3D.md; the v3/v4 sweep pipeline this doc
> describes was deleted from the branch (it lives on in git history and the
> explore/sweep-plane-3d-v3/-v4 branches).  Kept as the design record only.


STATUS: DESIGN DRAFT under crucible review. Supersedes the stage-D/E
architecture of docs/SweepPlane3D.md; stages A-C's proofed mechanics
carry over (softened where stated). Motivation is empirical, from
two full implementations (see that doc's closing records): the
sweep/section/consistency core implemented cleanly twice; the
classify-after-transfer bridge (stage D/E) sprouted independent
strata both times - all of v2's fidelity BREAKs and most of both
simplicity-audit garbage lists live on the bridge. The 2D engine's
own shape is "emit what the sweep retained"; stage D had no 2D
analog. Emmett's #1707 status is "an active cross section of
closed, non-overlapping polygons" - the retained solid section.
This design emits it.

## Architecture

Stage names map to code entry points: CANONICALIZE (Canonicalize),
SEAMS (FindSeams), SLABS (BuildSlabs), STRIPS (EmitStrips), CAPS
(EmitCaps/ComputeCap); the pipeline driver is SweepEmit.  Stage
LETTERS appearing in this doc and in SweepPlane3D.md refer to the
superseded v1/v2 architecture; for reading them against this design:
A = CANONICALIZE, B' = SEAMS, C' = SLABS, D' = STRIPS, E' = CAPS.

CANONICALIZE (formerly A - MERGE; unchanged from SweepPlane3D.md
stage A): eps vert merge, identical-tri multiplicity, zero-mult
drop.

SEAMS (formerly B' - criticals and seams), SOFTENED. The critical
set = sorted x of:
merged verts; edge-face intersection points; seam-seam crossing
points (triple points). Seams (tri-pair intersection segments) are
computed as 3D SEGMENTS for two purposes only: their endpoint x's
and crossing x's feed the critical set, and they serve as TRACKS
for section verts (below). NOTHING ELSE from old stage B survives:
no shared polyline objects, no global triple unification, no
identity write-back, no on-edge subdivision - those existed to give
face PSLGs shared ids, and there are no face PSLGs. Over-inclusion
of criticals is harmless (an extra critical = an extra slab = finer
strips); the diameter guard becomes unnecessary - a fuzzy triple
cluster yields several nearby criticals and thin slabs, handled
uniformly below. Degenerate tri-pair contacts (point/tangent/
sub-eps) contribute their x's as criticals and nothing more. (The
draft's "the chain hazard dissolves" claim was RETRACTED in round
1: features interior to a merged critical run are invisible to both
flanking sections - the explicit SubEpsFeature guard in the
degenerate-slab rule below owns that hazard.)
Coplanar face overlap and edge-in-plane remain OUT OF SCOPE
(detected exactly as today, fail closed): a coplanar pair has no
transversal section story at their shared plane.

SLABS (formerly C' - sections; unchanged mechanics): per slab wider
than eps, at
mid-x every straddling face contributes one directed segment
(sweepDir x outwardNormal convention, signed multiplicity), into
the 2D engine's arrangement + winding passes. The engine returns
the RETAINED boundary pieces (its native output - the piece capture
with (below, above) is no longer needed for emission; retained
pieces suffice) with source ids.

STRIPS (formerly D'). Every retained piece endpoint tracks a 3D
segment:
piece endpoints are (i) face-edge crossings of the section plane -
tracking the input EDGE - or (ii) crossings of two faces' section
segments - tracking the SEAM (both faces' planes contain the seam;
the section crossing at any x IS the seam's point at x). The
enumeration covers classes (i) and (ii); [R1-fold] BLOCK-RULE
FORCED-THROUGH verts are a confirmed THIRD class (empirical: 1% of
endpoints on a triple-point fixture) - a weld places an endpoint on
no linear track. The invariant therefore weakens honestly:
AGREEMENT-UP-TO-CAPS. Class i/ii endpoints extend by evaluating
their tracks at xLo/xHi via the shared Interpolate kernel and agree
BITWISE across the shared critical (empirically confirmed,
maxULP=0). [R2-fold, DEFINED] A forced-through endpoint extends CONSTANT in
(y, z) across its slab - the weld point does not move. Trivially
well-defined; the deviation from the true limit is confined to the
welded slab at the block's eps-scale spread; and BOTH bounding caps
are computed FROM these constant-extended limits, so the deviation
region is cap-covered exactly (the cap arithmetic sees exactly what
the strips emit). Strips are emitted as quads:
with piece endpoints (a0, a1) at xLo and (b0, b1) at xHi in the
piece's directed order, the quad is (a0, a1, b1, b0) split as
triangles (a0, a1, b1), (a0, b1, b0) [R3-fold: vertex order
stated]; orientation is ENGINE-NATIVE: the winding pass emits
retained pieces interior-on-left in section space, which together
with the sweep direction determines the material side of the strip
- no per-face re-derivation (transfer-in-disguise is forbidden).

CAPS (formerly E'). At each critical x = c, the retained REGIONS of
the two
adjacent slabs' sections (each a set of closed 2D loops - the
engine's retained boundary bounds them), both evaluated AT c via
the track extension, generally differ; the difference IS the output
surface content in the plane x = c. Cap = the signed 2D region
difference: material present on the left limit and absent on the
right bounds a cap facing +x, and vice versa. Computed as a 2D
boolean of the two extended retained boundaries - BY THE ENGINE
(seed left loops with +1, right loops with -1, WindRule::Add on
each sign class... precisely: cap_plus = region(L) minus region(R),
cap_minus = region(R) minus region(L), each a standard engine
Subtract). Cap regions triangulate via Triangulate and emit facing
+x / -x respectively. Axis-parallel input faces, the caps of
boxes, winding jumps at shared planes - ALL cap content, no
special path. Closure argument: across slab i's interior the
retained boundary sweeps the strips; at c the symmetric difference
of the two limits is exactly covered by the caps; hence every
strip edge at c is either shared with the neighbor strip (equal
limit locally) or bounded by cap boundary (differing limit).
[R1-fold, ONE SOURCE OF TRUTH PER CRITICAL] The cap's arrangement
pass subdivides the extended limits at their mutual crossings
(empirical: naive edge pairing was 5/8 because cap boundaries are
SUB-SEGMENTS of strip edges) - so the cap arrangement's output
vertex set at c DEFINES the boundary subdivision for BOTH the cap
triangles and the adjacent strips' edges at c: strips take their
c-side polyline from the cap arrangement, not from their own
unsplit extension. [R2-fold, the data flow and lift specified] The
order is acyclic: (1) both slabs extend their retained pieces to c
(tracks for class i/ii, constant for welds); (2) ONE cap
arrangement runs over both extended limits; (3) its output verts -
2D points (y, z) at the plane x = c, i.e. 3D points (c, y, z), no
lift ambiguity - subdivide the cap loops AND the strip edges on BOTH sides: slab
i's right-edge intervals and slab i+1's left-edge intervals are
each replaced by the SAME single arrangement's subdivision of them
- one arrangement at c, three consumers (cap, left strips, right
strips). A constructed subdivision vert lies within alpha of its
host segment (Smith 8.2, the shared kernels), hence within alpha of
the strip's source-face plane - the standard eps-validity, not a
new error class. Closure is then by shared construction, not
assertion. Both cap_plus and cap_minus come from the ONE
arrangement [M4-close]: the engine's winding measure runs twice
over the same collected arrangement (the second with negated
multiplicities - the winding pass constructs no points, so both
measures share the arrangement's geometry exactly); cap regions
triangulate via Triangulate (holes CW per its contract - never
fans); the exterior limit beyond the first/last critical is the
empty region, making the outermost caps ordinary. [M4-close] Two
adjudications from implementation: (a) "its output verts" means the
arrangement's FULL vert set (merged input verts + constructed
crossings), not the retained-edge subset - where the two limits
coincide their edges annihilate and retention is empty, yet the
strips on both sides still need the same subdivision to pair;
the merge/incidence machinery is what guarantees agreement, and
retention measures winding, not incidence. (b) Strip-edge chains
bind per adjacent-built-slab PAIR at the pair's canonical critical
- the FIRST critical of the gap between them (the critical itself
under direct adjacency, the generic case). One arrangement must own
each strip-strip seam: the two arrangements of a sub-eps critical
pair can disagree macroscopically about the subdivision of a shared
(cancelled) edge, since geometry changes discontinuously at a
critical by definition.

DEGENERATE SLABS (width <= eps): no section is built; the slab
contributes no strips. (SUPERSEDED in threshold by STRICT-FP SLAB
BUILDING: unbuilt slabs are now only adjacent-double critical
pairs, so "degenerate" reads ulp-scale, not eps-scale; the cap and
chain rules below are unchanged.) [R3-fold, AMENDED at the emission-closure
crucible] The original rule - every critical gets its own cap,
in-run copies collapsing by exact post-weld identity - was
corpus-falsified: track slopes amplify the in-run x-offsets past
eps and the copies stop being identical.  The standing rule: ONE
cap per adjacent-built-slab pair, at the pair-canonical critical
(ci == li+1, where the strip chains bind); criticals interior to
an unbuilt run emit nothing, since their content is the canonical
cap's L-R difference re-derived at a nearby x.  A macro geometry
change interior to a run remains the recorded fail-closed dead
zone (SubEpsFeature territory), not silently-covered content. One rule, no anchor propagation, no starvation class. [R1-fold - the dissolution claim was WRONG, convergent
finding] A feature living ENTIRELY inside a merged run appears in
neither flanking section and would vanish silently - macroscopic
sub-eps-thin chains are real (the proofed R1'-2 class). ONE guard
returns, explicit and named: any canonical face whose whole
x-extent lies inside a merged critical run and whose area exceeds
its perimeter * eps fails closed as SubEpsFeature (the house
dimensionally-correct threshold: area ~ length * eps separates
band-like sub-resolution geometry, droppable, from macroscopic
sheets, reportable). Smaller faces drop and are counted, exactly as
the eps contract allows, consistently on both sides by
construction.

## What dies / what is born

DIES: face PSLG walk, region-piece correspondence + clearance +
ambiguity, axis-parallel probes and neighbor-slab selection, anchor
components, span guard, seam balance, island/hole routing, global
triple unification + diameter guard, on-edge subdivision, piece
(below,above) capture for emission, the whole stage-D failure
taxonomy (ClassificationAmbiguity, AnchorConflict,
UnclassifiableComponent, Starvation, BalanceViolation, PSLGInvalid).
Failure contract shrinks to: CoplanarOverlap, EdgeInPlane,
SubEpsInput (bbox/eps degeneracy), SubEpsFeature (the merged-run
guard above), NonManifoldEmission (the output gate, kept as the
final honesty check).
[R2-fold] THE ENGINE CONTRACT, stated: SweepWinding gains a
retained-piece out-channel - for every RETAINED boundary piece,
(from, to, sourceId) in section coordinates, emission-oriented
(interior-on-left); this is the existing capture machinery
restricted to retained pieces (the id plumbing is already built and
2D-fence-verified); the (below, above) fields are not consumed by
this design. conflictCount semantics: nonfatal inside the 2D
engine (sourceId = -1 + counter), any nonzero count fatal to the
3D pipeline (EngineIdConflict) - the established rule. [M4-close]
SweepWinding/RemoveOverlaps2D additionally gain an optional
negated-measure out-param (negEdges/edgesNeg): one collected
arrangement, a second winding measure over it with negated
multiplicities, both edge sets materialized against the same vert
list. Null default; 2D callers unchanged.
BORN: track extension (Interpolate at two x's per piece endpoint),
cap construction (one arrangement + two signed measures per
critical), strip/cap assembly. Expected net: a large deletion.

## Eps posture

Unchanged: one absolute eps (EpsilonFromScale budget 1000), input
quantization at CANONICALIZE, slab-width gate, the merged-cap rule
for
sub-eps critical runs. The engine inside slabs and inside cap
booleans runs its standard posture. Alpha covers every constructed
point via the shared kernels.

## Output and tests

Output tessellation is per-slab strips + per-critical caps - finer
than input faces. Strips bounded by class-i/ii tracks lie
geometrically ON their source faces; strips at forced-through welds
deviate within the block's eps-scale spread (eps-valid, Smith 7.7)
- the on-face property is eps-valid, not exact, at welds;
Simplify owns coarsening (performance campaign later). TEST
EVOLUTION: oracle gates 1-5, the fixtures, fences, and pins that
assert PUBLIC behavior (inverted cube, islands P12, touching
coplanar, nested shells) carry over VERBATIM - they are
architecture-blind. White-box pins of deleted mechanisms retire
with their mechanisms; new white-box pins: track-extension bitwise
agreement, cap-boolean correctness on a known critical, merged-cap
rule on a sub-eps cluster. Gate 2 becomes: per slab, engine ingests
the section; per critical, the cap boolean closes against the
adjacent strips (edge-pairing check at the critical plane).

## RISKS (review lanes)

R1 (load-bearing, empirical): bitwise strip agreement across
criticals - track evaluation at shared x from both sides must be
exact; any FP asymmetry breaks closure. And cap booleans on dense/
tangent criticals - does the engine's output over extended limits
actually close against strips (edge-for-edge)?
R2: the track enumeration ("every retained piece endpoint is class
i or ii") - adversarial hunt for a third class (block-rule
constructions, MergeVerticals1D output, engine-constructed points
that track NOTHING linear). Also seam tracks at their ENDPOINTS
(a seam ends inside a slab? no - seam endpoint x's are criticals;
verify no seam-interior slab boundary cases are missed).
R3: cap orientation/multiplicity rigor (nested loops, mult > 1
sections, inverted shells in caps); the closure argument's rigor at
a critical where BOTH topology and geometry change; degenerate-slab
merged-cap soundness (two builds' limits within eps but not equal -
does the cap boolean absorb the eps gap or leak slivers?).
R4: simplicity audit - is anything here transfer-in-disguise; is
the failure contract really this small; test-evolution honesty.

## Crucible round 1 record (synthesis; folds pending)

Empirical (L1): bitwise track agreement CONFIRMED (maxULP=0, 100%
attribution without triples; trackMisses=0 end-to-end on 2-box);
cap-as-engine-Subtract CONFIRMED (exact area at a born critical;
signed seeding works with an aboveInside direction flag).
Convergent fractures to fold in round 2:
1. FORCED-THROUGH TRACKS (L2 + L1: 16/1544 at triple verts): block
   welds create endpoints on no linear track. Fold direction: weaken
   universal bitwise agreement to agreement-up-to-caps - welded
   mismatch regions are eps-thin limit differences the cap boolean
   already sees; state the invariant as "strip boundaries match
   bitwise OR their mismatch is cap-covered", with welded slivers
   dropped/emitted consistently by the cap's own arithmetic.
2. SUB-EPS CHAIN GUARD RETURNS (L2 + L3): features entirely interior
   to a merged critical run appear in neither flanking section and
   vanish silently. Fold: one explicit named guard - any face whose
   whole x-extent lies inside a merged run and whose area exceeds
   the eps band reports SubEpsFeature (fail closed). The dissolution
   claim was wrong; the guard is small and principled.
3. ONE SOURCE OF TRUTH PER CRITICAL (L2 + L1: 62.5% naive edge
   pairing): the cap arrangement's subdivision defines the boundary
   vertex set for BOTH the cap and the adjacent strips' edges at
   that critical - no post-hoc welding, engine-native, closure by
   shared construction.
Also fold: cap_plus/cap_minus both computed; non-convex cap
triangulation via Triangulate (not fans); exterior limits at
first/last criticals explicit; strip orientation stated
engine-native (L3's transfer-in-disguise note); the retained-output
API extension specified against the actual engine; test-evolution
table corrected per L3; L3's implementer-temptation list becomes
the impl prompt's forbidden list.

## Crucible close (round-3 cap, 2026-07-10)

Three rounds, five lanes (one empirical build). SURVIVED with
empirical validation: the architecture's core - per-slab sections
through the engine, class-i/ii track extension (bitwise, maxULP=0),
cap-as-engine-Subtract (exact area), zero track misses end-to-end
on 2-box, acyclic cap dataflow, the (c,y,z) lift, the weld
extension's closure story including neighbor-disappearance and
first/last-slab cases (all HELD on re-attack). OPEN AT CLOSE (the
round-3 list, all bounded):
1. Weld strips are NOT on-face (constant-(y,z) extension); the
   on-face claim must weaken to within-block-spread eps-validity at
   welds, everywhere the doc asserts it.
2. The merged-run rule's "criticals all within eps" is false for
   chained runs (gaps <= eps, span unbounded); needs per-critical
   caps inside runs or an explicit run-span bound feeding
   SubEpsFeature.
3. One-source-per-critical prose still has one both-directions
   sentence to make unambiguous.
4. The engine contract needs conflictCount semantics.
5. Strip quad vertex order unstated.
RECOMMENDATION recorded for the owner: fold these five (half-day),
skip a fourth review round - design review has reached diminishing
returns relative to the oracle-anchored implementation gates, which
are the stronger fence and already exist. Then implementation
crucible on this branch, with the L3 implementer-temptation list as
the forbidden catalog.

## Implementation crucible close (round 3 cap, 2026-07-10, 25f7c84a)

Round trajectory - simplicity judge: BREAK/21 (v1) -> BREAK/21 (v2)
-> BREAK/13 -> NEED-CHANGE/13 (first non-rewrite verdict across
four implementation generations); fidelity: 7 findings -> 2. Green:
31 + 1 scope-skip gates incl. the new triple-critical and 3-solid
oracle fixtures, fences 111, full suite 577, all four round-2
mechanism BREAKs implemented and verified. REMAINING (both in the
cap construction, one afternoon's scope): the emittedPairs
duplicate-skip is a run-merge in disguise (replace with the spec's
per-critical caps computed at each critical's own x); M4 still runs
two arrangements over shared input and strips still self-extend -
the spec's literal one-arrangement-three-consumers remains
unimplemented. Judge's 13 items stand listed in the audit answers.
The architecture verdict after three implementations: the
sweep-native emission core is sound and each round's findings
narrowed monotonically; what remains is finishing one function
cluster to its spec, not another rewrite.

## M4 close (main-agent, 2026-07-10)

The one-arrangement-three-consumers cluster is implemented to spec
after escaping three subagent lanes. Engine: optional negated
second winding measure over the one collected arrangement (nulled
for 2D callers; negation commutes with collect - verified by
bitwise-identical cap output against the old two-call scheme on the
box+rotated fixture). ComputeCap: ONE RemoveOverlaps2D call;
cap_plus from the positive measure, cap_minus from the negated one;
per-piece strip chains carry the arrangement's vert positions
bitwise (snapped endpoints via the input-vert map, interior splits
by projection onto the snapped chord). Strips stop self-extending:
extension now happens once, as cap input; EmitStrips zips
cap-produced position polylines. Died with the two-call scheme: the
capVerts bag + eps-dedup collection, the fuzzy x-keyed cap lookup,
the eps*0.1 neighbor-slab tolerance (slab lookup is index-exact:
slab bounds ARE the criticals), the param-space zipper's 2.0/1e-10
sentinels, the duplicated linear face-track searches (one map +
assert), and the fake criticals-as-vertices pattern (degenerate
contacts and seam-seam crossing x's now live in
ArrangementGeometry.criticalXs; a critical is not a vertex).

The two seam adjudications recorded in-line above ([M4-close] tags)
came from a real red: per-critical chain binding broke the
strip-strip weld across sub-eps critical pairs, and the
retained-subset scan emptied exactly at cancellation seams. The
invariant that makes per-critical caps survive sub-eps pairs: a
non-canonical cap of a run is sliver-only and collapses in the
assembly weld.

KNOWN DEAD ZONES at close (recorded, not regressions - all
inherited by construction and unreachable in the suite): (a) a
macro-scale geometry change at a non-canonical critical of a
sub-eps run (a face starting mid-run with a full yz edge) emits a
macro cap whose corners no strip carries; candidate SubEpsFeature
tightening. (b) A bare merged vert (every incident edge cancelled)
lying eps-on a retained cap edge is not incidence-split into that
edge but does subdivide the strips. (c) Multi-sliver runs (total
gap > eps) weld strip-to-strip across more than eps. Gate5
box+rotated is the de-facto pin for the canonical-binding rule (it
reds under per-critical binding); the M1 pin now asserts the
criticalXs mechanism and is mutation-verified against a stubbed M1
loop.

## Post-M4 audit round (main-agent, 2026-07-10)

Both closing audits returned NEED-CHANGE with converging short
lists (simplicity 13 -> 9, "no longer reads like a rewrite
candidate... recognizably one design"; fidelity: M4 open on one
point). Fixed after my own validation of each finding:
- Engine: `verts` under the negEdges contract now includes EVERY
  collected-arrangement vertex (retained-only materialization
  missed crossings whose four quadrant windings share one sign) -
  the spec sentence "its output verts" is now literal.
- Fail-closed release arms: missing piece-track attribution fatals
  as EngineIdConflict; a bad triangulator index fatals as
  NonManifoldEmission (both were assert-then-continue).
- The coplanar-overlap threshold is dimensional (area >
  perimeter * eps, the house rule) - the old eps^2 was FP-noise
  fragile for shared-vertex coplanar pairs.
- SeamSeamCrossX's near-parallel gate is dimensional (|cross| >
  eps * (lenA + lenB)) - the old eps^4 admitted pseudo-crossings
  of near-parallel seams at meaningless x (over-inclusive but
  noisy).
- Seam-track slab membership is exact (endpoint x's are criticals;
  xMid is strictly interior); InterpolateSafe documents its named
  extrapolation role and asserts x-nondegeneracy.
- Both degenerate-contact arms record endpoint x's in criticalXs.
- The M1 pin cross-checks criticalXs against brute-force seam-seam
  crossings; Pin_OneArrangementPerCritical pins the M4 count
  structurally via the capArrangements counter (a two-call
  regression doubles it).

DEFENDED against the audits (recorded disagreements):
- The shared-edge exemption cannot move after coplanar detection:
  opposite-diagonal triangulations of two solids touching on a
  common plane are legal (TouchingDisjoint pin) and locally
  indistinguishable from a folded flap; the hazardous same-winding
  case fails closed downstream (coincident section edges with two
  source ids and nonzero net multiplicity -> EngineIdConflict; the
  legal case cancels to net zero before any conflict). Empirical:
  the suggested reorder reds TouchingDisjoint.
- BuildImpl's exact-duplicate drop is the [R3] "their differences
  are empty" sentence realized (per-critical caps of a sub-eps pair
  computing the same macro difference twice), not a run-merge.
- The class-ii/weld eps-inference stands per the [R2-fold]
  adjudication.
- ChainSplitVerts remains a geometric on-chord scan (the
  tolerance-model reading of "verts subdivide the strip edges");
  per-input-edge chains from the engine would upgrade it to exact
  identity and are recorded as future work.

CONVERGENCE ROUND (same day): fidelity SURVIVE - M4 RESOLVED with
probe evidence (internal crossing verts materialize under negEdges;
default path unchanged), both new pins verified structural, no new
findings. Simplicity NEED-CHANGE narrowed to three release/drop
paths ("this is one design now, not patch strata"), all three
fixed: TriTriSeam no longer conflates point/tangent contact with
disjoint (an eps-admitted contact returns a collapsed seam and its
x enters criticalXs through the ordinary degenerate arm);
OutEdgesToPolygons gained a checked-extraction out-flag and a
non-closing cap boundary walk fatals as NonManifoldEmission instead
of silently dropping a loop; EmitStrips validates the
one-non-empty-chain-per-piece contract and fatals on violation
(ZipperEmit's empty-chain tolerance deleted with it).

ARC CLOSE AT ITERATION CAP (three revise rounds). The judge's
round-3 list was ONE item - EmitStrips validated chain sides
against each other but not against the piece count - fixed
post-verdict (both sides now check pieces.size() directly; equal
truncation cannot slip through). Verdict trajectory across the
implementation generations: BREAK/21 -> BREAK/21 -> BREAK/13 ->
NEED-CHANGE/13 -> NEED-CHANGE/9 -> NEED-CHANGE/3 -> NEED-CHANGE/1
(fixed). Fidelity: closed at SURVIVE. The adversarial reviewer is
prompted to always find something; the honest close criterion here
is the trajectory's convergence and the fidelity SURVIVE, with the
round-3 residual fixed and every defended item carrying its
recorded rebuttal above. Owner review is the next gate.

## Style crucible (main-agent, 2026-07-10)

A behavior-preserving simplicity/style/naming pass, two review
rounds (style judge, a dedicated linalg/reuse hunter, a naming
lane).  The stages got their names (the letters were sequential
noise); house style per the manifold-style-review skill landed:
linalg does its job (length2, lerp, 2D cross/perp forms, swizzle
blends - each verified bitwise against linalg.h before folding),
manifold::Box replaced a hand-rolled AABB, the seam-cluster helpers
return optionals instead of trailing out-refs, casts and
anonymous-namespace statics normalized, stale comments (including
an M1 pin banner describing the pre-criticalXs mechanism) rewritten,
and the Gate4e/4f names stopped lying (FailClosed; assertions
untouched).  The reuse lane converged to SURVIVE.  Residuals, both
deferred as upstream follow-ups rather than churned here: the
upstream boolean2.cpp's reopened anonymous namespaces, and its Box2
vs the public Rect.

## COPLANAR (extension design, round 2 - post D1/D2 fold)

The coplanar family comes IN SCOPE: coplanar face overlap and
edge-in-plane contacts resolve instead of failing closed.  CAD
inputs hit this constantly (stacked solids, shared walls).

PROBE EVIDENCE (scratch driver, fatals bypassed, before design):
identical-face stacking already resolves exactly (CANONICALIZE
annihilates coincident opposite faces; caps do the rest).  A shared
plane PERPENDICULAR to the sweep resolves exactly with no new
mechanism (pure cap arithmetic).  A shared plane PARALLEL to the
sweep with ANTI-oriented overlap - stacking and walls, the dominant
class - resolves exactly and conflict-free: +1/-1 collinear section
segments annihilate in PolySetAdd before any srcId merge.  The
mechanism gaps: SAME-oriented overlap (+2) loses attribution, and -
the round-1 BREAK, validated by instrumentation - junctions bounded
by in-plane structure have CRITICALS but NO TRACK, so their
extensions weld-freeze at xMid while the far side lands exactly on
the 3D vertex (hull fixture: an 85-eps open cap walk, CLASS-I vs
WELD across one critical).

DESIGN, five mechanisms:

1. COPLANAR GROUPS (SEAMS).  A grouping PRE-PASS unions canonical
   faces by plane-within-eps (either orientation; the existing
   detector's test is the membership test), INCLUDING shared-edge
   pairs - the flap's pair must join its group before any exemption
   skips it [D1-3].  Grouped pairs skip seam creation; the
   CoplanarOverlap fatal retires.  Group ids occupy nFaces + k.

2. GROUP-ID SEEDING (SLABS; zero engine change).  A grouped face
   seeds its section edges with its GROUP id; multiplicities stay
   per-face.  Coincident group segments merge under one id - no
   conflict exists to record.  Anti-oriented content cancels;
   same-oriented content sums and the winding pass REGULARIZES it:
   a 0 -> +2 boundary emits ONE unit-multiplicity retained piece
   (empirically verified; the |m| materialization loop never sees
   the 2) [D2-1].  Strips remain boundaries of the positive region
   only; the negated measure exists for caps alone [D2-2].
   A conflict that still reaches -1 means a face pair coincided in
   section WITHOUT sharing a plane group: geometrically impossible
   for genuinely non-coplanar planes (they meet in a line, sections
   cross at a point), so it indicates near-coplanar geometry the
   grouping eps-test missed - a named guard (EngineIdConflict
   narrows to this), not an ordinary residual class [D1-4].

3. PER-VERTEX EXTENSION RESOLUTION (STRIPS/CAPS; replaces the
   per-piece cascade - the round-1 BREAK's fix).  Extensions resolve
   ONCE PER SECTION VERTEX; every incident piece endpoint reuses
   the resolved point, so same-side closure holds by construction.
   The track candidate classes, in order:
   (i) face cutting edges (class-i, existing);
   (ii) seam tracks (class-ii, existing);
   (iii) IN-PLANE BOUNDARY EDGES of coplanar-group members - the 3D
   edges (lying in the shared plane) whose section crossings create
   the group's interior junctions; ordinary edge interpolations, new
   candidate class;
   (iv) weld constant - shrinks to genuine block-rule artifacts;
   the cap-coverage argument for weld deviation is unchanged.
   3D-IDENTITY PREFERENCE: when the governing track terminates at
   an arr.verts vertex at the target critical, the extension IS
   that vertex's (y,z) - exact identity, never an eps snap (the
   85-eps gap is far beyond snapping) [D1-1].  Ties across
   candidates resolve deterministically (lowest class, then lowest
   member id); an unresolvable vertex fails closed [D2-3].

4. IN-PLANE SKELETON CRITICALS (SEAMS; the M1 analog, widened in
   round 2).  Each group's IN-PLANE SKELETON = its members' boundary
   edges PLUS the in-plane seams of transversal faces with members
   (a seam with a member lies in the shared plane by construction).
   Pairwise skeleton-crossing x's enter criticalXs: edge x edge
   crossings (the original mechanism) AND seam x edge crossings -
   the event where a transversal junction slides across a group
   breakpoint, which neither M1 (seam x seam sharing a face) nor
   edge x edge covered [D3-R4].  Axis-aligned fixtures mask the
   class (crossings sit at vert x's); a rotated coplanar pair plus
   a rotated transversal pin it red-first.  Benign non-events,
   enumerated: in-plane tangency (the endpoint is a vertex, already
   critical); collinear-overlapping member edges (combinatorics
   change only at their endpoint verts; between them the coincident
   breakpoints share one 3D line and the tie rule is deterministic);
   group membership along x (per-face and x-independent).

5. CHECKED CLOSURE (already landed, the backstop): open retained
   walks in any cap fail closed as NonManifoldEmission.  Gate4c's
   hulls become MustResolve only when they empirically resolve
   under mechanisms 1-4; until then the open-walk guard is the
   honest boundary and the gate's skip narrows to it.

EDGE-IN-PLANE: the fatal dies with the family.  T-contacts section
to T-junctions the engine's incidence pre-split owns; the F-G seam
is an ordinary class-ii track; probe evidence shows the stacking
family resolving exactly through such contacts.

RE-ADJUDICATION (supersedes the style-round note): the same-winding
shared-edge flap becomes ordinary +2 content and RESOLVES to a
regularized single cover.  The fails-closed-via-EngineIdConflict
story retires with the conflict it relied on.  Requires a REAL
closed-shell folded-flap fixture (a single solid whose shell folds
back over itself), asserting single-cover output [D2-4].

GATE EVOLUTION: Gate4d's near-parallel plates group and RESOLVE
(adversarial attack failed to break this: dv = 3.5e-15 against the
oracle); the gate evolves to resolve-with-oracle-or-named-guard - a
strictly stronger accepted outcome.  Gate4c per mechanism 5 above.

NEW TEST SURFACE: oracle gates for stacked-perpendicular,
stacked-parallel, shared-wall, same-oriented partial overlap; a
rotated coplanar pair (mechanism 4, red-first); a three-face group;
inverted-shell and mixed-orientation groups (+ + -, - - +, and
subtract-encoded inverted stacking) [D2-2]; the closed folded-flap
shell; unit-boundary regularization of +2 content pinned
explicitly [D2-1].

EPS BOUNDARY [D3-R5]: group membership is the symmetric OR of the
two plane-direction tests (the pairwise detector's short-circuit
restructures into the grouping pre-pass).  Geometric coherence of
eps-CHAINED groups (pairwise within eps, wider in total) is
DELEGATED to the engine's vert merge - itself an eps-union-find over
the same geometry, so section verts chain-weld exactly where faces
chain-group; class-iii tracks are per-member true 3D edges, so no
common-plane projection distorts a wide group.  No diameter guard
(the old design's served triple-point unification, a dead
mechanism); a three-face eps-chain fixture pins the delegation.
The band just OUTSIDE eps (razor wedges, ~2-10 eps separation) is
recorded as UNTESTED: not grouped, no exact merge, no conflict
guard fires - thin-sliver behavior is empirical; a band fixture
joins the test surface and its outcome decides whether the
grouping threshold widens or a named band guard lands.

TEST SURFACE ADDITIONS [D3-R6]: a coplanar overlap corner within
eps of a vertex plane (group x sub-eps critical pair x
pair-canonical binding); the rotated skeleton pin (mechanism 4);
the razor-band fixture; the eps-chain fixture.  Exterior-cap
groups are covered by the stacked-at-extreme gates; coplanar + M1
interplay is subsumed by the skeleton rule.

DESIGN CRUCIBLE CLOSE (round 2): D1 (engine reality) SURVIVE after
the fold - its own instrumented hull junction resolves through
class-iii + 3D-identity to one arr.verts vertex from both sides;
class-iii candidates enumerable by member-edge scan, no hidden
in-plane arrangement.  D2 (winding) folded: unit-boundary
regularization stated as the invariant, strips positive-only,
mixed-orientation and folded-flap fixtures required, Gate4d
evolution attack failed (dv 3.5e-15 vs oracle).  D3 (criticals/eps
boundary; executed by the orchestrator after two lane deaths at
the output cap): the skeleton-crossing gap folded into mechanism
4, eps-boundary items above, no BREAK.  No unrefuted BREAK
remains; implementation is the next phase, with Gate4c's hulls as
the empirical acceptance fixture and every mechanism pinned
red-first per house discipline.

## COPLANAR implementation close (main-agent, 2026-07-10)

All five mechanisms landed in one pass; the four core oracle gates
went green immediately (including same-oriented +2 and the +3
triple group), and the extended surface followed: mixed-orientation
planes, the eps-chain (engine-merge delegation held), inverted
stacking, and the razor band resolving under its recorded contract.
TouchingDisjoint and the inverted-cube pins held throughout;
Gate4d resolves with oracle agreement as designed.

FINDINGS DURING IMPLEMENTATION, recorded:
- Mechanism 4 is the belt, seam endpoints are the suspenders: for
  closed solids, every in-plane boundary-edge crossing is also the
  endpoint of a side-face seam, hence already a vert.  The skeleton
  pin asserts the property (crossings are criticals however
  delivered) rather than a vertex-free entry; a mechanism-4 stub
  would survive it on closed-solid fixtures - recorded as a pin
  limitation, with the degenerate-adjacent construction (where the
  belt matters) as future test surface.
- Edge-on-face TOUCHING (the former EdgeInPlane P4 pins) welds into
  a genuinely non-manifold union - four faces around the contact
  line.  BuildImpl gates on the full 2-manifold check (Is2Manifold,
  vertex links included).  SUPERSEDED same day: the sheet splitter
  below resolves the class; the gate remains for what it cannot
  pair.
- The in-run macro-change dead zone from the M4 close is now
  REACHABLE: two perpendicular faces sub-eps apart produce macro
  cap content at both criticals of a sub-eps run and fail closed at
  the output gate.  Recorded-contract fixture pins named-guard-or-
  oracle; SubEpsFeature tightening moves up the priority list.
- Gate4c hulls fail closed at the near-coplanar guard: real
  geometry has faces whose planes cross at shallow angles -
  sections coincide within eps locally while the global
  three-verts-within-eps grouping test misses them.  The gate's
  skip narrows to {EngineIdConflict, NonManifoldEmission}; closing
  the hulls needs either local (per-section) conflict tolerance or
  overlap-region-scoped grouping - the next arc's opening problem.
- CoplanarOverlap and EdgeInPlane retired from the taxonomy; the
  detectors and their helpers (the coplanar clip, the edge-interior
  length test) deleted; the same-group skip and the grouping
  pre-pass replaced them.  The closed-shell folded-flap fixture
  remains open test surface (the +2 semantics are covered by the
  same-oriented and triple-group gates; a valid hand-authored
  folded shell is still owed).

## Touching contacts resolved (main-agent, 2026-07-11)

The "topological welding" note is implemented as a SHEET SPLITTER in
assembly, before topology construction.  Target semantics probed
first: Boolean3's own union of edge-touching cubes keeps the
coincident verts as separate topological copies (Decompose() = 2) -
the epsilon-valid posture; the regularized union of solids touching
on a measure-zero set is the solids, separate.

Mechanism (SplitTouchingSheets): (1) radial pairing per fan edge -
sort incident faces by angle around the edge axis; from the
outward-normal convention a forward halfedge carries material just
below its angle and a backward one just above, so material wedges
alternate with empty ones and each backward halfedge pairs with the
next forward one CCW; (2) vertex split - corners union through
PAIRED halfedges only, one output vert per connected component.
Covers edge-on-face, edge-edge, vertex contact, and self-touch
(same-component arms) uniformly.  Fail-closed arms (returning to
NonManifoldEmission): unbalanced halfedge counts, sliver third-vert
on the edge line, radial ties (tangent sheets: either pairing is a
coin flip, and the WRONG pairing preserves volume while garbling
topology - which is why the fixtures assert component COUNT, not
just oracle agreement), and non-alternating patterns (overlapping
material).

P4/P4b evolve from pinning the fatal to pinning RESOLUTION (oracle
vs Boolean3, two components); new edge-edge and vertex-only cube
fixtures.  Vertex-only contact needed no splitting (a shared vert
with disjoint fans and no shared edge is already tolerated
topologically) but the splitter separates it anyway - matching
Boolean3's output shape.  Gate4c is unaffected (its guard fires
before emission).

## Corpus measurement (2026-07-11)

The real-mesh corpus (test/models: four left/right operand pairs,
the hull pair, four Offset singles, openscad-nonmanifold-crash,
self_intersect A/B) through RemoveOverlaps3D, one subprocess per
case, timeout-bounded, oracle = Boolean3 where it round-trips.

RESOLVED, ORACLE-TRUE: the Cray pair (the first real operand pair
end-to-end) and three of four Offset singles.

FAIL-CLOSED, by class:
- NEAR-COPLANAR GUARD (EngineIdConflict): the hull pair and one
  Generic_Twin pair.  The known class: planes crossing at shallow
  angles coincide in section locally where the global grouping
  test cannot see them.
- EMISSION-CLOSURE DEFECT (NonManifoldEmission, unbalanced fan):
  one Generic_Twin pair and Havocglass8.  Instrumentation shows
  the sheet splitter's UNBALANCED arm (1F/0B, 1F/2B): the emitted
  triangulation itself has an unpaired edge - a missing or extra
  strip/cap triangle upstream, not a splitter limitation.  A new,
  distinct correctness residual on real geometry.

TIMEOUT (the performance wall, deferred by design): Offset1,
openscad-nonmanifold-crash, self_intersect A/B.  Diagnosed, not
guessed: the openscad case builds slabs in seconds but produces
half a million retained pieces across ~9k slabs (24k criticalXs
from CSG-dense in-plane structure), and assembly's quadratic
weld does the rest.  The perf campaign's first targets are known:
spatial-hash the assembly weld, bound the chain-split scans, and
consider criticals thinning.

The histogram's reading: the pipeline is correct-or-honest on
every case - nothing resolved wrong (zero ORACLE-OFF) - and the
residual work splits three ways: the near-coplanar guard (known),
one new emission-closure defect class (new, needs a minimized
fixture), and the deferred performance campaign.

## Emission-closure crucible (main-agent, 2026-07-11)

The corpus' unbalanced-fan class, root-caused through five instrumented
iterations on the Havocglass8 fixture (each layer peeled revealed
the next):

1. DOUBLED IN-RUN CAPS (fixed): per-critical caps inside a sub-eps
   run re-derive the same L-R difference; the [R3] premise that the
   copies collapse by exact post-weld identity is FALSIFIED when
   track slopes amplify the sub-eps x-offset past eps.  Fix, kept:
   caps emit only at PAIR-CANONICAL criticals (ci == li+1) - the
   in-run cap's content is the canonical cap's, re-derived; one cap
   per adjacent-built-slab pair.  Pin_OneArrangementPerCritical
   pins the rule.
2. RESOLVER TWIN DIVERGENCE (fixed): the winding pass constructs
   T-junction verts exactly and never merges them, so sections
   carry sub-eps-adjacent twins; the per-exact-vertex resolver
   could match twins to DIFFERENT tracks, and divergent steep
   tracks stretch a sub-eps section gap to many eps at the cap.
   Fix, kept: the resolver resolves per EPS-CLUSTER (one ball, one
   track - the same one-entity rule as the vert merges).
3. CROSS-DERIVATION NOISE FLOOR (kept, with honest scope): cap
   arrangements and chain subdivision run at capEps = 8 eps - the
   un-amplified kernel-noise floor for twice-constructed inputs.
   This is a floor, not a bound: it does NOT close the class.
4. THE IRREDUCIBLE REMAINDER (diagnosed, named, deferred): on real
   steep-track geometry the L and R limits can DISAGREE about a
   junction's position by 10-100 eps (observed: two retained loops
   splitting hairs across a 10-eps micro-edge; L's boundary through
   one twin, R's through the other).  The closure argument's "equal
   limit locally" premise fails at sub-cluster scale, and no
   constant radius fixes it - the amplification is unbounded.
   Closing it requires PROVENANCE-EXACT closure: strips consuming
   the retained cap graph's edges per input piece (with
   partial-retention semantics) instead of geometric chain
   reconstruction - the recorded per-input-edge upgrade, now
   promoted from future work to the demonstrated requirement, and
   the natural companion of the hulls' near-coplanar arc.

The corpus fixtures stand as recorded contracts
(Corpus_*_Recorded): a named guard or a true resolve, never silent
garbage.  Full suite green with all three landed changes.

## PROVENANCE CHAINS (main-agent, 2026-07-11)

The per-input-edge subdivision upgrade named at the emission-closure
crucible is implemented.  The mechanism and its adjudications:

THE CHANNEL (engine).  RemoveOverlaps2D gains an optional out-param
`edgeSubdiv`: one ordered vertex-POSITION polyline per INPUT edge -
that edge's exact arrangement subdivision, from its start merged vert
through every interior arrangement vertex on it to its end merged
vert, oriented v0->v1.  It is keyed by a COINCIDENCE CLASS = the
unordered merged-endpoint pair of the original input edge (a pure
label `classId` riding PolyVal/SweepEdge exactly as srcId does; -1 for
all 2D callers; null-defaulted; the arrangement geometry is BITWISE
unchanged - classId never touches m, erase, or conflict).  Interior
verts come from two sources: the incidence pre-split (seeded in the
driver) and the arrangement pass's own crossings + block-rule
forced-through points (recorded in SplitAt/ProcessEvent against the
edge's class; MergeVerticals1D preserves the class through the
vertical resolve).  The 3D cap (ComputeCap) consumes edgeSubdiv
directly as the strip chains (L pieces then R pieces, index-aligned),
retiring the geometric on-chord projection (ChainSplitVerts, deleted).
Pin: Boolean2.EdgeSubdivProvenanceChannel (coincident-class agreement
+ crossing capture, mutation-verified red against an endpoints-only
stub).

WHY CLASS-KEYING answers the design questions.  (a) partial
retention: the subdivision covers the WHOLE input edge (incidence
breakpoints + crossings), independent of which sub-edges are retained,
so a strip subdivided where a cancelled neighbor breaks it still pairs
with that neighbor.  (b) coincident chords (L + R anti-oriented, same
merged endpoints): identical class -> identical interior sequence by
construction (the twin-selection ambiguity that plagued the geometric
projection cannot arise - a twin on a DIFFERENT line has a different
class and is excluded; a real crossing is included by provenance).
(c) a vanished piece (merged endpoints equal) yields a single-vert
chain.  (d) determinism: std::map ordering + deterministic seed order;
2D fences (111 + the new pin) unchanged, default paths zero-cost.

WHAT IT FIXED, AND THE WALL IT EXPOSED.  The chains close the
UNPAIRED-EDGE class (Havocglass8's first-hit failure: strips that
mis-subdivided relative to the retained cap graph now consume it
exactly).  They do NOT close the corpus gates, because the true
blocker is one layer deeper and ORTHOGONAL to subdivision: at a
critical the L and R extensions of ONE junction can land 60-200 eps
apart (instrumented on Havocglass8: a forced-through WELD frozen at
its section position vs its track-extended twin landing on the
junction; the section itself already carries the pair ~60 eps apart,
below the resolver's eps cluster, so one finds a track and the other
welds, and the extension amplifies the gap).  Design question (b)'s
premise - "coincident L/R map to ONE arrangement sub-edge" - REQUIRES
those endpoints to merge; at 60-200 eps > capEps they do not, so L/R
become two sub-edges with a MICRO-EDGE between them, and no
subdivision scheme repairs a POSITION divergence.  The cap arrangement
keeps the micro-edge -> sliver cap triangles (Havocglass8, k=4
material-overlap) or a strip-less cap edge (GenericTwin7863, 1F/0B),
each failing closed at the sheet splitter.  This is exactly the
crucible's "irreducible remainder": no constant radius fixes it
(capEps is a floor, established prior), so closing it needs a
NON-constant-radius provenance junction unification - a true
3D-IDENTITY extension that snaps L and R endpoints of one 3D junction
to the same point (the coplanar mechanism-3 preference, unimplemented
in Extend, which currently interpolates and never snaps).  That is the
research-grade next step and the natural companion of the hulls'
near-coplanar arc; the provenance chains are its necessary
foundation, landed and pinned.  Corpus gates remain recorded
contracts, with the residual now named as twin-position divergence
rather than mis-subdivision.

## 3D-IDENTITY EXTENSION (main-agent, 2026-07-11)

The twin-position divergence named at the PROVENANCE CHAINS close was
attacked with the coplanar mechanism-3 "3D-IDENTITY PREFERENCE": snap the
divergent extensions of one junction to their SHARED CANONICAL 3D point.
This section records the instrumented finding, the mechanism, the A/B/C
adjudication, and the honest residual - the mechanism did NOT close the
corpus gates; the residual is deeper than divergence and stays recorded.

THE INSTRUMENTED JUNCTION (Havocglass8 + GenericTwin7863, env-gated dumps).
Every failing junction is a NEAR-DEGENERATE arr.verts cluster:
- The cap critical c sits SUB-EPS from the true 3D vertex V (V.x - c is
  0.3-0.5 eps - a sub-eps critical pair, so the slab between c and V.x is
  unbuilt).  A section vert on a STEEP track (observed slopes ~6 to ~125)
  amplifies that sub-eps x-offset into 60-200 eps at c: some incident tracks
  land on V, others extrapolate 60 eps past it, and forced-through welds
  freeze at their section position.  The cap keeps the spread as a micro-edge
  -> sliver/unbalanced fan -> fail closed.
- The SHARED CANONICAL is V itself: an arr.verts vertex, computed identically
  by L and R (same arr.verts).  Snapping every image of V to V.yz collapses
  the micro-edge by construction.

THE MECHANISM EXPLORED (candidate A, cap-input snap).  For each cap at x=c,
snap each extended cap-input vertex to the nearest arr.verts vertex V with
|V.x - c| <= eps (V on the cap plane), using a VORONOI-SAFE test: snap only
when the nearest such V is at least twice as close as the second-nearest
(dist < 0.5 * secondDist).  The Voronoi ratio is a PER-JUNCTION radius (half
the local vertex spacing), never a global constant - it absorbs an unbounded
junction spread where V is isolated and REFUSES to snap in dense areas.
Empirically this is SAFE: the full non-corpus suite stays green at an
effectively unbounded snap radius (zero regressions on 47 fixtures), because
the ratio guard alone prevents over-merge - the established constant-radius
dead end is avoided.

WHY IT DOES NOT CLOSE THE GATES (the refined wall).  Input-snapping resolves
the PRIMARY divergence but exposes two DEEPER residuals:
1. CLUSTER-COLLAPSE RE-EMISSION.  Collapsing a junction's images onto one V
   makes incident pieces VANISH (both endpoints -> V) and can FLIP the cap
   arrangement's topology (moving a vertex 60 eps crosses a neighboring edge,
   changing which retained region an edge bounds).  The cap arrangement runs
   over ALL cap input, so a local snap perturbs the boundary GLOBALLY -
   GenericTwin7863's snap fixes its steep-track micro-edge but opens a macro
   hole (an unbalanced fan) elsewhere on the same plane.  Load-bearing cap
   slivers cannot simply be dropped (that opens macro holes directly).
   Closing this needs COORDINATED cap+strip re-emission around a collapsed
   junction - a local re-mesh, not an input snap.
2. DISTINCT-VS-SPREAD AMBIGUITY.  Dense clusters carry MULTIPLE true vertices
   within a few hundred eps on one (sub-eps-adjacent) cap plane
   (Havocglass8: two canonical vertices ~300 eps apart, plus a near-x-
   degenerate track).  Their tracks already land on their own vertices - the
   defect is the emission AROUND two near-coincident-but-distinct vertices,
   and "which vertex does this image belong to" is ill-posed at the
   construction-noise scale.  This is the same eps-scale disambiguation the
   near-coplanar arc faces.

ADJUDICATION (A vs B vs C).  (A) 3D-identity snap is the correct shape for the
PRIMARY divergence and the shared canonical is unambiguous (arr.verts V), but
input-snapping alone is insufficient (residual 1).  (B) cross-side cap
unification shares A's input-snapping limits and additionally cannot
distinguish spread-of-one from distinct-near-degenerate vertices (residual
2 - the core hard problem).  (C) extend-through-the-junction has no single
constructible crossing in a dense cluster, and welds/degenerate tracks
self-locate nothing.  All three reduce to the same eps-scale re-emission and
disambiguation problem.

RESIDUAL, HONESTLY.  The Voronoi-safe 3D-identity snap is a SAFE, principled
partial that resolves the primary steep-track divergence, but the corpus
gates need coordinated re-emission around collapsed/near-degenerate junction
clusters and eps-scale spread-vs-distinct disambiguation - the research-grade
remainder, the companion of the near-coplanar arc.  The exploratory snap is
NOT landed (it trades sliver-overlaps for holes and its correctness pin
cannot be met via the gates); the mechanism and its instrumented evidence are
recorded here so the re-emission successor can rebuild on it.  Corpus gates
stay recorded contracts, with the residual refined from "twin-position
divergence" to "near-degenerate cluster re-emission".

## Corpus fail-closed attribution (main-agent, 2026-07-11)

NOTE (superseded by STRICT-FP SLAB BUILDING): the dominant dead-zone-c class
below is a consequence of the eps-width slab gate merging sub-eps runs.  Under
the strict-FP gate those runs build, so dead-zone-c does not arise; the residuals
move (openscad -> NonManifoldEmission, self_intersect A/B -> ArrangementBudget).
The steep-track wall-A fan and the arrangement-density budget remain the true
residuals.  Read the STRICT-FP section for the current class map.

The four SINGLE-mesh corpus fixtures (Offset1, openscad-nonmanifold-crash,
self_intersect A/B) all fail closed at NonManifoldEmission / "unresolvable
sheet contact" (the sheet splitter's arm).  The performance campaign, which
made them terminate, ASSUMED this was the steep-track junction-spread residual
(the "3D-IDENTITY EXTENSION" wall).  Instrumentation refutes that: the
dominant defect is a DIFFERENT, already-recorded class.  An adversarial
verification lane independently re-instrumented two of the four cases
(2026-07-12) and confirmed that class, with two corrections folded in below:
the first-failure nuance and the control scoping.

WHAT THE INSTRUMENTATION SHOWED.  Audited over the whole EMITTED SURFACE, the
defective edges are overwhelmingly BOUNDARY edges (a single halfedge = an open
hole), and those holes are MACRO (length orders of magnitude above eps,
whole-face scale), not the micro (eps-to-hundreds-of-eps) edges of the
steep-track wall - on the order of a thousand cap-plane holes per re-run case.
Every sampled hole carries the same displaced-partner signature: an open
cap-plane endpoint with a PARTNER vertex at the same (y, z) but a few eps away
in x, across a skipped run.  That partner is the far-side strip endpoint of a
SUB-EPS RUN whose TOTAL x-width exceeds eps: the assembly weld (3D distance
within eps) cannot fuse two points more than eps apart in x, so the shared
cap-plane loop tears into open boundary edges.  This dominant emitted-surface
defect is distinct from the literal FIRST edge the splitter returns (below).

THE CLASS: dead-zone (c), CHAINED-RUN STRIP DISPLACEMENT.  This is the
M4-close KNOWN DEAD ZONE (c) - "multi-sliver runs (total gap > eps) weld
strip-to-strip across more than eps" - recorded there as inherited and
unreachable in the suite.  On dense near-x-perpendicular real geometry (offset
walls, self-intersection sheets) it is PERVASIVELY reachable: the vast
majority of criticals collapse into sub-eps runs and many runs sum to more
than eps.  The canonical cap cannot bridge the gap - it fills only where the
two limits' regions DIFFER, but here the two flanking sections are the SAME
loop at two x-planes a few eps apart, so the cap is (correctly) empty.
CONTROL, AS CORRECTED BY THE VERIFICATION LANE.  The first-pass control - the
resolving Offset singles also have empty caps at essentially every critical -
refutes "empty caps == failure" but was CONFOUNDED as a discriminator, since
both populations share the empty-cap trait.  The real discriminator is run
width: the resolving controls have ZERO skipped runs wider than eps, while
each of the failing four has several.  Run width, not cap emptiness, separates
resolve from fail.

WALL A co-occurs as a MINORITY tail.  A small set of micro cap edges
(tangent-sheet ties, non-alternating fans, a few micro holes) is the
steep-track junction family.  The literal FIRST edge the splitter trips on is
NOT diagnostic of the dominant class - it may belong to EITHER class; on the
re-instrumented cases it was a wall-A-shaped micro fan (non-alternating), not
a macro hole - and even absent that first edge the surface carries macro
dead-zone-(c) holes throughout.  openscad is a richer variant of the same
class (part of its caps emit a small partial fill; its holes come in both
cap-present and strip-present flavors); its input imports as a VALID manifold
despite the "nonmanifold" filename, so the failure is ours, not the input's.
None of the four is WALL B (no engine id-conflict fires).

FIX DIRECTION (recorded, not implemented).  Two shapes: (1) BOUND THE RUN -
when a sub-eps run's total width exceeds eps, stop treating it as one gap; the
SubEpsFeature guard today catches only a single face living ENTIRELY inside a
run, missing content distributed across straddling faces, which is why these
report NonManifoldEmission rather than SubEpsFeature.  (2) WELD ACROSS THE RUN
by run-provenance - the far-side strip endpoints share (y, z) with the near
side (the partner search finds them exactly), so bind them by provenance rather
than by eps-distance.  Dead-zone (c) is the dominant blocker; the co-located
WALL A micro tail still needs the 3D-identity work afterward.

FIX-SHAPE ADJUDICATION (verification lane; recorded, not landed).  The leading
fix - place the strips' c-side chain verts at the pair-canonical critical's x,
so the cap and both adjacent strip boundaries share one exact plane and
closure becomes constructional rather than weld-dependent - is assessed SOUND
and BOUNDED (mechanism-sized, not research-grade), with one named scope
requirement: the far slab's strip then spans the skipped run, so a
run-width/SubEpsFeature guard (shape (1)'s bound) must own genuine macro
geometry changes inside the skipped interval.  This is the adjudicated next
mechanism, not landed work.

## CHAIN-PLANE RULE (main-agent, 2026-07-12)

NOTE (superseded in part by STRICT-FP SLAB BUILDING below): the wide-run guard
this section adds is UNREACHABLE under the strict-FP slab gate (runs shrink to
ulp scale), and the recorded PerpFaces fixture RESOLVES rather than fails closed.
The chain-plane rule's placement itself is UNCHANGED and remains load-bearing.
Read the STRICT-FP section for the current regime.

The dead-zone-c fix adjudicated above is implemented.  MECHANISM: a strip's
c-side chain emits at the PAIR-CANONICAL critical's x - the plane where its cap
bound it - not the strip's own slab boundary.  StripChains carries loX/hiX per
side; EmitCaps records crits[ci] at bind; ZipperEmit places each side there.  A
slab's hi side is always its own xHi (a slab is the left member of its pair), so
only the lo side of a slab following a skipped run moves: it spans back across
the run to the run's canonical critical, where the cap and the left slab's hi
edge already live.  Cap and both adjacent strip boundaries then share ONE exact
plane, so the shared cap loop closes CONSTRUCTIONALLY, not by an eps-weld the
run width can exceed.

GUARD ADJUDICATION.  The rule spans a skipped run by linear interpolation of the
two flanking sections.  That is EXACT only when they coincide (region(L) ==
region(R), empty cap - interpolating one loop reproduces it, at any width).  The
corpus dead-zone-c runs are exactly this (the same loop at two x-planes), so
they close with no guard.  But a run wider than eps whose flanks DIFFER
macroscopically (non-empty cap) is a real macro geometry change squeezed into a
sub-eps-per-slab interval that no linear span can carry - the slip-through the
single-face SubEpsFeature guard (BuildSlabs) misses, because its content is
distributed across straddling faces, not one face living wholly in the run.  ONE
guard closes it: ComputeCap fails closed SubEpsFeature when the cap is non-empty
over a run wider than eps (both flanks built; exterior end caps exempt).  The
run-width gate is load-bearing - without it the guard over-fires on DUST runs
(criticals sub-eps apart from vertex jitter) that carry legitimate box-end caps
and bridge eps-validly.  The recorded Coplanar_PerpFacesSubEpsApart fixture is
the constructed slip-through (offset squares across an ~eps run): before the
rule it fails closed by weld tear; with the rule alone it resolves oracle-WRONG
(an eps-displaced end face the winding grid samples); the guard restores a named
fail-closed.

RESIDUALS, honest.  The four single-mesh corpus fixtures no longer emit
dead-zone-c holes (instrumented to zero on Offset1 and self_intersectA/B, the
dominant class from the diagnosis).  Offset1 and self_intersectA/B still fail
closed at the sheet splitter on a co-occurring steep-track wall-A FAN (no open
boundary holes remain; out of scope, the 3D-identity research arc); openscad's
richer variant trips the new wide-run guard (SubEpsFeature).  GenericTwin7863
likewise trips it early (its pair gate now accepts that named guard alongside
NonManifoldEmission); Havocglass8 is unchanged.  The Cray pair and Offset2/3/4
keep resolving at unchanged volumes (zero runs wider than eps).  A red-first pin
(Pin_ChainPlaneRule_WideRunResolves: a ringed box whose long edges carry a
chained run wider than eps with coincident flanks) fails closed without the rule
and resolves oracle-true with it, mutation-verified against slab-boundary strip
placement.

## WALL-B: EngineIdConflict retired (main-agent, 2026-07-12)

NOTE (superseded in part by STRICT-FP SLAB BUILDING below): the residuals
recorded here moved under the strict-FP gate.  The hull no longer lands on the
chain-plane wide-run guard: its verified-real in-run content now SECTIONS, the
sections are dense enough to retain past the piece ceiling, and it refuses as
ArrangementBudget (Gate4c's skip is pinned to that detail).  The real-content
finding below stays accurate for the eps-gate regime it was measured under.
Read the STRICT-FP section for the current class map.

DECISION.  The EngineIdConflict fatal is RETIRED: BuildSlabs now counts a
section source-id conflict into cnt.engineIdConflicts (a benign diagnostic) and
CONTINUES instead of failing closed, and FatalReason::EngineIdConflict is
deleted (it had a single fire site).  This supersedes every earlier section that
treats EngineIdConflict as a live near-coplanar guard (Corpus measurement,
COPLANAR mechanism-5 skip narrowing, the touching-flap defense).

WHY IT IS SOUND - attribution is UNCONSUMED.  The conflict is raised by the 2D
engine when two distinct valid source ids coincide on one directed section edge
with nonzero net multiplicity: MergeSrcId sets srcId to -1 and reports the
event, but PolySetAdd sums the multiplicities and keeps the geometry REGARDLESS
(cancellation to net zero erases with no conflict; same-orientation content sums
and survives).  The 3D pipeline reads srcId NOWHERE: the SlabResolver keys and
matches pieces by POSITION, and winding retention is multiplicity-based.  So the
emitted boundary is byte-for-byte identical whether or not the conflict is
counted - the fatal guarded a label no consumer reads.  Grep-verified: the only
srcId/sourceId writer is the group-id seeding in BuildSlabs; there is no reader.

WHAT WOULD MAKE IT UNSOUND.  If a future change gives srcId a 3D CONSUMER (e.g.
per-source property transfer or provenance-attributed emission), a -1 conflicted
label would corrupt that consumer's input, and this demotion must be
re-adjudicated - the guard would again be load-bearing, and both the flap
re-adjudication and the residuals below would need review.

THE FLAP.  The touching-flap defense (audit round) rested on the hazardous
same-winding shared-edge flap "failing closed downstream via EngineIdConflict".
That backstop was already SUPERSEDED by the coplanar grouping pre-pass, which
unions coplanar faces INCLUDING shared-edge pairs into one group before any
exemption: a flap's two coincident faces then carry the SAME group source id, so
they never conflict - anti-oriented cancels, same-oriented sums to +2 and the
winding regularizes it to a single unit-boundary cover.  The demotion is
therefore ORTHOGONAL to the flap.  Pins: Coplanar_SameOriented_Oracle (two
same-oriented cubes composed into ONE input solid - a single self-overlapping
shell whose coplanar faces coincide same-winding over the overlap, resolving to
the union: the closed-shell folded-flap realization for the earlier [D2-4]
requirement) and Pin_TouchingDisjoint (the legal opposite-diagonal shared-edge
arm that cancels to net zero).  A bespoke doubled-coplanar-face fixture was tried
and REJECTED: a zero-volume doubled sheet is a degenerate input with no
well-defined winding oracle (it resolves empty), so it tests nothing.

RESIDUALS, honest - siblings to the wall-A arrangement arc:
- HULL (Gate4c).  With the fatal gone the hull runs past the former conflict and
  lands on the CHAIN-PLANE wide-run guard (SubEpsFeature) - the SAME honest
  boundary GenericTwin7863 trips.  Verified firing on REAL in-run content: the
  cap arrangement (run at 8*eps, so sub-8*eps construction noise has already
  annihilated) carries nonempty minus-side macro cap edges over a run wider than
  eps.  Gate4c's skip narrows from {EngineIdConflict, NonManifoldEmission} to
  {SubEpsFeature, NonManifoldEmission} and is pinned to the chain-plane guard's
  detail (any other SubEpsFeature must fail, not skip); red-first verified (the
  swapped skip reds while the fatal is armed, greens once it is demoted).
- GENERIC_TWIN_7081.  Its section conflicts are whole-segment near-degenerate
  junctions (three faces meeting on a near-point edge across several built-slab
  bands, correctly ungrouped - the dihedrals are macroscopic).  DIAGNOSIS
  (instrumented, correcting an earlier mislocation to "downstream emission"): the
  blowup is in the SLABS stage, not emission.  FindSeams emits an arrangement
  explosion on this near-degenerate geometry - critical x's numbering tens of
  times the face count out of only a few thousand seams - which dedups into tens
  of thousands of thin slabs.  BuildSlabs then both loops O(nSlabs*nFaces)
  straddle checks (the multi-minute cost) and RETAINS every built slab's section
  data (segments, edges, verts, pieces) simultaneously for the caps+strips
  stages, so cumulative retained pieces climb into the tens of millions and it
  dies by allocation roughly two-fifths through the slab loop; EmitCaps /
  EmitStrips / BuildImpl are NEVER reached.  HONESTY: the unbounded-retention +
  arrangement-explosion STRUCTURE is pre-existing (wall-A territory), but the
  EXPOSURE is demotion-created - the old fatal aborted this input in seconds, and
  removing it lets BuildSlabs run into the dense bands.  std::bad_alloc
  propagates as a clean catchable exception, but the swap-thrash before it is a
  de-facto hang, a regression against "never hang".  GUARD (spec [WALL-B]): the
  slabs stage now budgets cumulative retained pieces (FatalReason::
  ArrangementBudget, a refusal that the arrangement is too dense/degenerate to
  section within a sane resource budget - see kPieceBudget for the calibration).
  7081 fails closed fast (seconds, ~a gigabyte) instead of swap-thrashing, and is
  PINNED as a corpus fixture (Corpus_GenericTwin7081_Recorded) asserting that
  refusal.  The budget converts a de-facto hang into a recorded refusal; the real
  fix - arrangement robustness so the section does not explode - remains the
  wall-A arc.

Blast radius: the full Overlap3 suite stays green (Gate4c skips on the
chain-plane guard); the other former EngineIdConflict OR-list fixtures (Gate4d/e
/f, Coplanar_RazorBand, Coplanar_PerpFacesSubEpsApart) resolve or hit an accepted
named guard with EngineIdConflict dropped from their contracts; Boolean2 (2D
always seeds srcId 0, never conflicts) is unchanged.

## STRICT-FP SLAB BUILDING (main-agent, 2026-07-12)

DECISION.  The BuildSlabs slab-width gate lowers from the eps-width test
(xHi - xLo <= eps, which skipped whole runs of sub-eps slabs) to the strict-FP
floor: build every slab whose midpoint is strictly between its bounds in double
arithmetic (xLo < xMid < xHi).  An unbuilt slab is then only where the two
criticals are one ulp apart (their midpoint rounds onto a bound), so a section
plane never coincides with a critical.  The prior EngineIdConflict demotion
([WALL-B]) and the retained-section budget ([WALL-B]) are the two arcs that
unblocked this: the former near-coplanar attribution fatal (which the sub-eps
slabs surface pervasively) is now benign, and the budget bounds the denser
sections that building every slab produces.

WHY - the PRIZE is correctness, not a smaller gate.  The eps gate MERGED sub-eps
critical runs into single unbuilt gaps and refused their content (SubEpsFeature)
or collapsed it.  Strict-FP SECTIONS that content faithfully.  The money fixture
Coplanar_PerpFacesSubEps (two cubes ~0.5 eps apart in x, a macro a-square ->
b-square transition the eps gate merged) flips from a SubEpsFeature refusal to an
oracle-true RESOLVE - the principled cure (remove the CAUSE, build the slabs) as
opposed to the chain-plane arc's finding that merely DISABLING the wide-run guard
resolved the same fixture WRONG.  Fidelity gate held everywhere: zero oracle-wrong
resolves across the corpus and synthetic suite in either mode.  Cray resolves
bitwise-identically with 2-2.5x more built slabs; Offset2/3/4 shift by a few ulps
(deterministic, oracle-correct within the eps bound - the extra slabs re-associate
the arithmetic).  (Verification-round correction: the original "bitwise" claim
here and in the sub-eps probe came from a fixed-precision print; only Cray
survives an exact comparison.)

EVIDENCE (per-case class map, old eps gate -> strict-FP):
- Cray / Offset2 / Offset3 / Offset4: RESOLVE -> RESOLVE (Cray bitwise; the
  Offsets shift by ulps, oracle-correct).
- Coplanar_PerpFacesSubEps: SubEpsFeature -> RESOLVE oracle-true [STRENGTHEN].
- Pin_InPlaneSkeleton / Gate4a_Wedges8: RESOLVE -> RESOLVE (the probe-era
  EngineIdConflict breakage is gone with the demotion).
- GenericTwin7863 / openscad: SubEpsFeature (wide-run guard) -> NonManifoldEmission
  (the formerly-refused in-run content now processes to the wall-A sheet fan).
- hull: SubEpsFeature (wide-run guard) -> ArrangementBudget (its dense
  near-coplanar bands retain past the ceiling once every ulp slab builds).
- self_intersectA / self_intersectB: NonManifoldEmission -> ArrangementBudget
  (same density reason; FASTER - fails closed during slab build, not after a full
  section run).
- Offset1 / Havocglass8: NonManifoldEmission -> NonManifoldEmission (unchanged).
- GenericTwin7081: ArrangementBudget -> ArrangementBudget (unchanged).

WHAT DISSOLVES.  Nothing, and this CORRECTS the sub-eps probe's prediction that
the run machinery would delete.  Unbuilt runs do NOT vanish - an unbuilt slab is
exactly an adjacent-double critical pair, and runs CHAIN: hundreds of
non-canonical criticals on Offset2, a run of thousands of ulp slabs spanning
half an eps on degenerate corpus (GT7081), and a valid-manifold construction
pushes a run past eps.  (Verification-round correction: the earlier "widest run
well under eps" figure was measured only on cases that complete.)  So:
- CHAIN-PLANE RULE placement: STILL LOAD-BEARING.  Empirically deleting it (emit
  strips at slab bounds) regresses PerpFaces from RESOLVE to a fail-closed
  "unresolvable sheet contact" - a 1-ulp displacement between the post-run strip
  corner and the cap corner does NOT reliably weld into a paired sheet.  Root
  cause (verification round): the slab-bound placement truncates a
  constant-(y,z) sliver strip's x-width across the eps weld threshold - the weld
  collapses the sliver's verts, the degenerate triangles drop, and boundary
  holes open.  The chain-plane placement is the correct connect-to-the-cap-plane
  choice, not a workaround.  Kept; the probe assumed loX/hiX degenerate to slab
  bounds (a no-op), which is FALSE.
- WIDE-RUN GUARD: REACHABLE in precondition, budget-shadowed in fire
  (verification-round correction of "unreachable": a run wider than eps needs
  tens of thousands of consecutive adjacent-double criticals, and degenerate
  clusters DO produce them - observed at half an eps on real corpus, past eps by
  construction; but the same clustering trips the retained-piece budget or the
  flank cap's arrangement before the guard's critical is reached).  Kept as the
  chain-plane rule's fidelity BACKSTOP - deleting it, with the interpolating
  rule retained, would let a pathological wide run resolve silently oracle-wrong
  if the shadowing ever thins.
- NON-CANONICAL CAP SKIP (ci != li+1) and the single-face COVERAGE guard: both
  still exercised / reachable (multi-slab ulp runs; a near-x-perpendicular macro
  face landing wholly in a ulp run).  Kept.

WHAT REMAINS AND WHY (residuals, honest).  The former wide-run-guard cases
split: GT7863 and openscad now fail closed DEEPER - having processed their
in-run content, they hit the steep-track wall-A sheet fan (NonManifoldEmission).
The hull and the self-intersection pair fail closed EARLIER in the pipeline but
for an honest reason: their sections are now dense enough that the retained-piece
ceiling refuses during slab construction (ArrangementBudget).  Both are honest refusals in the same wall-A /
arrangement-robustness territory the earlier arcs named; strict-FP exposed what
the eps gate's early refusal masked, it did not create them.  The BUDGET is NOT
recalibrated up for the new ArrangementBudget cases: relaxing it would run them to
the emission blowup and STILL fail (NonManifold), strictly worse.  PERF: net
suite time is near parity, but that is a REDISTRIBUTION, not a free lunch -
clean cases build 2-2.5x more slabs at proportional O(nSlabs*nFaces) cost, dense
cases slow several-fold before their refusal (the hull roughly triples), and the
self-intersection pair's early budget refusal offsets both.  The sub-eps probe's
order-of-magnitude projection assumed those cases ran to completion; the budget
truncates them.  No single test approaches the suite bound.

Tests: Coplanar_PerpFacesSubEpsApart promoted to MUST-RESOLVE (red-first: it
fatals under the eps gate; mutation-verified as the chain-plane rule's pin -
reverting ZipperEmit to slab bounds reds it).  Pin_ChainPlaneRule_WideRunResolves
RETIRED: its RingedBox's eps-scale ring gaps all BUILD under strict-FP, so it can
no longer construct a wide interior run; the rule is now pinned by PerpFaces at
ulp scale.  Gate4c accepts ArrangementBudget (pinned to its detail).  The corpus
single-gate residual comments updated (self-intersection pair -> ArrangementBudget,
openscad -> NonManifoldEmission).

## WALL-A: the last correctness wall (anatomy, killed probes, the real fix)

The strict-FP class map leaves two fail-closed clusters: the sheet-fan
NonManifoldEmission cases (GT7863, openscad, Offset1, Havocglass8) and the
ArrangementBudget cases (hull, self_intersectA/B, GT7081).  This section records
the full three-part anatomy of that residual, every bounded shortcut probed and
killed with its evidence, the irreducible coupling the shortcuts all reduce to,
and the shape and preconditions of the real fix.  The prior adjudication (the
3D-IDENTITY EXTENSION section, the identity-extension notebook) named this
research-grade; the record below EARNS that verdict by killed probes rather than
assuming it, and refines the wall's anatomy under the strict-FP gate.

### ANATOMY, part 1: the sheet fan is the near-degenerate cluster, not a hole.

Under strict-FP the emitted surface of the four NonManifoldEmission cases carries
NO open boundary holes - the dead-zone-c macro holes the eps gate produced (the
old failclosed-diagnosis) are GONE, building the former sub-eps runs cured them.
The residual is FAN-INTERNAL: SplitTouchingSheets trips on radial TIES and
NON-ALTERNATING material, overwhelmingly on CAP-plane triangles (same-x), at a
small set of distinct junction clusters per case (order tens to hundreds by
eps-ball, not thousands).  The tie gaps are EXACTLY zero, not near the kAngleTie
threshold, and the non-tie fans sit far above it: the ties are two ANGULARLY
COINCIDENT sheets, not a marginal near-tangency.  The geometry: at a failing fan,
one 3D junction appears as a TWIN PAIR of emitted vertices that landed JUST over
the assembly weld eps and so did NOT merge; the caps of two ADJACENT slabs each
emitted the junction at its own extended position, producing a near-coincident
twin cap sheet (the non-alternating overlap) bounded by a micro-edge between the
twins (the gap-zero tie).  The twin separation is TRANSVERSE-DOMINATED, and its
axis is load-bearing for the kills below.  Measured on the failing Havocglass8 fan
(twice, independently, by both verification lanes): the two images sit SUB-eps in x
(order a hundredth of eps) and order a half eps to an eps-and-a-half apart in the
two transverse coordinates - a few eps in 3D, just over the weld.  This is NOT a
macro x-separation.  It is STEEP-TRACK amplification (order a hundredfold) of the
sub-eps x-offset between the two adjacent cap planes: TriangulateCap emits each cap
vertex at exactly its cap-plane x, so two criticals a sub-eps apart in x place two
cap images of the ONE junction at two distinct x's, and the steep extension track
magnifies that sub-eps x-gap into a super-eps transverse divergence.  The images
are kept distinct by the ASSEMBLY eps-weld - they live in SEPARATE cap arrangements
and meet only at assembly - NOT by Canonicalize.  This is the assembly-gate face of
the 3D-IDENTITY EXTENSION wall's near-degenerate cluster, in its case-(1) SPREAD
form (one arr.vert imaged twice), DISTINCT from the case-(2) DISTINCT pair (two
canonical verts kept apart by Canonicalize at a macro x-separation) that the
anisotropic-canonicalize kill below addresses.  Confusing the two axes mis-aims
that kill; the split is drawn explicitly there.

### ANATOMY, part 2: the budget cases are the same wall, plus over-inclusion.

The critical set is dominated by seam-seam CROSSING x's, not vertex endpoints -
crossings are the large majority of criticals on every case, resolving and failing
alike.  Crossings form dense BUNDLES (the 2D near-concurrence lifted to 3D: many
near-concurrent seams crossing pairwise pile their x's into a tiny band).  The
discriminator between resolve and fail is NOT bundle existence (the resolving
controls have bundles too) but bundle SPREAD: the resolving controls' crossing
bundles span well under eps, while every failing case has bundles wider than eps.
GT7081 is this at explosion scale - a single near-concurrent seam bundle piling on
the order of a million pairwise crossings (tens of times the face count) into a
band spanning a hundred-plus eps, from roughly a thousand near-concurrent seams;
the section cannot retain that many slabs and the retained-piece budget refuses
before emission.  hull and self_intersectA/B are NOT explosions - their crossing
ratio matches the resolving controls; they are merely LARGE, and the budget
refuses on cumulative retained pieces.  CROSS-ANATOMY verdict: wall A is ONE
degeneracy (near-concurrent-seam / near-degenerate-vertex clusters) on a SCALE
CONTINUUM.  hull, given enough budget to section, lands on the SAME cap-plane
sheet fan.  GT7081 is the fan degeneracy at a density that refuses before
emission.  self_intersectA/B are the far end: pure crossing over-inclusion with no
bad vertex twins underneath.

### KILLED PROBES (bounded shortcuts, each with its killing evidence).

- RADIAL-TIE RESOLUTION / kAngleTie tuning.  KILLED by the anatomy: the tie gaps
  are exactly zero (two coincident sheets = the twin over-emission), not a
  near-threshold spread, and the surviving fans are far above the threshold.  No
  maintained-order pairing exists for two coincident sheets - the fix is to remove
  the twin, not to pair it.  Tuning or replacing kAngleTie changes nothing.

- OUTPUT WELD-RADIUS BUMP / twin merge (BuildImpl).  The twins are only a few eps
  apart, tantalisingly close to the weld radius.  Bumping the assembly weld to
  merge them was probed to several times eps and does NOT resolve the sheet-fan
  case: merging the twin images trades the micro-edge for degenerate cap
  triangles that drop and reopen the fan elsewhere - the 3D-IDENTITY EXTENSION
  "cluster collapse creates vanishing pieces and flips topology", now confirmed at
  the OUTPUT (assembly) stage, not only at cap input.  A bounded output weld is
  the naked snap by another name.  KILLED.  SCOPE (verification round): the
  mechanism was measured on 2 of the 4 fan carriers - Havocglass8 directly (no
  resolve at 1.8x/2.5x/4x; cap triangles drop and the fan shrinks but never closes)
  and GT7863 via the identity-extension input-snap corroboration ("fixes the
  micro-edge, opens a macro hole elsewhere"); Offset1 and openscad were not
  separately weld-bumped, so the kill is a sample of the fan population, not a
  census.

- ANISOTROPIC CANONICALIZE (merge near-degenerate verts in the transverse plane,
  ignoring x).  KILLED by geometry - but for the case-(2) DISTINCT pair ONLY: that
  near-degenerate arr.verts pair is a real MACRO separation in x (order tens of eps)
  on a steep track - two distinct points, not one coincident junction.  Merging them
  in the transverse plane collapses real x-extent and moves a critical, which is
  exactly the topology flip the identity-extension snap hit.  SCOPE (verification
  round): this kill does NOT reach the case-(1) fan twin of ANATOMY part 1.  That
  twin is SUB-eps in x, not tens-of-eps, so a yz-merge there collapses no macro
  x-extent - the kill's premise is absent.  The fan twin's transverse merge is
  killed instead by the output weld-bump above (a yz-merge is a subset of the
  any-axis weld-bump that failed on Havocglass8) and by the cluster-collapse
  relocation in THE IRREDUCIBLE COUPLING below; the sub-eps-x cap-plane unification
  candidate below routes the same family through the cap-placement stage.

- PROVENANCE-IDENTITY TWIN MERGE (thread a junction label through emission, merge
  twins by label-equality).  The natural next idea after the provenance chains
  landed: if the two cap images "know" they are the SAME arr.vert, weld them by
  identity rather than by distance, sidestepping the >eps gap.  KILLED by
  RELOCATION, established by a code trace: the emission path is POSITIONS-ONLY
  end-to-end.  OutTri3D carries a bare vec3 triple (no id); the strip-chain
  edgeSubdiv is positions (vector<vector<vec2>>); the cap arrangement's
  RemoveOverlaps2D takes no edge-class/id argument (its public signature does not
  even accept one - class data lives on the internal winding pass); BuildSlabs'
  section call passes no class subdivision; and the assembly weld is an eps hash
  grid over vec3 (positions, radius eps).  So merging by identity requires threading
  a new label through slab -> cap -> edgeSubdiv -> OutTri3D -> weld, touching every
  emission structure.  And even granting the label, it changes only the SELECTION
  criterion (which images to merge, now robust to the >eps gap) - the merge ACTION
  still collapses a super-eps transverse divergence to one point, which doubles or
  degenerates the two cap sheets and reopens the fan exactly as the output weld-bump
  did.  This is the RSI-#3 relocation test: the shortcut moves the coupling into
  the emission plumbing without dissolving it.  Recorded as killed-by-relocation;
  the coordinated-re-emission fix below carries the identity already, by construction.

- BUNDLE COLLAPSE vs THINNING (a distinction, not a new lever).  Crossing-bundle
  thinning (below) removes bundle-INTERIOR criticals while keeping representatives
  at least eps apart, so no bounding vertex - hence no macro feature - moves or
  vanishes; it is over-inclusion pruning.  COLLAPSE-TO-ONE-CANONICAL of a wider-
  than-eps VERTEX cluster is a DIFFERENT shape: it MOVES a load-bearing vertex to a
  single canonical position, which is the 2D block-rule collapse lifted to 3D and
  the fan's actual root.  The two are easy to conflate ("just merge the cluster")
  but only thinning is bounded-and-safe; collapse is folded into the coordinated
  re-emission family in THE IRREDUCIBLE COUPLING, not into the thinning lever.

- SUB-EPS-X CAP-PLANE UNIFICATION (the correctly-aimed candidate for the case-(1)
  fan twin - PROBED, KILLED).  The twin's two images come from two cap planes only
  sub-eps apart in x (ANATOMY part 1), so the natural bounded shortcut is to let
  ONE cap plane own the junction: merge two criticals whose x-separation is sub-eps
  for cap-emission, emit the pair-canonical cap at one x, and let both flanks' strips
  bind there through the EXISTING chain-plane machinery (StripChains loX/hiX already
  place strip corners at cap planes).  Probed minimally and env-gated (demote a
  built slab narrower than a threshold-times-eps to unbuilt, so EmitCaps merges its
  two bounding criticals into one canonical cap and the chain-plane rule spans the
  sub-eps gap - composing with that rule, not rewriting it; thresholds a tenth-eps
  and one-eps).  RESULT: the fan does NOT close on ANY of the four fan carriers -
  each still fails NonManifoldEmission ("unresolvable sheet contact"), even at the
  threshold that demotes roughly half the slabs (the ~hundredth-eps twin slab among
  them).  KILLED, by two independently measured legs:
  (1) Cap-plane unification is INSUFFICIENT.  The merged cap runs ONE arrangement
  over both flanks at the construction-noise radius (well above the twin's few-eps
  spread), yet the twin survives.  Unifying the cap PLANE (a shared x) does not
  unify the junction POSITION: the two flanks still extend the junction
  independently and self-locate at different transverse points at that one x.  This
  is THE IRREDUCIBLE COUPLING below, now shown at the cap-PLACEMENT stage - a shared
  x is not the shared constructible point the fix needs.
  (2) The demoted BUILT slabs carry REAL content - the M4 disagreement, now with a
  BUILT slab.  The M4-close adjudication kept sub-eps critical pairs distinct
  because their arrangements disagree macroscopically about cancelled-edge
  subdivision; that was decided across UNBUILT gaps under the retired eps gate.
  Under strict-FP there is a BUILT slab between these criticals, and the
  disagreement STILL manifests: at the wider threshold a carrier trips the
  chain-plane wide-run FIDELITY BACKSTOP ("macro cap content over a skipped run
  wider than eps"), i.e. the flanks disagree macroscopically about the demoted run's
  content; and the resolving controls (the Offset trio) shift volume by ULPs when
  their sub-eps slabs are demoted - dropping that content MOVES the surface.  So the
  M4 adjudication EXTENDS to strict-FP and kills this family too: a built sub-eps
  slab is not over-inclusion to drop, it is real geometry.  This is exactly the
  BUNDLE-COLLAPSE-vs-THINNING line above - C2 thinning is resolve-preserving because
  it removes only CROSSINGS (over-inclusion); demoting built slabs moves geometry.
  The correctly-aimed fix remains the coordinated re-emission below (one shared
  junction vertex by construction), not any cap-placement merge.

### VALIDATED-SAFE, NOT LANDED: crossing-bundle thinning (the density sub-class).

Because crossings are the "over-inclusion is harmless" set (SEAMS: only their x is
consumed) and are the density's dominant source, thinning them is a bounded,
correctness-preserving lever - PROVIDED it touches only crossing x's (never
vertex endpoints, which bound every macro face) and keeps representatives at least
eps apart (so section resolution is never coarser than the fundamental tolerance,
and no macro feature can vanish - its bounding vertices survive).  Probed as a
sort-and-collapse of the crossing criticals to eps-spaced representatives:

- CORRECTNESS FENCE (strong): the full synthetic oracle suite stays green; the
  oracle pair (Cray) resolves ORACLE-BITWISE at every thinning tolerance - identical
  volume bits with and without thinning, verified by exact uint64 comparison; the
  resolving single meshes (the Offset trio) resolve to volumes identical TO ULPs
  with and without thinning.  Note the Offset trio is NOT uniformly bitwise:
  Offset3 is bitwise, but Offset2 shifts by ULPs (and drops one vertex) and Offset4
  by one ULP - all still oracle-correct, none moving a macro feature.  So on every
  case with a reference, thinning is resolve-preserving; the removed crossings
  genuinely contributed nothing (pure over-inclusion, as the design claims).  (A
  STANDING CAUTION recorded here because this is its second recurrence on the
  branch: a "bitwise" claim requires an EXACT bit comparison, never a fixed-precision
  print - fixed-precision output hid the Offset2/4 ULP shifts until an exact
  comparator was run, exactly as it did once before on the strict-FP arc.)
- EFFECT: self_intersectA and self_intersectB flip from the retained-piece budget
  refusal to a VALID resolve, TOL-INVARIANT (identical volume across a wide
  thinning range) and matching each other (two self-intersection meshes of the
  same object resolving to matching volumes) - strong corroboration the resolves
  are geometrically real, not thinning artifacts.  GT7081's crossings thin by
  nearly two orders of magnitude, but its ENDPOINT clusters (the near-degenerate
  arr.verts - the sheet-fan root itself) remain dense, so it still refuses.  hull
  sections once thinned and lands on its cap-plane sheet fan.

This is RECORDED, NOT LANDED, on the Voronoi-snap precedent (validated-safe,
documented, deferred).  Reasons: it flips NO oracle-bearing carrier - the only
flips are the single meshes, which have no a+b oracle (validity + tol-invariance
only), and asserting a resolve that cannot be oracle-checked violates the absolute
zero-oracle-wrong posture; and a core-stage change with recorded-contract churn
belongs to owner review.  (An earlier draft also cited hull turning a fast budget
refusal into a SLOWER emission-stage refusal.  The verification round did not
reproduce the slowdown - wall-time parity on repeated runs - so that leg is
dropped.  Thinning does still change hull's refusal from the retained-piece budget
to the cap-plane sheet fan, but it is a DIFFERENT fail-closed refusal at the same
cost, not a worse one.)  It is the density sub-class's landable mechanism for a
successor - once an oracle for the single meshes exists, or once it is combined
with the fix below that closes hull's fan.  (The perf-campaign journal already
named "criticals thinning" as a target; this is its correctness-fenced instance.)

### THE IRREDUCIBLE COUPLING, and the real fix.

Every bounded shortcut reduces to the same problem: ONE 3D junction is emitted
INDEPENDENTLY by two adjacent slabs' cap computations, at extended positions that
diverge past the weld radius; there is no single constructible point both sides
compute identically (welds and degenerate tracks self-locate nothing), and no
input snap or output weld can merge the images without collapsing load-bearing
pieces or flipping the cap arrangement's topology globally.  Closing it requires
the two caps to SHARE the junction vertex BY CONSTRUCTION - a cross-critical
identity - which means COORDINATED cap+strip re-emission around a near-degenerate
junction cluster: recognize the cluster as ONE entity, and rebuild its caps and
all incident strips together so both flanks consume one shared subdivision AND one
shared junction vertex.  This is the 2D law ("a dense near-concurrence collapses
to one shared vertex, and everything downstream re-emits through the maintained
order") lifted to 3D - and the same epistemic class as the June RSI-#3 arrangement
completion.  Its precondition, and the reason it is more than a local re-mesh, is
the spread-vs-distinct disambiguation the identity-extension arc named: a dense
cluster can carry MULTIPLE true vertices within a few hundred eps of construction
noise, so "which junction does this image belong to" must be decided coherently
for the whole cluster before re-emission, not per image.  The provenance chains
(landed) are the necessary subdivision foundation; the coordinated re-emission and
its cluster disambiguation are the research remainder, the companion of the
near-coplanar arrangement arc.  Corpus gates stay recorded contracts; the residual
is now anatomised as far as bounded mechanisms reach.

### REPRODUCTION NOTES (for a successor rebuilding the probes).

The probes above were env-gated and reverted; the notebooks record them by
mechanism and location, not exact diff.  Two details a reconstructor must recover
(a verification lane had to recover both by output-matching - recorded here so the
next one does not repeat the archaeology):

- The fan histogram (SplitTouchingSheets, whole-surface classification) counts TIE
  and NON-ALTERNATING as INDEPENDENT arms: a single edge can be both, so the two
  counters are not a partition and must be tallied separately - only then do the
  per-case histogram totals reproduce.
- The output weld-bump (BuildImpl getVertIdx) must scale the grid CELL along with
  the weld radius.  The weld is an eps hash grid with a 3x3x3 neighbor sweep; if
  only the radius grows while the cell stays at eps, the sweep misses candidates
  now inside the enlarged radius and the bump under-merges.  Scale both together.
