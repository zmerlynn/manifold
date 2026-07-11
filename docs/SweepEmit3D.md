# Sweep-native emission for 3D overlap removal (design, v3)

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
contributes no strips. [R3-fold, AMENDED at the emission-closure
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
