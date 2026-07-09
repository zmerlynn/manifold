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

A - MERGE (unchanged from SweepPlane3D.md stage A): eps vert merge,
identical-tri multiplicity, zero-mult drop.

B' - CRITICALS AND SEAMS, SOFTENED. The critical set = sorted x of:
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

C' - SECTIONS (unchanged mechanics): per slab wider than eps, at
mid-x every straddling face contributes one directed segment
(sweepDir x outwardNormal convention, signed multiplicity), into
the 2D engine's arrangement + winding passes. The engine returns
the RETAINED boundary pieces (its native output - the piece capture
with (below, above) is no longer needed for emission; retained
pieces suffice) with source ids.

D' - STRIPS. Every retained piece endpoint tracks a 3D segment:
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
the strips emit). Strips are emitted as quads
(two tris); orientation is ENGINE-NATIVE: the winding pass emits
retained pieces interior-on-left in section space, which together
with the sweep direction determines the material side of the strip
- no per-face re-derivation (transfer-in-disguise is forbidden).

E' - CAPS. At each critical x = c, the retained REGIONS of the two
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
lift ambiguity - subdivide the cap loops AND replace each adjacent
strip's c-side interval with the arrangement's subdivision of that
interval. A constructed subdivision vert lies within alpha of its
host segment (Smith 8.2, the shared kernels), hence within alpha of
the strip's source-face plane - the standard eps-validity, not a
new error class. Closure is then by shared construction, not
assertion. Both cap_plus and cap_minus are computed (one engine
Subtract each); cap regions triangulate via Triangulate (holes CW
per its contract - never fans); the exterior limit beyond the
first/last critical is the empty region, making the outermost caps
ordinary.

DEGENERATE SLABS (width <= eps): no section is built; the slab
contributes no strips; the caps at its two bounding criticals are
computed from the NEAREST BUILT slabs on each side extended to the
respective critical - equivalently, consecutive criticals with no
built slab between them merge into one cap plane evaluated once
(at the first critical of the run, all of them within eps of each
other; positions within eps are one eps-valid plane). One rule, no
anchor propagation, no span guard, no starvation class: a dense
cluster of criticals = one merged cap between its flanking built
slabs. [R1-fold - the dissolution claim was WRONG, convergent
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
this design.
BORN: track extension (Interpolate at two x's per piece endpoint),
cap construction (one engine Subtract per critical), strip/cap
assembly. Expected net: a large deletion.

## Eps posture

Unchanged: one absolute eps (EpsilonFromScale budget 1000), input
quantization at stage A, slab-width gate, the merged-cap rule for
sub-eps critical runs. The engine inside slabs and inside cap
booleans runs its standard posture. Alpha covers every constructed
point via the shared kernels.

## Output and tests

Output tessellation is per-slab strips + per-critical caps - finer
than input faces (geometrically ON input faces for strips);
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
