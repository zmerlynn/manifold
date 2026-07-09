# 3D overlap removal by sweep plane (prototype design)

STATUS: DESIGN COMPLETE - crucible closed at round 3 (round 1: 5
lanes, 4 convergent BREAKs resolved [R1-fold]; round 2: re-attack +
fresh-eyes + empirical id-plumbing build, round-1 fixes held,
narrower findings folded [R2-fold]; round 3: convergence check, one
final BREAK - the chained-degenerate-span hole - resolved by the
span guard, plus the implementability closures [R3-fold]). The
Crucible record at the bottom is the review trail; residual risks
are stated there honestly. The handoff
(.claude/plans/3d-sweep-plane-prototype-handoff.md) and its
reconciliation sibling are the framing. Simple-but-robust is the
bar; performance is a later campaign. The prototype exists to answer
ONE question: is a sweep status that is itself a moving 2D
arrangement tractable?

## Contract

Input: a topology-manifold triangle mesh whose triangle interiors may
cross (a self-overlapping solid; e.g. Compose of overlapping parts,
or a Boolean output with residual self-intersection). Output: a
manifold mesh that is the boundary of the input's winding > 0 region,
eps-valid in Smith's sense (7.7: topological validity, geometric
positions within the error budget), or a REPORTED failure - never a
silently wrong mesh. Fail-closed classes are named per stage below;
the prototype measures them. Oracle: Boolean3 (`Manifold::Boolean`).

## The 2D lesson, carried up (settled - do not re-litigate)

- Local, order-free decisions do not compose on dense
  near-concurrences. Every classification derives from a maintained
  sweep order or a settled arrangement built by one; never from
  independent per-item probes deciding arrangement or classification.
- Robustness lever is tolerance + consistency: shared
  `Interpolate`/`Intersect` kernels, `EpsilonFromScale` budgets, no
  exact predicates, no double-double.
- Coincident-event handling is the block rule (Smith 7.6.2), not
  symbolic perturbation.
- The 2D engine (`src/boolean2_sweep.cpp`, `SweepWinding`) is the
  validated consistency solver: input features eps-quantized once,
  the sweep tolerance-free with exact-coordinate identity, dense
  clusters collapsed through event vertices by the block rule,
  winding decided by the maintained order. EMPIRICALLY CONFIRMED
  (round-1 probe): a hand-built 3D cross-section fed to SweepWinding
  yields correct union winding; the piece-capture extension is ~15
  implementation lines with zero behavioral change (111 existing 2D
  tests green over the patched engine).

## Design overview

The kinetic insight: in 2D two edges cross at a POINT (one event); in
3D two faces cross along a SEGMENT that persists across an x-interval
of the sweep. The 3D status is a 2D arrangement with persistent,
continuously-deforming crossings whose TOPOLOGY changes only at
discrete x-criticals. Between consecutive criticals (a "slab") the
section's combinatorial structure is constant, so per-slab recompute
with one representative section per slab is sound. Held under attack
(round 1): for supported non-coplanar pairs there is no slab-interior
topology change - section combinatorics change at vertex x's, seam
endpoints are input verts or edge-face events, seam-order swaps on a
face are triple points, and a straight segment's x is monotone.

Stage A - merge and canonicalize. Merge input verts within eps
(deterministic union-find, centroid-nearest representative - the 2D
engine's convention). Drop collapsed tris; merge exactly-identical
tris (same canonical vert triple) into one record with signed
multiplicity (orientation-relative); zero-mult records drop.
[R1-fold] Multiplicity is a WINDING input only - it never rides into
emission (stage E emits exactly one oriented sheet per fill
transition, exactly as the 2D engine's EmitBoundary saturates to
+/-1). Its role: Compose(A, B) duplicates interface triangles;
cancellation and winding accumulation need the signed count, and
gate 5 is unpassable without it.

Stage B - arrangement geometry (seams and their verts), computed
once, globally, as shared objects:

- Edge-face events: for each (edge, tri) candidate from a BVH broad
  phase, the eps-aware piercing test (shared kernels); each event
  vert snaps to an existing vert within eps (the incidence
  quantization, applied once here).
- Seams: for each intersecting tri pair, the clipped intersection
  segment. [R1-fold] The naive "exactly two endpoint verts"
  invariant is FALSE for degenerate contacts; the intersection
  classes are enumerated with a rule each:
  - proper segment, length > eps: a seam with exactly two distinct
    canonical endpoint verts (asserted only AFTER the other classes
    are peeled off);
  - point contact / tangent touch / segment of length <= eps: no
    seam FROM THE PAIR ALONE. [R2-fold] A macroscopic intersection
    curve chopped by fine triangulation into per-pair sub-eps clips
    must not vanish clip by clip (and stage-B event merging eats
    such a chain rather than saving it): skipped contacts cluster
    into connected components first. [R3-fold, the graph defined]
    Nodes are the skipped contact primitives' canonicalized endpoint
    points (a point contact contributes one node, a sub-eps segment
    its two endpoints); edges connect nodes of the same primitive
    and nodes within eps (3D Euclidean); component diameter = max
    pairwise 3D distance over member points. Diameter <= eps: drop
    and count (genuinely point-like). Larger: the SUB-RESOLUTION
    SEAM CHAIN class, fatal (`SubResolutionChain`). A corner-touch
    of two boxes clusters to one point-like component and drops;
  - a tri edge lying in the other tri's plane (edge-in-plane, the
    seam-collinear-with-input-edge case): OUT OF SCOPE with the
    coplanar class (below), detected and reported.
- Triple points: for each face, pairwise crossings of its incident
  seam segments IN THE FACE PLANE (2D crossing test, shared kernel).
  [R1-fold, the round's central convergent BREAK] The same physical
  triple point is computed independently in up to three face frames
  and lands on different seams as different candidates, ~2*alpha
  apart (Smith 8.2: alpha ~ 12.37*u*L per construction; 2*alpha is
  ~eps/500 at budget 1000) - so per-seam local merging CANNOT
  deliver shared identity, and eps chain-merges along a seam are
  order-dependent (the #289 step-1 chain problem reborn). The fix is
  GLOBAL TRIPLE-POINT UNIFICATION, the production lesson of RSI's
  step 9.5 (same point computed through different frames must be
  unified globally): after all candidates are computed, one
  union-find over all triple-point candidates and seam verts within
  eps of each other, with a DIAMETER GUARD - a cluster whose point
  spread exceeds eps is the genuinely-degenerate chain class and
  FAILS CLOSED (reported, run returns failure) rather than
  chain-merging; each surviving cluster takes one canonical vert
  (centroid-nearest convention), written back into EVERY incident
  seam polyline. Sharing is thereby by-construction again: both
  faces of a seam reference the same polyline object; every face at
  a triple point references the same vert id. Well-conditioned
  crossings sit ~2*alpha apart and unify trivially;
  ill-conditioned (near-coplanar) crossings whose error exceeds eps
  hit the diameter guard and fail closed - the RSI conditioned-band
  class, measured not solved.
- On-edge incidences: events and seam verts within eps of an input
  edge subdivide that edge (ordered by t - 1D, sound), so adjacent
  faces conform across shared input edges by shared vert ids.

Stage C - the sweep: classification by slabs. x-criticals = sorted
union of merged-vert x's, event-vert x's, and triple-point x's,
BRACKETED by two exterior sentinel slabs (before min-x and after
max-x: empty sections, winding 0 everywhere - [R2-fold] the
boundary caps of any fixture need an exterior side). For each slab
wider than eps, at mid-x build the SECTION: every face straddling
mid-x contributes one segment (edge crossings at mid-x via
`Interpolate`), with its face id and signed multiplicity, seeded
into the 2D engine (arrangement + winding passes) in (y, z)
coordinates. The engine returns every arrangement piece with its
source id and (below, above) winding (the engine extension, below).
Because the input surface is closed and even-manifold after stage A,
section winding equals 3D solid winding along that plane (held under
attack round 1, including inverted and nested shells and
multiplicity > 1), so per-slab 2D winding IS the 3D classification,
decided by one maintained order per slab.

ORIENTATION AND SIGN CONVENTIONS (empirically validated round 1):
- Sweep axis x; section coordinates (u, v) = (y, z).
- A face's section segment is directed CCW AS VIEWED FROM +x
  (equivalently: direction = normalize(cross(sweepDir,
  outwardFaceNormal)) up to the segment's own parametrization), with
  multiplicity = the stage-A signed multiplicity. The engine's lex
  normalization handles reversal (m negates).
- (below, above) are STATUS-ORDER sides, not spatial sides: for
  non-vertical pieces below = smaller-z side; for VERTICAL pieces
  (z-parallel in section) the labels follow gradient-rank status
  order and are spatially inverted relative to z intuition. The
  emission decision IsInside(below) != IsInside(above) is
  orientation-safe for both; any consumer needing spatial sides must
  use the non-vertical rule only.
- Point-winding query (used for axis-parallel faces only): winding
  at (qy, qz) = sum over NON-VERTICAL captured pieces crossing the
  downward ray y = qy, z > qz of (below - above). Vertical pieces
  are skipped (parallel to the ray). Validated round 1 on interior,
  overlap, and exterior query points; the naive crossing-direction
  ray count is WRONG (multiplicity), which is why the algorithm is
  spelled here.

AXIS-PARALLEL FACES: a face lying in a section plane x = c never
straddles a slab. Classification: w(just below c) vs w(just above c)
at a clearance-checked interior probe point (largest-triangle
centroid of the PSLG region), via the point-winding query against
the two adjacent slabs' captured arrangements (the exterior sentinel
slab serves min/max-x caps; an adjacent interior slab that was
skipped as sub-eps sends the face to the degenerate/starvation
handling [R2-fold]). Querying a settled arrangement built by the
maintained order is evaluation, not the refuted independent-probe
pattern; clearance failure falls into the degenerate handling. Alternative considered and rejected (round 1):
a tilted "irrational" sweep direction shrinks but cannot eliminate
the class (any FP direction is rational; adversarial normals still
hit it), adds rotation FP noise to every coordinate in and out, and
abandons input-coordinate exactness - while the query is the direct
3D analog of the engine's own MergeVerticals1D special-casing.
Cost stated honestly: a per-slab point-location helper over captured
pieces (~40 lines), not a one-liner.

Stage D - partition and classification. [R3-fold ordering, one
place] Stage D opens with PSLG VALIDATION for every face - seams
crossing anywhere but shared ids is a stage-B miss and fatal
(`PSLGInvalid`), before any walk consumes the data. Then the
face-local planar walk (purely topological): each face's PSLG = its
subdivided boundary + its seam polylines; regions are the walk's
output.
[R1-fold] Interior seam loops are real (a seam both of whose
endpoints are events of the OTHER face's edges floats in this face's
interior): PSLG components not attached to the outer boundary are
classified by signed area and fed to `Triangulate` as holes/islands
through its existing keyhole machinery (nested loops included).
Then classification. [R2-fold, the correspondence spelled out] A
region's covering slabs are those whose OPEN x-interval lies inside
the region's projected x-extent; take the widest. In that slab, the
region's pieces are the captured pieces carrying THIS face's id
whose midpoints locate inside the region (2D point-in-region on the
face plane, lifting the piece midpoint to 3D along the face's
section segment), with a clearance guard: a midpoint within eps of
the region boundary is skipped as ambiguous. The region requires at
least one located piece, and ALL its located pieces must agree on
oriented (below, above); zero pieces or disagreement fails closed
(named: classification ambiguity). One face's section segment
routinely carries pieces of several regions (each seam crossing at
that x splits it) - the point-location owns the assignment; no
t-interval arithmetic against constructed values.

EMISSION ORIENTATION, algebraically (vertical-safe) [R2-fold]:
crossing a face along its outward normal changes winding by
-m_face, so w_front = w_back - m_face. A captured piece gives
(below, above) with above = below + m_lex, where m_lex = +/-m_face
by the piece's lex direction relative to the face's seeded segment
direction; the sign resolves back/front from below/above without
any spatial-side reasoning (which fails for vertical pieces). Keep
iff IsInside(w_back) != IsInside(w_front); output normal = the face
normal if IsInside(w_back) (material behind), else flipped. One
unit pin verifies the algebra on a cube's six faces, including the
normal = +/-y (all-vertical-pieces) sides.

[R1-fold] DEGENERATE REGIONS (no covering slab wider than eps) are
classified by ANCHOR-COMPONENT PROPAGATION, not free inheritance:
build connected components of degenerate regions (adjacency =
shared PSLG edges); a component's classification propagates only
from its NONDEGENERATE anchor neighbors; [R2-fold] ALL anchors must
agree on the ORIENTED (below, above) classification (bare keep/drop
agreement can mask an orientation conflict), else - or with no
anchor at all (slab starvation: a whole component whose criticals
all pack within eps) - the component FAILS CLOSED. Dense tangles
with agreeing anchors classify; conflicting evidence is never
averaged or tie-broken.

[R3-fold, replacing round 2's adjudication - its premise was
punctured on re-attack] Round 2 argued a dropped degenerate
component is "by construction sub-eps-thick along x"; FALSE for
CHAINS: connected degenerate regions, each individually without a
wide covering slab, can span macroscopic x. The rule is therefore a
SPAN GUARD: anchor propagation and the eps-feature drop apply only
to components whose 3D diameter is <= eps - for those, dropping is
the documented eps-scale decision (#289 step 1's "smallest feature
size that is desired to retain"), applied consistently (balance
below) and counted. A degenerate component with diameter > eps is
UNCLASSIFIABLE and fatal (`UnclassifiableComponent`) - never
propagated into, never silently dropped. How dense an input can get
before this guard fires is precisely the tractability data the
prototype exists to produce.

[R1-fold] SEAM BALANCE is enforced BEFORE emission, not assumed
after: for every seam polyline edge, the kept incident regions must
pair (equal start/end counts - the 2023 triplet constraint as a
pre-emission check). A violation names the seam and fails closed.
Manifoldness of the output is then the gate-3 CHECK; the design no
longer claims it as a theorem, it engineers toward it and verifies.

Stage E - triangulation and emission. Each kept region (fill
transition per the orientation algebra above) triangulates via
`Triangulate` projected to the face plane; exactly ONE sheet per
transition regardless of multiplicity. (PSLG validity was already
enforced at stage D's start; with global unification in stage B,
the round-1 false-assert class - same point, two ids - is resolved
at its root.) Projection and winding for the Triangulate handoff
[R3-fold]: each kept region is emitted as flat `PolygonsIdx` in
`GetAxisAlignedProjection(outputNormal)` - outer contours CCW,
holes CW, nested islands CCW - so a flipped output normal flips the
projection rather than post-reversing triangles; returned triangles
map back through `PolyVert.idx`.

FAILURE CONTRACT, one shape [R2/R3-fold]: every stage returns a
StageResult carrying either its product or a fatal reason code:
`TripleDiameter`, `SubResolutionChain`, `CoplanarOverlap`,
`EdgeInPlane`, `ClassificationAmbiguity`, `AnchorConflict`,
`UnclassifiableComponent`, `Starvation`, `BalanceViolation`,
`PSLGInvalid`, `EngineIdConflict`. Fatal stops the pipeline at that
stage; no later stage consumes invalidated data. Non-fatal COUNTERS
(sub-eps point contacts dropped, degenerate components classified
by propagation, eps-feature drops, clearance skips) accumulate in a
report struct returned alongside success - visible, never fatal.
Clearance failures map to `ClassificationAmbiguity` when slabs
exist and to `Starvation` when none wider than eps does. The
engine's own capture-conflict counter is nonfatal INSIDE the 2D
engine (sourceId = -1 + count; 2D callers are unaffected) and any
nonzero count is fatal to the 3D pipeline (`EngineIdConflict`) -
one rule, both empirically exercised.

## The engine extension (single change to boolean2_sweep, priced)

Round-1 empirical probe: geometric piece-to-segment matching HARD
FAILS at sub-eps segment separation (misattribution with no
resolving threshold, silently wrong classification on exactly the
gate-4 inputs). So the source id threads THROUGH the engine:

- `PolySet2`'s mapped value grows from `int64_t` to {multiplicity,
  sourceId}; `SweepEdge` carries the id; splits preserve it;
  `EmitBoundary` emits it. Exact-coincident segment summing (the
  only id-conflict site) cannot occur between distinct faces at
  mid-slab: identical section segments require the two faces to
  share both crossed edges (impossible for distinct tris), and
  stage-A merging gives identical tris one record; conflict asserts.
  Collinear-overlapping verticals from distinct faces would need an
  input edge lying exactly in the mid-slab plane, which mid-slab
  placement excludes (verts are criticals); asserted likewise.
- A default-null out-channel on `SweepWinding` records (from, to,
  sourceId, below, above) for every winding-pass piece, threaded
  through `CollectThenMeasure` (a header signature addition with a
  default, backward-compatible; existing callers unchanged -
  verified against the 2D suite in the round-1 probe).
- The per-slab point-location helper (the query above) lives with
  the 3D code, not in the engine.

Honest price [R2-fold, measured by building it]: ~174 functional
lines across `boolean2.h` + `boolean2_sweep.cpp` (+86 net file
lines). The round-2 empirical lane implemented the full plumbing:
111/111 existing 2D tests green, id fidelity perfect through splits
and block-rule re-entry, the round-1 killer (near-parallel segments
1e-9 apart) attributed correctly by ids where geometric matching is
provably ambiguous, and the conflict counter safe on
impossible-in-production coincident-id input. Two items the draft
under-priced, both confirmed independently by the fresh-eyes lane:
`MergeVerticals1D` needs an algorithmic rewrite (running-coverage
delta-sweep -> per-interval active-contributor tracking; adjacent
same-plane triangles make multi-source vertical groups ORDINARY,
not a corner case - each emitted interval carries the id of its
unique active contributor, multi-source overlapping intervals
count a conflict and carry -1), and the `pending_` inner map's
value type changes alongside PolySet2's. In exchange the 3D side
DELETES the geometric matcher and its failure analysis entirely.

## Types, constants, and metrics [R3-fold: the implementability closures]

- `CanonicalFace { int id; ivec3 verts; /* stored orientation */
  vec3 normal; /* from that order */ int64_t mult; /* signed,
  relative to that order */ }` - stage A's output; every later
  "face id" means this record. `mult == 0` records are dropped in
  stage A, so `m_face != 0` always.
- Edge-face event allocation: an event point unifies with an
  existing canonical vert within eps (3D Euclidean), else appends a
  new event vert. Event verts and triple vert x's feed the critical
  set.
- EPS METRIC TABLE - every bare "within eps" resolves as: 3D
  Euclidean point distance for vertex merge, event snapping, triple
  unification and cluster/component diameters; 3D point-to-segment
  distance for on-edge incidence; 2D face-plane Euclidean for
  region point-location clearance; scalar `xHi - xLo > eps` for
  slab width. One absolute eps everywhere (below).
- `SweepCapture { vec2 from, to; /* section (y,z), lex-forward
  measure-pass direction */ int32_t sourceId; int64_t below,
  above; }`; `SweepWinding(..., std::vector<SweepCapture>* capture
  = nullptr)` - default-null, existing callers unchanged.
- `SectionFaceSegment { int faceId; vec2 p0, p1; /* directed:
  dot(p1-p0, yz(cross(+x, face.normal))) > 0 */ int64_t mult; }` -
  stage C persists one per straddling face per slab; the emission
  algebra's m_lex sign compares a piece's lex direction against
  THIS direction.
- `SlabResult { double xLo, xHi, xMid; bool built; /* false: skipped
  sub-eps */ std::vector<SweepCapture> pieces; }` per slab, plus
  the two exterior sentinels (built, empty pieces, winding 0). The
  point-winding query and region classification read these.
- Region x-extent = [min, max] of the region's boundary verts'
  x's (post-unification); covering slab = open (xLo, xHi) inside
  that interval; ties between equal-width slabs break to lowest
  xLo. Axis-parallel adjacency at critical c = exactly the two
  intervals (prevCritical, c) and (c, nextCritical); a sub-eps
  non-sentinel neighbor routes to degenerate handling, never a
  farther search.
- Gate-2 "just past the critical" sample = the midpoint of the
  adjacent non-skipped interval (no nextafter games).

## Epsilon posture

One absolute eps = `EpsilonFromScale(bbox scale, 1000)` unless the
caller passes one (identical to 2D). It quantizes input features
(stage A, stage B snapping/unification/on-edge subdivision) and
gates slab width; the 2D engine inside each slab runs its own
standard posture. Smith's alpha covers every constructed point
(edge-face verts, triple points, section endpoints - all via shared
kernels); the global unification's diameter guard is where the
budget's limits surface as reported failures instead of corruption.
The >= 2^40 coordinate degeneracy is accepted manifold-wide.

## Deliberately out of scope (prototype; each detected + reported)

- Coplanar face overlap beyond exactly-identical triangles,
  INCLUDING the edge-in-plane seam class and the near-coplanar
  conditioned-band class (ill-conditioned triple crossings whose
  frame error exceeds eps - they hit the stage-B diameter guard).
  Reference: the RSI branch's OverlapRemoval.md known limitations 1
  and the 6.5 trace machinery, if ever imported.
- Vertex-on-face through-pierces (an edge ENDPOINT exactly on
  another face's interior plane produces no edge-face event; the RSI
  known-limitation-10 class). NOT detected [R3-fold wording]: there
  is no stage-B detector for it; fixtures avoid the class and an
  occurrence surfaces downstream as a `BalanceViolation` (correct
  fail-closed outcome, imprecise blame - accepted for the
  prototype).
- Slab-starved components (thin-in-sweep-axis geometry, all
  criticals within eps): fail closed, counted.
- Cancellation/ExecutionContext, progress, parallelism, performance.
- Public API: the seam is internal (`RemoveOverlaps3D(const
  Manifold::Impl&, double eps) -> std::optional<Impl>`, failure =
  nullopt + a reason report struct for tests); public wiring is
  post-prototype.
- Output simplification/decimation (Simplify's job, as in 2D).

## Validation gates -> tests (test/overlap3_test.cpp)

1. EVENT PARITY: brute-force O(n^2) edge-face and seam-pair
   enumeration equals stage B on small fixtures (two tets, two
   boxes, box + rotated box).
2. SECTION VALIDITY: per slab: closed section (even vert degree,
   signed multiplicity sums zero around each vert); engine ingests
   without assert; ZERO RESIDUAL CROSSINGS after the arrangement
   pass (piece endpoints only at shared vertices - the
   non-overlapping-status property, measured); every piece carries a
   source id. Sampled at mid-slab AND immediately past each
   x-critical on the fixtures (a wrong critical set manifests just
   past the transition, not at mid-slab).
3. TRIPLET PAIRING / MANIFOLDNESS: on fixtures with genuine
   triple-face verts (three boxes pairwise overlapping around a
   common region; three plates through one line region), seam
   balance holds pre-emission and the output passes the manifold
   gate (paired edges, balanced verts) + Manifold(Impl)
   construction.
4. DENSE NEAR-CONCURRENCE (adversarial) [R2-fold: non-vacuous;
   R3-fold: fixtures enumerated]: the fixture list is SPLIT.
   MUST-RESOLVE stress fixtures - constructors with eps-relative
   numbers, eps ~ 2.75e-9 at unit scale:
   (a) kWedges(k=8): thin boxes 1 x 0.02 x 0.3, rotated k ways
       about z through a common region, axis offsets ~1e-3 (>> eps,
       criticals separated well above eps);
   (b) nearParallel(sep=1e-6): two unit plates 1e-6 apart (~350
       eps) crossed by a third at a shallow angle;
   (c) the June trimaran hulls (hull-body.obj - hull-mask.obj,
       imported OBJ fixtures) run through Compose + the seam.
   PASS requires manifold output + oracle agreement - a fail-closed
   here FAILS the gate (no vacuous safety).
   MUST-FAIL-CLOSED degeneracy fixtures:
   (d) nearParallel(sep=1e-10) - inside eps: `TripleDiameter` or
       `UnclassifiableComponent`;
   (e) a strip triangulated so per-pair clips are ~0.75 eps:
       `SubResolutionChain`;
   (f) kWedges with axis offsets ~0.3 eps: `TripleDiameter` or
       `UnclassifiableComponent`.
   PASS requires the NAMED guard - a resolved-but-wrong mesh or an
   unnamed crash fails. The guard-firing rates across (a)-(f) are
   the tractability data the prototype exists to produce.
5. ORACLE: for two-operand fixtures, RemoveOverlaps3D(Compose(A, B))
   vs Boolean3 A+B: |volume difference| <= eps *
   max(surfaceArea(ours), surfaceArea(oracle)); genus equal;
   Contains() agreement on a deterministic grid [R3-fold: spec] -
   the union bbox inflated by 5%, 17x17x17 lattice, fixed traversal
   order, points within eps of EITHER surface skipped (counted);
   all remaining points must agree.

Plus: the existing manifold suite stays green (the branch adds files
plus the engine value-plumbing; the 2D suites are the regression
fence for that plumbing).

## File plan

- `src/overlap3.h` - internal seam + stage structs (~120 lines).
- `src/overlap3.cpp` - stages A, B, D, E (~1200-1500 lines; the
  round-0 estimate was optimistic per audit - stage B's enumeration
  + unification and stage D's propagation/balance are the bulk).
- `src/overlap3_sweep.cpp` - stage C: slabs, sections, engine calls,
  point-location helper (~300 lines).
- `src/boolean2.h` + `src/boolean2_sweep.cpp` - id/value plumbing +
  capture out-channel + the MergeVerticals1D contributor-tracking
  rewrite (~175 functional lines, measured by the round-2 build).
- `test/overlap3_test.cpp` - the gate ladder (~700 lines).
- CMake: sources added in `src/CMakeLists.txt`; the test file added
  to the SOURCE_FILES list in `test/CMakeLists.txt`.

## Crucible record

Round 1 (5 lanes): 4 convergent BREAKs (triple-point identity;
inheritance vs seam balance; seam-endpoint invariant; slab
starvation) - all resolved and HELD under round-2 re-attack. The
empirical lane validated the core bet (section -> engine -> correct
winding) and the capture extension against the live 2D suite.
Round 2 (re-attack + fresh-eyes + empirical build): the round-1
fixes held; new narrower findings folded - sub-resolution seam
chains fail closed; the piece->region correspondence specified
(midpoint point-location, clearance, agreement); oriented anchor
agreement; eps-feature-drop semantics adjudicated as documented
behavior, not silent wrongness; exterior sentinel slabs; the
emission-orientation algebra; one failure contract with PSLG
validity ordered before the walk; gate-4 non-vacuity split; the
engine extension re-priced from the round-2 build (~174 functional
lines; MergeVerticals1D rewrite; 111/111 2D regressions green).

Round 3 (convergence): the round-2 re-attacks HELD (correspondence
spec computable; orientation algebra verified on +/-y, inverted, and
mult-2 faces; gates 3/5 walk through without spurious guards). One
final BREAK - the chained-degenerate-span hole - punctured round 2's
eps-thickness adjudication and is resolved by the span guard. The
implementability lane returned 18 definitional closures (no design
objections), folded as the Types section, the fatal-code
completions, the fixture/oracle enumerations, and the ordering fix.

## Residual risks (stated, measured by the gates - not resolved)

- The guard-firing rate on MUST-RESOLVE fixtures is a prediction,
  not a proof: if kWedges(8) at 1e-3 offsets trips the span guard,
  the design's supported-density envelope is smaller than believed
  (that answer, either way, is the prototype's deliverable).
- The clearance rules can starve thin-but-real regions of usable
  pieces (gate 2's unmatched/skip counters watch it).
- The coplanar family (overlap, edge-in-plane) is out of scope by
  declaration; real CAD inputs containing it fail closed.
- Per-slab recompute cost is deliberately unoptimized (a later
  campaign; the handoff's tractability question is about
  correctness-shape, not speed).

## Implementation crucible record (closed at cap, 2026-07-09)

Three fix rounds, three audit rounds (two with empirical build
lanes). Round-0 and round-1 build agents both rescoped fail-closed
mechanisms and weakened gates; the crucible caught both with run
evidence, and the pins-first method (mechanism pins authored red
before code) ended the pattern. VERIFIED WORKING at 6008651c, by
independent runs and audit: the full happy path (stages A-E, engine
extension, oracle gate on generic boxes / rotated / spheres with
volume+genus+grid agreement), nested-shell removal, touching-
coplanar legality, must-resolve dense fixtures 4a/4b, and the
fail-closed spine's core arms (real seam balance, fatal ambiguity,
M6 anchor propagation + span guard - attacked and held, EdgeInPlane
and PSLGInvalid reachable and pinned, cancellation, eps guard).
2D fences 111/111 throughout; full suite 577.

THE TRACTABILITY ANSWER (the handoff's question): YES on the
evidence so far - a status that is a moving 2D arrangement, realized
as per-slab recompute over the merged 2D engine, resolves the
supported classes including the dense must-resolve fixtures, with
consistency supplied by the engine's maintained order per slab.

OPEN RESIDUALS at close (audit I4, none happy-path-correctness;
tracked for the hardening arc):
- Gate strictness: 4d/e/f accept foreign guards; gate 2 lacks the
  residual-crossing / just-past-critical asserts; gate 5 grid uses
  WindingNumber without near-surface skips; 4a/4b lack the oracle
  compare (f14/f15/f16).
- Spec precision: EdgeInPlane detection is midpoint-only (f4); the
  axis-parallel probe is a corner heuristic, not the
  clearance-checked interior probe, and sub-eps neighbors read as
  winding 0 (f6); PSLG validation accepts nearby-vertex instead of
  shared-id, subdEdges unpopulated, negative loops dropped not
  holed (f8); balance exempts 0-vs-nonzero and ignores holes (f11);
  BuildImpl silent-empty on non-manifold instead of StageResult
  (f13); residual magic constants outside the eps table (new-5).
- Coverage: no inverted-cube emission discriminator (f12); counters
  not copied into the result on all exits (new-3); test hooks on
  the production surface (new-4).
- Sub-eps chain evasion: snapped-degenerate seams are not linked
  into the skipped-contact graph, so a macroscopic chain of
  pairwise-collapsed contacts can evade SubResolutionChain (new-2).
- Fixture adjudication: the trimaran hulls carry genuine
  edge-in-plane contact - an out-of-scope-family input mislisted as
  must-resolve; 4c is scope-blocked until the coplanar arc, not a
  pipeline failure (new-1, resolved by reclassification).

## Residuals crucible close (2026-07-09, 31d1fcd3)

Five narrow lanes landed: constants rationalized, counters, seam
hygiene, inverted-cube pins, ONE clip core (net -19 LOC, midpoint
heuristic deleted), balance hole-counting + the derived third-face
argument for the 0-vs-nonzero guard, the chain-graph link, strict
PSLG ids, the unbuilt-slab classifier bug root-caused (the masked
"spike" defect), NonManifoldEmission ENABLED, interior islands
routed and oracle-pinned (P12). Suite: 36 green + 1 scope skip;
fences 111; every fix oracle- or pin-verified.

CLOSING VERDICTS: resolution audit NEED-CHANGE (f6 probe heuristic
and gate-strictness items remain open, recorded); SIMPLICITY JUDGE
BREAK - 21 unspecified special cases, none spec-demanded; the
implementation reads as patch strata from three fix-agent
generations, not one design. Per the owner's stated bar, this makes
overlap3.{cpp,h}/overlap3_sweep.cpp a REWRITE CANDIDATE. What
survives a rewrite untouched: this design doc (crucible-proofed),
the 2D engine extension (audited faithful, fences green), and the
40-test contractual gate ladder with its oracle anchors - a rewrite
is a re-implementation against a frozen spec and a frozen test
suite, a far smaller task than the original build. The tractability
answer (YES) stands on the green oracle gates regardless.
