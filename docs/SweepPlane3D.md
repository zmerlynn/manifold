# 3D overlap removal by sweep plane (prototype design)

STATUS: DESIGN, crucible rounds 1-2 folded (round 1: 5 lanes, 4
convergent BREAKs resolved, marked [R1-fold]; round 2: re-attack +
fresh-eyes + empirical id-plumbing build, all round-1 constructions
HELD-resolved, narrower findings folded, marked [R2-fold]). The handoff
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
    into connected components first; a component of total 3D
    diameter <= eps drops and is counted (a genuinely point-like
    contact); a larger component is the SUB-RESOLUTION SEAM CHAIN
    class and fails closed, named;
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

Stage D - partition and classification. The face-local planar walk
(purely topological) runs FIRST: each face's PSLG = its subdivided
boundary + its seam polylines; regions are the walk's output.
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

[R2-fold, adjudicated] A dropped degenerate component whose interior
carried a fill transition is NOT the silent-wrong class: such a
component is by construction sub-eps-thick along the sweep axis, and
dropping sub-eps features is the documented eps-scale decision (the
2D engine's own input quantization; #289 step 1's "smallest feature
size that is desired to retain"), applied CONSISTENTLY (both sides
of every seam, enforced by the balance check below) and COUNTED in
the report. Eps-validity admits it; the report makes it visible.

[R1-fold] SEAM BALANCE is enforced BEFORE emission, not assumed
after: for every seam polyline edge, the kept incident regions must
pair (equal start/end counts - the 2023 triplet constraint as a
pre-emission check). A violation names the seam and fails closed.
Manifoldness of the output is then the gate-3 CHECK; the design no
longer claims it as a theorem, it engineers toward it and verifies.

Stage E - triangulation and emission. Each kept region (fill
transition per the orientation algebra above) triangulates via
`Triangulate` projected to the face plane; exactly ONE sheet per
transition regardless of multiplicity. PSLG VALIDITY IS CHECKED
BEFORE THE STAGE-D WALK RUNS [R2-fold ordering]: seams crossing
anywhere but shared ids is a stage-B miss and a hard failure, never
a local repair, and no later stage consumes invalidated data. (With
global unification in stage B, the round-1 false-assert class -
same point, two ids - is resolved at its root.)

FAILURE CONTRACT, one shape [R2-fold]: every stage returns a
StageResult carrying either its product or a fatal reason code
(diameter guard, sub-resolution chain, classification ambiguity,
anchor conflict, starvation, balance violation, PSLG invalidity);
fatal stops the pipeline at that stage. Non-fatal COUNTERS
(sub-eps contacts dropped, degenerate components classified by
propagation, eps-feature drops) accumulate in a report struct
returned alongside success - visible, never fatal.

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
  known-limitation-10 class). Fixtures avoid it; occurrences surface
  as seam-balance failures.
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
4. DENSE NEAR-CONCURRENCE (adversarial) [R2-fold: non-vacuous]: the
   fixture list is SPLIT. MUST-RESOLVE stress fixtures (k thin
   wedges at feature scale through a common region; moderately
   near-parallel bundles well above eps; the June trimaran hulls as
   subtraction leftovers): PASS requires manifold output + oracle
   agreement - a fail-closed here FAILS the gate (no vacuous
   safety). MUST-FAIL-CLOSED degeneracy fixtures (bundles inside
   eps; sub-resolution chains; authored tolerance-scale features):
   PASS requires the NAMED guard firing - a resolved-but-wrong mesh
   or an unnamed crash fails. The anchor-propagation and
   diameter-guard rates on both lists are the tractability data the
   prototype exists to produce.
5. ORACLE: for two-operand fixtures, RemoveOverlaps3D(Compose(A, B))
   vs Boolean3 A+B: |volume difference| <= eps *
   max(surfaceArea(ours), surfaceArea(oracle)); genus equal;
   Contains() agreement on a deterministic point grid.

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

## RISKS (round-3 review lanes - the convergence check)

R1'' Re-attack the round-2 constructions verbatim against the folds
(the sub-eps chain cluster rule; the oriented-agreement rule; the
correspondence spec; the orientation algebra on a +/-y face) and
audit the folded design for internal contradictions introduced by
two rounds of edits.
R2'' Fresh-eyes full read: is the design now implementable as
written by an engineer who has seen none of the review history -
every stage contract stated, every constant named, every failure
path reachable and typed?
