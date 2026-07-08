# 3D overlap removal by sweep plane (prototype design)

STATUS: DESIGN DRAFT under adversarial review (crucible). The handoff
(.claude/plans/3d-sweep-plane-prototype-handoff.md) and its
reconciliation sibling are the framing; this doc is the concrete
design. Simple-but-robust is the bar; performance is a later
campaign. The prototype exists to answer ONE question: is a sweep
status that is itself a moving 2D arrangement tractable?

## Contract

Input: a topology-manifold triangle mesh whose triangle interiors may
cross (a self-overlapping solid; e.g. Compose of overlapping parts,
or a Boolean output with residual self-intersection). Output: a
manifold mesh that is the boundary of the input's winding > 0 region,
eps-valid in Smith's sense (7.7: topological validity, geometric
positions within the error budget), or a reported failure - never a
silently wrong mesh. Oracle: Boolean3 (`Manifold::Boolean`) computes
the same region for two-operand cases.

## The 2D lesson, carried up (settled - do not re-litigate)

- Local, order-free decisions do not compose on dense
  near-concurrences. Every classification in this design derives from
  a maintained sweep order, never from an independent per-item probe
  (no per-face ray casts; the June RSI pipeline's non-manifold
  fallbacks are the 3D empirical record of the same failure).
- Robustness lever is tolerance + consistency, not precision: shared
  `Interpolate`/`Intersect` kernels, `EpsilonFromScale` budgets, no
  exact predicates, no double-double.
- Coincident-event handling is the block rule (Smith 7.6.2), not
  symbolic perturbation.
- The 2D engine (`src/boolean2_sweep.cpp`, `SweepWinding`) is the
  validated consistency solver for a 2D arrangement: input features
  are eps-quantized once (vertex merge + incidence pre-split), the
  sweep itself is tolerance-free with exact-coordinate identity, the
  block rule collapses dense clusters through event vertices, and the
  winding pass decides fill by the maintained order.

## Design overview

Five stages. The kinetic insight that shapes them: in 2D, two edges
cross at a POINT (one sweep event); in 3D, two faces cross along a
SEGMENT that persists across an x-interval of the sweep. So the 3D
status is a 2D arrangement with PERSISTENT crossings that deform
continuously, and its TOPOLOGY changes only at discrete x-criticals.
Between consecutive criticals (a "slab"), the cross-section's
combinatorial structure is constant. That makes per-slab recompute
sound: one representative section per slab captures the whole slab.

Stage A - merge and canonicalize (the #289 steps 1-2 that both
approaches share). Merge input verts within eps (deterministic
union-find, centroid-nearest representative - the 2D engine's own
convention). Drop collapsed tris; merge exactly-identical tris
(same canonical vert triple) into one record with a signed
multiplicity (orientation-relative); zero-mult records drop.

Stage B - arrangement geometry (seams and their verts), computed
ONCE, globally, as shared objects:
- Edge-face events: for each (edge, tri) candidate from a BVH broad
  phase, the eps-aware piercing test (shared kernels); each event is
  one new vert, snapped to an existing vert within eps of it (the
  incidence quantization, applied once here - the 2D pre-split's
  analog).
- Seams: for each intersecting tri pair, the clipped intersection
  segment; its endpoints are exactly two events/verts (per #289 step
  7's invariant, which the merge is believed to guarantee - asserted,
  not assumed).
- Triple points: for each face, pairwise crossings of its incident
  seam segments IN THE FACE PLANE (2D crossing test, shared kernel);
  a crossing is one new vert shared by the three faces' seam
  structure; crossings within eps of each other or of a seam vert
  merge into it (1D sort along each seam, neighbor-merge - #289 step
  4's rule, which IS sound in 1D). Every seam becomes a polyline
  split at its triple points, and both incident faces reference the
  SAME polyline object with the same vert ids. This kills the join
  problem by construction: there is nothing to re-join later.
- On-edge incidences: events and seam verts within eps of an input
  edge subdivide that edge (ordered by t - sound in 1D), so adjacent
  faces conform across shared input edges by shared vert ids.

Stage C - the sweep: classification by slabs. x-criticals = the
sorted union of merged-vert x's, event-vert x's, and triple-point
x's. For each slab (an interval between consecutive criticals) that
is wider than eps, take the mid-x and build the SECTION: every face
straddling mid-x contributes one segment (its two edge crossings at
mid-x via `Interpolate`), tagged with (face id, multiplicity, side
orientation). Feed the section's segments to the 2D engine's
arrangement + winding passes in (y, z) coordinates. The engine
returns, for every arrangement piece, the winding below and above it
in section space (a small engine extension - see "engine extension"
below). Because the input surface is closed and even-manifold after
stage A, section winding at a point equals the 3D solid winding
along that plane, so this per-slab 2D winding IS the 3D
classification, decided by one maintained order per slab.

Stage D - classification transfer. Each face's sub-face (stage-E
region) intersects some set of slabs; it is classified from ONE slab
- the widest slab its x-extent covers - by locating its section
segment piece and reading (below, above). A sub-face whose every
covering slab is thinner than eps is a degenerate piece: it inherits
the classification of a neighboring piece across its widest seam
(recorded, counted, and reported - this is where a genuine 3D block
rule would act; the prototype measures the class instead of solving
it).

AXIS-PARALLEL FACES (found during drafting; in-scope, explicit): a
face lying in a section plane x = c never straddles any slab and has
no section segment anywhere. Its two sides face along the sweep
axis, so its classification is w(just below c) vs w(just above c) at
an interior point: probe the two ADJACENT slabs' finished 2D
arrangements with a point-winding query at a clearance-checked
interior point of the sub-face (largest-triangle centroid of the
region). Querying a SETTLED, valid arrangement is evaluation, not
the refuted pattern (which was deciding arrangement/classification
by independent probes); the query point is chosen with clearance
from all section edges, and a clearance failure falls into the
degenerate-piece rule. Nearly-axis-parallel faces section normally
but may have sub-eps x-extent; they fall under the thin-slab /
degenerate-piece rules and are counted in the same report.

Stage E - per-face partition and emission. Each face's PSLG = its
(subdivided) boundary + its seam polylines; regions are found by the
face-local planar walk and triangulated (`Triangulate`, projected to
the face plane). A region is EMITTED iff IsInside(below) !=
IsInside(above) (positive fill: w > 0 on exactly one side), oriented
with the kept material on the winding>0 side; region multiplicity
from stage A rides along (opposite coincident faces cancel). The
PSLG must already be a valid arrangement - seams crossing anywhere
but shared triple-point ids means stage B missed an event, which is
a HARD FAILURE reported by an assert-and-fallback, never a local
repair (repairs are the refuted pattern).

Manifoldness is then a THEOREM to check, not a mechanism to add:
seam ids are shared by construction (B), keep decisions are
winding-consistent by the maintained order (C), so every output edge
must pair. The output gate verifies it (every edge shared by exactly
two output triangles with opposite orientation, start/end balance at
triplet verts - the 2023 constraint as a check).

## The engine extension (minimal, the only 2D change)

`SweepWinding` today emits only the retained boundary. The winding
pass already computes (below, above) for every piece inside its
block loop (`EmitBoundary(from, to, m, below, above)`). Extension:
an optional out-channel that records (from, to, below, above) for
EVERY piece instead of only retained ones - a flag and one push_back,
no behavioral change to existing callers. Piece->face matching is
geometric: a piece lies on its source segment up to the engine's
constructed-point error; match by point-to-segment distance at eps
against the slab's tagged input segments, nearest wins, unmatched
pieces reported. (Threading ids through PolySet2 would be exact but
touches the engine's core key type; rejected for the prototype -
the matching is measured by gate 2 instead.)

## Epsilon posture

One absolute eps = `EpsilonFromScale(bbox scale, 1000)` unless the
caller passes one (identical to 2D). It quantizes input features
(stage A merge, stage B snapping/merging, on-edge subdivision) and
bounds the stage-D matching; the 2D engine inside each slab runs its
own standard posture on section coordinates (same scale). Slabs
thinner than eps are never classified from. The >= 2^40 coordinate
degeneracy is accepted manifold-wide, per the handoff. Smith's alpha
applies to every constructed point (edge-face verts, triple points,
section endpoints - all via the shared kernels); the budget-1000 eps
covers them exactly as in 2D.

## Deliberately out of scope (prototype)

- Coplanar face overlap beyond exactly-identical triangles (the
  hardest 2D-in-3D class; no seams exist between coplanar faces, so
  their mutual imprint is a separate mechanism - the June 6.5 trace
  machinery is the reference if it is ever imported). Detected
  (eps-coplanar overlapping pair with non-identical verts) and
  reported as unsupported, oracle fixtures avoid it.
- Cancellation/ExecutionContext, progress, parallelism, performance
  (per-slab recompute is deliberately quadratic-ish).
- Public API: the seam is internal (`RemoveOverlaps3D(const
  Manifold::Impl&, double eps) -> std::optional<Impl>`); tests reach
  it directly. Public wiring is post-prototype.
- Simplification/decimation of the output (Simplify's job, as in 2D).

## Validation gates -> tests (test/overlap3_test.cpp)

1. EVENT PARITY: brute-force O(n^2) edge-face and seam-pair
   enumeration equals stage B's output on small fixtures (two tets,
   two boxes, box+rotated box).
2. SECTION VALIDITY: per slab, the section is closed (every 2D vert
   has even degree; signed multiplicity sums zero around it) and the
   engine ingests it without assert; piece->face matching resolves
   every piece (unmatched count == 0 on the fixtures).
3. TRIPLET PAIRING / MANIFOLDNESS: on fixtures with genuine
   triple-face verts (three boxes pairwise overlapping around a
   common region; three plates through one line region), the output
   passes the manifold gate (paired edges, balanced verts) and
   Manifold(Impl) construction.
4. DENSE NEAR-CONCURRENCE (the adversarial gate): k thin wedges
   rotated about a near-common axis through one near-point (the 2D
   killer lifted); near-parallel face bundles eps apart; the June
   trimaran hulls (OBJ fixtures imported from the old branch) as
   subtraction leftovers. PASS = manifold output + winding oracle
   agreement; measured, not assumed.
5. ORACLE: for two-operand fixtures, RemoveOverlaps3D(Compose(A, B))
   vs Boolean3 A+B: volume within eps-derived bound, genus equal,
   Contains() agreement on a deterministic point grid (the
   Monte-Carlo winding oracle's 3D analog, deterministic seed).

Plus: the existing manifold suite stays green (the branch adds files
and one flagged engine out-channel; nothing else changes).

## File plan

- `src/overlap3.h` - internal seam + stage structs (~100 lines).
- `src/overlap3.cpp` - stages A, B, D, E (~600-900 lines).
- `src/overlap3_sweep.cpp` - stage C: slabs, section build, engine
  calls (~250 lines).
- `src/boolean2_sweep.cpp` - the piece-capture out-channel (~15
  lines).
- `test/overlap3_test.cpp` - the gate ladder (~600 lines).
- CMake: add the sources to the manifold target, the test to
  manifold_test.

## RISKS (the review lanes)

R1 CROSS-SLAB CONSISTENCY (the core bet). Sub-faces are classified
from single slabs; manifoldness needs keep-decisions consistent
across every shared edge whose incident pieces may classify from
DIFFERENT slabs. The design's claim: both slabs' sections are built
from the same stage-B objects evaluated by the same kernels, and a
correct 2D winding of a valid section is unique - so two slabs can
only disagree if a section is built WRONG (a missed critical, a
face-straddling misjudgment at a slab boundary). Attack: construct a
disagreement; find a topology change inside a "slab" that stage C's
critical set misses (tangencies? a seam whose x-extent endpoints are
not events? vertical faces?).

R2 DENSE 3D CONCURRENCE. There is no true 3D block rule here: stage
B merges triple points at eps (1D along seams) and the 2D block rule
acts only INSIDE each slab's section. A dense 3D tangle (many faces
near-concurrent through a point) produces: many near-coincident
triple points (merged at eps - order-dependence of that merge?),
degenerate slabs (skipped), and degenerate sub-faces (classification
inherited). Attack: break the inheritance rule; construct a tangle
where eps-merged triple points produce an invalid face PSLG or an
unpaired output edge.

R3 PER-FACE PSLG VALIDITY. Stage E asserts seams cross only at
shared ids; FP residual crossings between near-parallel seams within
a face are exactly the class the 2D engine solves with forced-through
- and the design REFUSES to run an independent per-face arrangement
(it would invent per-face topology and break cross-face sharing).
Attack: is assert-and-fail acceptable (how often does it fire on
gate-4 inputs), or does the design need a principled in-face
resolution that preserves sharing?

R4 SECTION ORIENTATION AND WINDING SIGNS. Section winding == solid
winding needs the right sign conventions (face normal -> section
segment orientation in (y,z); Add rule w>0; emission orientation).
Attack: work a box and an inverted box through the full sign chain
on paper; find a flipped case (vertical faces, faces parallel to the
sweep axis are IN the section plane - how are they sectioned?).

R5 PIECE->FACE MATCHING. Geometric matching at eps after the engine
splits pieces at constructed points. Attack: a slab where two faces'
section segments are within eps (near-parallel bundle) - nearest-
segment matching misattributes; quantify and decide whether the
id-threading alternative is required after all.

R6 STAGE-B INVARIANTS. #289 step 7's "exactly two shared verts per
seam" is believed-with-proof-wanted; tangent and edge-through-edge
cases produce 1 or 0. Attack: enumerate the degenerate seam
configurations and check each resolves to a stated rule (merge, skip,
or event) rather than an unstated assumption.

R7 SIMPLICITY AUDIT. Anything here not needed to answer the
tractability question? Anything manifold already provides being
reinvented? Is the engine extension truly minimal?
