# Plan: resolve interior-island stamps (faces with holes)

STATUS: IMPLEMENTED (see `src/overlap_removal.cpp`, `src/overlap_removal_internal.h`, `test/manifold_test.cpp`)

Working plan for the highest-ranked CAD limitation: a shell stamping
through the INTERIOR of a single face (the boss-through-plate class,
Known limitations item 8). Today the per-face partition cannot express
a face with a hole, so the island gate fails closed and the run
returns its input. This plan makes the partition hole-aware; the gate
remains as the fail-closed second line.

## Current behavior (ground truth)

- `PartitionFace` builds the face's sub-edge graph (boundary sub-edges
  + chord sub-edges), then BEFORE walking, detects island components:
  union-find over chord sub-edges (`seenSub`); a component with cycle
  rank (`compEdges >= compVerts`) and fewer than two distinct
  boundary-vert attachments sets `interiorIslandVerts` and returns;
  the driver gates (`return nullopt` -> caller returns its input).
- Pre-gate history (the TP5-era finding that motivated the gate): when
  islands were NOT gated, the angular walk traversed the island loop
  in BOTH orientations and emitted the outer boundary plus the loop
  twice; step 12's signed-multiplicity merge cancelled the
  opposite-orientation pair, silently erasing the cut. So the walk
  DOES reach island components (cycle extraction seeds from every
  directed sub-edge, not only boundary halfedges) - CONFIRMED by the
  round-1 plan review against the current walk (chords insert both
  directions; every directed sub-edge seeds; visited is
  per-direction).

## Design

### 1. Component reclassification (in PartitionFace)

Cycle-bearing chord components split by attachment count:
- `compAttach == 0`: hole CANDIDATE, but only after proving it is a
  CLEAN SIMPLE LOOP - cycle rank alone (the current detector's test)
  does not establish that. Ungate the component only when ALL hold:
  every component vert has degree exactly 2 within the component (a
  connected degree-2 component is exactly one closed loop); the loop
  is SIMPLE in the face frame - the pairwise non-adjacent segment
  check rejects every non-simple contact class (strict crossings,
  vertex-on-nonincident-edge, endpoint touches beyond shared loop
  adjacency, collinear overlaps), the same PSLG rules the
  triangulation validation applies, because degree-2 proves closure,
  not embedding; the loop
  has nonzero projected area in the face frame; no component vert is
  shared with any other island component; and no component vert is
  incident to a boundary-rider sub-edge that the walk SKIPPED before
  detection (a skipped rider means the component is attached after
  all - track the skipped set and check). Anything else (branchy
  detached graphs, zero-area loops, vert-sharing island clusters):
  GATE as today.
- `compAttach == 1`: PINCHED - stays GATED exactly as today (the walk's
  behavior at the pinch vert was the TP6-era false-negative class;
  resolving it is out of scope).
- `compAttach >= 2` or tree components (no cycle rank): partitionable
  as today.

### 2. Walk unchanged, then post-walk decomposition

SIGNATURE CHANGE: `PartitionFace` gains an explicit `double eps`
parameter (no default - the severance-proofing convention); the driver
passes the pipeline eps it already holds, and the seam tests pass it
explicitly. The decomposition's triEps derives from it as
`max(impl.tolerance_, impl.epsilon_, eps)`, the emit's existing seed.
`PartitionFace` also gains the driver-computed COPLANAR-PARTNER flag
described under failure arms.

Run the angular walk as today (no gate for 0-attachment islands). The
raw output for an island face contains, per island loop: the loop in
both orientations (signed areas +A and -A in the face frame), plus the
outer region cycle(s).

Decomposition pass (only when islands were detected):
- Compute each walked cycle's signed area in the face's existing 2D
  frame (the projected basis the walk already uses).
- POSITIVE-area cycles are region boundaries. They include both the
  outer region(s) and, for each island, the CCW orientation of the
  loop - which legitimately bounds the disk INSIDE the stamp footprint
  (that piece of the face exists and is later classified/dropped by
  winding like any other polygon). These pass through as ordinary
  simple polygons.
- Sub-resolution holes: the triangulator's own hole classification
  (`FindStart`) involves more than area - degenerate clipping and
  reflex-visibility rules - so the plan does NOT claim a verbatim
  precheck mirror. The DOUBLED-area threshold (origin-relative
  determinant sum without the 1/2 factor, compared against
  epsilon * max bbox side at the exact triEps) is a CONSERVATIVE
  fast-fail only: clearly-sub-resolution islands gate early. The
  AUTHORITY on whether the triangulator actually honored the hole is
  the post-triangulation validation below - a hole misclassified as a
  non-hole cannot satisfy boundary coverage (its sub-edges would not
  appear correctly oriented exactly once), so it gates there.
- Each island's NEGATIVE-area cycle is a HOLE. Candidate regions for
  it are positive-area cycles EXCLUDING its own orientation twin (the
  positive cycle over the same undirected canonical vert sequence -
  the twin shares every vert with the hole and would otherwise win a
  naive smallest-containing test) and excluding any cycle that shares
  ANY vert with the hole (a clean detached island per section 1
  shares verts with nothing but its twin). Among the remaining
  candidates, require STRICT point-in-polygon containment of a
  representative hole vert - a representative landing exactly ON a
  candidate's edge is treated as NOT contained (borderline incidence
  gates rather than guesses; the on-edge machinery exists upstream
  precisely because these cases are real). Pick the smallest area;
  break exact area ties by the candidates' deterministic walk order.
  No candidate -> gate. Nesting (stamp within a stamp) falls out: the
  inner island assigns to the inner disk region, which does not share
  its verts.
- For each region with holes: build a `PolygonsIdx` {outer
  positive-area, holes negative-area - the triangulator's confirmed
  convention: `FindStart` classifies `area < -minArea` contours as
  holes} carrying vert ids, call `TriangulateIdx` at the same epsilon
  ladder seed the emit uses (`max(impl.tolerance_, impl.epsilon_,
  eps)`), and append the resulting triangles to the face's polygon
  list IN FACE-WINDING ORIENTATION as ordinary walked cycles -
  `MergePolygons` computes canonical signed multiplicity itself; the
  decomposition synthesizes no `MergedPolygon` and no multiplicity.
  The returned indices are our vert ids (confirmed: `PolyVert.idx`
  rides through to the output triangles); NO new verts are
  introduced.
- Regions without holes are untouched.

### 3. Failure arms (the gate's new role)

Release `TriangulateIdx` does not throw and may return a non-empty
but invalid triangulation for bad input (its documented posture), so
the decomposition VALIDATES its output in every build config:
- per-triangle: three DISTINCT vert ids, STRICTLY POSITIVE signed
  area in the face frame - a NEW decomposition invariant (stricter
  than the partition's own exact-zero drop, which deliberately admits
  negative island walks and tiny slivers; triangles, unlike walked
  cycles, must all be face-winding regions) - and not canonical-equal
  to any hole contour (the 3-vert hole case);
- signed-area preservation: sum(triangle areas) == outerArea +
  sum(signed hole areas) - the holes are NEGATIVE, so this is
  outerArea minus their magnitudes; spelled signed to prevent the
  double-negation implementation bug. Tolerance (area-dimensioned, the
  only implementable form): |difference| <= triEps * regionBboxScale +
  kAreaRelTol * |expected|, with kAreaRelTol a small relative constant
  (1e-9 order) - NOT the pierce predicates' kPipelineRelTol, which is
  scoped to dimensionless ratios;
- boundary coverage: every outer and hole sub-edge appears exactly
  once, correctly oriented, among the triangle edges; every triangle
  edge is either a boundary/hole sub-edge or shared by exactly two
  triangles of this region;
- planar embedding, full PSLG rules (strict crossings alone are not
  enough: a T-junction, collinear overlap, duplicate edge, or
  vertex-on-edge passes a strict-crossing test yet breaks
  BuildCellComplex, which keys fans by exact edge endpoints): no two
  non-adjacent triangle edges strictly cross; no region vert lies on
  a non-incident triangle edge; no two distinct undirected edges
  overlap collinearly; no undirected edge appears beyond the
  boundary-once / interior-twice coverage rule. All tests in the face
  frame. Holed regions are small in practice (the boss class: a
  triangle outer and a handful of hole verts), so the quadratic
  checks are negligible; their cost is bounded by gating regions
  above a generous edge-count cap rather than skipping them.
COPLANAR-COINCIDENCE SCOPE CUT: step 12's cancellation matches
canonical cycles exactly, and two eps-coincident OPPOSITE holed faces
triangulated independently (opposite winding, mirrored frames) are not
guaranteed to choose matching diagonals - unmatched triangles would
defeat the cancellation that step 6.5 exists to feed, resurfacing
exact angular ties in BuildCellComplex. Canonical
winding-and-mirror-independent triangulation is research territory:
OUT OF SCOPE. Instead `CoplanarTraceChords` EXPOSES a per-face
CANCELLATION-HAZARD flag: among its plane-gated pairs (taken BEFORE
interval filtering, so equal-boundary coincidences whose
boundary-riding intervals are rejected still count), flag ONLY pairs
whose face WINDING normals - computed locally as normalized
cross(v1 - v0, v2 - v0), the 6.5 pass's existing convention, never the
stored faceNormal_ (which can oppose winding) - are ANTI-ALIGNED and
whose in-plane bounding boxes overlap. Same-oriented coplanar neighbors - every triangulated flat
face's own tris, which the plane gate also pairs - can never form a
step-12 cancellation pair (cancellation requires opposite
orientations) and MUST NOT be flagged: flagging them would re-gate
the boss-through-plate fixture itself, whose stamped face has a
same-oriented coplanar neighbor by construction. The driver passes
the flag set into `PartitionFace`; a flagged face with a detected
island gates exactly as today. The headline boss-through-plate class has no
coplanar partner on the stamped face and is unaffected; the cut is
recorded in Known limitations.

The `TriangulateIdx` call gets its own LOCAL #ifdef MANIFOLD_DEBUG
try/catch inside PartitionFace (the existing broad catch sits at the
RemoveOverlaps entry and can neither set `interiorIslandVerts` nor
serve the direct seam tests; the emit's retry catch is local to the
final output triangulation). On a caught debug throw, a validation
failure, a sub-minArea hole,
or a hole with no candidate region: set `interiorIslandVerts` and
return - exactly the current gate; the driver falls back
bit-identically. ONE rule for emptiness: ANY empty triangulation for
a region-with-holes gates (it fails the validation trivially; no
separate vert-count arm exists). Fail-closed posture unchanged; the gate covers the
pinched class, unclean detached components, and decomposition
failure.

### 4. Why downstream stages need no changes (to be re-verified by review)

- Triangles are simple cycles like every other partition polygon,
  appended in face-winding orientation; `MergePolygons` canonicalizes
  and signs them itself, exactly as it does the cycles they replace.
- The TP5-era cancellation hazard dissolves structurally: the CW
  island cycle is consumed INTO the triangulation (its sub-edges
  become triangle edges); it never reaches step 12 as a standalone
  polygon, so there is no opposite-orientation pair to cancel.
- Face-local diagonals introduced by triangulation are sub-edges
  shared by exactly 2 same-plane polygons on OPPOSITE sides with
  antiparallel in-face directions - the precise invariant
  `BuildCellComplex` digests on every flat quad's interior mesh edge
  today (exact angular ties arise only for SAME-direction pairs,
  which a valid triangulation cannot produce across a shared
  diagonal).
- The island boundary sub-edges themselves are chord sub-edges shared
  with the stamping shell's wall faces - the same 4-side radial fans
  any chord produces.
- No eps semantics change: island verts are ordinary arrangement verts
  (their conditioned radii were already handled at 9.5); the
  triangulation epsilon is the emit's existing ladder seed.
- Determinism: `TriangulateIdx` is deterministic for identical input;
  cycle and hole ordering derive from the walk's deterministic order.

## Tests

Unit (PartitionFace seams):
- `Step10InteriorIslandDetectedAndFailsClosed` SPLITS: the free-island
  half flips to expect hole-aware decomposition (disk polygon +
  annulus triangles, `interiorIslandVerts == 0`, signed area
  preserved); the pinched half becomes its own fail-closed pin; the
  proper-double-crossing pass case stays as is.
- Unclean detached components (branchy loop-plus-spur; two loops
  sharing a vert; zero-area loop): all gate.
- Two disjoint islands in one face; island-within-island (nested) if
  the fixture is cheap.
- Decomposition-failure arm: forced validation failure falls back to
  the gate (witnessed).

Feature:
- `RemoveSelfIntersectionsInteriorIslandFallsBack` FLIPS to a resolve
  pin: pierces == 0, one component, volume == the union (cube 1.0 +
  stamp 0.2 x 0.2 x 0.6 minus the embedded 0.2 x 0.2 x 0.3 overlap =
  1.012), idempotent on re-run. Witnessed red against the
  decomposition being disabled (gate restored).
- The fallback CLASS keeps coverage through the pinched unit pin (a
  feature-level pinched fixture is optional, not required).

Docs: Known limitations item 8 rewrites to "resolved for detached
islands; pinched (single-attachment) loops and
coplanar-cancellation-hazard faces remain gated"; pipeline overview
row for steps 10-11 mentions the hole decomposition; and the PUBLIC
doc comment in include/manifold/manifold.h, which today names "an
interior-island chord loop" among the fail-closed gates, updates to
the narrowed residual (pinched / hazard-flagged / validation-failure
arms).

## Stop-rule tripwires (escalate to the user, do not improvise)

- Hole-to-region assignment turning out to need exact predicates at
  conditioned geometry (point-in-polygon going borderline in real
  fixtures).
- `BuildCellComplex` exact angular ties on the k=2 coplanar diagonals
  that the flat-quad precedent does NOT already cover.
- Step-12 / 6.5 interaction when stamp walls are eps-coplanar with the
  stamped face (the pancake class meeting the island class).
- Any need to modify polygon.cpp / the triangulator itself.

## Out of scope

Pinched (1-attachment) loops; properties; the fold-gate class; any
global arrangement change. No new verts, no new eps, no new gates.

## Simplification addendum

After landing, a proportionality pass over the validation weight was
approved (docs/SimplificationPlan.md). The PSLG output-validation suite
(strict crossings, vertex-on-nonincident-edge, collinear overlap,
duplicate-edge analysis) and the edge-count cap were removed. The
surviving triad (degenerate triangle, signed-area preservation, boundary
coverage) remains. The rationale: manifold extends `TriangulateIdx`
unchecked trust at every other call site on weaker input guarantees;
the island call site's input is fully proven before the call, making
the PSLG suite disproportionate. The edge-count cap was deleted with the
suite it bounded; large valid holed regions now decompose instead of
gating. A sub-resolution island pin (Step10SubResolutionIslandGates)
was added to cover the one deterministically forceable fast-fail path
previously unpinned.
