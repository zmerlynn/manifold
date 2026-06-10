# Overlap Removal (RemoveSelfIntersections)

`Manifold::RemoveSelfIntersections()` returns a new manifold whose surface is
the boundary of the input's winding > 0 region, with no more geometric
self-intersections than the input. Boolean operations produce
topology-manifold output (every edge shared by exactly two triangles), but the
result may still contain triangle pairs whose interiors cross. This pass
resolves the full surface arrangement and reclassifies it by winding, or
returns the input unchanged when it cannot.

It is a faithful implementation of Emmett Lalish's 13-step sketch
(#289, the 2024-05-14 comment), with no recovery scaffold: no cap walker, no
pierce reducers, no post-hoc repair. Every kept output triangle exists because
the global winding classification kept it. Implemented in
`src/overlap_removal.{h,cpp}` + `src/overlap_removal_internal.h`; the only
production caller is `Manifold::RemoveSelfIntersections()`.

Two passes (6.5 and 9.5 below) are not in the sketch. Both exist because the
sketch is written as if arithmetic were exact; they make its implicit
assumptions hold in doubles, and add no new algorithmic ideas.

## Pipeline overview

| step | function | role |
|---|---|---|
| 1 | `MergeVertsEps` | eps-merge near-coincident verts (cluster centroids, MeshGL64 merge hints) |
| 2 | `EnumerateEdges` | canonical undirected edges with halfedge pairs |
| 4 | `BuildOnEdgeVertLists` | verts within eps of an edge's interior, sorted by t |
| 5 | `BuildOnTriVertLists` | verts within eps of a tri's interior (strict barycentric) |
| 6 | `FindEdgeTriIntersections` | transversal edge-pierces-tri events (BVH + Moller-Trumbore), snapped to the nearest existing vert within eps |
| 7 | `GenerateChordEdges` | one chord per tri-tri pair with exactly 2 pierce endpoints; event-vert resolution |
| 6.5 | `CoplanarTraceChords` | in-plane conformance cuts between coplanar overlapping faces (below) |
| - | `AddVertsToOnEdgeLists`, `PropagateNewVertsToOnEdgeLists` | pierce + trace verts subdivide the edges they lie on |
| 8 | `AddInteriorVertsToNewEdges` | on-tri verts threaded onto chord interiors |
| 9 | `FindOnChordEndpointContacts`, `FindChordChordCrossings`, `MergeAndPropagateCrossings`, `ResolveAndThreadClusters` | chord-chord crossings within each face (below) |
| 9.5 | `UnifyArrangementVerts` | arrangement-wide new-vert unification (below) |
| 10-11 | `PartitionFace` | per-face simple-cycle partition by angular walk |
| 12 | `MergePolygons` | canonical-cycle merge with signed multiplicity; coincident opposite pairs cancel |
| 13 | `BuildCellComplex`, `ClassifyCells` | radial fans -> volume cells -> seed cast -> winding BFS -> keep |
| emit | `BuildEmitTopology` + driver | inside-wedge twins, vertex rings, triangulation, MeshGL64 |
| gate | driver | folded-shell volume (pre-emit) + status / volume / pierce-monotonicity, else return input |

The driver (`RunOverlapRemovalImpl`) composes these serially (per-face
parallelization is an open follow-up; nothing is parallel yet). Exceptions
exist only in MANIFOLD_DEBUG builds (optional_assert.h defines the error
types there; release manifold is exception-free and all RSI failure arms
are status-based) - under MANIFOLD_DEBUG the entry point wraps the body in
try/catch so a throwing assertion anywhere falls back to the input, the
polygon.cpp guard pattern (two inner debug-only catches have local roles:
the emit triangulation retry ladder and the seed-cast target skip). EARLY-EXIT: if the combined chord list
(transversal + trace) is empty, the input is returned bit-identical - this
covers clean inputs, the all-pairs-dropped case, and pancake-free coplanar
contact.

## House terminology and style (from the boolean2 review logs)

- halfedge structure (not DCEL); near-line sliver (not T-junction);
  nearby-crossing merge (not duplicate-crossing); on-edge / on-tri / on-chord
  vert lists for verts lying on a primitive within eps.
- Determinism is hard: any `std::` trig deciding output geometry is a bug -
  use `math::` (vendored deterministic); the current design needs no trig at
  all (half-plane bucket + cross sign ordering throughout).
- Return result structs, never out-pointers (mutating-reference parameters
  only where the sibling precedent does, e.g. on-edge list insertion).
  `DEBUG_ASSERT` not `assert` (compiled out unless MANIFOLD_ASSERT &&
  MANIFOLD_DEBUG - not a production guard). `la::cross`/`la::dot` directly.
  `static_cast`, no C-casts. Dimensionally-correct thresholds: eps is a
  LENGTH; area ~ L^2 compares against length * eps.
- TDD: red first, failing for the stated reason.

## The eps contract

One absolute pipeline epsilon `eps` (from `InferEps` =
`AlphaBudgetEpsilon(bbox scale, 1000)` unless the caller passes one), used as:

- 1x eps: kernel acceptance, on-chord propagation, trace-interval length and
  interior-margin qualification, grazing guards, new-to-new event dedup, AND
  step 6's pierce-event vert snap. The step-6 snap is EVENT IDENTITY among
  one mesh's stored coordinates - the step-1 old-old scale, absorbing only
  the pipeline's own error - not an allocation snap: events in the
  (eps, 10 eps] band are unified onto originals by step 9.5 anyway, and a
  tolerance-scale radius on a tolerance-inflated input would drag pierce
  events onto far verts and deform the arrangement.
- `tolerance + eps` (with `tolerance = max(impl.tolerance_, eps)` from the
  post-merge impl): snaps that resolve NEW verts against existing ids at
  ALLOCATION time - step 6.5's corner-snap base, step 9's crossing
  resolution and on-chord contacts - matching boolean2's `newToOldThresh`
  (prior drift plus current-op error).
- 10x eps: the nearby-crossing merge radius (new-to-new for the SAME point
  computed twice), matching boolean2's `kIntersectionMergeEpsFactor`; reused
  by step 9.5's sweep.
- eps / sin(angle), capped at `kCondSnapCapEps` (128) x eps: the CONDITIONED
  radius of a computed crossing whose defining lines are near-parallel (the
  lever arm). Used by step 6.5's corner snap and new-to-new dedup and by
  step 9's crossing records, and carried per-vert into step 9.5's
  original-vert snap. Conditioned radii apply only
  where the conditioning is computable at the source; blanket widening was
  tried and rejected (it moves real geometry and re-pierces - see
  Known limitations).
- eps / len: per-chord t-space guards and dedup backstops.
- `ClassifyCells` receives the driver's eps as a hint and runs at
  max(hint, impl.epsilon_, machine eps at the arrangement scale) - without
  the hint it would probe at `impl.epsilon_`, which can be tighter than the
  epsilon the arrangement was built with.
- OUTPUT tolerance: `max(tolerance, 10 * eps, measured step-1 merge
  displacement, measured step-9.5 remap displacement)` - the pipeline
  deliberately moves verts by up to the merge radius, an eps-pair CHAIN or a
  conditioned snap can move one further, and claiming the input tolerance
  would overstate the output's precision. Ill-conditioned shallow-incidence
  corners can carry residual error beyond this, up to the conditioned band -
  a documented limitation, not part of the tolerance claim.

## Step 6.5: coplanar trace chords

Step 1's eps-merge flattens Boolean SoS slivers into zero-volume pancakes -
coplanar overlapping faces with opposite orientations - by design: step 12's
signed-multiplicity cancellation consumes the coincidences. But step 6 finds
only TRANSVERSAL intersections; coplanar pairs produce no chords, the sheets'
partitions never conform in-plane, cancellation cannot fire, and
`BuildCellComplex` hits exact angular ties. This pass does the in-plane
cutting the sketch implicitly assumes:

1. **Detect.** BVH tri-tri broad phase; plane gate: all six verts within eps
   of the LARGER face's plane (a near-zero-area sliver's own plane is noise).
2. **Single-frame clip.** ALL geometry for a pair is computed in ONE frame
   (the gate face's `FaceBasisFromNormal` - the larger face, whose plane the
   gate already trusts; ties to the lower face id); each face's 3 edges clip
   against the other's projected triangle (convex clip, one interval per
   edge). The same geometric crossing is computed once and shares its id
   across both clip directions by construction.
3. **Interval qualification.** An interval becomes a chord iff (a) longer
   than eps AND (b) its MIDPOINT is interior to the other face by > eps.
   Midpoint, not endpoints: a full-through cut - the generic overlap case -
   has BOTH endpoints on the boundary yet an interior midpoint; an
   endpoint-margin predicate would reject exactly the cuts conformance needs.
   Distance-to-boundary is concave along a segment in a convex face, so the
   midpoint margin >= half the deepest penetration: only dips shallower than
   2 * eps are rejected, the same FP-degenerate band the grazing guard skips.
   Boundary-riding intervals (coplanar NEIGHBORS - any flat region of any
   mesh) never qualify, so clean flat meshes keep the early-exit, and
   equal-size face-glued solids pass through bit-identical.
4. **Endpoint ids.** t = 0/1 endpoints use the source edge's vert id.
   Crossing endpoints snap to the nearest of the pair's six corners within
   max(tolerance + eps, the CAPPED conditioned radius eps/sin(angle)) -
   the cap applies to the conditioned term only, so an inflated mesh
   tolerance still widens the base - nearest, ties to smallest id (the
   step-9 convention); else allocate,
   deduping new-to-new over the whole pool at the SOURCE-GATED radius
   min(this crossing's conditioned radius, the pool entry's recorded radius),
   eps-floored - an ill-conditioned crossing must not claim an unrelated
   well-conditioned vert. The allocation's conditioned radius is recorded
   per vert (`TraceChordResult::newVertSnapR`, widened when a dedup hit
   claims more) for step 9.5. A
   lifted endpoint farther than eps from either original 3D edge rejects its
   interval (the near-grazing guard - conformance we provably cannot compute
   is skipped, not corrupted).
5. **Chords + boundary conformance.** Each qualifying interval emits a
   `PiercedNewEdge{v0, v1, faceA, faceB}` - it lies on BOTH coplanar faces,
   so the existing chord model applies verbatim (dedup within the pair; the
   partition's per-face dedup absorbs cross-pair repeats). New crossing verts
   on original mesh edges return as explicit (edge, vert, t) additions - an
   X crossing gets one record per edge with one shared vert id - applied by
   `AddVertsToOnEdgeLists` (id-dedup + per-edge t re-sort).
6. **Compose.** Trace chords append to the chord list BEFORE the early-exit,
   before `GroupChordsByFace`, and before step 8 - one grouping then feeds
   step 8 extras, step 9 crossings (trace-vs-regular and trace-vs-trace),
   and the partition with no further changes.

Semantics note: a coincident interior wall separating winding 1|1
(differently-sized solids glued face to face) DROPS by the step-13 keep rule -
the output is the winding-faithful welded solid. Correct #289 behavior, not a
regression: the input's w > 0 region IS one solid. Equal-size glued faces
never reach the pipeline (no qualifying intervals; early-exit).

## Step 9: chord-chord crossings within each face

Completes the per-face arrangement: after step 9 the chord set is CONFORMING
(no two sub-edges cross except at shared vert ids). Geometry only - no
classification.

0. **On-chord endpoint pass** (`FindOnChordEndpointContacts`): a chord
   endpoint lying on another same-face chord's interior, within
   tolerance + eps and inside the t-guard
   (t in ((tol+eps)/len, 1 - (tol+eps)/len) - the endpoint-proximity zone is
   excluded in t-space), recorded into the explicit `OnChordContact`
   accumulator. Required because the crossing kernel REJECTS near-endpoint
   crossings (`AwayFromEndpoints`); without this pass those contacts are
   silently lost as near-line slivers.
1. **Group** (`GroupChordsByFace`): chord indices with triA==t or triB==t.
2. **Pairwise crossings** (`FindChordChordCrossings`): per face, per sorted
   chord pair: re-project endpoints onto the face plane (drift from prior
   merges is zeroed by construction), map to 2D via the face's orthonormal
   in-plane basis (`FaceBasisFromNormal` - a true isometry, eps2d == eps3d),
   call `boolean2::IntersectSegments` with stableEdgeId = the chord's global
   index, lift accepted crossings back to 3D.
3. **Nearby-crossing merge** (`MergeAndPropagateCrossings`): union-find over
   raw crossings; unite when they share an incident face (hosting face plus
   both chords' face pairs) AND lie within 10 * eps. Cluster position =
   member centroid (ascending order), re-projected onto the HOST face plane
   via a member hosted there. The face gate (not a chord gate) is what
   unifies a 4-chord concurrence whose two crossings share no chord.
4. **Eager propagation**: each cluster position is tested against EVERY chord
   incident to any involved face (point-to-segment <= eps, the same t-guard),
   so a k-fold point lands on all k chords even when a pairwise intersection
   was missed.
5. **Resolve-then-allocate** (`ResolveAndThreadClusters`): per cluster, ONCE,
   symmetric across all incident chords: candidates within tolerance + eps =
   every incident chord endpoint, their existing extras, and the pass-0
   accumulator; nearest wins (ties to smallest id); else allocate. A crossing
   can therefore never thread as an endpoint id on one chord and a fresh id
   on another (the split-identity class). Each raw crossing records its
   conditioned radius (eps / sin(angle of the two chord lines), capped);
   clusters take the max over members; an allocation carries it - and a snap
   onto an existing new vert widens it - in
   `Step9Threading::newVertSnapR`, so step 9.5's original-vert snap covers
   ill-conditioned step-9 crossings exactly as it covers step 6.5's.
6. **Threading**: per chord, unify pass-0 records and resolved clusters;
   RECOMPUTE every t from the resolved position; drop step-9-added records
   whose recomputed t leaves the guarded range (pre-existing step-8 extras
   survive - they were admitted under step 8's weaker guard); id-dedup, sort
   by t (ties by id), eps/len t-dedup backstop.

Known accepted hole: when tolerance > 9 * eps, two clusters can sit within
tolerance + eps of each other yet beyond the 10 * eps merge radius and
allocate two near-coincident fresh ids; the t-dedup backstop collapses them
only when both land on one chord. Revisit if the tolerance-inflated regime
becomes a target.

## Step 9.5: arrangement-wide new-vert unification

Steps 6.5, 7, and 9 each dedup their own allocations, but the same geometric
point computed through two different frames lands up to ~10 * eps apart, and
a pair of such twins subdivides a shared sub-edge inconsistently across faces.
The unpaired sub-edges then read as open rims, whose ambient unification
collapses the cell complex (observed: rim folds merged nearly every
cell). One union-find sweep (`UnifyArrangementVerts`):

- new-new pairs unite within 10 * eps (the merge-radius philosophy);
- new verts snap onto ORIGINAL verts within max(10 * eps, the vert's recorded
  conditioned radius from step 6.5 or 9) - NEAREST wins, ties to the smallest
  id, originals before new (an id-priority pick could jump past the adjacent
  corner under a wide conditioned radius);
- consumers remap in place: chord endpoints, extras, and on-edge lists, with
  ids deduped, endpoint entries dropped, and ts RECOMPUTED from the remapped
  positions then re-sorted (a remap moves the consumed position; a stale
  order would hand the partition a crossed sub-edge sequence).

## Steps 10-11: per-face partition

`PartitionFace` partitions one face of the conforming arrangement into simple
sub-polygon cycles, CCW with respect to the face's OWN HALFEDGE WINDING:

- **The walk frame comes from the winding, not the stored faceNormal_**: on
  folded self-intersecting sheets the stored normal can be bit-exactly
  OPPOSITE the winding, which mirrors the projection and turns the
  face-on-left walk into a boundary hugger. The winding is the orientation
  the multiplicities mean.
- The face's three original edges contribute one halfedge per sub-edge
  (subdivided by the on-edge lists); incident chords contribute BOTH
  directions per sub-edge, deduped per face by undirected vert pair
  (coincident chords otherwise create exact angular ties) AND against the
  face's own boundary sub-edges (a trace chord rides its source edge by
  construction; on its host the doubled directed edge made the walk's
  exact-tie handling hes-order-sensitive - skipped, counted).
- Next-pointer rule: the smallest left turn among unvisited outgoing
  halfedges, skipping the immediate reverse UNLESS it is the sole candidate
  (the U-turn that traverses dangling-chord spurs instead of stalling).
  Closure is VERTEX ARRIVAL (the boolean2 OutEdgesToPolygons pattern).
- Each closed cycle splits at repeated vert ids (the PushSimpleLoops
  pattern); sub-3-vert loops are spurs - dropped, counted. A >= 3-vert cycle
  with EXACTLY zero projected area is a flattened spur (coincident post-merge
  positions under distinct ids, or an exactly-collinear out-and-back) -
  dropped, counted; its Newell normal would be undefined downstream.
  Tiny-but-nonzero areas are REAL slivers and pass.
- Zero-length chords (step-9 snapping can collapse v0 == v1): skipped,
  counted.

## Step 12: canonical polygon merge

`MergePolygons`: the canonical key of a cycle is the lexicographically
smallest rotation among all rotations of the cycle AND of its reversal; the
sign is +1 when the canonical form comes from the cycle as walked (CCW by the
face winding), -1 from the reversal. A simple cycle is never
rotation-equivalent to its own reversal, so the sign is well-defined. Equal
keys sum multiplicities; zero sums drop (coincident opposite-facing surfaces
cancel - the pancake killer, fed by step 6.5's conformance). Output ordered
by canonical key.

## Step 13: cells, winding, keep

The 3D lift of boolean2's twin-coupled winding filter: keep exactly the
polygons separating winding <= 0 from > 0, by GLOBAL propagation over volume
cells - never per-face classification plus repair.

SIGN CONVENTION (pinned; used identically by the seed cast and the BFS):
crossing a polygon front-to-back (front = +canonical-Newell side) changes w
by +mult - moving AGAINST the normal enters what the surface element wraps.
Pinned by `Step13CubeClassifyKeepsAllFaces`.

1. **Radial fans** (`BuildCellComplex`): per arrangement edge, incident
   polygons sorted CCW about the edge axis by their in-face direction
   (cross(normal, walk) - exact locally, convex or not), using the atan2-free
   comparator. An exact angular tie is a step-12 invariant failure
   (DEBUG_ASSERT; with conformance in place none remain on the fixtures).
   Newell normals are computed RELATIVE to each cycle's first vert - absolute
   positions cancel catastrophically for eps-thin slivers far from the origin.
2. **Wedges -> cells**: union-find over polygon sides (2p + side); between
   angularly-consecutive fan entries, unite the CCW-facing side of the
   earlier with the CW-facing side of the later. A k = 1 fan is an open
   sheet's rim: its single wedge wraps and unites the polygon's own front and
   back, as the ambient space does. Cells renumbered by smallest member key.
3. **Seed cast** (`ClassifyCells`): per connected component of the cell graph,
   a segment from P0 = arrangement bbox center + `kSeedCastDir` * 2x the bbox
   diagonal to an interior point of a component polygon, counting SIGNED
   crossings against ALL other separating polygons (other components'
   included - the true ambient winding is what a nested component cannot
   learn from its own polygons; open sheets are skipped, matching the BFS).
   Crossing tests run against a SIGNED FAN decomposition of each cycle (fan
   ears from the first vert tile any simple polygon with signed coverage; the
   kernel's direction-based steps cancel opposite-orientation overlap
   exactly; exact-zero ears contribute exactly nothing and are skipped).
   Any contact within eps of degenerate - endpoint on a surface, ear-edge
   graze, near-in-plane segment - invalidates the WHOLE cast (never skip one
   polygon and keep counting); retried on the component's next target.
   Targets are tried in DESCENDING-AREA order (ties ascending id) up to
   `kSeedCastMaxTargets` (8): big polygons' interior points sit far from
   their boundaries; slivers sort last. Triangle targets use their centroid;
   longer cycles take the largest ear of a REAL triangulation (a fan ear of a
   concave cycle can sit outside it); a triangulator throw skips to the next
   target. Exhaustion fails the classification (driver falls back).
4. **BFS**: w(back) = w(front) + mult across each separating polygon; a
   disagreement fails the classification in release too (ok = false) - the
   arrangement is not the closed surface the propagation assumes.
5. **Keep**: keep p iff IsInside(w_front) != IsInside(w_back), IsInside =
   w > 0; flip marks kept polygons whose canonical normal faces the inside
   (the emit reverses them so normals face outside). A polygon whose sides
   landed in one cell (a membrane) separates nothing and is never kept.

## Emit: twins, rings, triangulation

`CreateHalfedges` pairs halfedges by sort order, which mis-pairs when more
than two kept polygons meet at an edge, so `BuildEmitTopology` assigns the
topology explicitly:

1. **Inside-wedge twins**: at each radial fan restricted to kept polygons,
   the two flanking the same INSIDE (w > 0) wedge are twins - bare fan
   adjacency would pair across an outside wedge and weld solids that merely
   share an edge. Dropped polygons between consecutive kept ones cannot
   change the wedge's inside-ness. Odd kept fans, twin conflicts, and
   unpaired halfedges fail the topology (ok = false; driver falls back).
2. **Vertex rings**: orbits of nextAroundVert(h) = nextInPolygonCycle(twin(h))
   on the kept-polygon graph (pre-triangulation); one output vert per orbit -
   two solids touching at a vert or edge get distinct output verts (subsumes
   SplitPinchedVerts). Ring-separation argument: an orbit traces a link
   circle of the abstract surface; inside-wedge twins pair non-interleaved
   germ-sides around a fan, so one circle cannot pass through an edge's
   direction twice - every output edge carries exactly 2 halfedges
   (release-checked into ok). Rings numbered by (geometric vert, smallest
   incident kept polygon) - deterministic.
3. **Triangulation**: each kept cycle in its outward frame (its own
   relative-origin Newell basis), `manifold::Triangulate` with
   epsilon-doubling retries capped at 64x (a kept cycle can carry micro-tails
   of original verts clustered above the merge radius but below triangulable
   resolution; their ring ids are topologically pinned, so the cycle cannot
   be simplified - widening epsilon moves the tail into the triangulator's
   own degenerate class, and the zero-area output tris collapse at Manifold
   construction; far beyond 64x the CCW check would pass CW triangles over
   real geometry). Output MeshGL64: numProp = 3 - non-position properties
   are NOT preserved.

## Driver gate

Return the input unless (a) `out.Status() == NoError`, (b) `out.Volume() > 0`
for non-empty input (NaN fails too), (c) pierce-monotonicity:
`CheckSelfIntersection(out) <= CheckSelfIntersection(input)` (input count
computed once, early), and (d) the folded-shell volume gate below. NO
volume-ratio tripwire: the input volume is the winding-WEIGHTED integral - an
overlap lobe at w = k counts k times - so the legitimate output/input ratio
is (2-f)/(2+f) for overlap fraction f (1/3 at full overlap); any constant
bound vetoes correct outputs on exactly the heavy-overlap inputs the feature
targets. Gate trips are the expected fallback for adversarial inputs, not
asserts; internal invariants (BFS disagreement, odd kept fans, unpaired
halfedges) keep their DEBUG_ASSERTs and fail closed in release.

**Folded-shell volume gate** (`FoldedCellsEncloseVolume`, pre-emit): polygons
whose front and back cells united (front == back, e.g. across a k = 1 rim or
through a mis-ordered near-tangent radial fan) are dropped by the keep rule.
That is correct for membranes - enclosed volume below area x thickness, and
post-merge thickness is at most the 10 eps unification radius - but a
tangent-degenerate contact can fold a CLOSED shell's two cells together, and
dropping that fold silently deletes the shell (observed on the hull fixture:
whole disjoint hulls, most of the material). Per EDGE-CONNECTED component
of folded polygons (grouped by fold cell + shared undirected edge - per
cell, or vert-connected grouping, would let a positive shell and an
inverted twin that share a cell or merely touch at a snapped vert net
their signed volumes to nothing), sum the mult-weighted signed volume
(tetra fan anchored
at the component's own centroid: origin-independent for a closed set, near
zero for an open sheet) and fall back if any component exceeds
`area x kFoldedVolumePerAreaEps (= 100) x eps` - 10x headroom over the
thickest legitimate membrane, orders of magnitude below a real shell.

## Determinism constraints (pinned)

- Chord enumeration: `newEdges` order (map keyed by sorted tri pair; serial
  collider traversal). Trace pairs processed in ascending (faceA, faceB).
- Pair iteration: sorted within each face; faces ascending. Union-find unites
  in sorted order; centroids sum ascending.
- Candidate resolution: nearest, tie smallest id. Rings: (vert, smallest
  kept polygon). Cells: smallest member key. Seed targets: descending area,
  tie ascending id.
- `kSeedCastDir`: a fixed generic unit vector (no axis alignment) so casts
  into typical axis-aligned inputs avoid grazes on the first try.
- No `std::` trig in any decision path (none needed - bucket + cross sign).

## Validation

(Coverage is described by KIND, not count - counts and exact observed
numbers go stale instantly; the suite is the source of truth.)

`test/manifold_test.cpp`:
- `OverlapRemoval.*` unit tests per stage: step-9 kernel/merge/resolution/
  threading (including the split-identity, face-gate, conditioned-radius,
  and re-sort pins); partition (X-crossing, dangling-spur, coincident-dedup,
  zero-length, boundary-riding, chord-order invariance, stored-normal
  inversion); step-12 canonicalization and cancellation; cell complex
  (tetra, bipyramid-with-internal-face); winding classification (reversed
  representation, nested cubes, membrane-across-the-cast, concave seed
  target, seed-cast exhaustion); emit topology (the book fixture pinning
  twin pairing + ring splitting, odd-fan fail-closed); step 6.5 (pancake
  quad, coplanar neighbors, plane-gate offset, full-through cuts, on-edge
  additions, snapped-endpoint edge subdivision); step 9.5 (unification,
  nearest-original, per-vert radius carry including the composed
  step-9 -> 9.5 handoff, t recompute + re-sort); step-6 event-identity
  snap band; step-1 merge displacement and fixed-point convergence
  ordering; fold-gate membrane-pass / bent-open-trip / opposite-shell
  non-cancellation arms; interior-island detection (the free-island,
  pinched-loop, and proper-double-crossing cases).
- `Manifold.RemoveSelfIntersections*` feature tests: API smoke; clean-input
  and Boolean-result passthrough (bit-identical); the hull fixture (the
  trimaran fold class - pins the folded-shell gate's bit-identical
  fallback, see Known limitations); the ovoid dense-sliver fixture (falls
  back bit-identically - the outcome is pinned, not which internal
  guard fires); the interior-island stamp (bit-identical fallback; the
  identity assert discriminates gate removal, which would emit a
  wrong-but-valid mesh); the merge-displacement tolerance fixture (the
  eps-chain-strip weld; mutation-checked against the formula); empty
  input; idempotence (fallback fixed point + success-path monotonicity);
  determinism; far-from-origin (the same trimaran at scale: strict pierce
  reduction, all components preserved, volume preserved to the test's
  bar); glued boxes (equal-face early-exit bit-identical;
  smaller-on-larger welds, winding-faithfully, to one component).

## Known limitations

1. **The tangent-degenerate contact class** (the hull fixture). A
   shallow-incidence edge piercing two eps-SEPARATED coplanar sheets produces
   twin step-7 events that are geometrically REAL distinct points
   ~eps/sin(incidence) apart (observed at tens of eps), beside an original
   corner.
   The exact arrangement has a micro-triangle facet there that per-face FP
   partitions cannot consistently produce. Every snap policy beyond ~10 eps
   (fixed wide anchors and conditioned isotropic/anisotropic step-7 snap
   variants alike) traded the twin-rim holes for MORE eps-overlap pierces
   and was reverted: moving geometry tens of eps deforms kept triangles whose
   neighbors did not move with them. The class shows up at two severities
   on the hull fixture (multiple disjoint hulls grazed by one mask):
   - **Micro-facet pierce residue**: input pierces run thousands of eps
     deep; a successful rebuild leaves residual pierces tens of eps
     deep - above the 10x-eps output tolerance, inside the conditioned band
     of that corner.
   - **Folded shells**: the grazing contacts on the outrigger hulls leave
     vert clusters spread a few times past the 10 eps unification radius
     (k = 1 rim chains) and near-tangent k = 4 radial fans; either folds the
     entire shell's front cell onto its back cell, so every polygon of that
     hull reads front == back and the keep rule would silently delete the
     whole component. The folded-shell volume gate (Driver gate) detects
     this and falls back to the input bit-identically.
   Far from the origin the same geometry succeeds (strict pierce
   reduction, every component preserved): the scale-derived eps absorbs the
   clusters and the residual sits within a few eps - INSIDE its working
   band. Closing the class soundly (resolving the fixture at origin scale
   to zero pierces with every component preserved) needs
   exact/extended-precision
   local predicates (the family Emmett deferred).
2. **Dense slivers** (the ovoid class): the arrangement is not a closed
   surface after FP partitioning; the BFS disagreement guard detects it and
   the pipeline falls back to the input, monotonic by construction.
3. **Collinear-overlapping chords**: the partition's per-face dedup absorbs
   exact-id duplicates; eps-distinct near-collinear duplicates remain (the
   doubled-cut class; measure-zero for generic inputs).
4. **tolerance > 9 eps**: the step-9 fresh-id hole above.
5. **Welding**: coincident interior walls separating w = 1|1 drop (see step
   6.5's semantics note) - winding-faithful, documented, pinned by fixture.
6. **Properties**: non-position properties are not preserved (numProp = 3).
7. Output tolerance widens to cover the pipeline's applied movements (the
   10-eps floor plus the measured merge/unification displacements); see the
   eps contract's OUTPUT tolerance bullet for the exact formula.
8. **Interior-island stamps** (gated fail-closed): a shell whose
   intersection curve with a face does not properly cross its boundary - a
   "stamp" footprint strictly interior to the face, or one PINCHED onto a
   single boundary vert (reachable when a corner snap lands a loop endpoint
   on a boundary vert) - creates a chord loop whose surrounding region is a
   (possibly pinched) annulus, not representable as simple cycles. The
   partition would emit the loop in both orientations (canceled at step 12)
   and the bare boundary - silently erasing the cut and misclassifying the
   stamping shell as nested. The partition gates this per chord-only
   component: a cycle-bearing component with fewer than two distinct
   boundary attachments (`FacePartition::interiorIslandVerts`) fails the
   run closed, bit-identically. Deeper pinched compositions (e.g. nested
   loops bridged through one attachment) may still pass the count and fall
   to the downstream gates. Resolving the class needs hole-aware faces
   (bridge edges), a known arrangement technique deliberately out of
   scope.
9. **Fold-gate residual**: two opposite-orientation folded shells that share
   a REAL undirected edge (same two snapped vert ids, without the coincident
   opposite polygons step 12 would cancel) are edge-connected into one
   component and can net their signed volumes under the threshold. Requires
   eps-coincident input edge geometry between oppositely wound shells that
   BOTH fold into one cell - adversarial-construction territory, documented
   rather than guarded.

## Relationship to #289 and design history

Steps 1-13 map onto Emmett Lalish's 13-step sketch (the 2024-05-14 comment;
NOT the abandoned 2023 serial-seam sketch). The winding classification is
step 13; the in-plane crossing kernel is boolean2's `IntersectSegments`
(production since #1722/#1751); the angular ordering, the canonical merge,
and the winding filter are the boolean2 patterns lifted to 3D.

The design went through staged adversarial review, multiple rounds per
stage: step 9 (architecture; kernel contract; the split-identity
resolve-then-allocate rule; t-recompute ordering; post-implementation code
bugs fixed red-first). Steps 10-13 (the trim pre-pass replaced by the
U-turn + split-at-repeated-vert production pattern; per-(edge,pair) vert
duplication replaced by twin assignment + ring extraction; the volume
tripwire dropped with worked math). Step 6.5 (driver wiring, two-frame id
divergence, snap-rule attribution, early-exit and welding semantics, grazing
inflation; a max-endpoint-margin counterproposal was declined on the
concavity argument - it would reject full-through cuts). A whole-branch
post-implementation review found further fixes (stale ts after remap; the
cast counting membrane crossings; the cross-pair conditioned gap; the retry
cap; release-mode classification failure; an inverted sign statement in this
doc's ancestor), the worst red-first. An eps-propagation audit then plumbed
the driver's eps into the cast and made the output tolerance claim honest.

The empirical record behind the residue analysis (the snap-radius
experiments and their outcomes) is in the git history of
`docs/Steps10to13Design.md` - a working doc deleted when its content was
consolidated here; the history remains reachable at the deleting commit.

Repeated whole-branch dual-track review passes (Claude + Codex lanes:
correctness, numerical robustness, tests/docs/API, style and refactoring)
after consolidation produced the remaining hardening: the conditioned-radius
carry through step 9, source-gated trace dedup, nearest-original
unification, the measured-displacement tolerance claim, the single-leaf BVH
guard, the snapped-endpoint on-edge propagation, the position-fixed-point
merge convergence, and - found while chasing the hull fixture's volume red -
the folded-shell volume gate above.
