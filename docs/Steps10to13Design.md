# Steps 10-13 design v5: partition, merge, global winding, emit

v5 = v4 + the increment-v driver findings: a new COPLANAR CONFORMANCE pass
(step 6.5) and seed-target hardening - see "v5 delta" at the end. Both came
out of wiring the full driver: the hull fixture exposed coplanar pancake
regions that steps 4-12 as sketched cannot cancel, and the far-from-origin
fixture exhausted the 3-target seed retry on sliver ears.

v4 = v3 + round 3: the closure guard corrected to vertex-arrival (the
boolean2 production pattern); the ring orbit operator pinned on the POLYGON
graph; the book-case ring-separation argument written down (the round-3
HIGH) + an emit assert; three small framing fixes. Round-3 findings were
narrowing spec-wording items - design declared CLEAR for implementation.

Working doc for the faithful #289 rebuild (companion to docs/Step9Design.md,
complete and implemented). Not for upstream; lives in docs/ so review
subagents can read it. Clean up before any PR.

v3 folds in review round 2. Major changes: the trim pre-pass is REPLACED by
boolean2's production pattern (U-turn at sole-candidate verts + split-at-
repeated-vert + drop zero-area loops), which dissolves round 2's ordering bug
and the boundary-test spec gap; the emit's vert duplication is REPLACED by
explicit twin assignment + vertex-ring extraction (round 2 showed per-
(edge,pair) duplication does not compose at verts); the seed cast aborts (not
skips) on graze; the volume tripwire is DROPPED with a worked rationale (any
constant ratio bound against the winding-weighted input volume rejects
legitimate heavy-overlap outputs).

## Scope

Rebuild Emmett's #289 steps 10-13 on top of the completed step-9 arrangement,
then wire the full driver into `RunOverlapRemovalImpl`. NO recovery scaffold.
The manifold result comes from GLOBAL, twin-coupled winding propagation (the
boolean2 winding filter lifted to 3D), never from independent per-face
classification plus repair. Internal invariant violations are DEBUG_ASSERTs.

## Inputs (what steps 1-9 deliver)

- `onEdgeLists` (each original edge's on-edge verts by t). ORDERING PIN:
  `PropagateNewVertsToOnEdgeLists` runs right after `GenerateChordEdges`
  (step 7), before step 8. (Its header comment says "step 10" - stale; fix.)
- `chordsWithExtras` threaded by step 9: a CONFORMING arrangement (no two
  sub-edges cross except at shared vert ids).
- `newVertPositions` + `GetPos3`; `chordsByFace` from the pre-threading
  `newEdges`.
- `Step9Threading::crossings` is diagnostics only; 10-13 read
  `chordsWithExtras`.
- Tolerance PIN: `max(impl.tolerance_, eps)` from the post-`MergeVertsEps`
  impl.
- Zero-length chords (step-9 snapping can yield v0 == v1): filtered at
  partition entry (skipped, counted).

## Step 10-11: per-triangle halfedge partition

Per Emmett: per face, the three original edges (subdivided by onEdgeLists)
contribute one halfedge per sub-edge (CCW); every incident chord contributes
BOTH directions per sub-edge. Join into polygons by cyclic ordering around
each vert. Reimplementation of the deleted walker (reference:
`backup/rsi-recovery-scaffold`) with these changes:

a. **Sub-edge dedup.** Per face, chord sub-edges dedup by undirected
   (min,max) vert pair: coincident chords (the collinear-overlap class)
   otherwise create exact angular ties that mangle the sort. One undirected
   sub-edge -> one forward+backward halfedge pair.
b. **U-turn at sole-candidate verts (replaces v2's trim pass).** The
   next-pointer rule skips the immediate reverse halfedge UNLESS it is the
   only candidate at the vert (chord-degree-1 interior endpoints: dangling
   chords from dropped sibling pairs, deduped slit ends). The walk then
   traverses the spur out-and-back inside the surrounding cycle and CANNOT
   stall; a stall (next == -1) is a hard DEBUG_ASSERT.
c. **Split at repeated verts (production boolean2's `PushSimpleLoops`
   pattern).** Post-walk, split each closed cycle at repeated vert ids into
   simple sub-loops; loops with < 3 verts or zero area (the spur 2-gons)
   are dropped and counted. Output cycles are SIMPLE by construction -
   the emit's Triangulate precondition.
   KNOWN LIMITATION (documented, accepted): an exactly-coincident
   same-orientation doubled cut (two distinct shells crossing a face along
   the identical segment) survives as a spur, so that face is not split
   there and the cells across that segment are not separated. Requires
   measure-zero geometry that step-1 vert merging almost always collapses;
   the step-12 cancellation handles the opposite-orientation twin case.
d. **2D frame**: the step-9 orthonormal basis (true isometry, re-projected
   endpoints); angular ordering via boolean2's atan2-free comparator
   (half-plane bucket + cross sign). No trig (the old `std::atan2` violated
   the determinism rule).
e. **Closure guard: vertex arrival** (round 3): the walk closes when the
   next destination vert equals the start vert (`destV == startV`, the
   boolean2 `OutEdgesToPolygons` pattern) - NOT "re-select the start
   halfedge", which is impossible once visited[] marks it.

Output: per face, simple sub-polygon cycles, CCW w.r.t. the face normal.

## Step 12: canonical polygon merge with signed multiplicity

A NEW cycle canonicalizer (boolean2 `Canonicalize` is the key+sum PATTERN
only):

- Canonical key: lexicographically-smallest rotation among all rotations of
  the cycle AND of its reversal.
- SIGN RULE (independent of step 13): +1 if the canonical form is a rotation
  of the cycle as walked (CCW w.r.t. its face normal), -1 if of the
  reversal.
- Merge: map<key, summed multiplicity>; zero-sum entries drop (doubled
  surfaces with opposite orientation cancel - the projection convention
  makes their cycles reversals of each other by construction). This
  cancellation is also what grounds the emit's antiparallel-pairing proof:
  after step 12, two coplanar kept polygons cannot share an edge.
- Cycle keying (not Emmett's vert-set) is deliberate: distinct same-vert-set
  polygons must not merge; identical coplanar polygons share threaded vert
  ids and so share the key.

Output: `vector<MergedPolygon { cycle; mult; face }>`.

## Step 13: global winding classification (the heart)

Keep exactly the polygons separating winding <= 0 from > 0, by global
propagation over the arrangement's volume cells - the 3D lift of boolean2's
`FilterByWindingHalfedges`.

1. **Half-faces.** Each merged polygon has front (+normal) and back. SIGN
   CONVENTION (pinned, used identically by the seed cast and the BFS):
   crossing front-to-back changes w by -mult; back-to-front by +mult.
2. **Radial edge fans.** Per arrangement edge (undirected vert pair):
   collect incident polygon-sides; sort radially by the in-face direction
   (edge-perpendicular into the polygon; locally well-defined regardless of
   convexity) projected perpendicular to the edge, ordered by the atan2-free
   comparator in a deterministic frame (step-9 basis seeded from the edge
   direction). An EXACT angular tie (same polygon plane through the edge -
   not merely near-parallel distinct planes, which sort normally by their
   distinct angles) = step-12 invariant failure: DEBUG_ASSERT
   (ascending-polygon-id fallback only as release determinism insurance).
   Fan-loop check: signed mult changes around a full fan sum to zero, else
   DEBUG_ASSERT.
3. **Wedges -> cells.** Wedge = the region between angularly-consecutive
   half-faces around an edge. Unite: (a) wedges around each fan with the
   half-faces they touch; (b) all wedges along one SIDE of one polygon are
   one cell - sound because in a conforming arrangement no polygon
   subdivides the open region over another polygon's interior without
   sharing a chord edge (an intersecting plane's trace IS a chord; a
   parallel polygon does not subdivide the slab). Union-find over fan
   slots, deterministically keyed; cell ids = smallest member slot key.
4. **Seed + propagate.** Per connected component of the cell graph,
   SEGMENT-CAST seed (an extreme-vertex seed is wrong for nested
   components):
   - Target: an interior point of a chosen polygon Q = the centroid of the
     first ear of Q's triangulation (a concave polygon's vert-centroid can
     fall outside it).
   - Source: P0 = arrangement bbox center + kSeedCastDir * 2x bbox diagonal
     (a fixed generic unit direction; it is a segment cast, not a ray).
   - Count signed crossings of the open segment (P0 -> target) against ALL
     merged polygons EXCLUDING Q (counting Q would measure the far side).
     The arrival winding seeds the cell adjacent to Q's P0-FACING half-face:
     front if dot(target - P0, normal(Q)) < 0, else back.
   - GRAZE RULE (round 2): if ANY polygon's crossing test is within eps of
     degenerate (edge/vert graze, near-parallel plane), the ENTIRE cast is
     invalid - abort and retry with the next polygon of the component (up
     to 3 targets), then DEBUG_ASSERT + return input. Never skip a single
     polygon and keep counting: a skipped +-mult corrupts the seed.
     Crossing tests are Moller-Trumbore against the polygons' ear
     triangulations with the eps graze-reject; retry replaces SoS for the
     measure-zero cases.
   Then BFS across half-faces with the pinned sign convention; a BFS
   disagreement is a DEBUG_ASSERT.
5. **Keep.** Keep p iff IsInside(w_front) != IsInside(w_back), IsInside(w) =
   w > 0. Orient each kept polygon with its normal toward the outside
   (w <= 0) cell.

## Emit (round 2: twins + vertex rings replace vert duplication)

`CreateHalfedges` pairs halfedges by sort order; >2 kept polygons at one
edge would pair wrongly, and per-(edge,pair) vert duplication does not
compose at verts (round 2). So the emit assigns topology EXPLICITLY:

1. **Twin assignment.** For every kept polygon edge: at the edge's radial
   fan, restrict to KEPT polygons and pair each with the kept polygon
   flanking the SAME INSIDE (w > 0) WEDGE - precisely: walking the fan, an
   inside wedge of some region lies between two kept half-faces; those two
   are twins. (Bare fan-adjacency is ambiguous - adjacent-around-an-OUTSIDE-
   wedge pairs would invert orientation.) With exactly 2 kept this
   degenerates to the normal case. Two regions meeting 2+2 at an edge pair
   within their own regions automatically (each inside wedge belongs to one
   region). Outward orientation makes twins traverse the shared edge
   antiparallel (the standard closed-surface argument; grounded by step-12
   coplanar cancellation).
2. **Vertex rings.** The twin assignment defines an abstract closed
   surface over the kept POLYGONS (pre-triangulation). Orbit operator
   (round-3 pin): `nextAroundVert(h) = nextInPolygonCycle(twin(h))`, where
   `nextInPolygonCycle` is the CCW successor within the kept polygon's
   cycle - this runs on the polygon graph, NOT the post-triangulation
   halfedge mesh; ring ids are assigned BEFORE triangulation and the
   triangulation inherits them. Output vert ids = one id per orbit (ring)
   at each geometric vert, positions from GetPos3.
   RING-SEPARATION ARGUMENT (round-3 HIGH): two twin-pairs at one edge E
   cannot share rings at an endpoint u. An orbit traces a link circle of
   the abstract surface around u; each passage through E's direction
   consumes one germ-side pair, and inside-wedge twins pair
   NON-INTERLEAVED germ-sides around the fan (each pair flanks its own
   disjoint inside wedge), so one circle cannot pass through E's direction
   twice - the two pairs lie on two distinct rings, and every output edge
   carries exactly 2 halfedges. Belt-and-suspenders: the emit
   DEBUG_ASSERTs (and release-checks, feeding the gate) that no output
   (vertA, vertB) edge exceeds 2 halfedges.
   This subsumes SplitPinchedVerts (edge- and vertex-level pinches
   uniformly) and makes `CreateHalfedges`' sort-order pairing reconstruct
   exactly the assigned twins. Determinism: rings ordered by (geometric
   vert id, smallest incident kept-polygon id); regions/cells by smallest
   member key.
3. **Triangulation.** Project each kept polygon in its KEPT orientation's
   frame (reverse the cycle first if step 13 flipped it) using its face's
   step-9 basis; `manifold::Triangulate` on the projected SimplePolygon
   (simple by construction - partition rule c); emit triangles
   winding-consistent with the kept normal, vert ids from the ring map.
4. Build MeshGL64 (numProp = 3; non-position properties are NOT preserved -
   document), construct `Manifold`.

## Driver + contract gate

`RunOverlapRemovalImpl`:
- steps 1-7 -> `PropagateNewVertsToOnEdgeLists` -> step 8 -> EARLY-EXIT: if
  `ChordEdges::newEdges` is empty return `input` (covers both the
  no-intersections case and the all-pairs-dropped case; the doc note "may
  be bit-identical" covers both) -> step 9 -> 10-11 -> 12 -> 13 -> emit.
- GATE (thin, final): return `input` unless (a) `out.Status() == NoError`,
  (b) out is non-empty with `outVol > 0` for non-empty input (NaN volume
  fails this comparison too; note the acknowledged semantic edge: a
  legitimately zero-volume flat-sheet result also falls back - accepted),
  and (c) pierce-monotonicity `CheckSelfIntersection(out) <=
  CheckSelfIntersection(input)` (compute the input count once, early).
  NO volume-ratio tripwire (round 2, worked math): `input.Volume()` is the
  winding-WEIGHTED integral - an overlap lobe at w=k counts k times - while
  the output measures the w>0 region once, so the legitimate ratio is
  (2-f)/(2+f) for overlap fraction f (-> 1/3 at full overlap, -> 1/k for
  k-fold), and inverted lobes can push the ratio above 1; any constant bound
  vetoes correct outputs on exactly the heavy-overlap inputs the feature
  targets. Wrong-seed risk is owned by the debug-mode BFS/fan asserts and
  the nested/book fixtures, not a release heuristic.
  The gate arms do NOT DEBUG_ASSERT (changed from round 3 at
  implementation): a gate trip is the expected fallback for adversarial
  inputs (the dense-sliver class), not an invariant violation - asserting
  would fire on known-fallback fixtures in every debug run. Internal
  invariants (BFS disagreement, odd kept fans, unpaired halfedges, ...)
  keep their DEBUG_ASSERTs.
- Serial-first; `// TODO: parallelize` markers only.

## Doc debt at landing

- `manifold.h` RemoveSelfIntersections doc: drop sign-flip/merged-form/drift
  language; describe the thin gate; clean-input early-exit may be
  bit-identical; non-position properties not preserved; epsilon-merged
  positions caveat stays.
- `docs/OverlapRemoval.md`: rewrite the Algorithm Outline (old steps 11-13),
  the Guarantee list, and DELETE the stall/cap-walker Known-Limitations text;
  add the coincident-doubled-cut limitation (partition rule c).
- `overlap_removal_internal.h`: fix the stale "Step 10" comment on
  `PropagateNewVertsToOnEdgeLists`; remove still-dangling old-walker decls
  alongside increment i.

## Increments (TDD, in order)

i.   Partition (10-11) with dedup, U-turn, split-at-repeated-vert. Fixtures:
     real `Manifold::Impl` from MeshGL64 + synthetic threaded chords. Pin:
     X-crossing face -> 4 polygons; dangling chord -> spur split out and
     dropped (counted), surrounding polygon simple; coincident duplicate
     sub-edge dedups; zero-length chord skipped; no stalls.
ii.  Step 12: canonicalization, the sign rule, opposite-pair cancellation,
     multiplicity sums.
iii. Radial fans + wedge/cell union on a hand-built two-cell arrangement:
     fan order, fan-loop zero-sum, cell count, deterministic keys.
iv.  Seed + propagate + keep end-to-end: single cube (all faces kept,
     w 0|1 - the clean-passthrough pin); nested cubes (inner face separates
     w=1|2 -> NOT kept; seed must measure the inner component's true outer
     winding via the cast); a 4-kept-at-an-edge book fixture pinning the
     inside-wedge twin pairing + ring extraction.
v.   Driver wiring + emit; restore the two relaxed feature tests (HullMask
     EXPECT_EQ 0; FarFromOrigin EXPECT_LT) - the acceptance gate for the
     whole rebuild.

## Reuse

- Deleted walker (backup branch) as reference only.
- boolean2 `winding_filter.cpp` (seed/propagate/keep + the atan2-free
  comparator) and `PushSimpleLoops` (`boolean2.cpp`) for the repeated-vert
  split pattern.
- `manifold::Triangulate`, `GetPos3`; hoist the step-9 basis helper.
- `AnalyzeSelfMesh` is NOT used by 10-13.

## Known risks (current)

- The twin/ring emit is new ground; increment iv's book fixture is its proof
  obligation.
- Seed retry exhaustion returns input (counted).
- The coincident-doubled-cut limitation (partition rule c) leaves rule
  3(b)'s premise unmet on measure-zero inputs; documented.
- Scale: FarFromOrigin bar is strict reduction; plane re-projection
  cancellation at 1e4 is the known precision tax.
- Performance unmeasured; near-linear except the per-component seed cast.

## v5 delta: coplanar conformance (step 6.5) + seed-target hardening

### What the driver run exposed

Wiring the full driver (increment v) and running the hull fixture
(hull-body Subtract hull-mask, ~31 pierces) surfaced two design gaps.
Empirical landscape with the gaps waved through: every other feature test
passes (clean inputs early-exit bit-identical; the ovoid dense-sliver
fixture falls back cleanly, as documented); the hull fixture and its 1e4
translation are the only failures, and each isolates one gap.

GAP 1 - coplanar pancakes. Step 1's eps-merge flattens Boolean3's SoS
sliver wedges into zero-volume pancakes: coplanar overlapping face pairs
with opposite orientations. That is BY DESIGN - the merge manufactures the
coincidences step 12's signed-multiplicity cancellation consumes. But the
sketch's step 6 finds only TRANSVERSAL edge-tri pierces: coplanar pairs
produce no chords, the two sheets' partitions never conform in-plane,
step 12's equal-cycle cancellation cannot fire, and BuildCellComplex hits
EXACT angular ties (bit-identical in-face directions at shared fan edges -
observed directly on the fixture: a big tri on one sheet against a
3-sub-tri Steiner triangulation of the same region on the other, sharing
all four verts, surviving with unbalanced mults). Letting the tie through
on polygon-id order corrupts the wedge unions: downstream the winding BFS
disagrees and kept cycles emit non-simple (Triangulate's CCW check
throws). The tie assert is correct; the arrangement is what is incomplete.

GAP 2 - seed-target slivers. At 1e4 the seed cast exhausted its 3 targets:
ascending-polygon-id order picked eps-thin slivers whose first-ear
centroids sit within the graze margin of their own boundaries, so every
cast aborted. (Distinct from the documented re-projection precision tax;
this is target SELECTION, not cast math.)

### Step 6.5: coplanar trace chords (v6: post-review revision)

Two adversarial review lanes (geometry; integration/numerics) returned
NOT CLEAR on the v5 sketch. v6 pins every blocked rule. Review traceback:
[G1] driver wiring, [G2]/[I-A4] two-frame id divergence, [I-A1] snap-rule
misattribution, [I-A2] on-edge injection path, [I-B1] early-exit/welding
semantics + coplanar-neighbor noise, [I-C1] grazing-clip inflation,
[G4] detect thresholds, [G7] rule-3 endpoint ambiguity, [G9] rim fans.

After step 7 (GenerateChordEdges), BEFORE PropagateNewVertsToOnEdgeLists
and step 8, run a conformance pass over coplanar overlapping face pairs:

1. **Detect (broad + plane gate).** Fresh BVH over tri boxes (the
   SortedBVH ritual; serial, once per run); for each candidate pair, the
   plane gate: ALL SIX verts within eps of the LARGER-AREA face's plane
   (one-sided against the larger face only - a near-zero-area sliver's
   own plane is noise [G4]; eps is a length). Near-coplanar-but-tilted
   pairs fail and keep their transversal step-6 chords - no double
   handling.
2. **Single-frame clip [G2, I-A4].** ALL geometry for a pair is computed
   in ONE frame: the lower-face-id face's FaceBasisFromNormal frame.
   Clip each of B's 3 edges against A's projected triangle AND each of
   A's 3 edges against B's projected triangle in that same frame (convex
   clip: at most one interval per edge). The same geometric crossing is
   therefore computed exactly once and shared by both clip directions -
   the two-frame id-divergence failure mode is removed by construction.
3. **Interval qualification [I-B1, G4].** An interval becomes a trace
   chord only if (a) its length > eps AND (b) its midpoint lies strictly
   interior to the OTHER face with in-plane edge-distance margin > eps.
   (b) rejects boundary-coincident intervals: coplanar NEIGHBORS (any
   flat region of any mesh - e.g. a cube face's two tris) produce only
   boundary-riding intervals and emit NOTHING, so clean flat meshes
   still take the early-exit; equal-size face-glued solids likewise
   emit nothing and pass through bit-identical. A qualifying interval
   exists only where one face's boundary genuinely crosses the other's
   interior - the pancake class.
   MIDPOINT, not endpoints (re-review round 2 adjudication): a
   full-through cut - B's edge entering AND exiting A, the generic
   overlap case - has BOTH endpoints on A's boundary (endpoint margins
   ~0) yet an interior midpoint, so an endpoint-margin predicate would
   reject exactly the cuts conformance needs most. Distance-to-boundary
   is CONCAVE along a segment inside a convex face, so the midpoint
   margin >= half the deepest penetration: the midpoint rule keeps
   every cut deeper than 2 * eps and rejects only the FP-degenerate
   dip band the rule-4 grazing guard already skips.
4. **Endpoint ids [I-A1, G7].** Pin: interval endpoints that are
   original verts (a vert of one tri inside or on the other) use their
   ids directly. A genuine crossing endpoint snaps to the nearest of
   the pair's six corner verts within tolerance + eps (nearest, ties to
   smallest id) - this is the STEP-9 resolve-then-allocate convention,
   chosen deliberately: trace verts feed straight into step-9
   resolution, so they must obey its radius, NOT step 7's bare-eps
   first-found (which is the convention only for step-7's own
   pierce-event dedup). Unsnapped endpoints allocate new verts, then
   dedup new-to-new across ALL step-6.5 verts at eps, first-found, in
   deterministic pair order (the step-7 new-to-new convention, named as
   such). Lift back to 3D in the shared frame; clamp the crossing t and
   REJECT (skip the interval) any endpoint whose lifted position is
   farther than eps from either original 3D edge - the near-grazing
   inflation guard [I-C1]; a skipped degenerate interval costs only
   conformance we provably cannot compute.
5. **Trace chords.** Each qualifying interval emits
   `PiercedNewEdge{v0, v1, A, B}` (it lies on both coplanar faces; the
   existing chord model applies verbatim). Dedup by sorted endpoint ids
   within the pair; ACROSS pairs the partition's per-face undirected
   dedup (rule a) absorbs geometric repeats - there is NO pre-partition
   global chord dedup [G3].
6. **Boundary conformance [I-A2, G7].** Only NEW crossing endpoints
   that lie on an original mesh edge (within eps, by construction of
   the clip) generate on-edge insertions; original-vert endpoints are
   already endpoints and are NOT inserted. A crossing landing on edges
   of BOTH faces (the X case: B's edge crossing A's edge) gets one
   record PER EDGE - the two clip directions contribute one each, with
   the same vert id and per-edge t. The pass returns explicit
   (edgeIdx, vertId, t) additions - t recomputed from 3D positions -
   and the driver applies them via a NEW sibling helper
   `AddVertsToOnEdgeLists` (id-dedup + per-edge t re-sort, the same
   internals as PropagateNewVertsToOnEdgeLists, which structurally
   cannot carry them: its interface is parallel to etIsects).
7. **Compose [G1].** The driver appends trace chords to
   `ChordEdges::newEdges` and the new positions to
   `ChordEdges::newVertPositions` BEFORE the early-exit, BEFORE
   GroupChordsByFace, and BEFORE step 8:
   - EARLY-EXIT moves after this pass: return input iff the COMBINED
     newEdges is empty (a pancake-only defect has no transversal
     chords at all).
   - `chordsByFace` = GroupChordsByFace(COMBINED list) - one grouping
     feeds step 8 extras, step 9 (trace-vs-regular and trace-vs-trace
     crossings resolve there [G5]), and the partition. No other
     consumer changes.

After conformance both sheets partition into IDENTICAL sub-polygon
geometry over the overlap region (same ids by rules 2+4), step 12
cancels the doubled area exactly (equal canonical cycles, opposite
signs), and the angular-tie DEBUG_ASSERT stays - a remaining tie is
again a real invariant failure.

Semantics pin [I-B1]: a coincident interior wall separating winding
1|1 (differently-sized solids glued face to face, after conformance)
DROPS by the step-13 keep rule - the output is the winding-faithful
welded solid. This is correct #289 behavior, not a regression: the
input's w > 0 region IS one solid. Equal-size glued faces never reach
the pipeline (rule 3 emits nothing; early-exit). Documented with a
fixture.

Rim note [G9]: cancelling a pancake leaves its rim edges with a
single surviving non-coplanar polygon - the k = 1 open-sheet fan whose
existing rim rule (unite own front and back) is exactly the ambient
behavior wanted: the BFS crosses the survivor with its own mult and
the vanished pancake imposes no constraint.

Absorbed edge classes (unchanged from v5, restated tighter):
- Boundary-coincident intervals: now REJECTED at rule 3 (not emitted,
  rather than emitted-and-absorbed).
- A fully inside B: A's edges qualify (interior midpoints), B's edges
  do not reach A's interior; B's partition carves the A-shaped hole;
  cancellation proceeds.
- 3+ stacked sheets: pairwise traces compose; cross-pair crossings are
  step-9 chord-chord crossings on the shared face; multiplicities sum
  per step 12. A trace endpoint allocated near a THIRD sheet's corner
  (outside its own pair's six-corner snap set) resolves via step 8 ->
  step 9: the corner is an on-tri vert of the host face, step 8
  threads it onto the trace chord as an extra, and step-9 resolution
  snaps the eps-close endpoint to it - this RELIES on the rule-7
  ordering (step 8 runs over trace chords before step 9).
- Collinear overlapping edges: intervals ride the boundary - rejected
  by rule 3; eps-distinct near-collinear duplicates remain the
  documented doubled-cut limitation.

Cost: detect is one BVH pass + plane gates; clips run only on pairs
surviving the gate; everything is serial with `// TODO: parallelize`.

### Seed-target hardening (v6)

Replace ascending-polygon-id target order with DESCENDING-AREA order
(area = |relative-origin Newell| of the canonical cycle, computed in
deterministic cycle order; ties by ascending polygon id), and raise the
retry budget from 3 to `kSeedCastMaxTargets = 8` (named constant in the
.cpp numeric-defaults block [I-D1]). Big polygons have ear centroids far
from their boundaries, so the first target almost always casts cleanly;
slivers sort last instead of first. Deterministic given identical input
order [I-D2]; no cast-math change. All-8-graze still falls back to the
input (counted) - accepted.

### Increment plan (red-first; v6 reviewed)

v-2a. CoplanarTraceChords unit: overlapping coplanar tri pair (opposite
      orientation) -> expected trace chords with single-frame shared
      ids; the contained-vert, shared-edge (emits nothing), coplanar-
      neighbor (emits nothing), and near-grazing-reject cases.
v-2b. Driver pancake fixture: Compose(cube, zero-volume pancake whose
      two sheets triangulate the same quad with DIFFERENT diagonals) ->
      output is the cube alone (pancake cancels in step 12), 0 pierces.
      Plus the glued-boxes fixtures: equal faces -> bit-identical
      early-exit; smaller-on-larger -> welded, volume = sum, 0 pierces.
v-2c. Seed-target hardening + FarFromOrigin (EXPECT_LT restored).
v-2d. Hull fixture integration (EXPECT_EQ 0, tie assert restored,
      TEMP DEBUG instrumentation removed).
