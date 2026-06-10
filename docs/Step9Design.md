# Step 9 design v6: chord-chord crossings within each triangle

Working doc for the faithful #289 rebuild. Not for upstream; lives in docs/ so
review subagents can read it. Clean up before any PR.

v6 = v5 + round 5 (the post-implementation review): two code bugs fixed (the
host-plane point rule and the guard-scope rule below), two clarifications, one
documented known-hole. v5 incorporates four adversarial review rounds. Round 1: architecture (output
struct, kernel choice, k>=3 merge). Round 2: contract (kernel drops
near-endpoint crossings; id-dedup; tolerance+eps; re-projection; stableEdgeId).
Round 3: identity (the split-identity band -> the resolve-then-allocate rule;
pre-pass t-guard; 10x-eps 3D justification; centroid re-projection). Round 4:
ordering (t-recompute from the resolved position; the pass-0 accumulator made
explicit; the propagation t-guard). Round 4 confirmed the architecture sound;
design considered CLEAR - implementation proceeds.

## House terminology and style (from the boolean2 review logs, all of them)

- halfedge structure (not DCEL); near-line sliver (not T-junction);
  nearby-crossing merge (not duplicate-crossing); on-edge / on-tri / on-chord
  vert lists for verts lying on a primitive within eps.
- Determinism is hard: any `std::` trig deciding output geometry is a bug -
  use `math::cos/sin/acos` (vendored deterministic).
- Return result structs, never out-pointers. `DEBUG_ASSERT` not `assert` (and
  note `DEBUG_ASSERT` is compiled out unless MANIFOLD_ASSERT && MANIFOLD_DEBUG -
  not a production guard). `la::cross`/`la::dot` directly. `static_cast<int>`,
  no C-casts. Explicit types over compact `auto`. Named locals over range-for
  on temporaries. File signpost below includes, above namespace.
- Dimensionally-correct thresholds: area ~ L^2 compares against `length*eps`.
  `Rect::Scale()` is a position scale, not extent - pick deliberately.
- TDD: red first, failing for the stated reason; characterization tests for
  behavior-preserving changes.

## State

`overlap_removal.cpp` demolished to steps 1-8 + a no-op `RunOverlapRemovalImpl`
stub (2529 -> 839 lines; steps 9-13 and the recovery scaffold deleted; builds
green). Backup of the full pre-demolish pipeline: branch
`backup/rsi-recovery-scaffold` (f9412101).

Steps 1-8 produce:
- `ChordEdges { newVertPositions, newEdges, resolvedIds }` from
  `GenerateChordEdges`; `newEdges` is `vector<PiercedNewEdge{v0, v1, triA,
  triB}>` (one chord = the segment where faces triA and triB cross; endpoint
  ids existing (< baseId) or new (>= baseId, into newVertPositions)).
- `vector<NewEdgeWithExtras { edge, extraVerts, extraTs }>` from
  `AddInteriorVertsToNewEdges`: each chord plus the existing on-tri verts on
  it, sorted by t. Today extraVerts only holds ids < baseId.

Missing (Emmett #289 step 9): nothing finds where two *chords* cross inside a
triangle, so crossing chords are never split and any later per-triangle walk
cannot partition that face.

## Goal

Per triangle T: find all chord-chord crossings interior to T, merge concurrent
crossings to a single vertex, split every involved chord there - including
endpoint-on-chord contacts (a chord endpoint lying on another chord's
interior). Output a conforming per-triangle chord set (no two sub-edges cross
except at shared verts). Step 9 completes arrangement geometry only - no
classification, no manifold build.

## Output

```cpp
struct ChordCrossing {
  vec3 pos;          // merged 3D position, re-projected onto the face plane
  int id;            // CANONICAL id: an existing endpoint id when within
                     // tolerance+eps (resolve-then-allocate), else a fresh
                     // id in the baseId+ newVertPositions space
  // incident chords and per-chord t's, from the producing pairs plus eager
  // propagation:
  std::vector<int> chords;
  std::vector<double> ts;     // parallel to chords
};
```

Plus: `newVertPositions` grown by allocated verts; per-chord
`extraVerts`/`extraTs` threading DERIVED from the resolved records and the
pass-0 insertions, id-deduped over the unified list. Contract change:
extraVerts may now contain ids >= baseId. `GetPos3` handles both ranges; add
`DEBUG_ASSERT(v < baseId)` at `AddInteriorVertsToNewEdges`'s direct
`impl.vertPos_[v]` read (it must never see a new id).

## Algorithm

0. **On-chord endpoint pass.** The on-chord analog of step 3's on-edge vert
   lists: for each face T and each ordered chord pair (c, d) in
   `chordsByFace[T]`, test each endpoint e of c against segment d:
   - point-to-segment 3D distance <= `tolerance + eps`, AND
   - t-guard: t in ((tolerance+eps)/len(d), 1 - (tolerance+eps)/len(d)) -
     the endpoint-proximity zone is EXCLUDED in t-space (round-3 finding 1b:
     a bare t in (0,1) lets an endpoint-near-endpoint contact insert as a
     spurious interior vert on a short chord).
   If both hold, record (e on d at t) into the NAMED pass-0 accumulator
   `std::vector<OnChordContact{chord, vertId, t}>` (round 4: these records are
   NOT yet threaded when resolution runs; resolution must consult this
   accumulator explicitly, so it is a first-class structure, not an implicit
   "already threaded" set). Avoids near-line slivers; rationale:
   `IntersectSegments` REJECTS crossings within eps of either endpoint
   (`AwayFromEndpoints`), so without this pass those contacts are silently
   lost.
1. **Group**: `chordsByFace[t]` = chord indices with `triA==t || triB==t`.
2. **Pairwise interior crossings**: per face T, per unordered chord pair
   (ci, cj):
   - Re-project both chords' endpoints onto T's plane first
     (`p -= la::dot(p - v0, n) * n`, unit n) so out-of-plane drift from prior
     merges is zero by construction; then map to 2D via T's orthonormal
     in-plane basis (a true isometry, so eps2d == eps3d).
   - Call `boolean2::IntersectSegments(GraphSegment2D a, b, eps, &p2)` with
     `stableEdgeId` = the chord's global index in `newEdges` (deterministic:
     assigned by the `GenerateChordEdges` map walk; NOT BVH pair order). Lift
     accepted p2 back to 3D via the basis.
   - Near-endpoint contacts are handled by pass 0; the kernel's rejection of
     them here is correct, not a loss.
3. **Collect raw crossings** (no ids yet): `{pos, ci, cj, tA, tB}`.
4. **Nearby-crossing merge** (boolean2's nearby-intersection merge analog):
   union-find over raw crossings; unite when (structural) the two crossings
   share an incident face AND (geometric) within 10*eps. A crossing's
   incident faces = its hosting face plus both chords' face pairs (round-5
   clarification). Cluster position = centroid (summed in ascending member
   order), then RE-PROJECTED onto the HOST face plane (round 3: centroid
   drift must not leave the plane) - host = the lowest hosting face among
   members, and the plane POINT must come from a member hosted on that face
   (round-5 code bug: raw input need not arrive in face order, so the first
   member's face can differ from host; pairing host's normal with another
   face's point projects onto neither plane). Face-gate rationale: a
   4-chord concurrence yields crossings (c1,c2) and (c3,c4) sharing no chord
   but sharing the face; the 10*eps geometric gate still bounds what can
   merge. 10*eps note for 3D: inherited from boolean2's measured constant
   (covers shallow crossings to ~6 degrees, structural gate prevents
   over-merge); unlike 2D, two genuinely-distinct same-face crossings within
   10*eps WILL merge - within the error budget this is the correct collapse,
   and a pinned test documents it.
5. **Eager propagation (position only, before any id exists)**: for each
   merged cluster, test its position against EVERY chord incident to any face
   its members involved; record (chord, t) where point-to-segment 3D distance
   <= eps and t passes the SAME guard as pass 0:
   t in ((tolerance+eps)/len(chord), 1 - (tolerance+eps)/len(chord))
   (round 4: a bare (0,1) here re-admits the spurious near-endpoint
   insertions pass 0's guard exists to block). The cluster's incident-chord
   set is now complete (k-fold points land on all k chords even when a
   pairwise intersection was missed).
6. **Canonical id resolution - resolve-then-allocate (round 3, the unifying
   rule; replaces all per-chord endpoint snapping)**: per cluster, ONCE,
   symmetric across all incident chords:
   - Gather candidates within `tolerance + eps` of the cluster position:
     every endpoint id of every incident chord; the existing on-chord verts in
     their extraVerts (the step-8 on-tri verts); AND the pass-0 accumulator's
     vertIds on those chords (round 4: pass-0 records are not yet threaded -
     without consulting the accumulator, a pass-2 crossing at an
     endpoint-contact point re-creates the split identity).
   - If any: canonical id = the nearest (tie: smallest id). NO new vert is
     allocated. (If two distinct existing ids are both within tolerance+eps -
     possible in the (eps, tolerance+eps] band where steps 1-8 kept them
     distinct - pick nearest/smallest deterministically and record the pair;
     full unification is the later global pass's job. Counted, not silent.)
   - Else: allocate a fresh id into newVertPositions.
   This is per-crossing, BEFORE threading, so a crossing can never thread as
   an endpoint id on one chord and a fresh id on another (the round-3
   split-identity bug). It also unifies with pass-0: a pass-2 crossing at an
   endpoint-contact point resolves to the same existing id pass 0 inserted.
7. **Threading**: per chord, gather ALL (id, t) records - pass-0 insertions
   and resolved cluster ids together - then:
   - RECOMPUTE every record's t from the RESOLVED position `GetPos3(id)`
     (round 4: t's computed against the pre-resolution cluster centroid can
     disagree with the snapped vertex by up to (tolerance+eps)/len, enough to
     reorder the t-sort against geometric order). A record whose recomputed t
     falls outside pass-0's guarded range is DROPPED for that chord - the
     vertex sits at/beyond the endpoint, whose id is already present.
   - id-dedup over the UNIFIED list (keep the first occurrence; round-3
     finding 1c),
   - sort by recomputed t (equal t's break by ascending id; the t-dedup keeps
     the first - round-5 clarification),
   - t-dedup at eps/len(chord) as the backstop,
   - write into extraVerts/extraTs.
   GUARD SCOPE (round-5 code bug): the endpoint-zone guard applies to
   step-9-ADDED records (contacts + crossings) only. Pre-existing step-8
   extras were admitted under step 8's weaker (0,1) guard and must survive a
   rebuild - dropping them loses legitimate arrangement verts.

## eps contract

- Base eps: the single absolute pipeline eps (steps 1-8's `eps` param, from
  `InferEps`). NOT kPipelineRelTol, NOT the probe eps.
- 1x eps: kernel acceptance, on-chord propagation test.
- `tolerance + eps` (the mesh's actual `tolerance_`): all snaps onto EXISTING
  ids (pass 0, canonical-id resolution) - new-to-old, matching boolean2's
  `newToOldThresh`.
- 10x eps: nearby-crossing merge (new-to-new), matching
  `kIntersectionMergeEpsFactor`.
- eps/len(chord): per-chord t-dedup backstop (after id-dedup).

## Determinism constraints (pin these)

- Chord enumeration: `newEdges` order (stable: `GenerateChordEdges` walks a
  `std::map` keyed by sorted (triA,triB); `etIsects` from a non-parallel
  collider traversal).
- Pair iteration: sorted (ci < cj) within each face; faces ascending.
- Union-find: unite in sorted pair order; centroid summed ascending.
- Candidate resolution: nearest, tie smallest id.
- Propagation and threading: ascending chord index, then ascending t.
- stableEdgeId: chord index in newEdges.
- No `std::` trig anywhere in the decision path (use `math::` if trig is ever
  needed; current design needs none).

## Degenerate cases

- Crossing at/near an endpoint: pass 0 + canonical-id resolution give it the
  existing id everywhere (no loss, no spurious vertex, no split identity).
- Endpoint near endpoint ((eps, tolerance+eps] band): excluded from pass-0
  insertion by the t-guard; canonical resolution picks one id
  deterministically and counts the residual pair.
- k>=3 concurrent chords: nearby-crossing merge (face-gated) + propagation.
- Collinear/overlapping chords: kernel yields no single crossing; pass 0
  snaps endpoint contacts; the overlap interval itself is deferred to the
  step-12 multiplicity merge. Counted, not silently dropped.
- Near-tangent crossings: position uncertainty grows as eps/theta; the 10x
  merge absorbs moderate cases; genuinely tangent chords resolve via pass 0 /
  collinear handling. Accepted limitation, documented.
- Disjoint chords: nothing inserted.
- KNOWN HOLE (round 5, accepted): when tolerance > 9*eps, two clusters can
  sit within tolerance+eps of each other yet beyond the 10*eps merge radius;
  the later cluster cannot snap to the earlier one's fresh id (fresh ids are
  not candidates until threading), so two near-coincident verts can be
  allocated. The t-dedup backstop collapses them only when both land on one
  chord within eps/len. Revisit if the tolerance-inflated regime becomes a
  target.

## Increments (implement + test in this order)

i.   chordsByFace grouping + on-chord endpoint pass (bookkeeping + the
     point-to-segment test with the t-guard; literal-fixture tests).
ii.  pairwise interior crossing + canonical-id resolution + threading, single
     crossings only. Safe to land before (iii) ONLY with strictly
     non-concurrent fixtures; the concurrent case must not be constructible
     in (ii)'s tests.
iii. the nearby-crossing merge (face-gated union-find) + eager propagation.

Note (ii) now includes resolution: even a single crossing must
resolve-then-allocate (the split-identity band exists with k=1).

## Tests

- Synthetic-literal fixtures: hand-built PiercedNewEdge/NewEdgeWithExtras with
  a consistent baseId/newVertPositions namespace. CHARACTERIZATION tests
  (prove the crossing math + resolution + threading + dedup), not
  reachability - reachability from a real mesh stays open until steps 10-13
  exist. TDD: each behavioral case red-first, failing for its stated reason.
- Specific cases pinned: simple X crossing; endpoint-on-chord; the
  (eps, tolerance+eps] band crossing (canonical id == endpoint id on BOTH
  chords - the split-identity regression); endpoint-near-endpoint (pass-0
  t-guard exclusion); 3 concurrent chords (one shared id on all three);
  4-chord concurrence across two disjoint pairs (the face-gate case);
  pass-0-vs-pass-2 same-point (one id after unified dedup); same-id-twice-
  on-one-chord (id-dedup); two-distinct-crossings-within-10eps-same-face
  (documented collapse); collinear overlap (counted, endpoints snapped).
- Feature tests: keep the 7 no-op-compatible tests enabled; HullMask and
  FarFromOrigin strict assertions become EXPECT_LE / DISABLED_ with a
  tracking comment.
- Cleanup alongside: remove dangling back-half decls + dead constants from
  overlap_removal_internal.h (keep PiercedNewEdge/NewEdgeWithExtras/
  ChordEdges; keep PropagateNewVertsToOnEdgeLists as prior art for
  propagation).

## Resolved review findings (audit trail)

- R1: output struct replaces extraVerts-only; kernel = IntersectSegments (a
  new 3D coplanar intersector has its own unguarded near-parallel
  degeneracy); merge gate by FACE not chord (4-chord concurrence).
- R2 HIGH: on-chord endpoint pass (kernel rejects near-endpoint crossings;
  a post-hoc snap could never fire). R2: id-dedup before t-dedup;
  tolerance+eps for new-to-old; re-project endpoints onto T's plane;
  stableEdgeId pinned; determinism pinned.
- R3 HIGH (split identity): canonical id resolution is per-crossing,
  resolve-then-allocate, before threading - never per-chord.
- R3 HIGH (pass-0/pass-2 unification): same rule + unified id-dedup in
  threading.
- R3 MEDIUM: pass-0 t-guard excludes the endpoint-proximity zone in t-space.
- R3 MEDIUM: 10x-eps face-gate documented + pinned for 3D.
- R3 LOW: merged centroid re-projected onto the face plane.
- R4 HIGH (t consistency): threading recomputes every t from the resolved
  GetPos3(id); out-of-guard recomputed t's drop.
- R4 HIGH (propagation guard): pass 5 uses pass-0's t-guard, not bare (0,1).
- R4 MEDIUM (pass-0 visibility): the OnChordContact accumulator is explicit
  and is a required candidate source in resolution. (Increment (ii) inherits
  this; R4 confirmed cascading snaps bounded and (ii) otherwise consistent.)
- R1 minors: DEBUG_ASSERT in AddInteriorVertsToNewEdges; increment-(ii)
  fixture constraint (now also: (ii) includes resolution).
- R5 (post-implementation review, 3 lanes): spec-conformance lane found full
  conformance; correctness lane found the two code bugs above (host-plane
  point; guard scope), both fixed red-first and pinned
  (Step9MergeReprojectsOntoHostFacePlane,
  Step9ThreadingPreservesStepEightExtras); the tolerance>9eps fresh-id hole
  documented as accepted. Test-adequacy lane: 7/10 pinned cases covered;
  remaining additions TODO before any PR - collinear-overlap (HIGH), mixed id
  space with a real Impl / baseId>0 (MEDIUM-HIGH), explicit id-dedup,
  t-reorder-after-snap, end-to-end five-function composition, 9x-eps
  documented-collapse.
