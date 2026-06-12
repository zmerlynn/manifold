// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

// Internal data structures + helper-function declarations for the
// overlap-removal pipeline. Used by src/overlap_removal.cpp itself
// and by no production caller - kept out of src/overlap_removal.h so
// the src-internal seam is just RemoveOverlaps.

#include <set>
#include <utility>  // for std::pair
#include <vector>

#include "impl.h"               // Manifold::Impl (held by value below)
#include "manifold/common.h"    // vec3 alias
#include "manifold/manifold.h"  // for Manifold

namespace manifold {
namespace overlap_removal {

// Iteration-cap default for MergeVertsEps (a parameter default, so it
// lives here rather than in overlap_removal.cpp's anonymous
// namespace). Each pass moves merged verts to centroid and repeats
// until positions stabilize; converges in a few passes for working
// fixtures. The cap is a tripwire against pathological inputs.
constexpr int kMergeVertsMaxIter = 8;

// Canonical edge of an input mesh: an unordered pair of vert indices
// (v0 < v1 by convention) with the two halfedge ids that span it
// (forward = startVert == v0; paired is its twin or -1 on a non-
// manifold edge). Built once by EnumerateEdges() at pipeline start.
struct Edge {
  int v0, v1;           // vert indices, v0 < v1
  int halfedgeForward;  // halfedge id with startVert == v0
  int halfedgePaired;   // its pair (= -1 only on non-manifold edges)
};

// Verts that lie within eps of an Edge (strictly between its endpoints).
// Built by BuildOnEdgeVertLists in step 4. Sorted by parametric `t`
// along the edge, t in (0, 1). Used by PartitionFace to subdivide
// the edge's halfedges into sub-edges.
struct EdgeVertList {
  std::vector<int> verts;  // sorted by t along the edge
  std::vector<double> ts;  // parametric position 0 < t < 1
};

// Verts that lie within eps of a triangle's interior (strictly inside
// the barycentric simplex). Built by BuildOnTriVertLists in step 5.
// Used by step 6 (FindEdgeTriIntersections) to skip events at
// already-snapped verts, and by step 8 (AddInteriorVertsToNewEdges)
// to propagate verts onto chord edges.
struct TriVertList {
  std::vector<int> verts;  // verts within eps of interior
  std::vector<vec3> bary;  // barycentric (b0, b1, b2) - all > 0
};

// One edge-pierces-triangle event. Step 6 (FindEdgeTriIntersections)
// emits these via BVH broad phase + Moller-Trumbore narrow phase.
// `s` is the parameter along the edge (0 < s < 1, strict interior),
// `bary` is the barycentric on the triangle (all > 0). `snapTo` is
// the vert id this event resolves to if within eps of an existing
// vert; -1 if a fresh new vert needs allocating.
struct EdgeTriIntersection {
  int edgeIdx;    // index into edges[]
  int triIdx;     // triangle index in halfedge_/faceNormal_
  vec3 position;  // intersection point in 3D
  double s;       // parameter along the edge (0 < s < 1)
  vec3 bary;      // barycentric on the triangle (all > 0)
  int snapTo;     // -1 = new vert, >=0 = snap to existing vert
};

// A "chord" edge between two intersecting triangles. Step 7
// (GenerateChordEdges) generates one PiercedNewEdge per tri-tri
// pair that has exactly 2 distinct endpoints from EdgeTriIntersection
// events. The chord lies on both triA and triB and splits each into
// two sub-polygons.
struct PiercedNewEdge {
  int v0, v1;      // vert ids (sorted, v0 < v1; existing if < NumVert,
                   // else new with index `id - NumVert`)
  int triA, triB;  // the two tris this edge lies on (sorted)
};

// Output of GenerateChordEdges: the new vertex positions, chord
// edges, per-event resolved vert ids, and counters used by
// downstream stages (partition, classifier, emit).
struct ChordEdges {
  std::vector<vec3> newVertPositions;
  std::vector<PiercedNewEdge> newEdges;
  // Parallel to the EdgeTriIntersection input passed to
  // GenerateChordEdges. etIsect2Vert[i] is the vert id that
  // event i became - `snapTo` if >=0 in the input, otherwise a
  // freshly allocated id (>= NumVert) post-dedup. Used by
  // PropagateNewVertsToOnEdgeLists to add the new verts to the
  // on-edge lists of their piercing edges.
  std::vector<int> etIsect2Vert;
  int droppedTriTriPairsWithBadEndpointCount = 0;
};

// Step 8 (AddInteriorVertsToNewEdges) augments each PiercedNewEdge
// with the verts that lie on its interior (= snap-merged on-tri
// verts of either triA or triB that fall on the chord segment).
// Used by PartitionFace to subdivide the chord into sub-edges.
struct NewEdgeWithExtras {
  PiercedNewEdge edge;
  std::vector<int> extraVerts;  // sorted along the edge
  std::vector<double> extraTs;
};

// Setup helper: scale-invariant eps derived from an impl's bounding-
// box scale via AlphaBudgetEpsilon (in src/shared.h). Larger meshes get
// larger eps.
double InferEps(const Manifold::Impl& m);

// Result of MergeVertsEps below. `maxMove` is the largest total
// displacement any input vert received (final cluster centroid vs its
// input position - a CHAIN of eps-pairs can move a member well beyond
// eps across passes). The driver folds it into the output tolerance
// claim. When `mergedCount == 0`, `impl` is left default-constructed
// and the caller proceeds on its own input unchanged - nothing was
// rebuilt, so nothing can have drifted.
struct MergeVertsResult {
  Manifold::Impl impl;
  int mergedCount = 0;
  double maxMove = 0.0;
};

// Step 1 of the overlap-removal pipeline: merges all verts within eps
// of each other. Iterates broad-phase Collider self-collisions +
// DisjointSets union, applies cluster-centroid positions, then
// rebuilds an Impl directly: tri verts remapped to cluster
// representatives, collapsed tris dropped, followed by the house
// construction sweep (CreateHalfedges, CleanupTopology,
// SetNormalsAndCoplanar, RemoveDegenerates, RemoveUnreferencedVerts,
// SortGeometry) - the same invariant chain the MeshGL ctor provides,
// owned here because no ctor runs. The input's tolerance_ carries
// into the rebuild (SetEpsilon floors, never lowers it).
//
// `mergedCount` is the authoritative answer to "did anything get
// merged?" - counted directly from the union-find, not inferred from
// vert counts (the rebuild's sweep compacts merged-away verts, but
// counting by construction beats reverse-engineering the sweep).
// `ctx` (nullable) threads into the rebuild's finalize sweep
// (SortGeometry); a mid-sweep cancel comes back as an Impl with
// Error::Cancelled, which the driver propagates as observable status.
MergeVertsResult MergeVertsEps(const Manifold::Impl& in, double eps,
                               ExecutionContext::Impl* ctx,
                               int maxIter = kMergeVertsMaxIter);

// Step 2 of the pipeline (the sketch's step 3 is subsumed by the
// canonical-edge representation; numbering follows the sketch):
// enumerate canonical edges of a Manifold's
// Impl. Each undirected edge is represented once, with the canonical
// "forward" halfedge (= the one whose startVert < endVert) and its
// pair (-1 only on non-manifold edges). Order is deterministic.
std::vector<Edge> EnumerateEdges(const Manifold::Impl& impl);

// Step 4 of the pipeline: for each edge, find verts that lie within
// eps of its interior segment. Uses BVH broad phase + closest-point
// narrow phase. Skips edge endpoints, and skips verts that are
// neighbors of BOTH endpoints (= thin-tri apex case where the
// "near" vert is just a structural neighbor, not a real overlap).
std::vector<EdgeVertList> BuildOnEdgeVertLists(const Manifold::Impl& impl,
                                               const std::vector<Edge>& edges,
                                               double eps);

// Step 5 of the pipeline: for each triangle, find verts that lie
// strictly inside its barycentric interior (within eps of the plane,
// all barycentrics > 0). Uses BVH broad phase + plane-distance +
// strict barycentric narrow phase. Skips the triangle's own 3 verts.
std::vector<TriVertList> BuildOnTriVertLists(const Manifold::Impl& impl,
                                             double eps);

// Step 6 of the pipeline: edge-pierces-triangle events via BVH +
// Moller-Trumbore narrow phase. Strict-interior gates on segment
// parameter (0 < s < 1) AND barycentric (all > 0). Snaps the pierce
// point to the nearest existing vert within eps - event identity at
// the computational scale (the step-1 old-old convention;
// arrangement-wide identification happens at step 9.5), ties to
// smallest id.
//
// PAIR inclusion: only the edge's own two faces are skipped. In
// particular, a pair whose TRIANGLES share a vert (the post-Boolean-
// merge adjacency a naive #289 step 6 would skip) still yields events
// from the non-shared edges - an edge with an ENDPOINT on the tri
// cannot event anyway (a segment crosses the plane once, at s == 0).
// EVENTS stay strict-interior on both parameters regardless: an edge
// ENDPOINT lying on the tri plane (s == 0/1) generates no event by
// design - a pair left with one event drops as a tip touch (see
// GenerateChordEdges below), and a vertex-on-face through-pierce is
// part of the documented arrangement-incompleteness class
// (docs/OverlapRemoval.md, Known limitations), caught by the
// fail-closed gates.
std::vector<EdgeTriIntersection> FindEdgeTriIntersections(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<TriVertList>& onTriLists, double eps);

// Step 7 of the pipeline: resolve etIsect events to vert
// ids (snapping or allocating fresh), group by tri-tri pair, and
// emit chord edges. Pairs with:
//   - 2 endpoints -> new chord edge between the two pierce points.
//   - 1 endpoint  -> edge-tip touch, no through-pierce: nothing is
//     emitted or recorded.
//   - 0 or >= 3   -> dropped (counted).
ChordEdges GenerateChordEdges(const Manifold::Impl& impl,
                              const std::vector<Edge>& edges,
                              const std::vector<EdgeTriIntersection>& etIsects,
                              double eps);

// Step 8 of the pipeline: per chord edge, find original-mesh verts
// that lie within eps of its interior segment. Candidates are the
// union of in-tri verts of triA and triB (= step 5 onTriList
// outputs). Emits NewEdgeWithExtras (chord + sorted interior verts).
std::vector<NewEdgeWithExtras> AddInteriorVertsToNewEdges(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<PiercedNewEdge>& newEdges,
    const std::vector<TriVertList>& onTriLists, double eps);

// ---- Step 6.5: coplanar trace chords (docs/OverlapRemoval.md) ----

// In-plane conformance for coplanar overlapping face pairs. Step 1's
// eps-merge flattens Boolean SoS slivers into zero-volume pancakes -
// coplanar overlapping faces with opposite orientations - by design:
// step 12's signed-multiplicity cancellation consumes the
// coincidences. But transversal step-6 events cannot cut coplanar
// pairs, so without this pass the sheets' partitions never conform,
// cancellation cannot fire, and BuildCellComplex hits exact angular
// ties. For each BVH pair surviving the plane gate (every vert within
// eps of the LARGER face's plane), both faces' edges are clipped
// against the other face in ONE shared frame (the gate face's - the
// larger, ties to the lower id - so a crossing is computed once and
// shares its id across both clip directions). An interval qualifies
// iff it is longer than eps AND its midpoint is interior to the other
// face by > eps (full-through cuts qualify - their midpoints are
// interior; boundary-riding intervals from coplanar neighbors never
// do, so clean flat meshes emit nothing). Crossing endpoints snap to
// the pair's six corners at max(2 * eps, the conditioned radius)
// (nearest, ties to smallest - the step-9 convention), then
// dedup new-to-new over the pool at the SOURCE-GATED radius
// min(this crossing's conditioned radius, the entry's recorded
// radius), eps-floored; an endpoint farther than eps from either
// original 3D edge rejects its interval (the near-grazing guard).
// Each qualifying interval emits a chord lying on BOTH faces;
// interior crossing endpoints - allocated AND snapped - return as
// explicit (edge, vert, t) on-edge additions for every claimed edge
// whose endpoint they are not (an X crossing gets two records with
// one vert id).
struct OnEdgeAddition {
  int edge;  // index into edges[]
  int vertId;
  double t;
};
struct TraceChordResult {
  std::vector<PiercedNewEdge> chords;
  std::vector<vec3> newVertPositions;  // input extended (taken by value)
  std::vector<OnEdgeAddition> onEdgeAdditions;  // deduped (edge, vertId)
  // Parallel to newVertPositions: each new vert's conditioned snap
  // radius (eps for pre-existing pool entries; the capped
  // eps / sin(angle) of the crossing for step-6.5 allocations,
  // widened when a dedup hit claims more). Step 9.5 uses it for the
  // new-onto-original snap, closing the gap where one
  // pair corner-snaps an ill-conditioned crossing while the partner
  // pair allocates - leaving twins (10 eps, condR] apart.
  std::vector<double> newVertSnapR;
  int intervalsRejected = 0;  // boundary-riding / grazing / sub-eps
  // Faces that are CANCELLATION-HAZARD: among plane-gated pairs (taken
  // BEFORE interval filtering, so equal-boundary coincidences whose
  // boundary-riding intervals are rejected still count), faces whose
  // winding normals - normalized cross(v1-v0, v2-v0), never the stored
  // faceNormal_ - are ANTI-ALIGNED and whose in-plane bounding boxes
  // overlap. Same-oriented coplanar neighbors (every flat face's own
  // triangulation) are NOT flagged: they cannot form a step-12
  // cancellation pair, and flagging them would re-gate the
  // boss-through-plate fixture. A flagged face with a detected island
  // gates as today; an unflagged island face uses the decomposition.
  std::set<int> coplanarHazardFaces;
};
TraceChordResult CoplanarTraceChords(const Manifold::Impl& impl,
                                     const std::vector<Edge>& edges,
                                     const std::vector<int>& halfedge2Edge,
                                     std::vector<vec3> newVertPositions,
                                     double eps);

// Apply explicit on-edge additions: id-dedup against the existing
// list, then per-edge re-sort by t. The trace-chord sibling of
// PropagateNewVertsToOnEdgeLists (whose interface is parallel to
// etIsects and structurally cannot carry these).
void AddVertsToOnEdgeLists(const std::vector<OnEdgeAddition>& additions,
                           std::vector<EdgeVertList>& onEdgeLists);

// ---- Step 9: chord-chord crossings within each triangle ----
// (docs/OverlapRemoval.md.)

// Step 9 grouping: chord indices incident to each face. A chord lies
// on both its triA and triB, so it appears under both faces.
std::vector<std::vector<int>> GroupChordsByFace(
    const std::vector<PiercedNewEdge>& newEdges, int numTri);

// Step 9 pass 0: endpoint-on-chord contacts - the on-chord analog of
// the step-4 on-edge vert lists. A chord endpoint lying on another
// same-face chord's interior, within max(2*eps, snapR(endpoint)) and
// outside the endpoint-proximity zone in t-space (t in
// (contactR/len, 1 - contactR/len)), is recorded for threading.
// These records are consulted by the step-9 canonical-id resolution
// before any new vert is allocated; without this pass,
// IntersectSegments' near-endpoint rejection would silently drop
// these contacts (near-line slivers otherwise).
// snapR is the newVertSnapR view (parallel to newVertPositions):
// original endpoints carry zero conditioned radius (their 2*eps base
// covers them); new endpoints carry their recorded conditioned radius.
struct OnChordContact {
  int chord;   // chord index gaining the vert
  int vertId;  // existing vert id inserted onto it
  double t;    // parameter along that chord
};
std::vector<OnChordContact> FindOnChordEndpointContacts(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<double>& newVertSnapR,
    const std::vector<std::vector<int>>& face2Chords, double eps);

// Step 9 passes 2-3: one raw chord-chord crossing found in a face's
// plane (boolean2::IntersectSegments on the chords re-projected onto
// that plane), before the nearby-crossing merge and canonical-id
// resolution. A pair sharing two faces is recorded once (lowest face
// wins); chordA < chordB.
struct ChordChordCrossing {
  vec3 pos;            // crossing position, lifted back to 3D
  int chordA, chordB;  // indices into the chord vector
  double tA, tB;       // parameter along each chord
  int face;            // face whose plane hosted the kernel call
  // Conditioned radius of this crossing: eps / sin(angle between the
  // chord lines), eps-floored, capped at kCondSnapCapEps * eps - the
  // same conditioning as step 6.5's trace crossings. Carried through
  // clusters and allocation so step 9.5's new-onto-original snap can
  // widen for ill-conditioned step-9 verts too.
  double snapR;
};
std::vector<ChordChordCrossing> FindChordChordCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& face2Chords,
    VecView<const vec3> faceNormals, double eps);

// Step 9 singleton-cluster form (each raw crossing its own cluster;
// the test seam for resolution/threading without the merge -
// production goes through MergeAndPropagateCrossings +
// ResolveAndThreadClusters). Canonical-id resolution is
// resolve-then-allocate: snap to any existing endpoint / on-chord
// vert / pass-0 contact within the candidate-class radius BEFORE
// allocating, symmetric across both chords - a crossing must never
// thread as an endpoint id on one chord and a fresh id on the other.
// Candidate class split (the #1757 source-gate convention):
//   ORIGINAL candidates: max(2*eps, cl.snapR)
//   NEW candidates:      max(2*eps, min(cl.snapR, snapR(candidate)))
// Threading recomputes every t from the resolved position, re-applies
// the pass-0 endpoint-zone guard, id-dedups over the unified pass-0
// + crossing list, then t-sorts with an eps/len dedup backstop.
struct ChordCrossing {
  vec3 pos;                 // crossing position (face plane)
  int id;                   // canonical vert id (existing or fresh)
  std::vector<int> chords;  // incident chords
  std::vector<double> ts;   // parallel to chords (pre-recompute)
  double snapR = 0.0;       // max member conditioned radius (see above)
};
struct Step9Threading {
  std::vector<NewEdgeWithExtras> chords;
  std::vector<vec3> newVertPositions;
  // Parallel to newVertPositions: per-vert conditioned snap radius
  // (the incoming TraceChordResult::newVertSnapR entries, extended
  // with each fresh crossing vert's cluster radius, widened when a
  // snap onto an existing new vert claims more). Step 9.5 consumes it.
  std::vector<double> newVertSnapR;
  std::vector<ChordCrossing> crossings;
};
Step9Threading ResolveAndThreadCrossings(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions, std::vector<double> newVertSnapR,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<OnChordContact>& contacts, double eps);

// Step 9 nearby-crossing merge and eager propagation (runs between
// FindChordChordCrossings and ResolveAndThreadClusters). Raw crossings unite
// (union-find, sorted pair order) when they share an incident face AND lie
// within eps - the face gate, not a chord gate, so a 4-chord concurrence
// whose two crossings share no chord still merges. Two distinct crossings in
// the (eps, 10*eps] band are NOT merged here but remain distinct; step 9.5's
// frame-spread unification still unites them onto a REAL member position
// (not a manufactured centroid). Cluster position is the member centroid
// (ascending member order), re-projected onto the hosting face plane.
// Propagation then tests the cluster position against every chord incident
// to any involved face (point-to-segment <= eps, the pass-0 endpoint-zone
// t-guard re-applied) so a k-fold point lands on all k chords even when a
// pairwise crossing was missed. Output clusters carry id == -1;
// ResolveAndThreadClusters assigns canonical ids.
std::vector<ChordCrossing> MergeAndPropagateCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<std::vector<int>>& face2Chords,
    VecView<const vec3> faceNormals, double eps);

// Cluster form of ResolveAndThreadCrossings: resolution + threading
// over merged clusters (k incident chords). ResolveAndThreadCrossings
// delegates here with singleton clusters.
Step9Threading ResolveAndThreadClusters(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions, std::vector<double> newVertSnapR,
    const std::vector<ChordCrossing>& clusters,
    const std::vector<OnChordContact>& contacts, double eps);

// ---- Step 9.5: arrangement-wide new-vert unification ----

// Backstop dedup ACROSS allocation paths. Steps 6.5, 7, and 9 each
// dedup their own allocations, but the same geometric point computed
// through two different in-plane frames lands up to ~10 * eps apart
// (the crossing's lever arm on long near-parallel edges), and a pair
// of such twins subdivides a shared sub-edge inconsistently across
// faces - the unpaired sub-edges then read as open rims, whose
// ambient unification collapses the cell complex (observed: a handful
// of rims merged nearly every cell). One union-find sweep at the
// 9.5 unification radius (10 * eps - the "same point computed twice"
// constant; step 9's merge is eps, not this) unites new-new pairs and
// snaps new verts onto originals within max(that radius, the vert's
// recorded conditioned radius) - NEAREST wins across the cluster,
// ties to the smallest id, originals before new. Consumers are
// remapped in place:
// chord endpoints, chord extras (id-dedup, endpoint drops), and
// on-edge lists (id-dedup, endpoint drops, t re-sort). Returns the
// number of ids remapped plus the largest position displacement any
// remap applied (|pos(old) - pos(representative)|) - the driver folds
// that into the output tolerance claim.
// `perVertSnapR` (parallel prefix of newVertPositions; see
// TraceChordResult::newVertSnapR) widens the new-onto-original snap
// for verts whose allocation was ill-conditioned. Deliberately NO
// default argument: the driver must pass the threaded radii, and a
// default would let that handoff sever silently -
// callers with no radii pass {} explicitly.
struct UnifyResult {
  int changed = 0;
  double maxMove = 0.0;
};
UnifyResult UnifyArrangementVerts(const Manifold::Impl& impl,
                                  const std::vector<vec3>& newVertPositions,
                                  const std::vector<Edge>& edges,
                                  std::vector<EdgeVertList>& onEdgeLists,
                                  std::vector<NewEdgeWithExtras>& chords,
                                  double eps,
                                  const std::vector<double>& perVertSnapR);

// ---- Steps 10-11: per-face partition (docs/OverlapRemoval.md) ----

// Halfedge-id -> index into edges[] (both directions of an edge map
// to the same index; -1 where no canonical edge covers a halfedge).
// Built ONCE per pipeline run and shared by every PartitionFace call
// - rebuilding an all-edges lookup inside the per-face partition made
// it quadratic over the mesh.
std::vector<int> BuildHalfedgeToEdgeIndex(const Manifold::Impl& impl,
                                          const std::vector<Edge>& edges);

// Partition one face of the conforming step-9 arrangement into simple
// sub-polygon cycles (CCW w.r.t. the face's HALFEDGE WINDING - never
// the stored faceNormal_, which folded sheets invert). The face's three
// original edges contribute one halfedge per sub-edge (subdivided by
// onEdgeLists); incident chords contribute BOTH directions per
// sub-edge, deduped per face by undirected vert pair (coincident
// chords otherwise create exact angular ties). The angle-sorted walk
// skips the immediate reverse halfedge UNLESS it is the sole candidate
// (the U-turn that traverses dangling-chord spurs instead of
// stalling), closes on vertex arrival, and splits each closed cycle at
// repeated vert ids, dropping sub-3-vert spur loops (counted).
// Zero-length chords (step-9 snapping can collapse v0 == v1) are
// skipped at entry (counted). A stall is a DEBUG_ASSERT - with the
// conforming arrangement the walk is total.
struct FacePartition {
  std::vector<std::vector<int>> polygons;  // simple cycles of vert ids
  int spursDropped = 0;
  int zeroLengthChordsSkipped = 0;
  // >= 3-vert cycles with exactly-zero projected area (flattened
  // spurs: coincident post-merge positions under distinct ids, or an
  // exactly-collinear out-and-back walk). Dropped - their Newell
  // normal is undefined downstream.
  int degenerateCyclesDropped = 0;
  // Chord sub-edges coinciding with the face's own (subdivided)
  // boundary - step 6.5 trace chords riding their host's boundary
  // while cutting the partner face. Skipped: doubling a directed
  // edge makes the walk's exact-tie handling order-sensitive.
  int boundaryRidingSubEdgesSkipped = 0;
  // Verts of chord-only components that contain a CYCLE but attach
  // to the face boundary at fewer than two distinct verts (the
  // "stamp" class: a free interior loop, or a loop pinched onto one
  // boundary vert). Such a hole's surrounding region is a (possibly
  // pinched) annulus - not representable as simple cycles - so the
  // partition would silently erase the cut; a positive count means
  // the polygons are NOT a partition of this face, and the driver
  // fails the run closed.
  int interiorIslandVerts = 0;
};
// The walk frame derives from the face's OWN halfedge winding (NOT
// the stored faceNormal_, which can oppose the winding on folded
// self-intersecting sheets - that mirror turns the face-on-left walk
// into a boundary-hugging walk).
// `eps`: pipeline epsilon for hole-decomposition triangulation and
//   sub-resolution hole fast-fail. No default: the severance-proofing
//   convention requires callers to pass it explicitly.
// `coplanarHazard`: this face has an anti-aligned coplanar partner whose
//   in-plane bbox overlaps - triangulating its islands independently
//   could break step-12 cancellation (unmatched diagonals). A flagged
//   face with a detected clean island gates exactly as today.
FacePartition PartitionFace(const Manifold::Impl& impl, int face,
                            const std::vector<Edge>& edges,
                            const std::vector<int>& halfedge2Edge,
                            const std::vector<EdgeVertList>& onEdgeLists,
                            const std::vector<NewEdgeWithExtras>& chords,
                            const std::vector<int>& faceChords,
                            const std::vector<vec3>& newVertPositions,
                            double eps, bool coplanarHazard);

// ---- Step 12: canonical polygon merge (docs/OverlapRemoval.md) ----

// Merge equivalent sub-polygons with signed multiplicity. The
// canonical key of a cycle is the lexicographically-smallest rotation
// among all rotations of the cycle AND of its reversal; the sign is
// +1 when the canonical form is a rotation of the cycle as walked
// (CCW w.r.t. its face's halfedge winding - the partition's frame),
// -1 when it is a rotation of the reversal. A simple cycle with distinct verts
// is never rotation-equivalent to its own reversal, so the sign is always
// well-defined. Entries summing to zero drop (coincident
// opposite-facing surfaces cancel). `face` is the first contributor's
// (the plane/normal source); output is ordered by canonical key.
struct MergedPolygon {
  std::vector<int> cycle;  // the canonical rotation
  int mult;
  // First contributor's source face. Unused by the positions-only v1
  // emit, but it is the designed hook for property interpolation (the
  // boolean_result pattern: barycentric on the source face for new
  // verts, pass-through for originals - resolved PER VERT via
  // EmitTopology::ring2Vert, since a post-cancellation polygon can
  // span several source faces). Do not remove as dead.
  int face;
};
std::vector<MergedPolygon> MergePolygons(
    const std::vector<std::pair<int, std::vector<int>>>& facePolygons);

// ---- Step 13.2-13.3: radial fans + cells (docs/OverlapRemoval.md) ----

// The radial fan of one arrangement edge: incident polygons sorted
// CCW about the a -> b axis by their in-face direction, with each
// polygon's facing recorded (frontCcw: the front (+normal) side faces
// the CCW-adjacent wedge). An exact angular tie is a step-12
// invariant failure (DEBUG_ASSERT; ascending-polygon-id fallback only
// as release determinism insurance).
struct EdgeFan {
  int a, b;                   // undirected edge verts, a < b
  std::vector<int> polygons;  // incident polygon indices, radial order
  std::vector<bool> frontCcw;
};
// The volume-cell structure of the merged arrangement: cells are the
// equivalence classes of polygon SIDES (2 * polygon + side, side 0 =
// front = +Newell-normal of the canonical cycle) united through the
// wedges between angularly-consecutive fan entries - the 3D twin
// structure that makes winding propagation global. Cell ids are
// renumbered by smallest member key. A polygon's boundary edge with a
// single incident polygon (an open sheet's rim) unites its own front
// and back, as the ambient space does.
struct CellComplex {
  std::vector<EdgeFan> fans;       // ordered by (a, b)
  std::vector<int> polySide2Cell;  // [2 * polygon + side] -> cell id
  int numCells = 0;
};
CellComplex BuildCellComplex(const Manifold::Impl& impl,
                             const std::vector<MergedPolygon>& polygons,
                             const std::vector<vec3>& newVertPositions);

// ---- Step 13.4-13.5: winding seed + propagation + keep ----

// Winding classification of the cell complex. Per connected component
// of the cell graph (cells joined where a polygon separates them), a
// segment cast from outside the arrangement's bbox to an interior
// point of one of the component's polygons - targets in DESCENDING
// area (ties ascending id); triangles use their centroid, longer
// cycles the largest ear of a REAL triangulation - seeds the
// arrival-side cell with the true ambient winding:
// signed crossings are counted against ALL other polygons, including
// other components' - which is why an extreme-vertex seed is wrong
// for nested components. BFS then propagates windings through the
// component (crossing a polygon front to back adds its signed
// multiplicity). A cast within eps of any degenerate contact
// (endpoint on a surface, edge/vert graze, near-parallel plane) is
// invalid in FULL - never skip one polygon and keep counting - and is
// retried on the component's next polygon (kSeedCastMaxTargets); a
// component exhausting its targets fails the classification
// (ok = false, the driver falls back to the input). A polygon whose
// two sides landed in one cell (an open sheet united through its rim)
// separates nothing: it propagates no winding and is never kept.
// Keep a polygon iff exactly one side is inside (winding > 0); flip
// marks kept polygons whose canonical cycle's normal points toward
// the inside cell (the emit reverses those so normals face outside).
struct CellWinding {
  std::vector<int> winding;  // [cell] -> winding number
  std::vector<bool> keep;    // [polygon]
  std::vector<bool> flip;    // [polygon] -> reverse cycle on emit
  int seedCasts = 0;         // casts attempted, retries included
  bool ok = false;
};
// `epsHint`: the pipeline's working epsilon (the driver's eps); the
// effective graze margin is max(epsHint, impl.epsilon_, machine eps
// at the arrangement scale) - without the hint the cast would run at
// impl.epsilon_, which can be tighter than the epsilon the
// arrangement itself was built with. Deliberately NO default: a
// caller omitting it would silently reintroduce that skew; unit
// fixtures with no pipeline eps pass a non-positive hint.
CellWinding ClassifyCells(const Manifold::Impl& impl,
                          const std::vector<MergedPolygon>& polygons,
                          const std::vector<vec3>& newVertPositions,
                          const CellComplex& cells, double epsHint);

// ---- Step 13 emit topology: inside-wedge twins + vertex rings ----

// Explicit output topology for the kept polygons. CreateHalfedges
// pairs halfedges by sort order, which mis-pairs when more than two
// kept polygons meet at an arrangement edge, so the twins are
// assigned here first: at each radial fan, restricted to kept
// polygons, the two flanking the same INSIDE (winding > 0) wedge are
// twins - bare fan adjacency would pair across an outside wedge and
// invert orientation (e.g. weld two solids that share the edge).
// Dropped polygons between two consecutive kept ones cannot change
// the wedge's inside-ness (a dropped polygon has equal inside-ness on
// both sides), so consecutive-kept wedges are well-defined.
// The twin assignment defines an abstract closed surface over the
// kept polygons (pre-triangulation); vertex rings are its orbits
// nextAroundVert(h) = nextInPolygonCycle(twin(h)), one output vert id
// per ring - two solids touching at a vert or edge get distinct
// output verts (subsumes SplitPinchedVerts), and by the ring-
// separation argument every output edge carries exactly 2 halfedges
// (release-checked into ok, DEBUG_ASSERTed). Rings are numbered by
// (geometric vert id, smallest incident kept polygon) -
// deterministic. Cycles are emitted in OUTWARD orientation
// (CellWinding::flip applied) over ring ids.
struct EmitTopology {
  std::vector<int> keptPolygons;            // ascending polygon ids
  std::vector<std::vector<int>> outCycles;  // [kept idx] -> ring-id cycle
  std::vector<int> ring2Vert;               // [ring id] -> arrangement vert
  bool ok = false;  // closed + exactly-2-halfedges checks passed
};
EmitTopology BuildEmitTopology(const std::vector<MergedPolygon>& polygons,
                               const CellComplex& cells,
                               const CellWinding& winding);

// Post-step-7 on-edge propagation (before step 8 and the partition):
// etIsect resolved verts onto their piercing edges' on-edge lists, so
// the partition subdivides those halfedges at the pierce points. t is
// recomputed from the RESOLVED vert's position (a snapped event's
// raw parameter can order against the geometry) and must land in
// (0, 1) - a snapped vert projecting outside a short edge's interior
// subdivides nothing. Skips verts that are already edge endpoints or
// already in the list.
void PropagateNewVertsToOnEdgeLists(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<int>& etIsect2Vert, const std::vector<Edge>& edges,
    std::vector<EdgeVertList>& onEdgeLists);

// Diagnostic for self-intersections in a Manifold. Counts piercing
// triangle pairs via BVH broad phase + SegmentPiercesTriInterior
// narrow phase. Pairs sharing 2+ verts (= adjacent across an edge)
// are skipped as legal-by-construction.
//
// Used both as the pipeline's internal pierce guard (= fall back to
// merged input if pipeline output has more pierces than input) and
// as the user-visible check that overlap-removal achieved its goal.
struct SelfIntersectionResult {
  int interiorPierces;        // strict interior pierce count
  double maxPierceMagnitude;  // perpendicular depth of deepest pierce
  int candidatesChecked;      // pairs that passed broad phase
  int adjacentPairsSkipped;   // pairs sharing 2+ verts (legal)
  int trianglesTotal;
};

SelfIntersectionResult CheckSelfIntersection(const Manifold::Impl& m,
                                             double relTol = 1e-12 /* =
                                                 kPipelineRelTol */);

// Final-gate helper (pre-emit): true if any EDGE-CONNECTED component
// of folded polygons (front cell == back cell, grouped by fold cell +
// shared undirected edge) encloses mult-weighted signed volume beyond
// the membrane bound area x kFoldedVolumePerAreaEps x eps. Folded
// membranes legitimately drop; a folded CLOSED shell means the
// arrangement failed to embed it, and emitting would silently delete
// its material - the driver falls back. Per component, not per cell
// (and edge-connected, not vert-connected): opposite-orientation
// shells sharing a cell - or merely touching at a snapped vert - must
// not net their signed volumes.
bool FoldedCellsEncloseVolume(const Manifold::Impl& impl,
                              const std::vector<MergedPolygon>& polys,
                              const std::vector<vec3>& newVertPositions,
                              const CellComplex& cells, double eps);

}  // namespace overlap_removal
}  // namespace manifold
