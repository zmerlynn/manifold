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
// the public-to-src interface is just RunOverlapRemoval.

#include <functional>  // for PolygonClassifierFn
#include <map>
#include <set>
#include <vector>

#include "collider.h"           // for Collider, Box
#include "manifold/common.h"    // vec3 alias
#include "manifold/manifold.h"  // for Manifold
#include "overlap_removal.h"    // for RunOverlapRemoval

namespace manifold {
namespace overlap_removal {

// Iteration-cap defaults for parametric pipeline functions. Each
// loop has an early-exit on convergence; these caps are tripwires
// against pathological inputs. Companion caps for non-parametric
// loops live in overlap_removal.cpp's anonymous namespace.
//
// kMergeVertsMaxIter: eps-merge passes. Each pass moves merged verts
//   to centroid, iterates if any pair moved or unioned. Converges
//   in 1-3 passes for working fixtures.
// kPierceReducerMaxIter, kDropExcessOuterMax, kTrimOrphansMaxRounds:
//   each iter strictly drops / changes something or breaks; bounded
//   by underlying element count. 32 is conservative headroom.
constexpr int kMergeVertsMaxIter = 8;
constexpr int kPierceReducerMaxIter = 32;
constexpr int kDropExcessOuterMax = 32;
constexpr int kTrimOrphansMaxRounds = 32;

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
// along the edge, t in (0, 1). Used by the polygon walker to subdivide
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

// A "chord" edge between two intersecting triangles. Step 7 phase 2
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
// downstream stages (polygon walker, classifier, triangulator).
struct ChordEdges {
  std::vector<vec3> newVertPositions;
  std::vector<PiercedNewEdge> newEdges;
  // Parallel to the EdgeTriIntersection input passed to
  // GenerateChordEdges. resolvedIds[i] is the vert id that
  // event i became - `snapTo` if >=0 in the input, otherwise a
  // freshly allocated id (>= NumVert) post-dedup. Used by
  // PropagateNewVertsToOnEdgeLists to add the new verts to the
  // on-edge lists of their piercing edges.
  std::vector<int> resolvedIds;
  int droppedTriTriPairsWithBadEndpointCount = 0;
};

// Step 8 (AddInteriorVertsToNewEdges) augments each PiercedNewEdge
// with the verts that lie on its interior (= snap-merged on-tri
// verts of either triA or triB that fall on the chord segment).
// Used by the polygon walker to subdivide the chord into sub-edges.
struct NewEdgeWithExtras {
  PiercedNewEdge edge;
  std::vector<int> extraVerts;  // sorted along the edge
  std::vector<double> extraTs;
};

// One half-edge in a per-triangle graph. Step 11 phase 1
// (BuildPerTriHalfedgeGraphs) emits these for the three input edges
// (subdivided by onEdgeLists) plus all new-chord edges (both
// directions, each subdivided by extraVerts).
struct PerTriHalfedge {
  int startVert, endVert;  // vert ids (>= NumVert means new vert)
  bool isFromNewEdge;      // false = original edge, true = new edge
};

// The per-triangle halfedge graph that the polygon walker traverses.
// `nextHalfedge` is filled by step 11 phase 2 (AddNextPointers) via
// 2D angle-sorted next-around-face rotation.
struct PerTriHalfedgeGraph {
  int triId;
  std::vector<PerTriHalfedge> halfedges;
  std::set<int> verts;
  // Filled by AddNextPointers (step 11 phase 2). Parallel to
  // halfedges[]: nextHalfedge[i] is the index of the halfedge
  // that follows halfedges[i] in a polygon walk on the left side.
  // -1 = unset (no pairing computed).
  std::vector<int> nextHalfedge;
};

// Output of the polygon walker (step 11 phase 3, WalkPolygons): the
// sub-polygons of one tri after splitting by chord edges. Degenerate
// 2-vert "sub-polygons" (= chord pairs with both endpoints interior
// to the parent triangle) are kept separate so the classifier and
// triangulator only see >= 3-vert polygons; pair-sym Phase 1 still
// consults degenerates to enforce chord-pair constraints when one
// corner is degenerate.
struct PolygonWalkResult {
  std::vector<std::vector<int>> polygons;  // each = sequence of vert ids
  std::vector<std::vector<int>> degeneratePolygons;
  int stalledHalfedges = 0;  // hit an unset next pointer
};

// Result of TriangulateAndEmit: the final Manifold output plus
// per-stage counters useful for tests + diagnostics.
struct TriangulationResult {
  Manifold output;
  int polygonsTriangulated;
  int trianglesEmitted;
  int trisDroppedTooSmall;
  int polygonsKept = 0;
  int polygonsDropped = 0;
  int polygonsAutoKept = 0;  // skipped classifier (1-poly tris)
};

// Morton-sorted BVH builder shared across pipeline stages. Each
// stage that does broad-phase BVH overlap (MergeVertsEps,
// BuildOnEdgeVertLists, BuildOnTriVertLists, FindEdgeTriIntersections,
// CheckSelfIntersection, ScanForPiercingPairs, DoCapPass's
// rebuildCapBVH) used to inline the same 15-line ritual:
//   1. compute bbox = Union of leaf boxes
//   2. compute Morton codes for each leaf
//   3. build perm + stable_sort by Morton code
//   4. permute boxes / morton codes into sorted order
//   5. construct Collider from sorted views
// Hoisted here so all callers share one convention. Returns a
// struct that owns the sorted storage; the contained Collider
// has already copied the leaf boxes internally, so it remains
// valid after move/copy of the wrapper.
struct SortedBVH {
  Collider collider;
  std::vector<Box> boxes;        // boxes in Morton-sorted order
  std::vector<uint32_t> morton;  // sorted Morton codes (parallel to boxes)
  std::vector<size_t> perm;      // perm[sortedIdx] = origIdx
};

// Build a Morton-sorted BVH from a list of leaf boxes. The bbox
// covering all leaves is computed internally (= Union of inputs).
// Empty input returns a SortedBVH with an empty Collider.
SortedBVH BuildSortedBVH(VecView<const Box> leafBoxes);

// Setup helper: scale-invariant eps derived from a manifold's bounding-
// box scale via AlphaBudgetEpsilon (in src/shared.h). Larger meshes get
// larger eps.
double InferEps(const Manifold& m);

// Result of MergeVertsEps below.
struct MergeVertsResult {
  Manifold manifold;
  int mergedCount = 0;
};

// Step 1 of the overlap-removal pipeline: merges all verts within eps
// of each other. Iterates broad-phase Collider self-collisions +
// DisjointSets union, applies cluster-centroid positions, emits the
// merge hints via MeshGL64 mergeFromVert/mergeToVert. Returns the
// merged manifold and the count of merged pairs.
//
// `mergedCount` is the authoritative answer to "did anything get
// merged?" - `Manifold::NumVert()` may not reflect the merge if the
// merged verts didn't cause any tri collapse (Manifold's
// RemoveUnreferencedVerts sets unreferenced positions to NaN
// without compacting vertPos_).
MergeVertsResult MergeVertsEps(const Manifold& in, double eps,
                               int maxIter = kMergeVertsMaxIter);

// Step 2 of the pipeline: enumerate canonical edges of a Manifold's
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
// point to an existing vert if within eps.
//
// Includes pierces where the edge endpoint coincides with a tri vert
// (= the shared-vert pierce case post-Boolean merge), which the classic
// Emmett #289 step 6 behavior would miss.
std::vector<EdgeTriIntersection> FindEdgeTriIntersections(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<TriVertList>& onTriLists, double eps);

// Step 7 phase 2 of the pipeline: resolve etIsect events to vert
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

// ---- Step 9: chord-chord crossings within each triangle ----
// (docs/Step9Design.md; implemented incrementally.)

// Step 9 grouping: chord indices incident to each face. A chord lies
// on both its triA and triB, so it appears under both faces.
std::vector<std::vector<int>> GroupChordsByFace(
    const std::vector<PiercedNewEdge>& newEdges, int numTri);

// Step 9 pass 0: endpoint-on-chord contacts - the on-chord analog of
// the step-3 on-edge vert lists. A chord endpoint lying on another
// same-face chord's interior, within tolerance + eps and outside the
// endpoint-proximity zone in t-space (t in (snap/len, 1 - snap/len)),
// is recorded for threading. These records are consulted by the
// step-9 canonical-id resolution before any new vert is allocated;
// without this pass, IntersectSegments' near-endpoint rejection would
// silently drop these contacts (near-line slivers otherwise).
struct OnChordContact {
  int chord;   // chord index gaining the vert
  int vertId;  // existing vert id inserted onto it
  double t;    // parameter along that chord
};
std::vector<OnChordContact> FindOnChordEndpointContacts(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& chordsByFace, double tolerance,
    double eps);

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
};
std::vector<ChordChordCrossing> FindChordChordCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& chordsByFace,
    VecView<const vec3> faceNormals, double eps);

// Step 9 passes 5-7 (increment (ii): single crossings, no merge yet).
// Canonical-id resolution is resolve-then-allocate: snap to any
// existing endpoint / on-chord vert / pass-0 contact within
// tolerance + eps BEFORE allocating, symmetric across both chords -
// a crossing must never thread as an endpoint id on one chord and a
// fresh id on the other. Threading recomputes every t from the
// resolved position, re-applies the pass-0 endpoint-zone guard,
// id-dedups over the unified pass-0 + crossing list, then t-sorts
// with an eps/len dedup backstop.
struct ChordCrossing {
  vec3 pos;                 // crossing position (face plane)
  int id;                   // canonical vert id (existing or fresh)
  std::vector<int> chords;  // incident chords
  std::vector<double> ts;   // parallel to chords (pre-recompute)
};
struct Step9Threading {
  std::vector<NewEdgeWithExtras> chords;
  std::vector<vec3> newVertPositions;
  std::vector<ChordCrossing> crossings;
};
Step9Threading ResolveAndThreadCrossings(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<OnChordContact>& contacts, double tolerance, double eps);

// Step 10 of the pipeline: propagate etIsect resolved verts onto
// their piercing edges' on-edge lists, so the polygon walker
// subdivides those halfedges at the new pierce points. Skips verts
// that are already edge endpoints or already in the list.
void PropagateNewVertsToOnEdgeLists(
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<int>& resolvedIds, const std::vector<Edge>& edges,
    std::vector<EdgeVertList>& onEdgeLists);

// Position lookup helper that handles both original-mesh verts (id
// < baseId, into impl.vertPos_) and step-7 chord verts (id >=
// baseId, into newVertPositions). Used in many pipeline functions.
vec3 GetPos3(int id, int baseId, const Manifold::Impl& impl,
             const std::vector<vec3>& newVertPositions);

// Step 11 phase 1 of the pipeline: build per-tri halfedge graphs.
// For each tri T, emits one halfedge per sub-edge (along T's CCW
// direction) for each of T's 3 input edges (subdivided by
// onEdgeLists), plus BOTH directions of every new chord edge that
// touches T (subdivided by NewEdgeWithExtras::extraVerts).
std::vector<PerTriHalfedgeGraph> BuildPerTriHalfedgeGraphs(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<NewEdgeWithExtras>& newEdgesWithExtras);

// Step 11 phase 2 of the pipeline: compute next-around-face pointers
// for each per-tri graph. 2D project (drop dominant normal axis),
// atan2-sort outgoing halfedges per vert, find the next-clockwise
// outgoing as the "face on the left" walk's next pointer.
//
// Per-graph and batch overloads.
void AddNextPointers(const Manifold::Impl& impl,
                     const std::vector<vec3>& newVertPositions,
                     PerTriHalfedgeGraph& g);
void AddNextPointers(const Manifold::Impl& impl,
                     const std::vector<vec3>& newVertPositions,
                     std::vector<PerTriHalfedgeGraph>& graphs);

// Step 11 phase 3 of the pipeline: walk the polygon cycles in each
// per-tri graph using nextHalfedge. Each cycle = one sub-polygon
// (with >= 3 verts) or a degenerate 2-vert cycle (= chord pair with
// both endpoints interior to the parent triangle).
//
// Per-graph and batch overloads.
PolygonWalkResult WalkPolygons(const PerTriHalfedgeGraph& g);
std::vector<PolygonWalkResult> WalkPolygons(
    const std::vector<PerTriHalfedgeGraph>& graphs);

// Chord-partner lookup: for each (sorted-vert-pair, owning triId)
// of a chord edge, the partner triId on the other side. Built once
// from step 7's newEdges and consulted during pair-sym Phase 1 +
// the polygon classifier.
struct ChordPartnerMap {
  // key: (sorted v0, v1, owning triId), value: partner triId
  std::map<std::tuple<int, int, int>, int> partnerOf;
};

ChordPartnerMap BuildChordPartnerMap(
    const std::vector<PiercedNewEdge>& newEdges);

// Conversion helper: Manifold -> Impl via the public GetMeshGL64() API,
// to access internal halfedge / face-normal data.
Manifold::Impl ImplFromManifold(const Manifold& m);

// Geometric pierce predicate: does segment a-b strictly pierce
// triangle interior (v0, v1, v2)? Returns the pierce magnitude
// (perpendicular distance from nearer endpoint to plane) if yes,
// 0 if no.
//
// "Strict" excludes:
//   - Endpoint exactly on the plane (within FP tolerance).
//   - Intersection at a triangle edge or vertex.
// Used by both the in-pipeline pierce-aware cap walker and the
// post-pipeline CheckSelfIntersection diagnostic.
//
// `relTol` is the FP-noise threshold below which an endpoint is
// considered "on the plane" (returns 0). It is NOT a tolerance for
// filtering "small" pierces; it prevents zero-by-zero in the t
// computation.
double SegmentPiercesTriInterior(vec3 a, vec3 b, vec3 v0, vec3 v1, vec3 v2,
                                 double relTol = 1e-12);

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

SelfIntersectionResult CheckSelfIntersection(const Manifold& m,
                                             double relTol = 1e-12);

// Per-polygon classifier helper: decides whether to KEEP a polygon
// based on chord-bounded outward-side test. For each chord on the
// polygon's perimeter, find its partner triangle (= the tri on the
// other side of the chord) via chordPartners; if ALL non-chord
// polygon verts are on the OUTWARD side of that partner, drop the
// polygon (= it's in the carved-out region). Otherwise keep.
//
// "Outward" = dot(testPt - vertB, faceNormal_B) > threshold (where
// threshold is scale-relative to ignore borderline cases).
//
// Used as a fallback inside the production classifier when the
// per-vert sma-based classification (= AnalyzeSelfMesh) doesn't
// give a clear answer for an all-chord-vert polygon.
bool AnalyticalKeep(int triId, const std::vector<int>& polyVerts,
                    const ChordPartnerMap& chordPartners,
                    const std::vector<vec3>& positions3D, int baseId,
                    const Manifold::Impl& impl);

// Surface-cap walker: closes k=1 cycles in the current `out` mesh
// by fan-triangulation (with conflict-aware fan-apex selection) +
// greedy ear-clip fallback. Pierce-aware: refuses fan/ear tris
// that would pierce existing geometry (Moller-Trumbore against a
// BVH built per call). Forbidden-triple-aware: refuses tris whose
// vert triplet is in `forbiddenTriples` (= dropped by the post-cap
// pierce reducer; prevents drop+re-cap loops).
//
// Mutates `out.triVerts` (appends fan/ear tris). Returns
// (closed cycles, total tris added).
//
// Known limitation: k=1 cycles are non-planar space polygons, filled by
// pure combinatorial fan/ear-clip selected on edge-incidence and
// pierce-vs-existing-geometry only - there is no triangle orientation,
// planarity, or area check, so a badly-shaped cycle can over-inflate the
// surface (the documented cray volume blowup) or leave slivers. The
// post-pipeline volume-drift gate catches the gross case and falls back;
// sub-threshold inflation is the residual. A best-fit-plane projection +
// projected-simplicity check is the principled fix (tracked under #289).
std::pair<int, int> DoCapPass(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples);

// Pierce-aware reducer (PRE-cap): drops classifier-output tris that
// pierce each other. Catches "overlapping kept polygons" pierces
// the pierce-aware cap doesn't see (cap only checks NEW cap tris
// vs existing geometry).
//
// Algorithm: build BVH, find piercing pairs via Moller-Trumbore,
// drop one tri per pair (heuristic: higher pierce-count wins =
// more "load-bearing" to remove; ties -> lower triId). Iterate up
// to maxIters or until no pierces. Mutates out.triVerts in place.
// Returns total tris dropped.
int PierceAwareReducer(MeshGL64& out, int maxIters = kPierceReducerMaxIter);

// Pierce-aware reducer (POST-cap, with re-cap loop): same algorithm
// as PierceAwareReducer but alternates pierce-drop and cap-close.
// After dropping piercing tris, calls DoCapPass to close the
// resulting k=1 cycles. Iterate until fixed point. When the same
// pair set persists across iters (= drop+re-cap stuck), starts
// adding dropped tris to forbiddenTriples to break the cycle.
//
// Mutates out.triVerts AND forbiddenTriples (output param). Returns
// total tris dropped.
int PostCapPierceReducer(MeshGL64& out,
                         std::set<std::array<int, 3>>& forbiddenTriples);

// Directional k>2 reducer: for each edge with > 2 incidences, count
// contributions per direction. To make manifold (= 1 halfedge per
// direction), drop excess halfedges per direction; drop the entire
// tri containing the chosen halfedge. Then re-cap any new k=1
// cycles created by the drop. Iterate up to outerMax rounds.
//
// Drop heuristic: score-based. For each candidate tri, compute
// "drop cost" = number of its OTHER 2 edges currently at k=2 (=
// each becomes k=1 if dropped). Keep the candidate with HIGHEST
// cost (= most "load-bearing"); drop the rest.
//
// Mutates out.triVerts (drops + re-cap append). Returns
// (totalDropped, totalCapped, totalCapTris).
struct EdgeReducerResult {
  int totalDropped = 0;
  int totalCapped = 0;
  int totalCapTris = 0;
};
EdgeReducerResult DropExcessHalfedgeContributors(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples,
    int outerMax = kDropExcessOuterMax);

// Trim-orphans pass: drops tris with >= 2 k=1 edges (round 0) or
// >= 1 k=1 edges (subsequent rounds). Iterates with re-cap until
// no k=1 left or convergence. Effective on dangling chains the
// cap walker can't close.
//
// Mutates out.triVerts. Returns (totalDropped, totalCapped, totalCapTris).
EdgeReducerResult TrimOrphans(
    MeshGL64& out, const std::set<std::array<int, 3>>& forbiddenTriples,
    int maxRounds = kTrimOrphansMaxRounds);

// Per-polygon classifier function type. Returns a keep flag plus
// optional winding numbers (used in some classifier variants).
// Called by TriangulateAndEmit for polygons of multi-
// poly tris (single-poly tris are auto-kept).
struct PolygonClassification {
  bool keep;
  int windingUp;
  int windingDown;
};
using PolygonClassifierFn = std::function<PolygonClassification(
    const std::vector<int>&, const vec3&, int triId)>;

// Triangulate the kept polygons from each tri's PolygonWalkResult
// and emit them into a Manifold via MeshGL64. Composes:
//   1. Phase 1: tentative classification (= invoke classifier, or
//      auto-keep for single-poly tris).
//   2. Phase 2: cascade-drop forward (= drop auto-kept tris whose
//      sub-edges are entirely on the k=1 boundary).
//   3. Emit loop: for each kept polygon, fan-triangulate from on-
//      edge collinear apex if applicable, else manifold::Triangulate.
//   4. PierceAwareReducer (drops post-classifier overlapping tris).
//   5. DoCapPass (close k=1 cycles).
//   6. DropExcessHalfedgeContributors (drop excess directional contributors).
//   7. TrimOrphans (default-on; drop tris orphaned by the reducers).
//   8. PostCapPierceReducer (with re-cap loop, forbidden tracking).
//   9. Construct Manifold from out MeshGL64.
//
// classifier: optional. If null, all polygons are kept.
TriangulationResult TriangulateAndEmit(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<PolygonWalkResult>& walks,
    PolygonClassifierFn classifier = nullptr);

}  // namespace overlap_removal
}  // namespace manifold
