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

#include "overlap_removal.h"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

#include "collider.h"
#include "cross_section/boolean2/predicates.h"  // IntersectSegments (step 9)
#include "disjoint_sets.h"
#include "impl.h"
#include "manifold/polygon.h"  // for Triangulate
#include "overlap_removal_internal.h"
#include "shared.h"  // for AlphaBudgetEpsilon

// Internal pipeline implementation for Manifold::RemoveSelfIntersections().
// Implements Emmett Lalish's #289 13-step sketch: resolve the surface
// arrangement (chords, crossings, per-face partition, canonical merge),
// build the volume cell complex, classify cells by winding (seed cast +
// BFS), and keep the boundary between winding <= 0 and winding > 0.
// Fail-closed: any internal failure or gate trip returns the input
// unchanged. Design + known limitations: docs/OverlapRemoval.md.

namespace manifold {
namespace overlap_removal {

namespace {
// Numeric defaults shared across the pipeline. Centralized here so
// tuning + audits can change one value, not grep-and-replace.
//
// kPipelineRelTol: relative tolerance for SegmentPiercesTriInterior's
//   boundary-graze gate. 1e-12 lets eps-perturbed verts on tri edges
//   (within FP noise) NOT count as pierces, while keeping real
//   intersections detectable.
// kBarycentricFloor: minimum allowed barycentric coordinate when
//   classifying a point as strictly inside a triangle. 1e-12 rejects
//   points "on the boundary" (= one barycentric coord ~= 0) without
//   bumping into FP noise. Used by BuildOnTriVertLists and
//   FindEdgeTriIntersections.
constexpr double kPipelineRelTol = 1e-12;
constexpr double kBarycentricFloor = 1e-12;

// kSeedCastDir: direction from the arrangement bbox center to the
//   seed cast's source point P0 (step 13.4), at 2x the bbox diagonal.
//   A fixed generic unit vector - no axis alignment, no rational
//   component ratios - so casts into typical axis-aligned inputs
//   avoid edge/vert grazes on the first try. Unit to ~4 digits;
//   only its genericity and ~1 magnitude matter.
// kSeedCastMaxTargets: seed-cast target retries per cell-graph
//   component. Targets are tried in descending-area order, so the
//   budget only matters when the biggest polygons all graze - 8 is
//   generous headroom over typical usage.
// kCondSnapCapEps: cap, in eps units, on the condition-aware
//   corner-snap radius for trace-chord crossings (step 6.5). A
//   crossing's position error scales as eps / sin(angle) between the
//   two edge lines; 128 covers the near-parallel corner class with
//   headroom while keeping the snap local.
// kFoldedVolumePerAreaEps: final-gate threshold, in eps units, on a
//   folded cell's enclosed volume PER unit folded area. A fold set
//   (polygons whose front and back cells united, e.g. across a k = 1
//   rim) is dropped by the keep rule, which is correct for membranes:
//   their enclosed volume is below area x thickness, and post-merge
//   thickness is at most the 10 eps unification radius. A folded
//   CLOSED SHELL instead encloses volume ~ area x solid depth - tens
//   of orders above eps - and dropping it silently deletes real
//   material (observed: a tangent-degenerate mask contact folding an
//   entire disjoint hull). 100 gives 10x headroom over the thickest
//   legitimate membrane while sitting ~1e7 below any real shell.
const vec3 kSeedCastDir(0.278773, 0.581753, 0.764101);
constexpr int kSeedCastMaxTargets = 8;
constexpr double kCondSnapCapEps = 128.0;
constexpr double kFoldedVolumePerAreaEps = 100.0;

// Deterministic orthonormal in-plane basis for a unit face normal - a
// true isometry, so a kernel's 2D eps equals the pipeline's 3D eps
// (an axis-drop projection would not be an isometry and would
// contract distances). cross(u, v) == n, so CCW in (u, v)
// is CCW about the normal. Shared by step 6.5's clip frame, step 9's
// projection, and the step 10-11 angular ordering.
struct InPlaneBasis {
  vec3 u, v;
};
InPlaneBasis FaceBasisFromNormal(const vec3& n) {
  using la::cross;
  using la::dot;
  const double ax = std::fabs(n.x);
  const double ay = std::fabs(n.y);
  const double az = std::fabs(n.z);
  const vec3 ref = (ax <= ay && ax <= az) ? vec3(1.0, 0.0, 0.0)
                   : (ay <= az)           ? vec3(0.0, 1.0, 0.0)
                                          : vec3(0.0, 0.0, 1.0);
  vec3 u = cross(n, ref);
  u = u / std::sqrt(dot(u, u));
  return {u, cross(n, u)};
}

// Projection of p onto the line through a along ab (abLen2 =
// dot(ab, ab) > 0): the parameter t and the squared distance from
// the line. The shared kernel of the on-edge collection, the
// chord-interior threading, and the step-9 contact/propagation
// passes - callers apply their own t-range guards and radii.
struct LineProj {
  double t;
  double distSq;
};
LineProj ProjectToLine(const vec3& p, const vec3& a, const vec3& ab,
                       double abLen2) {
  using la::dot;
  const double t = dot(p - a, ab) / abLen2;
  const vec3 d = p - (a + ab * t);
  return {t, dot(d, d)};
}

// Conditioned snap radius of a crossing of two 2D directions: the
// position error scales as eps / sin(angle) (the near-parallel lever
// arm), eps-floored and capped at kCondSnapCapEps * eps; degenerate
// directions take the cap. Shared by step 6.5's trace crossings and
// step 9's chord-chord crossings.
double ConditionedSnapRadius(const vec2& d1, const vec2& d2, double eps) {
  const double l1 = std::sqrt(d1.x * d1.x + d1.y * d1.y);
  const double l2 = std::sqrt(d2.x * d2.x + d2.y * d2.y);
  const double crossAbs = std::fabs(d1.x * d2.y - d1.y * d2.x);
  return (l1 > 0 && l2 > 0 && crossAbs > 0)
             ? std::max(eps, std::min(eps * l1 * l2 / crossAbs,
                                      kCondSnapCapEps * eps))
             : kCondSnapCapEps * eps;
}

// Sort a parallel (vert id, edge parameter t) pair by t, breaking
// exact-t ties by ascending id so the order is independent of
// insertion order. Steps 8-11 consume these lists as monotone-t
// sub-edge sequences; shared by the on-edge/extra-vert sort sites.
void SortVertsByT(std::vector<int>& verts, std::vector<double>& ts) {
  std::vector<size_t> perm(verts.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(), [&](size_t i, size_t j) {
    if (ts[i] != ts[j]) return ts[i] < ts[j];
    return verts[i] < verts[j];
  });
  std::vector<int> sortedV(verts.size());
  std::vector<double> sortedT(ts.size());
  for (size_t k = 0; k < perm.size(); ++k) {
    sortedV[k] = verts[perm[k]];
    sortedT[k] = ts[perm[k]];
  }
  verts = std::move(sortedV);
  ts = std::move(sortedT);
}

// Seed-cast segment-vs-ear-triangle kernel (step 13.4). Transversal
// crossings only: anything within eps (a length) of a degenerate
// contact - a segment endpoint on the triangle's plane near the
// triangle, a crossing within eps of an ear edge (polygon boundary
// or interior diagonal alike), a near-in-plane segment - reports
// kGraze, which invalidates the WHOLE cast. The retry-on-next-target
// loop replaces SoS here; a per-polygon skip would corrupt the seed
// by a silent +-mult, so there is none.
enum class CastHit { kMiss, kHit, kGraze };
struct CastResult {
  CastHit kind;
  int step;  // on kHit: the winding increment per unit multiplicity,
             // -sign(dot(p1 - p0, ear normal))
};
CastResult CastSegmentAtEar(const vec3& p0, const vec3& p1, const vec3& a,
                            const vec3& b, const vec3& c, double eps) {
  using la::cross;
  using la::dot;
  const vec3 nRaw = cross(b - a, c - a);
  const double nLen = std::sqrt(dot(nRaw, nRaw));
  if (nLen == 0) return {CastHit::kGraze, 0};  // degenerate ear
  const vec3 n = nRaw / nLen;
  const double s0 = dot(p0 - a, n);
  const double s1 = dot(p1 - a, n);
  // Signed distance of an in-plane point from the ear's boundary:
  // the min over edges of the inward edge-line distance, positive
  // strictly interior. n perp each edge, so cross(n, dir) is unit.
  auto edgeMargin = [&](const vec3& x) {
    double dMin = std::numeric_limits<double>::infinity();
    const vec3 tri[3] = {a, b, c};
    for (int e = 0; e < 3; ++e) {
      const vec3 from = tri[e];
      vec3 dir = tri[(e + 1) % 3] - from;
      dir = dir / std::sqrt(dot(dir, dir));
      dMin = std::min(dMin, dot(cross(n, dir), x - from));
    }
    return dMin;
  };
  const bool near0 = std::fabs(s0) <= eps;
  const bool near1 = std::fabs(s1) <= eps;
  if (near0 && near1) return {CastHit::kGraze, 0};  // nearly in-plane
  if (near0 || near1) {
    // One endpoint within eps of the plane: degenerate only if its
    // plane contact is at/near the triangle itself - a coplanar
    // polygon far from the contact is a clean miss.
    const vec3 xNear = near0 ? p0 - n * s0 : p1 - n * s1;
    return {edgeMargin(xNear) <= -eps ? CastHit::kMiss : CastHit::kGraze, 0};
  }
  if ((s0 > 0) == (s1 > 0)) return {CastHit::kMiss, 0};
  const vec3 x = p0 + (p1 - p0) * (s0 / (s0 - s1));
  const double margin = edgeMargin(x);
  if (std::fabs(margin) <= eps) return {CastHit::kGraze, 0};
  if (margin < 0) return {CastHit::kMiss, 0};
  return {CastHit::kHit, s1 > s0 ? -1 : 1};
}

// Morton-sorted BVH wrapper shared by every broad-phase stage
// (MergeVertsEps, the on-edge/on-tri builders,
// FindEdgeTriIntersections, CheckSelfIntersection, step 9.5) - one
// bbox/Morton/sort/permute convention for all callers. Owns the
// sorted storage; the contained Collider has already copied the leaf
// boxes internally, so it remains valid after move/copy of the
// wrapper.
struct SortedBVH {
  Collider collider;
  std::vector<Box> boxes;         // boxes in Morton-sorted order
  std::vector<uint32_t> morton;   // sorted Morton codes (parallel to boxes)
  std::vector<size_t> leaf2Orig;  // [sorted leaf idx] -> input idx

  // Sequential broad phase: recorder.record(queryIdx, leafIdx) for
  // each query box overlapping a leaf box. Sequential deliberately:
  // every consumer feeds order-sensitive downstream stages, so the
  // deterministic visit order is part of the contract.
  template <typename Recorder, typename F>
  void Collisions(Recorder& recorder, F queryBox, int nQueries) const {
    collider.Collisions<false>(recorder, queryBox, nQueries,
                               /*parallel=*/false);
  }
};

SortedBVH BuildSortedBVH(VecView<const Box> leafBoxes) {
  SortedBVH out;
  const size_t n = leafBoxes.size();
  if (n == 0) return out;
  Box bbox;
  for (const auto& b : leafBoxes) bbox = bbox.Union(b);
  // MortonCode divides by the bbox extent per axis; a zero-extent axis
  // (e.g. all leaf centers coplanar) would yield NaN and an undefined
  // uint cast. Pad such axes with a scale-aware delta - a constant pad
  // is a no-op past 2^53 ULP scale (1e20 + 1.0 == 1e20), so size the
  // pad to the coordinate's own magnitude. Any positive extent maps
  // all centers to bucket 0 on that axis consistently.
  for (int k : {0, 1, 2}) {
    if (!(bbox.max[k] - bbox.min[k] > 0.0)) {
      const double pad =
          std::max(1.0, std::fabs(bbox.min[k]) * 4.0 *
                            std::numeric_limits<double>::epsilon());
      bbox.max[k] = bbox.min[k] + pad;
      DEBUG_ASSERT(bbox.max[k] - bbox.min[k] > 0.0, logicErr,
                   "BuildSortedBVH: degenerate axis pad failed");
    }
  }
  std::vector<uint32_t> rawMorton(n);
  for (size_t i = 0; i < n; ++i)
    rawMorton[i] = Collider::MortonCode(leafBoxes[i].Center(), bbox);
  out.leaf2Orig.resize(n);
  std::iota(out.leaf2Orig.begin(), out.leaf2Orig.end(), size_t{0});
  std::stable_sort(
      out.leaf2Orig.begin(), out.leaf2Orig.end(),
      [&](size_t a, size_t b) { return rawMorton[a] < rawMorton[b]; });
  out.boxes.resize(n);
  out.morton.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.boxes[i] = leafBoxes[out.leaf2Orig[i]];
    out.morton[i] = rawMorton[out.leaf2Orig[i]];
  }
  out.collider =
      Collider(VecView<const Box>(out.boxes.data(), out.boxes.size()),
               VecView<const uint32_t>(out.morton.data(), out.morton.size()));
  return out;
}

// Position lookup helper that handles both original-mesh verts (id
// < baseId, into impl.vertPos_) and step-7 chord verts (id >=
// baseId, into newVertPositions). Used in many pipeline functions.
vec3 GetPos3(int id, int baseId, const Manifold::Impl& impl,
             const std::vector<vec3>& newVertPositions) {
  if (id < baseId) {
    DEBUG_ASSERT(id >= 0 && id < static_cast<int>(impl.vertPos_.size()),
                 logicErr, "GetPos3: original vert id out of range");
    return impl.vertPos_[id];
  }
  const int j = id - baseId;
  DEBUG_ASSERT(j >= 0 && j < static_cast<int>(newVertPositions.size()),
               logicErr, "GetPos3: chord vert id beyond newVertPositions");
  return newVertPositions[j];
}

// Geometric pierce predicate: does segment a-b strictly pierce
// triangle interior (v0, v1, v2)? Returns the pierce magnitude
// (perpendicular distance from nearer endpoint to plane) if yes,
// 0 if no.
//
// "Strict" excludes:
//   - Endpoint exactly on the plane (within FP tolerance).
//   - Intersection at a triangle edge or vertex.
// Used by the CheckSelfIntersection diagnostic (the driver's
// pierce-monotonicity gate and the test-side pierce counter).
//
// `relTol` is the FP-noise threshold below which an endpoint is
// considered "on the plane" (returns 0). It is NOT a tolerance for
// filtering "small" pierces; it prevents zero-by-zero in the t
// computation.
double SegmentPiercesTriInterior(vec3 a, vec3 b, vec3 v0, vec3 v1, vec3 v2,
                                 double relTol = kPipelineRelTol) {
  using la::cross;
  using la::dot;
  const vec3 e1 = v1 - v0;
  const vec3 e2 = v2 - v0;
  const vec3 n = cross(e1, e2);
  const double nMag = std::sqrt(dot(n, n));
  if (nMag == 0) return 0.0;
  const double scaleA = std::sqrt(dot(a - v0, a - v0));
  const double scaleB = std::sqrt(dot(b - v0, b - v0));
  const double dTol = relTol * nMag * std::max(scaleA, scaleB);
  const double dA = dot(a - v0, n);
  const double dB = dot(b - v0, n);
  if (std::fabs(dA) < dTol || std::fabs(dB) < dTol) return 0.0;
  if (dA * dB >= 0) return 0.0;
  const double t = dA / (dA - dB);
  const vec3 p = a + t * (b - a);
  const vec3 c0 = cross(v1 - v0, p - v0);
  const vec3 c1 = cross(v2 - v1, p - v1);
  const vec3 c2 = cross(v0 - v2, p - v2);
  const double d0 = dot(c0, n);
  const double d1 = dot(c1, n);
  const double d2 = dot(c2, n);
  const double bTol = relTol * nMag * nMag;
  if (!(d0 > bTol && d1 > bTol && d2 > bTol)) return 0.0;
  return std::min(std::fabs(dA), std::fabs(dB)) / nMag;
}

// Shared construction sweep for pipeline-built Impls (the
// MergeVertsEps rebuild and the final emit): the same invariant chain
// the MeshGL ctor runs - CreateHalfedges, manifoldness check,
// CleanupTopology (pinched verts), SetNormalsAndCoplanar,
// RemoveDegenerates, RemoveUnreferencedVerts, SortGeometry - stated
// call by call because no ctor runs here, plus the explicit DERIVED
// metadata posture: one fresh reserved meshID on every triRef,
// originalID = -1. A rebuilt mesh is NOT an original (it matches what
// Manifold(MeshGL64) produces, and OriginalID() stays -1); callers
// wanting a provenance root use AsOriginal(). `toleranceSeed` lands
// in tolerance_ BEFORE SetEpsilon, which floors it at the
// bbox-derived epsilon_ and never lowers it - the working eps never
// overwrites epsilon_. Invariant failures come back as an error
// status (MakeEmpty), exactly as the ctor reports them. The sweep is
// ctx-aware at its heaviest stage: SortGeometry(ctx) can return early
// on cancel with partial state, so it is followed by an IsCancelled
// conversion to MakeEmpty(Error::Cancelled) - the Hull/LevelSet
// pattern.
Manifold::Impl BuildImplFromTris(std::vector<vec3>&& positions,
                                 const std::vector<ivec3>& tris,
                                 double toleranceSeed,
                                 ExecutionContext::Impl* ctx) {
  Manifold::Impl out;
  out.vertPos_.resize_nofill(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) out.vertPos_[i] = positions[i];
  Vec<ivec3> triProp;
  triProp.reserve(tris.size());
  for (const ivec3& t : tris) triProp.push_back(t);
  out.CreateHalfedges(triProp);
  if (!out.IsManifold()) {
    out.MakeEmpty(Manifold::Error::NotManifold);
    return out;
  }
  const int meshID = static_cast<int>(Manifold::Impl::ReserveIDs(1));
  auto& triRef = out.meshRelation_.triRef;
  triRef.resize_nofill(out.NumTri());
  for (size_t tri = 0; tri < triRef.size(); ++tri) {
    triRef[tri] = {meshID, meshID, -1, static_cast<int>(tri)};
  }
  out.meshRelation_.meshIDtransform[meshID] = {meshID, la::identity, false,
                                               false};
  out.meshRelation_.originalID = -1;
  out.CalculateBBox();
  out.tolerance_ = toleranceSeed;
  out.SetEpsilon();
  out.CleanupTopology();
  out.SetNormalsAndCoplanar();
  out.RemoveDegenerates();
  out.RemoveUnreferencedVerts();
  out.SortGeometry(ctx);
  if (IsCancelled(ctx)) {
    out.MakeEmpty(Manifold::Error::Cancelled);
    return out;
  }
  if (!out.IsFinite()) out.MakeEmpty(Manifold::Error::NonFiniteVertex);
  return out;
}

// An Impl carrying cancellation as observable status - the
// ADVANCE_PHASE_OR_RETURN idiom's return value. Never a silent
// fallback: the member wraps it like any errored result.
Manifold::Impl CancelledImpl() {
  Manifold::Impl impl;
  impl.MakeEmpty(Manifold::Error::Cancelled);
  return impl;
}

// Pipeline body - file-local so the public RemoveOverlaps below
// stays a thin fail-closed wrapper (every stage it composes is
// declared in overlap_removal_internal.h, so the definition can
// live here ahead of the stage implementations). Returns nullopt on
// the early-exit and every fallback arm: the CALLER owns "return the
// input bit-identically", which is a wrapper-identity concern, not a
// pipeline concern.
std::optional<Manifold::Impl> RemoveOverlapsImpl(const Manifold::Impl& input,
                                                 double eps,
                                                 ExecutionContext::Impl* ctx) {
  using la::cross;
  using la::dot;
  // Entry poll FIRST - before even the empty-input exit: a
  // pre-cancelled run is observable as Cancelled for every input
  // class, and the input pierce count below is the most expensive
  // pre-step, which a pre-cancelled caller should not pay.
  if (IsCancelled(ctx)) return CancelledImpl();
  if (input.IsEmpty()) return std::nullopt;
  if (eps <= 0) eps = InferEps(input);
  // The input's pierce count, computed once for the gate's
  // monotonicity arm.
  const int inputPierces = CheckSelfIntersection(input).interiorPierces;
  if (IsCancelled(ctx)) return CancelledImpl();

  // Step 1: merge verts within eps. Zero merges leaves the input
  // untouched (nothing rebuilt, nothing can drift).
  const MergeVertsResult merged = MergeVertsEps(input, eps, ctx);
  const Manifold::Impl& impl = merged.mergedCount > 0 ? merged.impl : input;
  // A cancel inside the merge rebuild's finalize sweep is observable
  // status, NOT a fallback - distinguish it before the generic arm.
  if (impl.status_ == Manifold::Error::Cancelled) return CancelledImpl();
  if (impl.IsEmpty() || impl.status_ != Manifold::Error::NoError) {
    return std::nullopt;
  }
  if (IsCancelled(ctx)) return CancelledImpl();
  const double tolerance = std::max(impl.tolerance_, eps);
  const int baseId = static_cast<int>(impl.NumVert());
  const int numTri = static_cast<int>(impl.NumTri());

  // Steps 2, 4-5 (no step 3 - see EnumerateEdges' note): canonical
  // edges, on-edge vert lists, on-tri vert lists.
  const std::vector<Edge> edges = EnumerateEdges(impl);
  std::vector<EdgeVertList> onEdgeLists =
      BuildOnEdgeVertLists(impl, edges, eps);
  const std::vector<TriVertList> onTriLists = BuildOnTriVertLists(impl, eps);

  // Step 6: edge-pierces-tri events.
  const std::vector<EdgeTriIntersection> etIsects =
      FindEdgeTriIntersections(impl, edges, onEdgeLists, onTriLists, eps);

  // Step 7: chords per tri-tri pair.
  ChordEdges chordEdges = GenerateChordEdges(impl, edges, etIsects, eps);

  // Step 6.5: coplanar trace chords - in-plane conformance for the
  // eps-merge-flattened pancake class. Appended to the chord list
  // BEFORE the early-exit, the per-face grouping, and step 8, so the
  // whole arrangement machinery consumes them unchanged.
  const std::vector<int> halfedge2Edge = BuildHalfedgeToEdgeIndex(impl, edges);
  TraceChordResult trace = CoplanarTraceChords(
      impl, edges, halfedge2Edge, std::move(chordEdges.newVertPositions),
      tolerance, eps);
  chordEdges.newVertPositions = std::move(trace.newVertPositions);
  chordEdges.newEdges.insert(chordEdges.newEdges.end(), trace.chords.begin(),
                             trace.chords.end());

  // EARLY-EXIT when the COMBINED chord list is empty: covers the
  // clean-input case (bit-identical return at the caller), the
  // all-pairs-dropped case, and pancake-free coplanar contact.
  // PRECEDENCE: a cancel RACING in since the last poll loses to this
  // exit - deliberately. The operation is effectively complete ("no
  // work to do"), and returning Cancelled here would make a clean
  // input's result depend on cancel timing; the FromMeshGL contract
  // sets the same precedent (validation outcomes win over a racing
  // cancel). Pre-cancelled callers were already caught at entry.
  if (chordEdges.newEdges.empty()) return std::nullopt;
  if (IsCancelled(ctx)) return CancelledImpl();

  // Boundary conformance for the trace crossings, then pierce verts
  // onto their piercing edges' on-edge lists, so the partition
  // subdivides those halfedges at the new verts.
  AddVertsToOnEdgeLists(trace.onEdgeAdditions, onEdgeLists);
  PropagateNewVertsToOnEdgeLists(impl, chordEdges.newVertPositions, etIsects,
                                 chordEdges.etIsect2Vert, edges, onEdgeLists);

  // Step 8: on-tri verts onto chord interiors.
  std::vector<NewEdgeWithExtras> chords = AddInteriorVertsToNewEdges(
      impl, chordEdges.newVertPositions, chordEdges.newEdges, onTriLists, eps);

  // Step 9: chord-chord crossings within each face - contacts,
  // pairwise crossings, nearby-crossing merge + propagation,
  // resolve-then-allocate threading.
  const std::vector<std::vector<int>> face2Chords =
      GroupChordsByFace(chordEdges.newEdges, numTri);
  const std::vector<OnChordContact> contacts = FindOnChordEndpointContacts(
      impl, chords, chordEdges.newVertPositions, face2Chords, tolerance, eps);
  const std::vector<ChordChordCrossing> rawCrossings =
      FindChordChordCrossings(impl, chords, chordEdges.newVertPositions,
                              face2Chords, impl.faceNormal_, eps);
  const std::vector<ChordCrossing> clusters = MergeAndPropagateCrossings(
      impl, chords, chordEdges.newVertPositions, rawCrossings, face2Chords,
      impl.faceNormal_, tolerance, eps);
  Step9Threading threaded = ResolveAndThreadClusters(
      impl, std::move(chords), std::move(chordEdges.newVertPositions),
      std::move(trace.newVertSnapR), clusters, contacts, tolerance, eps);

  // Step 9.5: unify new verts across allocation paths (the same
  // geometric point computed through two frames lands up to ~10 * eps
  // apart; unpaired twin sub-edges would read as open rims and
  // collapse the cell complex). The conditioned radii now cover the
  // step-9 allocations too (threaded, not just the trace pass's).
  const UnifyResult unified =
      UnifyArrangementVerts(impl, threaded.newVertPositions, edges, onEdgeLists,
                            threaded.chords, eps, threaded.newVertSnapR);
  if (IsCancelled(ctx)) return CancelledImpl();

  // Steps 10-11: partition every face (chordless faces still pick up
  // on-edge subdivision, so the arrangement conforms across shared
  // edges).
  std::vector<std::pair<int, std::vector<int>>> facePolygons;
  for (int f = 0; f < numTri; ++f) {
    if ((f & 0xFF) == 0 && IsCancelled(ctx)) return CancelledImpl();
    const bool hazard = trace.coplanarHazardFaces.count(f) > 0;
    FacePartition part = PartitionFace(
        impl, f, edges, halfedge2Edge, onEdgeLists, threaded.chords,
        face2Chords[f], threaded.newVertPositions, eps, hazard);
    // GATE (fail closed): an interior chord island means this face's
    // cycles are NOT a partition (the stamp class - see
    // FacePartition::interiorIslandVerts); emitting would silently
    // erase the cut and misclassify the stamping shell as nested.
    if (part.interiorIslandVerts > 0) return std::nullopt;
    for (std::vector<int>& cyc : part.polygons) {
      facePolygons.push_back({f, std::move(cyc)});
    }
  }

  // Step 12: canonical merge with signed multiplicity.
  const std::vector<MergedPolygon> polys = MergePolygons(facePolygons);
  if (polys.empty()) return std::nullopt;

  // Step 13: cells, winding, keep, emit topology. Classification or
  // topology failures (seed retries exhausted, unpaired halfedges,
  // ...) fall back to the input.
  const CellComplex cellCx =
      BuildCellComplex(impl, polys, threaded.newVertPositions);
  if (IsCancelled(ctx)) return CancelledImpl();
  const CellWinding winding =
      ClassifyCells(impl, polys, threaded.newVertPositions, cellCx, eps);
  if (!winding.ok) return std::nullopt;
  if (IsCancelled(ctx)) return CancelledImpl();
  // GATE (fail closed): a folded cell that encloses real volume means
  // a tangent-degenerate contact folded a closed shell's two sides
  // into one cell (its polygons all read front == back and the keep
  // rule drops them) - emitting would silently delete that shell.
  // Membranes legitimately fold flat and pass the area-relative
  // threshold; see FoldedCellsEncloseVolume.
  if (FoldedCellsEncloseVolume(impl, polys, threaded.newVertPositions, cellCx,
                               eps)) {
    return std::nullopt;
  }
  const EmitTopology topo = BuildEmitTopology(polys, cellCx, winding);
  if (!topo.ok || topo.keptPolygons.empty()) return std::nullopt;

  // Emit: one output vert per ring; triangulate each kept cycle in
  // its outward frame (Triangulate's CCW triangles project back
  // winding-consistent because the basis is right-handed).
  std::vector<vec3> ringPos(topo.ring2Vert.size());
  for (size_t r = 0; r < topo.ring2Vert.size(); ++r) {
    ringPos[r] =
        GetPos3(topo.ring2Vert[r], baseId, impl, threaded.newVertPositions);
  }
  std::vector<ivec3> outTris;
  for (const std::vector<int>& cyc : topo.outCycles) {
    if (cyc.size() == 3) {
      outTris.push_back(ivec3(cyc[0], cyc[1], cyc[2]));
      continue;
    }
    // Relative-origin Newell sum - see BuildCellComplex.
    const vec3 cycOrigin = ringPos[cyc[0]];
    vec3 nsum(0.0, 0.0, 0.0);
    for (size_t i = 0; i < cyc.size(); ++i) {
      nsum = nsum + cross(ringPos[cyc[i]] - cycOrigin,
                          ringPos[cyc[(i + 1) % cyc.size()]] - cycOrigin);
    }
    const double len2 = dot(nsum, nsum);
    DEBUG_ASSERT(len2 > 0, logicErr, "RemoveOverlaps: degenerate kept cycle");
    if (len2 <= 0) return std::nullopt;
    const InPlaneBasis basis = FaceBasisFromNormal(nsum / std::sqrt(len2));
    SimplePolygon poly2;
    poly2.reserve(cyc.size());
    const vec3 origin = ringPos[cyc[0]];
    for (const int r : cyc) {
      const vec3 d = ringPos[r] - origin;
      poly2.push_back(vec2(dot(d, basis.u), dot(d, basis.v)));
    }
    // Triangulate with epsilon-doubling retries: a kept cycle can
    // carry micro-tails of original verts clustered above the step-1
    // merge radius but below triangulable resolution (their ring ids
    // are topologically pinned, so the cycle cannot be simplified).
    // Widening epsilon moves the tail into the triangulator's own
    // degenerate class; the resulting zero-area tris collapse at
    // Manifold construction. Release builds return the same
    // triangulation without the debug CCW check, so behavior matches.
    std::vector<ivec3> tris;
    // Start from the pipeline eps, not just the mesh epsilon: kept
    // cycles carry pipeline-scale jitter (10 * eps merges), so seeding
    // the retry ladder below it just burns doubling attempts.
    double triEps = std::max({impl.tolerance_, impl.epsilon_, eps});
    // The retry ladder exists for MANIFOLD_DEBUG, where Triangulate's
    // CCW check can throw on micro-tail cycles; widening epsilon moves
    // the tail into the triangulator's own degenerate class. Capped at
    // 64x: far beyond that the check passes CW triangles spanning REAL
    // geometry, and the gate cannot be relied on to catch a mix.
    // Release Triangulate does not throw (it returns the same
    // triangulation unchecked), so behavior matches.
#ifdef MANIFOLD_DEBUG
    bool triangulated = false;
    for (int attempt = 0; attempt < 7 && !triangulated; ++attempt) {
      try {
        tris = Triangulate({poly2}, triEps, true);
        triangulated = true;
      } catch (...) {
        triEps *= 2.0;
      }
    }
    if (!triangulated) return std::nullopt;  // give the gate its fallback
#else
    tris = Triangulate({poly2}, triEps, true);
#endif
    // An empty triangulation of a >= 4-vert cycle would leave the
    // output topologically open: explicit fallback, matching the
    // debug ladder's, rather than relying on the volume gate.
    if (tris.empty()) return std::nullopt;
    for (const ivec3& t : tris) {
      outTris.push_back(ivec3(cyc[t[0]], cyc[t[1]], cyc[t[2]]));
    }
  }

  // Output: positions only - non-position properties are not
  // preserved (documented; MergedPolygon::face is the designed v2
  // hook). Built as an Impl directly, no Manifold(MeshGL64) in the
  // loop: the invariant sweep and the metadata policy are deliberate
  // calls in BuildImplFromTris, not ctor side effects. The tolerance
  // claim propagates the pipeline's MEASURED applied movements: the
  // 10 * eps floor covers the nearby-crossing merge radius (step-9
  // crossing merges, step-9.5 new-new unification), and the measured
  // step-1 cluster and step-9.5 remap displacements widen it when a
  // chain or a conditioned snap moved a vert further. Ill-conditioned
  // shallow-incidence corners can carry residual error beyond this, up
  // to the conditioned band (eps / sin(incidence), capped at
  // kCondSnapCapEps * eps) - a documented limitation, not part of the
  // tolerance claim.
  Manifold::Impl out = BuildImplFromTris(
      std::move(ringPos), outTris,
      std::max({tolerance, 10.0 * eps, merged.maxMove, unified.maxMove}), ctx);

  // GATE (thin, final): construction status, positive volume for a
  // non-empty input (NaN fails the comparison too), and pierce-
  // monotonicity. No volume-ratio tripwire: the input volume is the
  // winding-WEIGHTED integral, so heavy overlap legitimately shrinks
  // the measured volume (a full overlap reads 1/3). Gate failures
  // are an expected fallback for adversarial inputs, not asserts.
  // A cancel inside the emit's finalize sweep is observable status,
  // NOT a fallback - distinguish it before the generic arm.
  if (out.status_ == Manifold::Error::Cancelled) return out;
  if (out.status_ != Manifold::Error::NoError) return std::nullopt;
  if (!(out.GetProperty(Manifold::Impl::Property::Volume) > 0)) {
    return std::nullopt;
  }
  if (CheckSelfIntersection(out).interiorPierces > inputPierces) {
    return std::nullopt;
  }
  return out;
}
}  // namespace

double InferEps(const Manifold::Impl& m) {
  return AlphaBudgetEpsilon(m.bBox_.Scale(), 1000);
}

MergeVertsResult MergeVertsEps(const Manifold::Impl& in, double eps,
                               ExecutionContext::Impl* ctx, int maxIter) {
  // A nonpositive iteration cap is a degenerate INPUT, not the logic
  // regression the convergence tripwire below guards: fail closed
  // here, in every build config.
  if (in.IsEmpty() || maxIter <= 0) return {};

  // Cluster on a copy of the positions; `in.vertPos_` stays intact as
  // the displacement baseline. NaN positions cannot occur here:
  // RemoveUnreferencedVerts only MARKS unreferenced verts NaN, but
  // SortGeometry (the last stage of every construction sweep,
  // including BuildImplFromTris below) compacts them away, and this
  // pass only ever receives post-sweep impls. A NaN vert would NOT be
  // inert - it would poison the broad-phase bbox union and trip the
  // degenerate-axis pad assert in BuildSortedBVH.
  const size_t n = in.NumVert();
  std::vector<vec3> verts(n);
  for (size_t i = 0; i < n; ++i) verts[i] = in.vertPos_[i];

  // Per-pass: build eps/2-padded boxes, run Collider self-collision,
  // narrow-phase dist^2 < eps^2, unite via DisjointSets. Move each
  // cluster to its centroid and repeat until stable.
  int iter = 0;
  bool converged = false;
  std::vector<int> componentLabel(n);
  for (; iter < maxIter; ++iter) {
    std::vector<Box> boxes(n);
    const double halfEps = 0.5 * eps;
    for (size_t i = 0; i < n; ++i) {
      const vec3 lo(verts[i].x - halfEps, verts[i].y - halfEps,
                    verts[i].z - halfEps);
      const vec3 hi(verts[i].x + halfEps, verts[i].y + halfEps,
                    verts[i].z + halfEps);
      boxes[i] = Box(lo, hi);
    }
    SortedBVH bvh =
        BuildSortedBVH(VecView<const Box>(boxes.data(), boxes.size()));
    DisjointSets uf(static_cast<uint32_t>(n));
    const double eps2 = eps * eps;
    auto checkPair = [&](size_t qi, size_t li) {
      if (qi >= li) return;
      const size_t va = bvh.leaf2Orig[qi];
      const size_t vb = bvh.leaf2Orig[li];
      const vec3 d = verts[va] - verts[vb];
      const double d2 = la::dot(d, d);
      if (d2 > eps2) return;
      uf.unite(static_cast<uint32_t>(va), static_cast<uint32_t>(vb));
    };
    auto recorder = MakeSimpleRecorder(checkPair);
    auto qf = [&](int i) { return bvh.boxes[i]; };
    bvh.Collisions(recorder, qf, static_cast<int>(n));

    // Update positions to per-cluster centroid. A cluster whose
    // members are already bit-identical KEEPS that position: summing
    // n equal doubles and dividing is not bit-idempotent at large n
    // (iterated-addition rounding), and a drifting "centroid" of an
    // already-collapsed cluster could oscillate to the iteration
    // cap. With this skip, every cluster is
    // bit-identical one pass after it last grows, so convergence is
    // structural.
    int nComp = uf.connectedComponents(componentLabel);
    std::vector<vec3> sumByComp(nComp, vec3(0, 0, 0));
    std::vector<int> countByComp(nComp, 0);
    std::vector<vec3> firstByComp(nComp, vec3(0, 0, 0));
    std::vector<bool> allEqualByComp(nComp, true);
    for (size_t i = 0; i < n; ++i) {
      const int c = componentLabel[i];
      if (countByComp[c] == 0) {
        firstByComp[c] = verts[i];
      } else if (verts[i].x != firstByComp[c].x ||
                 verts[i].y != firstByComp[c].y ||
                 verts[i].z != firstByComp[c].z) {
        allEqualByComp[c] = false;
      }
      sumByComp[c] += verts[i];
      ++countByComp[c];
    }
    bool moved = false;
    for (size_t i = 0; i < n; ++i) {
      const int c = componentLabel[i];
      if (allEqualByComp[c]) continue;  // collapsed: keep, bit-idempotent
      const vec3 newPos = sumByComp[c] / static_cast<double>(countByComp[c]);
      if (newPos.x != verts[i].x || newPos.y != verts[i].y ||
          newPos.z != verts[i].z) {
        verts[i] = newPos;
        moved = true;
      }
    }
    // Convergence = a fixed point of POSITIONS: nothing moved this
    // pass, so this pass's component labels are final. A union count
    // is the wrong test - an already-merged coincident cluster
    // re-unites in every pass's fresh union-find, and whether that
    // reads as a "new" union depends on the union-find's internal
    // attachment order.
    if (!moved) {
      converged = true;
      break;
    }
  }
  DEBUG_ASSERT(converged, logicErr,
               "MergeVertsEps: hit kMergeVertsMaxIter without converging");
  // Convergence is structural (the bit-idempotence skip above), so a
  // non-converged exit means a logic regression or a nonpositive
  // maxIter: fail CLOSED in release rather than merging with the last
  // pass's labels - zero merges hands the driver its input unchanged.
  if (!converged) return {};

  // Total applied displacement (final centroid vs the INPUT position,
  // still intact in in.vertPos_) - exact, not a per-pass bound. The
  // driver folds it into the output tolerance claim. Singleton
  // clusters never moved (the centroid of one vert is itself).
  double maxMove = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const vec3 d = verts[i] - in.vertPos_[i];
    const double d2 = la::dot(d, d);
    if (d2 > 0.0) maxMove = std::max(maxMove, std::sqrt(d2));
  }

  std::vector<int> compRep(componentLabel.size(), -1);
  for (size_t i = 0; i < n; ++i) {
    const int comp = componentLabel[i];
    if (compRep[comp] == -1) compRep[comp] = static_cast<int>(i);
  }
  int mergedCount = 0;
  for (size_t i = 0; i < n; ++i) {
    if (compRep[componentLabel[i]] != static_cast<int>(i)) ++mergedCount;
  }
  // Zero merges: the caller proceeds on its own input - nothing was
  // rebuilt, so nothing can have drifted.
  if (mergedCount == 0) return {};

  // Rebuild directly: every tri's verts remapped to their cluster
  // representative (at the converged centroid positions held in
  // `verts`); tris that collapse to fewer than three distinct verts
  // drop here, the same dropping Manifold(MeshGL64) does with its
  // distinct-vert check before CreateHalfedges. BuildImplFromTris then
  // owns the construction sweep. The input's tolerance_ carries into
  // the rebuild (SetEpsilon floors, never lowers it).
  std::vector<ivec3> tris;
  tris.reserve(in.NumTri());
  for (size_t t = 0; t < in.NumTri(); ++t) {
    ivec3 tv;
    for (const int k : {0, 1, 2}) {
      tv[k] = compRep[componentLabel[in.halfedge_.Start(3 * t + k)]];
    }
    if (tv[0] != tv[1] && tv[1] != tv[2] && tv[2] != tv[0]) {
      tris.push_back(tv);
    }
  }
  return {BuildImplFromTris(std::move(verts), tris, in.tolerance_, ctx),
          mergedCount, maxMove};
}

std::vector<Edge> EnumerateEdges(const Manifold::Impl& impl) {
  // Each manifold edge is two halfedges; the canonical (forward) one
  // is the side where startVert < endVert. Walk all halfedges and
  // collect those.
  std::vector<Edge> out;
  out.reserve(impl.halfedge_.size() / 2);
  for (size_t h = 0; h < impl.halfedge_.size(); ++h) {
    const Halfedge he = impl.halfedge_.Get(h);
    if (!he.IsForward()) continue;
    out.push_back(
        {he.startVert, he.endVert, static_cast<int>(h), he.pairedHalfedge});
  }
  return out;
}

std::vector<EdgeVertList> BuildOnEdgeVertLists(const Manifold::Impl& impl,
                                               const std::vector<Edge>& edges,
                                               double eps) {
  using la::dot;
  const size_t nE = edges.size();
  const size_t nV = impl.vertPos_.size();
  const double eps2 = eps * eps;
  std::vector<EdgeVertList> out(nE);
  if (nE == 0 || nV == 0) return out;

  // vert->neighbor adjacency for the thin-tri-apex skip below (a
  // vert neighboring BOTH edge endpoints is a structural apex, not a
  // real on-edge overlap).
  std::vector<std::set<int>> adj(nV);
  for (size_t i = 0; i < impl.halfedge_.size(); ++i) {
    const int s = impl.halfedge_.Start(i);
    const int e = impl.halfedge_.End(i);
    adj[s].insert(e);
    adj[e].insert(s);
  }

  // Per-edge eps-padded AABB (BVH leaves) and per-vert eps-padded AABB
  // (queries).
  const vec3 pad(eps, eps, eps);
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.min -= pad;
    b.max += pad;
    edgeBoxes[i] = b;
  }
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }

  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(edgeBoxes.data(), edgeBoxes.size()));

  std::vector<std::vector<std::pair<double, int>>> hitsByEdge(nE);
  auto onCollision = [&](size_t vertIdx, size_t edgeIdxL) {
    const size_t edgeIdx = bvh.leaf2Orig[edgeIdxL];
    const auto& edge = edges[edgeIdx];
    if (static_cast<int>(vertIdx) == edge.v0 ||
        static_cast<int>(vertIdx) == edge.v1)
      return;
    if (adj[vertIdx].count(edge.v0) && adj[vertIdx].count(edge.v1)) return;
    const vec3 a = impl.vertPos_[edge.v0];
    const vec3 ab = impl.vertPos_[edge.v1] - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) return;
    const LineProj pr = ProjectToLine(impl.vertPos_[vertIdx], a, ab, abLen2);
    if (pr.t <= 0.0 || pr.t >= 1.0) return;
    if (pr.distSq <= eps2)
      hitsByEdge[edgeIdx].emplace_back(pr.t, static_cast<int>(vertIdx));
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  bvh.Collisions(recorder, qf, static_cast<int>(nV));

  for (size_t e = 0; e < nE; ++e) {
    auto& hits = hitsByEdge[e];
    std::sort(hits.begin(), hits.end());
    out[e].verts.reserve(hits.size());
    out[e].ts.reserve(hits.size());
    for (const auto& [t, v] : hits) {
      out[e].verts.push_back(v);
      out[e].ts.push_back(t);
    }
  }
  return out;
}

std::vector<TriVertList> BuildOnTriVertLists(const Manifold::Impl& impl,
                                             double eps) {
  using la::cross;
  using la::dot;
  const size_t nT = impl.NumTri();
  const size_t nV = impl.vertPos_.size();
  std::vector<TriVertList> out(nT);
  if (nT == 0 || nV == 0) return out;

  const vec3 pad(eps, eps, eps);
  std::vector<Box> triBoxes(nT);
  for (size_t t = 0; t < nT; ++t) {
    const int v0 = impl.halfedge_.Start(3 * t + 0);
    const int v1 = impl.halfedge_.Start(3 * t + 1);
    const int v2 = impl.halfedge_.Start(3 * t + 2);
    Box b(impl.vertPos_[v0], impl.vertPos_[v1]);
    b.Union(impl.vertPos_[v2]);
    b.min -= pad;
    b.max += pad;
    triBoxes[t] = b;
  }
  std::vector<Box> vertBoxes(nV);
  for (size_t i = 0; i < nV; ++i) {
    const vec3& p = impl.vertPos_[i];
    vertBoxes[i] = Box(vec3(p.x - eps, p.y - eps, p.z - eps),
                       vec3(p.x + eps, p.y + eps, p.z + eps));
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));

  auto onCollision = [&](size_t vertIdx, size_t triIdxL) {
    const size_t triIdx = bvh.leaf2Orig[triIdxL];
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    if (static_cast<int>(vertIdx) == t0 || static_cast<int>(vertIdx) == t1 ||
        static_cast<int>(vertIdx) == t2)
      return;
    const vec3 a = impl.vertPos_[t0];
    const vec3 b = impl.vertPos_[t1];
    const vec3 c = impl.vertPos_[t2];
    const vec3 e1 = b - a;
    const vec3 e2 = c - a;
    const vec3 n = cross(e1, e2);
    const double nMag2 = dot(n, n);
    if (nMag2 == 0) return;
    const double nMag = std::sqrt(nMag2);
    const vec3 p = impl.vertPos_[vertIdx];
    const double dist = dot(p - a, n) / nMag;
    if (std::fabs(dist) > eps) return;
    const vec3 pProj = p - dist * (n / nMag);
    const vec3 v0 = b - a, v1 = c - a, v2 = pProj - a;
    const double d00 = dot(v0, v0);
    const double d01 = dot(v0, v1);
    const double d11 = dot(v1, v1);
    const double d20 = dot(v2, v0);
    const double d21 = dot(v2, v1);
    const double denom = d00 * d11 - d01 * d01;
    if (denom == 0) return;
    const double bv = (d11 * d20 - d01 * d21) / denom;
    const double bw = (d00 * d21 - d01 * d20) / denom;
    const double bu = 1.0 - bv - bw;
    if (bu <= kBarycentricFloor || bv <= kBarycentricFloor ||
        bw <= kBarycentricFloor)
      return;
    out[triIdx].verts.push_back(static_cast<int>(vertIdx));
    out[triIdx].bary.push_back(vec3(bu, bv, bw));
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return vertBoxes[i]; };
  bvh.Collisions(recorder, qf, static_cast<int>(nV));
  return out;
}

std::vector<EdgeTriIntersection> FindEdgeTriIntersections(
    const Manifold::Impl& impl, const std::vector<Edge>& edges,
    const std::vector<EdgeVertList>& onEdgeLists,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using la::cross;
  using la::dot;
  std::vector<EdgeTriIntersection> out;
  const size_t nE = edges.size();
  const size_t nT = impl.NumTri();
  if (nE == 0 || nT == 0) return out;
  // Snap radius: eps - EVENT IDENTITY at the computational scale, the
  // same scale as the step-1 old-old merge (stored positions are the
  // geometry; only the pipeline's own error is absorbed). NOT the
  // tolerance + eps allocation radius: events in the (eps, 10 eps]
  // band are unified onto originals by step 9.5 anyway, and a
  // tolerance-scale radius on a tolerance-inflated input would drag
  // pierce events onto far verts and deform the arrangement (the
  // documented re-pierce class).
  const double snapR = eps;
  const double snapR2 = snapR * snapR;

  const vec3 pad(eps, eps, eps);
  std::vector<Box> edgeBoxes(nE);
  for (size_t i = 0; i < nE; ++i) {
    const vec3 v0 = impl.vertPos_[edges[i].v0];
    const vec3 v1 = impl.vertPos_[edges[i].v1];
    Box b(v0, v1);
    b.min -= pad;
    b.max += pad;
    edgeBoxes[i] = b;
  }
  std::vector<Box> triBoxes(nT);
  for (size_t t = 0; t < nT; ++t) {
    const int v0 = impl.halfedge_.Start(3 * t + 0);
    const int v1 = impl.halfedge_.Start(3 * t + 1);
    const int v2 = impl.halfedge_.Start(3 * t + 2);
    Box b(impl.vertPos_[v0], impl.vertPos_[v1]);
    b.Union(impl.vertPos_[v2]);
    b.min -= pad;
    b.max += pad;
    triBoxes[t] = b;
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(edgeBoxes.data(), edgeBoxes.size()));

  auto onCollision = [&](size_t triIdx, size_t edgeIdxL) {
    const size_t edgeIdx = bvh.leaf2Orig[edgeIdxL];
    const Edge& edge = edges[edgeIdx];
    if (static_cast<int>(triIdx) == edge.halfedgeForward / 3 ||
        (edge.halfedgePaired >= 0 &&
         static_cast<int>(triIdx) == edge.halfedgePaired / 3))
      return;
    const int t0 = impl.halfedge_.Start(3 * triIdx + 0);
    const int t1 = impl.halfedge_.Start(3 * triIdx + 1);
    const int t2 = impl.halfedge_.Start(3 * triIdx + 2);
    const vec3 a = impl.vertPos_[edge.v0];
    const vec3 b = impl.vertPos_[edge.v1];
    const vec3 p0 = impl.vertPos_[t0];
    const vec3 p1 = impl.vertPos_[t1];
    const vec3 p2 = impl.vertPos_[t2];
    const vec3 e1 = p1 - p0;
    const vec3 e2 = p2 - p0;
    const vec3 n = cross(e1, e2);
    const vec3 d = b - a;
    const double denom = dot(n, d);
    const double nMag2 = dot(n, n);
    if (nMag2 == 0) return;
    // Scale-relative parallel-edge reject (denom has units length^3 =
    // |n| * |d|), matching the rest of the pipeline's relative tolerances.
    if (std::fabs(denom) <=
        kPipelineRelTol * std::sqrt(nMag2) * std::sqrt(dot(d, d)))
      return;
    const double s = -dot(n, a - p0) / denom;
    if (s <= 0.0 || s >= 1.0) return;
    const vec3 pos = a + s * d;
    const vec3 v0 = e1, v1 = e2, v2 = pos - p0;
    const double d00 = dot(v0, v0);
    const double d01 = dot(v0, v1);
    const double d11 = dot(v1, v1);
    const double d20 = dot(v2, v0);
    const double d21 = dot(v2, v1);
    const double bDenom = d00 * d11 - d01 * d01;
    if (bDenom == 0) return;
    const double bv = (d11 * d20 - d01 * d21) / bDenom;
    const double bw = (d00 * d21 - d01 * d20) / bDenom;
    const double bu = 1.0 - bv - bw;
    if (bu <= kBarycentricFloor || bv <= kBarycentricFloor ||
        bw <= kBarycentricFloor)
      return;

    // Snap to the NEAREST existing vert within eps (ties to smallest
    // id - order-independent, unlike first-found). Event identity at
    // the computational scale only: wider radii deform the
    // arrangement and re-pierce (docs/OverlapRemoval.md, Known
    // limitations); arrangement-wide identification is step 9.5's.
    int snapTo = -1;
    double snapBest = std::numeric_limits<double>::infinity();
    auto trySnap = [&](int v) {
      const vec3 dd = pos - impl.vertPos_[v];
      const double d2 = dot(dd, dd);
      if (d2 > snapR2) return;
      if (d2 < snapBest || (d2 == snapBest && snapTo >= 0 && v < snapTo)) {
        snapBest = d2;
        snapTo = v;
      }
    };
    trySnap(edge.v0);
    trySnap(edge.v1);
    trySnap(t0);
    trySnap(t1);
    trySnap(t2);
    for (int v : onEdgeLists[edgeIdx].verts) trySnap(v);
    for (int v : onTriLists[triIdx].verts) trySnap(v);

    out.push_back({static_cast<int>(edgeIdx), static_cast<int>(triIdx), pos, s,
                   vec3(bu, bv, bw), snapTo});
  };
  auto recorder = MakeSimpleRecorder(onCollision);
  auto qf = [&](int i) { return triBoxes[i]; };
  bvh.Collisions(recorder, qf, static_cast<int>(nT));
  return out;
}

ChordEdges GenerateChordEdges(const Manifold::Impl& impl,
                              const std::vector<Edge>& edges,
                              const std::vector<EdgeTriIntersection>& etIsects,
                              double eps) {
  using la::dot;
  ChordEdges r;
  if (etIsects.empty()) return r;
  const double eps2 = eps * eps;
  const int baseId = static_cast<int>(impl.NumVert());

  // Resolve each etIsect to a vert id, deduping new positions
  // against each other within eps.
  r.etIsect2Vert.assign(etIsects.size(), -1);
  std::vector<int>& etIsect2Vert = r.etIsect2Vert;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    if (x.snapTo >= 0) {
      etIsect2Vert[i] = x.snapTo;
      continue;
    }
    int found = -1;
    // First match within eps, not nearest: which of two eps-coincident
    // events allocates the vert is insertion-order-dependent, but the
    // choice is identity-only - step 9.5's arrangement-wide
    // unification (nearest-wins) owns the final representative.
    for (size_t j = 0; j < r.newVertPositions.size(); ++j) {
      const vec3 d = x.position - r.newVertPositions[j];
      if (dot(d, d) <= eps2) {
        found = baseId + static_cast<int>(j);
        break;
      }
    }
    if (found >= 0) {
      etIsect2Vert[i] = found;
    } else {
      etIsect2Vert[i] = baseId + static_cast<int>(r.newVertPositions.size());
      r.newVertPositions.push_back(x.position);
    }
  }

  // Group resolved endpoints by tri-tri pair.
  std::map<std::pair<int, int>, std::set<int>> pairEndpoints;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    const Edge& e = edges[x.edgeIdx];
    const int triA1 = e.halfedgeForward / 3;
    const int triA2 = e.halfedgePaired >= 0 ? e.halfedgePaired / 3 : -1;
    auto add = [&](int t1, int t2) {
      auto key = (t1 < t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
      pairEndpoints[key].insert(etIsect2Vert[i]);
    };
    add(triA1, x.triIdx);
    if (triA2 >= 0) add(triA2, x.triIdx);
  }

  // Emit chords (= 2-endpoint pairs) or interior vert records (= 1-
  // endpoint pairs); drop 0 / >= 3.
  for (auto& [key, endpoints] : pairEndpoints) {
    if (endpoints.size() == 2) {
      auto it = endpoints.begin();
      const int a = *it++;
      const int b = *it;
      const int v0 = std::min(a, b);
      const int v1 = std::max(a, b);
      r.newEdges.push_back({v0, v1, key.first, key.second});
    } else if (endpoints.size() == 1) {
      // Single shared endpoint (edge-tip touch, no through-pierce): no
      // chord edge is emitted and nothing further is recorded.
    } else {
      ++r.droppedTriTriPairsWithBadEndpointCount;
    }
  }
  return r;
}

std::vector<NewEdgeWithExtras> AddInteriorVertsToNewEdges(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<PiercedNewEdge>& newEdges,
    const std::vector<TriVertList>& onTriLists, double eps) {
  using la::dot;
  std::vector<NewEdgeWithExtras> out;
  out.reserve(newEdges.size());
  const int baseId = static_cast<int>(impl.NumVert());
  const double eps2 = eps * eps;

  auto getPos = [&](int id) -> vec3 {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  for (const auto& edge : newEdges) {
    NewEdgeWithExtras nwe{edge, {}, {}};
    const vec3 a = getPos(edge.v0);
    const vec3 b = getPos(edge.v1);
    const vec3 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) {
      out.push_back(std::move(nwe));
      continue;
    }
    std::set<int> seen;
    auto check = [&](int v) {
      if (!seen.insert(v).second) return;
      if (v == edge.v0 || v == edge.v1) return;
      // Candidates come from onTriLists, which hold only original-mesh
      // verts; a chord vert id (>= baseId) here would read out of
      // bounds below (step 9 threads those separately).
      DEBUG_ASSERT(v >= 0 && v < baseId, logicErr,
                   "AddInteriorVertsToNewEdges: candidate must be an "
                   "original-mesh vert");
      const LineProj pr = ProjectToLine(impl.vertPos_[v], a, ab, abLen2);
      if (pr.t <= 0.0 || pr.t >= 1.0) return;
      if (pr.distSq > eps2) return;
      nwe.extraVerts.push_back(v);
      nwe.extraTs.push_back(pr.t);
    };
    for (int v : onTriLists[edge.triA].verts) check(v);
    for (int v : onTriLists[edge.triB].verts) check(v);
    SortVertsByT(nwe.extraVerts, nwe.extraTs);
    out.push_back(std::move(nwe));
  }
  return out;
}

std::vector<std::vector<int>> GroupChordsByFace(
    const std::vector<PiercedNewEdge>& newEdges, int numTri) {
  std::vector<std::vector<int>> byFace(numTri);
  for (size_t i = 0; i < newEdges.size(); ++i) {
    const PiercedNewEdge& e = newEdges[i];
    const int chord = static_cast<int>(i);
    if (e.triA >= 0 && e.triA < numTri) byFace[e.triA].push_back(chord);
    if (e.triB >= 0 && e.triB < numTri && e.triB != e.triA) {
      byFace[e.triB].push_back(chord);
    }
  }
  return byFace;
}

TraceChordResult CoplanarTraceChords(const Manifold::Impl& impl,
                                     const std::vector<Edge>& edges,
                                     const std::vector<int>& halfedge2Edge,
                                     std::vector<vec3> newVertPositions,
                                     double tolerance, double eps) {
  using la::cross;
  using la::dot;
  TraceChordResult out;
  out.newVertPositions = std::move(newVertPositions);
  out.newVertSnapR.assign(out.newVertPositions.size(), eps);
  const int numTri = static_cast<int>(impl.NumTri());
  const int baseId = static_cast<int>(impl.NumVert());
  if (numTri < 2) return out;
  auto triVert = [&](int t, int k) { return impl.halfedge_.Start(3 * t + k); };

  // Broad phase: tri-box self-collisions; pairs processed in
  // ascending (a, b) order so new-vert allocation is deterministic.
  // Boxes are eps-padded: a coplanar pair offset by up to eps along
  // the normal passes the plane gate below, but their unpadded
  // (zero-thickness) boxes would never overlap.
  std::vector<Box> triBoxes(numTri);
  const vec3 pad(eps, eps, eps);
  for (int t = 0; t < numTri; ++t) {
    Box b(impl.vertPos_[triVert(t, 0)], impl.vertPos_[triVert(t, 1)]);
    b.Union(impl.vertPos_[triVert(t, 2)]);
    b.min -= pad;
    b.max += pad;
    triBoxes[t] = b;
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));
  std::vector<std::pair<int, int>> pairs;
  auto recordPair = [&](size_t qi, size_t li) {
    if (qi >= li) return;
    const int ta = static_cast<int>(bvh.leaf2Orig[qi]);
    const int tb = static_cast<int>(bvh.leaf2Orig[li]);
    pairs.push_back({std::min(ta, tb), std::max(ta, tb)});
  };
  auto recorder = MakeSimpleRecorder(recordPair);
  auto qf = [&](int i) { return bvh.boxes[i]; };
  bvh.Collisions(recorder, qf, numTri);
  std::sort(pairs.begin(), pairs.end());
  pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

  std::set<std::pair<int, int>> additionSeen;  // (edge, vertId)
  for (const std::pair<int, int>& facePair : pairs) {
    // Locals, not structured bindings: the lambdas below capture them.
    const int fa = facePair.first;
    const int fb = facePair.second;
    const vec3 av[3] = {impl.vertPos_[triVert(fa, 0)],
                        impl.vertPos_[triVert(fa, 1)],
                        impl.vertPos_[triVert(fa, 2)]};
    const vec3 bv[3] = {impl.vertPos_[triVert(fb, 0)],
                        impl.vertPos_[triVert(fb, 1)],
                        impl.vertPos_[triVert(fb, 2)]};
    // Plane gate: every vert of the smaller face within eps of the
    // LARGER face's plane (a near-zero-area sliver's own plane is
    // noise; the larger plane is the trustworthy one).
    const vec3 nA = cross(av[1] - av[0], av[2] - av[0]);
    const vec3 nB = cross(bv[1] - bv[0], bv[2] - bv[0]);
    const double areaA2 = dot(nA, nA);
    const double areaB2 = dot(nB, nB);
    const bool aLarger = areaA2 >= areaB2;
    const vec3 nGRaw = aLarger ? nA : nB;
    const double nGLen = std::sqrt(dot(nGRaw, nGRaw));
    if (nGLen == 0) continue;  // degenerate face: nothing to trace
    const vec3 nG = nGRaw / nGLen;
    const vec3 gOrigin = aLarger ? av[0] : bv[0];
    const vec3* other = aLarger ? bv : av;
    bool coplanar = true;
    for (int k = 0; k < 3 && coplanar; ++k) {
      coplanar = std::fabs(dot(other[k] - gOrigin, nG)) <= eps;
    }
    if (!coplanar) continue;

    // Shared 2D frame: the LARGER face's basis (the gate already
    // trusts its plane; a sliver's own normal is noise), tie to fa
    // (the lower id) on equal areas - deterministic. One frame means
    // a crossing is computed once and shares its id across both clip
    // directions.
    const InPlaneBasis basis = FaceBasisFromNormal(nG);
    const vec3 origin = gOrigin;
    auto to2d = [&](const vec3& p) {
      const vec3 d = p - origin;
      return vec2(dot(d, basis.u), dot(d, basis.v));
    };
    vec2 a2[3], b2[3];
    for (int k = 0; k < 3; ++k) {
      a2[k] = to2d(av[k]);
      b2[k] = to2d(bv[k]);
    }

    // CANCELLATION-HAZARD flag (taken BEFORE interval filtering so even
    // equal-boundary coincidences whose intervals all ride the boundary
    // still count). Two conditions must both hold:
    // 1. Winding normals ANTI-ALIGNED: dot(normalize(nA), normalize(nB)) < 0.
    //    Same-oriented coplanar neighbors (every flat face's own tris, which
    //    the plane gate also pairs) have dot > 0 and are NOT flagged:
    //    flagging them would re-gate the boss-through-plate fixture.
    // 2. In-plane bounding boxes OVERLAP: if the bbox of fa's projected
    //    triangle and the bbox of fb's projected triangle do not overlap in
    //    2D, there can be no interior coincidence.
    if (areaA2 > 0 && areaB2 > 0 && dot(nA, nB) < 0) {
      // Check 2D bbox overlap (separating-axis test on the two projected
      // triangles' bounding boxes).
      double aMinU = a2[0].x, aMaxU = a2[0].x;
      double aMinV = a2[0].y, aMaxV = a2[0].y;
      double bMinU = b2[0].x, bMaxU = b2[0].x;
      double bMinV = b2[0].y, bMaxV = b2[0].y;
      for (int k = 1; k < 3; ++k) {
        aMinU = std::min(aMinU, a2[k].x);
        aMaxU = std::max(aMaxU, a2[k].x);
        aMinV = std::min(aMinV, a2[k].y);
        aMaxV = std::max(aMaxV, a2[k].y);
        bMinU = std::min(bMinU, b2[k].x);
        bMaxU = std::max(bMaxU, b2[k].x);
        bMinV = std::min(bMinV, b2[k].y);
        bMaxV = std::max(bMaxV, b2[k].y);
      }
      if (aMinU <= bMaxU && bMinU <= aMaxU && aMinV <= bMaxV &&
          bMinV <= aMaxV) {
        out.coplanarHazardFaces.insert(fa);
        out.coplanarHazardFaces.insert(fb);
      }
    }

    // Orientation sign per projected tri (the opposite sheet projects
    // CW), so the clip's inside test works for both.
    auto orient2 = [](const vec2* t3) {
      return (t3[1].x - t3[0].x) * (t3[2].y - t3[0].y) -
             (t3[1].y - t3[0].y) * (t3[2].x - t3[0].x);
    };
    // In-plane signed distance of x from tri edge k, positive inward.
    auto edgeDist = [&](const vec2* t3, double orientSign, int k,
                        const vec2& x) {
      const vec2 e0 = t3[k];
      const vec2 e1 = t3[(k + 1) % 3];
      const vec2 d(e1.x - e0.x, e1.y - e0.y);
      const double len = std::sqrt(d.x * d.x + d.y * d.y);
      if (len == 0) return 0.0;
      const double s = d.x * (x.y - e0.y) - d.y * (x.x - e0.x);
      return orientSign * s / len;
    };
    auto interiorMargin = [&](const vec2* t3, double orientSign,
                              const vec2& x) {
      double m = std::numeric_limits<double>::infinity();
      for (int k = 0; k < 3; ++k) {
        m = std::min(m, edgeDist(t3, orientSign, k, x));
      }
      return m;
    };

    // Clip one segment (2D p -> q) against one projected tri; on
    // success fills [tEnter, tExit] plus the clipping tri edge index
    // per end (-1 = the segment's own endpoint).
    auto clip = [&](const vec2& p, const vec2& q, const vec2* t3,
                    double orientSign, double& tEnter, double& tExit,
                    int& kEnter, int& kExit) {
      tEnter = 0.0;
      tExit = 1.0;
      kEnter = -1;
      kExit = -1;
      for (int k = 0; k < 3; ++k) {
        const double d0 = edgeDist(t3, orientSign, k, p);
        const double d1 = edgeDist(t3, orientSign, k, q);
        if (d0 < 0 && d1 < 0) return false;
        if (d0 < 0 || d1 < 0) {
          const double tc = d0 / (d0 - d1);
          if (d0 < 0) {  // entering
            if (tc > tEnter) {
              tEnter = tc;
              kEnter = k;
            }
          } else {  // exiting
            if (tc < tExit) {
              tExit = tc;
              kExit = k;
            }
          }
        }
      }
      return tEnter < tExit;
    };

    const double orientA = orient2(a2) >= 0 ? 1.0 : -1.0;
    const double orientB = orient2(b2) >= 0 ? 1.0 : -1.0;
    std::set<std::pair<int, int>> pairChords;  // sorted (v0, v1)

    // One clip direction: edges of `src` face against the `dst` tri.
    auto traceDirection = [&](int srcFace, const vec3* srcV, const vec2* src2,
                              int dstFace, const vec3* dstV, const vec2* dst2,
                              double dstOrient) {
      for (int k = 0; k < 3; ++k) {
        const vec2 p = src2[k];
        const vec2 q = src2[(k + 1) % 3];
        double t0, t1;
        int k0, k1;
        if (!clip(p, q, dst2, dstOrient, t0, t1, k0, k1)) continue;
        const vec3 P0 = srcV[k] + t0 * (srcV[(k + 1) % 3] - srcV[k]);
        const vec3 P1 = srcV[k] + t1 * (srcV[(k + 1) % 3] - srcV[k]);
        // Qualification: length > eps; midpoint interior to the dst
        // face by > eps (full-through cuts qualify - their midpoints
        // are interior; boundary-riding intervals never do).
        const vec3 d3 = P1 - P0;
        if (std::sqrt(dot(d3, d3)) <= eps) {
          ++out.intervalsRejected;
          continue;
        }
        const vec2 mid2((p.x + (q.x - p.x) * 0.5 * (t0 + t1)),
                        (p.y + (q.y - p.y) * 0.5 * (t0 + t1)));
        if (interiorMargin(dst2, dstOrient, mid2) <= eps) {
          ++out.intervalsRejected;
          continue;
        }
        // Resolve the two endpoint ids. t == 0/1: the src edge's own
        // vert. Crossings: snap to the nearest of the pair's six
        // corners at max(tolerance + eps, the conditioned radius)
        // (ties to smallest id - the step-9 convention), else
        // allocate, deduping new-to-new over the whole pool at the
        // SOURCE-GATED radius (see below).
        int ids[2];
        int clipK[2] = {k0, k1};
        vec3 pos3[2] = {P0, P1};
        bool grazeReject = false;
        for (int e = 0; e < 2 && !grazeReject; ++e) {
          const double t = e == 0 ? t0 : t1;
          if (t <= 0.0) {
            ids[e] = triVert(srcFace, k);
            continue;
          }
          if (t >= 1.0) {
            ids[e] = triVert(srcFace, (k + 1) % 3);
            continue;
          }
          // Grazing guard: the lifted crossing must sit within eps of
          // BOTH original 3D edges it claims to lie on.
          const vec3 x3 = pos3[e];
          auto distToSeg = [&](const vec3& s0, const vec3& s1) {
            const vec3 d = s1 - s0;
            const double len2 = dot(d, d);
            const double tt =
                len2 > 0 ? std::clamp(dot(x3 - s0, d) / len2, 0.0, 1.0) : 0.0;
            const vec3 c = s0 + tt * d - x3;
            return std::sqrt(dot(c, c));
          };
          if (distToSeg(srcV[k], srcV[(k + 1) % 3]) > eps ||
              (clipK[e] >= 0 &&
               distToSeg(dstV[clipK[e]], dstV[(clipK[e] + 1) % 3]) > eps)) {
            grazeReject = true;
            break;
          }
          // Conditioning of THIS crossing (see ConditionedSnapRadius).
          // Kept SEPARATE from the snap base: the tolerance + eps
          // new-to-old term must never widen the new-to-new dedup
          // below.
          double condOnly = eps;
          if (clipK[e] >= 0) {
            const vec2 d1(q.x - p.x, q.y - p.y);
            const vec2 d2(dst2[(clipK[e] + 1) % 3].x - dst2[clipK[e]].x,
                          dst2[(clipK[e] + 1) % 3].y - dst2[clipK[e]].y);
            condOnly = ConditionedSnapRadius(d1, d2, eps);
          }
          // Corner snap (new-to-old): the step-9 base radius, widened
          // by the conditioning; nearest, ties to smallest id.
          int best = -1;
          double bestD = std::max(tolerance + eps, condOnly);
          for (int c = 0; c < 6; ++c) {
            const int vid = c < 3 ? triVert(fa, c) : triVert(fb, c - 3);
            const vec3 dv = (c < 3 ? av[c] : bv[c - 3]) - x3;
            const double dd = std::sqrt(dot(dv, dv));
            // <= on the first hit: "within" the radius is inclusive
            // (a corner exactly at the boundary still snaps).
            if (dd < bestD || (dd == bestD && (best < 0 || vid < best))) {
              bestD = dd;
              best = vid;
            }
          }
          if (best >= 0) {
            ids[e] = best;
            continue;
          }
          // New-to-new dedup, SOURCE-GATED: the conditioned radius
          // applies only between entries that are themselves
          // conditioned (the min of the two claims) - a nearby
          // well-conditioned vert (a step-7 pierce, recorded at eps)
          // is a geometrically DISTINCT point, and absorbing it would
          // weld unrelated arrangement features. A missed twin fails
          // safe as a rim; a wrong weld fails silent.
          int found = -1;
          for (size_t j = 0; j < out.newVertPositions.size(); ++j) {
            const double matchR =
                std::max(eps, std::min(condOnly, out.newVertSnapR[j]));
            const vec3 dv = out.newVertPositions[j] - x3;
            if (std::sqrt(dot(dv, dv)) <= matchR) {
              found = baseId + static_cast<int>(j);
              break;
            }
          }
          if (found >= 0) {
            ids[e] = found;
            // The widest conditioning claim wins (feeds step 9.5's
            // new-onto-original snap).
            out.newVertSnapR[found - baseId] =
                std::max(out.newVertSnapR[found - baseId], condOnly);
          } else {
            ids[e] = baseId + static_cast<int>(out.newVertPositions.size());
            out.newVertPositions.push_back(x3);
            out.newVertSnapR.push_back(condOnly);
          }
        }
        if (grazeReject) {
          ++out.intervalsRejected;
          continue;
        }
        if (ids[0] == ids[1]) {  // both snapped to one corner
          ++out.intervalsRejected;
          continue;
        }
        // On-edge additions: an interior crossing endpoint lies on the
        // src edge, and - when a tri edge clipped it - on the dst
        // face's edge too (the X case: one record per edge, one vert
        // id). Emitted for SNAPPED endpoints as well as allocated ones
        // (a corner-snapped crossing still
        // subdivides the edges it crossed, or the claiming faces'
        // partitions never see the cut) - unless the resolved id IS
        // that edge's endpoint, where no subdivision is needed. The
        // resolved POSITION is used for t: a snap moves the point, and
        // a projection landing outside (0, 1) means the snap target
        // sits past the edge's end - the endpoint case again. Emitted
        // BEFORE the both-directions chord dedup below: each direction
        // claims a DIFFERENT src edge, and a snapped endpoint's only
        // subdividing record can come from the second direction
        // (additionSeen dedups per (edge, id)).
        for (int e = 0; e < 2; ++e) {
          if (ids[e] == triVert(srcFace, k) ||
              ids[e] == triVert(srcFace, (k + 1) % 3)) {
            continue;  // t == 0/1 endpoints resolve to the edge's verts
          }
          const vec3 resolvedPos = ids[e] < baseId
                                       ? impl.vertPos_[ids[e]]
                                       : out.newVertPositions[ids[e] - baseId];
          // A corner-snapped id can sit up to the conditioned corner
          // radius OFF this edge; threading it here moves the edge's
          // polyline by that much - the same displacement the snap
          // itself accepted (the conditioned band of Known
          // limitations #1), with ids consistent on both faces. A
          // snap to a face's OPPOSITE corner (not on the claimed
          // edge) is excluded: on an obtuse near-degenerate face its
          // projection can land in (0, 1) yet the vert is a whole
          // edge away.
          auto addOn = [&](int face, int kk, const vec3& s0, const vec3& s1) {
            if (kk < 0) return;
            if (ids[e] == triVert(face, (kk + 2) % 3)) return;
            const int edgeIdx = halfedge2Edge[3 * face + kk];
            if (edgeIdx < 0) return;
            if (ids[e] == edges[edgeIdx].v0 || ids[e] == edges[edgeIdx].v1) {
              return;
            }
            const vec3 d = s1 - s0;
            const double len2 = dot(d, d);
            const double tt = len2 > 0 ? dot(resolvedPos - s0, d) / len2 : 0.0;
            if (tt <= 0.0 || tt >= 1.0) return;
            if (!additionSeen.insert({edgeIdx, ids[e]}).second) return;
            // t along the canonical edge direction (v0 -> v1).
            const double tEdge =
                edges[edgeIdx].v0 == impl.halfedge_.Start(3 * face + kk)
                    ? tt
                    : 1.0 - tt;
            out.onEdgeAdditions.push_back({edgeIdx, ids[e], tEdge});
          };
          addOn(srcFace, k, srcV[k], srcV[(k + 1) % 3]);
          if (clipK[e] >= 0) {  // arg exprs index dstV: guard before call
            addOn(dstFace, clipK[e], dstV[clipK[e]], dstV[(clipK[e] + 1) % 3]);
          }
        }
        const std::pair<int, int> key{std::min(ids[0], ids[1]),
                                      std::max(ids[0], ids[1])};
        if (!pairChords.insert(key).second) continue;  // both directions
        PiercedNewEdge chord;
        chord.v0 = key.first;
        chord.v1 = key.second;
        chord.triA = fa;
        chord.triB = fb;
        out.chords.push_back(chord);
      }
    };
    traceDirection(fb, bv, b2, fa, av, a2, orientA);
    traceDirection(fa, av, a2, fb, bv, b2, orientB);
  }
  return out;
}

void AddVertsToOnEdgeLists(const std::vector<OnEdgeAddition>& additions,
                           std::vector<EdgeVertList>& onEdgeLists) {
  // Id-dedup against the existing list, then re-sort each touched
  // edge by t (the PropagateNewVertsToOnEdgeLists internals; that
  // sibling's interface is parallel to etIsects and cannot carry
  // these explicit additions).
  std::vector<int> touched;
  for (const OnEdgeAddition& a : additions) {
    EdgeVertList& list = onEdgeLists[a.edge];
    if (std::find(list.verts.begin(), list.verts.end(), a.vertId) !=
        list.verts.end()) {
      continue;
    }
    list.verts.push_back(a.vertId);
    list.ts.push_back(a.t);
    touched.push_back(a.edge);
  }
  std::sort(touched.begin(), touched.end());
  touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
  for (const int e : touched) {
    SortVertsByT(onEdgeLists[e].verts, onEdgeLists[e].ts);
  }
}

UnifyResult UnifyArrangementVerts(const Manifold::Impl& impl,
                                  const std::vector<vec3>& newVertPositions,
                                  const std::vector<Edge>& edges,
                                  std::vector<EdgeVertList>& onEdgeLists,
                                  std::vector<NewEdgeWithExtras>& chords,
                                  double eps,
                                  const std::vector<double>& perVertSnapR) {
  using la::dot;
  const int baseId = static_cast<int>(impl.NumVert());
  const int nNew = static_cast<int>(newVertPositions.size());
  if (nNew == 0) return {};
  // New-new pairs unite at the nearby-crossing merge radius (10 *
  // eps - frame-to-frame spread of one computed point). New verts
  // snap onto nearby originals at the same radius, widened PER VERT
  // by its conditioned allocation radius (TraceChordResult::
  // newVertSnapR) - blanket widening was tried and rejected: it
  // rounded real geometry into corners and re-pierced.
  const double radius = 10.0 * eps;
  const double radius2 = radius * radius;
  auto vertSnapR = [&](int j) {
    const double cond =
        static_cast<size_t>(j) < perVertSnapR.size() ? perVertSnapR[j] : 0.0;
    return std::max(radius, cond);
  };
  double maxSnapR = radius;
  for (int j = 0; j < nNew; ++j) maxSnapR = std::max(maxSnapR, vertSnapR(j));

  // Broad phase over the new verts; new-new self-collisions unite,
  // and each ORIGINAL vert within a new vert's snap radius becomes a
  // snap candidate for its cluster (nearest wins, ties to smallest).
  // One BVH pass padded at the GLOBAL maxSnapR serves both phases: for
  // ill-conditioned inputs (snapR up to kCondSnapCapEps * eps) the
  // new-new narrow phase (10 * eps) sees over-broad candidates - a
  // bounded, accepted perf cost; splitting into per-phase BVHs would
  // re-arm verification of this pass for zero behavioral gain.
  std::vector<Box> newBoxes(nNew);
  const double half = 0.5 * maxSnapR;
  for (int i = 0; i < nNew; ++i) {
    const vec3& p = newVertPositions[i];
    newBoxes[i] = Box(vec3(p.x - half, p.y - half, p.z - half),
                      vec3(p.x + half, p.y + half, p.z + half));
  }
  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(newBoxes.data(), newBoxes.size()));
  DisjointSets uf(static_cast<uint32_t>(nNew));
  auto unitePair = [&](size_t qi, size_t li) {
    if (qi >= li) return;
    const int a = static_cast<int>(bvh.leaf2Orig[qi]);
    const int b = static_cast<int>(bvh.leaf2Orig[li]);
    const vec3 d = newVertPositions[a] - newVertPositions[b];
    if (dot(d, d) <= radius2) {
      uf.unite(static_cast<uint32_t>(a), static_cast<uint32_t>(b));
    }
  };
  auto recorder = MakeSimpleRecorder(unitePair);
  auto qf = [&](int i) { return bvh.boxes[i]; };
  bvh.Collisions(recorder, qf, nNew);

  // [new idx] -> nearest original id within the vert's snap radius
  // (ties to smallest id - the pipeline's snap convention; an
  // id-priority pick could jump past the adjacent corner to a far
  // small-id vert under a wide conditioned radius).
  std::vector<int> snapTo(nNew, -1);
  std::vector<double> snapD2(nNew, std::numeric_limits<double>::infinity());
  {
    std::vector<Box> origBoxes(baseId);  // empty boxes overlap nothing
    for (int v = 0; v < baseId; ++v) {
      const vec3& p = impl.vertPos_[v];
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) {
        continue;  // unreferenced post-merge slots
      }
      origBoxes[v] = Box(vec3(p.x - half, p.y - half, p.z - half),
                         vec3(p.x + half, p.y + half, p.z + half));
    }
    auto hit = [&](size_t qi, size_t li) {
      const int v = static_cast<int>(qi);
      const int j = static_cast<int>(bvh.leaf2Orig[li]);
      const vec3 d = newVertPositions[j] - impl.vertPos_[v];
      const double r = vertSnapR(j);
      const double d2 = dot(d, d);
      if (d2 > r * r) return;
      if (d2 < snapD2[j] || (d2 == snapD2[j] && v < snapTo[j])) {
        snapD2[j] = d2;
        snapTo[j] = v;
      }
    };
    auto rec = MakeSimpleRecorder(hit);
    auto qfOrig = [&](int i) { return origBoxes[i]; };
    bvh.Collisions(rec, qfOrig, baseId);
  }

  // Cluster reps: the NEAREST snap-original anywhere in the cluster
  // (each member's snapTo is already its own nearest; across members
  // compare those distances, ties to the smaller id - an id-priority
  // pick could choose a far small-id original over an adjacent
  // corner), else baseId + the smallest member index.
  std::vector<int> repOriginal(nNew, -1);
  std::vector<double> repOriginalD2(nNew,
                                    std::numeric_limits<double>::infinity());
  std::vector<int> repMember(nNew, -1);
  for (int i = 0; i < nNew; ++i) {
    const int r = static_cast<int>(uf.find(static_cast<uint32_t>(i)));
    if (repMember[r] < 0 || i < repMember[r]) repMember[r] = i;
    if (snapTo[i] >= 0 &&
        (snapD2[i] < repOriginalD2[r] ||
         (snapD2[i] == repOriginalD2[r] && snapTo[i] < repOriginal[r]))) {
      repOriginal[r] = snapTo[i];
      repOriginalD2[r] = snapD2[i];
    }
  }
  std::vector<int> remap(baseId + nNew);
  std::iota(remap.begin(), remap.end(), 0);
  int changed = 0;
  double maxMove = 0.0;
  for (int i = 0; i < nNew; ++i) {
    const int r = static_cast<int>(uf.find(static_cast<uint32_t>(i)));
    const int rep =
        repOriginal[r] >= 0 ? repOriginal[r] : baseId + repMember[r];
    if (rep != baseId + i) {
      ++changed;
      // The consumers of id baseId + i now read rep's position: the
      // arrangement point moved by this distance (tolerance claim).
      const vec3 repPos =
          rep < baseId ? impl.vertPos_[rep] : newVertPositions[rep - baseId];
      const vec3 d = newVertPositions[i] - repPos;
      maxMove = std::max(maxMove, std::sqrt(dot(d, d)));
    }
    remap[baseId + i] = rep;
  }
  if (changed == 0) return {};

  // Remap consumers. Chord endpoints first; extras then dedup by id
  // and drop ids that became an endpoint. Ts are RECOMPUTED from the
  // remapped positions and re-sorted (a remap moves
  // the consumed position by up to the merge radius, and a stale t
  // order would make the partition build a crossed sub-edge
  // sequence).
  auto posOf = [&](int id) {
    return id < baseId ? impl.vertPos_[id] : newVertPositions[id - baseId];
  };
  for (NewEdgeWithExtras& nwe : chords) {
    nwe.edge.v0 = remap[nwe.edge.v0];
    nwe.edge.v1 = remap[nwe.edge.v1];
    std::vector<int> vs;
    std::vector<double> ts;
    std::set<int> seen;
    const vec3 c0 = posOf(nwe.edge.v0);
    const vec3 cd = posOf(nwe.edge.v1) - c0;
    const double cLen2 = dot(cd, cd);
    for (size_t i = 0; i < nwe.extraVerts.size(); ++i) {
      const int v = remap[nwe.extraVerts[i]];
      if (v == nwe.edge.v0 || v == nwe.edge.v1) continue;
      if (!seen.insert(v).second) continue;
      const double t =
          cLen2 > 0 ? dot(posOf(v) - c0, cd) / cLen2 : nwe.extraTs[i];
      // A representative can land PAST a chord endpoint (the nearest
      // original beyond it); outside (0, 1) it subdivides nothing -
      // the interior rule every on-edge/on-chord builder applies.
      if (t <= 0.0 || t >= 1.0) continue;
      vs.push_back(v);
      ts.push_back(t);
    }
    SortVertsByT(vs, ts);
    nwe.extraVerts = std::move(vs);
    nwe.extraTs = std::move(ts);
  }
  for (size_t e = 0; e < onEdgeLists.size(); ++e) {
    EdgeVertList& list = onEdgeLists[e];
    if (list.verts.empty()) continue;
    std::vector<int> vs;
    std::vector<double> ts;
    std::set<int> seen;
    const vec3 e0 = impl.vertPos_[edges[e].v0];
    const vec3 ed = impl.vertPos_[edges[e].v1] - e0;
    const double eLen2 = dot(ed, ed);
    for (size_t i = 0; i < list.verts.size(); ++i) {
      const int v = remap[list.verts[i]];
      if (v == edges[e].v0 || v == edges[e].v1) continue;
      if (!seen.insert(v).second) continue;
      const double t = eLen2 > 0 ? dot(posOf(v) - e0, ed) / eLen2 : list.ts[i];
      if (t <= 0.0 || t >= 1.0) continue;  // representative past an endpoint
      vs.push_back(v);
      ts.push_back(t);
    }
    SortVertsByT(vs, ts);
    list.verts = std::move(vs);
    list.ts = std::move(ts);
  }
  return {changed, maxMove};
}

std::vector<OnChordContact> FindOnChordEndpointContacts(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& face2Chords, double tolerance,
    double eps) {
  using la::dot;
  std::vector<OnChordContact> out;
  const int baseId = static_cast<int>(impl.NumVert());
  // New-to-old snap: prior drift plus current-op error, matching
  // boolean2's newToOldThresh (tolerance + eps), not bare eps.
  const double snap = tolerance + eps;
  const double snap2 = snap * snap;
  // A chord pair shares up to two faces; record each (chord, vert)
  // contact once.
  std::set<std::pair<int, int>> seen;
  for (const std::vector<int>& faceChords : face2Chords) {
    for (int ci : faceChords) {
      const PiercedNewEdge& c = chords[ci].edge;
      for (int di : faceChords) {
        if (di == ci) continue;
        const PiercedNewEdge& d = chords[di].edge;
        const vec3 a = GetPos3(d.v0, baseId, impl, newVertPositions);
        const vec3 b = GetPos3(d.v1, baseId, impl, newVertPositions);
        const vec3 ab = b - a;
        const double abLen2 = dot(ab, ab);
        if (abLen2 == 0) continue;
        // Endpoint-proximity zone in t-space: a contact within snap of
        // d's own endpoints is an endpoint-near-endpoint case, not an
        // interior vert (it would re-create the near-line sliver the
        // guard exists to block). snap/len >= 0.5 excludes the whole
        // chord (shorter than 2*snap).
        const double tGuard = snap / std::sqrt(abLen2);
        for (const int e : {c.v0, c.v1}) {
          if (e == d.v0 || e == d.v1) continue;
          const vec3 p = GetPos3(e, baseId, impl, newVertPositions);
          const LineProj pr = ProjectToLine(p, a, ab, abLen2);
          if (pr.t <= tGuard || pr.t >= 1.0 - tGuard) continue;
          if (pr.distSq > snap2) continue;
          if (!seen.insert({di, e}).second) continue;
          out.push_back({di, e, pr.t});
        }
      }
    }
  }
  // Deterministic output order: ascending chord, then ascending t.
  std::sort(out.begin(), out.end(),
            [](const OnChordContact& x, const OnChordContact& y) {
              if (x.chord != y.chord) return x.chord < y.chord;
              return x.t < y.t;
            });
  return out;
}

std::vector<ChordChordCrossing> FindChordChordCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<std::vector<int>>& face2Chords,
    VecView<const vec3> faceNormals, double eps) {
  using la::dot;
  std::vector<ChordChordCrossing> out;
  const int baseId = static_cast<int>(impl.NumVert());
  // A chord pair shares up to two faces; its crossing is recorded once
  // (lowest face wins by iteration order).
  std::set<std::pair<int, int>> seenPairs;
  for (size_t face = 0; face < face2Chords.size(); ++face) {
    const std::vector<int>& faceChords = face2Chords[face];
    if (faceChords.size() < 2) continue;
    DEBUG_ASSERT(face < faceNormals.size(), logicErr,
                 "FindChordChordCrossings: face normal missing");
    if (face >= faceNormals.size()) continue;
    const vec3 nRaw = faceNormals[face];
    const double nLen2 = dot(nRaw, nRaw);
    if (nLen2 == 0) continue;
    const vec3 n = nRaw / std::sqrt(nLen2);
    const InPlaneBasis basis = FaceBasisFromNormal(n);
    const vec3 u = basis.u;
    const vec3 v = basis.v;
    for (size_t i = 0; i < faceChords.size(); ++i) {
      for (size_t j = i + 1; j < faceChords.size(); ++j) {
        const int ci = faceChords[i];
        const int cj = faceChords[j];
        if (!seenPairs.insert({std::min(ci, cj), std::max(ci, cj)}).second) {
          continue;
        }
        const PiercedNewEdge& ea = chords[ci].edge;
        const PiercedNewEdge& eb = chords[cj].edge;
        // Re-project endpoints onto the face plane (through ea.v0)
        // before the 2D mapping, so out-of-plane drift from prior
        // merges is zero by construction.
        const vec3 planePt = GetPos3(ea.v0, baseId, impl, newVertPositions);
        auto to2d = [&](int id) -> vec2 {
          const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
          const vec3 inPlane = p - dot(p - planePt, n) * n;
          return vec2(dot(inPlane - planePt, u), dot(inPlane - planePt, v));
        };
        const vec2 a0 = to2d(ea.v0);
        const vec2 a1 = to2d(ea.v1);
        const vec2 b0 = to2d(eb.v0);
        const vec2 b1 = to2d(eb.v1);
        // stableEdgeId: the chord's index - deterministic (newEdges
        // order comes from GenerateChordEdges' sorted-pair map walk),
        // NOT broad-phase pair order.
        const boolean2::GraphSegment2D segA{a0, a1, ci};
        const boolean2::GraphSegment2D segB{b0, b1, cj};
        vec2 p2;
        if (!boolean2::IntersectSegments(segA, segB, eps, &p2)) continue;
        const vec2 da = a1 - a0;
        const vec2 db = b1 - b0;
        const double tA = dot(p2 - a0, da) / dot(da, da);
        const double tB = dot(p2 - b0, db) / dot(db, db);
        const vec3 pos = planePt + p2.x * u + p2.y * v;
        out.push_back({pos, ci, cj, tA, tB, static_cast<int>(face),
                       ConditionedSnapRadius(da, db, eps)});
      }
    }
  }
  return out;
}

Step9Threading ResolveAndThreadClusters(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions, std::vector<double> newVertSnapR,
    const std::vector<ChordCrossing>& clusters,
    const std::vector<OnChordContact>& contacts, double tolerance, double eps) {
  using la::dot;
  const int baseId = static_cast<int>(impl.NumVert());
  const double snap = tolerance + eps;
  const double snap2 = snap * snap;
  // Keep the per-vert conditioned radii parallel to the pool: eps for
  // any pool entry that arrived without one (defensive; the trace pass
  // emits a full-length vector).
  newVertSnapR.resize(newVertPositions.size(), eps);
  // Canonical-id resolution: resolve-then-allocate, per cluster,
  // symmetric across ALL incident chords. Candidates are every
  // existing vert the crossing could BE - chord endpoints,
  // already-threaded on-chord verts, and the pass-0 contacts (not yet
  // threaded at this point) - within tolerance + eps. Nearest wins;
  // ties take the smallest id. Only when no candidate exists is a
  // fresh vert allocated. This is what prevents a crossing threading
  // as an endpoint id on one chord and a fresh id on another (the
  // split-identity bug).
  std::vector<ChordCrossing> crossings;
  crossings.reserve(clusters.size());
  for (const ChordCrossing& cl : clusters) {
    int best = -1;
    double bestD2 = 0.0;
    auto consider = [&](int id) {
      const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
      const vec3 dv = cl.pos - p;
      const double d2 = dot(dv, dv);
      if (d2 > snap2) return;
      if (best < 0 || d2 < bestD2 || (d2 == bestD2 && id < best)) {
        best = id;
        bestD2 = d2;
      }
    };
    for (const int side : cl.chords) {
      const PiercedNewEdge& e = chords[side].edge;
      consider(e.v0);
      consider(e.v1);
      for (const int v : chords[side].extraVerts) consider(v);
    }
    for (const OnChordContact& c : contacts) {
      if (std::find(cl.chords.begin(), cl.chords.end(), c.chord) !=
          cl.chords.end()) {
        consider(c.vertId);
      }
    }
    int id = best;
    if (id < 0) {
      id = baseId + static_cast<int>(newVertPositions.size());
      newVertPositions.push_back(cl.pos);
      // A fresh crossing vert carries its cluster's conditioned
      // radius (eps-floored) into step 9.5's new-onto-original snap.
      newVertSnapR.push_back(std::max(eps, cl.snapR));
    } else if (id >= baseId) {
      // Snapping onto an existing NEW vert: that vert now also
      // stands for this ill-conditioned crossing - widen (the trace
      // pass's dedup rule, step-9 side).
      newVertSnapR[id - baseId] = std::max(newVertSnapR[id - baseId], cl.snapR);
    }
    ChordCrossing rec = cl;
    rec.id = id;
    crossings.push_back(std::move(rec));
  }
  // Threading: per chord, the unified record list (existing extras +
  // pass-0 contacts + resolved crossings), with every t RECOMPUTED
  // from the resolved position (a snap can move the vertex by up to
  // tolerance + eps - enough to reorder a stale t-sort), the pass-0
  // endpoint-zone guard re-applied, id-dedup over the unified list,
  // then t-sort with an eps/len dedup backstop.
  std::vector<std::vector<int>> pending(chords.size());
  for (const OnChordContact& c : contacts) {
    pending[c.chord].push_back(c.vertId);
  }
  for (const ChordCrossing& cc : crossings) {
    for (const int ch : cc.chords) pending[ch].push_back(cc.id);
  }
  for (size_t ci = 0; ci < chords.size(); ++ci) {
    if (pending[ci].empty()) continue;
    NewEdgeWithExtras& nwe = chords[ci];
    std::vector<int> ids = nwe.extraVerts;
    ids.insert(ids.end(), pending[ci].begin(), pending[ci].end());
    const vec3 a = GetPos3(nwe.edge.v0, baseId, impl, newVertPositions);
    const vec3 b = GetPos3(nwe.edge.v1, baseId, impl, newVertPositions);
    const vec3 ab = b - a;
    const double abLen2 = dot(ab, ab);
    if (abLen2 == 0) continue;
    const double len = std::sqrt(abLen2);
    const double tGuard = snap / len;
    const double tDedup = eps / len;
    // Pre-existing step-8 extras were admitted under step 8's weaker
    // t-guard; the step-9 endpoint-zone guard applies only to NEWLY
    // added records (contacts + crossings), so a rebuild does not
    // silently drop legitimate near-endpoint on-tri verts.
    const std::set<int> preexisting(nwe.extraVerts.begin(),
                                    nwe.extraVerts.end());
    std::vector<std::pair<double, int>> recs;
    recs.reserve(ids.size());
    std::set<int> seenIds;
    for (const int id : ids) {
      if (id == nwe.edge.v0 || id == nwe.edge.v1) continue;
      if (!seenIds.insert(id).second) continue;
      const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
      const double t = dot(p - a, ab) / abLen2;
      const bool isNew = preexisting.count(id) == 0;
      if (isNew && (t <= tGuard || t >= 1.0 - tGuard)) continue;
      recs.push_back({t, id});
    }
    std::sort(recs.begin(), recs.end());
    std::vector<int> outV;
    std::vector<double> outT;
    outV.reserve(recs.size());
    outT.reserve(recs.size());
    for (const auto& [t, id] : recs) {
      if (!outT.empty() && t - outT.back() <= tDedup) continue;
      outT.push_back(t);
      outV.push_back(id);
    }
    nwe.extraVerts = std::move(outV);
    nwe.extraTs = std::move(outT);
  }
  return {std::move(chords), std::move(newVertPositions),
          std::move(newVertSnapR), std::move(crossings)};
}

Step9Threading ResolveAndThreadCrossings(
    const Manifold::Impl& impl, std::vector<NewEdgeWithExtras> chords,
    std::vector<vec3> newVertPositions, std::vector<double> newVertSnapR,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<OnChordContact>& contacts, double tolerance, double eps) {
  // Singleton-cluster delegation: each raw crossing is its own
  // cluster (the test seam; MergeAndPropagateCrossings supplies
  // real clusters in the full pipeline).
  std::vector<ChordCrossing> clusters;
  clusters.reserve(raw.size());
  for (const ChordChordCrossing& rc : raw) {
    clusters.push_back(
        {rc.pos, -1, {rc.chordA, rc.chordB}, {rc.tA, rc.tB}, rc.snapR});
  }
  return ResolveAndThreadClusters(
      impl, std::move(chords), std::move(newVertPositions),
      std::move(newVertSnapR), clusters, contacts, tolerance, eps);
}

std::vector<ChordCrossing> MergeAndPropagateCrossings(
    const Manifold::Impl& impl, const std::vector<NewEdgeWithExtras>& chords,
    const std::vector<vec3>& newVertPositions,
    const std::vector<ChordChordCrossing>& raw,
    const std::vector<std::vector<int>>& face2Chords,
    VecView<const vec3> faceNormals, double tolerance, double eps) {
  using la::dot;
  std::vector<ChordCrossing> out;
  if (raw.empty()) return out;
  const int baseId = static_cast<int>(impl.NumVert());
  const double snap = tolerance + eps;
  // New-to-new merge radius, matching boolean2's
  // kIntersectionMergeEpsFactor (covers shallow crossings to ~6 deg).
  const double mergeR = 10.0 * eps;
  const double mergeR2 = mergeR * mergeR;

  // Incident faces of a raw crossing: the hosting face plus both
  // chords' face pairs.
  const uint32_t n = static_cast<uint32_t>(raw.size());
  std::vector<std::set<int>> rcFaces(n);
  for (uint32_t i = 0; i < n; ++i) {
    rcFaces[i].insert(raw[i].face);
    for (const int ch : {raw[i].chordA, raw[i].chordB}) {
      rcFaces[i].insert(chords[ch].edge.triA);
      rcFaces[i].insert(chords[ch].edge.triB);
    }
  }

  // Face-gated union-find in sorted pair order: unite when the two
  // crossings share an incident face AND lie within 10 * eps. The
  // FACE gate (not a chord gate) is what unites a 4-chord concurrence
  // whose two crossings share no chord.
  DisjointSets uf(n);
  for (uint32_t i = 0; i < n; ++i) {
    for (uint32_t j = i + 1; j < n; ++j) {
      const vec3 d = raw[i].pos - raw[j].pos;
      if (dot(d, d) > mergeR2) continue;
      bool shareFace = false;
      for (const int f : rcFaces[i]) {
        if (rcFaces[j].count(f) > 0) {
          shareFace = true;
          break;
        }
      }
      if (!shareFace) continue;
      uf.unite(i, j);
    }
  }

  // Clusters keyed by their smallest member index, members ascending,
  // so output order and centroid summation are deterministic.
  std::map<uint32_t, uint32_t> minOfRoot;
  for (uint32_t i = 0; i < n; ++i) {
    const uint32_t r = uf.find(i);
    auto it = minOfRoot.find(r);
    if (it == minOfRoot.end()) minOfRoot[r] = i;
  }
  std::map<uint32_t, std::vector<uint32_t>> byMin;
  for (uint32_t i = 0; i < n; ++i) {
    byMin[minOfRoot[uf.find(i)]].push_back(i);
  }

  for (const auto& [minIdx, members] : byMin) {
    // Centroid in ascending member order, then re-projected onto the
    // hosting face plane (lowest hosting face) so the merged position
    // does not drift out of plane.
    vec3 centroid(0.0, 0.0, 0.0);
    int host = raw[members[0]].face;
    for (const uint32_t m : members) {
      centroid = centroid + raw[m].pos;
      host = std::min(host, raw[m].face);
    }
    centroid = centroid / static_cast<double>(members.size());
    if (host >= 0 && static_cast<size_t>(host) < faceNormals.size()) {
      const vec3 nRaw = faceNormals[host];
      const double nLen2 = dot(nRaw, nRaw);
      if (nLen2 > 0) {
        const vec3 nrm = nRaw / std::sqrt(nLen2);
        // The plane point must lie on the HOST face's plane: take it
        // from the first member hosted there. raw need not arrive in
        // face order, so members[0]'s hosting face can differ from
        // host - pairing host's normal with another face's point would
        // project the centroid onto neither plane.
        uint32_t hostMember = members[0];
        for (const uint32_t m : members) {
          if (raw[m].face == host) {
            hostMember = m;
            break;
          }
        }
        const vec3 planePt = GetPos3(chords[raw[hostMember].chordA].edge.v0,
                                     baseId, impl, newVertPositions);
        centroid = centroid - dot(centroid - planePt, nrm) * nrm;
      }
    }

    // Incident chords from the members' producing pairs...
    ChordCrossing cluster{centroid, -1, {}, {}, 0.0};
    for (const uint32_t m : members) {
      cluster.snapR = std::max(cluster.snapR, raw[m].snapR);
    }
    std::set<int> seenChords;
    auto addChord = [&](int ch, double t) {
      if (!seenChords.insert(ch).second) return;
      cluster.chords.push_back(ch);
      cluster.ts.push_back(t);
    };
    for (const uint32_t m : members) {
      addChord(raw[m].chordA, raw[m].tA);
      addChord(raw[m].chordB, raw[m].tB);
    }
    // ...plus eager propagation: every chord incident to any involved
    // face that the merged position lies on (point-to-segment <= eps,
    // the pass-0 endpoint-zone t-guard re-applied), so a k-fold point
    // lands on all k chords even when a pairwise crossing was missed.
    std::set<int> faces;
    for (const uint32_t m : members) {
      faces.insert(rcFaces[m].begin(), rcFaces[m].end());
    }
    const double eps2 = eps * eps;
    for (const int f : faces) {
      if (f < 0 || static_cast<size_t>(f) >= face2Chords.size()) continue;
      for (const int ch : face2Chords[f]) {
        if (seenChords.count(ch) > 0) continue;
        const PiercedNewEdge& e = chords[ch].edge;
        const vec3 a = GetPos3(e.v0, baseId, impl, newVertPositions);
        const vec3 b = GetPos3(e.v1, baseId, impl, newVertPositions);
        const vec3 ab = b - a;
        const double abLen2 = dot(ab, ab);
        if (abLen2 == 0) continue;
        const double tGuard = snap / std::sqrt(abLen2);
        const LineProj pr = ProjectToLine(centroid, a, ab, abLen2);
        if (pr.t <= tGuard || pr.t >= 1.0 - tGuard) continue;
        if (pr.distSq > eps2) continue;
        addChord(ch, pr.t);
      }
    }
    out.push_back(std::move(cluster));
  }
  return out;
}

std::vector<int> BuildHalfedgeToEdgeIndex(const Manifold::Impl& impl,
                                          const std::vector<Edge>& edges) {
  std::vector<int> out(impl.halfedge_.size(), -1);
  for (size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].halfedgeForward >= 0) {
      out[edges[i].halfedgeForward] = static_cast<int>(i);
    }
    if (edges[i].halfedgePaired >= 0) {
      out[edges[i].halfedgePaired] = static_cast<int>(i);
    }
  }
  return out;
}

FacePartition PartitionFace(const Manifold::Impl& impl, int face,
                            const std::vector<Edge>& edges,
                            const std::vector<int>& halfedge2Edge,
                            const std::vector<EdgeVertList>& onEdgeLists,
                            const std::vector<NewEdgeWithExtras>& chords,
                            const std::vector<int>& faceChords,
                            const std::vector<vec3>& newVertPositions,
                            double eps, bool coplanarHazard) {
  using la::cross;
  using la::dot;
  FacePartition out;
  const int baseId = static_cast<int>(impl.NumVert());
  DEBUG_ASSERT(face >= 0 && static_cast<size_t>(face) < impl.NumTri(), logicErr,
               "PartitionFace: face out of range");
  if (face < 0 || static_cast<size_t>(face) >= impl.NumTri()) return out;
  // The walk frame comes from the face's OWN halfedge winding, not
  // the stored faceNormal_: on self-intersecting inputs (folded
  // sheets) the stored normal can be OPPOSITE the winding, which
  // mirrors the projection and turns the face-on-left walk into a
  // boundary-hugging face-on-right walk (observed on the hull
  // fixture). CCW about this normal IS the halfedge order, by
  // construction.
  const vec3 fp0 = impl.vertPos_[impl.halfedge_.Start(3 * face)];
  const vec3 fp1 = impl.vertPos_[impl.halfedge_.Start(3 * face + 1)];
  const vec3 fp2 = impl.vertPos_[impl.halfedge_.Start(3 * face + 2)];
  const vec3 nRaw = cross(fp1 - fp0, fp2 - fp0);
  const double nLen2 = dot(nRaw, nRaw);
  if (nLen2 == 0) return out;
  const vec3 n = nRaw / std::sqrt(nLen2);
  const InPlaneBasis basis = FaceBasisFromNormal(n);
  const vec3 planePt = impl.vertPos_[impl.halfedge_.Start(3 * face)];
  auto to2d = [&](int id) -> vec2 {
    const vec3 p = GetPos3(id, baseId, impl, newVertPositions);
    const vec3 inPlane = p - dot(p - planePt, n) * n;
    return vec2(dot(inPlane - planePt, basis.u),
                dot(inPlane - planePt, basis.v));
  };

  // Directed sub-edges of the face graph. Original edges contribute
  // one halfedge per sub-edge along the face's CCW winding; chords
  // contribute BOTH directions, deduped per face by undirected vert
  // pair (coincident chords otherwise create exact angular ties) AND
  // against the face's own boundary sub-edges: a trace chord riding a
  // subdivided boundary (step 6.5 cuts the PARTNER face; on its host
  // the segment coincides with the boundary) would double a directed
  // edge, and the walk's exact-tie handling of the doubles is
  // hes-order-sensitive - skipping riders keeps the walk total and
  // order-independent.
  struct SubHalfedge {
    int start, end;
  };
  std::vector<SubHalfedge> hes;
  std::set<std::pair<int, int>> boundarySub;
  for (int k = 0; k < 3; ++k) {
    const Halfedge he = impl.halfedge_.Get(3 * face + k);
    const int a = he.startVert;
    const int b = he.endVert;
    const int ei = halfedge2Edge[3 * face + k];
    DEBUG_ASSERT(ei >= 0, logicErr,
                 "PartitionFace: face edge missing from EnumerateEdges");
    if (ei < 0) continue;
    const EdgeVertList& evl = onEdgeLists[ei];
    std::vector<int> seq;
    seq.push_back(a);
    if (a == edges[ei].v0) {
      seq.insert(seq.end(), evl.verts.begin(), evl.verts.end());
    } else {
      seq.insert(seq.end(), evl.verts.rbegin(), evl.verts.rend());
    }
    seq.push_back(b);
    for (size_t i = 0; i + 1 < seq.size(); ++i) {
      if (seq[i] == seq[i + 1]) continue;  // snapped duplicates
      hes.push_back({seq[i], seq[i + 1]});
      boundarySub.insert(
          {std::min(seq[i], seq[i + 1]), std::max(seq[i], seq[i + 1])});
    }
  }
  std::set<std::pair<int, int>> seenSub;
  for (const int ci : faceChords) {
    const NewEdgeWithExtras& nwe = chords[ci];
    if (nwe.edge.v0 == nwe.edge.v1) {
      ++out.zeroLengthChordsSkipped;
      continue;
    }
    std::vector<int> seq;
    seq.push_back(nwe.edge.v0);
    seq.insert(seq.end(), nwe.extraVerts.begin(), nwe.extraVerts.end());
    seq.push_back(nwe.edge.v1);
    for (size_t i = 0; i + 1 < seq.size(); ++i) {
      const int a = seq[i];
      const int b = seq[i + 1];
      if (a == b) continue;
      const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
      if (boundarySub.count(key)) {
        ++out.boundaryRidingSubEdgesSkipped;
        continue;
      }
      if (!seenSub.insert(key).second) continue;
      hes.push_back({a, b});
      hes.push_back({b, a});
    }
  }
  if (hes.empty()) return out;

  // `skippedRiderVerts`: verts from chord sub-edges that were skipped
  // because they coincide with the face boundary (boundary riders).
  // A free-island vert in this set has a hidden boundary attachment - it
  // must gate.
  std::set<int> skippedRiderVerts;
  for (const int ci : faceChords) {
    const NewEdgeWithExtras& nwe = chords[ci];
    if (nwe.edge.v0 == nwe.edge.v1) continue;
    std::vector<int> riderSeq;
    riderSeq.push_back(nwe.edge.v0);
    riderSeq.insert(riderSeq.end(), nwe.extraVerts.begin(),
                    nwe.extraVerts.end());
    riderSeq.push_back(nwe.edge.v1);
    for (size_t i = 0; i + 1 < riderSeq.size(); ++i) {
      const int a = riderSeq[i];
      const int b = riderSeq[i + 1];
      if (a == b) continue;
      const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
      if (boundarySub.count(key) && !seenSub.count(key)) {
        skippedRiderVerts.insert(a);
        skippedRiderVerts.insert(b);
      }
    }
  }

  // INTERIOR-ISLAND GATE (reclassified): chord-only components containing
  // a cycle (compEdges >= compVerts) with fewer than 2 boundary attachments.
  //
  //  compAttach == 0: free island - hole CANDIDATE. Ungate only when:
  //    (a) every vert has degree exactly 2 in the component (one loop);
  //    (b) no component vert is in skippedRiderVerts (hidden attachment);
  //    (c) no component vert is shared with any other candidate island;
  //    (d) loop is SIMPLE (checked via 2D projections after to2d is set);
  //    (e) loop has nonzero projected area;
  //    (f) this face is NOT coplanar-hazard-flagged.
  //    Failure of any condition -> gate as today.
  //  compAttach == 1: PINCHED - gate as today.
  //  compAttach >= 2 or tree: pass as before.
  //
  // Collect clean-island loop sequences here (for decomposition later).
  // Each entry is the undirected edge set of one clean island loop.
  struct IslandLoop {
    std::vector<int> verts;  // sequence from the degree-2 walk
    std::set<int> vertSet;   // for fast containment checks
  };
  std::vector<IslandLoop> cleanIslandLoops;

  {
    std::map<int, int> vert2Idx;
    auto idxOf = [&](int v) {
      const auto [it, fresh] =
          vert2Idx.insert({v, static_cast<int>(vert2Idx.size())});
      return it->second;
    };
    for (const auto& [a, b] : seenSub) {
      idxOf(a);
      idxOf(b);
    }
    if (!vert2Idx.empty()) {
      std::set<int> boundaryVerts;
      for (const auto& [a, b] : boundarySub) {
        boundaryVerts.insert(a);
        boundaryVerts.insert(b);
      }
      DisjointSets uf(static_cast<uint32_t>(vert2Idx.size()));
      for (const auto& [a, b] : seenSub) {
        uf.unite(static_cast<uint32_t>(vert2Idx[a]),
                 static_cast<uint32_t>(vert2Idx[b]));
      }
      std::map<uint32_t, int> compVerts, compEdges, compAttach;
      std::map<uint32_t, std::vector<int>> compVertList;
      // Per-component, per-vert degree within the component.
      std::map<uint32_t, std::map<int, std::vector<int>>> compAdj;
      for (const auto& [v, idx] : vert2Idx) {
        const uint32_t r = uf.find(static_cast<uint32_t>(idx));
        ++compVerts[r];
        compVertList[r].push_back(v);
        if (boundaryVerts.count(v)) ++compAttach[r];
      }
      for (const auto& [a, b] : seenSub) {
        const uint32_t r = uf.find(static_cast<uint32_t>(vert2Idx[a]));
        ++compEdges[r];
        compAdj[r][a].push_back(b);
        compAdj[r][b].push_back(a);
      }

      // Collect candidate island roots (compAttach == 0, has cycle,
      // degree-2, no rider verts). Gate everything else with cycles.
      // Must also check for cross-component vert sharing among candidates.
      std::set<uint32_t> candidateRoots;
      for (const auto& [r, nV] : compVerts) {
        if (compEdges[r] < nV) continue;   // tree: no cycle, passes
        if (compAttach[r] >= 2) continue;  // properly attached: passes
        if (compAttach[r] == 1) {
          out.interiorIslandVerts += nV;  // pinched: gate
          continue;
        }
        // Free island (compAttach == 0). Check degree-2.
        bool allDeg2 = true;
        for (const int v : compVertList[r]) {
          if (static_cast<int>(compAdj[r][v].size()) != 2) {
            allDeg2 = false;
            break;
          }
        }
        if (!allDeg2) {
          out.interiorIslandVerts += nV;
          continue;
        }
        // Check no rider verts.
        bool hasRider = false;
        for (const int v : compVertList[r]) {
          if (skippedRiderVerts.count(v)) {
            hasRider = true;
            break;
          }
        }
        if (hasRider) {
          out.interiorIslandVerts += nV;
          continue;
        }
        candidateRoots.insert(r);
      }

      // Cross-component vert sharing among candidates: gate all.
      if (candidateRoots.size() > 1) {
        std::map<int, int> vertOwner;
        bool hasConflict = false;
        for (const uint32_t r : candidateRoots) {
          for (const int v : compVertList[r]) {
            if (!vertOwner.insert({v, static_cast<int>(r)}).second) {
              hasConflict = true;
              break;
            }
          }
          if (hasConflict) break;
        }
        if (hasConflict) {
          for (const uint32_t r : candidateRoots) {
            out.interiorIslandVerts += compVerts[r];
          }
          candidateRoots.clear();
        }
      }

      // Build the loop sequence for each surviving candidate (degree-2
      // walk: pick any start vert, follow the unique next neighbor not
      // yet visited).
      for (const uint32_t r : candidateRoots) {
        const std::vector<int>& vlist = compVertList[r];
        if (vlist.empty()) continue;
        // Walk the degree-2 graph to get an ordered loop sequence.
        std::set<int> loopVisited;
        std::vector<int> loopSeq;
        int cur = vlist[0];
        int prev = -1;
        for (;;) {
          if (!loopVisited.insert(cur).second) break;
          loopSeq.push_back(cur);
          int nxt = -1;
          for (const int nb : compAdj[r][cur]) {
            if (nb != prev) {
              nxt = nb;
              break;
            }
          }
          if (nxt < 0 || loopVisited.count(nxt)) break;
          prev = cur;
          cur = nxt;
        }
        if (loopSeq.size() < 3) {
          // Degenerate: gate.
          out.interiorIslandVerts += compVerts[r];
          continue;
        }
        IslandLoop isl;
        isl.verts = std::move(loopSeq);
        for (const int v : isl.verts) isl.vertSet.insert(v);
        cleanIslandLoops.push_back(std::move(isl));
      }

      // If any components were gated AND there are candidates, gate all
      // candidates too (conservative: mixed-state face is ambiguous).
      if (out.interiorIslandVerts > 0 && !cleanIslandLoops.empty()) {
        for (const IslandLoop& isl : cleanIslandLoops) {
          out.interiorIslandVerts += static_cast<int>(isl.verts.size());
        }
        cleanIslandLoops.clear();
        return out;
      }
      if (out.interiorIslandVerts > 0) return out;
    }
  }

  // Outgoing lists per vert, ordered CCW about the face normal by the
  // atan2-free comparator (half-plane bucket + cross sign - the
  // boolean2 winding_filter pattern; no trig in the decision path).
  std::map<int, vec2> pos2;
  auto p2 = [&](int id) -> vec2 {
    auto it = pos2.find(id);
    if (it == pos2.end()) it = pos2.insert({id, to2d(id)}).first;
    return it->second;
  };

  // ISLAND GEOMETRY CHECKS: for each candidate clean island loop,
  // verify simplicity and nonzero area. These need to2d (which uses the
  // already-computed n, basis, planePt above) via p2. Done here so the
  // same lazy pos2 cache is shared with the walk.
  //
  // The `triEps` ladder seed for TriangulateIdx.
  const double triEps = std::max({impl.tolerance_, impl.epsilon_, eps});
  if (!cleanIslandLoops.empty()) {
    // Geometry check helper: is loop simple in 2D?
    // Full PSLG rules: no strict crossings, no endpoint-on-nonincident-edge.
    // Collinear overlap and duplicate edges cannot arise from degree-2 loops
    // with distinct verts (checked below via the simple-segment pairwise
    // test). Returns false if any violation found.
    auto isLoopSimple2D = [&](const std::vector<int>& loopV) -> bool {
      const int nL = static_cast<int>(loopV.size());
      if (nL < 3) return false;
      // Collect 2D positions.
      std::vector<vec2> pts(nL);
      for (int i = 0; i < nL; ++i) pts[i] = p2(loopV[i]);
      // Pairwise non-adjacent segment test.
      // Segments: i -> (i+1)%nL. Non-adjacent means they share no vert.
      // Rejects ALL non-incident contact: strict crossings,
      // endpoint-on-interior of the other segment, endpoint-on-endpoint between
      // non-adjacent segments, and collinear overlap. Zero-determinant cases
      // (degenerate/touching) are NOT skipped - any contact between
      // non-incident segments is a violation.
      auto cross2 = [](const vec2& a, const vec2& b) {
        return a.x * b.y - a.y * b.x;
      };
      // pointOnSeg: does P lie on segment AB (including endpoints)?
      // Uses exact arithmetic (no tolerance) since these are 2D projected
      // positions that will be exact at this level.
      auto pointOnSeg = [&](const vec2& P, const vec2& A, const vec2& B) {
        const vec2 ab = {B.x - A.x, B.y - A.y};
        const vec2 ap = {P.x - A.x, P.y - A.y};
        // Must be collinear.
        if (cross2(ab, ap) != 0.0) return false;
        // Must be in [0,1] parametrically: dot in [0, len2].
        const double dot = ab.x * ap.x + ab.y * ap.y;
        const double len2 = ab.x * ab.x + ab.y * ab.y;
        return dot >= 0.0 && dot <= len2;
      };
      for (int i = 0; i < nL; ++i) {
        const vec2& a = pts[i];
        const vec2& b = pts[(i + 1) % nL];
        for (int j = i + 1; j < nL; ++j) {
          const vec2& c = pts[j];
          const vec2& d = pts[(j + 1) % nL];
          // Determine adjacency: adjacent pairs share exactly one vertex
          // (consecutive segments in the loop). The wrap-around pair
          // (i==0, j==nL-1) shares vertex a==d (index 0).
          const bool isNext = (j == i + 1);             // shares b == c
          const bool isWrap = (i == 0 && j == nL - 1);  // shares a == d
          const bool adjacent = isNext || isWrap;
          if (!adjacent) {
            // Non-adjacent: any contact is a simplicity violation.
            if (pointOnSeg(a, c, d)) return false;  // endpoint a on cd
            if (pointOnSeg(b, c, d)) return false;  // endpoint b on cd
            if (pointOnSeg(c, a, b)) return false;  // endpoint c on ab
            if (pointOnSeg(d, a, b)) return false;  // endpoint d on ab
            // Strict crossing (both pairs of opposite signs).
            const vec2 ab = {b.x - a.x, b.y - a.y};
            const vec2 ac = {c.x - a.x, c.y - a.y};
            const vec2 ad_v = {d.x - a.x, d.y - a.y};
            const double t1 = cross2(ab, ac);
            const double t2 = cross2(ab, ad_v);
            if (t1 * t2 < 0) {
              const vec2 cd = {d.x - c.x, d.y - c.y};
              const vec2 ca = {a.x - c.x, a.y - c.y};
              const vec2 cb = {b.x - c.x, b.y - c.y};
              const double t3 = cross2(cd, ca);
              const double t4 = cross2(cd, cb);
              if (t3 * t4 < 0) return false;  // strict crossing
            }
          } else {
            // Adjacent: one vertex is legitimately shared. Reject only if
            // the NON-shared endpoint of either segment lies on the other
            // segment (point-on-seg incl. endpoints), or if the segments
            // overlap collinearly with positive length beyond the shared
            // point (which the non-shared endpoint checks subsume).
            if (isNext) {
              // Shared vertex: b == c. Non-shared: a (of ab) and d (of cd).
              if (pointOnSeg(a, c, d)) return false;  // a on cd (beyond b)
              if (pointOnSeg(d, a, b)) return false;  // d on ab (beyond b)
            } else {
              // isWrap: shared vertex a == d (index 0). Non-shared: b and c.
              if (pointOnSeg(b, c, d)) return false;  // b on cd (beyond a)
              if (pointOnSeg(c, a, b)) return false;  // c on ab (beyond a)
            }
          }
        }
      }
      return true;
    };

    // Compute signed area (doubled) of a loop.
    auto signedArea2 = [&](const std::vector<int>& loopV) -> double {
      double a = 0.0;
      const int nL = static_cast<int>(loopV.size());
      for (int i = 0; i < nL; ++i) {
        const vec2 pi = p2(loopV[i]);
        const vec2 pj = p2(loopV[(i + 1) % nL]);
        a += pi.x * pj.y - pi.y * pj.x;
      }
      return a;
    };

    // Filter candidate loops through geometry checks.
    std::vector<IslandLoop> passingLoops;
    for (const IslandLoop& isl : cleanIslandLoops) {
      const double sa2 = signedArea2(isl.verts);
      if (sa2 == 0.0) {
        // Zero area: gate.
        out.interiorIslandVerts += static_cast<int>(isl.verts.size());
        continue;
      }
      if (!isLoopSimple2D(isl.verts)) {
        // Not simple: gate.
        out.interiorIslandVerts += static_cast<int>(isl.verts.size());
        continue;
      }
      passingLoops.push_back(isl);
    }
    cleanIslandLoops = std::move(passingLoops);

    // If any loop failed geometry or the face is coplanar-hazard:
    // gate all remaining.
    if (out.interiorIslandVerts > 0 || coplanarHazard) {
      for (const IslandLoop& isl : cleanIslandLoops) {
        out.interiorIslandVerts += static_cast<int>(isl.verts.size());
      }
      cleanIslandLoops.clear();
      return out;
    }
  }

  auto dirOf = [&](int he) -> vec2 {
    return p2(hes[he].end) - p2(hes[he].start);
  };
  auto bucketOf = [](const vec2& d) {
    return (d.y > 0 || (d.y == 0 && d.x > 0)) ? 0 : 1;
  };
  auto angleStrictLess = [&](const vec2& dA, const vec2& dB) {
    const int bA = bucketOf(dA);
    const int bB = bucketOf(dB);
    if (bA != bB) return bA < bB;
    return dA.x * dB.y - dA.y * dB.x > 0;
  };
  std::map<int, std::vector<int>> outgoing;
  for (size_t i = 0; i < hes.size(); ++i) {
    outgoing[hes[i].start].push_back(static_cast<int>(i));
  }
  for (auto& [vtx, list] : outgoing) {
    std::sort(list.begin(), list.end(), [&](int x, int y) {
      const vec2 dx = dirOf(x);
      const vec2 dy = dirOf(y);
      if (angleStrictLess(dx, dy)) return true;
      if (angleStrictLess(dy, dx)) return false;
      return x < y;  // deterministic exact-tie fallback
    });
  }

  // Successor of h (a -> b): among UNVISITED outgoing halfedges at b,
  // excluding the immediate reverse (b -> a) UNLESS it is the sole
  // candidate (the U-turn that traverses dangling-chord spurs), the
  // entry whose direction is the cyclic predecessor of the reverse
  // direction in CCW order - the smallest left turn.
  std::vector<bool> visited(hes.size(), false);
  auto successor = [&](int h) -> int {
    const int b = hes[h].end;
    const auto it = outgoing.find(b);
    if (it == outgoing.end()) return -1;
    const vec2 q = p2(hes[h].start) - p2(b);  // reverse direction
    int reverseHe = -1;
    int best = -1;      // max angle among directions strictly below q
    int bestWrap = -1;  // max angle overall (cyclic wrap)
    for (const int e : it->second) {
      if (visited[e]) continue;
      if (hes[e].end == hes[h].start) {
        // At most one reverse exists: seenSub dedups undirected chord
        // pairs and the boundary contributes one direction per
        // sub-edge - keep that invariant locally visible.
        DEBUG_ASSERT(reverseHe < 0, logicErr,
                     "PartitionFace: duplicate reverse halfedge");
        reverseHe = e;
        continue;
      }
      const vec2 d = dirOf(e);
      if (bestWrap < 0 || angleStrictLess(dirOf(bestWrap), d)) bestWrap = e;
      if (angleStrictLess(d, q) &&
          (best < 0 || angleStrictLess(dirOf(best), d))) {
        best = e;
      }
    }
    if (best >= 0) return best;
    if (bestWrap >= 0) return bestWrap;
    return reverseHe;  // sole candidate -> U-turn; -1 if none at all
  };

  // Walk every halfedge once; closure is VERTEX ARRIVAL (destination
  // equals the walk's start vert - the boolean2 OutEdgesToPolygons
  // pattern). Then split each closed cycle at repeated vert ids
  // (PushSimpleLoops pattern); sub-3-vert loops are spurs - dropped
  // and counted.
  for (size_t h0 = 0; h0 < hes.size(); ++h0) {
    if (visited[h0]) continue;
    const int startV = hes[h0].start;
    std::vector<int> loop;
    int cur = static_cast<int>(h0);
    bool closed = false;
    const size_t maxSteps = hes.size() + 1;
    for (size_t step = 0; step < maxSteps && cur >= 0; ++step) {
      visited[cur] = true;
      loop.push_back(hes[cur].start);
      if (hes[cur].end == startV) {
        closed = true;
        break;
      }
      cur = successor(cur);
    }
    DEBUG_ASSERT(closed, logicErr,
                 "PartitionFace: open walk in a conforming arrangement");
    if (!closed) continue;
    // Emit gate: sub-3-vert loops are spurs; a >= 3-vert simple cycle
    // whose projected signed area is EXACTLY zero is a flattened spur
    // (the near-line sliver artifact: post-merge coincident positions
    // under distinct ids contribute exactly-cancelling cross terms,
    // and an out-and-back walk along a line bounds nothing) - its
    // Newell normal is undefined downstream, so it is dropped and
    // counted. Tiny-but-nonzero areas are REAL slivers and pass.
    auto emit = [&](std::vector<int>&& cyc) {
      if (cyc.size() < 3) {
        ++out.spursDropped;
        return;
      }
      double area2 = 0.0;
      for (size_t i = 0; i < cyc.size(); ++i) {
        const vec2 pa = p2(cyc[i]);
        const vec2 pb = p2(cyc[(i + 1) % cyc.size()]);
        area2 += pa.x * pb.y - pa.y * pb.x;
      }
      if (area2 == 0.0) {
        ++out.degenerateCyclesDropped;
        return;
      }
      out.polygons.push_back(std::move(cyc));
    };
    for (;;) {
      bool split = false;
      for (size_t i = 1; i < loop.size() && !split; ++i) {
        for (size_t j = 0; j < i; ++j) {
          if (loop[i] != loop[j]) continue;
          std::vector<int> sub(loop.begin() + j, loop.begin() + i);
          emit(std::move(sub));
          loop.erase(loop.begin() + j + 1, loop.begin() + i + 1);
          split = true;
          break;
        }
      }
      if (!split) break;
    }
    emit(std::move(loop));
  }

  // POST-WALK ISLAND DECOMPOSITION: only when clean island loops were
  // detected. Split each region-with-holes (positive-area outer cycle +
  // one or more negative-area island cycles) into triangles via
  // TriangulateIdx. Regions without holes pass through untouched.
  //
  // All polygons emitted by the walk are in face-winding orientation
  // (CCW about the face normal). Island loops were inserted as chord
  // sub-edges in BOTH directions so the walk emits them twice: one
  // positive-area cycle (the disk inside the stamp footprint) and one
  // negative-area cycle (the hole in the outer region). The positive
  // orientation is the island loop's disk polygon and passes through
  // as an ordinary region. The negative orientation is the hole.
  if (!cleanIslandLoops.empty()) {
    // Identify island vert sets for twin matching.
    // A cycle whose vert set equals any clean island loop's vert set is
    // that island's orientation twin (the one we treat as the hole).
    // Match by vert set equality (canonical: sorted).
    auto cycleVertKey = [](const std::vector<int>& cyc) {
      std::vector<int> s = cyc;
      std::sort(s.begin(), s.end());
      return s;
    };
    std::set<std::vector<int>> islandVertKeys;
    for (const IslandLoop& isl : cleanIslandLoops) {
      std::vector<int> key = isl.verts;
      std::sort(key.begin(), key.end());
      islandVertKeys.insert(std::move(key));
    }

    // Partition walked polygons into positive (regions) and negative
    // (holes). Compute signed area for each.
    struct CycleInfo {
      std::vector<int> cyc;
      double signedArea2;  // doubled signed area (positive = CCW = region)
      std::vector<int> vertKey;  // sorted, for island twin matching
    };
    std::vector<CycleInfo> regions, holes;
    for (const std::vector<int>& cyc : out.polygons) {
      double sa2 = 0.0;
      for (size_t i = 0; i < cyc.size(); ++i) {
        const vec2 pi = p2(cyc[i]);
        const vec2 pj = p2(cyc[(i + 1) % cyc.size()]);
        sa2 += pi.x * pj.y - pi.y * pj.x;
      }
      std::vector<int> key = cycleVertKey(cyc);
      if (sa2 >= 0) {
        regions.push_back({cyc, sa2, std::move(key)});
      } else {
        holes.push_back({cyc, sa2, std::move(key)});
      }
    }

    // For each hole, find its assigned region: the smallest-area
    // strictly-containing positive cycle (excluding its orientation twin
    // and any cycle sharing a vert with it; border -> not contained;
    // tie -> walk order i.e. index in regions).
    struct HoleAssignment {
      size_t regionIdx;
      size_t holeIdx;
    };
    std::vector<HoleAssignment> assignments;

    // Build hole-vert sets for sharing check.
    std::vector<std::set<int>> holeVertSets(holes.size());
    for (size_t hi = 0; hi < holes.size(); ++hi) {
      for (const int v : holes[hi].cyc) holeVertSets[hi].insert(v);
    }

    // Strict point-in-polygon test (2D): is test point P strictly inside
    // polygon Q? Returns false if P lies on any boundary edge (including
    // at vertices). Explicit point-on-segment rejection runs first, before
    // the ray-crossing count, so horizontal edges and vertex-at-exact-y
    // cases do not yield false positives.
    auto pointOnSeg2D = [](const vec2& P, const vec2& A, const vec2& B) {
      const vec2 ab = {B.x - A.x, B.y - A.y};
      const vec2 ap = {P.x - A.x, P.y - A.y};
      if (ab.x * ap.y - ab.y * ap.x != 0.0) return false;
      const double dot = ab.x * ap.x + ab.y * ap.y;
      const double len2 = ab.x * ab.x + ab.y * ab.y;
      return dot >= 0.0 && dot <= len2;
    };
    auto strictPIP = [&](const vec2& P, const std::vector<int>& Q) -> bool {
      const int nQ = static_cast<int>(Q.size());
      // First: reject if P lies on any boundary edge (incl. at vertices).
      for (int i = 0; i < nQ; ++i) {
        const vec2 a = p2(Q[i]);
        const vec2 b = p2(Q[(i + 1) % nQ]);
        if (pointOnSeg2D(P, a, b)) return false;
      }
      // Ray-casting count: ray from P in +x direction.
      int crossings = 0;
      for (int i = 0; i < nQ; ++i) {
        const vec2 a = p2(Q[i]);
        const vec2 b = p2(Q[(i + 1) % nQ]);
        if ((a.y <= P.y) == (b.y <= P.y)) continue;  // same side
        const double x = a.x + (P.y - a.y) * (b.x - a.x) / (b.y - a.y);
        if (x > P.x) ++crossings;
      }
      return (crossings & 1) == 1;
    };

    bool decompositionFailed = false;
    for (size_t hi = 0; hi < holes.size(); ++hi) {
      const std::vector<int>& holeCyc = holes[hi].cyc;
      const std::vector<int>& holeKey = holes[hi].vertKey;
      // The hole's vert set (from holeVertSets[hi]).

      // Representative vert for point-in-polygon test.
      const vec2 testPt = p2(holeCyc[0]);

      double bestArea2 = std::numeric_limits<double>::infinity();
      int bestRegion = -1;
      for (size_t ri = 0; ri < regions.size(); ++ri) {
        const std::vector<int>& regionKey = regions[ri].vertKey;
        // Exclude the orientation twin: same vert set.
        if (regionKey == holeKey) continue;
        // Exclude any cycle sharing a vert with this hole.
        bool shares = false;
        for (const int v : holeVertSets[hi]) {
          if (std::binary_search(regionKey.begin(), regionKey.end(), v)) {
            shares = true;
            break;
          }
        }
        if (shares) continue;
        // Strict point-in-polygon.
        if (!strictPIP(testPt, regions[ri].cyc)) continue;
        // Smallest area; tie by walk order (ri).
        if (regions[ri].signedArea2 < bestArea2 ||
            (regions[ri].signedArea2 == bestArea2 && bestRegion >= 0 &&
             ri < static_cast<size_t>(bestRegion))) {
          bestArea2 = regions[ri].signedArea2;
          bestRegion = static_cast<int>(ri);
        }
      }
      if (bestRegion < 0) {
        // No candidate: gate.
        decompositionFailed = true;
        break;
      }
      assignments.push_back({static_cast<size_t>(bestRegion), hi});
    }

    if (decompositionFailed) {
      out.interiorIslandVerts = 1;  // signal gate to driver
      out.polygons.clear();
      return out;
    }

    // Group holes by region.
    std::map<size_t, std::vector<size_t>> regionHoles;
    for (const HoleAssignment& ha : assignments) {
      regionHoles[ha.regionIdx].push_back(ha.holeIdx);
    }

    // Rebuild out.polygons with decompositions for holed regions.
    std::vector<std::vector<int>> newPolygons;

    // Process ALL regions: any region with no assigned holes passes through
    // as a plain polygon (including island-disk cycles). An island-disk WITH
    // assigned holes (nested island: B inside A's disk) must be triangulated
    // to cut the nested hole out of the disk - passing it through uncut would
    // leave the nested island's negative cycle unconsumed.
    for (size_t ri = 0; ri < regions.size(); ++ri) {
      const CycleInfo& reg = regions[ri];

      auto itH = regionHoles.find(ri);
      if (itH == regionHoles.end()) {
        // No holes: pass through (island-disk or ordinary region).
        newPolygons.push_back(reg.cyc);
        continue;
      }

      // Holed region: triangulate.
      // Sub-resolution hole fast-fail: doubled area < triEps * maxBbox.
      // Compute 2D bbox of the region.
      // Build PolygonsIdx: outer (positive-area) first, then holes
      // (negative-area - the triangulator's convention).
      PolygonsIdx polysIdx;
      // Outer: winding as walked (CCW = positive).
      {
        SimplePolygonIdx outer;
        outer.reserve(reg.cyc.size());
        for (const int v : reg.cyc) {
          outer.push_back({p2(v), v});
        }
        polysIdx.push_back(std::move(outer));
      }
      // Compute the outer bbox once for the area-preservation tolerance.
      double rMinU = std::numeric_limits<double>::infinity();
      double rMaxU = -std::numeric_limits<double>::infinity();
      double rMinV = std::numeric_limits<double>::infinity();
      double rMaxV = -std::numeric_limits<double>::infinity();
      for (const int v : reg.cyc) {
        const vec2 pv = p2(v);
        rMinU = std::min(rMinU, pv.x);
        rMaxU = std::max(rMaxU, pv.x);
        rMinV = std::min(rMinV, pv.y);
        rMaxV = std::max(rMaxV, pv.y);
      }
      const double bboxMax = std::max(rMaxU - rMinU, rMaxV - rMinV);
      bool subResGate = false;
      for (const size_t hi : itH->second) {
        const CycleInfo& hci = holes[hi];
        // Sub-resolution fast-fail: conservative check using the hole
        // contour's own 2D bbox (FindStart uses the contour's own extent,
        // so per-hole bbox is the correct scale for this fast-fail).
        double hMinU = std::numeric_limits<double>::infinity();
        double hMaxU = -std::numeric_limits<double>::infinity();
        double hMinV = std::numeric_limits<double>::infinity();
        double hMaxV = -std::numeric_limits<double>::infinity();
        for (const int v : hci.cyc) {
          const vec2 pv = p2(v);
          hMinU = std::min(hMinU, pv.x);
          hMaxU = std::max(hMaxU, pv.x);
          hMinV = std::min(hMinV, pv.y);
          hMaxV = std::max(hMaxV, pv.y);
        }
        const double holeBboxMax = std::max(hMaxU - hMinU, hMaxV - hMinV);
        if (std::fabs(hci.signedArea2) < triEps * holeBboxMax) {
          subResGate = true;
          break;
        }
        // Hole cycle: appended AS WALKED - the negative-area
        // orientation is what the walk already produced for the
        // island's hole twin; reversing here would double-flip it
        // positive and force the validation fallback.
        SimplePolygonIdx hole;
        hole.reserve(hci.cyc.size());
        for (const int v : hci.cyc) {
          hole.push_back({p2(v), v});
        }
        polysIdx.push_back(std::move(hole));
      }
      if (subResGate) {
        out.interiorIslandVerts = 1;
        out.polygons.clear();
        return out;
      }

      // TriangulateIdx call with local debug try/catch.
      std::vector<ivec3> tris;
      bool triOk = true;
#ifdef MANIFOLD_DEBUG
      try {
        tris = TriangulateIdx(polysIdx, triEps);
      } catch (...) {
        triOk = false;
      }
#else
      tris = TriangulateIdx(polysIdx, triEps);
#endif
      if (!triOk || tris.empty()) {
        out.interiorIslandVerts = 1;
        out.polygons.clear();
        return out;
      }

      // VALIDATION (every build config, as required by the plan):
      // 1. Per-triangle: 3 distinct ids, strictly positive area.
      // 2. Signed-area preservation.
      // 3. Boundary coverage.
      // 4. Full PSLG embedding (no strict crossings, no
      //    vertex-on-nonincident-edge, no collinear overlaps, no
      //    duplicate edges) - quadratic but bounded by edge-count cap.

      // Build expected boundary edges from outer + holes.
      struct EdgeDir {
        int v0, v1;
      };
      std::vector<EdgeDir> boundaryEdges;
      for (const SimplePolygonIdx& cont : polysIdx) {
        const int nc = static_cast<int>(cont.size());
        for (int i = 0; i < nc; ++i) {
          boundaryEdges.push_back({cont[i].idx, cont[(i + 1) % nc].idx});
        }
      }
      const int totalBoundaryEdges = static_cast<int>(boundaryEdges.size());
      // Edge-count cap: regions exceeding this gate (fail closed) rather than
      // skipping validation. In practice the boss class is a triangle outer
      // + a handful of hole verts, so this cap is never hit on valid input.
      const int kMaxEdgesForValidation = 2000;
      const bool exceedsCap =
          totalBoundaryEdges + 3 * static_cast<int>(tris.size()) >
          kMaxEdgesForValidation;
      bool validFail = exceedsCap;  // gate if over cap

      // Build the set of canonical hole cycle keys for the 3-vert hole check:
      // a triangle whose vert set matches a hole contour would be filling in
      // the hole rather than cutting it, which is wrong.
      std::set<std::vector<int>> holeCycleKeys;
      for (size_t hIdx = 1; hIdx < polysIdx.size(); ++hIdx) {
        const SimplePolygonIdx& hcont = polysIdx[hIdx];
        if (hcont.size() == 3) {
          std::vector<int> hkey = {hcont[0].idx, hcont[1].idx, hcont[2].idx};
          std::sort(hkey.begin(), hkey.end());
          holeCycleKeys.insert(std::move(hkey));
        }
      }

      auto area2OfTri = [&](const ivec3& t) -> double {
        const vec2 pa = p2(t[0]);
        const vec2 pb = p2(t[1]);
        const vec2 pc = p2(t[2]);
        return (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);
      };

      // Per-triangle checks.
      if (!validFail) {
        for (const ivec3& t : tris) {
          if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) {
            validFail = true;
            break;
          }
          if (area2OfTri(t) <= 0.0) {
            validFail = true;
            break;
          }
          // Check not canonical-equal to any 3-vert hole contour.
          if (!holeCycleKeys.empty()) {
            std::vector<int> tkey = {t[0], t[1], t[2]};
            std::sort(tkey.begin(), tkey.end());
            if (holeCycleKeys.count(tkey)) {
              validFail = true;
              break;
            }
          }
        }
      }

      // Signed-area preservation:
      // expected = sum of signed areas of polysIdx contours (holes are
      // negative). Doubled areas (no 0.5 factor).
      if (!validFail) {
        double expectedArea2 = 0.0;
        for (const SimplePolygonIdx& cont : polysIdx) {
          const int nc = static_cast<int>(cont.size());
          for (int i = 0; i < nc; ++i) {
            const vec2 pi2 = cont[i].pos;
            const vec2 pj2 = cont[(i + 1) % nc].pos;
            expectedArea2 += pi2.x * pj2.y - pi2.y * pj2.x;
          }
        }
        double triArea2 = 0.0;
        for (const ivec3& t : tris) triArea2 += area2OfTri(t);
        const double kAreaRelTol = 1e-9;
        const double areaAbsTol =
            triEps * bboxMax + kAreaRelTol * std::fabs(expectedArea2);
        if (std::fabs(triArea2 - expectedArea2) > areaAbsTol) validFail = true;
      }

      // Boundary coverage (undirected semantics):
      // - Each boundary edge {v0,v1} appears exactly once in the required
      //   direction in the triangle fan.
      // - The reverse {v1,v0} of a boundary edge must NOT appear as a
      //   boundary-role edge (i.e., must not be in the boundary set itself).
      // - Each interior (non-boundary) undirected edge appears exactly twice
      //   in opposite directions: {u,v} with count 1 and {v,u} with count 1.
      // Uses .find()/.count() throughout - no operator[] phantom entries.
      if (!validFail) {
        std::map<std::pair<int, int>, int> edgeCount;
        for (const ivec3& t : tris) {
          for (int k = 0; k < 3; ++k) {
            const auto key = std::make_pair(t[k], t[(k + 1) % 3]);
            auto it = edgeCount.find(key);
            if (it == edgeCount.end())
              edgeCount.insert({key, 1});
            else
              ++it->second;
          }
        }
        // Build boundary set (required-direction directed edges).
        std::set<std::pair<int, int>> boundaryEdgeSet;
        for (const EdgeDir& be : boundaryEdges) {
          boundaryEdgeSet.insert({be.v0, be.v1});
        }
        // Check 1: each boundary edge appears exactly once in its direction.
        for (const EdgeDir& be : boundaryEdges) {
          auto it = edgeCount.find({be.v0, be.v1});
          if (it == edgeCount.end() || it->second != 1) {
            validFail = true;
            break;
          }
        }
        // Check 2: reverse of any boundary edge must not be a boundary edge
        // (directedness: {v1,v0} appearing in boundaryEdgeSet means a
        // boundary edge is present in the wrong orientation somewhere).
        if (!validFail) {
          for (const EdgeDir& be : boundaryEdges) {
            if (boundaryEdgeSet.count({be.v1, be.v0})) {
              validFail = true;
              break;
            }
          }
        }
        // Build undirected boundary edge set for reverse-occurrence check.
        // A boundary edge {v0,v1} appearing in the REQUIRED direction is
        // already checked in Check 1. Its reverse {v1,v0} must not appear
        // in the triangle fan at all (it would mean the boundary was
        // traversed in the wrong orientation by some triangle).
        std::set<std::pair<int, int>> boundaryReverseSet;
        for (const EdgeDir& be : boundaryEdges) {
          boundaryReverseSet.insert({be.v1, be.v0});
        }
        // Check 3: every directed edge in the triangle fan that is NOT a
        // required-direction boundary edge must:
        //   (a) not be the reverse of any boundary edge, and
        //   (b) have count exactly 1, and
        //   (c) its reverse must also have count exactly 1
        //       (the two directions of an interior undirected edge).
        if (!validFail) {
          for (const auto& [ep, cnt] : edgeCount) {
            if (boundaryEdgeSet.count(ep)) continue;
            // ep is not in the required-direction boundary set.
            // Reject if it is the reverse of a boundary edge.
            if (boundaryReverseSet.count(ep)) {
              validFail = true;
              break;
            }
            // ep is an interior directed edge; must appear exactly once.
            if (cnt != 1) {
              validFail = true;
              break;
            }
            // Its reverse must also appear exactly once.
            auto revIt = edgeCount.find({ep.second, ep.first});
            if (revIt == edgeCount.end() || revIt->second != 1) {
              validFail = true;
              break;
            }
          }
        }
      }

      // Full PSLG embedding checks (all four, in the face frame).
      // Regions exceeding the cap already gated above; this block only runs
      // when validFail is false and the cap was not exceeded.
      if (!validFail) {
        // Collect undirected triangle edges (use canonical {min,max} ordering
        // so each undirected edge is represented once for the duplicate check).
        struct Seg2D {
          vec2 a, b;
          int va, vb;
        };
        // Directed edges for crossing/point-on-edge checks (all 3*|T| of them).
        std::vector<Seg2D> allEdgesDir;
        allEdgesDir.reserve(3 * tris.size());
        for (const ivec3& t : tris) {
          for (int k = 0; k < 3; ++k) {
            allEdgesDir.push_back(
                {p2(t[k]), p2(t[(k + 1) % 3]), t[k], t[(k + 1) % 3]});
          }
        }
        // Undirected edge set for checks (c) no collinear overlap and
        // (d) no duplicate undirected edges beyond coverage.
        // Key: {min(va,vb), max(va,vb)}.
        std::map<std::pair<int, int>, int> undirEdgeCount;
        for (const Seg2D& s : allEdgesDir) {
          const auto key =
              std::make_pair(std::min(s.va, s.vb), std::max(s.va, s.vb));
          auto it = undirEdgeCount.find(key);
          if (it == undirEdgeCount.end())
            undirEdgeCount.insert({key, 1});
          else
            ++it->second;
        }

        auto cross2d = [](const vec2& u, const vec2& v) {
          return u.x * v.y - u.y * v.x;
        };
        // pointOnSegExact: does P lie on segment AB (incl. endpoints)?
        // Exact (zero-tolerance) cross product, parametric range check.
        auto pointOnSegExact = [&](const vec2& P, const vec2& A,
                                   const vec2& B) {
          const vec2 ab = {B.x - A.x, B.y - A.y};
          const vec2 ap = {P.x - A.x, P.y - A.y};
          if (cross2d(ab, ap) != 0.0) return false;
          const double dot = ab.x * ap.x + ab.y * ap.y;
          const double len2 = ab.x * ab.x + ab.y * ab.y;
          return dot >= 0.0 && dot <= len2;
        };
        // collinearOverlap: do two collinear segments [A,B] and [C,D]
        // overlap (share more than a single point)?
        // Assumes segments are already known to be collinear.
        auto collinearOverlap = [&](const vec2& A, const vec2& B, const vec2& C,
                                    const vec2& D) {
          const vec2 ab = {B.x - A.x, B.y - A.y};
          const double len2 = ab.x * ab.x + ab.y * ab.y;
          if (len2 == 0.0) return false;
          // Project C and D onto AB; overlap iff intervals [0,1] and
          // [tC,tD] share an interior point (tC != tD to exclude
          // single-point touch at an endpoint).
          const double tC = (ab.x * (C.x - A.x) + ab.y * (C.y - A.y)) / len2;
          const double tD = (ab.x * (D.x - A.x) + ab.y * (D.y - A.y)) / len2;
          const double lo = std::min(tC, tD);
          const double hi = std::max(tC, tD);
          // Overlap of [0,1] and [lo,hi] has nonzero length iff lo < 1 && hi >
          // 0 AND the overlap interval has positive measure (hi - lo > 0 or the
          // overlap [max(0,lo), min(1,hi)] has length > 0).
          return lo < 1.0 && hi > 0.0 && std::max(0.0, lo) < std::min(1.0, hi);
        };

        const int nSeg = static_cast<int>(allEdgesDir.size());
        for (int i = 0; i < nSeg && !validFail; ++i) {
          const Seg2D& si = allEdgesDir[i];
          for (int j = i + 1; j < nSeg && !validFail; ++j) {
            const Seg2D& sj = allEdgesDir[j];
            // Shared-vertex flags for adjacency detection.
            const bool shareAA = (si.va == sj.va);
            const bool shareAB = (si.va == sj.vb);
            const bool shareBA = (si.vb == sj.va);
            const bool shareBB = (si.vb == sj.vb);
            const bool adjacent = shareAA || shareAB || shareBA || shareBB;
            if (!adjacent) {
              // Non-adjacent: any contact is a PSLG violation.
              if (pointOnSegExact(si.a, sj.a, sj.b)) {
                validFail = true;  // (b) endpoint si.a on sj
                break;
              }
              if (pointOnSegExact(si.b, sj.a, sj.b)) {
                validFail = true;  // (b) endpoint si.b on sj
                break;
              }
              if (pointOnSegExact(sj.a, si.a, si.b)) {
                validFail = true;  // (b) endpoint sj.a on si
                break;
              }
              if (pointOnSegExact(sj.b, si.a, si.b)) {
                validFail = true;  // (b) endpoint sj.b on si
                break;
              }
              // (a) strict crossing.
              const vec2 ab = {si.b.x - si.a.x, si.b.y - si.a.y};
              const vec2 ac = {sj.a.x - si.a.x, sj.a.y - si.a.y};
              const vec2 ad = {sj.b.x - si.a.x, sj.b.y - si.a.y};
              const double t1 = cross2d(ab, ac);
              const double t2 = cross2d(ab, ad);
              if (t1 * t2 < 0) {
                const vec2 cd = {sj.b.x - sj.a.x, sj.b.y - sj.a.y};
                const vec2 ca = {si.a.x - sj.a.x, si.a.y - sj.a.y};
                const vec2 cb = {si.b.x - sj.a.x, si.b.y - sj.a.y};
                const double t3 = cross2d(cd, ca);
                const double t4 = cross2d(cd, cb);
                if (t3 * t4 < 0) {
                  validFail = true;  // (a) strict crossing
                  break;
                }
              }
              // (c) collinear overlap: both cross products zero.
              if (t1 == 0.0 && t2 == 0.0 &&
                  collinearOverlap(si.a, si.b, sj.a, sj.b)) {
                validFail = true;
                break;
              }
            } else {
              // Adjacent (share >= 1 vertex). The shared-endpoint contact
              // is legitimate. Reject only if a NON-shared endpoint of
              // either segment lies on the other segment. This detects
              // vertex-on-nonincident-edge violations that are hidden
              // behind a shared vertex (e.g., A-B-C collinear with C on
              // A->B, or A lying strictly on B->C).
              //
              // For nShared == 2 (same undirected edge, opposite directions),
              // all four endpoints are shared; we skip all endpoint checks.
              // Boundary coverage Check 3 catches true duplicates (same
              // directed edge appearing twice with cnt != 1), and the
              // collinear check is NOT applied here: two opposite-direction
              // copies of an interior edge are valid in a triangulation.
              if (!shareAA && !shareAB) {
                // si.a is not a shared vertex; check si.a on sj.
                if (pointOnSegExact(si.a, sj.a, sj.b)) {
                  validFail = true;
                  break;
                }
              }
              if (!shareBA && !shareBB) {
                // si.b is not a shared vertex; check si.b on sj.
                if (pointOnSegExact(si.b, sj.a, sj.b)) {
                  validFail = true;
                  break;
                }
              }
              if (!shareAA && !shareBA) {
                // sj.a is not a shared vertex; check sj.a on si.
                if (pointOnSegExact(sj.a, si.a, si.b)) {
                  validFail = true;
                  break;
                }
              }
              if (!shareAB && !shareBB) {
                // sj.b is not a shared vertex; check sj.b on si.
                if (pointOnSegExact(sj.b, si.a, si.b)) {
                  validFail = true;
                  break;
                }
              }
            }
          }
        }
        // (d) no duplicate undirected edges beyond the coverage rule:
        // boundary edges appear once (one directed copy) => undirected
        // count = 2 at most (one for each direction), but a boundary edge
        // has exactly one direction in the triangulation, so its undirected
        // count is 1. Interior edges have both directions, count = 2.
        // Any undirected edge with count > 2 is a duplicate.
        if (!validFail) {
          for (const auto& [key, cnt] : undirEdgeCount) {
            if (cnt > 2) {
              validFail = true;
              break;
            }
          }
        }
      }

      if (validFail) {
        out.interiorIslandVerts = 1;
        out.polygons.clear();
        return out;
      }

      // Append triangles as ordinary polygons in face-winding orientation
      // (MergePolygons computes signed multiplicity itself).
      for (const ivec3& t : tris) {
        newPolygons.push_back({t[0], t[1], t[2]});
      }
    }

    // Pass through non-holed negative cycles (should be empty for clean
    // islands, but defensive).
    for (size_t hi = 0; hi < holes.size(); ++hi) {
      bool assigned = false;
      for (const HoleAssignment& ha : assignments) {
        if (ha.holeIdx == hi) {
          assigned = true;
          break;
        }
      }
      if (!assigned) newPolygons.push_back(holes[hi].cyc);
    }

    out.polygons = std::move(newPolygons);
  }

  return out;
}

std::vector<MergedPolygon> MergePolygons(
    const std::vector<std::pair<int, std::vector<int>>>& facePolygons) {
  // Canonical form: the lexicographically-smallest rotation. O(n^2)
  // per cycle - sub-polygon cycles are small.
  auto minRotation = [](const std::vector<int>& c) {
    const size_t n = c.size();
    std::vector<int> best;
    std::vector<int> rot(n);
    for (size_t s = 0; s < n; ++s) {
      for (size_t i = 0; i < n; ++i) rot[i] = c[(s + i) % n];
      if (best.empty() || rot < best) best = rot;
    }
    return best;
  };
  struct Entry {
    int mult = 0;
    int face = -1;
  };
  std::map<std::vector<int>, Entry> merged;
  for (const auto& [face, cycle] : facePolygons) {
    if (cycle.size() < 3) continue;  // defensive; the partition emits >= 3
    const std::vector<int> fwd = minRotation(cycle);
    const std::vector<int> rev(cycle.rbegin(), cycle.rend());
    const std::vector<int> bwd = minRotation(rev);
    // A simple cycle with distinct verts is never rotation-equivalent
    // to its own reversal, so fwd != bwd and the sign is well-defined:
    // +1 when the canonical form comes from the cycle as walked, -1
    // from the reversal. A palindromic equality would mean the
    // partition emitted a non-simple cycle - assert, do not miscount.
    DEBUG_ASSERT(fwd != bwd, logicErr,
                 "MergePolygons: cycle equals its own reversal");
    const bool forward = fwd < bwd;
    Entry& e = merged[forward ? fwd : bwd];
    if (e.face < 0) e.face = face;
    e.mult += forward ? 1 : -1;
  }
  std::vector<MergedPolygon> out;
  out.reserve(merged.size());
  for (const auto& [key, e] : merged) {
    if (e.mult == 0) continue;  // coincident opposite-facing pair cancels
    out.push_back({key, e.mult, e.face});
  }
  return out;
}

CellComplex BuildCellComplex(const Manifold::Impl& impl,
                             const std::vector<MergedPolygon>& polygons,
                             const std::vector<vec3>& newVertPositions) {
  using la::cross;
  using la::dot;
  CellComplex out;
  const int nP = static_cast<int>(polygons.size());
  if (nP == 0) return out;
  const int baseId = static_cast<int>(impl.NumVert());
  auto posOf = [&](int id) {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Newell normal per polygon: the canonical cycle's own orientation -
  // the frame its multiplicity is signed against (the polygon's FRONT
  // is the +normal side). Computed over positions RELATIVE to the
  // cycle's first vert: the cyclic sum is mathematically identical,
  // but absolute positions cancel catastrophically for an eps-thin
  // sliver far from the origin (cross terms ~ R * L with signal
  // ~ L^2), flattening real normals to exact zero.
  std::vector<vec3> normal(nP);
  for (int p = 0; p < nP; ++p) {
    const std::vector<int>& c = polygons[p].cycle;
    const vec3 origin = posOf(c[0]);
    vec3 nsum(0.0, 0.0, 0.0);
    for (size_t i = 0; i < c.size(); ++i) {
      nsum = nsum +
             cross(posOf(c[i]) - origin, posOf(c[(i + 1) % c.size()]) - origin);
    }
    const double len2 = dot(nsum, nsum);
    DEBUG_ASSERT(len2 > 0, logicErr, "BuildCellComplex: degenerate polygon");
    normal[p] = len2 > 0 ? nsum / std::sqrt(len2) : vec3(0.0, 0.0, 1.0);
  }

  // Fan entries per undirected arrangement edge. The in-face direction
  // at a cycle edge a -> b is cross(normal, walkDir): perpendicular to
  // the edge, pointing into the polygon (exact locally, convex or
  // not). frontCcw: the front side faces the CCW-adjacent wedge iff
  // rotating the in-face direction +90 degrees about the edge axis
  // lands on the +normal side.
  struct FanEntry {
    int polygon;
    vec2 dir2;
    bool frontCcw;
  };
  std::map<std::pair<int, int>, std::vector<FanEntry>> fans;
  for (int p = 0; p < nP; ++p) {
    const std::vector<int>& c = polygons[p].cycle;
    for (size_t i = 0; i < c.size(); ++i) {
      const int a = c[i];
      const int b = c[(i + 1) % c.size()];
      const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
      const vec3 pa = posOf(key.first);
      const vec3 pb = posOf(key.second);
      const vec3 axisRaw = pb - pa;
      const double axisLen2 = dot(axisRaw, axisRaw);
      DEBUG_ASSERT(axisLen2 > 0, logicErr,
                   "BuildCellComplex: zero-length arrangement edge");
      if (axisLen2 == 0) continue;
      const vec3 axis = axisRaw / std::sqrt(axisLen2);
      const vec3 walkRaw = posOf(b) - posOf(a);
      const vec3 walk = walkRaw / std::sqrt(dot(walkRaw, walkRaw));
      const vec3 inFace = cross(normal[p], walk);
      const InPlaneBasis fanBasis = FaceBasisFromNormal(axis);
      const vec2 dir2(dot(inFace, fanBasis.u), dot(inFace, fanBasis.v));
      const bool frontCcw = dot(cross(axis, inFace), normal[p]) > 0;
      fans[key].push_back({p, dir2, frontCcw});
    }
  }

  // Radial sort per fan (atan2-free comparator); an exact angular tie
  // is a step-12 invariant failure - DEBUG_ASSERT, with ascending
  // polygon id only as release determinism insurance.
  auto bucketOf = [](const vec2& d) {
    return (d.y > 0 || (d.y == 0 && d.x > 0)) ? 0 : 1;
  };
  for (auto& [key, entries] : fans) {
    std::sort(entries.begin(), entries.end(),
              [&](const FanEntry& x, const FanEntry& y) {
                const int bx = bucketOf(x.dir2);
                const int by = bucketOf(y.dir2);
                if (bx != by) return bx < by;
                const double c = x.dir2.x * y.dir2.y - x.dir2.y * y.dir2.x;
                if (c != 0) return c > 0;
                DEBUG_ASSERT(false, logicErr,
                             "BuildCellComplex: exact angular tie (step-12 "
                             "invariant failure)");
                return x.polygon < y.polygon;
              });
  }

  // Wedges -> cells: between angularly-consecutive entries lies one
  // wedge; unite the CCW-facing side of the earlier entry with the
  // CW-facing side of the later. A k = 1 fan is an open sheet's rim:
  // its single wedge wraps around and unites the polygon's own front
  // and back, as the ambient space does.
  DisjointSets uf(static_cast<uint32_t>(2 * nP));
  for (const auto& [key, entries] : fans) {
    const size_t k = entries.size();
    for (size_t i = 0; i < k; ++i) {
      const size_t j = (i + 1) % k;
      const int sideI = entries[i].frontCcw ? 0 : 1;  // faces CCW
      const int sideJ = entries[j].frontCcw ? 1 : 0;  // faces CW
      uf.unite(static_cast<uint32_t>(2 * entries[i].polygon + sideI),
               static_cast<uint32_t>(2 * entries[j].polygon + sideJ));
    }
  }

  // Renumber cells by smallest member key; emit fans ordered by edge.
  out.polySide2Cell.assign(2 * nP, -1);
  std::map<uint32_t, int> cellIdOf;
  for (int s = 0; s < 2 * nP; ++s) {
    const uint32_t r = uf.find(static_cast<uint32_t>(s));
    auto it = cellIdOf.find(r);
    if (it == cellIdOf.end()) {
      it = cellIdOf.insert({r, out.numCells++}).first;
    }
    out.polySide2Cell[s] = it->second;
  }
  out.fans.reserve(fans.size());
  for (const auto& [key, entries] : fans) {
    EdgeFan fan;
    fan.a = key.first;
    fan.b = key.second;
    fan.polygons.reserve(entries.size());
    fan.frontCcw.reserve(entries.size());
    for (const FanEntry& e : entries) {
      fan.polygons.push_back(e.polygon);
      fan.frontCcw.push_back(e.frontCcw);
    }
    out.fans.push_back(std::move(fan));
  }
  return out;
}

CellWinding ClassifyCells(const Manifold::Impl& impl,
                          const std::vector<MergedPolygon>& polygons,
                          const std::vector<vec3>& newVertPositions,
                          const CellComplex& cells, double epsHint) {
  using la::cross;
  using la::dot;
  CellWinding out;
  const int nP = static_cast<int>(polygons.size());
  out.winding.assign(cells.numCells, 0);
  out.keep.assign(nP, false);
  out.flip.assign(nP, false);
  if (nP == 0) {
    out.ok = true;
    return out;
  }
  const int baseId = static_cast<int>(impl.NumVert());
  auto posOf = [&](int id) {
    return GetPos3(id, baseId, impl, newVertPositions);
  };

  // Canonical Newell normals (the frame each signed multiplicity is
  // measured against; relative-origin sum - see BuildCellComplex),
  // the arrangement bbox, and a length-correct graze margin: the
  // impl's epsilon when it has one, else machine eps at the
  // arrangement's own scale.
  std::vector<vec3> normal(nP);
  std::vector<double> polyArea2(nP, 0.0);  // (2 * area)^2, target ordering
  vec3 bbMin = posOf(polygons[0].cycle[0]);
  vec3 bbMax = bbMin;
  for (int p = 0; p < nP; ++p) {
    const std::vector<int>& cyc = polygons[p].cycle;
    const vec3 origin = posOf(cyc[0]);
    vec3 nsum(0.0, 0.0, 0.0);
    for (size_t i = 0; i < cyc.size(); ++i) {
      const vec3 pa = posOf(cyc[i]);
      nsum =
          nsum + cross(pa - origin, posOf(cyc[(i + 1) % cyc.size()]) - origin);
      bbMin = la::min(bbMin, pa);
      bbMax = la::max(bbMax, pa);
    }
    const double len2 = dot(nsum, nsum);
    DEBUG_ASSERT(len2 > 0, logicErr, "ClassifyCells: degenerate polygon");
    normal[p] = len2 > 0 ? nsum / std::sqrt(len2) : vec3(0.0, 0.0, 1.0);
    polyArea2[p] = len2;
  }
  const double scale = la::length(bbMax - bbMin);
  const double eps = std::max(
      {epsHint, impl.epsilon_, std::numeric_limits<double>::epsilon() * scale});
  const vec3 p0 = 0.5 * (bbMin + bbMax) + kSeedCastDir * (2.0 * scale);

  // Ear decomposition for the cast's crossing tests: a FAN from each
  // cycle's first vert. Fan ears tile any simple polygon as a SIGNED
  // winding decomposition (negative where a concavity puts the fan
  // outside), and the cast kernel's direction-based step makes
  // opposite-sign coverage cancel exactly - so the per-polygon
  // crossing sum is right with no triangulator and no degeneracy
  // checks (near-line sliver cycles made manifold::Triangulate's CCW
  // check throw). Exactly-zero ears contribute exactly nothing and
  // are skipped; near-zero ears stay and at worst graze a cast into
  // its retry.
  std::vector<std::vector<ivec3>> ears(nP);
  for (int p = 0; p < nP; ++p) {
    const std::vector<int>& cyc = polygons[p].cycle;
    const int n = static_cast<int>(cyc.size());
    const vec3 origin = posOf(cyc[0]);
    for (int i = 1; i + 1 < n; ++i) {
      const vec3 e1 = posOf(cyc[i]) - origin;
      const vec3 e2 = posOf(cyc[i + 1]) - origin;
      const vec3 c = cross(e1, e2);
      if (dot(c, c) == 0.0) continue;  // exact-zero ear: contributes nothing
      ears[p].push_back(ivec3(0, i, i + 1));
    }
  }

  // Cell graph: polygon p joins its front cell (the +canonical-normal
  // side, key 2p) to its back (2p + 1); crossing front to back moves
  // AGAINST the normal - entering what the surface element wraps - so
  // the winding gains the signed multiplicity. A polygon whose two
  // sides united (an open sheet, k = 1 rim) separates nothing and
  // propagates nothing; the keep rule below never keeps it.
  std::vector<std::vector<std::pair<int, int>>> adj(cells.numCells);
  for (int p = 0; p < nP; ++p) {
    const int cF = cells.polySide2Cell[2 * p];
    const int cB = cells.polySide2Cell[2 * p + 1];
    if (cF == cB) continue;
    adj[cF].push_back({cB, polygons[p].mult});
    adj[cB].push_back({cF, -polygons[p].mult});
  }

  // Per connected component of the cell graph: seed one cell by a
  // segment cast, then propagate by BFS. Components are discovered in
  // ascending cell order; targets are tried in DESCENDING area, ties
  // ascending polygon id - deterministic.
  std::vector<int> compOf(cells.numCells, -1);
  std::vector<bool> seen(cells.numCells, false);
  int numComps = 0;
  for (int c0 = 0; c0 < cells.numCells; ++c0) {
    if (compOf[c0] >= 0) continue;
    const int comp = numComps++;
    std::vector<int> compCells{c0};
    compOf[c0] = comp;
    for (size_t i = 0; i < compCells.size(); ++i) {
      for (const auto& [d, m] : adj[compCells[i]]) {
        if (compOf[d] < 0) {
          compOf[d] = comp;
          compCells.push_back(d);
        }
      }
    }
    // Separating polygons of this component, as cast targets, in
    // DESCENDING-AREA order (ties to ascending polygon id): a big
    // polygon's interior point sits far from its boundary, so the
    // first target almost always casts cleanly; slivers sort last. A
    // component with no separating polygon (a lone open sheet) needs
    // no seed: its polygons are dropped whatever the winding.
    std::vector<int> candidates;
    for (int p = 0; p < nP; ++p) {
      if (compOf[cells.polySide2Cell[2 * p]] == comp &&
          cells.polySide2Cell[2 * p] != cells.polySide2Cell[2 * p + 1]) {
        candidates.push_back(p);
      }
    }
    if (candidates.empty()) continue;
    std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
      if (polyArea2[a] != polyArea2[b]) return polyArea2[a] > polyArea2[b];
      return a < b;
    });
    bool seeded = false;
    int tried = 0;
    for (const int q : candidates) {
      if (tried == kSeedCastMaxTargets) break;
      ++tried;
      ++out.seedCasts;
      // Target: an interior point of q. Triangles use their centroid
      // directly; longer cycles take the largest ear's centroid from
      // a REAL triangulation (a concave polygon's vert-centroid - and
      // a fan ear's - can fall outside it). A triangulator throw on a
      // near-line sliver cycle just skips to the next target.
      const std::vector<int>& cyc = polygons[q].cycle;
      vec3 target;
      if (cyc.size() == 3) {
        target = (posOf(cyc[0]) + posOf(cyc[1]) + posOf(cyc[2])) / 3.0;
      } else {
        const InPlaneBasis basis = FaceBasisFromNormal(normal[q]);
        const vec3 origin = posOf(cyc[0]);
        SimplePolygon poly2;
        poly2.reserve(cyc.size());
        for (const int v : cyc) {
          const vec3 d = posOf(v) - origin;
          poly2.push_back(vec2(dot(d, basis.u), dot(d, basis.v)));
        }
        std::vector<ivec3> tris;
        // eps is the effective classifier epsilon (>= impl.epsilon_ and
        // the pipeline epsHint): cycles here carry pipeline-eps-scale
        // jitter, so triangulating at the raw mesh epsilon could fail
        // on cycles the pipeline considers clean. A MANIFOLD_DEBUG
        // throw on a near-line sliver skips to the next target;
        // release Triangulate does not throw.
#ifdef MANIFOLD_DEBUG
        try {
          tris = Triangulate({poly2}, std::max(impl.tolerance_, eps), true);
        } catch (...) {
          continue;
        }
#else
        tris = Triangulate({poly2}, std::max(impl.tolerance_, eps), true);
#endif
        if (tris.empty()) continue;
        double bestA = -1.0;
        ivec3 best = tris[0];
        for (const ivec3& t : tris) {
          const vec2 a2 = poly2[t[1]] - poly2[t[0]];
          const vec2 b2 = poly2[t[2]] - poly2[t[0]];
          const double a = std::fabs(a2.x * b2.y - a2.y * b2.x);
          if (a > bestA) {
            bestA = a;
            best = t;
          }
        }
        target =
            (posOf(cyc[best[0]]) + posOf(cyc[best[1]]) + posOf(cyc[best[2]])) /
            3.0;
      }
      // The cast must arrive transversally: P0 within eps of q's own
      // plane is a tangential arrival - retry.
      if (std::fabs(dot(p0 - target, normal[q])) <= eps) continue;
      // Count signed crossings of the open segment P0 -> target
      // against ALL other polygons (other components' included - the
      // true ambient winding is exactly what a nested component
      // cannot learn from its own polygons). Q itself is excluded:
      // counting it would measure the far side.
      int wArr = 0;
      bool graze = false;
      for (int p = 0; p < nP && !graze; ++p) {
        if (p == q) continue;
        // An open sheet (front == back cell) separates nothing and is
        // excluded from the BFS - the cast must skip it too, or the
        // seeded winding disagrees with what propagates from it.
        if (cells.polySide2Cell[2 * p] == cells.polySide2Cell[2 * p + 1])
          continue;
        // SIGNED sum over the fan ears: where a concavity makes fan
        // ears overlap, the opposite-orientation hits cancel exactly
        // (the winding-decomposition argument behind the fan choice).
        int stepSum = 0;
        for (const ivec3& e : ears[p]) {
          const CastResult r =
              CastSegmentAtEar(p0, target, posOf(polygons[p].cycle[e[0]]),
                               posOf(polygons[p].cycle[e[1]]),
                               posOf(polygons[p].cycle[e[2]]), eps);
          if (r.kind == CastHit::kGraze) {
            graze = true;
            break;
          }
          if (r.kind == CastHit::kHit) stepSum += r.step;
        }
        if (!graze && stepSum != 0) wArr += stepSum * polygons[p].mult;
      }
      if (graze) continue;
      // Seed the side the segment arrives through: the front iff the
      // canonical normal points back along the arrival direction.
      const int side = dot(target - p0, normal[q]) < 0 ? 0 : 1;
      const int seedCell = cells.polySide2Cell[2 * q + side];
      out.winding[seedCell] = wArr;
      seen[seedCell] = true;
      std::vector<int> frontier{seedCell};
      for (size_t i = 0; i < frontier.size(); ++i) {
        const int c = frontier[i];
        for (const auto& [d, m] : adj[c]) {
          if (!seen[d]) {
            seen[d] = true;
            out.winding[d] = out.winding[c] + m;
            frontier.push_back(d);
          } else if (out.winding[d] != out.winding[c] + m) {
            // A disagreement means the arrangement is not the closed
            // surface the propagation assumes - fail the
            // classification in release too, so the driver falls
            // back (matches BuildEmitTopology).
            DEBUG_ASSERT(false, logicErr,
                         "ClassifyCells: winding propagation disagreement");
            out.ok = false;
            return out;
          }
        }
      }
      seeded = true;
      break;
    }
    if (!seeded) {
      // Every target in the retry budget grazed: a degenerate
      // configuration the caller handles by falling back to the input.
      DEBUG_ASSERT(false, logicErr,
                   "ClassifyCells: seed cast retries exhausted");
      out.ok = false;
      return out;
    }
  }

  // Keep a polygon iff exactly one side is inside (winding > 0);
  // orient kept polygons with the normal toward the outside cell, so
  // flip exactly those whose front is the inside.
  for (int p = 0; p < nP; ++p) {
    const bool inFront = out.winding[cells.polySide2Cell[2 * p]] > 0;
    const bool inBack = out.winding[cells.polySide2Cell[2 * p + 1]] > 0;
    out.keep[p] = inFront != inBack;
    out.flip[p] = out.keep[p] && inFront;
  }
  out.ok = true;
  return out;
}

EmitTopology BuildEmitTopology(const std::vector<MergedPolygon>& polygons,
                               const CellComplex& cells,
                               const CellWinding& winding) {
  EmitTopology out;
  const int nP = static_cast<int>(polygons.size());
  // Kept polygons in output (outward) orientation, arrangement ids.
  std::vector<int> keptIdxOf(nP, -1);
  for (int p = 0; p < nP; ++p) {
    if (!winding.keep[p]) continue;
    keptIdxOf[p] = static_cast<int>(out.keptPolygons.size());
    out.keptPolygons.push_back(p);
    std::vector<int> cyc = polygons[p].cycle;
    if (winding.flip[p]) std::reverse(cyc.begin(), cyc.end());
    out.outCycles.push_back(std::move(cyc));
  }
  const int nK = static_cast<int>(out.keptPolygons.size());
  if (nK == 0) {
    out.ok = true;  // nothing kept: an empty emit is consistent
    return out;
  }

  // Halfedge ids: cycleStart[k] + cycle position. A simple cycle
  // visits an undirected edge at most once, so (kept polygon, edge)
  // names a halfedge uniquely.
  std::vector<int> cycleStart(nK + 1, 0);
  for (int k = 0; k < nK; ++k) {
    cycleStart[k + 1] =
        cycleStart[k] + static_cast<int>(out.outCycles[k].size());
  }
  const int nH = cycleStart[nK];
  std::map<std::tuple<int, int, int>, int> posOfEdge;  // (k, lo, hi) -> pos
  std::vector<int> hStart(nH), hNextInCycle(nH), hKept(nH);
  for (int k = 0; k < nK; ++k) {
    const std::vector<int>& cyc = out.outCycles[k];
    const int len = static_cast<int>(cyc.size());
    for (int s = 0; s < len; ++s) {
      const int u = cyc[s];
      const int v = cyc[(s + 1) % len];
      posOfEdge[{k, std::min(u, v), std::max(u, v)}] = s;
      hStart[cycleStart[k] + s] = u;
      hNextInCycle[cycleStart[k] + s] = cycleStart[k] + (s + 1) % len;
      hKept[cycleStart[k] + s] = k;
    }
  }

  // Twin assignment per radial fan, restricted to kept polygons:
  // consecutive kept entries flanking an INSIDE wedge are twins. The
  // wedge CCW of fan entry x is entry x's CCW-facing side's cell;
  // dropped entries inside the kept-to-kept span cannot change
  // inside-ness (equal on both of their sides), so that first cell
  // speaks for the merged wedge. Alternation (each kept polygon
  // flips inside-ness) gives each kept entry exactly one inside-
  // wedge partner; failures mark the topology inconsistent and ok
  // stays false for the caller's fallback.
  std::vector<int> twin(nH, -1);
  bool consistent = true;
  for (const EdgeFan& fan : cells.fans) {
    std::vector<int> keptPos;
    for (size_t i = 0; i < fan.polygons.size(); ++i) {
      if (winding.keep[fan.polygons[i]]) {
        keptPos.push_back(static_cast<int>(i));
      }
    }
    if (keptPos.empty()) continue;
    if (keptPos.size() % 2 != 0) {
      DEBUG_ASSERT(false, logicErr,
                   "BuildEmitTopology: odd kept count at a fan");
      consistent = false;
      continue;
    }
    const size_t kc = keptPos.size();
    for (size_t i = 0; i < kc; ++i) {
      const int x = keptPos[i];
      const int y = keptPos[(i + 1) % kc];
      const int wedgeCell =
          cells.polySide2Cell[2 * fan.polygons[x] + (fan.frontCcw[x] ? 0 : 1)];
      if (winding.winding[wedgeCell] <= 0) continue;  // outside wedge
      const int ka = keptIdxOf[fan.polygons[x]];
      const int kb = keptIdxOf[fan.polygons[y]];
      const auto ita = posOfEdge.find({ka, fan.a, fan.b});
      const auto itb = posOfEdge.find({kb, fan.a, fan.b});
      if (ita == posOfEdge.end() || itb == posOfEdge.end()) {
        DEBUG_ASSERT(false, logicErr,
                     "BuildEmitTopology: fan polygon missing its edge");
        consistent = false;
        continue;
      }
      const int ha = cycleStart[ka] + ita->second;
      const int hb = cycleStart[kb] + itb->second;
      // Twins bound the same inside region with outward normals, so
      // they traverse the shared edge antiparallel.
      const bool aForward = out.outCycles[ka][ita->second] == fan.a;
      const bool bForward = out.outCycles[kb][itb->second] == fan.a;
      if (aForward == bForward || twin[ha] != -1 || twin[hb] != -1) {
        DEBUG_ASSERT(false, logicErr, "BuildEmitTopology: twin conflict");
        consistent = false;
        continue;
      }
      twin[ha] = hb;
      twin[hb] = ha;
    }
  }
  for (int h = 0; h < nH && consistent; ++h) {
    if (twin[h] < 0) {
      DEBUG_ASSERT(false, logicErr, "BuildEmitTopology: unpaired halfedge");
      consistent = false;
    }
  }
  if (!consistent) return out;

  // Vertex rings: orbits of nextAroundVert(h) =
  // nextInPolygonCycle(twin(h)), walked on the kept-polygon graph
  // (pre-triangulation). One output vert per orbit; the ring-
  // separation argument makes every output edge carry exactly 2
  // halfedges (release-checked below).
  std::vector<int> orbitOf(nH, -1);
  struct Orbit {
    int vert;
    int minPoly;
  };
  std::vector<Orbit> orbits;
  for (int h0 = 0; h0 < nH; ++h0) {
    if (orbitOf[h0] >= 0) continue;
    const int id = static_cast<int>(orbits.size());
    Orbit orb{hStart[h0], std::numeric_limits<int>::max()};
    int h = h0;
    do {
      orbitOf[h] = id;
      orb.minPoly = std::min(orb.minPoly, out.keptPolygons[hKept[h]]);
      DEBUG_ASSERT(hStart[h] == orb.vert, logicErr,
                   "BuildEmitTopology: orbit left its vert");
      h = hNextInCycle[twin[h]];
    } while (h != h0);
    orbits.push_back(orb);
  }
  // Deterministic ring numbering by (geometric vert, smallest
  // incident kept polygon) - unique, since a simple cycle starts at
  // a vert once, putting each polygon in one orbit per vert.
  std::vector<int> order(orbits.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    if (orbits[a].vert != orbits[b].vert) {
      return orbits[a].vert < orbits[b].vert;
    }
    return orbits[a].minPoly < orbits[b].minPoly;
  });
  std::vector<int> ringOfOrbit(orbits.size());
  out.ring2Vert.resize(orbits.size());
  for (size_t r = 0; r < order.size(); ++r) {
    DEBUG_ASSERT(r == 0 || orbits[order[r - 1]].vert != orbits[order[r]].vert ||
                     orbits[order[r - 1]].minPoly != orbits[order[r]].minPoly,
                 logicErr, "BuildEmitTopology: duplicate ring key");
    ringOfOrbit[order[r]] = static_cast<int>(r);
    out.ring2Vert[r] = orbits[order[r]].vert;
  }
  for (int k = 0; k < nK; ++k) {
    std::vector<int>& cyc = out.outCycles[k];
    for (size_t s = 0; s < cyc.size(); ++s) {
      cyc[s] = ringOfOrbit[orbitOf[cycleStart[k] + static_cast<int>(s)]];
    }
  }

  // Belt-and-suspenders: exactly 2 halfedges per output (ring, ring)
  // edge, in release too - a violation feeds the driver's gate.
  std::map<std::pair<int, int>, int> edgeCount;
  for (const std::vector<int>& cyc : out.outCycles) {
    for (size_t s = 0; s < cyc.size(); ++s) {
      const int u = cyc[s];
      const int v = cyc[(s + 1) % cyc.size()];
      ++edgeCount[{std::min(u, v), std::max(u, v)}];
    }
  }
  for (const auto& [e, n] : edgeCount) {
    if (n != 2) {
      DEBUG_ASSERT(false, logicErr,
                   "BuildEmitTopology: output edge without exactly 2 "
                   "halfedges");
      return out;
    }
  }
  out.ok = true;
  return out;
}

void PropagateNewVertsToOnEdgeLists(
    const Manifold::Impl& impl, const std::vector<vec3>& newVertPositions,
    const std::vector<EdgeTriIntersection>& etIsects,
    const std::vector<int>& etIsect2Vert, const std::vector<Edge>& edges,
    std::vector<EdgeVertList>& onEdgeLists) {
  using la::dot;
  const int baseId = static_cast<int>(impl.NumVert());
  // For each etIsect, add the resolved vert id to the on-edge list of
  // the piercing edge. t is recomputed from the RESOLVED position, not
  // taken from the raw event parameter x.s: a snapped event's vert
  // sits up to the snap radius from the event point, and two nearby
  // events' stale parameters could order against the geometry (the
  // partition assumes monotone t). Skip if the vert is the edge's
  // endpoint or already in the list.
  std::vector<int> touched;
  for (size_t i = 0; i < etIsects.size(); ++i) {
    const auto& x = etIsects[i];
    auto& list = onEdgeLists[x.edgeIdx];
    const int v = etIsect2Vert[i];
    if (v == edges[x.edgeIdx].v0 || v == edges[x.edgeIdx].v1) continue;
    if (std::find(list.verts.begin(), list.verts.end(), v) != list.verts.end())
      continue;
    const vec3 a = impl.vertPos_[edges[x.edgeIdx].v0];
    const vec3 ab = impl.vertPos_[edges[x.edgeIdx].v1] - a;
    const double abLen2 = dot(ab, ab);
    const vec3 p = GetPos3(v, baseId, impl, newVertPositions);
    const double t = abLen2 > 0 ? dot(p - a, ab) / abLen2 : x.s;
    // A snapped vert can project outside the edge's interior when the
    // edge is short (~the snap radius): it then subdivides nothing
    // here - the same (0, 1) interior rule the on-edge builders use.
    if (t <= 0.0 || t >= 1.0) continue;
    list.verts.push_back(v);
    list.ts.push_back(t);
    touched.push_back(x.edgeIdx);
  }
  // BuildOnEdgeVertLists left each list sorted by t, but the appends above
  // are unordered. PartitionFace consumes verts in stored order to
  // build consecutive sub-edges assuming monotone t, so re-sort each
  // touched list by t (carrying verts along).
  std::sort(touched.begin(), touched.end());
  touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
  for (int e : touched) {
    SortVertsByT(onEdgeLists[e].verts, onEdgeLists[e].ts);
  }
}

SelfIntersectionResult CheckSelfIntersection(const Manifold::Impl& m,
                                             double relTol) {
  SelfIntersectionResult r{};
  if (m.IsEmpty()) return r;
  const size_t nTri = m.NumTri();
  r.trianglesTotal = static_cast<int>(nTri);
  if (nTri < 2) return r;

  auto vp = [&](int idx) { return m.vertPos_[idx]; };
  std::vector<Box> triBoxes(nTri);
  std::vector<std::array<int, 3>> triIdx(nTri);
  for (size_t t = 0; t < nTri; ++t) {
    const int i0 = m.halfedge_.Start(3 * t + 0);
    const int i1 = m.halfedge_.Start(3 * t + 1);
    const int i2 = m.halfedge_.Start(3 * t + 2);
    triIdx[t] = {i0, i1, i2};
    Box b(vp(i0), vp(i1));
    b.Union(vp(i2));
    triBoxes[t] = b;
  }

  SortedBVH bvh =
      BuildSortedBVH(VecView<const Box>(triBoxes.data(), triBoxes.size()));

  auto checkPair = [&](size_t qi, size_t li) {
    if (qi >= li) return;
    const size_t ta = bvh.leaf2Orig[qi];
    const size_t tb = bvh.leaf2Orig[li];
    int shared = 0;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        if (triIdx[ta][a] == triIdx[tb][b]) ++shared;
    if (shared >= 2) {
      ++r.adjacentPairsSkipped;
      return;
    }
    ++r.candidatesChecked;
    const vec3 a0 = vp(triIdx[ta][0]), a1 = vp(triIdx[ta][1]),
               a2 = vp(triIdx[ta][2]);
    const vec3 b0 = vp(triIdx[tb][0]), b1 = vp(triIdx[tb][1]),
               b2 = vp(triIdx[tb][2]);
    double maxMag = 0.0;
    auto take = [&](double m) {
      if (m > maxMag) maxMag = m;
    };
    take(SegmentPiercesTriInterior(a0, a1, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(a1, a2, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(a2, a0, b0, b1, b2, relTol));
    take(SegmentPiercesTriInterior(b0, b1, a0, a1, a2, relTol));
    take(SegmentPiercesTriInterior(b1, b2, a0, a1, a2, relTol));
    take(SegmentPiercesTriInterior(b2, b0, a0, a1, a2, relTol));
    if (maxMag > 0.0) {
      ++r.interiorPierces;
      if (maxMag > r.maxPierceMagnitude) r.maxPierceMagnitude = maxMag;
    }
  };

  auto recorder = MakeSimpleRecorder(checkPair);
  auto qf = [&](int i) { return bvh.boxes[i]; };
  bvh.Collisions(recorder, qf, static_cast<int>(nTri));
  return r;
}

// Final-gate helper: does any folded surface enclose real volume?
// Polygons whose two sides united (front cell == back cell, e.g.
// across a k = 1 rim) are dropped by the keep rule. That is correct
// for flat membranes - enclosed volume below area x thickness, and
// post-merge thickness is at most the 10 eps unification radius -
// but a tangent-degenerate contact can fold a CLOSED shell's two
// cells together, and dropping that fold silently deletes the
// shell's material. Folded polygons group into EDGE-CONNECTED
// components (same fold cell + shared undirected edge - the surface
// notion: a vert-touch pinch between two shells must NOT merge them,
// or a positive shell and an inverted twin touching at one snapped
// vert would net their signed volumes to nothing). Per component,
// sum the mult-weighted signed volume (tetra fan anchored at the
// component's own centroid: origin-independent for a closed set,
// near zero for a FLAT open sheet) and compare against the membrane
// bound area x kFoldedVolumePerAreaEps x eps. A macro-BENT open fold
// also trips the gate (its anchored cone volume is O(area x bend
// depth)) - deliberate: an open fold of real extent is arrangement
// damage, and falling back is the safe side.
bool FoldedCellsEncloseVolume(const Manifold::Impl& impl,
                              const std::vector<MergedPolygon>& polys,
                              const std::vector<vec3>& newVertPositions,
                              const CellComplex& cells, double eps) {
  const int baseId = static_cast<int>(impl.NumVert());
  auto posOf = [&](int id) {
    return GetPos3(id, baseId, impl, newVertPositions);
  };
  std::vector<int> folded;  // polygon ids with front == back
  for (size_t p = 0; p < polys.size(); ++p) {
    if (cells.polySide2Cell[2 * p] == cells.polySide2Cell[2 * p + 1]) {
      folded.push_back(static_cast<int>(p));
    }
  }
  if (folded.empty()) return false;
  const uint32_t nF = static_cast<uint32_t>(folded.size());
  DisjointSets uf(nF);
  std::map<std::tuple<int, int, int>, uint32_t> cellEdge2First;
  for (uint32_t k = 0; k < nF; ++k) {
    const int c = cells.polySide2Cell[2 * folded[k]];
    const std::vector<int>& cyc = polys[folded[k]].cycle;
    for (size_t i = 0; i < cyc.size(); ++i) {
      const int u = cyc[i];
      const int v = cyc[(i + 1) % cyc.size()];
      const auto [it, fresh] =
          cellEdge2First.insert({{c, std::min(u, v), std::max(u, v)}, k});
      if (!fresh) uf.unite(it->second, k);
    }
  }
  std::vector<uint32_t> root(nF);
  std::map<uint32_t, uint32_t> root2Comp;
  for (uint32_t k = 0; k < nF; ++k) {
    root[k] = uf.find(k);
    root2Comp.insert({root[k], static_cast<uint32_t>(root2Comp.size())});
  }
  const size_t nComp = root2Comp.size();
  std::vector<vec3> centroid(nComp, vec3(0.0, 0.0, 0.0));
  std::vector<int> centroidVerts(nComp, 0);
  for (uint32_t k = 0; k < nF; ++k) {
    const uint32_t comp = root2Comp[root[k]];
    for (const int v : polys[folded[k]].cycle) {
      centroid[comp] = centroid[comp] + posOf(v);
      ++centroidVerts[comp];
    }
  }
  for (size_t c = 0; c < nComp; ++c) {
    if (centroidVerts[c] > 0) centroid[c] = centroid[c] / centroidVerts[c];
  }
  std::vector<double> volume(nComp, 0.0);
  std::vector<double> area(nComp, 0.0);
  for (uint32_t k = 0; k < nF; ++k) {
    const uint32_t comp = root2Comp[root[k]];
    const std::vector<int>& cyc = polys[folded[k]].cycle;
    const vec3 a0 = posOf(cyc[0]);
    vec3 nsum(0.0, 0.0, 0.0);
    double vol = 0.0;
    for (size_t i = 1; i + 1 < cyc.size(); ++i) {
      const vec3 a = posOf(cyc[i]);
      const vec3 b = posOf(cyc[i + 1]);
      const vec3 cr = la::cross(a - a0, b - a0);
      nsum = nsum + cr;
      vol += la::dot(a0 - centroid[comp], cr) / 6.0;
    }
    volume[comp] += polys[folded[k]].mult * vol;
    area[comp] += 0.5 * std::sqrt(la::dot(nsum, nsum));
  }
  for (size_t c = 0; c < nComp; ++c) {
    if (std::fabs(volume[c]) > area[c] * kFoldedVolumePerAreaEps * eps) {
      return true;
    }
  }
  return false;
}

std::optional<Manifold::Impl> RemoveOverlaps(const Manifold::Impl& input,
                                             double eps,
                                             ExecutionContext::Impl* ctx) {
  // Throwing assertions exist only when BOTH MANIFOLD_ASSERT and
  // MANIFOLD_DEBUG are set (optional_assert.h compiles DEBUG_ASSERT /
  // ASSERT to throws under that conjunction and defines the error
  // types under MANIFOLD_DEBUG; release manifold is exception-free and
  // errors are status enums). The guard keys on MANIFOLD_DEBUG - the
  // superset where the types exist; without MANIFOLD_ASSERT it is dead
  // but harmless. A throwing assertion anywhere in the pipeline -
  // including inside manifold's own Triangulate checks - becomes the
  // nullopt fallback arm (the caller returns its input), preserving
  // pierce-monotonicity; the guard pattern matches polygon.cpp's
  // TriangulateIdxHalfedges.
#ifdef MANIFOLD_DEBUG
  try {
    return RemoveOverlapsImpl(input, eps, ctx);
  } catch (...) {
    return std::nullopt;
  }
#else
  return RemoveOverlapsImpl(input, eps, ctx);
#endif
}

}  // namespace overlap_removal
}  // namespace manifold
