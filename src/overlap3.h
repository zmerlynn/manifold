// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// 3D overlap removal via sweep-native emission (strips + caps).
// Design: docs/SweepEmit3D.md (three crucible rounds + empirical validation).
// Internal seam only; public wiring is post-prototype.

#pragma once

#include <climits>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "boolean2.h"
#include "impl.h"
#include "manifold/common.h"

namespace manifold {

// ---------------------------------------------------------------------------
// Canonicalize-stage output.
// ---------------------------------------------------------------------------

// One canonical face after vert-merge and multiplicity accumulation.
// mult == 0 records are dropped before reaching the pipeline.
struct CanonicalFace {
  int id;        // representative face index in the original Impl
  ivec3 verts;   // canonical merged-vert ids
  vec3 normal;   // from stored vert order
  int64_t mult;  // signed multiplicity (>= 1 or <= -1)
};

// ---------------------------------------------------------------------------
// Fatal-reason taxonomy (spec "FAILURE CONTRACT" - shrunken contract).
// ---------------------------------------------------------------------------

enum class FatalReason {
  SubEpsInput,          // eps <= 0 or degenerate input geometry
  SubEpsFeature,        // macro-scale face in merged sub-eps critical run
  NonManifoldEmission,  // emitted triangulation is not 2-manifold
  // Arrangement too dense/degenerate to section within a sane resource budget -
  // a refusal, not an error.  The whole slab decomposition is held in memory at
  // once, so a near-degenerate input packing thousands of faces into each of
  // tens of thousands of thin slabs would swap-thrash into bad_alloc; the slabs
  // stage fails closed here instead (spec [WALL-B]).
  ArrangementBudget,
  // Regularization operator (docs/Regularize3D.md): a dirty (self-intersecting)
  // component reached candidate B, but B's production dirty-core resolver is
  // not yet built.  Stage-1 fail-closed STUB SENTINEL - the honest outcome
  // while the dirty core is a stub, never a silent wrong result.  When B lands
  // it is replaced by either a clean re-gated resolve or a real re-gate fatal
  // (NonManifoldEmission / a self-intersection re-gate failure).
  DirtyComponentUnresolved,
};

// Non-fatal counter accumulator.
struct Overlap3Counters {
  int subEpsContactsDropped = 0;  // point-like skipped contacts
  int engineIdConflicts = 0;      // total engine id-conflict events; benign
                                  // diagnostic (srcId is unconsumed in 3D) -
                                  // see spec [WALL-B]
  int capArrangements = 0;        // 2D arrangements run by the caps stage
                                  // (one per critical with cap input - the
                                  // M4 pin)
};

// Per-stage result: either a product or a fatal reason.
template <typename T>
struct StageResult {
  std::optional<T> value;
  std::optional<FatalReason> fatal;
  std::string detail;

  bool ok() const { return value.has_value() && !fatal.has_value(); }

  static StageResult<T> Ok(T v) {
    StageResult<T> r;
    r.value = std::move(v);
    return r;
  }
  static StageResult<T> Fatal(FatalReason reason, std::string msg = {}) {
    StageResult<T> r;
    r.fatal = reason;
    r.detail = std::move(msg);
    return r;
  }
};

// ---------------------------------------------------------------------------
// Slabs-stage types.
// ---------------------------------------------------------------------------

// Directed section segment for one straddling face in a slab.
// p1 - p0 has positive dot with yz(cross(+x, face.normal)).
struct SectionFaceSegment {
  int faceId;
  vec2 p0, p1;   // section (y,z) directed as above
  int64_t mult;  // signed multiplicity from CanonicalFace
};

// Per-face track for strip generation (strips stage).
// Records the 3D edge pairs whose interpolations define the face's section
// segment at any x in the slab:
//   at x: p0(x) = Interpolate(va0, vb0, x),  p1(x) = Interpolate(va1, vb1, x)
// xMid evaluation matches (p0, p1) stored in SectionFaceSegment.
struct FaceTrack {
  int faceId;
  vec2 p0, p1;  // section (y,z) at xMid, directed per spec
  vec3 va0,
      vb0;  // 3D edge pair for p0: Interpolate(va0, vb0, x) at any x in slab
  vec3 va1, vb1;  // 3D edge pair for p1
};

// Linearly interpolate (y,z) of segment va-vb at x=xTarget.  This is the
// named EXTRAPOLATION path (unlike shared.h:Interpolate): extending a piece
// backward/forward across an unbuilt run evaluates a track outside
// [va.x, vb.x] by design.  x-degenerate tracks cannot reach here: face tracks
// exclude section-parallel edges (ComputeSectionSegment) and seam tracks use
// exact span membership with xMid strictly between criticals (BuildSlabs);
// the guard is release-safety only.
inline vec2 InterpolateSafe(vec3 va, vec3 vb, double xTarget) {
  const double dx = vb.x - va.x;
  DEBUG_ASSERT(dx != 0.0, logicErr, "x-degenerate track in InterpolateSafe");
  if (dx == 0.0) return {va.y, va.z};
  const double t = (xTarget - va.x) / dx;
  return va.yz() + t * (vb.yz() - va.yz());
}

// Per-seam track for cap/strip extension of class-ii endpoints (spec
// STRIPS/CAPS).
// A seam crossing at yzMid lies on the 3D seam segment [vA, vB]; the correct
// extension to any xTarget is InterpolateSafe(vA, vB, xTarget).yz.
struct SeamTrackEntry {
  vec2 yzMid;   // seam's (y,z) at this slab's xMid (for endpoint lookup)
  vec3 vA, vB;  // 3D seam segment endpoints
};

// Per-slab output.
struct SlabResult {
  double xLo, xHi, xMid;
  bool built;                        // false = sub-eps width, skipped
  std::vector<SweepCapture> pieces;  // retained boundary pieces from engine
  std::vector<SectionFaceSegment> segments;  // directed section segments
  std::vector<FaceTrack> faceTracks;  // per-face tracks for strip extension
  std::vector<SeamTrackEntry> seamTracks;  // per-seam tracks for class-ii ext
  // Test-hook: raw section edges and verts before arrangement.
  std::vector<EdgeM> sectionEdges;
  std::vector<vec2> sectionVerts;
};

// ---------------------------------------------------------------------------
// Seams-stage types.
// ---------------------------------------------------------------------------

// A face-pair seam: 3D segment [vertId0, vertId1] at the intersection of
// two canonical faces.
struct Seam {
  int faceId0, faceId1;
  int vertId0, vertId1;  // canonical merged-vert ids
};

// ---------------------------------------------------------------------------
// Arrangement geometry (canonicalize + seams output; slabs, caps, and
// strips input).
// ---------------------------------------------------------------------------

struct ArrangementGeometry {
  std::vector<vec3> verts;           // all canonical 3D verts
  std::vector<CanonicalFace> faces;  // canonical faces
  std::vector<Seam> seams;           // face-pair seam segments
  // Coplanar groups (spec COPLANAR mechanism 1): face2Group[fi] is the
  // group index of face fi, or -1 when ungrouped.  Grouped faces seed their
  // section edges with id = faces.size() + group (mechanism 2), so
  // coincident in-plane content merges under one source id.
  std::vector<int> face2Group;
  int numGroups = 0;
  // Extra x-criticals with no vert identity: degenerate-contact endpoints and
  // seam-seam crossing x's (spec SEAMS: only the x is consumed;
  // over-inclusion is harmless).  Kept separate from verts - a critical is
  // not a vertex.
  std::vector<double> criticalXs;
};

// ---------------------------------------------------------------------------
// Public entry point.
// ---------------------------------------------------------------------------

struct Overlap3Result {
  std::optional<Manifold::Impl> impl;
  Overlap3Counters counters;
  std::optional<FatalReason> fatal;
  std::string detail;
};

// eps = 0 -> compute from bounding-box scale.
Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps = 0.0);

// ---------------------------------------------------------------------------
// Regularization operator (docs/Regularize3D.md) - parallel entry point.
//
// RegularizeImpl maps a valid oriented face soup to the boundary of the solid
// {p : w_S(p) >= 1}, PER CONNECTED COMPONENT (it never fuses separate DISJOINT
// components - fusion is the Boolean's job, already done upstream).  The one
// exception is a genuine COPLANAR self-overlap that connectivity would split (a
// buried plug with a coincident cap): those components are UNITED at decompose
// time so the fold sees the overlap as internal (touching/disjoint objects with
// no 2D-area overlap never merge).  The pipeline: DECOMPOSE by connectivity (+
// coplanar-overlap merge) -> per-component GATE (validity + IsSelfIntersecting
// + coplanar overlap) -> EARLY-EXIT clean components -> route DIRTY components
// to candidate B -> RE-GATE B's output -> COMPOSE BACK by concatenation.  The
// v3 sweep entry point RemoveOverlaps3D is untouched; this is an additive
// second entry point, not a rewrite.
// ---------------------------------------------------------------------------

// White-box dispatch counters (the Stage-1 pins read these directly).
struct RegularizeCounters {
  int components = 0;  // components after decompose + coplanar-overlap merge
  int clean = 0;       // passed the gate; early-exit copied through
  int dirty = 0;       // failed the gate (self-intersecting OR coplanar
                       // overlap); routed to candidate B
  int regularized =
      0;               // B produced a clean re-gated output (0 while B is stub)
  int failClosed = 0;  // components that fail-closed (dirty stub, re-gate, or
                       // an unexpected non-manifold input component)
};

struct RegularizeResult {
  std::optional<Manifold::Impl> impl;  // composed output; absent on any fatal
  std::optional<FatalReason> fatal;
  std::string detail;
  RegularizeCounters counters;
};

// eps = 0 -> compute from bounding-box scale.
RegularizeResult RegularizeImpl(const Manifold::Impl& in, double eps = 0.0);

// ---------------------------------------------------------------------------
// Test hooks (overlap3_test.cpp only).
// ---------------------------------------------------------------------------

struct Overlap3Internals {
  ArrangementGeometry arr;
  std::vector<SlabResult> slabs;
  std::optional<FatalReason> fatal;
  std::string detail;
  Overlap3Counters counters;
};

// Run canonicalize + seams + slabs; slabs include sectionEdges/sectionVerts
// for gate-2.
Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps = 0.0);

// Candidate B mechanism probe (docs/Regularize3D.md "B's mechanism"), exposed
// so the port of the validated fragment (enumeration + coupled winding) is
// tested directly against the fragment's recorded numbers.  seamCount = genuine
// non-adjacent self-crossings; boundaryTouchPairs = pairs whose deciding
// predicate hit an exact-zero / static-filter-uncertain boundary (the
// single-global-SoS axis); probeWinding[i] = coupled soup winding w_S at
// probes[i] cast to `seed` (kWindingUncertain if a deciding predicate was
// filter-uncertain).
struct CandidateBProbe {
  int seamCount = 0;
  int boundaryTouchPairs = 0;
  int coplanarClusterFaces = 0;  // faces in an exact-coplanar overlap cluster
  std::vector<int> probeWinding;
};
constexpr int kWindingUncertain = INT_MIN;
CandidateBProbe RegularizeB_Probe(const Manifold::Impl& dirty,
                                  const std::vector<vec3>& probes,
                                  const vec3& seed);

// Test hook: run candidate B (fold + build + re-gate) directly on a soup
// treated as ONE dirty component, bypassing decompose and the
// IsSelfIntersecting gate. The self-intersection gate does not flag a pure
// coplanar overlap, so this hook is the way to exercise the exact-coplanar
// fold's RESOLVE path on an isolated coplanar cluster (no transversal
// entanglement).
RegularizeResult RegularizeDirtyDirect(const Manifold::Impl& soup, double eps);

}  // namespace manifold
