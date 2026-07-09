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

// 3D overlap removal by sweep-plane classification, prototype.
// Design: docs/SweepPlane3D.md (crucible round 3, impl crucible closed).
// Internal seam only; public wiring is post-prototype.

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "boolean2.h"
#include "impl.h"
#include "manifold/common.h"

namespace manifold {

// ---------------------------------------------------------------------------
// Stage-A output (spec "Types, constants, and metrics").
// ---------------------------------------------------------------------------

// One canonical face after vert-merge and multiplicity accumulation.
// `id` is the representative triangle index in the original Impl.
// `mult == 0` records are dropped before this reaches the pipeline.
struct CanonicalFace {
  int id;        // representative face index
  ivec3 verts;   // canonical merged-vert ids
  vec3 normal;   // from stored vert order
  int64_t mult;  // signed multiplicity (>= 1 or <= -1)
};

// ---------------------------------------------------------------------------
// Fatal-reason taxonomy (spec "FAILURE CONTRACT").
// ---------------------------------------------------------------------------

enum class FatalReason {
  TripleDiameter,           // triple-point cluster diameter > eps
  SubResolutionChain,       // skipped-contact cluster diameter > eps
  CoplanarOverlap,          // coplanar face overlap (out of scope)
  EdgeInPlane,              // edge lying in another face's plane (out of scope)
  ClassificationAmbiguity,  // no located piece, or pieces disagree
  AnchorConflict,           // degenerate component: conflicting anchors
  UnclassifiableComponent,  // degenerate component: no anchors or span guard
  Starvation,               // eps <= 0 or no usable slab
  BalanceViolation,         // seam edge: kept-region counts unequal
  PSLGInvalid,              // seams cross without shared vert ids
  EngineIdConflict,         // 2D engine conflict counter nonzero
  NonManifoldEmission,      // emitted triangulation is not 2-manifold
};

// Non-fatal counter accumulator.
struct Overlap3Counters {
  int subEpsContactsDropped = 0;  // point-like skipped contacts
  int degenerateClassified = 0;   // degenerate components via anchor
  int epsFeaturesDropped = 0;     // eps-scale degenerate components dropped
  int clearanceSkips = 0;         // pieces skipped for clearance
  int engineIdConflicts = 0;      // total engine id conflict events
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
// Stage-C types (spec "SlabResult", "SectionFaceSegment").
// ---------------------------------------------------------------------------

// Directed section segment for one straddling face in a slab.
// p1 - p0 has positive dot with yz(cross(+x, face.normal)).
struct SectionFaceSegment {
  int faceId;
  vec2 p0, p1;   // section (y,z) directed as above
  int64_t mult;  // signed multiplicity from CanonicalFace
};

// Per-slab output.
struct SlabResult {
  double xLo, xHi, xMid;
  bool built;  // false = sub-eps, skipped
  std::vector<SweepCapture> pieces;
  std::vector<SectionFaceSegment> segments;
  // Test-hook: raw section edges and verts before arrangement.
  std::vector<EdgeM> sectionEdges;
  std::vector<vec2> sectionVerts;
};

// ---------------------------------------------------------------------------
// Stage-B types (seam geometry).
// ---------------------------------------------------------------------------

struct MergedVert {
  vec3 pos;
};

// A seam: the intersection segment of two canonical faces.
// vertIds holds the canonical merged-vert sequence (2 verts for a plain seam,
// more if triple-point insertion added interior verts).
struct Seam {
  int faceId0, faceId1;
  std::vector<int> vertIds;
};

// A subdivided boundary edge for one face: the halfedge's vert sequence
// (including any event/seam verts inserted during stage B).
struct SubdividedEdge {
  int faceId;
  int edgeIdx;               // halfedge index in original mesh
  std::vector<int> vertIds;  // from start to end, both endpoints included
};

// ---------------------------------------------------------------------------
// Arrangement geometry (stage A+B output, stage C+D+E input).
// ---------------------------------------------------------------------------

struct ArrangementGeometry {
  std::vector<MergedVert> verts;            // all canonical 3D verts
  std::vector<CanonicalFace> faces;         // canonical faces
  std::vector<Seam> seams;                  // face-pair seam polylines
  std::vector<SubdividedEdge> subdEdges;    // per halfedge, with inserted verts
  std::vector<std::vector<int>> faceSeams;  // face -> list of seam indices
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
// Test seams (overlap3_test.cpp only).
// ---------------------------------------------------------------------------

// PSLGRegion: one face sub-region produced by the PSLG walk.
// Exposed so tests can synthesise inputs for the white-box functions below.
struct PSLGRegion {
  std::vector<int> loopVerts;               // outer boundary (CCW)
  std::vector<std::vector<int>> holeVerts;  // interior hole loops
  int64_t below = 0, above = 0;             // status-order windings
  bool classified = false;
  bool degenerate = false;  // no covering slab wider than eps
};

struct ClassifyRegionResult {
  std::optional<FatalReason> fatal;
  int64_t below = 0, above = 0;
  bool classified = false;
};

struct Overlap3Internals {
  ArrangementGeometry arr;
  std::vector<SlabResult> slabs;
  std::optional<FatalReason> fatal;
  std::string detail;
  Overlap3Counters counters;
};

// Run stages A+B+C; slabs include sectionEdges/sectionVerts for gate-2.
Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps = 0.0);

// White-box: seam-balance check on synthetic data (P1/P11 pins).
std::optional<FatalReason> CheckSeamBalance_Test(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions);

// White-box: classify one region against given slabs (P2/P3 pins).
ClassifyRegionResult ClassifyRegion_Test(const PSLGRegion& region, int faceId,
                                         const std::vector<SlabResult>& slabs,
                                         const std::vector<MergedVert>& verts,
                                         const CanonicalFace& face, double eps);

// Run stages C+D+E from a pre-built ArrangementGeometry (P5 PSLGInvalid pin).
Overlap3Result RemoveOverlaps3D_FromArr(const ArrangementGeometry& arr,
                                        double eps);

// White-box M6: anchor-component propagation on synthetic data (P8-P10 pins).
std::optional<FatalReason> PropagateAnchorComponents_Test(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt);

}  // namespace manifold
