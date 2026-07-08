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

// 3D overlap removal by sweep-plane classification, prototype implementation.
// Design: docs/SweepPlane3D.md (crucible round 3).
// Internal seam (`RemoveOverlaps3D`) only; public wiring is post-prototype.

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
// Types (spec "Types, constants, and metrics" section)
// ---------------------------------------------------------------------------

// Stage-A output: one canonical face record per unique triangle after merge.
// `id` is the face index in the original Impl (the canonical one chosen by
// centroid-nearest after vert merge). `verts` are the canonical vert ids after
// merge. `normal` is computed from the stored `verts` order. `mult` is the
// signed multiplicity (same orientation = +1 per appearance, opposite = -1).
// Stage A drops records where mult == 0.
struct CanonicalFace {
  int id;  // original face index (canonical representative)
  ivec3 verts;
  vec3 normal;
  int64_t mult;
};

// ---------------------------------------------------------------------------
// Failure-mode taxonomy (spec "FAILURE CONTRACT")
// ---------------------------------------------------------------------------

enum class FatalReason {
  TripleDiameter,           // triple-point cluster diameter > eps
  SubResolutionChain,       // skipped-contact cluster diameter > eps
  CoplanarOverlap,          // coplanar triangle overlap (out of scope)
  EdgeInPlane,              // edge lying in another tri's plane (out of scope)
  ClassificationAmbiguity,  // no located piece, or pieces disagree
  AnchorConflict,           // degenerate component with conflicting anchors
  UnclassifiableComponent,  // degenerate component diameter > eps, no anchors
  Starvation,               // all slabs for a region are sub-eps wide
  BalanceViolation,         // seam edge has unbalanced kept-region counts
  PSLGInvalid,              // seams cross without shared vert ids
  EngineIdConflict,         // 2D engine conflict counter nonzero
};

// Non-fatal counter accumulator (one per pipeline run).
struct Overlap3Counters {
  int subEpsContactsDropped = 0;  // point-like skipped contacts
  int degenerateClassified = 0;   // degenerate components classified by anchor
  int epsFeaturesDropped = 0;     // eps-scale degenerate components dropped
  int clearanceSkips = 0;         // pieces skipped for clearance
  int engineIdConflicts = 0;      // total engine id conflict events
};

// Per-stage result carrying either the product or a fatal reason.
// Failure stops the pipeline; counters are accumulated throughout.
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
  static StageResult<T> Fatal(FatalReason reason, std::string detail = {}) {
    StageResult<T> r;
    r.fatal = reason;
    r.detail = std::move(detail);
    return r;
  }
};

// ---------------------------------------------------------------------------
// Stage-C types
// ---------------------------------------------------------------------------

// One section segment for a straddling face in a slab. Directed: the vector
// p1-p0 has positive dot product with yz(cross(+x, face.normal)), per the
// orientation spec (section "ORIENTATION AND SIGN CONVENTIONS"). Used by stage
// E to determine m_lex sign relative to the seeded direction.
struct SectionFaceSegment {
  int faceId;
  vec2 p0, p1;   // section (y,z) directed per above convention
  int64_t mult;  // signed multiplicity (from CanonicalFace)
};

// Per-slab result from stage C.
struct SlabResult {
  double xLo, xHi, xMid;
  bool built;  // false = skipped (sub-eps slab)
  std::vector<SweepCapture> pieces;
  std::vector<SectionFaceSegment> segments;  // all straddling faces this slab
  // Test-hook: raw section edges and verts fed to the 2D engine (before
  // arrangement). Populated only when running via RemoveOverlaps3D_TestHooks.
  std::vector<EdgeM> sectionEdges;
  std::vector<vec2> sectionVerts;
};

// ---------------------------------------------------------------------------
// Stage-B internal types (seam geometry)
// ---------------------------------------------------------------------------

// A 3D vert in the merged vert pool (shared by stage A verts + event verts +
// triple-point verts).
struct MergedVert {
  vec3 pos;
};

// A seam between two canonical faces: a sequence of merged-vert ids from one
// face-face pair. Both faces share this polyline object; vert ids are globally
// canonical (from the unified triple-point pass).
struct Seam {
  int faceId0, faceId1;      // the two faces sharing this seam
  std::vector<int> vertIds;  // merged-vert ids, in segment order
};

// Per-face subdivided boundary edge: a list of vert ids spanning the original
// halfedge, including any event/seam/triple verts inserted in stage B.
struct SubdividedEdge {
  int faceId;
  int edgeIdx;  // halfedge index in the original mesh (for pairing)
  std::vector<int> vertIds;  // from start to end, including both endpoints
};

// ---------------------------------------------------------------------------
// Full pipeline input (after stage A+B)
// ---------------------------------------------------------------------------

struct ArrangementGeometry {
  std::vector<MergedVert> verts;          // all canonical 3D verts
  std::vector<CanonicalFace> faces;       // stage-A canonical faces
  std::vector<Seam> seams;                // face-face seam polylines
  std::vector<SubdividedEdge> subdEdges;  // per halfedge, with inserted verts
  // For each face: list of seam indices incident to that face.
  std::vector<std::vector<int>> faceSeams;
};

// ---------------------------------------------------------------------------
// Public seam
// ---------------------------------------------------------------------------

struct Overlap3Result {
  std::optional<Manifold::Impl> impl;
  Overlap3Counters counters;
  std::optional<FatalReason> fatal;
  std::string detail;
};

// Main entry point. `eps` = 0 -> compute from bounding box scale
// (EpsilonFromScale).
Overlap3Result RemoveOverlaps3D(const Manifold::Impl& in, double eps = 0.0);

// ---------------------------------------------------------------------------
// Test seams (overlap3_test.cpp only)
// ---------------------------------------------------------------------------

// Intermediate pipeline state exposed for gate-2 section validity checks and
// gate-1 seam parity assertions.
struct Overlap3Internals {
  ArrangementGeometry arr;  // after stages A+B
  std::vector<SlabResult>
      slabs;  // after stage C (sectionEdges/sectionVerts populated)
  std::optional<FatalReason> fatal;
  std::string detail;
  Overlap3Counters counters;
};

// Run stages A+B+C and return internals. Slabs include
// sectionEdges/sectionVerts.
Overlap3Internals RemoveOverlaps3D_TestHooks(const Manifold::Impl& in,
                                             double eps = 0.0);

// A region in the face PSLG, exposed for white-box balance and classification
// pins.  Stage D builds these internally; they are exposed here so test code
// can synthesize inputs for CheckSeamBalance_Test and ClassifyRegion_Test.
struct PSLGRegion {
  std::vector<int> loopVerts;               // boundary vert ids (CCW)
  std::vector<std::vector<int>> holeVerts;  // hole loops
  int64_t below = 0, above = 0;             // status-order windings
  bool classified = false;
  bool degenerate = false;  // no covering slab > eps
};

// Return value from ClassifyRegion_Test.
struct ClassifyRegionResult {
  std::optional<FatalReason> fatal;
  int64_t below = 0, above = 0;
  bool classified = false;
};

// White-box: call the seam-balance check on synthetic region data.  Returns
// the FatalReason if balance is violated, nullopt if the check passes.
std::optional<FatalReason> CheckSeamBalance_Test(
    const ArrangementGeometry& arr,
    const std::vector<std::vector<PSLGRegion>>& faceRegions);

// White-box: call ClassifyRegion on a single synthetic PSLGRegion with the
// given slab data.
ClassifyRegionResult ClassifyRegion_Test(const PSLGRegion& region, int faceId,
                                         const std::vector<SlabResult>& slabs,
                                         const std::vector<MergedVert>& verts,
                                         const CanonicalFace& face, double eps);

// Run stages C+D+E starting from a pre-built ArrangementGeometry.  Used by
// PSLGInvalid pin (P5) to inject synthetic bad-seam data.
Overlap3Result RemoveOverlaps3D_FromArr(const ArrangementGeometry& arr,
                                        double eps);

// White-box M6: run anchor-component propagation on synthetic region data.
// Returns nullopt on success (faceRegions updated), or a FatalReason on
// failure (AnchorConflict, UnclassifiableComponent).
std::optional<FatalReason> PropagateAnchorComponents_Test(
    const ArrangementGeometry& arr,
    std::vector<std::vector<PSLGRegion>>& faceRegions, double eps,
    Overlap3Counters& cnt);

}  // namespace manifold
