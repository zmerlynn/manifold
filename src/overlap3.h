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
  CoplanarOverlap,      // coplanar face overlap (out of scope)
  EdgeInPlane,          // edge lying in another face's plane (out of scope)
  SubEpsInput,          // eps <= 0 or degenerate input geometry
  SubEpsFeature,        // macro-scale face in merged sub-eps critical run
  EngineIdConflict,     // 2D engine source-id conflict in a slab
  NonManifoldEmission,  // emitted triangulation is not 2-manifold
};

// Non-fatal counter accumulator.
struct Overlap3Counters {
  int subEpsContactsDropped = 0;  // point-like skipped contacts
  int engineIdConflicts = 0;      // total engine id-conflict events
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

struct MergedVert {
  vec3 pos;
};

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
  std::vector<MergedVert> verts;     // all canonical 3D verts
  std::vector<CanonicalFace> faces;  // canonical faces
  std::vector<Seam> seams;           // face-pair seam segments
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

}  // namespace manifold
