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

// 3D overlap removal as a per-component regularization operator.
// Design: docs/Regularize3D.md.  Maps a valid oriented (possibly
// self-overlapping) face soup to the boundary of {p : w_S(p) >= 1}, per
// connected component, never fusing separate components (fusion is the
// Boolean's job).  Internal intersection segments only; public wiring is
// post-prototype.

#pragma once

#include <climits>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "boolean2.h"  // EpsilonFromScale, EdgeM, RemoveOverlaps2D
#include "impl.h"
#include "manifold/common.h"

namespace manifold {

// ---------------------------------------------------------------------------
// Fatal-reason taxonomy.
// ---------------------------------------------------------------------------

enum class FatalReason {
  SubEpsInput,          // eps <= 0 or degenerate input geometry
  NonManifoldEmission,  // emitted triangulation is not 2-manifold
  // Regularization operator (docs/Regularize3D.md): regularization RAN on a
  // component that needed it but could not complete exactly - a non-coplanar
  // exact-zero SoS residue, a >2-sheet
  // triple point, a coplanar/transversal entanglement, a near-coplanar
  // planarity-guard failure, or a filter-uncertain winding probe.  The honest
  // fail-closed (recorded reason, no output), never a silent wrong result.
  RegularizationIncomplete,
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
// Regularization operator (docs/Regularize3D.md) - public entry point.
//
// RemoveOverlaps3D maps a valid oriented face soup to the boundary of the solid
// {p : w_S(p) >= 1}, PER CONNECTED COMPONENT (it never fuses separate
// components - fusion is the Boolean's job, already done upstream).
// Cross-component overlap of any kind - touching, coplanar, or transversal - is
// out of scope: a buried plug with a coincident cap is two components, each
// regularized on its own and concatenated back unchanged (docs/Regularize3D.md
// non-fusion contract).  A coplanar self-overlap is resolved only when it is
// INTERNAL to one connected component, where component classification detects
// it and routes it to regularization.  The pipeline: DECOMPOSE by connectivity
// -> CLASSIFY (validity + IsSelfIntersecting + within-component coplanar
// overlap) -> EARLY-EXIT unchanged components -> regularize the remaining
// components -> VALIDATE the output -> COMPOSE BACK by concatenation.
// ---------------------------------------------------------------------------

// White-box dispatch counters (the dispatch pins read these directly).
struct RegularizeCounters {
  int components = 0;           // components after decompose by connectivity
  int unchanged = 0;            // unchanged by regularization
  int needsRegularization = 0;  // self-intersecting OR coplanar overlap
  int regularized = 0;          // regularization produced valid output
  int failClosed =
      0;  // components that fail closed (regularization incomplete, output
          // validation failure, or unexpected non-manifold input)
};

struct RegularizeResult {
  std::optional<Manifold::Impl> impl;  // composed output; absent on any fatal
  std::optional<FatalReason> fatal;
  std::string detail;
  RegularizeCounters counters;
};

// eps = 0 -> compute from bounding-box scale.
RegularizeResult RemoveOverlaps3D(const Manifold::Impl& in, double eps = 0.0);

// ---------------------------------------------------------------------------
// Test hooks (overlap3_test.cpp only).
// ---------------------------------------------------------------------------

// The regularization mechanism probe (docs/Regularize3D.md
// "the regularizer's
// mechanism"), exposed so the port of the validated fragment (enumeration +
// coupled winding) is tested directly against the fragment's recorded numbers.
// intersectionSegmentCount = genuine non-adjacent self-crossings;
// boundaryTouchPairs = pairs whose deciding predicate hit an exact-zero /
// static-filter-uncertain boundary (the single-global-SoS case);
// probeWinding[i] = coupled soup winding w_S at probes[i] cast to `seed`
// (kWindingUncertain if a deciding predicate was filter-uncertain).
struct ComponentEnumProbe {
  int intersectionSegmentCount = 0;
  int boundaryTouchPairs = 0;
  int coplanarClusterFaces = 0;  // faces in an exact-coplanar overlap cluster
  std::vector<int> probeWinding;
};
constexpr int kWindingUncertain = INT_MIN;
ComponentEnumProbe EnumerateComponent_Probe(
    const Manifold::Impl& inputComponent, const std::vector<vec3>& probes,
    const vec3& seed);

// Test hook: run regularization (coplanar arrangement + emission + output
// validation) directly on a soup treated as ONE input component, bypassing
// decomposition and classification. The self-intersection check does not flag a
// pure coplanar overlap, so this hook is the way to exercise the exact-coplanar
// arrangement's path on an isolated coplanar cluster (no transversal
// entanglement).
RegularizeResult RegularizeComponentDirect(const Manifold::Impl& soup,
                                           double eps);

// Test hook: the exact-kernel tie test behind the stage-6 SoS - the EXACT
// orient3d sign (0 iff the four points are exactly coplanar).  The w==1
// instantiation of the ONE homogeneous predicate form
// (sos::HomogOrient3DSign<true>), adaptive-width two's-complement accumulator,
// TOTAL for every finite-double input (no window-fail, no expansion fallback).
// Exposed so the property pin can grade it directly (filter agreement,
// antisymmetry, constructed exact zeros, scaling invariance).
int Orient3DExactSignProbe(const vec3& a, const vec3& b, const vec3& c,
                           const vec3& d);

}  // namespace manifold
