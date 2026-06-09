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

// Public-to-src interface for the overlap-removal pipeline that backs
// Manifold::RemoveSelfIntersections() (declared in
// include/manifold/manifold.h, implemented in src/manifold.cpp via
// RunOverlapRemoval below).
//
// Internal types + helpers — the per-stage data structures, the BVH
// query helpers, the chord-pair-sym phases, the cap walker, etc. —
// live in src/overlap_removal_internal.h. That header is for
// overlap_removal.cpp's own use and is not part of any installed
// public API.

#include <utility>  // for std::pair

#include "manifold/manifold.h"  // for Manifold

namespace manifold {
namespace overlap_removal {

// Per-stage diagnostic counters returned alongside the output mesh.
// All fields zero-initialized; populated by RunOverlapRemoval as the
// pipeline executes. Step numbering tracks Emmett Lalish's #289
// 13-step sketch (see docs/Overlap3D.md).
struct RemoveSelfIntersectionsStats {
  int mergedVerts;               // step 1: ε-merged vert pairs
  int vertsOnEdges;              // step 4: total verts on edge interiors
  int edgeEdgeIntersections;     // step 3: skew-edge crossings
  int vertsInsideTris;           // step 5: total verts in tri interiors
  int edgeTriIntersections;      // step 6: edge-pierces-tri events
  int chordVerts;                // step 7p2: new vert positions added
  int chordEdges;                // step 7p2: chord edges generated
  int vertsPropagatedToChords;   // step 8: in-tri verts on chord segments
  int newEdgeIntersections;      // step 9: chord-chord crossings
  int totalSubdividedHalfedges;  // step 10: per-tri halfedge count
};

// Top-level entry point composing all 13 #289 pipeline stages plus
// pair-symmetric chord enforcement, pierce-aware cap walker,
// pre/post-cap pierce reducers, and a pierce/drift gate with
// sign-flip recovery.
//
// Backs Manifold::RemoveSelfIntersections() (= one and only caller
// in production source). Returns {output, stats}.
//
// `eps` is optional: defaults to AlphaBudgetEpsilon(input.bbox_scale)
// when 0 (= the production pipeline's choice). Callers that already
// know their tolerance can pass it explicitly.
std::pair<Manifold, RemoveSelfIntersectionsStats> RunOverlapRemoval(
    const Manifold& input, double eps = 0.0);

}  // namespace overlap_removal
}  // namespace manifold
