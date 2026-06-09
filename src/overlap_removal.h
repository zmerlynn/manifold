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
// Internal types + helpers - the per-stage data structures, the BVH
// query helpers, the chord-pair-sym phases, the cap walker, etc. -
// live in src/overlap_removal_internal.h. That header is for
// overlap_removal.cpp's own use and is not part of any installed
// public API.

#include "manifold/manifold.h"  // for Manifold

namespace manifold {
namespace overlap_removal {

// Top-level entry point composing all 13 #289 pipeline stages plus
// pair-symmetric chord enforcement, pierce-aware cap walker,
// pre/post-cap pierce reducers, and a pierce/drift gate with
// sign-flip recovery.
//
// Backs Manifold::RemoveSelfIntersections() (= one and only caller in
// production source). Returns the cleaned manifold, or the input
// (unchanged or eps-merged) on the fallback paths.
//
// `eps` is optional: defaults to AlphaBudgetEpsilon(input.bbox_scale)
// when 0 (= the production pipeline's choice). Callers that already
// know their tolerance can pass it explicitly.
Manifold RunOverlapRemoval(const Manifold& input, double eps = 0.0);

}  // namespace overlap_removal
}  // namespace manifold
