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
// Internal types + the per-stage function seams live in
// src/overlap_removal_internal.h - a test-only header, not part of
// any installed public API; file-local kernels stay in
// overlap_removal.cpp's anonymous namespace.

#include "manifold/manifold.h"  // for Manifold

namespace manifold {
namespace overlap_removal {

// Top-level entry point composing the 13 #289 pipeline stages
// (arrangement -> cell complex -> winding classification -> emit)
// behind a fail-closed gate.
//
// Backs Manifold::RemoveSelfIntersections() (= one and only caller in
// production source). Returns the rebuilt manifold, or the input
// BIT-IDENTICALLY on the early-exit and every fallback path.
//
// `eps` is optional: defaults to
// AlphaBudgetEpsilon(input.BoundingBox().Scale(), 1000)
// when 0 (= the production pipeline's choice). Callers that already
// know their WORKING EPSILON can pass it explicitly (this is the
// pipeline's computational scale, not the mesh tolerance - see the
// eps contract in docs/OverlapRemoval.md).
Manifold RunOverlapRemoval(const Manifold& input, double eps = 0.0);

}  // namespace overlap_removal
}  // namespace manifold
