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

// Src-internal interface for the overlap-removal pipeline that backs
// Manifold::RemoveSelfIntersections() (declared in
// include/manifold/manifold.h, implemented in src/manifold.cpp via
// RemoveOverlaps below).
//
// Internal types + the per-stage function seams live in
// src/overlap_removal_internal.h - a test-only header, not part of
// any installed public API; file-local kernels stay in
// overlap_removal.cpp's anonymous namespace.

#include <optional>

#include "execution_impl.h"  // ExecutionContext::Impl, IsCancelled
#include "impl.h"            // Manifold::Impl

namespace manifold {
namespace overlap_removal {

// Top-level entry point composing the 13 #289 pipeline stages
// (arrangement -> cell complex -> winding classification -> emit)
// behind a fail-closed gate. Impl-to-Impl, like Boolean3: no Manifold
// is held anywhere inside the pipeline.
//
// Backs Manifold::RemoveSelfIntersections() (= one and only caller in
// production source; the member owns status propagation and wrapper
// identity). Three outcomes:
//  - a rebuilt Impl on success;
//  - nullopt on the early-exit and EVERY fallback arm - the caller
//    returns its own input bit-identically (wrapper identity is a
//    Manifold concern, not a pipeline concern);
//  - an Impl made empty with Error::Cancelled when `ctx` reports
//    cancellation at a stage boundary (the ADVANCE_PHASE_OR_RETURN
//    idiom: cancellation is observable status, never a silent
//    fallback).
//
// `eps` is the WORKING EPSILON (the pipeline's computational scale,
// not the mesh tolerance - see the eps contract in
// docs/OverlapRemoval.md). eps <= 0 (the inferred-eps path) derives
// AlphaBudgetEpsilon(input.bBox_.Scale(), 1000) for the first attempt,
// then retries a small fixed ladder of wider multiples of that base
// before returning nullopt. eps > 0 is a single attempt at exactly
// that eps (no ladder). `ctx` may be null (no cancellation checks).
std::optional<Manifold::Impl> RemoveOverlaps(const Manifold::Impl& input,
                                             double eps,
                                             ExecutionContext::Impl* ctx);

}  // namespace overlap_removal
}  // namespace manifold
