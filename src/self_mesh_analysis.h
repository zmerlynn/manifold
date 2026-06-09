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

// =============================================================================
// SelfMeshAnalysis: per-vert above/below winding for a single self-
// intersecting input Manifold::Impl.
//
// Two-sided winding analysis for the polygon-keep classifier in
// src/overlap_removal.cpp. Reuses the production Winding03_ flood-fill
// structure (DisjointSets + component decomposition) but substitutes a
// literal geometric eps-offset ray-cast for the per-component
// classification step. Reasons:
//
//   1. Boolean3 itself can't be called with M==M - its SoS perturbation
//      collapses at shared verts (= why `m.Boolean(m, OpType::Add)`
//      doesn't work as a self-union shortcut).
//   2. The SoS-perturbed Winding03_ in boolean3.cpp gives one
//      z-projection winding number per vert; the `expandP` template
//      parameter controls a tiebreaker, not a geometric "above" vs
//      "below" probe direction. Confirmed empirically on the .obj
//      battery: Winding03_<true> and Winding03_<false> produce nearly
//      identical results for self-mesh.
//   3. A polygon-keep classifier needs both above and below windings
//      to decide "is this polygon on the boundary of the (winding >= 1)
//      region of M?" The clean answer is to probe at v +/- eps * n(v) and
//      ray-cast.
//
// Flood-fill validity: verts connected by intact halfedges share the
// same outward direction near the surface (intact halfedges preserve
// normal orientation between adjacent faces). One ray-cast per
// component labels every vert in that component.
// =============================================================================

#include <array>

#include "impl.h"
#include "vec.h"

namespace manifold {

struct SelfMeshAnalysis {
  // Per-vert winding number of M evaluated at v +/- eps * n(v). Length =
  // M.NumVert(). For a clean (non-self-intersecting) closed manifold,
  // w_above=1 / w_below=0 (or the reverse) for every vert. For
  // self-intersecting M, both may be >= 1 inside overlap regions; the
  // surface separates regions where |w_above - w_below| = 1 from the
  // exterior (w=0).
  Vec<int> w_above;
  Vec<int> w_below;
};

// Run the self-mesh winding-classify pass on M, given a pre-computed
// edge-face pierce list (= the "broken halfedge" set).
//
// p1q2[i] = {forward halfedge index in M, pierced face index in M},
// sorted by halfedge index ascending - same shape Winding03_ expects
// for its flood-fill component decomposition.
//
// Per-component flood-fill (DisjointSets) + per-component representative
// ray-cast for the above and below windings.  No mutation of M.
SelfMeshAnalysis AnalyzeSelfMesh(const Manifold::Impl& M,
                                 VecView<const std::array<int, 2>> p1q2);

// Compute the signed winding number of M at a single 3D point `origin`
// by casting a ray in `direction` (must be unit length) of length
// `length` (should comfortably exceed the mesh diameter), summing
// signed face crossings via Moller-Trumbore. Convention: "inside = +1",
// "outside = 0" (= matches the standard outward-normal-oriented closed
// manifold convention).
//
// For the polygon-keep classifier in overlap_removal: probe at
// (poly_centroid +/- eps * n_T) with ray direction picked to avoid
// axis-grazing, and read both windings to decide if the polygon is
// on the boundary of M's winding >= 1 region.
int WindingAt(const Manifold::Impl& M, vec3 origin, vec3 direction,
              double length);

// Shared probe scale: Box::Scale() (absolute-largest coordinate, the
// manifold convention) when bBox_ is finite, else the vertPos_ diagonal;
// 0 for degenerate input. Used by both AnalyzeSelfMesh and
// overlap_removal's polygon-keep classifier so the per-vert and
// per-polygon probes derive identical eps and ray length.
double ProbeMeshScale(const Manifold::Impl& M);

// Probe-direction + scaling constants shared between AnalyzeSelfMesh
// (per-vert eps-offset ray-cast) and overlap_removal's polygon-keep
// classifier (per-polygon centroid eps-offset ray-cast). Keeping a
// single source-of-truth so the two probes always agree.
//
// kProbeRayDir: irrational unit-vector chosen to minimize the chance
//   of axis-grazing through the mesh's BVH (= each component is a
//   distinct non-rational fraction; pre-normalized to unit length to
//   avoid an FP normalize() call in the hot loop).
// kProbeOffsetCoeff: eps-offset distance from probe origin to ray
//   start, expressed as a fraction of the mesh's bbox scale. Small
//   enough that the offset doesn't cross any mesh feature for normal
//   inputs; large enough to clear FP noise around exact-on-surface
//   verts.
// kProbeRayLengthCoeff: ray length, as a fraction of bbox scale.
//   Must comfortably exceed the mesh diameter (= 4x bbox is
//   conservative since mesh diameter <= sqrt(3)*bbox edge ~= 1.73x).
// kProbeRayDir is `inline const` (not constexpr) because la::normalize
// is not constexpr in linalg.h. Computed once at startup, no hot-path
// cost. The direction is arbitrary but off-axis to avoid axis-aligned
// ray degeneracies.
inline const vec3 kProbeRayDir = la::normalize(vec3(0.7234, 0.4567, 0.5191));
constexpr double kProbeOffsetCoeff = 1e-9;
constexpr double kProbeRayLengthCoeff = 4.0;

}  // namespace manifold
