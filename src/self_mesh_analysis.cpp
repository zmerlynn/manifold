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

#include "self_mesh_analysis.h"

#include <algorithm>
#include <cmath>
#include <set>

#include "disjoint_sets.h"
#include "impl.h"
#include "parallel.h"
#include "vec.h"

// =============================================================================
// AnalyzeSelfMesh implementation.
//
// Algorithm:
//   1. DisjointSets flood-fill over M's verts. Halfedges in p1q2 are
//      "broken" (their edge is pierced by a non-incident face). Verts
//      connected by intact halfedges land in the same component;
//      intact halfedges preserve outward-normal orientation, so the
//      "outward" region near each component's surface is well-defined.
//   2. Per component, pick a representative vert v_rep. Compute the
//      winding number of M at v_rep + eps * n(v_rep) (= "above") and
//      v_rep - eps * n(v_rep) (= "below") via Manifold::Impl::RayCast.
//   3. Propagate: for any vert v in component c, w_above[v] =
//      w_above[c_rep] and similarly for below.
//
// Why eps-offset ray-cast instead of pure SoS? The production
// `Winding03_<expandP, true>(M, M, p1q2)` in boolean3.cpp
// computes a single z-projection winding via Kernel02 BVH-collide;
// the `expandP` SoS direction does not correspond to "above" vs
// "below" the surface (it perturbs the comparison's tiebreaker, not
// the geometric probe direction). For a polygon-keep classifier we
// need both sides, so we substitute a literal geometric eps offset
// along v_rep's outward normal.
//
// Ray-cast narrow phase: Manifold::Impl::RayCast (Kernel02/11/12 with
// symbolic perturbation), BVH broad phase via M.collider_. SoS
// resolves on-edge / on-vert grazes deterministically. Signed
// contribution per hit is dot(direction, hit.normal): exit-with-flow
// (s>0) is +1, enter-against-flow (s<0) is -1.
// =============================================================================

namespace {

int RayCastWindingShared(const manifold::Manifold::Impl& M,
                         const manifold::vec3& origin,
                         const manifold::vec3& direction, double length) {
  using manifold::la::dot;
  const manifold::vec3 endpoint = origin + direction * length;
  const auto hits = M.RayCast(origin, endpoint);
  int w = 0;
  for (const auto& hit : hits) {
    // Convention: "inside = +1." A ray traveling in `direction`
    // exits the solid through a face whose outward normal aligns
    // with `direction` (dot > 0); each such hit means we were inside
    // before that exit, so it contributes +1 to the winding count
    // of the origin.
    const double s = dot(direction, hit.normal);
    if (s > 0)
      w += 1;
    else if (s < 0)
      w -= 1;
  }
  return w;
}

}  // namespace

namespace manifold {

int WindingAt(const Manifold::Impl& M, vec3 origin, vec3 direction,
              double length) {
  return RayCastWindingShared(M, origin, direction, length);
}

double ProbeMeshScale(const Manifold::Impl& M) {
  // Box::Scale() (absolute-largest coordinate, the manifold convention)
  // when bBox_ is finite, else the diagonal from vertPos_; 0 for a
  // degenerate or single-point input.
  if (M.bBox_.IsFinite()) return M.bBox_.Scale();
  Box bb;
  for (size_t v = 0; v < M.NumVert(); ++v) bb.Union(M.vertPos_[v]);
  return bb.IsFinite() ? bb.Scale() : 0.0;
}

SelfMeshAnalysis AnalyzeSelfMesh(const Manifold::Impl& M,
                                 VecView<const std::array<int, 2>> p1q2) {
  using manifold::la::length;
  using manifold::la::normalize;
  SelfMeshAnalysis r;
  const size_t nVert = M.NumVert();
  r.w_above.resize(nVert, 0);
  r.w_below.resize(nVert, 0);
  if (nVert == 0) return r;

  // A degenerate (single-point) input yields meshScale=0 -> eps=0 ->
  // all-zero windings (the trivial case handled below).
  const double meshScale = ProbeMeshScale(M);
  // 1. Flood-fill components via DisjointSets, breaking at p1q2 edges.
  DisjointSets uA(nVert);
  for_each(autoPolicy(M.halfedge_.size()), countAt(0),
           countAt(static_cast<int>(M.halfedge_.size())), [&](int edge) {
             const Halfedge he = M.halfedge_.Get(edge);
             if (!he.IsForward()) return;
             auto it = std::lower_bound(
                 p1q2.begin(), p1q2.end(), edge,
                 [](const std::array<int, 2>& cp, int e) { return cp[0] < e; });
             if (it == p1q2.end() || (*it)[0] != edge)
               uA.unite(he.startVert, he.endVert);
           });

  // 2. Collect representatives. std::set (sorted) keeps iteration
  // order deterministic across runs - we hit the unordered_set
  // determinism class once already (commit 5acc6ab5).
  //
  // Determinism fix (2026-05-10): DisjointSets::unite is thread-safe
  // (CAS-based) but the chosen REP per component depends on the
  // parallel union order - same connectivity, different rep across
  // runs. Since the rep's vertNormal_ becomes the eps-offset probe
  // direction, a different rep gives a different probe -> different
  // ray-cast winding -> different keep decisions. Symptom: 1/5 runs
  // on self-intersect produced k=1=645 instead of k=1=379 (with
  // pierce 2061 vs 312). Fix: re-anchor each component to its
  // smallest-id member (= deterministic regardless of union order).
  std::map<int, int> compMinId;  // raw uA.find(v) -> smallest v in component
  for (size_t v = 0; v < nVert; ++v) {
    const int c = static_cast<int>(uA.find(v));
    auto it = compMinId.find(c);
    if (it == compMinId.end())
      compMinId[c] = static_cast<int>(v);
    else if (static_cast<int>(v) < it->second)
      it->second = static_cast<int>(v);
  }
  // Build deterministic rep map: rawRep -> minId.
  // Use sorted output set to preserve deterministic iteration.
  std::set<int> compSet;
  for (const auto& [_, minId] : compMinId) compSet.insert(minId);
  Vec<int> reps;
  reps.reserve(compSet.size());
  for (int c : compSet) reps.push_back(c);

  // 3. Per-rep geometric eps-offset ray-cast (above + below).
  // Degenerate-input guard: if mesh has zero scale (all coincident
  // verts), skip the ray-casts and leave w_above/w_below at zero.
  if (meshScale == 0.0) return r;
  const double eps = meshScale * kProbeOffsetCoeff;
  const double rayLen = meshScale * kProbeRayLengthCoeff;
  const vec3& rayDir = kProbeRayDir;
  const bool hasVertNormals = (M.vertNormal_.size() == nVert);

  Vec<int> wa(reps.size(), 0), wb(reps.size(), 0);
  for_each(autoPolicy(reps.size()), countAt(0),
           countAt(static_cast<int>(reps.size())), [&](int i) {
             const int v = reps[i];
             // Defensive: degenerate or missing vert normal -> use
             // rayDir as offset. This still gives a valid winding
             // (the offset just isn't perfectly aligned with the
             // local outward direction); the per-component flood-fill
             // is still consistent because all verts in a component
             // share the same probe.
             vec3 offset = rayDir;
             if (hasVertNormals) {
               const vec3 n = M.vertNormal_[v];
               const double nLen = length(n);
               if (nLen > 0.0) offset = n / nLen;
             }
             const vec3 originAbove = M.vertPos_[v] + offset * eps;
             const vec3 originBelow = M.vertPos_[v] - offset * eps;
             wa[i] = RayCastWindingShared(M, originAbove, rayDir, rayLen);
             wb[i] = RayCastWindingShared(M, originBelow, rayDir, rayLen);
           });

  // 4. Propagate via DisjointSets: every vert inherits its rep's
  // winding. Both sides labelled "outward of the local surface
  // patch" via the rep's outward normal - well-defined per
  // component because intact halfedges preserve normal orientation.
  // Indexed by deterministic min-id rep (= via compMinId).
  Vec<int> repAbove(nVert, 0), repBelow(nVert, 0);
  for (size_t i = 0; i < reps.size(); ++i) {
    repAbove[reps[i]] = wa[i];
    repBelow[reps[i]] = wb[i];
  }
  for_each(autoPolicy(nVert), countAt(0_uz), countAt(nVert), [&](size_t v) {
    const int rawRep = static_cast<int>(uA.find(v));
    const int minRep = compMinId.at(rawRep);
    r.w_above[v] = repAbove[minRep];
    r.w_below[v] = repBelow[minRep];
  });
  return r;
}

}  // namespace manifold
