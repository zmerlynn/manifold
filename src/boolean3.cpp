// Copyright 2021 The Manifold Authors.
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

#include "boolean3.h"

#include <limits>

#include "parallel.h"
#include "shared.h"
#include "winding03.h"

#if (MANIFOLD_PAR == 1)
#include <tbb/combinable.h>
#endif

// =============================================================================
// 3D boolean implementation: Smith framework mapping.
//
// This file is the BVH-on-Smith 3D analog of the prototype in
// extras/overlap2d_proto.cpp. It implements the kernel cascade that
// finds geometric intersections (Kernel11 = edge-edge, Kernel02 =
// vert-in-face, Kernel12 = edge-in-face) and the winding-classification
// step (Winding03). See docs/Overlap3D.md for the mapping between
// Emmett Lalish's 13-step #289 sketch and the actual entry points
// here. Brief vocabulary table:
//
//   `Shadows` (in shared.h)            symbolic-perturbation orientation
//                                      predicate (Edelsbrunner-Mucke SoS
//                                      via withSign-flagged signs).
//   `Interpolate` (in shared.h)        FP-stable axis-overlap interpolation
//                                      with smaller-|dy| basepoint, used
//                                      by Shadow01 / Kernel11 / etc.
//   `Shadow01` / `Kernel11`            edge-edge (1D-1D) intersection;
//                                      Smith chapter 8 §8.1.
//   `Kernel02` / `Shadow02`            vert-in-face (0D-2D) classification.
//   `Kernel12` / `Intersect12`         edge-face (1D-2D) intersection;
//                                      the main producer of new verts.
//   `AddNewEdgeVerts` (boolean_result) eager propagation of new
//                                      intersection verts to all
//                                      relevant halfedge lists. 2D
//                                      analog: extras/overlap2d_proto.cpp
//                                      step 4 propagation phase.
//   `Winding03`                        per-vert winding-number
//                                      classification (the "isInside"
//                                      decision per face). Standard
//                                      arrangement-classify pattern.
//   `CollapseShortEdges` (edge_op)     post-boolean cleanup of near-
//                                      duplicate verts. 2D analog:
//                                      extras/overlap2d_proto.cpp
//                                      step 4b structural re-merge.
//
// Eps story: Manifold::Impl::epsilon_ and tolerance_ track per-instance
// position-precision bounds. Smith's alpha-budget formula (eps = (k+1) *
// sqrt(153) * u * L) lives in shared.h as `EpsilonFromScale`; this file
// uses whatever epsilon the Impl already carries. See docs/Overlap3D.md.
// =============================================================================

using namespace manifold;

namespace {

// `withSign`, `Interpolate`, `Intersect`, and `Shadows` are in shared.h
// (upstream #1708 lifted them for Boolean2 reuse).
// `Shadow01`, `Kernel02`, `Kernel11`, `Kernel12`, `Winding03_`, `Winding03`
// are in winding03.h (branch lift so overlap_removal can reuse them; see
// issues #1640, #289, #1445).

template <bool expandP, bool forward>
struct Kernel12Recorder {
  using Local = Intersections;
  Kernel12<expandP, forward>& k12;

#if MANIFOLD_PAR == 1
  tbb::combinable<Intersections> store;
  Local& local() { return store.local(); }
#else
  Intersections localStore;
  Local& local() { return localStore; }
#endif

  void record(int queryIdx, int leafIdx, Local& tmp) {
    const auto [x12, v12] = k12(queryIdx, leafIdx);
    if (std::isfinite(v12[0])) {
      if (forward)
        tmp.p1q2.push_back({queryIdx, leafIdx});
      else
        tmp.p1q2.push_back({leafIdx, queryIdx});
      tmp.x12.push_back(x12);
      tmp.v12.push_back(v12);
    }
  }

  Intersections get() {
#if MANIFOLD_PAR == 1
    Intersections result;
    std::vector<Intersections> tmps;
    store.combine_each(
        [&](Intersections& data) { tmps.emplace_back(std::move(data)); });
    std::vector<size_t> sizes;
    size_t total_size = 0;
    for (const auto& tmp : tmps) {
      sizes.push_back(total_size);
      total_size += tmp.x12.size();
    }
    result.p1q2.resize(total_size);
    result.x12.resize(total_size);
    result.v12.resize(total_size);
    for_each_n(ExecutionPolicy::Seq, countAt(0), tmps.size(), [&](size_t i) {
      std::copy(tmps[i].p1q2.begin(), tmps[i].p1q2.end(),
                result.p1q2.begin() + sizes[i]);
      std::copy(tmps[i].x12.begin(), tmps[i].x12.end(),
                result.x12.begin() + sizes[i]);
      std::copy(tmps[i].v12.begin(), tmps[i].v12.end(),
                result.v12.begin() + sizes[i]);
    });
    return result;
#else
    return localStore;
#endif
  }
};

// Run all edge-face intersection queries between inP's edges and
// inQ's faces (or vice versa, when forward = false), via the BVH
// `Collider` broad phase. Returns the per-pair intersection results
// permuted into edge-major order so downstream `AddNewEdgeVerts` can
// scan contiguously per edge. This is the BVH-based 3D analog of the
// 2D prototype's step 4 pair iteration.
template <bool expandP, bool forward>
Intersections Intersect12_(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                           ExecutionContext::Impl* ctx) {
  ZoneScoped;
  // Invariant: every ctx-passing parallel op is followed by IsCancelled to
  // keep partial output from feeding unconditional downstream consumers.
  // a: 1 (edge), b: 2 (face)
  const Manifold::Impl& a = forward ? inP : inQ;
  const Manifold::Impl& b = forward ? inQ : inP;

  Kernel02<expandP, forward> k02{a, b};
  Kernel11<expandP> k11{inP, inQ};

  Kernel12<expandP, forward> k12{a, b, k02, k11};
  Kernel12Recorder<expandP, forward> recorder{k12, {}};
  auto f = [&a](int i) {
    const int start = a.halfedge_.Start(i);
    const int end = a.halfedge_.End(i);
    return start < end ? Box(a.vertPos_[start], a.vertPos_[end]) : Box();
  };
  b.collider_.Collisions<false>(recorder, f, a.halfedge_.size(), true, ctx);
  if (IsCancelled(ctx)) return Intersections{};

  Intersections result = recorder.get();
  auto& p1q2 = result.p1q2;
  // sort p1q2 according to edges
  Vec<size_t> i12(p1q2.size());
  sequence(i12.begin(), i12.end());

  int index = forward ? 0 : 1;
  stable_sort(i12.begin(), i12.end(), [&](auto a, auto b) {
    return p1q2[a][index] < p1q2[b][index] ||
           (p1q2[a][index] == p1q2[b][index] &&
            p1q2[a][1 - index] < p1q2[b][1 - index]);
  });
  Permute(p1q2, i12);
  Permute(result.x12, i12);
  Permute(result.v12, i12);
  return result;
};

template <bool forward>
Intersections Intersect12(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                          bool expandP, ExecutionContext::Impl* ctx) {
  if (expandP)
    return Intersect12_<true, forward>(inP, inQ, ctx);
  else
    return Intersect12_<false, forward>(inP, inQ, ctx);
}

// boolean3.cpp keeps its own Winding03_/Winding03 here (vs the
// branch's winding03.h copy) because this version threads
// ExecutionContext for mid-Boolean cancellation. The winding03.h
// copy is the one overlap_removal / self_mesh_analysis use.
template <bool expandP, bool forward>
Vec<int> Winding03_(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                    const VecView<std::array<int, 2>> p1q2,
                    ExecutionContext::Impl* ctx) {
  ZoneScoped;
  // a: 0 (vert), b: 2 (face)
  const Manifold::Impl& a = forward ? inP : inQ;
  const Manifold::Impl& b = forward ? inQ : inP;
  Vec<int> brokenHalfedges;
  int index = forward ? 0 : 1;

  // Invariant: every ctx-passing parallel op is followed by IsCancelled to
  // keep partial output from feeding unconditional downstream consumers.
  DisjointSets uA(a.vertPos_.size());
  for_each(autoPolicy(a.halfedge_.size()), countAt(0),
           countAt(static_cast<int>(a.halfedge_.size())), ctx, [&](int edge) {
             const int start = a.halfedge_.Start(edge);
             const int end = a.halfedge_.End(edge);
             if (start >= end) return;
             // check if the edge is broken
             auto it = std::lower_bound(
                 p1q2.begin(), p1q2.end(), edge,
                 [index](const std::array<int, 2>& collisionPair, int e) {
                   return collisionPair[index] < e;
                 });
             if (it == p1q2.end() || (*it)[index] != edge) uA.unite(start, end);
           });
  if (IsCancelled(ctx)) return Vec<int>{};

  // find components, the hope is the number of components should be small
  std::unordered_set<size_t> components;
#if (MANIFOLD_PAR == 1)
  if (a.vertPos_.size() > 1e5) {
    tbb::combinable<std::unordered_set<size_t>> componentsShared;
    for_each(autoPolicy(a.vertPos_.size()), countAt(0_uz),
             countAt(a.vertPos_.size()), ctx,
             [&](size_t v) { componentsShared.local().insert(uA.find(v)); });
    componentsShared.combine_each([&](const std::unordered_set<size_t>& data) {
      components.insert(data.begin(), data.end());
    });
  } else
#endif
  {
    for (size_t v = 0; v < a.vertPos_.size(); v++)
      components.insert(uA.find(v));
  }
  if (IsCancelled(ctx)) return Vec<int>{};
  Vec<int> verts;
  verts.reserve(components.size());
  for (size_t c : components) verts.push_back(static_cast<int>(c));

  Vec<int> w03(a.NumVert(), 0);
  Kernel02<expandP, forward> k02{a, b};
  auto recorderf = [&](int i, int b) {
    const auto [s02, z02] = k02(verts[i], b);
    // note that i is distinct on each thread, and verts contains unique
    // elements, so this does not require atomics
    if (std::isfinite(z02)) w03[verts[i]] += s02 * (forward ? 1 : -1);
  };
  auto recorder = MakeSimpleRecorder(recorderf);
  auto f = [&](int i) { return a.vertPos_[verts[i]]; };
  b.collider_.Collisions<false>(recorder, f, verts.size(), true, ctx);
  if (IsCancelled(ctx)) return Vec<int>{};
  // flood fill
  for_each(autoPolicy(w03.size()), countAt(0_uz), countAt(w03.size()), ctx,
           [&](size_t i) {
             size_t root = uA.find(i);
             if (root == i) return;
             w03[i] = w03[root];
           });
  if (IsCancelled(ctx)) return Vec<int>{};
  return w03;
}

template <bool forward>
Vec<int> Winding03(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                   const VecView<std::array<int, 2>> p1q2, bool expandP,
                   ExecutionContext::Impl* ctx) {
  if (expandP)
    return Winding03_<true, forward>(inP, inQ, p1q2, ctx);
  else
    return Winding03_<false, forward>(inP, inQ, p1q2, ctx);
}
}  // namespace

namespace manifold {
Boolean3::Boolean3(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                   OpType op, ExecutionContext::Impl* ctx)
    : inP_(inP), inQ_(inQ), expandP_(op == OpType::Add), ctx_(ctx) {
  ZoneScoped;
  // Symbolic perturbation:
  // Union -> expand inP, expand inQ
  // Difference, Intersection -> contract inP, expand inQ
  // Technically Intersection should contract inQ, but doing it this way makes
  // Split faster and any suboptimal cases seem pretty rare.

  constexpr size_t INT_MAX_SZ =
      static_cast<size_t>(std::numeric_limits<int>::max());

  if (inP.IsEmpty() || inQ.IsEmpty() || !inP.bBox_.DoesOverlap(inQ.bBox_)) {
    PRINT("No overlap, early out");
    w03_.resize(inP.NumVert(), 0);
    w30_.resize(inQ.NumVert(), 0);
    return;
  }

  // Phase-boundary fast-path: skip launching the next phase if cancel fired
  // between phases. The per-phase invariant is enforced inside the called
  // functions (Intersect12_, Winding03_).
#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  Timer intersections;
  intersections.Start();
  Timer intersect12P, intersect12Q, winding03P, winding03Q;
  intersect12P.Start();
#endif

  // Level 3
  // Build up the intersection of the edges and triangles, keeping only those
  // that intersect, and record the direction the edge is passing through the
  // triangle.
  if (IsCancelled(ctx_)) return;
  xv12_ = Intersect12<true>(inP, inQ, expandP_, ctx_);
#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  intersect12P.Stop();
  intersect12Q.Start();
#endif
  if (IsCancelled(ctx_)) return;
  xv21_ = Intersect12<false>(inP, inQ, expandP_, ctx_);
#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  intersect12Q.Stop();
#endif

  if (xv12_.x12.size() > INT_MAX_SZ || xv21_.x12.size() > INT_MAX_SZ) {
    valid = false;
    return;
  }

  // Compute winding numbers of all vertices using flood fill
  // Vertices on the same connected component have the same winding number
#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  winding03P.Start();
#endif
  if (IsCancelled(ctx_)) return;
  w03_ = Winding03<true>(inP, inQ, xv12_.p1q2, expandP_, ctx_);
#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  winding03P.Stop();
  winding03Q.Start();
#endif
  if (IsCancelled(ctx_)) return;
  w30_ = Winding03<false>(inP, inQ, xv21_.p1q2, expandP_, ctx_);
  // No trailing check: Winding03_ already returns empty on cancel and
  // Boolean3::Result re-checks on entry.

#if defined(MANIFOLD_DEBUG) || defined(MANIFOLD_TIMING)
  winding03Q.Stop();
  intersections.Stop();

  if (ManifoldParams().verbose >= 2) {
    intersect12P.Print("  Intersect12 P->Q");
    intersect12Q.Print("  Intersect12 Q->P");
    winding03P.Print("  Winding03 P");
    winding03Q.Print("  Winding03 Q");
    intersections.Print("Intersections (total)");
  }
#endif
}
Vec<int> Manifold::Impl::PointWinding(VecView<const vec3> points) const {
  ZoneScoped;
  Vec<int> winding(points.size(), 0);
  if (points.empty() || IsEmpty()) return winding;

  // Points outside the bounding box have winding 0 for any closed manifold.
  Vec<int> query2input;
  query2input.reserve(points.size());
  for (size_t i = 0; i < points.size(); ++i)
    if (bBox_.Contains(points[i])) query2input.push_back(static_cast<int>(i));
  if (query2input.empty()) return winding;

  // Build a minimal Impl for the query points. Kernel02 only reads vertPos_
  // and vertNormal_ from inA; halfedges and faces are not accessed.
  Impl pointImpl;
  pointImpl.vertPos_.resize(query2input.size());
  pointImpl.vertNormal_.resize(query2input.size(), vec3(0.0));
  for (size_t i = 0; i < query2input.size(); ++i)
    pointImpl.vertPos_[i] = points[query2input[i]];

  // expandP=false: no symbolic perturbation from the query side (zero normals).
  // forward=true: project along +Z to count signed face crossings above each
  // point. f returns vec3 → collider uses DoesOverlap(vec3), which is
  // XY-projected, so all faces with XY overlap are returned regardless of Z
  // (correct for +Z winding).
  Kernel02<false, true> k02{pointImpl, *this};
  auto recorderf = [&](int localIdx, int tri) {
    const auto [s02, z02] = k02(localIdx, tri);
    if (std::isfinite(z02)) winding[query2input[localIdx]] += s02;
  };
  auto recorder = MakeSimpleRecorder(recorderf);
  auto f = [&pointImpl](int i) { return pointImpl.vertPos_[i]; };
  // Each queryIdx is processed by at most one thread (for_each_n), so the
  // per-queryIdx accumulation into winding[] is race-free without atomics.
  collider_.Collisions<false>(recorder, f, static_cast<int>(query2input.size()),
                              true);
  return winding;
}

std::vector<RayHit> Manifold::Impl::RayCast(vec3 origin, vec3 endpoint) const {
  ZoneScoped;
  if (IsEmpty()) return {};
  const vec3 dir = endpoint - origin;
  if (la::dot(dir, dir) == 0.0) return {};

  // Build a minimal single-edge Impl representing the ray segment.
  // Kernel12 treats inA as an edge mesh. Halfedges derives endVert from the
  // next edge in a triangle, so this helper uses a padded degenerate face:
  // edge 0 is the forward ray, edge 1 is its reverse, and edge 2 only closes
  // the local storage loop. Zero vertex normals and face normal mean the ray
  // contributes nothing to perturbation tiebreakers; consistency at shared
  // edges/vertices depends entirely on the mesh's own normals.
  Impl rayImpl;
  rayImpl.vertPos_.resize(2);
  rayImpl.vertPos_[0] = origin;
  rayImpl.vertPos_[1] = endpoint;
  rayImpl.vertNormal_.resize(2);
  rayImpl.vertNormal_[0] = vec3(0.0);
  rayImpl.vertNormal_[1] = vec3(0.0);
  rayImpl.halfedge_.resize(3);
  rayImpl.halfedge_.Set(0, 0, 1, 0);  // forward: vert 0 → 1
  rayImpl.halfedge_.Set(1, 1, 0, 0);  // backward: vert 1 → 0
  rayImpl.halfedge_.Set(2, -1, -1, 0);
  rayImpl.faceNormal_.resize(1);
  rayImpl.faceNormal_[0] = vec3(0.0);

  // expandP=false with zero vertNormal means the ray-side perturbation is
  // zero. forward=true means we project along +Z for the lower-dimensional
  // kernel cascade (Shadow01 → Kernel02 → Kernel11 → Kernel12).
  Kernel02<false, true> k02{rayImpl, *this};
  Kernel11<false> k11{rayImpl, *this};
  Kernel12<false, true> k12{rayImpl, *this, k02, k11};

  // Use the component with largest magnitude for stable t computation.
  const vec3 absDir = la::abs(dir);
  const int tAxis = absDir.x > absDir.y && absDir.x > absDir.z ? 0
                    : absDir.y > absDir.z                      ? 1
                                                               : 2;

  std::vector<RayHit> hits;
  // Query the BVH with the ray's AABB.
  const Box rayBox(la::min(origin, endpoint), la::max(origin, endpoint));
  auto recorderf = [&](int /*queryIdx*/, int tri) {
    const auto [s, v] = k12(0, tri);  // halfedge 0 vs triangle tri
    if (s != 0 && std::isfinite(v.x)) {
      // v is the 3D intersection point computed by Kernel12.
      // Compute parametric t ∈ [0,1] along the ray segment.
      const double t = (v[tAxis] - origin[tAxis]) / dir[tAxis];
      if (t >= 0.0 && t <= 1.0) {
        hits.push_back({static_cast<uint64_t>(tri), t, v, faceNormal_[tri]});
      }
    }
  };
  auto recorder = MakeSimpleRecorder(recorderf);
  auto f = [&rayBox](int) { return rayBox; };
  collider_.Collisions<false>(recorder, f, 1, false);

  std::sort(hits.begin(), hits.end(), [](const RayHit& a, const RayHit& b) {
    return a.distance < b.distance;
  });
  return hits;
}

}  // namespace manifold
