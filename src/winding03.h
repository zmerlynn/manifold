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
// SoS-clean kernels lifted from the anonymous namespace of `boolean3.cpp`
// so non-Boolean callers can reuse them: `Shadow01` / `Kernel02` (#289
// step 5), `Kernel11` / `Kernel12` (#289 steps 4 / 6), `Winding03` (#289
// step 13), and the `Intersect` helper they share. Companion exposures
// of `Shadows`, `Interpolate`, `withSign` live in `shared.h`. Overlap-
// removal (`src/overlap_removal.cpp`) uses Kernel12 for SoS-clean edge-
// face piercing; the 3D spike (`extras/overlap3d_proto.cpp`) uses
// `Winding03_` to label per-vert winding inside a self-intersecting
// input mesh. See issues #1640, #289, #1445.
// =============================================================================

#include <unordered_set>

#include "collider.h"
#include "disjoint_sets.h"
#include "impl.h"
#include "parallel.h"
#include "shared.h"
#include "vec.h"

#if (MANIFOLD_PAR == 1)
#include <tbb/combinable.h>
#endif

namespace manifold {

// Edge-vert (0D-1D) symbolic-overlap test. Reports whether vert `a0` of
// inA shadows halfedge `b1` of inB along the y axis, with FP-stable
// `Interpolate` for the (y, z) of the projected intersection.
// Symbolic-perturbation tiebreaker via `Shadows` (in shared.h).
template <bool expandP, bool forward>
inline std::pair<int, vec2> Shadow01(const int a0, const int b1,
                                     const Manifold::Impl& inA,
                                     const Manifold::Impl& inB) {
  const int b1s = inB.halfedge_.Start(b1);
  const int b1e = inB.halfedge_.End(b1);
  const double a0x = inA.vertPos_[a0].x;
  const double b1sx = inB.vertPos_[b1s].x;
  const double b1ex = inB.vertPos_[b1e].x;
  const double a0xp = inA.vertNormal_[a0].x;
  const double b1sxp = inB.vertNormal_[b1s].x;
  const double b1exp = inB.vertNormal_[b1e].x;
  int s01 = forward ? Shadows(a0x, b1ex, withSign(expandP, a0xp) - b1exp) -
                          Shadows(a0x, b1sx, withSign(expandP, a0xp) - b1sxp)
                    : Shadows(b1sx, a0x, withSign(expandP, b1sxp) - a0xp) -
                          Shadows(b1ex, a0x, withSign(expandP, b1exp) - a0xp);
  vec2 yz01(NAN);

  if (s01 != 0) {
    yz01 =
        Interpolate(inB.vertPos_[b1s], inB.vertPos_[b1e], inA.vertPos_[a0].x);
    const int b1pair = inB.halfedge_.Pair(b1);
    const double dir =
        inB.faceNormal_[b1 / 3].y + inB.faceNormal_[b1pair / 3].y;
    if (forward) {
      if (!Shadows(inA.vertPos_[a0].y, yz01[0], -dir)) s01 = 0;
    } else {
      if (!Shadows(yz01[0], inA.vertPos_[a0].y, withSign(expandP, dir)))
        s01 = 0;
    }
  }
  return std::make_pair(s01, yz01);
}

// Vert-in-face (0D-2D) classification kernel. Tells whether a vertex
// of one input lies on the plane of the other input's triangle, with
// symbolic resolution of the on-edge / on-vert degeneracies via
// `Shadow01`. Emmett's #289 step 5.
template <bool expandP, bool forward>
struct Kernel02 {
  const Manifold::Impl& inA;
  const Manifold::Impl& inB;

  std::pair<int, double> operator()(int a0, int b2) {
    int s02 = 0;
    double z02 = 0.0;

    // For yzzLR[k], k==0 is the left and k==1 is the right.
    int k = 0;
    vec3 yzzRL[2];
    // Either the left or right must shadow, but not both. This ensures the
    // intersection is between the left and right.
    bool shadows = false;

    for (const int i : {0, 1, 2}) {
      const int b1 = 3 * b2 + i;
      const Halfedge edgeB = inB.halfedge_.Get(b1);
      const int b1F = edgeB.IsForward() ? b1 : edgeB.pairedHalfedge;

      const auto syz01 = Shadow01<expandP, forward>(a0, b1F, inA, inB);
      const int s01 = syz01.first;
      const vec2 yz01 = syz01.second;
      // If the value is NaN, then these do not overlap.
      if (std::isfinite(yz01[0])) {
        s02 += s01 * (forward == edgeB.IsForward() ? -1 : 1);
        if (k < 2 && (k == 0 || (s01 != 0) != shadows)) {
          shadows = s01 != 0;
          yzzRL[k++] = vec3(yz01[0], yz01[1], yz01[1]);
        }
      }
    }

    if (s02 == 0) {  // No intersection
      z02 = NAN;
    } else {
      DEBUG_ASSERT(k == 2, logicErr, "Boolean manifold error: s02");
      vec3 vertPosA = inA.vertPos_[a0];
      z02 = Interpolate(yzzRL[0], yzzRL[1], vertPosA.y)[1];
      if (forward) {
        if (!Shadows(vertPosA.z, z02, -inB.faceNormal_[b2].z)) s02 = 0;
      } else {
        if (!Shadows(z02, vertPosA.z, withSign(expandP, inB.faceNormal_[b2].z)))
          s02 = 0;
      }
    }
    return std::make_pair(s02, z02);
  }
};

// `Intersect` (used by Kernel11/Kernel12 below) lives in shared.h
// (upstream #1708 lifted it for Boolean2 reuse).

// Edge-edge (1D-1D) intersection kernel. Smith ch.8 §8.1; Emmett's
// #289 step 4. Returns (sign, xyzz4) where xyzz4 is the 4D-embedded
// intersection point and sign is 0/+1/-1 (no-cross / direction).
// Symbolic perturbation via `Shadow01`/`Shadows` ensures the cross-
// or-not decision is decided non-ambiguously even when the four edge
// endpoints are coplanar in the projection.
template <bool expandP>
struct Kernel11 {
  const Manifold::Impl& inP;
  const Manifold::Impl& inQ;

  std::pair<int, vec4> operator()(int p1, int q1) {
    vec4 xyzz11 = vec4(NAN);
    int s11 = 0;

    int k = 0;
    vec3 pRL[2], qRL[2];
    bool shadows = false;
    s11 = 0;

    const int p0[2] = {inP.halfedge_.Start(p1), inP.halfedge_.End(p1)};
    for (int i : {0, 1}) {
      const auto [s01, yz01] = Shadow01<expandP, true>(p0[i], q1, inP, inQ);
      if (std::isfinite(yz01[0])) {
        s11 += s01 * (i == 0 ? -1 : 1);
        if (k < 2 && (k == 0 || (s01 != 0) != shadows)) {
          shadows = s01 != 0;
          pRL[k] = inP.vertPos_[p0[i]];
          qRL[k] = vec3(pRL[k].x, yz01.x, yz01.y);
          ++k;
        }
      }
    }

    const int q0[2] = {inQ.halfedge_.Start(q1), inQ.halfedge_.End(q1)};
    for (int i : {0, 1}) {
      const auto [s10, yz10] = Shadow01<expandP, false>(q0[i], p1, inQ, inP);
      if (std::isfinite(yz10[0])) {
        s11 += s10 * (i == 0 ? -1 : 1);
        if (k < 2 && (k == 0 || (s10 != 0) != shadows)) {
          shadows = s10 != 0;
          qRL[k] = inQ.vertPos_[q0[i]];
          pRL[k] = vec3(qRL[k].x, yz10.x, yz10.y);
          ++k;
        }
      }
    }

    if (s11 == 0) {
      xyzz11 = vec4(NAN);
    } else {
      DEBUG_ASSERT(k == 2, logicErr, "Boolean manifold error: s11");
      xyzz11 = Intersect(pRL[0], pRL[1], qRL[0], qRL[1]);
      const int p1pair = inP.halfedge_.Pair(p1);
      const double dirP =
          inP.faceNormal_[p1 / 3].z + inP.faceNormal_[p1pair / 3].z;
      const int q1pair = inQ.halfedge_.Pair(q1);
      const double dirQ =
          inQ.faceNormal_[q1 / 3].z + inQ.faceNormal_[q1pair / 3].z;
      if (!Shadows(xyzz11.z, xyzz11.w, withSign(expandP, dirP) - dirQ)) s11 = 0;
    }
    return std::make_pair(s11, xyzz11);
  }
};

// Edge-face (1D-2D) intersection kernel — the main producer of new
// vertices in the 3D pipeline. Composes `Kernel02` (vert-in-face) and
// `Kernel11` (edge-edge) to handle the geometric and degenerate cases
// uniformly. Emmett's #289 step 6.
template <bool expandP, bool forward>
struct Kernel12 {
  const Manifold::Impl& inA;
  const Manifold::Impl& inB;
  Kernel02<expandP, forward> k02;
  Kernel11<expandP> k11;

  std::pair<int, vec3> operator()(int a1, int b2) {
    int x12 = 0;
    vec3 v12 = vec3(NAN);

    int k = 0;
    vec3 xzyLR0[2];
    vec3 xzyLR1[2];
    bool shadows = false;
    x12 = 0;

    const Halfedge edgeA = inA.halfedge_.Get(a1);

    for (int vertA : {edgeA.startVert, edgeA.endVert}) {
      const auto [s, z] = k02(vertA, b2);
      if (std::isfinite(z)) {
        x12 += s * ((vertA == edgeA.startVert) == forward ? 1 : -1);
        if (k < 2 && (k == 0 || (s != 0) != shadows)) {
          shadows = s != 0;
          xzyLR0[k] = inA.vertPos_[vertA];
          std::swap(xzyLR0[k].y, xzyLR0[k].z);
          xzyLR1[k] = xzyLR0[k];
          xzyLR1[k][1] = z;
          k++;
        }
      }
    }

    for (const int i : {0, 1, 2}) {
      const int b1 = 3 * b2 + i;
      const Halfedge edgeB = inB.halfedge_.Get(b1);
      const int b1F = edgeB.IsForward() ? b1 : edgeB.pairedHalfedge;
      const auto [s, xyzz] = forward ? k11(a1, b1F) : k11(b1F, a1);
      if (std::isfinite(xyzz[0])) {
        x12 -= s * (edgeB.IsForward() ? 1 : -1);
        if (k < 2 && (k == 0 || (s != 0) != shadows)) {
          shadows = s != 0;
          xzyLR0[k][0] = xyzz.x;
          xzyLR0[k][1] = xyzz.z;
          xzyLR0[k][2] = xyzz.y;
          xzyLR1[k] = xzyLR0[k];
          xzyLR1[k][1] = xyzz.w;
          if (!forward) std::swap(xzyLR0[k][1], xzyLR1[k][1]);
          k++;
        }
      }
    }

    if (x12 == 0) {
      v12 = vec3(NAN);
    } else {
      DEBUG_ASSERT(k == 2, logicErr, "Boolean manifold error: v12");
      const vec4 xzyy = Intersect(xzyLR0[0], xzyLR0[1], xzyLR1[0], xzyLR1[1]);
      v12.x = xzyy[0];
      v12.y = xzyy[2];
      v12.z = xzyy[1];
    }
    return std::make_pair(x12, v12);
  }
};

// Per-vert winding-number classification (the "isInside" decision for
// each face of inP relative to inQ, or vice versa). Standard
// arrangement-and-classify pattern (de Berg ch.8). Uses `DisjointSets`
// to propagate inclusion across uninterrupted halfedge runs, so each
// component only needs one explicit ray-cast / classification result.
// Emmett's #289 step 13. The 2D analog is `FilterByWindingDCEL` in
// extras/overlap2d_proto.cpp.
template <bool expandP, bool forward>
Vec<int> Winding03_(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                    const VecView<std::array<int, 2>> p1q2) {
  ZoneScoped;
  // a: 0 (vert), b: 2 (face)
  const Manifold::Impl& a = forward ? inP : inQ;
  const Manifold::Impl& b = forward ? inQ : inP;
  Vec<int> brokenHalfedges;
  int index = forward ? 0 : 1;

  DisjointSets uA(a.vertPos_.size());
  for_each(autoPolicy(a.halfedge_.size()), countAt(0),
           countAt(a.halfedge_.size()), [&](int edge) {
             const Halfedge he = a.halfedge_.Get(edge);
             if (!he.IsForward()) return;
             // check if the edge is broken
             auto it = std::lower_bound(
                 p1q2.begin(), p1q2.end(), edge,
                 [index](const std::array<int, 2>& collisionPair, int e) {
                   return collisionPair[index] < e;
                 });
             if (it == p1q2.end() || (*it)[index] != edge)
               uA.unite(he.startVert, he.endVert);
           });

  // find components, the hope is the number of components should be small
  std::unordered_set<int> components;
#if (MANIFOLD_PAR == 1)
  if (a.vertPos_.size() > 1e5) {
    tbb::combinable<std::unordered_set<int>> componentsShared;
    for_each(autoPolicy(a.vertPos_.size()), countAt(0),
             countAt(a.vertPos_.size()),
             [&](int v) { componentsShared.local().insert(uA.find(v)); });
    componentsShared.combine_each([&](const std::unordered_set<int>& data) {
      components.insert(data.begin(), data.end());
    });
  } else
#endif
  {
    for (size_t v = 0; v < a.vertPos_.size(); v++)
      components.insert(uA.find(v));
  }
  Vec<int> verts;
  verts.reserve(components.size());
  for (int c : components) verts.push_back(c);

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
  b.collider_.Collisions<false>(recorder, f, verts.size());
  // flood fill
  for_each(autoPolicy(w03.size()), countAt(0), countAt(w03.size()),
           [&](size_t i) {
             size_t root = uA.find(i);
             if (root == i) return;
             w03[i] = w03[root];
           });
  return w03;
}

template <bool forward>
Vec<int> Winding03(const Manifold::Impl& inP, const Manifold::Impl& inQ,
                   const VecView<std::array<int, 2>> p1q2, bool expandP) {
  if (expandP)
    return Winding03_<true, forward>(inP, inQ, p1q2);
  else
    return Winding03_<false, forward>(inP, inQ, p1q2);
}

}  // namespace manifold
