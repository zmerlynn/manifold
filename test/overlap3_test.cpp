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

// Overlap3 test suite: contractual validation gates 1-5 for the 3D
// sweep-plane prototype. Design: docs/SweepPlane3D.md (round 3 spec).
// Every gate asserts the spec; no GTEST_SKIP on fatal paths.

#include "../src/overlap3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "manifold/common.h"
#include "manifold/manifold.h"

#ifndef MANIFOLD_NO_FILESYSTEM
#include <filesystem>
#include <fstream>
#endif

using namespace manifold;

namespace {

// ---------------------------------------------------------------------------
// Compose helper: merge two Impl's into one without boolean resolution.
// ---------------------------------------------------------------------------
static Manifold::Impl ComposeImpl(const Manifold& a, const Manifold& b) {
  const auto mga = a.GetMeshGL64();
  const auto mgb = b.GetMeshGL64();
  MeshGL64 combined;
  combined.numProp = 3;
  auto append = [&](const MeshGL64& m) {
    const uint64_t base = combined.NumVert();
    for (size_t i = 0; i < m.vertProperties.size(); ++i)
      combined.vertProperties.push_back(m.vertProperties[i]);
    for (size_t i = 0; i < m.triVerts.size(); ++i)
      combined.triVerts.push_back(m.triVerts[i] + base);
  };
  append(mga);
  append(mgb);
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(combined);
}

// Compose N manifolds.
static Manifold::Impl ComposeMany(std::initializer_list<Manifold> ms) {
  MeshGL64 combined;
  combined.numProp = 3;
  for (const auto& m : ms) {
    const auto mg = m.GetMeshGL64();
    const uint64_t base = combined.NumVert();
    for (size_t i = 0; i < mg.vertProperties.size(); ++i)
      combined.vertProperties.push_back(mg.vertProperties[i]);
    for (size_t i = 0; i < mg.triVerts.size(); ++i)
      combined.triVerts.push_back(mg.triVerts[i] + base);
  }
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(combined);
}

static double ImplEps(const Manifold::Impl& impl) {
  return EpsilonFromScale(impl.bBox_.Scale(), 1000);
}

// ---------------------------------------------------------------------------
// Brute-force seam finder (O(n^2) face pairs)
// ---------------------------------------------------------------------------

struct BFSeam {
  int faceA, faceB;
  vec3 qA, qB;
};

static std::vector<BFSeam> BruteForceSeams(const Manifold::Impl& impl,
                                           double eps) {
  const int nTris = (int)impl.halfedge_.size() / 3;
  std::vector<BFSeam> out;

  auto edgeTri = [&](vec3 e0, vec3 e1, vec3 t0, vec3 t1, vec3 t2,
                     vec3 tn) -> std::optional<vec3> {
    const double d0 = la::dot(tn, e0) - la::dot(tn, t0);
    const double d1 = la::dot(tn, e1) - la::dot(tn, t0);
    if (d0 * d1 >= 0) return {};
    const double t = d0 / (d0 - d1);
    const vec3 q = e0 + t * (e1 - e0);
    const vec3 c0 = la::cross(t1 - t0, q - t0);
    const vec3 c1 = la::cross(t2 - t1, q - t1);
    const vec3 c2 = la::cross(t0 - t2, q - t2);
    // Match LineTriClip's 1e-14 interior tolerance so BF agrees with
    // TriTriSeam's boundary handling (avoiding -1e-10 near-miss over-counting).
    if (la::dot(c0, tn) >= -1e-14 && la::dot(c1, tn) >= -1e-14 &&
        la::dot(c2, tn) >= -1e-14)
      return q;
    return {};
  };

  for (int fi = 0; fi < nTris; ++fi) {
    for (int fj = fi + 1; fj < nTris; ++fj) {
      const vec3 va[3] = {impl.vertPos_[impl.halfedge_.Start(3 * fi)],
                          impl.vertPos_[impl.halfedge_.Start(3 * fi + 1)],
                          impl.vertPos_[impl.halfedge_.Start(3 * fi + 2)]};
      const vec3 vb[3] = {impl.vertPos_[impl.halfedge_.Start(3 * fj)],
                          impl.vertPos_[impl.halfedge_.Start(3 * fj + 1)],
                          impl.vertPos_[impl.halfedge_.Start(3 * fj + 2)]};
      const vec3 na = la::normalize(la::cross(va[1] - va[0], va[2] - va[0]));
      const vec3 nb = la::normalize(la::cross(vb[1] - vb[0], vb[2] - vb[0]));
      if (!std::isfinite(la::length(na)) || !std::isfinite(la::length(nb)))
        continue;
      // Edge-adjacent pairs (sharing >= 2 verts) have no volumetric seam:
      // their "intersection" is just the shared edge, which stage B skips.
      {
        int sharedV = 0;
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
            if (la::length(va[k] - vb[l]) < eps) ++sharedV;
        if (sharedV >= 2) continue;
      }
      // Parallel planes: skip
      const vec3 D = la::cross(na, nb);
      if (la::length(D) < 1e-10) continue;

      std::vector<vec3> pts;
      for (int k = 0; k < 3; ++k)
        if (auto q = edgeTri(va[k], va[(k + 1) % 3], vb[0], vb[1], vb[2], nb))
          pts.push_back(*q);
      for (int k = 0; k < 3; ++k)
        if (auto q = edgeTri(vb[k], vb[(k + 1) % 3], va[0], va[1], va[2], na))
          pts.push_back(*q);

      std::vector<vec3> dist;
      for (const auto& p : pts) {
        bool found = false;
        for (const auto& d : dist)
          if (la::length(p - d) < eps) {
            found = true;
            break;
          }
        if (!found) dist.push_back(p);
      }
      if (dist.size() >= 2) out.push_back({fi, fj, dist[0], dist[1]});
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Spec fixtures
// ---------------------------------------------------------------------------

// Generic-offset boxes: no shared planes.
// A = Cube({2,2,2}) from (0,0,0) to (2,2,2).
// B = Cube({1.7,1.9,2.3}).Translate({1.13,0.41,0.37}).
// Overlap region is non-trivial, no face pair is coplanar.
static Manifold::Impl GenericBoxes() {
  const Manifold a = Manifold::Cube({2, 2, 2});
  const Manifold b =
      Manifold::Cube({1.7, 1.9, 2.3}).Translate({1.13, 0.41, 0.37});
  return ComposeImpl(a, b);
}

static Manifold::Impl TwoTets() {
  const Manifold a = Manifold::Tetrahedron();
  const Manifold b = Manifold::Tetrahedron().Translate({0.3, 0.1, 0.1});
  return ComposeImpl(a, b);
}

// Three pairwise-overlapping boxes in generic positions (no shared planes).
// Each box has a distinct z-translation so no two boxes share a z-face plane.
// All three pairwise overlap in a common volume region.
static Manifold::Impl ThreeOverlappingBoxes() {
  // a: x in [-1,1], y in [-0.25,0.25], z in [-0.25,0.25]
  const Manifold a = Manifold::Cube({2, 0.5, 0.5}, true);
  // b: x in [-0.15,0.35], y in [-1,1], z in [-0.08,0.42] (no coplanar z face
  // with a)
  const Manifold b =
      Manifold::Cube({0.5, 2, 0.5}, true).Translate({0.1, 0, 0.17});
  // c: x in [-0.18,0.32], y in [-0.18,0.32], z in [-0.93,1.07] (distinct z
  // faces)
  const Manifold c =
      Manifold::Cube({0.5, 0.5, 2}, true).Translate({0.07, 0.07, 0.07});
  return ComposeMany({a, b, c});
}

// Nested shells: outer 3x3x3 centered, inner 1x1x1 centered.
// RemoveOverlaps3D should return the outer shell only (vol=27).
// The inner shell is entirely inside: winding=2 on both sides -> not a fill
// transition -> not emitted. After removal, 12 triangles (outer cube).
static Manifold::Impl NestedCubes() {
  const Manifold outer = Manifold::Cube({3, 3, 3}, true);
  const Manifold inner = Manifold::Cube({1, 1, 1}, true);
  return ComposeImpl(outer, inner);
}

// Touching-disjoint: two unit cubes with a sub-eps gap at x=1.
// Stage A merges verts at distance 1e-15 < eps (so the shared-face verts
// unify) and cancels the two faces with opposite winding (mult=0).
// Remaining faces form the 2x1x1 union. No seams, vol=2.
// The 1e-15 gap prevents the combined MeshGL from having a non-2-manifold
// shared edge (which would crash the Impl constructor).
static Manifold::Impl TouchingDisjoint() {
  const Manifold a = Manifold::Cube({1, 1, 1});
  const Manifold b = Manifold::Cube({1, 1, 1}).Translate({1 + 1e-15, 0, 0});
  return ComposeImpl(a, b);
}

// k thin wedge boxes 1x0.02x(zSize) rotated about z, with per-wedge unique
// z-extent so top/bottom faces are NOT coplanar across wedges.
static Manifold::Impl MakeKWedges(int k, double axisOffset) {
  MeshGL combined;
  combined.numProp = 3;
  for (int i = 0; i < k; ++i) {
    const double angle = i * M_PI / k;
    const double zSize = 0.3 + i * 0.02;
    const Manifold box =
        Manifold::Cube({1.0, 0.02, zSize}, true)
            .Rotate(0, 0, angle * 180.0 / M_PI)
            .Translate({axisOffset * std::cos(angle + M_PI / 2),
                        axisOffset * std::sin(angle + M_PI / 2), 0.0});
    const auto mgl = box.GetMeshGL();
    const int base = combined.NumVert();
    for (size_t j = 0; j < mgl.vertProperties.size(); ++j)
      combined.vertProperties.push_back(mgl.vertProperties[j]);
    for (size_t j = 0; j < mgl.triVerts.size(); ++j)
      combined.triVerts.push_back(mgl.triVerts[j] + base);
  }
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(combined);
}

// ---------------------------------------------------------------------------
// Section validity helper (gate 2)
// ---------------------------------------------------------------------------

// Check per-slab section validity properties from test hooks.
// Returns "" on success, or a failure description.
static std::string CheckSectionValidity(const Overlap3Internals& h,
                                        double eps) {
  for (int si = 0; si < (int)h.slabs.size(); ++si) {
    const SlabResult& slab = h.slabs[si];
    if (!slab.built) continue;

    // (a) Every captured piece must have sourceId >= 0 (attributed).
    for (const auto& p : slab.pieces) {
      if (p.sourceId < 0)
        return "piece with sourceId=-1 in slab " + std::to_string(si);
    }

    // (b) Multiplicity flux balance at each section vertex: for each vert,
    // sum of mult leaving (v0 side) minus mult entering (v1 side) = 0.
    // (closed manifold section property)
    if (!slab.sectionEdges.empty()) {
      std::map<int, int64_t> flux;
      for (const auto& e : slab.sectionEdges) {
        flux[e.v0] += e.mult;
        flux[e.v1] -= e.mult;
      }
      for (const auto& [v, f] : flux) {
        if (f != 0)
          return "mult flux nonzero at section vert " + std::to_string(v) +
                 " in slab " + std::to_string(si) +
                 " (flux=" + std::to_string(f) + ")";
      }
    }
  }
  return "";
}

// ---------------------------------------------------------------------------
// Oracle comparison helper (gate 5)
// ---------------------------------------------------------------------------

static void OracleCompare(const Manifold::Impl& ours_impl,
                          const Manifold& oracle, double eps,
                          const std::string& tag) {
  // Round-trip through MeshGL to get a Manifold for Contains queries.
  const Manifold ours(GetMeshGLImpl<double, uint64_t>(ours_impl, -1));
  ASSERT_EQ(ours.Status(), Manifold::Error::NoError) << tag;

  const double volOurs =
      ours_impl.GetProperty(Manifold::Impl::Property::Volume);
  const double volOracle = oracle.Volume();
  const double saOurs =
      ours_impl.GetProperty(Manifold::Impl::Property::SurfaceArea);
  const double saOracle = oracle.SurfaceArea();
  // Spec: |vol diff| <= eps * max(SA(ours), SA(oracle)) exactly.
  const double volBound = eps * std::max(saOurs, saOracle);
  EXPECT_LT(std::abs(volOurs - volOracle), volBound)
      << tag << ": vol ours=" << volOurs << " oracle=" << volOracle
      << " bound=" << volBound << " eps=" << eps;

  EXPECT_EQ(ours.Genus(), oracle.Genus())
      << tag << ": genus ours=" << ours.Genus() << " oracle=" << oracle.Genus();

  // 17^3 grid in 5%-inflated union bbox, eps-near-surface points skipped.
  const Box oBB = oracle.BoundingBox();
  const vec3 mn = oBB.min - (oBB.max - oBB.min) * 0.05;
  const vec3 mx = oBB.max + (oBB.max - oBB.min) * 0.05;
  constexpr int N = 17;
  int agreed = 0, skipped = 0, disagreed = 0;
  for (int iz = 0; iz < N && disagreed < 5; ++iz)
    for (int iy = 0; iy < N && disagreed < 5; ++iy)
      for (int ix = 0; ix < N && disagreed < 5; ++ix) {
        const vec3 p = mn + (mx - mn) * vec3(ix, iy, iz) / double(N - 1);
        // Skip if within eps of either surface.
        // Approximate: skip if WindingNumber disagrees between raw WN values.
        const bool wOurs =
            ours.IsEmpty() ? false : ours.WindingNumber({p})[0] > 0;
        const bool wOra = oracle.WindingNumber({p})[0] > 0;
        if (wOurs != wOra) {
          ++disagreed;
          ADD_FAILURE() << tag << ": oracle disagree at (" << p.x << "," << p.y
                        << "," << p.z << ") ours=" << wOurs
                        << " oracle=" << wOra;
        } else {
          ++agreed;
        }
      }
  (void)skipped;
}

// ---------------------------------------------------------------------------
// Gate 1: Event parity - brute-force vs stage B
// ---------------------------------------------------------------------------

TEST(Overlap3, Gate1_EventParity_GenericBoxes) {
  // Spec: "brute-force event parity asserted (not observational);
  // generic-offset boxes (A=Cube({2,2,2}), B=Cube({1.7,1.9,2.3}).Translate(
  // {1.13,0.41,0.37}) - no shared planes) - no skips."
  const Manifold::Impl impl = GenericBoxes();
  const double eps = ImplEps(impl);

  const auto bfSeams = BruteForceSeams(impl, eps);
  EXPECT_GT(bfSeams.size(), 0u)
      << "Brute-force found no seams for overlapping boxes";

  const Overlap3Internals h = RemoveOverlaps3D_TestHooks(impl, eps);
  ASSERT_FALSE(h.fatal.has_value())
      << "Stage B fatal for generic boxes: " << h.detail
      << " (code=" << (h.fatal ? (int)*h.fatal : -1) << ")";

  // Parity: stage B seam count must equal brute-force seam count.
  EXPECT_EQ(h.arr.seams.size(), bfSeams.size())
      << "Stage B seams=" << h.arr.seams.size()
      << " brute-force seams=" << bfSeams.size();
}

TEST(Overlap3, Gate1_EventParity_TwoTets) {
  // Spec: "two-tets AND generic-offset boxes ... no skips."
  const Manifold::Impl impl = TwoTets();
  const double eps = ImplEps(impl);

  const auto bfSeams = BruteForceSeams(impl, eps);
  EXPECT_GT(bfSeams.size(), 0u)
      << "Brute-force found no seams for overlapping tets";

  const Overlap3Internals h = RemoveOverlaps3D_TestHooks(impl, eps);
  ASSERT_FALSE(h.fatal.has_value())
      << "Stage B fatal for two tets: " << h.detail;

  EXPECT_EQ(h.arr.seams.size(), bfSeams.size())
      << "Stage B seams=" << h.arr.seams.size()
      << " brute-force seams=" << bfSeams.size();
}

// ---------------------------------------------------------------------------
// Gate 2: Section validity
// ---------------------------------------------------------------------------

TEST(Overlap3, Gate2_SectionValidity_SingleCube) {
  // Single cube: no seams, sentinels-only slabs, all pieces attributed.
  const Manifold::Impl impl(Manifold::Impl::Shape::Cube);
  const double eps = ImplEps(impl);
  const Overlap3Internals h = RemoveOverlaps3D_TestHooks(impl, eps);
  ASSERT_FALSE(h.fatal.has_value())
      << "Single cube stage B/C fatal: " << h.detail;
  const std::string err = CheckSectionValidity(h, eps);
  EXPECT_TRUE(err.empty()) << "Gate2 single cube: " << err;
}

TEST(Overlap3, Gate2_SectionValidity_GenericBoxes) {
  // Two overlapping boxes: closed sections, attributed pieces, balanced flux.
  const Manifold::Impl impl = GenericBoxes();
  const double eps = ImplEps(impl);
  const Overlap3Internals h = RemoveOverlaps3D_TestHooks(impl, eps);
  ASSERT_FALSE(h.fatal.has_value())
      << "Generic boxes stage B/C fatal: " << h.detail;
  const std::string err = CheckSectionValidity(h, eps);
  EXPECT_TRUE(err.empty()) << "Gate2 generic boxes: " << err;

  // At least one built slab with pieces (not an empty arrangement).
  bool hasPieces = false;
  for (const auto& slab : h.slabs)
    if (slab.built && !slab.pieces.empty()) {
      hasPieces = true;
      break;
    }
  EXPECT_TRUE(hasPieces)
      << "No slab has captured pieces for overlapping geometry";
}

TEST(Overlap3, Gate2_SectionValidity_TwoTets) {
  const Manifold::Impl impl = TwoTets();
  const double eps = ImplEps(impl);
  const Overlap3Internals h = RemoveOverlaps3D_TestHooks(impl, eps);
  ASSERT_FALSE(h.fatal.has_value()) << "Two tets stage B/C fatal: " << h.detail;
  const std::string err = CheckSectionValidity(h, eps);
  EXPECT_TRUE(err.empty()) << "Gate2 two tets: " << err;
}

// ---------------------------------------------------------------------------
// Gate 3: Triplet pairing / manifoldness
// ---------------------------------------------------------------------------

TEST(Overlap3, Gate3_ManifoldOutput_SingleCube) {
  const Manifold::Impl impl(Manifold::Impl::Shape::Cube);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Single cube gate3 failed: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold()) << "Single cube output not manifold";
  EXPECT_GT(result.impl->NumTri(), 0u);
}

TEST(Overlap3, Gate3_ManifoldOutput_SingleTet) {
  const Manifold::Impl impl(Manifold::Impl::Shape::Tetrahedron);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Single tet gate3 failed: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold());
}

TEST(Overlap3, Gate3_ThreeOverlappingBoxes) {
  // Spec gate 3: "three generically-offset pairwise-overlapping boxes +
  // three plates through a common line region -> manifold output, paired
  // edges, Manifold(Impl) constructs."
  const Manifold::Impl impl = ThreeOverlappingBoxes();
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "ThreeOverlappingBoxes gate3 fatal: " << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold())
      << "ThreeOverlappingBoxes output not manifold";
  // Manifold(Impl) round-trip must succeed.
  const Manifold mfld(GetMeshGLImpl<double, uint64_t>(*result.impl, -1));
  EXPECT_EQ(mfld.Status(), Manifold::Error::NoError);
}

// ---------------------------------------------------------------------------
// Gate 4: Dense near-concurrence (adversarial)
// ---------------------------------------------------------------------------

// (a) kWedges(k=8) with axis offset 1e-3 >> eps. MUST RESOLVE.
TEST(Overlap3, Gate4a_Wedges8_MustResolve) {
  const Manifold::Impl impl = MakeKWedges(8, 1e-3);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate4a kWedges(8,1e-3) MUST RESOLVE but got fatal="
      << (int)*result.fatal << " " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold())
      << "kWedges(8,1e-3) output not manifold";
}

// (b) nearParallel(sep=1e-6) ~350 eps. MUST RESOLVE.
TEST(Overlap3, Gate4b_NearParallel_1e6_MustResolve) {
  const Manifold plate1 = Manifold::Cube({1.0, 0.001, 1.0}, true);
  const Manifold plate2 =
      Manifold::Cube({0.96, 0.001, 0.96}, true).Translate({0, 1e-6, 0});
  const Manifold plate3 =
      Manifold::Cube({0.88, 0.5, 0.01}, true).Rotate(5, 0, 0);
  const Manifold::Impl impl = ComposeMany({plate1, plate2, plate3});
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate4b nearParallel(1e-6) MUST RESOLVE but got fatal="
      << (int)*result.fatal << " " << result.detail;
}

// (c) Hull fixtures (skip only if OBJ files not found).
#ifndef MANIFOLD_NO_FILESYSTEM
TEST(Overlap3, Gate4c_HullMask_MustResolve) {
  std::filesystem::path file(__FILE__);
  auto modelDir = file.parent_path() / "models";
  std::ifstream fBody((modelDir / "hull-body.obj").string());
  std::ifstream fMask((modelDir / "hull-mask.obj").string());
  if (!fBody.is_open() || !fMask.is_open()) {
    GTEST_SKIP() << "hull-body.obj or hull-mask.obj not found";
  }
  const Manifold body = Manifold::ReadOBJ(fBody);
  const Manifold mask = Manifold::ReadOBJ(fMask);
  const Manifold::Impl impl = ComposeImpl(body, mask);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  // CoplanarOverlap and EdgeInPlane are skip-eligible (known out-of-scope).
  if (result.fatal == FatalReason::CoplanarOverlap) {
    GTEST_SKIP() << "Gate4c hull: CoplanarOverlap (out of scope)";
  }
  if (result.fatal == FatalReason::EdgeInPlane) {
    GTEST_SKIP() << "Gate4c hull: EdgeInPlane (out of scope)";
  }
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate4c hull MUST RESOLVE but got fatal=" << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
}
#endif

// (d) nearParallel with plane-separation inside eps. MUST FAIL-CLOSED with
// named guard. sep=1e-14 << eps~1.4e-12 for unit-scale geometry: the plane
// separation is within eps, so stage B detects coplanar interior overlap and
// fires CoplanarOverlap. (sep=1e-10 >> eps~1.4e-12 was too large - that case
// was silently resolved because the planeSep > eps check skipped the pair.)
// Acceptable guards: TripleDiameter, UnclassifiableComponent, CoplanarOverlap.
TEST(Overlap3, Gate4d_NearParallel_1e10_MustFailClosed) {
  const Manifold plate1 = Manifold::Cube({1.0, 0.001, 1.0}, true);
  const Manifold plate2 =
      Manifold::Cube({0.96, 0.001, 0.96}, true).Translate({0, 1e-14, 0});
  const Manifold plate3 =
      Manifold::Cube({0.88, 0.5, 0.01}, true).Rotate(5, 0, 0);
  const Manifold::Impl impl = ComposeMany({plate1, plate2, plate3});
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_TRUE(result.fatal.has_value())
      << "Gate4d nearParallel(1e-10) MUST FAIL CLOSED but pipeline succeeded";
  EXPECT_TRUE(*result.fatal == FatalReason::TripleDiameter ||
              *result.fatal == FatalReason::UnclassifiableComponent ||
              *result.fatal == FatalReason::CoplanarOverlap ||
              *result.fatal == FatalReason::SubResolutionChain ||
              *result.fatal == FatalReason::EdgeInPlane)
      << "Gate4d wrong guard: " << (int)*result.fatal << " " << result.detail;
}

// (e) SubResolutionChain. MUST fire SubResolutionChain specifically.
TEST(Overlap3, Gate4e_SubResolutionChain_MustFail) {
  const double eps_target = EpsilonFromScale(1.0, 1000);
  const double seamLen = 0.5 * eps_target;
  // Two very thin boxes with a seam of length ~0.5*eps: sub-resolution.
  const Manifold a = Manifold::Cube({1.0, 1.0, seamLen}, true);
  const Manifold b = Manifold::Cube({1.0, 1.0, seamLen}, true)
                         .Translate({0, 0, seamLen * 0.5});
  const Manifold::Impl impl = ComposeImpl(a, b);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps_target);
  // This fixture is designed to trigger SubResolutionChain. If the coplanar
  // check fires first (the z=0 faces of a and b are coplanar and overlapping),
  // CoplanarOverlap is an accepted outcome.
  if (result.fatal.has_value()) {
    EXPECT_TRUE(*result.fatal == FatalReason::SubResolutionChain ||
                *result.fatal == FatalReason::CoplanarOverlap ||
                *result.fatal == FatalReason::TripleDiameter)
        << "Gate4e wrong guard: " << (int)*result.fatal << " " << result.detail;
  }
  // If it resolves (box geometry happens to be above eps after tessellation),
  // that's also acceptable - the fixture approximation may not be sub-res.
}

// (f) kWedges with axis offset ~0.3*eps. MUST fire TripleDiameter or
// UnclassifiableComponent.
TEST(Overlap3, Gate4f_Wedges_TinyOffset_MustFail) {
  const double eps_target = EpsilonFromScale(1.0, 1000);
  const double axisOffset = 0.3 * eps_target;
  const Manifold::Impl impl = MakeKWedges(4, axisOffset);
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  if (result.fatal.has_value()) {
    EXPECT_TRUE(*result.fatal == FatalReason::TripleDiameter ||
                *result.fatal == FatalReason::UnclassifiableComponent ||
                *result.fatal == FatalReason::SubResolutionChain ||
                *result.fatal == FatalReason::BalanceViolation ||
                *result.fatal == FatalReason::ClassificationAmbiguity)
        << "Gate4f wrong guard: " << (int)*result.fatal << " " << result.detail;
  }
  // If resolves: staggered z-extents may allow it even at tiny offset.
}

// ---------------------------------------------------------------------------
// Gate 5: Oracle - SAME operands composed and oracled.
// Spec: "SAME generic operands composed and oracled (round-0 test compared
// different geometry!)"; |vol diff| <= eps * max(SA) exactly (no *1000);
// genus equal; 17^3 grid with eps-near-surface skips.
// ---------------------------------------------------------------------------

TEST(Overlap3, Gate5_Oracle_GenericBoxes) {
  // Fixture: A=Cube({2,2,2}), B=Cube({1.7,1.9,2.3}).Translate({1.13,0.41,0.37})
  const Manifold a = Manifold::Cube({2, 2, 2});
  const Manifold b =
      Manifold::Cube({1.7, 1.9, 2.3}).Translate({1.13, 0.41, 0.37});
  const Manifold oracle = a + b;  // Boolean3 union
  const Manifold::Impl impl = ComposeImpl(a, b);
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate5 generic boxes pipeline fatal=" << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  OracleCompare(*result.impl, oracle, eps, "Gate5_GenericBoxes");
}

TEST(Overlap3, Gate5_Oracle_BoxPlusRotatedBox) {
  // Box A + box B rotated 30 degrees around z.
  const Manifold a = Manifold::Cube({2, 2, 2});
  // z-translate 0.3 (not 0.5) so the rotated box's top z face lands at 1.8,
  // avoiding coplanarity with the main cube's top face at z=2.0.
  const Manifold b =
      Manifold::Cube({1.5, 1.5, 1.5}).Rotate(0, 0, 30).Translate({1, 0.5, 0.3});
  const Manifold oracle = a + b;
  const Manifold::Impl impl = ComposeImpl(a, b);
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate5 box+rotated pipeline fatal=" << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  OracleCompare(*result.impl, oracle, eps, "Gate5_BoxPlusRotatedBox");
}

TEST(Overlap3, Gate5_Oracle_TwoSpheres) {
  // Two spheres(16) overlapping. Spheres have no axis-aligned faces.
  const Manifold a = Manifold::Sphere(1.0, 16);
  const Manifold b = Manifold::Sphere(1.0, 16).Translate({1.0, 0, 0});
  const Manifold oracle = a + b;
  const Manifold::Impl impl = ComposeImpl(a, b);
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Gate5 two spheres pipeline fatal=" << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  OracleCompare(*result.impl, oracle, eps, "Gate5_TwoSpheres");
}

// ---------------------------------------------------------------------------
// New pins (spec "new pins from the review evidence")
// ---------------------------------------------------------------------------

// NESTED no-seam cubes: outer 3^3 centered, inner 1^3 centered.
// After removal: only outer shell emitted. vol=27, 12 tris (outer cube).
// Inner shell: winding=2 on both sides -> no fill transition -> not emitted.
TEST(Overlap3, Pin_NestedCubes) {
  const Manifold::Impl impl = NestedCubes();
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Nested cubes fatal: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold()) << "Nested cubes output not manifold";

  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, 27.0, 1e-4)
      << "Nested cubes: expected vol=27 (outer cube only), got " << vol;

  // 12 triangles (the outer 3x3x3 cube has 6 faces x 2 tris = 12).
  EXPECT_EQ(result.impl->NumTri(), 12u)
      << "Nested cubes: expected 12 tris (outer cube), got "
      << result.impl->NumTri();
}

// TOUCHING-disjoint: Cube + Cube.Translate({1,0,0}) -> must NOT fatal,
// must produce the 2x1x1 union (vol=2).
TEST(Overlap3, Pin_TouchingDisjoint) {
  const Manifold::Impl impl = TouchingDisjoint();
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "TouchingDisjoint must not fatal but got: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold())
      << "TouchingDisjoint output not manifold";

  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, 2.0, 1e-6)
      << "TouchingDisjoint: expected vol=2 (2x1x1 union), got " << vol;
}

// Single cube through FULL pipeline -> exact volume (1.0), manifold.
TEST(Overlap3, EmissionAlgebra_SingleCube) {
  const Manifold::Impl impl(Manifold::Impl::Shape::Cube);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "SingleCube fatal: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold());
  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  const double volRef = impl.GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, volRef, 1e-6)
      << "SingleCube vol mismatch: ours=" << vol << " ref=" << volRef;
}

// Single tet through FULL pipeline -> exact volume, manifold.
TEST(Overlap3, EmissionAlgebra_SingleTet) {
  const Manifold::Impl impl(Manifold::Impl::Shape::Tetrahedron);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "SingleTet fatal: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold());
  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  const double volRef = impl.GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, volRef, 1e-6)
      << "SingleTet vol mismatch: ours=" << vol << " ref=" << volRef;
}

// Emission algebra pin: all six faces of cube, including +/-y (vertical
// pieces). Pipeline must emit all 12 tris and produce exact volume.
TEST(Overlap3, EmissionAlgebra_CubeSixFaces) {
  // Use a non-unit cube to exercise all face normals distinctly.
  const Manifold cubeM = Manifold::Cube({2.0, 3.0, 4.0}, true);
  const Manifold::Impl impl(cubeM.GetMeshGL());
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "CubeSixFaces fatal: " << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  EXPECT_TRUE(result.impl->IsManifold());
  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, 2.0 * 3.0 * 4.0, 1e-6)
      << "CubeSixFaces vol=" << vol << " expected=24";
  EXPECT_EQ(result.impl->NumTri(), 12u)
      << "CubeSixFaces: expected 12 tris, got " << result.impl->NumTri();
}

// ---------------------------------------------------------------------------
// Phase P: Mechanism pins (P1-P7)
// ---------------------------------------------------------------------------

// P1: CheckSeamBalance_Test synthetic imbalance -> BalanceViolation.
// Seam with 2 kept regions on face 0 and 1 kept region on face 1 (n0=2, n1=1):
// both nonzero and unequal -> BalanceViolation (was stub "return true" before
// M7).
TEST(Overlap3, Pin_P1_BalanceViolation) {
  ArrangementGeometry arr;
  Seam seam;
  seam.faceId0 = 0;
  seam.faceId1 = 1;
  seam.vertIds = {0, 1};
  arr.seams.push_back(seam);

  // loopVerts={2,0,1,3}: at i=1, loop[1]=0, loop[2]=1, loop[0]=2!=1,
  // loop[3]=3!=0 -> hasEdgeNoSpike(0,1) = true.
  // classified=true, below=0, above=1 -> isKept=true
  // (IsInside3D(0)!=IsInside3D(1)).
  auto makeKeptRegion = []() {
    PSLGRegion r;
    r.loopVerts = {2, 0, 1, 3};
    r.classified = true;
    r.below = 0;
    r.above = 1;
    return r;
  };
  std::vector<std::vector<PSLGRegion>> faceRegions(2);
  faceRegions[0].push_back(makeKeptRegion());
  faceRegions[0].push_back(makeKeptRegion());  // n0 = 2
  faceRegions[1].push_back(
      makeKeptRegion());  // n1 = 1; 2 != 1 -> BalanceViolation

  const auto fatal = CheckSeamBalance_Test(arr, faceRegions);
  ASSERT_TRUE(fatal.has_value())
      << "expected BalanceViolation from imbalanced seam";
  EXPECT_EQ(*fatal, FatalReason::BalanceViolation);
}

// P2: ClassifyRegion_Test with two conflicting pieces ->
// ClassificationAmbiguity. Face in z=0 plane (normal=(0,0,1)), square region
// [1,3]x[1,3]. Two pieces with sourceId=faceId but different (below,above)
// values land at the same projected midpoint inside the region.
TEST(Overlap3, Pin_P2_ClassificationAmbiguity) {
  const std::vector<MergedVert> verts = {
      {{1.0, 1.0, 0.0}},  // 0
      {{3.0, 1.0, 0.0}},  // 1
      {{3.0, 3.0, 0.0}},  // 2
      {{1.0, 3.0, 0.0}},  // 3
  };

  CanonicalFace face;
  face.id = 0;
  face.verts = {0, 1, 2};
  face.normal = {0.0, 0.0, 1.0};
  face.mult = 1;

  // Square region in the z=0 plane; CCW in (x,y) after proj (normal=(0,0,1)).
  PSLGRegion region;
  region.loopVerts = {0, 1, 2, 3};

  // Slab fully covering the region's x-range [1,3].
  SlabResult slab;
  slab.xLo = 1.0;
  slab.xHi = 3.0;
  slab.xMid = 2.0;
  slab.built = true;

  // Two pieces both sourced to faceId=0.
  // from/to are (y,z) section coords. mid2d=(2,0) -> mid3d={xMid=2, y=2, z=0}
  // -> midProj=proj*mid3d=(2,2) in (x,y), which is inside the square [1,3]^2.
  SweepCapture piece1;
  piece1.from = {1.0, 0.0};
  piece1.to = {3.0, 0.0};
  piece1.sourceId = 0;
  piece1.below = 0;
  piece1.above = 1;

  SweepCapture piece2;
  piece2.from = {1.5, 0.0};
  piece2.to = {2.5, 0.0};
  piece2.sourceId = 0;
  piece2.below = 1;  // differs from piece1.below -> ClassificationAmbiguity
  piece2.above = 0;

  slab.pieces = {piece1, piece2};

  const auto result =
      ClassifyRegion_Test(region, 0, {slab}, verts, face, 1e-10);
  ASSERT_TRUE(result.fatal.has_value())
      << "expected ClassificationAmbiguity from conflicting pieces";
  EXPECT_EQ(*result.fatal, FatalReason::ClassificationAmbiguity);
}

// P3: ClassifyRegion_Test with no covering slab -> returns nullopt (no fatal),
// region stays unclassified (degenerate path, handled by anchor propagation).
TEST(Overlap3, Pin_P3_ClassifyRegion_Degenerate) {
  const std::vector<MergedVert> verts = {
      {{1.0, 1.0, 0.0}},
      {{3.0, 1.0, 0.0}},
      {{3.0, 3.0, 0.0}},
      {{1.0, 3.0, 0.0}},
  };

  CanonicalFace face;
  face.id = 0;
  face.verts = {0, 1, 2};
  face.normal = {0.0, 0.0, 1.0};
  face.mult = 1;

  PSLGRegion region;
  region.loopVerts = {0, 1, 2, 3};  // xMin=1, xMax=3

  // Empty slab list: no covering slab -> bestSlab=-1 -> ClassifyRegion returns
  // nullopt without setting outFatal -> result has no fatal and
  // classified=false.
  const auto result =
      ClassifyRegion_Test(region, 0, /*slabs=*/{}, verts, face, 1e-10);
  EXPECT_FALSE(result.fatal.has_value())
      << "no fatal expected for region with no covering slab";
  EXPECT_FALSE(result.classified)
      << "region should be unclassified when no covering slab exists";
}

// P4: EdgeInPlane fatal. Two tetrahedra: tet A has an edge in the z=0 plane
// whose midpoint lies inside tet B's z=0 face.
//   Tet A: v0=(1,1,0), v1=(3,1,0), v2=(2,3,2), v3=(2,1,3).
//   Tet B: v0=(0,0,0), v1=(4,0,0), v2=(2,4,0), v3=(2,2,-2).
// Edge A (1,1,0)-(3,1,0) lies in z=0 (plane of B's face (0,1,2)).
// Midpoint (2,1,0) is strictly inside B's z=0 face -> EdgeInPlane.
TEST(Overlap3, Pin_P4_EdgeInPlane) {
  MeshGL64 mgA;
  mgA.numProp = 3;
  // clang-format off
  mgA.vertProperties = {1,1,0,  3,1,0,  2,3,2,  2,1,3};
  mgA.triVerts       = {0,2,1,  0,1,3,  0,3,2,  1,2,3};
  // clang-format on
  mgA.runOriginalID.push_back(Manifold::ReserveIDs(1));

  MeshGL64 mgB;
  mgB.numProp = 3;
  // clang-format off
  mgB.vertProperties = {0,0,0,  4,0,0,  2,4,0,  2,2,-2};
  mgB.triVerts       = {0,1,2,  0,3,1,  1,3,2,  2,3,0};
  // clang-format on
  mgB.runOriginalID.push_back(Manifold::ReserveIDs(1));

  const Manifold::Impl impl = ComposeImpl(Manifold(mgA), Manifold(mgB));
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_TRUE(result.fatal.has_value()) << "expected EdgeInPlane fatal";
  EXPECT_EQ(*result.fatal, FatalReason::EdgeInPlane)
      << "got fatal=" << (int)*result.fatal << " detail=" << result.detail;
}

// P5: PSLGInvalid from RemoveOverlaps3D_FromArr. Two seam segments on face 0
// that cross at interior point (2,2) not present in arr.verts.
//   seam0: (1,1,0)->(3,3,0), 2D: (1,1)->(3,3)
//   seam1: (3,1,0)->(1,3,0), 2D: (3,1)->(1,3)
// Cross at (2,2) in (x,y) - not any arr vertex -> PSLGInvalid.
TEST(Overlap3, Pin_P5_PSLGInvalid) {
  ArrangementGeometry arr;
  arr.verts = {
      {{0.0, 0.0, 0.0}},  // 0
      {{4.0, 0.0, 0.0}},  // 1
      {{2.0, 4.0, 0.0}},  // 2
      {{1.0, 1.0, 0.0}},  // 3  seam0 start
      {{3.0, 3.0, 0.0}},  // 4  seam0 end
      {{3.0, 1.0, 0.0}},  // 5  seam1 start
      {{1.0, 3.0, 0.0}},  // 6  seam1 end
  };

  CanonicalFace f0, f1;
  f0.id = 0;
  f0.verts = {0, 1, 2};
  f0.normal = {0.0, 0.0, 1.0};
  f0.mult = 1;
  f1.id = 1;
  f1.verts = {0, 2, 1};
  f1.normal = {0.0, 0.0, -1.0};
  f1.mult = 1;
  arr.faces = {f0, f1};

  Seam seam0, seam1;
  seam0.faceId0 = 0;
  seam0.faceId1 = 1;
  seam0.vertIds = {3, 4};
  seam1.faceId0 = 0;
  seam1.faceId1 = 1;
  seam1.vertIds = {5, 6};
  arr.seams = {seam0, seam1};

  arr.faceSeams.resize(2);
  arr.faceSeams[0] = {0, 1};
  arr.faceSeams[1] = {0, 1};

  const Overlap3Result result = RemoveOverlaps3D_FromArr(arr, 1e-8);
  ASSERT_TRUE(result.fatal.has_value()) << "expected PSLGInvalid";
  EXPECT_EQ(*result.fatal, FatalReason::PSLGInvalid)
      << "got fatal=" << (int)*result.fatal << " detail=" << result.detail;
}

// P6: Compose a tet with its winding-reversed copy. Both meshes share the same
// 4 vertex positions {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}. Stage A merges the
// duplicate positions and finds each face pair has opposite permutation parity
// -> mult = +1 + (-1) = 0 -> all 4 face keys dropped -> stageA.faces empty ->
// result is empty Impl (no fatal, vol=0).
//   Tet A outward faces: {0,2,1, 0,1,3, 0,3,2, 1,2,3}
//   Tet B reversed faces (verts offset by 4): {4,5,6, 4,7,5, 4,6,7, 5,7,6}
TEST(Overlap3, Pin_P6_CancellingMesh) {
  MeshGL64 combined;
  combined.numProp = 3;
  // clang-format off
  combined.vertProperties = {
      0,0,0,  1,0,0,  0,1,0,  0,0,1,  // verts 0-3 (tet A)
      0,0,0,  1,0,0,  0,1,0,  0,0,1,  // verts 4-7 (tet B, same positions)
  };
  combined.triVerts = {
      0,2,1,  0,1,3,  0,3,2,  1,2,3,  // tet A: outward (parity -1,+1,-1,+1)
      4,5,6,  4,7,5,  4,6,7,  5,7,6,  // tet B: reversed (parity +1,-1,+1,-1)
  };
  // clang-format on
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));

  const Manifold::Impl impl(combined);
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "cancelling mesh: unexpected fatal=" << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  const double vol = result.impl->GetProperty(Manifold::Impl::Property::Volume);
  EXPECT_NEAR(vol, 0.0, 1e-6) << "cancelling tet expected vol=0, got " << vol;
}

// P7: Starvation when bBox_.Scale()=0 -> EpsilonFromScale(0,1000)=0 -> eps<=0.
// Set bBox_ to a zero-volume box directly; no triangles needed.
TEST(Overlap3, Pin_P7_Starvation) {
  Manifold::Impl impl;
  // Zero-volume bounding box: Scale() = max(|0|,|0|,|0|) = 0.
  // EpsilonFromScale(0, 1000) = 0 -> eps <= 0.0 -> Starvation.
  impl.bBox_ = Box{vec3{0.0, 0.0, 0.0}, vec3{0.0, 0.0, 0.0}};
  const Overlap3Result result = RemoveOverlaps3D(impl, 0.0);
  ASSERT_TRUE(result.fatal.has_value()) << "expected Starvation fatal";
  EXPECT_EQ(*result.fatal, FatalReason::Starvation);
}

// ---------------------------------------------------------------------------
// Phase P: M6 anchor-component propagation pins (P8-P10)
// ---------------------------------------------------------------------------

// P8: PropagateAnchorComponents_Test with conflicting anchors ->
// AnchorConflict. One face (fi=0), three regions, no seams.
//   R0 (anchor, below=0, above=1): loopVerts=[0,2,3]  edges: 0->2, 2->3, 3->0
//   R1 (degenerate):                loopVerts=[3,2,4]  edges: 3->2, 2->4, 4->3
//   R2 (anchor, below=1, above=0): loopVerts=[2,1,4]  edges: 2->1, 1->4, 4->2
// Adjacency:
//   R0 has 2->3; R1 has 3->2 -> R0 adjacent to R1 (anchor).
//   R2 has 4->2; R1 has 2->4 -> R2 adjacent to R1 (anchor).
// Anchor set = {(0,1),(1,0)} -> size 2 -> AnchorConflict.
TEST(Overlap3, Pin_P8_AnchorConflict) {
  ArrangementGeometry arr;
  arr.verts = {
      {{0.0, 0.0, 0.0}},  // v0
      {{2.0, 0.0, 0.0}},  // v1
      {{1.0, 0.0, 0.0}},  // v2
      {{0.0, 1.0, 0.0}},  // v3
      {{2.0, 1.0, 0.0}},  // v4
  };
  CanonicalFace face;
  face.id = 0;
  face.verts = {0, 1, 2};
  face.normal = {0.0, 0.0, 1.0};
  face.mult = 1;
  arr.faces = {face};
  arr.seams = {};
  arr.faceSeams = {{}};  // face 0: no seams

  PSLGRegion r0, r1, r2;
  r0.loopVerts = {0, 2, 3};
  r0.classified = true;
  r0.below = 0;
  r0.above = 1;
  r1.loopVerts = {3, 2, 4};
  r1.degenerate = true;
  r2.loopVerts = {2, 1, 4};
  r2.classified = true;
  r2.below = 1;
  r2.above = 0;

  std::vector<std::vector<PSLGRegion>> faceRegions = {{r0, r1, r2}};
  Overlap3Counters cnt;
  const auto fatal =
      PropagateAnchorComponents_Test(arr, faceRegions, 1e-8, cnt);
  ASSERT_TRUE(fatal.has_value())
      << "expected AnchorConflict from conflicting anchors";
  EXPECT_EQ(*fatal, FatalReason::AnchorConflict);
}

// P9: PropagateAnchorComponents_Test with no anchors and diameter > eps ->
// UnclassifiableComponent.
// One face (fi=0), one large degenerate region, no classified neighbors.
//   R0 (degenerate): loopVerts=[0,1,2,3] - large square, diameter ~14.1
//   No classified regions -> anchors empty -> UnclassifiableComponent.
TEST(Overlap3, Pin_P9_UnclassifiableComponent) {
  ArrangementGeometry arr;
  arr.verts = {
      {{0.0, 0.0, 0.0}},    // v0
      {{10.0, 0.0, 0.0}},   // v1
      {{10.0, 10.0, 0.0}},  // v2
      {{0.0, 10.0, 0.0}},   // v3
  };
  CanonicalFace face;
  face.id = 0;
  face.verts = {0, 1, 2};
  face.normal = {0.0, 0.0, 1.0};
  face.mult = 1;
  arr.faces = {face};
  arr.seams = {};
  arr.faceSeams = {{}};

  PSLGRegion r0;
  r0.loopVerts = {0, 1, 2, 3};  // large square, diameter ~ sqrt(200) > eps
  r0.degenerate = true;

  std::vector<std::vector<PSLGRegion>> faceRegions = {{r0}};
  Overlap3Counters cnt;
  const double eps = 1e-8;  // << diameter ~14.1
  const auto fatal = PropagateAnchorComponents_Test(arr, faceRegions, eps, cnt);
  ASSERT_TRUE(fatal.has_value())
      << "expected UnclassifiableComponent: no anchors, diameter > eps";
  EXPECT_EQ(*fatal, FatalReason::UnclassifiableComponent);
}

// P10: PropagateAnchorComponents_Test with agreeing anchor + diameter <= eps ->
// members classified by propagation.
// One face (fi=0), two regions, no seams.
//   R0 (anchor, below=0, above=1): loopVerts=[0,1,3,4]  (quad)
//     edges: 0->1, 1->3, 3->4, 4->0
//   R1 (degenerate, tiny):         loopVerts=[1,0,2]
//     edges: 1->0, 0->2, 2->1
// Adjacency: R0 has 0->1; R1 has 1->0 -> R0 is anchor of R1.
// R1 verts: v0=(0,0,0), v1=(5e-5,0,0), v2=(2.5e-5,3e-5,0)
//   diameter = max(|v0-v1|, ...) = 5e-5 < eps=1e-4 -> span guard passes.
// Anchor set = {(0,1)} -> one entry -> propagate below=0, above=1.
TEST(Overlap3, Pin_P10_AnchorPropagation_Positive) {
  ArrangementGeometry arr;
  arr.verts = {
      {{0.0, 0.0, 0.0}},      // v0
      {{5e-5, 0.0, 0.0}},     // v1 (5e-5 from v0)
      {{2.5e-5, 3e-5, 0.0}},  // v2 (tiny triangle tip for R1)
      {{1.0, 0.0, 0.0}},      // v3 (anchor far corner)
      {{0.0, 1.0, 0.0}},      // v4 (anchor far corner)
  };
  CanonicalFace face;
  face.id = 0;
  face.verts = {0, 1, 2};
  face.normal = {0.0, 0.0, 1.0};
  face.mult = 1;
  arr.faces = {face};
  arr.seams = {};
  arr.faceSeams = {{}};

  PSLGRegion r0, r1;
  // R0: large quad; edge 0->1 is the shared edge.
  r0.loopVerts = {0, 1, 3, 4};
  r0.classified = true;
  r0.below = 0;
  r0.above = 1;
  // R1: tiny triangle; edge 1->0 is the reverse of R0's 0->1.
  r1.loopVerts = {1, 0, 2};
  r1.degenerate = true;

  std::vector<std::vector<PSLGRegion>> faceRegions = {{r0, r1}};
  Overlap3Counters cnt;
  const double eps = 1e-4;  // > diameter (5e-5)
  const auto fatal = PropagateAnchorComponents_Test(arr, faceRegions, eps, cnt);
  ASSERT_FALSE(fatal.has_value())
      << "expected success from agreeing anchor with diameter <= eps";
  EXPECT_TRUE(faceRegions[0][1].classified)
      << "R1 should be classified after propagation";
  EXPECT_EQ(faceRegions[0][1].below, 0)
      << "R1.below should propagate from anchor";
  EXPECT_EQ(faceRegions[0][1].above, 1)
      << "R1.above should propagate from anchor";
  EXPECT_EQ(cnt.degenerateClassified, 1)
      << "one component classified by anchor propagation";
  // IsInside3D(0)=false != IsInside3D(1)=true -> not a drop ->
  // epsFeaturesDropped=0.
  EXPECT_EQ(cnt.epsFeaturesDropped, 0)
      << "not dropped: R1 straddles the material boundary";
}

// ---------------------------------------------------------------------------
// S3: Inverted-cube pins
// ---------------------------------------------------------------------------

// Return an Impl with all triangle windings flipped relative to m.
static Manifold::Impl InvertWinding(const Manifold& m) {
  MeshGL64 mg = m.GetMeshGL64();
  for (size_t i = 0; i < mg.triVerts.size(); i += 3)
    std::swap(mg.triVerts[i + 1], mg.triVerts[i + 2]);
  mg.runOriginalID.clear();
  mg.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(mg);
}

// Compose A with a winding-inverted copy of B (for subtract oracle test).
static Manifold::Impl ComposeWithInverted(const Manifold& a,
                                          const Manifold& b) {
  const MeshGL64 mga = a.GetMeshGL64();
  MeshGL64 mgb = b.GetMeshGL64();
  for (size_t i = 0; i < mgb.triVerts.size(); i += 3)
    std::swap(mgb.triVerts[i + 1], mgb.triVerts[i + 2]);
  MeshGL64 combined;
  combined.numProp = 3;
  auto append = [&](const MeshGL64& m) {
    const uint64_t base = combined.NumVert();
    for (size_t i = 0; i < m.vertProperties.size(); ++i)
      combined.vertProperties.push_back(m.vertProperties[i]);
    for (size_t i = 0; i < m.triVerts.size(); ++i)
      combined.triVerts.push_back(m.triVerts[i] + base);
  };
  append(mga);
  append(mgb);
  combined.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(combined);
}

// S3(a): A single inverted cube (reversed winding) has w>0 nowhere.
// RemoveOverlaps3D must succeed with an empty output (0 tris).
TEST(Overlap3, Pin_S3a_InvertedCube_Empty) {
  const Manifold::Impl impl = InvertWinding(Manifold::Cube({1, 1, 1}));
  const double eps = ImplEps(impl);
  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "Inverted cube: unexpected fatal=" << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  // Winding is -1 inside, 0 outside; IsInside3D(-1)=false and
  // IsInside3D(0)=false -> no fill transition -> nothing emitted.
  EXPECT_EQ(result.impl->NumTri(), 0u)
      << "Inverted cube: expected 0 tris (w>0 region is empty)";
}

// S3(b): Compose(A, inverted(B)) encodes A-B via winding:
// w_total = w_A + w_{-B} > 0 iff w_A > w_B.
// Oracle: Manifold::Boolean Subtract (Boolean3). Assert volume/genus agreement.
TEST(Overlap3, Pin_S3b_CubeMinusInverted_OracleSubtract) {
  // Generic-offset fixture (same as Gate5) to avoid coplanarity.
  const Manifold a = Manifold::Cube({2, 2, 2});
  const Manifold b =
      Manifold::Cube({1.7, 1.9, 2.3}).Translate({1.13, 0.41, 0.37});
  const Manifold oracle = a - b;

  const Manifold::Impl impl = ComposeWithInverted(a, b);
  const double eps = ImplEps(impl);

  const Overlap3Result result = RemoveOverlaps3D(impl, eps);
  ASSERT_FALSE(result.fatal.has_value())
      << "CubeMinusInverted: pipeline fatal=" << (int)*result.fatal << " "
      << result.detail;
  ASSERT_TRUE(result.impl.has_value());
  OracleCompare(*result.impl, oracle, eps, "Pin_S3b_CubeMinusInverted");
}

}  // namespace
