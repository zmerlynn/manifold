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

// Overlap3 test suite: contract pins for RemoveOverlaps3D, the per-component
// regularization operator. Design: docs/Regularize3D.md.  Pins are authored
// red-first and mutation-verified; skips exist only for missing model files.

#include "../src/overlap3.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <optional>
#include <random>
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
Manifold::Impl ComposeImpl(const Manifold& a, const Manifold& b) {
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
Manifold::Impl ComposeMany(std::initializer_list<Manifold> ms) {
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

double ImplEps(const Manifold::Impl& impl) {
  return EpsilonFromScale(impl.bBox_.Scale(), 1000);
}

}  // namespace

// ===========================================================================
// RETIRED: v3 sweep-mechanism pins (one line each; mechanism deleted).
//
// When RemoveOverlaps3D became the per-component regularization operator
// (docs/Regularize3D.md), the v3 sweep pipeline - Canonicalize/FindSeams/
// BuildSlabs/SlabResolver/EmitCaps/EmitStrips/SweepEmit and its Gate1-4
// section gates plus the ArrangementBudget/SubEpsFeature refusals - was
// deleted.  The pins that asserted that mechanism retire with it; their
// history lives on the v3/v4 branches (explore/sweep-plane-3d-v3, -v4).
// The corpus FIXTURES are re-pinned against the sole-impl contracts in the
// Regularize section below (clean pairs = bitwise pass-through + Decompose
// count; dirty singles = resolve oracle-true / named fail-closed).
//
// Retired pins:
//   Gate1_EventParity_GenericBoxes, Gate1_EventParity_TwoTets,
//   Gate2_SectionValidity_SingleCube, Gate2_SectionValidity_GenericBoxes,
//   Gate2_SectionValidity_TwoTets,
//   Gate2_SectionValidity_ThreeOverlappingBoxes,
//   Gate3_ManifoldOutput_SingleCube, Gate3_ManifoldOutput_SingleTet,
//   Gate3_ThreeOverlappingBoxes, Gate4a_Wedges8_MustResolve,
//   Gate4b_NearParallel_1e6_MustResolve, Gate4c_HullMask_MustResolve,
//   Gate4d_NearParallel_ResolveOrFailClosed,
//   Gate4e_SubResolutionChain_FailClosed, Gate4f_Wedges_TinyOffset_FailClosed,
//   Gate5_Oracle_GenericBoxes, Gate5_Oracle_BoxPlusRotatedBox,
//   Gate5_Oracle_TwoSpheres, Gate5_Oracle_ThreeOverlappingBoxes,
//   Pin_M1_TripleCritical, Pin_NestedCubes, Pin_TouchingDisjoint,
//   Pin_P4_EdgeOnFace_Touching, Pin_P4b_EdgeOnFace_OffMidpoint,
//   Pin_P6_CancellingMesh, Pin_P7_SubEpsInput, Pin_S3a_InvertedCube_Empty,
//   Pin_S3b_CubeMinusInverted_OracleSubtract,
//   Pin_P12_InteriorIsland_StampThroughFace, Pin_PerCriticalCaps,
//   Pin_StripSubdivFromCap, Pin_OneArrangementPerCritical,
//   Pin_InPlaneSkeletonCriticals, EmissionAlgebra_SingleCube,
//   EmissionAlgebra_SingleTet, EmissionAlgebra_CubeSixFaces,
//   Coplanar_StackedPerp_Oracle, Coplanar_StackedParallel_Oracle,
//   Coplanar_SharedWall_Oracle, Coplanar_SameOriented_Oracle,
//   Coplanar_ThreeFaceGroup_Oracle, Coplanar_MixedOrientation_Oracle,
//   Coplanar_EpsChain_Oracle, Coplanar_RazorBand_Recorded,
//   Coplanar_PerpFacesSubEpsApart_Resolves, Coplanar_InvertedStacking,
//   Touch_EdgeEdge_Cubes, Touch_VertexOnly_Cubes, Corpus_Havocglass8_Recorded,
//   Corpus_GenericTwin7863_Recorded, Corpus_GenericTwin7081_Recorded,
//   Corpus_Offset1_Recorded, Corpus_OpenscadNonmanifold_Recorded,
//   Corpus_SelfIntersectA_Recorded, Corpus_SelfIntersectB_Recorded.
// ===========================================================================

// ===========================================================================
// Regularization operator (docs/Regularize3D.md) - GATE + DISPATCH.
// RemoveOverlaps3D: decompose by connectivity -> per-component gate (validity +
// IsSelfIntersecting + within-component coplanar overlap) -> early-exit clean
// -> route dirty to the resolver -> re-gate -> compose back by concatenation.
// These pins are mutation-verified in the lane notebook.
// ===========================================================================

// Bit-pattern mesh identity: vert positions (raw doubles) and triangle vertex
// ids, compared exactly.  Returns a bool; never prints the coordinates (house
// discipline: equality claims by bit-pattern compare, no coordinate dumps).
static bool BitIdenticalMesh(const Manifold::Impl& a, const Manifold::Impl& b) {
  if (a.vertPos_.size() != b.vertPos_.size()) return false;
  if (a.halfedge_.size() != b.halfedge_.size()) return false;
  for (size_t i = 0; i < a.vertPos_.size(); ++i)
    if (std::memcmp(&a.vertPos_[i], &b.vertPos_[i], sizeof(vec3)) != 0)
      return false;
  for (size_t h = 0; h < a.halfedge_.size(); ++h)
    if (a.halfedge_.Start(static_cast<int>(h)) !=
        b.halfedge_.Start(static_cast<int>(h)))
      return false;
  return true;
}

// A single connected VALID 2-manifold that SELF-INTERSECTS: a centered unit
// cube with one corner dragged PAST the opposite corner, so the three tris
// incident to that corner sweep across the whole body and pierce the far
// faces.  The warp leaves topology unchanged (still IsManifold && Is2Manifold);
// the pierce makes IsSelfIntersecting true (verified: a shallow poke exits
// through a side and is NOT flagged, so the target is past the far corner).
// Cheap CI-safe stand-in for the heavy corpus self-intersectors (siA/siB are
// ~17k tris).
static Manifold::Impl PokedCube() {
  const Manifold cube = Manifold::Cube({1, 1, 1}, true);  // [-0.5, 0.5]^3
  const Manifold poked = cube.Warp([](vec3& v) {
    if (v.x > 0 && v.y > 0 && v.z > 0) v = vec3(-1.0, -1.0, -1.0);
  });
  return Manifold::Impl(poked.GetMeshGL64());
}

// Independent winding-number volume oracle (Van Oosterom-Strackee solid
// angle over the raw soup; grid flood of round(GWN)>=1) for the resolver
// acceptance battery (reg3d-wjump).
static double GwnVolumeOracle(const Manifold::Impl& in, vec3 lo, vec3 hi,
                              int G) {
  const int nTri = static_cast<int>(in.NumTri());
  std::vector<std::array<vec3, 3>> T(nTri);
  for (int t = 0; t < nTri; ++t)
    for (int k = 0; k < 3; ++k)
      T[t][k] = in.vertPos_[in.halfedge_.Start(3 * t + k)];
  auto gwn = [&](vec3 p) {
    double s = 0.0;
    for (int t = 0; t < nTri; ++t) {
      const vec3 a = T[t][0] - p, b = T[t][1] - p, c = T[t][2] - p;
      const double la0 = la::length(a), lb = la::length(b), lc = la::length(c);
      const double num = la::dot(a, la::cross(b, c));
      const double den = la0 * lb * lc + la::dot(a, b) * lc +
                         la::dot(b, c) * la0 + la::dot(c, a) * lb;
      s += std::atan2(num, den);
    }
    return s / (2.0 * 3.14159265358979323846);  // winding number
  };
  long inCnt = 0, tot = 0;
  const vec3 d = hi - lo;
  for (int ix = 0; ix < G; ++ix)
    for (int iy = 0; iy < G; ++iy)
      for (int iz = 0; iz < G; ++iz) {
        vec3 p(lo.x + d.x * (ix + 0.5) / G, lo.y + d.y * (iy + 0.5) / G,
               lo.z + d.z * (iz + 0.5) / G);
        if (std::lround(gwn(p)) >= 1) ++inCnt;
        ++tot;
      }
  return d.x * d.y * d.z * (double)inCnt / (double)tot;
}

// Pin 1: clean-input identity.  A clean single connected component passes the
// gate and is copied through BITWISE-unchanged (mesh geometry + topology).
TEST(Overlap3, Regularize_CleanSingleComponent_BitwisePassThrough) {
  const Manifold::Impl in(Manifold::Impl::Shape::Cube);
  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value())
      << "clean cube must not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 1);
  EXPECT_EQ(r.counters.clean, 1);
  EXPECT_EQ(r.counters.dirty, 0);
  EXPECT_EQ(r.counters.failClosed, 0);
  EXPECT_TRUE(BitIdenticalMesh(*r.impl, in))
      << "clean single component must pass through bitwise-unchanged";
}

// Pin 2: multi-component dispatch, all clean.  Two disjoint cubes decompose
// into two components, both early-exit; composed back without fusion, the
// output is still two topological components.
TEST(Overlap3, Regularize_MultiComponent_AllClean_DispatchCounts) {
  const Manifold a = Manifold::Cube({1, 1, 1});
  const Manifold b = Manifold::Cube({1, 1, 1}).Translate({3, 0, 0});
  const Manifold::Impl in = ComposeImpl(a, b);
  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value()) << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 2);
  EXPECT_EQ(r.counters.clean, 2);
  EXPECT_EQ(r.counters.dirty, 0);
  EXPECT_EQ(r.counters.failClosed, 0);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Decompose().size(), 2u)
      << "compose-back must not fuse the two components";
}

// Pin 3: multi-component dispatch, one clean + one dirty.  The clean cube
// early-exits; the poked cube routes to the resolver and now RESOLVES (the
// everted-spike arrangement closes once the shares-vertex genuine-crossing
// recovery records the crossings the broadphase skip dropped, reg3d-wjump). The
// dispatch counts show one regularized, none fail-closed; compose-back keeps
// the two components separate (non-fusion).
TEST(Overlap3, Regularize_MultiComponent_CleanPlusDirty_DispatchCounts) {
  const Manifold clean = Manifold::Cube({1, 1, 1}).Translate({3, 0, 0});
  const Manifold dirtyM(GetMeshGLImpl<double, uint64_t>(PokedCube(), -1));
  const Manifold::Impl in = ComposeImpl(clean, dirtyM);
  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value()) << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 2);
  EXPECT_EQ(r.counters.clean, 1);
  EXPECT_EQ(r.counters.dirty, 1);
  EXPECT_EQ(r.counters.regularized, 1);
  EXPECT_EQ(r.counters.failClosed, 0);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "resolved output must be clean";
}

// Pin 4: a dirty single component (the poked cube) routes to the resolver and
// now RESOLVES to the boundary of {w_S>=1}.  The everted +++->--- corner is a
// near-triple-point where three faces sharing the collapsed spike vertex cross
// transversally OFF the shared vertex; the old shares-vertex broadphase skip
// dropped those crossings, leaving an incomplete arrangement whose emission
// opened (unbalanced fans - the 8-edge-hole family, reg3d-c1a).  The
// reg3d-wjump shares-vertex genuine-crossing recovery records them (SoS-decided
// pierces, seam [V, offVertexP]), completing the arrangement so the per-face
// witness rule classifies the everted (w=-1, dropped) region apart from the
// {w>=1} boundary. Oracle-graded: the emitted volume matches an INDEPENDENT
// winding-number MC oracle, is tol-invariant, non-self-intersecting, and a
// valid manifold.
TEST(Overlap3, Regularize_DirtySingleComponent_RoutesToResolver) {
  const Manifold::Impl dirty = PokedCube();
  ASSERT_TRUE(dirty.IsManifold() && dirty.Is2Manifold())
      << "fixture must be a valid 2-manifold";
  ASSERT_TRUE(dirty.IsSelfIntersecting())
      << "fixture must self-intersect (else it is not a dirty component)";
  const double eps = ImplEps(dirty);
  const RegularizeResult r = RemoveOverlaps3D(dirty, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << "must resolve, not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 1);
  EXPECT_EQ(r.counters.dirty, 1);
  EXPECT_EQ(r.counters.regularized, 1);
  EXPECT_EQ(r.counters.failClosed, 0);

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError) << "output not a manifold";
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output still self-intersects";

  // GEOMETRIC CORRECTNESS: independent GWN volume oracle over the raw soup.
  const double oracle = GwnVolumeOracle(dirty, vec3(-1.05), vec3(0.55), 120);
  const double vol = out.Volume();
  EXPECT_NEAR(vol, oracle, 0.02)
      << "resolved {w>=1} volume off the winding oracle";

  // TOL-INVARIANCE: the topology is decided from input data, not rounded
  // coords.
  const RegularizeResult r2 = RemoveOverlaps3D(dirty, eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * std::abs(vol) + 1e-9)
      << "enclosed volume is not tol-invariant";
}

// White-box classification pin (reg3d-s7b + reg3d-arr): the clean-face
// retention inside the resolver's EmitCleanFaces must be PER-FACE, not
// per-patch.  The everted-corner carrier (PokedCube: the +++ corner collapsed
// onto ---) puts genuine {w_S>=1} boundary clean faces (own +n winding == 0)
// and exterior clean faces (own winding == -1) in ONE clean-clean-connected
// patch.  The old flood classified the whole patch by one representative; here
// the representative is exterior (w=-1) -> it DROPPED the boundary faces,
// leaving the open fan the carrier fails on - AND, in the mirror config where
// the representative is a boundary face, it would EMIT the exterior faces = a
// silent wrong retain (this is the latent-wrongness class, not merely a
// fail-closed).  Per-face is uniformly correct: keep iff the face's OWN winding
// == 0.
//
// This pins the FIX at the classification level even though the whole carrier
// still fails downstream (the unmatched pierce points at the everted spike are
// a separate arrangement-incompleteness wall, reg3d Stage 2/3).  MUTATION:
// revert EmitCleanFaces to the per-patch representative flood -> the boundary
// faces (own winding == 0) are dropped with their patch -> `kept != (ownWinding
// == 0)`
// -> this pin REDs.  Bitwise-safe on the resolving fixtures (siA/siB: every
// face in a patch reaches the same decision, so the emitted set is unchanged).
TEST(Overlap3, Regularize_CleanFacePerFace_RetainsEvertedBoundary) {
  const Manifold::Impl in = PokedCube();
  const CleanFaceProbe p = ClassifyCleanFaces_Probe(in);
  ASSERT_FALSE(p.faceIdx.empty()) << "PokedCube must reach EmitCleanFaces with "
                                     "clean faces (arrangement prefix must not "
                                     "fail before the clean pass)";
  int nBoundary = 0, nExterior = 0, nUncertain = 0;
  for (size_t i = 0; i < p.faceIdx.size(); ++i) {
    // Per-face invariant: a clean face is retained IFF its OWN winding == 0.
    // This is the mutation-sensitive assertion - the flood breaks it on the
    // mixed everted-corner patch.
    const bool ownBoundary = p.ownWinding[i] == 0;
    EXPECT_EQ(static_cast<bool>(p.kept[i]), ownBoundary)
        << "clean face " << p.faceIdx[i] << " ownWinding=" << p.ownWinding[i]
        << " kept=" << static_cast<int>(p.kept[i])
        << " - retention must follow the face's OWN winding, not a patch "
           "representative";
    if (p.ownWinding[i] == kWindingUncertain)
      ++nUncertain;
    else if (ownBoundary)
      ++nBoundary;
    else if (p.ownWinding[i] < 0)
      ++nExterior;
  }
  EXPECT_EQ(nUncertain, 0) << "no clean face's winding may be filter-uncertain "
                              "on this fixture (else the pin is vacuous)";
  // Non-vacuity: the everted corner is a GENUINELY MIXED patch (>=1 boundary
  // AND
  // >=1 exterior clean face), which is exactly the config the flood mislabels.
  EXPECT_GT(nBoundary, 0) << "everted corner must present >=1 genuine boundary "
                             "clean face (own winding == 0)";
  EXPECT_GT(nExterior, 0) << "and >=1 exterior clean face (own winding < 0) - "
                             "the mixed patch the flood dropped wholesale";
}

// Acceptance: the corpus single-
// shell self-intersectors self_intersectA/B - genuine w_S in {0,1,2} dirty
// components with NO negative winding (the clean, safe-by-margin resolver
// target; PokedCube reaches w_S=-1, an openscad-class negative-winding
// KNOWN-OPEN, so it is deliberately NOT the resolve fixture) - must be
// REGULARIZED to the boundary of {w_S >= 1}.  The acceptance battery is
// GEOMETRIC, not just "some clean shape" (an earlier verify lane found the
// prior pin greened against a stub returning an unrelated clean cube): the
// emitted boundary must enclose the {w_S>=1} volume within an INDEPENDENT
// reference band, be one connected solid, and be tol-invariant.
//
// Reference figures (independent MC winding-integration oracle, N=4e5, three
// cocycle-checked seeds; and a grid flood-fill component oracle stable over
// G=40/56/72 - see reg3d-s2 lane notebook):
//   self_intersectA: vol(w>=1) = 0.1438 +- 1e-3, {w>=1} = 1 component
//   self_intersectB: vol(w>=1) = 0.1446 +- 1e-3, {w>=1} = 1 component
// The volume bands below are wide enough for eps-noise + oracle variance yet
// reject any clean-but-WRONG B (a unit cube = 1.0, a bbox cube ~3.2 are both
// decisively outside).  Scalars compared, no coordinate dumps (house
// discipline).
static void ExpectSelfIntersectorRegularizes(const char* name, double volLo,
                                             double volHi) {
  std::filesystem::path file(__FILE__);
  std::ifstream fin((file.parent_path() / "models" / name).string());
  if (!fin.is_open()) GTEST_SKIP() << "model not found";
  const MeshGL64 mesh = ReadOBJ(fin);
  const Manifold::Impl in(mesh);
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold())
      << name << " fixture must be a valid 2-manifold";
  ASSERT_TRUE(in.IsSelfIntersecting())
      << name << " fixture must be a dirty single component";
  const double eps = ImplEps(in);

  const RegularizeResult r = RemoveOverlaps3D(in, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << name << " B must resolve, not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 1);
  EXPECT_EQ(r.counters.dirty, 1);
  EXPECT_EQ(r.counters.regularized, 1);
  EXPECT_EQ(r.counters.failClosed, 0);

  // Output re-gate: a valid closed 2-manifold with the self-overlap removed.
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError)
      << name << " output is not a valid manifold";
  EXPECT_FALSE(r.impl->IsSelfIntersecting())
      << name << " output still self-intersects";

  // GEOMETRIC CORRECTNESS: the emitted boundary encloses the {w_S>=1} volume
  // (independent MC oracle band) - a clean-but-WRONG B lands outside and reds.
  const double vol = out.Volume();
  EXPECT_GT(vol, volLo) << name << " {w>=1} volume below the reference band";
  EXPECT_LT(vol, volHi) << name << " {w>=1} volume above the reference band";

  // TOPOLOGY SANITY: {w>=1} is one connected solid (grid oracle); genus finite.
  EXPECT_EQ(out.Decompose().size(), 1u)
      << name << " {w>=1} must be exactly one solid";
  EXPECT_GE(out.Genus(), 0) << name << " degenerate genus";

  // TOL-INVARIANCE: the output TOPOLOGY is decided from input data, not from
  // rounded positions, so the enclosed volume is invariant to the weld radius.
  // A second run at a tightened eps must reproduce it (a tol-dependent resolve
  // would drift).
  const RegularizeResult r2 = RemoveOverlaps3D(in, eps * 0.5);
  ASSERT_FALSE(r2.fatal.has_value())
      << name << " tol-variant run fatal: " << r2.detail;
  ASSERT_TRUE(r2.impl.has_value());
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-3 * vol)
      << name << " enclosed volume is not tol-invariant";
}

TEST(Overlap3, Regularize_SelfIntersectA_Regularized) {
  ExpectSelfIntersectorRegularizes("self_intersectA.obj", 0.130, 0.158);
}

TEST(Overlap3, Regularize_SelfIntersectB_Regularized) {
  ExpectSelfIntersectorRegularizes("self_intersectB.obj", 0.130, 0.160);
}

// ---------------------------------------------------------------------------
// The resolver mechanism (docs/Regularize3D.md "the resolver's mechanism") -
// white-box port verification against the FRAGMENT-VALIDATED numbers (v5b-r3/r4
// notebooks). These pins are NOT disabled: they prove the ported ENUMERATION
// and coupled WINDING (the substrate B runs today, on top of which THE BUILD is
// unbuilt) are correct, independent of the boundary-emission wall.  Anchor
// points carry their winding from the independent MC oracle (reg3d-s2
// notebook), cocycle-stable across two unrelated seeds.
// ---------------------------------------------------------------------------
static void ExpectBMechanism(const char* name, int expectSeams, vec3 w1,
                             vec3 w2, vec3 far, vec3 seed, vec3 seed2) {
  std::filesystem::path file(__FILE__);
  std::ifstream fin((file.parent_path() / "models" / name).string());
  if (!fin.is_open()) GTEST_SKIP() << "model not found";
  const MeshGL64 mesh = ReadOBJ(fin);
  const Manifold::Impl in(mesh);

  // ENUMERATION: the ported level-0 pierce enumeration reproduces the
  // fragment's genuine self-crossing count exactly, with no exact-zero ties
  // (safe-by-margin, so the single-global-SoS axis is not exercised on the
  // corpus).
  const std::vector<vec3> probes = {far, w1, w2};
  const ComponentEnumProbe p = EnumerateComponent_Probe(in, probes, seed);
  EXPECT_EQ(p.seamCount, expectSeams) << name << " enumeration seam count";
  EXPECT_EQ(p.boundaryTouchPairs, 0)
      << name << " must be safe-by-margin (no exact-zero tie)";

  // WINDING: exterior (0), single-cover (1), and the w=2 self-overlap stratum.
  ASSERT_EQ(p.probeWinding.size(), 3u);
  EXPECT_EQ(p.probeWinding[0], 0) << name << " exterior winding must be 0";
  EXPECT_EQ(p.probeWinding[1], 1) << name << " single-cover winding must be 1";
  EXPECT_EQ(p.probeWinding[2], 2) << name << " double-cover (w=2) winding";

  // PATH-INDEPENDENCE (cocycle): an unrelated second seed reproduces the w=2
  // classification (the coupled winding is single-valued off-surface).
  const ComponentEnumProbe p2 = EnumerateComponent_Probe(in, {w2}, seed2);
  ASSERT_EQ(p2.probeWinding.size(), 1u);
  EXPECT_EQ(p2.probeWinding[0], 2) << name << " winding not path-independent";
}

TEST(Overlap3, Regularize_BMechanism_SelfIntersectA) {
  ExpectBMechanism(
      "self_intersectA.obj", 338,
      {-0.84405223113934791, 0.50434676505537546, 1.3952491390912483},
      {0.088618132029630078, 0.23160096832022406, 0.59518845825120714},
      {8.8552219880000003, 5.7315517699999994, 11.819220435},
      {226.70071018299998, -89.319086551699996, 109.69480060789999},
      {-347.737426747, 170.0367771507, -195.02711216070003});
}

TEST(Overlap3, Regularize_BMechanism_SelfIntersectB) {
  ExpectBMechanism(
      "self_intersectB.obj", 338,
      {-1.5288847346363044, 0.55692625122266082, 1.4831033934933748},
      {-0.60008906037949228, 0.26029227470132027, 0.74852618084744027},
      {8.1568456297999994, 7.8764903549999996, 10.74691606},
      {225.69857201970999, -125.15321454389999, 97.788207489599998},
      {-347.93884238262996, 237.8741813569, -173.12363447680002});
}

// ===========================================================================
// Regularization axis: NEGATIVE WINDING (openscad / subtraction class).
// docs/Regularize3D.md open item "negative winding / subtraction, untested".
// Witness theorem (general form): for an oriented mult-1 face w_below =
// w_above + 1 UNIVERSALLY, so a cell is on d{w_S>=1} iff EXACTLY ONE side has
// w>=1, which reduces to w_above==0; a negative w_above means BOTH sides are
// exterior (w_below = w_above+1 <= 0), so the cell is DROPPED, never
// fail-closed.  B previously fail-closed on w_above<0; this axis removes that.
// ===========================================================================

// Independent generalized winding number (Van Oosterom-Strackee solid-angle
// sum) over an oriented triangle soup - a genuinely different algorithm from
// B's ray-crossing WindingAt, so it is a real cross-check, not a tautology.
static double GWN(const std::vector<std::array<vec3, 3>>& tris, const vec3& p) {
  double sum = 0.0;
  for (const auto& t : tris) {
    const vec3 a = t[0] - p, b = t[1] - p, c = t[2] - p;
    const double la_ = la::length(a), lb = la::length(b), lc = la::length(c);
    const double det = la::dot(a, la::cross(b, c));
    const double den = la_ * lb * lc + la::dot(a, b) * lc +
                       la::dot(b, c) * la_ + la::dot(c, a) * lb;
    sum += 2.0 * std::atan2(det, den);
  }
  return sum / (4.0 * 3.14159265358979323846);
}

static std::vector<std::array<vec3, 3>> SoupTris(const Manifold::Impl& in) {
  std::vector<std::array<vec3, 3>> out(in.NumTri());
  for (int t = 0; t < static_cast<int>(in.NumTri()); ++t)
    for (int k = 0; k < 3; ++k)
      out[t][k] = in.vertPos_[in.halfedge_.Start(3 * t + k)];
  return out;
}

// A single connected VALID 2-manifold that self-intersects AND reaches negative
// soup winding: a sphere with its top cap pushed straight DOWN through the body
// and out the bottom (an inverted cap = a w_S<0 region), then a small (2%)
// smooth generic warp to put it in GENERAL POSITION.  The warp is load-bearing:
// without it the UV-sphere's shared meridian/parallel planes trip the
// exact-zero (SoS) pierce-tie gate BEFORE the winding classifier - negative
// winding and cross-operand exact-zero ties are ENTANGLED in structured
// fixtures, so a generic-position carrier is needed to exercise the
// negative-winding axis alone (reg3d-s3 lane; the axis-aligned openscad carrier
// stays fail-closed at SoS).
static Manifold::Impl PushedCapSphere() {
  const Manifold s = Manifold::Sphere(1.0, 16);
  Manifold p = s.Warp([](vec3& v) {
    if (v.z > 0.5) v.z -= 3.0;
  });
  p = p.Warp([](vec3& v) {
    v += 0.02 * vec3(std::sin(3.1 * v.y + 1.2 * v.z),
                     std::sin(2.7 * v.z + 0.9 * v.x),
                     std::sin(3.3 * v.x + 1.7 * v.y));
  });
  return Manifold::Impl(p.GetMeshGL64());
}

TEST(Overlap3, Regularize_NegativeWinding_PushedCapSphere) {
  const Manifold::Impl in = PushedCapSphere();
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold())
      << "fixture must be a valid 2-manifold";
  ASSERT_TRUE(in.IsSelfIntersecting())
      << "fixture must self-intersect (a dirty component)";
  const auto inTris = SoupTris(in);

  // Sample the bbox once; store input GWN (double + rounded).  Require the
  // fixture to REACH negative soup winding - else it would not exercise the
  // axis (a mutation that removes the inverted cap would red here).
  const Box bb = in.bBox_;
  const vec3 lo = bb.min, hi = bb.max;
  std::mt19937 rng(20260713u);
  std::uniform_real_distribution<double> ux(lo.x, hi.x), uy(lo.y, hi.y),
      uz(lo.z, hi.z);
  const int N = 40000;
  std::vector<vec3> qs;
  std::vector<double> giD;
  qs.reserve(N);
  giD.reserve(N);
  int neg = 0, nIn = 0;
  for (int i = 0; i < N; ++i) {
    const vec3 q(ux(rng), uy(rng), uz(rng));
    const double g = GWN(inTris, q);
    qs.push_back(q);
    giD.push_back(g);
    const long gr = std::lround(g);
    if (gr < 0) ++neg;
    if (gr >= 1) ++nIn;
  }
  ASSERT_GT(neg, N / 100)
      << "fixture must reach negative winding (the axis under test)";

  const double eps = ImplEps(in);
  const RegularizeResult r = RemoveOverlaps3D(in, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << "negative winding must be absorbed, not fail-closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.dirty, 1);
  EXPECT_EQ(r.counters.regularized, 1);
  EXPECT_EQ(r.counters.failClosed, 0);

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError)
      << "output is not a valid manifold";
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output still self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << "{w>=1} must be one solid";

  // INDEPENDENT ORACLE: emitted volume matches vol{w_S>=1}, and MEMBERSHIP
  // (round(GWN_input)>=1) == (round(GWN_output)>=1) at every unambiguous point
  // (a subtly-wrong resolve that passes the wide volume band still disagrees).
  const auto outTris = SoupTris(*r.impl);
  const double bbVol = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
  const double gwnVol = bbVol * nIn / N;
  EXPECT_NEAR(out.Volume(), gwnVol, 0.04 * gwnVol)
      << "emitted volume must match the independent {w>=1} GWN volume";
  int checked = 0, disagree = 0;
  for (int i = 0; i < N; ++i) {
    if (std::abs(giD[i] - std::round(giD[i])) > 0.05) continue;  // near in-surf
    const double go = GWN(outTris, qs[i]);
    if (std::abs(go - std::round(go)) > 0.05) continue;  // near out-surf
    ++checked;
    if ((std::lround(giD[i]) >= 1) != (std::lround(go) >= 1)) ++disagree;
  }
  EXPECT_EQ(disagree, 0) << "input {w>=1} vs output solid disagree at "
                         << disagree << " of " << checked
                         << " unambiguous points";

  // TOL-INVARIANCE: topology decided from input data, so volume is weld-radius
  // invariant.
  const RegularizeResult r2 = RemoveOverlaps3D(in, eps * 0.5);
  ASSERT_FALSE(r2.fatal.has_value()) << "tol-variant run fatal: " << r2.detail;
  ASSERT_TRUE(r2.impl.has_value());
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(out.Volume(), out2.Volume(), 1e-3 * out.Volume())
      << "emitted volume is not tol-invariant";
}

// ===========================================================================
// Regularization axis: EXACT-ZERO TIES (single-global SoS),
// docs/Regularize3D.md stage 6, LANDED.  A level-0 pierce predicate that is an
// EXACT ZERO (a vertex exactly on a face, or an edge grazing a triangle edge -
// a NON-coplanar transversal tie) is now DECIDED by the single-global
// symbolic-perturbation convention (Orient3DSoS: filter fast-path, else the
// exact e^0 sign, else the Edelsbrunner-Mucke cascade) instead of failing
// closed.  The coplanar family stays the FOLD's (coplanar pairs are never
// SoS-perturbed - the s3/s4adj sliver rail).  BridgedCaps (above) is the
// oracle-true RESOLVE of this tie family through the genus-handle junction. The
// carriers below are the RESIDUE: they PASS the SoS gate (no longer
// boundaryTouch) and now fail NARROWER, at emission
// - the coplanar-DOMINATED soups (GT7863) and the collapsed-vertex spike
// (PokedCube) reduce to sub-eps / touching-sheet slivers with no representable
// double-precision manifold boundary (NonManifoldEmission, a hard fail-closed =
// no output; the doc's stage-7 thin-cell territory, honestly named).  These
// pins assert the residue is precise and load-bearing: B never silently
// resolves an exact-tie carrier to (possibly wrong) geometry - the SoS gate
// NARROWED, the thin-cell wall remains.
// ===========================================================================

// PROPERTY PIN for the micro exact tie-test (Orient3DExactSignProbe, the
// owner-contract exact orient3d sign - ONE integer path, adaptive-width, TOTAL:
// no window-fail refusal, no expansion fallback).  Graded properties, each on
// randomized configs:
//  P1 FILTER AGREEMENT: whenever the replicated static Shewchuk filter
//     certifies a sign, the exact sign matches it (the filter is sound, so a
//     disagreement would indict the exact path).
//  P2 ANTISYMMETRY: swapping any two of the four points flips the sign
//     (including on exactly-degenerate configs, where 0 stays 0).
//  P3 CONSTRUCTED EXACT ZEROS: four points on an exact plane -> sign 0.
//  P4 SCALING INVARIANCE: scaling all coordinates by a power of two multiplies
//     the determinant by a POSITIVE factor, so the sign is invariant.  NOTE
//     (reg3d-s6r): uniform pow2 scaling shifts every term's exponent by the
//     SAME amount, so it PRESERVES the term-exponent spread - it does NOT walk
//     any accumulator-width boundary (that is P5's job, below).
//  P5 MIXED-MAGNITUDE / ZERO-STRADDLING: coordinates spanning a huge exponent
//     range in one predicate (fine scale ~2^-996 next to coarse ~2^-10) - a
//     term-exponent spread far past the retired fixed 90-bit window.  This is
//     the wide-spread path that SURVIVES the rework (the adaptive-width
//     accumulator).  A specific adversary is pinned to its arbitrary-precision
//     oracle sign, exact coplanarity at wide spread must still decide 0, and
//     antisymmetry must hold - the path is total, never refusing.
TEST(Overlap3, Orient3DExactSign_PropertyPin) {
  // Replicated static filter (Orient3DFilterSign's math, kept in lockstep).
  auto filter = [](const vec3& a, const vec3& b, const vec3& c, const vec3& d) {
    const vec3 ad = a - d, bd = b - d, cd = c - d;
    const double bc = bd.y * cd.z, cb = cd.y * bd.z;
    const double ca = cd.y * ad.z, ac = ad.y * cd.z;
    const double ab = ad.y * bd.z, ba = bd.y * ad.z;
    const double det = ad.x * (bc - cb) + bd.x * (ca - ac) + cd.x * (ab - ba);
    const double perm = (std::abs(bc) + std::abs(cb)) * std::abs(ad.x) +
                        (std::abs(ca) + std::abs(ac)) * std::abs(bd.x) +
                        (std::abs(ab) + std::abs(ba)) * std::abs(cd.x);
    constexpr double u = 0x1p-53;
    const double errb = (7.0 + 56.0 * u) * u * perm;
    if (errb > 0.0 && std::abs(det) > errb) return det > 0.0 ? 1 : -1;
    return 0;
  };
  std::mt19937_64 rng(20260714u);
  std::uniform_real_distribution<double> U(-3, 3);
  std::uniform_int_distribution<int> Sc(-120, 120);
  int p1 = 0, p2 = 0, p3 = 0, p4 = 0;
  for (int t = 0; t < 20000; ++t) {
    vec3 q[4];
    const bool coplanar = (t % 4 == 0);
    for (auto& v : q) {
      v.x = std::round(U(rng) * 4);
      v.y = std::round(U(rng) * 4);
      v.z = coplanar ? 1.0 : std::round(U(rng) * 4);
    }
    const int s = Orient3DExactSignProbe(q[0], q[1], q[2], q[3]);
    // P1: filter agreement.
    const int f = filter(q[0], q[1], q[2], q[3]);
    if (f != 0 && s != f) ++p1;
    // P2: antisymmetry across all 6 transpositions.
    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j) {
        vec3 w[4] = {q[0], q[1], q[2], q[3]};
        std::swap(w[i], w[j]);
        if (Orient3DExactSignProbe(w[0], w[1], w[2], w[3]) != -s) ++p2;
      }
    // P3: constructed exact zero.
    if (coplanar && s != 0) ++p3;
    // P4: power-of-two scaling invariance (uniform scaling preserves spread).
    const double k = std::ldexp(1.0, Sc(rng));
    if (Orient3DExactSignProbe(q[0] * k, q[1] * k, q[2] * k, q[3] * k) != s)
      ++p4;
  }
  EXPECT_EQ(p1, 0) << "exact sign must agree with every certified filter sign";
  EXPECT_EQ(p2, 0) << "exact sign must be antisymmetric";
  EXPECT_EQ(p3, 0) << "exactly-coplanar points must give sign 0";
  EXPECT_EQ(p4, 0) << "sign must be invariant to power-of-two scaling";

  // P5 MIXED-MAGNITUDE / ZERO-STRADDLING: the wide-spread path the rework
  // keeps. (a) A predicate straddling coordinate zero at fine scale: coords
  // mixing 2^-996 (~1e-300) with 2^-10 (~1e-3) give a term-exponent spread of
  // ~2963 bits, far past the retired fixed 90-bit window.  The exact sign is
  // pinned to its independent arbitrary-precision (Python Fraction) oracle: -1
  // (reg3d-s6r).  The old int256 fixed window REFUSED this config; the adaptive
  // path decides it.
  const vec3 adv[4] = {
      {std::ldexp(1.0, -996), std::ldexp(3.0, -10), 0.0},
      {std::ldexp(5.0, -10), std::ldexp(1.0, -996), std::ldexp(2.0, -10)},
      {std::ldexp(7.0, -10), std::ldexp(2.0, -10), std::ldexp(1.0, -996)},
      {std::ldexp(9.0, -10), std::ldexp(4.0, -10), std::ldexp(3.0, -10)}};
  EXPECT_EQ(Orient3DExactSignProbe(adv[0], adv[1], adv[2], adv[3]), -1)
      << "wide-spread adversary must match its arbitrary-precision oracle";
  EXPECT_EQ(Orient3DExactSignProbe(adv[1], adv[0], adv[2], adv[3]), 1)
      << "and be antisymmetric at wide spread";
  // (b) Exact coplanarity at wide spread must still decide 0 (never a spurious
  // nonzero, never a refusal).  Four points share an exact z, with x,y spanning
  // ~2^-500..2^40 - a spread of ~1000 bits - so the determinant is exactly 0.
  int p5 = 0;
  std::uniform_int_distribution<int> Wide(-500, 40);
  for (int t = 0; t < 4000; ++t) {
    const double z = std::ldexp(1.0, Sc(rng));  // common z -> exactly coplanar
    vec3 w[4];
    for (auto& v : w) {
      v.x = std::ldexp(std::round(U(rng) * 8), Wide(rng));
      v.y = std::ldexp(std::round(U(rng) * 8), Wide(rng));
      v.z = z;
    }
    if (Orient3DExactSignProbe(w[0], w[1], w[2], w[3]) != 0) ++p5;
    // antisymmetry at mixed magnitude on a generic (non-coplanar) config
    vec3 g[4];
    for (auto& v : g) {
      v.x = std::ldexp(std::round(U(rng) * 8), Wide(rng));
      v.y = std::ldexp(std::round(U(rng) * 8), Wide(rng));
      v.z = std::ldexp(std::round(U(rng) * 8), Wide(rng));
    }
    const int sg = Orient3DExactSignProbe(g[0], g[1], g[2], g[3]);
    if (Orient3DExactSignProbe(g[1], g[0], g[2], g[3]) != -sg) ++p5;
  }
  EXPECT_EQ(p5, 0)
      << "wide-spread exact coplanarity must give 0 and stay antisymmetric";
}

// Constructed carrier: a SINGLE self-intersecting component with ISOLATED
// exact-zero ties (vertex-on-plane + pierce-on-edge, NO coplanar overlap - the
// PokedCube's axis-aligned spike pierces the far faces at exact configs).  The
// SoS DECIDES those ties (the carrier passes the exact-zero tie gate) AND the
// shares-vertex genuine-crossing recovery (reg3d-wjump) completes the
// everted-spike arrangement, so the carrier now RESOLVES oracle-true instead of
// opening at emission.  Both axes are exercised: the exact-zero tie is reached
// (boundaryTouchPairs>0) and the resolve is winding-oracle-correct.
TEST(Overlap3, Regularize_ExactZeroTie_Constructed_Resolves) {
  const Manifold::Impl in = PokedCube();
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  ASSERT_TRUE(in.IsSelfIntersecting()) << "carrier must be a dirty component";
  // The carrier genuinely REACHES an exact-zero pierce tie (else it would not
  // exercise the axis): B's enumeration reports filter-tie pairs.
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.boundaryTouchPairs, 0)
      << "carrier must reach an exact-zero tie (the axis under test)";

  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value())
      << "must resolve, not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.regularized, 1);
  EXPECT_EQ(r.counters.failClosed, 0);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting());
  const double oracle = GwnVolumeOracle(in, vec3(-1.05), vec3(0.55), 120);
  EXPECT_NEAR(out.Volume(), oracle, 0.02)
      << "resolved {w>=1} volume off the winding oracle";
}

// Real carrier: the GT7863 twin pair composed as ONE soup.  Connectivity splits
// it into 4 pieces (non-fusion contract: they stay separate).  TWO route dirty:
// one is SELF-INTERSECTING (the 8-edge-hole emission carrier), and one carries
// a WITHIN-component COPLANAR OVERLAP (a double sheet).  NARROWED by
// reg3d-wjump: the shares-vertex genuine-crossing recovery CLOSES the
// self-intersecting component's 8-edge hole - measured white-box below:
// decomposed, that component RESOLVES oracle-true (its {w>=1} volume is
// preserved to the input's signed volume, tol-invariant,
// non-self-intersecting).  The compose still fails closed on the OTHER dirty
// component: the coplanar double sheet fails the fold self-check (in-plane
// cover m=2, but the 3D winding jumps 1 - GWN-verified a single-sheet {w>=1}
// exists, but the doubled connectivity is the Cluster-1a emission
// REPRESENTABILITY wall, reg3d-wjump honest wall).  So the terminal fatal moves
// from NonManifoldEmission to DirtyComponentUnresolved (narrower;
// mutation-verified: reverting the recovery reopens the 8-edge hole and the
// fatal returns to NonManifoldEmission).  Any component fail-closed suppresses
// output, so the whole compose fails closed - never a silent wrong resolve.
TEST(Overlap3, Regularize_ExactZeroTie_GT7863_FailClosed) {
  std::filesystem::path file(__FILE__);
  auto load = [&](const char* n) -> std::optional<MeshGL64> {
    std::ifstream fin((file.parent_path() / "models" / n).string());
    if (!fin.is_open()) return std::nullopt;
    return ReadOBJ(fin);
  };
  const auto mL = load("Generic_Twin_7863.1.t0_left.obj");
  const auto mR = load("Generic_Twin_7863.1.t0_right.obj");
  if (!mL || !mR) GTEST_SKIP() << "model not found";
  MeshGL64 comb;
  comb.numProp = 3;
  auto app = [&](const MeshGL64& m) {
    const uint64_t base = comb.NumVert();
    for (double d : m.vertProperties) comb.vertProperties.push_back(d);
    for (uint64_t t : m.triVerts) comb.triVerts.push_back(t + base);
  };
  app(*mL);
  app(*mR);
  comb.runOriginalID.push_back(Manifold::ReserveIDs(1));
  const Manifold::Impl in(comb);

  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  EXPECT_EQ(r.counters.components, 4) << "four components stay separate";
  EXPECT_EQ(r.counters.clean, 2) << "2 components early-exit clean";
  EXPECT_EQ(r.counters.dirty, 2)
      << "1 self-intersecting + 1 within-component coplanar overlap";
  ASSERT_TRUE(r.fatal.has_value())
      << "the coplanar double-sheet residue must fail closed";
  // NARROWED (reg3d-wjump): the 8-edge-hole emission carrier now resolves; the
  // terminal wall is the coplanar double-sheet fold self-check decline.
  EXPECT_EQ(*r.fatal, FatalReason::DirtyComponentUnresolved) << r.detail;
  EXPECT_FALSE(r.impl.has_value())
      << "any component fail-closed suppresses the whole compose (no partial)";

  // WHITE-BOX: the SELF-INTERSECTING component (the 8-edge hole) now RESOLVES
  // oracle-true - the reg3d-wjump recovery closes the hole.  Decompose, find
  // it, resolve it in isolation, and grade its volume against the input soup's
  // signed volume (the near-tangent self-overlap is measure-~0 so {w>=1} volume
  // is preserved) + tol-invariance + non-self-intersection.
  const Manifold M(GetMeshGLImpl<double, uint64_t>(in, -1));
  int resolvedSI = 0;
  for (const Manifold& c : M.Decompose()) {
    Manifold::Impl ci(c.GetMeshGL64());
    if (!ci.IsSelfIntersecting()) continue;
    const double ceps = ImplEps(ci);
    const RegularizeResult rc = ResolveComponentDirect(ci, ceps);
    ASSERT_FALSE(rc.fatal.has_value())
        << "the 8-edge-hole component must resolve: " << rc.detail;
    ASSERT_TRUE(rc.impl.has_value());
    EXPECT_FALSE(rc.impl->IsSelfIntersecting());
    const double inVol =
        Manifold(GetMeshGLImpl<double, uint64_t>(ci, -1)).Volume();
    const double outVol =
        Manifold(GetMeshGLImpl<double, uint64_t>(*rc.impl, -1)).Volume();
    EXPECT_NEAR(outVol, inVol, 1e-3 * std::abs(inVol))
        << "resolved {w>=1} volume must preserve the near-tangent soup volume";
    const RegularizeResult rc2 = ResolveComponentDirect(ci, ceps * 0.5);
    ASSERT_TRUE(rc2.impl.has_value());
    EXPECT_NEAR(
        outVol,
        Manifold(GetMeshGLImpl<double, uint64_t>(*rc2.impl, -1)).Volume(),
        1e-6 * std::abs(outVol) + 1e-9)
        << "resolved volume not tol-invariant";
    ++resolvedSI;
  }
  EXPECT_EQ(resolvedSI, 1)
      << "exactly one self-intersecting component resolves";
}

// ===========================================================================
// Regularization axis: EXACT-COPLANAR IN-PLANE FOLD (docs/Regularize3D.md
// coplanar axis).  Exactly-coplanar overlapping faces are folded in their
// shared plane with a per-cell signed multiplicity; the generalized retention
// w_below = w_above + m emits the {w_S>=1} boundary (mult-1 stays the existing
// w_above==0 rule).  The self-intersection gate does NOT flag a pure coplanar
// overlap (the doc's R2(i) blind spot), so these carriers are exercised through
// the ResolveComponentDirect hook (the resolver on a soup treated as one dirty
// component).  Each resolve is checked against an INDEPENDENT generalized
// winding-number oracle (Van Oosterom-Strackee over the input soup) - a wrong
// fold lands outside the band or disagrees on membership - plus tol-invariance.
// ===========================================================================

namespace coplanarfold {
struct MB {
  std::vector<double> vp;
  std::vector<uint64_t> tv;
  std::map<std::tuple<double, double, double>, uint64_t> idx;
  uint64_t V(double x, double y, double z) {
    auto key = std::make_tuple(x, y, z);
    auto it = idx.find(key);
    if (it != idx.end()) return it->second;
    uint64_t id = vp.size() / 3;
    vp.push_back(x);
    vp.push_back(y);
    vp.push_back(z);
    idx.emplace(key, id);
    return id;
  }
  void Quad(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    tv.push_back(a);
    tv.push_back(b);
    tv.push_back(c);
    tv.push_back(a);
    tv.push_back(c);
    tv.push_back(d);
  }
  MeshGL64 Mesh() {
    MeshGL64 m;
    m.numProp = 3;
    m.vertProperties = vp;
    m.triVerts = tv;
    m.runOriginalID.push_back(Manifold::ReserveIDs(1));
    return m;
  }
};

// Generalized winding number of an oriented triangle SOUP (Van Oosterom-
// Strackee) - independent of the resolver's ray-crossing winding.
double GWNsoup(const Manifold::Impl& in, const vec3& p) {
  double sum = 0.0;
  const int nTri = static_cast<int>(in.NumTri());
  for (int t = 0; t < nTri; ++t) {
    const vec3 A = in.vertPos_[in.halfedge_.Start(3 * t)] - p;
    const vec3 B = in.vertPos_[in.halfedge_.Start(3 * t + 1)] - p;
    const vec3 C = in.vertPos_[in.halfedge_.Start(3 * t + 2)] - p;
    const double la_ = la::length(A), lb = la::length(B), lc = la::length(C);
    const double det = la::dot(A, la::cross(B, C));
    const double den = la_ * lb * lc + la::dot(A, B) * lc +
                       la::dot(B, C) * la_ + la::dot(C, A) * lb;
    sum += 2.0 * std::atan2(det, den);
  }
  return sum / (4.0 * 3.14159265358979323846);
}

// SLANT-PLUG carrier: box B's footprint STRICTLY inside box A's, both sharing
// A's flat bottom (z=0) AND slanted top plane zt(x,y) EXACTLY (dyadic 0.5/0.25
// coefficients over dyadic coords -> bit-exact coplanarity).  Both caps
// coincide as mult-2 overlaps; B's walls are interior risers (no wall crosses
// another face), so it is a PURE coplanar overlap - no transversal
// entanglement. flipB inverts B (anti-oriented, mult-0 cancellation).
MeshGL64 SlantPlug(bool flipB) {
  MB b;
  auto zt = [](double x, double y) { return 1.0 + 0.5 * x + 0.25 * y; };
  auto sbox = [&](double x0, double x1, double y0, double y1, bool flip) {
    auto T = [&](double x, double y) { return b.V(x, y, zt(x, y)); };
    auto B = [&](double x, double y) { return b.V(x, y, 0); };
    auto Q = [&](uint64_t p, uint64_t q, uint64_t r, uint64_t s) {
      if (!flip)
        b.Quad(p, q, r, s);
      else
        b.Quad(p, s, r, q);
    };
    Q(T(x0, y0), T(x1, y0), T(x1, y1), T(x0, y1));  // slanted top +z
    Q(B(x0, y0), B(x0, y1), B(x1, y1), B(x1, y0));  // flat bottom -z
    Q(B(x0, y0), B(x1, y0), T(x1, y0), T(x0, y0));  // -y wall
    Q(B(x1, y0), B(x1, y1), T(x1, y1), T(x1, y0));  // +x wall
    Q(B(x1, y1), B(x0, y1), T(x0, y1), T(x1, y1));  // +y wall
    Q(B(x0, y1), B(x0, y0), T(x0, y0), T(x0, y1));  // -x wall
  };
  sbox(0, 3, 0, 2, false);          // A (outer)
  sbox(0.5, 1.5, 0.5, 1.5, flipB);  // B (inner plug)
  return b.Mesh();
}

// N boxes sharing A's slanted top + flat bottom.  boxes[0] is the outer box A;
// the rest are strictly-inside plugs.  All the caps coincide in the slant plane
// (coplanar clusters).  Lets a plug's cap edge land on another plug's cap edge
// interior to force a T-JUNCTION in the cap-cluster arrangement.
MeshGL64 MultiPlug(const std::vector<std::array<double, 4>>& boxes) {
  MB b;
  auto zt = [](double x, double y) { return 1.0 + 0.5 * x + 0.25 * y; };
  auto sbox = [&](double x0, double x1, double y0, double y1) {
    auto T = [&](double x, double y) { return b.V(x, y, zt(x, y)); };
    auto B = [&](double x, double y) { return b.V(x, y, 0); };
    b.Quad(T(x0, y0), T(x1, y0), T(x1, y1), T(x0, y1));
    b.Quad(B(x0, y0), B(x0, y1), B(x1, y1), B(x1, y0));
    b.Quad(B(x0, y0), B(x1, y0), T(x1, y0), T(x0, y0));
    b.Quad(B(x1, y0), B(x1, y1), T(x1, y1), T(x1, y0));
    b.Quad(B(x1, y1), B(x0, y1), T(x0, y1), T(x1, y1));
    b.Quad(B(x0, y1), B(x0, y0), T(x0, y0), T(x0, y1));
  };
  for (const auto& bx : boxes) sbox(bx[0], bx[1], bx[2], bx[3]);
  return b.Mesh();
}

// Run the fold and grade the resolve against the GWN oracle + tol-invariance.
// The carriers below are MULTI-BOX soups (a cross-component coplanar overlap);
// under the non-fusion contract RemoveOverlaps3D decomposes them into separate
// clean components, so the fold ARRANGEMENT (Fix 1) is exercised through the
// ResolveComponentDirect hook, which treats the whole soup as one dirty
// component and runs the resolver on it directly (bypassing decompose + the
// gate).
static void ExpectFoldResolves(const char* tag, const MeshGL64& mesh,
                               double volLo, double volHi) {
  const Manifold::Impl in(mesh);
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold()) << tag;
  // The pure coplanar overlap is NOT flagged by the self-intersection gate
  // (R2(i)); the carrier reaches the fold via its coplanar cap clusters.
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.coplanarClusterFaces, 0) << tag << " must reach the fold";
  const double eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  auto run = [&](double e) { return ResolveComponentDirect(in, e); };

  const RegularizeResult r = run(eps);
  ASSERT_FALSE(r.fatal.has_value())
      << tag << " fold must resolve, not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.regularized, 1) << tag;

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError) << tag;
  EXPECT_FALSE(r.impl->IsSelfIntersecting())
      << tag << " output self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << tag << " must be one solid";
  const double vol = out.Volume();
  EXPECT_GT(vol, volLo) << tag << " volume below band";
  EXPECT_LT(vol, volHi) << tag << " volume above band";

  // INDEPENDENT GWN oracle: (w_soup>=1) == (w_out>0) at every unambiguous
  // point.
  std::mt19937 rng(0xC0FFEE);
  const Box bb = out.BoundingBox();
  const vec3 mn = bb.min - (bb.max - bb.min) * 0.05;
  const vec3 mx = bb.max + (bb.max - bb.min) * 0.05;
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<vec3> qs;
  for (int k = 0; k < 8000; ++k)
    qs.push_back(mn + (mx - mn) * vec3(U(rng), U(rng), U(rng)));
  const auto wOut = out.WindingNumber(qs);
  int disagree = 0, checked = 0;
  for (int k = 0; k < static_cast<int>(qs.size()) && disagree < 5; ++k) {
    const double g = GWNsoup(in, qs[k]);
    if (std::abs(g - std::round(g)) > 0.15) continue;  // near-surface skip
    ++checked;
    if ((std::lround(g) >= 1) != (wOut[k] > 0.5)) {
      ++disagree;
      ADD_FAILURE() << tag << " GWN membership disagreement at (" << qs[k].x
                    << "," << qs[k].y << "," << qs[k].z << ")";
    }
  }
  EXPECT_GT(checked, 1000) << tag << " oracle undersampled";

  // TOL-INVARIANCE: the retained topology is decided from input data, so the
  // enclosed volume is invariant to the weld radius.
  const RegularizeResult r2 = run(eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value())
      << tag << " tol-variant fatal: " << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * vol) << tag << " not tol-invariant";
}
}  // namespace coplanarfold

// Same-oriented mult-2 buried plug (B strictly inside A, caps coincident):
// {w_S>=1} = A (B is buried), volume = 12 exactly.  SlantPlug is TWO disjoint
// boxes with a CROSS-component coplanar overlap - under the non-fusion contract
// RemoveOverlaps3D passes both through (see the SlantPlug pass-through pin
// below); the FOLD ARRANGEMENT (Fix 1) is exercised on the whole soup via the
// direct hook.  A no-op / wrong fold fails the 11.9-12.1 band and the GWN
// oracle.
TEST(Overlap3, Regularize_CoplanarFold_Mult2_Resolves) {
  coplanarfold::ExpectFoldResolves("mult2", coplanarfold::SlantPlug(false),
                                   11.9, 12.1);
}

// Anti-oriented cancellation (B inverted inside A): the coincident caps cancel
// (mult 0 -> dropped), {w_S>=1} = A minus B, volume = 12 - 1.75 = 10.25.
TEST(Overlap3, Regularize_CoplanarFold_AntiCancellation_Resolves) {
  coplanarfold::ExpectFoldResolves("anti", coplanarfold::SlantPlug(true), 10.15,
                                   10.35);
}

// T-JUNCTION cap cluster (Fix-1 load-bearing): A + two adjacent plugs B,C
// sharing wall y=1.0, with C's cap bottom edge landing ON B's cap top edge
// interior -> T-junctions at (1,1) and (1.5,1) that split B's cap edge.  The
// rebuilt edgeSubdiv arrangement makes those splits; the retired hand-roll
// (proper crossings only) would leave B's cap edge unsplit and miswalk the
// cells.  B and C are buried in A -> {w_S>=1} = A, volume = 12.
TEST(Overlap3, Regularize_CoplanarFold_TJunction_Resolves) {
  coplanarfold::ExpectFoldResolves(
      "tjunc",
      coplanarfold::MultiPlug(
          {{0, 3, 0, 2}, {0.5, 2.5, 0.5, 1.0}, {1.0, 1.5, 1.0, 1.5}}),
      11.9, 12.1);
}

// mult-3 nested (A superset B superset C, all caps coincident, dyadic): the
// central cap region is covered three times; all buried, {w_S>=1} = A, the 4x4
// slanted box volume = 40.
TEST(Overlap3, Regularize_CoplanarFold_Mult3Nested_Resolves) {
  coplanarfold::ExpectFoldResolves(
      "mult3",
      coplanarfold::MultiPlug(
          {{0, 4, 0, 4}, {1, 3, 1, 3}, {1.5, 2.5, 1.5, 2.5}}),
      39.9, 40.1);
}

// Real carrier through the REAL entry: the openscad soup imports as a valid
// manifold with a WITHIN-component pure COPLANAR overlap
// (IsSelfIntersecting=0), so a single connectivity component routes DIRTY on
// the re-scoped coplanar gate (no cross-component merge - the non-fusion
// contract).  Post stage-6 SoS its non-coplanar ties DECIDE (it passes the
// exact-zero tie gate); then MANY of its truncated seams are cap-plane
// entanglement (a wall-wall seam endpoint on a coplanar cap cluster vertex)
// that the reg3d-ent junction completion now RECORDS - measured strictly
// narrower: the RecordSeams truncation count drops from 496 to 344 with the
// recovery active (152 seams completed).  The residue that REMAINS (~344) is a
// DIFFERENT wall: near-coplanar slivers (the GT7863-class thin seam sub-face)
// and cap- INTERIOR seam endpoints (the seam pierces a cap face interior, not
// an overlay vertex - would need the pierce injected as new fold input), both
// separate research axes.  So the component still fails NARROWER, a hard
// fail-closed = no output.  Still DirtyComponentUnresolved, a strictly narrower
// named reason than before, never a silent wrong resolve, no OOM
// (bbox-prefiltered on this 1442-face model).
TEST(Overlap3, Regularize_ExactZeroTie_Openscad_FailClosed) {
  std::filesystem::path file(__FILE__);
  std::ifstream fin(
      (file.parent_path() / "models" / "openscad-nonmanifold-crash.obj")
          .string());
  if (!fin.is_open()) GTEST_SKIP() << "model not found";
  const Manifold::Impl in(ReadOBJ(fin));
  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  EXPECT_GE(r.counters.dirty, 1) << "coplanar overlap must route to B";
  ASSERT_TRUE(r.fatal.has_value())
      << "the narrowed residue must fail closed, never silently resolve";
  EXPECT_EQ(*r.fatal, FatalReason::DirtyComponentUnresolved) << r.detail;
  // NARROWED: past the SoS gate AND past the cap-plane entanglement (152 seams
  // completed); the residue is the near-coplanar-sliver / cap-interior wall.
  EXPECT_NE(r.detail.find("degenerate incidence"), std::string::npos)
      << "residue must name the degenerate-seam wall: " << r.detail;
  EXPECT_FALSE(r.impl.has_value()) << "fail-closed yields no partial output";
}

// WITHIN-component coplanar carrier (the re-scoped GateComponent axis).
// BridgedCaps: box A [0,6]x[0,5]x[0,2] and box B [2,4]x[2,3]x[2,4] STACKED, B's
// bottom cap coincident with A's top over B's footprint (a DOUBLED internal
// wall = pure coplanar overlap).  A single containing quad over B keeps the
// caps detectable (they share NO verts) while B's walls stay interior to it.  A
// and B are joined into ONE connected component by a small solid L-ROD far from
// the caps (a hole in A's top annulus + a hole in B's +x wall, stitched through
// exterior) - a genus handle joining A-interior and B-interior, both w=1.  This
// makes the coplanar overlap INTERNAL to one component: IsSelfIntersecting
// misses it (R2(i)), so the WITHIN-component coplanar gate is what routes it
// DIRTY (the point of the re-scope).  The resolver then FAILS CLOSED on the
// bridge-junction SoS residue: the rod's faces are coplanar with the
// axis-normal wall/frame planes they connect to, an exact-zero orient3d that
// RecordSeams refuses (the single-global-SoS axis, PokedCube-class).  Building
// a genus handle is the ONLY way to make an internal coplanar overlap one
// connected 2-manifold (sharing the coincidence boundary is a non-manifold
// pinch or a transversal crossing; a buried plug's 3 winding levels can't be
// bridged consistently), and every such junction is on that SoS axis - so no
// within-component coplanar carrier RESOLVES through B today (a completeness
// counterexample, docs open list).  The FOLD ARRANGEMENT itself resolves the
// same defect on the direct hook (the CoplanarFold_* pins above).
static MeshGL64 BridgedCaps() {
  coplanarfold::MB b;
  auto V = [&](double x, double y, double z) { return b.V(x, y, z); };
  auto qz = [&](double x0, double x1, double y0, double y1, double z, bool up) {
    if (up)
      b.Quad(V(x0, y0, z), V(x1, y0, z), V(x1, y1, z), V(x0, y1, z));  // +z
    else
      b.Quad(V(x0, y0, z), V(x0, y1, z), V(x1, y1, z), V(x1, y0, z));  // -z
  };
  auto qx = [&](double y0, double y1, double z0, double z1, double x, bool pl) {
    if (pl)
      b.Quad(V(x, y0, z0), V(x, y1, z0), V(x, y1, z1), V(x, y0, z1));  // +x
    else
      b.Quad(V(x, y0, z0), V(x, y0, z1), V(x, y1, z1), V(x, y1, z0));  // -x
  };
  auto qy = [&](double x0, double x1, double z0, double z1, double y, bool pl) {
    if (pl)
      b.Quad(V(x0, y, z0), V(x0, y, z1), V(x1, y, z1), V(x1, y, z0));  // +y
    else
      b.Quad(V(x0, y, z0), V(x1, y, z0), V(x1, y, z1), V(x0, y, z1));  // -y
  };
  // frameZ / frameX: a rectangle in the z=const / x=const plane MINUS a
  // rectangular hole, tiled by 4 MITER trapezoids (outer edges are the full box
  // edges, so they match the plain-quad neighbours; the hole is interior).
  auto frameZ = [&](double X0, double X1, double Y0, double Y1, double hx0,
                    double hx1, double hy0, double hy1, double z, bool up) {
    auto q = [&](double ax, double ay, double bx, double by, double cx,
                 double cy, double dx, double dy) {
      if (up)
        b.Quad(V(ax, ay, z), V(bx, by, z), V(cx, cy, z), V(dx, dy, z));
      else
        b.Quad(V(ax, ay, z), V(dx, dy, z), V(cx, cy, z), V(bx, by, z));
    };
    q(X0, Y0, X1, Y0, hx1, hy0, hx0, hy0);
    q(X1, Y0, X1, Y1, hx1, hy1, hx1, hy0);
    q(X1, Y1, X0, Y1, hx0, hy1, hx1, hy1);
    q(X0, Y1, X0, Y0, hx0, hy0, hx0, hy1);
  };
  auto frameX = [&](double Y0, double Y1, double Z0, double Z1, double hy0,
                    double hy1, double hz0, double hz1, double x, bool pl) {
    auto q = [&](double ay, double az, double by, double bz, double cy,
                 double cz, double dy, double dz) {
      if (pl)
        b.Quad(V(x, ay, az), V(x, by, bz), V(x, cy, cz), V(x, dy, dz));
      else
        b.Quad(V(x, ay, az), V(x, dy, dz), V(x, cy, cz), V(x, by, bz));
    };
    q(Y0, Z0, Y1, Z0, hy1, hz0, hy0, hz0);
    q(Y1, Z0, Y1, Z1, hy1, hz1, hy1, hz0);
    q(Y1, Z1, Y0, Z1, hy0, hz1, hy1, hz1);
    q(Y0, Z1, Y0, Z0, hy0, hz0, hy0, hz1);
  };
  // A = [0,6]x[0,5]x[0,2].  A's top (z=2) is split at the grid lines
  // x{0,1,5,6}, y{0,1,4,5} so that:
  //  - a single CENTRAL quad [1,5]x[1,4] fully CONTAINS B [2,4]x[2,3], so B's
  //    walls are interior to it (no edge crosses a wall plane) and the coplanar
  //    overlap is DETECTED (the central quad shares no verts with B);
  //  - the surrounding FRAME cells are all bbox-disjoint (in x or y) from B's
  //    walls, so they are never compared with them - killing the spurious
  //    coplanar wall/annulus exact-zero ties that fail the fold closed.
  // The rod-hole sits in the +x frame cell [5,6]x[1,4] (miter frameZ so its
  // outer edges stay full and match the central quad + wall).
  {  // A bottom: 3x3 grid (lines x{0,1,5,6}, y{0,1,4,5}) to match the split
     // walls
    const double gx[] = {0, 1, 5, 6}, gy[] = {0, 1, 4, 5};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        qz(gx[i], gx[i + 1], gy[j], gy[j + 1], 0, /*up=*/false);
  }
  qz(1, 5, 1, 4, 2, /*up=*/true);  // A top CENTRAL quad (contains B)
  qz(0, 1, 0, 1, 2, /*up=*/true);  // A top frame cells (bbox-far from B)
  qz(1, 5, 0, 1, 2, /*up=*/true);
  qz(5, 6, 0, 1, 2, /*up=*/true);
  qz(0, 1, 1, 4, 2, /*up=*/true);
  qz(0, 1, 4, 5, 2, /*up=*/true);
  qz(1, 5, 4, 5, 2, /*up=*/true);
  qz(5, 6, 4, 5, 2, /*up=*/true);
  frameZ(5, 6, 1, 4, 5.2, 5.6, 2.35, 2.65, 2,
         /*up=*/true);  // +x frame+rod-hole
  // A walls, split at the frame grid lines so the top edges match.
  const double xs[][2] = {{0, 1}, {1, 5}, {5, 6}};
  for (auto& s : xs) {
    qy(s[0], s[1], 0, 2, 0, /*pl=*/false);  // -y wall segment
    qy(s[0], s[1], 0, 2, 5, /*pl=*/true);   // +y wall segment
  }
  const double ys[][2] = {{0, 1}, {1, 4}, {4, 5}};
  for (auto& s : ys) {
    qx(s[0], s[1], 0, 2, 0, /*pl=*/false);  // -x wall segment
    qx(s[0], s[1], 0, 2, 6, /*pl=*/true);   // +x wall segment
  }
  // B = [2,4]x[2,3]x[2,4].  Clean caps (2 tris each -> detectable coplanar
  // overlap with A's central quad).  +x wall (x=4) minus the rod-hole.
  qz(2, 4, 2, 3, 4, /*up=*/true);   // B top
  qz(2, 4, 2, 3, 2, /*up=*/false);  // B bottom (coincident cap)
  qy(2, 4, 2, 4, 2, /*pl=*/false);  // B -y
  qy(2, 4, 2, 4, 3, /*pl=*/true);   // B +y
  qx(2, 3, 2, 4, 2, /*pl=*/false);  // B -x
  frameX(2, 3, 2, 4, 2.35, 2.65, 2.8, 3, 4, /*pl=*/true);  // B +x (rod-hole)
  // Solid L-rod through EXTERIOR (above the cap, bbox-clear of the cap plane):
  // horizontal seg x[4,5.6] z[2.8,3] into B's +x wall, vertical seg x[5.2,5.6]
  // z[2,3] down into A's frame hole.  thickness y[2.35,2.65].
  const double ry0 = 2.35, ry1 = 2.65;
  qz(4, 5.2, ry0, ry1, 2.8, /*up=*/false);  // bottom of horizontal seg (-z)
  qx(ry0, ry1, 2, 2.8, 5.2, /*pl=*/false);  // inner step (-x)
  qx(ry0, ry1, 2, 3, 5.6, /*pl=*/true);     // right of vertical seg (+x)
  qz(4, 5.6, ry0, ry1, 3, /*up=*/true);     // top (+z)
  auto tri = [&](double ax, double az, double bx, double bz, double cx,
                 double cz, double y, bool plusY) {
    if (plusY) {
      b.tv.push_back(V(ax, y, az));
      b.tv.push_back(V(cx, y, cz));
      b.tv.push_back(V(bx, y, bz));
    } else {
      b.tv.push_back(V(ax, y, az));
      b.tv.push_back(V(bx, y, bz));
      b.tv.push_back(V(cx, y, cz));
    }
  };
  // y-faces = the L-hexagon, fanned from the reentrant corner (5.2,2.8).
  for (double y : {ry0, ry1}) {
    const bool pl = (y == ry1);
    tri(5.2, 2.8, 5.2, 2, 5.6, 2, y, pl);
    tri(5.2, 2.8, 5.6, 2, 5.6, 3, y, pl);
    tri(5.2, 2.8, 5.6, 3, 4, 3, y, pl);
    tri(5.2, 2.8, 4, 3, 4, 2.8, y, pl);
  }
  return b.Mesh();
}

// POINT-3 (owner adjudication): SlantPlug is TWO disjoint boxes with a
// CROSS-component coplanar overlap.  The non-fusion contract is UNIFORM -
// cross-component interaction (overlapping / touching / coplanar or
// transversal) is NEVER our job; fusion is the Boolean's.  So both boxes
// early-exit CLEAN and the operator emits them UNCHANGED (bitwise), Decompose
// == 2, volume 13.75. Mutation guard: a gate that wrongly merged the
// coplanar-overlapping components or routed either dirty would resolve the
// buried plug to 12.0 and fail the bitwise / count / volume checks below.
TEST(Overlap3, Regularize_SlantPlug_CrossComponent_PassThrough) {
  const Manifold::Impl in(coplanarfold::SlantPlug(false));
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  EXPECT_FALSE(in.IsSelfIntersecting())
      << "the coplanar overlap is cross-component and not self-intersecting";

  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value()) << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 2) << "two disjoint boxes stay separate";
  EXPECT_EQ(r.counters.clean, 2) << "both early-exit clean (no merge)";
  EXPECT_EQ(r.counters.dirty, 0)
      << "cross-component coplanar overlap is invisible"
         " to the per-component gate by design";
  EXPECT_EQ(r.counters.regularized, 0);

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_NEAR(out.Volume(), 13.75, 1e-9) << "un-regularized (A + protruding B)";
  EXPECT_EQ(out.Decompose().size(), 2u) << "still two solids";

  // BITWISE pass-through: the output vertex positions are exactly the input's
  // (concatenation of the two unchanged components, no fold, no weld shift).
  std::multiset<std::tuple<double, double, double>> inV, outV;
  for (const vec3& p : in.vertPos_) inV.emplace(p.x, p.y, p.z);
  const MeshGL64 og = GetMeshGLImpl<double, uint64_t>(*r.impl, -1);
  for (size_t i = 0; i + 2 < og.vertProperties.size(); i += 3)
    outV.emplace(og.vertProperties[i], og.vertProperties[i + 1],
                 og.vertProperties[i + 2]);
  EXPECT_EQ(inV, outV) << "pass-through must be bit-identical vertex positions";
}

// POINT-4 (owner adjudication): the WITHIN-component coplanar gate reached
// through the REAL RemoveOverlaps3D entry.  BridgedCaps is ONE connected
// 2-manifold (a stacked pair joined by a solid rod) whose only defect is an
// INTERNAL coplanar overlap - IsSelfIntersecting does NOT flag it (R2(i)), so
// the re-scoped GateComponent coplanar check is the ONLY thing that can route
// it to the resolver.  The resolver then fails closed on the bridge-junction
// single-global-SoS residue (the connection needed to make an internal coplanar
// overlap one component is itself on the unbuilt SoS axis; see the BridgedCaps
// header + the lane notebook).  Mutation guard: disable the
// DetectCoplanarClusters check in GateComponent and this component early-exits
// CLEAN (a silent wrong pass-through of an un-regularized doubled wall) - so
// the check is load-bearing.
// TARGET (a), docs/Regularize3D.md stage-6 SoS: BridgedCaps - a
// within-component coplanar overlap joined into ONE connected 2-manifold by a
// solid L-rod (a genus handle).  Its bridge-junction is the vertex-on-face /
// edge-on-edge exact-zero orient3d tie the SINGLE-GLOBAL SoS now DECIDES: the
// coplanar caps are consumed by the fold, and the rod junction (rod faces
// coplanar/ perpendicular with the wall & frame planes) is resolved by the SoS
// instead of failing closed.  RESOLVES oracle-true through the REAL
// RemoveOverlaps3D entry: {w_S>=1} = A[0,6]x[0,5]x[0,2] (60) +
// B[2,4]x[2,3]x[2,4] (4) + the L-rod
// (~0.192) = 64.192, one solid, the doubled cap interior.  Graded by the
// INDEPENDENT GWN solid-angle oracle (membership + volume band +
// tol-invariance). MUTATION-VERIFIED in-lane (reg3d-s6 notebook): disabling the
// SoS reverts this to the old fail-closed at the boundaryTouch gate (SoS is
// load-bearing); flipping the global perturbation direction still resolves
// oracle-true (the exact e^0 sign is unaffected, only sub-eps ties flip, not
// the {w>=1} topology).
// Compose two closed manifolds into one position-only (numProp=3) soup.
static Manifold::Impl TwoBoxSoup(const Manifold& a, const Manifold& b) {
  MeshGL64 out;
  out.numProp = 3;
  auto app = [&](const Manifold& m) {
    const MeshGL64 g = m.GetMeshGL64();
    const uint64_t base = out.NumVert();
    const size_t np = g.numProp;
    for (uint64_t v = 0; v < g.NumVert(); ++v)
      for (int k = 0; k < 3; ++k)
        out.vertProperties.push_back(g.vertProperties[v * np + k]);
    for (uint64_t t : g.triVerts) out.triVerts.push_back(t + base);
  };
  app(a);
  app(b);
  out.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(out);
}

// TARGET (e), docs/Regularize3D.md: ENTANGLED PRISMS.  Two axis-aligned bars
// crossing at right angles and sharing z[2,4] - a 3D plus.  The coincident
// z-caps over the crossing are a COPLANAR overlap (the FOLD's), and the
// transversal x/y walls cross at the vertical reentrant edges (x,y) in
// {(2,2),(2,4),(4,2),(4,4)}, z[2,4].  Each such wall-wall seam has ONE ordinary
// interior endpoint (on the wall diagonal) and ONE endpoint sitting on a cap
// cluster plane at an overlap-corner vertex OF THE FOLD'S IN-PLANE ARRANGEMENT
// (e.g. (2,2,2)).  The cap-plane endpoint was formerly DROPPED by
// edgeInClusterPlane's blanket suppression, truncating the seam to nPts==1 and
// failing closed ("degenerate incidence") - the COPLANAR/TRANSVERSAL
// ENTANGLEMENT wall.  The junction completion (reg3d-ent) gives that endpoint
// the fold arrangement's identity: RecordSeams records the cap-plane crossing
// (a genuine transversal endpoint on the cap fold) and dedups the symmetric
// double-pierce, so nPts==2, the seamed wall drops its interior span, and the
// folded cap + seamed walls weld shut at the reentrant corner.  RESOLVES
// oracle-true: {w_S>=1} is the exact union of the two bars (24 + 24 - 8 = 40).
// The junction is NOT a >2-sheet triple point (measured: every arrangement edge
// is exactly two walls; the coincident caps are the fold's, not a radial sheet)
// so NO radial rule / new predicate is needed (reg3d-radial STEP 1).
TEST(Overlap3, Regularize_ExactZeroTie_EntangledBars_Resolves) {
  const Manifold::Impl in =
      TwoBoxSoup(Manifold::Cube({6, 2, 2}).Translate({0, 2, 2}),
                 Manifold::Cube({2, 6, 2}).Translate({2, 0, 2}));
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.seamCount, 0) << "must reach a transversal crossing";
  EXPECT_GT(p.coplanarClusterFaces, 0) << "must reach the coplanar caps (fold)";

  const double eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  const RegularizeResult r = ResolveComponentDirect(in, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << "the cap-plane seam endpoint must resolve, not truncate: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.regularized, 1);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << "must be one solid";
  // Exact union of the two bars: 24 + 24 - 8 (the 2x2x2 overlap) = 40.
  const double vol = out.Volume();
  EXPECT_GT(vol, 39.999) << "volume below the exact-union band";
  EXPECT_LT(vol, 40.001) << "volume above the exact-union band";

  // INDEPENDENT GWN oracle: (w_soup>=1) == (w_out>0) at every unambiguous
  // point.
  const auto inTris = SoupTris(in);
  std::mt19937 rng(0xE47A0);
  const Box bb = out.BoundingBox();
  const vec3 mn = bb.min - (bb.max - bb.min) * 0.05;
  const vec3 mx = bb.max + (bb.max - bb.min) * 0.05;
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<vec3> qs;
  for (int k = 0; k < 20000; ++k)
    qs.push_back(mn + (mx - mn) * vec3(U(rng), U(rng), U(rng)));
  const auto wOut = out.WindingNumber(qs);
  int disagree = 0, checked = 0;
  for (int k = 0; k < static_cast<int>(qs.size()) && disagree < 5; ++k) {
    const double g = GWN(inTris, qs[k]);
    if (std::abs(g - std::round(g)) > 0.15) continue;  // near-surface skip
    ++checked;
    if ((std::lround(g) >= 1) != (wOut[k] > 0.5)) {
      ++disagree;
      ADD_FAILURE() << "GWN membership disagreement at (" << qs[k].x << ","
                    << qs[k].y << "," << qs[k].z << ")";
    }
  }
  EXPECT_GT(checked, 3000) << "oracle undersampled";

  // TOL-INVARIANCE: the retained topology is decided from input data.
  const RegularizeResult r2 = ResolveComponentDirect(in, eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << "tol-variant fatal: " << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * vol) << "not tol-invariant";
}

// ADVERSARIAL variant of the entanglement completion: the same two crossing
// bars but B ROTATED 50 degrees about z, so the coincident z-caps still cluster
// (rotation about z preserves the z=2/z=4 planes) yet the wall-wall seams cross
// at IRRATIONAL points and the recovered cap-plane endpoints are NOT bit-
// identical to the fold's overlay vertices (they differ sub-eps).  This
// stresses the closure: the seamed wall and the folded cap agree at the
// reentrant corner only through the emission WELD, not by exact construction.
// It RESOLVES oracle-true - volume matches the library's INDEPENDENT boolean
// union (a different algorithm), GWN membership agrees, one solid, not
// self-intersecting.
TEST(Overlap3, Regularize_ExactZeroTie_EntangledBarsRotated_Resolves) {
  const Manifold A = Manifold::Cube({6, 2, 2}).Translate({-3, -1, 2});
  const Manifold B = Manifold::Cube({6, 2, 2})
                         .Translate({-3, -1, 0})
                         .Rotate(0, 0, 50)
                         .Translate({0, 0, 2});
  const Manifold::Impl in = TwoBoxSoup(A, B);
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.seamCount, 0);
  EXPECT_GT(p.coplanarClusterFaces, 0);

  const double eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  const RegularizeResult r = ResolveComponentDirect(in, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << "rotated cap-plane junction must resolve: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << "must be one solid";
  // INDEPENDENT oracle #1: the library's boolean union (a different algorithm)
  // gives the exact union volume; the regularized output must match it.
  const double volUnion = (A + B).Volume();
  const double vol = out.Volume();
  EXPECT_NEAR(vol, volUnion, 1e-6 * volUnion) << "volume off the union";

  // INDEPENDENT oracle #2: GWN membership (solid-angle) at unambiguous points.
  const auto inTris = SoupTris(in);
  std::mt19937 rng(0xE47A1);
  const Box bb = out.BoundingBox();
  const vec3 mn = bb.min - (bb.max - bb.min) * 0.05;
  const vec3 mx = bb.max + (bb.max - bb.min) * 0.05;
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<vec3> qs;
  for (int k = 0; k < 20000; ++k)
    qs.push_back(mn + (mx - mn) * vec3(U(rng), U(rng), U(rng)));
  const auto wOut = out.WindingNumber(qs);
  int disagree = 0, checked = 0;
  for (int k = 0; k < static_cast<int>(qs.size()) && disagree < 5; ++k) {
    const double g = GWN(inTris, qs[k]);
    if (std::abs(g - std::round(g)) > 0.15) continue;  // near-surface skip
    ++checked;
    if ((std::lround(g) >= 1) != (wOut[k] > 0.5)) {
      ++disagree;
      ADD_FAILURE() << "GWN membership disagreement at (" << qs[k].x << ","
                    << qs[k].y << "," << qs[k].z << ")";
    }
  }
  EXPECT_GT(checked, 3000) << "oracle undersampled";

  // TOL-INVARIANCE.
  const RegularizeResult r2 = ResolveComponentDirect(in, eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << "tol-variant fatal: " << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * vol) << "not tol-invariant";
}

// TARGET (e) variant: the same two bars OFFSET in z (no coincident caps - the
// fold is NOT reached, coplanarClusterFaces == 0), walls crossing edge-on-edge
// at exact ties (the SoS axis, decided).  The seamed faces build; the
// clean-face winding classify then GRAZED every seed at the triangle centroid
// on this axis-aligned integer geometry - the old flood, probing one centroid
// per patch, fail-closed on the COMPONENT-LOCAL SEED graze.  The per-face clean
// rule (reg3d-arr) re-probes OTHER interior points of the same clean triangle,
// which sample the SAME (constant) winding cell above an uncrossed face,
// dodging the graze SOUNDLY (every emitted face's winding is directly measured,
// no uniformity/borrow assumption).  So this now RESOLVES oracle-true: the
// {w_S>=1} boundary is the exact union of the two bars.  This CLOSES the
// component-local seed-policy open for this carrier (a bounded
// decision-completion: sample the constant cell, not a global seed policy).
TEST(Overlap3, Regularize_ExactZeroTie_BarsCrossZ_Resolves) {
  const Manifold::Impl in =
      TwoBoxSoup(Manifold::Cube({6, 2, 4}).Translate({0, 2, 0}),
                 Manifold::Cube({2, 6, 2}).Translate({2, 0, 1}));
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.seamCount, 0) << "must reach transversal crossings";
  EXPECT_EQ(p.coplanarClusterFaces, 0) << "no coplanar caps (fold not reached)";

  const double eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  const RegularizeResult r = ResolveComponentDirect(in, eps);
  ASSERT_FALSE(r.fatal.has_value()) << "the per-face clean rule must resolve "
                                       "the seed-graze, not fail closed: "
                                    << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.regularized, 1);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << "must be one solid";
  // Exact union of the two bars: 48 + 24 - 8 (the 2x2x2 overlap) = 64.
  const double vol = out.Volume();
  EXPECT_GT(vol, 63.999) << "volume below the exact-union band";
  EXPECT_LT(vol, 64.001) << "volume above the exact-union band";

  // INDEPENDENT GWN oracle: (w_soup>=1) == (w_out>0) at every unambiguous
  // point.
  const auto inTris = SoupTris(in);
  std::mt19937 rng(0xBC0FF);
  const Box bb = out.BoundingBox();
  const vec3 mn = bb.min - (bb.max - bb.min) * 0.05;
  const vec3 mx = bb.max + (bb.max - bb.min) * 0.05;
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<vec3> qs;
  for (int k = 0; k < 20000; ++k)
    qs.push_back(mn + (mx - mn) * vec3(U(rng), U(rng), U(rng)));
  const auto wOut = out.WindingNumber(qs);
  int disagree = 0, checked = 0;
  for (int k = 0; k < static_cast<int>(qs.size()) && disagree < 5; ++k) {
    const double g = GWN(inTris, qs[k]);
    if (std::abs(g - std::round(g)) > 0.15) continue;  // near-surface skip
    ++checked;
    if ((std::lround(g) >= 1) != (wOut[k] > 0.5)) {
      ++disagree;
      ADD_FAILURE() << "GWN membership disagreement at (" << qs[k].x << ","
                    << qs[k].y << "," << qs[k].z << ")";
    }
  }
  EXPECT_GT(checked, 3000) << "oracle undersampled";

  // TOL-INVARIANCE: the retained topology is decided from input data.
  const RegularizeResult r2 = ResolveComponentDirect(in, eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << "tol-variant fatal: " << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * vol) << "not tol-invariant";
}

static Manifold::Impl NBoxSoup(const std::vector<Manifold>& ms) {
  MeshGL64 out;
  out.numProp = 3;
  for (const Manifold& m : ms) {
    const MeshGL64 g = m.GetMeshGL64();
    const uint64_t base = out.NumVert();
    const size_t np = g.numProp;
    for (uint64_t v = 0; v < g.NumVert(); ++v)
      for (int k = 0; k < 3; ++k)
        out.vertProperties.push_back(g.vertProperties[v * np + k]);
    for (uint64_t t : g.triVerts) out.triVerts.push_back(t + base);
  }
  out.runOriginalID.push_back(Manifold::ReserveIDs(1));
  return Manifold::Impl(out);
}

// COPLANAR/TRANSVERSAL CAP-SEAMING ENTANGLEMENT (census Cluster 2a, ledger
// C-2a). A base coplanar cap CLUSTER at z=2 - box A [0,6]x[0,6]x[0,2] top cap
// coincident with box B [1,5]x[1,5]x[0,2] top cap (a mult-2 overlap the fold
// owns) - is PIERCED by a small box C [2,4]x[2,4] straddling z[1,3].  C's four
// walls cross the z=2 plane strictly INSIDE the cap faces, so the coplanar cap
// faces are TRANSVERSALLY SEAMED (A.seamed set on a cluster face) and the
// wall-cap seam endpoints land in the cap INTERIOR (C's footprint corners /
// cap-triangulation diagonals - NOT an overlay vertex the fold owns).
//
// This is the C-2a target the census named "cap-INTERIOR pierce injection" and
// classified a bounded decision-completion.  This lane's per-pair measurement
// (reg3d-c2a) REFUTES the bounded classification: a wall-wall seam endpoint on
// a cap-cluster plane is where two walls' bottom edges (both in the plane)
// meet, and a wall bottom edge in the cluster plane is a cap BOUNDARY edge - so
// their meeting is a fold overlay vertex (capvert), never strictly
// cap-interior.  A cap-INTERIOR endpoint therefore only arises when a wall
// PIERCES a cap FACE, and that pierce SEAMS the cap = the coplanar/transversal
// entanglement the fold declines BY DESIGN (a transversal seam splits a fold
// cell with a 3D winding jump the coplanar mult does not carry).  So the F11
// seam truncation here is only a SYMPTOM; the fold-decline (F3) sits underneath
// (verified in-lane: measured seamedCluster>0 with the seams recorded).
// Injecting the pierce and arranging the cap around the transversal seam =
// unifying the coplanar fold with the transversal seam machinery + classifying
// each split sub-cell by the real 3D coupled winding = RESEARCH-GRADE (the
// Cluster-1 coordinated-emission wall), not a bounded endpoint injection.  So
// this FAILS CLOSED (narrower than a silent wrong resolve), documenting the
// wall.  The oracle exists (the union is a valid solid, volume A + C-above-cap
// = 72 + 4 = 76) - the geometry is resolvable in principle, just not by a
// bounded completion.
//
// CONTRAST (the delta is the PIERCE): the same A+B cluster with C sitting ON
// the cap (bottom at z=2, no transversal pierce - a buried plug of the
// coplanar-fold family) RESOLVES oracle-true (volume 72 + 8 = 80).  Piercing vs
// resting on the cap is the entire boundary between the fold family (resolves)
// and the cap-seaming entanglement (research-grade fail-closed).
TEST(Overlap3, Regularize_CapSeamingEntanglement_CapInteriorPierce_FailClosed) {
  const Manifold A = Manifold::Cube({6, 6, 2});                       // [0,6]^2
  const Manifold B = Manifold::Cube({4, 4, 2}).Translate({1, 1, 0});  // [1,5]^2
  const Manifold Cpierce =
      Manifold::Cube({2, 2, 2}).Translate({2, 2, 1});  // [2,4]^2 z[1,3] PIERCES
  const Manifold::Impl in = NBoxSoup({A, B, Cpierce});
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.coplanarClusterFaces, 0) << "must reach the coplanar cap cluster";
  EXPECT_GT(p.seamCount, 0) << "the pierce must transversally seam the cap";

  const double eps = EpsilonFromScale(in.bBox_.Scale(), 1000);
  const RegularizeResult r = ResolveComponentDirect(in, eps);
  // FAIL CLOSED: the cap-seaming entanglement is research-grade, never a silent
  // wrong resolve.  Today the RecordSeams truncation (F11 "degenerate
  // incidence") is the first gate; the fold-decline (F3) sits underneath it.
  ASSERT_TRUE(r.fatal.has_value())
      << "the cap-seaming entanglement must fail closed, not resolve: "
      << r.detail;
  EXPECT_EQ(*r.fatal, FatalReason::DirtyComponentUnresolved) << r.detail;
  EXPECT_FALSE(r.impl.has_value()) << "fail-closed yields no partial output";

  // CONTRAST: C sitting ON the cap (no pierce) is the resolving buried-plug
  // case.
  const Manifold Crest =
      Manifold::Cube({2, 2, 2}).Translate({2, 2, 2});  // [2,4]^2 z[2,4] RESTS
  const Manifold::Impl inRest = NBoxSoup({A, B, Crest});
  const RegularizeResult rRest = ResolveComponentDirect(inRest, eps);
  ASSERT_FALSE(rRest.fatal.has_value())
      << "the non-piercing plug (fold family) must resolve: " << rRest.detail;
  ASSERT_TRUE(rRest.impl.has_value());
  const Manifold outRest(GetMeshGLImpl<double, uint64_t>(*rRest.impl, -1));
  EXPECT_EQ(outRest.Decompose().size(), 1u) << "rest: must be one solid";
  EXPECT_FALSE(rRest.impl->IsSelfIntersecting()) << "rest: self-intersects";
  const double volRest = outRest.Volume();
  EXPECT_GT(volRest, 79.99) << "rest: volume below the union band";
  EXPECT_LT(volRest, 80.01) << "rest: volume above the union band";
}

TEST(Overlap3, Regularize_WithinComponentCoplanar_BridgedCaps_Resolves) {
  const Manifold::Impl in(BridgedCaps());
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold())
      << "a valid single connected 2-manifold";
  EXPECT_EQ(
      Manifold(GetMeshGLImpl<double, uint64_t>(in, -1)).Decompose().size(), 1u)
      << "ONE connected component (the rod joins A and B)";
  EXPECT_FALSE(in.IsSelfIntersecting())
      << "the internal coplanar overlap is NOT self-intersecting (R2(i))";
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(97.1, 33.7, 51.3));
  EXPECT_GT(p.coplanarClusterFaces, 0)
      << "the doubled cap is a genuine within-component coplanar overlap";

  const double eps = ImplEps(in);
  const RegularizeResult r = RemoveOverlaps3D(in, eps);
  ASSERT_FALSE(r.fatal.has_value())
      << "the SoS must resolve the bridge junction, not fail closed: "
      << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 1);
  EXPECT_EQ(r.counters.clean, 0)
      << "the within-component coplanar gate routes it DIRTY";
  EXPECT_EQ(r.counters.dirty, 1);
  EXPECT_EQ(r.counters.regularized, 1);

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(r.impl->IsSelfIntersecting()) << "output self-intersects";
  EXPECT_EQ(out.Decompose().size(), 1u) << "must be one solid";
  const double vol = out.Volume();
  EXPECT_GT(vol, 64.0) << "volume below band";
  EXPECT_LT(vol, 64.4) << "volume above band";

  // INDEPENDENT GWN oracle: (w_soup>=1) == (w_out>0) at every unambiguous
  // point.
  const auto inTris = SoupTris(in);
  std::mt19937 rng(0xB1D9E);
  const Box bb = out.BoundingBox();
  const vec3 mn = bb.min - (bb.max - bb.min) * 0.05;
  const vec3 mx = bb.max + (bb.max - bb.min) * 0.05;
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<vec3> qs;
  for (int k = 0; k < 12000; ++k)
    qs.push_back(mn + (mx - mn) * vec3(U(rng), U(rng), U(rng)));
  const auto wOut = out.WindingNumber(qs);
  int disagree = 0, checked = 0;
  for (int k = 0; k < static_cast<int>(qs.size()) && disagree < 5; ++k) {
    const double g = GWN(inTris, qs[k]);
    if (std::abs(g - std::round(g)) > 0.15) continue;  // near-surface skip
    ++checked;
    if ((std::lround(g) >= 1) != (wOut[k] > 0.5)) {
      ++disagree;
      ADD_FAILURE() << "GWN membership disagreement at (" << qs[k].x << ","
                    << qs[k].y << "," << qs[k].z << ")";
    }
  }
  EXPECT_GT(checked, 2000) << "oracle undersampled";

  // TOL-INVARIANCE: the retained topology is decided from input data.
  const RegularizeResult r2 = RemoveOverlaps3D(in, eps * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << "tol-variant fatal: " << r2.detail;
  const Manifold out2(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1));
  EXPECT_NEAR(vol, out2.Volume(), 1e-6 * vol) << "not tol-invariant";
}

// TARGET (c), docs/Regularize3D.md stage-5 + non-fusion contract: the hull
// (body + mask) carries the corpus's real near-coplanar overlap geometry - two
// large flat facets within eps of coplanar (research memo: 14 decidable-thin
// near-coplanar pairs).  But that overlap is CROSS-component (between the body
// and mask solids), so RemoveOverlaps3D decomposes into separate clean
// components and passes them through UNCHANGED under the uniform non-fusion
// contract - the per-component near-coplanar widen never sees a cross-component
// cluster (and no component is dirty).  A regression guard: a gate that wrongly
// fused or routed the cross-component overlap dirty would change the component
// count / output.
#ifndef MANIFOLD_NO_FILESYSTEM
TEST(Overlap3, Regularize_Hull_CrossComponent_PassThrough) {
  std::filesystem::path file(__FILE__);
  auto modelDir = file.parent_path() / "models";
  std::ifstream fBody((modelDir / "hull-body.obj").string());
  std::ifstream fMask((modelDir / "hull-mask.obj").string());
  if (!fBody.is_open() || !fMask.is_open()) GTEST_SKIP() << "hull model absent";
  const Manifold::Impl impl =
      ComposeImpl(Manifold::ReadOBJ(fBody), Manifold::ReadOBJ(fMask));
  const RegularizeResult r = RemoveOverlaps3D(impl, ImplEps(impl));
  ASSERT_FALSE(r.fatal.has_value()) << r.detail;
  EXPECT_GT(r.counters.components, 1) << "body + mask are distinct components";
  EXPECT_EQ(r.counters.clean, r.counters.components)
      << "every component is clean; the near-coplanar overlap is "
         "cross-component";
  EXPECT_EQ(r.counters.dirty, 0);
  EXPECT_EQ(r.counters.regularized, 0);
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.impl->NumTri(), impl.NumTri()) << "cross-component pass-through";
  // The output remains self-intersecting BY CONTRACT: the body/mask overlap is
  // cross-component, which this operator never resolves (fusion is the
  // Boolean's job) - the near-coplanar widen is within-component only.
}
#endif

// ===========================================================================
// STAGE-5: NEAR-COPLANAR widen + global-planarity guard (reg3d-s5,
// docs/Regularize3D.md; reg3d-nearcoplanar-research candidate (a)).  A face
// pair within eps of coplanar but ABOVE the filter's ~1-ULP bound is DECIDABLE
// non-coplanar, so the exact fold does not cluster it; enumerated transversally
// its sub-eps-thin cell double-rounds to a sliver ("unresolvable sheet
// contact").  The input-side planarization (SnapNearCoplanarClusters) snaps
// such a cluster onto its fitted plane so B re-derives the arrangement from an
// exactly-coplanar input (the landed exact fold verbatim).  A curved chain (fit
// deviation > eps) FAILS the global-planarity guard, named.
//
// NOTE on reachability: a PURE near-coplanar overlap (deviation < 2eps) is
// invisible to IsSelfIntersecting (its 2*eps normal-nudge always separates two
// near-coplanar faces), so - exactly like the exact CoplanarFold_* carriers -
// the fold is exercised through the ResolveComponentDirect hook (the resolver
// on the soup as one dirty component), not the gate.
// ===========================================================================

// SlantPlug with B's TOP cap tilted a sub-eps `delta` off A's slant plane along
// +x (one-sided, `cross`=false) or straddling it (`cross`=true) -
// NEAR-coplanar, filter-decidable non-coplanar.  The z=0 bottom caps stay
// EXACTLY coplanar.
static MeshGL64 NearSlantPlug(double delta, bool cross, bool flipB) {
  MeshGL64 m = coplanarfold::SlantPlug(flipB);
  auto zt = [](double x, double y) { return 1.0 + 0.5 * x + 0.25 * y; };
  for (uint64_t v = 0; v < m.NumVert(); ++v) {
    const double x = m.vertProperties[v * 3], y = m.vertProperties[v * 3 + 1],
                 z = m.vertProperties[v * 3 + 2];
    const bool bTop = (x == 0.5 || x == 1.5) && (y == 0.5 || y == 1.5) &&
                      std::abs(z - zt(x, y)) < 1e-9;
    if (bTop)
      m.vertProperties[v * 3 + 2] = z + delta * (cross ? x - 1.0 : x - 0.5);
  }
  return m;
}

// K nested unit-thick boxes (each footprint strictly inside the larger) whose
// TOP caps fan through z=0 with tilt k*dm about the y-axis (top corner
// z=k*dm*x). Consecutive tops overlap in footprint and share a z-band
// (bbox-detected) with gap ~ dm*s < eps so they chain in the widened
// union-find; the fan's tilt spread makes the cluster deviate > eps from any
// plane, so it fails the guard.
static MeshGL64 CurvedChainSoup(double dm, int K) {
  coplanarfold::MB b;
  const double H = 1.0;
  for (int k = 0; k < K; ++k) {
    const double s = 3.0 - 0.3 * k;
    auto V = [&](double x, double y, double z) { return b.V(x, y, z); };
    auto T = [&](double x, double y) { return V(x, y, k * dm * x); };
    auto Bo = [&](double x, double y) { return V(x, y, -H); };
    b.Quad(T(-s, -s), T(s, -s), T(s, s), T(-s, s));      // top (tilted)
    b.Quad(Bo(-s, -s), Bo(-s, s), Bo(s, s), Bo(s, -s));  // bottom -z
    b.Quad(Bo(-s, -s), Bo(s, -s), T(s, -s), T(-s, -s));  // -y wall
    b.Quad(Bo(s, -s), Bo(s, s), T(s, s), T(s, -s));      // +x wall
    b.Quad(Bo(s, s), Bo(-s, s), T(-s, s), T(s, s));      // +y wall
    b.Quad(Bo(-s, s), Bo(-s, -s), T(-s, -s), T(-s, s));  // -x wall
  }
  return b.Mesh();
}

// Base eps for the SlantPlug family (bbox scale ~3).
static double SlantPlugEps() {
  return EpsilonFromScale(
      Manifold::Impl(coplanarfold::SlantPlug(false)).bBox_.Scale(), 1000);
}

// TARGET (a): same-oriented mult-2 buried plug whose TOP caps are NEAR-coplanar
// (tilted 0.3 eps).  Without the widen this fails closed at the sliver
// emission; the snap folds it to {w_S>=1} = A, volume 12 (identical to the
// exact carrier). coplanarClusterFaces == 4 pins that only the z=0 bottoms
// cluster EXACTLY, so the resolve to 12 requires the NEAR top caps to be folded
// by stage-5.
TEST(Overlap3, Regularize_NearCoplanarFold_Mult2_Resolves) {
  const double delta =
      0.3 * SlantPlugEps();  // near band; tol-invariant to eps/2
  const MeshGL64 mesh = NearSlantPlug(delta, /*cross=*/false, /*flipB=*/false);
  const Manifold::Impl in(mesh);
  const ComponentEnumProbe p = EnumerateComponent_Probe(
      in, {}, in.bBox_.Center() + vec3(9.71, 3.37, 5.13));
  EXPECT_EQ(p.coplanarClusterFaces, 4)
      << "only the exact z=0 bottom caps cluster; the tilted top caps are the "
         "near band (would be 8 if the tops were exactly coplanar)";
  coplanarfold::ExpectFoldResolves("nearmult2", mesh, 11.9, 12.1);
}

// TARGET (a) variant: anti-oriented cancellation with NEAR-coplanar caps (B
// inverted, top tilted).  The snapped coincident caps cancel (mult 0), {w_S>=1}
// = A minus B, volume 10.25 - the mult algebra stays the stage-4 form.
TEST(Overlap3, Regularize_NearCoplanarFold_AntiCancellation_Resolves) {
  const double delta = 0.3 * SlantPlugEps();
  coplanarfold::ExpectFoldResolves(
      "nearanti", NearSlantPlug(delta, /*cross=*/false, /*flipB=*/true), 10.15,
      10.35);
}

// TARGET (b): a CURVED near-coplanar chain (tilted-fan nested boxes) whose
// consecutive tops are within eps (they chain) but whose global fit deviation
// exceeds eps.  The global-planarity guard REFUSES it with a distinct named
// reason - a strictly-narrower fail-closed than today's blanket refusal, never
// a silent wrong resolve.  (dm chosen so the consecutive gap ~0.8 eps < eps but
// the fan deviation ~1.6 eps > eps.)
TEST(Overlap3, Regularize_NearCoplanarChain_GuardFailClosed) {
  const int K = 5;
  const double s0 = 3.0;
  const Manifold::Impl probe(CurvedChainSoup(1e-13, K));
  const double eps = EpsilonFromScale(probe.bBox_.Scale(), 1000);
  const Manifold::Impl in(CurvedChainSoup(0.8 * eps / s0, K));
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold());
  const RegularizeResult r = ResolveComponentDirect(in, eps);
  ASSERT_TRUE(r.fatal.has_value())
      << "a curved near-coplanar chain must fail closed, never fold to a wrong "
         "plane";
  EXPECT_EQ(*r.fatal, FatalReason::DirtyComponentUnresolved) << r.detail;
  EXPECT_NE(r.detail.find("global-planarity guard"), std::string::npos)
      << "residue must name the guard: " << r.detail;
  EXPECT_FALSE(r.impl.has_value());
}

// ===========================================================================
// Corpus re-pin: the sole-impl (RemoveOverlaps3D / RemoveOverlaps3D) contract
// on the file fixtures.  The v3 sweep Corpus_* pins retired with the sweep;
// these MEASURE and pin what the regularization operator does today.  Under the
// non-fusion contract an overlapping pair / multi-shell input decomposes into
// per-shell components, each gated independently: a clean shell early-exits
// (bitwise pass-through), a self-intersecting or within-component-coplanar
// shell routes to the resolver.  Figures are the measured dispatch (reg3d-sole
// notebook); a mutation that fused components or misrouted a clean shell reds
// the count / vertex-set check.
// ===========================================================================
#ifndef MANIFOLD_NO_FILESYSTEM

// Clean pass-through: every component gates clean, so the output is the input's
// components concatenated unchanged - bitwise-identical vertex positions, no
// fold, no weld shift.  `nComp` is the measured Decompose count.
static void ExpectCorpusCleanPassThrough(const char* tag,
                                         const Manifold::Impl& in, int nComp) {
  ASSERT_TRUE(in.IsManifold() && in.Is2Manifold())
      << tag << " fixture must be a valid 2-manifold";
  const RegularizeResult r = RemoveOverlaps3D(in, ImplEps(in));
  ASSERT_FALSE(r.fatal.has_value())
      << tag << " must not fail closed: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, nComp) << tag << " decompose count";
  EXPECT_EQ(r.counters.clean, nComp) << tag << " every component gates clean";
  EXPECT_EQ(r.counters.dirty, 0) << tag << " no within-component defect";
  EXPECT_EQ(r.counters.regularized, 0);
  EXPECT_EQ(r.counters.failClosed, 0);

  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError) << tag << " output invalid";
  EXPECT_EQ(out.Decompose().size(), static_cast<size_t>(nComp))
      << tag << " components stay separate (no fusion)";

  // BITWISE pass-through: the output vertex positions are exactly the input's
  // (concatenation of the unchanged components; no fold, no weld shift).  A
  // gate that misrouted a clean shell to B, or fused components, would perturb
  // them.
  std::multiset<std::tuple<double, double, double>> inV, outV;
  for (const vec3& p : in.vertPos_) inV.emplace(p.x, p.y, p.z);
  const MeshGL64 og = GetMeshGLImpl<double, uint64_t>(*r.impl, -1);
  for (size_t i = 0; i + 2 < og.vertProperties.size(); i += 3)
    outV.emplace(og.vertProperties[i], og.vertProperties[i + 1],
                 og.vertProperties[i + 2]);
  EXPECT_EQ(inV, outV) << tag << " pass-through must be bit-identical verts";
}

static std::optional<Manifold::Impl> LoadCorpusPair(const char* l,
                                                    const char* rr) {
  std::filesystem::path f(__FILE__);
  std::ifstream fl((f.parent_path() / "models" / l).string());
  std::ifstream fr((f.parent_path() / "models" / rr).string());
  if (!fl.is_open() || !fr.is_open()) return std::nullopt;
  return ComposeImpl(Manifold::ReadOBJ(fl), Manifold::ReadOBJ(fr));
}
static std::optional<Manifold::Impl> LoadCorpusSingle(const char* n) {
  std::filesystem::path f(__FILE__);
  std::ifstream fin((f.parent_path() / "models" / n).string());
  if (!fin.is_open()) return std::nullopt;
  return Manifold::Impl(ReadOBJ(fin));
}

// Cray pair: joins into ONE connected component that gates clean.
TEST(Overlap3, Corpus_Cray_CleanPassThrough) {
  const auto in = LoadCorpusPair("Cray_left.obj", "Cray_right.obj");
  if (!in) GTEST_SKIP() << "model not found";
  ExpectCorpusCleanPassThrough("Cray", *in, 1);
}

// Havocglass8 pair: two overlapping shells.  Cross-component overlap is out of
// scope (non-fusion), each shell gates clean -> pass-through.
TEST(Overlap3, Corpus_Havocglass8_CleanPassThrough) {
  const auto in =
      LoadCorpusPair("Havocglass8_left.obj", "Havocglass8_right.obj");
  if (!in) GTEST_SKIP() << "model not found";
  ExpectCorpusCleanPassThrough("Havocglass8", *in, 2);
}

// Offset meshes: multi-shell solids whose shells are each a clean solid.
TEST(Overlap3, Corpus_Offsets_CleanPassThrough) {
  const struct {
    const char* name;
    int nComp;
  } cases[] = {{"Offset1.obj", 39},
               {"Offset2.obj", 45},
               {"Offset3.obj", 1},
               {"Offset4.obj", 1}};
  for (const auto& c : cases) {
    const auto in = LoadCorpusSingle(c.name);
    if (!in) GTEST_SKIP() << c.name << " not found";
    ExpectCorpusCleanPassThrough(c.name, *in, c.nComp);
  }
}

// GenericTwin7081: the pair decomposes into 13 components; 11 gate clean and 2
// carry a within-component defect - two near-flat SHELLS with a 0.002deg
// near-tangent self-overlap.  History: census -> Cluster-2 seam sub-face;
// reg3d- c2b -> winding-probe filter-precision (closed by the reg3d-c2bx
// exact-classify escalation); then the pre-existing Cluster-1 emission wall
// (unbalanced fans from crossings the shares-vertex broadphase skip DROPPED).
// reg3d-wjump CLOSES it: the shares-vertex genuine-crossing recovery records
// those crossings, so both shells complete their arrangement and RESOLVE.
// Oracle: the near-tangent overlap is measure-~0, so each shell's {w>=1} volume
// is preserved to its input signed volume (a GWN grid oracle is unusable at
// these near-flat sheets at scale ~20000 - one cell exceeds the whole solid);
// the resolve is non-self-intersecting, a valid manifold, and tol-invariant.
// MUTATION-VERIFIED: reverting the recovery reopens the fans ->
// NonManifoldEmission.  HEAVY (~45s); run under the corpus resource cap (ulimit
// -v 4000000; timeout 900).
TEST(Overlap3, Corpus_GenericTwin7081_Resolves) {
  const auto in = LoadCorpusPair("Generic_Twin_7081.1.t0_left.obj",
                                 "Generic_Twin_7081.1.t0_right.obj");
  if (!in) GTEST_SKIP() << "model not found";
  const double inVol =
      Manifold(GetMeshGLImpl<double, uint64_t>(*in, -1)).Volume();
  const RegularizeResult r = RemoveOverlaps3D(*in, ImplEps(*in));
  ASSERT_FALSE(r.fatal.has_value()) << "GT7081 must resolve: " << r.detail;
  ASSERT_TRUE(r.impl.has_value());
  EXPECT_EQ(r.counters.components, 13) << "decompose count";
  EXPECT_EQ(r.counters.clean, 11);
  EXPECT_EQ(r.counters.dirty, 2)
      << "two shells carry a within-component defect";
  EXPECT_EQ(r.counters.regularized, 2);
  EXPECT_EQ(r.counters.failClosed, 0);
  const Manifold out(GetMeshGLImpl<double, uint64_t>(*r.impl, -1));
  EXPECT_EQ(out.Status(), Manifold::Error::NoError) << "output not a manifold";
  // The composed whole self-intersects by design (the two twins overlap in
  // space - a CROSS-component overlap the non-fusion contract keeps separate,
  // never this operator's job), so check per-component: every resolved shell is
  // individually non-self-intersecting (the within-component defect removed).
  for (const Manifold& c : out.Decompose())
    EXPECT_FALSE(Manifold::Impl(c.GetMeshGL64()).IsSelfIntersecting())
        << "a resolved component still self-intersects";
  // Measure-~0 near-tangent overlap: {w>=1} volume preserved.
  EXPECT_NEAR(out.Volume(), inVol, 1e-3 * std::abs(inVol))
      << "resolved volume must preserve the near-tangent soup volume";
  const RegularizeResult r2 = RemoveOverlaps3D(*in, ImplEps(*in) * 0.5);
  ASSERT_TRUE(r2.impl.has_value()) << r2.detail;
  EXPECT_NEAR(out.Volume(),
              Manifold(GetMeshGLImpl<double, uint64_t>(*r2.impl, -1)).Volume(),
              1e-6 * std::abs(out.Volume()) + 1e-9)
      << "enclosed volume is not tol-invariant";
}
#endif
