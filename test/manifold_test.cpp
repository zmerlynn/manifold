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

#include "manifold/manifold.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <thread>

#include "../src/cross_section/boolean2/predicates.h"  // EpsilonFromScale
#include "../src/execution_impl.h"
#include "../src/impl.h"
#include "../src/overlap_removal.h"  // for RemoveOverlaps (explicit eps)
#include "../src/overlap_removal_internal.h"
#include "manifold/cross_section.h"
#include "test.h"

namespace {

using namespace manifold;

template <typename T>
int NumUnique(const std::vector<T>& in) {
  std::set<int> unique;
  for (const T& v : in) {
    unique.emplace(v);
  }
  return unique.size();
}

// A disjoint Union of a CubeSTL (no hasNormals) and a Sphere with
// CalculateNormals. The result has mixed meshIDs - some with hasNormals,
// some without - the exact shape the per-meshID handling in
// Impl::Transform / Compose / CreateProperties was added to support.
Manifold MixedNormalsCubePlusSphere(double sphereRadius = 5.0) {
  MeshGL cubeGL = CubeSTL();
  cubeGL.Merge();
  return Manifold(cubeGL).Translate({20, 0, 0}) +
         Manifold::Sphere(sphereRadius, 32).CalculateNormals(0);
}

// Count verts on the sphere surface (|pos| ~ sphereRadius) whose stored
// normal at slot 3..5+offset aligns with `expected(pos)` (dot > 0.9).
// Returns (good, bad).
template <typename ExpectedFn>
std::pair<int, int> CountSphereNormalAlignment(const MeshGL& gl,
                                               double sphereRadius,
                                               int normalSlotOffset,
                                               ExpectedFn expected) {
  int good = 0, bad = 0;
  for (size_t v = 0; v < gl.NumVert(); ++v) {
    const vec3 pos(gl.vertProperties[v * gl.numProp + 0],
                   gl.vertProperties[v * gl.numProp + 1],
                   gl.vertProperties[v * gl.numProp + 2]);
    if (std::abs(la::length(pos) - sphereRadius) > 0.1) continue;
    const vec3 n(gl.vertProperties[v * gl.numProp + normalSlotOffset + 0],
                 gl.vertProperties[v * gl.numProp + normalSlotOffset + 1],
                 gl.vertProperties[v * gl.numProp + normalSlotOffset + 2]);
    if (la::dot(n, expected(pos)) > 0.9)
      ++good;
    else
      ++bad;
  }
  return {good, bad};
}

}  // namespace

/**
 * This tests that turning a mesh into a manifold and returning it to a mesh
 * produces a consistent result.
 */
TEST(Manifold, GetMeshGL) {
  Manifold manifold = Manifold::Sphere(0.01);
  auto mesh_out = manifold.GetMeshGL();
  Manifold manifold2(mesh_out);
  auto mesh_out2 = manifold2.GetMeshGL();
  Identical(mesh_out, mesh_out2);
}

TEST(Manifold, MeshDeterminism) {
  Manifold cube1 = Manifold::Cube(vec3(2.0, 2.0, 2.0), true);
  Manifold cube2 = Manifold::Cube(vec3(2.0, 2.0, 2.0), true)
                       .Translate(vec3(-1.1091, 0.88509, 1.3099));

  Manifold result = cube1 - cube2;
  MeshGL out = result.GetMeshGL();

  uint32_t triVerts[]{0,  2,  7,  0,  10, 1,  0,  6,  10, 0, 1,  2,  1, 3,  2,
                      1,  5,  3,  1,  11, 5,  0,  7,  6,  6, 7,  8,  6, 8,  13,
                      10, 12, 11, 1,  10, 11, 11, 13, 5,  6, 12, 10, 6, 13, 12,
                      13, 9,  5,  13, 8,  9,  11, 12, 13, 4, 2,  3,  4, 3,  5,
                      4,  7,  2,  4,  5,  8,  4,  8,  7,  9, 8,  5};

  float vertProperties[]{-1,      -1,       -1,     -1,      -1,       1,
                         -1,      -0.11491, 0.3099, -1,      -0.11491, 1,
                         -0.1091, -0.11491, 0.3099, -0.1091, -0.11491, 1,
                         -1,      1,        -1,     -1,      1,        0.3099,
                         -0.1091, 1,        0.3099, -0.1091, 1,        1,
                         1,       -1,       -1,     1,       -1,       1,
                         1,       1,        -1,     1,       1,        1};

  bool flag = true;
  for (size_t i = 0; i != out.triVerts.size(); i++) {
    if (out.triVerts[i] != triVerts[i]) {
      flag = false;
      break;
    }
  }

  for (size_t i = 0; flag && i != out.vertProperties.size(); i++) {
    if (out.vertProperties[i] != vertProperties[i]) {
      flag = false;
      break;
    }
  }

  EXPECT_TRUE(flag);
}

TEST(Manifold, Empty) {
  MeshGL emptyMesh;
  Manifold empty(emptyMesh);

  EXPECT_TRUE(empty.IsEmpty());
  EXPECT_EQ(empty.Status(), Manifold::Error::NoError);
}

TEST(Manifold, ValidInput) {
  MeshGL tetGL = TetGL();
  Manifold tet(tetGL);
  EXPECT_FALSE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::NoError);
}

TEST(Manifold, ValidInputOneRunIndex) {
  MeshGL emptyMesh;
  emptyMesh.runIndex = {0};
  Manifold empty(emptyMesh);
  EXPECT_TRUE(empty.IsEmpty());
  EXPECT_EQ(empty.Status(), Manifold::Error::NoError);
}

TEST(Manifold, InvalidInput1) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold tet(in);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, InvalidInput2) {
  MeshGL in = TetGL();
  std::swap(in.triVerts[2 * 3 + 1], in.triVerts[2 * 3 + 2]);
  Manifold tet(in);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::NotManifold);
}

TEST(Manifold, InvalidInput3) {
  MeshGL in = TetGL();
  for (uint32_t& triVert : in.triVerts) {
    if (triVert == 2) triVert = -2;
  }
  Manifold tet(in);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::VertexOutOfBounds);
}

TEST(Manifold, InvalidInput4) {
  MeshGL in = TetGL();
  for (uint32_t& triVert : in.triVerts) {
    if (triVert == 2) triVert = 4;
  }
  Manifold tet(in);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::NotManifold);
}

TEST(Manifold, InvalidInput5) {
  MeshGL tetGL = TetGL();
  tetGL.mergeFromVert[tetGL.mergeFromVert.size() - 1] = 7;
  Manifold tet(tetGL);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::MergeIndexOutOfBounds);
}

TEST(Manifold, InvalidInput6) {
  MeshGL tetGL = TetGL();
  tetGL.triVerts[tetGL.triVerts.size() - 1] = 7;
  Manifold tet(tetGL);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::VertexOutOfBounds);
}

TEST(Manifold, InvalidInput7) {
  MeshGL cube = CubeUV();
  cube.runIndex = {0, 1, static_cast<uint32_t>(cube.triVerts.size())};
  Manifold tet(cube);
  EXPECT_TRUE(tet.IsEmpty());
  EXPECT_EQ(tet.Status(), Manifold::Error::RunIndexWrongLength);
}

TEST(Manifold, ErrorPropagationDecompose) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  auto parts = errored.Decompose();
  ASSERT_EQ(parts.size(), 1);
  EXPECT_EQ(parts[0].Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationHull) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Hull().Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationHullMulti) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  Manifold good = Manifold::Cube();
  EXPECT_EQ(Manifold::Hull({good, errored}).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationSetProperties) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.SetProperties(1, nullptr).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationCalculateCurvature) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.CalculateCurvature(0, 1).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationCalculateNormals) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.CalculateNormals(0).Status(),
            Manifold::Error::NonFiniteVertex);
}

// CalculateNormals(idx) followed by GetMeshGL() (no idx) used to drop
// the transform-on-export step, returning input-frame data.
TEST(Manifold, NormalsCavity) {
  // The #1712 repro: inner-sphere normals from a Boolean diff should point
  // toward the origin (outward from the surrounding solid).
  MeshGL mesh = (Manifold::Sphere(10.0, 32) - Manifold::Sphere(3.0, 32))
                    .CalculateNormals(0)
                    .GetMeshGL();
  ASSERT_GE(mesh.numProp, 6);
  auto [good, bad] = CountSphereNormalAlignment(
      mesh, 3.0, 3, [](vec3 pos) { return la::normalize(-pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, NormalsRotateBeforeCalc) {
  // Rotation before CalculateNormals: SetNormals computes from already-
  // rotated faceNormal_ and stores world-frame at slot 0.
  MeshGL mesh = Manifold::Sphere(10.0, 32)
                    .Rotate(45, 0, 0)
                    .CalculateNormals(0)
                    .GetMeshGL();
  auto [_, bad] = CountSphereNormalAlignment(
      mesh, 10.0, 3, [](vec3 pos) { return la::normalize(pos); });
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, NormalsRotateAfterCalc) {
  // Rotation *after* CalculateNormals: Impl::Transform eager-transforms the
  // stored slot 0..2 so it tracks the new orientation.
  MeshGL mesh = Manifold::Sphere(10.0, 32)
                    .CalculateNormals(0)
                    .Rotate(45, 0, 0)
                    .GetMeshGL();
  auto [_, bad] = CountSphereNormalAlignment(
      mesh, 10.0, 3, [](vec3 pos) { return la::normalize(pos); });
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, NormalsAutoSubstitute) {
  // No-arg invocation: defaults to slot 0 and sets the per-run hasNormals
  // bit on every output run.
  MeshGL mesh = Manifold::Sphere(10.0, 32).CalculateNormals().GetMeshGL();
  ASSERT_GE(mesh.numProp, 6);
  ASSERT_GT(mesh.NumRun(), 0u);
  EXPECT_TRUE(mesh.HasNormals(0));
}

TEST(Manifold, NormalsRoundTrip) {
  // getMesh -> ofMesh -> getMesh preserves the per-run flag, so the second
  // getMesh still emits world-frame normals.
  Manifold round = (Manifold::Sphere(10.0, 32) - Manifold::Sphere(3.0, 32))
                       .CalculateNormals();
  MeshGL out1 = round.GetMeshGL();
  EXPECT_TRUE(out1.HasNormals(0));
  MeshGL out2 = Manifold(out1).GetMeshGL();
  EXPECT_TRUE(out2.HasNormals(0));
  auto [good, bad] = CountSphereNormalAlignment(
      out2, 3.0, 3, [](vec3 pos) { return la::normalize(-pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, NormalsRefinePreserved) {
  // Refine keeps the recording: linearly-interpolated normals at the new
  // verts are less precise than recomputed ones but still meaningful.
  MeshGL mesh =
      Manifold::Sphere(10.0, 32).CalculateNormals().Refine(2).GetMeshGL();
  ASSERT_GT(mesh.NumRun(), 0u);
  EXPECT_TRUE(mesh.HasNormals(0));
}

TEST(Manifold, NormalsSmoothByNormalsNoArg) {
  // Smoke test that the no-arg SmoothByNormals reads from the recorded
  // slot 0 and produces a valid manifold.
  Manifold smoothed =
      Manifold::Sphere(10.0, 32).CalculateNormals().SmoothByNormals();
  EXPECT_EQ(smoothed.Status(), Manifold::Error::NoError);
}

TEST(Manifold, NormalsNonStandardSlotNotRecorded) {
  // CalculateNormals(non-zero) does NOT set the per-run recording, since a
  // non-standard slot can't be safely auto-substituted on GetMeshGL(-1).
  MeshGL mesh = Manifold::Sphere(10.0, 32).CalculateNormals(3).GetMeshGL();
  ASSERT_GT(mesh.NumRun(), 0u);
  EXPECT_FALSE(mesh.HasNormals(0));
}

TEST(Manifold, NormalsSharedPropVertMixedFlagsUndefined) {
  // hasNormals is per-run, but a single propVert holds one slot 0..2
  // value. If a hand-built MeshGL shares a propVert between a
  // hasNormals=true run and a hasNormals=false run, a Transform rotates
  // the slot on behalf of the hasNormals=true camp; the other camp's
  // interpretation (e.g. color) is collateral damage. Standard
  // CalculateNormals / Boolean / Compose outputs never share propVerts
  // this way - this test pins the documented behaviour for the only
  // input shape that can produce it.
  MeshGL gl;
  gl.numProp = 6;
  const float pts[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int v = 0; v < 4; ++v) {
    for (int j : {0, 1, 2}) gl.vertProperties.push_back(pts[v][j]);
    // Slot 0..2: +Z. Run 0 treats as normal, Run 1 treats as color.
    gl.vertProperties.push_back(0);
    gl.vertProperties.push_back(0);
    gl.vertProperties.push_back(1);
  }
  gl.triVerts = {0, 1, 3, 0, 2, 1, 0, 3, 2, 1, 2, 3};
  gl.runOriginalID = {Manifold::ReserveIDs(1), Manifold::ReserveIDs(1)};
  gl.runIndex = {0, 6, 12};
  gl.runFlags = {0x02, 0x00};  // run 0 hasNormals, run 1 plain color

  Manifold m(gl);
  ASSERT_EQ(m.Status(), Manifold::Error::NoError);

  MeshGL out = m.Rotate(90, 0, 0).GetMeshGL();

  int run0Bad = 0;  // hasNormals camp: slot rotates as expected.
  int run1Bad = 0;  // no-normals camp: collateral, slot is now rotated.
  for (size_t run = 0; run < out.NumRun(); ++run) {
    const bool isNormalsRun = out.HasNormals(run);
    for (uint32_t i = out.runIndex[run]; i < out.runIndex[run + 1]; ++i) {
      const uint32_t v = out.triVerts[i];
      const vec3 val(out.vertProperties[v * out.numProp + 3],
                     out.vertProperties[v * out.numProp + 4],
                     out.vertProperties[v * out.numProp + 5]);
      if (isNormalsRun) {
        const vec3 expected(0, -1, 0);
        if (la::length(val - expected) > 0.01) ++run0Bad;
      } else {
        const vec3 expected(0, 0, 1);
        if (la::length(val - expected) > 0.01) ++run1Bad;
      }
    }
  }
  EXPECT_EQ(run0Bad, 0);
  EXPECT_GT(run1Bad, 0);
}

TEST(Manifold, GetNormalLegacyContract) {
  // Pre-#1718, slot N normals were stored in per-mesh frame and runTransform
  // had to be applied on read - there was no runFlags bit 1 to mark
  // world-frame storage. GetNormal must still honour that contract when
  // reading a MeshGL without the bit set, or SmoothByNormals on legacy data
  // produces wrong tangents.
  //
  // Build twins: take a CalculateNormals'd rotated cube, emit both as a
  // modern MeshGL (world-frame + bit 1 set) and as a legacy MeshGL (the same
  // normals inverse-rotated into per-mesh frame, bit 1 cleared). Both should
  // produce the same SmoothByNormals.Refine output if GetNormal recomposes
  // correctly.
  const Manifold rotated =
      Manifold::Cube({1, 1, 1}, true).Rotate(30, 45, 0).CalculateNormals(0);
  MeshGL gl_modern = rotated.GetMeshGL();
  ASSERT_GT(gl_modern.NumRun(), 0u);
  ASSERT_TRUE(gl_modern.HasNormals(0));

  MeshGL gl_legacy = gl_modern;
  std::vector<bool> visited(gl_legacy.NumVert(), false);
  for (size_t run = 0; run < gl_legacy.NumRun(); ++run) {
    // Invert the per-run normal transform that the legacy GetNormal will
    // apply on read, so the stored values look "per-mesh frame" again.
    const mat3 fwd =
        la::inverse(la::transpose(mat3(gl_legacy.GetRunTransform(run)))) *
        (gl_legacy.Backside(run) ? -1.0 : 1.0);
    const mat3 inv = la::inverse(fwd);
    // `data() + i` rather than `&v[i]`: the latter is UB when `i == size()`
    // (one-past-the-end of the final run) and traps under libc++ fast
    // hardening - see #1735.
    for (uint32_t* itr = gl_legacy.triVerts.data() + gl_legacy.runIndex[run];
         itr < gl_legacy.triVerts.data() + gl_legacy.runIndex[run + 1]; ++itr) {
      const uint32_t v = *itr;
      if (visited[v]) continue;
      visited[v] = true;
      vec3 n(gl_legacy.vertProperties[v * gl_legacy.numProp + 3],
             gl_legacy.vertProperties[v * gl_legacy.numProp + 4],
             gl_legacy.vertProperties[v * gl_legacy.numProp + 5]);
      n = inv * n;
      gl_legacy.vertProperties[v * gl_legacy.numProp + 3] = n.x;
      gl_legacy.vertProperties[v * gl_legacy.numProp + 4] = n.y;
      gl_legacy.vertProperties[v * gl_legacy.numProp + 5] = n.z;
    }
    gl_legacy.runFlags[run] &= ~uint8_t(2);
  }
  ASSERT_FALSE(gl_legacy.HasNormals(0));

  Manifold m_modern(gl_modern);
  Manifold m_legacy(gl_legacy);

  Manifold sm_modern = m_modern.SmoothByNormals(0).Refine(4);
  Manifold sm_legacy = m_legacy.SmoothByNormals(0).Refine(4);
  EXPECT_NEAR(sm_modern.Volume(), sm_legacy.Volume(), 1e-4);
  EXPECT_NEAR(sm_modern.SurfaceArea(), sm_legacy.SurfaceArea(), 1e-4);
}

TEST(Manifold, TransformMixedNormalsPerMeshID) {
  // Mixed Boolean output (no-normals + with-normals): the result's
  // AllHaveNormals() is false (AND across meshIDs), but the with-normals
  // meshIDs still hold world-frame normals at slot 0..2 that must rotate
  // with a subsequent Transform. Impl::Transform must per-meshID iterate,
  // not skip on the whole-impl flag.
  MeshGL gl = MixedNormalsCubePlusSphere().Rotate(90, 0, 0).GetMeshGL();
  auto [good, bad] = CountSphereNormalAlignment(
      gl, 5.0, 3, [](vec3 pos) { return la::normalize(pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, ComposeMixedNormalsPerMeshID) {
  // Compose's per-node eager-transform check is per-Relation, not whole-node.
  // Pass a mixed Manifold with a pending Rotate into a 3-input BatchBoolean,
  // forcing the disjoint Compose path where node->transform_ is non-identity
  // and node->pImpl_->AllHaveNormals() is false (mixed). Sphere's normals
  // within the mixed input must still rotate.
  const Manifold mixed_rot = MixedNormalsCubePlusSphere().Rotate(90, 0, 0);
  const Manifold t1 = Manifold::Tetrahedron().Translate({-50, 0, 0});
  const Manifold t2 = Manifold::Tetrahedron().Translate({-100, 0, 0});
  MeshGL gl =
      Manifold::BatchBoolean({mixed_rot, t1, t2}, OpType::Add).GetMeshGL();
  auto [good, bad] = CountSphereNormalAlignment(
      gl, 5.0, 3, [](vec3 pos) { return la::normalize(pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, BooleanSubtractMixedQPerMeshIDNegation) {
  // CreateProperties' cavity sign-flip must be per-source-meshID. Subtract a
  // mixed Q (some meshIDs with hasNormals, some without). Cavity verts from
  // the hasNormals meshIDs need their slot 0..2 flipped to point outward
  // from the result solid (into the cavity = toward sphere center).
  const Manifold A = Manifold::Cube({100, 100, 100}, true);
  MeshGL gl = (A - MixedNormalsCubePlusSphere()).GetMeshGL();
  auto [good, bad] = CountSphereNormalAlignment(
      gl, 5.0, 3, [](vec3 pos) { return la::normalize(-pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, CalculateNormalsNonZeroIdxSurvivesTransform) {
  // Non-zero normalIdx is the legacy deferred-transform path: SetNormals
  // stores slot N in per-mesh frame, and GetMeshGL(N)'s legacy export
  // applies the per-run runTransform to recover world-frame. This must
  // survive transforms applied between CalculateNormals and GetMeshGL.
  const int idx = 3;
  MeshGL gl = Manifold::Sphere(5.0, 32)
                  .Rotate(30)
                  .CalculateNormals(idx)
                  .Rotate(60)
                  .GetMeshGL(idx);
  ASSERT_GE(gl.numProp, 3 + idx + 3);
  auto [good, bad] = CountSphereNormalAlignment(
      gl, 5.0, 3 + idx, [](vec3 pos) { return la::normalize(pos); });
  EXPECT_GT(good, 0);
  EXPECT_EQ(bad, 0);
}

TEST(Manifold, ErrorPropagationSmoothByNormals) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.SmoothByNormals(0).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationSmoothOut) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.SmoothOut().Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationRefine) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Refine(2).Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.RefineToLength(0.1).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.RefineToTolerance(0.1).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationSetTolerance) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.SetTolerance(0.1).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationAsOriginal) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.AsOriginal().Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationWarp) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Warp([](vec3& v) {}).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.WarpBatch([](VecView<vec3>) {}).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationMinkowski) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  Manifold good = Manifold::Cube();
  EXPECT_EQ(errored.MinkowskiSum(good).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(good.MinkowskiSum(errored).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.MinkowskiDifference(good).Status(),
            Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(good.MinkowskiDifference(errored).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationSplitByPlane) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  auto parts = errored.SplitByPlane(vec3(0, 0, 1), 0);
  EXPECT_EQ(parts.first.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(parts.second.Status(), Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationMirror) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Mirror(vec3(1, 0, 0)).Status(),
            Manifold::Error::NonFiniteVertex);
  // Degenerate normal (zero vector) on errored input should still propagate.
  EXPECT_EQ(errored.Mirror(vec3(0, 0, 0)).Status(),
            Manifold::Error::NonFiniteVertex);
}

TEST(Manifold, ErrorPropagationSimplify) {
  MeshGL in = TetGL();
  in.vertProperties[2 * 3 + 1] = NAN;
  Manifold errored(in);
  ASSERT_EQ(errored.Status(), Manifold::Error::NonFiniteVertex);
  EXPECT_EQ(errored.Simplify().Status(), Manifold::Error::NonFiniteVertex);
}

#ifndef MANIFOLD_NO_IOSTREAM
TEST(Manifold, ObjRoundTrip) {
  Manifold m = Manifold::Cube();
  std::stringstream ss;
  m.WriteOBJ(ss);
  ss.seekg(0);
  Manifold m2 = Manifold::ReadOBJ(ss);
  EXPECT_EQ(m2.Status(), Manifold::Error::NoError);
  EXPECT_EQ(m2.Volume(), 1);
}
#endif

TEST(Manifold, OppositeFace) {
  MeshGL gl;
  gl.vertProperties = {
      0, 0, 0,  //
      1, 0, 0,  //
      0, 1, 0,  //
      1, 1, 0,  //
      0, 0, 1,  //
      1, 0, 1,  //
      0, 1, 1,  //
      1, 1, 1,  //
      2, 0, 0,  //
      2, 1, 0,  //
      2, 0, 1,  //
      2, 1, 1,  //
  };
  gl.triVerts = {
      0, 1,  4,   //
      0, 2,  3,   //
      0, 3,  1,   //
      0, 4,  2,   //
      1, 3,  5,   //
      1, 3,  9,   //
      1, 5,  3,   //
      1, 5,  4,   //
      1, 8,  5,   //
      1, 9,  8,   //
      2, 4,  6,   //
      2, 6,  7,   //
      2, 7,  3,   //
      3, 5,  7,   //
      3, 7,  5,   //
      3, 7,  11,  //
      3, 11, 9,   //
      4, 5,  6,   //
      5, 7,  6,   //
      5, 8,  10,  //
      5, 10, 7,   //
      7, 10, 11,  //
      8, 9,  10,  //
      9, 11, 10,  //
  };
  Manifold man(gl);
  EXPECT_EQ(man.Status(), Manifold::Error::NoError);
  EXPECT_EQ(man.NumVert(), 12);
  EXPECT_FLOAT_EQ(man.Volume(), 2);
}

/**
 * ExpectMeshes performs a decomposition, so this test ensures that compose and
 * decompose are inverse operations.
 */
TEST(Manifold, Decompose) {
  std::vector<Manifold> manifoldList;
  manifoldList.emplace_back(Manifold::Tetrahedron());
  manifoldList.emplace_back(Manifold::Cube().Translate({2, 0, 0}).AsOriginal());
  manifoldList.emplace_back(
      Manifold::Sphere(1, 4).Translate({4, 0, 0}).AsOriginal());
  Manifold manifolds = Manifold::BatchBoolean(manifoldList, OpType::Add);

  ExpectMeshes(manifolds, {{8, 12}, {6, 8}, {4, 4}});

  std::vector<MeshGL> input;

  for (const Manifold& manifold : manifoldList) {
    EXPECT_GE(manifold.OriginalID(), 0);
    input.emplace_back(manifold.GetMeshGL());
  }

  RelatedGL(manifolds, input);
}

TEST(Manifold, DecomposeProps) {
  std::vector<MeshGL> input;
  std::vector<Manifold> manifoldList;
  auto tet = WithPositionColors(Manifold::Tetrahedron());
  manifoldList.emplace_back(tet);
  input.emplace_back(tet.GetMeshGL());
  auto cube =
      WithPositionColors(Manifold::Cube().Translate({2, 0, 0}).AsOriginal());
  manifoldList.emplace_back(cube);
  input.emplace_back(cube.GetMeshGL());
  auto sphere = WithPositionColors(
      Manifold::Sphere(1, 4).Translate({4, 0, 0}).AsOriginal());
  manifoldList.emplace_back(sphere);
  input.emplace_back(sphere.GetMeshGL());
  Manifold manifolds = Manifold::BatchBoolean(manifoldList, OpType::Add);

  ExpectMeshes(manifolds, {{8, 12, 3}, {6, 8, 3}, {4, 4, 3}});

  RelatedGL(manifolds, input);

  for (const Manifold& manifold : manifolds.Decompose()) {
    RelatedGL(manifold, input);
  }
}

/**
 * These tests check the various manifold constructors.
 */
TEST(Manifold, Sphere) {
  int n = 25;
  Manifold sphere = Manifold::Sphere(1.0, 4 * n);
  EXPECT_EQ(sphere.NumTri(), n * n * 8);
}

TEST(Manifold, Cylinder) {
  int n = 10000;
  Manifold cylinder = Manifold::Cylinder(2, 2, 2, n);
  EXPECT_EQ(cylinder.NumTri(), 4 * n - 4);
}

TEST(Manifold, CylinderZeroRadiusLow) {
  // cylinder(h, 0, r) should produce a cone with apex at z=0 and base at z=h.
  // Use enough segments for accurate volume.
  const int n = 256;
  const double h = 5.0, r = 3.0;
  Manifold coneApexBottom = Manifold::Cylinder(h, 0.0, r, n);
  Manifold coneApexTop = Manifold::Cylinder(h, r, 0.0, n);

  EXPECT_EQ(coneApexBottom.Status(), Manifold::Error::NoError);
  EXPECT_FALSE(coneApexBottom.IsEmpty());

  // Both cones must have equal total volume.
  const double totalVol = coneApexTop.Volume();
  EXPECT_NEAR(coneApexBottom.Volume(), totalVol, 1e-6);

  // Differentiate orientation by intersecting with the bottom half (z in
  // [0, h/2]). A sub-cone scaled by 1/2 in all dimensions has 1/8 the volume,
  // so apex-at-bottom gives V/8 and apex-at-top gives 7*V/8.
  Manifold slicer = Manifold::Cube(vec3(2 * r + 1, 2 * r + 1, h / 2))
                        .Translate(vec3(-(r + 0.5), -(r + 0.5), 0.0));
  EXPECT_NEAR((coneApexBottom ^ slicer).Volume(), totalVol / 8.0, 0.01);
  EXPECT_NEAR((coneApexTop ^ slicer).Volume(), 7.0 * totalVol / 8.0, 0.01);
}

TEST(Manifold, Extrude) {
  Polygons polys = SquareHole();
  Manifold donut = Manifold::Extrude(polys, 1.0, 3);
  EXPECT_EQ(donut.Genus(), 1);
  EXPECT_FLOAT_EQ(donut.Volume(), 12.0);
  EXPECT_FLOAT_EQ(donut.SurfaceArea(), 48.0);
}

TEST(Manifold, ExtrudeCone) {
  Polygons polys = SquareHole();
  Manifold donut = Manifold::Extrude(polys, 1.0, 0, 0, vec2(0.0));
  EXPECT_EQ(donut.Genus(), 0);
  EXPECT_FLOAT_EQ(donut.Volume(), 4.0);
}

Polygons RotatePolygons(Polygons polys, const int index) {
  Polygons rotatedPolys;
  for (auto& polygon : polys) {
    auto rotatedPolygon = polygon;
    std::rotate(rotatedPolygon.begin(), rotatedPolygon.begin() + index,
                rotatedPolygon.end());
    rotatedPolys.push_back(rotatedPolygon);
  }
  return rotatedPolys;
}

TEST(Manifold, Revolve) {
  Polygons polys = SquareHole();
  Manifold vug;
  for (size_t i = 0; i < polys[0].size(); i++) {
    Polygons rotatedPolys = RotatePolygons(polys, i);
    vug = Manifold::Revolve(rotatedPolys, 48);
    EXPECT_EQ(vug.Genus(), -1);
    EXPECT_NEAR(vug.Volume(), 14.0 * kPi, 0.2);
    EXPECT_NEAR(vug.SurfaceArea(), 30.0 * kPi, 0.2);
  }
}

TEST(Manifold, Revolve2) {
  Polygons polys = SquareHole(2.0);
  Manifold donutHole = Manifold::Revolve(polys, 48);
  EXPECT_EQ(donutHole.Genus(), 0);
  EXPECT_NEAR(donutHole.Volume(), 48.0 * kPi, 1.0);
  EXPECT_NEAR(donutHole.SurfaceArea(), 96.0 * kPi, 1.0);
}

TEST(Manifold, Revolve3) {
  CrossSection circle = CrossSection::Circle(1, 32);
  Manifold sphere = Manifold::Revolve(circle.ToPolygons(), 32);
  EXPECT_NEAR(sphere.Volume(), 4.0 / 3.0 * kPi, 0.1);
  EXPECT_NEAR(sphere.SurfaceArea(), 4 * kPi, 0.15);
}

TEST(Manifold, RevolveClip) {
  Polygons polys = {{{-5, -10}, {5, 0}, {-5, 10}}};
  Polygons clipped = {{{0, -5}, {5, 0}, {0, 5}}};
  Manifold first = Manifold::Revolve(polys, 48);
  Manifold second = Manifold::Revolve(clipped, 48);
  EXPECT_EQ(first.Genus(), second.Genus());
  EXPECT_EQ(first.Volume(), second.Volume());
  EXPECT_EQ(first.SurfaceArea(), second.SurfaceArea());
}

TEST(Manifold, PartialRevolveOnYAxis) {
  Polygons polys = SquareHole(2.0);
  Polygons offsetPolys = SquareHole(10.0);

  Manifold revolute;
  for (size_t i = 0; i < polys[0].size(); i++) {
    Polygons rotatedPolys = RotatePolygons(polys, i);
    revolute = Manifold::Revolve(rotatedPolys, 48, 180);
    EXPECT_EQ(revolute.Genus(), 1);
    EXPECT_NEAR(revolute.Volume(), 24.0 * kPi, 1.0);
    EXPECT_NEAR(revolute.SurfaceArea(),
                48.0 * kPi + 4.0 * 4.0 * 2.0 - 2.0 * 2.0 * 2.0, 1.0);
  }
}

TEST(Manifold, PartialRevolveOffset) {
  Polygons polys = SquareHole(10.0);

  Manifold revolute;
  for (size_t i = 0; i < polys[0].size(); i++) {
    Polygons rotatedPolys = RotatePolygons(polys, i);
    revolute = Manifold::Revolve(rotatedPolys, 48, 180);
    EXPECT_EQ(revolute.Genus(), 1);
    EXPECT_NEAR(revolute.SurfaceArea(), 777.0, 1.0);
    EXPECT_NEAR(revolute.Volume(), 376.0, 1.0);
  }
}

TEST(Manifold, Warp) {
  CrossSection square = CrossSection::Square({1, 1});
  Manifold shape =
      Manifold::Extrude(square.ToPolygons(), 2, 10).Warp([](vec3& v) {
        v.x += v.z * v.z;
      });

  Manifold simplified = Manifold::BatchBoolean({shape}, OpType::Add);

  EXPECT_NEAR(shape.Volume(), simplified.Volume(), 0.0001);
  EXPECT_NEAR(shape.SurfaceArea(), simplified.SurfaceArea(), 0.0001);
  EXPECT_NEAR(shape.Volume(), 2, 0.0001);
}

TEST(Manifold, Warp2) {
  CrossSection circle = CrossSection::Circle(5, 20).Translate(vec2(10.0, 10.0));

  Manifold shape =
      Manifold::Extrude(circle.ToPolygons(), 2, 10).Warp([](vec3& v) {
        int nSegments = 10;
        double angleStep = 2.0 / 3.0 * kPi / nSegments;
        int zIndex = nSegments - 1 - std::round(v.z);
        double angle = zIndex * angleStep;
        v.z = v.y;
        v.y = v.x * sin(angle);
        v.x = v.x * cos(angle);
      });

  Manifold simplified = Manifold::BatchBoolean({shape}, OpType::Add);

  EXPECT_NEAR(shape.Volume(), simplified.Volume(), 0.0001);
  EXPECT_NEAR(shape.SurfaceArea(), simplified.SurfaceArea(), 0.0001);
  EXPECT_NEAR(shape.Volume(), 321, 1);
}

TEST(Manifold, WarpBatch) {
  Manifold cube = Manifold::Cube({2, 3, 4});
  const int id = cube.OriginalID();

  Manifold shape1 = cube.Warp([](vec3& v) { v.x += v.z * v.z; });
  Manifold shape2 = cube.WarpBatch([](VecView<vec3> vecs) {
    for (vec3& v : vecs) {
      v.x += v.z * v.z;
    }
  });

  EXPECT_GE(id, 0);
  EXPECT_EQ(shape1.OriginalID(), -1);
  EXPECT_EQ(shape2.OriginalID(), -1);
  std::vector<uint32_t> runOriginalID1 = shape1.GetMeshGL().runOriginalID;
  EXPECT_EQ(runOriginalID1.size(), 1);
  EXPECT_EQ(runOriginalID1[0], id);
  std::vector<uint32_t> runOriginalID2 = shape2.GetMeshGL().runOriginalID;
  EXPECT_EQ(runOriginalID2.size(), 1);
  EXPECT_EQ(runOriginalID2[0], id);
  EXPECT_EQ(shape1.Volume(), shape2.Volume());
  EXPECT_EQ(shape1.SurfaceArea(), shape2.SurfaceArea());
}

TEST(Manifold, Project) {
  MeshGL input;
  input.numProp = 3;
  input.vertProperties = {0,    0,       0,     //
                          -2,   -0.7,    -0.1,  //
                          -2,   -0.7,    0,     //
                          -1.9, -0.7,    -0.1,  //
                          -1.9, -0.6901, -0.1,  //
                          -1.9, -0.7,    0,     //
                          -1.9, -0.6901, 0,     //
                          -2,   -1,      3,     //
                          -1.9, -1,      3,     //
                          -2,   -1,      4,     //
                          -1.9, -1,      4,     //
                          -1.9, -0.6901, 3,     //
                          -1.9, -0.6901, 4,     //
                          -1.7, -0.6901, 3,     //
                          -1.7, -0.6901, 3.2,   //
                          -2,   0,       -0.1,  //
                          -2,   0,       0,     //
                          -2,   0,       3,     //
                          -2,   0,       4,     //
                          -1.7, 0,       3,     //
                          -1.7, 0,       3.2,   //
                          -1,   -0.6901, -0.1,  //
                          -1,   -0.6901, 0,     //
                          -1,   -0.6901, 3.2,   //
                          -1,   -0.6901, 4,     //
                          -1,   0,       -0.1,  //
                          -1,   0,       0,     //
                          -1,   0,       3.2,   //
                          -1,   0,       4};
  input.triVerts = {1,  3,  2,   //
                    1,  4,  3,   //
                    2,  3,  5,   //
                    5,  6,  2,   //
                    3,  4,  6,   //
                    5,  3,  6,   //
                    6,  4,  21,  //
                    26, 22, 25,  //
                    21, 25, 22,  //
                    25, 15, 26,  //
                    26, 6,  22,  //
                    21, 4,  25,  //
                    21, 22, 6,   //
                    16, 26, 15,  //
                    16, 6,  26,  //
                    4,  15, 25,  //
                    15, 1,  16,  //
                    16, 2,  6,   //
                    4,  1,  15,  //
                    1,  2,  16,  //
                    12, 14, 23,  //
                    12, 13, 14,  //
                    12, 11, 13,  //
                    18, 9,  12,  //
                    11, 7,  17,  //
                    7,  9,  18,  //
                    17, 7,  18,  //
                    13, 11, 19,  //
                    17, 18, 20,  //
                    19, 11, 17,  //
                    19, 17, 20,  //
                    14, 13, 20,  //
                    18, 12, 24,  //
                    20, 13, 19,  //
                    20, 18, 27,  //
                    12, 10, 11,  //
                    24, 12, 23,  //
                    9,  10, 12,  //
                    9,  8,  10,  //
                    8,  11, 10,  //
                    8,  7,  11,  //
                    8,  9,  7,   //
                    14, 20, 27,  //
                    24, 28, 18,  //
                    27, 18, 28,  //
                    23, 14, 27,  //
                    24, 23, 28,  //
                    28, 23, 27};
  Manifold in(input);
  CrossSection projected = in.Project();
  EXPECT_NEAR(projected.Area(), 0.72, 0.01);
}

/**
 * Testing more advanced Manifold operations.
 */

TEST(Manifold, Transform) {
  Manifold cube = Manifold::Cube({1, 2, 3});
  Manifold cube2 = cube;
  cube = cube.Rotate(30, 40, 50).Scale({6, 5, 4}).Translate({1, 2, 3});

  mat3 rX({1.0, 0.0, 0.0},            //
          {0.0, cosd(30), sind(30)},  //
          {0.0, -sind(30), cosd(30)});
  mat3 rY({cosd(40), 0.0, -sind(40)},  //
          {0.0, 1.0, 0.0},             //
          {sind(40), 0.0, cosd(40)});
  mat3 rZ({cosd(50), sind(50), 0.0},   //
          {-sind(50), cosd(50), 0.0},  //
          {0.0, 0.0, 1.0});
  mat3 s;
  s[0][0] = 6;
  s[1][1] = 5;
  s[2][2] = 4;
  mat3x4 transform = mat3x4(s * rZ * rY * rX, vec3(0.0));
  transform[3] = vec3(1, 2, 3);
  cube2 = cube2.Transform(transform);

  Identical(cube.GetMeshGL(), cube2.GetMeshGL());
}

TEST(Manifold, Slice) {
  Manifold cube = Manifold::Cube();
  CrossSection bottom = cube.Slice();
  CrossSection top = cube.Slice(1);
  EXPECT_EQ(bottom.Area(), 1);
  EXPECT_EQ(top.Area(), 0);
}

TEST(Manifold, SliceEmptyObject) {
  Manifold empty;
  EXPECT_TRUE(empty.IsEmpty());
  CrossSection bottom = empty.Slice();
}

TEST(Manifold, Simplify) {
  Polygons polyCircle =
      CrossSection::Circle(1, 20).Translate({10, 0}).ToPolygons();
  Manifold torus = Manifold::Revolve(polyCircle, 100);
  Manifold simplified = torus.Simplify(0.4);
  EXPECT_NEAR(torus.Volume(), simplified.Volume(), 25);
  EXPECT_NEAR(torus.SurfaceArea(), simplified.SurfaceArea(), 10);

  if (options.exportModels) WriteTestOBJ("torus.obj", simplified);
}

// ---- Step 9 arrangement tests (docs/OverlapRemoval.md) ----
// Characterization tests on synthetic literal chords: they prove the
// step-9 bookkeeping/geometry, not reachability from a real mesh
// (open until steps 10-13 exist). All ids >= baseId == 0 against an
// empty Impl, so positions come from newVertPositions directly.

TEST(OverlapRemoval, Step9GroupChordsByFace) {
  using overlap_removal::PiercedNewEdge;
  const std::vector<PiercedNewEdge> chords = {
      {0, 1, 2, 5},  // chord 0 on faces 2 and 5
      {2, 3, 5, 7},  // chord 1 on faces 5 and 7
  };
  const std::vector<std::vector<int>> byFace =
      overlap_removal::GroupChordsByFace(chords, 8);
  ASSERT_EQ(byFace.size(), 8u);
  EXPECT_TRUE(byFace[0].empty());
  ASSERT_EQ(byFace[2].size(), 1u);
  EXPECT_EQ(byFace[2][0], 0);
  ASSERT_EQ(byFace[5].size(), 2u);
  EXPECT_EQ(byFace[5][0], 0);
  EXPECT_EQ(byFace[5][1], 1);
  ASSERT_EQ(byFace[7].size(), 1u);
  EXPECT_EQ(byFace[7][0], 1);
}

namespace {
// Two same-face chords: d (chord 0) spans ids 0->1 along x; c (chord 1)
// has endpoint id 2 at (0.5, h, 0) - distance h from d's interior at
// t = 0.5 - and far endpoint id 3. h is the knob the cases turn.
overlap_removal::NewEdgeWithExtras MakeChord(int v0, int v1, int triA,
                                             int triB) {
  return {{v0, v1, triA, triB}, {}, {}};
}
std::vector<manifold::vec3> EndpointOnChordPositions(double h) {
  return {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, h, 0.0}, {0.5, 1.0, 0.0}};
}
}  // namespace

TEST(OverlapRemoval, Step9EndpointOnChordContact) {
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;  // empty: baseId == 0
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  const std::vector<manifold::vec3> pos =
      EndpointOnChordPositions(/*h=*/1e-9);  // within tolerance + eps
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const std::vector<overlap_removal::OnChordContact> contacts =
      overlap_removal::FindOnChordEndpointContacts(impl, chords, pos, byFace,
                                                   tolerance, eps);
  ASSERT_EQ(contacts.size(), 1u);
  EXPECT_EQ(contacts[0].chord, 0);
  EXPECT_EQ(contacts[0].vertId, 2);
  EXPECT_NEAR(contacts[0].t, 0.5, 1e-12);
}

TEST(OverlapRemoval, Step9EndpointOnChordRespectsDistance) {
  // h far beyond tolerance + eps: no contact.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  const std::vector<manifold::vec3> pos = EndpointOnChordPositions(1e-6);
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  EXPECT_TRUE(overlap_removal::FindOnChordEndpointContacts(
                  impl, chords, pos, byFace, tolerance, eps)
                  .empty());
}

TEST(OverlapRemoval, Step9EndpointOnChordExcludesEndpointZone) {
  // The contact point projects to t ~= 1e-12 on d - inside the
  // endpoint-proximity zone (snap/len = 2e-9) - so it must NOT be
  // inserted as an interior vert (round-3 finding 1b).
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  std::vector<manifold::vec3> pos = EndpointOnChordPositions(1e-10);
  pos[2] = {1e-12, 1e-10, 0.0};  // near d's v0 endpoint
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  EXPECT_TRUE(overlap_removal::FindOnChordEndpointContacts(
                  impl, chords, pos, byFace, tolerance, eps)
                  .empty());
}

TEST(OverlapRemoval, Step9EndpointOnChordRequiresSharedFace) {
  // Same geometry as the hit case, but the chords share no face.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 4, 11)};
  const std::vector<manifold::vec3> pos = EndpointOnChordPositions(1e-9);
  const std::vector<std::vector<int>> byFace = {{0}, {}, {}, {}, {1}};
  EXPECT_TRUE(overlap_removal::FindOnChordEndpointContacts(
                  impl, chords, pos, byFace, tolerance, eps)
                  .empty());
}

namespace {
const manifold::vec3 kStep9FaceNormal{0.0, 0.0, 1.0};
}  // namespace

TEST(OverlapRemoval, Step9SimpleCrossing) {
  // Two chords forming an X in face 0's plane (normal z): one raw
  // crossing at (0.5, 0, 0), t = 0.5 on both; resolution finds no
  // existing vert within tolerance + eps, so a fresh id is allocated
  // and threaded onto both chords.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(
          impl, chords, pos, byFace,
          manifold::VecView<const manifold::vec3>(&kStep9FaceNormal, 1), eps);
  ASSERT_EQ(raw.size(), 1u);
  EXPECT_EQ(raw[0].chordA, 0);
  EXPECT_EQ(raw[0].chordB, 1);
  EXPECT_NEAR(raw[0].tA, 0.5, 1e-12);
  EXPECT_NEAR(raw[0].tB, 0.5, 1e-12);
  EXPECT_NEAR(raw[0].pos.x, 0.5, 1e-12);
  EXPECT_NEAR(raw[0].pos.y, 0.0, 1e-12);
  // Perpendicular crossing: the conditioned radius sits at its floor.
  EXPECT_NEAR(raw[0].snapR, eps, 1e-12 * eps);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 4);  // fresh: baseId 0 + 4 existing
  ASSERT_EQ(threaded.newVertPositions.size(), 5u);
  ASSERT_EQ(threaded.newVertSnapR.size(), 5u);  // parallel to the pool
  EXPECT_NEAR(threaded.newVertSnapR[4], eps, 1e-12 * eps);
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);
  EXPECT_NEAR(threaded.chords[0].extraTs[0], 0.5, 1e-12);
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 4);
  EXPECT_NEAR(threaded.chords[1].extraTs[0], 0.5, 1e-12);
}

TEST(OverlapRemoval, Step9CrossingResolvesToEndpointBand) {
  // The split-identity regression (design round 3): a crossing in the
  // (eps, tolerance+eps] band of chord 0's endpoint id 1. The kernel
  // accepts it (> eps from endpoints), but resolution must give the
  // ENDPOINT id to the crossing on BOTH chords - never a fresh id on
  // one and the endpoint on the other. On chord 0 itself the
  // recomputed t for id 1 is 1.0, outside the endpoint-zone guard, so
  // chord 0 threads nothing.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  const double x0 = 1.0 - 1.5e-9;  // in (eps, tolerance+eps] of id 1
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {x0, -0.5, 0.0}, {x0, 0.5, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(
          impl, chords, pos, byFace,
          manifold::VecView<const manifold::vec3>(&kStep9FaceNormal, 1), eps);
  ASSERT_EQ(raw.size(), 1u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 1);           // snapped, not allocated
  EXPECT_EQ(threaded.newVertPositions.size(), 4u);  // no fresh vert
  EXPECT_TRUE(threaded.chords[0].extraVerts.empty());
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 1);
  EXPECT_NEAR(threaded.chords[1].extraTs[0], 0.5, 1e-6);
}

TEST(OverlapRemoval, Step9CrossingSeesPassZeroContacts) {
  // A third chord's endpoint (id 4) rests on both crossing chords
  // within tolerance + eps; pass 0 records it onto both. The c-x-d
  // crossing lands within tolerance + eps of id 4, so resolution must
  // pick id 4 (consulting the pass-0 accumulator - design round 4),
  // allocate nothing, and the unified id-dedup must leave exactly one
  // entry per chord.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // d: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // c: diagonal through (0.5, 0, 0)
      MakeChord(4, 5, 0, 13)};  // e: endpoint 4 hovers 1e-9 above (0.5,0,0)
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},  {1.0, 0.0, 0.0},
                                     {0.0, -0.5, 0.0}, {1.0, 0.5, 0.0},
                                     {0.5, 1e-9, 0.0}, {0.9, 1.0, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1, 2}};
  const std::vector<overlap_removal::OnChordContact> contacts =
      overlap_removal::FindOnChordEndpointContacts(impl, chords, pos, byFace,
                                                   tolerance, eps);
  ASSERT_FALSE(contacts.empty());  // id 4 rests on d (and on c)
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(
          impl, chords, pos, byFace,
          manifold::VecView<const manifold::vec3>(&kStep9FaceNormal, 1), eps);
  ASSERT_EQ(raw.size(), 1u);  // only c x d properly cross
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(impl, std::move(chords),
                                                 std::move(pos), {}, raw,
                                                 contacts, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 4);           // pass-0 vert, not fresh
  EXPECT_EQ(threaded.newVertPositions.size(), 6u);  // no allocation
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 1u);  // d: one entry
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);  // c: one entry
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 4);
}

TEST(OverlapRemoval, Step9ThreeConcurrentChordsShareOneVert) {
  // Three chords concurrent at (0.5, 0, 0): the three pairwise
  // crossings must merge into ONE cluster, resolve to ONE fresh id,
  // and thread exactly once onto each of the three chords.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // (0.5,-0.5,0)-(0.5,0.5,0)
      MakeChord(4, 5, 0, 13)};  // (0,-0.5,0)-(1,0.5,0): y = x - 0.5
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},  {1.0, 0.0, 0.0},
                                     {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0},
                                     {0.0, -0.5, 0.0}, {1.0, 0.5, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1, 2}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 3u);  // all three pairs cross
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normals, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);
  EXPECT_EQ(merged[0].chords.size(), 3u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadClusters(impl, std::move(chords),
                                                std::move(pos), {}, merged, {},
                                                tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 6);           // one fresh vert
  EXPECT_EQ(threaded.newVertPositions.size(), 7u);  // exactly one allocation
  for (int ci : {0, 1, 2}) {
    ASSERT_EQ(threaded.chords[ci].extraVerts.size(), 1u) << "chord " << ci;
    EXPECT_EQ(threaded.chords[ci].extraVerts[0], 6) << "chord " << ci;
    EXPECT_NEAR(threaded.chords[ci].extraTs[0], 0.5, 1e-9) << "chord " << ci;
  }
}

TEST(OverlapRemoval, Step9FaceGateMergesDisjointPairs) {
  // The face-gate regression (design round 3): crossings (c1,c2) and
  // (c3,c4) share NO chord - a chord-gated merge would leave two
  // distinct verts 5e-12 apart at a genuine 4-chord concurrence. The
  // face gate + 10x-eps radius must unite them; resolution then snaps
  // to the strictly-nearest existing endpoint, c2's top (id 3, 0.5e-12
  // from the centroid vs 1e-12 for c4's bottom).
  const double eps = 1e-12;
  const double tolerance = 1e-12;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // c1: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // c2: (0.5,-0.1,0)-(0.5,2e-12,0)
      MakeChord(4, 5, 0, 13),   // c3: (0,5e-12,0)-(1,5e-12,0)
      MakeChord(6, 7, 0, 15)};  // c4: (0.5,3.5e-12,0)-(0.5,0.6,0)
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},     {1.0, 0.0, 0.0},
                                     {0.5, -0.1, 0.0},    {0.5, 2e-12, 0.0},
                                     {0.0, 5e-12, 0.0},   {1.0, 5e-12, 0.0},
                                     {0.5, 3.5e-12, 0.0}, {0.5, 0.6, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1, 2, 3}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::OnChordContact> contacts =
      overlap_removal::FindOnChordEndpointContacts(impl, chords, pos, byFace,
                                                   tolerance, eps);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 2u);  // only (c1,c2) and (c3,c4) properly cross
  EXPECT_TRUE(raw[0].chordA != raw[1].chordA && raw[0].chordB != raw[1].chordB);
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normals, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);  // face gate united the disjoint pairs
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadClusters(impl, std::move(chords),
                                                std::move(pos), {}, merged,
                                                contacts, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 3);           // snapped: tie -> id 3
  EXPECT_EQ(threaded.newVertPositions.size(), 8u);  // no allocation
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 3);
  EXPECT_TRUE(threaded.chords[1].extraVerts.empty());  // id 3 is c2's endpoint
  ASSERT_EQ(threaded.chords[2].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[2].extraVerts[0], 3);
  EXPECT_TRUE(threaded.chords[3].extraVerts.empty());  // outside c4's span
}

TEST(OverlapRemoval, Step9MergeReprojectsOntoHostFacePlane) {
  // Multi-face cluster: crossing r0 hosts on face 7 (plane z = 0) and
  // r1 on face 2 (plane x = 0.5); they unite via shared incident
  // faces. host = min(7, 2) = 2, but members[0] is r0 - the plane
  // POINT for re-projection must come from a member whose hosting
  // face IS the host (r1), or the centroid is projected onto a plane
  // through face 7's vertex with face 2's normal, landing on neither
  // face (here x would collapse to 0 instead of staying 0.5).
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 7, 9),   // c1 in z=0:    (0,0,0)-(1,0,0)
      MakeChord(2, 3, 7, 11),  // c2 in z=0:    (0.5,-0.5,0)-(0.5,0.5,0)
      MakeChord(4, 5, 2, 7),   // c3 in x=0.5:  (0.5,-0.5,5e-9)-(0.5,0.5,5e-9)
      MakeChord(6, 7, 2, 9)};  // c4 in x=0.5:  varies z through y=0
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0},   {1.0, 0.0, 0.0},  {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0},
      {0.5, -0.5, 5e-9}, {0.5, 0.5, 5e-9}, {0.5, 0.0, -0.5}, {0.5, 0.0, 0.5}};
  std::vector<std::vector<int>> byFace(8);
  byFace[2] = {2, 3};
  byFace[7] = {0, 1};
  std::vector<manifold::vec3> normals(8, manifold::vec3(0.0, 0.0, 1.0));
  normals[2] = manifold::vec3(1.0, 0.0, 0.0);
  const manifold::VecView<const manifold::vec3> normalsView(normals.data(),
                                                            normals.size());
  std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normalsView, eps);
  ASSERT_EQ(raw.size(), 2u);  // (c3,c4) at (0.5,0,5e-9); (c1,c2) at (0.5,0,0)
  // FindChordChordCrossings emits raw ascending by face, which makes
  // the first cluster member's face coincide with the host face. The
  // merge must not DEPEND on that caller ordering - reverse it so the
  // first member hosts on face 7 while the host face is 2.
  std::reverse(raw.begin(), raw.end());
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normalsView, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);              // shared incident faces unite
  EXPECT_NEAR(merged[0].pos.x, 0.5, 1e-12);  // on the host plane, not x=0
  EXPECT_NEAR(merged[0].pos.y, 0.0, 1e-12);
}

TEST(OverlapRemoval, Step9ThreadingPreservesStepEightExtras) {
  // A pre-existing step-8 on-tri vert (id 4, t = 1e-9 - inside the
  // endpoint-proximity zone the step-9 guard excludes) must survive a
  // threading rebuild triggered by an unrelated crossing: the
  // endpoint-zone guard applies to step-9-ADDED records only, not to
  // verts step 8 already admitted under its weaker guard.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  chords[0].extraVerts = {4};
  chords[0].extraTs = {1e-9};
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},
                                     {1.0, 0.0, 0.0},
                                     {0.5, -0.5, 0.0},
                                     {0.5, 0.5, 0.0},
                                     {1e-9, 0.0, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 1u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 2u);  // id 4 preserved
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);
  EXPECT_EQ(threaded.chords[0].extraVerts[1], 5);  // the fresh crossing id
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 5);
}

TEST(OverlapRemoval, Step9CollinearOverlapSnapsEndpointsNoCrossing) {
  // Collinear overlapping chords: the kernel yields NO single crossing
  // (the overlap interval defers to the step-12 multiplicity merge),
  // but pass 0 must snap the overlap endpoints onto the other chord -
  // counted, not silently dropped.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // d: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11)};  // e: (0.3,0,0)-(1.3,0,0), overlap [0.3,1]
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.3, 0.0, 0.0}, {1.3, 0.0, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  EXPECT_TRUE(overlap_removal::FindChordChordCrossings(impl, chords, pos,
                                                       byFace, normals, eps)
                  .empty());
  const std::vector<overlap_removal::OnChordContact> contacts =
      overlap_removal::FindOnChordEndpointContacts(impl, chords, pos, byFace,
                                                   tolerance, eps);
  ASSERT_EQ(contacts.size(), 2u);
  EXPECT_EQ(contacts[0].chord, 0);  // e's start rests on d at t = 0.3
  EXPECT_EQ(contacts[0].vertId, 2);
  EXPECT_NEAR(contacts[0].t, 0.3, 1e-12);
  EXPECT_EQ(contacts[1].chord, 1);  // d's end rests on e at t = 0.7
  EXPECT_EQ(contacts[1].vertId, 1);
  EXPECT_NEAR(contacts[1].t, 0.7, 1e-12);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(impl, std::move(chords),
                                                 std::move(pos), {}, {},
                                                 contacts, tolerance, eps);
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 2);
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 1);
}

TEST(OverlapRemoval, Step9ResolutionSnapsToRealMeshVert) {
  // Mixed id space: a real Impl (tetrahedron, baseId == 4) supplies a
  // chord endpoint with id < baseId. A crossing in the
  // (eps, tolerance+eps] band of that REAL vert must resolve to its id
  // through GetPos3's impl.vertPos_ branch - the path every
  // empty-Impl test misses.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl(Manifold::Tetrahedron().GetMeshGL64());
  ASSERT_EQ(impl.NumVert(), 4u);  // baseId == 4
  const manifold::vec3 v0 = impl.vertPos_[0];
  // Chord d: real id 0 -> new id 4 along +x from v0; chord c crosses d
  // at v0 + (1.5e-9, 0, 0), inside the band of real vert 0.
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 4, 0, 9), MakeChord(5, 6, 0, 11)};
  std::vector<manifold::vec3> newPos = {
      v0 + manifold::vec3(1.0, 0.0, 0.0),      // id 4
      v0 + manifold::vec3(1.5e-9, -0.5, 0.0),  // id 5
      v0 + manifold::vec3(1.5e-9, 0.5, 0.0)};  // id 6
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, newPos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 1u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(impl, std::move(chords),
                                                 std::move(newPos), {}, raw, {},
                                                 tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 0);              // the REAL mesh vert
  EXPECT_EQ(threaded.newVertPositions.size(), 3u);     // no allocation
  EXPECT_TRUE(threaded.chords[0].extraVerts.empty());  // id 0 is d's endpoint
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 0);
  EXPECT_NEAR(threaded.chords[1].extraTs[0], 0.5, 1e-6);
}

TEST(OverlapRemoval, Step9IdDedupAcrossExtrasAndCrossing) {
  // Explicit id-dedup: a step-8 extraVert already threads id 4 on the
  // chord, and a crossing RESOLVES to the same id 4 - the unified
  // dedup must leave exactly one entry, not two.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  chords[0].extraVerts = {4};
  chords[0].extraTs = {0.5};
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},
                                     {1.0, 0.0, 0.0},
                                     {0.5, -0.5, 0.0},
                                     {0.5, 0.5, 0.0},
                                     {0.5, 1e-9, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 1u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 4);  // snapped to the step-8 vert
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 1u);  // deduped, not two
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 4);
}

TEST(OverlapRemoval, Step9ThreadingRecomputeReordersAfterSnap) {
  // The round-4 t-recompute rule, pinned against its exact regression:
  // two crossings on one chord whose SNAPPED positions invert their
  // pre-resolution t-order. tolerance is deliberately large (0.03) so
  // each crossing snaps to its slanted chord's near endpoint, moving
  // ALONG the host chord: pre t = (0.45, 0.46), post t = (0.475,
  // 0.435). Threading must sort by the recomputed t's - stale ts give
  // [id2, id4]; recomputed give [id4, id2].
  const double eps = 1e-9;
  const double tolerance = 0.03;
  Manifold::Impl impl;
  // c1 crosses d at 0.45 from ABOVE with its near endpoint at x=0.475;
  // c2 crosses d at 0.46 from BELOW with its near endpoint at x=0.435.
  // Approaching from opposite sides keeps c1 strictly left of c2 over
  // their shared y-band, so c1 and c2 themselves never cross.
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // d: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // c1: snap target (0.475, +0.001)
      MakeChord(4, 5, 0, 13)};  // c2: snap target (0.435, -0.001)
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},      {1.0, 0.0, 0.0},
                                     {0.475, 0.001, 0.0},  {0.39, -0.0024, 0.0},
                                     {0.435, -0.001, 0.0}, {0.52, 0.0024, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1, 2}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 2u);  // (d,c1) and (d,c2) only
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 2u);
  EXPECT_EQ(threaded.crossings[0].id, 2);  // snapped to c1's near endpoint
  EXPECT_EQ(threaded.crossings[1].id, 4);  // snapped to c2's near endpoint
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 2u);
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);  // 0.435 first
  EXPECT_EQ(threaded.chords[0].extraVerts[1], 2);  // 0.475 second
  EXPECT_LT(threaded.chords[0].extraTs[0], threaded.chords[0].extraTs[1]);
}

TEST(OverlapRemoval, Step9MergeRadiusCollapsesNearbyDistinctCrossings) {
  // Documented collapse: two genuinely distinct crossings 9 * eps
  // apart in one face merge by the 10x-eps radius alone (shrinking the
  // radius below 9x would break this pin). The off-chord cost (the
  // merged vert sits 4.5e-9 off each vertical chord) is the accepted
  // error-budget trade.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // d: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // c1: vertical at x = 0.5
      MakeChord(4, 5, 0, 13)};  // c2: vertical at x = 0.5 + 9e-9
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},         {0.5, -0.5, 0.0},
      {0.5, 0.5, 0.0}, {0.5 + 9e-9, -0.5, 0.0}, {0.5 + 9e-9, 0.5, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1, 2}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 2u);  // (d,c1) and (d,c2); c1 || c2
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normals, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);  // collapsed by the radius, not structure
  EXPECT_NEAR(merged[0].pos.x, 0.5 + 4.5e-9, 1e-12);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadClusters(impl, std::move(chords),
                                                std::move(pos), {}, merged, {},
                                                tolerance, eps);
  EXPECT_EQ(threaded.newVertPositions.size(), 7u);  // one fresh vert
  for (int ci : {0, 1, 2}) {
    ASSERT_EQ(threaded.chords[ci].extraVerts.size(), 1u) << "chord " << ci;
    EXPECT_EQ(threaded.chords[ci].extraVerts[0], 6) << "chord " << ci;
  }
}

TEST(OverlapRemoval, Step9ShallowCrossingCarriesConditionedRadius) {
  // Review finding (cross-pair conditioning, step-9 side): a shallow
  // crossing's position error scales as eps / sin(angle); the radius
  // must ride the crossing record into the cluster and the allocated
  // vert so step 9.5 can widen its new-onto-original snap. Chord 1
  // rises 0.01 per unit x -> sin ~ 0.01 -> condR ~ 100 eps (inside
  // the 128-eps cap).
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -0.005, 0.0}, {1.0, 0.005, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 1u);
  EXPECT_GT(raw[0].snapR, 90.0 * eps);
  EXPECT_LT(raw[0].snapR, 110.0 * eps);
  const double rawSnapR = raw[0].snapR;
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadCrossings(
          impl, std::move(chords), std::move(pos), {}, raw, {}, tolerance, eps);
  ASSERT_EQ(threaded.newVertPositions.size(), 5u);  // fresh vert id 4
  ASSERT_EQ(threaded.newVertSnapR.size(), 5u);
  EXPECT_EQ(threaded.newVertSnapR[4], rawSnapR);  // carried, not reset
  for (int j : {0, 1, 2, 3}) {  // pre-existing pool entries default eps
    EXPECT_EQ(threaded.newVertSnapR[j], eps) << "pool vert " << j;
  }
}

TEST(OverlapRemoval, Step9ClusterCarriesMaxConditionedRadius) {
  // The nearby-crossing merge takes the WIDEST member radius: the
  // merged point stands for every member crossing, including the
  // worst-conditioned one.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9), MakeChord(2, 3, 0, 11)};
  const std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0}};
  const std::vector<std::vector<int>> byFace = {{0, 1}};
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  // Two hand-built raw records for the same pair, 5 eps apart (within
  // the 10 eps merge radius, sharing face 0), radii 5 eps and 30 eps.
  const std::vector<overlap_removal::ChordChordCrossing> raw = {
      {{0.5, 0.0, 0.0}, 0, 1, 0.5, 0.5, 0, 5.0 * eps},
      {{0.5 + 5e-9, 0.0, 0.0}, 0, 1, 0.5, 0.5, 0, 30.0 * eps}};
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normals, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);
  EXPECT_EQ(merged[0].snapR, 30.0 * eps);
}

TEST(OverlapRemoval, Step9EndToEndComposition) {
  // The five step-9 functions composed the way the production driver
  // will call them: grouping -> pass 0 -> crossings -> merge ->
  // resolve/thread, on the 3-concurrent fixture plus a resting
  // endpoint, all in one flow.
  const double eps = 1e-9;
  const double tolerance = 1e-9;
  Manifold::Impl impl;
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      MakeChord(0, 1, 0, 9),    // d: (0,0,0)-(1,0,0)
      MakeChord(2, 3, 0, 11),   // c: vertical through (0.5, 0)
      MakeChord(4, 5, 0, 13)};  // e: endpoint 4 rests on d at t = 0.25
  std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},   {1.0, 0.0, 0.0},
                                     {0.5, -0.5, 0.0},  {0.5, 0.5, 0.0},
                                     {0.25, 1e-9, 0.0}, {0.25, 1.0, 0.0}};
  std::vector<overlap_removal::PiercedNewEdge> newEdges;
  for (const overlap_removal::NewEdgeWithExtras& nwe : chords) {
    newEdges.push_back(nwe.edge);
  }
  const std::vector<std::vector<int>> byFace =
      overlap_removal::GroupChordsByFace(newEdges, 1);
  ASSERT_EQ(byFace.size(), 1u);
  ASSERT_EQ(byFace[0].size(), 3u);
  const manifold::VecView<const manifold::vec3> normals(&kStep9FaceNormal, 1);
  const std::vector<overlap_removal::OnChordContact> contacts =
      overlap_removal::FindOnChordEndpointContacts(impl, chords, pos, byFace,
                                                   tolerance, eps);
  ASSERT_EQ(contacts.size(), 1u);  // e's endpoint on d
  EXPECT_EQ(contacts[0].vertId, 4);
  const std::vector<overlap_removal::ChordChordCrossing> raw =
      overlap_removal::FindChordChordCrossings(impl, chords, pos, byFace,
                                               normals, eps);
  ASSERT_EQ(raw.size(), 1u);  // only d x c properly cross
  const std::vector<overlap_removal::ChordCrossing> merged =
      overlap_removal::MergeAndPropagateCrossings(
          impl, chords, pos, raw, byFace, normals, tolerance, eps);
  ASSERT_EQ(merged.size(), 1u);
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadClusters(impl, std::move(chords),
                                                std::move(pos), {}, merged,
                                                contacts, tolerance, eps);
  ASSERT_EQ(threaded.crossings.size(), 1u);
  EXPECT_EQ(threaded.crossings[0].id, 6);               // fresh: far from id 4
  EXPECT_EQ(threaded.newVertPositions.size(), 7u);      // one allocation
  ASSERT_EQ(threaded.chords[0].extraVerts.size(), 2u);  // d: contact + cross
  EXPECT_EQ(threaded.chords[0].extraVerts[0], 4);       // t = 0.25
  EXPECT_EQ(threaded.chords[0].extraVerts[1], 6);       // t = 0.5
  ASSERT_EQ(threaded.chords[1].extraVerts.size(), 1u);  // c: the crossing
  EXPECT_EQ(threaded.chords[1].extraVerts[0], 6);
  EXPECT_TRUE(threaded.chords[2].extraVerts.empty());  // e: id 4 is its own
}

// ---- Steps 10-11 partition tests (docs/OverlapRemoval.md) ----
// Fixtures run on a real tetrahedron Impl. Manifold::Impl(MeshGL64)
// runs SortGeometry, which permutes vert and face ids, so the fixture
// recovers the z = 0 face A(0,0,0), B(1,0,0), C(0,1,0) and its corner
// ids FROM POSITIONS. Chords and on-edge verts are synthetic
// step-9-style literals in that plane; ids 4+ index newVertPositions.

namespace {
int Step10EdgeIndex(const std::vector<overlap_removal::Edge>& edges, int a,
                    int b) {
  const int v0 = std::min(a, b);
  const int v1 = std::max(a, b);
  for (size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].v0 == v0 && edges[i].v1 == v1) return static_cast<int>(i);
  }
  return -1;
}
bool Step10CycleIsSimple(const std::vector<int>& poly) {
  std::set<int> s(poly.begin(), poly.end());
  return s.size() == poly.size();
}
struct Step10Fixture {
  Manifold::Impl impl;
  std::vector<overlap_removal::Edge> edges;
  std::vector<int> he2e;       // halfedge id -> edge index
  int face = -1;               // the z = 0 face, post-SortGeometry
  int A = -1, B = -1, C = -1;  // ids of (0,0,0), (1,0,0), (0,1,0)
};
Step10Fixture MakeStep10Fixture() {
  Step10Fixture f;
  MeshGL64 m;
  m.numProp = 3;
  m.vertProperties = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  m.triVerts = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
  f.impl = Manifold::Impl(m);
  f.edges = overlap_removal::EnumerateEdges(f.impl);
  f.he2e = overlap_removal::BuildHalfedgeToEdgeIndex(f.impl, f.edges);
  auto idAt = [&](double x, double y, double z) {
    for (size_t i = 0; i < f.impl.NumVert(); ++i) {
      const manifold::vec3 d = f.impl.vertPos_[i] - manifold::vec3(x, y, z);
      if (la::dot(d, d) < 1e-24) return static_cast<int>(i);
    }
    return -1;
  };
  f.A = idAt(0, 0, 0);
  f.B = idAt(1, 0, 0);
  f.C = idAt(0, 1, 0);
  const std::set<int> want = {f.A, f.B, f.C};
  for (size_t t = 0; t < f.impl.NumTri(); ++t) {
    std::set<int> got;
    for (int k = 0; k < 3; ++k) {
      got.insert(f.impl.halfedge_.Start(3 * static_cast<int>(t) + k));
    }
    if (got == want) {
      f.face = static_cast<int>(t);
      break;
    }
  }
  return f;
}
// Insert vertId (at position p) into the on-edge list of edge (a, b),
// keeping the list sorted by t along the edge's v0 -> v1.
void Step10AddOnEdge(const Step10Fixture& f,
                     std::vector<overlap_removal::EdgeVertList>& lists, int a,
                     int b, int vertId, manifold::vec3 p) {
  const int ei = Step10EdgeIndex(f.edges, a, b);
  ASSERT_GE(ei, 0);
  const manifold::vec3 v0 = f.impl.vertPos_[f.edges[ei].v0];
  const manifold::vec3 v1 = f.impl.vertPos_[f.edges[ei].v1];
  const double t = la::dot(p - v0, v1 - v0) / la::dot(v1 - v0, v1 - v0);
  overlap_removal::EdgeVertList& l = lists[ei];
  size_t pos = 0;
  while (pos < l.ts.size() && l.ts[pos] < t) ++pos;
  l.verts.insert(l.verts.begin() + pos, vertId);
  l.ts.insert(l.ts.begin() + pos, t);
}
}  // namespace

TEST(OverlapRemoval, Step10XCrossingPartitionsIntoFour) {
  // Two chords crossing at id 8 split the z = 0 face into four simple
  // polygons (three quads around the corners, one triangle in the
  // middle), every one containing the crossing vert exactly once.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // p1 = id 4 (0.5,0,0) and p3 = id 6 (0.25,0,0) on edge (A,B);
  // p2 = id 5 (0,0.5,0) on (A,C); p4 = id 7 (0.5,0.5,0) on (B,C).
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 6, {0.25, 0.0, 0.0});
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.C, 5, {0.0, 0.5, 0.0});
  Step10AddOnEdge(fx, onEdgeLists, fx.B, fx.C, 7, {0.5, 0.5, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {8}, {1.0 / 3.0}},    // p1 -> p2 through x
      {{6, 7, fx.face, 101}, {8}, {1.0 / 3.0}}};  // p3 -> p4 through x
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0},
                                              {0.0, 0.5, 0.0},
                                              {0.25, 0.0, 0.0},
                                              {0.5, 0.5, 0.0},
                                              {1.0 / 3.0, 1.0 / 6.0, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1}, newPos,
      1e-10, false);
  ASSERT_EQ(part.polygons.size(), 4u);
  EXPECT_EQ(part.spursDropped, 0);
  std::multiset<size_t> sizes;
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_TRUE(Step10CycleIsSimple(poly));
    EXPECT_EQ(std::count(poly.begin(), poly.end(), 8), 1);
    sizes.insert(poly.size());
  }
  EXPECT_EQ(sizes, (std::multiset<size_t>{3, 4, 4, 4}));
}

TEST(OverlapRemoval, Step10DanglingChordSpurDropped) {
  // A chord from the boundary to an interior dead end (its sibling
  // pair was dropped upstream): the walk U-turns at the dead end, the
  // spur splits out as a sub-3-vert loop and is dropped, and the one
  // remaining polygon is simple and excludes the interior vert.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}};  // p1 -> interior d, dangling
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}, {0.3, 0.3, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0}, newPos,
      1e-10, false);
  ASSERT_EQ(part.polygons.size(), 1u);
  EXPECT_GE(part.spursDropped, 1);
  EXPECT_TRUE(Step10CycleIsSimple(part.polygons[0]));
  EXPECT_EQ(std::count(part.polygons[0].begin(), part.polygons[0].end(), 5), 0);
  EXPECT_EQ(part.polygons[0].size(), 4u);  // corners + the on-edge vert
}

TEST(OverlapRemoval, Step10CoincidentChordsDedup) {
  // Two chords spanning the same vert pair (the collinear-overlap
  // class) dedup to one cut: the partition equals the single-chord
  // result - two simple polygons.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.C, 5, {0.0, 0.5, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}, {{4, 5, fx.face, 101}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1}, newPos,
      1e-10, false);
  ASSERT_EQ(part.polygons.size(), 2u);
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_TRUE(Step10CycleIsSimple(poly));
  }
}

TEST(OverlapRemoval, Step10BoundaryRidingChordsSkipped) {
  // The hull pancake's host-side configuration: the face's boundary
  // edge (A,B) is subdivided at vert 4, and three chords arrive as
  // [rider {A,4}, interior {4,C}, rider {4,B}] - the order that made
  // the doubled-directed-edge walk emit an uncut quad. Riders
  // coincide with boundary sub-edges and are skipped (counted); the
  // interior chord cuts the face into exactly two triangles
  // regardless of chord order.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{std::min(fx.A, 4), std::max(fx.A, 4), fx.face, 99}, {}, {}},
      {{std::min(4, fx.C), std::max(4, fx.C), fx.face, 99}, {}, {}},
      {{std::min(4, fx.B), std::max(4, fx.B), fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
      newPos, 1e-10, false);
  EXPECT_EQ(part.boundaryRidingSubEdgesSkipped, 2);
  ASSERT_EQ(part.polygons.size(), 2u);
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_EQ(poly.size(), 3u);
    EXPECT_TRUE(Step10CycleIsSimple(poly));
    EXPECT_EQ(std::count(poly.begin(), poly.end(), 4), 1) << "pierce vert";
    EXPECT_EQ(std::count(poly.begin(), poly.end(), fx.C), 1) << "apex";
  }
}

TEST(OverlapRemoval, Step10WalkFrameIgnoresStoredNormal) {
  // The walk frame derives from the halfedge WINDING, not the stored
  // faceNormal_ (folded-sheet faces arrive with the two disagreeing -
  // SortGeometry can deliver either). Inverting the stored normal must
  // not change the partition, and - the discriminating assert, since
  // cycle SIZES are topological - every output cycle must stay CCW
  // about the WINDING normal: a regression to the stored normal
  // mirrors the projection and flips the output orientation.
  Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  fx.impl.faceNormal_[fx.face] = -fx.impl.faceNormal_[fx.face];
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.C, 5, {0.0, 0.5, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0}, newPos,
      1e-10, false);
  ASSERT_EQ(part.polygons.size(), 2u);
  const int baseId = static_cast<int>(fx.impl.NumVert());
  auto posOf = [&](int id) {
    return id < baseId ? fx.impl.vertPos_[id] : newPos[id - baseId];
  };
  const manifold::vec3 p0 =
      fx.impl.vertPos_[fx.impl.halfedge_.Start(3 * fx.face)];
  const manifold::vec3 p1 =
      fx.impl.vertPos_[fx.impl.halfedge_.Start(3 * fx.face + 1)];
  const manifold::vec3 p2 =
      fx.impl.vertPos_[fx.impl.halfedge_.Start(3 * fx.face + 2)];
  const manifold::vec3 windingN = la::cross(p1 - p0, p2 - p0);
  std::multiset<size_t> sizes;
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_TRUE(Step10CycleIsSimple(poly));
    sizes.insert(poly.size());
    const manifold::vec3 origin = posOf(poly[0]);
    manifold::vec3 nsum(0, 0, 0);
    for (size_t i = 0; i < poly.size(); ++i) {
      nsum = nsum + la::cross(posOf(poly[i]) - origin,
                              posOf(poly[(i + 1) % poly.size()]) - origin);
    }
    EXPECT_GT(la::dot(nsum, windingN), 0.0) << "cycle wound against winding";
  }
  EXPECT_EQ(sizes, (std::multiset<size_t>{3, 4}));
}

TEST(OverlapRemoval, Step10BoundaryRiderResultIsChordOrderInvariant) {
  // The rider configuration again, with the chord list presented in
  // all rotations: the partition must not depend on chord order (the
  // hes-order sensitivity the rider skip was built to remove).
  Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  Step10AddOnEdge(fx, onEdgeLists, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> base = {
      {{std::min(fx.A, 4), std::max(fx.A, 4), fx.face, 99}, {}, {}},
      {{std::min(4, fx.C), std::max(4, fx.C), fx.face, 99}, {}, {}},
      {{std::min(4, fx.B), std::max(4, fx.B), fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}};
  std::vector<std::vector<std::vector<int>>> results;
  for (int rot = 0; rot < 3; ++rot) {
    std::vector<overlap_removal::NewEdgeWithExtras> chords;
    for (int k = 0; k < 3; ++k) chords.push_back(base[(k + rot) % 3]);
    const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
        fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
        newPos, 1e-10, false);
    std::vector<std::vector<int>> canon;
    for (const std::vector<int>& poly : part.polygons) {
      // Canonical rotation: smallest id first (orientation preserved).
      const auto mn = std::min_element(poly.begin(), poly.end());
      std::vector<int> c(mn, poly.end());
      c.insert(c.end(), poly.begin(), mn);
      canon.push_back(std::move(c));
    }
    std::sort(canon.begin(), canon.end());
    results.push_back(std::move(canon));
  }
  EXPECT_EQ(results[0], results[1]);
  EXPECT_EQ(results[0], results[2]);
}

TEST(OverlapRemoval, Step10FreeIslandDecomposesToAnnulusTriangles) {
  // The stamp class (free-island): a closed chord loop strictly interior
  // to the face, no connection to its boundary. After reclassification,
  // the clean free island is now hole-aware: the walk runs normally, the
  // positive-area island cycle (disk polygon) passes through, and the
  // outer region's negative-area island cycle is the hole. TriangulateIdx
  // decomposes the annulus into triangles.
  //
  // Expected: interiorIslandVerts == 0, polygons non-empty (disk polygon
  // + annulus triangles), every polygon is a simple triangle or the disk,
  // and the total signed area (outer face minus hole) is preserved.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // A triangle of chords between three interior verts (ids 4-6).
  // The face is the z=0 triangle A(0,0,0), B(1,0,0), C(0,1,0).
  // Island at (0.2,0.2), (0.5,0.2), (0.2,0.5).
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {
      {0.2, 0.2, 0.0}, {0.5, 0.2, 0.0}, {0.2, 0.5, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
      newPos, 1e-10, false);
  // The clean free island is now decomposed, not gated.
  EXPECT_EQ(part.interiorIslandVerts, 0);
  EXPECT_FALSE(part.polygons.empty());
  // All output polygons are simple and have >= 3 verts.
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_GE(poly.size(), 3u);
    EXPECT_TRUE(Step10CycleIsSimple(poly));
  }

  // Build a full position table (impl verts + new verts) for area checks.
  // Vert ids 0..NumVert()-1 come from impl; ids >= NumVert() are newPos.
  const int baseId = static_cast<int>(fx.impl.NumVert());
  auto posOf2D = [&](int id) -> manifold::vec2 {
    // Project onto xy-plane (the face is z=0).
    const manifold::vec3 p = id < baseId
                                 ? fx.impl.vertPos_[id]
                                 : newPos[static_cast<size_t>(id - baseId)];
    return {p.x, p.y};
  };
  auto signedArea2D = [&](const std::vector<int>& poly) {
    double a = 0.0;
    const int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
      const manifold::vec2 pi = posOf2D(poly[i]);
      const manifold::vec2 pj = posOf2D(poly[(i + 1) % n]);
      a += pi.x * pj.y - pi.y * pj.x;
    }
    return a * 0.5;
  };

  // The face is stored as triVerts {0,2,1} = {A,C,B}, which in xy is CW
  // (normal points -z). The walk emits all cycles in face-winding orientation
  // (CCW about face normal = CW in +z view), so every output polygon has
  // negative signed area in global xy projection.
  // Face outer cycle signed area using halfedge vertex order {A,C,B}:
  const double faceOuterArea = signedArea2D({fx.A, fx.C, fx.B});
  // Sanity: face is CW in global xy (normal -z).
  ASSERT_LT(faceOuterArea, 0.0);

  // No standalone positive-area cycle in the output: every polygon must
  // share the face winding (negative in global xy). A positive cycle would
  // be a hole that escaped the decomposition uncut.
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_LE(signedArea2D(poly), 0.0);
  }

  // Total signed area of all output polygons equals the face outer area
  // (disk + annulus triangles together tile the face exactly once).
  double totalArea = 0.0;
  for (const std::vector<int>& poly : part.polygons) {
    totalArea += signedArea2D(poly);
  }
  EXPECT_NEAR(totalArea, faceOuterArea, 1e-10);

  // At least one triangle (3-vert polygon) is present beyond the disk
  // polygon (which is itself 3 verts; the annulus triangles are additional).
  int triCount = 0;
  for (const std::vector<int>& poly : part.polygons) {
    if (poly.size() == 3) ++triCount;
  }
  // The disk is one triangle; the annulus triangulation adds more.
  EXPECT_GT(triCount, 1);
}

TEST(OverlapRemoval, Step10PinchedIslandGatesFail) {
  // The PINCHED variant: a chord loop attached to the boundary at
  // exactly ONE vert (corner A). The walk would split at A and the loop
  // still cancels (pinched annulus). Must gate: compAttach == 1.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> pinched = {
      {{std::min(fx.A, 4), std::max(fx.A, 4), fx.face, 99}, {}, {}},
      {{4, 5, fx.face, 99}, {}, {}},
      {{std::min(fx.A, 5), std::max(fx.A, 5), fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> pinchedPos = {{0.3, 0.2, 0.0},
                                                  {0.2, 0.3, 0.0}};
  const overlap_removal::FacePartition pinchedPart =
      overlap_removal::PartitionFace(fx.impl, fx.face, fx.edges, fx.he2e,
                                     onEdgeLists, pinched, {0, 1, 2},
                                     pinchedPos, 1e-10, false);
  EXPECT_GT(pinchedPart.interiorIslandVerts, 0);
  EXPECT_TRUE(pinchedPart.polygons.empty());
}

TEST(OverlapRemoval, Step10ProperDoubleCrossingPasses) {
  // A proper DOUBLE crossing - a CYCLE-bearing component attached at
  // two distinct boundary verts - splits into representable regions
  // and must PASS (the gate must not over-fire; tree components are
  // covered by the X-crossing fixture, this one pins cycles).
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  std::vector<overlap_removal::EdgeVertList> onEdge2(fx.edges.size());
  Step10AddOnEdge(fx, onEdge2, fx.A, fx.B, 4, {0.5, 0.0, 0.0});
  Step10AddOnEdge(fx, onEdge2, fx.A, fx.C, 5, {0.0, 0.5, 0.0});
  const std::vector<overlap_removal::NewEdgeWithExtras> through = {
      {{4, 6, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 5, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> throughPos = {
      {0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.3, 0.3, 0.0}};
  const overlap_removal::FacePartition throughPart =
      overlap_removal::PartitionFace(fx.impl, fx.face, fx.edges, fx.he2e,
                                     onEdge2, through, {0, 1, 2}, throughPos,
                                     1e-10, false);
  EXPECT_EQ(throughPart.interiorIslandVerts, 0);
  EXPECT_GE(throughPart.polygons.size(), 3u);
}

TEST(OverlapRemoval, Step10BranchyLoopPlusSpurGates) {
  // An unclean detached component: the island loop plus a dangling spur
  // from one of its verts. The spur makes one vert degree-3 (not all
  // degree-2), so the degree-2 proof fails and the component gates.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // Triangle island 4-5-6 plus spur from 4 to interior vert 7.
  // newPos: 4=(0.2,0.2), 5=(0.5,0.2), 6=(0.2,0.5), 7=(0.1,0.1) (spur).
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}},
      {{4, 7, fx.face, 99}, {}, {}}};  // spur from island vert
  const std::vector<manifold::vec3> newPos = {
      {0.2, 0.2, 0.0}, {0.5, 0.2, 0.0}, {0.2, 0.5, 0.0}, {0.1, 0.1, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2, 3},
      newPos, 1e-10, false);
  EXPECT_GT(part.interiorIslandVerts, 0);
}

TEST(OverlapRemoval, Step10TwoLoopsSharingVertGates) {
  // Two island loops sharing a vert: cross-component vert sharing
  // means neither is clean. Both gate.
  // Loop1: 4-5-6 triangle. Loop2: 6-7-8 triangle (shares vert 6).
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}, {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}, {{6, 7, fx.face, 99}, {}, {}},
      {{7, 8, fx.face, 99}, {}, {}}, {{6, 8, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.15, 0.15, 0.0},
                                              {0.30, 0.15, 0.0},
                                              {0.15, 0.30, 0.0},
                                              {0.45, 0.15, 0.0},
                                              {0.30, 0.30, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords,
      {0, 1, 2, 3, 4, 5}, newPos, 1e-10, false);
  EXPECT_GT(part.interiorIslandVerts, 0);
}

TEST(OverlapRemoval, Step10ZeroAreaIslandLoopGates) {
  // An island loop whose 2D projected area is exactly zero (a degenerate
  // collinear loop). Must gate even though degree-2 holds.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // A "triangle" with all three verts collinear along y=0.2.
  // 2D projected signed area will be ~0.
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {
      {0.1, 0.2, 0.0}, {0.3, 0.2, 0.0}, {0.2, 0.2, 0.0}};  // collinear
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
      newPos, 1e-10, false);
  EXPECT_GT(part.interiorIslandVerts, 0);
}

TEST(OverlapRemoval, Step10HazardFlaggedFaceWithIslandGates) {
  // A coplanar-hazard-flagged face with a clean free island must gate
  // (the hazard flag signals that triangulating independently could break
  // step-12 cancellation of anti-aligned coplanar partners).
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {
      {0.2, 0.2, 0.0}, {0.5, 0.2, 0.0}, {0.2, 0.5, 0.0}};
  // Same island as FreeIslandDecomposesToAnnulusTriangles, but hazard=true.
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
      newPos, 1e-10, true);  // coplanarHazard = true
  EXPECT_GT(part.interiorIslandVerts, 0);
  EXPECT_TRUE(part.polygons.empty());
}

TEST(OverlapRemoval, Step10TwoDisjointIslandsInOneFace) {
  // Two disjoint clean free islands in one face. Both should be
  // decomposed independently (separate hole assignments to the outer
  // region), producing annulus triangles for each.
  // Face: A(0,0,0) B(1,0,0) C(0,1,0). Two small triangles inside.
  // Island1: verts 4,5,6 at (0.1,0.1),(0.2,0.1),(0.1,0.2).
  // Island2: verts 7,8,9 at (0.4,0.1),(0.6,0.1),(0.4,0.3).
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}, {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}, {{7, 8, fx.face, 99}, {}, {}},
      {{8, 9, fx.face, 99}, {}, {}}, {{7, 9, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.1, 0.1, 0.0}, {0.2, 0.1, 0.0},
                                              {0.1, 0.2, 0.0}, {0.4, 0.1, 0.0},
                                              {0.6, 0.1, 0.0}, {0.4, 0.3, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords,
      {0, 1, 2, 3, 4, 5}, newPos, 1e-10, false);
  EXPECT_EQ(part.interiorIslandVerts, 0);
  EXPECT_FALSE(part.polygons.empty());
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_GE(poly.size(), 3u);
    EXPECT_TRUE(Step10CycleIsSimple(poly));
  }
}

TEST(OverlapRemoval, Step10NestedFreeIslandSeamPin) {
  // Nested free-island seam pin: an inner island loop (7,8,9) strictly
  // inside the disk region of an outer island loop (4,5,6), both free
  // (no boundary attachment).
  //
  // Face: A(0,0,0), B(1,0,0), C(0,1,0).
  // Outer island: 4(0.2,0.2), 5(0.5,0.2), 6(0.2,0.5) - "outer stamp".
  // Inner island: 7(0.25,0.25), 8(0.35,0.25), 9(0.25,0.35) - strictly
  //   inside the outer island triangle.
  //
  // Expected:
  //   - interiorIslandVerts == 0 (both islands processed cleanly).
  //   - The outer face annulus triangulates (triangles containing face
  //     boundary verts and outer island verts).
  //   - The OUTER ISLAND DISK region (4,5,6) is also triangulated (NOT
  //     passed through as a plain polygon) because it carries the inner
  //     island as a hole. The inner annulus produces triangles mixing
  //     outer island verts {4,5,6} and inner island verts {7,8,9}.
  //   - The inner island disk (7,8,9) passes through as a polygon.
  //   - No standalone negative cycle in the output.
  //   - Total signed area is preserved.
  //
  // This pin goes RED if process-ALL-regions is reverted to emit-then-skip:
  // in that case the outer disk (4,5,6) is emitted uncut and no polygon
  // bridges {4,5,6} and {7,8,9}.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // Outer island chords.
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}}, {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}, {{7, 8, fx.face, 99}, {}, {}},
      {{8, 9, fx.face, 99}, {}, {}}, {{7, 9, fx.face, 99}, {}, {}}};
  // newPos[0..2] = outer island verts 4,5,6; [3..5] = inner island verts 7,8,9.
  const std::vector<manifold::vec3> newPos = {
      {0.2, 0.2, 0.0},   {0.5, 0.2, 0.0},
      {0.2, 0.5, 0.0},  // outer (4,5,6)
      {0.25, 0.25, 0.0}, {0.35, 0.25, 0.0},
      {0.25, 0.35, 0.0}};  // inner (7,8,9)
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords,
      {0, 1, 2, 3, 4, 5}, newPos, 1e-10, false);

  EXPECT_EQ(part.interiorIslandVerts, 0);
  EXPECT_FALSE(part.polygons.empty());

  // Build position table for area computation.
  const int baseId = static_cast<int>(fx.impl.NumVert());
  auto posOf2D = [&](int id) -> manifold::vec2 {
    const manifold::vec3 p = id < baseId
                                 ? fx.impl.vertPos_[id]
                                 : newPos[static_cast<size_t>(id - baseId)];
    return {p.x, p.y};
  };
  auto signedArea2D = [&](const std::vector<int>& poly) {
    double a = 0.0;
    const int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
      const manifold::vec2 pi = posOf2D(poly[i]);
      const manifold::vec2 pj = posOf2D(poly[(i + 1) % n]);
      a += pi.x * pj.y - pi.y * pj.x;
    }
    return a * 0.5;
  };

  // All output polygons are simple and at least triangles.
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_GE(poly.size(), 3u);
    EXPECT_TRUE(Step10CycleIsSimple(poly));
  }

  // No standalone negative cycle: face winding is CW in global xy (normal
  // -z), so every output polygon must have negative global-xy signed area.
  for (const std::vector<int>& poly : part.polygons) {
    EXPECT_LE(signedArea2D(poly), 0.0);
  }

  // Total signed area preserved: face outer area (signed, includes the two
  // island disks minus the two holes = net face area).
  // The face outer cycle {A,C,B} in global xy is CW (negative).
  const double faceOuterArea = signedArea2D({fx.A, fx.C, fx.B});
  ASSERT_LT(faceOuterArea, 0.0);
  double totalArea = 0.0;
  for (const std::vector<int>& poly : part.polygons) {
    totalArea += signedArea2D(poly);
  }
  EXPECT_NEAR(totalArea, faceOuterArea, 1e-10);

  // Discriminating assertion: the outer island disk (verts 4,5,6) must NOT
  // appear as an uncut plain polygon. The disk region carries the inner
  // island as a hole, so it must be triangulated. In the correctly processed
  // output, some triangle has verts from both the outer island set {4,5,6}
  // and the inner island set {7,8,9} (the inner-annulus triangulation
  // bridges the two boundaries). In the emit-then-skip bug, the outer disk
  // emits as {4,5,6} with no cross-boundary triangles.
  const std::set<int> outerIslandVerts = {4, 5, 6};
  const std::set<int> innerIslandVerts = {7, 8, 9};
  bool foundCrossTri = false;
  for (const std::vector<int>& poly : part.polygons) {
    bool hasOuter = false, hasInner = false;
    for (int v : poly) {
      if (outerIslandVerts.count(v)) hasOuter = true;
      if (innerIslandVerts.count(v)) hasInner = true;
    }
    if (hasOuter && hasInner) {
      foundCrossTri = true;
      break;
    }
  }
  // Must find at least one polygon bridging outer and inner island boundaries.
  EXPECT_TRUE(foundCrossTri)
      << "no polygon bridges outer-island {4,5,6} and inner-island {7,8,9}; "
         "the outer disk was not triangulated (emit-then-skip regression)";
}

TEST(OverlapRemoval, Step10ZeroLengthChordSkippedAndCleanFace) {
  // A zero-length chord (step-9 snapping collapsed it) is skipped and
  // counted; with no effective cuts the face partitions into its own
  // boundary cycle. Also pins the chord-free clean-face path.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 4, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {{0.5, 0.0, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0}, newPos,
      1e-10, false);
  ASSERT_EQ(part.polygons.size(), 1u);
  EXPECT_EQ(part.zeroLengthChordsSkipped, 1);
  EXPECT_EQ(part.polygons[0].size(), 3u);  // the bare corner cycle
  EXPECT_TRUE(Step10CycleIsSimple(part.polygons[0]));

  const overlap_removal::FacePartition clean =
      overlap_removal::PartitionFace(fx.impl, fx.face, fx.edges, fx.he2e,
                                     onEdgeLists, {}, {}, {}, 1e-10, false);
  ASSERT_EQ(clean.polygons.size(), 1u);
  EXPECT_EQ(clean.polygons[0].size(), 3u);
}

TEST(OverlapRemoval, Step10SubResolutionIslandGates) {
  // A clean degree-2 island loop whose 2D projected area is nonzero but
  // below the sub-resolution fast-fail threshold (|signedArea2| <
  // triEps * holeBboxMax) must gate rather than decompose. Pins the
  // CLASS OUTCOME, not the fast-fail arm itself: the arm is redundant
  // by construction (its hole-own-bbox threshold is always <= the
  // triangulator's contour-set threshold, so anything it gates would
  // otherwise be classified non-hole and rejected by the validation
  // triad) - it is a short-circuit, not independently observable from
  // the outcome.
  //
  // The face is A(0,0,0), B(1,0,0), C(0,1,0) (z=0). The island is a
  // tiny equilateral-like right triangle strictly inside the face at
  // (0.2, 0.2, 0), with leg length 5e-12 - well below triEps (~1e-10)
  // so the sub-resolution condition fires deterministically.
  // signedArea2 ~ (5e-12)^2 = 25e-24; holeBboxMax = 5e-12;
  // triEps * holeBboxMax ~ 1e-10 * 5e-12 = 5e-22 >> 25e-24: gates.
  const Step10Fixture fx = MakeStep10Fixture();
  ASSERT_GE(fx.face, 0);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  // Island verts 4,5,6 at eps-scale near (0.2,0.2,0).
  const double leg = 5e-12;
  const std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{4, 5, fx.face, 99}, {}, {}},
      {{5, 6, fx.face, 99}, {}, {}},
      {{4, 6, fx.face, 99}, {}, {}}};
  const std::vector<manifold::vec3> newPos = {
      {0.2, 0.2, 0.0}, {0.2 + leg, 0.2, 0.0}, {0.2, 0.2 + leg, 0.0}};
  const overlap_removal::FacePartition part = overlap_removal::PartitionFace(
      fx.impl, fx.face, fx.edges, fx.he2e, onEdgeLists, chords, {0, 1, 2},
      newPos, 1e-10, false);
  // Must gate: sub-resolution island cannot be triangulated reliably.
  EXPECT_GT(part.interiorIslandVerts, 0);
  EXPECT_TRUE(part.polygons.empty());
}

// ---- Step 12 canonical merge tests (docs/OverlapRemoval.md) ----

TEST(OverlapRemoval, Step12RotationsMergeAndSum) {
  // The same cycle in two rotations (coplanar same-orientation
  // duplicates from two faces) merges to one entry with summed
  // multiplicity; the face of the FIRST contributor wins.
  const std::vector<std::pair<int, std::vector<int>>> in = {{7, {2, 3, 1}},
                                                            {9, {1, 2, 3}}};
  const std::vector<overlap_removal::MergedPolygon> out =
      overlap_removal::MergePolygons(in);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].cycle, (std::vector<int>{1, 2, 3}));
  EXPECT_EQ(out[0].mult, 2);
  EXPECT_EQ(out[0].face, 7);
}

TEST(OverlapRemoval, Step12OppositePairCancels) {
  // A cycle and its reversal (a doubled surface's two faces) sum to
  // zero and drop.
  const std::vector<std::pair<int, std::vector<int>>> in = {{7, {1, 2, 3, 4}},
                                                            {9, {4, 3, 2, 1}}};
  EXPECT_TRUE(overlap_removal::MergePolygons(in).empty());
}

TEST(OverlapRemoval, Step12SignRuleReversedCycle) {
  // A cycle whose canonical form is a rotation of its REVERSAL enters
  // with multiplicity -1, stored under the canonical (reversed)
  // rotation.
  const std::vector<std::pair<int, std::vector<int>>> in = {{7, {3, 2, 1}}};
  const std::vector<overlap_removal::MergedPolygon> out =
      overlap_removal::MergePolygons(in);
  ASSERT_EQ(out.size(), 1u);
  EXPECT_EQ(out[0].cycle, (std::vector<int>{1, 2, 3}));
  EXPECT_EQ(out[0].mult, -1);
}

TEST(OverlapRemoval, Step12DistinctSameVertSetStaySeparate) {
  // Two cycles over the same vert set but different cyclic order are
  // DIFFERENT polygons and must not merge (the reason step 12 keys on
  // the cycle, not Emmett's vert set).
  const std::vector<std::pair<int, std::vector<int>>> in = {{7, {1, 2, 3, 4}},
                                                            {9, {1, 3, 2, 4}}};
  const std::vector<overlap_removal::MergedPolygon> out =
      overlap_removal::MergePolygons(in);
  ASSERT_EQ(out.size(), 2u);
  EXPECT_NE(out[0].cycle, out[1].cycle);
}

TEST(OverlapRemoval, Step12DeterministicKeyOrder) {
  // Output ordered by canonical key regardless of input order.
  const std::vector<std::pair<int, std::vector<int>>> a = {{7, {5, 6, 7}},
                                                           {9, {1, 2, 3}}};
  const std::vector<std::pair<int, std::vector<int>>> b = {{9, {1, 2, 3}},
                                                           {7, {5, 6, 7}}};
  const std::vector<overlap_removal::MergedPolygon> outA =
      overlap_removal::MergePolygons(a);
  const std::vector<overlap_removal::MergedPolygon> outB =
      overlap_removal::MergePolygons(b);
  ASSERT_EQ(outA.size(), 2u);
  ASSERT_EQ(outB.size(), 2u);
  EXPECT_EQ(outA[0].cycle, outB[0].cycle);
  EXPECT_EQ(outA[1].cycle, outB[1].cycle);
  EXPECT_EQ(outA[0].cycle, (std::vector<int>{1, 2, 3}));
}

// ---- Step 13.2-13.3 cell-complex tests (docs/OverlapRemoval.md) ----
// Closed synthetic complexes with all ids in newVertPositions (empty
// Impl, baseId == 0). Side key: 2 * polygon + side, side 0 = front.

namespace {
bool Step13SameCell(const overlap_removal::CellComplex& cc, int sideA,
                    int sideB) {
  return cc.polySide2Cell[sideA] == cc.polySide2Cell[sideB];
}

// Appends the 12 outward-wound triangles of an axis-aligned box over
// the given 8 vert ids as merged polygons (mult +1), where vert i is
// the corner with x = hi.x iff i & 1, y = hi.y iff i & 2, z = hi.z
// iff i & 4.
void AppendBoxTris(const int (&ids)[8],
                   std::vector<overlap_removal::MergedPolygon>& polys) {
  const int tris[12][3] = {{0, 2, 3}, {0, 3, 1},   // -z
                           {4, 5, 7}, {4, 7, 6},   // +z
                           {0, 1, 5}, {0, 5, 4},   // -y
                           {2, 7, 3}, {2, 6, 7},   // +y
                           {0, 4, 6}, {0, 6, 2},   // -x
                           {1, 7, 5}, {1, 3, 7}};  // +x
  for (const auto& t : tris) {
    polys.push_back({{ids[t[0]], ids[t[1]], ids[t[2]]}, 1, 0});
  }
}

// Appends the 8 verts and 12 outward-wound triangles of the
// axis-aligned cube [lo, hi]^3.
void AppendCubePolys(double lo, double hi, std::vector<manifold::vec3>& pos,
                     std::vector<overlap_removal::MergedPolygon>& polys) {
  const int b = static_cast<int>(pos.size());
  int ids[8];
  for (int i = 0; i < 8; ++i) {
    pos.push_back({(i & 1) ? hi : lo, (i & 2) ? hi : lo, (i & 4) ? hi : lo});
    ids[i] = b + i;
  }
  AppendBoxTris(ids, polys);
}
}  // namespace

TEST(OverlapRemoval, Step13TetraSurfaceHasTwoCells) {
  // A closed tetra surface (4 outward-wound faces): every edge fan has
  // k = 2, and the cells are exactly inside and outside - all fronts
  // (outward) one cell, all backs the other.
  Manifold::Impl impl;
  const std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  const std::vector<overlap_removal::MergedPolygon> polys = {
      {{0, 2, 1}, 1, 0},   // base, outward -z
      {{0, 1, 3}, 1, 0},   // outward -y
      {{1, 2, 3}, 1, 0},   // outward +x+y+z
      {{2, 0, 3}, 1, 0}};  // outward -x
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.polySide2Cell.size(), 8u);
  EXPECT_EQ(cc.numCells, 2);
  ASSERT_EQ(cc.fans.size(), 6u);
  for (const overlap_removal::EdgeFan& fan : cc.fans) {
    EXPECT_EQ(fan.polygons.size(), 2u);
  }
  for (int p = 1; p < 4; ++p) {
    EXPECT_TRUE(Step13SameCell(cc, 0, 2 * p));      // fronts: outside
    EXPECT_TRUE(Step13SameCell(cc, 1, 2 * p + 1));  // backs: inside
  }
  EXPECT_FALSE(Step13SameCell(cc, 0, 1));
}

TEST(OverlapRemoval, Step13BipyramidWithInternalFaceHasThreeCells) {
  // Triangular bipyramid (apexes above and below the shared base
  // triangle) with the base itself included as an internal polygon
  // (normal +z, toward the upper solid): three cells - outside,
  // upper inside, lower inside - and the base edges carry k = 3 fans.
  Manifold::Impl impl;
  const std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},
                                           {1.0, 0.0, 0.0},
                                           {0.0, 1.0, 0.0},
                                           {0.3, 0.3, 1.0},
                                           {0.3, 0.3, -1.0}};
  // A=0 B=1 C=2 U=3 D=4. Outer faces wound outward; base (A,B,C) CCW
  // from above (front = +z = upper inside).
  const std::vector<overlap_removal::MergedPolygon> polys = {
      {{0, 1, 3}, 1, 0},   // 0: ABU upper side
      {{1, 2, 3}, 1, 0},   // 1: BCU upper side
      {{2, 0, 3}, 1, 0},   // 2: CAU upper side
      {{1, 0, 4}, 1, 0},   // 3: BAD lower side
      {{2, 1, 4}, 1, 0},   // 4: CBD lower side
      {{0, 2, 4}, 1, 0},   // 5: ACD lower side
      {{0, 1, 2}, 1, 0}};  // 6: base, front toward U
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.polySide2Cell.size(), 14u);
  EXPECT_EQ(cc.numCells, 3);
  // Base edges (A,B), (B,C), (A,C) have three incident polygons.
  int k3Fans = 0;
  for (const overlap_removal::EdgeFan& fan : cc.fans) {
    if (fan.polygons.size() == 3u) ++k3Fans;
  }
  EXPECT_EQ(k3Fans, 3);
  // Outside: all six outer fronts together.
  for (int p = 1; p < 6; ++p) {
    EXPECT_TRUE(Step13SameCell(cc, 0, 2 * p)) << "outer front " << p;
  }
  // Upper inside: the upper faces' backs + the base's front.
  EXPECT_TRUE(Step13SameCell(cc, 1, 3));
  EXPECT_TRUE(Step13SameCell(cc, 1, 5));
  EXPECT_TRUE(Step13SameCell(cc, 1, 2 * 6));
  // Lower inside: the lower faces' backs + the base's back.
  EXPECT_TRUE(Step13SameCell(cc, 7, 9));
  EXPECT_TRUE(Step13SameCell(cc, 7, 11));
  EXPECT_TRUE(Step13SameCell(cc, 7, 2 * 6 + 1));
  // The three cells are distinct.
  EXPECT_FALSE(Step13SameCell(cc, 0, 1));
  EXPECT_FALSE(Step13SameCell(cc, 0, 7));
  EXPECT_FALSE(Step13SameCell(cc, 1, 7));
}

TEST(OverlapRemoval, Step13CubeClassifyKeepsAllFaces) {
  // A single outward cube: two cells with windings 0 (outside) and 1
  // (inside); every face separates them and is kept. The bottom face
  // is one quad polygon passed in its REVERSED representation (cycle
  // reversed, mult -1): classification must be representation-
  // invariant, and only that polygon - whose canonical normal points
  // inside - gets the emit-orientation flip. The quad is the first
  // polygon, so the seed cast targets its first Triangulate ear.
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  std::vector<overlap_removal::MergedPolygon> polys;
  AppendCubePolys(0.0, 1.0, pos, polys);
  // Replace the two bottom tris with the reversed quad: outward is
  // {0, 2, 3, 1} (Newell -z), so the reversed form carries mult -1.
  polys.erase(polys.begin(), polys.begin() + 2);
  polys.insert(polys.begin(), {{1, 3, 2, 0}, -1, 0});
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.numCells, 2);
  const overlap_removal::CellWinding cw =
      overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/-1.0);
  ASSERT_TRUE(cw.ok);
  EXPECT_EQ(cw.seedCasts, 1);
  ASSERT_EQ(cw.winding.size(), 2u);
  ASSERT_EQ(cw.keep.size(), 11u);
  ASSERT_EQ(cw.flip.size(), 11u);
  for (int p = 0; p < 11; ++p) {
    EXPECT_TRUE(cw.keep[p]) << "polygon " << p;
    const int wFront = cw.winding[cc.polySide2Cell[2 * p]];
    const int wBack = cw.winding[cc.polySide2Cell[2 * p + 1]];
    if (p == 0) {
      // Reversed representation: front (+canonical normal) is inside.
      EXPECT_EQ(wFront, 1);
      EXPECT_EQ(wBack, 0);
      EXPECT_TRUE(cw.flip[p]);
    } else {
      EXPECT_EQ(wFront, 0) << "polygon " << p;
      EXPECT_EQ(wBack, 1) << "polygon " << p;
      EXPECT_FALSE(cw.flip[p]) << "polygon " << p;
    }
  }
}

TEST(OverlapRemoval, Step13NestedCubesInnerFacesNotKept) {
  // Two disjoint outward cubes, one inside the other: two cell-graph
  // components (no shared arrangement edges, so the between-region is
  // represented by two cells, one per component). The inner
  // component's seed cast passes through the outer wall, measuring
  // its true ambient winding 1; the inner faces then separate w 1|2 -
  // both inside - and are NOT kept, while the outer faces (0|1) are.
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  std::vector<overlap_removal::MergedPolygon> polys;
  AppendCubePolys(-2.0, 2.0, pos, polys);
  AppendCubePolys(-0.5, 0.5, pos, polys);
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.numCells, 4);
  const overlap_removal::CellWinding cw =
      overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/-1.0);
  ASSERT_TRUE(cw.ok);
  EXPECT_EQ(cw.seedCasts, 2);
  ASSERT_EQ(cw.winding.size(), 4u);
  ASSERT_EQ(cw.keep.size(), 24u);
  for (int p = 0; p < 12; ++p) {
    EXPECT_TRUE(cw.keep[p]) << "outer " << p;
    EXPECT_FALSE(cw.flip[p]) << "outer " << p;
    EXPECT_EQ(cw.winding[cc.polySide2Cell[2 * p]], 0) << "outer front " << p;
    EXPECT_EQ(cw.winding[cc.polySide2Cell[2 * p + 1]], 1) << "outer back " << p;
  }
  for (int p = 12; p < 24; ++p) {
    EXPECT_FALSE(cw.keep[p]) << "inner " << p;
    EXPECT_FALSE(cw.flip[p]) << "inner " << p;
    EXPECT_EQ(cw.winding[cc.polySide2Cell[2 * p]], 1) << "inner front " << p;
    EXPECT_EQ(cw.winding[cc.polySide2Cell[2 * p + 1]], 2) << "inner back " << p;
  }
}

TEST(OverlapRemoval, Step13SeedCastSkipsMembranes) {
  // Review finding: an open sheet (k = 1 rims unite its front and
  // back - a membrane) separates nothing and is excluded from the
  // winding BFS, so the seed cast must not count crossings of it
  // either - else the seeded winding disagrees with what the BFS
  // propagates from it. A huge membrane hangs between the cast
  // source and the cube: with the bug the cube seeds at w = 2 and
  // nothing is kept.
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  std::vector<overlap_removal::MergedPolygon> polys;
  AppendCubePolys(0.0, 1.0, pos, polys);  // polys 0-11
  pos.push_back({-10.0, -10.0, 2.0});     // 8
  pos.push_back({20.0, -10.0, 2.0});      // 9
  pos.push_back({-10.0, 20.0, 2.0});      // 10
  polys.push_back({{8, 9, 10}, 1, 0});    // poly 12: the membrane
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.numCells, 3);  // outside, inside, membrane (united)
  ASSERT_EQ(cc.polySide2Cell[2 * 12], cc.polySide2Cell[2 * 12 + 1]);  // premise
  const overlap_removal::CellWinding cw =
      overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/-1.0);
  ASSERT_TRUE(cw.ok);
  for (int p = 0; p < 12; ++p) {
    EXPECT_TRUE(cw.keep[p]) << "cube poly " << p;
  }
  EXPECT_FALSE(cw.keep[12]);  // the membrane separates nothing
}

TEST(OverlapRemoval, Step13BookTwinPairingSplitsSharedEdge) {
  // Two boxes sharing exactly one arrangement edge (verts 3 and 7 at
  // x = y = 1): the shared fan carries 4 kept faces. Inside-wedge
  // twin pairing must pair each box's own two faces (bare fan
  // adjacency would pair across an OUTSIDE wedge and weld the
  // solids), so ring extraction yields TWO output verts at each
  // shared vert - 16 rings over 14 geometric verts - and a closed
  // surface: every directed ring-id edge appears once, with its
  // antiparallel twin. One tri is passed in reversed representation
  // (cycle reversed, mult -1) to pin flip handling; all output cycles
  // must come out wound outward.
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  std::vector<overlap_removal::MergedPolygon> polys;
  AppendCubePolys(0.0, 1.0, pos, polys);  // box A: verts 0-7
  // Box B = [1,2] x [1,2] x [0,1] shares A's verts 3 = (1,1,0) and
  // 7 = (1,1,1) (its corners 0 and 4).
  pos.push_back({2.0, 1.0, 0.0});  // 8
  pos.push_back({1.0, 2.0, 0.0});  // 9
  pos.push_back({2.0, 2.0, 0.0});  // 10
  pos.push_back({2.0, 1.0, 1.0});  // 11
  pos.push_back({1.0, 2.0, 1.0});  // 12
  pos.push_back({2.0, 2.0, 1.0});  // 13
  const int bIds[8] = {3, 8, 9, 10, 7, 11, 12, 13};
  AppendBoxTris(bIds, polys);
  // Reversed representation for one of A's bottom tris.
  ASSERT_EQ(polys[0].cycle, (std::vector<int>{0, 2, 3}));
  polys[0] = {{3, 2, 0}, -1, 0};
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  // Outside, inside-A, inside-B; one component (joined via outside).
  ASSERT_EQ(cc.numCells, 3);
  const overlap_removal::CellWinding cw =
      overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/-1.0);
  ASSERT_TRUE(cw.ok);
  EXPECT_EQ(cw.seedCasts, 1);
  const overlap_removal::EmitTopology et =
      overlap_removal::BuildEmitTopology(polys, cc, cw);
  ASSERT_TRUE(et.ok);
  ASSERT_EQ(et.keptPolygons.size(), 24u);
  ASSERT_EQ(et.outCycles.size(), 24u);
  // Rings: one per vert except the shared-edge verts, which carry one
  // ring per solid.
  ASSERT_EQ(et.ring2Vert.size(), 16u);
  int rings3 = 0;
  int rings7 = 0;
  for (const int v : et.ring2Vert) {
    rings3 += v == 3;
    rings7 += v == 7;
  }
  EXPECT_EQ(rings3, 2);
  EXPECT_EQ(rings7, 2);
  // Closed in ring space: each directed edge once, twin antiparallel.
  std::map<std::pair<int, int>, int> directed;
  for (const std::vector<int>& cyc : et.outCycles) {
    ASSERT_EQ(cyc.size(), 3u);
    for (size_t i = 0; i < cyc.size(); ++i) {
      ++directed[{cyc[i], cyc[(i + 1) % cyc.size()]}];
    }
  }
  for (const auto& [e, n] : directed) {
    EXPECT_EQ(n, 1) << e.first << "->" << e.second;
    const auto rev = directed.find({e.second, e.first});
    EXPECT_TRUE(rev != directed.end() && rev->second == 1)
        << e.first << "->" << e.second << " unpaired";
  }
  // Outward orientation, reversed-representation tri included: each
  // output tri's normal points away from its box's center.
  for (size_t k = 0; k < et.keptPolygons.size(); ++k) {
    const manifold::vec3 center = et.keptPolygons[k] < 12
                                      ? manifold::vec3(0.5, 0.5, 0.5)
                                      : manifold::vec3(1.5, 1.5, 0.5);
    const std::vector<int>& cyc = et.outCycles[k];
    const manifold::vec3 a = pos[et.ring2Vert[cyc[0]]];
    const manifold::vec3 b = pos[et.ring2Vert[cyc[1]]];
    const manifold::vec3 c = pos[et.ring2Vert[cyc[2]]];
    const manifold::vec3 n = la::cross(b - a, c - a);
    EXPECT_GT(la::dot(n, (a + b + c) / 3.0 - center), 0.0)
        << "kept polygon " << et.keptPolygons[k];
  }
}

TEST(OverlapRemoval, Step13OddKeptFanFailsEmitTopology) {
  // Three kept pages sharing one edge: inside-wedge twin pairing
  // cannot release the shared edge with an odd page count, and the
  // emit must fail CLOSED rather than emit an unpaired halfedge - in
  // release as ok = false, in MANIFOLD_DEBUG as the corresponding
  // assert (the driver catches it and falls back either way).
  Manifold::Impl impl;
  const std::vector<manifold::vec3> pos = {{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 1.0},
                                           {1.0, 0.0, 0.5},
                                           {-0.5, 0.5, 0.5},
                                           {-0.5, -0.5, 0.5}};
  const std::vector<overlap_removal::MergedPolygon> polys = {
      {{0, 1, 2}, 1, 0}, {{0, 1, 3}, 1, 1}, {{0, 1, 4}, 1, 2}};
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  overlap_removal::CellWinding cw;
  cw.ok = true;
  cw.seedCasts = 0;
  cw.winding.assign(cc.numCells, 1);  // every wedge "inside"
  cw.keep.assign(3, true);
  cw.flip.assign(3, false);
  bool failedClosed = false;
  try {
    const overlap_removal::EmitTopology et =
        overlap_removal::BuildEmitTopology(polys, cc, cw);
    failedClosed = !et.ok;
  } catch (const std::exception& e) {
    failedClosed = true;
    EXPECT_NE(std::string(e.what()).find("odd kept count"), std::string::npos)
        << e.what();
  }
  EXPECT_TRUE(failedClosed);
}

TEST(OverlapRemoval, Step13ConcaveSeedTargetUsesRealTriangulationEar) {
  // Seed targets for > 3-vert cycles come from the largest ear of a
  // REAL triangulation: a concave L-prism's hex faces (the largest
  // polygons, so they are tried first) have their vert-centroid
  // EXACTLY on the reflex corner - a guaranteed graze for a centroid
  // regression (and fan ears from the first vert can land outside a
  // concavity in general). The real-ear target is interior: one cast,
  // every face kept at winding 0|1.
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  const double xy[6][2] = {{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}};
  for (const auto& p : xy) pos.push_back({p[0], p[1], 0.0});  // 0-5 bottom
  for (const auto& p : xy) pos.push_back({p[0], p[1], 0.2});  // 6-11 top
  std::vector<overlap_removal::MergedPolygon> polys;
  polys.push_back({{0, 5, 4, 3, 2, 1}, 1, 0});    // bottom, outward -z
  polys.push_back({{6, 7, 8, 9, 10, 11}, 1, 1});  // top, outward +z
  for (int i = 0; i < 6; ++i) {                   // sides, outward
    const int j = (i + 1) % 6;
    polys.push_back({{i, j, 6 + j, 6 + i}, 1, 2 + i});
  }
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.numCells, 2);
  const overlap_removal::CellWinding cw =
      overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/-1.0);
  ASSERT_TRUE(cw.ok);
  EXPECT_EQ(cw.seedCasts, 1);  // the first (largest) target cast cleanly
  for (size_t p = 0; p < polys.size(); ++p) {
    EXPECT_TRUE(cw.keep[p]) << "polygon " << p;
    const int wF = cw.winding[cc.polySide2Cell[2 * p]];
    const int wB = cw.winding[cc.polySide2Cell[2 * p + 1]];
    EXPECT_EQ(std::min(wF, wB), 0) << "polygon " << p;
    EXPECT_EQ(std::max(wF, wB), 1) << "polygon " << p;
  }
}

TEST(OverlapRemoval, Step13SeedCastExhaustionFailsClosed) {
  // Every seed target grazes (an oversized epsilon hint makes each
  // arrival read as tangential): the cast budget caps at
  // kSeedCastMaxTargets and classification fails CLOSED - in release
  // as ok = false with seedCasts == 8, in MANIFOLD_DEBUG as the
  // exhaustion assert (the driver catches it and falls back).
  Manifold::Impl impl;
  std::vector<manifold::vec3> pos;
  std::vector<overlap_removal::MergedPolygon> polys;
  AppendCubePolys(0.0, 1.0, pos, polys);  // 12 tris > the 8-cast budget
  const overlap_removal::CellComplex cc =
      overlap_removal::BuildCellComplex(impl, polys, pos);
  ASSERT_EQ(cc.numCells, 2);
  bool failedClosed = false;
  try {
    const overlap_removal::CellWinding cw =
        overlap_removal::ClassifyCells(impl, polys, pos, cc, /*epsHint=*/1e6);
    failedClosed = !cw.ok;
    EXPECT_EQ(cw.seedCasts, 8);  // the kSeedCastMaxTargets cap
  } catch (const std::exception& e) {
    failedClosed = true;
    EXPECT_NE(std::string(e.what()).find("seed cast retries exhausted"),
              std::string::npos)
        << e.what();
  }
  EXPECT_TRUE(failedClosed);
}

TEST(OverlapRemoval, Step95UnifyArrangementVerts) {
  // Twins from two allocation paths 5 * eps apart, plus a new vert
  // 3 * eps from an original corner: one sweep at the nearby-crossing
  // radius (10 * eps) unifies both - smallest id wins, originals
  // before new - and remaps chords (endpoint + extras with id-dedup)
  // and on-edge lists (endpoint entries dropped).
  const Step10Fixture fx = MakeStep10Fixture();
  const double eps = 1e-6;
  // baseId = 4 verts: ids 4, 5, 6.
  const std::vector<manifold::vec3> newPos = {
      {0.5, 0.25, 0.0},         // 4
      {0.5 + 5e-6, 0.25, 0.0},  // 5: within 10 eps of 4 -> rep 4
      {3e-6, 0.0, 0.0}};        // 6: within 10 eps of corner A -> A
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const int eAB = Step10EdgeIndex(fx.edges, fx.A, fx.B);
  ASSERT_GE(eAB, 0);
  onEdgeLists[eAB].verts = {6, 5};
  onEdgeLists[eAB].ts = {0.1, 0.5};
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{std::min(5, fx.C), std::max(5, fx.C), fx.face, 99}, {4}, {0.4}}};
  const int changed =
      overlap_removal::UnifyArrangementVerts(fx.impl, newPos, fx.edges,
                                             onEdgeLists, chords, eps, {})
          .changed;
  EXPECT_EQ(changed, 2);  // 5 -> 4 and 6 -> A
  EXPECT_EQ(std::min(chords[0].edge.v0, chords[0].edge.v1), fx.C);
  EXPECT_EQ(std::max(chords[0].edge.v0, chords[0].edge.v1), 4);
  // The extra that became an endpoint id is dropped.
  EXPECT_TRUE(chords[0].extraVerts.empty());
  // On-edge list: 6 remapped to corner A (an endpoint: dropped);
  // 5 remapped to 4 (kept).
  ASSERT_EQ(onEdgeLists[eAB].verts.size(), 1u);
  EXPECT_EQ(onEdgeLists[eAB].verts[0], 4);
}

TEST(OverlapRemoval, Step95UnifyRecomputesAndResortsTs) {
  // Review finding: remapping an extra's id moves its consumed
  // POSITION by up to the merge radius, so stored ts go stale and
  // two adjacent extras can invert their true order - the partition
  // would then build a crossed sub-edge sequence. After unification,
  // extras must carry ts recomputed from the remapped positions and
  // be re-sorted.
  const Step10Fixture fx = MakeStep10Fixture();
  const double eps = 1e-2;  // merge radius 10 eps = 0.1
  const int v0 = std::min(fx.A, fx.B);
  const int v1 = std::max(fx.A, fx.B);
  const manifold::vec3 a = fx.impl.vertPos_[v0];
  const manifold::vec3 b = fx.impl.vertPos_[v1];
  auto along = [&](double t, double off) {
    manifold::vec3 p = a + t * (b - a);
    p.y += off;  // transverse to the x-axis edge, in the z = 0 face
    return p;
  };
  // baseId = 4: ids 4, 5, 6. Extra 5 (t 0.30) unifies with 4
  // (t 0.345, smaller id wins), jumping PAST extra 6 (t 0.32, held
  // 0.15 off-axis - outside both clusters in 3D, same projected t).
  const std::vector<manifold::vec3> newPos = {
      along(0.345, 0.0), along(0.30, 0.0), along(0.32, 0.15)};
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{v0, v1, fx.face, 99}, {5, 6}, {0.30, 0.32}}};
  const int changed =
      overlap_removal::UnifyArrangementVerts(fx.impl, newPos, fx.edges,
                                             onEdgeLists, chords, eps, {})
          .changed;
  EXPECT_EQ(changed, 1);  // 5 -> 4 only
  ASSERT_EQ(chords[0].extraVerts.size(), 2u);
  EXPECT_EQ(chords[0].extraVerts[0], 6);  // t 0.32 now precedes
  EXPECT_EQ(chords[0].extraVerts[1], 4);  // t 0.345
  ASSERT_EQ(chords[0].extraTs.size(), 2u);
  EXPECT_LT(chords[0].extraTs[0], chords[0].extraTs[1]);
  EXPECT_NEAR(chords[0].extraTs[1], 0.345, 1e-9);
}

TEST(OverlapRemoval, Step95UnifySnapsToNearestOriginal) {
  // Review finding: the original-vert snap picked the SMALLEST id
  // among in-radius candidates instead of the nearest - the
  // pipeline's convention everywhere else is nearest, ties to
  // smallest. A new vert sitting between two corners, closer to the
  // LARGER-id one, must snap to that nearer corner.
  const Step10Fixture fx = MakeStep10Fixture();
  const double eps = 0.08;  // snap radius 10 eps = 0.8 covers both
  const int lo = std::min(fx.A, fx.B);
  const int hi = std::max(fx.A, fx.B);
  // On the A-B edge (unit length): 0.3 from hi, 0.7 from lo.
  const manifold::vec3 p = fx.impl.vertPos_[hi] +
                           0.3 * (fx.impl.vertPos_[lo] - fx.impl.vertPos_[hi]);
  const std::vector<manifold::vec3> newPos = {p};  // id 4
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{std::min(4, fx.C), std::max(4, fx.C), fx.face, 99}, {}, {}}};
  const int changed =
      overlap_removal::UnifyArrangementVerts(fx.impl, newPos, fx.edges,
                                             onEdgeLists, chords, eps, {})
          .changed;
  EXPECT_EQ(changed, 1);
  EXPECT_EQ(std::min(chords[0].edge.v0, chords[0].edge.v1), std::min(hi, fx.C));
  EXPECT_EQ(std::max(chords[0].edge.v0, chords[0].edge.v1), std::max(hi, fx.C));
}

TEST(OverlapRemoval, Step95UnifyWidensSnapByPerVertRadius) {
  // The conditioned-radius carry (steps 6.5/9 -> 9.5): a vert whose
  // allocation was ill-conditioned widens ITS new-onto-original snap
  // beyond the 10 eps default; an unconditioned vert at the same
  // distance stays put. The applied displacement is reported for the
  // driver's output tolerance claim.
  const Step10Fixture fx = MakeStep10Fixture();
  const double eps = 0.01;  // 10 eps = 0.1 < the 0.3 separation
  const int lo = std::min(fx.A, fx.B);
  const int hi = std::max(fx.A, fx.B);
  const manifold::vec3 p = fx.impl.vertPos_[hi] +
                           0.3 * (fx.impl.vertPos_[lo] - fx.impl.vertPos_[hi]);
  const std::vector<manifold::vec3> newPos = {p};  // id 4
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(fx.edges.size());
  const std::vector<overlap_removal::NewEdgeWithExtras> chordsProto = {
      {{std::min(4, fx.C), std::max(4, fx.C), fx.face, 99}, {}, {}}};
  {  // Control: no conditioned radius - 0.3 is outside 10 eps.
    std::vector<overlap_removal::NewEdgeWithExtras> chords = chordsProto;
    const overlap_removal::UnifyResult r =
        overlap_removal::UnifyArrangementVerts(fx.impl, newPos, fx.edges,
                                               onEdgeLists, chords, eps, {});
    EXPECT_EQ(r.changed, 0);
    EXPECT_EQ(r.maxMove, 0.0);
  }
  std::vector<overlap_removal::NewEdgeWithExtras> chords = chordsProto;
  const overlap_removal::UnifyResult r = overlap_removal::UnifyArrangementVerts(
      fx.impl, newPos, fx.edges, onEdgeLists, chords, eps, {0.35});
  EXPECT_EQ(r.changed, 1);
  EXPECT_NEAR(r.maxMove, 0.3, 1e-12);  // |p - hi|, the applied snap
  EXPECT_EQ(std::min(chords[0].edge.v0, chords[0].edge.v1), std::min(hi, fx.C));
  EXPECT_EQ(std::max(chords[0].edge.v0, chords[0].edge.v1), std::max(hi, fx.C));
}

TEST(OverlapRemoval, Step9SnapRadiusReachesStep95) {
  // The driver handoff the review flagged as severable: a conditioned
  // step-9 allocation's radius must ride Step9Threading::newVertSnapR
  // into step 9.5 and widen its new-onto-original snap, composed the
  // way the driver composes them. An original vert sits 0.06 from the
  // crossing - outside the 10 eps = 0.05 default, inside the cluster's
  // 0.1 conditioned radius - so the snap happens ONLY if the radii
  // flow through.
  const double eps = 0.005;
  const double tolerance = eps;
  // A real impl supplies the original verts (id 4 at x = 0.56 is the
  // snap target; the rest are far).
  Manifold::Impl impl;
  const manifold::vec3 verts[6] = {{0.0, -5.0, 0.0}, {10.0, 0.0, 0.0},
                                   {0.0, 10.0, 0.0}, {-10.0, 0.0, 0.0},
                                   {0.56, 0.0, 0.0}, {-10.0, 5.0, 0.0}};
  for (const auto& v : verts) impl.vertPos_.push_back(v);
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const int baseId = static_cast<int>(impl.NumVert());  // 6
  // Two chords (new endpoints 6-9) crossing at (0.5, 0, 0).
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{baseId + 0, baseId + 1, 0, 9}, {}, {}},
      {{baseId + 2, baseId + 3, 0, 11}, {}, {}}};
  std::vector<manifold::vec3> pos = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.5, -0.5, 0.0}, {0.5, 0.5, 0.0}};
  // Hand-built cluster at the crossing with a conditioned radius of
  // 0.1 (production computation of this radius is pinned by
  // Step9ShallowCrossingCarriesConditionedRadius).
  const std::vector<overlap_removal::ChordCrossing> clusters = {
      {{0.5, 0.0, 0.0}, -1, {0, 1}, {0.5, 0.5}, 0.1}};
  const overlap_removal::Step9Threading threaded =
      overlap_removal::ResolveAndThreadClusters(impl, std::move(chords),
                                                std::move(pos), {}, clusters,
                                                {}, tolerance, eps);
  ASSERT_EQ(threaded.newVertPositions.size(), 5u);  // fresh crossing
  ASSERT_EQ(threaded.newVertSnapR.size(), 5u);
  EXPECT_EQ(threaded.newVertSnapR[4], 0.1);  // carried
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
  std::vector<overlap_removal::NewEdgeWithExtras> chords2 = threaded.chords;
  const overlap_removal::UnifyResult r = overlap_removal::UnifyArrangementVerts(
      impl, threaded.newVertPositions, edges, onEdgeLists, chords2, eps,
      threaded.newVertSnapR);
  // The crossing vert snapped onto original 4: both chords now thread
  // vert 4 where they threaded the fresh crossing id.
  EXPECT_GE(r.changed, 1);
  // The conditioned snap moved the crossing 0.06 onto the original -
  // PAST the 10 eps = 0.05 floor; maxMove must report it (the driver
  // folds this into the exported tolerance as the step-9.5 term).
  EXPECT_NEAR(r.maxMove, 0.06, 1e-9);
  bool threads4 = false;
  for (const overlap_removal::NewEdgeWithExtras& c : chords2) {
    for (const int v : c.extraVerts) {
      EXPECT_NE(v, baseId + 4);  // the fresh id must be gone
      if (v == 4) threads4 = true;
    }
  }
  EXPECT_TRUE(threads4) << "conditioned radius did not reach step 9.5";
}

TEST(OverlapRemoval, Step95RemapPastEndpointDropsFromLists) {
  // A unification representative can land PAST an edge endpoint. The
  // on-edge member at t = 0.98 unifies with a twin hanging beyond the
  // edge's far end; the globally nearest (member, original) pair is
  // (the twin, an original at x = 2.08) - NOT the edge endpoint, which
  // would drop via the endpoint check and never reach the range
  // filter (the first version of this fixture made that mistake and
  // discriminated nothing). The remapped entry's recomputed t = 1.04
  // must DROP from the on-edge list (it subdivides nothing) -
  // retaining it would hand PartitionFace a boundary sequence
  // stepping outside the face.
  Manifold::Impl impl;
  const manifold::vec3 verts[6] = {{0.0, 0.0, 0.0},   {2.0, 0.0, 0.0},
                                   {0.0, 10.0, 0.0},  {2.08, 0.0, 0.0},
                                   {10.0, 10.0, 0.0}, {-10.0, 5.0, 0.0}};
  for (const auto& v : verts) impl.vertPos_.push_back(v);
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  int e01 = -1;
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].v0 == 0 && edges[e].v1 == 1) e01 = static_cast<int>(e);
  }
  ASSERT_GE(e01, 0);
  const double eps = 0.01;  // 10 eps = 0.1
  // New vert 6 ON edge (0,1) at x = 1.96 (t = 0.98; 0.04 from endpoint
  // v1, 0.12 from v3); new vert 7 at x = 2.05 (0.05 from v1, 0.03 from
  // v3). 6 and 7 unify (0.09 apart, inside the 0.1 radius). Nearest
  // across members: (7 -> v3, 0.03) beats (6 -> v1, 0.04) ->
  // representative v3 at x = 2.08, PAST the edge's far end.
  const std::vector<manifold::vec3> newPos = {{1.96, 0.0, 0.0},
                                              {2.05, 0.0, 0.0}};
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
  onEdgeLists[e01].verts = {6};
  onEdgeLists[e01].ts = {0.98};
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{6, 7, 0, 1}, {}, {}}};
  const overlap_removal::UnifyResult r = overlap_removal::UnifyArrangementVerts(
      impl, newPos, edges, onEdgeLists, chords, eps, {});
  EXPECT_GE(r.changed, 1);
  EXPECT_EQ(chords[0].edge.v0, 3);  // premise: representative is v3
  EXPECT_EQ(chords[0].edge.v1, 3);
  // The remapped entry (now original 3, t = 1.04) is OFF the edge.
  EXPECT_TRUE(onEdgeLists[e01].verts.empty())
      << "representative past the endpoint stayed in the on-edge list";

  // MIRROR (chord-extras site): a long chord from v0 to new vert 7
  // threads 6 as an extra at t ~ 0.956; after 6 and 7 unify to v3,
  // the extra's recomputed t exceeds 1 (v3 sits past the remapped far
  // endpoint... the endpoint ALSO remapped to v3, so the extra drops
  // by the endpoint rule first). Pin the genuinely-past case instead:
  // chord v0 -> 6 with extra 7 (t ~ 1.046 against the remapped chord
  // v0 -> v3): extra 7 remaps to v3 = the chord's own new endpoint -
  // again endpoint-dropped. The chord-extras range filter is only
  // reachable when the EXTRA's representative differs from both
  // remapped endpoints yet projects outside; build exactly that: a
  // SECOND original v6 past v3, a third new vert 8 near it, chord
  // v0 -> 7 with extra 8.
  Manifold::Impl impl2;
  const manifold::vec3 verts2[7] = {
      {0.0, 0.0, 0.0},   {2.0, 0.0, 0.0},   {0.0, 10.0, 0.0}, {2.08, 0.0, 0.0},
      {10.0, 10.0, 0.0}, {-10.0, 5.0, 0.0}, {2.30, 0.0, 0.0}};
  for (const auto& v : verts2) impl2.vertPos_.push_back(v);
  // Two open tris referencing all 7 verts (6 used via tri2's slot).
  impl2.halfedge_.push_back(0, -1, -1);
  impl2.halfedge_.push_back(1, -1, -1);
  impl2.halfedge_.push_back(2, -1, -1);
  impl2.halfedge_.push_back(3, -1, -1);
  impl2.halfedge_.push_back(4, -1, -1);
  impl2.halfedge_.push_back(6, -1, -1);
  const std::vector<overlap_removal::Edge> edges2 =
      overlap_removal::EnumerateEdges(impl2);
  // New verts: 7 at x = 2.20 (the chord's far endpoint, between v3
  // and v6, snapping to NEITHER: 0.12 from v3, 0.10 from v6 - both
  // outside 10 eps with eps = 0.005, radius 0.05); 8 at x = 2.26 and
  // 9 at x = 2.305: 8 and 9 unify (0.045 < 0.05) and 9's nearest
  // original v6 (0.005) wins across members -> 8 remaps to v6 at
  // x = 2.30, past the chord's far endpoint (t = 2.30 / 2.20 > 1).
  const double eps2 = 0.005;
  const std::vector<manifold::vec3> newPos2 = {
      {2.20, 0.0, 0.0}, {2.26, 0.0, 0.0}, {2.305, 0.0, 0.0}};
  std::vector<overlap_removal::EdgeVertList> onEdge2(edges2.size());
  std::vector<overlap_removal::NewEdgeWithExtras> chords2 = {
      {{0, 7, 0, 1}, {8}, {2.26 / 2.20}}};
  const overlap_removal::UnifyResult r2 =
      overlap_removal::UnifyArrangementVerts(impl2, newPos2, edges2, onEdge2,
                                             chords2, eps2, {});
  EXPECT_GE(r2.changed, 1);
  EXPECT_EQ(chords2[0].edge.v1, 7);  // chord endpoint did NOT remap
  EXPECT_TRUE(chords2[0].extraVerts.empty())
      << "extra's representative past the chord endpoint stayed threaded";
}

TEST(OverlapRemoval, Step95ClusterSnapsToNearestOriginalAcrossMembers) {
  // Review finding (round 2): each member's snapTo is its own nearest
  // original, but the CLUSTER representative was picked by smallest
  // id across members - a farther small-id original could beat the
  // adjacent corner. Fixture: originals at x = 0 (id 0) and x = 0.08
  // (id 3); new verts at x = 0.035 (member nearest = id 0, d = 0.035)
  // and x = 0.07 (member nearest = id 3, d = 0.01) unify into one
  // cluster (0.035 apart < 10 eps = 0.05). Nearest across members is
  // id 3; the old id-priority rule picked id 0.
  Manifold::Impl impl;
  const manifold::vec3 verts[6] = {{0.0, 0.0, 0.0},   {10.0, 0.0, 0.0},
                                   {0.0, 10.0, 0.0},  {0.08, 0.0, 0.0},
                                   {10.0, 10.0, 0.0}, {-10.0, 5.0, 0.0}};
  for (const auto& v : verts) impl.vertPos_.push_back(v);
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const double eps = 0.005;
  const std::vector<manifold::vec3> newPos = {{0.035, 0.0, 0.0},
                                              {0.07, 0.0, 0.0}};  // ids 6, 7
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
  std::vector<overlap_removal::NewEdgeWithExtras> chords = {
      {{6, 7, 0, 1}, {}, {}}};
  const overlap_removal::UnifyResult r = overlap_removal::UnifyArrangementVerts(
      impl, newPos, edges, onEdgeLists, chords, eps, {});
  EXPECT_EQ(r.changed, 2);          // both members remap to the cluster rep
  EXPECT_EQ(chords[0].edge.v0, 3);  // nearest across members, NOT id 0
  EXPECT_EQ(chords[0].edge.v1, 3);  // (a collapsed chord: both ends remap)
  EXPECT_NEAR(r.maxMove, 0.045, 1e-12);  // member at 0.035 -> rep at 0.08
}

TEST(OverlapRemoval, Step13FoldedOppositeShellsDoNotCancel) {
  // Review findings (rounds 2-3): the folded-shell gate must evaluate
  // EDGE-CONNECTED components of folded polygons, not whole cells and
  // not vert-connected groups - a positive shell and an inverted twin
  // folded into the SAME cell (or merely touching at one snapped
  // vert) would otherwise net to zero signed volume and slip under
  // the area threshold, silently deleting the positive shell.
  // (a) Disjoint twin: two unit cubes, one reversed, one cell.
  {
    Manifold::Impl impl;
    std::vector<manifold::vec3> pos;
    std::vector<overlap_removal::MergedPolygon> polys;
    AppendCubePolys(0.0, 1.0, pos, polys);  // +1 volume
    const size_t firstReversed = polys.size();
    AppendCubePolys(3.0, 4.0, pos, polys);  // disjoint twin...
    for (size_t p = firstReversed; p < polys.size(); ++p) {
      std::reverse(polys[p].cycle.begin(), polys[p].cycle.end());
    }
    overlap_removal::CellComplex cells;
    cells.numCells = 1;
    cells.polySide2Cell.assign(2 * polys.size(), 0);  // everything folded
    EXPECT_TRUE(overlap_removal::FoldedCellsEncloseVolume(impl, polys, pos,
                                                          cells,
                                                          /*eps=*/1e-9));
  }
  // (b) Vertex-touch twin: the reversed cube's corner vert id is
  // REPLACED by the positive cube's coincident corner id (the
  // post-step-9.5 snapped-vert configuration), pinching the two
  // shells at one vert. Vert-connected grouping would merge them and
  // cancel; edge-connected grouping must keep them apart and trip.
  {
    Manifold::Impl impl;
    std::vector<manifold::vec3> pos;
    std::vector<overlap_removal::MergedPolygon> polys;
    AppendCubePolys(0.0, 1.0, pos, polys);
    const size_t firstReversed = polys.size();
    AppendCubePolys(1.0, 2.0, pos, polys);  // shares corner (1,1,1)
    int sharedA = -1, sharedB = -1;
    for (int v = 0; v < 8; ++v) {
      if (pos[v] == manifold::vec3(1.0, 1.0, 1.0)) sharedA = v;
      if (pos[8 + v] == manifold::vec3(1.0, 1.0, 1.0)) sharedB = 8 + v;
    }
    ASSERT_GE(sharedA, 0);
    ASSERT_GE(sharedB, 0);
    for (size_t p = firstReversed; p < polys.size(); ++p) {
      std::reverse(polys[p].cycle.begin(), polys[p].cycle.end());
      for (int& v : polys[p].cycle) {
        if (v == sharedB) v = sharedA;  // pinch at the shared corner
      }
    }
    overlap_removal::CellComplex cells;
    cells.numCells = 1;
    cells.polySide2Cell.assign(2 * polys.size(), 0);
    EXPECT_TRUE(overlap_removal::FoldedCellsEncloseVolume(impl, polys, pos,
                                                          cells,
                                                          /*eps=*/1e-9));
  }
}

TEST(OverlapRemoval, EpsilonFromScaleDelegatesToAlphaBudget) {
  // Pins exact behavioral EQUALITY between boolean2's EpsilonFromScale
  // and shared.h's AlphaBudgetEpsilon at every scale class (zero,
  // sub-unit, exact power of two, non-power, large) and budget. This
  // guards DRIFT (a diverging formula or changed default budget fails
  // here); the one-copy-in-the-tree rule itself is a review contract -
  // a reintroduced bit-identical local copy would pass this test.
  for (const double L : {0.0, 0.37, 1.0, 1.5, 1024.0, 7.3e15}) {
    EXPECT_EQ(boolean2::EpsilonFromScale(L), AlphaBudgetEpsilon(L));
    EXPECT_EQ(boolean2::EpsilonFromScale(L, 0), AlphaBudgetEpsilon(L, 0));
    EXPECT_EQ(boolean2::EpsilonFromScale(L, 7), AlphaBudgetEpsilon(L, 7));
  }
}

// Impl-typed view of a fixture Manifold for the Impl-to-Impl pipeline
// seams (RemoveOverlaps, MergeVertsEps, CheckSelfIntersection). Goes
// through the public MeshGL64 boundary + the Impl ctor - for an
// already-constructed Manifold whose construction sweep is complete
// and stable that round-trip is a fixed point (sorted stays sorted),
// pinned by the fresh-rebuild determinism test. The qualifier
// matters: the ctor re-runs CleanupTopology / RemoveDegenerates, so a
// fixture with collapsible micro-features (sub-epsilon slivers, the
// 775k ball below) can come back ALTERED - such fixtures must
// hand-build their Impl instead.
Manifold::Impl MakeImpl(const Manifold& m) {
  return Manifold::Impl(m.GetMeshGL64());
}

TEST(OverlapRemoval, Step1MergeLargeClusterConverges) {
  // A whole high-resolution sphere inside one eps-cluster at a large
  // coordinate: thousands of members collapse to one centroid. Sums
  // of many near-equal doubles are not bit-idempotent, so without the
  // collapsed-cluster skip the recomputed centroid can keep drifting
  // an ULP per pass toward the iteration cap (drift length depends on
  // the exact (count, value) pair, so this pins the CLASS - large
  // cluster, high coordinate - rather than one drift trace; the
  // structural guarantee is the skip itself).
  //
  Manifold ball =
      Manifold::Sphere(1e-9, 64).Translate({774996.8, 774996.8, 774996.8});
  ASSERT_EQ(ball.Status(), Manifold::Error::NoError);
  const MeshGL64 mesh = ball.GetMeshGL64();
  const size_t n = mesh.NumVert();
  ASSERT_GT(n, 1000u);
  // Hand-build the Impl from the exported geometry: at this
  // coordinate scale a mesh-ctor round-trip (MakeImpl and the public
  // Manifold(MeshGL64) alike) collapses the whole sub-epsilon sphere
  // during construction - the documented export-reconstruction
  // lossiness class - while the production driver receives the
  // evaluated leaf Impl with every vert intact. CreateHalfedges plus
  // bbox/epsilon is all the merge reads.
  Manifold::Impl ballImpl;
  ballImpl.vertPos_.resize_nofill(n);
  for (size_t i = 0; i < n; ++i) {
    ballImpl.vertPos_[i] = vec3(mesh.vertProperties[mesh.numProp * i + 0],
                                mesh.vertProperties[mesh.numProp * i + 1],
                                mesh.vertProperties[mesh.numProp * i + 2]);
  }
  Vec<ivec3> triProp;
  triProp.reserve(mesh.NumTri());
  for (size_t t = 0; t < mesh.NumTri(); ++t) {
    triProp.push_back(ivec3(mesh.triVerts[3 * t + 0], mesh.triVerts[3 * t + 1],
                            mesh.triVerts[3 * t + 2]));
  }
  ballImpl.CreateHalfedges(triProp);
  ASSERT_TRUE(ballImpl.IsManifold());
  ballImpl.CalculateBBox();
  ballImpl.SetEpsilon();
  const overlap_removal::MergeVertsResult r =
      overlap_removal::MergeVertsEps(ballImpl, 1e-3, nullptr);  // no throw
  EXPECT_EQ(r.mergedCount, static_cast<int>(n) - 1);
  EXPECT_LT(r.maxMove, 1e-6);  // members moved at most ~the sphere size
}

TEST(OverlapRemoval, Step13FoldedMembranePassesBentOpenFoldTrips) {
  // The gate's other two documented arms, directly. (a) A FLAT folded
  // membrane (coplanar open sheet) encloses no volume about its own
  // centroid and must PASS - tripping it would veto the legitimate
  // membrane drops the keep rule performs. (b) A macro-BENT open fold
  // has centroid-anchored cone volume O(area x bend depth) and must
  // TRIP - an open fold of real extent is arrangement damage, and
  // falling back is the safe side.
  const double eps = 1e-9;
  {  // (a) flat quad sheet, two coplanar triangles, folded
    Manifold::Impl impl;
    const std::vector<manifold::vec3> pos = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    const std::vector<overlap_removal::MergedPolygon> polys = {
        {{0, 1, 2}, 1, 0}, {{0, 2, 3}, 1, 1}};
    overlap_removal::CellComplex cells;
    cells.numCells = 1;
    cells.polySide2Cell.assign(2 * polys.size(), 0);
    EXPECT_FALSE(overlap_removal::FoldedCellsEncloseVolume(impl, polys, pos,
                                                           cells, eps));
  }
  {  // (b) the same sheet folded 90 degrees along the diagonal
    Manifold::Impl impl;
    const std::vector<manifold::vec3> pos = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 1}};  // bent flap
    const std::vector<overlap_removal::MergedPolygon> polys = {
        {{0, 1, 2}, 1, 0}, {{0, 2, 3}, 1, 1}};
    overlap_removal::CellComplex cells;
    cells.numCells = 1;
    cells.polySide2Cell.assign(2 * polys.size(), 0);
    EXPECT_TRUE(overlap_removal::FoldedCellsEncloseVolume(impl, polys, pos,
                                                          cells, eps));
  }
}

TEST(OverlapRemoval, Step1MergeReportsMaxMove) {
  // The step-1 merge's applied displacement feeds the driver's output
  // tolerance claim. One eps-pair (verts 0 and 3, 2^-11 apart so the
  // halving is exact in FP) merges to its centroid: each member moves
  // half the separation.
  const double h = 0.00048828125;  // 2^-11
  MeshGL64 m;
  m.numProp = 3;
  m.vertProperties = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, h};
  m.triVerts = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
  Manifold tet((MeshGL64(m)));
  ASSERT_EQ(tet.Status(), Manifold::Error::NoError);
  const overlap_removal::MergeVertsResult r =
      overlap_removal::MergeVertsEps(MakeImpl(tet), 1e-3, nullptr);
  EXPECT_EQ(r.mergedCount, 1);
  EXPECT_NEAR(r.maxMove, h / 2, 1e-15);
}

TEST(OverlapRemoval, Step1MergeConvergesWhenHigherIdSortsFirst) {
  // Convergence regression premise: a per-pass union COUNT is the
  // wrong convergence test.
  // A coincident (already-merged) cluster re-unites in every pass's
  // fresh union-find, and DisjointSets attaches the LARGER id under
  // the smaller - so whenever the cluster's higher-id vert
  // Morton-sorts FIRST, find(va) changed every pass, the count never
  // hit zero, and an assert-enabled build threw the iteration-cap
  // tripwire on a perfectly ordinary merge. Convergence is now a
  // position fixed point. This fixture puts the HIGHER id (3) at the
  // bbox-min corner so it sorts first: under the old rule it never
  // converged.
  const double h = 0.00048828125;  // 2^-11, below eps = 1e-3
  MeshGL64 m;
  m.numProp = 3;
  // vert 3 at the origin (Morton-first), vert 0 at (0,0,h); verts 1,2
  // swapped vs the canonical tetra so orientation stays positive.
  m.vertProperties = {0, 0, h, 0, 1, 0, 1, 0, 0, 0, 0, 0};
  m.triVerts = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
  Manifold tet((MeshGL64(m)));
  ASSERT_EQ(tet.Status(), Manifold::Error::NoError);
  const overlap_removal::MergeVertsResult r =
      overlap_removal::MergeVertsEps(MakeImpl(tet), 1e-3,
                                     nullptr);  // must not throw
  EXPECT_EQ(r.mergedCount, 1);
  EXPECT_NEAR(r.maxMove, h / 2, 1e-15);
}

// ---- Step 6.5 coplanar trace chords (docs/OverlapRemoval.md) ----

// Hand-built Impl holding two disjoint open triangles (halfedges
// unpaired; EnumerateEdges keeps only forward halfedges, and the
// trace pass guards missing edge indices), for coplanar-pair unit
// fixtures that the public constructor would reject or strip.
Manifold::Impl MakeTwoTriImpl(const manifold::vec3 t1[3],
                              const manifold::vec3 t2[3]) {
  Manifold::Impl impl;
  for (int i = 0; i < 3; ++i) impl.vertPos_.push_back(t1[i]);
  for (int i = 0; i < 3; ++i) impl.vertPos_.push_back(t2[i]);
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  return impl;
}

TEST(OverlapRemoval, Step1MergeNonPositiveMaxIterFailsClosed) {
  // A nonpositive iteration cap is handled as a degenerate input at
  // entry: zero merges in EVERY build config (previously an
  // assert-enabled build threw the convergence tripwire and release
  // would have merged with uninitialized labels). Genuine
  // non-convergence keeps the debug tripwire plus a release
  // fail-closed arm.
  const double h = 0.00048828125;  // 2^-11, below eps = 1e-3
  MeshGL64 m;
  m.numProp = 3;
  m.vertProperties = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, h};
  m.triVerts = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
  Manifold tet((MeshGL64(m)));
  ASSERT_EQ(tet.Status(), Manifold::Error::NoError);
  const overlap_removal::MergeVertsResult r = overlap_removal::MergeVertsEps(
      MakeImpl(tet), 1e-3, nullptr, /*maxIter=*/0);
  EXPECT_EQ(r.mergedCount, 0);
  EXPECT_EQ(r.maxMove, 0.0);
}

TEST(OverlapRemoval, Step6SharedVertPairStillPierces) {
  // PAIR inclusion: an edge-tri pair sharing a vert is NOT skipped
  // (the post-Boolean-merge shared-vert pierce case classic #289
  // step 6 would miss) - only the edge's own two faces skip. T2
  // shares vert 0 with T1 while its opposite edge pierces T1's
  // interior; a reintroduced shared-vert pair skip loses the event.
  Manifold::Impl impl;
  impl.vertPos_.push_back({0, 0, 0});       // v0, shared corner
  impl.vertPos_.push_back({1, 0, 0});       // v1 (T1)
  impl.vertPos_.push_back({0, 1, 0});       // v2 (T1)
  impl.vertPos_.push_back({0.3, 0.3, -1});  // v3 (T2)
  impl.vertPos_.push_back({0.3, 0.3, 1});   // v4 (T2)
  const int tris[2][3] = {{0, 1, 2}, {0, 3, 4}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const double eps = 1e-6;
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
  const std::vector<overlap_removal::TriVertList> onTriLists(impl.NumTri());
  const std::vector<overlap_removal::EdgeTriIntersection> events =
      overlap_removal::FindEdgeTriIntersections(impl, edges, onEdgeLists,
                                                onTriLists, eps);
  // Exactly the v3-v4 edge pierces T1's interior at (0.3, 0.3, 0).
  ASSERT_EQ(events.size(), 1u);
  EXPECT_EQ(events[0].triIdx, 0);
  const overlap_removal::Edge& e = edges[events[0].edgeIdx];
  EXPECT_EQ(e.v0, 3);
  EXPECT_EQ(e.v1, 4);
  EXPECT_NEAR(events[0].position.x, 0.3, 1e-12);
  EXPECT_NEAR(events[0].position.y, 0.3, 1e-12);
  EXPECT_NEAR(events[0].position.z, 0.0, 1e-12);
}

TEST(OverlapRemoval, Step4OnEdgeListStrictInteriorSortedApexSkipped) {
  // BuildOnEdgeVertLists' three contracted behaviors on one edge:
  // verts within eps of the OPEN segment list (sorted by t), verts
  // beyond eps do not, and a vert that neighbors BOTH endpoints (the
  // thin-tri apex) is skipped even when within eps.
  const double eps = 1e-3;
  Manifold::Impl impl;
  impl.vertPos_.push_back({0, 0, 0});            // v0 - edge start
  impl.vertPos_.push_back({1, 0, 0});            // v1 - edge end
  impl.vertPos_.push_back({0.5, 0.4 * eps, 0});  // v2 - APEX of (v0,v1,v2)
  // A second, distant triangle contributes the probe verts.
  impl.vertPos_.push_back({0.7, 0.5 * eps, 0});  // v3 - within eps, t=0.7
  impl.vertPos_.push_back({0.3, 0.5 * eps, 0});  // v4 - within eps, t=0.3
  impl.vertPos_.push_back({0.5, 2.0 * eps, 0});  // v5 - OUTSIDE the band
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  int e01 = -1;
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].v0 == 0 && edges[e].v1 == 1) e01 = static_cast<int>(e);
  }
  ASSERT_GE(e01, 0);
  const std::vector<overlap_removal::EdgeVertList> lists =
      overlap_removal::BuildOnEdgeVertLists(impl, edges, eps);
  // v4 (t=0.3) then v3 (t=0.7), sorted ascending; apex v2 skipped
  // despite sitting closer to the edge than either probe; v5 excluded.
  ASSERT_EQ(lists[e01].verts.size(), 2u);
  EXPECT_EQ(lists[e01].verts[0], 4);
  EXPECT_EQ(lists[e01].verts[1], 3);
  EXPECT_NEAR(lists[e01].ts[0], 0.3, 1e-9);
  EXPECT_NEAR(lists[e01].ts[1], 0.7, 1e-9);
}

TEST(OverlapRemoval, Step5OnTriListStrictInterior) {
  // BuildOnTriVertLists: a vert within eps of the plane and strictly
  // inside the barycentric simplex lists (all-positive bary); a vert
  // past the eps band does not; a vert over an EDGE of the triangle
  // (one barycentric pinned at zero) does not - strict interior only.
  const double eps = 1e-3;
  Manifold::Impl impl;
  impl.vertPos_.push_back({0, 0, 0});                // T1
  impl.vertPos_.push_back({1, 0, 0});                // T1
  impl.vertPos_.push_back({0, 1, 0});                // T1
  impl.vertPos_.push_back({0.3, 0.3, 0.5 * eps});    // v3 - interior, near
  impl.vertPos_.push_back({0.25, 0.25, 2.0 * eps});  // v4 - too far
  impl.vertPos_.push_back({0.5, 0.0, 0.5 * eps});    // v5 - over the edge
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const std::vector<overlap_removal::TriVertList> lists =
      overlap_removal::BuildOnTriVertLists(impl, eps);
  ASSERT_EQ(lists.size(), 2u);
  ASSERT_EQ(lists[0].verts.size(), 1u);
  EXPECT_EQ(lists[0].verts[0], 3);
  EXPECT_GT(lists[0].bary[0].x, 0.0);
  EXPECT_GT(lists[0].bary[0].y, 0.0);
  EXPECT_GT(lists[0].bary[0].z, 0.0);
}

TEST(OverlapRemoval, Step8ExtrasThreadedSortedAlongChord) {
  // AddInteriorVertsToNewEdges: on-tri verts lying on a chord's
  // interior thread onto it as extras, sorted by t along the chord
  // regardless of their order in the on-tri list; off-chord verts do
  // not thread.
  const double eps = 1e-3;
  Manifold::Impl impl;
  impl.vertPos_.push_back({0, 0, 0});      // T1
  impl.vertPos_.push_back({2, 0, 0});      // T1
  impl.vertPos_.push_back({0, 2, 0});      // T1
  impl.vertPos_.push_back({1.2, 0.1, 0});  // v3 - ON the chord, t~0.69
  impl.vertPos_.push_back({0.5, 0.1, 0});  // v4 - ON the chord, t~0.25
  impl.vertPos_.push_back({0.8, 0.4, 0});  // v5 - off the chord
  const int tris[2][3] = {{0, 1, 2}, {3, 4, 5}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  const int baseId = static_cast<int>(impl.NumVert());
  // A chord across T1 between two NEW verts at y = 0.1.
  const std::vector<manifold::vec3> newVertPositions = {{0.1, 0.1, 0},
                                                        {1.7, 0.1, 0}};
  std::vector<overlap_removal::PiercedNewEdge> chords(1);
  chords[0].v0 = baseId + 0;
  chords[0].v1 = baseId + 1;
  chords[0].triA = 0;
  chords[0].triB = 1;
  // On-tri list for T1 deliberately ordered far-then-near (v3 before
  // v4) to discriminate the sort.
  std::vector<overlap_removal::TriVertList> onTriLists(impl.NumTri());
  onTriLists[0].verts = {3, 4, 5};
  onTriLists[0].bary = {{0.35, 0.6, 0.05}, {0.7, 0.25, 0.05}, {0.4, 0.4, 0.2}};
  const std::vector<overlap_removal::NewEdgeWithExtras> out =
      overlap_removal::AddInteriorVertsToNewEdges(impl, newVertPositions,
                                                  chords, onTriLists, eps);
  ASSERT_EQ(out.size(), 1u);
  ASSERT_EQ(out[0].extraVerts.size(), 2u);
  EXPECT_EQ(out[0].extraVerts[0], 4);  // t ~ 0.25 first
  EXPECT_EQ(out[0].extraVerts[1], 3);  // t ~ 0.69 second
  EXPECT_LT(out[0].extraTs[0], out[0].extraTs[1]);
}

TEST(OverlapRemoval, Step65OppositeCornerSnapDoesNotSubdivideEdge) {
  // The opposite-corner exclusion: a crossing on T1's bottom edge that
  // corner-snaps to T1's own APEX (the obtuse face's third vert, 0.45
  // off the edge yet projecting to t ~ 0.475) must NOT thread the apex
  // onto the bottom edge - the apex is a whole altitude away; only the
  // chord uses the snapped id. Tolerance 0.5 makes the apex the
  // nearest in-radius corner for both crossings.
  const double eps = 1e-6;
  const double tolerance = 0.5;
  const manifold::vec3 t1[3] = {{0, 0, 0}, {4, 0, 0}, {1.9, 0.45, 0}};
  const manifold::vec3 t2[3] = {
      {1.5, -1, 0}, {1.9, 1, 0}, {2.5, -1, 0}};  // opposite winding
  Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
  const int apexId = 2;  // t1[2]
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {}, tolerance,
                                           eps);
  // Premise: the crossings snapped to the apex (a chord uses its id).
  bool chordUsesApex = false;
  for (const overlap_removal::PiercedNewEdge& ch : res.chords) {
    if (ch.v0 == apexId || ch.v1 == apexId) chordUsesApex = true;
  }
  ASSERT_TRUE(chordUsesApex) << "premise: crossing snapped to the apex";
  int bottomEdge = -1;
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].v0 == 0 && edges[e].v1 == 1) bottomEdge = static_cast<int>(e);
  }
  ASSERT_GE(bottomEdge, 0);
  for (const overlap_removal::OnEdgeAddition& a : res.onEdgeAdditions) {
    EXPECT_FALSE(a.edge == bottomEdge && a.vertId == apexId)
        << "apex threaded onto the edge it is the opposite corner of";
  }
}

TEST(OverlapRemoval, Step65SourceGatedDedupKeepsDistinctPoolVert) {
  // The source-gated trace dedup, pinned: an ill-conditioned crossing
  // (condR ~ 50 eps from a near-parallel pair) must NOT absorb a
  // pre-existing well-conditioned pool vert in the (eps, condR] band -
  // the pool entry's recorded radius (eps) gates the match to
  // min(condR, eps). A regression to a bare condR radius would weld
  // them. T2's long edge runs nearly parallel to T1's bottom edge
  // (sin ~ 0.02), crossing it at ~(1.04, 0); the preseeded pool vert
  // sits 2e-5 away (20 eps).
  const double eps = 1e-6;
  const double tolerance = eps;
  const manifold::vec3 t1[3] = {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
  const manifold::vec3 t2[3] = {{1, 0.0008, 0}, {5, -0.0792, 0}, {0.9, -2, 0}};
  Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
  const int baseId = static_cast<int>(impl.NumVert());  // 6: preseed id
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(
          impl, edges, he2e, {{1.04 + 2e-5, 0.0, 0.0}}, tolerance, eps);
  // The pool keeps the preseed plus TWO fresh crossings (near ~1.0 and
  // far ~1.04) - the far crossing allocates instead of welding.
  ASSERT_EQ(res.newVertPositions.size(), 3u);
  for (const overlap_removal::PiercedNewEdge& ch : res.chords) {
    EXPECT_NE(ch.v0, baseId) << "chord welded to the preseeded pool vert";
    EXPECT_NE(ch.v1, baseId) << "chord welded to the preseeded pool vert";
  }
  // The far crossing's conditioned radius rode along (~eps / 0.02).
  int farIdx = -1;
  for (size_t j = 0; j < res.newVertPositions.size(); ++j) {
    if (std::abs(res.newVertPositions[j].x - 1.04) < 1e-3 &&
        j > 0) {  // skip the preseed itself
      farIdx = static_cast<int>(j);
    }
  }
  ASSERT_GE(farIdx, 1);
  EXPECT_GT(res.newVertSnapR[farIdx], 4e-5);
  EXPECT_LT(res.newVertSnapR[farIdx], 6e-5);
}

TEST(OverlapRemoval, Step6SnapBandIsBareEps) {
  // The eps contract's event-identity arm: a pierce event snaps to an
  // existing vert only within BARE eps (the step-1 old-old scale) -
  // NOT the tolerance + eps allocation radius. An event 2 eps from a
  // corner must allocate (snapTo == -1); half an eps away must snap.
  // (Events in the (eps, 10 eps] band are unified by step 9.5; a
  // tolerance-scale radius here would drag pierce events onto far
  // verts and deform the arrangement.)
  const double eps = 1e-3;
  auto pierceAt = [&](double x0, double y0) {
    const manifold::vec3 t1[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const manifold::vec3 t2[3] = {{x0, y0, -1}, {x0, y0, 1}, {5, 5, 5}};
    Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
    const std::vector<overlap_removal::Edge> edges =
        overlap_removal::EnumerateEdges(impl);
    const std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
    const std::vector<overlap_removal::TriVertList> onTriLists(impl.NumTri());
    const std::vector<overlap_removal::EdgeTriIntersection> events =
        overlap_removal::FindEdgeTriIntersections(impl, edges, onEdgeLists,
                                                  onTriLists, eps);
    for (const overlap_removal::EdgeTriIntersection& x : events) {
      if (x.triIdx == 0) return x.snapTo;  // the pierce into T1
    }
    return -2;  // no event found
  };
  const double inv = 1.0 / std::sqrt(2.0);
  // 2 eps from corner (0,0,0): outside the band - allocate.
  EXPECT_EQ(pierceAt(2.0 * eps * inv, 2.0 * eps * inv), -1);
  // 0.5 eps from the corner: inside - snap to vert id 0.
  EXPECT_EQ(pierceAt(0.5 * eps * inv, 0.5 * eps * inv), 0);
}

TEST(OverlapRemoval, Step65PlaneGatePairWithinEpsOffsetIsFound) {
  // Review finding: the broad-phase boxes were built unpadded, so two
  // coplanar tris offset by half an eps along the normal - which PASS
  // the plane gate - had non-overlapping zero-thickness boxes and
  // were never tested. The pair must produce trace chords.
  const double eps = 1e-6;
  const double d = 0.5 * eps;
  const manifold::vec3 t1[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  const manifold::vec3 t2[3] = {
      {0.2, 0.2, d}, {0.2, 1.2, d}, {1.2, 0.2, d}};  // opposite winding
  Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {}, eps, eps);
  EXPECT_FALSE(res.chords.empty());
}

TEST(OverlapRemoval, Step65FullThroughCutQualifies) {
  // Pins the midpoint (not endpoint) interval qualification: a
  // full-through cut - the source edge entering AND exiting the other
  // face - has BOTH endpoints on that face's boundary, so any
  // endpoint-margin predicate would reject it; its midpoint is
  // interior and it must emit a chord between two NEW verts at the
  // boundary crossings.
  const double eps = 1e-9;
  const manifold::vec3 t1[3] = {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
  const manifold::vec3 t2[3] = {
      {-1, 1, 0}, {2, 5, 0}, {5, 1, 0}};  // opposite winding
  Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
  const int baseId = static_cast<int>(impl.NumVert());
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {}, eps, eps);
  // T2's edge (-1,1)->(5,1) crosses T1's left edge at (0,1) and its
  // hypotenuse at (3,1): both chord endpoints are NEW verts there.
  bool found = false;
  for (const overlap_removal::PiercedNewEdge& ch : res.chords) {
    if (ch.v0 < baseId || ch.v1 < baseId) continue;
    const manifold::vec3 p0 = res.newVertPositions[ch.v0 - baseId];
    const manifold::vec3 p1 = res.newVertPositions[ch.v1 - baseId];
    const manifold::vec3 lo = p0.x < p1.x ? p0 : p1;
    const manifold::vec3 hi = p0.x < p1.x ? p1 : p0;
    if (std::abs(lo.x - 0.0) < 1e-9 && std::abs(lo.y - 1.0) < 1e-9 &&
        std::abs(hi.x - 3.0) < 1e-9 && std::abs(hi.y - 1.0) < 1e-9) {
      found = true;
    }
  }
  EXPECT_TRUE(found) << "full-through chord (0,1)-(3,1) missing";
}

TEST(OverlapRemoval, Step65SnappedEndpointStillSubdividesSourceEdge) {
  // Regression premise: on-edge additions were emitted only for
  // ALLOCATED crossing endpoints. A crossing on a source edge that
  // SNAPS to a corner of the OTHER face (within the tolerance + eps /
  // conditioned corner radius, but farther than eps from the edge -
  // so step 4 never listed it) must still subdivide that source edge
  // at the snapped id, or the partition of the source face never sees
  // the cut and the sheets do not conform. T2's apex (1, 0.03) sits
  // near the interior of T1's bottom edge; the near crossing
  // (~x = 0.9985) snaps to the apex with tolerance 0.05, the far
  // crossing (~x = 1.059) allocates.
  const double eps = 1e-6;
  const double tolerance = 0.05;
  const manifold::vec3 t1[3] = {{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
  const manifold::vec3 t2[3] = {
      {1, 0.03, 0}, {5, -2, 0}, {0.9, -2, 0}};  // opposite winding
  Manifold::Impl impl = MakeTwoTriImpl(t1, t2);
  const int apexId = 3;  // t2[0]
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {}, tolerance,
                                           eps);
  // Premise: a chord with the apex as one endpoint exists.
  bool chordUsesApex = false;
  for (const overlap_removal::PiercedNewEdge& ch : res.chords) {
    if (ch.v0 == apexId || ch.v1 == apexId) chordUsesApex = true;
  }
  ASSERT_TRUE(chordUsesApex) << "premise: near crossing snapped to the apex";
  // The bug: no on-edge addition threading the apex onto T1's bottom
  // edge (verts 0 and 1).
  int bottomEdge = -1;
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].v0 == 0 && edges[e].v1 == 1) bottomEdge = static_cast<int>(e);
  }
  ASSERT_GE(bottomEdge, 0);
  bool apexOnBottomEdge = false;
  for (const overlap_removal::OnEdgeAddition& a : res.onEdgeAdditions) {
    if (a.edge == bottomEdge && a.vertId == apexId) {
      apexOnBottomEdge = true;
      EXPECT_NEAR(a.t, 0.9985 / 4.0, 0.01);
    }
  }
  EXPECT_TRUE(apexOnBottomEdge)
      << "snapped endpoint not threaded onto its source edge";
  // Conditioned-radius production (the eps contract's trace arm):
  // newVertSnapR is parallel to the pool, eps-floored, and the
  // allocated far crossing records its eps / sin(angle) conditioning
  // (the T2 edge meets T1's bottom edge at sin ~ 0.45 -> ~2.2 eps).
  ASSERT_EQ(res.newVertSnapR.size(), res.newVertPositions.size());
  for (const double r : res.newVertSnapR) EXPECT_GE(r, eps);
  ASSERT_GE(res.newVertPositions.size(), 1u);
  EXPECT_GT(res.newVertSnapR[0], 2.0 * eps);
  EXPECT_LT(res.newVertSnapR[0], 2.5 * eps);
}

TEST(OverlapRemoval, Step65PancakeQuadTraceChords) {
  // Zero-volume pancake: two sheets over the unit quad, triangulated
  // with DIFFERENT diagonals (top 0-2 at +z winding, bottom 1-3 at
  // -z), a valid closed manifold. Every cross-sheet pair overlaps;
  // the only genuine crossing is diag x anti-diag at the quad center.
  // Expect: ONE new vert at (0.5, 0.5, 0); trace chords spanning
  // center-to-corner, on cross-sheet pairs only (same-sheet
  // neighbors ride boundaries and emit nothing); and exactly two
  // on-edge additions - the center vert at t = 0.5 on BOTH diagonal
  // edges (the X case: one record per edge, one shared vert id).
  // Hand-built Impl: the public constructor strips zero-volume
  // components, but hull-class pancakes live ON a larger connected
  // surface; the unit test needs only vertPos_ and paired halfedges.
  Manifold::Impl impl;
  impl.vertPos_.push_back({0.0, 0.0, 0.0});
  impl.vertPos_.push_back({1.0, 0.0, 0.0});
  impl.vertPos_.push_back({1.0, 1.0, 0.0});
  impl.vertPos_.push_back({0.0, 1.0, 0.0});
  const int tris[4][3] = {{0, 1, 2}, {0, 2, 3}, {3, 1, 0}, {3, 2, 1}};
  for (const auto& t : tris) {
    for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
  }
  std::map<std::pair<int, int>, int> directedHe;
  for (int h = 0; h < 12; ++h) {
    directedHe[{impl.halfedge_.Start(h), impl.halfedge_.End(h)}] = h;
  }
  for (const auto& [uv, h] : directedHe) {
    impl.halfedge_.SetPair(h, directedHe.at({uv.second, uv.first}));
  }
  ASSERT_EQ(impl.NumTri(), 4u);
  const int c00 = 0;
  const int c10 = 1;
  const int c11 = 2;
  const int c01 = 3;
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const double eps = std::max(impl.epsilon_, 1e-12);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {},
                                           std::max(impl.tolerance_, eps), eps);
  const int baseId = static_cast<int>(impl.NumVert());
  ASSERT_EQ(res.newVertPositions.size(), 1u);
  const manifold::vec3 x = res.newVertPositions[0];
  EXPECT_NEAR(x.x, 0.5, 1e-9);
  EXPECT_NEAR(x.y, 0.5, 1e-9);
  EXPECT_NEAR(x.z, 0.0, 1e-9);
  const int center = baseId;
  // Chords: center-to-corner segments, each from two cross-sheet
  // pairs; never within one sheet.
  EXPECT_EQ(res.chords.size(), 8u);
  std::map<std::pair<int, int>, int> segCount;
  auto sheetZ = [&](int t) {  // winding sign (faceNormal_ not built)
    const manifold::vec3 a = impl.vertPos_[impl.halfedge_.Start(3 * t)];
    const manifold::vec3 b = impl.vertPos_[impl.halfedge_.Start(3 * t + 1)];
    const manifold::vec3 c = impl.vertPos_[impl.halfedge_.Start(3 * t + 2)];
    return la::cross(b - a, c - a).z;
  };
  for (const overlap_removal::PiercedNewEdge& ch : res.chords) {
    segCount[{std::min(ch.v0, ch.v1), std::max(ch.v0, ch.v1)}]++;
    EXPECT_LT(sheetZ(ch.triA) * sheetZ(ch.triB), 0.0)
        << "chord faces must be cross-sheet";
  }
  for (const int corner : {c00, c10, c11, c01}) {
    const std::pair<int, int> key{std::min(corner, center),
                                  std::max(corner, center)};
    EXPECT_EQ(segCount[key], 2) << "corner " << corner;
  }
  // On-edge additions: the center vert on both diagonal edges.
  ASSERT_EQ(res.onEdgeAdditions.size(), 2u);
  std::set<int> additionEdges;
  for (const overlap_removal::OnEdgeAddition& a : res.onEdgeAdditions) {
    EXPECT_EQ(a.vertId, center);
    EXPECT_NEAR(a.t, 0.5, 1e-9);
    additionEdges.insert(a.edge);
  }
  auto edgeIndexOf = [&](int u, int v) {
    for (size_t i = 0; i < edges.size(); ++i) {
      if (edges[i].v0 == std::min(u, v) && edges[i].v1 == std::max(u, v)) {
        return static_cast<int>(i);
      }
    }
    return -1;
  };
  EXPECT_TRUE(additionEdges.count(edgeIndexOf(c00, c11)));
  EXPECT_TRUE(additionEdges.count(edgeIndexOf(c10, c01)));
  EXPECT_GT(res.intervalsRejected, 0);
}

TEST(OverlapRemoval, Step65CoplanarNeighborsEmitNothing) {
  // A clean cube: every face's two tris are coplanar neighbors whose
  // clip intervals all ride the shared boundary - no trace chords, no
  // additions, so clean flat meshes keep the driver's early-exit.
  Manifold::Impl impl(Manifold::Cube({1, 1, 1}).GetMeshGL64());
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  const std::vector<int> he2e =
      overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
  const double eps = std::max(impl.epsilon_, 1e-12);
  const overlap_removal::TraceChordResult res =
      overlap_removal::CoplanarTraceChords(impl, edges, he2e, {},
                                           std::max(impl.tolerance_, eps), eps);
  EXPECT_TRUE(res.chords.empty());
  EXPECT_TRUE(res.onEdgeAdditions.empty());
  EXPECT_TRUE(res.newVertPositions.empty());
}

TEST(OverlapRemoval, Step65HazardFlagAntiAlignedOnly) {
  // The cancellation-hazard flag's PRODUCTION path (the gating pin
  // injects the flag directly; this pins where it comes from):
  // anti-aligned overlapping coplanar pairs flag their faces, while
  // same-oriented coplanar neighbors - every flat face's own tris -
  // must NOT flag, or the hole decomposition would re-gate ordinary
  // stamped faces (the plan-review round-5 class).
  {  // Pancake: two opposite-winding sheets over the unit quad.
    Manifold::Impl impl;
    impl.vertPos_.push_back({0.0, 0.0, 0.0});
    impl.vertPos_.push_back({1.0, 0.0, 0.0});
    impl.vertPos_.push_back({1.0, 1.0, 0.0});
    impl.vertPos_.push_back({0.0, 1.0, 0.0});
    const int tris[4][3] = {{0, 1, 2}, {0, 2, 3}, {3, 1, 0}, {3, 2, 1}};
    for (const auto& t : tris) {
      for (int k = 0; k < 3; ++k) impl.halfedge_.push_back(t[k], -1, -1);
    }
    const std::vector<overlap_removal::Edge> edges =
        overlap_removal::EnumerateEdges(impl);
    const std::vector<int> he2e =
        overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
    const double eps = 1e-9;
    const overlap_removal::TraceChordResult res =
        overlap_removal::CoplanarTraceChords(impl, edges, he2e, {}, eps, eps);
    // Every face overlaps an anti-aligned partner: all four flag.
    EXPECT_EQ(res.coplanarHazardFaces.size(), 4u);
  }
  {  // Clean cube: same-oriented coplanar neighbors only - no flags.
    Manifold::Impl impl(Manifold::Cube({1, 1, 1}).GetMeshGL64());
    const std::vector<overlap_removal::Edge> edges =
        overlap_removal::EnumerateEdges(impl);
    const std::vector<int> he2e =
        overlap_removal::BuildHalfedgeToEdgeIndex(impl, edges);
    const double eps = std::max(impl.epsilon_, 1e-12);
    const overlap_removal::TraceChordResult res =
        overlap_removal::CoplanarTraceChords(
            impl, edges, he2e, {}, std::max(impl.tolerance_, eps), eps);
    EXPECT_TRUE(res.coplanarHazardFaces.empty());
  }
}

TEST(OverlapRemoval, Step65AddVertsToOnEdgeLists) {
  // The trace-chord sibling of PropagateNewVertsToOnEdgeLists:
  // id-dedup against the existing list, then per-edge t re-sort.
  std::vector<overlap_removal::EdgeVertList> lists(2);
  lists[0].verts = {7};
  lists[0].ts = {0.5};
  const std::vector<overlap_removal::OnEdgeAddition> additions = {
      {0, 9, 0.8},
      {0, 8, 0.2},
      {0, 7, 0.5},  // 7 duplicates: skipped
      {1, 9, 0.4}};
  overlap_removal::AddVertsToOnEdgeLists(additions, lists);
  ASSERT_EQ(lists[0].verts.size(), 3u);
  EXPECT_EQ(lists[0].verts, (std::vector<int>{8, 7, 9}));
  EXPECT_EQ(lists[0].ts, (std::vector<double>{0.2, 0.5, 0.8}));
  ASSERT_EQ(lists[1].verts.size(), 1u);
  EXPECT_EQ(lists[1].verts[0], 9);
}

TEST(OverlapRemoval, Step7PropagateDropsOutOfRangeSnappedT) {
  // The (0, 1) guard on the recomputed propagate t: a pierce event
  // whose resolved vert SNAPPED far enough that its projection lands
  // outside the piercing edge's interior must not subdivide that edge
  // (the vert is still the chord's endpoint; only the subdivision is
  // skipped). Edge (0, 1) is short; the resolved vert sits past its
  // far end (t ~ 1.3).
  Manifold::Impl impl;
  const manifold::vec3 verts[4] = {
      {0.0, 0.0, 0.0}, {0.1, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.13, 0.0, 0.0}};
  for (const auto& v : verts) impl.vertPos_.push_back(v);
  impl.halfedge_.push_back(0, -1, -1);
  impl.halfedge_.push_back(1, -1, -1);
  impl.halfedge_.push_back(2, -1, -1);
  const std::vector<overlap_removal::Edge> edges =
      overlap_removal::EnumerateEdges(impl);
  int e01 = -1;
  for (size_t e = 0; e < edges.size(); ++e) {
    if (edges[e].v0 == 0 && edges[e].v1 == 1) e01 = static_cast<int>(e);
  }
  ASSERT_GE(e01, 0);
  // One event on edge (0,1) at s = 0.5, resolved (snapped) to vert 3,
  // which projects to t = 1.3 on that edge.
  const std::vector<overlap_removal::EdgeTriIntersection> ets = {
      {e01, 0, {0.05, 0.0, 0.0}, 0.5, {0.3, 0.3, 0.4}, 3}};
  const std::vector<int> et2v = {3};
  std::vector<overlap_removal::EdgeVertList> onEdgeLists(edges.size());
  overlap_removal::PropagateNewVertsToOnEdgeLists(impl, {}, ets, et2v, edges,
                                                  onEdgeLists);
  EXPECT_TRUE(onEdgeLists[e01].verts.empty())
      << "out-of-range snapped t subdivided the edge";
}

// White-box interior-pierce count via the internal checker (external
// linkage in the linked manifold library), used to assert the
// pierce-monotonicity contract that the public API does not expose.
// relTol = 1e-12 is tight for origin-scale fixtures and tolerated for
// the far-from-origin ones: at 1e4 coordinates the FP noise floor is
// of the same order, but no far-scale assertion depends on an exact
// count (they bound with GT 0 / LT / LE).
int InteriorPierces(const Manifold& m) {
  return overlap_removal::CheckSelfIntersection(MakeImpl(m), 1e-12)
      .interiorPierces;
}

// Bit-identical passthrough: the early-exit (empty chord list) and
// every fallback arm return the input Manifold exactly - every
// exported MeshGL64 field, not just volume/count agreement (review
// finding: loose bounds let a remesh, a tolerance change, or a
// metadata drop slip through).
void ExpectMeshGL64Identical(const Manifold& got, const Manifold& want) {
  const MeshGL64 g = got.GetMeshGL64();
  const MeshGL64 w = want.GetMeshGL64();
  EXPECT_EQ(g.numProp, w.numProp);
  EXPECT_EQ(g.NumVert(), w.NumVert());
  EXPECT_EQ(g.NumTri(), w.NumTri());
  EXPECT_EQ(g.vertProperties, w.vertProperties);
  EXPECT_EQ(g.triVerts, w.triVerts);
  EXPECT_EQ(g.tolerance, w.tolerance);
  EXPECT_EQ(g.mergeFromVert, w.mergeFromVert);
  EXPECT_EQ(g.mergeToVert, w.mergeToVert);
  EXPECT_EQ(g.runIndex, w.runIndex);
  EXPECT_EQ(g.runOriginalID, w.runOriginalID);
  EXPECT_EQ(g.runTransform, w.runTransform);
  EXPECT_EQ(g.runFlags, w.runFlags);
  EXPECT_EQ(g.faceID, w.faceID);
  EXPECT_EQ(g.halfedgeTangent, w.halfedgeTangent);
}

// Determinism of a fresh REBUILD: geometric output equality. Identity
// metadata (runOriginalID and friends) is construction-unique by
// design (see Manifold.MeshID), so the full-field helper above only
// applies to passthrough/fallback, where the SAME object comes back.
void ExpectMeshGL64GeometryIdentical(const Manifold& got,
                                     const Manifold& want) {
  const MeshGL64 g = got.GetMeshGL64();
  const MeshGL64 w = want.GetMeshGL64();
  EXPECT_EQ(g.numProp, w.numProp);
  EXPECT_EQ(g.vertProperties, w.vertProperties);
  EXPECT_EQ(g.triVerts, w.triVerts);
  EXPECT_EQ(g.tolerance, w.tolerance);
  EXPECT_EQ(g.mergeFromVert, w.mergeFromVert);
  EXPECT_EQ(g.mergeToVert, w.mergeToVert);
}

TEST(Manifold, RemoveSelfIntersectionsApi) {
  // API smoke: a clean cube has no self-intersections; output is the
  // input bit-identically (the empty-chord early-exit).
  Manifold cube = Manifold::Cube({1, 1, 1});
  Manifold cleaned = cube.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(cleaned, cube);
}

TEST(Manifold, RemoveSelfIntersectionsPropagatesErrorStatus) {
  // An errored input comes back with its status and contents intact.
  // This pins the OBSERVABLE contract only: manifold's invariant is
  // errored => empty, so the dedicated status short-circuit in
  // RemoveSelfIntersections is observationally equivalent to the
  // pipeline's own empty-input early-exit - no public fixture can
  // discriminate the arm (the guard's value is consistency with
  // every other member function).
  MeshGL64 open;  // a single triangle: not a closed manifold
  open.numProp = 3;
  open.vertProperties = {0, 0, 0, 1, 0, 0, 0, 1, 0};
  open.triVerts = {0, 1, 2};
  Manifold invalid((MeshGL64(open)));
  ASSERT_NE(invalid.Status(), Manifold::Error::NoError);
  ASSERT_TRUE(invalid.IsEmpty());  // the errored => empty invariant
  Manifold cleaned = invalid.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), invalid.Status());
  ExpectMeshGL64Identical(cleaned, invalid);
}

TEST(Manifold, RemoveSelfIntersectionsCancelBeforeRun) {
  // A pre-cancelled ExecutionContext stops the pipeline at the ENTRY
  // poll - before even the input pierce count - and the cancellation
  // is OBSERVABLE: an empty result carrying Error::Cancelled, exactly
  // like a cancelled boolean - never a silent input-return, which
  // would be indistinguishable from "nothing to do". Discriminates
  // poll removal: without the IsCancelled checks this pierce-free
  // input early-exits to a bit-identical NoError passthrough.
  Manifold cube = Manifold::Cube({1, 1, 1});
  ExecutionContext ctx;
  ctx.Cancel();
  Manifold cleaned = cube.WithContext(ctx).RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::Cancelled);
  EXPECT_TRUE(cleaned.IsEmpty());
}

TEST(Manifold, RemoveSelfIntersectionsCancelBeforeRunEmptyInput) {
  // The entry poll precedes even the empty-input exit: a pre-cancelled
  // run is observable as Cancelled for EVERY input class. Discriminates
  // re-ordering the empty exit ahead of the poll (which would return
  // the empty input as NoError).
  Manifold empty;
  ASSERT_TRUE(empty.IsEmpty());
  ASSERT_EQ(empty.Status(), Manifold::Error::NoError);
  ExecutionContext ctx;
  ctx.Cancel();
  Manifold cleaned = empty.WithContext(ctx).RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::Cancelled);
  EXPECT_TRUE(cleaned.IsEmpty());
}

TEST(Manifold, RemoveSelfIntersectionsUncancelledContextRuns) {
  // The ctx plumbing must not disturb an uncancelled run: same
  // bit-identical clean-input passthrough as without a context.
  Manifold cube = Manifold::Cube({1, 1, 1});
  ExecutionContext ctx;
  Manifold cleaned = cube.WithContext(ctx).RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(cleaned, cube);
  // The D6 contract's public face: RSI schedules no phases, so the
  // ctx reads as trivially complete. (Progress() also reads 1.0 for
  // COMPLETED accounting, so this pins the convention, not full D6
  // discrimination - recorded in the architecture review.)
  EXPECT_EQ(ctx.Progress(), 1.0);
}

TEST(Manifold, RemoveSelfIntersectionsCancelBeforeLazyEval) {
  // A pre-cancelled ctx on a LAZY CSG input: GetCsgLeafNode(ctx)
  // evaluates the boolean under the cancelled ctx, so the member's
  // status-propagation arm surfaces Cancelled before the pipeline
  // ever runs. The pipeline's own entry poll is a CO-witness (with
  // the GetCsgLeafNode handoff removed, the entry poll still yields
  // Cancelled), so this pins the lazy-input cancellation OUTCOME
  // end-to-end, not the CSG arm specifically - removing all ctx
  // plumbing makes this a NoError weld.
  Manifold lazy = Manifold::Cube({1, 1, 1}) +
                  Manifold::Cube({1, 1, 1}).Translate({0.5, 0.5, 0.5});
  ExecutionContext ctx;
  ctx.Cancel();
  Manifold cleaned = lazy.WithContext(ctx).RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::Cancelled);
  EXPECT_TRUE(cleaned.IsEmpty());
}

TEST(Manifold, RemoveSelfIntersectionsCleanInputUnchanged) {
  // A mesh with no self-intersections passes through bit-identically.
  Manifold sphere = Manifold::Sphere(1.0, 32);
  Manifold cleaned = sphere.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(cleaned, sphere);
}

TEST(Manifold, RemoveSelfIntersectionsBooleanResult) {
  // The common case: a clean Boolean output round-trips bit-
  // identically (no pierce events, no coplanar trace chords - the
  // empty-chord early-exit). (HullMask and SelfIntersect below cover
  // the actually-piercing paths.)
  Manifold a = Manifold::Cube({2, 2, 2}, true);
  Manifold b = Manifold::Cube({2, 2, 2}, true)
                   .Translate({1, 0.5, 0.3})
                   .Rotate(15, 30, 7);
  Manifold result = a + b;
  EXPECT_EQ(result.Status(), Manifold::Error::NoError);
  ASSERT_EQ(InteriorPierces(result), 0);  // premise: clean input
  Manifold cleaned = result.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(cleaned, result);
}

TEST(Manifold, RemoveSelfIntersectionsHullMaskFixture) {
  // Real-world adversarial fixture: hull-body Subtract hull-mask. The
  // body is THREE disjoint hulls (a trimaran); the mask grazes all of
  // them tangentially and Boolean3 produces a self-intersecting
  // manifold (~31 interior pierces). The tangent-degenerate contacts
  // fold two of the three untouched-shell cell complexes (a k = 1 rim
  // chain on one, near-tangent radial fans on the other), which would
  // silently delete those hulls (volume -> 1/3) - the folded-shell
  // volume gate detects this and falls back to the input unchanged.
  // Actually resolving this fixture (0 pierces, all three hulls kept)
  // needs exact predicates at the degenerate contacts; see
  // docs/OverlapRemoval.md "Known limitations". The same geometry
  // far from the origin resolves cleanly (its scale-derived eps
  // absorbs the cluster) - see RemoveSelfIntersectionsFarFromOrigin.
  Manifold body = Manifold(ReadTestMeshGL64OBJ("hull-body.obj"));
  Manifold mask = Manifold(ReadTestMeshGL64OBJ("hull-mask.obj"));
  Manifold result = body - mask;
  EXPECT_EQ(result.Status(), Manifold::Error::NoError);
  ASSERT_GT(InteriorPierces(result), 0);  // premise: input actually pierces
  Manifold cleaned = result.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  // Fail-closed fallback: bit-identical input passthrough.
  ExpectMeshGL64Identical(cleaned, result);
  EXPECT_EQ(InteriorPierces(cleaned), InteriorPierces(result));
}

TEST(Manifold, RemoveSelfIntersectionsSelfIntersectFixture) {
  // self_intersect: Add of two interpenetrating ovoids, dense interior
  // pierces. This is the dense-sliver fallback class - the pipeline
  // cannot cleanly reduce it, so it must fall back to the input
  // BIT-IDENTICALLY (the fail-closed contract; a partial rebuild
  // slipping out with merely-monotonic pierces would break it).
  //
  // test_main.cpp sets ManifoldParams().processOverlaps = false for
  // stricter validation; that enables a CCW-check assertion in Boolean3's
  // Triangulate that fires on this fixture. Restore the default here.
  ManifoldParamGuard guard;
  ManifoldParams().processOverlaps = true;
  std::filesystem::path file(__FILE__);
  std::filesystem::path modelsDir = file.parent_path() / "models";
  std::ifstream fa(modelsDir / "self_intersectA.obj");
  std::ifstream fb(modelsDir / "self_intersectB.obj");
  Manifold a = Manifold::ReadOBJ(fa);
  Manifold b = Manifold::ReadOBJ(fb);
  ASSERT_EQ(a.Status(), Manifold::Error::NoError);
  ASSERT_EQ(b.Status(), Manifold::Error::NoError);
  Manifold result = a + b;
  EXPECT_EQ(result.Status(), Manifold::Error::NoError);
  ASSERT_GT(InteriorPierces(result), 0);  // premise: input genuinely pierces
  Manifold cleaned = result.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(cleaned, result);
  EXPECT_EQ(InteriorPierces(cleaned), InteriorPierces(result));
}

// Appends an axis-aligned box [lo, hi] as 8 fresh verts + 12 tris to
// a raw MeshGL64 - a pure concatenation. (Manifold::Compose runs a
// real boolean union here and would resolve touching faces before
// RemoveSelfIntersections ever sees them; the constructor likewise
// strips standalone zero-volume components, which is why the pancake
// class has no standalone driver fixture - the hull fixture covers it
// end to end.)
void AppendBoxToMesh(MeshGL64& m, const vec3& lo, const vec3& hi) {
  const uint64_t b = m.NumVert();
  for (int i = 0; i < 8; ++i) {
    m.vertProperties.push_back((i & 1) ? hi.x : lo.x);
    m.vertProperties.push_back((i & 2) ? hi.y : lo.y);
    m.vertProperties.push_back((i & 4) ? hi.z : lo.z);
  }
  const uint64_t tris[12][3] = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6},
                                {0, 1, 5}, {0, 5, 4}, {2, 7, 3}, {2, 6, 7},
                                {0, 4, 6}, {0, 6, 2}, {1, 7, 5}, {1, 3, 7}};
  for (const auto& t : tris) {
    m.triVerts.push_back(b + t[0]);
    m.triVerts.push_back(b + t[1]);
    m.triVerts.push_back(b + t[2]);
  }
}

TEST(Manifold, RemoveSelfIntersectionsInteriorIslandResolves) {
  // A shell stamping through the INTERIOR of a single large face (its
  // footprint touching no face boundary) is the free-island class. After
  // the hole-aware partition, the annulus is decomposed into triangles and
  // the pipeline resolves correctly. The stamp box straddles the bottom
  // face of a unit cube: the intersection curve is a chord loop strictly
  // inside that face.
  //
  // Cube 1.0 + stamp_outside 0.2x0.2x0.3 = 1.012 (winding-union volume).
  MeshGL64 m;
  m.numProp = 3;
  AppendBoxToMesh(m, {0, 0, 0}, {1, 1, 1});
  // Strictly inside the bottom face's y > x tri, poking through.
  AppendBoxToMesh(m, {0.15, 0.55, -0.3}, {0.35, 0.75, 0.3});
  Manifold input((MeshGL64(m)));
  ASSERT_EQ(input.Status(), Manifold::Error::NoError);
  ASSERT_GT(InteriorPierces(input), 0);  // premise: genuinely pierces
  const std::optional<Manifold::Impl> out =
      overlap_removal::RemoveOverlaps(MakeImpl(input), 1e-3, nullptr);
  ASSERT_TRUE(out.has_value());
  ASSERT_EQ(out->status_, Manifold::Error::NoError);
  // Geometric checks via the public mesh boundary.
  Manifold cleaned(GetMeshGLImpl<double, uint64_t>(*out, -1));
  ASSERT_EQ(cleaned.Status(), Manifold::Error::NoError);
  EXPECT_EQ(InteriorPierces(cleaned), 0);
  EXPECT_EQ(cleaned.Decompose().size(), 1u);  // one connected component
  EXPECT_NEAR(cleaned.Volume(), 1.012, 1e-4);
  // Idempotence: a second RemoveOverlaps on the rebuilt result returns
  // nullopt (clean input early-exit) or the same volume.
  const std::optional<Manifold::Impl> out2 =
      overlap_removal::RemoveOverlaps(*out, 1e-3, nullptr);
  if (out2.has_value()) {
    Manifold cleaned2(GetMeshGLImpl<double, uint64_t>(*out2, -1));
    EXPECT_NEAR(cleaned2.Volume(), 1.012, 1e-4);
  }
  // No pierces in either case.
}

TEST(Manifold, RemoveSelfIntersectionsInteriorIslandPublicAPI) {
  // Public-API arm for the interior island fixture:
  // input.RemoveSelfIntersections() uses inferred eps
  // (AlphaBudgetEpsilon(bbox_scale, 1000) ~2.75e-12 for a unit-scale mesh). The
  // chord intersection points for this fixture land at machine-precision
  // coordinates, which the pipeline can detect at any positive eps, so the
  // inferred eps is sufficient.
  MeshGL64 m;
  m.numProp = 3;
  AppendBoxToMesh(m, {0, 0, 0}, {1, 1, 1});
  AppendBoxToMesh(m, {0.15, 0.55, -0.3}, {0.35, 0.75, 0.3});
  Manifold input((MeshGL64(m)));
  ASSERT_EQ(input.Status(), Manifold::Error::NoError);
  ASSERT_GT(InteriorPierces(input), 0);  // premise: genuinely pierces
  Manifold cleaned = input.RemoveSelfIntersections();
  ASSERT_EQ(cleaned.Status(), Manifold::Error::NoError);
  EXPECT_EQ(InteriorPierces(cleaned), 0);
  EXPECT_NEAR(cleaned.Volume(), 1.012, 1e-4);
}

TEST(Manifold, RemoveSelfIntersectionsToleranceCoversMergeDisplacement) {
  // The one output-tolerance formula term nothing else discriminated:
  // the driver folding the MEASURED step-1 merge displacement into the
  // exported tolerance. A 25-vert eps-chain strip on the cube's top
  // face (consecutive spacing 0.9 eps, so the whole strip transitively
  // merges to its centroid) moves its END verts ~10.8 eps - past the
  // 10 eps floor - while a smaller box glued onto the x = 1 face far
  // away welds through the step-6.5 trace path (the proven success
  // class; transversal raw-shell placements all land on documented
  // drop/gate edges). Dropping merged.maxMove from the formula exports
  // exactly the 10 eps floor and fails the bound below.
  const double eps = 1e-3;
  const int kChain = 25;
  const double spacing = 0.9 * eps;
  const double x0 = 0.5 - 0.5 * (kChain - 1) * spacing;
  MeshGL64 m;
  m.numProp = 3;
  AppendBoxToMesh(m, {0, 0, 0}, {1, 1, 1});  // verts 0-7
  // Replace the top face (z = 1: verts 4,5,7,6; table tris {4,5,7},
  // {4,7,6} = entries 2,3) with fans around the chain strip.
  m.triVerts.erase(m.triVerts.begin() + 6, m.triVerts.begin() + 12);
  const uint64_t c0 = m.NumVert();
  for (int i = 0; i < kChain; ++i) {
    m.vertProperties.push_back(x0 + i * spacing);
    m.vertProperties.push_back(0.5);
    m.vertProperties.push_back(1.0);
  }
  auto tri = [&](uint64_t a, uint64_t b, uint64_t c) {
    m.triVerts.push_back(a);
    m.triVerts.push_back(b);
    m.triVerts.push_back(c);
  };
  // South of the strip: fan from corner 4 (CCW seen from +z).
  tri(4, 5, c0 + kChain - 1);
  for (int i = kChain - 1; i > 0; --i) tri(4, c0 + i, c0 + i - 1);
  // North: fan from corner 6.
  for (int i = 0; i + 1 < kChain; ++i) tri(6, c0 + i, c0 + i + 1);
  tri(6, c0 + kChain - 1, 7);
  tri(4, c0, 6);               // west cap
  tri(5, 7, c0 + kChain - 1);  // east cap
  AppendBoxToMesh(m, {1.0, 0.25, 0.25}, {1.5, 0.75, 0.75});
  Manifold input((MeshGL64(m)));
  ASSERT_EQ(input.Status(), Manifold::Error::NoError);
  ASSERT_EQ(input.Decompose().size(), 2u);  // premise: two raw shells
  const std::optional<Manifold::Impl> out =
      overlap_removal::RemoveOverlaps(MakeImpl(input), eps, nullptr);
  // Premise: the SUCCESS path - a rebuilt Impl, not the nullopt
  // fallback (early-exit and every fallback would hand the caller
  // the 2-component input back bit-identically).
  ASSERT_TRUE(out.has_value());
  ASSERT_EQ(out->status_, Manifold::Error::NoError);
  // The tolerance claim is asserted on the Impl field DIRECTLY - the
  // exact value the pipeline computed, before any export flooring.
  // The chain ends moved (kChain - 1) / 2 * spacing = 10.8 eps to the
  // strip centroid; the claimed tolerance must cover it (the 10 eps
  // floor alone is 0.010).
  const double expectedMove = 0.5 * (kChain - 1) * spacing;
  ASSERT_GT(expectedMove, 10.0 * eps);  // fixture premise
  EXPECT_GE(out->tolerance_, expectedMove * 0.97);
  // ...and no wider: the formula takes the MAX of its terms, all of
  // which this fixture bounds (an over-wide claim - e.g. a blanket
  // conditioned-band term - would over-weld downstream consumers).
  EXPECT_LE(out->tolerance_, expectedMove * 1.03);
  // Geometric assertions go through the public mesh boundary, as a
  // user would see the result.
  Manifold cleaned(GetMeshGLImpl<double, uint64_t>(*out, -1));
  ASSERT_EQ(cleaned.Status(), Manifold::Error::NoError);
  ASSERT_EQ(cleaned.Decompose().size(), 1u);  // the glued shells welded
  EXPECT_EQ(InteriorPierces(cleaned), 0);
  EXPECT_NEAR(cleaned.Volume(), 1.0 + 0.5 * 0.5 * 0.5, 1e-6);
}

TEST(Manifold, RemoveSelfIntersectionsGluedBoxes) {
  // Two boxes glued face to face, built by raw mesh concatenation.
  // (a) Equal faces: step 1's eps-merge unifies the duplicated corner
  //     verts, the touching faces become exactly-equal opposite
  //     cycles, but every clip interval rides the boundary - no
  //     chords at all, and the early-exit returns the input
  //     unchanged.
  MeshGL64 mEq;
  mEq.numProp = 3;
  AppendBoxToMesh(mEq, {0, 0, 0}, {1, 1, 1});
  AppendBoxToMesh(mEq, {1, 0, 0}, {2, 1, 1});
  Manifold glued(mEq);
  ASSERT_EQ(glued.Status(), Manifold::Error::NoError);
  Manifold cleanedEq = glued.RemoveSelfIntersections();
  EXPECT_EQ(cleanedEq.Status(), Manifold::Error::NoError);
  EXPECT_EQ(InteriorPierces(cleanedEq), 0);
  EXPECT_NEAR(cleanedEq.Volume(), 2.0, 1e-9);
  ExpectMeshGL64Identical(cleanedEq, glued);  // early-exit, bit-identical
  // (b) A smaller box glued onto a larger face: the coincident
  //     interior wall region separates winding 1|1 and DROPS - the
  //     output is the winding-faithful welded solid (one component).
  MeshGL64 mMix;
  mMix.numProp = 3;
  AppendBoxToMesh(mMix, {0, 0, 0}, {1, 1, 1});
  AppendBoxToMesh(mMix, {1.0, 0.25, 0.25}, {1.5, 0.75, 0.75});
  Manifold mixed(mMix);
  ASSERT_EQ(mixed.Status(), Manifold::Error::NoError);
  ASSERT_EQ(mixed.Decompose().size(), 2u);  // premise: two components in
  Manifold cleanedMix = mixed.RemoveSelfIntersections();
  EXPECT_EQ(cleanedMix.Status(), Manifold::Error::NoError);
  EXPECT_EQ(InteriorPierces(cleanedMix), 0);
  EXPECT_NEAR(cleanedMix.Volume(), 1.125, 1e-9);
  EXPECT_EQ(cleanedMix.Decompose().size(), 1u);  // welded
}

TEST(Manifold, RemoveSelfIntersectionsEmptyInput) {
  // Empty input round-trips to empty without error.
  Manifold empty;
  Manifold cleaned = empty.RemoveSelfIntersections();
  EXPECT_TRUE(cleaned.IsEmpty());
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
}

TEST(Manifold, RemoveSelfIntersectionsIdempotent) {
  // Repeated passes never regress. At origin scale the hull fixture
  // takes the folded-shell fallback, so pass 2 of an identical input
  // must be the identical fallback (fixed point, bit-identical). On
  // the success path (the same fixture at 1e4) a second pass must
  // hold the monotonicity contract: no new pierces, valid status.
  Manifold body = Manifold(ReadTestMeshGL64OBJ("hull-body.obj"));
  Manifold mask = Manifold(ReadTestMeshGL64OBJ("hull-mask.obj"));
  Manifold input = body - mask;
  ASSERT_GT(InteriorPierces(input), 0);
  Manifold once = input.RemoveSelfIntersections();
  ASSERT_EQ(once.Status(), Manifold::Error::NoError);
  EXPECT_LE(InteriorPierces(once), InteriorPierces(input));
  Manifold twice = once.RemoveSelfIntersections();
  EXPECT_EQ(twice.Status(), Manifold::Error::NoError);
  ExpectMeshGL64Identical(twice, once);

  Manifold far = input.Translate({1e4, 1e4, 1e4});
  Manifold farOnce = far.RemoveSelfIntersections();
  ASSERT_EQ(farOnce.Status(), Manifold::Error::NoError);
  ASSERT_LT(InteriorPierces(farOnce), InteriorPierces(far));  // success path
  Manifold farTwice = farOnce.RemoveSelfIntersections();
  EXPECT_EQ(farTwice.Status(), Manifold::Error::NoError);
  EXPECT_LE(InteriorPierces(farTwice), InteriorPierces(farOnce));
  EXPECT_GT(farTwice.Volume(), 0);
}

TEST(Manifold, RemoveSelfIntersectionsDeterministic) {
  // The pipeline is deterministic: repeated runs on the same input
  // produce identical output. The origin-scale hull exercises the
  // fallback path; the far-from-origin variant exercises the FULL
  // rebuild (the fallback runs alone would pin only
  // trivial identity, not pipeline-ordering determinism).
  Manifold body = Manifold(ReadTestMeshGL64OBJ("hull-body.obj"));
  Manifold mask = Manifold(ReadTestMeshGL64OBJ("hull-mask.obj"));
  Manifold result = body - mask;
  ASSERT_GT(InteriorPierces(result), 0);  // premise: pipeline does real work
  Manifold first = result.RemoveSelfIntersections();
  ASSERT_EQ(first.Status(), Manifold::Error::NoError);
  for (int i = 0; i < 4; ++i) {
    Manifold again = result.RemoveSelfIntersections();
    ExpectMeshGL64Identical(again, first);
  }
  Manifold far = result.Translate({1e4, 1e4, 1e4});
  Manifold farFirst = far.RemoveSelfIntersections();
  ASSERT_EQ(farFirst.Status(), Manifold::Error::NoError);
  ASSERT_LT(InteriorPierces(farFirst), InteriorPierces(far));  // success path
  for (int i = 0; i < 2; ++i) {
    Manifold again = far.RemoveSelfIntersections();
    ExpectMeshGL64GeometryIdentical(again, farFirst);
  }
}

TEST(Manifold, RemoveSelfIntersectionsFarFromOrigin) {
  // Scale-robustness pin: the same genuinely self-intersecting input
  // (hull Subtract, ~31 pierces) translated to 1e4 must still see a
  // strict pierce reduction. At 1e4 the absolute FP grid is ~1e-12,
  // so the arrangement's plane re-projection cancellation is the
  // known precision tax - this asserts reduction, not zero.
  Manifold body = Manifold(ReadTestMeshGL64OBJ("hull-body.obj"));
  Manifold mask = Manifold(ReadTestMeshGL64OBJ("hull-mask.obj"));
  Manifold result = (body - mask).Translate({1e4, 1e4, 1e4});
  ASSERT_EQ(result.Status(), Manifold::Error::NoError);
  ASSERT_GT(InteriorPierces(result), 0);
  Manifold cleaned = result.RemoveSelfIntersections();
  EXPECT_EQ(cleaned.Status(), Manifold::Error::NoError);
  EXPECT_LT(InteriorPierces(cleaned), InteriorPierces(result));
  EXPECT_GT(cleaned.Volume(), 0);
  // All three disjoint hulls survive (the fold class resolves here:
  // the larger scale-derived eps absorbs the degenerate clusters) and
  // the volume holds to 0.1% - the success-path complement of the
  // origin-scale folded-shell fallback.
  EXPECT_EQ(cleaned.Decompose().size(), 3u);
  EXPECT_NEAR(cleaned.Volume(), result.Volume(), result.Volume() * 1e-3);
  // The derived metadata posture: a rebuilt mesh is NOT an original
  // (OriginalID() == -1, matching Manifold(MeshGL64) construction).
  // An emit that called InitializeOriginal() instead would flip this
  // to a fresh nonnegative id - an observable public change.
  EXPECT_EQ(cleaned.OriginalID(), -1);
  // The documented output-tolerance floor: a successful rebuild's
  // tolerance covers at least the 10 x working-eps merge radius (the
  // measured-displacement terms can only widen it further). The bound
  // spells out the pipeline's default eps (= InferEps on this input)
  // through the public bbox, since the member infers it internally.
  EXPECT_GE(cleaned.GetMeshGL64().tolerance,
            10.0 * AlphaBudgetEpsilon(result.BoundingBox().Scale(), 1000));
  // The documented positions-only rebuild contract: a successful
  // rebuild drops non-position properties to numProp == 3.
  Manifold propped = result.SetProperties(
      2, [](double* prop, manifold::vec3 p, const double*) {
        prop[0] = p.x;
        prop[1] = p.y;
      });
  ASSERT_EQ(propped.GetMeshGL64().numProp, 5u);  // 3 position + 2 props
  Manifold cleanedProps = propped.RemoveSelfIntersections();
  EXPECT_EQ(cleanedProps.Status(), Manifold::Error::NoError);
  ASSERT_LT(InteriorPierces(cleanedProps), InteriorPierces(propped));
  const MeshGL64 freshM = cleanedProps.GetMeshGL64();
  EXPECT_EQ(freshM.numProp, 3u);
  // Rebuilt metadata is construction-fresh, never inherited: a new
  // originalID and no tangents (faceID is populated by construction
  // itself, so emptiness is not the contract there; the passthrough/
  // fallback tests pin the opposite - full-field identity).
  EXPECT_NE(freshM.runOriginalID, propped.GetMeshGL64().runOriginalID);
  EXPECT_TRUE(freshM.halfedgeTangent.empty());
}

TEST(Manifold, MeshID) {
  const Manifold cube = Manifold::Cube();
  MeshGL cubeGL = cube.GetMeshGL();
  cubeGL.runIndex.clear();
  cubeGL.runOriginalID.clear();
  Manifold cube1 = Manifold(cubeGL);
  Manifold cube2 = Manifold(cubeGL);
  EXPECT_NE(cube1.GetMeshGL().runOriginalID[0],
            cube2.GetMeshGL().runOriginalID[0]);
}

TEST(Manifold, MeshRelation) {
  Manifold gyroid = WithPositionColors(Gyroid());
  MeshGL gyroidMeshGL = gyroid.GetMeshGL();
  gyroid = gyroid.Simplify();

  if (options.exportModels) WriteTestOBJ("gyroid.obj", gyroid);

  RelatedGL(gyroid, {gyroidMeshGL});
}

TEST(Manifold, MeshRelationTransform) {
  const Manifold cube = Manifold::Cube();
  const MeshGL cubeGL = cube.GetMeshGL();
  const Manifold turned = cube.Rotate(45, 90);

  RelatedGL(turned, {cubeGL});
}

TEST(Manifold, MeshRelationRefine) {
  Manifold csaszar = WithPositionColors(Csaszar()).AsOriginal();
  MeshGL inGL = csaszar.GetMeshGL();

  RelatedGL(csaszar, {inGL});
  csaszar = csaszar.RefineToLength(1);
  ExpectMeshes(csaszar, {{9019, 18038, 3}});
  RelatedGL(csaszar, {inGL});

  if (options.exportModels) WriteTestOBJ("csaszar.obj", csaszar);
}

TEST(Manifold, MeshRelationRefinePrecision) {
  MeshGL inGL = WithPositionColors(Csaszar()).GetMeshGL();
  const int id = inGL.runOriginalID[0];
  Manifold csaszar = Manifold::Smooth(inGL);

  csaszar = csaszar.RefineToTolerance(0.05);
  ExpectMeshes(csaszar, {{2343, 4686, 3}});
  std::vector<uint32_t> runOriginalID = csaszar.GetMeshGL().runOriginalID;
  EXPECT_EQ(runOriginalID.size(), 1);
  EXPECT_EQ(runOriginalID[0], id);

  if (options.exportModels) WriteTestOBJ("csaszarSmooth.obj", csaszar);
}

TEST(Manifold, MeshGLRoundTrip) {
  const Manifold cylinder = Manifold::Cylinder(2, 1);
  EXPECT_GE(cylinder.OriginalID(), 0);
  MeshGL inGL = cylinder.GetMeshGL();
  const Manifold cylinder2(inGL);
  const MeshGL outGL = cylinder2.GetMeshGL();

  EXPECT_EQ(inGL.runOriginalID.size(), 1);
  EXPECT_EQ(outGL.runOriginalID.size(), 1);
  EXPECT_EQ(outGL.runOriginalID[0], inGL.runOriginalID[0]);

  RelatedGL(cylinder2, {inGL});
}

void CheckCube(const MeshGL& cubeSTL) {
  Manifold cube(cubeSTL);
  cube = cube.AsOriginal();
  EXPECT_EQ(cube.NumTri(), 12);
  EXPECT_EQ(cube.NumVert(), 8);
  EXPECT_EQ(cube.NumPropVert(), 24);

  EXPECT_FLOAT_EQ(cube.Volume(), 1.0);
  EXPECT_FLOAT_EQ(cube.SurfaceArea(), 6.0);
}

TEST(Manifold, Merge) {
  MeshGL cubeSTL = CubeSTL();
  EXPECT_EQ(cubeSTL.NumTri(), 12);
  EXPECT_EQ(cubeSTL.NumVert(), 36);

  Manifold cubeBad(cubeSTL);
  EXPECT_TRUE(cubeBad.IsEmpty());
  EXPECT_EQ(cubeBad.Status(), Manifold::Error::NotManifold);

  EXPECT_TRUE(cubeSTL.Merge());
  EXPECT_EQ(cubeSTL.mergeFromVert.size(), 28);
  CheckCube(cubeSTL);

  EXPECT_FALSE(cubeSTL.Merge());
  EXPECT_EQ(cubeSTL.mergeFromVert.size(), 28);
  cubeSTL.mergeFromVert.resize(14);
  cubeSTL.mergeToVert.resize(14);

  EXPECT_TRUE(cubeSTL.Merge());
  EXPECT_EQ(cubeSTL.mergeFromVert.size(), 28);
  CheckCube(cubeSTL);
}

TEST(Manifold, MergeEmpty) {
  MeshGL shape;
  shape.numProp = 7;
  shape.triVerts = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                    24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
  shape.vertProperties = {0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  0.5,  -0.43450000882149, 0.0, 0.0, 0.0, 1.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          0.0,  -0.5, 0.434500008821487, 0.0, 0.0, 1.0, 0.0,
                          0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          -0.0, 0.5,  0.434500008821487, 0.0, 0.0, 0.0, 1.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          0.0,  0.5,  -0.43450000882149, 0.0, 0.0, 1.0, 0.0,
                          0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 1.0, 1.0,
                          0.0,  -0.5, 0.434500008821487, 0.0, 0.0, 0.0, 1.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 1.0, 1.0,
                          0.0,  0.5,  0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          -0.0, 0.5,  0.434500008821487, 0.0, 0.0, 1.0, 0.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 0.0, 0.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 1.0, 1.0,
                          -0.0, 0.5,  0.434500008821487, 0.0, 0.0, 0.0, 1.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 1.0, 1.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 0.0, 0.0,
                          -0.0, -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 0.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  -0.5, 0.434500008821487, 0.0, 0.0, 0.0, 1.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          -0.0, -0.5, 0.434500008821487, 0.0, 0.0, 0.0, 0.0,
                          -0.0, -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 0.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 0.0, 0.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  0.5,  -0.43450000882149, 0.0, 0.0, 0.0, 1.0,
                          -0.0, 0.5,  -0.43450000882149, 0.0, 0.0, 1.0, 1.0,
                          0.0,  -0.5, -0.43450000882149, 0.0, 0.0, 0.0, 0.0,
                          -0.0, -0.5, -0.43450000882149, 0.0, 0.0, 1.0, 0.0};
  EXPECT_TRUE(shape.Merge());
  Manifold man(shape);
  EXPECT_EQ(man.Status(), Manifold::Error::NoError);
  EXPECT_TRUE(man.IsEmpty());
}

TEST(Manifold, PinchedVert) {
  MeshGL shape;
  shape.numProp = 3;
  shape.vertProperties = {0,        0,  0,   //
                          1,        1,  0,   //
                          1,        -1, 0,   //
                          -0.00001, 0,  0,   //
                          -1,       -1, -0,  //
                          -1,       1,  0,   //
                          0,        0,  2,   //
                          0,        0,  -2};
  shape.triVerts = {0, 2, 6,  //
                    2, 1, 6,  //
                    1, 0, 6,  //
                    4, 3, 6,  //
                    3, 5, 6,  //
                    5, 4, 6,  //
                    2, 0, 4,  //
                    0, 3, 4,  //
                    3, 0, 1,  //
                    3, 1, 5,  //
                    7, 2, 4,  //
                    7, 4, 5,  //
                    7, 5, 1,  //
                    7, 1, 2};
  Manifold touch(shape);
  EXPECT_FALSE(touch.IsEmpty());
  EXPECT_EQ(touch.Status(), Manifold::Error::NoError);
  EXPECT_EQ(touch.Genus(), 0);
}

TEST(Manifold, FaceIDRoundTrip) {
  const Manifold cube = Manifold::Cube();
  EXPECT_GE(cube.OriginalID(), 0);
  MeshGL inGL = cube.GetMeshGL();
  EXPECT_EQ(NumUnique(inGL.faceID), 6);
  inGL.faceID = {3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5};

  const Manifold cube2(inGL);
  const MeshGL outGL = cube2.GetMeshGL();
  EXPECT_EQ(NumUnique(outGL.faceID), 2);
}

TEST(Manifold, MirrorUnion) {
  auto a = Manifold::Cube({5., 5., 5.}, true);
  auto b = a.Translate({2.5, 2.5, 2.5});
  auto result = a + b + b.Mirror({1, 1, 0});

  if (options.exportModels) WriteTestOBJ("manifold_mirror_union.obj", result);

  auto vol_a = a.Volume();
  EXPECT_FLOAT_EQ(vol_a * 2.75, result.Volume());
  EXPECT_TRUE(a.Mirror(vec3(0.0)).IsEmpty());
}

TEST(Manifold, MirrorUnion2) {
  auto a = Manifold::Cube();
  auto result = Manifold::BatchBoolean({a.Mirror({1, 0, 0})}, OpType::Add);
  EXPECT_TRUE(result.MatchesTriNormals());
}

TEST(Manifold, Invalid) {
  auto invalid = Manifold::Error::InvalidConstruction;
  auto circ = CrossSection::Circle(10.);
  auto empty_circ = CrossSection::Circle(-2.);
  auto empty_sq = CrossSection::Square(vec2(0.0));

  EXPECT_EQ(Manifold::Sphere(0).Status(), invalid);
  EXPECT_EQ(Manifold::Cylinder(0, 5).Status(), invalid);
  EXPECT_EQ(Manifold::Cylinder(2, -5).Status(), invalid);
  EXPECT_EQ(Manifold::Cylinder(2, 0).Status(), invalid);
  EXPECT_EQ(Manifold::Cylinder(2, 0, 0).Status(), invalid);
  EXPECT_EQ(Manifold::Cube(vec3(0.0)).Status(), invalid);
  EXPECT_EQ(Manifold::Cube({-1, 1, 1}).Status(), invalid);
  EXPECT_EQ(Manifold::Extrude(circ.ToPolygons(), 0.).Status(), invalid);
  EXPECT_EQ(Manifold::Extrude(empty_circ.ToPolygons(), 10.).Status(), invalid);
  EXPECT_EQ(Manifold::Revolve(empty_sq.ToPolygons()).Status(), invalid);
}

TEST(Manifold, MergeDegenerates) {
  MeshGL cube = Manifold::Cube(vec3(1), true).GetMeshGL();
  MeshGL squash;
  squash.vertProperties = cube.vertProperties;
  squash.triVerts = cube.triVerts;
  // Move one vert to the position of its neighbor and remove one triangle
  // linking them to break the manifold.
  squash.vertProperties[squash.vertProperties.size() - 1] *= -1;
  squash.triVerts.resize(squash.triVerts.size() - 3);
  // Rotate the degenerate triangle to the middle to catch more problems.
  std::rotate(squash.triVerts.begin(), squash.triVerts.begin() + 3 * 5,
              squash.triVerts.end());
  // Merge should remove the now duplicate vertex.
  EXPECT_TRUE(squash.Merge());
  // Manifold should remove the triangle with two references to the same vert.
  Manifold squashed = Manifold(squash);
  EXPECT_FALSE(squashed.IsEmpty());
  EXPECT_EQ(squashed.Status(), Manifold::Error::NoError);
}

TEST(Manifold, MergeRefine) {
  MeshGL mesh;
  mesh.tolerance = 1e-5;
  mesh.vertProperties = {
      0.01383194,  -1.88699961, -1.09618223, -0.72694844, -1.30568409,
      -1.09618223, -0.72694844, -1.30568397, -1.09618223, -1.14769018,
      0.00177073,  -1.09618223, -0.72777849, 1.30731535,  -1.09618223,
      -0.72777843, 1.30731547,  -1.09618223, 0.01344040,  1.88713050,
      -1.09618223, 0.75422078,  1.30581498,  -1.09618223, 0.75422555,
      1.30580890,  -1.09618223, 1.14608598,  -0.00000003, -1.09618223,
      0.75343317,  -1.30844986, -1.09618223, 0.00032274,  2.45551205,
      0.70305431,  -0.00045545, 0.49932927,  1.63340366,  -0.00045485,
      0.49933004,  1.63340306,  0.59210259,  2.54137087,  0.46526274,
      -0.00083422, -0.50377685, 1.63340366,  -0.00153223, -2.45697832,
      0.70305431,  0.59018338,  -2.54327321, 0.46526241,  0.59018356,
      -2.54327297, 0.46526247,  -0.00083402, -0.50377715, 1.63340342,
      -0.00153225, -2.45697832, 0.70305425,  -0.00083453, -0.50377727,
      1.63340330,  -0.59329146, -2.54207492, 0.46526280,  -0.00045549,
      0.49932927,  1.63340366,  -0.59137100, 2.54256892,  0.46526301,
      -0.59137243, 2.54256916,  0.46526241,  -1.09916806, 1.80846310,
      0.62867594,  -1.09916782, 1.80846262,  0.62867624,  1.09932983,
      1.80623639,  0.62867647,  1.09933031,  1.80623722,  0.62867594,
      0.59210330,  2.54137087,  0.46526241,  1.09796500,  -1.80969083,
      0.62867600,  -1.10053313, -1.80746412, 0.62867641,  -1.10053349,
      -1.80746484, 0.62867594,  -0.59329236, -2.54207492, 0.46526238,
      -0.00045545, 0.49932924,  1.63340366,  -0.45734894, -0.00169915,
      1.59454787,  -0.00065268, -0.00251516, 1.81522143,  0.45606264,
      -0.00262397, 1.59454787,  0.45606312,  -0.00262326, 1.59454739,
      0.45606264,  -0.00262397, 1.59454787,  -0.00065267, -0.00251517,
      1.81522143,  -0.00083422, -0.50377685, 1.63340366,  1.09796488,
      -1.80969071, 0.62867606,  0.45606279,  -0.00262419, 1.59454775,
      -0.45734891, -0.00169916, 1.59454787,  -0.45734927, -0.00169968,
      1.59454751,  -1.56936228, 1.35052788,  0.52379429,  -0.45734900,
      -0.00169913, 1.59454787,  1.56918728,  1.34735024,  0.52379429,
      1.56811047,  -1.35147583, 0.52380723,  1.56808507,  -1.35131860,
      0.52387440,  -1.57038140, -1.34824181, 0.52379429,  -1.56936240,
      1.35052788,  0.52379423,  -1.44328797, 0.00054518,  0.88277382,
      -1.44327283, 0.00043958,  0.88280576,  1.44195509,  -0.00248179,
      0.88287276,  1.44223559,  -0.00054497, 0.88228679,  1.56918740,
      1.34735012,  0.52379423,  -1.57040441, -1.34818101, 0.52378970,
      -1.57038152, -1.34824169, 0.52379423,  0.00060985,  3.12254715,
      -0.12592883, 0.00060586,  3.12254739,  -0.12593007, -0.21464077,
      3.10757208,  -0.15445986, -0.42353559, 3.05566168,  -0.13569236,
      -0.62559998, 2.98231292,  -0.08878353, 0.62670720,  2.98100901,
      -0.08873826, 0.62661505,  2.98107791,  -0.08880345, 0.42469645,
      3.05480289,  -0.13569044, 0.21584108,  3.10713625,  -0.15445758,
      -0.00174801, -3.12165236, -0.12629783, 0.21322219,  -3.10669637,
      -0.15490441, 0.42281535,  -3.05461216, -0.13599920, 0.62440765,
      -2.98143482, -0.08901373, 0.62447023,  -2.98138809, -0.08896933,
      -0.62786478, -2.98014331, -0.08899876, -0.62786037, -2.98014665,
      -0.08900189, -0.42584255, -3.05390787, -0.13610038, -0.21698712,
      -3.10624146, -0.15494294, -0.00175192, -3.12165260, -0.12629905,
      -0.62560898, 2.98230958,  -0.08879187, -0.78726906, 2.86185288,
      -0.16757384, -0.93646783, 2.72652268,  -0.21237959, -1.07508647,
      2.58084154,  -0.23096853, -1.20526373, 2.42826319,  -0.23021236,
      -1.32975245, 2.27157402,  -0.21704006, -1.45062196, 2.05095077,
      -0.12514453, -1.45065486, 2.05087924,  -0.12510653, 1.45106840,
      2.04791498,  -0.12508209, 1.45102346,  2.04801273,  -0.12513401,
      1.33032060,  2.26888084,  -0.21701908, 1.20595014,  2.42582202,
      -0.23018987, 1.07588828,  2.57866383,  -0.23094597, 0.93737990,
      2.72462535,  -0.21235918, 0.78851759,  2.86004448,  -0.16762893,
      0.78625345,  -2.86083961, -0.16809504, 0.93532181,  -2.72562766,
      -0.21302269, 1.07394040,  -2.57994676, -0.23167852, 1.20411777,
      -2.42736840, -0.23091961, 1.32860637,  -2.27067924, -0.21769993,
      1.44944537,  -2.05011153, -0.12549706, 1.56815219,  -1.35151458,
      0.52370584,  -1.45212758, -2.04719448, -0.12551625, -1.33146667,
      -2.26798582, -0.21772107, -1.20709634, -2.42492723, -0.23094229,
      -1.07703435, -2.57776904, -0.23170128, -0.93852592, -2.72373033,
      -0.21304333, -0.78942937, -2.85936260, -0.16807115, -1.80910003,
      1.01058388,  -0.03365797, -1.84168148, 0.80862963,  -0.10271203,
      -1.86654651, 0.60610521,  -0.15332751, -1.88623869, 0.40378988,
      -0.19160835, -1.90002000, 0.20241861,  -0.21582428, -1.90505385,
      0.00237573,  -0.21922338, -1.90501344, 0.00071962,  -0.21905561,
      -1.90017104, -0.19767725, -0.21568473, -1.88654184, -0.39907673,
      -0.19146930, -1.86700249, -0.60143071, -0.15318951, -1.84229040,
      -0.80400622, -0.10257494, -1.80987811, -1.00592279, -0.03355779,
      1.90388918,  -0.00071938, -0.21917287, 1.89902496,  0.19857220,
      -0.21578738, 1.88539588,  0.39997119,  -0.19157638, 1.86585641,
      0.60232592,  -0.15330335, 1.84114432,  0.80490118,  -0.10269797,
      1.80871546,  1.00692093,  -0.03365800, 1.80794466,  -1.00973761,
      -0.03369782, 1.80795383,  -1.00968957, -0.03371784, 1.84053540,
      -0.80773503, -0.10275956, 1.86540043,  -0.60521013, -0.15336604,
      1.88509274,  -0.40289518, -0.19163992, 1.89881432,  -0.20239551,
      -0.21574667, -1.54677641, 1.84270406,  -0.22365679, -1.63208210,
      1.62775970,  -0.25862685, -1.70706439, 1.40805590,  -0.23586856,
      -1.77522779, 1.18613660,  -0.17659405, -1.80909991, 1.01058459,
      -0.03365869, 1.80871534,  1.00692177,  -0.03365877, 1.77497590,
      1.18254220,  -0.17653592, 1.70699000,  1.40456688,  -0.23577766,
      1.63230419,  1.62443995,  -0.25918254, 1.54714429,  1.83953130,
      -0.22410297, -1.80986142, -1.00602686, -0.03365880, -1.77612185,
      -1.18164730, -0.17739406, -1.70812643, -1.40370321, -0.23700003,
      -1.63330925, -1.62356102, -0.25988522, -1.54816461, -1.83868098,
      -0.22471850, -1.45216942, -2.04711795, -0.12558441, 1.44947588,
      -2.05005574, -0.12554680, 1.54517174,  -1.84302914, -0.22456542,
      1.63157570,  -1.62538552, -0.26028806, 1.70592797,  -1.40712976,
      -0.23690872, 1.77408183,  -1.18524182, -0.17733525, 0.00060586,
      3.12254739,  -0.46959737, -0.21464077, 3.10757208,  -0.45468479,
      0.00060591,  3.12254739,  -0.46959740, 0.21584108,  3.10713625,
      -0.45459676, -0.00175192, -3.12165260, -0.46959740, 0.21322219,
      -3.10669637, -0.45459723, -0.00175193, -3.12165260, -0.46959740,
      -0.21698712, -3.10624146, -0.45467842, -0.42353559, 3.05566168,
      -0.46107247, 0.42469645,  3.05480289,  -0.46099803, 0.42281535,
      -3.05461216, -0.46102279, -0.42584255, -3.05390787, -0.46106711,
      -0.62560898, 2.98230958,  -0.48016420, 0.62440765,  -2.98143482,
      -0.48015898, 0.62661505,  2.98107791,  -0.48015478, -0.62786037,
      -2.98014665, -0.48016420, -0.78726906, 2.86185288,  -0.45703074,
      0.62445766,  -2.98139739, -0.48017445, 0.78625345,  -2.86083961,
      -0.45658660, 0.78851759,  2.86004448,  -0.45661709, 0.62670475,
      2.98101068,  -0.48018256, -0.78942937, -2.85936260, -0.45704728,
      -0.62786043, -2.98014665, -0.48016423, -0.93646783, 2.72652268,
      -0.44622171, 0.93532181,  -2.72562766, -0.44563186, -0.93852592,
      -2.72373033, -0.44624415, 0.93737990,  2.72462535,  -0.44568357,
      -1.07508647, 2.58084154,  -0.44488293, 1.07394040,  -2.57994676,
      -0.44435173, -1.07703435, -2.57776904, -0.44490317, 1.07588828,
      2.57866383,  -0.44439828, -1.20526373, 2.42826319,  -0.45050618,
      1.20411777,  -2.42736840, -0.45019045, -1.20709634, -2.42492723,
      -0.45051819, 1.20595014,  2.42582202,  -0.45021811, -1.32975245,
      2.27157426,  -0.46058667, -1.32975245, 2.27157402,  -0.46058667,
      1.32860637,  -2.27067924, -0.46058667, -1.33146667, -2.26798582,
      -0.46058667, -1.33146656, -2.26798606, -0.46058667, 1.33032060,
      2.26888084,  -0.46058667, 1.33032048,  2.26888084,  -0.46058667,
      -1.45062196, 2.05095077,  -0.45238319, 1.32860649,  -2.27067924,
      -0.46058670, 1.44947588,  -2.05005574, -0.45129201, -1.45216942,
      -2.04711795, -0.45243445, 1.45102346,  2.04801273,  -0.45137730,
      -1.54677641, 1.84270406,  -0.45860177, 1.54517174,  -1.84302914,
      -0.45671692, -1.54816461, -1.83868098, -0.45867091, 1.54714429,
      1.83953130,  -0.45686778, -1.63208210, 1.62775970,  -0.47587365,
      -1.63330925, -1.62356102, -0.47593805, 1.63157570,  -1.62538552,
      -0.47375324, 1.63230419,  1.62443995,  -0.47368667, -1.70706439,
      1.40805590,  -0.50303531, -1.70812643, -1.40370321, -0.50307423,
      1.70592797,  -1.40712976, -0.50065041, 1.70699000,  1.40456688,
      -0.50072700, -1.77522779, 1.18613660,  -0.53636789, -1.77612185,
      -1.18164730, -0.53636789, 1.77308404,  -1.18848991, -0.53344285,
      1.77408183,  -1.18524182, -0.53272378, 1.77497590,  1.18254220,
      -0.53206521, 1.77341759,  1.18763101,  -0.53318721, -1.80909991,
      1.01058459,  -0.51588774, -1.80986142, -1.00602686, -0.51592720,
      -1.77612185, -1.18164718, -0.53636789, 1.80795383,  -1.00968957,
      -0.51267529, 1.80871534,  1.00692177,  -0.51214570, -1.84168148,
      0.80862963,  -0.49827772, -1.84229040, -0.80400622, -0.49834538,
      1.84053540,  -0.80773503, -0.49560285, 1.84114432,  0.80490118,
      -0.49522084, -1.86654651, 0.60610521,  -0.48789516, -1.86700249,
      -0.60143071, -0.48796996, 1.86540043,  -0.60521013, -0.48580948,
      1.86585641,  0.60232592,  -0.48557436, -1.88623869, 0.40378988,
      -0.48232013, -1.88654184, -0.39907673, -0.48238823, 1.88509274,
      -0.40289518, -0.48085603, 1.88539588,  0.39997119,  -0.48076698,
      -1.90002000, 0.20241861,  -0.48219424, -1.90017104, -0.19767725,
      -0.48223990, 1.89881432,  -0.20239551, -0.48138389, 1.89902496,
      0.19857220,  -0.48143867, 1.89892781,  0.20000805,  -0.48143691,
      -1.90505385, 0.00237582,  -0.49013662, -1.90505385, 0.00237573,
      -0.49013659, 1.90388918,  -0.00071975, -0.21917289, 1.90388918,
      -0.00071975, -0.49010351, 1.90388906,  -0.00072469, -0.49010390,
      1.90128529,  0.00611682,  -0.49218670, 1.40738547,  1.02835560,
      -0.88730782, 1.40482545,  -1.03119016, -0.88952076, 1.40470743,
      -1.02958739, -0.88973862, -1.41726339, 1.01772320,  -0.88060957,
      -1.90505385, 0.00237566,  -0.49013662, -1.41691077, -1.01434207,
      -0.88060957, 1.32860637,  -2.27067900, -0.46058673, 1.40398669,
      -1.03246403, -0.88968676, 1.40373409,  1.02831805,  -0.89023006,
      -1.32975245, 2.27157402,  -0.46058670, -1.41726327, 1.01772332,
      -0.88060963, -1.41726327, 1.01772320,  -0.88060963, -1.41691077,
      -1.01434207, -0.88060963, -1.14769018, 0.00177076,  -1.09618223,
      1.14696825,  -0.00000003, -1.09586132, 1.40373087,  1.02831531,
      -0.89023262, 0.40889478,  -2.28727150, -0.89325231, 0.40938482,
      2.28683925,  -0.89325225, -0.62560904, 2.98230958,  -0.48016423,
      -0.39072043, 2.28764248,  -0.89325225, -0.39121491, -2.28645349,
      -0.89325231, 0.40889484,  -2.28727150, -0.89325231, 0.75434840,
      -1.30917358, -1.09577036, 0.40938476,  2.28683901,  -0.89325231,
      -0.39072037, 2.28764224,  -0.89325231, -0.00175192, -3.12165236,
      -0.46959746, 0.40889472,  -2.28727150, -0.89325231, 0.00060591,
      3.12254691,  -0.46959764, -0.39072025, 2.28764224,  -0.89325237,
      0.40938473,  2.28683901,  -0.89325237, 0.01344044,  1.88713050,
      -1.09618223, 0.00032271,  2.45551205,  0.70305431};
  mesh.triVerts = {
      0,   1,   2,   0,   2,   3,   0,   3,   4,   0,   4,   5,   0,   5,   6,
      0,   6,   7,   0,   7,   8,   0,   8,   9,   9,   10,  0,   11,  12,  13,
      13,  14,  11,  15,  16,  17,  15,  17,  18,  18,  19,  15,  20,  15,  21,
      21,  22,  20,  23,  24,  25,  23,  25,  26,  26,  27,  23,  14,  13,  28,
      14,  28,  29,  29,  30,  14,  19,  18,  31,  22,  21,  32,  22,  32,  33,
      33,  34,  22,  35,  27,  36,  36,  37,  35,  12,  37,  38,  12,  38,  39,
      12,  39,  28,  40,  41,  42,  40,  42,  19,  40,  19,  43,  43,  44,  40,
      42,  41,  45,  42,  45,  46,  42,  46,  32,  32,  21,  42,  47,  48,  27,
      27,  26,  47,  28,  39,  49,  49,  29,  28,  43,  31,  50,  43,  50,  51,
      51,  44,  43,  32,  46,  52,  52,  33,  32,  48,  47,  53,  48,  53,  54,
      54,  55,  48,  44,  51,  56,  56,  40,  44,  38,  56,  57,  38,  57,  58,
      38,  58,  49,  49,  39,  38,  45,  55,  59,  45,  59,  60,  45,  60,  52,
      52,  46,  45,  24,  11,  61,  24,  61,  62,  24,  62,  63,  24,  63,  64,
      24,  64,  65,  65,  25,  24,  61,  11,  14,  61,  14,  30,  61,  30,  66,
      61,  66,  67,  61,  67,  68,  68,  69,  61,  17,  20,  70,  17,  70,  71,
      17,  71,  72,  17,  72,  73,  73,  74,  17,  70,  20,  22,  70,  22,  34,
      70,  34,  75,  70,  75,  76,  70,  76,  77,  70,  77,  78,  78,  79,  70,
      25,  65,  80,  25,  80,  81,  25,  81,  82,  25,  82,  83,  25,  83,  84,
      25,  84,  85,  25,  85,  86,  25,  86,  87,  25,  87,  53,  25,  53,  47,
      47,  26,  25,  66,  30,  29,  66,  29,  49,  66,  49,  58,  66,  58,  88,
      66,  88,  89,  66,  89,  90,  66,  90,  91,  66,  91,  92,  66,  92,  93,
      93,  94,  66,  17,  74,  95,  17,  95,  96,  17,  96,  97,  17,  97,  98,
      17,  98,  99,  17,  99,  100, 17,  100, 101, 17,  101, 50,  17,  50,  31,
      52,  60,  102, 52,  102, 103, 52,  103, 104, 52,  104, 105, 52,  105, 106,
      52,  106, 107, 52,  107, 75,  52,  75,  34,  34,  33,  52,  54,  53,  108,
      54,  108, 109, 54,  109, 110, 54,  110, 111, 54,  111, 112, 54,  112, 113,
      113, 114, 54,  54,  114, 115, 54,  115, 116, 54,  116, 117, 54,  117, 118,
      54,  118, 119, 54,  119, 59,  59,  55,  54,  58,  57,  120, 58,  120, 121,
      58,  121, 122, 58,  122, 123, 58,  123, 124, 124, 125, 58,  120, 57,  56,
      120, 56,  51,  120, 51,  50,  120, 50,  101, 120, 101, 126, 120, 126, 127,
      120, 127, 128, 120, 128, 129, 120, 129, 130, 120, 130, 131, 53,  87,  132,
      53,  132, 133, 53,  133, 134, 53,  134, 135, 53,  135, 136, 136, 108, 53,
      58,  125, 137, 58,  137, 138, 58,  138, 139, 58,  139, 140, 58,  140, 141,
      141, 88,  58,  59,  119, 142, 59,  142, 143, 59,  143, 144, 59,  144, 145,
      59,  145, 146, 59,  146, 147, 59,  147, 102, 102, 60,  59,  101, 100, 148,
      101, 148, 149, 101, 149, 150, 101, 150, 151, 101, 151, 152, 152, 126, 101,
      62,  153, 154, 154, 63,  62,  155, 62,  61,  155, 61,  69,  69,  156, 155,
      79,  157, 158, 79,  158, 71,  71,  70,  79,  159, 79,  78,  78,  160, 159,
      63,  154, 161, 161, 64,  63,  69,  68,  162, 162, 156, 69,  71,  158, 163,
      163, 72,  71,  78,  77,  164, 164, 160, 78,  64,  161, 165, 64,  165, 80,
      80,  65,  64,  72,  163, 166, 166, 73,  72,  68,  67,  167, 167, 162, 68,
      77,  76,  168, 168, 164, 77,  80,  165, 169, 169, 81,  80,  73,  166, 170,
      73,  170, 171, 73,  171, 95,  95,  74,  73,  66,  94,  172, 66,  172, 173,
      66,  173, 167, 167, 67,  66,  75,  107, 174, 75,  174, 175, 175, 76,  75,
      81,  169, 176, 176, 82,  81,  95,  171, 177, 177, 96,  95,  107, 106, 178,
      178, 174, 107, 94,  93,  179, 179, 172, 94,  82,  176, 180, 180, 83,  82,
      96,  177, 181, 181, 97,  96,  106, 105, 182, 182, 178, 106, 93,  92,  183,
      183, 179, 93,  83,  180, 184, 184, 84,  83,  97,  181, 185, 185, 98,  97,
      105, 104, 186, 186, 182, 105, 92,  91,  187, 187, 183, 92,  84,  184, 188,
      189, 85,  84,  98,  185, 190, 190, 99,  98,  104, 103, 191, 192, 186, 104,
      91,  90,  193, 194, 187, 91,  85,  189, 195, 195, 86,  85,  99,  196, 197,
      99,  197, 148, 148, 100, 99,  102, 147, 198, 102, 198, 191, 191, 103, 102,
      90,  89,  199, 199, 193, 90,  86,  195, 200, 86,  200, 132, 132, 87,  86,
      148, 197, 201, 201, 149, 148, 147, 146, 202, 202, 198, 147, 88,  141, 203,
      88,  203, 199, 199, 89,  88,  132, 200, 204, 204, 133, 132, 146, 145, 205,
      205, 202, 146, 149, 201, 206, 206, 150, 149, 141, 140, 207, 207, 203, 141,
      133, 204, 208, 208, 134, 133, 145, 144, 209, 209, 205, 145, 150, 206, 210,
      210, 151, 150, 140, 139, 211, 211, 207, 140, 134, 208, 212, 212, 135, 134,
      144, 143, 213, 213, 209, 144, 151, 210, 214, 151, 214, 215, 215, 152, 151,
      139, 138, 216, 139, 216, 217, 217, 211, 139, 135, 212, 218, 218, 136, 135,
      143, 142, 219, 143, 219, 220, 152, 215, 221, 152, 221, 127, 127, 126, 152,
      138, 137, 222, 222, 216, 138, 136, 218, 223, 136, 223, 109, 109, 108, 136,
      119, 118, 224, 119, 224, 219, 219, 142, 119, 127, 221, 225, 225, 128, 127,
      125, 124, 226, 125, 226, 222, 222, 137, 125, 109, 223, 227, 227, 110, 109,
      118, 117, 228, 228, 224, 118, 128, 225, 229, 229, 129, 128, 124, 123, 230,
      230, 226, 124, 110, 227, 231, 231, 111, 110, 117, 116, 232, 232, 228, 117,
      129, 229, 233, 233, 130, 129, 123, 122, 234, 234, 230, 123, 111, 231, 235,
      235, 112, 111, 116, 115, 236, 236, 232, 116, 130, 233, 237, 237, 131, 130,
      122, 121, 238, 122, 238, 239, 239, 234, 122, 112, 235, 240, 240, 113, 112,
      113, 241, 236, 113, 236, 115, 115, 114, 113, 238, 121, 120, 238, 120, 242,
      242, 243, 238, 243, 242, 131, 243, 131, 237, 243, 237, 244, 243, 245, 239,
      239, 238, 243, 239, 245, 246, 239, 246, 217, 239, 217, 216, 239, 216, 222,
      239, 222, 226, 239, 226, 230, 230, 234, 239, 214, 247, 248, 214, 248, 244,
      214, 244, 237, 214, 237, 233, 214, 233, 229, 214, 229, 225, 214, 225, 221,
      221, 215, 214, 212, 249, 240, 212, 240, 235, 212, 235, 231, 212, 231, 227,
      212, 227, 223, 223, 218, 212, 236, 250, 251, 236, 251, 220, 236, 220, 219,
      236, 219, 224, 236, 224, 228, 228, 232, 236, 197, 252, 253, 197, 253, 247,
      197, 247, 214, 197, 214, 210, 197, 210, 206, 206, 201, 197, 217, 246, 254,
      217, 254, 194, 217, 194, 193, 217, 193, 199, 217, 199, 203, 217, 203, 207,
      207, 211, 217, 195, 255, 256, 195, 257, 249, 195, 249, 212, 195, 212, 208,
      195, 208, 204, 204, 200, 195, 220, 258, 192, 220, 191, 198, 220, 198, 202,
      220, 202, 205, 220, 205, 209, 259, 251, 250, 259, 250, 240, 259, 240, 249,
      249, 257, 259, 244, 248, 260, 244, 260, 261, 244, 261, 254, 244, 254, 246,
      244, 246, 245, 171, 170, 262, 171, 262, 252, 171, 190, 185, 171, 185, 181,
      181, 177, 171, 187, 194, 263, 187, 263, 173, 187, 173, 172, 187, 172, 179,
      179, 183, 187, 169, 264, 265, 169, 265, 255, 169, 255, 188, 169, 188, 184,
      169, 184, 180, 180, 176, 169, 192, 266, 175, 192, 175, 174, 192, 174, 178,
      192, 178, 182, 182, 186, 192, 252, 267, 268, 268, 253, 252, 254, 261, 8,
      254, 8,   7,   254, 7,   269, 254, 269, 263, 263, 194, 254, 265, 270, 5,
      265, 4,   256, 256, 255, 265, 192, 258, 2,   192, 2,   1,   1,   266, 192,
      158, 271, 272, 158, 272, 170, 158, 170, 166, 166, 163, 158, 263, 269, 273,
      263, 273, 155, 263, 155, 156, 263, 156, 162, 263, 162, 167, 167, 173, 263,
      271, 160, 164, 271, 164, 168, 168, 266, 271, 153, 273, 274, 153, 274, 270,
      153, 270, 265, 153, 265, 264, 153, 264, 161, 161, 154, 153, 271, 266, 0,
      0,   272, 271, 274, 273, 275, 275, 276, 274, 260, 9,   8,   8,   261, 260,
      9,   260, 248, 9,   248, 247, 9,   247, 253, 9,   253, 268, 268, 10,  9,
      256, 4,   259, 3,   2,   258, 272, 0,   10,  272, 10,  268, 268, 267, 272,
      7,   276, 275, 274, 6,   5,   5,   270, 274, 0,   266, 1,   277, 24,  23};
  mesh.Merge();
  Manifold manifold(mesh);
  manifold = manifold.RefineToLength(1.0);
  EXPECT_NEAR(manifold.Volume(), 31.21, 0.01);
}

#ifdef MANIFOLD_DEBUG
TEST(Manifold, OpenscadCrash) {
  ManifoldParamGuard guard;
  ManifoldParams().processOverlaps = true;
  Manifold m = ReadTestOBJ("openscad-nonmanifold-crash.obj");
  // m is not empty
  EXPECT_EQ(m.IsEmpty(), false);
  Manifold m2 = m + m.Translate({0, 0.6, 0});
  EXPECT_EQ(m2.IsEmpty(), false);
}
#endif

// Deeply-nested CsgOpNode chain (e.g. repeated `+=` in a loop) must not
// stack-overflow in the leaf-counting pre-pass. Cancel up front so we only
// exercise NumLeaves, not the full boolean evaluation.
//
// 300k depth SIGSEGVs on the recursive NumLeaves (8MB default stack on
// macOS/Linux); runs in ~1s / ~1.2GB peak RSS on the iterative version.
TEST(Manifold, DeepChainDoesNotOverflowNumLeaves) {
  constexpr int kDepth = 300000;
  Manifold m = Manifold::Cube(vec3(1), true);
  for (int i = 0; i < kDepth; ++i) {
    m = m + Manifold::Cube(vec3(1), true).Translate(vec3(i * 2.0, 0, 0));
  }
  ExecutionContext ctx;
  ctx.Cancel();
  EXPECT_EQ(m.WithContext(ctx).Status(), Manifold::Error::Cancelled);
  // kDepth + 1 leaves in the chain → kDepth booleans to reduce.
  auto& privateCtx = *ctx.impl_;
  EXPECT_EQ(privateCtx.totalBooleans.load(), kDepth);
}
