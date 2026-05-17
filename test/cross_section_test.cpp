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

#include "manifold/cross_section.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <vector>

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
#include "../src/cross_section/boolean2/boolean2.h"
#include "../src/cross_section/boolean2/intersections.h"
#endif
#include "manifold/common.h"
#include "manifold/manifold.h"
#include "test.h"

using namespace manifold;

namespace {

SimplePolygon RandomTopologicalRing(int n, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<vec2> verts(n);
  for (int i = 0; i < n; ++i) {
    const double theta = 2.0 * kPi * dist(rng);
    verts[i] = {std::cos(theta), std::sin(theta)};
  }

  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::shuffle(order.begin(), order.end(), rng);

  SimplePolygon ring;
  ring.reserve(n);
  for (int i : order) ring.push_back(verts[i]);
  return ring;
}

template <typename Edge>
std::map<int, int> ComputeBalance(const std::vector<Edge>& edges) {
  std::map<int, int> balance;
  for (const auto& edge : edges) {
    balance[edge.v0] += edge.mult;
    balance[edge.v1] -= edge.mult;
  }
  return balance;
}

bool CheckTopologicalValidity(const boolean2::OverlapResult& result,
                              const std::vector<boolean2::EdgeM>& inputEdges,
                              const std::vector<int>& inputRemap,
                              int numMergedVerts) {
  std::vector<boolean2::EdgeM> remapped;
  remapped.reserve(inputEdges.size());
  for (const auto& edge : inputEdges) {
    const int a = inputRemap[edge.v0];
    const int b = inputRemap[edge.v1];
    if (a != b) remapped.push_back({a, b, edge.mult});
  }

  const auto expected = ComputeBalance(remapped);
  const auto actual = ComputeBalance(result.edges);
  for (const auto& [v, actualBalance] : actual) {
    if (v < 0 || v >= static_cast<int>(result.verts.size())) {
      return false;
    }
    (void)actualBalance;
  }
  for (int v = 0; v < static_cast<int>(result.verts.size()); ++v) {
    const int expectedBalance =
        expected.count(v) ? expected.find(v)->second : 0;
    const int actualBalance = actual.count(v) ? actual.find(v)->second : 0;
    const int target = (v < numMergedVerts) ? expectedBalance : 0;
    if (actualBalance != target) return false;
  }
  return true;
}

std::pair<std::vector<vec2>, std::vector<boolean2::EdgeM>> CombinedInput(
    const Polygons& a, const Polygons& b, int bMult) {
  auto [verts, edges] = boolean2::PolygonsToInput(a);
  auto [bVerts, bEdges] = boolean2::PolygonsToInput(b);
  const int base = static_cast<int>(verts.size());
  verts.insert(verts.end(), bVerts.begin(), bVerts.end());
  for (auto edge : bEdges) {
    edge.v0 += base;
    edge.v1 += base;
    edge.mult *= bMult;
    edges.push_back(edge);
  }
  return {std::move(verts), std::move(edges)};
}

}  // namespace

TEST(CrossSection, Square) {
  auto a = Manifold::Cube({5, 5, 5});
  auto b = Manifold::Extrude(CrossSection::Square({5, 5}).ToPolygons(), 5);

  EXPECT_FLOAT_EQ((a - b).Volume(), 0.);
}

TEST(CrossSection, MirrorUnion) {
  auto a = CrossSection::Square({5., 5.}, true);
  auto b = a.Translate({2.5, 2.5});
  auto cross = a + b + b.Mirror({1, 1});
  auto result = Manifold::Extrude(cross.ToPolygons(), 5.);

  if (options.exportModels)
    WriteTestOBJ("cross_section_mirror_union.obj", result);

  EXPECT_FLOAT_EQ(2.5 * a.Area(), cross.Area());
  EXPECT_TRUE(a.Mirror(vec2(0.0)).IsEmpty());
}

TEST(CrossSection, MirrorCheckAxis) {
  auto tri = CrossSection({{0., 0.}, {5., 5.}, {0., 10.}});

  auto a = tri.Mirror({1., 1.}).Bounds();
  auto a_expected = CrossSection({{0., 0.}, {-10., 0.}, {-5., -5.}}).Bounds();

  EXPECT_NEAR(a.min.x, a_expected.min.x, 0.001);
  EXPECT_NEAR(a.min.y, a_expected.min.y, 0.001);
  EXPECT_NEAR(a.max.x, a_expected.max.x, 0.001);
  EXPECT_NEAR(a.max.y, a_expected.max.y, 0.001);

  auto b = tri.Mirror({-1., 1.}).Bounds();
  auto b_expected = CrossSection({{0., 0.}, {10., 0.}, {5., 5.}}).Bounds();

  EXPECT_NEAR(b.min.x, b_expected.min.x, 0.001);
  EXPECT_NEAR(b.min.y, b_expected.min.y, 0.001);
  EXPECT_NEAR(b.max.x, b_expected.max.x, 0.001);
  EXPECT_NEAR(b.max.y, b_expected.max.y, 0.001);
}

TEST(CrossSection, RoundOffset) {
  auto a = CrossSection::Square({20., 20.}, true);
  int segments = 20;
  auto rounded = a.Offset(5., CrossSection::JoinType::Round, 2, segments);
  auto result = Manifold::Extrude(rounded.ToPolygons(), 5.);

  if (options.exportModels)
    WriteTestOBJ("cross_section_round_offset.obj", result);

  EXPECT_EQ(result.Genus(), 0);
  EXPECT_NEAR(result.Volume(), 4386, 1);
  EXPECT_EQ(rounded.NumVert(), segments + 4);
}

TEST(CrossSection, BevelOffset) {
  auto a = CrossSection::Square({20., 20.}, true);
  int segments = 20;
  auto rounded = a.Offset(5., CrossSection::JoinType::Bevel, 2, segments);
  auto result = Manifold::Extrude(rounded.ToPolygons(), 5.);

  if (options.exportModels)
    WriteTestOBJ("cross_section_bevel_offset.obj", result);

  EXPECT_EQ(result.Genus(), 0);
  EXPECT_NEAR(result.Volume(),
              5 * (((20. + (2 * 5.)) * (20. + (2 * 5.))) - (2 * 5. * 5)), 1);
  EXPECT_EQ(rounded.NumVert(), 4 + 4);
}

TEST(CrossSection, MiterOffset) {
  auto square = CrossSection::Square({20., 20.}, true);
  auto offset = square.Offset(5., CrossSection::JoinType::Miter);
  auto result = Manifold::Extrude(offset.ToPolygons(), 1.);

  EXPECT_EQ(result.Genus(), 0);
  EXPECT_NEAR(result.Volume(), 30. * 30., 0.01);
  EXPECT_EQ(offset.NumVert(), 4);
}

TEST(CrossSection, OffsetWithHole) {
  SimplePolygon outer = {{-10, -10}, {10, -10}, {10, 10}, {-10, 10}};
  SimplePolygon hole = {{-2, -2}, {-2, 2}, {2, 2}, {2, -2}};
  CrossSection cs({outer, hole}, CrossSection::FillRule::EvenOdd);

  auto inflated = cs.Offset(1., CrossSection::JoinType::Miter);
  EXPECT_NEAR(inflated.Area(), 484. - 4., 0.01);

  auto inset = cs.Offset(-1., CrossSection::JoinType::Miter);
  EXPECT_NEAR(inset.Area(), 324. - 36., 0.01);
}

TEST(CrossSection, OffsetThinPinch) {
  CrossSection thin = CrossSection::Square({100., 2.}, true);
  auto inset = thin.Offset(-1., CrossSection::JoinType::Miter);

  EXPECT_TRUE(inset.IsEmpty() || inset.Area() < 0.01);
}

TEST(CrossSection, OffsetMultipleDisjointOuters) {
  SimplePolygon a = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  SimplePolygon b = {{10, 0}, {11, 0}, {11, 1}, {10, 1}};
  CrossSection cs({a, b});

  auto separated = cs.Offset(0.5, CrossSection::JoinType::Miter);
  EXPECT_EQ(separated.NumContour(), 2);
  EXPECT_NEAR(separated.Area(), 8., 0.01);

  auto merged = cs.Offset(5., CrossSection::JoinType::Miter);
  EXPECT_EQ(merged.NumContour(), 1);
}

TEST(CrossSection, FourConcurrentEdges) {
  auto rhomb = [](double angleDegrees) {
    const double angle = angleDegrees * kPi / 180.;
    const double cosA = std::cos(angle);
    const double sinA = std::sin(angle);
    constexpr double radius = 1.0;
    constexpr double width = 0.1;
    return SimplePolygon{{radius * cosA, radius * sinA},
                         {-width * sinA, width * cosA},
                         {-radius * cosA, -radius * sinA},
                         {width * sinA, -width * cosA}};
  };

  Polygons polys;
  for (double angle : {0., 45., 90., 135.}) polys.push_back(rhomb(angle));
  CrossSection cs(polys, CrossSection::FillRule::NonZero);

  EXPECT_EQ(cs.NumContour(), 1);
  EXPECT_NEAR(cs.Area(), 0.644423, 1e-4);
  for (const auto& loop : cs.ToPolygons()) {
    for (const vec2& v : loop) {
      EXPECT_TRUE(std::isfinite(v.x));
      EXPECT_TRUE(std::isfinite(v.y));
      EXPECT_LE(la::length(v), 1.05);
    }
  }
}

TEST(CrossSection, ConcurrentIndependentEdgePairs) {
  auto rhomb = [](double angleDegrees) {
    const double angle = angleDegrees * kPi / 180.;
    const double cosA = std::cos(angle);
    const double sinA = std::sin(angle);
    constexpr double radius = 1.0;
    constexpr double width = 0.08;
    return SimplePolygon{{radius * cosA, radius * sinA},
                         {-width * sinA, width * cosA},
                         {-radius * cosA, -radius * sinA},
                         {width * sinA, -width * cosA}};
  };

  Polygons polys;
  for (double angle : {0., 30., 90., 120.}) polys.push_back(rhomb(angle));
  CrossSection cs(polys, CrossSection::FillRule::NonZero);

  EXPECT_EQ(cs.NumContour(), 1);
  EXPECT_NEAR(cs.Area(), 0.527482, 1e-4);
  for (const auto& loop : cs.ToPolygons()) {
    for (const vec2& v : loop) {
      EXPECT_TRUE(std::isfinite(v.x));
      EXPECT_TRUE(std::isfinite(v.y));
      EXPECT_LE(la::length(v), 1.05);
    }
  }
}

TEST(CrossSection, PropagatesShallowIndependentIntersections) {
  using boolean2::Box2;
  using boolean2::BoxOf2DEdge;
  using boolean2::BVH;
  using boolean2::EdgeM;
  using boolean2::FindAndInsertIntersections;

  constexpr double eps = 1e-3;
  std::vector<vec2> verts = {{-1., 0.},
                             {1., 0.},
                             {0., -1.},
                             {0., 1.},
                             {-1., 5. * eps},
                             {1., 5. * eps},
                             {0., -1. + 5. * eps},
                             {0., 1. + 5. * eps},
                             {0., -1.},
                             {0., 1.}};
  std::vector<EdgeM> edges = {
      {0, 1, 1}, {2, 3, 1}, {4, 5, 1}, {6, 7, 1}, {8, 9, 1}};
  std::vector<Box2> edgeBoxes;
  edgeBoxes.reserve(edges.size());
  for (const EdgeM& edge : edges) {
    edgeBoxes.push_back(BoxOf2DEdge(verts[edge.v0], verts[edge.v1], eps));
  }
  std::vector<std::vector<int>> lists(edges.size());
  std::vector<std::vector<int>> vertEdges;
  const std::vector<std::pair<int, int>> pairs = {{0, 1}, {2, 3}};

  FindAndInsertIntersections(edges, &verts, &lists, &vertEdges, eps, edgeBoxes,
                             BVH{}, pairs);

  ASSERT_EQ(verts.size(), 12);
  EXPECT_EQ(lists[4].size(), 2);
}

TEST(CrossSection, TranslatedShallowConcurrentEdges) {
  auto rhomb = [](double angleDegrees, vec2 offset) {
    const double angle = angleDegrees * kPi / 180.;
    const double cosA = std::cos(angle);
    const double sinA = std::sin(angle);
    constexpr double radius = 1.0;
    constexpr double width = 0.08;
    return SimplePolygon{{offset.x + radius * cosA, offset.y + radius * sinA},
                         {offset.x - width * sinA, offset.y + width * cosA},
                         {offset.x - radius * cosA, offset.y - radius * sinA},
                         {offset.x + width * sinA, offset.y - width * cosA}};
  };

  auto polysAt = [&](vec2 offset) {
    Polygons polys;
    for (double angle : {0., 6., 90., 96.})
      polys.push_back(rhomb(angle, offset));
    return polys;
  };

  const double base = std::ldexp(1.0, 40);
  CrossSection origin(polysAt({0., 0.}), CrossSection::FillRule::NonZero);
  CrossSection shifted(polysAt({base, -base}), CrossSection::FillRule::NonZero);
  CrossSection shiftedBack = shifted.Translate({-base, base});

  EXPECT_EQ(origin.NumContour(), 1);
  EXPECT_EQ(shifted.NumContour(), 1);
  EXPECT_EQ(shiftedBack.NumContour(), origin.NumContour());
  EXPECT_NEAR(shiftedBack.Area(), origin.Area(), 1e-4);
  EXPECT_NEAR(shiftedBack.Bounds().Size().x, origin.Bounds().Size().x, 1e-4);
  EXPECT_NEAR(shiftedBack.Bounds().Size().y, origin.Bounds().Size().y, 1e-4);
}

TEST(CrossSection, TranslatedSmallPolygonKeepsFeatures) {
  const double base = std::ldexp(1.0, 49) * 1.5;
  SimplePolygon square = {{base, -base},
                          {base + 10.0, -base},
                          {base + 10.0, -base + 10.0},
                          {base, -base + 10.0}};

  CrossSection cs(square);

  EXPECT_EQ(cs.NumContour(), 1);
  EXPECT_EQ(cs.NumVert(), 4);
  const vec2 size = cs.Bounds().Size();
  EXPECT_NEAR(size.x, 10.0, 1e-9);
  EXPECT_NEAR(size.y, 10.0, 1e-9);
}

// Regression test for the BR-cell hole pattern from Samples.Sponge4. Two
// CCW polygons that share an endpoint AND form a T-junction at the
// non-shared endpoint of one edge. Before the broad-phase / fused-pass
// shared-endpoint filter fix, the (long, short) edge pair was dropped at
// the broad phase, so the narrow phase never inserted the T-junction
// vertex on the long edge. Canonical sub-edges came out with the wrong
// multiplicities, the DCEL face traversal merged faces of different
// windings, and small CW holes were silently dropped from the output.
TEST(CrossSection, TJunctionAtSharedEndpoint) {
  // Outer CCW unit square plus a smaller CCW square sharing the (0,0)
  // corner. The outer's bottom edge (0,0)->(1,0) has the inner's
  // (0.5,0) vertex on its interior; both share endpoint (0,0).
  SimplePolygon outer = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
  SimplePolygon inner = {{0.0, 0.0}, {0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}};

  CrossSection cs(Polygons{outer, inner});
  // Union with Positive fill rule: the inner CCW square is fully inside
  // the outer, so the result is just the outer's area.
  EXPECT_EQ(cs.NumContour(), 1);
  EXPECT_NEAR(cs.Area(), 1.0, 1e-9);

  // A more demanding case from the same family: a butterfly polygon with
  // cancel-pair retraces (one CW lobe around a unit cell) combined with
  // an adjacent CCW staircase polygon whose edge runs along the cell's
  // boundary, creating a T-junction at a vertex shared with the
  // butterfly. This is the minimal three-poly pattern that produced a
  // missing hole in Samples.Sponge4 before the fix.
  SimplePolygon hex = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
  // CW butterfly enclosing the cell (-0.2, -0.2) .. (0.2, 0.2) with a
  // cancel-pair retrace at (-0.5, 0.0)/(0.0, -0.5) to mimic the Sponge16
  // shape. Net signed area is negative (a hole-bearing loop).
  SimplePolygon butterfly = {
      {-0.2, 0.0},  {-0.2, 0.2},  {-0.5, 0.2},  {-0.2, 0.2}, {0.0, 0.2},
      {0.2, 0.2},   {0.2, 0.0},   {0.5, -0.2},  {0.2, -0.2}, {0.2, -0.5},
      {0.2, -0.2},  {0.0, -0.2},  {-0.2, -0.2}};
  // CCW staircase whose v6->v7 edge ((0.2, -0.2) -> (-0.2, -0.2)) is the
  // butterfly's bottom edge, sharing both endpoints with butterfly
  // sub-edges (-0.2, -0.2)..(0.0, -0.2)..(0.2, -0.2) and creating a
  // T-junction at (0.0, -0.2).
  SimplePolygon staircase = {{0.2, -0.2}, {-0.2, -0.2}, {-0.2, -0.4},
                             {0.2, -0.4}};

  CrossSection csBR(Polygons{hex, butterfly, staircase});
  // The butterfly's full interior (cell 0.16 + SE detour triangle 0.03 =
  // 0.19) is a connected CW hole. The CCW staircase below the cell is
  // double-covered (winding 2: hex + staircase), so it stays filled.
  // Expected: outer hex CCW + one CW hole at the butterfly's interior.
  // Without the fix, the butterfly's edges along y=-0.2 (which share
  // endpoints with the staircase's edge) miss the T-junction split at
  // (0.0, -0.2), the canonical mults come out wrong, and the hole
  // collapses or merges with neighboring faces.
  EXPECT_EQ(csBR.NumContour(), 2);
  EXPECT_NEAR(csBR.Area(), 4.0 - 0.19, 1e-9);
}

// Audit follow-up: regression tests for boolean2 filters the post-
// Sponge4 audit argued were correct but didn't have a targeted case.
// Each test verifies a specific topology against a closed-form expected
// result.

// Audit target: canonical sub-edge multiplicity summing for collinear
// overlapping edges. Two CCW rectangles share a collinear segment along
// their boundary but not endpoints; the narrow phase must split A's
// edge at B's endpoints, then Finalize must sum the contributing mults
// correctly so the shared interior segment cancels and only the outer
// "T" outline remains. Without correct mult summing the result would
// drop the shared bottom strip entirely or carry a spurious interior
// edge.
TEST(CrossSection, CollinearSegmentOverlap) {
  SimplePolygon A = {{0.0, 0.0}, {10.0, 0.0}, {10.0, 1.0}, {0.0, 1.0}};
  SimplePolygon B = {{3.0, 0.0}, {7.0, 0.0}, {7.0, 2.0}, {3.0, 2.0}};

  CrossSection cs(Polygons{A, B}, CrossSection::FillRule::NonZero);

  EXPECT_EQ(cs.NumContour(), 1);
  // A = 10, B = 8, overlap (4x1 strip on shared bottom) = 4, union = 14.
  EXPECT_NEAR(cs.Area(), 14.0, 1e-9);
}

// Audit target: collinearity gate in CollectIntersectionPairs. N polygon
// wedges all share the origin as a vertex; every pair shares that
// endpoint. The gate must drop pairs that fan out at the center (no
// T-junction possible) while letting through pairs that need their
// non-shared endpoints checked. Result must be the inscribed regular
// N-gon: one contour, area = N/2 * sin(2pi/N). If the gate over-drops
// (filters something needed), interior wedge boundaries leak through;
// if it under-drops (filters nothing), the test still passes but the
// pre-1a057638 perf gain is gone.
TEST(CrossSection, ManyPolygonsShareCenterVertex) {
  constexpr int N = 8;
  Polygons polys;
  for (int i = 0; i < N; ++i) {
    const double a1 = 2.0 * kPi * i / N;
    const double a2 = 2.0 * kPi * (i + 1) / N;
    polys.push_back({{0.0, 0.0},
                     {std::cos(a1), std::sin(a1)},
                     {std::cos(a2), std::sin(a2)}});
  }
  CrossSection cs(polys, CrossSection::FillRule::NonZero);

  EXPECT_EQ(cs.NumContour(), 1);
  const double expectedArea = 0.5 * N * std::sin(2.0 * kPi / N);
  EXPECT_NEAR(cs.Area(), expectedArea, 1e-9);
}

// Seed: SimplePositiveOffset (2026-05-16 iteration #3)
// Counterexample-hash: 50ede5b9d980d52c
// Suspected owner: pr/boolean2-core (20-gon with extreme-magnitude
//   radii alternating between O(0.1) and O(8.9); Offset(7.21, Bevel)
//   returns a polygon with area ~ input.Area() instead of expanding -
//   effectively a no-op. Likely b2::Offset's normal/miter calc breaks
//   on the near-zero edges produced by the +0.1 floor on tiny radii).
TEST(CrossSection, OffsetPositiveOnExtremeRadiusStar) {
  const std::vector<double> radii = {
      0.,
      0.40098345505108085,
      4.7498113621644447e-99,
      3.7120810186334277,
      8.8389852354367608,
      7.8626648130962875e-111,
      0.37816850000826657,
      0.,
      2.0448856906274785e-158,
      0.,
      0.2582179499017509,
      0.,
      7.224596115948677e-174,
      0.,
      0.21952411055214244,
      0.,
      9.550952284653982e-128,
      0.18006645730017631,
      0.,
      3.9220883255587997e-118};
  SimplePolygon ring;
  ring.reserve(radii.size());
  const int n = static_cast<int>(radii.size());
  for (int i = 0; i < n; ++i) {
    const double r = 0.1 + std::fabs(radii[i]);
    const double theta = 2.0 * kPi * i / n;
    ring.push_back({r * std::cos(theta), r * std::sin(theta)});
  }
  const CrossSection input(ring);
  ASSERT_FALSE(input.IsEmpty());
  ASSERT_GT(std::fabs(input.Area()), 1e-9);

  const double delta = 7.2097955766145416;
  const auto output =
      input.Offset(delta, CrossSection::JoinType::Bevel,
                   /*miter_limit=*/2.0, /*circularSegments=*/0);
  EXPECT_FALSE(output.IsEmpty());
  // A positive offset on a non-self-intersecting ring should always
  // grow the area. Observed locally: output.Area() ~= input.Area(),
  // i.e. Offset is effectively a no-op for this input.
  EXPECT_GE(output.Area(), input.Area() - 1e-6 * (1.0 + input.Area()));
}

// Seed: SubtractInvariants (2026-05-16 iteration #6)
// Counterexample-hash: a68ac82747b27394
// Suspected owner: pr/boolean2-core (extreme magnitude mismatch
//   between two star polygons - A is ~1e-9 scale, B is ~100 scale.
//   At least one of the boolean algebraic invariants
//   `area(A-B)+area(A∩B)=area(A)` etc. fails. Likely an eps-inference
//   or vertex-merge cliff where A's tiny scale gets crushed by B's
//   eps).
TEST(CrossSection, SubtractInvariantsTinyVsLargeStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    ring.reserve(radii.size());
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };

  const std::vector<double> radiiA = {2.0371964671064575e-14,
                                      1.814002251251902e-10, 0.,
                                      2.1825088357767143e-09};
  const std::vector<double> radiiB = {113.5978182908662, 0.,
                                      114.34968677141997,
                                      6.5076626333939721e-10};
  const CrossSection a(star(radiiA));
  const CrossSection b = CrossSection(star(radiiB))
                             .Translate({0., 4.667562921730494e-33});

  const auto aMinusB = a - b;
  const auto bMinusA = b - a;
  const auto aIntersectB = a.Boolean(b, OpType::Intersect);
  const auto aUnionB = a + b;

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aMinusB.Area() + aIntersectB.Area(), a.Area(), tol)
      << "area(A - B) + area(A ∩ B) != area(A)";
  EXPECT_NEAR(bMinusA.Area() + aIntersectB.Area(), b.Area(), tol)
      << "area(B - A) + area(A ∩ B) != area(B)";
  EXPECT_NEAR(aUnionB.Area(), a.Area() + b.Area() - aIntersectB.Area(), tol)
      << "inclusion-exclusion violated";
}

// Seed: SubtractInvariants (2026-05-16 iteration #14)
// Counterexample-hash: a7c02d027b57bf97
// Suspected owner: pr/boolean2-core (two spiky stars at origin - A
//   is 5-pointed with mixed magnitudes including two ~1000-radius
//   spikes and a near-zero vertex; B is 4-pointed with three
//   ~1000-scale spikes plus a near-zero vertex. One of the boolean
//   algebraic invariants fails. Likely a different code path than
//   the iter#6 tiny-vs-large needle case - this is two large spiky
//   shapes, not a tiny shape vs a needle. May share the off-axis
//   T-junction root cause from iter#9 diagnosis, but the spike-
//   collision geometry could be its own failure mode).
TEST(CrossSection, SubtractInvariantsSpikyStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    ring.reserve(radii.size());
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };

  const std::vector<double> radiiA = {1000., 886.10628147264833, 1000., 1., 0.};
  const std::vector<double> radiiB = {1000., 827.10387617193078,
                                      548.20533789242359, 0.};
  const CrossSection a(star(radiiA));
  const CrossSection b(star(radiiB));

  const auto aMinusB = a - b;
  const auto bMinusA = b - a;
  const auto aIntersectB = a.Boolean(b, OpType::Intersect);
  const auto aUnionB = a + b;

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aMinusB.Area() + aIntersectB.Area(), a.Area(), tol)
      << "area(A - B) + area(A ∩ B) != area(A)";
  EXPECT_NEAR(bMinusA.Area() + aIntersectB.Area(), b.Area(), tol)
      << "area(B - A) + area(A ∩ B) != area(B)";
  EXPECT_NEAR(aUnionB.Area(), a.Area() + b.Area() - aIntersectB.Area(), tol)
      << "inclusion-exclusion violated";
}

// Seed: SubtractInvariants (2026-05-17 local sweep)
// Counterexample-hash: a1de309b23c81d31
// Suspected owner: pr/boolean2-core (17-vertex star with most
//   radii zero collapsed to 0.1, plus a 4-pointed star with one
//   ~112 spike). area(b - a) + area(a ∩ b) = 402.04, b.Area() =
//   634.18 - 232 units of area leak. Different geometry from the
//   iter#6/#14 needle/spiky-star seeds; the consumer's `bb4533c2`
//   winding-seed fix didn't cover this case. Reproduces via direct
//   call with these exact args).
TEST(CrossSection, SubtractInvariantsLeakySparseStar) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    ring.reserve(radii.size());
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };

  const std::vector<double> radiiA = {
      627.11994612906153, 0.11978140887764367, 601.04982226530046,
      102.31660501134466, 2.1549912703437601,  0.,
      0.,                 0.,                  0.,
      0.,                 0.,                  0.,
      0.,                 0.,                  0.,
      0.,                 0.};
  const std::vector<double> radiiB = {112.00504648760347,
                                      4.3823848227606392,
                                      1.2795858846899296e-07,
                                      6.7216166802304542};
  const CrossSection a(star(radiiA));
  const CrossSection b = CrossSection(star(radiiB))
                             .Translate({-2.7193614970785894e-29, 0.});

  const auto aMinusB = a - b;
  const auto bMinusA = b - a;
  const auto aIntersectB = a.Boolean(b, OpType::Intersect);
  const auto aUnionB = a + b;

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aMinusB.Area() + aIntersectB.Area(), a.Area(), tol)
      << "area(A - B) + area(A ∩ B) != area(A)";
  EXPECT_NEAR(bMinusA.Area() + aIntersectB.Area(), b.Area(), tol)
      << "area(B - A) + area(A ∩ B) != area(B)";
  EXPECT_NEAR(aUnionB.Area(), a.Area() + b.Area() - aIntersectB.Area(), tol)
      << "inclusion-exclusion violated";
}

// Seed: BooleanCommutativity (2026-05-17 CI fuzz, run 25975407735)
// Counterexample-hash: c4c1f5a6ca197fa8
// Suspected owner: pr/boolean2-core (asymmetric handling of two
//   inputs in the boolean engine - A+B and B+A produce different
//   results. A is an 11-pointed star with mixed magnitudes
//   (~115 down to ~1e-32); B is a 6-pointed star with mixed
//   magnitudes (~96 down to ~1e-53); translation (-0.72, -4.58).
//   The boolean engine treats one input as "subject" and one as
//   "clip"; an order-dependent bug in vertex merge or winding sign
//   would surface here).
TEST(CrossSection, BooleanCommutativityMixedScaleStars) {
  const std::vector<double> radiiA = {
      3.0969192681814191e-32, 2.3071236813515518e-31, 2.353005480384586e-24,
      115.24729490924352,     83.850539440722784,     4.3348536506605848,
      3.0129653594607548,     3.1208956318387746,     4.4747952332988792,
      51.973983183682101,     5.99364665010221e-15};
  const std::vector<double> radiiB = {3.9905161863280855e-53,
                                      0.,
                                      96.331001430191975,
                                      0.,
                                      0.43952001502476951,
                                      48.319452987958854};
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    ring.reserve(radii.size());
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };

  const CrossSection a(star(radiiA));
  const CrossSection b = CrossSection(star(radiiB))
                             .Translate({-0.72134106089064531,
                                         -4.5808240251267858});

  const auto aPlusB = a + b;
  const auto bPlusA = b + a;
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto bIntA = b.Boolean(a, OpType::Intersect);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aPlusB.Area(), bPlusA.Area(), tol) << "A + B != B + A";
  EXPECT_NEAR(aIntB.Area(), bIntA.Area(), tol) << "A ∩ B != B ∩ A";
  EXPECT_EQ(aPlusB.NumContour(), bPlusA.NumContour());
  EXPECT_EQ(aIntB.NumContour(), bIntA.NumContour());
}

// Seed: SubtractInvariants (2026-05-17 CI run 25976076740)
// Counterexample-hash: 4e5ca8dc9060f53e
// Suspected owner: pr/boolean2-core (19-vertex sparse star with
//   only 2 nonzero radii (~39, ~21) vs 4-pointed star with one
//   ~73 spike and tiny radii). area(b - a) + area(a ∩ b) is off
//   by ~234 from b.Area() (~788). Different geometry from prior
//   Subtract seeds.
TEST(CrossSection, SubtractInvariantsTwoNonzeroVsNeedle) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {
      0., 39.150613710409736, 0., 21.472413919643575, 0., 0.,
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
      3.237760857312173e-24};
  const std::vector<double> radiiB = {72.814105930277677,
                                      21.384583693573447,
                                      1.6644635543899127e-243,
                                      1.6688953794582528e-229};
  const CrossSection a(star(radiiA));
  const CrossSection b(star(radiiB));
  const auto aMinusB = a - b;
  const auto bMinusA = b - a;
  const auto aIntersectB = a.Boolean(b, OpType::Intersect);
  const auto aUnionB = a + b;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aMinusB.Area() + aIntersectB.Area(), a.Area(), tol);
  EXPECT_NEAR(bMinusA.Area() + aIntersectB.Area(), b.Area(), tol);
  EXPECT_NEAR(aUnionB.Area(), a.Area() + b.Area() - aIntersectB.Area(), tol);
}

// Seed: PrismBooleanMatchesCrossSection (2026-05-17 CI run 25976076740)
// Counterexample-hash: 74a5b06eb5c583d8
// Suspected owner: pr/boolean2-core (two equilateral triangles with
//   circumradii 0.1 and 0.1+6.88e-13, op=Add. Two near-IDENTICAL
//   triangles - the union should be essentially one triangle.
//   Volume check fails by ~0.26 absolute on h=5; the iter#2
//   narrowing of Prism only skipped Project/Slice for Subtract,
//   not the Volume check which catches this).
TEST(CrossSection, PrismNearIdenticalTrianglesAdd) {
  auto regular = [](int sides, double radius) {
    SimplePolygon ring;
    const double r = 0.1 + std::fabs(radius);
    for (int i = 0; i < sides; ++i) {
      const double th = 2.0 * kPi * i / sides;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const CrossSection a(regular(3, 0.10000000000000001));
  const CrossSection b(regular(3, 0.10000000000068791));
  const auto expected = a + b;
  const double h = 5.0;
  const auto solidA = Manifold::Extrude(a.ToPolygons(), h);
  const auto solidB = Manifold::Extrude(b.ToPolygons(), h);
  const auto result = solidA + solidB;
  const double tolScale = 1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                          std::fabs(expected.Area());
  EXPECT_NEAR(result.Volume(), expected.Area() * h, 1e-5 * tolScale * h);
}

// Seed: BooleanAssociativity (2026-05-17 local sweep on post-743a75b7
// binary)
// Counterexample-hash: 39c5c204f2936291
// Suspected owner: pr/boolean2-core (three 4-vertex stars with
//   mixed magnitudes; (A∪B)∪C = 2.994 but A∪(B∪C) = 5.472 - off
//   by 2.48. Intersection associativity is fine (matches to 3e-17),
//   so the asymmetry is in the union's binary-vs-batch path.
//   Translations are all near-zero (~1e-35 to 1e-28); the trigger
//   is the radii distribution, not translation).
TEST(CrossSection, BooleanAssociativityUnionMixedTriples) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> rA = {24.575460587300253,
                                  2.4572617728922851e-10,
                                  3.3952718785303559e-06, 0.};
  const std::vector<double> rB = {29.731318514644453,
                                  1.5729051003875837e-06,
                                  0.0082858661009423962, 0.};
  const std::vector<double> rC = {0., 0., 6.5474426075871467e-16,
                                  2.9296729904240054e-06};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB))
                             .Translate({-9.1863209962415243e-35,
                                         -9.8444830049208569e-28});
  const CrossSection c = CrossSection(star(rC))
                             .Translate({-6.0489837474564476e-29, 0.});

  const auto ab_c = (a + b) + c;
  const auto a_bc = a + (b + c);
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(ab_c.Area(), a_bc.Area(), tol)
      << "(A ∪ B) ∪ C != A ∪ (B ∪ C)";
}

// Seed: BooleanCommutativity (2026-05-17 CI run 25976076740)
// Counterexample-hash: 570252c8cf569aa4
// Suspected owner: pr/boolean2-core (11-vertex star vs 12-vertex
//   star both with extreme-magnitude radii ranging from O(1e-40)
//   to O(55), translated by (-0.36, 2.83). A+B = 11.12 but B+A =
//   72.75 - off by ~62. Different inputs from
//   DISABLED_BooleanCommutativityMixedScaleStars).
TEST(CrossSection, BooleanCommutativityVeryMixedStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {
      5.5755616189085049e-09, 55.713117722084,
      0.24425911685471227,    12.355613962816753,
      3.2696445771376686e-21, 7.7223346956456043e-27,
      9.1566312909847612e-29, 1.1757300951916251e-32,
      7.5589164688107744e-40, 2.1073720366622818e-34,
      4.7311936525798036e-23};
  const std::vector<double> radiiB = {
      0., 55.263016827250659, 0., 2.1568234275576556e-26,
      4.3751784176513804e-26, 3.0217717245002794e-34, 0., 0., 0.,
      5.2142580450244983e-22, 6.1161686141794186e-15,
      2.8398350775198264};
  const CrossSection a(star(radiiA));
  const CrossSection b = CrossSection(star(radiiB))
                             .Translate({-0.35827067719250716,
                                         2.8298129605573439});
  const auto aPlusB = a + b;
  const auto bPlusA = b + a;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aPlusB.Area(), bPlusA.Area(), tol) << "A + B != B + A";
  EXPECT_EQ(aPlusB.NumContour(), bPlusA.NumContour());
}

// Seed: SubtractInvariants (2026-05-17 CI run 25979228195)
// Counterexample-hash: 340510df5dd776fd
// Suspected owner: pr/boolean2-core (34-vert star vs 20-vert star
//   each with one dominant radius (~33 and ~36) and the rest very
//   small. Inclusion-exclusion violated by ~8 absolute on a∪b≈10.
//   Verified against post-68cbade7 binary).
TEST(CrossSection, SubtractInvariantsDominantSpikeStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> rA = {
      0.15588462038906103, 0.00015361176400222398, 5.8511202761938412e-16,
      3.3009196314697135e-10, 2.7199740700330746e-09, 4.690497141269639e-14,
      8.9922571322790269e-10, 6.3553011716417212e-16, 4.2872009761335588e-18,
      5.6483346599537904e-17, 3.479258338980731e-12, 9.9500479914355979e-14,
      1.822663758661907e-10, 3.7366088533519133e-05, 7.6137290475474603e-12,
      1.5917667211519553e-08, 1.3293803091881353e-06, 6.5804393417497265e-09,
      0.0015777855119048808, 3.5317048943785664e-07, 1.353313563107637e-09,
      1.1788667553671715e-05, 0.00014847370669519226, 1.509263403954591e-11,
      6.5516568176598574e-10, 3.8749326674543751e-07, 3.450117114134831e-06,
      8.3840226251340948e-08, 3.0175612629932874e-07, 33.778786246753299,
      0.033251289695133752, 0., 0., 1.5817601023378236};
  const std::vector<double> rB = {
      3.910135632588194e-07, 1.0064704530951155e-13, 2.0158855952179934e-14,
      1.4165380430901347e-21, 3.6148752792123783e-11, 1.3306453783491879e-10,
      2.3794160990202857e-13, 4.0821728496033454e-14, 4.6816579341479393e-18,
      1.7855060689360453e-17, 2.4858744170474308e-10, 1.5313588823631113e-11,
      6.6960358868004941e-08, 2.1802650917628242e-10, 5.5611296920831909e-13,
      1.9586931478995027e-08, 0., 36.09032519236748, 0.,
      7.5123873966542927e-05};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB))
                             .Translate({0.72140998591309213, 0.});
  const auto aIb = a.Boolean(b, OpType::Intersect);
  const auto aUb = a + b;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aUb.Area(), a.Area() + b.Area() - aIb.Area(), tol)
      << "inclusion-exclusion violated";
}

// Seed: PrismBooleanMatchesCrossSection (2026-05-17 CI run 25979228195)
// Counterexample-hash: 2bda468cc57be581
// Suspected owner: pr/boolean2-core (two equilateral triangles
//   with circumradii 1.675 vs 1.669, op=Add. Volume check fails by
//   ~20 absolute on h=5. Different from the iter#28
//   PrismNearIdenticalTrianglesAdd seed which had radii differing
//   by 6.88e-13; this one has a ~0.4% radius difference. Verified
//   against post-68cbade7 binary.).
TEST(CrossSection, PrismCloseRadiiTrianglesAdd) {
  auto regular = [](int sides, double radius) {
    SimplePolygon ring;
    const double r = 0.1 + std::fabs(radius);
    for (int i = 0; i < sides; ++i) {
      const double th = 2.0 * kPi * i / sides;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const CrossSection a(regular(3, 1.6750962039134867));
  const CrossSection b(regular(3, 1.6691810932888278));
  const auto expected = a + b;
  const double h = 5.0;
  const auto solidA = Manifold::Extrude(a.ToPolygons(), h);
  const auto solidB = Manifold::Extrude(b.ToPolygons(), h);
  const auto result = solidA + solidB;
  const double tolScale = 1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                          std::fabs(expected.Area());
  EXPECT_NEAR(result.Volume(), expected.Area() * h, 1e-5 * tolScale * h);
}

// Seed: BooleanDistributivity (2026-05-17 CI run 25979228195)
// Counterexample-hash: 7532e050f1386752
// Suspected owner: pr/boolean2-core (three 4-vertex stars with all
//   near-zero radii. A∩(B∪C) returns 0 but (A∩B)∪(A∩C) returns
//   0.02 - off by 0.02. Tiny inputs near eps; the distributivity
//   path may have lost a contour during simplification. Verified
//   against post-68cbade7 binary).
TEST(CrossSection, BooleanDistributivityTinyStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> rA = {0., 0., 0., 0.};
  const std::vector<double> rB = {0., 0., 0., 2.6526376418441693e-13};
  const std::vector<double> rC = {5.8033140339376039e-12,
                                  4.8677534136980153e-13, 0., 0.};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB)).Translate({0., 0.});
  const CrossSection c = CrossSection(star(rC))
                             .Translate({-2.4305117516652688e-13, 0.});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIb = a.Boolean(b, OpType::Intersect);
  const auto aIc = a.Boolean(c, OpType::Intersect);
  const auto right = aIb + aIc;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
}

// Seed: BooleanAssociativity (2026-05-17 CI run 25979228195)
// Counterexample-hash: 4b2f604bfb3900e9
// Suspected owner: pr/boolean2-core (three 4-vertex stars with
//   near-zero radii. (A∪B)∪C = 0.04 but A∪(B∪C) = 0.02 - the
//   second form drops one unit of area, likely a contour-merge
//   bug. Distinct from
//   iter#30 BooleanAssociativityUnionMixedTriples seed which had
//   larger radii. Verified against post-68cbade7 binary).
TEST(CrossSection, BooleanAssociativityTinyStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double theta = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(theta), r * std::sin(theta)});
    }
    return ring;
  };
  const std::vector<double> rA = {0., 0., 0., 2.9078938473355343e-13};
  const std::vector<double> rB = {1.4795645772678251e-12,
                                  2.4342935254638573e-13,
                                  1.4800031437954507e-12,
                                  1.1368244372782033e-305};
  const std::vector<double> rC = {0., 0., 0., 0.};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB)).Translate({0., 0.});
  const CrossSection c(star(rC));
  const auto ab_c = (a + b) + c;
  const auto a_bc = a + (b + c);
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(ab_c.Area(), a_bc.Area(), tol)
      << "(A ∪ B) ∪ C != A ∪ (B ∪ C)";
}

TEST(CrossSection, NonFiniteInputReturnsEmpty) {
  const double inf = std::numeric_limits<double>::infinity();
  SimplePolygon bad = {{0.0, 0.0}, {1.0, 0.0}, {inf, 1.0}, {0.0, 1.0}};

  CrossSection constructed(bad);
  EXPECT_TRUE(constructed.IsEmpty());

  Polygons finite = {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
  EXPECT_TRUE(boolean2::Boolean2D(Polygons{bad}, finite, OpType::Add).empty());
}

TEST(CrossSection, SimplifyUsesFixedPointWrapper) {
  Polygons polys{RandomTopologicalRing(8, 618)};
  const double eps = boolean2::InferEps(polys, {});

  const auto [verts, edges] = boolean2::PolygonsToInput(polys);
  auto single = boolean2::RemoveOverlaps2D(verts, edges, eps);
  auto singlePolys = boolean2::OutEdgesToPolygons(single.verts, single.edges);
  ASSERT_EQ(singlePolys.size(), 2);

  auto simplified = boolean2::Simplify(polys, eps);
  ASSERT_EQ(simplified.size(), 2);
  EXPECT_EQ(simplified[0].size(), 11);
  EXPECT_EQ(simplified[1].size(), 3);
  EXPECT_NEAR(boolean2::TotalSignedArea(simplified), 1.7657076120973501,
              1e-12);
}

TEST(CrossSection, ConstructorUsesFixedPointWrapper) {
  Polygons polys{RandomTopologicalRing(8, 618)};
  const double eps = boolean2::InferEps(polys, {});
  const auto [verts, edges] = boolean2::PolygonsToInput(polys);
  auto single = boolean2::RemoveOverlaps2D(verts, edges, eps);
  auto singlePolys = boolean2::OutEdgesToPolygons(single.verts, single.edges);
  ASSERT_EQ(singlePolys.size(), 2);

  CrossSection constructed(polys, CrossSection::FillRule::NonZero);

  EXPECT_EQ(constructed.NumContour(), 1);
  EXPECT_EQ(constructed.NumVert(), 17);
  EXPECT_NEAR(constructed.Area(), 1.7657086786950753, 1e-9);
}

TEST(CrossSection, Empty) {
  Polygons polys(2);
  auto e = CrossSection(polys);
  EXPECT_TRUE(e.IsEmpty());
}

TEST(CrossSection, Rect) {
  double w = 10;
  double h = 5;
  auto rect = Rect({0, 0}, {w, h});
  CrossSection cross(rect);
  auto area = rect.Area();

  EXPECT_FLOAT_EQ(area, w * h);
  EXPECT_FLOAT_EQ(area, cross.Area());
  EXPECT_TRUE(rect.Contains({5, 5}));
  EXPECT_TRUE(rect.Contains(cross.Bounds()));
  EXPECT_TRUE(rect.Contains(Rect()));
  EXPECT_TRUE(rect.DoesOverlap(Rect({5, 5}, {15, 15})));
  EXPECT_TRUE(Rect().IsEmpty());
}

TEST(CrossSection, Transform) {
  auto sq = CrossSection::Square({10., 10.});
  auto a = sq.Rotate(45).Scale({2, 3}).Translate({4, 5});

  mat3 trans({1.0, 0.0, 0.0},  //
             {0.0, 1.0, 0.0},  //
             {4.0, 5.0, 1.0});
  mat3 rot({cosd(45), sind(45), 0.0},   //
           {-sind(45), cosd(45), 0.0},  //
           {0.0, 0.0, 1.0});
  mat3 scale({2.0, 0.0, 0.0},  //
             {0.0, 3.0, 0.0},  //
             {0.0, 0.0, 1.0});

  auto b = sq.Transform(mat2x3(trans * scale * rot));
  auto b_copy = CrossSection(b);

  auto ex_b = Manifold::Extrude(b.ToPolygons(), 1.).GetMeshGL();
  Identical(Manifold::Extrude(a.ToPolygons(), 1.).GetMeshGL(), ex_b);

  // same transformations are applied in b_copy (giving same result)
  Identical(ex_b, Manifold::Extrude(b_copy.ToPolygons(), 1.).GetMeshGL());
}

TEST(CrossSection, Warp) {
  auto sq = CrossSection::Square({10., 10.});
  auto a = sq.Scale({2, 3}).Translate({4, 5});
  auto b = sq.Warp([](vec2& v) {
    v.x = v.x * 2 + 4;
    v.y = v.y * 3 + 5;
  });

  EXPECT_EQ(sq.NumVert(), 4);
  EXPECT_EQ(sq.NumContour(), 1);
}

TEST(CrossSection, Decompose) {
  auto a = CrossSection::Square({2., 2.}, true) -
           CrossSection::Square({1., 1.}, true);
  auto b = a.Translate({4, 4});
  auto ab = a + b;
  auto decomp = ab.Decompose();
  auto recomp = CrossSection::Compose(decomp);

  EXPECT_EQ(decomp.size(), 2);
  EXPECT_EQ(decomp[0].NumContour(), 2);
  EXPECT_EQ(decomp[1].NumContour(), 2);

  Identical(Manifold::Extrude(a.ToPolygons(), 1.).GetMeshGL(),
            Manifold::Extrude(decomp[0].ToPolygons(), 1.).GetMeshGL());
  Identical(Manifold::Extrude(b.ToPolygons(), 1.).GetMeshGL(),
            Manifold::Extrude(decomp[1].ToPolygons(), 1.).GetMeshGL());
  Identical(Manifold::Extrude(ab.ToPolygons(), 1.).GetMeshGL(),
            Manifold::Extrude(recomp.ToPolygons(), 1.).GetMeshGL());
}

TEST(CrossSection, DecomposeNestedHoleAndIsland) {
  SimplePolygon outer = {{-5, -5}, {5, -5}, {5, 5}, {-5, 5}};
  SimplePolygon hole = {{-3, -3}, {-3, 3}, {3, 3}, {3, -3}};
  SimplePolygon island = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  CrossSection input({outer, hole, island}, CrossSection::FillRule::EvenOdd);

  auto components = input.Decompose();
  ASSERT_EQ(components.size(), 2);

  double totalArea = 0.0;
  size_t donutCount = 0;
  size_t islandCount = 0;
  for (const auto& component : components) {
    totalArea += component.Area();
    if (component.NumContour() == 2) {
      ++donutCount;
    } else if (component.NumContour() == 1) {
      ++islandCount;
    }
  }

  EXPECT_EQ(donutCount, 1);
  EXPECT_EQ(islandCount, 1);
  EXPECT_NEAR(totalArea, input.Area(), 1e-9);
  EXPECT_NEAR(CrossSection::Compose(components).Area(), input.Area(), 1e-9);
}

TEST(CrossSection, FillRule) {
  SimplePolygon polygon = {
      {-7, 13},   //
      {-7, 12},   //
      {-5, 9},    //
      {-5, 8.1},  //
      {-4.8, 8},  //
  };

  CrossSection positive(polygon);
  EXPECT_NEAR(positive.Area(), 0.683, 0.001);

  CrossSection negative(polygon, CrossSection::FillRule::Negative);
  EXPECT_NEAR(negative.Area(), 0.193, 0.001);

  CrossSection evenOdd(polygon, CrossSection::FillRule::EvenOdd);
  EXPECT_NEAR(evenOdd.Area(), 0.875, 0.001);

  CrossSection nonZero(polygon, CrossSection::FillRule::NonZero);
  EXPECT_NEAR(nonZero.Area(), 0.875, 0.001);
}

TEST(CrossSection, Hull) {
  auto circ = CrossSection::Circle(10, 360);
  auto circs = std::vector<CrossSection>{circ, circ.Translate({0, 30}),
                                         circ.Translate({30, 0})};
  auto circ_tri = CrossSection::Hull(circs);
  auto centres = SimplePolygon{{0, 0}, {0, 30}, {30, 0}, {15, 5}};
  auto tri = CrossSection::Hull(centres);

  if (options.exportModels) {
    auto circ_tri_ex = Manifold::Extrude(circ_tri.ToPolygons(), 10);
    WriteTestOBJ("cross_section_hull_circ_tri.obj", circ_tri_ex);
  }

  auto circ_area = circ.Area();
  EXPECT_FLOAT_EQ(circ_area, (circ - circ.Scale({0.8, 0.8})).Hull().Area());
  EXPECT_FLOAT_EQ(
      circ_area * 2.5,
      (CrossSection::BatchBoolean(circs, OpType::Add) - tri).Area());
}

TEST(CrossSection, HullError) {
  auto rounded_rectangle = [](double x, double y, double radius, int segments) {
    auto circ = CrossSection::Circle(radius, segments);
    std::vector<CrossSection> vl{};
    vl.push_back(circ.Translate(vec2{radius, radius}));
    vl.push_back(circ.Translate(vec2{x - radius, radius}));
    vl.push_back(circ.Translate(vec2{x - radius, y - radius}));
    vl.push_back(circ.Translate(vec2{radius, y - radius}));
    return CrossSection::Hull(vl);
  };
  auto rr = rounded_rectangle(51, 36, 9.0, 36);

  auto rr_area = rr.Area();
  auto rr_verts = rr.NumVert();
  EXPECT_FLOAT_EQ(rr_area, 1765.1790375559026);
  EXPECT_FLOAT_EQ(rr_verts, 40);
}

TEST(CrossSection, BatchBoolean) {
  CrossSection square = CrossSection::Square({100, 100});
  CrossSection circle1 = CrossSection::Circle(30, 30).Translate({-10, 30});
  CrossSection circle2 = CrossSection::Circle(20, 30).Translate({110, 20});
  CrossSection circle3 = CrossSection::Circle(40, 30).Translate({50, 110});

  CrossSection intersect = CrossSection::BatchBoolean(
      {square, circle1, circle2, circle3}, OpType::Intersect);

  EXPECT_FLOAT_EQ(intersect.Area(), 0);
  EXPECT_FLOAT_EQ(intersect.NumVert(), 0);

  CrossSection add = CrossSection::BatchBoolean(
      {square, circle1, circle2, circle3}, OpType::Add);

  CrossSection subtract = CrossSection::BatchBoolean(
      {square, circle1, circle2, circle3}, OpType::Subtract);

  EXPECT_FLOAT_EQ(add.Area(), 16278.637002);
  EXPECT_FLOAT_EQ(add.NumVert(), 66);

  EXPECT_FLOAT_EQ(subtract.Area(), 7234.478452);
  EXPECT_FLOAT_EQ(subtract.NumVert(), 42);
}

// Seed: BooleanRobustness (2026-05-23 daemon find)
// Counterexample-hash: pending minimization
// Suspected owner: pr/boolean2-core (output topology balance check fails
// after Add on raw multi-contour, near-duplicate-heavy polygons).
TEST(CrossSection, BooleanRobustnessMergeTopologyBalance) {
  const Polygons a = {{{1024., 1024.},
                       {1., 1.},
                       {9.9999999999999995e-07, -9.9999999999999995e-07}}};
  const Polygons b = {{{0., 1.},
                       {9.9999999999999995e-07, 1024.},
                       {-1., -1024.},
                       {-0., 1024.},
                       {1., -1.}},
                      {{1., 1.},
                       {-9.9999999999999995e-07, -0.},
                       {1., 1024.},
                       {9.9999999999999995e-07, 0.}},
                      {{-111.03407576117854, 0.},
                       {560.59976308273758, 2.9313131393310714},
                       {-1024., 1.},
                       {89.678187020143696, 0.},
                       {1024., 1.}},
                      {{0., 0.},
                       {3.2408667776584608, 1024.},
                       {-1000., 9.9999999999999995e-07},
                       {1024., -1.},
                       {-1., -0.}}};
  const auto [verts, edges] = CombinedInput(a, b, /*bMult=*/1);
  ASSERT_FALSE(verts.empty());
  const double eps = boolean2::InferEps(a, b);
  const auto overlap = boolean2::RemoveOverlaps2D(
      verts, edges, eps, /*debug=*/false, boolean2::WindRule::Add);
  EXPECT_TRUE(CheckTopologicalValidity(overlap, edges, overlap.inputRemap,
                                       overlap.numMergedVerts));
}

// Seed: BooleanRobustness (2026-05-23 daemon find)
// Counterexample-hash: pending minimization
// Suspected owner: pr/boolean2-core (direct-cast disconnected-component
// winding fallback emitted an open spur after Add).
TEST(CrossSection, BooleanRobustnessDirectCastKeepsExpectedArea) {
  const Polygons a = {{{-9.9999999999999995e-07, 9.9999999999999995e-07},
                       {-9.9999999999999995e-07, -9.9999999999999995e-07},
                       {-0., 0.},
                       {1., 1024.}},
                      {{1024., 1024.},
                       {1., 1.},
                       {9.9999999999999995e-07, -9.9999999999999995e-07}}};
  const Polygons b = {{{0., 1.},
                       {9.9999999999999995e-07, 1024.},
                       {-1., -1024.},
                       {-0., 1024.},
                       {1., -1.}},
                      {{1., 1.},
                       {-9.9999999999999995e-07, -0.},
                       {1., 1024.},
                       {9.9999999999999995e-07, 0.}},
                      {{-111.03407576117854, 0.},
                       {560.59976308273758, 2.9313131393310714},
                       {-1024., 1.},
                       {89.678187020143696, 0.},
                       {1024., 1.}}};
  const auto [verts, edges] = CombinedInput(a, b, /*bMult=*/1);
  ASSERT_FALSE(verts.empty());
  const double eps = boolean2::InferEps(a, b);
  const auto overlap = boolean2::RemoveOverlaps2D(
      verts, edges, eps, /*debug=*/false, boolean2::WindRule::Add);
  EXPECT_TRUE(CheckTopologicalValidity(overlap, edges, overlap.inputRemap,
                                       overlap.numMergedVerts));

  const auto polys = boolean2::OutEdgesToPolygons(overlap.verts, overlap.edges);
  EXPECT_EQ(polys.size(), 11);
  EXPECT_NEAR(boolean2::TotalSignedArea(polys), 1678.2538553263785,
              1e-10 * (1.0 + 1678.2538553263785));
}

TEST(CrossSection, BooleanOperatorAssignments) {
  CrossSection a = CrossSection::Square({10, 10});
  CrossSection b = CrossSection::Square({10, 10}).Translate({5, 0});

  EXPECT_NEAR((a + b).Area(), 150.0, 1e-9);
  EXPECT_NEAR((a - b).Area(), 50.0, 1e-9);
  EXPECT_NEAR((a ^ b).Area(), 50.0, 1e-9);

  CrossSection add = a;
  add += b;
  EXPECT_NEAR(add.Area(), 150.0, 1e-9);

  CrossSection subtract = a;
  subtract -= b;
  EXPECT_NEAR(subtract.Area(), 50.0, 1e-9);

  CrossSection intersect = a;
  intersect ^= b;
  EXPECT_NEAR(intersect.Area(), 50.0, 1e-9);
}

TEST(CrossSection, NegativeOffset) {
  CrossSection plusSign = CrossSection::Square({30, 50}, true) +
                          CrossSection::Square({50, 30}, true);
  CrossSection dilated =
      plusSign.Offset(-10, CrossSection::JoinType::Round, 2.0, 1024);
  EXPECT_NEAR(dilated.Area(), 30 * 30 - 10 * 10 * kPi, 0.01);
}
