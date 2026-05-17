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

// ---------------------------------------------------------------------
// Build-config guardrails. These exist because we wasted an iteration
// on 2026-05-16 chasing a "boolean2 Simplify idempotence cliff" that
// was actually a Clipper2 cliff: the local fuzz binary had been built
// with `-DMANIFOLD_CROSS_SECTION_BACKEND=clipper2` (the CMake default)
// instead of boolean2. The fuzz targets in this TU assume:
//   - boolean2 backend (cross_section_boolean2.cpp linked),
//   - AddressSanitizer enabled,
//   - UndefinedBehaviorSanitizer enabled.
// If any of these is wrong, the build fails here. The boolean2 check
// is a link-time fingerprint (see cross_section_boolean2.cpp).
// ---------------------------------------------------------------------

#if defined(__has_feature)
#  if !__has_feature(address_sanitizer)
#    error \
        "cross_section_fuzz requires AddressSanitizer (-fsanitize=address). " \
        "Reconfigure with MANIFOLD_FUZZ=ON or rebuild with the fuzz-asan preset."
#  endif
#  if !__has_feature(undefined_behavior_sanitizer)
#    error \
        "cross_section_fuzz requires UndefinedBehaviorSanitizer " \
        "(-fsanitize=undefined). Rebuild with the fuzz-asan preset."
#  endif
#else
#  if !defined(__SANITIZE_ADDRESS__)
#    error \
        "cross_section_fuzz requires AddressSanitizer; no sanitizer feature " \
        "macros detected. Rebuild with the fuzz-asan preset."
#  endif
#endif

// Backend fingerprint: defined only by cross_section_boolean2.cpp.
// A Clipper2 build of the Manifold library has no definition and the
// link below fails with `undefined reference to
// ManifoldCrossSectionBackendIsBoolean2`. A misconfigured-but-linked
// build (e.g. weak-symbol stub) still aborts at process start via the
// constructor below.
extern "C" int ManifoldCrossSectionBackendIsBoolean2();

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "../src/cross_section/boolean2/boolean2.h"
#include "fuzztest/fuzztest.h"
#include "gtest/gtest.h"
#include "manifold/common.h"
#include "manifold/cross_section.h"
#include "manifold/manifold.h"

using namespace fuzztest;

namespace {

[[maybe_unused]] const int kCrossSectionFuzzBackendCheck = [] {
  if (ManifoldCrossSectionBackendIsBoolean2() != 1) {
    std::fprintf(stderr,
                 "FATAL: cross_section_fuzz built against a non-boolean2 "
                 "CrossSection backend (fingerprint != 1). Reconfigure with "
                 "-DMANIFOLD_CROSS_SECTION_BACKEND=boolean2 and rebuild.\n");
    std::abort();
  }
  return 0;
}();

using RawPoint = std::pair<double, double>;
using RawRing = std::vector<RawPoint>;
using RawPolygons = std::vector<RawRing>;
using BooleanSeed = std::tuple<RawPolygons, RawPolygons, manifold::OpType>;
using OffsetSeed =
    std::tuple<RawPolygons, double, manifold::CrossSection::JoinType, int>;
using ExtrudeSeed = std::tuple<RawPolygons, double, int>;
using StarSeed = std::tuple<std::vector<double>>;
using StarExtrudeSeed = std::tuple<std::vector<double>, double, int>;
using BooleanExtrudeSeed =
    std::tuple<std::vector<double>, std::vector<double>, double, double,
               manifold::OpType, double, int>;
using PrismBooleanSeed =
    std::tuple<int, double, int, double, double, double, manifold::OpType>;
using TranslatedSliceSeed =
    std::tuple<int, double, double, double, double, double>;
using RotatedSliceSeed = std::tuple<int, double, double, double>;
using HoledBooleanExtrudeSeed =
    std::tuple<std::vector<double>, std::vector<double>, double, double, double,
               double, manifold::OpType, double, int>;
using TransformExtrudeSeed =
    std::tuple<std::vector<double>, double, double, double, double, int>;
using MultiComponentSeed = std::tuple<std::vector<double>, int, double>;
using BatchBooleanSeed = std::tuple<int, double, int, double>;
using WarpAffineSeed =
    std::tuple<int, double, double, double, double, double, int>;
using MirrorExtrudeSeed = std::tuple<int, double, double, double, double, int>;
using SimpleOffsetSeed = std::tuple<std::vector<double>, double,
                                    manifold::CrossSection::JoinType, int>;
using OffsetExtrudeSeed =
    std::tuple<std::vector<double>, double, manifold::CrossSection::JoinType,
               int, double, int>;
using HoleExtrudeSeed = std::tuple<std::vector<double>, double, double, int>;

manifold::Polygons ToPolygons(const RawPolygons& raw) {
  manifold::Polygons polys;
  polys.reserve(raw.size());
  for (const auto& rawRing : raw) {
    manifold::SimplePolygon ring;
    ring.reserve(rawRing.size());
    for (const auto& [x, y] : rawRing) {
      if (std::isfinite(x) && std::isfinite(y)) ring.push_back({x, y});
    }
    polys.push_back(std::move(ring));
  }
  return polys;
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

bool CheckTopologicalValidity(
    const manifold::boolean2::OverlapResult& result,
    const std::vector<manifold::boolean2::EdgeM>& inputEdges,
    const std::vector<int>& inputRemap, int numMergedVerts) {
  std::vector<manifold::boolean2::EdgeM> remapped;
  remapped.reserve(inputEdges.size());
  for (const auto& edge : inputEdges) {
    const int a = inputRemap[edge.v0];
    const int b = inputRemap[edge.v1];
    if (a != b) remapped.push_back({a, b, edge.mult});
  }

  const auto expected = ComputeBalance(remapped);
  const auto actual = ComputeBalance(result.edges);
  for (int v = 0; v < static_cast<int>(result.verts.size()); ++v) {
    const int expectedBalance =
        expected.count(v) ? expected.find(v)->second : 0;
    const int actualBalance = actual.count(v) ? actual.find(v)->second : 0;
    const int target = (v < numMergedVerts) ? expectedBalance : 0;
    if (actualBalance != target) return false;
  }
  return true;
}

void ExpectFinite(const manifold::Polygons& polys) {
  for (const auto& ring : polys) {
    for (const auto& point : ring) {
      EXPECT_TRUE(std::isfinite(point.x));
      EXPECT_TRUE(std::isfinite(point.y));
    }
  }
}

void ExpectBoolean2TopologyValid(const manifold::Polygons& polys) {
  const auto [verts, edges] = manifold::boolean2::PolygonsToInput(polys);
  if (verts.empty()) return;
  const double eps = manifold::boolean2::InferEps(polys, {});
  const auto result = manifold::boolean2::RemoveOverlaps2D(verts, edges, eps);
  EXPECT_TRUE(CheckTopologicalValidity(result, edges, result.inputRemap,
                                       result.numMergedVerts));
}

void ExpectCrossSectionValid(const manifold::CrossSection& crossSection) {
  EXPECT_TRUE(std::isfinite(crossSection.Area()));
  const auto polys = crossSection.ToPolygons();
  ExpectFinite(polys);
  ExpectBoolean2TopologyValid(polys);

  const auto simplified = crossSection.Simplify();
  EXPECT_TRUE(std::isfinite(simplified.Area()));
  const auto simplifiedAgain = simplified.Simplify();
  EXPECT_NEAR(simplified.Area(), simplifiedAgain.Area(),
              1e-8 * (1.0 + std::fabs(simplified.Area())));
}

manifold::CrossSection ApplyBoolean(const manifold::CrossSection& a,
                                    const manifold::CrossSection& b,
                                    manifold::OpType op) {
  switch (op) {
    case manifold::OpType::Add:
      return a + b;
    case manifold::OpType::Subtract:
      return a - b;
    case manifold::OpType::Intersect:
      return a.Boolean(b, manifold::OpType::Intersect);
  }
  return a + b;
}

manifold::Manifold ApplyBoolean(const manifold::Manifold& a,
                                const manifold::Manifold& b,
                                manifold::OpType op) {
  switch (op) {
    case manifold::OpType::Add:
      return a + b;
    case manifold::OpType::Subtract:
      return a - b;
    case manifold::OpType::Intersect:
      return a ^ b;
  }
  return a + b;
}

// WARNING: this target has been observed to *hang* (infinite loop in
// boolean2 on a specific mutated input). Producer's sweep on
// 2026-05-17 against the post-743a75b7-pre-68cbade7 binary froze
// after corpus_size=256 / fuzz_time=1m22s with CPU pinned but no
// log progress. Could not isolate the exact input (PRNG-driven
// mutation, not in saved corpus); the subsequent fixes (68cbade7
// "Keep winding seed rays inside narrow faces" and later) may have
// resolved it. If you see this target hang again:
//   1. Note the wall time / corpus_size at the freeze.
//   2. `gdb -p <fuzz pid> --batch -ex "thread apply all bt"` to
//      capture the stack.
//   3. Re-run with smaller `--fuzz_for` to narrow the input.
//   4. Once isolated, seed as a regression on pr/boolean2-tests
//      and remove this warning.
void BooleanRobustness(const RawPolygons& rawA, const RawPolygons& rawB,
                       manifold::OpType op) {
  const manifold::Polygons aPolys = ToPolygons(rawA);
  const manifold::Polygons bPolys = ToPolygons(rawB);
  const manifold::CrossSection a(aPolys);
  const manifold::CrossSection b(bPolys);
  const auto result = ApplyBoolean(a, b, op);
  ExpectCrossSectionValid(result);
}

void OffsetRobustness(const RawPolygons& raw, double delta,
                      manifold::CrossSection::JoinType joinType,
                      int circularSegments) {
  const manifold::Polygons polys = ToPolygons(raw);
  const manifold::CrossSection input(polys);
  const auto output =
      input.Offset(delta, joinType, /*miter_limit=*/2.0, circularSegments);
  ExpectCrossSectionValid(output);
}

void ManifoldExtrudeRoundTrip(const RawPolygons& raw, double height,
                              int nDivisions) {
  const manifold::CrossSection input(ToPolygons(raw));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(input.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
}

manifold::SimplePolygon StarPolygon(const std::vector<double>& radii) {
  manifold::SimplePolygon ring;
  ring.reserve(radii.size());
  const double kPi = std::acos(-1.0);
  const int n = static_cast<int>(radii.size());
  for (int i = 0; i < n; ++i) {
    const double r = 0.1 + std::fabs(radii[i]);
    const double theta = 2.0 * kPi * i / n;
    ring.push_back({r * std::cos(theta), r * std::sin(theta)});
  }
  return ring;
}

manifold::SimplePolygon RegularPolygon(int sides, double radius) {
  manifold::SimplePolygon ring;
  ring.reserve(sides);
  const double kPi = std::acos(-1.0);
  const double r = 0.1 + std::fabs(radius);
  for (int i = 0; i < sides; ++i) {
    const double theta = 2.0 * kPi * i / sides;
    ring.push_back({r * std::cos(theta), r * std::sin(theta)});
  }
  return ring;
}

manifold::Polygons StarWithHole(const std::vector<double>& radii,
                                double innerScale) {
  manifold::Polygons polys;
  manifold::SimplePolygon outer;
  outer.reserve(radii.size());
  manifold::SimplePolygon hole;
  hole.reserve(radii.size());
  const double kPi = std::acos(-1.0);
  const int n = static_cast<int>(radii.size());
  for (int i = 0; i < n; ++i) {
    const double theta = 2.0 * kPi * i / n;
    outer.push_back({(10.0 + std::fabs(radii[i])) * std::cos(theta),
                     (10.0 + std::fabs(radii[i])) * std::sin(theta)});
    hole.push_back({(5.0 * innerScale) * std::cos(theta),
                    (5.0 * innerScale) * std::sin(theta)});
  }
  polys.push_back(std::move(outer));
  std::reverse(hole.begin(), hole.end());
  polys.push_back(std::move(hole));
  return polys;
}

std::vector<manifold::CrossSection> SeparatedStars(
    const std::vector<double>& radii, int copies, double spacing) {
  std::vector<manifold::CrossSection> stars;
  stars.reserve(copies);
  const manifold::CrossSection base(StarPolygon(radii));
  const double offset = 2500.0 + spacing;
  for (int i = 0; i < copies; ++i) {
    stars.push_back(base.Translate({offset * i, 0.0}));
  }
  return stars;
}

std::vector<manifold::CrossSection> SeparatedRegulars(int sides, double radius,
                                                      int copies,
                                                      double spacing) {
  std::vector<manifold::CrossSection> sections;
  sections.reserve(copies);
  const double r = 0.1 + std::fabs(radius);
  const manifold::CrossSection base(RegularPolygon(sides, r));
  const double offset = 100.0 + spacing + 3.0 * r;
  for (int i = 0; i < copies; ++i) {
    sections.push_back(base.Translate({offset * i, 0.0}));
  }
  return sections;
}

void ManifoldSimpleExtrudeRoundTrip(const std::vector<double>& radii,
                                    double height, int nDivisions) {
  const manifold::CrossSection input(StarPolygon(radii));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(input.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), input.Area() * height,
              1e-6 * (1.0 + std::fabs(input.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));
}

void BooleanExtrudeRoundTrip(const std::vector<double>& radiiA,
                             const std::vector<double>& radiiB,
                             double translateX, double translateY,
                             manifold::OpType op, double height,
                             int nDivisions) {
  const manifold::CrossSection a(StarPolygon(radiiA));
  const manifold::CrossSection b = manifold::CrossSection(StarPolygon(radiiB))
                                       .Translate({translateX, translateY});
  const auto result = ApplyBoolean(a, b, op);
  ExpectCrossSectionValid(result);
  if (result.IsEmpty() || std::fabs(result.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(result.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), result.Area() * height,
              1e-6 * (1.0 + std::fabs(result.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
}

void PrismBooleanMatchesCrossSection(int sidesA, double radiusA, int sidesB,
                                     double radiusB, double translateX,
                                     double translateY, manifold::OpType op) {
  const manifold::CrossSection a(RegularPolygon(sidesA, radiusA));
  const manifold::CrossSection b =
      manifold::CrossSection(RegularPolygon(sidesB, radiusB))
          .Translate({translateX, translateY});
  const auto expected = ApplyBoolean(a, b, op);
  ExpectCrossSectionValid(expected);

  const double height = 5.0;
  const double areaToleranceScale = 1.0 + std::fabs(a.Area()) +
                                    std::fabs(b.Area()) +
                                    std::fabs(expected.Area());
  const auto solidA = manifold::Manifold::Extrude(a.ToPolygons(), height);
  const auto solidB = manifold::Manifold::Extrude(b.ToPolygons(), height);
  const auto result = ApplyBoolean(solidA, solidB, op);
  EXPECT_EQ(result.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(result.Volume()));
  EXPECT_NEAR(result.Volume(), expected.Area() * height,
              1e-5 * areaToleranceScale * height);

  // Project()/Slice() area checks only run for Add and Intersect.
  // Subtract can produce hole-containing 3D manifolds, and
  // Manifold::Project() returns polygons whose winding is
  // incompatible with the Positive fill rule that CrossSection's
  // default constructor applies - the hole's interior gets counted
  // as filled, yielding the outer hull's area instead of the
  // annulus. That's an upstream Manifold::Project() concern, not a
  // boolean2 concern, so skip the projection assertions for the op
  // that can produce holes.
  if (op != manifold::OpType::Subtract) {
    const manifold::CrossSection projected(result.Project());
    ExpectCrossSectionValid(projected);
    EXPECT_NEAR(projected.Area(), expected.Area(), 1e-5 * areaToleranceScale);

    const manifold::CrossSection middle(result.Slice(height * 0.5));
    ExpectCrossSectionValid(middle);
    EXPECT_NEAR(middle.Area(), expected.Area(), 1e-5 * areaToleranceScale);
  }
}

void DecomposeComposeAndHull(const std::vector<double>& radii, int copies,
                             double spacing) {
  const auto stars = SeparatedStars(radii, copies, spacing);
  const auto input = manifold::CrossSection::Compose(stars);
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto components = input.Decompose();
  EXPECT_FALSE(components.empty());
  double componentArea = 0.0;
  for (const auto& component : components) {
    ExpectCrossSectionValid(component);
    componentArea += component.Area();
  }
  EXPECT_NEAR(componentArea, input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const auto recomposed = manifold::CrossSection::Compose(components);
  ExpectCrossSectionValid(recomposed);
  EXPECT_NEAR(recomposed.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const auto hull = manifold::CrossSection::Hull(components);
  ExpectCrossSectionValid(hull);
  EXPECT_GE(hull.Area(), input.Area() - 1e-6 * (1.0 + input.Area()));
  const auto solid = manifold::Manifold::Extrude(hull.ToPolygons(), 1.0);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
}

void BatchBooleanSeparated(int sides, double radius, int copies,
                           double spacing) {
  const auto sections = SeparatedRegulars(sides, radius, copies, spacing);
  ASSERT_FALSE(sections.empty());

  double sumArea = 0.0;
  for (const auto& section : sections) {
    ExpectCrossSectionValid(section);
    sumArea += section.Area();
  }

  const auto added =
      manifold::CrossSection::BatchBoolean(sections, manifold::OpType::Add);
  ExpectCrossSectionValid(added);
  EXPECT_NEAR(added.Area(), sumArea, 1e-6 * (1.0 + std::fabs(sumArea)));

  const auto intersected = manifold::CrossSection::BatchBoolean(
      sections, manifold::OpType::Intersect);
  ExpectCrossSectionValid(intersected);
  if (sections.size() == 1) {
    EXPECT_NEAR(intersected.Area(), sections.front().Area(),
                1e-6 * (1.0 + std::fabs(sections.front().Area())));
  } else {
    EXPECT_TRUE(intersected.IsEmpty());
  }

  const auto subtracted = manifold::CrossSection::BatchBoolean(
      sections, manifold::OpType::Subtract);
  ExpectCrossSectionValid(subtracted);
  EXPECT_NEAR(subtracted.Area(), sections.front().Area(),
              1e-6 * (1.0 + std::fabs(sections.front().Area())));

  const auto solid = manifold::Manifold::Extrude(added.ToPolygons(), 1.0);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
}

void DecomposedExtrusionsRecompose(int sides, double radius, int copies,
                                   double spacing) {
  const auto sections = SeparatedRegulars(sides, radius, copies, spacing);
  const auto input = manifold::CrossSection::Compose(sections);
  ExpectCrossSectionValid(input);

  const auto components = input.Decompose();
  ASSERT_FALSE(components.empty());

  const double height = 3.0;
  std::vector<manifold::Manifold> solids;
  solids.reserve(components.size());
  for (const auto& component : components) {
    ExpectCrossSectionValid(component);
    solids.push_back(
        manifold::Manifold::Extrude(component.ToPolygons(), height));
    EXPECT_EQ(solids.back().Status(), manifold::Manifold::Error::NoError);
  }

  const auto recomposed =
      manifold::Manifold::BatchBoolean(solids, manifold::OpType::Add);
  EXPECT_EQ(recomposed.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(recomposed.Volume()));
  EXPECT_NEAR(recomposed.Volume(), input.Area() * height,
              1e-6 * (1.0 + std::fabs(input.Area() * height)));

  const manifold::CrossSection projected(recomposed.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const manifold::CrossSection middle(recomposed.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));
}

void TranslatedExtrusionSliceMatchesCrossSection(int sides, double radius,
                                                 double translateX,
                                                 double translateY,
                                                 double translateZ,
                                                 double height) {
  const manifold::CrossSection input(RegularPolygon(sides, radius));
  const auto expected = input.Translate({translateX, translateY});
  ExpectCrossSectionValid(expected);

  const auto solid = manifold::Manifold::Extrude(input.ToPolygons(), height)
                         .Translate({translateX, translateY, translateZ});
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), input.Area() * height,
              1e-6 * (1.0 + std::fabs(input.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), expected.Area(),
              1e-6 * (1.0 + std::fabs(expected.Area())));

  const manifold::CrossSection middle(solid.Slice(translateZ + height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), expected.Area(),
              1e-6 * (1.0 + std::fabs(expected.Area())));
}

void RotatedExtrusionSliceMatchesCrossSection(int sides, double radius,
                                              double rotation, double height) {
  const manifold::CrossSection input(RegularPolygon(sides, radius));
  const auto expected = input.Rotate(rotation);
  ExpectCrossSectionValid(expected);

  const auto solid = manifold::Manifold::Extrude(input.ToPolygons(), height)
                         .Rotate(0.0, 0.0, rotation);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), expected.Area() * height,
              1e-6 * (1.0 + std::fabs(expected.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), expected.Area(),
              1e-6 * (1.0 + std::fabs(expected.Area())));

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), expected.Area(),
              1e-6 * (1.0 + std::fabs(expected.Area())));
}

void HoledBooleanExtrudeRoundTrip(const std::vector<double>& radiiA,
                                  const std::vector<double>& radiiB,
                                  double innerScaleA, double innerScaleB,
                                  double translateX, double translateY,
                                  manifold::OpType op, double height,
                                  int nDivisions) {
  const manifold::CrossSection a(StarWithHole(radiiA, innerScaleA));
  const manifold::CrossSection b =
      manifold::CrossSection(StarWithHole(radiiB, innerScaleB))
          .Translate({translateX, translateY});
  const auto result = ApplyBoolean(a, b, op);
  ExpectCrossSectionValid(result);
  if (result.IsEmpty() || std::fabs(result.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(result.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), result.Area() * height,
              1e-6 * (1.0 + std::fabs(result.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
}

void WarpAffineEquivalence(int sides, double radius, double scaleX,
                           double scaleY, double translateX, double translateY,
                           int nDivisions) {
  const manifold::CrossSection base(RegularPolygon(sides, radius));
  const auto transformed =
      base.Scale({scaleX, scaleY}).Translate({translateX, translateY});
  const auto warped = base.Warp([=](manifold::vec2& v) {
    v.x = v.x * scaleX + translateX;
    v.y = v.y * scaleY + translateY;
  });

  ExpectCrossSectionValid(transformed);
  ExpectCrossSectionValid(warped);
  EXPECT_NEAR(warped.Area(), transformed.Area(),
              1e-6 * (1.0 + std::fabs(transformed.Area())));

  const auto transformedSolid =
      manifold::Manifold::Extrude(transformed.ToPolygons(), 1.0, nDivisions);
  const auto warpedSolid =
      manifold::Manifold::Extrude(warped.ToPolygons(), 1.0, nDivisions);
  EXPECT_EQ(transformedSolid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_EQ(warpedSolid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(transformedSolid.Volume()));
  EXPECT_TRUE(std::isfinite(warpedSolid.Volume()));
  EXPECT_NEAR(warpedSolid.Volume(), transformedSolid.Volume(),
              1e-6 * (1.0 + std::fabs(transformedSolid.Volume())));
}

void MirrorExtrudeRoundTrip(int sides, double radius, double axisX,
                            double axisY, double height, int nDivisions) {
  const manifold::CrossSection base(RegularPolygon(sides, radius));
  const auto mirrored = base.Mirror({axisX, axisY});
  if (axisX == 0.0 && axisY == 0.0) {
    EXPECT_TRUE(mirrored.IsEmpty());
    return;
  }

  ExpectCrossSectionValid(base);
  ExpectCrossSectionValid(mirrored);
  EXPECT_NEAR(mirrored.Area(), base.Area(),
              1e-6 * (1.0 + std::fabs(base.Area())));

  const auto solid =
      manifold::Manifold::Extrude(mirrored.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), mirrored.Area() * height,
              1e-6 * (1.0 + std::fabs(mirrored.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), mirrored.Area(),
              1e-6 * (1.0 + std::fabs(mirrored.Area())));

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), mirrored.Area(),
              1e-6 * (1.0 + std::fabs(mirrored.Area())));
}

void ManifoldTransformedExtrudeRoundTrip(const std::vector<double>& radii,
                                         double rotation, double scaleX,
                                         double scaleY, double height,
                                         int nDivisions) {
  const manifold::CrossSection input =
      manifold::CrossSection(StarPolygon(radii))
          .Rotate(rotation)
          .Scale({scaleX, scaleY});
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(input.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), input.Area() * height,
              1e-6 * (1.0 + std::fabs(input.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));
}

void SimpleBooleanIdentities(const std::vector<double>& radii) {
  const manifold::CrossSection input(StarPolygon(radii));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto unioned = input + input;
  ExpectCrossSectionValid(unioned);
  EXPECT_NEAR(unioned.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const auto intersected = input.Boolean(input, manifold::OpType::Intersect);
  ExpectCrossSectionValid(intersected);
  EXPECT_NEAR(intersected.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const auto subtracted = input - input;
  ExpectCrossSectionValid(subtracted);
  EXPECT_TRUE(subtracted.IsEmpty());
}

// Stresses edge_vert_lists.h:177 apex-skip with `1e-4` collinear
// escape. A triangle whose apex sits at perp distance `apexPerpDist`
// to its base is combined with a thin crossing quad. If the apex-skip
// path drops a legitimate apex/quad intersection on the base edge,
// the union is no longer idempotent under duplicating the triangle.
void ApexSkipNearLine(double apexPerpDist, double crossOffset) {
  if (!std::isfinite(apexPerpDist) || !std::isfinite(crossOffset)) return;
  if (apexPerpDist <= 0) return;

  const manifold::SimplePolygon tri = {
      {-1.0, 0.0}, {1.0, 0.0}, {0.0, apexPerpDist}};
  const manifold::SimplePolygon quad = {{crossOffset - 0.05, -1.0},
                                        {crossOffset + 0.05, -1.0},
                                        {crossOffset + 0.05, 2.0},
                                        {crossOffset - 0.05, 2.0}};

  const manifold::Polygons inputAB{tri, quad};
  const manifold::CrossSection cs(inputAB,
                                  manifold::CrossSection::FillRule::NonZero);
  ExpectCrossSectionValid(cs);

  const manifold::Polygons inputABA{tri, quad, tri};
  const manifold::CrossSection csDup(
      inputABA, manifold::CrossSection::FillRule::NonZero);
  ExpectCrossSectionValid(csDup);
  EXPECT_NEAR(cs.Area(), csDup.Area(),
              1e-6 * (1.0 + std::fabs(cs.Area())));
}

// Large-coordinate translation invariance: a CrossSection
// constructed at the origin and one constructed at large
// translation t should agree on area and contour count. Catches FP-
// precision cliffs in the boolean2 pipeline (eps inference, vertex
// merge, intersection compute) when input coordinates are far from
// the origin.
void TranslationInvariance(const std::vector<double>& radii, double translateX,
                           double translateY) {
  if (!std::isfinite(translateX) || !std::isfinite(translateY)) return;

  const manifold::SimplePolygon ring = StarPolygon(radii);
  const manifold::CrossSection base(ring);
  ExpectCrossSectionValid(base);
  if (base.IsEmpty() || std::fabs(base.Area()) <= 1e-9) return;

  manifold::SimplePolygon translated;
  translated.reserve(ring.size());
  for (const auto& v : ring) {
    translated.push_back({v.x + translateX, v.y + translateY});
  }
  const manifold::CrossSection shifted(translated);
  ExpectCrossSectionValid(shifted);

  // Area and contour count are invariants of translation.
  const double tol = 1e-6 * (1.0 + std::fabs(base.Area()) +
                             std::fabs(translateX) + std::fabs(translateY));
  EXPECT_NEAR(shifted.Area(), base.Area(), tol);
  EXPECT_EQ(shifted.NumContour(), base.NumContour());
}

// Boolean commutativity: A + B == B + A and A ∩ B == B ∩ A. The
// op is mathematically symmetric, but boolean2's internal pipeline
// distinguishes the two inputs (one is the "subject", one is the
// "clip"). Order-dependent bugs in vertex merge, edge collapse, or
// winding sign would show up as area or contour-count mismatches.
void BooleanCommutativity(const std::vector<double>& radiiA,
                          const std::vector<double>& radiiB,
                          double translateX, double translateY) {
  if (!std::isfinite(translateX) || !std::isfinite(translateY)) return;

  const manifold::CrossSection a(StarPolygon(radiiA));
  const manifold::CrossSection b =
      manifold::CrossSection(StarPolygon(radiiB)).Translate({translateX, translateY});
  if (a.IsEmpty() || b.IsEmpty()) return;

  const auto unionAB = a + b;
  const auto unionBA = b + a;
  const auto intersectAB = a.Boolean(b, manifold::OpType::Intersect);
  const auto intersectBA = b.Boolean(a, manifold::OpType::Intersect);
  ExpectCrossSectionValid(unionAB);
  ExpectCrossSectionValid(unionBA);
  ExpectCrossSectionValid(intersectAB);
  ExpectCrossSectionValid(intersectBA);

  const double tol =
      1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(unionAB.Area(), unionBA.Area(), tol) << "A + B != B + A";
  EXPECT_NEAR(intersectAB.Area(), intersectBA.Area(), tol)
      << "A ∩ B != B ∩ A";
  EXPECT_EQ(unionAB.NumContour(), unionBA.NumContour())
      << "A + B and B + A produce different contour counts";
  EXPECT_EQ(intersectAB.NumContour(), intersectBA.NumContour())
      << "A ∩ B and B ∩ A produce different contour counts";
}

// Subtract invariants: for any two 2D regions A and B,
//   area(A - B) + area(A ∩ B) = area(A)
//   area(B - A) + area(A ∩ B) = area(B)
//   area(A ∪ B) = area(A) + area(B) - area(A ∩ B)
// A pure algebra check for the 2D backend's boolean ops. Recovers
// the subtract-side correctness coverage that PrismBooleanMatchesCrossSection
// gave up when it had to skip 3D Project()/Slice() assertions for the
// hole-producing Subtract op.
// Boolean associativity: (A ∪ B) ∪ C == A ∪ (B ∪ C),
// (A ∩ B) ∩ C == A ∩ (B ∩ C). Both ops are associative; reordering
// the bracketing shouldn't change area or contour count. Catches
// order-dependent bugs in the two boolean entry points (binary vs
// batch).
void BooleanAssociativity(const std::vector<double>& radiiA,
                          const std::vector<double>& radiiB,
                          const std::vector<double>& radiiC, double tBx,
                          double tBy, double tCx, double tCy) {
  if (!std::isfinite(tBx) || !std::isfinite(tBy)) return;
  if (!std::isfinite(tCx) || !std::isfinite(tCy)) return;

  const manifold::CrossSection a(StarPolygon(radiiA));
  const manifold::CrossSection b =
      manifold::CrossSection(StarPolygon(radiiB)).Translate({tBx, tBy});
  const manifold::CrossSection c =
      manifold::CrossSection(StarPolygon(radiiC)).Translate({tCx, tCy});
  if (a.IsEmpty() || b.IsEmpty() || c.IsEmpty()) return;

  const auto ab = a + b;
  const auto ab_c = ab + c;
  const auto bc = b + c;
  const auto a_bc = a + bc;
  ExpectCrossSectionValid(ab_c);
  ExpectCrossSectionValid(a_bc);

  const auto aIntB = a.Boolean(b, manifold::OpType::Intersect);
  const auto aIntB_C = aIntB.Boolean(c, manifold::OpType::Intersect);
  const auto bIntC = b.Boolean(c, manifold::OpType::Intersect);
  const auto a_IntBC = a.Boolean(bIntC, manifold::OpType::Intersect);
  ExpectCrossSectionValid(aIntB_C);
  ExpectCrossSectionValid(a_IntBC);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(ab_c.Area(), a_bc.Area(), tol)
      << "(A ∪ B) ∪ C != A ∪ (B ∪ C)";
  EXPECT_NEAR(aIntB_C.Area(), a_IntBC.Area(), tol)
      << "(A ∩ B) ∩ C != A ∩ (B ∩ C)";
  EXPECT_EQ(ab_c.NumContour(), a_bc.NumContour())
      << "union associativity: contour count differs";
}

// Boolean distributivity: A ∩ (B ∪ C) == (A ∩ B) ∪ (A ∩ C).
// Connects union and intersection through a non-trivial three-way
// boolean expression. Order matters here (vs commutativity /
// associativity which only reorder same-op terms) - distributivity
// stresses the engine's handling of interleaved op types.
void BooleanDistributivity(const std::vector<double>& radiiA,
                           const std::vector<double>& radiiB,
                           const std::vector<double>& radiiC, double tBx,
                           double tBy, double tCx, double tCy) {
  if (!std::isfinite(tBx) || !std::isfinite(tBy)) return;
  if (!std::isfinite(tCx) || !std::isfinite(tCy)) return;

  const manifold::CrossSection a(StarPolygon(radiiA));
  const manifold::CrossSection b =
      manifold::CrossSection(StarPolygon(radiiB)).Translate({tBx, tBy});
  const manifold::CrossSection c =
      manifold::CrossSection(StarPolygon(radiiC)).Translate({tCx, tCy});
  if (a.IsEmpty() || b.IsEmpty() || c.IsEmpty()) return;

  // Left side: A ∩ (B ∪ C)
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, manifold::OpType::Intersect);

  // Right side: (A ∩ B) ∪ (A ∩ C)
  const auto aIntB = a.Boolean(b, manifold::OpType::Intersect);
  const auto aIntC = a.Boolean(c, manifold::OpType::Intersect);
  const auto right = aIntB + aIntC;

  ExpectCrossSectionValid(left);
  ExpectCrossSectionValid(right);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_EQ(left.NumContour(), right.NumContour())
      << "distributivity: contour count differs";
}

// Scale invariance: area(f(scale(P, k))) == k^2 * area(f(P)).
// Multiplying all coordinates by a positive scalar k scales area
// quadratically. The boolean engine should produce a scaled result
// for scaled input. Catches FP-precision cliffs at extreme scales
// (eps inference, vertex merge, intersection compute).
void ScaleInvariance(const std::vector<double>& radii, double scale) {
  if (!std::isfinite(scale) || scale <= 0) return;
  if (scale < 1e-4 || scale > 1e4) return;  // bound the domain

  const manifold::SimplePolygon ring = StarPolygon(radii);
  const manifold::CrossSection base(ring);
  ExpectCrossSectionValid(base);
  if (base.IsEmpty() || std::fabs(base.Area()) <= 1e-9) return;

  manifold::SimplePolygon scaled;
  scaled.reserve(ring.size());
  for (const auto& v : ring) scaled.push_back({v.x * scale, v.y * scale});
  const manifold::CrossSection scaledCs(scaled);
  ExpectCrossSectionValid(scaledCs);

  const double expectedArea = base.Area() * scale * scale;
  const double tol = 1e-6 * (1.0 + std::fabs(expectedArea));
  EXPECT_NEAR(scaledCs.Area(), expectedArea, tol)
      << "area(scale(P, " << scale << ")) != " << scale << "^2 * area(P)";
  EXPECT_EQ(scaledCs.NumContour(), base.NumContour())
      << "contour count changed under scaling";
}

// Rotation invariance: area(rotate(P, theta)) == area(P) for any
// angle theta. Rotation preserves area exactly in continuous math;
// FP precision should hold the area to within a tight relative
// tolerance. Catches angle-sensitive choices in eps inference or
// vertex merge that bias one orientation over another.
void RotationInvariance(const std::vector<double>& radii, double thetaRadians) {
  if (!std::isfinite(thetaRadians)) return;

  const manifold::SimplePolygon ring = StarPolygon(radii);
  const manifold::CrossSection base(ring);
  ExpectCrossSectionValid(base);
  if (base.IsEmpty() || std::fabs(base.Area()) <= 1e-9) return;

  const double c = std::cos(thetaRadians);
  const double s = std::sin(thetaRadians);
  manifold::SimplePolygon rotated;
  rotated.reserve(ring.size());
  for (const auto& v : ring) {
    rotated.push_back({c * v.x - s * v.y, s * v.x + c * v.y});
  }
  const manifold::CrossSection rotatedCs(rotated);
  ExpectCrossSectionValid(rotatedCs);

  const double tol = 1e-6 * (1.0 + std::fabs(base.Area()));
  EXPECT_NEAR(rotatedCs.Area(), base.Area(), tol)
      << "area changed under rotation by " << thetaRadians << " radians";
  EXPECT_EQ(rotatedCs.NumContour(), base.NumContour())
      << "contour count changed under rotation";
}

// Offset(0) is identity: input.Offset(0.0, anyJoin) == input.
// A trivial-delta sanity check. Offset's polygon-vert-walk should
// produce identical output when no expansion is requested. Catches
// edge-case bugs in the offset entry where delta=0 takes a wrong
// branch (e.g. dividing by delta, or skipping the join-corner
// generation entirely).
void OffsetIdentityAtZero(const std::vector<double>& radii,
                          manifold::CrossSection::JoinType joinType) {
  const manifold::CrossSection input(StarPolygon(radii));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto output =
      input.Offset(0.0, joinType, /*miter_limit=*/2.0, /*segments=*/0);
  ExpectCrossSectionValid(output);

  const double tol = 1e-6 * (1.0 + std::fabs(input.Area()));
  EXPECT_NEAR(output.Area(), input.Area(), tol)
      << "Offset(0, " << static_cast<int>(joinType) << ") changed area";
  EXPECT_EQ(output.NumContour(), input.NumContour())
      << "Offset(0) changed contour count";
}

// Empty identities: A + empty == A, A ∩ empty == empty, A - empty
// == A. Exercises the boolean engine's handling of empty inputs.
// Catches null-pointer / empty-vector edge cases that don't surface
// in the standard fuzz targets (which generate non-trivial inputs).
void EmptyIdentities(const std::vector<double>& radii) {
  const manifold::CrossSection a(StarPolygon(radii));
  ExpectCrossSectionValid(a);
  if (a.IsEmpty() || std::fabs(a.Area()) <= 1e-9) return;
  const manifold::CrossSection e;
  ASSERT_TRUE(e.IsEmpty());

  const auto aUnionE = a + e;
  const auto eUnionA = e + a;
  const auto aIntE = a.Boolean(e, manifold::OpType::Intersect);
  const auto aMinusE = a - e;
  ExpectCrossSectionValid(aUnionE);
  ExpectCrossSectionValid(eUnionA);
  ExpectCrossSectionValid(aIntE);
  ExpectCrossSectionValid(aMinusE);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()));
  EXPECT_NEAR(aUnionE.Area(), a.Area(), tol) << "A + empty != A";
  EXPECT_NEAR(eUnionA.Area(), a.Area(), tol) << "empty + A != A";
  EXPECT_TRUE(aIntE.IsEmpty()) << "A ∩ empty is non-empty";
  EXPECT_NEAR(aMinusE.Area(), a.Area(), tol) << "A - empty != A";
}

// Double mirror: A.Mirror(axis).Mirror(axis) == A. Mirror is its
// own inverse along any axis. Catches sign-flip or axis-handling
// bugs in the transform path.
void DoubleMirrorIdentity(const std::vector<double>& radii, double axisX,
                          double axisY) {
  if (!std::isfinite(axisX) || !std::isfinite(axisY)) return;
  if (std::fabs(axisX) + std::fabs(axisY) < 1e-9) return;  // zero axis

  const manifold::CrossSection a(StarPolygon(radii));
  ExpectCrossSectionValid(a);
  if (a.IsEmpty() || std::fabs(a.Area()) <= 1e-9) return;

  const manifold::vec2 axis{axisX, axisY};
  const auto twice = a.Mirror(axis).Mirror(axis);
  ExpectCrossSectionValid(twice);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()));
  EXPECT_NEAR(twice.Area(), a.Area(), tol)
      << "Mirror(axis).Mirror(axis) changed area";
  EXPECT_EQ(twice.NumContour(), a.NumContour())
      << "Mirror.Mirror changed contour count";
}

void SubtractInvariants(const std::vector<double>& radiiA,
                        const std::vector<double>& radiiB, double translateX,
                        double translateY) {
  if (!std::isfinite(translateX) || !std::isfinite(translateY)) return;

  const manifold::CrossSection a(StarPolygon(radiiA));
  const manifold::CrossSection b =
      manifold::CrossSection(StarPolygon(radiiB)).Translate({translateX, translateY});
  ExpectCrossSectionValid(a);
  ExpectCrossSectionValid(b);
  if (a.IsEmpty() || b.IsEmpty()) return;

  const auto aMinusB = a - b;
  const auto bMinusA = b - a;
  const auto aIntersectB = a.Boolean(b, manifold::OpType::Intersect);
  const auto aUnionB = a + b;
  ExpectCrossSectionValid(aMinusB);
  ExpectCrossSectionValid(bMinusA);
  ExpectCrossSectionValid(aIntersectB);
  ExpectCrossSectionValid(aUnionB);

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));

  // area(A - B) + area(A ∩ B) == area(A)
  EXPECT_NEAR(aMinusB.Area() + aIntersectB.Area(), a.Area(), tol)
      << "area(A - B) + area(A ∩ B) != area(A)";

  // area(B - A) + area(A ∩ B) == area(B)
  EXPECT_NEAR(bMinusA.Area() + aIntersectB.Area(), b.Area(), tol)
      << "area(B - A) + area(A ∩ B) != area(B)";

  // area(A ∪ B) == area(A) + area(B) - area(A ∩ B)
  EXPECT_NEAR(aUnionB.Area(), a.Area() + b.Area() - aIntersectB.Area(), tol)
      << "inclusion-exclusion violated";
}

// Iterate-to-fixed-point convergence: a random star polygon should
// converge under iterate.h's fingerprint loop within a few iterations.
// Smith's alpha-budget proves <=2 iterations with symbolic
// intersections; production passes maxIter=2 and ignores the status.
// This dimension exercises maxIter=10 and asserts the result is
// either Converged or Cycled, never MaxedOut. A MaxedOut counterexample
// is itself a regression seed - it means topology drifts past
// iteration 10, which today would only show up as off-by-eps cliffs.
void IterateToFixedPointConverges(const std::vector<double>& radii) {
  const manifold::Polygons polys{StarPolygon(radii)};
  const auto [verts, edges] = manifold::boolean2::PolygonsToInput(polys);
  if (verts.empty()) return;
  const double eps = manifold::boolean2::InferEps(polys, {});

  int iters = -1;
  manifold::boolean2::IterStatus status = manifold::boolean2::IterStatus::MaxedOut;
  manifold::boolean2::IterateToFixedPoint(verts, edges, eps, /*maxIter=*/10,
                                          &iters, &status);

  EXPECT_NE(static_cast<int>(status),
            static_cast<int>(manifold::boolean2::IterStatus::MaxedOut))
      << "IterateToFixedPoint hit MaxedOut at iter=10 (production maxIter=2 "
      << "would have produced different topology). Input is a regression "
      << "seed for boolean2 iteration tail behavior.";
}

// Decompose/Compose round-trip on a HOLED CrossSection. The existing
// DecomposeComposeAndHull covers separated stars (no negative-orientation
// rings), which doesn't exercise hole containment in the decompose
// path. Build a holed shape via outer - inner_translated_subtract, then
// Decompose -> Compose and assert area + NumContour preservation.
void DecomposeRecomposeWithHoles(const std::vector<double>& outerRadii,
                                 const std::vector<double>& holeRadii,
                                 double holeOffsetX, double holeOffsetY) {
  if (!std::isfinite(holeOffsetX) || !std::isfinite(holeOffsetY)) return;

  const manifold::CrossSection outer(StarPolygon(outerRadii));
  const manifold::CrossSection hole =
      manifold::CrossSection(StarPolygon(holeRadii))
          .Translate({holeOffsetX, holeOffsetY});
  ExpectCrossSectionValid(outer);
  ExpectCrossSectionValid(hole);
  if (outer.IsEmpty() || hole.IsEmpty()) return;
  if (std::fabs(outer.Area()) <= 1e-9) return;

  const auto holed = outer - hole;
  ExpectCrossSectionValid(holed);
  if (holed.IsEmpty()) return;
  // Need actual multi-ring topology for this dim. If the hole missed
  // or completely swallowed the outer, skip.
  if (holed.NumContour() < 2) return;

  const auto components = holed.Decompose();
  EXPECT_FALSE(components.empty());

  double componentArea = 0.0;
  size_t componentContourSum = 0;
  for (const auto& component : components) {
    ExpectCrossSectionValid(component);
    componentArea += component.Area();
    componentContourSum += component.NumContour();
  }
  const double tol = 1e-6 * (1.0 + std::fabs(holed.Area()));
  EXPECT_NEAR(componentArea, holed.Area(), tol)
      << "Decompose lost area on holed input";
  EXPECT_EQ(componentContourSum, holed.NumContour())
      << "Decompose split or merged contours unexpectedly";

  const auto recomposed = manifold::CrossSection::Compose(components);
  ExpectCrossSectionValid(recomposed);
  EXPECT_NEAR(recomposed.Area(), holed.Area(), tol)
      << "Compose(Decompose(holed)) changed area";
  EXPECT_EQ(recomposed.NumContour(), holed.NumContour())
      << "Compose(Decompose(holed)) changed contour count";
}

// Offset round-trip on convex inputs: input.Offset(d, Miter).Offset(-d,
// Miter) should return to the original area (and same contour count) for
// delta small enough that the negative offset doesn't collapse features.
// Miter-only because Square clips sharp corners flat (legitimately losing
// area in the round-trip for triangles / sharp n-gons) and Round
// introduces curvature-dependent area drift; both are properties of the
// join type, not regressions. With miter_limit=2.0 and a regular n-gon
// (exterior angle 2*pi/n), the miter never clips because
// delta/sin(pi/n - pi/2) = delta/cos(pi/n) <= 2*delta when cos(pi/n) >=
// 0.5, i.e. n >= 3. Bounds delta to a fraction of the inscribed radius
// to keep the un-offset polygon from vanishing.

// Hull idempotence: Hull(Hull(X)) == Hull(X). The convex hull of any
// already-convex shape is itself, so applying Hull twice should be
// indistinguishable from once. Catches vertex jitter, contour drift,
// or area regression on the second hull pass.
void HullIdempotence(const std::vector<double>& radii) {
  const manifold::CrossSection input(StarPolygon(radii));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto hull1 = input.Hull();
  ExpectCrossSectionValid(hull1);
  if (hull1.IsEmpty() || std::fabs(hull1.Area()) <= 1e-9) return;

  const auto hull2 = hull1.Hull();
  ExpectCrossSectionValid(hull2);

  const double tol = 1e-6 * (1.0 + std::fabs(hull1.Area()));
  EXPECT_NEAR(hull2.Area(), hull1.Area(), tol)
      << "Hull(Hull(X)).Area() != Hull(X).Area()";
  EXPECT_EQ(hull2.NumContour(), hull1.NumContour())
      << "Hull(Hull(X)).NumContour() != Hull(X).NumContour()";
  EXPECT_EQ(hull2.NumVert(), hull1.NumVert())
      << "Hull(Hull(X)).NumVert() != Hull(X).NumVert()";
}

void OffsetInverseConvex(int sides, double radius, double delta) {
  if (!std::isfinite(delta)) return;
  // n=3 (equilateral triangle) drifts ~0.35% on the boolean2 Offset
  // round-trip even with Miter and well under miter_limit; captured as
  // a separate regression seed on pr/boolean2-tests. Restrict the
  // fuzz dim to n>=4 so it catches future regressions in less-sharp
  // n-gon round-trips, not just rediscover the triangle case.
  if (sides < 4 || sides > 32) return;

  const double kPi = std::acos(-1.0);
  const double r = 0.1 + std::fabs(radius);
  // Inscribed radius for a regular n-gon = r * cos(pi/n). Keep |delta|
  // safely below half the inscribed radius so the offset triangle
  // doesn't collapse to a point.
  const double inscribed = r * std::cos(kPi / sides);
  if (std::fabs(delta) <= 1e-6) return;
  if (std::fabs(delta) > 0.4 * inscribed) return;

  const manifold::CrossSection input(RegularPolygon(sides, radius));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto expanded =
      input.Offset(delta, manifold::CrossSection::JoinType::Miter,
                   /*miter_limit=*/2.0, /*circular=*/0);
  ExpectCrossSectionValid(expanded);
  if (expanded.IsEmpty()) return;

  const auto roundTrip =
      expanded.Offset(-delta, manifold::CrossSection::JoinType::Miter,
                      /*miter_limit=*/2.0, /*circular=*/0);
  ExpectCrossSectionValid(roundTrip);

  // For convex Miter-join inputs the round-trip is exact in principle;
  // FP drift accumulates so 1e-4 relative tol catches real bugs without
  // false positives from arithmetic noise.
  const double tol = 1e-4 * (1.0 + std::fabs(input.Area()));
  EXPECT_NEAR(roundTrip.Area(), input.Area(), tol)
      << "Offset(" << delta << ", Miter).Offset(" << -delta
      << ", Miter) on " << sides << "-gon r=" << r
      << " didn't round-trip area";
  EXPECT_EQ(roundTrip.NumContour(), input.NumContour())
      << "Offset round-trip changed contour count";
}

void SimplePositiveOffset(const std::vector<double>& radii, double delta,
                          manifold::CrossSection::JoinType joinType,
                          int circularSegments) {
  const manifold::CrossSection input(StarPolygon(radii));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto output =
      input.Offset(delta, joinType, /*miter_limit=*/2.0, circularSegments);
  ExpectCrossSectionValid(output);
  EXPECT_FALSE(output.IsEmpty());
  EXPECT_GE(output.Area(), input.Area() - 1e-6 * (1.0 + input.Area()));
}

void OffsetExtrudeRoundTrip(const std::vector<double>& radii, double delta,
                            manifold::CrossSection::JoinType joinType,
                            int circularSegments, double height,
                            int nDivisions) {
  const manifold::CrossSection input(StarPolygon(radii));
  const auto output =
      input.Offset(delta, joinType, /*miter_limit=*/2.0, circularSegments);
  ExpectCrossSectionValid(output);
  if (output.IsEmpty() || std::fabs(output.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(output.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_NEAR(solid.Volume(), output.Area() * height,
              1e-6 * (1.0 + std::fabs(output.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
}

void ManifoldHoledExtrudeRoundTrip(const std::vector<double>& radii,
                                   double innerScale, double height,
                                   int nDivisions) {
  const manifold::CrossSection input(StarWithHole(radii, innerScale));
  ExpectCrossSectionValid(input);
  if (input.IsEmpty() || std::fabs(input.Area()) <= 1e-9) return;

  const auto solid =
      manifold::Manifold::Extrude(input.ToPolygons(), height, nDivisions);
  EXPECT_EQ(solid.Status(), manifold::Manifold::Error::NoError);
  EXPECT_TRUE(std::isfinite(solid.Volume()));
  EXPECT_EQ(solid.Genus(), 1);
  EXPECT_NEAR(solid.Volume(), input.Area() * height,
              1e-6 * (1.0 + std::fabs(input.Area() * height)));

  const manifold::CrossSection projected(solid.Project());
  ExpectCrossSectionValid(projected);
  EXPECT_NEAR(projected.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));

  const manifold::CrossSection middle(solid.Slice(height * 0.5));
  ExpectCrossSectionValid(middle);
  EXPECT_NEAR(middle.Area(), input.Area(),
              1e-6 * (1.0 + std::fabs(input.Area())));
}

std::vector<BooleanSeed> BooleanSeeds() {
  return {
      {{{{{0, 0}, {10, 0}, {10, 10}, {0, 10}}}},
       {{{{5, -1}, {11, 5}, {5, 11}, {-1, 5}}}},
       manifold::OpType::Add},
      {{{{{0, 0}, {10, 0}, {10, 1}, {0, 1}}}},
       {{{{5, 0}, {15, 0}, {15, 1}, {5, 1}}}},
       manifold::OpType::Intersect},
      {{{{{0, 0}, {8, 0}, {8, 8}, {0, 8}}, {{2, 6}, {6, 6}, {6, 2}, {2, 2}}}},
       {{{{4, -2}, {10, 4}, {4, 10}, {-2, 4}}}},
       manifold::OpType::Subtract},
  };
}

std::vector<OffsetSeed> OffsetSeeds() {
  return {
      {{{{{0, 0}, {20, 0}, {20, 20}, {0, 20}}}},
       3.0,
       manifold::CrossSection::JoinType::Round,
       16},
      {{{{{0, 0}, {20, 0}, {20, 2}, {2, 2}, {2, 20}, {0, 20}}}},
       -0.9,
       manifold::CrossSection::JoinType::Miter,
       0},
      {{{{{0, 0}, {12, 0}, {6, 0.01}, {12, 8}, {0, 8}}}},
       2.0,
       manifold::CrossSection::JoinType::Square,
       0},
  };
}

std::vector<ExtrudeSeed> ExtrudeSeeds() {
  return {
      {{{{{0, 0}, {10, 0}, {10, 10}, {0, 10}}}}, 5.0, 0},
      {{{{{0, 0}, {20, 0}, {20, 20}, {0, 20}},
         {{5, 15}, {15, 15}, {15, 5}, {5, 5}}}},
       3.0,
       2},
      {{{{{0, 0}, {12, 0}, {6, 0.01}, {12, 8}, {0, 8}}}}, 2.0, 1},
  };
}

std::vector<StarExtrudeSeed> StarExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}, 5.0, 0},
      {{1.0, 2.0, 0.5, 2.0, 1.0}, 3.0, 2},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0}, 2.0, 1},
  };
}

std::vector<BooleanExtrudeSeed> BooleanExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0},
       {1.0, 1.0, 1.0, 1.0},
       0.5,
       0.0,
       manifold::OpType::Add,
       5.0,
       0},
      {{1.0, 2.0, 0.5, 2.0, 1.0},
       {0.5, 1.5, 0.5, 1.5, 0.5},
       0.25,
       -0.25,
       manifold::OpType::Subtract,
       3.0,
       1},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0},
       {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
       0.0,
       0.0,
       manifold::OpType::Intersect,
       2.0,
       2},
  };
}

std::vector<PrismBooleanSeed> PrismBooleanSeeds() {
  return {
      {4, 1.0, 4, 1.0, 0.5, 0.0, manifold::OpType::Add},
      {5, 2.0, 6, 1.5, 0.25, -0.25, manifold::OpType::Subtract},
      {8, 3.0, 7, 2.0, 1.0, 0.0, manifold::OpType::Intersect},
      {17, 35.115498790314518, 7, 38.316457438114696, -0.0,
       -0.30479173448248414, manifold::OpType::Subtract},
  };
}

std::vector<TranslatedSliceSeed> TranslatedSliceSeeds() {
  return {
      {4, 1.0, 0.0, 0.0, 0.0, 5.0},
      {5, 2.0, 10.0, -20.0, 3.0, 2.0},
      {12, 3.0, -100.0, 50.0, -25.0, 10.0},
  };
}

std::vector<RotatedSliceSeed> RotatedSliceSeeds() {
  return {
      {4, 1.0, 0.0, 5.0},
      {5, 2.0, 45.0, 2.0},
      {12, 3.0, -90.0, 10.0},
  };
}

std::vector<HoledBooleanExtrudeSeed> HoledBooleanExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0},
       {1.0, 1.0, 1.0, 1.0},
       0.25,
       0.2,
       2.0,
       0.0,
       manifold::OpType::Add,
       5.0,
       0},
      {{1.0, 2.0, 0.5, 2.0, 1.0},
       {0.5, 1.5, 0.5, 1.5, 0.5},
       0.2,
       0.15,
       0.25,
       -0.25,
       manifold::OpType::Subtract,
       3.0,
       1},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0},
       {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
       0.1,
       0.3,
       0.0,
       0.0,
       manifold::OpType::Intersect,
       2.0,
       2},
      {{1.0, 4.0, 1.0, 4.0, 1.0, 4.0, 1.0, 4.0},
       {4.0, 1.0, 4.0, 1.0, 4.0, 1.0, 4.0, 1.0},
       0.1,
       0.1,
       1.0,
       1.0,
       manifold::OpType::Add,
       4.0,
       4},
      {{3.0, 0.5, 2.0, 0.5, 3.0, 0.5},
       {2.0, 0.5, 3.0, 0.5, 2.0, 0.5},
       0.15,
       0.15,
       -0.5,
       0.5,
       manifold::OpType::Subtract,
       1.5,
       2},
  };
}

std::vector<TransformExtrudeSeed> TransformExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}, 45.0, 2.0, 0.5, 5.0, 0},
      {{1.0, 2.0, 0.5, 2.0, 1.0}, -30.0, 1.5, 3.0, 3.0, 2},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0}, 90.0, 0.25, 4.0, 2.0, 1},
  };
}

std::vector<MultiComponentSeed> MultiComponentSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}, 2, 0.0},
      {{1.0, 2.0, 0.5, 2.0, 1.0}, 3, 500.0},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0}, 2, 1000.0},
  };
}

std::vector<BatchBooleanSeed> BatchBooleanSeeds() {
  return {
      {4, 1.0, 1, 0.0},
      {5, 2.0, 2, 50.0},
      {12, 3.0, 3, 100.0},
  };
}

std::vector<WarpAffineSeed> WarpAffineSeeds() {
  return {
      {4, 1.0, 2.0, 3.0, 4.0, 5.0, 0},
      {5, 2.0, 0.5, 4.0, -10.0, 20.0, 1},
      {12, 3.0, 10.0, 0.25, 100.0, -50.0, 2},
  };
}

std::vector<MirrorExtrudeSeed> MirrorExtrudeSeeds() {
  return {
      {4, 1.0, 1.0, 0.0, 5.0, 0},
      {5, 2.0, 1.0, 1.0, 3.0, 1},
      {12, 3.0, -1.0, 1.0, 2.0, 2},
      {6, 1.0, 0.0, 0.0, 1.0, 0},
  };
}

std::vector<StarSeed> StarSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}},
      {{1.0, 2.0, 0.5, 2.0, 1.0}},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0}},
  };
}

std::vector<SimpleOffsetSeed> SimpleOffsetSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}, 0.5, manifold::CrossSection::JoinType::Round, 16},
      {{1.0, 2.0, 0.5, 2.0, 1.0},
       0.25,
       manifold::CrossSection::JoinType::Miter,
       0},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0},
       0.1,
       manifold::CrossSection::JoinType::Bevel,
       0},
  };
}

std::vector<OffsetExtrudeSeed> OffsetExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0},
       0.5,
       manifold::CrossSection::JoinType::Round,
       16,
       5.0,
       0},
      {{1.0, 2.0, 0.5, 2.0, 1.0},
       -0.1,
       manifold::CrossSection::JoinType::Miter,
       0,
       3.0,
       1},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0},
       0.25,
       manifold::CrossSection::JoinType::Bevel,
       0,
       2.0,
       2},
  };
}

std::vector<HoleExtrudeSeed> HoleExtrudeSeeds() {
  return {
      {{1.0, 1.0, 1.0, 1.0}, 0.25, 5.0, 0},
      {{1.0, 2.0, 0.5, 2.0, 1.0}, 0.2, 3.0, 2},
      {{0.1, 3.0, 0.1, 3.0, 0.1, 3.0}, 0.1, 2.0, 1},
  };
}

auto CoordinateDomain() {
  return OneOf(InRange(-1000.0, 1000.0),
               ElementOf({-1024.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 1024.0}));
}

auto RawRingDomain() {
  return VectorOf(PairOf(CoordinateDomain(), CoordinateDomain()))
      .WithMinSize(3)
      .WithMaxSize(48);
}

auto RawPolygonsDomain() {
  return VectorOf(RawRingDomain()).WithMinSize(1).WithMaxSize(4);
}

auto StarRadiiDomain() {
  return VectorOf(InRange(0.0, 1000.0)).WithMinSize(4).WithMaxSize(48);
}

auto SmallStarRadiiDomain() {
  return VectorOf(InRange(0.0, 50.0)).WithMinSize(4).WithMaxSize(24);
}

}  // namespace

FUZZ_TEST(CrossSectionFuzz, BooleanRobustness)
    .WithDomains(RawPolygonsDomain(), RawPolygonsDomain(),
                 ElementOf({manifold::OpType::Add, manifold::OpType::Subtract,
                            manifold::OpType::Intersect}))
    .WithSeeds(BooleanSeeds);

FUZZ_TEST(CrossSectionFuzz, OffsetRobustness)
    .WithDomains(RawPolygonsDomain(), InRange(-100.0, 100.0),
                 ElementOf({manifold::CrossSection::JoinType::Square,
                            manifold::CrossSection::JoinType::Round,
                            manifold::CrossSection::JoinType::Miter,
                            manifold::CrossSection::JoinType::Bevel}),
                 ElementOf({0, 4, 8, 16, 32}))
    .WithSeeds(OffsetSeeds);

FUZZ_TEST(CrossSectionFuzz, ManifoldExtrudeRoundTrip)
    .WithDomains(RawPolygonsDomain(), InRange(0.1, 50.0),
                 ElementOf({0, 1, 2, 4}))
    .WithSeeds(ExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, ManifoldSimpleExtrudeRoundTrip)
    .WithDomains(StarRadiiDomain(), InRange(0.1, 50.0), ElementOf({0, 1, 2, 4}))
    .WithSeeds(StarExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, BooleanExtrudeRoundTrip)
    .WithDomains(SmallStarRadiiDomain(), SmallStarRadiiDomain(),
                 InRange(-100.0, 100.0), InRange(-100.0, 100.0),
                 ElementOf({manifold::OpType::Add, manifold::OpType::Subtract,
                            manifold::OpType::Intersect}),
                 InRange(0.1, 50.0), ElementOf({0, 1, 2, 4}))
    .WithSeeds(BooleanExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, PrismBooleanMatchesCrossSection)
    .WithDomains(InRange(3, 32), InRange(0.1, 50.0), InRange(3, 32),
                 InRange(0.1, 50.0), InRange(-75.0, 75.0), InRange(-75.0, 75.0),
                 ElementOf({manifold::OpType::Add, manifold::OpType::Subtract,
                            manifold::OpType::Intersect}))
    .WithSeeds(PrismBooleanSeeds);

FUZZ_TEST(CrossSectionFuzz, ManifoldTransformedExtrudeRoundTrip)
    .WithDomains(StarRadiiDomain(), InRange(-180.0, 180.0), InRange(0.1, 10.0),
                 InRange(0.1, 10.0), InRange(0.1, 50.0),
                 ElementOf({0, 1, 2, 4}))
    .WithSeeds(TransformExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, DecomposeComposeAndHull)
    .WithDomains(StarRadiiDomain(), ElementOf({1, 2, 3}), InRange(0.0, 2000.0))
    .WithSeeds(MultiComponentSeeds);

FUZZ_TEST(CrossSectionFuzz, BatchBooleanSeparated)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), ElementOf({1, 2, 3}),
                 InRange(0.0, 2000.0))
    .WithSeeds(BatchBooleanSeeds);

FUZZ_TEST(CrossSectionFuzz, DecomposedExtrusionsRecompose)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), ElementOf({1, 2, 3}),
                 InRange(0.0, 2000.0))
    .WithSeeds(BatchBooleanSeeds);

FUZZ_TEST(CrossSectionFuzz, TranslatedExtrusionSliceMatchesCrossSection)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), InRange(-1000.0, 1000.0),
                 InRange(-1000.0, 1000.0), InRange(-1000.0, 1000.0),
                 InRange(0.1, 50.0))
    .WithSeeds(TranslatedSliceSeeds);

FUZZ_TEST(CrossSectionFuzz, RotatedExtrusionSliceMatchesCrossSection)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), InRange(-180.0, 180.0),
                 InRange(0.1, 50.0))
    .WithSeeds(RotatedSliceSeeds);

FUZZ_TEST(CrossSectionFuzz, HoledBooleanExtrudeRoundTrip)
    .WithDomains(SmallStarRadiiDomain(), SmallStarRadiiDomain(),
                 InRange(0.05, 0.45), InRange(0.05, 0.45), InRange(-25.0, 25.0),
                 InRange(-25.0, 25.0),
                 ElementOf({manifold::OpType::Add, manifold::OpType::Subtract,
                            manifold::OpType::Intersect}),
                 InRange(0.1, 50.0), ElementOf({0, 1, 2, 4}))
    .WithSeeds(HoledBooleanExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, WarpAffineEquivalence)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), InRange(0.1, 10.0),
                 InRange(0.1, 10.0), InRange(-1000.0, 1000.0),
                 InRange(-1000.0, 1000.0), ElementOf({0, 1, 2}))
    .WithSeeds(WarpAffineSeeds);

FUZZ_TEST(CrossSectionFuzz, MirrorExtrudeRoundTrip)
    .WithDomains(InRange(3, 64), InRange(0.1, 100.0), InRange(-10.0, 10.0),
                 InRange(-10.0, 10.0), InRange(0.1, 50.0),
                 ElementOf({0, 1, 2, 4}))
    .WithSeeds(MirrorExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, SimpleBooleanIdentities)
    .WithDomains(StarRadiiDomain())
    .WithSeeds(StarSeeds);

FUZZ_TEST(CrossSectionFuzz, ApexSkipNearLine)
    .WithDomains(InRange(1e-15, 1e-1), InRange(-1.5, 1.5));

FUZZ_TEST(CrossSectionFuzz, TranslationInvariance)
    .WithDomains(StarRadiiDomain(), InRange(1e3, 1e9), InRange(1e3, 1e9));

FUZZ_TEST(CrossSectionFuzz, IterateToFixedPointConverges)
    .WithDomains(StarRadiiDomain());

FUZZ_TEST(CrossSectionFuzz, SubtractInvariants)
    .WithDomains(StarRadiiDomain(), StarRadiiDomain(), InRange(-5.0, 5.0),
                 InRange(-5.0, 5.0));

FUZZ_TEST(CrossSectionFuzz, BooleanCommutativity)
    .WithDomains(StarRadiiDomain(), StarRadiiDomain(), InRange(-5.0, 5.0),
                 InRange(-5.0, 5.0));

FUZZ_TEST(CrossSectionFuzz, BooleanAssociativity)
    .WithDomains(StarRadiiDomain(), StarRadiiDomain(), StarRadiiDomain(),
                 InRange(-5.0, 5.0), InRange(-5.0, 5.0), InRange(-5.0, 5.0),
                 InRange(-5.0, 5.0));

FUZZ_TEST(CrossSectionFuzz, BooleanDistributivity)
    .WithDomains(StarRadiiDomain(), StarRadiiDomain(), StarRadiiDomain(),
                 InRange(-5.0, 5.0), InRange(-5.0, 5.0), InRange(-5.0, 5.0),
                 InRange(-5.0, 5.0));

FUZZ_TEST(CrossSectionFuzz, ScaleInvariance)
    .WithDomains(StarRadiiDomain(), InRange(1e-4, 1e4));

FUZZ_TEST(CrossSectionFuzz, RotationInvariance)
    .WithDomains(StarRadiiDomain(), InRange(0.0, 6.2831853071795862));

FUZZ_TEST(CrossSectionFuzz, OffsetIdentityAtZero)
    .WithDomains(StarRadiiDomain(),
                 ElementOf({manifold::CrossSection::JoinType::Square,
                            manifold::CrossSection::JoinType::Round,
                            manifold::CrossSection::JoinType::Miter,
                            manifold::CrossSection::JoinType::Bevel}));

FUZZ_TEST(CrossSectionFuzz, EmptyIdentities)
    .WithDomains(StarRadiiDomain());

FUZZ_TEST(CrossSectionFuzz, DoubleMirrorIdentity)
    .WithDomains(StarRadiiDomain(), InRange(-10.0, 10.0),
                 InRange(-10.0, 10.0));

FUZZ_TEST(CrossSectionFuzz, DecomposeRecomposeWithHoles)
    .WithDomains(SmallStarRadiiDomain(), SmallStarRadiiDomain(),
                 InRange(-5.0, 5.0), InRange(-5.0, 5.0));

FUZZ_TEST(CrossSectionFuzz, OffsetInverseConvex)
    .WithDomains(InRange(4, 32), InRange(0.05, 25.0), InRange(-10.0, 10.0));

FUZZ_TEST(CrossSectionFuzz, HullIdempotence).WithDomains(StarRadiiDomain());

FUZZ_TEST(CrossSectionFuzz, SimplePositiveOffset)
    .WithDomains(SmallStarRadiiDomain(), InRange(0.05, 25.0),
                 ElementOf({manifold::CrossSection::JoinType::Square,
                            manifold::CrossSection::JoinType::Round,
                            manifold::CrossSection::JoinType::Miter,
                            manifold::CrossSection::JoinType::Bevel}),
                 ElementOf({0, 4, 8, 16, 32}))
    .WithSeeds(SimpleOffsetSeeds);

FUZZ_TEST(CrossSectionFuzz, OffsetExtrudeRoundTrip)
    .WithDomains(SmallStarRadiiDomain(), InRange(-5.0, 10.0),
                 ElementOf({manifold::CrossSection::JoinType::Square,
                            manifold::CrossSection::JoinType::Round,
                            manifold::CrossSection::JoinType::Miter,
                            manifold::CrossSection::JoinType::Bevel}),
                 ElementOf({0, 4, 8, 16, 32}), InRange(0.1, 50.0),
                 ElementOf({0, 1, 2, 4}))
    .WithSeeds(OffsetExtrudeSeeds);

FUZZ_TEST(CrossSectionFuzz, ManifoldHoledExtrudeRoundTrip)
    .WithDomains(StarRadiiDomain(), InRange(0.05, 0.45), InRange(0.1, 50.0),
                 ElementOf({0, 1, 2, 4}))
    .WithSeeds(HoleExtrudeSeeds);
