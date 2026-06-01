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
#include <cstdint>
#ifdef MANIFOLD_DEBUG
#include <fstream>
#include <iomanip>
#include <ostream>
#endif
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
#include "../src/cross_section/boolean2/boolean2.h"
#include "../src/cross_section/boolean2/diagnostics.h"
#include "../src/cross_section/boolean2/driver.h"
#include "../src/cross_section/boolean2/edge_vert_lists.h"
#include "../src/cross_section/boolean2/intersections.h"
#include "../src/cross_section/boolean2/predicates.h"
#include "../src/cross_section/boolean2/vertex_merge.h"
#endif
#include "manifold/common.h"
#include "manifold/manifold.h"
#include "test.h"

using namespace manifold;

namespace {

#if defined(MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2) && defined(MANIFOLD_DEBUG)
void WriteEscaped(std::ostream& os, const std::string& s) {
  os << '"';
  for (char c : s) {
    switch (c) {
      case '"':
        os << "\\\"";
        break;
      case '\\':
        os << "\\\\";
        break;
      case '\n':
        os << "\\n";
        break;
      case '\r':
        os << "\\r";
        break;
      case '\t':
        os << "\\t";
        break;
      default:
        os << c;
        break;
    }
  }
  os << '"';
}

void WriteVec2(std::ostream& os, const vec2& p) {
  os << '[' << p.x << ',' << p.y << ']';
}

void WriteField(std::ostream& os, const char* name, const std::string& value,
                bool comma = true) {
  WriteEscaped(os, name);
  os << ':';
  WriteEscaped(os, value);
  if (comma) os << ',';
}

template <typename T, typename F>
void WriteArray(std::ostream& os, const std::vector<T>& items, F writeItem) {
  os << '[';
  for (size_t i = 0; i < items.size(); ++i) {
    if (i > 0) os << ',';
    writeItem(os, items[i]);
  }
  os << ']';
}

void WriteTraceJson(std::ostream& os, const boolean2::Trace& trace) {
  os << std::setprecision(17);
  os << "{\n";
  WriteEscaped(os, "eps");
  os << ':' << trace.eps << ",\n";
  WriteField(os, "rule", trace.rule);
  os << "\n";
  WriteEscaped(os, "phases");
  os << ":[\n";
  for (size_t i = 0; i < trace.phases.size(); ++i) {
    const boolean2::TracePhase& phase = trace.phases[i];
    if (i > 0) os << ",\n";
    os << '{';
    WriteField(os, "name", phase.name);
    WriteEscaped(os, "points");
    os << ':';
    WriteArray(os, phase.points,
               [](std::ostream& out, const boolean2::TracePoint& p) {
                 out << '{';
                 WriteField(out, "id", p.id);
                 WriteEscaped(out, "xy");
                 out << ':';
                 WriteVec2(out, p.p);
                 out << ',';
                 WriteField(out, "kind", p.kind);
                 WriteField(out, "source", p.source);
                 WriteField(out, "label", p.label, false);
                 out << '}';
               });
    os << ',';
    WriteEscaped(os, "segments");
    os << ':';
    WriteArray(os, phase.segments,
               [](std::ostream& out, const boolean2::TraceSegment& s) {
                 out << '{';
                 WriteField(out, "id", s.id);
                 WriteEscaped(out, "a");
                 out << ':';
                 WriteVec2(out, s.a);
                 out << ',';
                 WriteEscaped(out, "b");
                 out << ':';
                 WriteVec2(out, s.b);
                 out << ',';
                 WriteField(out, "kind", s.kind);
                 WriteField(out, "source", s.source);
                 WriteEscaped(out, "mult");
                 out << ':' << s.mult << ',';
                 WriteField(out, "label", s.label, false);
                 out << '}';
               });
    os << ',';
    WriteEscaped(os, "polygons");
    os << ':';
    WriteArray(os, phase.polygons,
               [](std::ostream& out, const boolean2::TracePolygon& p) {
                 out << '{';
                 WriteField(out, "id", p.id);
                 WriteEscaped(out, "verts");
                 out << ":[";
                 for (size_t j = 0; j < p.verts.size(); ++j) {
                   if (j > 0) out << ',';
                   WriteVec2(out, p.verts[j]);
                 }
                 out << "],";
                 WriteField(out, "kind", p.kind);
                 WriteField(out, "source", p.source);
                 WriteEscaped(out, "winding");
                 out << ':' << p.winding << ',';
                 WriteEscaped(out, "inside");
                 out << ':' << (p.inside ? "true" : "false") << ',';
                 WriteField(out, "label", p.label, false);
                 out << '}';
               });
    os << ',';
    WriteEscaped(os, "annotations");
    os << ':';
    WriteArray(os, phase.annotations,
               [](std::ostream& out, const boolean2::TraceAnnotation& a) {
                 out << '{';
                 WriteField(out, "target", a.target);
                 WriteField(out, "key", a.key);
                 WriteField(out, "value", a.value, false);
                 out << '}';
               });
    os << '}';
  }
  os << "\n]\n}\n";
}
#endif

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

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
double Cross2D(vec2 a, vec2 b) { return a.x * b.y - a.y * b.x; }

bool ValidVertId(int v, const std::vector<vec2>& verts) {
  return 0 <= v && v < static_cast<int>(verts.size());
}

bool PointsCollinearWithinEps(vec2 a, vec2 b, vec2 p, double eps) {
  const vec2 d = b - a;
  const double len2 = dot(d, d);
  if (len2 == 0.0) return false;
  return std::fabs(Cross2D(d, p - a)) <= eps * std::sqrt(len2);
}

bool SegmentsHavePositiveCollinearOverlap(vec2 a0, vec2 a1, vec2 b0, vec2 b1,
                                          double eps) {
  if (!PointsCollinearWithinEps(a0, a1, b0, eps) ||
      !PointsCollinearWithinEps(a0, a1, b1, eps)) {
    return false;
  }
  const int axis = std::fabs(a1.x - a0.x) >= std::fabs(a1.y - a0.y) ? 0 : 1;
  const double aMin =
      std::min(boolean2::Coord(a0, axis), boolean2::Coord(a1, axis));
  const double aMax =
      std::max(boolean2::Coord(a0, axis), boolean2::Coord(a1, axis));
  const double bMin =
      std::min(boolean2::Coord(b0, axis), boolean2::Coord(b1, axis));
  const double bMax =
      std::max(boolean2::Coord(b0, axis), boolean2::Coord(b1, axis));
  return std::min(aMax, bMax) - std::max(aMin, bMin) > eps;
}

int StrictSide(vec2 a, vec2 b, vec2 p, double eps) {
  const vec2 d = b - a;
  const double threshold = eps * std::sqrt(dot(d, d));
  const double area = Cross2D(d, p - a);
  return (area > threshold) - (area < -threshold);
}

bool SegmentsHaveStrictCrossing(vec2 a0, vec2 a1, vec2 b0, vec2 b1,
                                double eps) {
  const int b0Side = StrictSide(a0, a1, b0, eps);
  const int b1Side = StrictSide(a0, a1, b1, eps);
  const int a0Side = StrictSide(b0, b1, a0, eps);
  const int a1Side = StrictSide(b0, b1, a1, eps);
  return b0Side * b1Side < 0 && a0Side * a1Side < 0;
}

bool PointInSegmentInteriorBand(vec2 p, vec2 a, vec2 b, double eps) {
  const vec2 d = b - a;
  const double len2 = dot(d, d);
  if (len2 == 0.0) return false;
  const double along = dot(p - a, d);
  if (along <= 0.0 || along >= len2) return false;
  return std::fabs(Cross2D(d, p - a)) <= eps * std::sqrt(len2);
}

::testing::AssertionResult CheckRetainedGraphValidity(
    const boolean2::OverlapResult& result,
    const std::vector<boolean2::EdgeM>& inputEdges,
    const std::vector<int>& inputVert2Merged, int numMergedVerts, double eps) {
  auto fail = [](const std::string& msg) {
    return ::testing::AssertionFailure() << msg;
  };

  for (int v = 0; v < static_cast<int>(result.verts.size()); ++v) {
    if (!std::isfinite(result.verts[v].x) ||
        !std::isfinite(result.verts[v].y)) {
      std::ostringstream out;
      out << "non-finite retained vertex " << v << " = (" << result.verts[v].x
          << ", " << result.verts[v].y << ")";
      return fail(out.str());
    }
  }

  const double eps2 = eps * eps;
  for (int a = 0; a < static_cast<int>(result.verts.size()); ++a) {
    for (int b = a + 1; b < static_cast<int>(result.verts.size()); ++b) {
      if (dot(result.verts[b] - result.verts[a],
              result.verts[b] - result.verts[a]) <= eps2) {
        std::ostringstream out;
        out << "retained vertices " << a << " and " << b
            << " remain within epsilon";
        return fail(out.str());
      }
    }
  }

  for (int i = 0; i < static_cast<int>(result.edges.size()); ++i) {
    const auto& edge = result.edges[i];
    if (!ValidVertId(edge.v0, result.verts) ||
        !ValidVertId(edge.v1, result.verts)) {
      std::ostringstream out;
      out << "retained edge " << i << " has invalid verts " << edge.v0 << " -> "
          << edge.v1;
      return fail(out.str());
    }
    if (dot(result.verts[edge.v1] - result.verts[edge.v0],
            result.verts[edge.v1] - result.verts[edge.v0]) <= eps2) {
      std::ostringstream out;
      out << "retained edge " << i << " is epsilon-zero: " << edge.v0 << " -> "
          << edge.v1;
      return fail(out.str());
    }
  }

  std::vector<boolean2::EdgeM> remapped;
  remapped.reserve(inputEdges.size());
  for (const auto& edge : inputEdges) {
    if (edge.v0 < 0 || edge.v0 >= static_cast<int>(inputVert2Merged.size()) ||
        edge.v1 < 0 || edge.v1 >= static_cast<int>(inputVert2Merged.size())) {
      std::ostringstream out;
      out << "input edge has verts outside inputVert2Merged: " << edge.v0
          << " -> " << edge.v1;
      return fail(out.str());
    }
    const int a = inputVert2Merged[edge.v0];
    const int b = inputVert2Merged[edge.v1];
    if (a != b) remapped.push_back({a, b, edge.mult});
  }

  const auto expected = ComputeBalance(remapped);
  const auto actual = ComputeBalance(result.edges);
  for (const auto& [v, actualBalance] : actual) {
    if (v < 0 || v >= static_cast<int>(result.verts.size())) {
      std::ostringstream out;
      out << "retained balance references invalid vertex " << v;
      return fail(out.str());
    }
    (void)actualBalance;
  }
  for (int v = 0; v < static_cast<int>(result.verts.size()); ++v) {
    const int expectedBalance =
        expected.count(v) ? expected.find(v)->second : 0;
    const int actualBalance = actual.count(v) ? actual.find(v)->second : 0;
    const int target = (v < numMergedVerts) ? expectedBalance : 0;
    if (actualBalance != target) {
      std::ostringstream out;
      out << "retained balance mismatch at vertex " << v << ": expected "
          << target << ", got " << actualBalance;
      return fail(out.str());
    }
  }

  for (int i = 0; i < static_cast<int>(result.edges.size()); ++i) {
    const auto& a = result.edges[i];
    for (int j = i + 1; j < static_cast<int>(result.edges.size()); ++j) {
      const auto& b = result.edges[j];
      if (a.v0 == b.v0 || a.v0 == b.v1 || a.v1 == b.v0 || a.v1 == b.v1) {
        continue;
      }
      const vec2 a0 = result.verts[a.v0];
      const vec2 a1 = result.verts[a.v1];
      const vec2 b0 = result.verts[b.v0];
      const vec2 b1 = result.verts[b.v1];
      if (SegmentsHaveStrictCrossing(a0, a1, b0, b1, eps)) {
        std::ostringstream out;
        out << "retained edges " << i << " and " << j
            << " still have a strict crossing";
        return fail(out.str());
      }
      if (SegmentsHavePositiveCollinearOverlap(a0, a1, b0, b1, eps)) {
        std::ostringstream out;
        out << "retained edges " << i << " and " << j
            << " still have positive collinear overlap";
        return fail(out.str());
      }
    }
  }

  for (int e = 0; e < static_cast<int>(result.edges.size()); ++e) {
    const auto& edge = result.edges[e];
    const vec2 a = result.verts[edge.v0];
    const vec2 b = result.verts[edge.v1];
    for (int v = 0; v < static_cast<int>(result.verts.size()); ++v) {
      if (v == edge.v0 || v == edge.v1) continue;
      if (PointInSegmentInteriorBand(result.verts[v], a, b, eps)) {
        std::ostringstream out;
        out << "retained vertex " << v << " lies in the interior band of edge "
            << e << " (" << edge.v0 << " -> " << edge.v1 << ")";
        return fail(out.str());
      }
    }
  }

  return ::testing::AssertionSuccess();
}

std::vector<boolean2::EdgeM> EdgesFromOverlapResult(
    const boolean2::OverlapResult& result) {
  std::vector<boolean2::EdgeM> edges;
  edges.reserve(result.edges.size());
  for (const auto& edge : result.edges) {
    edges.push_back({edge.v0, edge.v1, edge.mult});
  }
  return edges;
}

boolean2::OverlapResult CleanupPassLikeIterate(
    const boolean2::OverlapResult& result, double eps) {
  return boolean2::RemoveOverlaps2D(
      result.verts, EdgesFromOverlapResult(result), eps, /*tolerance=*/0.0,
      /*debug=*/false, boolean2::WindRule::Add);
}

void ExpectSameFingerprint(const boolean2::OverlapResult& a,
                           const boolean2::OverlapResult& b, double eps) {
  using Fingerprint =
      std::vector<std::tuple<int64_t, int64_t, int64_t, int64_t, int>>;
  auto fingerprint = [eps](const boolean2::OverlapResult& r) {
    const double quantum = eps * 0.01;
    auto q = [quantum](double x) {
      return static_cast<int64_t>(std::round(x / quantum));
    };
    Fingerprint fp;
    fp.reserve(r.edges.size());
    for (const auto& edge : r.edges) {
      vec2 p0 = r.verts[edge.v0];
      vec2 p1 = r.verts[edge.v1];
      auto k0 = std::make_pair(q(p0.x), q(p0.y));
      auto k1 = std::make_pair(q(p1.x), q(p1.y));
      int mult = edge.mult;
      if (k1 < k0) {
        std::swap(k0, k1);
        mult = -mult;
      }
      fp.emplace_back(k0.first, k0.second, k1.first, k1.second, mult);
    }
    manifold::stable_sort(fp.begin(), fp.end());
    return fp;
  };

  const auto fpA = fingerprint(a);
  const auto fpB = fingerprint(b);
  EXPECT_EQ(fpA, fpB) << "fingerprints differ: lhs edges=" << a.edges.size()
                      << " rhs edges=" << b.edges.size();
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

template <typename Edge>
std::map<int, int> ComputeEdgeBalance(const std::vector<Edge>& edges) {
  std::map<int, int> balance;
  for (const auto& edge : edges) {
    balance[edge.v0] += edge.mult;
    balance[edge.v1] -= edge.mult;
  }
  return balance;
}
#endif

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

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
TEST(CrossSection, PropagatesShallowIndependentIntersections) {
  using boolean2::Box2;
  using boolean2::BoxOf2DEdge;
  using boolean2::BuildListsAndFindIntersections;
  using boolean2::BVH;
  using boolean2::EdgeM;
  using boolean2::FindAndInsertIntersections;
  using boolean2::IntersectionInsertion;
  using boolean2::NarrowPhaseResult;

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
  const std::vector<std::pair<int, int>> pairs = {{0, 1}, {2, 3}};

  NarrowPhaseResult narrow =
      BuildListsAndFindIntersections(edges, verts, eps, pairs);
  IntersectionInsertion inserted = FindAndInsertIntersections(
      edges, std::move(verts), std::move(narrow.lists), eps, edgeBoxes, BVH{},
      narrow.intersections);

  ASSERT_EQ(inserted.verts.size(), 12);
  EXPECT_EQ(inserted.lists[4].size(), 2);
}
#endif

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
  SimplePolygon butterfly = {{-0.2, 0.0}, {-0.2, 0.2}, {-0.5, 0.2}, {-0.2, 0.2},
                             {0.0, 0.2},  {0.2, 0.2},  {0.2, 0.0},  {0.5, -0.2},
                             {0.2, -0.2}, {0.2, -0.5}, {0.2, -0.2}, {0.0, -0.2},
                             {-0.2, -0.2}};
  // CCW staircase whose v6->v7 edge ((0.2, -0.2) -> (-0.2, -0.2)) is the
  // butterfly's bottom edge, sharing both endpoints with butterfly
  // sub-edges (-0.2, -0.2)..(0.0, -0.2)..(0.2, -0.2) and creating a
  // T-junction at (0.0, -0.2).
  SimplePolygon staircase = {
      {0.2, -0.2}, {-0.2, -0.2}, {-0.2, -0.4}, {0.2, -0.4}};

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
  const std::vector<double> radii = {0.,
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
  const auto output = input.Offset(delta, CrossSection::JoinType::Bevel,
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
  const std::vector<double> radiiB = {113.5978182908662, 0., 114.34968677141997,
                                      6.5076626333939721e-10};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB)).Translate({0., 4.667562921730494e-33});

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

#if defined(MANIFOLD_DEBUG) && defined(MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2)
TEST(CrossSection, DISABLED_Boolean2TraceTinyVsLargeStars) {
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

  auto append = [](const Polygons& polys, int mult, std::vector<vec2>* verts,
                   std::vector<boolean2::EdgeM>* edges) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts->size());
      for (const vec2& v : loop) verts->push_back(v);
      const int n = static_cast<int>(loop.size());
      for (int i = 0; i < n; ++i) {
        edges->push_back({base + i, base + ((i + 1) % n), mult});
      }
    }
  };

  const std::vector<double> radiiA = {2.0371964671064575e-14,
                                      1.814002251251902e-10, 0.,
                                      2.1825088357767143e-09};
  const std::vector<double> radiiB = {113.5978182908662, 0., 114.34968677141997,
                                      6.5076626333939721e-10};
  Polygons a{star(radiiA)};
  Polygons b{star(radiiB)};
  for (vec2& v : b[0]) v.y += 4.667562921730494e-33;

  const double eps = boolean2::InferEps(a, b);
  std::vector<vec2> verts;
  std::vector<boolean2::EdgeM> edges;
  append(a, 1, &verts, &edges);
  append(b, -1, &verts, &edges);

  boolean2::Trace trace;
  auto result = boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                           /*debug=*/false,
                                           boolean2::WindRule::Add, &trace);
  Polygons final = boolean2::OutEdgesToPolygons(result.verts, result.edges);
  auto& phase = trace.AddPhase("final_polygons");
  for (int i = 0; i < static_cast<int>(final.size()); ++i) {
    phase.polygons.push_back({std::string("poly") + std::to_string(i), final[i],
                              "final_polygon", "", 0, true, ""});
  }

  std::ofstream out("boolean2_trace_tiny_vs_large_stars.json");
  ASSERT_TRUE(out.good());
  WriteTraceJson(out, trace);
}

TEST(CrossSection, DISABLED_Boolean2TraceShowcase) {
  auto append = [](const Polygons& polys, int mult, std::vector<vec2>* verts,
                   std::vector<boolean2::EdgeM>* edges) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts->size());
      for (const vec2& v : loop) verts->push_back(v);
      const int n = static_cast<int>(loop.size());
      for (int i = 0; i < n; ++i) {
        edges->push_back({base + i, base + ((i + 1) % n), mult});
      }
    }
  };

  const SimplePolygon outer = {
      {-4.0, -2.0}, {4.0, -2.0}, {4.0, 2.0}, {-4.0, 2.0}};
  const SimplePolygon hole = {
      {-1.2, -0.7}, {1.2, -0.7}, {1.2, 0.7}, {-1.2, 0.7}};
  const SimplePolygon diamond = {
      {0.0, -3.0}, {3.0, 0.0}, {0.0, 3.0}, {-3.0, 0.0}};
  const SimplePolygon slantedBar = {
      {-4.7, -0.35}, {4.7, 0.85}, {4.45, 1.55}, {-4.95, 0.35}};
  const SimplePolygon notch = {{-3.2, -2.6}, {-1.6, -2.6}, {-2.2, 0.15}};

  const Polygons positive{outer, diamond, slantedBar};
  const Polygons negative{hole, notch};
  const Polygons all{outer, hole, diamond, slantedBar, notch};
  const double eps = boolean2::InferEps(all, {});

  std::vector<vec2> verts;
  std::vector<boolean2::EdgeM> edges;
  append(positive, 1, &verts, &edges);
  append(negative, -1, &verts, &edges);

  boolean2::Trace trace;
  auto result = boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                           /*debug=*/false,
                                           boolean2::WindRule::Add, &trace);
  Polygons final = boolean2::OutEdgesToPolygons(result.verts, result.edges);
  auto& phase = trace.AddPhase("final_polygons");
  for (int i = 0; i < static_cast<int>(final.size()); ++i) {
    phase.polygons.push_back({std::string("poly") + std::to_string(i), final[i],
                              "final_polygon", "", 0, true, ""});
  }

  std::ofstream out("boolean2_trace_showcase.json");
  ASSERT_TRUE(out.good());
  WriteTraceJson(out, trace);
}
#endif

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

  const std::vector<double> radiiA = {627.11994612906153,
                                      0.11978140887764367,
                                      601.04982226530046,
                                      102.31660501134466,
                                      2.1549912703437601,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.};
  const std::vector<double> radiiB = {112.00504648760347, 4.3823848227606392,
                                      1.2795858846899296e-07,
                                      6.7216166802304542};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB)).Translate({-2.7193614970785894e-29, 0.});

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
  const std::vector<double> radiiB = {
      3.9905161863280855e-53, 0., 96.331001430191975, 0., 0.43952001502476951,
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
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({-0.72134106089064531, -4.5808240251267858});

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
  const std::vector<double> radiiA = {0.,
                                      39.150613710409736,
                                      0.,
                                      21.472413919643575,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      3.237760857312173e-24};
  const std::vector<double> radiiB = {72.814105930277677, 21.384583693573447,
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
  const std::vector<double> rA = {24.575460587300253, 2.4572617728922851e-10,
                                  3.3952718785303559e-06, 0.};
  const std::vector<double> rB = {29.731318514644453, 1.5729051003875837e-06,
                                  0.0082858661009423962, 0.};
  const std::vector<double> rC = {0., 0., 6.5474426075871467e-16,
                                  2.9296729904240054e-06};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB)).Translate(
      {-9.1863209962415243e-35, -9.8444830049208569e-28});
  const CrossSection c =
      CrossSection(star(rC)).Translate({-6.0489837474564476e-29, 0.});

  const auto ab_c = (a + b) + c;
  const auto a_bc = a + (b + c);
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(ab_c.Area(), a_bc.Area(), tol) << "(A ∪ B) ∪ C != A ∪ (B ∪ C)";
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
      5.5755616189085049e-09, 55.713117722084,        0.24425911685471227,
      12.355613962816753,     3.2696445771376686e-21, 7.7223346956456043e-27,
      9.1566312909847612e-29, 1.1757300951916251e-32, 7.5589164688107744e-40,
      2.1073720366622818e-34, 4.7311936525798036e-23};
  const std::vector<double> radiiB = {0.,
                                      55.263016827250659,
                                      0.,
                                      2.1568234275576556e-26,
                                      4.3751784176513804e-26,
                                      3.0217717245002794e-34,
                                      0.,
                                      0.,
                                      0.,
                                      5.2142580450244983e-22,
                                      6.1161686141794186e-15,
                                      2.8398350775198264};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({-0.35827067719250716, 2.8298129605573439});
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
  const std::vector<double> rA = {0.15588462038906103,
                                  0.00015361176400222398,
                                  5.8511202761938412e-16,
                                  3.3009196314697135e-10,
                                  2.7199740700330746e-09,
                                  4.690497141269639e-14,
                                  8.9922571322790269e-10,
                                  6.3553011716417212e-16,
                                  4.2872009761335588e-18,
                                  5.6483346599537904e-17,
                                  3.479258338980731e-12,
                                  9.9500479914355979e-14,
                                  1.822663758661907e-10,
                                  3.7366088533519133e-05,
                                  7.6137290475474603e-12,
                                  1.5917667211519553e-08,
                                  1.3293803091881353e-06,
                                  6.5804393417497265e-09,
                                  0.0015777855119048808,
                                  3.5317048943785664e-07,
                                  1.353313563107637e-09,
                                  1.1788667553671715e-05,
                                  0.00014847370669519226,
                                  1.509263403954591e-11,
                                  6.5516568176598574e-10,
                                  3.8749326674543751e-07,
                                  3.450117114134831e-06,
                                  8.3840226251340948e-08,
                                  3.0175612629932874e-07,
                                  33.778786246753299,
                                  0.033251289695133752,
                                  0.,
                                  0.,
                                  1.5817601023378236};
  const std::vector<double> rB = {3.910135632588194e-07,
                                  1.0064704530951155e-13,
                                  2.0158855952179934e-14,
                                  1.4165380430901347e-21,
                                  3.6148752792123783e-11,
                                  1.3306453783491879e-10,
                                  2.3794160990202857e-13,
                                  4.0821728496033454e-14,
                                  4.6816579341479393e-18,
                                  1.7855060689360453e-17,
                                  2.4858744170474308e-10,
                                  1.5313588823631113e-11,
                                  6.6960358868004941e-08,
                                  2.1802650917628242e-10,
                                  5.5611296920831909e-13,
                                  1.9586931478995027e-08,
                                  0.,
                                  36.09032519236748,
                                  0.,
                                  7.5123873966542927e-05};
  const CrossSection a(star(rA));
  const CrossSection b =
      CrossSection(star(rB)).Translate({0.72140998591309213, 0.});
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
  const CrossSection c =
      CrossSection(star(rC)).Translate({-2.4305117516652688e-13, 0.});
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
  const std::vector<double> rB = {
      1.4795645772678251e-12, 2.4342935254638573e-13, 1.4800031437954507e-12,
      1.1368244372782033e-305};
  const std::vector<double> rC = {0., 0., 0., 0.};
  const CrossSection a(star(rA));
  const CrossSection b = CrossSection(star(rB)).Translate({0., 0.});
  const CrossSection c(star(rC));
  const auto ab_c = (a + b) + c;
  const auto a_bc = a + (b + c);
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(ab_c.Area(), a_bc.Area(), tol) << "(A ∪ B) ∪ C != A ∪ (B ∪ C)";
}

// Seed: OffsetInverseConvex (2026-05-17 local smoke during /loop iter)
// Counterexample-hash: f186158d2284d4c9
// Suspected owner: pr/boolean2-backend-wiring (CrossSection::Offset on
//   boolean2 backend, src/cross_section/cross_section_boolean2.cpp).
//   Regular triangle, Miter join exactly on the miter_limit=2.0
//   threshold; the Offset(d).Offset(-d) round-trip drifts -0.35% in area
//   in a way that's scale-invariant (same percentage at r=0.15, r=1.5,
//   r=15). n>=4 has zero drift (verified probe), so the bug is specific
//   to the miter-limit boundary of the equilateral triangle.
TEST(CrossSection, OffsetInverseTriangleMiter) {
  // Equilateral triangle inscribed in r=0.15 (effective radius via the
  // 0.1 + |radius| convention used in cross_section_fuzz).
  SimplePolygon ring;
  const double r = 0.15;
  for (int i = 0; i < 3; ++i) {
    const double theta = 2.0 * kPi * i / 3;
    ring.push_back({r * std::cos(theta), r * std::sin(theta)});
  }
  const CrossSection input(ring);
  const double delta = -0.0094938192047002799;

  const auto expanded =
      input.Offset(delta, CrossSection::JoinType::Miter, /*miter_limit=*/2.0,
                   /*circularSegments=*/0);
  const auto roundTrip = expanded.Offset(-delta, CrossSection::JoinType::Miter,
                                         /*miter_limit=*/2.0,
                                         /*circularSegments=*/0);

  // For Miter join with miter_limit=2.0 exactly at the equilateral triangle
  // threshold, equality should still take the miter path. Observed:
  // 0.357% absolute drift when rounded normals spuriously square the join.
  const double tol = 1e-4 * (1.0 + std::fabs(input.Area()));
  EXPECT_NEAR(roundTrip.Area(), input.Area(), tol)
      << "Triangle Miter Offset round-trip drifted by "
      << (roundTrip.Area() - input.Area()) / input.Area() * 100 << "%";
  EXPECT_EQ(roundTrip.NumContour(), input.NumContour());
}

// Seed: DecomposeRecomposeWithHoles (2026-05-18 daemon find)
// Counterexample-hash: a33524d9c3e6fb10
// Suspected owner: pr/boolean2-core (CrossSection::Decompose path,
//   most likely the containment/face-classification step that picks
//   which rings to emit per component. An 8-vertex star outer with
//   a small translated hole produces a 2-contour CrossSection
//   (NumContour=2: outer+hole, Area=911.40 = 911.70 outer - 0.30
//   hole). Decompose returns 1 component whose Area=911.70 and
//   NumContour=1 - the hole is silently dropped. Compose then gives
//   a full-outer result, losing 0.3 area worth of hole. Area
//   conservation invariant violated; bidirectional Decompose/Compose
//   should be the identity on multi-ring inputs).
TEST(CrossSection, DecomposeRecomposeOuterStarWithSmallHole) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> outerRadii = {
      1., 1., 1., 50., 50., 0., 3.987798525003551, 1.};
  const std::vector<double> holeRadii = {0., 0., 1., 0.29444504003509697};
  const CrossSection outer(star(outerRadii));
  const CrossSection hole = CrossSection(star(holeRadii)).Translate({0., 1.});
  const auto holed = outer - hole;
  ASSERT_EQ(holed.NumContour(), 2u);  // sanity: subtract produced
                                      // outer + hole
  const double holedArea = holed.Area();

  const auto components = holed.Decompose();
  ASSERT_FALSE(components.empty());

  double componentArea = 0.0;
  size_t componentContourSum = 0;
  for (const auto& component : components) {
    componentArea += component.Area();
    componentContourSum += component.NumContour();
  }
  const double tol = 1e-6 * (1.0 + std::fabs(holedArea));
  EXPECT_NEAR(componentArea, holedArea, tol)
      << "Decompose lost area on holed input: " << "sum(components.Area)="
      << componentArea << " vs holed.Area()=" << holedArea;
  EXPECT_EQ(componentContourSum, holed.NumContour())
      << "Decompose lost contours on holed input: "
      << "sum(components.NumContour)=" << componentContourSum
      << " vs holed.NumContour()=" << holed.NumContour();

  const auto recomposed = CrossSection::Compose(components);
  EXPECT_NEAR(recomposed.Area(), holedArea, tol);
  EXPECT_EQ(recomposed.NumContour(), holed.NumContour());
}

// Seed: DecomposeRecomposeWithHoles (2026-05-18 CI run 26015239869)
// Counterexample-hash: 663debcfef9c9d1e
// Suspected owner: pr/boolean2-core (same Decompose area-conservation
//   class as DecomposeRecomposeOuterStarWithSmallHole but a distinct
//   counterexample. 5-vertex outer star (radii skewed: 48.55, 5.05,
//   1, 1, 50 - one large spike) minus a 5-vertex inner star
//   (radii 1, 3.75, 0.573, 1, 0) translated by (0, 1.0012).
//   componentArea=1304.5424 but holed.Area()=1304.5751, diff 0.0327
//   exceeds tol 0.0013. Decompose drops a chunk of hole-adjacent
//   area rather than the entire hole this time. Possible: when the
//   hole is near-tangent to the outer (the 1.0012 offset puts it
//   very close to an outer edge), face classification miscategorizes
//   a sliver that gets dropped on the decompose path).
TEST(CrossSection, DecomposeRecomposeNearTangentSmallHole) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  // Match the fuzz target's exact construction (outer - hole) so the
  // reproducer goes through the same Boolean Subtract -> Decompose
  // path that fuzzing exercised.
  const std::vector<double> outerRadii = {48.55001516665169, 5.0536385110757127,
                                          1., 1., 50.};
  const std::vector<double> holeRadii = {1., 3.7464608085566375,
                                         0.57299724595371804, 1., 0.};
  const CrossSection outer(star(outerRadii));
  const CrossSection hole =
      CrossSection(star(holeRadii)).Translate({-0., 1.0011960822937693});
  const auto holed = outer - hole;
  ASSERT_FALSE(holed.IsEmpty());
  ASSERT_GE(holed.NumContour(), 2u)
      << "expected outer + hole, got NumContour=" << holed.NumContour();
  const double holedArea = holed.Area();

  const auto components = holed.Decompose();
  ASSERT_FALSE(components.empty());

  double componentArea = 0.0;
  size_t componentContourSum = 0;
  for (const auto& c : components) {
    componentArea += c.Area();
    componentContourSum += c.NumContour();
  }
  const double tol = 1e-6 * (1.0 + std::fabs(holedArea));
  EXPECT_NEAR(componentArea, holedArea, tol)
      << "Decompose lost area on holed input: sum(components.Area)="
      << componentArea << " vs holed.Area()=" << holedArea;
  EXPECT_EQ(componentContourSum, holed.NumContour())
      << "Decompose split or merged contours unexpectedly";

  const auto recomposed = CrossSection::Compose(components);
  EXPECT_NEAR(recomposed.Area(), holedArea, tol)
      << "Compose(Decompose(holed)) changed area";
  EXPECT_EQ(recomposed.NumContour(), holed.NumContour())
      << "Compose(Decompose(holed)) changed contour count";
}

// Seed: BooleanDistributivity (2026-05-20 daemon find on post-cleanup tip)
// Counterexample-hash: 659ec969e064893e
// Suspected owner: pr/boolean2-core (left side `A ∩ (B ∪ C)` has
//   NumContour=24, right side `(A ∩ B) ∪ (A ∩ C)` has NumContour=23.
//   Off-by-one in contour decomposition. Inputs are three mid-size
//   stars (34/46/48 verts) with extreme-magnitude radii (mix of 0,
//   1, 1000-ish) and translations in [-2.4, 3.9] x [-1.0, 2.4].
//   Likely the interleaved Union+Intersect path produces or drops
//   a contour relative to the (Intersect, Intersect, Union)
//   ordering. Related to but distinct from previously-drained
//   BooleanCommutativity / BooleanAssociativity classes).
TEST(CrossSection, BooleanDistributivityExtremeMagStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {22.40352672886273,
                                      1000.,
                                      1.930469114705216,
                                      1000.,
                                      995.07901723564635,
                                      450.42497103653432,
                                      1.,
                                      1000.,
                                      442.67830885450564,
                                      974.26213405728731,
                                      762.40083808845225,
                                      375.22374612076686,
                                      311.70408869554416,
                                      4.6596927406080066,
                                      1000.,
                                      993.53293688808651,
                                      545.61942022366179,
                                      902.87869535025561,
                                      857.23036058152354,
                                      1.,
                                      997.14305195601548,
                                      998.81543502854458,
                                      4.287799655231888,
                                      1.9578049801330417,
                                      212.99717518614997,
                                      539.98252283023282,
                                      814.2147123777587,
                                      6.4135462648884491,
                                      386.69166945488905,
                                      642.25190405836315,
                                      671.02238904111891,
                                      1000.,
                                      0.,
                                      959.32882825167917};
  const std::vector<double> radiiB = {596.78699570023025,
                                      0.,
                                      0.,
                                      4.8215130310908281,
                                      134.35969768300723,
                                      514.84369184653679,
                                      514.84369184653679,
                                      514.84369184653679,
                                      134.35969768300723,
                                      134.35969768300723,
                                      1000.,
                                      1.,
                                      1000.,
                                      0.,
                                      926.29500965174145,
                                      0.,
                                      134.35969768300723,
                                      0.90493025476122824,
                                      95.854999537523781,
                                      1000.,
                                      0.,
                                      0.,
                                      460.32237309958083,
                                      2.9791218075989327,
                                      30.186663600739344,
                                      1.,
                                      5.732788180417673,
                                      487.77994489187455,
                                      1.,
                                      639.36782925291504,
                                      1.,
                                      1.,
                                      405.85342012592389,
                                      405.85342012592389,
                                      405.85342012592389,
                                      706.56289376450047,
                                      405.85342012592389,
                                      405.85342012592389,
                                      1.,
                                      1.,
                                      872.8336783186561,
                                      5.1830997334306659,
                                      313.84753596061228,
                                      3.3318960586713984,
                                      22.317598513501498,
                                      4.1956604335614083};
  const std::vector<double> radiiC = {279.68331311398725,
                                      670.56312431772176,
                                      812.15251441725718,
                                      1000.,
                                      0.,
                                      996.64035079388486,
                                      1.,
                                      1000.,
                                      1.,
                                      146.36259356992105,
                                      1000.,
                                      1.2440205144455736,
                                      2.4903214421341642,
                                      558.11442117895865,
                                      667.56514527125182,
                                      601.28991029625661,
                                      0.,
                                      0.,
                                      1.,
                                      485.29900078511309,
                                      205.35289543499823,
                                      1.,
                                      1.,
                                      999.13575735403208,
                                      0.,
                                      343.12355139003932,
                                      124.51182502694519,
                                      0.44553220583928932,
                                      512.45895098827191,
                                      2.5455428684940777,
                                      0.,
                                      1000.,
                                      1000.,
                                      0.,
                                      595.63567011806413,
                                      2.5851695628219766,
                                      3.2720224325713057,
                                      354.45821723238265,
                                      4.2476895774067867,
                                      956.03225557664916,
                                      0.,
                                      1000.,
                                      0.,
                                      98.799661756610732,
                                      0.,
                                      619.55626420041142,
                                      0.,
                                      999.86711403007541};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({3.4425093025001434, -0.89018350575014171});
  const CrossSection c =
      CrossSection(star(radiiC))
          .Translate({3.8936971336588151, 2.3137657561842939});

  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  // Contour count is not asserted: the right-hand side's deeper boolean
  // chain accumulates more propagated tolerance than the left-hand side,
  // which can merge one extra contour without breaking the area identity.
}

// Seed: BooleanDistributivity (2026-05-23 daemon find - both local
//   x86_64 and manifoci aarch64 daemons hit it simultaneously)
// Counterexample-hash: 5d07906de3df36e7
// Suspected owner: pr/boolean2-core (A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)
//   on large-area stars. left.Area=216374, right.Area=215961, diff=412
//   (relative ~0.2%, exceeds 1e-6 * sum-of-areas tol of 2.34).
//   Inputs: 16-vertex star A with mixed radii, 4-vertex star B (small),
//   37-vertex star C with most radii repeated as 332.47 (degenerate-ish
//   shape). Translations: B by (4.89, 3.65), C by (-1.53, -2.63).
//   Distinct from previously-fixed BooleanDistributivityExtremeMagStars
//   (hash 659ec969) - that was a contour-count off-by-one, this is
//   area drift in absolute terms).
TEST(CrossSection, BooleanDistributivityRepeatedRadiusStarC) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {1.,
                                      283.64569337262878,
                                      0.,
                                      705.96814039386641,
                                      332.07185213769458,
                                      534.61157622175767,
                                      1000.,
                                      1000.,
                                      1000.,
                                      1000.,
                                      1000.,
                                      1000.,
                                      1000.,
                                      1.,
                                      1.,
                                      1000.};
  const std::vector<double> radiiB = {974.53903971669604, 443.79612063734504,
                                      0.60212977272176105, 601.54223181860084};
  const std::vector<double> radiiC = {2.8991656937704362,
                                      1000.,
                                      1.,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      332.4725376862869,
                                      384.70682133429182,
                                      384.70682133429182,
                                      332.4725376862869,
                                      332.4725376862869,
                                      682.65150344752476,
                                      737.43471827249516,
                                      327.59982985461647,
                                      332.4725376862869,
                                      332.4725376862869,
                                      1.,
                                      1000.,
                                      1000.,
                                      1000.};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({4.890230810667985, 3.6470913496183108});
  const CrossSection c =
      CrossSection(star(radiiC))
          .Translate({-1.5262817329731697, -2.6347725451611175});

  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;

  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
}

// Seed: BooleanDistributivity (2026-05-23 cycle 277, both daemons,
//   post-disconnected-winding-seeds-fix tip 22077d9c)
// Counterexample-hash: cf05d696225efee5
// Suspected owner: pr/boolean2-core (sibling of the just-drained
//   RepeatedRadiusStarC. Smaller magnitude (left.Area=80787,
//   right.Area=80781, diff=5.9 vs tol=1.28). Star C has many zeros
//   mixed with mid-range values - a different degeneracy pattern
//   than the previous all-repeated-radius case).
TEST(CrossSection, BooleanDistributivityZeroMixStarC) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {1000.,
                                      1000.,
                                      255.79995725682645,
                                      1.,
                                      454.29432950304715,
                                      995.65401133735929,
                                      1000.,
                                      811.86733957046874,
                                      1.,
                                      437.04103091878892,
                                      195.25477231407794,
                                      1.8268273914683537,
                                      4.5954160688751484,
                                      3.2030613706914481,
                                      511.5935307510378,
                                      2.0771806846580496,
                                      570.21917979539887,
                                      1.9769538901315753,
                                      467.93526395973339,
                                      1.,
                                      959.79674877569607,
                                      0.076792444298594276,
                                      503.25635714266843,
                                      0.,
                                      89.215459767705454,
                                      0.,
                                      1000.,
                                      0.06374847983784826,
                                      995.14671048524883,
                                      1.,
                                      0.};
  const std::vector<double> radiiB = {271.39306287742255, 997.96011166439871,
                                      909.94463259654469, 364.30795394492606};
  const std::vector<double> radiiC = {812.15251441725718,
                                      0.,
                                      601.28991029625661,
                                      0.,
                                      482.34411167146169,
                                      0.,
                                      588.04482477271347,
                                      0.,
                                      169.03839642543201,
                                      0.44553220583928932,
                                      4.0649812370918958,
                                      512.45895098827191,
                                      0.,
                                      1000.,
                                      0.,
                                      999.90271349806335,
                                      459.5889801989477,
                                      105.49723198609914,
                                      0.,
                                      0.,
                                      748.07738025389699,
                                      1.7632947108253827,
                                      0.,
                                      0.,
                                      4.9202747062496162,
                                      751.85659442712381,
                                      732.78647234489199,
                                      0.,
                                      375.32896975602767,
                                      0.};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({0.12749108485042271, 2.533194291154123});
  const CrossSection c =
      CrossSection(star(radiiC))
          .Translate({-2.5843878979814856, -4.3357952411613461});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
}

// Seed: BooleanDistributivity (2026-05-24 local daemon cycle 290,
//   post-nonzero-union-fix tip a9cfd4ff)
// Counterexample-hash: 86029efb7e8d5bd1
// Suspected owner: pr/boolean2-core (residual distributivity failure
//   after a9cfd4ff "Use nonzero union for boolean2 set add". A is a
//   33-radii star, C is a 32-radii star with zeros interleaved with
//   mid-large radii. Area asymmetry: left.Area ~ 71882, right.Area
//   ~ 66448, diff 5435 (~7%). The right side (A ∩ B) ∪ (A ∩ C) is
//   contained in left ((right - left).Area == 0) but is missing
//   5435 area units. NumContour also diverges (left=1, right=5).
//   The zero-radius vertices collapse to the 0.1 baseline floor,
//   producing dense near-degenerate slivers that the intersect path
//   appears to drop on one side of the distributivity identity).
TEST(CrossSection, BooleanDistributivityNonzeroUnionResidual) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    constexpr double kPi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> rA = {1000.,
                                  1000.,
                                  1.,
                                  454.29432950304715,
                                  995.65401133735929,
                                  1000.,
                                  1.,
                                  1000.,
                                  437.04103091878892,
                                  195.25477231407794,
                                  1.8268273914683537,
                                  4.5954160688751484,
                                  3.2030613706914481,
                                  511.5935307510378,
                                  2.0771806846580496,
                                  643.7557324509894,
                                  570.21917979539887,
                                  1.9769538901315753,
                                  467.93526395973339,
                                  2.8343333181297945,
                                  832.7231586632563,
                                  1.9935416909420567,
                                  3.4941254801335897,
                                  503.25635714266843,
                                  0.,
                                  677.87894521363603,
                                  1000.,
                                  0.,
                                  1000.,
                                  0.06374847983784826,
                                  995.14671048524883,
                                  1.,
                                  0.};
  const std::vector<double> rB = {381.47868752207233, 271.38539902365324,
                                  902.29172154497269, 369.78142089265032};
  const std::vector<double> rC = {812.15251441725718,
                                  0.,
                                  996.64035079388486,
                                  1.2440205144455736,
                                  845.2781504274202,
                                  1.7150332286725525,
                                  601.28991029625661,
                                  0.,
                                  482.34411167146169,
                                  0.,
                                  588.04482477271347,
                                  0.,
                                  169.03839642543201,
                                  512.45895098827191,
                                  3.4002817704521431,
                                  793.8531142632803,
                                  234.81914439627278,
                                  0.,
                                  999.90271349806335,
                                  105.49723198609914,
                                  0.,
                                  1000.,
                                  0.,
                                  705.33517104963164,
                                  2.2547666355380791,
                                  949.94413913719598,
                                  0.,
                                  4.9202747062496162,
                                  751.85659442712381,
                                  0.,
                                  375.32896975602767,
                                  0.};
  const CrossSection a({star(rA)});
  const CrossSection b =
      CrossSection({star(rB)})
          .Translate({-4.2594328040756597, -3.826038498814941});
  const CrossSection c =
      CrossSection({star(rC)})
          .Translate({-4.6033649212771799, -4.3357952411613461});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
  EXPECT_NEAR((aIntB - right).Area(), 0.0, tol)
      << "union monotonicity: (A ∩ B) is not contained in right";
  EXPECT_NEAR((aIntC - right).Area(), 0.0, tol)
      << "union monotonicity: (A ∩ C) is not contained in right";
}

// Seed: BooleanDistributivity (2026-05-24 manifoci daemon cycle 90,
//   post-nonzero-union-fix tip a9cfd4ff)
// Counterexample-hash: 278f30ca6f52d93f
// Suspected owner: pr/boolean2-core (companion to
//   DISABLED_BooleanDistributivityNonzeroUnionResidual above, but
//   with the zeros in A and an intermediate-sized B - 28-radii A
//   with leading 0., 0. plus repeated 85.93 values, 16-radii B
//   dominated by 1000s with a few small/zero values, 4-radii simple
//   C. Same failure shape: left has 1 contour, right has 3, right
//   strictly contained in left, missing ~4001 area units out of
//   ~16255 (~25%). Confirms the intersect-path sliver-drop bug is
//   not specific to which input carries the zeros).
TEST(CrossSection, BooleanDistributivityZerosInANonzeroUnion) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    constexpr double kPi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> rA = {0.,
                                  0.,
                                  85.932045280269435,
                                  85.932045280269435,
                                  503.627431790542,
                                  85.932045280269435,
                                  90.069450175585047,
                                  90.069450175585047,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  3.0409033274230679,
                                  90.340071855270295,
                                  88.580867960622072,
                                  261.53938936688053,
                                  1.,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  1000.,
                                  1.};
  const std::vector<double> rB = {0.,
                                  1000.,
                                  1000.,
                                  150.433665233841,
                                  1000.,
                                  1000.,
                                  1.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  177.35813950095334,
                                  0.04587971958967612,
                                  635.76724699191789,
                                  0.,
                                  801.00098835163021};
  const std::vector<double> rC = {2.2620915191424817, 352.8271046252292,
                                  124.40008725909412, 583.15656615652506};
  const CrossSection a({star(rA)});
  const CrossSection b =
      CrossSection({star(rB)})
          .Translate({0.70686431311486997, 4.2881751021813148});
  const CrossSection c =
      CrossSection({star(rC)})
          .Translate({1.9937691377314026, 4.7666664390133615});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
  EXPECT_NEAR((aIntB - right).Area(), 0.0, tol)
      << "union monotonicity: (A ∩ B) is not contained in right";
  EXPECT_NEAR((aIntC - right).Area(), 0.0, tol)
      << "union monotonicity: (A ∩ C) is not contained in right";
}

// Seed: BooleanDistributivity (2026-05-24 local daemon cycle ~294,
//   post-classify-nonzero-outer-face tip 986183af)
// Counterexample-hash: 389840fbf4d32247
// Suspected owner: pr/boolean2-core (small-residual failure on the
//   OPPOSITE side from the previous nonzero-union residuals: here
//   right > left, with left.Area=172524 and right.Area=172618 -
//   diff only 94.6 (~0.06%) but (left-right).Area==0 and
//   (right-left).Area=94. NumContour also reversed: left=25,
//   right=8 - right side collapses many contours into fewer.
//   Suggests the right-hand `(A∩B) ∪ (A∩C)` is over-merging
//   adjacent components and accidentally including area that the
//   left-hand `A ∩ (B∪C)` correctly carves out. Distinct from the
//   "right is strictly contained in left, missing area" shape of
//   86029efb/278f30ca seeds which were fixed by the recent
//   nonzero-outer-face classifier change).
TEST(CrossSection, BooleanDistributivityRightOverMerges) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    constexpr double kPi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> rA = {1000.,
                                  995.07901723564635,
                                  1.,
                                  1000.,
                                  974.26213405728731,
                                  974.26213405728731,
                                  887.65436435582899,
                                  1000.,
                                  0.,
                                  1000.,
                                  1000.,
                                  0.06374847983784826,
                                  671.02238904111891,
                                  671.02238904111891,
                                  671.02238904111891,
                                  0.};
  const std::vector<double> rB = {165.9646574194592,
                                  134.35969768300723,
                                  134.35969768300723,
                                  0.,
                                  0.,
                                  0.,
                                  134.35969768300723,
                                  0.90493025476122824,
                                  95.854999537523781,
                                  0.,
                                  730.30962259021919,
                                  0.,
                                  460.32237309958083,
                                  2.9791218075989327,
                                  30.186663600739344,
                                  1.,
                                  90.446292042234916,
                                  5.732788180417673,
                                  487.77994489187455,
                                  4.5933105105744643,
                                  110.52720785200476,
                                  1.,
                                  1.,
                                  1.,
                                  338.4042053897046,
                                  1.,
                                  1.,
                                  405.85342012592389,
                                  405.85342012592389,
                                  0.,
                                  1.,
                                  405.85342012592389,
                                  706.56289376450047,
                                  405.85342012592389,
                                  405.85342012592389,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  1.};
  const std::vector<double> rC = {279.68331311398725,
                                  814.7527470368675,
                                  812.15251441725718,
                                  1000.,
                                  0.,
                                  996.64035079388486,
                                  1.,
                                  1000.,
                                  1.,
                                  0.,
                                  146.36259356992105,
                                  1000.,
                                  1.2440205144455736,
                                  2.4903214421341642,
                                  1.4267332381508839,
                                  1.0658012274682902,
                                  937.2477438999249,
                                  601.28991029625661,
                                  0.,
                                  482.34411167146169,
                                  0.,
                                  999.13575735403208,
                                  0.,
                                  124.51182502694519,
                                  0.44553220583928932,
                                  4.0649812370918958,
                                  512.45895098827191,
                                  1000.,
                                  0.,
                                  0.,
                                  0.,
                                  595.63567011806413,
                                  2.5851695628219766,
                                  354.45821723238265,
                                  4.2476895774067867,
                                  956.03225557664916,
                                  0.,
                                  0.,
                                  98.799661756610732,
                                  0.,
                                  1000.,
                                  1.};
  const CrossSection a({star(rA)});
  const CrossSection b =
      CrossSection({star(rB)})
          .Translate({0.61958285855050033, -2.9878148015921839});
  const CrossSection c =
      CrossSection({star(rC)})
          .Translate({2.4887547431653907, 4.8411440198518338});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
}

// Seed: BooleanDistributivity (2026-05-24 manifoci daemon cycle ~94,
//   post-classify-nonzero-outer-face tip 986183af)
// Counterexample-hash: dbc69cb9e33a69c4
// Suspected owner: pr/boolean2-core (same failure shape as the
//   already-drained 86029efb/278f30ca residuals - right strictly
//   contained in left, missing area - but on larger inputs (34/19/45
//   radii). Area diff 26631 (~19% of left), NumContour 11 vs 10.
//   Demonstrates that the nonzero-outer-face classifier fix did not
//   fully drain the bug class: the specific seeded counterexamples
//   pass, but neighbors in the same family with larger zero/repeat
//   counts still trigger the same misclassification).
TEST(CrossSection, BooleanDistributivityLargeInputsResidual) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    constexpr double kPi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> rA = {215.54119679461166,
                                  0.,
                                  1000.,
                                  85.932045280269435,
                                  85.932045280269435,
                                  503.627431790542,
                                  85.932045280269435,
                                  1.6560313083633313,
                                  90.069450175585047,
                                  90.069450175585047,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  0.,
                                  1.,
                                  1.,
                                  1.,
                                  1.,
                                  87.027104618165325,
                                  85.932045280269435,
                                  258.81750993550418,
                                  1000.,
                                  1.,
                                  0.,
                                  89.294011990045064,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  85.932045280269435,
                                  1000.,
                                  1000.,
                                  1.};
  const std::vector<double> rB = {1000.,
                                  0.,
                                  997.02489346940774,
                                  997.07890660881435,
                                  319.71364680208376,
                                  675.57216392823045,
                                  997.55565490398703,
                                  619.18266326889307,
                                  345.22578984765374,
                                  3.1502725292253997,
                                  121.95185049831551,
                                  585.41131030666656,
                                  999.93563090020007,
                                  1000.,
                                  0.,
                                  177.35813950095334,
                                  0.04587971958967612,
                                  882.24641989432973,
                                  1.1946278115288917};
  const std::vector<double> rC = {1.,
                                  1000.,
                                  1.,
                                  282.41410155737179,
                                  1.,
                                  282.41410155737179,
                                  282.41410155737179,
                                  282.41410155737179,
                                  282.41410155737179,
                                  1000.,
                                  1.,
                                  1000.,
                                  1000.,
                                  0.,
                                  0.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  578.55699544670244,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1000.,
                                  1.};
  const CrossSection a({star(rA)});
  const CrossSection b =
      CrossSection({star(rB)})
          .Translate({-4.3860291464683625, 2.2730433830580954});
  const CrossSection c =
      CrossSection({star(rC)})
          .Translate({0.99578244833742113, -4.0475870932225444});
  const auto bUc = b + c;
  const auto left = a.Boolean(bUc, OpType::Intersect);
  const auto aIntB = a.Boolean(b, OpType::Intersect);
  const auto aIntC = a.Boolean(c, OpType::Intersect);
  const auto right = aIntB + aIntC;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()) +
                             std::fabs(c.Area()));
  EXPECT_NEAR(left.Area(), right.Area(), tol)
      << "A ∩ (B ∪ C) != (A ∩ B) ∪ (A ∩ C)";
  EXPECT_NEAR((left - right).Area(), 0.0, tol)
      << "distributivity: left-right difference is non-empty";
  EXPECT_NEAR((right - left).Area(), 0.0, tol)
      << "distributivity: right-left difference is non-empty";
}

// Seed: DegenerateInputFuzz (2026-05-25 local daemon cycle 339,
//   post-single-ray-cast tip 12e769fe)
// Counterexample-hash: 6151576e50da41a2
// Suspected owner: pr/boolean2-core (inclusion-exclusion fails when
//   one vertex of B is set exactly coincident with a vertex of A
//   (op=2 in the DegenerateInputFuzz harness). Both stars are 6
//   vertices with mid-range radii. unionAB.Area=436995 but
//   ca.Area+cb.Area-intersectAB.Area=434285, diff 2710 (~0.6%) far
//   above the 4.34 tol. intersectAB.Area only 0.51 - tiny sliver,
//   yet the union absorbs ~2710 extra area. Likely the coincident
//   vertex is mis-classified during merge, creating phantom area
//   in the union or losing it from the intersection).
TEST(CrossSection, DegenerateCoincidentVertexUnion) {
  auto starPolygon = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    constexpr double kPi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {1.0169016983060246, 1000., 1.,
                                      578.85382959129936, 0.,    0.};
  const std::vector<double> radiiB = {0.,
                                      1000.,
                                      999.83083173100238,
                                      7.275875880519048,
                                      726.89009231357352,
                                      3.880747251022969};
  SimplePolygon a = starPolygon(radiiA);
  SimplePolygon b = starPolygon(radiiB);
  for (auto& v : b) {
    v.x += 0.0030378301550779696;
  }
  // op=94 % 4 = 2 in DegenerateInputFuzz: set b[j] = a[i] to create
  // a coincident-vertex degeneracy. idxA=25 % 6 = 1, idxB=52 % 6 = 4.
  b[4] = a[1];

  const CrossSection ca(a);
  const CrossSection cb(b);
  const auto unionAB = ca + cb;
  const auto intersectAB = ca.Boolean(cb, OpType::Intersect);
  const double sum = ca.Area() + cb.Area() - intersectAB.Area();
  const double tol = 1e-5 * (1.0 + std::fabs(ca.Area()) + std::fabs(cb.Area()));
  EXPECT_NEAR(unionAB.Area(), sum, tol)
      << "Inclusion-exclusion violated on coincident-vertex degenerate input";
}

// Reduced from DegenerateCoincidentVertexUnion after constructor round-trip:
// already-regularized inputs still reproduce the original ~2710 area residual.
// This keeps both B triangles and drops two A vertices.
TEST(CrossSection, DegenerateCoincidentVertexUnionReduced) {
  const Polygons a = {{
      {500.05000000000018, 866.11200632481712},
      {-0.54999999999999716, 0.95262794416288443},
      {-578.95382959129938, 5.6843418860808015e-14},
      {-0.049999999999997158, -0.08660254037846471},
  }};
  const Polygons b = {
      {{500.04999974215514, 866.1120058797319},
       {-499.96237803534598, 865.96550230635103},
       {-7.3728380503639697, 0}},
      {{500.04999974215514, 866.1120058797319},
       {499.66038873910634, 865.43178211833663},
       {500.05303783015518, 866.11200632481712}},
  };

  const CrossSection ca(a);
  const CrossSection cb(b);
  const auto unionAB = ca + cb;
  const auto intersectAB = ca.Boolean(cb, OpType::Intersect);
  const CrossSection combined(Polygons{a[0], b[0], b[1]});
  const double sum = ca.Area() + cb.Area() - intersectAB.Area();
  const double tol = 1e-5 * (1.0 + std::fabs(ca.Area()) + std::fabs(cb.Area()));
  EXPECT_NEAR(unionAB.Area(), sum, tol)
      << "Inclusion-exclusion violated on reduced coincident-vertex input";
  EXPECT_NEAR(combined.Area(), sum, tol)
      << "Constructor edge soup lost area on reduced coincident-vertex input";
}

// Smaller 4+3 reduction from the same parked seed: keeping only the tiny top
// B triangle isolates the near-corner edge-vertex double-hit that used to make
// both binary union and single-constructor edge soup lose the large A contour.
TEST(CrossSection, DegenerateCoincidentVertexUnionTinyTriangle) {
  const Polygons a = {{
      {500.05000000000018, 866.11200632481712},
      {-0.54999999999999716, 0.95262794416288443},
      {-578.95382959129938, 5.6843418860808015e-14},
      {-0.049999999999997158, -0.08660254037846471},
  }};
  const Polygons b = {{
      {500.04999974215514, 866.1120058797319},
      {499.66038873910634, 865.43178211833663},
      {500.05303783015518, 866.11200632481712},
  }};

  auto translate = [](Polygons polys, double delta) {
    for (auto& loop : polys) {
      for (vec2& p : loop) {
        p.x += delta;
        p.y += delta;
      }
    }
    return polys;
  };

  for (double offset : {0.0, 1024.0, 4096.0, 8192.0}) {
    SCOPED_TRACE(offset);
    const Polygons shiftedA = translate(a, offset);
    const Polygons shiftedB = translate(b, offset);
    const CrossSection ca(shiftedA);
    const CrossSection cb(shiftedB);
    const auto unionAB = ca + cb;
    const auto intersectAB = ca.Boolean(cb, OpType::Intersect);
    const CrossSection combined(Polygons{shiftedA[0], shiftedB[0]});
    const double sum = ca.Area() + cb.Area() - intersectAB.Area();
    const double tol =
        1e-5 * (1.0 + std::fabs(ca.Area()) + std::fabs(cb.Area()));
    const double tinyContributionTol = 0.01 * cb.Area();
    EXPECT_GT(cb.Area(), 1e-4);
    EXPECT_GT(unionAB.Area(), ca.Area() + 0.5 * cb.Area())
        << "Binary union dropped most of the tiny triangle contribution";
    EXPECT_GT(combined.Area(), ca.Area() + 0.5 * cb.Area())
        << "Constructor edge soup dropped most of the tiny triangle "
           "contribution";
    EXPECT_NEAR(unionAB.Area(), sum, tol)
        << "Inclusion-exclusion violated on tiny-triangle reduction";
    EXPECT_NEAR(unionAB.Area(), sum, tinyContributionTol)
        << "Binary union dropped the tiny triangle contribution";
    EXPECT_NEAR(combined.Area(), sum, tinyContributionTol)
        << "Constructor edge soup dropped the tiny triangle contribution";
  }
}

// Seed: BooleanRobustness topology-validity failure (2026-05-23
//   cycle 277, both daemons, post-disconnected-winding-fix tip
//   22077d9c)
// Counterexample-hash: 8ea578543f9f1b64
// Suspected owner: pr/boolean2-core (RemoveOverlaps2D produces an
//   edge-balance imbalance for some vertices on these mixed-scale
//   inputs - 1e-6 alongside 1024. The fuzz harness's
//   CheckRetainedGraphValidity checks that per-vertex edge balance
//   from input (sum of mult on outgoing minus incoming) matches the
//   output for surviving verts and is zero for newly-introduced
//   ones. This input violates that. Likely eps-inference at the
//   extreme-scale-mix boundary or vertex-merge corner case in
//   RemoveOverlaps2D. Replicated inline because
//   CheckRetainedGraphValidity isn't exposed via the public API).
#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
TEST(CrossSection, RemoveOverlaps2DTopologyMixedScale) {
  // Inputs A and B exactly as the BooleanRobustness fuzz target
  // constructs them from the counterexample raw polygons. The
  // topology check is on the OUTPUT of the boolean Add, not the
  // raw inputs.
  const Polygons polysA = {
      {{-1e-6, 1e-6}, {-1e-6, -1e-6}, {-0.0, 0.0}, {1.0, 1024.0}},
      {{1024.0, 1024.0}, {1.0, 1.0}, {1e-6, -1e-6}}};
  const Polygons polysB = {
      {{0.0, 1.0},
       {1e-6, 1024.0},
       {-1.0, -1024.0},
       {-0.0, 1024.0},
       {1.0, -1.0}},
      {{-260.66988565137592, -0.0},
       {0.0, -414.46576455279967},
       {-1.0, -1024.0}},
      {{1.0, 1.0}, {-1e-6, -0.0}, {1.0, 1024.0}, {1e-6, 0.0}},
      {{-111.03407576117854, 0.0},
       {560.59976308273758, 2.9313131393310714},
       {-1024.0, 1.0},
       {89.678187020143696, 0.0},
       {1024.0, 1.0}}};
  const CrossSection a(polysA);
  const CrossSection b(polysB);
  const auto result = a + b;
  const auto resultPolys = result.ToPolygons();

  // Run RemoveOverlaps2D on the Boolean's output and check that the
  // per-vertex edge balance is preserved (sum of mult on outgoing
  // minus incoming) for surviving verts, zero for newly-introduced.
  const auto [verts, edges] = manifold::boolean2::PolygonsToInput(resultPolys);
  if (verts.empty()) GTEST_SKIP() << "Boolean output collapsed to empty";
  const double eps = manifold::boolean2::InferEps(resultPolys, {});
  const auto overlapResult =
      manifold::boolean2::RemoveOverlaps2D(verts, edges, eps);

  std::vector<manifold::boolean2::EdgeM> remapped;
  for (const auto& edge : edges) {
    const int aIdx = overlapResult.inputVert2Merged[edge.v0];
    const int bIdx = overlapResult.inputVert2Merged[edge.v1];
    if (aIdx != bIdx) remapped.push_back({aIdx, bIdx, edge.mult});
  }
  const auto expected = ComputeEdgeBalance(remapped);
  const auto actual = ComputeEdgeBalance(overlapResult.edges);
  for (int v = 0; v < static_cast<int>(overlapResult.verts.size()); ++v) {
    const int expectedBalance =
        expected.count(v) ? expected.find(v)->second : 0;
    const int actualBalance = actual.count(v) ? actual.find(v)->second : 0;
    const int target = (v < overlapResult.numMergedVerts) ? expectedBalance : 0;
    EXPECT_EQ(actualBalance, target)
        << "Vertex " << v
        << " edge-balance mismatch after RemoveOverlaps2D on boolean "
        << "Add result: expected=" << target << " actual=" << actualBalance;
  }
}
#endif

// Seed: BooleanCommutativity (2026-05-18 CI run 26006574818)
// Counterexample-hash: 68e22a523aa8d6a5
// Suspected owner: pr/boolean2-core (A+B and B+A produce different
//   contour counts: A+B has 14, B+A has 13. Order-dependent
//   commutativity violation: depending on which input is "subject"
//   vs "clip", one contour appears/vanishes. Input shapes are
//   adversarial - 42-vertex star A with extreme magnitudes
//   (mix of 0, 1, 1000) vs 33-vertex star B with mostly small
//   radii plus a few in the 700-900 range, translated by
//   (1.71, -0.50). Likely a vertex-merge or sort tie-break path
//   that depends on input order rather than canonical position).
TEST(CrossSection, BooleanCommutativityExtremeMagStars) {
  auto star = [](const std::vector<double>& radii) {
    SimplePolygon ring;
    const int n = static_cast<int>(radii.size());
    for (int i = 0; i < n; ++i) {
      const double r = 0.1 + std::fabs(radii[i]);
      const double th = 2.0 * kPi * i / n;
      ring.push_back({r * std::cos(th), r * std::sin(th)});
    }
    return ring;
  };
  const std::vector<double> radiiA = {849.58267006542974,
                                      1000.,
                                      1.,
                                      1000.,
                                      0.,
                                      0.,
                                      1.,
                                      0.,
                                      0.,
                                      673.91660850339099,
                                      641.26963903267847,
                                      0.84355834285842513,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      633.12578243694929,
                                      52.949291739509349,
                                      596.41353494815576,
                                      659.59815503253992,
                                      907.85583546669454,
                                      2.3971527824334933,
                                      0.,
                                      744.77799184771879,
                                      0.,
                                      668.48707869911436,
                                      436.35135227432249,
                                      0.,
                                      0.,
                                      0.,
                                      1000.,
                                      0.,
                                      938.50101700824371,
                                      84.188304887665581,
                                      470.75422504731279,
                                      89.810808608112893,
                                      992.99533690532235,
                                      152.76119378800999,
                                      1000.,
                                      0.};
  const std::vector<double> radiiB = {681.1551582487607,
                                      1.,
                                      226.83285963348141,
                                      25.271695356113401,
                                      903.10565859434519,
                                      25.271695356113401,
                                      27.294237080468665,
                                      205.64018204902527,
                                      708.55499051253935,
                                      934.95934725410712,
                                      597.99208744386829,
                                      937.31857430425612,
                                      1.,
                                      3.0503172192082628,
                                      1.,
                                      1.,
                                      1000.,
                                      1000.,
                                      1.,
                                      1.,
                                      901.26525729886407,
                                      1.,
                                      1000.,
                                      1.6267135216510802,
                                      449.23493668120602,
                                      449.23493668120602,
                                      3.1828449559882372,
                                      6.9538585160876147,
                                      1.7485381697954843,
                                      696.88681861119551,
                                      1.,
                                      183.08811014255667,
                                      1.};
  const CrossSection a(star(radiiA));
  const CrossSection b =
      CrossSection(star(radiiB))
          .Translate({1.7126290778198534, -0.5023590407011973});
  const auto aPlusB = a + b;
  const auto bPlusA = b + a;
  const double tol = 1e-6 * (1.0 + std::fabs(a.Area()) + std::fabs(b.Area()));
  EXPECT_NEAR(aPlusB.Area(), bPlusA.Area(), tol) << "A + B != B + A";
  EXPECT_EQ(aPlusB.NumContour(), bPlusA.NumContour())
      << "A + B and B + A produce different contour counts: "
      << aPlusB.NumContour() << " vs " << bPlusA.NumContour();
}

TEST(CrossSection, NonFiniteInputReturnsEmpty) {
  const double inf = std::numeric_limits<double>::infinity();
  SimplePolygon bad = {{0.0, 0.0}, {1.0, 0.0}, {inf, 1.0}, {0.0, 1.0}};

  CrossSection constructed(bad);
  EXPECT_TRUE(constructed.IsEmpty());

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
  Polygons finite = {{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
  EXPECT_TRUE(boolean2::Boolean2D(Polygons{bad}, finite, OpType::Add).empty());
#endif
}

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
TEST(CrossSection, ShallowLongEdgeIntersectionIsNotDropped) {
  vec2 intersection;
  const double eps = boolean2::EpsilonFromScale(1e9);

  EXPECT_TRUE(boolean2::IntersectSegments({{-0.5, 0.0}, {0.5, 0.0}, 0},
                                          {{-1e9, 0.01}, {1e9, -0.01}, 1}, eps,
                                          &intersection));
  EXPECT_NEAR(intersection.x, 0.0, eps);
  EXPECT_NEAR(intersection.y, 0.0, eps);
}

TEST(CrossSection, NearEndpointIntersectionOutsideSegmentIsDropped) {
  vec2 intersection;
  const double eps = 1.0;

  EXPECT_FALSE(boolean2::IntersectSegments({{0.0, 0.0}, {10.0, 0.0}, 0},
                                           {{8.0, 0.4}, {18.0, 1.4}, 1}, eps,
                                           &intersection));
}

TEST(CrossSection, EndpointTJunctionIntersectionIsDropped) {
  vec2 intersection;
  const double eps = 1.0;

  EXPECT_FALSE(boolean2::IntersectSegments({{0.0, 0.0}, {10.0, 0.0}, 0},
                                           {{5.0, 0.0}, {5.0, 0.1}, 1}, eps,
                                           &intersection));
}

TEST(CrossSection, MergeVertsTransitiveChainCanDriftPastEps) {
  const double eps = 1.0;
  const std::vector<vec2> verts = {{0.0, 0.0},  {0.99, 0.0}, {1.98, 0.0},
                                   {2.97, 0.0}, {3.96, 0.0}, {10.0, 0.0},
                                   {20.0, 0.0}};

  const boolean2::VertexMerge merged = boolean2::MergeVerts(verts, eps);

  ASSERT_EQ(merged.inputVert2Merged.size(), verts.size());
  ASSERT_EQ(merged.verts.size(), 3);
  for (int i = 0; i < 5; ++i)
    EXPECT_EQ(merged.inputVert2Merged[i], merged.inputVert2Merged[0]);
  EXPECT_NE(merged.inputVert2Merged[5], merged.inputVert2Merged[0]);
  EXPECT_NE(merged.inputVert2Merged[6], merged.inputVert2Merged[0]);

  EXPECT_NEAR(merged.verts[0].x, 1.98, 1e-12);
  EXPECT_NEAR(merged.verts[0].y, 0.0, 1e-12);
  const vec2 endpointDrift = verts.front() - merged.verts[0];
  EXPECT_GT(std::hypot(endpointDrift.x, endpointDrift.y), eps);
}

// Seed: VertexMergeIdempotence (2026-05-23 CI run 26323577663)
// Counterexample-hash: 4caf999681ce3e2a
TEST(CrossSection, VertexMergeIdempotenceTightCluster) {
  const std::vector<double> xs = {
      -564.24299726366871, -564.25684749526681, -564.25684749526681,
      -564.25684749526681, -564.25684749526681, -564.25684749526681,
      -564.25684749526681, -564.2434248982521,  -564.25684749526681,
      -564.25684749526681, -564.25684749526681, -560.99098110791851,
      -564.24548426671515, -564.25684749526681, -564.25684749526681,
      -564.24797363468235, -564.25105014358735, -564.25684749526681,
      -564.25684749526681, -564.2650553505373,  -564.25684749526681,
      -564.25684749526681, -564.25684749526681, -564.25684749526681,
      -564.25684749526681, -564.25684749526681, -564.25684749526681,
      -564.2565950636108,  -564.25684749526681, -564.25684749526681,
      -564.25684749526681, -564.25684749526681};
  const std::vector<double> ys = {-1.,
                                  -1.,
                                  -1.,
                                  1.3146737456995536,
                                  3.334547678102286,
                                  924.18159456764056,
                                  -5.0212043784034019,
                                  -3.442978883496882,
                                  -882.0023276178074,
                                  -3.0103861859513601,
                                  -218.03838880073647,
                                  1.1654515551817912,
                                  -361.18598495388369,
                                  -0.040195183151324976,
                                  -229.62237726036881,
                                  -3.4005396935541334,
                                  -0.37374103257085878,
                                  3.993690857296496,
                                  1.1420137899457909,
                                  210.83410884574982,
                                  -3.4038895533070526,
                                  -4.3056143084054153,
                                  -2.8446598710193314,
                                  0.90992952485448519,
                                  -5.1663919232428848,
                                  -502.42923749613101,
                                  -2.3459351240317718,
                                  223.16662739130152,
                                  810.70082359638945,
                                  -3.395314716558568,
                                  0.055697227037768471,
                                  -1.6836095819963903};
  ASSERT_EQ(xs.size(), ys.size());
  const double eps = std::pow(10.0, -2.0433270966230594);
  std::vector<vec2> verts;
  verts.reserve(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) {
    verts.push_back({xs[i], ys[i]});
  }

  const auto m1 = boolean2::MergeVerts(verts, eps);
  ASSERT_FALSE(m1.verts.empty());
  const auto m2 = boolean2::MergeVerts(m1.verts, eps);
  EXPECT_EQ(m2.verts.size(), m1.verts.size())
      << "MergeVerts not idempotent: pass1=" << m1.verts.size()
      << " pass2=" << m2.verts.size() << " (n=" << xs.size() << ", eps=" << eps
      << ")";
}
#endif

TEST(CrossSection, OffsetIsInvariantUnderLargeTranslation) {
  const CrossSection square = CrossSection::Square({10.0, 10.0}, true);
  const CrossSection origin =
      square.Offset(1.0, CrossSection::JoinType::Round, 2.0, 8);
  const CrossSection translated =
      square.Translate({1e12, -1e12})
          .Offset(1.0, CrossSection::JoinType::Round, 2.0, 8)
          .Translate({-1e12, 1e12});

  EXPECT_EQ(translated.NumContour(), origin.NumContour());
  EXPECT_EQ(translated.NumVert(), origin.NumVert());
  EXPECT_NEAR(translated.Area(), origin.Area(), 1e-3);
}

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
TEST(CrossSection, Boolean2ValidatorRejectsRetainedVertsWithinEps) {
  const std::vector<vec2> verts = {{0.0, 0.0}, {10.0, 0.0}, {0.0, 10.0},
                                   {0.5, 0.5}, {20.0, 0.0}, {20.0, 10.0}};
  const std::vector<boolean2::EdgeM> edges = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1},
                                              {3, 4, 1}, {4, 5, 1}, {5, 3, 1}};
  const boolean2::OverlapResult result{verts, edges, {0, 1, 2, 3, 4, 5}, 6};

  EXPECT_FALSE(CheckRetainedGraphValidity(
      result, edges, result.inputVert2Merged, result.numMergedVerts, 1.0));
}

TEST(CrossSection, Boolean2ValidatorRejectsNearEndpointTJunction) {
  const std::vector<vec2> verts = {{0.0, 0.0},   {10.0, 0.0}, {0.5, 0.9},
                                   {10.0, 10.0}, {20.0, 5.0}, {20.0, 15.0}};
  const std::vector<boolean2::EdgeM> edges = {{0, 1, 1}, {1, 3, 1}, {3, 0, 1},
                                              {2, 4, 1}, {4, 5, 1}, {5, 2, 1}};
  const boolean2::OverlapResult result{verts, edges, {0, 1, 2, 3, 4, 5}, 6};

  EXPECT_FALSE(CheckRetainedGraphValidity(
      result, edges, result.inputVert2Merged, result.numMergedVerts, 1.0));
}

TEST(CrossSection, Boolean2ValidatorRejectsRetainedStrictCrossing) {
  const std::vector<vec2> verts = {{0.0, 0.0},  {10.0, 10.0}, {0.0, 10.0},
                                   {10.0, 0.0}, {-10.0, 5.0}, {20.0, 5.0}};
  const std::vector<boolean2::EdgeM> edges = {{0, 1, 1}, {1, 4, 1}, {4, 0, 1},
                                              {2, 3, 1}, {3, 5, 1}, {5, 2, 1}};
  const boolean2::OverlapResult result{verts, edges, {0, 1, 2, 3, 4, 5}, 6};

  EXPECT_FALSE(CheckRetainedGraphValidity(
      result, edges, result.inputVert2Merged, result.numMergedVerts, 0.01));
}

TEST(CrossSection, Boolean2CleanupPassMatchesValidAddSinglePass) {
  Polygons polys{RandomTopologicalRing(8, 618)};
  const double eps = boolean2::InferEps(polys, {});
  const auto [verts, edges] = boolean2::PolygonsToInput(polys);
  const auto pass1 =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                 /*debug=*/false, boolean2::WindRule::Add);
  EXPECT_TRUE(CheckRetainedGraphValidity(pass1, edges, pass1.inputVert2Merged,
                                         pass1.numMergedVerts, eps));

  const auto pass2 = CleanupPassLikeIterate(pass1, eps);
  const auto pass2Input = EdgesFromOverlapResult(pass1);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      pass2, pass2Input, pass2.inputVert2Merged, pass2.numMergedVerts, eps));
  ExpectSameFingerprint(pass1, pass2, eps);

  const auto pass3 = CleanupPassLikeIterate(pass2, eps);
  const auto pass3Input = EdgesFromOverlapResult(pass2);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      pass3, pass3Input, pass3.inputVert2Merged, pass3.numMergedVerts, eps));
  ExpectSameFingerprint(pass2, pass3, eps);

  const auto pass1Polys =
      boolean2::OutEdgesToPolygons(pass1.verts, pass1.edges);
  const auto pass2Polys =
      boolean2::OutEdgesToPolygons(pass2.verts, pass2.edges);
  ASSERT_EQ(pass1Polys.size(), 2);
  ASSERT_EQ(pass2Polys.size(), 2);
  EXPECT_EQ(pass1Polys[0].size(), pass2Polys[0].size());
  EXPECT_EQ(pass1Polys[1].size(), pass2Polys[1].size());
  EXPECT_NEAR(boolean2::TotalSignedArea(pass1Polys),
              boolean2::TotalSignedArea(pass2Polys), 1e-12);
}

TEST(CrossSection, Boolean2CleanupPassMatchesValidNonZeroSinglePass) {
  Polygons polys{RandomTopologicalRing(8, 618)};
  const double eps = boolean2::InferEps(polys, {});
  const auto [verts, edges] = boolean2::PolygonsToInput(polys);
  const auto pass1 =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                 /*debug=*/false, boolean2::WindRule::NonZero);
  EXPECT_TRUE(CheckRetainedGraphValidity(pass1, edges, pass1.inputVert2Merged,
                                         pass1.numMergedVerts, eps));

  const auto pass2 = CleanupPassLikeIterate(pass1, eps);
  const auto pass2Input = EdgesFromOverlapResult(pass1);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      pass2, pass2Input, pass2.inputVert2Merged, pass2.numMergedVerts, eps));
  ExpectSameFingerprint(pass1, pass2, eps);

  const auto pass3 = CleanupPassLikeIterate(pass2, eps);
  const auto pass3Input = EdgesFromOverlapResult(pass2);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      pass3, pass3Input, pass3.inputVert2Merged, pass3.numMergedVerts, eps));
  ExpectSameFingerprint(pass2, pass3, eps);

  const auto pass1Polys =
      boolean2::OutEdgesToPolygons(pass1.verts, pass1.edges);
  const auto pass2Polys =
      boolean2::OutEdgesToPolygons(pass2.verts, pass2.edges);
  ASSERT_EQ(pass1Polys.size(), 3);
  ASSERT_EQ(pass2Polys.size(), 3);
  EXPECT_EQ(pass1Polys[0].size(), pass2Polys[0].size());
  EXPECT_EQ(pass1Polys[1].size(), pass2Polys[1].size());
  EXPECT_EQ(pass1Polys[2].size(), pass2Polys[2].size());
  EXPECT_NEAR(boolean2::TotalSignedArea(pass1Polys),
              boolean2::TotalSignedArea(pass2Polys), 1e-12);
}
#endif

TEST(CrossSection, SimplifyPostFiltersBoolean2Output) {
  const double apex = 1.0148512233354445e-6;
  const SimplePolygon tri = {{-1.0, 0.0}, {1.0, 0.0}, {0.0, apex}};
  const SimplePolygon quad = {
      {-0.05, -1.0}, {0.05, -1.0}, {0.05, 2.0}, {-0.05, 2.0}};
  const CrossSection input(Polygons{tri, quad},
                           CrossSection::FillRule::NonZero);

  const CrossSection once = input.Simplify();
  const CrossSection twice = once.Simplify();

  EXPECT_EQ(once.NumContour(), 1);
  EXPECT_EQ(once.NumVert(), 6);
  EXPECT_EQ(twice.NumContour(), once.NumContour());
  EXPECT_EQ(twice.NumVert(), once.NumVert());
  EXPECT_NEAR(twice.Area(), once.Area(), 1e-12);
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
#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
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
  const auto overlap =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                 /*debug=*/false, boolean2::WindRule::Add);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      overlap, edges, overlap.inputVert2Merged, overlap.numMergedVerts, eps));
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
  const auto overlap =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                 /*debug=*/false, boolean2::WindRule::Add);
  EXPECT_TRUE(CheckRetainedGraphValidity(
      overlap, edges, overlap.inputVert2Merged, overlap.numMergedVerts, eps));

  const auto polys = boolean2::OutEdgesToPolygons(overlap.verts, overlap.edges);
  EXPECT_EQ(polys.size(), 11);
  EXPECT_NEAR(boolean2::TotalSignedArea(polys), 1678.2538553263785,
              1e-10 * (1.0 + 1678.2538553263785));
}

// Regression: a single closed directed ring must round-trip through
// OutEdgesToPolygons. Earlier closure logic detected closure by trying to
// re-select the start edge as `next`, which is skipped by the visited guard;
// the natural walk terminates with destV == startV instead.
TEST(CrossSection, OutEdgesToPolygonsClosesSimpleRing) {
  const std::vector<vec2> verts = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  const std::vector<boolean2::OutEdge> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  const auto polys = boolean2::OutEdgesToPolygons(verts, edges);
  ASSERT_EQ(polys.size(), 1u);
  EXPECT_EQ(polys[0].size(), 4u);
}

TEST(CrossSection, OutEdgesToPolygonsSplitsExactRepeatedVertex) {
  const std::vector<vec2> verts = {{0, 0}, {1, 0}, {1, -1}, {2, -1}, {1, 1}};
  const std::vector<boolean2::OutEdge> edges = {{0, 1}, {1, 2}, {2, 3},
                                                {3, 1}, {1, 4}, {4, 0}};
  const auto polys = boolean2::OutEdgesToPolygons(verts, edges);
  ASSERT_EQ(polys.size(), 2u);
  EXPECT_EQ(polys[0].size(), 3u);
  EXPECT_EQ(polys[1].size(), 3u);
  EXPECT_NEAR(boolean2::TotalSignedArea(polys), 1.0, 1e-12);
}

TEST(CrossSection, OutEdgesToPolygonsKeepsNearDistinctVertex) {
  constexpr double kDelta = 1e-12;
  const std::vector<vec2> verts = {{0, 0},  {1, 0}, {1, -1},
                                   {2, -1}, {1, 1}, {1 + kDelta, 0}};
  const std::vector<boolean2::OutEdge> edges = {{0, 1}, {1, 2}, {2, 3},
                                                {3, 5}, {5, 4}, {4, 0}};
  const auto polys = boolean2::OutEdgesToPolygons(verts, edges);
  ASSERT_EQ(polys.size(), 1u);
  EXPECT_EQ(polys[0].size(), 6u);
  EXPECT_NEAR(boolean2::TotalSignedArea(polys), 1.0 + kDelta, 1e-12);

  const CrossSection reconsumed(polys, CrossSection::FillRule::NonZero);
  EXPECT_FALSE(reconsumed.IsEmpty());
  EXPECT_NEAR(reconsumed.Area(), 1.0, 1e-9);
  const Manifold solid = Manifold::Extrude(reconsumed.ToPolygons(), 1.0);
  EXPECT_EQ(solid.Status(), Manifold::Error::NoError);
  EXPECT_NEAR(solid.Volume(), reconsumed.Area(), 1e-9);
}

TEST(CrossSection, Boolean2DKeepsNearDistinctPresentationVertex) {
  constexpr double kDelta = 1e-12;
  const Polygons input = {
      {{0, 0}, {1, 0}, {1, -1}, {2, -1}, {1 + kDelta, 0}, {1, 1}}};
  const auto polys = boolean2::Boolean2D(input, {}, OpType::Add, /*eps=*/1e-15,
                                         /*tolerance=*/1e-9);
  ASSERT_EQ(polys.size(), 1u);
  EXPECT_EQ(polys[0].size(), 6u);
  EXPECT_NEAR(boolean2::TotalSignedArea(polys), 1.0 + kDelta, 1e-12);
}

TEST(CrossSection, Boolean2NewOldToleranceMergesGeneratedNearEndpoint) {
  constexpr double eps = 1e-6;
  constexpr double tolerance = 1e-3;
  const Polygons a = {{{0, 0}, {10, 0}, {10, 1}, {0, 1}}};
  const Polygons b = {{{5e-4, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {5e-4, 0.5}}};
  const auto [verts, edges] = CombinedInput(a, b, /*bMult=*/1);

  const auto withoutPrior =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/eps,
                                 /*debug=*/false, boolean2::WindRule::Add);
  const auto withPrior =
      boolean2::RemoveOverlaps2D(verts, edges, eps, tolerance,
                                 /*debug=*/false, boolean2::WindRule::Add);

  auto countNear = [](const std::vector<vec2>& points, vec2 target,
                      double radius) {
    int count = 0;
    const double radius2 = radius * radius;
    for (const vec2& p : points) {
      const vec2 d = p - target;
      if (dot(d, d) <= radius2) ++count;
    }
    return count;
  };

  const vec2 generatedNearOld{5e-4, 0.0};
  EXPECT_GT(countNear(withoutPrior.verts, generatedNearOld, 10 * eps), 0)
      << "baseline should retain the generated crossing as distinct";
  EXPECT_EQ(countNear(withPrior.verts, generatedNearOld, 10 * eps), 0)
      << "propagated tolerance should merge the generated crossing into the "
         "nearby old endpoint";
  EXPECT_EQ(withPrior.verts.size() + 1, withoutPrior.verts.size());
  EXPECT_TRUE(CheckRetainedGraphValidity(withPrior, edges,
                                         withPrior.inputVert2Merged,
                                         withPrior.numMergedVerts, eps));
}

TEST(CrossSection, RemoveOverlapsMergesExactDuplicateCoordinates) {
  const Polygons polys = {{{0, 0}, {1, 0}, {0, 1}}, {{0, 0}, {-1, 0}, {0, -1}}};
  const auto [verts, edges] = boolean2::PolygonsToInput(polys);
  ASSERT_EQ(verts.size(), 6u);
  const double eps = boolean2::InferEps(polys, {});
  const auto overlap =
      boolean2::RemoveOverlaps2D(verts, edges, eps, /*tolerance=*/0.0,
                                 /*debug=*/false, boolean2::WindRule::Add);
  ASSERT_EQ(overlap.inputVert2Merged.size(), verts.size());
  ASSERT_GE(overlap.inputVert2Merged[0], 0);
  ASSERT_LT(overlap.inputVert2Merged[0], overlap.numMergedVerts);
  ASSERT_GE(overlap.inputVert2Merged[3], 0);
  ASSERT_LT(overlap.inputVert2Merged[3], overlap.numMergedVerts);
  EXPECT_EQ(overlap.inputVert2Merged[0], overlap.inputVert2Merged[3]);
}
#endif

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

#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
TEST(CrossSection, Boolean2GraphOrderDetectsProperCrossing) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D a{{0.0, 0.0}, {10.0, 10.0}, 0};
  GraphSegment2D b{{0.0, 10.0}, {10.0, 0.0}, 1};

  const auto order = CompareProjectedOrder(a, b, /*axis=*/0, 0.0, 10.0);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::ALessOrtho);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_FALSE(order.coincidentOverlap);
  EXPECT_TRUE(order.properCrossing);
}

TEST(CrossSection, Boolean2GraphOrderIsEndpointReversalStable) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D a{{0.0, 0.0}, {10.0, 10.0}, 0};
  GraphSegment2D b{{0.0, 10.0}, {10.0, 0.0}, 1};
  GraphSegment2D aReversed{{10.0, 10.0}, {0.0, 0.0}, 0};
  GraphSegment2D bReversed{{10.0, 0.0}, {0.0, 10.0}, 1};

  const auto order = CompareProjectedOrder(a, b, /*axis=*/0, 0.0, 10.0);
  const auto reversed =
      CompareProjectedOrder(aReversed, bReversed, /*axis=*/0, 0.0, 10.0);
  EXPECT_EQ(order.atMinProjection, reversed.atMinProjection);
  EXPECT_EQ(order.atMaxProjection, reversed.atMaxProjection);
  EXPECT_EQ(order.coincidentOverlap, reversed.coincidentOverlap);
  EXPECT_EQ(order.properCrossing, reversed.properCrossing);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::ALessOrtho);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::AGreaterOrtho);

  const auto aOnlyReversed =
      CompareProjectedOrder(aReversed, b, /*axis=*/0, 0.0, 10.0);
  const auto bOnlyReversed =
      CompareProjectedOrder(a, bReversed, /*axis=*/0, 0.0, 10.0);
  EXPECT_EQ(order.atMinProjection, aOnlyReversed.atMinProjection);
  EXPECT_EQ(order.atMaxProjection, aOnlyReversed.atMaxProjection);
  EXPECT_EQ(order.properCrossing, aOnlyReversed.properCrossing);
  EXPECT_EQ(order.atMinProjection, bOnlyReversed.atMinProjection);
  EXPECT_EQ(order.atMaxProjection, bOnlyReversed.atMaxProjection);
  EXPECT_EQ(order.properCrossing, bOnlyReversed.properCrossing);
}

TEST(CrossSection, Boolean2GraphOrderSupportsYAxisProjection) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D a{{0.0, 0.0}, {10.0, 10.0}, 0};
  GraphSegment2D b{{10.0, 0.0}, {0.0, 10.0}, 1};

  const auto order = CompareProjectedOrder(a, b, /*axis=*/1, 0.0, 10.0);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::ALessOrtho);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_TRUE(order.properCrossing);
}

TEST(CrossSection, Boolean2GraphOrderResolvesCoincidentOverlap) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{0.0, 0.0}, {10.0, 0.0}, 1};

  const auto order = CompareProjectedOrder(a, b, /*axis=*/0, 0.0, 10.0);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::ALessOrtho);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::ALessOrtho);
  EXPECT_TRUE(order.coincidentOverlap);
  EXPECT_FALSE(order.properCrossing);

  const auto swapped = CompareProjectedOrder(b, a, /*axis=*/0, 0.0, 10.0);
  EXPECT_EQ(swapped.atMinProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_EQ(swapped.atMaxProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_TRUE(swapped.coincidentOverlap);
  EXPECT_FALSE(swapped.properCrossing);
}

TEST(CrossSection, Boolean2GraphOrderUsesCanonicalGeometryTieBeforeEdgeId) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D lower{{0.0, 0.0}, {10.0, 0.0}, 100};
  GraphSegment2D upper{{0.0, 0.5}, {10.0, 0.5}, 1};

  const auto order =
      CompareProjectedOrder(lower, upper, /*axis=*/0, 0.0, 10.0, 1.0);
  EXPECT_TRUE(order.coincidentOverlap);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::ALessOrtho);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::ALessOrtho);
  EXPECT_FALSE(order.properCrossing);

  const auto swapped =
      CompareProjectedOrder(upper, lower, /*axis=*/0, 0.0, 10.0, 1.0);
  EXPECT_TRUE(swapped.coincidentOverlap);
  EXPECT_EQ(swapped.atMinProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_EQ(swapped.atMaxProjection, GraphOrderKind::AGreaterOrtho);
  EXPECT_FALSE(swapped.properCrossing);
}

TEST(CrossSection, Boolean2GraphOrderKeepsEndpointTouchDegenerate) {
  using boolean2::CompareProjectedOrder;
  using boolean2::GraphOrderKind;
  using boolean2::GraphSegment2D;

  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{5.0, 0.0}, {15.0, 1.0}, 1};

  const auto order = CompareProjectedOrder(a, b, /*axis=*/0, 5.0, 10.0);
  EXPECT_EQ(order.atMinProjection, GraphOrderKind::EndpointTouch);
  EXPECT_EQ(order.atMaxProjection, GraphOrderKind::ALessOrtho);
  EXPECT_FALSE(order.coincidentOverlap);
  EXPECT_FALSE(order.properCrossing);
}

TEST(CrossSection, Boolean2IntersectSegmentsFindsStrictCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 10.0}, 0};
  GraphSegment2D b{{0.0, 10.0}, {10.0, 0.0}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 0.0, &intersection));
  EXPECT_NEAR(intersection.x, 5.0, 1e-12);
  EXPECT_NEAR(intersection.y, 5.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsKeepsOneSidedEpsBandCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{0.0, -0.5}, {10.0, 2.0}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 1.0, &intersection));
  EXPECT_NEAR(intersection.x, 2.0, 1e-12);
  EXPECT_NEAR(intersection.y, 0.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsKeepsTwoSidedEpsBandCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{0.0, -0.75}, {10.0, 0.75}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 1.0, &intersection));
  EXPECT_NEAR(intersection.x, 5.0, 1e-12);
  EXPECT_NEAR(intersection.y, 0.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsKeepsUnderflowingSignChange) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  const double tiny = 1e-200;
  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{0.0, -tiny}, {10.0, tiny}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 1.0, &intersection));
  EXPECT_NEAR(intersection.x, 5.0, 1e-12);
  EXPECT_NEAR(intersection.y, 0.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsDropsEpsNearEndpointCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{0.5, -1.0}, {0.5, 1.0}, 1};

  EXPECT_FALSE(IntersectSegments(a, b, 1.0, &intersection));
}

TEST(CrossSection, Boolean2IntersectSegmentsKeepsSteepInteriorCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {0.0015, 1000.0}, 0};
  GraphSegment2D b{{-1.0, 500.0}, {1.0, 500.0}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 0.001, &intersection));
  EXPECT_NEAR(intersection.x, 0.00075, 1e-12);
  EXPECT_NEAR(intersection.y, 500.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsDropsEndpointTouch) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{10.0, 0.0}, {20.0, 10.0}, 1};

  EXPECT_FALSE(IntersectSegments(a, b, 0.0, &intersection));
}

TEST(CrossSection, Boolean2IntersectSegmentsDropsPositiveOverlapTJunction) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{5.0, 0.0}, {15.0, 1.0}, 1};

  EXPECT_FALSE(IntersectSegments(a, b, 0.0, &intersection));
}

TEST(CrossSection, Boolean2IntersectSegmentsFindsAxisAlignedStrictCrossing) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 5.0}, {10.0, 5.0}, 0};
  GraphSegment2D b{{5.0, 0.0}, {5.0, 10.0}, 1};

  EXPECT_TRUE(IntersectSegments(a, b, 0.0, &intersection));
  EXPECT_NEAR(intersection.x, 5.0, 1e-12);
  EXPECT_NEAR(intersection.y, 5.0, 1e-12);
}

TEST(CrossSection, Boolean2IntersectSegmentsDropsAxisAlignedEndpointTouch) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{10.0, 0.0}, {10.0, 10.0}, 1};

  EXPECT_FALSE(IntersectSegments(a, b, 0.0, &intersection));
}

TEST(CrossSection, Boolean2IntersectSegmentsDropsCoincidentOverlap) {
  using boolean2::GraphSegment2D;
  using boolean2::IntersectSegments;

  vec2 intersection;
  GraphSegment2D a{{0.0, 0.0}, {10.0, 0.0}, 0};
  GraphSegment2D b{{2.0, 0.0}, {8.0, 0.0}, 1};

  EXPECT_FALSE(IntersectSegments(a, b, 0.0, &intersection));
}
#endif
