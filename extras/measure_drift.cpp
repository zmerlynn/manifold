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
//
// Reports max cross-input vertex-pair distance among non-collapsed verts in
// units of eps, across a corpus of polygon pairs. Mirrors BinaryOpByRule's
// pipeline.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

#include "cross_section/boolean2/boolean2.h"
#include "cross_section/boolean2/vertex_merge.h"
#include "manifold/cross_section.h"

namespace mb = manifold::boolean2;

namespace {

using manifold::CrossSection;
using manifold::Polygons;
using manifold::Rect;
using manifold::SimplePolygon;
using manifold::vec2;

// Cutoff: pairs farther than this are unrelated geometry, not drift.
constexpr double kNearCoincidenceMaxRatio = 1e5;

struct DriftSample {
  std::string name;
  double maxRatio = 0.0;  // largest dist/eps among near-coincidence pairs
  double eps = 0.0;
  int nearPairs = 0;  // pairs with dist/eps < kNearCoincidenceMaxRatio
};

Polygons Translate(const Polygons& src, vec2 d) {
  Polygons out = src;
  for (auto& loop : out)
    for (vec2& v : loop) v = v + d;
  return out;
}

vec2 AabbCenter(const Polygons& a, const Polygons& b) {
  Rect box;
  for (const auto* src : {&a, &b})
    for (const auto& loop : *src)
      for (const vec2& v : loop) box.Union(v);
  return box.IsEmpty() ? vec2(0.0) : box.Center();
}

DriftSample Measure(const std::string& name, const CrossSection& a,
                    const CrossSection& b) {
  DriftSample s;
  s.name = name;
  const Polygons aPolys = a.ToPolygons();
  const Polygons bPolys = b.ToPolygons();
  if (aPolys.empty() || bPolys.empty()) return s;

  const vec2 origin = AabbCenter(aPolys, bPolys);
  const Polygons localA = Translate(aPolys, -origin);
  const Polygons localB = Translate(bPolys, -origin);
  s.eps = mb::InferEps(localA, localB);
  if (s.eps <= 0.0) return s;

  std::vector<vec2> verts;
  for (const auto& loop : localA)
    for (const vec2& v : loop) verts.push_back(v);
  const int sourceBoundary = static_cast<int>(verts.size());
  for (const auto& loop : localB)
    for (const vec2& v : loop) verts.push_back(v);
  if (sourceBoundary == 0 || sourceBoundary == (int)verts.size()) return s;

  const mb::VertexMerge merge = mb::MergeVerts(verts, s.eps);

  const double maxAcceptedDist2 =
      (kNearCoincidenceMaxRatio * s.eps) * (kNearCoincidenceMaxRatio * s.eps);
  for (int i = 0; i < sourceBoundary; ++i) {
    const int mi = merge.inputVert2Merged[i];
    if (mi < 0) continue;
    const vec2 pi = merge.verts[mi];
    for (int j = sourceBoundary; j < (int)verts.size(); ++j) {
      const int mj = merge.inputVert2Merged[j];
      if (mj < 0 || mj == mi) continue;
      const vec2 d = merge.verts[mj] - pi;
      const double d2 = d.x * d.x + d.y * d.y;
      if (d2 >= maxAcceptedDist2) continue;
      ++s.nearPairs;
      const double ratio = std::sqrt(d2) / s.eps;
      if (ratio > s.maxRatio) s.maxRatio = ratio;
    }
  }
  return s;
}

SimplePolygon StarPolygon(const std::vector<double>& radii) {
  SimplePolygon ring;
  const int n = static_cast<int>(radii.size());
  constexpr double kPi = 3.14159265358979323846;
  for (int i = 0; i < n; ++i) {
    const double r = 0.1 + std::fabs(radii[i]);
    const double th = 2.0 * kPi * i / n;
    ring.push_back({r * std::cos(th), r * std::sin(th)});
  }
  return ring;
}

std::vector<DriftSample> RunCorpus() {
  std::vector<DriftSample> out;

  // Parked seed: forced coincident vertex (B[4] = A[1]) drifted by FillByRule.
  {
    SimplePolygon a = StarPolygon(
        {1.0169016983060246, 1000.0, 1.0, 578.85382959129936, 0.0, 0.0});
    SimplePolygon b =
        StarPolygon({0.0, 1000.0, 999.83083173100238, 7.275875880519048,
                     726.89009231357352, 3.880747251022969});
    for (auto& v : b) v.x += 0.0030378301550779696;
    b[4] = a[1];
    out.push_back(Measure("DegenerateCoincidentVertexUnion seed",
                          CrossSection(a), CrossSection(b)));
  }

  // Disjoint squares.
  out.push_back(Measure("Disjoint squares (unit, far apart)",
                        CrossSection::Square({1, 1}, true),
                        CrossSection::Square({1, 1}, true).Translate({5, 0})));

  // Overlapping squares.
  out.push_back(
      Measure("Overlapping squares (unit, half-shift)",
              CrossSection::Square({1, 1}, true),
              CrossSection::Square({1, 1}, true).Translate({0.5, 0.3})));

  // Mixed scales.
  out.push_back(Measure(
      "Tiny vs large squares (1 vs 1e6)", CrossSection::Square({1, 1}, true),
      CrossSection::Square({1e6, 1e6}, true).Translate({100, 100})));

  // Concentric circles.
  out.push_back(Measure("Concentric circles (r=1 vs r=2)",
                        CrossSection::Circle(1.0, 64),
                        CrossSection::Circle(2.0, 64)));

  // Touching squares share a vertex by construction.
  {
    SimplePolygon a = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    SimplePolygon b = {{1, 1}, {2, 1}, {2, 2}, {1, 2}};
    out.push_back(
        Measure("Squares touching at (1,1)", CrossSection(a), CrossSection(b)));
  }

  // Random star pairs at varied scales/coincidence.
  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> rad(0.0, 1.0);
  std::uniform_real_distribution<double> scale(0.0, 3.0);
  std::uniform_real_distribution<double> shift(-0.01, 0.01);
  for (int i = 0; i < 50; ++i) {
    const int n = 4 + (i % 5);
    const double s = std::pow(10.0, scale(rng));
    std::vector<double> ra(n), rb(n);
    for (int k = 0; k < n; ++k) {
      ra[k] = s * rad(rng);
      rb[k] = s * rad(rng);
    }
    SimplePolygon a = StarPolygon(ra);
    SimplePolygon b = StarPolygon(rb);
    for (auto& v : b) v.x += s * shift(rng);
    if (i % 7 == 0) b[n / 2] = a[n / 3];
    char name[64];
    std::snprintf(name, sizeof(name), "Random star pair #%02d (n=%d, s=%.2g)",
                  i, n, s);
    out.push_back(Measure(name, CrossSection(a), CrossSection(b)));
  }

  // Chained ops: A.Offset.Boolean(B) exercises propagated tolerance + Round
  // arcTolerance inflation.
  {
    auto a = CrossSection::Square({10, 10}, true)
                 .Offset(0.5, CrossSection::JoinType::Round);
    auto b = CrossSection::Circle(0.01, 16).Translate({5, 5});
    out.push_back(Measure("Square.Offset(0.5,Round) vs tiny circle", a, b));
  }

  // Decompose+recompose stress: many components, each inherits parent tol.
  {
    std::vector<CrossSection> pieces;
    for (int i = 0; i < 16; ++i) {
      pieces.push_back(
          CrossSection::Square({0.5, 0.5}, true).Translate({1.0 * i, 0.0}));
    }
    auto big = CrossSection::BatchBoolean(pieces, manifold::OpType::Add);
    out.push_back(Measure("16-square BatchBoolean vs unit circle", big,
                          CrossSection::Circle(0.3).Translate({4, 0})));
  }

  return out;
}

}  // namespace

int main() {
  const auto samples = RunCorpus();
  std::vector<double> ratios;
  ratios.reserve(samples.size());
  double maxRatio = 0.0;
  std::string maxName = "(none)";
  for (const auto& s : samples) {
    if (s.eps <= 0.0 || s.nearPairs == 0) continue;
    ratios.push_back(s.maxRatio);
    if (s.maxRatio > maxRatio) {
      maxRatio = s.maxRatio;
      maxName = s.name;
    }
  }
  std::sort(ratios.begin(), ratios.end());

  std::printf("%-60s %12s %8s %12s\n", "case", "eps", "near#", "maxRatio");
  std::printf("%-60s %12s %8s %12s\n", "----", "---", "-----", "--------");
  for (const auto& s : samples) {
    std::printf("%-60s %12.3e %8d %12.3f\n", s.name.c_str(), s.eps, s.nearPairs,
                s.maxRatio);
  }
  std::printf("\n%zu samples with near-coincidence pairs (cutoff %.0f eps)\n",
              ratios.size(), kNearCoincidenceMaxRatio);
  std::printf("max ratio = %.3f (%s)\n", maxRatio, maxName.c_str());
  if (!ratios.empty()) {
    const auto pct = [&](double q) {
      const size_t i = std::min<size_t>(ratios.size() - 1,
                                        static_cast<size_t>(q * ratios.size()));
      return ratios[i];
    };
    std::printf("p50=%.3f  p95=%.3f  p99=%.3f\n", pct(0.50), pct(0.95),
                pct(0.99));
  }
  return 0;
}
