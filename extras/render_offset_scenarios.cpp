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
// Diagnostic only (not shipped). Emits per-scenario input+output polygon JSON
// for whichever CrossSection backend's libmanifold this links against. Build it
// in a clipper2 and a boolean2 configuration, then merge with
// compose_offset_renders.py; render_offset_c2_b2.sh drives both.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "manifold/cross_section.h"

using namespace manifold;
using JT = CrossSection::JoinType;

static const char* kBackend =
#ifdef MANIFOLD_CROSS_SECTION_BACKEND_BOOLEAN2
    "boolean2";
#else
    "clipper2";
#endif

static void emitRing(FILE* f, const SimplePolygon& r) {
  fputc('[', f);
  for (size_t i = 0; i < r.size(); ++i)
    fprintf(f, "%s[%.9g,%.9g]", i ? "," : "", r[i].x, r[i].y);
  fputc(']', f);
}
static void emitPolys(FILE* f, const Polygons& p) {
  fputc('[', f);
  for (size_t i = 0; i < p.size(); ++i) {
    if (i) fputc(',', f);
    emitRing(f, p[i]);
  }
  fputc(']', f);
}
static double areaOf(const Polygons& p) {
  double a = 0;
  for (const auto& r : p) {
    if (r.size() < 3) continue;
    double s = 0;
    const vec2 o = r[0];
    for (size_t i = 0; i < r.size(); ++i) {
      const vec2 u = r[i], v = r[(i + 1) % r.size()];
      s += (u.x - o.x) * (v.y - o.y) - (v.x - o.x) * (u.y - o.y);
    }
    a += 0.5 * s;
  }
  return a;
}
static size_t numVert(const Polygons& p) {
  size_t n = 0;
  for (const auto& r : p) n += r.size();
  return n;
}

static SimplePolygon Star(int pts, double rOuter, double rInner) {
  SimplePolygon s;
  const double kPi = 3.14159265358979323846;
  for (int i = 0; i < 2 * pts; ++i) {
    const double r = (i % 2) ? rInner : rOuter;
    const double th = kPi * i / pts;
    s.push_back({r * std::cos(th), r * std::sin(th)});
  }
  return s;
}

struct Scenario {
  std::string name;
  Polygons input;
  CrossSection result;
};

int main() {
  std::vector<Scenario> scen;

  // 1. Notched square ("U" slot) + Round: reflex-corner heal + arc
  // tessellation.
  {
    SimplePolygon p = {{0, 0},     {4, 0},     {4, 4},   {2.5, 4},
                       {2.5, 1.5}, {1.5, 1.5}, {1.5, 4}, {0, 4}};
    scen.push_back({"Notched square - Round +0.5",
                    {p},
                    CrossSection(p).Offset(0.5, JT::Round, 2.0, 32)});
  }
  // 2. 5-point star + Miter: sharp convex points (miter limit) + reflex heal.
  {
    auto s = Star(5, 2.0, 0.8);
    scen.push_back({"5-pt star - Miter +0.3 (limit 4)",
                    {s},
                    CrossSection(s).Offset(0.3, JT::Miter, 4.0)});
  }
  // 3. Same star + Bevel: flat-cut corners (contrast with miter).
  {
    auto s = Star(5, 2.0, 0.8);
    scen.push_back({"5-pt star - Bevel +0.3",
                    {s},
                    CrossSection(s).Offset(0.3, JT::Bevel)});
  }
  // 4. Square with a square hole + Round dilate: hole shrinks, outer gains
  // arcs.
  {
    Polygons ph = {{{0, 0}, {5, 0}, {5, 5}, {0, 5}},
                   {{1.5, 1.5}, {1.5, 3.5}, {3.5, 3.5}, {3.5, 1.5}}};
    scen.push_back({"Square+hole - Round +0.4", ph,
                    CrossSection(ph).Offset(0.4, JT::Round, 2.0, 32)});
  }
  // 5. Thin spike + Square: sharp base corners get square caps.
  {
    SimplePolygon p = {{0, 0}, {6, 0}, {3, 0.04}};
    scen.push_back({"Thin spike - Square +0.3",
                    {p},
                    CrossSection(p).Offset(0.3, JT::Square)});
  }
  // 6. Thin spike + Round: same shape, arc caps (contrast with square).
  {
    SimplePolygon p = {{0, 0}, {6, 0}, {3, 0.04}};
    scen.push_back({"Thin spike - Round +0.3",
                    {p},
                    CrossSection(p).Offset(0.3, JT::Round, 2.0, 32)});
  }
  // 7. Razor spike + Square: half-angle near 90 deg, where the square cap is
  // numerically fragile (cap-retention case).
  {
    SimplePolygon p = {{0, 0}, {6, 0}, {3, 0.0015}};
    scen.push_back({"Razor spike (h=1.5e-3) - Square +0.3",
                    {p},
                    CrossSection(p).Offset(0.3, JT::Square)});
  }
  // 8. Union of 3 circles, then Round offset: boolean junctions + arc joins.
  {
    auto c = CrossSection::Circle(1.0, 48);
    auto u = c + c.Translate({1.4, 0.0}) + c.Translate({0.7, 1.2});
    scen.push_back({"3-circle union - Round +0.3", u.ToPolygons(),
                    u.Offset(0.3, JT::Round, 2.0, 48)});
  }
  // 9. Two squares with a 0.4 gap + Round +0.3: offset closes the gap, merging
  // two rings into one (topology change).
  {
    Polygons two = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}},
                    {{1.4, 0}, {2.4, 0}, {2.4, 1}, {1.4, 1}}};
    scen.push_back({"Two squares (gap 0.4) - Round +0.3 merges", two,
                    CrossSection(two).Offset(0.3, JT::Round, 2.0, 32)});
  }
  // 10. 12-point star + Round: many convex+reflex joins at fine arc resolution.
  {
    auto s = Star(12, 2.0, 1.3);
    scen.push_back({"12-pt star - Round +0.25",
                    {s},
                    CrossSection(s).Offset(0.25, JT::Round, 2.0, 48)});
  }

  FILE* f = stdout;
  fprintf(f, "{\"backend\":\"%s\",\"scenarios\":[", kBackend);
  for (size_t i = 0; i < scen.size(); ++i) {
    const Polygons out = scen[i].result.ToPolygons();
    if (i) fputc(',', f);
    fprintf(f,
            "{\"name\":\"%s\",\"numVert\":%zu,\"area\":%.6g,\"rings\":%zu,"
            "\"input\":",
            scen[i].name.c_str(), numVert(out), areaOf(out), out.size());
    emitPolys(f, scen[i].input);
    fprintf(f, ",\"output\":");
    emitPolys(f, out);
    fputc('}', f);
  }
  fprintf(f, "]}\n");
  return 0;
}
