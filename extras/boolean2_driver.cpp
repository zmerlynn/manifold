// Copyright 2026 The Manifold Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Standalone boolean2 driver for 2D overlap-removal verification,
// diagnostics, corpus oracles, fuzzing, and benchmarking.
//
// Single-translation-unit driver following Julian Smith,
// "Towards robust inexact geometric computation", UCAM-CL-TR-766
// (2009), Chapter 7, with the BVH adaptation sketched by elalish in
// github.com/elalish/manifold/issues/289.
//
// References (cited inline as "Smith §X" / "Smith fig X" below):
//   [Smith09] Julian M. Smith, "Towards robust inexact geometric
//             computation", PhD dissertation, University of Cambridge,
//             Technical Report UCAM-CL-TR-766, 2009.
//             https://www.cl.cam.ac.uk/techreports/UCAM-CL-TR-766.pdf
//
//             - Chapter 6: figs 6.4 and 6.5 (degenerate-vertex test
//               cases used by `Smith fig 6.5(d)` in the test battery).
//             - Chapter 7: overall algorithm.
//             - §7.4 onwards: robust Bentley-Ottmann.
//             - §7.7 / fig 7.16: iterate-to-fixed-point bound (≤2
//               iterations under symbolic positions).
//             - Table 7.3: PolySet2 (canonical sub-edges with signed
//               multiplicity; see `Canonicalize` and `OutEdge` here).
//             - Chapter 8: floating-point error analysis;
//               u = 2⁻⁵³ for IEEE 754 doubles, per-intersection
//               bound α = √153·u·L (used by `EpsilonFromScale`).
//
//   [Issue289] https://github.com/elalish/manifold/issues/289
//             elalish's BVH-adapted sketch (no sweep line),
//             which this driver implements with Smith's symbolic
//             perturbation as the FP-robust core.
//
// Build (from the manifold repo root):
//   g++ -std=c++17 -O2 -I include -I src -DMANIFOLD_PAR=-1 \
//     -ffp-contract=off -fexcess-precision=standard \
//     extras/boolean2_driver.cpp -o boolean2_driver
// Run:
//   ./boolean2_driver                  # full test battery
//   ./boolean2_driver diagnose 0       # diagnostic dump for one case
//   ./boolean2_driver deepfuzz 100     # broader randomized verification
//   ./boolean2_driver offsetfuzz 50    # offset-path adversarial coverage
//   ./boolean2_driver time --repeat 5 clipper2corpus
//                                     # stable aggregate corpus timing
//   ./boolean2_driver vsclipper2 [all|clipper2|offsets|jts|mfogel]
//                                     # head-to-head bench (needs
//                                     # -DBOOLEAN2_WITH_CLIPPER2 +
//                                     # scripts/fetch_clipper2_test_data.sh)
//
// Algorithm shape (matches manifold internals where possible):
//   - Spatial queries via boolean2's 2D BVH (Morton-sorted leaves).
//   - Edge-edge intersection via a trim-and-`Interpolate` symbolic kernel
//     atop the `Shadows` orientation primitive from src/shared.h.
//   - Vertex equality via `manifold::DisjointSets` lock-free union-find.
//   - Face traversal via DCEL (the same structure manifold's 3D
//     `Manifold::Impl::halfedge_` uses), winding-rayed per face from a
//     point on the LEFT side of any boundary half-edge.
//   - Iterate to fixed point per Smith §7.7 (default maxIter=2, his bound).
//
// Build can be single-threaded with MANIFOLD_PAR=-1 or use manifold's TBB
// parallel helpers with MANIFOLD_PAR=1. Optional comparison and corpus modes
// link Clipper2 and/or libmanifold.
//
// Deferred for graduation to manifold's mainline build:
//   - Mechanical std::vector → manifold::Vec rename (Vec is the project's
//     CPU/GPU-portable vector). Drop-in API-compatible for our usage; the
//     conversion is ~135 declarations and is best done in the same patch
//     that wires the driver into the build system, so the type churn
//     and the build-graph churn land together.
//   - ZoneScoped Tracy markers at phase boundaries.
//   - ExecutionContext::Impl* ctx threading for parallelism dispatch.
//   - Internal namespace (boolean2::detail) hiding everything except
//     the public Simplify / Boolean2D entry points.

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <climits>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include "../src/collider.h"
#include "../src/disjoint_sets.h"
#include "../src/iters.h"
#include "../src/parallel.h"
#include "../src/shared.h"
#include "../src/utils.h"
#include "manifold/common.h"
#ifdef BOOLEAN2_WITH_MANIFOLD
#include "manifold/manifold.h"
#endif

#include "cross_section/boolean2/predicates.h"
#include "cross_section/boolean2/bvh.h"
#include "cross_section/boolean2/vertex_merge.h"
#include "cross_section/boolean2/edge_vert_lists.h"
#include "cross_section/boolean2/intersections.h"
#include "cross_section/boolean2/canonicalize.h"
#include "cross_section/boolean2/winding_filter.h"
#include "cross_section/boolean2/driver.h"
#include "cross_section/boolean2/boolean2.h"
#include "cross_section/boolean2/iterate.h"
#include "cross_section/boolean2/offset.h"

namespace manifold {
namespace boolean2 {

using manifold::Box;
using manifold::CCW;
using manifold::Collider;
using manifold::MakeSimpleRecorder;
using manifold::OpType;
using manifold::vec2;
using manifold::vec3;
using manifold::VecView;
using manifold::la::dot;
using manifold::la::length;

// =============================================================================
// Property checks.
// =============================================================================

// Compute per-vertex (out_mult - in_mult) balance for an edge set.
template <typename Edge>
std::map<int, int> ComputeBalance(const std::vector<Edge>& edges) {
  std::map<int, int> b;
  for (const auto& e : edges) {
    b[e.v0] += e.mult;
    b[e.v1] -= e.mult;
  }
  return b;
}

// Topological invariant: for vertices that exist in BOTH input and output
// (i.e. original input verts surviving the merge), the output balance must
// equal the input balance. New intersection-created verts must have balance
// zero. This generalizes "sum in == sum out" to handle open inputs.
//
// Note: empty output trivially passes; closed input edges contribute
// equal in/out at every vertex, expected balance is 0 everywhere, and an
// empty edge list also has 0 balance everywhere. Callers that care about
// "did we drop everything?" must additionally check `r.edges` non-empty.
//
// numMergedVerts = number of verts emitted by vertex merging (the first N
// in r.verts are either original or merged-from-original; intersection
// verts come after).
bool CheckTopologicalValidity(const OverlapResult& r,
                              const std::vector<EdgeM>& inputEdges,
                              const std::vector<int>& inputRemap,
                              int numMergedVerts) {
  // Expected balance: remap input edges through `inputRemap`, compute
  // per-output-vertex balance contribution from the input.
  std::vector<EdgeM> remapped;
  for (const auto& e : inputEdges) {
    int a = inputRemap[e.v0];
    int b = inputRemap[e.v1];
    if (a != b) remapped.push_back({a, b, e.mult});
  }
  auto expected = ComputeBalance(remapped);
  auto actual = ComputeBalance(r.edges);
  bool ok = true;
  for (int v = 0; v < static_cast<int>(r.verts.size()); ++v) {
    int exp = expected.count(v) ? expected[v] : 0;
    int act = actual.count(v) ? actual[v] : 0;
    // Original-input verts (idx < numMergedVerts): output balance must
    // match input. Intersection verts (idx >= numMergedVerts): must be 0.
    int target = (v < numMergedVerts) ? exp : 0;
    if (act != target) {
      std::cerr << "[FAIL] vertex " << v << " (orig=" << (v < numMergedVerts)
                << ") expected balance " << target << ", got " << act
                << std::endl;
      ok = false;
    }
  }
  return ok;
}

// =============================================================================
// Test patterns from Smith §7.7.
// =============================================================================

// Random topological polygon: N points on unit circle in random angular
// order. Produces a self-intersecting closed curve.
std::pair<std::vector<vec2>, std::vector<EdgeM>> RandomTopologicalPolygon(
    int n, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u01(0, 1);
  std::vector<vec2> verts(n);
  for (int i = 0; i < n; ++i) {
    double theta = 2 * M_PI * u01(rng);
    verts[i] = {std::cos(theta), std::sin(theta)};
  }
  // Random permutation for connectivity.
  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::shuffle(order.begin(), order.end(), rng);
  std::vector<EdgeM> edges;
  edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    edges.push_back({order[i], order[(i + 1) % n], 1});
  }
  return {std::move(verts), std::move(edges)};
}

// Polygonal star: n equispaced points on unit circle, connected by skipping
// k vertices. Highly self-intersecting and symmetric.
std::pair<std::vector<vec2>, std::vector<EdgeM>> PolygonalStar(int n,
                                                               int skip) {
  std::vector<vec2> verts(n);
  for (int i = 0; i < n; ++i) {
    double theta = 2 * M_PI * i / n;
    verts[i] = {std::cos(theta), std::sin(theta)};
  }
  std::vector<EdgeM> edges;
  edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    edges.push_back({i, (i + skip) % n, 1});
  }
  return {std::move(verts), std::move(edges)};
}

// Apply a translation that brings the bounding box near a power of 2.
// Smith §7.7's "displacement attack" makes representable grid spacing
// comparable to feature sizes.
void Displace(std::vector<vec2>* verts, double offset) {
  for (auto& v : *verts) {
    v.x += offset;
    v.y += offset;
  }
}

// L-shaped concave polygon (CCW). 6 verts, no self-intersection.
std::pair<std::vector<vec2>, std::vector<EdgeM>> LShape() {
  std::vector<vec2> v = {
      {0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2},
  };
  std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1},
                          {3, 4, 1}, {4, 5, 1}, {5, 0, 1}};
  return {std::move(v), std::move(e)};
}

// Square (CCW) with a square hole (CW). Two disjoint loops.
std::pair<std::vector<vec2>, std::vector<EdgeM>> SquareWithHole() {
  std::vector<vec2> v = {
      // Outer CCW square
      {0, 0},
      {4, 0},
      {4, 4},
      {0, 4},
      // Inner CW square (hole)
      {1, 1},
      {1, 3},
      {3, 3},
      {3, 1},
  };
  std::vector<EdgeM> e = {
      {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // outer CCW
      {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},  // inner CW
  };
  return {std::move(v), std::move(e)};
}

// Two disjoint CCW squares.
std::pair<std::vector<vec2>, std::vector<EdgeM>> TwoSquares() {
  std::vector<vec2> v = {
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {3, 0}, {3, 1}, {2, 1},
  };
  std::vector<EdgeM> e = {
      {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},
      {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},
  };
  return {std::move(v), std::move(e)};
}

// =============================================================================
// Main: run a battery of tests and report.
// =============================================================================
struct TestCase {
  std::string name;
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  double eps;
};

// Run RemoveOverlaps2D and verify topological invariant. Returns the result
// for further property checks.
OverlapResult RunCase(const TestCase& tc, bool* allPass) {
  std::cout << "=== " << tc.name << " ===" << std::endl;
  std::cout << "  Input: " << tc.verts.size() << " verts, " << tc.edges.size()
            << " edges, eps=" << tc.eps << std::endl;
  auto r = RemoveOverlaps2D(tc.verts, tc.edges, tc.eps);
  std::cout << "  Output: " << r.verts.size() << " verts, " << r.edges.size()
            << " edges" << std::endl;
  bool topo =
      CheckTopologicalValidity(r, tc.edges, r.inputRemap, r.numMergedVerts);
  std::cout << "  Topological validity: " << (topo ? "PASS" : "FAIL")
            << std::endl;
  if (!topo) *allPass = false;
  return r;
}

// Idempotence: feeding the output back through RemoveOverlaps2D should
// produce the same result. Uses Fingerprint (position-based, eps/100
// quantized) for cross-iteration vert-renumbering tolerance, matching
// the displacement fuzz's check.
bool CheckIdempotence(const OverlapResult& first, double eps) {
  std::vector<EdgeM> asEdgeM;
  for (const auto& e : first.edges) {
    asEdgeM.push_back({e.v0, e.v1, e.mult});
  }
  auto second = RemoveOverlaps2D(first.verts, asEdgeM, eps);
  bool ok = (Fingerprint(first, eps) == Fingerprint(second, eps));
  std::cout << "  Idempotence: " << (ok ? "PASS" : "FAIL") << std::endl;
  return ok;
}

}  // namespace boolean2
}  // namespace manifold

// Diagnostic: run ONE failing case and dump intermediate state.
// Invoked via `./boolean2_driver diagnose <seed> [kPow] [n]`. Default
// kPow=30, n=50, but any DeepFuzz parameter combo can be targeted.
// IMPORTANT: the seed→RNG mapping must match DeepFuzz exactly.
namespace manifold {
namespace boolean2 {
void Diagnose(uint64_t seed, int kPow = 30, int n = 50) {
  const double offset = std::ldexp(1.5, kPow);
  const double eps = EpsilonFromScale(offset);
  std::cerr << "=== Diagnose: kPow=" << kPow << " n=" << n << " seed=" << seed
            << " offset=" << offset << " eps=" << eps << " ===\n";

  // Match the seed mapping used by both main()'s displacement fuzz and
  // DeepFuzz: `seed + 1000 * kPow`.
  auto [vIn, eIn] = RandomTopologicalPolygon(n, seed + 1000ull * kPow);
  Displace(&vIn, offset);
  std::cerr << "Input: " << vIn.size() << " verts, " << eIn.size()
            << " edges\n";

  // Run pipeline manually to inspect intermediate state.
  auto merge = MergeVerts(vIn, eps);
  std::cerr << "After vertex merge: " << merge.verts.size() << " verts (was "
            << vIn.size() << ")\n";

  auto edges = RemapAndCollapse(eIn, merge.remap);
  std::cerr << "After edge collapse: " << edges.size() << " edges (was "
            << eIn.size() << ")\n";

  std::vector<Box2> diagEdgeBoxes(edges.size());
  for (size_t i = 0; i < edges.size(); ++i)
    diagEdgeBoxes[i] = BoxOf2DEdge(merge.verts[edges[i].v0],
                                   merge.verts[edges[i].v1], eps);
  BVH diagBvh = BVHBuildFromBoxes(diagEdgeBoxes);
  auto lists =
      BuildEdgeVertLists(edges, merge.verts, eps, diagEdgeBoxes, diagBvh);
  size_t totalListSize = 0;
  for (const auto& l : lists) totalListSize += l.size();
  std::cerr << "After near-vertex indexing: " << totalListSize
            << " total list entries across " << edges.size() << " edges\n";

  const int beforeIntersections = static_cast<int>(merge.verts.size());
  std::vector<std::vector<int>> vertEdges;
  auto diagPairs = CollectIntersectionPairs(edges, diagEdgeBoxes, diagBvh);
  FindAndInsertIntersections(edges, &merge.verts, &lists, &vertEdges, eps,
                             diagEdgeBoxes, diagBvh, diagPairs);
  const int afterIntersections = static_cast<int>(merge.verts.size());
  std::cerr << "After intersections: " << merge.verts.size()
            << " verts (added " << (afterIntersections - beforeIntersections)
            << ")\n";

  // Structural merge is omitted here: production
  // RemoveOverlaps2D uses union-find over verts that share a parent edge
  // and fall within 10*eps. The diagnostic shows the post-intersection
  // (pre-structural-merge) state so any residual near-coincident
  // intersection clusters surface in the dump below.
  std::cerr << "After structural re-merge: skipped (production uses "
               "structural merge; see RemoveOverlaps2D)\n";

  auto canon = Canonicalize(edges, lists);
  std::cerr << "After canonicalization: " << canon.edges.size()
            << " sub-edges (after multiplicity collapse)\n";

  // Vertex balance from canonicalized sub-edges.
  std::map<int, int> canonBalance;
  for (const auto& edge : canon.edges) {
    canonBalance[edge.vMin] += edge.mult;
    canonBalance[edge.vMax] -= edge.mult;
  }

  // Expected balance from remapped input.
  std::vector<EdgeM> remappedInput;
  for (const auto& e : eIn) {
    int a = merge.remap[e.v0];
    int b = merge.remap[e.v1];
    if (a != b) remappedInput.push_back({a, b, e.mult});
  }
  auto expBalance = ComputeBalance(remappedInput);

  // Find imbalanced vertices.
  const int numMergedVerts = beforeIntersections;
  std::cerr << "\n=== Per-vertex balance check ===\n";
  int badCount = 0;
  for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
    int target =
        (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
    int actual = canonBalance.count(v) ? canonBalance[v] : 0;
    if (actual != target) {
      ++badCount;
      if (badCount <= 8) {
        std::cerr << "  v" << v << " (orig=" << (v < numMergedVerts)
                  << ") pos=(" << merge.verts[v].x << "," << merge.verts[v].y
                  << ") expected=" << target << " actual=" << actual << "\n";
        // Dump sub-edges touching this vertex.
        std::cerr << "    sub-edges touching v" << v << ":\n";
        int touched = 0;
        for (const auto& edge : canon.edges) {
          if (edge.vMin == v || edge.vMax == v) {
            std::cerr << "      (" << edge.vMin << "," << edge.vMax
                      << ") mult=" << edge.mult << "\n";
            if (++touched >= 12) {
              std::cerr << "      ... (more truncated)\n";
              break;
            }
          }
        }
      }
    }
  }
  std::cerr << "Total imbalanced vertices after canonicalization: "
            << badCount << "\n";

  // Now also run the winding filter and check balance again.
  auto out = FilterByWindingDCEL(canon, merge.verts);
  std::cerr << "\nAfter winding filter: " << out.size()
            << " output edges (from " << canon.edges.size()
            << " canonical sub-edges)\n";
  std::map<int, int> outBalance;
  for (const auto& e : out) {
    outBalance[e.v0] += e.mult;
    outBalance[e.v1] -= e.mult;
  }
  int outBadCount = 0;
  for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
    int target =
        (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
    int actual = outBalance.count(v) ? outBalance[v] : 0;
    if (actual != target) {
      ++outBadCount;
      if (outBadCount <= 8) {
        std::cerr << "  v" << v << " (orig=" << (v < numMergedVerts)
                  << ") pos=(" << merge.verts[v].x << "," << merge.verts[v].y
                  << ") expected=" << target << " actual=" << actual << "\n";
        std::cerr << "    output edges touching v" << v << ":\n";
        int touched = 0;
        for (const auto& e : out) {
          if (e.v0 == v || e.v1 == v) {
            std::cerr << "      (" << e.v0 << "->" << e.v1
                      << ") mult=" << e.mult << "\n";
            if (++touched >= 12) {
              std::cerr << "      ... (more truncated)\n";
              break;
            }
          }
        }
      }
    }
  }
  std::cerr << "Total imbalanced vertices after winding filter: "
            << outBadCount << "\n";

  // Deep dive on first imbalanced vertex: dump winding decisions for every
  // canonical sub-edge touching it, including the wL/wR values.
  if (outBadCount > 0) {
    int targetV = -1;
    for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
      int target =
          (v < numMergedVerts) ? (expBalance.count(v) ? expBalance[v] : 0) : 0;
      int actual = outBalance.count(v) ? outBalance[v] : 0;
      if (actual != target) {
        targetV = v;
        break;
      }
    }
    std::cerr << "\n=== Deep dive on v" << targetV << " ===\n";
    // Also dump nearby vertices and check distances for §7.6-style grouping.
    std::cerr << "  v" << targetV << " neighbors (within 10*eps):\n";
    const double thresh2 = (10 * eps) * (10 * eps);
    for (int v = 0; v < static_cast<int>(merge.verts.size()); ++v) {
      if (v == targetV) continue;
      vec2 d = merge.verts[v] - merge.verts[targetV];
      double dist2 = dot(d, d);
      if (dist2 <= thresh2) {
        std::cerr << "    v" << v << " dist=" << std::sqrt(dist2) << " ("
                  << (std::sqrt(dist2) / eps) << "*eps)\n";
      }
    }
    // Per-edge perpendicular-offset ray-cast for diagnostic visualization
    // only (the actual winding filter uses DCEL face traversal).
    const double ofs = eps * 1e-1;
    for (const auto& edge : canon.edges) {
      if (edge.vMin != targetV && edge.vMax != targetV) continue;
      vec2 p0 = merge.verts[edge.vMin];
      vec2 p1 = merge.verts[edge.vMax];
      vec2 mid = (p0 + p1) * 0.5;
      vec2 d = p1 - p0;
      double len = length(d);
      vec2 perp = {-d.y / len, d.x / len};
      vec2 leftPt = mid + perp * ofs;
      vec2 rightPt = mid + perp * -ofs;
      int wL = CastWindingRay(leftPt, canon, merge.verts);
      int wR = CastWindingRay(rightPt, canon, merge.verts);
      std::cerr << "  sub-edge (" << edge.vMin << "," << edge.vMax
                << ") mult=" << edge.mult << " | mid=(" << mid.x - offset
                << "," << mid.y - offset << ") perp=(" << perp.x << ","
                << perp.y << ") | wL=" << wL << " wR=" << wR
                << (((wL > 0) != (wR > 0)) ? " => KEPT" : " => DROPPED")
                << "\n";
    }
  }

  // Cluster proximity: how many distinct verts are within eps of each other?
  int closeCount = 0;
  const double eps2 = eps * eps;
  for (int i = 0; i < static_cast<int>(merge.verts.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(merge.verts.size()); ++j) {
      vec2 d = merge.verts[i] - merge.verts[j];
      if (dot(d, d) <= eps2) ++closeCount;
    }
  }
  std::cerr << "After structural re-merge: " << closeCount
            << " pairs of distinct verts within eps of each other\n";

  // Now run iterate-to-fixed-point and report on the iterated result.
  std::cerr << "\n=== IterateToFixedPoint ===\n";
  int iters = 0;
  IterStatus iterStatus = IterStatus::Converged;
  auto rIter = IterateToFixedPoint(vIn, eIn, eps, /*maxIter=*/8, &iters,
                                   &iterStatus);
  const char* statusStr =
      iterStatus == IterStatus::Converged
          ? "Converged"
          : (iterStatus == IterStatus::Cycled ? "Cycled" : "MaxedOut");
  std::cerr << "Iters: " << iters << " status: " << statusStr << "\n";
  std::cerr << "Iterated result: " << rIter.verts.size() << " verts, "
            << rIter.edges.size() << " edges, numMergedVerts="
            << rIter.numMergedVerts << "\n";

  // Check topology validity on the iterated result with composed remap.
  bool iterValid = CheckTopologicalValidity(rIter, eIn, rIter.inputRemap,
                                            rIter.numMergedVerts);
  std::cerr << "Iterated topology valid: " << (iterValid ? "YES" : "NO")
            << "\n";

  // Sanity-check pass 1 the way DeepFuzz does it (the production
  // RemoveOverlaps2D pipeline, NOT the legacy diagnostic above).
  std::cerr << "\n=== Pass 1 (RemoveOverlaps2D, the production path) ===\n";
  auto pass1 = RemoveOverlaps2D(vIn, eIn, eps);
  bool pass1Valid = CheckTopologicalValidity(pass1, eIn, pass1.inputRemap,
                                             pass1.numMergedVerts);
  std::cerr << "Pass 1 verts=" << pass1.verts.size()
            << " edges=" << pass1.edges.size()
            << " numMergedVerts=" << pass1.numMergedVerts << " valid="
            << (pass1Valid ? "YES" : "NO") << "\n";

  // Find imbalanced verts and dump their canonical sub-edge environment.
  if (!pass1Valid) {
    std::cerr << "\n=== Imbalanced vert analysis ===\n";
    // Compute expected balance from input edges remapped through
    // pass1.inputRemap (same as CheckTopologicalValidity).
    std::vector<EdgeM> remapped;
    for (const auto& ie : eIn) {
      int a = pass1.inputRemap[ie.v0], b = pass1.inputRemap[ie.v1];
      if (a != b) remapped.push_back({a, b, ie.mult});
    }
    auto exp = ComputeBalance(remapped);
    std::map<int, int> act;
    for (const auto& oe : pass1.edges) {
      act[oe.v0] += oe.mult;
      act[oe.v1] -= oe.mult;
    }
    for (int v = 0; v < (int)pass1.verts.size(); ++v) {
      int target = (v < pass1.numMergedVerts)
                       ? (exp.count(v) ? exp[v] : 0)
                       : 0;
      int actual = act.count(v) ? act[v] : 0;
      if (actual != target) {
        const vec2 p = pass1.verts[v];
        std::cerr << "  v" << v << " (orig=" << (v < pass1.numMergedVerts)
                  << ") pos=(" << (p.x - offset) << "," << (p.y - offset)
                  << ") expected=" << target << " got=" << actual << "\n";
        // Edges in pass1.edges touching this vert.
        std::cerr << "    output edges touching:\n";
        int tch = 0;
        for (const auto& oe : pass1.edges) {
          if (oe.v0 == v || oe.v1 == v) {
            std::cerr << "      " << oe.v0 << "→" << oe.v1
                      << " mult=" << oe.mult << "\n";
            if (++tch >= 6) break;
          }
        }
      }
    }
    // For each imbalanced vert, find canonical sub-edges incident to it,
    // and check how the winding filter voted on each.
    std::cerr << "\n=== Winding-filter verdicts for imbalanced verts ===\n";
    // Re-run the production pipeline to canonical sub-edges (structural
    // structural re-merge included; same as RemoveOverlaps2D through
    // canonicalization.
    // We want canon as the production pipeline produces it.
    // Note: pass1.edges is after winding filtering; we need the
    // canonicalized edges. So re-run.
    auto m3 = MergeVerts(vIn, eps);
    auto e3 = RemapAndCollapse(eIn, m3.remap);
    std::vector<Box2> e3Boxes(e3.size());
    for (size_t i = 0; i < e3.size(); ++i)
      e3Boxes[i] =
          BoxOf2DEdge(m3.verts[e3[i].v0], m3.verts[e3[i].v1], eps);
    BVH e3Bvh = BVHBuildFromBoxes(e3Boxes);
    auto l3 = BuildEdgeVertLists(e3, m3.verts, eps, e3Boxes, e3Bvh);
    std::vector<std::vector<int>> ve3;
    auto e3Pairs = CollectIntersectionPairs(e3, e3Boxes, e3Bvh);
    FindAndInsertIntersections(e3, &m3.verts, &l3, &ve3, eps, e3Boxes, e3Bvh,
                               e3Pairs);
    // Structural re-merge (copy of the production code).
    {
      DisjointSets uf(static_cast<int>(m3.verts.size()));
      const double mergeThresh2 = (10.0 * eps) * (10.0 * eps);
      for (size_t a = 0; a < m3.verts.size(); ++a) {
        if (a >= ve3.size() || ve3[a].empty()) continue;
        for (size_t b = a + 1; b < m3.verts.size(); ++b) {
          if (b >= ve3.size() || ve3[b].empty()) continue;
          bool shared = false;
          for (int e : ve3[a]) if (VESetContains(ve3[b], e)) { shared = true; break; }
          if (!shared) continue;
          vec2 d = m3.verts[b] - m3.verts[a];
          if (dot(d, d) > mergeThresh2) continue;
          uf.unite((int)a, (int)b);
        }
      }
      std::map<int, std::pair<vec2, int>> sums;
      for (size_t i = 0; i < m3.verts.size(); ++i) {
        int r = uf.find((int)i);
        auto& s = sums[r];
        s.first = s.first + m3.verts[i];
        s.second += 1;
      }
      if (sums.size() < m3.verts.size()) {
        std::map<int, int> rootToNew;
        std::vector<vec2> nv;
        for (const auto& [root, sc] : sums) {
          rootToNew[root] = (int)nv.size();
          nv.push_back(sc.first * (1.0 / sc.second));
        }
        std::vector<int> rmp(m3.verts.size());
        for (size_t i = 0; i < m3.verts.size(); ++i)
          rmp[i] = rootToNew[uf.find((int)i)];
        for (auto& ed : e3) { ed.v0 = rmp[ed.v0]; ed.v1 = rmp[ed.v1]; }
        for (auto& list : l3) {
          for (auto& v : list) v = rmp[v];
          list.erase(std::unique(list.begin(), list.end()), list.end());
        }
        m3.verts = std::move(nv);
      }
    }
    auto canon = Canonicalize(e3, l3);
    std::cerr << "  canonical sub-edges: "
              << canon.edges.size() << "\n";

    // For each imbalanced vert, dump every canonical edge touching it and
    // the winding filter's verdict.
    std::set<int> imbalanced{101, 110, 112, 195};
    const double ofs = eps * 1e-1;
    for (int target : imbalanced) {
      if (target >= (int)m3.verts.size()) continue;
      vec2 vp = m3.verts[target];
      std::cerr << "  v" << target << " pos=(" << (vp.x - offset) << ","
                << (vp.y - offset) << ")\n";
      int sub_count = 0, kept = 0, dropped = 0;
      for (const auto& cedge : canon.edges) {
        if (cedge.vMin != target && cedge.vMax != target) continue;
        ++sub_count;
        const int other = (cedge.vMin == target) ? cedge.vMax : cedge.vMin;
        const vec2 op = m3.verts[other];
        vec2 p0 = vp, p1 = op;
        if (cedge.vMin != target) std::swap(p0, p1);
        vec2 mid = (p0 + p1) * 0.5;
        vec2 d = p1 - p0;
        const double len = length(d);
        if (len == 0) continue;
        vec2 perp = {-d.y / len, d.x / len};
        vec2 lpt = mid + perp * ofs, rpt = mid - perp * ofs;
        int wL = CastWindingRay(lpt, canon, m3.verts);
        int wR = CastWindingRay(rpt, canon, m3.verts);
        bool keep = (wL > 0) != (wR > 0);
        if (keep) ++kept; else ++dropped;
        std::cerr << "    (" << cedge.vMin << "↔" << cedge.vMax
                  << ") mult=" << cedge.mult << " wL=" << wL << " wR=" << wR
                  << (keep ? " KEEP" : " DROP") << " other.pos=("
                  << (op.x - offset) << "," << (op.y - offset) << ")\n";
        if (sub_count >= 8) {
          std::cerr << "    ...\n";
          break;
        }
      }
      std::cerr << "    " << kept << " kept, " << dropped
                << " dropped at this vert\n";
    }
  }
  // Also report what the test is expecting per-vertex.
  std::cerr << "Composed remap (original input vert → final vert idx):\n";
  for (size_t i = 0; i < rIter.inputRemap.size(); ++i) {
    if (i < 10 || i == rIter.inputRemap.size() - 1)
      std::cerr << "  v" << i << " → " << rIter.inputRemap[i] << "\n";
  }
}
}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// Polygon corpus runner.
//
// Reads the manifold project's existing curated test corpus at
// test/polygons/polygon_corpus.txt - 100 named adversarial polygons
// (Eberly, Sliver, Comb, KissingZigzag, BarelyValid, Tricky, Sponge4a,
// CondensedMatter64, etc.) used in the project's triangulation tests.
// Each entry has form (matches test/polygon_test.cpp's
// RegisterPolygonTestsFile):
//
//   <name> <expected_tri_count> <eps> <num_loops>
//   <loop1_vert_count>
//   <x> <y>     (loop1_vert_count times)
//   <loop2_vert_count>   (only if multi-loop)
//   <x> <y>     (loop2_vert_count times)
//   ...
//
// eps = -1.0 means "use a sensible default for this scale". For the
// overlap-removal tests we infer one from the bounding box.
//
// Runs each polygon through Simplify, checks topology
// validity, and verifies signed area is preserved within a tolerance
// (since these are mostly already-valid polygons - any large area drift
// indicates the algorithm modified something it shouldn't have, or
// failed to clean up a self-intersection it should have). The
// expected_tri_count is informational for now (it's a triangulation
// metric, not directly an overlap-removal one).
// =============================================================================
namespace manifold {
namespace boolean2 {

struct CorpusEntry {
  std::string name;
  manifold::Polygons polys;
  double eps;          // -1 if not specified in file (caller infers)
  int expected_tris;   // triangulation triangle count (informational)
};

inline std::vector<CorpusEntry> LoadCorpus(const std::string& path) {
  std::vector<CorpusEntry> entries;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return entries;
  }
  while (in.good()) {
    CorpusEntry entry;
    int num_loops = 0;
    if (!(in >> entry.name >> entry.expected_tris >> entry.eps >> num_loops))
      break;
    for (int i = 0; i < num_loops; ++i) {
      int loop_size = 0;
      if (!(in >> loop_size)) break;
      manifold::SimplePolygon loop;
      loop.reserve(loop_size);
      for (int j = 0; j < loop_size; ++j) {
        double x, y;
        if (!(in >> x >> y)) break;
        loop.push_back({x, y});
      }
      entry.polys.push_back(std::move(loop));
    }
    entries.push_back(std::move(entry));
  }
  return entries;
}

inline void RunCorpus(const std::string& path) {
  auto entries = LoadCorpus(path);
  std::cout << "=== Polygon corpus: " << entries.size() << " entries from "
            << path << " ===\n";
  int total = 0;
  int topo_ok = 0;
  int topo_fail = 0;
  int area_preserved = 0;     // |drift| < 1%
  int area_drifted = 0;       // 1% <= |drift| < 10%
  int area_collapsed = 0;     // |drift| >= 10%
  int area_input_zero = 0;    // input had zero net area
  std::vector<std::tuple<std::string, double>> driftList;
  for (const auto& entry : entries) {
    ++total;
    double eps = entry.eps > 0 ? entry.eps : InferEps(entry.polys, {});
    if (eps <= 0) eps = EpsilonFromScale(1.0);
    auto out = Simplify(entry.polys, eps);
    // Topology: convert input polys to (verts, edges) form, run pipeline,
    // run the project's CheckTopologicalValidity.
    auto [v, e] = PolygonsToInput(entry.polys);
    auto r = RemoveOverlaps2D(v, e, eps);
    bool topoOk = CheckTopologicalValidity(r, e, r.inputRemap,
                                           r.numMergedVerts);
    if (topoOk) ++topo_ok;
    else ++topo_fail;
    const double inArea = std::fabs(TotalSignedArea(entry.polys));
    const double outArea = std::fabs(TotalSignedArea(out));
    if (inArea < 1e-12) {
      ++area_input_zero;
    } else {
      const double drift = std::fabs(outArea - inArea) / inArea;
      if (drift >= 0.10) {
        ++area_collapsed;
        driftList.emplace_back(entry.name, drift);
      } else if (drift >= 0.01) {
        ++area_drifted;
        driftList.emplace_back(entry.name, drift);
      } else {
        ++area_preserved;
      }
    }
  }
  std::cout << "  Total: " << total << "\n";
  std::cout << "  Topology valid:   " << topo_ok << " / " << total << "\n";
  std::cout << "  Topology INVALID: " << topo_fail << " / " << total << "\n";
  std::cout << "  Area preserved (drift < 1%):       " << area_preserved
            << "\n";
  std::cout << "  Area drifted   (1% <= drift < 10%): " << area_drifted
            << "\n";
  std::cout << "  Area collapsed (drift >= 10%):       " << area_collapsed
            << "\n";
  std::cout << "  Input had zero net area: " << area_input_zero
            << " (degenerate or pure-CW polygons; not an oracle case)\n";
  if (!driftList.empty()) {
    std::sort(driftList.begin(), driftList.end(),
              [](const auto& a, const auto& b) {
                return std::get<1>(a) > std::get<1>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 20):\n";
    for (size_t i = 0; i < driftList.size() && i < 20; ++i) {
      auto& [name, dr] = driftList[i];
      std::cout << "    " << name << "  drift=" << (dr * 100.0) << "%\n";
    }
  }
}

}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// Clipper2 Polygons.txt corpus runner.
//
// Independent third-party oracle from Clipper2's test suite (195 numbered
// cases at build/_deps/clipper2-src/Tests/Polygons.txt). Each case has:
//
//   CAPTION: <n>. [<note>]
//   CLIPTYPE: UNION | INTERSECTION | DIFFERENCE
//   FILLRULE: NONZERO | EVENODD | NEGATIVE | POSITIVE
//   SOL_AREA: <expected absolute area, or -1 if not validated>
//   SOL_COUNT: <expected loop count, or -1>
//   SUBJECTS
//   x,y, x,y, ...        (one polygon per line; integer coords)
//   ...
//   [CLIPS
//   x,y, x,y, ...]
//   <blank line>
//
// Dispatch model: Clipper2's CLIPTYPE × FILLRULE semantics is
//
//   result = boolean_op(fill_rule(subjects), fill_rule(clips))
//
// i.e. each input is independently canonicalized to a CCW arrangement
// under the declared fill rule, then the boolean op combines them. The
// boolean2's public Simplify/Boolean2D/Xor APIs assume CCW-canonical
// input (predicates use Smith's w>0 union convention), so they don't
// match Clipper2's semantics directly when inputs are self-intersecting
// or in screen-coord (y-down → CW in y-up math) orientation.
//
// To fairly compare, we run each case through a pre-fill stage:
//
//   - FILLRULE: NONZERO → RemoveOverlaps2D with predicate w != 0
//   - FILLRULE: EVENODD → RemoveOverlaps2D with predicate w & 1
//   - FILLRULE: NEGATIVE/POSITIVE → unsupported (6 cases skipped).
//
// After pre-fill, both sides are CCW-canonical, and we combine via:
//
//   - UNION        → Boolean2D::Add        (w>0 on combined CCW)
//   - INTERSECTION → Boolean2D::Intersect  (w>1 on combined CCW)
//   - DIFFERENCE   → Boolean2D::Subtract   (B mult flipped, then w>0)
//
// The eps inferred from the bounding box (Smith α-budget at L = max
// |coord|) lands at ~1e-12 for these int-coordinate inputs - well below
// the unit grid the corpus uses, so iteration shouldn't be triggered
// by snap-induced changes.
// =============================================================================
namespace manifold {
namespace boolean2 {

struct Clipper2Case {
  int n = 0;
  std::string caption;
  std::string cliptype;
  std::string fillrule;
  double solArea = -1.0;
  int solCount = -1;
  manifold::Polygons subjects;
  manifold::Polygons clips;
};

// Parse one polygon line "x,y, x,y, x,y" into a SimplePolygon. Coordinates
// are comma-separated; pairs are flat in the list. Tolerates trailing
// whitespace/commas.
inline manifold::SimplePolygon ParseClipper2PolyLine(const std::string& line) {
  manifold::SimplePolygon poly;
  std::vector<double> nums;
  std::string buf;
  for (char c : line) {
    if (c == ',' || std::isspace(static_cast<unsigned char>(c))) {
      if (!buf.empty()) {
        try { nums.push_back(std::stod(buf)); } catch (...) {}
        buf.clear();
      }
    } else {
      buf.push_back(c);
    }
  }
  if (!buf.empty()) {
    try { nums.push_back(std::stod(buf)); } catch (...) {}
  }
  for (size_t i = 0; i + 1 < nums.size(); i += 2) {
    poly.push_back({nums[i], nums[i + 1]});
  }
  return poly;
}

inline std::vector<Clipper2Case> LoadClipper2Cases(const std::string& path) {
  std::vector<Clipper2Case> cases;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    if (path.find("clipper2_") != std::string::npos) {
      std::cerr << "  hint: scripts/fetch_clipper2_test_data.sh\n";
    }
    return cases;
  }
  Clipper2Case cur;
  bool inCase = false;
  enum class Section { None, Subjects, Clips };
  Section sec = Section::None;
  std::string line;
  auto flush = [&]() {
    if (inCase && !cur.caption.empty()) cases.push_back(std::move(cur));
    cur = Clipper2Case{};
    inCase = false;
    sec = Section::None;
  };
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
      line.pop_back();
    if (line.empty()) {
      flush();
      continue;
    }
    auto starts = [&](const char* k) {
      size_t n = std::strlen(k);
      return line.size() >= n && line.compare(0, n, k) == 0;
    };
    if (starts("CAPTION:")) {
      flush();
      inCase = true;
      cur.caption = line.substr(8);
      // Strip leading whitespace and trailing dot.
      size_t p = 0;
      while (p < cur.caption.size() && std::isspace(static_cast<unsigned char>(cur.caption[p]))) ++p;
      cur.caption = cur.caption.substr(p);
      cur.n = std::atoi(cur.caption.c_str());
      continue;
    }
    if (!inCase) continue;
    if (starts("CLIPTYPE:")) { cur.cliptype = line.substr(9); auto p = cur.cliptype.find_first_not_of(" \t"); cur.cliptype = (p == std::string::npos) ? "" : cur.cliptype.substr(p); continue; }
    if (starts("FILLRULE:")) { cur.fillrule = line.substr(9); auto p = cur.fillrule.find_first_not_of(" \t"); cur.fillrule = (p == std::string::npos) ? "" : cur.fillrule.substr(p); continue; }
    if (starts("SOL_AREA:")) { cur.solArea = std::stod(line.substr(9)); continue; }
    if (starts("SOL_COUNT:")) { cur.solCount = std::atoi(line.c_str() + 10); continue; }
    if (line == "SUBJECTS" || line == "SUBJECTS_OPEN") { sec = Section::Subjects; continue; }
    if (line == "CLIPS" || line == "CLIPS_OPEN") { sec = Section::Clips; continue; }
    auto poly = ParseClipper2PolyLine(line);
    if (poly.size() >= 3) {
      if (sec == Section::Clips) cur.clips.push_back(std::move(poly));
      else cur.subjects.push_back(std::move(poly));
    }
  }
  flush();
  return cases;
}

enum class Clipper2Reach { Supported, Unsupported };

inline Clipper2Reach ClassifyCase(const Clipper2Case& c) {
  if (c.fillrule == "NEGATIVE" || c.fillrule == "POSITIVE")
    return Clipper2Reach::Unsupported;
  if (c.fillrule != "NONZERO" && c.fillrule != "EVENODD")
    return Clipper2Reach::Unsupported;
  if (c.cliptype != "UNION" && c.cliptype != "INTERSECTION" &&
      c.cliptype != "DIFFERENCE")
    return Clipper2Reach::Unsupported;
  return Clipper2Reach::Supported;
}

// Pre-fill polys under the declared fill rule. Output is a CCW-canonical
// arrangement (no self-overlap). Reuses RemoveOverlaps2D with a custom
// predicate so we don't perturb the public API surface.
inline manifold::Polygons FillUnderRule(const manifold::Polygons& polys,
                                        const std::string& rule, double eps) {
  if (polys.empty()) return {};
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  for (const auto& loop : polys) {
    if (loop.size() < 3) continue;
    const int base = static_cast<int>(verts.size());
    const int n = static_cast<int>(loop.size());
    for (const auto& v : loop) verts.push_back(v);
    for (int i = 0; i < n; ++i)
      edges.push_back({base + i, base + ((i + 1) % n), 1});
  }
  if (verts.empty()) return {};
  const WindRule windRule =
      (rule == "EVENODD") ? WindRule::EvenOdd : WindRule::NonZero;
  auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, windRule);
  return OutEdgesToPolygons(r.verts, r.edges);
}

// Run one case under the universal pre-fill model. Returns the output
// polygons; sets `reach` to indicate whether the case was dispatched.
inline manifold::Polygons RunClipper2Case(const Clipper2Case& c,
                                          Clipper2Reach* reach,
                                          double* epsOut = nullptr) {
  *reach = ClassifyCase(c);
  if (*reach == Clipper2Reach::Unsupported) return {};
  const double eps = InferEps(c.subjects, c.clips);
  if (epsOut) *epsOut = eps;
  auto a = FillUnderRule(c.subjects, c.fillrule, eps);
  auto b = FillUnderRule(c.clips, c.fillrule, eps);
  if (c.cliptype == "UNION") {
    if (b.empty()) return a;
    if (a.empty()) return b;
    return Boolean2D(a, b, OpType::Add, eps);
  }
  if (c.cliptype == "INTERSECTION") {
    if (a.empty() || b.empty()) return {};
    return Boolean2D(a, b, OpType::Intersect, eps);
  }
  if (c.cliptype == "DIFFERENCE") {
    if (a.empty()) return {};
    if (b.empty()) return a;
    return Boolean2D(a, b, OpType::Subtract, eps);
  }
  return {};
}

inline void RunClipper2Corpus(const std::string& path) {
  auto cases = LoadClipper2Cases(path);
  std::cout << "=== Clipper2 corpus: " << cases.size() << " cases from "
            << path << " ===\n\n";
  // Per-case TSV (case_n, op, fillrule, structure, sol_area, output_area)
  // is written to /tmp/clipper2_proto_areas.tsv so the per-case area can
  // be diffed against Clipper2's own current output. SOL_AREA in the
  // corpus is older Clipper2 stored output and may differ from what
  // current Clipper2 produces for the same input.
  std::ofstream tsv("/tmp/clipper2_proto_areas.tsv");
  tsv << "n\top\tfillrule\tstructure\tsol_area\tproto_area\n";
  int total = 0;
  int run = 0, unsupported = 0;
  int topoOk = 0, topoFail = 0;
  int areaPreserved = 0;     // < 1% of SOL_AREA
  int areaDrifted = 0;       // 1% .. 10%
  int areaCollapsed = 0;     // > 10% (or zero output where SOL_AREA > 0)
  int solSkipped = 0;        // SOL_AREA <= 0 (no oracle)
  int outZeroOk = 0;         // SOL_AREA == 0 and we returned ~0 area
  std::vector<std::tuple<int, std::string, double, double, double>> drifters;
  // (case_n, descriptor, drift_pct, outArea, solArea)
  for (const auto& c : cases) {
    ++total;
    Clipper2Reach reach;
    double eps = 0.0;
    manifold::Polygons out;
    try {
      out = RunClipper2Case(c, &reach, &eps);
    } catch (const std::exception& ex) {
      std::cerr << "  case " << c.n << " (" << c.cliptype << "+" << c.fillrule
                << "): EXCEPTION " << ex.what() << "\n";
      ++topoFail;
      continue;
    }
    std::string desc = c.cliptype + "+" + c.fillrule +
                       (c.clips.empty() ? "/sub-only" : "/sub+clip");
    if (reach == Clipper2Reach::Unsupported) {
      ++unsupported;
      continue;
    }
    ++run;
    {
      const double outArea = std::fabs(TotalSignedArea(out));
      tsv << c.n << '\t' << c.cliptype << '\t' << c.fillrule << '\t'
          << (c.clips.empty() ? "subonly" : "subclip") << '\t'
          << c.solArea << '\t' << outArea << '\n';
    }

    // Topology check on the output polygons (treat as a fresh input).
    auto [v, e] = PolygonsToInput(out);
    if (!v.empty()) {
      auto r2 = RemoveOverlaps2D(v, e, eps);
      bool topo = CheckTopologicalValidity(r2, e, r2.inputRemap,
                                           r2.numMergedVerts);
      // For a valid output, re-running RemoveOverlaps2D should be a no-op
      // (idempotent under sub-eps drift). Topology failure here flags
      // either an invalid output or a corner case the pipeline mishandles
      // when fed back its own result.
      if (topo) ++topoOk;
      else ++topoFail;
    } else {
      // Empty output: only "valid" when SOL_AREA == 0 or unsupported.
      ++topoOk;
    }

    const double outArea = std::fabs(TotalSignedArea(out));
    if (c.solArea < 0) {
      ++solSkipped;
    } else if (c.solArea == 0.0) {
      if (outArea < 1.0) ++outZeroOk;
      else { ++areaCollapsed; drifters.emplace_back(c.n, desc, outArea, outArea, c.solArea); }
    } else {
      const double drift = std::fabs(outArea - c.solArea) / c.solArea;
      if (drift < 0.01) ++areaPreserved;
      else if (drift < 0.10) {
        ++areaDrifted;
        drifters.emplace_back(c.n, desc, drift, outArea, c.solArea);
      } else {
        ++areaCollapsed;
        drifters.emplace_back(c.n, desc, drift, outArea, c.solArea);
      }
    }
  }
  std::cout << "  Total cases:     " << total << "\n";
  std::cout << "  Run:             " << run << "\n";
  std::cout << "  Unsupported:     " << unsupported
            << " (NEGATIVE/POSITIVE fill rule)\n";
  std::cout << "\n  Of " << run << " cases run:\n";
  std::cout << "    Topology valid:   " << topoOk << " / " << run << "\n";
  std::cout << "    Topology INVALID: " << topoFail << " / " << run << "\n";
  std::cout << "    Area within  1% of SOL_AREA: " << areaPreserved << "\n";
  std::cout << "    Area within  1-10% of SOL_AREA: " << areaDrifted << "\n";
  std::cout << "    Area > 10% from SOL_AREA:     " << areaCollapsed << "\n";
  std::cout << "    SOL_AREA == 0, returned ~0:    " << outZeroOk << "\n";
  std::cout << "    SOL_AREA < 0 (no oracle):      " << solSkipped << "\n";
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 30):\n";
    for (size_t i = 0; i < drifters.size() && i < 30; ++i) {
      auto& [n, desc, drift, outA, solA] = drifters[i];
      std::cout << "    #" << n << "  " << desc
                << "  drift=" << (drift * 100.0) << "%"
                << "  out=" << outA << "  sol=" << solA << "\n";
    }
  }
}

}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// mfogel/polygon-clipping end-to-end fixtures.
//
// Independent third-party corpus harvested from real bug reports against
// the JS polygon-clipping library. Each fixture lives in its own dir
// under build/_deps/mfogel-polygon-clipping-src/test/end-to-end/<name>/:
//
//   args.geojson   - GeoJSON FeatureCollection of N input operands.
//                    Each Feature is a Polygon or MultiPolygon.
//   union.geojson, intersection.geojson, difference.geojson,
//   xor.geojson    - expected outputs for each op (at most one per file).
//   all.geojson    - expected output that's the same for every op
//                    (typical for single-input self-cleanup cases).
//   broken-issue-* - args only; no expected output. Skip these.
//
// GeoJSON convention (RFC 7946): exterior rings CCW, interior (hole)
// rings CW. NONZERO fill rule under that orientation isolates outer−hole.
// The runner pre-fills each feature under NONZERO, then chains the
// boolean op left-to-right across features.
//
// The corpus is small (~120 dirs) but its naming is itself a degeneracy
// taxonomy: triple-coincident-segments, almost-colinear-segments,
// infinitely-thin-polygon, self-intersects-but-doesnt-cross, ...
// This is exactly the layer Clipper2's int-coord corpus doesn't cover.
// =============================================================================
namespace manifold {
namespace boolean2 {

// ---------------------------------------------------------------------------
// Minimal JSON value + parser. Hand-rolled to avoid a third-party
// dependency for this single-translation-unit diagnostic driver.
// Supports only what GeoJSON needs: object, array, string, number, bool,
// null. No escape parsing in strings beyond \" and \\ - none of the
// fixture files use richer escapes.
// ---------------------------------------------------------------------------
struct JsonValue {
  enum Kind { Null, Bool, Number, String, Array, Object };
  Kind kind = Null;
  bool b = false;
  double num = 0.0;
  std::string str;
  std::vector<JsonValue> arr;
  std::map<std::string, JsonValue> obj;
};

struct JsonParser {
  const std::string& s;
  size_t i = 0;
  explicit JsonParser(const std::string& src) : s(src) {}
  void skipWs() {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  }
  bool match(char c) {
    skipWs();
    if (i < s.size() && s[i] == c) { ++i; return true; }
    return false;
  }
  void expect(char c) {
    if (!match(c)) throw std::runtime_error(
        std::string("expected '") + c + "' at " + std::to_string(i));
  }
  JsonValue parseValue() {
    skipWs();
    if (i >= s.size()) throw std::runtime_error("unexpected eof");
    char c = s[i];
    if (c == '{') return parseObject();
    if (c == '[') return parseArray();
    if (c == '"') return parseString();
    if (c == 't' || c == 'f') return parseBool();
    if (c == 'n') return parseNull();
    return parseNumber();
  }
  JsonValue parseObject() {
    JsonValue v; v.kind = JsonValue::Object;
    expect('{');
    skipWs();
    if (match('}')) return v;
    for (;;) {
      JsonValue key = parseString();
      expect(':');
      v.obj[key.str] = parseValue();
      skipWs();
      if (match(',')) continue;
      expect('}');
      break;
    }
    return v;
  }
  JsonValue parseArray() {
    JsonValue v; v.kind = JsonValue::Array;
    expect('[');
    skipWs();
    if (match(']')) return v;
    for (;;) {
      v.arr.push_back(parseValue());
      skipWs();
      if (match(',')) continue;
      expect(']');
      break;
    }
    return v;
  }
  JsonValue parseString() {
    JsonValue v; v.kind = JsonValue::String;
    expect('"');
    while (i < s.size() && s[i] != '"') {
      if (s[i] == '\\' && i + 1 < s.size()) {
        char n = s[i + 1];
        if (n == '"' || n == '\\') { v.str.push_back(n); i += 2; continue; }
        v.str.push_back(s[i]); ++i;
      } else {
        v.str.push_back(s[i]); ++i;
      }
    }
    expect('"');
    return v;
  }
  JsonValue parseBool() {
    JsonValue v; v.kind = JsonValue::Bool;
    if (s.compare(i, 4, "true") == 0) { v.b = true; i += 4; }
    else if (s.compare(i, 5, "false") == 0) { v.b = false; i += 5; }
    else throw std::runtime_error("bad bool");
    return v;
  }
  JsonValue parseNull() {
    JsonValue v; v.kind = JsonValue::Null;
    if (s.compare(i, 4, "null") == 0) { i += 4; return v; }
    throw std::runtime_error("bad null");
  }
  JsonValue parseNumber() {
    JsonValue v; v.kind = JsonValue::Number;
    size_t start = i;
    if (s[i] == '-' || s[i] == '+') ++i;
    while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) ||
                            s[i] == '.' || s[i] == 'e' || s[i] == 'E' ||
                            s[i] == '-' || s[i] == '+'))
      ++i;
    v.num = std::stod(s.substr(start, i - start));
    return v;
  }
};

inline JsonValue ParseJsonFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("open: " + path);
  std::ostringstream oss; oss << in.rdbuf();
  std::string text = oss.str();
  JsonParser p(text);
  return p.parseValue();
}

// mfogel's polygon-clipping library normalizes ring winding by position:
// the first ring of each polygon is outer (CCW), subsequent rings are
// holes (CW). It does this regardless of the input's actual winding -
// the corpus has explicit "backward-ring-winding-order" test cases that
// rely on this auto-normalization. Strict GeoJSON conformance (RFC 7946)
// agrees with this convention but doesn't mandate enforcement; mfogel
// enforces. We do the same here so NONZERO pre-fill on the combined
// rings reproduces outer−holes correctly.
inline void NormalizeMfogelWinding(manifold::Polygons& rings) {
  for (size_t i = 0; i < rings.size(); ++i) {
    const double a = SignedArea(rings[i]);
    const bool wantCCW = (i == 0);
    if ((wantCCW && a < 0) || (!wantCCW && a > 0))
      std::reverse(rings[i].begin(), rings[i].end());
  }
}

// Extract a list of polygons from one Feature's geometry.
// Polygon: coordinates = [ring1, ring2, ...] → one polygon with multiple
//          rings (exterior + holes).
// MultiPolygon: coordinates = [[ring1, ring2], [ring1], ...] → one
//          polygon per outer entry.
// Returns vector<Polygons> (each entry = one polygon, possibly with holes).
inline std::vector<manifold::Polygons> ExtractFeatureGeometry(
    const JsonValue& feature) {
  std::vector<manifold::Polygons> polys;
  auto gIt = feature.obj.find("geometry");
  if (gIt == feature.obj.end() || gIt->second.kind != JsonValue::Object)
    return polys;
  const auto& geom = gIt->second;
  auto tIt = geom.obj.find("type");
  auto cIt = geom.obj.find("coordinates");
  if (tIt == geom.obj.end() || cIt == geom.obj.end()) return polys;
  const std::string& gtype = tIt->second.str;
  auto parseRing = [](const JsonValue& ring) {
    manifold::SimplePolygon out;
    for (const auto& pt : ring.arr) {
      if (pt.arr.size() >= 2)
        out.push_back({pt.arr[0].num, pt.arr[1].num});
    }
    // Drop the closing duplicate vertex if present (RFC 7946 requires it).
    if (out.size() >= 2 && out.front().x == out.back().x &&
        out.front().y == out.back().y)
      out.pop_back();
    return out;
  };
  // mfogel drops the entire polygon if its outer (first) ring is
  // degenerate (< 3 verts after closure-removal). The corpus has the
  // explicit `rings-with-no-area` case for this. Subsequent (hole)
  // rings can be silently dropped if degenerate; the polygon survives.
  auto collect = [&](const JsonValue& ringList,
                     manifold::Polygons& dest) -> bool {
    bool first = true;
    for (const auto& ring : ringList.arr) {
      auto r = parseRing(ring);
      if (first) {
        first = false;
        if (r.size() < 3) return false;
        dest.push_back(std::move(r));
      } else {
        if (r.size() >= 3) dest.push_back(std::move(r));
      }
    }
    return !dest.empty();
  };
  if (gtype == "Polygon") {
    manifold::Polygons one;
    if (collect(cIt->second, one)) {
      NormalizeMfogelWinding(one);
      polys.push_back(std::move(one));
    }
  } else if (gtype == "MultiPolygon") {
    for (const auto& poly : cIt->second.arr) {
      manifold::Polygons one;
      if (collect(poly, one)) {
        NormalizeMfogelWinding(one);
        polys.push_back(std::move(one));
      }
    }
  }
  return polys;
}

// Flatten a list-of-Polygons (one per outer polygon in a MultiPolygon)
// into a single Polygons (concat all rings). NONZERO pre-fill on the
// concatenation reproduces the multipolygon's covered region under
// GeoJSON's CCW-outer / CW-hole convention.
inline manifold::Polygons FlattenMulti(
    const std::vector<manifold::Polygons>& multi) {
  manifold::Polygons flat;
  for (const auto& p : multi)
    for (const auto& r : p) flat.push_back(r);
  return flat;
}

// Top-level extractor: returns one polygon-group per Feature in the file.
// Handles either a single Feature or a FeatureCollection at the root.
inline std::vector<std::vector<manifold::Polygons>> ExtractFeatures(
    const std::string& path) {
  JsonValue root = ParseJsonFile(path);
  std::vector<std::vector<manifold::Polygons>> out;
  if (root.kind != JsonValue::Object) return out;
  auto tIt = root.obj.find("type");
  if (tIt == root.obj.end()) return out;
  if (tIt->second.str == "FeatureCollection") {
    auto fIt = root.obj.find("features");
    if (fIt == root.obj.end()) return out;
    for (const auto& f : fIt->second.arr)
      out.push_back(ExtractFeatureGeometry(f));
  } else if (tIt->second.str == "Feature") {
    out.push_back(ExtractFeatureGeometry(root));
  }
  return out;
}

struct MfogelCase {
  std::string name;
  std::vector<manifold::Polygons> features;  // one entry per arg feature
  std::map<std::string, manifold::Polygons> expected;  // op -> expected (flat)
};

inline std::vector<MfogelCase> LoadMfogelCorpus(const std::string& dir) {
  std::vector<MfogelCase> cases;
  // Walk subdirs without <filesystem> dependency: ls via popen.
  std::string cmd = "ls -1 " + dir;
  FILE* fp = popen(cmd.c_str(), "r");
  if (!fp) return cases;
  char buf[1024];
  std::vector<std::string> names;
  while (std::fgets(buf, sizeof(buf), fp)) {
    std::string n = buf;
    if (!n.empty() && n.back() == '\n') n.pop_back();
    if (!n.empty()) names.push_back(n);
  }
  pclose(fp);
  for (const auto& name : names) {
    const std::string caseDir = dir + "/" + name;
    const std::string args = caseDir + "/args.geojson";
    std::ifstream chk(args);
    if (!chk) continue;
    MfogelCase c;
    c.name = name;
    try {
      auto feats = ExtractFeatures(args);
      for (auto& f : feats) {
        // Collapse each feature's polygon list into one rings-list. We
        // pre-fill under NONZERO so the orientation-encoded holes work.
        c.features.push_back(FlattenMulti(f));
      }
    } catch (const std::exception& e) {
      std::cerr << "  parse failed for " << name << ": " << e.what() << "\n";
      continue;
    }
    for (const char* op :
         {"union", "intersection", "difference", "xor", "all"}) {
      const std::string p = caseDir + "/" + op + ".geojson";
      std::ifstream chk2(p);
      if (!chk2) continue;
      try {
        auto feats = ExtractFeatures(p);
        manifold::Polygons exp;
        for (auto& f : feats) {
          auto flat = FlattenMulti(f);
          for (auto& r : flat) exp.push_back(std::move(r));
        }
        c.expected[op] = std::move(exp);
      } catch (...) {}
    }
    cases.push_back(std::move(c));
  }
  return cases;
}

// Run one (case, op) pair through boolean2 as a single-pass
// combined-input pipeline.
//
// Earlier this chained `Boolean2D(acc, filled[i], op)` left-to-right.
// The chain breaks on N>=3-feature cases like `issue-44` because the
// intermediate `acc` accumulates T-junction vertices (extra colinear
// verts on shared edges) that the next feature doesn't have, and the
// pipeline mishandles those collinear runs at the next merge.
//
// Single-pass dispatch sidesteps the chain by feeding all features into
// one RemoveOverlaps2D call with appropriate multiplicities + predicate:
//
//   union          → mult=+1 each, predicate w != 0
//   xor            → mult=+1 each, predicate (w & 1)
//   intersection   → mult=+1 each, predicate w >= N (covered by every input)
//   difference     → A mult=+1, others mult=-1, predicate w > 0
//
// Each feature is pre-filled under NONZERO first so its rings are
// CCW-canonical before being merged into the combined input.
inline manifold::Polygons RunMfogelOp(const MfogelCase& c,
                                      const std::string& op, double eps) {
  if (c.features.empty()) return {};
  std::vector<manifold::Polygons> filled;
  for (const auto& f : c.features) {
    if (f.empty()) { filled.push_back({}); continue; }
    filled.push_back(FillUnderRule(f, "NONZERO", eps));
  }
  if (op == "all") return filled[0];

  // Combine all features' edges into one input list with appropriate
  // per-feature multiplicities.
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  size_t nNonEmpty = 0;
  auto append = [&](const manifold::Polygons& polys, int mult) {
    if (polys.empty()) return;
    ++nNonEmpty;
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i)
        edges.push_back({base + i, base + ((i + 1) % n), mult});
    }
  };
  if (op == "difference") {
    append(filled[0], +1);
    for (size_t i = 1; i < filled.size(); ++i) append(filled[i], -1);
  } else {
    for (auto& f : filled) append(f, +1);
  }
  if (verts.empty()) return {};
  OverlapResult r;
  if (op == "union") {
    r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, WindRule::NonZero);
  } else if (op == "xor") {
    r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, WindRule::EvenOdd);
  } else if (op == "intersection") {
    // N-ary intersection: keep faces with w >= N. Parameterized
    // predicate, so we go through RemoveOverlaps2D's templated
    // overload with a lambda (only the diagnostic driver pays the
    // extra instantiation).
    const int N = static_cast<int>(filled.size());
    r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                         [N](int w) { return w >= N; });
  } else if (op == "difference") {
    r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false, WindRule::Add);
  } else {
    return {};
  }
  return OutEdgesToPolygons(r.verts, r.edges);
}

inline double InferEpsForFeatures(
    const std::vector<manifold::Polygons>& features) {
  manifold::Polygons all;
  for (const auto& f : features)
    for (const auto& r : f) all.push_back(r);
  return InferEps(all, {});
}

inline void RunMfogelCorpus(const std::string& dir) {
  auto cases = LoadMfogelCorpus(dir);
  std::cout << "=== mfogel/polygon-clipping corpus: " << cases.size()
            << " dirs from " << dir << " ===\n\n";
  int casesWithExpected = 0;
  int totalOps = 0;
  int topoOk = 0, topoFail = 0;
  int areaPreserved = 0, areaDrifted = 0, areaCollapsed = 0;
  int oracleEmptyOk = 0;          // expected = 0 area, we returned ~0
  int oracleEmptyFail = 0;        // expected = 0, we returned non-zero
  std::vector<std::tuple<std::string, std::string, double, double, double>>
      drifters;
  for (const auto& c : cases) {
    if (c.expected.empty()) continue;
    ++casesWithExpected;
    const double eps = InferEpsForFeatures(c.features);
    for (const auto& [op, expected] : c.expected) {
      ++totalOps;
      manifold::Polygons out;
      try {
        out = RunMfogelOp(c, op, eps);
      } catch (const std::exception& e) {
        std::cerr << "  EXC " << c.name << "/" << op << ": " << e.what()
                  << "\n";
        ++topoFail;
        continue;
      }
      // Topology check by re-feeding output into the pipeline.
      auto [v, e] = PolygonsToInput(out);
      if (!v.empty()) {
        auto r2 = RemoveOverlaps2D(v, e, eps);
        bool topo = CheckTopologicalValidity(r2, e, r2.inputRemap,
                                             r2.numMergedVerts);
        if (topo) ++topoOk;
        else ++topoFail;
      } else {
        ++topoOk;
      }
      const double outA = std::fabs(TotalSignedArea(out));
      const double expA = std::fabs(TotalSignedArea(expected));
      if (expA < 1e-9) {
        if (outA < 1e-9) ++oracleEmptyOk;
        else { ++oracleEmptyFail; drifters.emplace_back(c.name, op, outA, outA, 0.0); }
      } else {
        const double drift = std::fabs(outA - expA) / expA;
        if (drift < 0.01) ++areaPreserved;
        else if (drift < 0.10) {
          ++areaDrifted;
          drifters.emplace_back(c.name, op, drift, outA, expA);
        } else {
          ++areaCollapsed;
          drifters.emplace_back(c.name, op, drift, outA, expA);
        }
      }
    }
  }
  std::cout << "  Cases with expected-output files: " << casesWithExpected
            << " / " << cases.size() << "\n";
  std::cout << "  (Skipped: cases with only args.geojson - known broken in"
               " mfogel itself, e.g. broken-issue-*.)\n";
  std::cout << "  Total (case, op) pairs run:       " << totalOps << "\n\n";
  std::cout << "  Topology valid:   " << topoOk << " / " << totalOps << "\n";
  std::cout << "  Topology INVALID: " << topoFail << " / " << totalOps << "\n";
  std::cout << "  Area within  1%:  " << areaPreserved << "\n";
  std::cout << "  Area 1-10%:        " << areaDrifted << "\n";
  std::cout << "  Area >= 10%:       " << areaCollapsed << "\n";
  std::cout << "  Expected empty, returned ~0:  " << oracleEmptyOk << "\n";
  std::cout << "  Expected empty, returned >0: " << oracleEmptyFail << "\n";
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<2>(a) > std::get<2>(b);
              });
    std::cout << "\n  Drift cases (worst first, up to 30):\n";
    for (size_t i = 0; i < drifters.size() && i < 30; ++i) {
      auto& [name, op, drift, outA, expA] = drifters[i];
      std::cout << "    " << name << " / " << op
                << "  drift=" << (drift * 100.0) << "%"
                << "  out=" << outA << "  exp=" << expA << "\n";
    }
  }
}

}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// JTS robust/overlay corpus runner.
//
// JTS Topology Suite (LocationTech) ships an XML+WKB test format with
// cases harvested from real GIS bug reports (GEOS, PostGIS, QGIS, JTS
// itself, Shapely, OSGeo). The dominant op is `overlayAreaTest`, which
// checks the invariant
//
//     |A∪B| + |A∩B| = |A| + |B|
//
// This is a self-checking oracle - no expected-output WKB needed -
// and it's the strongest single cross-check available: a bug in either
// union or intersection that creates or loses area shows up immediately.
// Real-bug-derived inputs from production GIS systems mean the invariant
// is the right shape of test for them.
//
// XML+WKB parsing lives outside this driver in test/polygons/
// jts_to_corpus.py, which produces test/polygons/jts_overlay_corpus.txt.
// That keeps the WKB/XML/EWKB/WKT plumbing out of the algorithm code
// and makes the corpus inspectable as plain text. The same converter
// handles other JTS-format-compatible corpora (GEOS testxml, etc.) if
// added later.
//
// Format consumed here:
//
//   CASE <n> <op> <source>[#<idx>] [<expected_value>]
//   A <num_rings>
//   <num_pts> x y x y ...    (one ring per line; coords space-separated)
//   B <num_rings>            (omit if no B)
//   <num_pts> ...
//   <blank line>
//
// Supported ops: `overlayAreaTest` (invariant) and `unionArea` (oracle).
// =============================================================================
namespace manifold {
namespace boolean2 {

struct JtsCase {
  int n = 0;
  std::string op;
  std::string source;
  double expected = 0.0;
  bool hasExpected = false;
  manifold::Polygons a;  // flat ring list (orientation preserved)
  manifold::Polygons b;
};

inline std::vector<JtsCase> LoadJtsCorpus(const std::string& path) {
  std::vector<JtsCase> cases;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Failed to open " << path << "\n";
    return cases;
  }
  std::string line;
  JtsCase cur;
  bool inCase = false;
  // When we see "A <n>" or "B <n>" we then expect <n> ring lines.
  int ringsRemaining = 0;
  bool readingA = false;
  auto flush = [&]() {
    if (inCase && cur.n > 0) cases.push_back(std::move(cur));
    cur = JtsCase{};
    inCase = false;
    ringsRemaining = 0;
    readingA = false;
  };
  while (std::getline(in, line)) {
    while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
      line.pop_back();
    if (line.empty() || line[0] == '#') {
      if (line.empty()) flush();
      continue;
    }
    std::istringstream ss(line);
    if (ringsRemaining > 0) {
      // Ring line: <num_pts> x y x y ...
      int npts;
      ss >> npts;
      manifold::SimplePolygon ring;
      ring.reserve(npts);
      for (int i = 0; i < npts; ++i) {
        double x, y;
        if (!(ss >> x >> y)) break;
        ring.push_back({x, y});
      }
      if (ring.size() >= 3) {
        if (readingA) cur.a.push_back(std::move(ring));
        else cur.b.push_back(std::move(ring));
      }
      --ringsRemaining;
      continue;
    }
    std::string token;
    ss >> token;
    if (token == "CASE") {
      flush();
      inCase = true;
      ss >> cur.n >> cur.op >> cur.source;
      double v;
      if (ss >> v) { cur.expected = v; cur.hasExpected = true; }
    } else if (token == "A") {
      ss >> ringsRemaining;
      readingA = true;
    } else if (token == "B") {
      ss >> ringsRemaining;
      readingA = false;
    }
    // Other lines silently ignored (forward-compat with format extensions).
  }
  flush();
  return cases;
}

inline void RunJtsCorpus(const std::string& path) {
  auto cases = LoadJtsCorpus(path);
  std::cout << "=== JTS overlay corpus: " << cases.size()
            << " cases from " << path << " ===\n\n";
  size_t areaTests = 0, unionAreaTests = 0, skipped = 0;
  size_t topoOk = 0, topoFail = 0;
  size_t invariantHeld = 0;
  size_t invariantDriftSmall = 0;  // < 1%
  size_t invariantDriftLarge = 0;  // >= 1%
  size_t unionAreaOk = 0, unionAreaDrift = 0;
  std::vector<std::tuple<std::string, double, double, double, double>>
      drifters;  // (source, drift_pct, A∪B, A∩B, |A|+|B|)
  for (const auto& c : cases) {
    if (c.op != "overlayAreaTest" && c.op != "overlayareatest" &&
        c.op != "unionArea") {
      ++skipped;
      continue;
    }
    try {
      const double eps = InferEps(c.a, c.b);
      auto fa = FillUnderRule(c.a, "NONZERO", eps);
      if (c.op == "unionArea") {
        ++unionAreaTests;
        auto [v, e] = PolygonsToInput(fa);
        if (!v.empty()) {
          auto r2 = RemoveOverlaps2D(v, e, eps);
          if (CheckTopologicalValidity(r2, e, r2.inputRemap,
                                       r2.numMergedVerts)) ++topoOk;
          else ++topoFail;
        } else ++topoOk;
        const double a = std::fabs(TotalSignedArea(fa));
        if (c.hasExpected && c.expected > 0) {
          double drift = std::fabs(a - c.expected) / c.expected;
          if (drift < 0.01) ++unionAreaOk; else ++unionAreaDrift;
        } else {
          ++unionAreaOk;
        }
        continue;
      }
      ++areaTests;
      auto fb = FillUnderRule(c.b, "NONZERO", eps);
      const double Aa = std::fabs(TotalSignedArea(fa));
      const double Ab = std::fabs(TotalSignedArea(fb));
      manifold::Polygons unionAB, interAB;
      if (fa.empty()) unionAB = fb;
      else if (fb.empty()) unionAB = fa;
      else unionAB = Boolean2D(fa, fb, OpType::Add, eps);
      if (!fa.empty() && !fb.empty())
        interAB = Boolean2D(fa, fb, OpType::Intersect, eps);
      const double Au = std::fabs(TotalSignedArea(unionAB));
      const double Ai = std::fabs(TotalSignedArea(interAB));
      auto [v, e] = PolygonsToInput(unionAB);
      if (!v.empty()) {
        auto r2 = RemoveOverlaps2D(v, e, eps);
        if (CheckTopologicalValidity(r2, e, r2.inputRemap,
                                     r2.numMergedVerts)) ++topoOk;
        else ++topoFail;
      } else ++topoOk;
      const double lhs = Au + Ai;
      const double rhs = Aa + Ab;
      const double scale = std::max(rhs, 1e-12);
      const double drift = std::fabs(lhs - rhs) / scale;
      if (drift < 1e-9) ++invariantHeld;
      else if (drift < 0.01) ++invariantDriftSmall;
      else {
        ++invariantDriftLarge;
        drifters.emplace_back(c.source, drift, Au, Ai, rhs);
      }
    } catch (const std::exception& ex) {
      std::cerr << "  EXC " << c.source << ": " << ex.what() << "\n";
      ++topoFail;
    }
  }
  std::cout << "  overlayAreaTest run:   " << areaTests << "\n";
  std::cout << "  unionArea run:         " << unionAreaTests << "\n";
  std::cout << "  Skipped (other ops):   " << skipped << "\n\n";
  std::cout << "  Topology valid:   " << topoOk << "\n";
  std::cout << "  Topology INVALID: " << topoFail << "\n\n";
  std::cout << "  Invariant holds (drift < 1e-9):  " << invariantHeld
            << " / " << areaTests << "\n";
  std::cout << "  Invariant near-holds (< 1%):     " << invariantDriftSmall
            << " / " << areaTests << "\n";
  std::cout << "  Invariant violated (>= 1%):      " << invariantDriftLarge
            << " / " << areaTests << "\n";
  if (unionAreaTests > 0) {
    std::cout << "  unionArea matches expected:   " << unionAreaOk << " / "
              << unionAreaTests << "\n";
    std::cout << "  unionArea drifts:             " << unionAreaDrift << " / "
              << unionAreaTests << "\n";
  }
  if (!drifters.empty()) {
    std::sort(drifters.begin(), drifters.end(),
              [](const auto& a, const auto& b) {
                return std::get<1>(a) > std::get<1>(b);
              });
    std::cout << "\n  Invariant violators (worst first, up to 20):\n";
    for (size_t i = 0; i < drifters.size() && i < 20; ++i) {
      auto& [src, drift, Au, Ai, rhs] = drifters[i];
      std::cout << "    " << src
                << "  drift=" << (drift * 100.0) << "%"
                << "  A∪B=" << Au << "  A∩B=" << Ai << "  |A|+|B|=" << rhs
                << "\n";
    }
  }
}

// =============================================================================
// Generated CAD-ish corpus.
//
// JTS is intentionally GIS-shaped. This lane gives the perf harness a compact,
// deterministic mix closer to Manifold CrossSection inputs: clean overlapping
// primitives, smooth polygon approximations, orthogonal notches, radial parts,
// and a frame with a hole. It is not a correctness oracle; it is a timing lane
// that exercises cleaner modeling-style geometry at several vertex counts.
// =============================================================================
struct CadCase {
  std::string name;
  manifold::Polygons a;
  manifold::Polygons b;
};

inline manifold::SimplePolygon Rect(double x0, double y0, double x1,
                                    double y1) {
  return {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
}

inline manifold::SimplePolygon RegularPoly(int n, double cx, double cy,
                                           double r, double phase = 0.0) {
  manifold::SimplePolygon out;
  out.reserve(n);
  for (int i = 0; i < n; ++i) {
    const double a = phase + 2.0 * M_PI * i / n;
    out.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
  }
  return out;
}

inline manifold::SimplePolygon RadialPart(int n, double cx, double cy,
                                          double r0, double r1,
                                          double phase = 0.0) {
  manifold::SimplePolygon out;
  out.reserve(n);
  for (int i = 0; i < n; ++i) {
    const double a = phase + 2.0 * M_PI * i / n;
    const double t = 0.5 + 0.5 * std::sin(5.0 * a) * std::cos(3.0 * a);
    const double r = r0 + (r1 - r0) * t;
    out.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
  }
  return out;
}

inline manifold::SimplePolygon NotchedRect(double x0, double y0, double x1,
                                           double y1, double notch) {
  const double mx = 0.5 * (x0 + x1);
  return {{x0, y0}, {x1, y0}, {x1, y1}, {mx + notch, y1},
          {mx + notch, y1 - notch}, {mx - notch, y1 - notch},
          {mx - notch, y1}, {x0, y1}};
}

inline manifold::Polygons Frame(double outer, double inner) {
  auto out = RegularPoly(4, 0.0, 0.0, outer, M_PI / 4.0);
  auto hole = RegularPoly(4, 0.0, 0.0, inner, M_PI / 4.0);
  std::reverse(hole.begin(), hole.end());
  return {out, hole};
}

#ifdef BOOLEAN2_WITH_MANIFOLD
inline bool HasAnyLoop(const manifold::Polygons& polygons) {
  for (const auto& loop : polygons) {
    if (loop.size() >= 3) return true;
  }
  return false;
}

inline void AddCadCase(std::vector<CadCase>* cases, std::string name,
                       manifold::Polygons a, manifold::Polygons b) {
  if (!HasAnyLoop(a) || !HasAnyLoop(b)) return;
  cases->push_back({std::move(name), std::move(a), std::move(b)});
}

inline void AddProjectPair(std::vector<CadCase>* cases, const std::string& name,
                           const Manifold& a, const Manifold& b) {
  AddCadCase(cases, name + " project", a.Project(), b.Project());
}

inline void AddSlicePair(std::vector<CadCase>* cases, const std::string& name,
                         const Manifold& a, const Manifold& b,
                         std::initializer_list<double> heights) {
  for (double z : heights) {
    std::ostringstream label;
    label << name << " slice z=" << z;
    AddCadCase(cases, label.str(), a.Slice(z), b.Slice(z));
  }
}

inline Manifold SlottedBlock() {
  Manifold block = Manifold::Cube({1.45, 1.2, 1.0}, true);
  const Manifold slotX =
      Manifold::Cylinder(1.8, 0.18, -1.0, 48, true).Rotate(0, 90, 0);
  const Manifold slotY =
      Manifold::Cylinder(1.6, 0.13, -1.0, 40, true).Rotate(90, 0, 0);
  block -= slotX.Translate({0.0, 0.28, 0.0});
  block -= slotX.Translate({0.0, -0.28, 0.0});
  block -= slotY.Translate({0.32, 0.0, 0.18});
  block -= slotY.Translate({-0.32, 0.0, -0.18});
  return block;
}

inline Manifold TwistedStarPrism(int n, double phase) {
  return Manifold::Extrude({RadialPart(n, 0.0, 0.0, 0.42, 0.78, phase)},
                           1.2, 10, 35.0, {0.72, 1.08})
      .Translate({0.0, 0.0, -0.6});
}

inline Manifold GyroidPatch() {
  return Manifold::LevelSet(
      [](vec3 p) {
        return std::sin(4.0 * p.x) * std::cos(4.0 * p.y) +
               std::sin(4.0 * p.y) * std::cos(4.0 * p.z) +
               std::sin(4.0 * p.z) * std::cos(4.0 * p.x) - 0.15;
      },
      Box({-0.85, -0.85, -0.85}, {0.85, 0.85, 0.85}), 0.14);
}

inline void AppendManifoldDerivedCadCases(std::vector<CadCase>* cases) {
  const Manifold sphere = Manifold::Sphere(0.72, 48);
  const Manifold drilled =
      Manifold::Cube(vec3(1.4), true) -
      Manifold::Cylinder(1.7, 0.23, -1.0, 40, true).Rotate(90, 0, 0) -
      Manifold::Cylinder(1.7, 0.18, -1.0, 36, true).Rotate(0, 90, 0);
  AddCadCase(cases, "drilled cube project vs circle", drilled.Project(),
             {RegularPoly(96, 0.12, -0.04, 0.58)});
  AddSlicePair(cases, "drilled cube vs sphere", drilled,
               sphere.Translate({0.08, -0.06, 0.0}),
               {-0.34, -0.16, 0.0, 0.16, 0.34});

  const Manifold slotted = SlottedBlock();
  const Manifold crossingCylinder =
      Manifold::Cylinder(1.7, 0.16, -1.0, 64, true).Rotate(90, 0, 0);
  AddProjectPair(cases, "slotted block vs crossing cylinder", slotted,
                 crossingCylinder.Translate({0.18, 0.0, 0.05}));
  AddSlicePair(cases, "slotted block vs crossing cylinder", slotted,
               crossingCylinder.Translate({0.18, 0.0, 0.05}),
               {-0.32, -0.12, 0.12, 0.32});

  const Manifold starA = TwistedStarPrism(72, 0.0);
  const Manifold starB =
      TwistedStarPrism(96, M_PI / 11.0).Rotate(0, 0, 18).Translate(
          {0.14, -0.08, 0.0});
  AddProjectPair(cases, "twisted radial prisms", starA, starB);
  AddSlicePair(cases, "twisted radial prisms", starA, starB,
               {-0.42, -0.18, 0.0, 0.18, 0.42});

  const Manifold cylinder =
      Manifold::Cylinder(1.8, 0.46, 0.22, 96, true).Rotate(18, 0, 0);
  AddProjectPair(cases, "cone-cylinder vs sphere", cylinder,
                 sphere.Translate({0.12, 0.02, 0.0}));
  AddSlicePair(cases, "cone-cylinder vs slotted block", cylinder, slotted,
               {-0.36, -0.12, 0.12, 0.36});

  const Manifold gyroid = GyroidPatch();
  AddProjectPair(cases, "gyroid patch vs sphere", gyroid,
                 sphere.Translate({0.05, -0.04, 0.0}));
  AddSlicePair(cases, "gyroid patch vs crossing cylinder", gyroid,
               crossingCylinder,
               {-0.36, -0.18, 0.0, 0.18, 0.36});
}
#else
inline void AppendManifoldDerivedCadCases(std::vector<CadCase>*) {}
#endif

inline std::vector<CadCase> GenerateCadCases() {
  std::vector<CadCase> cases;
  cases.push_back({"overlapping boxes",
                   {Rect(-1.0, -0.7, 1.0, 0.7)},
                   {Rect(-0.2, -1.0, 1.4, 0.55)}});
  cases.push_back({"orthogonal notches",
                   {NotchedRect(-1.4, -1.0, 1.2, 1.0, 0.35)},
                   {NotchedRect(-0.7, -1.2, 1.5, 0.9, 0.25)}});
  for (int n : {24, 96, 256}) {
    cases.push_back({"regular " + std::to_string(n),
                     {RegularPoly(n, -0.15, 0.0, 1.0)},
                     {RegularPoly(n, 0.22, 0.08, 0.92, M_PI / n)}});
    cases.push_back({"radial " + std::to_string(n),
                     {RadialPart(n, -0.1, 0.0, 0.65, 1.15)},
                     {RadialPart(n, 0.18, 0.05, 0.55, 1.05, M_PI / n)}});
  }
  cases.push_back({"frame vs bar", Frame(1.2, 0.45),
                   {Rect(-1.5, -0.25, 1.5, 0.25)}});
  AppendManifoldDerivedCadCases(&cases);
  return cases;
}

inline void RunCadCorpus() {
  const auto cases = GenerateCadCases();
  std::cout << "=== Generated CAD corpus: " << cases.size() << " cases ===\n";
  size_t outputs = 0;
  double totalArea = 0.0;
  for (const auto& c : cases) {
    const double eps = InferEps(c.a, c.b);
    auto add = Boolean2D(c.a, c.b, OpType::Add, eps);
    auto sub = Boolean2D(c.a, c.b, OpType::Subtract, eps);
    auto isec = Boolean2D(c.a, c.b, OpType::Intersect, eps);
    outputs += add.size() + sub.size() + isec.size();
    totalArea += std::fabs(TotalSignedArea(add)) +
                 std::fabs(TotalSignedArea(sub)) +
                 std::fabs(TotalSignedArea(isec));
  }
  std::cout << "  Boolean ops: " << (cases.size() * 3) << "\n";
  std::cout << "  Output loops: " << outputs << "\n";
  std::cout << "  Aggregate abs area: " << totalArea << "\n";
}

}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// Head-to-head benchmark vs Clipper2 (the library boolean2 is intended to
// eventually replace as `CrossSection`'s boolean/Simplify backend). Built
// only when -DBOOLEAN2_WITH_CLIPPER2 is defined and Clipper2 is linked.
//
// Compile/link:
//   g++ ... -DBOOLEAN2_WITH_CLIPPER2 -DMANIFOLD_PAR=1 \
//     -I build/_deps/clipper2-src/CPP/Clipper2Lib/include \
//     extras/boolean2_driver.cpp \
//     build/_deps/clipper2-build/libClipper2.a -ltbb \
//     -o boolean2_driver
//
// What's timed: the boolean call itself, on already-loaded inputs. We
// exclude geometry conversion (manifold::Polygons → Clipper2 PathsD)
// because that's a one-shot cost both sides would amortize out under
// any real-world wiring.
//
// Workloads: each of the three corpora we already validate against,
// run a UNION on every (sub, clip) pair (or every subject-only case).
// Clipper2's NONZERO is the cross-comparable fill rule; cases with
// EVENODD or NEGATIVE/POSITIVE are skipped here so both libraries do
// the same work. (Correctness comparison happens in the corpus runs;
// this benchmark is purely about wall-clock at parity.)
//
// Reference: `Clipper64` (native int64 API), not `ClipperD` (real-coord
// wrapper). Per elalish/manifold#1683 the planned post-precision-bug
// `CrossSection` wiring is `Paths64`-based, so `Clipper64` is the right
// baseline; `ClipperD`'s per-op rescaling would flatter our number
// against an end state nobody plans to ship.
// =============================================================================
#ifdef BOOLEAN2_WITH_CLIPPER2
#include "clipper2/clipper.h"
namespace manifold {
namespace boolean2 {

inline Clipper2Lib::PathsD ToPathsD(const manifold::Polygons& p) {
  Clipper2Lib::PathsD out;
  out.reserve(p.size());
  for (const auto& loop : p) {
    Clipper2Lib::PathD path;
    path.reserve(loop.size());
    for (const auto& v : loop) path.emplace_back(v.x, v.y);
    out.push_back(std::move(path));
  }
  return out;
}

// Scale doubles to int64 with a fixed factor and convert to Paths64 (the
// native Clipper2 API; ClipperD wraps this with rescaling at every op).
// Per elalish/manifold#1683 the planned post-Clipper2-precision-bug
// CrossSection wiring is Paths64-based, so the head-to-head we want
// reports against this faster, full-int-precision path. Scale chosen
// so a 1e6-coordinate input still has ~1 unit precision in int64.
constexpr double kPaths64Scale = 1e9;
inline Clipper2Lib::Paths64 ToPaths64(const manifold::Polygons& p) {
  Clipper2Lib::Paths64 out;
  out.reserve(p.size());
  for (const auto& loop : p) {
    Clipper2Lib::Path64 path;
    path.reserve(loop.size());
    for (const auto& v : loop) {
      path.emplace_back(static_cast<int64_t>(v.x * kPaths64Scale),
                        static_cast<int64_t>(v.y * kPaths64Scale));
    }
    out.push_back(std::move(path));
  }
  return out;
}

struct BenchTotals {
  int64_t oursNs = 0;
  int64_t clipperNs = 0;
  int cases = 0;
  int oursDrops = 0;     // we returned empty
  int clipperDrops = 0;  // clipper2 returned empty
};

inline void TimeOneUnion(const manifold::Polygons& a,
                         const manifold::Polygons& b, BenchTotals* totals) {
  using Clock = std::chrono::steady_clock;
  const double eps = InferEps(a, b);
  // Ours: single-pass dispatch (mfogel-runner shape). Combine subject
  // and clip edges with mult=+1 each, run RemoveOverlaps2D once with
  // a NONZERO `w != 0` predicate - matches Clipper2's UNION+NONZERO
  // semantics directly without paying for two extra FillUnderRule
  // pipeline runs that the multi-pass Boolean2D path would incur.
  std::vector<vec2> verts;
  std::vector<EdgeM> edges;
  auto append = [&](const manifold::Polygons& polys) {
    for (const auto& loop : polys) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i)
        edges.push_back({base + i, base + ((i + 1) % n), 1});
    }
  };
  manifold::Polygons ours;
  auto t0 = Clock::now();
  append(a);
  append(b);
  if (!verts.empty()) {
    auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                              [](int w) { return w != 0; });
    ours = OutEdgesToPolygons(r.verts, r.edges);
  }
  auto t1 = Clock::now();
  totals->oursNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        t1 - t0).count();
  if (ours.empty()) ++totals->oursDrops;

  // Clipper2 Clipper64 (native int64 API; what CrossSection is targeting
  // post-#1683 to avoid PathsD's precision cap and rescale overhead).
  // Conversion to Paths64 is excluded from the timed region; we only
  // time Execute, the same as for boolean2.
  auto a2 = ToPaths64(a);
  auto b2 = ToPaths64(b);
  Clipper2Lib::Clipper64 clip;
  clip.AddSubject(a2);
  if (!b.empty()) clip.AddClip(b2);
  Clipper2Lib::Paths64 sol;
  auto t2 = Clock::now();
  clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero,
               sol);
  auto t3 = Clock::now();
  totals->clipperNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           t3 - t2).count();
  if (sol.empty()) ++totals->clipperDrops;
  ++totals->cases;
}

inline void ReportBench(const std::string& label, const BenchTotals& t) {
  auto fmt = [](int64_t ns) {
    std::ostringstream o; o.setf(std::ios::fixed); o.precision(3);
    o << ns * 1e-6 << " ms";
    return o.str();
  };
  std::cout << "  " << label << ":\n";
  std::cout << "    cases:        " << t.cases << "\n";
  std::cout << "    boolean2:    " << fmt(t.oursNs)
            << " (" << (t.cases ? t.oursNs / t.cases / 1000.0 : 0)
            << " µs/case avg)";
  if (t.oursDrops) std::cout << "  drops=" << t.oursDrops;
  std::cout << "\n";
  std::cout << "    Clipper2 :    " << fmt(t.clipperNs)
            << " (" << (t.cases ? t.clipperNs / t.cases / 1000.0 : 0)
            << " µs/case avg)";
  if (t.clipperDrops) std::cout << "  drops=" << t.clipperDrops;
  std::cout << "\n";
  if (t.clipperNs > 0) {
    const double ratio = static_cast<double>(t.oursNs) / t.clipperNs;
    std::cout << "    ratio:        ";
    std::cout.setf(std::ios::fixed); std::cout.precision(2);
    std::cout << ratio << "x  ("
              << (ratio < 1 ? "boolean2 faster" : "Clipper2 faster")
              << ")\n";
    std::cout.unsetf(std::ios::fixed);
  }
}

inline void RunVsClipper2_Clipper2Corpus(const std::string& path) {
  auto cases = LoadClipper2Cases(path);
  BenchTotals t;
  for (const auto& c : cases) {
    if (c.fillrule != "NONZERO" || c.cliptype != "UNION") continue;
    TimeOneUnion(c.subjects, c.clips, &t);
  }
  ReportBench("Clipper2 corpus / UNION+NONZERO", t);
}

inline void RunVsClipper2_CadCorpus() {
  auto cases = GenerateCadCases();
  BenchTotals t;
  for (const auto& c : cases) {
    TimeOneUnion(c.a, c.b, &t);
  }
  ReportBench("Generated CAD corpus / UNION", t);
}

// Per-case offset benchmark + comparison. Clipper2's Offsets.txt doesn't
// encode per-case delta / jointype, so we pick a sensible default: 1%
// inset (negative delta) with Round join. Inset is the harder direction
// (concave features can collapse), and the relative delta scales with
// input bbox so we exercise both the unit-scale and the GIS-scale (~1e3)
// cases in the same way.
inline void TimeOneOffset(const manifold::Polygons& subjects, double delta,
                          manifold::boolean2::JoinType jt, double miterLimit,
                          double arcTol, BenchTotals* totals) {
  using Clock = std::chrono::steady_clock;
  auto t0 = Clock::now();
  auto ours = manifold::boolean2::Offset(subjects, delta, jt, miterLimit,
                                         arcTol);
  auto t1 = Clock::now();
  totals->oursNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        t1 - t0).count();
  if (ours.empty()) ++totals->oursDrops;
#ifdef BOOLEAN2_WITH_CLIPPER2
  auto s64 = ToPaths64(subjects);
  Clipper2Lib::JoinType cjt;
  switch (jt) {
    case manifold::boolean2::JoinType::Round:
      cjt = Clipper2Lib::JoinType::Round;
      break;
    case manifold::boolean2::JoinType::Miter:
      cjt = Clipper2Lib::JoinType::Miter;
      break;
    case manifold::boolean2::JoinType::Bevel:
      cjt = Clipper2Lib::JoinType::Bevel;
      break;
    case manifold::boolean2::JoinType::Square:
      cjt = Clipper2Lib::JoinType::Square;
      break;
  }
  auto t2 = Clock::now();
  auto sol =
      Clipper2Lib::InflatePaths(s64, delta, cjt, Clipper2Lib::EndType::Polygon,
                                miterLimit, arcTol);
  auto t3 = Clock::now();
  totals->clipperNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
                           t3 - t2).count();
  if (sol.empty()) ++totals->clipperDrops;
#endif
  ++totals->cases;
}

inline void RunVsClipper2_OffsetsCorpus(const std::string& path) {
  auto cases = LoadClipper2Cases(path);
  if (cases.empty()) return;
  BenchTotals t;
  for (const auto& c : cases) {
    if (c.cliptype != "OFFSET") continue;
    double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300;
    for (const auto& loop : c.subjects) {
      for (const auto& v : loop) {
        xmin = std::min(xmin, v.x);
        ymin = std::min(ymin, v.y);
        xmax = std::max(xmax, v.x);
        ymax = std::max(ymax, v.y);
      }
    }
    if (xmax < xmin) continue;  // empty case
    const double extent = std::max(xmax - xmin, ymax - ymin);
    if (extent <= 0) continue;
    const double delta = -extent / 100.0;
    const double absDelta = std::fabs(delta);
    TimeOneOffset(c.subjects, delta, manifold::boolean2::JoinType::Round,
                  /*miterLimit=*/2.0, /*arcTol=*/absDelta / 100.0, &t);
  }
  ReportBench("Clipper2 offsets corpus / Round / -1% inset", t);
}

// Case-insensitive ASCII compare for op-tag dispatch.
inline bool OpEq(const std::string& a, const char* lit) {
  size_t n = std::strlen(lit);
  if (a.size() != n) return false;
  for (size_t i = 0; i < n; ++i) {
    char ca = a[i];
    char cb = lit[i];
    if (ca >= 'A' && ca <= 'Z') ca = static_cast<char>(ca - 'A' + 'a');
    if (cb >= 'A' && cb <= 'Z') cb = static_cast<char>(cb - 'A' + 'a');
    if (ca != cb) return false;
  }
  return true;
}

inline void RunVsClipper2_JtsCorpus(const std::string& path,
                                    const char* label = "JTS overlay") {
  auto cases = LoadJtsCorpus(path);
  BenchTotals t;
  const bool perCase = std::getenv("BOOLEAN2_VERBOSE_CORPUS") != nullptr;
  for (const auto& c : cases) {
    if (perCase) {
      std::fprintf(stderr, "  TimeOneUnion: %s\n", c.source.c_str());
      std::fflush(stderr);
    }
    if (OpEq(c.op, "overlayAreaTest")) {
      TimeOneUnion(c.a, c.b, &t);
    } else if (OpEq(c.op, "unionArea")) {
      TimeOneUnion(c.a, {}, &t);
    }
  }
  std::string heading = std::string(label) + " / UNION";
  ReportBench(heading, t);
}

// Manifold-derived corpus (Project / Slice outputs from representative
// 3D meshes; capture lives in extras/capture_manifold_corpus.cpp). Same
// file format as the JTS corpus, so we re-use the loader. Used to
// regression-test boolean2 against the polygon distributions that
// arise from manifold's 3D-to-2D operations.
//
// We also run a *verbose* pass under WindRule::Add / Clipper2's
// Positive fill rule, mirroring what CrossSection's constructors
// (default `FillRule::Positive`) actually call. The driver's
// canonical `TimeOneUnion` uses NonZero, which doesn't expose the
// known-failure cases where Positive returns empty while Clipper2
// produces a non-empty boundary.
inline void RunVsClipper2_ManifoldCorpus(const std::string& path) {
  RunVsClipper2_JtsCorpus(path, "Manifold corpus");

  auto cases = LoadJtsCorpus(path);
  std::cout << "\n  Manifold corpus / Positive (CrossSection ctor "
               "semantics):\n";
  size_t pass = 0, oursEmptyTheyNot = 0, theirsEmptyOursNot = 0,
         areaMismatch = 0;
  for (const auto& c : cases) {
    if (!OpEq(c.op, "unionArea")) continue;
    const double eps = InferEps(c.a, {});
    // Run boolean2 (FillByRule equivalent: WindRule::Add).
    std::vector<vec2> verts;
    std::vector<EdgeM> edges;
    for (const auto& loop : c.a) {
      if (loop.size() < 3) continue;
      const int base = static_cast<int>(verts.size());
      const int n = static_cast<int>(loop.size());
      for (const auto& v : loop) verts.push_back(v);
      for (int i = 0; i < n; ++i)
        edges.push_back({base + i, base + ((i + 1) % n), 1});
    }
    manifold::Polygons ours;
    if (!verts.empty()) {
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                                [](int w) { return w > 0; });
      ours = OutEdgesToPolygons(r.verts, r.edges);
    }
    auto a2 = ToPaths64(c.a);
    Clipper2Lib::Clipper64 clip;
    clip.AddSubject(a2);
    Clipper2Lib::Paths64 sol;
    clip.Execute(Clipper2Lib::ClipType::Union,
                 Clipper2Lib::FillRule::Positive, sol);
    const bool oursEmpty = ours.empty();
    const bool theirsEmpty = sol.empty();
    double oursArea = 0.0;
    for (const auto& l : ours) {
      const auto& r = l.front();
      for (size_t i = 0; i < l.size(); ++i) {
        const auto& v0 = l[i];
        const auto& v1 = l[(i + 1) % l.size()];
        oursArea +=
            0.5 * ((v0.x - r.x) * (v1.y - r.y) - (v1.x - r.x) * (v0.y - r.y));
      }
    }
    double theirsArea = 0.0;
    for (const auto& p : sol) {
      theirsArea += Clipper2Lib::Area(p);
    }
    // Clipper2 area is in scaled-int64 units; undo the scale^2.
    theirsArea /= (kPaths64Scale * kPaths64Scale);
    const char* status = "PASS";
    if (oursEmpty && !theirsEmpty) {
      status = "FAIL ours_empty";
      ++oursEmptyTheyNot;
    } else if (!oursEmpty && theirsEmpty) {
      status = "FAIL theirs_empty";
      ++theirsEmptyOursNot;
    } else if (std::fabs(oursArea - theirsArea) >
               1e-6 * std::max(std::fabs(oursArea), std::fabs(theirsArea))) {
      status = "FAIL area_mismatch";
      ++areaMismatch;
    } else {
      ++pass;
    }
    std::cout << "    [" << status << "] " << c.source
              << "  ours=" << ours.size() << "/" << oursArea
              << " theirs=" << sol.size() << "/" << theirsArea << "\n";
  }
  std::cout << "  pass=" << pass << " fail=" << (oursEmptyTheyNot + theirsEmptyOursNot + areaMismatch)
            << " (ours_empty=" << oursEmptyTheyNot
            << " theirs_empty=" << theirsEmptyOursNot
            << " area_mismatch=" << areaMismatch << ")\n";
}

inline void RunVsClipper2_MfogelCorpus(const std::string& dir) {
  auto cases = LoadMfogelCorpus(dir);
  BenchTotals t;
  for (const auto& c : cases) {
    // Run a UNION on the first two features (the most common 2-input
    // shape; subject-only single-feature cases are skipped here since
    // that's a self-clean and Clipper2's fill-rule handling differs).
    if (c.features.size() >= 2 && !c.features[0].empty() &&
        !c.features[1].empty()) {
      TimeOneUnion(c.features[0], c.features[1], &t);
    }
  }
  ReportBench("mfogel end-to-end / UNION", t);
}

// Return `primary` if it exists, otherwise `fallback`. The Clipper2 test
// corpora live in `test/polygons/clipper2_*.txt` after running
// `scripts/fetch_clipper2_test_data.sh`; the fallback path is the
// historical location under a sister Clipper2 build.
inline std::string CorpusPath(const std::string& primary,
                              const std::string& fallback) {
  std::ifstream a(primary);
  if (a) return primary;
  std::ifstream b(fallback);
  if (b) return fallback;
  return primary;  // surface the canonical path in the error message
}

inline void RunVsClipper2(const std::string& which) {
  std::cout << "=== boolean2 vs Clipper2 head-to-head ===\n\n";
  if (which == "all" || which == "clipper2") {
    RunVsClipper2_Clipper2Corpus(CorpusPath(
        "test/polygons/clipper2_polygons.txt",
        "build/_deps/clipper2-src/Tests/Polygons.txt"));
  }
  if (which == "all" || which == "cad") {
    RunVsClipper2_CadCorpus();
  }
  if (which == "all" || which == "offsets") {
    RunVsClipper2_OffsetsCorpus(CorpusPath(
        "test/polygons/clipper2_offsets.txt",
        "build/_deps/clipper2-src/Tests/Offsets.txt"));
  }
  if (which == "all" || which == "jts") {
    RunVsClipper2_JtsCorpus("test/polygons/jts_overlay_corpus.txt");
  }
  if (which == "all" || which == "mfogel") {
    RunVsClipper2_MfogelCorpus(
        "build/_deps/mfogel-polygon-clipping-src/test/end-to-end");
  }
  if (which == "all" || which == "manifoldcorpus" || which == "manifold") {
    RunVsClipper2_ManifoldCorpus("test/polygons/manifold_corpus.txt");
  }
}

// Diagnostic: walk the JTS corpus and print every case where boolean2's
// output diverges from Clipper2's (one drops while the other doesn't, or
// area differs by more than a small tolerance). Used to triage drop cases
// without polluting the timed benchmark path.
inline void RunJtsDropDiag(const std::string& path) {
  auto cases = LoadJtsCorpus(path);
  std::cout << "=== JTS drop / divergence diagnostic ===\n";
  int oursDrops = 0, clipperDrops = 0, bothEmpty = 0, bothNonEmpty = 0;
  int ourEmptyTheyNot = 0, theirEmptyOursNot = 0;
  for (const auto& c : cases) {
    bool isUnion = (c.op == "overlayAreaTest" || c.op == "overlayareatest");
    bool isUnionArea = (c.op == "unionArea");
    if (!isUnion && !isUnionArea) continue;
    const manifold::Polygons& A = c.a;
    const manifold::Polygons& B = isUnion ? c.b : manifold::Polygons{};
    const double eps = InferEps(A, B);
    // Build inputs for boolean2.
    std::vector<vec2> verts;
    std::vector<EdgeM> edges;
    auto append = [&](const manifold::Polygons& polys) {
      for (const auto& loop : polys) {
        if (loop.size() < 3) continue;
        const int base = static_cast<int>(verts.size());
        const int n = static_cast<int>(loop.size());
        for (const auto& v : loop) verts.push_back(v);
        for (int i = 0; i < n; ++i)
          edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
    };
    append(A); append(B);
    manifold::Polygons ours;
    if (!verts.empty()) {
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                                [](int w) { return w != 0; });
      ours = OutEdgesToPolygons(r.verts, r.edges);
    }
    // Clipper2 reference.
    auto a2 = ToPaths64(A);
    auto b2 = ToPaths64(B);
    Clipper2Lib::Clipper64 clip;
    clip.AddSubject(a2);
    if (!B.empty()) clip.AddClip(b2);
    Clipper2Lib::Paths64 sol;
    clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero,
                 sol);
    const bool oursEmpty = ours.empty();
    const bool theirsEmpty = sol.empty();
    if (oursEmpty) ++oursDrops;
    if (theirsEmpty) ++clipperDrops;
    if (oursEmpty && theirsEmpty) ++bothEmpty;
    if (!oursEmpty && !theirsEmpty) ++bothNonEmpty;
    if (oursEmpty && !theirsEmpty) {
      ++ourEmptyTheyNot;
      // Clipper2 area for context.
      double clipArea = 0;
      for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      // Total input vert/edge count for size context.
      int totalVerts = 0, totalRings = 0;
      for (const auto& p : A) { totalVerts += p.size(); ++totalRings; }
      for (const auto& p : B) { totalVerts += p.size(); ++totalRings; }
      std::cout << "OURS-DROP case=" << c.n << " op=" << c.op
                << " src=" << c.source
                << " A.rings=" << A.size() << " B.rings=" << B.size()
                << " verts=" << totalVerts
                << " eps=" << eps
                << " clipperArea=" << clipArea << "\n";
    }
    if (!oursEmpty && theirsEmpty) {
      ++theirEmptyOursNot;
      std::cout << "THEY-DROP case=" << c.n << " op=" << c.op
                << " src=" << c.source << "\n";
    }
  }
  std::cout << "\nSummary:\n";
  std::cout << "  total cases:           " << cases.size() << "\n";
  std::cout << "  both nonempty:         " << bothNonEmpty << "\n";
  std::cout << "  both empty:            " << bothEmpty << "\n";
  std::cout << "  ours drops, theirs ok: " << ourEmptyTheyNot << "\n";
  std::cout << "  theirs drops, ours ok: " << theirEmptyOursNot << "\n";
  std::cout << "  total ours drops:      " << oursDrops << "\n";
  std::cout << "  total theirs drops:    " << clipperDrops << "\n";
}

}  // namespace boolean2
}  // namespace manifold
#endif  // BOOLEAN2_WITH_CLIPPER2

// =============================================================================
// Deep fuzz: wider seed sweep + iter-vs-topology verification.
//
// The standard displacement fuzz uses 30 seeds per (kPow, n) cell. This
// mode runs `seedsPerCell` seeds per cell across a wider parameter grid
// and, for every case the iterate-to-fixed-point reports as iter ≥ 2,
// compares pass-1 vs pass-2 outputs at:
//
//   - eps/100 quantum (the standard fingerprint, what triggers iter ≥ 2)
//   - eps quantum (topology-only; captures sub-edges modulo sub-eps drift)
//
// If pass-1 and pass-2 differ at eps/100 but match at eps for every iter≥2
// case, the iteration is shifting positions only, not topology. That's the
// cut-corner answer for whether iter=2 is doing real work.
// =============================================================================
namespace manifold {
namespace boolean2 {
void DeepFuzz(int seedsPerCell) {
  const std::vector<int> kPows = {10, 20, 30, 35, 40, 45, 49};
  const std::vector<int> sizes = {8, 20, 50, 100};
  std::cout << "=== DeepFuzz: " << seedsPerCell << " seeds × " << sizes.size()
            << " sizes × " << kPows.size() << " kPow = "
            << (seedsPerCell * sizes.size() * kPows.size()) << " cases ===\n";

  int total = 0;
  int firstPassFail = 0;
  int convergedPassValid = 0;     // converged result: topo valid
  int convergedPassInvalid = 0;   // converged result: topo invalid
  // Polygons-API regularization counter. Each fuzz case is run through
  // both the lower-level pipeline (counts every retained sub-edge) and
  // the public Simplify (regularizes via OutEdgesToPolygons, dropping
  // zero-area lens loops). When the counts disagree, the difference is
  // exactly the regularization drop, which is correct CGAL/Clipper2/SVG
  // behavior. We track the rate as a regularization-frequency stat
  // rather than as a "mismatch" (it isn't a bug), and as a tripwire
  // for any future drift in regularization behavior.
  int regUnchanged = 0;       // post-regularization edge count == lower-level
  int regDropped = 0;         // dropped one or more zero-area lens loops
  std::map<int, int> iterDist;
  int iterGE2 = 0;
  int iterGE2_topo_match = 0;     // pass1 == pass2 at eps quantum
  int iterGE2_fine_match = 0;     // pass1 == pass2 at eps/100 quantum
  int iterGE2_topo_diff = 0;
  int iterGE2_pass1_invalid = 0;  // iter≥2 cases where pass 1 was wrong
  int iterGE2_pass1_valid = 0;    // iter≥2 cases where pass 1 was already valid
  int iter_repaired = 0;          // pass-1 invalid → converged valid
  int iter_degraded = 0;          // pass-1 valid → converged INVALID
  int pass1_invalid_unfixed = 0;  // pass-1 invalid AND converged invalid
  // Geometric-emptiness check: counts cases where pass 2 outputs zero edges
  // even though pass 1 output non-zero. Indicates pass-1 → pass-2 winding-
  // sign disagreement (output orientation not Smith-convention-compatible).
  int pass2_collapsed = 0;
  int pass1_nonempty = 0;
  // Area-drift tracking: a quantitative oracle on top of the topology
  // check. Pass 1 produces an output area; iterate-to-fixed-point
  // produces another. If they differ by > 1% we flag it as a regression
  // (per the design doc's quantitative oracle recommendation). Catches
  // thin-polygon-style collapses where iteration empties a non-empty
  // pass-1 output, which topology balance alone misses.
  int area_pass1_nonzero = 0;       // pass-1 had a measurable area
  int area_drift_over_1pct = 0;     // > 1% drift between pass 1 and converged
  int area_drift_over_10pct = 0;    // > 10% drift (typically a collapse)
  double area_drift_max = 0.0;      // largest fractional drift seen
  std::vector<std::tuple<int, int, uint64_t, double>> areaDriftList;
  std::vector<std::tuple<int, int, uint64_t>> topoMismatch;
  std::vector<std::tuple<int, int, uint64_t, int>> convergedInvalidList;
  std::vector<std::tuple<int, int, uint64_t, int>> degradedList;
  // (kPow, n, seed, iters) for both
  std::vector<std::tuple<int, int, uint64_t>> geomCollapseList;

  for (int kPow : kPows) {
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    for (int n : sizes) {
      for (int seed = 0; seed < seedsPerCell; ++seed) {
        ++total;
        // Seed mapping matches the standard displacement fuzz in main()
        // so any failing (kPow, n, seed) here can be passed verbatim to
        // `./boolean2_driver diagnose <seed> <kPow> <n>`.
        auto [v, e] = RandomTopologicalPolygon(
            n, static_cast<uint64_t>(seed) + 1000ull * kPow);
        Displace(&v, offset);

        auto pass1 = RemoveOverlaps2D(v, e, eps);
        const bool pass1Valid = CheckTopologicalValidity(
            pass1, e, pass1.inputRemap, pass1.numMergedVerts);
        if (!pass1Valid) ++firstPassFail;

        // Polygons round-trip: convert (v, e) (which is a single closed
        // cycle from RandomTopologicalPolygon) into manifold::Polygons,
        // run the public Simplify API, and check the
        // boundary-edge count matches the lower-level RemoveOverlaps2D
        // output. Catches regressions in PolygonsToInput's drop-degenerate
        // logic and OutEdgesToPolygons's loop-walking.
        {
          // Walk the edge cycle (v0 -> v1 chain) into a single SimplePolygon.
          std::vector<int> nextV(v.size(), -1);
          for (const auto& edge : e) nextV[edge.v0] = edge.v1;
          manifold::SimplePolygon loop;
          int cur = e.empty() ? -1 : e[0].v0;
          for (size_t walk = 0; walk < v.size() && cur >= 0; ++walk) {
            loop.push_back(v[cur]);
            cur = nextV[cur];
            if (cur == e[0].v0) break;
          }
          manifold::Polygons polys = {std::move(loop)};
          auto out = Simplify(polys, eps);
          // Compare lower-level edge count to post-regularization count.
          // A drop indicates OutEdgesToPolygons regularized one or more
          // zero-area lens loops out of the result, which is correct
          // CGAL/Clipper2/SVG behavior; track the rate as a stat, not
          // a bug.
          size_t outEdges = 0;
          for (const auto& l : out) outEdges += l.size();
          if (outEdges == pass1.edges.size())
            ++regUnchanged;
          else
            ++regDropped;
        }

        // Run iterate-to-fixed-point to classify.
        int iters = 0;
        IterStatus status = IterStatus::Converged;
        auto rIter =
            IterateToFixedPoint(v, e, eps, /*maxIter=*/8, &iters, &status);
        ++iterDist[iters];
        // Note: rIter.inputRemap is the LAST iteration's remap, not the
        // composed original->final. The fuzz inputs are random topological
        // polygons where vertex merging produces an identity remap (no input verts
        // merge), so this is OK here.
        const bool convergedValid = CheckTopologicalValidity(
            rIter, e, rIter.inputRemap, rIter.numMergedVerts);
        if (convergedValid) {
          ++convergedPassValid;
          if (!pass1Valid) ++iter_repaired;
        } else {
          ++convergedPassInvalid;
          convergedInvalidList.emplace_back(kPow, n,
                                            static_cast<uint64_t>(seed),
                                            iters);
          if (pass1Valid) {
            ++iter_degraded;
            degradedList.emplace_back(kPow, n,
                                      static_cast<uint64_t>(seed), iters);
          } else {
            ++pass1_invalid_unfixed;
          }
        }

        // Area-drift check: compare pass 1's signed area against the
        // converged result's signed area. If pass 1 was meaningful
        // (non-empty, non-zero area) and the iteration changed it by
        // more than 1%, flag as a regression candidate. Special case:
        // pass-1-nonempty -> converged-empty is the thin-polygon
        // collapse mode (drift = 100%, area_drift_over_10pct picks it up).
        const double pass1Area =
            SignedAreaFromOutEdges(pass1.verts, pass1.edges);
        const double convArea =
            SignedAreaFromOutEdges(rIter.verts, rIter.edges);
        if (std::fabs(pass1Area) > 0) {
          ++area_pass1_nonzero;
          const double drift =
              std::fabs(convArea - pass1Area) / std::fabs(pass1Area);
          if (drift > area_drift_max) area_drift_max = drift;
          if (drift > 0.01) {
            ++area_drift_over_1pct;
            areaDriftList.emplace_back(kPow, n,
                                       static_cast<uint64_t>(seed), drift);
          }
          if (drift > 0.10) ++area_drift_over_10pct;
        }

        if (iters >= 2) {
          ++iterGE2;
          if (pass1Valid)
            ++iterGE2_pass1_valid;
          else
            ++iterGE2_pass1_invalid;
          // Hypothesis check: pass1 vs pass2 differ at eps/100 but match
          // at eps quantum (positions drift, topology doesn't)?
          std::vector<EdgeM> p2edges;
          for (const auto& oe : pass1.edges)
            p2edges.push_back({oe.v0, oe.v1, oe.mult});
          auto pass2 = RemoveOverlaps2D(pass1.verts, p2edges, eps);
          if (!pass1.edges.empty()) ++pass1_nonempty;
          if (!pass1.edges.empty() && pass2.edges.empty()) {
            ++pass2_collapsed;
            geomCollapseList.emplace_back(kPow, n,
                                          static_cast<uint64_t>(seed));
          }
          const auto fp1_fine = FingerprintAt(pass1, eps * 0.01);
          const auto fp2_fine = FingerprintAt(pass2, eps * 0.01);
          const auto fp1_topo = CoarseFingerprint(pass1, eps);
          const auto fp2_topo = CoarseFingerprint(pass2, eps);
          if (fp1_fine == fp2_fine) ++iterGE2_fine_match;
          if (fp1_topo == fp2_topo)
            ++iterGE2_topo_match;
          else {
            ++iterGE2_topo_diff;
            topoMismatch.emplace_back(kPow, n,
                                      static_cast<uint64_t>(seed));
          }
        }
      }
    }
  }

  std::cout << "  Total cases: " << total << "\n";
  std::cout << "  First-pass topology FAIL: " << firstPassFail << "\n";
  std::cout << "  Converged result topology valid:   " << convergedPassValid
            << "\n";
  std::cout << "  Converged result topology INVALID: " << convergedPassInvalid
            << "\n";
  std::cout << "    breakdown:\n";
  std::cout << "      iteration REPAIRED (invalid→valid):   " << iter_repaired
            << "\n";
  std::cout << "      iteration DEGRADED (valid→invalid):   " << iter_degraded
            << "\n";
  std::cout << "      pass-1 invalid + iteration didn't fix: "
            << pass1_invalid_unfixed << "\n";
  std::cout << "  Iter distribution: ";
  for (auto [k, c] : iterDist) std::cout << k << ":" << c << " ";
  std::cout << "\n";
  std::cout << "  iter≥2 cases: " << iterGE2 << "\n";
  std::cout << "    of which pass 1 was already topo-valid: "
            << iterGE2_pass1_valid << "\n";
  std::cout << "    of which pass 1 was topo-INVALID (iteration repaired): "
            << iterGE2_pass1_invalid << "\n";
  std::cout << "    where pass1 == pass2 at eps quantum (TOPO match): "
            << iterGE2_topo_match << "\n";
  std::cout << "    where pass1 == pass2 at eps/100 quantum (FINE match): "
            << iterGE2_fine_match << "\n";
  std::cout << "    where TOPO differs (real algorithm change): "
            << iterGE2_topo_diff << "\n";
  std::cout << "  Geometric collapse: pass1 nonempty=" << pass1_nonempty
            << " of which pass2=empty: " << pass2_collapsed
            << " (output orientation bug)\n";
  std::cout << "  Polygons regularization (zero-area lens drops): "
            << regDropped << " of " << total
            << " cases regularized away one or more zero-area loops "
               "(correct CGAL/Clipper2/SVG behavior; not a bug)\n";
  std::cout << "  Area drift (pass-1 vs converged, |delta|/|pass1|):\n";
  std::cout << "    pass-1 with measurable area: " << area_pass1_nonzero
            << "\n";
  std::cout << "    drift > 1%:  " << area_drift_over_1pct << "\n";
  std::cout << "    drift > 10%: " << area_drift_over_10pct
            << "  (typically thin-polygon collapse)\n";
  std::cout << "    max drift seen: " << (area_drift_max * 100.0) << "%\n";
  // For idempotence diagnosis: dump first 5 iter=2 cases (pass 1 valid).
  if (!topoMismatch.empty()) {
    std::cout << "\n  First 5 iter≥2 cases (for idempotence probe):\n";
    for (size_t i = 0; i < topoMismatch.size() && i < 5; ++i) {
      auto [kp, nn, sd] = topoMismatch[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd << "\n";
    }
  }
  // Geometric collapse cases: pass 1 nonempty, pass 2 empty. These are
  // the genuine bugs where iteration is changing the algorithm output.
  // Track them by re-running and checking. (Slow: requires re-running
  // pass 2 for each case. Done only for first 30 cases.)
  if (!geomCollapseList.empty()) {
    std::cout << "\n  Geometric COLLAPSE cases (pass 2 empty when pass 1 "
                 "nonempty):\n";
    for (size_t i = 0; i < geomCollapseList.size() && i < 20; ++i) {
      auto [kp, nn, sd] = geomCollapseList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd << "\n";
    }
  }
  if (!degradedList.empty()) {
    std::cout << "\n  DEGRADED cases (pass 1 valid, iteration broke it; first 20):\n";
    for (size_t i = 0; i < degradedList.size() && i < 20; ++i) {
      auto [kp, nn, sd, it] = degradedList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << " iters=" << it << "\n";
    }
  }
  if (!areaDriftList.empty()) {
    std::cout << "\n  Area-drift cases > 1% (first 10):\n";
    // Sort by drift magnitude so the worst cases surface first.
    std::sort(areaDriftList.begin(), areaDriftList.end(),
              [](const auto& a, const auto& b) {
                return std::get<3>(a) > std::get<3>(b);
              });
    for (size_t i = 0; i < areaDriftList.size() && i < 10; ++i) {
      auto [kp, nn, sd, dr] = areaDriftList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << "  drift=" << (dr * 100.0) << "%\n";
    }
  }
  if (!convergedInvalidList.empty()) {
    std::cout
        << "\n  CONVERGED-INVALID cases (iteration didn't fix; first 30):\n";
    for (size_t i = 0; i < convergedInvalidList.size() && i < 30; ++i) {
      auto [kp, nn, sd, it] = convergedInvalidList[i];
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << " iters=" << it << "\n";
    }
  }
}
}  // namespace boolean2
}  // namespace manifold

// =============================================================================
// OffsetFuzz: parameterized adversarial offset coverage.
//
// Mirror of DeepFuzz for the offset path. DeepFuzz exercises boolean's
// alpha-budget arithmetic against displacement-scaled random self-
// intersecting polygons; offsets need their own coverage because their
// failure modes (miter clamps, arc-chord subdivision, concave-side
// pinches, sub-eps thin-feature collapse) don't show up in boolean fuzz.
//
// Per-cell parameters: jointype x deltaFrac x size x kPow. Per-seed
// random simple polygon (star-shaped on a unit circle, displaced and
// scaled by 2^kPow), then run CrossSection::Offset and validate the
// output.
//
// Validators (any failure prints a diagnostic line and bumps a counter):
//   - Output coords all finite.
//   - Output area bounded by the actual input bbox inflated by |delta|.
//   - Output topology valid when re-run through RemoveOverlaps2D.
//   - For positive delta on non-degenerate input, output must be
//     non-empty (the inflated polygon must contain the original).
//
// Outputs a histogram of failure modes plus the first 10 failing
// (kPow, n, seed, jointype, deltaFrac) tuples for each failure category.
namespace manifold {
namespace boolean2 {

// Star-polygon generator: n verts at increasing angles around origin
// with random radii. Always simple (non-self-intersecting) because the
// boundary is monotonic in angle. Convex iff all radii equal; concave
// otherwise. Used as offset-fuzz input because offset on self-
// intersecting input is undefined.
inline manifold::SimplePolygon RandomSimpleStar(int n, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u01(0.0, 1.0);
  // Use equispaced angles + radial jitter to guarantee monotonic-in-angle
  // (and therefore simple, non-self-intersecting) polygons. An earlier
  // attempt that sampled random angles + sorted + enforced a minimum gap
  // could push cumulative angles past 2*pi when many gaps got bumped,
  // wrapping vertices past index 0 and producing a self-intersecting
  // polygon. Offset on a self-intersecting input is undefined, so all
  // those cases came back as area-runaway garbage.
  std::vector<double> angles(n);
  const double dTheta = 2.0 * M_PI / n;
  for (int i = 0; i < n; ++i) {
    // Jitter each angle within +/- 40% of half its slot so order is
    // preserved (no overtaking the next slot).
    const double jitter = 0.4 * dTheta * (2.0 * u01(rng) - 1.0);
    angles[i] = i * dTheta + jitter;
  }
  manifold::SimplePolygon ring(n);
  for (int i = 0; i < n; ++i) {
    // Radius in [0.5, 1.5]: gentle concavity.
    const double r = 0.5 + u01(rng);
    ring[i] = {r * std::cos(angles[i]), r * std::sin(angles[i])};
  }
  return ring;
}

inline const char* JtName(JoinType jt) {
  switch (jt) {
    case JoinType::Round: return "Round";
    case JoinType::Miter: return "Miter";
    case JoinType::Bevel: return "Bevel";
    case JoinType::Square: return "Square";
  }
  return "?";
}

inline JoinType ParseJt(const std::string& s) {
  if (s == "round" || s == "Round") return JoinType::Round;
  if (s == "miter" || s == "Miter") return JoinType::Miter;
  if (s == "bevel" || s == "Bevel") return JoinType::Bevel;
  if (s == "square" || s == "Square") return JoinType::Square;
  return JoinType::Round;
}

inline std::pair<vec2, vec2> BoundsOf(const Polygons& polys) {
  vec2 lo(std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity());
  vec2 hi(-std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity());
  for (const auto& loop : polys) {
    for (const auto& v : loop) {
      lo.x = std::min(lo.x, v.x);
      lo.y = std::min(lo.y, v.y);
      hi.x = std::max(hi.x, v.x);
      hi.y = std::max(hi.y, v.y);
    }
  }
  return {lo, hi};
}

#ifdef BOOLEAN2_WITH_CLIPPER2
inline Clipper2Lib::JoinType ToClipperJoin(JoinType jt) {
  switch (jt) {
    case JoinType::Round:
      return Clipper2Lib::JoinType::Round;
    case JoinType::Miter:
      return Clipper2Lib::JoinType::Miter;
    case JoinType::Bevel:
      return Clipper2Lib::JoinType::Bevel;
    case JoinType::Square:
      return Clipper2Lib::JoinType::Square;
  }
  return Clipper2Lib::JoinType::Round;
}

inline double ClipperAreaD(const Clipper2Lib::PathsD& paths) {
  double area = 0.0;
  for (const auto& path : paths) area += Clipper2Lib::Area(path);
  return area;
}
#endif

void OffsetCase(int kPow, int n, uint64_t seed64, JoinType jt,
                double deltaFrac) {
  const double scale = std::ldexp(1.0, kPow);
  auto ring = RandomSimpleStar(n, seed64);
  for (auto& v : ring) {
    v.x *= scale;
    v.y *= scale;
  }
  Polygons in = {ring};
  const auto [lo, hi] = BoundsOf(in);
  const double bboxExtent = std::max(hi.x - lo.x, hi.y - lo.y);
  const double delta = deltaFrac * bboxExtent;
  auto raw = offset_detail::OffsetContour(ring, delta, jt, 2.0, 0.0);
  Polygons rawPolys;
  if (raw.size() >= 3) rawPolys.push_back(raw);
  const double eps = rawPolys.empty() ? EpsilonFromScale(scale)
                                      : InferEps(rawPolys, {});
  auto rawNonZero = FillByRule(rawPolys, WindRule::NonZero, eps);
  auto rawPositive = FillByRule(rawPolys, WindRule::Add, eps);
  auto rawNegative = FillByRule(rawPolys, WindRule::Negative, eps);
  auto out = Offset(in, delta, jt, /*miterLimit=*/2.0, /*arcTol=*/0.0);
  std::cout << "=== OffsetCase ===\n";
  std::cout << "  kPow=" << kPow << " n=" << n << " seed=" << seed64
            << " jt=" << JtName(jt) << " deltaFrac=" << deltaFrac
            << " delta=" << delta << "\n";
  std::cout << "  input bbox: [" << lo.x << "," << lo.y << "] - [" << hi.x
            << "," << hi.y << "], extent=" << bboxExtent << "\n";
  std::cout << "  input area: " << TotalSignedArea(in) << "\n";
  std::cout << "  raw contour: verts=" << raw.size()
            << " area=" << SignedArea(raw) << "\n";
  std::cout << "  raw fill NonZero: paths=" << rawNonZero.size()
            << " area=" << TotalSignedArea(rawNonZero) << "\n";
  std::cout << "  raw fill Positive: paths=" << rawPositive.size()
            << " area=" << TotalSignedArea(rawPositive) << "\n";
  std::cout << "  raw fill Negative: paths=" << rawNegative.size()
            << " area=" << TotalSignedArea(rawNegative) << "\n";
  std::cout << "  boolean2: paths=" << out.size()
            << " verts=" << PolygonsToInput(out).first.size()
            << " area=" << TotalSignedArea(out) << "\n";
#ifdef BOOLEAN2_WITH_CLIPPER2
  auto clip = Clipper2Lib::InflatePaths(
      ToPathsD(in), delta, ToClipperJoin(jt), Clipper2Lib::EndType::Polygon,
      2.0, 0.0);
  std::cout << "  Clipper2D: paths=" << clip.size()
            << " area=" << ClipperAreaD(clip) << "\n";
#endif
}

void OffsetFuzz(int seedsPerCell) {
  using JT = JoinType;
  const std::vector<int> kPows = {0, 10, 30};
  const std::vector<int> sizes = {6, 20, 100};
  const std::vector<double> deltaFracs = {-0.5, -0.1, -0.01, 0.01, 0.1, 0.5};
  const std::vector<JT> jts = {JT::Round, JT::Miter, JT::Bevel, JT::Square};

  const int totalCells = static_cast<int>(kPows.size() * sizes.size() *
                                          deltaFracs.size() * jts.size());
  const int totalCases = totalCells * seedsPerCell;
  std::cout << "=== OffsetFuzz: " << seedsPerCell << " seeds x "
            << kPows.size() << " kPow x " << sizes.size() << " size x "
            << deltaFracs.size() << " deltaFrac x " << jts.size()
            << " jt = " << totalCases << " cases ===\n";

  int total = 0;
  int nonFinite = 0;
  int areaRunaway = 0;
  int topologyInvalid = 0;
  int emptyPositive = 0;
  int emptyNegative = 0;
#ifdef BOOLEAN2_WITH_CLIPPER2
  int clipperCollapsedOursKept = 0;
#endif
  using FailKey = std::tuple<int, int, uint64_t, int, double>;
  std::vector<FailKey> nonFiniteList, areaList, topoList, emptyPosList;
#ifdef BOOLEAN2_WITH_CLIPPER2
  std::vector<FailKey> collapseMismatchList;
#endif

  for (int kPow : kPows) {
    const double scale = std::ldexp(1.0, kPow);
    for (int n : sizes) {
      for (double deltaFrac : deltaFracs) {
        for (JT jt : jts) {
          for (int seed = 0; seed < seedsPerCell; ++seed) {
            ++total;
            const uint64_t seed64 =
                static_cast<uint64_t>(seed) + 1000ull * kPow + 10000ull * n;
            auto ring = RandomSimpleStar(n, seed64);
            for (auto& v : ring) {
              v.x *= scale;
              v.y *= scale;
            }
            manifold::Polygons in = {ring};
            const auto [lo, hi] = BoundsOf(in);
            const double bboxWidth = hi.x - lo.x;
            const double bboxHeight = hi.y - lo.y;
            const double bboxExtent = std::max(bboxWidth, bboxHeight);
            const double delta = deltaFrac * bboxExtent;
            auto out = Offset(in, delta, jt, /*miterLimit=*/2.0,
                              /*arcTol=*/0.0);
            bool anyNonFinite = false;
            double outArea = 0.0;
            for (const auto& loop : out) {
              outArea += SignedArea(loop);
              for (const auto& vv : loop) {
                if (!std::isfinite(vv.x) || !std::isfinite(vv.y))
                  anyNonFinite = true;
              }
            }
            if (anyNonFinite) {
              ++nonFinite;
              if (nonFiniteList.size() < 10)
                nonFiniteList.emplace_back(kPow, n, seed64, static_cast<int>(jt),
                                           deltaFrac);
              continue;
            }
            // Area-bound: post-offset bbox dimensions are bounded by the
            // actual input bbox inflated by |delta| in both directions. Allow
            // 4x slack for join geometry and FP jitter; this catches wild
            // runaway without rejecting valid large round/miter offsets.
            const double absDelta = std::fabs(delta);
            const double maxExpected =
                4.0 * (bboxWidth + 2.0 * absDelta) *
                (bboxHeight + 2.0 * absDelta);
            if (std::fabs(outArea) > maxExpected) {
              ++areaRunaway;
              if (areaList.size() < 10)
                areaList.emplace_back(kPow, n, seed64, static_cast<int>(jt),
                                      deltaFrac);
              continue;
            }
            if (!out.empty()) {
              const double eps = InferEps(out, {});
              auto [ov, oe] = PolygonsToInput(out);
              auto r2 = RemoveOverlaps2D(ov, oe, eps);
              if (!CheckTopologicalValidity(r2, oe, r2.inputRemap,
                                            r2.numMergedVerts)) {
                ++topologyInvalid;
                if (topoList.size() < 10)
                  topoList.emplace_back(kPow, n, seed64,
                                        static_cast<int>(jt), deltaFrac);
                continue;
              }
            }
#ifdef BOOLEAN2_WITH_CLIPPER2
            if (delta < 0) {
              auto clip = Clipper2Lib::InflatePaths(
                  ToPathsD(in), delta, ToClipperJoin(jt),
                  Clipper2Lib::EndType::Polygon, 2.0, 0.0);
              if (clip.empty() && !out.empty()) {
                ++clipperCollapsedOursKept;
                if (collapseMismatchList.size() < 10)
                  collapseMismatchList.emplace_back(
                      kPow, n, seed64, static_cast<int>(jt), deltaFrac);
              }
            }
#endif
            if (out.empty()) {
              if (delta > 0) {
                ++emptyPositive;
                if (emptyPosList.size() < 10)
                  emptyPosList.emplace_back(kPow, n, seed64,
                                            static_cast<int>(jt), deltaFrac);
              } else {
                ++emptyNegative;
              }
            }
          }
        }
      }
    }
  }
  std::cout << "  total cases:       " << total << "\n";
  std::cout << "  non-finite output: " << nonFinite << "\n";
  std::cout << "  area run-away:     " << areaRunaway << "\n";
  std::cout << "  topology invalid:  " << topologyInvalid << "\n";
  std::cout << "  empty (delta>0):   " << emptyPositive
            << "  (regression: positive delta should never collapse)\n";
  std::cout << "  empty (delta<0):   " << emptyNegative
            << "  (expected for inset-eats-shape cases)\n";
#ifdef BOOLEAN2_WITH_CLIPPER2
  std::cout << "  Clipper empty, boolean2 non-empty (delta<0): "
            << clipperCollapsedOursKept << "\n";
#endif
  auto dumpList = [&](const char* label, const std::vector<FailKey>& list) {
    if (list.empty()) return;
    std::cout << "\n  " << label << " (first " << list.size() << "):\n";
    for (const auto& [kp, nn, sd, jti, df] : list) {
      std::cout << "    kPow=" << kp << " n=" << nn << " seed=" << sd
                << " jt=" << JtName(static_cast<JT>(jti))
                << " deltaFrac=" << df << "\n";
    }
  };
  dumpList("non-finite cases", nonFiniteList);
  dumpList("area-runaway cases", areaList);
  dumpList("topology-invalid cases", topoList);
  dumpList("unexpected-empty (positive delta) cases", emptyPosList);
#ifdef BOOLEAN2_WITH_CLIPPER2
  dumpList("Clipper-collapsed / boolean2-kept cases", collapseMismatchList);
#endif
}
}  // namespace boolean2
}  // namespace manifold

namespace {

void PrintUsage(const char* argv0) {
  std::cout
      << "Usage: " << argv0 << " [mode] [args]\n\n"
      << "Core modes:\n"
      << "  <none>                    run the built-in smoke/fuzz battery\n"
      << "  diagnose [seed [kPow [n]]]\n"
      << "  deepfuzz [seedsPerCell]\n"
      << "  offsetfuzz [seedsPerCell]\n"
      << "  offsetcase [kPow n seed jointype deltaFrac]\n"
      << "  corpus [path]\n"
      << "  clipper2corpus [path]\n"
      << "  mfogelcorpus [dir]\n"
      << "  jtscorpus [path]\n"
      << "  cadcorpus\n"
      << "  time [--repeat N] [mode [args...]]\n"
#ifdef BOOLEAN2_WITH_CLIPPER2
      << "\nClipper2 comparison/diagnostic modes:\n"
      << "  vsclipper2 [all|clipper2|cad|offsets|jts|mfogel|manifoldcorpus]\n"
      << "  jtsdrops [path]\n"
      << "  jtsdiag <case_n>\n"
      << "  manifolddiag <case_n>\n"
      << "  manifoldpick <case_n> <ring_idx> ...\n"
      << "  manifoldsubrun <case_n> <num_rings> [skip]\n"
#endif
      << "\nDebug modes:\n"
      << "  pentagon [kPow]\n"
      << "  idempotence <seed> [kPow [n]]\n";
}

}  // namespace

int main(int argc, char** argv) {
  using namespace manifold::boolean2;
  if (argc > 1) {
    const std::string mode = argv[1];
    if (mode == "-h" || mode == "--help" || mode == "help") {
      PrintUsage(argv[0]);
      return 0;
    }
  }
  if (argc > 1 && std::string(argv[1]) == "diagnose") {
    uint64_t seed = (argc > 2) ? std::stoull(argv[2]) : 0;
    int kPow = (argc > 3) ? std::atoi(argv[3]) : 30;
    int n = (argc > 4) ? std::atoi(argv[4]) : 50;
    Diagnose(seed, kPow, n);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "deepfuzz") {
    int seedsPerCell = (argc > 2) ? std::atoi(argv[2]) : 200;
    DeepFuzz(seedsPerCell);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "offsetfuzz") {
    int seedsPerCell = (argc > 2) ? std::atoi(argv[2]) : 50;
    OffsetFuzz(seedsPerCell);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "offsetcase") {
    int kPow = (argc > 2) ? std::atoi(argv[2]) : 0;
    int n = (argc > 3) ? std::atoi(argv[3]) : 100;
    uint64_t seed = (argc > 4) ? std::stoull(argv[4]) : 1000000ull;
    JoinType jt = (argc > 5) ? ParseJt(argv[5]) : JoinType::Round;
    double deltaFrac = (argc > 6) ? std::atof(argv[6]) : -0.5;
    OffsetCase(kPow, n, seed, jt, deltaFrac);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "corpus") {
    const std::string path = (argc > 2)
                                  ? argv[2]
                                  : "test/polygons/polygon_corpus.txt";
    RunCorpus(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "clipper2corpus") {
    std::string path =
        "build/_deps/clipper2-src/Tests/Polygons.txt";
    if (argc > 2) path = argv[2];
    RunClipper2Corpus(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "mfogelcorpus") {
    std::string dir =
        "build/_deps/mfogel-polygon-clipping-src/test/end-to-end";
    if (argc > 2) dir = argv[2];
    RunMfogelCorpus(dir);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtscorpus") {
    std::string path = "test/polygons/jts_overlay_corpus.txt";
    if (argc > 2) path = argv[2];
    RunJtsCorpus(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "cadcorpus") {
    RunCadCorpus();
    return 0;
  }
#ifdef BOOLEAN2_WITH_CLIPPER2
  if (argc > 1 && std::string(argv[1]) == "vsclipper2") {
    std::string which = (argc > 2) ? argv[2] : "all";
    RunVsClipper2(which);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtsdrops") {
    std::string path = (argc > 2) ? argv[2] : "test/polygons/jts_overlay_corpus.txt";
    RunJtsDropDiag(path);
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "jtsdiag") {
    if (argc < 3) {
      std::cerr << "Usage: jtsdiag <case_n>\n";
      return 1;
    }
    int targetCase = std::atoi(argv[2]);
    auto cases = LoadJtsCorpus("test/polygons/jts_overlay_corpus.txt");
    for (const auto& c : cases) {
      if (c.n != targetCase) continue;
      std::cout << "=== JTS case " << c.n << " ===\n";
      std::cout << "  op=" << c.op << " src=" << c.source << "\n";
      std::cout << "  A.rings=" << c.a.size() << " B.rings=" << c.b.size() << "\n";
      double xMin = 1e300, yMin = 1e300, xMax = -1e300, yMax = -1e300;
      int totalVerts = 0;
      auto bound = [&](const manifold::Polygons& polys) {
        for (const auto& loop : polys) {
          totalVerts += loop.size();
          for (const auto& v : loop) {
            xMin = std::min(xMin, v.x); yMin = std::min(yMin, v.y);
            xMax = std::max(xMax, v.x); yMax = std::max(yMax, v.y);
          }
        }
      };
      bound(c.a); bound(c.b);
      const double eps = InferEps(c.a, c.b);
      std::cout << "  bbox: x=[" << xMin << ".." << xMax << "] y=[" << yMin << ".." << yMax << "]\n";
      std::cout << "  total verts: " << totalVerts << " eps: " << eps << "\n";
      // Run our pipeline with debug=true.
      const manifold::Polygons& A = c.a;
      const manifold::Polygons& B = (c.op == "overlayAreaTest" || c.op == "overlayareatest") ? c.b : manifold::Polygons{};
      std::vector<vec2> verts;
      std::vector<EdgeM> edges;
      auto append = [&](const manifold::Polygons& polys) {
        for (const auto& loop : polys) {
          if (loop.size() < 3) continue;
          const int base = static_cast<int>(verts.size());
          const int n = static_cast<int>(loop.size());
          for (const auto& v : loop) verts.push_back(v);
          for (int i = 0; i < n; ++i)
            edges.push_back({base + i, base + ((i + 1) % n), 1});
        }
      };
      append(A); append(B);
      std::cout << "  Pipeline input: " << verts.size() << " verts, " << edges.size() << " edges\n";
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/true,
                                [](int w) { return w != 0; });
      auto out = OutEdgesToPolygons(r.verts, r.edges);
      std::cout << "  Pipeline output: " << r.verts.size() << " verts, " << r.edges.size() << " edges, " << out.size() << " polygons\n";
      // Clipper2 reference for comparison.
      auto a2 = ToPaths64(A);
      auto b2 = ToPaths64(B);
      Clipper2Lib::Clipper64 clip;
      clip.AddSubject(a2);
      if (!B.empty()) clip.AddClip(b2);
      Clipper2Lib::Paths64 sol;
      clip.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero, sol);
      double clipArea = 0; for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      std::cout << "  Clipper2 output: " << sol.size() << " paths, area=" << clipArea << "\n";
      return 0;
    }
    std::cerr << "Case " << targetCase << " not found\n";
    return 1;
  }
  // manifoldpick <case_n> <ring_idx_1> <ring_idx_2> ... - runs a custom
  // subset by listing exact ring indices (0-based).
  if (argc > 1 && std::string(argv[1]) == "manifoldpick") {
    if (argc < 4) {
      std::cerr << "Usage: manifoldpick <case_n> <ring_idx> ...\n";
      return 1;
    }
    int targetCase = std::atoi(argv[2]);
    std::vector<int> picks;
    for (int i = 3; i < argc; ++i) picks.push_back(std::atoi(argv[i]));
    auto cases = LoadJtsCorpus("test/polygons/manifold_corpus.txt");
    for (const auto& c : cases) {
      if (c.n != targetCase) continue;
      manifold::Polygons sub;
      for (int i : picks) {
        if (i < 0 || i >= (int)c.a.size()) {
          std::cerr << "ring " << i << " out of range\n";
          return 1;
        }
        sub.push_back(c.a[i]);
      }
      const double eps = InferEps(sub, {});
      std::vector<vec2> verts;
      std::vector<EdgeM> edges;
      for (const auto& loop : sub) {
        if (loop.size() < 3) continue;
        const int base = static_cast<int>(verts.size());
        const int n = static_cast<int>(loop.size());
        for (const auto& v : loop) verts.push_back(v);
        for (int i = 0; i < n; ++i)
          edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
      const bool debug = std::getenv("BOOLEAN2_DEBUG") != nullptr;
      auto r = RemoveOverlaps2D(verts, edges, eps, debug,
                                [](int w) { return w > 0; });
      auto out = OutEdgesToPolygons(r.verts, r.edges);
      double oursArea = 0;
      for (const auto& l : out) {
        if (l.empty()) continue;
        const auto& r0 = l.front();
        for (size_t i = 0; i < l.size(); ++i) {
          const auto& a = l[i];
          const auto& b = l[(i + 1) % l.size()];
          oursArea += 0.5 *
              ((a.x - r0.x) * (b.y - r0.y) - (b.x - r0.x) * (a.y - r0.y));
        }
      }
      auto a2 = ToPaths64(sub);
      Clipper2Lib::Clipper64 clip;
      clip.AddSubject(a2);
      Clipper2Lib::Paths64 sol;
      clip.Execute(Clipper2Lib::ClipType::Union,
                   Clipper2Lib::FillRule::Positive, sol);
      double clipArea = 0;
      for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      clipArea /= (kPaths64Scale * kPaths64Scale);
      std::cout << "picks=" << picks.size() << " ours=" << out.size()
                << "/" << std::setprecision(8) << oursArea
                << " theirs=" << sol.size() << "/" << clipArea
                << " diff=" << (oursArea - clipArea)
                << (std::fabs(oursArea - clipArea) > 1e-6 ? " DIFFER" : "")
                << "\n";
      if (std::getenv("BOOLEAN2_DUMP_OUT")) {
        auto ringArea = [](const manifold::SimplePolygon& l) {
          double a = 0;
          if (l.empty()) return 0.0;
          const auto& r0 = l.front();
          for (size_t i = 0; i < l.size(); ++i) {
            const auto& p = l[i];
            const auto& q = l[(i + 1) % l.size()];
            a += 0.5 *
                ((p.x - r0.x) * (q.y - r0.y) - (q.x - r0.x) * (p.y - r0.y));
          }
          return a;
        };
        std::cout << "OURS (" << out.size() << "):\n";
        for (size_t i = 0; i < out.size(); ++i) {
          std::cout << "  poly " << i << " verts=" << out[i].size()
                    << " area=" << ringArea(out[i]) << "\n";
        }
        std::cout << "THEIRS (" << sol.size() << "):\n";
        for (size_t i = 0; i < sol.size(); ++i) {
          double a = Clipper2Lib::Area(sol[i]) / (kPaths64Scale * kPaths64Scale);
          std::cout << "  poly " << i << " verts=" << sol[i].size()
                    << " area=" << a << "\n";
        }
      }
      return 0;
    }
    std::cerr << "Case " << targetCase << " not found\n";
    return 1;
  }
  // manifoldsubrun <case_n> <num_rings> [skip]
  // Runs boolean2 + Clipper2 on the first `num_rings` of case_n. If
  // `skip` is provided, skips that many rings at the start (so you can
  // test [skip..skip+num_rings)). Reports area diff. Used for bisecting
  // which input subset triggers the discrepancy.
  if (argc > 1 && std::string(argv[1]) == "manifoldsubrun") {
    if (argc < 4) {
      std::cerr << "Usage: manifoldsubrun <case_n> <num_rings> [skip]\n";
      return 1;
    }
    int targetCase = std::atoi(argv[2]);
    int numRings = std::atoi(argv[3]);
    int skip = argc > 4 ? std::atoi(argv[4]) : 0;
    auto cases = LoadJtsCorpus("test/polygons/manifold_corpus.txt");
    for (const auto& c : cases) {
      if (c.n != targetCase) continue;
      // Slice the input rings to [skip..skip+numRings).
      manifold::Polygons sub;
      int end = std::min((int)c.a.size(), skip + numRings);
      for (int i = skip; i < end; ++i) sub.push_back(c.a[i]);
      const double eps = InferEps(sub, {});
      std::vector<vec2> verts;
      std::vector<EdgeM> edges;
      for (const auto& loop : sub) {
        if (loop.size() < 3) continue;
        const int base = static_cast<int>(verts.size());
        const int n = static_cast<int>(loop.size());
        for (const auto& v : loop) verts.push_back(v);
        for (int i = 0; i < n; ++i)
          edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/false,
                                [](int w) { return w > 0; });
      auto out = OutEdgesToPolygons(r.verts, r.edges);
      double oursArea = 0;
      for (const auto& l : out) {
        if (l.empty()) continue;
        const auto& r0 = l.front();
        for (size_t i = 0; i < l.size(); ++i) {
          const auto& a = l[i];
          const auto& b = l[(i + 1) % l.size()];
          oursArea += 0.5 *
              ((a.x - r0.x) * (b.y - r0.y) - (b.x - r0.x) * (a.y - r0.y));
        }
      }
      auto a2 = ToPaths64(sub);
      Clipper2Lib::Clipper64 clip;
      clip.AddSubject(a2);
      Clipper2Lib::Paths64 sol;
      clip.Execute(Clipper2Lib::ClipType::Union,
                   Clipper2Lib::FillRule::Positive, sol);
      double clipArea = 0;
      for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      clipArea /= (kPaths64Scale * kPaths64Scale);
      double diff = oursArea - clipArea;
      std::cout << "rings=[" << skip << "," << end << ") n=" << (end - skip)
                << " ours=" << out.size() << "/" << std::setprecision(8)
                << oursArea << " theirs=" << sol.size() << "/" << clipArea
                << " diff=" << diff
                << (std::fabs(diff) > 1e-6 ? " DIFFER" : "") << "\n";
      return 0;
    }
    std::cerr << "Case " << targetCase << " not found\n";
    return 1;
  }
  // Manifold-corpus diagnostic: like jtsdiag but uses WindRule::Add
  // (Positive fillrule, what CrossSection's ctors call) instead of
  // NonZero. Targets the manifold corpus.
  if (argc > 1 && std::string(argv[1]) == "manifolddiag") {
    if (argc < 3) {
      std::cerr << "Usage: manifolddiag <case_n>\n";
      return 1;
    }
    int targetCase = std::atoi(argv[2]);
    auto cases = LoadJtsCorpus("test/polygons/manifold_corpus.txt");
    for (const auto& c : cases) {
      if (c.n != targetCase) continue;
      std::cout << "=== Manifold corpus case " << c.n << " ===\n";
      std::cout << "  op=" << c.op << " src=" << c.source << "\n";
      std::cout << "  A.rings=" << c.a.size() << "\n";
      for (size_t i = 0; i < c.a.size(); ++i) {
        double area = 0;
        const auto& r = c.a[i].front();
        for (size_t k = 0; k < c.a[i].size(); ++k) {
          const auto& a = c.a[i][k];
          const auto& b = c.a[i][(k + 1) % c.a[i].size()];
          area += 0.5 * ((a.x - r.x) * (b.y - r.y) - (b.x - r.x) * (a.y - r.y));
        }
        std::cout << "    A[" << i << "] verts=" << c.a[i].size()
                  << " signedArea=" << area << "\n";
      }
      double eps = InferEps(c.a, {});
      // Optional eps override via BOOLEAN2_EPS env var (absolute) or
      // BOOLEAN2_EPS_MULT (multiplier on inferred eps). Useful for
      // testing sensitivity on the manifold corpus.
      if (const char* m = std::getenv("BOOLEAN2_EPS_MULT")) {
        eps *= std::atof(m);
      }
      if (const char* e = std::getenv("BOOLEAN2_EPS")) {
        eps = std::atof(e);
      }
      std::cout << "  eps=" << eps << "\n";
      std::vector<vec2> verts;
      std::vector<EdgeM> edges;
      for (const auto& loop : c.a) {
        if (loop.size() < 3) continue;
        const int base = static_cast<int>(verts.size());
        const int n = static_cast<int>(loop.size());
        for (const auto& v : loop) verts.push_back(v);
        for (int i = 0; i < n; ++i)
          edges.push_back({base + i, base + ((i + 1) % n), 1});
      }
      std::cout << "  Pipeline input: " << verts.size() << " verts, "
                << edges.size() << " edges\n";
      const char* ruleEnv = std::getenv("BOOLEAN2_RULE");
      std::string rule = ruleEnv ? std::string(ruleEnv) : "positive";
      auto rPred = [&](int w) {
        if (rule == "nonzero") return w != 0;
        if (rule == "evenodd") return (w & 1) != 0;
        if (rule == "negative") return w < 0;
        return w > 0;  // positive (default)
      };
      std::cout << "  rule=" << rule << "\n";
      auto r = RemoveOverlaps2D(verts, edges, eps, /*debug=*/true, rPred);
      auto out = OutEdgesToPolygons(r.verts, r.edges);
      std::cout << "  Pipeline output: " << r.verts.size() << " verts, "
                << r.edges.size() << " edges, " << out.size() << " polygons\n";
      // Also try IterateToFixedPoint with more passes.
      auto iter = IterateToFixedPoint(verts, edges, eps, /*maxIter=*/5);
      auto outIter = OutEdgesToPolygons(iter.verts, iter.edges);
      double iterArea = 0;
      for (const auto& l : outIter) {
        const auto& r0 = l.front();
        for (size_t i = 0; i < l.size(); ++i) {
          const auto& a = l[i];
          const auto& b = l[(i + 1) % l.size()];
          iterArea += 0.5 *
              ((a.x - r0.x) * (b.y - r0.y) - (b.x - r0.x) * (a.y - r0.y));
        }
      }
      std::cout << "  IterateToFixedPoint(maxIter=5): " << outIter.size()
                << " polygons, area=" << iterArea << "\n";
      auto a2 = ToPaths64(c.a);
      Clipper2Lib::Clipper64 clip;
      clip.AddSubject(a2);
      Clipper2Lib::Paths64 sol;
      Clipper2Lib::FillRule cFr = Clipper2Lib::FillRule::Positive;
      if (rule == "nonzero") cFr = Clipper2Lib::FillRule::NonZero;
      if (rule == "evenodd") cFr = Clipper2Lib::FillRule::EvenOdd;
      if (rule == "negative") cFr = Clipper2Lib::FillRule::Negative;
      clip.Execute(Clipper2Lib::ClipType::Union, cFr, sol);
      double clipArea = 0;
      for (const auto& p : sol) clipArea += Clipper2Lib::Area(p);
      clipArea /= (kPaths64Scale * kPaths64Scale);
      std::cout << "  Clipper2 (" << rule << "): " << sol.size()
                << " paths, area=" << clipArea << "\n";
      // Per-polygon area histogram + sorted diff (boolean2 vs Clipper2).
      auto ringArea = [](const manifold::SimplePolygon& l) {
        double a = 0;
        if (l.empty()) return 0.0;
        const auto& r0 = l.front();
        for (size_t i = 0; i < l.size(); ++i) {
          const auto& p = l[i];
          const auto& q = l[(i + 1) % l.size()];
          a += 0.5 *
              ((p.x - r0.x) * (q.y - r0.y) - (q.x - r0.x) * (p.y - r0.y));
        }
        return a;
      };
      std::vector<double> oursAreas, theirsAreas;
      for (const auto& l : out) oursAreas.push_back(ringArea(l));
      for (const auto& p : sol) {
        double a = Clipper2Lib::Area(p) / (kPaths64Scale * kPaths64Scale);
        theirsAreas.push_back(a);
      }
      std::sort(oursAreas.begin(), oursAreas.end());
      std::sort(theirsAreas.begin(), theirsAreas.end());
      std::cout << "  --- per-polygon areas (sorted) ---\n";
      std::cout << "  ours    counts: ";
      auto histBin = [](double a) {
        if (a == 0) return std::string("zero");
        double ab = std::fabs(a);
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.0e", ab);
        return std::string(buf);
      };
      auto histogram = [&](const std::vector<double>& v) {
        std::map<std::string, std::pair<int, int>> b;  // (cw, ccw)
        for (double a : v) {
          auto& p = b[histBin(a)];
          if (a < 0) ++p.first;
          else ++p.second;
        }
        for (auto& [k, p] : b) {
          std::cout << "    bin=" << k << " cw=" << p.first
                    << " ccw=" << p.second << "\n";
        }
      };
      std::cout << "  --- ours histogram by |area| ---\n";
      histogram(oursAreas);
      std::cout << "  --- theirs histogram by |area| ---\n";
      histogram(theirsAreas);
      // Print outer ring area to ~17 digits to see if it diverges.
      double oursOuter = 0, theirsOuter = 0;
      for (double a : oursAreas) if (a > 0) oursOuter += a;
      for (double a : theirsAreas) if (a > 0) theirsOuter += a;
      double oursHoles = 0, theirsHoles = 0;
      for (double a : oursAreas) if (a < 0) oursHoles += a;
      for (double a : theirsAreas) if (a < 0) theirsHoles += a;
      std::cout << "  positives sum  ours=" << std::setprecision(17)
                << oursOuter << "  theirs=" << theirsOuter << "\n";
      std::cout << "  negatives sum  ours=" << oursHoles
                << "  theirs=" << theirsHoles << "\n";
      std::cout << "  net area       ours=" << (oursOuter + oursHoles)
                << "  theirs=" << (theirsOuter + theirsHoles) << "\n";
      return 0;
    }
    std::cerr << "Case " << targetCase << " not found\n";
    return 1;
  }
#endif
  if (argc > 1 && std::string(argv[1]) == "time") {
    // Per-phase wall-clock timing on a representative workload. Default
    // workload is DeepFuzz with seedsPerCell=200 (5,600 cases × ~3
    // RemoveOverlaps2D calls each ≈ 17k pipeline runs); selectable by a
    // second positional arg matching another mode name. Examples:
    //
    //   ./boolean2_driver time                 # deepfuzz 200
    //   ./boolean2_driver time deepfuzz 1000   # 28k cases
    //   ./boolean2_driver time offsetfuzz 50   # offset-path coverage
    //   ./boolean2_driver time --repeat 5 clipper2corpus
    //   ./boolean2_driver time jtscorpus       # JTS overlay corpus
    //   ./boolean2_driver time clipper2corpus
    //   ./boolean2_driver time mfogelcorpus
    //   ./boolean2_driver time cadcorpus
    //   ./boolean2_driver time corpus          # polygon_corpus.txt
    //
    // Output: total wall-clock + per-phase ns + percentage breakdown.
    SetTimingEnabled(true);
    GlobalPhases().Reset();
    int arg = 2;
    int repeats = 1;
    if (argc > arg && std::string(argv[arg]) == "--repeat") {
      if (argc <= arg + 1) {
        std::cerr << "missing repeat count\n";
        return 2;
      }
      repeats = std::max(1, std::atoi(argv[arg + 1]));
      arg += 2;
    }
    const auto wallStart = std::chrono::steady_clock::now();
    std::string sub = (argc > arg) ? argv[arg] : "deepfuzz";
    // Silence per-mode stdout so timing summary is uncluttered.
    std::ostringstream sink;
    std::streambuf* oldCout = std::cout.rdbuf(sink.rdbuf());
    for (int r = 0; r < repeats; ++r) {
      if (sub == "deepfuzz") {
        int seeds = (argc > arg + 1) ? std::atoi(argv[arg + 1]) : 200;
        DeepFuzz(seeds);
      } else if (sub == "offsetfuzz") {
        int seeds = (argc > arg + 1) ? std::atoi(argv[arg + 1]) : 50;
        OffsetFuzz(seeds);
      } else if (sub == "corpus") {
        const std::string p = (argc > arg + 1) ? argv[arg + 1] : "test/polygons/polygon_corpus.txt";
        RunCorpus(p);
      } else if (sub == "clipper2corpus") {
        const std::string p = (argc > arg + 1) ? argv[arg + 1]
          : std::string("build/_deps/clipper2-src/Tests/Polygons.txt");
        RunClipper2Corpus(p);
      } else if (sub == "mfogelcorpus") {
        const std::string p = (argc > arg + 1) ? argv[arg + 1]
          : std::string("build/_deps/mfogel-polygon-clipping-src/test/end-to-end");
        RunMfogelCorpus(p);
      } else if (sub == "jtscorpus") {
        const std::string p = (argc > arg + 1) ? argv[arg + 1] : "test/polygons/jts_overlay_corpus.txt";
        RunJtsCorpus(p);
      } else if (sub == "cadcorpus") {
        RunCadCorpus();
      } else {
        std::cout.rdbuf(oldCout);
        std::cerr << "unknown subcommand: " << sub << "\n";
        return 2;
      }
    }
    std::cout.rdbuf(oldCout);
    const auto wallEnd = std::chrono::steady_clock::now();
    const int64_t wallNs = std::chrono::duration_cast<std::chrono::nanoseconds>(
        wallEnd - wallStart).count();
    const auto& P = GlobalPhases();
    const int64_t total = P.totalNs.load();
    const int64_t cases = P.cases.load();
    auto pct = [&](int64_t ns) {
      return total > 0 ? (ns * 100.0 / total) : 0.0;
    };
    auto fmt = [](int64_t ns) {
      double ms = ns * 1e-6;
      std::ostringstream o;
      o.setf(std::ios::fixed);
      o.precision(3);
      o << ms << " ms";
      return o.str();
    };
    std::cout << "=== timing: " << sub << " ===\n";
    if (repeats > 1) std::cout << "  Repeats:                 " << repeats << "\n";
    std::cout << "  Wall (driver+pipeline): " << fmt(wallNs) << "\n";
    std::cout << "  Pipeline only (sum):    " << fmt(total) << "\n";
    std::cout << "  RemoveOverlaps2D calls: " << cases << "\n";
    if (cases > 0) {
      std::cout << "  Avg per call:           " << fmt(total / cases) << "\n";
    }
    std::cout << "\n  Per-phase breakdown (% of pipeline total):\n";
    auto row = [&](const char* name, int64_t ns) {
      std::cout << "    " << std::left << std::setw(28) << name
                << std::right << std::setw(12) << fmt(ns)
                << "  " << std::setw(6);
      std::cout.setf(std::ios::fixed);
      std::cout.precision(2);
      std::cout << pct(ns) << "%\n";
      std::cout.unsetf(std::ios::fixed);
    };
    row("vertex merge",                P.mergeNs.load());
    row("edge collapse",               P.remapNs.load());
    row("near-vertex indexing",        P.buildListsNs.load());
    row("intersection insertion",      P.findIxNs.load());
    row("structural re-merge",         P.restructNs.load());
    row("canonicalization",            P.canonNs.load());
    row("winding filter",              P.filterDcelNs.load());
    std::cout << "\n  Sub-phase breakdown (% of pipeline total):\n";
    row("  merge BVH build",            P.mergeBvhBuildNs.load());
    row("  merge BVH self-collide",     P.mergeCollideNs.load());
    row("  shared edge BVH build",      P.bvhBuildNs.load());
    row("  broad candidate work",       P.broadPairWorkNs.load());
    row("  edge-vertex lists",          P.edgeVertListsNs.load());
    row("  intersection broad",         P.intersectionBroadNs.load());
    row("  intersection narrow",        P.intersectionNarrowNs.load());
    row("  intersection propagation",   P.intersectionPropagationNs.load());
    const int64_t evCand = P.edgeVertCandidates.load();
    if (evCand > 0) {
      std::cout << "\n  Edge-vertex candidate breakdown:\n";
      auto crow = [&](const char* name, int64_t count) {
        std::cout << "    " << std::left << std::setw(28) << name
                  << std::right << std::setw(12) << count
                  << "  " << std::setw(6);
        std::cout.setf(std::ios::fixed);
        std::cout.precision(2);
        std::cout << (count * 100.0 / evCand) << "%\n";
        std::cout.unsetf(std::ios::fixed);
      };
      crow("candidates", evCand);
      crow("endpoint rejects", P.edgeVertEndpointRejects.load());
      crow("degenerate edge rejects", P.edgeVertDegenerateRejects.load());
      crow("t-range rejects", P.edgeVertTRangeRejects.load());
      crow("distance rejects", P.edgeVertDistanceRejects.load());
      crow("apex rejects", P.edgeVertApexRejects.load());
      crow("hits", P.edgeVertHits.load());
      const int64_t evCalls = P.edgeVertCalls.load();
      const int64_t evEdges = P.edgeVertTotalEdges.load();
      const int64_t evVerts = P.edgeVertTotalVerts.load();
      if (evCalls > 0) {
        std::cout << "\n  Edge-vertex shape breakdown:\n";
        std::cout << "    calls=" << evCalls
                  << " brute=" << P.edgeVertBruteCalls.load()
                  << " bvh=" << P.edgeVertBvhCalls.load()
                  << " avgE=" << (evEdges * 1.0 / evCalls)
                  << " avgV=" << (evVerts * 1.0 / evCalls)
                  << " flatHits=" << P.edgeVertHitsFlat.load() << "\n";
        std::cout << "    edge buckets: <64=" << P.edgeVertBucketLt64.load()
                  << " <256=" << P.edgeVertBucketLt256.load()
                  << " <1024=" << P.edgeVertBucketLt1024.load()
                  << " >=1024=" << P.edgeVertBucketGe1024.load() << "\n";
      }
    }
    const int64_t propCalls = P.propagationCalls.load();
    if (propCalls > 0) {
      std::cout << "\n  Intersection propagation shape:\n";
      std::cout << "    calls=" << propCalls
                << " skippedNoNearDup="
                << P.propagationSkippedNoNearDup.load() << "\n";
    }
    return 0;
  }
  SetTimingEnabled(false);
  if (argc > 1 && std::string(argv[1]) == "pentagon") {
    // Sanity: feed a clean pentagon (displaced if requested), see what DCEL
    // does. With kPow arg, displaces to that scale to mimic re-feed conditions.
    int kPow = (argc > 2) ? std::atoi(argv[2]) : 0;
    const double offset = (kPow > 0) ? std::ldexp(1.5, kPow) : 0.0;
    const double scale = (kPow > 0) ? offset : 1.0;
    std::vector<vec2> v = {
        {std::cos(0.0), std::sin(0.0)},
        {std::cos(2 * M_PI / 5), std::sin(2 * M_PI / 5)},
        {std::cos(4 * M_PI / 5), std::sin(4 * M_PI / 5)},
        {std::cos(6 * M_PI / 5), std::sin(6 * M_PI / 5)},
        {std::cos(8 * M_PI / 5), std::sin(8 * M_PI / 5)},
    };
    if (kPow > 0)
      for (auto& p : v) {
        p.x += offset;
        p.y += offset;
      }
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1},
                            {3, 4, 1}, {4, 0, 1}};
    auto r = RemoveOverlaps2D(v, e, EpsilonFromScale(scale));
    std::cerr << "Pentagon (kPow=" << kPow << "): " << r.verts.size()
              << " verts, " << r.edges.size() << " edges\n";
    for (const auto& oe : r.edges) {
      std::cerr << "  " << oe.v0 << "→" << oe.v1 << " mult=" << oe.mult
                << "\n";
    }
    return 0;
  }
  if (argc > 1 && std::string(argv[1]) == "idempotence") {
    // ./boolean2_driver idempotence <seed> [kPow] [n]
    // Compares pass 1 vs pass 2 sub-edge sets to identify what's drifting.
    uint64_t seed = (argc > 2) ? std::stoull(argv[2]) : 0;
    int kPow = (argc > 3) ? std::atoi(argv[3]) : 10;
    int n = (argc > 4) ? std::atoi(argv[4]) : 8;
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    auto [vIn, eIn] = RandomTopologicalPolygon(n, seed + 1000ull * kPow);
    Displace(&vIn, offset);
    std::cerr << "=== idempotence: kPow=" << kPow << " n=" << n
              << " seed=" << seed << " eps=" << eps << " ===\n";
    std::cerr << "--- pass 1 ---\n";
    auto pass1 = RemoveOverlaps2D(vIn, eIn, eps, /*debug=*/true);
    std::vector<EdgeM> p2in;
    for (const auto& oe : pass1.edges)
      p2in.push_back({oe.v0, oe.v1, oe.mult});
    std::cerr << "--- pass 2 (input: " << pass1.verts.size() << " verts, "
              << p2in.size() << " edges) ---\n";
    // Manually run pass 2's pipeline to see where edges disappear.
    auto p2_mrg = MergeVerts(pass1.verts, eps);
    std::cerr << "  vertex merge: " << pass1.verts.size() << "→"
              << p2_mrg.verts.size() << "\n";
    auto p2_e = RemapAndCollapse(p2in, p2_mrg.remap);
    std::cerr << "  edge collapse: " << p2in.size() << "→" << p2_e.size()
              << " edges\n";
    std::vector<Box2> p2_eBoxes(p2_e.size());
    for (size_t i = 0; i < p2_e.size(); ++i)
      p2_eBoxes[i] = BoxOf2DEdge(p2_mrg.verts[p2_e[i].v0],
                                 p2_mrg.verts[p2_e[i].v1], eps);
    BVH p2_bvh = BVHBuildFromBoxes(p2_eBoxes);
    auto p2_l = BuildEdgeVertLists(p2_e, p2_mrg.verts, eps, p2_eBoxes, p2_bvh);
    int p2_totalList = 0;
    for (auto& l : p2_l) p2_totalList += l.size();
    std::cerr << "  near-vertex lists: " << p2_totalList
              << " total entries\n";
    std::vector<std::vector<int>> p2_ve;
    auto p2_pairs = CollectIntersectionPairs(p2_e, p2_eBoxes, p2_bvh);
    FindAndInsertIntersections(p2_e, &p2_mrg.verts, &p2_l, &p2_ve, eps,
                               p2_eBoxes, p2_bvh, p2_pairs);
    std::cerr << "  intersections: " << p2_mrg.verts.size()
              << " verts after\n";
    auto p2_canon = Canonicalize(p2_e, p2_l);
    std::cerr << "  canonicalization: " << p2_canon.edges.size()
              << " sub-edges\n";

    auto pass2 = RemoveOverlaps2D(pass1.verts, p2in, eps, /*debug=*/true);
    std::cerr << "pass1: " << pass1.verts.size() << " verts, "
              << pass1.edges.size() << " edges\n";
    std::cerr << "pass2: " << pass2.verts.size() << " verts, "
              << pass2.edges.size() << " edges\n";
    // Dump pass 1 edges as (qx, qy)→(qx, qy) at eps/100 quantum.
    auto qedges = [&](const OverlapResult& r) {
      const double q = eps * 0.01;
      auto qq = [q](double x) { return (int64_t)std::round(x / q); };
      std::set<std::tuple<int64_t, int64_t, int64_t, int64_t, int>> s;
      for (const auto& oe : r.edges) {
        const vec2 p0 = r.verts[oe.v0], p1 = r.verts[oe.v1];
        auto a = std::make_pair(qq(p0.x), qq(p0.y));
        auto b = std::make_pair(qq(p1.x), qq(p1.y));
        int m = oe.mult;
        if (b < a) {
          std::swap(a, b);
          m = -m;
        }
        s.emplace(a.first, a.second, b.first, b.second, m);
      }
      return s;
    };
    auto e1 = qedges(pass1), e2 = qedges(pass2);
    std::set<std::tuple<int64_t, int64_t, int64_t, int64_t, int>> only1, only2;
    std::set_difference(e1.begin(), e1.end(), e2.begin(), e2.end(),
                        std::inserter(only1, only1.end()));
    std::set_difference(e2.begin(), e2.end(), e1.begin(), e1.end(),
                        std::inserter(only2, only2.end()));
    std::cerr << "Edges in pass1 but not pass2: " << only1.size() << "\n";
    for (const auto& [ax, ay, bx, by, m] : only1) {
      std::cerr << "  (" << ax << "," << ay << ")→(" << bx << "," << by
                << ") mult=" << m << "\n";
    }
    std::cerr << "Edges in pass2 but not pass1: " << only2.size() << "\n";
    for (const auto& [ax, ay, bx, by, m] : only2) {
      std::cerr << "  (" << ax << "," << ay << ")→(" << bx << "," << by
                << ") mult=" << m << "\n";
    }
    return 0;
  }
  bool allPass = true;

  // ---- Simple sanity ----
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
    auto r = RunCase({"CCW square (no overlap)", v, e, EpsilonFromScale(1.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Boolean2D entry point: Add / Subtract / Intersect ----
  // Two CCW unit squares overlapping by a 0.5×1 strip in the middle.
  {
    manifold::Polygons a = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    manifold::Polygons b = {{{0.5, 0}, {1.5, 0}, {1.5, 1}, {0.5, 1}}};
    auto add = Boolean2D(a, b, OpType::Add);
    auto sub = Boolean2D(a, b, OpType::Subtract);
    auto isec = Boolean2D(a, b, OpType::Intersect);
    // Expected boundary edge counts (each square is 4 edges; expected
    // post-overlap: union 4 edges (rectangle 0..1.5 by 0..1 with 4 sides
    // but each side has a midpoint vert from the cross-cut so actually
    // 6 boundary edges); difference 4 edges (square 0..0.5 by 0..1 with
    // 4 sides), intersection 4 edges (square 0.5..1 by 0..1 with 4
    // sides). Allow some flex from intersection-vert insertion.
    auto edgeCount = [](const manifold::Polygons& p) {
      size_t s = 0;
      for (const auto& l : p) s += l.size();
      return s;
    };
    std::cout << "=== Boolean2D: two overlapping squares ===\n";
    std::cout << "  Add:       " << add.size() << " loop(s), "
              << edgeCount(add) << " edges\n";
    std::cout << "  Subtract:  " << sub.size() << " loop(s), "
              << edgeCount(sub) << " edges\n";
    std::cout << "  Intersect: " << isec.size() << " loop(s), "
              << edgeCount(isec) << " edges\n";
    bool ok = !add.empty() && !sub.empty() && !isec.empty() &&
              edgeCount(add) >= 4 && edgeCount(sub) >= 4 &&
              edgeCount(isec) >= 4;
    std::cout << "  Smoke: " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  // Two axis-aligned CCW squares overlapping at an L-corner. Their
  // boundaries cross at TWO perpendicular axis-aligned edge pairs
  // (A's right edge × B's bottom edge at (2,1); A's top edge × B's
  // left edge at (1,2)). The Boolean2D smoke test above offsets the
  // squares horizontally so their crossing edges are *collinear*
  // (both top edges at y=1), not perpendicular - it doesn't exercise
  // axis-aligned perpendicular intersection at all. This test does.
  {
    std::vector<vec2> v = {{0, 0}, {2, 0}, {2, 2}, {0, 2},   // square A CCW
                           {1, 1}, {3, 1}, {3, 3}, {1, 3}};  // square B CCW
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // A
        {4, 5, 1}, {5, 6, 1}, {6, 7, 1}, {7, 4, 1},  // B
    };
    auto r = RunCase(
        {"Two axis-aligned squares overlapping (L-shape union)", v, e,
         EpsilonFromScale(3.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(3.0))) allPass = false;
    // Expected union: single L-shaped loop with 8 boundary verts and
    // total area 4 + 4 − 1 = 7. Intersection insertion must detect the two
    // perpendicular axis-aligned intersections.
    manifold::Polygons in = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}},
                             {{1, 1}, {3, 1}, {3, 3}, {1, 3}}};
    auto out = Simplify(in, EpsilonFromScale(3.0));
    size_t totalEdges = 0;
    for (const auto& l : out) totalEdges += l.size();
    bool shapeOk = (out.size() == 1 && totalEdges == 8);
    std::cout << "  Polygons union: " << out.size() << " loop(s), "
              << totalEdges << " edges, "
              << (shapeOk ? "PASS" : "FAIL") << "\n";
    if (!shapeOk) allPass = false;
    std::cout << std::endl;
  }

  // Two CCW unit squares overlapping diagonally, run through all three
  // Boolean2D ops with quantitative area checks. Smith fig 6.2 in
  // spirit: a single configuration that exercises every winding rule.
  // Critically, Intersect drives the `w > 1` predicate (the central
  // 1x1 overlap region has winding 2 because both inputs cover it).
  // The existing Boolean2D smoke test (two squares offset *horizontally*)
  // produces a wind=2 strip but its boundaries are collinear, not
  // perpendicular, and the Intersect winding rule was untested against
  // axis-aligned perpendicular crossings until this case.
  {
    manifold::Polygons a = {{{0, 0}, {2, 0}, {2, 2}, {0, 2}}};   // CCW
    manifold::Polygons b = {{{1, 1}, {3, 1}, {3, 3}, {1, 3}}};   // CCW
    const double eps = EpsilonFromScale(3.0);
    auto add = Boolean2D(a, b, OpType::Add, eps);
    auto sub = Boolean2D(a, b, OpType::Subtract, eps);
    auto isec = Boolean2D(a, b, OpType::Intersect, eps);
    auto xorOp = Xor(a, b, eps);
    const double areaAdd = TotalSignedArea(add);
    const double areaSub = TotalSignedArea(sub);
    const double areaIsec = TotalSignedArea(isec);
    const double areaXor = TotalSignedArea(xorOp);
    // Expected: A area 4, B area 4, overlap 1.
    // Add: 4 + 4 - 1 = 7.  Subtract: 4 - 1 = 3.  Intersect: 1.
    // Xor (symmetric difference): (A | B) - (A & B) = 7 - 1 = 6.
    auto near = [](double a, double b) { return std::fabs(a - b) < 1e-9; };
    bool ok = near(areaAdd, 7.0) && near(areaSub, 3.0) &&
              near(areaIsec, 1.0) && near(areaXor, 6.0);
    std::cout << "=== Boolean2D area regression: diagonal squares ===\n";
    std::cout << "  Add area:       " << areaAdd << " (expect 7)\n";
    std::cout << "  Subtract area:  " << areaSub << " (expect 3)\n";
    std::cout << "  Intersect area: " << areaIsec << " (expect 1)\n";
    std::cout << "  Xor area:       " << areaXor << " (expect 6)\n";
    std::cout << "  " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = LShape();
    auto r = RunCase({"L-shape (concave, CCW)", v, e, EpsilonFromScale(2.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(2.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = SquareWithHole();
    auto r = RunCase(
        {"Square with hole (outer CCW, inner CW)", v, e, EpsilonFromScale(4.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(4.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = TwoSquares();
    auto r = RunCase({"Two disjoint CCW squares", v, e, EpsilonFromScale(3.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(3.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Polygons-API smoke tests ----
  // Verify the public Simplify wrapper round-trips simple
  // closed loops through the manifold::Polygons type.
  {
    std::cout << "=== Polygons API: square ===" << std::endl;
    manifold::Polygons in = {{{0, 0}, {1, 0}, {1, 1}, {0, 1}}};
    auto out = Simplify(in, EpsilonFromScale(1.0));
    bool ok = (out.size() == 1 && out[0].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    if (!out.empty()) std::cout << ", " << out[0].size() << " verts in loop 0";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: square with hole ===" << std::endl;
    manifold::Polygons in = {
        {{0, 0}, {4, 0}, {4, 4}, {0, 4}},          // outer CCW
        {{1, 1}, {1, 3}, {3, 3}, {3, 1}},          // inner CW (hole)
    };
    auto out = Simplify(in, EpsilonFromScale(4.0));
    // Expect two loops (outer + hole) with 4 verts each.
    bool ok = (out.size() == 2 && out[0].size() == 4 && out[1].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: bowtie keeps the wind>0 lobe ==="
              << std::endl;
    // Self-intersecting bowtie input: under Smith's `wind > 0`
    // convention, only the CCW lobe is in the boolean union; the
    // other lobe has wind = −1. So output = 1 triangle.
    manifold::Polygons in = {{{0, 0}, {2, 2}, {2, 0}, {0, 2}}};
    auto out = Simplify(in, EpsilonFromScale(2.0));
    bool ok = (out.size() == 1 && out[0].size() == 3);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }
  {
    std::cout << "=== Polygons API: two disjoint squares ===" << std::endl;
    manifold::Polygons in = {
        {{0, 0}, {1, 0}, {1, 1}, {0, 1}},          // square A (CCW)
        {{2, 0}, {3, 0}, {3, 1}, {2, 1}},          // square B (CCW)
    };
    auto out = Simplify(in, EpsilonFromScale(3.0));
    bool ok = (out.size() == 2 && out[0].size() == 4 && out[1].size() == 4);
    std::cout << "  Output: " << out.size() << " loop(s)";
    for (size_t i = 0; i < out.size(); ++i)
      std::cout << " [loop " << i << ": " << out[i].size() << " verts]";
    std::cout << ": " << (ok ? "PASS" : "FAIL") << std::endl;
    if (!ok) allPass = false;
    std::cout << std::endl;
  }

  // ---- Non-closed / non-manifold input tests ----
  // The boolean2 pipeline is designed for closed polygons. These tests probe
  // what happens at the edges of that assumption: open polylines,
  // shared verts/edges between polygons, T-junctions, etc.
  //
  // For *open* inputs (polylines, T-junctions) the per-vertex balance
  // check `CheckTopologicalValidity` doesn't apply: it compares input
  // edge contributions per vertex against output edge contributions,
  // and for open inputs the input has end-vertex balance ±1 while the
  // algorithm's correct output is empty (balance 0). So we accept
  // those as observational and only fail closed-input cases.
  enum class NMExpect { ClosedTopo, OpenObservational };
  auto runNonManifold = [&allPass](
                            const char* name, const std::vector<vec2>& v,
                            const std::vector<EdgeM>& e, double eps,
                            const char* expectation, NMExpect mode) {
    std::cout << "=== Non-manifold: " << name << " ===" << std::endl;
    std::cout << "  Input: " << v.size() << " verts, " << e.size() << " edges, "
              << "expectation: " << expectation << std::endl;
    auto r = RemoveOverlaps2D(v, e, eps);
    std::cout << "  Output: " << r.verts.size() << " verts, " << r.edges.size()
              << " edges" << std::endl;
    bool topoOk =
        CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
    if (mode == NMExpect::ClosedTopo) {
      std::cout << "  Topological validity: " << (topoOk ? "PASS" : "FAIL")
                << std::endl;
      if (!topoOk) allPass = false;
    } else {
      std::cout << "  Topological validity: " << (topoOk ? "PASS" : "FAIL")
                << " (observational; open input, balance check N/A)"
                << std::endl;
    }
    std::cout << std::endl;
    return r;
  };

  // (1) Open polyline (zigzag): no closure, no enclosed region.
  // Expected: empty output (wind > 0 is empty everywhere).
  {
    std::vector<vec2> v = {{0, 0}, {1, 1}, {2, 0}, {3, 1}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
    runNonManifold("open polyline (Z-zigzag)", v, e, EpsilonFromScale(3.0),
                   "empty (no enclosed region)",
                   NMExpect::OpenObservational);
  }

  // (2) Two polygons touching at a single vertex (kissing corner).
  // Squares share v=(1,1); each is otherwise CCW and disjoint.
  // Expected: 2 disjoint loops (8 output edges total) + the shared
  // vertex appears in both.
  {
    std::vector<vec2> v = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1},  // square A
        {2, 1}, {2, 2}, {1, 2},          // square B (shares (1,1))
    };
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // square A CCW
        {2, 4, 1}, {4, 5, 1}, {5, 6, 1}, {6, 2, 1},  // square B CCW (via (1,1))
    };
    runNonManifold("two polys touching at a vertex", v, e,
                   EpsilonFromScale(2.0), "2 disjoint loops",
                   NMExpect::ClosedTopo);
  }

  // (3) Two polygons sharing an edge (non-manifold along the shared
  // segment, two CCW squares glued along (1,0)→(1,1)). The shared
  // edge has canonical mult = +2 in the combined input. Expected: 1
  // outer loop = the rectangle [0,2]×[0,1] (the two squares fused).
  {
    std::vector<vec2> v = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1},  // square A
        {2, 0}, {2, 1},                  // square B's other corners
    };
    std::vector<EdgeM> e = {
        {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},  // square A CCW
        {1, 4, 1}, {4, 5, 1}, {5, 2, 1}, {2, 1, 1},  // square B CCW
    };
    runNonManifold("two polys sharing an edge", v, e, EpsilonFromScale(2.0),
                   "1 fused outer loop", NMExpect::ClosedTopo);
  }

  // (4) T-junction: vertex with degree 3, no closure. Three edges
  // dangling from a common point. No enclosed region.
  // Expected: empty output.
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {0, 1}, {-1, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {0, 2, 1}, {0, 3, 1}};
    runNonManifold("T-junction (degree-3 dangling)", v, e,
                   EpsilonFromScale(1.0), "empty",
                   NMExpect::OpenObservational);
  }

  // (5) Open polyline crossing itself (figure-8 missing the closing
  // edge). Topology-wise: 5 verts, 4 edges, one self-crossing.
  // Expected: empty output (no closed boundary).
  {
    std::vector<vec2> v = {{0, 0}, {2, 2}, {2, 0}, {0, 2}, {0, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}};
    // Note: although the input has 5 verts forming an *open* zigzag,
    // verts 0 and 4 share the same coordinate (0,0). Vertex merging
    // collapses them, turning the input into a *closed* bowtie. So
    // this case actually exercises the closed bowtie path, not an
    // open-polyline path. The topology check applies.
    runNonManifold("polyline w/ duplicate-position endpoints (=> closed bowtie)",
                   v, e, EpsilonFromScale(2.0), "1 triangle (= bowtie lobe)",
                   NMExpect::ClosedTopo);
  }

  // (6) Degenerate triangle: three collinear verts (zero-area input).
  // Expected: empty (no enclosed region).
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {2, 0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    runNonManifold("collinear triangle (zero area)", v, e,
                   EpsilonFromScale(2.0), "empty", NMExpect::ClosedTopo);
  }

  // (7) Duplicate-vertex polygon: same coordinate listed twice. Vertex
  // merging should fold them; edge collapse should drop the resulting
  // self-loop edge.
  // Expected: 1 valid loop (the remaining triangle).
  {
    std::vector<vec2> v = {{0, 0}, {1, 0}, {0, 1}, {0, 1}};  // v2 == v3
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1}};
    runNonManifold("polygon with duplicate vertex", v, e, EpsilonFromScale(1.0),
                   "1 triangle (after merge)", NMExpect::ClosedTopo);
  }

  // (8) Adversarial 4+ concurrent edges: 4 line segments all passing
  // through origin at different orientations. Intersection insertion produces
  // C(4, 2) = 6 pairwise intersections that should all snap to one
  // point in the structural merge. Some pairs share no edge in their incidence
  // sets (e.g. seg AB×CD vs seg EF×GH share no input edge), so that merge's
  // structural gate cannot merge them directly. The iterate-to-
  // fixed-point pass picks up any leftover near-duplicates. Open
  // input → empty output (no enclosed region).
  {
    std::vector<vec2> v = {{-1, 0},     {1, 0},      {0, -1},     {0, 1},
                           {-0.7, -0.7}, {0.7, 0.7},  {-0.7, 0.7}, {0.7, -0.7}};
    std::vector<EdgeM> e = {{0, 1, 1}, {2, 3, 1}, {4, 5, 1}, {6, 7, 1}};
    runNonManifold("4 segments concurrent at origin (asterisk)", v, e,
                   EpsilonFromScale(1.0), "empty (open input)",
                   NMExpect::OpenObservational);
  }

  // (9) Near-equal-coords at displaced scale: two triangle verts
  // separated by a few ULPs at 2^49 displacement. Eps-radius vertex merging
  // should fold them, collapsing the triangle to empty output
  // (degenerate after merge).
  {
    const double offset = std::ldexp(1.5, 49);
    const double ulp = std::ldexp(1.0, 49 - 52);  // ~1 ULP at offset
    std::vector<vec2> v = {{offset + 1.0, offset},
                           {offset + 1.0 + ulp, offset},
                           {offset, offset + 1.0}};
    std::vector<EdgeM> e = {{0, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    runNonManifold("near-equal coords at 2^49 displacement", v, e,
                   EpsilonFromScale(offset), "empty (degenerate after merge)",
                   NMExpect::ClosedTopo);
  }

  // ---- Self-intersecting cases ----
  {
    std::vector<vec2> v = {{0, 0}, {2, 0}, {0, 2}, {2, 2}};
    std::vector<EdgeM> e = {{0, 3, 1}, {3, 1, 1}, {1, 2, 1}, {2, 0, 1}};
    auto r = RunCase(
        {"Bowtie (one interior intersection)", v, e, EpsilonFromScale(2.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(2.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = RandomTopologicalPolygon(8, 42);
    auto r = RunCase(
        {"Random topological polygon (n=8)", v, e, EpsilonFromScale(1.0)},
        &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = PolygonalStar(7, 3);
    auto r = RunCase({"Heptagram (n=7, skip=3)", v, e, EpsilonFromScale(1.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(1.0))) allPass = false;
    std::cout << std::endl;
  }

  {
    auto [v, e] = PolygonalStar(7, 3);
    const double offset = std::ldexp(1.5, 49);
    Displace(&v, offset);
    RunCase({"Heptagram displaced by 1.5*2^49", v, e, EpsilonFromScale(offset)},
            &allPass);
    std::cout << std::endl;
  }

  // ---- Smith §7.2 figure 6.5(d) variants ----
  // Smith catalogs two inconsistent-edge sets that exhibit the cyclic-
  // ordering failure: {AF, BG, CE} and {AG, BF, CE}. Both embedded in
  // closed self-intersecting hexagons.
  const double u = std::ldexp(1.0, -53);
  std::vector<vec2> smithVerts = {
      {0.5 + 9 * u, 0.5 + 23 * u},   // A=0
      {0.5 + 23 * u, 0.5 + 25 * u},  // B=1
      {9.0, 0.5},                    // C=2
      {9.0, 32.0},                   // E=3
      {32.0, 32.0 - 32 * u},         // F=4
      {32.0, 32.0},                  // G=5
  };
  {
    // A→F→C→E→B→G→A (uses AF, BG, with CE inside)
    std::vector<EdgeM> e = {{0, 4, 1}, {4, 2, 1}, {2, 3, 1},
                            {3, 1, 1}, {1, 5, 1}, {5, 0, 1}};
    auto r = RunCase({"Smith fig 6.5(d) hexagon {AF, BG, CE}", smithVerts, e,
                      EpsilonFromScale(32.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(32.0))) allPass = false;
    std::cout << std::endl;
  }
  {
    // A→G→C→E→B→F→A (uses AG, BF, with CE inside)
    std::vector<EdgeM> e = {{0, 5, 1}, {5, 2, 1}, {2, 3, 1},
                            {3, 1, 1}, {1, 4, 1}, {4, 0, 1}};
    auto r = RunCase({"Smith fig 6.5(d) hexagon {AG, BF, CE}", smithVerts, e,
                      EpsilonFromScale(32.0)},
                     &allPass);
    if (!CheckIdempotence(r, EpsilonFromScale(32.0))) allPass = false;
    std::cout << std::endl;
  }

  // ---- Fuzz: random topological polygons at multiple seeds and sizes ----
  std::cout << "=== Fuzz: random topological polygons ===" << std::endl;
  int fuzzPass = 0, fuzzFail = 0;
  for (int n : {5, 8, 12, 20, 50, 100, 200}) {
    for (uint64_t seed = 0; seed < 50; ++seed) {
      auto [v, e] = RandomTopologicalPolygon(n, seed);
      auto r = RemoveOverlaps2D(v, e, EpsilonFromScale(1.0));
      bool ok = CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
      if (ok) {
        ++fuzzPass;
      } else {
        ++fuzzFail;
        std::cout << "  FAIL: n=" << n << " seed=" << seed << std::endl;
      }
    }
  }
  std::cout << "  Topological invariant: " << fuzzPass << " pass, " << fuzzFail
            << " fail (7 sizes × 50 seeds = 350 cases)" << std::endl;
  if (fuzzFail > 0) allPass = false;

  // ---- Fuzz: displaced random polygons (Smith §7.7 displacement attack) ----
  std::cout << "=== Fuzz: random polygons displaced near 2^k boundaries ==="
            << std::endl;
  int dispPass = 0, dispFail = 0;
  int singlePassIdem = 0, iterConverged = 0, iterCycled = 0, iterMaxedOut = 0;
  std::map<int, int> iterDistribution;  // iter count -> count
  for (int kPow : {10, 20, 30, 40, 49}) {
    const double offset = std::ldexp(1.5, kPow);
    const double eps = EpsilonFromScale(offset);
    for (int n : {8, 20, 50}) {
      for (uint64_t seed = 0; seed < 30; ++seed) {
        auto [v, e] = RandomTopologicalPolygon(n, seed + 1000 * kPow);
        Displace(&v, offset);
        // Single-pass: topological validity only.
        auto r = RemoveOverlaps2D(v, e, eps);
        bool ok =
            CheckTopologicalValidity(r, e, r.inputRemap, r.numMergedVerts);
        if (ok) {
          ++dispPass;
        } else {
          ++dispFail;
          std::cout << "  FAIL: kPow=" << kPow << " n=" << n << " seed=" << seed
                    << std::endl;
        }
        // Single-pass idempotence: re-run and compare fingerprint.
        std::vector<EdgeM> rEdges;
        for (const auto& oe : r.edges)
          rEdges.push_back({oe.v0, oe.v1, oe.mult});
        auto r2 = RemoveOverlaps2D(r.verts, rEdges, eps);
        if (Fingerprint(r, eps) == Fingerprint(r2, eps)) ++singlePassIdem;
        // Iterate-to-fixed-point.
        int iters = 0;
        IterStatus status = IterStatus::Converged;
        auto rIter =
            IterateToFixedPoint(v, e, eps, /*maxIter=*/8, &iters, &status);
        ++iterDistribution[iters];
        switch (status) {
          case IterStatus::Converged:
            ++iterConverged;
            break;
          case IterStatus::Cycled:
            ++iterCycled;
            break;
          case IterStatus::MaxedOut:
            ++iterMaxedOut;
            break;
        }
        // Topological validity on the iterated result. NOTE: rIter.inputRemap
        // is the LAST iteration's input->output remap, not the composed
        // original->final. Works for the random-topology displaced inputs
        // here because vertex merging produces an identity remap when no input verts
        // merge (true for these inputs); a production version of
        // IterateToFixedPoint should compose remaps across iterations.
        if (!CheckTopologicalValidity(rIter, e, rIter.inputRemap,
                                      rIter.numMergedVerts)) {
          std::cout << "  ITER FAIL: kPow=" << kPow << " n=" << n
                    << " seed=" << seed << std::endl;
        }
      }
    }
  }
  std::cout << "  Single-pass idempotence: " << singlePassIdem << "/450"
            << std::endl;
  std::cout << "  Iterate-to-fixed-point: " << iterConverged << " converged, "
            << iterCycled << " cycled, " << iterMaxedOut << " hit max-iter (8)"
            << std::endl;
  std::cout << "  Iter distribution: ";
  for (const auto& [k, v] : iterDistribution) {
    std::cout << k << ":" << v << " ";
  }
  std::cout << std::endl;
  std::cout << "  Displacement invariant: " << dispPass << " pass, " << dispFail
            << " fail (5 scales × 3 sizes × 30 seeds = 450 cases)" << std::endl;
  if (dispFail > 0) allPass = false;
  std::cout << std::endl;

  std::cout << "==== OVERALL: " << (allPass ? "PASS" : "FAIL")
            << " ====" << std::endl;
  return allPass ? 0 : 1;
}
