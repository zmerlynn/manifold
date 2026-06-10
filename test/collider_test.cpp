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

#include "../src/collider.h"

#include <limits>
#include <vector>

#include "test.h"

using namespace manifold;

namespace {

// Collects (queryIdx, leafIdx) hits sequentially.
struct HitCollector {
  std::vector<std::pair<int, int>> hits;

  void operator()(int queryIdx, int leafIdx) {
    hits.push_back({queryIdx, leafIdx});
  }
};

Collider MakeCollider(const std::vector<Box>& boxes) {
  Box outer;
  for (const Box& b : boxes) outer = outer.Union(b);
  Vec<Box> leafBB;
  Vec<uint32_t> leafMorton;
  for (const Box& b : boxes) {
    leafBB.push_back(b);
    leafMorton.push_back(Collider::MortonCode(b.Center(), outer));
  }
  return Collider(leafBB, leafMorton);
}

}  // namespace

// Trees with zero and one leaf have no internal nodes; the radix
// traversal, box propagation, and bounding-box root lookup all need the
// degenerate forms. These pin them at the Collider level (the house
// mesh paths never construct fewer than 4 leaves).

TEST(Collider, EmptyTree) {
  Collider collider;
  HitCollector f;
  auto recorder = MakeSimpleRecorder(f);
  Vec<Box> queries;
  queries.push_back(Box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  collider.Collisions<false>(recorder, queries.cview(), /*parallel=*/false);
  EXPECT_TRUE(f.hits.empty());
  const Box bb = collider.GetBoundingBox();
  EXPECT_EQ(bb.min.x, std::numeric_limits<double>::infinity());
}

TEST(Collider, SingleLeafCollisions) {
  Collider collider = MakeCollider({Box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})});

  HitCollector f;
  auto recorder = MakeSimpleRecorder(f);
  Vec<Box> queries;
  queries.push_back(Box({2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}));     // miss
  queries.push_back(Box({0.5, 0.5, 0.5}, {2.0, 2.0, 2.0}));     // hit
  queries.push_back(Box({-1.0, -1.0, -1.0}, {0.0, 0.0, 0.0}));  // touch = hit
  collider.Collisions<false>(recorder, queries.cview(), /*parallel=*/false);

  ASSERT_EQ(f.hits.size(), 2u);
  EXPECT_EQ(f.hits[0], (std::pair<int, int>(1, 0)));
  EXPECT_EQ(f.hits[1], (std::pair<int, int>(2, 0)));
}

TEST(Collider, SingleLeafBoundingBox) {
  const Box leaf({-2.0, 1.0, 0.5}, {3.0, 4.0, 2.5});
  Collider collider = MakeCollider({leaf});
  const Box bb = collider.GetBoundingBox();
  EXPECT_EQ(bb.min, leaf.min);
  EXPECT_EQ(bb.max, leaf.max);
}

TEST(Collider, SingleLeafSelfCollision) {
  const Box leaf({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  Collider collider = MakeCollider({leaf});
  HitCollector f;
  auto recorder = MakeSimpleRecorder(f);
  Vec<Box> queries;
  queries.push_back(leaf);
  collider.Collisions<true>(recorder, queries.cview(), /*parallel=*/false);
  // selfCollision skips queryIdx == leafIdx: the lone leaf against
  // itself records nothing.
  EXPECT_TRUE(f.hits.empty());
}

TEST(Collider, SingleLeafUpdateBoxes) {
  Collider collider = MakeCollider({Box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0})});
  Vec<Box> moved;
  moved.push_back(Box({10.0, 10.0, 10.0}, {11.0, 11.0, 11.0}));
  collider.UpdateBoxes(moved.cview());

  const Box bb = collider.GetBoundingBox();
  EXPECT_EQ(bb.min, moved[0].min);
  EXPECT_EQ(bb.max, moved[0].max);

  HitCollector f;
  auto recorder = MakeSimpleRecorder(f);
  Vec<Box> queries;
  queries.push_back(Box({0.5, 0.5, 0.5}, {2.0, 2.0, 2.0}));        // old spot
  queries.push_back(Box({10.5, 10.5, 10.5}, {12.0, 12.0, 12.0}));  // new spot
  collider.Collisions<false>(recorder, queries.cview(), /*parallel=*/false);
  ASSERT_EQ(f.hits.size(), 1u);
  EXPECT_EQ(f.hits[0], (std::pair<int, int>(1, 0)));
}

TEST(Collider, TwoLeavesStillTraverse) {
  // Smallest tree with an internal node: the radix path, not the
  // single-leaf arm.
  Collider collider = MakeCollider({Box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}),
                                    Box({5.0, 0.0, 0.0}, {6.0, 1.0, 1.0})});
  HitCollector f;
  auto recorder = MakeSimpleRecorder(f);
  Vec<Box> queries;
  queries.push_back(Box({5.5, 0.5, 0.5}, {7.0, 2.0, 2.0}));  // second leaf only
  collider.Collisions<false>(recorder, queries.cview(), /*parallel=*/false);
  ASSERT_EQ(f.hits.size(), 1u);
  EXPECT_EQ(f.hits[0].first, 0);
  EXPECT_EQ(f.hits[0].second, 1);
}
