// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <numeric>
#include "latch/LatchDependencyGraph.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

Primitive cell(std::vector<size_t> inputs, std::vector<size_t> outputs) {
  Primitive result;
  result.inputs = std::move(inputs);
  result.outputs = std::move(outputs);
  return result;
}

TEST(LatchDependencyGraphTests, EmptyNetworkHasNoComponents) {
  const auto graph = analyzeDependencies({});
  EXPECT_TRUE(graph.components.empty());
  EXPECT_TRUE(graph.feedbackComponents.empty());
  EXPECT_TRUE(graph.componentOf.empty());
}

TEST(LatchDependencyGraphTests, AcyclicChainIsOneSchedulingIslandWithoutFeedback) {
  const auto graph = analyzeDependencies({4, {0},
      {cell({0}, {1}), cell({1}, {2}), cell({2}, {3})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1, 2}}));
  EXPECT_TRUE(graph.feedbackComponents.empty());
  EXPECT_EQ(graph.componentOf, (std::vector<size_t>{0, 0, 0}));
}

TEST(LatchDependencyGraphTests, InterconnectedLoopsRemainDistinctInsideOneIsland) {
  const auto graph = analyzeDependencies({5, {}, {cell({1}, {0}), cell({0}, {1}),
      cell({1, 3}, {2}), cell({2}, {3}), cell({3}, {4})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1, 2, 3, 4}}));
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0, 1}, {2, 3}}));
}

TEST(LatchDependencyGraphTests, SelfFeedbackIsAnExplicitFeedbackComponent) {
  const auto graph = analyzeDependencies({2, {}, {cell({0}, {0}), cell({}, {1})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0}, {1}}));
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0}}));
}

TEST(LatchDependencyGraphTests, EveryInputIncludingControlClosesTheGraph) {
  // The second pin may be a clock, enable or async reset. Its role is irrelevant
  // to event connectivity; cutting it would lose transient control propagation.
  const auto graph = analyzeDependencies({4, {0, 1},
      {cell({3}, {2}), cell({0, 2, 1}, {3})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1}}));
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0, 1}}));
}

TEST(LatchDependencyGraphTests, SequentialStorageIsNotAnImplicitGraphCut) {
  auto first = cell({1}, {0});
  auto second = cell({0}, {1});
  first.storageBits = second.storageBits = 1;
  const auto graph = analyzeDependencies({2, {}, {first, second}});
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0, 1}}));
}

TEST(LatchDependencyGraphTests, SharedReadOnlyInputsDoNotJoinIndependentIslands) {
  const auto graph = analyzeDependencies({4, {0, 1},
      {cell({0, 1}, {2}), cell({0, 1}, {3})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0}, {1}}));
  EXPECT_TRUE(graph.feedbackComponents.empty());
  EXPECT_EQ(graph.componentOf, (std::vector<size_t>{0, 1}));
}

TEST(LatchDependencyGraphTests, MultipleWritersJoinWithoutInventingDirectedFeedback) {
  const auto graph = analyzeDependencies({4, {0, 1},
      {cell({0}, {2}), cell({1}, {2}), cell({2}, {3})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1, 2}}));
  EXPECT_TRUE(graph.feedbackComponents.empty());
}

TEST(LatchDependencyGraphTests, MultipleUnconsumedWritersStillShareAnIsland) {
  const auto graph = analyzeDependencies({3, {0, 1},
      {cell({0}, {2}), cell({1}, {2})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1}}));
  EXPECT_TRUE(graph.feedbackComponents.empty());
}

TEST(LatchDependencyGraphTests, DuplicatePinArcsDoNotChangeComponents) {
  const auto graph = analyzeDependencies({2, {},
      {cell({1, 1}, {0, 0}), cell({0, 0}, {1, 1})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 1}}));
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0, 1}}));
}

TEST(LatchDependencyGraphTests, IsolatedAndOutputlessCellsAreRetained) {
  const auto graph = analyzeDependencies({2, {0},
      {cell({}, {}), cell({0}, {1}), cell({1}, {})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0}, {1, 2}}));
  EXPECT_EQ(graph.componentOf, (std::vector<size_t>{0, 1, 1}));
}

TEST(LatchDependencyGraphTests, ComponentsAreSortedByLowestPrimitiveIndex) {
  const auto graph = analyzeDependencies({5, {}, {cell({4}, {0}), cell({3}, {1}),
      cell({}, {2}), cell({1}, {3}), cell({0}, {4})}});
  EXPECT_EQ(graph.components, (std::vector<std::vector<size_t>>{{0, 4}, {1, 3}, {2}}));
  EXPECT_EQ(graph.feedbackComponents, (std::vector<std::vector<size_t>>{{0, 4}, {1, 3}}));
  EXPECT_EQ(graph.componentOf, (std::vector<size_t>{0, 1, 2, 1, 0}));
}

TEST(LatchDependencyGraphTests, InvalidNetIndicesAreRejected) {
  EXPECT_THROW(analyzeDependencies(Network{1, {1}, {}}), std::invalid_argument);
  EXPECT_THROW(analyzeDependencies(Network{1, {}, {cell({1}, {})}}), std::invalid_argument);
  EXPECT_THROW(analyzeDependencies(Network{1, {}, {cell({}, {1})}}), std::invalid_argument);
}

TEST(LatchDependencyGraphTests, TenThousandCellChainUsesNoRecursiveTraversal) {
  constexpr size_t length = 10000;
  Network network;
  network.netCount = length + 1;
  network.externalInputs = {0};
  for (size_t i = 0; i < length; ++i) network.primitives.push_back(cell({i}, {i + 1}));
  auto graph = analyzeDependencies(network);
  ASSERT_EQ(graph.components.size(), 1);
  EXPECT_EQ(graph.components.front().size(), length);
  EXPECT_TRUE(graph.feedbackComponents.empty());
  // Close the same long chain to exercise the reverse SCC traversal as well.
  network.primitives.front().inputs = {length};
  graph = analyzeDependencies(network);
  ASSERT_EQ(graph.feedbackComponents.size(), 1);
  EXPECT_EQ(graph.feedbackComponents.front().size(), length);
}

TEST(LatchDependencyGraphTests, AllThreeNodeDirectedGraphsMatchTransitiveClosure) {
  constexpr size_t size = 3;
  for (size_t mask = 0; mask < (size_t{1} << (size * size)); ++mask) {
    SCOPED_TRACE(mask);
    Network network;
    network.netCount = size;
    bool directed[size][size] = {}, weak[size][size] = {};
    for (size_t target = 0; target < size; ++target) {
      auto primitive = cell({}, {target});
      for (size_t source = 0; source < size; ++source) {
        if (!(mask & (size_t{1} << (source * size + target)))) continue;
        primitive.inputs.push_back(source);
        directed[source][target] = true;
        weak[source][target] = weak[target][source] = true;
      }
      network.primitives.push_back(std::move(primitive));
    }
    for (size_t k = 0; k < size; ++k)
      for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j) {
          directed[i][j] |= directed[i][k] && directed[k][j];
          weak[i][j] |= weak[i][k] && weak[k][j];
        }
    std::vector<std::vector<size_t>> expectedFeedback, expectedIslands;
    bool assignedFeedback[size] = {}, assignedIsland[size] = {};
    for (size_t i = 0; i < size; ++i) {
      if (!assignedIsland[i]) {
        std::vector<size_t> group;
        for (size_t j = i; j < size; ++j)
          if (i == j || weak[i][j]) { group.push_back(j); assignedIsland[j] = true; }
        expectedIslands.push_back(std::move(group));
      }
      if (directed[i][i] && !assignedFeedback[i]) {
        std::vector<size_t> group;
        for (size_t j = i; j < size; ++j)
          if (directed[i][j] && directed[j][i]) { group.push_back(j); assignedFeedback[j] = true; }
        expectedFeedback.push_back(std::move(group));
      }
    }
    const auto graph = analyzeDependencies(network);
    EXPECT_EQ(graph.components, expectedIslands);
    EXPECT_EQ(graph.feedbackComponents, expectedFeedback);
  }
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC::LATCH
