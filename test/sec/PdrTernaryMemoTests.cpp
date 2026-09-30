// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "BoolExprCache.h"
#include "pdr/PDREngine.h"

namespace KEPLER_FORMAL::SEC {
namespace {

class PdrTernaryMemoTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  KInductionProblem twoMappedRoots(bool largerFirst) {
    KInductionProblem problem;
    problem.state0Symbols = {2, 3};
    problem.state1Symbols = {4, 5};
    problem.inputSymbols = {6, 7};
    problem.allSymbols = {2, 3, 4, 5, 6, 7};
    problem.initialStateAssignments = {{2, false}, {3, true}, {4, false}, {5, true}};
    problem.initialCondition = BoolExpr::createTrue();
    problem.initializedStateCount = problem.totalStateCount = 4;
    problem.bad = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(4));
    problem.property = BoolExpr::Not(problem.bad);
    problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
    auto& lazy = *problem.lazyTransitions;
    lazy.localToCombinedByDesign[0] = {{10, 3}, {11, 6}, {12, 7}};
    lazy.localToCombinedByDesign[1] = {{10, 5}, {11, 6}, {12, 7}};
    auto* small = BoolExpr::Or(BoolExpr::Var(10), BoolExpr::Var(11));
    auto* large = BoolExpr::Or(BoolExpr::Var(10),
        BoolExpr::And(BoolExpr::Var(11), BoolExpr::Var(12)));
    lazy.sourceByStateSymbol.emplace(2, LazyTransitionSource{
        0, largerFirst ? large : small, LazyTransitionRail::Binary});
    lazy.sourceByStateSymbol.emplace(4, LazyTransitionSource{
        1, largerFirst ? small : large, LazyTransitionRail::Binary});
    lazy.sourceByStateSymbol.emplace(3, LazyTransitionSource{
        0, BoolExpr::Var(10), LazyTransitionRail::Binary});
    lazy.sourceByStateSymbol.emplace(5, LazyTransitionSource{
        1, BoolExpr::Var(10), LazyTransitionRail::Binary});
    return problem;
  }
};

TEST_F(PdrTernaryMemoTests, LaterMappedRootCanAddParentsBeyondEarlierMemoSize) {
  // Both design-local DAGs share variable nodes but use different symbol maps.
  // Compiling the second root adds parents to those shared nodes after the
  // first map's memo was allocated. Probing a predecessor literal must never
  // index that shorter memo with a later root's parent index (ASan regression).
  for (bool largerFirst : {false, true}) {
    SCOPED_TRACE(largerFirst);
    const auto problem = twoMappedRoots(largerFirst);
    auto cache = std::make_shared<PDRExactInitCache>(problem, Config::SolverType::KISSAT);
    PDREngine engine(problem, Config::SolverType::KISSAT, 0, cache);
    for (size_t repetition = 0; repetition < 3; ++repetition) {
      const auto result = engine.run(2);
      EXPECT_EQ(result.status, PDRStatus::Different);
      EXPECT_EQ(result.bound, 1u);
    }
  }
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
