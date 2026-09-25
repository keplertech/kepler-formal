// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BoolExprCache.h"
#include "export/SecBtor2Exporter.h"
#include "imc/CraigInterpolatingModelChecker.h"
#include "imc/IMCEngine.h"
#include "kinduction/KInductionEngine.h"
#include "pdr/PDREngine.h"
#include "proof/ProofEngineShared.h"

namespace KEPLER_FORMAL::SEC {
namespace {

BoolExpr* equal(size_t lhs, size_t rhs) {
  return BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs)));
}

void setBad(KInductionProblem& problem, BoolExpr* bad) {
  problem.bad = problem.inductionBad = bad;
  problem.property = problem.inductionProperty = BoolExpr::Not(bad);
  problem.observedOutputNames = {"result"};
  problem.observedOutputExprs0 = {bad};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
}

KInductionProblem relationalProblem() {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.auxiliaryStateSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  problem.transitions1 = {{3, BoolExpr::Var(3)}};
  problem.auxiliaryTransitions = {{4, BoolExpr::Var(4)}};
  // A unit fact must not replace the independently supplied relation. Counts
  // describe known bits, not the completeness of this set-valued BOOT state.
  problem.initialStateAssignments = {{4, false}};
  problem.initializedStateCount = 1;
  problem.totalStateCount = 3;
  problem.initialCondition = equal(2, 3);
  problem.hasExactRelationalInitialState = true;
  setBad(problem, BoolExpr::Or(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3)),
                              BoolExpr::Var(4)));
  return problem;
}

class RelationalInitializationTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
};

class RelationalInitializationEngineTests : public RelationalInitializationTests,
    public ::testing::WithParamInterface<int> {
 protected:
  void check(const KInductionProblem& problem, bool different, size_t bound = 0) {
    constexpr auto solver = Config::SolverType::KISSAT;
    if (GetParam() == 0) {
      const auto result = KInductionEngine(problem, solver).run(5);
      EXPECT_EQ(result.status, different ? KInductionStatus::Different : KInductionStatus::Equivalent);
      if (different) EXPECT_EQ(result.bound, bound);
    } else if (GetParam() == 1) {
      const auto result = PDREngine(problem, solver).run(5);
      EXPECT_EQ(result.status, different ? PDRStatus::Different : PDRStatus::Equivalent);
      if (different) EXPECT_EQ(result.bound, bound);
    } else {
      const auto result = IMCEngine(problem, solver).run(5);
      EXPECT_EQ(result.status, different ? IMCStatus::Different : IMCStatus::Equivalent);
      if (different) EXPECT_EQ(result.bound, bound);
    }
  }
};

TEST_P(RelationalInitializationEngineTests, RelationAndUnitFactsBothConstrainBoot) {
  check(relationalProblem(), false);
}

TEST_P(RelationalInitializationEngineTests, MismatchingBootIsNotAssumedAway) {
  auto problem = relationalProblem();
  problem.initialCondition = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  check(problem, true, 0);
}

TEST_P(RelationalInitializationEngineTests, ExactUnconstrainedBootDoesNotAssumeOutputEquality) {
  auto problem = relationalProblem();
  problem.initialCondition = BoolExpr::createTrue();
  problem.initialStateAssignments.clear();
  problem.initializedStateCount = 0;
  check(problem, true, 0);
}

TEST_P(RelationalInitializationEngineTests, InitializationRelationIsNotAnInvariant) {
  auto problem = relationalProblem();
  problem.transitions1 = {{3, BoolExpr::Not(BoolExpr::Var(3))}};
  check(problem, true, 1);
}

TEST_P(RelationalInitializationEngineTests, InitialPredicateSupportExtendsOutputCoi) {
  auto problem = relationalProblem();
  // The output depends only on 2; the initial relation must bring 4 into the
  // initial COI so its separate unit fact can constrain 2.
  problem.initialCondition = equal(2, 4);
  setBad(problem, BoolExpr::Var(2));
  check(problem, false);
}

TEST_P(RelationalInitializationEngineTests, OutputSubsetsRetainBootRelation) {
  auto problem = relationalProblem();
  auto* different = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.observedOutputNames = {"relation", "unit"};
  problem.observedOutputExprs0 = {different, BoolExpr::Var(4)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse(), BoolExpr::createFalse()};
  check(problem, false);
  problem.transitions1 = {{3, BoolExpr::Not(BoolExpr::Var(3))}};
  check(problem, true, 1);
}

TEST_P(RelationalInitializationEngineTests, SharedInputOriginAndIndependentStorageOriginsDiffer) {
  for (bool shared : {true, false}) {
    SCOPED_TRACE(shared);
    auto problem = relationalProblem();
    problem.auxiliaryStateSymbols = {4, 5, 6};
    problem.auxiliaryTransitions = {{4, BoolExpr::Var(4)}, {5, BoolExpr::Var(5)},
                                    {6, BoolExpr::Var(6)}};
    problem.allSymbols = {2, 3, 4, 5, 6};
    problem.totalStateCount = 5;
    problem.initialStateAssignments = {{6, false}};
    problem.initialCondition = BoolExpr::And(equal(2, 4), equal(3, shared ? 4 : 5));
    setBad(problem, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3)));
    check(problem, !shared, 0);
  }
}

INSTANTIATE_TEST_SUITE_P(KiPdrImc, RelationalInitializationEngineTests,
                        ::testing::Values(0, 1, 2));

TEST_F(RelationalInitializationTests, SharedProofInitKeepsUnitsWithoutInventingPointState) {
  auto problem = relationalProblem();
  auto* initial = buildProofInitFormula(problem);
  ASSERT_NE(initial, nullptr);
  for (bool value : {false, true}) {
    EXPECT_TRUE(initial->evaluate({{2, value}, {3, value}, {4, false}}));
    EXPECT_FALSE(initial->evaluate({{2, value}, {3, !value}, {4, false}}));
    EXPECT_FALSE(initial->evaluate({{2, value}, {3, value}, {4, true}}));
  }
}

TEST_F(RelationalInitializationTests, CraigProjectionPreservesRelationAndItsHiddenOrigin) {
  auto problem = relationalProblem();
  problem.initialCondition = equal(2, 4);
  setBad(problem, BoolExpr::Var(2));
  const auto result = CraigInterpolatingModelChecker(problem).run(4);
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(RelationalInitializationTests, ImcCounterexampleStartsAtExactBootRelation) {
  auto problem = relationalProblem();
  problem.initialCondition = BoolExpr::And(equal(2, 4), equal(3, 4));
  problem.transitions0 = {{2, BoolExpr::Var(3)}};
  problem.transitions1 = {{3, BoolExpr::createTrue()}};
  setBad(problem, BoolExpr::Var(2));
  const auto result = IMCEngine(problem, Config::SolverType::KISSAT).run(4);
  EXPECT_EQ(result.status, IMCStatus::Different);
  EXPECT_EQ(result.bound, 2);
}

TEST_F(RelationalInitializationTests, ImcExactReachabilityIncludesAuxiliaryTransitions) {
  auto problem = relationalProblem();
  // Three-bit rotating one-hot state. The public property alone is not
  // inductive; neither is the one-step image. The exact reachability path must
  // retain the auxiliary transition to close the three-state orbit.
  problem.initialStateAssignments = {{2, true}, {3, false}, {4, false}};
  problem.initializedStateCount = 3;
  problem.initialCondition = BoolExpr::createTrue();
  problem.transitions0 = {{2, BoolExpr::Var(3)}};
  problem.transitions1 = {{3, BoolExpr::Var(4)}};
  problem.auxiliaryTransitions = {{4, BoolExpr::Var(2)}};
  setBad(problem, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  const auto result = IMCEngine(problem, Config::SolverType::KISSAT).run(4);
  EXPECT_EQ(result.status, IMCStatus::Equivalent);
}

TEST_F(RelationalInitializationTests, LargeDualRailImcKeepsBooleanRelations) {
  auto problem = relationalProblem();
  problem.usesDualRailStateEncoding = true;
  // Select Craig IMC's large-state path; the arbitrary origins are Boolean,
  // not X=11 values or fabricated fixed initial assignments.
  for (size_t symbol = 5; symbol != 18; ++symbol) {
    problem.auxiliaryStateSymbols.push_back(symbol);
    problem.auxiliaryTransitions.emplace_back(symbol, BoolExpr::Var(symbol));
    problem.allSymbols.push_back(symbol);
  }
  problem.totalStateCount = problem.combinedStateSymbols().size();
  const auto result = IMCEngine(problem, Config::SolverType::KISSAT).run(4);
  EXPECT_EQ(result.status, IMCStatus::Equivalent);
}

// Independent exhaustive Boolean BTOR2 execution for these tiny no-input
// systems. Checks semantics, including the exporter's first-frame monitor.
std::vector<bool> exportedBadFrames(const KInductionProblem& problem, size_t frames) {
  std::ostringstream output;
  exportSecBtor2(problem, output);
  struct Node { size_t id; std::string op; std::vector<std::string> args; };
  std::vector<Node> nodes;
  std::map<size_t, size_t> offsets, initial, next;
  std::istringstream source(output.str());
  std::string line;
  while (std::getline(source, line)) {
    std::istringstream fields(line.substr(0, line.find(';')));
    Node node{};
    if (!(fields >> node.id >> node.op)) continue;
    std::string argument;
    while (fields >> argument) node.args.push_back(argument);
    if (node.op == "state") offsets.emplace(node.id, offsets.size());
    if (node.op == "init" || node.op == "next") {
      (node.op == "init" ? initial : next).emplace(std::stoul(node.args.at(1)),
                                                  std::stoul(node.args.at(2)));
    }
    nodes.push_back(std::move(node));
  }
  if (offsets.size() > 10) throw std::runtime_error("Tiny BTOR2 test state limit");
  std::set<size_t> reachable;
  for (size_t value = 0; value != (size_t{1} << offsets.size()); ++value) reachable.insert(value);
  std::vector<bool> result(frames, false);
  for (size_t frame = 0; frame != frames; ++frame) {
    std::set<size_t> successors;
    for (size_t state : reachable) {
      std::map<size_t, bool> values;
      bool legal = true, bad = false;
      for (const auto& node : nodes) {
        const auto operand = [&](size_t index) { return values.at(std::stoul(node.args.at(index))); };
        if (node.op == "state") values[node.id] = (state >> offsets.at(node.id)) & 1;
        else if (node.op == "const") values[node.id] = node.args.at(1) == "1";
        else if (node.op == "zero") values[node.id] = false;
        else if (node.op == "one" || node.op == "ones") values[node.id] = true;
        else if (node.op == "not") values[node.id] = !operand(1);
        else if (node.op == "and") values[node.id] = operand(1) && operand(2);
        else if (node.op == "or") values[node.id] = operand(1) || operand(2);
        else if (node.op == "xor") values[node.id] = operand(1) != operand(2);
        else if (node.op == "eq" || node.op == "xnor") values[node.id] = operand(1) == operand(2);
        else if (node.op == "ite") values[node.id] = operand(1) ? operand(2) : operand(3);
        else if (node.op == "constraint") legal &= operand(0);
        else if (node.op == "bad") bad |= operand(0);
        else if (node.op != "sort" && node.op != "init" && node.op != "next" && node.op != "output")
          throw std::runtime_error("Unexpected BTOR2 Boolean opcode: " + node.op);
      }
      if (frame == 0) for (const auto& [symbol, expression] : initial)
        legal &= values.at(symbol) == values.at(expression);
      if (!legal) continue;
      result[frame] = result[frame] || bad;
      size_t successor = 0;
      for (const auto& [symbol, expression] : next)
        successor |= size_t(values.at(expression)) << offsets.at(symbol);
      successors.insert(successor);
    }
    EXPECT_FALSE(successors.empty()) << "Exporter removed all behavior at frame " << frame;
    reachable = std::move(successors);
  }
  return result;
}

TEST_F(RelationalInitializationTests, ExportPreservesRelationAlongsideUnitsAndOnlyAtBoot) {
  auto problem = relationalProblem();
  EXPECT_EQ(exportedBadFrames(problem, 3), std::vector<bool>({false, false, false}));
  problem.transitions1 = {{3, BoolExpr::Not(BoolExpr::Var(3))}};
  EXPECT_EQ(exportedBadFrames(problem, 3), std::vector<bool>({false, true, false}));
}

TEST_F(RelationalInitializationTests, ExportKeepsInitialMismatchWithoutObservationMask) {
  auto problem = relationalProblem();
  problem.initialCondition = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  EXPECT_EQ(exportedBadFrames(problem, 2), std::vector<bool>({true, true}));
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
