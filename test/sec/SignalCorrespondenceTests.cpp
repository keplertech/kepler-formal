// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "BoolExprCache.h"
#include "proof/InternalRelations.h"
#include "proof/SignalCorrespondence.h"

namespace KEPLER_FORMAL::SEC {
namespace {

using ::testing::IsEmpty;
using ::testing::UnorderedElementsAre;

constexpr auto kRegisters = SignalCorrespondenceScope::RegistersAndOutputs;
constexpr auto kAllSignals = SignalCorrespondenceScope::AllSignals;

BoolExpr* mux(BoolExpr* select, BoolExpr* ifTrue, BoolExpr* ifFalse) {
  return BoolExpr::Or(BoolExpr::And(select, ifTrue),
                      BoolExpr::And(BoolExpr::Not(select), ifFalse));
}

class SignalCorrespondenceTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  // Design 1 keeps one bit `m` and outputs it. Design 2 keeps the same bit
  // and, beside it, a copy `buf` that it outputs while `valid` is set. A write
  // stores the complement of the output, so the registers are equal only
  // while the outputs are, and the outputs only while the registers are: the
  // shape of issue 317. `buf` boots unlike `m`, so no register of design 2 is
  // always equal to it. With `staleCopy` the copy is never written and the
  // designs differ after the first write.
  //   2: m (design 1)   3: m, 4: buf, 5: valid (design 2)   6: write, 7: clear
  KInductionProblem shadowBit(bool staleCopy, bool resetPort) {
    KInductionProblem problem;
    problem.allSymbols = {2, 3, 4, 5, 6, 7};
    problem.state0Symbols = {2};
    problem.state1Symbols = {3, 4, 5};
    problem.inputSymbols = {6, 7};
    BoolExpr* write = BoolExpr::Var(6);
    BoolExpr* clear = BoolExpr::Var(7);
    BoolExpr* out0 = BoolExpr::Var(2);
    BoolExpr* out1 = mux(BoolExpr::Var(5), BoolExpr::Var(4), BoolExpr::Var(3));
    BoolExpr* next0 = mux(write, BoolExpr::Not(out0), BoolExpr::Var(2));
    BoolExpr* next1 = mux(write, BoolExpr::Not(out1), BoolExpr::Var(3));
    BoolExpr* nextValid =
        BoolExpr::Or(write, BoolExpr::And(BoolExpr::Not(clear), BoolExpr::Var(5)));
    if (resetPort) {
      // 8: reset. It clears `m` and `valid` and leaves `buf` alone; nothing
      // has a boot value.
      problem.allSymbols.push_back(8);
      problem.inputSymbols.push_back(8);
      BoolExpr* idle = BoolExpr::Not(BoolExpr::Var(8));
      next0 = BoolExpr::And(idle, next0);
      next1 = BoolExpr::And(idle, next1);
      nextValid = BoolExpr::And(idle, nextValid);
      problem.resetBootstrapCycles = 1;
      problem.resetBootstrapInputs = {{8, true}};
    } else {
      problem.initialStateAssignments = {{2, false}, {3, false}, {4, true}, {5, false}};
    }
    problem.transitions0 = {{2, next0}};
    problem.transitions1 = {
        {3, next1},
        {4, staleCopy ? BoolExpr::Var(4) : mux(write, BoolExpr::Not(out1), BoolExpr::Var(4))},
        {5, nextValid}};
    problem.observedOutputExprs0 = {out0};
    problem.observedOutputExprs1 = {out1};
    return problem;
  }

  // One register per design in the dual-rail encoding, reset to 0 and held.
  //   2, 3: rails of design 1's register   4, 5: of design 2's   6: reset
  KInductionProblem dualRailResetRegisters() {
    KInductionProblem problem;
    problem.usesDualRailStateEncoding = true;
    problem.allSymbols = {2, 3, 4, 5, 6};
    problem.state0Symbols = {2, 3};
    problem.state1Symbols = {4, 5};
    problem.inputSymbols = {6};
    problem.dualRailStatePairs = {{2, 3}, {4, 5}};
    // Both registers boot unknown: both rails set.
    problem.initialStateAssignments = {{2, true}, {3, true}, {4, true}, {5, true}};
    problem.resetBootstrapCycles = 2;
    problem.resetBootstrapInputs = {{6, true}};
    BoolExpr* reset = BoolExpr::Var(6);
    const auto clearedOne = [reset](size_t rail) {
      return BoolExpr::And(BoolExpr::Not(reset), BoolExpr::Var(rail));
    };
    const auto clearedZero = [reset](size_t rail) {
      return BoolExpr::Or(reset, BoolExpr::Var(rail));
    };
    problem.transitions0 = {{2, clearedOne(2)}, {3, clearedZero(3)}};
    problem.transitions1 = {{4, clearedOne(4)}, {5, clearedZero(5)}};
    problem.dualRailObservedOutputs0 = {{BoolExpr::Var(2), BoolExpr::Var(3)}};
    problem.dualRailObservedOutputs1 = {{BoolExpr::Var(4), BoolExpr::Var(5)}};
    // The engines' obligation; the pass reads the rail cones above.
    problem.observedOutputExprs0 = {BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4)))};
    problem.observedOutputExprs1 = {BoolExpr::createTrue()};
    return problem;
  }
};

TEST_F(SignalCorrespondenceTests, ProvesRegistersAndOutputsThatOnlyHoldTogether) {
  for (auto solver : {Config::SolverType::KISSAT, Config::SolverType::CADICAL,
                      Config::SolverType::GLUCOSE}) {
    const auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/false);
    // The register learner has nothing to prove here: `m` of design 2 is not
    // inductive against `m` of design 1 without the outputs, and `buf` is
    // equal to neither at boot.
    const InternalRelationCandidate sameName{{{2, 3}}, {}};
    const InternalRelationCandidate copy{{{3, 4}}, {}};
    EXPECT_THAT(proveInternalRelations(problem, {sameName, copy}, {}, solver), IsEmpty());

    for (const auto scope : {kRegisters, kAllSignals}) {
      const auto result = proveSignalCorrespondence(problem, scope, solver);
      ASSERT_TRUE(result.completed);
      EXPECT_THAT(result.provedOutputs, testing::ElementsAre(true));
      EXPECT_THAT(result.stateEqualities, UnorderedElementsAre(std::pair<size_t, size_t>{2, 3}));
      EXPECT_GE(result.rounds, 1u);
    }
    EXPECT_GT(proveSignalCorrespondence(problem, kAllSignals, solver).signals,
              proveSignalCorrespondence(problem, kRegisters, solver).signals);
  }
}

TEST_F(SignalCorrespondenceTests, KeepsDifferentDesignsUnproved) {
  const auto problem = shadowBit(/*staleCopy=*/true, /*resetPort=*/false);
  for (const auto scope : {kRegisters, kAllSignals}) {
    const auto result = proveSignalCorrespondence(problem, scope, Config::SolverType::KISSAT);
    ASSERT_TRUE(result.completed);
    EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
    EXPECT_THAT(result.stateEqualities, IsEmpty());
  }
}

TEST_F(SignalCorrespondenceTests, StartsFromTheStateTheResetWindowLeaves) {
  auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/true);
  const auto afterReset =
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(afterReset.completed);
  // Equal from the first compared frame on, with `buf` still unknown there.
  EXPECT_THAT(afterReset.provedOutputs, testing::ElementsAre(true));
  // Not equal in every frame: the two `m` are unknown at boot, so the pair is
  // not handed to the engines as an invariant of the reset window.
  EXPECT_THAT(afterReset.stateEqualities, IsEmpty());

  // Without the reset window the start state is the boot state, where the
  // two `m` are unrelated.
  problem.resetBootstrapCycles = 0;
  problem.resetBootstrapInputs.clear();
  const auto fromBoot =
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(fromBoot.completed);
  EXPECT_THAT(fromBoot.provedOutputs, testing::ElementsAre(false));
}

TEST_F(SignalCorrespondenceTests, RelatesSignalsThatHoldOppositeValues) {
  // Each design toggles one register; design 2 keeps the complement and
  // outputs its negation.
  KInductionProblem problem;
  problem.allSymbols = {2, 3};
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.initialStateAssignments = {{2, false}, {3, true}};
  problem.transitions0 = {{2, BoolExpr::Not(BoolExpr::Var(2))}};
  problem.transitions1 = {{3, BoolExpr::Not(BoolExpr::Var(3))}};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Not(BoolExpr::Var(3))};

  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(true));
  // The registers are complements, not equal: no equality is reported.
  EXPECT_THAT(result.stateEqualities, IsEmpty());

  // The same registers with equal boot values are not complements.
  problem.initialStateAssignments = {{2, false}, {3, false}};
  const auto equalBoot =
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(equalBoot.completed);
  EXPECT_THAT(equalBoot.provedOutputs, testing::ElementsAre(false));
}

TEST_F(SignalCorrespondenceTests, DualRailReportsCompleteRegistersEqualInEveryFrame) {
  auto problem = dualRailResetRegisters();
  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(true));
  // Unknown at boot on both sides, 0 on both sides from the first reset cycle.
  EXPECT_THAT(result.stateEqualities,
              UnorderedElementsAre(std::pair<size_t, size_t>{2, 4},
                                   std::pair<size_t, size_t>{3, 5}));

  // Design 2's register boots at 1: equal after the reset window, so the
  // output is proved, but not at boot, so neither rail is reported.
  problem.initialStateAssignments = {{2, true}, {3, true}, {4, true}, {5, false}};
  const auto differentBoot =
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(differentBoot.completed);
  EXPECT_THAT(differentBoot.provedOutputs, testing::ElementsAre(true));
  EXPECT_THAT(differentBoot.stateEqualities, IsEmpty());

  // A pair the register learner certified is not reported again.
  problem = dualRailResetRegisters();
  problem.sameFrameStateEqualityPairs0 = {{2, 4}, {3, 5}};
  EXPECT_THAT(proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT)
                  .stateEqualities,
              IsEmpty());
}

TEST_F(SignalCorrespondenceTests, LeavesOutputsProvedEarlierAlone) {
  auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/false);
  problem.structurallyProvedOutputs = {true};
  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));

  const auto learned =
      learnSignalCorrespondence(problem, Config::SolverType::KISSAT, /*diagnostics=*/false);
  EXPECT_THAT(learned.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(learned.stateEqualities, IsEmpty());
}

TEST_F(SignalCorrespondenceTests, LearningTriesAllSignalsWhenOutputsRemain) {
  testing::internal::CaptureStdout();
  const auto proved = learnSignalCorrespondence(
      shadowBit(/*staleCopy=*/false, /*resetPort=*/false), Config::SolverType::KISSAT,
      /*diagnostics=*/true);
  const std::string provedLog = testing::internal::GetCapturedStdout();
  EXPECT_THAT(proved.provedOutputs, testing::ElementsAre(true));
  EXPECT_THAT(proved.stateEqualities, UnorderedElementsAre(std::pair<size_t, size_t>{2, 3}));
  EXPECT_NE(provedLog.find("scope=registers_and_outputs"), std::string::npos) << provedLog;
  // Everything was proved on registers and outputs; the second scope is not run.
  EXPECT_EQ(provedLog.find("scope=all_signals"), std::string::npos) << provedLog;

  testing::internal::CaptureStdout();
  const auto open = learnSignalCorrespondence(
      shadowBit(/*staleCopy=*/true, /*resetPort=*/false), Config::SolverType::KISSAT,
      /*diagnostics=*/true);
  const std::string openLog = testing::internal::GetCapturedStdout();
  EXPECT_THAT(open.provedOutputs, testing::ElementsAre(false));
  EXPECT_NE(openLog.find("scope=all_signals"), std::string::npos) << openLog;
}

TEST_F(SignalCorrespondenceTests, GivesUpOnAProblemWithoutUsableCones) {
  auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/false);
  problem.observedOutputExprs1 = {BoolExpr::createInvalid()};
  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  EXPECT_FALSE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(result.stateEqualities, IsEmpty());
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
