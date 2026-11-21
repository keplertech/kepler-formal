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

TEST_F(SignalCorrespondenceTests, LearningStopsWhenNothingMoreCanBeProved) {
  testing::internal::CaptureStdout();
  const auto proved = learnSignalCorrespondence(
      shadowBit(/*staleCopy=*/false, /*resetPort=*/false), Config::SolverType::KISSAT,
      /*diagnostics=*/true);
  const std::string provedLog = testing::internal::GetCapturedStdout();
  EXPECT_THAT(proved.provedOutputs, testing::ElementsAre(true));
  EXPECT_THAT(proved.stateEqualities, UnorderedElementsAre(std::pair<size_t, size_t>{2, 3}));
  // Everything was proved by the first run; nothing deeper or larger is tried.
  EXPECT_NE(provedLog.find("scope=registers_and_outputs frames=1"), std::string::npos)
      << provedLog;
  EXPECT_EQ(provedLog.find("frames=2"), std::string::npos) << provedLog;
  EXPECT_EQ(provedLog.find("scope=all_signals"), std::string::npos) << provedLog;

  // The stale copy makes the outputs differ in a state simulation reaches.
  // No depth and no scope can prove such an output, so the learner stops too.
  const auto stale = shadowBit(/*staleCopy=*/true, /*resetPort=*/false);
  const auto run = proveSignalCorrespondence(stale, kRegisters, Config::SolverType::KISSAT);
  EXPECT_THAT(run.differingOutputs, testing::ElementsAre(true));
  testing::internal::CaptureStdout();
  const auto open =
      learnSignalCorrespondence(stale, Config::SolverType::KISSAT, /*diagnostics=*/true);
  const std::string openLog = testing::internal::GetCapturedStdout();
  EXPECT_THAT(open.provedOutputs, testing::ElementsAre(false));
  EXPECT_NE(openLog.find("frames=1"), std::string::npos) << openLog;
  EXPECT_EQ(openLog.find("frames=2"), std::string::npos) << openLog;
  EXPECT_EQ(openLog.find("scope=all_signals"), std::string::npos) << openLog;
}

TEST_F(SignalCorrespondenceTests, GivesUpOnAProblemWithoutUsableCones) {
  auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/false);
  problem.observedOutputExprs1 = {BoolExpr::createInvalid()};
  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  EXPECT_FALSE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(result.stateEqualities, IsEmpty());
}

TEST_F(SignalCorrespondenceTests, CarriesXorAndConstantsThroughTheResetWindow) {
  // Each design toggles a register while `enable` is set and clears it on
  // reset. Design 2 also has `zero`, whose next state is the constant 0 and
  // which is unknown at boot, and `mark`, the exclusive or of reset and
  // `zero`; it outputs its register ORed with `zero`.
  //   2: q (design 1)   3: p, 4: zero, 5: mark (design 2)   6: enable, 7: reset
  KInductionProblem problem;
  problem.allSymbols = {2, 3, 4, 5, 6, 7};
  problem.state0Symbols = {2};
  problem.state1Symbols = {3, 4, 5};
  problem.inputSymbols = {6, 7};
  problem.resetBootstrapCycles = 2;
  problem.resetBootstrapInputs = {{7, true}};
  BoolExpr* idle = BoolExpr::Not(BoolExpr::Var(7));
  const auto toggled = [idle](size_t state) {
    return BoolExpr::And(idle, BoolExpr::Xor(BoolExpr::Var(state), BoolExpr::Var(6)));
  };
  problem.transitions0 = {{2, toggled(2)}};
  problem.transitions1 = {{3, toggled(3)},
                          {4, BoolExpr::createFalse()},
                          {5, BoolExpr::Xor(BoolExpr::Var(7), BoolExpr::Var(4))}};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4))};

  // Nothing has a boot value. The window clears both toggles and settles
  // `zero`, so from the first compared frame on the outputs are equal.
  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(true));
  EXPECT_THAT(result.stateEqualities, IsEmpty());

  // One reset cycle is enough for the toggles and for `zero` as well.
  problem.resetBootstrapCycles = 1;
  EXPECT_THAT(proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT)
                  .provedOutputs,
              testing::ElementsAre(true));
  // Without the window `zero` is unknown in the start state, and so are the
  // toggles.
  problem.resetBootstrapCycles = 0;
  problem.resetBootstrapInputs.clear();
  EXPECT_THAT(proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT)
                  .provedOutputs,
              testing::ElementsAre(false));
}

TEST_F(SignalCorrespondenceTests, UnknownStartRailsAreNotAssumedEqual) {
  // Dual rail, no reset. Design 1's register has no boot value at all;
  // design 2's has only its may-be-zero rail set. 90 and 91 are the rails of
  // a register without a transition, which the pass reads as free inputs.
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.allSymbols = {2, 3, 4, 5, 90, 91};
  problem.state0Symbols = {2, 3, 90, 91};
  problem.state1Symbols = {4, 5};
  problem.dualRailStatePairs = {{2, 3}, {4, 5}, {90, 91}};
  problem.initialStateAssignments = {{5, true}};
  problem.transitions0 = {{2, BoolExpr::Var(2)}, {3, BoolExpr::Var(3)}};
  problem.transitions1 = {{4, BoolExpr::Var(4)}, {5, BoolExpr::Var(5)}};
  problem.dualRailObservedOutputs0 = {{BoolExpr::Var(2), BoolExpr::Var(3)}};
  problem.dualRailObservedOutputs1 = {{BoolExpr::Var(4), BoolExpr::Var(5)}};
  problem.observedOutputExprs0 = {BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4)))};
  problem.observedOutputExprs1 = {BoolExpr::createTrue()};

  for (const auto scope : {kRegisters, kAllSignals}) {
    const auto result = proveSignalCorrespondence(problem, scope, Config::SolverType::KISSAT);
    ASSERT_TRUE(result.completed);
    // Two unknown registers are two different unknowns.
    EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
    EXPECT_THAT(result.stateEqualities, IsEmpty());
  }
}

TEST_F(SignalCorrespondenceTests, CounterexamplesAtTheStartAndInTheStepSplitClasses) {
  // Design 1's register loads the AND of forty inputs, design 2's stays 0.
  // Random simulation never sets all forty, so it groups the two registers,
  // the AND and the constant; only the solver separates them: the AND in the
  // start state, the registers one step later.
  //   2: design 1's register   3: design 2's   10..49: inputs
  KInductionProblem problem;
  problem.allSymbols = {2, 3};
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.initialStateAssignments = {{2, false}, {3, false}};
  BoolExpr* all = BoolExpr::createTrue();
  for (size_t input = 10; input < 50; ++input) {
    problem.allSymbols.push_back(input);
    problem.inputSymbols.push_back(input);
    all = BoolExpr::And(all, BoolExpr::Var(input));
  }
  problem.transitions0 = {{2, all}};
  problem.transitions1 = {{3, BoolExpr::createFalse()}};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(3)};

  const auto result = proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT);
  ASSERT_TRUE(result.completed);
  EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(result.stateEqualities, IsEmpty());
  // The round that split the registers is followed by one that splits nothing.
  EXPECT_GE(result.rounds, 2u);

  // Stopped before that second round, the classes are not a certificate.
  SignalCorrespondenceLimits oneRound;
  oneRound.rounds = 1;
  const auto unfinished =
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT, oneRound);
  EXPECT_FALSE(unfinished.completed);
  EXPECT_THAT(unfinished.stateEqualities, IsEmpty());
}

TEST_F(SignalCorrespondenceTests, GivesUpAtItsLimitsWithoutReportingAnything) {
  const auto problem = shadowBit(/*staleCopy=*/false, /*resetPort=*/false);
  const auto expectNothing = [&problem](SignalCorrespondenceScope scope,
                                        const SignalCorrespondenceLimits& limits) {
    const auto result =
        proveSignalCorrespondence(problem, scope, Config::SolverType::KISSAT, limits);
    EXPECT_FALSE(result.completed);
    EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
    EXPECT_THAT(result.stateEqualities, IsEmpty());
  };
  // The problem has four state symbols.
  SignalCorrespondenceLimits limits;
  limits.stateSymbols = 3;
  expectNothing(kRegisters, limits);
  // Not even the constant fits, then the registers do not.
  limits = {};
  limits.nodes = 0;
  expectNothing(kRegisters, limits);
  limits.nodes = 4;
  expectNothing(kRegisters, limits);
  // Too large for all signals only.
  limits = {};
  limits.allSignalsNodes = 1;
  expectNothing(kAllSignals, limits);
  EXPECT_TRUE(
      proveSignalCorrespondence(problem, kRegisters, Config::SolverType::KISSAT, limits).completed);
  // No solver work allowed.
  limits = {};
  limits.budgetConflicts = 0;
  expectNothing(kRegisters, limits);
  // A cone with an operator the pass does not read.
  auto unreadable = problem;
  BoolExpr notAFormula;
  unreadable.observedOutputExprs1 = {&notAFormula};
  EXPECT_FALSE(
      proveSignalCorrespondence(unreadable, kRegisters, Config::SolverType::KISSAT).completed);
  // Learning reports the run that gave up and proves nothing.
  testing::internal::CaptureStdout();
  const auto learned = learnSignalCorrespondence(unreadable, Config::SolverType::KISSAT,
                                                 /*diagnostics=*/true);
  const std::string log = testing::internal::GetCapturedStdout();
  EXPECT_THAT(learned.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(learned.stateEqualities, IsEmpty());
  EXPECT_NE(log.find("scope=registers_and_outputs frames=1 completed=0"), std::string::npos)
      << log;
  // A depth that gave up is not followed by a deeper one.
  EXPECT_EQ(log.find("frames=2"), std::string::npos) << log;
}

TEST_F(SignalCorrespondenceTests, UndecidedSignalsAreDroppedNotAssumed) {
  // With no conflict allowed a query is decided by propagation or not at all.
  SignalCorrespondenceLimits noSearch;
  noSearch.groupConflicts = 0;
  noSearch.classConflicts = 0;
  noSearch.groupClasses = 4;
  for (const bool staleCopy : {false, true}) {
    const auto problem = shadowBit(staleCopy, /*resetPort=*/false);
    const auto full = proveSignalCorrespondence(problem, kAllSignals, Config::SolverType::KISSAT);
    const auto limited =
        proveSignalCorrespondence(problem, kAllSignals, Config::SolverType::KISSAT, noSearch);
    ASSERT_TRUE(full.completed);
    ASSERT_TRUE(limited.completed);
    // An undecided signal leaves its class, so the relation can only shrink
    // and a design that differs stays unproved. Groups are halved and classes
    // asked member by member on the way, which costs more queries.
    EXPECT_LE(limited.classes, full.classes);
    EXPECT_GT(limited.solverCalls, full.solverCalls);
    if (staleCopy) {
      EXPECT_THAT(limited.provedOutputs, testing::ElementsAre(false));
    }
  }
}

// Two phases per operation. In phase 0 each design outputs its bit and
// latches it; in phase 1 a write stores the complement of the latched bit.
// Design 2 also writes a copy `buf` and reads it in place of the bit while
// `valid` is set. In phase 1 nothing shows whether the copy is right, so the
// classes are not inductive over one time frame; over two, the frame of
// phase 0 is always among the assumed ones. With `staleCopy` the copy is
// never written.
//   2: phase, 3: m, 4: latch (design 1)
//   5: phase, 6: m, 7: latch, 8: buf, 9: valid (design 2)   10: write, 11: clear
KInductionProblem twoPhaseShadowBit(bool staleCopy) {
  KInductionProblem problem;
  problem.allSymbols = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  problem.state0Symbols = {2, 3, 4};
  problem.state1Symbols = {5, 6, 7, 8, 9};
  problem.inputSymbols = {10, 11};
  problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}, {5, false},
                                     {6, false}, {7, false}, {8, true},  {9, false}};
  BoolExpr* write = BoolExpr::Var(10);
  BoolExpr* clear = BoolExpr::Var(11);
  const auto v = [](size_t symbol) { return BoolExpr::Var(symbol); };
  BoolExpr* read1 = mux(v(9), v(8), v(6));
  BoolExpr* writing1 = BoolExpr::And(v(5), write);
  problem.transitions0 = {
      {2, BoolExpr::Not(v(2))},
      {3, mux(BoolExpr::And(v(2), write), BoolExpr::Not(v(4)), v(3))},
      {4, mux(v(2), v(4), v(3))}};
  problem.transitions1 = {
      {5, BoolExpr::Not(v(5))},
      {6, mux(writing1, BoolExpr::Not(v(7)), v(6))},
      {7, mux(v(5), v(7), read1)},
      {8, staleCopy ? v(8) : mux(writing1, BoolExpr::Not(v(7)), v(8))},
      {9, BoolExpr::Or(writing1,
                       BoolExpr::And(BoolExpr::Not(BoolExpr::And(v(5), clear)), v(9)))}};
  problem.observedOutputExprs0 = {BoolExpr::And(BoolExpr::Not(v(2)), v(3))};
  problem.observedOutputExprs1 = {BoolExpr::And(BoolExpr::Not(v(5)), read1)};
  return problem;
}

TEST_F(SignalCorrespondenceTests, RelationOverAMultiPhaseOperationNeedsDeeperInduction) {
  const auto problem = twoPhaseShadowBit(/*staleCopy=*/false);
  const auto solver = Config::SolverType::KISSAT;
  // Registers, outputs and next-state functions say nothing about the copy
  // in phase 1: one frame proves the phase registers only.
  const auto oneFrame = proveSignalCorrespondence(problem, kRegisters, solver, {}, 1);
  ASSERT_TRUE(oneFrame.completed);
  EXPECT_THAT(oneFrame.provedOutputs, testing::ElementsAre(false));
  EXPECT_THAT(oneFrame.stateEqualities, UnorderedElementsAre(std::pair<size_t, size_t>{2, 5}));

  for (const size_t frames : {2, 4}) {
    const auto deeper = proveSignalCorrespondence(problem, kRegisters, solver, {}, frames);
    ASSERT_TRUE(deeper.completed);
    EXPECT_THAT(deeper.provedOutputs, testing::ElementsAre(true));
    EXPECT_THAT(deeper.stateEqualities,
                UnorderedElementsAre(std::pair<size_t, size_t>{2, 5},
                                     std::pair<size_t, size_t>{3, 6},
                                     std::pair<size_t, size_t>{4, 7}));
  }
  // Among all signals is the read of design 2 itself, the copy or the bit,
  // which equals the bit in every phase: one frame is enough there.
  EXPECT_THAT(proveSignalCorrespondence(problem, kAllSignals, solver, {}, 1).provedOutputs,
              testing::ElementsAre(true));

  // One frame leaves the output open without showing a difference, so the
  // learner deepens the induction by itself.
  EXPECT_THAT(oneFrame.differingOutputs, testing::ElementsAre(false));
  testing::internal::CaptureStdout();
  const auto learned = learnSignalCorrespondence(problem, solver, /*diagnostics=*/true);
  const std::string log = testing::internal::GetCapturedStdout();
  EXPECT_THAT(learned.provedOutputs, testing::ElementsAre(true));
  EXPECT_NE(log.find("frames=1 completed=1"), std::string::npos) << log;
  EXPECT_NE(log.find("frames=2 completed=1"), std::string::npos) << log;
  EXPECT_EQ(log.find("frames=4"), std::string::npos) << log;
  EXPECT_EQ(log.find("scope=all_signals"), std::string::npos) << log;

  // Held to one frame, it goes on to all signals, where one frame is enough.
  SignalCorrespondenceLimits shallow;
  shallow.inductionFrames = 1;
  testing::internal::CaptureStdout();
  const auto wide = learnSignalCorrespondence(problem, solver, /*diagnostics=*/true, shallow);
  const std::string wideLog = testing::internal::GetCapturedStdout();
  EXPECT_THAT(wide.provedOutputs, testing::ElementsAre(true));
  EXPECT_NE(wideLog.find("scope=all_signals frames=1 completed=1"), std::string::npos) << wideLog;

  // A copy that is never written is a real difference at every depth.
  const auto stale = twoPhaseShadowBit(/*staleCopy=*/true);
  for (const size_t frames : {1, 2, 4, 8}) {
    const auto result = proveSignalCorrespondence(stale, kAllSignals, solver, {}, frames);
    ASSERT_TRUE(result.completed);
    EXPECT_THAT(result.provedOutputs, testing::ElementsAre(false));
  }
  EXPECT_THAT(learnSignalCorrespondence(stale, solver, /*diagnostics=*/false).provedOutputs,
              testing::ElementsAre(false));

  // The unrolled frames count against the node limit.
  SignalCorrespondenceLimits small;
  small.nodes = 100;
  EXPECT_TRUE(proveSignalCorrespondence(problem, kRegisters, solver, small, 1).completed);
  EXPECT_FALSE(proveSignalCorrespondence(problem, kRegisters, solver, small, 8).completed);
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
