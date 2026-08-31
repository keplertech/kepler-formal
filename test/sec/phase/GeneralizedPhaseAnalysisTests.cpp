// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "phase/GeneralizedPhaseAnalysis.h"

using namespace KEPLER_FORMAL::SEC;
using namespace KEPLER_FORMAL::SEC::PHASE;

namespace {

SignalKey makeSignalKey(size_t ordinal) {
  SignalKey key;
  key.first.push_back(ordinal);
  return key;
}

NormalizedTransitionSystem makePeriodTwoSystem() {
  NormalizedTransitionSystem system;
  const size_t state =
      system.addState(makeSignalKey(1), "period2", TernaryValue::Zero);
  const ExprID current =
      system.variableExpression(system.states()[state].variable);
  system.setNextState(state, system.expressions().makeNot(current));
  return system;
}

NormalizedTransitionSystem makePeriodTwoAndThreeSystem() {
  NormalizedTransitionSystem system;
  const size_t periodTwo =
      system.addState(makeSignalKey(1), "period2", TernaryValue::Zero);
  const size_t periodThreeLow =
      system.addState(makeSignalKey(2), "period3_low", TernaryValue::Zero);
  const size_t periodThreeHigh =
      system.addState(makeSignalKey(3), "period3_high", TernaryValue::Zero);

  const ExprID two =
      system.variableExpression(system.states()[periodTwo].variable);
  const ExprID low =
      system.variableExpression(system.states()[periodThreeLow].variable);
  const ExprID high =
      system.variableExpression(system.states()[periodThreeHigh].variable);
  system.setNextState(periodTwo, system.expressions().makeNot(two));
  // 00 -> 01 -> 10 -> 00 (low is the least-significant state bit).
  system.setNextState(
      periodThreeLow,
      system.expressions().makeNot(system.expressions().makeOr(low, high)));
  system.setNextState(periodThreeHigh, low);
  return system;
}

NormalizedTransitionSystem makeInputDependentSystem() {
  NormalizedTransitionSystem system;
  const size_t input = system.addInput(makeSignalKey(1), "input");
  const size_t state =
      system.addState(makeSignalKey(2), "state", TernaryValue::Zero);
  system.setNextState(state, system.variableExpression(input));
  return system;
}

NormalizedTransitionSystem makeTransientConstantSystem(
    TernaryValue initialValue) {
  NormalizedTransitionSystem system;
  const size_t state =
      system.addState(makeSignalKey(1), "transient", initialValue);
  system.setNextState(state, system.expressions().constant(true));
  return system;
}

struct GatedLatchSystem {
  NormalizedTransitionSystem system;
  size_t gateVariable = 0;
};

GatedLatchSystem makeGatedLatchSystem() {
  GatedLatchSystem result;
  auto& system = result.system;
  result.gateVariable = system.addInput(makeSignalKey(1), "gate");
  const size_t dataVariable = system.addInput(makeSignalKey(2), "data");
  const size_t carrier =
      system.addState(makeSignalKey(3), "carrier", TernaryValue::Zero);
  const size_t latchState =
      system.addState(makeSignalKey(4), "latch_q", TernaryValue::Unknown);

  const ExprID gate = system.variableExpression(result.gateVariable);
  const ExprID data = system.variableExpression(dataVariable);
  const ExprID carrierValue =
      system.variableExpression(system.states()[carrier].variable);
  const ExprID currentLatch =
      system.variableExpression(system.states()[latchState].variable);
  const ExprID transparency = system.expressions().makeAnd(carrierValue, gate);

  system.setNextState(carrier, system.expressions().makeNot(carrierValue));
  system.setNextState(latchState, system.expressions().makeIte(
                                      transparency, data, currentLatch));
  system.addLatch({latchState, "latch0", "synthetic_latch", data, transparency,
                   std::nullopt, std::nullopt, ClearPresetValue::Unknown});
  return result;
}

NormalizedTransitionSystem makeManyGatedLatches(size_t latchCount) {
  NormalizedTransitionSystem system;
  const size_t carrier =
      system.addState(makeSignalKey(1), "carrier", TernaryValue::Zero);
  const ExprID carrierValue =
      system.variableExpression(system.states()[carrier].variable);
  system.setNextState(carrier, system.expressions().makeNot(carrierValue));

  for (size_t index = 0; index < latchCount; ++index) {
    const size_t gateVariable = system.addInput(
        makeSignalKey(2 + index * 2), "gate_" + std::to_string(index));
    const size_t latchState = system.addState(
        makeSignalKey(3 + index * 2), "latch_q_" + std::to_string(index),
        TernaryValue::Unknown);
    const ExprID gate = system.variableExpression(gateVariable);
    const ExprID current =
        system.variableExpression(system.states()[latchState].variable);
    const ExprID transparency =
        system.expressions().makeAnd(carrierValue, gate);
    const ExprID data = system.expressions().constant(false);
    system.setNextState(
        latchState, system.expressions().makeIte(transparency, data, current));
    system.addLatch({latchState, "latch_" + std::to_string(index),
                     "synthetic_latch", data, transparency, std::nullopt,
                     std::nullopt, ClearPresetValue::Unknown});
  }
  return system;
}

GeneralizedPhaseAnalysisOptions analysisOptions(size_t concurrency) {
  GeneralizedPhaseAnalysisOptions options;
  options.maxSimulationSteps = 64;
  options.maxStoredStateBytes = 1024 * 1024;
  options.maxPhases = 8;
  options.maxConcurrency = concurrency;
  return options;
}

void expectSamePeriodicResult(const PeriodicSignalAnalysisResult& lhs,
                              const PeriodicSignalAnalysisResult& rhs) {
  EXPECT_EQ(lhs.phaseCount, rhs.phaseCount);
  EXPECT_EQ(lhs.complete, rhs.complete);
  EXPECT_EQ(lhs.diagnostic, rhs.diagnostic);
  EXPECT_EQ(lhs.statistics.candidateStateCount,
            rhs.statistics.candidateStateCount);
  EXPECT_EQ(lhs.statistics.resetPeriodicStateCount,
            rhs.statistics.resetPeriodicStateCount);
  EXPECT_EQ(lhs.statistics.usableGeneratorCount,
            rhs.statistics.usableGeneratorCount);
  ASSERT_EQ(lhs.generators.size(), rhs.generators.size());
  for (size_t index = 0; index < lhs.generators.size(); ++index) {
    const auto& lhsGenerator = lhs.generators[index];
    const auto& rhsGenerator = rhs.generators[index];
    EXPECT_EQ(lhsGenerator.stateIndex, rhsGenerator.stateIndex);
    EXPECT_EQ(lhsGenerator.minimumPeriod, rhsGenerator.minimumPeriod);
    EXPECT_EQ(lhsGenerator.generatorWord, rhsGenerator.generatorWord);
    EXPECT_EQ(lhsGenerator.expandedWord, rhsGenerator.expandedWord);
    EXPECT_EQ(lhsGenerator.usableAtSelectedPhaseCount,
              rhsGenerator.usableAtSelectedPhaseCount);
  }
}

void expectSameProfileResult(const LatchPhaseProfileResult& lhs,
                             const LatchPhaseProfileResult& rhs) {
  EXPECT_EQ(lhs.complete, rhs.complete);
  EXPECT_EQ(lhs.diagnostic, rhs.diagnostic);
  ASSERT_EQ(lhs.profiles.size(), rhs.profiles.size());
  for (size_t profileIndex = 0; profileIndex < lhs.profiles.size();
       ++profileIndex) {
    const auto& lhsProfile = lhs.profiles[profileIndex];
    const auto& rhsProfile = rhs.profiles[profileIndex];
    EXPECT_EQ(lhsProfile.latchIndex, rhsProfile.latchIndex);
    ASSERT_EQ(lhsProfile.phases.size(), rhsProfile.phases.size());
    for (size_t phase = 0; phase < lhsProfile.phases.size(); ++phase) {
      EXPECT_EQ(lhsProfile.phases[phase].kind, rhsProfile.phases[phase].kind);
      EXPECT_EQ(lhsProfile.phases[phase].residualExpression,
                rhsProfile.phases[phase].residualExpression);
    }
  }
  ASSERT_EQ(lhs.residualExpressions.nodes().size(),
            rhs.residualExpressions.nodes().size());
  for (size_t index = 0; index < lhs.residualExpressions.nodes().size();
       ++index) {
    const auto& lhsNode = lhs.residualExpressions.nodes()[index];
    const auto& rhsNode = rhs.residualExpressions.nodes()[index];
    EXPECT_EQ(lhsNode.op, rhsNode.op);
    EXPECT_EQ(lhsNode.variable, rhsNode.variable);
    EXPECT_EQ(lhsNode.constant, rhsNode.constant);
    EXPECT_EQ(lhsNode.left, rhsNode.left);
    EXPECT_EQ(lhsNode.right, rhsNode.right);
  }
}

void expectSameAnalysis(const GeneralizedPhaseAnalysisResult& lhs,
                        const GeneralizedPhaseAnalysisResult& rhs) {
  EXPECT_EQ(lhs.complete, rhs.complete);
  EXPECT_EQ(lhs.diagnostics, rhs.diagnostics);
  EXPECT_EQ(lhs.simulation.status, rhs.simulation.status);
  EXPECT_EQ(lhs.simulation.trace, rhs.simulation.trace);
  EXPECT_EQ(lhs.simulation.stemLength, rhs.simulation.stemLength);
  EXPECT_EQ(lhs.simulation.cycleLength, rhs.simulation.cycleLength);
  EXPECT_EQ(lhs.simulation.diagnostic, rhs.simulation.diagnostic);
  EXPECT_EQ(lhs.simulation.statistics.simulatedSteps,
            rhs.simulation.statistics.simulatedSteps);
  EXPECT_EQ(lhs.simulation.statistics.storedStates,
            rhs.simulation.statistics.storedStates);
  EXPECT_EQ(lhs.simulation.statistics.stateWidth,
            rhs.simulation.statistics.stateWidth);
  EXPECT_EQ(lhs.simulation.statistics.expressionNodes,
            rhs.simulation.statistics.expressionNodes);
  expectSamePeriodicResult(lhs.periodicSignals, rhs.periodicSignals);
  expectSameProfileResult(lhs.latchProfiles, rhs.latchProfiles);
}

TEST(GeneralizedPhaseAnalysisTests, DiscoversPeriodTwoGenerator) {
  const auto system = makePeriodTwoSystem();
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.simulation.stemLength, 0u);
  EXPECT_EQ(result.simulation.cycleLength, 2u);
  EXPECT_EQ(result.periodicSignals.phaseCount, 2u);
  ASSERT_EQ(result.periodicSignals.generators.size(), 1u);
  const auto& generator = result.periodicSignals.generators.front();
  EXPECT_EQ(generator.stateIndex, 0u);
  EXPECT_EQ(generator.minimumPeriod, 2u);
  EXPECT_EQ(generator.generatorWord,
            (std::vector<TernaryValue>{TernaryValue::Zero, TernaryValue::One}));
  EXPECT_EQ(generator.expandedWord, generator.generatorWord);
}

TEST(GeneralizedPhaseAnalysisTests,
     SelectsSixPhasesForPeriodTwoAndPeriodThreeGenerators) {
  const auto system = makePeriodTwoAndThreeSystem();
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.simulation.cycleLength, 6u);
  EXPECT_EQ(result.periodicSignals.phaseCount, 6u);
  ASSERT_EQ(result.periodicSignals.generators.size(), 3u);
  EXPECT_EQ(result.periodicSignals.generators[0].minimumPeriod, 2u);
  EXPECT_EQ(result.periodicSignals.generators[1].minimumPeriod, 3u);
  EXPECT_EQ(result.periodicSignals.generators[2].minimumPeriod, 3u);
  EXPECT_EQ(result.periodicSignals.statistics.usableGeneratorCount, 3u);
}

TEST(GeneralizedPhaseAnalysisTests,
     DoesNotInventGeneratorForInputDependentUnknownState) {
  const auto system = makeInputDependentSystem();
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.simulation.stemLength, 1u);
  EXPECT_EQ(result.simulation.cycleLength, 1u);
  EXPECT_EQ(result.periodicSignals.phaseCount, 1u);
  EXPECT_TRUE(result.periodicSignals.generators.empty());
  EXPECT_EQ(result.periodicSignals.statistics.resetPeriodicStateCount, 0u);
}

TEST(GeneralizedPhaseAnalysisTests,
     RejectsEventuallyConstantSignalThatIsNotPeriodicFromReset) {
  const auto system = makeTransientConstantSystem(TernaryValue::Zero);
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.simulation.stemLength, 1u);
  EXPECT_EQ(result.simulation.cycleLength, 1u);
  EXPECT_TRUE(result.periodicSignals.generators.empty());
}

TEST(GeneralizedPhaseAnalysisTests,
     RejectsGeneratorCandidateThatIsUnknownInTheStem) {
  const auto system = makeTransientConstantSystem(TernaryValue::Unknown);
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.simulation.stemLength, 1u);
  EXPECT_EQ(result.simulation.cycleLength, 1u);
  EXPECT_TRUE(result.periodicSignals.generators.empty());
}

TEST(GeneralizedPhaseAnalysisTests, PreservesResidualGateInLatchPhaseProfile) {
  const auto fixture = makeGatedLatchSystem();
  const auto result =
      analyzeGeneralizedPhases(fixture.system, analysisOptions(1));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.periodicSignals.phaseCount, 2u);
  ASSERT_EQ(result.latchProfiles.profiles.size(), 1u);
  const auto& profile = result.latchProfiles.profiles.front();
  ASSERT_EQ(profile.phases.size(), 2u);
  EXPECT_EQ(profile.phases[0].kind, TransparencyKind::Closed);
  EXPECT_EQ(profile.phases[1].kind, TransparencyKind::Conditional);
  ASSERT_TRUE(result.latchProfiles.residualExpressions.isValid(
      profile.phases[1].residualExpression));
  const auto& residual = result.latchProfiles.residualExpressions.node(
      profile.phases[1].residualExpression);
  EXPECT_EQ(residual.op, ExpressionOp::Variable);
  EXPECT_EQ(residual.variable, fixture.gateVariable);
}

TEST(GeneralizedPhaseAnalysisTests,
     ParallelAndSerialAnalysisAreDeterministicallyIdentical) {
  const auto fixture = makeGatedLatchSystem();
  const auto serial =
      analyzeGeneralizedPhases(fixture.system, analysisOptions(1));
  ASSERT_TRUE(serial.complete);

  for (size_t iteration = 0; iteration < 20; ++iteration) {
    const auto parallel =
        analyzeGeneralizedPhases(fixture.system, analysisOptions(4));
    expectSameAnalysis(serial, parallel);
  }
}

TEST(GeneralizedPhaseAnalysisTests, EvaluatesStrongKleeneBooleanOperators) {
  NormalizedTransitionSystem system;
  const size_t lhsVariable = system.addInput(makeSignalKey(1), "lhs");
  const size_t rhsVariable = system.addInput(makeSignalKey(2), "rhs");
  const ExprID lhs = system.variableExpression(lhsVariable);
  const ExprID rhs = system.variableExpression(rhsVariable);
  const ExprID notLhs = system.expressions().makeNot(lhs);
  const ExprID conjunction = system.expressions().makeAnd(lhs, rhs);
  const ExprID disjunction = system.expressions().makeOr(lhs, rhs);
  const ExprID exclusiveOr = system.expressions().makeXor(lhs, rhs);

  constexpr std::array values{TernaryValue::Zero, TernaryValue::One,
                              TernaryValue::Unknown};
  constexpr std::array negated{TernaryValue::One, TernaryValue::Zero,
                               TernaryValue::Unknown};
  constexpr TernaryValue andExpected[3][3] = {
      {TernaryValue::Zero, TernaryValue::Zero, TernaryValue::Zero},
      {TernaryValue::Zero, TernaryValue::One, TernaryValue::Unknown},
      {TernaryValue::Zero, TernaryValue::Unknown, TernaryValue::Unknown}};
  constexpr TernaryValue orExpected[3][3] = {
      {TernaryValue::Zero, TernaryValue::One, TernaryValue::Unknown},
      {TernaryValue::One, TernaryValue::One, TernaryValue::One},
      {TernaryValue::Unknown, TernaryValue::One, TernaryValue::Unknown}};
  constexpr TernaryValue xorExpected[3][3] = {
      {TernaryValue::Zero, TernaryValue::One, TernaryValue::Unknown},
      {TernaryValue::One, TernaryValue::Zero, TernaryValue::Unknown},
      {TernaryValue::Unknown, TernaryValue::Unknown, TernaryValue::Unknown}};

  std::vector<TernaryValue> environment(system.variables().size(),
                                        TernaryValue::Unknown);
  for (size_t lhsIndex = 0; lhsIndex < values.size(); ++lhsIndex) {
    environment[lhsVariable] = values[lhsIndex];
    EXPECT_EQ(evaluateTernaryExpression(system, notLhs, environment),
              negated[lhsIndex]);
    for (size_t rhsIndex = 0; rhsIndex < values.size(); ++rhsIndex) {
      SCOPED_TRACE("lhs=" + std::to_string(lhsIndex) +
                   " rhs=" + std::to_string(rhsIndex));
      environment[rhsVariable] = values[rhsIndex];
      EXPECT_EQ(evaluateTernaryExpression(system, conjunction, environment),
                andExpected[lhsIndex][rhsIndex]);
      EXPECT_EQ(evaluateTernaryExpression(system, disjunction, environment),
                orExpected[lhsIndex][rhsIndex]);
      EXPECT_EQ(evaluateTernaryExpression(system, exclusiveOr, environment),
                xorExpected[lhsIndex][rhsIndex]);
    }
  }
}

TEST(GeneralizedPhaseAnalysisTests,
     ReportsDeterministicSimulationBudgetStatuses) {
  const auto system = makePeriodTwoSystem();
  TernarySimulationOptions options;
  options.maxConcurrency = 1;
  options.maxSteps = 1;
  options.maxStoredStateBytes = 0;

  const auto stepLimited = simulateTernary(system, options);
  EXPECT_EQ(stepLimited.status, TernarySimulationStatus::StepLimitExceeded);
  EXPECT_EQ(stepLimited.statistics.simulatedSteps, 1u);
  EXPECT_EQ(stepLimited.statistics.storedStates, 2u);

  options.maxSteps = 4;
  options.maxStoredStateBytes = 1;
  const auto memoryLimited = simulateTernary(system, options);
  EXPECT_EQ(memoryLimited.status, TernarySimulationStatus::MemoryLimitExceeded);
  EXPECT_EQ(memoryLimited.statistics.simulatedSteps, 1u);
  EXPECT_EQ(memoryLimited.statistics.storedStates, 1u);
}

TEST(GeneralizedPhaseAnalysisTests,
     ProfilesTwoThousandLatchesSharingOneCarrier) {
  constexpr size_t kLatchCount = 2000;
  const auto system = makeManyGatedLatches(kLatchCount);
  const auto result = analyzeGeneralizedPhases(system, analysisOptions(4));

  ASSERT_TRUE(result.complete);
  EXPECT_EQ(result.periodicSignals.phaseCount, 2u);
  EXPECT_EQ(result.periodicSignals.statistics.usableGeneratorCount, 1u);
  ASSERT_EQ(result.latchProfiles.profiles.size(), kLatchCount);
  for (const auto& profile : result.latchProfiles.profiles) {
    ASSERT_EQ(profile.phases.size(), 2u);
    EXPECT_EQ(profile.phases[0].kind, TransparencyKind::Closed);
    EXPECT_EQ(profile.phases[1].kind, TransparencyKind::Conditional);
  }
  // Constants plus one residual variable per gate keeps the result linear in
  // the latch count; the analysis never constructs pairwise latch relations.
  EXPECT_LE(result.latchProfiles.residualExpressions.nodes().size(),
            kLatchCount + 2);
}

}  // namespace
