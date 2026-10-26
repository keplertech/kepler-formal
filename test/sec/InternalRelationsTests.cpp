// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cstdlib>
#include <gtest/gtest.h>
#include <sstream>

#include "export/SecBtor2Exporter.h"

#include "BoolExprCache.h"
#include "common/BoolExprUtils.h"
#include "model/SequentialDesignModel.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "imc/IMCEngine.h"
#include "proof/InternalRelations.h"
#include "proof/ProofEngineShared.h"
#include "strategy/SequentialEquivalenceStrategy.h"

// Set when this test binary runs under AddressSanitizer.
#if defined(__SANITIZE_ADDRESS__)
#define KEPLER_TEST_ADDRESS_SANITIZER 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define KEPLER_TEST_ADDRESS_SANITIZER 1
#endif
#endif

namespace KEPLER_FORMAL::SEC {
namespace {

class InternalRelationsTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  KInductionProblem binaryRing() {
    KInductionProblem problem;
    problem.allSymbols = {2, 3, 4, 5};
    problem.state0Symbols = {2, 3};
    problem.state1Symbols = {4, 5};
    problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}, {5, false}};
    problem.transitions0 = {{2, BoolExpr::Var(3)}, {3, BoolExpr::Not(BoolExpr::Var(2))}};
    problem.transitions1 = {{4, BoolExpr::Var(5)}, {5, BoolExpr::Not(BoolExpr::Var(4))}};
    return problem;
  }

  KInductionProblem heldX() {
    auto problem = binaryRing();
    problem.usesDualRailStateEncoding = true;
    problem.dualRailStatePairs = {{2, 3}, {4, 5}};
    problem.initialStateAssignments = {{2, true}, {3, true}, {4, true}, {5, true}};
    problem.transitions0 = {{2, BoolExpr::Var(2)}, {3, BoolExpr::Var(3)}};
    problem.transitions1 = {{4, BoolExpr::Var(4)}, {5, BoolExpr::Var(5)}};
    return problem;
  }

  SequentialDesignModel heldModel(bool initialized, bool value, bool extraKnownOutput = false) {
    const SignalKey state{{1}, {1}};
    const SignalKey output{{2}, {2}};
    SequentialDesignModel model;
    model.stateBits = {state};
    model.inputVarByKey[state] = 2;
    model.displayNameByKey[state] = "q";
    model.nextStateExprByStateKey[state] = BoolExpr::Var(2);
    model.allObservedOutputs = model.observedOutputs = {output};
    model.displayNameByKey[output] = "y";
    model.observedOutputExprByKey[output] = BoolExpr::Var(2);
    if (initialized) {
      model.initialStateValueByKey[state] = value;
    }
    if (extraKnownOutput) {
      const SignalKey known{{3}, {3}};
      model.allObservedOutputs.push_back(known);
      model.observedOutputs.push_back(known);
      model.displayNameByKey[known] = "known";
      model.observedOutputExprByKey[known] = BoolExpr::createTrue();
    }
    return model;
  }
};

TEST_F(InternalRelationsTests, MergingProvedRelationsFoldsIdenticalLogicToOneNode) {
  // Two designs with one register each, proved equal, reading shared input 4.
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.transitions0 = {{2, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(4))}};
  problem.transitions1 = {{3, BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4))}};
  problem.observedOutputExprs0 = {BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(4))};
  problem.observedOutputExprs1 = {BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4))};
  problem.sameFrameStateEqualityPairs0 = {{2, 3}};
  ASSERT_NE(problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);

  mergeProvedStateRelations(problem, /*diagnostics=*/false);

  // Design 2 now reads the representative register, so identical logic is the
  // same hash-consed node and the output comparison needs no search.
  EXPECT_EQ(problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  EXPECT_EQ(problem.transitions0[0].second, problem.transitions1[0].second);
  const auto support = problem.observedOutputExprs1[0]->getSupportVars();
  EXPECT_TRUE(support.contains(2));
  EXPECT_FALSE(support.contains(3));
  // Both registers stay in the state; only reads of the merged one changed.
  EXPECT_EQ(problem.transitions1[0].first, 3u);
}

TEST_F(InternalRelationsTests, JointInductionProvesMutuallySupportingRegisters) {
  const auto problem = binaryRing();
  const InternalRelationCandidate first{{{2, 4}}, {}};
  const InternalRelationCandidate second{{{3, 5}}, {}};
  for (auto solver : {Config::SolverType::KISSAT, Config::SolverType::CADICAL,
                      Config::SolverType::GLUCOSE}) {
    EXPECT_TRUE(proveInternalRelations(problem, {first}, {}, solver).empty());
    EXPECT_EQ(proveInternalRelations(problem, {first, second}, {}, solver).size(), 2u);
    EXPECT_TRUE(proveInternalRelations(problem, {first, second}, {false, true}, solver).empty());
  }
}

TEST_F(InternalRelationsTests, WideChainWithSwappedPinsIsProvedStructurally) {
  // A 32-input AND reduction feeds one register in each design; the second
  // design has the pins swapped, as physical optimization does. The pair is
  // proved by hashing the sorted chain, with no solver call.
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.allSymbols = {2, 3};
  for (size_t input = 4; input < 36; ++input) {
    problem.inputSymbols.push_back(input);
    problem.allSymbols.push_back(input);
  }
  problem.initialStateAssignments = {{2, false}, {3, false}};
  BoolExpr* ascending = BoolExpr::Var(4);
  BoolExpr* descending = BoolExpr::Var(35);
  for (size_t i = 1; i < 32; ++i) {
    ascending = BoolExpr::And(ascending, BoolExpr::Var(4 + i));
    descending = BoolExpr::And(descending, BoolExpr::Var(35 - i));
  }
  problem.transitions0 = {{2, ascending}};
  problem.transitions1 = {{3, descending}};

  setenv("KEPLER_SEC_DIAG", "1", 1);
  testing::internal::CaptureStdout();
  const auto proved = proveInternalRelations(problem, {{{{2, 3}}, {}}}, {},
                                             Config::SolverType::KISSAT);
  const std::string diag = testing::internal::GetCapturedStdout();
  unsetenv("KEPLER_SEC_DIAG");

  EXPECT_EQ(proved.size(), 1u);
  EXPECT_NE(diag.find("round=0 active=1 structural=1 refuted=0"), std::string::npos)
      << diag;
}

TEST_F(InternalRelationsTests, RefutedHypothesesAreRemovedBeforeRecheckingDependents) {
  auto problem = binaryRing();
  problem.transitions1[1].second = BoolExpr::Var(4);
  EXPECT_TRUE(proveInternalRelations(problem, {{{{2, 4}}, {}}, {{{3, 5}}, {}}}, {},
                                    Config::SolverType::KISSAT).empty());
}

TEST_F(InternalRelationsTests, IndependentUninitializedBooleanRegistersAreNotEqual) {
  auto problem = binaryRing();
  problem.initialStateAssignments.clear();
  EXPECT_TRUE(proveInternalRelations(problem, {{{{2, 4}}, {}}, {{{3, 5}}, {}}}, {},
                                    Config::SolverType::KISSAT).empty());
}

TEST_F(InternalRelationsTests, XEqualityUsesCompleteTernaryValuesAndHonorsSwitch) {
  auto problem = heldX();
  const InternalRelationCandidate candidate{{{2, 4}, {3, 5}}, {{2, 3}, {4, 5}}};
  EXPECT_EQ(proveInternalRelations(problem, {candidate}, {}, Config::SolverType::KISSAT).size(), 2u);
  EXPECT_TRUE(proveInternalRelations(problem, {candidate}, {true, false},
                                    Config::SolverType::KISSAT).empty());
  problem.initialStateAssignments[3].second = false;
  EXPECT_TRUE(proveInternalRelations(problem, {candidate}, {}, Config::SolverType::KISSAT).empty());
}

TEST_F(InternalRelationsTests, MatchingXAtBootDoesNotCertifyDivergingTransitions) {
  auto problem = heldX();
  problem.transitions0 = {{2, BoolExpr::createTrue()}, {3, BoolExpr::createFalse()}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}, {5, BoolExpr::createTrue()}};
  EXPECT_TRUE(proveInternalRelations(problem, {{{{2, 4}, {3, 5}}, {{2, 3}, {4, 5}}}}, {},
                                    Config::SolverType::KISSAT).empty());
}

TEST_F(InternalRelationsTests, KnownValuesCanBeLearnedWithXEqualityDisabled) {
  auto problem = heldX();
  problem.initialStateAssignments = {{2, false}, {3, true}, {4, false}, {5, true}};
  EXPECT_EQ(proveInternalRelations(problem, {{{{2, 4}, {3, 5}}, {{2, 3}, {4, 5}}}},
                                  {true, false}, Config::SolverType::KISSAT).size(), 2u);
}

TEST_F(InternalRelationsTests, DisallowingXChecksDefinednessInTheInductionStepToo) {
  auto problem = heldX();
  problem.initialStateAssignments = {{2, false}, {3, true}, {4, false}, {5, true}};
  problem.transitions0 = {{2, BoolExpr::createTrue()}, {3, BoolExpr::createTrue()}};
  problem.transitions1 = {{4, BoolExpr::createTrue()}, {5, BoolExpr::createTrue()}};
  const InternalRelationCandidate candidate{{{2, 4}, {3, 5}}, {{2, 3}, {4, 5}}};
  EXPECT_EQ(proveInternalRelations(problem, {candidate}, {}, Config::SolverType::KISSAT).size(), 2u);
  EXPECT_TRUE(proveInternalRelations(problem, {candidate}, {true, false},
                                    Config::SolverType::KISSAT).empty());
}

TEST_F(InternalRelationsTests, OptionsLeaveTheOutputProofObligationUnchanged) {
  const auto model = heldModel(false, false);
  for (bool learn : {false, true}) {
    for (bool allowX : {false, true}) {
      auto problem = heldX();
      problem.property = BoolExpr::createFalse();
      problem.bad = BoolExpr::createTrue();
      problem.inductionProperty = BoolExpr::Var(2);
      problem.observedOutputExprs0 = {BoolExpr::Var(2)};
      problem.observedOutputExprs1 = {BoolExpr::Var(4)};
      learnInternalStateRelations(model, model, problem, {learn, allowX},
                                  Config::SolverType::KISSAT, false);
      EXPECT_EQ(problem.property, BoolExpr::createFalse());
      EXPECT_EQ(problem.bad, BoolExpr::createTrue());
      EXPECT_EQ(problem.inductionProperty, BoolExpr::Var(2));
      EXPECT_EQ(problem.observedOutputExprs0, std::vector<BoolExpr*>({BoolExpr::Var(2)}));
      EXPECT_EQ(problem.observedOutputExprs1, std::vector<BoolExpr*>({BoolExpr::Var(4)}));
      EXPECT_EQ(problem.sameFrameStateEqualityPairs0.size(), learn && allowX ? 2u : 0u);
    }
  }
}

TEST_F(InternalRelationsTests, DisablingBothSwitchesPreservesThePreparedProblemExactly) {
  const auto model = heldModel(false, false);
  auto problem = heldX();
  problem.property = BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4)));
  problem.bad = BoolExpr::Not(problem.property);
  problem.observedOutputNames = {"y"};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  std::ostringstream before, after;
  exportSecBtor2(problem, before);
  learnInternalStateRelations(model, model, problem, {false, false},
                              Config::SolverType::KISSAT, false);
  exportSecBtor2(problem, after);
  EXPECT_EQ(before.str(), after.str());
  EXPECT_TRUE(problem.sameFrameStateEqualityPairs0.empty());
}

TEST_F(InternalRelationsTests, OutputAssumptionsCannotCertifyAnInvalidInternalRelation) {
  auto problem = heldX();
  problem.property = problem.inductionProperty = BoolExpr::createFalse();
  problem.transitions0 = {{2, BoolExpr::createTrue()}, {3, BoolExpr::createFalse()}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}, {5, BoolExpr::createTrue()}};
  const auto model = heldModel(false, false);
  learnInternalStateRelations(model, model, problem, {}, Config::SolverType::KISSAT, false);
  // The same-name guess proposes equality, which the induction step refutes:
  // one register becomes constant 1, the other constant 0. What survives is
  // the valid complement found by rehashing (rails crosswise), never the
  // equality the false property would have "justified".
  EXPECT_EQ(problem.sameFrameStateEqualityPairs0,
            (std::vector<std::pair<size_t, size_t>>{{2, 5}, {3, 4}}));
}

TEST_F(InternalRelationsTests, SharedDriversSuggestOnlyInitiallyCompatibleSameDesignRelations) {
  auto model = heldModel(true, false);
  const SignalKey second{{5}, {5}};
  model.stateBits.push_back(second);
  model.inputVarByKey[second] = 3;
  model.displayNameByKey[second] = "other_name";
  for (const auto& key : model.stateBits) {
    model.nextStateExprByStateKey[key] = BoolExpr::createFalse();
  }
  auto problem = binaryRing();
  problem.transitions0 = {{2, BoolExpr::createFalse()}, {3, BoolExpr::createFalse()}};
  problem.state1Symbols.clear();
  problem.transitions1.clear();
  problem.allSymbols = {2, 3};
  problem.initialStateAssignments = {{2, false}, {3, false}};
  learnInternalStateRelations(model, {}, problem, {}, Config::SolverType::KISSAT, false);
  EXPECT_EQ(problem.sameFrameStateEqualityPairs0,
            (std::vector<std::pair<size_t, size_t>>{{2, 3}}));
  problem.sameFrameStateEqualityPairs0.clear();
  problem.initialStateAssignments[1].second = true;
  learnInternalStateRelations(model, {}, problem, {}, Config::SolverType::KISSAT, false);
  EXPECT_TRUE(problem.sameFrameStateEqualityPairs0.empty());
}

// Two designs with a two-stage pipeline fed by input symbol 6. Only the first
// stage reads the input directly; the second stage reads the first.
struct TwoStage {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  KInductionProblem problem;
};

TwoStage twoStagePipelines(bool sameNetNames) {
  TwoStage t;
  const SignalKey a0{{1}, {0}}, a1{{1}, {1}}, b0{{2}, {0}}, b1{{2}, {1}};
  t.model0.stateBits = {a0, a1};
  t.model1.stateBits = {b0, b1};
  t.model0.displayNameByKey = {{a0, "a0"}, {a1, "a1"}};
  t.model1.displayNameByKey = {{b0, "b0"}, {b1, "b1"}};
  t.model0.nextStateExprByStateKey = {{a0, BoolExpr::Var(6)}, {a1, BoolExpr::Var(2)}};
  t.model1.nextStateExprByStateKey = {{b0, BoolExpr::Var(6)}, {b1, BoolExpr::Var(4)}};
  if (sameNetNames) {
    t.model0.outputNetNameByKey = {{a0, "n0"}, {a1, "n1"}};
    t.model1.outputNetNameByKey = {{b0, "n0"}, {b1, "n1"}};
  }
  t.problem.allSymbols = {2, 3, 4, 5, 6};
  t.problem.inputSymbols = {6};
  t.problem.state0Symbols = {2, 3};
  t.problem.state1Symbols = {4, 5};
  t.problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}, {5, false}};
  t.problem.transitions0 = {{2, BoolExpr::Var(6)}, {3, BoolExpr::Var(2)}};
  t.problem.transitions1 = {{4, BoolExpr::Var(6)}, {5, BoolExpr::Var(4)}};
  return t;
}

TEST_F(InternalRelationsTests, PropagationReachesRegistersBehindProvedOnes) {
  // No names match. Stage 0 of both designs reads the same input and is found
  // structurally; stage 1 reads different state symbols and is found only
  // after stage 0 is merged and the transitions rehashed (van Eijk, 2000).
  auto t = twoStagePipelines(false);
  learnInternalStateRelations(t.model0, t.model1, t.problem, {}, Config::SolverType::KISSAT, false);
  EXPECT_EQ(t.problem.sameFrameStateEqualityPairs0,
            (std::vector<std::pair<size_t, size_t>>{{2, 4}, {3, 5}}));
}

TEST_F(InternalRelationsTests, OutputNetNamesProposeWhatRehashingCannotSee) {
  // Stage 1 of design 1 computes the same function with a different shape, so
  // it never hashes equal to its counterpart. The shared output net name
  // proposes the pair and the solver proves it.
  for (bool sameNetNames : {false, true}) {
    auto t = twoStagePipelines(sameNetNames);
    // (s & x) | (s & !x) == s, a shape neither the factory nor the rewrite
    // table folds; only the solver sees through it.
    BoolExpr* reshaped = BoolExpr::Or(
        BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(6)),
        BoolExpr::And(BoolExpr::Var(4), BoolExpr::Not(BoolExpr::Var(6))));
    t.model1.nextStateExprByStateKey[t.model1.stateBits[1]] = reshaped;
    t.problem.transitions1[1].second = reshaped;
    learnInternalStateRelations(t.model0, t.model1, t.problem, {}, Config::SolverType::KISSAT, false);
    const std::vector<std::pair<size_t, size_t>> expected =
        sameNetNames ? std::vector<std::pair<size_t, size_t>>{{2, 4}, {3, 5}}
                     : std::vector<std::pair<size_t, size_t>>{{2, 4}};
    EXPECT_EQ(t.problem.sameFrameStateEqualityPairs0, expected) << "sameNetNames=" << sameNetNames;
  }
}

TEST_F(InternalRelationsTests, DualRailComplementsArePairedAcrossRails) {
  // Design 0 exposes Q, design 1 exposes QN of the same register: the rails
  // of one are the rails of the other exchanged. Rehashing pairs them with the
  // rails crosswise, which needs no negation in the dual-rail encoding.
  const SignalKey q{{1}, {0}}, qn{{2}, {0}};
  SequentialDesignModel model0, model1;
  model0.stateBits = {q};
  model1.stateBits = {qn};
  model0.displayNameByKey = {{q, "q"}};
  model1.displayNameByKey = {{qn, "qn"}};
  model0.nextStateExprByStateKey = {{q, BoolExpr::Var(10)}};
  model1.nextStateExprByStateKey = {{qn, BoolExpr::Var(11)}};
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.allSymbols = {2, 3, 4, 5, 10, 11};
  problem.inputSymbols = {10, 11};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.dualRailStatePairs = {{2, 3}, {4, 5}};
  problem.initialStateAssignments = {{2, true}, {3, true}, {4, true}, {5, true}};
  // Q may be one when D may be one; QN may be one when D may be zero.
  problem.transitions0 = {{2, BoolExpr::Var(10)}, {3, BoolExpr::Var(11)}};
  problem.transitions1 = {{4, BoolExpr::Var(11)}, {5, BoolExpr::Var(10)}};
  learnInternalStateRelations(model0, model1, problem, {}, Config::SolverType::KISSAT, false);
  EXPECT_EQ(problem.sameFrameStateEqualityPairs0,
            (std::vector<std::pair<size_t, size_t>>{{2, 5}, {3, 4}}));
}

TEST_F(InternalRelationsTests, SameNamesAndMatchingStartupXCannotHideBinaryMismatch) {
  EXPECT_TRUE(InternalRelationOptions{}.learnInternalRelations);
  EXPECT_TRUE(InternalRelationOptions{}.allowXEqualityInInternalRelations);
  for (bool initialized : {false, true}) {
    auto zero = heldModel(initialized, false);
    zero.nextStateExprByStateKey.begin()->second = BoolExpr::createFalse();
    auto diverging = zero;
    diverging.nextStateExprByStateKey.begin()->second = BoolExpr::createTrue();
    for (auto engine : {SecEngine::Pdr, SecEngine::KInduction, SecEngine::Imc}) {
      for (bool learn : {false, true}) {
        for (bool allowX : {false, true}) {
          SCOPED_TRACE(::testing::Message() << "engine=" << static_cast<int>(engine)
              << " initialized=" << initialized << " learn=" << learn << " allowX=" << allowX);
          SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT, engine);
          strategy.setInternalRelationOptions({learn, allowX});
          EXPECT_EQ(strategy.runExtractedModels(zero, diverging, 2).status,
                    SequentialEquivalenceStatus::Different);
          if (initialized) {
            EXPECT_EQ(strategy.runExtractedModels(zero, zero, 2).status,
                      SequentialEquivalenceStatus::Equivalent);
          }
        }
      }
    }
  }
}

TEST_F(InternalRelationsTests, BinaryRoundSharesTransitionsThatDifferOnlyInShape) {
  // Both designs hold q <= in, but design 2 writes it as a mux of in by a
  // register that is never defined; the per-gate lifting makes that X
  // whenever in is 1, so the rail check refutes the pair.
  const SignalKey q{{1}, {1}};
  const SignalKey y{{2}, {2}};
  const SignalKey sel{{3}, {3}};
  const SignalKey in{{4}, {4}};
  const auto design = [&](bool restructured) {
    SequentialDesignModel model;
    model.stateBits = {q, sel};
    model.inputVarByKey[q] = 2;
    model.inputVarByKey[sel] = 3;
    model.inputVarByKey[in] = 4;
    model.environmentInputs = {in};
    model.displayNameByKey[q] = "q";
    model.displayNameByKey[sel] = "sel";
    model.displayNameByKey[in] = "in";
    model.displayNameByKey[y] = "y";
    model.nextStateExprByStateKey[sel] = BoolExpr::Var(3);
    model.nextStateExprByStateKey[q] = restructured
        ? BoolExpr::Or(BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)),
                       BoolExpr::And(BoolExpr::Not(BoolExpr::Var(3)), BoolExpr::Var(4)))
        : BoolExpr::Var(4);
    model.allObservedOutputs = model.observedOutputs = {y};
    model.observedOutputExprByKey[y] = BoolExpr::Var(2);
    return model;
  };
  const auto model0 = design(false);
  const auto model1 = design(true);
  // The lazily remapped dual-rail problem: rails 10..17 for q and sel of each
  // design, symbol 18 for the shared input.
  const auto problemFor = [&]() {
    KInductionProblem problem;
    problem.usesDualRailStateEncoding = true;
    problem.state0Symbols = {10, 11, 12, 13};
    problem.state1Symbols = {14, 15, 16, 17};
    problem.inputSymbols = {18};
    problem.allSymbols = {10, 11, 12, 13, 14, 15, 16, 17, 18};
    problem.dualRailStatePairs = {{10, 11}, {12, 13}, {14, 15}, {16, 17}};
    for (size_t symbol = 10; symbol < 18; ++symbol) {
      problem.initialStateAssignments.emplace_back(symbol, true);
    }
    problem.observedOutputNames = {"y"};
    problem.observedOutputExprs0 = {BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::Var(14)))};
    problem.observedOutputExprs1 = {BoolExpr::createTrue()};
    problem.property = problem.observedOutputExprs0[0];
    auto store = std::make_shared<LazyTransitionStore>();
    const std::array<const SequentialDesignModel*, 2> models{&model0, &model1};
    for (size_t design = 0; design < 2; ++design) {
      store->localToCombinedByDesign[design] = {{4, 18}};
      const size_t base = 10 + 4 * design;
      store->dualRailStateByLocalSymbolByDesign[design] = {
          {2, {base, base + 1}}, {3, {base + 2, base + 3}}};
      for (size_t bit = 0; bit < 2; ++bit) {
        BoolExpr* next = models[design]->nextStateExprByStateKey.at(bit == 0 ? q : sel);
        store->sourceByStateSymbol[base + 2 * bit] = {design, next, LazyTransitionRail::DualRailOne};
        store->sourceByStateSymbol[base + 2 * bit + 1] = {design, next, LazyTransitionRail::DualRailZero};
      }
    }
    problem.lazyTransitions = store;
    return problem;
  };
  for (bool binaryRound : {false, true}) {
    auto problem = problemFor();
    learnInternalStateRelations(model0, model1, problem, {true, true, binaryRound},
                                Config::SolverType::KISSAT, false);
    // sel is X in both designs and holds; the rail check proves it either way.
    const auto& pairs = problem.sameFrameStateEqualityPairs0;
    EXPECT_TRUE(std::find(pairs.begin(), pairs.end(), std::make_pair(12ul, 16ul)) != pairs.end());
    EXPECT_EQ(std::find(pairs.begin(), pairs.end(), std::make_pair(10ul, 14ul)) != pairs.end(),
              binaryRound) << "binary_relation_round=" << binaryRound;
    EXPECT_EQ(problem.sharedTransitionPairs.size(), binaryRound ? 2u : 0u);
    EXPECT_EQ(problem.observedOutputExprs0[0] == BoolExpr::createTrue(), binaryRound);
    EXPECT_EQ(problem.property == BoolExpr::createTrue(), binaryRound);
  }
}

TEST_F(InternalRelationsTests, BinaryDischargeSkipsOutputsWithoutAFunctionOnBothSides) {
  // Same designs as BinaryRoundSharesTransitionsThatDifferOnlyInShape, plus
  // two observed output names that have no function to compare: "w" is a
  // name design 1 knows without an expression, "z" is unknown to both.
  const SignalKey q{{1}, {1}};
  const SignalKey y{{2}, {2}};
  const SignalKey sel{{3}, {3}};
  const SignalKey in{{4}, {4}};
  const SignalKey w{{5}, {5}};
  const auto design = [&](bool restructured) {
    SequentialDesignModel model;
    model.stateBits = {q, sel};
    model.inputVarByKey[q] = 2;
    model.inputVarByKey[sel] = 3;
    model.inputVarByKey[in] = 4;
    model.environmentInputs = {in};
    model.displayNameByKey[q] = "q";
    model.displayNameByKey[sel] = "sel";
    model.displayNameByKey[in] = "in";
    model.displayNameByKey[y] = "y";
    model.nextStateExprByStateKey[sel] = BoolExpr::Var(3);
    model.nextStateExprByStateKey[q] = restructured
        ? BoolExpr::Or(BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)),
                       BoolExpr::And(BoolExpr::Not(BoolExpr::Var(3)), BoolExpr::Var(4)))
        : BoolExpr::Var(4);
    model.allObservedOutputs = model.observedOutputs = {y};
    model.observedOutputExprByKey[y] = BoolExpr::Var(2);
    if (!restructured) {
      model.allObservedOutputs.push_back(w);
      model.displayNameByKey[w] = "w";
    }
    return model;
  };
  const auto model0 = design(false);
  const auto model1 = design(true);
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {10, 11, 12, 13};
  problem.state1Symbols = {14, 15, 16, 17};
  problem.inputSymbols = {18};
  problem.allSymbols = {10, 11, 12, 13, 14, 15, 16, 17, 18};
  problem.dualRailStatePairs = {{10, 11}, {12, 13}, {14, 15}, {16, 17}};
  for (size_t symbol = 10; symbol < 18; ++symbol) {
    problem.initialStateAssignments.emplace_back(symbol, true);
  }
  problem.observedOutputNames = {"y", "w", "z"};
  BoolExpr* yEquality = BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::Var(14)));
  BoolExpr* wEquality = BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(12), BoolExpr::Var(16)));
  problem.observedOutputExprs0 = {yEquality, wEquality, wEquality};
  problem.observedOutputExprs1 = {BoolExpr::createTrue(), BoolExpr::createTrue(),
                                  BoolExpr::createTrue()};
  problem.property = BoolExpr::And(yEquality, wEquality);
  auto store = std::make_shared<LazyTransitionStore>();
  const std::array<const SequentialDesignModel*, 2> models{&model0, &model1};
  for (size_t designIndex = 0; designIndex < 2; ++designIndex) {
    store->localToCombinedByDesign[designIndex] = {{4, 18}};
    const size_t base = 10 + 4 * designIndex;
    store->dualRailStateByLocalSymbolByDesign[designIndex] = {
        {2, {base, base + 1}}, {3, {base + 2, base + 3}}};
    for (size_t bit = 0; bit < 2; ++bit) {
      BoolExpr* next = models[designIndex]->nextStateExprByStateKey.at(bit == 0 ? q : sel);
      store->sourceByStateSymbol[base + 2 * bit] = {designIndex, next, LazyTransitionRail::DualRailOne};
      store->sourceByStateSymbol[base + 2 * bit + 1] = {designIndex, next, LazyTransitionRail::DualRailZero};
    }
  }
  problem.lazyTransitions = store;

  learnInternalStateRelations(model0, model1, problem, {true, true, true},
                              Config::SolverType::KISSAT, false);

  // "y" is one function on both sides under the proved merge; the other two
  // outputs keep their obligations.
  EXPECT_EQ(problem.observedOutputExprs0[0], BoolExpr::createTrue());
  EXPECT_EQ(problem.observedOutputExprs0[1], wEquality);
  EXPECT_EQ(problem.observedOutputExprs0[2], wEquality);
}

TEST_F(InternalRelationsTests, BinaryRetryWithoutLazyTransitionsProvesNothing) {
  const auto problem = binaryRing();
  EXPECT_EQ(problem.lazyTransitions, nullptr);
  EXPECT_TRUE(proveInternalRelations(problem, {{{{2, 4}}, {}}}, {},
                                    Config::SolverType::KISSAT,
                                    /*commonStartBinary=*/true)
                  .empty());
}

// The most significant bit of the product of two `a.size()`-bit operands,
// as an array multiplier of AND, XOR and OR nodes.
BoolExpr* productMsb(const std::vector<BoolExpr*>& a,
                     const std::vector<BoolExpr*>& b) {
  const size_t n = a.size();
  std::vector<BoolExpr*> acc(2 * n, BoolExpr::createFalse());
  for (size_t i = 0; i < n; ++i) {
    BoolExpr* carry = BoolExpr::createFalse();
    for (size_t j = 0; j < n; ++j) {
      BoolExpr* partial = BoolExpr::And(a[j], b[i]);
      BoolExpr* half = BoolExpr::Xor(acc[i + j], partial);
      BoolExpr* sum = BoolExpr::Xor(half, carry);
      carry = BoolExpr::Or(BoolExpr::And(acc[i + j], partial),
                           BoolExpr::And(carry, half));
      acc[i + j] = sum;
    }
    for (size_t k = i + n; k < 2 * n; ++k) {
      BoolExpr* sum = BoolExpr::Xor(acc[k], carry);
      carry = BoolExpr::And(acc[k], carry);
      acc[k] = sum;
    }
  }
  return acc[2 * n - 1];
}

TEST_F(InternalRelationsTests, UndecidedPartitionGivesUpEachPairUnderItsOwnBudget) {
#if defined(KEPLER_TEST_ADDRESS_SANITIZER)
  // The hard instance drives CaDiCaL 3.0.0 into its transitive-reduction
  // inprocessing, whose memcpy AddressSanitizer reports as reading past a
  // buffer (seen with libc++ on macOS). That is the solver's own code, so
  // the budget path is exercised in the unsanitized builds only.
  GTEST_SKIP() << "CaDiCaL transitive reduction trips AddressSanitizer";
#endif
  // Two registers holding the top bit of a*b and of b*a: equal, but the
  // commutativity of a 14-bit multiplier is far beyond the partition's
  // conflict budget, so the pair is retried alone and then given up.
  constexpr size_t kWidth = 14;
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.allSymbols = {2, 3};
  std::vector<BoolExpr*> a;
  std::vector<BoolExpr*> b;
  for (size_t bit = 0; bit < kWidth; ++bit) {
    for (auto* operand : {&a, &b}) {
      const size_t symbol = problem.allSymbols.size() + 10;
      problem.inputSymbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
      operand->push_back(BoolExpr::Var(symbol));
    }
  }
  problem.initialStateAssignments = {{2, false}, {3, false}};
  problem.transitions0 = {{2, productMsb(a, b)}};
  problem.transitions1 = {{3, productMsb(b, a)}};

  setenv("KEPLER_SEC_DIAG", "1", 1);
  testing::internal::CaptureStdout();
  const auto proved = proveInternalRelations(problem, {{{{2, 3}}, {}}}, {},
                                             Config::SolverType::KISSAT);
  const std::string diag = testing::internal::GetCapturedStdout();
  unsetenv("KEPLER_SEC_DIAG");

  EXPECT_TRUE(proved.empty());
  EXPECT_NE(diag.find("status=unknown"), std::string::npos) << diag;
  EXPECT_NE(diag.find("SEC diag: internal relation given up"), std::string::npos)
      << diag;
}

TEST_F(InternalRelationsTests, RefutedRelationsAreReportedWithDiagnostics) {
  auto problem = binaryRing();
  problem.transitions1[1].second = BoolExpr::Var(4);

  setenv("KEPLER_SEC_DIAG", "1", 1);
  testing::internal::CaptureStdout();
  const auto proved = proveInternalRelations(
      problem, {{{{2, 4}}, {}}, {{{3, 5}}, {}}}, {}, Config::SolverType::KISSAT);
  const std::string diag = testing::internal::GetCapturedStdout();
  unsetenv("KEPLER_SEC_DIAG");

  EXPECT_TRUE(proved.empty());
  EXPECT_NE(diag.find("SEC diag: internal relation refuted:"), std::string::npos)
      << diag;
}

TEST_F(InternalRelationsTests, HashingGuessSkipsRegistersWithoutATransition) {
  // The second register of the model has no transition in the problem, as
  // an opaque state has: the rehashing guess leaves it out without proposing
  // anything for it.
  auto model = heldModel(true, false);
  const SignalKey second{{5}, {5}};
  model.stateBits.push_back(second);
  model.inputVarByKey[second] = 3;
  model.displayNameByKey[second] = "other_name";
  model.nextStateExprByStateKey[second] = BoolExpr::Var(3);
  auto problem = binaryRing();
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  problem.state1Symbols.clear();
  problem.transitions1.clear();
  problem.allSymbols = {2, 3};
  problem.initialStateAssignments = {{2, false}, {3, false}};

  learnInternalStateRelations(model, {}, problem, {}, Config::SolverType::KISSAT, false);

  EXPECT_TRUE(problem.sameFrameStateEqualityPairs0.empty());
}

TEST_F(InternalRelationsTests, MergeSkipsVeryLargeLazilyRemappedProblems) {
  KInductionProblem problem;
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  problem.state0Symbols.assign((size_t{1} << 17) + 1, 2);
  problem.sameFrameStateEqualityPairs0 = {{2, 3}};

  testing::internal::CaptureStdout();
  mergeProvedStateRelations(problem, /*diagnostics=*/true);
  const std::string diag = testing::internal::GetCapturedStdout();

  EXPECT_NE(diag.find("SEC diag: merged state symbols skipped: 131073 state symbols"),
            std::string::npos)
      << diag;
  EXPECT_EQ(problem.sameFrameStateEqualityPairs0.size(), 1u);
}

TEST_F(InternalRelationsTests, MergeWithOnlySelfPairsLeavesTheProblemUnchanged) {
  auto problem = binaryRing();
  problem.sameFrameStateEqualityPairs0 = {{2, 2}};
  problem.sameFrameStateEqualityPairs1 = {{4, 4}};
  const auto transitions0 = problem.transitions0;
  const auto transitions1 = problem.transitions1;

  mergeProvedStateRelations(problem, /*diagnostics=*/false);

  EXPECT_EQ(problem.transitions0, transitions0);
  EXPECT_EQ(problem.transitions1, transitions1);
}

TEST_F(InternalRelationsTests, ResetProofBehaviorIsPreservedWithLearning) {
  auto model = heldModel(false, false);
  const SignalKey reset{{4}, {4}};
  model.environmentInputs = {reset};
  model.inputVarByKey[reset] = 3;
  model.displayNameByKey[reset] = "rst";
  model.nextStateExprByStateKey.begin()->second = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(3)), BoolExpr::Var(2));
  for (auto engine : {SecEngine::Pdr, SecEngine::KInduction, SecEngine::Imc}) {
    SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
        engine, SecEncoding::DualRailSteady, {1, {{"rst", true}}});
    strategy.setInternalRelationOptions({false, false});
    const auto baseline = strategy.runExtractedModels(model, model, 2);
    for (bool allowX : {false, true}) {
      strategy.setInternalRelationOptions({true, allowX});
      const auto result = strategy.runExtractedModels(model, model, 2);
      EXPECT_EQ(result.status, baseline.status);
      EXPECT_EQ(result.coveredOutputs, baseline.coveredOutputs);
    }
  }
}

TEST_F(InternalRelationsTests, UnfinishedRefinementPublishesNoSurvivingGuesses) {
  KInductionProblem problem;
  std::vector<InternalRelationCandidate> candidates;
  for (size_t i = 0; i < 70; ++i) {
    const size_t lhs = 2 + 2 * i;
    const size_t rhs = lhs + 1;
    problem.allSymbols.insert(problem.allSymbols.end(), {lhs, rhs});
    problem.state0Symbols.push_back(lhs);
    problem.state1Symbols.push_back(rhs);
    problem.initialStateAssignments.emplace_back(lhs, false);
    problem.initialStateAssignments.emplace_back(rhs, false);
    problem.transitions0.emplace_back(lhs,
        i == 69 ? BoolExpr::createTrue() : BoolExpr::Var(lhs + 2));
    problem.transitions1.emplace_back(rhs,
        i == 69 ? BoolExpr::createFalse() : BoolExpr::Var(rhs + 2));
    candidates.push_back({{{lhs, rhs}}, {}});
  }
  EXPECT_TRUE(proveInternalRelations(problem, candidates, {}, Config::SolverType::KISSAT).empty());
}

SequentialDesignModel internalRelationRingModel(const std::string& prefix = "") {
  SequentialDesignModel model;
  const SignalKey q{{1}, {1}}, r{{2}, {2}}, output{{3}, {3}};
  model.stateBits = {q, r};
  model.inputVarByKey = {{q, 2}, {r, 3}};
  model.displayNameByKey = {{q, prefix + "q"}, {r, prefix + "r"}, {output, "y"}};
  model.initialStateValueByKey = {{q, false}, {r, false}};
  model.nextStateExprByStateKey = {{q, BoolExpr::Var(3)}, {r, BoolExpr::Not(BoolExpr::Var(2))}};
  model.allObservedOutputs = model.observedOutputs = {output};
  model.observedOutputExprByKey[output] = BoolExpr::Var(2);
  return model;
}

TEST_F(InternalRelationsTests, ExactImcUsesCertifiedRelationsToCloseTheInductionGap) {
  const auto model = internalRelationRingModel();
  for (auto solver : {Config::SolverType::KISSAT, Config::SolverType::CADICAL,
                      Config::SolverType::GLUCOSE}) {
    for (bool learn : {false, true}) {
      for (bool allowX : {false, true}) {
        SCOPED_TRACE(::testing::Message() << "solver=" << static_cast<int>(solver)
                     << " learn=" << learn << " allowX=" << allowX);
        SequentialEquivalenceStrategy strategy(nullptr, nullptr, solver,
                                               SecEngine::Imc, SecEncoding::Binary);
        strategy.setInternalRelationOptions({learn, allowX});
        // The one-step frontier is not closed, and output equality alone is
        // not inductive. The two jointly certified register equalities close it.
        const auto result = strategy.runExtractedModels(model, model, 0);
        EXPECT_EQ(result.status, learn ? SequentialEquivalenceStatus::Equivalent
                                       : SequentialEquivalenceStatus::Inconclusive);
      }
    }
  }
}

TEST_F(InternalRelationsTests, ExactImcConstraintsRespectBothSwitchesAndPreserveTheProperty) {
  const auto model = heldModel(false, false);
  for (bool learn : {false, true}) {
    for (bool allowX : {false, true}) {
      auto problem = heldX();
      problem.auxiliaryStateSymbols = {6};
      problem.allSymbols.push_back(6);
      problem.initialStateAssignments.emplace_back(6, false);
      problem.auxiliaryTransitions = {{6, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4))}};
      problem.property = BoolExpr::Not(BoolExpr::Var(6));
      problem.bad = BoolExpr::Var(6);
      learnInternalStateRelations(model, model, problem, {learn, allowX},
                                  Config::SolverType::KISSAT, false);
      const auto constraints = problem.learnedInternalRelationInvariant;
      const bool enabled = learn && allowX;
      EXPECT_EQ(constraints != nullptr, enabled);
      EXPECT_EQ(isInductiveInvariant(problem, problem.property, Config::SolverType::KISSAT,
                                     constraints), enabled);
      // Existing callers that omit the optional constraints keep the old query.
      EXPECT_FALSE(isInductiveInvariant(problem, problem.property, Config::SolverType::KISSAT));
      EXPECT_FALSE(invariantExcludesBadStates(problem, BoolExpr::createTrue(),
                                             Config::SolverType::KISSAT, constraints));
      EXPECT_EQ(problem.property, BoolExpr::Not(BoolExpr::Var(6)));
      EXPECT_EQ(problem.bad, BoolExpr::Var(6));

      auto relationOutput = problem;
      relationOutput.bad = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4));
      EXPECT_EQ(invariantExcludesBadStates(relationOutput, BoolExpr::createTrue(),
                                           Config::SolverType::KISSAT, constraints), enabled);
    }
  }
}

TEST_F(InternalRelationsTests, ExactInterpolationConstrainsBothTransitionFramesOnlyWhenEnabled) {
  const auto model = internalRelationRingModel();
  const std::unordered_map<size_t, size_t> nextSymbols{{2, 6}, {3, 7}, {4, 8}, {5, 9}};
  for (bool learn : {false, true}) {
    auto problem = binaryRing();
    learnInternalStateRelations(model, model, problem, {learn, false},
                                Config::SolverType::KISSAT, false);
    const auto original = buildOneStepTransitionFormula(problem, nextSymbols);
    const auto constrained = buildOneStepTransitionFormula(
        problem, nextSymbols, problem.learnedInternalRelationInvariant);
    if (!learn) {
      EXPECT_EQ(original, constrained);
    }
    for (const auto& [lhs, rhs] : {std::pair<size_t, size_t>{2, 4}, {6, 8}}) {
      const auto unequal = BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs));
      EXPECT_TRUE(isProofFormulaSatisfiable(BoolExpr::And(original, unequal), Config::SolverType::KISSAT));
      EXPECT_EQ(isProofFormulaSatisfiable(BoolExpr::And(constrained, unequal),
                                         Config::SolverType::KISSAT), !learn);
    }
  }
}

TEST_F(InternalRelationsTests, ExactImcReachableFrontierStillClosesWithLearningEnabledOrDisabled) {
  auto model0 = internalRelationRingModel("left_");
  auto model1 = internalRelationRingModel("right_");
  const SignalKey helper{{4}, {4}};
  for (auto* model : {&model0, &model1}) {
    model->stateBits.push_back(helper);
    model->inputVarByKey[helper] = 4;
    model->displayNameByKey[helper] = "helper";
    model->initialStateValueByKey[helper] = false;
    model->nextStateExprByStateKey[helper] = BoolExpr::createFalse();
  }
  for (bool learn : {false, true}) {
    SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
                                           SecEngine::Imc, SecEncoding::Binary);
    strategy.setInternalRelationOptions({learn, false});
    // Only the unrelated helper matches by name. The public ring equality
    // still needs the exact reachable frontier to cover its four-state cycle.
    const auto result = strategy.runExtractedModels(model0, model1, 3);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_EQ(result.bound, 3u);
  }
}

TEST_F(InternalRelationsTests, ExactImcRetainsCounterexamplesWithValidLearnedRelations) {
  const auto model0 = internalRelationRingModel();
  auto model1 = model0;
  model1.observedOutputExprByKey.begin()->second = BoolExpr::Not(BoolExpr::Var(2));
  for (bool learn : {false, true}) {
    for (bool allowX : {false, true}) {
      SequentialEquivalenceStrategy strategy(nullptr, nullptr, Config::SolverType::KISSAT,
                                             SecEngine::Imc, SecEncoding::Binary);
      strategy.setInternalRelationOptions({learn, allowX});
      EXPECT_EQ(strategy.runExtractedModels(model0, model1, 0).status,
                SequentialEquivalenceStatus::Different);
    }
  }
}

TEST_F(InternalRelationsTests, RejectedCandidatesNeverReachExactImcConstraints) {
  auto problem = binaryRing();
  problem.transitions1[1].second = BoolExpr::Var(4);
  const auto model = internalRelationRingModel();
  learnInternalStateRelations(model, model, problem, {}, Config::SolverType::KISSAT, false);
  EXPECT_EQ(problem.learnedInternalRelationInvariant, nullptr);
  EXPECT_TRUE(problem.sameFrameStateEqualityPairs0.empty());
}

// Two designs of one output register q over held registers. Design 0 uses the
// symbols 2 (q) and 3..3+width-1, design 1 the same count after them.
struct GatePair {
  KInductionProblem problem;
  std::vector<InternalRelationCandidate> candidates;
  std::vector<BoolExpr*> in0;
  std::vector<BoolExpr*> in1;
};

GatePair gatePair(size_t width) {
  GatePair pair;
  const size_t stride = width + 1;
  for (size_t i = 0; i < stride; ++i) {
    const size_t lhs = 2 + i;
    const size_t rhs = lhs + stride;
    pair.problem.state0Symbols.push_back(lhs);
    pair.problem.state1Symbols.push_back(rhs);
    pair.problem.initialStateAssignments.emplace_back(lhs, false);
    pair.problem.initialStateAssignments.emplace_back(rhs, false);
    pair.candidates.push_back({{{lhs, rhs}}, {}});
    if (i != 0) {
      pair.problem.transitions0.emplace_back(lhs, BoolExpr::Var(lhs));
      pair.problem.transitions1.emplace_back(rhs, BoolExpr::Var(rhs));
      pair.in0.push_back(BoolExpr::Var(lhs));
      pair.in1.push_back(BoolExpr::Var(rhs));
    }
  }
  pair.problem.allSymbols = pair.problem.state0Symbols;
  pair.problem.allSymbols.insert(pair.problem.allSymbols.end(),
                                 pair.problem.state1Symbols.begin(),
                                 pair.problem.state1Symbols.end());
  return pair;
}

size_t provedPairs(GatePair pair, BoolExpr* next0, BoolExpr* next1) {
  pair.problem.transitions0.emplace_back(2, next0);
  pair.problem.transitions1.emplace_back(2 + pair.in0.size() + 1, next1);
  return proveInternalRelations(pair.problem, pair.candidates, {},
                                Config::SolverType::KISSAT).size();
}

TEST_F(InternalRelationsTests, SwappedPinsOfASymmetricGateAreTheSameTransition) {
  const auto pair = gatePair(3);
  const auto& a = pair.in0;
  const auto& b = pair.in1;
  EXPECT_EQ(provedPairs(pair, BoolExpr::And(BoolExpr::And(a[0], a[1]), a[2]),
                        BoolExpr::And(BoolExpr::And(b[2], b[1]), b[0])), 4u);
  EXPECT_EQ(provedPairs(pair, BoolExpr::Not(BoolExpr::Or(BoolExpr::Or(a[0], a[1]), a[2])),
                        BoolExpr::Not(BoolExpr::Or(b[1], BoolExpr::Or(b[2], b[0])))), 4u);
  EXPECT_EQ(provedPairs(pair, BoolExpr::Xor(a[0], BoolExpr::Not(a[1])),
                        BoolExpr::Xor(BoolExpr::Not(b[1]), b[0])), 4u);
}

TEST_F(InternalRelationsTests, SwappedDataPinsOfAMuxAreADifferentTransition) {
  const auto pair = gatePair(3);
  const auto mux = [](BoolExpr* select, BoolExpr* one, BoolExpr* zero) {
    return BoolExpr::Or(BoolExpr::And(select, one),
                        BoolExpr::And(BoolExpr::Not(select), zero));
  };
  const auto& a = pair.in0;
  const auto& b = pair.in1;
  EXPECT_EQ(provedPairs(pair, mux(a[0], a[1], a[2]), mux(b[0], b[1], b[2])), 4u);
  // The three held registers stay proved; only q is refuted.
  EXPECT_EQ(provedPairs(pair, mux(a[0], a[1], a[2]), mux(b[0], b[2], b[1])), 3u);
}

TEST_F(InternalRelationsTests, EqualFunctionsOfDifferentShapeAreLeftToTheSolver) {
  const auto pair = gatePair(2);
  const auto& a = pair.in0;
  const auto& b = pair.in1;
  EXPECT_EQ(provedPairs(pair, BoolExpr::Not(BoolExpr::And(a[0], a[1])),
                        BoolExpr::Or(BoolExpr::Not(b[0]), BoolExpr::Not(b[1]))), 3u);
  EXPECT_EQ(provedPairs(pair, BoolExpr::And(a[0], a[1]), BoolExpr::Or(b[0], b[1])), 2u);
}

TEST_F(InternalRelationsTests, LongChainsKeepTheirShapeAndStillCompare) {
  const auto pair = gatePair(20);
  BoolExpr* same0 = pair.in0[0];
  BoolExpr* same1 = pair.in1[0];
  BoolExpr* reversed1 = pair.in1[19];
  for (size_t i = 1; i < 20; ++i) {
    same0 = BoolExpr::And(same0, pair.in0[i]);
    same1 = BoolExpr::And(same1, pair.in1[i]);
    reversed1 = BoolExpr::And(reversed1, pair.in1[19 - i]);
  }
  EXPECT_EQ(provedPairs(pair, same0, same1), 21u);
  EXPECT_EQ(provedPairs(pair, same0, reversed1), 21u);
}

// Van Eijk's merge step: the same cone over each design's registers becomes
// one formula once the proved register pairs are merged, and the output is
// discharged without an engine.  A cone that differs keeps its obligation.
TEST_F(InternalRelationsTests, StructuralDischargeMergesProvedRegistersInBinary) {
  KInductionProblem problem;
  problem.allSymbols = {2, 3, 4, 5, 6};
  problem.inputSymbols = {6};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.sameFrameStateEqualityPairs0 = {{2, 4}, {3, 5}};
  problem.observedOutputExprs0 = {
      BoolExpr::Xor(BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(6)), BoolExpr::Var(3)),
      BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3))};
  problem.observedOutputExprs1 = {
      BoolExpr::Xor(BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(6)), BoolExpr::Var(5)),
      BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(5))};
  BoolExpr* property = BoolExpr::createTrue();
  for (size_t i = 0; i < 2; ++i) {
    property = BoolExpr::And(property,
        makeEqualityExpr(problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]));
  }
  problem.property = problem.inductionProperty = property;
  problem.bad = problem.inductionBad = BoolExpr::Not(property);

  EXPECT_EQ(dischargeStructurallyEquivalentOutputs(problem, false), 1u);

  ASSERT_EQ(problem.structurallyProvedOutputs.size(), 2u);
  EXPECT_TRUE(problem.structurallyProvedOutputs[0]);
  EXPECT_FALSE(problem.structurallyProvedOutputs[1]);
  EXPECT_EQ(problem.observedOutputExprs0[0], BoolExpr::createTrue());
  EXPECT_EQ(problem.observedOutputExprs1[0], BoolExpr::createTrue());
  EXPECT_NE(problem.observedOutputExprs0[1], BoolExpr::createTrue());
  // The rebuilt property is the remaining obligation only: it fails exactly
  // where the second output pair differs.
  EXPECT_NE(problem.property, BoolExpr::createTrue());
  EXPECT_TRUE(problem.property->evaluate(
      {{2, true}, {3, true}, {4, true}, {5, true}, {6, false}}));
  EXPECT_FALSE(problem.property->evaluate(
      {{2, true}, {3, false}, {4, true}, {5, false}, {6, false}}));
  EXPECT_EQ(problem.inductionProperty, problem.property);
  EXPECT_EQ(problem.bad, problem.inductionBad);
}

TEST_F(InternalRelationsTests, StructuralDischargeRequiresTheProvedEqualities) {
  KInductionProblem problem;
  problem.allSymbols = {2, 3, 4, 5};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.observedOutputExprs0 = {BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3))};
  problem.observedOutputExprs1 = {BoolExpr::Xor(BoolExpr::Var(4), BoolExpr::Var(5))};
  BoolExpr* const original = problem.observedOutputExprs0[0];
  problem.property = problem.inductionProperty =
      makeEqualityExpr(problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = problem.inductionBad = BoolExpr::Not(problem.property);
  BoolExpr* const originalProperty = problem.property;

  // Same shape, different registers, nothing proved: no merge, no discharge.
  EXPECT_EQ(dischargeStructurallyEquivalentOutputs(problem, false), 0u);
  EXPECT_EQ(problem.observedOutputExprs0[0], original);
  EXPECT_EQ(problem.property, originalProperty);
  ASSERT_EQ(problem.structurallyProvedOutputs.size(), 1u);
  EXPECT_FALSE(problem.structurallyProvedOutputs[0]);

  // Chained equalities merge through their representative.
  problem.sameFrameStateEqualityPairs0 = {{2, 4}};
  problem.sameFrameStateEqualityPairs1 = {{5, 3}};
  EXPECT_EQ(dischargeStructurallyEquivalentOutputs(problem, false), 1u);
  EXPECT_EQ(problem.property, BoolExpr::createTrue());
}

TEST_F(InternalRelationsTests, StructuralDischargeComparesBothRailsInDualRail) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.allSymbols = {2, 3, 4, 5};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.sameFrameStateEqualityPairs0 = {{2, 4}, {3, 5}};
  // Output 0: both rails identical after the merge.  Output 1: the may-be-zero
  // rail differs, so the ternary values are not the same formula.
  problem.dualRailObservedOutputs0 = {{BoolExpr::Var(2), BoolExpr::Var(3)},
                                      {BoolExpr::Var(2), BoolExpr::Var(3)}};
  problem.dualRailObservedOutputs1 = {{BoolExpr::Var(4), BoolExpr::Var(5)},
                                      {BoolExpr::Var(4), BoolExpr::Not(BoolExpr::Var(5))}};
  problem.observedOutputExprs0 = {makeEqualityExpr(BoolExpr::Var(2), BoolExpr::Var(4)),
                                  makeEqualityExpr(BoolExpr::Var(3), BoolExpr::Var(5))};
  problem.observedOutputExprs1 = {BoolExpr::createTrue(), BoolExpr::createTrue()};
  problem.dualRailOutputStrictEqualityExprs = problem.observedOutputExprs0;
  BoolExpr* const guard1 = problem.observedOutputExprs0[1];
  problem.property = problem.inductionProperty =
      BoolExpr::And(problem.observedOutputExprs0[0], problem.observedOutputExprs0[1]);
  problem.bad = problem.inductionBad = BoolExpr::Not(problem.property);

  EXPECT_EQ(dischargeStructurallyEquivalentOutputs(problem, false), 1u);

  EXPECT_TRUE(problem.structurallyProvedOutputs[0]);
  EXPECT_FALSE(problem.structurallyProvedOutputs[1]);
  EXPECT_EQ(problem.observedOutputExprs0[0], BoolExpr::createTrue());
  EXPECT_EQ(problem.dualRailOutputStrictEqualityExprs[0], BoolExpr::createTrue());
  EXPECT_EQ(problem.observedOutputExprs0[1], guard1);
  EXPECT_NE(problem.dualRailOutputStrictEqualityExprs[1], BoolExpr::createTrue());
  // Only the second guard remains in the property.
  EXPECT_TRUE(problem.property->evaluate({{2, true}, {3, true}, {4, false}, {5, true}}));
  EXPECT_FALSE(problem.property->evaluate({{2, true}, {3, true}, {4, true}, {5, false}}));
}

// End to end: once the learner has certified the register pairs of two
// identical rings, their identical output cones are discharged and nothing is
// left for an engine.
TEST_F(InternalRelationsTests, LearningDischargesIdenticalOutputsAfterProvingRegisters) {
  auto problem = binaryRing();
  problem.observedOutputNames = {"y"};
  problem.observedOutputExprs0 = {BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3))};
  problem.observedOutputExprs1 = {BoolExpr::Xor(BoolExpr::Var(4), BoolExpr::Var(5))};
  problem.property = problem.inductionProperty =
      makeEqualityExpr(problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = problem.inductionBad = BoolExpr::Not(problem.property);

  const SignalKey first{{1}, {1}};
  const SignalKey second{{2}, {2}};
  auto makeModel = [&](size_t firstSymbol, size_t secondSymbol) {
    SequentialDesignModel model;
    model.stateBits = {first, second};
    model.displayNameByKey[first] = "r0";
    model.displayNameByKey[second] = "r1";
    model.inputVarByKey[first] = firstSymbol;
    model.inputVarByKey[second] = secondSymbol;
    model.nextStateExprByStateKey[first] = BoolExpr::Var(secondSymbol);
    model.nextStateExprByStateKey[second] = BoolExpr::Not(BoolExpr::Var(firstSymbol));
    return model;
  };
  const auto model0 = makeModel(2, 3);
  const auto model1 = makeModel(4, 5);

  learnInternalStateRelations(model0, model1, problem, {}, Config::SolverType::KISSAT, false);

  EXPECT_EQ(problem.sameFrameStateEqualityPairs0.size(), 2u);
  ASSERT_EQ(problem.structurallyProvedOutputs.size(), 1u);
  EXPECT_TRUE(problem.structurallyProvedOutputs[0]);
  EXPECT_EQ(problem.property, BoolExpr::createTrue());
  EXPECT_EQ(problem.bad, BoolExpr::createFalse());
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
