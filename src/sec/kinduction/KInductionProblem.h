// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BoolExpr.h"
#include "common/SignalKey.h"
#include "proof/DualRailEncoding.h"

namespace KEPLER_FORMAL::SEC {

enum class LazyTransitionRail {
  Binary,
  DualRailOne,
  DualRailZero,
};

struct DualRailSymbolPair {
  size_t mayBeOne = 0;
  size_t mayBeZero = 0;
};

struct LazyTransitionSource {
  size_t designIndex = 0;
  BoolExpr* localExpr = nullptr;
  LazyTransitionRail rail = LazyTransitionRail::Binary;
};

struct InductionTransitionSupportCache;

struct DualRailResidualPublicKiAttempt {
  BoolExpr* outputExpr0 = nullptr;
  BoolExpr* outputExpr1 = nullptr;
  BoolExpr* inductionProperty = nullptr;
  BoolExpr* inductionBad = nullptr;
  size_t inductionOutputCount = 0;
  size_t solverType = 0;
  size_t attemptedMaxK = 0;
  size_t resultBound = 0;

  bool matches(BoolExpr* candidateOutputExpr0,
               BoolExpr* candidateOutputExpr1,
               BoolExpr* candidateInductionProperty,
               BoolExpr* candidateInductionBad,
               size_t candidateInductionOutputCount,
               size_t candidateSolverType) const {
    return outputExpr0 == candidateOutputExpr0 &&
           outputExpr1 == candidateOutputExpr1 &&
           inductionProperty == candidateInductionProperty &&
           inductionBad == candidateInductionBad &&
           inductionOutputCount == candidateInductionOutputCount &&
           solverType == candidateSolverType;
  }
};

struct LazyTransitionStore {
  // Large SEC designs can have hundreds of thousands of modeled state bits.
  // K-induction proves one output cone at a time, so eagerly remapping every
  // state update into the shared symbol space wastes time and memory before
  // COI reduction can remove most of it. This store keeps the original local
  // next-state expressions plus the per-design symbol remap tables; the
  // k-induction encoders remap only transition equations that are actually
  // pulled into the current proof cone.
  std::unordered_map<size_t, LazyTransitionSource> sourceByStateSymbol;
  std::array<std::unordered_map<size_t, size_t>, 2> localToCombinedByDesign;
  std::array<std::unordered_map<size_t, DualRailSymbolPair>, 2>
      dualRailStateByLocalSymbolByDesign;
  mutable std::array<std::unordered_map<BoolExpr*, BoolExpr*>, 2> remapMemoByDesign;
  mutable std::array<std::unordered_map<BoolExpr*, DualRailBoolExpr>, 2>
      dualRailRemapMemoByDesign;
  mutable std::unordered_map<size_t, BoolExpr*> remappedByStateSymbol;
  // Output-batched SEC creates a fresh transition resolver for each PDR slice.
  // Keep lazy support and size metadata with the shared transition store so
  // reset-frontier COI rebuilding does not repeatedly walk the same large
  // source BoolExpr DAGs across batches.
  mutable std::unordered_map<size_t, std::set<size_t>> supportByStateSymbol;
  mutable std::unordered_map<size_t, size_t> nodeCountByStateSymbol;
  // Residual dual-rail KI slices are created outside the KI engine.  Remember
  // the full residual public-output conjunction here so later one-output leaves
  // can still use the same public induction hypothesis without adding any
  // cross-design internal-state relation.
  BoolExpr* dualRailResidualPublicProperty = nullptr;
  BoolExpr* dualRailResidualPublicBad = nullptr;
  size_t dualRailResidualPublicOutputCount = 0;
  // A cached UNKNOWN is never accepted as coverage.  It only avoids rerunning
  // the same strict one-output KI attempt after a residual batch has already
  // failed under the same public-output induction hypothesis.
  std::vector<DualRailResidualPublicKiAttempt>
      inconclusiveDualRailResidualPublicKiAttempts;

  bool findInconclusiveDualRailResidualPublicKiAttempt(
      BoolExpr* outputExpr0,
      BoolExpr* outputExpr1,
      BoolExpr* inductionProperty,
      BoolExpr* inductionBad,
      size_t inductionOutputCount,
      size_t solverType,
      size_t requestedMaxK,
      size_t& resultBound) const {
    for (const DualRailResidualPublicKiAttempt& attempt :
         inconclusiveDualRailResidualPublicKiAttempts) {
      if (attempt.matches(
              outputExpr0,
              outputExpr1,
              inductionProperty,
              inductionBad,
              inductionOutputCount,
              solverType) &&
          attempt.attemptedMaxK >= requestedMaxK) {
        resultBound = attempt.resultBound;
        return true;
      }
    }
    return false;
  }

  bool hasInconclusiveDualRailResidualPublicKiAttempt(
      BoolExpr* outputExpr0,
      BoolExpr* outputExpr1,
      BoolExpr* inductionProperty,
      BoolExpr* inductionBad,
      size_t inductionOutputCount,
      size_t solverType,
      size_t requestedMaxK) const {
    size_t ignoredBound = 0;
    return findInconclusiveDualRailResidualPublicKiAttempt(
        outputExpr0,
        outputExpr1,
        inductionProperty,
        inductionBad,
        inductionOutputCount,
        solverType,
        requestedMaxK,
        ignoredBound);
  }

  void rememberInconclusiveDualRailResidualPublicKiAttempt(
      BoolExpr* outputExpr0,
      BoolExpr* outputExpr1,
      BoolExpr* inductionProperty,
      BoolExpr* inductionBad,
      size_t inductionOutputCount,
      size_t solverType,
      size_t attemptedMaxK,
      size_t resultBound) {
    for (DualRailResidualPublicKiAttempt& attempt :
         inconclusiveDualRailResidualPublicKiAttempts) {
      if (attempt.matches(
              outputExpr0,
              outputExpr1,
              inductionProperty,
              inductionBad,
              inductionOutputCount,
              solverType)) {
        attempt.attemptedMaxK = // LCOV_EXCL_LINE
            std::max(attempt.attemptedMaxK, attemptedMaxK); // LCOV_EXCL_LINE
        attempt.resultBound = std::max(attempt.resultBound, resultBound); // LCOV_EXCL_LINE
        return; // LCOV_EXCL_LINE
      }
    }
    inconclusiveDualRailResidualPublicKiAttempts.push_back(
        DualRailResidualPublicKiAttempt{
            outputExpr0,
            outputExpr1,
            inductionProperty,
            inductionBad,
            inductionOutputCount,
            solverType,
            attemptedMaxK,
            resultBound});
  }
};

struct KInductionProblem {
  std::vector<SignalKey> environmentInputs;
  std::vector<SignalKey> observedOutputs;
  std::vector<std::string> environmentInputNames;
  std::vector<std::string> observedOutputNames;
  // Preserve the top-level SEC output width after batching/slicing. Some PDR
  // heuristics must size themselves by the original property, not by the
  // currently selected one-output leaf.
  size_t originalObservedOutputCount = 0;
  std::vector<size_t> inputSymbols;
  size_t resetBootstrapCycles = 0;
  std::vector<std::pair<size_t, bool>> resetBootstrapInputs;
  std::vector<std::pair<size_t, bool>> initialStateAssignments;
  std::vector<std::pair<size_t, bool>> bootstrapStateAssignments;
  std::vector<size_t> state0Symbols;
  std::vector<size_t> state1Symbols;
  // Verifier-owned monitor state is part of the proof transition system but
  // belongs to neither design.  Keeping it separate prevents accidental
  // cross-design state matching by name or position.
  std::vector<size_t> auxiliaryStateSymbols;
  std::vector<size_t> allSymbols;
  std::vector<std::pair<size_t, size_t>> complementedStatePairs0;
  std::vector<std::pair<size_t, size_t>> complementedStatePairs1;
  // Same-design state equalities that hold in every frame. Dual-rail SEC uses
  // this for Q/QN complemented state outputs, where the structural relation is
  // cross-rail equality rather than Boolean complement on one rail.
  std::vector<std::pair<size_t, size_t>> sameFrameStateEqualityPairs0;
  std::vector<std::pair<size_t, size_t>> sameFrameStateEqualityPairs1;
  std::vector<DualRailSymbolPair> dualRailStatePairs;
  std::vector<BoolExpr*> observedOutputExprs0;
  std::vector<BoolExpr*> observedOutputExprs1;
  // Exact rail equality is retained for shared SAT query surfaces. It is not
  // an equivalence criterion because matching 11 rails are still X.
  std::vector<BoolExpr*> dualRailOutputStrictEqualityExprs;
  std::vector<std::string> dualRailOutputSkipReasons;
  std::vector<std::pair<size_t, BoolExpr*>> transitions0;
  std::vector<std::pair<size_t, BoolExpr*>> transitions1;
  std::vector<std::pair<size_t, BoolExpr*>> auxiliaryTransitions;
  std::shared_ptr<LazyTransitionStore> lazyTransitions;
  BoolExpr* initialCondition = nullptr;
  size_t initializedStateCount = 0;
  size_t totalStateCount = 0;
  BoolExpr* property = nullptr;
  BoolExpr* bad = nullptr;
  BoolExpr* inductionProperty = nullptr;
  BoolExpr* inductionBad = nullptr;
  // Dual-rail SEC has a complete rail-valued boot state, but it still needs
  // the normal reset-bootstrap prefix so reset controls are driven exactly as
  // they are in the binary SEC flow.
  bool usesDualRailStateEncoding = false;
  // The second dual-rail SEC round proves strict equality of both rails. Its
  // recursive output splits use path-local incremental PDR solver contexts.
  bool usesStrictDualRailEqualityProperty = false;
  // Output-batched dual-rail KI proves each output slice independently.  When
  // this flag is set, the slice skips local base checks because the caller will
  // validate the shared full-output base prefix once after all slices prove.
  bool deferBaseCaseChecks = false;
  // KI output batching and increasing-k retries ask for the same transition
  // target supports many times.  Keep an exact per-problem cache of those DAG
  // walks so every retry still builds the same strict KI formula without
  // re-traversing unchanged transition cones.
  mutable std::shared_ptr<InductionTransitionSupportCache>
      inductionTransitionSupportCache;
  std::string description;

  bool hasSequentialState() const {
    return !state0Symbols.empty() || !state1Symbols.empty() ||
           !auxiliaryStateSymbols.empty();
  }

  bool hasExplicitInitialState() const {
    return initializedStateCount != 0;
  }

  bool hasCompleteInitialState() const {
    return initializedStateCount != 0 && initializedStateCount == totalStateCount;
  }

  bool hasResetBootstrap() const {
    return !resetBootstrapInputs.empty();
  }

  size_t effectiveTotalStateCount() const {
    return totalStateCount != 0 ? totalStateCount
                                : state0Symbols.size() + state1Symbols.size() +
                                      auxiliaryStateSymbols.size(); // LCOV_EXCL_LINE
  }

  bool hasCompleteBootstrapStateAssignments() const {
    const size_t stateCount = effectiveTotalStateCount();
    return stateCount != 0 && bootstrapStateAssignments.size() >= stateCount;
  }

  bool usesResetBootstrapObservationFrontier() const {
    // Binary SEC cannot compare resetless internal state across designs unless
    // that relation was proved. Dual rail already carries unknown startup state
    // on value/known rails, so forcing it onto this frontier can turn an
    // over-approximate startup state into a false counterexample.
    return !usesDualRailStateEncoding && hasSequentialState() &&
           hasResetBootstrap() && resetBootstrapCycles != 0 &&
           property != nullptr && !hasCompleteBootstrapStateAssignments();
  }

  bool canReportSteadyFrontierMismatchAsCounterexample() const { // LCOV_EXCL_LINE
    // With reset/bootstrap startup, the steady-frontier SAT model can be an
    // over-approximate startup assignment rather than a concrete design trace,
    // even when the rail-state bootstrap map is complete. Resetless sequential
    // SEC can still report a real top-output frontier mismatch through the
    // selected engine path.
    return !(hasSequentialState() && hasResetBootstrap() && // LCOV_EXCL_LINE
             resetBootstrapCycles != 0); // LCOV_EXCL_LINE
  }

  std::vector<size_t> combinedStateSymbols() const {
    std::vector<size_t> combined = state0Symbols;
    combined.insert(combined.end(), state1Symbols.begin(), state1Symbols.end());
    combined.insert(
        combined.end(),
        auxiliaryStateSymbols.begin(),
        auxiliaryStateSymbols.end());
    return combined;
  }
};

}  // namespace KEPLER_FORMAL::SEC
