// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "BoolExpr.h"
#include "common/SignalKey.h"

namespace KEPLER_FORMAL::SEC {

struct KInductionProblem {
  std::vector<SignalKey> environmentInputs;
  std::vector<SignalKey> observedOutputs;
  std::vector<std::string> environmentInputNames;
  std::vector<std::string> observedOutputNames;
  std::vector<size_t> inputSymbols;
  size_t resetBootstrapCycles = 0;
  std::vector<std::pair<size_t, bool>> resetBootstrapInputs;
  std::vector<std::pair<size_t, size_t>> initialStateEqualityPairs;
  std::vector<std::pair<size_t, bool>> bootstrapStateAssignments;
  std::vector<std::pair<size_t, size_t>> bootstrapStateEqualityPairs;
  std::vector<std::pair<size_t, size_t>> inductiveStateEqualityPairs;
  std::vector<size_t> state0Symbols;
  std::vector<size_t> state1Symbols;
  std::vector<size_t> allSymbols;
  std::vector<std::pair<size_t, size_t>> complementedStatePairs0;
  std::vector<std::pair<size_t, size_t>> complementedStatePairs1;
  std::vector<BoolExpr*> observedOutputExprs0;
  std::vector<BoolExpr*> observedOutputExprs1;
  std::vector<std::pair<size_t, BoolExpr*>> transitions0;
  std::vector<std::pair<size_t, BoolExpr*>> transitions1;
  BoolExpr* initialCondition = nullptr;
  size_t initializedStateCount = 0;
  size_t totalStateCount = 0;
  BoolExpr* property = nullptr;
  BoolExpr* bad = nullptr;
  BoolExpr* inductionProperty = nullptr;
  BoolExpr* inductionBad = nullptr;
  std::string description;

  bool hasSequentialState() const {
    return !state0Symbols.empty() || !state1Symbols.empty();
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

  std::vector<size_t> combinedStateSymbols() const {
    std::vector<size_t> combined = state0Symbols;
    combined.insert(combined.end(), state1Symbols.begin(), state1Symbols.end());
    return combined;
  }
};

}  // namespace KEPLER_FORMAL::SEC
