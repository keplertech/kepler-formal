// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "BoolExpr.h"
#include "common/SignalKey.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL::SEC {

struct ComplementedStateRelation {
  SignalKey primaryKey;
  SignalKey complementedKey;
};

enum class ConnectivitySkipOrigin {
  NoDriver,
  MultiDriver,
};

struct ConnectivitySkipInfo {
  ConnectivitySkipOrigin origin = ConnectivitySkipOrigin::NoDriver;
  std::string detail;
};

// Normalized view of a sequential design after extracting the interface we
// need for SEC: environment inputs, current-state bits, observed outputs, and
// the Boolean formulas that describe outputs and next-state updates.
struct SequentialDesignModel {
  std::vector<SignalKey> environmentInputs;
  std::vector<SignalKey> stateBits;
  std::vector<SignalKey> allObservedOutputs;
  std::vector<SignalKey> observedOutputs;
  std::vector<SignalKey> skippedStateBits;
  std::vector<SignalKey> skippedObservedOutputs;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputVarByKey;
  std::unordered_map<SignalKey, std::string, SignalKeyHash> displayNameByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> observedOutputExprByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> nextStateExprByStateKey;
  std::unordered_map<SignalKey, bool, SignalKeyHash> initialStateValueByKey;
  std::unordered_map<SignalKey, ConnectivitySkipInfo, SignalKeyHash>
      connectivitySkipInfoByKey;
  std::vector<ComplementedStateRelation> complementedStateRelations;
  std::vector<std::string> abstractedSequentialBoundaries;
  std::vector<std::string> unsupportedReasons;

  // Extract the model from the given top design. Unsupported sequential
  // structures are recorded in unsupportedReasons instead of being guessed.
  static SequentialDesignModel extract(naja::NL::SNLDesign* top);

  bool hasUnsupportedFeatures() const {
    return !unsupportedReasons.empty();
  }

  size_t coveredObservedOutputCount() const { return observedOutputs.size(); }
  size_t totalObservedOutputCount() const { return allObservedOutputs.size(); }
};

}  // namespace KEPLER_FORMAL::SEC
