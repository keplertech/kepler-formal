// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <vector>

#include "BoolExpr.h"
#include "clocks/SecClockModel.h"
#include "common/SignalKey.h"
#include "../../utils/DesignBoundary.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL::SEC {
namespace LATCH { struct EventResetInterface; }

struct ComplementedStateRelation {  // LCOV_EXCL_LINE
  SignalKey primaryKey;
  SignalKey complementedKey;
};

enum class ConnectivitySkipOrigin {
  NoDriver,
  MultiDriver,
  LogicalLoop,
  MultiClockDomain,
  OpaqueInternal,
  UnknownConstant,
};

struct ConnectivitySkipInfo {  // LCOV_EXCL_LINE
  ConnectivitySkipOrigin origin = ConnectivitySkipOrigin::NoDriver;
  std::string detail;
};

// Normalized view of a sequential design after extracting the interface we
// need for SEC: environment inputs, current-state bits, observed outputs, and
// the Boolean formulas that describe outputs and next-state updates.
struct SequentialDesignModel {  // LCOV_EXCL_LINE
  std::vector<SignalKey> environmentInputs;
  std::vector<SignalKey> stateBits;
  std::vector<SignalKey> topInputKeys;
  std::vector<SignalKey> topOutputKeys;
  std::vector<SignalKey> allObservedOutputs;
  std::vector<SignalKey> observedOutputs;
  std::vector<SignalKey> skippedStateBits;
  std::vector<SignalKey> skippedObservedOutputs;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputVarByKey;
  std::unordered_map<SignalKey, std::string, SignalKeyHash> displayNameByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> observedOutputExprByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> nextStateExprByStateKey;
  std::unordered_map<SignalKey, bool, SignalKeyHash> initialStateValueByKey;
  // Variables proven during extraction to be pure routed clock carriers.
  // Downstream SEC matching can classify them with the top clock without
  // making any name-based assumption about internal sequential state.
  std::vector<size_t> clockCarrierVarIDs;
  std::vector<ClockCarrierClass> clockCarrierClasses;
  std::unordered_map<SignalKey, ClockEvent, SignalKeyHash> clockEventByStateKey;
  std::unordered_map<SignalKey, ConnectivitySkipInfo, SignalKeyHash>
      connectivitySkipInfoByKey;
  std::vector<ComplementedStateRelation> complementedStateRelations;
  std::vector<std::string> unsupportedReasons;
  // Empty for legacy clock-cycle extraction; event models retain their contract
  // after compact mode releases the netlists.
  std::string eventContract;
  // Copied input-level/clock metadata for reset-cycle expansion after compact
  // extraction; never retains pointers into the source netlist.
  std::shared_ptr<const LATCH::EventResetInterface> eventResetInterface;
  size_t eventResetCycles = 0;

  // Extract the model from the given top design. Opaque per-output cones are
  // skipped; globally unsupported structures are recorded in unsupportedReasons.
  static SequentialDesignModel extract(naja::NL::SNLDesign* top,
                                       const BoundaryPairs& pairs = {},
                                       size_t side = 0);

  bool hasUnsupportedFeatures() const {
    return !unsupportedReasons.empty();
  }

  size_t coveredObservedOutputCount() const { return observedOutputs.size(); }
  size_t totalObservedOutputCount() const { return allObservedOutputs.size(); }
};

}  // namespace KEPLER_FORMAL::SEC
