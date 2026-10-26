// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
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
// Word-level view of one structured memory instance, kept next to the
// bit-level cell transitions. Engines that search bounded traces can model
// the array through its ports (a read returns the latest earlier write to the
// same address, read-over-write after Burch and Dill, CAV 1994) instead of
// carrying every cell through every frame.
struct StructuredMemoryReadPortModel {  // LCOV_EXCL_LINE
  std::vector<BoolExpr*> addressExprs;  // least significant bit first
  std::vector<SignalKey> outputKeys;    // registered read data, per bit
};

struct StructuredMemoryWritePortModel {  // LCOV_EXCL_LINE
  std::vector<BoolExpr*> addressExprs;
  std::vector<BoolExpr*> dataExprs;
  std::vector<BoolExpr*> maskExprs;    // empty: every bit is written
  std::vector<BoolExpr*> enableExprs;  // all must hold for a write
};

struct StructuredMemoryModel {  // LCOV_EXCL_LINE
  size_t width = 0;
  size_t depth = 0;
  size_t abits = 0;
  std::vector<std::vector<SignalKey>> cellKeys;  // [cell][bit]
  std::vector<StructuredMemoryReadPortModel> readPorts;
  std::vector<StructuredMemoryWritePortModel> writePorts;
  // Asserted reset clears every cell and read output after the writes.
  BoolExpr* resetAssertedExpr = nullptr;
};

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
  // Name of the net on a state bit's output pin, when that net is named: the
  // parent instance path, the net name and the bit. Only a guess source for
  // register correspondence (netlist flows keep net names where they rename
  // or re-cluster instances); it never takes part in a proof.
  std::unordered_map<SignalKey, std::string, SignalKeyHash> outputNetNameByKey;
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
  // Port-level descriptions of the structured memories whose cells are in
  // stateBits; skipped ports and disabled write ports are left out.
  std::vector<StructuredMemoryModel> structuredMemories;

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
