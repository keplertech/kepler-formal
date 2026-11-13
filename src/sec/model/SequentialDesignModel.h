// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <set>
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

// The clock of a memory as the netlist shows it: the top-level inputs that
// drive its clock pins through combinational cells. `internalSource` is set
// when a pin is driven some other way, for example from a register, a
// constant or a net without a single driver.
struct MemoryClockSources {  // LCOV_EXCL_LINE
  std::vector<SignalKey> topInputs;
  bool internalSource = false;

  bool operator==(const MemoryClockSources&) const = default;
};

// A structured memory compared at its interface instead of modeled cell by
// cell: its connected input pins are observed outputs and its read-data pins
// are environment inputs, all keyed by the memory's hierarchical path. Two
// designs that drive such a memory identically read identical data from it.
struct AbstractedMemory {  // LCOV_EXCL_LINE
  SignalKey instanceKey;
  std::string displayName;
  // Primitive name and instance parameters: equal for identical memories.
  std::string signature;
  size_t width = 0;
  size_t depth = 0;
  std::vector<SignalKey> interfaceKeys;
  std::vector<SignalKey> readDataKeys;
  std::vector<SignalKey> writeEnableKeys;
  MemoryClockSources clockSources;
  // Why comparing the pins cannot stand in for this memory's cells, when it
  // cannot: the memory then has to be extracted again with its cells.
  std::string obstacle;
  // The ports whose data follows their address within the step. A path from
  // such data back to its own address is a combinational loop.
  struct AsynchronousReadPort {  // LCOV_EXCL_LINE
    std::vector<SignalKey> addressKeys;
    std::vector<SignalKey> dataKeys;
  };
  std::vector<AsynchronousReadPort> asynchronousReadPorts;
};

// A structured memory modeled with one state bit per cell bit.
struct BitLevelMemory {  // LCOV_EXCL_LINE
  std::string displayName;
  size_t width = 0;
  size_t depth = 0;
  // Why the extraction did not abstract this memory although it was asked to,
  // when the memory itself is the reason. Empty otherwise.
  std::string abstractionObstacle;
  MemoryClockSources clockSources;
  // The states of this memory: its cells, whose keys share this path part,
  // and its clocked read outputs.
  decltype(SignalKey::first) cellKeyPath;
  std::vector<SignalKey> readOutputKeys;
};

struct SequentialExtractOptions {  // LCOV_EXCL_LINE
  // Abstract the structured memories as described on AbstractedMemory.
  bool abstractMemories = false;
  // When set, only the memories with these instance keys are abstracted.
  std::optional<std::set<SignalKey>> abstractOnly;

  bool operator==(const SequentialExtractOptions&) const = default;
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
  // Name of the net on a state bit's output pin, when that net is named: the
  // parent instance path, the net name and the bit. Only a guess source for
  // register correspondence (netlist flows keep net names where they rename
  // or re-cluster instances); it never takes part in a proof.
  std::unordered_map<SignalKey, std::string, SignalKeyHash> outputNetNameByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> observedOutputExprByKey;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> nextStateExprByStateKey;
  // Read data of asynchronously read memories: the cell addressed in the current
  // step. The formula replaces the read-data variable everywhere and the
  // variable leaves the state.
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> memoryReadExprByKey;
  std::unordered_map<SignalKey, bool, SignalKeyHash> initialStateValueByKey;
  // Variables proven during extraction to be pure routed clock carriers.
  // Downstream SEC matching can classify them with the top clock without
  // making any name-based assumption about internal sequential state.
  std::vector<size_t> clockCarrierVarIDs;
  std::vector<ClockCarrierClass> clockCarrierClasses;
  std::unordered_map<SignalKey, ClockEvent, SignalKeyHash> clockEventByStateKey;
  // The clock domain of the registers an observed output reads, for an output
  // whose cone stays in one domain. The model has one step for all clocks, so
  // SEC compares an output across two designs only when both clock it from
  // the same domain.
  std::unordered_map<SignalKey, SignalKey, SignalKeyHash>
      clockDomainByObservedOutputKey;
  std::unordered_map<SignalKey, ConnectivitySkipInfo, SignalKeyHash>
      connectivitySkipInfoByKey;
  std::vector<ComplementedStateRelation> complementedStateRelations;
  std::vector<AbstractedMemory> abstractedMemories;
  std::vector<BitLevelMemory> bitLevelMemories;
  // An observed output that matters only on the frames where every listed
  // observed output is 1: the address, mask and data pins of an abstracted
  // memory write port have no effect while the port is not writing. SEC
  // compares such an output only on those frames.
  std::unordered_map<SignalKey, std::vector<SignalKey>, SignalKeyHash>
      observedOutputGuardKeysByKey;
  std::vector<std::string> unsupportedReasons;

  // Extract the model from the given top design. Opaque per-output cones are
  // skipped; globally unsupported structures are recorded in unsupportedReasons.
  static SequentialDesignModel extract(naja::NL::SNLDesign* top,
                                       const BoundaryPairs& pairs = {},
                                       size_t side = 0,
                                       const SequentialExtractOptions& options = {});

  bool hasUnsupportedFeatures() const {
    return !unsupportedReasons.empty();
  }

  size_t coveredObservedOutputCount() const { return observedOutputs.size(); }
  size_t totalObservedOutputCount() const { return allObservedOutputs.size(); }
};

}  // namespace KEPLER_FORMAL::SEC
