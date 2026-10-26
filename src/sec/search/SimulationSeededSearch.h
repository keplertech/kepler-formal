// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

#include "../../config/Config.h"
#include "common/AlignedSignals.h"
#include "common/SignalKey.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/KInductionProblem.h"
#include "model/SequentialDesignModel.h"
#include "search/SimulationSeededSearchOptions.h"

namespace KEPLER_FORMAL::SEC {

// Three-valued simulation value in the two-rail encoding used by the proof
// engines: bit 0 is "may be 1", bit 1 is "may be 0".
enum class TernaryValue : uint8_t {
  One = 1,
  Zero = 2,
  Unknown = 3,
};

inline bool isKnownTernaryValue(TernaryValue value) {
  return value == TernaryValue::One || value == TernaryValue::Zero;
}

// Evaluates one extracted design's next-state functions on three-valued
// state with (optionally unknown) inputs. Every gate is lifted exactly as the
// dual-rail proof encoding lifts it, so a simulated state is a state of the
// dual-rail transition system itself (ternary simulation after Bryant,
// IEEE Trans. Computers 1987; Seger and Bryant, FMSD 1995).
class TernaryModelSimulator {
 public:
  explicit TernaryModelSimulator(const SequentialDesignModel& model);

  const SequentialDesignModel& model() const { return model_; }
  size_t stateCount() const { return state_.size(); }

  // Initial values from the model; every other state bit is unknown.
  void resetToInitialState();
  void setInput(size_t localVar, TernaryValue value);
  void clearInputs();
  // One clock cycle: next state from the current state and inputs.
  // Returns whether any state bit changed.
  bool step();

  TernaryValue stateValue(size_t stateIndex) const {
    return state_[stateIndex];
  }
  const std::vector<TernaryValue>& state() const { return state_; }
  void setState(const std::vector<TernaryValue>& state);
  // State bits that changed between two known values in any step since the
  // last clearToggles().
  const std::vector<bool>& toggled() const { return toggled_; }
  void clearToggles();
  size_t knownStateBitCount() const;
  // Observed output value for the current state and inputs.
  TernaryValue evaluateOutput(const SignalKey& outputKey);

 private:
  struct Node {
    uint8_t op = 0;
    uint32_t left = 0;
    uint32_t right = 0;
  };

  uint32_t compile(BoolExpr* expr);
  void ensureVarCapacity(size_t localVar);
  void evaluateCombinational();

  const SequentialDesignModel& model_;
  std::vector<Node> nodes_;
  std::vector<uint8_t> values_;
  std::unordered_map<BoolExpr*, uint32_t> nodeIndexByExpr_;
  std::vector<uint32_t> leafNodeByLocalVar_;
  std::vector<size_t> stateIndexByLocalVar_;
  std::vector<uint8_t> inputValueByLocalVar_;
  std::vector<uint32_t> nextRootByState_;
  std::unordered_map<SignalKey, uint32_t, SignalKeyHash> outputRootByKey_;
  std::vector<TernaryValue> state_;
  std::vector<bool> toggled_;
  bool combinationalValid_ = false;
};

// Runs the simulation-seeded search on a dual-rail SEC problem. The rails
// maps give each design's state keys their proof symbols. A returned witness
// is a confirmed trace: the simulated prefix is replayed together with the
// bounded-search suffix on the exact models before it is reported.
std::optional<KInductionResult::CounterexampleWitness>
runSimulationSeededSearch(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>&
        rails0ByKey,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>&
        rails1ByKey,
    const KInductionProblem& problem,
    const SimulationSeededSearchOptions& options,
    KEPLER_FORMAL::Config::SolverType solverType,
    bool diagnostics);

}  // namespace KEPLER_FORMAL::SEC
