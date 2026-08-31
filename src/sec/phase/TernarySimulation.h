// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "phase/NormalizedTransitionSystem.h"

namespace KEPLER_FORMAL::SEC::PHASE {

using TernaryState = std::vector<TernaryValue>;

enum class TernarySimulationStatus : uint8_t {
  Complete,
  InvalidSystem,
  StepLimitExceeded,
  MemoryLimitExceeded,
  EvaluationFailed,
};

struct TernarySimulationOptions {
  // Maximum number of transitions evaluated while looking for a repeated
  // ternary state.  The initial state does not consume a step.
  size_t maxSteps = 4096;
  // Approximate bytes available for stored state-vector payloads.  Zero means
  // no explicit memory limit.
  size_t maxStoredStateBytes = 256u * 1024u * 1024u;
  // Zero uses TBB's automatic concurrency; one is serial.  KEPLER_NO_MT always
  // forces one worker.
  size_t maxConcurrency = 0;
};

struct TernarySimulationStatistics {
  size_t simulatedSteps = 0;
  size_t storedStates = 0;
  size_t stateWidth = 0;
  size_t expressionNodes = 0;
  size_t effectiveMaxConcurrency = 1;
};

struct TernarySimulationResult {
  TernarySimulationStatus status = TernarySimulationStatus::InvalidSystem;
  // Contains the initial state and every subsequently reached unique state.
  // The repeated state that closes the cycle is not appended.
  std::vector<TernaryState> trace;
  size_t stemLength = 0;
  size_t cycleLength = 0;
  std::string diagnostic;
  TernarySimulationStatistics statistics;

  bool complete() const { return status == TernarySimulationStatus::Complete; }
};

// Evaluate one owned expression using strong Kleene three-valued logic.
// Throws std::invalid_argument for a malformed environment or expression.
TernaryValue evaluateTernaryExpression(
    const NormalizedTransitionSystem& system, ExprID expression,
    const std::vector<TernaryValue>& variableValues);

// Reset-relative deterministic X simulation.  Time steps are sequential;
// independent next-state roots are evaluated in parallel.
TernarySimulationResult simulateTernary(
    const NormalizedTransitionSystem& system,
    const TernarySimulationOptions& options = {});

const char* toString(TernarySimulationStatus status);

}  // namespace KEPLER_FORMAL::SEC::PHASE
