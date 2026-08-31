// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "phase/TernarySimulation.h"

namespace KEPLER_FORMAL::SEC::PHASE {

struct PeriodicSignalAnalysisOptions {
  // Bjesse and Kukula use a small global unfolding factor.  Eight is their
  // implementation heuristic; callers may deliberately choose another cap.
  size_t maxPhases = 8;
  // 0 means the TBB default, 1 forces deterministic serial execution.
  size_t maxConcurrency = 0;
};

struct PeriodicGenerator {
  size_t stateIndex = 0;
  size_t minimumPeriod = 0;
  std::vector<TernaryValue> generatorWord;
  std::vector<TernaryValue> expandedWord;
  bool usableAtSelectedPhaseCount = false;
};

struct PeriodicSignalStatistics {
  size_t candidateStateCount = 0;
  size_t resetPeriodicStateCount = 0;
  size_t usableGeneratorCount = 0;
};

struct PeriodicSignalAnalysisResult {
  // The selected global unfolding factor N. It is always at least one.
  size_t phaseCount = 1;
  std::vector<PeriodicGenerator> generators;
  PeriodicSignalStatistics statistics;
  bool complete = true;
  std::string diagnostic;
};

// Extract reset-anchored minimum generator words from a completed ternary
// trace. A candidate must remain Boolean and match its word across both the
// stem and cycle. Then select the smallest N <= maxPhases that maximizes the
// number of usable generators (a generator is usable when its minimum period
// divides N).
PeriodicSignalAnalysisResult analyzePeriodicSignals(
    const NormalizedTransitionSystem& system,
    const TernarySimulationResult& simulation,
    const PeriodicSignalAnalysisOptions& options = {});

}  // namespace KEPLER_FORMAL::SEC::PHASE
