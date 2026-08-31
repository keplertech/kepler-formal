// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "phase/LatchPhaseProfile.h"

namespace KEPLER_FORMAL::SEC::PHASE {

struct GeneralizedPhaseAnalysisOptions {
  size_t maxSimulationSteps = 4096;
  size_t maxStoredStateBytes = 256u * 1024u * 1024u;
  size_t maxPhases = 8;
  // 0 means TBB default concurrency. KEPLER_NO_MT always forces one worker.
  size_t maxConcurrency = 0;
};

struct GeneralizedPhaseAnalysisResult {
  TernarySimulationResult simulation;
  PeriodicSignalAnalysisResult periodicSignals;
  LatchPhaseProfileResult latchProfiles;
  bool complete = true;
  std::vector<std::string> diagnostics;
};

// Run the paper-aligned reset simulation and periodic generator analysis,
// followed by Kepler's latch-specific phase interpretation.
GeneralizedPhaseAnalysisResult analyzeGeneralizedPhases(
    const NormalizedTransitionSystem& system,
    const GeneralizedPhaseAnalysisOptions& options = {});

struct CollectedGeneralizedPhaseAnalysis {
  NormalizedTransitionSystem transitionSystem;
  GeneralizedPhaseAnalysisResult analysis;
};

// Convenience frontend. Netlist collection is serialized; the returned owned
// transition system and all subsequent analysis are independent of DNL.
CollectedGeneralizedPhaseAnalysis analyzeLatchPhases(
    naja::NL::SNLDesign* top,
    const NormalizationOptions& normalizationOptions = {},
    const GeneralizedPhaseAnalysisOptions& analysisOptions = {});

}  // namespace KEPLER_FORMAL::SEC::PHASE
