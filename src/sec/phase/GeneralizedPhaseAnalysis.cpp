// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "phase/GeneralizedPhaseAnalysis.h"

#include <sstream>
#include <utility>

namespace KEPLER_FORMAL::SEC::PHASE {

namespace {

std::string formatNormalizationDiagnostic(
    const NormalizationDiagnostic& diagnostic) {
  std::ostringstream out;
  out << (diagnostic.fatal ? "fatal" : "warning");
  if (!diagnostic.object.empty()) {
    out << " [" << diagnostic.object << "]";
  }
  out << ": " << diagnostic.message;
  return out.str();
}

}  // namespace

GeneralizedPhaseAnalysisResult analyzeGeneralizedPhases(
    const NormalizedTransitionSystem& system,
    const GeneralizedPhaseAnalysisOptions& options) {
  GeneralizedPhaseAnalysisResult result;
  for (const auto& diagnostic : system.diagnostics()) {
    result.diagnostics.push_back(formatNormalizationDiagnostic(diagnostic));
  }

  TernarySimulationOptions simulationOptions;
  simulationOptions.maxSteps = options.maxSimulationSteps;
  simulationOptions.maxStoredStateBytes = options.maxStoredStateBytes;
  simulationOptions.maxConcurrency = options.maxConcurrency;
  result.simulation = simulateTernary(system, simulationOptions);
  if (!result.simulation.complete()) {
    result.complete = false;
    if (!result.simulation.diagnostic.empty()) {
      result.diagnostics.push_back(result.simulation.diagnostic);
    }
    result.periodicSignals.complete = false;
    result.periodicSignals.diagnostic =
        "periodic-signal analysis was not run because simulation did not "
        "complete";
    result.latchProfiles.complete = false;
    result.latchProfiles.diagnostic =
        "latch profiling was not run because simulation did not complete";
    return result;
  }

  PeriodicSignalAnalysisOptions periodicOptions;
  periodicOptions.maxPhases = options.maxPhases;
  periodicOptions.maxConcurrency = options.maxConcurrency;
  result.periodicSignals =
      analyzePeriodicSignals(system, result.simulation, periodicOptions);
  if (!result.periodicSignals.complete) {
    result.complete = false;
    if (!result.periodicSignals.diagnostic.empty()) {
      result.diagnostics.push_back(result.periodicSignals.diagnostic);
    }
    result.latchProfiles.complete = false;
    result.latchProfiles.diagnostic =
        "latch profiling was not run because generator analysis did not "
        "complete";
    return result;
  }

  LatchPhaseProfileOptions profileOptions;
  profileOptions.maxConcurrency = options.maxConcurrency;
  result.latchProfiles =
      buildLatchPhaseProfiles(system, result.periodicSignals, profileOptions);
  if (!result.latchProfiles.complete) {
    result.complete = false;
    if (!result.latchProfiles.diagnostic.empty()) {
      result.diagnostics.push_back(result.latchProfiles.diagnostic);
    }
  }
  return result;
}

CollectedGeneralizedPhaseAnalysis analyzeLatchPhases(
    naja::NL::SNLDesign* top, const NormalizationOptions& normalizationOptions,
    const GeneralizedPhaseAnalysisOptions& analysisOptions) {
  CollectedGeneralizedPhaseAnalysis result;
  result.transitionSystem =
      collectNormalizedTransitionSystem(top, normalizationOptions);
  result.analysis =
      analyzeGeneralizedPhases(result.transitionSystem, analysisOptions);
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::PHASE
