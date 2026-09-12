// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "../../config/Config.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/KInductionProblem.h"
#include "kinduction/OutputBatching.h"

namespace KEPLER_FORMAL::SEC {

enum class IMCStatus {
  Equivalent,
  Different,
  Inconclusive,
};

struct IMCResult {
  IMCStatus status = IMCStatus::Inconclusive;
  size_t bound = 0;
  std::optional<KInductionResult::CounterexampleWitness> witness;
  // Large dual-rail Craig IMC records each proved output independently so an
  // inconclusive output does not hide later successful proofs.
  std::vector<bool> coveredOutputs;
};

// IMC-specific batching for large dual-rail Craig proofs. Unlike the shared
// KI batcher, this also checks marginal support growth so wide unrelated cones
// are split before launching an expensive interpolation query.
std::vector<std::pair<size_t, size_t>> buildLargeDualRailCraigImcOutputBatches(
    const KInductionProblem& problem,
    const OutputBatchingLimits& limits);

// Test-visible strict projection cap used by the large dual-rail Craig path.
size_t largeDualRailCraigImcProjectionStateLimit();

// Test-visible policy for combining two previously proved Craig helpers.  The
// merged helper is still strict IMC: both inputs are Craig-derived inductive
// invariants from the same SEC problem.
bool shouldCombineCraigHelpersForSmallRawSingleton(
    bool useSmallRawSingletonInvariant,
    bool reusableInvariantHasRegions);

// Interpolation-Based Model Checking over the extracted SEC problem. It keeps
// counterexample discovery on the shared BMC path, then grows a proof frontier
// with validated interpolation/reachability invariants until that frontier
// becomes inductive and excludes the bad states.
class IMCEngine {
 public:
  IMCEngine(const KInductionProblem& problem,
            KEPLER_FORMAL::Config::SolverType solverType);

  IMCResult run(size_t maxK) const;

 private:
  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
};

}  // namespace KEPLER_FORMAL::SEC
