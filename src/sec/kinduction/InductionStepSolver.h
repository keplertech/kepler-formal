// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <optional>

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

enum class InductionProofStatus {
  Proved,
  NotProved,
  Unknown,
};

InductionProofStatus proveByInductionStatus(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t k,
    std::optional<unsigned> kissatDecisionLimit = std::nullopt);

// Solves the k-induction step over a simple path of length k.
bool provesByInduction(const KInductionProblem& problem,
                       KEPLER_FORMAL::Config::SolverType solverType,
                       size_t k);

}  // namespace KEPLER_FORMAL::SEC
