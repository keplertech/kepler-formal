// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <optional>

#include "../../config/Config.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/KInductionProblem.h"

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
};

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
