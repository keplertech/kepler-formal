// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <optional>

#include "../../config/Config.h"
#include "BoolExpr.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

// A small-state exact interpolant helper used by IMC and the legacy SEC flow
// to strengthen induction. It enumerates shared-state valuations exactly, so
// it is only enabled for compact transition systems where the symbolic support
// stays manageable.
class ExactInterpolantSynthesizer {
 public:
  ExactInterpolantSynthesizer(const KInductionProblem& problem,
                              KEPLER_FORMAL::Config::SolverType solverType);

  std::optional<BoolExpr*> deriveOneStepReachableStateInvariant(
      size_t maxSharedStateBits = 12) const;

 private:
  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
};

}  // namespace KEPLER_FORMAL::SEC
