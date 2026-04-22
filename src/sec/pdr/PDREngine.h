// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

enum class PDRStatus {
  Equivalent,
  Different,
  Inconclusive,
};

struct PDRResult {
  PDRStatus status = PDRStatus::Inconclusive;
  size_t bound = 0;
};

// Top-level clause-based Property Directed Reachability strategy for SEC. It
// follows the classic proof-obligation/blocking loop over the already-built
// SEC transition system.
class PDREngine {
 public:
  PDREngine(const KInductionProblem& problem,
            KEPLER_FORMAL::Config::SolverType solverType);

  PDRResult run(size_t maxFrames) const;

 private:
  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
};

}  // namespace KEPLER_FORMAL::SEC
