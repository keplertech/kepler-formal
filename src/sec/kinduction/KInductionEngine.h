// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

enum class KInductionStatus {
  Equivalent,
  Different,
  Inconclusive,
};

struct KInductionResult {
  KInductionStatus status = KInductionStatus::Inconclusive;
  size_t bound = 0;
  struct SignalMismatch {  // LCOV_EXCL_LINE
    std::string signal;
    bool design0Value = false;
    bool design1Value = false;
  };
  struct FrameInputAssignments {  // LCOV_EXCL_LINE
    struct Assignment {  // LCOV_EXCL_LINE
      std::string signal;
      bool value = false;
    };

    size_t frame = 0;
    std::vector<Assignment> assignments;
  };
  struct CounterexampleWitness {  // LCOV_EXCL_LINE
    size_t badFrame = 0;
    std::vector<FrameInputAssignments> inputTrace;
    std::vector<SignalMismatch> outputMismatches;
  };

  std::optional<CounterexampleWitness> witness;
};

// Top-level k-induction strategy for SEC: a BMC-style base search for concrete
// counterexamples plus an induction step over simple paths.
class KInductionEngine {
 public:
  KInductionEngine(const KInductionProblem& problem,
                   KEPLER_FORMAL::Config::SolverType solverType);

  KInductionResult run(size_t maxK) const;

 private:
  const KInductionProblem& problem_;
  KEPLER_FORMAL::Config::SolverType solverType_;
};

}  // namespace KEPLER_FORMAL::SEC
