// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

struct OutputBatchingLimits {
  size_t maxOutputBatchSize = 32;
  size_t outputBatchSupportLimit = 512;
};

// Build small SEC output slices that can be proved independently and combined.
// Each slice is still a real SEC property: the full conjunction is valid iff
// every sliced conjunction is valid.  The split only gives the SAT encoders a
// tighter cone of influence for each query.
std::vector<std::pair<size_t, size_t>> buildSupportBoundedOutputBatches(
    const KInductionProblem& problem);

std::vector<std::pair<size_t, size_t>> buildSupportBoundedOutputBatches(
    const KInductionProblem& problem,
    const OutputBatchingLimits& limits);

OutputBatchingLimits defaultOutputBatchingLimitsForProblem(
    const KInductionProblem& problem);

// Reuse the large shared SEC problem object while replacing only the observed
// output/property slice.  This avoids copying all state, transition, and memory
// metadata per output batch.
void configureOutputBatchProblem(KInductionProblem& batch,
                                 const KInductionProblem& source,
                                 size_t firstOutput,
                                 size_t endOutput);

}  // namespace KEPLER_FORMAL::SEC
