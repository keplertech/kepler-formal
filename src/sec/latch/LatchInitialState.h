// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "common/AlignedSignals.h"
#include "kinduction/KInductionProblem.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {

// At BOOT, an input's remembered level is its origin. Existentially eliminate
// the redundant fixed-once origin, without identifying future input levels.
void reuseInitialInputHistory(SequentialDesignModel& model,
    const std::vector<SignalKey>& inputs, const std::vector<BoolExpr*>& history);

// Preserve exact BOOT relations while sharing only initial external levels,
// never unrelated hardware storage, between the two designs.
void integrateEventInitialState(
    const SequentialDesignModel& first, const SequentialDesignModel& second,
    const AlignedSignals& inputs,
    const std::unordered_map<size_t, size_t>& symbols0,
    const std::unordered_map<size_t, size_t>& symbols1,
    KInductionProblem& problem);

// Event certification is Boolean. An unspecified Boolean origin is either 0
// or 1, not the ternary value X (both rails asserted).
void constrainEventInitialRails(
    const SequentialDesignModel& model,
    const std::unordered_map<size_t, DualRailSymbolPair>& rails,
    KInductionProblem& problem, bool secondDesign = false);

}  // namespace KEPLER_FORMAL::SEC::LATCH
