// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <unordered_map>

#include "common/AlignedSignals.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC {

// The reachable-state strengthening that SEC derives before proving anything:
// which state pairs are compatible at startup, which ones stay anchored after
// reset/bootstrap, and which concrete state values become known after reset
// propagation.
struct ReachableStateInvariant {
  size_t bootstrapCycles = 0;
  AlignedSignals initialStateCorrespondence;
  AlignedSignals anchoredStateEqualities;
  std::unordered_map<SignalKey, bool, SignalKeyHash> bootstrapValues0;
  std::unordered_map<SignalKey, bool, SignalKeyHash> bootstrapValues1;
};

ReachableStateInvariant buildReachableStateInvariant(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& inductiveStateEqualities,
    bool secDiagEnabled = false);

}  // namespace KEPLER_FORMAL::SEC
