// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
// Empty when disabled, or when a latch-free design has no explicitly requested
// event configuration. Omitted BOOT values stay symbolic, never default to zero.
std::optional<SequentialDesignModel> extractEventDesign(
    naja::NL::SNLDesign* top, const BoundaryPairs& pairs, size_t side);
}  // namespace KEPLER_FORMAL::SEC::LATCH
