// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {
// Empty only when the explicit event contract is disabled. In event mode even
// a latch-free comparison side must use the same external-event semantics.
std::optional<SequentialDesignModel> extractEventDesign(
    naja::NL::SNLDesign* top, const BoundaryPairs& pairs, size_t side);
}  // namespace KEPLER_FORMAL::SEC::LATCH
