// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <string>
#include "DNL.h"

namespace KEPLER_FORMAL {
class LeafBoundary;
}

namespace KEPLER_FORMAL::SEC {

struct SequentialDesignModel;

// Outputless primitive/black-box instances have no output key for the normal
// boundary collector. Scan the explicitly supplied graph only in strict mode.
// Top/ordinary empty hierarchy and explicitly selected boundaries are not
// implicitly reclassified as opaque cells.
void recordOpaqueOutputlessCells(SequentialDesignModel& model,
                                 const naja::DNL::DNLFull& dnl,
                                 const LeafBoundary* boundary = nullptr);

// Promote already-classified opacity to a hard unsupported result only when
// the opt-in policy is enabled. Call before pruning and after late extraction
// classifications: disconnected cells must not escape the policy.
void applyOpaquePolicy(SequentialDesignModel& model,
                       const std::string& topName, size_t side);

}  // namespace KEPLER_FORMAL::SEC
