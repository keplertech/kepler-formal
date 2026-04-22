// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>
#include <vector>

#include "common/SignalKey.h"

namespace KEPLER_FORMAL::SEC {

// Represents a paired SEC interface between the two designs. The labels are
// only for diagnostics; the real matching lives in the paired SignalKeys.
struct AlignedSignals {
  std::vector<std::string> names;
  std::vector<SignalKey> keys0;
  std::vector<SignalKey> keys1;
};

}  // namespace KEPLER_FORMAL::SEC
