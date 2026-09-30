// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchSupportOptions.h"
#include <utility>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace { thread_local SupportOptions currentOptions; }
const SupportOptions& supportOptions() { return currentOptions; }
ScopedSupportOptions::ScopedSupportOptions(SupportOptions options)
    : previous_(currentOptions) { currentOptions = std::move(options); }
ScopedSupportOptions::~ScopedSupportOptions() { currentOptions = std::move(previous_); }
}  // namespace KEPLER_FORMAL::SEC::LATCH
