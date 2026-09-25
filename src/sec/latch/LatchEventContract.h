// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <optional>
#include <string>

namespace KEPLER_FORMAL::SEC::LATCH {
// Also used for already-extracted/borrowed models: CLI validation alone cannot
// protect those paths from mixing cycles, events, or initial-state contracts.
inline std::optional<std::string> eventContractError(
    const std::string& left, const std::string& right, bool resetCycles) {
  if (left != right)
    return "SEC models use different clock/event or Boolean initialization contracts";
  if (!left.empty() && resetCycles)
    return "SEC external-event models cannot use clock-cycle reset bootstrap";
  return {};
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
