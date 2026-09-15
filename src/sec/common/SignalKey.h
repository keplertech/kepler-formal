// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <sstream>
#include <string>

#include "../../strategies/miter/BuildPrimaryOutputClauses.h"

namespace KEPLER_FORMAL::SEC {

using SignalKey = KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey;
using SignalKeyHash = KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash;

struct SignalKeyLess {
  bool operator()(const SignalKey& lhs, const SignalKey& rhs) const {
    return lhs < rhs;
  }
};

// LCOV_EXCL_START
inline std::string signalKeyToString(const SignalKey& key) {
  std::ostringstream oss;
  for (const auto& nameID : key.first) {
    oss << nameID << ".";
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  for (const auto& objectID : key.second) {
    oss << objectID << ".";
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return oss.str();
}
// LCOV_EXCL_STOP

}  // namespace KEPLER_FORMAL::SEC
