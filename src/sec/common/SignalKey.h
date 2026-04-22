// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstdint>
#include <sstream>
#include <string>
#include <string_view>

#include "../../strategies/miter/BuildPrimaryOutputClauses.h"

namespace KEPLER_FORMAL::SEC {

using SignalKey = KEPLER_FORMAL::BuildPrimaryOutputClauses::PathKey;
using SignalKeyHash = KEPLER_FORMAL::BuildPrimaryOutputClauses::KeyHash;

struct SignalKeyLess {
  bool operator()(const SignalKey& lhs, const SignalKey& rhs) const {
    return lhs < rhs;
  }
};

inline naja::NL::NLID::DesignObjectID stableSignalKeyNameID(
    std::string_view name) {
  // Use a deterministic string hash so SEC aligns signals by hierarchical
  // path names instead of parser-local object IDs, which can drift after
  // structural edits even when the user-visible signal names stay the same.
  constexpr uint64_t kFnvOffset = 1469598103934665603ull;
  constexpr uint64_t kFnvPrime = 1099511628211ull;

  uint64_t hash = kFnvOffset;
  for (const unsigned char ch : name) {
    hash ^= ch;
    hash *= kFnvPrime;
  }
  return static_cast<naja::NL::NLID::DesignObjectID>(hash);
}

inline std::string signalKeyToString(const SignalKey& key) {
  std::ostringstream oss;
  for (const auto& nameID : key.first) {
    oss << nameID << ".";
  }
  for (const auto& objectID : key.second) {
    oss << objectID << ".";
  }
  return oss.str();
}

}  // namespace KEPLER_FORMAL::SEC
