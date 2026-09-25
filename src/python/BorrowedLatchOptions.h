// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>

#include "latch/LatchSupportOptions.h"

namespace KEPLER_FORMAL {

enum class LatchInputChanges { Any, Single };

// Optional members distinguish an explicitly stated event contract from an
// omitted setting. Omitted initialization stays symbolic; omitted inputChanges
// permits any input changes. Overrides never become implicit assumptions.
struct BorrowedLatchOptions {
  bool enabled = true;
  std::optional<LatchInputChanges> inputChanges;
  std::optional<bool> initialInputs;
  std::optional<bool> initialStorage;
  std::optional<size_t> workers;
  std::optional<size_t> maxWaves;
  std::optional<size_t> maxStates;
  std::optional<size_t> maxTransactions;
  std::optional<size_t> maxSymbolicNodes;
  std::optional<size_t> maxSatConflicts;
  std::optional<size_t> maxSatDecisions;

  // Throws invalid_argument without changing any process or netlist state.
  SEC::LATCH::SupportOptions validated(bool isSec, bool hasLeafBoundaries) const;
};

}  // namespace KEPLER_FORMAL
