// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC {
struct SecResetSpec;
namespace LATCH {

// Owned expressions survive compact extraction. currentInputs describes the
// remembered external levels at a settled boundary, never fresh environment PIs.
struct EventResetInterface {
  std::vector<SignalKey> inputKeys;
  std::vector<std::string> inputNames;
  std::vector<BoolExpr*> currentInputs;
  bool singleInputChange = false;
  std::vector<size_t> selectorSymbols;
  std::optional<size_t> valueSymbol;
  std::optional<size_t> clockInputIndex;
  std::string clockError;
  size_t maxCompositionNodes = 2000000;
};

struct ResetCycleAdaptation {
  std::optional<SequentialDesignModel> model;
  std::string error;
};

// Compile a finite reset-cycle prefix by composing already-certified event
// transitions. Original model and ambient configuration remain unchanged.
ResetCycleAdaptation adaptResetCycles(const SequentialDesignModel& model,
                                     const SecResetSpec& reset);

}  // namespace LATCH
}  // namespace KEPLER_FORMAL::SEC
