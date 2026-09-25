// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <optional>

#include "latch/LatchSettlingCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {

// This is an explicit semantic contract, not an inference from cell names or
// clock carriers. A step is one external transaction followed by full settling.
struct SupportOptions {
  bool enabled = true;
  bool singleInputChange = false;
  std::optional<bool> initialInputs;
  std::optional<bool> initialStorage;
  size_t workers = 0;
  CompilerLimits limits;
  size_t maxSymbolicNodes = 2000000;
  unsigned maxSatConflicts = 500000;
  unsigned maxSatDecisions = 5000000;

  // Enabling support never invents power-up values. Without an explicit
  // contract, extraction retains the legacy model and opaque latch cones.
  bool hasEventContract() const {
    return enabled && initialInputs.has_value() && initialStorage.has_value();
  }
};

const SupportOptions& supportOptions();

// Extraction is serial; workers receive immutable copies, never this context.
// A scoped option also keeps embedded invocations from leaking semantics.
class ScopedSupportOptions {
 public:
  explicit ScopedSupportOptions(SupportOptions options);
  ~ScopedSupportOptions();
  ScopedSupportOptions(const ScopedSupportOptions&) = delete;
  ScopedSupportOptions& operator=(const ScopedSupportOptions&) = delete;
 private:
  SupportOptions previous_;
};

}  // namespace KEPLER_FORMAL::SEC::LATCH
