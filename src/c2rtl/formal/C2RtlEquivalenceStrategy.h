// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "../../config/Config.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL::C2RTL {

struct C2RtlEquivalenceOptions {
  // Delay zero compares the current combinational reference value. Delay N
  // compares against the reference value sampled N active clock events ago.
  std::unordered_map<std::string, size_t> outputDelays;
};

enum class C2RtlEquivalenceStatus {
  Equivalent,
  Different,
  Inconclusive,
  Unsupported,
};

struct C2RtlEquivalenceResult {
  C2RtlEquivalenceStatus status = C2RtlEquivalenceStatus::Unsupported;
  size_t bound = 0;
  size_t comparedBits = 0;
  size_t comparedOutputs = 0;
  std::string reason;
  std::string clock;
  std::string reset;
  bool resetActiveHigh = true;
  std::vector<std::pair<std::string, size_t>> outputDelays;
};

// Builds an isolated temporal miter around an XLS-generated combinational RTL
// reference and a user-provided clocked RTL implementation. The resulting
// transition system is handed to the existing PDR API without changing PDR or
// the normal sequential-equivalence construction path.
class C2RtlEquivalenceStrategy {
public:
  C2RtlEquivalenceStrategy(naja::NL::SNLDesign *combinationalReference,
                           naja::NL::SNLDesign *clockedImplementation,
                           KEPLER_FORMAL::Config::SolverType solverType,
                           C2RtlEquivalenceOptions options);

  C2RtlEquivalenceResult run(size_t maxFrames) const;

private:
  naja::NL::SNLDesign *combinationalReference_ = nullptr;
  naja::NL::SNLDesign *clockedImplementation_ = nullptr;
  KEPLER_FORMAL::Config::SolverType solverType_;
  C2RtlEquivalenceOptions options_;
};

} // namespace KEPLER_FORMAL::C2RTL
