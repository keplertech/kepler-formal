// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace KEPLER_FORMAL::SEC {

struct KInductionProblem;

struct SecBtor2Metadata {
  size_t totalOutputCount = 0;
  std::vector<std::string> skippedOutputs;
};

// Export the prepared SEC obligation with the concrete base-case startup and
// observation semantics. The artifact is bit-level and describes the selected
// output coverage, not a proof result or the engine's search heuristics.
void exportSecBtor2(const KInductionProblem& problem,
                   std::ostream& output,
                   const SecBtor2Metadata& metadata = {});

// Replace the destination only after a complete export has been written.
void exportSecBtor2File(const KInductionProblem& problem,
                       const std::string& path,
                       const SecBtor2Metadata& metadata = {});

}  // namespace KEPLER_FORMAL::SEC
