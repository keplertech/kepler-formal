// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "latch/LatchEventModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {

struct DependencyGraph {
  // Event-connected scheduling islands. All producer/consumer arcs count,
  // including clocks, enables and asynchronous controls; no implicit FF cut.
  // Multiple writers share an island so an invalid driver cannot be isolated
  // from another writer or its consumers. Read-only shared PIs are not edges.
  std::vector<std::vector<size_t>> components;
  // Directed strongly connected components with a cycle, including self-loops.
  // A structural cycle is not evidence of an active loop or a settling proof.
  std::vector<std::vector<size_t>> feedbackComponents;
  std::vector<size_t> componentOf;
};

// Indices refer to Network::primitives. Members and groups are ordered by their
// smallest primitive index. Analysis is iterative (including both SCC passes),
// invokes no primitive callbacks, and throws invalid_argument for bad net IDs.
DependencyGraph analyzeDependencies(const Network& network);

}  // namespace KEPLER_FORMAL::SEC::LATCH
