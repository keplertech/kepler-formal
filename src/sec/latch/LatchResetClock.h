// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>
#include <string>
#include <vector>
#include "BoolExpr.h"
#include "latch/LatchEventModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {

// Copied semantic metadata, parallel to Network::primitives. Expressions use
// local input symbols pinIndex + 2; symbols 0 and 1 are Boolean constants.
// Latch enables are deliberately not clock roots. This discovery schedules
// reset cycles; it does not classify or replace latch event semantics.
struct ResetClockPrimitive {
  enum class Kind { Combinational, Latch, FlipFlop, Unsupported };
  Kind kind = Kind::Unsupported;
  BoolExpr* clock = nullptr;
  std::vector<BoolExpr*> outputs;
};

struct ResetClockDiscovery {
  enum class Status { Resolved, NoEdgeClock, Unsupported };
  Status status = Status::Unsupported;
  std::optional<size_t> rootNet;
  std::optional<size_t> externalInputIndex;
  std::string detail;

  bool resolved() const { return status == Status::Resolved; }
};

// Finds one common external carrier for EVERY explicit flip-flop clock. Only
// exact constant/literal reductions through combinational gates are accepted.
// An arbitrary multi-input gate, state-generated clock, or second root cannot
// be guessed away. Both edge polarities on the same root are supported: the
// reset sequencer emits complete source-clock cycles, settling after each edge.
// The implementation is iterative, including for cyclic or very deep routing.
ResetClockDiscovery discoverResetClock(
    const Network& network, const std::vector<ResetClockPrimitive>& metadata);

}  // namespace KEPLER_FORMAL::SEC::LATCH
