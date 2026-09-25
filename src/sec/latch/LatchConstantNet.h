// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <optional>

#include "DNL.h"

namespace KEPLER_FORMAL::SEC::LATCH {

struct DriverlessConstant {
  std::optional<bool> value;
  bool conflictingOrUnknown = false;
};

// Some driverless terminals have no flattened iso even when connected to an
// explicitly constant net. Walk the entire connected hierarchy, not just its
// local net annotation: a second annotation can conflict with that constant.
// Scratch connectivity and default no-op callbacks leave the caller's DNL and
// source models unchanged. Ordinary driven/floating nets are never constants.
class DriverlessConstantResolver {
 public:
  explicit DriverlessConstantResolver(const naja::DNL::DNLFull& dnl)
      : builder_(scratch_, dnl) {}
  DriverlessConstant resolve(const naja::DNL::DNLTerminalFull& term) {
    naja::DNL::DNLComplexIso connected;
    builder_.treatDriver(term, connected, visited_);
    if (!connected.getDrivers().empty()) return {};
    DriverlessConstant result;
    for (const auto* net : connected.getNets()) {
      if (net->isConstantX() || net->isConstantZ()) return {{}, true};
      if (!net->isConstant0() && !net->isConstant1()) continue;
      const bool value = net->isConstant1();
      if (result.value && *result.value != value) return {{}, true};
      result.value = value;
    }
    return result;
  }
 private:
  naja::DNL::DNLIsoDB scratch_;
  naja::DNL::DNLIsoDBBuilder<naja::DNL::DNLInstanceFull, naja::DNL::DNLTerminalFull> builder_;
  // Reuse the traversal bitmap instead of allocating one
  // graph-sized visited set for each missing-iso terminal.
  naja::DNL::visited visited_;
};

}  // namespace KEPLER_FORMAL::SEC::LATCH
