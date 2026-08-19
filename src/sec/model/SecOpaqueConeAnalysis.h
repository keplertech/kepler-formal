// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>
#include <vector>

#include "DNL.h"

namespace KEPLER_FORMAL::SEC {

struct OpaqueTerminalSeed {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  std::string detail;
};

struct SecLocalTerminalDependency {
  naja::DNL::DNLID sourceTermID = naja::DNL::DNLID_MAX;
  naja::DNL::DNLID outputTermID = naja::DNL::DNLID_MAX;
};

struct OpaqueReachedTopOutput {
  naja::DNL::DNLID topOutputTermID = naja::DNL::DNLID_MAX;
  OpaqueTerminalSeed source;
};

std::vector<OpaqueReachedTopOutput> findTopOutputsReachedByOpaqueTerminals(
    naja::DNL::DNLFull* dnl, std::vector<OpaqueTerminalSeed> opaqueSeeds,
    const std::vector<SecLocalTerminalDependency>& extraDependencies = {});

}  // namespace KEPLER_FORMAL::SEC
