// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>

#include "DNL.h"

namespace KEPLER_FORMAL::SEC {

struct OpaqueTerminalSeed {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  std::string detail;
};

struct OpaqueReachedTopOutput {
  naja::DNL::DNLID topOutputTermID = naja::DNL::DNLID_MAX;
  OpaqueTerminalSeed source;
};

class SecNetlistChecks {
 public:
  explicit SecNetlistChecks(naja::DNL::DNLFull* dnl) : dnl_(dnl) {}

  std::vector<OpaqueReachedTopOutput> findTopOutputsReachedByOpaqueTerminals(
      std::vector<OpaqueTerminalSeed> opaqueSeeds) const;

 private:
  naja::DNL::DNLFull* dnl_ = nullptr;
};

}  // namespace KEPLER_FORMAL::SEC
