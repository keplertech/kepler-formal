// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "RunResult.h"

namespace KEPLER_FORMAL {

const char* runStatusName(RunStatus status) {
  switch (status) {
    case RunStatus::NoResult:
      return "no_result";
    case RunStatus::Equivalent:
      return "equivalent";
    case RunStatus::Different:
      return "different";
    case RunStatus::PartiallyProved:
      return "partially_proved";
    case RunStatus::Inconclusive:
      return "inconclusive";
    case RunStatus::Unsupported:
      return "unsupported";
    case RunStatus::Exported:
      return "exported";
    case RunStatus::Error:
    default:
      return "error";
  }
}

}  // namespace KEPLER_FORMAL
