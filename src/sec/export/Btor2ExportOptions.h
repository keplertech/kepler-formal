// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <string>

namespace KEPLER_FORMAL::SEC {

struct Btor2ExportOptions {
  std::string path;
  bool dumpOnly = false;

  bool enabled() const { return !path.empty(); }
};

}  // namespace KEPLER_FORMAL::SEC
