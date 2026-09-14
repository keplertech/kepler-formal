// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include "KeplerFormalDriver.h"

namespace KEPLER_FORMAL {

// The binding owns Python argument conversion, the GIL and Python exceptions.
// This adapter supplies the in-process policy to the common native run API.
int runPythonVerification(int argc, char** argv, RunResult& result);

}  // namespace KEPLER_FORMAL
