// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <filesystem>

namespace naja::NL {
class NLLibrary;
}

namespace KEPLER_FORMAL {

// Ordinary Liberty loading plus explicit scalar latch semantics not yet
// supplied by the Naja frontend. Intended only for the opt-in latch-event path.
// Unsupported definitions remain unmodeled, with a warning; existing cells
// retain the ordinary constructor's first-definition-wins behavior.
void constructLibertyWithLatchModels(naja::NL::NLLibrary* library,
                                     const std::filesystem::path& path);

}  // namespace KEPLER_FORMAL
