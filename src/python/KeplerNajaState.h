// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "DNL.h"

namespace KEPLER_FORMAL {

// Detach/restore the caller's DNL without constructing or deleting a graph.
// The caller must serialize access to Naja's process-global state.
naja::DNL::DNLFull* exchangeNajaDNL(naja::DNL::DNLFull* replacement) noexcept;

}  // namespace KEPLER_FORMAL
