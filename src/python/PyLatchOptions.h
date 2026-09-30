// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <Python.h>
#include <string_view>

#include "BorrowedLatchOptions.h"

namespace KEPLER_FORMAL {

bool isBorrowedLatchOption(std::string_view key);
bool parseBorrowedLatchOptions(PyObject* dictionary, BorrowedLatchOptions& options,
                              bool isSec, bool hasLeafBoundaries);

}  // namespace KEPLER_FORMAL
