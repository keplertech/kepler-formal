// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "property/PropertyCheckStrategy.h"

namespace KEPLER_FORMAL::SEC {

// Translator for the assertions the SystemVerilog frontend captures as
// naja_formal_property instances.
//
// The frontend stores an assertion as source text plus the signals it reads;
// it builds no logic for it. This translator reads the elaborated assertion
// the frontend kept for the instance, so the design must have been loaded
// with both preserveAssertions and keepASTLink. It compiles each Boolean
// operand into logic over the instance's DATA bits and writes the temporal
// structure in the form PropertyExpression.h describes.
PropertyTranslator makeSvaPropertyTranslator();

}  // namespace KEPLER_FORMAL::SEC
