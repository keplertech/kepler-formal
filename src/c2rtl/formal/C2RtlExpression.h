// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace KEPLER_FORMAL {
class BoolExpr;
}

namespace KEPLER_FORMAL::C2RTL {

// Unsigned values, least significant bit first. Comparisons zero-extend the
// shorter operand; Boolean operators treat any nonzero value as true.
using ExpressionBits = std::vector<BoolExpr*>;
using ExpressionSignalResolver = std::function<ExpressionBits(
    const std::string& name, std::optional<size_t> bitIndex)>;

// Supported syntax: signal names, optional [decimal bit index], nonnegative
// decimal/hexadecimal/binary literals, true/false, !, &&, ||, ==, !=, <, >, <=,
// >= and parentheses, with C operator precedence. The resolver receives the
// complete name (including model./rtl.) and enforces the caller's signal scope.
// All operands are resolved, including those behind a constant true/false.
// Invalid syntax or signal resolution throws std::invalid_argument with a byte
// offset. Expressions are limited to 64 KiB, 128 nested operands and 4096-bit
// numeric literals; literals are never truncated to a signal's width.
BoolExpr* compileC2RtlExpression(
    std::string_view text, const ExpressionSignalResolver& resolve);

// Checks the same syntax without resolving signals or allocating BoolExprs.
void validateC2RtlExpression(std::string_view text);

}  // namespace KEPLER_FORMAL::C2RTL
