// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "property/SvaPropertyTranslator.h"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "NLDB0.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLSVConstructor.h"
#include "property/BoolExprOps.h"

#include "slang/ast/ASTContext.h"
#include "slang/ast/Compilation.h"
#include "slang/ast/EvalContext.h"
#include "slang/ast/Expression.h"
#include "slang/ast/Scope.h"
#include "slang/ast/Statement.h"
#include "slang/ast/TimingControl.h"
#include "slang/ast/expressions/AssertionExpr.h"
#include "slang/ast/expressions/CallExpression.h"
#include "slang/ast/expressions/ConversionExpression.h"
#include "slang/ast/expressions/LiteralExpressions.h"
#include "slang/ast/expressions/MiscExpressions.h"
#include "slang/ast/expressions/OperatorExpressions.h"
#include "slang/ast/expressions/SelectExpressions.h"
#include "slang/ast/statements/MiscStatements.h"
#include "slang/ast/symbols/BlockSymbols.h"
#include "slang/ast/symbols/MemberSymbols.h"
#include "slang/ast/symbols/VariableSymbols.h"
#include "slang/ast/types/Type.h"
#include "slang/numeric/ConstantValue.h"
#include "slang/numeric/SVInt.h"

namespace KEPLER_FORMAL::SEC {

namespace {

using namespace PropertyLogic;
using slang::ast::AssertionExpr;
using slang::ast::AssertionExprKind;
using slang::ast::BinaryOperator;
using slang::ast::Expression;
using slang::ast::ExpressionKind;
using slang::ast::SymbolKind;
using slang::ast::UnaryOperator;

// A packed value in the current step, least significant bit first.
using Bits = std::vector<BoolExpr*>;

struct Unsupported {
  std::string reason;
};

// Bounds that keep one property from allocating an unreasonable amount of
// logic or history.
constexpr size_t kMaxWidth = 4096;
constexpr int64_t kMaxPastTicks = 4096;

BoolExpr* orReduce(const Bits& bits) {
  BoolExpr* result = constFalse();
  for (auto* bit : bits) {
    result = orExpr(result, bit);
  }
  return result;
}

BoolExpr* andReduce(const Bits& bits) {
  BoolExpr* result = constTrue();
  for (auto* bit : bits) {
    result = andExpr(result, bit);
  }
  return result;
}

BoolExpr* xorReduce(const Bits& bits) {
  BoolExpr* result = constFalse();
  for (auto* bit : bits) {
    result = xorExpr(result, bit);
  }
  return result;
}

BoolExpr* muxBit(BoolExpr* select, BoolExpr* whenTrue, BoolExpr* whenFalse) {
  if (whenTrue == whenFalse) {
    return whenTrue;
  }
  return orExpr(andExpr(select, whenTrue), andExpr(notExpr(select), whenFalse));
}

// Truncates, or extends with the sign bit or with zeros.
Bits resize(Bits bits, size_t width, bool signExtend) {
  BoolExpr* fill = signExtend && !bits.empty() ? bits.back() : constFalse();
  bits.resize(width, fill);
  return bits;
}

Bits addBits(const Bits& a, const Bits& b, BoolExpr* carry) {
  Bits sum(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    BoolExpr* partial = xorExpr(a[i], b[i]);
    sum[i] = xorExpr(partial, carry);
    carry = orExpr(andExpr(a[i], b[i]), andExpr(carry, partial));
  }
  return sum;
}

Bits invert(Bits bits) {
  for (auto& bit : bits) {
    bit = notExpr(bit);
  }
  return bits;
}

BoolExpr* equalBits(const Bits& a, const Bits& b) {
  BoolExpr* result = constTrue();
  for (size_t i = 0; i < a.size(); ++i) {
    result = andExpr(result, notExpr(xorExpr(a[i], b[i])));
  }
  return result;
}

// a < b. For signed operands the sign bits swap roles: a negative value is
// below a non-negative one.
BoolExpr* lessThan(const Bits& a, const Bits& b, bool isSigned) {
  BoolExpr* less = constFalse();
  for (size_t i = 0; i < a.size(); ++i) {
    BoolExpr* ai = a[i];
    BoolExpr* bi = b[i];
    if (isSigned && i + 1 == a.size()) {
      std::swap(ai, bi);
    }
    less = orExpr(andExpr(notExpr(ai), bi),
                  andExpr(notExpr(xorExpr(ai, bi)), less));
  }
  return less;
}

Bits shiftBy(const Bits& value, size_t distance, bool left, BoolExpr* fill) {
  Bits shifted(value.size(), fill);
  for (size_t i = 0; i < value.size(); ++i) {
    if (left) {
      if (i >= distance) {
        shifted[i] = value[i - distance];
      }
    } else if (i + distance < value.size()) {
      shifted[i] = value[i + distance];
    }
  }
  return shifted;
}

// Compiles the Boolean and packed operands of one assertion into logic over
// the signals the frontend bound for it. Widths, signs and implicit
// conversions are the ones slang elaborated; nothing is re-derived here.
class ExpressionCompiler {
 public:
  ExpressionCompiler(const slang::ast::Scope& scope,
                     std::string bodyPrefix,
                     std::unordered_map<std::string, Bits> signals,
                     const std::function<BoolExpr*(BoolExpr*)>& previous)
      : astContext_(scope, slang::ast::LookupLocation::max),
        bodyPrefix_(std::move(bodyPrefix)),
        signals_(std::move(signals)),
        previous_(previous) {}

  Bits compile(const Expression& expr) {
    if (expr.kind == ExpressionKind::MinTypMax) {
      return compile(expr.as<slang::ast::MinTypMaxExpression>().selected());
    }
    if (!expr.type->isIntegral()) {
      throw Unsupported{"an operand that is not a packed value"};
    }
    const size_t width = expr.type->getBitWidth();
    if (width == 0 || width > kMaxWidth) {
      throw Unsupported{"an operand of unsupported width"};
    }
    if (const auto value = constant(expr)) {
      return fromConstant(*value, width);
    }
    Bits bits = compileNonConstant(expr, width);
    if (bits.size() != width) {
      throw Unsupported{"an operand whose width could not be resolved"};  // LCOV_EXCL_LINE
    }
    return bits;
  }

  // A Boolean operand is true when its value is not zero.
  BoolExpr* truth(const Expression& expr) { return orReduce(compile(expr)); }

 private:
  // The value of an expression that does not depend on any signal.
  std::optional<slang::SVInt> constant(const Expression& expr) const {
    slang::ast::EvalContext context(astContext_);
    const auto value = expr.eval(context);
    if (!value.isInteger()) {
      return std::nullopt;
    }
    return value.integer();
  }

  std::optional<int64_t> constantInt(const Expression& expr) const {
    const auto value = constant(expr);
    if (!value || value->hasUnknown()) {
      return std::nullopt;
    }
    return value->as<int64_t>();
  }

  static Bits fromConstant(const slang::SVInt& value, size_t width) {
    if (value.hasUnknown()) {
      throw Unsupported{"a constant with unknown bits"};
    }
    Bits bits;
    bits.reserve(width);
    const size_t valueWidth = value.getBitWidth();
    for (size_t i = 0; i < width && i < valueWidth; ++i) {
      bits.push_back(value[static_cast<int32_t>(i)].value == 1 ? constTrue()
                                                               : constFalse());
    }
    return resize(std::move(bits), width, value.isSigned());
  }

  Bits compileNonConstant(const Expression& expr, size_t width) {
    switch (expr.kind) {
      case ExpressionKind::NamedValue:
        return signal(expr, width);
      case ExpressionKind::Conversion: {
        const auto& conversion = expr.as<slang::ast::ConversionExpression>();
        const auto& operand = conversion.operand();
        if (!operand.type->isIntegral()) {
          throw Unsupported{"a conversion from a value that is not packed"};
        }
        // A propagated conversion extends by the sign of the type it
        // propagates; every other conversion by the sign of its operand.
        const bool signExtend =
            conversion.conversionKind == slang::ast::ConversionKind::Propagated
                ? expr.type->isSigned()
                : operand.type->isSigned();
        return resize(compile(operand), width, signExtend);
      }
      case ExpressionKind::UnaryOp:
        return unary(expr.as<slang::ast::UnaryExpression>(), width);
      case ExpressionKind::BinaryOp:
        return binary(expr.as<slang::ast::BinaryExpression>(), width);
      case ExpressionKind::ConditionalOp: {
        const auto& conditional = expr.as<slang::ast::ConditionalExpression>();
        if (conditional.conditions.size() != 1 ||
            conditional.conditions.front().pattern != nullptr) {
          throw Unsupported{"a conditional operator with a pattern"};
        }
        BoolExpr* select = truth(*conditional.conditions.front().expr);
        const Bits whenTrue = compile(conditional.left());
        const Bits whenFalse = compile(conditional.right());
        Bits result(width);
        for (size_t i = 0; i < width; ++i) {
          result[i] = muxBit(select, whenTrue[i], whenFalse[i]);
        }
        return result;
      }
      case ExpressionKind::Inside:
        return resize({inside(expr.as<slang::ast::InsideExpression>())}, width,
                      false);
      case ExpressionKind::Concatenation: {
        // Operands are written most significant first.
        const auto operands =
            expr.as<slang::ast::ConcatenationExpression>().operands();
        Bits result;
        for (auto it = operands.rbegin(); it != operands.rend(); ++it) {
          if ((*it)->type->isVoid()) {
            continue;
          }
          const Bits part = compile(**it);
          result.insert(result.end(), part.begin(), part.end());
        }
        return result;
      }
      case ExpressionKind::Replication: {
        const auto& replication = expr.as<slang::ast::ReplicationExpression>();
        const auto count = constantInt(replication.count());
        if (!count || *count < 0) {
          throw Unsupported{"a replication with a variable count"};
        }
        const Bits part = compile(replication.concat());
        Bits result;
        for (int64_t i = 0; i < *count; ++i) {
          result.insert(result.end(), part.begin(), part.end());
        }
        return result;
      }
      case ExpressionKind::ElementSelect:
        return elementSelect(expr.as<slang::ast::ElementSelectExpression>(),
                             width);
      case ExpressionKind::RangeSelect:
        return rangeSelect(expr.as<slang::ast::RangeSelectExpression>(), width);
      case ExpressionKind::MemberAccess: {
        const auto& access = expr.as<slang::ast::MemberAccessExpression>();
        if (access.member.kind != SymbolKind::Field ||
            !access.value().type->isIntegral()) {
          throw Unsupported{"a member access outside a packed structure"};
        }
        const Bits value = compile(access.value());
        const auto offset = static_cast<size_t>(
            access.member.as<slang::ast::FieldSymbol>().bitOffset);
        if (offset + width > value.size()) {
          throw Unsupported{"a member access outside its structure"};  // LCOV_EXCL_LINE
        }
        return Bits(value.begin() + offset, value.begin() + offset + width);
      }
      case ExpressionKind::Call:
        return call(expr.as<slang::ast::CallExpression>(), width);
      case ExpressionKind::HierarchicalValue:
        throw Unsupported{"a hierarchical reference"};
      default:
        throw Unsupported{"an expression form that is not supported yet"};
    }
  }

  Bits signal(const Expression& expr, size_t width) const {
    const auto& symbol = expr.as<slang::ast::ValueExpressionBase>().symbol;
    // The frontend names a bound signal by its path below the module.
    std::string name = symbol.getHierarchicalPath();
    if (name.starts_with(bodyPrefix_)) {
      name.erase(0, bodyPrefix_.size());
    }
    if (name.empty()) {
      name = std::string(symbol.name);  // LCOV_EXCL_LINE
    }
    const auto found = signals_.find(name);
    if (found == signals_.end() || found->second.size() != width) {
      throw Unsupported{"signal `" + name + "` is not bound by the frontend"};
    }
    return found->second;
  }

  Bits unary(const slang::ast::UnaryExpression& expr, size_t width) {
    const Bits operand = compile(expr.operand());
    switch (expr.op) {
      case UnaryOperator::Plus:
        return operand;
      case UnaryOperator::Minus:
        return addBits(invert(operand), Bits(operand.size(), constFalse()),
                       constTrue());
      case UnaryOperator::BitwiseNot:
        return invert(operand);
      case UnaryOperator::LogicalNot:
        return resize({notExpr(orReduce(operand))}, width, false);
      case UnaryOperator::BitwiseAnd:
        return resize({andReduce(operand)}, width, false);
      case UnaryOperator::BitwiseNand:
        return resize({notExpr(andReduce(operand))}, width, false);
      case UnaryOperator::BitwiseOr:
        return resize({orReduce(operand)}, width, false);
      case UnaryOperator::BitwiseNor:
        return resize({notExpr(orReduce(operand))}, width, false);
      case UnaryOperator::BitwiseXor:
        return resize({xorReduce(operand)}, width, false);
      case UnaryOperator::BitwiseXnor:
        return resize({notExpr(xorReduce(operand))}, width, false);
      default:
        throw Unsupported{"an increment or decrement operator"};
    }
  }

  Bits binary(const slang::ast::BinaryExpression& expr, size_t width) {
    const auto bit = [&](BoolExpr* value) { return resize({value}, width, false); };
    switch (expr.op) {
      case BinaryOperator::LogicalAnd:
        return bit(andExpr(truth(expr.left()), truth(expr.right())));
      case BinaryOperator::LogicalOr:
        return bit(orExpr(truth(expr.left()), truth(expr.right())));
      case BinaryOperator::LogicalImplication:
        return bit(orExpr(notExpr(truth(expr.left())), truth(expr.right())));
      case BinaryOperator::LogicalEquivalence:
        return bit(notExpr(xorExpr(truth(expr.left()), truth(expr.right()))));
      case BinaryOperator::LogicalShiftLeft:
      case BinaryOperator::ArithmeticShiftLeft:
        return shift(expr, true, constFalse());
      case BinaryOperator::LogicalShiftRight:
        return shift(expr, false, constFalse());
      case BinaryOperator::ArithmeticShiftRight: {
        const Bits value = compile(expr.left());
        return shift(expr, false,
                     expr.type->isSigned() ? value.back() : constFalse());
      }
      default:
        break;
    }

    // Every remaining operator works on two operands of one common width,
    // signed only when both are.
    const bool isSigned =
        expr.left().type->isSigned() && expr.right().type->isSigned();
    Bits left = compile(expr.left());
    Bits right = compile(expr.right());
    const size_t common = std::max(left.size(), right.size());
    left = resize(std::move(left), common, isSigned);
    right = resize(std::move(right), common, isSigned);
    switch (expr.op) {
      case BinaryOperator::Add:
        return resize(addBits(left, right, constFalse()), width, isSigned);
      case BinaryOperator::Subtract:
        return resize(addBits(left, invert(right), constTrue()), width, isSigned);
      case BinaryOperator::BinaryAnd:
      case BinaryOperator::BinaryOr:
      case BinaryOperator::BinaryXor:
      case BinaryOperator::BinaryXnor: {
        Bits result(common);
        for (size_t i = 0; i < common; ++i) {
          switch (expr.op) {
            case BinaryOperator::BinaryAnd:
              result[i] = andExpr(left[i], right[i]);
              break;
            case BinaryOperator::BinaryOr:
              result[i] = orExpr(left[i], right[i]);
              break;
            case BinaryOperator::BinaryXor:
              result[i] = xorExpr(left[i], right[i]);
              break;
            default:
              result[i] = notExpr(xorExpr(left[i], right[i]));
              break;
          }
        }
        return resize(std::move(result), width, isSigned);
      }
      // Values are two-valued, so case equality is plain equality.
      case BinaryOperator::Equality:
      case BinaryOperator::CaseEquality:
        return bit(equalBits(left, right));
      case BinaryOperator::Inequality:
      case BinaryOperator::CaseInequality:
        return bit(notExpr(equalBits(left, right)));
      case BinaryOperator::LessThan:
        return bit(lessThan(left, right, isSigned));
      case BinaryOperator::GreaterThan:
        return bit(lessThan(right, left, isSigned));
      case BinaryOperator::LessThanEqual:
        return bit(notExpr(lessThan(right, left, isSigned)));
      case BinaryOperator::GreaterThanEqual:
        return bit(notExpr(lessThan(left, right, isSigned)));
      default:
        throw Unsupported{"multiplication, division, modulo, power or a "
                          "wildcard comparison of signals"};
    }
  }

  Bits shift(const slang::ast::BinaryExpression& expr, bool left, BoolExpr* fill) {
    Bits value = compile(expr.left());
    if (const auto amount = constantInt(expr.right())) {
      if (*amount < 0) {
        throw Unsupported{"a negative shift amount"};  // LCOV_EXCL_LINE
      }
      const auto distance = static_cast<uint64_t>(*amount);
      return distance >= value.size()
                 ? Bits(value.size(), fill)
                 : shiftBy(value, static_cast<size_t>(distance), left, fill);
    }
    // One stage per bit of the amount: stage k shifts by 2^k when that bit is
    // set, or shifts everything out once 2^k reaches the width.
    const Bits amount = compile(expr.right());
    for (size_t stage = 0; stage < amount.size(); ++stage) {
      const bool beyond = stage >= 63 || (uint64_t{1} << stage) >= value.size();
      const Bits shifted =
          beyond ? Bits(value.size(), fill)
                 : shiftBy(value, size_t{1} << stage, left, fill);
      for (size_t i = 0; i < value.size(); ++i) {
        value[i] = muxBit(amount[stage], shifted[i], value[i]);
      }
    }
    return value;
  }

  BoolExpr* inside(const slang::ast::InsideExpression& expr) {
    const bool leftSigned = expr.left().type->isSigned();
    const Bits left = compile(expr.left());
    // Each comparison is made at the common width of its two operands, signed
    // only when both are.
    struct Pair {
      Bits value;
      Bits other;
      bool isSigned;
    };
    const auto widen = [&](const Expression& other) {
      const Bits bits = compile(other);
      const size_t common = std::max(left.size(), bits.size());
      const bool isSigned = leftSigned && other.type->isSigned();
      return Pair{resize(left, common, isSigned), resize(bits, common, isSigned),
                  isSigned};
    };
    BoolExpr* result = constFalse();
    for (const auto* item : expr.rangeList()) {
      if (item->kind != ExpressionKind::ValueRange) {
        const Pair pair = widen(*item);
        result = orExpr(result, equalBits(pair.value, pair.other));
        continue;
      }
      const auto& range = item->as<slang::ast::ValueRangeExpression>();
      if (range.rangeKind != slang::ast::ValueRangeKind::Simple) {
        throw Unsupported{"a tolerance range"};
      }
      const Pair low = widen(range.left());
      const Pair high = widen(range.right());
      result = orExpr(
          result,
          andExpr(notExpr(lessThan(low.value, low.other, low.isSigned)),
                  notExpr(lessThan(high.other, high.value, high.isSigned))));
    }
    return result;
  }

  // Offset of element `index` of a packed value, in elements from its least
  // significant end.
  static std::optional<size_t> elementOffset(const slang::ConstantRange& range,
                                             int64_t index) {
    if (index < range.lower() || index > range.upper()) {
      return std::nullopt;
    }
    return static_cast<size_t>(range.translateIndex(static_cast<int32_t>(index)));
  }

  Bits elementSelect(const slang::ast::ElementSelectExpression& expr,
                     size_t width) {
    const auto& valueType = *expr.value().type;
    if (!valueType.isIntegral() || !valueType.hasFixedRange()) {
      throw Unsupported{"an element select on a value that is not packed"};
    }
    const Bits value = compile(expr.value());
    const auto range = valueType.getFixedRange();
    // A read outside the vector is unknown: it is reported, not given a value.
    const auto element = [&](int64_t index) {
      const auto offset = elementOffset(range, index);
      if (!offset || (*offset + 1) * width > value.size()) {
        throw Unsupported{"an index outside its vector"};
      }
      return Bits(value.begin() + *offset * width,
                  value.begin() + (*offset + 1) * width);
    };
    if (const auto index = constantInt(expr.selector())) {
      return element(*index);
    }
    const Bits selector = compile(expr.selector());
    const bool selectorSigned = expr.selector().type->isSigned();
    // Every value the index can take must name an element.
    const size_t indexBits = selector.size();
    if (indexBits > 30 ||
        (selectorSigned ? -(int64_t{1} << (indexBits - 1)) : 0) < range.lower() ||
        (int64_t{1} << (selectorSigned ? indexBits - 1 : indexBits)) - 1 >
            range.upper()) {
      throw Unsupported{"an index that can fall outside its vector"};
    }
    Bits result(width, constFalse());
    for (int64_t index = range.lower(); index <= range.upper(); ++index) {
      BoolExpr* selected = equalsConstant(selector, index, selectorSigned);
      if (selected == constFalse()) {
        continue;
      }
      const Bits bits = element(index);
      for (size_t i = 0; i < width; ++i) {
        result[i] = orExpr(result[i], andExpr(selected, bits[i]));
      }
    }
    return result;
  }

  static BoolExpr* equalsConstant(const Bits& bits, int64_t value, bool isSigned) {
    const size_t width = bits.size();
    if (width < 63) {
      const int64_t low = isSigned ? -(int64_t{1} << (width - 1)) : 0;
      const int64_t high = isSigned ? (int64_t{1} << (width - 1)) - 1
                                    : (int64_t{1} << width) - 1;
      if (value < low || value > high) {
        return constFalse();
      }
    } else if (!isSigned && value < 0) {
      return constFalse();  // LCOV_EXCL_LINE
    }
    BoolExpr* result = constTrue();
    for (size_t i = 0; i < width; ++i) {
      const bool one = i < 63 ? ((value >> i) & 1) != 0 : value < 0;
      result = andExpr(result, one ? bits[i] : notExpr(bits[i]));
    }
    return result;
  }

  Bits rangeSelect(const slang::ast::RangeSelectExpression& expr, size_t width) {
    const auto& valueType = *expr.value().type;
    if (!valueType.isIntegral() || !valueType.hasFixedRange()) {
      throw Unsupported{"a part select on a value that is not packed"};
    }
    const auto first = constantInt(expr.left());
    const auto second = constantInt(expr.right());
    if (!first || !second) {
      throw Unsupported{"a part select with a variable bound"};
    }
    int64_t last = *second;
    if (expr.getSelectionKind() == slang::ast::RangeSelectionKind::IndexedUp) {
      last = *first + *second - 1;
    } else if (expr.getSelectionKind() ==
               slang::ast::RangeSelectionKind::IndexedDown) {
      last = *first - *second + 1;
    }
    const Bits value = compile(expr.value());
    const auto range = valueType.getFixedRange();
    const size_t elementWidth = value.size() / range.width();
    const auto a = elementOffset(range, *first);
    const auto b = elementOffset(range, last);
    if (!a || !b || std::min(*a, *b) * elementWidth + width > value.size()) {
      throw Unsupported{"a part select outside its vector"};
    }
    const size_t low = std::min(*a, *b) * elementWidth;
    return Bits(value.begin() + low, value.begin() + low + width);
  }

  // The value `bit` had at the previous tick of the property clock.
  BoolExpr* previous(BoolExpr* bit) {
    const auto found = previousByBit_.find(bit);
    if (found != previousByBit_.end()) {
      return found->second;
    }
    return previousByBit_.emplace(bit, previous_(bit)).first->second;
  }

  Bits call(const slang::ast::CallExpression& expr, size_t width) {
    const std::string name(expr.getSubroutineName());
    if (!expr.isSystemCall()) {
      throw Unsupported{"a call to function `" + name + "`"};
    }
    const auto arguments = expr.arguments();
    const auto present = [&](size_t index) {
      return index < arguments.size() && arguments[index] != nullptr &&
             arguments[index]->kind != ExpressionKind::EmptyArgument;
    };
    const auto bit = [&](BoolExpr* value) { return resize({value}, width, false); };
    if (!present(0)) {
      throw Unsupported{"`" + name + "` without an argument"};
    }
    const Expression& argument = *arguments[0];

    if (name == "$past") {
      if (present(2) || present(3)) {
        throw Unsupported{"$past with a gating expression or its own clock"};
      }
      int64_t ticks = 1;
      if (present(1)) {
        const auto count = constantInt(*arguments[1]);
        if (!count || *count < 1 || *count > kMaxPastTicks) {
          throw Unsupported{"$past with an unsupported number of ticks"};
        }
        ticks = *count;
      }
      Bits bits = compile(argument);
      for (int64_t tick = 0; tick < ticks; ++tick) {
        for (auto& value : bits) {
          value = previous(value);
        }
      }
      return resize(std::move(bits), width, argument.type->isSigned());
    }
    if (name == "$rose" || name == "$fell" || name == "$stable" ||
        name == "$changed") {
      if (present(1)) {
        throw Unsupported{"`" + name + "` with its own clock"};
      }
      const Bits now = compile(argument);
      // $rose and $fell look at the least significant bit only.
      if (name == "$rose") {
        return bit(andExpr(now.front(), notExpr(previous(now.front()))));
      }
      if (name == "$fell") {
        return bit(andExpr(notExpr(now.front()), previous(now.front())));
      }
      Bits before = now;
      for (auto& value : before) {
        value = previous(value);
      }
      BoolExpr* stable = equalBits(now, before);
      return bit(name == "$stable" ? stable : notExpr(stable));
    }
    if (name == "$onehot" || name == "$onehot0") {
      BoolExpr* one = constFalse();
      BoolExpr* several = constFalse();
      for (auto* value : compile(argument)) {
        several = orExpr(several, andExpr(one, value));
        one = orExpr(one, value);
      }
      return bit(name == "$onehot" ? andExpr(one, notExpr(several))
                                   : notExpr(several));
    }
    if (name == "$countones") {
      Bits count(width, constFalse());
      for (auto* value : compile(argument)) {
        Bits addend(width, constFalse());
        addend.front() = value;
        count = addBits(count, addend, constFalse());
      }
      return count;
    }
    if (name == "$isunknown") {
      // Values are two-valued: nothing is ever unknown.
      return bit(constFalse());
    }
    // In a concurrent assertion every value is already the sampled one.
    if (name == "$signed" || name == "$unsigned" || name == "$sampled") {
      return resize(compile(argument), width, argument.type->isSigned());
    }
    throw Unsupported{"the system function `" + name + "`"};
  }

  slang::ast::ASTContext astContext_;
  std::string bodyPrefix_;
  std::unordered_map<std::string, Bits> signals_;
  const std::function<BoolExpr*(BoolExpr*)>& previous_;
  std::unordered_map<BoolExpr*, BoolExpr*> previousByBit_;
};

std::string formatRange(const slang::ast::SequenceRange& range) {
  return std::to_string(range.min) + " " +
         (range.max ? std::to_string(*range.max) : std::string("$"));
}

std::string wrapRepetition(
    const std::optional<slang::ast::SequenceRepetition>& repetition,
    const std::string& operand) {
  if (!repetition) {
    return operand;
  }
  const char* kind = "consecutive";
  if (repetition->kind == slang::ast::SequenceRepetition::Nonconsecutive) {
    kind = "nonconsecutive";
  } else if (repetition->kind == slang::ast::SequenceRepetition::GoTo) {
    kind = "goto";
  }
  return std::string("(repeat ") + kind + " " + formatRange(repetition->range) +
         " " + operand + ")";
}

// Writes the temporal structure of a property as the prefix expression
// PropertyCompiler reads, compiling each Boolean operand into one leaf.
class Transcriber {
 public:
  explicit Transcriber(ExpressionCompiler& compiler) : compiler_(compiler) {}

  std::vector<BoolExpr*> leaves;

  std::string leaf(const Expression& expr) {
    BoolExpr* value = compiler_.truth(expr);
    const auto found = std::find(leaves.begin(), leaves.end(), value);
    const auto index = static_cast<size_t>(found - leaves.begin());
    if (found == leaves.end()) {
      leaves.push_back(value);
    }
    return "(leaf " + std::to_string(index) + ")";
  }

  std::string transcribe(const AssertionExpr& expr) {
    switch (expr.kind) {
      case AssertionExprKind::Clocking:
        throw Unsupported{"a property with more than one clocking event"};
      case AssertionExprKind::DisableIff:
        throw Unsupported{"a nested disable iff"};
      case AssertionExprKind::Simple: {
        const auto& simple = expr.as<slang::ast::SimpleAssertionExpr>();
        if (simple.expr.kind != ExpressionKind::AssertionInstance) {
          return wrapRepetition(simple.repetition, leaf(simple.expr));
        }
        // A named sequence or property: slang already substituted the actual
        // arguments, so its body is transcribed in place.
        const auto& instance =
            simple.expr.as<slang::ast::AssertionInstanceExpression>();
        if (instance.isRecursiveProperty) {
          throw Unsupported{"a recursive property"};
        }
        if (!instance.localVars.empty()) {
          throw Unsupported{"local variables in sequences and properties"};
        }
        return wrapRepetition(simple.repetition, transcribe(instance.body));
      }
      case AssertionExprKind::SequenceConcat: {
        std::string result = "(concat";
        for (const auto& element :
             expr.as<slang::ast::SequenceConcatExpr>().elements) {
          const std::string sequence = transcribe(*element.sequence);
          result += " (" + formatRange(element.delay) + " " + sequence + ")";
        }
        return result + ")";
      }
      case AssertionExprKind::SequenceWithMatch: {
        const auto& withMatch = expr.as<slang::ast::SequenceWithMatchExpr>();
        if (!withMatch.matchItems.empty()) {
          throw Unsupported{"local variable assignments in sequences"};
        }
        return wrapRepetition(withMatch.repetition, transcribe(withMatch.expr));
      }
      case AssertionExprKind::Unary: {
        const auto& unary = expr.as<slang::ast::UnaryAssertionExpr>();
        const char* name = "not";
        switch (unary.op) {
          case slang::ast::UnaryAssertionOperator::Not: name = "not"; break;
          case slang::ast::UnaryAssertionOperator::NextTime: name = "nexttime"; break;
          case slang::ast::UnaryAssertionOperator::SNextTime: name = "s_nexttime"; break;
          case slang::ast::UnaryAssertionOperator::Always: name = "always"; break;
          case slang::ast::UnaryAssertionOperator::SAlways: name = "s_always"; break;
          case slang::ast::UnaryAssertionOperator::Eventually: name = "eventually"; break;
          case slang::ast::UnaryAssertionOperator::SEventually: name = "s_eventually"; break;
        }
        std::string result = std::string("(") + name;
        if (unary.range) {
          result += " " + formatRange(*unary.range);
        }
        return result + " " + transcribe(unary.expr) + ")";
      }
      case AssertionExprKind::Binary: {
        const auto& binary = expr.as<slang::ast::BinaryAssertionExpr>();
        const char* name = "and";
        switch (binary.op) {
          case slang::ast::BinaryAssertionOperator::And: name = "and"; break;
          case slang::ast::BinaryAssertionOperator::Or: name = "or"; break;
          case slang::ast::BinaryAssertionOperator::Intersect: name = "intersect"; break;
          case slang::ast::BinaryAssertionOperator::Throughout: name = "throughout"; break;
          case slang::ast::BinaryAssertionOperator::Within: name = "within"; break;
          case slang::ast::BinaryAssertionOperator::Iff: name = "iff"; break;
          case slang::ast::BinaryAssertionOperator::Until: name = "until"; break;
          case slang::ast::BinaryAssertionOperator::SUntil: name = "s_until"; break;
          case slang::ast::BinaryAssertionOperator::UntilWith: name = "until_with"; break;
          case slang::ast::BinaryAssertionOperator::SUntilWith: name = "s_until_with"; break;
          case slang::ast::BinaryAssertionOperator::Implies: name = "property_implies"; break;
          case slang::ast::BinaryAssertionOperator::OverlappedImplication: name = "implies"; break;
          case slang::ast::BinaryAssertionOperator::NonOverlappedImplication: name = "implies_next"; break;
          case slang::ast::BinaryAssertionOperator::OverlappedFollowedBy: name = "followed_by"; break;
          case slang::ast::BinaryAssertionOperator::NonOverlappedFollowedBy: name = "followed_by_next"; break;
        }
        // Leaves are numbered in source order: left before right.
        const std::string left = transcribe(binary.left);
        const std::string right = transcribe(binary.right);
        return std::string("(") + name + " " + left + " " + right + ")";
      }
      case AssertionExprKind::FirstMatch: {
        const auto& firstMatch = expr.as<slang::ast::FirstMatchAssertionExpr>();
        if (!firstMatch.matchItems.empty()) {
          throw Unsupported{"local variable assignments in sequences"};
        }
        return "(first_match " + transcribe(firstMatch.seq) + ")";
      }
      case AssertionExprKind::StrongWeak: {
        const auto& strongWeak = expr.as<slang::ast::StrongWeakAssertionExpr>();
        return std::string(strongWeak.strength ==
                                   slang::ast::StrongWeakAssertionExpr::Strong
                               ? "(strong "
                               : "(weak ") +
               transcribe(strongWeak.expr) + ")";
      }
      case AssertionExprKind::Abort: {
        const auto& abort = expr.as<slang::ast::AbortAssertionExpr>();
        std::string name = abort.isSync ? "sync_" : "";
        name += abort.action == slang::ast::AbortAssertionExpr::Accept
                    ? "accept_on"
                    : "reject_on";
        const std::string condition = leaf(abort.condition);
        return "(" + name + " " + condition + " " + transcribe(abort.expr) + ")";
      }
      case AssertionExprKind::Conditional: {
        const auto& conditional = expr.as<slang::ast::ConditionalAssertionExpr>();
        std::string result = "(if " + leaf(conditional.condition);
        result += " " + transcribe(conditional.ifExpr);
        if (conditional.elseExpr) {
          result += " " + transcribe(*conditional.elseExpr);
        }
        return result + ")";
      }
      case AssertionExprKind::Case:
        throw Unsupported{"a case property"};
      case AssertionExprKind::Invalid:
        break;
    }
    throw Unsupported{"an invalid property expression"};  // LCOV_EXCL_LINE
  }

 private:
  ExpressionCompiler& compiler_;
};

// The assertion statement of a module-level assertion item.
const slang::ast::ConcurrentAssertionStatement* findAssertion(
    const slang::ast::Statement& body) {
  const slang::ast::Statement* current = &body;
  for (;;) {
    if (current->kind == slang::ast::StatementKind::Block) {
      current = &current->as<slang::ast::BlockStatement>().body;
    } else if (current->kind == slang::ast::StatementKind::List &&
               current->as<slang::ast::StatementList>().list.size() == 1) {
      current = current->as<slang::ast::StatementList>().list.front();
    } else {
      break;
    }
  }
  return current->kind == slang::ast::StatementKind::ConcurrentAssertion
             ? &current->as<slang::ast::ConcurrentAssertionStatement>()
             : nullptr;
}

TranslatedProperty translate(const PropertyTranslationRequest& request) {
  TranslatedProperty result;
  try {
    const auto* instance = request.instance;
    const auto* link = naja::NL::SNLSVLiveASTLinkRegistry::get(
        instance->getDesign()->getDB());
    const auto* symbol = link != nullptr ? link->getSymbol(instance) : nullptr;
    if (symbol == nullptr || symbol->kind != SymbolKind::ProceduralBlock) {
      throw Unsupported{"the parsed assertion is not available: the design "
                        "must be loaded with keepASTLink"};
    }
    const auto& block = symbol->as<slang::ast::ProceduralBlockSymbol>();
    const auto* assertion = findAssertion(block.getBody());
    const auto* scope = block.getParentScope();
    if (assertion == nullptr || scope == nullptr) {
      throw Unsupported{"no concurrent assertion was found for the instance"};  // LCOV_EXCL_LINE
    }

    // Signal names are paths below the module that holds the assertion.
    std::string bodyPrefix;
    for (const auto* outer = scope; outer != nullptr;
         outer = outer->asSymbol().getParentScope()) {
      if (outer->asSymbol().kind == SymbolKind::InstanceBody) {
        bodyPrefix = outer->asSymbol().getHierarchicalPath() + ".";
        break;
      }
    }
    std::unordered_map<std::string, Bits> signals;
    for (const auto& binding :
         naja::NL::NLDB0::getFormalProperty(instance).bindings) {
      if (binding.offset + binding.width > request.data.size()) {
        throw Unsupported{"a signal binding outside the captured data"};  // LCOV_EXCL_LINE
      }
      signals.emplace(binding.name,
                      Bits(request.data.begin() + binding.offset,
                           request.data.begin() + binding.offset + binding.width));
    }
    ExpressionCompiler compiler(*scope, std::move(bodyPrefix), std::move(signals),
                                request.previous);

    // Peel what only selects the clock and the disable condition, through
    // named properties.
    const AssertionExpr* body = &assertion->propertySpec;
    const slang::ast::TimingControl* clocking = nullptr;
    const Expression* disable = nullptr;
    for (;;) {
      if (body->kind == AssertionExprKind::Clocking && clocking == nullptr) {
        const auto& clocked = body->as<slang::ast::ClockingAssertionExpr>();
        clocking = &clocked.clocking;
        body = &clocked.expr;
      } else if (body->kind == AssertionExprKind::DisableIff &&
                 disable == nullptr) {
        const auto& disabled = body->as<slang::ast::DisableIffAssertionExpr>();
        disable = &disabled.condition;
        body = &disabled.expr;
      } else if (body->kind == AssertionExprKind::Simple &&
                 !body->as<slang::ast::SimpleAssertionExpr>().repetition &&
                 body->as<slang::ast::SimpleAssertionExpr>().expr.kind ==
                     ExpressionKind::AssertionInstance) {
        const auto& named = body->as<slang::ast::SimpleAssertionExpr>()
                                .expr.as<slang::ast::AssertionInstanceExpression>();
        if (named.isRecursiveProperty) {
          throw Unsupported{"a recursive property"};
        }
        if (!named.localVars.empty()) {
          throw Unsupported{"local variables in sequences and properties"};
        }
        body = &named.body;
      } else {
        break;
      }
    }

    auto& compilation = scope->getCompilation();
    if (clocking == nullptr) {
      if (const auto* inherited = compilation.getDefaultClocking(*scope);
          inherited != nullptr && inherited->kind == SymbolKind::ClockingBlock) {
        clocking =
            &inherited->as<slang::ast::ClockingBlockSymbol>().getEvent();
      }
    }
    if (clocking == nullptr) {
      throw Unsupported{"the property has no clocking event"};
    }
    if (clocking->kind != slang::ast::TimingControlKind::SignalEvent) {
      throw Unsupported{"a clocking event that is not one posedge or negedge"};
    }
    const auto& event = clocking->as<slang::ast::SignalEventControl>();
    if (event.iffCondition != nullptr) {
      throw Unsupported{"a clocking event with an iff condition"};
    }
    if (event.edge != slang::ast::EdgeKind::PosEdge &&
        event.edge != slang::ast::EdgeKind::NegEdge) {
      throw Unsupported{"a clocking event that is not one posedge or negedge"};
    }
    const Bits clock = compiler.compile(event.expr);
    if (clock.size() != 1) {
      throw Unsupported{"a property clock that is not a single bit"};
    }
    // The property advances on the rising edge of this signal.
    result.clock = event.edge == slang::ast::EdgeKind::PosEdge
                       ? clock.front()
                       : notExpr(clock.front());

    if (disable == nullptr) {
      disable = compilation.getDefaultDisable(*scope);
    }
    result.disable = disable != nullptr ? compiler.truth(*disable) : constFalse();

    Transcriber transcriber(compiler);
    result.expression = transcriber.transcribe(*body);
    result.leaves = std::move(transcriber.leaves);
  } catch (const Unsupported& unsupported) {
    result = TranslatedProperty{};
    result.unsupportedReason = unsupported.reason;
  }
  return result;
}

}  // namespace

PropertyTranslator makeSvaPropertyTranslator() {
  return translate;
}

}  // namespace KEPLER_FORMAL::SEC
