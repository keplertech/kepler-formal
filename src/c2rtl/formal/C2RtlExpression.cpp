// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "C2RtlExpression.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>

#include "BoolExpr.h"

namespace KEPLER_FORMAL::C2RTL {
namespace {

class Parser {
 public:
  Parser(std::string_view text, const ExpressionSignalResolver* resolve)
      : text_(text), resolve_(resolve) {
    if (text.size() > 65536) {
      fail(0, "expression exceeds 64 KiB");
    }
  }

  BoolExpr* parse() {
    auto value = logicalOr();
    whitespace();
    if (position_ != text_.size()) {
      fail(position_, "unexpected token");
    }
    return resolve_ ? truth(value) : nullptr;
  }

 private:
  [[noreturn]] void fail(size_t offset, const std::string& message) const {
    throw std::invalid_argument("expression at offset " +
                                std::to_string(offset) + ": " + message);
  }

  void whitespace() {
    while (position_ < text_.size() &&
           std::isspace(static_cast<unsigned char>(text_[position_]))) {
      ++position_;
    }
  }

  bool take(std::string_view token) {
    whitespace();
    if (text_.substr(position_, token.size()) != token) {
      return false;
    }
    position_ += token.size();
    return true;
  }

  static bool nameStart(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '$';
  }

  static bool namePart(char c) {
    return nameStart(c) || std::isdigit(static_cast<unsigned char>(c));
  }

  static BoolExpr* truth(const ExpressionBits& value) {
    auto* result = BoolExpr::createFalse();
    for (auto* bit : value) {
      result = BoolExpr::Or(result, bit);
    }
    return result;
  }

  static BoolExpr* bit(const ExpressionBits& value, size_t index) {
    return index < value.size() ? value[index] : BoolExpr::createFalse();
  }

  static BoolExpr* equal(const ExpressionBits& lhs, const ExpressionBits& rhs) {
    auto* result = BoolExpr::createTrue();
    for (size_t i = 0; i < std::max(lhs.size(), rhs.size()); ++i) {
      result = BoolExpr::And(
          result, BoolExpr::Not(BoolExpr::Xor(bit(lhs, i), bit(rhs, i))));
    }
    return result;
  }

  static BoolExpr* less(const ExpressionBits& lhs, const ExpressionBits& rhs) {
    // A more significant differing bit overrides the accumulated lower bits.
    auto* result = BoolExpr::createFalse();
    for (size_t i = 0; i < std::max(lhs.size(), rhs.size()); ++i) {
      auto* a = bit(lhs, i);
      auto* b = bit(rhs, i);
      result = BoolExpr::Or(
          BoolExpr::And(BoolExpr::Not(a), b),
          BoolExpr::And(BoolExpr::Not(BoolExpr::Xor(a, b)), result));
    }
    return result;
  }

  ExpressionBits logicalOr() {
    auto lhs = logicalAnd();
    while (take("||")) {
      auto rhs = logicalAnd();
      if (resolve_) {
        lhs = {BoolExpr::Or(truth(lhs), truth(rhs))};
      }
    }
    return lhs;
  }

  ExpressionBits logicalAnd() {
    auto lhs = equality();
    while (take("&&")) {
      auto rhs = equality();
      if (resolve_) {
        lhs = {BoolExpr::And(truth(lhs), truth(rhs))};
      }
    }
    return lhs;
  }

  ExpressionBits equality() {
    auto lhs = comparison();
    while (true) {
      const bool isEqual = take("==");
      if (!isEqual && !take("!=")) {
        return lhs;
      }
      auto rhs = comparison();
      if (resolve_) {
        auto* result = equal(lhs, rhs);
        lhs = {isEqual ? result : BoolExpr::Not(result)};
      }
    }
  }

  ExpressionBits comparison() {
    auto lhs = unary();
    while (true) {
      std::string_view operation;
      for (auto candidate : {"<=", ">=", "<", ">"}) {
        if (take(candidate)) {
          operation = candidate;
          break;
        }
      }
      if (operation.empty()) {
        return lhs;
      }
      auto rhs = unary();
      if (resolve_) {
        auto* result = (operation == "<" || operation == ">=")
                           ? less(lhs, rhs)
                           : less(rhs, lhs);
        lhs = {(operation == "<=" || operation == ">=")
                   ? BoolExpr::Not(result)
                   : result};
      }
    }
  }

  ExpressionBits unary() {
    whitespace();
    if (++depth_ > 128) {
      fail(position_, "nesting exceeds 128 operands");
    }
    ExpressionBits result;
    if (take("!")) {
      auto operand = unary();
      if (resolve_) {
        result = {BoolExpr::Not(truth(operand))};
      }
    } else if (take("(")) {
      result = logicalOr();
      if (!take(")")) {
        fail(position_, "expected ')'");
      }
    } else if (position_ < text_.size() && nameStart(text_[position_])) {
      result = signal();
    } else if (position_ < text_.size() &&
               std::isdigit(static_cast<unsigned char>(text_[position_]))) {
      result = number();
    } else {
      fail(position_, "expected a signal, numeric constant or '('");
    }
    --depth_;
    return result;
  }

  ExpressionBits signal() {
    const size_t start = position_;
    do {
      ++position_;
      while (position_ < text_.size() && namePart(text_[position_])) {
        ++position_;
      }
      if (position_ == text_.size() || text_[position_] != '.') {
        break;
      }
      ++position_;
      if (position_ == text_.size() || !nameStart(text_[position_])) {
        fail(position_, "expected signal name after '.'");
      }
    } while (true);
    const std::string name(text_.substr(start, position_ - start));
    if (name == "true" || name == "false") {
      return resolve_ ? ExpressionBits{BoolExpr::Var(name == "true")}
                      : ExpressionBits{};
    }

    std::optional<size_t> index;
    if (take("[")) {
      whitespace();
      const size_t indexStart = position_;
      size_t value = 0;
      while (position_ < text_.size() &&
             std::isdigit(static_cast<unsigned char>(text_[position_]))) {
        const auto digit = static_cast<size_t>(text_[position_] - '0');
        if (value > (std::numeric_limits<size_t>::max() - digit) / 10) {
          fail(indexStart, "bit index is out of range");
        }
        value = value * 10 + digit;
        ++position_;
      }
      if (position_ == indexStart) {
        fail(position_, "expected a nonnegative decimal bit index");
      }
      if (!take("]")) {
        fail(position_, "expected ']' after bit index");
      }
      index = value;
    }
    if (!resolve_) {
      return {};
    }
    ExpressionBits value;
    try {
      value = (*resolve_)(name, index);
    } catch (const std::exception& error) {
      fail(start, error.what());
    }
    if (value.empty() || std::any_of(value.begin(), value.end(), [](auto* node) {
          return node == nullptr || !node->isValid();
        })) {
      fail(start, "signal '" + name + "' has no valid value");
    }
    if (index && value.size() != 1) {
      fail(start, "bit selection must resolve to exactly one bit");
    }
    return value;
  }

  ExpressionBits number() {
    const size_t start = position_;
    unsigned base = 10;
    if (text_.substr(position_, 2) == "0x" ||
        text_.substr(position_, 2) == "0X") {
      base = 16;
      position_ += 2;
    } else if (text_.substr(position_, 2) == "0b" ||
               text_.substr(position_, 2) == "0B") {
      base = 2;
      position_ += 2;
    }
    const size_t digitsStart = position_;
    std::vector<unsigned char> bits(1, 0);
    while (position_ < text_.size() && namePart(text_[position_])) {
      const char c = text_[position_];
      unsigned digit = c >= '0' && c <= '9' ? c - '0'
                       : c >= 'a' && c <= 'f' ? c - 'a' + 10
                       : c >= 'A' && c <= 'F' ? c - 'A' + 10
                                            : base;
      if (digit >= base) {
        fail(position_, "invalid digit in numeric constant");
      }
      unsigned carry = digit;
      for (auto& bitValue : bits) {
        carry += bitValue * base;
        bitValue = carry & 1;
        carry >>= 1;
      }
      while (carry) {
        bits.push_back(carry & 1);
        carry >>= 1;
      }
      if (bits.size() > 4096) {
        fail(start, "numeric constant exceeds 4096 bits");
      }
      ++position_;
    }
    if (position_ == digitsStart) {
      fail(position_, "expected digits in numeric constant");
    }
    ExpressionBits value;
    if (resolve_) {
      value.reserve(bits.size());
      for (auto bitValue : bits) {
        value.push_back(BoolExpr::Var(bitValue));
      }
    }
    return value;
  }

  std::string_view text_;
  const ExpressionSignalResolver* resolve_;
  size_t position_ = 0;
  size_t depth_ = 0;
};

}  // namespace

BoolExpr* compileC2RtlExpression(
    std::string_view text, const ExpressionSignalResolver& resolve) {
  return Parser(text, &resolve).parse();
}

void validateC2RtlExpression(std::string_view text) {
  Parser(text, nullptr).parse();
}

}  // namespace KEPLER_FORMAL::C2RTL
