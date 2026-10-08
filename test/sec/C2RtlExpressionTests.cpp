// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "formal/C2RtlExpression.h"

using namespace KEPLER_FORMAL;
using namespace KEPLER_FORMAL::C2RTL;

namespace {

class C2RtlExpressionTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  void addSignal(const std::string& name, size_t width) {
    auto& bits = signals_[name];
    for (size_t bit = 0; bit < width; ++bit) {
      bits.push_back(BoolExpr::Var(nextId_++));
    }
  }

  void setValue(const std::string& name, uint64_t value) {
    size_t bit = 0;
    for (auto* node : signals_.at(name)) {
      values_[node->getId()] = bit < 64 && ((value >> bit) & 1);
      ++bit;
    }
  }

  BoolExpr* compile(const std::string& text) {
    return compileC2RtlExpression(text, [&](const std::string& name,
                                          std::optional<size_t> index) {
      const auto found = signals_.find(name);
      if (found == signals_.end()) {
        throw std::invalid_argument("unknown signal '" + name + "'");
      }
      if (!index) {
        return found->second;
      }
      if (*index >= found->second.size()) {
        throw std::invalid_argument("bit index outside signal '" + name + "'");
      }
      return ExpressionBits{found->second[*index]};
    });
  }

  std::unordered_map<std::string, ExpressionBits> signals_;
  std::unordered_map<size_t, bool> values_;
  size_t nextId_ = 2;
};

TEST_F(C2RtlExpressionTests, UnsignedComparisonsZeroExtendOperands) {
  addSignal("a", 3);
  addSignal("b", 5);
  auto* eq = compile("a == b");
  auto* ne = compile("a != b");
  auto* lt = compile("a < b");
  auto* le = compile("a <= b");
  auto* gt = compile("a > b");
  auto* ge = compile("a >= b");
  for (unsigned a = 0; a < 8; ++a) {
    for (unsigned b = 0; b < 32; ++b) {
      SCOPED_TRACE("a=" + std::to_string(a) + ", b=" + std::to_string(b));
      setValue("a", a);
      setValue("b", b);
      EXPECT_EQ(eq->evaluate(values_), a == b);
      EXPECT_EQ(ne->evaluate(values_), a != b);
      EXPECT_EQ(lt->evaluate(values_), a < b);
      EXPECT_EQ(le->evaluate(values_), a <= b);
      EXPECT_EQ(gt->evaluate(values_), a > b);
      EXPECT_EQ(ge->evaluate(values_), a >= b);
    }
  }
}

TEST_F(C2RtlExpressionTests, BooleanPrecedenceAndBusTruthiness) {
  addSignal("a", 2);
  addSignal("b", 2);
  auto* expr = compile("!a || b && a == 3");
  auto* grouped = compile("(!a || b) && a == 3");
  auto* comparison = compile("a == b < 2");
  for (unsigned a = 0; a < 4; ++a) {
    for (unsigned b = 0; b < 4; ++b) {
      setValue("a", a);
      setValue("b", b);
      EXPECT_EQ(expr->evaluate(values_), !a || (b && a == 3));
      EXPECT_EQ(grouped->evaluate(values_), (!a || b) && a == 3);
      EXPECT_EQ(comparison->evaluate(values_), a == (b < 2));
    }
  }
}

TEST_F(C2RtlExpressionTests, ConstraintsPreserveStrictNumericBoundaries) {
  addSignal("A_in", 32);
  auto* expr = compile("A_in > 0x10000000 && A_in < 0x7FFFF000");
  for (uint64_t value : {0ULL, 0x10000000ULL, 0x10000001ULL,
                         0x7FFFEFFFULL, 0x7FFFF000ULL, 0xFFFFFFFFULL}) {
    setValue("A_in", value);
    EXPECT_EQ(expr->evaluate(values_), value > 0x10000000 && value < 0x7FFFF000);
  }
}

TEST_F(C2RtlExpressionTests, NumericLiteralsAreNotTruncatedToSignalWidth) {
  addSignal("byte", 8);
  auto* lt = compile("byte < 256");
  auto* eq = compile("byte == 256");
  auto* gt = compile("0x100 > byte");
  for (unsigned value = 0; value < 256; ++value) {
    setValue("byte", value);
    EXPECT_TRUE(lt->evaluate(values_));
    EXPECT_FALSE(eq->evaluate(values_));
    EXPECT_TRUE(gt->evaluate(values_));
  }
}

TEST_F(C2RtlExpressionTests, ConstantsSupportMoreThan64BitsAndAllRadices) {
  EXPECT_TRUE(compile("1208925819614629174706176 == 0x100000000000000000000")
                  ->evaluate({}));
  EXPECT_TRUE(compile("0x100000000000000000000 > 18446744073709551615")
                  ->evaluate({}));
  EXPECT_TRUE(compile("0b101011 == 43 && 0X2B == 0B101011")->evaluate({}));
  EXPECT_TRUE(compile("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF == "
                      "340282366920938463463374607431768211455")
                  ->evaluate({}));
  EXPECT_TRUE(compile("!0 && !!0x100000000000000000000 && true && !false")
                  ->evaluate({}));
}

TEST_F(C2RtlExpressionTests, WideSignalsRetainBitsAbove64) {
  addSignal("wide", 81);
  auto* eq = compile("wide == 1208925819614629174706176");
  auto* lt = compile("wide < 0x100000000000000000000");
  setValue("wide", 0);
  EXPECT_FALSE(eq->evaluate(values_));
  EXPECT_TRUE(lt->evaluate(values_));
  values_[signals_.at("wide")[80]->getId()] = true;
  EXPECT_TRUE(eq->evaluate(values_));
  EXPECT_FALSE(lt->evaluate(values_));
  values_[signals_.at("wide")[0]->getId()] = true;
  EXPECT_FALSE(eq->evaluate(values_));
  EXPECT_FALSE(lt->evaluate(values_));
}

TEST_F(C2RtlExpressionTests, ContradictoryConstraintsRejectEveryInputValue) {
  addSignal("a", 4);
  auto* expr = compile("a > 4 && a < 3");
  for (unsigned value = 0; value < 16; ++value) {
    setValue("a", value);
    EXPECT_FALSE(expr->evaluate(values_));
  }
}

TEST_F(C2RtlExpressionTests, QualifiedSignalsAndBit40SelectAreResolved) {
  addSignal("model.mantissa", 48);
  addSignal("rtl.mantissa", 48);
  auto* bitEq = compile("model.mantissa[40] == rtl.mantissa[40]");
  auto* busEq = compile("model.mantissa == rtl.mantissa");
  setValue("model.mantissa", uint64_t{1} << 40);
  setValue("rtl.mantissa", (uint64_t{1} << 40) | 1);
  EXPECT_TRUE(bitEq->evaluate(values_));
  EXPECT_FALSE(busEq->evaluate(values_));
  setValue("rtl.mantissa", 1);
  EXPECT_FALSE(bitEq->evaluate(values_));
}

TEST_F(C2RtlExpressionTests, ConditionalNanInfEqualityIgnoresInactiveFields) {
  for (const auto& side : {"model.", "rtl."}) {
    addSignal(std::string(side) + "nan", 1);
    addSignal(std::string(side) + "inf", 1);
    addSignal(std::string(side) + "mantissa", 4);
    addSignal(std::string(side) + "exponent", 4);
  }
  auto* expr = compile(
      "model.nan == rtl.nan && "
      "(model.nan || model.inf == rtl.inf) && "
      "(model.nan || model.inf || model.mantissa == rtl.mantissa) && "
      "(model.nan || model.inf || !model.mantissa || "
      "model.exponent == rtl.exponent)");
  for (const auto& [name, bits] : signals_) {
    setValue(name, 0);
  }
  setValue("rtl.exponent", 15);
  EXPECT_TRUE(expr->evaluate(values_));  // Zero ignores exponent.
  setValue("model.mantissa", 2);
  setValue("rtl.mantissa", 2);
  EXPECT_FALSE(expr->evaluate(values_));
  setValue("model.inf", 1);
  setValue("rtl.inf", 1);
  setValue("rtl.mantissa", 3);
  EXPECT_TRUE(expr->evaluate(values_));  // Infinity ignores mantissa/exponent.
  setValue("rtl.inf", 0);
  EXPECT_FALSE(expr->evaluate(values_));
  setValue("model.nan", 1);
  setValue("rtl.nan", 1);
  EXPECT_TRUE(expr->evaluate(values_));  // NaN ignores the other outputs.
  setValue("rtl.nan", 0);
  EXPECT_FALSE(expr->evaluate(values_));
}

TEST_F(C2RtlExpressionTests, SyntaxValidationDoesNotRequireSignalResolution) {
  EXPECT_NO_THROW(validateC2RtlExpression(
      "unknown.input > 0x100 && (model.mantissa[40] == rtl.mantissa[40])"));
  EXPECT_NO_THROW(validateC2RtlExpression("true || unknown"));
  EXPECT_NO_THROW(validateC2RtlExpression("false && unknown"));
  EXPECT_NO_THROW(validateC2RtlExpression("0x" + std::string(1024, 'f')));
}

TEST_F(C2RtlExpressionTests, UnknownSignalsAreNotHiddenByConstantBranches) {
  EXPECT_THROW(compile("true || missing"), std::invalid_argument);
  EXPECT_THROW(compile("false && missing"), std::invalid_argument);
  try {
    compile("true || missing");
    FAIL() << "expected unknown signal error";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string(error.what()).find("offset 8"), std::string::npos);
    EXPECT_NE(std::string(error.what()).find("unknown signal 'missing'"),
              std::string::npos);
  }
}

TEST_F(C2RtlExpressionTests, RejectsMalformedExpressionsAndIndices) {
  for (const std::string text : {
           "", " ", "()", "a &&", "a ==", "a <", "a b", "a & b",
           "a | b", "a = b", "a + b", "a ^ b", "-1", "0x", "0b2",
           "0xg", "12abc", "a.", "a..b", "(a", "a)", "a[]", "a[-1]",
           "a[1.5]", "a[0x1]", "a[1", "a[1][2]", "true[0]",
           "a[9999999999999999999999999999999999999999]"}) {
    SCOPED_TRACE(text);
    EXPECT_THROW(validateC2RtlExpression(text), std::invalid_argument);
  }
  addSignal("a", 8);
  EXPECT_THROW(compile("a[8]"), std::invalid_argument);
  try {
    validateC2RtlExpression("a && @");
    FAIL() << "expected token error";
  } catch (const std::invalid_argument& error) {
    EXPECT_NE(std::string(error.what()).find("offset 5"), std::string::npos);
  }
}

TEST_F(C2RtlExpressionTests, RejectsInvalidResolvedSignalValues) {
  for (ExpressionBits bits : {ExpressionBits{}, ExpressionBits{nullptr},
                              ExpressionBits{BoolExpr::createInvalid()}}) {
    EXPECT_THROW(compileC2RtlExpression("a", [&](const auto&, auto) {
      return bits;
    }), std::invalid_argument);
  }
  EXPECT_THROW(compileC2RtlExpression("a[0]", [](const auto&, auto) {
    return ExpressionBits{BoolExpr::createFalse(), BoolExpr::createTrue()};
  }), std::invalid_argument);
}

TEST_F(C2RtlExpressionTests, RejectsExcessiveNestingLengthAndLiteralWidth) {
  EXPECT_THROW(validateC2RtlExpression(std::string(129, '!') + "true"),
               std::invalid_argument);
  EXPECT_THROW(validateC2RtlExpression(std::string(65537, ' ')),
               std::invalid_argument);
  EXPECT_THROW(validateC2RtlExpression("0x" + std::string(1025, 'f')),
               std::invalid_argument);
}

}  // namespace
