// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "BoolExprCnfWriter.h"

using namespace KEPLER_FORMAL;

namespace {

class BoolExprCnfWriterTests : public ::testing::Test {
 protected:
  void TearDown() override {
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
};

BoolExpr* maybeNegated(BoolExpr* expr, bool positive) {
  return positive ? expr : BoolExpr::Not(expr);
}

BoolExpr* makeTwoLiteralTerm(BoolExpr* a,
                             bool aPositive,
                             BoolExpr* b,
                             bool bPositive) {
  return BoolExpr::And(
      maybeNegated(a, aPositive),
      maybeNegated(b, bPositive));
}

BoolExpr* makeThreeLiteralTerm(BoolExpr* a,
                               bool aPositive,
                               BoolExpr* b,
                               bool bPositive,
                               BoolExpr* c,
                               bool cPositive) {
  return BoolExpr::And(
      makeTwoLiteralTerm(a, aPositive, b, bPositive),
      maybeNegated(c, cPositive));
}

BoolExpr* makeOr3(BoolExpr* first, BoolExpr* second, BoolExpr* third) {
  return BoolExpr::Or(BoolExpr::Or(first, second), third);
}

BoolExpr* makeOr4(BoolExpr* first,
                  BoolExpr* second,
                  BoolExpr* third,
                  BoolExpr* fourth) {
  return BoolExpr::Or(makeOr3(first, second, third), fourth);
}

bool evaluateThreeInputExpr(BoolExpr* expr, bool a, bool b, bool c) {
  return expr->evaluate({{2, a}, {3, b}, {4, c}});
}

void expectEquivalentForThreeInputs(BoolExpr* actual, BoolExpr* expected) {
  for (size_t values = 0; values < 8; ++values) {
    const bool a = (values & 1u) != 0;
    const bool b = (values & 2u) != 0;
    const bool c = (values & 4u) != 0;
    EXPECT_EQ(
        evaluateThreeInputExpr(actual, a, b, c),
        evaluateThreeInputExpr(expected, a, b, c));
  }
}

std::string makeAllocationOrderSensitiveExprText(bool reverseOrder) {
  BoolExprCache::destroy();
  BoolExpr* lowVars = nullptr;
  BoolExpr* highVars = nullptr;
  if (reverseOrder) {
    highVars = BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(5));
    lowVars = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  } else {
    lowVars = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
    highVars = BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(5));
  }
  const std::string text = BoolExpr::Or(lowVars, highVars)->toString();
  BoolExprCache::destroy();
  return text;
}

}  // namespace

TEST_F(BoolExprCnfWriterTests, EncodeAndWriteAndExpression) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* expr = BoolExpr::And(a, b);

  CnfFormula cnf = encodeBoolExprToCnf(expr);

  EXPECT_EQ(cnf.numVars, 3);
  EXPECT_EQ(cnf.clauses.size(), 3u);
  EXPECT_TRUE(cnf.varNameToDimacs.count("x2") > 0);
  EXPECT_TRUE(cnf.varNameToDimacs.count("x3") > 0);
  EXPECT_GT(cnf.rootLit, 0);

  std::stringstream out;
  EXPECT_TRUE(writeDimacsCnf(cnf, out, true));

  const std::string text = out.str();
  EXPECT_NE(text.find("p cnf 3 4\n"), std::string::npos);
  EXPECT_NE(text.find("\n3 0\n"), std::string::npos);
}

TEST_F(BoolExprCnfWriterTests, WriteWithoutRootClause) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* expr = BoolExpr::And(a, b);

  CnfFormula cnf = encodeBoolExprToCnf(expr);

  std::stringstream out;
  EXPECT_TRUE(writeDimacsCnf(cnf, out, false));

  const std::string text = out.str();
  EXPECT_NE(text.find("p cnf 3 3\n"), std::string::npos);
  EXPECT_EQ(text.find("\n3 0\n"), std::string::npos);
}

TEST_F(BoolExprCnfWriterTests, EncodeConstantTrue) {
  BoolExpr* expr = BoolExpr::Var(1);
  CnfFormula cnf = encodeBoolExprToCnf(expr);

  EXPECT_EQ(cnf.numVars, 1);
  EXPECT_EQ(cnf.clauses.size(), 1u);
  EXPECT_EQ(cnf.rootLit, 1);
}

TEST_F(BoolExprCnfWriterTests, BoolExprCacheDestroyKeepsCacheReusable) {
  BoolExpr* beforeDestroy = BoolExpr::And(BoolExpr::Var(55), BoolExpr::Var(56));
  ASSERT_NE(beforeDestroy, nullptr);
  EXPECT_EQ(beforeDestroy->getOp(), Op::AND);

  BoolExprCache::destroy();

  // The production cache intentionally lives until process exit to avoid a
  // large shutdown-time destructor, but explicit cleanup must still support the
  // unit-test and embedded re-use paths.
  BoolExpr* afterDestroy = BoolExpr::And(BoolExpr::Var(55), BoolExpr::Var(56));
  ASSERT_NE(afterDestroy, nullptr);
  EXPECT_EQ(afterDestroy->getOp(), Op::AND);
  EXPECT_EQ(afterDestroy, BoolExpr::And(BoolExpr::Var(55), BoolExpr::Var(56)));
}

TEST_F(BoolExprCnfWriterTests,
       CommutativeExprOrderingIgnoresAllocationOrder) {
  EXPECT_EQ(
      makeAllocationOrderSensitiveExprText(/*reverseOrder=*/false),
      makeAllocationOrderSensitiveExprText(/*reverseOrder=*/true));
}

TEST_F(BoolExprCnfWriterTests, NormalizeArithmeticRewritesFullAdderSumDnf) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* c = BoolExpr::Var(4);
  BoolExpr* dnf = makeOr4(
      makeThreeLiteralTerm(a, true, b, false, c, false),
      makeThreeLiteralTerm(a, false, b, true, c, false),
      makeThreeLiteralTerm(a, false, b, false, c, true),
      makeThreeLiteralTerm(a, true, b, true, c, true));
  BoolExpr* expected = BoolExpr::Xor(BoolExpr::Xor(a, b), c);

  BoolExpr* normalized = BoolExpr::normalizeArithmetic(dnf);

  EXPECT_NE(normalized, dnf);
  expectEquivalentForThreeInputs(normalized, expected);
}

TEST_F(BoolExprCnfWriterTests, NormalizeArithmeticRewritesFullAdderCarryDnf) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* c = BoolExpr::Var(4);
  BoolExpr* dnf = makeOr4(
      makeThreeLiteralTerm(a, true, b, true, c, false),
      makeThreeLiteralTerm(a, true, b, false, c, true),
      makeThreeLiteralTerm(a, false, b, true, c, true),
      makeThreeLiteralTerm(a, true, b, true, c, true));
  BoolExpr* expected = BoolExpr::Or(
      BoolExpr::And(a, b),
      BoolExpr::And(c, BoolExpr::Xor(a, b)));

  BoolExpr* normalized = BoolExpr::normalizeArithmetic(dnf);

  EXPECT_NE(normalized, dnf);
  expectEquivalentForThreeInputs(normalized, expected);
}

TEST_F(BoolExprCnfWriterTests,
       NormalizeArithmeticRewritesMinimizedFullAdderCarry) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* c = BoolExpr::Var(4);
  BoolExpr* majority = makeOr3(
      makeTwoLiteralTerm(a, true, b, true),
      makeTwoLiteralTerm(a, true, c, true),
      makeTwoLiteralTerm(b, true, c, true));
  BoolExpr* expected = BoolExpr::Or(
      BoolExpr::And(a, b),
      BoolExpr::And(c, BoolExpr::Xor(a, b)));

  BoolExpr* normalized = BoolExpr::normalizeArithmetic(majority);

  EXPECT_NE(normalized, majority);
  expectEquivalentForThreeInputs(normalized, expected);
}

TEST_F(BoolExprCnfWriterTests,
       NormalizeArithmeticRejectsIncompleteSumLookalike) {
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* c = BoolExpr::Var(4);
  BoolExpr* lookalike = makeOr4(
      a,
      b,
      c,
      makeThreeLiteralTerm(a, true, b, true, c, true));

  BoolExpr* normalized = BoolExpr::normalizeArithmetic(lookalike);

  EXPECT_EQ(normalized, lookalike);
}

TEST_F(BoolExprCnfWriterTests, DumpToFileAndInvalidPath) {
  BoolExpr* expr = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));

  std::filesystem::path tmpDir = std::filesystem::temp_directory_path();
  std::filesystem::path filePath = tmpDir / "bool_expr_cnf_test.cnf";
  std::filesystem::path dirPath = tmpDir / "bool_expr_cnf_dir";

  if (std::filesystem::exists(filePath)) {
    std::filesystem::remove(filePath);
  }
  if (std::filesystem::exists(dirPath)) {
    std::filesystem::remove_all(dirPath);
  }
  std::filesystem::create_directory(dirPath);

  EXPECT_TRUE(dumpBoolExprToDimacs(expr, filePath.string()));
  EXPECT_TRUE(std::filesystem::exists(filePath));

  std::ifstream in(filePath);
  std::string header;
  std::getline(in, header);
  EXPECT_TRUE(header.rfind("p cnf", 0) == 0);

  CnfFormula cnf = encodeBoolExprToCnf(expr);
  EXPECT_FALSE(dumpDimacsCnf(cnf, dirPath.string()));

  std::filesystem::remove(filePath);
  std::filesystem::remove_all(dirPath);
}
