// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "Btor2Writer.h"

using namespace KEPLER_FORMAL;

namespace {

class Btor2WriterTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }
};

// A deliberately small independent evaluator for the writer's combinational
// subset. It checks declaration order while interpreting the emitted operators.
bool evaluateBtor2(const std::string& text,
                   const std::unordered_map<std::string, bool>& inputs,
                   Btor2Writer::NodeId root) {
  std::unordered_map<Btor2Writer::NodeId, bool> values;
  std::istringstream lines(text);
  std::string line;
  Btor2Writer::NodeId previous = 0;
  while (std::getline(lines, line)) {
    if (line.empty() || line.front() == ';') {
      continue;
    }
    std::istringstream fields(line);
    Btor2Writer::NodeId id = 0;
    std::string op;
    fields >> id >> op;
    if (!fields || id <= previous) {
      throw std::runtime_error("Invalid BTOR2 declaration order");
    }
    previous = id;
    if (op == "sort") {
      std::string kind;
      unsigned width = 0;
      fields >> kind >> width;
      if (id != 1 || kind != "bitvec" || width != 1) {
        throw std::runtime_error("Unexpected BTOR2 sort");
      }
      continue;
    }
    Btor2Writer::NodeId sort = 0;
    fields >> sort;
    if (sort != 1) {
      throw std::runtime_error("Unexpected BTOR2 value sort");
    }
    if (op == "input") {
      std::string symbol;
      fields >> symbol;
      values.emplace(id, inputs.at(symbol));
    } else if (op == "const") {
      std::string bits;
      fields >> bits;
      if (bits != "0" && bits != "1") {
        throw std::runtime_error("Invalid BTOR2 constant");
      }
      values.emplace(id, bits == "1");
    } else {
      Btor2Writer::NodeId left = 0;
      Btor2Writer::NodeId right = 0;
      fields >> left;
      const bool a = values.at(left);
      if (op == "not") {
        values.emplace(id, !a);
        continue;
      }
      fields >> right;
      const bool b = values.at(right);
      if (op == "and") values.emplace(id, a && b);
      else if (op == "or") values.emplace(id, a || b);
      else if (op == "xor") values.emplace(id, a != b);
      else if (op == "eq") values.emplace(id, a == b);
      else throw std::runtime_error("Unexpected BTOR2 operator");
    }
  }
  return values.at(root);
}

}  // namespace

TEST_F(Btor2WriterTests, EmitsStateWiringAndPropertiesWithStandardSyntax) {
  std::ostringstream out;
  Btor2Writer writer(out);
  const auto input = writer.input("in");
  const auto zero = writer.constant(false);
  const auto state = writer.state("q");
  EXPECT_EQ(writer.constant(false), zero);
  writer.init(state, zero);
  writer.next(state, input);
  const auto mismatch = writer.logicalXor(state, input);
  writer.bad(mismatch, "mismatch");
  writer.constraint(input);
  writer.output(state, "out");

  EXPECT_EQ(out.str(),
      "1 sort bitvec 1\n"
      "2 input 1 in\n"
      "3 const 1 0\n"
      "4 state 1 q\n"
      "5 init 1 4 3\n"
      "6 next 1 4 2\n"
      "7 xor 1 4 2\n"
      "8 bad 7 mismatch\n"
      "9 constraint 2\n"
      "10 output 4 out\n");
}

TEST_F(Btor2WriterTests, PreservesBooleanSemanticsAndSharedCones) {
  std::ostringstream out;
  Btor2Writer writer(out);
  writer.bindVariable(2, writer.input("a"));
  writer.bindVariable(3, writer.input("b"));
  writer.bindVariable(4, writer.input("c"));
  BoolExpr* a = BoolExpr::Var(2);
  BoolExpr* b = BoolExpr::Var(3);
  BoolExpr* c = BoolExpr::Var(4);
  BoolExpr* shared = BoolExpr::And(a, b);
  BoolExpr* expr = BoolExpr::Or(
      BoolExpr::Xor(shared, c), BoolExpr::And(shared, BoolExpr::Not(c)));
  const auto root = writer.expression(expr);
  const std::string emitted = out.str();
  EXPECT_EQ(writer.expression(expr), root);
  writer.expression(shared);
  EXPECT_EQ(out.str(), emitted);
  // One sort, three inputs, and the five unique operation nodes.
  EXPECT_EQ(std::count(emitted.begin(), emitted.end(), '\n'), 9);
  for (size_t assignment = 0; assignment < 8; ++assignment) {
    const bool av = (assignment & 1) != 0;
    const bool bv = (assignment & 2) != 0;
    const bool cv = (assignment & 4) != 0;
    EXPECT_EQ(evaluateBtor2(emitted, {{"a", av}, {"b", bv}, {"c", cv}}, root),
              expr->evaluate({{2, av}, {3, bv}, {4, cv}}));
  }
}

TEST_F(Btor2WriterTests, WritesConstantsAndMonitorOperations) {
  std::ostringstream out;
  Btor2Writer writer(out);
  const auto a = writer.input("a");
  const auto b = writer.input("b");
  const auto falseNode = writer.expression(BoolExpr::createFalse());
  const auto trueNode = writer.expression(BoolExpr::createTrue());
  EXPECT_NE(falseNode, trueNode);
  EXPECT_EQ(falseNode, writer.constant(false));
  EXPECT_EQ(trueNode, writer.constant(true));
  const auto conjunction = writer.logicalAnd(a, writer.logicalNot(b));
  const auto eq = writer.equal(a, b);
  const auto root = writer.logicalOr(conjunction, eq);
  for (bool av : {false, true}) {
    for (bool bv : {false, true}) {
      const std::unordered_map<std::string, bool> inputs{{"a", av}, {"b", bv}};
      EXPECT_FALSE(evaluateBtor2(out.str(), inputs, falseNode));
      EXPECT_TRUE(evaluateBtor2(out.str(), inputs, trueNode));
      EXPECT_EQ(evaluateBtor2(out.str(), inputs, root), av || !bv);
    }
  }
}

TEST_F(Btor2WriterTests, ExportsDeepDagWithoutRecursion) {
  std::ostringstream out;
  Btor2Writer writer(out);
  writer.bindVariable(2, writer.input("a"));
  writer.bindVariable(3, writer.input("b"));
  BoolExpr* expr = BoolExpr::Var(2);
  constexpr size_t depth = 30000;
  for (size_t i = 0; i < depth; ++i) {
    expr = BoolExpr::Xor(expr, BoolExpr::Var(3));
  }
  const auto root = writer.expression(expr);
  const std::string emitted = out.str();
  EXPECT_EQ(static_cast<size_t>(std::count(emitted.begin(), emitted.end(), '\n')),
            depth + 3);
  EXPECT_TRUE(evaluateBtor2(emitted, {{"a", true}, {"b", true}}, root));
  EXPECT_FALSE(evaluateBtor2(emitted, {{"a", false}, {"b", true}}, root));
}

TEST_F(Btor2WriterTests, EncodesUnsafeSymbolsAndFlattensComments) {
  std::ostringstream out;
  Btor2Writer writer(out);
  writer.input("a b;\n\t%");
  writer.input("a%20b");
  writer.comment("first\n999 bad 1\r\tend");
  writer.state();
  EXPECT_EQ(out.str(),
      "1 sort bitvec 1\n"
      "2 input 1 a%20b%3B%0A%09%25\n"
      "3 input 1 a%2520b\n"
      "; first 999 bad 1  end\n"
      "4 state 1\n");
}

TEST_F(Btor2WriterTests, RejectsInvalidNodesAndVariableBindings) {
  std::ostringstream out;
  Btor2Writer writer(out);
  const auto input = writer.input("a");
  const auto zero = writer.constant(false);
  const auto state = writer.state("q");
  writer.bindVariable(2, input);
  EXPECT_NO_THROW(writer.bindVariable(2, input));
  EXPECT_THROW(writer.bindVariable(2, state), std::invalid_argument);
  EXPECT_THROW(writer.bindVariable(0, input), std::invalid_argument);
  EXPECT_THROW(writer.bindVariable(1, input), std::invalid_argument);
  EXPECT_THROW(writer.bindVariable(std::numeric_limits<size_t>::max(), input),
               std::invalid_argument);
  EXPECT_THROW(writer.bindVariable(3, 999), std::invalid_argument);
  EXPECT_THROW(writer.logicalNot(0), std::invalid_argument);
  EXPECT_THROW(writer.logicalAnd(input, 1), std::invalid_argument);
  EXPECT_THROW(writer.bad(999), std::invalid_argument);
  EXPECT_THROW(writer.init(input, state), std::invalid_argument);
  EXPECT_THROW(writer.next(state, 999), std::invalid_argument);
  const auto init = writer.init(state, zero);
  EXPECT_THROW(writer.output(init), std::invalid_argument);
  EXPECT_THROW(writer.init(state, input), std::invalid_argument);
  writer.next(state, input);
  EXPECT_THROW(writer.next(state, input), std::invalid_argument);
}

TEST_F(Btor2WriterTests, RejectsInitializationFromLaterNodesOrInputs) {
  std::ostringstream out;
  Btor2Writer writer(out);
  const auto input = writer.input("a");
  const auto inputCone = writer.logicalNot(input);
  const auto state = writer.state("q");
  const auto lateConstant = writer.constant(false);
  EXPECT_THROW(writer.init(state, input), std::invalid_argument);
  EXPECT_THROW(writer.init(state, inputCone), std::invalid_argument);
  EXPECT_THROW(writer.init(state, state), std::invalid_argument);
  EXPECT_THROW(writer.init(state, lateConstant), std::invalid_argument);
  const auto initialized = writer.state("r");
  EXPECT_NO_THROW(writer.init(initialized, state));
}

TEST_F(Btor2WriterTests, RejectsInvalidAndUnboundExpressions) {
  std::ostringstream out;
  Btor2Writer writer(out);
  BoolExpr invalid;
  EXPECT_THROW(writer.expression(nullptr), std::invalid_argument);
  EXPECT_THROW(writer.expression(&invalid), std::invalid_argument);
  EXPECT_THROW(writer.expression(BoolExpr::createInvalid()), std::invalid_argument);
  EXPECT_THROW(writer.expression(BoolExpr::Var(2)), std::invalid_argument);
  writer.bindVariable(2, writer.input("a"));
  EXPECT_NO_THROW(writer.expression(BoolExpr::Var(2)));
  BoolExpr* malformed = BoolExprCache::getExpression(
      {Op::AND, 0, BoolExpr::Var(2), nullptr});
  EXPECT_THROW(writer.expression(malformed), std::invalid_argument);
}

TEST_F(Btor2WriterTests, IgnoresStreamNumberFormattingAndReportsWriteErrors) {
  std::ostringstream out;
  out << std::hex << std::setw(100);
  Btor2Writer writer(out);
  for (unsigned i = 0; i < 10; ++i) {
    writer.input();
  }
  EXPECT_EQ(out.str().find("1 sort bitvec 1\n"), 0u);
  EXPECT_NE(out.str().find("11 input 1\n"), std::string::npos);
  out.setstate(std::ios::badbit);
  EXPECT_THROW(writer.input("fails"), std::runtime_error);
  EXPECT_THROW(writer.comment("fails"), std::runtime_error);

  std::ostringstream broken;
  broken.setstate(std::ios::badbit);
  EXPECT_THROW(Btor2Writer{broken}, std::runtime_error);
}
