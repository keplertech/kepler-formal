// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

#include "BoolExpr.h"
#include "BoolExprCache.h"
#include "export/SecBtor2Exporter.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {
namespace {

class SecBtor2ExportValidationTests : public ::testing::Test {
 protected:
  void TearDown() override { BoolExprCache::destroy(); }

  static KInductionProblem preparedProblem() {
    KInductionProblem problem;
    problem.property = BoolExpr::createTrue();
    problem.bad = BoolExpr::createFalse();
    problem.state0Symbols = {2};
    problem.transitions0 = {{2, BoolExpr::Var(2)}};
    problem.inputSymbols = {3};
    problem.allSymbols = {2, 3};
    problem.totalStateCount = 1;
    return problem;
  }

  static std::string dump(const KInductionProblem& problem) {
    std::ostringstream output;
    exportSecBtor2(problem, output);
    return output.str();
  }

  static void expectError(const KInductionProblem& problem,
                          const std::string& message) {
    try {
      dump(problem);
      FAIL() << "Expected export to reject the malformed SEC problem";
    } catch (const std::runtime_error& error) {
      EXPECT_EQ(error.what(), message);
    }
  }
};

TEST_F(SecBtor2ExportValidationTests, RejectsMissingPreparedPropertyBeforeWriting) {
  for (bool missingBad : {false, true}) {
    auto problem = preparedProblem();
    (missingBad ? problem.bad : problem.property) = nullptr;
    std::ostringstream output;
    EXPECT_THROW(exportSecBtor2(problem, output), std::runtime_error);
    EXPECT_TRUE(output.str().empty());
  }
}

TEST_F(SecBtor2ExportValidationTests, ReportsCoveredAndSkippedOutputMetadata) {
  auto problem = preparedProblem();
  problem.observedOutputNames = {"data", "ready"};
  std::ostringstream output;
  exportSecBtor2(problem, output, {4, {"opaque_output", "undriven_output"}});
  const auto text = output.str();
  EXPECT_NE(text.find("covered_outputs=2 total_outputs=4"), std::string::npos);
  EXPECT_NE(text.find("covered_output=data"), std::string::npos);
  EXPECT_NE(text.find("covered_output=ready"), std::string::npos);
  EXPECT_NE(text.find("skipped_output=opaque_output"), std::string::npos);
  EXPECT_NE(text.find("skipped_output=undriven_output"), std::string::npos);
  EXPECT_NE(dump(problem).find("covered_outputs=2 total_outputs=2"),
            std::string::npos);
}

TEST_F(SecBtor2ExportValidationTests, RejectsReservedAndStateAliasingInputs) {
  for (size_t symbol : {0, 1, 2}) {
    auto problem = preparedProblem();
    problem.inputSymbols = {symbol};
    expectError(problem, "Invalid SEC input symbol " + std::to_string(symbol));
  }
}

TEST_F(SecBtor2ExportValidationTests, RejectsReservedStateSymbols) {
  for (size_t symbol : {0, 1}) {
    auto problem = preparedProblem();
    problem.state0Symbols = {symbol};
    expectError(problem,
                "Invalid or duplicate SEC state symbol " + std::to_string(symbol));
  }
}

TEST_F(SecBtor2ExportValidationTests, RejectsRelationsToInputsOrUndeclaredStates) {
  for (size_t symbol : {3, 99}) {
    auto problem = preparedProblem();
    problem.sameFrameStateEqualityPairs0 = {{2, symbol}};
    expectError(problem,
                "BTOR2 export references undeclared state " + std::to_string(symbol));
  }
}

TEST_F(SecBtor2ExportValidationTests, RejectsInitializationOfAnUndeclaredState) {
  auto problem = preparedProblem();
  problem.initializedStateCount = 1;
  problem.initialCondition = BoolExpr::createTrue();
  problem.initialStateAssignments = {{99, false}};
  expectError(problem, "BTOR2 export references undeclared state 99");
}

TEST_F(SecBtor2ExportValidationTests, ChecksResetCounterOverflowBoundary) {
  auto problem = preparedProblem();
  problem.resetBootstrapInputs = {{3, true}};
  problem.resetBootstrapCycles = std::numeric_limits<size_t>::max();
  expectError(problem, "BTOR2 export reset cycle count is too large");
  --problem.resetBootstrapCycles;
  EXPECT_NO_THROW(dump(problem));
}

TEST_F(SecBtor2ExportValidationTests, ChecksConsistencyOfRepeatedInitialFacts) {
  auto problem = preparedProblem();
  problem.initializedStateCount = 1;
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initialStateAssignments = {{2, false}};
  const auto singleAssignment = dump(problem);
  problem.initialStateAssignments.push_back({2, false});
  EXPECT_EQ(dump(problem), singleAssignment);
  problem.initialStateAssignments.back().second = true;
  expectError(problem, "Conflicting SEC initial assignments");
}

TEST_F(SecBtor2ExportValidationTests, RejectsResetPortsThatAreNotInputs) {
  for (size_t symbol : {2, 99}) {
    auto problem = preparedProblem();
    problem.resetBootstrapInputs = {{symbol, true}};
    problem.resetBootstrapCycles = 1;
    expectError(problem, "BTOR2 export reset port is not an input");
  }
}

TEST_F(SecBtor2ExportValidationTests, ReportsStreamFailureDuringExport) {
  class ShortWriteBuffer : public std::streambuf {
   public:
    explicit ShortWriteBuffer(std::streamsize budget) : budget_(budget) {}

   protected:
    std::streamsize xsputn(const char*, std::streamsize count) override {
      const auto written = std::min(count, budget_);
      budget_ -= written;
      return written;
    }

   private:
    std::streamsize budget_;
  } buffer(static_cast<std::streamsize>(dump(preparedProblem()).size() / 2));
  std::ostream output(&buffer);
  EXPECT_THROW(exportSecBtor2(preparedProblem(), output), std::runtime_error);
  EXPECT_TRUE(output.bad());
}

TEST_F(SecBtor2ExportValidationTests, RejectsEmptyFilePath) {
  try {
    exportSecBtor2File(preparedProblem(), "");
    FAIL() << "Expected export to reject the empty destination";
  } catch (const std::runtime_error& error) {
    EXPECT_STREQ(error.what(), "BTOR2 export path must not be empty");
  }
}

}  // namespace
}  // namespace KEPLER_FORMAL::SEC
