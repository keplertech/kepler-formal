// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

#include "LibertyLatchModels.h"
#include "NLUniverse.h"
#include "NLLibrary.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLLibertyConstructor.h"
#include "SNLScalarTerm.h"

namespace {

using namespace naja::NL;
using Modeling = SNLDesignModeling;
using Expression = Modeling::BooleanExpression;
using Conflict = Modeling::SequentialState::ClearPresetValue;

bool evaluate(const Expression& expression,
              const std::map<std::string, bool>& inputs,
              const std::vector<bool>& states = {}) {
  const auto visit = [&](const auto& self, size_t id) -> bool {
    const auto& node = expression.nodes.at(id);
    switch (node.operation) {
      case Expression::Operator::Constant: return node.constant;
      case Expression::Operator::Term: return inputs.at(node.term->getName().getString());
      case Expression::Operator::State: return states.at(node.state);
      case Expression::Operator::Not: return !self(self, node.operands.at(0));
      case Expression::Operator::And: {
        bool value = true;
        for (const auto operand : node.operands) value &= self(self, operand);
        return value;
      }
      case Expression::Operator::Or: {
        bool value = false;
        for (const auto operand : node.operands) value |= self(self, operand);
        return value;
      }
      case Expression::Operator::Xor: {
        bool value = false;
        for (const auto operand : node.operands) value ^= self(self, operand);
        return value;
      }
    }
    throw std::runtime_error("unknown expression operator");
  };
  return visit(visit, expression.root);
}

class LibertyLatchModelsTests : public ::testing::Test {
 protected:
  void SetUp() override {
    auto* db = NLDB::create(NLUniverse::create());
    library_ = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("cells"));
    directory_ = std::filesystem::temp_directory_path() /
        ("kepler_latch_liberty_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory_);
  }
  void TearDown() override {
    if (auto* universe = NLUniverse::get()) universe->destroy();
    std::filesystem::remove_all(directory_);
  }
  std::string cell(const std::string& behavior,
                   const std::string& function = "IQ",
                   const std::string& extra = "",
                   const std::string& name = "LATCH") {
    return "cell(" + name + ") { " + behavior +
        " pin(D) { direction : input; } pin(E) { direction : input; }"
        " pin(CLK) { direction : input; } pin(R) { direction : input; }"
        " pin(S) { direction : input; } pin(Q) { direction : output; function : \"" +
        function + "\"; } " + extra + " } ";
  }
  SNLDesign* load(const std::string& cells, bool gzip = false) {
    const auto path = directory_ / ("cells" + std::to_string(nextFile_++) + ".lib");
    const auto contents = "library(test) { " + cells + " }";
    if (gzip) {
      auto* file = gzopen(path.string().c_str(), "wb");
      if (!file) throw std::runtime_error("cannot create gzip fixture");
      const auto count = gzwrite(file, contents.data(), contents.size());
      const auto status = gzclose(file);
      if (count != contents.size() || status != Z_OK) throw std::runtime_error("gzip write failed");
    } else {
      std::ofstream(path) << contents;
    }
    KEPLER_FORMAL::constructLibertyWithLatchModels(library_, path);
    return library_->getSNLDesign(NLName("LATCH"));
  }
  const Modeling::SequentialModel& model(SNLDesign* design) {
    return Modeling::getSequentialModel(design);
  }
  NLLibrary* library_ = nullptr;
  std::filesystem::path directory_;
  unsigned nextFile_ = 0;
};

TEST_F(LibertyLatchModelsTests, ParsesActiveLowLatchDataAndComplementedOutput) {
  auto* design = load(cell("latch(IQ, IQN) { enable : \"!E\"; data_in : \"D\"; }",
                          "IQN"));
  ASSERT_TRUE(Modeling::hasSequentialModel(design));
  const auto& latch = model(design);
  EXPECT_EQ(latch.kind, Modeling::SequentialModel::Kind::Latch);
  ASSERT_EQ(latch.states.size(), 1u);
  EXPECT_TRUE(evaluate(latch.clockedOn, {{"E", false}}));
  EXPECT_FALSE(evaluate(latch.clockedOn, {{"E", true}}));
  EXPECT_TRUE(evaluate(latch.states[0].nextState, {{"D", true}}));
  EXPECT_FALSE(evaluate(latch.outputs[0].function, {}, {true}));
  EXPECT_FALSE(Modeling::getOutputRelatedClocks(design->getScalarTerm(NLName("Q"))).empty());
}

TEST_F(LibertyLatchModelsTests, ClockGateUsesActualLatchAndOutputExpressions) {
  auto* design = load(cell("latch(IQ, IQN) { enable : \"!CLK\"; data_in : \"D | E\"; }",
                          "IQ & CLK", "clock_gating_integrated_cell : latch_posedge;"));
  ASSERT_TRUE(Modeling::hasSequentialModel(design));
  const auto& latch = model(design);
  EXPECT_TRUE(evaluate(latch.clockedOn, {{"CLK", false}}));
  EXPECT_TRUE(evaluate(latch.states[0].nextState, {{"D", false}, {"E", true}}));
  EXPECT_FALSE(evaluate(latch.outputs[0].function, {{"CLK", false}}, {true}));
  EXPECT_TRUE(evaluate(latch.outputs[0].function, {{"CLK", true}}, {true}));
}

TEST_F(LibertyLatchModelsTests, ResetPresetAndConflictModesAreRetained) {
  const std::vector<std::pair<std::string, Conflict>> cases = {
      {"L", Conflict::Zero}, {"H", Conflict::One}, {"N", Conflict::Hold},
      {"T", Conflict::Toggle}, {"X", Conflict::Unknown}};
  for (size_t i = 0; i < cases.size(); ++i) {
    const auto& [value, expected] = cases[i];
    const auto second = value == "L" ? "H" : value == "H" ? "L" : value;
    const auto name = "LATCH" + std::to_string(i);
    load(cell("latch(IQ, IQN) { enable : E; data_in : D; clear : R; preset : S; "
              "clear_preset_var1 : " + value + "; clear_preset_var2 : " + second + "; }",
              "IQ", "", name));
    auto* design = library_->getSNLDesign(NLName(name));
    ASSERT_TRUE(Modeling::hasSequentialModel(design));
    const auto& state = model(design).states[0];
    EXPECT_EQ(state.clearPresetValue, expected);
    ASSERT_TRUE(state.clear.has_value());
    ASSERT_TRUE(state.preset.has_value());
    EXPECT_TRUE(evaluate(*state.clear, {{"R", true}}));
    EXPECT_FALSE(evaluate(*state.preset, {{"S", false}}));
  }
}

TEST_F(LibertyLatchModelsTests, MultipleGroupsRequireSharedEnable) {
  auto* design = load(cell("latch(IQ, IQN) { enable : E; data_in : D; }"
                          "latch(JQ, JQN) { enable : E; data_in : !D; }", "IQ ^ JQ"));
  ASSERT_TRUE(Modeling::hasSequentialModel(design));
  ASSERT_EQ(model(design).states.size(), 2u);
  EXPECT_TRUE(evaluate(model(design).outputs[0].function, {}, {true, false}));
  EXPECT_FALSE(evaluate(model(design).states[1].nextState, {{"D", true}}));
}

TEST_F(LibertyLatchModelsTests, UnsupportedLatchDescriptionsStayUnmodeled) {
  const std::vector<std::string> unsupported = {
      "latch(IQ, IQN) { data_in : D; }",
      "latch(IQ, IQN) { enable : E; }",
      "latch(IQ, IQN) { enable : E; data_in : MISSING; }",
      "latch(IQ, IQN) { enable : E; data_in : Q; }",
      "latch(IQ, IQN) { enable : E; data_in : D; enable_also : CLK; }",
      "latch(IQ) { enable : E; data_in : D; clear_preset_var2 : H; }",
      "latch(IQ, IQN) { enable : E; data_in : D; clear : R; preset : S; "
          "clear_preset_var1 : L; clear_preset_var2 : L; }",
      "latch(IQ, IQN) { enable : E; data_in : D; clear : R; preset : S; "
          "clear_preset_var1 : H; }",
      "latch(IQ, IQN) { enable : E; data_in : D; }"
          "latch(JQ, JQN) { enable : CLK; data_in : D; }",
      "latch(IQ, IQN) { enable : E; data_in : D; }"
          "latch(IQ, JQN) { enable : E; data_in : D; }",
      "latch(D, IQN) { enable : E; data_in : D; }",
      "latch(IQ, IQN) { enable : E; data_in : D; }"
          "ff(FQ, FQN) { clocked_on : CLK; next_state : D; }",
      "latch(IQ, IQN) { enable : E; data_in : D; }"
          "statetable(\"D E\", \"IQ\") { table : \"- - : - : -\"; }",
      "latch(IQ, IQN) { enable : E; data_in : D; } power_down_function : R;",
  };
  for (size_t i = 0; i < unsupported.size(); ++i) {
    const auto name = "UNSUPPORTED" + std::to_string(i);
    SCOPED_TRACE(unsupported[i]);
    load(cell(unsupported[i], "IQ", "", name));
    auto* design = library_->getSNLDesign(NLName(name));
    ASSERT_NE(design, nullptr);
    EXPECT_FALSE(Modeling::hasSequentialModel(design));
  }
}

TEST_F(LibertyLatchModelsTests, DoesNotInferLatchFromClockGateMetadataOrNames) {
  auto* design = load(cell("", "D & CLK",
                          "clock_gating_integrated_cell : latch_posedge;"));
  EXPECT_FALSE(Modeling::hasSequentialModel(design));
}

TEST_F(LibertyLatchModelsTests, UnsupportedOutputPinsAndBusRemainUnmodeled) {
  const std::vector<std::string> extraPins = {
      "pin(Z) { direction : output; }",
      "pin(Z) { direction : output; function : IQ; three_state : S; }",
      "pin(Z) { direction : inout; }",
      "bus(B) { direction : input; bus_type : two_bits; }",
  };
  for (size_t i = 0; i < extraPins.size(); ++i) {
    const auto name = "PINS" + std::to_string(i);
    load("type(two_bits) { base_type : array; data_type : bit; bit_width : 2; "
         "bit_from : 1; bit_to : 0; downto : true; } " +
         cell("latch(IQ, IQN) { enable : E; data_in : D; }", "IQ", extraPins[i], name));
    auto* design = library_->getSNLDesign(NLName(name));
    ASSERT_NE(design, nullptr);
    EXPECT_FALSE(Modeling::hasSequentialModel(design));
  }
}

TEST_F(LibertyLatchModelsTests, ExistingDefinitionWinsAcrossFiles) {
  auto* first = load(cell("latch(IQ, IQN) { enable : E; data_in : D; }"));
  const auto before = model(first).clockedOn;
  auto* second = load(cell("latch(IQ, IQN) { enable : !E; data_in : !D; }"));
  EXPECT_EQ(first, second);
  EXPECT_TRUE(evaluate(model(second).clockedOn, {{"E", true}}));
  EXPECT_TRUE(evaluate(before, {{"E", true}}));
}

TEST_F(LibertyLatchModelsTests, EarlierUnsupportedDuplicateIsNotUpgraded) {
  auto* design = load(cell("latch(IQ, IQN) { enable : E; }") +
                      cell("latch(IQ, IQN) { enable : E; data_in : D; }"));
  EXPECT_FALSE(Modeling::hasSequentialModel(design));
  load(cell("latch(IQ, IQN) { enable : E; data_in : D; }"));
  EXPECT_FALSE(Modeling::hasSequentialModel(design));
}

TEST_F(LibertyLatchModelsTests, GzipSignatureIsSupportedWithoutGzipExtension) {
  auto* design = load(cell("latch(IQ, IQN) { enable : E; data_in : D; }"), true);
  EXPECT_TRUE(Modeling::hasSequentialModel(design));
}

TEST_F(LibertyLatchModelsTests, OrdinaryFrontendBehaviorRemainsUnchanged) {
  const auto path = directory_ / "ordinary.lib";
  std::ofstream(path) << "library(test) { "
      << cell("latch(IQ, IQN) { enable : E; data_in : D; }") << " }";
  SNLLibertyConstructor(library_).construct(path);
  auto* design = library_->getSNLDesign(NLName("LATCH"));
  EXPECT_FALSE(Modeling::hasSequentialModel(design));
  KEPLER_FORMAL::constructLibertyWithLatchModels(library_, path);
  EXPECT_FALSE(Modeling::hasSequentialModel(design));
}

TEST_F(LibertyLatchModelsTests, RejectsMissingInputAndNullLibrary) {
  EXPECT_THROW(KEPLER_FORMAL::constructLibertyWithLatchModels(
      library_, directory_ / "missing.lib"), std::exception);
  EXPECT_THROW(KEPLER_FORMAL::constructLibertyWithLatchModels(
      nullptr, directory_ / "missing.lib"), std::invalid_argument);
}

}  // namespace
