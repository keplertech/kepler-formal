// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <functional>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "Tree2BoolExpr.h"
#include "phase/NormalizedTransitionSystem.h"

using namespace naja::NL;
using namespace KEPLER_FORMAL::SEC::PHASE;

namespace {

class NormalizedTransitionSystemTests : public ::testing::Test {
 protected:
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) {
      universe->destroy();
    }
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
};

SNLDesign* createDirectLatchTop(NLLibrary* library) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName("top"));
  auto* data =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("data"));
  auto* enable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("enable"));
  auto* output =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* latchModel = NLDB0::getDLatch();
  if (latchModel == nullptr) {
    throw std::runtime_error("DB0 latch primitive is unavailable");
  }
  auto* latch = SNLInstance::create(top, latchModel, NLName("latch0"));

  auto* dataNet = SNLScalarNet::create(top, NLName("data_net"));
  auto* enableNet = SNLScalarNet::create(top, NLName("enable_net"));
  auto* outputNet = SNLScalarNet::create(top, NLName("output_net"));
  data->setNet(dataNet);
  enable->setNet(enableNet);
  output->setNet(outputNet);
  latch->getInstTerm(NLDB0::getDLatchData())->setNet(dataNet);
  latch->getInstTerm(NLDB0::getDLatchEnable())->setNet(enableNet);
  latch->getInstTerm(NLDB0::getDLatchOutput())->setNet(outputNet);
  return top;
}

std::string formatDiagnostics(const NormalizedTransitionSystem& system) {
  std::ostringstream out;
  for (const auto& diagnostic : system.diagnostics()) {
    out << (diagnostic.fatal ? "fatal: " : "warning: ") << diagnostic.object
        << ": " << diagnostic.message << '\n';
  }
  return out.str();
}

size_t findVariable(const NormalizedTransitionSystem& system,
                    const std::string& displayName) {
  for (size_t index = 0; index < system.variables().size(); ++index) {
    if (system.variables()[index].displayName == displayName) {
      return index;
    }
  }
  throw std::runtime_error("Missing normalized variable `" + displayName + "`");
}

bool evaluateExpression(const ExpressionArena& arena, ExprID root,
                        const std::unordered_map<size_t, bool>& values) {
  if (!arena.isValid(root)) {
    throw std::runtime_error("Cannot evaluate an invalid phase expression");
  }
  std::vector<std::optional<bool>> evaluated(arena.nodes().size());
  std::function<bool(ExprID)> evaluate = [&](ExprID expression) {
    if (evaluated[expression].has_value()) {
      return *evaluated[expression];
    }
    const auto& node = arena.node(expression);
    bool result = false;
    switch (node.op) {
      case ExpressionOp::Constant:
        result = node.constant;
        break;
      case ExpressionOp::Variable: {
        const auto value = values.find(node.variable);
        if (value == values.end()) {
          throw std::runtime_error("Missing value for normalized variable " +
                                   std::to_string(node.variable));
        }
        result = value->second;
        break;
      }
      case ExpressionOp::Not:
        result = !evaluate(node.left);
        break;
      case ExpressionOp::And:
        result = evaluate(node.left) && evaluate(node.right);
        break;
      case ExpressionOp::Or:
        result = evaluate(node.left) || evaluate(node.right);
        break;
      case ExpressionOp::Xor:
        result = evaluate(node.left) != evaluate(node.right);
        break;
    }
    evaluated[expression] = result;
    return result;
  };
  return evaluate(root);
}

TEST_F(NormalizedTransitionSystemTests,
       CollectsDb0LatchWithOwnedHoldAndTransparencySemantics) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createDirectLatchTop(library);

  NormalizationOptions options;
  options.collectObservedOutputs = false;
  const auto system = collectNormalizedTransitionSystem(top, options);

  ASSERT_FALSE(system.hasFatalDiagnostics()) << formatDiagnostics(system);
  std::string validationReason;
  ASSERT_TRUE(system.validate(&validationReason)) << validationReason;
  ASSERT_EQ(system.states().size(), 1u);
  ASSERT_EQ(system.latches().size(), 1u);

  const auto& state = system.states().front();
  const auto& latch = system.latches().front();
  EXPECT_EQ(system.stateVariable(0).displayName, "latch0.Q[0]");
  EXPECT_EQ(state.initialValue, TernaryValue::Unknown);
  EXPECT_EQ(latch.stateIndex, 0u);
  EXPECT_EQ(latch.instancePath, "latch0");
  EXPECT_EQ(latch.modelName, "naja_dlatch");
  EXPECT_TRUE(system.expressions().isValid(latch.data));
  EXPECT_TRUE(system.expressions().isValid(latch.transparency));
  ASSERT_EQ(latch.outputs.size(), 1u);
  EXPECT_EQ(latch.outputs.front().displayName, "latch0.Q[0]");
  EXPECT_FALSE(latch.outputs.front().complemented);

  const size_t dataVariable = findVariable(system, "data[0]");
  const size_t enableVariable = findVariable(system, "enable[0]");
  const size_t stateVariable = state.variable;
  const auto evaluateNext = [&](bool data, bool enable, bool current) {
    return evaluateExpression(system.expressions(), state.nextState,
                              {{dataVariable, data},
                               {enableVariable, enable},
                               {stateVariable, current}});
  };

  EXPECT_FALSE(evaluateNext(true, false, false));
  EXPECT_TRUE(evaluateNext(false, false, true));
  EXPECT_FALSE(evaluateNext(false, true, true));
  EXPECT_TRUE(evaluateNext(true, true, false));
  EXPECT_FALSE(evaluateExpression(system.expressions(), latch.transparency,
                                  {{enableVariable, false}}));
  EXPECT_TRUE(evaluateExpression(system.expressions(), latch.transparency,
                                 {{enableVariable, true}}));
}

TEST_F(NormalizedTransitionSystemTests,
       AppliesExplicitResetStateOverrideWithoutInferringOneFromLatchModel) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createDirectLatchTop(library);

  NormalizationOptions options;
  options.collectObservedOutputs = false;
  options.initialValueByDisplayName.emplace("latch0.Q[0]", TernaryValue::One);
  const auto system = collectNormalizedTransitionSystem(top, options);

  ASSERT_FALSE(system.hasFatalDiagnostics()) << formatDiagnostics(system);
  ASSERT_EQ(system.states().size(), 1u);
  EXPECT_EQ(system.states().front().initialValue, TernaryValue::One);
}

}  // namespace
