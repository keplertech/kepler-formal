// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <string>
#include <utility>

#include "BoolExprCache.h"
#include "formal/C2RtlEquivalenceStrategy.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "Tree2BoolExpr.h"

namespace KEPLER_FORMAL::C2RTL {
namespace {

using namespace naja::NL;

class C2RtlEquivalenceTests : public ::testing::Test {
 protected:
  enum class Source { A, B, Zero, InvertedA };

  void SetUp() override {
    auto* universe = NLUniverse::create();
    auto* db = NLDB::create(universe);
    designs_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    auto* primitives =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
    inverter_ = SNLDesign::create(primitives, SNLDesign::Type::Primitive,
                                  NLName("INV"));
    auto* a = SNLScalarTerm::create(inverter_, SNLTerm::Direction::Input,
                                   NLName("A"));
    auto* y = SNLScalarTerm::create(inverter_, SNLTerm::Direction::Output,
                                   NLName("Y"));
    SNLDesignModeling::addCombinatorialArcs({a}, {y});
    SNLDesignModeling::setTruthTable(inverter_, SNLTruthTable::Inv());
  }

  void TearDown() override {
    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    Tree2BoolExpr::iso2boolExpr_.clear();
    BoolExprCache::destroy();
  }

  SNLScalarNet* input(SNLDesign* top, const char* name) {
    auto* net = SNLScalarNet::create(top, NLName(name));
    SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName(name))
        ->setNet(net);
    return net;
  }

  SNLDesign* design(const std::string& name, Source source = Source::A,
                    size_t stages = 0, bool resettable = false) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* a = input(top, "a");
    auto* b = input(top, "b");
    auto* value = source == Source::B ? b : a;
    if (source == Source::Zero) {
      value = SNLScalarNet::create(top, NLName("zero"));
      value->setType(SNLNet::Type::Assign0);
    } else if (source == Source::InvertedA) {
      value = SNLScalarNet::create(top, NLName("inverted_a"));
      auto* inv = SNLInstance::create(top, inverter_, NLName("inv"));
      inv->getInstTerm(inverter_->getScalarTerm(NLName("A")))->setNet(a);
      inv->getInstTerm(inverter_->getScalarTerm(NLName("Y")))->setNet(value);
    }
    if (stages) {
      auto* clk = input(top, "clk");
      auto* rst = resettable ? input(top, "rst") : nullptr;
      for (size_t stage = 0; stage < stages; ++stage) {
        const auto suffix = std::to_string(stage);
        auto* ff = SNLInstance::create(
            top, resettable ? NLDB0::getDFFR() : NLDB0::getDFF(),
            NLName("ff" + suffix));
        auto* q = SNLScalarNet::create(top, NLName("q" + suffix));
        if (resettable) {
          ff->getInstTerm(NLDB0::getDFFRClock())->setNet(clk);
          ff->getInstTerm(NLDB0::getDFFRReset())->setNet(rst);
          ff->getInstTerm(NLDB0::getDFFRData())->setNet(value);
          ff->getInstTerm(NLDB0::getDFFROutput())->setNet(q);
        } else {
          ff->getInstTerm(NLDB0::getDFFClock())->setNet(clk);
          ff->getInstTerm(NLDB0::getDFFData())->setNet(value);
          ff->getInstTerm(NLDB0::getDFFOutput())->setNet(q);
        }
        value = q;
      }
    }
    SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"))
        ->setNet(value);
    return top;
  }

  static C2RtlEquivalenceOptions delayed(size_t cycles = 1) {
    C2RtlEquivalenceOptions options;
    options.outputDelays.emplace("y", cycles);
    return options;
  }

  SNLDesign* ascendingBusDesign(const std::string& name, bool registeredZero) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* bus = SNLBusTerm::create(top, SNLTerm::Direction::Input, 0, 1,
                                  NLName("a"));
    auto* msb = SNLScalarNet::create(top, NLName("msb"));
    auto* lsb = SNLScalarNet::create(top, NLName("lsb"));
    bus->getBit(0)->setNet(msb);
    bus->getBit(1)->setNet(lsb);
    auto* value = msb;
    if (registeredZero) {
      auto* zero = SNLScalarNet::create(top, NLName("zero"));
      zero->setType(SNLNet::Type::Assign0);
      auto* clk = input(top, "clk");
      auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff"));
      value = SNLScalarNet::create(top, NLName("q"));
      ff->getInstTerm(NLDB0::getDFFClock())->setNet(clk);
      ff->getInstTerm(NLDB0::getDFFData())->setNet(zero);
      ff->getInstTerm(NLDB0::getDFFOutput())->setNet(value);
    }
    SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("y"))
        ->setNet(value);
    return top;
  }

  SNLScalarNet* registerValue(SNLDesign* top, SNLScalarNet* value,
                              SNLScalarNet* clk, const std::string& name) {
    auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName(name + "_ff"));
    auto* q = SNLScalarNet::create(top, NLName(name + "_q"));
    ff->getInstTerm(NLDB0::getDFFClock())->setNet(clk);
    ff->getInstTerm(NLDB0::getDFFData())->setNet(value);
    ff->getInstTerm(NLDB0::getDFFOutput())->setNet(q);
    return q;
  }

  SNLScalarNet* invert(SNLDesign* top, SNLScalarNet* value,
                       const std::string& name) {
    auto* result = SNLScalarNet::create(top, NLName(name));
    auto* inv = SNLInstance::create(top, inverter_, NLName(name + "_inv"));
    inv->getInstTerm(inverter_->getScalarTerm(NLName("A")))->setNet(value);
    inv->getInstTerm(inverter_->getScalarTerm(NLName("Y")))->setNet(result);
    return result;
  }

  SNLDesign* wideOutputDesign(const std::string& name, bool registered) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* value = input(top, "a");
    if (registered) {
      value = registerValue(top, value, input(top, "clk"), "data");
    }
    auto* lowBit = registered ? invert(top, value, "wrong_low_bit") : value;
    auto* output = SNLBusTerm::create(top, SNLTerm::Direction::Output,
                                     40, 0, NLName("mantissa"));
    for (int bit = 0; bit <= 40; ++bit) {
      output->getBit(bit)->setNet(bit == 0 ? lowBit : value);
    }
    return top;
  }

  SNLDesign* classifiedDesign(const std::string& name, bool registered) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* clk = registered ? input(top, "clk") : nullptr;
    for (const auto& output : {"nan", "inf", "mantissa", "exponent"}) {
      auto* value = input(top, (std::string(output) + "_in").c_str());
      if (registered) {
        if (std::string(output) == "exponent") {
          value = invert(top, value, "wrong_exponent");
        }
        value = registerValue(top, value, clk, output);
      }
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(output))
          ->setNet(value);
    }
    return top;
  }

  static C2RtlEquivalenceOptions eventual(size_t cycle,
                                        const std::string& condition = "true",
                                        const std::string& equality =
                                            "model.y == rtl.y") {
    C2RtlEquivalenceOptions options;
    options.eventuals.push_back({cycle, condition, equality});
    return options;
  }

  static C2RtlEquivalenceResult prove(SNLDesign* model, SNLDesign* rtl,
                                     C2RtlEquivalenceOptions options) {
    return C2RtlEquivalenceStrategy(model, rtl, Config::SolverType::KISSAT,
                                    std::move(options)).run(12);
  }

  NLLibrary* designs_ = nullptr;
  SNLDesign* inverter_ = nullptr;
};

TEST_F(C2RtlEquivalenceTests, LegacyOutputDelayStillProvesRegisteredOutput) {
  const auto result = prove(design("model"), design("rtl", Source::A, 2),
                            delayed(2));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.comparedOutputs, 1u);
  EXPECT_EQ(result.comparedBits, 1u);
  EXPECT_EQ(result.outputDelays,
            (std::vector<std::pair<std::string, size_t>>{{"y", 2}}));
}

TEST_F(C2RtlEquivalenceTests, LegacyIncorrectOutputDelayStillFindsMismatch) {
  const auto result = prove(design("model"), design("rtl", Source::A, 2),
                            delayed(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
  EXPECT_FALSE(result.counterexampleTrace.empty());
}

TEST_F(C2RtlEquivalenceTests, LegacyResettableOutputStillProves) {
  const auto result = prove(design("model"), design("rtl", Source::A, 1, true),
                            delayed());
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.reset, "rst[0]");
}

TEST_F(C2RtlEquivalenceTests, LegacyMissingOutputDelayRemainsUnsupported) {
  const auto result = prove(design("model"), design("rtl", Source::A, 1), {});
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("delay"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, ConstraintRestrictsProofToLegalInputDomain) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::Zero, 1);
  EXPECT_EQ(prove(model, rtl, delayed()).status, C2RtlEquivalenceStatus::Different);
  auto options = delayed();
  options.constraints = {"a == 0"};
  const auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, ConstraintDoesNotHideMismatchInsideLegalDomain) {
  auto options = delayed();
  options.constraints = {"a == 0"};
  const auto result = prove(design("model"),
                            design("rtl", Source::InvertedA, 1), options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
  EXPECT_FALSE(result.counterexampleTrace.empty());
}

TEST_F(C2RtlEquivalenceTests, MultipleConstraintsAreConjoined) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::Zero, 1);
  auto options = delayed();
  options.constraints = {"model.a == rtl.b", "b < 1"};
  const auto constrained = prove(model, rtl, options);
  EXPECT_EQ(constrained.status, C2RtlEquivalenceStatus::Equivalent)
      << constrained.reason;
  for (const auto& single : {"model.a == rtl.b", "b < 1"}) {
    options.constraints = {single};
    const auto result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
  }
}

TEST_F(C2RtlEquivalenceTests, ContradictoryConstraintsRejectVacuousProofByDefault) {
  auto options = delayed();
  ASSERT_TRUE(options.checkReachability);
  options.constraints = {"a == 0", "a != 0"};
  const auto result = prove(design("model"), design("rtl", Source::A, 1), options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("unsatisfiable"), std::string::npos);
  EXPECT_NE(result.reason.find("vacuous"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, ReachabilityCheckCanBeExplicitlyDisabled) {
  auto options = delayed();
  options.constraints = {"a == 0", "a != 0"};
  options.checkReachability = false;
  const auto result = prove(design("model"),
                            design("rtl", Source::InvertedA, 1), options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, InvalidEarlierInputsInvalidateTheEntireTracePrefix) {
  // An earlier a != b can persist in the pipelines after the inputs become
  // equal. A current-frame-only assumption would incorrectly report a mismatch.
  auto options = delayed(2);
  options.constraints = {"a == b"};
  const auto result = prove(design("model"), design("rtl", Source::B, 2), options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, ConstraintsCannotReferenceOutputsEvenWhenInactive) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  for (const auto& expression : {"model.y == 0", "rtl.y == 0",
                                  "false && model.y"}) {
    auto options = delayed();
    options.constraints = {expression};
    const auto result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
    EXPECT_NE(result.reason.find("only inputs"), std::string::npos);
  }
}

TEST_F(C2RtlEquivalenceTests, ConstraintReportsUnknownTerminalAndExpressionIndex) {
  auto options = delayed();
  options.constraints = {"true", "missing_input > 0"};
  const auto result = prove(design("model"), design("rtl", Source::A, 1), options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
  EXPECT_NE(result.reason.find("constraint[1]"), std::string::npos);
  EXPECT_NE(result.reason.find("missing_input"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, AscendingBusKeepsNumericSignificanceAndBitIndices) {
  auto* model = ascendingBusDesign("model", false);
  auto* rtl = ascendingBusDesign("rtl", true);
  for (const auto& expression : {"a < 2", "model.a < 2", "rtl.a < 2",
                                  "a[0] == 0"}) {
    auto options = delayed();
    options.constraints = {expression};
    const auto result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent)
        << expression << ": " << result.reason;
  }
  // In [0:1], index 1 is the low bit. Constraining it leaves the high bit free.
  auto options = delayed();
  options.constraints = {"a[1] == 0"};
  const auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualCycleZeroObservesInitialState) {
  const auto result = prove(design("model"), design("rtl", Source::A, 1),
                            eventual(0));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualCycleOneObservesFirstClockTransition) {
  const auto result = prove(design("model"), design("rtl", Source::A, 1),
                            eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualCycleFourUsesFixedCycleZeroInputs) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 4);
  const auto afterPropagation = prove(model, rtl, eventual(4));
  EXPECT_EQ(afterPropagation.status, C2RtlEquivalenceStatus::Equivalent)
      << afterPropagation.reason;
  const auto tooEarly = prove(model, rtl, eventual(3));
  EXPECT_EQ(tooEarly.status, C2RtlEquivalenceStatus::Different) << tooEarly.reason;
}

TEST_F(C2RtlEquivalenceTests, MultipleEventualsAtOneCycleAreAllRequired) {
  auto options = eventual(1);
  options.eventuals.push_back({1, "false", "model.y != rtl.y"});
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  EXPECT_EQ(prove(model, rtl, options).status, C2RtlEquivalenceStatus::Equivalent);
  options.eventuals.push_back({1, "true", "model.y != rtl.y"});
  const auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, FalseEventualConditionDisablesOnlyItsEquality) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::InvertedA, 1);
  const auto ignored = prove(model, rtl, eventual(1, "false"));
  EXPECT_EQ(ignored.status, C2RtlEquivalenceStatus::Equivalent) << ignored.reason;
  const auto required = prove(model, rtl, eventual(1, "true"));
  EXPECT_EQ(required.status, C2RtlEquivalenceStatus::Different) << required.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualConditionsCanUseInputsAndOutputsAtCycle) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::Zero, 1);
  for (const auto& condition : {"a == 0 && model.b == rtl.b",
                                "model.y == 0 && rtl.y == 0"}) {
    const auto result = prove(model, rtl, eventual(1, condition));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent)
        << condition << ": " << result.reason;
  }
  const auto result = prove(model, rtl, eventual(1, "rtl.y == 0 && rtl.a != 0"));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualEqualityRejectsInputTerminalsEvenIfInactive) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  for (const auto& equality : {"model.a == rtl.a", "a == 0",
                               "true || rtl.b"}) {
    const auto result = prove(model, rtl, eventual(1, "false", equality));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
    EXPECT_NE(result.reason.find("only outputs"), std::string::npos);
  }
}

TEST_F(C2RtlEquivalenceTests, EventualReportsUnknownOrUnqualifiedOutputs) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  for (const auto& equality : {"model.missing == rtl.y", "y == rtl.y"}) {
    const auto result = prove(model, rtl, eventual(1, "true", equality));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
  }
  const auto condition = prove(model, rtl, eventual(1, "rtl.missing"));
  EXPECT_EQ(condition.status, C2RtlEquivalenceStatus::Unsupported)
      << condition.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualCanCompareOutputsWithDifferentNames) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  rtl->getScalarTerm(NLName("y"))->setName(NLName("renamed"));
  const auto result = prove(model, rtl,
                            eventual(1, "true", "model.y == rtl.renamed"));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, EventualBit40EqualityLeavesOtherBitsUnconstrained) {
  auto* model = wideOutputDesign("model", false);
  auto* rtl = wideOutputDesign("rtl", true);
  const auto bit = prove(model, rtl, eventual(
      1, "true", "model.mantissa[40] == rtl.mantissa[40]"));
  EXPECT_EQ(bit.status, C2RtlEquivalenceStatus::Equivalent) << bit.reason;
  const auto bus = prove(model, rtl, eventual(
      1, "true", "model.mantissa == rtl.mantissa"));
  EXPECT_EQ(bus.status, C2RtlEquivalenceStatus::Different) << bus.reason;
}

TEST_F(C2RtlEquivalenceTests, ConditionalNanInfAndZeroIgnoreOnlyInactiveFields) {
  auto* model = classifiedDesign("model", false);
  auto* rtl = classifiedDesign("rtl", true);
  C2RtlEquivalenceOptions options;
  options.eventuals = {
      {1, "true", "model.nan == rtl.nan"},
      {1, "!model.nan", "model.inf == rtl.inf"},
      {1, "!model.nan && !model.inf", "model.mantissa == rtl.mantissa"},
      {1, "!model.nan && !model.inf && model.mantissa != 0",
       "model.exponent == rtl.exponent"}};
  for (const auto& domain : {"nan_in", "!nan_in && inf_in",
                             "!nan_in && !inf_in && !mantissa_in"}) {
    options.constraints = {domain};
    const auto result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent)
        << domain << ": " << result.reason;
  }
  options.constraints = {"!nan_in && !inf_in && mantissa_in"};
  const auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
}

}  // namespace
}  // namespace KEPLER_FORMAL::C2RTL
