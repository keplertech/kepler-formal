// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <limits>
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

  SNLScalarNet* binaryGate(SNLDesign* top, SNLScalarNet* a, SNLScalarNet* b,
                           uint64_t truthTable, const std::string& name) {
    auto* primitive = SNLDesign::create(inverter_->getLibrary(),
        SNLDesign::Type::Primitive, NLName(name + "_primitive"));
    auto* pinA = SNLScalarTerm::create(primitive, SNLTerm::Direction::Input,
                                      NLName("A"));
    auto* pinB = SNLScalarTerm::create(primitive, SNLTerm::Direction::Input,
                                      NLName("B"));
    auto* pinY = SNLScalarTerm::create(primitive, SNLTerm::Direction::Output,
                                      NLName("Y"));
    SNLDesignModeling::addCombinatorialArcs({pinA, pinB}, {pinY});
    SNLDesignModeling::setTruthTable(primitive,
        SNLTruthTable(2, truthTable, SNLTruthTable::fullDependencies(2)));
    auto* instance = SNLInstance::create(top, primitive, NLName(name));
    auto* value = SNLScalarNet::create(top, NLName(name + "_out"));
    instance->getInstTerm(pinA)->setNet(a);
    instance->getInstTerm(pinB)->setNet(b);
    instance->getInstTerm(pinY)->setNet(value);
    return value;
  }

  SNLDesign* wideOutputDesign(const std::string& name, bool registered,
                              int low = 0) {
    auto* top = SNLDesign::create(designs_, NLName(name));
    auto* value = input(top, "a");
    if (registered) {
      value = registerValue(top, value, input(top, "clk"), "data");
    }
    auto* lowBit = registered ? invert(top, value, "wrong_low_bit") : value;
    auto* output = SNLBusTerm::create(top, SNLTerm::Direction::Output,
                                     40, low, NLName("mantissa"));
    for (int bit = low; bit <= 40; ++bit) {
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

TEST_F(C2RtlEquivalenceTests, RejectsMissingDesignsAndIncorrectStatePlacement) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  for (const auto pair : {std::make_pair(static_cast<SNLDesign*>(nullptr), rtl),
                          std::make_pair(model, static_cast<SNLDesign*>(nullptr))}) {
    const auto result = prove(pair.first, pair.second, eventual(1));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
    EXPECT_NE(result.reason.find("requires both parsed RTL designs"),
              std::string::npos);
  }
  auto result = prove(rtl, model, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("purely combinational"), std::string::npos);
  result = prove(model, design("unclocked"), eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("requires clocked RTL state"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, InputMappingListsMissingAndUnexpectedTerminals) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  rtl->getScalarTerm(NLName("a"))->setName(NLName("alpha"));
  rtl->getScalarTerm(NLName("b"))->setName(NLName("beta"));
  const auto result = prove(model, rtl, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("missing in RTL: a, b"), std::string::npos);
  EXPECT_NE(result.reason.find("unmatched RTL inputs: alpha, beta"),
            std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, RejectsDifferentDataInputWidths) {
  auto* model = ascendingBusDesign("model", false);
  auto* rtl = design("rtl", Source::A, 1);
  rtl->getScalarTerm(NLName("b"))->destroy();
  const auto result = prove(model, rtl, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("input width mismatch for `a`"), std::string::npos);
  EXPECT_NE(result.reason.find("reference=2, RTL=1"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, LegacyOutputMappingRejectsNamesAndUnknownDelays) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  auto options = delayed();
  options.outputDelays.emplace("unknown_z", 1);
  options.outputDelays.emplace("unknown_a", 1);
  auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("unknown outputs: unknown_a, unknown_z"),
            std::string::npos);

  rtl->getScalarTerm(NLName("y"))->setName(NLName("actual_y"));
  result = prove(model, rtl, delayed());
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("missing in RTL: y"), std::string::npos);
  EXPECT_NE(result.reason.find("extra RTL outputs: actual_y"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, UnobservableOutputCannotProduceAnEquivalenceProof) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  rtl->getScalarTerm(NLName("y"))->setNet(
      SNLScalarNet::create(rtl, NLName("undriven")));
  const auto result = prove(model, rtl, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("Missing observed output expression"),
            std::string::npos) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, DifferentClockDomainsAreUnsupported) {
  auto* rtl = design("rtl", Source::A, 2);
  rtl->getInstance(NLName("ff1"))->getInstTerm(NLDB0::getDFFClock())
      ->setNet(input(rtl, "other_clock"));
  const auto result = prove(design("model"), rtl, eventual(2));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("complete coverage"), std::string::npos)
      << result.reason;
}

TEST_F(C2RtlEquivalenceTests, IndependentOutputClocksMustStillShareOneDomain) {
  auto* rtl = design("rtl", Source::A, 2);
  auto* first = rtl->getInstance(NLName("ff0"));
  auto* second = rtl->getInstance(NLName("ff1"));
  second->getInstTerm(NLDB0::getDFFClock())->setNet(input(rtl, "other_clock"));
  second->getInstTerm(NLDB0::getDFFData())
      ->setNet(rtl->getScalarTerm(NLName("b"))->getNet());
  SNLScalarTerm::create(rtl, SNLTerm::Direction::Output, NLName("first_y"))
      ->setNet(first->getInstTerm(NLDB0::getDFFOutput())->getNet());
  const auto result = prove(design("model"), rtl, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("one clock domain"), std::string::npos)
      << result.reason;
}

TEST_F(C2RtlEquivalenceTests, ClockBusAndNegativeInputIndicesAreRejected) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  auto* clock = rtl->getScalarTerm(NLName("clk"));
  auto* clockNet = clock->getNet();
  clock->destroy();
  auto* bus = SNLBusTerm::create(rtl, SNLTerm::Direction::Input, 1, 0,
                                NLName("clk"));
  bus->getBit(0)->setNet(clockNet);
  bus->getBit(1)->setNet(SNLScalarNet::create(rtl, NLName("unused_clock")));
  auto result = prove(model, rtl, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("clock must be a scalar"), std::string::npos)
      << result.reason;

  auto* offset = SNLBusTerm::create(model, SNLTerm::Direction::Input, -1, -2,
                                   NLName("negative_index"));
  for (const int bit : {-1, -2}) {
    offset->getBit(bit)->setNet(SNLScalarNet::create(model,
        NLName("negative" + std::to_string(-bit))));
  }
  result = prove(model, design("scalar_clock", Source::A, 1), eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("expected a non-negative bit index"),
            std::string::npos) << result.reason;
}

TEST_F(C2RtlEquivalenceTests, FrameLimitReportsInconclusiveRatherThanEquivalent) {
  const auto result = C2RtlEquivalenceStrategy(
      design("model"), design("rtl", Source::A, 4),
      Config::SolverType::KISSAT, eventual(4)).run(1);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Inconclusive) << result.reason;
  EXPECT_NE(result.reason.find("configured frame bound"), std::string::npos);
  EXPECT_TRUE(result.counterexampleTrace.empty());
}

TEST_F(C2RtlEquivalenceTests, CompoundDataExpressionsRemainEquivalentAfterDelay) {
  for (const uint64_t table : {0b1110u, 0b0110u}) {
    const auto suffix = std::to_string(table);
    auto* model = design("model" + suffix);
    auto* rtl = design("rtl" + suffix, Source::A, 1);
    for (auto* top : {model, rtl}) {
      auto* a = static_cast<SNLScalarNet*>(top->getScalarTerm(NLName("a"))->getNet());
      auto* b = static_cast<SNLScalarNet*>(top->getScalarTerm(NLName("b"))->getNet());
      auto* expression = binaryGate(top, a, b, table,
          top->getName().getString() + "_logic");
      if (top == model) {
        top->getScalarTerm(NLName("y"))->setNet(expression);
      } else {
        top->getInstance(NLName("ff0"))->getInstTerm(NLDB0::getDFFData())
            ->setNet(expression);
      }
    }
    const auto result = prove(model, rtl, delayed());
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  }
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

TEST_F(C2RtlEquivalenceTests, EventualRejectsAmbiguousTimingAndOverflowingCycle) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  auto mixed = eventual(1);
  mixed.outputDelays.emplace("y", 1);
  const auto ambiguous = prove(model, rtl, mixed);
  EXPECT_EQ(ambiguous.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(ambiguous.reason.find("cannot be combined"), std::string::npos);

  const auto overflow = prove(model, rtl,
      eventual(std::numeric_limits<size_t>::max()));
  EXPECT_EQ(overflow.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(overflow.reason.find("cycle is too large"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, MultipleCyclesFindEarliestFailedCycle) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 2);
  auto options = eventual(2);
  options.eventuals.push_back({3, "true", "model.y == rtl.y"});
  auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.comparedOutputs, 2u);
  EXPECT_EQ(result.comparedBits, 2u);

  options.eventuals.push_back({1, "true", "model.y == rtl.y"});
  result = prove(model, rtl, options);
  ASSERT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;
  EXPECT_NE(result.counterexampleTrace.find("first bad frame at cycle 1"),
            std::string::npos);
  EXPECT_NE(result.counterexampleTrace.find("Failed eventual checks at cycle 1"),
            std::string::npos);
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

TEST_F(C2RtlEquivalenceTests, ClockTerminalUsesTheActiveEdgeValue) {
  auto* model = design("model");
  for (const bool negative : {false, true}) {
    auto* rtl = design(negative ? "rtl_negative" : "rtl_positive",
                       Source::InvertedA, 1);
    if (negative) {
      auto* clk = rtl->getScalarTerm(NLName("clk"))->getNet();
      rtl->getInstance(NLName("ff0"))->getInstTerm(NLDB0::getDFFClock())
          ->setNet(invert(rtl, static_cast<SNLScalarNet*>(clk), "clock_inverted"));
    }
    const auto inactive = negative ? "rtl.clk" : "!clk";
    auto result = prove(model, rtl, eventual(1, inactive));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
    const auto active = negative ? "!rtl.clk" : "clk";
    result = prove(model, rtl, eventual(1, active));
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Different) << result.reason;

    auto options = eventual(1);
    options.constraints = {inactive};
    result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported) << result.reason;
    EXPECT_NE(result.reason.find("unsatisfiable"), std::string::npos);
  }
}

TEST_F(C2RtlEquivalenceTests, ResetConstraintAppliesBeforeTheEventualCycle) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1, true);
  // A reset at cycle 0 can clear the result even if reset is low at cycle 1.
  auto options = eventual(1, "!rtl.rst");
  EXPECT_EQ(prove(model, rtl, options).status, C2RtlEquivalenceStatus::Different);
  options.constraints = {"!rst"};
  const auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_EQ(result.reset, "rst[0]");
  EXPECT_TRUE(result.resetActiveHigh);
}

TEST_F(C2RtlEquivalenceTests, InvertedResetKeepsItsExternalPolarity) {
  auto* rtl = design("rtl", Source::A, 1, true);
  auto* reset = static_cast<SNLScalarNet*>(rtl->getScalarTerm(NLName("rst"))->getNet());
  rtl->getInstance(NLName("ff0"))->getInstTerm(NLDB0::getDFFRReset())
      ->setNet(invert(rtl, reset, "active_low_reset"));
  auto options = eventual(1);
  options.constraints = {"rst"};
  const auto result = prove(design("model"), rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Equivalent) << result.reason;
  EXPECT_FALSE(result.resetActiveHigh);
}

TEST_F(C2RtlEquivalenceTests, CompoundResetAndDifferentResetSourcesAreRejected) {
  auto* model = design("model");
  auto* compound = design("compound", Source::A, 1, true);
  auto* reset = static_cast<SNLScalarNet*>(compound->getScalarTerm(NLName("rst"))->getNet());
  auto* data = static_cast<SNLScalarNet*>(compound->getScalarTerm(NLName("a"))->getNet());
  compound->getInstance(NLName("ff0"))->getInstTerm(NLDB0::getDFFRReset())
      ->setNet(binaryGate(compound, reset, data, 0b1110, "compound_reset"));
  auto result = prove(model, compound, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("unsupported combinational logic"), std::string::npos);

  auto* separate = design("separate", Source::A, 2, true);
  separate->getInstance(NLName("ff1"))->getInstTerm(NLDB0::getDFFRReset())
      ->setNet(input(separate, "other_reset"));
  result = prove(model, separate, eventual(2));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("one reset source"), std::string::npos);
}

TEST_F(C2RtlEquivalenceTests, ConstantAndPartiallyResetPipelinesAreRejected) {
  auto* model = design("model");
  auto* tied = design("tied", Source::A, 1, true);
  auto* zero = SNLScalarNet::create(tied, NLName("zero"));
  zero->setType(SNLNet::Type::Assign0);
  tied->getInstance(NLName("ff0"))->getInstTerm(NLDB0::getDFFRReset())
      ->setNet(zero);
  auto result = prove(model, tied, eventual(1));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("constant reset"), std::string::npos)
      << result.reason;

  auto* partial = design("partial", Source::A, 1, true);
  auto* output = partial->getScalarTerm(NLName("y"));
  output->setNet(registerValue(partial,
      static_cast<SNLScalarNet*>(output->getNet()),
      static_cast<SNLScalarNet*>(partial->getScalarTerm(NLName("clk"))->getNet()),
      "resetless"));
  result = prove(model, partial, eventual(2));
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("reset that covers only part"), std::string::npos)
      << result.reason;
}

TEST_F(C2RtlEquivalenceTests, InvalidBitSelectionsIdentifyTheConfigurationField) {
  auto* model = design("model");
  auto* rtl = design("rtl", Source::A, 1);
  auto options = eventual(1);
  options.constraints = {"a[1] == 0"};
  auto result = prove(model, rtl, options);
  EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("constraint[0]"), std::string::npos);
  EXPECT_NE(result.reason.find("bit index out of range"), std::string::npos);

  for (const bool condition : {false, true}) {
    options = eventual(1);
    auto& field = condition ? options.eventuals[0].condition
                            : options.eventuals[0].equality;
    field = condition ? "rtl.b[1]" : "model.y[1] == rtl.y";
    result = prove(model, rtl, options);
    EXPECT_EQ(result.status, C2RtlEquivalenceStatus::Unsupported);
    EXPECT_NE(result.reason.find(condition ? "eventual[0].condition"
                                          : "eventual[0].equality"),
              std::string::npos);
    EXPECT_NE(result.reason.find("bit index out of range"), std::string::npos);
  }
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

TEST_F(C2RtlEquivalenceTests, NonZeroBasedOutputSupportsExplicitBitSelection) {
  auto* model = wideOutputDesign("model", false, 39);
  auto* rtl = wideOutputDesign("rtl", true, 39);
  const auto bit = prove(model, rtl, eventual(
      1, "true", "model.mantissa[40] == rtl.mantissa[40]"));
  EXPECT_EQ(bit.status, C2RtlEquivalenceStatus::Equivalent) << bit.reason;
  const auto bus = prove(model, rtl, eventual(
      1, "true", "model.mantissa == rtl.mantissa"));
  EXPECT_EQ(bus.status, C2RtlEquivalenceStatus::Unsupported);
  EXPECT_NE(bus.reason.find("contiguous zero-based expression terminal"),
            std::string::npos);
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
