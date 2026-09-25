// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "BoolExprCache.h"
#include "Config.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "Tree2BoolExpr.h"
#include "model/OpaquePolicy.h"
#include "model/SequentialDesignModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace {

using namespace naja::NL;
using namespace KEPLER_FORMAL::SEC;
using KEPLER_FORMAL::Config;

class OpaquePolicyTests : public ::testing::Test {
 protected:
  void SetUp() override {
    oldPolicy_ = Config::getErrorOnOpaque();
    Config::setErrorOnOpaque(false);
  }
  void TearDown() override {
    Config::setErrorOnOpaque(oldPolicy_);
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) universe->destroy();
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
  }

  SequentialDesignModel classifiedModel(
      ConnectivitySkipOrigin origin = ConnectivitySkipOrigin::OpaqueInternal) {
    SequentialDesignModel model;
    const SignalKey key{{1}, {2}};
    model.displayNameByKey.emplace(key, "island.cell.Q[0]");
    model.connectivitySkipInfoByKey.emplace(
        key, ConnectivitySkipInfo{origin, "missing primitive model"});
    return model;
  }

  void createLibraries() {
    auto* db = NLDB::create(NLUniverse::create());
    library_ = NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
    auto* primitives =
        NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("primitives"));
    opaque_ = SNLDesign::create(
        primitives, SNLDesign::Type::Primitive, NLName("UNMODELED"));
    SNLScalarTerm::create(opaque_, SNLTerm::Direction::Input, NLName("D"));
    SNLScalarTerm::create(opaque_, SNLTerm::Direction::Output, NLName("Q"));
  }

  SNLDesign* top(const std::string& name, SNLDesign* cell = nullptr,
                 bool connected = false) {
    auto* design = SNLDesign::create(library_, SNLDesign::Type::Standard, NLName(name));
    auto* data = SNLScalarNet::create(design, NLName("data"));
    auto* enable = SNLScalarNet::create(design, NLName("enable"));
    SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("a"))->setNet(data);
    SNLScalarTerm::create(design, SNLTerm::Direction::Input, NLName("e"))->setNet(enable);
    SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("good"))->setNet(data);
    if (cell) {
      auto* instance = SNLInstance::create(design, cell, NLName("unused_cell"));
      instance->getInstTerm(cell->getScalarTerm(NLName("D")))->setNet(data);
      if (auto* gate = cell->getScalarTerm(NLName("E"))) {
        instance->getInstTerm(gate)->setNet(enable);
      }
      auto* output = SNLScalarNet::create(design, NLName("cell_output"));
      instance->getInstTerm(cell->getScalarTerm(NLName("Q")))->setNet(output);
      if (connected) {
        SNLScalarTerm::create(design, SNLTerm::Direction::Output, NLName("bad"))->setNet(output);
      }
    }
    return design;
  }

  bool oldPolicy_ = false;
  NLLibrary* library_ = nullptr;
  SNLDesign* opaque_ = nullptr;
};

TEST_F(OpaquePolicyTests, DisabledPreservesClassificationAndReasons) {
  auto model = classifiedModel();
  applyOpaquePolicy(model, "top", 0);
  EXPECT_FALSE(model.hasUnsupportedFeatures());
  ASSERT_EQ(model.connectivitySkipInfoByKey.size(), 1u);
  EXPECT_EQ(model.connectivitySkipInfoByKey.begin()->second.detail,
            "missing primitive model");
}

TEST_F(OpaquePolicyTests, EnabledIdentifiesDesignSignalAndReason) {
  Config::setErrorOnOpaque(true);
  auto model = classifiedModel();
  applyOpaquePolicy(model, "candidate", 1);
  ASSERT_EQ(model.unsupportedReasons.size(), 1u);
  EXPECT_EQ(model.unsupportedReasons.front(),
            "SEC error-on-opaque: design 2 (`candidate`), signal "
            "`island.cell.Q[0]`: missing primitive model");
  applyOpaquePolicy(model, "candidate", 1);
  EXPECT_EQ(model.unsupportedReasons.size(), 1u);
}

TEST_F(OpaquePolicyTests, OtherConnectivitySkipKindsDoNotBecomeOpaque) {
  Config::setErrorOnOpaque(true);
  for (const auto origin : {ConnectivitySkipOrigin::NoDriver,
                            ConnectivitySkipOrigin::MultiDriver,
                            ConnectivitySkipOrigin::LogicalLoop,
                            ConnectivitySkipOrigin::MultiClockDomain,
                            ConnectivitySkipOrigin::UnknownConstant}) {
    auto model = classifiedModel(origin);
    applyOpaquePolicy(model, "top", 0);
    EXPECT_FALSE(model.hasUnsupportedFeatures());
  }
}

TEST_F(OpaquePolicyTests, DiagnosticOrderIsStableAndEarlierErrorsArePreserved) {
  Config::setErrorOnOpaque(true);
  auto model = classifiedModel();
  model.unsupportedReasons.push_back("earlier error");
  const SignalKey earlier{{3}, {4}};
  model.displayNameByKey.emplace(earlier, "a.Q[0]");
  model.connectivitySkipInfoByKey.emplace(
      earlier, ConnectivitySkipInfo{ConnectivitySkipOrigin::OpaqueInternal, "other reason"});
  applyOpaquePolicy(model, "top", 0);
  ASSERT_EQ(model.unsupportedReasons.size(), 2u);
  EXPECT_EQ(model.unsupportedReasons.front(), "earlier error");
  EXPECT_NE(model.unsupportedReasons.back().find("a.Q[0]"), std::string::npos);
}

TEST_F(OpaquePolicyTests, MissingDisplayNameHasSignalKeyFallback) {
  Config::setErrorOnOpaque(true);
  auto model = classifiedModel();
  model.displayNameByKey.clear();
  applyOpaquePolicy(model, "top", 0);
  ASSERT_EQ(model.unsupportedReasons.size(), 1u);
  EXPECT_NE(model.unsupportedReasons.front().find("signal `1.2.`"), std::string::npos);
}

TEST_F(OpaquePolicyTests, DefaultExtractionStillChecksSupportedOutput) {
  createLibraries();
  const auto model = SequentialDesignModel::extract(top("top", opaque_, true));
  EXPECT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_EQ(model.coveredObservedOutputCount(), 1u);
  EXPECT_EQ(model.totalObservedOutputCount(), 2u);
}

TEST_F(OpaquePolicyTests, EnabledExtractionRejectsDisconnectedNonLatchCell) {
  createLibraries();
  Config::setErrorOnOpaque(true);
  const auto model = SequentialDesignModel::extract(top("candidate", opaque_));
  ASSERT_TRUE(model.hasUnsupportedFeatures());
  EXPECT_NE(model.unsupportedReasons.front().find("unused_cell"), std::string::npos);
  EXPECT_NE(model.unsupportedReasons.front().find("UNMODELED"), std::string::npos);
  EXPECT_FALSE(naja::DNL::isCreated());
}

TEST_F(OpaquePolicyTests, EnabledExtractionRejectsDisconnectedLatch) {
  createLibraries();
  Config::setErrorOnOpaque(true);
  const auto model = SequentialDesignModel::extract(top("candidate", NLDB0::getDLatch()));
  ASSERT_TRUE(model.hasUnsupportedFeatures());
  EXPECT_NE(model.unsupportedReasons.front().find("latch"), std::string::npos);
}

TEST_F(OpaquePolicyTests, SupportedDesignIsUnchangedWhenEnabled) {
  createLibraries();
  Config::setErrorOnOpaque(true);
  const auto model = SequentialDesignModel::extract(top("supported"));
  EXPECT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_EQ(model.coveredObservedOutputCount(), 1u);
}

TEST_F(OpaquePolicyTests, RejectsOpacityInEitherDesign) {
  createLibraries();
  Config::setErrorOnOpaque(true);
  auto* supported = top("supported");
  auto* unsupported = top("unsupported", opaque_);
  for (const bool first : {true, false}) {
    const SequentialEquivalenceStrategy strategy(
        first ? unsupported : supported, first ? supported : unsupported,
        Config::KISSAT, SecEngine::Pdr, SecEncoding::Binary);
    const auto result = strategy.run(1);
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
    EXPECT_NE(result.reason.find(first ? "design 1" : "design 2"), std::string::npos);
  }
}

TEST_F(OpaquePolicyTests, StrictModeDetectsInputOnlyOpaqueLeafCells) {
  createLibraries();
  auto* sink = SNLDesign::create(opaque_->getLibrary(), SNLDesign::Type::Primitive,
                                NLName("OPAQUE_SINK"));
  auto* input = SNLScalarTerm::create(sink, SNLTerm::Direction::Input, NLName("D"));
  auto* design = top("with_sink");
  auto* instance = SNLInstance::create(design, sink, NLName("unknown_sink"));
  instance->getInstTerm(input)->setNet(design->getScalarTerm(NLName("a"))->getNet());
  EXPECT_FALSE(SequentialDesignModel::extract(design).hasUnsupportedFeatures());
  Config::setErrorOnOpaque(true);
  const auto model = SequentialDesignModel::extract(design);
  ASSERT_TRUE(model.hasUnsupportedFeatures());
  EXPECT_NE(model.unsupportedReasons.front().find("unknown_sink"), std::string::npos);
  EXPECT_NE(model.unsupportedReasons.front().find("outputless"), std::string::npos);
}

TEST_F(OpaquePolicyTests, OutputlessScanDoesNotTreatEmptyStandardHierarchyAsOpaque) {
  createLibraries();
  auto* design = top("empty_hierarchy");
  auto* empty = SNLDesign::create(library_, SNLDesign::Type::Standard, NLName("empty"));
  SNLInstance::create(design, empty, NLName("empty_child"));
  Config::setErrorOnOpaque(true);
  const auto model = SequentialDesignModel::extract(design);
  EXPECT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_EQ(model.coveredObservedOutputCount(), 1u);
}

TEST_F(OpaquePolicyTests, OutputlessScanHandlesBlackBoxesAndNoPortPrimitives) {
  createLibraries();
  for (const auto type : {SNLDesign::Type::UserBlackBox, SNLDesign::Type::Primitive}) {
    const auto suffix = type == SNLDesign::Type::Primitive ? "primitive" : "blackbox";
    auto* design = top(std::string("top_") + suffix);
    auto* unknown = SNLDesign::create(
        type == SNLDesign::Type::Primitive ? opaque_->getLibrary() : library_,
        type, NLName(suffix));
    SNLInstance::create(design, unknown, NLName("unknown"));
    Config::setErrorOnOpaque(true);
    const auto model = SequentialDesignModel::extract(design);
    ASSERT_TRUE(model.hasUnsupportedFeatures());
    EXPECT_NE(model.unsupportedReasons.front().find("unknown"), std::string::npos);
  }
}

}  // namespace
