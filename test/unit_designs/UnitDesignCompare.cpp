// Copyright 2024-2025 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "gtest/gtest.h"

#include <filesystem>
#include <fstream>

#include "NLUniverse.h"

#include "SNLDesign.h"
#include "SNLScalarTerm.h"
#include "SNLScalarNet.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLTruthTable.h"
#include "SNLUtils.h"
#include "SNLVRLConstructor.h"
#include "MiterStrategy.h"

using namespace naja::NL;

#ifndef BENCHMARKS_PATH
#define BENCHMARKS_PATH "Undefined"
#endif

class UnitDesignCompare: public ::testing::Test {
  protected:
    void SetUp() override {
      NLUniverse* universe = NLUniverse::create();
      auto db = NLDB::create(universe);
      library0_ = NLLibrary::create(db, NLName("LIB0"));
      library1_ = NLLibrary::create(db, NLName("LIB1"));
    }
    void TearDown() override {
      NLUniverse::get()->destroy();
    }
  protected:
    NLLibrary*  library0_;
    NLLibrary*  library1_;
};

TEST_F(UnitDesignCompare, testSameDesigns) {
  SNLVRLConstructor constructor0(library0_);
  std::filesystem::path benchmarksPath(BENCHMARKS_PATH);
  constructor0.construct(benchmarksPath/"simple0.v");
  auto top = SNLUtils::findTop(library0_);
  auto halfadder0 = library0_->getSNLDesign(NLName("halfadder"));
  ASSERT_NE(nullptr, halfadder0);
  auto sumXor0 = halfadder0->getInstance(NLName("sum_xor"));
  auto ttSum0 = SNLDesignModeling::getTruthTable(sumXor0->getModel());
  auto or2 = top->getInstance(NLName("cout_or"));
  printf("or2: %s\n", or2->getDescription().c_str());
  auto ttOr2 = SNLDesignModeling::getTruthTable(or2->getModel());
  printf("TT or2: %s\n", ttOr2.getString().c_str());
  ASSERT_NE(nullptr, or2);
  
  auto sum = top->getScalarTerm(NLName("sum"));
  ASSERT_NE(nullptr, sum);
  EXPECT_EQ(SNLTerm::Direction::Output, sum->getDirection());
  auto sumNet = sum->getNet();
  ASSERT_NE(nullptr, sumNet);
  EXPECT_EQ(2, sumNet->getComponents().size());

  auto cout = top->getScalarTerm(NLName("cout"));
  ASSERT_NE(nullptr, cout);
  EXPECT_EQ(SNLTerm::Direction::Output, cout->getDirection());
  auto coutNet = cout->getNet();
  ASSERT_NE(nullptr, coutNet);
  EXPECT_EQ(2, coutNet->getComponents().size());
  auto topClone = top->clone(NLName("topClone"));
  KEPLER_FORMAL::MiterStrategy miterS(top, topClone);
  miterS.init();
  EXPECT_TRUE(miterS.run());
}

TEST_F(UnitDesignCompare, testDifferentDesigns) {
  std::filesystem::path benchmarksPath(BENCHMARKS_PATH);
  SNLVRLConstructor constructor0(library0_);
  constructor0.construct(benchmarksPath/"simple0.v");
  auto top0 = SNLUtils::findTop(library0_);

  SNLVRLConstructor constructor1(library1_);
  constructor1.construct(benchmarksPath/"simple1.v");
  auto top1 = SNLUtils::findTop(library1_);

  auto halfadder0 = library0_->getSNLDesign(NLName("halfadder"));
  ASSERT_NE(nullptr, halfadder0);
  auto sumXor0 = halfadder0->getInstance(NLName("sum_xor"));
  ASSERT_NE(nullptr, sumXor0);
  auto ttSum0 = SNLDesignModeling::getTruthTable(sumXor0->getModel());
  printf("TT sum0: %s\n", ttSum0.getString().c_str());
  for (auto term : sumXor0->getModel()->getBitTerms()) {
    printf("Term in sum0: %s\n", term->getDescription().c_str());
    // print direction 
    printf(" Direction: %d\n", static_cast<int>(term->getDirection()));
    if (term->getDirection() == SNLBitTerm::Direction::Input) {
      continue;
    }
    auto ttSum0 = SNLDesignModeling::getTruthTable(sumXor0->getModel(), term->getOrderID());
  }

  auto halfadder1 = library1_->getSNLDesign(NLName("halfadder"));
  ASSERT_NE(nullptr, halfadder1);
  auto sumXor1 = halfadder1->getInstance(NLName("sum_xor"));
  ASSERT_NE(nullptr, sumXor1);
  auto ttSum1 = SNLDesignModeling::getTruthTable(sumXor1->getModel());
  printf("TT sum1: %s\n", ttSum1.getString().c_str());

  EXPECT_TRUE(ttSum0.isGeneric());
  EXPECT_EQ(ttSum0.getGenericType(), SNLTruthTable::GenericType::XOR);
  EXPECT_TRUE(ttSum1.isGeneric());
  EXPECT_EQ(ttSum1.getGenericType(), SNLTruthTable::GenericType::OR);

  const auto logPath = std::filesystem::current_path() / "logB";
  std::error_code ec;
  std::filesystem::remove(logPath, ec);
  KEPLER_FORMAL::MiterStrategy miterS(top0, top1, logPath.string());
  miterS.init();
  EXPECT_FALSE(miterS.run());
  //should be different
  //here the issue comes from missing truth table but nothing is reported
}

TEST_F(UnitDesignCompare, testDiffWithConstants0) {
  std::filesystem::path benchmarksPath(BENCHMARKS_PATH);
  SNLVRLConstructor constructor0(library0_);
  constructor0.construct(benchmarksPath/"simple1.v");
  auto top0 = SNLUtils::findTop(library0_);

  SNLVRLConstructor constructor1(library1_);
  constructor1.construct(benchmarksPath/"simple2.v");
  auto top1 = SNLUtils::findTop(library1_);

  KEPLER_FORMAL::MiterStrategy miterS(top0, top1);
  miterS.init();
  EXPECT_FALSE(miterS.run());
}

TEST_F(UnitDesignCompare, testDiffWithConstan1) {
  std::filesystem::path benchmarksPath(BENCHMARKS_PATH);
  SNLVRLConstructor constructor0(library0_);
  constructor0.construct(benchmarksPath/"simple1.v");
  auto top0 = SNLUtils::findTop(library0_);

  SNLVRLConstructor constructor1(library1_);
  constructor1.construct(benchmarksPath/"simple3.v");
  auto top1 = SNLUtils::findTop(library1_);

  KEPLER_FORMAL::MiterStrategy miterS(top0, top1);
  miterS.init();
  EXPECT_FALSE(miterS.run());
}

TEST_F(UnitDesignCompare, testManualAndDesigns) {
  auto primitives = NLLibrary::create(
      library0_->getDB(), NLLibrary::Type::Primitives, NLName("PRIMITIVES"));
  const char* andName = "AND";
  auto andModel =
    SNLDesign::create(primitives, SNLDesign::Type::Primitive, NLName(andName));
  auto andIn1 = SNLScalarTerm::create(andModel, SNLTerm::Direction::Input,
                                        NLName("a"));
  auto andIn2 = SNLScalarTerm::create(andModel, SNLTerm::Direction::Input,
                                        NLName("b"));
  auto andOut = SNLScalarTerm::create(andModel, SNLTerm::Direction::Output,
                                        NLName("y"));

  auto makeAndTop = [&](NLLibrary* library, const char* topName,
                       const char* andName) -> SNLDesign* {

    auto top =
        SNLDesign::create(library, SNLDesign::Type::Standard, NLName(topName));
    auto topIn1 = SNLScalarTerm::create(top, SNLTerm::Direction::Input,
                                        NLName("in1"));
    auto topIn2 = SNLScalarTerm::create(top, SNLTerm::Direction::Input,
                                        NLName("in2"));
    auto topOut = SNLScalarTerm::create(top, SNLTerm::Direction::Output,
                                        NLName("out"));

    auto andInst = SNLInstance::create(top, andModel, NLName("and0"));
    auto netIn1 = SNLScalarNet::create(top, NLName("net_in1"));
    auto netIn2 = SNLScalarNet::create(top, NLName("net_in2"));
    auto netOut = SNLScalarNet::create(top, NLName("net_out"));

    topIn1->setNet(netIn1);
    andInst->getInstTerm(andIn1)->setNet(netIn1);
    topIn2->setNet(netIn2);
    andInst->getInstTerm(andIn2)->setNet(netIn2);
    topOut->setNet(netOut);
    andInst->getInstTerm(andOut)->setNet(netOut);

    return top;
  };

  auto top0 = makeAndTop(library0_, "top0", "AND");
  auto top1 = makeAndTop(library1_, "top1", "AND");

  ASSERT_NE(nullptr, top0);
  ASSERT_NE(nullptr, top1);
  EXPECT_EQ(1, top0->getInstances().size());
  EXPECT_EQ(1, top1->getInstances().size());

  {
    //should be the same designs but no truth table
    KEPLER_FORMAL::MiterStrategy miterS(top0, top1);
    miterS.init();
    EXPECT_TRUE(miterS.run());
  }

  SNLDesignModeling::setTruthTable(andModel, SNLTruthTable(2, 8, SNLTruthTable::fullDependencies(2)));
  EXPECT_TRUE(SNLDesignModeling::hasModeling(andModel));

  {
    //should be the same designs with same truth table
    KEPLER_FORMAL::MiterStrategy miterS(top0, top1);
    miterS.init();
    EXPECT_TRUE(miterS.run());
  }
}
