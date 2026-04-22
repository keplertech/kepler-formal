// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <array>
#include <cstdlib>
#include <optional>
#include <string>
#include <unordered_set>
#include <unistd.h>

#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLBusNet.h"
#include "SNLBusNetBit.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "common/BoolExprUtils.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "imc/IMCEngine.h"
#include "kinduction/KInductionEngine.h"
#include "pdr/PDREngine.h"
#include "kinduction/BaseCaseSolver.h"
#include "kinduction/SatEncoding.h"
#include "kinduction/InductionStepSolver.h"
#include "model/SequentialDesignModel.h"
#include "strategy/ReachableStateInvariant.h"
#include "strategy/SequentialEquivalenceStrategy.h"
#include "strategy/StructuralStateInvariant.h"

using namespace naja::NL;
using namespace KEPLER_FORMAL::SEC;
using KEPLER_FORMAL::BoolExpr;

namespace KEPLER_FORMAL::SEC::detail {

BoolExpr* buildNextStateExprForTest(
    size_t stateTermID,
    const std::unordered_map<std::string, naja::DNL::DNLID>& pinTermIDs,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm);

std::optional<bool> detectInitialStateValueForTest(
    const std::unordered_map<std::string, naja::DNL::DNLID>& pinTermIDs);

std::optional<bool> evaluateConstantUnderAssignmentsForTest(
    BoolExpr* expr,
    const std::unordered_map<size_t, bool>& assignments);

void inferSynthesizedResetInitialStateValuesForTest(SequentialDesignModel& model);

std::optional<bool> getResetAssertionValueForTest(const std::string& displayName);

std::unordered_map<SignalKey, bool, SignalKeyHash>
deriveResetBootstrapStateValuesForTest(
    const SequentialDesignModel& model,
    size_t cycles);

AlignedSignals filterStateEqualitiesByInitialValueForTest(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& candidateStates);

std::string formatStringListForTest(const std::vector<std::string>& values,
                                    size_t limit);

std::string formatConeLevelsForTest(
    const std::vector<std::vector<std::string>>& levels);

std::string formatCounterexampleWitnessForTest(
    const KInductionResult& result,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1);

AlignedSignals alignSignalsByNameForTest(
    const std::vector<SignalKey>& keys0,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames0,
    const std::vector<SignalKey>& keys1,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames1,
    const char* label);

SignalKey getTerminalPathKeyForTest(
    const naja::DNL::DNLTerminalFull& terminal);

std::optional<naja::DNL::DNLID> findTermByKeyForTest(
    naja::DNL::DNLFull* dnl,
    const SignalKey& key);

std::vector<naja::DNL::DNLID> selectRequiredBuilderOutputsForTest(
    const std::vector<naja::DNL::DNLID>& collectedOutputs,
    const std::unordered_set<naja::DNL::DNLID>& topOutputTerms,
    const std::vector<naja::DNL::DNLID>& sequentialDependencyTerms,
    const std::unordered_set<naja::DNL::DNLID>& prunedBuilderOutputTerms);

}  // namespace KEPLER_FORMAL::SEC::detail

namespace {

class SequentialEquivalenceStrategyTests : public ::testing::Test {
 protected:
  void TearDown() override {
    naja::DNL::destroy();
    if (auto* universe = NLUniverse::get()) {
      universe->destroy();
    }
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
};

class ScopedEnvVar {
 public:
  ScopedEnvVar(const char* name, const char* value)
      : name_(name) {
    if (const char* current = std::getenv(name_); current != nullptr) {
      previousValue_ = current;
    }
    setenv(name_, value, 1);
  }

  ~ScopedEnvVar() {
    if (previousValue_.has_value()) {
      setenv(name_, previousValue_->c_str(), 1);
    } else {
      unsetenv(name_);
    }
  }

 private:
  const char* name_;
  std::optional<std::string> previousValue_;
};

class ScopedSecBoundaryAbstraction {
 public:
  explicit ScopedSecBoundaryAbstraction(bool enabled)
      : previousValue_(
            KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary()) {
    KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(enabled);
  }

  ~ScopedSecBoundaryAbstraction() {
    KEPLER_FORMAL::Config::setSecTreatUncomputableSeqAsBoundary(previousValue_);
  }

 private:
  bool previousValue_;
};

class ScopedFdCapture {
 public:
  explicit ScopedFdCapture(int fd)
      : fd_(fd) {
    if (::pipe(pipeFds_) != 0) {
      throw std::runtime_error("Failed to create capture pipe");
    }
    savedFd_ = ::dup(fd_);
    if (savedFd_ < 0) {
      ::close(pipeFds_[0]);
      ::close(pipeFds_[1]);
      throw std::runtime_error("Failed to duplicate capture fd");
    }
    if (::dup2(pipeFds_[1], fd_) < 0) {
      ::close(savedFd_);
      ::close(pipeFds_[0]);
      ::close(pipeFds_[1]);
      throw std::runtime_error("Failed to redirect capture fd");
    }
    ::close(pipeFds_[1]);
    pipeFds_[1] = -1;
  }

  ~ScopedFdCapture() {
    if (savedFd_ >= 0) {
      finish();
    }
  }

  std::string finish() {
    if (savedFd_ < 0) {
      return captured_;
    }
    ::dup2(savedFd_, fd_);
    ::close(savedFd_);
    savedFd_ = -1;

    char buffer[4096];
    while (true) {
      const ssize_t readSize = ::read(pipeFds_[0], buffer, sizeof(buffer));
      if (readSize <= 0) {
        break;
      }
      captured_.append(buffer, static_cast<size_t>(readSize));
    }
    ::close(pipeFds_[0]);
    pipeFds_[0] = -1;
    return captured_;
  }

 private:
  int fd_;
  int savedFd_ = -1;
  int pipeFds_[2] = {-1, -1};
  std::string captured_;
};

SignalKey makeSignalKey(const std::string& name) {
  SignalKey key;
  key.first.push_back(stableSignalKeyNameID(name));
  key.second.push_back(stableSignalKeyNameID(name));
  return key;
}

KInductionProblem buildLinearChainSecProblem(size_t logicalStateCount) {
  const auto bitCount = [logicalStateCount]() {
    size_t bits = 0;
    size_t encodedStates = 1;
    while (encodedStates < logicalStateCount) {
      encodedStates <<= 1;
      ++bits;
    }
    return std::max<size_t>(bits, 1);
  }();

  const auto buildStateExpr = [](const std::vector<size_t>& symbols, size_t value) {
    BoolExpr* expr = BoolExpr::createTrue();
    for (size_t bit = 0; bit < symbols.size(); ++bit) {
      expr = BoolExpr::And(
          expr,
          (value & (size_t{1} << bit)) != 0 ? BoolExpr::Var(symbols[bit])
                                            : BoolExpr::Not(BoolExpr::Var(symbols[bit])));
    }
    return BoolExpr::simplify(expr);
  };

  const auto buildNextBitExpr =
      [&](const std::vector<size_t>& symbols, size_t bitIndex) {
        BoolExpr* expr = BoolExpr::createFalse();
        for (size_t logicalState = 0; logicalState < logicalStateCount; ++logicalState) {
          const size_t nextLogicalState =
              logicalState + 1 < logicalStateCount ? logicalState + 1 : logicalState;
          if ((nextLogicalState & (size_t{1} << bitIndex)) == 0) {
            continue;
          }
          expr = BoolExpr::Or(expr, buildStateExpr(symbols, logicalState));
        }
        return BoolExpr::simplify(expr);
      };

  KInductionProblem problem;
  problem.state0Symbols.reserve(bitCount);
  problem.state1Symbols.reserve(bitCount);
  problem.allSymbols.reserve(bitCount * 2);

  size_t nextSymbol = 2;
  for (size_t bit = 0; bit < bitCount; ++bit) {
    problem.state0Symbols.push_back(nextSymbol++);
  }
  for (size_t bit = 0; bit < bitCount; ++bit) {
    problem.state1Symbols.push_back(nextSymbol++);
  }
  problem.allSymbols.insert(
      problem.allSymbols.end(), problem.state0Symbols.begin(), problem.state0Symbols.end());
  problem.allSymbols.insert(
      problem.allSymbols.end(), problem.state1Symbols.begin(), problem.state1Symbols.end());

  for (size_t bit = 0; bit < bitCount; ++bit) {
    problem.transitions0.emplace_back(
        problem.state0Symbols[bit], buildNextBitExpr(problem.state0Symbols, bit));
    problem.transitions1.emplace_back(
        problem.state1Symbols[bit], buildNextBitExpr(problem.state1Symbols, bit));
  }

  problem.initialCondition = BoolExpr::And(
      buildStateExpr(problem.state0Symbols, 0), buildStateExpr(problem.state1Symbols, 0));
  problem.initializedStateCount = problem.allSymbols.size();
  problem.totalStateCount = problem.allSymbols.size();
  problem.observedOutputExprs0 = {
      buildStateExpr(problem.state0Symbols, logicalStateCount - 1)};
  problem.observedOutputExprs1 = {
      buildStateExpr(problem.state1Symbols, logicalStateCount - 1)};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0.front(), problem.observedOutputExprs1.front());
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

SNLDesignModeling::BitTerms collectBitTerms(SNLBusTerm* bus) {
  SNLDesignModeling::BitTerms bits;
  for (auto* bit : bus->getBits()) {
    bits.push_back(bit);
  }
  return bits;
}

SNLDesign* createInvModel(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("INV"));
  auto* input =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  SNLDesignModeling::addCombinatorialArcs({input}, {output});
  SNLDesignModeling::setTruthTable(model, SNLTruthTable::Inv());
  return model;
}

SNLDesign* createAnd2Model(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("AND2"));
  auto* input0 =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* input1 =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("B"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  SNLDesignModeling::addCombinatorialArcs({input0, input1}, {output});
  SNLDesignModeling::setTruthTable(
      model,
      SNLTruthTable(
          2,
          SNLTruthTable::GenericType::AND,
          SNLTruthTable::fullDependencies(2)));
  return model;
}

SNLDesign* createOr2Model(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("OR2"));
  auto* input0 =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* input1 =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("B"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  SNLDesignModeling::addCombinatorialArcs({input0, input1}, {output});
  SNLDesignModeling::setTruthTable(
      model,
      SNLTruthTable(
          2,
          SNLTruthTable::GenericType::OR,
          SNLTruthTable::fullDependencies(2)));
  return model;
}

SNLDesign* createDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    bool invertData,
    bool invertOutput,
    const std::string& inputName,
    const std::string& outputName,
    const std::string& ffName) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName(inputName));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName(outputName));

  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName(ffName));
  SNLInstance* dataInv = nullptr;
  SNLInstance* outputInv = nullptr;
  if (invertData) {
    dataInv = SNLInstance::create(top, invModel, NLName("inv_data"));
  }
  if (invertOutput) {
    outputInv = SNLInstance::create(top, invModel, NLName("inv_out"));
  }

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netData = SNLScalarNet::create(top, NLName("net_data"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);

  if (invertData) {
    dataInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netIn);
    dataInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netData);
  }

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(invertData ? netData : netIn);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  if (invertOutput) {
    topOut->setNet(netOut);
    outputInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netQ);
    outputInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netOut);
  } else {
    topOut->setNet(netQ);
  }

  return top;
}

SNLDesign* createDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    bool invertData,
    bool invertOutput,
    const std::string& ffName = "ff0") {
  return createDffTop(
      library, name, invModel, invertData, invertOutput, "in", "out", ffName);
}

SNLDesign* createNamedOutputDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    const std::string& outputName) {
  return createDffTop(
      library, name, invModel, false, false, "in", outputName, "ff0");
}

SNLDesign* createNamedInputDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    const std::string& inputName) {
  return createDffTop(
      library, name, invModel, false, false, inputName, "out", "ff0");
}

SNLDesign* createExtraOutputDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));
  auto* topExtra =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out_extra"));

  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  auto* extraInv = SNLInstance::create(top, invModel, NLName("inv_extra"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));
  auto* netExtra = SNLScalarNet::create(top, NLName("net_extra"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topOut->setNet(netQ);
  topExtra->setNet(netExtra);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);
  extraInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netQ);
  extraInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netExtra);

  return top;
}

SNLDesign* createExtraInputDffTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topExtra =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in_extra"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netExtra = SNLScalarNet::create(top, NLName("net_extra"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topExtra->setNet(netExtra);
  topClock->setNet(netClock);
  topOut->setNet(netQ);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createDffeTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topEnable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("en"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* ff = SNLInstance::create(top, NLDB0::getDFFE(), NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netEnable = SNLScalarNet::create(top, NLName("net_en"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topEnable->setNet(netEnable);
  topClock->setNet(netClock);
  topOut->setNet(netQ);

  ff->getInstTerm(NLDB0::getDFFEClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFEData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFEEnable())->setNet(netEnable);
  ff->getInstTerm(NLDB0::getDFFEOutput())->setNet(netQ);

  return top;
}

SNLDesign* createResetInitializedPipelineTop(
    NLLibrary* library,
    const std::string& name,
    bool driveLastStageFromReset,
    const std::vector<std::string>& ffNames);

SNLDesign* createResetInitializedShiftPipelineTopWithStages(
    NLLibrary* library,
    const std::string& name,
    size_t stages);

SNLDesign* createResetInitializedPipelineTop(
    NLLibrary* library,
    const std::string& name,
    bool driveLastStageFromReset) {
  return createResetInitializedPipelineTop(
      library,
      name,
      driveLastStageFromReset,
      {"ff0", "ff1", "ff2"});
}

SNLDesign* createResetInitializedPipelineTop(
    NLLibrary* library,
    const std::string& name,
    bool driveLastStageFromReset,
    const std::vector<std::string>& ffNames) {
  if (ffNames.size() != 3) {
    throw std::invalid_argument(
        "createResetInitializedPipelineTop expects exactly three flop names");
  }

  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topResetN =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst_n"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* ff0 = SNLInstance::create(top, NLDB0::getDFFRN(), NLName(ffNames[0]));
  auto* ff1 = SNLInstance::create(top, NLDB0::getDFFRN(), NLName(ffNames[1]));
  auto* ff2 = SNLInstance::create(top, NLDB0::getDFFRN(), NLName(ffNames[2]));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netQ0 = SNLScalarNet::create(top, NLName("net_q0"));
  auto* netQ1 = SNLScalarNet::create(top, NLName("net_q1"));
  auto* netQ2 = SNLScalarNet::create(top, NLName("net_q2"));

  topIn->setNet(netIn);
  topResetN->setNet(netResetN);
  topClock->setNet(netClock);
  topOut->setNet(netQ0);

  for (auto* ff : {ff0, ff1, ff2}) {
    ff->getInstTerm(NLDB0::getDFFRNClock())->setNet(netClock);
    ff->getInstTerm(NLDB0::getDFFRNResetN())->setNet(netResetN);
  }

  ff0->getInstTerm(NLDB0::getDFFRNData())->setNet(netQ1);
  ff0->getInstTerm(NLDB0::getDFFRNOutput())->setNet(netQ0);
  ff1->getInstTerm(NLDB0::getDFFRNData())->setNet(netQ2);
  ff1->getInstTerm(NLDB0::getDFFRNOutput())->setNet(netQ1);
  ff2->getInstTerm(NLDB0::getDFFRNData())->setNet(
      driveLastStageFromReset ? netResetN : netIn);
  ff2->getInstTerm(NLDB0::getDFFRNOutput())->setNet(netQ2);

  return top;
}

SNLDesign* createResetInitializedShiftPipelineTopWithStages(
    NLLibrary* library,
    const std::string& name,
    size_t stages) {
  if (stages == 0) {
    throw std::invalid_argument(
        "createResetInitializedShiftPipelineTopWithStages expects at least one stage");
  }

  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topResetN =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst_n"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  std::vector<SNLScalarNet*> stageNets;
  stageNets.reserve(stages);
  for (size_t i = 0; i < stages; ++i) {
    stageNets.push_back(
        SNLScalarNet::create(top, NLName("net_q" + std::to_string(i))));
  }

  topIn->setNet(netIn);
  topResetN->setNet(netResetN);
  topClock->setNet(netClock);
  topOut->setNet(stageNets.front());

  for (size_t i = 0; i < stages; ++i) {
    auto* ff = SNLInstance::create(
        top, NLDB0::getDFFRN(), NLName("ff" + std::to_string(i)));
    ff->getInstTerm(NLDB0::getDFFRNClock())->setNet(netClock);
    ff->getInstTerm(NLDB0::getDFFRNResetN())->setNet(netResetN);
    ff->getInstTerm(NLDB0::getDFFRNData())->setNet(
        i + 1 == stages ? netIn : stageNets[i + 1]);
    ff->getInstTerm(NLDB0::getDFFRNOutput())->setNet(stageNets[i]);
  }

  return top;
}

SNLDesign* createBootstrapPipelineTopWithStages(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    SNLDesign* andModel,
    size_t stages) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* resetInv = SNLInstance::create(top, invModel, NLName("reset_inv"));
  std::vector<SNLInstance*> gates;
  std::vector<SNLInstance*> flops;
  gates.reserve(stages);
  flops.reserve(stages);
  for (size_t i = 0; i < stages; ++i) {
    gates.push_back(
        SNLInstance::create(top, andModel, NLName("gate" + std::to_string(i))));
    flops.push_back(
        SNLInstance::create(top, NLDB0::getDFF(), NLName("ff" + std::to_string(i))));
  }

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  std::vector<SNLScalarNet*> dataNets;
  std::vector<SNLScalarNet*> stateNets;
  dataNets.reserve(stages);
  stateNets.reserve(stages);
  for (size_t i = 0; i < stages; ++i) {
    dataNets.push_back(
        SNLScalarNet::create(top, NLName("net_d" + std::to_string(i))));
    stateNets.push_back(
        SNLScalarNet::create(top, NLName("net_q" + std::to_string(i))));
  }

  topIn->setNet(netIn);
  topReset->setNet(netReset);
  topClock->setNet(netClock);
  topOut->setNet(stateNets.front());

  resetInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netReset);
  resetInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netResetN);

  for (size_t i = 0; i < stages; ++i) {
    gates[i]->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(
        i + 1 == stages ? netIn : stateNets[i + 1]);
    gates[i]->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netResetN);
    gates[i]->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(dataNets[i]);
  }

  for (auto* ff : flops) {
    ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  }
  for (size_t i = 0; i < stages; ++i) {
    flops[i]->getInstTerm(NLDB0::getDFFData())->setNet(dataNets[i]);
    flops[i]->getInstTerm(NLDB0::getDFFOutput())->setNet(stateNets[i]);
  }

  return top;
}

SNLDesign* createBootstrapPipelineTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    SNLDesign* andModel) {
  return createBootstrapPipelineTopWithStages(library, name, invModel, andModel, 3);
}

SNLDesign* createResetLoadsInputTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    SNLDesign* andModel,
    SNLDesign* orModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* resetInv = SNLInstance::create(top, invModel, NLName("reset_inv"));
  auto* loadData = SNLInstance::create(top, andModel, NLName("load_data"));
  auto* holdData = SNLInstance::create(top, andModel, NLName("hold_data"));
  auto* muxOut = SNLInstance::create(top, orModel, NLName("mux_out"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLoad = SNLScalarNet::create(top, NLName("net_load"));
  auto* netHold = SNLScalarNet::create(top, NLName("net_hold"));
  auto* netD = SNLScalarNet::create(top, NLName("net_d"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topReset->setNet(netReset);
  topClock->setNet(netClock);
  topOut->setNet(netQ);

  resetInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netReset);
  resetInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netResetN);

  loadData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netReset);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netIn);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netLoad);

  holdData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netResetN);
  holdData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netQ);
  holdData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netHold);

  muxOut->getInstTerm(orModel->getScalarTerm(NLName("A")))->setNet(netLoad);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("B")))->setNet(netHold);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("Y")))->setNet(netD);

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netD);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createResetLoadsInputTwoStageTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    SNLDesign* andModel,
    SNLDesign* orModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* resetInv = SNLInstance::create(top, invModel, NLName("reset_inv"));
  auto* loadData = SNLInstance::create(top, andModel, NLName("load_data"));
  auto* holdData = SNLInstance::create(top, andModel, NLName("hold_data"));
  auto* muxOut = SNLInstance::create(top, orModel, NLName("mux_out"));
  auto* ffHidden = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff_hidden"));
  auto* ffOut = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff_out"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLoad = SNLScalarNet::create(top, NLName("net_load"));
  auto* netHold = SNLScalarNet::create(top, NLName("net_hold"));
  auto* netHiddenD = SNLScalarNet::create(top, NLName("net_hidden_d"));
  auto* netHiddenQ = SNLScalarNet::create(top, NLName("net_hidden_q"));
  auto* netOutQ = SNLScalarNet::create(top, NLName("net_out_q"));

  topIn->setNet(netIn);
  topReset->setNet(netReset);
  topClock->setNet(netClock);
  topOut->setNet(netOutQ);

  resetInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netReset);
  resetInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netResetN);

  loadData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netReset);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netIn);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netLoad);

  holdData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netResetN);
  holdData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netHiddenQ);
  holdData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netHold);

  muxOut->getInstTerm(orModel->getScalarTerm(NLName("A")))->setNet(netLoad);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("B")))->setNet(netHold);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("Y")))->setNet(netHiddenD);

  for (auto* ff : {ffHidden, ffOut}) {
    ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  }
  ffHidden->getInstTerm(NLDB0::getDFFData())->setNet(netHiddenD);
  ffHidden->getInstTerm(NLDB0::getDFFOutput())->setNet(netHiddenQ);
  ffOut->getInstTerm(NLDB0::getDFFData())->setNet(netHiddenQ);
  ffOut->getInstTerm(NLDB0::getDFFOutput())->setNet(netOutQ);

  return top;
}

SNLDesign* createResetLoadsInputShiftPipelineTopWithStages(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    SNLDesign* andModel,
    SNLDesign* orModel,
    size_t stages) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* resetInv = SNLInstance::create(top, invModel, NLName("reset_inv"));
  auto* loadData = SNLInstance::create(top, andModel, NLName("load_data"));
  auto* holdData = SNLInstance::create(top, andModel, NLName("hold_data"));
  auto* muxOut = SNLInstance::create(top, orModel, NLName("mux_out"));

  std::vector<SNLInstance*> flops;
  flops.reserve(stages);
  for (size_t i = 0; i < stages; ++i) {
    flops.push_back(
        SNLInstance::create(top, NLDB0::getDFF(), NLName("ff" + std::to_string(i))));
  }

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netResetN = SNLScalarNet::create(top, NLName("net_rst_n"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLoad = SNLScalarNet::create(top, NLName("net_load"));
  auto* netHold = SNLScalarNet::create(top, NLName("net_hold"));
  auto* netLastD = SNLScalarNet::create(top, NLName("net_last_d"));
  std::vector<SNLScalarNet*> stateNets;
  stateNets.reserve(stages);
  for (size_t i = 0; i < stages; ++i) {
    stateNets.push_back(
        SNLScalarNet::create(top, NLName("net_q" + std::to_string(i))));
  }

  topIn->setNet(netIn);
  topReset->setNet(netReset);
  topClock->setNet(netClock);
  topOut->setNet(stateNets.front());

  resetInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netReset);
  resetInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netResetN);

  loadData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netReset);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netIn);
  loadData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netLoad);

  holdData->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netResetN);
  holdData->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(stateNets.back());
  holdData->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netHold);

  muxOut->getInstTerm(orModel->getScalarTerm(NLName("A")))->setNet(netLoad);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("B")))->setNet(netHold);
  muxOut->getInstTerm(orModel->getScalarTerm(NLName("Y")))->setNet(netLastD);

  for (auto* ff : flops) {
    ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  }
  for (size_t i = 0; i + 1 < stages; ++i) {
    flops[i]->getInstTerm(NLDB0::getDFFData())->setNet(stateNets[i + 1]);
    flops[i]->getInstTerm(NLDB0::getDFFOutput())->setNet(stateNets[i]);
  }
  flops.back()->getInstTerm(NLDB0::getDFFData())->setNet(netLastD);
  flops.back()->getInstTerm(NLDB0::getDFFOutput())->setNet(stateNets.back());

  return top;
}

SNLDesign* createDffQnModel(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("DFF_Q_QN"));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* q =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  auto* qn =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("QN"));
  SNLDesignModeling::addInputsToClockArcs({data}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {q, qn});
  return model;
}

SNLDesign* createNamedComplementSequentialModel(
    NLLibrary* library,
    const std::string& name,
    const std::string& primaryPinName,
    const std::string& complementPinName) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* primary = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(primaryPinName));
  auto* complement = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(complementPinName));
  SNLDesignModeling::addInputsToClockArcs({data}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {primary, complement});
  return model;
}

SNLDesign* createComplementFirstSequentialModel(
    NLLibrary* library,
    const std::string& name,
    const std::string& primaryPinName,
    const std::string& complementPinName) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* complement = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(complementPinName));
  auto* primary = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(primaryPinName));
  SNLDesignModeling::addInputsToClockArcs({data}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {primary, complement});
  return model;
}

SNLDesign* createSetOnlySequentialModel(NLLibrary* library,
                                        const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* set =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("S"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  SNLDesignModeling::addInputsToClockArcs({data, set}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {output});
  return model;
}

SNLDesign* createBusSequentialModel(NLLibrary* library,
                                    const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data = SNLBusTerm::create(
      model, SNLTerm::Direction::Input, 1, 0, NLName("D"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* output = SNLBusTerm::create(
      model, SNLTerm::Direction::Output, 1, 0, NLName("Q"));
  SNLDesignModeling::addInputsToClockArcs(collectBitTerms(data), clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, collectBitTerms(output));
  return model;
}

SNLDesign* createNoDataSequentialModel(NLLibrary* library,
                                       const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  SNLDesignModeling::addClockToOutputsArcs(clock, {output});
  return model;
}

SNLDesign* createExtraUpdatePinSequentialModel(NLLibrary* library,
                                               const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* address =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  SNLDesignModeling::addInputsToClockArcs({data, address}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {output});
  return model;
}

SNLDesign* createResetSetSequentialModel(NLLibrary* library,
                                         const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* reset =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("R"));
  auto* set =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("S"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  SNLDesignModeling::addInputsToClockArcs({data, reset, set}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {output});
  return model;
}

SNLDesign* createNamedComplementSetSequentialModel(
    NLLibrary* library,
    const std::string& name,
    const std::string& primaryPinName,
    const std::string& complementPinName) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* data =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  auto* set =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("S"));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CK"));
  auto* primary = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(primaryPinName));
  auto* complement = SNLScalarTerm::create(
      model, SNLTerm::Direction::Output, NLName(complementPinName));
  SNLDesignModeling::addInputsToClockArcs({data, set}, clock);
  SNLDesignModeling::addClockToOutputsArcs(clock, {primary, complement});
  return model;
}

SNLDesign* createSequentialOutputPairTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel,
    const std::string& primaryPinName,
    const std::string& secondaryPinName) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topPrimary =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out_primary"));
  auto* topSecondary = SNLScalarTerm::create(
      top, SNLTerm::Direction::Output, NLName("out_secondary"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netPrimary = SNLScalarNet::create(top, NLName("net_primary"));
  auto* netSecondary = SNLScalarNet::create(top, NLName("net_secondary"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topPrimary->setNet(netPrimary);
  topSecondary->setNet(netSecondary);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("D")))->setNet(netIn);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName(primaryPinName)))->setNet(
      netPrimary);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName(secondaryPinName)))->setNet(
      netSecondary);

  return top;
}

SNLDesign* createSetOnlySequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topSet =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("set"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netSet = SNLScalarNet::create(top, NLName("net_set"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topSet->setNet(netSet);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("D")))->setNet(netIn);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("S")))->setNet(netSet);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("Q")))->setNet(netOut);

  return top;
}

SNLDesign* createBusSequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn = SNLBusTerm::create(
      top, SNLTerm::Direction::Input, 1, 0, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut = SNLBusTerm::create(
      top, SNLTerm::Direction::Output, 1, 0, NLName("out"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLBusNet::create(top, 1, 0, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLBusNet::create(top, 1, 0, NLName("net_out"));

  topClock->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);

  auto* modelData = sequentialModel->getBusTerm(NLName("D"));
  auto* modelOutput = sequentialModel->getBusTerm(NLName("Q"));
  for (int bit = 0; bit <= 1; ++bit) {
    topIn->getBit(bit)->setNet(netIn->getBit(bit));
    topOut->getBit(bit)->setNet(netOut->getBit(bit));
    seq->getInstTerm(modelData->getBit(bit))->setNet(netIn->getBit(bit));
    seq->getInstTerm(modelOutput->getBit(bit))->setNet(netOut->getBit(bit));
  }

  return top;
}

SNLDesign* createComplementedSetSequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel,
    const std::string& primaryPinName,
    const std::string& complementPinName) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topSet =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("set"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topPrimary =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out_primary"));
  auto* topSecondary = SNLScalarTerm::create(
      top, SNLTerm::Direction::Output, NLName("out_secondary"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netSet = SNLScalarNet::create(top, NLName("net_set"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netPrimary = SNLScalarNet::create(top, NLName("net_primary"));
  auto* netSecondary = SNLScalarNet::create(top, NLName("net_secondary"));

  topIn->setNet(netIn);
  topSet->setNet(netSet);
  topClock->setNet(netClock);
  topPrimary->setNet(netPrimary);
  topSecondary->setNet(netSecondary);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("D")))->setNet(netIn);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("S")))->setNet(netSet);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName(primaryPinName)))->setNet(
      netPrimary);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName(complementPinName)))->setNet(
      netSecondary);

  return top;
}

SNLDesign* createNoDataSequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topClock->setNet(netClock);
  topOut->setNet(netOut);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("Q")))->setNet(netOut);

  return top;
}

SNLDesign* createExtraUpdatePinSequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topAddr =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("addr"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netAddr = SNLScalarNet::create(top, NLName("net_addr"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topAddr->setNet(netAddr);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("D")))->setNet(netIn);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("A")))->setNet(netAddr);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("Q")))->setNet(netOut);

  return top;
}

SNLDesign* createPartialCoverageNoDriverTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topGood =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("good"));
  auto* topBad =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("bad"));

  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netData = SNLScalarNet::create(top, NLName("net_data"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topGood->setNet(netIn);
  topBad->setNet(netQ);

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netData);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createPartialCoverageMultiDriverTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topInA =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in_a"));
  auto* topInB =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in_b"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topGood =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("good"));
  auto* topBad =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("bad"));

  auto* inv0 = SNLInstance::create(top, invModel, NLName("inv0"));
  auto* inv1 = SNLInstance::create(top, invModel, NLName("inv1"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  auto* netInA = SNLScalarNet::create(top, NLName("net_in_a"));
  auto* netInB = SNLScalarNet::create(top, NLName("net_in_b"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netMulti = SNLScalarNet::create(top, NLName("net_multi"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topInA->setNet(netInA);
  topInB->setNet(netInB);
  topClock->setNet(netClock);
  topGood->setNet(netInA);
  topBad->setNet(netQ);

  inv0->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netInA);
  inv0->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netMulti);
  inv1->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netInB);
  inv1->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netMulti);

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netMulti);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createUnsupportedPrimitiveCoverageTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topGood =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("good"));
  auto* topBad =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("bad"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topGood->setNet(netIn);
  topBad->setNet(netOut);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("Q")))->setNet(netOut);

  return top;
}

SNLDesign* createCombinationalInvTop(NLLibrary* library,
                                     const std::string& name,
                                     SNLDesign* invModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* inv = SNLInstance::create(top, invModel, NLName("inv0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topOut->setNet(netOut);
  inv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netIn);
  inv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netOut);

  return top;
}

SNLDesign* createHierarchicalCombinationalInvChild(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel) {
  auto* child =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* childIn =
      SNLScalarTerm::create(child, SNLTerm::Direction::Input, NLName("in"));
  auto* childOut =
      SNLScalarTerm::create(child, SNLTerm::Direction::Output, NLName("out"));

  auto* inv = SNLInstance::create(child, invModel, NLName("inv0"));
  auto* netIn = SNLScalarNet::create(child, NLName("net_in"));
  auto* netOut = SNLScalarNet::create(child, NLName("net_out"));

  childIn->setNet(netIn);
  childOut->setNet(netOut);
  inv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netIn);
  inv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netOut);

  return child;
}

SNLDesign* createHierarchicalCombinationalInvTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* childModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* child = SNLInstance::create(top, childModel, NLName("child0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topOut->setNet(netOut);
  child->getInstTerm(childModel->getScalarTerm(NLName("in")))->setNet(netIn);
  child->getInstTerm(childModel->getScalarTerm(NLName("out")))->setNet(netOut);

  return top;
}

SNLDesign* createResetSetSequentialTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* sequentialModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topSet =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("set"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* seq = SNLInstance::create(top, sequentialModel, NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netSet = SNLScalarNet::create(top, NLName("net_set"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topReset->setNet(netReset);
  topSet->setNet(netSet);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("D")))->setNet(netIn);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("R")))->setNet(netReset);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("S")))->setNet(netSet);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  seq->getInstTerm(sequentialModel->getScalarTerm(NLName("Q")))->setNet(netOut);

  return top;
}

SNLDesign* createDffreTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topEnable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("en"));
  auto* topReset =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* ff = SNLInstance::create(top, NLDB0::getDFFRE(), NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netEnable = SNLScalarNet::create(top, NLName("net_en"));
  auto* netReset = SNLScalarNet::create(top, NLName("net_rst"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topEnable->setNet(netEnable);
  topReset->setNet(netReset);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  ff->getInstTerm(NLDB0::getDFFREData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFREEnable())->setNet(netEnable);
  ff->getInstTerm(NLDB0::getDFFREReset())->setNet(netReset);
  ff->getInstTerm(NLDB0::getDFFREClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFREOutput())->setNet(netOut);

  return top;
}

SNLDesign* createComplementedOutputTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* ffModel,
    SNLDesign* invModel,
    bool rebuildOutputsFromComplements) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOutQ =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out_q"));
  auto* topOutQn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out_qn"));

  auto* ff = SNLInstance::create(top, ffModel, NLName("ff0"));
  SNLInstance* qnToQInv = nullptr;
  SNLInstance* qToQnInv = nullptr;
  if (rebuildOutputsFromComplements) {
    qnToQInv = SNLInstance::create(top, invModel, NLName("inv_qn_to_q"));
    qToQnInv = SNLInstance::create(top, invModel, NLName("inv_q_to_qn"));
  }

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));
  auto* netQn = SNLScalarNet::create(top, NLName("net_qn"));
  auto* netOutQ = SNLScalarNet::create(top, NLName("net_out_q"));
  auto* netOutQn = SNLScalarNet::create(top, NLName("net_out_qn"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);

  ff->getInstTerm(ffModel->getScalarTerm(NLName("CK")))->setNet(netClock);
  ff->getInstTerm(ffModel->getScalarTerm(NLName("D")))->setNet(netIn);
  ff->getInstTerm(ffModel->getScalarTerm(NLName("Q")))->setNet(netQ);
  ff->getInstTerm(ffModel->getScalarTerm(NLName("QN")))->setNet(netQn);

  if (rebuildOutputsFromComplements) {
    topOutQ->setNet(netOutQ);
    topOutQn->setNet(netOutQn);
    qnToQInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netQn);
    qnToQInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netOutQ);
    qToQnInv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netQ);
    qToQnInv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netOutQn);
  } else {
    topOutQ->setNet(netQ);
    topOutQn->setNet(netQn);
  }

  return top;
}

SignalKey findKeyByDisplayName(const SequentialDesignModel& model,
                               const std::string& displayName) {
  for (const auto& [key, currentName] : model.displayNameByKey) {
    if (currentName == displayName) {
      return key;
    }
  }
  throw std::runtime_error("Missing display name in extracted model: " + displayName);
}

}  // namespace

TEST_F(SequentialEquivalenceStrategyTests, IdenticalDffDesignsAreEquivalent) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, IdenticalDffDesignsAreEquivalentWithPdrEngine) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library = NLLibrary::create(db, NLName("LIB"));
  auto* primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
  auto* invModel = createInvModel(primitives);

  auto* top0 =
      createDffTop(library, "top0", invModel, false, false, "in", "out", "ff0");
  auto* top1 =
      createDffTop(library, "top1", invModel, false, false, "in", "out", "ff1");

  SequentialEquivalenceStrategy strategy(
      top0,
      top1,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IdenticalDffDesignsAreEquivalentWithKInductionEngine) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library = NLLibrary::create(db, NLName("LIB"));
  auto* primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
  auto* invModel = createInvModel(primitives);

  auto* top0 =
      createDffTop(library, "top0", invModel, false, false, "in", "out", "ff0");
  auto* top1 =
      createDffTop(library, "top1", invModel, false, false, "in", "out", "ff1");

  SequentialEquivalenceStrategy strategy(
      top0,
      top1,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, IdenticalDffDesignsAreEquivalentWithImcEngine) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library = NLLibrary::create(db, NLName("LIB"));
  auto* primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
  auto* invModel = createInvModel(primitives);

  auto* top0 =
      createDffTop(library, "top0", invModel, false, false, "in", "out", "ff0");
  auto* top1 =
      createDffTop(library, "top1", invModel, false, false, "in", "out", "ff1");

  SequentialEquivalenceStrategy strategy(
      top0,
      top1,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, OutputMismatchFailsAfterInitialObservation) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, true);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests, NextStateMismatchFailsAtOneStep) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, true, false);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, DffeHoldSemanticsAreProved) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createDffeTop(library, "top0");
  auto* top1 = createDffeTop(library, "top1");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, ComplementedStateOutputsRemainConsistent) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* dffQnModel = createDffQnModel(primitives);
  auto* top0 =
      createComplementedOutputTop(library, "top0", dffQnModel, invModel, false);
  auto* top1 =
      createComplementedOutputTop(library, "top1", dffQnModel, invModel, true);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests, EquivalentDesignsWithRenamedStateAreAccepted) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false, "state_a");
  auto* top1 = createDffTop(library, "top1", invModel, false, false, "state_b");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RenamedStatePipelineIsProvedWithoutNameBasedStateMatching) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(
      library, "top0", false, {"left0", "left1", "left2"});
  auto* top1 = createResetInitializedPipelineTop(
      library, "top1", false, {"right0", "right1", "right2"});

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetInitializedThreeStagePipelineFailsAtThreeSteps) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(library, "top0", false);
  auto* top1 = createResetInitializedPipelineTop(library, "top1", true);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(4);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetInitializedEquivalentPipelineIsProved) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(library, "top0", false);
  auto* top1 = createResetInitializedPipelineTop(library, "top1", false);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetInitializedRenamedPipelineClosesWithinThreeStepSecProof) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(
      library, "top0", false, {"ff0", "ff1", "ff2"});
  auto* top1 = createResetInitializedPipelineTop(
      library, "top1", false, {"state_a", "state_b", "state_c"});

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(4);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapEquivalentPipelineIsProved) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* top0 =
      createBootstrapPipelineTop(library, "top0", invModel, andModel);
  auto* top1 =
      createBootstrapPipelineTop(library, "top1", invModel, andModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapCanAnchorEqualStatesWithoutConstantValues) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* orModel = createOr2Model(primitives);
  auto* top0 =
      createResetLoadsInputTop(library, "top0", invModel, andModel, orModel);
  auto* top1 =
      createResetLoadsInputTop(library, "top1", invModel, andModel, orModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapCanAnchorHiddenEqualStatesWithoutConstantValues) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* orModel = createOr2Model(primitives);
  auto* top0 = createResetLoadsInputTwoStageTop(
      library, "top0", invModel, andModel, orModel);
  auto* top1 = createResetLoadsInputTwoStageTop(
      library, "top1", invModel, andModel, orModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapAutomaticallyExtendsForHiddenShiftPipelines) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* orModel = createOr2Model(primitives);
  auto* top0 = createResetLoadsInputShiftPipelineTopWithStages(
      library, "top0", invModel, andModel, orModel, 20);
  auto* top1 = createResetLoadsInputShiftPipelineTopWithStages(
      library, "top1", invModel, andModel, orModel, 20);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapLongEquivalentPipelineStillClosesAtSmallK) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* top0 =
      createBootstrapPipelineTopWithStages(library, "top0", invModel, andModel, 12);
  auto* top1 =
      createBootstrapPipelineTopWithStages(library, "top1", invModel, andModel, 12);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       StructuralInvariantHandlesMismatchedStateCountsWithoutOscillation) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedShiftPipelineTopWithStages(
      library, "top0", 5);
  auto* top1 = createResetInitializedShiftPipelineTopWithStages(
      library, "top1", 1);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(6);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_LE(result.bound, 6u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantUsesExplicitInitialCompatibilityWithoutReset) {
  const SignalKey state0A = makeSignalKey("state0A");
  const SignalKey state0B = makeSignalKey("state0B");
  const SignalKey state1A = makeSignalKey("state1A");
  const SignalKey state1B = makeSignalKey("state1B");

  SequentialDesignModel model0;
  model0.stateBits = {state0A, state0B};
  model0.initialStateValueByKey.emplace(state0A, false);
  model0.initialStateValueByKey.emplace(state0B, true);

  SequentialDesignModel model1;
  model1.stateBits = {state1A, state1B};
  model1.initialStateValueByKey.emplace(state1A, false);
  model1.initialStateValueByKey.emplace(state1B, false);

  AlignedSignals candidateStates;
  candidateStates.names = {"state_a", "state_b"};
  candidateStates.keys0 = {state0A, state0B};
  candidateStates.keys1 = {state1A, state1B};

  const auto invariant =
      buildReachableStateInvariant(model0, model1, AlignedSignals{}, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 0u);
  ASSERT_EQ(invariant.initialStateCorrespondence.names.size(), 1u);
  EXPECT_EQ(invariant.initialStateCorrespondence.names[0], "state_a");
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 1u);
  EXPECT_EQ(invariant.anchoredStateEqualities.names[0], "state_a");
  EXPECT_TRUE(invariant.bootstrapValues0.empty());
  EXPECT_TRUE(invariant.bootstrapValues1.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantSkipsBootstrapWhenResetAndInitialStateAreComplete) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.initialStateValueByKey.emplace(state0, false);

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.initialStateValueByKey.emplace(state1, false);

  AlignedSignals alignedInputs;
  alignedInputs.names = {"rst"};
  alignedInputs.keys0 = {rst0};
  alignedInputs.keys1 = {rst1};

  AlignedSignals candidateStates;
  candidateStates.names = {"state"};
  candidateStates.keys0 = {state0};
  candidateStates.keys1 = {state1};

  const auto invariant =
      buildReachableStateInvariant(model0, model1, alignedInputs, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 0u);
  ASSERT_EQ(invariant.initialStateCorrespondence.names.size(), 1u);
  EXPECT_EQ(invariant.initialStateCorrespondence.names[0], "state");
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 1u);
  EXPECT_EQ(invariant.anchoredStateEqualities.names[0], "state");
  EXPECT_TRUE(invariant.bootstrapValues0.empty());
  EXPECT_TRUE(invariant.bootstrapValues1.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantDerivesBootstrapValuesAndAnchorsFromReset) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));

  AlignedSignals alignedInputs;
  alignedInputs.names = {"rst"};
  alignedInputs.keys0 = {rst0};
  alignedInputs.keys1 = {rst1};

  AlignedSignals candidateStates;
  candidateStates.names = {"state"};
  candidateStates.keys0 = {state0};
  candidateStates.keys1 = {state1};

  const auto invariant =
      buildReachableStateInvariant(model0, model1, alignedInputs, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.initialStateCorrespondence.names.size(), 1u);
  EXPECT_EQ(invariant.initialStateCorrespondence.names[0], "state");
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 1u);
  EXPECT_EQ(invariant.anchoredStateEqualities.names[0], "state");
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantBootstrapRecoversEqualitiesAfterMismatchedInitialValues) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");
  const SignalKey shadow0 = makeSignalKey("shadow0");
  const SignalKey shadow1 = makeSignalKey("shadow1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {state0, shadow0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(state0, 3);
  model0.inputVarByKey.emplace(shadow0, 6);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.initialStateValueByKey.emplace(state0, false);
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));
  model0.nextStateExprByStateKey.emplace(shadow0, BoolExpr::Var(6));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {state1, shadow1};
  model1.inputVarByKey.emplace(rst1, 4);
  model1.inputVarByKey.emplace(state1, 5);
  model1.inputVarByKey.emplace(shadow1, 7);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.initialStateValueByKey.emplace(state1, true);
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(4)));
  model1.nextStateExprByStateKey.emplace(shadow1, BoolExpr::Var(7));

  AlignedSignals alignedInputs;
  alignedInputs.names = {"rst"};
  alignedInputs.keys0 = {rst0};
  alignedInputs.keys1 = {rst1};

  AlignedSignals candidateStates;
  candidateStates.names = {"state"};
  candidateStates.keys0 = {state0};
  candidateStates.keys1 = {state1};

  ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto invariant = buildReachableStateInvariant(
      model0, model1, alignedInputs, candidateStates, true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  EXPECT_TRUE(invariant.initialStateCorrespondence.names.empty());
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 1u);
  EXPECT_EQ(invariant.anchoredStateEqualities.names[0], "state");
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
  EXPECT_NE(stderrOutput.find("SEC diag: bootstrap step 1"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantFallsBackWhenBootstrapHasNoCandidates) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(state0, 4);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(state1, 5);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));

  const auto invariant = buildReachableStateInvariant(
      model0, model1, AlignedSignals{}, AlignedSignals{});

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  EXPECT_TRUE(invariant.anchoredStateEqualities.names.empty());
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantFallsBackWhenOnlyOneSideHasResetAssignments) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(state0, 3);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.initialStateValueByKey.emplace(state0, false);

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(state1, 4);
  model1.initialStateValueByKey.emplace(state1, false);

  AlignedSignals candidateStates;
  candidateStates.names = {"state"};
  candidateStates.keys0 = {state0};
  candidateStates.keys1 = {state1};

  const auto invariant = buildReachableStateInvariant(
      model0, model1, AlignedSignals{}, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 0u);
  ASSERT_EQ(invariant.initialStateCorrespondence.names.size(), 1u);
  EXPECT_EQ(invariant.initialStateCorrespondence.names[0], "state");
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 1u);
  EXPECT_EQ(invariant.anchoredStateEqualities.names[0], "state");
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantCoversBootstrapValuePropagationEdgeCases) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey truth0 = makeSignalKey("truth0");
  const SignalKey truth1 = makeSignalKey("truth1");
  const SignalKey invalid0 = makeSignalKey("invalid0");
  const SignalKey invalid1 = makeSignalKey("invalid1");
  const SignalKey null0 = makeSignalKey("null0");
  const SignalKey null1 = makeSignalKey("null1");
  const SignalKey diff0 = makeSignalKey("diff0");
  const SignalKey diff1 = makeSignalKey("diff1");
  const SignalKey hidden0 = makeSignalKey("hidden0");
  const SignalKey hidden1 = makeSignalKey("hidden1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {truth0, invalid0, null0, diff0, hidden0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(truth0, 3);
  model0.inputVarByKey.emplace(invalid0, 4);
  model0.inputVarByKey.emplace(null0, 5);
  model0.inputVarByKey.emplace(diff0, 6);
  model0.inputVarByKey.emplace(hidden0, 7);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.nextStateExprByStateKey.emplace(
      truth0,
      BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::createFalse()));
  model0.nextStateExprByStateKey.emplace(invalid0, BoolExpr::Var(7));
  model0.nextStateExprByStateKey.emplace(null0, nullptr);
  model0.nextStateExprByStateKey.emplace(diff0, BoolExpr::Var(2));
  model0.nextStateExprByStateKey.emplace(hidden0, BoolExpr::Var(7));
  model0.initialStateValueByKey.emplace(truth0, false);
  model0.initialStateValueByKey.emplace(invalid0, false);
  model0.initialStateValueByKey.emplace(null0, false);
  model0.initialStateValueByKey.emplace(diff0, false);

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {truth1, invalid1, null1, diff1, hidden1};
  model1.inputVarByKey.emplace(rst1, 10);
  model1.inputVarByKey.emplace(truth1, 11);
  model1.inputVarByKey.emplace(invalid1, 12);
  model1.inputVarByKey.emplace(null1, 13);
  model1.inputVarByKey.emplace(diff1, 14);
  model1.inputVarByKey.emplace(hidden1, 15);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(
      truth1,
      BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::createFalse()));
  model1.nextStateExprByStateKey.emplace(invalid1, BoolExpr::Var(15));
  model1.nextStateExprByStateKey.emplace(null1, nullptr);
  model1.nextStateExprByStateKey.emplace(diff1, BoolExpr::Not(BoolExpr::Var(10)));
  model1.nextStateExprByStateKey.emplace(hidden1, BoolExpr::Var(15));
  model1.initialStateValueByKey.emplace(truth1, true);
  model1.initialStateValueByKey.emplace(invalid1, true);
  model1.initialStateValueByKey.emplace(null1, true);
  model1.initialStateValueByKey.emplace(diff1, true);

  AlignedSignals alignedInputs;
  alignedInputs.names = {"rst"};
  alignedInputs.keys0 = {rst0};
  alignedInputs.keys1 = {rst1};

  AlignedSignals candidateStates;
  candidateStates.names = {"truth", "invalid", "null", "diff"};
  candidateStates.keys0 = {truth0, invalid0, null0, diff0};
  candidateStates.keys1 = {truth1, invalid1, null1, diff1};

  const auto invariant = buildReachableStateInvariant(
      model0, model1, alignedInputs, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  EXPECT_TRUE(invariant.initialStateCorrespondence.names.empty());
  EXPECT_TRUE(invariant.bootstrapValues0.at(truth0));
  EXPECT_TRUE(invariant.bootstrapValues1.at(truth1));
  EXPECT_TRUE(invariant.bootstrapValues0.at(diff0));
  EXPECT_FALSE(invariant.bootstrapValues1.at(diff1));
  EXPECT_TRUE(
      std::find(
          invariant.anchoredStateEqualities.names.begin(),
          invariant.anchoredStateEqualities.names.end(),
          "truth") != invariant.anchoredStateEqualities.names.end());
  EXPECT_TRUE(
      std::find(
          invariant.anchoredStateEqualities.names.begin(),
          invariant.anchoredStateEqualities.names.end(),
          "diff") == invariant.anchoredStateEqualities.names.end());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantBootstrapValuesEvaluateConstTrueXorAndInvalidExprStates) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey const0 = makeSignalKey("const0");
  const SignalKey const1 = makeSignalKey("const1");
  const SignalKey xor0 = makeSignalKey("xor0");
  const SignalKey xor1 = makeSignalKey("xor1");
  const SignalKey invalid0 = makeSignalKey("invalid0");
  const SignalKey invalid1 = makeSignalKey("invalid1");
  const SignalKey hidden0 = makeSignalKey("hidden0");
  const SignalKey hidden1 = makeSignalKey("hidden1");
  BoolExpr invalidExpr0;
  BoolExpr invalidExpr1;

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {const0, xor0, invalid0, hidden0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(const0, 3);
  model0.inputVarByKey.emplace(xor0, 4);
  model0.inputVarByKey.emplace(invalid0, 5);
  model0.inputVarByKey.emplace(hidden0, 6);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.initialStateValueByKey.emplace(const0, false);
  model0.initialStateValueByKey.emplace(xor0, false);
  model0.initialStateValueByKey.emplace(invalid0, false);
  model0.nextStateExprByStateKey.emplace(const0, BoolExpr::createTrue());
  model0.nextStateExprByStateKey.emplace(
      xor0, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::createFalse()));
  model0.nextStateExprByStateKey.emplace(invalid0, &invalidExpr0);
  model0.nextStateExprByStateKey.emplace(hidden0, BoolExpr::Var(6));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {const1, xor1, invalid1, hidden1};
  model1.inputVarByKey.emplace(rst1, 10);
  model1.inputVarByKey.emplace(const1, 11);
  model1.inputVarByKey.emplace(xor1, 12);
  model1.inputVarByKey.emplace(invalid1, 13);
  model1.inputVarByKey.emplace(hidden1, 14);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.initialStateValueByKey.emplace(const1, false);
  model1.initialStateValueByKey.emplace(xor1, false);
  model1.initialStateValueByKey.emplace(invalid1, false);
  model1.nextStateExprByStateKey.emplace(const1, BoolExpr::createTrue());
  model1.nextStateExprByStateKey.emplace(
      xor1, BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::createFalse()));
  model1.nextStateExprByStateKey.emplace(invalid1, &invalidExpr1);
  model1.nextStateExprByStateKey.emplace(hidden1, BoolExpr::Var(14));

  AlignedSignals candidateStates;
  candidateStates.names = {"const", "xor", "invalid"};
  candidateStates.keys0 = {const0, xor0, invalid0};
  candidateStates.keys1 = {const1, xor1, invalid1};

  const auto invariant = buildReachableStateInvariant(
      model0, model1, AlignedSignals{}, candidateStates);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.anchoredStateEqualities.names.size(), 3u);
  EXPECT_TRUE(invariant.bootstrapValues0.at(const0));
  EXPECT_TRUE(invariant.bootstrapValues1.at(const1));
  EXPECT_TRUE(invariant.bootstrapValues0.at(xor0));
  EXPECT_TRUE(invariant.bootstrapValues1.at(xor1));
  EXPECT_EQ(invariant.bootstrapValues0.count(invalid0), 0u);
  EXPECT_EQ(invariant.bootstrapValues1.count(invalid1), 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprRemapThrowsOnMissingVariableMapping) {
  auto* expr = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));

  EXPECT_THROW(
      static_cast<void>(remapBoolExprVariables(expr, {{2, 10}})),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprHelpersCoverNullXorAndInvalidOperators) {
  EXPECT_EQ(remapBoolExprVariables(nullptr, {}), nullptr);
  EXPECT_EQ(substituteBoolExprVariables(nullptr, {}), nullptr);

  BoolExpr invalid;
  EXPECT_THROW(
      static_cast<void>(remapBoolExprVariables(&invalid, {})),
      std::runtime_error);

  auto* xorExpr = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  auto* substituted =
      substituteBoolExprVariables(xorExpr, {{2, true}, {3, false}});
  EXPECT_TRUE(substituted->evaluate({}));

  EXPECT_THROW(
      static_cast<void>(substituteBoolExprVariables(&invalid, {})),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprSubstitutionRewritesAssignedVariablesAndKeepsOthers) {
  auto* expr = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Not(BoolExpr::Var(3)));
  auto* substituted = substituteBoolExprVariables(expr, {{2, true}, {3, false}});

  EXPECT_TRUE(substituted->evaluate({}));
  EXPECT_EQ(substituted->getOp(), KEPLER_FORMAL::Op::VAR);
  EXPECT_EQ(substituted->getId(), 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SatEncodingHelpersCoverConstantCachingAndErrorBranches) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  FrameVariableStore variables(solver, {2, 3}, 2);

  EXPECT_THROW(
      static_cast<void>(variables.getLiteral(99, 0)),
      std::runtime_error);
  EXPECT_THROW(
      static_cast<void>(variables.makeLeafLits(3)),
      std::runtime_error);

  FrameFormulaEncoder encoder(solver, variables.makeLeafLits(0));
  EXPECT_THROW(static_cast<void>(encoder.encode(nullptr)), std::invalid_argument);
  EXPECT_THROW(
      static_cast<void>(encoder.encode(BoolExpr::Var(99))),
      std::runtime_error);

  BoolExpr invalid;
  EXPECT_THROW(static_cast<void>(encoder.encode(&invalid)), std::runtime_error);

  const int trueLit = encoder.encode(BoolExpr::createTrue());
  EXPECT_EQ(trueLit, encoder.encode(BoolExpr::createTrue()));

  addSimplePathConstraint(solver, variables, {}, 2);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SatSolverWrapperGetLiteralValueHandlesConstantsUnknownModelsAndErrors) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::GLUCOSE);
  const int symbol = solver.newVar() + 2;
  EXPECT_TRUE(solver.solve());

  EXPECT_FALSE(solver.getLiteralValue(0));
  EXPECT_TRUE(solver.getLiteralValue(1));
  EXPECT_FALSE(solver.getLiteralValue(symbol));
  EXPECT_THROW(static_cast<void>(solver.getLiteralValue(-1)), std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverFindsCombinationalCounterexampleAtFrameZero) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.allSymbols = {2};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Not(BoolExpr::Var(2))};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 0u);
  ASSERT_EQ(witness->inputTrace.size(), 1u);
  EXPECT_EQ(witness->inputTrace[0].frame, 0u);
  ASSERT_EQ(witness->outputMismatches.size(), 1u);
  EXPECT_EQ(witness->outputMismatches[0].signal, "out");
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverObservationOnlyStartsSearchingAtFrameOne) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(2)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
  ASSERT_EQ(witness->inputTrace.size(), 2u);
  EXPECT_EQ(witness->inputTrace.front().frame, 0u);
  EXPECT_EQ(witness->inputTrace.back().frame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverOffsetsWitnessAfterResetBootstrap) {
  KInductionProblem problem;
  problem.environmentInputNames = {"rst", "in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2, 5};
  problem.resetBootstrapCycles = 2;
  problem.resetBootstrapInputs = {{2, true}};
  problem.bootstrapStateAssignments = {{3, false}, {4, false}};
  problem.bootstrapStateEqualityPairs = {{3, 4}};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4, 5};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(5)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
  ASSERT_EQ(witness->inputTrace.size(), 2u);
  EXPECT_EQ(witness->inputTrace[0].frame, 0u);
  ASSERT_EQ(witness->inputTrace[0].assignments.size(), 2u);
  EXPECT_EQ(witness->inputTrace[0].assignments[0].signal, "rst");
  EXPECT_FALSE(witness->inputTrace[0].assignments[0].value);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverHandlesActiveLowResetBootstrapInputs) {
  KInductionProblem problem;
  problem.environmentInputNames = {"rst_n", "in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2, 5};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{2, false}};
  problem.bootstrapStateAssignments = {{3, false}, {4, false}};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4, 5};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(5)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
  ASSERT_EQ(witness->inputTrace.size(), 2u);
  // The reported witness trace is offset past the hidden bootstrap frame,
  // so an active-low reset is already deasserted in the visible input trace.
  EXPECT_TRUE(witness->inputTrace[0].assignments[0].value);
  EXPECT_TRUE(witness->inputTrace[1].assignments[0].value);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverPartialInitWithoutStateRelationUsesObservationFallback) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(2)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(3));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 2;
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverPartialInitWithStateRelationKeepsFrameZeroAligned) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(2)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(3));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 2;
  problem.initialStateEqualityPairs = {{3, 4}};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ExactInterpolantSynthesizerDerivesOneStepReachableStateInvariant) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2, 3};
  problem.inputSymbols = {3};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ExactInterpolantSynthesizer engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto interpolant = engine.deriveOneStepReachableStateInvariant(4);

  ASSERT_TRUE(interpolant.has_value());
  EXPECT_TRUE((*interpolant)->evaluate({{2, false}}));
  EXPECT_FALSE((*interpolant)->evaluate({{2, true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ExactInterpolantSynthesizerReturnsNulloptWhenStateBudgetIsExceeded) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.transitions0.emplace_back(2, BoolExpr::Var(2));
  problem.transitions0.emplace_back(3, BoolExpr::Var(3));
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ExactInterpolantSynthesizer engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto interpolant = engine.deriveOneStepReachableStateInvariant(1);

  EXPECT_FALSE(interpolant.has_value());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ExactInterpolantSynthesizerReturnsNulloptWhenBadIsReachableInOneStep) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ExactInterpolantSynthesizer engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto interpolant = engine.deriveOneStepReachableStateInvariant(4);

  EXPECT_FALSE(interpolant.has_value());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ExactInterpolantSynthesizerRejectsNonInductiveInterpolant) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.transitions0.emplace_back(2, BoolExpr::Var(3));
  problem.transitions0.emplace_back(3, BoolExpr::createTrue());
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ExactInterpolantSynthesizer engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto interpolant = engine.deriveOneStepReachableStateInvariant(4);

  EXPECT_FALSE(interpolant.has_value());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ExactInterpolantSynthesizerUsesBootstrapAssignmentsAndComplementedStatePairs) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {{2, false}, {3, true}};
  problem.bootstrapStateEqualityPairs = {{2, 2}};
  problem.complementedStatePairs0 = {{2, 3}};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.transitions0.emplace_back(3, BoolExpr::createTrue());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ExactInterpolantSynthesizer engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto interpolant = engine.deriveOneStepReachableStateInvariant(4);

  ASSERT_TRUE(interpolant.has_value());
  EXPECT_TRUE((*interpolant)->evaluate({{2, false}, {3, true}}));
  EXPECT_FALSE((*interpolant)->evaluate({{2, true}, {3, false}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineProvesEquivalentSmallTransitionSystem) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFindsReachableBadState) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDoesNotConvergeAtFourFramesWithoutInitialConstraint) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(4);

  EXPECT_EQ(result.status, PDRStatus::Inconclusive);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineStillAcceptsImmediateProofWhenFrameBudgetIsZero) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesBootstrapAssignmentsAndComplementedStatePairs) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {{2, false}, {3, true}};
  problem.bootstrapStateEqualityPairs = {{2, 2}};
  problem.complementedStatePairs0 = {{2, 3}};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.transitions0.emplace_back(3, BoolExpr::createTrue());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDoesNotReuseNonInductiveStrengtheningAsFrameInvariant) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::Var(2));
  problem.bad = BoolExpr::Var(3);
  problem.property = BoolExpr::Not(problem.bad);
  // Init implies !x, but it is not inductive: x becomes true in the next step.
  // Reusing it as a frame fact would incorrectly hide the real bad state 11.
  problem.inductionProperty = BoolExpr::Not(BoolExpr::Var(2));
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesPropertyAsFallbackImmediateProof) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  // This strengthening is not implied by Init, so PDR must fall back to the
  // checked SEC property instead of dropping straight into the full clause loop.
  problem.inductionProperty = BoolExpr::Var(2);
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineProvesEquivalentExactlyAtThreeFrames) {
  const auto problem = buildLinearChainSecProblem(4);

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineProvesEquivalentSmallTransitionSystem) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineFindsReachableBadState) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, KInductionStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.bound, 1u);
  EXPECT_EQ(result.witness->badFrame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineProvesEquivalentWithExactInterpolant) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.transitions0.emplace_back(2, BoolExpr::Var(3));
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineUsesValidatedSharedStrengtheningInvariant) {
  KInductionProblem problem;
  problem.observedOutputNames = {"out"};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{2, BoolExpr::Var(3)}, {3, BoolExpr::Var(3)}};
  problem.transitions1 = {{4, BoolExpr::Var(5)}, {5, BoolExpr::Var(5)}};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)), BoolExpr::Not(BoolExpr::Var(5))));
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.property =
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4)));
  problem.bad = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4));
  problem.inductionProperty = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4))),
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(5))));
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineFindsReachableBadState) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, IMCStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.bound, 1u);
  EXPECT_EQ(result.witness->badFrame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineProvesEquivalentExactlyAtThreeFrames) {
  const auto problem = buildLinearChainSecProblem(4);

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  EXPECT_EQ(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineRemainsInconclusiveAtFourFramesWhenFiveAreNeeded) {
  const auto problem = buildLinearChainSecProblem(6);

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(4);

  EXPECT_EQ(result.status, IMCStatus::Inconclusive);
  EXPECT_EQ(result.bound, 4u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineProvesEquivalentExactlyAtThreeFrames) {
  const auto problem = buildLinearChainSecProblem(4);

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineRemainsInconclusiveAtFourFramesWhenFiveAreNeeded) {
  const auto problem = buildLinearChainSecProblem(6);

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(4);

  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(result.bound, 4u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       InductionStepSolverUsesExplicitInvariantWhenProvided) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.allSymbols = {2, 3};
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  problem.transitions1 = {{3, BoolExpr::Var(3)}};
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.inductionProperty =
      makeEqualityExpr(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  EXPECT_TRUE(
      provesByInduction(problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineCombinationalProblemReturnsImmediately) {
  KInductionProblem problem;
  problem.allSymbols = {2};
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SynthesizedResetInferencePropagatesThroughLongBootstrapPipeline) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* top = createBootstrapPipelineTopWithStages(
      library, "top", invModel, andModel, 12);

  const auto model = SequentialDesignModel::extract(top);

  EXPECT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_EQ(model.initialStateValueByKey.size(), model.stateBits.size());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SynthesizedResetInferenceScalesPastLargeStateCutoff) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* andModel = createAnd2Model(primitives);
  auto* top = createBootstrapPipelineTopWithStages(
      library, "top", invModel, andModel, 2200);

  const auto model = SequentialDesignModel::extract(top);

  EXPECT_FALSE(model.hasUnsupportedFeatures());
  EXPECT_EQ(model.initialStateValueByKey.size(), model.stateBits.size());
}

TEST_F(SequentialEquivalenceStrategyTests,
       StructuralStateInvariantRefinesConstantsUnknownVarsAndInvalidNodes) {
  const SignalKey a0 = makeSignalKey("a0");
  const SignalKey b0 = makeSignalKey("b0");
  const SignalKey c0 = makeSignalKey("c0");
  const SignalKey a1 = makeSignalKey("a1");
  const SignalKey b1 = makeSignalKey("b1");
  const SignalKey c1 = makeSignalKey("c1");
  BoolExpr invalid0;
  BoolExpr invalid1;

  SequentialDesignModel model0;
  model0.stateBits = {a0, b0, c0};
  model0.inputVarByKey.emplace(a0, 2);
  model0.inputVarByKey.emplace(b0, 3);
  model0.inputVarByKey.emplace(c0, 4);
  model0.initialStateValueByKey.emplace(a0, false);
  model0.initialStateValueByKey.emplace(b0, true);
  model0.initialStateValueByKey.emplace(c0, false);
  model0.nextStateExprByStateKey.emplace(a0, BoolExpr::createFalse());
  model0.nextStateExprByStateKey.emplace(b0, BoolExpr::Var(99));
  model0.nextStateExprByStateKey.emplace(c0, &invalid0);
  model0.complementedStateRelations.push_back({a0, b0});

  SequentialDesignModel model1;
  model1.stateBits = {b1, c1, a1};
  model1.inputVarByKey.emplace(a1, 5);
  model1.inputVarByKey.emplace(b1, 6);
  model1.inputVarByKey.emplace(c1, 7);
  model1.initialStateValueByKey.emplace(a1, false);
  model1.initialStateValueByKey.emplace(b1, true);
  model1.initialStateValueByKey.emplace(c1, false);
  model1.nextStateExprByStateKey.emplace(a1, BoolExpr::createFalse());
  model1.nextStateExprByStateKey.emplace(b1, BoolExpr::Var(123));
  model1.nextStateExprByStateKey.emplace(c1, &invalid1);
  model1.complementedStateRelations.push_back({a1, b1});

  const auto aligned = inferStructurallyEquivalentStatePairs(
      model0, model1, AlignedSignals{});

  EXPECT_EQ(aligned.names.size(), 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       StructuralStateInvariantReturnsEmptyWhenOneSideHasNoState) {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  model1.stateBits = {makeSignalKey("s1")};
  model1.inputVarByKey.emplace(model1.stateBits.front(), 2);
  model1.nextStateExprByStateKey.emplace(model1.stateBits.front(), BoolExpr::createFalse());

  const auto aligned = inferStructurallyEquivalentStatePairs(
      model0, model1, AlignedSignals{});

  EXPECT_TRUE(aligned.names.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       StructuralStateInvariantReturnsEmptyWhenBothSidesHaveNoState) {
  SequentialDesignModel model0;
  SequentialDesignModel model1;

  const auto aligned = inferStructurallyEquivalentStatePairs(
      model0, model1, AlignedSignals{});

  EXPECT_TRUE(aligned.names.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       StructuralStateInvariantHandlesNullNextStateExpressions) {
  const SignalKey null0 = makeSignalKey("null0");
  const SignalKey false0 = makeSignalKey("false0");
  const SignalKey null1 = makeSignalKey("null1");
  const SignalKey false1 = makeSignalKey("false1");

  SequentialDesignModel model0;
  model0.stateBits = {null0, false0};
  model0.inputVarByKey.emplace(null0, 2);
  model0.inputVarByKey.emplace(false0, 3);
  model0.nextStateExprByStateKey.emplace(null0, nullptr);
  model0.nextStateExprByStateKey.emplace(false0, BoolExpr::createFalse());

  SequentialDesignModel model1;
  model1.stateBits = {false1, null1};
  model1.inputVarByKey.emplace(false1, 4);
  model1.inputVarByKey.emplace(null1, 5);
  model1.nextStateExprByStateKey.emplace(false1, BoolExpr::createFalse());
  model1.nextStateExprByStateKey.emplace(null1, nullptr);

  const auto aligned = inferStructurallyEquivalentStatePairs(
      model0, model1, AlignedSignals{});

  EXPECT_EQ(aligned.names.size(), 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSupportsGenericComplementedStateNames) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNamedComplementSequentialModel(
      primitives, "DFF_STATE_STATEN", "STATE", "STATEN");
  auto* top = createSequentialOutputPairTop(
      library, "top", model, "STATE", "STATEN");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.complementedStateRelations.size(), 1u);
  EXPECT_EQ(
      extracted.complementedStateRelations.front().primaryKey,
      extracted.stateBits.front());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractTracksVectorStateBitsPerOutputTerm) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createBusSequentialModel(primitives, "DFF_BUS");
  auto* top = createBusSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.stateBits.size(), 2u);

  const auto in0Key = findKeyByDisplayName(extracted, "in[0]");
  const auto in1Key = findKeyByDisplayName(extracted, "in[1]");
  const auto q0Key = findKeyByDisplayName(extracted, "ff0.Q[0]");
  const auto q1Key = findKeyByDisplayName(extracted, "ff0.Q[1]");

  auto* q0Expr = extracted.nextStateExprByStateKey.at(q0Key);
  auto* q1Expr = extracted.nextStateExprByStateKey.at(q1Key);
  EXPECT_TRUE(q0Expr->evaluate({{extracted.inputVarByKey.at(in0Key), true},
                                {extracted.inputVarByKey.at(in1Key), false}}));
  EXPECT_FALSE(q1Expr->evaluate({{extracted.inputVarByKey.at(in0Key), true},
                                 {extracted.inputVarByKey.at(in1Key), false}}));
  EXPECT_FALSE(q0Expr->evaluate({{extracted.inputVarByKey.at(in0Key), false},
                                 {extracted.inputVarByKey.at(in1Key), true}}));
  EXPECT_TRUE(q1Expr->evaluate({{extracted.inputVarByKey.at(in0Key), false},
                                {extracted.inputVarByKey.at(in1Key), true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSkipsCombinationalInstancesWithoutStateOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top = createCombinationalInvTop(library, "top", invModel);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.stateBits.empty());
  EXPECT_EQ(extracted.observedOutputs.size(), 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractResolvesHierarchicalTopOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* childModel =
      createHierarchicalCombinationalInvChild(library, "child", invModel);
  auto* top =
      createHierarchicalCombinationalInvTop(library, "top", childModel);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.observedOutputs.size(), 1u);
  EXPECT_TRUE(extracted.skippedObservedOutputs.empty());

  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto outKey = findKeyByDisplayName(extracted, "out[0]");
  auto* expr = extracted.observedOutputExprByKey.at(outKey);

  EXPECT_TRUE(expr->evaluate({{extracted.inputVarByKey.at(inKey), false}}));
  EXPECT_FALSE(expr->evaluate({{extracted.inputVarByKey.at(inKey), true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractResolvesNestedHierarchicalTopOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* childModel =
      createHierarchicalCombinationalInvChild(library, "child", invModel);
  auto* midModel =
      createHierarchicalCombinationalInvTop(library, "mid", childModel);
  auto* top =
      createHierarchicalCombinationalInvTop(library, "top", midModel);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.observedOutputs.size(), 1u);
  EXPECT_TRUE(extracted.skippedObservedOutputs.empty());

  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto outKey = findKeyByDisplayName(extracted, "out[0]");
  auto* expr = extracted.observedOutputExprByKey.at(outKey);

  EXPECT_TRUE(expr->evaluate({{extracted.inputVarByKey.at(inKey), false}}));
  EXPECT_FALSE(expr->evaluate({{extracted.inputVarByKey.at(inKey), true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractPropagatesNoDriverSkipsToStateAndOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createPartialCoverageNoDriverTop(library, "top");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.totalObservedOutputCount(), 2u);
  EXPECT_EQ(extracted.coveredObservedOutputCount(), 1u);
  EXPECT_EQ(extracted.skippedObservedOutputs.size(), 1u);
  EXPECT_EQ(extracted.skippedStateBits.size(), 1u);
  EXPECT_NE(
      std::find(
          extracted.observedOutputs.begin(),
          extracted.observedOutputs.end(),
          findKeyByDisplayName(extracted, "good[0]")),
      extracted.observedOutputs.end());
  EXPECT_EQ(
      extracted.skippedObservedOutputs.front(),
      findKeyByDisplayName(extracted, "bad[0]"));
  const auto badKey = findKeyByDisplayName(extracted, "bad[0]");
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  ASSERT_NE(extracted.connectivitySkipInfoByKey.find(badKey),
            extracted.connectivitySkipInfoByKey.end());
  ASSERT_NE(extracted.connectivitySkipInfoByKey.find(stateKey),
            extracted.connectivitySkipInfoByKey.end());
  EXPECT_EQ(
      extracted.connectivitySkipInfoByKey.at(stateKey).origin,
      ConnectivitySkipOrigin::NoDriver);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractPropagatesMultiDriverSkipsToStateAndOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top = createPartialCoverageMultiDriverTop(library, "top", invModel);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.totalObservedOutputCount(), 2u);
  EXPECT_EQ(extracted.coveredObservedOutputCount(), 1u);
  EXPECT_EQ(extracted.skippedObservedOutputs.size(), 1u);
  EXPECT_EQ(extracted.skippedStateBits.size(), 1u);
  const auto badKey = findKeyByDisplayName(extracted, "bad[0]");
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  ASSERT_NE(extracted.connectivitySkipInfoByKey.find(badKey),
            extracted.connectivitySkipInfoByKey.end());
  ASSERT_NE(extracted.connectivitySkipInfoByKey.find(stateKey),
            extracted.connectivitySkipInfoByKey.end());
  EXPECT_EQ(
      extracted.connectivitySkipInfoByKey.at(stateKey).origin,
      ConnectivitySkipOrigin::MultiDriver);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSupportsSetHighInitialState) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createSetOnlySequentialModel(primitives, "DFF_SET");
  auto* top = createSetOnlySequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.stateBits.size(), 1u);
  EXPECT_TRUE(extracted.initialStateValueByKey.at(extracted.stateBits.front()));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractPreservesSetHighControlSemantics) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createSetOnlySequentialModel(primitives, "DFF_SET");
  auto* top = createSetOnlySequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);
  const auto stateKey = extracted.stateBits.front();
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto setKey = findKeyByDisplayName(extracted, "set[0]");
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(setKey), true}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(setKey), false}}));
  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(setKey), false}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSupportsResetHighInitialState) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createDffreTop(library, "top");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.stateBits.size(), 1u);
  EXPECT_FALSE(extracted.initialStateValueByKey.at(extracted.stateBits.front()));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractPreservesEnableAndResetControlSemantics) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createDffreTop(library, "top");

  const auto extracted = SequentialDesignModel::extract(top);
  const auto stateKey = extracted.stateBits.front();
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto enKey = findKeyByDisplayName(extracted, "en[0]");
  const auto rstKey = findKeyByDisplayName(extracted, "rst[0]");
  const size_t stateVar = extracted.inputVarByKey.at(stateKey);
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(enKey), true},
       {extracted.inputVarByKey.at(rstKey), true},
       {stateVar, true}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(enKey), true},
       {extracted.inputVarByKey.at(rstKey), false},
       {stateVar, false}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(enKey), false},
       {extracted.inputVarByKey.at(rstKey), false},
       {stateVar, true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractRejectsNullTop) {
  EXPECT_THROW(
      static_cast<void>(SequentialDesignModel::extract(nullptr)),
      std::invalid_argument);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractRejectsMissingUniverse) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top = createDffTop(library, "top", invModel, false, false);
  NLUniverse::get()->destroy();

  EXPECT_THROW(
      static_cast<void>(SequentialDesignModel::extract(top)),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractAbstractsUncomputableSequentialAsBoundaryByDefault) {
  ScopedSecBoundaryAbstraction boundaryAbstraction(true);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNoDataSequentialModel(primitives, "SEQ_NO_D");
  auto* top = createNoDataSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.stateBits.empty());
  EXPECT_NE(
      std::find(
          extracted.environmentInputs.begin(),
          extracted.environmentInputs.end(),
          findKeyByDisplayName(extracted, "ff0.Q[0]")),
      extracted.environmentInputs.end());
  EXPECT_NE(
      std::find(
          extracted.observedOutputs.begin(),
          extracted.observedOutputs.end(),
          findKeyByDisplayName(extracted, "ff0.CK[0]")),
      extracted.observedOutputs.end());
  ASSERT_EQ(extracted.abstractedSequentialBoundaries.size(), 1u);
  EXPECT_NE(
      extracted.abstractedSequentialBoundaries.front().find("ff0"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractAbstractsSequentialWithUnsupportedUpdatePinsAsBoundaryByDefault) {
  ScopedSecBoundaryAbstraction boundaryAbstraction(true);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createExtraUpdatePinSequentialModel(primitives, "SEQ_ADDR");
  auto* top = createExtraUpdatePinSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.stateBits.empty());
  EXPECT_FALSE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_NE(
      extracted.abstractedSequentialBoundaries.front().find("update pin `A`"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractReportsUnsupportedSequentialWithoutDInput) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNoDataSequentialModel(primitives, "SEQ_NO_D");
  auto* top = createNoDataSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  EXPECT_FALSE(extracted.unsupportedReasons.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractRejectsSequentialWithUnsupportedUpdatePinsInStrictMode) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createExtraUpdatePinSequentialModel(primitives, "SEQ_ADDR");
  auto* top = createExtraUpdatePinSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  ASSERT_FALSE(extracted.unsupportedReasons.empty());
  EXPECT_NE(
      extracted.unsupportedReasons.front().find("update pin `A`"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractRejectsMultipleControlStyles) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createResetSetSequentialModel(primitives, "SEQ_RST_SET");
  auto* top = createResetSetSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  ASSERT_FALSE(extracted.unsupportedReasons.empty());
  EXPECT_NE(
      extracted.unsupportedReasons.front().find("multiple control styles"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractMirrorsComplementedInitialStateValue) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNamedComplementSetSequentialModel(
      primitives, "DFF_STATE_SET", "STATE", "STATEN");
  auto* top = createComplementedSetSequentialTop(
      library, "top", model, "STATE", "STATEN");

  const auto extracted = SequentialDesignModel::extract(top);

  ASSERT_EQ(extracted.stateBits.size(), 2u);
  ASSERT_EQ(extracted.initialStateValueByKey.size(), 2u);
  const auto& relation = extracted.complementedStateRelations.front();
  EXPECT_TRUE(extracted.initialStateValueByKey.at(relation.primaryKey));
  EXPECT_FALSE(extracted.initialStateValueByKey.at(relation.complementedKey));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractReportsSharedScalarDataForMultiOutputPrimitive) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNamedComplementSequentialModel(
      primitives, "DFF_STATE_ALT", "STATE", "ALT");
  auto* top = createSequentialOutputPairTop(
      library, "top", model, "STATE", "ALT");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  ASSERT_FALSE(extracted.unsupportedReasons.empty());
  EXPECT_NE(
      extracted.unsupportedReasons.front().find("multiple independent state outputs"),
      std::string::npos);
  for (const auto& reason : extracted.unsupportedReasons) {
    EXPECT_EQ(reason.find("Missing next-state relation"), std::string::npos);
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractReportsSharedScalarDataForQAndUnrelatedOutput) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNamedComplementSequentialModel(
      primitives, "DFF_Q_ALT", "Q", "ALT");
  auto* top = createSequentialOutputPairTop(
      library, "top", model, "Q", "ALT");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  ASSERT_FALSE(extracted.unsupportedReasons.empty());
  EXPECT_NE(
      extracted.unsupportedReasons.front().find("multiple independent state outputs"),
      std::string::npos);
  for (const auto& reason : extracted.unsupportedReasons) {
    EXPECT_EQ(reason.find("Missing next-state relation"), std::string::npos);
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractStopsBeforeConeBuildForUnsupportedPrimitiveInfo) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createNamedComplementSequentialModel(
      primitives, "DFF_STATE_ALT", "STATE", "ALT");
  auto* top = createSequentialOutputPairTop(
      library, "top", model, "STATE", "ALT");

  ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto extracted = SequentialDesignModel::extract(top);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_TRUE(extracted.hasUnsupportedFeatures());
  EXPECT_NE(
      stderrOutput.find("SEC diag: extract(top) collect begin"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: extract(top) early unsupported exit before build"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("SEC diag: extract(top) build begin"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractFindsPrimaryStateOutputWhenComplementAppearsFirst) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createComplementFirstSequentialModel(
      primitives, "DFF_STATEN_STATE", "STATE", "STATEN");
  auto* top = createSequentialOutputPairTop(
      library, "top", model, "STATE", "STATEN");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  ASSERT_EQ(extracted.complementedStateRelations.size(), 1u);
  EXPECT_EQ(
      extracted.complementedStateRelations.front().primaryKey,
      extracted.stateBits.front());
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedObservedOutputNamesAreReportedAsUnsupported) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedOutputDffTop(library, "top0", invModel, "out0");
  auto* top1 = createNamedOutputDffTop(library, "top1", invModel, "out1");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("Mismatched observed output sets"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedObservedOutputCountsAreReportedAsUnsupported) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedOutputDffTop(library, "top0", invModel, "out");
  auto* top1 = createExtraOutputDffTop(library, "top1", invModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("Mismatched observed output sets"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedInputNamesAreReportedAsUnsupported) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedInputDffTop(library, "top0", invModel, "in0");
  auto* top1 = createNamedInputDffTop(library, "top1", invModel, "in1");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(
      result.reason.find("Mismatched environment input sets"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedInputCountsAreReportedAsUnsupported) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedInputDffTop(library, "top0", invModel, "in");
  auto* top1 = createExtraInputDffTop(library, "top1");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(
      result.reason.find("Mismatched environment input sets"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       TooSmallBoundRemainsInconclusiveBeforeCounterexampleDepth) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(library, "top0", false);
  auto* top1 = createResetInitializedPipelineTop(library, "top1", true);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.bound, 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ZeroBoundRemainsInconclusiveForEquivalentSequentialDesigns) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(0);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DifferentResultIncludesCounterexampleTracebackDetails) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, true, false);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_NE(result.reason.find("Input trace:"), std::string::npos);
  EXPECT_NE(
      result.reason.find("Observed output mismatches at cycle"),
      std::string::npos);
  EXPECT_NE(
      result.reason.find("Traceback for first differing point"),
      std::string::npos);
  EXPECT_NE(
      result.reason.find("design0 cone to environment inputs"),
      std::string::npos);
  EXPECT_NE(result.reason.find("cone terms only in design1"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       UnsupportedReasonsFromBothDesignsAreJoined) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* unsupportedModel = createNamedComplementSequentialModel(
      primitives, "DFF_STATE_ALT", "STATE", "ALT");
  auto* top0 = createSequentialOutputPairTop(
      library, "top0", unsupportedModel, "STATE", "ALT");
  auto* top1 = createSequentialOutputPairTop(
      library, "top1", unsupportedModel, "STATE", "ALT");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_FALSE(result.reason.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       StrategyStopsBeforeSecondExtractionAndProofOnUnsupportedPrimitiveInfo) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* unsupportedModel = createNamedComplementSequentialModel(
      primitives, "DFF_STATE_ALT", "STATE", "ALT");
  auto* invModel = createInvModel(primitives);
  auto* top0 = createSequentialOutputPairTop(
      library, "top0", unsupportedModel, "STATE", "ALT");
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(stderrOutput.find("SEC diag: extracted design0"), std::string::npos);
  EXPECT_EQ(stderrOutput.find("SEC diag: extracted design1"), std::string::npos);
  EXPECT_EQ(stderrOutput.find("SEC diag: aligning inputs/outputs"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentDesignsCanProveSecOnCoveredOutputsOnlyAfterNoDriverSkipping) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createPartialCoverageNoDriverTop(library, "top0");
  auto* top1 = createPartialCoverageNoDriverTop(library, "top1");

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(result.skippedObservedOutputs.front().find("no-driver"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentDesignsCanProveSecOnCoveredOutputsOnlyAfterMultiDriverSkipping) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createPartialCoverageMultiDriverTop(library, "top0", invModel);
  auto* top1 = createPartialCoverageMultiDriverTop(library, "top1", invModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(result.skippedObservedOutputs.front().find("multi-driver"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       UnsupportedSequentialInterfacesCanBeAbstractedAsSecBoundariesByDefault) {
  ScopedSecBoundaryAbstraction boundaryAbstraction(true);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* unsupportedModel = createNoDataSequentialModel(primitives, "SEQ_NO_D");
  auto* top0 =
      createUnsupportedPrimitiveCoverageTop(library, "top0", unsupportedModel);
  auto* top1 =
      createUnsupportedPrimitiveCoverageTop(library, "top1", unsupportedModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_FALSE(result.abstractedSequentialBoundaries.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       UnsupportedPrimitiveInformationStillFailsEvenWithOtherCoveredOutputs) {
  ScopedSecBoundaryAbstraction strictSequentialModeling(false);
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* unsupportedModel = createNoDataSequentialModel(primitives, "SEQ_NO_D");
  auto* top0 =
      createUnsupportedPrimitiveCoverageTop(library, "top0", unsupportedModel);
  auto* top1 =
      createUnsupportedPrimitiveCoverageTop(library, "top1", unsupportedModel);

  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_TRUE(result.skippedObservedOutputs.empty());
  EXPECT_FALSE(result.reason.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       DISABLED_DiagnosticModePrintsStrategyAndExtractionProgress) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  // This regression is about the diagnostic text itself, not the parallel
  // cloud builder. Keep the run single-threaded so stderr capture stays
  // deterministic under GoogleTest on macOS.
  ScopedEnvVar noMt("KEPLER_NO_MT", "1");
  ScopedFdCapture stderrCapture(STDERR_FILENO);
  SequentialEquivalenceStrategy strategy(top0, top1);
  const auto result = strategy.run(3);
  const std::string stderrOutput = stderrCapture.finish();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_NE(stderrOutput.find("SEC diag: start run"), std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: extract(top0) collect begin"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: remapped next-state formulas"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: entering legacy engine"),
      std::string::npos);
  EXPECT_NE(stderrOutput.find("SEC diag: aligned_inputs="), std::string::npos);
  EXPECT_NE(stderrOutput.find("SEC diag: property_is_true="), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapInductionProvesPostResetInvariant) {
  KInductionProblem problem;
  problem.environmentInputNames = {"rst"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.resetBootstrapInputs = {{2, true}};
  problem.bootstrapStateEqualityPairs = {{3, 4}};
  problem.inductiveStateEqualityPairs = {{3, 4}};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(3)}};
  problem.transitions1 = {{4, BoolExpr::Var(4)}};
  problem.property =
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(4)));
  problem.bad = BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(4));
  problem.description = "bootstrap induction regression";

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       StrongerInductionInvariantClosesOutputOnlySecAtOneStep) {
  KInductionProblem problem;
  problem.observedOutputNames = {"out"};
  problem.state0Symbols = {2, 3};
  problem.state1Symbols = {4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{2, BoolExpr::Var(3)}, {3, BoolExpr::Var(3)}};
  problem.transitions1 = {{4, BoolExpr::Var(5)}, {5, BoolExpr::Var(5)}};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)), BoolExpr::Not(BoolExpr::Var(5))));
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.property =
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4)));
  problem.bad = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4));
  problem.description = "output-only SEC needs a stronger invariant";

  KInductionEngine withoutInvariant(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto withoutInvariantResult = withoutInvariant.run(1);
  EXPECT_EQ(withoutInvariantResult.status, KInductionStatus::Inconclusive);

  problem.inductionProperty = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(4))),
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(5))));
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  KInductionEngine withInvariant(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto withInvariantResult = withInvariant.run(1);
  EXPECT_EQ(withInvariantResult.status, KInductionStatus::Equivalent);
  EXPECT_EQ(withInvariantResult.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelDetailHelpersCoverNextStateAndInitErrors) {
  const std::unordered_map<naja::DNL::DNLID, BoolExpr*> outputExprByTerm = {
      {11, BoolExpr::Var(7)},
      {12, BoolExpr::Var(8)},
      {13, BoolExpr::Var(9)},
  };

  EXPECT_THROW(
      detail::buildNextStateExprForTest(5, {{"D", 11}}, {2, 3}, outputExprByTerm),
      std::runtime_error);
  EXPECT_THROW(
      detail::buildNextStateExprForTest(0, {{"D", 11}}, {1}, outputExprByTerm),
      std::runtime_error);
  EXPECT_THROW(
      detail::buildNextStateExprForTest(0, {}, {2}, outputExprByTerm),
      std::runtime_error);
  EXPECT_THROW(
      detail::buildNextStateExprForTest(
          0, {{"D", 11}, {"R", 12}, {"S", 13}}, {2}, outputExprByTerm),
      std::runtime_error);
  auto* holdExpr = detail::buildNextStateExprForTest(
      0, {{"D", 11}, {"E", 12}, {"RN", 13}}, {2}, outputExprByTerm);
  EXPECT_FALSE(holdExpr->evaluate({{2, true}, {7, true}, {8, true}, {9, false}}));
  EXPECT_TRUE(holdExpr->evaluate({{2, false}, {7, true}, {8, true}, {9, true}}));
  EXPECT_TRUE(holdExpr->evaluate({{2, true}, {7, false}, {8, false}, {9, true}}));
  auto* setExpr = detail::buildNextStateExprForTest(
      0, {{"D", 11}, {"S", 12}}, {2}, outputExprByTerm);
  EXPECT_TRUE(setExpr->evaluate({{2, false}, {7, false}, {8, true}}));
  EXPECT_FALSE(setExpr->evaluate({{2, false}, {7, false}, {8, false}}));

  EXPECT_EQ(
      detail::detectInitialStateValueForTest({{"R", 11}}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::detectInitialStateValueForTest({{"RN", 11}}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::detectInitialStateValueForTest({{"S", 11}}),
      std::optional<bool>(true));
  EXPECT_EQ(detail::detectInitialStateValueForTest({}), std::nullopt);
  EXPECT_THROW(
      detail::detectInitialStateValueForTest({{"R", 11}, {"S", 12}}),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelDetailResetInferenceAndReachableStateHelpersCoverBranches) {
  const auto requiredOutputs = detail::selectRequiredBuilderOutputsForTest(
      {10, 11, 12, 13, 14},
      {10, 14},
      {12, 13, 13},
      {14});
  EXPECT_EQ(
      requiredOutputs,
      (std::vector<naja::DNL::DNLID>{10, 12, 13}));

  EXPECT_EQ(
      detail::getResetAssertionValueForTest("rst[0]"),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("rst_n[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(detail::getResetAssertionValueForTest("enable[0]"), std::nullopt);

  const auto shared = BoolExpr::Not(BoolExpr::Var(3));
  EXPECT_EQ(detail::evaluateConstantUnderAssignmentsForTest(nullptr, {}), std::nullopt);
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(BoolExpr::Var(1), {}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(BoolExpr::Var(0), {}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::And(shared, shared), {{3, false}}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::And(BoolExpr::createFalse(), BoolExpr::Var(99)), {}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, false}, {4, true}}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, true}, {4, false}}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, true}, {4, true}}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Or(BoolExpr::createTrue(), BoolExpr::Var(99)), {}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, true}, {4, false}}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, false}, {4, true}}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, false}, {4, false}}),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, true}, {4, false}}),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(
          BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::Var(4)), {{3, true}}),
      std::nullopt);
  BoolExpr invalidExpr;
  EXPECT_EQ(
      detail::evaluateConstantUnderAssignmentsForTest(&invalidExpr, {}),
      std::nullopt);

  const auto rstKey = makeSignalKey("rst");
  const auto stateAKey = makeSignalKey("state_a");
  const auto stateBKey = makeSignalKey("state_b");
  const auto stateAComplementKey = makeSignalKey("state_a_n");

  SequentialDesignModel inferredModel;
  inferredModel.environmentInputs = {rstKey};
  inferredModel.stateBits = {stateAKey, stateBKey, stateAComplementKey};
  inferredModel.displayNameByKey[rstKey] = "rst[0]";
  inferredModel.inputVarByKey[rstKey] = 10;
  inferredModel.inputVarByKey[stateAKey] = 2;
  inferredModel.inputVarByKey[stateBKey] = 3;
  inferredModel.inputVarByKey[stateAComplementKey] = 4;
  inferredModel.nextStateExprByStateKey[stateAKey] = BoolExpr::Var(10);
  inferredModel.nextStateExprByStateKey[stateBKey] =
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::createTrue());
  inferredModel.nextStateExprByStateKey[stateAComplementKey] =
      BoolExpr::Not(BoolExpr::Var(2));
  inferredModel.complementedStateRelations.push_back(
      {stateAKey, stateAComplementKey});

  detail::inferSynthesizedResetInitialStateValuesForTest(inferredModel);
  EXPECT_EQ(
      inferredModel.initialStateValueByKey.at(stateAKey),
      true);
  EXPECT_EQ(
      inferredModel.initialStateValueByKey.at(stateBKey),
      true);
  EXPECT_EQ(
      inferredModel.initialStateValueByKey.at(stateAComplementKey),
      false);

  const auto missingDisplayResetKey = makeSignalKey("rst_missing_display");
  const auto missingVarResetKey = makeSignalKey("rst_missing_var");
  const auto nullStateKey = makeSignalKey("null_state");
  const auto derivedStateKey = makeSignalKey("derived_state");
  const auto partnerPrimaryKey = makeSignalKey("partner_primary");
  const auto partnerComplementKey = makeSignalKey("partner_complement");

  SequentialDesignModel edgeCaseModel;
  edgeCaseModel.environmentInputs = {missingDisplayResetKey, missingVarResetKey, rstKey};
  edgeCaseModel.stateBits = {
      nullStateKey, derivedStateKey, partnerPrimaryKey, partnerComplementKey};
  edgeCaseModel.displayNameByKey[missingVarResetKey] = "rst[0]";
  edgeCaseModel.displayNameByKey[rstKey] = "rst[0]";
  edgeCaseModel.inputVarByKey[missingDisplayResetKey] = 30;
  edgeCaseModel.inputVarByKey[rstKey] = 31;
  edgeCaseModel.inputVarByKey[nullStateKey] = 2;
  edgeCaseModel.inputVarByKey[derivedStateKey] = 3;
  edgeCaseModel.inputVarByKey[partnerPrimaryKey] = 4;
  edgeCaseModel.inputVarByKey[partnerComplementKey] = 5;
  auto* sharedResetVar = BoolExpr::Var(31);
  edgeCaseModel.nextStateExprByStateKey[nullStateKey] = nullptr;
  edgeCaseModel.nextStateExprByStateKey[derivedStateKey] = BoolExpr::And(
      sharedResetVar, BoolExpr::Or(BoolExpr::Var(99), sharedResetVar));
  edgeCaseModel.nextStateExprByStateKey[partnerPrimaryKey] = BoolExpr::createFalse();
  edgeCaseModel.nextStateExprByStateKey[partnerComplementKey] = BoolExpr::createTrue();
  edgeCaseModel.initialStateValueByKey[partnerPrimaryKey] = false;
  edgeCaseModel.complementedStateRelations.push_back(
      {partnerPrimaryKey, partnerComplementKey});

  detail::inferSynthesizedResetInitialStateValuesForTest(edgeCaseModel);
  EXPECT_TRUE(edgeCaseModel.initialStateValueByKey.at(derivedStateKey));
  EXPECT_TRUE(edgeCaseModel.initialStateValueByKey.at(partnerComplementKey));

  const auto dependencyKnownKey = makeSignalKey("dependency_known");
  const auto dependencyDerivedKey = makeSignalKey("dependency_derived");
  SequentialDesignModel dependencyModel;
  dependencyModel.environmentInputs = {rstKey};
  dependencyModel.stateBits = {dependencyKnownKey, dependencyDerivedKey};
  dependencyModel.displayNameByKey[rstKey] = "rst[0]";
  dependencyModel.inputVarByKey[rstKey] = 40;
  dependencyModel.inputVarByKey[dependencyKnownKey] = 2;
  dependencyModel.inputVarByKey[dependencyDerivedKey] = 3;
  dependencyModel.initialStateValueByKey[dependencyKnownKey] = true;
  auto* sharedStateExpr = BoolExpr::Var(2);
  dependencyModel.nextStateExprByStateKey[dependencyKnownKey] = sharedStateExpr;
  dependencyModel.nextStateExprByStateKey[dependencyDerivedKey] = BoolExpr::And(
      sharedStateExpr,
      BoolExpr::Or(BoolExpr::Var(99), sharedStateExpr));

  detail::inferSynthesizedResetInitialStateValuesForTest(dependencyModel);
  EXPECT_TRUE(dependencyModel.initialStateValueByKey.at(dependencyDerivedKey));

  const auto derivedKey0 = makeSignalKey("derived0");
  const auto derivedKey1 = makeSignalKey("derived1");
  const auto xorKey = makeSignalKey("derived_xor");
  SequentialDesignModel bootstrapModel0;
  bootstrapModel0.environmentInputs = {rstKey};
  bootstrapModel0.stateBits = {derivedKey0, derivedKey1, xorKey};
  bootstrapModel0.displayNameByKey[rstKey] = "rst[0]";
  bootstrapModel0.inputVarByKey[rstKey] = 20;
  bootstrapModel0.inputVarByKey[derivedKey0] = 2;
  bootstrapModel0.inputVarByKey[derivedKey1] = 3;
  bootstrapModel0.inputVarByKey[xorKey] = 4;
  bootstrapModel0.initialStateValueByKey[derivedKey0] = true;
  bootstrapModel0.initialStateValueByKey[derivedKey1] = false;
  bootstrapModel0.nextStateExprByStateKey[derivedKey0] = BoolExpr::Var(2);
  bootstrapModel0.nextStateExprByStateKey[derivedKey1] = BoolExpr::Var(3);
  bootstrapModel0.nextStateExprByStateKey[xorKey] =
      BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));

  const auto bootstrapValues =
      detail::deriveResetBootstrapStateValuesForTest(bootstrapModel0, 1);
  EXPECT_EQ(bootstrapValues.at(xorKey), true);

  SequentialDesignModel bootstrapModel1 = bootstrapModel0;
  bootstrapModel1.initialStateValueByKey[derivedKey1] = true;

  AlignedSignals candidateStates;
  candidateStates.names = {"state_a", "state_b"};
  candidateStates.keys0 = {derivedKey0, derivedKey1};
  candidateStates.keys1 = {derivedKey0, derivedKey1};
  const auto anchoredStates = detail::filterStateEqualitiesByInitialValueForTest(
      bootstrapModel0, bootstrapModel1, candidateStates);
  ASSERT_EQ(anchoredStates.names.size(), 1u);
  EXPECT_EQ(anchoredStates.names.front(), "state_a");
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialEquivalenceStrategyDetailFormattingHelpersCoverFallbackPaths) {
  EXPECT_EQ(detail::formatStringListForTest({}, 3), "<none>");
  EXPECT_EQ(
      detail::formatStringListForTest({"a", "b", "c"}, 2),
      "a, b, ... +1 more");

  EXPECT_NE(
      detail::formatConeLevelsForTest({}).find("<no traced cone terms>"),
      std::string::npos);

  std::vector<std::vector<std::string>> levels;
  for (size_t level = 0; level < 13; ++level) {
    std::vector<std::string> levelTerms;
    for (size_t term = 0; term < 13; ++term) {
      levelTerms.push_back(
          "term_" + std::to_string(level) + "_" + std::to_string(term));
    }
    levels.push_back(std::move(levelTerms));
  }
  const auto formattedLevels = detail::formatConeLevelsForTest(levels);
  EXPECT_NE(formattedLevels.find("... +1 more trace steps"), std::string::npos);
  EXPECT_NE(formattedLevels.find("... +1 more"), std::string::npos);

  KInductionResult noWitnessResult;
  EXPECT_TRUE(detail::formatCounterexampleWitnessForTest(
                  noWitnessResult, {}, {}, nullptr, nullptr)
                  .empty());

  KInductionResult emptyTraceResult;
  emptyTraceResult.witness = KInductionResult::CounterexampleWitness{
      .badFrame = 2,
      .inputTrace = {},
      .outputMismatches = {{"ghost[0]", true, false}},
  };
  const auto emptyTraceText = detail::formatCounterexampleWitnessForTest(
      emptyTraceResult, {}, {}, nullptr, nullptr);
  EXPECT_NE(emptyTraceText.find("Input trace: <none>"), std::string::npos);
  EXPECT_NE(emptyTraceText.find("Cone traceback unavailable:"), std::string::npos);

  KInductionResult noEnvTraceResult;
  noEnvTraceResult.witness = KInductionResult::CounterexampleWitness{
      .badFrame = 1,
      .inputTrace = {{{0, {}}}},
      .outputMismatches = {{"ghost[0]", false, true}},
  };
  const auto noEnvTraceText = detail::formatCounterexampleWitnessForTest(
      noEnvTraceResult, {}, {}, nullptr, nullptr);
  EXPECT_NE(
      noEnvTraceText.find("<no environment inputs>"),
      std::string::npos);

  KInductionResult noMismatchTraceResult;
  noMismatchTraceResult.witness = KInductionResult::CounterexampleWitness{
      .badFrame = 1,
      .inputTrace = {{{0, {{"in[0]", true}}}}},
      .outputMismatches = {},
  };
  const auto noMismatchTraceText = detail::formatCounterexampleWitnessForTest(
      noMismatchTraceResult, {}, {}, nullptr, nullptr);
  EXPECT_EQ(
      noMismatchTraceText.find("Traceback for first differing point"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialEquivalenceStrategyDetailAlignmentAndLookupExposeErrors) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  const auto model0 = SequentialDesignModel::extract(top0);
  const auto model1 = SequentialDesignModel::extract(top1);

  KInductionResult missingSignalResult;
  missingSignalResult.witness = KInductionResult::CounterexampleWitness{
      .badFrame = 1,
      .inputTrace = {{{0, {{"in[0]", true}}}}},
      .outputMismatches = {{"ghost[0]", false, true}},
  };
  const auto missingSignalText = detail::formatCounterexampleWitnessForTest(
      missingSignalResult, model0, model1, top0, top1);
  EXPECT_NE(
      missingSignalText.find("design0 could not resolve the differing SEC signal back into the DNL"),
      std::string::npos);
  EXPECT_NE(
      missingSignalText.find("design1 could not resolve the differing SEC signal back into the DNL"),
      std::string::npos);

  KInductionResult topInputTraceResult;
  topInputTraceResult.witness = KInductionResult::CounterexampleWitness{
      .badFrame = 1,
      .inputTrace = {{{0, {{"in[0]", true}}}}},
      .outputMismatches = {{"in[0]", false, true}},
  };
  const auto topInputTraceText = detail::formatCounterexampleWitnessForTest(
      topInputTraceResult, model0, model1, top0, top1);
  EXPECT_NE(
      topInputTraceText.find("Observed output mismatches at cycle 1:"),
      std::string::npos);

  const auto inputKey = makeSignalKey("in");
  const auto outputKey = makeSignalKey("out");
  std::unordered_map<SignalKey, std::string, SignalKeyHash> displayNames0 = {
      {inputKey, "in[0]"},
  };
  std::unordered_map<SignalKey, std::string, SignalKeyHash> displayNames1 = {
      {inputKey, "in[0]"},
  };

  EXPECT_THROW(
      detail::alignSignalsByNameForTest(
          {inputKey}, displayNames0, {outputKey}, displayNames1, "inputs"),
      std::runtime_error);
  EXPECT_THROW(
      detail::alignSignalsByNameForTest(
          {inputKey, inputKey}, displayNames0, {inputKey}, displayNames1, "inputs"),
      std::runtime_error);

  auto* universe = NLUniverse::get();
  ASSERT_NE(universe, nullptr);
  universe->setTopDesign(top0);
  auto* dnl = naja::DNL::get();
  ASSERT_NE(dnl, nullptr);

  std::optional<naja::DNL::DNLID> outTermID;
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull() || !term.isTopPort()) {
      continue;
    }
    if (term.getSnlBitTerm()->getName().getString() == "out") {
      outTermID = termID;
      break;
    }
  }
  ASSERT_TRUE(outTermID.has_value());
  const auto outKey = detail::getTerminalPathKeyForTest(
      dnl->getDNLTerminalFromID(*outTermID));
  EXPECT_EQ(detail::findTermByKeyForTest(dnl, outKey), outTermID);

  SignalKey missingKey = outKey;
  ++missingKey.first.front();
  EXPECT_FALSE(detail::findTermByKeyForTest(dnl, missingKey).has_value());
}
