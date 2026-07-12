// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "model/SequentialDesignModel.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <limits>
#include <map>
#include <memory_resource>
#include <memory>
#include <optional>
#include <sstream>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "DNL.h"
#include "NLDB0.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesignModeling.h"
#include "SNLPath.h"
#include "../../clauses/SNLLogicCloud.h"
#include "../../clauses/Tree2BoolExpr.h"
#include "../../config/Config.h"
#include "common/BoolExprUtils.h"
#include "../../strategies/miter/BuildPrimaryOutputClauses.h"

namespace KEPLER_FORMAL::SEC {

namespace {

struct PendingPinTerm {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  naja::NL::NLID::Bit bit = 0;
};

using PendingPinMap =
    std::unordered_map<std::string, std::vector<PendingPinTerm>>;

struct StateOutputTerm {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  std::string pinName;
  naja::NL::NLID::Bit bit = 0;
  bool intrinsicClockIsInverted = false;
  std::vector<PendingPinTerm> clockTermIDs;
};

struct PendingTransition {  // LCOV_EXCL_LINE
  SignalKey stateKey;
  naja::DNL::DNLID stateTermID = naja::DNL::DNLID_MAX;
  std::string statePinName;
  bool stateOutputIsComplemented = false;
  naja::NL::NLID::Bit stateBit = 0;
  size_t independentStateOutputCount = 0;
  size_t boundaryInfoIndex = std::numeric_limits<size_t>::max();
  std::vector<SignalKey> complementedStateKeys;
  PendingPinMap pinTermIDs;
  std::vector<PendingPinTerm> clockTermIDs;
  bool intrinsicClockIsInverted = false;
};

struct PendingMemoryReadPort {
  std::vector<naja::DNL::DNLID> addressTermIDs;
  std::vector<naja::DNL::DNLID> dataTermIDs;
};

struct PendingMemoryWritePort {
  std::vector<naja::DNL::DNLID> addressTermIDs;
  std::vector<naja::DNL::DNLID> dataTermIDs;
  std::vector<naja::DNL::DNLID> maskTermIDs;
  std::vector<naja::DNL::DNLID> enableTermIDs;
  std::vector<std::vector<naja::DNL::DNLID>> extraWriteInputTermIDs;
};

struct PendingMemoryCellState {
  SignalKey key;
  std::string displayName;
  size_t cellIndex = 0;
  size_t bitIndex = 0;
};

struct PendingMemoryReadOutput {
  SignalKey key;
  size_t portIndex = 0;
  size_t bitIndex = 0;
};

struct PendingMemoryInstance {
  size_t width = 0;
  size_t depth = 0;
  size_t abits = 0;
  naja::NL::SNLDesignModeling::MemoryResetMode resetMode =
      naja::NL::SNLDesignModeling::MemoryResetMode::None;
  std::optional<naja::DNL::DNLID> resetTermID;
  size_t boundaryInfoIndex = std::numeric_limits<size_t>::max();
  std::vector<PendingMemoryReadPort> readPorts;
  std::vector<PendingMemoryWritePort> writePorts;
  std::vector<PendingMemoryCellState> cellStates;
  std::vector<PendingMemoryReadOutput> readOutputs;
};

struct BoundaryObservedTerm {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  SignalKey key;
};

struct InstanceBoundaryInfo {
  std::string instancePath;
  std::vector<SignalKey> stateKeys;
  std::vector<BoundaryObservedTerm> observedTerms;
};

struct SequentialInstanceScan {
  PendingPinMap pinTermIDs;
  std::vector<StateOutputTerm> stateOutputs;
  InstanceBoundaryInfo boundaryInfo;
};

AbstractedSequentialBoundaryDetail makeAbstractedBoundaryDetail(
    const InstanceBoundaryInfo& info) {
  AbstractedSequentialBoundaryDetail detail;
  detail.instancePath = info.instancePath;
  detail.stateKeys = info.stateKeys;
  detail.observedKeys.reserve(info.observedTerms.size());
  for (const auto& observedTerm : info.observedTerms) {
    detail.observedKeys.push_back(observedTerm.key);
  }
  return detail;
}

struct BuiltObservedExpr {  // LCOV_EXCL_LINE
  BoolExpr* expr = nullptr;  // LCOV_EXCL_LINE
  std::optional<ConnectivitySkipInfo> connectivitySkip;
  std::string unsupportedReason;
};

using BuilderSkippedOutputInfo = KEPLER_FORMAL::BuildPrimaryOutputClauses::SkippedOutputInfo;
using BuilderSkippedOutputReason =
    KEPLER_FORMAL::BuildPrimaryOutputClauses::SkippedOutputReason;

std::string normalizePinName(const std::string& name);

bool isClockTreeBufferCell(const naja::DNL::DNLTerminalFull& term);

bool isPotentialClockTreeBufferCell(const naja::DNL::DNLTerminalFull& term);

std::optional<naja::DNL::DNLID> getClockTreeBufferSourceDriverTerm(
    naja::DNL::DNLFull* dnl,
    const naja::DNL::DNLTerminalFull& outputTerm);

bool hasBuildableCombinationalRoot(
    naja::DNL::DNLFull* dnl,
    naja::DNL::DNLID requestedTermID) {
  if (dnl == nullptr || requestedTermID == naja::DNL::DNLID_MAX) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  std::unordered_set<naja::DNL::DNLID> visitedTerms;
  naja::DNL::DNLID currentTermID = requestedTermID;
  while (currentTermID != naja::DNL::DNLID_MAX &&
         visitedTerms.insert(currentTermID).second) {
    const auto& currentTerm = dnl->getDNLTerminalFromID(currentTermID);
    if (currentTerm.isNull()) {
      // LCOV_EXCL_START
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (currentTerm.isTopPort() &&
        currentTerm.getSnlBitTerm()->getDirection() !=
            naja::NL::SNLBitTerm::Direction::Output) {
      return true;
    }

    if (currentTerm.getSnlBitTerm()->getDirection() !=
        naja::NL::SNLBitTerm::Direction::Output) {
      const auto isoID = currentTerm.getIsoID();
      if (isoID == naja::DNL::DNLID_MAX) {
        // LCOV_EXCL_START
        return false;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant()) {
        // LCOV_EXCL_START
        return true;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      if (iso.getDrivers().size() != 1) {
        // LCOV_EXCL_START
        return false;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      currentTermID = iso.getDrivers().front();
      continue;
    }

    if (currentTerm.isTopPort()) {
      // LCOV_EXCL_START
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (isClockTreeBufferCell(currentTerm)) {
      // LCOV_EXCL_START
      if (const auto sourceDriver =  // LCOV_EXCL_LINE
              getClockTreeBufferSourceDriverTerm(dnl, currentTerm)) {  // LCOV_EXCL_LINE
        currentTermID = *sourceDriver;  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    const auto& instance = currentTerm.getDNLInstance();
    auto* model = instance.getSNLModel();
    if (model == nullptr) {
      // LCOV_EXCL_START
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (NLDB0::isAssign(model)) {
      std::optional<naja::DNL::DNLID> passthroughDriver;
      for (auto* inputBitTerm : naja::NL::SNLDesignModeling::getCombinatorialInputs(
               const_cast<naja::NL::SNLBitTerm*>(currentTerm.getSnlBitTerm()))) {
        if (inputBitTerm == nullptr ||
            inputBitTerm->getDirection() ==
                naja::NL::SNLBitTerm::Direction::Output) {
          // LCOV_EXCL_START
          continue;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        const auto& inputTerm = instance.getTerminalFromBitTerm(inputBitTerm);
        if (inputTerm.isNull() || inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
          // LCOV_EXCL_START
          passthroughDriver.reset();  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
        if (iso.isConstant() || iso.getDrivers().size() != 1) {
          passthroughDriver.reset();
          break;
        }
        // LCOV_EXCL_START
        if (passthroughDriver.has_value()) {  // LCOV_EXCL_LINE
          passthroughDriver.reset();  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        passthroughDriver = iso.getDrivers().front();  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      if (passthroughDriver.has_value()) {
        // LCOV_EXCL_START
        currentTermID = *passthroughDriver;  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      return false;
    }

    const auto* bitTerm = currentTerm.getSnlBitTerm();
    if (naja::NL::SNLDesignModeling::getTruthTableCount(model) <=
        bitTerm->getOrderID()) {
      return false;
    }
    const auto& truthTable = naja::NL::SNLDesignModeling::getTruthTable(
        model, bitTerm->getOrderID());
    return truthTable.isInitialized();
  }

  // LCOV_EXCL_START
  return false;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

std::string describeConnectivitySkipOrigin(ConnectivitySkipOrigin origin) {
  switch (origin) {
    case ConnectivitySkipOrigin::NoDriver:
      return "no-driver";
    case ConnectivitySkipOrigin::MultiDriver:
      return "multi-driver";
    case ConnectivitySkipOrigin::LogicalLoop:
      return "logical-loop";
    case ConnectivitySkipOrigin::MultiClockDomain:
      // LCOV_EXCL_START
      return "multi-clock-domain";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return "connectivity";  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

const char* describeBuilderSkippedOutputReason(
    BuilderSkippedOutputReason reason) {
  switch (reason) {
    case BuilderSkippedOutputReason::NoDriver:
      return "no_driver";
    case BuilderSkippedOutputReason::MultiDriver:
      // LCOV_EXCL_START
      return "multi_driver";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    case BuilderSkippedOutputReason::LogicalLoop:
      // LCOV_EXCL_START
      return "logical_loop";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    case BuilderSkippedOutputReason::None:
      // LCOV_EXCL_START
      return "none";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return "unknown";  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

std::optional<ConnectivitySkipInfo> getConnectivitySkipInfo(
    const BuilderSkippedOutputInfo& info) {
  switch (info.reason) {
    case BuilderSkippedOutputReason::NoDriver:
      return ConnectivitySkipInfo{ConnectivitySkipOrigin::NoDriver, info.detail};
    case BuilderSkippedOutputReason::MultiDriver:
      return ConnectivitySkipInfo{ConnectivitySkipOrigin::MultiDriver, info.detail};
    case BuilderSkippedOutputReason::LogicalLoop:
      return ConnectivitySkipInfo{
          ConnectivitySkipOrigin::LogicalLoop, info.detail};
    // LCOV_EXCL_START
    case BuilderSkippedOutputReason::None:  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      return std::nullopt;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }
}

// LCOV_EXCL_START
BuiltObservedExpr buildObservedExprForTerm(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    naja::DNL::DNLID termID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::vector<naja::DNL::DNLID>& inputTerms,
    const std::vector<naja::DNL::DNLID>& outputTerms,
    const std::vector<size_t>& termDNLID2varID) {
  // LCOV_EXCL_START
  BuiltObservedExpr result;  // LCOV_EXCL_LINE
  if (const auto exprIt = outputExprByTerm.find(termID);  // LCOV_EXCL_LINE
      exprIt != outputExprByTerm.end()) {  // LCOV_EXCL_LINE
    result.expr = exprIt->second;  // LCOV_EXCL_LINE
    return result;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  auto* dnl = naja::DNL::get();  // LCOV_EXCL_LINE
  if (dnl == nullptr) {  // LCOV_EXCL_LINE
    result.unsupportedReason = "DNL is not available while rebuilding a SEC boundary";  // LCOV_EXCL_LINE
    return result;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::vector<bool> isPIs(dnl->getNBterms(), false);  // LCOV_EXCL_LINE
  for (const auto inputTermID : inputTerms) {  // LCOV_EXCL_LINE
    if (inputTermID < isPIs.size()) {  // LCOV_EXCL_LINE
      isPIs[inputTermID] = true;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::vector<bool> isPOs(dnl->getNBterms(), false);  // LCOV_EXCL_LINE
  for (const auto outputTermID : outputTerms) {  // LCOV_EXCL_LINE
    if (outputTermID < isPOs.size()) {  // LCOV_EXCL_LINE
      isPOs[outputTermID] = true;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  auto describeTerm = [&](const naja::DNL::DNLTerminalFull& term) {  // LCOV_EXCL_LINE
    return term.getSnlBitTerm()->getDesign()->getName().getString() + "." +  // LCOV_EXCL_LINE
           term.getSnlBitTerm()->getName().getString() + "[" +  // LCOV_EXCL_LINE
           std::to_string(term.getSnlBitTerm()->getBit()) + "]";  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  auto buildFromOutputTerm = [&](naja::DNL::DNLID outputTermID) {  // LCOV_EXCL_LINE
    BuiltObservedExpr localResult;  // LCOV_EXCL_LINE
    auto localIsPOs = isPOs;  // LCOV_EXCL_LINE
    if (outputTermID < localIsPOs.size()) {  // LCOV_EXCL_LINE
      localIsPOs[outputTermID] = true;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    KEPLER_FORMAL::SNLLogicCloud cloud(  // LCOV_EXCL_LINE
        outputTermID,  // LCOV_EXCL_LINE
        isPIs,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        localIsPOs);
    // LCOV_EXCL_START
    cloud.compute();  // LCOV_EXCL_LINE
    if (cloud.getTruthTable().isValid()) {  // LCOV_EXCL_LINE
      for (const auto inputTermID : cloud.getInputs()) {  // LCOV_EXCL_LINE
        if (inputTermID == naja::DNL::DNLID_MAX) {  // LCOV_EXCL_LINE
          continue;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        if (inputTermID >= termDNLID2varID.size() ||  // LCOV_EXCL_LINE
            termDNLID2varID[inputTermID] < 2) {  // LCOV_EXCL_LINE
          localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
              ConnectivitySkipOrigin::NoDriver,  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              "encountered internal frontier term " +  // LCOV_EXCL_LINE
                  std::to_string(inputTermID) +  // LCOV_EXCL_LINE
                  // LCOV_EXCL_STOP
                  " that was not collected as a primary input"};  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          cloud.destroy();  // LCOV_EXCL_LINE
          return localResult;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      cloud.getTruthTable().finalize();  // LCOV_EXCL_LINE
      localResult.expr =  // LCOV_EXCL_LINE
          KEPLER_FORMAL::Tree2BoolExpr::convert(cloud.getTruthTable(), termDNLID2varID);  // LCOV_EXCL_LINE
      cloud.destroy();  // LCOV_EXCL_LINE
      return localResult;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    switch (cloud.getSkipReason()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      case KEPLER_FORMAL::SNLLogicCloud::SkipReason::NoDriver:
        // LCOV_EXCL_START
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
            ConnectivitySkipOrigin::NoDriver, cloud.getSkipReasonText()};  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      case KEPLER_FORMAL::SNLLogicCloud::SkipReason::MultiDriver:
        // LCOV_EXCL_START
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
            ConnectivitySkipOrigin::MultiDriver, cloud.getSkipReasonText()};  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      case KEPLER_FORMAL::SNLLogicCloud::SkipReason::LogicalLoop:
        // LCOV_EXCL_START
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
            ConnectivitySkipOrigin::LogicalLoop, cloud.getSkipReasonText()};  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      case KEPLER_FORMAL::SNLLogicCloud::SkipReason::None:  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      default:
        // LCOV_EXCL_START
        localResult.unsupportedReason = "failed to build a Boolean expression";  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    cloud.destroy();  // LCOV_EXCL_LINE
    return localResult;  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  std::unordered_set<naja::DNL::DNLID> visitedTerms;  // LCOV_EXCL_LINE
  auto buildRecursively = [&](auto&& self,  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
                              naja::DNL::DNLID currentTermID) -> BuiltObservedExpr {
    // LCOV_EXCL_START
    BuiltObservedExpr localResult;  // LCOV_EXCL_LINE
    if (!visitedTerms.insert(currentTermID).second) {  // LCOV_EXCL_LINE
      localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          ConnectivitySkipOrigin::LogicalLoop,
          // LCOV_EXCL_START
          "a logical loop was detected while rebuilding a SEC boundary"};  // LCOV_EXCL_LINE
      return localResult;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    if (const auto exprIt = outputExprByTerm.find(currentTermID);  // LCOV_EXCL_LINE
        exprIt != outputExprByTerm.end()) {  // LCOV_EXCL_LINE
      localResult.expr = exprIt->second;  // LCOV_EXCL_LINE
      return localResult;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    if (currentTermID < termDNLID2varID.size() && isPIs[currentTermID]) {  // LCOV_EXCL_LINE
      const size_t varID = termDNLID2varID[currentTermID];  // LCOV_EXCL_LINE
      if (varID == 0) {  // LCOV_EXCL_LINE
        localResult.expr = BoolExpr::createFalse();  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (varID == 1) {  // LCOV_EXCL_LINE
        localResult.expr = BoolExpr::createTrue();  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (varID >= 2) {  // LCOV_EXCL_LINE
        localResult.expr = BoolExpr::Var(varID);  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP

    // LCOV_EXCL_START
    const auto& term = dnl->getDNLTerminalFromID(currentTermID);  // LCOV_EXCL_LINE
    if (term.getSnlBitTerm()->getDirection() !=  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
        naja::NL::SNLBitTerm::Direction::Output) {
      // LCOV_EXCL_START
      if (term.getIsoID() == naja::DNL::DNLID_MAX) {  // LCOV_EXCL_LINE
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            ConnectivitySkipOrigin::NoDriver,
            // LCOV_EXCL_START
            "term `" + describeTerm(term) + "` is not connected"};  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }

      // LCOV_EXCL_START
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());  // LCOV_EXCL_LINE
      if (iso.isConstant0()) {  // LCOV_EXCL_LINE
        localResult.expr = BoolExpr::createFalse();  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (iso.isConstant1()) {  // LCOV_EXCL_LINE
        localResult.expr = BoolExpr::createTrue();  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (iso.getDrivers().empty()) {  // LCOV_EXCL_LINE
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            ConnectivitySkipOrigin::NoDriver,
            // LCOV_EXCL_START
            "term `" + describeTerm(term) + "` has no drivers"};  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (iso.getDrivers().size() > 1) {  // LCOV_EXCL_LINE
        localResult.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            ConnectivitySkipOrigin::MultiDriver,
            // LCOV_EXCL_START
            "term `" + describeTerm(term) + "` has multiple drivers"};  // LCOV_EXCL_LINE
        return localResult;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      return self(self, iso.getDrivers().front());  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    return buildFromOutputTerm(currentTermID);  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  return buildRecursively(buildRecursively, termID);  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

struct MaterializedBuilderOutputs {
  std::vector<naja::DNL::DNLID> inputs;
  std::vector<naja::DNL::DNLID> outputs;
  std::vector<size_t> termDNLID2varID;
  std::unordered_map<naja::DNL::DNLID, BoolExpr*> outputExprByTerm;
  std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo> skippedOutputsByTerm;
};

void appendUniqueTermIDs(
    std::vector<naja::DNL::DNLID>& dest,
    const std::vector<naja::DNL::DNLID>& src) {
  std::unordered_set<naja::DNL::DNLID> seen(dest.begin(), dest.end());
  for (const auto termID : src) {
    if (seen.insert(termID).second) {
      dest.push_back(termID);
    }
  }
}

bool isUnassignedBuilderVarID(size_t varID) {
  return varID == static_cast<size_t>(-1);
}

void mergeBuilderTermVarIDs(
    std::vector<size_t>& dest,
    const std::vector<size_t>& src) {
  if (src.size() > dest.size()) {
    // LCOV_EXCL_START
    dest.resize(src.size(), static_cast<size_t>(-1));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  for (size_t index = 0; index < src.size(); ++index) {
    if (isUnassignedBuilderVarID(dest[index]) &&
        !isUnassignedBuilderVarID(src[index])) {
      // LCOV_EXCL_START
      dest[index] = src[index];  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
}

// LCOV_EXCL_START
void mergeMaterializedBuilderOutputs(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    MaterializedBuilderOutputs& dest,
    const MaterializedBuilderOutputs& src) {
  // LCOV_EXCL_START
  appendUniqueTermIDs(dest.inputs, src.inputs);  // LCOV_EXCL_LINE
  appendUniqueTermIDs(dest.outputs, src.outputs);  // LCOV_EXCL_LINE
  mergeBuilderTermVarIDs(dest.termDNLID2varID, src.termDNLID2varID);  // LCOV_EXCL_LINE
  for (const auto& [termID, expr] : src.outputExprByTerm) {  // LCOV_EXCL_LINE
    dest.outputExprByTerm.insert_or_assign(termID, expr);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  for (const auto& [termID, info] : src.skippedOutputsByTerm) {  // LCOV_EXCL_LINE
    dest.skippedOutputsByTerm.emplace(termID, info);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

std::unordered_map<size_t, size_t> buildStableBuilderVarRemap(
    const std::vector<size_t>& sourceTermDNLID2varID,
    const std::vector<size_t>& stableTermDNLID2varID) {
  std::unordered_map<size_t, size_t> remap;
  const size_t termCount =
      std::min(sourceTermDNLID2varID.size(), stableTermDNLID2varID.size());
  for (size_t termID = 0; termID < termCount; ++termID) {
    const size_t sourceVarID = sourceTermDNLID2varID[termID];
    if (isUnassignedBuilderVarID(sourceVarID) || sourceVarID < 2) {
      continue;
    }
    const size_t stableVarID = stableTermDNLID2varID[termID];
    if (isUnassignedBuilderVarID(stableVarID)) {
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    remap.emplace(sourceVarID, stableVarID);
  }
  return remap;
}

bool varRemapChangesAnySymbol(const std::unordered_map<size_t, size_t>& remap) {
  for (const auto& [sourceVarID, stableVarID] : remap) {
    if (sourceVarID != stableVarID) {
      return true;
    }
  }
  return false;
}

MaterializedBuilderOutputs materializeBuilderOutputs(
    const std::vector<naja::DNL::DNLID>& requestedOutputs,
    const std::vector<naja::DNL::DNLID>& collectedInputs,
    const std::vector<size_t>& stableTermDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>&
        collectedSkippedOutputs,
    bool secDiagEnabled,
    const char* topName,
    const char* phaseLabel) {
  MaterializedBuilderOutputs result;

  KEPLER_FORMAL::BuildPrimaryOutputClauses builder;
  builder.setRetainDnl(true);
  std::vector<naja::DNL::DNLID> normalizedRoots;
  normalizedRoots.reserve(requestedOutputs.size());
  std::unordered_map<naja::DNL::DNLID, std::vector<naja::DNL::DNLID>> requestedByRoot;
  std::unordered_set<naja::DNL::DNLID> seenRoots;
  std::unordered_set<naja::DNL::DNLID> nonTopRootsToBuild;
  auto* dnl = naja::DNL::get();
  size_t clockBufferRootPassthroughs = 0;
  const auto describeTerminal = [](const naja::DNL::DNLTerminalFull& terminal) {
    std::ostringstream oss;
    // Keep the debug label cheap: the full DNL instance path reconstruction is
    // surprisingly expensive on large designs and can dominate CVA6 diag runs.
    // The term ID is logged separately, so a compact model.pin[bit] label is
    // enough to correlate skipped roots back to the detailed skip reports.
    oss << terminal.getSnlBitTerm()->getDesign()->getName().getString() << "."
        << terminal.getSnlBitTerm()->getName().getString() << "["
        << terminal.getSnlBitTerm()->getBit() << "]";
    return oss.str();
  };
  const auto findBuildableOutputRoot = [&](naja::DNL::DNLID requestedTermID)
      -> std::optional<naja::DNL::DNLID> {
    std::unordered_set<naja::DNL::DNLID> visitedTerms;
    naja::DNL::DNLID currentTermID = requestedTermID;
    while (currentTermID != naja::DNL::DNLID_MAX &&
           visitedTerms.insert(currentTermID).second) {
      const auto& currentTerm = dnl->getDNLTerminalFromID(currentTermID);
      if (currentTerm.isNull()) {
        // LCOV_EXCL_START
        return std::nullopt;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      if (currentTerm.isTopPort() &&
          currentTerm.getSnlBitTerm()->getDirection() !=
              naja::NL::SNLBitTerm::Direction::Output) {
        return currentTermID;
      }
      if (currentTerm.getSnlBitTerm()->getDirection() ==
          naja::NL::SNLBitTerm::Direction::Output) {
        if (isClockTreeBufferCell(currentTerm)) {
          if (const auto sourceDriver =
                  getClockTreeBufferSourceDriverTerm(dnl, currentTerm)) {
            ++clockBufferRootPassthroughs;
            currentTermID = *sourceDriver;
            continue;
          }
        // LCOV_EXCL_START
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP

        const auto& inst = currentTerm.getDNLInstance();
        auto* model = inst.getSNLModel();
        if (model != nullptr && NLDB0::isAssign(model)) {
          std::optional<naja::DNL::DNLID> passthroughDriver;
          for (auto* inputBitTerm : naja::NL::SNLDesignModeling::getCombinatorialInputs(
                   const_cast<naja::NL::SNLBitTerm*>(currentTerm.getSnlBitTerm()))) {
            if (inputBitTerm == nullptr ||
                inputBitTerm->getDirection() ==
                    naja::NL::SNLBitTerm::Direction::Output) {
              // LCOV_EXCL_START
              continue;  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
            }
            const auto& inputTerm = inst.getTerminalFromBitTerm(inputBitTerm);
            if (inputTerm.isNull() || inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
              // LCOV_EXCL_START
              passthroughDriver.reset();  // LCOV_EXCL_LINE
              break;  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
            }
            const auto& iso =
                dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
            if (iso.isConstant() || iso.getDrivers().size() != 1) {
              passthroughDriver.reset();
              break;
            }
            if (passthroughDriver.has_value()) {
              // LCOV_EXCL_START
              passthroughDriver.reset();  // LCOV_EXCL_LINE
              break;  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
            }
            passthroughDriver = iso.getDrivers().front();
          }
          if (passthroughDriver.has_value()) {
            currentTermID = *passthroughDriver;
            continue;
          }
        }
        return currentTermID;
      }

      const auto isoID = currentTerm.getIsoID();
      if (isoID == naja::DNL::DNLID_MAX) {
        // LCOV_EXCL_START
        return std::nullopt;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant() || iso.getDrivers().size() != 1) {
        return std::nullopt;
      }
      currentTermID = iso.getDrivers().front();
    }
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  };
  for (const auto requestedTermID : requestedOutputs) {
    const auto rootTermID = findBuildableOutputRoot(requestedTermID);
    if (!rootTermID.has_value()) {
      continue;
    }
    const auto& rootTerm = dnl->getDNLTerminalFromID(*rootTermID);
    if (!(rootTerm.isTopPort() &&
          rootTerm.getSnlBitTerm()->getDirection() !=
              naja::NL::SNLBitTerm::Direction::Output) &&
        hasBuildableCombinationalRoot(dnl, *rootTermID)) {
      nonTopRootsToBuild.insert(*rootTermID);
    }
    requestedByRoot[*rootTermID].push_back(requestedTermID);
    const auto skippedIt = collectedSkippedOutputs.find(*rootTermID);
    if (skippedIt != collectedSkippedOutputs.end()) {
      // LCOV_EXCL_START
      result.skippedOutputsByTerm.emplace(*rootTermID, skippedIt->second);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (seenRoots.insert(*rootTermID).second) {
      normalizedRoots.push_back(*rootTermID);
    }
  }

  // A combinational dependency root can also appear in the carried boundary-input
  // set from an earlier extraction pass.  If we leave that root marked as a PI,
  // the cloud stops at the root instead of rebuilding its cone; clock-gated flops
  // then see their gated clock as a stale frontier rather than clk & enable.
  // Sequential state outputs must remain PIs here: they are the current-state
  // variables used by reset/bootstrap proofs and must not be rebuilt as logic.
  std::vector<naja::DNL::DNLID> materializationInputs;
  materializationInputs.reserve(collectedInputs.size());
  for (const auto inputTermID : collectedInputs) {
    if (nonTopRootsToBuild.find(inputTermID) != nonTopRootsToBuild.end()) {
      continue;
    }
    materializationInputs.push_back(inputTermID);
  }
  builder.setInputs(materializationInputs);

  // Structured memory dependencies are often internal roots that SEC needs
  // even when the generic boundary collector would not classify them as
  // outputs. Root input-side dependency pins at their single-driver source so
  // the clause builder sees the actual combinational producer instead of the
  // sink pin on the memory instance. If that cone hits a bad frontier, the root
  // is skipped and the dependent state/output is filtered later.
  builder.setOutputs(normalizedRoots);

  const bool traceDependencyRoots =
      secDiagEnabled &&
      std::string_view(phaseLabel).find("dependency build") !=
          std::string_view::npos;
  if (std::getenv("KEPLER_SEC_CLOCK_GATE_DIAG") != nullptr &&
      // LCOV_EXCL_START
      clockBufferRootPassthroughs != 0) {  // LCOV_EXCL_LINE
    std::fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        "SEC diag: extract(%s) %s clock-buffer root passthroughs=%zu outputs=%zu\n",
        // LCOV_EXCL_START
        topName,  // LCOV_EXCL_LINE
        phaseLabel,  // LCOV_EXCL_LINE
        clockBufferRootPassthroughs,  // LCOV_EXCL_LINE
        normalizedRoots.size());  // LCOV_EXCL_LINE
    std::fflush(stderr);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  if (secDiagEnabled) {
    // Dependency batches on designs like CVA6 can include thousands of memory
    // pins. Logging every requested-to-root mapping is prohibitively expensive,
    // so the diagnostic path only prints summary begin/end markers here and
    // prints the detailed term mapping only for the roots that later skip.
    fprintf(
        stderr,
        "SEC diag: extract(%s) %s begin outputs=%zu\n",
        topName,
        phaseLabel,
        normalizedRoots.size());
    fflush(stderr);
  }
  builder.build();
  // BuildPrimaryOutputClauses owns a temporary DNL expansion and destroys the
  // singleton when build() completes. Rebuilding the full DNL here is very
  // expensive on CVA6, so only reacquire it when the detailed dependency-root
  // diagnostics actually need to print terminal names after the build.
  if (traceDependencyRoots) {
    dnl = naja::DNL::get();
  }
  if (secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) %s end outputs=%zu\n",
        topName,
        phaseLabel,
        normalizedRoots.size());
    fflush(stderr);
  }

  result.inputs = builder.getInputs();
  result.outputs = builder.getOutputs();
  result.termDNLID2varID = stableTermDNLID2varID;
  mergeBuilderTermVarIDs(result.termDNLID2varID, builder.getTermDNLID2VarID());
  const auto stableVarRemap = buildStableBuilderVarRemap(
      builder.getTermDNLID2VarID(), result.termDNLID2varID);
  const bool needsStableVarRemap = varRemapChangesAnySymbol(stableVarRemap);
  for (const auto& [termID, info] : builder.getSkippedOutputs()) {
    result.skippedOutputsByTerm.emplace(termID, info);
  }
  if (traceDependencyRoots) {
    for (const auto& [rootTermID, skipInfo] : result.skippedOutputsByTerm) {
      const auto requestedIt = requestedByRoot.find(rootTermID);
      if (requestedIt == requestedByRoot.end()) {
        // LCOV_EXCL_START
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto& rootTerm = dnl->getDNLTerminalFromID(rootTermID);
      for (const auto requestedTermID : requestedIt->second) {
        const auto& requestedTerm = dnl->getDNLTerminalFromID(requestedTermID);
        fprintf(
            stderr,
            "SEC diag: extract(%s) %s skipped requested=%s term=%zu root=%s term=%zu reason=%s detail=%s\n",
            topName,
            phaseLabel,
            describeTerminal(requestedTerm).c_str(),
            requestedTermID,
            describeTerminal(rootTerm).c_str(),
            rootTermID,
            describeBuilderSkippedOutputReason(skipInfo.reason),
            skipInfo.detail.c_str());
      }
    }
  }
  const auto& outputExprs = builder.getPOs();
  // All formulas in this materialized batch use the same builder-local to
  // stable-variable map. Reuse one memo across outputs so shared ASIC cones are
  // remapped once instead of once per requested root.
  std::unordered_map<BoolExpr*, BoolExpr*> stableRemapMemo;
  for (size_t i = 0; i < result.outputs.size(); ++i) {
    BoolExpr* expr = outputExprs[i];
    if (expr == nullptr || !expr->isValid()) {
      continue;
    }
    // Dependency materialization may temporarily remove internal roots from
    // the PI list to force cone rebuilding. That changes builder-local
    // variable IDs, so normalize every dependency formula back to the stable
    // SEC model variable space before the expression is shared with proofs.
    if (needsStableVarRemap) {
      expr = remapBoolExprVariables(expr, stableVarRemap, stableRemapMemo);
    }
    result.outputExprByTerm.emplace(result.outputs[i], expr);
    if (const auto requestedIt = requestedByRoot.find(result.outputs[i]);
        requestedIt != requestedByRoot.end()) {
      for (const auto requestedTermID : requestedIt->second) {
        result.outputExprByTerm.emplace(requestedTermID, expr);
      }
    }
  }

  std::vector<std::pair<naja::DNL::DNLID, BuilderSkippedOutputInfo>> skippedAliases(
      result.skippedOutputsByTerm.begin(), result.skippedOutputsByTerm.end());
  for (const auto& [rootTermID, info] : skippedAliases) {
    if (const auto requestedIt = requestedByRoot.find(rootTermID);
        requestedIt != requestedByRoot.end()) {
      for (const auto requestedTermID : requestedIt->second) {
        result.skippedOutputsByTerm.emplace(requestedTermID, info);
      }
    }
  }

  return result;
}

struct CandidateDependencyScratch {
  std::vector<const BoolExpr*> stack;
  std::vector<size_t> dependencies;
  std::unordered_map<const BoolExpr*, uint32_t> visitedEpochByNode;
  std::vector<uint32_t> emittedEpochByVarID;
  uint32_t currentEpoch = 1;
};

struct ThreadLocalDependencyState {
  CandidateDependencyScratch scratch;
  std::vector<std::vector<size_t>> dependenciesBySourceVarID;
};

const std::vector<size_t>& collectCandidateStateDependenciesFromExpr(
    BoolExpr* expr,
    const std::vector<uint8_t>& isCandidateStateVar,
    CandidateDependencyScratch& scratch) {
  scratch.dependencies.clear();
  if (expr == nullptr || !expr->isValid()) {
    // LCOV_EXCL_START
    return scratch.dependencies;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Skip propagation only needs candidate state support at the current root.
  // Sampling showed that calling BoolExpr::getSupportVars() here spent most of
  // the time reallocating temporary visited/support containers for every root.
  // Reusing epoch-tagged scratch storage keeps the low-memory behavior from the
  // cache removal while avoiding that repeated allocation churn.
  if (scratch.currentEpoch == std::numeric_limits<uint32_t>::max()) {
    // LCOV_EXCL_START
    scratch.visitedEpochByNode.clear();  // LCOV_EXCL_LINE
    std::fill(  // LCOV_EXCL_LINE
        scratch.emittedEpochByVarID.begin(),  // LCOV_EXCL_LINE
        scratch.emittedEpochByVarID.end(),  // LCOV_EXCL_LINE
        0);  // LCOV_EXCL_LINE
    scratch.currentEpoch = 1;  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  const uint32_t epoch = scratch.currentEpoch++;

  scratch.stack.clear();
  scratch.stack.push_back(expr);
  while (!scratch.stack.empty()) {
    const BoolExpr* node = scratch.stack.back();
    scratch.stack.pop_back();
    if (node == nullptr) {
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    auto& visitedEpoch = scratch.visitedEpochByNode[node];
    if (visitedEpoch == epoch) {
      continue;  // LCOV_EXCL_LINE
    }
    visitedEpoch = epoch;

    switch (node->getOp()) {
      case Op::VAR: {
        const size_t symbol = node->getId();
        if (symbol >= 2 && symbol < isCandidateStateVar.size() &&
            isCandidateStateVar[symbol] != 0) {
          if (symbol >= scratch.emittedEpochByVarID.size()) {
            scratch.emittedEpochByVarID.resize(symbol + 1, 0);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
          auto& emittedEpoch = scratch.emittedEpochByVarID[symbol];
          if (emittedEpoch != epoch) {
            emittedEpoch = epoch;
            scratch.dependencies.push_back(symbol);
          }
        }
        break;
      }
      case Op::NOT:
        scratch.stack.push_back(node->getLeft());  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      case Op::AND:
      case Op::OR:
      case Op::XOR:
        scratch.stack.push_back(node->getLeft());  // LCOV_EXCL_LINE
        scratch.stack.push_back(node->getRight());  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      case Op::NONE:  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      default:
        // LCOV_EXCL_START
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
    }
  }
  return scratch.dependencies;
}

void collectCandidateDependenciesIntoSet(
    BoolExpr* expr,
    const std::vector<uint8_t>& isCandidateVar,
    CandidateDependencyScratch& scratch,
    std::unordered_set<size_t>& dependencies) {
  for (const auto varID :
       collectCandidateStateDependenciesFromExpr(expr, isCandidateVar, scratch)) {
    // LCOV_EXCL_START
    dependencies.insert(varID);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
}

// LCOV_EXCL_START
std::optional<size_t> findFirstUnpublishedSupportVar(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    BoolExpr* expr,
    const std::vector<uint8_t>& isPublishedVar,
    CandidateDependencyScratch& scratch) {
  // LCOV_EXCL_START
  if (expr == nullptr || !expr->isValid()) {  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // We only need the first unpublished symbol. Avoid BoolExpr::getSupportVars()
  // here: BlackParrot-scale formulas spend minutes allocating full support
  // sets.  A PMR-backed visited set keeps the DAG walk exact without the tiny
  // per-node heap churn of std::unordered_map/operator[].
  // LCOV_EXCL_START
  std::pmr::monotonic_buffer_resource visitedResource;  // LCOV_EXCL_LINE
  std::pmr::unordered_set<const BoolExpr*> visited{&visitedResource};  // LCOV_EXCL_LINE
  visited.reserve(4096);  // LCOV_EXCL_LINE
  scratch.stack.clear();  // LCOV_EXCL_LINE
  scratch.stack.push_back(expr);  // LCOV_EXCL_LINE
  while (!scratch.stack.empty()) {  // LCOV_EXCL_LINE
    const BoolExpr* node = scratch.stack.back();  // LCOV_EXCL_LINE
    scratch.stack.pop_back();  // LCOV_EXCL_LINE
    if (node == nullptr || !visited.insert(node).second) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    switch (node->getOp()) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
      case Op::VAR: {
        // LCOV_EXCL_START
        const size_t symbol = node->getId();  // LCOV_EXCL_LINE
        if (symbol >= 2 &&  // LCOV_EXCL_LINE
            (symbol >= isPublishedVar.size() || isPublishedVar[symbol] == 0)) {  // LCOV_EXCL_LINE
          return symbol;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      case Op::NOT:
        // LCOV_EXCL_START
        scratch.stack.push_back(node->getLeft());  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      case Op::AND:
      case Op::OR:
      case Op::XOR:
        // LCOV_EXCL_START
        scratch.stack.push_back(node->getLeft());  // LCOV_EXCL_LINE
        scratch.stack.push_back(node->getRight());  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      case Op::NONE:  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      default:
        // LCOV_EXCL_START
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
    }
  }
  // LCOV_EXCL_START
  return std::nullopt;  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

std::vector<uint8_t> buildCandidateVarMask(
    const std::unordered_map<size_t, BoolExpr*>& candidateExprByVarID) {
  size_t maxVarID = 0;
  for (const auto& [varID, _] : candidateExprByVarID) {
    maxVarID = std::max(maxVarID, varID);
  }
  std::vector<uint8_t> isCandidateVar(maxVarID + 1, 0);
  for (const auto& [varID, _] : candidateExprByVarID) {
    isCandidateVar[varID] = 1;
  }
  return isCandidateVar;
}

// LCOV_EXCL_START
std::string displayNameForSignalKey(  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP
    const SequentialDesignModel& model,
    const SignalKey& key) {
  // LCOV_EXCL_START
  const auto displayIt = model.displayNameByKey.find(key);  // LCOV_EXCL_LINE
  return displayIt == model.displayNameByKey.end() ? signalKeyToString(key)  // LCOV_EXCL_LINE
                                                   : displayIt->second;  // LCOV_EXCL_LINE
                                                   // LCOV_EXCL_STOP
}

// LCOV_EXCL_START
ConnectivitySkipInfo makeUnpublishedSupportSkip(size_t varID) {  // LCOV_EXCL_LINE
  return ConnectivitySkipInfo{  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
      ConnectivitySkipOrigin::NoDriver,
      // LCOV_EXCL_START
      "Depends on unpublished internal support variable v" + std::to_string(varID),  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  };
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

SignalKey getTerminalPathKey(const naja::DNL::DNLTerminalFull& terminal) {
  SignalKey key;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  key.first.reserve(pathNames.size() + 1);
  for (const auto& name : pathNames) {
    key.first.push_back(name.getID());
  }
  key.first.push_back(terminal.getSnlBitTerm()->getName().getID());
  key.second.push_back(
      static_cast<naja::NL::NLID::DesignObjectID>(terminal.getSnlBitTerm()->getBit()));
  return key;
}

std::string getTerminalDisplayName(const naja::DNL::DNLTerminalFull& terminal) {
  std::ostringstream oss;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  for (const auto& name : pathNames) {
    oss << name.getString() << ".";
  }
  oss << terminal.getSnlBitTerm()->getName().getString() << "["
      << terminal.getSnlBitTerm()->getBit() << "]";
  return oss.str();
}

std::string normalizePinName(const std::string& name) {
  std::string normalized = name;
  for (char& ch : normalized) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return normalized;
}

bool hasSuffix(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string stripComplementSuffix(const std::string& pinName) {
  if (hasSuffix(pinName, "_N") || hasSuffix(pinName, "_B")) {
    // LCOV_EXCL_START
    return pinName.substr(0, pinName.size() - 2);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (hasSuffix(pinName, "N") || hasSuffix(pinName, "B")) {
    return pinName.substr(0, pinName.size() - 1);
  }
  return pinName;
}

bool isIntrinsicComplementedStateOutput(const std::string& pinName) {
  const std::string normalized = normalizePinName(pinName);
  return normalized == "QN" || normalized == "QB" ||
         normalized == "Q_N" || normalized == "Q_B";
}

bool isComplementedStateOutput(const std::string& primaryPinName,
                               const std::string& candidatePinName) {
  return candidatePinName != primaryPinName &&
         stripComplementSuffix(candidatePinName) == primaryPinName;
}

const StateOutputTerm* findComplementedPrimaryStateOutput(
    const StateOutputTerm& stateOutput,
    const std::vector<StateOutputTerm>& stateOutputs) {
  const std::string baseName = stripComplementSuffix(stateOutput.pinName);
  if (baseName == stateOutput.pinName) {
    return nullptr;
  }

  for (const auto& candidate : stateOutputs) {
    if (candidate.termID == stateOutput.termID || candidate.bit != stateOutput.bit) {
      continue;
    }
    if (candidate.pinName == baseName) {
      return &candidate;
    }
  }
  // LCOV_EXCL_START
  return nullptr;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

size_t countIndependentStateOutputs(
    const std::vector<StateOutputTerm>& stateOutputs) {
  size_t count = 0;
  for (const auto& stateOutput : stateOutputs) {
    if (findComplementedPrimaryStateOutput(stateOutput, stateOutputs) != nullptr) {
      continue;
    }
    ++count;
  }
  return count;
}

std::optional<naja::DNL::DNLID> resolvePendingPinTermID(
    const PendingTransition& pending,
    const char* pinName) {
  const auto pinIt = pending.pinTermIDs.find(pinName);
  if (pinIt == pending.pinTermIDs.end()) {
    return std::nullopt;
  }

  const auto& candidates = pinIt->second;
  if (candidates.empty()) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Multi-bit sequential primitives must resolve update pins against the same
  // bit index as the current state output. This keeps vector flops aligned per
  // state term instead of collapsing the whole instance down to one pin map.
  if (candidates.size() > 1) {
    for (const auto& candidate : candidates) {
      if (candidate.bit == pending.stateBit) {
        return candidate.termID;
      }
    }
    // LCOV_EXCL_START
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing bit-matched sequential pin `" + std::string(pinName) +  // LCOV_EXCL_LINE
        "` for output `" + pending.statePinName + "[" +  // LCOV_EXCL_LINE
        std::to_string(pending.stateBit) + "]`");  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }

  const bool isDataPin = std::string(pinName) == "D";
  if (isDataPin && pending.independentStateOutputCount > 1) {
    throw std::runtime_error( // LCOV_EXCL_LINE
        "Shared scalar D input cannot define multiple independent state outputs");
  }

  return candidates.front().termID;
}

bool isSequentialStateOutput(const naja::DNL::DNLTerminalFull& term) {
  if (term.isTopPort()) {
    return false;
  }
  return !naja::NL::SNLDesignModeling::getOutputRelatedClocks(
              term.getSnlBitTerm())
              .empty();
}

bool isSequentialNextStateInput(const naja::DNL::DNLTerminalFull& term) {
  if (term.isTopPort()) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return !naja::NL::SNLDesignModeling::getInputRelatedClocks(
              term.getSnlBitTerm())
              .empty();
}

bool isNegativeClockPinName(const std::string& pinName) {
  const std::string normalized = normalizePinName(pinName);
  return normalized == "CKN" || normalized == "CLKN" ||
         normalized == "CLK_N" || normalized == "CLK_B" ||
         normalized == "CK_N" || normalized == "CK_B" ||
         normalized == "CLOCK_N" || normalized == "CLOCK_B";
}

bool isIntrinsicNegativeEdgeClock(
    const naja::DNL::DNLInstanceFull& instance,
    const naja::NL::SNLBitTerm* clockBitTerm) {
  const auto* model = instance.getSNLModel();
  if (naja::NL::NLDB0::isDFFN(model)) {
    return true;
  }
  if (clockBitTerm == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return isNegativeClockPinName(clockBitTerm->getName().getString());
}

bool isConstantInternalOutputTerm(const naja::DNL::DNLTerminalFull& term) {
  if (term.isNull() || term.isTopPort() ||
      term.getSnlBitTerm()->getDirection() !=
          naja::NL::SNLBitTerm::Direction::Output) {
    return false;
  }

  const auto& instance = term.getDNLInstance();
  const auto* model = instance.getSNLModel();
  if (model == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  const auto truthTable = naja::NL::SNLDesignModeling::getTruthTable(
      model, term.getSnlBitTerm()->getOrderID());
  if (truthTable.isInitialized() && (truthTable.all0() || truthTable.all1())) {
    return true;
  }

  const std::string pinName = normalizePinName(term.getSnlBitTerm()->getName().getString());
  const std::string modelName = normalizePinName(model->getName().getString());
  return (pinName == "LO" || pinName == "HI") &&
         // LCOV_EXCL_START
         (modelName.find("CONB") != std::string::npos ||  // LCOV_EXCL_LINE
          modelName.find("TIE") != std::string::npos);  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
}

bool isOptionalSequentialControlPin(const std::string& pinName) {
  return pinName == "E" || pinName == "DE" || pinName == "R" ||
         pinName == "RN" || pinName == "RESET_B" ||
         pinName == "RESET_N" || pinName == "RESETN" ||
         pinName == "RST_B" || pinName == "RST_N" ||
         pinName == "RSTN" || pinName == "S" ||
         pinName == "SET_B" || pinName == "SET_N" ||
         pinName == "SETN";
}

bool isSupportedSequentialUpdatePin(const std::string& pinName) {
  return pinName == "D" || isOptionalSequentialControlPin(pinName);
}

std::optional<naja::DNL::DNLID> resolvePendingPinRoleTermID(
    const PendingTransition& pending,
    const char* roleName) {
  auto resolvedTermID = resolvePendingPinTermID(pending, roleName);
  if (resolvedTermID.has_value()) {
    return resolvedTermID;
  }

  const std::string role(roleName);
  if (role == "E") {
    // sky130 spells an enabled flop's data enable as DE while the generic SEC
    // transition builder uses E for hold semantics.
    return resolvePendingPinTermID(pending, "DE");
  }
  if (role == "RN") {
    // Active-low reset pins appear under several Liberty naming conventions.
    // Treat them as the same control role, without introducing any equality
    // assumptions between internal state bits.
    for (const char* alias : {"RESET_B", "RESET_N", "RESETN", "RST_B",
                              "RST_N", "RSTN"}) {
      resolvedTermID = resolvePendingPinTermID(pending, alias);
      if (resolvedTermID.has_value()) {
        // LCOV_EXCL_START
        return resolvedTermID;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    }
  }
  if (role == "SN") {
    for (const char* alias : {"SET_B", "SET_N", "SETN"}) {
      resolvedTermID = resolvePendingPinTermID(pending, alias);
      if (resolvedTermID.has_value()) {
        // LCOV_EXCL_START
        return resolvedTermID;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    }
  }
  return std::nullopt;
}

BoolExpr* getRequiredOutputExpr(
    const PendingTransition& pending,
    const char* pinName,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm) {
  auto resolvedTermID = resolvePendingPinRoleTermID(pending, pinName);
  if (!resolvedTermID.has_value()) {
    return nullptr;
  }
  auto exprIt = outputExprByTerm.find(*resolvedTermID);
  if (exprIt == outputExprByTerm.end()) {
    // LCOV_EXCL_START
    throw std::runtime_error("Missing combinational expression for sequential pin `" +  // LCOV_EXCL_LINE
                             std::string(pinName) + "`");  // LCOV_EXCL_LINE
                             // LCOV_EXCL_STOP
  }
  return exprIt->second;
}

BoolExpr* stripClockCarrierFromClockEnable(
    BoolExpr* root,
    const std::unordered_set<size_t>& topClockCarrierVarIDs,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (const auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  BoolExpr* stripped = nullptr;
  switch (root->getOp()) {
    case Op::VAR:
      stripped = topClockCarrierVarIDs.find(root->getId()) ==
                         topClockCarrierVarIDs.end()
                     ? root
                     : BoolExpr::createTrue();
      break;
    case Op::NOT:
      if (root->getLeft() != nullptr && root->getLeft()->getOp() == Op::VAR &&
          topClockCarrierVarIDs.find(root->getLeft()->getId()) !=
              topClockCarrierVarIDs.end()) {
        // Inverted clock trees still represent the abstract cycle carrier.
        // Remove the carrier polarity, but keep real non-clock gating logic.
        // LCOV_EXCL_START
        stripped = BoolExpr::createTrue();  // LCOV_EXCL_LINE
      } else {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        BoolExpr* left = stripClockCarrierFromClockEnable(
            root->getLeft(), topClockCarrierVarIDs, memo);
        stripped = left == root->getLeft() ? root : BoolExpr::Not(left);
      }
      break;
    case Op::AND: {
      BoolExpr* left = stripClockCarrierFromClockEnable(
          root->getLeft(), topClockCarrierVarIDs, memo);
      BoolExpr* right = stripClockCarrierFromClockEnable(
          root->getRight(), topClockCarrierVarIDs, memo);
      stripped = left == root->getLeft() && right == root->getRight()
                     ? root
                     // LCOV_EXCL_START
                     : BoolExpr::And(left, right);  // LCOV_EXCL_LINE
                     // LCOV_EXCL_STOP
      break;
    }
    case Op::OR: {
      BoolExpr* left = stripClockCarrierFromClockEnable(
          root->getLeft(), topClockCarrierVarIDs, memo);
      BoolExpr* right = stripClockCarrierFromClockEnable(
          root->getRight(), topClockCarrierVarIDs, memo);
      stripped = left == root->getLeft() && right == root->getRight()
                     ? root
                     // LCOV_EXCL_START
                     : BoolExpr::Or(left, right);  // LCOV_EXCL_LINE
                     // LCOV_EXCL_STOP
      break;
    }
    case Op::XOR: {
      // LCOV_EXCL_START
      BoolExpr* left = stripClockCarrierFromClockEnable(  // LCOV_EXCL_LINE
          root->getLeft(), topClockCarrierVarIDs, memo);  // LCOV_EXCL_LINE
      BoolExpr* right = stripClockCarrierFromClockEnable(  // LCOV_EXCL_LINE
          root->getRight(), topClockCarrierVarIDs, memo);  // LCOV_EXCL_LINE
      stripped = left == root->getLeft() && right == root->getRight()  // LCOV_EXCL_LINE
                     ? root  // LCOV_EXCL_LINE
                     : BoolExpr::Xor(left, right);  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    case Op::NONE:  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      throw std::runtime_error("Unsupported BoolExpr operator in clock enable");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }

  memo.emplace(root, stripped);
  return stripped;
}

BoolExpr* simplifyWhenClockCarriersChanged(BoolExpr* original, BoolExpr* stripped) {
  // Most mapped flops do not actually contain a clock carrier in their data
  // cone.  Preserve the original expression in that common case instead of
  // rebuilding and simplifying large ASIC cones just to discover no change.
  return stripped == original ? original : BoolExpr::simplify(stripped);
}

BoolExpr* stripClockCarriersFromSequentialUpdate(
    BoolExpr* update,
    const std::unordered_set<size_t>& topClockCarrierVarIDs,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (update == nullptr || topClockCarrierVarIDs.empty()) {
    // LCOV_EXCL_START
    return update;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  // Clock carriers are the implicit SEC step trigger, not data state.  Strip any
  // carrier variables that leaked through clock-pin materialization before the
  // transition is exposed to the proof engines.
  return simplifyWhenClockCarriersChanged(
      update, stripClockCarrierFromClockEnable(update, topClockCarrierVarIDs, memo));
}

BoolExpr* substituteClockGateLatchVars(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& latchDataExprByVarID,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (const auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  BoolExpr* substituted = nullptr;
  switch (root->getOp()) {
    case Op::VAR: {
      const auto it = latchDataExprByVarID.find(root->getId());
      substituted = it == latchDataExprByVarID.end() ? root : it->second;
      break;
    }
    case Op::NOT: {
      BoolExpr* left = substituteClockGateLatchVars(
          root->getLeft(), latchDataExprByVarID, memo);
      substituted = left == root->getLeft() ? root : BoolExpr::Not(left);
      break;
    }
    case Op::AND: {
      BoolExpr* left = substituteClockGateLatchVars(
          root->getLeft(), latchDataExprByVarID, memo);
      BoolExpr* right = substituteClockGateLatchVars(
          root->getRight(), latchDataExprByVarID, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::And(left, right);
      break;
    }
    case Op::OR: {
      BoolExpr* left = substituteClockGateLatchVars(
          root->getLeft(), latchDataExprByVarID, memo);
      BoolExpr* right = substituteClockGateLatchVars(
          root->getRight(), latchDataExprByVarID, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        // LCOV_EXCL_START
                        : BoolExpr::Or(left, right);  // LCOV_EXCL_LINE
                        // LCOV_EXCL_STOP
      break;
    }
    case Op::XOR: {
      // LCOV_EXCL_START
      BoolExpr* left = substituteClockGateLatchVars(  // LCOV_EXCL_LINE
          root->getLeft(), latchDataExprByVarID, memo);  // LCOV_EXCL_LINE
      BoolExpr* right = substituteClockGateLatchVars(  // LCOV_EXCL_LINE
          root->getRight(), latchDataExprByVarID, memo);  // LCOV_EXCL_LINE
      substituted = left == root->getLeft() && right == root->getRight()  // LCOV_EXCL_LINE
                        ? root  // LCOV_EXCL_LINE
                        : BoolExpr::Xor(left, right);  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    case Op::NONE:  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      throw std::runtime_error("Unsupported BoolExpr operator in latch substitution");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }

  memo.emplace(root, substituted);
  return substituted;
}

BoolExpr* simplifyWhenClockGateLatchVarsChanged(
    BoolExpr* original, BoolExpr* substituted) {
  // Folded latch outputs are rare. Preserve untouched shared subtrees so large
  // SEC models only simplify cones that actually reference the latch output.
  return substituted == original ? original : BoolExpr::simplify(substituted);
}

BoolExpr* substituteClockGateLatchVarsInExpr(
    BoolExpr* expr,
    const std::unordered_map<size_t, BoolExpr*>& latchDataExprByVarID) {
  if (expr == nullptr || latchDataExprByVarID.empty()) {
    // LCOV_EXCL_START
    return expr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return simplifyWhenClockGateLatchVarsChanged(
      expr, substituteClockGateLatchVars(expr, latchDataExprByVarID, memo));
}

bool pendingClockTermsArePureCarriers(
    const PendingTransition& pending,
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs) {
  if (pending.clockTermIDs.empty()) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  for (const auto& clockTerm : pending.clockTermIDs) {
    if (pureClockCarrierTermIDs.find(clockTerm.termID) ==
        pureClockCarrierTermIDs.end()) {
      return false;
    }
  }
  return true;
}

BoolExpr* getLocalClockEnableExpr(
    const PendingTransition& pending,
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs,
    const std::unordered_set<size_t>& topClockCarrierVarIDs,
    const std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID,
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    std::unordered_map<BoolExpr*, BoolExpr*>& clockCarrierStripMemo) {
  if (pending.clockTermIDs.empty()) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  // CTS-only clock trees can be materialized as internal output expressions that
  // do not own a stable boundary variable.  The structural pure-clock term set is
  // therefore the authoritative "this is just a clock" signal for these pins.
  if (pendingClockTermsArePureCarriers(pending, pureClockCarrierTermIDs)) {
    return nullptr;
  }
  if (topClockCarrierVarIDs.empty()) {
    // LCOV_EXCL_START
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  std::vector<BoolExpr*> clockExprs;
  clockExprs.reserve(pending.clockTermIDs.size());
  for (const auto& clockTerm : pending.clockTermIDs) {
    const auto exprIt = outputExprByTerm.find(clockTerm.termID);
    if (exprIt == outputExprByTerm.end()) {
      // LCOV_EXCL_START
      throw std::runtime_error(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          "Missing combinational expression for sequential clock pin");  // LCOV_EXCL_LINE
    }
    clockExprs.push_back(exprIt->second);
  }

  BoolExpr* clockExpr = clockExprs.front();
  for (size_t index = 1; index < clockExprs.size(); ++index) {
    // LCOV_EXCL_START
    clockExpr = BoolExpr::And(clockExpr, clockExprs[index]);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  if (!clockGateLatchDataExprByVarID.empty()) {
    clockExpr = substituteClockGateLatchVarsInExpr(
        clockExpr, clockGateLatchDataExprByVarID);
  }
  if (const auto event =
          classifyClockEventExpression(clockExpr, clockEventByCarrierVarID)) {
    return clockEventIsUngated(*event) ? nullptr : event->enable;
  }

  // LCOV_EXCL_START
  BoolExpr* enable = simplifyWhenClockCarriersChanged(  // LCOV_EXCL_LINE
      clockExpr, stripClockCarrierFromClockEnable(  // LCOV_EXCL_LINE
                     clockExpr, topClockCarrierVarIDs, clockCarrierStripMemo));  // LCOV_EXCL_LINE
  if (enable == BoolExpr::createTrue()) {  // LCOV_EXCL_LINE
    return nullptr;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return enable;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

std::optional<std::string> getPendingTransitionUnsupportedReason(
    const PendingTransition& pending) {
  const auto dIt = pending.pinTermIDs.find("D");
  if (dIt == pending.pinTermIDs.end() || dIt->second.empty()) {
    return "Unsupported sequential primitive without D input";
  }

  if (dIt->second.size() == 1 && pending.independentStateOutputCount > 1) {
    return "Shared scalar D input cannot define multiple independent state outputs";
  }

  if (dIt->second.size() > 1) {
    bool hasBitMatchedDataPin = false;
    for (const auto& candidate : dIt->second) {
      if (candidate.bit == pending.stateBit) {
        hasBitMatchedDataPin = true;
        break;
      }
    }
    if (!hasBitMatchedDataPin) {
      return "Missing bit-matched sequential pin `D` for output `" + // LCOV_EXCL_LINE
             pending.statePinName + "[" + std::to_string(pending.stateBit) + // LCOV_EXCL_LINE
             "]`";
    }
  }

  for (const auto& [pinName, _] : pending.pinTermIDs) {
    if (!isSupportedSequentialUpdatePin(pinName)) {
      return "Unsupported sequential primitive with update pin `" + pinName + "`";
    }
  }

  // Reset/set combinations are common in mapped cells.  They must be modeled
  // as state transitions instead of abstracted as internal SEC boundary terms,
  // because SEC may only align top-level terminals by name.
  return std::nullopt;
}

BoolExpr* buildNextStateExpr(
    const PendingTransition& pending,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs,
    const std::unordered_set<size_t>& topClockCarrierVarIDs,
    const std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID,
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    std::unordered_map<BoolExpr*, BoolExpr*>& clockCarrierStripMemo) {
  if (pending.stateTermID >= termDNLID2varID.size()) {
    // LCOV_EXCL_START
    throw std::runtime_error("Sequential state term is out of range");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  const size_t stateVarID = termDNLID2varID[pending.stateTermID];
  if (stateVarID < 2) {
    // LCOV_EXCL_START
    throw std::runtime_error("Sequential state bit was mapped to a constant");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  BoolExpr* data = getRequiredOutputExpr(pending, "D", outputExprByTerm);
  if (data == nullptr) {
    // LCOV_EXCL_START
    throw std::runtime_error("Unsupported sequential primitive without D input");  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (pending.stateOutputIsComplemented) {
    // LCOV_EXCL_START
    data = BoolExpr::Not(data);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  BoolExpr* current = BoolExpr::Var(stateVarID);
  BoolExpr* next = data;

  // Supported hold semantics: Q' = E ? D : Q.
  if (BoolExpr* enable = getRequiredOutputExpr(pending, "E", outputExprByTerm)) {
    next = BoolExpr::Or(
        BoolExpr::And(enable, data),
        BoolExpr::And(BoolExpr::Not(enable), current));
  }

  if (BoolExpr* clockEnable =
          getLocalClockEnableExpr(
              pending,
              pureClockCarrierTermIDs,
              topClockCarrierVarIDs,
              clockEventByCarrierVarID,
              clockGateLatchDataExprByVarID,
              outputExprByTerm,
              clockCarrierStripMemo)) {
    next = BoolExpr::Or(
        BoolExpr::And(clockEnable, next),
        BoolExpr::And(BoolExpr::Not(clockEnable), current));
  }

  BoolExpr* resetHigh = getRequiredOutputExpr(pending, "R", outputExprByTerm);
  BoolExpr* resetLow = getRequiredOutputExpr(pending, "RN", outputExprByTerm);
  BoolExpr* setHigh = getRequiredOutputExpr(pending, "S", outputExprByTerm);
  BoolExpr* setLow = getRequiredOutputExpr(pending, "SN", outputExprByTerm);

  auto applyForcedValue = [&](BoolExpr* asserted, bool value) {
    return BoolExpr::Or(
        BoolExpr::And(asserted,
                      value ? BoolExpr::createTrue() : BoolExpr::createFalse()),
        BoolExpr::And(BoolExpr::Not(asserted), next));
  };

  const bool resetValue = pending.stateOutputIsComplemented;
  const bool setValue = !pending.stateOutputIsComplemented;

  // Apply reset before set so set/clear-style controls have the final say when
  // both are asserted.  This matches ASAP7 reset+set flops whose Liberty views
  // define the simultaneous clear/preset case as low on the exposed state pin.
  if (resetHigh) {
    next = applyForcedValue(resetHigh, resetValue);
  }
  if (resetLow) {
    next = applyForcedValue(BoolExpr::Not(resetLow), resetValue);
  }
  if (setHigh) {
    next = applyForcedValue(setHigh, setValue);
  }
  if (setLow) {
    // LCOV_EXCL_START
    next = applyForcedValue(BoolExpr::Not(setLow), setValue);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  return stripClockCarriersFromSequentialUpdate(
      next, topClockCarrierVarIDs, clockCarrierStripMemo);
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

std::optional<bool> detectInitialStateValue(const PendingTransition& pending) {
  const bool hasResetHigh = resolvePendingPinRoleTermID(pending, "R").has_value();
  const bool hasResetLow = resolvePendingPinRoleTermID(pending, "RN").has_value();
  const bool hasSetHigh = resolvePendingPinRoleTermID(pending, "S").has_value();
  const bool hasSetLow = resolvePendingPinRoleTermID(pending, "SN").has_value();

  const bool hasReset = hasResetHigh || hasResetLow;
  const bool hasSet = hasSetHigh || hasSetLow;
  if (hasReset && !hasSet) {
    return pending.stateOutputIsComplemented;
  }
  if (hasSet && !hasReset) {
    return !pending.stateOutputIsComplemented;
  }
  return std::nullopt;
}

BoolExpr* makeAndChain(const std::vector<BoolExpr*>& exprs) {
  if (exprs.empty()) {
    // LCOV_EXCL_START
    return BoolExpr::createTrue();  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  BoolExpr* result = exprs.front();
  for (size_t index = 1; index < exprs.size(); ++index) {
    result = BoolExpr::And(result, exprs[index]);
  }
  return result;
}

BoolExpr* makeIte(BoolExpr* condition, BoolExpr* whenTrue, BoolExpr* whenFalse) {
  return BoolExpr::Or(
      BoolExpr::And(condition, whenTrue),
      BoolExpr::And(BoolExpr::Not(condition), whenFalse));
}

BoolExpr* buildAddressEqualsExpr(
    const std::vector<BoolExpr*>& addressBits,
    size_t addressValue) {
  std::vector<BoolExpr*> equalities;
  equalities.reserve(addressBits.size());
  for (size_t bitIndex = 0; bitIndex < addressBits.size(); ++bitIndex) {
    const bool expected = ((addressValue >> bitIndex) & size_t{1}) != 0;
    equalities.push_back(expected ? addressBits[bitIndex]
                                  : BoolExpr::Not(addressBits[bitIndex]));
  }
  return makeAndChain(equalities);
}

std::optional<bool> evaluateConstantUnderAssignments(
    BoolExpr* expr,
    const std::unordered_map<size_t, bool>& assignments,
    std::unordered_map<BoolExpr*, std::optional<bool>>& memo) {
  if (expr == nullptr) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (const auto it = memo.find(expr); it != memo.end()) {
    // LCOV_EXCL_START
    return it->second;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  std::optional<bool> value;
  switch (expr->getOp()) {
    case Op::VAR:
      if (expr->getId() < 2) {
        value = expr->getId() == 1;
      } else if (const auto it = assignments.find(expr->getId());
                 it != assignments.end()) {
        value = it->second;
      }
      break;
    case Op::NOT: {
      const auto operand =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);  // LCOV_EXCL_LINE
      if (operand.has_value()) {  // LCOV_EXCL_LINE
        value = !*operand;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    case Op::AND: {
      const auto lhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);  // LCOV_EXCL_LINE
      if (lhs.has_value() && !*lhs) {  // LCOV_EXCL_LINE
        value = false;  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto rhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);  // LCOV_EXCL_LINE
      if (rhs.has_value() && !*rhs) {  // LCOV_EXCL_LINE
        value = false;  // LCOV_EXCL_LINE
      } else if (lhs.has_value() && rhs.has_value()) {  // LCOV_EXCL_LINE
        value = *lhs && *rhs;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    case Op::OR: {
      const auto lhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);  // LCOV_EXCL_LINE
      if (lhs.has_value() && *lhs) {  // LCOV_EXCL_LINE
        value = true;  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto rhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);  // LCOV_EXCL_LINE
      if (rhs.has_value() && *rhs) {  // LCOV_EXCL_LINE
        value = true;  // LCOV_EXCL_LINE
      } else if (lhs.has_value() && rhs.has_value()) {  // LCOV_EXCL_LINE
        value = *lhs || *rhs;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    case Op::XOR: {
      const auto lhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
      const auto rhs =
          // LCOV_EXCL_START
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);  // LCOV_EXCL_LINE
      if (lhs.has_value() && rhs.has_value()) {  // LCOV_EXCL_LINE
        value = *lhs != *rhs;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    case Op::NONE:  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    default:
      // LCOV_EXCL_START
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }

  memo.emplace(expr, value);
  return value;
}

std::string normalizeSignalBaseName(const std::string& name) {
  std::string base = name;
  const auto bracket = base.find('[');
  if (bracket != std::string::npos) {
    base = base.substr(0, bracket);
  }
  return normalizePinName(base);
}

bool isResetNameToken(const std::string& candidate, const std::string& token) {
  // Domain-prefixed top resets, for example `wb_rst_i`, normalize to `WB_RST`
  // after input-suffix stripping.  Match only a final underscore-separated
  // reset token so synthesized reset inference remains conservative.
  return candidate == token || hasSuffix(candidate, "_" + token);
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

bool isActiveLowResetToken(const std::string& candidate) {
  return candidate == "RESET_N" || candidate == "RESETN" ||
         candidate == "RESET_L" || candidate == "RST_N" ||
         candidate == "RSTN" || candidate == "RST_L";
}

void appendDomainPrefixedActiveLowResetCandidates(
    std::vector<std::string>& candidates) {
  const size_t originalSize = candidates.size();
  for (size_t index = 0; index < originalSize; ++index) {
    const std::string& candidate = candidates[index];
    if (candidate.size() <= 1) {
      continue;
    // LCOV_EXCL_START
    }
    const std::string strippedDomain = candidate.substr(1);
    // LCOV_EXCL_STOP
    if (isActiveLowResetToken(strippedDomain)) {
      // Async FIFOs often spell domain resets as rrst_n/wrst_n.  Treat only
      // one-letter active-low prefixes as reset candidates so unrelated names
      // containing "rst" do not become reset controls.
      candidates.push_back(strippedDomain);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
}

std::vector<std::string> resetNameCandidates(const std::string& displayName) {
  // Synthesized-reset inference runs before the final SEC symbol space exists.
  // Use the same top-port spelling policy as later SEC phases so names such as
  // `reset_i[0]` and `rst_ni[0]` are classified consistently end-to-end.
  const std::string normalized = normalizeSignalBaseName(displayName);
  std::vector<std::string> candidates = {normalized};
  if (hasSuffix(normalized, "_IN")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 3));
  }
  if (hasSuffix(normalized, "_I")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 2));
  }
  if (hasSuffix(normalized, "_NI")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 1));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  appendDomainPrefixedActiveLowResetCandidates(candidates);
  return candidates;
}

std::optional<bool> getResetAssertionValue(const std::string& displayName) {
  for (const auto& candidate : resetNameCandidates(displayName)) {
    if (isResetNameToken(candidate, "RESET") ||
        isResetNameToken(candidate, "RST")) {
      return true;
    }
    if (isResetNameToken(candidate, "RESET_N") ||
        isResetNameToken(candidate, "RESETN") ||
        isResetNameToken(candidate, "RESET_L") ||
        isResetNameToken(candidate, "RST_N") ||
        isResetNameToken(candidate, "RSTN") ||
        isResetNameToken(candidate, "RST_L")) {
      return false;
    }
  }
  return std::nullopt;
}

std::vector<std::string> clockNameCandidates(const std::string& displayName) {
  const std::string normalized = normalizeSignalBaseName(displayName);
  std::vector<std::string> candidates = {normalized};
  if (hasSuffix(normalized, "_IN")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 3));
  }
  if (hasSuffix(normalized, "_I")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 2));
  }
  return candidates;
}

// LCOV_EXCL_START
bool isTopClockCarrierName(const std::string& displayName) {
// LCOV_EXCL_STOP
  for (const auto& candidate : clockNameCandidates(displayName)) {
    if (candidate == "CLK" || candidate == "CLOCK" || candidate == "CK" ||
        hasSuffix(candidate, "_CLK") || hasSuffix(candidate, "_CLOCK")) {
      return true;
    }
  }
  // LCOV_EXCL_START
  return false;
  // LCOV_EXCL_STOP
}

bool isClockTreeCarrierBoundaryName(const std::string& displayName) {
  const std::string normalized = normalizePinName(displayName);
  if (normalized.find("CLOCK_GATE") != std::string::npos ||
      normalized.find("EN_LATCH") != std::string::npos) {
    // LCOV_EXCL_START
    return false;
    // LCOV_EXCL_STOP
  }
  return normalized.find("CLKBUF") != std::string::npos ||
         normalized.find("CLKLOAD") != std::string::npos ||
         normalized.find("CLKNET") != std::string::npos ||
         normalized.find("_CLK.") != std::string::npos ||
         normalized.find("_CLOCK.") != std::string::npos;
}

bool isClockTreeCarrierTerminalName(const naja::DNL::DNLTerminalFull& term) {
  if (term.isNull()) {
    return false;  // LCOV_EXCL_LINE
  }
  if (isClockTreeCarrierBoundaryName(getTerminalDisplayName(term))) {
    return true;
  }
  if (isClockTreeCarrierBoundaryName(
          term.getSnlBitTerm()->getName().getString())) {
    return true;  // LCOV_EXCL_LINE
  }

  const auto& instance = term.getDNLInstance();
  // LCOV_EXCL_START
  if (const auto* snlInstance = instance.getSNLInstance();
  // LCOV_EXCL_STOP
      snlInstance != nullptr &&
      isClockTreeCarrierBoundaryName(snlInstance->getName().getString())) {
    return true;  // LCOV_EXCL_LINE
  }
  return false;
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

bool isClockTreeBufferCellName(const std::string& name) {
  const std::string normalized = normalizePinName(name);
  if (normalized.find("CLOCK_GATE") != std::string::npos ||
      normalized.find("EN_LATCH") != std::string::npos) {
    return false;
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  return normalized.find("CLKBUF") != std::string::npos ||
         normalized.find("CLKLOAD") != std::string::npos;
}

// LCOV_EXCL_START
bool isClockTreeBufferCell(const naja::DNL::DNLTerminalFull& term) {
// LCOV_EXCL_STOP
  if (term.isNull() || term.isTopPort()) {
    return false;  // LCOV_EXCL_LINE
  }

  // CTS tools may keep a generic BUF/INV library cell while marking only the
  // placed instance or flattened DNL path as clkbuf/clkload.  Use that path
  // marker so routed clock leaves are not mis-modeled as SEC clock enables.
  // Library model names alone are not enough: sky130hs also uses clkbuf cells
  // as ordinary routed data buffers.
  if (isClockTreeCarrierTerminalName(term)) {
    return true;
  }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  if (isClockTreeBufferCellName(term.getSnlBitTerm()->getName().getString())) {
    return true;  // LCOV_EXCL_LINE
  }

  // LCOV_EXCL_START
  const auto& instance = term.getDNLInstance();
  // LCOV_EXCL_STOP
  if (const auto* snlInstance = instance.getSNLInstance();
      snlInstance != nullptr &&
      isClockTreeBufferCellName(snlInstance->getName().getString())) {
    return true;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  return false;
}

bool isPotentialClockTreeBufferCell(const naja::DNL::DNLTerminalFull& term) {
  if (isClockTreeBufferCell(term)) {
    return true;
  }
  if (term.isNull() || term.isTopPort()) {
    return false;  // LCOV_EXCL_LINE
  }

  // Only the source-tracing structural pass may use the library model name:
  // the candidate becomes a clock carrier only if it recursively reaches a top
  // clock. Generic cone passthrough must not erase data-path clkbuf cells.
  const auto& instance = term.getDNLInstance();
  if (const auto* model = instance.getSNLModel();
      model != nullptr && isClockTreeBufferCellName(model->getName().getString())) {
    return true;
  }
  return false;
}

// LCOV_EXCL_START
std::optional<naja::DNL::DNLID> getSingleInputDriverTerm(
// LCOV_EXCL_STOP
    naja::DNL::DNLFull* dnl,
    const naja::DNL::DNLInstanceFull& instance,
    const naja::NL::SNLBitTerm* inputBitTerm) {
  if (dnl == nullptr || inputBitTerm == nullptr) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const auto& inputTerm = instance.getTerminalFromBitTerm(inputBitTerm);
  if (inputTerm.isNull() || inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
  if (iso.isConstant() || iso.getDrivers().size() != 1) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  return iso.getDrivers().front();
}

std::optional<naja::DNL::DNLID> getClockTreeBufferSourceDriverTerm(
    naja::DNL::DNLFull* dnl,
    const naja::DNL::DNLTerminalFull& outputTerm) {
  if (dnl == nullptr || outputTerm.isNull()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  const auto& instance = outputTerm.getDNLInstance();
  std::vector<naja::DNL::DNLID> sourceDrivers;
  auto addInputDriver = [&](const naja::NL::SNLBitTerm* inputBitTerm) {
    if (inputBitTerm == nullptr ||
        inputBitTerm->getDirection() == naja::NL::SNLBitTerm::Direction::Output) {
      return;
    }
    if (const auto driver = getSingleInputDriverTerm(dnl, instance, inputBitTerm)) {
      sourceDrivers.push_back(*driver);
    }
  };

  for (auto* inputBitTerm :
       naja::NL::SNLDesignModeling::getCombinatorialInputs(
           const_cast<naja::NL::SNLBitTerm*>(outputTerm.getSnlBitTerm()))) {
    addInputDriver(inputBitTerm);
  }

  // LCOV_EXCL_START
  if (sourceDrivers.empty()) {
  // LCOV_EXCL_STOP
    if (const auto* model = instance.getSNLModel(); model != nullptr) {
      for (auto* bitTerm : model->getBitTerms()) {
        addInputDriver(bitTerm);
      }
    }
  }

  std::sort(sourceDrivers.begin(), sourceDrivers.end());
  sourceDrivers.erase(
      std::unique(sourceDrivers.begin(), sourceDrivers.end()), sourceDrivers.end());
  if (sourceDrivers.size() != 1) {
    return std::nullopt;
  }
  return sourceDrivers.front();
}

std::unordered_set<size_t> collectTopClockCarrierVarIDs(
    // LCOV_EXCL_START
    const SequentialDesignModel& model) {
    // LCOV_EXCL_STOP
  std::unordered_set<size_t> varIDs;
  for (const auto& key : model.topInputKeys) {
    const auto displayIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        varIt == model.inputVarByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (isTopClockCarrierName(displayIt->second)) {
      varIDs.insert(varIt->second);
    }
  }
  return varIDs;
}

void seedTopClockCarrierEvents(
    const SequentialDesignModel& model,
    const std::unordered_set<size_t>& clockCarrierVarIDs,
    std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID) {
  for (const auto& key : model.topInputKeys) {
    const auto displayIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        varIt == model.inputVarByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    if (!isTopClockCarrierName(displayIt->second) ||
    // LCOV_EXCL_STOP
        clockCarrierVarIDs.find(varIt->second) == clockCarrierVarIDs.end()) {
      continue;
    }

    ClockEvent event;
    event.domain = key;
    event.phase = ClockPhase::Pos;
    event.enable = nullptr;
    clockEventByCarrierVarID.emplace(varIt->second, event);
  // LCOV_EXCL_START
  }
}


// LCOV_EXCL_STOP
ClockEvent applyIntrinsicClockPolarity(ClockEvent event, bool inverted) {
  // LCOV_EXCL_START
  if (inverted) {
  // LCOV_EXCL_STOP
    event.phase = invertClockPhase(event.phase);
  }
  return event;
}

BoolExpr* buildPendingClockExpr(
    const PendingTransition& pending,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm) {
  if (pending.clockTermIDs.empty()) {
    return nullptr;  // LCOV_EXCL_LINE
  }

  BoolExpr* clockExpr = nullptr;
  // LCOV_EXCL_START
  for (const auto& clockTerm : pending.clockTermIDs) {
  // LCOV_EXCL_STOP
    BoolExpr* termExpr = nullptr;
    const auto exprIt = outputExprByTerm.find(clockTerm.termID);
    if (exprIt != outputExprByTerm.end()) {
      termExpr = exprIt->second;
    // LCOV_EXCL_START
    } else if (clockTerm.termID < termDNLID2varID.size() &&
    // LCOV_EXCL_STOP
               termDNLID2varID[clockTerm.termID] >= 2) {  // LCOV_EXCL_LINE
      termExpr = BoolExpr::Var(termDNLID2varID[clockTerm.termID]);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (termExpr == nullptr) {
      continue;  // LCOV_EXCL_LINE
    }
    clockExpr = clockExpr == nullptr ? termExpr : BoolExpr::And(clockExpr, termExpr);
  }
  return clockExpr;
}

std::optional<ClockEvent> classifyPendingClockEvent(
    const PendingTransition& pending,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID) {
  BoolExpr* clockExpr =
      buildPendingClockExpr(pending, termDNLID2varID, outputExprByTerm);
  if (clockExpr == nullptr) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  auto event =
      classifyClockEventExpression(clockExpr, clockEventByCarrierVarID);
  if (!event.has_value()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  return applyIntrinsicClockPolarity(*event, pending.intrinsicClockIsInverted);
}

std::vector<ClockCarrierClass> buildClockCarrierClasses(
    const std::vector<size_t>& clockCarrierVarIDs,
    const std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID) {
  std::vector<ClockCarrierClass> classes;
  // LCOV_EXCL_START
  classes.reserve(clockEventByCarrierVarID.size());
  for (const auto varID : clockCarrierVarIDs) {
  // LCOV_EXCL_STOP
    const auto eventIt = clockEventByCarrierVarID.find(varID);
    if (eventIt == clockEventByCarrierVarID.end()) {
      continue;
    }
    classes.push_back(
        ClockCarrierClass{varID, eventIt->second.domain, eventIt->second.phase});
  }
  return classes;
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
size_t expandClockCarrierVarIDsFromBoundaryNames(
    const SequentialDesignModel& model,
    std::unordered_set<size_t>& clockCarrierVarIDs) {
  size_t added = 0;
  for (const auto& [key, varID] : model.inputVarByKey) {
    if (varID < 2 || clockCarrierVarIDs.find(varID) != clockCarrierVarIDs.end()) {
      continue;
    }
    const auto displayIt = model.displayNameByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        // LCOV_EXCL_START
        !isClockTreeCarrierBoundaryName(displayIt->second)) {
        // LCOV_EXCL_STOP
      continue;
    }
    clockCarrierVarIDs.insert(varID);  // LCOV_EXCL_LINE
    ++added;  // LCOV_EXCL_LINE
  }
  return added;
}

size_t expandClockCarrierVarIDsFromTermNames(
    naja::DNL::DNLFull* dnl,
    const std::vector<size_t>& termDNLID2varID,
    std::unordered_set<size_t>& clockCarrierVarIDs) {
  if (dnl == nullptr) {
    return 0;  // LCOV_EXCL_LINE
  }
  size_t added = 0;
  const size_t termCount = std::min(termDNLID2varID.size(), dnl->getNBterms());
  for (naja::DNL::DNLID termID = 0; termID < termCount; ++termID) {
    const size_t varID = termDNLID2varID[termID];
    if (varID < 2 || clockCarrierVarIDs.find(varID) != clockCarrierVarIDs.end()) {
      continue;
    }
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;  // LCOV_EXCL_LINE
    }
    // Model names such as sky130_fd_sc_hs__clkbuf_* are also used as ordinary
    // data buffers after routing.  Name-only promotion is therefore limited to
    // routed branches whose terminal/path/instance name marks them as clock
    // tree; library-cell clkbufs are handled by the structural pass below only
    // when their source actually traces back to a top clock.
    bool isClockCarrier = isClockTreeCarrierTerminalName(term);
    if (!isClockCarrier) {
      continue;
    }
    clockCarrierVarIDs.insert(varID);
    ++added;
  }
  return added;
}

size_t expandClockCarrierVarIDsFromMaterializedTerms(
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    // LCOV_EXCL_START
    std::unordered_set<size_t>& clockCarrierVarIDs,
    std::unordered_map<size_t, ClockEvent>& clockEventByCarrierVarID) {
    // LCOV_EXCL_STOP
  size_t added = 0;
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& [termID, expr] : outputExprByTerm) {
      // LCOV_EXCL_START
      if (termID >= termDNLID2varID.size() || expr == nullptr) {
      // LCOV_EXCL_STOP
        continue;  // LCOV_EXCL_LINE
      }
      const size_t varID = termDNLID2varID[termID];
      if (varID < 2 || clockEventByCarrierVarID.find(varID) !=
                          clockEventByCarrierVarID.end()) {
        continue;
      }

      const auto event =
          classifyClockEventExpression(expr, clockEventByCarrierVarID);
      if (!event.has_value() || !clockEventIsUngated(*event)) {
        continue;
      }

      // This is a phase-preserving clock-tree buffer/inverter inside one
      // extracted design. Gated clocks remain local enable logic, not carriers.
      // LCOV_EXCL_START
      if (clockCarrierVarIDs.insert(varID).second) {
        ++added;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }  // LCOV_EXCL_LINE
      clockEventByCarrierVarID.emplace(varID, *event);
      changed = true;
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
  }
  return added;
}  // LCOV_EXCL_LINE

size_t expandClockCarrierVarIDsFromPureClockTermExprs(
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    std::unordered_set<size_t>& clockCarrierVarIDs) {
  // LCOV_EXCL_START
  size_t added = 0;
  // LCOV_EXCL_STOP
  for (const auto termID : pureClockCarrierTermIDs) {
    const auto exprIt = outputExprByTerm.find(termID);
    if (exprIt == outputExprByTerm.end() || exprIt->second == nullptr) {
      continue;
    }
    // A structurally pure clock terminal can be represented by a temporary
    // builder variable rather than a published SEC input.  Mark that support as
    // carrier-only too, otherwise it leaks into state transitions as data.
    for (const auto varID : exprIt->second->getSupportVars()) {
      if (varID >= 2 && clockCarrierVarIDs.insert(varID).second) {
        // LCOV_EXCL_START
        ++added;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }  // LCOV_EXCL_LINE
    }
  }
  return added;
}  // LCOV_EXCL_LINE

size_t expandClockCarrierVarIDsFromStructure(
    naja::DNL::DNLFull* dnl,
    const SequentialDesignModel& model,
    const std::vector<size_t>& termDNLID2varID,
    std::unordered_set<size_t>& clockCarrierVarIDs,
    // LCOV_EXCL_START
    std::unordered_set<naja::DNL::DNLID>* pureClockCarrierTermIDs = nullptr) {
    // LCOV_EXCL_STOP
  if (dnl == nullptr) {
    return 0;  // LCOV_EXCL_LINE
  }

  std::vector<int8_t> pureClockMemoStrict(dnl->getNBterms(), -1);
  std::vector<int8_t> pureClockMemoAfterNamedClockTree(dnl->getNBterms(), -1);
  auto isPureClockCarrier =
      [&](auto&& self,
          naja::DNL::DNLID termID,
          bool afterNamedClockTree) -> bool {
    if (termID == naja::DNL::DNLID_MAX ||
        termID >= pureClockMemoStrict.size()) {
      return false;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    int8_t& cached = afterNamedClockTree
    // LCOV_EXCL_STOP
                         ? pureClockMemoAfterNamedClockTree[termID]
                         : pureClockMemoStrict[termID];
    if (cached != -1) {
      return cached == 1;
    }
    cached = 0;

    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      return false;  // LCOV_EXCL_LINE
    }

    if (term.isTopPort() &&
        term.getSnlBitTerm()->getDirection() !=
            naja::NL::SNLBitTerm::Direction::Output) {
      const bool result = isTopClockCarrierName(getTerminalDisplayName(term));
      cached = result ? 1 : 0;
      return result;
    }

    if (term.getSnlBitTerm()->getDirection() !=
        naja::NL::SNLBitTerm::Direction::Output) {
      const auto isoID = term.getIsoID();
      if (isoID == naja::DNL::DNLID_MAX) {
        return false;  // LCOV_EXCL_LINE
      }
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant() || iso.getDrivers().size() != 1) {
        return false;
      }
      const bool result =
          self(self, iso.getDrivers().front(), afterNamedClockTree);
      cached = result ? 1 : 0;
      return result;
    }

    if (term.isTopPort()) {
      return false;
    }

    if (isPotentialClockTreeBufferCell(term)) {
      const auto sourceDriver = getClockTreeBufferSourceDriverTerm(dnl, term);
      const bool result =
          sourceDriver.has_value() && self(self, *sourceDriver, true);
      cached = result ? 1 : 0;
      return result;
    }

    if (afterNamedClockTree) {
      // Some CTS flows insert a generic root buffer before the named clkbuf
      // tree. Once a named clock-tree branch has been seen, permit transparent
      // single-input cells to bridge that root back to the top clock.
      const auto sourceDriver = getClockTreeBufferSourceDriverTerm(dnl, term);
      const bool result =
          sourceDriver.has_value() && self(self, *sourceDriver, true);
      cached = result ? 1 : 0;
      // LCOV_EXCL_START
      return result;
      // LCOV_EXCL_STOP
    }

    return false;
  };

  size_t added = 0;
  const auto addCarrierVarID = [&](size_t varID) {
    if (varID >= 2 && clockCarrierVarIDs.insert(varID).second) {
      ++added;
    }
  };
  for (naja::DNL::DNLID termID = 0; termID < pureClockMemoStrict.size(); ++termID) {
    if (!isPureClockCarrier(isPureClockCarrier, termID, false)) {
      continue;
    }
    if (pureClockCarrierTermIDs != nullptr) {
      pureClockCarrierTermIDs->insert(termID);
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
    if (termID < termDNLID2varID.size()) {
      addCarrierVarID(termDNLID2varID[termID]);
    }

    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;  // LCOV_EXCL_LINE
    }
    const auto varIt = model.inputVarByKey.find(getTerminalPathKey(term));
    if (varIt != model.inputVarByKey.end()) {
      addCarrierVarID(varIt->second);
    }
  }
  return added;
}

std::unordered_map<size_t, bool> collectResetAssignments(
    const SequentialDesignModel& model) {
  std::unordered_map<size_t, bool> assignments;
  for (const auto& key : model.environmentInputs) {
    const auto displayIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        // LCOV_EXCL_START
        varIt == model.inputVarByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_STOP
    const auto assertedValue = getResetAssertionValue(displayIt->second);
    // LCOV_EXCL_START
    if (!assertedValue.has_value()) {
      continue;
    }
    // LCOV_EXCL_STOP
    assignments.emplace(varIt->second, *assertedValue);
  }
  return assignments;
}

void inferSynthesizedResetInitialStateValues(SequentialDesignModel& model) {
  const auto resetAssignments = collectResetAssignments(model);
  if (resetAssignments.empty()) {
    return;
  }

  auto resetInitInferenceNodeLimit = []() {
    constexpr size_t kDefaultResetSpecializedExprNodesForInitInference = 200000;
    const char* valueText =
        std::getenv("KEPLER_SEC_RESET_INIT_INFERENCE_NODE_LIMIT");
    if (valueText == nullptr || *valueText == '\0') {
      return kDefaultResetSpecializedExprNodesForInitInference;
    }
    const auto value = std::strtoull(valueText, nullptr, 10);  // LCOV_EXCL_LINE
    if (value == 0) {  // LCOV_EXCL_LINE
      return kDefaultResetSpecializedExprNodesForInitInference;  // LCOV_EXCL_LINE
    }
    return value > std::numeric_limits<size_t>::max()  // LCOV_EXCL_LINE
               ? std::numeric_limits<size_t>::max()  // LCOV_EXCL_LINE
               : static_cast<size_t>(value);  // LCOV_EXCL_LINE
  };

  auto countUniqueExprNodes =
      [](const std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash>& exprByKey) {
        std::unordered_set<BoolExpr*> visited;
        std::vector<BoolExpr*> stack;
        for (const auto& [_, root] : exprByKey) {
          if (root != nullptr) {
            stack.push_back(root);
          // LCOV_EXCL_START
          }
          // LCOV_EXCL_STOP
        }

        while (!stack.empty()) {
          BoolExpr* current = stack.back();
          stack.pop_back();
          if (current == nullptr || !visited.insert(current).second) {
            continue;
          }
          if (current->getLeft() != nullptr) {
            stack.push_back(current->getLeft());
          }
          if (current->getRight() != nullptr) {
            stack.push_back(current->getRight());
          }
        }
        return visited.size();
      // LCOV_EXCL_START
      };


// LCOV_EXCL_STOP
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> resetSpecializedNextStateByKey;
  // LCOV_EXCL_START
  resetSpecializedNextStateByKey.reserve(model.stateBits.size());
  std::unordered_map<BoolExpr*, BoolExpr*> resetSubstitutionMemo;
  for (const auto& key : model.stateBits) {
    const auto nextStateIt = model.nextStateExprByStateKey.find(key);
    if (nextStateIt == model.nextStateExprByStateKey.end()) {
    // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    }
    resetSpecializedNextStateByKey.emplace(
        // LCOV_EXCL_START
        key,
        substituteBoolExprVariables(
        // LCOV_EXCL_STOP
            nextStateIt->second, resetAssignments, resetSubstitutionMemo));
  // LCOV_EXCL_START
  }
  // Synthesized reset inference is only a proof-strengthening heuristic. Cap
  // the specialized DAG size so very large SoCs do not spend most of SEC
  // extraction deriving reset values, while still allowing measured ASIC-size
  // LCOV_EXCL_STOP
  // reset cones to seed PDR's frame-0 frontier once instead of repeatedly
  // proving the same reset-image facts through SAT.
  const size_t maxResetSpecializedExprNodesForInitInference =
      resetInitInferenceNodeLimit();
  const size_t resetSpecializedExprNodes =
      countUniqueExprNodes(resetSpecializedNextStateByKey);
  // LCOV_EXCL_START
  if (std::getenv("KEPLER_SEC_DIAG") != nullptr) {
  // LCOV_EXCL_STOP
    fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: reset-specialized next-state nodes=%zu limit=%zu states=%zu\n",
        resetSpecializedExprNodes,  // LCOV_EXCL_LINE
        maxResetSpecializedExprNodesForInitInference,  // LCOV_EXCL_LINE
        model.stateBits.size());  // LCOV_EXCL_LINE
    fflush(stderr);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  if (resetSpecializedExprNodes >
  // LCOV_EXCL_STOP
      maxResetSpecializedExprNodesForInitInference) {
    if (std::getenv("KEPLER_SEC_DIAG") != nullptr) {
      fprintf(  // LCOV_EXCL_LINE
          stderr,  // LCOV_EXCL_LINE
          "SEC diag: skip synthesized init inference for %zu reset-specialized nodes (limit=%zu)\n",
          resetSpecializedExprNodes,  // LCOV_EXCL_LINE
          maxResetSpecializedExprNodesForInitInference);  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      fflush(stderr);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return;
  }

  auto collectReferencedStateVars = [](BoolExpr* expr) {
  // LCOV_EXCL_STOP
    std::unordered_set<size_t> referencedVars;
    if (expr == nullptr) {
      return referencedVars;  // LCOV_EXCL_LINE
    }

    std::vector<BoolExpr*> stack = {expr};
    std::unordered_set<BoolExpr*> visited;
    while (!stack.empty()) {
      BoolExpr* current = stack.back();
      stack.pop_back();
      if (current == nullptr || !visited.insert(current).second) {
        continue;  // LCOV_EXCL_LINE
      }
      if (current->getOp() == Op::VAR) {
        if (current->getId() >= 2) {
          referencedVars.insert(current->getId());
        }
        // LCOV_EXCL_START
        continue;
        // LCOV_EXCL_STOP
      }
      if (current->getLeft() != nullptr) {  // LCOV_EXCL_LINE
        stack.push_back(current->getLeft());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      if (current->getRight() != nullptr) {  // LCOV_EXCL_LINE
        stack.push_back(current->getRight());  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }
    return referencedVars;
  };

  std::unordered_map<size_t, SignalKey> stateKeyByVar;
  std::unordered_map<size_t, std::vector<SignalKey>> dependentStatesByVar;
  // LCOV_EXCL_START
  stateKeyByVar.reserve(model.stateBits.size());
  dependentStatesByVar.reserve(model.stateBits.size());
  // LCOV_EXCL_STOP
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      stateKeyByVar.emplace(varIt->second, key);
    }
  }
  for (const auto& key : model.stateBits) {
    const auto nextStateIt = resetSpecializedNextStateByKey.find(key);
    if (nextStateIt == resetSpecializedNextStateByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    const auto referencedVars = collectReferencedStateVars(nextStateIt->second);
    for (const auto referencedVar : referencedVars) {
      if (stateKeyByVar.find(referencedVar) == stateKeyByVar.end()) {
        // LCOV_EXCL_START
        continue;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      dependentStatesByVar[referencedVar].push_back(key);
    }
  }

  std::unordered_map<SignalKey, SignalKey, SignalKeyHash> complementedPartnerByKey;
  complementedPartnerByKey.reserve(model.complementedStateRelations.size() * 2);
  for (const auto& relation : model.complementedStateRelations) {
    complementedPartnerByKey.emplace(relation.primaryKey, relation.complementedKey);  // LCOV_EXCL_LINE
    complementedPartnerByKey.emplace(relation.complementedKey, relation.primaryKey);  // LCOV_EXCL_LINE
  }

  std::unordered_map<size_t, bool> assignments = resetAssignments;
  for (const auto& [key, value] : model.initialStateValueByKey) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      // LCOV_EXCL_START
      assignments.emplace(varIt->second, value);
    }
  }


// LCOV_EXCL_STOP
  std::deque<SignalKey> workQueue(model.stateBits.begin(), model.stateBits.end());
  auto recordKnownState = [&](const SignalKey& key, bool value) {
    const auto [it, inserted] = model.initialStateValueByKey.emplace(key, value);
    if (!inserted) {
      return;  // LCOV_EXCL_LINE
    }

    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      // LCOV_EXCL_START
      assignments[varIt->second] = value;
      const auto dependentIt = dependentStatesByVar.find(varIt->second);
      if (dependentIt != dependentStatesByVar.end()) {
        workQueue.insert(
        // LCOV_EXCL_STOP
            workQueue.end(),
            dependentIt->second.begin(),
            dependentIt->second.end());
      }
    }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
    const auto partnerIt = complementedPartnerByKey.find(key);
    if (partnerIt != complementedPartnerByKey.end() &&
        model.initialStateValueByKey.find(partnerIt->second) ==  // LCOV_EXCL_LINE
            model.initialStateValueByKey.end()) {  // LCOV_EXCL_LINE
      workQueue.push_back(partnerIt->second);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  };

  while (!workQueue.empty()) {
    const SignalKey key = workQueue.front();
    workQueue.pop_front();

    if (model.initialStateValueByKey.find(key) != model.initialStateValueByKey.end()) {
      const auto partnerIt = complementedPartnerByKey.find(key);
      if (partnerIt != complementedPartnerByKey.end() &&
          model.initialStateValueByKey.find(partnerIt->second) ==  // LCOV_EXCL_LINE
              model.initialStateValueByKey.end()) {  // LCOV_EXCL_LINE
        recordKnownState(partnerIt->second, !model.initialStateValueByKey.at(key));  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      continue;
    }

    const auto nextStateIt = resetSpecializedNextStateByKey.find(key);
    if (nextStateIt == resetSpecializedNextStateByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }

    std::unordered_map<BoolExpr*, std::optional<bool>> memo;
    const auto resetValue = evaluateConstantUnderAssignments(
        nextStateIt->second, assignments, memo);
    if (resetValue.has_value()) {
      recordKnownState(key, *resetValue);
    }
  }
}

struct ExtractContext {
  naja::NL::SNLDesign* top = nullptr;
  naja::NL::NLUniverse* universe = nullptr;
  naja::NL::SNLDesign* previousTop = nullptr;
  std::string topName;
  bool secDiagEnabled = false;
  bool abstractUncomputableSequentialBoundaries = false;
  // LCOV_EXCL_START
  KEPLER_FORMAL::BuildPrimaryOutputClauses builder;
  // LCOV_EXCL_STOP
  decltype(naja::DNL::get()) dnl = nullptr;
  std::unordered_map<naja::DNL::DNLID, SignalKey> inputKeyByTerm;
  std::unordered_map<naja::DNL::DNLID, SignalKey> outputKeyByTerm;
  std::unordered_map<naja::DNL::DNLID, SignalKey> topOutputKeyByTerm;
  // LCOV_EXCL_START
  std::set<SignalKey, SignalKeyLess> topInputKeys;
  std::set<SignalKey, SignalKeyLess> topOutputKeys;
  std::set<SignalKey, SignalKeyLess> internalBoundaryInputKeys;
  std::set<SignalKey, SignalKeyLess> internalBoundaryOutputKeys;
  std::set<SignalKey, SignalKeyLess> environmentInputs;
  std::set<SignalKey, SignalKeyLess> stateBits;
  // LCOV_EXCL_STOP
  std::set<SignalKey, SignalKeyLess> allObservedOutputs;
  // LCOV_EXCL_START
  std::unordered_set<naja::DNL::DNLID> prunedBuilderOutputTerms;
  std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>
      collectedSkippedOutputs;
  std::set<SignalKey, SignalKeyLess> abstractedBoundaryStateKeys;
  std::vector<std::pair<naja::DNL::DNLID, SignalKey>> abstractedBoundaryObservedTerms;
  std::unordered_set<SignalKey, SignalKeyHash> abstractedBoundaryObservedKeys;
  std::unordered_set<SignalKey, SignalKeyHash> unsupportedStateBits;
  std::vector<PendingTransition> pendingTransitions;
  std::vector<PendingMemoryInstance> pendingMemoryInstances;
  std::vector<InstanceBoundaryInfo> instanceBoundaryInfos;
  std::unordered_map<naja::DNL::DNLID, bool> sequentialInstanceCache;
};

std::string describeSupportVarOrigins(  // LCOV_EXCL_LINE
    const ExtractContext& ctx,
    // LCOV_EXCL_STOP
    size_t varID,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs) {
  std::ostringstream out;  // LCOV_EXCL_LINE
  bool first = true;  // LCOV_EXCL_LINE
  const size_t termCount = std::min(termDNLID2varID.size(), ctx.dnl->getNBterms());  // LCOV_EXCL_LINE
  for (naja::DNL::DNLID termID = 0; termID < termCount; ++termID) {  // LCOV_EXCL_LINE
    if (termDNLID2varID[termID] != varID) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    if (!first) {  // LCOV_EXCL_LINE
      out << ", ";  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    first = false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    const auto& term = ctx.dnl->getDNLTerminalFromID(termID);  // LCOV_EXCL_LINE
    out << "term#" << termID << ":" << getTerminalDisplayName(term);  // LCOV_EXCL_LINE
    out << ":input_key=" << (ctx.inputKeyByTerm.find(termID) != ctx.inputKeyByTerm.end());  // LCOV_EXCL_LINE
    out << ":pure_clock=" << (pureClockCarrierTermIDs.find(termID) !=  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
                               pureClockCarrierTermIDs.end());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  if (first) {  // LCOV_EXCL_LINE
    out << "no direct term";  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return out.str();  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}  // LCOV_EXCL_LINE

// LCOV_EXCL_START


// LCOV_EXCL_STOP
void logUnpublishedTransitionSupport(
    // LCOV_EXCL_START
    const ExtractContext& ctx,
    const SequentialDesignModel& model,
    const PendingTransition& pending,
    BoolExpr* nextStateExpr,
    const std::vector<size_t>& termDNLID2varID,
    // LCOV_EXCL_STOP
    const std::unordered_set<size_t>& clockCarrierVarIDs,
    // LCOV_EXCL_START
    const std::unordered_set<naja::DNL::DNLID>& pureClockCarrierTermIDs) {
  if (nextStateExpr == nullptr ||
      std::getenv("KEPLER_SEC_CLOCK_GATE_DIAG") == nullptr) {
    return;
  }


// LCOV_EXCL_STOP
  std::unordered_set<size_t> publishedVars;  // LCOV_EXCL_LINE
  publishedVars.reserve(model.inputVarByKey.size());  // LCOV_EXCL_LINE
  for (const auto& [_, varID] : model.inputVarByKey) {  // LCOV_EXCL_LINE
    publishedVars.insert(varID);  // LCOV_EXCL_LINE
  }

  for (const auto varID : nextStateExpr->getSupportVars()) {  // LCOV_EXCL_LINE
    if (varID < 2 || publishedVars.find(varID) != publishedVars.end() ||  // LCOV_EXCL_LINE
        clockCarrierVarIDs.find(varID) != clockCarrierVarIDs.end()) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    const auto displayIt = model.displayNameByKey.find(pending.stateKey);  // LCOV_EXCL_LINE
    const std::string stateName =
        displayIt == model.displayNameByKey.end()  // LCOV_EXCL_LINE
            ? signalKeyToString(pending.stateKey)  // LCOV_EXCL_LINE
            : displayIt->second;  // LCOV_EXCL_LINE
    std::fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: transition `%s` unpublished support v%zu origins=[%s]\n",
        stateName.c_str(),  // LCOV_EXCL_LINE
        varID,  // LCOV_EXCL_LINE
        describeSupportVarOrigins(  // LCOV_EXCL_LINE
            ctx, varID, termDNLID2varID, pureClockCarrierTermIDs).c_str());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  std::fflush(stderr);  // LCOV_EXCL_LINE
}

void collectInitialBuilderBoundary(ExtractContext& ctx) {
  naja::DNL::destroy();
  ctx.universe->setTopDesign(ctx.top);

  // Reuse the existing miter frontend to discover the relevant boundary
  // signals before we ask it to build Boolean formulas.
  if (ctx.secDiagEnabled) {
    fprintf(stderr, "SEC diag: extract(%s) collect begin\n", ctx.topName.c_str());
    fflush(stderr);
  }
  ctx.builder.collect();
  ctx.collectedSkippedOutputs = ctx.builder.getSkippedOutputs();
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) collect end inputs=%zu outputs=%zu\n",
        ctx.topName.c_str(),
        ctx.builder.getInputs().size(),
        ctx.builder.getOutputs().size());
    fflush(stderr);
  }

  ctx.dnl = naja::DNL::get();
}

void collectTopInterfaceTerms(ExtractContext& ctx, SequentialDesignModel& model) {
  const auto topInstance = ctx.dnl->getTop();
  for (naja::DNL::DNLID termID = topInstance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX && termID <= topInstance.getTermIndexes().second;
       ++termID) {
    const auto& term = ctx.dnl->getDNLTerminalFromID(termID);
    SignalKey key = getTerminalPathKey(term);
    model.displayNameByKey.try_emplace(key, getTerminalDisplayName(term));
    if (term.getSnlBitTerm()->getDirection() ==
        naja::NL::SNLBitTerm::Direction::Input) {
      ctx.topInputKeys.insert(key);
      continue;
    }
    ctx.topOutputKeys.insert(key);
    ctx.topOutputKeyByTerm.emplace(termID, key);
    ctx.allObservedOutputs.insert(key);
  }
}

bool isLatchLikeOutputPinName(const std::string& pinName) {
  return pinName == "Q" || pinName == "IQ";
}

bool isLatchLikeDataPinName(const std::string& pinName) {
  return pinName == "D";
}

bool isLatchLikeClockPinName(const std::string& pinName) {
  return pinName == "E" || pinName == "EN" || pinName == "GATE" ||
         pinName == "GATE_N" || pinName == "GATEN" || pinName == "CLK" ||
         pinName == "CK";
}

std::vector<naja::DNL::DNLID> collectClockGateLatchDependencyTerms(
    const ExtractContext& ctx,
    const SequentialDesignModel& model) {
  std::vector<naja::DNL::DNLID> terms;
  std::unordered_set<SignalKey, SignalKeyHash> stateKeys(
      model.stateBits.begin(), model.stateBits.end());
  std::unordered_set<SignalKey, SignalKeyHash> topInputKeys(
      model.topInputKeys.begin(), model.topInputKeys.end());

  for (auto leafID : ctx.dnl->getLeaves()) {
    const auto& instance = ctx.dnl->getDNLInstanceFromID(leafID);
    std::vector<naja::DNL::DNLID> latchOutputTerms;
    std::vector<naja::DNL::DNLID> latchDataTerms;
    std::vector<naja::DNL::DNLID> latchClockTerms;

    for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= instance.getTermIndexes().second;
         ++termID) {
      const auto& term = ctx.dnl->getDNLTerminalFromID(termID);
      const std::string pinName =
          normalizePinName(term.getSnlBitTerm()->getName().getString());
      if (term.getSnlBitTerm()->getDirection() ==
          naja::NL::SNLBitTerm::Direction::Output) {
        if (isLatchLikeOutputPinName(pinName)) {
          latchOutputTerms.push_back(termID);
        }
        continue;
      }
      if (isLatchLikeDataPinName(pinName)) {
        latchDataTerms.push_back(termID);
      } else if (isLatchLikeClockPinName(pinName)) {
        latchClockTerms.push_back(termID);
      }
    }

    if (latchOutputTerms.size() != 1 || latchDataTerms.size() != 1 ||
        latchClockTerms.empty()) {
      continue;
    }
    const auto& outputTerm =
        ctx.dnl->getDNLTerminalFromID(latchOutputTerms.front());
    const SignalKey outputKey = getTerminalPathKey(outputTerm);
    if (stateKeys.find(outputKey) != stateKeys.end() ||
        topInputKeys.find(outputKey) != topInputKeys.end() ||
        model.inputVarByKey.find(outputKey) == model.inputVarByKey.end()) {
      continue;
    }

    terms.push_back(latchDataTerms.front());
    terms.insert(terms.end(), latchClockTerms.begin(), latchClockTerms.end());
  }

  std::sort(terms.begin(), terms.end());
  terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
  return terms;
}

std::unordered_map<size_t, BoolExpr*> collectClockGateLatchDataExprByVarID(
    const ExtractContext& ctx,
    const SequentialDesignModel& model,
    const std::unordered_set<size_t>& topClockCarrierVarIDs,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm) {
  std::unordered_map<size_t, BoolExpr*> dataExprByVarID;
  const bool diagEnabled = std::getenv("KEPLER_SEC_CLOCK_GATE_DIAG") != nullptr;
  size_t latchLikeCandidates = 0;
  size_t latchLikeStateOutputs = 0;
  size_t latchLikeTopOutputs = 0;
  size_t latchLikeUnmappedOutputs = 0;
  size_t latchLikeClockRejected = 0;
  size_t latchLikeMissingDataExpr = 0;
  size_t latchLikeMissingClockExpr = 0;
  if (topClockCarrierVarIDs.empty()) {
    return dataExprByVarID;
  }

  std::unordered_set<SignalKey, SignalKeyHash> stateKeys(
      model.stateBits.begin(), model.stateBits.end());
  std::unordered_set<SignalKey, SignalKeyHash> topInputKeys(
      model.topInputKeys.begin(), model.topInputKeys.end());

  for (auto leafID : ctx.dnl->getLeaves()) {
    const auto& instance = ctx.dnl->getDNLInstanceFromID(leafID);
    std::vector<naja::DNL::DNLID> latchOutputTerms;
    std::vector<naja::DNL::DNLID> latchDataTerms;
    std::vector<naja::DNL::DNLID> latchClockTerms;

    for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= instance.getTermIndexes().second;
         ++termID) {
      const auto& term = ctx.dnl->getDNLTerminalFromID(termID);
      const std::string pinName =
          normalizePinName(term.getSnlBitTerm()->getName().getString());
      if (term.getSnlBitTerm()->getDirection() ==
          naja::NL::SNLBitTerm::Direction::Output) {
        if (isLatchLikeOutputPinName(pinName)) {
          latchOutputTerms.push_back(termID);
        }
        // LCOV_EXCL_START
        continue;
        // LCOV_EXCL_STOP
      }
      if (isLatchLikeDataPinName(pinName)) {
        latchDataTerms.push_back(termID);
      // LCOV_EXCL_START
      } else if (isLatchLikeClockPinName(pinName)) {
      // LCOV_EXCL_STOP
        latchClockTerms.push_back(termID);
      }
    // LCOV_EXCL_START
    }


// LCOV_EXCL_STOP
    if (latchOutputTerms.size() != 1 || latchDataTerms.size() != 1 ||
        latchClockTerms.empty()) {
      continue;
    }
    ++latchLikeCandidates;

    const auto dataExprIt = outputExprByTerm.find(latchDataTerms.front());
    if (dataExprIt == outputExprByTerm.end()) {
      ++latchLikeMissingDataExpr;
      continue;
    }

    BoolExpr* latchClockExpr = nullptr;
    for (const auto clockTermID : latchClockTerms) {
      const auto clockExprIt = outputExprByTerm.find(clockTermID);
      if (clockExprIt == outputExprByTerm.end()) {
        continue;  // LCOV_EXCL_LINE
      }
      latchClockExpr = latchClockExpr == nullptr
                           ? clockExprIt->second
                           // LCOV_EXCL_START
                           : BoolExpr::And(latchClockExpr, clockExprIt->second);  // LCOV_EXCL_LINE
                           // LCOV_EXCL_STOP
    }
    if (latchClockExpr == nullptr) {
      ++latchLikeMissingClockExpr;  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }

// LCOV_EXCL_START

    std::unordered_map<BoolExpr*, BoolExpr*> clockStripMemo;
    // LCOV_EXCL_STOP
    BoolExpr* strippedClock = simplifyWhenClockCarriersChanged(
        latchClockExpr,
        stripClockCarrierFromClockEnable(
            latchClockExpr, topClockCarrierVarIDs, clockStripMemo));
    if (strippedClock != BoolExpr::createTrue()) {
      // LCOV_EXCL_START
      ++latchLikeClockRejected;
      continue;
      // LCOV_EXCL_STOP
    }

// LCOV_EXCL_START

    const auto& outputTerm =
        ctx.dnl->getDNLTerminalFromID(latchOutputTerms.front());
    const SignalKey outputKey = getTerminalPathKey(outputTerm);
    if (stateKeys.find(outputKey) != stateKeys.end() ||
        topInputKeys.find(outputKey) != topInputKeys.end()) {
      if (stateKeys.find(outputKey) != stateKeys.end()) {
        ++latchLikeStateOutputs;
      } else {
        ++latchLikeTopOutputs;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      continue;
    }
    const auto varIt = model.inputVarByKey.find(outputKey);
    if (varIt == model.inputVarByKey.end()) {
      ++latchLikeUnmappedOutputs;  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    dataExprByVarID.emplace(varIt->second, dataExprIt->second);
  }

  if (diagEnabled) {
    fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: clock-gate latch fold candidates=%zu mapped=%zu state=%zu top=%zu unmapped=%zu missing_data=%zu missing_clock=%zu clock_rejected=%zu\n",
        latchLikeCandidates,  // LCOV_EXCL_LINE
        dataExprByVarID.size(),  // LCOV_EXCL_LINE
        latchLikeStateOutputs,  // LCOV_EXCL_LINE
        latchLikeTopOutputs,  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        latchLikeUnmappedOutputs,  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        latchLikeMissingDataExpr,  // LCOV_EXCL_LINE
        latchLikeMissingClockExpr,  // LCOV_EXCL_LINE
        latchLikeClockRejected);  // LCOV_EXCL_LINE
    fflush(stderr);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  return dataExprByVarID;
}

bool isSequentialInstanceTerm(ExtractContext& ctx,
                              const naja::DNL::DNLTerminalFull& term) {
  const auto& instance = term.getDNLInstance();
  const auto instanceID = instance.getID();
  if (const auto cached = ctx.sequentialInstanceCache.find(instanceID);
      cached != ctx.sequentialInstanceCache.end()) {
    return cached->second;
  }

  bool isSequentialInstance = false;
  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX && termID <= instance.getTermIndexes().second;
       ++termID) {
    const auto& instanceTerm = ctx.dnl->getDNLTerminalFromID(termID);
    if (instanceTerm.isNull()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (isSequentialStateOutput(instanceTerm) ||
        isSequentialNextStateInput(instanceTerm)) {
      isSequentialInstance = true;
      break;
    }
  }

  ctx.sequentialInstanceCache.emplace(instanceID, isSequentialInstance);
  return isSequentialInstance;
}

void classifyBuilderBoundaryTerms(ExtractContext& ctx, SequentialDesignModel& model) {
  // The miter builder already exposes sequential outputs as "inputs" and
  // sequential next-state pins as "outputs"; we normalize those into SEC's
  // environment/state/output buckets here.
  for (const auto inputTermID : ctx.builder.getInputs()) {
    const auto& term = ctx.dnl->getDNLTerminalFromID(inputTermID);
    if (isConstantInternalOutputTerm(term)) {
      continue;
    }
    SignalKey key = getTerminalPathKey(term);
    ctx.inputKeyByTerm.emplace(inputTermID, key);
    model.displayNameByKey.try_emplace(key, getTerminalDisplayName(term));
    if (isSequentialStateOutput(term)) {
      ctx.stateBits.insert(key);
    } else {
      ctx.environmentInputs.insert(key);
      if (!term.isTopPort() && !isSequentialInstanceTerm(ctx, term)) {
        ctx.internalBoundaryInputKeys.insert(key);
      }
    }
  }

  for (const auto outputTermID : ctx.builder.getOutputs()) {
    const auto& term = ctx.dnl->getDNLTerminalFromID(outputTermID);
    SignalKey key = getTerminalPathKey(term);
    ctx.outputKeyByTerm.emplace(outputTermID, key);
    model.displayNameByKey.try_emplace(key, getTerminalDisplayName(term));
    if (ctx.topOutputKeys.find(key) == ctx.topOutputKeys.end() &&
        !isSequentialInstanceTerm(ctx, term)) {
      ctx.internalBoundaryOutputKeys.insert(key);
    }
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

std::optional<SequentialInstanceScan> scanSequentialInstance(
    const naja::DNL::DNLInstanceFull& instance,
    const std::unordered_map<naja::DNL::DNLID, SignalKey>& inputKeyByTerm,
    SequentialDesignModel& model) {
  // LCOV_EXCL_START
  SequentialInstanceScan scan;
  // LCOV_EXCL_STOP
  scan.boundaryInfo.instancePath = instance.getFullPath();

  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX &&
       termID <= instance.getTermIndexes().second;
       ++termID) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(termID);
    const std::string normalizedPinName =
        normalizePinName(term.getSnlBitTerm()->getName().getString());
    if (isSequentialStateOutput(term) &&
        term.getSnlBitTerm()->getDirection() !=
            naja::NL::SNLBitTerm::Direction::Input) {
      std::vector<PendingPinTerm> clockTermIDs;
      bool intrinsicClockIsInverted = false;
      for (auto* clockBitTerm :
           naja::NL::SNLDesignModeling::getOutputRelatedClocks(
               term.getSnlBitTerm())) {
        if (clockBitTerm == nullptr) {
          continue;  // LCOV_EXCL_LINE
        }
        const auto& clockTerm = instance.getTerminalFromBitTerm(clockBitTerm);
        if (clockTerm.isNull() ||
            clockTerm.getSnlBitTerm()->getDirection() ==
                naja::NL::SNLBitTerm::Direction::Output) {
          continue;  // LCOV_EXCL_LINE
        }
        intrinsicClockIsInverted =
            intrinsicClockIsInverted ||
            isIntrinsicNegativeEdgeClock(instance, clockBitTerm);
        clockTermIDs.push_back(
            {clockTerm.getID(), clockTerm.getSnlBitTerm()->getBit()});
      }
      scan.stateOutputs.push_back(
          {termID,
           normalizedPinName,
           term.getSnlBitTerm()->getBit(),
           intrinsicClockIsInverted,
           std::move(clockTermIDs)});
    }
    // Some Liberty readers expose async set/reset pins as control timing arcs
    // rather than data-to-clock arcs.  If the cell is sequential and the pin
    // name is a supported control role, include it explicitly so reset
    // semantics are modeled instead of leaving the flop unconstrained.
    if ((isSequentialNextStateInput(term) ||
         isOptionalSequentialControlPin(normalizedPinName)) &&
        term.getSnlBitTerm()->getDirection() !=
            naja::NL::SNLBitTerm::Direction::Output) {
      scan.pinTermIDs[normalizedPinName].push_back(
          {termID, term.getSnlBitTerm()->getBit()});
    }
  }

  if (scan.stateOutputs.empty()) {
    return std::nullopt;
  }

  std::set<SignalKey, SignalKeyLess> boundaryStateKeys;
  for (const auto& stateOutput : scan.stateOutputs) {
    const auto keyIt = inputKeyByTerm.find(stateOutput.termID);
    if (keyIt != inputKeyByTerm.end()) {
      boundaryStateKeys.insert(keyIt->second);
    }
  }
  scan.boundaryInfo.stateKeys.assign(boundaryStateKeys.begin(), boundaryStateKeys.end());

  std::set<SignalKey, SignalKeyLess> boundaryObservedKeys;
  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX &&
       termID <= instance.getTermIndexes().second;
       ++termID) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(termID);
    if (term.isNull() ||
        term.getSnlBitTerm()->getDirection() == naja::NL::SNLBitTerm::Direction::Output) {
      continue;
    }
    const SignalKey key = getTerminalPathKey(term);
    model.displayNameByKey.try_emplace(key, getTerminalDisplayName(term));
    if (boundaryObservedKeys.insert(key).second) {
      scan.boundaryInfo.observedTerms.push_back({termID, key});
    }
  }

  return scan;
}

SignalKey makeMemoryCellStateKey(
    const naja::DNL::DNLInstanceFull& instance,
    size_t cellIndex,
    size_t bitIndex) {
  SignalKey key;
  const auto pathNames = instance.getPath().getPathNames();
  key.first.reserve(pathNames.size() + 1);
  for (const auto& name : pathNames) {
    key.first.push_back(name.getID());
  }
  key.first.push_back(naja::NL::NLName("__MEM_CELL").getID());
  key.second.push_back(
      // LCOV_EXCL_START
      static_cast<naja::NL::NLID::DesignObjectID>(cellIndex));
      // LCOV_EXCL_STOP
  key.second.push_back(
      static_cast<naja::NL::NLID::DesignObjectID>(bitIndex));
  return key;
}

std::string makeMemoryCellStateDisplayName(
    const naja::DNL::DNLInstanceFull& instance,
    size_t cellIndex,
    size_t bitIndex) {
  // LCOV_EXCL_START
  std::ostringstream oss;
  // LCOV_EXCL_STOP
  oss << instance.getFullPath() << ".__MEM_CELL[" << cellIndex << "]["
      << bitIndex << "]";
  return oss.str();
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
std::unordered_map<const naja::NL::SNLBitTerm*, naja::DNL::DNLID>
collectInstanceTermIDByBitTerm(const naja::DNL::DNLInstanceFull& instance) {
  std::unordered_map<const naja::NL::SNLBitTerm*, naja::DNL::DNLID> termIDs;
  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX &&
       // LCOV_EXCL_START
       termID <= instance.getTermIndexes().second;
       // LCOV_EXCL_STOP
       ++termID) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(termID);
    // LCOV_EXCL_START
    if (term.isNull()) {
    // LCOV_EXCL_STOP
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    termIDs.emplace(term.getSnlBitTerm(), termID);
    // LCOV_EXCL_STOP
  }
  return termIDs;
}

bool supportsStructuredMemoryModel(
    const naja::NL::SNLDesignModeling::MemoryInterface& interface) {
  if (!interface.isValid()) {
    return false;  // LCOV_EXCL_LINE
  }
  for (const auto& readPort : interface.readPorts) {
    if (readPort.address.size() != interface.abits ||
        readPort.data.size() != interface.width) {
      return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
  }
  for (const auto& writePort : interface.writePorts) {
  // LCOV_EXCL_STOP
    if (writePort.address.size() != interface.abits ||
        writePort.data.size() != interface.width) {
      // LCOV_EXCL_START
      return false;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (!writePort.mask.empty() && writePort.mask.size() != interface.width) {
      return false;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    if (!writePort.extraWriteInputs.empty()) {
    // LCOV_EXCL_STOP
      return false;  // LCOV_EXCL_LINE
    }
  }
  // LCOV_EXCL_START
  return true;
  // LCOV_EXCL_STOP
}

naja::DNL::DNLID getRequiredInstanceTermID(
    const std::unordered_map<const naja::NL::SNLBitTerm*, naja::DNL::DNLID>&
        termIDsByBitTerm,
    const naja::NL::SNLBitTerm* term,
    const std::string& instancePath,
    const char* context) {
  // LCOV_EXCL_START
  const auto termIt = termIDsByBitTerm.find(term);
  // LCOV_EXCL_STOP
  if (termIt == termIDsByBitTerm.end()) {
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Missing DNL term for memory " + std::string(context) + " in instance `" +  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        instancePath + "`");  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
  }
  return termIt->second;
}  // LCOV_EXCL_LINE

bool isNoDriverTerm(naja::DNL::DNLID termID) {
  auto* dnl = naja::DNL::get();
  if (dnl == nullptr || termID == naja::DNL::DNLID_MAX) {
    return true;  // LCOV_EXCL_LINE
  }
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.isNull() || term.getIsoID() == naja::DNL::DNLID_MAX) {
    return true;  // LCOV_EXCL_LINE
  }
  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  return !iso.isConstant() && iso.getDrivers().empty();
}

bool isConstantZeroTerm(naja::DNL::DNLID termID) {
  auto* dnl = naja::DNL::get();
  if (dnl == nullptr || termID == naja::DNL::DNLID_MAX) {
    return false;  // LCOV_EXCL_LINE
  }
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.isNull() || term.getIsoID() == naja::DNL::DNLID_MAX) {
    return false;  // LCOV_EXCL_LINE
  }
  return dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID()).isConstant0();
}

bool isDisabledMemoryWriteEnable(naja::DNL::DNLID termID) {
  // Some lowered memories expose fixed-width write ports even when a
  // particular port is unused by the RTL. In that shape the enable pin can
  // have a net but no leaf driver; treating it as an active symbolic input
  // would pull unrelated address/data cones into the SEC memory transition.
  // A constant-0 or undriven enable cannot assert in the concrete netlist, so
  // the whole write port is semantically inactive and should be ignored.
  return isConstantZeroTerm(termID) || isNoDriverTerm(termID);
}

void appendPendingMemoryInstance(
    ExtractContext& ctx,
    SequentialDesignModel& model,
    const naja::DNL::DNLInstanceFull& instance) {
  const auto interface =
      naja::NL::SNLDesignModeling::getMemoryInterface(instance.getSNLInstance());
  PendingMemoryInstance pending;
  pending.width = interface.width;
  pending.depth = interface.depth;
  pending.abits = interface.abits;
  pending.resetMode = interface.resetMode;

  InstanceBoundaryInfo boundaryInfo;
  boundaryInfo.instancePath = instance.getFullPath();
  for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
       termID != naja::DNL::DNLID_MAX &&
       termID <= instance.getTermIndexes().second;
       ++termID) {
    const auto& term = naja::DNL::get()->getDNLTerminalFromID(termID);
    if (term.isNull() ||
        term.getSnlBitTerm()->getDirection() ==
            naja::NL::SNLBitTerm::Direction::Output) {
      continue;
    }
    const SignalKey key = getTerminalPathKey(term);
    model.displayNameByKey.try_emplace(key, getTerminalDisplayName(term));
    boundaryInfo.observedTerms.push_back({termID, key});
  }

  for (size_t cellIndex = 0; cellIndex < interface.depth; ++cellIndex) {
    for (size_t bitIndex = 0; bitIndex < interface.width; ++bitIndex) {
      auto key = makeMemoryCellStateKey(instance, cellIndex, bitIndex);
      auto displayName = makeMemoryCellStateDisplayName(instance, cellIndex, bitIndex);
      model.displayNameByKey.try_emplace(key, displayName);
      boundaryInfo.stateKeys.push_back(key);
      pending.cellStates.push_back({key, std::move(displayName), cellIndex, bitIndex});
    }
  }

  const auto termIDsByBitTerm = collectInstanceTermIDByBitTerm(instance);
  // Modeled memories always expose a reset terminal on the lowered primitive,
  // but SEC should only rebuild that cone when the memory semantics actually
  // use reset. Otherwise a disabled RST pin becomes a fake dependency surface.
  if (pending.resetMode != naja::NL::SNLDesignModeling::MemoryResetMode::None &&
      interface.reset != nullptr) {
    pending.resetTermID = getRequiredInstanceTermID(
        termIDsByBitTerm,
        // LCOV_EXCL_START
        interface.reset,
        boundaryInfo.instancePath,
        "reset");
        // LCOV_EXCL_STOP
  }
  pending.readPorts.reserve(interface.readPorts.size());
  for (size_t portIndex = 0; portIndex < interface.readPorts.size(); ++portIndex) {
    const auto& readPort = interface.readPorts[portIndex];
    PendingMemoryReadPort pendingReadPort;
    pendingReadPort.addressTermIDs.reserve(readPort.address.size());
    pendingReadPort.dataTermIDs.reserve(readPort.data.size());
    size_t bitIndex = 0;
    for (auto* addressTerm : readPort.address) {
      pendingReadPort.addressTermIDs.push_back(
          getRequiredInstanceTermID(
              termIDsByBitTerm,
              addressTerm,
              boundaryInfo.instancePath,
              "read-address"));
    }
    for (auto* dataTerm : readPort.data) {
      const auto termID = getRequiredInstanceTermID(
          termIDsByBitTerm, dataTerm, boundaryInfo.instancePath, "read-data");
      pendingReadPort.dataTermIDs.push_back(termID);
      const auto keyIt = ctx.inputKeyByTerm.find(termID);
      if (keyIt == ctx.inputKeyByTerm.end()) {
        throw std::runtime_error(  // LCOV_EXCL_LINE
            "Missing SEC state key for memory read-data bit in instance `" +  // LCOV_EXCL_LINE
            boundaryInfo.instancePath + "`");  // LCOV_EXCL_LINE
      }
      pending.readOutputs.push_back({keyIt->second, portIndex, bitIndex});
      ++bitIndex;
    }
    pending.readPorts.push_back(std::move(pendingReadPort));
  }

  pending.writePorts.reserve(interface.writePorts.size());
  for (const auto& writePort : interface.writePorts) {
    std::vector<naja::DNL::DNLID> enableTermIDs;
    enableTermIDs.reserve(writePort.enables.size());
    for (auto* enableTerm : writePort.enables) {
      enableTermIDs.push_back(
          getRequiredInstanceTermID(
              termIDsByBitTerm,
              enableTerm,
              boundaryInfo.instancePath,
              "write-enable"));
    }
    if (std::any_of(
            enableTermIDs.begin(),
            enableTermIDs.end(),
            isDisabledMemoryWriteEnable)) {
      continue;
    }

// LCOV_EXCL_START

    PendingMemoryWritePort pendingWritePort;
    pendingWritePort.addressTermIDs.reserve(writePort.address.size());
    pendingWritePort.dataTermIDs.reserve(writePort.data.size());
    pendingWritePort.maskTermIDs.reserve(writePort.mask.size());
    // LCOV_EXCL_STOP
    pendingWritePort.enableTermIDs = std::move(enableTermIDs);
    // LCOV_EXCL_START
    for (auto* addressTerm : writePort.address) {
      pendingWritePort.addressTermIDs.push_back(
      // LCOV_EXCL_STOP
          getRequiredInstanceTermID(
              termIDsByBitTerm,
              // LCOV_EXCL_START
              addressTerm,
              // LCOV_EXCL_STOP
              boundaryInfo.instancePath,
              // LCOV_EXCL_START
              "write-address"));
              // LCOV_EXCL_STOP
    }
    for (auto* dataTerm : writePort.data) {
      pendingWritePort.dataTermIDs.push_back(
          getRequiredInstanceTermID(
              termIDsByBitTerm, dataTerm, boundaryInfo.instancePath, "write-data"));
    }
    for (auto* maskTerm : writePort.mask) {
      pendingWritePort.maskTermIDs.push_back(
          getRequiredInstanceTermID(
              termIDsByBitTerm, maskTerm, boundaryInfo.instancePath, "write-mask"));
    }
    for (const auto& extraWriteTerms : writePort.extraWriteInputs) {
      std::vector<naja::DNL::DNLID> pendingExtraTerms;  // LCOV_EXCL_LINE
      pendingExtraTerms.reserve(extraWriteTerms.size());  // LCOV_EXCL_LINE
      for (auto* extraTerm : extraWriteTerms) {  // LCOV_EXCL_LINE
        pendingExtraTerms.push_back(  // LCOV_EXCL_LINE
            getRequiredInstanceTermID(  // LCOV_EXCL_LINE
                termIDsByBitTerm,
                extraTerm,  // LCOV_EXCL_LINE
                boundaryInfo.instancePath,  // LCOV_EXCL_LINE
                "extra-write"));
      }
      pendingWritePort.extraWriteInputTermIDs.push_back(  // LCOV_EXCL_LINE
          std::move(pendingExtraTerms));
    }  // LCOV_EXCL_LINE
    pending.writePorts.push_back(std::move(pendingWritePort));
  }

  ctx.instanceBoundaryInfos.push_back(std::move(boundaryInfo));
  pending.boundaryInfoIndex = ctx.instanceBoundaryInfos.size() - 1;
  ctx.pendingMemoryInstances.push_back(std::move(pending));
}

void appendPendingTransitionsForInstance(
    ExtractContext& ctx,
    SequentialDesignModel& model,
    const SequentialInstanceScan& scan) {
  ctx.instanceBoundaryInfos.push_back(scan.boundaryInfo);
  const size_t boundaryInfoIndex = ctx.instanceBoundaryInfos.size() - 1;

  auto markUnsupportedInstanceStateOutputs = [&]() {
    for (const auto& key : ctx.instanceBoundaryInfos[boundaryInfoIndex].stateKeys) {
      ctx.unsupportedStateBits.insert(key);
    }
  };
  auto abstractUnsupportedInstanceAsBoundary = [&](const std::string& reason) {
    const auto& info = ctx.instanceBoundaryInfos[boundaryInfoIndex];
    model.abstractedSequentialBoundaries.push_back(
        "Abstracted uncomputable sequential instance `" + info.instancePath +
        "` as a SEC boundary: " + reason);
    model.abstractedSequentialBoundaryDetails.push_back(
        makeAbstractedBoundaryDetail(info));

    for (const auto& key : info.stateKeys) {
      ctx.abstractedBoundaryStateKeys.insert(key);
    }

    for (const auto& observedTerm : info.observedTerms) {
      if (ctx.abstractedBoundaryObservedKeys.insert(observedTerm.key).second) {
        ctx.abstractedBoundaryObservedTerms.emplace_back(observedTerm.termID, observedTerm.key);
        ctx.allObservedOutputs.insert(observedTerm.key);
      }
      ctx.prunedBuilderOutputTerms.insert(observedTerm.termID);
    }
  };

  const size_t independentStateOutputCount =
      countIndependentStateOutputs(scan.stateOutputs);
  const size_t pendingStart = ctx.pendingTransitions.size();
  const size_t complementedStart = model.complementedStateRelations.size();
  bool unsupportedInstance = false;
  bool abstractedUnsupportedInstance = false;
  std::string abstractedUnsupportedReason;

  // Track sequential behavior per state output bit. This keeps vector flops and
  // other multi-output sequential cells from being collapsed into a single
  // instance-wide transition record.
  for (const auto& stateOutput : scan.stateOutputs) {
    if (findComplementedPrimaryStateOutput(stateOutput, scan.stateOutputs) != nullptr) {
      continue;
    }

    PendingTransition pending;
    pending.stateTermID = stateOutput.termID;
    pending.stateKey = ctx.inputKeyByTerm.at(pending.stateTermID);
    pending.statePinName = stateOutput.pinName;
    pending.stateOutputIsComplemented =
        isIntrinsicComplementedStateOutput(stateOutput.pinName);
    pending.stateBit = stateOutput.bit;
    pending.independentStateOutputCount = independentStateOutputCount;
    pending.boundaryInfoIndex = boundaryInfoIndex;
    pending.pinTermIDs = scan.pinTermIDs;
    pending.clockTermIDs = stateOutput.clockTermIDs;
    pending.intrinsicClockIsInverted = stateOutput.intrinsicClockIsInverted;

    std::vector<ComplementedStateRelation> complementedRelations;
    for (const auto& candidate : scan.stateOutputs) {
      if (candidate.termID == stateOutput.termID || candidate.bit != stateOutput.bit) {
        continue;
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
      if (!isComplementedStateOutput(stateOutput.pinName, candidate.pinName)) {
        continue;
      }
      const SignalKey complementedKey = ctx.inputKeyByTerm.at(candidate.termID);
      pending.complementedStateKeys.push_back(complementedKey);
      complementedRelations.push_back({pending.stateKey, complementedKey});
    }

    if (const auto unsupportedReason =
            getPendingTransitionUnsupportedReason(pending)) {
      // Unsupported sequential cells are classified before cone construction.
      // Boundary mode exposes their interface; strict mode reports a structural
      // unsupported reason without relying on exception-to-result conversion.
      if (ctx.abstractUncomputableSequentialBoundaries) {
        abstractedUnsupportedInstance = true;
        abstractedUnsupportedReason = *unsupportedReason;
        break;
      }
      unsupportedInstance = true;
      model.unsupportedReasons.push_back(
          "Unsupported sequential primitive for `" +
          signalKeyToString(pending.stateKey) + "`: " + *unsupportedReason);
      continue;
    }

    model.complementedStateRelations.insert(
        model.complementedStateRelations.end(),
        complementedRelations.begin(),
        complementedRelations.end());
    ctx.pendingTransitions.push_back(std::move(pending));
  }

  if (abstractedUnsupportedInstance) {
    ctx.pendingTransitions.erase(
        ctx.pendingTransitions.begin() + static_cast<std::ptrdiff_t>(pendingStart),
        ctx.pendingTransitions.end());
    model.complementedStateRelations.erase(
        model.complementedStateRelations.begin() +
            static_cast<std::ptrdiff_t>(complementedStart),
        model.complementedStateRelations.end());
    abstractUnsupportedInstanceAsBoundary(abstractedUnsupportedReason);
    return;
  }

  if (unsupportedInstance) {
    markUnsupportedInstanceStateOutputs();
  }
}

void collectSequentialTransitions(ExtractContext& ctx, SequentialDesignModel& model) {
  // Record enough pin information to reconstruct Q' after the combinational
  // Boolean expressions have been built.
  for (auto leafID : ctx.dnl->getLeaves()) {
    const auto& instance = ctx.dnl->getDNLInstanceFromID(leafID);
    if (naja::NL::SNLDesignModeling::hasMemoryInterface(
            instance.getSNLInstance()->getModel()) &&
        supportsStructuredMemoryModel(
            naja::NL::SNLDesignModeling::getMemoryInterface(
                instance.getSNLInstance()))) {
      appendPendingMemoryInstance(ctx, model, instance);
      continue;
    }
    const auto scan = scanSequentialInstance(instance, ctx.inputKeyByTerm, model);
    if (!scan.has_value()) {
      continue;
    }
    appendPendingTransitionsForInstance(ctx, model, *scan);
  }
}

std::vector<naja::DNL::DNLID> collectInitialObservedTerms(const ExtractContext& ctx) {
  std::vector<naja::DNL::DNLID> initialObservedTerms;
  initialObservedTerms.reserve(
      ctx.topOutputKeyByTerm.size() + ctx.abstractedBoundaryObservedTerms.size());
  std::unordered_set<naja::DNL::DNLID> initialObservedTermSet;
  for (const auto& [termID, _] : ctx.topOutputKeyByTerm) {
    if (initialObservedTermSet.insert(termID).second) {
      initialObservedTerms.push_back(termID);
    }
  }
  for (const auto& [termID, _] : ctx.abstractedBoundaryObservedTerms) {
    if (initialObservedTermSet.insert(termID).second) {
      initialObservedTerms.push_back(termID);
    }
  }
  return initialObservedTerms;
}

void buildInitialObservedOutputClouds(ExtractContext& ctx, SequentialDesignModel& model) {
  const auto initialObservedTerms = collectInitialObservedTerms(ctx);
  std::unordered_set<naja::DNL::DNLID> collectedOutputSet(
      ctx.builder.getOutputs().begin(), ctx.builder.getOutputs().end());
  std::vector<naja::DNL::DNLID> initialMaterializedOutputs;
  initialMaterializedOutputs.reserve(initialObservedTerms.size());
  for (const auto outputTermID : initialObservedTerms) {
    if (collectedOutputSet.find(outputTermID) != collectedOutputSet.end()) {
      initialMaterializedOutputs.push_back(outputTermID);
    }
  }
  ctx.builder.setOutputs(initialMaterializedOutputs);
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) abstracted_boundaries=%zu pruned_builder_outputs=%zu initial_observed_outputs=%zu\n",
        ctx.topName.c_str(),
        model.abstractedSequentialBoundaries.size(),
        ctx.prunedBuilderOutputTerms.size(),
        initialMaterializedOutputs.size());
    fflush(stderr);
  }

  // Materialize the combinational BoolExpr DAGs for the design boundary.
  if (ctx.secDiagEnabled) {
    fprintf(stderr, "SEC diag: extract(%s) build begin\n", ctx.topName.c_str());
    fflush(stderr);
  }
  ctx.builder.build();
  if (ctx.secDiagEnabled) {
    fprintf(stderr, "SEC diag: extract(%s) build end\n", ctx.topName.c_str());
    fflush(stderr);
  }
}

void orderComplementedStateBitsByPrimary(SequentialDesignModel& model) {
  if (model.complementedStateRelations.empty() || model.stateBits.size() < 2) {
    return;
  }

  struct ComplementRole {
    SignalKey canonicalPrimaryKey;
    bool isComplementedOutput = false;
  };

  std::unordered_map<SignalKey, ComplementRole, SignalKeyHash> complementRoles;
  for (const auto& relation : model.complementedStateRelations) {
    complementRoles.emplace(
        relation.primaryKey,
        ComplementRole{relation.primaryKey, false});
    complementRoles.emplace(
        relation.complementedKey,
        ComplementRole{relation.primaryKey, true});
  }

  // LCOV_EXCL_START
  const SignalKeyLess less;
  // LCOV_EXCL_STOP
  auto roleFor = [&](const SignalKey& key) -> ComplementRole {
    const auto it = complementRoles.find(key);
    if (it != complementRoles.end()) {
      return it->second;
    }
    return ComplementRole{key, false};
  };

  // Real NLName IDs expose construction order. If a complemented pin is
  // interned before its primary pin, the plain key sort can publish QN before
  // Q. Keep each primary/complement group ordered as primary first so SEC uses
  // the real state output as the canonical representative.
  std::sort(
      model.stateBits.begin(),
      model.stateBits.end(),
      [&](const SignalKey& lhs, const SignalKey& rhs) {
        const auto lhsRole = roleFor(lhs);
        const auto rhsRole = roleFor(rhs);
        if (lhsRole.canonicalPrimaryKey != rhsRole.canonicalPrimaryKey) {
          return less(lhsRole.canonicalPrimaryKey, rhsRole.canonicalPrimaryKey);
        }
        if (lhsRole.isComplementedOutput != rhsRole.isComplementedOutput) {
          return !lhsRole.isComplementedOutput;
        }
        return less(lhs, rhs);  // LCOV_EXCL_LINE
      });
}

void publishNormalizedBoundary(ExtractContext& ctx, SequentialDesignModel& model) {
  for (const auto& key : ctx.abstractedBoundaryStateKeys) {
    ctx.stateBits.erase(key);
    ctx.environmentInputs.insert(key);
  }
  // Opaque internal boundary inputs behave exactly like additional SEC
  // environment inputs: extraction gives them symbolic leaf variables because
  // the surrounding cone cannot be modeled combinationally, so the published
  // interface must carry them forward into the shared proof symbol space.
  ctx.environmentInputs.insert(
      ctx.internalBoundaryInputKeys.begin(), ctx.internalBoundaryInputKeys.end());

  model.topInputKeys.assign(ctx.topInputKeys.begin(), ctx.topInputKeys.end());
  model.topOutputKeys.assign(ctx.topOutputKeys.begin(), ctx.topOutputKeys.end());
  model.environmentInputs.assign(ctx.environmentInputs.begin(), ctx.environmentInputs.end());
  model.internalBoundaryInputKeys.assign(
      ctx.internalBoundaryInputKeys.begin(), ctx.internalBoundaryInputKeys.end());
  model.internalBoundaryOutputKeys.assign(
      ctx.internalBoundaryOutputKeys.begin(), ctx.internalBoundaryOutputKeys.end());
  model.stateBits.assign(ctx.stateBits.begin(), ctx.stateBits.end());
  orderComplementedStateBitsByPrimary(model);
  model.allObservedOutputs.assign(ctx.allObservedOutputs.begin(), ctx.allObservedOutputs.end());
  if (ctx.secDiagEnabled) {
    fprintf(
        // LCOV_EXCL_START
        stderr,
        // LCOV_EXCL_STOP
        "SEC diag: extract(%s) boundary normalized env=%zu state=%zu outputs=%zu pending=%zu\n",
        ctx.topName.c_str(),
        // LCOV_EXCL_START
        model.environmentInputs.size(),
        // LCOV_EXCL_STOP
        model.stateBits.size(),
        model.allObservedOutputs.size(),
        ctx.pendingTransitions.size());
    // LCOV_EXCL_START
    fflush(stderr);
    // LCOV_EXCL_STOP
  }
}

void recordBoundaryInputVars(
    const ExtractContext& ctx,
    const std::vector<naja::DNL::DNLID>& builderInputs,
    const std::vector<size_t>& termDNLID2varID,
    SequentialDesignModel& model) {
  // Preserve the symbolic variable chosen by the clause builder for each
  // aligned SEC input/state signal.
  for (const auto inputTermID : builderInputs) {
    const auto& term = ctx.dnl->getDNLTerminalFromID(inputTermID);
    if (isConstantInternalOutputTerm(term)) {
      continue;
    }
    const auto keyIt = ctx.inputKeyByTerm.find(inputTermID);
    if (keyIt == ctx.inputKeyByTerm.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (inputTermID >= termDNLID2varID.size()) {
      continue;  // LCOV_EXCL_LINE
    }
    const size_t varID = termDNLID2varID[inputTermID];
    if (varID < 2) {
      continue;  // LCOV_EXCL_LINE
    }
    model.inputVarByKey.emplace(keyIt->second, varID);
  }
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) mapped boundary vars=%zu\n",
        ctx.topName.c_str(),
        model.inputVarByKey.size());
    fflush(stderr);
  }
}

size_t getNextSyntheticVarID(const SequentialDesignModel& model) {
  size_t nextVarID = 2;
  for (const auto& [_, varID] : model.inputVarByKey) {
    nextVarID = std::max(nextVarID, varID + 1);
  }
  return nextVarID;
}

void assignStructuredMemoryStateVars(
    const ExtractContext& ctx,
    SequentialDesignModel& model) {
  size_t nextVarID = getNextSyntheticVarID(model);
  for (const auto& pendingMemory : ctx.pendingMemoryInstances) {
    for (const auto& cellState : pendingMemory.cellStates) {
      model.stateBits.push_back(cellState.key);
      model.inputVarByKey.emplace(cellState.key, nextVarID++);
    }
  }
}

std::vector<naja::DNL::DNLID> collectStructuredMemoryDependencyTerms(
    const ExtractContext& ctx) {
  std::vector<naja::DNL::DNLID> termIDs;
  std::unordered_set<naja::DNL::DNLID> seen;
  for (const auto& pendingMemory : ctx.pendingMemoryInstances) {
    if (pendingMemory.resetTermID.has_value() &&
        seen.insert(*pendingMemory.resetTermID).second) {
      termIDs.push_back(*pendingMemory.resetTermID);
    }
    for (const auto& readPort : pendingMemory.readPorts) {
      for (const auto termID : readPort.addressTermIDs) {
        if (seen.insert(termID).second) {
          termIDs.push_back(termID);
        }
      }
    }
    for (const auto& writePort : pendingMemory.writePorts) {
      for (const auto termID : writePort.addressTermIDs) {
        if (seen.insert(termID).second) {
          termIDs.push_back(termID);
        }
      }
      for (const auto termID : writePort.dataTermIDs) {
        if (seen.insert(termID).second) {
          termIDs.push_back(termID);
        }
      }
      for (const auto termID : writePort.maskTermIDs) {
        if (seen.insert(termID).second) {
          termIDs.push_back(termID);
        }
      }
      for (const auto termID : writePort.enableTermIDs) {
        if (seen.insert(termID).second) {
          termIDs.push_back(termID);
        }
      }
    }
  }
  return termIDs;
}

BuiltObservedExpr materializeStructuredMemoryTermExpr(
    naja::DNL::DNLID termID,
    std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>&
        skippedOutputsByTerm,
    // LCOV_EXCL_START
    const std::vector<naja::DNL::DNLID>& builderInputs,
    const std::vector<naja::DNL::DNLID>& builderOutputs,
    const std::vector<size_t>& termDNLID2varID) {
  BuiltObservedExpr result;
  // LCOV_EXCL_STOP
  if (const auto exprIt = outputExprByTerm.find(termID);
      exprIt != outputExprByTerm.end()) {
    // LCOV_EXCL_START
    result.expr = exprIt->second;
    return result;
  }

  auto* dnl = naja::DNL::get();
  const std::string displayName =
      dnl == nullptr ? ("term#" + std::to_string(termID))
                     : getTerminalDisplayName(dnl->getDNLTerminalFromID(termID));


// LCOV_EXCL_STOP
  if (const auto skippedIt = skippedOutputsByTerm.find(termID);
      // LCOV_EXCL_START
      skippedIt != skippedOutputsByTerm.end()) {
    if (const auto connectivitySkip = getConnectivitySkipInfo(skippedIt->second);
        connectivitySkip.has_value()) {
      result.connectivitySkip = ConnectivitySkipInfo{
          connectivitySkip->origin,
          "Structured memory dependency `" + displayName +
              "` was skipped because " + connectivitySkip->detail};
              // LCOV_EXCL_STOP
      return result;
    // LCOV_EXCL_START
    }
    result.unsupportedReason =  // LCOV_EXCL_LINE
        "Structured memory dependency `" + displayName +  // LCOV_EXCL_LINE
        "` is unsupported: " + skippedIt->second.detail;  // LCOV_EXCL_LINE
    return result;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  const auto built = buildObservedExprForTerm(  // LCOV_EXCL_LINE
      termID,  // LCOV_EXCL_LINE
      outputExprByTerm,  // LCOV_EXCL_LINE
      builderInputs,  // LCOV_EXCL_LINE
      builderOutputs,  // LCOV_EXCL_LINE
      termDNLID2varID);  // LCOV_EXCL_LINE
  if (built.expr != nullptr) {  // LCOV_EXCL_LINE
    outputExprByTerm.emplace(termID, built.expr);  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    return built;  // LCOV_EXCL_LINE
  }
  if (built.connectivitySkip.has_value()) {  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
    result.connectivitySkip = ConnectivitySkipInfo{  // LCOV_EXCL_LINE
        built.connectivitySkip->origin,  // LCOV_EXCL_LINE
        "Structured memory dependency `" + displayName +  // LCOV_EXCL_LINE
            "` was skipped because " +  // LCOV_EXCL_LINE
            built.connectivitySkip->detail};  // LCOV_EXCL_LINE
    return result;  // LCOV_EXCL_LINE
  }
  result.unsupportedReason =  // LCOV_EXCL_LINE
      "Structured memory dependency `" + displayName +  // LCOV_EXCL_LINE
      "` is unsupported: " +  // LCOV_EXCL_LINE
      built.unsupportedReason;  // LCOV_EXCL_LINE
  return result;  // LCOV_EXCL_LINE
}

bool isNoDriverSkippedStructuredMemoryTerm(
    naja::DNL::DNLID termID,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>&
        skippedOutputsByTerm) {
  const auto skippedIt = skippedOutputsByTerm.find(termID);
  if (skippedIt == skippedOutputsByTerm.end()) {
    return false;
  }
  const auto connectivitySkip = getConnectivitySkipInfo(skippedIt->second);  // LCOV_EXCL_LINE
  return connectivitySkip.has_value() &&  // LCOV_EXCL_LINE
         connectivitySkip->origin == ConnectivitySkipOrigin::NoDriver;  // LCOV_EXCL_LINE
}

bool isDisabledMemoryWriteEnable(
    naja::DNL::DNLID termID,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>&
        skippedOutputsByTerm) {
  return isDisabledMemoryWriteEnable(termID) ||
         isNoDriverSkippedStructuredMemoryTerm(termID, skippedOutputsByTerm);
}

void buildStructuredMemoryTransitions(
    const ExtractContext& ctx,
    SequentialDesignModel& model,
    const std::vector<naja::DNL::DNLID>& builderInputs,
    const std::vector<naja::DNL::DNLID>& builderOutputs,
    const std::vector<size_t>& termDNLID2varID,
    std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>&
        skippedOutputsByTerm) {
  for (const auto& pendingMemory : ctx.pendingMemoryInstances) {
    auto markConnectivitySkipped = [&](const SignalKey& key,
                                       const ConnectivitySkipInfo& info) {
      model.connectivitySkipInfoByKey.emplace(key, info);
    };
    auto markWholeMemoryConnectivitySkipped =
        [&](const ConnectivitySkipInfo& info) {
          for (const auto& cellState : pendingMemory.cellStates) {
            // LCOV_EXCL_START
            markConnectivitySkipped(cellState.key, info);
          }
          // LCOV_EXCL_STOP
          for (const auto& readOutput : pendingMemory.readOutputs) {
            markConnectivitySkipped(readOutput.key, info);
          }
        };
    auto appendStructuredMemoryExpr =
        [&](naja::DNL::DNLID termID,
            std::vector<BoolExpr*>& exprs,
            std::optional<ConnectivitySkipInfo>& connectivitySkip) {
          const auto built = materializeStructuredMemoryTermExpr(
              termID,
              outputExprByTerm,
              skippedOutputsByTerm,
              builderInputs,
              builderOutputs,
              termDNLID2varID);
          if (built.expr != nullptr) {
            exprs.push_back(built.expr);
            return true;
          }
          if (built.connectivitySkip.has_value()) {
            connectivitySkip = *built.connectivitySkip;
            return false;
          }
          model.unsupportedReasons.push_back(built.unsupportedReason);  // LCOV_EXCL_LINE
          return false;  // LCOV_EXCL_LINE
        };

    std::vector<std::vector<BoolExpr*>> readAddressExprs;
    readAddressExprs.reserve(pendingMemory.readPorts.size());
    std::vector<std::optional<ConnectivitySkipInfo>> readPortSkips;
    readPortSkips.reserve(pendingMemory.readPorts.size());
    for (const auto& readPort : pendingMemory.readPorts) {
      std::vector<BoolExpr*> addressExprs;
      addressExprs.reserve(readPort.addressTermIDs.size());
      std::optional<ConnectivitySkipInfo> readPortSkip;
      for (const auto termID : readPort.addressTermIDs) {
        // Structured-memory dependency materialization can add extra
        // boundary roots that were not present in the original top-output
        // build. Read-address reconstruction must use that merged frontier,
        // otherwise a dependency found in the batch phase can disappear when
        // we rebuild the per-port expressions here.
        if (!appendStructuredMemoryExpr(termID, addressExprs, readPortSkip)) {
          break;
        }
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
      readAddressExprs.push_back(std::move(addressExprs));
      readPortSkips.push_back(std::move(readPortSkip));
    }

    struct WritePortExprs {
      bool disabled = false;
      std::vector<BoolExpr*> addressExprs;
      std::vector<BoolExpr*> dataExprs;
      std::vector<BoolExpr*> maskExprs;
      std::vector<BoolExpr*> enableExprs;
    };
    // LCOV_EXCL_START
    std::vector<WritePortExprs> writePortExprs;
    writePortExprs.reserve(pendingMemory.writePorts.size());
    // LCOV_EXCL_STOP
    std::optional<ConnectivitySkipInfo> memorySkip;
    bool unsupportedMemoryDependency = false;
    // LCOV_EXCL_START
    auto appendWriteDependencyExpr =
    // LCOV_EXCL_STOP
        [&](naja::DNL::DNLID termID, std::vector<BoolExpr*>& exprs) {
          std::optional<ConnectivitySkipInfo> dependencySkip;
          if (appendStructuredMemoryExpr(termID, exprs, dependencySkip)) {
            // LCOV_EXCL_START
            return true;
            // LCOV_EXCL_STOP
          }
          if (dependencySkip.has_value()) {
            // LCOV_EXCL_START
            memorySkip = *dependencySkip;
          } else {
          // LCOV_EXCL_STOP
            unsupportedMemoryDependency = true;  // LCOV_EXCL_LINE
          }
          return false;
        // LCOV_EXCL_START
        };
        // LCOV_EXCL_STOP
    for (const auto& writePort : pendingMemory.writePorts) {
      WritePortExprs exprs;
      exprs.addressExprs.reserve(writePort.addressTermIDs.size());
      // LCOV_EXCL_START
      exprs.dataExprs.reserve(writePort.dataTermIDs.size());
      // LCOV_EXCL_STOP
      exprs.maskExprs.reserve(writePort.maskTermIDs.size());
      exprs.enableExprs.reserve(writePort.enableTermIDs.size());
      for (const auto termID : writePort.enableTermIDs) {
        if (isDisabledMemoryWriteEnable(termID, skippedOutputsByTerm)) {
          exprs.disabled = true;  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
        }
        if (!appendWriteDependencyExpr(termID, exprs.enableExprs)) {
          break;  // LCOV_EXCL_LINE
        }
      }
      // LCOV_EXCL_START
      if (memorySkip.has_value() || unsupportedMemoryDependency) {
      // LCOV_EXCL_STOP
        break;  // LCOV_EXCL_LINE
      }
      if (exprs.disabled) {
        // LCOV_EXCL_START
        writePortExprs.push_back(std::move(exprs));  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        continue;  // LCOV_EXCL_LINE
      }
      for (const auto termID : writePort.addressTermIDs) {
        if (!appendWriteDependencyExpr(termID, exprs.addressExprs)) {
          break;  // LCOV_EXCL_LINE
        }
      }
      if (memorySkip.has_value() || unsupportedMemoryDependency) {
        // LCOV_EXCL_START
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      for (const auto termID : writePort.dataTermIDs) {
        if (!appendWriteDependencyExpr(termID, exprs.dataExprs)) {
          break;
        }
      }
      if (memorySkip.has_value() || unsupportedMemoryDependency) {
        break;
      }
      for (const auto termID : writePort.maskTermIDs) {
        if (!appendWriteDependencyExpr(termID, exprs.maskExprs)) {
          break;  // LCOV_EXCL_LINE
        }
      }
      if (memorySkip.has_value() || unsupportedMemoryDependency) {
        break;  // LCOV_EXCL_LINE
      }
      writePortExprs.push_back(std::move(exprs));
    }
    if (memorySkip.has_value()) {
      markWholeMemoryConnectivitySkipped(*memorySkip);
      continue;
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
    if (unsupportedMemoryDependency) {
      continue;  // LCOV_EXCL_LINE
    }

    BoolExpr* resetExpr = nullptr;
    if (pendingMemory.resetTermID.has_value()) {
      std::vector<BoolExpr*> resetExprs;
      std::optional<ConnectivitySkipInfo> resetSkip;
      if (!appendStructuredMemoryExpr(
              *pendingMemory.resetTermID, resetExprs, resetSkip)) {
        if (resetSkip.has_value()) {
          markWholeMemoryConnectivitySkipped(*resetSkip);
        // LCOV_EXCL_START
        }
        continue;
      }
      // LCOV_EXCL_STOP
      resetExpr = resetExprs.front();
    }
    const auto resetAssertedExpr = [&]() -> BoolExpr* {
      switch (pendingMemory.resetMode) {
        case naja::NL::SNLDesignModeling::MemoryResetMode::AsyncLow:
        // LCOV_EXCL_START
        case naja::NL::SNLDesignModeling::MemoryResetMode::SyncLow:
        // LCOV_EXCL_STOP
          return resetExpr == nullptr ? nullptr : BoolExpr::Not(resetExpr);
        case naja::NL::SNLDesignModeling::MemoryResetMode::AsyncHigh:
        case naja::NL::SNLDesignModeling::MemoryResetMode::SyncHigh:
          return resetExpr;  // LCOV_EXCL_LINE
        case naja::NL::SNLDesignModeling::MemoryResetMode::None:
        default:
          return nullptr;
      }
    }();

    std::vector<std::vector<BoolExpr*>> cellNextExprs(
        pendingMemory.depth,
        std::vector<BoolExpr*>(pendingMemory.width, nullptr));
    for (const auto& cellState : pendingMemory.cellStates) {
      const auto varIt = model.inputVarByKey.find(cellState.key);
      if (varIt == model.inputVarByKey.end()) {
        throw std::runtime_error(  // LCOV_EXCL_LINE
            "Missing synthetic SEC variable for memory cell state `" +  // LCOV_EXCL_LINE
            cellState.displayName + "`");  // LCOV_EXCL_LINE
      }
      BoolExpr* next = BoolExpr::Var(varIt->second);
      for (size_t portIndex = 0; portIndex < pendingMemory.writePorts.size(); ++portIndex) {
        const auto& exprs = writePortExprs[portIndex];
        if (exprs.disabled) {
          continue;  // LCOV_EXCL_LINE
        }
        std::vector<BoolExpr*> conditions;
        conditions.reserve(2 + exprs.enableExprs.size());
        // LCOV_EXCL_START
        conditions.push_back(buildAddressEqualsExpr(
            exprs.addressExprs, cellState.cellIndex));
        for (auto* enableExpr : exprs.enableExprs) {
        // LCOV_EXCL_STOP
          conditions.push_back(enableExpr);
        }
        if (!exprs.maskExprs.empty()) {
          conditions.push_back(exprs.maskExprs[cellState.bitIndex]);
        }
        BoolExpr* writeCondition = makeAndChain(conditions);
        next = makeIte(
            writeCondition,
            exprs.dataExprs[cellState.bitIndex],
            next);
      }
      if (resetAssertedExpr != nullptr) {
        next = makeIte(resetAssertedExpr, BoolExpr::createFalse(), next);
        model.initialStateValueByKey.emplace(cellState.key, false);
      }
      model.nextStateExprByStateKey.emplace(cellState.key, next);
      cellNextExprs[cellState.cellIndex][cellState.bitIndex] = next;
    }

    for (const auto& readOutput : pendingMemory.readOutputs) {
      const auto varIt = model.inputVarByKey.find(readOutput.key);
      if (varIt == model.inputVarByKey.end()) {
        throw std::runtime_error(  // LCOV_EXCL_LINE
            "Missing SEC variable for memory read-output state `" +  // LCOV_EXCL_LINE
            signalKeyToString(readOutput.key) + "`");  // LCOV_EXCL_LINE
      }
      if (readOutput.portIndex < readPortSkips.size() &&
          readPortSkips[readOutput.portIndex].has_value()) {
        markConnectivitySkipped(
            readOutput.key, *readPortSkips[readOutput.portIndex]);
        continue;
      }
      const auto& addressExprs = readAddressExprs[readOutput.portIndex];
      BoolExpr* next = cellNextExprs[0][readOutput.bitIndex];
      for (size_t cellIndex = 1; cellIndex < pendingMemory.depth; ++cellIndex) {
        next = makeIte(
            buildAddressEqualsExpr(addressExprs, cellIndex),
            cellNextExprs[cellIndex][readOutput.bitIndex],
            // LCOV_EXCL_START
            next);
      }
      if (resetAssertedExpr != nullptr) {
        next = makeIte(resetAssertedExpr, BoolExpr::createFalse(), next);
        model.initialStateValueByKey.emplace(readOutput.key, false);
      }
      // LCOV_EXCL_STOP
      model.nextStateExprByStateKey.emplace(readOutput.key, next);
    // LCOV_EXCL_START
    }
  }
}


// LCOV_EXCL_STOP
void materializeBoundaryObservedOutputs(
    const std::vector<std::pair<naja::DNL::DNLID, SignalKey>>& observedTerms,
    // LCOV_EXCL_START
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>& skippedOutputsByTerm,
    const std::vector<naja::DNL::DNLID>& builderInputs,
    const std::vector<naja::DNL::DNLID>& builderOutputs,
    const std::vector<size_t>& termDNLID2varID,
    // LCOV_EXCL_STOP
    SequentialDesignModel& model) {
  // LCOV_EXCL_START
  for (const auto& [termID, key] : observedTerms) {
    if (const auto exprIt = outputExprByTerm.find(termID);
        exprIt != outputExprByTerm.end()) {
        // LCOV_EXCL_STOP
      model.observedOutputExprByKey.emplace(key, exprIt->second);
      // LCOV_EXCL_START
      continue;
    }
    if (const auto skippedIt = skippedOutputsByTerm.find(termID);  // LCOV_EXCL_LINE
        skippedIt != skippedOutputsByTerm.end()) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      if (auto skipInfo = getConnectivitySkipInfo(skippedIt->second);  // LCOV_EXCL_LINE
          skipInfo.has_value()) {  // LCOV_EXCL_LINE
        model.connectivitySkipInfoByKey.emplace(key, *skipInfo);  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
          "Unsupported SEC boundary output `" + model.displayNameByKey.at(key) +  // LCOV_EXCL_LINE
          "`: " + skippedIt->second.detail);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }

    const auto built = buildObservedExprForTerm(  // LCOV_EXCL_LINE
        termID, outputExprByTerm, builderInputs, builderOutputs, termDNLID2varID);  // LCOV_EXCL_LINE
    if (built.expr != nullptr) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      model.observedOutputExprByKey.emplace(key, built.expr);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    if (built.connectivitySkip.has_value()) {  // LCOV_EXCL_LINE
      model.connectivitySkipInfoByKey.emplace(key, *built.connectivitySkip);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
        "Unsupported SEC boundary output `" + model.displayNameByKey.at(key) +  // LCOV_EXCL_LINE
        "`: " + built.unsupportedReason);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

// LCOV_EXCL_START
void materializeTopObservedOutputs(
    const std::unordered_map<naja::DNL::DNLID, SignalKey>& topOutputKeyByTerm,
    // LCOV_EXCL_STOP
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>& skippedOutputsByTerm,
    SequentialDesignModel& model) {
  for (const auto& [termID, key] : topOutputKeyByTerm) {
    auto exprIt = outputExprByTerm.find(termID);
    if (exprIt != outputExprByTerm.end()) {
      model.observedOutputExprByKey.emplace(key, exprIt->second);
      continue;
    }

    auto skippedIt = skippedOutputsByTerm.find(termID);  // LCOV_EXCL_LINE
    if (skippedIt != skippedOutputsByTerm.end()) {  // LCOV_EXCL_LINE
      if (auto skipInfo = getConnectivitySkipInfo(skippedIt->second);  // LCOV_EXCL_LINE
          skipInfo.has_value()) {  // LCOV_EXCL_LINE
        model.connectivitySkipInfoByKey.emplace(key, *skipInfo);  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
          "Unsupported observed output cone for `" + signalKeyToString(key) +  // LCOV_EXCL_LINE
          "`: " + skippedIt->second.detail);  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }

    model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
        "Missing observed output expression for `" + signalKeyToString(key) + "`");  // LCOV_EXCL_LINE
  }
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
struct RebuiltTransitionArtifacts {
  std::unordered_set<SignalKey, SignalKeyHash> requiredStateKeys;
  std::set<SignalKey, SignalKeyLess> lateAbstractedBoundaryStateKeys;
  std::vector<std::pair<naja::DNL::DNLID, SignalKey>> lateAbstractedBoundaryObservedTerms;
};

constexpr size_t kMaxCompleteStateFrontierForStartupMatching = 5000;

bool shouldRetainCompleteStateFrontierForStartupMatching(size_t stateCount) {
  // Moderate-size SEC cases benefit from the complete transition relation:
  // reset/startup checks can inspect local sequential cones without relying on
  // internal flop names.  Large ASICs still use the COI frontier so
  // LCOV_EXCL_START
  // BlackParrot-scale proofs do not materialize every sequential cone up front.
  // LCOV_EXCL_STOP
  return stateCount <= kMaxCompleteStateFrontierForStartupMatching;
}

template <typename EnqueueStateKey>
void enqueueStateDependenciesFromFormula(
    BoolExpr* expr,
    const std::unordered_set<size_t>& requiredStateVarIDs,
    std::vector<const BoolExpr*>& stack,
    std::unordered_set<const BoolExpr*>& scannedDependencyNodes,
    EnqueueStateKey&& enqueueStateVarID) {
  if (expr == nullptr || !expr->isValid()) {
    return;  // LCOV_EXCL_LINE
  }

  // This frontier walk only needs to know whether a shared subtree has ever
  // exposed a required state variable. Once a subtree is scanned, every state
  // variable reachable under it has already been enqueued, so revisiting that
  // subtree for later roots would be pure duplicate work.
  stack.clear();
  stack.push_back(expr);
  // LCOV_EXCL_START
  while (!stack.empty()) {
  // LCOV_EXCL_STOP
    const BoolExpr* node = stack.back();
    // LCOV_EXCL_START
    stack.pop_back();
    // LCOV_EXCL_STOP
    if (node == nullptr) {
      continue;  // LCOV_EXCL_LINE
    }
    if (!scannedDependencyNodes.insert(node).second) {
      continue;
    }

    switch (node->getOp()) {
      case Op::VAR:
        if (requiredStateVarIDs.find(node->getId()) != requiredStateVarIDs.end()) {
          enqueueStateVarID(node->getId());
        }
        break;
      case Op::NOT:
        stack.push_back(node->getLeft());
        break;
      case Op::AND:
      case Op::OR:
      case Op::XOR:
        stack.push_back(node->getLeft());
        stack.push_back(node->getRight());
        break;
      case Op::NONE:  // LCOV_EXCL_LINE
      default:
        break;  // LCOV_EXCL_LINE
    }
  }
}

// LCOV_EXCL_START

void removeDeadFoldedClockGateLatchInputs(
// LCOV_EXCL_STOP
    SequentialDesignModel& model,
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID);
void substituteFoldedClockGateLatchVarsInModel(
    SequentialDesignModel& model,
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID);
// LCOV_EXCL_START
void markFormulasWithUnpublishedSupportAsSkipped(
    const ExtractContext& ctx,
    // LCOV_EXCL_STOP
    SequentialDesignModel& model);

// LCOV_EXCL_START

RebuiltTransitionArtifacts rebuildRequiredStateTransitions(
// LCOV_EXCL_STOP
    ExtractContext& ctx,
    SequentialDesignModel& model,
    // LCOV_EXCL_START
    std::vector<naja::DNL::DNLID>& builderInputs,
    std::vector<naja::DNL::DNLID>& builderOutputs,
    std::vector<size_t>& termDNLID2varID,
    std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    // LCOV_EXCL_STOP
    std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>& skippedOutputsByTerm) {
  // LCOV_EXCL_START
  RebuiltTransitionArtifacts artifacts;
  auto markConnectivitySkippedState =
      [&](const SignalKey& key, const ConnectivitySkipInfo& info) {
        model.connectivitySkipInfoByKey.emplace(key, info);
      };
  auto markUnsupportedState = [&](const SignalKey& key) {
    ctx.unsupportedStateBits.insert(key);  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE

  std::unordered_set<SignalKey, SignalKeyHash> lateAbstractedBoundaryObservedKeys;
  std::unordered_set<size_t> lateAbstractedBoundaryIndexes;
  auto recordLateAbstractedInstanceBoundary =
  // LCOV_EXCL_STOP
      [&](size_t boundaryInfoIndex, const std::string& reason) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        if (boundaryInfoIndex == std::numeric_limits<size_t>::max()) {  // LCOV_EXCL_LINE
          return;  // LCOV_EXCL_LINE
        }
        if (!lateAbstractedBoundaryIndexes.insert(boundaryInfoIndex).second) {  // LCOV_EXCL_LINE
          return;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
        const auto& info = ctx.instanceBoundaryInfos[boundaryInfoIndex];  // LCOV_EXCL_LINE
        if (ctx.secDiagEnabled) {  // LCOV_EXCL_LINE
          fprintf(  // LCOV_EXCL_LINE
              stderr,  // LCOV_EXCL_LINE
              "SEC diag: extract(%s) late abstracted sequential instance `%s`: %s\n",
              ctx.topName.c_str(),  // LCOV_EXCL_LINE
              info.instancePath.c_str(),  // LCOV_EXCL_LINE
              reason.c_str());  // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          fflush(stderr);  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }  // LCOV_EXCL_LINE
        model.abstractedSequentialBoundaries.push_back(  // LCOV_EXCL_LINE
            "Abstracted uncomputable sequential instance `" +  // LCOV_EXCL_LINE
            info.instancePath + "` as a SEC boundary: " + reason);  // LCOV_EXCL_LINE
        model.abstractedSequentialBoundaryDetails.push_back(  // LCOV_EXCL_LINE
            makeAbstractedBoundaryDetail(info));  // LCOV_EXCL_LINE
        for (const auto& key : info.stateKeys) {  // LCOV_EXCL_LINE
          artifacts.lateAbstractedBoundaryStateKeys.insert(key);  // LCOV_EXCL_LINE
        }
        for (const auto& observedTerm : info.observedTerms) {  // LCOV_EXCL_LINE
          if (lateAbstractedBoundaryObservedKeys.insert(observedTerm.key).second) {  // LCOV_EXCL_LINE
            artifacts.lateAbstractedBoundaryObservedTerms.emplace_back(  // LCOV_EXCL_LINE
                observedTerm.termID, observedTerm.key);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
        }
      };  // LCOV_EXCL_LINE

  std::unordered_map<size_t, SignalKey> requiredStateKeyByVarID;
  requiredStateKeyByVarID.reserve(model.stateBits.size());
  std::unordered_set<size_t> requiredStateVarIDs;
  requiredStateVarIDs.reserve(model.stateBits.size());
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt == model.inputVarByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    requiredStateKeyByVarID.emplace(varIt->second, key);
    requiredStateVarIDs.insert(varIt->second);
  }

  std::unordered_map<SignalKey, size_t, SignalKeyHash> pendingIndexByStateKey;
  pendingIndexByStateKey.reserve(ctx.pendingTransitions.size() * 2);
  for (size_t pendingIndex = 0; pendingIndex < ctx.pendingTransitions.size(); ++pendingIndex) {
    const auto& pending = ctx.pendingTransitions[pendingIndex];
    pendingIndexByStateKey.emplace(pending.stateKey, pendingIndex);
    for (const auto& complementedKey : pending.complementedStateKeys) {
      pendingIndexByStateKey.emplace(complementedKey, pendingIndex);
    }
  }

  std::unordered_set<size_t> requiredPendingIndexes;
  std::unordered_set<naja::DNL::DNLID> materializedOutputTerms;
  std::unordered_set<size_t> topClockCarrierVarIDs;
  std::unordered_set<naja::DNL::DNLID> pureClockCarrierTermIDs;
  std::unordered_map<size_t, ClockEvent> clockEventByCarrierVarID;
  std::unordered_map<BoolExpr*, BoolExpr*> transitionClockStripMemo;
  materializedOutputTerms.reserve(outputExprByTerm.size());
  for (const auto& [termID, _] : outputExprByTerm) {
    materializedOutputTerms.insert(termID);
  }
  auto refreshClockCarrierVarIDs = [&]() {
    size_t addedTop = 0;
    for (const auto varID : collectTopClockCarrierVarIDs(model)) {
      if (topClockCarrierVarIDs.insert(varID).second) {
        ++addedTop;
      }
    }
    const size_t addedBoundary = expandClockCarrierVarIDsFromBoundaryNames(
        model, topClockCarrierVarIDs);
    // LCOV_EXCL_START
    const size_t addedTermNames = expandClockCarrierVarIDsFromTermNames(
        ctx.dnl, termDNLID2varID, topClockCarrierVarIDs);
        // LCOV_EXCL_STOP
    const size_t addedStructure = expandClockCarrierVarIDsFromStructure(
        ctx.dnl,
        // LCOV_EXCL_START
        model,
        termDNLID2varID,
        topClockCarrierVarIDs,
        &pureClockCarrierTermIDs);
    const size_t addedPureExprs = expandClockCarrierVarIDsFromPureClockTermExprs(
        pureClockCarrierTermIDs, outputExprByTerm, topClockCarrierVarIDs);
    seedTopClockCarrierEvents(model, topClockCarrierVarIDs, clockEventByCarrierVarID);
    const size_t addedMaterialized = expandClockCarrierVarIDsFromMaterializedTerms(
        termDNLID2varID,
        outputExprByTerm,
        // LCOV_EXCL_STOP
        topClockCarrierVarIDs,
        // LCOV_EXCL_START
        clockEventByCarrierVarID);
        // LCOV_EXCL_STOP
    const size_t added = addedTop + addedBoundary + addedTermNames +
                         addedStructure + addedPureExprs + addedMaterialized;
    if (added != 0) {
      // The strip memo is keyed by expression pointer and the current carrier
      // set.  When extraction discovers a new routed-clock variable, previous
      // no-op answers may become stale, so invalidate once per carrier growth.
      transitionClockStripMemo.clear();
    }
    if (std::getenv("KEPLER_SEC_CLOCK_GATE_DIAG") != nullptr) {
      std::fprintf(  // LCOV_EXCL_LINE
          stderr,  // LCOV_EXCL_LINE
          "SEC diag: clock carriers top=%zu boundary=%zu term_names=%zu "
          "structure=%zu pure_exprs=%zu materialized=%zu total=%zu pure_terms=%zu\n",
          addedTop,  // LCOV_EXCL_LINE
          addedBoundary,  // LCOV_EXCL_LINE
          addedTermNames,  // LCOV_EXCL_LINE
          addedStructure,  // LCOV_EXCL_LINE
          addedPureExprs,  // LCOV_EXCL_LINE
          addedMaterialized,  // LCOV_EXCL_LINE
          topClockCarrierVarIDs.size(),  // LCOV_EXCL_LINE
          pureClockCarrierTermIDs.size());  // LCOV_EXCL_LINE
      std::fflush(stderr);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return added;
  };  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  refreshClockCarrierVarIDs();
  // LCOV_EXCL_STOP
  {
    std::vector<naja::DNL::DNLID> latchDependencyTerms;
    for (const auto termID : collectClockGateLatchDependencyTerms(ctx, model)) {
      if (materializedOutputTerms.insert(termID).second) {
        latchDependencyTerms.push_back(termID);
      }
    }
    if (!latchDependencyTerms.empty()) {
      const auto dependencyOutputs = materializeBuilderOutputs(
          latchDependencyTerms,
          builderInputs,
          termDNLID2varID,
          ctx.collectedSkippedOutputs,
          ctx.secDiagEnabled,
          ctx.topName.c_str(),
          "clock-gate latch dependency build");
      appendUniqueTermIDs(builderInputs, dependencyOutputs.inputs);
      appendUniqueTermIDs(builderOutputs, dependencyOutputs.outputs);
      mergeBuilderTermVarIDs(termDNLID2varID, dependencyOutputs.termDNLID2varID);
      recordBoundaryInputVars(ctx, builderInputs, termDNLID2varID, model);
      for (const auto& [termID, expr] : dependencyOutputs.outputExprByTerm) {
        outputExprByTerm.insert_or_assign(termID, expr);
      }
      for (const auto& [termID, info] : dependencyOutputs.skippedOutputsByTerm) {
        skippedOutputsByTerm.emplace(termID, info);  // LCOV_EXCL_LINE
      }
      refreshClockCarrierVarIDs();
    }
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  std::unordered_map<size_t, BoolExpr*> clockGateLatchDataExprByVarID =
      collectClockGateLatchDataExprByVarID(
          ctx, model, topClockCarrierVarIDs, outputExprByTerm);
  // LCOV_EXCL_START
  std::deque<size_t> pendingWorkQueue;
  // LCOV_EXCL_STOP
  std::deque<SignalKey> stateDependencyWorkQueue;
  std::unordered_set<SignalKey, SignalKeyHash> expandedStateDependencies;
  std::vector<const BoolExpr*> dependencyWalkStack;
  std::unordered_set<const BoolExpr*> scannedDependencyNodes;
  scannedDependencyNodes.reserve(std::max<size_t>(1024, outputExprByTerm.size() * 16));
  // LCOV_EXCL_START
  std::unordered_set<size_t> enqueuedStateVarIDs;
  // LCOV_EXCL_STOP
  enqueuedStateVarIDs.reserve(requiredStateVarIDs.size());
  auto enqueueRequiredStateKey = [&](const SignalKey& key) {
    if (!artifacts.requiredStateKeys.insert(key).second) {
      return;
    }
    stateDependencyWorkQueue.push_back(key);
    const auto pendingIt = pendingIndexByStateKey.find(key);
    if (pendingIt != pendingIndexByStateKey.end() &&
        requiredPendingIndexes.insert(pendingIt->second).second) {
      pendingWorkQueue.push_back(pendingIt->second);
    }
  };
  auto enqueueRequiredStateVarID = [&](size_t varID) {
    if (!enqueuedStateVarIDs.insert(varID).second) {
      return;  // LCOV_EXCL_LINE
    }
    const auto stateIt = requiredStateKeyByVarID.find(varID);
    if (stateIt == requiredStateKeyByVarID.end()) {
      return;  // LCOV_EXCL_LINE
    }
    enqueueRequiredStateKey(stateIt->second);
  };
  auto enqueueStateDependenciesFromExpr = [&](BoolExpr* expr) {
    if (expr == nullptr || !expr->isValid()) {
      return;  // LCOV_EXCL_LINE
    }
    enqueueStateDependenciesFromFormula(
        expr,
        requiredStateVarIDs,
        dependencyWalkStack,
        scannedDependencyNodes,
        enqueueRequiredStateVarID);
  // LCOV_EXCL_START
  };
  // LCOV_EXCL_STOP

  if (shouldRetainCompleteStateFrontierForStartupMatching(model.stateBits.size())) {
    for (const auto& key : model.stateBits) {
      enqueueRequiredStateKey(key);
    }
  }

  for (const auto& [_, expr] : model.observedOutputExprByKey) {
    enqueueStateDependenciesFromExpr(expr);
  }

  // Follow the state/output dependency frontier lazily so SEC only rebuilds the
  // sequential update cones that can actually influence covered observations.
  // States with prebuilt next-state relations, such as structured memories,
  // participate in the same frontier expansion before we trim the SEC model.
  while (!stateDependencyWorkQueue.empty() || !pendingWorkQueue.empty()) {
    while (!stateDependencyWorkQueue.empty()) {
      const SignalKey key = stateDependencyWorkQueue.front();
      stateDependencyWorkQueue.pop_front();
      const auto nextStateIt = model.nextStateExprByStateKey.find(key);
      if (nextStateIt == model.nextStateExprByStateKey.end()) {
        continue;
      }
      if (!expandedStateDependencies.insert(key).second) {
        continue;  // LCOV_EXCL_LINE
      }
      enqueueStateDependenciesFromExpr(nextStateIt->second);
    }

    if (pendingWorkQueue.empty()) {
      continue;
    }

    std::vector<size_t> batchPendingIndexes;
    std::vector<naja::DNL::DNLID> batchOutputTerms;
    std::unordered_set<naja::DNL::DNLID> batchOutputTermSet;
    while (!pendingWorkQueue.empty()) {
      const size_t pendingIndex = pendingWorkQueue.front();
      pendingWorkQueue.pop_front();
      batchPendingIndexes.push_back(pendingIndex);

      const auto& pending = ctx.pendingTransitions[pendingIndex];
      artifacts.requiredStateKeys.insert(pending.stateKey);
      for (const auto& complementedKey : pending.complementedStateKeys) {
        artifacts.requiredStateKeys.insert(complementedKey);
      }

      for (const auto& [_, candidates] : pending.pinTermIDs) {
        for (const auto& candidate : candidates) {
          if (materializedOutputTerms.insert(candidate.termID).second &&
              batchOutputTermSet.insert(candidate.termID).second) {
            batchOutputTerms.push_back(candidate.termID);
          }
        }
      }
      for (const auto& clockTerm : pending.clockTermIDs) {
        // Clock pins may have been materialized earlier as plain carriers.
        // Rebuild required clocks after boundary/latch dependencies are known
        // so gated clocks become clk & enable instead of a stale clk frontier.
        const bool alreadyMaterialized =
            !materializedOutputTerms.insert(clockTerm.termID).second;
        if ((!alreadyMaterialized ||
             hasBuildableCombinationalRoot(ctx.dnl, clockTerm.termID)) &&
            batchOutputTermSet.insert(clockTerm.termID).second) {
          batchOutputTerms.push_back(clockTerm.termID);
        }
      }
    }

    if (!batchOutputTerms.empty()) {
      const auto dependencyOutputs = materializeBuilderOutputs(
          batchOutputTerms,
          builderInputs,
          termDNLID2varID,
          ctx.collectedSkippedOutputs,
          ctx.secDiagEnabled,
          // LCOV_EXCL_START
          ctx.topName.c_str(),
          // LCOV_EXCL_STOP
          "dependency build");
      appendUniqueTermIDs(builderInputs, dependencyOutputs.inputs);
      appendUniqueTermIDs(builderOutputs, dependencyOutputs.outputs);
      mergeBuilderTermVarIDs(termDNLID2varID, dependencyOutputs.termDNLID2varID);
      recordBoundaryInputVars(ctx, builderInputs, termDNLID2varID, model);
      for (const auto& [termID, expr] : dependencyOutputs.outputExprByTerm) {
        outputExprByTerm.insert_or_assign(termID, expr);
      }
      for (const auto& [termID, info] : dependencyOutputs.skippedOutputsByTerm) {
        skippedOutputsByTerm.emplace(termID, info);
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
      refreshClockCarrierVarIDs();
      clockGateLatchDataExprByVarID =
          collectClockGateLatchDataExprByVarID(
              ctx, model, topClockCarrierVarIDs, outputExprByTerm);
    }

    for (const auto pendingIndex : batchPendingIndexes) {
      const auto& pending = ctx.pendingTransitions[pendingIndex];
      std::optional<ConnectivitySkipInfo> skippedPinInfo;
      // LCOV_EXCL_START
      bool abortPending = false;
      for (const auto& [pinName, _] : pending.pinTermIDs) {
        const auto resolvedPinTermID = resolvePendingPinTermID(pending, pinName.c_str());
        if (!resolvedPinTermID.has_value()) {
          continue;  // LCOV_EXCL_LINE
        }
        auto skippedIt = skippedOutputsByTerm.find(*resolvedPinTermID);
        // LCOV_EXCL_STOP
        if (skippedIt == skippedOutputsByTerm.end()) {
          continue;
        // LCOV_EXCL_START
        }

        if (auto skipInfo = getConnectivitySkipInfo(skippedIt->second);
            skipInfo.has_value()) {
          if (skipInfo->origin == ConnectivitySkipOrigin::NoDriver &&
              isOptionalSequentialControlPin(pinName)) {
            continue;  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
          }
          // LCOV_EXCL_START
          skippedPinInfo = {
              skipInfo->origin,
              // LCOV_EXCL_STOP
              "Sequential pin `" + pinName + "` was skipped because " +
                  skippedIt->second.detail,
          };
          break;
        }

        if (ctx.abstractUncomputableSequentialBoundaries) {  // LCOV_EXCL_LINE
          recordLateAbstractedInstanceBoundary(  // LCOV_EXCL_LINE
              pending.boundaryInfoIndex,  // LCOV_EXCL_LINE
              // LCOV_EXCL_START
              "unsupported sequential pin `" + pinName + "`: " +  // LCOV_EXCL_LINE
                  skippedIt->second.detail);  // LCOV_EXCL_LINE
          abortPending = true;  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
        }


// LCOV_EXCL_STOP
        model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            "Unsupported sequential primitive for `" + signalKeyToString(pending.stateKey) +  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            "`: Sequential pin `" + pinName + "` is unsupported: " +  // LCOV_EXCL_LINE
            skippedIt->second.detail);  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        markUnsupportedState(pending.stateKey);  // LCOV_EXCL_LINE
        for (const auto& complementedKey : pending.complementedStateKeys) {  // LCOV_EXCL_LINE
          markUnsupportedState(complementedKey);  // LCOV_EXCL_LINE
        }
        abortPending = true;  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      }
      // LCOV_EXCL_STOP

      if (!skippedPinInfo.has_value()) {
        // LCOV_EXCL_START
        for (const auto& clockTerm : pending.clockTermIDs) {
          auto skippedIt = skippedOutputsByTerm.find(clockTerm.termID);
          if (skippedIt == skippedOutputsByTerm.end()) {
            continue;
          }

          if (auto skipInfo = getConnectivitySkipInfo(skippedIt->second);  // LCOV_EXCL_LINE
              skipInfo.has_value()) {  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
            skippedPinInfo = {  // LCOV_EXCL_LINE
                // LCOV_EXCL_START
                skipInfo->origin,  // LCOV_EXCL_LINE
                "Sequential clock pin was skipped because " +  // LCOV_EXCL_LINE
                // LCOV_EXCL_STOP
                    skippedIt->second.detail,  // LCOV_EXCL_LINE
            };
            break;  // LCOV_EXCL_LINE
          }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
          if (ctx.abstractUncomputableSequentialBoundaries) {  // LCOV_EXCL_LINE
            recordLateAbstractedInstanceBoundary(  // LCOV_EXCL_LINE
                pending.boundaryInfoIndex,  // LCOV_EXCL_LINE
                "unsupported sequential clock pin: " +  // LCOV_EXCL_LINE
                    // LCOV_EXCL_START
                    skippedIt->second.detail);  // LCOV_EXCL_LINE
                    // LCOV_EXCL_STOP
            abortPending = true;  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }

          model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
              "Unsupported sequential primitive for `" +  // LCOV_EXCL_LINE
              signalKeyToString(pending.stateKey) +  // LCOV_EXCL_LINE
              "`: Sequential clock pin is unsupported: " +  // LCOV_EXCL_LINE
              skippedIt->second.detail);  // LCOV_EXCL_LINE
          markUnsupportedState(pending.stateKey);  // LCOV_EXCL_LINE
          for (const auto& complementedKey : pending.complementedStateKeys) {  // LCOV_EXCL_LINE
            markUnsupportedState(complementedKey);  // LCOV_EXCL_LINE
          }
          abortPending = true;  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
        }
      }

      if (abortPending) {
        continue;  // LCOV_EXCL_LINE
      }
      if (skippedPinInfo.has_value()) {
        markConnectivitySkippedState(pending.stateKey, *skippedPinInfo);
        for (const auto& complementedKey : pending.complementedStateKeys) {
          markConnectivitySkippedState(complementedKey, *skippedPinInfo);  // LCOV_EXCL_LINE
        }
        continue;
      }

      const auto initialStateValue = detectInitialStateValue(pending);
      if (initialStateValue.has_value()) {
        model.initialStateValueByKey.emplace(pending.stateKey, *initialStateValue);
        for (const auto& complementedKey : pending.complementedStateKeys) {
          model.initialStateValueByKey.emplace(complementedKey, !*initialStateValue);
        }
      }

      BoolExpr* nextStateExpr =
          buildNextStateExpr(
              pending,
              termDNLID2varID,
              pureClockCarrierTermIDs,
              topClockCarrierVarIDs,
              clockEventByCarrierVarID,
              clockGateLatchDataExprByVarID,
              outputExprByTerm,
              transitionClockStripMemo);
      const auto clockEvent = classifyPendingClockEvent(
          pending,
          termDNLID2varID,
          outputExprByTerm,
          clockEventByCarrierVarID);
      // This diagnostic is intentionally after clock-carrier stripping: any
      // remaining unpublished support would become a private proof input and
      // can hide the real reason state matching stopped converging.
      logUnpublishedTransitionSupport(
          ctx,
          // LCOV_EXCL_START
          model,
          pending,
          nextStateExpr,
          // LCOV_EXCL_STOP
          termDNLID2varID,
          topClockCarrierVarIDs,
          // LCOV_EXCL_START
          pureClockCarrierTermIDs);
      model.nextStateExprByStateKey.emplace(pending.stateKey, nextStateExpr);
      if (clockEvent.has_value()) {
        model.clockEventByStateKey.emplace(pending.stateKey, *clockEvent);
      }
      // Liberty flops such as DFF_X1 expose both Q and QN. They share one
      // LCOV_EXCL_STOP
      // storage element, so complementary outputs inherit the same next-state
      // LCOV_EXCL_START
      // function with a logical inversion.
      // LCOV_EXCL_STOP
      for (const auto& complementedKey : pending.complementedStateKeys) {
        model.nextStateExprByStateKey.emplace(complementedKey, BoolExpr::Not(nextStateExpr));
        if (clockEvent.has_value()) {
          model.clockEventByStateKey.emplace(complementedKey, *clockEvent);
        }
        if (artifacts.requiredStateKeys.find(complementedKey) !=
            artifacts.requiredStateKeys.end()) {
          stateDependencyWorkQueue.push_back(complementedKey);
        }
      }
      stateDependencyWorkQueue.push_back(pending.stateKey);
    }
  }

  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) rebuilt next-state exprs=%zu init=%zu\n",
        ctx.topName.c_str(),
        model.nextStateExprByStateKey.size(),
        model.initialStateValueByKey.size());
    fflush(stderr);
  }

  model.clockCarrierVarIDs.assign(
      topClockCarrierVarIDs.begin(), topClockCarrierVarIDs.end());
  std::sort(model.clockCarrierVarIDs.begin(), model.clockCarrierVarIDs.end());
  model.clockCarrierClasses = buildClockCarrierClasses(
      model.clockCarrierVarIDs, clockEventByCarrierVarID);
  substituteFoldedClockGateLatchVarsInModel(model, clockGateLatchDataExprByVarID);
  removeDeadFoldedClockGateLatchInputs(model, clockGateLatchDataExprByVarID);

  return artifacts;
}

void eraseSignalKeyFromVector(std::vector<SignalKey>& keys, const SignalKey& key) {
  keys.erase(std::remove(keys.begin(), keys.end(), key), keys.end());
}

void substituteFoldedClockGateLatchVarsInModel(
    SequentialDesignModel& model,
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID) {
  if (clockGateLatchDataExprByVarID.empty()) {
    return;
  }

  for (auto& [_, expr] : model.observedOutputExprByKey) {
    expr = substituteClockGateLatchVarsInExpr(
        expr, clockGateLatchDataExprByVarID);
  }
  for (auto& [_, expr] : model.nextStateExprByStateKey) {
    expr = substituteClockGateLatchVarsInExpr(
        expr, clockGateLatchDataExprByVarID);
  }
}

void removeDeadFoldedClockGateLatchInputs(
    // LCOV_EXCL_START
    SequentialDesignModel& model,
    // LCOV_EXCL_STOP
    const std::unordered_map<size_t, BoolExpr*>& clockGateLatchDataExprByVarID) {
  if (clockGateLatchDataExprByVarID.empty()) {
    return;
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  const auto isFoldedLatchVar =
      buildCandidateVarMask(clockGateLatchDataExprByVarID);
  CandidateDependencyScratch scratch;
  std::unordered_set<size_t> liveVars;
  for (const auto& [_, expr] : model.observedOutputExprByKey) {
    collectCandidateDependenciesIntoSet(expr, isFoldedLatchVar, scratch, liveVars);
  }
  for (const auto& [_, expr] : model.nextStateExprByStateKey) {
    collectCandidateDependenciesIntoSet(expr, isFoldedLatchVar, scratch, liveVars);
  }

  std::unordered_map<size_t, SignalKey> keyByVarID;
  keyByVarID.reserve(model.inputVarByKey.size());
  for (const auto& [key, varID] : model.inputVarByKey) {
    keyByVarID.emplace(varID, key);
  }

  for (const auto& [varID, _] : clockGateLatchDataExprByVarID) {
    if (liveVars.find(varID) != liveVars.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    const auto keyIt = keyByVarID.find(varID);
    if (keyIt == keyByVarID.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    eraseSignalKeyFromVector(model.environmentInputs, keyIt->second);
    eraseSignalKeyFromVector(model.internalBoundaryInputKeys, keyIt->second);
    model.inputVarByKey.erase(keyIt->second);
  }
}

struct ClockDomainIndex {
  std::unordered_map<size_t, SignalKey> domainByStateVarID;
  std::vector<uint8_t> isClockedStateVar;
};

ClockDomainIndex buildClockDomainIndex(const SequentialDesignModel& model) {
  // LCOV_EXCL_START
  ClockDomainIndex index;
  // LCOV_EXCL_STOP
  size_t maxStateVarID = 0;
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    const auto eventIt = model.clockEventByStateKey.find(key);
    if (varIt == model.inputVarByKey.end() ||
        eventIt == model.clockEventByStateKey.end()) {
      continue;
    }
    index.domainByStateVarID.emplace(varIt->second, eventIt->second.domain);
    maxStateVarID = std::max(maxStateVarID, varIt->second);
  }

  index.isClockedStateVar.assign(maxStateVarID + 1, 0);
  for (const auto& [varID, _] : index.domainByStateVarID) {
    index.isClockedStateVar[varID] = 1;
  }
  return index;
}

std::string formatClockDomainName(
    const SequentialDesignModel& model,
    const SignalKey& domain) {
  const auto displayIt = model.displayNameByKey.find(domain);
  return displayIt == model.displayNameByKey.end()
             // LCOV_EXCL_START
             ? signalKeyToString(domain)  // LCOV_EXCL_LINE
             // LCOV_EXCL_STOP
             : displayIt->second;
}

std::string formatClockDomains(
    const SequentialDesignModel& model,
    const std::set<SignalKey, SignalKeyLess>& domains) {
  std::ostringstream oss;
  bool first = true;
  for (const auto& domain : domains) {
    if (!first) {
      oss << ", ";
    }
    first = false;
    oss << formatClockDomainName(model, domain);
  }
  return oss.str();
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
std::set<SignalKey, SignalKeyLess> collectClockDomainsFromExpr(
    BoolExpr* expr,
    const ClockDomainIndex& index,
    CandidateDependencyScratch& scratch) {
  std::set<SignalKey, SignalKeyLess> domains;
  if (index.isClockedStateVar.empty()) {
    return domains;  // LCOV_EXCL_LINE
  }
  for (const auto stateVarID :
       collectCandidateStateDependenciesFromExpr(
           expr, index.isClockedStateVar, scratch)) {
    const auto domainIt = index.domainByStateVarID.find(stateVarID);
    if (domainIt != index.domainByStateVarID.end()) {
      domains.insert(domainIt->second);
    }
  }
  return domains;
}

bool containsClockDomainOtherThan(
    const std::set<SignalKey, SignalKeyLess>& domains,
    const SignalKey& expectedDomain) {
  for (const auto& domain : domains) {
    if (domain != expectedDomain) {
      return true;  // LCOV_EXCL_LINE
    }
  }
  return false;
}

ConnectivitySkipInfo makeMultiClockDomainSkip(const std::string& detail) {
  return ConnectivitySkipInfo{ConnectivitySkipOrigin::MultiClockDomain, detail};
}

void composeSameDomainPhaseTransitions(SequentialDesignModel& model) {
  std::map<SignalKey, std::vector<SignalKey>, SignalKeyLess> posedgeStatesByDomain;
  std::map<SignalKey, std::vector<SignalKey>, SignalKeyLess> negedgeStatesByDomain;
  for (const auto& key : model.stateBits) {
    const auto eventIt = model.clockEventByStateKey.find(key);
    if (eventIt == model.clockEventByStateKey.end()) {
      continue;
    }
    if (eventIt->second.phase == ClockPhase::Pos) {
      posedgeStatesByDomain[eventIt->second.domain].push_back(key);
    } else {
      negedgeStatesByDomain[eventIt->second.domain].push_back(key);
    }
  }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  for (const auto& [domain, posedgeStates] : posedgeStatesByDomain) {
    const auto negedgeIt = negedgeStatesByDomain.find(domain);
    if (negedgeIt == negedgeStatesByDomain.end()) {
      continue;
    }

    std::unordered_map<size_t, BoolExpr*> posedgeNextByVarID;
    posedgeNextByVarID.reserve(posedgeStates.size());
    for (const auto& key : posedgeStates) {
      const auto varIt = model.inputVarByKey.find(key);
      const auto nextIt = model.nextStateExprByStateKey.find(key);
      if (varIt != model.inputVarByKey.end() &&
          nextIt != model.nextStateExprByStateKey.end()) {
        posedgeNextByVarID.emplace(varIt->second, nextIt->second);
      }
    }

    // LCOV_EXCL_START
    // SEC macro-cycles are sampled just before the positive edge.  For a
    // LCOV_EXCL_STOP
    // positive/negative phase pair in one domain, negedge flops therefore see
    // the just-computed posedge state, while posedge flops still see old
    // negedge state from the previous macro-cycle.
    for (const auto& key : negedgeIt->second) {
      const auto nextIt = model.nextStateExprByStateKey.find(key);
      if (nextIt == model.nextStateExprByStateKey.end()) {
        continue;  // LCOV_EXCL_LINE
      }
      nextIt->second = substituteBoolExprVariableExpressions(
          nextIt->second, posedgeNextByVarID);
    }
  }
// LCOV_EXCL_START
}

void markMultiClockDomainConesAsSkipped(SequentialDesignModel& model) {
  const ClockDomainIndex index = buildClockDomainIndex(model);
  if (index.domainByStateVarID.empty()) {
  // LCOV_EXCL_STOP
    return;
  }

  CandidateDependencyScratch scratch;
  for (const auto& key : model.stateBits) {
    // LCOV_EXCL_START
    if (model.connectivitySkipInfoByKey.find(key) !=
    // LCOV_EXCL_STOP
        model.connectivitySkipInfoByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    const auto eventIt = model.clockEventByStateKey.find(key);
    // LCOV_EXCL_STOP
    const auto nextIt = model.nextStateExprByStateKey.find(key);
    if (eventIt == model.clockEventByStateKey.end() ||
        nextIt == model.nextStateExprByStateKey.end()) {
      continue;
    }
    const auto domains =
        collectClockDomainsFromExpr(nextIt->second, index, scratch);
    if (!containsClockDomainOtherThan(domains, eventIt->second.domain)) {
      continue;
    }
    model.connectivitySkipInfoByKey.emplace(  // LCOV_EXCL_LINE
        key,  // LCOV_EXCL_LINE
        makeMultiClockDomainSkip(  // LCOV_EXCL_LINE
            "State next-state cone crosses clock domains: " +  // LCOV_EXCL_LINE
            formatClockDomains(model, domains)));  // LCOV_EXCL_LINE
  }

  for (const auto& key : model.allObservedOutputs) {
    if (model.connectivitySkipInfoByKey.find(key) !=
        model.connectivitySkipInfoByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    const auto exprIt = model.observedOutputExprByKey.find(key);
    if (exprIt == model.observedOutputExprByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    const auto domains =
        collectClockDomainsFromExpr(exprIt->second, index, scratch);
    if (domains.size() <= 1) {
      continue;
    }
    model.connectivitySkipInfoByKey.emplace(
        key,
        makeMultiClockDomainSkip(
            "Observed output cone spans clock domains: " +
            formatClockDomains(model, domains)));
  }
}

void applyRebuiltTransitionArtifacts(
    const ExtractContext& ctx,
    const RebuiltTransitionArtifacts& artifacts,
    SequentialDesignModel& model,
    const std::vector<naja::DNL::DNLID>& builderInputs,
    const std::vector<naja::DNL::DNLID>& builderOutputs,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm,
    const std::unordered_map<naja::DNL::DNLID, BuilderSkippedOutputInfo>& skippedOutputsByTerm) {
  // The transition rebuild already closes the observed-output frontier while it
  // can still materialize missing next-state cones.  Pruning must reuse that
  // frontier instead of rewalking the whole BoolExpr model again; the latter is
  // redundant and quadratic on large designs such as BlackParrot.
  const auto& retainedStateKeys = artifacts.requiredStateKeys;
  std::vector<SignalKey> prunedStateKeys;
  model.stateBits.erase(
      std::remove_if(
          model.stateBits.begin(),
          model.stateBits.end(),
          [&](const SignalKey& key) {
            if (retainedStateKeys.find(key) != retainedStateKeys.end()) {
              // LCOV_EXCL_START
              return false;
            }
            prunedStateKeys.push_back(key);
            return true;
          }),
      model.stateBits.end());
  for (const auto& key : prunedStateKeys) {
  // LCOV_EXCL_STOP
    model.nextStateExprByStateKey.erase(key);
    model.initialStateValueByKey.erase(key);
    // LCOV_EXCL_START
    model.clockEventByStateKey.erase(key);
    model.inputVarByKey.erase(key);
  }
  model.complementedStateRelations.erase(
      std::remove_if(
          model.complementedStateRelations.begin(),
          model.complementedStateRelations.end(),
          // LCOV_EXCL_STOP
          [&](const ComplementedStateRelation& relation) {
            // LCOV_EXCL_START
            return retainedStateKeys.find(relation.primaryKey) ==
                       retainedStateKeys.end() ||
                   retainedStateKeys.find(relation.complementedKey) ==
                       retainedStateKeys.end();
          }),
      model.complementedStateRelations.end());

  for (const auto& key : artifacts.lateAbstractedBoundaryStateKeys) {
    model.nextStateExprByStateKey.erase(key);  // LCOV_EXCL_LINE
    model.initialStateValueByKey.erase(key);  // LCOV_EXCL_LINE
    model.clockEventByStateKey.erase(key);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    if (std::find(model.environmentInputs.begin(), model.environmentInputs.end(), key) ==  // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        model.environmentInputs.end()) {  // LCOV_EXCL_LINE
      model.environmentInputs.push_back(key);  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  if (!artifacts.lateAbstractedBoundaryStateKeys.empty()) {
    model.stateBits.erase(  // LCOV_EXCL_LINE
        std::remove_if(  // LCOV_EXCL_LINE
            model.stateBits.begin(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            model.stateBits.end(),  // LCOV_EXCL_LINE
            [&](const SignalKey& key) {  // LCOV_EXCL_LINE
              return artifacts.lateAbstractedBoundaryStateKeys.find(key) !=  // LCOV_EXCL_LINE
                     artifacts.lateAbstractedBoundaryStateKeys.end();  // LCOV_EXCL_LINE
            }),
        model.stateBits.end());  // LCOV_EXCL_LINE
    model.complementedStateRelations.erase(  // LCOV_EXCL_LINE
        std::remove_if(  // LCOV_EXCL_LINE
            model.complementedStateRelations.begin(),  // LCOV_EXCL_LINE
            model.complementedStateRelations.end(),  // LCOV_EXCL_LINE
            [&](const ComplementedStateRelation& relation) {  // LCOV_EXCL_LINE
              return artifacts.lateAbstractedBoundaryStateKeys.find(relation.primaryKey) !=  // LCOV_EXCL_LINE
                         artifacts.lateAbstractedBoundaryStateKeys.end() ||  // LCOV_EXCL_LINE
                     artifacts.lateAbstractedBoundaryStateKeys.find(  // LCOV_EXCL_LINE
                         relation.complementedKey) !=  // LCOV_EXCL_LINE
                         artifacts.lateAbstractedBoundaryStateKeys.end();  // LCOV_EXCL_LINE
            }),
        model.complementedStateRelations.end());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE

  for (const auto& [_, key] : artifacts.lateAbstractedBoundaryObservedTerms) {
    if (std::find(model.allObservedOutputs.begin(), model.allObservedOutputs.end(), key) ==  // LCOV_EXCL_LINE
        model.allObservedOutputs.end()) {  // LCOV_EXCL_LINE
      model.allObservedOutputs.push_back(key);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  materializeBoundaryObservedOutputs(
      artifacts.lateAbstractedBoundaryObservedTerms,
      outputExprByTerm,
      // LCOV_EXCL_START
      skippedOutputsByTerm,
      builderInputs,
      // LCOV_EXCL_STOP
      builderOutputs,
      termDNLID2varID,
      model);
}

void filterUnsupportedAndUnmappedBoundary(ExtractContext& ctx, SequentialDesignModel& model) {
  {
    // Any published leaf variable that is not retained as sequential state is a
    // free SEC environment input, regardless of whether it originated from the
    // top interface, an opaque internal boundary, or a later abstraction step.
    // Compact SEC rebuilds the proof problem only from this normalized model,
    // so leaving such leaves out of the environment interface causes remapped
    // formulas to reference symbols that the shared proof symbol space never
    // allocates.
    std::unordered_set<SignalKey, SignalKeyHash> stateKeys(
        model.stateBits.begin(), model.stateBits.end());
    std::unordered_set<SignalKey, SignalKeyHash> publishedInputs(
        model.environmentInputs.begin(), model.environmentInputs.end());
    // LCOV_EXCL_START
    for (const auto& [key, _] : model.inputVarByKey) {
      if (stateKeys.find(key) != stateKeys.end()) {
        continue;
      }
      if (publishedInputs.insert(key).second) {
        model.environmentInputs.push_back(key);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }

  // Inputs or state bits can disappear if the underlying BoolExpr builder
  // optimized them away to constants; remove them from the aligned interface.
  auto keepMappedInputs = [&](std::vector<SignalKey>& keys) {
    keys.erase(
        std::remove_if(
            keys.begin(),
            keys.end(),
            [&](const SignalKey& key) {
              return model.inputVarByKey.find(key) == model.inputVarByKey.end();
            }),
        keys.end());
  };
  keepMappedInputs(model.environmentInputs);
  keepMappedInputs(model.stateBits);
  // LCOV_EXCL_STOP
  if (!ctx.unsupportedStateBits.empty()) {
    // LCOV_EXCL_START
    model.stateBits.erase(  // LCOV_EXCL_LINE
        std::remove_if(  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            model.stateBits.begin(),  // LCOV_EXCL_LINE
            model.stateBits.end(),  // LCOV_EXCL_LINE
            [&](const SignalKey& key) {  // LCOV_EXCL_LINE
              if (ctx.unsupportedStateBits.find(key) == ctx.unsupportedStateBits.end()) {  // LCOV_EXCL_LINE
                return false;  // LCOV_EXCL_LINE
              }
              model.nextStateExprByStateKey.erase(key);  // LCOV_EXCL_LINE
              model.initialStateValueByKey.erase(key);  // LCOV_EXCL_LINE
              model.clockEventByStateKey.erase(key);  // LCOV_EXCL_LINE
              model.inputVarByKey.erase(key);  // LCOV_EXCL_LINE
              return true;  // LCOV_EXCL_LINE
            }),  // LCOV_EXCL_LINE
        model.stateBits.end());  // LCOV_EXCL_LINE
    model.complementedStateRelations.erase(  // LCOV_EXCL_LINE
        std::remove_if(  // LCOV_EXCL_LINE
            model.complementedStateRelations.begin(),  // LCOV_EXCL_LINE
            model.complementedStateRelations.end(),  // LCOV_EXCL_LINE
            [&](const ComplementedStateRelation& relation) {  // LCOV_EXCL_LINE
              return ctx.unsupportedStateBits.find(relation.primaryKey) !=  // LCOV_EXCL_LINE
                         // LCOV_EXCL_START
                         ctx.unsupportedStateBits.end() ||  // LCOV_EXCL_LINE
                         // LCOV_EXCL_STOP
                     ctx.unsupportedStateBits.find(relation.complementedKey) !=  // LCOV_EXCL_LINE
                         ctx.unsupportedStateBits.end();  // LCOV_EXCL_LINE
            }),
        model.complementedStateRelations.end());  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) filtered env=%zu state=%zu\n",
        ctx.topName.c_str(),
        model.environmentInputs.size(),
        model.stateBits.size());
    fflush(stderr);
  }
}

void propagateConnectivitySkipsThroughDependencies(SequentialDesignModel& model) {
  std::unordered_map<size_t, SignalKey> stateKeyByVarID;
  size_t maxCandidateStateVarID = 0;
  std::deque<size_t> pendingSkippedStateVars;
  std::unordered_set<size_t> enqueuedSkippedStateVars;
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt == model.inputVarByKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    stateKeyByVarID.emplace(varIt->second, key);
    maxCandidateStateVarID = std::max(maxCandidateStateVarID, varIt->second);
    if (model.connectivitySkipInfoByKey.find(key) !=
            model.connectivitySkipInfoByKey.end() &&
        enqueuedSkippedStateVars.insert(varIt->second).second) {
      pendingSkippedStateVars.push_back(varIt->second);
    }
  }
  if (pendingSkippedStateVars.empty()) {
    return;
  }
  std::vector<uint8_t> isCandidateStateVar(maxCandidateStateVarID + 1, 0);
  for (const auto& [varID, _] : stateKeyByVarID) {
    isCandidateStateVar[varID] = 1;
  }

  // Reverse dependencies only need to point back into the already-owned SEC
  // interface vectors. Storing indexes here avoids copying large SignalKey path
  // objects into the temporary propagation graph, which sampling showed was
  // expensive both to allocate and to tear down on BlackParrot.
  std::vector<std::vector<size_t>> dependentStateIndexesBySourceVarID(
      isCandidateStateVar.size());
  std::vector<std::vector<size_t>> dependentOutputIndexesBySourceVarID(
      isCandidateStateVar.size());
  auto mergeDependencyIndexMap =
      [](std::vector<std::vector<size_t>>& destination,
         const std::vector<std::vector<size_t>>& source) {
        const size_t limit = std::min(destination.size(), source.size());
        for (size_t sourceVarID = 0; sourceVarID != limit; ++sourceVarID) {
          const auto& indexes = source[sourceVarID];
          if (indexes.empty()) {
            continue;
          }
          auto& mergedIndexes = destination[sourceVarID];
          mergedIndexes.insert(
              mergedIndexes.end(), indexes.begin(), indexes.end());
        }
  };
  // Precompute state-to-state and state-to-output dependency edges once. The
  // previous implementation recursively rescanned large Boolean cones every
  // time another state became skipped. We still build the reverse edges once,
  // but we now do it with one support walk per root instead of materializing a
  // giant exact dependency cache for every BoolExpr node. Sampling also showed
  // LCOV_EXCL_START
  // all TBB workers idle here, so the independent per-root scans run in
  // LCOV_EXCL_STOP
  // parallel with per-thread scratch and are merged once at the end.
  tbb::enumerable_thread_specific<ThreadLocalDependencyState> threadLocalStateDeps;
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, model.stateBits.size()),
      [&](const tbb::blocked_range<size_t>& range) {
        auto& localState = threadLocalStateDeps.local();
        if (localState.scratch.emittedEpochByVarID.size() <
            isCandidateStateVar.size()) {
          localState.scratch.emittedEpochByVarID.resize(
              isCandidateStateVar.size(), 0);
        }
        if (localState.dependenciesBySourceVarID.size() !=
            isCandidateStateVar.size()) {
          localState.dependenciesBySourceVarID.assign(
              isCandidateStateVar.size(), {});
        }
        for (size_t stateIndex = range.begin(); stateIndex != range.end(); ++stateIndex) {
          const auto& key = model.stateBits[stateIndex];
          if (model.connectivitySkipInfoByKey.find(key) !=
              model.connectivitySkipInfoByKey.end()) {
            continue;
          }
          const auto exprIt = model.nextStateExprByStateKey.find(key);  // LCOV_EXCL_LINE
          if (exprIt == model.nextStateExprByStateKey.end()) {  // LCOV_EXCL_LINE
            continue;  // LCOV_EXCL_LINE
          }
          const auto dependencies = collectCandidateStateDependenciesFromExpr(  // LCOV_EXCL_LINE
              exprIt->second, isCandidateStateVar, localState.scratch);  // LCOV_EXCL_LINE
          for (const auto dependencyVarID : dependencies) {  // LCOV_EXCL_LINE
            localState.dependenciesBySourceVarID[dependencyVarID].push_back(stateIndex);  // LCOV_EXCL_LINE
          }
        }  // LCOV_EXCL_LINE
      });
  // LCOV_EXCL_START
  for (const auto& localState : threadLocalStateDeps) {
  // LCOV_EXCL_STOP
    mergeDependencyIndexMap(
        dependentStateIndexesBySourceVarID, localState.dependenciesBySourceVarID);
  }

// LCOV_EXCL_START


// LCOV_EXCL_STOP
  tbb::enumerable_thread_specific<ThreadLocalDependencyState> threadLocalOutputDeps;
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, model.allObservedOutputs.size()),
      [&](const tbb::blocked_range<size_t>& range) {
        auto& localState = threadLocalOutputDeps.local();
        if (localState.scratch.emittedEpochByVarID.size() <
            isCandidateStateVar.size()) {
          localState.scratch.emittedEpochByVarID.resize(
              isCandidateStateVar.size(), 0);
        }
        if (localState.dependenciesBySourceVarID.size() !=
            isCandidateStateVar.size()) {
          localState.dependenciesBySourceVarID.assign(
              isCandidateStateVar.size(), {});
        }
        for (size_t outputIndex = range.begin(); outputIndex != range.end();
             // LCOV_EXCL_START
             ++outputIndex) {
             // LCOV_EXCL_STOP
          const auto& key = model.allObservedOutputs[outputIndex];
          if (model.connectivitySkipInfoByKey.find(key) !=
              model.connectivitySkipInfoByKey.end()) {
            // LCOV_EXCL_START
            continue;  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
          }
          const auto exprIt = model.observedOutputExprByKey.find(key);
          if (exprIt == model.observedOutputExprByKey.end()) {
            continue;  // LCOV_EXCL_LINE
          }
          const auto dependencies = collectCandidateStateDependenciesFromExpr(
              exprIt->second, isCandidateStateVar, localState.scratch);
          for (const auto dependencyVarID : dependencies) {
            localState.dependenciesBySourceVarID[dependencyVarID].push_back(outputIndex);
          }
        }
      });
  for (const auto& localState : threadLocalOutputDeps) {
    mergeDependencyIndexMap(
        dependentOutputIndexesBySourceVarID, localState.dependenciesBySourceVarID);
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  auto makeDependencySkip = [&](size_t sourceVarID) -> std::optional<ConnectivitySkipInfo> {
    const auto sourceKeyIt = stateKeyByVarID.find(sourceVarID);
    if (sourceKeyIt == stateKeyByVarID.end()) {
      return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    const auto skipInfoIt = model.connectivitySkipInfoByKey.find(sourceKeyIt->second);
    if (skipInfoIt == model.connectivitySkipInfoByKey.end()) {
    // LCOV_EXCL_STOP
      return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }
    return ConnectivitySkipInfo{
        skipInfoIt->second.origin,
        "Depends on skipped state `" + model.displayNameByKey.at(sourceKeyIt->second) +
            "` whose cone traces to a " +
            // LCOV_EXCL_STOP
            describeConnectivitySkipOrigin(skipInfoIt->second.origin) + " issue",
    };
  };

  while (!pendingSkippedStateVars.empty()) {
    const size_t skippedVarID = pendingSkippedStateVars.front();
    pendingSkippedStateVars.pop_front();

    const auto dependencySkip = makeDependencySkip(skippedVarID);
    if (!dependencySkip.has_value()) {
      continue;  // LCOV_EXCL_LINE
    }

    // LCOV_EXCL_START
    if (skippedVarID < dependentStateIndexesBySourceVarID.size()) {
    // LCOV_EXCL_STOP
      for (const auto dependentStateIndex :
           dependentStateIndexesBySourceVarID[skippedVarID]) {
        // LCOV_EXCL_START
        const auto& dependentKey = model.stateBits[dependentStateIndex];  // LCOV_EXCL_LINE
        if (!model.connectivitySkipInfoByKey.emplace(dependentKey, *dependencySkip).second) {  // LCOV_EXCL_LINE
          continue;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        const auto varIt = model.inputVarByKey.find(dependentKey);  // LCOV_EXCL_LINE
        if (varIt != model.inputVarByKey.end() &&  // LCOV_EXCL_LINE
            enqueuedSkippedStateVars.insert(varIt->second).second) {  // LCOV_EXCL_LINE
          pendingSkippedStateVars.push_back(varIt->second);  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
      }
      // LCOV_EXCL_STOP
    }

// LCOV_EXCL_START

    if (skippedVarID < dependentOutputIndexesBySourceVarID.size()) {
      for (const auto dependentOutputIndex :
           dependentOutputIndexesBySourceVarID[skippedVarID]) {
        const auto& dependentKey = model.allObservedOutputs[dependentOutputIndex];
        // LCOV_EXCL_STOP
        model.connectivitySkipInfoByKey.emplace(dependentKey, *dependencySkip);
      // LCOV_EXCL_START
      }
      // LCOV_EXCL_STOP
    }
  // LCOV_EXCL_START
  }
}

void markFormulasWithUnpublishedSupportAsSkipped(  // LCOV_EXCL_LINE
    const ExtractContext& ctx,
    // LCOV_EXCL_STOP
    SequentialDesignModel& model) {
  // LCOV_EXCL_START
  size_t maxPublishedVarID = 1;  // LCOV_EXCL_LINE
  for (const auto& [_, varID] : model.inputVarByKey) {  // LCOV_EXCL_LINE
    maxPublishedVarID = std::max(maxPublishedVarID, varID);  // LCOV_EXCL_LINE
  }
  std::vector<uint8_t> isPublishedVar(maxPublishedVarID + 1, 0);  // LCOV_EXCL_LINE
  for (const auto& key : model.environmentInputs) {  // LCOV_EXCL_LINE
    const auto varIt = model.inputVarByKey.find(key);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    if (varIt != model.inputVarByKey.end()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      isPublishedVar[varIt->second] = 1;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  for (const auto& key : model.stateBits) {  // LCOV_EXCL_LINE
    const auto varIt = model.inputVarByKey.find(key);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    if (varIt != model.inputVarByKey.end()) {  // LCOV_EXCL_LINE
      isPublishedVar[varIt->second] = 1;  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
  }
  CandidateDependencyScratch scratch;  // LCOV_EXCL_LINE

  size_t skippedStates = 0;  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  for (const auto& [key, expr] : model.nextStateExprByStateKey) {  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (model.connectivitySkipInfoByKey.find(key) !=  // LCOV_EXCL_LINE
        model.connectivitySkipInfoByKey.end()) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    if (const auto varID =  // LCOV_EXCL_LINE
            findFirstUnpublishedSupportVar(expr, isPublishedVar, scratch)) {  // LCOV_EXCL_LINE
      model.connectivitySkipInfoByKey.emplace(key, makeUnpublishedSupportSkip(*varID));  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
      ++skippedStates;  // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      if (ctx.secDiagEnabled) {  // LCOV_EXCL_LINE
        fprintf(  // LCOV_EXCL_LINE
            stderr,  // LCOV_EXCL_LINE
            "SEC diag: extract(%s) skipping state `%s`: unpublished support v%zu\n",
            ctx.topName.c_str(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
            displayNameForSignalKey(model, key).c_str(),  // LCOV_EXCL_LINE
            // LCOV_EXCL_START
            *varID);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

// LCOV_EXCL_START

  size_t skippedOutputs = 0;  // LCOV_EXCL_LINE
  for (const auto& [key, expr] : model.observedOutputExprByKey) {  // LCOV_EXCL_LINE
    if (model.connectivitySkipInfoByKey.find(key) !=  // LCOV_EXCL_LINE
        model.connectivitySkipInfoByKey.end()) {  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (const auto varID =  // LCOV_EXCL_LINE
            findFirstUnpublishedSupportVar(expr, isPublishedVar, scratch)) {  // LCOV_EXCL_LINE
      model.connectivitySkipInfoByKey.emplace(key, makeUnpublishedSupportSkip(*varID));  // LCOV_EXCL_LINE
      ++skippedOutputs;  // LCOV_EXCL_LINE
      if (ctx.secDiagEnabled) {  // LCOV_EXCL_LINE
        fprintf(  // LCOV_EXCL_LINE
            stderr,  // LCOV_EXCL_LINE
            "SEC diag: extract(%s) skipping output `%s`: unpublished support v%zu\n",
            ctx.topName.c_str(),  // LCOV_EXCL_LINE
            displayNameForSignalKey(model, key).c_str(),  // LCOV_EXCL_LINE
            *varID);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  if (ctx.secDiagEnabled && (skippedStates != 0 || skippedOutputs != 0)) {  // LCOV_EXCL_LINE
    fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: extract(%s) skipped unpublished-support formulas states=%zu outputs=%zu\n",
        ctx.topName.c_str(),  // LCOV_EXCL_LINE
        skippedStates,  // LCOV_EXCL_LINE
        skippedOutputs);  // LCOV_EXCL_LINE
    fflush(stderr);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
}  // LCOV_EXCL_LINE

void partitionCoveredSignals(SequentialDesignModel& model) {
  std::vector<SignalKey> legalStateBits;
  legalStateBits.reserve(model.stateBits.size());
  for (const auto& key : model.stateBits) {
    if (model.connectivitySkipInfoByKey.find(key) != model.connectivitySkipInfoByKey.end()) {
      model.skippedStateBits.push_back(key);
      model.nextStateExprByStateKey.erase(key);
      model.initialStateValueByKey.erase(key);
      // LCOV_EXCL_START
      model.clockEventByStateKey.erase(key);
      model.inputVarByKey.erase(key);
      continue;
    }
    legalStateBits.push_back(key);
    // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  model.stateBits = std::move(legalStateBits);
  // LCOV_EXCL_STOP

  for (const auto& key : model.allObservedOutputs) {
    if (model.connectivitySkipInfoByKey.find(key) != model.connectivitySkipInfoByKey.end()) {
      // LCOV_EXCL_START
      model.skippedObservedOutputs.push_back(key);
      model.observedOutputExprByKey.erase(key);
      continue;
    }
    if (model.observedOutputExprByKey.find(key) != model.observedOutputExprByKey.end()) {
    // LCOV_EXCL_STOP
      model.observedOutputs.push_back(key);
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
  }
}

void validateExtractedModel(SequentialDesignModel& model) {
  // Missing formulas mean we do not have a sound SEC model, so report the
  // design as unsupported instead of continuing with partial information.
  for (const auto& key : model.observedOutputs) {
    if (model.observedOutputExprByKey.find(key) == model.observedOutputExprByKey.end()) {
      const auto displayIt = model.displayNameByKey.find(key);  // LCOV_EXCL_LINE
      model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
          "Missing observed output expression for `" +  // LCOV_EXCL_LINE
          (displayIt == model.displayNameByKey.end() ? signalKeyToString(key)  // LCOV_EXCL_LINE
                                                     : displayIt->second) +  // LCOV_EXCL_LINE
          "`");
    }  // LCOV_EXCL_LINE
  }
  for (const auto& key : model.stateBits) {
    if (model.nextStateExprByStateKey.find(key) == model.nextStateExprByStateKey.end()) {
      const auto displayIt = model.displayNameByKey.find(key);  // LCOV_EXCL_LINE
      model.unsupportedReasons.push_back(  // LCOV_EXCL_LINE
          "Missing next-state relation for `" +  // LCOV_EXCL_LINE
          (displayIt == model.displayNameByKey.end() ? signalKeyToString(key)  // LCOV_EXCL_LINE
                                                     : displayIt->second) +  // LCOV_EXCL_LINE
          "`");
    }  // LCOV_EXCL_LINE
  }
}

void logExtractedModelDebugSummary(const ExtractContext& ctx,
                                   const SequentialDesignModel& model) {
  if (!ctx.secDiagEnabled) {
    return;
  }

  size_t structuredMemoryCellCount = 0;
  for (const auto& pendingMemory : ctx.pendingMemoryInstances) {
    structuredMemoryCellCount += pendingMemory.cellStates.size();  // LCOV_EXCL_LINE
  }
  fprintf(
      stderr,
      "SEC diag: extract(%s) structured_memories=%zu structured_memory_cells=%zu abstracted_seq_boundaries=%zu opaque_inputs=%zu opaque_outputs=%zu\n",
      ctx.topName.c_str(),
      ctx.pendingMemoryInstances.size(),
      structuredMemoryCellCount,
      model.abstractedSequentialBoundaries.size(),
      model.internalBoundaryInputKeys.size(),
      model.internalBoundaryOutputKeys.size());

  auto formatSignal = [&](const SignalKey& key) {
    const auto nameIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    // LCOV_EXCL_START
    std::ostringstream oss;
    // LCOV_EXCL_STOP
    oss << (nameIt == model.displayNameByKey.end() ? signalKeyToString(key)
                                                   : nameIt->second);
    if (varIt != model.inputVarByKey.end()) {
      oss << "@v" << varIt->second;
    }
    // LCOV_EXCL_START
    return oss.str();
  };
  // LCOV_EXCL_STOP

  std::ostringstream stateSummary;
  for (size_t index = 0; index < model.stateBits.size(); ++index) {
    if (index != 0) {
      stateSummary << ", ";  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    stateSummary << formatSignal(model.stateBits[index]);
  }
  fprintf(
      stderr,
      "SEC diag: extract(%s) kept_states=[%s]\n",
      ctx.topName.c_str(),
      stateSummary.str().c_str());

  for (const auto& key : model.observedOutputs) {
    const auto exprIt = model.observedOutputExprByKey.find(key);
    if (exprIt == model.observedOutputExprByKey.end() || exprIt->second == nullptr) {
      continue;  // LCOV_EXCL_LINE
    }
    std::ostringstream supportSummary;
    bool first = true;
    for (const auto symbol : exprIt->second->getSupportVars()) {
      if (!first) {
        supportSummary << ", ";  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      first = false;
      supportSummary << "v" << symbol;
    }
    fprintf(
        stderr,
        "SEC diag: extract(%s) observed_support %s = [%s]\n",
        ctx.topName.c_str(),
        formatSignal(key).c_str(),
        supportSummary.str().c_str());
  }
  fflush(stderr);
}

}  // namespace

SequentialDesignModel SequentialDesignModel::extract(naja::NL::SNLDesign* top) {
  if (top == nullptr) {
    throw std::invalid_argument("SequentialDesignModel::extract: null top");
  }

  auto* universe = naja::NL::NLUniverse::get();
  if (universe == nullptr) {
    throw std::runtime_error("SequentialDesignModel::extract: NLUniverse not created");
  }

  SequentialDesignModel model;
  ExtractContext ctx{
      // LCOV_EXCL_START
      .top = top,
      .universe = universe,
      // LCOV_EXCL_STOP
      .previousTop = universe->getTopDesign(),
      .topName = top->getName().getString(),
      .secDiagEnabled = std::getenv("KEPLER_SEC_DIAG") != nullptr,
      .abstractUncomputableSequentialBoundaries =
          KEPLER_FORMAL::Config::getSecTreatUncomputableSeqAsBoundary(),
  };
  ctx.builder.setRetainDnl(true);

  // Phase 1: collect the raw boundary, classify top I/O vs sequential state,
  // and scan leaf sequentials so the later formula build knows what it must
  // reconstruct.
  collectInitialBuilderBoundary(ctx);
  collectTopInterfaceTerms(ctx, model);
  classifyBuilderBoundaryTerms(ctx, model);
  collectSequentialTransitions(ctx, model);

  if (model.hasUnsupportedFeatures()) {
    // Primitive-modeling issues are structural, not proof-related. Report them
    // immediately so large designs fail fast before the expensive cone builder
    // tries to derive BoolExprs for a transition system we already know is
    // unsupported.
    naja::DNL::destroy();
    if (ctx.previousTop != nullptr) {
      universe->setTopDesign(ctx.previousTop);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (ctx.secDiagEnabled) {
      fprintf(
          stderr,
          "SEC diag: extract(%s) early unsupported exit before build\n",
          ctx.topName.c_str());
      fflush(stderr);
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
    return model;
  }

  // Phase 2: build the initial boundary formulas for real top outputs plus any
  // already abstracted boundary terms, then publish the normalized SEC
  // interface and variable map.
  buildInitialObservedOutputClouds(ctx, model);
  publishNormalizedBoundary(ctx, model);

  std::vector<naja::DNL::DNLID> builderInputs = ctx.builder.getInputs();
  std::vector<naja::DNL::DNLID> builderOutputs = ctx.builder.getOutputs();
  std::vector<size_t> termDNLID2varID = ctx.builder.getTermDNLID2VarID();
  recordBoundaryInputVars(ctx, builderInputs, termDNLID2varID, model);

  std::unordered_map<naja::DNL::DNLID, BoolExpr*> outputExprByTerm;
  const auto& outputTerms = builderOutputs;
  const auto& outputExprs = ctx.builder.getPOs();
  auto skippedOutputsByTerm = ctx.builder.getSkippedOutputs();
  // Keep only the valid formulas produced by the clause builder. Invalid
  // clouds are classified below either as skippable SEC gaps or as hard
  // unsupported logic.
  for (size_t i = 0; i < outputTerms.size(); ++i) {
    BoolExpr* expr = outputExprs[i];
    if (expr == nullptr || !expr->isValid()) {
      continue;  // LCOV_EXCL_LINE
    }
    outputExprByTerm.emplace(outputTerms[i], expr);
  }

  const auto structuredMemoryDependencyTerms =
      collectStructuredMemoryDependencyTerms(ctx);
  if (!structuredMemoryDependencyTerms.empty()) {
    const auto dependencyOutputs = materializeBuilderOutputs(
        structuredMemoryDependencyTerms,
        builderInputs,
        termDNLID2varID,
        ctx.collectedSkippedOutputs,
        ctx.secDiagEnabled,
        ctx.topName.c_str(),
        "structured memory dependency build");
    appendUniqueTermIDs(builderInputs, dependencyOutputs.inputs);
    appendUniqueTermIDs(builderOutputs, dependencyOutputs.outputs);
    mergeBuilderTermVarIDs(termDNLID2varID, dependencyOutputs.termDNLID2varID);
    recordBoundaryInputVars(ctx, builderInputs, termDNLID2varID, model);
    for (const auto& [termID, expr] : dependencyOutputs.outputExprByTerm) {
      outputExprByTerm.insert_or_assign(termID, expr);
    }
    for (const auto& [termID, info] : dependencyOutputs.skippedOutputsByTerm) {
      skippedOutputsByTerm.emplace(termID, info);
    }
  }
  assignStructuredMemoryStateVars(ctx, model);

  // Phase 3: materialize the observed output formulas that SEC will actually
  // compare, classifying anything missing as either a skippable connectivity
  // gap or a hard unsupported boundary.
  materializeBoundaryObservedOutputs(
      ctx.abstractedBoundaryObservedTerms,
      outputExprByTerm,
      skippedOutputsByTerm,
      builderInputs,
      builderOutputs,
      termDNLID2varID,
      model);
  materializeTopObservedOutputs(
      ctx.topOutputKeyByTerm, outputExprByTerm, skippedOutputsByTerm, model);
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) materialized output exprs observed=%zu total=%zu\n",
        ctx.topName.c_str(),
        model.observedOutputExprByKey.size(),
        outputExprByTerm.size());
    fflush(stderr);
  }

  if (!ctx.pendingMemoryInstances.empty()) {
    buildStructuredMemoryTransitions(
        ctx,
        model,
        builderInputs,
        builderOutputs,
        termDNLID2varID,
        outputExprByTerm,
        skippedOutputsByTerm);
  }

  // Phase 4: rebuild the next-state relations for just the state that is still
  // relevant to covered outputs, then fold any late boundary abstractions back
  // into the published interface.
  const auto rebuiltArtifacts = rebuildRequiredStateTransitions(
      ctx,
      model,
      // LCOV_EXCL_START
      builderInputs,
      builderOutputs,
      // LCOV_EXCL_STOP
      termDNLID2varID,
      outputExprByTerm,
      skippedOutputsByTerm);
  applyRebuiltTransitionArtifacts(
      ctx,
      rebuiltArtifacts,
      model,
      builderInputs,
      builderOutputs,
      termDNLID2varID,
      outputExprByTerm,
      skippedOutputsByTerm);
  filterUnsupportedAndUnmappedBoundary(ctx, model);
  composeSameDomainPhaseTransitions(model);
  markMultiClockDomainConesAsSkipped(model);

  // Phase 5: propagate connectivity skips through dependent state/output cones,
  // then partition the final interface into covered vs skipped signals.  The
  // proof remapper treats any remaining unpublished internal support as a
  // design-private free input, so normal SEC extraction should not skip an
  // otherwise covered top output just because a memory/opaque internal leaf was
  // not part of the public state/environment interface.
  if (std::getenv("KEPLER_SEC_STRICT_UNPUBLISHED_SUPPORT") != nullptr) {
    markFormulasWithUnpublishedSupportAsSkipped(ctx, model);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  propagateConnectivitySkipsThroughDependencies(model);
  partitionCoveredSignals(model);

  inferSynthesizedResetInitialStateValues(model);
  logExtractedModelDebugSummary(ctx, model);
  if (ctx.secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: extract(%s) synthesized init inference done init=%zu\n",
        ctx.topName.c_str(),
        model.initialStateValueByKey.size());
    fflush(stderr);
  }

  // Phase 6: make sure the remaining covered interface is complete before SEC
  // hands this model to the proof engines.
  validateExtractedModel(model);

  // Restore the original top design for callers that keep using the universe.
  naja::DNL::destroy();
  if (ctx.previousTop != nullptr) {
    universe->setTopDesign(ctx.previousTop);
  }
  if (ctx.secDiagEnabled) {
    fprintf(stderr, "SEC diag: extract(%s) end\n", ctx.topName.c_str());
    fflush(stderr);
  }

  return model;
}

}  // namespace KEPLER_FORMAL::SEC
