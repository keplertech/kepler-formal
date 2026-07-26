// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BoolExprCache.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLLibertyConstructor.h"
#include "SNLSVConstructor.h"
#include "SNLPath.h"
#include "SNLBusNet.h"
#include "SNLBusNetBit.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLInstance.h"
#include "SNLScalarNet.h"
#include "SNLScalarTerm.h"
#include "common/AlignedSignals.h"
#include "common/BoolExprUtils.h"
#include "common/PrivateProofSymbol.h"
#include "common/ProofProblemDebug.h"
#include "common/SecDiag.h"
#include "imc/CraigInterpolatingModelChecker.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "imc/IMCEngine.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/OutputBatching.h"
#include "pdr/PDREngine.h"
#include "kinduction/BaseCaseSolver.h"
#include "kinduction/SatEncoding.h"
#include "kinduction/InductionStepSolver.h"
#include "model/SequentialDesignModel.h"
#include "proof/TransitionExprResolver.h"
#include "BuildPrimaryOutputClauses.h"
#include "Tree2BoolExpr.h"
#include "clocks/SecClockModel.h"
#include "strategy/SequentialEquivalenceStrategy.h"

using namespace naja::NL;
using namespace KEPLER_FORMAL::SEC;
using KEPLER_FORMAL::BoolExpr;

namespace KEPLER_FORMAL::SEC::detail {

namespace {

struct PendingPinTermForTest {
  naja::DNL::DNLID termID = naja::DNL::DNLID_MAX;
  naja::NL::NLID::Bit bit = 0;
};

struct PendingTransitionForTest {
  naja::DNL::DNLID stateTermID = naja::DNL::DNLID_MAX;
  naja::NL::NLID::Bit stateBit = 0;
  size_t independentStateOutputCount = 0;
  std::unordered_map<std::string, std::vector<PendingPinTermForTest>> pinTermIDs;
};

size_t countTextOccurrences(
    const std::string& text,
    const std::string& needle) {
  size_t count = 0;
  size_t pos = 0;
  while ((pos = text.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

struct ConeTraceForTest {
  std::vector<std::vector<std::string>> levels;
  std::set<std::string> allTerms;
};

struct ConeDiffReportForTest {
  ConeTraceForTest trace;
  std::string error;
};

struct ScopedDnlContextForTest {
  explicit ScopedDnlContextForTest(naja::NL::SNLDesign* top)
      : universe_(naja::NL::NLUniverse::get()),
        previousTop_(universe_ ? universe_->getTopDesign() : nullptr) {
    if (universe_ == nullptr) {
      throw std::runtime_error("NLUniverse not created for SEC cone tracing");
    }

    naja::DNL::destroy();
    universe_->setTopDesign(top);
    dnl_ = naja::DNL::get();
  }

  ~ScopedDnlContextForTest() {
    naja::DNL::destroy();
    if (universe_ != nullptr && previousTop_ != nullptr) {
      universe_->setTopDesign(previousTop_);
    }
  }

  naja::DNL::DNLFull* dnl() const {
    return dnl_;
  }

 private:
  naja::NL::NLUniverse* universe_ = nullptr;
  naja::NL::SNLDesign* previousTop_ = nullptr;
  naja::DNL::DNLFull* dnl_ = nullptr;
};

std::string formatBoolValueForTest(bool value) {
  return value ? "1" : "0";
}


std::optional<naja::DNL::DNLID> resolvePendingPinTermIDForTest(
    const PendingTransitionForTest& pending,
    const char* pinName) {
  const auto pinIt = pending.pinTermIDs.find(pinName);
  if (pinIt == pending.pinTermIDs.end()) {
    return std::nullopt;
  }

  const auto& candidates = pinIt->second;
  if (candidates.empty()) {
    return std::nullopt;
  }

  if (candidates.size() > 1) {
    for (const auto& candidate : candidates) {
      if (candidate.bit == pending.stateBit) {
        return candidate.termID;
      }
    }
    throw std::runtime_error(
        "Missing bit-matched sequential pin `" + std::string(pinName) + "`");
  }

  const bool isDataPin = std::string(pinName) == "D";
  if (isDataPin && pending.independentStateOutputCount > 1) {
    throw std::runtime_error(
        "Shared scalar D input cannot define multiple independent state outputs");
  }

  return candidates.front().termID;
}

BoolExpr* getRequiredOutputExprForTest(
    const PendingTransitionForTest& pending,
    const char* pinName,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm) {
  const auto resolvedTermID = resolvePendingPinTermIDForTest(pending, pinName);
  if (!resolvedTermID.has_value()) {
    return nullptr;
  }
  const auto exprIt = outputExprByTerm.find(*resolvedTermID);
  if (exprIt == outputExprByTerm.end()) {
    throw std::runtime_error(
        "Missing combinational expression for sequential pin `" +
        std::string(pinName) + "`");
  }
  return exprIt->second;
}

std::vector<std::string> setDifferenceForTest(const std::set<std::string>& lhs,
                                              const std::set<std::string>& rhs) {
  std::vector<std::string> diff;
  std::set_difference(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(diff));
  return diff;
}

std::string describeMismatchedNamesForTest(const std::vector<std::string>& lhs,
                                           const std::vector<std::string>& rhs,
                                           const char* label) {
  std::ostringstream oss;
  oss << "Mismatched " << label << " sets";
  if (!lhs.empty()) {
    oss << " lhs=[";
    for (size_t i = 0; i < lhs.size(); ++i) {
      if (i) {
        oss << ", ";
      }
      oss << lhs[i];
    }
    oss << "]";
  }
  if (!rhs.empty()) {
    oss << " rhs=[";
    for (size_t i = 0; i < rhs.size(); ++i) {
      if (i) {
        oss << ", ";
      }
      oss << rhs[i];
    }
    oss << "]";
  }
  return oss.str();
}

std::map<SignalKey, std::string, SignalKeyLess> buildKeyToNameMapForTest(
    const std::vector<SignalKey>& keys,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames,
    const char* label) {
  std::map<SignalKey, std::string, SignalKeyLess> byKey;
  for (const auto& key : keys) {
    const auto nameIt = displayNames.find(key);
    if (nameIt == displayNames.end()) {
      throw std::runtime_error(
          std::string("Missing display name for SEC ") + label);
    }
    const auto [_, inserted] = byKey.emplace(key, nameIt->second);
    if (!inserted) {
      throw std::runtime_error(
          std::string("Duplicate SEC ") + label + " key `" +
          signalKeyToString(key) + "`");
    }
  }
  return byKey;
}

std::vector<std::string> sortedNamesForTest(
    const std::map<SignalKey, std::string, SignalKeyLess>& byKey) {
  std::vector<std::string> names;
  names.reserve(byKey.size());
  for (const auto& [_, name] : byKey) {
    names.push_back(name);
  }
  return names;
}

std::optional<naja::DNL::DNLID> findTermByDisplayNameForTest(
    naja::DNL::DNLFull* dnl,
    const std::string& signalName);

std::string getTerminalDisplayNameForTest(
    const naja::DNL::DNLTerminalFull& terminal);

std::vector<naja::DNL::DNLID> resolveTermsByKeyForTest(
    naja::DNL::DNLFull* dnl,
    const std::vector<SignalKey>& keys);

std::string formatConeTermForTest(naja::DNL::DNLFull* dnl,
                                  naja::DNL::DNLID termID);

ConeTraceForTest buildConeTraceForTest(
    naja::DNL::DNLFull* dnl,
    naja::DNL::DNLID seedTermID,
    const std::vector<naja::DNL::DNLID>& environmentInputs);

std::string formatConeTracebackForTest(
    const KInductionResult::CounterexampleWitness& witness,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1);

}  // namespace

BoolExpr* buildNextStateExprForTest(
    size_t stateTermID,
    const std::unordered_map<std::string, naja::DNL::DNLID>& pinTermIDs,
    const std::vector<size_t>& termDNLID2varID,
    const std::unordered_map<naja::DNL::DNLID, BoolExpr*>& outputExprByTerm) {
  PendingTransitionForTest pending;
  pending.stateTermID = stateTermID;
  pending.independentStateOutputCount = 1;
  for (const auto& [pinName, termID] : pinTermIDs) {
    pending.pinTermIDs[pinName].push_back({termID, 0});
  }

  if (pending.stateTermID >= termDNLID2varID.size()) {
    throw std::runtime_error("Sequential state term is out of range");
  }

  const size_t stateVarID = termDNLID2varID[pending.stateTermID];
  if (stateVarID < 2) {
    throw std::runtime_error("Sequential state bit was mapped to a constant");
  }

  BoolExpr* data = getRequiredOutputExprForTest(pending, "D", outputExprByTerm);
  if (data == nullptr) {
    throw std::runtime_error("Unsupported sequential primitive without D input");
  }

  BoolExpr* current = BoolExpr::Var(stateVarID);
  BoolExpr* next = data;

  if (BoolExpr* enable = getRequiredOutputExprForTest(pending, "E", outputExprByTerm)) {
    next = BoolExpr::Or(
        BoolExpr::And(enable, data),
        BoolExpr::And(BoolExpr::Not(enable), current));
  }

  BoolExpr* resetHigh =
      getRequiredOutputExprForTest(pending, "R", outputExprByTerm);
  BoolExpr* resetLow =
      getRequiredOutputExprForTest(pending, "RN", outputExprByTerm);
  BoolExpr* setHigh =
      getRequiredOutputExprForTest(pending, "S", outputExprByTerm);
  BoolExpr* setLow =
      getRequiredOutputExprForTest(pending, "SN", outputExprByTerm);

  auto applyForcedValue = [&](BoolExpr* asserted, bool value) {
    return BoolExpr::Or(
        BoolExpr::And(asserted,
                      value ? BoolExpr::createTrue() : BoolExpr::createFalse()),
        BoolExpr::And(BoolExpr::Not(asserted), next));
  };

  // Match production ordering: set/clear controls override reset/preset when
  // both are asserted, which keeps ASAP7 reset+set flops modeled in SEC.
  if (resetHigh) {
    next = applyForcedValue(resetHigh, false);
  }
  if (resetLow) {
    next = applyForcedValue(BoolExpr::Not(resetLow), false);
  }
  if (setHigh) {
    next = applyForcedValue(setHigh, true);
  }
  if (setLow) {
    next = applyForcedValue(BoolExpr::Not(setLow), true);
  }

  return next;
}

std::string formatStringListForTest(const std::vector<std::string>& values,
                                    size_t limit) {
  if (values.empty()) {
    return "<none>";
  }

  std::ostringstream oss;
  const size_t printed = std::min(values.size(), limit);
  for (size_t i = 0; i < printed; ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  if (values.size() > printed) {
    oss << ", ... +" << (values.size() - printed) << " more";
  }
  return oss.str();
}

std::string formatConeLevelsForTest(
    const std::vector<std::vector<std::string>>& levels) {
  constexpr size_t kMaxLevels = 12;
  constexpr size_t kMaxTermsPerLevel = 12;

  if (levels.empty()) {
    return "    <no traced cone terms>\n";
  }

  std::ostringstream oss;
  const size_t printedLevels = std::min(levels.size(), kMaxLevels);
  for (size_t level = 0; level < printedLevels; ++level) {
    oss << "    step " << level << ": "
        << formatStringListForTest(levels[level], kMaxTermsPerLevel) << "\n";
  }
  if (levels.size() > printedLevels) {
    oss << "    ... +" << (levels.size() - printedLevels)
        << " more trace steps\n";
  }
  return oss.str();
}

std::string formatCounterexampleWitnessForTest(
    const KInductionResult& result,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1) {
  if (!result.witness.has_value()) {
    return "";
  }

  const auto& witness = *result.witness;
  std::ostringstream oss;
  oss << "Counterexample reaches the first bad frame at cycle "
      << witness.badFrame << ".\n";

  if (witness.inputTrace.empty()) {
    oss << "Input trace: <none>\n";
  } else {
    oss << "Input trace:\n";
    for (const auto& frame : witness.inputTrace) {
      oss << "  cycle " << frame.frame << ": ";
      if (frame.assignments.empty()) {
        oss << "<no environment inputs>";
      } else {
        for (size_t i = 0; i < frame.assignments.size(); ++i) {
          if (i) {
            oss << ", ";
          }
          oss << frame.assignments[i].signal << "="
              << formatBoolValueForTest(frame.assignments[i].value);
        }
      }
      oss << "\n";
    }
  }

  if (!witness.outputMismatches.empty()) {
    oss << "Observed output mismatches at cycle " << witness.badFrame << ":\n";
    for (const auto& mismatch : witness.outputMismatches) {
      oss << "  " << mismatch.signal << ": design0="
          << formatBoolValueForTest(mismatch.design0Value)
          << ", design1=" << formatBoolValueForTest(mismatch.design1Value) << "\n";
    }
  }

  oss << formatConeTracebackForTest(witness, model0, model1, top0, top1);
  return oss.str();
}

AlignedSignals alignSignalsByNameForTest(
    const std::vector<SignalKey>& keys0,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames0,
    const std::vector<SignalKey>& keys1,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames1,
    const char* label) {
  const auto byKey0 = buildKeyToNameMapForTest(keys0, displayNames0, label);
  const auto byKey1 = buildKeyToNameMapForTest(keys1, displayNames1, label);
  if (byKey0.size() != byKey1.size()) {
    throw std::runtime_error(describeMismatchedNamesForTest(
        sortedNamesForTest(byKey0), sortedNamesForTest(byKey1), label));
  }

  auto it0 = byKey0.begin();
  auto it1 = byKey1.begin();
  for (; it0 != byKey0.end() && it1 != byKey1.end(); ++it0, ++it1) {
    if (it0->first != it1->first) {
      throw std::runtime_error(describeMismatchedNamesForTest(
          sortedNamesForTest(byKey0), sortedNamesForTest(byKey1), label));
    }
  }

  AlignedSignals aligned;
  aligned.names.reserve(byKey0.size());
  aligned.keys0.reserve(byKey0.size());
  aligned.keys1.reserve(byKey0.size());
  for (const auto& [key, displayName] : byKey0) {
    aligned.names.push_back(displayName);
    aligned.keys0.push_back(key);
    aligned.keys1.push_back(key);
  }
  return aligned;
}

SignalKey getTerminalPathKeyForTest(
    const naja::DNL::DNLTerminalFull& terminal) {
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

std::optional<naja::DNL::DNLID> findTermByKeyForTest(
    naja::DNL::DNLFull* dnl,
    const SignalKey& key) {
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    if (getTerminalPathKeyForTest(term) == key) {
      return termID;
    }
  }
  return std::nullopt;
}

std::vector<naja::DNL::DNLID> selectRequiredBuilderOutputsForTest(
    const std::vector<naja::DNL::DNLID>& collectedOutputs,
    const std::unordered_set<naja::DNL::DNLID>& topOutputTerms,
    const std::vector<naja::DNL::DNLID>& sequentialDependencyTerms,
    const std::unordered_set<naja::DNL::DNLID>& prunedBuilderOutputTerms) {
  const std::unordered_set<naja::DNL::DNLID> sequentialDependencySet(
      sequentialDependencyTerms.begin(), sequentialDependencyTerms.end());
  std::vector<naja::DNL::DNLID> filteredOutputs;
  filteredOutputs.reserve(collectedOutputs.size());

  for (const auto outputTermID : collectedOutputs) {
    if (prunedBuilderOutputTerms.find(outputTermID) !=
        prunedBuilderOutputTerms.end()) {
      continue;
    }
    if (topOutputTerms.find(outputTermID) != topOutputTerms.end() ||
        sequentialDependencySet.find(outputTermID) != sequentialDependencySet.end()) {
      filteredOutputs.push_back(outputTermID);
    }
  }

  return filteredOutputs;
}

namespace {

std::string getTerminalDisplayNameForTest(
    const naja::DNL::DNLTerminalFull& terminal) {
  std::ostringstream oss;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  for (const auto& name : pathNames) {
    oss << name.getString() << ".";
  }
  oss << terminal.getSnlBitTerm()->getName().getString() << "["
      << terminal.getSnlBitTerm()->getBit() << "]";
  return oss.str();
}

struct DisplayNameQueryForTest {
  static DisplayNameQueryForTest exact(const std::string& displayName) {
    DisplayNameQueryForTest query;
    query.parse(displayName, /*requiresBit=*/true);
    return query;
  }

  static DisplayNameQueryForTest prefix(const std::string& displayPrefix) {
    DisplayNameQueryForTest query;
    query.parse(displayPrefix, /*requiresBit=*/false);
    return query;
  }

  bool canPrefilter() const { return valid; }

  bool matchesLeaf(const naja::DNL::DNLTerminalFull& terminal) const {
    if (!valid) {
      return true;
    }
    auto* bitTerm = terminal.getSnlBitTerm();
    if (bitTerm == nullptr) {
      return false;
    }
    if (bitTerm->getName().getString() != leafName) {
      return false;
    }
    if (!bitText.empty() && std::to_string(bitTerm->getBit()) != bitText) {
      return false;
    }
    return true;
  }

 private:
  void parse(const std::string& displayText, bool requiresBit) {
    const size_t bracket = displayText.rfind('[');
    if (bracket == std::string::npos) {
      return;
    }
    const size_t dot = displayText.rfind('.', bracket);
    if (dot == std::string::npos || dot + 1 >= bracket) {
      return;
    }
    if (requiresBit) {
      const size_t close = displayText.rfind(']');
      if (close == std::string::npos || close + 1 != displayText.size() ||
          close <= bracket + 1) {
        return;
      }
      bitText = displayText.substr(bracket + 1, close - bracket - 1);
    }
    leafName = displayText.substr(dot + 1, bracket - dot - 1);
    valid = !leafName.empty();
  }

  bool valid = false;
  std::string leafName;
  std::string bitText;
};

std::optional<naja::DNL::DNLID> findTermByDisplayNameForTest(
    naja::DNL::DNLFull* dnl,
    const std::string& signalName) {
  const DisplayNameQueryForTest query =
      DisplayNameQueryForTest::exact(signalName);
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    // CVA6-scale tests can contain millions of terminals.  Avoid building a
    // hierarchical DNL path unless the cheap leaf-name/bit filter matches.
    if (query.canPrefilter() && !query.matchesLeaf(term)) {
      continue;
    }
    if (getTerminalDisplayNameForTest(term) == signalName) {
      return termID;
    }
  }
  return std::nullopt;
}

std::optional<naja::DNL::DNLID> findFirstTermByDisplayPrefixForTest(
    naja::DNL::DNLFull* dnl,
    const std::string& signalPrefix) {
  const DisplayNameQueryForTest query =
      DisplayNameQueryForTest::prefix(signalPrefix);
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    if (query.canPrefilter() && !query.matchesLeaf(term)) {
      continue;
    }
    if (getTerminalDisplayNameForTest(term).rfind(signalPrefix, 0) == 0) {
      return termID;
    }
  }
  return std::nullopt;
}

std::optional<naja::DNL::DNLID> findBuildableOutputRootForTest(
    naja::DNL::DNLFull* dnl,
    naja::DNL::DNLID requestedTermID,
    std::vector<std::string>* chain = nullptr) {
  std::unordered_set<naja::DNL::DNLID> visitedTerms;
  naja::DNL::DNLID currentTermID = requestedTermID;
  while (currentTermID != naja::DNL::DNLID_MAX &&
         visitedTerms.insert(currentTermID).second) {
    const auto& currentTerm = dnl->getDNLTerminalFromID(currentTermID);
    if (currentTerm.isNull()) {
      return std::nullopt;
    }
    if (chain != nullptr) {
      chain->push_back(getTerminalDisplayNameForTest(currentTerm));
    }
    if (currentTerm.isTopPort() &&
        currentTerm.getSnlBitTerm()->getDirection() !=
            naja::NL::SNLBitTerm::Direction::Output) {
      return currentTermID;
    }
    if (currentTerm.getSnlBitTerm()->getDirection() ==
        naja::NL::SNLBitTerm::Direction::Output) {
      const auto& inst = currentTerm.getDNLInstance();
      auto* model = inst.getSNLModel();
      if (model != nullptr && naja::NL::NLDB0::isAssign(model)) {
        std::optional<naja::DNL::DNLID> passthroughDriver;
        for (auto* inputBitTerm :
             naja::NL::SNLDesignModeling::getCombinatorialInputs(
                 const_cast<naja::NL::SNLBitTerm*>(currentTerm.getSnlBitTerm()))) {
          if (inputBitTerm == nullptr ||
              inputBitTerm->getDirection() ==
                  naja::NL::SNLBitTerm::Direction::Output) {
            continue;
          }
          const auto& inputTerm = inst.getTerminalFromBitTerm(inputBitTerm);
          if (inputTerm.isNull() || inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
            passthroughDriver.reset();
            break;
          }
          const auto& iso =
              dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
          if (iso.isConstant() || iso.getDrivers().size() != 1) {
            passthroughDriver.reset();
            break;
          }
          if (passthroughDriver.has_value()) {
            passthroughDriver.reset();
            break;
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
      return std::nullopt;
    }
    const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
    if (iso.isConstant() || iso.getDrivers().size() != 1) {
      return std::nullopt;
    }
    currentTermID = iso.getDrivers().front();
  }
  return std::nullopt;
}

struct BuilderOutputProbeForTest {
  std::optional<naja::DNL::DNLID> normalizedRoot;
  std::string normalizedRootName;
  std::string normalizedRootModelName;
  bool hasBuiltExpr = false;
  bool hasSkip = false;
  std::string skipDetail;
  std::vector<std::string> normalizationChain;
  std::vector<std::string> rootSupportTerms;
  std::vector<std::string> rootCombinationalInputs;
  std::vector<std::string> driverSpine;
};

BuilderOutputProbeForTest probeRequestedBuilderOutputForTest(
    naja::DNL::DNLFull* dnl,
    naja::DNL::DNLID requestedTermID) {
  BuilderOutputProbeForTest probe;
  probe.normalizedRoot =
      findBuildableOutputRootForTest(dnl, requestedTermID, &probe.normalizationChain);
  if (!probe.normalizedRoot.has_value()) {
    return probe;
  }
  probe.normalizedRootName =
      getTerminalDisplayNameForTest(dnl->getDNLTerminalFromID(*probe.normalizedRoot));
  const auto& rootTerm = dnl->getDNLTerminalFromID(*probe.normalizedRoot);
  const auto& rootInst = rootTerm.getDNLInstance();
  probe.normalizedRootModelName = rootInst.getSNLModel()->getName().getString();
  for (naja::DNL::DNLID termID = rootInst.getTermIndexes().first;
       termID <= rootInst.getTermIndexes().second; ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull() || term.getSnlBitTerm()->getDirection() ==
                             naja::NL::SNLBitTerm::Direction::Output) {
      continue;
    }
    std::ostringstream support;
    support << getTerminalDisplayNameForTest(term);
    if (term.getIsoID() == naja::DNL::DNLID_MAX) {
      support << " iso=<invalid>";
    } else {
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
      support << " iso=" << iso.getIsoID()
              << " drivers=" << iso.getDrivers().size()
              << " readers=" << iso.getReaders().size()
              << " const0=" << (iso.isConstant0() ? "true" : "false")
              << " const1=" << (iso.isConstant1() ? "true" : "false")
              << " driver_terms=[";
      for (size_t index = 0; index < iso.getDrivers().size(); ++index) {
        if (index != 0) {
          support << ", ";
        }
        const auto& driverTerm = dnl->getDNLTerminalFromID(iso.getDrivers()[index]);
        support << getTerminalDisplayNameForTest(driverTerm)
                << "{dir="
                << static_cast<int>(driverTerm.getSnlBitTerm()->getDirection())
                << ",model="
                << driverTerm.getDNLInstance().getSNLModel()->getName().getString()
                << "}";
      }
      support << "]";
    }
    probe.rootSupportTerms.push_back(support.str());
  }
  for (auto* bitTerm :
       naja::NL::SNLDesignModeling::getCombinatorialInputs(rootTerm.getSnlBitTerm())) {
    if (bitTerm == nullptr) {
      continue;
    }
    const auto& inputTerm = rootInst.getTerminalFromBitTerm(bitTerm);
    if (inputTerm.isNull()) {
      continue;
    }
    std::ostringstream support;
    support << getTerminalDisplayNameForTest(inputTerm);
    if (inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
      support << " iso=<invalid>";
    } else {
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
      support << " iso=" << iso.getIsoID()
              << " drivers=" << iso.getDrivers().size()
              << " readers=" << iso.getReaders().size()
              << " const0=" << (iso.isConstant0() ? "true" : "false")
              << " const1=" << (iso.isConstant1() ? "true" : "false")
              << " driver_terms=[";
      for (size_t index = 0; index < iso.getDrivers().size(); ++index) {
        if (index != 0) {
          support << ", ";
        }
        const auto& driverTerm = dnl->getDNLTerminalFromID(iso.getDrivers()[index]);
        support << getTerminalDisplayNameForTest(driverTerm)
                << "{dir="
                << static_cast<int>(driverTerm.getSnlBitTerm()->getDirection())
                << ",model="
                << driverTerm.getDNLInstance().getSNLModel()->getName().getString()
                << "}";
      }
      support << "]";
    }
    probe.rootCombinationalInputs.push_back(support.str());
  }
  {
    std::unordered_set<naja::DNL::DNLID> visited;
    naja::DNL::DNLID currentTermID = *probe.normalizedRoot;
    while (visited.insert(currentTermID).second) {
      const auto& currentTerm = dnl->getDNLTerminalFromID(currentTermID);
      std::ostringstream step;
      step << getTerminalDisplayNameForTest(currentTerm)
           << "{model="
           << currentTerm.getDNLInstance().getSNLModel()->getName().getString()
           << "}";
      probe.driverSpine.push_back(step.str());

      std::optional<naja::DNL::DNLID> nextDriver;
      bool ambiguousNext = false;
      const auto& currentInst = currentTerm.getDNLInstance();
      for (auto* bitTerm : naja::NL::SNLDesignModeling::getCombinatorialInputs(
               const_cast<naja::NL::SNLBitTerm*>(currentTerm.getSnlBitTerm()))) {
        if (bitTerm == nullptr ||
            bitTerm->getDirection() ==
                naja::NL::SNLBitTerm::Direction::Output) {
          continue;
        }
        const auto& inputTerm = currentInst.getTerminalFromBitTerm(bitTerm);
        if (inputTerm.isNull() || inputTerm.getIsoID() == naja::DNL::DNLID_MAX) {
          ambiguousNext = true;
          nextDriver.reset();
          break;
        }
        const auto& iso =
            dnl->getDNLIsoDB().getIsoFromIsoIDconst(inputTerm.getIsoID());
        if (iso.isConstant()) {
          continue;
        }
        if (iso.getDrivers().size() != 1) {
          ambiguousNext = true;
          nextDriver.reset();
          break;
        }
        if (nextDriver.has_value()) {
          ambiguousNext = true;
          nextDriver.reset();
          break;
        }
        nextDriver = iso.getDrivers().front();
      }
      if (ambiguousNext || !nextDriver.has_value()) {
        break;
      }
      currentTermID = *nextDriver;
    }
  }

  KEPLER_FORMAL::BuildPrimaryOutputClauses builder;
  builder.collect();
  builder.setOutputs({*probe.normalizedRoot});
  builder.build();

  const auto& outputs = builder.getOutputs();
  const auto& exprs = builder.getPOs();
  for (size_t index = 0; index < outputs.size(); ++index) {
    if (outputs[index] != *probe.normalizedRoot) {
      continue;
    }
    probe.hasBuiltExpr =
        exprs[index] != nullptr && exprs[index]->isValid();
    break;
  }
  if (const auto skippedIt = builder.getSkippedOutputs().find(*probe.normalizedRoot);
      skippedIt != builder.getSkippedOutputs().end()) {
    probe.hasSkip = true;
    probe.skipDetail = skippedIt->second.detail;
  }
  return probe;
}

std::optional<BuilderOutputProbeForTest>
findFirstBuildableOutputProbeByDisplayPrefixForTest(
    naja::DNL::DNLFull* dnl,
    const std::string& signalPrefix) {
  const DisplayNameQueryForTest query =
      DisplayNameQueryForTest::prefix(signalPrefix);
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    // Generated CVA6 memory ports can contain multiple WE bits.  Newer Naja
    // frontend lowering may make an earlier bit normalize through a skipped
    // temporary, so pick the first matching bit that the real clause builder
    // can materialize instead of tying the regression to terminal ordering.
    if (query.canPrefilter() && !query.matchesLeaf(term)) {
      continue;
    }
    if (getTerminalDisplayNameForTest(term).rfind(signalPrefix, 0) != 0) {
      continue;
    }
    BuilderOutputProbeForTest probe =
        probeRequestedBuilderOutputForTest(dnl, termID);
    if (probe.normalizedRoot.has_value() && probe.hasBuiltExpr &&
        !probe.hasSkip) {
      return probe;
    }
  }
  return std::nullopt;
}

std::vector<naja::DNL::DNLID> resolveTermsByKeyForTest(
    naja::DNL::DNLFull* dnl,
    const std::vector<SignalKey>& keys) {
  std::vector<naja::DNL::DNLID> resolved;
  resolved.reserve(keys.size());
  for (const auto& key : keys) {
    if (auto termID = findTermByKeyForTest(dnl, key); termID.has_value()) {
      resolved.push_back(*termID);
    }
  }
  return resolved;
}

std::string formatConeTermForTest(naja::DNL::DNLFull* dnl,
                                  naja::DNL::DNLID termID) {
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.isNull()) {
    return "<null>";
  }
  if (term.getIsoID() == naja::DNL::DNLID_MAX) {
    return getTerminalDisplayNameForTest(term);
  }

  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  if (iso.isConstant0()) {
    return "Constant 0";
  }
  if (iso.isConstant1()) {
    return "Constant 1";
  }
  return getTerminalDisplayNameForTest(term);
}

ConeTraceForTest buildConeTraceForTest(
    naja::DNL::DNLFull* dnl,
    naja::DNL::DNLID seedTermID,
    const std::vector<naja::DNL::DNLID>& environmentInputs) {
  ConeTraceForTest trace;
  std::vector<bool> isEnvironmentInput(dnl->getDNLTerms().size(), false);
  for (const auto termID : environmentInputs) {
    if (termID < isEnvironmentInput.size()) {
      isEnvironmentInput[termID] = true;
    }
  }

  const auto seedIsoID = dnl->getDNLTerminalFromID(seedTermID).getIsoID();
  if (seedIsoID == naja::DNL::DNLID_MAX) {
    return trace;
  }

  std::vector<naja::DNL::DNLID> currentIsos = {seedIsoID};
  std::unordered_set<naja::DNL::DNLID> visitedIsos;
  while (!currentIsos.empty()) {
    std::set<std::string> levelTerms;
    std::vector<naja::DNL::DNLID> nextIsos;

    for (const auto isoID : currentIsos) {
      if (isoID == naja::DNL::DNLID_MAX || !visitedIsos.insert(isoID).second) {
        continue;
      }

      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant0()) {
        levelTerms.insert("Constant 0");
        continue;
      }
      if (iso.isConstant1()) {
        levelTerms.insert("Constant 1");
        continue;
      }

      for (const auto driver : iso.getDrivers()) {
        if (driver == naja::DNL::DNLID_MAX) {
          continue;
        }

        const auto& driverTerm = dnl->getDNLTerminalFromID(driver);
        if (driverTerm.isNull()) {
          continue;
        }

        levelTerms.insert(formatConeTermForTest(dnl, driver));
        if (driver < isEnvironmentInput.size() && isEnvironmentInput[driver]) {
          continue;
        }

        const auto& inst = driverTerm.getDNLInstance();
        for (naja::DNL::DNLID termID = inst.getTermIndexes().first;
             termID != naja::DNL::DNLID_MAX && termID <= inst.getTermIndexes().second;
             ++termID) {
          const auto& term = dnl->getDNLTerminalFromID(termID);
          if (term.isNull()) {
            continue;
          }
          if (term.getSnlBitTerm()->getDirection() ==
              naja::NL::SNLBitTerm::Direction::Output) {
            continue;
          }
          if (term.getIsoID() != naja::DNL::DNLID_MAX) {
            nextIsos.push_back(term.getIsoID());
          }
        }
      }
    }

    if (!levelTerms.empty()) {
      std::vector<std::string> orderedTerms(levelTerms.begin(), levelTerms.end());
      trace.allTerms.insert(orderedTerms.begin(), orderedTerms.end());
      trace.levels.push_back(std::move(orderedTerms));
    }

    std::sort(nextIsos.begin(), nextIsos.end());
    nextIsos.erase(std::unique(nextIsos.begin(), nextIsos.end()), nextIsos.end());
    currentIsos = std::move(nextIsos);
  }

  return trace;
}

ConeDiffReportForTest buildConeDiffReportForTest(
    naja::NL::SNLDesign* top,
    const std::string& differenceSignal,
    const std::vector<SignalKey>& environmentInputs) {
  ConeDiffReportForTest report;
  ScopedDnlContextForTest dnlContext(top);
  auto* dnl = dnlContext.dnl();

  const auto seedTermID = findTermByDisplayNameForTest(dnl, differenceSignal);
  if (!seedTermID.has_value()) {
    report.error =
        "could not resolve the differing SEC signal back into the DNL";
    return report;
  }

  report.trace = buildConeTraceForTest(
      dnl, *seedTermID, resolveTermsByKeyForTest(dnl, environmentInputs));
  return report;
}

std::string formatConeTracebackForTest(
    const KInductionResult::CounterexampleWitness& witness,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1) {
  if (witness.outputMismatches.empty()) {
    return "";
  }
  const auto& differencePoint = witness.outputMismatches.front();

  std::ostringstream oss;
  oss << "Traceback for first differing point `" << differencePoint.signal
      << "` at cycle " << witness.badFrame << ":\n";

  try {
    const auto report0 = buildConeDiffReportForTest(
        top0, differencePoint.signal, model0.environmentInputs);
    const auto report1 = buildConeDiffReportForTest(
        top1, differencePoint.signal, model1.environmentInputs);

    if (!report0.error.empty() || !report1.error.empty()) {
      oss << "  Cone traceback unavailable: ";
      if (!report0.error.empty()) {
        oss << "design0 " << report0.error;
      }
      if (!report0.error.empty() && !report1.error.empty()) {
        oss << "; ";
      }
      if (!report1.error.empty()) {
        oss << "design1 " << report1.error;
      }
      oss << "\n";
      return oss.str();
    }

    oss << "  design0 cone to environment inputs:\n"
        << formatConeLevelsForTest(report0.trace.levels);
    oss << "  design1 cone to environment inputs:\n"
        << formatConeLevelsForTest(report1.trace.levels);

    constexpr size_t kMaxDiffTerms = 20;
    const auto onlyInDesign0 =
        setDifferenceForTest(report0.trace.allTerms, report1.trace.allTerms);
    const auto onlyInDesign1 =
        setDifferenceForTest(report1.trace.allTerms, report0.trace.allTerms);
    oss << "  cone terms only in design0: "
        << formatStringListForTest(onlyInDesign0, kMaxDiffTerms) << "\n";
    oss << "  cone terms only in design1: "
        << formatStringListForTest(onlyInDesign1, kMaxDiffTerms) << "\n";
  } catch (const std::exception& e) {
    oss << "  Cone traceback unavailable: " << e.what() << "\n";
  }

  return oss.str();
}

}  // namespace

}  // namespace KEPLER_FORMAL::SEC::detail

namespace {

class SequentialEquivalenceStrategyTests : public ::testing::Test {
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

SequentialEquivalenceStrategy makeBinarySecStrategy(
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    SecEngine engine = SecEngine::Pdr) {
  return SequentialEquivalenceStrategy(
      top0,
      top1,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      engine,
      SecEncoding::Binary);
}

SequentialEquivalenceStrategy makeBinaryExtractedSecStrategy(
    SecEngine engine = SecEngine::Pdr) {
  return makeBinarySecStrategy(nullptr, nullptr, engine);
}

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

class ScopedUnsetEnvVar {
 public:
  explicit ScopedUnsetEnvVar(const char* name)
      : name_(name) {
    if (const char* current = std::getenv(name_); current != nullptr) {
      previousValue_ = current;
    }
    unsetenv(name_);
  }

  ~ScopedUnsetEnvVar() {
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

// Synthetic tests below do not always build a Naja universe. Production SEC
// keys come from NLName::getID() on real terminals; this local allocator only
// gives those unit-only artificial keys stable, collision-free identities.
naja::NL::NLID::DesignObjectID makeSyntheticSignalNameID(
    const std::string& name) {
  static std::unordered_map<std::string, naja::NL::NLID::DesignObjectID> ids;
  const auto nextID =
      static_cast<naja::NL::NLID::DesignObjectID>(ids.size() + 1);
  return ids.emplace(name, nextID).first->second;
}

SignalKey makeSignalKey(const std::string& name) {
  SignalKey key;
  const auto nameID = makeSyntheticSignalNameID(name);
  key.first.push_back(nameID);
  key.second.push_back(nameID);
  return key;
}

BoolExpr* makeOrChain(const std::vector<size_t>& vars) {
  BoolExpr* expr = BoolExpr::createFalse();
  for (const auto var : vars) {
    expr = BoolExpr::Or(expr, BoolExpr::Var(var));
  }
  return expr;
}

BoolExpr* makeRepeatedSmallSupportCone(size_t varA, size_t varB, size_t depth) {
  BoolExpr* toggle = BoolExpr::And(BoolExpr::Var(varA), BoolExpr::Var(varB));
  BoolExpr* expr = BoolExpr::Var(varA);
  for (size_t i = 0; i < depth; ++i) {
    expr = BoolExpr::Xor(expr, toggle);
  }
  return expr;
}

SequentialDesignModel makeCombinationalExtractedModel(BoolExpr* outputExpr) {
  SequentialDesignModel model;
  const SignalKey inputKey = makeSignalKey("in");
  const SignalKey outputKey = makeSignalKey("out");
  model.environmentInputs = {inputKey};
  model.topInputKeys = {inputKey};
  model.topOutputKeys = {outputKey};
  model.allObservedOutputs = {outputKey};
  model.observedOutputs = {outputKey};
  model.inputVarByKey.emplace(inputKey, 2);
  model.displayNameByKey.emplace(inputKey, "in[0]");
  model.displayNameByKey.emplace(outputKey, "out[0]");
  model.observedOutputExprByKey.emplace(outputKey, outputExpr);
  return model;
}

KInductionProblem makeDualRailComplementedOutputProblemForTest() {
  constexpr size_t p0One = 2;
  constexpr size_t p0Zero = 3;
  constexpr size_t c0One = 4;
  constexpr size_t c0Zero = 5;
  constexpr size_t p1One = 6;
  constexpr size_t p1Zero = 7;
  constexpr size_t c1One = 8;
  constexpr size_t c1Zero = 9;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {p0One, p0Zero, c0One, c0Zero};
  problem.state1Symbols = {p1One, p1Zero, c1One, c1Zero};
  problem.allSymbols = {
      p0One, p0Zero, c0One, c0Zero, p1One, p1Zero, c1One, c1Zero};
  problem.dualRailStatePairs = {
      {p0One, p0Zero}, {c0One, c0Zero}, {p1One, p1Zero}, {c1One, c1Zero}};
  problem.sameFrameStateEqualityPairs0 = {{c0One, p0Zero}, {c0Zero, p0One}};
  problem.sameFrameStateEqualityPairs1 = {{c1One, p1Zero}, {c1Zero, p1One}};

  for (const auto symbol : problem.state0Symbols) {
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  for (const auto symbol : problem.state1Symbols) {
    problem.transitions1.emplace_back(symbol, BoolExpr::Var(symbol));
  }

  // This fixture checks the allowed same-design rail relation only. It must
  // not assume any by-name relation between design0 and design1 internals.
  BoolExpr* complementedOutputEquality = BoolExpr::And(
      makeEqualityExpr(BoolExpr::Var(c0One), BoolExpr::Var(p0Zero)),
      makeEqualityExpr(BoolExpr::Var(c0Zero), BoolExpr::Var(p0One)));
  problem.observedOutputNames = {"complemented_rail_output"};
  problem.observedOutputExprs0 = {complementedOutputEquality};
  problem.observedOutputExprs1 = {BoolExpr::createTrue()};
  problem.property = complementedOutputEquality;
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = complementedOutputEquality;
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);
  problem.totalStateCount = problem.state0Symbols.size() + problem.state1Symbols.size();
  return problem;
}

KInductionProblem makeDualRailLocalFalseInvariantProblemForTest(
    size_t outputCount = 1) {
  constexpr size_t state = 2;
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(state));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{state, BoolExpr::createFalse()}};

  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < outputCount; ++output) {
    problem.observedOutputNames.push_back("out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(state));
    problem.observedOutputExprs1.push_back(BoolExpr::createFalse());
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(
            problem.observedOutputExprs0.back(),
            problem.observedOutputExprs1.back()));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

void addStateBitForTest(SequentialDesignModel& model,
                        const SignalKey& key,
                        size_t localVar,
                        const std::string& displayName,
                        BoolExpr* nextState) {
  model.stateBits.push_back(key);
  model.inputVarByKey.emplace(key, localVar);
  model.displayNameByKey.emplace(key, displayName);
  model.nextStateExprByStateKey.emplace(key, nextState);
}

struct LargeDualRailResidualCase {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  size_t residualOutputs = 0;
};

LargeDualRailResidualCase makeLargeDualRailResidualCaseForTest(
    const std::string& prefix,
    size_t residualOutputs,
    size_t stateBitsPerDesign = 1,
    bool includeImpliedOutput = true) {
  const SignalKey implied = makeSignalKey(prefix + "ImpliedOut");

  LargeDualRailResidualCase testCase;
  testCase.residualOutputs = residualOutputs;
  size_t nextLocalVar = 2;
  for (size_t i = 0; i < stateBitsPerDesign; ++i) {
    const SignalKey state0 =
        makeSignalKey(prefix + "State0_" + std::to_string(i));
    const SignalKey state1 =
        makeSignalKey(prefix + "State1_" + std::to_string(i));
    const size_t stateVar = nextLocalVar++;
    addStateBitForTest(
        testCase.model0,
        state0,
        stateVar,
        prefix + "_large_residual_state0[" + std::to_string(i) + "]",
        BoolExpr::Var(stateVar));
    addStateBitForTest(
        testCase.model1,
        state1,
        stateVar,
        prefix + "_large_residual_state1[" + std::to_string(i) + "]",
        BoolExpr::Var(stateVar));
  }
  if (includeImpliedOutput) {
    testCase.model0.allObservedOutputs.push_back(implied);
    testCase.model0.observedOutputs.push_back(implied);
    testCase.model0.displayNameByKey.emplace(implied, "implied_out[0]");
    testCase.model0.observedOutputExprByKey.emplace(
        implied, BoolExpr::createTrue());

    testCase.model1.allObservedOutputs.push_back(implied);
    testCase.model1.observedOutputs.push_back(implied);
    testCase.model1.displayNameByKey.emplace(implied, "implied_out[0]");
    testCase.model1.observedOutputExprByKey.emplace(
        implied, BoolExpr::createTrue());
  }

  for (size_t i = 0; i < residualOutputs; ++i) {
    const SignalKey out =
        makeSignalKey(prefix + "ResidualOut" + std::to_string(i));
    const std::string outputName =
        "large_residual_out[" + std::to_string(i) + "]";
    testCase.model0.allObservedOutputs.push_back(out);
    testCase.model0.observedOutputs.push_back(out);
    testCase.model0.displayNameByKey.emplace(out, outputName);
    testCase.model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(2));

    testCase.model1.allObservedOutputs.push_back(out);
    testCase.model1.observedOutputs.push_back(out);
    testCase.model1.displayNameByKey.emplace(out, outputName);
    testCase.model1.observedOutputExprByKey.emplace(
        out, BoolExpr::createFalse());
  }

  return testCase;
}

struct DelayedRailMismatchModels {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
};

struct WideFrameZeroMismatchModels {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
};

WideFrameZeroMismatchModels makeWideFrameZeroMismatchModelsForTest(
    size_t outputCount,
    size_t dummyStateCount) {
  const SignalKey input = makeSignalKey("wideFrameZeroInput");
  WideFrameZeroMismatchModels models;
  models.model0.environmentInputs = {input};
  models.model1.environmentInputs = {input};
  models.model0.inputVarByKey.emplace(input, 2);
  models.model1.inputVarByKey.emplace(input, 2);
  models.model0.displayNameByKey.emplace(input, "wide_frame_zero_in[0]");
  models.model1.displayNameByKey.emplace(input, "wide_frame_zero_in[0]");

  size_t nextLocalVar = 3;
  for (size_t i = 0; i < dummyStateCount; ++i) {
    const SignalKey state0 =
        makeSignalKey("wideFrameZeroState0_" + std::to_string(i));
    const SignalKey state1 =
        makeSignalKey("wideFrameZeroState1_" + std::to_string(i));
    addStateBitForTest(
        models.model0,
        state0,
        nextLocalVar,
        "wide_frame_zero_state[" + std::to_string(i) + "]",
        BoolExpr::Var(nextLocalVar));
    addStateBitForTest(
        models.model1,
        state1,
        nextLocalVar,
        "wide_frame_zero_state[" + std::to_string(i) + "]",
        BoolExpr::Var(nextLocalVar));
    ++nextLocalVar;
  }

  for (size_t i = 0; i + 1 < outputCount; ++i) {
    const SignalKey output =
        makeSignalKey("wideFrameZeroStableOut" + std::to_string(i));
    const std::string outputName =
        "wide_frame_zero_stable[" + std::to_string(i) + "]";
    models.model0.allObservedOutputs.push_back(output);
    models.model0.observedOutputs.push_back(output);
    models.model0.displayNameByKey.emplace(output, outputName);
    models.model0.observedOutputExprByKey.emplace(output, BoolExpr::createFalse());

    models.model1.allObservedOutputs.push_back(output);
    models.model1.observedOutputs.push_back(output);
    models.model1.displayNameByKey.emplace(output, outputName);
    models.model1.observedOutputExprByKey.emplace(output, BoolExpr::createFalse());
  }

  const SignalKey probe = makeSignalKey("wideFrameZeroProbe");
  models.model0.allObservedOutputs.push_back(probe);
  models.model0.observedOutputs.push_back(probe);
  models.model0.displayNameByKey.emplace(probe, "wide_frame_zero_probe[0]");
  models.model0.observedOutputExprByKey.emplace(probe, BoolExpr::Var(2));

  models.model1.allObservedOutputs.push_back(probe);
  models.model1.observedOutputs.push_back(probe);
  models.model1.displayNameByKey.emplace(probe, "wide_frame_zero_probe[0]");
  models.model1.observedOutputExprByKey.emplace(
      probe, BoolExpr::Not(BoolExpr::Var(2)));
  return models;
}

SequentialDesignModel makeDelayedRailMismatchModelForTest(
    const std::string& prefix,
    const SignalKey& output,
    const std::string& outputName,
    bool seedNextValue,
    size_t dummyStateCount) {
  SequentialDesignModel model;
  model.allObservedOutputs = {output};
  model.observedOutputs = {output};
  model.displayNameByKey.emplace(output, outputName);

  size_t nextLocalVar = 2;
  const SignalKey seed = makeSignalKey(prefix + "Seed");
  const SignalKey stage0 = makeSignalKey(prefix + "Stage0");
  const SignalKey stage1 = makeSignalKey(prefix + "Stage1");

  const size_t seedVar = nextLocalVar++;
  addStateBitForTest(
      model,
      seed,
      seedVar,
      prefix + ".seed_q[0]",
      seedNextValue ? BoolExpr::createTrue() : BoolExpr::createFalse());

  const size_t stage0Var = nextLocalVar++;
  addStateBitForTest(
      model,
      stage0,
      stage0Var,
      prefix + ".stage0_q[0]",
      BoolExpr::Var(seedVar));

  const size_t stage1Var = nextLocalVar++;
  addStateBitForTest(
      model,
      stage1,
      stage1Var,
      prefix + ".stage1_q[0]",
      BoolExpr::Var(stage0Var));
  model.observedOutputExprByKey.emplace(output, BoolExpr::Var(seedVar));

  for (size_t i = 0; i < dummyStateCount; ++i) {
    const SignalKey dummy = makeSignalKey(
        prefix + "Dummy" + std::to_string(i));
    const size_t dummyVar = nextLocalVar++;
    // The dummies keep the test above the old small-PDR certificate limit while
    // staying outside the output COI, matching the compact CSR regression shape.
    addStateBitForTest(
        model,
        dummy,
        dummyVar,
        prefix + ".dummy_q[" + std::to_string(i) + "]",
        BoolExpr::Var(dummyVar));
  }
  return model;
}

DelayedRailMismatchModels makeDelayedRailMismatchModelsForTest(
    size_t dummyStateCount) {
  const SignalKey output = makeSignalKey("dualRailDelayedPdrOut");
  DelayedRailMismatchModels models;
  models.model0 = makeDelayedRailMismatchModelForTest(
      "left",
      output,
      "delayed_out[0]",
      true,
      dummyStateCount);
  models.model1 = makeDelayedRailMismatchModelForTest(
      "right",
      output,
      "delayed_out[0]",
      false,
      dummyStateCount);
  return models;
}

DelayedRailMismatchModels makeHeldRailModelsForTest(
    const std::string& prefix,
    std::optional<bool> initialValue0,
    std::optional<bool> initialValue1) {
  const SignalKey output = makeSignalKey(prefix + "Output");
  const SignalKey state0 = makeSignalKey(prefix + "State0");
  const SignalKey state1 = makeSignalKey(prefix + "State1");
  DelayedRailMismatchModels models;

  addStateBitForTest(
      models.model0,
      state0,
      /*symbol=*/2,
      prefix + ".left_q[0]",
      BoolExpr::Var(2));
  models.model0.allObservedOutputs = {output};
  models.model0.observedOutputs = {output};
  models.model0.displayNameByKey.emplace(output, prefix + "_out[0]");
  models.model0.observedOutputExprByKey.emplace(output, BoolExpr::Var(2));
  if (initialValue0.has_value()) {
    models.model0.initialStateValueByKey.emplace(state0, *initialValue0);
  }

  addStateBitForTest(
      models.model1,
      state1,
      /*symbol=*/2,
      prefix + ".right_q[0]",
      BoolExpr::Var(2));
  models.model1.allObservedOutputs = {output};
  models.model1.observedOutputs = {output};
  models.model1.displayNameByKey.emplace(output, prefix + "_out[0]");
  models.model1.observedOutputExprByKey.emplace(output, BoolExpr::Var(2));
  if (initialValue1.has_value()) {
    models.model1.initialStateValueByKey.emplace(state1, *initialValue1);
  }
  return models;
}

size_t bitCountForPdrChainStateCount(size_t logicalStateCount) {
  size_t bits = 0;
  size_t encodedStates = 1;
  while (encodedStates < logicalStateCount) {
    encodedStates <<= 1;
    ++bits;
  }
  return std::max<size_t>(bits, 1);
}

BoolExpr* makePdrChainStateExpr(const std::vector<size_t>& symbols,
                                size_t value) {
  BoolExpr* expr = BoolExpr::createTrue();
  for (size_t bit = 0; bit < symbols.size(); ++bit) {
    expr = BoolExpr::And(
        expr,
        (value & (size_t{1} << bit)) != 0
            ? BoolExpr::Var(symbols[bit])
            : BoolExpr::Not(BoolExpr::Var(symbols[bit])));
  }
  return BoolExpr::simplify(expr);
}

BoolExpr* makeSaturatingPdrChainNextBitExpr(
    const std::vector<size_t>& symbols,
    size_t logicalStateCount,
    size_t bitIndex) {
  BoolExpr* expr = BoolExpr::createFalse();
  for (size_t logicalState = 0; logicalState < logicalStateCount;
       ++logicalState) {
    const size_t nextLogicalState =
        logicalState + 1 < logicalStateCount ? logicalState + 1
                                             : logicalState;
    if ((nextLogicalState & (size_t{1} << bitIndex)) == 0) {
      continue;
    }
    expr = BoolExpr::Or(expr, makePdrChainStateExpr(symbols, logicalState));
  }
  return BoolExpr::simplify(expr);
}

KInductionProblem buildLinearChainSecProblem(size_t logicalStateCount) {
  const size_t bitCount = bitCountForPdrChainStateCount(logicalStateCount);

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
        problem.state0Symbols[bit],
        makeSaturatingPdrChainNextBitExpr(
            problem.state0Symbols, logicalStateCount, bit));
    problem.transitions1.emplace_back(
        problem.state1Symbols[bit],
        makeSaturatingPdrChainNextBitExpr(
            problem.state1Symbols, logicalStateCount, bit));
  }

  problem.initialCondition = BoolExpr::And(
      makePdrChainStateExpr(problem.state0Symbols, 0),
      makePdrChainStateExpr(problem.state1Symbols, 0));
  problem.initializedStateCount = problem.allSymbols.size();
  problem.totalStateCount = problem.allSymbols.size();
  problem.observedOutputExprs0 = {
      makePdrChainStateExpr(problem.state0Symbols, logicalStateCount - 1)};
  problem.observedOutputExprs1 = {
      makePdrChainStateExpr(problem.state1Symbols, logicalStateCount - 1)};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0.front(), problem.observedOutputExprs1.front());
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

KInductionProblem buildSharedPdrNarrowProbeProblem() {
  KInductionProblem problem = buildLinearChainSecProblem(4);
  const size_t decoyFirst = problem.allSymbols.back() + 1;
  const size_t decoyDelayed = decoyFirst + 1;
  problem.state0Symbols.push_back(decoyFirst);
  problem.state0Symbols.push_back(decoyDelayed);
  problem.allSymbols.push_back(decoyFirst);
  problem.allSymbols.push_back(decoyDelayed);
  problem.transitions0.emplace_back(decoyFirst, BoolExpr::createTrue());
  problem.transitions0.emplace_back(decoyDelayed, BoolExpr::Var(decoyFirst));
  problem.initialCondition = BoolExpr::And(
      problem.initialCondition,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(decoyFirst)),
          BoolExpr::Not(BoolExpr::Var(decoyDelayed))));
  problem.initializedStateCount += 2;
  problem.totalStateCount += 2;
  problem.usesDualRailStateEncoding = true;
  problem.usesStrictDualRailEqualityProperty = true;

  // The broad parent fails when the independent two-cycle decoy rises. It
  // populates the shared higher-frame context without certifying an invariant
  // that could bypass the narrower child run.
  problem.property = BoolExpr::And(
      problem.property, BoolExpr::Not(BoolExpr::Var(decoyDelayed)));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

KInductionProblem buildCraigResetSecProblem(bool equivalent) {
  constexpr size_t reset = 2;
  constexpr size_t data = 3;
  constexpr size_t state0 = 4;
  constexpr size_t state1 = 5;

  KInductionProblem problem;
  problem.environmentInputNames = {"reset", "data"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {reset, data};
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.allSymbols = {reset, data, state0, state1};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.observedOutputExprs0 = {BoolExpr::Var(state0)};
  problem.observedOutputExprs1 = {BoolExpr::Var(state1)};
  problem.transitions0 = {{
      state0,
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(reset)),
                    BoolExpr::Var(data))}};
  problem.transitions1 = {{
      state1,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          equivalent ? BoolExpr::Var(data)
                     : BoolExpr::Not(BoolExpr::Var(data)))}};
  problem.property =
      makeEqualityExpr(BoolExpr::Var(state0), BoolExpr::Var(state1));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.usesDualRailStateEncoding = true;
  problem.totalStateCount = 13;
  return problem;
}

InterpolantRegion makeStateLiteralCraigRegion(size_t symbol, bool positive) {
  InterpolantRegion region;
  region.type = InterpolantRegion::Type::Normal;
  region.root = {true, symbol, positive};
  return region;
}

BoolExpr* makeConjunctionOfVars(const std::vector<size_t>& symbols) {
  BoolExpr* expr = BoolExpr::createTrue();
  for (const size_t symbol : symbols) {
    expr = BoolExpr::And(expr, BoolExpr::Var(symbol));
  }
  return BoolExpr::simplify(expr);
}

KInductionProblem buildWideSharedConeImcProblem(size_t outputCount,
                                               size_t sharedSupportCount) {
  constexpr size_t firstSharedSymbol = 2;
  const size_t firstOutputSymbol = firstSharedSymbol + sharedSupportCount;

  KInductionProblem problem;
  std::vector<size_t> sharedSupport;
  sharedSupport.reserve(sharedSupportCount);
  problem.state0Symbols.reserve(sharedSupportCount + outputCount);
  problem.allSymbols.reserve(sharedSupportCount + outputCount);
  BoolExpr* init = BoolExpr::createTrue();
  for (size_t index = 0; index < sharedSupportCount + outputCount; ++index) {
    const size_t symbol = firstSharedSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
    if (index < sharedSupportCount) {
      sharedSupport.push_back(symbol);
    }
  }

  BoolExpr* sharedCone = makeConjunctionOfVars(sharedSupport);
  BoolExpr* property = BoolExpr::createTrue();
  for (size_t output = 0; output < outputCount; ++output) {
    const size_t outputSymbol = firstOutputSymbol + output;
    problem.observedOutputNames.push_back("wide_out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(
        BoolExpr::And(sharedCone, BoolExpr::Var(outputSymbol)));
    problem.observedOutputExprs1.push_back(
        BoolExpr::And(sharedCone, BoolExpr::Not(BoolExpr::Var(outputSymbol))));
    property = BoolExpr::And(
        property,
        makeEqualityExpr(
            problem.observedOutputExprs0.back(),
            problem.observedOutputExprs1.back()));
  }

  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.allSymbols.size();
  problem.totalStateCount = problem.allSymbols.size();
  problem.usesDualRailStateEncoding = true;
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

KInductionProblem buildWideSharedConeImcProblem(size_t outputCount) {
  return buildWideSharedConeImcProblem(outputCount, 129);
}

KInductionProblem buildDisjointWideConeImcBatchProblem() {
  constexpr size_t kSupportPerOutput = 129;
  constexpr size_t firstSymbol = 2;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.totalStateCount = kSupportPerOutput * 2;
  problem.state0Symbols.reserve(kSupportPerOutput * 2);
  problem.allSymbols.reserve(kSupportPerOutput * 2);

  for (size_t output = 0; output < 2; ++output) {
    std::vector<size_t> support;
    support.reserve(kSupportPerOutput);
    for (size_t index = 0; index < kSupportPerOutput; ++index) {
      const size_t symbol =
          firstSymbol + output * kSupportPerOutput + index;
      support.push_back(symbol);
      problem.state0Symbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
    }
    problem.observedOutputNames.push_back(
        "disjoint_wide_out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(makeConjunctionOfVars(support));
    problem.observedOutputExprs1.push_back(BoolExpr::createFalse());
  }

  return problem;
}

KInductionProblem buildProjectionSharedImcBatchProblem(
    size_t sharedTransitionStates = 24,
    size_t outputCount = 2) {
  constexpr size_t firstSymbol = 2;
  const size_t firstOutputState = firstSymbol + sharedTransitionStates;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols.reserve(sharedTransitionStates + outputCount);
  problem.allSymbols.reserve(sharedTransitionStates + outputCount);

  std::vector<size_t> sharedStates;
  sharedStates.reserve(sharedTransitionStates);
  BoolExpr* init = BoolExpr::createTrue();
  for (size_t index = 0; index < sharedTransitionStates; ++index) {
    const size_t symbol = firstSymbol + index;
    sharedStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }

  BoolExpr* sharedCone = makeConjunctionOfVars(sharedStates);
  BoolExpr* property = BoolExpr::createTrue();
  for (size_t output = 0; output < outputCount; ++output) {
    const size_t outputState = firstOutputState + output;
    problem.state0Symbols.push_back(outputState);
    problem.allSymbols.push_back(outputState);
    problem.initialStateAssignments.push_back({outputState, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(outputState)));
    // Raw output supports are disjoint, but both next-state functions pull in
    // the same transition cone. This mirrors AES text_out bits, where a
    // one-bit bad predicate expands to a much larger shared Craig projection.
    problem.transitions0.emplace_back(
        outputState, sharedCone);
    problem.observedOutputNames.push_back(
        "projection_shared_out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(outputState));
    problem.observedOutputExprs1.push_back(BoolExpr::createFalse());
    property = BoolExpr::And(
        property,
        makeEqualityExpr(
            problem.observedOutputExprs0.back(),
            problem.observedOutputExprs1.back()));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.allSymbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

KInductionProblem buildBootstrapModelGuidedCraigProjectionProblem(
    size_t supportCount = 96,
    bool assignSupportBootstrap = true,
    bool addDualRailPartners = false) {
  constexpr size_t outputState = 2;
  constexpr size_t firstSupportState = 3;
  const size_t firstPartnerState = firstSupportState + supportCount;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols.push_back(outputState);
  problem.allSymbols.push_back(outputState);
  problem.bootstrapStateAssignments.push_back({outputState, false});

  std::vector<size_t> supportStates;
  supportStates.reserve(supportCount);
  for (size_t index = 0; index < supportCount; ++index) {
    const size_t symbol = firstSupportState + index;
    supportStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    if (assignSupportBootstrap) {
      problem.bootstrapStateAssignments.push_back({symbol, false});
    }
    if (addDualRailPartners) {
      const size_t partnerSymbol = firstPartnerState + index;
      problem.state0Symbols.push_back(partnerSymbol);
      problem.allSymbols.push_back(partnerSymbol);
      if (assignSupportBootstrap) {
        problem.bootstrapStateAssignments.push_back({partnerSymbol, true});
      }
      problem.dualRailStatePairs.push_back(
          DualRailSymbolPair{symbol, partnerSymbol});
    }
  }

  problem.transitions0.emplace_back(
      outputState, makeConjunctionOfVars(supportStates));
  problem.observedOutputNames = {"model_guided_out"};
  problem.observedOutputExprs0 = {BoolExpr::Var(outputState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.property = BoolExpr::Not(BoolExpr::Var(outputState));
  problem.bad = BoolExpr::Var(outputState);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildHighFaninCraigProjectionProblemForTest(
    size_t helperTrackedCount = 4,
    size_t decoySupportCount = 96) {
  constexpr size_t outputState = 2;
  constexpr size_t firstHelperState = 3;
  constexpr size_t firstDecoyState = 100;
  const size_t sharedControlState = firstDecoyState + decoySupportCount + 100;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols.push_back(outputState);
  problem.allSymbols.push_back(outputState);
  problem.bootstrapStateAssignments.push_back({outputState, false});

  std::vector<size_t> decoyStates;
  decoyStates.reserve(decoySupportCount);
  for (size_t index = 0; index < decoySupportCount; ++index) {
    const size_t symbol = firstDecoyState + index;
    decoyStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }
  problem.state0Symbols.push_back(sharedControlState);
  problem.allSymbols.push_back(sharedControlState);

  BoolExpr* outputNext = BoolExpr::Var(sharedControlState);
  for (const size_t symbol : decoyStates) {
    outputNext = BoolExpr::Or(outputNext, BoolExpr::Var(symbol));
  }
  problem.transitions0.emplace_back(outputState, outputNext);

  for (size_t index = 0; index < helperTrackedCount; ++index) {
    const size_t symbol = firstHelperState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
    // sharedControlState is intentionally a high-id symbol.  A sorted fallback
    // would spend the first slice on decoys; fan-in scoring should rank this
    // shared dependency first.
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(sharedControlState));
  }

  problem.bad = BoolExpr::Var(outputState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildTinyModelGuidedBackfillCraigProjectionProblemForTest(
    size_t decoyPairCount,
    size_t controlCount = 1) {
  constexpr size_t outputState = 2;
  constexpr size_t firstControlState = 3;
  const size_t firstDecoyState = firstControlState + controlCount + 100;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {outputState};
  problem.allSymbols = {outputState};
  problem.bootstrapStateAssignments = {{outputState, false}};

  BoolExpr* controlCone = BoolExpr::createTrue();
  for (size_t index = 0; index < controlCount; ++index) {
    const size_t symbol = firstControlState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
    controlCone = BoolExpr::And(controlCone, BoolExpr::Var(symbol));
  }

  BoolExpr* impossibleDecoyCone = BoolExpr::createFalse();
  for (size_t index = 0; index < decoyPairCount; ++index) {
    const size_t symbol = firstDecoyState + index * 2;
    const size_t partnerSymbol = symbol + 1;
    problem.state0Symbols.push_back(symbol);
    problem.state0Symbols.push_back(partnerSymbol);
    problem.allSymbols.push_back(symbol);
    problem.allSymbols.push_back(partnerSymbol);
    problem.complementedStatePairs0.push_back({symbol, partnerSymbol});
    // Same-design complement semantics makes symbol && partnerSymbol
    // unreachable, but both rails remain in the transition support.
    BoolExpr* contradiction =
        BoolExpr::And(BoolExpr::Var(symbol), BoolExpr::Var(partnerSymbol));
    impossibleDecoyCone = BoolExpr::Or(impossibleDecoyCone, contradiction);
  }

  // The full transition support is large, but the SAT witness can only make
  // progress through controlState.  Tiny model-guided slices should therefore
  // be backfilled by the scored bounded selector.
  problem.transitions0.emplace_back(
      outputState,
      BoolExpr::Or(controlCone, impossibleDecoyCone));
  problem.bad = BoolExpr::Var(outputState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

void addOneHotPdrStateSymbols(KInductionProblem& problem, size_t stateCount) {
  problem.state0Symbols.clear();
  problem.allSymbols.clear();
  problem.state0Symbols.reserve(stateCount);
  problem.allSymbols.reserve(stateCount);
  size_t nextSymbol = 2;
  for (size_t state = 0; state < stateCount; ++state) {
    const size_t symbol = nextSymbol++;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }
}

BoolExpr* makeOneHotStateExpr(const std::vector<size_t>& symbols,
                              size_t hotIndex) {
  BoolExpr* expr = BoolExpr::createTrue();
  for (size_t index = 0; index < symbols.size(); ++index) {
    expr = BoolExpr::And(
        expr,
        index == hotIndex ? BoolExpr::Var(symbols[index])
                          : BoolExpr::Not(BoolExpr::Var(symbols[index])));
  }
  return BoolExpr::simplify(expr);
}

void finishOneHotPdrProblem(KInductionProblem& problem, size_t badIndex) {
  problem.initialCondition = makeOneHotStateExpr(problem.state0Symbols, 0);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::Var(problem.state0Symbols[badIndex]);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
}

KInductionProblem buildClassicPdrOneHotUnreachableBadChainProblem(
    size_t proofDepth) {
  KInductionProblem problem;
  addOneHotPdrStateSymbols(problem, proofDepth + 1);

  for (size_t state = 0; state < problem.state0Symbols.size(); ++state) {
    BoolExpr* nextExpr = BoolExpr::createFalse();
    if (state == 0) {
      nextExpr = BoolExpr::Var(problem.state0Symbols[0]);
    } else if (state > 1) {
      nextExpr = BoolExpr::Var(problem.state0Symbols[state - 1]);
    }
    problem.transitions0.emplace_back(problem.state0Symbols[state], nextExpr);
  }

  // Classic safe PDR shape: state 0 is the only reachable state, while the bad
  // state is fed by an unreachable predecessor chain 1 -> 2 -> ... -> bad.
  finishOneHotPdrProblem(problem, proofDepth);
  return problem;
}

KInductionProblem buildClassicPdrOneHotReachableBadChainProblem(
    size_t badDepth) {
  KInductionProblem problem;
  addOneHotPdrStateSymbols(problem, badDepth + 1);

  for (size_t state = 0; state < problem.state0Symbols.size(); ++state) {
    BoolExpr* nextExpr = BoolExpr::createFalse();
    if (state > 0) {
      nextExpr = BoolExpr::Var(problem.state0Symbols[state - 1]);
    }
    if (state + 1 == problem.state0Symbols.size()) {
      nextExpr = BoolExpr::Or(nextExpr, BoolExpr::Var(problem.state0Symbols[state]));
    }
    problem.transitions0.emplace_back(
        problem.state0Symbols[state], BoolExpr::simplify(nextExpr));
  }

  // Classic bug PDR shape: 0 -> 1 -> ... -> bad, so the first bad state appears
  // exactly at badDepth transitions.
  finishOneHotPdrProblem(problem, badDepth);
  return problem;
}


KInductionProblem buildDocumentedBooleanPdrCounterexampleProblem() {
  KInductionProblem problem;
  problem.description = "documented Boolean PDR counterexample miter";
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"q_miter"};
  problem.inputSymbols = {4};
  problem.state0Symbols = {2};
  problem.state1Symbols = {3};
  problem.allSymbols = {2, 3, 4};
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Not(BoolExpr::Var(3)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(3)};
  problem.transitions0.emplace_back(2, BoolExpr::Not(BoolExpr::Var(2)));
  problem.transitions1.emplace_back(3, BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)));
  problem.property = makeEqualityExpr(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.bad = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
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

SNLDesign* createBufModelWithName(NLLibrary* library, const std::string& name) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* input =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  auto* output =
      SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  SNLDesignModeling::addCombinatorialArcs({input}, {output});
  SNLDesignModeling::setTruthTable(
      model, SNLTruthTable(1, 0b10, SNLTruthTable::fullDependencies(1)));
  return model;
}

SNLDesign* createBufModel(NLLibrary* library) {
  return createBufModelWithName(library, "BUF");
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
      SNLTruthTable(2, 0b1000, SNLTruthTable::fullDependencies(2)));
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

SNLDesign* createOpaqueLeafModel(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("OPAQUE"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("A"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Y"));
  return model;
}

SNLDesign* createSinglePortMemoryModel(
    NLLibrary* library,
    const std::string& name,
    bool withReset = false) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName(name));
  auto* clock =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CLK"));
  auto* reset = withReset
      ? SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("RST"))
      : nullptr;
  auto* chipEnable =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("CE"));
  auto* writeEnable =
      SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("WE"));
  auto* address =
      SNLBusTerm::create(model, SNLTerm::Direction::Input, 1, 0, NLName("ADDR"));
  auto* writeData =
      SNLBusTerm::create(model, SNLTerm::Direction::Input, 3, 0, NLName("WDATA"));
  auto* writeMask =
      SNLBusTerm::create(model, SNLTerm::Direction::Input, 3, 0, NLName("WMASK"));
  auto* readData =
      SNLBusTerm::create(model, SNLTerm::Direction::Output, 3, 0, NLName("RDATA"));

  SNLDesignModeling::BitTerms readDataBits;
  SNLDesignModeling::BitTerms readAddressBits;
  SNLDesignModeling::BitTerms writeDataBits;
  SNLDesignModeling::BitTerms writeMaskBits;
  for (int bit = 0; bit <= 3; ++bit) {
    readDataBits.push_back(readData->getBit(bit));
    writeDataBits.push_back(writeData->getBit(bit));
    writeMaskBits.push_back(writeMask->getBit(bit));
    if (bit <= 1) {
      readAddressBits.push_back(address->getBit(bit));
    }
  }
  SNLDesignModeling::addClockToOutputsArcs(clock, readDataBits);
  SNLDesignModeling::BitTerms clockInputs = {
      address->getBit(0),
      address->getBit(1),
      writeData->getBit(0),
      writeData->getBit(1),
      writeData->getBit(2),
      writeData->getBit(3),
      writeMask->getBit(0),
      writeMask->getBit(1),
      writeMask->getBit(2),
      writeMask->getBit(3),
      chipEnable,
      writeEnable};
  if (reset != nullptr) {
    clockInputs.push_back(reset);
  }
  SNLDesignModeling::addInputsToClockArcs(clockInputs, clock);
  SNLDesignModeling::addCombinatorialArcs(readAddressBits, readDataBits);

  SNLDesignModeling::MemoryInterface interface;
  interface.width = 4;
  interface.depth = 4;
  interface.abits = 2;
  interface.clock = clock;
  if (reset != nullptr) {
    interface.resetMode = SNLDesignModeling::MemoryResetMode::AsyncHigh;
    interface.reset = reset;
  }
  interface.readPorts.push_back(
      {.address = {address->getBit(0), address->getBit(1)},
       .data = {readData->getBit(0),
                readData->getBit(1),
                readData->getBit(2),
                readData->getBit(3)}});
  interface.writePorts.push_back(
      {.address = {address->getBit(0), address->getBit(1)},
       .data = {writeData->getBit(0),
                writeData->getBit(1),
                writeData->getBit(2),
                writeData->getBit(3)},
       .mask = {writeMask->getBit(0),
                writeMask->getBit(1),
                writeMask->getBit(2),
                writeMask->getBit(3)},
       .enables = {chipEnable, writeEnable}});
  SNLDesignModeling::setMemoryInterface(model, interface);
  return model;
}

SNLDesign* createSinglePortMemoryTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* memoryModel,
    std::optional<int> floatingWriteDataBit = std::nullopt,
    bool floatingReset = false) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topChipEnable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ce"));
  auto* topWriteEnable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("we"));
  auto* modelReset = memoryModel->getScalarTerm(NLName("RST"));
  auto* topReset = modelReset == nullptr || floatingReset
      ? nullptr
      : SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("rst"));
  auto* topAddress =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 1, 0, NLName("addr"));
  auto* topWriteData =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 3, 0, NLName("wdata"));
  auto* topWriteMask =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 3, 0, NLName("wmask"));
  auto* topOut =
      SNLBusTerm::create(top, SNLTerm::Direction::Output, 3, 0, NLName("out"));

  auto* memory = SNLInstance::create(top, memoryModel, NLName("mem0"));
  auto* clockNet = SNLScalarNet::create(top, NLName("clk_net"));
  auto* chipEnableNet = SNLScalarNet::create(top, NLName("ce_net"));
  auto* writeEnableNet = SNLScalarNet::create(top, NLName("we_net"));
  auto* resetNet = modelReset == nullptr
      ? nullptr
      : SNLScalarNet::create(top, NLName("rst_net"));
  auto* addressNet = SNLBusNet::create(top, 1, 0, NLName("addr_net"));
  auto* writeDataNet = SNLBusNet::create(top, 3, 0, NLName("wdata_net"));
  auto* writeMaskNet = SNLBusNet::create(top, 3, 0, NLName("wmask_net"));
  auto* outNet = SNLBusNet::create(top, 3, 0, NLName("out_net"));

  topClock->setNet(clockNet);
  topChipEnable->setNet(chipEnableNet);
  topWriteEnable->setNet(writeEnableNet);
  if (topReset != nullptr) {
    topReset->setNet(resetNet);
  }
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("CLK")))->setNet(clockNet);
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("CE")))->setNet(chipEnableNet);
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("WE")))->setNet(writeEnableNet);
  if (modelReset != nullptr) {
    memory->getInstTerm(modelReset)->setNet(resetNet);
  }

  auto* modelAddress = memoryModel->getBusTerm(NLName("ADDR"));
  auto* modelWriteData = memoryModel->getBusTerm(NLName("WDATA"));
  auto* modelWriteMask = memoryModel->getBusTerm(NLName("WMASK"));
  auto* modelReadData = memoryModel->getBusTerm(NLName("RDATA"));
  for (int bit = 0; bit <= 1; ++bit) {
    topAddress->getBit(bit)->setNet(addressNet->getBit(bit));
    memory->getInstTerm(modelAddress->getBit(bit))->setNet(addressNet->getBit(bit));
  }
  for (int bit = 0; bit <= 3; ++bit) {
    if (floatingWriteDataBit.has_value() && bit == *floatingWriteDataBit) {
      auto* floatingWriteDataNet = SNLScalarNet::create(
          top, NLName("floating_wdata" + std::to_string(bit) + "_net"));
      memory->getInstTerm(modelWriteData->getBit(bit))->setNet(floatingWriteDataNet);
    } else {
      topWriteData->getBit(bit)->setNet(writeDataNet->getBit(bit));
      memory->getInstTerm(modelWriteData->getBit(bit))->setNet(writeDataNet->getBit(bit));
    }
    topWriteMask->getBit(bit)->setNet(writeMaskNet->getBit(bit));
    topOut->getBit(bit)->setNet(outNet->getBit(bit));
    memory->getInstTerm(modelWriteMask->getBit(bit))->setNet(writeMaskNet->getBit(bit));
    memory->getInstTerm(modelReadData->getBit(bit))->setNet(outNet->getBit(bit));
  }

  return top;
}

std::filesystem::path repoRootForSecTests() {
  if (const char* prefix = std::getenv("TEST_DATA_PREFIX");
      prefix != nullptr) {
    return std::filesystem::path(prefix);
  }
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

struct Cva6SourceContextForSecTests {
  std::filesystem::path cva6RepoDir;
  std::filesystem::path hpdcacheDir;
  std::string targetCfg;
};

std::optional<Cva6SourceContextForSecTests> resolveCva6SourceContextForSecTests() {
  const char* cva6RepoDirEnv = std::getenv("CVA6_REPO_DIR");
  const char* hpdcacheDirEnv = std::getenv("HPDCACHE_DIR");
  const char* targetCfgEnv = std::getenv("TARGET_CFG");

  std::filesystem::path cva6RepoDir;
  if (cva6RepoDirEnv != nullptr && *cva6RepoDirEnv != '\0') {
    cva6RepoDir = cva6RepoDirEnv;
  } else {
    const auto fallbackRepoDir = std::filesystem::path("/Users/noamcohen/dev/CVA6/cva6");
    if (std::filesystem::exists(fallbackRepoDir)) {
      cva6RepoDir = fallbackRepoDir;
    }
  }
  if (cva6RepoDir.empty() || !std::filesystem::exists(cva6RepoDir)) {
    return std::nullopt;
  }

  std::filesystem::path hpdcacheDir;
  if (hpdcacheDirEnv != nullptr && *hpdcacheDirEnv != '\0') {
    hpdcacheDir = hpdcacheDirEnv;
  } else {
    const auto fallbackHpdcacheDir =
        cva6RepoDir / "core" / "cache_subsystem" / "hpdcache";
    if (std::filesystem::exists(fallbackHpdcacheDir)) {
      hpdcacheDir = fallbackHpdcacheDir;
    }
  }
  if (hpdcacheDir.empty() || !std::filesystem::exists(hpdcacheDir)) {
    return std::nullopt;
  }

  std::string targetCfg = "cv64a6_imafdc_sv39";
  if (targetCfgEnv != nullptr && *targetCfgEnv != '\0') {
    targetCfg = targetCfgEnv;
  }

  return Cva6SourceContextForSecTests{
      std::move(cva6RepoDir), std::move(hpdcacheDir), std::move(targetCfg)};
}

std::string substituteCva6FlistVariablesForSecTests(
    std::string text,
    const Cva6SourceContextForSecTests& context) {
  const std::array<std::pair<std::string_view, std::string>, 3> substitutions{{
      {"${CVA6_REPO_DIR}", context.cva6RepoDir.string()},
      {"${HPDCACHE_DIR}", context.hpdcacheDir.string()},
      {"${TARGET_CFG}", context.targetCfg},
  }};

  for (const auto& [needle, replacement] : substitutions) {
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != std::string::npos) {
      text.replace(pos, needle.size(), replacement);
      pos += replacement.size();
    }
  }
  return text;
}

void collectExpandedSlangArgsFromCommandFileForSecTests(
    const std::filesystem::path& commandFile,
    const Cva6SourceContextForSecTests& context,
    std::unordered_set<std::string>& visitedFiles,
    std::vector<std::string>& args) {
  const auto normalizedPath = std::filesystem::weakly_canonical(commandFile);
  if (!visitedFiles.insert(normalizedPath.string()).second) {
    return;
  }

  std::ifstream input(normalizedPath);
  ASSERT_TRUE(input.good()) << "Failed to read command file: " << normalizedPath.string();

  std::string line;
  while (std::getline(input, line)) {
    line = substituteCva6FlistVariablesForSecTests(line, context);
    const auto commentPos = line.find("//");
    if (commentPos != std::string::npos) {
      line.erase(commentPos);
    }
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
      continue;
    }
    const auto last = line.find_last_not_of(" \t\r\n");
    line = line.substr(first, last - first + 1);
    if (line.empty()) {
      continue;
    }

    if (line.rfind("-F ", 0) == 0 || line.rfind("-f ", 0) == 0) {
      const auto nestedPath = normalizedPath.parent_path() / line.substr(3);
      collectExpandedSlangArgsFromCommandFileForSecTests(
          nestedPath, context, visitedFiles, args);
      continue;
    }
    args.push_back(std::move(line));
  }
}

SNLSVConstructor::Paths buildExpandedCva6SlangArgsForSecTests(
    const Cva6SourceContextForSecTests& context,
    const std::string& topName,
    const std::vector<std::filesystem::path>& extraSources = {}) {
  const auto flistPath = context.cva6RepoDir / "core" / "Flist.cva6";
  std::unordered_set<std::string> visitedFiles;
  std::vector<std::string> args;
  collectExpandedSlangArgsFromCommandFileForSecTests(
      flistPath, context, visitedFiles, args);
  for (const auto& extraSource : extraSources) {
    args.push_back(extraSource.string());
  }
  args.push_back("--top");
  args.push_back(topName);

  SNLSVConstructor::Paths paths;
  paths.reserve(args.size());
  for (const auto& arg : args) {
    paths.emplace_back(arg);
  }
  return paths;
}

SNLDesign* loadLibertyMemoryModel(
    NLLibrary* primitivesLibrary,
    const std::string& libertyFileName,
    const std::string& cellName) {
  SNLLibertyConstructor constructor(primitivesLibrary);
  constructor.construct(repoRootForSecTests() / "example" / libertyFileName);
  auto* model = primitivesLibrary->getSNLDesign(NLName(cellName));
  if (model == nullptr) {
    throw std::runtime_error("Failed to load Liberty memory model `" + cellName + "`");
  }
  return model;
}

SNLDesign* loadSystemVerilogTopFromSource(
    NLLibrary* designLibrary,
    const std::string& moduleName,
    const std::string& sourceText) {
  const auto svDir =
      std::filesystem::temp_directory_path() / ("sec_sv_" + moduleName);
  std::filesystem::remove_all(svDir);
  std::filesystem::create_directories(svDir);

  const auto svPath = svDir / (moduleName + ".sv");
  std::ofstream svFile(svPath);
  if (!svFile.good()) {
    throw std::runtime_error(
        "Failed to create temporary SystemVerilog source `" +
        svPath.string() + "`");
  }
  svFile << sourceText;
  svFile.close();

  SNLSVConstructor constructor(designLibrary);
  constructor.construct(svPath);
  auto* top = designLibrary->getSNLDesign(NLName(moduleName));
  if (top == nullptr) {
    throw std::runtime_error(
        "Failed to construct SystemVerilog top `" + moduleName + "`");
  }
  return top;
}

SNLDesign* loadSystemVerilogTopFromPaths(
    NLLibrary* designLibrary,
    const std::string& moduleName,
    const SNLSVConstructor::Paths& paths) {
  SNLSVConstructor constructor(designLibrary);
  constructor.construct(paths);
  auto* top = designLibrary->getSNLDesign(NLName(moduleName));
  if (top == nullptr) {
    throw std::runtime_error(
        "Failed to construct SystemVerilog top `" + moduleName + "`");
  }
  return top;
}

SNLDesign* loadRealCva6PerfCountersTargetConfigTopForSecTests(
    NLLibrary* designLibrary,
    const Cva6SourceContextForSecTests& context,
    const std::string& moduleName) {
  const auto svDir = std::filesystem::temp_directory_path() / moduleName;
  std::filesystem::remove_all(svDir);
  std::filesystem::create_directories(svDir);
  const auto wrapperPath = svDir / (moduleName + ".sv");
  std::ofstream wrapperFile(wrapperPath);
  if (!wrapperFile.good()) {
    throw std::runtime_error(
        "Failed to create SEC CVA6 wrapper `" + wrapperPath.string() + "`");
  }
  wrapperFile
      << R"(module )"
      << moduleName
      << R"(
  import ariane_pkg::*;
  import cva6_config_pkg::*;
#(
  parameter config_pkg::cva6_cfg_t CVA6Cfg =
      build_config_pkg::build_config(cva6_cfg)
) ();
  localparam type branchpredict_sbe_t = struct packed {
    cf_t                     cf;
    logic [CVA6Cfg.VLEN-1:0] predict_address;
  };

  localparam type exception_t = struct packed {
    logic [CVA6Cfg.XLEN-1:0]  cause;
    logic [CVA6Cfg.XLEN-1:0]  tval;
    logic [CVA6Cfg.GPLEN-1:0] tval2;
    logic [31:0]              tinst;
    logic                     gva;
    logic                     valid;
  };

  localparam type bp_resolve_t = struct packed {
    logic                    valid;
    logic [CVA6Cfg.VLEN-1:0] pc;
    logic [CVA6Cfg.VLEN-1:0] target_address;
    logic                    is_mispredict;
    logic                    is_taken;
    cf_t                     cf_type;
  };

  localparam type icache_dreq_t = struct packed {
    logic                    req;
    logic                    kill_s1;
    logic                    kill_s2;
    logic                    spec;
    logic [CVA6Cfg.VLEN-1:0] vaddr;
  };

  localparam type cbo_t = logic [7:0];

  localparam type dcache_req_i_t = struct packed {
    logic [CVA6Cfg.DCACHE_INDEX_WIDTH-1:0] address_index;
    logic [CVA6Cfg.DCACHE_TAG_WIDTH-1:0]   address_tag;
    logic [CVA6Cfg.XLEN-1:0]               data_wdata;
    logic [CVA6Cfg.DCACHE_USER_WIDTH-1:0]  data_wuser;
    logic                                  data_req;
    logic                                  data_we;
    logic [(CVA6Cfg.XLEN/8)-1:0]           data_be;
    logic [1:0]                            data_size;
    logic [CVA6Cfg.DcacheIdWidth-1:0]      data_id;
    logic                                  kill_req;
    logic                                  tag_valid;
    cbo_t                                  cbo_op;
  };

  localparam type scoreboard_entry_t = struct packed {
    logic [CVA6Cfg.VLEN-1:0]              pc;
    logic [CVA6Cfg.TRANS_ID_BITS-1:0]     trans_id;
    fu_t                                  fu;
    fu_op                                 op;
    logic [REG_ADDR_SIZE-1:0]             rs1;
    logic [REG_ADDR_SIZE-1:0]             rs2;
    logic [REG_ADDR_SIZE-1:0]             rd;
    logic [CVA6Cfg.XLEN-1:0]              result;
    logic                                 valid;
    logic                                 use_imm;
    logic                                 use_zimm;
    logic                                 use_pc;
    exception_t                           ex;
    branchpredict_sbe_t                   bp;
    logic                                 is_compressed;
    logic                                 is_macro_instr;
    logic                                 is_last_macro_instr;
    logic                                 is_double_rd_macro_instr;
    logic                                 vfp;
    logic                                 is_zcmt;
  };

  logic clk_i;
  logic rst_ni;
  logic debug_mode_i;
  logic [11:0] addr_i;
  logic we_i;
  logic [CVA6Cfg.XLEN-1:0] data_i;
  logic [CVA6Cfg.XLEN-1:0] data_o;
  scoreboard_entry_t [CVA6Cfg.NrCommitPorts-1:0] commit_instr_i;
  logic [CVA6Cfg.NrCommitPorts-1:0] commit_ack_i;
  logic l1_icache_miss_i;
  logic l1_dcache_miss_i;
  logic itlb_miss_i;
  logic dtlb_miss_i;
  logic sb_full_i;
  logic if_empty_i;
  exception_t ex_i;
  logic eret_i;
  bp_resolve_t resolved_branch_i;
  exception_t branch_exceptions_i;
  icache_dreq_t l1_icache_access_i;
  dcache_req_i_t [2:0] l1_dcache_access_i;
  logic [2:0][CVA6Cfg.DCACHE_SET_ASSOC-1:0] miss_vld_bits_i;
  logic i_tlb_flush_i;
  logic stall_issue_i;
  logic [31:0] mcountinhibit_i;

  assign clk_i = 1'b0;
  assign rst_ni = 1'b1;
  assign debug_mode_i = 1'b0;
  assign addr_i = '0;
  assign we_i = 1'b0;
  assign data_i = '0;
  assign commit_instr_i = '0;
  assign commit_ack_i = '0;
  assign l1_icache_miss_i = 1'b0;
  assign l1_dcache_miss_i = 1'b0;
  assign itlb_miss_i = 1'b0;
  assign dtlb_miss_i = 1'b0;
  assign sb_full_i = 1'b0;
  assign if_empty_i = 1'b0;
  assign ex_i = '0;
  assign eret_i = 1'b0;
  assign resolved_branch_i = '0;
  assign branch_exceptions_i = '0;
  assign l1_icache_access_i = '0;
  assign l1_dcache_access_i = '0;
  assign miss_vld_bits_i = '0;
  assign i_tlb_flush_i = 1'b0;
  assign stall_issue_i = 1'b0;
  assign mcountinhibit_i = '0;

  perf_counters #(
    .CVA6Cfg(CVA6Cfg),
    .bp_resolve_t(bp_resolve_t),
    .dcache_req_i_t(dcache_req_i_t),
    .dcache_req_o_t(dcache_req_i_t),
    .exception_t(exception_t),
    .icache_dreq_t(icache_dreq_t),
    .scoreboard_entry_t(scoreboard_entry_t)
  ) dut (
    .clk_i(clk_i),
    .rst_ni(rst_ni),
    .debug_mode_i(debug_mode_i),
    .addr_i(addr_i),
    .we_i(we_i),
    .data_i(data_i),
    .data_o(data_o),
    .commit_instr_i(commit_instr_i),
    .commit_ack_i(commit_ack_i),
    .l1_icache_miss_i(l1_icache_miss_i),
    .l1_dcache_miss_i(l1_dcache_miss_i),
    .itlb_miss_i(itlb_miss_i),
    .dtlb_miss_i(dtlb_miss_i),
    .sb_full_i(sb_full_i),
    .if_empty_i(if_empty_i),
    .ex_i(ex_i),
    .eret_i(eret_i),
    .resolved_branch_i(resolved_branch_i),
    .branch_exceptions_i(branch_exceptions_i),
    .l1_icache_access_i(l1_icache_access_i),
    .l1_dcache_access_i(l1_dcache_access_i),
    .miss_vld_bits_i(miss_vld_bits_i),
    .i_tlb_flush_i(i_tlb_flush_i),
    .stall_issue_i(stall_issue_i),
    .mcountinhibit_i(mcountinhibit_i)
  );
endmodule
)";
  wrapperFile.close();

  const auto args = buildExpandedCva6SlangArgsForSecTests(
      context, moduleName, {wrapperPath});
  return loadSystemVerilogTopFromPaths(designLibrary, moduleName, args);
}

SNLDesign* createMirroredInstanceTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* model) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* instance = SNLInstance::create(top, model, NLName("mem0"));

  for (auto* scalarTerm : model->getScalarTerms()) {
    auto* topTerm = SNLScalarTerm::create(
        top, scalarTerm->getDirection(), scalarTerm->getName());
    auto* net = SNLScalarNet::create(
        top, NLName(scalarTerm->getName().getString() + "_net"));
    topTerm->setNet(net);
    instance->getInstTerm(scalarTerm)->setNet(net);
  }

  for (auto* busTerm : model->getBusTerms()) {
    auto* topTerm = SNLBusTerm::create(
        top,
        busTerm->getDirection(),
        busTerm->getMSB(),
        busTerm->getLSB(),
        busTerm->getName());
    auto* net = SNLBusNet::create(
        top,
        busTerm->getMSB(),
        busTerm->getLSB(),
        NLName(busTerm->getName().getString() + "_net"));
    for (int bit = busTerm->getLSB(); bit <= busTerm->getMSB(); ++bit) {
      topTerm->getBit(bit)->setNet(net->getBit(bit));
      instance->getInstTerm(busTerm->getBit(bit))->setNet(net->getBit(bit));
    }
  }

  return top;
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

SNLDesign* createOpaqueBoundaryTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* opaqueModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* opaque = SNLInstance::create(top, opaqueModel, NLName("opaque0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topOut->setNet(netOut);
  opaque->getInstTerm(opaqueModel->getScalarTerm(NLName("A")))->setNet(netIn);
  opaque->getInstTerm(opaqueModel->getScalarTerm(NLName("Y")))->setNet(netOut);

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

SNLDesign* createPartialCoverageNoDriverDataConeTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* andModel) {
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
  auto* andGate = SNLInstance::create(top, andModel, NLName("data_gate"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netFloating = SNLScalarNet::create(top, NLName("net_floating"));
  auto* netData = SNLScalarNet::create(top, NLName("net_data"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topGood->setNet(netIn);
  topBad->setNet(netQ);

  andGate->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netIn);
  andGate->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netFloating);
  andGate->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netData);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netData);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createPartialCoverageDrivenTop(
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
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topGood->setNet(netIn);
  topBad->setNet(netQ);

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
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

SNLDesign* createPartialCoverageLogicalLoopTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topSel =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("sel"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topGood =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("good"));
  auto* topBad =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("bad"));

  auto* assign = SNLInstance::create(top, NLDB0::getAssign(), NLName("assign0"));
  auto* mux = SNLInstance::create(top, NLDB0::getMux2(), NLName("mux0"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netSel = SNLScalarNet::create(top, NLName("net_sel"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLoopSeed = SNLScalarNet::create(top, NLName("net_loop_seed"));
  auto* netLoopIn = SNLScalarNet::create(top, NLName("net_loop_in"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topSel->setNet(netSel);
  topClock->setNet(netClock);
  topGood->setNet(netIn);
  topBad->setNet(netQ);

  assign->getInstTerm(NLDB0::getAssignInput())->setNet(netLoopIn);
  assign->getInstTerm(NLDB0::getAssignOutput())->setNet(netLoopSeed);

  mux->getInstTerm(NLDB0::getMux2InputA()->getBit(0))->setNet(netLoopSeed);
  mux->getInstTerm(NLDB0::getMux2InputB()->getBit(0))->setNet(netIn);
  mux->getInstTerm(NLDB0::getMux2Select())->setNet(netSel);
  mux->getInstTerm(NLDB0::getMux2Output()->getBit(0))->setNet(netLoopIn);

  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netLoopSeed);
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

SNLDesign* createOpaqueClockGateLatchModel(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("DLATCH_N"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("D"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Input, NLName("GATE"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("Q"));
  return model;
}

SNLDesign* createConstantLowModel(NLLibrary* library) {
  auto* model =
      SNLDesign::create(library, SNLDesign::Type::Primitive, NLName("CONB"));
  SNLScalarTerm::create(model, SNLTerm::Direction::Output, NLName("LO"));
  SNLDesignModeling::setTruthTable(
      model, SNLTruthTable(0, 0, SNLTruthTable::fullDependencies(0)));
  return model;
}

SNLDesign* createClockGateLatchDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* andModel,
    SNLDesign* latchModel) {
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

  auto* latch = SNLInstance::create(top, latchModel, NLName("clock_gate_i.en_latch"));
  auto* gateAnd = SNLInstance::create(top, andModel, NLName("clock_gate_i.and_clk"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netEnable = SNLScalarNet::create(top, NLName("net_en"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLatchQ = SNLScalarNet::create(top, NLName("net_latch_q"));
  auto* netGatedClock = SNLScalarNet::create(top, NLName("net_gated_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topEnable->setNet(netEnable);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  latch->getInstTerm(latchModel->getScalarTerm(NLName("D")))->setNet(netEnable);
  latch->getInstTerm(latchModel->getScalarTerm(NLName("GATE")))->setNet(netClock);
  latch->getInstTerm(latchModel->getScalarTerm(NLName("Q")))->setNet(netLatchQ);

  gateAnd->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netClock);
  gateAnd->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netLatchQ);
  gateAnd->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netGatedClock);

  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netGatedClock);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netOut);

  return top;
}

SNLDesign* createClockTreeBufferedDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* bufferModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* rootBuffer =
      SNLInstance::create(top, bufferModel, NLName("wire4069"));
  auto* clockBuffer =
      SNLInstance::create(top, bufferModel, NLName("clkbuf_leaf_0_clk"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netClockRoot = SNLScalarNet::create(top, NLName("net4068"));
  auto* netLeafClock = SNLScalarNet::create(top, NLName("clknet_leaf_0_clk"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  rootBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("A")))->setNet(
      netClock);
  rootBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("Y")))->setNet(
      netClockRoot);

  clockBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("A")))->setNet(
      netClockRoot);
  clockBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("Y")))->setNet(
      netLeafClock);

  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netLeafClock);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netOut);

  return top;
}

SNLDesign* createDataBufferedDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* bufferModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* dataBuffer =
      SNLInstance::create(top, bufferModel, NLName("data_buffer"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netData = SNLScalarNet::create(top, NLName("net_data"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topOut->setNet(netOut);

  dataBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("A")))->setNet(
      netIn);
  dataBuffer->getInstTerm(bufferModel->getScalarTerm(NLName("Y")))->setNet(
      netData);

  ff->getInstTerm(NLDB0::getDFFData())->setNet(netData);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netOut);

  return top;
}

SNLDesign* createInvertedClockDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* invModel,
    const std::string& invInstanceName = "inv0") {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* inv = SNLInstance::create(top, invModel, NLName(invInstanceName));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netClockN = SNLScalarNet::create(top, NLName("net_clk_n"));
  auto* netQ = SNLScalarNet::create(top, NLName("net_q"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topOut->setNet(netQ);

  inv->getInstTerm(invModel->getScalarTerm(NLName("A")))->setNet(netClock);
  inv->getInstTerm(invModel->getScalarTerm(NLName("Y")))->setNet(netClockN);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClockN);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netQ);

  return top;
}

SNLDesign* createPosToNegSameDomainTop(
    NLLibrary* library,
    const std::string& name) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topIn =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* posFf = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff_pos"));
  auto* negFf = SNLInstance::create(top, NLDB0::getDFFN(), NLName("ff_neg"));

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netPosQ = SNLScalarNet::create(top, NLName("net_pos_q"));
  auto* netNegQ = SNLScalarNet::create(top, NLName("net_neg_q"));

  topIn->setNet(netIn);
  topClock->setNet(netClock);
  topOut->setNet(netNegQ);

  posFf->getInstTerm(NLDB0::getDFFData())->setNet(netIn);
  posFf->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  posFf->getInstTerm(NLDB0::getDFFOutput())->setNet(netPosQ);

  negFf->getInstTerm(NLDB0::getDFFNData())->setNet(netPosQ);
  negFf->getInstTerm(NLDB0::getDFFNClock())->setNet(netClock);
  negFf->getInstTerm(NLDB0::getDFFNOutput())->setNet(netNegQ);

  return top;
}

SNLDesign* createMultiClockDomainOutputTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* andModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topInA =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in_a"));
  auto* topInB =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("in_b"));
  auto* topClockA =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("a_clk"));
  auto* topClockB =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("b_clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* ffA = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff_a"));
  auto* ffB = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff_b"));
  auto* andInst = SNLInstance::create(top, andModel, NLName("and_domains"));

  auto* netInA = SNLScalarNet::create(top, NLName("net_in_a"));
  auto* netInB = SNLScalarNet::create(top, NLName("net_in_b"));
  auto* netClockA = SNLScalarNet::create(top, NLName("net_a_clk"));
  auto* netClockB = SNLScalarNet::create(top, NLName("net_b_clk"));
  auto* netQa = SNLScalarNet::create(top, NLName("net_qa"));
  auto* netQb = SNLScalarNet::create(top, NLName("net_qb"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topInA->setNet(netInA);
  topInB->setNet(netInB);
  topClockA->setNet(netClockA);
  topClockB->setNet(netClockB);
  topOut->setNet(netOut);

  ffA->getInstTerm(NLDB0::getDFFData())->setNet(netInA);
  ffA->getInstTerm(NLDB0::getDFFClock())->setNet(netClockA);
  ffA->getInstTerm(NLDB0::getDFFOutput())->setNet(netQa);
  ffB->getInstTerm(NLDB0::getDFFData())->setNet(netInB);
  ffB->getInstTerm(NLDB0::getDFFClock())->setNet(netClockB);
  ffB->getInstTerm(NLDB0::getDFFOutput())->setNet(netQb);
  andInst->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netQa);
  andInst->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netQb);
  andInst->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netOut);

  return top;
}

SNLDesign* createClockGateLatchDataDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* andModel,
    SNLDesign* latchModel,
    bool includeIndependentDff = false) {
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
  SNLScalarTerm* topIndependentIn = nullptr;
  SNLScalarTerm* topIndependentOut = nullptr;
  if (includeIndependentDff) {
    topIndependentIn = SNLScalarTerm::create(
        top, SNLTerm::Direction::Input, NLName("independent_in"));
    topIndependentOut = SNLScalarTerm::create(
        top, SNLTerm::Direction::Output, NLName("independent_out"));
  }

  auto* latch = SNLInstance::create(top, latchModel, NLName("clock_gate_i.en_latch"));
  auto* dataAnd = SNLInstance::create(top, andModel, NLName("data_and"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));
  SNLInstance* independentFf = nullptr;
  if (includeIndependentDff) {
    independentFf = SNLInstance::create(
        top, NLDB0::getDFF(), NLName("independent_ff"));
  }

  auto* netIn = SNLScalarNet::create(top, NLName("net_in"));
  auto* netEnable = SNLScalarNet::create(top, NLName("net_en"));
  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netLatchQ = SNLScalarNet::create(top, NLName("net_latch_q"));
  auto* netData = SNLScalarNet::create(top, NLName("net_data"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));
  SNLScalarNet* netIndependentIn = nullptr;
  SNLScalarNet* netIndependentOut = nullptr;
  if (includeIndependentDff) {
    netIndependentIn = SNLScalarNet::create(top, NLName("net_independent_in"));
    netIndependentOut =
        SNLScalarNet::create(top, NLName("net_independent_out"));
  }

  topIn->setNet(netIn);
  topEnable->setNet(netEnable);
  topClock->setNet(netClock);
  topOut->setNet(netOut);
  if (includeIndependentDff) {
    topIndependentIn->setNet(netIndependentIn);
    topIndependentOut->setNet(netIndependentOut);
  }

  latch->getInstTerm(latchModel->getScalarTerm(NLName("D")))->setNet(netEnable);
  latch->getInstTerm(latchModel->getScalarTerm(NLName("GATE")))->setNet(netClock);
  latch->getInstTerm(latchModel->getScalarTerm(NLName("Q")))->setNet(netLatchQ);

  dataAnd->getInstTerm(andModel->getScalarTerm(NLName("A")))->setNet(netIn);
  dataAnd->getInstTerm(andModel->getScalarTerm(NLName("B")))->setNet(netLatchQ);
  dataAnd->getInstTerm(andModel->getScalarTerm(NLName("Y")))->setNet(netData);

  ff->getInstTerm(NLDB0::getDFFData())->setNet(netData);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netOut);
  if (includeIndependentDff) {
    // This cone intentionally does not reference the folded latch output. It
    // catches regressions where latch substitution rebuilds unrelated SEC
    // state expressions instead of preserving no-op subtrees.
    independentFf->getInstTerm(NLDB0::getDFFData())->setNet(netIndependentIn);
    independentFf->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
    independentFf->getInstTerm(NLDB0::getDFFOutput())->setNet(netIndependentOut);
  }

  return top;
}

SNLDesign* createConstantDrivenDffTop(
    NLLibrary* library,
    const std::string& name,
    SNLDesign* constantModel) {
  auto* top =
      SNLDesign::create(library, SNLDesign::Type::Standard, NLName(name));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topOut =
      SNLScalarTerm::create(top, SNLTerm::Direction::Output, NLName("out"));

  auto* constant = SNLInstance::create(top, constantModel, NLName("tie0"));
  auto* ff = SNLInstance::create(top, NLDB0::getDFF(), NLName("ff0"));

  auto* netClock = SNLScalarNet::create(top, NLName("net_clk"));
  auto* netConstant = SNLScalarNet::create(top, NLName("net_const"));
  auto* netOut = SNLScalarNet::create(top, NLName("net_out"));

  topClock->setNet(netClock);
  topOut->setNet(netOut);

  constant->getInstTerm(constantModel->getScalarTerm(NLName("LO")))->setNet(netConstant);
  ff->getInstTerm(NLDB0::getDFFData())->setNet(netConstant);
  ff->getInstTerm(NLDB0::getDFFClock())->setNet(netClock);
  ff->getInstTerm(NLDB0::getDFFOutput())->setNet(netOut);

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

void expectAllExpressionSupportIsPublished(const SequentialDesignModel& model) {
  std::unordered_set<size_t> publishedVars;
  for (const auto& key : model.environmentInputs) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      publishedVars.insert(varIt->second);
    }
  }
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      publishedVars.insert(varIt->second);
    }
  }

  auto checkExpr = [&](const char* expressionKind,
                       const SignalKey& key,
                       BoolExpr* expr) {
    ASSERT_NE(expr, nullptr);
    const auto nameIt = model.displayNameByKey.find(key);
    const std::string displayName =
        nameIt == model.displayNameByKey.end() ? signalKeyToString(key)
                                               : nameIt->second;
    for (const auto varID : expr->getSupportVars()) {
      if (varID < 2) {
        continue;
      }
      EXPECT_NE(publishedVars.find(varID), publishedVars.end())
          << expressionKind << " `" << displayName
          << "` references unpublished variable " << varID;
    }
  };

  for (const auto& [key, expr] : model.observedOutputExprByKey) {
    checkExpr("observed output", key, expr);
  }
  for (const auto& [key, expr] : model.nextStateExprByStateKey) {
    checkExpr("next-state expression", key, expr);
  }
}

}  // namespace

TEST_F(SequentialEquivalenceStrategyTests,
       EmitSecDiagIsQuietWithoutDiagnosticEnvironment) {
  const ScopedUnsetEnvVar secDiag("KEPLER_SEC_DIAG");
  const ScopedUnsetEnvVar kiDiag("KEPLER_SEC_KI_DIAG");
  const ScopedUnsetEnvVar pdrStats("KEPLER_SEC_PDR_STATS");
  const ScopedUnsetEnvVar pdrTrace("KEPLER_SEC_PDR_TRACE");
  const ScopedUnsetEnvVar summaryStats("KEPLER_SEC_SUMMARY_STATS");

  testing::internal::CaptureStderr();
  emitSecDiag("SEC diag: should stay quiet");
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_TRUE(stderrOutput.empty()) << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       EmitSecDiagWritesWithDiagnosticEnvironment) {
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  emitSecDiag("SEC diag: visible ", 42);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(stderrOutput.find("SEC diag: visible 42"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecEngineProofProgressListsUnprovenOutputs) {
  const SequentialEquivalenceProofProgress progress =
      detail::buildSecEngineProofProgress(
          "IMC",
          {"wide_out0", "wide_out1", "wide_out2"},
          /*totalOutputCount=*/3,
          /*provenOutputCount=*/1);
  const std::vector<std::string> lines =
      detail::buildSecEngineProofProgressDiagLines(
          "IMC",
          {"wide_out0", "wide_out1", "wide_out2"},
          /*totalOutputCount=*/3,
          /*provenOutputCount=*/1);

  EXPECT_EQ(progress.engineLabel, "IMC");
  EXPECT_EQ(progress.provenOutputs, 1u);
  EXPECT_EQ(progress.totalOutputs, 3u);
  ASSERT_EQ(progress.unprovenOutputs.size(), 2u);
  EXPECT_EQ(progress.unprovenOutputs[0].index, 1u);
  EXPECT_EQ(progress.unprovenOutputs[0].name, "wide_out1");
  EXPECT_EQ(progress.unprovenOutputs[1].index, 2u);
  EXPECT_EQ(progress.unprovenOutputs[1].name, "wide_out2");
  ASSERT_EQ(lines.size(), 3u);
  EXPECT_EQ(lines[0], "SEC diag: SEC IMC proven outputs: 1/3");
  EXPECT_EQ(
      lines[1], "SEC diag: SEC IMC not proven output[1]=wide_out1");
  EXPECT_EQ(
      lines[2], "SEC diag: SEC IMC not proven output[2]=wide_out2");
}

TEST_F(SequentialEquivalenceStrategyTests,
       BinarySecSkipsUninitializedDffOutput) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library = NLLibrary::create(db, NLName("LIB"));
  auto* primitives = NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("PRIMS"));
  auto* invModel = createInvModel(primitives);

  auto* top0 =
      createDffTop(library, "top0", invModel, false, false, "in", "out", "ff0");
  auto* top1 =
      createDffTop(library, "top1", invModel, false, false, "in", "out", "ff1");

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::Pdr);
  const auto result = strategy.run(2);

  // Binary SEC cannot compare independently initialized internal state. The
  // state-dependent output is skipped instead of reporting a false mismatch.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, 0u);
}



TEST_F(SequentialEquivalenceStrategyTests,
       BoolFormulaImplicationProvesCommutedConeUnderStateEquality) {
  BoolExpr* stateEquality = makeEqualityExpr(BoolExpr::Var(2), BoolExpr::Var(4));
  BoolExpr* outputEquality = makeEqualityExpr(
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
      BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)));
  BoolExpr* unrelatedEquality =
      makeEqualityExpr(BoolExpr::Var(2), BoolExpr::Var(3));

  EXPECT_TRUE(boolFormulaImplies(
      stateEquality,
      outputEquality,
      KEPLER_FORMAL::Config::getSolverType()));
  EXPECT_FALSE(boolFormulaImplies(
      stateEquality,
      unrelatedEquality,
      KEPLER_FORMAL::Config::getSolverType()));
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
  EXPECT_FALSE(isBoolFormulaSatisfiable(
      nullptr, KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_FALSE(isBoolFormulaSatisfiable(
      BoolExpr::createFalse(), KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_TRUE(isBoolFormulaSatisfiable(
      BoolExpr::createTrue(), KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_FALSE(boolFormulaImplies(
      BoolExpr::createTrue(), nullptr, KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_TRUE(boolFormulaImplies(
      BoolExpr::createFalse(),
      BoolExpr::Var(2),
      KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_TRUE(boolFormulaImplies(
      BoolExpr::Var(2),
      BoolExpr::createTrue(),
      KEPLER_FORMAL::Config::SolverType::KISSAT));

  BoolExpr invalid;
  EXPECT_THROW(
      static_cast<void>(remapBoolExprVariables(&invalid, {})),
      std::runtime_error);

  auto* xorExpr = BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::Var(3));
  auto* substituted =
      substituteBoolExprVariables(xorExpr, {{2, true}, {3, false}});
  EXPECT_TRUE(substituted->evaluate({}));

  KEPLER_FORMAL::BoolExprCache::Key rawAndWithFalse{
      KEPLER_FORMAL::Op::AND, 0, BoolExpr::createFalse(), BoolExpr::Var(4)};
  EXPECT_FALSE(isBoolFormulaSatisfiable(
      KEPLER_FORMAL::BoolExprCache::getExpression(rawAndWithFalse),
      KEPLER_FORMAL::Config::SolverType::KISSAT));

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
       BoolExprSubstitutionPreservesUnchangedBootstrapCones) {
  auto* preservedDataCone =
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::Not(BoolExpr::Var(3)));
  auto* expr = BoolExpr::Or(preservedDataCone, BoolExpr::Var(4));

  // Reset/bootstrap specialization may visit large ASIC data cones with an
  // assignment frontier that does not touch them. Keeping those nodes by pointer
  // prevents a runtime regression where BlackParrot rebuilt equivalent cones.
  EXPECT_EQ(substituteBoolExprVariables(expr, {{99, true}}), expr);
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

TEST_F(SequentialEquivalenceStrategyTests, SatEncodingFlatCacheGrows) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  std::unordered_map<size_t, int> leafLits;
  leafLits.emplace(2, solver.newVar() + 2);
  leafLits.emplace(3, solver.newVar() + 2);

  BoolExpr* expr = BoolExpr::Var(2);
  for (size_t index = 0; index < 1500; ++index) {
    expr = BoolExpr::Xor(expr, BoolExpr::Var(3));
  }

  FrameFormulaEncoder encoder(solver, std::move(leafLits));
  EXPECT_NE(encoder.encode(expr), 0);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SatEncodingCachedPostorderMatchesNormalDagEncoding) {
  constexpr size_t a = 2;
  constexpr size_t b = 3;
  constexpr size_t c = 4;
  constexpr size_t state = 5;
  constexpr size_t secondState = 6;
  BoolExpr* shared = BoolExpr::And(BoolExpr::Var(a), BoolExpr::Var(b));
  BoolExpr* transition =
      BoolExpr::Or(shared, BoolExpr::Xor(shared, BoolExpr::Var(c)));
  KInductionProblem problem;
  problem.state0Symbols = {state, secondState};
  problem.inputSymbols = {a, b, c};
  problem.allSymbols = {a, b, c, state, secondState};
  problem.transitions0 = {{state, transition}, {secondState, transition}};

  const TransitionExprResolver resolver(problem);
  const std::vector<BoolExpr*>& postorder = resolver.encodingPostorder(state);
  EXPECT_EQ(postorder,
            (std::vector<BoolExpr*>{BoolExpr::Var(a), BoolExpr::Var(b), shared,
                                    BoolExpr::Var(c), transition->getRight(),
                                    transition}));
  EXPECT_EQ(&postorder, &resolver.encodingPostorder(state));
  EXPECT_EQ(&postorder, &resolver.encodingPostorder(secondState));

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  FrameVariableStore variables(solver, {a, b, c}, 1);
  const auto leafLits = variables.makeLeafLits(0);
  FrameFormulaEncoder normalEncoder(solver, leafLits);
  FrameFormulaEncoder plannedEncoder(solver, leafLits);
  const int normalRoot = normalEncoder.encode(transition);
  const int plannedRoot = plannedEncoder.encode(transition, postorder);

  // Both encoders share the same leaves but allocate independent Tseitin
  // literals. Opposite root assumptions must therefore be UNSAT in both
  // directions, proving that the cached traversal changes preparation only.
  EXPECT_EQ(solver.solveWithAssumptionsStatus({normalRoot, -plannedRoot}),
            SATSolverWrapper::SolveStatus::Unsat);
  EXPECT_EQ(solver.solveWithAssumptionsStatus({-normalRoot, plannedRoot}),
            SATSolverWrapper::SolveStatus::Unsat);

  SATSolverWrapper invalidSolver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  FrameVariableStore invalidVariables(invalidSolver, {a, b, c}, 1);
  FrameFormulaEncoder invalidEncoder(invalidSolver,
                                     invalidVariables.makeLeafLits(0));
  EXPECT_THROW(
      static_cast<void>(invalidEncoder.encode(transition, {transition})),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SatEncodingHonorsLargeHintedFrameFormulaReserve) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  std::unordered_map<size_t, int> leafLits;
  leafLits.emplace(2, solver.newVar() + 2);
  leafLits.emplace(3, solver.newVar() + 2);

  BoolExpr* expr = BoolExpr::Var(2);
  constexpr size_t kDepthAboveOldSolverReserve = 70000;
  for (size_t index = 0; index < kDepthAboveOldSolverReserve; ++index) {
    expr = BoolExpr::Xor(expr, BoolExpr::Var(3));
  }

  // BP-sized strict KI frames provide an exact transition-DAG hint. Encoding a
  // cone above the old 64K solver reserve keeps this optimization covered
  // without changing any Tseitin clauses or SAT result.
  FrameFormulaEncoder encoder(
      solver, std::move(leafLits), kDepthAboveOldSolverReserve + 1);
  EXPECT_NE(encoder.encode(expr), 0);
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
       KissatResourceLimitedSolveReportsUnknownInsteadOfUnsat) {
  SATSolverWrapper limitedSolver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  limitedSolver.configureForSecResetExpressionProof();
  const int x = limitedSolver.newVar() + 2;
  const int y = limitedSolver.newVar() + 2;
  limitedSolver.addClause({x, y});
  limitedSolver.addClause({-x, y});

  // Optional PDR shortcuts may cap Kissat work. A limit hit must be observable
  // as UNKNOWN so callers do not accidentally learn a bogus UNSAT cube.
  EXPECT_EQ(
      limitedSolver.solveWithKissatResourceLimits(
          std::numeric_limits<unsigned>::max(),
          /*decisionLimit=*/0),
      SATSolverWrapper::SolveStatus::Unknown);

  SATSolverWrapper unboundedSolver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  unboundedSolver.configureForSecResetExpressionProof();
  const int ux = unboundedSolver.newVar() + 2;
  const int uy = unboundedSolver.newVar() + 2;
  unboundedSolver.addClause({ux, uy});
  unboundedSolver.addClause({-ux, uy});
  EXPECT_EQ(unboundedSolver.solveStatus(), SATSolverWrapper::SolveStatus::Sat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalResourceLimitedSolveReportsUnknownOnDecisionBudget) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecPdrQuery();
  const int x = solver.newVar() + 2;
  const int y = solver.newVar() + 2;
  solver.addClause({x, y});
  solver.addClause({-x, y});

  // Dual-rail PDR bad-cube repair can be decision-heavy without quickly
  // generating conflicts. The generic limit wrapper must expose a decision-cap
  // hit as UNKNOWN so callers skip the residual output instead of waiting.
  EXPECT_EQ(
      solver.solveWithResourceLimits(
          std::numeric_limits<unsigned>::max(),
          /*decisionLimit=*/0),
      SATSolverWrapper::SolveStatus::Unknown);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalCumulativeBudgetLeavesExactQueriesUnboundedAndAssumptionsLocal) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecPdrQuery();
  const int x = solver.newVar() + 2;
  solver.addClause({x});

  SATSolverWrapper::CadicalWorkBudget budget(
      /*conflictLimit=*/0,
      /*decisionLimit=*/0,
      /*tickLimit=*/0);
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
    // Exact Init queries do not opt into the output budget and must retain a
    // definite result even after property-local SAT work is exhausted.
    EXPECT_EQ(
        solver.solveWithAssumptionsStatus({-x}),
        SATSolverWrapper::SolveStatus::Unsat);
    EXPECT_EQ(budget.conflictsUsed(), 0u);
    EXPECT_EQ(budget.decisionsUsed(), 0u);
    EXPECT_EQ(budget.ticksUsed(), 0u);
    // The bounded query returns before enqueueing its assumption.
    EXPECT_EQ(
        solver.solveWithAssumptionsStatus(
            {-x},
            /*conflictLimit=*/100,
            /*propagationLimit=*/100),
        SATSolverWrapper::SolveStatus::Unknown);
  }

  EXPECT_TRUE(budget.exhausted());
  SATSolverWrapper::CadicalWorkBudget nextOutputBudget(
      /*conflictLimit=*/100,
      /*decisionLimit=*/100,
      /*tickLimit=*/100);
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(nextOutputBudget);
    EXPECT_EQ(
        solver.solveWithAssumptionsStatus(
            {x},
            /*conflictLimit=*/100,
            /*propagationLimit=*/100),
        SATSolverWrapper::SolveStatus::Sat);
  }
  EXPECT_FALSE(nextOutputBudget.exhausted());
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalCumulativeBudgetChargesActualWorkAcrossSolverInstances) {
  SATSolverWrapper firstSolver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  firstSolver.configureForSecPdrQuery();
  const int x = firstSolver.newVar() + 2;
  const int y = firstSolver.newVar() + 2;
  firstSolver.addClause({x, y});
  firstSolver.addClause({-x, y});

  SATSolverWrapper secondSolver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  secondSolver.configureForSecPdrQuery();
  const int secondX = secondSolver.newVar() + 2;
  secondSolver.addClause({secondX});

  SATSolverWrapper::CadicalWorkBudget budget(
      /*conflictLimit=*/100,
      /*decisionLimit=*/1,
      /*tickLimit=*/100);
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
    EXPECT_EQ(
        firstSolver.solveWithResourceLimits(
            /*conflictLimit=*/100,
            /*decisionLimit=*/100),
        SATSolverWrapper::SolveStatus::Unknown);
    EXPECT_EQ(budget.decisionsUsed(), 1u);
    EXPECT_TRUE(budget.exhausted());
    // A different incremental solver still belongs to the same output.
    EXPECT_EQ(
        secondSolver.solveWithResourceLimits(
            /*conflictLimit=*/100,
            /*decisionLimit=*/100),
        SATSolverWrapper::SolveStatus::Unknown);
  }

  EXPECT_EQ(secondSolver.solveStatus(), SATSolverWrapper::SolveStatus::Sat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalCumulativeBudgetBoundsPropagationHeavyWorkWithTicks) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecPdrQuery();
  const int x = solver.newVar() + 2;
  const int y = solver.newVar() + 2;
  solver.addClause({x, y});
  solver.addClause({-x, y});

  SATSolverWrapper::CadicalWorkBudget budget(
      /*conflictLimit=*/100,
      /*decisionLimit=*/100,
      /*tickLimit=*/1);
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
    EXPECT_EQ(
        solver.solveWithResourceLimits(
            /*conflictLimit=*/100,
            /*decisionLimit=*/100),
        SATSolverWrapper::SolveStatus::Unknown);
  }

  EXPECT_EQ(budget.ticksUsed(), 1u);
  EXPECT_TRUE(budget.exhausted());
  EXPECT_EQ(solver.solveStatus(), SATSolverWrapper::SolveStatus::Sat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalAssumptionTickLimitReturnsUnknownWithoutLeakingAssumptions) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecPdrQuery();
  const int x = solver.newVar() + 2;
  const int y = solver.newVar() + 2;
  solver.addClause({x, y});
  solver.addClause({-x, y});

  EXPECT_EQ(
      solver.solveWithAssumptionsStatus(
          {-y},
          /*conflictLimit=*/100,
          /*propagationLimit=*/100,
          /*tickLimit=*/0),
      SATSolverWrapper::SolveStatus::Unknown);
  // UNKNOWN clears the bounded query's assumptions and resets all limits on
  // the following exact call.
  EXPECT_EQ(
      solver.solveWithAssumptionsStatus({-y}),
      SATSolverWrapper::SolveStatus::Unsat);
  EXPECT_EQ(
      solver.solveWithAssumptionsStatus({y}),
      SATSolverWrapper::SolveStatus::Sat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalAssumptionQueriesDoNotPersistAssumptions) {
  EXPECT_EQ(
      SATSolverWrapper::assumptionSolverTypeFor(
          KEPLER_FORMAL::Config::SolverType::KISSAT),
      KEPLER_FORMAL::Config::SolverType::CADICAL);

  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.configureForSecPdrQuery();
  const int x = solver.newVar() + 2;
  solver.addClause({x});

  // PDR reuses one CaDiCaL context while changing only the target cube. A
  // failed target assumption must not become a permanent unit clause.
  EXPECT_EQ(
      solver.solveWithAssumptionsStatus({-x}),
      SATSolverWrapper::SolveStatus::Unsat);
  EXPECT_EQ(solver.failedAssumptions(), std::vector<int>{-x});
  EXPECT_EQ(
      solver.solveWithAssumptionsStatus({}),
      SATSolverWrapper::SolveStatus::Sat);
  EXPECT_EQ(
      solver.solveWithAssumptionsStatus({x}),
      SATSolverWrapper::SolveStatus::Sat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CadicalCraigInterpolationReturnsProofDerivedGlobalClause) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  solver.enableCraigInterpolation();
  solver.setCraigVariablePartition(
      SATSolverWrapper::CraigVariablePartition::Global);
  const int shared = solver.newVar() + 2;

  solver.setCraigClausePartition(
      SATSolverWrapper::CraigClausePartition::A);
  solver.addClause({-shared});
  solver.setCraigClausePartition(
      SATSolverWrapper::CraigClausePartition::B);
  solver.addClause({shared});

  ASSERT_EQ(
      solver.solveStatus(), SATSolverWrapper::SolveStatus::Unsat);
  const auto interpolant = solver.createCraigInterpolant();
  EXPECT_EQ(
      interpolant.type,
      SATSolverWrapper::CraigInterpolantCnf::Type::Normal);
  EXPECT_EQ(interpolant.clauses, std::vector<std::vector<int>>{{-shared}});
}

TEST_F(SequentialEquivalenceStrategyTests,
       KissatLargeSecConeProofProfileRemainsUsable) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  solver.configureForSecConeProof(/*coneSymbols=*/32768);
  const int x = solver.newVar() + 2;
  solver.addClause({x});
  solver.addClause({-x});

  // The medium-large SEC profile disables expensive speculative preprocessing.
  // Keep a direct regression guard that the selected Kissat options remain
  // compatible with the embedded solver used by CMake/Bazel.
  EXPECT_EQ(solver.solveStatus(), SATSolverWrapper::SolveStatus::Unsat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KissatDualRailMediumConeProofProfileRemainsUsable) {
  SATSolverWrapper solver(KEPLER_FORMAL::Config::SolverType::KISSAT);
  solver.configureForSecDualRailConeProof(/*coneSymbols=*/4096);
  const int x = solver.newVar() + 2;
  solver.addClause({x});
  solver.addClause({-x});

  // Dynamic-node dual-rail KI reaches medium-sized rail cones that are too
  // small for the generic large-cone cutoff but still suffer from speculative
  // Kissat preprocessing.  The dual-rail profile must remain a valid UNSAT
  // proof path.
  EXPECT_EQ(solver.solveStatus(), SATSolverWrapper::SolveStatus::Unsat);

  SATSolverWrapper cadicalSolver(KEPLER_FORMAL::Config::SolverType::CADICAL);
  cadicalSolver.configureForSecDualRailConeProof(/*coneSymbols=*/4096);
  const int y = cadicalSolver.newVar() + 2;
  cadicalSolver.addClause({y});
  cadicalSolver.addClause({-y});
  EXPECT_EQ(cadicalSolver.solveStatus(), SATSolverWrapper::SolveStatus::Unsat);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseValidationKeepsRequestedSolverAndUsesFastLocalProfile) {
  KInductionProblem problem;

  // KI base-case queries are one-shot validation probes, not incremental
  // assumption queries. Keep the selected solver for these BMC checks, but use
  // the fast local profile for embedded CDCL solvers so ASIC-sized SEC cases do
  // not spend minutes in standalone preprocessing.
  EXPECT_EQ(baseCaseValidationSolverType(
                problem, KEPLER_FORMAL::Config::SolverType::KISSAT),
            KEPLER_FORMAL::Config::SolverType::KISSAT);
  EXPECT_EQ(baseCaseValidationSolverType(
                problem, KEPLER_FORMAL::Config::SolverType::CADICAL),
            KEPLER_FORMAL::Config::SolverType::CADICAL);
  EXPECT_EQ(baseCaseValidationSolverType(
                problem, KEPLER_FORMAL::Config::SolverType::GLUCOSE),
            KEPLER_FORMAL::Config::SolverType::GLUCOSE);
  EXPECT_TRUE(baseCaseValidationUsesLocalQueryProfile(
      KEPLER_FORMAL::Config::SolverType::KISSAT));
  EXPECT_TRUE(baseCaseValidationUsesLocalQueryProfile(
      KEPLER_FORMAL::Config::SolverType::CADICAL));
  EXPECT_FALSE(baseCaseValidationUsesLocalQueryProfile(
      KEPLER_FORMAL::Config::SolverType::GLUCOSE));
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
       LocalBaseCaseCachePreservesBatchedNewestFrontierWitness) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"stable", "out"};
  problem.usesDualRailStateEncoding = true;
  problem.inputSymbols = {2};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::createFalse(), BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse(), BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(2)}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.property = BoolExpr::And(
      makeEqualityExpr(
          problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]),
      makeEqualityExpr(
          problem.observedOutputExprs0[1], problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const auto cache = makeImcBaseCounterexampleCache(problem);
  const auto kiCache = makeKInductionBaseCounterexampleCache(problem);

  // KI and IMC reuse this cache while sweeping depths and while localizing a
  // dual-rail residual output batch. Depth 0 is still before the
  // observation-only bad frontier, while depth 1 must match the public base
  // validator's witness.
  EXPECT_FALSE(findImcBaseCounterexampleAtFrontier(
      *cache, KEPLER_FORMAL::Config::SolverType::KISSAT, 0));
  EXPECT_FALSE(findBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0));
  const auto cachedWitness = findImcBaseCounterexampleAtFrontier(
      *cache, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);
  const auto kiCachedWitness = findKInductionBaseCounterexampleAtFrontier(
      *kiCache, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);
  const auto uncachedWitness = findBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  ASSERT_TRUE(cachedWitness.has_value());
  ASSERT_TRUE(kiCachedWitness.has_value());
  ASSERT_TRUE(uncachedWitness.has_value());
  EXPECT_EQ(cachedWitness->badFrame, uncachedWitness->badFrame);
  EXPECT_EQ(kiCachedWitness->badFrame, uncachedWitness->badFrame);
  ASSERT_EQ(cachedWitness->outputMismatches.size(), 1u);
  ASSERT_EQ(kiCachedWitness->outputMismatches.size(), 1u);
  ASSERT_EQ(uncachedWitness->outputMismatches.size(), 1u);
  EXPECT_EQ(cachedWitness->outputMismatches[0].signal, "out");
  EXPECT_EQ(kiCachedWitness->outputMismatches[0].signal, "out");
  EXPECT_EQ(cachedWitness->outputMismatches[0].signal,
            uncachedWitness->outputMismatches[0].signal);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PublicBaseCaseFrontierDoesNotRequireEarlierSafeFrames) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  problem.initialCondition = BoolExpr::Var(state);
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{state, BoolExpr::createTrue()}};
  problem.observedOutputNames = {"out"};
  problem.observedOutputExprs0 = {BoolExpr::Var(state)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  // Frontier BMC asks whether bad is reachable at this exact frame.  Requiring
  // earlier frames to be safe would incorrectly reject this valid frame-1
  // counterexample and would also widen dual-rail residual sweeps.
  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
  ASSERT_EQ(witness->outputMismatches.size(), 1u);
  EXPECT_EQ(witness->outputMismatches[0].signal, "out");
}

TEST_F(SequentialEquivalenceStrategyTests,
       LocalBaseCaseCacheProvesSafeMultiOutputFrontier) {
  KInductionProblem problem;
  problem.environmentInputNames = {"in"};
  problem.observedOutputNames = {"stable", "tracked"};
  problem.usesDualRailStateEncoding = true;
  problem.inputSymbols = {2};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.observedOutputExprs0 = {BoolExpr::createFalse(), BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse(), BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::Var(2)}};
  problem.transitions1 = {{4, BoolExpr::Var(2)}};
  problem.property = BoolExpr::And(
      makeEqualityExpr(
          problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]),
      makeEqualityExpr(
          problem.observedOutputExprs0[1], problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const auto cache = makeImcBaseCounterexampleCache(problem);
  const auto kiCache = makeKInductionBaseCounterexampleCache(problem);

  // Multi-output residual batches should prove the whole newest frontier safe
  // before splitting into per-output witness localization. This is the exact
  // fast path dual-rail KI residual runs need, and it must still agree with the
  // public base validator.
  EXPECT_FALSE(findImcBaseCounterexampleAtFrontier(
      *cache, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
  EXPECT_FALSE(findKInductionBaseCounterexampleAtFrontier(
      *kiCache, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
  EXPECT_FALSE(findBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0));
  EXPECT_FALSE(findBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailPublicBaseCacheCoversDefaultHorizon) {
  auto problem = std::make_unique<KInductionProblem>();
  problem->usesDualRailStateEncoding = true;
  problem->observedOutputNames = {"out0", "out1"};
  problem->state0Symbols = {2};
  problem->state1Symbols = {3};
  problem->allSymbols = {2, 3};
  problem->observedOutputExprs0 = {BoolExpr::Var(2), BoolExpr::Var(2)};
  problem->observedOutputExprs1 = {BoolExpr::Var(3), BoolExpr::Var(3)};
  problem->transitions0 = {{2, BoolExpr::Var(2)}};
  problem->transitions1 = {{3, BoolExpr::Var(3)}};
  for (size_t index = 0; index < 300; ++index) {
    problem->dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + index * 2, 1001 + index * 2});
  }
  problem->property = BoolExpr::And(
      makeEqualityExpr(
          problem->observedOutputExprs0[0], problem->observedOutputExprs1[0]),
      makeEqualityExpr(
          problem->observedOutputExprs0[1], problem->observedOutputExprs1[1]));
  problem->bad = BoolExpr::Not(problem->property);

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(findBaseCounterexampleAtFrontier(
      *problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0));
  EXPECT_FALSE(findBaseCounterexampleAtFrontier(
      *problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 31));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Large dual-rail residual sweeps should prove the normal default SEC
  // horizon once, then answer later safe frontiers from the exact prefix cache
  // instead of rebuilding the same AES-sized public base proof at 8,17,26,35.
  EXPECT_NE(
      stderrOutput.find("k-induction base coi k=32"),
      std::string::npos);
  EXPECT_EQ(
      detail::countTextOccurrences(stderrOutput, "k-induction base coi k=32"),
      1u);
  EXPECT_EQ(
      stderrOutput.find("k-induction base coi k=8"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LocalBaseCaseCacheChecksExactFrontierWithoutSafePrefixAssumption) {
  KInductionProblem problem;
  constexpr size_t stickyBadState = 2;
  problem.state0Symbols = {stickyBadState};
  problem.allSymbols = {stickyBadState};
  problem.initialCondition = BoolExpr::Var(stickyBadState);
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{stickyBadState, BoolExpr::createTrue()}};
  problem.bad = BoolExpr::Var(stickyBadState);
  problem.property = BoolExpr::Not(problem.bad);

  const auto cache = makeImcBaseCounterexampleCache(problem);
  const auto witness = findImcBaseCounterexampleAtFrontier(
      *cache, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  // Newest-frontier validation must not assume earlier frames are safe.  A bad
  // state at frame 1 is still a real counterexample even when frame 0 was also
  // bad; the caller's monotonic sweep decides which witness to report first.
  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverUsesFallbackWitnessNamesForUnnamedSignals) {
  KInductionProblem problem;
  constexpr size_t input = 2;
  problem.inputSymbols = {input};
  problem.allSymbols = {input};
  problem.observedOutputExprs0 = {BoolExpr::createFalse()};
  problem.observedOutputExprs1 = {BoolExpr::Var(input)};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0);

  ASSERT_TRUE(witness.has_value());
  ASSERT_EQ(witness->inputTrace.size(), 1u);
  ASSERT_EQ(witness->inputTrace[0].assignments.size(), 1u);
  EXPECT_EQ(witness->inputTrace[0].assignments[0].signal, "input_2");
  ASSERT_EQ(witness->outputMismatches.size(), 1u);
  EXPECT_EQ(witness->outputMismatches[0].signal, "output_0");
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverPdrProofOnlyRejectsReachableFrontierBad) {
  KInductionProblem problem;
  constexpr size_t badState = 2;
  problem.state0Symbols = {badState};
  problem.allSymbols = {badState};
  problem.initialStateAssignments = {{badState, true}};
  problem.initializedStateCount = 1;
  problem.initialCondition = BoolExpr::Var(badState);
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionBad = problem.bad;
  problem.inductionProperty = problem.property;
  problem.totalStateCount = problem.state0Symbols.size();

  // PDR callers use this helper as a proof certificate. A reachable bad
  // frontier is SAT and must never be collapsed into "no witness means UNSAT".
  EXPECT_FALSE(provesNoBaseCounterexampleAtFrontier(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0));
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
       BaseCaseSolverIncompleteResetBootstrapUsesObservationFrontier) {
  KInductionProblem problem;
  problem.environmentInputNames = {"rst"};
  problem.observedOutputNames = {"out"};
  problem.inputSymbols = {2};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{2, true}};
  problem.bootstrapStateAssignments = {{3, false}};
  problem.state0Symbols = {3};
  problem.state1Symbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.totalStateCount = 2;
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::Var(4)};
  problem.transitions0 = {{3, BoolExpr::createFalse()}};
  problem.transitions1 = {{4, BoolExpr::createFalse()}};
  problem.property = makeEqualityExpr(
      problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);

  // The reset summary leaves state1 arbitrary, so frame 0 is the top-level
  // observation frontier.  A mismatching internal value there is not a concrete
  // SEC counterexample unless it survives into a later visible cycle.
  EXPECT_FALSE(findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0).has_value());
  EXPECT_FALSE(findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1).has_value());

  KInductionProblem completeMismatch = problem;
  completeMismatch.bootstrapStateAssignments = {{3, false}, {4, true}};
  completeMismatch.transitions1 = {{4, BoolExpr::createTrue()}};
  const auto completeWitness = findBaseCounterexample(
      completeMismatch, KEPLER_FORMAL::Config::SolverType::KISSAT, 0);
  ASSERT_TRUE(completeWitness.has_value());
  EXPECT_EQ(completeWitness->badFrame, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResetBootstrapDoesNotUseBinaryObservationFrontier) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{2, true}};
  problem.bootstrapStateAssignments = {{3, true}, {4, true}};
  problem.state0Symbols = {3, 4};
  problem.state1Symbols = {5, 6};
  problem.totalStateCount = 4;
  problem.property = BoolExpr::Var(7);

  // Dual rail represents resetless startup values as unknown rails.  Reusing
  // the binary observation frontier would make PDR validate an over-approximate
  // startup mismatch instead of proving the selected dual-rail property.
  EXPECT_FALSE(problem.usesResetBootstrapObservationFrontier());

  KInductionProblem binaryProblem = problem;
  binaryProblem.usesDualRailStateEncoding = false;
  EXPECT_TRUE(binaryProblem.usesResetBootstrapObservationFrontier());
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverExactDualRailBootstrapDerivesConcreteRailValue) {
  KInductionProblem problem;
  constexpr size_t mayBeOne = 2;
  constexpr size_t mayBeZero = 3;
  constexpr size_t reset = 4;
  problem.usesDualRailStateEncoding = true;
  problem.inputSymbols = {reset};
  problem.state0Symbols = {mayBeOne, mayBeZero};
  problem.allSymbols = {mayBeOne, mayBeZero, reset};
  problem.totalStateCount = 2;
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.initialStateAssignments = {
      {mayBeOne, true}, {mayBeZero, true}};
  problem.initializedStateCount = 2;
  problem.bootstrapStateAssignments = {{mayBeOne, false}, {mayBeZero, true}};
  problem.dualRailStatePairs = {DualRailSymbolPair{mayBeOne, mayBeZero}};
  problem.transitions0 = {
      {mayBeOne, BoolExpr::createFalse()},
      {mayBeZero, BoolExpr::createTrue()}};
  problem.property = BoolExpr::Not(
      BoolExpr::And(BoolExpr::Var(mayBeOne), BoolExpr::Var(mayBeZero)));
  problem.bad = BoolExpr::Not(problem.property);

  EXPECT_FALSE(findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0).has_value());

  KInductionProblem initialProblem = problem;
  initialProblem.inputSymbols.clear();
  initialProblem.allSymbols = {mayBeOne, mayBeZero};
  initialProblem.resetBootstrapCycles = 0;
  initialProblem.resetBootstrapInputs.clear();
  initialProblem.bootstrapStateAssignments.clear();
  initialProblem.initialStateAssignments = {
      {mayBeOne, false}, {mayBeZero, true}};
  initialProblem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(mayBeOne)), BoolExpr::Var(mayBeZero));
  initialProblem.initializedStateCount = 2;

  EXPECT_FALSE(findBaseCounterexample(
      initialProblem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0).has_value());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialSteadyFrontierMismatchRequiresEngineValidation) {
  KInductionProblem combinationalProblem;
  EXPECT_TRUE(
      combinationalProblem.canReportSteadyFrontierMismatchAsCounterexample());

  KInductionProblem sequentialProblem;
  sequentialProblem.state0Symbols = {2};

  // A SAT steady-frontier result on sequential SEC can come from an arbitrary
  // reset-bootstrap startup assignment, so it is not a concrete counterexample
  // until the selected engine validates reachability.
  sequentialProblem.resetBootstrapCycles = 1;
  sequentialProblem.resetBootstrapInputs = {{3, true}};
  EXPECT_FALSE(
      sequentialProblem.canReportSteadyFrontierMismatchAsCounterexample());

  KInductionProblem completeBootstrapProblem = sequentialProblem;
  completeBootstrapProblem.bootstrapStateAssignments = {{2, true}};
  EXPECT_FALSE(
      completeBootstrapProblem.canReportSteadyFrontierMismatchAsCounterexample());

  KInductionProblem resetlessSequentialProblem;
  resetlessSequentialProblem.state0Symbols = {2};
  EXPECT_TRUE(
      resetlessSequentialProblem.canReportSteadyFrontierMismatchAsCounterexample());
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
       PDREngineFullyGeneralizesCheapConstantBlockedWideCubes) {
  KInductionProblem problem;
  // Keep this above PDR's per-group reserve-hint threshold. Missing a reserve
  // hint must not turn an otherwise exact, cheap proof into inconclusive.
  constexpr size_t kStateCount = 513;
  const size_t firstStateSymbol = 2;
  const size_t constantFalseSymbol = firstStateSymbol + kStateCount - 1;

  BoolExpr* bad = BoolExpr::createTrue();
  BoolExpr* init = BoolExpr::createTrue();
  problem.state0Symbols.reserve(kStateCount);
  problem.allSymbols.reserve(kStateCount);
  for (size_t index = 0; index < kStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    // Only the last target bit is impossible. The other bits share the opposite
    // SAT-root polarity, so this also checks that the failed-assumption core is
    // mapped back to the exact target literal rather than its first alias.
    problem.transitions0.emplace_back(
        symbol,
        symbol == constantFalseSymbol ? BoolExpr::createFalse()
                                      : BoolExpr::createTrue());
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kStateCount;
  problem.totalStateCount = kStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar secPdrTrace("KEPLER_SEC_PDR_TRACE", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("(!x" + std::to_string(constantFalseSymbol) + ")\n"),
      std::string::npos)
      << stderrOutput;

  // A broad scheduling probe may decline the same exact query before walking
  // more than 512 target transitions. UNKNOWN causes the caller to split; a
  // singleton/full run above still performs the exact proof.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  KInductionProblem probeProblem = problem;
  probeProblem.usesDualRailStateEncoding = true;
  testing::internal::CaptureStderr();
  PDREngine probeEngine(
      probeProblem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto probeResult = probeEngine.run(
      2, probeProblem.property,
      PDRQueryLimits{/*predecessorConflictLimit=*/10000,
                     /*predecessorDecisionLimit=*/150000,
                     /*blockingConflictLimit=*/10000,
                     /*blockingDecisionLimit=*/150000,
                     /*predecessorEncodingNodeLimit=*/5 * 1000 * 1000,
                     /*predecessorNodeHintTargetLimit=*/512});
  const std::string probeStderr = testing::internal::GetCapturedStderr();

  EXPECT_EQ(probeResult.status, PDRStatus::Inconclusive) << probeStderr;
  EXPECT_NE(
      probeStderr.find(
          "predecessor encoding budget exhausted targets=513 nodes=0 "
          "node_limit=5000000 node_hint_target_limit=512"),
      std::string::npos)
      << probeStderr;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesIncrementalCoreForExactWideBlockedCube) {
  KInductionProblem problem;
  constexpr size_t kStateCount = 96;
  constexpr size_t firstStateSymbol = 2;

  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();
  problem.state0Symbols.reserve(kStateCount);
  problem.allSymbols.reserve(kStateCount);
  for (size_t index = 0; index < kStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    problem.initialStateAssignments.push_back({symbol, false});
    // Keep the cheap 8-literal seed reachable, but make the full wide cube
    // unreachable through the remaining identity-held reset-low bits.
    problem.transitions0.emplace_back(
        symbol,
        index < 8 ? BoolExpr::createTrue() : BoolExpr::Var(symbol));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kStateCount;
  problem.totalStateCount = kStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // Section V returns the failed target assumptions from the same incremental
  // solveRelative query. Figure 7 should start from that core without a second
  // SAT query or a separate validation solver.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "generalized blocked cube level=1 size=96->1 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineKeepsFailedAssumptionCoreOutsideExactInit) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)),
      BoolExpr::Not(BoolExpr::Var(3)));
  problem.initialStateAssignments = {{2, false}, {3, false}};
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.bad = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)), BoolExpr::Var(3));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // The failed core is !x2, which overlaps Init. The authors' reduction keeps
  // an original target literal until the learned blocker is Init-disjoint.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "predecessor core kept outside init core=1->2 target=2"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesIncrementalCoreForWideCubeWithMediumSupport) {
  KInductionProblem problem;
  constexpr size_t kStateCount = 96;
  constexpr size_t kHeldStateCount = 16;
  constexpr size_t firstStateSymbol = 2;

  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();
  problem.state0Symbols.reserve(kStateCount);
  problem.allSymbols.reserve(kStateCount);
  for (size_t index = 0; index < kStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    problem.initialStateAssignments.push_back({symbol, false});
    const bool heldResetLow =
        index >= 8 && index < 8 + kHeldStateCount;
    problem.transitions0.emplace_back(
        symbol, heldResetLow ? BoolExpr::Var(symbol) : BoolExpr::createTrue());
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kStateCount;
  problem.totalStateCount = kStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // Failed assumptions must be used directly regardless of transition-cone
  // width; there is no support threshold in the paper's solveRelative flow.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "generalized blocked cube level=1 size=96->1 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesIncrementalCoreForMediumCubeWithWideSupport) {
  KInductionProblem problem;
  constexpr size_t kTargetStateCount = 12;
  constexpr size_t kSupportStateCount = 40;
  constexpr size_t firstStateSymbol = 2;
  constexpr size_t firstSupportSymbol = firstStateSymbol + kTargetStateCount;

  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();
  BoolExpr* wideSupport = BoolExpr::createTrue();
  problem.state0Symbols.reserve(kTargetStateCount + kSupportStateCount);
  problem.allSymbols.reserve(kTargetStateCount + kSupportStateCount);
  for (size_t index = 0; index < kSupportStateCount; ++index) {
    const size_t symbol = firstSupportSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    wideSupport = BoolExpr::And(wideSupport, BoolExpr::Var(symbol));
  }
  for (size_t index = 0; index < kTargetStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    problem.initialStateAssignments.push_back({symbol, false});
    // The first four target bits remain reachable. The remaining target bits
    // depend on a broad support cone that is false in the startup frontier,
    // matching the measured AES 12-literal, 113-support level-zero blockers.
    problem.transitions0.emplace_back(
        symbol,
        index < 4
            ? BoolExpr::createTrue()
            : BoolExpr::And(BoolExpr::Var(symbol), wideSupport));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kTargetStateCount + kSupportStateCount;
  problem.totalStateCount = kTargetStateCount + kSupportStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "generalized blocked cube level=1 size=12->1 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesIncrementalCoreForDualRailBlockedCube) {
  KInductionProblem problem;
  constexpr size_t kTargetStateCount = 12;
  constexpr size_t kSupportStateCount = 40;
  constexpr size_t firstStateSymbol = 2;
  constexpr size_t firstSupportSymbol = firstStateSymbol + kTargetStateCount;

  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();
  BoolExpr* wideSupport = BoolExpr::createTrue();
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols.reserve(kTargetStateCount + kSupportStateCount);
  problem.allSymbols.reserve(kTargetStateCount + kSupportStateCount);
  for (size_t index = 0; index < kSupportStateCount; ++index) {
    const size_t symbol = firstSupportSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    wideSupport = BoolExpr::And(wideSupport, BoolExpr::Var(symbol));
  }
  for (size_t index = 0; index < kTargetStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(
        symbol,
        index < 4
            ? BoolExpr::createTrue()
            : BoolExpr::And(BoolExpr::Var(symbol), wideSupport));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kTargetStateCount + kSupportStateCount;
  problem.totalStateCount = kTargetStateCount + kSupportStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "generalized blocked cube level=1 size=12->1 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesIncrementalCoreForBroadDualRailBlockedCube) {
  KInductionProblem problem;
  constexpr size_t kTargetStateCount = 12;
  constexpr size_t kSupportStateCount = 160;
  constexpr size_t firstStateSymbol = 2;
  constexpr size_t firstSupportSymbol = firstStateSymbol + kTargetStateCount;

  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();
  BoolExpr* wideSupport = BoolExpr::createTrue();
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols.reserve(kTargetStateCount + kSupportStateCount);
  problem.allSymbols.reserve(kTargetStateCount + kSupportStateCount);
  for (size_t index = 0; index < kSupportStateCount; ++index) {
    const size_t symbol = firstSupportSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    wideSupport = BoolExpr::And(wideSupport, BoolExpr::Var(symbol));
  }
  for (size_t index = 0; index < kTargetStateCount; ++index) {
    const size_t symbol = firstStateSymbol + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(
        symbol,
        index < 4
            ? BoolExpr::createTrue()
            : BoolExpr::And(BoolExpr::Var(symbol), wideSupport));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = kTargetStateCount + kSupportStateCount;
  problem.totalStateCount = kTargetStateCount + kSupportStateCount;
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "generalized blocked cube level=1 size=12->1 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineIndexedInitFactsPreserveWideExactFrameZero) {
  constexpr size_t kStateCount = 2048;
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols.reserve(kStateCount);
  problem.allSymbols.reserve(kStateCount);
  problem.initialStateAssignments.reserve(kStateCount);
  problem.transitions0.reserve(kStateCount);
  for (size_t index = 0; index < kStateCount; ++index) {
    const size_t symbol = index + 2;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.emplace_back(symbol, false);
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = kStateCount;
  problem.totalStateCount = kStateCount;
  problem.bad = BoolExpr::Var(problem.state0Symbols.back());
  problem.property = BoolExpr::Not(problem.bad);

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // Exact F[0] fixes every bit low. The immutable assignment index changes
  // only lookup cost; the final bad bit must remain unreachable.
  EXPECT_EQ(result.status, PDRStatus::Equivalent);
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
       PDREngineFindsOneStepCounterexampleForDocumentedBooleanMiter) {
  const auto problem = buildDocumentedBooleanPdrCounterexampleProblem();

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationRemovesIrrelevantBadStateLiterals) {
  auto problem = buildDocumentedBooleanPdrCounterexampleProblem();
  problem.state0Symbols.push_back(5);
  problem.allSymbols.push_back(5);
  problem.initialCondition = BoolExpr::And(
      problem.initialCondition, BoolExpr::Not(BoolExpr::Var(5)));
  problem.transitions0.emplace_back(
      5, BoolExpr::Xor(BoolExpr::Var(5), BoolExpr::Var(4)));

  const ScopedEnvVar secPdrTrace("KEPLER_SEC_PDR_TRACE", "1");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  ASSERT_EQ(result.status, PDRStatus::Different);
  const auto badCubePos = stderrOutput.find("SEC PDR trace: bad_cube@F1");
  ASSERT_NE(badCubePos, std::string::npos);
  const auto nextTracePos = stderrOutput.find("SEC PDR trace:", badCubePos + 1);
  const std::string badCubeTrace = stderrOutput.substr(
      badCubePos,
      nextTracePos == std::string::npos ? std::string::npos
                                        : nextTracePos - badCubePos);
  EXPECT_NE(badCubeTrace.find("x2="), std::string::npos);
  EXPECT_NE(badCubeTrace.find("x3="), std::string::npos);
  // Section III-B of the FMCAD'11 PDR paper probes each state bit with X.
  // x5 cannot affect the bad XOR, so it must not survive in the obligation.
  EXPECT_EQ(badCubeTrace.find("x5="), std::string::npos);
  // Section V's on-demand encoding also avoids materializing x5 before the
  // same ternary reduction; the exact bad predicate has two state inputs.
  EXPECT_NE(
      stderrOutput.find("state_symbols=2 model_cube=2"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationShrinksPredecessorCubes) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 5};
  problem.allSymbols = problem.state0Symbols;
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)),
      BoolExpr::And(BoolExpr::Var(3), BoolExpr::Not(BoolExpr::Var(5))));
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0 = {
      {2, BoolExpr::Var(3)},
      {3, BoolExpr::Var(3)},
      {5, BoolExpr::Var(5)},
  };
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  ASSERT_EQ(result.status, PDRStatus::Different);
  // Only x3 controls x2'. The current values of x2 and x5 must be removed by
  // the paper's ternary predecessor simulation before the cube is queued.
  EXPECT_NE(stderrOutput.find("predecessor_cube=1"), std::string::npos)
      << stderrOutput;
  // Exact F[0] still contributes all three initialized state constraints, but
  // only x3 is requested as a predecessor model value.
  EXPECT_NE(
      stderrOutput.find("predecessor_symbols=1 solver_symbols=3"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationStopsPropagationAtStableParent) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = problem.state0Symbols;
  problem.initialStateAssignments = {
      {2, false}, {3, true}, {4, true}};
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0 = {
      {2, BoolExpr::Or(BoolExpr::Var(3), BoolExpr::Var(4))},
      {3, BoolExpr::Var(3)},
      {4, BoolExpr::Var(4)}};
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  // Either true input controls the OR. Replacing the first input by X leaves
  // the exact ternary value stable, so incremental cache propagation stops at
  // that parent without changing the paper's literal-removal result.
  EXPECT_NE(
      stderrOutput.find("ternary incremental propagation changed_values="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("stable_parents="), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationKeepsRemovedLiteralsUnknown) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = problem.state0Symbols;
  problem.initialStateAssignments = {{2, false}, {3, false}};
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0 = {
      {2, BoolExpr::createTrue()}, {3, BoolExpr::createTrue()}};
  problem.bad = BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  ASSERT_EQ(result.status, PDRStatus::Different);
  ASSERT_EQ(result.bound, 1u);
  // With x2=x3=1, either literal alone initially appears removable from OR.
  // The paper keeps the first removal at X while testing the second, leaving
  // one literal in the generalized bad cube instead of unsoundly removing both.
  EXPECT_NE(
      stderrOutput.find(
          "bad cube level=1 source=ternary_simulation state_symbols=2 "
          "model_cube=2 cube=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationPropagatesRestoredLiteralDependencies) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = problem.state0Symbols;
  problem.initialStateAssignments = {{2, true}, {3, false}, {4, false}};
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0 = {
      {2, BoolExpr::Var(2)},
      {3, BoolExpr::Var(3)},
      {4,
       BoolExpr::And(
           BoolExpr::Var(2),
           BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3)))}};
  problem.bad = BoolExpr::Var(4);
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  ASSERT_EQ(result.status, PDRStatus::Different);
  ASSERT_EQ(result.bound, 1u);
  // Probing x2 with X fails, so x2 is restored before x3 is probed. With
  // x2=1, x3 is irrelevant and the exact paper reduction keeps only x2.
  // Reusing the failed probe's stale X value would incorrectly retain x3.
  EXPECT_NE(stderrOutput.find("predecessor_cube=1"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("ternary incremental propagation changed_values="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("ternary evaluation memo storage reused entries="),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineTernarySimulationPreservesSharedTransitionRoots) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = problem.state0Symbols;
  problem.initialStateAssignments = {{2, false}, {3, false}, {4, true}};
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;

  // Both target bits use the same transition DAG. The ternary reducer may
  // share evaluation work, but it must still retain the controlling x4 value.
  BoolExpr* sharedTransition = BoolExpr::Var(4);
  problem.transitions0 = {
      {2, sharedTransition}, {3, sharedTransition}, {4, BoolExpr::Var(4)}};
  problem.bad = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  auto exactCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine engine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0, exactCache);
  const auto result = engine.run(1);
  const auto repeatedResult = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  // Exact metadata and dense ternary memo allocations survive serial output
  // runs, but every reducer must still produce the same predecessor.
  EXPECT_EQ(repeatedResult.status, result.status);
  EXPECT_EQ(repeatedResult.bound, result.bound);
  EXPECT_NE(stderrOutput.find("predecessor_cube=1"), std::string::npos)
      << stderrOutput;
  // Each tentative X assignment propagates through computed ancestors while
  // unaffected values reuse the shared transition-DAG evaluation memo.
  EXPECT_NE(
      stderrOutput.find("ternary evaluation memo storage reused entries="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("ternary mapped support cache reused="),
      std::string::npos)
      << stderrOutput;
  // The support index must retain the one controlling symbol shared by both
  // roots; this avoids rescanning unrelated roots without changing reduction.
  EXPECT_NE(stderrOutput.find("root_index_symbols=1"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBlockingUsesPaperRelativeInductionQuery) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{2, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  ASSERT_EQ(result.status, PDRStatus::Equivalent);
  // Query Q2 in Figure 1 includes !s in the current frame.
  EXPECT_NE(stderrOutput.find("exclude_target=1"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesBootstrapAssignmentsAndComplementedStatePairs) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {{2, false}, {3, true}};
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
       PDREngineDualRailUsesSameFrameComplementRailEqualities) {
  KInductionProblem problem = makeDualRailComplementedOutputProblemForTest();

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRelationComponentCachePreservesTransitiveClosure) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = problem.state0Symbols;
  problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}};
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0 = {
      {2, BoolExpr::Var(2)},
      {3, BoolExpr::Var(3)},
      {4, BoolExpr::Var(4)},
  };
  problem.sameFrameStateEqualityPairs0 = {{2, 3}, {3, 4}};
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // The target cone initially contains only x2. Its immutable equality
  // component must still close transitively through x3 to x4 before the
  // target surface is retained.
  EXPECT_NE(stderrOutput.find("predecessor target surface cached target=1"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("cached solver created level=0 symbols=3"),
            std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineIndexedRelationsPreserveEverySameDesignConstraint) {
  constexpr size_t complement0 = 2;
  constexpr size_t complement0N = 3;
  constexpr size_t equality0 = 4;
  constexpr size_t equality0Peer = 5;
  constexpr size_t railOne = 6;
  constexpr size_t railZero = 7;
  constexpr size_t complement1 = 8;
  constexpr size_t complement1N = 9;
  constexpr size_t equality1 = 10;
  constexpr size_t equality1Peer = 11;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {complement0,   complement0N, equality0,
                           equality0Peer, railOne,      railZero};
  problem.state1Symbols = {complement1, complement1N, equality1, equality1Peer};
  problem.allSymbols = problem.state0Symbols;
  problem.allSymbols.insert(problem.allSymbols.end(),
                            problem.state1Symbols.begin(),
                            problem.state1Symbols.end());
  problem.complementedStatePairs0 = {{complement0, complement0N}};
  problem.complementedStatePairs1 = {{complement1, complement1N}};
  problem.sameFrameStateEqualityPairs0 = {{equality0, equality0Peer}};
  problem.sameFrameStateEqualityPairs1 = {{equality1, equality1Peer}};
  problem.dualRailStatePairs = {{railOne, railZero}};
  for (size_t index = 0; index < 32; ++index) {
    const size_t base = 1000 + index * 16;
    problem.complementedStatePairs0.emplace_back(base, base + 1);
    problem.complementedStatePairs1.emplace_back(base + 2, base + 3);
    problem.sameFrameStateEqualityPairs0.emplace_back(base + 4, base + 5);
    problem.sameFrameStateEqualityPairs1.emplace_back(base + 6, base + 7);
    problem.dualRailStatePairs.push_back({base + 8, base + 9});
  }
  problem.totalStateCount = problem.allSymbols.size();

  for (const size_t symbol : problem.state0Symbols) {
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  for (const size_t symbol : problem.state1Symbols) {
    problem.transitions1.emplace_back(symbol, BoolExpr::Var(symbol));
  }

  // Every disjunct violates one relation kind. The unrelated pairs force the
  // sparse index path; if it drops a relevant pair, F[0] admits a false CEX.
  BoolExpr* bad =
      makeEqualityExpr(BoolExpr::Var(complement0), BoolExpr::Var(complement0N));
  bad = BoolExpr::Or(bad, makeEqualityExpr(BoolExpr::Var(complement1),
                                           BoolExpr::Var(complement1N)));
  bad = BoolExpr::Or(bad, BoolExpr::Xor(BoolExpr::Var(equality0),
                                        BoolExpr::Var(equality0Peer)));
  bad = BoolExpr::Or(bad, BoolExpr::Xor(BoolExpr::Var(equality1),
                                        BoolExpr::Var(equality1Peer)));
  bad = BoolExpr::Or(bad, BoolExpr::Not(BoolExpr::Or(BoolExpr::Var(railOne),
                                                     BoolExpr::Var(railZero))));
  problem.bad = bad;
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionBad = problem.bad;
  problem.inductionProperty = problem.property;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDeterministicWorklistSortsHashSetSymbols) {
  std::unordered_set<size_t> symbols;
  std::vector<size_t> rawIteration;
  std::vector<size_t> expected;

  for (const size_t symbol : {
           2u, 3u, 5u, 7u, 11u, 17u, 31u, 64u, 127u, 257u, 1025u}) {
    symbols.insert(symbol);
    rawIteration.assign(symbols.begin(), symbols.end());
    expected = rawIteration;
    std::sort(expected.begin(), expected.end());
    if (rawIteration != expected) {
      break;
    }
  }
  ASSERT_NE(rawIteration, expected)
      << "test data must expose non-sorted unordered_set traversal";

  EXPECT_EQ(detail::makeDeterministicPdrWorklist(symbols), expected);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrProofObligationsUseLowestFrameThenLifoOrder) {
  // Figure 6 orders by the lowest frame and reports a small benefit from
  // stack-like behavior among obligations in the same frame.
  EXPECT_TRUE(detail::pdrProofObligationPriorityLess(
      /*lhsLevel=*/0, /*lhsSequence=*/1,
      /*rhsLevel=*/1, /*rhsSequence=*/100));
  EXPECT_TRUE(detail::pdrProofObligationPriorityLess(
      /*lhsLevel=*/2, /*lhsSequence=*/8,
      /*rhsLevel=*/2, /*rhsSequence=*/7));
  EXPECT_FALSE(detail::pdrProofObligationPriorityLess(
      /*lhsLevel=*/2, /*lhsSequence=*/7,
      /*rhsLevel=*/2, /*rhsSequence=*/8));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrMergeSortedSymbolVectorsKeepsStableSurfaceSortedUnique) {
  const std::vector<size_t> stableSurface = {2, 4, 8, 16};
  const std::vector<size_t> localSymbols = {3, 4, 16, 32};

  const std::vector<size_t> merged =
      detail::mergeSortedPdrSymbolVectors(stableSurface, localSymbols);

  // Incremental predecessor caching and bad-state queries use this instead of
  // rebuilding a large unordered_set on every query. The result must remain
  // the exact sorted union used by FrameVariableStore.
  const std::vector<size_t> expected = {2, 3, 4, 8, 16, 32};
  EXPECT_EQ(merged, expected);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrWidenSortedSymbolSurfaceOnlyAddsMissingSymbols) {
  std::vector<size_t> stableSurface = {2, 4, 8};

  EXPECT_FALSE(
      detail::widenSortedPdrSymbolSurface(stableSurface, {2, 8}));
  EXPECT_EQ(stableSurface, (std::vector<size_t>{2, 4, 8}));

  // Section V's incremental predecessor solver uses this to stay alive when
  // neighboring target cubes add support symbols. The widened surface must
  // stay sorted and unique because it becomes the SAT variable list.
  EXPECT_TRUE(
      detail::widenSortedPdrSymbolSurface(stableSurface, {3, 4, 9}));
  EXPECT_EQ(stableSurface, (std::vector<size_t>{2, 3, 4, 8, 9}));

  EXPECT_FALSE(
      detail::widenSortedPdrSymbolSurface(stableSurface, {2, 3, 9}));
  EXPECT_EQ(stableSurface, (std::vector<size_t>{2, 3, 4, 8, 9}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrFrameSymbolSurfaceCacheRetainsInterleavedFrames) {
  detail::PdrFrameSymbolSurfaceCache cache;
  const int firstModel = 0;
  const int secondModel = 0;

  EXPECT_EQ(cache.widen(&firstModel, 1, {2, 6}),
            (std::vector<size_t>{2, 6}));
  EXPECT_EQ(cache.widen(&firstModel, 2, {3, 7}),
            (std::vector<size_t>{3, 7}));
  // Returning to F[1] must retain its earlier surface instead of replacing it
  // with the most recently visited frame's cache entry.
  EXPECT_EQ(cache.widen(&firstModel, 1, {4}),
            (std::vector<size_t>{2, 4, 6}));
  EXPECT_EQ(cache.widen(&firstModel, 2, {3}),
            (std::vector<size_t>{3, 7}));

  // Symbol identities are local to one transition model.
  EXPECT_EQ(cache.widen(&secondModel, 1, {5}),
            (std::vector<size_t>{5}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrWeightedLruCacheEvictsOnlyLeastRecentEntries) {
  detail::PdrWeightedLruCache<int, std::string, std::hash<int>> cache(5);

  EXPECT_EQ(cache.insert(1, "one", 2).evictedEntries, 0u);
  EXPECT_EQ(cache.insert(2, "two", 3).evictedEntries, 0u);
  ASSERT_NE(cache.find(1), nullptr);

  // Touching key 1 makes key 2 least recent. Inserting key 3 should preserve
  // that hot exact value instead of clearing the complete cache at its bound.
  const auto inserted = cache.insert(3, "three", 3);
  ASSERT_NE(inserted.value, nullptr);
  EXPECT_EQ(inserted.evictedEntries, 1u);
  EXPECT_EQ(*inserted.value, "three");
  EXPECT_NE(cache.find(1), nullptr);
  EXPECT_EQ(cache.find(2), nullptr);
  EXPECT_NE(cache.find(3), nullptr);
  EXPECT_EQ(cache.size(), 2u);
  EXPECT_EQ(cache.retainedWeight(), 5u);

  // An entry larger than the complete byte budget remains uncached and must
  // not evict exact entries that are still reusable.
  const auto oversized = cache.insert(4, "four", 6);
  EXPECT_EQ(oversized.value, nullptr);
  EXPECT_EQ(oversized.evictedEntries, 0u);
  EXPECT_EQ(cache.size(), 2u);
  EXPECT_EQ(cache.retainedWeight(), 5u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrStableUnsatCacheUsesItsOwnEntryBound) {
  constexpr size_t stableUnsatLimit = 8;

  // Exact-result LRU churn is intentionally absent from this policy: stable
  // UNSAT facts remain valid as IC3 strengthens a frame.
  EXPECT_FALSE(detail::shouldResetPdrStableUnsatCache(
      stableUnsatLimit - 1, stableUnsatLimit));
  EXPECT_TRUE(detail::shouldResetPdrStableUnsatCache(
      stableUnsatLimit, stableUnsatLimit));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrPredecessorUnsatCoreSharingUsesBaseContextOnly) {
  // A predecessor UNSAT core can be reused for stronger target cubes only when
  // it came from the monotonic base frame context. Q2 selector assumptions stay
  // target-local.
  EXPECT_TRUE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/0,
      /*excludeTargetOnCurrentFrame=*/false));
  EXPECT_FALSE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/7,
      /*excludeTargetOnCurrentFrame=*/false));
  EXPECT_FALSE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/0,
      /*excludeTargetOnCurrentFrame=*/true));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDeterministicCubeOrderingSortsCandidates) {
  using CubeKey = std::vector<std::pair<size_t, bool>>;
  std::vector<CubeKey> cubes = {
      {{5, true}, {8, false}},
      {{3, true}},
      {{5, false}, {9, true}},
      {{5, false}, {8, true}},
      {{5, false}, {8, false}},
  };

  std::sort(
      cubes.begin(),
      cubes.end(),
      [](const CubeKey& lhs, const CubeKey& rhs) {
        return detail::pdrCubeAssignmentOrderLess(lhs, rhs);
      });

  const std::vector<CubeKey> expected = {
      {{3, true}},
      {{5, false}, {8, false}},
      {{5, false}, {8, true}},
      {{5, false}, {9, true}},
      {{5, true}, {8, false}},
  };
  EXPECT_EQ(cubes, expected);
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
       PDREngineTriesSharedStrengtheningBeforeWeakStateSubsetFallback) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4, 5, 6, 7};
  problem.allSymbols = problem.state0Symbols;
  BoolExpr* init = BoolExpr::createTrue();
  for (const auto symbol : problem.state0Symbols) {
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
  }
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createFalse());
  problem.transitions0.emplace_back(5, BoolExpr::createFalse());
  problem.transitions0.emplace_back(6, BoolExpr::createTrue());
  problem.transitions0.emplace_back(7, BoolExpr::createFalse());
  // PDR should prefer the validated shared strengthening below before falling
  // back to its regular frame-clause loop.
  problem.observedOutputExprs0 = {BoolExpr::Var(6)};
  problem.observedOutputExprs1 = {BoolExpr::Var(7)};
  problem.bad = BoolExpr::Var(5);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = BoolExpr::Not(BoolExpr::Var(4));
  problem.inductionBad = BoolExpr::Not(problem.inductionProperty);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("frame invariant shared_strengthening"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("PDR using validated SEC strengthening frame invariant"),
      std::string::npos);
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
       PDREngineProvesEquivalentWithinThreeFrames) {
  const auto problem = buildLinearChainSecProblem(4);

  // This is an engine-regression check for the current binary-chain model and
  // current clause-generalization behavior.  It is not a portable "classic PDR
  // must prove safe exactly at k=3" theorem: safe IC3/PDR proofs may converge
  // earlier whenever a stronger inductive invariant is learned.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
  // Figure 7 visits each original literal once. A four-literal cube must not
  // trigger a fifth Q2 probe by retrying an earlier SAT removal after the cube
  // changes, which is the stronger generalization rejected in Section VI-B.
  EXPECT_EQ(stderrOutput.find("size=4->2 checks=5"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("size=4->2 checks=4"), std::string::npos)
      << stderrOutput;
  // Section V's default solveRelative call returns only SAT/UNSAT. Figure 7
  // does not request EXTRACTMODEL, so generalization must not decode and
  // ternary-reduce a predecessor cube that its caller immediately discards.
  EXPECT_NE(
      stderrOutput.find(
          "result=sat cached_assumptions=1 model_extracted=0 "
          "purpose=generalize"),
      std::string::npos)
      << stderrOutput;
  // Once Figure 7 finishes, its temporary Q2 selectors no longer help
  // recursive blocking. Retire them so long exact runs do not accumulate
  // inactive generalization contexts in the shared SAT solver.
  EXPECT_NE(
      stderrOutput.find(
          "Q2 status selectors retired after generalization count="),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailGeneralizationUsesNarrowExactSolver) {
  const auto problem = buildSharedPdrNarrowProbeProblem();
  BoolExpr* narrowProperty = makeEqualityExpr(
      problem.observedOutputExprs0.front(),
      problem.observedOutputExprs1.front());
  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);
  const auto broadResult = engine.run(3, problem.property);
  const auto narrowResult = engine.run(3, narrowProperty);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(broadResult.status, PDRStatus::Different) << stderrOutput;
  EXPECT_EQ(narrowResult.status, PDRStatus::Equivalent) << stderrOutput;
  // Once the level-local persistent solver has widened, Figure 7 keeps this
  // cube's exact status-only queries in their own smaller SAT context.
  EXPECT_NE(
      stderrOutput.find("narrow_generalization_probe result="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("persistent_request_symbols="),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCumulativeCadicalBudgetExhaustionIsInconclusive) {
  constexpr size_t state = 2;
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  problem.totalStateCount = 1;
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(state));
  problem.initialStateAssignments = {{state, false}};
  problem.initializedStateCount = 1;
  problem.transitions0 = {{state, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(state);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(state)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"bounded_output"};

  SATSolverWrapper::CadicalWorkBudget exhaustedBudget(
      /*conflictLimit=*/0,
      /*decisionLimit=*/0,
      /*tickLimit=*/0);
  PDRResult boundedResult;
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(exhaustedBudget);
    PDREngine boundedEngine(
        problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    boundedResult = boundedEngine.run(2);
  }

  PDREngine unboundedEngine(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto unboundedResult = unboundedEngine.run(2);

  EXPECT_EQ(boundedResult.status, PDRStatus::Inconclusive);
  EXPECT_TRUE(exhaustedBudget.exhausted());
  EXPECT_EQ(unboundedResult.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineNarrowGeneralizationUnknownFallsBackToPersistentSolver) {
  const auto problem = buildSharedPdrNarrowProbeProblem();
  BoolExpr* narrowProperty = makeEqualityExpr(
      problem.observedOutputExprs0.front(),
      problem.observedOutputExprs1.front());
  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);
  const auto broadResult = engine.run(3, problem.property);
  const auto narrowResult = engine.run(
      3,
      narrowProperty,
      PDRQueryLimits{
          /*predecessorConflictLimit=*/0,
          /*predecessorDecisionLimit=*/1,
          /*blockingConflictLimit=*/250 * 1000,
          /*blockingDecisionLimit=*/10 * 1000 * 1000});
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(broadResult.status, PDRStatus::Different) << stderrOutput;
  EXPECT_EQ(narrowResult.status, PDRStatus::Inconclusive) << stderrOutput;
  // UNKNOWN from the bounded narrow attempt is not interpreted as blocked.
  // The same exact Q2 query must continue through the established solver path.
  EXPECT_NE(
      stderrOutput.find(
          "narrow_generalization_probe result=unknown fallback=persistent"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "predecessor query budget exhausted limit=0 decision_limit=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDoesNotSubstituteRetainedModelForExactPredecessorSolve) {
  const auto problem = buildLinearChainSecProblem(4);
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");

  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // Reusing the prepared CNF and assumptions is safe. Reusing the previous SAT
  // assignment as the answer changes finite-budget witness selection and once
  // collapsed strict dual-rail SEC coverage from 539/598 outputs to 0/598.
  EXPECT_NE(
      stderrOutput.find("predecessor target assumptions reused"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor cached SAT model reused"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineLiftsBlockedCubeThroughRelativeInductiveFrames) {
  const auto problem = buildLinearChainSecProblem(6);
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");

  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(5);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // Figure 6 repeatedly calls solveRelative(next(z)); this fixture has a
  // blocker first learned in F2 that is also relatively inductive in F3.
  EXPECT_NE(
      stderrOutput.find("blocked cube lifted level=2->3"),
      std::string::npos)
      << stderrOutput;
  // Reusing a predecessor solver at a stable symbol surface should consume
  // only clauses appended since its previous query.
  EXPECT_NE(
      stderrOutput.find(
          "predecessor cached solver frame sync source=frame_log"),
      std::string::npos)
      << stderrOutput;
  // Repeated solveRelative calls for one obligation must reuse only the
  // prepared assumption surface; each call still performs its exact SAT solve.
  EXPECT_NE(
      stderrOutput.find("predecessor target assumptions reused"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineProvesClassicSafeChainsWithinDepthsThreeFourFive) {
  for (const size_t proofDepth : {size_t{3}, size_t{4}, size_t{5}}) {
    const auto problem =
        buildClassicPdrOneHotUnreachableBadChainProblem(proofDepth);

    // Do not tighten this to bound == proofDepth.  For safe properties, exact
    // convergence depth is not a classic PDR/IC3 semantic contract: clause
    // generalization is allowed to learn the whole unreachable suffix earlier.
    // Exact depth is meaningful for the reachable counterexample test below.
    PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(proofDepth);

    EXPECT_EQ(result.status, PDRStatus::Equivalent)
        << "proofDepth=" << proofDepth;
    EXPECT_LE(result.bound, proofDepth) << "proofDepth=" << proofDepth;
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFindsClassicCounterexamplesAtDepthsThreeFourFive) {
  for (const size_t badDepth : {size_t{3}, size_t{4}, size_t{5}}) {
    const auto problem =
        buildClassicPdrOneHotReachableBadChainProblem(badDepth);

    PDREngine earlyEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto earlyResult = earlyEngine.run(badDepth - 1);

    // Figure 6 requeues next(s), specifically allowing a counterexample longer
    // than the current trace. A trace one frame short must still find this bug.
    EXPECT_EQ(earlyResult.status, PDRStatus::Different)
        << "badDepth=" << badDepth;
    EXPECT_EQ(earlyResult.bound, badDepth) << "badDepth=" << badDepth;

    PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(badDepth);

    EXPECT_EQ(result.status, PDRStatus::Different)
        << "badDepth=" << badDepth;
    EXPECT_EQ(result.bound, badDepth) << "badDepth=" << badDepth;
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDebugFormattingPrintsDocumentedBooleanMiterProblemAndFrames) {
  const auto problem = buildDocumentedBooleanPdrCounterexampleProblem();

  const std::string formattedProblem = formatKInductionProblemForDebug(problem);
  EXPECT_NE(formattedProblem.find("description: documented Boolean PDR counterexample miter"),
            std::string::npos);
  EXPECT_NE(formattedProblem.find("state0_symbols: [x2]"), std::string::npos);
  EXPECT_NE(formattedProblem.find("state1_symbols: [x3]"), std::string::npos);
  EXPECT_NE(formattedProblem.find("input_symbols: [x4]"), std::string::npos);
  EXPECT_NE(formattedProblem.find("transition_formula:"), std::string::npos);
  EXPECT_NE(formattedProblem.find("x2' = ~x2"), std::string::npos);
  EXPECT_TRUE(
      formattedProblem.find("x3' = x3 AND x4") != std::string::npos ||
      formattedProblem.find("x3' = x4 AND x3") != std::string::npos)
      << formattedProblem;
  EXPECT_NE(formattedProblem.find("bad: x2 XOR x3"), std::string::npos);

  const ScopedEnvVar secPdrTrace("KEPLER_SEC_PDR_TRACE", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_NE(stderrOutput.find("SEC PDR trace: problem"), std::string::npos);
  EXPECT_NE(stderrOutput.find("transition_formula:"), std::string::npos);
  EXPECT_NE(stderrOutput.find("SEC PDR trace: seeded_frames"), std::string::npos);
  EXPECT_NE(stderrOutput.find("F[1]"), std::string::npos);
  EXPECT_NE(stderrOutput.find("SEC PDR trace: bad_cube@F1"), std::string::npos);
  EXPECT_NE(stderrOutput.find("{x2=1, x3=0}"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ProofProblemDebugFormatsConstantsAndInvalidExpressions) {
  BoolExpr invalid;
  KInductionProblem problem;
  problem.description = "debug edge cases";
  problem.allSymbols = {0, 1, 2};
  problem.state0Symbols = {0, 1};
  problem.transitions0.emplace_back(0, &invalid);
  problem.property = BoolExpr::Not(BoolExpr::Var(0));
  problem.bad = BoolExpr::Var(1);

  const std::string formattedProblem = formatKInductionProblemForDebug(problem);

  EXPECT_NE(formattedProblem.find("state0_symbols: [FALSE, TRUE]"),
            std::string::npos);
  EXPECT_NE(formattedProblem.find("FALSE' = <invalid>"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ProofProblemDebugFormatsDeepExpressionsWithoutCallStackRecursion) {
  constexpr size_t kExpressionDepth = 60000;
  BoolExpr* property = BoolExpr::Var(2);
  for (size_t index = 0; index < kExpressionDepth; ++index) {
    property = BoolExpr::And(property, BoolExpr::Var(index + 3));
  }

  KInductionProblem problem;
  problem.description = "deep debug expression";
  problem.property = property;

  // Full PDR tracing formats production ASIC DAGs deeper than the native
  // thread stack. The formatter must traverse that exact DAG iteratively.
  const std::string formattedProblem = formatKInductionProblemForDebug(problem);

  EXPECT_NE(formattedProblem.find("property:"), std::string::npos);
  EXPECT_NE(formattedProblem.find("x2"), std::string::npos);
  EXPECT_NE(
      formattedProblem.find("x" + std::to_string(kExpressionDepth + 2)),
      std::string::npos);
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
       KInductionDualRailUsesSameFrameComplementRailEqualities) {
  KInductionProblem problem = makeDualRailComplementedOutputProblemForTest();

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionEngineBatchesWideOutputProofs) {
  KInductionProblem problem;
  for (size_t i = 0; i < 129; ++i) {
    const size_t symbol = 2 + i;
    problem.allSymbols.push_back(symbol);
    problem.observedOutputNames.push_back("out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(symbol));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(symbol));
  }

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  size_t baseChecks = 0;
  size_t pos = 0;
  const std::string needle = "SEC diag: k-induction base k=0 begin";
  while ((pos = stderrOutput.find(needle, pos)) != std::string::npos) {
    ++baseChecks;
    pos += needle.size();
  }

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(baseChecks, 5u);
  EXPECT_NE(stderrOutput.find("SEC diag: k-induction problem outputs=129"),
            std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionBatchDecisionLimitSettingKeepsProofSound) {
  KInductionProblem problem =
      makeDualRailLocalFalseInvariantProblemForTest(/*outputCount=*/2);

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  const ScopedEnvVar decisionLimit("KEPLER_SEC_KI_BATCH_DECISION_LIMIT", "0");
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // The batch decision limit is a performance guard for broad output slices.
  // Even an aggressively tiny limit must not change the proof result: the
  // engine can still split or prove smaller slices using only the encoded
  // transition relation.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionStepConstrainsResetInputsAfterBootstrap) {
  KInductionProblem problem;
  constexpr size_t state0 = 2;
  constexpr size_t state1 = 3;
  constexpr size_t reset = 4;
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.inputSymbols = {reset};
  problem.allSymbols = {state0, state1, reset};
  problem.transitions0 = {{state0, BoolExpr::Var(reset)}};
  problem.transitions1 = {{state1, BoolExpr::createFalse()}};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.property =
      makeEqualityExpr(BoolExpr::Var(state0), BoolExpr::Var(state1));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // The induction query starts after the boot/reset prefix.  If reset is left
  // unconstrained here, the proof can invent a later reset assertion and break
  // a property that is valid in the post-bootstrap SEC environment.
  EXPECT_TRUE(provesByInduction(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailLeafResourceLimitContinuesToNextFrontier) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  // The leaf limit bounds decisions, not propagation.  This tiny invariant
  // proof closes without decisions, so it must still be accepted instead of
  // being treated as a resource-limited UNKNOWN.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailDeferredBaseLeafContinuesAfterResourceLimitedStep) {
  KInductionProblem problem;
  constexpr size_t stateA = 2;
  constexpr size_t stateB = 3;
  BoolExpr* const badState =
      BoolExpr::And(BoolExpr::Var(stateA), BoolExpr::Not(BoolExpr::Var(stateB)));
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {stateA, stateB};
  problem.allSymbols = {stateA, stateB};
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(stateA)),
                    BoolExpr::Not(BoolExpr::Var(stateB)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0 = {
      {stateA, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateB, BoolExpr::createFalse()}};
  for (size_t i = 0; i < 2; ++i) {
    problem.observedOutputNames.push_back("out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(badState);
    problem.observedOutputExprs1.push_back(BoolExpr::createFalse());
  }
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  // Output slices share one concrete base check at the end.  If a capped leaf
  // step is UNKNOWN at k=1, KI must still try the k=2 obligation instead of
  // reporting the output uncovered; the final shared base check remains the
  // gate that makes the sliced proof a real SEC proof.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(result.bound, 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailDeferredBaseLeafContinuesAfterUnknownStep) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  constexpr size_t kXorStateSymbols = 4096;
  BoolExpr* xorCone = BoolExpr::Var(2);
  for (size_t symbol = 3; symbol < 2 + kXorStateSymbols; ++symbol) {
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    xorCone = BoolExpr::Xor(xorCone, BoolExpr::Var(symbol));
  }
  problem.state0Symbols.insert(problem.state0Symbols.begin(), 2);
  problem.allSymbols.insert(problem.allSymbols.begin(), 2);
  problem.observedOutputNames = {"out0", "out1"};
  problem.observedOutputExprs0 = {xorCone, xorCone};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "100");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "100");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(4);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // This intentionally hard residual has already lost the original output
  // count, as happens after SEC residual subsetting.  The rail-state surface
  // must still stop repeated resource-limited leaves, but only after allowing
  // later strict KI depths for outputs that first become inductive there.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(result.bound, 4u);
  EXPECT_NE(
      stderrOutput.find("resource-limited; deferred base continues"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find(
          "resource-limited deferred leaf limit reached; reporting inconclusive"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailDeferredBaseSmallLeafKeepsSearchingAfterUnknownStep) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  constexpr size_t kXorStateSymbols = 64;
  BoolExpr* xorCone = BoolExpr::Var(2);
  for (size_t symbol = 3; symbol < 2 + kXorStateSymbols; ++symbol) {
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    xorCone = BoolExpr::Xor(xorCone, BoolExpr::Var(symbol));
  }
  problem.state0Symbols.insert(problem.state0Symbols.begin(), 2);
  problem.allSymbols.insert(problem.allSymbols.begin(), 2);
  problem.observedOutputNames = {"out0", "out1"};
  problem.observedOutputExprs0 = {xorCone, xorCone};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "100");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "100");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Small residual surfaces are cheap enough to keep the full strict KI depth
  // search.  The wide-output early stop must not drop coverage for compact
  // designs that may close at a later k.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(result.bound, 3u);
  EXPECT_EQ(
      stderrOutput.find(
          "resource-limited deferred leaf limit reached; reporting inconclusive"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSmallLeafUsesNormalProofProfileByDefault) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.originalObservedOutputCount = 18;

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Compact designs such as GCD need an unbounded strict KI leaf by default;
  // otherwise every depth can become UNKNOWN before the useful induction
  // search.  They also keep the strict simple-path strengthening because the
  // rail-state surface is small enough for the loop-free clauses to pay off.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction simple path states="),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("profile_symbols=4096"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailInductionKeepsUnknownRailValuesInLegalStateDomain) {
  KInductionProblem problem;
  constexpr size_t mayBeOne = 2;
  constexpr size_t mayBeZero = 3;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {mayBeOne, mayBeZero};
  problem.allSymbols = {mayBeOne, mayBeZero};
  problem.totalStateCount = 2;
  problem.dualRailStatePairs = {DualRailSymbolPair{mayBeOne, mayBeZero}};
  problem.property = BoolExpr::Not(
      BoolExpr::And(BoolExpr::Var(mayBeOne), BoolExpr::Var(mayBeZero)));
  problem.bad = BoolExpr::Not(problem.property);

  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::NotProved);

  problem.bootstrapStateAssignments = {{mayBeOne, false}, {mayBeZero, true}};

  // Bootstrap assignments constrain only the reachable initial frontier.
  // Arbitrary induction frames still use the paper's legal ternary domain,
  // where 11 is X, so this property is not inductive.
  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::NotProved);

  problem.deferBaseCaseChecks = true;
  // Deferring a local base check changes neither that domain nor the result.
  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::NotProved);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSmallDeferredLeafUsesDirectProofProfileWithoutCdcl) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.deferBaseCaseChecks = true;
  problem.originalObservedOutputCount = 18;

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "100");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // A deferred small residual leaf may still use a decision cap when explicitly
  // requested.  Use the direct dual-rail proof profile, but reserve the
  // SAT-oriented direct-CDCL mix for genuinely large parent rail-state leaves.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("profile_symbols=4096"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("direct_cdcl=0"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailDeferredLeafUsesDirectProofProfileAfterSplitting) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.deferBaseCaseChecks = true;
  for (size_t index = 0; index < 300; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + index * 2, 1001 + index * 2});
  }

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "100");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Deferred leaves come from residual splitting; the large parent rail-state
  // surface should keep Kissat on the direct capped profile even when the
  // reduced strict-KI cone for this leaf is small.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("profile_symbols=4096"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("direct_cdcl=1"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DirectDualRailInductionAllowsKissatPropagationUnderLeafLimit) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");

  // The dual-rail leaf cap bounds decisions, but Kissat may still close tiny
  // UNSAT obligations by propagation before making any decision.
  EXPECT_TRUE(provesByInduction(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1));
  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::Proved);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SupportBoundedOutputBatchingKeepsModerateOutputSlicesTogether) {
  KInductionProblem problem;
  for (size_t i = 0; i < 40; ++i) {
    const size_t symbol = 2 + i;
    problem.allSymbols.push_back(symbol);
    problem.observedOutputNames.push_back("out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(symbol));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(symbol));
  }

  // PDR uses this moderate output-batch shape to avoid hundreds of identical
  // one-output proof attempts while still bounding each SAT cone by output
  // count and support.
  const auto batches =
      buildSupportBoundedOutputBatches(problem, OutputBatchingLimits{16, 512});

  ASSERT_EQ(batches.size(), 3u);
  EXPECT_EQ(batches[0], (std::pair<size_t, size_t>(0, 16)));
  EXPECT_EQ(batches[1], (std::pair<size_t, size_t>(16, 32)));
  EXPECT_EQ(batches[2], (std::pair<size_t, size_t>(32, 40)));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ConfigureOutputBatchPreservesImmutableSameDesignRelations) {
  KInductionProblem source;
  source.observedOutputNames = {"out0", "out1"};
  source.observedOutputExprs0 = {BoolExpr::Var(2), BoolExpr::Var(3)};
  source.observedOutputExprs1 = {BoolExpr::Var(4), BoolExpr::Var(5)};
  source.sameFrameStateEqualityPairs0 = {{10, 11}, {12, 13}};
  source.sameFrameStateEqualityPairs1 = {{20, 21}, {22, 23}};

  // A reusable batch starts as the source model. Selecting another output
  // range must replace only output-owned fields and retain model relations.
  KInductionProblem batch = source;
  configureOutputBatchProblem(batch, source, 0, 1);
  EXPECT_EQ(batch.sameFrameStateEqualityPairs0,
            source.sameFrameStateEqualityPairs0);
  EXPECT_EQ(batch.sameFrameStateEqualityPairs1,
            source.sameFrameStateEqualityPairs1);

  configureOutputBatchProblem(batch, source, 1, 2);
  EXPECT_EQ(batch.sameFrameStateEqualityPairs0,
            source.sameFrameStateEqualityPairs0);
  EXPECT_EQ(batch.sameFrameStateEqualityPairs1,
            source.sameFrameStateEqualityPairs1);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailOutputBatchingStartsWithModerateSharedConeSlices) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  for (size_t i = 0; i < 80; ++i) {
    problem.observedOutputNames.push_back("out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(10 + i));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(20 + i));
  }

  const auto batches = buildSupportBoundedOutputBatches(problem);

  ASSERT_EQ(batches.size(), 10u);
  // Dual-rail KI keeps small exact OR batches together, but avoids the old
  // 64-output starting slice that built a wide transition CNF before splitting.
  EXPECT_EQ(batches[0], (std::pair<size_t, size_t>(0, 8)));
  EXPECT_EQ(batches[9], (std::pair<size_t, size_t>(72, 80)));
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailHugeOutputBatchingStartsAtStrictLeaves) {
  KInductionProblem problem;
  constexpr size_t kSwervRailStatePairs = 45096;
  problem.usesDualRailStateEncoding = true;
  for (size_t i = 0; i < 16; ++i) {
    problem.observedOutputNames.push_back("huge_out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(10 + i));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(100 + i));
  }
  for (size_t i = 0; i < kSwervRailStatePairs; ++i) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + i * 2, 1001 + i * 2});
  }

  const auto batches = buildSupportBoundedOutputBatches(problem);

  // A Swerv/BP-sized rail surface should not first rebuild broad UNKNOWN
  // batches that immediately split.  Each returned slice is still a normal
  // one-output strict k-induction obligation.
  ASSERT_EQ(batches.size(), 16u);
  EXPECT_EQ(batches.front(), (std::pair<size_t, size_t>(0, 1)));
  EXPECT_EQ(batches.back(), (std::pair<size_t, size_t>(15, 16)));
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSmallOutputBatchingKeepsPublicConjunctionTogether) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  for (size_t i = 0; i < 18; ++i) {
    problem.observedOutputNames.push_back("out" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(10 + i));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(100 + i));
  }
  for (size_t i = 0; i < 70; ++i) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + i * 2, 1001 + i * 2});
  }

  const auto batches = buildSupportBoundedOutputBatches(problem);

  // GCD-sized dual-rail designs often close only when the public output
  // conjunction is the strict KI property.  Keep those compact surfaces as one
  // proof obligation instead of prematurely falling into residual leaves.
  ASSERT_EQ(batches.size(), 1u);
  EXPECT_EQ(batches[0], (std::pair<size_t, size_t>(0, 18)));
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualBatchingCarriesPublicHypothesisToSingleLeaf) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  residual.initialCondition = BoolExpr::simplify(initial);
  residual.initializedStateCount = kOutputs * 2;
  residual.totalStateCount = kOutputs * 2;
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);
  EXPECT_EQ(
      residual.lazyTransitions->dualRailResidualPublicOutputCount,
      kOutputs);
  EXPECT_EQ(
      residual.lazyTransitions->dualRailResidualPublicProperty,
      residual.property);

  KInductionProblem leaf = residual;
  configureOutputBatchProblem(leaf, residual, 0, 1);
  leaf.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Strategy residual batching may eventually call KI on one public output.  The
  // one-output leaf still proves a strict step from the remembered residual
  // public conjunction; no internal cross-design state equality is introduced.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("source_outputs=4"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailPublicHypothesisKeepsLeafBadPredicate) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  residual.state0Symbols = {10, 11};
  residual.state1Symbols = {20, 21};
  residual.allSymbols = {10, 11, 20, 21};
  residual.transitions0 = {
      {10, BoolExpr::Var(11)},
      {11, BoolExpr::createTrue()}};
  residual.transitions1 = {
      {20, BoolExpr::Var(21)},
      {21, BoolExpr::createFalse()}};
  residual.observedOutputNames = {"out0", "out1"};
  residual.observedOutputExprs0 = {BoolExpr::Var(10), BoolExpr::Var(11)};
  residual.observedOutputExprs1 = {BoolExpr::Var(20), BoolExpr::Var(21)};
  residual.property = BoolExpr::And(
      makeEqualityExpr(
          residual.observedOutputExprs0[0],
          residual.observedOutputExprs1[0]),
      makeEqualityExpr(
          residual.observedOutputExprs0[1],
          residual.observedOutputExprs1[1]));
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);

  KInductionProblem leaf = residual;
  configureOutputBatchProblem(leaf, residual, 0, 1);
  leaf.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The public conjunction is only the induction hypothesis.  The final-frame
  // bad predicate must remain the selected leaf output; otherwise this case
  // would try to prove the intentionally non-inductive second output too.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("source_outputs=2"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualBatchTriesStoredPublicHypothesisBeforeSplit) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
  }
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);

  KInductionProblem slice = residual;
  configureOutputBatchProblem(slice, residual, 0, 2);
  slice.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(slice, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The two-output slice needs the full four-output previous-frame public
  // hypothesis, but it can prove as a batch.  Trying that hypothesis before
  // recursive splitting avoids one-output proof repetition on AES-like groups.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("source_outputs=4"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("output range [0,1)"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSplitLeafSkipsLargerStoredFallbackAfterParentFailure) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
  }
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);

  KInductionProblem slice = residual;
  configureOutputBatchProblem(slice, residual, 0, 2);
  slice.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(slice, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The full residual public hypothesis is public-output-only and is tried at
  // the parent slice.  After that parent failure, children keep the smaller
  // parent public conjunction instead of replaying the wider stored fallback.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_NE(
      stderrOutput.find("range [0,2) source_outputs=4"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("range [0,1) source_outputs=2"),
      std::string::npos);
  EXPECT_EQ(
      detail::countTextOccurrences(stderrOutput, "source_outputs=4"),
      1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualBatchSkipsOversizedStoredPublicHypothesis) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(state0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(state1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
  }
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));
  residual.lazyTransitions->dualRailResidualPublicProperty = residual.property;
  residual.lazyTransitions->dualRailResidualPublicBad = residual.bad;
  residual.lazyTransitions->dualRailResidualPublicOutputCount = 129;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(residual, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Oversized remembered residual surfaces are not replayed at every recursive
  // AES-style batch.  The split still uses the smaller strict public
  // conjunction from the current four-output range.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(
      stderrOutput.find("source_outputs=129"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("range [0,2) source_outputs=4"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualLeafKeepsLocalPublicHypothesisWhenEnough) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(state0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(state1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
  }
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);

  KInductionProblem leaf = residual;
  configureOutputBatchProblem(leaf, residual, 0, 1);
  leaf.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // If the selected public output is already inductive by itself, keep that
  // smaller strict KI proof.  The full residual public conjunction remains only
  // a fallback for leaves that really need cross-output public history.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(
      stderrOutput.find("source_outputs=4"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualLeafCachesInconclusiveStoredPublicHypothesis) {
  KInductionProblem residual;
  residual.usesDualRailStateEncoding = true;
  residual.deferBaseCaseChecks = true;
  residual.lazyTransitions = std::make_shared<LazyTransitionStore>();
  constexpr size_t kOutputs = 4;
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    residual.state0Symbols.push_back(state0);
    residual.state1Symbols.push_back(state1);
    residual.allSymbols.push_back(state0);
    residual.allSymbols.push_back(state1);
    residual.transitions0.emplace_back(state0, BoolExpr::Var(state0));
    residual.transitions1.emplace_back(state1, BoolExpr::Var(state1));
    residual.observedOutputNames.push_back("out" + std::to_string(output));
    residual.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    residual.observedOutputExprs1.push_back(BoolExpr::Var(state1));
  }
  residual.property = BoolExpr::createTrue();
  for (size_t output = 0; output < residual.observedOutputExprs0.size(); ++output) {
    residual.property = BoolExpr::And(
        residual.property,
        makeEqualityExpr(
            residual.observedOutputExprs0[output],
            residual.observedOutputExprs1[output]));
  }
  residual.property = BoolExpr::simplify(residual.property);
  residual.bad = BoolExpr::simplify(BoolExpr::Not(residual.property));

  defaultOutputBatchingLimitsForProblem(residual);

  KInductionProblem leaf = residual;
  configureOutputBatchProblem(leaf, residual, 0, 1);
  leaf.deferBaseCaseChecks = true;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine firstEngine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto firstResult = firstEngine.run(0);
  KInductionEngine secondEngine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto secondResult = secondEngine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The cache is only for repeated UNKNOWN leaves.  It must suppress the second
  // identical strengthened KI attempt without marking the output equivalent.
  EXPECT_EQ(firstResult.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(secondResult.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(
      detail::countTextOccurrences(stderrOutput, "source_outputs=4"),
      1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailStandaloneLeafSkipsOversizedStoredPublicHypothesis) {
  KInductionProblem leaf;
  leaf.usesDualRailStateEncoding = true;
  leaf.deferBaseCaseChecks = true;
  leaf.lazyTransitions = std::make_shared<LazyTransitionStore>();
  leaf.state0Symbols = {10};
  leaf.state1Symbols = {20};
  leaf.allSymbols = {10, 20};
  leaf.transitions0 = {{10, BoolExpr::Not(BoolExpr::Var(10))}};
  leaf.transitions1 = {{20, BoolExpr::Var(20)}};
  leaf.observedOutputNames = {"out0"};
  leaf.observedOutputExprs0 = {BoolExpr::Var(10)};
  leaf.observedOutputExprs1 = {BoolExpr::Var(20)};
  leaf.property = makeEqualityExpr(
      leaf.observedOutputExprs0[0], leaf.observedOutputExprs1[0]);
  leaf.bad = BoolExpr::simplify(BoolExpr::Not(leaf.property));
  leaf.lazyTransitions->dualRailResidualPublicProperty = leaf.property;
  leaf.lazyTransitions->dualRailResidualPublicBad = leaf.bad;
  leaf.lazyTransitions->dualRailResidualPublicOutputCount = 129;

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The local leaf is still strict KI.  Once it is UNKNOWN, do not replay an
  // oversized remembered AES-style public conjunction for every standalone
  // residual bit; that only returns UNKNOWN and dominates the workflow.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(
      stderrOutput.find("source_outputs=129"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailResidualLeafCachesLocalInconclusiveAttempt) {
  KInductionProblem leaf;
  leaf.usesDualRailStateEncoding = true;
  leaf.deferBaseCaseChecks = true;
  leaf.lazyTransitions = std::make_shared<LazyTransitionStore>();
  leaf.state0Symbols = {10};
  leaf.state1Symbols = {20};
  leaf.allSymbols = {10, 20};
  leaf.transitions0 = {{10, BoolExpr::Not(BoolExpr::Var(10))}};
  leaf.transitions1 = {{20, BoolExpr::Var(20)}};
  leaf.observedOutputNames = {"out0"};
  leaf.observedOutputExprs0 = {BoolExpr::Var(10)};
  leaf.observedOutputExprs1 = {BoolExpr::Var(20)};
  leaf.property = makeEqualityExpr(
      leaf.observedOutputExprs0[0], leaf.observedOutputExprs1[0]);
  leaf.bad = BoolExpr::simplify(BoolExpr::Not(leaf.property));

  const ScopedEnvVar secDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine firstEngine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto firstResult = firstEngine.run(1);
  KInductionEngine secondEngine(leaf, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto secondResult = secondEngine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Reusing a cached inconclusive leaf is a runtime cache only.  The second
  // identical strict KI attempt should not rebuild the induction step, and it
  // must still report inconclusive rather than covered.
  EXPECT_EQ(firstResult.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(firstResult.bound, 1u);
  EXPECT_EQ(secondResult.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(secondResult.bound, 1u);
  EXPECT_EQ(
      detail::countTextOccurrences(stderrOutput, "k-induction step k=1 begin"),
      1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSmallOutputBatchKeepsSimplePathBeforeDirectCdclProfile) {
  KInductionProblem problem =
      makeDualRailLocalFalseInvariantProblemForTest(/*outputCount=*/2);
  problem.originalObservedOutputCount = 33;

  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto proofStatus = proveByInductionStatus(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      std::nullopt);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Tiny strict dual-rail batches still fit the simple-path strengthening. Do
  // not switch them to the direct-CDCL profile that is reserved for compact
  // cones after simple-path has been suppressed.
  EXPECT_EQ(proofStatus, InductionProofStatus::Proved);
  EXPECT_NE(
      stderrOutput.find("k-induction simple path states="),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("k-induction compact dual-rail direct profile"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("direct_cdcl=1"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailSmallDeferredLeafUsesSimplePathWithoutLimit) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.deferBaseCaseChecks = true;
  problem.originalObservedOutputCount = 33;

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto proofStatus = proveByInductionStatus(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      std::nullopt);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Shared-base output batching proves split leaves with the base check
  // deferred to the full output set.  When the reduced rail cone is small, keep
  // the strict simple-path strengthening instead of switching to the direct
  // no-simple-path CDCL profile.
  EXPECT_EQ(proofStatus, InductionProofStatus::Proved);
  EXPECT_NE(
      stderrOutput.find("k-induction simple path states="),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("k-induction compact dual-rail direct profile"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("direct_cdcl=1"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailWideOriginalSmallDeferredLeafUsesDefaultCap) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.deferBaseCaseChecks = true;
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES residual leaves inherit a wide original output count.  Keep those
  // leaves resource-bounded by default so a hard output cannot monopolize the
  // workflow, while still adding loop-free simple-path strengthening when the
  // reduced rail cone is below the state threshold.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("k-induction simple path states="),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailCompactOriginalLargeDeferredLeafKeepsDefaultCapForLargeRailSurface) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();
  problem.deferBaseCaseChecks = true;
  problem.originalObservedOutputCount = 34;
  for (size_t index = 0; index < 300; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + index * 2, 1001 + index * 2});
  }

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Keep the default cap once a deferred leaf inherits a large rail-state
  // surface.  The exact COI may be small, but removing this cap lets one hard
  // residual leaf dominate the dual-rail KI workflow.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_EQ(
      stderrOutput.find("compact deferred leaf exact coi unbounded"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailWideDeferredEightOutputBatchUsesDefaultCap) {
  KInductionProblem problem =
      makeDualRailLocalFalseInvariantProblemForTest(/*outputCount=*/8);
  problem.deferBaseCaseChecks = true;
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES residuals arrive as 8-output batches.  The first strict KI attempt is
  // deliberately bounded; if it is UNKNOWN, recursive splitting continues with
  // the parent public-output conjunction as the child induction hypothesis.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailWideSplitUsesParentPublicConjunctionHypothesis) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.originalObservedOutputCount = 129;
  constexpr size_t kOutputs = 8;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    problem.state0Symbols.push_back(state0);
    problem.state1Symbols.push_back(state1);
    problem.allSymbols.push_back(state0);
    problem.allSymbols.push_back(state1);
    problem.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    problem.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    problem.observedOutputNames.push_back("out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  problem.initialCondition = BoolExpr::simplify(initial);
  problem.initializedStateCount = kOutputs * 2;
  problem.totalStateCount = kOutputs * 2;
  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(
            problem.observedOutputExprs0[output],
            problem.observedOutputExprs1[output]));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Each output equality depends on another public output's previous equality.
  // Recursive splitting must therefore keep the parent public conjunction as
  // the induction hypothesis; no cross-design internal state equality is used.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("dual-rail public conjunction hypothesis"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailWideBatchesUseFullPublicConjunctionHypothesis) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  constexpr size_t kOutputs = 40;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    problem.state0Symbols.push_back(state0);
    problem.state1Symbols.push_back(state1);
    problem.allSymbols.push_back(state0);
    problem.allSymbols.push_back(state1);
    problem.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    problem.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    problem.observedOutputNames.push_back("out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  problem.initialCondition = BoolExpr::simplify(initial);
  problem.initializedStateCount = kOutputs * 2;
  problem.totalStateCount = kOutputs * 2;
  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(
            problem.observedOutputExprs0[output],
            problem.observedOutputExprs1[output]));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The first 8-output batch depends on equality of output 8 in the previous
  // frame, so per-batch history is not enough.  Use the full public SEC
  // conjunction as the strict KI hypothesis across initial batches.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("source_outputs=40"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailCompactPublicSurfaceRespectsLargeStateGuard) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  constexpr size_t kOutputs = 40;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    problem.state0Symbols.push_back(state0);
    problem.state1Symbols.push_back(state1);
    problem.allSymbols.push_back(state0);
    problem.allSymbols.push_back(state1);
    problem.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    problem.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    problem.observedOutputNames.push_back("out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  for (size_t index = 0; index < 3000; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + index * 2, 1001 + index * 2});
  }
  problem.initialCondition = BoolExpr::simplify(initial);
  problem.initializedStateCount = kOutputs * 2;
  problem.totalStateCount = kOutputs * 2;
  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(
            problem.observedOutputExprs0[output],
            problem.observedOutputExprs1[output]));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Large unrelated rail-state arrays make the full public hypothesis too
  // expensive to replay into every split.  Keep strict KI on the local batch
  // instead of rebuilding the 40-output source hypothesis.
  EXPECT_NE(result.status, KInductionStatus::Different);
  EXPECT_EQ(
      stderrOutput.find("source_outputs=40"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailCompactStoredPublicHypothesisRespectsLargeStateGuard) {
  KInductionProblem fullProblem;
  fullProblem.usesDualRailStateEncoding = true;
  constexpr size_t kOutputs = 40;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    fullProblem.state0Symbols.push_back(state0);
    fullProblem.state1Symbols.push_back(state1);
    fullProblem.allSymbols.push_back(state0);
    fullProblem.allSymbols.push_back(state1);
    fullProblem.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    fullProblem.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    fullProblem.observedOutputNames.push_back("out" + std::to_string(output));
    fullProblem.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    fullProblem.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  for (size_t index = 0; index < 3000; ++index) {
    fullProblem.dualRailStatePairs.push_back(
        DualRailSymbolPair{1000 + index * 2, 1001 + index * 2});
  }
  fullProblem.initialCondition = BoolExpr::simplify(initial);
  fullProblem.initializedStateCount = kOutputs * 2;
  fullProblem.totalStateCount = kOutputs * 2;
  fullProblem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < fullProblem.observedOutputExprs0.size(); ++output) {
    fullProblem.property = BoolExpr::And(
        fullProblem.property,
        makeEqualityExpr(
            fullProblem.observedOutputExprs0[output],
            fullProblem.observedOutputExprs1[output]));
  }
  fullProblem.property = BoolExpr::simplify(fullProblem.property);
  fullProblem.bad = BoolExpr::simplify(BoolExpr::Not(fullProblem.property));
  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->dualRailResidualPublicProperty = fullProblem.property;
  lazyTransitions->dualRailResidualPublicBad = fullProblem.bad;
  lazyTransitions->dualRailResidualPublicOutputCount = kOutputs;
  fullProblem.lazyTransitions = lazyTransitions;

  KInductionProblem batchProblem = fullProblem;
  configureOutputBatchProblem(batchProblem, fullProblem, 0, 8);
  batchProblem.lazyTransitions = lazyTransitions;

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(batchProblem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Strategy residual batching may remember the full public residual surface,
  // but large rail-state surfaces should not replay that stored hypothesis into
  // each split batch.  This preserves strict KI while avoiding AES-sized
  // rebuilds at every recursive level.
  EXPECT_NE(result.status, KInductionStatus::Different);
  EXPECT_EQ(
      stderrOutput.find("source_outputs=40"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailCompactSplitUsesPublicConjunctionHypothesis) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  constexpr size_t kOutputs = 2;
  BoolExpr* initial = BoolExpr::createTrue();
  for (size_t output = 0; output < kOutputs; ++output) {
    const size_t state0 = 10 + output;
    const size_t state1 = 100 + output;
    const size_t next0 = 10 + ((output + 1) % kOutputs);
    const size_t next1 = 100 + ((output + 1) % kOutputs);
    problem.state0Symbols.push_back(state0);
    problem.state1Symbols.push_back(state1);
    problem.allSymbols.push_back(state0);
    problem.allSymbols.push_back(state1);
    problem.transitions0.emplace_back(state0, BoolExpr::Var(next0));
    problem.transitions1.emplace_back(state1, BoolExpr::Var(next1));
    problem.observedOutputNames.push_back("out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(state0));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(state1));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state0)));
    initial = BoolExpr::And(initial, BoolExpr::Not(BoolExpr::Var(state1)));
  }
  problem.initialCondition = BoolExpr::simplify(initial);
  problem.initializedStateCount = kOutputs * 2;
  problem.totalStateCount = kOutputs * 2;
  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < problem.observedOutputExprs0.size(); ++output) {
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(
            problem.observedOutputExprs0[output],
            problem.observedOutputExprs1[output]));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Decomposed KI for a compact public conjunction is still strict KI: each
  // output depends on the other output's previous equality, so each split must
  // prove under the full public-output induction hypothesis and the engine
  // accepts only after all splits plus the shared base check.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("dual-rail public conjunction hypothesis"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailCompactBatchUsesBoundedStrictKInductionByDefault) {
  KInductionProblem problem =
      makeDualRailLocalFalseInvariantProblemForTest(/*outputCount=*/2);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "");
  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Compact multi-output batches may still prove by propagation, but their
  // first strict KI attempt must use the bounded batch path so hard conjunctions
  // split instead of monopolizing the workflow.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailBatchedKInductionChecksSharedBaseAfterSliceProofs) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  constexpr size_t invariantState = 3;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state, invariantState};
  problem.allSymbols = {state, invariantState};
  problem.initialCondition = BoolExpr::Var(state);
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{state, BoolExpr::Var(state)}};
  problem.observedOutputNames = {"out0", "out1"};
  problem.observedOutputExprs0 = {BoolExpr::Var(state), BoolExpr::Var(state)};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // The induction step can prove that "state is false" is invariant, but the
  // initial condition sets state true.  Batched dual-rail KI may defer per-slice
  // base checks for runtime, but it must still validate the shared full-output
  // base prefix before reporting equivalence.
  EXPECT_EQ(result.status, KInductionStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.witness->badFrame, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BaseCaseSolverFindsObservationOnlyProbeAtFrontier) {
  KInductionProblem problem;
  constexpr size_t input = 2;
  constexpr size_t state0 = 3;
  constexpr size_t state1 = 4;
  problem.environmentInputNames = {"in"};
  problem.inputSymbols = {input};
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.allSymbols = {input, state0, state1};
  problem.transitions0 = {{state0, BoolExpr::createFalse()}};
  problem.transitions1 = {{state1, BoolExpr::createFalse()}};
  problem.observedOutputNames = {"state_out", "probe"};
  problem.observedOutputExprs0 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::Var(input)};
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  EXPECT_FALSE(findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0).has_value());
  const auto witness = findBaseCounterexample(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1);

  // Shared KI base validation must include the proved frontier. This exact
  // top-output probe is invisible at k=0 under observation-only startup, but a
  // real SEC counterexample appears as soon as frame 1 is checked.
  ASSERT_TRUE(witness.has_value());
  EXPECT_EQ(witness->badFrame, 1u);
  ASSERT_EQ(witness->outputMismatches.size(), 1u);
  EXPECT_EQ(witness->outputMismatches[0].signal, "probe");
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionValidatesConcreteBaseAfterResetBootstrapProof) {
  KInductionProblem problem;
  constexpr size_t reset = 2;
  constexpr size_t state0 = 3;
  constexpr size_t state1 = 4;
  constexpr size_t probe = 5;
  problem.environmentInputNames = {"reset", "probe"};
  problem.inputSymbols = {reset, probe};
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.allSymbols = {reset, state0, state1, probe};
  problem.totalStateCount = 2;
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{state0, false}};
  problem.transitions0 = {{state0, BoolExpr::createFalse()}};
  problem.transitions1 = {{state1, BoolExpr::createFalse()}};
  problem.observedOutputNames = {"probe_out"};
  problem.observedOutputExprs0 = {BoolExpr::createFalse()};
  problem.observedOutputExprs1 = {BoolExpr::Var(probe)};
  problem.bad = BoolExpr::Var(probe);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // A strengthened induction certificate is not itself a SEC result. The
  // engine must still check the real top-output base predicate at the proved
  // frontier before reporting equivalence.
  EXPECT_EQ(result.status, KInductionStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.witness->badFrame, 1u);
  ASSERT_EQ(result.witness->outputMismatches.size(), 1u);
  EXPECT_EQ(result.witness->outputMismatches[0].signal, "probe_out");
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailBatchedKInductionAllowsKissatPropagationUnderBatchLimit) {
  KInductionProblem problem =
      makeDualRailLocalFalseInvariantProblemForTest(/*outputCount=*/2);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "100");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "100");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The cap is still applied to decisions, but this small invariant is solved
  // by propagation.  Keep the test as a guard that the KISSAT dual-rail path
  // remains a valid proof route after removing the old CaDiCaL detour.
  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("k-induction direct dual-rail capped proof profile"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailKInductionUsesBatchedLazyTransitionSupport) {
  KInductionProblem problem;
  constexpr size_t localState = 2;
  constexpr size_t railOne = 10;
  constexpr size_t railZero = 11;
  constexpr size_t targetBase = 100;
  constexpr size_t targetPairs = 12;
  BoolExpr* localNext = BoolExpr::Var(localState);

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->dualRailStateByLocalSymbolByDesign[0].emplace(
      localState, DualRailSymbolPair{railOne, railZero});
  lazyTransitions->sourceByStateSymbol.emplace(
      railOne, LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailOne});
  lazyTransitions->sourceByStateSymbol.emplace(
      railZero, LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailZero});

  problem.usesDualRailStateEncoding = true;
  problem.lazyTransitions = lazyTransitions;
  problem.state0Symbols = {railOne, railZero};
  problem.dualRailStatePairs.push_back(DualRailSymbolPair{railOne, railZero});
  problem.initialCondition =
      BoolExpr::And(BoolExpr::Var(railOne), BoolExpr::Not(BoolExpr::Var(railZero)));
  problem.property =
      BoolExpr::And(BoolExpr::Var(railOne), BoolExpr::Not(BoolExpr::Var(railZero)));
  for (size_t index = 0; index < targetPairs; ++index) {
    const size_t targetOne = targetBase + index * 2;
    const size_t targetZero = targetOne + 1;
    problem.state0Symbols.push_back(targetOne);
    problem.state0Symbols.push_back(targetZero);
    problem.dualRailStatePairs.push_back(DualRailSymbolPair{targetOne, targetZero});
    lazyTransitions->sourceByStateSymbol.emplace(
        targetOne,
        LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailOne});
    lazyTransitions->sourceByStateSymbol.emplace(
        targetZero,
        LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailZero});
    problem.initialCondition = BoolExpr::And(
        problem.initialCondition,
        BoolExpr::And(BoolExpr::Var(targetOne),
                      BoolExpr::Not(BoolExpr::Var(targetZero))));
    problem.property = BoolExpr::And(
        problem.property,
        BoolExpr::And(BoolExpr::Var(targetOne),
                      BoolExpr::Not(BoolExpr::Var(targetZero))));
  }
  problem.allSymbols = problem.state0Symbols;
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, KInductionStatus::Equivalent);
  // KI/IMC COI construction should use the resolver's frame-wide support
  // collector.  The old per-target path filled this lazy cache once for every
  // rail target, which dominated Ibex dual-rail IMC before SAT solving.
  EXPECT_TRUE(lazyTransitions->supportByStateSymbol.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionReusesTransitionSupportCacheAcrossDepths) {
  KInductionProblem problem = makeDualRailLocalFalseInvariantProblemForTest();

  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::Proved);
  const auto firstCache = problem.inductionTransitionSupportCache;
  ASSERT_TRUE(firstCache);

  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          2,
          std::nullopt),
      InductionProofStatus::Proved);

  // Increasing-k retries on the same strict KI problem should reuse the exact
  // transition-target support cache; only the SAT formula is rebuilt for the
  // new depth.
  EXPECT_EQ(problem.inductionTransitionSupportCache.get(), firstCache.get());
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionTracksExactTransitionFrameSupport) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  constexpr size_t complementedState = 3;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state, complementedState};
  problem.allSymbols = {state, complementedState};
  problem.complementedStatePairs0 = {{state, complementedState}};
  problem.transitions0 = {{state, BoolExpr::createFalse()}};
  problem.property = BoolExpr::Not(BoolExpr::Var(state));
  problem.bad = BoolExpr::Var(state);

  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto proofStatus = proveByInductionStatus(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      std::nullopt);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // A constant transition still has a target equality, but it reads no frame
  // leaves.  KI uses that exact frame support to avoid building wide leaf maps
  // for unrelated symbols while preserving the same transition obligation.
  EXPECT_EQ(proofStatus, InductionProofStatus::Proved);
  EXPECT_NE(
      stderrOutput.find("transition_targets=1 transition_support=0"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionUsesTransitionNodeCountReserveHint) {
  KInductionProblem problem;
  constexpr size_t state0 = 2;
  constexpr size_t state1 = 3;
  constexpr size_t inputA = 4;
  constexpr size_t inputB = 5;
  BoolExpr* sharedCone =
      makeRepeatedSmallSupportCone(inputA, inputB, /*depth=*/24);
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.inputSymbols = {inputA, inputB};
  problem.allSymbols = {state0, state1, inputA, inputB};
  problem.transitions0 = {{state0, sharedCone}};
  problem.transitions1 = {{state1, sharedCone}};
  problem.property = makeEqualityExpr(BoolExpr::Var(state0), BoolExpr::Var(state1));
  problem.bad = BoolExpr::Not(problem.property);

  const ScopedEnvVar coiDiag("KEPLER_SEC_KI_COI_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto proofStatus = proveByInductionStatus(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      std::nullopt);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The strict KI transition relation still uses the exact frame support for
  // leaves, but the encoder now reserves from the transition DAG node count so
  // deep small-support cones do not repeatedly grow their SAT-side cache.
  EXPECT_EQ(proofStatus, InductionProofStatus::Proved);
  EXPECT_NE(
      stderrOutput.find("k-induction transition node reserve"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("support=2"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailBatchedKInductionReturnsInconclusiveOnDecisionBudget) {
  KInductionProblem problem;
  constexpr size_t state0 = 2;
  constexpr size_t state1 = 3;
  constexpr size_t freeInput0 = 4;
  constexpr size_t freeInput1 = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.inputSymbols = {freeInput0, freeInput1};
  problem.allSymbols = {state0, state1, freeInput0, freeInput1};
  problem.transitions0 = {{state0, BoolExpr::Var(freeInput0)}};
  problem.transitions1 = {{state1, BoolExpr::Var(freeInput1)}};
  problem.observedOutputNames = {"out0", "out1"};
  problem.observedOutputExprs0 = {
      BoolExpr::Var(state0), BoolExpr::Not(BoolExpr::Var(state0))};
  problem.observedOutputExprs1 = {
      BoolExpr::Var(state1), BoolExpr::Not(BoolExpr::Var(state1))};
  problem.property = BoolExpr::And(
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]),
      makeEqualityExpr(problem.observedOutputExprs0[1],
                       problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // Hard residual dual-rail outputs are allowed to stay uncovered, but they
  // must not block the workflow once the configured SAT decision budget is hit.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
}

TEST_F(SequentialEquivalenceStrategyTests,
       DualRailHugeDeferredLeafStopsAtFirstResourceLimit) {
  KInductionProblem problem;
  constexpr size_t kStatePairs = 45096;
  constexpr size_t kFirstOneRail = 2;
  constexpr size_t kSecondOneRail = 4;
  problem.usesDualRailStateEncoding = true;
  problem.deferBaseCaseChecks = true;
  problem.state0Symbols.reserve(kStatePairs * 2);
  problem.allSymbols.reserve(kStatePairs * 2);
  for (size_t index = 0; index < kStatePairs; ++index) {
    const size_t oneRail = 2 + index * 2;
    const size_t zeroRail = oneRail + 1;
    problem.state0Symbols.push_back(oneRail);
    problem.state0Symbols.push_back(zeroRail);
    problem.allSymbols.push_back(oneRail);
    problem.allSymbols.push_back(zeroRail);
    problem.dualRailStatePairs.push_back(DualRailSymbolPair{oneRail, zeroRail});
    problem.transitions0.emplace_back(oneRail, BoolExpr::Var(oneRail));
    problem.transitions0.emplace_back(zeroRail, BoolExpr::Var(zeroRail));
  }
  problem.observedOutputExprs0 = {BoolExpr::Var(kFirstOneRail)};
  problem.observedOutputExprs1 = {BoolExpr::Var(kSecondOneRail)};
  problem.property =
      makeEqualityExpr(problem.observedOutputExprs0[0],
                       problem.observedOutputExprs1[0]);
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(8);

  // A capped UNKNOWN on a Swerv/BP-sized deferred residual leaf is not a proof.
  // It should be reported as uncovered immediately instead of rebuilding the
  // same huge strict-KI obligation at larger k values.
  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(result.bound, 1u);
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
       KInductionEngineDiagnosticsCoverCounterexampleAndInconclusivePaths) {
  {
    KInductionProblem problem;
    problem.allSymbols = {2};
    problem.bad = BoolExpr::createTrue();
    problem.property = BoolExpr::createFalse();
    problem.inductionProperty = problem.property;
    problem.inductionBad = problem.bad;

    const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
    testing::internal::CaptureStderr();
    KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(3);
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    EXPECT_EQ(result.status, KInductionStatus::Different);
    EXPECT_EQ(result.bound, 0u);
    EXPECT_NE(
        stderrOutput.find("SEC diag: k-induction base k=0 found cex"),
        std::string::npos);
  }

  {
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

    const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
    testing::internal::CaptureStderr();
    KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(3);
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    EXPECT_EQ(result.status, KInductionStatus::Different);
    EXPECT_EQ(result.bound, 1u);
    EXPECT_NE(
        stderrOutput.find("SEC diag: k-induction base k=1 found cex"),
        std::string::npos);
  }

  {
    const auto problem = buildLinearChainSecProblem(6);

    const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
    testing::internal::CaptureStderr();
    KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(4);
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
    EXPECT_EQ(result.bound, 4u);
    EXPECT_NE(
        stderrOutput.find("SEC diag: k-induction step k=4 inconclusive"),
        std::string::npos);
  }
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
       IMCEngineUsesObservationOnlyFrontierWithoutExplicitInit) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::Var(2));
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcDerivesRelationWithoutCrossDesignStateAssumptions) {
  const KInductionProblem problem = buildCraigResetSecProblem(true);

  CraigInterpolatingModelChecker checker(problem);
  const CraigImcResult result = checker.run(4);

  // The two state symbols are deliberately unrelated in the problem. Their
  // equality may appear only as a consequence of the reset and transition
  // clauses in CaDiCaL's UNSAT proof, never as an internal-name assumption.
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
  EXPECT_GE(result.iterations, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSimplifiesRedundantInterpolantClauses) {
  InterpolantRegion region;
  region.type = InterpolantRegion::Type::Normal;
  region.root = {true, 10, true};
  region.definitionLiterals = {
      {true, 10, true},
      {true, 11, true},
      {true, 10, true},
      {true, 10, true},
      {true, 12, true},
      {true, 12, false}};
  region.definitionClauseEnds = {2, 3, 4, 6};

  const InterpolantRegion simplified =
      simplifyCraigInterpolantRegion(std::move(region));

  // This is the IMC-local version of McMillan's redundant-interpolant cleanup:
  // remove duplicate clauses, tautologies, and clauses subsumed by a smaller
  // clause before future Craig iterations keep re-instantiating them.
  ASSERT_EQ(simplified.type, InterpolantRegion::Type::Normal);
  EXPECT_EQ(simplified.definitionClauseEnds.size(), 1u);
  ASSERT_EQ(simplified.definitionLiterals.size(), 1u);
  EXPECT_TRUE(simplified.definitionLiterals.front().isState);
  EXPECT_EQ(simplified.definitionLiterals.front().index, 10u);
  EXPECT_TRUE(simplified.definitionLiterals.front().positive);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcCompactsRedundantReachableRegions) {
  constexpr size_t state = 10;
  KInductionProblem problem;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  const std::unordered_set<size_t> trackedStates = {state};
  const std::vector<InterpolantRegion> helperRegions;
  std::vector<InterpolantRegion> reachableRegions = {
      makeStateLiteralCraigRegion(state, true),
      makeStateLiteralCraigRegion(state, true),
      makeStateLiteralCraigRegion(state, false)};

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  const size_t removed = compactCraigReachableRegions(
      problem,
      trackedStates,
      helperRegions,
      reachableRegions,
      /*compactionStart=*/1,
      /*candidateLimit=*/4);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The first x region is contained in the union of the duplicate x and !x.
  // Removing only one region per pass keeps duplicate regions from justifying
  // each other's deletion in the same compaction sweep.
  EXPECT_EQ(removed, 1u);
  ASSERT_EQ(reachableRegions.size(), 2u);
  EXPECT_NE(
      stderrOutput.find("imc Craig compacted reachable regions 3->2"),
      std::string::npos)
      << stderrOutput;
  testing::internal::CaptureStderr();
  EXPECT_EQ(
      compactCraigReachableRegions(
          problem,
          trackedStates,
          helperRegions,
          reachableRegions,
          /*compactionStart=*/1,
          /*candidateLimit=*/4),
      0u);
  testing::internal::GetCapturedStderr();
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesMcMillanFixedpointContainment) {
  const KInductionProblem problem = buildCraigResetSecProblem(true);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(problem);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
  // max_k is the Craig lookahead. The paper's single loop may still execute
  // several Q := Q OR I passes before that lookahead changes.
  EXPECT_NE(
      stderrOutput.find("lookahead=1 q_pass=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig fixedpoint containment"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig incremental inductiveness"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcReportsGrowthBudgetBeforeFixedpointContainment) {
  const KInductionProblem problem = buildCraigResetSecProblem(true);
  CraigImcOptions options;
  options.growthBudget.enabled = true;
  options.growthBudget.maxQExpansionPass = 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  const std::unordered_set<size_t> initialTrackedStates(
      problem.state0Symbols.begin(), problem.state0Symbols.end());

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The budget is not an equivalence result. It is a caller-visible signal that
  // a large dual-rail output batch should be split and retried as strict IMC.
  EXPECT_EQ(result.status, CraigImcStatus::BudgetExceeded);
  EXPECT_NE(
      stderrOutput.find("imc Craig growth budget exceeded reason=q_pass"),
      std::string::npos)
      << stderrOutput;
}

KInductionProblem buildCraigAuxiliaryConstantProblemForTest() {
  constexpr size_t stableState = 2;
  constexpr size_t unstableState = 3;
  constexpr size_t input = 4;

  KInductionProblem problem;
  problem.inputSymbols = {input};
  problem.state0Symbols = {stableState, unstableState};
  problem.allSymbols = {stableState, unstableState, input};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {
      {stableState, false},
      {unstableState, false}};
  problem.transitions0 = {
      {stableState, BoolExpr::createFalse()},
      {unstableState, BoolExpr::Var(input)}};
  problem.bad = BoolExpr::Var(stableState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigAuxiliaryConstantGatedConeProblemForTest(
    size_t coneStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t guardState = 3;
  constexpr size_t firstConeState = 4;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, guardState};
  problem.allSymbols = {badState, guardState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {guardState, false}};

  std::vector<size_t> coneStates;
  coneStates.reserve(coneStateCount);
  for (size_t index = 0; index < coneStateCount; ++index) {
    const size_t symbol = firstConeState + index;
    coneStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }

  problem.transitions0 = {
      {badState,
       BoolExpr::And(BoolExpr::Var(guardState),
                     makeConjunctionOfVars(coneStates))},
      {guardState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigResetGatedConeProblemForTest(
    size_t coneStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t resetInput = 3;
  constexpr size_t firstConeState = 4;

  KInductionProblem problem;
  problem.inputSymbols = {resetInput};
  problem.resetBootstrapInputs = {{resetInput, true}};
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState};
  problem.allSymbols = {badState, resetInput};
  problem.bootstrapStateAssignments = {{badState, false}};

  std::vector<size_t> coneStates;
  coneStates.reserve(coneStateCount);
  for (size_t index = 0; index < coneStateCount; ++index) {
    const size_t symbol = firstConeState + index;
    coneStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }

  problem.transitions0 = {
      {badState,
       BoolExpr::And(BoolExpr::Var(resetInput),
                     makeConjunctionOfVars(coneStates))}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigAuxiliaryEqualityProblemForTest() {
  constexpr size_t lhsState = 2;
  constexpr size_t rhsState = 3;
  constexpr size_t unrelatedState = 4;
  constexpr size_t input = 5;

  KInductionProblem problem;
  problem.inputSymbols = {input};
  problem.state0Symbols = {lhsState, rhsState, unrelatedState};
  problem.allSymbols = {lhsState, rhsState, unrelatedState, input};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {
      {lhsState, false},
      {rhsState, false},
      {unrelatedState, false}};
  problem.transitions0 = {
      {lhsState, BoolExpr::Var(input)},
      {rhsState, BoolExpr::Var(input)},
      {unrelatedState, BoolExpr::Not(BoolExpr::Var(input))}};
  problem.property =
      makeEqualityExpr(BoolExpr::Var(lhsState), BoolExpr::Var(rhsState));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigCrossDesignAuxiliaryEqualityProblemForTest() {
  constexpr size_t state0 = 2;
  constexpr size_t state1 = 3;
  constexpr size_t input = 4;

  KInductionProblem problem;
  problem.inputSymbols = {input};
  problem.state0Symbols = {state0};
  problem.state1Symbols = {state1};
  problem.allSymbols = {state0, state1, input};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {{state0, false}, {state1, false}};
  problem.transitions0 = {{state0, BoolExpr::Var(input)}};
  problem.transitions1 = {{state1, BoolExpr::Var(input)}};
  problem.property =
      makeEqualityExpr(BoolExpr::Var(state0), BoolExpr::Var(state1));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount =
      problem.state0Symbols.size() + problem.state1Symbols.size();
  return problem;
}

KInductionProblem buildCraigLocalAuxiliaryEqualityProblemForTest(
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t lhsState = 3;
  constexpr size_t rhsState = 4;
  constexpr size_t input = 5;
  constexpr size_t firstPaddingState = 6;

  KInductionProblem problem;
  problem.inputSymbols = {input};
  problem.state0Symbols = {badState, lhsState, rhsState};
  problem.allSymbols = {badState, lhsState, rhsState, input};
  problem.resetBootstrapCycles = 1;
  problem.bootstrapStateAssignments = {
      {badState, false},
      {lhsState, false},
      {rhsState, false}};
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  // lhsState and rhsState are not constants, but they start equal and share the
  // same next-state function.  With the oversized padding the global auxiliary
  // pass must skip, so only local equality mining can keep badState out without
  // importing the lhs/rhs support into the Craig projection.
  problem.transitions0 = {
      {badState,
       BoolExpr::Not(makeEqualityExpr(
           BoolExpr::Var(lhsState), BoolExpr::Var(rhsState)))},
      {lhsState, BoolExpr::Var(input)},
      {rhsState, BoolExpr::Var(input)}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigAuxiliaryCandidateGuardProblemForTest(
    size_t candidateCount) {
  constexpr size_t firstState = 2;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols.reserve(candidateCount);
  problem.allSymbols.reserve(candidateCount);
  problem.bootstrapStateAssignments.reserve(candidateCount);
  for (size_t index = 0; index < candidateCount; ++index) {
    const size_t symbol = firstState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  // Keep only the first state in the property so the test exercises the
  // auxiliary-candidate guard without building a large Craig proof.
  problem.transitions0 = {{firstState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(firstState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigLocalAuxiliaryInvariantProblemForTest(
    size_t noiseStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t guardState = 3;
  constexpr size_t firstNoiseState = 4;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, guardState};
  problem.allSymbols = {badState, guardState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {guardState, false}};
  for (size_t index = 0; index < noiseStateCount; ++index) {
    const size_t symbol = firstNoiseState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  // The property state can only become bad through guardState.  The oversized
  // bootstrap vector forces global auxiliary mining to skip, so the proof needs
  // the local Craig projection to prove guardState is a stable constant.
  problem.transitions0 = {
      {badState, BoolExpr::Var(guardState)},
      {guardState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigAuxiliaryReuseAcrossOutputsProblemForTest(
    size_t noiseStateCount) {
  constexpr size_t outputCount = 60;
  constexpr size_t firstBadState = 2;
  constexpr size_t guardState = firstBadState + outputCount;
  constexpr size_t firstNoiseState = guardState + 1;

  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  for (size_t output = 0; output < outputCount; ++output) {
    const size_t badState = firstBadState + output;
    problem.state0Symbols.push_back(badState);
    problem.allSymbols.push_back(badState);
    problem.bootstrapStateAssignments.push_back({badState, false});
    problem.transitions0.emplace_back(badState, BoolExpr::Var(guardState));
    problem.observedOutputNames.push_back(
        "aux_reuse_out" + std::to_string(output));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(badState));
    problem.observedOutputExprs1.push_back(BoolExpr::createFalse());
  }
  problem.state0Symbols.push_back(guardState);
  problem.allSymbols.push_back(guardState);
  problem.bootstrapStateAssignments.push_back({guardState, false});
  for (size_t index = 0; index < noiseStateCount; ++index) {
    const size_t symbol = firstNoiseState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  // All outputs depend on the same stable guard.  The first output proves the
  // guard as a local auxiliary constant after global mining skips; the second
  // output should seed that fact from the saved Craig helper instead of mining
  // it again.
  problem.transitions0.emplace_back(guardState, BoolExpr::createFalse());
  problem.property = BoolExpr::createTrue();
  for (size_t output = 0; output < outputCount; ++output) {
    problem.property = BoolExpr::And(
        problem.property,
        makeEqualityExpr(problem.observedOutputExprs0[output],
                         problem.observedOutputExprs1[output]));
  }
  problem.property = BoolExpr::simplify(problem.property);
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigComplementLocalAuxiliaryProblemForTest(
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t primaryRailState = 3;
  constexpr size_t complementRailState = 4;
  constexpr size_t firstPaddingState = 5;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, primaryRailState, complementRailState};
  problem.allSymbols = {badState, primaryRailState, complementRailState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {primaryRailState, false},
      {complementRailState, true}};
  problem.complementedStatePairs0.push_back(
      {primaryRailState, complementRailState});
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  // The bad transition only mentions the complement rail, while only the
  // primary rail has an explicit next-state function.  Local auxiliary mining
  // must use same-design Q/QN semantics to prove and apply complement constants
  // without growing the Craig projection.
  problem.transitions0 = {
      {badState, BoolExpr::Not(BoolExpr::Var(complementRailState))},
      {primaryRailState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigSelfPreservingLocalAuxiliaryPriorityProblemForTest(
    size_t noiseStateCount,
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t firstNoiseState = 3;

  const size_t guardState = firstNoiseState + noiseStateCount;
  const size_t firstPaddingState = guardState + 1;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState};
  problem.allSymbols = {badState};
  problem.bootstrapStateAssignments = {{badState, false}};

  std::vector<size_t> noiseStates;
  noiseStates.reserve(noiseStateCount);
  for (size_t index = 0; index < noiseStateCount; ++index) {
    const size_t symbol = firstNoiseState + index;
    noiseStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, index == 1});
  }
  problem.state0Symbols.push_back(guardState);
  problem.allSymbols.push_back(guardState);
  problem.bootstrapStateAssignments.push_back({guardState, false});
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  if (noiseStates.size() >= 2) {
    problem.complementedStatePairs0.push_back(
        {noiseStates[0], noiseStates[1]});
  }

  BoolExpr* unreachableNoiseCone = BoolExpr::createFalse();
  if (noiseStates.size() >= 2) {
    unreachableNoiseCone =
        BoolExpr::And(BoolExpr::Var(noiseStates[0]), BoolExpr::Var(noiseStates[1]));
    for (size_t index = 2; index < noiseStates.size(); ++index) {
      unreachableNoiseCone =
          BoolExpr::And(unreachableNoiseCone, BoolExpr::Var(noiseStates[index]));
    }
  }

  // guardState has a high symbol number and would be outside a simple
  // first-1024 slice.  The noise cone keeps the support huge but is impossible
  // under same-design complement semantics, so proving guardState=0 is the
  // useful local auxiliary fact.
  problem.transitions0 = {
      {badState, BoolExpr::Or(BoolExpr::Var(guardState), unreachableNoiseCone)},
      {guardState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigLargeLocalAuxiliaryInvariantProblemForTest(
    size_t supportNoiseCount,
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t guardState = 3;
  constexpr size_t firstNoiseState = 4;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, guardState};
  problem.allSymbols = {badState, guardState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {guardState, false}};

  std::vector<size_t> noiseStates;
  noiseStates.reserve(supportNoiseCount);
  for (size_t index = 0; index < supportNoiseCount; ++index) {
    const size_t symbol = firstNoiseState + index;
    noiseStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  const size_t firstPaddingState = firstNoiseState + supportNoiseCount;
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  BoolExpr* unreachableNoiseCone = BoolExpr::createFalse();
  if (noiseStates.size() >= 2) {
    problem.complementedStatePairs0.push_back(
        {noiseStates[0], noiseStates[1]});
    unreachableNoiseCone =
        BoolExpr::And(BoolExpr::Var(noiseStates[0]), BoolExpr::Var(noiseStates[1]));
    for (size_t index = 2; index < noiseStates.size(); ++index) {
      unreachableNoiseCone =
          BoolExpr::And(unreachableNoiseCone, BoolExpr::Var(noiseStates[index]));
    }
  }

  // The large support cone forces local auxiliary mining through its bounded
  // candidate selector.  Only guardState has a transition proof; the noise
  // terms keep the support large without needing thousands of useful constants.
  problem.transitions0 = {
      {badState, BoolExpr::Or(BoolExpr::Var(guardState), unreachableNoiseCone)},
      {guardState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigKnownSupportLocalAuxiliaryProblemForTest(
    size_t supportStateCount,
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t guardState = 3;
  constexpr size_t firstSupportState = 4;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, guardState};
  problem.allSymbols = {badState, guardState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {guardState, false}};

  std::vector<size_t> supportStates;
  supportStates.reserve(supportStateCount);
  for (size_t index = 0; index < supportStateCount; ++index) {
    const size_t symbol = firstSupportState + index;
    supportStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  const size_t firstPaddingState = firstSupportState + supportStateCount;
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  BoolExpr* knownFalseSupportCone = makeConjunctionOfVars(supportStates);
  problem.transitions0.emplace_back(
      badState,
      BoolExpr::Or(BoolExpr::Var(guardState), knownFalseSupportCone));
  problem.transitions0.emplace_back(guardState, knownFalseSupportCone);
  for (const size_t symbol : supportStates) {
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }

  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigLocalAuxiliaryRetryLimitProblemForTest(
    size_t supportStateCount,
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t firstSupportState = 3;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState};
  problem.allSymbols = {badState};
  problem.bootstrapStateAssignments = {{badState, false}};

  std::vector<size_t> supportStates;
  supportStates.reserve(supportStateCount);
  for (size_t index = 0; index < supportStateCount; ++index) {
    const size_t symbol = firstSupportState + index;
    supportStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }
  const size_t firstPaddingState = firstSupportState + supportStateCount;
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
  }

  BoolExpr* supportCanMakeBad = BoolExpr::createFalse();
  for (const size_t symbol : supportStates) {
    supportCanMakeBad =
        BoolExpr::Or(supportCanMakeBad, BoolExpr::Var(symbol));
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }
  problem.transitions0.emplace_back(
      badState, BoolExpr::simplify(supportCanMakeBad));

  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

KInductionProblem buildCraigLargeProjectionAuxiliarySuppressionProblemForTest(
    size_t paddingStateCount) {
  constexpr size_t badState = 2;
  constexpr size_t guardState = 3;
  constexpr size_t firstPaddingState = 4;

  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  problem.state0Symbols = {badState, guardState};
  problem.allSymbols = {badState, guardState};
  problem.bootstrapStateAssignments = {
      {badState, false},
      {guardState, false}};
  for (size_t index = 0; index < paddingStateCount; ++index) {
    const size_t symbol = firstPaddingState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }

  // Keep the useful proof tiny while the tracked projection itself is large.
  // The stable guard still creates a non-empty auxiliary set, while the
  // trivially safe property lets the test focus on query-clause suppression
  // instead of depending on auxiliaries for the proof.
  problem.transitions0 = {
      {badState, BoolExpr::Var(guardState)},
      {guardState, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::createFalse();
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();
  return problem;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcAuxiliaryConstantsRequireTransitionProof) {
  const KInductionProblem problem = buildCraigAuxiliaryConstantProblemForTest();

  const ScopedEnvVar enableAux("KEPLER_SEC_IMC_AUX_INVARIANTS", "1");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(problem);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Reset/bootstrap constants are not trusted blindly. Only stableState is
  // kept because its transition proves the same value; unstableState can follow
  // an unconstrained input and must not become an auxiliary invariant.
  EXPECT_NE(
      stderrOutput.find("imc Craig auxiliary constants=1 candidates=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcAuxiliaryEqualitiesStayWithinOneDesign) {
  const KInductionProblem problem =
      buildCraigCrossDesignAuxiliaryEqualityProblemForTest();
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The two states have the same reset value and next expression, but they are
  // in different designs. Craig IMC must not promote that into a relation.
  EXPECT_NE(
      stderrOutput.find("imc Craig auxiliary constants=0 candidates=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("equalities=0 equality_candidates=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcPromotesTransitionProvenSameDesignEqualities) {
  const KInductionProblem problem = buildCraigAuxiliaryEqualityProblemForTest();
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // lhsState and rhsState are not constants, but their equality is inductive.
  // unrelatedState shares the reset value but not the next function.
  EXPECT_NE(
      stderrOutput.find("imc Craig auxiliary constants=0 candidates=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("equalities=1 equality_candidates=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcOptionsPromoteAuxiliaryConstantsWithoutEnv) {
  const KInductionProblem problem = buildCraigAuxiliaryConstantProblemForTest();
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Large dual-rail IMC uses this option instead of relying on a process-wide
  // env var. The constant still has to be transition-preserved before use.
  EXPECT_NE(
      stderrOutput.find("imc Craig auxiliary constants=1 candidates=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSubstitutesAuxiliaryConstantsBeforeTransitionEncoding) {
  const KInductionProblem problem =
      buildCraigAuxiliaryConstantGatedConeProblemForTest(
          /*coneStateCount=*/4096);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The stable guard is a proven auxiliary constant.  Encode it as a fixed
  // transition leaf so the large gated cone never enters the Craig proof.
  EXPECT_NE(
      stderrOutput.find("imc Craig auxiliary constants=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("transition_states=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_FALSE(result.auxiliaryStateInvariants.empty());
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSubstitutesInactiveResetBeforeTransitionEncoding) {
  const KInductionProblem problem =
      buildCraigResetGatedConeProblemForTest(/*coneStateCount=*/4096);
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Reset is already asserted inactive in every post-bootstrap transition
  // query. Substitute that value before Tseitin encoding so reset-gated support
  // cones do not enter the Craig proof trace.
  EXPECT_NE(
      stderrOutput.find("transition_states=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSkipsAuxiliaryInvariantMiningForHugeBootstrap) {
  const KInductionProblem problem =
      buildCraigAuxiliaryCandidateGuardProblemForTest(10001);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  (void)checker.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Auxiliary constants are only a pruning aid.  Very large bootstrap maps
  // must not spend memory mining transition support before Craig IMC starts.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10001"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig auxiliary constants="),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesLocalAuxiliaryConstantsAfterGlobalMiningSkip) {
  const KInductionProblem problem =
      buildCraigLocalAuxiliaryInvariantProblemForTest(10001);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The global pass must still avoid the huge bootstrap map, but a SAT witness
  // over the tiny projected support can mine transition-proven constants
  // locally and retry the same strict Craig query.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10003"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig local auxiliary constants=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesLocalAuxiliaryEqualitiesAfterGlobalMiningSkip) {
  const KInductionProblem problem =
      buildCraigLocalAuxiliaryEqualityProblemForTest(10000);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The local pass should reuse the same same-design equality proof as the
  // global auxiliary path, but only over the selected Craig support.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10003"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig local auxiliary constants=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("equalities=1 equality_candidates=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcClosesLocalAuxiliaryConstantsOverComplementRails) {
  const KInductionProblem problem =
      buildCraigComplementLocalAuxiliaryProblemForTest(10000);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Same-design complement rails are already part of the Craig state
  // semantics.  Local auxiliary mining should close constants over those pairs
  // so support that only contains QN can still use a proven Q constant.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10003"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=3 candidates=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcPrioritizesSelfPreservingLocalAuxiliaryCandidates) {
  const KInductionProblem problem =
      buildCraigSelfPreservingLocalAuxiliaryPriorityProblemForTest(
          /*noiseStateCount=*/8193,
          /*paddingStateCount=*/2000);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The stable guard is above the first 1024 symbols in a huge support pool.
  // Candidate scoring should promote cheaply self-preserving constants before
  // fan-in-only noise, then the normal transition proof validates the choice.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10195"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=1 candidates=8195 "
          "selected=1025 candidate_limit=1024"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSelectsLocalAuxiliaryConstantsFromLargeSupport) {
  const KInductionProblem problem =
      buildCraigLargeLocalAuxiliaryInvariantProblemForTest(
          /*supportNoiseCount=*/4097,
          /*paddingStateCount=*/6000);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Local auxiliary mining must not give up when the current projection support
  // is above the candidate budget. It should select a bounded slice, prove the
  // stable guard constant, and retry before growing the projection.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10099"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=1 candidates=4099 "
          "selected=4096"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcShrinksLocalAuxiliarySliceForHugeSupport) {
  const KInductionProblem problem =
      buildCraigLargeLocalAuxiliaryInvariantProblemForTest(
          /*supportNoiseCount=*/8193,
          /*paddingStateCount=*/2000);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Huge local support pools are sorted by transition fanin score and then
  // capped before the transition-proof pass.  The useful guard is still in the
  // selected slice, but IMC avoids thousands of low-value SAT checks.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10195"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=1 candidates=8195 "
          "selected=1024 candidate_limit=1024"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcProvesLocalAuxiliaryConstantsFromKnownSupport) {
  const KInductionProblem problem =
      buildCraigKnownSupportLocalAuxiliaryProblemForTest(
          /*supportStateCount=*/300,
          /*paddingStateCount=*/9800);
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The guard's next-state support is larger than the SAT proof support limit,
  // but every support bit is already a local candidate constant.  IMC should
  // prove the whole local set by evaluating that known support instead of
  // expanding the projection.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10102"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=302 candidates=302 "
          "selected=302"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesSecondLocalAuxiliaryRetryBeforeProjection) {
  const KInductionProblem problem =
      buildCraigLocalAuxiliaryRetryLimitProblemForTest(
          /*supportStateCount=*/4097,
          /*paddingStateCount=*/6000);
  const std::unordered_set<size_t> initialTrackedStates = {2};
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.maxProjectionStates = 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The first bounded local-aux slice is useful but leaves two support bits
  // unconstrained.  A second bounded slice should prove the remaining local
  // constants before IMC grows the Craig projection.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips auxiliary invariants candidates=10098"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=4096 candidates=4098 "
          "selected=4096"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig local auxiliary constants=2 candidates=2 "
          "selected=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig local auxiliary retry limit reached"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig refines transition projection states=1->"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesDirectCadicalProfileForLargeProjectionQueries) {
  const KInductionProblem problem =
      buildCraigLargeProjectionAuxiliarySuppressionProblemForTest(
          /*paddingStateCount=*/6200);
  const std::unordered_set<size_t> initialTrackedStates(
      problem.state0Symbols.begin(), problem.state0Symbols.end());
  CraigImcOptions options;
  options.enableAuxiliaryInvariants = true;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // High-projection image queries keep the strict Craig formula, but must avoid
  // CaDiCaL variable elimination because Craig tracing records every derived
  // resolvent. The direct profile keeps that optimization local to IMC.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig uses direct CaDiCaL profile phase=image "
          "tracked_states=6202"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcOptionsUseDirectConcreteCubeSourceWithoutEnv) {
  const KInductionProblem problem = buildCraigAuxiliaryConstantProblemForTest();
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Large dual-rail IMC turns this on through options so the exact bootstrap
  // cube does not need to be synthesized as another Craig source region.
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("direct_cube_source=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcCarriesModerateHelperAuxiliarySeed) {
  constexpr size_t firstState = 2;
  constexpr size_t stateCount = 600;
  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  for (size_t index = 0; index < stateCount; ++index) {
    const size_t symbol = firstState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }
  problem.bad = BoolExpr::Var(firstState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();

  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  for (const size_t symbol : problem.state0Symbols) {
    options.helperAuxiliaryStateInvariants.push_back({symbol, false});
  }
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // BP's retained-helper tail carries roughly 3K transition-proven facts.
  // Keep moderate helper packets available so later singleton outputs do not
  // rebuild the same auxiliary pruning from scratch.
  EXPECT_NE(
      stderrOutput.find("imc Craig seeds helper auxiliary invariants"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig skips broad helper auxiliary invariants"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSkipsBroadHelperAuxiliarySeed) {
  constexpr size_t firstState = 2;
  constexpr size_t stateCount = 5000;
  KInductionProblem problem;
  problem.resetBootstrapCycles = 1;
  for (size_t index = 0; index < stateCount; ++index) {
    const size_t symbol = firstState + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }
  problem.bad = BoolExpr::Var(firstState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.totalStateCount = problem.state0Symbols.size();

  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  for (const size_t symbol : problem.state0Symbols) {
    options.helperAuxiliaryStateInvariants.push_back({symbol, false});
  }
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      /*initialTrackedStates=*/nullptr,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Very broad helper auxiliary packets are not a compact proof hint anymore;
  // they become a side condition on every image query.  Skip them and let
  // strict IMC rediscover only the local facts it actually needs.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips broad helper auxiliary invariants "
          "constants=5000 equalities=0 limit=4096"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig seeds helper auxiliary invariants"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcAdvancesFocusedLookaheadAfterSaturatedQBudget) {
  CraigImcGrowthBudget budget;
  budget.maxInterpolantClauses = 100000;
  budget.maxInterpolantLiterals = 250000;
  budget.maxInterpolantAuxiliaries = 50000;

  // The saturated focused q-pass limit is only a proof-size guard.  Strict IMC
  // can still try the next unroll depth from S0 instead of reporting the
  // singleton inconclusive or attempting the known-explosive next q pass.
  EXPECT_TRUE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/12,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/82359,
      /*interpolantLiterals=*/192171,
      /*interpolantAuxiliaries=*/27453));
  EXPECT_TRUE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/7,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/214323,
      /*interpolantLiterals=*/500087,
      /*interpolantAuxiliaries=*/71441));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/true,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/12,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/82359,
      /*interpolantLiterals=*/192171,
      /*interpolantAuxiliaries=*/27453));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/5,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/13443,
      /*interpolantLiterals=*/31367,
      /*interpolantAuxiliaries=*/4481));
  EXPECT_TRUE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/3,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/2541,
      /*interpolantLiterals=*/5929,
      /*interpolantAuxiliaries=*/847,
      /*qExpansionPassLimit=*/3));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/3,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/2541,
      /*interpolantLiterals=*/5929,
      /*interpolantAuxiliaries=*/847,
      /*qExpansionPassLimit=*/6));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"clauses",
      /*qExpansionPass=*/12,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/82359,
      /*interpolantLiterals=*/192171,
      /*interpolantAuxiliaries=*/27453));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/5,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/82359,
      /*interpolantLiterals=*/192171,
      /*interpolantAuxiliaries=*/27453));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/12,
      /*lookahead=*/32,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/82359,
      /*interpolantLiterals=*/192171,
      /*interpolantAuxiliaries=*/27453));
  EXPECT_TRUE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/12,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/567639,
      /*interpolantLiterals=*/1324491,
      /*interpolantAuxiliaries=*/189213));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterSaturatedFocusedQBudget(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*qExpansionPass=*/12,
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget,
      /*interpolantClauses=*/900000,
      /*interpolantLiterals=*/1900000,
      /*interpolantAuxiliaries=*/260000));
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcAdvancesFocusedLookaheadAfterBudgetedSatFrontier) {
  CraigImcGrowthBudget budget;
  budget.maxQExpansionPass = 4;

  // BP's retained singleton tail reached a fully tracked focused SAT frontier
  // at lookahead 3/q7.  That is strict IMC's cue to try k=4, not a reason to
  // stop merely because the q-pass growth guard had already fired.
  EXPECT_TRUE(shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/true,
      /*budgetReason=*/"q_pass",
      /*lookahead=*/3,
      /*maxLookahead=*/32,
      budget));
  EXPECT_FALSE(shouldAdvanceCraigLookaheadAfterBudgetedFocusedSat(
      /*focusedTransitionProjection=*/true,
      /*hasUntrackedTransitionSupport=*/false,
      /*budgetReason=*/"q_pass",
      /*lookahead=*/32,
      /*maxLookahead=*/32,
      budget));
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcBulkImportsModestFocusedProjectionSupport) {
  // BP's early single-output churn reached a focused 36K transition cone and
  // was paying one saturated Craig rebuild per 4K imported states.  Import that
  // modest focused support in one strict refinement step.
  EXPECT_EQ(
      craigBoundedProjectionRefinementLimit(
          /*candidateCount=*/26232,
          /*transitionSupportSize=*/36088,
          /*focusedTransitionProjection=*/true),
      26232u);
  EXPECT_EQ(
      craigBoundedProjectionRefinementLimit(
          /*candidateCount=*/21562,
          /*transitionSupportSize=*/36088,
          /*focusedTransitionProjection=*/true),
      21562u);

  // Larger focused pools still use the bounded stride that protects the
  // memory-heavy BP tail.
  EXPECT_EQ(
      craigBoundedProjectionRefinementLimit(
          /*candidateCount=*/32769,
          /*transitionSupportSize=*/36088,
          /*focusedTransitionProjection=*/true),
      4096u);
  EXPECT_EQ(
      craigBoundedProjectionRefinementLimit(
          /*candidateCount=*/80046,
          /*transitionSupportSize=*/117866,
          /*focusedTransitionProjection=*/true),
      512u);
  EXPECT_EQ(
      craigBoundedProjectionRefinementLimit(
          /*candidateCount=*/96,
          /*transitionSupportSize=*/36088,
          /*focusedTransitionProjection=*/false),
      64u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcSkipsLocalAuxiliaryMiningForLargeRetainedHelperTail) {
  // Helper-backed singleton outputs still need their small local-auxiliary
  // packets; those transition-proven facts closed BP's first hard singleton
  // before the projection grew into the broad retained-helper tail.
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/false,
      /*trackedStateCount=*/448,
      /*transitionSupportSize=*/736,
      /*helperInvariantRegionCount=*/3));
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/1878,
      /*transitionSupportSize=*/10640,
      /*helperInvariantRegionCount=*/3));

  EXPECT_TRUE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/42584,
      /*transitionSupportSize=*/84516,
      /*helperInvariantRegionCount=*/6));

  // Once helper auxiliaries prune BP's tail to 49K support, the existing
  // bounded local-auxiliary pass is cheaper than another deep lookahead walk.
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/49024,
      /*transitionSupportSize=*/49024,
      /*helperInvariantRegionCount=*/6));
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/false,
      /*trackedStateCount=*/36088,
      /*transitionSupportSize=*/36644,
      /*helperInvariantRegionCount=*/6));
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/36088,
      /*transitionSupportSize=*/36644,
      /*helperInvariantRegionCount=*/3));
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/9856,
      /*transitionSupportSize=*/36088,
      /*helperInvariantRegionCount=*/6));
  EXPECT_FALSE(shouldSkipCraigLocalAuxiliaryMiningForLargeRetainedHelper(
      /*focusedTransitionProjection=*/true,
      /*trackedStateCount=*/36088,
      /*transitionSupportSize=*/9856,
      /*helperInvariantRegionCount=*/0));
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcShortensSaturatedQReplayForLargeRetainedHelperTail) {
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      3u);

  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/49024,
          /*transitionSupportSize=*/49024,
          /*helperInvariantRegionCount=*/6),
      3u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/49682,
          /*transitionSupportSize=*/49682,
          /*helperInvariantRegionCount=*/3),
      6u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/39728,
          /*transitionSupportSize=*/39728,
          /*helperInvariantRegionCount=*/3),
      6u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/3),
      6u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/85498,
          /*transitionSupportSize=*/85498,
          /*helperInvariantRegionCount=*/3),
      6u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/9856,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      6u);
  EXPECT_EQ(
      craigFocusedSaturatedQExpansionPassLimit(
          /*focusedTransitionProjection=*/false,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      6u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcRefinesLargeRetainedHelperTailProjectionAfterQ3) {
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      3u);

  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/49024,
          /*transitionSupportSize=*/49024,
          /*helperInvariantRegionCount=*/6),
      3u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/49682,
          /*transitionSupportSize=*/49682,
          /*helperInvariantRegionCount=*/3),
      0u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/39728,
          /*transitionSupportSize=*/84336,
          /*helperInvariantRegionCount=*/3),
      0u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/3),
      0u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/82450,
          /*transitionSupportSize=*/85498,
          /*helperInvariantRegionCount=*/3),
      0u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/true,
          /*trackedStateCount=*/9856,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      0u);
  EXPECT_EQ(
      craigFocusedProjectionRefinementQExpansionPassLimit(
          /*focusedTransitionProjection=*/false,
          /*trackedStateCount=*/42584,
          /*transitionSupportSize=*/84516,
          /*helperInvariantRegionCount=*/6),
      0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcCapsExplosiveFocusedImageRequestExpansion) {
  // BP's lookahead-4 tail expanded from a useful 9,336-state suffix request to
  // a 37,984-state transition request and 130K local leaves.  Keep that query
  // as a weaker strict over-approximation instead of materializing the spike.
  EXPECT_FALSE(shouldCapCraigFocusedImageTransitionRequests(
      /*expandedRequestCount=*/9336));
  EXPECT_TRUE(shouldCapCraigFocusedImageTransitionRequests(
      /*expandedRequestCount=*/37984));
  EXPECT_EQ(
      cappedCraigFocusedImageTransitionRequestCount(
          /*currentRequestCount=*/9336,
          /*expandedRequestCount=*/37984),
      12000u);
  EXPECT_EQ(
      craigFocusedImageTransitionRequestLimit(
          /*trackedStateCount=*/49024,
          /*helperInvariantRegionCount=*/6),
      8192u);
  EXPECT_EQ(
      craigFocusedImageTransitionRequestLimit(
          /*trackedStateCount=*/36088,
          /*helperInvariantRegionCount=*/6),
      12000u);
  EXPECT_EQ(
      craigFocusedImageTransitionRequestLimit(
          /*trackedStateCount=*/49024,
          /*helperInvariantRegionCount=*/3),
      10000u);
  EXPECT_EQ(
      craigFocusedImageTransitionRequestLimit(
          /*trackedStateCount=*/58104,
          /*helperInvariantRegionCount=*/9),
      8192u);
  EXPECT_EQ(
      craigFocusedImageTransitionRequestLimit(
          /*trackedStateCount=*/9856,
          /*helperInvariantRegionCount=*/6),
      12000u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcCombinesRetainedSingletonWithReusableHelper) {
  EXPECT_TRUE(shouldCombineCraigHelpersForSmallRawSingleton(
      /*useSmallRawSingletonInvariant=*/true,
      /*reusableInvariantHasRegions=*/true));
  EXPECT_FALSE(shouldCombineCraigHelpersForSmallRawSingleton(
      /*useSmallRawSingletonInvariant=*/false,
      /*reusableInvariantHasRegions=*/true));
  EXPECT_FALSE(shouldCombineCraigHelpersForSmallRawSingleton(
      /*useSmallRawSingletonInvariant=*/true,
      /*reusableInvariantHasRegions=*/false));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcProjectionBudgetCoversRetainedHelperTail) {
  // BP's retained-helper tail needs the full 84,516-state focused projection
  // surface, but should stay below the broader >90K memory-risk regime.
  EXPECT_GE(largeDualRailCraigImcProjectionStateLimit(), 84516u);
  EXPECT_LT(largeDualRailCraigImcProjectionStateLimit(), 90000u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcUsesProofDerivedCraigInterpolation) {
  const KInductionProblem problem = buildCraigResetSecProblem(true);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(4);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("imc Craig projection round="), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcStartupProbeUsesRawSupportSizing) {
  const KInductionProblem problem = buildCraigResetSecProblem(true);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Startup witness probes are concrete bounded SAT checks.  Their ordering
  // should not spend a full Craig projection-closure walk before the first
  // proof attempt; RISC-V-sized output sets were dominated by that sizing pass.
  EXPECT_NE(
      stderrOutput.find("imc large dual-rail bounded witness probes="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("projection_sizing=raw_support"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcBatchesWideSharedOutputCones) {
  const KInductionProblem problem = buildWideSharedConeImcProblem(3);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Each output individually exceeds the old 128-symbol IMC Craig support
  // limit.  Since the support is almost identical, one classic IMC query over
  // the three-output conjunction avoids rebuilding the same wide cone per bit.
  EXPECT_NE(
      stderrOutput.find("imc Craig output batch first=0 end=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(result.status, IMCStatus::Different);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcPreSplitsDisjointWideOutputCones) {
  const KInductionProblem problem = buildDisjointWideConeImcBatchProblem();
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});

  // The two outputs are adjacent but share no support. IMC must split them
  // before SAT instead of spending minutes learning that the combined Craig
  // proof is too broad.
  ASSERT_EQ(batches.size(), 2u);
  EXPECT_EQ(batches[0], std::make_pair(size_t{0}, size_t{1}));
  EXPECT_EQ(batches[1], std::make_pair(size_t{1}, size_t{2}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcBatchesOutputsWithSharedProjectionSurface) {
  const KInductionProblem problem = buildProjectionSharedImcBatchProblem();
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});

  // Raw output support alone would split these outputs because the bad
  // predicates touch different output registers. Craig IMC pays for the
  // transition projection, though, and both outputs expand to the same
  // same-design state surface. Batch them so future AES-like cases do not
  // rediscover that projection one bit at a time.
  ASSERT_EQ(batches.size(), 1u);
  EXPECT_EQ(batches[0], std::make_pair(size_t{0}, size_t{2}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcAdaptsBatchSizeForHugeProjectionSurface) {
  const KInductionProblem problem = buildProjectionSharedImcBatchProblem(
      /*sharedTransitionStates=*/2050,
      /*outputCount=*/3);
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});

  // Each output has tiny raw support, but its Craig projection pulls a wide
  // same-design transition surface. Keep the batch small before the first
  // interpolating SAT query creates a wide OR'ed bad predicate.
  ASSERT_EQ(batches.size(), 2u);
  EXPECT_EQ(batches[0], std::make_pair(size_t{0}, size_t{2}));
  EXPECT_EQ(batches[1], std::make_pair(size_t{2}, size_t{3}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcSkipsProjectionCacheForHugeOutputStateProduct) {
  constexpr size_t kSharedTransitionStates = 1000;
  constexpr size_t kOutputCount = 400;
  const KInductionProblem problem = buildProjectionSharedImcBatchProblem(
      kSharedTransitionStates,
      kOutputCount);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // BP-like dual-rail problems have many outputs over a huge state surface.
  // Avoid materializing one transition-closed projection per output before
  // the first Craig query; raw support batching keeps the cache bounded.
  EXPECT_NE(
      stderrOutput.find("imc Craig skips projection support cache"),
      std::string::npos)
      << stderrOutput;
  ASSERT_EQ(batches.size(), kOutputCount);
  EXPECT_EQ(batches.front(), std::make_pair(size_t{0}, size_t{1}));
  EXPECT_EQ(
      batches.back(),
      std::make_pair(kOutputCount - 1, kOutputCount));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcSkipsProjectionCacheForDenseRiscvSurface) {
  constexpr size_t kSharedTransitionStates = 4096;
  constexpr size_t kOutputCount = 65;
  const KInductionProblem problem = buildProjectionSharedImcBatchProblem(
      kSharedTransitionStates,
      kOutputCount);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // sky130*_riscv32i has 4K+ dual-rail state bits and a 99-output bus
  // surface.  That sits below the old product cap but still spends minutes
  // computing transition closures before Craig IMC starts.
  EXPECT_NE(
      stderrOutput.find("imc Craig skips projection support cache"),
      std::string::npos)
      << stderrOutput;
  ASSERT_EQ(batches.size(), kOutputCount);
  EXPECT_EQ(batches.front(), std::make_pair(size_t{0}, size_t{1}));
  EXPECT_EQ(
      batches.back(),
      std::make_pair(kOutputCount - 1, kOutputCount));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcBudgetSplitReportsFirstUnprovenOutput) {
  KInductionProblem problem = buildCraigResetSecProblem(/*equivalent=*/true);
  problem.observedOutputNames = {"out0", "out1"};
  problem.observedOutputExprs0.push_back(problem.observedOutputExprs0.front());
  problem.observedOutputExprs1.push_back(problem.observedOutputExprs1.front());
  // The duplicate output keeps the proof tiny but still exercises the same
  // multi-output Craig batch split and first-unproven-output propagation.
  problem.property = BoolExpr::And(
      makeEqualityExpr(
          problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]),
      makeEqualityExpr(
          problem.observedOutputExprs0[1], problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  const ScopedEnvVar maxQPass("KEPLER_SEC_IMC_CRAIG_MAX_Q_PASS", "1");

  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The scheduled parent may seed the shared sibling surface, but after a
  // budget split each child must restart from only its own output range.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig output batch first=0 end=2 first_name=out0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig splitting output batch after growth budget first=0 end=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig output batch first=0 end=1 first_name=out0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig stopping after inconclusive output batch first=0 end=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, IMCStatus::Inconclusive);
  ASSERT_TRUE(result.firstUnprovenOutput.has_value());
  EXPECT_EQ(*result.firstUnprovenOutput, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesModelGuidedBootstrapProjectionRefinement) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildBootstrapModelGuidedCraigProjectionProblem();
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The first projected SAT model leaves only a small near-saturated remainder,
  // so Craig IMC imports that full remainder instead of spending another
  // bounded refinement pass on the same support.
  EXPECT_NE(
      stderrOutput.find("imc Craig projection round=0 states=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig imports near-saturated projection remainder "
          "selected=96 tracked_states=1 full=97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig model-guided projection refinement "
          "candidates=96 selected=64 full=97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig refines transition projection states=1->97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(result.status, CraigImcStatus::CounterexampleCandidate);
  EXPECT_NE(result.status, CraigImcStatus::BudgetExceeded);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcKeepsModelGuidedRefinementSemanticallyMinimal) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildBootstrapModelGuidedCraigProjectionProblem(
          /*supportCount=*/96,
          /*assignSupportBootstrap=*/true,
          /*addDualRailPartners=*/true);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Missing dual-rail partners make the projected query an over-approximation,
  // so proof soundness is preserved.  The near-saturated import pulls the whole
  // remaining transition surface, including partners, in one pass.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig imports near-saturated projection remainder "
          "selected=192 tracked_states=1 full=193"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig model-guided projection refinement "
          "candidates=96 selected=64 full=193"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig refines transition projection states=1->193"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(result.status, CraigImcStatus::CounterexampleCandidate);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcBoundsLargeProjectionRefinementWithoutModelHint) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildBootstrapModelGuidedCraigProjectionProblem(
          /*supportCount=*/96,
          /*assignSupportBootstrap=*/false);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // When the SAT model does not contradict reset-known support states, Craig
  // still recognizes the small near-saturated remainder and imports it in one
  // pass instead of iterating a bounded 64-state slice.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig imports near-saturated projection remainder "
          "selected=96 tracked_states=1 full=97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement "
          "candidates=96 selected=64 full=97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig refines transition projection states=1->97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(result.status, CraigImcStatus::CounterexampleCandidate);
  EXPECT_NE(result.status, CraigImcStatus::BudgetExceeded);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcScoresBoundedProjectionRefinementByFanin) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest();
  std::unordered_set<size_t> initialTrackedStates = {2, 3, 4, 5, 6};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = 600;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // When model-guided reset facts are unavailable, prioritize support states
  // that feed many tracked next-state functions.  This keeps high-fan-in
  // control symbols ahead of low-id decoy support.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=97"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=5"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcWidensBoundedProjectionRefinementForHugeSupport) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest(
          /*helperTrackedCount=*/4,
          /*decoySupportCount=*/8200);
  std::unordered_set<size_t> initialTrackedStates = {2, 3, 4, 5, 6};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = initialTrackedStates.size() + 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // BP-like support pools should move by a wider scored slice than the small
  // 64-symbol default, while still ranking the high-fan-in control first.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=8201 "
          "selected=512"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=5"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcWidensHighSupportRefinementWhileScoresAreStrong) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest(
          /*helperTrackedCount=*/100,
          /*decoySupportCount=*/84000);
  std::unordered_set<size_t> initialTrackedStates = {2};
  for (size_t symbol = 3; symbol < 103; ++symbol) {
    initialTrackedStates.insert(symbol);
  }
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = initialTrackedStates.size() + 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The partial focused image cap leaves BP's hard tail around 84K support.
  // While the fan-in score is still strong, use the wider strict slice instead
  // of rebuilding the same capped proof one 256-state step at a time.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=83998 "
          "selected=4096"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=101"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcUsesWideStrideForHighSupportLowScoreRefinement) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest(
          /*helperTrackedCount=*/4,
          /*decoySupportCount=*/84000);
  std::unordered_set<size_t> initialTrackedStates = {2, 3, 4, 5, 6};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = initialTrackedStates.size() + 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // BP's partial-cap tail enters a flat-score surface near 84K support.  Move
  // by the focused strict slice there; the >100K test below keeps the tighter
  // cap for truly huge low-score cones.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=84001 "
          "selected=4096"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=5"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcWidensVeryHighSupportRefinementWhileScoresAreStrong) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest(
          /*helperTrackedCount=*/100,
          /*decoySupportCount=*/100000);
  std::unordered_set<size_t> initialTrackedStates = {2};
  for (size_t symbol = 3; symbol < 103; ++symbol) {
    initialTrackedStates.insert(symbol);
  }
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = initialTrackedStates.size() + 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Very large BP tails should still move quickly while the fan-in scorer has a
  // strong signal.  The low-score plateau test below keeps the later cap.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=99998 "
          "selected=1024"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("top_score=101"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcFreezesAndCapsLowScoreFaninRefinementForHugeSupport) {
  const KInductionProblem problem =
      buildHighFaninCraigProjectionProblemForTest(
          /*helperTrackedCount=*/4,
          /*decoySupportCount=*/100000);
  std::unordered_set<size_t> initialTrackedStates = {2, 3, 4, 5, 6};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = initialTrackedStates.size() + 1;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Once a BP-sized support pool has only low fan-in scores left, do not keep
  // materializing enormous resolver support sets just to refresh the ranking.
  const std::string freezeDiag =
      "imc Craig freezes low-score fanin scoring "
      "candidates=100001 top_score=5 support=100006";
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig skips local auxiliary invariants support=100006 "
          "support_limit=65536"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig weakens projected local state semantics leaves=100006 "
          "leaf_limit=100000"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(freezeDiag), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      KEPLER_FORMAL::SEC::detail::countTextOccurrences(
          stderrOutput, freezeDiag),
      1u)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig caps low-score bounded refinement "
          "candidates=100001 selected_limit=128 top_score=5 "
          "support=100006 score_limit=64"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig bounded projection refinement candidates=100001 "
          "selected=128 full=100006 top_score=5"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcBackfillsTinyModelGuidedRefinementForLargeSupport) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildTinyModelGuidedBackfillCraigProjectionProblemForTest(
          /*decoyPairCount=*/4100);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // A large support cone with a one-bit SAT witness should not spend one Craig
  // rebuild per model-guided bit.  Backfill the tiny witness with a bounded
  // scored slice in the same strict IMC refinement step.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig model-guided projection refinement candidates=1 "
          "selected=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig bounded projection refinement candidates="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("selected=512"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcBackfillsMediumModelGuidedRefinementForLargeSupport) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildTinyModelGuidedBackfillCraigProjectionProblemForTest(
          /*decoyPairCount=*/4100,
          /*controlCount=*/128);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  (void)checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // On large supports the model-guided slice can be much larger than 64 and
  // still smaller than the adaptive 512-symbol stride.  Backfill those medium
  // witnesses too, so BP-sized cones do not spend one rebuild per model slice.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig model-guided projection refinement candidates=128 "
          "selected=128"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig bounded projection refinement candidates="),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcStopsProjectionRefinementAtGrowthBudget) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildBootstrapModelGuidedCraigProjectionProblem(
          /*supportCount=*/96,
          /*assignSupportBootstrap=*/false);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxProjectionStates = 32;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // A single hard output cannot be bisected further. The transition build can
  // complete, but the projection-state guard still reports the strict growth
  // budget before importing an oversized refinement slice.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig growth budget exceeded reason=projection_states"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(" state_limit=32"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::BudgetExceeded);
}

TEST_F(SequentialEquivalenceStrategyTests,
       CraigImcStopsOversizedImageBeforeCraigBuild) {
  constexpr size_t outputState = 2;
  const KInductionProblem problem =
      buildBootstrapModelGuidedCraigProjectionProblem(
          /*supportCount=*/96,
          /*assignSupportBootstrap=*/false);
  const std::unordered_set<size_t> initialTrackedStates = {outputState};
  CraigImcOptions options;
  options.enableDirectConcreteCubeSource = true;
  options.growthBudget.enabled = true;
  options.growthBudget.maxImageTransitionStates = 32;
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  CraigInterpolatingModelChecker checker(
      problem,
      /*helperInvariantRegions=*/nullptr,
      &initialTrackedStates,
      options);
  const CraigImcResult result = checker.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The transition planner can see the projected image surface before the
  // Craig solver is allocated.  Trip the strict growth budget there instead of
  // materializing a proof trace that the caller will discard as inconclusive.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig growth budget exceeded reason=transition_build"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(" transition_support=97"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(" support_limit=32"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("imc Craig image build begin"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, CraigImcStatus::BudgetExceeded);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcBudgetRetrySingleOutputWithSelfProjection) {
  const KInductionProblem problem = buildProjectionSharedImcBatchProblem(
      /*sharedTransitionStates=*/1100,
      /*outputCount=*/1);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  const ScopedEnvVar maxQPass("KEPLER_SEC_IMC_CRAIG_MAX_Q_PASS", "1");

  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES-like single-output slices can have tiny bad support but a huge
  // precomputed transition seed. After the first seeded Craig attempt exceeds
  // budget, retry strict IMC from the property's own support instead of
  // spending the rest of the run on the same oversized seed.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig output batch first=0 end=1 "
          "first_name=projection_shared_out0 bad_support=1 "
          "tracked_seed_states=1102"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig retrying single output with self projection first=0 "
          "end=1 dropped_seed_states=1102"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig projection round=0 states=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(result.status, IMCStatus::Different);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcPreSplitsHugeSharedOutputCones) {
  const KInductionProblem problem =
      buildWideSharedConeImcProblem(/*outputCount=*/2,
                                    /*sharedSupportCount=*/1100);
  const auto batches = buildLargeDualRailCraigImcOutputBatches(
      problem, OutputBatchingLimits{/*maxOutputBatchSize=*/8,
                                    /*outputBatchSupportLimit=*/8192});

  // Huge shared cones look batch-friendly by overlap alone, but the combined
  // bad predicate creates much larger Craig interpolants. Keep them single.
  ASSERT_EQ(batches.size(), 2u);
  EXPECT_EQ(batches[0], std::make_pair(size_t{0}, size_t{1}));
  EXPECT_EQ(batches[1], std::make_pair(size_t{1}, size_t{2}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcReusesCraigInvariantAcrossWideBatches) {
  const KInductionProblem problem = buildWideSharedConeImcProblem(9);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The first eight outputs reach the McMillan fixed point with one Craig
  // step. Later outputs have the same reachable-state surface, so IMC can
  // reuse the saved region against the new bad predicate.
  EXPECT_NE(
      stderrOutput.find("imc Craig output batch first=0 end=8"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig reused invariant for output batch first=8 end=9"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, IMCStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcReusesSmallCraigInvariantAcrossBatches) {
  const KInductionProblem problem =
      buildWideSharedConeImcProblem(/*outputCount=*/9,
                                    /*sharedSupportCount=*/8);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES proves a small control output before the text bus. Keep such Craig
  // invariants available as helpers even when they track far fewer than the
  // wide datapath-state threshold.
  EXPECT_NE(
      stderrOutput.find("imc Craig output batch first=0 end=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig reused invariant for output batch first=1 end=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig reused invariant for output batch first=8 end=9"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, IMCStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LargeDualRailImcSeedsAuxiliaryInvariantsFromReusableHelper) {
  const KInductionProblem problem =
      buildCraigAuxiliaryReuseAcrossOutputsProblemForTest(10001);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const IMCResult result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The saved helper carries transition-proven auxiliary constants into the
  // next strict Craig batch.  This keeps large dual-rail output sweeps from
  // paying the same local-mining cost and proof width repeatedly.
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig output batch first=0 end=1 first_name=aux_reuse_out0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "imc Craig output batch first=1 end=2 first_name=aux_reuse_out1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("imc Craig seeds helper auxiliary invariants"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(result.status, IMCStatus::Equivalent);
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
       IMCEngineFindsDelayedBadStateAfterFailedInductionStep) {
  KInductionProblem problem;
  constexpr size_t firstState = 2;
  constexpr size_t delayedBadState = 3;
  problem.state0Symbols = {firstState, delayedBadState};
  problem.allSymbols = {firstState, delayedBadState};
  problem.transitions0 = {
      {firstState, BoolExpr::createTrue()},
      {delayedBadState, BoolExpr::Var(firstState)}};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(firstState)),
      BoolExpr::Not(BoolExpr::Var(delayedBadState)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.bad = BoolExpr::Var(delayedBadState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  // IMC now tries the sound induction step before frontier witness search.  If
  // that proof does not close, it must still extend the concrete frontier and
  // report a real delayed counterexample.
  EXPECT_EQ(result.status, IMCStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.bound, 2u);
  EXPECT_EQ(result.witness->badFrame, 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineProofOnlyFrontierStillLocalizesReachableBad) {
  KInductionProblem problem;
  constexpr size_t firstState = 2;
  constexpr size_t delayedBadState = 3;
  problem.observedOutputNames = {"stable", "delayed"};
  problem.state0Symbols = {firstState, delayedBadState};
  problem.allSymbols = {firstState, delayedBadState};
  problem.observedOutputExprs0 = {
      BoolExpr::createFalse(), BoolExpr::Var(delayedBadState)};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.transitions0 = {
      {firstState, BoolExpr::createTrue()},
      {delayedBadState, BoolExpr::Var(firstState)}};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(firstState)),
      BoolExpr::Not(BoolExpr::Var(delayedBadState)));
  problem.initializedStateCount = 2;
  problem.totalStateCount = 13;
  problem.usesDualRailStateEncoding = true;
  problem.property = BoolExpr::And(
      makeEqualityExpr(
          problem.observedOutputExprs0[0], problem.observedOutputExprs1[0]),
      makeEqualityExpr(
          problem.observedOutputExprs0[1], problem.observedOutputExprs1[1]));
  problem.bad = BoolExpr::Not(problem.property);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  // Once a delayed bad state is reachable, IMC must follow the exact witness
  // path and report the mismatch instead of closing from a proof-only summary.
  EXPECT_EQ(result.status, IMCStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.bound, 2u);
  EXPECT_EQ(result.witness->badFrame, 2u);
  ASSERT_EQ(result.witness->outputMismatches.size(), 1u);
  EXPECT_EQ(result.witness->outputMismatches[0].signal, "delayed");
}

TEST_F(SequentialEquivalenceStrategyTests,
       IMCEngineValidatesConcreteBaseWhenInductionCloses) {
  KInductionProblem problem;
  constexpr size_t badState = 2;
  problem.state0Symbols = {badState};
  problem.allSymbols = {badState};
  problem.transitions0 = {{badState, BoolExpr::createTrue()}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(badState));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(badState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  // A closed induction step is only a proof once the concrete SEC base horizon
  // is also clean.  This catches accidental "step-only" equivalence results.
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
       IMCEngineBatchesMultiOutputProblemsWithoutFallbackInduction) {
  auto problem = buildLinearChainSecProblem(4);
  problem.observedOutputNames = {"terminal_state", "low_state_bit"};
  problem.observedOutputExprs0.push_back(BoolExpr::Var(problem.state0Symbols.front()));
  problem.observedOutputExprs1.push_back(BoolExpr::Var(problem.state1Symbols.front()));

  BoolExpr* property = BoolExpr::createTrue();
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    property = BoolExpr::And(
        property,
        makeEqualityExpr(problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]));
  }
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  IMCEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, IMCStatus::Equivalent);
  // Exact IMC should close on the real reachable frontier here. The old
  // fallback-induction shortcut closed one frame earlier, but IMC must not jump
  // into a different proof engine.
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
  problem.usesStrictDualRailEqualityProperty = true;
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
       SequentialEquivalenceResultReportsZeroCoverageWhenNoOutputsExist) {
  SequentialEquivalenceResult result;

  EXPECT_EQ(result.outputCoveragePercent(), 0.0);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsRequiresEngineProofWithoutCrossDesignStateCore) {
  const SignalKey out = makeSignalKey("out");
  const SignalKey stateA0 = makeSignalKey("stateA0");
  const SignalKey stateB0 = makeSignalKey("stateB0");
  const SignalKey stateA1 = makeSignalKey("stateA1");
  const SignalKey stateB1 = makeSignalKey("stateB1");

  SequentialDesignModel model0;
  model0.stateBits = {stateA0, stateB0};
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model0.inputVarByKey.emplace(stateA0, 2);
  model0.inputVarByKey.emplace(stateB0, 3);
  model0.displayNameByKey.emplace(out, "out[0]");
  model0.displayNameByKey.emplace(stateA0, "state_a[0]");
  model0.displayNameByKey.emplace(stateB0, "state_b[0]");
  model0.initialStateValueByKey.emplace(stateA0, false);
  model0.initialStateValueByKey.emplace(stateB0, false);
  model0.nextStateExprByStateKey.emplace(stateA0, BoolExpr::Var(2));
  model0.nextStateExprByStateKey.emplace(stateB0, BoolExpr::Var(3));
  model0.observedOutputExprByKey.emplace(
      out,
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));

  SequentialDesignModel model1;
  model1.stateBits = {stateA1, stateB1};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model1.inputVarByKey.emplace(stateA1, 4);
  model1.inputVarByKey.emplace(stateB1, 5);
  model1.displayNameByKey.emplace(out, "out[0]");
  model1.displayNameByKey.emplace(stateA1, "state_a[0]");
  model1.displayNameByKey.emplace(stateB1, "state_b[0]");
  model1.initialStateValueByKey.emplace(stateA1, false);
  model1.initialStateValueByKey.emplace(stateB1, false);
  model1.nextStateExprByStateKey.emplace(stateA1, BoolExpr::Var(4));
  model1.nextStateExprByStateKey.emplace(stateB1, BoolExpr::Var(5));
  model1.observedOutputExprByKey.emplace(
      out,
      BoolExpr::Not(BoolExpr::Or(
          BoolExpr::Not(BoolExpr::Var(4)),
          BoolExpr::Not(BoolExpr::Var(5)))));

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStdout();
  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::KInduction);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stdoutOutput = testing::internal::GetCapturedStdout();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
  // The output is proved by strict k-induction from per-design init and
  // transition facts; no strategy summary may advertise side-proof coverage.
  EXPECT_EQ(stdoutOutput.find("sat_implied_outputs"), std::string::npos);
  EXPECT_EQ(stdoutOutput.find("abstract_equiv_outputs"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsDualRailDoesNotAdvertisePreEngineCertificates) {
  const SignalKey out = makeSignalKey("pdrDualRailNoPrecheckOut");

  SequentialDesignModel model0;
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model0.displayNameByKey.emplace(out, "no_precheck_out[0]");
  model0.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());

  SequentialDesignModel model1;
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model1.displayNameByKey.emplace(out, "no_precheck_out[0]");
  model1.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stdoutOutput = testing::internal::GetCapturedStdout();
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  // These fields belonged to deleted pre-engine proof paths.  The selected PDR
  // engine may still prove the trivial property, but the strategy must not mark
  // coverage through side implication or flush certificates.
  EXPECT_EQ(stdoutOutput.find("sat_implied_outputs"), std::string::npos);
  EXPECT_EQ(stdoutOutput.find("flush_certified_outputs"), std::string::npos);
  EXPECT_EQ(stdoutOutput.find("implication_conflict_limit"), std::string::npos);
  EXPECT_EQ(stdoutOutput.find("flush_conflict_limit"), std::string::npos);
  EXPECT_EQ(stderrOutput.find("dual-rail flush certificate"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailValidatesRiscvSizedFrameZeroSurface) {
  constexpr size_t kOutputCount = 99;
  SequentialDesignModel model0;
  SequentialDesignModel model1;

  for (size_t index = 0; index < kOutputCount; ++index) {
    const SignalKey out =
        makeSignalKey("pdrDualRailRiscvSizedFrameZeroOut" +
                      std::to_string(index));
    const SignalKey state =
        makeSignalKey("pdrDualRailRiscvSizedFrameZeroState" +
                      std::to_string(index));
    const size_t stateVar = 2000 + index;
    const std::string outputName =
        "riscv_sized_frame_zero_out[" + std::to_string(index) + "]";

    model0.allObservedOutputs.push_back(out);
    model0.observedOutputs.push_back(out);
    model0.displayNameByKey.emplace(out, outputName);
    addStateBitForTest(
        model0,
        state,
        stateVar,
        "riscv_sized_frame_zero_state[" + std::to_string(index) + "]",
        BoolExpr::createFalse());
    model0.initialStateValueByKey.emplace(state, false);
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(stateVar));

    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.displayNameByKey.emplace(out, outputName);
    model1.observedOutputExprByKey.emplace(out, BoolExpr::createFalse());
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // RISC-V exposes 99 observed outputs.  Keep that medium surface inside the
  // exact frame-0 guard so PDR does not rely on abstract empty-cube witnesses
  // when reset/bootstrap facts are enough to prove the post-reset outputs.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, kOutputCount);
  EXPECT_EQ(result.totalOutputs, kOutputCount);
  EXPECT_EQ(
      stderrOutput.find("skipped dual-rail frame-0 validation"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailValidatesDynamicNodeSizedFrameZeroSurface) {
  constexpr size_t kOutputCount = 331;
  SequentialDesignModel model0;
  SequentialDesignModel model1;

  for (size_t index = 0; index < kOutputCount; ++index) {
    const SignalKey out =
        makeSignalKey("pdrDualRailDynamicNodeFrameZeroOut" +
                      std::to_string(index));
    const SignalKey state =
        makeSignalKey("pdrDualRailDynamicNodeFrameZeroState" +
                      std::to_string(index));
    const size_t stateVar = 5000 + index;
    const std::string outputName =
        "dynamic_node_frame_zero_out[" + std::to_string(index) + "]";

    model0.allObservedOutputs.push_back(out);
    model0.observedOutputs.push_back(out);
    model0.displayNameByKey.emplace(out, outputName);
    addStateBitForTest(
        model0,
        state,
        stateVar,
        "dynamic_node_frame_zero_state[" + std::to_string(index) + "]",
        BoolExpr::createFalse());
    model0.initialStateValueByKey.emplace(state, false);
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(stateVar));

    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.displayNameByKey.emplace(out, outputName);
    model1.observedOutputExprByKey.emplace(out, BoolExpr::createFalse());
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Ternary obligations remove unrelated state from each exact output slice,
  // so PDR can prove this initialized 331-output surface without reduction.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, kOutputCount);
  EXPECT_EQ(result.totalOutputs, kOutputCount);
  EXPECT_EQ(
      stderrOutput.find("skipped dual-rail frame-0 validation"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailProvesVeryWideEquivalentSurface) {
  constexpr size_t kOutputCount = 385;
  SequentialDesignModel model0;
  SequentialDesignModel model1;

  for (size_t index = 0; index < kOutputCount; ++index) {
    const SignalKey out =
        makeSignalKey("pdrDualRailVeryWideEquivalentOut" +
                      std::to_string(index));
    const SignalKey state =
        makeSignalKey("pdrDualRailVeryWideEquivalentState" +
                      std::to_string(index));
    const size_t stateVar = 9000 + index;
    const std::string outputName =
        "very_wide_equivalent_out[" + std::to_string(index) + "]";

    model0.allObservedOutputs.push_back(out);
    model0.observedOutputs.push_back(out);
    model0.displayNameByKey.emplace(out, outputName);
    addStateBitForTest(
        model0,
        state,
        stateVar,
        "very_wide_equivalent_state[" + std::to_string(index) + "]",
        BoolExpr::createFalse());
    model0.initialStateValueByKey.emplace(state, false);
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(stateVar));

    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.displayNameByKey.emplace(out, outputName);
    model1.observedOutputExprByKey.emplace(out, BoolExpr::createFalse());
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Ternary obligations let exact PDR prove every initialized output directly;
  // no validation/deferral layer or reduced transition system is involved.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, kOutputCount);
  EXPECT_EQ(result.totalOutputs, kOutputCount);
  EXPECT_EQ(
      stderrOutput.find("deferred wide dual-rail equivalent validation outputs=385"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsKiDualRailFindsProductionResidualMismatch) {
  constexpr size_t kResidualOutputs = 129;
  constexpr size_t kStateBitsPerDesign = 2049;
  auto testCase = makeLargeDualRailResidualCaseForTest(
      "kiDualRailProductionLargeResidual",
      kResidualOutputs,
      kStateBitsPerDesign);
  const SignalKey firstState0 =
      makeSignalKey("kiDualRailProductionLargeResidualState0_0");
  testCase.model0.initialStateValueByKey[firstState0] = true;

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(testCase.model0, testCase.model1, 1);

  // This exceeds the historical large-residual skip threshold.  The selected
  // KI path must still receive the real obligation and report a top-output
  // mismatch instead of pre-skipping the surface.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.coveredOutputs, kResidualOutputs + 1);
  EXPECT_EQ(result.totalOutputs, kResidualOutputs + 1);
  EXPECT_TRUE(result.skippedObservedOutputs.empty());
  EXPECT_NE(result.reason.find("large_residual_out[0]"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsImcDualRailRunsSelectedEngineDirectly) {
  constexpr size_t kResidualOutputs = 1;
  constexpr size_t kStateBitsPerDesign = 2049;
  const auto testCase = makeLargeDualRailResidualCaseForTest(
      "imcDualRailDirectEngine",
      kResidualOutputs,
      kStateBitsPerDesign,
      /*includeImpliedOutput=*/false);

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(testCase.model0, testCase.model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // IMC must enter its own interpolation engine directly. The held X versus
  // constant 0 is outside the steady-state mismatch predicate.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_TRUE(result.skippedObservedOutputs.empty());
  EXPECT_NE(
      stderrOutput.find("SEC diag: entering imc engine"),
      std::string::npos);
  EXPECT_EQ(stderrOutput.find("dual-rail imc"), std::string::npos);
  EXPECT_EQ(stderrOutput.find("k-induction problem"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsKiDualRailRemainsInconclusiveWhenNoOutputIsCovered) {
  constexpr size_t kStatePairs = 5;
  constexpr size_t kOutputs = kStatePairs * 2;
  SequentialDesignModel model0;
  SequentialDesignModel model1;

  for (size_t i = 0; i < kStatePairs; ++i) {
    const SignalKey out = makeSignalKey(
        "pdrDualRailZeroCoveredOutQ" + std::to_string(i));
    const SignalKey outN = makeSignalKey(
        "pdrDualRailZeroCoveredOutQN" + std::to_string(i));
    const SignalKey state0 = makeSignalKey(
        "pdrDualRailZeroCoveredState0" + std::to_string(i));
    const SignalKey state1 = makeSignalKey(
        "pdrDualRailZeroCoveredState1" + std::to_string(i));
    const size_t state0Var = 200 + i;
    const size_t state1Var = 300 + i;
    const std::string outputName =
        "zero_covered_out[" + std::to_string(i) + "]";
    const std::string outputNName =
        "zero_covered_out_n[" + std::to_string(i) + "]";

    model0.allObservedOutputs.push_back(out);
    model0.allObservedOutputs.push_back(outN);
    model0.observedOutputs.push_back(out);
    model0.observedOutputs.push_back(outN);
    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.allObservedOutputs.push_back(outN);
    model1.observedOutputs.push_back(outN);
    model0.displayNameByKey.emplace(out, outputName);
    model0.displayNameByKey.emplace(outN, outputNName);
    model1.displayNameByKey.emplace(out, outputName);
    model1.displayNameByKey.emplace(outN, outputNName);
    addStateBitForTest(
        model0,
        state0,
        state0Var,
        "left_zero_covered_q[" + std::to_string(i) + "]",
        BoolExpr::Var(state0Var));
    addStateBitForTest(
        model1,
        state1,
        state1Var,
        "right_zero_covered_q[" + std::to_string(i) + "]",
        BoolExpr::Not(BoolExpr::Var(state1Var)));
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(state0Var));
    model0.observedOutputExprByKey.emplace(
        outN, BoolExpr::Not(BoolExpr::Var(state0Var)));
    model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(state1Var));
    model1.observedOutputExprByKey.emplace(
        outN, BoolExpr::Not(BoolExpr::Var(state1Var)));
  }

  const ScopedEnvVar batchLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
  const ScopedEnvVar leafLimit(
      "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
  SequentialEquivalenceStrategy kiStrategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto kiResult = kiStrategy.runExtractedModels(model0, model1, 1);

  // A dual-rail residual engine that proves no top output must report zero
  // coverage, not reuse the original all-output coverage surface.
  EXPECT_EQ(kiResult.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(kiResult.coveredOutputs, 0u);
  EXPECT_EQ(kiResult.totalOutputs, kOutputs);
  ASSERT_EQ(kiResult.skippedObservedOutputs.size(), kOutputs);
  EXPECT_NE(
      kiResult.reason.find("Dual-rail k-induction did not prove any output"),
      std::string::npos);
  EXPECT_NE(
      kiResult.skippedObservedOutputs.front().find(
          "k-induction proof was inconclusive"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsSkipsResetUnanchoredStateDependentOutputs) {
  const SignalKey rst = makeSignalKey("skipResetUnanchoredRst");
  const SignalKey data = makeSignalKey("skipResetUnanchoredData");
  const SignalKey good = makeSignalKey("skipResetUnanchoredGood");
  const SignalKey bad = makeSignalKey("skipResetUnanchoredBad");
  const SignalKey state0 = makeSignalKey("skipResetUnanchoredState0");
  const SignalKey state1 = makeSignalKey("skipResetUnanchoredState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst, data};
  model0.stateBits = {state0};
  model0.allObservedOutputs = {good, bad};
  model0.observedOutputs = {good, bad};
  model0.inputVarByKey.emplace(rst, 2);
  model0.inputVarByKey.emplace(data, 4);
  model0.inputVarByKey.emplace(state0, 6);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(data, "data[0]");
  model0.displayNameByKey.emplace(good, "good[0]");
  model0.displayNameByKey.emplace(bad, "bad[0]");
  model0.displayNameByKey.emplace(state0, "u_left.q[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));
  model0.observedOutputExprByKey.emplace(good, BoolExpr::Var(4));
  model0.observedOutputExprByKey.emplace(bad, BoolExpr::Var(6));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst, data};
  model1.stateBits = {state1};
  model1.allObservedOutputs = {good, bad};
  model1.observedOutputs = {good, bad};
  model1.inputVarByKey.emplace(rst, 3);
  model1.inputVarByKey.emplace(data, 5);
  model1.inputVarByKey.emplace(state1, 7);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(data, "data[0]");
  model1.displayNameByKey.emplace(good, "good[0]");
  model1.displayNameByKey.emplace(bad, "bad[0]");
  model1.displayNameByKey.emplace(state1, "u_right.q[0]");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::createFalse());
  model1.observedOutputExprByKey.emplace(good, BoolExpr::Var(5));
  model1.observedOutputExprByKey.emplace(bad, BoolExpr::Var(7));

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::KInduction);
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  // The bad output is top-visible and both state bits have concrete reset
  // values, but its temporal equality would still require assuming an internal
  // flop correspondence.  SEC should report it as uncovered instead.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(
      result.skippedObservedOutputs.front().find("reset-unanchored"),
      std::string::npos);
  ASSERT_EQ(result.resetUnanchoredSkippedOutputs.size(), 1u);
  EXPECT_EQ(
      result.resetUnanchoredSkippedOutputs.front(),
      result.skippedObservedOutputs.front());
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsBinaryReportsWideResetUnanchoredSurfaceAsSkippedWithoutRecovery) {
  const SignalKey rst = makeSignalKey("binaryImcWideResetUnanchoredRst");
  const SignalKey data = makeSignalKey("binaryImcWideResetUnanchoredData");
  const SignalKey state0 = makeSignalKey("binaryImcWideResetUnanchoredState0");
  const SignalKey state1 = makeSignalKey("binaryImcWideResetUnanchoredState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst, data};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(rst, 2);
  model0.inputVarByKey.emplace(data, 4);
  model0.inputVarByKey.emplace(state0, 6);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(data, "data[0]");
  model0.displayNameByKey.emplace(state0, "u_left.q[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst, data};
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(rst, 3);
  model1.inputVarByKey.emplace(data, 5);
  model1.inputVarByKey.emplace(state1, 7);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(data, "data[0]");
  model1.displayNameByKey.emplace(state1, "u_right.q[0]");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::createFalse());

  constexpr size_t kOutputCount = 64;
  for (size_t i = 0; i < kOutputCount; ++i) {
    const SignalKey out =
        makeSignalKey("binaryImcWideResetUnanchoredOut" + std::to_string(i));
    const std::string outputName =
        "binary_imc_reset_unanchored_out[" + std::to_string(i) + "]";
    model0.allObservedOutputs.push_back(out);
    model0.observedOutputs.push_back(out);
    model0.displayNameByKey.emplace(out, outputName);
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(6));

    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.displayNameByKey.emplace(out, outputName);
    model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(7));
  }

  const std::vector<SecEngine> engines = {
      SecEngine::Pdr, SecEngine::KInduction, SecEngine::Imc};
  for (const SecEngine engine : engines) {
    auto strategy = makeBinaryExtractedSecStrategy(engine);
    const auto result = strategy.runExtractedModels(model0, model1, 1);

    // Binary SEC must not silently switch a fully reset-unanchored surface into
    // another engine or encoding, and zero covered outputs must not be reported
    // as a successful equivalence proof.
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
    EXPECT_EQ(result.coveredOutputs, 0u);
    EXPECT_EQ(result.totalOutputs, kOutputCount);
    ASSERT_EQ(result.skippedObservedOutputs.size(), kOutputCount);
    EXPECT_EQ(result.resetUnanchoredSkippedOutputs.size(), kOutputCount);
    EXPECT_NE(
        result.skippedObservedOutputs.front().find("reset-unanchored"),
        std::string::npos);
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsDualRailSteadyIgnoresXContainingCycles) {
  const auto models = makeHeldRailModelsForTest(
      "dualRailXVersusBinary", std::nullopt, false);

  for (const SecEngine engine :
       {SecEngine::Pdr, SecEngine::KInduction, SecEngine::Imc}) {
    SCOPED_TRACE(static_cast<int>(engine));
    SequentialEquivalenceStrategy strategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        engine,
        SecEncoding::DualRailSteady);
    const auto result =
        strategy.runExtractedModels(models.model0, models.model1, 2);

    // The steady-state property checks only cycles where both outputs are
    // binary-defined. An X-versus-binary cycle is outside that property.
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_EQ(result.coveredOutputs, 1u);
    EXPECT_EQ(result.totalOutputs, 1u);
    EXPECT_TRUE(result.skippedObservedOutputs.empty());
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsDualRailLeavesResidualsUncoveredWithoutPdrFallback) {
  const SignalKey good = makeSignalKey("dualRailResidualGood");
  const SignalKey residual0 = makeSignalKey("dualRailResidualState");
  const SignalKey residual1 = makeSignalKey("dualRailResidualStateN");
  const SignalKey data = makeSignalKey("dualRailResidualData");
  const SignalKey state0 = makeSignalKey("dualRailResidualState0");
  const SignalKey state1 = makeSignalKey("dualRailResidualState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {data};
  model0.stateBits = {state0};
  model0.allObservedOutputs = {good, residual0, residual1};
  model0.observedOutputs = {good, residual0, residual1};
  model0.inputVarByKey.emplace(data, 4);
  model0.inputVarByKey.emplace(state0, 2);
  model0.displayNameByKey.emplace(data, "data[0]");
  model0.displayNameByKey.emplace(good, "proved_const[0]");
  model0.displayNameByKey.emplace(residual0, "residual_state[0]");
  model0.displayNameByKey.emplace(residual1, "residual_state_n[0]");
  model0.displayNameByKey.emplace(state0, "left_state_q[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Var(2));
  model0.observedOutputExprByKey.emplace(good, BoolExpr::createTrue());
  model0.observedOutputExprByKey.emplace(residual0, BoolExpr::Var(2));
  model0.observedOutputExprByKey.emplace(residual1, BoolExpr::Not(BoolExpr::Var(2)));

  SequentialDesignModel model1;
  model1.environmentInputs = {data};
  model1.stateBits = {state1};
  model1.allObservedOutputs = {good, residual0, residual1};
  model1.observedOutputs = {good, residual0, residual1};
  model1.inputVarByKey.emplace(data, 5);
  model1.inputVarByKey.emplace(state1, 3);
  model1.displayNameByKey.emplace(data, "data[0]");
  model1.displayNameByKey.emplace(good, "proved_const[0]");
  model1.displayNameByKey.emplace(residual0, "residual_state[0]");
  model1.displayNameByKey.emplace(residual1, "residual_state_n[0]");
  model1.displayNameByKey.emplace(state1, "right_state_q[0]");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));
  model1.observedOutputExprByKey.emplace(good, BoolExpr::createTrue());
  model1.observedOutputExprByKey.emplace(residual0, BoolExpr::Var(3));
  model1.observedOutputExprByKey.emplace(residual1, BoolExpr::Not(BoolExpr::Var(3)));

  {
    const ScopedEnvVar batchLimit(
        "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
    const ScopedEnvVar leafLimit(
        "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
    SequentialEquivalenceStrategy strategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        SecEngine::KInduction,
        SecEncoding::DualRailSteady);
    const auto result = strategy.runExtractedModels(model0, model1, 1);

    // KI may report partial output coverage only for obligations it proved
    // itself; it must not invoke PDR behind the selected engine.
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
    EXPECT_EQ(result.coveredOutputs, 1u);
    EXPECT_EQ(result.totalOutputs, 3u);
    ASSERT_EQ(result.skippedObservedOutputs.size(), 2u);
    EXPECT_NE(
        result.skippedObservedOutputs[0].find("k-induction proof was inconclusive"),
        std::string::npos);
    EXPECT_NE(
        result.skippedObservedOutputs[1].find("k-induction proof was inconclusive"),
        std::string::npos);

    SequentialEquivalenceStrategy pdrStrategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        SecEngine::Pdr,
        SecEncoding::DualRailSteady);
    const auto pdrResult = pdrStrategy.runExtractedModels(model0, model1, 1);

    // PDR proves the steady-state property for all three outputs. Permanent X
    // values are outside the guarded binary-mismatch predicate.
    EXPECT_EQ(pdrResult.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_EQ(pdrResult.coveredOutputs, 3u);
    EXPECT_EQ(pdrResult.totalOutputs, 3u);
    EXPECT_TRUE(pdrResult.skippedObservedOutputs.empty());

    SequentialEquivalenceStrategy imcStrategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        SecEngine::Imc,
        SecEncoding::DualRailSteady);
    const auto imcResult = imcStrategy.runExtractedModels(model0, model1, 1);

    EXPECT_EQ(imcResult.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_EQ(imcResult.coveredOutputs, 3u);
    EXPECT_EQ(imcResult.totalOutputs, 3u);
    EXPECT_TRUE(imcResult.skippedObservedOutputs.empty());
  }

  {
    const ScopedEnvVar batchLimit(
        "KEPLER_SEC_KI_DUAL_RAIL_BATCH_DECISION_LIMIT", "0");
    const ScopedEnvVar leafLimit(
        "KEPLER_SEC_KI_DUAL_RAIL_LEAF_DECISION_LIMIT", "0");
    const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    SequentialEquivalenceStrategy strategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        SecEngine::KInduction,
        SecEncoding::DualRailSteady);
    const auto result = strategy.runExtractedModels(model0, model1, 1);
    const std::string stdoutOutput = testing::internal::GetCapturedStdout();
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    // KI must honor the selected engine: a resource-limited KI proof keeps only
    // the dual-rail implied coverage and leaves residuals uncovered rather than
    // launching a hidden PDR retry.
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
    EXPECT_EQ(result.coveredOutputs, 1u);
    EXPECT_EQ(result.totalOutputs, 3u);
    ASSERT_EQ(result.skippedObservedOutputs.size(), 2u);
    EXPECT_NE(
        result.skippedObservedOutputs[0].find("k-induction proof was inconclusive"),
        std::string::npos);
    EXPECT_EQ(stdoutOutput.find("PDR repair"), std::string::npos);
    EXPECT_EQ(stderrOutput.find("PDR repair"), std::string::npos);
  }

  {
    const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    SequentialEquivalenceStrategy imcStrategy(
        nullptr,
        nullptr,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        SecEngine::Imc,
        SecEncoding::DualRailSteady);
    const auto imcResult = imcStrategy.runExtractedModels(model0, model1, 0);
    const std::string stdoutOutput = testing::internal::GetCapturedStdout();
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    // IMC has the same contract: selecting IMC must not silently run PDR.  On
    // this small model IMC may reach its own result directly, so the assertion
    // is about engine isolation rather than a forced coverage shape.
    EXPECT_EQ(imcResult.totalOutputs, 3u);
    EXPECT_EQ(stdoutOutput.find("PDR repair"), std::string::npos);
    EXPECT_EQ(stderrOutput.find("PDR repair"), std::string::npos);
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailFindsResetlessTransitionMismatch) {
  constexpr size_t kDummyStatesPerDesign = 1024;
  const auto models =
      makeDelayedRailMismatchModelsForTest(kDummyStatesPerDesign);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 4);

  // Both arbitrary initial values reach opposite constants after one step, so
  // exact PDR must report the observable cycle-1 mismatch.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_FALSE(result.reason.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailReportsTransientFrameZeroMismatch) {
  constexpr const char* kPrefix = "dualRailTransientStartupMismatch";
  auto models = makeHeldRailModelsForTest(kPrefix, false, true);
  const SignalKey state0 = makeSignalKey(std::string(kPrefix) + "State0");
  const SignalKey state1 = makeSignalKey(std::string(kPrefix) + "State1");
  models.model0.nextStateExprByStateKey.at(state0) = BoolExpr::createFalse();
  models.model1.nextStateExprByStateKey.at(state1) = BoolExpr::createFalse();

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(
      models.model0, models.model1, 2);

  // The cycle-zero values are binary-defined and opposite, even though both
  // designs transition to the same value afterward.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailFindsMismatchAlongsideXTrace) {
  constexpr const char* kPrefix = "dualRailMixedStartup";
  auto models = makeHeldRailModelsForTest(kPrefix, std::nullopt, true);
  const SignalKey input = makeSignalKey("dualRailMixedStartupInput");
  const SignalKey output = makeSignalKey(std::string(kPrefix) + "Output");
  const SignalKey state0 = makeSignalKey(std::string(kPrefix) + "State0");
  const SignalKey state1 = makeSignalKey(std::string(kPrefix) + "State1");
  for (SequentialDesignModel* model : {&models.model0, &models.model1}) {
    model->environmentInputs = {input};
    model->inputVarByKey.emplace(input, 3);
    model->displayNameByKey.emplace(input, "select[0]");
  }
  models.model0.nextStateExprByStateKey.at(state0) = BoolExpr::createFalse();
  models.model1.nextStateExprByStateKey.at(state1) = BoolExpr::createFalse();
  models.model0.observedOutputExprByKey.at(output) =
      BoolExpr::Or(BoolExpr::Var(2), BoolExpr::Var(3));
  models.model1.observedOutputExprByKey.at(output) =
      BoolExpr::And(BoolExpr::Var(3), BoolExpr::Not(BoolExpr::Var(2)));

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 2);

  // At cycle zero select=0 exposes X/0, while select=1 produces the defined
  // mismatch 1/0. The X trace must not hide the real counterexample.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailReportsDefinedBinaryMismatch) {
  const auto models = makeHeldRailModelsForTest(
      "dualRailDefinedMismatch", false, true);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 2);

  // Both sides are binary at frame zero: 01 versus 10 is a real mismatch.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDoesNotMineCrossDesignStateEqualities) {
  constexpr size_t kStateBitsPerDesign = 4097;
  const SignalKey out = makeSignalKey("smallOutputLargeStateOut");
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model0.displayNameByKey.emplace(out, "small_output_large_state_out[0]");
  model1.displayNameByKey.emplace(out, "small_output_large_state_out[0]");
  model0.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());
  model1.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());

  size_t nextLocalVar = 2;
  for (size_t i = 0; i < kStateBitsPerDesign; ++i) {
    const SignalKey state0 =
        makeSignalKey("smallOutputLargeStateLeft" + std::to_string(i));
    const SignalKey state1 =
        makeSignalKey("smallOutputLargeStateRight" + std::to_string(i));
    const size_t var0 = nextLocalVar++;
    const size_t var1 = nextLocalVar++;
    addStateBitForTest(
        model0,
        state0,
        var0,
        "small_output_large_state_q[" + std::to_string(i) + "]",
        BoolExpr::Var(var0));
    addStateBitForTest(
        model1,
        state1,
        var1,
        "small_output_large_state_q[" + std::to_string(i) + "]",
        BoolExpr::Var(var1));
  }

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // SEC engines must not assume that internal state from different designs is
  // related. With precoverage removed, this max_k=0 PDR run may remain
  // inconclusive, but it must still avoid cross-design state mining.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_EQ(
      stderrOutput.find("cross-design internal state equality mining"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("SEC diag: inferring inductive state equalities"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsImcSkipsPreEngineStateMining) {
  const SignalKey out = makeSignalKey("imcNoPreMiningOut");
  const SignalKey state0 = makeSignalKey("imcNoPreMiningState0");
  const SignalKey state1 = makeSignalKey("imcNoPreMiningState1");

  SequentialDesignModel model0;
  SequentialDesignModel model1;
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model0.displayNameByKey.emplace(out, "imc_no_pre_mining_out[0]");
  model1.displayNameByKey.emplace(out, "imc_no_pre_mining_out[0]");
  model0.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());
  model1.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());
  addStateBitForTest(
      model0,
      state0,
      /*symbol=*/2,
      "imc_no_pre_mining_q0[0]",
      BoolExpr::Var(2));
  addStateBitForTest(
      model1,
      state1,
      /*symbol=*/3,
      "imc_no_pre_mining_q1[0]",
      BoolExpr::Var(3));

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Classic IMC should start from interpolation/reachability, never from
  // assumed cross-design internal-state relations before the engine runs.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(
      stderrOutput.find("cross-design internal state equality mining"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("SEC diag: inferring inductive state equalities"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsKInductionDoesNotMineCrossDesignStateEqualities) {
  const SignalKey out = makeSignalKey("kiNoCrossStateOut");
  const SignalKey state0 = makeSignalKey("kiNoCrossState0");
  const SignalKey state1 = makeSignalKey("kiNoCrossState1");

  SequentialDesignModel model0;
  SequentialDesignModel model1;
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model0.displayNameByKey.emplace(out, "ki_no_cross_state_out[0]");
  model1.displayNameByKey.emplace(out, "ki_no_cross_state_out[0]");
  model0.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());
  model1.observedOutputExprByKey.emplace(out, BoolExpr::createTrue());
  addStateBitForTest(
      model0, state0, /*symbol=*/2, "ki_no_cross_state_q0[0]", BoolExpr::Var(2));
  addStateBitForTest(
      model1, state1, /*symbol=*/3, "ki_no_cross_state_q1[0]", BoolExpr::Var(3));

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(
      stderrOutput.find("cross-design internal state equality mining"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find("SEC diag: inferring inductive state equalities"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailFindsWideFrameZeroCounterexampleBeforeSkip) {
  constexpr size_t kObservedOutputs = 133;
  constexpr size_t kDummyStatesPerDesign = 4100;
  const auto models = makeWideFrameZeroMismatchModelsForTest(
      kObservedOutputs, kDummyStatesPerDesign);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 1);

  // The unrelated dummy flops keep this out of the small leaf certificate path.
  // PDR must still honor a concrete top-output frame-0 mismatch before any
  // residual output can be skipped as an inconclusive hard leaf.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, kObservedOutputs);
  EXPECT_EQ(result.totalOutputs, kObservedOutputs);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsKInductionDualRailFindsWideFrameZeroCounterexampleBeforeResidualProof) {
  constexpr size_t kObservedOutputs = 133;
  constexpr size_t kDummyStatesPerDesign = 4100;
  const auto models = makeWideFrameZeroMismatchModelsForTest(
      kObservedOutputs, kDummyStatesPerDesign);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 1);

  // KI residual batching must not bury a concrete input-only top-output edit
  // behind an expensive residual proof.  The pre-batch witness check keeps the
  // counterexample in the selected SEC engine path without using LEC fallback.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, kObservedOutputs);
  EXPECT_EQ(result.totalOutputs, kObservedOutputs);
  EXPECT_NE(
      result.reason.find("wide_frame_zero_probe[0]"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsImcDualRailFindsWideFrameZeroCounterexampleBeforeCraigBatching) {
  constexpr size_t kObservedOutputs = 133;
  constexpr size_t kDummyStatesPerDesign = 4100;
  const auto models = makeWideFrameZeroMismatchModelsForTest(
      kObservedOutputs, kDummyStatesPerDesign);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 1);

  // Large dual-rail IMC uses Craig batching for proof work, but concrete
  // frame-0 edits still have to be found by the selected IMC base query before
  // an early inconclusive batch can hide a later edited output.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, kObservedOutputs);
  EXPECT_EQ(result.totalOutputs, kObservedOutputs);
  EXPECT_NE(
      result.reason.find("wide_frame_zero_probe[0]"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsImcDualRailLeavesResetlessStateLeafInconclusive) {
  constexpr size_t kDummyStatesPerDesign = 1024;
  const auto models =
      makeDelayedRailMismatchModelsForTest(kDummyStatesPerDesign);

  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc,
      SecEncoding::DualRailSteady);
  const auto result =
      strategy.runExtractedModels(models.model0, models.model1, 4);

  // Both all-X startup states reach opposite binary constants after one
  // transition. The guarded steady-state property therefore has a real public
  // output mismatch at cycle 1, independent of any internal-state relation.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_NE(result.reason.find("delayed_out[0]"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsSkipsOnlyResetUnanchoredOutputInsteadOfComparingState) {
  const SignalKey rst = makeSignalKey("keepOnlyResetUnanchoredRst");
  const SignalKey out = makeSignalKey("keepOnlyResetUnanchoredOut");
  const SignalKey state0 = makeSignalKey("keepOnlyResetUnanchoredState0");
  const SignalKey state1 = makeSignalKey("keepOnlyResetUnanchoredState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst};
  model0.stateBits = {state0};
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model0.inputVarByKey.emplace(rst, 2);
  model0.inputVarByKey.emplace(state0, 4);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(out, "only_out[0]");
  model0.displayNameByKey.emplace(state0, "u_left.q[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));
  model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(4));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst};
  model1.stateBits = {state1};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model1.inputVarByKey.emplace(rst, 3);
  model1.inputVarByKey.emplace(state1, 5);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(out, "only_out[0]");
  model1.displayNameByKey.emplace(state1, "u_right.q[0]");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::createFalse());
  model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(5));

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::KInduction);
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  // Binary SEC may only compare top outputs.  If a top output's first
  // post-reset value is still an arbitrary internal flop value in each design,
  // do not restore the output just to avoid zero coverage: that would compare an
  // internal design0 state bit with an internal design1 state bit by implication.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, 1u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("only_out[0]"), std::string::npos);
  EXPECT_NE(
      result.skippedObservedOutputs.front().find("reset-unanchored"),
      std::string::npos);
  ASSERT_EQ(result.resetUnanchoredSkippedOutputs.size(), 1u);
  EXPECT_EQ(
      result.resetUnanchoredSkippedOutputs.front(),
      result.skippedObservedOutputs.front());
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsRejectsWideStartupRelationWithoutStateRelations) {
  const SignalKey rst = makeSignalKey("wideStartupRelationRst");
  const SignalKey stateA0 = makeSignalKey("wideStartupRelationStateA0");
  const SignalKey stateB0 = makeSignalKey("wideStartupRelationStateB0");
  const SignalKey stateA1 = makeSignalKey("wideStartupRelationStateA1");
  const SignalKey stateB1 = makeSignalKey("wideStartupRelationStateB1");
  const SignalKey stateC1 = makeSignalKey("wideStartupRelationStateC1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst};
  model0.stateBits = {stateA0, stateB0};
  model0.inputVarByKey.emplace(rst, 2);
  model0.inputVarByKey.emplace(stateA0, 4);
  model0.inputVarByKey.emplace(stateB0, 6);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(stateA0, "u_left.a_q[0]");
  model0.displayNameByKey.emplace(stateB0, "u_left.b_q[0]");
  model0.nextStateExprByStateKey.emplace(stateA0, BoolExpr::Var(6));
  model0.nextStateExprByStateKey.emplace(
      stateB0,
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::createFalse()));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst};
  model1.stateBits = {stateA1, stateB1, stateC1};
  model1.inputVarByKey.emplace(rst, 3);
  model1.inputVarByKey.emplace(stateA1, 5);
  model1.inputVarByKey.emplace(stateB1, 7);
  model1.inputVarByKey.emplace(stateC1, 9);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(stateA1, "u_right.a_q[0]");
  model1.displayNameByKey.emplace(stateB1, "u_right.b_q[0]");
  model1.displayNameByKey.emplace(stateC1, "u_right.extra_q[0]");
  model1.nextStateExprByStateKey.emplace(stateA1, BoolExpr::Var(7));
  model1.nextStateExprByStateKey.emplace(stateB1, BoolExpr::createFalse());
  model1.nextStateExprByStateKey.emplace(stateC1, BoolExpr::createFalse());

  constexpr size_t kWideStartupRelationOutputs = 129;
  for (size_t i = 0; i < kWideStartupRelationOutputs; ++i) {
    const SignalKey out =
        makeSignalKey("wideStartupRelationOut" + std::to_string(i));
    const std::string outputName =
        "wide_startup_relation_out[" + std::to_string(i) + "]";
    model0.allObservedOutputs.push_back(out);
    model0.observedOutputs.push_back(out);
    model0.displayNameByKey.emplace(out, outputName);
    model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(4));
    model1.allObservedOutputs.push_back(out);
    model1.observedOutputs.push_back(out);
    model1.displayNameByKey.emplace(out, outputName);
    model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(5));
  }

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::Pdr);
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  // Wide reset/bootstrap state-dependent outputs cannot be proved by a
  // cross-design startup relation.  Binary SEC reports the conservative skip
  // surface instead of proving outputs from internal element names.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, kWideStartupRelationOutputs);
  EXPECT_EQ(result.skippedObservedOutputs.size(), kWideStartupRelationOutputs);
  EXPECT_NE(
      result.skippedObservedOutputs.front().find("reset-unanchored internal state"),
      std::string::npos);

  SequentialEquivalenceStrategy dualRailKiStrategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto dualRailKiResult =
      dualRailKiStrategy.runExtractedModels(model0, model1, 1);
  EXPECT_EQ(dualRailKiResult.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(dualRailKiResult.bound, 1u);
  EXPECT_EQ(dualRailKiResult.coveredOutputs, 0u);
  EXPECT_EQ(dualRailKiResult.totalOutputs, kWideStartupRelationOutputs);
  EXPECT_FALSE(dualRailKiResult.skippedObservedOutputs.empty());

  SequentialEquivalenceStrategy dualRailImcStrategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Imc,
      SecEncoding::DualRailSteady);
  const auto dualRailImcResult =
      dualRailImcStrategy.runExtractedModels(model0, model1, 1);
  // Without an initial predicate, IMC must not turn an X-only startup relation
  // into an equivalence proof or count guarded equality as output coverage.
  EXPECT_EQ(dualRailImcResult.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(dualRailImcResult.bound, 1u);
  EXPECT_EQ(dualRailImcResult.coveredOutputs, 0u);
  EXPECT_EQ(dualRailImcResult.totalOutputs, kWideStartupRelationOutputs);
  EXPECT_EQ(
      dualRailImcResult.skippedObservedOutputs.size(),
      kWideStartupRelationOutputs);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsKeepsResetCombinationalOutputsCovered) {
  const SignalKey rst = makeSignalKey("keepResetCombRst");
  const SignalKey data = makeSignalKey("keepResetCombData");
  const SignalKey out = makeSignalKey("keepResetCombOut");
  const SignalKey state0 = makeSignalKey("keepResetCombState0");
  const SignalKey state1 = makeSignalKey("keepResetCombState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst, data};
  model0.stateBits = {state0};
  model0.allObservedOutputs = {out};
  model0.observedOutputs = {out};
  model0.inputVarByKey.emplace(rst, 2);
  model0.inputVarByKey.emplace(data, 4);
  model0.inputVarByKey.emplace(state0, 6);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(data, "data[0]");
  model0.displayNameByKey.emplace(out, "out[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Var(6));
  model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(4));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst, data};
  model1.stateBits = {state1};
  model1.allObservedOutputs = {out};
  model1.observedOutputs = {out};
  model1.inputVarByKey.emplace(rst, 3);
  model1.inputVarByKey.emplace(data, 5);
  model1.inputVarByKey.emplace(state1, 7);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(data, "data[0]");
  model1.displayNameByKey.emplace(out, "out[0]");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(7)));
  model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(5));

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::Pdr);
  // Unconstrained startup state is represented exactly by F[0] = true. The
  // combinational output is independent of that state, so exact PDR proves it.
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_TRUE(result.skippedObservedOutputs.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsStopsOnUnsupportedFirstModelWithBoundaryReports) {
  auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  auto model1 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  const SignalKey stateKey = makeSignalKey("state");
  const SignalKey internalIn = makeSignalKey("internal_in");
  const SignalKey internalOut = makeSignalKey("internal_out");
  model0.unsupportedReasons = {"unsupported sequential state"};
  model0.abstractedSequentialBoundaries = {"abstracted cell u_ff"};
  model0.internalBoundaryInputKeys = {internalIn};
  model0.internalBoundaryOutputKeys = {internalOut};
  model0.displayNameByKey.emplace(stateKey, "u_ff.STATE[0]");
  model0.displayNameByKey.emplace(internalIn, "u_logic.A[0]");
  model0.displayNameByKey.emplace(internalOut, "u_logic.Y[0]");
  model0.abstractedSequentialBoundaryDetails.push_back(
      {"u_ff", {stateKey}, model0.allObservedOutputs});

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("unsupported sequential state"), std::string::npos);
  ASSERT_EQ(result.abstractedSequentialBoundaries.size(), 1u);
  EXPECT_EQ(result.abstractedSequentialBoundaries.front(),
            "design0 abstracted cell u_ff");
  EXPECT_GE(result.extractedBoundaryReports.size(), 4u);
  const auto stateReport = std::find_if(
      result.extractedBoundaryReports.begin(),
      result.extractedBoundaryReports.end(),
      [](const ExtractedBoundaryReportEntry& entry) {
        return entry.signal == "u_ff.STATE[0]";
      });
  ASSERT_NE(stateReport, result.extractedBoundaryReports.end());
  EXPECT_NE(
      std::find(
          stateReport->roles.begin(),
          stateReport->roles.end(),
          "abstracted_sequential_state"),
      stateReport->roles.end());
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsStopsOnUnsupportedSecondModelAfterFirstReports) {
  auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  auto model1 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  model0.abstractedSequentialBoundaries = {"kept first-side boundary"};
  model1.unsupportedReasons = {"unsupported second side"};

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_NE(result.reason.find("unsupported second side"), std::string::npos);
  ASSERT_EQ(result.abstractedSequentialBoundaries.size(), 1u);
  EXPECT_EQ(result.abstractedSequentialBoundaries.front(),
            "design0 kept first-side boundary");
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsReportsAllConnectivitySkippedOutputs) {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  std::array<ConnectivitySkipOrigin, 4> origins = {
      ConnectivitySkipOrigin::NoDriver,
      ConnectivitySkipOrigin::MultiDriver,
      ConnectivitySkipOrigin::LogicalLoop,
      ConnectivitySkipOrigin::MultiClockDomain};
  for (size_t i = 0; i < origins.size(); ++i) {
    const SignalKey key = makeSignalKey("skipped_out_" + std::to_string(i));
    const std::string name = "out" + std::to_string(i) + "[0]";
    model0.allObservedOutputs.push_back(key);
    model1.allObservedOutputs.push_back(key);
    model0.displayNameByKey.emplace(key, name);
    model1.displayNameByKey.emplace(key, name);
    model0.connectivitySkipInfoByKey.emplace(
        key, ConnectivitySkipInfo{origins[i], "left side"});
    model1.connectivitySkipInfoByKey.emplace(
        key, ConnectivitySkipInfo{origins[i], "right side"});
  }

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, origins.size());
  ASSERT_EQ(result.skippedObservedOutputs.size(), origins.size());
  const auto hasSkipText = [&](const char* text) {
    return std::any_of(
        result.skippedObservedOutputs.begin(),
        result.skippedObservedOutputs.end(),
        [&](const std::string& skipped) {
          return skipped.find(text) != std::string::npos;
        });
  };
  EXPECT_TRUE(hasSkipText("no-driver connectivity"));
  EXPECT_TRUE(hasSkipText("multi-driver connectivity"));
  EXPECT_TRUE(hasSkipText("logical-loop connectivity"));
  EXPECT_TRUE(hasSkipText("multi-clock-domain connectivity"));
  ASSERT_EQ(result.multiClockDomainSkippedOutputs.size(), 1u);
  EXPECT_NE(
      result.multiClockDomainSkippedOutputs.front().find("multi-clock-domain"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsThrowsOnMissingOutputExpression) {
  auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  auto model1 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  model1.observedOutputExprByKey.clear();

  auto strategy = makeBinaryExtractedSecStrategy();
  EXPECT_THROW(
      static_cast<void>(strategy.runExtractedModels(model0, model1, 1)),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsIgnoresStaleObservedOutputsOutsideCoverage) {
  auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  auto model1 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  const SignalKey extraKey = makeSignalKey("extra_observed");
  model0.observedOutputs.push_back(extraKey);
  model1.observedOutputs.push_back(extraKey);
  model0.displayNameByKey.emplace(extraKey, "extra[0]");
  model1.displayNameByKey.emplace(extraKey, "extra[0]");

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsAcceptsSameValueModelWithoutBuildingSatProblem) {
  const auto model = makeCombinationalExtractedModel(BoolExpr::Var(2));

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(model, model, 9);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_EQ(result.outputCoveragePercent(), 100.0);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsEmitsSelectedEngineDiagnostics) {
  const auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  const auto model1 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  const std::array<std::pair<SecEngine, const char*>, 3> expected = {{
      {SecEngine::Pdr, "pdr engine"},
      {SecEngine::Imc, "imc engine"},
      {SecEngine::KInduction, "classic k-induction engine"},
  }};

  for (const auto& [engine, label] : expected) {
    testing::internal::CaptureStderr();
    auto strategy = makeBinaryExtractedSecStrategy(engine);
    const auto result = strategy.runExtractedModels(model0, model1, 1);
    const std::string stderrOutput = testing::internal::GetCapturedStderr();

    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_NE(stderrOutput.find(label), std::string::npos);
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsFormatsCompactCounterexampleWithoutDnlTraceback) {
  const auto model0 = makeCombinationalExtractedModel(BoolExpr::Var(2));
  const auto model1 = makeCombinationalExtractedModel(BoolExpr::createTrue());

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::KInduction);
  const auto result = strategy.runExtractedModels(model0, model1, 0);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_NE(result.reason.find("Counterexample reaches"), std::string::npos);
  EXPECT_NE(result.reason.find("compact SEC released"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFindsInitialBadStateBeforeGrowingFrames) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::Var(2);
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar secPdrTrace("KEPLER_SEC_PDR_TRACE", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 0u);
  EXPECT_NE(stderrOutput.find("bad_cube@F0"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineNoCopyPathPreservesDefaultBadNormalization) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.initialCondition = BoolExpr::Var(2);
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  problem.property = BoolExpr::Not(BoolExpr::Var(2));
  // Internal callers historically had their bad root normalized by run().
  // An inconsistent stored root must therefore use the copy fallback.
  problem.bad = BoolExpr::createFalse();

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBuildsExactBootstrapFrameZeroFromResetPrefix) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  // F[0] is incomplete because the bootstrap summary does not assign y.
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(3);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_GE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineKeepsExactRelationalBootstrapState) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t reset = 4;
  constexpr size_t data = 5;
  problem.state0Symbols = {x, y};
  problem.inputSymbols = {reset, data};
  problem.allSymbols = {x, y, reset, data};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};

  // Reset loads both registers from the same arbitrary input. After reset the
  // registers swap, so the exact x == y frontier is invariant. Treating the
  // incomplete constant summary as F[0] would admit 01 and report a spurious
  // 01 -> 10 counterexample.
  BoolExpr* resetAndData =
      BoolExpr::And(BoolExpr::Var(reset), BoolExpr::Var(data));
  BoolExpr* resetInactive = BoolExpr::Not(BoolExpr::Var(reset));
  problem.transitions0.emplace_back(
      x,
      BoolExpr::Or(
          resetAndData,
          BoolExpr::And(resetInactive, BoolExpr::Var(y))));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::Or(
          resetAndData,
          BoolExpr::And(resetInactive, BoolExpr::Var(x))));
  problem.bad =
      BoolExpr::And(BoolExpr::Var(x), BoolExpr::Not(BoolExpr::Var(y)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_GE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesResetImageDespiteCompleteDualRailBootstrapSummary) {
  KInductionProblem problem;
  constexpr size_t xOne = 2;
  constexpr size_t xZero = 3;
  constexpr size_t yOne = 4;
  constexpr size_t yZero = 5;
  constexpr size_t reset = 6;
  constexpr size_t data = 7;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {xOne, xZero, yOne, yZero};
  problem.inputSymbols = {reset, data};
  problem.allSymbols = {xOne, xZero, yOne, yZero, reset, data};
  problem.totalStateCount = problem.state0Symbols.size();
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.dualRailStatePairs = {
      DualRailSymbolPair{xOne, xZero},
      DualRailSymbolPair{yOne, yZero}};

  // A full-size ternary summary may lose reset-created relations: X versus 0
  // appears different even though reset loads both registers from the same
  // public input. F[0] must therefore come from the reset transition image.
  problem.bootstrapStateAssignments = {
      {xOne, true}, {xZero, true}, {yOne, false}, {yZero, true}};
  BoolExpr* resetInactive = BoolExpr::Not(BoolExpr::Var(reset));
  problem.transitions0 = {
      {xOne,
       BoolExpr::Or(
           BoolExpr::And(BoolExpr::Var(reset), BoolExpr::Var(data)),
           BoolExpr::And(resetInactive, BoolExpr::Var(xOne)))},
      {xZero,
       BoolExpr::Or(
           BoolExpr::And(
               BoolExpr::Var(reset), BoolExpr::Not(BoolExpr::Var(data))),
           BoolExpr::And(resetInactive, BoolExpr::Var(xZero)))},
      {yOne,
       BoolExpr::Or(
           BoolExpr::And(BoolExpr::Var(reset), BoolExpr::Var(data)),
           BoolExpr::And(resetInactive, BoolExpr::Var(yOne)))},
      {yZero,
       BoolExpr::Or(
           BoolExpr::And(
               BoolExpr::Var(reset), BoolExpr::Not(BoolExpr::Var(data))),
           BoolExpr::And(resetInactive, BoolExpr::Var(yZero)))}};
  problem.bad = BoolExpr::Or(
      BoolExpr::Xor(BoolExpr::Var(xOne), BoolExpr::Var(yOne)),
      BoolExpr::Xor(BoolExpr::Var(xZero), BoolExpr::Var(yZero)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_GE(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineResetPrefixStartsFromDualRailXInitialization) {
  KInductionProblem problem;
  constexpr size_t xOne = 2;
  constexpr size_t xZero = 3;
  constexpr size_t yOne = 4;
  constexpr size_t yZero = 5;
  constexpr size_t reset = 6;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {xOne, xZero};
  problem.state1Symbols = {yOne, yZero};
  problem.inputSymbols = {reset};
  problem.allSymbols = {xOne, xZero, yOne, yZero, reset};
  problem.totalStateCount = 4;
  problem.initializedStateCount = 4;
  problem.initialStateAssignments = {
      {xOne, true}, {xZero, true}, {yOne, true}, {yZero, true}};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.dualRailStatePairs = {
      DualRailSymbolPair{xOne, xZero},
      DualRailSymbolPair{yOne, yZero}};
  problem.transitions0 = {
      {xOne, BoolExpr::Var(xOne)}, {xZero, BoolExpr::Var(xZero)}};
  problem.transitions1 = {
      {yOne, BoolExpr::Var(yOne)}, {yZero, BoolExpr::Var(yZero)}};

  // This per-register summary claims 1 versus 0, but the exact reset prefix
  // preserves the real X versus X initialization through both hold flops.
  problem.bootstrapStateAssignments = {
      {xOne, true}, {xZero, false}, {yOne, false}, {yZero, true}};
  BoolExpr* bothDefined = BoolExpr::And(
      BoolExpr::Xor(BoolExpr::Var(xOne), BoolExpr::Var(xZero)),
      BoolExpr::Xor(BoolExpr::Var(yOne), BoolExpr::Var(yZero)));
  problem.bad = BoolExpr::And(
      bothDefined,
      BoolExpr::Xor(BoolExpr::Var(xOne), BoolExpr::Var(yOne)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionBad = problem.bad;
  problem.inductionProperty = problem.property;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSharesExactDualRailInitWithoutSharingBatchResults) {
  KInductionProblem sourceProblem;
  constexpr size_t xOne = 2;
  constexpr size_t xZero = 3;
  constexpr size_t yOne = 4;
  constexpr size_t yZero = 5;
  constexpr size_t reset = 6;
  sourceProblem.usesDualRailStateEncoding = true;
  sourceProblem.state0Symbols = {xOne, xZero};
  sourceProblem.state1Symbols = {yOne, yZero};
  sourceProblem.inputSymbols = {reset};
  sourceProblem.allSymbols = {xOne, xZero, yOne, yZero, reset};
  sourceProblem.totalStateCount = 4;
  sourceProblem.initializedStateCount = 4;
  sourceProblem.initialStateAssignments = {
      {xOne, true}, {xZero, true}, {yOne, true}, {yZero, true}};
  sourceProblem.resetBootstrapCycles = 1;
  sourceProblem.resetBootstrapInputs = {{reset, true}};
  sourceProblem.dualRailStatePairs = {
      DualRailSymbolPair{xOne, xZero},
      DualRailSymbolPair{yOne, yZero}};
  sourceProblem.transitions0 = {
      {xOne, BoolExpr::Var(xOne)}, {xZero, BoolExpr::Var(xZero)}};
  sourceProblem.transitions1 = {
      {yOne, BoolExpr::Var(yOne)}, {yZero, BoolExpr::Var(yZero)}};

  BoolExpr* bothDefined = BoolExpr::And(
      BoolExpr::Xor(BoolExpr::Var(xOne), BoolExpr::Var(xZero)),
      BoolExpr::Xor(BoolExpr::Var(yOne), BoolExpr::Var(yZero)));
  BoolExpr* guardedMismatch = BoolExpr::And(
      bothDefined,
      BoolExpr::Xor(BoolExpr::Var(xOne), BoolExpr::Var(yOne)));
  // The full source property reserves the support needed by either batch.
  sourceProblem.bad = BoolExpr::Or(guardedMismatch, BoolExpr::Var(xOne));
  sourceProblem.property = BoolExpr::Not(sourceProblem.bad);
  sourceProblem.inductionBad = sourceProblem.bad;
  sourceProblem.inductionProperty = sourceProblem.property;
  sourceProblem.observedOutputExprs0 = {BoolExpr::Var(xOne)};
  sourceProblem.observedOutputExprs1 = {BoolExpr::Var(yOne)};
  sourceProblem.dualRailOutputStrictEqualityExprs = {BoolExpr::Var(xOne)};

  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      sourceProblem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  KInductionProblem equivalentBatch = sourceProblem;
  equivalentBatch.bad = guardedMismatch;
  equivalentBatch.property = BoolExpr::Not(equivalentBatch.bad);
  equivalentBatch.inductionBad = equivalentBatch.bad;
  equivalentBatch.inductionProperty = equivalentBatch.property;
  KInductionProblem secondEquivalentBatch = sourceProblem;
  secondEquivalentBatch.bad = BoolExpr::Not(BoolExpr::Var(xOne));
  secondEquivalentBatch.property = BoolExpr::Not(secondEquivalentBatch.bad);
  secondEquivalentBatch.inductionBad = secondEquivalentBatch.bad;
  secondEquivalentBatch.inductionProperty =
      secondEquivalentBatch.property;
  KInductionProblem differentBatch = sourceProblem;
  differentBatch.bad = BoolExpr::Var(xOne);
  differentBatch.property = BoolExpr::Not(differentBatch.bad);
  differentBatch.inductionBad = differentBatch.bad;
  differentBatch.inductionProperty = differentBatch.property;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine equivalentEngine(
      equivalentBatch,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      0,
      exactInitCache);
  const auto equivalentResult = equivalentEngine.run(2);
  PDREngine secondEquivalentEngine(
      secondEquivalentBatch,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      0,
      exactInitCache);
  const auto secondEquivalentResult = secondEquivalentEngine.run(2);
  PDREngine differentEngine(
      differentBatch,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      0,
      exactInitCache);
  const auto differentResult = differentEngine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(equivalentResult.status, PDRStatus::Equivalent);
  EXPECT_EQ(secondEquivalentResult.status, PDRStatus::Equivalent);
  EXPECT_EQ(differentResult.status, PDRStatus::Different);
  EXPECT_EQ(differentResult.bound, 0u);
  const std::string built = "shared exact F[0] cache built";
  const size_t firstBuild = stderrOutput.find(built);
  ASSERT_NE(firstBuild, std::string::npos) << stderrOutput;
  EXPECT_EQ(stderrOutput.find(built, firstBuild + built.size()),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("shared exact F[0] cache reused"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared exact F[0] solver used for bad cube"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared exact F[0] solver used for init intersection"),
      std::string::npos)
      << stderrOutput;
  const std::string metadataBuilt = "immutable model metadata built";
  const size_t firstMetadataBuild = stderrOutput.find(metadataBuilt);
  ASSERT_NE(firstMetadataBuild, std::string::npos) << stderrOutput;
  EXPECT_EQ(stderrOutput.find(metadataBuilt,
                              firstMetadataBuild + metadataBuilt.size()),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("immutable model metadata reused"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("reusable invariant candidates stored"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("reusable invariant clauses certified"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("reusable invariant clauses injected"),
            std::string::npos)
      << stderrOutput;
  const std::string predecessorCreated =
      "predecessor cached solver created level=0";
  const size_t firstPredecessorBuild = stderrOutput.find(predecessorCreated);
  ASSERT_NE(firstPredecessorBuild, std::string::npos) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          predecessorCreated,
          firstPredecessorBuild + predecessorCreated.size()),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared exact F[0] predecessor solver reused"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRunsIndependentPropertiesWithVerifierOwnedState) {
  KInductionProblem problem;
  constexpr size_t designState = 2;
  constexpr size_t monitorState = 3;
  problem.state0Symbols = {designState};
  problem.auxiliaryStateSymbols = {monitorState};
  problem.allSymbols = {designState, monitorState};
  problem.totalStateCount = 2;
  problem.initializedStateCount = 2;
  problem.initialStateAssignments = {
      {designState, false}, {monitorState, false}};
  problem.transitions0 = {{designState, BoolExpr::createTrue()}};
  problem.auxiliaryTransitions = {
      {monitorState, BoolExpr::createTrue()}};
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();

  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);

  BoolExpr* holdsAfterMonitor = BoolExpr::Or(
      BoolExpr::Not(BoolExpr::Var(monitorState)),
      BoolExpr::Var(designState));
  BoolExpr* failsAfterMonitor = BoolExpr::Or(
      BoolExpr::Not(BoolExpr::Var(monitorState)),
      BoolExpr::Not(BoolExpr::Var(designState)));

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  const auto proved = engine.run(2, holdsAfterMonitor);
  const auto different = engine.run(2, failsAfterMonitor);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Certified reachability invariants may be shared, but each supplied safety
  // property still gets an independent bad-state search and verdict.
  EXPECT_EQ(proved.status, PDRStatus::Equivalent);
  EXPECT_EQ(different.status, PDRStatus::Different);
  EXPECT_EQ(different.bound, 1u);
  EXPECT_NE(stderrOutput.find("reusable invariant clauses certified"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("reusable invariant clauses injected"),
            std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRejectsUnfinishedFrameClausesBeforeCrossPropertyReuse) {
  KInductionProblem problem;
  constexpr size_t firstState = 2;
  constexpr size_t delayedState = 3;
  constexpr size_t otherFirstState = 4;
  constexpr size_t otherDelayedState = 5;
  problem.state0Symbols = {
      firstState, delayedState, otherFirstState, otherDelayedState};
  problem.allSymbols = problem.state0Symbols;
  problem.totalStateCount = 4;
  problem.initializedStateCount = 4;
  problem.initialStateAssignments = {
      {firstState, false},
      {delayedState, false},
      {otherFirstState, false},
      {otherDelayedState, false}};
  problem.transitions0 = {
      {firstState, BoolExpr::createTrue()},
      {delayedState, BoolExpr::Var(firstState)},
      {otherFirstState, BoolExpr::createTrue()},
      {otherDelayedState, BoolExpr::Var(otherFirstState)}};
  BoolExpr* delayedIsFalse = BoolExpr::Not(BoolExpr::Var(delayedState));
  BoolExpr* otherDelayedIsFalse =
      BoolExpr::Not(BoolExpr::Var(otherDelayedState));
  problem.property = delayedIsFalse;
  problem.bad = BoolExpr::Not(delayedIsFalse);

  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  // !delayed is true initially but is not inductive: first=1, delayed=0
  // transitions to delayed'=1. The first run therefore leaves it only in F1.
  const auto unfinished = engine.run(1);
  const auto proved = engine.run(1, BoolExpr::createTrue());
  // A second independent unfinished clause must be certified alone. The
  // FMCAD'11 removal loop permanently discarded the first rejected candidate.
  const auto otherUnfinished = engine.run(1, otherDelayedIsFalse);
  const auto otherProved = engine.run(1, BoolExpr::createTrue());
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(unfinished.status, PDRStatus::Inconclusive);
  EXPECT_EQ(proved.status, PDRStatus::Equivalent);
  EXPECT_EQ(otherUnfinished.status, PDRStatus::Inconclusive);
  EXPECT_EQ(otherProved.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("reusable invariant candidates stored"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("reusable invariant clauses certified"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("inductive_rejected=1"),
      std::string::npos)
      << stderrOutput;
  const std::string singleRejectedCertification =
      "reusable invariant clauses certified candidates=1 retained=0 "
      "initial_rejected=0 inductive_rejected=1";
  const size_t firstCertification =
      stderrOutput.find(singleRejectedCertification);
  ASSERT_NE(firstCertification, std::string::npos) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          singleRejectedCertification,
          firstCertification + singleRejectedCertification.size()),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          "reusable invariant clauses certified candidates=2"),
      std::string::npos)
      << stderrOutput;
  // Direct Init contradictions do not need to enter the whole-model SAT owner.
  // Retaining those learned clauses across hundreds of batches once exhausted
  // CI memory on a large dual-rail design.
  EXPECT_NE(
      stderrOutput.find(
          "reusable invariant initial facts resolved=1 unresolved=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reusable invariant local init solver created"),
      std::string::npos)
      << stderrOutput;
  // Inductive certification is local to one candidate set. Retaining this SAT
  // owner across recursive batches caused whole-design transition cones to
  // accumulate beyond the CI memory limit.
  const std::string inductiveSolverCreated =
      "reusable invariant inductive local solver created";
  const size_t firstInductiveSolverCreation =
      stderrOutput.find(inductiveSolverCreated);
  ASSERT_NE(firstInductiveSolverCreation, std::string::npos) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          inductiveSolverCreated,
          firstInductiveSolverCreation + inductiveSolverCreated.size()),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reusable invariant inductive solver reused"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCumulativeCertificationBudgetDoesNotChangePropertyVerdict) {
  KInductionProblem problem;
  constexpr size_t firstState = 2;
  constexpr size_t delayedState = 3;
  problem.state0Symbols = {firstState, delayedState};
  problem.allSymbols = problem.state0Symbols;
  problem.totalStateCount = 2;
  problem.initializedStateCount = 2;
  problem.initialStateAssignments = {
      {firstState, false}, {delayedState, false}};
  problem.transitions0 = {
      {firstState, BoolExpr::createTrue()},
      {delayedState, BoolExpr::Var(firstState)}};
  BoolExpr* delayedIsFalse = BoolExpr::Not(BoolExpr::Var(delayedState));
  problem.property = delayedIsFalse;
  problem.bad = BoolExpr::Not(delayedIsFalse);
  problem.usesDualRailStateEncoding = true;

  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);

  EXPECT_EQ(engine.run(1).status, PDRStatus::Inconclusive);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDRQueryLimits limits{
      /*predecessorConflictLimit=*/100,
      /*predecessorDecisionLimit=*/100,
      /*blockingConflictLimit=*/100,
      /*blockingDecisionLimit=*/100,
      /*predecessorEncodingNodeLimit=*/0,
      /*predecessorNodeHintTargetLimit=*/0,
      /*invariantCertificationTotalTickLimit=*/0};
  const PDRResult boundedResult =
      engine.run(1, BoolExpr::createTrue(), limits);
  const auto retriedResult = engine.run(1, BoolExpr::createTrue());
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(boundedResult.status, PDRStatus::Equivalent);
  EXPECT_EQ(retriedResult.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "reusable invariant cumulative certification budget created "
          "ticks=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "reusable invariant cumulative certification budget exhausted; "
          "disabling optional certification"),
      std::string::npos)
      << stderrOutput;
  // The later run must not retry pending candidates after the shared budget
  // is exhausted, and no uncertified clause may become reusable.
  EXPECT_EQ(
      stderrOutput.find("reusable invariant clauses certified"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineGuardsSharedHigherFrameSolverAcrossProperties) {
  KInductionProblem problem;
  constexpr size_t firstState = 2;
  constexpr size_t delayedState = 3;
  constexpr size_t otherFirstState = 4;
  constexpr size_t otherDelayedState = 5;
  constexpr size_t independentFirstState = 6;
  constexpr size_t independentDelayedState = 7;
  problem.state0Symbols = {
      firstState,
      delayedState,
      otherFirstState,
      otherDelayedState,
      independentFirstState,
      independentDelayedState};
  problem.allSymbols = problem.state0Symbols;
  problem.totalStateCount = 6;
  problem.initializedStateCount = 6;
  problem.initialStateAssignments = {
      {firstState, false},
      {delayedState, false},
      {otherFirstState, false},
      {otherDelayedState, false},
      {independentFirstState, false},
      {independentDelayedState, false}};
  problem.transitions0 = {
      {firstState, BoolExpr::createTrue()},
      {delayedState, BoolExpr::Var(firstState)},
      {otherFirstState, BoolExpr::createTrue()},
      {otherDelayedState, BoolExpr::Var(otherFirstState)},
      {independentFirstState, BoolExpr::createTrue()},
      {independentDelayedState, BoolExpr::Var(independentFirstState)}};
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.usesStrictDualRailEqualityProperty = true;

  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      exactInitCache);

  BoolExpr* safeProperty = BoolExpr::Or(
      BoolExpr::Not(BoolExpr::Var(delayedState)),
      BoolExpr::Var(firstState));
  BoolExpr* otherSafeProperty = BoolExpr::Or(
      BoolExpr::Not(BoolExpr::Var(otherDelayedState)),
      BoolExpr::Var(otherFirstState));
  BoolExpr* combinedSafeProperty =
      BoolExpr::And(safeProperty, otherSafeProperty);
  BoolExpr* independentSafeProperty = BoolExpr::Or(
      BoolExpr::Not(BoolExpr::Var(independentDelayedState)),
      BoolExpr::Var(independentFirstState));
  BoolExpr* depthTwoFailure = BoolExpr::Not(BoolExpr::Var(delayedState));
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");

  testing::internal::CaptureStderr();
  // The first property creates a path-local solver. A disjoint property gets a
  // separate entry, while certified invariants can avoid later solver calls.
  const auto combinedProved = engine.run(3, combinedSafeProperty);
  const auto subsetProved = engine.run(3, safeProperty);
  const auto independentProved = engine.run(3, independentSafeProperty);
  const auto different = engine.run(3, depthTwoFailure);
  std::vector<PDRResult> repeatedSubsetResults;
  repeatedSubsetResults.reserve(20);
  // Repeated proved properties may finish from the certified invariant without
  // creating or retiring additional guarded SAT contexts.
  for (size_t iteration = 0; iteration < 20; ++iteration) {
    repeatedSubsetResults.push_back(engine.run(3, safeProperty));
  }
  // The combined parent's other child stays inside the same exact family
  // surface and reuses its learned clauses without touching the independent
  // entry.
  const auto siblingProved = engine.run(3, otherSafeProperty);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(combinedProved.status, PDRStatus::Equivalent);
  EXPECT_EQ(subsetProved.status, PDRStatus::Equivalent);
  EXPECT_EQ(independentProved.status, PDRStatus::Equivalent);
  EXPECT_EQ(different.status, PDRStatus::Different);
  EXPECT_EQ(different.bound, 2u);
  for (const PDRResult& repeatedResult : repeatedSubsetResults) {
    EXPECT_EQ(repeatedResult.status, PDRStatus::Equivalent);
  }
  EXPECT_EQ(siblingProved.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("shared predecessor context activated run=1 level=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared predecessor context activated run=3 level=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared predecessor context activated run=4 level=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "shared predecessor solver pool selected level=1 run=1 entry=0 "
          "cache_hit=0 evicted=0 family_symbols=4"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "shared predecessor solver pool selected level=1 run=3 entry=1 "
          "cache_hit=0 evicted=0 family_symbols=2"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "shared predecessor solver pool selected level=1 run=4 entry=0 "
          "cache_hit=1 evicted=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("restarted=1"),
      std::string::npos)
      << stderrOutput;
  const std::string createdLevelOne =
      "predecessor cached solver created level=1";
  const size_t firstCreation = stderrOutput.find(createdLevelOne);
  ASSERT_NE(firstCreation, std::string::npos) << stderrOutput;
  const size_t secondCreation = stderrOutput.find(
      createdLevelOne, firstCreation + createdLevelOne.size());
  ASSERT_NE(secondCreation, std::string::npos) << stderrOutput;
  EXPECT_EQ(stderrOutput.find(
                createdLevelOne, secondCreation + createdLevelOne.size()),
            std::string::npos)
      << stderrOutput;

  KInductionProblem guardedProblem = problem;
  guardedProblem.usesStrictDualRailEqualityProperty = false;
  auto guardedExactInitCache = std::make_shared<PDRExactInitCache>(
      guardedProblem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  PDREngine guardedEngine(
      guardedProblem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*maxPredecessorQueries=*/0,
      guardedExactInitCache);
  testing::internal::CaptureStderr();
  const auto guardedCombined = guardedEngine.run(3, combinedSafeProperty);
  const auto guardedSubset = guardedEngine.run(3, safeProperty);
  const auto guardedSibling = guardedEngine.run(3, otherSafeProperty);
  std::vector<PDRResult> guardedRepeatedResults;
  guardedRepeatedResults.reserve(16);
  // Repeated guarded properties may finish from the certified invariant
  // without creating or retiring additional guarded SAT contexts.
  for (size_t iteration = 0; iteration < 16; ++iteration) {
    guardedRepeatedResults.push_back(guardedEngine.run(3, safeProperty));
  }
  const std::string guardedStderr =
      testing::internal::GetCapturedStderr();
  EXPECT_EQ(guardedCombined.status, PDRStatus::Equivalent);
  EXPECT_EQ(guardedSubset.status, PDRStatus::Equivalent);
  EXPECT_EQ(guardedSibling.status, PDRStatus::Equivalent);
  for (const PDRResult& repeatedResult : guardedRepeatedResults) {
    EXPECT_EQ(repeatedResult.status, PDRStatus::Equivalent);
  }
  EXPECT_NE(
      guardedStderr.find("reusable invariant clauses certified"),
      std::string::npos)
      << guardedStderr;
  EXPECT_NE(
      guardedStderr.find("reusable invariant clauses injected"),
      std::string::npos)
      << guardedStderr;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFindsCounterexampleFromExactRelationalBootstrapState) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t reset = 4;
  constexpr size_t data = 5;
  problem.state0Symbols = {x, y};
  problem.inputSymbols = {reset, data};
  problem.allSymbols = {x, y, reset, data};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  // This complete but lossy summary excludes the concrete reset state that
  // reaches the mismatch. Neither F[0] nor its cube checks may rely on it.
  problem.bootstrapStateAssignments = {{x, true}, {y, true}};

  BoolExpr* resetAndData =
      BoolExpr::And(BoolExpr::Var(reset), BoolExpr::Var(data));
  BoolExpr* resetInactive = BoolExpr::Not(BoolExpr::Var(reset));
  problem.transitions0.emplace_back(
      x,
      BoolExpr::Or(
          resetAndData,
          BoolExpr::And(
              resetInactive, BoolExpr::Not(BoolExpr::Var(y)))));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::Or(
          resetAndData,
          BoolExpr::And(resetInactive, BoolExpr::Var(x))));
  problem.bad =
      BoolExpr::And(BoolExpr::Var(x), BoolExpr::Not(BoolExpr::Var(y)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBootstrapPrefixEnforcesHistoricalDualRailValidity) {
  KInductionProblem problem;
  constexpr size_t target = 2;
  constexpr size_t mayBeOne = 3;
  constexpr size_t mayBeZero = 4;
  constexpr size_t reset = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {target, mayBeOne, mayBeZero};
  problem.inputSymbols = {reset};
  problem.allSymbols = {target, mayBeOne, mayBeZero, reset};
  problem.totalStateCount = 3;
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{mayBeOne, true}, {mayBeZero, true}};
  problem.dualRailStatePairs = {
      DualRailSymbolPair{mayBeOne, mayBeZero}};

  // The target can become true only from the invalid rail encoding 00. The
  // exact reset prefix must exclude that historical pseudo-state just as every
  // ordinary PDR frame does.
  problem.transitions0 = {
      {target,
       BoolExpr::Not(BoolExpr::Or(
           BoolExpr::Var(mayBeOne), BoolExpr::Var(mayBeZero)))},
      {mayBeOne, BoolExpr::createTrue()},
      {mayBeZero, BoolExpr::createTrue()}};
  problem.bad = BoolExpr::Var(target);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       LazyTransitionSupportCacheIsSharedAcrossResolversWithoutDagRemap) {
  KInductionProblem problem;
  constexpr size_t combinedState = 10;
  constexpr size_t combinedInput = 11;
  constexpr size_t localState = 2;
  constexpr size_t localInput = 3;
  BoolExpr* localNext =
      BoolExpr::And(BoolExpr::Var(localState), BoolExpr::Var(localInput));

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->localToCombinedByDesign[0].emplace(localState, combinedState);
  lazyTransitions->localToCombinedByDesign[0].emplace(localInput, combinedInput);
  lazyTransitions->sourceByStateSymbol.emplace(
      combinedState, LazyTransitionSource{0, localNext});
  problem.lazyTransitions = lazyTransitions;
  problem.state0Symbols = {combinedState};
  problem.inputSymbols = {combinedInput};
  problem.allSymbols = {combinedState, combinedInput};

  {
    const TransitionExprResolver transitionByState(problem);
    const auto& support = transitionByState.support(combinedState);
    EXPECT_EQ(support, (std::set<size_t>{combinedState, combinedInput}));
    EXPECT_EQ(transitionByState.nodeCount(combinedState), 3u);
  }

  ASSERT_NE(
      lazyTransitions->supportByStateSymbol.find(combinedState),
      lazyTransitions->supportByStateSymbol.end());
  ASSERT_NE(
      lazyTransitions->nodeCountByStateSymbol.find(combinedState),
      lazyTransitions->nodeCountByStateSymbol.end());
  // Support and node-count queries must not force a lazy BoolExpr remap. In
  // BlackParrot PDR those queries happen across many output batches; sharing
  // this metadata avoids repeatedly
  // walking the same source DAGs before any transition needs SAT encoding.
  EXPECT_TRUE(lazyTransitions->remappedByStateSymbol.empty());

  const TransitionExprResolver secondTransitionByState(problem);
  EXPECT_EQ(
      secondTransitionByState.support(combinedState),
      (std::set<size_t>{combinedState, combinedInput}));
  EXPECT_EQ(secondTransitionByState.nodeCount(combinedState), 3u);
  EXPECT_TRUE(lazyTransitions->remappedByStateSymbol.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       LazyTransitionUnpublishedSupportStaysDesignPrivate) {
  KInductionProblem problem;
  constexpr size_t combinedState0 = 10;
  constexpr size_t combinedState1 = 20;
  constexpr size_t localState = 2;
  constexpr size_t unpublishedLocal = 42;
  BoolExpr* localNext =
      BoolExpr::Xor(BoolExpr::Var(localState), BoolExpr::Var(unpublishedLocal));

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->localToCombinedByDesign[0].emplace(localState, combinedState0);
  lazyTransitions->localToCombinedByDesign[1].emplace(localState, combinedState1);
  lazyTransitions->sourceByStateSymbol.emplace(
      combinedState0, LazyTransitionSource{0, localNext});
  lazyTransitions->sourceByStateSymbol.emplace(
      combinedState1, LazyTransitionSource{1, localNext});
  problem.lazyTransitions = lazyTransitions;
  problem.state0Symbols = {combinedState0};
  problem.state1Symbols = {combinedState1};
  problem.allSymbols = {combinedState0, combinedState1};

  const TransitionExprResolver transitionByState(problem);
  const size_t private0 = makePrivateProofLeafSymbol(0, unpublishedLocal);
  const size_t private1 = makePrivateProofLeafSymbol(1, unpublishedLocal);

  EXPECT_NE(private0, private1);
  EXPECT_EQ(
      transitionByState.support(combinedState0),
      (std::set<size_t>{combinedState0, private0}));
  EXPECT_EQ(
      transitionByState.support(combinedState1),
      (std::set<size_t>{combinedState1, private1}));
  EXPECT_EQ(
      lazyTransitions->localToCombinedByDesign[0].at(unpublishedLocal),
      private0);
  EXPECT_EQ(
      lazyTransitions->localToCombinedByDesign[1].at(unpublishedLocal),
      private1);

  EXPECT_EQ(
      transitionByState.at(combinedState0)->getSupportVars(),
      (std::set<size_t>{combinedState0, private0}));
  EXPECT_EQ(
      transitionByState.at(combinedState1)->getSupportVars(),
      (std::set<size_t>{combinedState1, private1}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       LazyDualRailTransitionSupportUsesBothRailsWithoutDagRemap) {
  KInductionProblem problem;
  constexpr size_t railOne = 10;
  constexpr size_t railZero = 11;
  constexpr size_t combinedInput = 12;
  constexpr size_t localState = 2;
  constexpr size_t localInput = 3;
  BoolExpr* localNext =
      BoolExpr::Xor(BoolExpr::Var(localState), BoolExpr::Var(localInput));

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->dualRailStateByLocalSymbolByDesign[0].emplace(
      localState, DualRailSymbolPair{railOne, railZero});
  lazyTransitions->localToCombinedByDesign[0].emplace(localInput, combinedInput);
  lazyTransitions->sourceByStateSymbol.emplace(
      railOne, LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailOne});
  problem.lazyTransitions = lazyTransitions;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {railOne, railZero};
  problem.inputSymbols = {combinedInput};
  problem.allSymbols = {railOne, railZero, combinedInput};

  const TransitionExprResolver transitionByState(problem);
  EXPECT_EQ(
      transitionByState.support(railOne),
      (std::set<size_t>{railOne, railZero, combinedInput}));
  // A support-only query must stay in the lazy source-expression layer; the
  // lifted dual-rail BoolExpr is materialized later only if SAT encoding needs
  // this transition.
  EXPECT_TRUE(lazyTransitions->remappedByStateSymbol.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       LazyDualRailMaterializationUsesRailsInsteadOfBinaryPrivateLeaves) {
  KInductionProblem problem;
  constexpr size_t railOne = 10;
  constexpr size_t railZero = 11;
  constexpr size_t combinedInput = 12;
  constexpr size_t localState = 2;
  constexpr size_t localInput = 3;
  BoolExpr* localNext =
      BoolExpr::Xor(BoolExpr::Var(localState), BoolExpr::Var(localInput));

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->dualRailStateByLocalSymbolByDesign[0].emplace(
      localState, DualRailSymbolPair{railOne, railZero});
  lazyTransitions->localToCombinedByDesign[0].emplace(localInput, combinedInput);
  lazyTransitions->sourceByStateSymbol.emplace(
      railOne, LazyTransitionSource{0, localNext, LazyTransitionRail::DualRailOne});
  problem.lazyTransitions = lazyTransitions;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {railOne, railZero};
  problem.inputSymbols = {combinedInput};
  problem.allSymbols = {railOne, railZero, combinedInput};

  const TransitionExprResolver transitionByState(problem);

  EXPECT_EQ(
      transitionByState.at(railOne)->getSupportVars(),
      (std::set<size_t>{railOne, railZero, combinedInput}));
  EXPECT_NE(
      lazyTransitions->remappedByStateSymbol.find(railOne),
      lazyTransitions->remappedByStateSymbol.end());
}

TEST_F(SequentialEquivalenceStrategyTests,
       LazyDualRailWideTargetSupportIsCollectedAsOneUnion) {
  KInductionProblem problem;
  constexpr size_t railOne = 10;
  constexpr size_t railZero = 11;
  constexpr size_t combinedInput = 12;
  constexpr size_t localState = 2;
  constexpr size_t localInput = 3;
  constexpr size_t targetBase = 100;
  constexpr size_t targetCount = 20;
  BoolExpr* localNext =
      BoolExpr::Xor(BoolExpr::Var(localState), BoolExpr::Var(localInput));

  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->dualRailStateByLocalSymbolByDesign[0].emplace(
      localState, DualRailSymbolPair{railOne, railZero});
  lazyTransitions->localToCombinedByDesign[0].emplace(localInput, combinedInput);
  problem.state0Symbols = {railOne, railZero};
  std::vector<size_t> targets;
  targets.reserve(targetCount);
  for (size_t index = 0; index < targetCount; ++index) {
    const size_t target = targetBase + index;
    const auto rail =
        (index % 2 == 0) ? LazyTransitionRail::DualRailOne
                         : LazyTransitionRail::DualRailZero;
    targets.push_back(target);
    problem.state0Symbols.push_back(target);
    lazyTransitions->sourceByStateSymbol.emplace(
        target, LazyTransitionSource{0, localNext, rail});
  }
  problem.lazyTransitions = lazyTransitions;
  problem.usesDualRailStateEncoding = true;
  problem.inputSymbols = {combinedInput};
  problem.allSymbols = problem.state0Symbols;
  problem.allSymbols.push_back(combinedInput);

  const TransitionExprResolver transitionByState(problem);
  const std::unordered_set<size_t> knownStateSymbols(
      problem.state0Symbols.begin(), problem.state0Symbols.end());
  std::unordered_set<size_t> stateSupport;
  std::unordered_set<size_t> allSupport;

  transitionByState.collectSupportForTargets(
      targets, knownStateSymbols, stateSupport, allSupport);

  EXPECT_EQ(stateSupport, (std::unordered_set<size_t>{railOne, railZero}));
  EXPECT_EQ(
      allSupport,
      (std::unordered_set<size_t>{railOne, railZero, combinedInput}));
  // This regression covers the dynamic-node dual-rail wall: wide lazy
  // dual-rail COIs should walk the shared source expression once as a union,
  // not populate one cached per-target support set before SAT even starts.
  EXPECT_TRUE(lazyTransitions->supportByStateSymbol.empty());
  EXPECT_TRUE(lazyTransitions->remappedByStateSymbol.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineConstrainsResetInputsOnFirstPostBootstrapStep) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.inputSymbols = {3, 4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{3, true}, {4, false}};
  // The concrete reset prefix drives x=0, then deasserts the reset controls as
  // r=0 and g=1.  The abstract PDR F[0] summary contains only x=0, so a
  // level-0 predecessor query that forgets reset-input deassertion can invent
  // r=1,g=1 on the first normal step and falsely reach x'=1.
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(
      2, BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4)));
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  EXPECT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1)
          .has_value());

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineConstrainsResetInputsInPostBootstrapBadQueries) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.inputSymbols = {3, 4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{3, true}, {4, false}};
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  // The bad predicate is input-only. PDR must still apply the post-reset
  // deasserted reset controls before deciding whether this is a real bad frame.
  problem.bad = BoolExpr::And(BoolExpr::Var(3), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  EXPECT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}






TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailBadCubeSkipsUnchangedFrameClauseSync) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  problem.totalStateCount = 1;
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(state));
  problem.initialStateAssignments = {{state, false}};
  problem.initializedStateCount = 1;
  problem.transitions0.emplace_back(state, BoolExpr::Var(state));
  problem.observedOutputExprs0 = {
      BoolExpr::Var(state),
      BoolExpr::Not(BoolExpr::Var(state))};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(),
      BoolExpr::createFalse()};
  problem.observedOutputNames = {"state_is_one", "state_is_zero"};
  problem.bad = BoolExpr::Or(
      BoolExpr::Var(state), BoolExpr::Not(BoolExpr::Var(state)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Two output-bad formulas share the same frame and symbol surface. The
  // cached bad-cube solver should not rescan already synchronized frame clauses
  // before asking the second formula.
  EXPECT_NE(
      stderrOutput.find("bad cube cached frame clauses unchanged"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailHugeStateSurfaceCachesNarrowBadQueries) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  constexpr size_t invariantState = 3;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state, invariantState};
  problem.allSymbols = {state, invariantState};
  // Ariane has a multi-million-bit rail surface. Model only the cache-policy
  // signal here so the unit test stays tiny while still protecting that shape.
  problem.totalStateCount = 300000;
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(state)),
      BoolExpr::Not(BoolExpr::Var(invariantState)));
  problem.initialStateAssignments = {
      {state, false}, {invariantState, false}};
  problem.initializedStateCount = 2;
  problem.transitions0.emplace_back(state, BoolExpr::Var(state));
  problem.transitions0.emplace_back(invariantState, BoolExpr::createFalse());
  problem.observedOutputExprs0 = {
      BoolExpr::Var(state), BoolExpr::Var(state)};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.observedOutputNames = {
      "huge_state_bad_cube_cache_0", "huge_state_bad_cube_cache_1"};
  problem.originalObservedOutputCount = 278;
  problem.bad = BoolExpr::Var(state);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = BoolExpr::Not(BoolExpr::Var(invariantState));
  problem.inductionBad = BoolExpr::Var(invariantState);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  testing::internal::CaptureStderr();
  PDREngine firstEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0,
                        exactInitCache);
  const auto firstResult = firstEngine.run(1);
  PDREngine secondEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0,
                         exactInitCache);
  const auto secondResult = secondEngine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(firstResult.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_EQ(secondResult.status, firstResult.status) << stderrOutput;
  EXPECT_EQ(secondResult.bound, firstResult.bound) << stderrOutput;
  // F[0] is immutable across output batches. Bad-state queries must reuse the
  // exact predecessor SAT instance instead of retaining a duplicate Init CNF.
  EXPECT_NE(
      stderrOutput.find("shared exact F[0] solver used for bad cube"),
      std::string::npos)
      << stderrOutput;
  // Higher-frame caching is bounded by the exact SAT surface, not unrelated
  // state elsewhere in the design. The narrow query retains its solver while
  // the newly learned IC3 clause is streamed into that exact frame context.
  EXPECT_NE(
      stderrOutput.find("bad cube cached frame clauses added=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("shared exact F[0] symbols reused for bad cube"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("formula closed support reused"),
            std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("bad cube cached solver disabled"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("predecessor cached solver created level=1"),
            std::string::npos)
      << stderrOutput;
  // Frame closure caching is exact and applies to this non-local two-output
  // surface as well as the historical one-output residual-leaf case.
  EXPECT_NE(
      stderrOutput.find("predecessor frame symbol cache built level=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("predecessor target surface cached"),
            std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("predecessor target surface reused"),
            std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("predecessor cached solver disabled"),
            std::string::npos)
      << stderrOutput;
}












TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailKeepsSelectedSolverForBadCubeQueries) {
  KInductionProblem problem;
  std::vector<size_t> symbols;
  constexpr size_t kSmallStateCount = 6;
  symbols.reserve(kSmallStateCount);

  for (size_t i = 0; i < kSmallStateCount; ++i) {
    const size_t symbol = i + 2;
    symbols.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }

  problem.usesDualRailStateEncoding = true;
  problem.initialCondition = BoolExpr::createTrue();
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::simplify(makeOrChain(symbols));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar proofConflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_BAD_CUBE_CONFLICT_LIMIT", "5000");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  // Bad-cube queries are core PDR obligations. They should stay on the selected
  // engine solver rather than silently switching to a separate backend.
  EXPECT_EQ(
      stderrOutput.find("SEC PDR stats: bad cube query solver=cadical"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailBadCubeBudgetReturnsInconclusive) {
  KInductionProblem problem;
  constexpr size_t stateA = 2;
  constexpr size_t stateB = 3;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {stateA, stateB};
  problem.allSymbols = {stateA, stateB};
  problem.totalStateCount = 2;
  problem.transitions0 = {
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};

  // This contradiction is not unit-propagation trivial in CNF form.  With a
  // one-conflict local budget, the bad-cube query must return UNKNOWN and make
  // PDR inconclusive rather than treating UNKNOWN as "no bad cube".
  BoolExpr* bad = BoolExpr::And(
      BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB)),
      BoolExpr::And(
          BoolExpr::Or(BoolExpr::Not(BoolExpr::Var(stateA)), BoolExpr::Var(stateB)),
          BoolExpr::And(
              BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Not(BoolExpr::Var(stateB))),
              BoolExpr::Or(
                  BoolExpr::Not(BoolExpr::Var(stateA)),
                  BoolExpr::Not(BoolExpr::Var(stateB))))));
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar conflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_BAD_CUBE_CONFLICT_LIMIT", "1");
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, PDRStatus::Inconclusive);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPredecessorZeroDecisionLimitIsUnbounded) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB};
  problem.allSymbols = {targetState, stateA, stateB};
  problem.totalStateCount = 3;

  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"single_output_unbounded"};
  problem.originalObservedOutputCount = 1266;

  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different) << stderrOutput;
  // Zero means unlimited. Section V's incremental query must answer directly,
  // rather than turning this reachable predecessor into UNKNOWN.
  EXPECT_NE(
      stderrOutput.find("result=sat cached_assumptions=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor query budget exhausted"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPredecessorResourceLimitReturnsInconclusive) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB};
  problem.allSymbols = {targetState, stateA, stateB};
  problem.totalStateCount = 3;

  // This four-clause formula excludes every assignment without simplifying to
  // a unit contradiction. It forces the predecessor SAT query to consume its
  // deliberately tiny local budget.
  BoolExpr* hardFalse = BoolExpr::And(
      BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB)),
      BoolExpr::And(
          BoolExpr::Or(BoolExpr::Not(BoolExpr::Var(stateA)),
                       BoolExpr::Var(stateB)),
          BoolExpr::And(
              BoolExpr::Or(BoolExpr::Var(stateA),
                           BoolExpr::Not(BoolExpr::Var(stateB))),
              BoolExpr::Or(BoolExpr::Not(BoolExpr::Var(stateA)),
                           BoolExpr::Not(BoolExpr::Var(stateB))))));
  problem.transitions0 = {
      {targetState, hardFalse},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"single_output_limited"};
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar conflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_CONFLICT_LIMIT", "1");
  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "1");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  // UNKNOWN is not UNSAT. The local limit skips this output as inconclusive;
  // it must not retry the same proof obligation through another solver path.
  EXPECT_NE(
      stderrOutput.find("predecessor query budget exhausted"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("cached_assumptions=1"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("target_cube="), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("observed_outputs=1"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("residual_budget="), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPropagationResourceLimitContinuesCurrentRun) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t gateState = 3;
  constexpr size_t firstRailState = 10;
  constexpr size_t railPairCount = 16;

  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, gateState};
  problem.allSymbols = problem.state0Symbols;
  problem.transitions0 = {
      {targetState, BoolExpr::Var(gateState)},
      {gateState, BoolExpr::Var(gateState)}};
  for (size_t pair = 0; pair < railPairCount; ++pair) {
    const size_t mayBeOne = firstRailState + pair * 2;
    const size_t mayBeZero = mayBeOne + 1;
    problem.state0Symbols.push_back(mayBeOne);
    problem.state0Symbols.push_back(mayBeZero);
    problem.allSymbols.push_back(mayBeOne);
    problem.allSymbols.push_back(mayBeZero);
    problem.transitions0.emplace_back(mayBeOne, BoolExpr::Var(mayBeOne));
    problem.transitions0.emplace_back(mayBeZero, BoolExpr::Var(mayBeZero));
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{mayBeOne, mayBeZero});
  }
  problem.totalStateCount = problem.state0Symbols.size();
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(targetState)),
      BoolExpr::Not(BoolExpr::Var(gateState)));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"propagation_budget"};

  // F[0] makes the first blocker immediate. Propagating that blocker from
  // F[1] has a satisfiable predecessor and deliberately exceeds this tiny
  // decision budget. Figure 9 leaves the clause in F[1] and continues.
  const ScopedEnvVar conflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_CONFLICT_LIMIT", "0");
  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "1");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("propagation left clause in frame"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("purpose=propagate"), std::string::npos)
      << stderrOutput;
  // Reaching the normal frame bound proves UNKNOWN did not abort the run.
  EXPECT_NE(
      stderrOutput.find("max frame budget exhausted"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailAesSizedSatLeafUsesIncrementalPredecessor) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB};
  problem.allSymbols = {targetState, stateA, stateB};
  problem.totalStateCount = 3;

  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"single_output_aes_sized_incremental"};
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different) << stderrOutput;
  // Section V reuses the exact incremental SAT model as the predecessor. A SAT
  // answer must not rebuild the same complete query in a fresh solver.
  EXPECT_NE(
      stderrOutput.find("result=sat cached_assumptions=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor_cube="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("cached_query_us="),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("fallback=exact"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailBlockingQueryUsesFiniteRoleBudget) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB};
  problem.allSymbols = {targetState, stateA, stateB};
  problem.totalStateCount = 3;

  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"single_output_blocking_budget"};
  problem.originalObservedOutputCount = 1266;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Mandatory Figure 6 blocking receives one finite role-based budget. The
  // enclosing design's original output count does not select this policy.
  EXPECT_NE(
      stderrOutput.find("conflict_limit=250000"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("decision_limit=10000000"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("purpose=block"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("residual_budget="), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineExplicitBatchProbeLimitsOverrideFullRoleBudget) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB};
  problem.allSymbols = problem.state0Symbols;
  problem.totalStateCount = problem.state0Symbols.size();
  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"explicit_batch_probe_budget"};

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(
      1, problem.property,
      PDRQueryLimits{/*predecessorConflictLimit=*/10000,
                     /*predecessorDecisionLimit=*/150000});
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(stderrOutput.find("conflict_limit=10000"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("decision_limit=150000"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("purpose=block"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("conflict_limit=250000"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPredecessorUsesExactOnDemandStateSurface) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  constexpr size_t decoyState = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB, decoyState};
  problem.allSymbols = {targetState, stateA, stateB, decoyState};
  problem.totalStateCount = 4;

  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)},
      {decoyState, BoolExpr::Var(decoyState)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"single_output_on_demand_cache"};
  problem.originalObservedOutputCount = 1266;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Section V encodes the exact frame and target-transition cone on demand.
  // The decoy is absent from both, so existentially omitting it changes no SAT
  // answer and avoids carrying it through every incremental query.
  EXPECT_NE(
      stderrOutput.find(
          "predecessor_symbols=2 solver_symbols=3 cached_solver_symbols=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "predecessor cached solver created level=0 symbols=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("solver_symbols=4 cached_solver_symbols=4"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor frame symbol cache built level=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor transition encoder cached"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor target surface cached"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineExtendsIncrementalPredecessorSolverForWiderCone) {
  KInductionProblem problem;
  constexpr size_t targetA = 2;
  constexpr size_t targetB = 3;
  constexpr size_t sourceA = 4;
  constexpr size_t sourceB = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetA, targetB, sourceA, sourceB};
  problem.allSymbols = problem.state0Symbols;
  problem.totalStateCount = problem.state0Symbols.size();
  problem.initialCondition = BoolExpr::createTrue();
  problem.initialStateAssignments = {
      {targetA, false},
      {targetB, false},
      {sourceA, false},
      {sourceB, false}};
  problem.transitions0 = {
      {targetA, BoolExpr::Var(sourceA)},
      {targetB, BoolExpr::Var(sourceB)},
      {sourceA, BoolExpr::createFalse()},
      {sourceB, BoolExpr::createFalse()}};
  problem.bad = BoolExpr::Or(
      BoolExpr::Var(targetA), BoolExpr::Var(targetB));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  const std::string frameZeroCreated =
      "predecessor cached solver created level=0";
  const size_t firstFrameZeroBuild = stderrOutput.find(frameZeroCreated);
  ASSERT_NE(firstFrameZeroBuild, std::string::npos) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          frameZeroCreated,
          firstFrameZeroBuild + frameZeroCreated.size()),
      std::string::npos)
      << stderrOutput;
  // A wider cone at the same PDR level extends the one incremental SAT
  // instance, while each distinct frame keeps its own exact solver context.
  EXPECT_NE(
      stderrOutput.find("predecessor cached solver surface extended"),
      std::string::npos)
      << stderrOutput;
  // Relation clauses already present in the incremental solver are retained;
  // extending the exact query adds only relations exposed by the new symbols.
  EXPECT_NE(
      stderrOutput.find(
          "predecessor relation surface extended added_symbols="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "predecessor transition encoder cache extended encoders="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "predecessor cached solver frame sync source=full_frame"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailHugeStateSurfaceCachesExactFrameZeroPredecessor) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  constexpr size_t decoyState = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB, decoyState};
  problem.allSymbols = {targetState, stateA, stateB, decoyState};
  // Model Ariane's multi-million-bit rail surface without allocating it in the
  // unit test. Exact F[0] and T are shared across output batches even at this
  // width; later learned PDR frames still retain their bounded cache policy.
  problem.totalStateCount = 300000;

  problem.transitions0 = {
      {targetState, BoolExpr::Or(BoolExpr::Var(stateA), BoolExpr::Var(stateB))},
      {stateA, BoolExpr::Var(stateA)},
      {stateB, BoolExpr::Var(stateB)},
      {decoyState, BoolExpr::Var(decoyState)}};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {BoolExpr::Var(targetState)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"huge_state_preparation_surface"};
  problem.originalObservedOutputCount = 278;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  auto exactInitCache = std::make_shared<PDRExactInitCache>(
      problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  testing::internal::CaptureStderr();
  PDREngine firstEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0,
                        exactInitCache);
  const auto firstResult = firstEngine.run(1);
  PDREngine secondEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0,
                         exactInitCache);
  const auto secondResult = secondEngine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(firstResult.status, secondResult.status);
  EXPECT_EQ(firstResult.bound, secondResult.bound);
  EXPECT_NE(stderrOutput.find("predecessor target surface cached"),
            std::string::npos)
      << stderrOutput;
  const std::string frameZeroCreated =
      "predecessor cached solver created level=0";
  const size_t firstFrameZeroBuild = stderrOutput.find(frameZeroCreated);
  ASSERT_NE(firstFrameZeroBuild, std::string::npos) << stderrOutput;
  EXPECT_EQ(stderrOutput.find(
                frameZeroCreated,
                firstFrameZeroBuild + frameZeroCreated.size()),
            std::string::npos)
      << stderrOutput;
  // The second output batch asks the same immutable level-zero predecessor
  // question. It should reuse the exact answer before rebuilding any target
  // preparation or invoking the already shared SAT solver.
  EXPECT_NE(
      stderrOutput.find(
          "predecessor result cache hit level=0 has_predecessor=1 "
          "has_model=1 shared_f0=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(
                "shared exact F[0] predecessor symbols reused"),
            std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("predecessor cached solver disabled"),
            std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPredecessorEncodingBudgetReturnsInconclusive) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  std::vector<size_t> sourceStates;
  sourceStates.reserve(12);
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols.push_back(targetState);
  problem.allSymbols.push_back(targetState);
  for (size_t offset = 0; offset < 12; ++offset) {
    const size_t symbol = 3 + offset;
    sourceStates.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  problem.totalStateCount = problem.state0Symbols.size();
  problem.transitions0.emplace_back(targetState, makeOrChain(sourceStates));
  // This proof would normally build the target transition cone before solving
  // the predecessor query. Keep the artificial encoding limit tiny so the test
  // verifies the local inconclusive exit before expensive clause generation.
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar nodeLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_ENCODING_NODE_LIMIT", "4");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor encoding budget exhausted"),
      std::string::npos)
      << stderrOutput;
}









TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineReturnsInconclusiveWhenZeroBudgetNeedsFrames) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);

  EXPECT_EQ(result.status, PDRStatus::Inconclusive);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractDoesNotTreatResetAsInitialState) {
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
  EXPECT_FALSE(model.stateBits.empty());
  // Reset controls the transition relation. Without a declared initializer it
  // must not constrain IC3's exact initial frame.
  EXPECT_TRUE(model.initialStateValueByKey.empty());
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
       SequentialDesignModelExtractClassifiesNegativeEdgePrimitiveClock) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createPosToNegSameDomainTop(library, "top");

  const auto extracted = SequentialDesignModel::extract(top);

  const auto posKey = findKeyByDisplayName(extracted, "ff_pos.Q[0]");
  const auto negKey = findKeyByDisplayName(extracted, "ff_neg.Q[0]");
  ASSERT_NE(
      extracted.clockEventByStateKey.find(posKey),
      extracted.clockEventByStateKey.end());
  ASSERT_NE(
      extracted.clockEventByStateKey.find(negKey),
      extracted.clockEventByStateKey.end());
  EXPECT_EQ(extracted.clockEventByStateKey.at(posKey).phase, ClockPhase::Pos);
  EXPECT_EQ(extracted.clockEventByStateKey.at(negKey).phase, ClockPhase::Neg);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractComposesPosedgeBeforeNegedgeSameDomain) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createPosToNegSameDomainTop(library, "top");

  const auto extracted = SequentialDesignModel::extract(top);

  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto posKey = findKeyByDisplayName(extracted, "ff_pos.Q[0]");
  const auto negKey = findKeyByDisplayName(extracted, "ff_neg.Q[0]");
  auto* negNext = extracted.nextStateExprByStateKey.at(negKey);

  EXPECT_TRUE(negNext->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(posKey), false},
       {extracted.inputVarByKey.at(negKey), false}}));
  EXPECT_FALSE(negNext->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(posKey), true},
       {extracted.inputVarByKey.at(negKey), true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractClassifiesInvertedClockTreePhase) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top = createInvertedClockDffTop(library, "top", invModel);

  const auto extracted = SequentialDesignModel::extract(top);

  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  ASSERT_NE(
      extracted.clockEventByStateKey.find(stateKey),
      extracted.clockEventByStateKey.end());
  EXPECT_EQ(extracted.clockEventByStateKey.at(stateKey).phase, ClockPhase::Neg);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractClassifiesNamedClockTreeInverterPhase) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top = createInvertedClockDffTop(
      library, "top", invModel, "clkinv_leaf_0");

  const auto extracted = SequentialDesignModel::extract(top);

  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  ASSERT_NE(
      extracted.clockEventByStateKey.find(stateKey),
      extracted.clockEventByStateKey.end());
  EXPECT_EQ(extracted.clockEventByStateKey.at(stateKey).phase, ClockPhase::Neg);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSkipsOutputConeSpanningClockDomains) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* andModel = createAnd2Model(primitives);
  auto* top = createMultiClockDomainOutputTop(library, "top", andModel);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.totalObservedOutputCount(), 1u);
  EXPECT_EQ(extracted.coveredObservedOutputCount(), 0u);
  ASSERT_EQ(extracted.skippedObservedOutputs.size(), 1u);
  const auto outputKey = findKeyByDisplayName(extracted, "out[0]");
  const auto skipIt = extracted.connectivitySkipInfoByKey.find(outputKey);
  ASSERT_NE(skipIt, extracted.connectivitySkipInfoByKey.end());
  EXPECT_EQ(skipIt->second.origin, ConnectivitySkipOrigin::MultiClockDomain);

  auto strategy = makeBinaryExtractedSecStrategy();
  const auto result = strategy.runExtractedModels(extracted, extracted, 1);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, 1u);
  ASSERT_EQ(result.multiClockDomainSkippedOutputs.size(), 1u);
  EXPECT_NE(
      result.multiClockDomainSkippedOutputs.front().find("multi-clock-domain"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractModelsStructuredMemoryWithoutBoundaryFallback) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createSinglePortMemoryModel(primitives, "MEM1P");
  auto* top = createSinglePortMemoryTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  auto hasBoundaryRoleName = [&](const std::vector<SignalKey>& keys,
                                 const std::string& prefix) {
    return std::any_of(keys.begin(), keys.end(), [&](const SignalKey& key) {
      const auto it = extracted.displayNameByKey.find(key);
      return it != extracted.displayNameByKey.end() &&
             it->second.rfind(prefix, 0) == 0;
    });
  };

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryInputKeys, "mem0."));
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryOutputKeys, "mem0."));
  EXPECT_FALSE(extracted.stateBits.empty());
  EXPECT_FALSE(extracted.nextStateExprByStateKey.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractModelsImportedLibertyMemoryWithoutOpaqueBoundaryTerms) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = loadLibertyMemoryModel(
      primitives, "fakeram45_64x32.lib", "fakeram45_64x32");
  auto* top = createMirroredInstanceTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  auto hasBoundaryRoleName = [&](const std::vector<SignalKey>& keys,
                                 const std::string& prefix) {
    return std::any_of(keys.begin(), keys.end(), [&](const SignalKey& key) {
      const auto it = extracted.displayNameByKey.find(key);
      return it != extracted.displayNameByKey.end() &&
             it->second.rfind(prefix, 0) == 0;
    });
  };

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryInputKeys, "mem0."));
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryOutputKeys, "mem0."));
  EXPECT_FALSE(extracted.stateBits.empty());
  EXPECT_FALSE(extracted.nextStateExprByStateKey.empty());
}

TEST_F(
    SequentialEquivalenceStrategyTests,
    StructuredMemoryDependencyBatchBuildsRealCva6PerfCounterReadAddressRoot) {
  const auto context = resolveCva6SourceContextForSecTests();
  if (!context.has_value()) {
    GTEST_SKIP()
        << "CVA6 source tree is not available for the real-source SEC memory regression test";
  }

  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadRealCva6PerfCountersTargetConfigTopForSecTests(
      library,
      *context,
      "real_cva6_perf_counters_module_with_target_config_sec_memory_probe");

  detail::ScopedDnlContextForTest dnlContext(top);
  auto* dnl = dnlContext.dnl();
  ASSERT_NE(dnl, nullptr);
  auto requestedTermID = detail::findTermByDisplayNameForTest(
      dnl, "dut.generic_counter_q_mem.RADDR[18]");
  if (!requestedTermID.has_value()) {
    // Different CVA6 target configs can infer a narrower generated memory
    // address port. The dependency regression only needs a real RADDR bit, not
    // a specific width-dependent index.
    requestedTermID = detail::findFirstTermByDisplayPrefixForTest(
        dnl, "dut.generic_counter_q_mem.RADDR[");
  }
  ASSERT_TRUE(requestedTermID.has_value());

  const auto probe =
      detail::probeRequestedBuilderOutputForTest(dnl, *requestedTermID);
  std::ostringstream normalizationChain;
  for (size_t index = 0; index < probe.normalizationChain.size(); ++index) {
    if (index != 0) {
      normalizationChain << " -> ";
    }
    normalizationChain << probe.normalizationChain[index];
  }
  std::ostringstream supportTerms;
  for (size_t index = 0; index < probe.rootSupportTerms.size(); ++index) {
    if (index != 0) {
      supportTerms << " | ";
    }
    supportTerms << probe.rootSupportTerms[index];
  }
  std::ostringstream combinationalInputs;
  for (size_t index = 0; index < probe.rootCombinationalInputs.size(); ++index) {
    if (index != 0) {
      combinationalInputs << " | ";
    }
    combinationalInputs << probe.rootCombinationalInputs[index];
  }
  std::ostringstream driverSpine;
  for (size_t index = 0; index < probe.driverSpine.size(); ++index) {
    if (index != 0) {
      driverSpine << " -> ";
    }
    driverSpine << probe.driverSpine[index];
  }

  ASSERT_TRUE(probe.normalizedRoot.has_value())
      << "normalization chain: " << normalizationChain.str();
  EXPECT_TRUE(probe.hasBuiltExpr)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " did not yield a valid clause-builder expression; support terms: "
      << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str();
  EXPECT_FALSE(probe.hasSkip)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " was skipped: " << probe.skipDetail
      << "; support terms: " << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str();
}

TEST_F(
    SequentialEquivalenceStrategyTests,
    StructuredMemoryDependencyBatchBuildsRealCva6TopPerfCounterWriteEnableRoot) {
  const auto context = resolveCva6SourceContextForSecTests();
  if (!context.has_value()) {
    GTEST_SKIP()
        << "CVA6 source tree is not available for the real-source SEC memory regression test";
  }

  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadRealCva6PerfCountersTargetConfigTopForSecTests(
      library,
      *context,
      "real_cva6_perf_counters_module_with_target_config_sec_we_probe");

  detail::ScopedDnlContextForTest dnlContext(top);
  auto* dnl = dnlContext.dnl();
  ASSERT_NE(dnl, nullptr);
  const auto skippedReportPath =
      std::filesystem::current_path() / "skipped_no_driver_pos.txt";
  std::filesystem::remove(skippedReportPath);
  const bool previousReportSkippedPOs =
      KEPLER_FORMAL::Config::getReportSkippedPOs();
  KEPLER_FORMAL::Config::setReportSkippedPOs(true);
  const auto selectedProbe =
      detail::findFirstBuildableOutputProbeByDisplayPrefixForTest(
          dnl,
          "dut.generic_counter_q_mem.WE[");
  KEPLER_FORMAL::Config::setReportSkippedPOs(previousReportSkippedPOs);
  ASSERT_TRUE(selectedProbe.has_value())
      << "No buildable generated CVA6 perf-counter memory WE bit was found";
  const auto& probe = *selectedProbe;
  std::ostringstream normalizationChain;
  for (size_t index = 0; index < probe.normalizationChain.size(); ++index) {
    if (index != 0) {
      normalizationChain << " -> ";
    }
    normalizationChain << probe.normalizationChain[index];
  }
  std::ostringstream supportTerms;
  for (size_t index = 0; index < probe.rootSupportTerms.size(); ++index) {
    if (index != 0) {
      supportTerms << " | ";
    }
    supportTerms << probe.rootSupportTerms[index];
  }
  std::ostringstream combinationalInputs;
  for (size_t index = 0; index < probe.rootCombinationalInputs.size(); ++index) {
    if (index != 0) {
      combinationalInputs << " | ";
    }
    combinationalInputs << probe.rootCombinationalInputs[index];
  }
  std::ostringstream driverSpine;
  for (size_t index = 0; index < probe.driverSpine.size(); ++index) {
    if (index != 0) {
      driverSpine << " -> ";
    }
    driverSpine << probe.driverSpine[index];
  }
  std::string skippedReport;
  if (std::ifstream skippedReportFile(skippedReportPath);
      skippedReportFile.good()) {
    std::ostringstream report;
    report << skippedReportFile.rdbuf();
    skippedReport = report.str();
  }

  ASSERT_TRUE(probe.normalizedRoot.has_value())
      << "normalization chain: " << normalizationChain.str();
  EXPECT_TRUE(probe.hasBuiltExpr)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " did not yield a valid clause-builder expression; support terms: "
      << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str()
      << "; skipped report: " << skippedReport;
  EXPECT_FALSE(probe.hasSkip)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " was skipped: " << probe.skipDetail
      << "; support terms: " << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str()
      << "; skipped report: " << skippedReport;
}

TEST_F(
    SequentialEquivalenceStrategyTests,
    StructuredMemoryDependencyBatchBuildsRealCva6TopPerfCounterWriteDataSliceRoot) {
  const auto context = resolveCva6SourceContextForSecTests();
  if (!context.has_value()) {
    GTEST_SKIP()
        << "CVA6 source tree is not available for the real-source SEC memory regression test";
  }

  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  const auto paths = buildExpandedCva6SlangArgsForSecTests(*context, "cva6");
  auto* top = loadSystemVerilogTopFromPaths(library, "cva6", paths);

  detail::ScopedDnlContextForTest dnlContext(top);
  auto* dnl = dnlContext.dnl();
  ASSERT_NE(dnl, nullptr);
  auto requestedTermID = detail::findTermByDisplayNameForTest(
      dnl,
      "cva6_gen_perf_counter_perf_counters_i.generic_counter_q_mem.WDATA[384]");
  if (!requestedTermID.has_value()) {
    // Keep the probe tied to the generated WDATA port while allowing the local
    // CVA6 config to choose a smaller packed memory width.
    requestedTermID = detail::findFirstTermByDisplayPrefixForTest(
        dnl,
        "cva6_gen_perf_counter_perf_counters_i.generic_counter_q_mem.WDATA[");
  }
  ASSERT_TRUE(requestedTermID.has_value());

  const auto skippedNoDriverReportPath =
      std::filesystem::current_path() / "skipped_no_driver_pos.txt";
  const auto skippedLogicalLoopReportPath =
      std::filesystem::current_path() / "skipped_logical_loop_pos.txt";
  std::filesystem::remove(skippedNoDriverReportPath);
  std::filesystem::remove(skippedLogicalLoopReportPath);
  const bool previousReportSkippedPOs =
      KEPLER_FORMAL::Config::getReportSkippedPOs();
  KEPLER_FORMAL::Config::setReportSkippedPOs(true);
  const auto probe =
      detail::probeRequestedBuilderOutputForTest(dnl, *requestedTermID);
  KEPLER_FORMAL::Config::setReportSkippedPOs(previousReportSkippedPOs);
  std::ostringstream normalizationChain;
  for (size_t index = 0; index < probe.normalizationChain.size(); ++index) {
    if (index != 0) {
      normalizationChain << " -> ";
    }
    normalizationChain << probe.normalizationChain[index];
  }
  std::ostringstream supportTerms;
  for (size_t index = 0; index < probe.rootSupportTerms.size(); ++index) {
    if (index != 0) {
      supportTerms << " | ";
    }
    supportTerms << probe.rootSupportTerms[index];
  }
  std::ostringstream combinationalInputs;
  for (size_t index = 0; index < probe.rootCombinationalInputs.size(); ++index) {
    if (index != 0) {
      combinationalInputs << " | ";
    }
    combinationalInputs << probe.rootCombinationalInputs[index];
  }
  std::ostringstream driverSpine;
  for (size_t index = 0; index < probe.driverSpine.size(); ++index) {
    if (index != 0) {
      driverSpine << " -> ";
    }
    driverSpine << probe.driverSpine[index];
  }
  auto readReportFile = [](const std::filesystem::path& path) {
    std::string contents;
    if (std::ifstream file(path); file.good()) {
      std::ostringstream buffer;
      buffer << file.rdbuf();
      contents = buffer.str();
    }
    return contents;
  };
  const auto skippedNoDriverReport = readReportFile(skippedNoDriverReportPath);
  const auto skippedLogicalLoopReport =
      readReportFile(skippedLogicalLoopReportPath);

  ASSERT_TRUE(probe.normalizedRoot.has_value())
      << "normalization chain: " << normalizationChain.str();
  EXPECT_TRUE(probe.hasBuiltExpr)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " did not yield a valid clause-builder expression; support terms: "
      << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str()
      << "; skipped no-driver report: " << skippedNoDriverReport
      << "; skipped logical-loop report: " << skippedLogicalLoopReport;
  EXPECT_FALSE(probe.hasSkip)
      << "normalized root " << probe.normalizedRootName
      << " model=" << probe.normalizedRootModelName
      << " was skipped: " << probe.skipDetail
      << "; support terms: " << supportTerms.str()
      << "; combinational inputs: " << combinationalInputs.str()
      << "; driver spine: " << driverSpine.str()
      << "; skipped no-driver report: " << skippedNoDriverReport
      << "; skipped logical-loop report: " << skippedLogicalLoopReport;
}

TEST_F(
    SequentialEquivalenceStrategyTests,
    SequentialDesignModelExtractModelsRealCva6PerfCountersTargetConfigMemoryWithoutBoundaryFallback) {
  const auto context = resolveCva6SourceContextForSecTests();
  if (!context.has_value()) {
    GTEST_SKIP()
        << "CVA6 source tree is not available for the real-source SEC memory regression test";
  }

  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadRealCva6PerfCountersTargetConfigTopForSecTests(
      library,
      *context,
      "real_cva6_perf_counters_module_with_target_config_sec_memory_supported");

  const auto extracted = SequentialDesignModel::extract(top);
  auto hasBoundaryRoleName = [&](const std::vector<SignalKey>& keys,
                                 const std::string& needle) {
    return std::any_of(keys.begin(), keys.end(), [&](const SignalKey& key) {
      const auto it = extracted.displayNameByKey.find(key);
      return it != extracted.displayNameByKey.end() &&
             it->second.find(needle) != std::string::npos;
    });
  };

  // Guard the exact configured CVA6 perf-counter memory path that fails in
  // full SEC runs: the inferred memory should stay inside the sequential model
  // instead of leaking back out as generic boundary terms.
  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_FALSE(
      hasBoundaryRoleName(extracted.internalBoundaryInputKeys, "generic_counter_q_mem"));
  EXPECT_FALSE(
      hasBoundaryRoleName(extracted.internalBoundaryOutputKeys, "generic_counter_q_mem"));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractModelsInferredMemoryWithConstantFalseCommitGuard) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadSystemVerilogTopFromSource(
      library,
      "qd_next_indexed_commit_constant_false_guard_skip_supported",
      R"(module qd_next_indexed_commit_constant_false_guard_skip_supported(
  input  logic       clk_i,
  input  logic [1:0] addr_i,
  input  logic [7:0] data_i,
  output logic [7:0] data_o
);
  logic [7:0] mem_q [0:3];
  logic [7:0] mem_d [0:3];
  logic [7:0] mem_next [0:3];

  always_comb begin
    mem_d = mem_q;
    mem_d[addr_i] = data_i;
  end

  always_comb begin
    for (int i = 0; i < 4; i++) begin
      mem_next[i] = mem_d[i];
      if (i == 0) begin
        mem_next[i] = mem_q[i];
      end
    end
  end

  always_ff @(posedge clk_i) begin
    mem_q <= mem_next;
  end

  assign data_o = mem_q[addr_i];
endmodule
)");

  const auto extracted = SequentialDesignModel::extract(top);

  auto hasBoundaryRoleName = [&](const std::vector<SignalKey>& keys,
                                 const std::string& prefix) {
    return std::any_of(keys.begin(), keys.end(), [&](const SignalKey& key) {
      const auto it = extracted.displayNameByKey.find(key);
      return it != extracted.displayNameByKey.end() &&
             it->second.rfind(prefix, 0) == 0;
    });
  };

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryInputKeys, "mem_q"));
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryOutputKeys, "mem_q"));
  EXPECT_FALSE(extracted.stateBits.empty());
  EXPECT_FALSE(extracted.nextStateExprByStateKey.empty());
  EXPECT_TRUE(extracted.skippedObservedOutputs.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractReportsStructuredMemoryDependencyNameForUndrivenAddressBit) {
  // Regression guard for CVA6 diag runs: structured-memory dependency
  // materialization may skip an internal root, and the skip diagnostic must
  // describe that root without using stale DNL terminals after clause building.
  ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadSystemVerilogTopFromSource(
      library,
      "structured_memory_undriven_address_dependency",
      R"(module structured_memory_undriven_address_dependency(
  input  logic       clk_i,
  input  logic       we_i,
  input  logic       addr_i,
  input  logic [7:0] data_i,
  output logic [7:0] data_o
);
  logic [7:0] mem_q [0:3];
  logic       bad_addr_bit;

  always_ff @(posedge clk_i) begin
    if (we_i) begin
      mem_q[{1'b0, addr_i}] <= data_i;
    end
  end

  assign data_o = mem_q[{bad_addr_bit, addr_i}];
endmodule
)");

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_FALSE(extracted.skippedObservedOutputs.empty());
  const auto hasStructuredMemorySkipDetail = std::any_of(
      extracted.connectivitySkipInfoByKey.begin(),
      extracted.connectivitySkipInfoByKey.end(),
      [](const auto& entry) {
        return entry.second.detail.find("Structured memory dependency") !=
                   std::string::npos &&
               entry.second.detail.find("mem_q_mem.RADDR[1]") !=
                   std::string::npos;
      });
  EXPECT_TRUE(hasStructuredMemorySkipDetail);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSkipsUndrivenMemoryWriteEnablePort) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* memoryModel = createSinglePortMemoryModel(primitives, "MEM_DISABLED_WE");
  auto* top = SNLDesign::create(
      library,
      SNLDesign::Type::Standard,
      NLName("structured_memory_undriven_write_enable_port"));
  auto* topClock =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("clk"));
  auto* topChipEnable =
      SNLScalarTerm::create(top, SNLTerm::Direction::Input, NLName("ce"));
  auto* topAddress =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 1, 0, NLName("addr"));
  auto* topWriteData =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 3, 0, NLName("wdata"));
  auto* topWriteMask =
      SNLBusTerm::create(top, SNLTerm::Direction::Input, 3, 0, NLName("wmask"));
  auto* topOut =
      SNLBusTerm::create(top, SNLTerm::Direction::Output, 3, 0, NLName("out"));

  auto* memory = SNLInstance::create(top, memoryModel, NLName("mem0"));
  auto* clockNet = SNLScalarNet::create(top, NLName("clk_net"));
  auto* chipEnableNet = SNLScalarNet::create(top, NLName("ce_net"));
  auto* undrivenWriteEnableNet = SNLScalarNet::create(top, NLName("we_net"));
  auto* addressNet = SNLBusNet::create(top, 1, 0, NLName("addr_net"));
  auto* writeDataNet = SNLBusNet::create(top, 3, 0, NLName("wdata_net"));
  auto* writeMaskNet = SNLBusNet::create(top, 3, 0, NLName("wmask_net"));
  auto* outNet = SNLBusNet::create(top, 3, 0, NLName("out_net"));

  topClock->setNet(clockNet);
  topChipEnable->setNet(chipEnableNet);
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("CLK")))->setNet(clockNet);
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("CE")))->setNet(chipEnableNet);
  memory->getInstTerm(memoryModel->getScalarTerm(NLName("WE")))
      ->setNet(undrivenWriteEnableNet);

  auto* modelAddress = memoryModel->getBusTerm(NLName("ADDR"));
  auto* modelWriteData = memoryModel->getBusTerm(NLName("WDATA"));
  auto* modelWriteMask = memoryModel->getBusTerm(NLName("WMASK"));
  auto* modelReadData = memoryModel->getBusTerm(NLName("RDATA"));
  for (int bit = 0; bit <= 1; ++bit) {
    topAddress->getBit(bit)->setNet(addressNet->getBit(bit));
    memory->getInstTerm(modelAddress->getBit(bit))->setNet(addressNet->getBit(bit));
  }
  for (int bit = 0; bit <= 3; ++bit) {
    topWriteData->getBit(bit)->setNet(writeDataNet->getBit(bit));
    topWriteMask->getBit(bit)->setNet(writeMaskNet->getBit(bit));
    topOut->getBit(bit)->setNet(outNet->getBit(bit));
    memory->getInstTerm(modelWriteData->getBit(bit))->setNet(writeDataNet->getBit(bit));
    memory->getInstTerm(modelWriteMask->getBit(bit))->setNet(writeMaskNet->getBit(bit));
    memory->getInstTerm(modelReadData->getBit(bit))->setNet(outNet->getBit(bit));
  }

  const auto extracted = SequentialDesignModel::extract(top);

  auto findInputVarContaining = [&](std::string_view needle) {
    std::optional<size_t> varID;
    for (const auto& key : extracted.environmentInputs) {
      const auto nameIt = extracted.displayNameByKey.find(key);
      const auto varIt = extracted.inputVarByKey.find(key);
      if (nameIt == extracted.displayNameByKey.end() ||
          varIt == extracted.inputVarByKey.end()) {
        continue;
      }
      if (nameIt->second.find(needle) != std::string::npos) {
        varID = varIt->second;
        break;
      }
    }
    return varID;
  };
  const auto disabledWriteDataVar = findInputVarContaining("wdata[0]");

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  if (disabledWriteDataVar.has_value()) {
    for (const auto& [_, nextStateExpr] : extracted.nextStateExprByStateKey) {
      EXPECT_EQ(nextStateExpr->getSupportVars().count(*disabledWriteDataVar), 0u)
          << "An undriven write-enable port must not pull its write-data cone "
             "into the modeled memory transition";
    }
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSkipsWholeMemoryForUndrivenWriteData) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* memoryModel = createSinglePortMemoryModel(primitives, "MEM_FLOAT_WDATA");
  auto* top = createSinglePortMemoryTop(
      library, "structured_memory_undriven_write_data", memoryModel, 0);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_FALSE(extracted.connectivitySkipInfoByKey.empty());
  EXPECT_FALSE(extracted.skippedObservedOutputs.empty());
  EXPECT_TRUE(std::any_of(
      extracted.connectivitySkipInfoByKey.begin(),
      extracted.connectivitySkipInfoByKey.end(),
      [](const auto& entry) {
        return entry.second.detail.find("Structured memory dependency") !=
               std::string::npos;
      }));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSkipsWholeMemoryForUndrivenReset) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* memoryModel =
      createSinglePortMemoryModel(primitives, "MEM_FLOAT_RST", true);
  auto* top = createSinglePortMemoryTop(
      library,
      "structured_memory_undriven_reset",
      memoryModel,
      std::nullopt,
      true);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_FALSE(extracted.connectivitySkipInfoByKey.empty());
  EXPECT_FALSE(extracted.skippedObservedOutputs.empty());
  EXPECT_TRUE(std::any_of(
      extracted.connectivitySkipInfoByKey.begin(),
      extracted.connectivitySkipInfoByKey.end(),
      [](const auto& entry) {
        return entry.second.detail.find("Structured memory dependency") !=
               std::string::npos;
      }));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractModelsInferredStructMemoryWithLogicalOrCommitGuard) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = loadSystemVerilogTopFromSource(
      library,
      "qd_next_local_commit_logical_or_supported",
      R"(module qd_next_local_commit_logical_or_supported(
  input  logic       clk_i,
  input  logic [1:0] addr_i,
  input  logic [1:0] mode_i,
  input  logic [1:0] access_i,
  input  logic [3:0] payload_i,
  output logic [7:0] data_o
);
  typedef struct packed {
    logic [1:0] mode;
    logic [1:0] access;
    logic [3:0] payload;
  } entry_t;

  entry_t mem_q [0:3];
  entry_t mem_d [0:3];
  entry_t mem_next [0:3];

  always_comb begin
    mem_d = mem_q;
    mem_d[addr_i].mode = mode_i;
    mem_d[addr_i].access = access_i;
    mem_d[addr_i].payload = payload_i;
  end

  always_comb begin
    for (int i = 0; i < 4; i++) begin
      mem_next[i] = mem_d[i];
      if ((mem_d[i].mode == 2'b11) || (mem_d[i].access == 2'b01)) begin
        mem_next[i] = mem_q[i];
      end
    end
  end

  always_ff @(posedge clk_i) begin
    mem_q <= mem_next;
  end

  assign data_o = {mem_q[addr_i].mode, mem_q[addr_i].access, mem_q[addr_i].payload};
endmodule
)");

  const auto extracted = SequentialDesignModel::extract(top);

  auto hasBoundaryRoleName = [&](const std::vector<SignalKey>& keys,
                                 const std::string& prefix) {
    return std::any_of(keys.begin(), keys.end(), [&](const SignalKey& key) {
      const auto it = extracted.displayNameByKey.find(key);
      return it != extracted.displayNameByKey.end() &&
             it->second.rfind(prefix, 0) == 0;
    });
  };

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryInputKeys, "mem_q"));
  EXPECT_FALSE(hasBoundaryRoleName(extracted.internalBoundaryOutputKeys, "mem_q"));
  EXPECT_FALSE(extracted.stateBits.empty());
  EXPECT_FALSE(extracted.nextStateExprByStateKey.empty());
  EXPECT_TRUE(extracted.skippedObservedOutputs.empty());
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
       SequentialDesignModelExtractPropagatesNestedNoDriverDataConeSkips) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* andModel = createAnd2Model(primitives);
  auto* top =
      createPartialCoverageNoDriverDataConeTop(library, "top", andModel);

  const auto extracted = SequentialDesignModel::extract(top);

  // A no-driver buried below a sequential data cone makes the state update
  // unmodelable. Outputs that depend on that state must be skipped instead of
  // comparing arbitrary private symbols between the two SEC sides.
  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.totalObservedOutputCount(), 2u);
  EXPECT_EQ(extracted.coveredObservedOutputCount(), 1u);
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
       SequentialDesignModelExtractPropagatesLogicalLoopSkipsToStateAndOutputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top = createPartialCoverageLogicalLoopTop(library, "top");

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
      ConnectivitySkipOrigin::LogicalLoop);
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
       SequentialDesignModelExtractFoldsOpaqueClockGateLatchEnable) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* andModel = createAnd2Model(primitives);
  auto* latchModel = createOpaqueClockGateLatchModel(primitives);
  auto* top = createClockGateLatchDffTop(
      library, "top", andModel, latchModel);

  const auto extracted = SequentialDesignModel::extract(top);
  expectAllExpressionSupportIsPublished(extracted);
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto enableKey = findKeyByDisplayName(extracted, "en[0]");
  const auto latchKey =
      findKeyByDisplayName(extracted, "clock_gate_i.en_latch.Q[0]");
  const size_t stateVar = extracted.inputVarByKey.at(stateKey);
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.inputVarByKey.find(latchKey), extracted.inputVarByKey.end());
  EXPECT_EQ(
      std::find(
          extracted.environmentInputs.begin(),
          extracted.environmentInputs.end(),
          latchKey),
      extracted.environmentInputs.end());
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(enableKey), false},
       {stateVar, true}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(enableKey), true},
       {stateVar, false}}));
  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(enableKey), true},
       {stateVar, true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractTreatsNamedClockTreeBufferAsCarrier) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* bufferModel = createBufModel(primitives);
  auto* top = createClockTreeBufferedDffTop(library, "top", bufferModel);

  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto extracted = SequentialDesignModel::extract(top);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();
  expectAllExpressionSupportIsPublished(extracted);
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto clockKey = findKeyByDisplayName(extracted, "clk[0]");
  const size_t stateVar = extracted.inputVarByKey.at(stateKey);
  const size_t clockVar = extracted.inputVarByKey.at(clockKey);
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(expr->getSupportVars().count(clockVar), 0u);
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true}, {stateVar, false}}));
  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false}, {stateVar, true}}));
  const std::string structureIndexBuilt =
      "immutable clock structure indexed terms=";
  const size_t firstStructureIndex = stderrOutput.find(structureIndexBuilt);
  ASSERT_NE(firstStructureIndex, std::string::npos) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          structureIndexBuilt,
          firstStructureIndex + structureIndexBuilt.size()),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractDoesNotTreatDataClkbufCellAsClockCarrier) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* dataBufferModel =
      createBufModelWithName(primitives, "sky130_fd_sc_hs__clkbuf_1");
  auto* top = createDataBufferedDffTop(library, "top", dataBufferModel);

  const auto extracted = SequentialDesignModel::extract(top);
  expectAllExpressionSupportIsPublished(extracted);
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto clockKey = findKeyByDisplayName(extracted, "clk[0]");
  const size_t inVar = extracted.inputVarByKey.at(inKey);
  const size_t clockVar = extracted.inputVarByKey.at(clockKey);
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  // sky130hs uses clkbuf cells for ordinary routed data too.  Only a branch
  // that structurally traces to the top clock may become a SEC clock carrier.
  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(expr->getSupportVars().count(inVar), 1u);
  EXPECT_EQ(expr->getSupportVars().count(clockVar), 0u);
  EXPECT_NE(
      std::find(
          extracted.clockCarrierVarIDs.begin(),
          extracted.clockCarrierVarIDs.end(),
          clockVar),
      extracted.clockCarrierVarIDs.end());
  EXPECT_EQ(
      std::find(
          extracted.clockCarrierVarIDs.begin(),
          extracted.clockCarrierVarIDs.end(),
          inVar),
      extracted.clockCarrierVarIDs.end());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractSubstitutesFoldedClockGateLatchDataUses) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* andModel = createAnd2Model(primitives);
  auto* latchModel = createOpaqueClockGateLatchModel(primitives);
  auto* top = createClockGateLatchDataDffTop(
      library, "top", andModel, latchModel);

  const auto extracted = SequentialDesignModel::extract(top);
  expectAllExpressionSupportIsPublished(extracted);
  const auto stateKey = findKeyByDisplayName(extracted, "ff0.Q[0]");
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto enableKey = findKeyByDisplayName(extracted, "en[0]");
  const auto latchKey =
      findKeyByDisplayName(extracted, "clock_gate_i.en_latch.Q[0]");
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.inputVarByKey.find(latchKey), extracted.inputVarByKey.end());
  EXPECT_EQ(
      std::find(
          extracted.environmentInputs.begin(),
          extracted.environmentInputs.end(),
          latchKey),
      extracted.environmentInputs.end());
  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(enableKey), false}}));
  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(enableKey), true}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(enableKey), true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractPreservesUnrelatedClockGateLatchCones) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* andModel = createAnd2Model(primitives);
  auto* latchModel = createOpaqueClockGateLatchModel(primitives);
  auto* top = createClockGateLatchDataDffTop(
      library, "top", andModel, latchModel, true);

  const auto extracted = SequentialDesignModel::extract(top);
  expectAllExpressionSupportIsPublished(extracted);
  const auto independentStateKey =
      findKeyByDisplayName(extracted, "independent_ff.Q[0]");
  const auto independentInKey =
      findKeyByDisplayName(extracted, "independent_in[0]");
  const auto enableKey = findKeyByDisplayName(extracted, "en[0]");
  const auto latchKey =
      findKeyByDisplayName(extracted, "clock_gate_i.en_latch.Q[0]");
  auto* independentExpr =
      extracted.nextStateExprByStateKey.at(independentStateKey);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_EQ(extracted.inputVarByKey.find(latchKey), extracted.inputVarByKey.end());
  EXPECT_EQ(independentExpr->getSupportVars().count(
                extracted.inputVarByKey.at(enableKey)),
            0u);
  EXPECT_TRUE(independentExpr->evaluate(
      {{extracted.inputVarByKey.at(independentInKey), true}}));
  EXPECT_FALSE(independentExpr->evaluate(
      {{extracted.inputVarByKey.at(independentInKey), false}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelExtractDoesNotPublishConstantInternalBoundaryInputs) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* constantModel = createConstantLowModel(primitives);
  auto* top = createConstantDrivenDffTop(library, "top", constantModel);

  const auto extracted = SequentialDesignModel::extract(top);
  expectAllExpressionSupportIsPublished(extracted);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  for (const auto& key : extracted.environmentInputs) {
    const auto nameIt = extracted.displayNameByKey.find(key);
    ASSERT_NE(nameIt, extracted.displayNameByKey.end());
    EXPECT_NE(nameIt->second, "tie0.LO[0]");
  }
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
       SequentialDesignModelExtractSupportsResetSetControlStyles) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* model = createResetSetSequentialModel(primitives, "SEQ_RST_SET");
  auto* top = createResetSetSequentialTop(library, "top", model);

  const auto extracted = SequentialDesignModel::extract(top);

  EXPECT_FALSE(extracted.hasUnsupportedFeatures());
  EXPECT_TRUE(extracted.abstractedSequentialBoundaries.empty());
  EXPECT_TRUE(extracted.internalBoundaryOutputKeys.empty());
  ASSERT_EQ(extracted.stateBits.size(), 1u);
  ASSERT_EQ(extracted.topOutputKeys.size(), 1u);
  EXPECT_EQ(extracted.observedOutputs.size(), extracted.topOutputKeys.size());
  EXPECT_TRUE(extracted.initialStateValueByKey.empty());

  const auto stateKey = extracted.stateBits.front();
  const auto inKey = findKeyByDisplayName(extracted, "in[0]");
  const auto rstKey = findKeyByDisplayName(extracted, "rst[0]");
  const auto setKey = findKeyByDisplayName(extracted, "set[0]");
  auto* expr = extracted.nextStateExprByStateKey.at(stateKey);

  EXPECT_FALSE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(rstKey), true},
       {extracted.inputVarByKey.at(setKey), false}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(rstKey), false},
       {extracted.inputVarByKey.at(setKey), true}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), true},
       {extracted.inputVarByKey.at(rstKey), false},
       {extracted.inputVarByKey.at(setKey), false}}));
  EXPECT_TRUE(expr->evaluate(
      {{extracted.inputVarByKey.at(inKey), false},
       {extracted.inputVarByKey.at(rstKey), true},
       {extracted.inputVarByKey.at(setKey), true}}));
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
       MismatchedObservedOutputNamesThrow) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedOutputDffTop(library, "top0", invModel, "out0");
  auto* top1 = createNamedOutputDffTop(library, "top1", invModel, "out1");

  auto strategy = makeBinarySecStrategy(top0, top1);
  EXPECT_THROW(static_cast<void>(strategy.run(1)), std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedObservedOutputCountsThrow) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedOutputDffTop(library, "top0", invModel, "out");
  auto* top1 = createExtraOutputDffTop(library, "top1", invModel);

  auto strategy = makeBinarySecStrategy(top0, top1);
  EXPECT_THROW(static_cast<void>(strategy.run(1)), std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedInputNamesThrow) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedInputDffTop(library, "top0", invModel, "in0");
  auto* top1 = createNamedInputDffTop(library, "top1", invModel, "in1");

  auto strategy = makeBinarySecStrategy(top0, top1);
  EXPECT_THROW(static_cast<void>(strategy.run(1)), std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       MismatchedInputCountsThrow) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createNamedInputDffTop(library, "top0", invModel, "in");
  auto* top1 = createExtraInputDffTop(library, "top1");

  auto strategy = makeBinarySecStrategy(top0, top1);
  EXPECT_THROW(static_cast<void>(strategy.run(1)), std::runtime_error);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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
  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(result.skippedObservedOutputs.front().find("no-driver"), std::string::npos);
  EXPECT_TRUE(result.resetUnanchoredSkippedOutputs.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentDesignsCanProveSecWhenOutputSkipOnlyAppearsOnOneSide) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createPartialCoverageNoDriverTop(library, "top0");
  auto* top1 = createPartialCoverageDrivenTop(library, "top1");

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(result.skippedObservedOutputs.front().find("design0"), std::string::npos);
  EXPECT_EQ(result.skippedObservedOutputs.front().find("design1"), std::string::npos);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(result.skippedObservedOutputs.front().find("multi-driver"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentDesignsCanProveSecOnCoveredOutputsOnlyAfterLogicalLoopSkipping) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createPartialCoverageLogicalLoopTop(library, "top0");
  auto* top1 = createPartialCoverageLogicalLoopTop(library, "top1");

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::PartiallyProved);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(
      result.skippedObservedOutputs.front().find("logical-loop"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       UninitializedDffDesignsReportTopBoundarySurface) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::KInduction);
  const auto result = strategy.run(2);

  auto hasRole = [&](const char* design, const char* signal, const char* role) {
    return std::any_of(
        result.extractedBoundaryReports.begin(),
        result.extractedBoundaryReports.end(),
        [&](const ExtractedBoundaryReportEntry& entry) {
          return entry.design == design && entry.signal == signal &&
                 std::find(entry.roles.begin(), entry.roles.end(), role) !=
                     entry.roles.end();
        });
  };

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_TRUE(hasRole("design0", "clk[0]", "top_input"));
  EXPECT_TRUE(hasRole("design0", "in[0]", "top_input"));
  EXPECT_TRUE(hasRole("design0", "out[0]", "top_output"));
  EXPECT_TRUE(hasRole("design1", "clk[0]", "top_input"));
  EXPECT_TRUE(hasRole("design1", "in[0]", "top_input"));
  EXPECT_TRUE(hasRole("design1", "out[0]", "top_output"));
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentOpaqueLeafDesignsReportOpaqueInternalBoundaryTerms) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* opaqueModel = createOpaqueLeafModel(primitives);
  auto* top0 = createOpaqueBoundaryTop(library, "top0", opaqueModel);
  auto* top1 = createOpaqueBoundaryTop(library, "top1", opaqueModel);

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  auto hasRole = [&](const char* design, const char* signal, const char* role) {
    return std::any_of(
        result.extractedBoundaryReports.begin(),
        result.extractedBoundaryReports.end(),
        [&](const ExtractedBoundaryReportEntry& entry) {
          return entry.design == design && entry.signal == signal &&
                 std::find(entry.roles.begin(), entry.roles.end(), role) !=
                     entry.roles.end();
        });
  };

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_TRUE(hasRole("design0", "in[0]", "top_input"));
  EXPECT_TRUE(hasRole("design0", "out[0]", "top_output"));
  EXPECT_TRUE(hasRole("design0", "opaque0.Y[0]", "opaque_internal_input"));
  EXPECT_TRUE(hasRole("design0", "opaque0.A[0]", "opaque_internal_output"));
  EXPECT_TRUE(hasRole("design1", "opaque0.Y[0]", "opaque_internal_input"));
  EXPECT_TRUE(hasRole("design1", "opaque0.A[0]", "opaque_internal_output"));
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  auto hasRole = [&](const char* design, const char* signal, const char* role) {
    return std::any_of(
        result.extractedBoundaryReports.begin(),
        result.extractedBoundaryReports.end(),
        [&](const ExtractedBoundaryReportEntry& entry) {
          return entry.design == design && entry.signal == signal &&
                 std::find(entry.roles.begin(), entry.roles.end(), role) !=
                     entry.roles.end();
        });
  };

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_FALSE(result.abstractedSequentialBoundaries.empty());
  EXPECT_TRUE(hasRole("design0", "clk[0]", "top_input"));
  EXPECT_TRUE(hasRole("design0", "in[0]", "top_input"));
  EXPECT_TRUE(hasRole("design0", "good[0]", "top_output"));
  EXPECT_TRUE(hasRole("design0", "bad[0]", "top_output"));
  EXPECT_FALSE(hasRole("design0", "ff0.Q[0]", "opaque_internal_input"));
  EXPECT_FALSE(hasRole("design0", "ff0.CK[0]", "opaque_internal_output"));
  EXPECT_TRUE(hasRole("design0", "ff0.Q[0]", "abstracted_sequential_state"));
  EXPECT_TRUE(hasRole("design0", "ff0.CK[0]", "abstracted_sequential_observed"));
  EXPECT_FALSE(hasRole("design1", "ff0.Q[0]", "opaque_internal_input"));
  EXPECT_FALSE(hasRole("design1", "ff0.CK[0]", "opaque_internal_output"));
  EXPECT_TRUE(hasRole("design1", "ff0.Q[0]", "abstracted_sequential_state"));
  EXPECT_TRUE(hasRole("design1", "ff0.CK[0]", "abstracted_sequential_observed"));
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(2);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Unsupported);
  EXPECT_TRUE(result.skippedObservedOutputs.empty());
  EXPECT_FALSE(result.reason.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       DiagnosticModePrintsStrategyAndExtractionProgress) {
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
  ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      top0,
      top1,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.run(3);
  const std::string stdoutOutput = testing::internal::GetCapturedStdout();
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_NE(stderrOutput.find("SEC diag: start run"), std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: extract(top0) collect begin"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: entering pdr engine"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("singleton CaDiCaL cumulative budget"),
      std::string::npos);
  EXPECT_NE(stdoutOutput.find("SEC diag: aligned_inputs="), std::string::npos);
  EXPECT_NE(
      stdoutOutput.find("SEC summary: encoding=dual_rail_steady"),
      std::string::npos);
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
       KInductionBatchedOutputsCombineInconclusiveSlices) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(2));
  problem.initialStateAssignments = {{2, false}};
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  problem.transitions0 = {{2, BoolExpr::Var(2)}};
  for (size_t i = 0; i < 33; ++i) {
    problem.observedOutputs.push_back(makeSignalKey("batched_inconclusive_" + std::to_string(i)));
    problem.observedOutputNames.push_back("out_" + std::to_string(i));
    problem.observedOutputExprs0.push_back(BoolExpr::Var(2));
    problem.observedOutputExprs1.push_back(BoolExpr::Var(2));
  }
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);

  EXPECT_EQ(result.status, KInductionStatus::Inconclusive);
  EXPECT_EQ(result.bound, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       KInductionBatchedOutputsReturnDifferentSlice) {
  KInductionProblem problem;
  for (size_t i = 0; i < 33; ++i) {
    problem.observedOutputs.push_back(makeSignalKey("batched_different_" + std::to_string(i)));
    problem.observedOutputNames.push_back("out_" + std::to_string(i));
    problem.observedOutputExprs0.push_back(i == 0 ? BoolExpr::createFalse()
                                                  : BoolExpr::createTrue());
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
  }
  problem.property = BoolExpr::createFalse();
  problem.bad = BoolExpr::createTrue();

  KInductionEngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);

  EXPECT_EQ(result.status, KInductionStatus::Different);
  ASSERT_TRUE(result.witness.has_value());
  EXPECT_EQ(result.witness->badFrame, 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       UnpublishedInternalSupportIsDesignPrivateDuringSecRemap) {
  SequentialDesignModel model0;
  SequentialDesignModel model1;
  const SignalKey out = makeSignalKey("privateSupportOut");
  const SignalKey out0 = out;
  const SignalKey out1 = out;
  model0.topOutputKeys = {out0};
  model0.allObservedOutputs = {out0};
  model0.observedOutputs = {out0};
  model0.displayNameByKey.emplace(out0, "out[0]");
  model0.observedOutputExprByKey.emplace(out0, BoolExpr::Var(42));

  model1.topOutputKeys = {out1};
  model1.allObservedOutputs = {out1};
  model1.observedOutputs = {out1};
  model1.displayNameByKey.emplace(out1, "out[0]");
  model1.observedOutputExprByKey.emplace(out1, BoolExpr::Var(42));

  auto strategy = makeBinaryExtractedSecStrategy(SecEngine::Imc);
  const auto result = strategy.runExtractedModels(model0, model1, 1);

  // Same local variable ID, different extracted designs: this support is not a
  // top terminal, so SEC must not equate it by name or by local BoolExpr ID.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelDetailHelpersCoverNextStateErrors) {
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
  auto* resetSetExpr = detail::buildNextStateExprForTest(
      0, {{"D", 11}, {"R", 12}, {"S", 13}}, {2}, outputExprByTerm);
  EXPECT_FALSE(
      resetSetExpr->evaluate({{2, true}, {7, true}, {8, true}, {9, false}}));
  EXPECT_TRUE(
      resetSetExpr->evaluate({{2, true}, {7, false}, {8, false}, {9, true}}));
  EXPECT_TRUE(
      resetSetExpr->evaluate({{2, false}, {7, false}, {8, true}, {9, true}}));
  auto* holdExpr = detail::buildNextStateExprForTest(
      0, {{"D", 11}, {"E", 12}, {"RN", 13}}, {2}, outputExprByTerm);
  EXPECT_FALSE(holdExpr->evaluate({{2, true}, {7, true}, {8, true}, {9, false}}));
  EXPECT_TRUE(holdExpr->evaluate({{2, false}, {7, true}, {8, true}, {9, true}}));
  EXPECT_TRUE(holdExpr->evaluate({{2, true}, {7, false}, {8, false}, {9, true}}));
  auto* setExpr = detail::buildNextStateExprForTest(
      0, {{"D", 11}, {"S", 12}}, {2}, outputExprByTerm);
  EXPECT_TRUE(setExpr->evaluate({{2, false}, {7, false}, {8, true}}));
  EXPECT_FALSE(setExpr->evaluate({{2, false}, {7, false}, {8, false}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprVariableRemapRejectsUnsupportedOperator) {
  BoolExpr invalidExpr;

  EXPECT_THROW(
      static_cast<void>(remapBoolExprVariables(&invalidExpr, {})),
      std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprVariableRemapPreservesIdentitySubtrees) {
  BoolExpr* left = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  BoolExpr* right = BoolExpr::Or(BoolExpr::Var(4), BoolExpr::Var(5));
  BoolExpr* root = BoolExpr::Xor(left, right);

  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  EXPECT_EQ(
      remapBoolExprVariables(
          root, {{2, 2}, {3, 3}, {4, 4}, {5, 5}}, memo),
      root);

  memo.clear();
  BoolExpr* remapped =
      remapBoolExprVariables(root, {{2, 2}, {3, 6}, {4, 4}, {5, 5}}, memo);
  EXPECT_NE(remapped, root);
  EXPECT_NE(remapped->getLeft(), left);
  // The rebuilt parent may be canonicalized by the global BoolExpr cache when
  // this test runs after other large SEC tests. The remapper contract is that
  // unchanged sub-DAGs are preserved in the memo before parent construction.
  ASSERT_NE(memo.find(right), memo.end());
  EXPECT_EQ(memo.at(right), right);
}

TEST_F(SequentialEquivalenceStrategyTests,
       BoolExprVariableRemapCanReuseBatchMemoForSharedCones) {
  BoolExpr* shared = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  BoolExpr* firstRoot = BoolExpr::Or(shared, BoolExpr::Var(4));
  BoolExpr* secondRoot = BoolExpr::Xor(shared, BoolExpr::Var(5));

  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  static_cast<void>(
      remapBoolExprVariables(firstRoot, {{2, 7}, {3, 8}, {4, 4}}, memo));
  ASSERT_NE(memo.find(shared), memo.end());
  BoolExpr* remappedShared = memo.at(shared);

  // Dependency materialization remaps many output formulas from one builder
  // batch. Reusing the memo across outputs is safe because the variable map is
  // fixed, and it keeps large shared cones from being rebuilt repeatedly.
  static_cast<void>(
      remapBoolExprVariables(secondRoot, {{2, 7}, {3, 8}, {5, 5}}, memo));
  ASSERT_NE(memo.find(shared), memo.end());
  EXPECT_EQ(memo.at(shared), remappedShared);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SequentialDesignModelDetailSelectsRequiredBuilderOutputs) {
  const auto requiredOutputs = detail::selectRequiredBuilderOutputsForTest(
      {10, 11, 12, 13, 14},
      {10, 14},
      {12, 13, 13},
      {14});
  EXPECT_EQ(
      requiredOutputs,
      (std::vector<naja::DNL::DNLID>{10, 12, 13}));

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

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelHelpersCoverPhaseAndUngatedEventPaths) {
  ClockEvent defaultEvent;
  EXPECT_EQ(defaultEvent.phase, ClockPhase::Pos);
  EXPECT_EQ(defaultEvent.enable, nullptr);

  const auto domain = makeSignalKey("clk");
  ClockEvent ungated{domain, ClockPhase::Pos, nullptr};
  EXPECT_EQ(invertClockPhase(ClockPhase::Pos), ClockPhase::Neg);
  EXPECT_EQ(invertClockPhase(ClockPhase::Neg), ClockPhase::Pos);
  EXPECT_STREQ(clockPhaseName(ClockPhase::Pos), "posedge");
  EXPECT_STREQ(clockPhaseName(ClockPhase::Neg), "negedge");
  EXPECT_TRUE(clockEventIsUngated(ungated));
  EXPECT_EQ(clockEventEnableOrTrue(ungated), BoolExpr::createTrue());

  ClockEvent explicitTrue{domain, ClockPhase::Pos, BoolExpr::createTrue()};
  EXPECT_TRUE(clockEventIsUngated(explicitTrue));
  EXPECT_EQ(clockEventEnableOrTrue(explicitTrue), BoolExpr::createTrue());

  ClockEvent gated{domain, ClockPhase::Neg, BoolExpr::Var(42)};
  EXPECT_FALSE(clockEventIsUngated(gated));
  EXPECT_EQ(clockEventEnableOrTrue(gated), BoolExpr::Var(42));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelClassifiesCarrierEdgesAndGates) {
  const auto domain = makeSignalKey("clk");
  const std::unordered_map<size_t, ClockEvent> carrierEvents = {
      {10, ClockEvent{domain, ClockPhase::Pos, nullptr}},
  };

  const auto posEdge = classifyClockEventExpression(BoolExpr::Var(10), carrierEvents);
  ASSERT_TRUE(posEdge.has_value());
  EXPECT_EQ(posEdge->domain, domain);
  EXPECT_EQ(posEdge->phase, ClockPhase::Pos);
  EXPECT_EQ(posEdge->enable, nullptr);

  const auto negEdge =
      classifyClockEventExpression(BoolExpr::Not(BoolExpr::Var(10)), carrierEvents);
  ASSERT_TRUE(negEdge.has_value());
  EXPECT_EQ(negEdge->domain, domain);
  EXPECT_EQ(negEdge->phase, ClockPhase::Neg);
  EXPECT_EQ(negEdge->enable, nullptr);

  const auto andGated = classifyClockEventExpression(
      BoolExpr::And(BoolExpr::Var(10), BoolExpr::Var(20)), carrierEvents);
  ASSERT_TRUE(andGated.has_value());
  EXPECT_EQ(andGated->phase, ClockPhase::Pos);
  ASSERT_NE(andGated->enable, nullptr);
  EXPECT_FALSE(andGated->enable->evaluate({{20, false}}));
  EXPECT_TRUE(andGated->enable->evaluate({{20, true}}));

  const auto orGated = classifyClockEventExpression(
      BoolExpr::Or(BoolExpr::Var(10), BoolExpr::Var(21)), carrierEvents);
  ASSERT_TRUE(orGated.has_value());
  EXPECT_EQ(orGated->phase, ClockPhase::Pos);
  ASSERT_NE(orGated->enable, nullptr);
  EXPECT_TRUE(orGated->enable->evaluate({{21, false}}));
  EXPECT_FALSE(orGated->enable->evaluate({{21, true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelRejectsAmbiguousAndMissingCarrierExpressions) {
  const auto domain = makeSignalKey("clk");
  const std::unordered_map<size_t, ClockEvent> carrierEvents = {
      {10, ClockEvent{domain, ClockPhase::Pos, nullptr}},
      {11, ClockEvent{makeSignalKey("clk2"), ClockPhase::Pos, nullptr}},
  };

  EXPECT_FALSE(classifyClockEventExpression(nullptr, carrierEvents).has_value());
  EXPECT_FALSE(classifyClockEventExpression(BoolExpr::Var(10), {}).has_value());
  EXPECT_FALSE(
      classifyClockEventExpression(
          BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::Var(11)),
          carrierEvents)
          .has_value());
  EXPECT_FALSE(
      classifyClockEventExpression(
          BoolExpr::Xor(BoolExpr::Var(10), BoolExpr::Var(20)),
          carrierEvents)
          .has_value());
  EXPECT_FALSE(
      classifyClockEventExpression(BoolExpr::Var(99), carrierEvents).has_value());
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelCoversRawDagCarrierSpecializationEdges) {
  const auto domain = makeSignalKey("clk");
  const std::unordered_map<size_t, ClockEvent> carrierEvents = {
      {10, ClockEvent{domain, ClockPhase::Pos, nullptr}},
  };

  const auto makeRaw = [](KEPLER_FORMAL::Op op,
                          BoolExpr* left,
                          BoolExpr* right) {
    return KEPLER_FORMAL::BoolExprCache::getExpression({op, 0, left, right});
  };
  auto* clock = BoolExpr::Var(10);
  auto* gateA = BoolExpr::Var(20);
  auto* gateB = BoolExpr::Var(21);
  auto* gatedAnd = makeRaw(KEPLER_FORMAL::Op::AND, gateA, gateB);
  auto* gatedOr = makeRaw(KEPLER_FORMAL::Op::OR, gateA, gateB);
  auto* gatedXor = makeRaw(KEPLER_FORMAL::Op::XOR, gateA, gateB);

  for (auto* gateExpr : {gatedAnd, gatedOr, gatedXor}) {
    const auto edge =
        classifyClockEventExpression(BoolExpr::And(clock, gateExpr), carrierEvents);
    ASSERT_TRUE(edge.has_value());
    EXPECT_EQ(edge->phase, ClockPhase::Pos);
    EXPECT_EQ(edge->enable, gateExpr);
  }

  auto* rawClockAndTrue =
      makeRaw(KEPLER_FORMAL::Op::AND, clock, BoolExpr::createTrue());
  auto* shared = makeRaw(KEPLER_FORMAL::Op::AND, rawClockAndTrue, gateA);
  auto* repeated = makeRaw(KEPLER_FORMAL::Op::OR, shared, shared);
  const auto repeatedEdge =
      classifyClockEventExpression(repeated, carrierEvents);
  ASSERT_TRUE(repeatedEdge.has_value());
  EXPECT_EQ(repeatedEdge->phase, ClockPhase::Pos);
  ASSERT_NE(repeatedEdge->enable, nullptr);
  EXPECT_FALSE(repeatedEdge->enable->evaluate({{20, false}}));
  EXPECT_TRUE(repeatedEdge->enable->evaluate({{20, true}}));

  auto* invalid =
      makeRaw(KEPLER_FORMAL::Op::NONE, BoolExpr::Var(10), nullptr);
  EXPECT_THROW(classifyClockEventExpression(invalid, carrierEvents),
               std::runtime_error);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelSubstitutesVariableExpressionsAcrossSharedDag) {
  BoolExpr* shared = BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3));
  BoolExpr* root = BoolExpr::Xor(
      BoolExpr::Or(BoolExpr::Not(shared), BoolExpr::Var(4)),
      BoolExpr::And(shared, BoolExpr::Var(5)));

  const auto unchanged = substituteBoolExprVariableExpressions(root, {});
  EXPECT_EQ(unchanged, root);
  EXPECT_EQ(substituteBoolExprVariableExpressions(nullptr, {{2, BoolExpr::Var(9)}}),
            nullptr);

  BoolExpr* replacement = BoolExpr::Or(BoolExpr::Var(8), BoolExpr::Var(9));
  BoolExpr* substituted = substituteBoolExprVariableExpressions(
      root, {{2, replacement}, {5, BoolExpr::createTrue()}});
  ASSERT_NE(substituted, nullptr);
  EXPECT_NE(substituted, root);
  EXPECT_EQ(
      substituted->evaluate({{3, true}, {4, false}, {8, false}, {9, false}}),
      root->evaluate({{2, false}, {3, true}, {4, false}, {5, true}}));
  EXPECT_EQ(
      substituted->evaluate({{3, true}, {4, true}, {8, true}, {9, false}}),
      root->evaluate({{2, true}, {3, true}, {4, true}, {5, true}}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelKeepsUnchangedCompositeSubstitutionNodes) {
  const std::unordered_map<size_t, BoolExpr*> replacements = {
      {2, BoolExpr::Var(9)},
  };

  auto* andExpr = BoolExpr::And(BoolExpr::Var(20), BoolExpr::Var(21));
  auto* orExpr = BoolExpr::Or(BoolExpr::Var(20), BoolExpr::Var(21));
  auto* xorExpr = BoolExpr::Xor(BoolExpr::Var(20), BoolExpr::Var(21));

  EXPECT_EQ(substituteBoolExprVariableExpressions(andExpr, replacements), andExpr);
  EXPECT_EQ(substituteBoolExprVariableExpressions(orExpr, replacements), orExpr);
  EXPECT_EQ(substituteBoolExprVariableExpressions(xorExpr, replacements), xorExpr);
}

TEST_F(SequentialEquivalenceStrategyTests,
       SecClockModelRejectsInvalidBoolExprOperators) {
  BoolExpr invalid;
  EXPECT_THROW(
      substituteBoolExprVariableExpressions(&invalid, {{2, BoolExpr::Var(9)}}),
      std::runtime_error);
}
