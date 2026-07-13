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
#include "strategy/ReachableStateInvariant.h"
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

std::string normalizePinNameForTest(const std::string& name) {
  std::string normalized = name;
  for (char& ch : normalized) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return normalized;
}

std::string normalizeSignalBaseNameForTest(const std::string& name) {
  std::string base = name;
  const auto bracket = base.find('[');
  if (bracket != std::string::npos) {
    base = base.substr(0, bracket);
  }
  return normalizePinNameForTest(base);
}

bool hasSuffixForTest(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isResetNameTokenForTest(
    const std::string& candidate,
    const std::string& token) {
  return candidate == token || hasSuffixForTest(candidate, "_" + token);
}

bool isActiveLowResetTokenForTest(const std::string& candidate) {
  return candidate == "RESET_N" || candidate == "RESETN" ||
         candidate == "RESET_L" || candidate == "RST_N" ||
         candidate == "RSTN" || candidate == "RST_L";
}

void appendDomainPrefixedActiveLowResetCandidatesForTest(
    std::vector<std::string>& candidates) {
  const size_t originalSize = candidates.size();
  for (size_t index = 0; index < originalSize; ++index) {
    const std::string& candidate = candidates[index];
    if (candidate.size() <= 1) {
      continue;
    }
    const std::string strippedDomain = candidate.substr(1);
    if (isActiveLowResetTokenForTest(strippedDomain)) {
      candidates.push_back(strippedDomain);
    }
  }
}

std::optional<bool> getResetAssertionValueFromDisplayNameForTest(
    const std::string& displayName);

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

std::optional<bool> evaluateConstantUnderAssignmentsImplForTest(
    BoolExpr* expr,
    const std::unordered_map<size_t, bool>& assignments,
    std::unordered_map<BoolExpr*, std::optional<bool>>& memo) {
  if (expr == nullptr) {
    return std::nullopt;
  }
  if (const auto it = memo.find(expr); it != memo.end()) {
    return it->second;
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
      const auto operand = evaluateConstantUnderAssignmentsImplForTest(
          expr->getLeft(), assignments, memo);
      if (operand.has_value()) {
        value = !*operand;
      }
      break;
    }
    case Op::AND: {
      const auto lhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getLeft(), assignments, memo);
      if (lhs.has_value() && !*lhs) {
        value = false;
        break;
      }
      const auto rhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getRight(), assignments, memo);
      if (rhs.has_value() && !*rhs) {
        value = false;
      } else if (lhs.has_value() && rhs.has_value()) {
        value = *lhs && *rhs;
      }
      break;
    }
    case Op::OR: {
      const auto lhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getLeft(), assignments, memo);
      if (lhs.has_value() && *lhs) {
        value = true;
        break;
      }
      const auto rhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getRight(), assignments, memo);
      if (rhs.has_value() && *rhs) {
        value = true;
      } else if (lhs.has_value() && rhs.has_value()) {
        value = *lhs || *rhs;
      }
      break;
    }
    case Op::XOR: {
      const auto lhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getLeft(), assignments, memo);
      const auto rhs = evaluateConstantUnderAssignmentsImplForTest(
          expr->getRight(), assignments, memo);
      if (lhs.has_value() && rhs.has_value()) {
        value = *lhs != *rhs;
      }
      break;
    }
    case Op::NONE:
    default:
      break;
  }

  memo.emplace(expr, value);
  return value;
}

std::unordered_map<size_t, bool> collectResetAssignmentsForTest(
    const SequentialDesignModel& model) {
  std::unordered_map<size_t, bool> assignments;
  for (const auto& key : model.environmentInputs) {
    const auto displayIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        varIt == model.inputVarByKey.end()) {
      continue;
    }
    const auto assertedValue =
        getResetAssertionValueFromDisplayNameForTest(displayIt->second);
    if (!assertedValue.has_value()) {
      continue;
    }
    assignments.emplace(varIt->second, *assertedValue);
  }
  return assignments;
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

std::optional<bool> detectInitialStateValueForTest(
    const std::unordered_map<std::string, naja::DNL::DNLID>& pinTermIDs) {
  PendingTransitionForTest pending;
  pending.independentStateOutputCount = 1;
  for (const auto& [pinName, termID] : pinTermIDs) {
    pending.pinTermIDs[pinName].push_back({termID, 0});
  }

  const bool hasResetHigh = resolvePendingPinTermIDForTest(pending, "R").has_value();
  const bool hasResetLow = resolvePendingPinTermIDForTest(pending, "RN").has_value();
  const bool hasSetHigh = resolvePendingPinTermIDForTest(pending, "S").has_value();
  const bool hasSetLow = resolvePendingPinTermIDForTest(pending, "SN").has_value();

  const bool hasReset = hasResetHigh || hasResetLow;
  const bool hasSet = hasSetHigh || hasSetLow;
  if (hasReset && !hasSet) {
    return false;
  }
  if (hasSet && !hasReset) {
    return true;
  }
  return std::nullopt;
}

std::optional<bool> evaluateConstantUnderAssignmentsForTest(
    BoolExpr* expr,
    const std::unordered_map<size_t, bool>& assignments) {
  std::unordered_map<BoolExpr*, std::optional<bool>> memo;
  return evaluateConstantUnderAssignmentsImplForTest(expr, assignments, memo);
}

void inferSynthesizedResetInitialStateValuesForTest(SequentialDesignModel& model) {
  const auto resetAssignments = collectResetAssignmentsForTest(model);
  if (resetAssignments.empty()) {
    return;
  }

  auto countUniqueExprNodes =
      [](const std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash>& exprByKey) {
        std::unordered_set<BoolExpr*> visited;
        std::vector<BoolExpr*> stack;
        for (const auto& [_, root] : exprByKey) {
          if (root != nullptr) {
            stack.push_back(root);
          }
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
      };

  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> resetSpecializedNextStateByKey;
  resetSpecializedNextStateByKey.reserve(model.stateBits.size());
  std::unordered_map<BoolExpr*, BoolExpr*> resetSubstitutionMemo;
  for (const auto& key : model.stateBits) {
    const auto nextStateIt = model.nextStateExprByStateKey.find(key);
    if (nextStateIt == model.nextStateExprByStateKey.end()) {
      continue;
    }
    resetSpecializedNextStateByKey.emplace(
        key,
        substituteBoolExprVariables(
            nextStateIt->second, resetAssignments, resetSubstitutionMemo));
  }

  constexpr size_t kMaxResetSpecializedExprNodesForInitInference = 50000;
  if (countUniqueExprNodes(resetSpecializedNextStateByKey) >
      kMaxResetSpecializedExprNodesForInitInference) {
    return;
  }

  auto collectReferencedStateVars = [](BoolExpr* expr) {
    std::unordered_set<size_t> referencedVars;
    if (expr == nullptr) {
      return referencedVars;
    }

    std::vector<BoolExpr*> stack = {expr};
    std::unordered_set<BoolExpr*> visited;
    while (!stack.empty()) {
      BoolExpr* current = stack.back();
      stack.pop_back();
      if (current == nullptr || !visited.insert(current).second) {
        continue;
      }
      if (current->getOp() == Op::VAR) {
        if (current->getId() >= 2) {
          referencedVars.insert(current->getId());
        }
        continue;
      }
      if (current->getLeft() != nullptr) {
        stack.push_back(current->getLeft());
      }
      if (current->getRight() != nullptr) {
        stack.push_back(current->getRight());
      }
    }
    return referencedVars;
  };

  std::unordered_map<size_t, SignalKey> stateKeyByVar;
  std::unordered_map<size_t, std::vector<SignalKey>> dependentStatesByVar;
  stateKeyByVar.reserve(model.stateBits.size());
  dependentStatesByVar.reserve(model.stateBits.size());
  for (const auto& key : model.stateBits) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      stateKeyByVar.emplace(varIt->second, key);
    }
  }
  for (const auto& key : model.stateBits) {
    const auto nextStateIt = resetSpecializedNextStateByKey.find(key);
    if (nextStateIt == resetSpecializedNextStateByKey.end()) {
      continue;
    }
    const auto referencedVars = collectReferencedStateVars(nextStateIt->second);
    for (const auto referencedVar : referencedVars) {
      if (stateKeyByVar.find(referencedVar) == stateKeyByVar.end()) {
        continue;
      }
      dependentStatesByVar[referencedVar].push_back(key);
    }
  }

  std::unordered_map<SignalKey, SignalKey, SignalKeyHash> complementedPartnerByKey;
  complementedPartnerByKey.reserve(model.complementedStateRelations.size() * 2);
  for (const auto& relation : model.complementedStateRelations) {
    complementedPartnerByKey.emplace(relation.primaryKey, relation.complementedKey);
    complementedPartnerByKey.emplace(relation.complementedKey, relation.primaryKey);
  }

  std::unordered_map<size_t, bool> assignments = resetAssignments;
  for (const auto& [key, value] : model.initialStateValueByKey) {
    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      assignments.emplace(varIt->second, value);
    }
  }

  std::deque<SignalKey> workQueue(model.stateBits.begin(), model.stateBits.end());
  auto recordKnownState = [&](const SignalKey& key, bool value) {
    const auto [it, inserted] = model.initialStateValueByKey.emplace(key, value);
    if (!inserted) {
      return;
    }

    const auto varIt = model.inputVarByKey.find(key);
    if (varIt != model.inputVarByKey.end()) {
      assignments[varIt->second] = value;
      const auto dependentIt = dependentStatesByVar.find(varIt->second);
      if (dependentIt != dependentStatesByVar.end()) {
        workQueue.insert(
            workQueue.end(),
            dependentIt->second.begin(),
            dependentIt->second.end());
      }
    }

    const auto partnerIt = complementedPartnerByKey.find(key);
    if (partnerIt != complementedPartnerByKey.end() &&
        model.initialStateValueByKey.find(partnerIt->second) ==
            model.initialStateValueByKey.end()) {
      workQueue.push_back(partnerIt->second);
    }
  };

  while (!workQueue.empty()) {
    const SignalKey key = workQueue.front();
    workQueue.pop_front();

    if (model.initialStateValueByKey.find(key) != model.initialStateValueByKey.end()) {
      const auto partnerIt = complementedPartnerByKey.find(key);
      if (partnerIt != complementedPartnerByKey.end() &&
          model.initialStateValueByKey.find(partnerIt->second) ==
              model.initialStateValueByKey.end()) {
        recordKnownState(partnerIt->second, !model.initialStateValueByKey.at(key));
      }
      continue;
    }

    const auto nextStateIt = resetSpecializedNextStateByKey.find(key);
    if (nextStateIt == resetSpecializedNextStateByKey.end()) {
      continue;
    }

    std::unordered_map<BoolExpr*, std::optional<bool>> memo;
    const auto resetValue = evaluateConstantUnderAssignmentsImplForTest(
        nextStateIt->second, assignments, memo);
    if (resetValue.has_value()) {
      recordKnownState(key, *resetValue);
    }
  }
}

std::optional<bool> getResetAssertionValueForTest(const std::string& displayName) {
  return getResetAssertionValueFromDisplayNameForTest(displayName);
}

namespace {

std::optional<bool> getResetAssertionValueFromDisplayNameForTest(
    const std::string& displayName) {
  const std::string normalized = normalizeSignalBaseNameForTest(displayName);
  std::vector<std::string> candidates = {normalized};
  if (hasSuffixForTest(normalized, "_I")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 2));
  }
  if (hasSuffixForTest(normalized, "_NI")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 1));
  }
  appendDomainPrefixedActiveLowResetCandidatesForTest(candidates);
  for (const auto& candidate : candidates) {
    if (isResetNameTokenForTest(candidate, "RESET") ||
        isResetNameTokenForTest(candidate, "RST")) {
      return true;
    }
    if (isResetNameTokenForTest(candidate, "RESET_N") ||
        isResetNameTokenForTest(candidate, "RESETN") ||
        isResetNameTokenForTest(candidate, "RESET_L") ||
        isResetNameTokenForTest(candidate, "RST_N") ||
        isResetNameTokenForTest(candidate, "RSTN") ||
        isResetNameTokenForTest(candidate, "RST_L")) {
      return false;
    }
  }
  return std::nullopt;
}
 
}  // namespace

std::unordered_map<SignalKey, bool, SignalKeyHash>
deriveResetBootstrapStateValuesForTest(
    const SequentialDesignModel& model,
    size_t cycles) {
  const auto resetAssignments = collectResetAssignmentsForTest(model);
  if (resetAssignments.empty() || cycles == 0) {
    return {};
  }

  std::unordered_map<SignalKey, bool, SignalKeyHash> knownStates =
      model.initialStateValueByKey;
  for (size_t step = 0; step < cycles; ++step) {
    std::unordered_map<size_t, bool> assignments = resetAssignments;
    for (const auto& [key, value] : knownStates) {
      const auto varIt = model.inputVarByKey.find(key);
      if (varIt != model.inputVarByKey.end()) {
        assignments.emplace(varIt->second, value);
      }
    }

    std::unordered_map<SignalKey, bool, SignalKeyHash> nextKnownStates;
    std::unordered_map<BoolExpr*, std::optional<bool>> memo;
    for (const auto& key : model.stateBits) {
      const auto value = evaluateConstantUnderAssignmentsImplForTest(
          model.nextStateExprByStateKey.at(key), assignments, memo);
      if (value.has_value()) {
        nextKnownStates.emplace(key, *value);
      }
    }
    knownStates = std::move(nextKnownStates);
  }

  return knownStates;
}

AlignedSignals filterStateEqualitiesByInitialValueForTest(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& candidateStates) {
  AlignedSignals anchoredStates;
  for (size_t i = 0; i < candidateStates.names.size(); ++i) {
    const auto initial0 = model0.initialStateValueByKey.find(candidateStates.keys0[i]);
    const auto initial1 = model1.initialStateValueByKey.find(candidateStates.keys1[i]);
    if (initial0 == model0.initialStateValueByKey.end() ||
        initial1 == model1.initialStateValueByKey.end() ||
        initial0->second != initial1->second) {
      continue;
    }

    anchoredStates.names.push_back(candidateStates.names[i]);
    anchoredStates.keys0.push_back(candidateStates.keys0[i]);
    anchoredStates.keys1.push_back(candidateStates.keys1[i]);
  }
  return anchoredStates;
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

std::string makeOneHotPdrFullFlowImplSource(const std::string& moduleName,
                                            size_t depth,
                                            bool reachableBad) {
  const size_t stateCount = depth + 1;
  std::ostringstream source;
  source << "module " << moduleName
         << "(input clk, input reset, output out);\n";
  for (size_t index = 0; index < stateCount; ++index) {
    source << "  reg s" << index << ";\n";
  }
  // Keep the parsed full-flow PDR fixture in one clocked process.  Newer
  // SystemVerilog frontend lowering can split independent procedural blocks in
  // a way that makes this tiny synthetic chain frontend-shape dependent, while
  // the intended SEC/PDR behavior is only the one-hot temporal chain below.
  source << "  always @(posedge clk) begin\n";
  source << "    if (reset) begin\n";
  for (size_t index = 0; index < stateCount; ++index) {
    source << "      s" << index << " <= "
           << (index == 0 ? "1'b1" : "1'b0") << ";\n";
  }
  source << "    end else begin\n";
  for (size_t index = 0; index < stateCount; ++index) {
    source << "      s" << index << " <= ";
    if (index == 0) {
      source << (reachableBad ? "1'b0" : "s0");
    } else if (reachableBad || index > 1) {
      source << "s" << (index - 1);
    } else {
      source << "1'b0";
    }
    source << ";\n";
  }
  source << "    end\n";
  source << "  end\n";
  source << "  assign out = s" << depth << ";\n";
  source << "endmodule\n";
  return source.str();
}

std::string makeOneHotPdrFullFlowReferenceSource(
    const std::string& moduleName) {
  std::ostringstream source;
  source << "module " << moduleName
         << "(input clk, input reset, output out);\n";
  source << "  assign out = 1'b0;\n";
  source << "endmodule\n";
  return source.str();
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
  const ScopedUnsetEnvVar resetShortcutDiag(
      "KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG");
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::Pdr);
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

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::KInduction);
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

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::Imc);
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");
  testing::internal::CaptureStderr();
  const auto result = strategy.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 1u);
  EXPECT_NE(
      stderrOutput.find("SEC diag: SEC IMC proven outputs: 1/1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("SEC diag: SEC IMC not proven output"),
      std::string::npos)
      << stderrOutput;
  ASSERT_TRUE(result.proofProgress.has_value());
  EXPECT_EQ(result.proofProgress->engineLabel, "IMC");
  EXPECT_EQ(result.proofProgress->provenOutputs, 1u);
  EXPECT_EQ(result.proofProgress->totalOutputs, 1u);
  EXPECT_TRUE(result.proofProgress->unprovenOutputs.empty());
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  // Without a cross-design state assumption, the inverted registered output is
  // first a concrete SEC mismatch after one transition.
  EXPECT_EQ(result.bound, 1u);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  // PDR can close the invariant before the visible output stage.
  EXPECT_LE(result.bound, 3u);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(3);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_LE(result.bound, 3u);
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
       ResetBootstrapHiddenShiftPipelineDoesNotCloseBelowDepth) {
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(1);

  // A one-step proof cannot justify equality of a hidden 20-stage shift chain
  // because SEC does not assume cross-design internal state correspondence.
  // Keep the result inconclusive until the caller supplies a sufficient KI
  // horizon.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetBootstrapLongEquivalentPipelineDoesNotCloseAtSmallK) {
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

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::KInduction);
  const auto result = strategy.run(3);

  // The removed startup fast path used internal cross-design state facts.
  // With strict top-output k-induction inputs only, this 12-stage pipe needs
  // a deeper caller horizon than k=3.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.bound, 3u);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(6);

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_LE(result.bound, 6u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantDoesNotRelateInternalStateWithoutReset) {
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.stateBits = {state0};
  model0.displayNameByKey.emplace(state0, "same_name_state[0]");
  model0.initialStateValueByKey.emplace(state0, false);

  SequentialDesignModel model1;
  model1.stateBits = {state1};
  model1.displayNameByKey.emplace(state1, "same_name_state[0]");
  model1.initialStateValueByKey.emplace(state1, false);

  const auto invariant = buildReachableStateInvariant(model0, model1);

  EXPECT_EQ(invariant.bootstrapCycles, 0u);
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
  model1.initialStateValueByKey.emplace(state1, true);

  const auto invariant = buildReachableStateInvariant(model0, model1);

  EXPECT_EQ(invariant.bootstrapCycles, 0u);
  EXPECT_TRUE(invariant.bootstrapValues0.empty());
  EXPECT_TRUE(invariant.bootstrapValues1.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantDerivesDesignLocalBootstrapValuesFromReset) {
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
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.inputVarByKey.emplace(state1, 5);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));

  const auto invariant = buildReachableStateInvariant(model0, model1);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantCanSkipDesignLocalBootstrapValueSweep) {
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
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.inputVarByKey.emplace(state1, 5);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));

  const auto invariant = buildReachableStateInvariant(
      model0, model1, /*deriveResetBootstrapStrengthening=*/false);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  EXPECT_TRUE(invariant.bootstrapValues0.empty());
  EXPECT_TRUE(invariant.bootstrapValues1.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantRecognizesInputSuffixedResetNames) {
  const SignalKey reset0 = makeSignalKey("reset0");
  const SignalKey reset1 = makeSignalKey("reset1");
  const SignalKey activeLowReset0 = makeSignalKey("activeLowReset0");
  const SignalKey activeLowReset1 = makeSignalKey("activeLowReset1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");
  const SignalKey lowState0 = makeSignalKey("lowState0");
  const SignalKey lowState1 = makeSignalKey("lowState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {reset0, activeLowReset0};
  model0.stateBits = {state0, lowState0};
  model0.inputVarByKey.emplace(reset0, 2);
  model0.inputVarByKey.emplace(activeLowReset0, 4);
  model0.inputVarByKey.emplace(state0, 6);
  model0.inputVarByKey.emplace(lowState0, 8);
  model0.displayNameByKey.emplace(reset0, "reset_i");
  model0.displayNameByKey.emplace(activeLowReset0, "rst_ni");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Not(BoolExpr::Var(2)));
  model0.nextStateExprByStateKey.emplace(lowState0, BoolExpr::Var(4));

  SequentialDesignModel model1;
  model1.environmentInputs = {reset1, activeLowReset1};
  model1.stateBits = {state1, lowState1};
  model1.inputVarByKey.emplace(reset1, 3);
  model1.inputVarByKey.emplace(activeLowReset1, 5);
  model1.inputVarByKey.emplace(state1, 7);
  model1.inputVarByKey.emplace(lowState1, 9);
  model1.displayNameByKey.emplace(reset1, "reset_i");
  model1.displayNameByKey.emplace(activeLowReset1, "rst_ni");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));
  model1.nextStateExprByStateKey.emplace(lowState1, BoolExpr::Var(5));

  const auto invariant = buildReachableStateInvariant(model0, model1);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.bootstrapValues0.size(), 2u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  EXPECT_FALSE(invariant.bootstrapValues0.at(lowState0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 2u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
  EXPECT_FALSE(invariant.bootstrapValues1.at(lowState1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantRecognizesDomainPrefixedActiveLowResetNames) {
  const SignalKey readReset0 = makeSignalKey("readReset0");
  const SignalKey writeReset1 = makeSignalKey("writeReset1");
  const SignalKey state0 = makeSignalKey("state0");
  const SignalKey state1 = makeSignalKey("state1");

  SequentialDesignModel model0;
  model0.environmentInputs = {readReset0};
  model0.stateBits = {state0};
  model0.inputVarByKey.emplace(readReset0, 2);
  model0.inputVarByKey.emplace(state0, 4);
  model0.displayNameByKey.emplace(readReset0, "rrst_n");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Var(2));

  SequentialDesignModel model1;
  model1.environmentInputs = {writeReset1};
  model1.stateBits = {state1};
  model1.inputVarByKey.emplace(writeReset1, 3);
  model1.inputVarByKey.emplace(state1, 5);
  model1.displayNameByKey.emplace(writeReset1, "wrst_n");
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Var(3));

  const auto invariant = buildReachableStateInvariant(model0, model1);

  // These resets are design-local active-low controls, not cross-design state
  // relations.
  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues0.at(state0));
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_FALSE(invariant.bootstrapValues1.at(state1));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantHandlesInputInSuffixAndMissingNextState) {
  const SignalKey reset0 = makeSignalKey("reset0");
  const SignalKey reset1 = makeSignalKey("reset1");
  const SignalKey driven0 = makeSignalKey("driven0");
  const SignalKey driven1 = makeSignalKey("driven1");
  const SignalKey missingNext0 = makeSignalKey("missingNext0");
  const SignalKey missingNext1 = makeSignalKey("missingNext1");

  SequentialDesignModel model0;
  model0.environmentInputs = {reset0};
  model0.stateBits = {driven0, missingNext0};
  model0.inputVarByKey.emplace(reset0, 2);
  model0.inputVarByKey.emplace(driven0, 4);
  model0.inputVarByKey.emplace(missingNext0, 6);
  model0.displayNameByKey.emplace(reset0, "reset_in");
  model0.nextStateExprByStateKey.emplace(driven0, BoolExpr::Var(2));

  SequentialDesignModel model1;
  model1.environmentInputs = {reset1};
  model1.stateBits = {driven1, missingNext1};
  model1.inputVarByKey.emplace(reset1, 3);
  model1.inputVarByKey.emplace(driven1, 5);
  model1.inputVarByKey.emplace(missingNext1, 7);
  model1.displayNameByKey.emplace(reset1, "reset_in");
  model1.nextStateExprByStateKey.emplace(driven1, BoolExpr::Var(3));

  const auto invariant = buildReachableStateInvariant(model0, model1);

  // States without a local next-state expression cannot receive a reset-derived
  // bootstrap value.
  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  ASSERT_EQ(invariant.bootstrapValues0.size(), 1u);
  EXPECT_TRUE(invariant.bootstrapValues0.at(driven0));
  EXPECT_EQ(invariant.bootstrapValues0.count(missingNext0), 0u);
  ASSERT_EQ(invariant.bootstrapValues1.size(), 1u);
  EXPECT_TRUE(invariant.bootstrapValues1.at(driven1));
  EXPECT_EQ(invariant.bootstrapValues1.count(missingNext1), 0u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       ReachableStateInvariantEvaluatesBootstrapValueOperatorsAndInvalidNodes) {
  const SignalKey rst0 = makeSignalKey("rst0");
  const SignalKey rst1 = makeSignalKey("rst1");
  const SignalKey const0 = makeSignalKey("const0");
  const SignalKey const1 = makeSignalKey("const1");
  const SignalKey xor0 = makeSignalKey("xor0");
  const SignalKey xor1 = makeSignalKey("xor1");
  const SignalKey diff0 = makeSignalKey("diff0");
  const SignalKey diff1 = makeSignalKey("diff1");
  const SignalKey invalid0 = makeSignalKey("invalid0");
  const SignalKey invalid1 = makeSignalKey("invalid1");
  BoolExpr invalidExpr0;
  BoolExpr invalidExpr1;

  SequentialDesignModel model0;
  model0.environmentInputs = {rst0};
  model0.stateBits = {const0, xor0, diff0, invalid0};
  model0.inputVarByKey.emplace(rst0, 2);
  model0.inputVarByKey.emplace(const0, 4);
  model0.inputVarByKey.emplace(xor0, 6);
  model0.inputVarByKey.emplace(diff0, 8);
  model0.inputVarByKey.emplace(invalid0, 10);
  model0.displayNameByKey.emplace(rst0, "rst");
  model0.nextStateExprByStateKey.emplace(const0, BoolExpr::createTrue());
  model0.nextStateExprByStateKey.emplace(
      xor0, BoolExpr::Xor(BoolExpr::Var(2), BoolExpr::createFalse()));
  model0.nextStateExprByStateKey.emplace(diff0, BoolExpr::Var(2));
  model0.nextStateExprByStateKey.emplace(invalid0, &invalidExpr0);

  SequentialDesignModel model1;
  model1.environmentInputs = {rst1};
  model1.stateBits = {const1, xor1, diff1, invalid1};
  model1.inputVarByKey.emplace(rst1, 3);
  model1.inputVarByKey.emplace(const1, 5);
  model1.inputVarByKey.emplace(xor1, 7);
  model1.inputVarByKey.emplace(diff1, 9);
  model1.inputVarByKey.emplace(invalid1, 11);
  model1.displayNameByKey.emplace(rst1, "rst");
  model1.nextStateExprByStateKey.emplace(const1, BoolExpr::createTrue());
  model1.nextStateExprByStateKey.emplace(
      xor1, BoolExpr::Xor(BoolExpr::Var(3), BoolExpr::createFalse()));
  model1.nextStateExprByStateKey.emplace(diff1, BoolExpr::Not(BoolExpr::Var(3)));
  model1.nextStateExprByStateKey.emplace(invalid1, &invalidExpr1);

  const auto invariant = buildReachableStateInvariant(model0, model1);

  EXPECT_EQ(invariant.bootstrapCycles, 3u);
  EXPECT_TRUE(invariant.bootstrapValues0.at(const0));
  EXPECT_TRUE(invariant.bootstrapValues1.at(const1));
  EXPECT_TRUE(invariant.bootstrapValues0.at(xor0));
  EXPECT_TRUE(invariant.bootstrapValues1.at(xor1));
  EXPECT_TRUE(invariant.bootstrapValues0.at(diff0));
  EXPECT_FALSE(invariant.bootstrapValues1.at(diff1));
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
       BaseCaseSolverConcreteDualRailBootstrapForbidsUnknownRailValues) {
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
  problem.bootstrapStateAssignments = {{mayBeOne, false}, {mayBeZero, true}};
  problem.dualRailStatePairs = {DualRailSymbolPair{mayBeOne, mayBeZero}};
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
  constexpr size_t kStateCount = 512;
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
    // Only the last target bit is impossible. The other bits are constants that
    // make the first bad cube wide, matching the ASIC case where the actual
    // predecessor surface was tiny but the blocked cube had many literals.
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
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kStateCount,
      /*preciseBadCubeStateLimit=*/kStateCount);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("(!x" + std::to_string(constantFalseSymbol) + ")\n"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesPredecessorCoresForProjectedWideBlockedCubes) {
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

  // BlackParrot sampling showed projected-frame stages learning many adjacent
  // wide blockers. The predecessor-core path is still sound in projected mode:
  // if a weaker frame query cannot reach the reduced cube, the complete frame
  // cannot reach it either.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kStateCount,
      /*preciseBadCubeStateLimit=*/kStateCount,
      /*useExactFrameClauses=*/false);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("predecessor core target=96->1 source_level=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesPredecessorCoresForMediumSupportWideBlockedCubes) {
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

  // The cheap 8-literal seed is reachable, while the full wide cube is blocked
  // by a medium-sized transition surface. BlackParrot produced thousands of
  // such 68/88-literal blockers, so wide cubes should try the predecessor-core
  // oracle before bounded chunk dropping even when their support is not huge.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kStateCount,
      /*preciseBadCubeStateLimit=*/kStateCount,
      /*useExactFrameClauses=*/false);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("predecessor core target=96->1 source_level=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesPredecessorCoresForMediumHighSupportBlockedCubes) {
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
    // The first four target bits keep the cheap seed reachable. The remaining
    // target bits depend on a broad support cone that is false in the startup
    // frontier, matching the measured AES 12-literal, 113-support level-zero
    // blockers that were too small for the old medium predecessor-core gate
    // but still expensive to learn one neighboring valuation at a time.
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
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kTargetStateCount + kSupportStateCount,
      /*preciseBadCubeStateLimit=*/kTargetStateCount,
      /*useExactFrameClauses=*/false);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("predecessor core target=12->1 source_level=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesPredecessorCoresForLocalBroadDualRailBlockedCubes) {
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
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kTargetStateCount + kSupportStateCount,
      /*preciseBadCubeStateLimit=*/kTargetStateCount,
      /*useExactFrameClauses=*/false);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("predecessor core target=12->1 source_level=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("skipped dual-rail predecessor core"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesCachedCoreForHugeBroadDualRailBlockedCubes) {
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
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/kTargetStateCount + kSupportStateCount,
      /*preciseBadCubeStateLimit=*/kTargetStateCount,
      /*useExactFrameClauses=*/false);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // The fresh predecessor-core oracle must still be skipped on this broad
  // dual-rail surface, but the already-run cached predecessor query now gives
  // us the failed-assumption core for free.
  EXPECT_NE(
      stderrOutput.find("skipped dual-rail predecessor core"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor cached core target=12->1 source_level=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("predecessor core target="), std::string::npos)
      << stderrOutput;
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
       PDREngineProjectsBadCubesToRelevantStateSupport) {
  auto problem = buildDocumentedBooleanPdrCounterexampleProblem();
  problem.state0Symbols.push_back(5);
  problem.allSymbols.push_back(5);
  problem.initialCondition = BoolExpr::And(
      problem.initialCondition, BoolExpr::Not(BoolExpr::Var(5)));
  problem.transitions0.emplace_back(
      5, BoolExpr::Xor(BoolExpr::Var(5), BoolExpr::Var(4)));

  const ScopedEnvVar secPdrTrace("KEPLER_SEC_PDR_TRACE", "1");
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
  EXPECT_EQ(badCubeTrace.find("x5"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesObservationOnlyFrontierWithoutExplicitInit) {
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

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesBadFormulaRepairOnResetObservationFrontier) {
  KInductionProblem problem;
  problem.inputSymbols = {2};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{2, true}};
  problem.state0Symbols = {3};
  problem.allSymbols = {2, 3};
  problem.totalStateCount = 1;
  problem.transitions0.emplace_back(3, BoolExpr::Var(3));
  problem.observedOutputExprs0 = {BoolExpr::Var(3)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"reset_observation_out"};
  problem.bad = BoolExpr::Var(3);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionBad = problem.bad;
  problem.inductionProperty = problem.property;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("refined projected counterexample with validated "
                        "bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDoesNotUseImmediateProofWhenFrameBudgetIsZero) {
  KInductionProblem problem;
  problem.state0Symbols = {2};
  problem.allSymbols = {2};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.usesDualRailStateEncoding = true;
  problem.totalStateCount = 13;
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);

  EXPECT_EQ(result.status, PDRStatus::Inconclusive);
  EXPECT_EQ(result.bound, 0u);
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
       PDREngineDualRailResetBootstrapUsesOriginalOutputSurface) {
  EXPECT_TRUE(detail::pdrResetBootstrapPrecheckTooLarge(
      /*usesDualRailStateEncoding=*/true,
      /*observedOutputCount=*/16,
      /*originalObservedOutputCount=*/99,
      /*transitionSources=*/8,
      /*transitionSourceLimit=*/1024,
      /*outputLimit=*/64));
  EXPECT_FALSE(detail::pdrResetBootstrapPrecheckTooLarge(
      /*usesDualRailStateEncoding=*/false,
      /*observedOutputCount=*/16,
      /*originalObservedOutputCount=*/99,
      /*transitionSources=*/8,
      /*transitionSourceLimit=*/1024,
      /*outputLimit=*/64));
  EXPECT_FALSE(detail::pdrResetBootstrapPrecheckTooLarge(
      /*usesDualRailStateEncoding=*/true,
      /*observedOutputCount=*/16,
      /*originalObservedOutputCount=*/16,
      /*transitionSources=*/8,
      /*transitionSourceLimit=*/1024,
      /*outputLimit=*/64));
  EXPECT_FALSE(detail::pdrResetBootstrapPrecheckTooLarge(
      /*usesDualRailStateEncoding=*/true,
      /*observedOutputCount=*/99,
      /*originalObservedOutputCount=*/99,
      /*transitionSources=*/4224,
      /*transitionSourceLimit=*/8192));
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
       PdrOrderedClauseFingerprintSeparatesProjectedRetries) {
  struct ClauseLiteralForFingerprintTest {
    size_t symbol = 0;
    bool positive = false;
  };
  using ClauseForFingerprintTest =
      std::vector<ClauseLiteralForFingerprintTest>;
  const std::vector<ClauseForFingerprintTest> emptyClauses;
  const std::vector<ClauseForFingerprintTest> retryClauses = {
      {{2, true}, {3, false}}, {{5, true}}};
  const std::vector<ClauseForFingerprintTest> sameClauses = {
      {{2, true}, {3, false}}, {{5, true}}};
  const std::vector<ClauseForFingerprintTest> reorderedClauses = {
      {{5, true}}, {{2, true}, {3, false}}};
  const std::vector<ClauseForFingerprintTest> polarityChangedClauses = {
      {{2, true}, {3, true}}, {{5, true}}};

  // Projected PDR retries use this fingerprint in the predecessor result-cache
  // key. Empty means "no extra retry clauses"; every real local refinement must
  // remain distinct so cached UNSAT/SAT answers cannot leak across retries.
  EXPECT_EQ(detail::pdrOrderedClauseFingerprint(emptyClauses), 0u);
  EXPECT_EQ(
      detail::pdrOrderedClauseFingerprint(retryClauses),
      detail::pdrOrderedClauseFingerprint(sameClauses));
  EXPECT_NE(
      detail::pdrOrderedClauseFingerprint(retryClauses),
      detail::pdrOrderedClauseFingerprint(reorderedClauses));
  EXPECT_NE(
      detail::pdrOrderedClauseFingerprint(retryClauses),
      detail::pdrOrderedClauseFingerprint(polarityChangedClauses));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrMergeSortedSymbolVectorsKeepsStableSurfaceSortedUnique) {
  const std::vector<size_t> stableSurface = {2, 4, 8, 16};
  const std::vector<size_t> localSymbols = {3, 4, 16, 32};

  const std::vector<size_t> merged =
      detail::mergeSortedPdrSymbolVectors(stableSurface, localSymbols);

  // Local dual-rail predecessor caching uses this instead of rebuilding a
  // large unordered_set and sorting it on every query. The result must still be
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

  // Local dual-rail predecessor caches use this to keep one solver alive when
  // neighboring target cubes add a few non-state support symbols. The widened
  // surface must stay sorted/unique because it becomes the SAT variable list.
  EXPECT_TRUE(
      detail::widenSortedPdrSymbolSurface(stableSurface, {3, 4, 9}));
  EXPECT_EQ(stableSurface, (std::vector<size_t>{2, 3, 4, 8, 9}));

  EXPECT_FALSE(
      detail::widenSortedPdrSymbolSurface(stableSurface, {2, 3, 9}));
  EXPECT_EQ(stableSurface, (std::vector<size_t>{2, 3, 4, 8, 9}));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrStableLocalPredecessorCacheSurfaceIsLevelZeroOnly) {
  // The stable solver surface is only a cache optimization for startup/frontier
  // local leaves. Level-1+ predecessor retries must stay on their exact local
  // symbol surface so Swerv does not spend budgets on broad SAT instances.
  EXPECT_TRUE(
      detail::shouldUseStableLocalPredecessorCacheSurface(true, true, 0));
  EXPECT_FALSE(
      detail::shouldUseStableLocalPredecessorCacheSurface(true, true, 1));
  EXPECT_FALSE(
      detail::shouldUseStableLocalPredecessorCacheSurface(true, false, 0));
  EXPECT_FALSE(
      detail::shouldUseStableLocalPredecessorCacheSurface(false, true, 0));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrPredecessorUnsatCoreSharingUsesBaseContextOnly) {
  // A predecessor UNSAT core can be reused for stronger target cubes only when
  // it came from the monotonic base frame context. Selector assumptions and
  // projected retry clauses stay target-local.
  EXPECT_TRUE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/0,
      /*extraFrameFingerprint=*/0,
      /*excludeTargetOnCurrentFrame=*/false));
  EXPECT_FALSE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/7,
      /*extraFrameFingerprint=*/0,
      /*excludeTargetOnCurrentFrame=*/false));
  EXPECT_FALSE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/0,
      /*extraFrameFingerprint=*/11,
      /*excludeTargetOnCurrentFrame=*/false));
  EXPECT_FALSE(detail::shouldSharePredecessorUnsatCore(
      /*frameFingerprint=*/0,
      /*extraFrameFingerprint=*/0,
      /*excludeTargetOnCurrentFrame=*/true));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrLargeDualRailPredecessorResetFrontierRepairModePolicy) {
  constexpr size_t kDefaultSupportLimit = 8192;
  constexpr size_t kLocalSupportLimit = 16 * 1024;
  EXPECT_TRUE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/false,
      /*observedOutputCount=*/1,
      /*level=*/0,
      /*targetCubeSize=*/32,
      /*transitionSupportSize=*/4096,
      /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_TRUE(detail::shouldUseOneShotLargeDualRailResetFrontierPredecessor(
      /*hasLargeDualRailResetFrontierSurface=*/true,
      /*hasLocalDualRailLeafRepairSurface=*/false));
  EXPECT_FALSE(detail::shouldUseOneShotLargeDualRailResetFrontierPredecessor(
      /*hasLargeDualRailResetFrontierSurface=*/true,
      /*hasLocalDualRailLeafRepairSurface=*/true));
  EXPECT_FALSE(detail::shouldUseOneShotLargeDualRailResetFrontierPredecessor(
      /*hasLargeDualRailResetFrontierSurface=*/false,
      /*hasLocalDualRailLeafRepairSurface=*/false));
  EXPECT_FALSE(detail::shouldRunLargeDualRailResetFrontierQuery(
      /*resetFrontierQueryAllowed=*/true,
      /*hasLargeDualRailResetFrontierSurface=*/true,
      /*hasLocalDualRailLeafRepairSurface=*/false));
  EXPECT_TRUE(detail::shouldRunLargeDualRailResetFrontierQuery(
      /*resetFrontierQueryAllowed=*/true,
      /*hasLargeDualRailResetFrontierSurface=*/true,
      /*hasLocalDualRailLeafRepairSurface=*/true));
  EXPECT_TRUE(detail::shouldRunLargeDualRailResetFrontierQuery(
      /*resetFrontierQueryAllowed=*/true,
      /*hasLargeDualRailResetFrontierSurface=*/false,
      /*hasLocalDualRailLeafRepairSurface=*/false));
  EXPECT_FALSE(detail::shouldRunLargeDualRailResetFrontierQuery(
      /*resetFrontierQueryAllowed=*/false,
      /*hasLargeDualRailResetFrontierSurface=*/false,
      /*hasLocalDualRailLeafRepairSurface=*/false));
  EXPECT_TRUE(
      detail::shouldPrecheckLargeDualRailPredecessorWithResetFrontier(
          /*usesDualRailStateEncoding=*/true,
          /*exactResetFrontierChecksEnabled=*/false,
          /*observedOutputCount=*/1,
          /*level=*/0,
          /*targetCubeSize=*/32,
          /*transitionSupportSize=*/4096,
          /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(
      detail::shouldPrecheckLargeDualRailPredecessorWithResetFrontier(
          /*usesDualRailStateEncoding=*/true,
          /*exactResetFrontierChecksEnabled=*/false,
          /*observedOutputCount=*/1,
          /*level=*/0,
          /*targetCubeSize=*/16,
          /*transitionSupportSize=*/4096,
          /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(
      detail::shouldPrecheckLargeDualRailPredecessorWithResetFrontier(
          /*usesDualRailStateEncoding=*/true,
          /*exactResetFrontierChecksEnabled=*/false,
          /*observedOutputCount=*/1,
          /*level=*/0,
          /*targetCubeSize=*/32,
          /*transitionSupportSize=*/3999,
          /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/true,
      /*observedOutputCount=*/1,
      /*level=*/0,
      /*targetCubeSize=*/32,
      /*transitionSupportSize=*/4096,
      /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/false,
      /*observedOutputCount=*/2,
      /*level=*/0,
      /*targetCubeSize=*/32,
      /*transitionSupportSize=*/4096,
      /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/false,
      /*observedOutputCount=*/1,
      /*level=*/1,
      /*targetCubeSize=*/32,
      /*transitionSupportSize=*/4096,
      /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/false,
      /*observedOutputCount=*/1,
      /*level=*/0,
      /*targetCubeSize=*/33,
      /*transitionSupportSize=*/4096,
      /*exactResetPrecheckSupportLimit=*/8192));
  EXPECT_FALSE(detail::shouldRetryLargeDualRailPredecessorWithResetFrontier(
      /*usesDualRailStateEncoding=*/true,
      /*exactResetFrontierChecksEnabled=*/false,
      /*observedOutputCount=*/1,
      /*level=*/0,
      /*targetCubeSize=*/32,
      /*transitionSupportSize=*/8193,
      /*exactResetPrecheckSupportLimit=*/8192));
  // Single-output local leaves can run the exact reset precheck before the
  // ordinary predecessor SAT query for slightly wider 8k-16k support cones.
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/1,
                /*level=*/0,
                /*targetCubeSize=*/32,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kLocalSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/1,
                /*level=*/0,
                /*targetCubeSize=*/27,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kDefaultSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/1,
                /*level=*/0,
                /*targetCubeSize=*/33,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kDefaultSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/2,
                /*level=*/0,
                /*targetCubeSize=*/32,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kDefaultSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/1,
                /*level=*/1,
                /*targetCubeSize=*/32,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kDefaultSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/false,
                /*observedOutputCount=*/1,
                /*level=*/0,
                /*targetCubeSize=*/32,
                kDefaultSupportLimit,
                kLocalSupportLimit),
            kDefaultSupportLimit);
  EXPECT_EQ(detail::effectiveLocalDualRailExactResetPrecheckSupportLimit(
                /*hasLocalDualRailLeafRepairSurface=*/true,
                /*observedOutputCount=*/1,
                /*level=*/0,
                /*targetCubeSize=*/32,
                /*configuredSupportLimit=*/0,
                kLocalSupportLimit),
            0);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrExactResetPredecessorSiblingSeedingCoversResidualCubeSize) {
  // Swerv residual dual-rail leaves produce 32-literal cubes whose exact reset
  // proof often minimizes to one singleton.  Seeding sibling singletons from
  // that same reset context avoids rediscovering the bus one full cube at a
  // time, while 33+ literal broad cubes stay outside this bounded shortcut.
  EXPECT_TRUE(
      detail::shouldSeedExactResetPredecessorSiblingCores(
          /*cubeSize=*/32,
          /*knownCoreSize=*/1));
  EXPECT_FALSE(
      detail::shouldSeedExactResetPredecessorSiblingCores(
          /*cubeSize=*/33,
          /*knownCoreSize=*/1));
  EXPECT_FALSE(
      detail::shouldSeedExactResetPredecessorSiblingCores(
          /*cubeSize=*/32,
          /*knownCoreSize=*/2));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrResidualDualRailPredecessorBudgetCoversLocalLeafShape) {
  // Swerv final dual-rail leaves can produce 28-32 literal residual targets
  // with a local 9k-15k symbol solver surface. Those are not broad batches, so
  // they should use the restored residual predecessor budget instead of the
  // lightweight batch retry cap.
  EXPECT_TRUE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/1, /*level=*/0,
          /*targetCubeSize=*/32,
          /*solverSymbolCount=*/16 * 1024));
  EXPECT_TRUE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/1, /*level=*/1,
          /*targetCubeSize=*/16,
          /*solverSymbolCount=*/8192));

  EXPECT_FALSE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/2, /*level=*/0,
          /*targetCubeSize=*/32,
          /*solverSymbolCount=*/16 * 1024));
  EXPECT_FALSE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/1, /*level=*/1,
          /*targetCubeSize=*/32,
          /*solverSymbolCount=*/16 * 1024));
  EXPECT_FALSE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/1, /*level=*/0,
          /*targetCubeSize=*/33,
          /*solverSymbolCount=*/16 * 1024));
  EXPECT_FALSE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          true, /*observedOutputCount=*/1, /*level=*/0,
          /*targetCubeSize=*/32,
          /*solverSymbolCount=*/16 * 1024 + 1));
  EXPECT_FALSE(
      detail::shouldUseResidualDualRailPredecessorBudget(
          false, /*observedOutputCount=*/1, /*level=*/0,
          /*targetCubeSize=*/32,
          /*solverSymbolCount=*/16 * 1024));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrBroadDualRailResidualSurfaceExcludesAesSizedLeaf) {
  constexpr size_t kDualRailMediumOutputLimit = 384;

  // AES-sized designs become one-output residual leaves after splitting, but
  // they must keep the 376a017 partial reset-conflict and fresh predecessor
  // route instead of broad-output residual shortcuts.
  EXPECT_FALSE(
      detail::isBroadDualRailResidualOutputSurface(
          /*usesDualRailStateEncoding=*/true,
          /*observedOutputCount=*/1,
          /*originalObservedOutputCount=*/129,
          kDualRailMediumOutputLimit));
  EXPECT_TRUE(
      detail::isBroadDualRailResidualOutputSurface(
          /*usesDualRailStateEncoding=*/true,
          /*observedOutputCount=*/1,
          /*originalObservedOutputCount=*/1266,
          kDualRailMediumOutputLimit));
  EXPECT_FALSE(
      detail::isBroadDualRailResidualOutputSurface(
          /*usesDualRailStateEncoding=*/true,
          /*observedOutputCount=*/128,
          /*originalObservedOutputCount=*/1266,
          kDualRailMediumOutputLimit));
  EXPECT_FALSE(
      detail::isBroadDualRailResidualOutputSurface(
          /*usesDualRailStateEncoding=*/false,
          /*observedOutputCount=*/1,
          /*originalObservedOutputCount=*/1266,
          kDualRailMediumOutputLimit));
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDeterministicCubeOrderingSortsResetCoreCandidates) {
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
       PdrDualRailFrameZeroValidationDefersHugeRailStateSurface) {
  // RISC-V HS has a smaller 99-output surface.  Keep it on the existing exact
  // validation path even if the rail-state surface is large, because the
  // Ariane runtime fix must not weaken that previous coverage fix.
  EXPECT_FALSE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/99,
      /*railStateSymbolSurface=*/2000000));
  EXPECT_FALSE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/133,
      /*railStateSymbolSurface=*/2000000));
  EXPECT_FALSE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/278,
      /*railStateSymbolSurface=*/1000000));

  // Ariane136 has only a mid-wide output bus, but compact dual-rail extraction
  // expands the rail state into a million-scale surface.  That shape should
  // enter PDR directly instead of spending minutes in the pre-PDR frame-0
  // validation pass.
  EXPECT_TRUE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/278,
      /*railStateSymbolSurface=*/1000001));

  // Dynamic-node has a mid-wide 331-output surface.  It should keep the exact
  // validation path unless compact extraction also creates the huge rail-state
  // shape seen in Ariane.
  EXPECT_FALSE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/331,
      /*railStateSymbolSurface=*/1000000));
  EXPECT_TRUE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/331,
      /*railStateSymbolSurface=*/1000001));

  EXPECT_TRUE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/385,
      /*railStateSymbolSurface=*/1));

  // BlackParrot-style wide output surfaces were already deferred by the old
  // output-count rule.  The new state-size rule should not change that behavior.
  EXPECT_TRUE(detail::shouldDeferPdrDualRailFrameZeroValidation(
      /*observedOutputSurface=*/598,
      /*railStateSymbolSurface=*/1));
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
       PDREngineProvesEquivalentExactlyAtThreeFrames) {
  const auto problem = buildLinearChainSecProblem(4);

  // This is an engine-regression check for the current binary-chain model and
  // current clause-generalization behavior.  It is not a portable "classic PDR
  // must prove safe exactly at k=3" theorem: safe IC3/PDR proofs may converge
  // earlier whenever a stronger inductive invariant is learned.
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(result.bound, 3u);
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

    EXPECT_EQ(earlyResult.status, PDRStatus::Inconclusive)
        << "badDepth=" << badDepth;
    EXPECT_EQ(earlyResult.bound, badDepth - 1) << "badDepth=" << badDepth;

    PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
    const auto result = engine.run(badDepth);

    EXPECT_EQ(result.status, PDRStatus::Different)
        << "badDepth=" << badDepth;
    EXPECT_EQ(result.bound, badDepth) << "badDepth=" << badDepth;
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrFullFlowProvesParsedVerilogSafeChainsWithinDepthsThreeFourFive) {
  for (const size_t proofDepth : {size_t{3}, size_t{4}, size_t{5}}) {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* library = NLLibrary::create(db, NLName("LIB"));
    const std::string suffix = std::to_string(proofDepth);
    auto* impl = loadSystemVerilogTopFromSource(
        library,
        "pdr_full_safe_impl_k" + suffix,
        makeOneHotPdrFullFlowImplSource(
            "pdr_full_safe_impl_k" + suffix,
            proofDepth,
            /*reachableBad=*/false));
    auto* reference = loadSystemVerilogTopFromSource(
        library,
        "pdr_full_safe_ref_k" + suffix,
        makeOneHotPdrFullFlowReferenceSource("pdr_full_safe_ref_k" + suffix));

    auto strategy = makeBinarySecStrategy(impl, reference, SecEngine::Pdr);
    const auto result = strategy.run(proofDepth);

    // Full-flow parsed-Verilog safe proofs have the same limitation as direct
    // PDR safe proofs: exact proof depth is implementation-dependent because
    // PDR may generalize to an invariant before the requested depth.
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent)
        << "proofDepth=" << proofDepth << " reason=" << result.reason;
    EXPECT_LE(result.bound, proofDepth) << "proofDepth=" << proofDepth;
    EXPECT_EQ(result.coveredOutputs, 1u) << "proofDepth=" << proofDepth;
    EXPECT_EQ(result.totalOutputs, 1u) << "proofDepth=" << proofDepth;

    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrFullFlowFindsParsedVerilogCounterexamplesAtDepthsThreeFourFive) {
  for (const size_t badDepth : {size_t{3}, size_t{4}, size_t{5}}) {
    NLUniverse::create();
    auto* db = NLDB::create(NLUniverse::get());
    auto* library = NLLibrary::create(db, NLName("LIB"));
    const std::string suffix = std::to_string(badDepth);
    auto* impl = loadSystemVerilogTopFromSource(
        library,
        "pdr_full_diff_impl_k" + suffix,
        makeOneHotPdrFullFlowImplSource(
            "pdr_full_diff_impl_k" + suffix,
            badDepth,
            /*reachableBad=*/true));
    auto* reference = loadSystemVerilogTopFromSource(
        library,
        "pdr_full_diff_ref_k" + suffix,
        makeOneHotPdrFullFlowReferenceSource("pdr_full_diff_ref_k" + suffix));

    // Exact early-depth PDR behavior is covered above by the direct PDREngine
    // test.  This full-flow test verifies the Verilog parser, SEC miter build,
    // and selected PDR engine agree on the reachable counterexample.
    auto strategy = makeBinarySecStrategy(impl, reference, SecEngine::Pdr);
    const auto result = strategy.run(badDepth);

    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Different)
        << "badDepth=" << badDepth << " reason=" << result.reason;
    EXPECT_EQ(result.bound, badDepth) << "badDepth=" << badDepth;
    EXPECT_EQ(result.coveredOutputs, 1u) << "badDepth=" << badDepth;
    EXPECT_EQ(result.totalOutputs, 1u) << "badDepth=" << badDepth;

    naja::DNL::destroy();
    NLUniverse::get()->destroy();
    KEPLER_FORMAL::Tree2BoolExpr::iso2boolExpr_.clear();
    KEPLER_FORMAL::BoolExprCache::destroy();
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
       DualRailConcreteBootstrapStateForbidsUnknownRailValues) {
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

  // Complete bootstrap assignments describe a concrete Boolean rail value.
  // Induction may forbid the synthetic "unknown" rail state (both rails true)
  // without assuming any relation between the two compared designs.
  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::Proved);

  problem.deferBaseCaseChecks = true;
  // Deferring the local base check for an output slice must not widen the
  // induction step beyond the same concrete dual-rail state domain.
  EXPECT_EQ(
      proveByInductionStatus(
          problem,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          1,
          std::nullopt),
      InductionProofStatus::Proved);
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
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
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

  // Dynamic-node has 331 observed outputs. Treat it as a medium-wide PDR
  // surface so frame-0 validation seeds the exact reset/bootstrap facts instead
  // of falling into all-output abstract cube repair.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, kOutputCount);
  EXPECT_EQ(result.totalOutputs, kOutputCount);
  EXPECT_EQ(
      stderrOutput.find("skipped dual-rail frame-0 validation"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailDefersVeryWideEquivalentValidation) {
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

  // A very wide dual-rail batch is already closed by PDR here.  Do not rebuild
  // the broad concrete BMC validator after the proof, because that is the
  // BlackParrot runtime wall this path is meant to avoid.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, kOutputCount);
  EXPECT_EQ(result.totalOutputs, kOutputCount);
  EXPECT_NE(
      stderrOutput.find("deferred wide dual-rail equivalent validation outputs=385"),
      std::string::npos)
      << stderrOutput;
}

KInductionProblem makeDualRailResetFrontierGuardProblemForTest(
    size_t railPairs,
    size_t transitionSources) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.totalStateCount = railPairs;
  problem.dualRailStatePairs.reserve(railPairs);
  for (size_t index = 0; index < railPairs; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{index * 2, index * 2 + 1});
  }
  problem.transitions0.reserve(transitionSources);
  for (size_t index = 0; index < transitionSources; ++index) {
    const size_t symbol = 100000 + index;
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  return problem;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierAllowsIbexSizedSurface) {
  // Nangate45 Ibex needs exact reset-frontier repair at this rail/transition
  // scale; otherwise dual-rail PDR reports abstract init-reaching roots.
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/7748,
          /*transitionSources=*/15496);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierAllowsRiscvSizedMediumOutputSurface) {
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/2112,
          /*transitionSources=*/4224);
  while (problem.observedOutputExprs0.size() < 99) {
    problem.observedOutputExprs0.push_back(BoolExpr::Var(7));
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
  }
  problem.originalObservedOutputCount = 99;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // sky130hs_riscv32i has a medium 99-output dual-rail bus surface.  Keeping it
  // eligible for exact reset-frontier repair avoids exhausting ordinary PDR
  // bad-cube budgets on its reset-unanchored datapath buses.
  EXPECT_EQ(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierAllowsDynamicNodeSizedMediumOutputSurface) {
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/9028,
          /*transitionSources=*/18056);
  while (problem.observedOutputExprs0.size() < 331) {
    problem.observedOutputExprs0.push_back(BoolExpr::Var(7));
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
  }
  problem.originalObservedOutputCount = 331;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Nangate45 dynamic-node sits inside the same medium-wide reset-frontier
  // envelope as RISC-V by output count, but has a much larger rail surface.
  EXPECT_EQ(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierBlocksAesSizedMediumOutputSurface) {
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/1124,
          /*transitionSources=*/2248);
  while (problem.observedOutputExprs0.size() < 129) {
    problem.observedOutputExprs0.push_back(BoolExpr::Var(7));
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
  }
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES-sized residual leaves are medium-width by original output count but
  // too small by rail-state surface to justify exact reset-frontier context.
  // Keep them on the lower-memory predecessor/reset-conflict path from the
  // known-good 376a017 behavior.
  EXPECT_NE(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
  EXPECT_NE(stderrOutput.find("medium_state_min=4096"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierStillBlocksTooWideOutputSurface) {
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/1124,
          /*transitionSources=*/2248);
  while (problem.observedOutputExprs0.size() < 385) {
    problem.observedOutputExprs0.push_back(BoolExpr::Var(7));
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
  }
  problem.originalObservedOutputCount = 385;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The broad-output guard still prevents SoC-scale dual-rail runs from
  // rebuilding exact reset-prefix solvers for every residual leaf.
  EXPECT_NE(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
  EXPECT_NE(stderrOutput.find("original_outputs=385"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PdrDualRailExactResetFrontierGuardCountsRailSymbols) {
  KInductionProblem problem =
      makeDualRailResetFrontierGuardProblemForTest(
          /*railPairs=*/10001,
          /*transitionSources=*/0);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The PDR reset-frontier repair pays for both rails.  Guard on rail symbols
  // so large dual-rail SEC runs do not re-enter the exact reset-prefix
  // validation wall.
  EXPECT_NE(
      stderrOutput.find(
          "exact reset-frontier checks disabled for large dual-rail problem"),
      std::string::npos);
  EXPECT_NE(stderrOutput.find("rail_state_symbols=20002"),
            std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesConservativeFrontierWhenResetBmcSkipsEmptyDualRailSlice) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.originalObservedOutputCount = 385;
  problem.observedOutputExprs0 = {BoolExpr::Var(2)};
  problem.observedOutputExprs1 = {BoolExpr::Var(2)};
  problem.property = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Empty pruned reset-bootstrap slices should still let PDR prove properties
  // that hold over all states.  Returning inconclusive before the PDR loop
  // leaves harmless one-output dual-rail leaves uncovered.
  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("skipped dual-rail reset-bootstrap BMC precheck"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("max frame budget exhausted"),
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

  // IMC must enter its own interpolation engine directly.  The old dual-rail
  // residual helper batched/skipped outputs before IMCEngine saw the problem;
  // strict IMC may return inconclusive on this large state-dependent surface.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_NE(
      stderrOutput.find("SEC diag: entering imc engine"),
      std::string::npos);
  EXPECT_EQ(stderrOutput.find("dual-rail imc"), std::string::npos);
  EXPECT_EQ(stderrOutput.find("k-induction problem"), std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       RunExtractedModelsPdrDualRailRemainsInconclusiveWhenNoOutputIsCovered) {
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

  const ScopedEnvVar pdrQueryBudget(
      "KEPLER_SEC_PDR_DUAL_RAIL_PROJECTED_QUERY_BUDGET", "1");
  const ScopedEnvVar pdrFinalBudget(
      "KEPLER_SEC_PDR_DUAL_RAIL_FINAL_QUERY_BUDGET", "1");
  const ScopedEnvVar pdrClosureLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_BATCH_CLOSURE_LIMIT", "1");
  const ScopedEnvVar predecessorDecisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
  const ScopedEnvVar secDiag("KEPLER_SEC_DIAG", "1");

  testing::internal::CaptureStderr();
  SequentialEquivalenceStrategy strategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::Pdr,
      SecEncoding::DualRailSteady);
  const auto result = strategy.runExtractedModels(model0, model1, 1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // PDR can be resource-limited on every resetless rail output.  It must stay
  // inconclusive instead of reporting a vacuous zero-output equivalence or
  // invoking another SEC engine behind the selected mode.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.coveredOutputs, 0u);
  EXPECT_EQ(result.totalOutputs, kOutputs);
  EXPECT_NE(
      result.reason.find("Dual-rail PDR exhausted repair/projection"),
      std::string::npos);
  EXPECT_EQ(stderrOutput.find("trying k-induction"), std::string::npos);

  {
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
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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
       RunExtractedModelsDualRailCoversMatchingResetlessResidualWithoutPdrFallback) {
  const SignalKey good = makeSignalKey("dualRailResetlessGood");
  const SignalKey out = makeSignalKey("dualRailResetlessOut");
  const SignalKey rst = makeSignalKey("dualRailResetlessRst");
  const SignalKey data = makeSignalKey("dualRailResetlessData");
  const SignalKey state0 = makeSignalKey("dualRailResetlessState0");
  const SignalKey state1 = makeSignalKey("dualRailResetlessState1");

  SequentialDesignModel model0;
  model0.environmentInputs = {rst, data};
  model0.stateBits = {state0};
  model0.allObservedOutputs = {good, out};
  model0.observedOutputs = {good, out};
  model0.inputVarByKey.emplace(rst, 4);
  model0.inputVarByKey.emplace(data, 6);
  model0.inputVarByKey.emplace(state0, 2);
  model0.displayNameByKey.emplace(rst, "rst");
  model0.displayNameByKey.emplace(data, "data[0]");
  model0.displayNameByKey.emplace(good, "good[0]");
  model0.displayNameByKey.emplace(out, "resetless_out[0]");
  model0.displayNameByKey.emplace(state0, "left_state_q[0]");
  model0.nextStateExprByStateKey.emplace(state0, BoolExpr::Var(2));
  model0.observedOutputExprByKey.emplace(good, BoolExpr::Var(6));
  model0.observedOutputExprByKey.emplace(out, BoolExpr::Var(2));

  SequentialDesignModel model1;
  model1.environmentInputs = {rst, data};
  model1.stateBits = {state1};
  model1.allObservedOutputs = {good, out};
  model1.observedOutputs = {good, out};
  model1.inputVarByKey.emplace(rst, 5);
  model1.inputVarByKey.emplace(data, 7);
  model1.inputVarByKey.emplace(state1, 3);
  model1.displayNameByKey.emplace(rst, "rst");
  model1.displayNameByKey.emplace(data, "data[0]");
  model1.displayNameByKey.emplace(good, "good[0]");
  model1.displayNameByKey.emplace(out, "resetless_out[0]");
  model1.displayNameByKey.emplace(state1, "right_state_q[0]");
  // Binary SEC cannot use a cross-design state equality here: one side holds
  // the resetless state while the other toggles it.  Dual-rail mode proves the
  // top-output rail equality directly with KI, without invoking PDR.
  model1.nextStateExprByStateKey.emplace(state1, BoolExpr::Not(BoolExpr::Var(3)));
  model1.observedOutputExprByKey.emplace(good, BoolExpr::Var(7));
  model1.observedOutputExprByKey.emplace(out, BoolExpr::Var(3));

  auto binaryStrategy = makeBinaryExtractedSecStrategy(SecEngine::KInduction);
  const auto binaryResult = binaryStrategy.runExtractedModels(model0, model1, 1);
  EXPECT_EQ(binaryResult.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(binaryResult.coveredOutputs, 1u);
  EXPECT_EQ(binaryResult.totalOutputs, 2u);
  ASSERT_EQ(binaryResult.resetUnanchoredSkippedOutputs.size(), 1u);
  EXPECT_NE(
      binaryResult.resetUnanchoredSkippedOutputs.front().find("resetless_out[0]"),
      std::string::npos);

  SequentialEquivalenceStrategy dualRailStrategy(
      nullptr,
      nullptr,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      SecEngine::KInduction,
      SecEncoding::DualRailSteady);
  const auto dualRailResult =
      dualRailStrategy.runExtractedModels(model0, model1, 2);

  EXPECT_EQ(dualRailResult.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(dualRailResult.coveredOutputs, 2u);
  EXPECT_EQ(dualRailResult.totalOutputs, 2u);
  EXPECT_TRUE(dualRailResult.skippedObservedOutputs.empty());
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
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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

    EXPECT_EQ(pdrResult.status, SequentialEquivalenceStatus::Equivalent);
    EXPECT_EQ(pdrResult.coveredOutputs, 3u);
    EXPECT_EQ(pdrResult.totalOutputs, 3u);
    EXPECT_TRUE(pdrResult.skippedObservedOutputs.empty());

    {
      const ScopedEnvVar pdrQueryBudget(
          "KEPLER_SEC_PDR_DUAL_RAIL_PROJECTED_QUERY_BUDGET", "1");
      const ScopedEnvVar pdrClosureLimit(
          "KEPLER_SEC_PDR_DUAL_RAIL_BATCH_CLOSURE_LIMIT", "1");
      const ScopedEnvVar predecessorDecisionLimit(
          "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
      const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
      testing::internal::CaptureStderr();
      SequentialEquivalenceStrategy budgetedPdrStrategy(
          nullptr,
          nullptr,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          SecEngine::Pdr,
          SecEncoding::DualRailSteady);
      const auto budgetedPdrResult =
          budgetedPdrStrategy.runExtractedModels(model0, model1, 1);
      const std::string stderrOutput =
          testing::internal::GetCapturedStderr();

      // Even with tight projected-query knobs, PDR must either prove its own
      // batches or stay in its own inconclusive path; it must not call KI as a
      // hidden fallback.
      EXPECT_EQ(
          budgetedPdrResult.status, SequentialEquivalenceStatus::Equivalent);
      // If the tightened budgets still let PDR prove all residual leaves, that
      // is the preferred result.  The regression guard is engine isolation, not
      // forcing PDR to stay partial.
      EXPECT_EQ(budgetedPdrResult.coveredOutputs, 3u);
      EXPECT_EQ(budgetedPdrResult.totalOutputs, 3u);
      EXPECT_TRUE(budgetedPdrResult.skippedObservedOutputs.empty());
      EXPECT_NE(stderrOutput.find("closure_limit=1"), std::string::npos);
      EXPECT_EQ(stderrOutput.find("trying k-induction"), std::string::npos);
    }

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
    EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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
       RunExtractedModelsPdrDualRailDoesNotReportResetlessStateLeafMismatch) {
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

  // The unrelated resetless dummies push the rail-state count above the small
  // PDR certificate fast path. The observable edit is driven by unanchored
  // internal state, so the direct dual-rail PDR path must not report it as a
  // concrete SEC counterexample without a public startup relation.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.bound, 1u);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_TRUE(result.reason.empty());
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
  EXPECT_NE(
      result.reason.find("wide_frame_zero_probe[0]"),
      std::string::npos);
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

  // Large dual-rail IMC residuals may use Craig/deferred proof work, but a
  // resetless state-dependent edit is not a concrete top-output mismatch unless
  // startup state is publicly anchored. Keep it inconclusive instead of
  // reporting a difference through an internal-state relation.
  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Inconclusive);
  EXPECT_EQ(result.bound, 4u);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 1u);
  EXPECT_EQ(result.reason.find("delayed_out[0]"), std::string::npos);
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
  EXPECT_EQ(dualRailImcResult.status, SequentialEquivalenceStatus::Equivalent);
  // Craig IMC now follows the paper-style loop, so the first valid fixed-point
  // proof is the k=1 interpolant containment check rather than the removed
  // immediate transition-closure shortcut.
  EXPECT_EQ(dualRailImcResult.bound, 1u);
  EXPECT_EQ(dualRailImcResult.coveredOutputs, kWideStartupRelationOutputs);
  EXPECT_EQ(dualRailImcResult.totalOutputs, kWideStartupRelationOutputs);
  EXPECT_TRUE(dualRailImcResult.skippedObservedOutputs.empty());
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
  // This output is purely combinational even though unrelated state nearby is
  // reset-unanchored.  The conservative state-dependent coverage filter must not
  // drop such top outputs.
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
       PDREngineDoesNotTreatPartialBootstrapSummaryAsExactInitialState) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  // The bootstrap summary can be partial: x is known at the post-reset
  // frontier, while y is only known by actually unrolling the reset transition.
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::createFalse());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(3);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  EXPECT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 0)
          .has_value());

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRejectsBootstrapPredecessorsOutsideConcreteResetImage) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  // The summary proves x=0 at the post-reset frontier but says nothing about y.
  // The concrete reset unroll also forces y=0.  Without the level-0 refinement,
  // PDR can invent the abstract post-reset state y=1 and use it to reach x'=1.
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(3)));
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
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
       PDREngineCanDeferExactResetFrontierChecksToCallerValidation) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(3)));
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1)
          .has_value());

  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/false,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);

  // Projected SEC stages can return this abstract trace quickly because the
  // top-level SEC strategy immediately validates every PDR difference with the
  // exact bounded base-case query above before accepting it.
  EXPECT_EQ(result.status, PDRStatus::Different);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailLocalF0SkipsResetFrontierPrecheckForMediumCube) {
  KInductionProblem problem;
  constexpr size_t reset = 100;
  problem.usesDualRailStateEncoding = true;
  std::vector<size_t> xs;
  std::vector<size_t> ys;
  xs.reserve(16);
  ys.reserve(16);
  for (size_t bit = 0; bit < 16; ++bit) {
    xs.push_back(2 + bit);
    ys.push_back(32 + bit);
  }
  problem.state0Symbols = xs;
  problem.state0Symbols.insert(
      problem.state0Symbols.end(), ys.begin(), ys.end());
  problem.inputSymbols = {reset};
  problem.allSymbols = problem.state0Symbols;
  problem.allSymbols.push_back(reset);
  problem.totalStateCount = problem.state0Symbols.size();
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  // The abstract bootstrap summary leaves every y bit unconstrained, but the
  // concrete reset unroll below forces all y bits to 0.  A local dual-rail F0
  // predecessor query can therefore invent the all-ones y vector. This medium
  // cube is still cheaper to try through ordinary PDR first; the early exact
  // reset-frontier repair is reserved for larger, high-support residual cubes.
  for (const size_t x : xs) {
    problem.bootstrapStateAssignments.emplace_back(x, false);
  }
  for (size_t bit = 0; bit < xs.size(); ++bit) {
    problem.transitions0.emplace_back(
        xs[bit],
        BoolExpr::And(
            BoolExpr::Not(BoolExpr::Var(reset)), BoolExpr::Var(ys[bit])));
    problem.transitions0.emplace_back(ys[bit], BoolExpr::createFalse());
  }
  BoolExpr* bad = BoolExpr::createTrue();
  for (const size_t x : xs) {
    bad = BoolExpr::And(bad, BoolExpr::Var(x));
  }
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {problem.bad};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  problem.observedOutputNames = {"dual_rail_local_reset_frontier"};

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  (void)result;
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(
      stderrOutput.find("predecessor reset-frontier precheck"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor query budget exhausted"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesCheapResetConstantFactsWhenExactResetChecksAreDisabled) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(3)));
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1)
          .has_value());

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(stderrOutput.find("post_bootstrap_steps=1"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesResetSpecializedRelationsBeforeExactRootResetFrontier) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {x, y, w};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Not(BoolExpr::Var(y)), BoolExpr::Var(w))));
  // The reset transition creates y == w at the F[0] frontier, but neither bit
  // is a reset constant.  This guards the sampled ASIC path where exact deeper
  // reset checks are disabled, yet PDR should still learn the abstract F[0]
  // predecessor is outside the concrete post-reset image before doing a wide
  // root-cube validation query.
  problem.transitions0.emplace_back(y, BoolExpr::Var(w));
  problem.transitions0.emplace_back(w, BoolExpr::Var(w));
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("post_bootstrap_steps=1"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesResetSpecializedExpressionSatBeforeExactRootResetFrontier) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t a = 5;
  constexpr size_t b = 6;
  constexpr size_t reset = 7;
  problem.state0Symbols = {x, y, w, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, a, b, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  // y and w reset to equivalent XNOR forms that are not reduced by the cheap
  // structural implication rules. This keeps the test focused on the bounded
  // expression-SAT shortcut instead of the faster syntactic reset proofs.
  problem.transitions0.emplace_back(
      y,
      BoolExpr::Or(
          BoolExpr::And(BoolExpr::Var(a), BoolExpr::Var(b)),
          BoolExpr::And(BoolExpr::Not(BoolExpr::Var(a)), BoolExpr::Not(BoolExpr::Var(b)))));
  problem.transitions0.emplace_back(
      w,
      BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(a), BoolExpr::Var(b))));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("reset-specialized expression conflict"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("reset-specialized expression solver_profile=reset_expression"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineMinimizesResetSpecializedExpressionSatConflictToPair) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t e0 = 5;
  constexpr size_t e1 = 6;
  constexpr size_t a = 7;
  constexpr size_t b = 8;
  constexpr size_t reset = 9;
  problem.state0Symbols = {x, y, w, e0, e1, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, e0, e1, a, b, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  BoolExpr* badDriver = BoolExpr::And(
      BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w))),
      BoolExpr::And(BoolExpr::Var(e0), BoolExpr::Var(e1)));
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(reset)), badDriver));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::And(
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Var(b)),
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Not(BoolExpr::Var(b)))));
  problem.transitions0.emplace_back(w, BoolExpr::Var(a));
  problem.transitions0.emplace_back(e0, BoolExpr::Var(e0));
  problem.transitions0.emplace_back(e1, BoolExpr::Var(e1));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/32,
      /*preciseBadCubeStateLimit=*/32,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("via=pair_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("->2 via=pair_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineOrdersResetExpressionPairProbesBySupport) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t e0 = 5;
  constexpr size_t a = 6;
  constexpr size_t b = 7;
  constexpr size_t reset = 8;
  constexpr size_t firstWideLeaf = 9;
  constexpr size_t wideLeafCount = 12;
  problem.state0Symbols = {x, y, w, e0, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, e0, a, b, reset};
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(
              BoolExpr::And(BoolExpr::Var(y), BoolExpr::Var(w)),
              BoolExpr::Not(BoolExpr::Var(e0)))));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::And(
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Var(b)),
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Not(BoolExpr::Var(b)))));
  problem.transitions0.emplace_back(w, BoolExpr::Not(BoolExpr::Var(a)));
  BoolExpr* wideExpr = BoolExpr::Var(firstWideLeaf);
  for (size_t index = 1; index < wideLeafCount; ++index) {
    wideExpr =
        BoolExpr::Xor(wideExpr, BoolExpr::Var(firstWideLeaf + index));
  }
  problem.transitions0.emplace_back(e0, wideExpr);
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/32,
      /*preciseBadCubeStateLimit=*/32,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // The conflicting pair is same-valued but has tiny support; the opposite
  // valued pairs pull in the wide e0 cone.  Probe support first so sampled AES
  // failures do not spend time proving wide SAT distractor pairs.
  EXPECT_NE(
      stderrOutput.find(
          "reset-specialized expression solve cube=2 target_step=1 support=2 "),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("via=pair_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression solve cube=2 support=13"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression solve cube=2 support=14"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineContinuesResetExpressionPairProbesPastSatDistractors) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t d0 = 3;
  constexpr size_t d1 = 4;
  constexpr size_t d2 = 5;
  constexpr size_t d3 = 6;
  constexpr size_t y = 7;
  constexpr size_t w = 8;
  constexpr size_t p0 = 9;
  constexpr size_t p1 = 10;
  constexpr size_t p2 = 11;
  constexpr size_t p3 = 12;
  constexpr size_t a = 13;
  constexpr size_t b = 14;
  constexpr size_t reset = 15;
  problem.state0Symbols = {x, d0, d1, d2, d3, y, w, p0, p1, p2, p3, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {
      x, d0, d1, d2, d3, y, w, p0, p1, p2, p3, a, b, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(
              BoolExpr::And(BoolExpr::Var(d0), BoolExpr::Var(d1)),
              BoolExpr::And(
                  BoolExpr::And(BoolExpr::Var(d2), BoolExpr::Var(d3)),
                  BoolExpr::And(BoolExpr::Var(y),
                                BoolExpr::Not(BoolExpr::Var(w)))))));
  problem.transitions0.emplace_back(d0, BoolExpr::Var(p0));
  problem.transitions0.emplace_back(d1, BoolExpr::Var(p1));
  problem.transitions0.emplace_back(d2, BoolExpr::Var(p2));
  problem.transitions0.emplace_back(d3, BoolExpr::Var(p3));
  problem.transitions0.emplace_back(y, BoolExpr::Xor(BoolExpr::Var(a), BoolExpr::Var(b)));
  problem.transitions0.emplace_back(
      w,
      BoolExpr::Or(
          BoolExpr::And(BoolExpr::Var(a), BoolExpr::Not(BoolExpr::Var(b))),
          BoolExpr::And(BoolExpr::Not(BoolExpr::Var(a)), BoolExpr::Var(b))));
  problem.transitions0.emplace_back(p0, BoolExpr::Var(p0));
  problem.transitions0.emplace_back(p1, BoolExpr::Var(p1));
  problem.transitions0.emplace_back(p2, BoolExpr::Var(p2));
  problem.transitions0.emplace_back(p3, BoolExpr::Var(p3));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/32,
      /*preciseBadCubeStateLimit=*/32,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // Boolean-equivalent but structurally different reset expressions should
  // still shrink a wide reset cube through the SAT pair-probe path before any
  // exact reset-frontier query is needed.
  EXPECT_NE(
      stderrOutput.find("reset-specialized expression conflict cube=7->2 via=pair_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFindsResetExpressionTripleConflictWhenPairsAreSat) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t z = 5;
  constexpr size_t d = 6;
  constexpr size_t a = 7;
  constexpr size_t b = 8;
  constexpr size_t reset = 9;
  problem.state0Symbols = {x, y, w, z, d, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, z, d, a, b, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(
              BoolExpr::And(BoolExpr::Var(y), BoolExpr::Var(w)),
              BoolExpr::And(BoolExpr::Var(z), BoolExpr::Var(d)))));
  problem.transitions0.emplace_back(y, BoolExpr::Var(a));
  problem.transitions0.emplace_back(w, BoolExpr::Var(b));
  problem.transitions0.emplace_back(z, BoolExpr::Xor(BoolExpr::Var(a), BoolExpr::Var(b)));
  problem.transitions0.emplace_back(d, BoolExpr::Var(d));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/32,
      /*preciseBadCubeStateLimit=*/32,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  // a=1, b=1, and a^b=1 is impossible, but every pair of those literals is
  // satisfiable.  The triple probe should learn that smaller reset-image
  // conflict before the optional full-cube SAT fallback.
  EXPECT_NE(
      stderrOutput.find(
          "reset-specialized expression conflict cube=5->3 via=triple_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

KInductionProblem makeWideResetExpressionSatShortcutProblem(
    size_t wideLeafCount) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t a = 5;
  constexpr size_t reset = 6;
  constexpr size_t firstWideLeaf = 7;
  problem.state0Symbols = {x, y, w, a};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, a, reset};
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Var(w))));
  BoolExpr* wideExpr = BoolExpr::Var(firstWideLeaf);
  for (size_t index = 1; index < wideLeafCount; ++index) {
    wideExpr =
        BoolExpr::Xor(wideExpr, BoolExpr::Var(firstWideLeaf + index));
  }
  problem.transitions0.emplace_back(y, wideExpr);
  problem.transitions0.emplace_back(w, BoolExpr::Var(a));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

KInductionProblem makeWideMultiLiteralResetExpressionSatShortcutProblem(
    size_t wideLeafCount) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t extra0 = 5;
  constexpr size_t extra1 = 6;
  constexpr size_t a = 7;
  constexpr size_t reset = 8;
  constexpr size_t firstWideLeaf = 9;
  problem.state0Symbols = {x, y, w, extra0, extra1, a};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, extra0, extra1, a, reset};
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  BoolExpr* rootCubeDriver = BoolExpr::And(
      BoolExpr::Var(y),
      BoolExpr::And(BoolExpr::Var(w), BoolExpr::Var(extra0)));
  rootCubeDriver = BoolExpr::And(rootCubeDriver, BoolExpr::Var(extra1));
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(reset)), rootCubeDriver));
  BoolExpr* wideExpr = BoolExpr::Var(firstWideLeaf);
  for (size_t index = 1; index < wideLeafCount; ++index) {
    wideExpr =
        BoolExpr::Xor(wideExpr, BoolExpr::Var(firstWideLeaf + index));
  }
  problem.transitions0.emplace_back(y, wideExpr);
  problem.transitions0.emplace_back(w, BoolExpr::Var(a));
  problem.transitions0.emplace_back(extra0, BoolExpr::Var(extra0));
  problem.transitions0.emplace_back(extra1, BoolExpr::Var(extra1));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  for (size_t index = 0; index < wideLeafCount; ++index) {
    const size_t symbol = firstWideLeaf + index;
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }
  problem.bad = BoolExpr::Var(x);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  return problem;
}

void addLargeDualRailResetFrontierSurfaceForTest(
    KInductionProblem& problem,
    size_t railPairCount = 10001) {
  problem.usesDualRailStateEncoding = true;
  for (size_t index = 0; index < railPairCount; ++index) {
    problem.dualRailStatePairs.push_back({1000 + index * 2, 1001 + index * 2});
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineAttemptsModerateWideResetExpressionSatShortcut) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/128);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // AES sampling kept useful support-129/135 reset-image pair proofs. This
  // guarded case keeps only local proof shapes eligible for the
  // reset-expression SAT path; whether the particular cube is SAT or UNSAT
  // remains a solver result.
  EXPECT_NE(
      stderrOutput.find("reset-specialized expression solve"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("support=129"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("miss reason=full_sat_support_cap"),
            std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsResetExpressionShortcutWhenResourceLimitHits) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/128);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar proofConflictLimit(
      "KEPLER_SEC_PDR_RESET_EXPRESSION_CONFLICT_LIMIT", "0");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Reset-expression SAT is an optional shortcut. If Kissat hits the local
  // resource cap, PDR must report a miss and continue through the normal
  // validation/refinement path instead of treating UNKNOWN as UNSAT.
  EXPECT_NE(
      stderrOutput.find("miss reason=solver_resource_limit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression conflict"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsBroadMultiLiteralResetExpressionSatShortcut) {
  KInductionProblem problem =
      makeWideMultiLiteralResetExpressionSatShortcutProblem(
          /*wideLeafCount=*/600);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Sampling showed that admitting broad four-literal reset-image cubes into
  // the full SAT shortcut simply moved the wall from exact reset-frontier BMC
  // into Kissat.  Pair/triple probes are still allowed above, but the complete
  // multi-literal SAT fallback must stay below the smaller support cap.
  EXPECT_NE(
      stderrOutput.find("miss reason=full_sat_support_cap"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression solve cube=5 support=603"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsHighSupportSmallResetExpressionSatShortcut) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/900);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The useful sampled AES shortcuts had support 129/135. Wider small cubes
  // still fall through before opening Kissat, so this optional proof path
  // cannot become a whole-chip SAT query.
  EXPECT_NE(
      stderrOutput.find("miss reason=full_sat_support_cap"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression solve cube=3 support=901"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsOneShotResetValidationForWideFinalRootCegar) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/900);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Final projected-counterexample repair can encounter a small root cube whose
  // reset image has a huge support.  Keep that shape behind the reset-expression
  // support cap instead of opening either the one-shot or cached exact reset
  // validator; BlackParrot/AES sampling showed those broad exact queries are the
  // runtime wall.
  EXPECT_NE(
      stderrOutput.find("miss reason=full_sat_support_cap"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("mode=cached_assumptions"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("mode=one_shot_unit_clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized concrete-frame conflict"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineUsesCachedResetValidationForWideDualRailRootCegar) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/900);
  problem.bad = BoolExpr::And(
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
      BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  addLargeDualRailResetFrontierSurfaceForTest(problem);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find(
          "concrete cube reachability begin cube=3 max_step=1 "
          "mode=cached_assumptions"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          "concrete cube reachability begin cube=3 max_step=1 "
          "mode=one_shot_unit_clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail PDR transient caches reason=pdr_run_exit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("bad_solver=1"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("predecessor_solver=1"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail reset-frontier memory reason=pdr_run_exit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("skipped large dual-rail reset-frontier precheck"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail reset-frontier memory "
          "reason=before_singleton_reset_frontier_core"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "reset frontier one-shot cube coi post_bootstrap_steps=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail reset-frontier memory "
          "reason=reachable_singleton_reset_frontier_probe"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail reset-frontier memory "
          "reason=cheap_concrete_frame_conflict"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineLargeDualRailRunExitReleasesTransientCaches) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  constexpr size_t reset = 3;
  problem.state0Symbols = {state};
  problem.inputSymbols = {reset};
  problem.allSymbols = {state, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{state, false}};
  problem.transitions0.emplace_back(state, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(state);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  addLargeDualRailResetFrontierSurfaceForTest(problem);
  auto lazyTransitions = std::make_shared<LazyTransitionStore>();
  lazyTransitions->remappedByStateSymbol.emplace(state, BoolExpr::Var(state));
  lazyTransitions->remapMemoByDesign[0].emplace(
      BoolExpr::Var(state), BoolExpr::Var(state));
  lazyTransitions->dualRailRemapMemoByDesign[0].emplace(
      BoolExpr::Var(state),
      DualRailBoolExpr{BoolExpr::Var(state), BoolExpr::Not(BoolExpr::Var(state))});
  lazyTransitions->supportByStateSymbol.emplace(
      state, std::set<size_t>{state, reset});
  lazyTransitions->nodeCountByStateSymbol.emplace(state, 2);
  problem.lazyTransitions = lazyTransitions;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result =
      engine.run(/*maxFrames=*/0, /*resetBootstrapFrameCheckedSafe=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail PDR transient caches reason=pdr_run_exit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "released large dual-rail reset-frontier memory reason=pdr_run_exit"),
      std::string::npos)
      << stderrOutput;
  // The release hook must drop materialized BoolExpr remaps, but keep compact
  // support metadata so sibling dual-rail PDR batches do not repeat the same
  // lazy transition DAG walks.
  EXPECT_TRUE(lazyTransitions->remappedByStateSymbol.empty());
  EXPECT_TRUE(lazyTransitions->remapMemoByDesign[0].empty());
  EXPECT_TRUE(lazyTransitions->dualRailRemapMemoByDesign[0].empty());
  EXPECT_FALSE(lazyTransitions->supportByStateSymbol.empty());
  EXPECT_FALSE(lazyTransitions->nodeCountByStateSymbol.empty());
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBudgetsWideDualRailConcreteRootValidation) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/900);
  constexpr size_t firstWideLeaf = 7;
  constexpr size_t wideBadLeafCount = 29;
  BoolExpr* bad = BoolExpr::And(
      BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
      BoolExpr::Var(4));
  for (size_t index = 0; index < wideBadLeafCount; ++index) {
    bad = BoolExpr::And(bad, BoolExpr::Var(firstWideLeaf + index));
  }
  problem.bad = bad;
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.usesDualRailStateEncoding = true;
  for (size_t index = 0; index < 4097; ++index) {
    problem.dualRailStatePairs.push_back({1000 + index * 2, 1001 + index * 2});
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar frontierStateLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_RESET_FRONTIER_STATE_SYMBOL_LIMIT", "8193");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true,
      /*maxProjectedCounterexampleRefinements=*/2);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("skipped large dual-rail concrete root validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("concrete cube reachability begin cube=32"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBudgetsDeepLargeDualRailExactResetFrontierQuery) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.dualRailStatePairs = {{20, 21}};
  problem.resetBootstrapCycles = 3;
  problem.state0Symbols = {2, 3, 4, 5, 6, 7};
  problem.allSymbols = {2, 3, 4, 5, 6, 7};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                        BoolExpr::Not(BoolExpr::Var(5))),
          BoolExpr::And(BoolExpr::Not(BoolExpr::Var(6)),
                        BoolExpr::Not(BoolExpr::Var(7)))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false},
      {5, false}, {6, false}, {7, false}};
  problem.initializedStateCount = 6;
  problem.totalStateCount = 6;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.transitions0.emplace_back(6, BoolExpr::Var(5));
  problem.transitions0.emplace_back(7, BoolExpr::Var(6));
  problem.bad = BoolExpr::And(BoolExpr::Var(7), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 5)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar frontierStateLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_RESET_FRONTIER_STATE_SYMBOL_LIMIT", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/4,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(5);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          "skipped fresh large dual-rail exact reset-frontier query "
          "reason=concrete_frame_reachability post_bootstrap_steps=4"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "reset-specialized expression relaxed_budget cube="),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset frontier cube coi post_bootstrap_steps=4"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsVeryWideResetExpressionSatShortcut) {
  KInductionProblem problem =
      makeWideResetExpressionSatShortcutProblem(/*wideLeafCount=*/1032);

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  (void)engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Reset-expression SAT is only a shortcut.  If the local reset support is
  // broader than the bounded ASIC-sized proof path and no cheap relation proof
  // applies, skip the optional SAT query and let the caller's exact
  // validation/refinement path decide.
  EXPECT_NE(
      stderrOutput.find("miss reason=full_sat_support_cap"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-specialized expression solve cube=2 support=1033"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCanonicalizesResetSpecializedExpressionsBeforeSat) {
  KInductionProblem problem;
  constexpr size_t x0 = 2;
  constexpr size_t x1 = 3;
  constexpr size_t y = 4;
  constexpr size_t w = 5;
  constexpr size_t e0 = 6;
  constexpr size_t e1 = 7;
  constexpr size_t a = 8;
  constexpr size_t b = 9;
  constexpr size_t reset = 10;

  problem.state0Symbols = {x0, x1, y, w, e0, e1, a, b};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x0, x1, y, w, e0, e1, a, b, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  const auto gatedBad =
      [&](BoolExpr* extra) {
        return BoolExpr::And(
            BoolExpr::Not(BoolExpr::Var(reset)),
            BoolExpr::And(
                BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w))),
                extra));
      };
  problem.transitions0.emplace_back(x0, gatedBad(BoolExpr::Var(e0)));
  problem.transitions0.emplace_back(x1, gatedBad(BoolExpr::Var(e1)));
  // Both bad predecessors are reset-unreachable for the same reason:
  // y' = a | (a & b) is Boolean-equivalent to w' = a. The sampled AES run was
  // spending time in the reset-specialized SAT fallback for this shape, so the
  // canonical pass should learn the conflict before invoking that solver.
  problem.transitions0.emplace_back(
      y,
      BoolExpr::Or(BoolExpr::Var(a),
                   BoolExpr::And(BoolExpr::Var(a), BoolExpr::Var(b))));
  problem.transitions0.emplace_back(w, BoolExpr::Var(a));
  problem.transitions0.emplace_back(e0, BoolExpr::Var(e0));
  problem.transitions0.emplace_back(e1, BoolExpr::Var(e1));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Or(BoolExpr::Var(x0), BoolExpr::Var(x1));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/32,
      /*preciseBadCubeStateLimit=*/32,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("via=canonical"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineValidatedBadFormulaLearningRepairsBeforePostBootstrapPrecheck) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.inputSymbols = {4};
  problem.allSymbols = {2, 3, 4};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{4, false}};
  problem.bootstrapStateAssignments = {{2, false}};
  problem.transitions0.emplace_back(2, BoolExpr::And(BoolExpr::Var(4), BoolExpr::Var(3)));
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(2);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 1)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/true,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("k-induction base coi"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineExactResetFrontierBlocksBeforeRootMinimization) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t z = 5;
  constexpr size_t reset = 6;
  problem.state0Symbols = {x, y, w, z};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, z, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  // The concrete reset step creates y == w, but that relation is intentionally
  // not summarized in F0. Exact reset-frontier predecessor checks should block
  // the abstract predecessor before PDR learns a root obligation and starts any
  // optional root-cube minimization work.
  problem.transitions0.emplace_back(y, BoolExpr::Var(w));
  problem.transitions0.emplace_back(w, BoolExpr::Var(w));
  problem.transitions0.emplace_back(z, BoolExpr::Var(z));
  problem.bad = BoolExpr::And(BoolExpr::Var(x), BoolExpr::Var(z));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(stderrOutput.find("post_bootstrap_steps=1"), std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "exact_reset_frontier=1 mode=cached_assumptions result=unsat"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset-frontier core"), std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("post_bootstrap_steps=0 frames=2 "
                        "solver_symbols=5 transition_targets=4 cube_literals=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsExactResetPrecheckForUnprojectedPredecessorQuery) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t z = 5;
  constexpr size_t reset = 6;
  problem.state0Symbols = {x, y, w, z};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, z, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  problem.transitions0.emplace_back(y, BoolExpr::Var(w));
  problem.transitions0.emplace_back(w, BoolExpr::Var(w));
  problem.transitions0.emplace_back(z, BoolExpr::Var(z));
  problem.bad = BoolExpr::And(BoolExpr::Var(x), BoolExpr::Var(z));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/0,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("exact_reset_frontier=skipped"),
      std::string::npos)
      << stderrOutput;
  // In unprojected mode the normal predecessor SAT query is already exact.
  // Do not spend the sampled AES wall time on an extra one-step reset-image
  // query before that exact predecessor query has a chance to run.
  EXPECT_EQ(
      stderrOutput.find("reset frontier cube coi post_bootstrap_steps=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineValidatedLearningKeepsRootResetFrontierRefinementDisabled) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t w = 4;
  constexpr size_t z = 5;
  constexpr size_t reset = 6;
  problem.state0Symbols = {x, y, w, z};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, w, z, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  problem.transitions0.emplace_back(y, BoolExpr::Var(w));
  problem.transitions0.emplace_back(w, BoolExpr::Var(w));
  problem.transitions0.emplace_back(z, BoolExpr::Var(z));
  problem.bad = BoolExpr::And(BoolExpr::Var(x), BoolExpr::Var(z));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/true,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset frontier cube coi post_bootstrap_steps=0"),
            std::string::npos)
      << stderrOutput;
  EXPECT_EQ(stderrOutput.find("reset-frontier core"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineMinimizesExactResetPredecessorCoreAfterRelaxedPrecheck) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t q = 3;
  constexpr size_t y = 4;
  constexpr size_t w = 5;
  constexpr size_t z = 6;
  constexpr size_t a = 7;
  constexpr size_t reset = 8;
  problem.state0Symbols = {x, q, y, w, z, a};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, q, y, w, z, a, reset};
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  problem.transitions0.emplace_back(
      q,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::And(BoolExpr::Var(y), BoolExpr::Not(BoolExpr::Var(w)))));
  // The reset/bootstrap frontier has y == w. Therefore target x == 1 is
  // impossible in one post-reset step, and q == 1 is an independent sibling
  // core worth seeding from the same exact reset-frontier context.
  problem.transitions0.emplace_back(y, BoolExpr::Var(w));
  problem.transitions0.emplace_back(w, BoolExpr::Var(w));
  problem.transitions0.emplace_back(z, BoolExpr::Var(z));
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.bad = BoolExpr::And(
      BoolExpr::And(BoolExpr::Var(x), BoolExpr::Var(q)),
      BoolExpr::And(BoolExpr::Var(z), BoolExpr::Var(a)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;
  problem.observedOutputExprs0 = {problem.bad};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("bad cube level=1 source=precise support=4 cube=4"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("limit=6"), std::string::npos) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("exact reset-predecessor core cube=4->1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("seeded exact reset-predecessor sibling cores cube=4 seeded=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "seeded exact reset-predecessor sibling cores cube=4 seeded=1 "
          "post_bootstrap_steps=1 cached=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "learned exact reset-predecessor singleton clauses level=1 added=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("generalized blocked cube level=1 size=4->2"),
      std::string::npos)
      << stderrOutput;

  testing::internal::CaptureStderr();
  PDREngine cachedEngine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      PDREngine::kDefaultBoundedRootGeneralizationAttempts,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/true);
  const auto cachedResult = cachedEngine.run(3);
  const std::string cachedStderrOutput =
      testing::internal::GetCapturedStderr();

  EXPECT_EQ(cachedResult.status, PDRStatus::Equivalent);
  EXPECT_NE(
      cachedStderrOutput.find("exact reset-predecessor core cube=4->1"),
      std::string::npos)
      << cachedStderrOutput;
  EXPECT_NE(
      cachedStderrOutput.find(
          "learned exact reset-predecessor singleton clauses level=1 added=1"),
      std::string::npos)
      << cachedStderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityExtractsSmallUnreachableCubeCore) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t freeState0 = 3;
  constexpr size_t freeState1 = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {resetForcedLow, freeState0, freeState1};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, freeState0, freeState1, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(freeState0, BoolExpr::Var(freeState0));
  problem.transitions0.emplace_back(freeState1, BoolExpr::Var(freeState1));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> wideUnreachableCube = {
      {resetForcedLow, true}, {freeState0, true}, {freeState1, false}};

  ASSERT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      wideUnreachableCube,
      0));
  ASSERT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      std::vector<std::pair<size_t, bool>>{{resetForcedLow, true}},
      0));
  const auto core = findResetFrontierUnreachableCubeCore(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      wideUnreachableCube,
      0);

  ASSERT_TRUE(core.has_value());
  EXPECT_LT(core->size(), wideUnreachableCube.size());
  EXPECT_FALSE(isStateCubeReachableAtResetFrontier(
      *context, KEPLER_FORMAL::Config::SolverType::KISSAT, *core, 0));
  EXPECT_NE(
      std::find(
          core->begin(), core->end(), std::pair<size_t, bool>{resetForcedLow, true}),
      core->end());
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityOneShotMatchesCachedAssumptionSolver) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t freeState = 3;
  constexpr size_t reset = 4;
  problem.state0Symbols = {resetForcedLow, freeState};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, freeState, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(freeState, BoolExpr::Var(freeState));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> unreachableCube = {
      {resetForcedLow, true}, {freeState, true}};
  const std::vector<std::pair<size_t, bool>> reachableCube = {
      {resetForcedLow, false}, {freeState, true}};

  // The one-shot path is the same exact bounded-prefix query as the cached
  // assumption solver, but it gives final PDR candidate validation a way to use
  // the selected SEC solver instead of a long-lived incremental Glucose query.
  EXPECT_EQ(
      isStateCubeReachableAtResetFrontier(
          *context,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          unreachableCube,
          0),
      isStateCubeReachableAtResetFrontierOneShot(
          *context,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          unreachableCube,
          0));
  EXPECT_EQ(
      isStateCubeReachableAtResetFrontier(
          *context,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          reachableCube,
          0),
      isStateCubeReachableAtResetFrontierOneShot(
          *context,
          KEPLER_FORMAL::Config::SolverType::KISSAT,
          reachableCube,
          0));
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilitySkipsBroadRelaxedCachedPrecheck) {
  KInductionProblem problem;
  constexpr size_t observed = 2;
  constexpr size_t reset = 3;
  constexpr size_t supportBase = 100;
  constexpr size_t supportCount = 300;
  problem.state0Symbols = {observed};
  problem.inputSymbols = {reset};
  problem.allSymbols = {observed, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};

  BoolExpr* observedNext = BoolExpr::createFalse();
  for (size_t offset = 0; offset < supportCount; ++offset) {
    const size_t symbol = supportBase + offset;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    observedNext = BoolExpr::Or(observedNext, BoolExpr::Var(symbol));
  }
  problem.transitions0.emplace_back(observed, observedNext);

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_TRUE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      std::vector<std::pair<size_t, bool>>{
          {observed, true},
          {supportBase, false},
          {supportBase + 1, false},
          {supportBase + 2, false},
          {supportBase + 3, false}},
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The relaxed precheck is only a local UNSAT shortcut. If it still pulls a
  // broad transition surface, skip solving it and fall through to the exact
  // cached reset-frontier query instead of creating an unbounded PDR wall.
  EXPECT_NE(
      stderrOutput.find(
          "reset frontier relaxed cached precheck skipped reason=coi_cap"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("reset frontier cube coi"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityAllowsTinyBroadRelaxedCachedPrecheck) {
  KInductionProblem problem;
  constexpr size_t observed = 2;
  constexpr size_t reset = 3;
  constexpr size_t supportBase = 100;
  constexpr size_t supportCount = 300;
  problem.state0Symbols = {observed};
  problem.inputSymbols = {reset};
  problem.allSymbols = {observed, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};

  BoolExpr* observedNext = BoolExpr::createFalse();
  for (size_t offset = 0; offset < supportCount; ++offset) {
    const size_t symbol = supportBase + offset;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    observedNext = BoolExpr::Or(observedNext, BoolExpr::Var(symbol));
  }
  problem.transitions0.emplace_back(observed, observedNext);

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_TRUE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      std::vector<std::pair<size_t, bool>>{{observed, true}},
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // MockAlu-like PDR leaves use tiny bad cubes whose relaxed reset-frontier COI
  // is larger than the broad-cube cap.  They should still get the bounded
  // relaxed exact precheck before falling back to the heavier cached solver.
  EXPECT_NE(
      stderrOutput.find("reset frontier relaxed cached cube coi"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          "reset frontier relaxed cached precheck skipped reason=coi_cap"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityUsesValidatedFrameInvariantAfterStartup) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t z = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {x, y, z};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, z, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      x,
      BoolExpr::Not(makeEqualityExpr(BoolExpr::Var(y), BoolExpr::Var(z))));
  problem.transitions0.emplace_back(y, BoolExpr::Var(y));
  problem.transitions0.emplace_back(z, BoolExpr::Var(z));

  const TransitionExprResolver transitionByState(problem);
  const std::vector<std::pair<size_t, bool>> targetCube = {{x, true}};
  const auto plainContext =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  ASSERT_TRUE(isStateCubeReachableAtResetFrontierOneShot(
      *plainContext,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      targetCube,
      1));

  // PDR validates the invariant separately before passing it into this helper.
  // The bounded transition prefix is unchanged, but from the startup frontier
  // onward y==z makes x unreachable one post-bootstrap step later.
  BoolExpr* frameInvariant =
      makeEqualityExpr(BoolExpr::Var(y), BoolExpr::Var(z));
  const auto invariantContext =
      makeResetFrontierReachabilityContext(
          problem, transitionByState, frameInvariant);

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(isStateCubeReachableAtResetFrontierOneShot(
      *invariantContext,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      targetCube,
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("frame_invariant_symbols=2"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityReusesCachedUnreachableCores) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t neighborState0 = 3;
  constexpr size_t neighborState1 = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {resetForcedLow, neighborState0, neighborState1};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, neighborState0, neighborState1, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(neighborState0, BoolExpr::Var(neighborState0));
  problem.transitions0.emplace_back(neighborState1, BoolExpr::Var(neighborState1));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> firstUnreachableCube = {
      {resetForcedLow, true}, {neighborState0, true}, {neighborState1, false}};
  const std::vector<std::pair<size_t, bool>> neighboringUnreachableCube = {
      {resetForcedLow, true}, {neighborState0, false}};

  ASSERT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      firstUnreachableCube,
      0));

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      neighboringUnreachableCube,
      0));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("reset frontier cached unreachable core hit"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityCachesPostBootstrapFailedAssumptionCores) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t neighborState0 = 3;
  constexpr size_t neighborState1 = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {resetForcedLow, neighborState0, neighborState1};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, neighborState0, neighborState1, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(neighborState0, BoolExpr::Var(neighborState0));
  problem.transitions0.emplace_back(neighborState1, BoolExpr::Var(neighborState1));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> firstUnreachableCube = {
      {resetForcedLow, true}, {neighborState0, true}, {neighborState1, false}};
  const std::vector<std::pair<size_t, bool>> neighboringUnreachableCube = {
      {resetForcedLow, true}, {neighborState0, false}};

  ASSERT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      firstUnreachableCube,
      1));

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      neighboringUnreachableCube,
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("reset frontier relaxed cached cube coi"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset frontier cube coi"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityUsesPriorCoreAsSafePrefixBlocker) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t reset = 3;
  problem.state0Symbols = {resetForcedLow};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> unreachableCube = {
      {resetForcedLow, true}};

  ASSERT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      unreachableCube,
      0));

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      unreachableCube,
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("reset frontier previous unreachable blockers=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       ResetFrontierReachabilityCachesPostBootstrapOneShotFailures) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t neighborState = 3;
  constexpr size_t reset = 4;
  problem.state0Symbols = {resetForcedLow, neighborState};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, neighborState, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(neighborState, BoolExpr::Var(neighborState));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> unreachableCube = {
      {resetForcedLow, true}, {neighborState, false}};

  // One-shot PDR prechecks at post-bootstrap depths must populate the shared
  // unreachable-core cache too; otherwise a repeated target rebuilds the same
  // reset COI instead of taking the cheap cache hit.
  ASSERT_FALSE(isStateCubeReachableAtResetFrontierOneShot(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      unreachableCube,
      1));

  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_FALSE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      unreachableCube,
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("reset frontier cached unreachable core hit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset frontier cube coi"),
      std::string::npos)
      << stderrOutput;
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
  // BlackParrot PDR those queries happen while rebuilding reset-frontier COIs
  // across many output batches; sharing this metadata avoids repeatedly
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
       ResetFrontierReachabilityReusesWiderCachedSolverForSubsetCube) {
  KInductionProblem problem;
  constexpr size_t resetForcedLow = 2;
  constexpr size_t neighborState0 = 3;
  constexpr size_t neighborState1 = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {resetForcedLow, neighborState0, neighborState1};
  problem.inputSymbols = {reset};
  problem.allSymbols = {resetForcedLow, neighborState0, neighborState1, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments = {{resetForcedLow, false}};
  problem.transitions0.emplace_back(
      resetForcedLow,
      BoolExpr::And(
          BoolExpr::Not(BoolExpr::Var(reset)),
          BoolExpr::Var(resetForcedLow)));
  problem.transitions0.emplace_back(neighborState0, BoolExpr::Var(neighborState0));
  problem.transitions0.emplace_back(neighborState1, BoolExpr::Var(neighborState1));

  const TransitionExprResolver transitionByState(problem);
  const auto context =
      makeResetFrontierReachabilityContext(problem, transitionByState);
  const std::vector<std::pair<size_t, bool>> wideReachableCube = {
      {resetForcedLow, false}, {neighborState0, true}, {neighborState1, false}};
  const std::vector<std::pair<size_t, bool>> subsetReachableCube = {
      {resetForcedLow, false}, {neighborState0, false}};

  ASSERT_TRUE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      wideReachableCube,
      1));

  // PDR often checks neighboring cubes where a previous reset-frontier solver
  // already covers a wider COI. Reusing that exact solver avoids rebuilding
  // transition support and clauses for every small cube variant.
  const ScopedEnvVar kiDiag("KEPLER_SEC_KI_DIAG", "1");
  testing::internal::CaptureStderr();
  EXPECT_TRUE(isStateCubeReachableAtResetFrontier(
      *context,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      subsetReachableCube,
      1));
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("reset frontier solver superset cache hit"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineSkipsExactResetPrecheckAboveConfiguredSupportLimit) {
  KInductionProblem problem;
  constexpr size_t y = 2;
  constexpr size_t reset = 100;
  constexpr size_t supportBase = 200;
  constexpr size_t supportCount = 300;
  problem.state0Symbols.push_back(y);
  problem.allSymbols.push_back(y);
  BoolExpr* nextY = BoolExpr::createFalse();
  problem.bootstrapStateAssignments.push_back({y, false});
  for (size_t index = 0; index < supportCount; ++index) {
    const size_t symbol = supportBase + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.bootstrapStateAssignments.push_back({symbol, false});
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    nextY = BoolExpr::Or(nextY, BoolExpr::Var(symbol));
  }
  problem.inputSymbols = {reset};
  problem.allSymbols.push_back(reset);
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(y, nextY);
  problem.bad = BoolExpr::Var(y);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar exactResetPrecheckLimit(
      "KEPLER_SEC_PDR_EXACT_RESET_PRECHECK_SUPPORT_LIMIT", "256");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("exact_reset_frontier=skipped"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineProjectsLevelZeroResetPredecessorsAfterConcretePrecheck) {
  KInductionProblem problem;
  constexpr size_t y = 2;
  constexpr size_t reset = 3;
  constexpr size_t supportBase = 100;
  constexpr size_t supportCount = 64;
  problem.state0Symbols.push_back(y);
  problem.inputSymbols = {reset};
  problem.allSymbols = {y, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.bootstrapStateAssignments.push_back({y, false});

  BoolExpr* nextY = BoolExpr::createFalse();
  for (size_t index = 0; index < supportCount; ++index) {
    const size_t symbol = supportBase + index;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    // Make the concrete reset predecessor real, but leave a wide support cone.
    // PDR should not carry the full 64-bit support cube after the exact
    // reset-frontier precheck already established concrete reachability.
    problem.bootstrapStateAssignments.push_back({symbol, index == 0});
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
    nextY = BoolExpr::Or(nextY, BoolExpr::Var(symbol));
  }
  problem.transitions0.emplace_back(y, nextY);
  problem.bad = BoolExpr::Var(y);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/8,
      /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_NE(
      stderrOutput.find(
          "exact_reset_frontier=1 mode=cached_assumptions result=sat"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find("predecessor_cube=1"), std::string::npos)
      << stderrOutput;
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
       PDREngineDoesNotLoopOnProjectedResetFrontierRefinements) {
  KInductionProblem problem;
  constexpr size_t a = 2;
  constexpr size_t b = 3;
  constexpr size_t y = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {a, b, y};
  problem.inputSymbols = {reset};
  problem.allSymbols = {a, b, y, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};

  BoolExpr* resetDeasserted = BoolExpr::Not(BoolExpr::Var(reset));
  problem.transitions0.emplace_back(
      a, BoolExpr::And(resetDeasserted, BoolExpr::Var(a)));
  problem.transitions0.emplace_back(
      b, BoolExpr::And(resetDeasserted, BoolExpr::Var(b)));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::And(
          resetDeasserted,
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Var(b))));
  problem.bad = BoolExpr::Var(y);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  // Force the same condition sampled on BlackParrot: projected F[0] encoding
  // can omit one reset-frontier refinement even though the full frame already
  // blocks the predecessor cube. PDR must not keep re-enqueuing that stale
  // projected predecessor.
  const ScopedEnvVar clauseLimit(
      "KEPLER_SEC_PDR_PROJECTED_FRAME_CLAUSE_LIMIT", "1");
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/100);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineProjectedFrameClosureIsStableUnderClauseCap) {
  auto makeProblem = [](bool reverseSeeds) {
    KInductionProblem problem;
    constexpr size_t a = 2;
    constexpr size_t b = 3;
    constexpr size_t y = 4;
    constexpr size_t reset = 5;
    problem.state0Symbols = reverseSeeds ? std::vector<size_t>{y, b, a}
                                         : std::vector<size_t>{a, b, y};
    problem.inputSymbols = {reset};
    problem.allSymbols = reverseSeeds ? std::vector<size_t>{reset, y, b, a}
                                      : std::vector<size_t>{a, b, y, reset};
    problem.resetBootstrapCycles = 1;
    problem.resetBootstrapInputs = {{reset, true}};

    BoolExpr* resetDeasserted = BoolExpr::Not(BoolExpr::Var(reset));
    const auto addA = [&]() {
      problem.transitions0.emplace_back(
          a, BoolExpr::And(resetDeasserted, BoolExpr::Var(a)));
    };
    const auto addB = [&]() {
      problem.transitions0.emplace_back(
          b, BoolExpr::And(resetDeasserted, BoolExpr::Var(b)));
    };
    const auto addY = [&]() {
      problem.transitions0.emplace_back(
          y,
          BoolExpr::And(
              resetDeasserted,
              BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Var(b))));
    };
    if (reverseSeeds) {
      addY();
      addB();
      addA();
    } else {
      addA();
      addB();
      addY();
    }
    problem.bad = BoolExpr::Var(y);
    problem.property = BoolExpr::Not(problem.bad);
    problem.inductionProperty = problem.property;
    problem.inductionBad = problem.bad;
    return problem;
  };

  const ScopedEnvVar clauseLimit(
      "KEPLER_SEC_PDR_PROJECTED_FRAME_CLAUSE_LIMIT", "1");
  constexpr size_t kDeterministicQueryBudget = 19;

  for (const bool reverseSeeds : {false, true}) {
    KInductionProblem problem = makeProblem(reverseSeeds);
    ASSERT_FALSE(
        findBaseCounterexample(
            problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
            .has_value());
    PDREngine engine(
        problem,
        KEPLER_FORMAL::Config::SolverType::KISSAT,
        /*predecessorProjectionLimit=*/PDREngine::kDefaultPredecessorProjectionLimit,
        /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit,
        /*useExactFrameClauses=*/false,
        /*maxPredecessorQueries=*/kDeterministicQueryBudget);
    const auto result = engine.run(3);

    EXPECT_EQ(result.status, PDRStatus::Equivalent)
        << "reverseSeeds=" << reverseSeeds;
  }
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineExactRetriesWhenProjectedPredecessorIsAlreadyBlocked) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3};
  problem.allSymbols = {2, 3};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(2)),
      BoolExpr::Not(BoolExpr::Var(3)));
  problem.initialStateAssignments = {{2, false}, {3, false}};
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0.emplace_back(2, BoolExpr::Var(2));
  problem.transitions0.emplace_back(3, BoolExpr::Var(2));
  problem.bad = BoolExpr::Var(3);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // With a one-literal predecessor projection, the level-2 bad obligation for
  // b=1 first projects to a=1.  PDR then learns !a in F1 while blocking that
  // predecessor.  Re-querying b=1 against a projected frame can rediscover the
  // now-blocked a=1 cube forever unless the engine retries the parent query
  // against the exact learned frame before re-enqueueing it.
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      1);
  const auto result = engine.run(4);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCapsProjectedFrameRefinementsBeforeExactRetry) {
  KInductionProblem problem;
  constexpr size_t a = 2;
  constexpr size_t b = 3;
  constexpr size_t y = 4;
  constexpr size_t reset = 5;
  problem.state0Symbols = {a, b, y};
  problem.inputSymbols = {reset};
  problem.allSymbols = {a, b, y, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};

  BoolExpr* resetDeasserted = BoolExpr::Not(BoolExpr::Var(reset));
  problem.transitions0.emplace_back(
      a, BoolExpr::And(resetDeasserted, BoolExpr::Var(a)));
  problem.transitions0.emplace_back(
      b, BoolExpr::And(resetDeasserted, BoolExpr::Var(b)));
  problem.transitions0.emplace_back(
      y,
      BoolExpr::And(
          resetDeasserted,
          BoolExpr::Or(BoolExpr::Var(a), BoolExpr::Var(b))));
  problem.bad = BoolExpr::Var(y);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  // Force projected frame repair to see one omitted blocker, then cap it so
  // the same obligation immediately retries with exact frame clauses. This
  // protects the BlackParrot case where projected repair kept adding many local
  // blockers for the same obligation before reaching the exact retry.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar clauseLimit(
      "KEPLER_SEC_PDR_PROJECTED_FRAME_CLAUSE_LIMIT", "1");
  const ScopedEnvVar refinementLimit(
      "KEPLER_SEC_PDR_PROJECTED_FRAME_REFINEMENT_LIMIT", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      1,
      1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/50);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("projected-frame refinement cap reached"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCachedFallbackAvoidsRepeatedDualRailProjectedBadCube) {
  KInductionProblem problem;
  constexpr size_t a = 2;
  constexpr size_t b = 3;
  constexpr size_t input = 4;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {a, b};
  problem.inputSymbols = {input};
  problem.allSymbols = {a, b, input};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(a)),
      BoolExpr::Not(BoolExpr::Var(b)));
  problem.initialStateAssignments = {{a, false}, {b, false}};
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0.emplace_back(a, BoolExpr::Var(input));
  problem.transitions0.emplace_back(b, BoolExpr::Not(BoolExpr::Var(input)));
  problem.bad = BoolExpr::And(BoolExpr::Var(a), BoolExpr::Var(b));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  // The bad cube (a=1,b=1) is unreachable because one input drives opposite
  // next-state values.  Each single literal is reachable, so PDR must learn the
  // two-literal blocker.  The cached frame-clause fallback should now block the
  // cube instead of letting the projected bad query rediscover it forever.
  const ScopedEnvVar literalLimit(
      "KEPLER_SEC_PDR_PROJECTED_FRAME_LITERAL_LIMIT", "1");
  const ScopedEnvVar repeatedBadCubeLimit(
      "KEPLER_SEC_PDR_REPEATED_PROJECTED_BAD_CUBE_LIMIT", "2");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("bad cube cached frame clauses added="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("source=frame_log"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("repeated projected bad cube exhausted"),
      std::string::npos)
      << stderrOutput;
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
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0);
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
       PDREngineDualRailHugeStateSurfaceAvoidsRetainedBadCubeCache) {
  KInductionProblem problem;
  constexpr size_t state = 2;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  // Ariane has a multi-million-bit rail surface. Model only the cache-policy
  // signal here so the unit test stays tiny while still protecting that shape.
  problem.totalStateCount = 300000;
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(state));
  problem.initialStateAssignments = {{state, false}};
  problem.initializedStateCount = 1;
  problem.transitions0.emplace_back(state, BoolExpr::Var(state));
  problem.observedOutputExprs0 = {BoolExpr::Var(state)};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"huge_state_bad_cube_cache"};
  problem.originalObservedOutputCount = 278;
  problem.bad = BoolExpr::Var(state);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/PDREngine::kDefaultPredecessorProjectionLimit,
      /*preciseBadCubeStateLimit=*/PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("bad cube cached solver disabled"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("bad cube cached frame clauses"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineFallsBackWhenStructuralBadCubeIsEmpty) {
  KInductionProblem problem;
  constexpr size_t a = 2;
  constexpr size_t b = 3;
  constexpr size_t input = 4;
  problem.state0Symbols = {a, b};
  problem.inputSymbols = {input};
  problem.allSymbols = {a, b, input};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::Not(BoolExpr::Var(a)),
      BoolExpr::Not(BoolExpr::Var(b)));
  problem.initialStateAssignments = {{a, false}, {b, false}};
  problem.initializedStateCount = 2;
  problem.totalStateCount = 2;
  problem.transitions0.emplace_back(a, BoolExpr::Var(a));
  problem.transitions0.emplace_back(b, BoolExpr::Var(b));
  problem.bad = BoolExpr::Or(
      BoolExpr::Var(input),
      BoolExpr::And(BoolExpr::Var(a), BoolExpr::Var(b)));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_NE(
      stderrOutput.find("source=structural_model_fallback cube=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("source=structural cube=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRefinesProjectedCounterexampleWithBoundedReachability) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = {2, 3, 4};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::Not(BoolExpr::Var(4)));
  problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}};
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(
      4, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  problem.bad = BoolExpr::Var(4);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  // The full predecessor of z=1 needs both x=1 and y=1, but a one-literal
  // projected obligation may keep only x=1.  Since x=1 is reachable while y=1
  // is not, PDR must refine the spurious bounded path instead of reporting a
  // counterexample for the projected cube.
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCanDeferProjectedCounterexampleValidationToCaller) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4};
  problem.allSymbols = {2, 3, 4};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::Not(BoolExpr::Var(4)));
  problem.initialStateAssignments = {{2, false}, {3, false}, {4, false}};
  problem.initializedStateCount = 3;
  problem.totalStateCount = 3;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(
      4, BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)));
  problem.bad = BoolExpr::Var(4);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  // SEC strategy validates every PDR "Different" result with concrete BMC.
  // Its projected precision stages can therefore return the abstract candidate
  // immediately instead of doing the same bounded-prefix validation inside PDR.
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/1,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/false);
  const auto result = engine.run(3);

  EXPECT_EQ(result.status, PDRStatus::Different);
  EXPECT_EQ(result.bound, 2u);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineGeneralizesUnreachableProjectedCounterexampleRoot) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                    BoolExpr::Not(BoolExpr::Var(5))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false}, {5, false}};
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.bad = BoolExpr::And(BoolExpr::Var(5), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  // The projected predecessor can keep only the reachable x=1 literal and
  // therefore reaches Init abstractly. The concrete bad root still includes
  // b=1, which is unreachable because y is permanently false. PDR should learn
  // the exact bounded generalization b=0 rather than refining just b=1,z=1.
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("refined projected counterexample bad_frame=2 root_cube=2->1"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailGeneralizesUnreachableProjectedCounterexampleRoot) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {2, 3, 4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                    BoolExpr::Not(BoolExpr::Var(5))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false}, {5, false}};
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.bad = BoolExpr::And(BoolExpr::Var(5), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/4,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Dual-rail final SEC/PDR validates projected roots exactly.  Keep bounded
  // root generalization available so the repair learns a useful exact clause
  // instead of enumerating every sibling full-cube root one by one.
  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("refined projected counterexample bad_frame=2 root_cube=2->1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("concrete cube reachability begin cube=0"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("refined projected counterexample bad_frame=2 root_cube=2->2 checks=0"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineBudgetsDualRailProjectedCounterexampleRepairs) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {2, 3, 4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                    BoolExpr::Not(BoolExpr::Var(5))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false}, {5, false}};
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.bad = BoolExpr::And(BoolExpr::Var(5), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/4,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false,
      /*maxProjectedCounterexampleRefinements=*/1);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive);
  EXPECT_NE(
      stderrOutput.find(
          "projected counterexample repair budget exhausted refinement_limit=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCanRefineProjectedCounterexampleWithoutRootGeneralization) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                    BoolExpr::Not(BoolExpr::Var(5))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false}, {5, false}};
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.bad = BoolExpr::And(BoolExpr::Var(5), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "concrete cube reachability begin cube=2 max_step=2 "
          "mode=one_shot_unit_clauses"),
      std::string::npos);
  EXPECT_EQ(
      stderrOutput.find(
          "concrete cube reachability begin cube=2 max_step=2 "
          "mode=cached_assumptions"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("refined projected counterexample bad_frame=2 root_cube=2->2 checks=0"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineAvoidsOneShotConcreteFrameBmcAfterResetShortcutSat) {
  KInductionProblem problem;
  constexpr size_t x = 2;
  constexpr size_t y = 3;
  constexpr size_t gate = 4;
  constexpr size_t badState = 5;
  constexpr size_t reset = 6;
  problem.state0Symbols = {x, y, gate, badState};
  problem.inputSymbols = {reset};
  problem.allSymbols = {x, y, gate, badState, reset};
  problem.resetBootstrapCycles = 1;
  problem.resetBootstrapInputs = {{reset, true}};
  problem.transitions0.emplace_back(x, BoolExpr::Not(BoolExpr::Var(reset)));
  problem.transitions0.emplace_back(y, BoolExpr::createFalse());
  problem.transitions0.emplace_back(gate, BoolExpr::Not(BoolExpr::Var(reset)));
  problem.transitions0.emplace_back(
      badState,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(x), BoolExpr::Var(y)),
          BoolExpr::Var(gate)));
  problem.bad = BoolExpr::And(BoolExpr::Var(badState), BoolExpr::Var(gate));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  const ScopedEnvVar resetDiag("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/false,
      /*useExactResetFrontierChecks=*/false);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "reset-specialized expression miss reason=sat"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("one-shot cube coi post_bootstrap_steps=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineCanLearnValidatedBadFormulaClausesAfterRejectedTrace) {
  KInductionProblem problem;
  problem.state0Symbols = {2, 3, 4, 5};
  problem.allSymbols = {2, 3, 4, 5};
  problem.initialCondition = BoolExpr::And(
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(2)),
                    BoolExpr::Not(BoolExpr::Var(3))),
      BoolExpr::And(BoolExpr::Not(BoolExpr::Var(4)),
                    BoolExpr::Not(BoolExpr::Var(5))));
  problem.initialStateAssignments = {
      {2, false}, {3, false}, {4, false}, {5, false}};
  problem.initializedStateCount = 4;
  problem.totalStateCount = 4;
  problem.transitions0.emplace_back(2, BoolExpr::createTrue());
  problem.transitions0.emplace_back(3, BoolExpr::createFalse());
  problem.transitions0.emplace_back(4, BoolExpr::createTrue());
  problem.transitions0.emplace_back(
      5,
      BoolExpr::And(
          BoolExpr::And(BoolExpr::Var(2), BoolExpr::Var(3)),
          BoolExpr::Var(4)));
  problem.bad = BoolExpr::And(BoolExpr::Var(5), BoolExpr::Var(4));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/1,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailLearnsValidatedBadFormulaWithRailExpandedSupport) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  BoolExpr* bad = BoolExpr::createTrue();

  for (size_t offset = 0; offset < 12; ++offset) {
    const size_t symbol = 2 + offset;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    bad = BoolExpr::And(bad, BoolExpr::Var(symbol));

    // Keep one rail permanently false.  The full bad predicate is unreachable,
    // but a projected PDR cube can still need the validated-bad-formula repair.
    problem.transitions0.emplace_back(
        symbol,
        offset == 1 ? BoolExpr::createFalse() : BoolExpr::createTrue());
  }

  problem.usesDualRailStateEncoding = true;
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.observedOutputExprs0 = {bad};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"rail_expanded_output"};
  problem.bad = BoolExpr::simplify(bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailUsesStructuredInitWhenInitialFormulaIsPlaceholder) {
  KInductionProblem problem;
  constexpr size_t state = 2;

  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {state};
  problem.allSymbols = {state};
  problem.initialStateAssignments = {{state, false}};
  problem.initializedStateCount = 1;
  problem.totalStateCount = 1;
  // Dual-rail extraction keeps the real startup facts structured and uses this
  // non-null formula only to select the structured-init encoder path.
  problem.initialCondition = BoolExpr::createTrue();
  problem.transitions0.emplace_back(state, BoolExpr::createFalse());
  problem.bad = BoolExpr::Var(state);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
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
       PDREngineDualRailPredecessorBudgetReturnsInconclusive) {
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
  // Init excludes the target bad state, while F1 can still propose it as a PDR
  // obligation. With a zero-decision cap, the otherwise SAT predecessor query
  // must report UNKNOWN and leave the local proof slice inconclusive.
  problem.initialCondition = BoolExpr::Not(BoolExpr::Var(targetState));
  problem.bad = BoolExpr::Var(targetState);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor query budget exhausted"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailPredecessorReusesCachedFallback) {
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
  problem.observedOutputNames = {"single_output_cached_retry"};
  problem.originalObservedOutputCount = 1266;

  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  // The fallback should spend its retry budget in the already encoded cached
  // predecessor solver instead of reconstructing the same frame and transition
  // cone as a fresh SAT instance.
  EXPECT_NE(
      stderrOutput.find("cached_assumptions=unknown retry=cached_solver"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("cached_solver_retry=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailAesSizedLeafUsesReferencePredecessorFallback) {
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
  problem.observedOutputNames = {"single_output_aes_sized_fallback"};
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar decisionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_PREDECESSOR_DECISION_LIMIT", "0");
  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true);
  const auto result = engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  // AES residual leaves have a one-output local shape after splitting, but the
  // original design is still inside the medium-output guard. Keep the
  // 376a017-style fresh predecessor fallback and do not retain the cached retry
  // solver that caused available-memory spikes on AES.
  EXPECT_NE(
      stderrOutput.find("cached_assumptions=unknown fallback=exact"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("cached_assumptions=unknown retry=cached_solver"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("cached_solver_retry=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor result cache hit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor cached core"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailAesSizedSatLeafUsesFreshPredecessorFallback) {
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
  problem.observedOutputNames = {"single_output_aes_sized_sat_fallback"};
  problem.originalObservedOutputCount = 129;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // The good AES path used the cached solver only as a cheap probe.  A cached
  // SAT answer still falls through to the ordinary exact predecessor solver so
  // AES-sized leaves do not keep the broad residual-leaf cached model path.
  EXPECT_NE(
      stderrOutput.find("cached_assumptions=sat fallback=exact"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("result=sat cached_assumptions=1"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor result cache hit"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor cached core"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailSingleOutputResidualRaisesPredecessorBudget) {
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
  problem.observedOutputNames = {"single_output_residual"};

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Final single-output dual-rail repairs are still ordinary PDR predecessor
  // checks, and the residual predecessor budget must stay at the original
  // proof-search bound. Runtime fixes should reduce rebuild cost instead of
  // shrinking this legal PDR search budget.
  EXPECT_NE(
      stderrOutput.find("conflict_limit=200000"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailSingleOutputUsesStableCachedPredecessorSurface) {
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
  problem.observedOutputNames = {"single_output_stable_cache"};
  problem.originalObservedOutputCount = 1266;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  // Single-output dual-rail leaves should build the reusable cached
  // predecessor solver on the stable local surface. The unrelated decoy state
  // must not be pulled in merely because this is a dual-rail leaf.
  EXPECT_NE(
      stderrOutput.find("solver_symbols=3 cached_solver_symbols=3"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "predecessor cached solver created level=0 symbols=3"),
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
      stderrOutput.find("predecessor closed symbol cache seed="),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor target surface cached"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailHugeStateSurfaceAvoidsRetainedPredecessorCaches) {
  KInductionProblem problem;
  constexpr size_t targetState = 2;
  constexpr size_t stateA = 3;
  constexpr size_t stateB = 4;
  constexpr size_t decoyState = 5;
  problem.usesDualRailStateEncoding = true;
  problem.state0Symbols = {targetState, stateA, stateB, decoyState};
  problem.allSymbols = {targetState, stateA, stateB, decoyState};
  // Model Ariane's multi-million-bit rail surface without allocating it in the
  // unit test. The real query surface stays tiny, but the PDR cache policy must
  // still choose the low-retention path for this shape.
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
  problem.observedOutputNames = {"huge_state_uncached_surface"};
  problem.originalObservedOutputCount = 278;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar pdrStatsInterval("KEPLER_SEC_PDR_STATS_INTERVAL", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      PDREngine::kDefaultPredecessorProjectionLimit,
      PDREngine::kDefaultPreciseBadCubeStateLimit,
      /*useExactFrameClauses=*/true);
  (void)engine.run(1);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("predecessor target surface uncached"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("predecessor cached solver disabled"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor target surface cached"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("predecessor cached solver created"),
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
       PDREngineDualRailLearnsPerOutputBadFormulaPastBinaryClauseLimit) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<BoolExpr*> outputs;
  size_t nextSymbol = 2;

  for (size_t output = 0; output < 2; ++output) {
    std::vector<size_t> symbols;
    symbols.reserve(8);
    for (size_t bit = 0; bit < 8; ++bit) {
      const size_t symbol = nextSymbol++;
      symbols.push_back(symbol);
      problem.state0Symbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
      problem.initialStateAssignments.push_back({symbol, false});
      init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
      problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
    }
    // OR over eight rails produces 255 state-only bad assignments.  Binary SEC
    // intentionally keeps that above its eager repair cap; dual rail needs it
    // because a small unknown/known output can enumerate many Boolean rail
    // combinations.
    outputs.push_back(BoolExpr::simplify(makeOrChain(symbols)));
    problem.observedOutputNames.push_back("rail_or_" + std::to_string(output));
  }

  problem.usesDualRailStateEncoding = true;
  problem.observedOutputExprs0 = outputs;
  problem.observedOutputExprs1.assign(outputs.size(), BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::createFalse();
  for (auto* outputExpr : outputs) {
    problem.bad = BoolExpr::Or(problem.bad, outputExpr);
  }
  problem.bad = BoolExpr::simplify(problem.bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 2)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar proofConflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_BAD_CUBE_CONFLICT_LIMIT", "5000");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "per-output validated bad-formula clauses bad_frame=1 output=0 "
          "outputs=1 "
          "clauses=255"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find(
          "skipped per-output bad-formula validation bad_frame=1 output=0 "
          "clauses=255 reason=clause_limit"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailKeepsPerOutputGroupsPastTotalClauseCap) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<BoolExpr*> outputs;
  size_t nextSymbol = 2;

  for (size_t output = 0; output < 10; ++output) {
    std::vector<size_t> symbols;
    symbols.reserve(9);
    for (size_t bit = 0; bit < 9; ++bit) {
      const size_t symbol = nextSymbol++;
      symbols.push_back(symbol);
      problem.state0Symbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
      problem.initialStateAssignments.push_back({symbol, false});
      init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
      problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
    }
    outputs.push_back(BoolExpr::simplify(makeOrChain(symbols)));
    problem.observedOutputNames.push_back("rail_group_" + std::to_string(output));
  }

  problem.usesDualRailStateEncoding = true;
  problem.observedOutputExprs0 = outputs;
  problem.observedOutputExprs1.assign(outputs.size(), BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::createFalse();
  for (auto* outputExpr : outputs) {
    problem.bad = BoolExpr::Or(problem.bad, outputExpr);
  }
  problem.bad = BoolExpr::simplify(problem.bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar proofConflictLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_BAD_CUBE_CONFLICT_LIMIT", "5000");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "per-output validated bad-formula clauses bad_frame=1 output=9 "
          "outputs=10"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineLearnsValidatedBadFormulaClausesPerOutputInBatchedSecSlice) {
  KInductionProblem problem;
  const std::vector<std::vector<size_t>> outputStateGroups = {
      {2, 4, 5, 6, 7, 8},
      {9, 11, 12, 13, 14, 15}};

  BoolExpr* init = BoolExpr::createTrue();
  auto makeConjunction = [](const std::vector<size_t>& symbols) {
    BoolExpr* expr = BoolExpr::createTrue();
    for (const auto symbol : symbols) {
      expr = BoolExpr::And(expr, BoolExpr::Var(symbol));
    }
    return BoolExpr::simplify(expr);
  };

  auto addState = [&](size_t symbol, BoolExpr* next) {
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, next);
  };
  auto addOutputGroup = [&](size_t base) {
    addState(base + 0, BoolExpr::createTrue());
    addState(base + 1, BoolExpr::createFalse());
    addState(base + 2, BoolExpr::createTrue());
    addState(
        base + 3,
        BoolExpr::And(
            BoolExpr::And(BoolExpr::Var(base + 0), BoolExpr::Var(base + 1)),
            BoolExpr::Var(base + 2)));
    addState(base + 4, BoolExpr::createTrue());
    addState(base + 5, BoolExpr::createTrue());
    addState(base + 6, BoolExpr::createTrue());
  };

  addOutputGroup(2);
  addOutputGroup(9);
  for (const auto& group : outputStateGroups) {
    for (const auto symbol : group) {
      // Each per-output bad predicate is small enough to learn directly, but
      // the combined batched bad predicate has 12 state symbols. This guards
      // the BlackParrot case where the useful refinement is per observed output
      // rather than over the whole output-batch support union.
      ASSERT_NE(
          std::find(
              problem.state0Symbols.begin(), problem.state0Symbols.end(), symbol),
          problem.state0Symbols.end());
    }
  }

  BoolExpr* output0 = makeConjunction(outputStateGroups[0]);
  BoolExpr* output1 = makeConjunction(outputStateGroups[1]);
  problem.observedOutputExprs0 = {output0, output1};
  problem.observedOutputExprs1 = {
      BoolExpr::createFalse(), BoolExpr::createFalse()};
  problem.observedOutputNames = {"o0", "o1"};
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = 14;
  problem.totalStateCount = 14;
  problem.bad = BoolExpr::simplify(BoolExpr::Or(output0, output1));
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(stderrOutput.find(" clauses=2"), std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineRepairsBroadValidatedBadFormulaLearningPerOutput) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  auto makeConjunction = [](const std::vector<size_t>& symbols) {
    BoolExpr* expr = BoolExpr::createTrue();
    for (const auto symbol : symbols) {
      expr = BoolExpr::And(expr, BoolExpr::Var(symbol));
    }
    return BoolExpr::simplify(expr);
  };

  std::vector<BoolExpr*> outputs;
  size_t nextSymbol = 2;
  for (size_t output = 0; output < 10; ++output) {
    const size_t base = nextSymbol;
    nextSymbol += 7;
    std::vector<size_t> group;
    group.reserve(6);
    for (size_t offset = 0; offset < 6; ++offset) {
      const size_t symbol = base + offset;
      group.push_back(symbol);
      problem.state0Symbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
      problem.initialStateAssignments.push_back({symbol, false});
      init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    }

    // One permanently false bit keeps every per-output conjunction unreachable,
    // while the projected bad cube remains small enough to exercise PDR's
    // concrete root-cube refinement after broad bad-formula validation is
    // deliberately skipped.
    problem.transitions0.emplace_back(base + 0, BoolExpr::createTrue());
    problem.transitions0.emplace_back(base + 1, BoolExpr::createFalse());
    problem.transitions0.emplace_back(base + 2, BoolExpr::createTrue());
    problem.transitions0.emplace_back(base + 3, BoolExpr::createTrue());
    problem.transitions0.emplace_back(base + 4, BoolExpr::createTrue());
    problem.transitions0.emplace_back(base + 5, BoolExpr::createTrue());

    outputs.push_back(makeConjunction(group));
    problem.observedOutputNames.push_back("o" + std::to_string(output));
  }

  problem.observedOutputExprs0 = outputs;
  problem.observedOutputExprs1.assign(outputs.size(), BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::createFalse();
  for (auto* outputExpr : outputs) {
    problem.bad = BoolExpr::Or(problem.bad, outputExpr);
  }
  problem.bad = BoolExpr::simplify(problem.bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  ASSERT_FALSE(
      findBaseCounterexample(
          problem, KEPLER_FORMAL::Config::SolverType::KISSAT, 3)
          .has_value());

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/false,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent);
  EXPECT_NE(
      stderrOutput.find("skipped broad bad-formula validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "per-output validated bad-formula clauses bad_frame=1 output=0 "
          "outputs=1 clauses=1"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailExactMultiOutputAvoidsWholeBatchBadFormulaRepair) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<BoolExpr*> outputs;
  size_t nextSymbol = 2;

  for (size_t output = 0; output < 17; ++output) {
    const size_t base = nextSymbol;
    nextSymbol += 2;
    problem.state0Symbols.push_back(base);
    problem.state0Symbols.push_back(base + 1);
    problem.allSymbols.push_back(base);
    problem.allSymbols.push_back(base + 1);
    problem.initialStateAssignments.push_back({base, false});
    problem.initialStateAssignments.push_back({base + 1, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(base)));
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(base + 1)));
    problem.transitions0.emplace_back(base, BoolExpr::createTrue());
    problem.transitions0.emplace_back(base + 1, BoolExpr::createFalse());
    outputs.push_back(BoolExpr::And(BoolExpr::Var(base), BoolExpr::Var(base + 1)));
    problem.observedOutputNames.push_back("rail_batch_" + std::to_string(output));
  }

  // Keep the state surface above the multi-output defer threshold. The test is
  // about policy selection, so these filler states are intentionally outside
  // the output cones.
  while (problem.state0Symbols.size() < 520) {
    const size_t symbol = nextSymbol++;
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, BoolExpr::Var(symbol));
  }

  problem.usesDualRailStateEncoding = true;
  problem.observedOutputExprs0 = outputs;
  problem.observedOutputExprs1.assign(outputs.size(), BoolExpr::createFalse());
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = BoolExpr::createFalse();
  for (auto* outputExpr : outputs) {
    problem.bad = BoolExpr::Or(problem.bad, outputExpr);
  }
  problem.bad = BoolExpr::simplify(problem.bad);
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(3);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("skipped broad bad-formula validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("batched reset-cube validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("trying deep whole bad-formula validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("refined projected counterexample with validated "
                        "bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailLargeFrontierSkipsResetCubeBadFormulaRepair) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<size_t> outputSymbols;
  outputSymbols.reserve(8);

  for (size_t symbol = 2; symbol < 10; ++symbol) {
    outputSymbols.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }
  // Keep this above the production dual-rail reset-frontier cap.  The default
  // was raised for Ibex-sized proofs, so this regression must scale with it
  // instead of forcing the real flow back to the smaller historical limit.
  for (size_t index = 0; index < 10001; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{10000 + index * 2, 10001 + index * 2});
  }

  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.observedOutputExprs0 = {BoolExpr::simplify(makeOrChain(outputSymbols))};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"wide_rail_leaf"};
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = problem.observedOutputExprs0.front();
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("skipped deep bad-formula base validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("validated bad-formula clauses with reset cubes"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-frontier bad-formula proof"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailLargeTransitionSurfaceSkipsResetCubeRepair) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<size_t> outputSymbols;
  outputSymbols.reserve(8);

  for (size_t symbol = 2; symbol < 10; ++symbol) {
    outputSymbols.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  for (size_t index = 0; index < 20001; ++index) {
    problem.lazyTransitions->sourceByStateSymbol.emplace(
        20000 + index,
        LazyTransitionSource{0, BoolExpr::createFalse(), LazyTransitionRail::Binary});
  }

  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.observedOutputExprs0 = {BoolExpr::simplify(makeOrChain(outputSymbols))};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"wide_lazy_surface_leaf"};
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = problem.observedOutputExprs0.front();
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2, /*resetBootstrapFrameCheckedSafe=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("skipped deep bad-formula base validation"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("validated bad-formula clauses with reset cubes"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailSkipsResetCubeBadFormulaRepair) {
  KInductionProblem problem;
  BoolExpr* init = BoolExpr::createTrue();
  std::vector<size_t> outputSymbols;
  outputSymbols.reserve(8);

  for (size_t symbol = 2; symbol < 10; ++symbol) {
    outputSymbols.push_back(symbol);
    problem.state0Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.initialStateAssignments.push_back({symbol, false});
    init = BoolExpr::And(init, BoolExpr::Not(BoolExpr::Var(symbol)));
    problem.transitions0.emplace_back(symbol, BoolExpr::createFalse());
  }

  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.observedOutputExprs0 = {BoolExpr::simplify(makeOrChain(outputSymbols))};
  problem.observedOutputExprs1 = {BoolExpr::createFalse()};
  problem.observedOutputNames = {"small_rail_leaf"};
  problem.initialCondition = BoolExpr::simplify(init);
  problem.initializedStateCount = problem.state0Symbols.size();
  problem.totalStateCount = problem.state0Symbols.size();
  problem.bad = problem.observedOutputExprs0.front();
  problem.property = BoolExpr::Not(problem.bad);
  problem.inductionProperty = problem.property;
  problem.inductionBad = problem.bad;

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(
      problem,
      KEPLER_FORMAL::Config::SolverType::KISSAT,
      /*predecessorProjectionLimit=*/2,
      /*preciseBadCubeStateLimit=*/2,
      /*useExactFrameClauses=*/true,
      /*maxPredecessorQueries=*/0,
      /*refineProjectedCounterexamples=*/true,
      /*maxBoundedRootGeneralizationAttempts=*/0,
      /*learnValidatedBadFormulaClauses=*/true);
  const auto result = engine.run(2, /*resetBootstrapFrameCheckedSafe=*/true);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Equivalent) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find(
          "refined projected counterexample with validated bad-formula clauses"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("validated bad-formula clauses with reset cubes"),
      std::string::npos)
      << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("reset-frontier bad-formula proof"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailLargeResetBootstrapSkipsBmcPrecheck) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.initialCondition = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.property = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();
  problem.inductionProperty = BoolExpr::createTrue();
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  for (size_t index = 0; index < 8193; ++index) {
    problem.lazyTransitions->sourceByStateSymbol.emplace(
        20000 + index,
        LazyTransitionSource{0, BoolExpr::createFalse(), LazyTransitionRail::Binary});
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_NE(
      stderrOutput.find("skipped dual-rail reset-bootstrap BMC precheck"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailResetBootstrapPrecheckHonorsTransitionLimitOverride) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.initialCondition = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.property = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();
  problem.inductionProperty = BoolExpr::createTrue();
  problem.lazyTransitions = std::make_shared<LazyTransitionStore>();
  for (size_t index = 0; index < 10001; ++index) {
    problem.lazyTransitions->sourceByStateSymbol.emplace(
        20000 + index,
        LazyTransitionSource{
            0, BoolExpr::createFalse(), LazyTransitionRail::Binary});
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  const ScopedEnvVar transitionLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_RESET_BMC_TRANSITION_SOURCE_LIMIT", "10001");
  testing::internal::CaptureStderr();
  PDREngine engine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  const auto result = engine.run(0);
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, PDRStatus::Inconclusive) << stderrOutput;
  EXPECT_EQ(
      stderrOutput.find("skipped dual-rail reset-bootstrap BMC precheck"),
      std::string::npos)
      << stderrOutput;
}

TEST_F(SequentialEquivalenceStrategyTests,
       PDREngineDualRailResetFrontierHonorsStateSymbolLimitOverride) {
  KInductionProblem problem;
  problem.usesDualRailStateEncoding = true;
  problem.resetBootstrapCycles = 1;
  problem.initialCondition = BoolExpr::createTrue();
  problem.bad = BoolExpr::createFalse();
  problem.property = BoolExpr::createTrue();
  problem.inductionBad = BoolExpr::createFalse();
  problem.inductionProperty = BoolExpr::createTrue();
  for (size_t index = 0; index < 10001; ++index) {
    problem.dualRailStatePairs.push_back(
        DualRailSymbolPair{20000 + index * 2, 20001 + index * 2});
  }

  const ScopedEnvVar pdrStats("KEPLER_SEC_PDR_STATS", "1");
  testing::internal::CaptureStderr();
  PDREngine defaultEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_NE(
      stderrOutput.find("exact reset-frontier checks disabled"),
      std::string::npos)
      << stderrOutput;

  const ScopedEnvVar stateLimit(
      "KEPLER_SEC_PDR_DUAL_RAIL_RESET_FRONTIER_STATE_SYMBOL_LIMIT", "20002");
  testing::internal::CaptureStderr();
  PDREngine overriddenEngine(problem, KEPLER_FORMAL::Config::SolverType::KISSAT);
  stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(
      stderrOutput.find("exact reset-frontier checks disabled"),
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

  const auto extracted = SequentialDesignModel::extract(top);
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
  const size_t stateVar = extracted.inputVarByKey.at(stateKey);
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
       TooSmallBoundRemainsInconclusiveBeforeCounterexampleDepth) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* top0 = createResetInitializedPipelineTop(library, "top0", false);
  auto* top1 = createResetInitializedPipelineTop(library, "top1", true);

  auto strategy = makeBinarySecStrategy(top0, top1, SecEngine::KInduction);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  auto strategy = makeBinarySecStrategy(top0, top1);
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

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
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

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_EQ(result.coveredOutputs, 1u);
  EXPECT_EQ(result.totalOutputs, 2u);
  ASSERT_EQ(result.skippedObservedOutputs.size(), 1u);
  EXPECT_NE(result.skippedObservedOutputs.front().find("bad[0]"), std::string::npos);
  EXPECT_NE(
      result.skippedObservedOutputs.front().find("logical-loop"),
      std::string::npos);
}

TEST_F(SequentialEquivalenceStrategyTests,
       EquivalentDffDesignsReportTopBoundarySurface) {
  NLUniverse::create();
  auto* db = NLDB::create(NLUniverse::get());
  auto* primitives =
      NLLibrary::create(db, NLLibrary::Type::Primitives, NLName("prims"));
  auto* library =
      NLLibrary::create(db, NLLibrary::Type::Standard, NLName("designs"));
  auto* invModel = createInvModel(primitives);
  auto* top0 = createDffTop(library, "top0", invModel, false, false);
  auto* top1 = createDffTop(library, "top1", invModel, false, false);

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
  testing::internal::CaptureStdout();
  testing::internal::CaptureStderr();
  auto strategy = makeBinarySecStrategy(top0, top1);
  const auto result = strategy.run(3);
  const std::string stdoutOutput = testing::internal::GetCapturedStdout();
  const std::string stderrOutput = testing::internal::GetCapturedStderr();

  EXPECT_EQ(result.status, SequentialEquivalenceStatus::Equivalent);
  EXPECT_NE(stderrOutput.find("SEC diag: start run"), std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: extract(top0) collect begin"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: deferred next-state formula remapping"),
      std::string::npos);
  EXPECT_NE(
      stderrOutput.find("SEC diag: entering pdr engine"),
      std::string::npos);
  EXPECT_NE(stdoutOutput.find("SEC diag: aligned_inputs="), std::string::npos);
  EXPECT_NE(stdoutOutput.find("SEC summary: property_is_true="), std::string::npos);
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
  EXPECT_EQ(
      detail::detectInitialStateValueForTest({{"R", 11}, {"S", 12}}),
      std::nullopt);
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
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("reset_i[0]"),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("rst_ni[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("rst_l[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("reset_l[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("wb_rst_i[0]"),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("wb_reset_i[0]"),
      std::optional<bool>(true));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("wb_rst_ni[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("rrst_n[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("wrst_n[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(
      detail::getResetAssertionValueForTest("aresetn[0]"),
      std::optional<bool>(false));
  EXPECT_EQ(detail::getResetAssertionValueForTest("burst_n[0]"), std::nullopt);
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
