// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "strategy/StructuralStateInvariant.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BoolExpr.h"

namespace KEPLER_FORMAL::SEC {

// Overall structural-state matching algorithm:
// 1. Put aligned SEC inputs into a shared abstract symbol space.
// 2. Seed each state bit with coarse classes from init/complement information.
// 3. Repeatedly fingerprint next-state functions under the current classes.
// 4. Refine those classes to a fixed point.
// 5. Pair states whose final structural fingerprints match across the designs.
// 6. Use an ordered fast path only when every state already matches in order.

namespace {

using KEPLER_FORMAL::BoolExpr;

struct ExprPairHash {
  size_t operator()(const std::pair<BoolExpr*, BoolExpr*>& pair) const noexcept {
    size_t seed = std::hash<const void*>()(pair.first);
    seed ^= std::hash<const void*>()(pair.second) + 0x9e3779b9 + (seed << 6) +
            (seed >> 2);
    return seed;
  }
};

bool areEquivalentUnderAbstractMapsImpl(
    BoolExpr* expr0,
    BoolExpr* expr1,
    const LocalToAbstractVarMap& abstractMap0,
    const LocalToAbstractVarMap& abstractMap1,
    std::unordered_map<std::pair<BoolExpr*, BoolExpr*>, bool, ExprPairHash>& memo) {
  if (expr0 == nullptr || expr1 == nullptr) {
    return expr0 == expr1;
  }

  const auto key = std::make_pair(expr0, expr1);
  if (const auto it = memo.find(key); it != memo.end()) {
    return it->second;
  }

  bool equivalent = false;
  if (expr0->getOp() == Op::VAR && expr1->getOp() == Op::VAR) {
    if (expr0->getId() < 2 || expr1->getId() < 2) {
      equivalent = expr0->getId() == expr1->getId();
    } else {
      const auto it0 = abstractMap0.find(expr0->getId());
      const auto it1 = abstractMap1.find(expr1->getId());
      equivalent = it0 != abstractMap0.end() && it1 != abstractMap1.end() &&
                   it0->second == it1->second;
    }
  } else if (expr0->getOp() == expr1->getOp()) {
    equivalent = areEquivalentUnderAbstractMapsImpl(
                     expr0->getLeft(),
                     expr1->getLeft(),
                     abstractMap0,
                     abstractMap1,
                     memo) &&
                 areEquivalentUnderAbstractMapsImpl(
                     expr0->getRight(),
                     expr1->getRight(),
                     abstractMap0,
                     abstractMap1,
                     memo);
  }

  memo.emplace(key, equivalent);
  return equivalent;
}

std::pair<LocalToAbstractVarMap, LocalToAbstractVarMap> buildAbstractTransitionMapsImpl(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedStates) {
  // Matched SEC inputs and already-correlated state bits get the same abstract
  // IDs on both sides. Everything else stays private so equivalence checks only
  // identify what has really been aligned.
  LocalToAbstractVarMap abstractMap0;
  LocalToAbstractVarMap abstractMap1;
  size_t nextAbstractSymbol = 2;

  for (size_t i = 0; i < alignedInputs.names.size(); ++i) {
    const size_t symbol = nextAbstractSymbol++;
    abstractMap0.emplace(model0.inputVarByKey.at(alignedInputs.keys0[i]), symbol);
    abstractMap1.emplace(model1.inputVarByKey.at(alignedInputs.keys1[i]), symbol);
  }
  for (size_t i = 0; i < alignedStates.names.size(); ++i) {
    const size_t symbol = nextAbstractSymbol++;
    abstractMap0.emplace(model0.inputVarByKey.at(alignedStates.keys0[i]), symbol);
    abstractMap1.emplace(model1.inputVarByKey.at(alignedStates.keys1[i]), symbol);
  }

  auto assignPrivateStateSymbols = [&](const SequentialDesignModel& model,
                                       LocalToAbstractVarMap& abstractMap) {
    for (const auto& key : model.stateBits) {
      const size_t localVar = model.inputVarByKey.at(key);
      if (abstractMap.find(localVar) != abstractMap.end()) {
        continue;
      }
      abstractMap.emplace(localVar, nextAbstractSymbol++);
    }
  };
  assignPrivateStateSymbols(model0, abstractMap0);
  assignPrivateStateSymbols(model1, abstractMap1);

  return {std::move(abstractMap0), std::move(abstractMap1)};
}

AlignedSignals buildOrderedStatePairs(const SequentialDesignModel& model0,
                                      const SequentialDesignModel& model1) {
  if (model0.stateBits.size() != model1.stateBits.size()) {
    return {};
  }

  AlignedSignals aligned;
  aligned.names.reserve(model0.stateBits.size());
  aligned.keys0.reserve(model0.stateBits.size());
  aligned.keys1.reserve(model0.stateBits.size());
  for (size_t i = 0; i < model0.stateBits.size(); ++i) {
    aligned.names.push_back("ordered_state_" + std::to_string(i));
    aligned.keys0.push_back(model0.stateBits[i]);
    aligned.keys1.push_back(model1.stateBits[i]);
  }
  return aligned;
}

bool areAllOrderedStatesEquivalent(const SequentialDesignModel& model0,
                                   const SequentialDesignModel& model1,
                                   const AlignedSignals& alignedInputs,
                                   const AlignedSignals& orderedStates) {
  if (orderedStates.names.empty()) {
    return false;
  }

  const auto [abstractMap0, abstractMap1] =
      buildAbstractTransitionMapsImpl(model0, model1, alignedInputs, orderedStates);
  for (size_t i = 0; i < orderedStates.names.size(); ++i) {
    const auto& key0 = orderedStates.keys0[i];
    const auto& key1 = orderedStates.keys1[i];
    std::unordered_map<std::pair<BoolExpr*, BoolExpr*>, bool, ExprPairHash> memo;
    if (!areEquivalentUnderAbstractMapsImpl(
            model0.nextStateExprByStateKey.at(key0),
            model1.nextStateExprByStateKey.at(key1),
            abstractMap0,
            abstractMap1,
            memo)) {
      return false;
    }
  }
  return true;
}

std::unordered_map<SignalKey, size_t, SignalKeyHash> buildStateIndexMap(
    const std::vector<SignalKey>& stateBits) {
  std::unordered_map<SignalKey, size_t, SignalKeyHash> indices;
  indices.reserve(stateBits.size());
  for (size_t i = 0; i < stateBits.size(); ++i) {
    indices.emplace(stateBits[i], i);
  }
  return indices;
}

std::unordered_map<size_t, size_t> buildInputClassMap(
    const SequentialDesignModel& model,
    const std::vector<SignalKey>& alignedInputKeys) {
  std::unordered_map<size_t, size_t> classes;
  classes.reserve(alignedInputKeys.size());
  for (size_t i = 0; i < alignedInputKeys.size(); ++i) {
    classes.emplace(model.inputVarByKey.at(alignedInputKeys[i]), i);
  }
  return classes;
}

std::unordered_map<SignalKey, char, SignalKeyHash> buildComplementRoleMap(
    const SequentialDesignModel& model) {
  std::unordered_map<SignalKey, char, SignalKeyHash> roles;
  for (const auto& key : model.stateBits) {
    roles.emplace(key, 'N');
  }
  for (const auto& relation : model.complementedStateRelations) {
    roles[relation.primaryKey] = 'P';
    roles[relation.complementedKey] = 'C';
  }
  return roles;
}

uint64_t mixHash(uint64_t value) {
  value += 0x9e3779b97f4a7c15ull;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
  return value ^ (value >> 31);
}

uint64_t combineHashes(std::initializer_list<uint64_t> parts) {
  uint64_t hash = 0x243f6a8885a308d3ull;
  for (const auto part : parts) {
    hash = mixHash(hash ^ mixHash(part));
  }
  return hash;
}

uint64_t initialStateSignature(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, char, SignalKeyHash>& complementRoles,
    const SignalKey& key) {
  uint64_t initHash = 0;
  if (const auto initIt = model.initialStateValueByKey.find(key);
      initIt != model.initialStateValueByKey.end()) {
    initHash = initIt->second ? 1 : 2;
  }
  uint64_t roleHash = 0;
  if (const auto roleIt = complementRoles.find(key);
      roleIt != complementRoles.end()) {
    roleHash = static_cast<unsigned char>(roleIt->second);
  }
  return combineHashes({initHash, roleHash});
}

uint64_t fingerprintExpr(
    BoolExpr* expr,
    const std::unordered_map<size_t, size_t>& inputClasses,
    const std::unordered_map<size_t, size_t>& stateClasses,
    std::unordered_map<BoolExpr*, uint64_t>& memo) {
  if (expr == nullptr) {
    return 0;
  }
  if (const auto it = memo.find(expr); it != memo.end()) {
    return it->second;
  }

  uint64_t fingerprint = 0;
  switch (expr->getOp()) {
    case Op::VAR: {
      if (expr->getId() < 2) {
        fingerprint = expr->getId() == 0 ? 7 : 11;
        break;
      }
      if (const auto inputIt = inputClasses.find(expr->getId());
          inputIt != inputClasses.end()) {
        fingerprint = combineHashes({13, inputIt->second});
      } else if (const auto stateIt = stateClasses.find(expr->getId());
                 stateIt != stateClasses.end()) {
        fingerprint = combineHashes({17, stateIt->second});
      } else {
        fingerprint = 19;
      }
      break;
    }
    case Op::NOT:
      fingerprint = combineHashes(
          {23, fingerprintExpr(expr->getLeft(), inputClasses, stateClasses, memo)});
      break;
    case Op::AND:
    case Op::OR:
    case Op::XOR: {
      uint64_t lhs =
          fingerprintExpr(expr->getLeft(), inputClasses, stateClasses, memo);
      uint64_t rhs =
          fingerprintExpr(expr->getRight(), inputClasses, stateClasses, memo);
      if (lhs > rhs) {
        std::swap(lhs, rhs);
      }
      const uint64_t opTag =
          expr->getOp() == Op::AND ? 29 : (expr->getOp() == Op::OR ? 31 : 37);
      fingerprint = combineHashes({opTag, lhs, rhs});
      break;
    }
    case Op::NONE:
    default:
      fingerprint = 41;
      break;
  }

  memo.emplace(expr, fingerprint);
  return fingerprint;
}

std::vector<size_t> refineClasses(
    const SequentialDesignModel& model,
    const std::unordered_map<size_t, size_t>& inputClasses,
    const std::vector<size_t>& currentClasses) {
  // Replace each state's current class by a finer signature that combines the
  // previous class and the structure of its next-state expression.
  std::unordered_map<size_t, size_t> stateClasses;
  stateClasses.reserve(model.stateBits.size());
  for (size_t i = 0; i < model.stateBits.size(); ++i) {
    stateClasses.emplace(model.inputVarByKey.at(model.stateBits[i]), currentClasses[i]);
  }

  using RefinementSignature = std::pair<size_t, uint64_t>;
  std::vector<RefinementSignature> signatures(model.stateBits.size());
  std::unordered_map<BoolExpr*, uint64_t> memo;
  for (size_t i = 0; i < model.stateBits.size(); ++i) {
    signatures[i] = {
        currentClasses[i],
        fingerprintExpr(
            model.nextStateExprByStateKey.at(model.stateBits[i]),
            inputClasses,
            stateClasses,
            memo)};
  }

  std::vector<RefinementSignature> unique = signatures;
  std::sort(unique.begin(), unique.end());
  unique.erase(std::unique(unique.begin(), unique.end()), unique.end());

  std::map<RefinementSignature, size_t> classBySignature;
  for (size_t i = 0; i < unique.size(); ++i) {
    classBySignature.emplace(unique[i], i);
  }

  std::vector<size_t> refined(model.stateBits.size(), 0);
  for (size_t i = 0; i < signatures.size(); ++i) {
    refined[i] = classBySignature.at(signatures[i]);
  }
  return refined;
}

std::vector<size_t> seedClasses(const SequentialDesignModel& model) {
  const auto complementRoles = buildComplementRoleMap(model);
  std::vector<uint64_t> signatures;
  signatures.reserve(model.stateBits.size());
  for (const auto& key : model.stateBits) {
    signatures.push_back(initialStateSignature(model, complementRoles, key));
  }

  std::vector<uint64_t> unique = signatures;
  std::sort(unique.begin(), unique.end());
  unique.erase(std::unique(unique.begin(), unique.end()), unique.end());

  std::unordered_map<uint64_t, size_t> classBySignature;
  classBySignature.reserve(unique.size());
  for (size_t i = 0; i < unique.size(); ++i) {
    classBySignature.emplace(unique[i], i);
  }

  std::vector<size_t> seeded(signatures.size(), 0);
  for (size_t i = 0; i < signatures.size(); ++i) {
    seeded[i] = classBySignature.at(signatures[i]);
  }
  return seeded;
}

struct StateClassFingerprint {
  uint64_t seedSignature = 0;
  uint64_t transitionFingerprint = 0;
};

std::vector<StateClassFingerprint> computeFinalFingerprints(
    const SequentialDesignModel& model,
    const std::unordered_map<size_t, size_t>& inputClasses,
    const std::vector<size_t>& finalClasses) {
  const auto complementRoles = buildComplementRoleMap(model);

  std::unordered_map<size_t, size_t> stateClasses;
  stateClasses.reserve(model.stateBits.size());
  for (size_t i = 0; i < model.stateBits.size(); ++i) {
    stateClasses.emplace(model.inputVarByKey.at(model.stateBits[i]), finalClasses[i]);
  }

  std::vector<StateClassFingerprint> fingerprints;
  fingerprints.reserve(model.stateBits.size());
  std::unordered_map<BoolExpr*, uint64_t> memo;
  for (const auto& key : model.stateBits) {
    fingerprints.push_back(
        {initialStateSignature(model, complementRoles, key),
         fingerprintExpr(
             model.nextStateExprByStateKey.at(key),
             inputClasses,
             stateClasses,
             memo)});
  }
  return fingerprints;
}

std::string stateFingerprintKey(const StateClassFingerprint& fingerprint) {
  return std::to_string(fingerprint.seedSignature) + ":" +
      std::to_string(fingerprint.transitionFingerprint);
}

}  // namespace

std::pair<LocalToAbstractVarMap, LocalToAbstractVarMap> buildAbstractTransitionMaps(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedStates) {
  return buildAbstractTransitionMapsImpl(
      model0, model1, alignedInputs, alignedStates);
}

bool areEquivalentUnderAbstractMaps(
    BoolExpr* expr0,
    BoolExpr* expr1,
    const LocalToAbstractVarMap& abstractMap0,
    const LocalToAbstractVarMap& abstractMap1) {
  std::unordered_map<std::pair<BoolExpr*, BoolExpr*>, bool, ExprPairHash> memo;
  return areEquivalentUnderAbstractMapsImpl(
      expr0, expr1, abstractMap0, abstractMap1, memo);
}

AlignedSignals inferStructurallyEquivalentStatePairs(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs) {
  if (model0.stateBits.empty() || model1.stateBits.empty()) {
    return {};
  }

  const AlignedSignals orderedStates = buildOrderedStatePairs(model0, model1);
  if (!orderedStates.names.empty() &&
      areAllOrderedStatesEquivalent(model0, model1, alignedInputs, orderedStates)) {
    // This fast path is purely structural: same order, same transition shape.
    return orderedStates;
  }

  const auto inputClasses0 = buildInputClassMap(model0, alignedInputs.keys0);
  const auto inputClasses1 = buildInputClassMap(model1, alignedInputs.keys1);

  std::vector<size_t> classes0 = seedClasses(model0);
  std::vector<size_t> classes1 = seedClasses(model1);

  while (true) {
    // Fixed-point refinement: keep splitting classes until another pass no
    // longer learns anything new from the transition structure.
    const std::vector<size_t> refined0 = refineClasses(model0, inputClasses0, classes0);
    const std::vector<size_t> refined1 = refineClasses(model1, inputClasses1, classes1);
    if (refined0 == classes0 && refined1 == classes1) {
      break;
    }
    classes0 = refined0;
    classes1 = refined1;
  }

  const auto fingerprints0 = computeFinalFingerprints(model0, inputClasses0, classes0);
  const auto fingerprints1 = computeFinalFingerprints(model1, inputClasses1, classes1);

  std::map<std::string, std::vector<size_t>> indicesByFingerprint0;
  std::map<std::string, std::vector<size_t>> indicesByFingerprint1;
  for (size_t i = 0; i < fingerprints0.size(); ++i) {
    indicesByFingerprint0[stateFingerprintKey(fingerprints0[i])].push_back(i);
  }
  for (size_t i = 0; i < fingerprints1.size(); ++i) {
    indicesByFingerprint1[stateFingerprintKey(fingerprints1[i])].push_back(i);
  }

  AlignedSignals aligned;
  size_t pairIndex = 0;
  for (const auto& [fingerprint, indices0] : indicesByFingerprint0) {
    const auto it1 = indicesByFingerprint1.find(fingerprint);
    if (it1 == indicesByFingerprint1.end()) {
      continue;
    }
    const size_t matchedCount = std::min(indices0.size(), it1->second.size());
    for (size_t i = 0; i < matchedCount; ++i) {
      aligned.names.push_back("structural_state_" + std::to_string(pairIndex++));
      aligned.keys0.push_back(model0.stateBits[indices0[i]]);
      aligned.keys1.push_back(model1.stateBits[it1->second[i]]);
    }
  }
  return aligned;
}

}  // namespace KEPLER_FORMAL::SEC
