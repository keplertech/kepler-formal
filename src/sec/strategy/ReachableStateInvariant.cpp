// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "strategy/ReachableStateInvariant.h"

#include <algorithm>
#include <cctype>
#include <optional>
#include <unordered_map>

#include "common/BoolExprUtils.h"
#include "common/SecDiag.h"
#include "strategy/StructuralStateInvariant.h"

namespace KEPLER_FORMAL::SEC {

// Overall reachable-state strengthening algorithm:
// 1. Discover any reset inputs and explicit initial-state information.
// 2. Filter structural state pairs down to ones compatible with startup.
// 3. If reset/bootstrap is available, specialize next-state logic under reset.
// 4. Symbolically push known reset-state values forward for a few cycles.
// 5. Keep only the state equalities that survive each bootstrap step.
// 6. Return both the startup correspondence and the anchored equalities/value
//    facts that later proof engines can rely on.

namespace {

bool isConstBoolExpr(BoolExpr* expr, bool value) {
  return expr != nullptr && expr->getOp() == Op::VAR &&
         expr->getId() == static_cast<size_t>(value ? 1 : 0);
}

std::string normalizeSignalBaseName(const std::string& name) {
  std::string base = name;
  const auto bracket = base.find('[');
  if (bracket != std::string::npos) {
    base = base.substr(0, bracket);
  }
  std::transform(base.begin(), base.end(), base.begin(), [](unsigned char ch) {
    return static_cast<char>(std::toupper(ch));
  });
  return base;
}

std::optional<bool> getResetAssertionValue(const std::string& displayName) {
  const std::string normalized = normalizeSignalBaseName(displayName);
  if (normalized == "RESET" || normalized == "RST") {
    return true;
  }
  if (normalized == "RESET_N" || normalized == "RESETN" ||
      normalized == "RST_N" || normalized == "RSTN") {
    return false;
  }
  return std::nullopt;
}

std::unordered_map<size_t, bool> collectResetAssignments(
    const SequentialDesignModel& model) {
  // Reset controls are identified from the aligned user-visible signal names
  // and then converted into the local BoolExpr variable IDs used by the model.
  std::unordered_map<size_t, bool> assignments;
  for (const auto& key : model.environmentInputs) {
    const auto displayIt = model.displayNameByKey.find(key);
    const auto varIt = model.inputVarByKey.find(key);
    if (displayIt == model.displayNameByKey.end() ||
        varIt == model.inputVarByKey.end()) {
      continue;
    }
    const auto assertedValue = getResetAssertionValue(displayIt->second);
    if (!assertedValue.has_value()) {
      continue;
    }
    assignments.emplace(varIt->second, *assertedValue);
  }
  return assignments;
}

AlignedSignals filterStateEqualitiesByInitialValue(
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

AlignedSignals filterStateEqualitiesByInitialCompatibility(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& candidateStates) {
  AlignedSignals compatibleStates;
  for (size_t i = 0; i < candidateStates.names.size(); ++i) {
    const auto initial0 = model0.initialStateValueByKey.find(candidateStates.keys0[i]);
    const auto initial1 = model1.initialStateValueByKey.find(candidateStates.keys1[i]);
    if (initial0 != model0.initialStateValueByKey.end() &&
        initial1 != model1.initialStateValueByKey.end() &&
        initial0->second != initial1->second) {
      continue;
    }

    compatibleStates.names.push_back(candidateStates.names[i]);
    compatibleStates.keys0.push_back(candidateStates.keys0[i]);
    compatibleStates.keys1.push_back(candidateStates.keys1[i]);
  }
  return compatibleStates;
}

size_t defaultResetBootstrapCycles(bool hasResetBootstrap, bool hasCompleteInitialState) {
  return (hasResetBootstrap && !hasCompleteInitialState) ? 3u : 0u;
}

std::optional<bool> evaluateConstantUnderAssignments(
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
      const auto operand =
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);
      if (operand.has_value()) {
        value = !*operand;
      }
      break;
    }
    case Op::AND: {
      const auto lhs =
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);
      if (lhs.has_value() && !*lhs) {
        value = false;
        break;
      }
      const auto rhs =
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);
      if (rhs.has_value() && !*rhs) {
        value = false;
      } else if (lhs.has_value() && rhs.has_value()) {
        value = *lhs && *rhs;
      }
      break;
    }
    case Op::OR: {
      const auto lhs =
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);
      if (lhs.has_value() && *lhs) {
        value = true;
        break;
      }
      const auto rhs =
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);
      if (rhs.has_value() && *rhs) {
        value = true;
      } else if (lhs.has_value() && rhs.has_value()) {
        value = *lhs || *rhs;
      }
      break;
    }
    case Op::XOR: {
      const auto lhs =
          evaluateConstantUnderAssignments(expr->getLeft(), assignments, memo);
      const auto rhs =
          evaluateConstantUnderAssignments(expr->getRight(), assignments, memo);
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

std::unordered_map<SignalKey, bool, SignalKeyHash> deriveResetBootstrapStateValues(
    const SequentialDesignModel& model,
    size_t cycles) {
  // This is a small symbolic reset simulation. We keep only states whose value
  // becomes constant under the asserted-reset environment.
  const auto resetAssignments = collectResetAssignments(model);
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
      const auto value = evaluateConstantUnderAssignments(
          model.nextStateExprByStateKey.at(key), assignments, memo);
      if (value.has_value()) {
        nextKnownStates.emplace(key, *value);
      }
    }
    knownStates = std::move(nextKnownStates);
  }

  return knownStates;
}

AlignedSignals deriveResetBootstrapStateEqualities(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& candidateStates,
    size_t cycles,
    bool secDiagEnabled) {
  // Push candidate state equalities through the reset-specialized next-state
  // logic. A pair survives only if both sides either collapse to the same
  // constant or stay structurally equivalent after each bootstrap step.
  if (cycles == 0 || candidateStates.names.empty()) {
    return filterStateEqualitiesByInitialValue(model0, model1, candidateStates);
  }

  const auto resetAssignments0 = collectResetAssignments(model0);
  const auto resetAssignments1 = collectResetAssignments(model1);
  if (resetAssignments0.empty() || resetAssignments1.empty()) {
    return filterStateEqualitiesByInitialValue(model0, model1, candidateStates);
  }

  auto specializeForReset = [](const SequentialDesignModel& model,
                               const std::unordered_map<size_t, bool>& resetAssignments) {
    std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> specialized;
    specialized.reserve(model.stateBits.size());
    std::unordered_map<BoolExpr*, BoolExpr*> memo;
    for (const auto& key : model.stateBits) {
      specialized.emplace(
          key,
          substituteBoolExprVariables(
              model.nextStateExprByStateKey.at(key), resetAssignments, memo));
    }
    return specialized;
  };

  const auto resetSpecializedNext0 = specializeForReset(model0, resetAssignments0);
  const auto resetSpecializedNext1 = specializeForReset(model1, resetAssignments1);

  AlignedSignals currentEqualities = filterStateEqualitiesByInitialValue(
      model0, model1, candidateStates);
  std::unordered_map<SignalKey, bool, SignalKeyHash> currentKnownValues0 =
      model0.initialStateValueByKey;
  std::unordered_map<SignalKey, bool, SignalKeyHash> currentKnownValues1 =
      model1.initialStateValueByKey;

  for (size_t step = 0; step < cycles; ++step) {
    std::unordered_map<size_t, bool> stateAssignments0;
    std::unordered_map<size_t, bool> stateAssignments1;
    stateAssignments0.reserve(currentKnownValues0.size());
    stateAssignments1.reserve(currentKnownValues1.size());
    for (const auto& [key, value] : currentKnownValues0) {
      if (const auto it = model0.inputVarByKey.find(key); it != model0.inputVarByKey.end()) {
        stateAssignments0.emplace(it->second, value);
      }
    }
    for (const auto& [key, value] : currentKnownValues1) {
      if (const auto it = model1.inputVarByKey.find(key); it != model1.inputVarByKey.end()) {
        stateAssignments1.emplace(it->second, value);
      }
    }

    std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> specializedNext0;
    std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> specializedNext1;
    specializedNext0.reserve(model0.stateBits.size());
    specializedNext1.reserve(model1.stateBits.size());
    std::unordered_map<BoolExpr*, BoolExpr*> stateSubMemo0;
    std::unordered_map<BoolExpr*, BoolExpr*> stateSubMemo1;
    for (const auto& key : model0.stateBits) {
      specializedNext0.emplace(
          key,
          substituteBoolExprVariables(
              resetSpecializedNext0.at(key), stateAssignments0, stateSubMemo0));
    }
    for (const auto& key : model1.stateBits) {
      specializedNext1.emplace(
          key,
          substituteBoolExprVariables(
              resetSpecializedNext1.at(key), stateAssignments1, stateSubMemo1));
    }

    std::unordered_map<SignalKey, bool, SignalKeyHash> nextKnownValues0;
    std::unordered_map<SignalKey, bool, SignalKeyHash> nextKnownValues1;
    nextKnownValues0.reserve(model0.stateBits.size());
    nextKnownValues1.reserve(model1.stateBits.size());
    for (const auto& key : model0.stateBits) {
      if (isConstBoolExpr(specializedNext0.at(key), false)) {
        nextKnownValues0.emplace(key, false);
      } else if (isConstBoolExpr(specializedNext0.at(key), true)) {
        nextKnownValues0.emplace(key, true);
      }
    }
    for (const auto& key : model1.stateBits) {
      if (isConstBoolExpr(specializedNext1.at(key), false)) {
        nextKnownValues1.emplace(key, false);
      } else if (isConstBoolExpr(specializedNext1.at(key), true)) {
        nextKnownValues1.emplace(key, true);
      }
    }

    const auto [abstractMap0, abstractMap1] = buildAbstractTransitionMaps(
        model0, model1, alignedInputs, currentEqualities);
    AlignedSignals nextEqualities;
    for (size_t i = 0; i < candidateStates.names.size(); ++i) {
      const auto& key0 = candidateStates.keys0[i];
      const auto& key1 = candidateStates.keys1[i];

      bool equalAfterStep = false;
      const auto known0 = nextKnownValues0.find(key0);
      const auto known1 = nextKnownValues1.find(key1);
      if (known0 != nextKnownValues0.end() && known1 != nextKnownValues1.end() &&
          known0->second == known1->second) {
        equalAfterStep = true;
      } else {
        equalAfterStep = areEquivalentUnderAbstractMaps(
            specializedNext0.at(key0),
            specializedNext1.at(key1),
            abstractMap0,
            abstractMap1);
      }

      if (!equalAfterStep) {
        continue;
      }

      nextEqualities.names.push_back(candidateStates.names[i]);
      nextEqualities.keys0.push_back(key0);
      nextEqualities.keys1.push_back(key1);
    }

    currentEqualities = std::move(nextEqualities);
    currentKnownValues0 = std::move(nextKnownValues0);
    currentKnownValues1 = std::move(nextKnownValues1);
    if (secDiagEnabled) {
      emitSecDiag(
          "SEC diag: bootstrap step ",
          step + 1,
          " equalities=",
          currentEqualities.names.size(),
          " known0=",
          currentKnownValues0.size(),
          " known1=",
          currentKnownValues1.size());
    }
  }

  return currentEqualities;
}

bool hasCompleteInitialState(const SequentialDesignModel& model0,
                             const SequentialDesignModel& model1) {
  return model0.initialStateValueByKey.size() == model0.stateBits.size() &&
         model1.initialStateValueByKey.size() == model1.stateBits.size();
}

bool hasExplicitInitialState(const SequentialDesignModel& model0,
                             const SequentialDesignModel& model1) {
  return !model0.initialStateValueByKey.empty() ||
         !model1.initialStateValueByKey.empty();
}

}  // namespace

ReachableStateInvariant buildReachableStateInvariant(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& inductiveStateEqualities,
    bool secDiagEnabled) {
  ReachableStateInvariant invariant;
  // First decide which startup model we have: explicit init, reset bootstrap,
  // both, or neither. That determines how strong the frame-0 correspondence
  // may safely be.
  const bool hasResetBootstrap = !collectResetAssignments(model0).empty() &&
      !collectResetAssignments(model1).empty();

  invariant.bootstrapCycles = defaultResetBootstrapCycles(
      hasResetBootstrap, hasCompleteInitialState(model0, model1));
  invariant.initialStateCorrespondence = filterStateEqualitiesByInitialCompatibility(
      model0, model1, inductiveStateEqualities);

  if (hasResetBootstrap) {
    // Prefer the strongest already-justified correspondence. If none exists at
    // frame 0, derive one by walking the reset window forward.
    if (!invariant.initialStateCorrespondence.names.empty()) {
      invariant.anchoredStateEqualities = invariant.initialStateCorrespondence;
    } else {
      invariant.anchoredStateEqualities = deriveResetBootstrapStateEqualities(
          model0,
          model1,
          alignedInputs,
          inductiveStateEqualities,
          invariant.bootstrapCycles,
          secDiagEnabled);
    }

    invariant.bootstrapValues0 =
        deriveResetBootstrapStateValues(model0, invariant.bootstrapCycles);
    invariant.bootstrapValues1 =
        deriveResetBootstrapStateValues(model1, invariant.bootstrapCycles);
  } else if (hasExplicitInitialState(model0, model1)) {
    // Without reset, we can only anchor the state pairs whose explicit init
    // values agree on both sides.
    invariant.anchoredStateEqualities = filterStateEqualitiesByInitialValue(
        model0, model1, inductiveStateEqualities);
  }

  return invariant;
}

namespace detail {

std::unordered_map<SignalKey, bool, SignalKeyHash>
deriveResetBootstrapStateValuesForTest(
    const SequentialDesignModel& model,
    size_t cycles) {
  return deriveResetBootstrapStateValues(model, cycles);
}

AlignedSignals filterStateEqualitiesByInitialValueForTest(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& candidateStates) {
  return filterStateEqualitiesByInitialValue(model0, model1, candidateStates);
}

}  // namespace detail

}  // namespace KEPLER_FORMAL::SEC
