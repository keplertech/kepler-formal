// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "clocks/SecClockModel.h"

#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace KEPLER_FORMAL::SEC {

namespace {

bool isBoolConst(BoolExpr* expr, bool value) {
  return expr != nullptr && expr->getOp() == Op::VAR &&
         expr->getId() == static_cast<size_t>(value ? 1 : 0);
}

BoolExpr* simplifyWhenChanged(BoolExpr* original, BoolExpr* changed) {
  return changed == original ? original : BoolExpr::simplify(changed);
}

BoolExpr* substituteBoolExprAssignments(
    BoolExpr* root,
    const std::unordered_map<size_t, bool>& assignments,
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
      const size_t id = root->getId();
      if (id < 2) {
        substituted = root;
        break;
      }
      const auto assignmentIt = assignments.find(id);
      substituted = assignmentIt == assignments.end()
                        ? root
                        : (assignmentIt->second ? BoolExpr::createTrue()
                                                : BoolExpr::createFalse());
      break;
    }
    case Op::NOT: {
      BoolExpr* left =
          substituteBoolExprAssignments(root->getLeft(), assignments, memo);
      substituted = left == root->getLeft() ? root : BoolExpr::Not(left);
      break;
    }
    case Op::AND: {
      BoolExpr* left =
          substituteBoolExprAssignments(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprAssignments(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::And(left, right);
      break;
    }
    case Op::OR: {
      BoolExpr* left =
          substituteBoolExprAssignments(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprAssignments(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Or(left, right);
      break;
    }
    case Op::XOR: {
      BoolExpr* left =
          substituteBoolExprAssignments(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprAssignments(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Xor(left, right);
      break;
    }
    case Op::NONE:
    default:
      throw std::runtime_error("Unsupported BoolExpr operator in clock analysis");
  }

  memo.emplace(root, substituted);
  return substituted;
}

BoolExpr* specializeClockCarrier(BoolExpr* expr, size_t carrierVarID, bool value) {
  std::unordered_map<size_t, bool> assignments{{carrierVarID, value}};
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return simplifyWhenChanged(
      expr, substituteBoolExprAssignments(expr, assignments, memo));
}

BoolExpr* makeRisingEnable(BoolExpr* whenClockLow, BoolExpr* whenClockHigh) {
  return BoolExpr::simplify(
      BoolExpr::And(BoolExpr::Not(whenClockLow), whenClockHigh));
}

BoolExpr* makeFallingEnable(BoolExpr* whenClockLow, BoolExpr* whenClockHigh) {
  return BoolExpr::simplify(
      BoolExpr::And(whenClockLow, BoolExpr::Not(whenClockHigh)));
}

std::optional<size_t> singleCarrierInSupport(
    BoolExpr* expr,
    const std::unordered_map<size_t, ClockEvent>& carrierEvents) {
  std::optional<size_t> carrierVarID;
  for (const auto varID : expr->getSupportVars()) {
    if (carrierEvents.find(varID) == carrierEvents.end()) {
      continue;
    }
    if (carrierVarID.has_value() && *carrierVarID != varID) {
      return std::nullopt;
    }
    carrierVarID = varID;
  }
  return carrierVarID;
}

BoolExpr* normalizeEventEnable(BoolExpr* enable) {
  return isBoolConst(enable, true) ? nullptr : enable;
}

}  // namespace

ClockPhase invertClockPhase(ClockPhase phase) {
  return phase == ClockPhase::Pos ? ClockPhase::Neg : ClockPhase::Pos;
}

const char* clockPhaseName(ClockPhase phase) {
  return phase == ClockPhase::Pos ? "posedge" : "negedge";
}

BoolExpr* clockEventEnableOrTrue(const ClockEvent& event) {
  return event.enable == nullptr ? BoolExpr::createTrue() : event.enable;
}

bool clockEventIsUngated(const ClockEvent& event) {
  return event.enable == nullptr || isBoolConst(event.enable, true);
}

std::optional<ClockEvent> classifyClockEventExpression(
    BoolExpr* expr,
    const std::unordered_map<size_t, ClockEvent>& carrierEvents) {
  if (expr == nullptr || carrierEvents.empty()) {
    return std::nullopt;
  }

  const auto carrierVarID = singleCarrierInSupport(expr, carrierEvents);
  if (!carrierVarID.has_value()) {
    return std::nullopt;
  }

  const auto carrierIt = carrierEvents.find(*carrierVarID);
  if (carrierIt == carrierEvents.end()) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  BoolExpr* whenClockLow = specializeClockCarrier(expr, *carrierVarID, false);
  BoolExpr* whenClockHigh = specializeClockCarrier(expr, *carrierVarID, true);
  BoolExpr* risingEnable = makeRisingEnable(whenClockLow, whenClockHigh);
  BoolExpr* fallingEnable = makeFallingEnable(whenClockLow, whenClockHigh);
  const bool hasRising = !isBoolConst(risingEnable, false);
  const bool hasFalling = !isBoolConst(fallingEnable, false);

  // A dynamically selected phase, such as clk XOR sel, can create either edge.
  // SEC does not model CDC/glitch assumptions, so fail classification closed.
  if (hasRising == hasFalling) {
    return std::nullopt;
  }

  ClockEvent event = carrierIt->second;
  if (hasRising) {
    event.enable = normalizeEventEnable(risingEnable);
  } else {
    event.phase = invertClockPhase(event.phase);
    event.enable = normalizeEventEnable(fallingEnable);
  }
  return event;
}

BoolExpr* substituteBoolExprVariableExpressionsImpl(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements,
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
      const auto replacementIt = replacements.find(root->getId());
      substituted = replacementIt == replacements.end() ? root : replacementIt->second;
      break;
    }
    case Op::NOT: {
      BoolExpr* left =
          substituteBoolExprVariableExpressionsImpl(root->getLeft(), replacements, memo);
      substituted = left == root->getLeft() ? root : BoolExpr::Not(left);
      break;
    }
    case Op::AND: {
      BoolExpr* left =
          substituteBoolExprVariableExpressionsImpl(root->getLeft(), replacements, memo);
      BoolExpr* right =
          substituteBoolExprVariableExpressionsImpl(root->getRight(), replacements, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::And(left, right);
      break;
    }
    case Op::OR: {
      BoolExpr* left =
          substituteBoolExprVariableExpressionsImpl(root->getLeft(), replacements, memo);
      BoolExpr* right =
          substituteBoolExprVariableExpressionsImpl(root->getRight(), replacements, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Or(left, right);
      break;
    }
    case Op::XOR: {
      BoolExpr* left =
          substituteBoolExprVariableExpressionsImpl(root->getLeft(), replacements, memo);
      BoolExpr* right =
          substituteBoolExprVariableExpressionsImpl(root->getRight(), replacements, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Xor(left, right);
      break;
    }
    case Op::NONE:
    default:
      throw std::runtime_error("Unsupported BoolExpr operator in expression substitution");
  }

  memo.emplace(root, substituted);
  return substituted;
}

BoolExpr* substituteBoolExprVariableExpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements) {
  if (root == nullptr || replacements.empty()) {
    return root;
  }
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return simplifyWhenChanged(
      root, substituteBoolExprVariableExpressionsImpl(root, replacements, memo));
}

}  // namespace KEPLER_FORMAL::SEC
