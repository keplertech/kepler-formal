// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "common/BoolExprUtils.h"

#include <stdexcept>

namespace KEPLER_FORMAL::SEC {

BoolExpr* remapBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, size_t>& varMap,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    return nullptr;
  }

  if (auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  BoolExpr* remapped = nullptr;
  switch (root->getOp()) {
    case Op::VAR: {
      const size_t id = root->getId();
      if (id < 2) {
        remapped = BoolExpr::Var(id);
        break;
      }
      auto it = varMap.find(id);
      if (it == varMap.end()) {
        throw std::runtime_error("Missing BoolExpr remap for variable " +
                                 std::to_string(id));
      }
      remapped = BoolExpr::Var(it->second);
      break;
    }
    case Op::NOT:
      remapped = BoolExpr::Not(
          remapBoolExprVariables(root->getLeft(), varMap, memo));
      break;
    case Op::AND:
      remapped = BoolExpr::And(
          remapBoolExprVariables(root->getLeft(), varMap, memo),
          remapBoolExprVariables(root->getRight(), varMap, memo));
      break;
    case Op::OR:
      remapped = BoolExpr::Or(
          remapBoolExprVariables(root->getLeft(), varMap, memo),
          remapBoolExprVariables(root->getRight(), varMap, memo));
      break;
    case Op::XOR:
      remapped = BoolExpr::Xor(
          remapBoolExprVariables(root->getLeft(), varMap, memo),
          remapBoolExprVariables(root->getRight(), varMap, memo));
      break;
    case Op::NONE:
    default:
      throw std::runtime_error("Unsupported BoolExpr operator in remap");
  }

  memo.emplace(root, remapped);
  return remapped;
}

BoolExpr* remapBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, size_t>& varMap) {
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return remapBoolExprVariables(root, varMap, memo);
}

BoolExpr* substituteBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, bool>& assignments,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    return nullptr;
  }

  if (auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  BoolExpr* substituted = nullptr;
  switch (root->getOp()) {
    case Op::VAR: {
      const size_t id = root->getId();
      if (id < 2) {
        substituted = BoolExpr::Var(id);
        break;
      }
      auto it = assignments.find(id);
      substituted =
          it == assignments.end() ? BoolExpr::Var(id)
                                  : (it->second ? BoolExpr::createTrue()
                                                : BoolExpr::createFalse());
      break;
    }
    case Op::NOT:
      substituted = BoolExpr::Not(
          substituteBoolExprVariables(root->getLeft(), assignments, memo));
      break;
    case Op::AND:
      substituted = BoolExpr::And(
          substituteBoolExprVariables(root->getLeft(), assignments, memo),
          substituteBoolExprVariables(root->getRight(), assignments, memo));
      break;
    case Op::OR:
      substituted = BoolExpr::Or(
          substituteBoolExprVariables(root->getLeft(), assignments, memo),
          substituteBoolExprVariables(root->getRight(), assignments, memo));
      break;
    case Op::XOR:
      substituted = BoolExpr::Xor(
          substituteBoolExprVariables(root->getLeft(), assignments, memo),
          substituteBoolExprVariables(root->getRight(), assignments, memo));
      break;
    case Op::NONE:
    default:
      throw std::runtime_error("Unsupported BoolExpr operator in substitution");
  }

  memo.emplace(root, substituted);
  return substituted;
}

BoolExpr* substituteBoolExprVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, bool>& assignments) {
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return substituteBoolExprVariables(root, assignments, memo);
}

BoolExpr* substituteBoolExprSubexpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    return nullptr;
  }

  if (auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  BoolExpr* substituted = nullptr;
  switch (root->getOp()) {
    case Op::VAR: {
      const size_t id = root->getId();
      if (id < 2) {
        substituted = BoolExpr::Var(id);
        break;
      }
      auto it = replacements.find(id);
      substituted = it == replacements.end() ? BoolExpr::Var(id) : it->second;
      break;
    }
    case Op::NOT:
      substituted = BoolExpr::Not(
          substituteBoolExprSubexpressions(root->getLeft(), replacements, memo));
      break;
    case Op::AND:
      substituted = BoolExpr::And(
          substituteBoolExprSubexpressions(root->getLeft(), replacements, memo),
          substituteBoolExprSubexpressions(root->getRight(), replacements, memo));
      break;
    case Op::OR:
      substituted = BoolExpr::Or(
          substituteBoolExprSubexpressions(root->getLeft(), replacements, memo),
          substituteBoolExprSubexpressions(root->getRight(), replacements, memo));
      break;
    case Op::XOR:
      substituted = BoolExpr::Xor(
          substituteBoolExprSubexpressions(root->getLeft(), replacements, memo),
          substituteBoolExprSubexpressions(root->getRight(), replacements, memo));
      break;
    case Op::NONE:
    default:
      throw std::runtime_error(
          "Unsupported BoolExpr operator in subexpression substitution");
  }

  memo.emplace(root, substituted);
  return substituted;
}

BoolExpr* substituteBoolExprSubexpressions(
    BoolExpr* root,
    const std::unordered_map<size_t, BoolExpr*>& replacements) {
  std::unordered_map<BoolExpr*, BoolExpr*> memo;
  return substituteBoolExprSubexpressions(root, replacements, memo);
}

}  // namespace KEPLER_FORMAL::SEC
