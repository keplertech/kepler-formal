// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "common/BoolExprUtils.h"

#include <optional>
#include <stdexcept>
#include <vector>

#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC {

namespace {

bool isBoolConst(BoolExpr* expr, bool value) {
  return expr != nullptr && expr->getOp() == Op::VAR &&
         expr->getId() == static_cast<size_t>(value ? 1 : 0);
}

SATSolverWrapper::SolveStatus solveBoolFormulaStatus(
    BoolExpr* formula,
    KEPLER_FORMAL::Config::SolverType solverType,
    std::optional<unsigned> conflictLimit) {
  if (formula == nullptr || isBoolConst(formula, false)) {
    return SATSolverWrapper::SolveStatus::Unsat;
  }
  if (isBoolConst(formula, true)) {
    return SATSolverWrapper::SolveStatus::Sat;
  }

  SATSolverWrapper solver(solverType);
  const auto support = formula->getSupportVars();
  solver.configureForSecLocalBooleanCheck(support.size());
  std::unordered_map<size_t, int> leafLits;
  leafLits.reserve(support.size());
  for (const auto symbol : support) {
    if (symbol < 2) {
      continue;
    }
    leafLits.emplace(symbol, solver.newVar() + 2);
  }

  FrameFormulaEncoder encoder(solver, std::move(leafLits));
  solver.addClause({encoder.encode(formula)});
  if (conflictLimit.has_value()) {
    if (solverType == KEPLER_FORMAL::Config::SolverType::KISSAT) { // LCOV_EXCL_LINE
      return solver.solveWithKissatResourceLimits( // LCOV_EXCL_LINE
          *conflictLimit, *conflictLimit); // LCOV_EXCL_LINE
    }
    return solver.solveWithAssumptionsStatus( // LCOV_EXCL_LINE
        {}, // LCOV_EXCL_LINE
        *conflictLimit, // LCOV_EXCL_LINE
        /*propagationLimit=*/-1);
  }
  return solver.solveStatus();
}

}  // namespace

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

  struct StackFrame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };

  // Large gate-level SEC formulas are often long, hash-consed DAGs.  A
  // recursive remap is semantically simple, but profiling on BlackParrot-sized
  // designs showed it consuming substantial time and stack depth before the
  // SAT engine even sees a query.  This iterative post-order traversal keeps
  // the same memoized remapping while making the pass linear in the reached
  // DAG and robust for very deep cones.
  std::vector<StackFrame> stack;
  stack.push_back({root, false});
  while (!stack.empty()) {
    const StackFrame current = stack.back();
    stack.pop_back();
    BoolExpr* node = current.expr;
    // LCOV_EXCL_START
    if (node == nullptr || memo.find(node) != memo.end()) {
      continue;
    }
    // LCOV_EXCL_STOP

    if (node->getOp() == Op::VAR) {
      const size_t id = node->getId();
      if (id < 2) {
        memo.emplace(node, node);
      } else {
        auto it = varMap.find(id);
        if (it == varMap.end()) {
          throw std::runtime_error("Missing BoolExpr remap for variable " +
                                   std::to_string(id));
        }
        memo.emplace(node, it->second == id ? node : BoolExpr::Var(it->second));
      }
      continue;
    }

    if (node->getOp() == Op::NONE) {
      throw std::runtime_error("Unsupported BoolExpr operator in remap");
    }

    if (!current.visited) {
      stack.push_back({node, true});
      if (node->getRight() != nullptr &&
          memo.find(node->getRight()) == memo.end()) {
        stack.push_back({node->getRight(), false});
      }
      if (node->getLeft() != nullptr &&
          memo.find(node->getLeft()) == memo.end()) {
        stack.push_back({node->getLeft(), false});
      }
      continue;
    }

    BoolExpr* remapped = nullptr;
    switch (node->getOp()) {
      case Op::NOT: {
        BoolExpr* left = memo.at(node->getLeft());
        // Stable-variable remapping is often identity on most of a large SEC
        // cone. Preserve unchanged sub-DAGs instead of rebuilding them through
        // the global BoolExpr cache.
        remapped = left == node->getLeft() ? node : BoolExpr::Not(left);
        break;
      }
      case Op::AND: {
        BoolExpr* left = memo.at(node->getLeft());
        BoolExpr* right = memo.at(node->getRight());
        remapped = left == node->getLeft() && right == node->getRight()
                       ? node
                       : BoolExpr::And(left, right);
        break;
      }
      case Op::OR: {
        BoolExpr* left = memo.at(node->getLeft());
        BoolExpr* right = memo.at(node->getRight());
        remapped = left == node->getLeft() && right == node->getRight()
                       ? node
                       : BoolExpr::Or(left, right);
        break;
      }
      case Op::XOR: {
        BoolExpr* left = memo.at(node->getLeft());
        BoolExpr* right = memo.at(node->getRight());
        remapped = left == node->getLeft() && right == node->getRight()
                       ? node
                       : BoolExpr::Xor(left, right);
        break;
      }
      // LCOV_EXCL_START
      case Op::VAR:
      // LCOV_DISABLED_START
      case Op::NONE:
      // LCOV_DISABLED_STOP
      default:
        // LCOV_DISABLED_START
        throw std::runtime_error("Unsupported BoolExpr operator in remap");
        // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
    memo.emplace(node, remapped);
  }

  return memo.at(root);
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
        substituted = root;
        break;
      }
      auto it = assignments.find(id);
      substituted =
          it == assignments.end() ? root
                                  : (it->second ? BoolExpr::createTrue()
                                                : BoolExpr::createFalse());
      break;
    }
    case Op::NOT: {
      BoolExpr* left =
          substituteBoolExprVariables(root->getLeft(), assignments, memo);
      // Preserve no-op subtrees. Reset/bootstrap specialization touches only a
      // small frontier on large ASIC cones; rebuilding untouched cones made
      // BlackParrot spend minutes constructing equivalent BoolExpr nodes.
      substituted = left == root->getLeft() ? root : BoolExpr::Not(left);
      break;
    }
    case Op::AND: {
      BoolExpr* left =
          substituteBoolExprVariables(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprVariables(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::And(left, right);
      break;
    }
    case Op::OR: {
      BoolExpr* left =
          substituteBoolExprVariables(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprVariables(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Or(left, right);
      break;
    }
    case Op::XOR: {
      BoolExpr* left =
          substituteBoolExprVariables(root->getLeft(), assignments, memo);
      BoolExpr* right =
          substituteBoolExprVariables(root->getRight(), assignments, memo);
      substituted = left == root->getLeft() && right == root->getRight()
                        ? root
                        : BoolExpr::Xor(left, right);
      break;
    }
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

bool isBoolFormulaSatisfiable(
    BoolExpr* formula,
    KEPLER_FORMAL::Config::SolverType solverType) {
  return solveBoolFormulaStatus(formula, solverType, std::nullopt) ==
         SATSolverWrapper::SolveStatus::Sat;
}

bool boolFormulaImplies(
    BoolExpr* assumptions,
    BoolExpr* conclusion,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (conclusion == nullptr) {
    return false;  // LCOV_EXCL_LINE
  }
  if (isBoolConst(conclusion, true) || isBoolConst(assumptions, false)) {
    return true;
  }

  // Ask SAT for an assignment that satisfies the assumptions and violates the
  // conclusion. If no such assignment exists, the implication is a theorem.
  return !isBoolFormulaSatisfiable(
      BoolExpr::And(assumptions, BoolExpr::Not(conclusion)), solverType);
}

std::optional<bool> boolFormulaImpliesWithConflictLimit(
    BoolExpr* assumptions,
    BoolExpr* conclusion,
    KEPLER_FORMAL::Config::SolverType solverType,
    unsigned conflictLimit) {
  if (conclusion == nullptr) {
    // LCOV_EXCL_START
    return false;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (isBoolConst(conclusion, true) || isBoolConst(assumptions, false)) {
    return true; // LCOV_EXCL_LINE
  }

  const auto status = solveBoolFormulaStatus(
      BoolExpr::And(assumptions, BoolExpr::Not(conclusion)),
      solverType,
      conflictLimit);
  if (status == SATSolverWrapper::SolveStatus::Unknown) {
    // LCOV_EXCL_START
    return std::nullopt;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  return status == SATSolverWrapper::SolveStatus::Unsat;
}

}  // namespace KEPLER_FORMAL::SEC
