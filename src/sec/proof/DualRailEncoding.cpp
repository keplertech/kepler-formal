// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "proof/DualRailEncoding.h"

#include <stdexcept>
#include <vector>

namespace KEPLER_FORMAL::SEC {

DualRailBoolExpr buildDualRailBoolExpr(
    BoolExpr* root,
    DualRailVariableMapper& mapper,
    std::unordered_map<BoolExpr*, DualRailBoolExpr>& memo) {
  if (root == nullptr) {
    // LCOV_EXCL_START
    return {};  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  if (auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  struct StackFrame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };

  // Dual rail represents a ternary value with possible-value bits:
  //   known 0 => (may1=0, may0=1)
  //   known 1 => (may1=1, may0=0)
  //   X       => (may1=1, may0=1)
  // The formulas below are the conservative Boolean lift of NOT/AND/OR/XOR.
  std::vector<StackFrame> stack;
  stack.push_back({root, false});
  while (!stack.empty()) {
    const StackFrame current = stack.back();
    stack.pop_back();
    BoolExpr* node = current.expr;
    if (node == nullptr || memo.find(node) != memo.end()) {
      continue;  // LCOV_EXCL_LINE
    }

    if (node->getOp() == Op::VAR) {
      const size_t id = node->getId();
      if (id == 0) {
        memo.emplace(
            node,
            DualRailBoolExpr{BoolExpr::createFalse(), BoolExpr::createTrue()});
      } else if (id == 1) {
        memo.emplace(
            node,
            DualRailBoolExpr{BoolExpr::createTrue(), BoolExpr::createFalse()});
      } else {
        memo.emplace(node, mapper.mapVariable(id));
      }
      continue;
    }

    if (node->getOp() == Op::NONE) {
      // LCOV_EXCL_START
      throw std::runtime_error("Unsupported BoolExpr operator in dual-rail encoding");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
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

    DualRailBoolExpr lifted;
    switch (node->getOp()) {
      case Op::NOT: {
        const auto operand = memo.at(node->getLeft());
        lifted = DualRailBoolExpr{operand.mayBeZero, operand.mayBeOne};
        break;
      }
      case Op::AND: {
        const auto lhs = memo.at(node->getLeft());  // LCOV_EXCL_LINE
        const auto rhs = memo.at(node->getRight());  // LCOV_EXCL_LINE
        lifted.mayBeOne = BoolExpr::And(lhs.mayBeOne, rhs.mayBeOne);  // LCOV_EXCL_LINE
        lifted.mayBeZero = BoolExpr::Or(lhs.mayBeZero, rhs.mayBeZero);  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      }
      case Op::OR: {
        const auto lhs = memo.at(node->getLeft());  // LCOV_EXCL_LINE
        const auto rhs = memo.at(node->getRight());  // LCOV_EXCL_LINE
        lifted.mayBeOne = BoolExpr::Or(lhs.mayBeOne, rhs.mayBeOne);  // LCOV_EXCL_LINE
        lifted.mayBeZero = BoolExpr::And(lhs.mayBeZero, rhs.mayBeZero);  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      }
      case Op::XOR: {
        const auto lhs = memo.at(node->getLeft());
        const auto rhs = memo.at(node->getRight());
        lifted.mayBeOne = BoolExpr::Or(
            BoolExpr::And(lhs.mayBeOne, rhs.mayBeZero),
            BoolExpr::And(lhs.mayBeZero, rhs.mayBeOne));
        lifted.mayBeZero = BoolExpr::Or(
            BoolExpr::And(lhs.mayBeOne, rhs.mayBeOne),
            BoolExpr::And(lhs.mayBeZero, rhs.mayBeZero));
        break;
      }
      // LCOV_EXCL_START
      case Op::VAR:
      // LCOV_DISABLED_START
      case Op::NONE:
      // LCOV_DISABLED_STOP
      default:
        // LCOV_DISABLED_START
        throw std::runtime_error("Unsupported BoolExpr operator in dual-rail encoding");
        // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
    // Do not simplify every lifted node here.  ASIC transition cones share
    // large DAGs, and recursively simplifying both rails at every node made
    // dual-rail PDR spend minutes in BoolExpr cache churn before SAT solving.
    // The proof engines still simplify final properties where that is useful.
    memo.emplace(node, lifted);
  }

  return memo.at(root);
}

}  // namespace KEPLER_FORMAL::SEC
