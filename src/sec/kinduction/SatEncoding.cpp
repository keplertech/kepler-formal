// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "kinduction/SatEncoding.h"

#include <stack>
#include <stdexcept>

namespace KEPLER_FORMAL::SEC {

namespace {

int newSolverLiteral(SATSolverWrapper& solver) {
  // BoolExpr reserves 0/1 for false/true, so fresh SAT literals start above
  // those special ids.
  return solver.newVar() + 2;
}

}  // namespace

FrameVariableStore::FrameVariableStore(SATSolverWrapper& solver,
                                       const std::vector<size_t>& symbols,
                                       size_t numFrames) {
  // Every symbolic SEC variable gets one SAT literal per time frame.
  for (const auto symbol : symbols) {
    auto& frameLits = symbolFrameLits_[symbol];
    frameLits.reserve(numFrames);
    for (size_t frame = 0; frame < numFrames; ++frame) {
      frameLits.push_back(newSolverLiteral(solver));
    }
  }
}

int FrameVariableStore::getLiteral(size_t symbol, size_t frame) const {
  auto it = symbolFrameLits_.find(symbol);
  if (it == symbolFrameLits_.end() || frame >= it->second.size()) {
    throw std::runtime_error("Missing frame variable for symbol " +
                             std::to_string(symbol));
  }
  return it->second[frame];
}

std::unordered_map<size_t, int> FrameVariableStore::makeLeafLits(
    size_t frame) const {
  std::unordered_map<size_t, int> leafLits;
  leafLits.reserve(symbolFrameLits_.size());
  for (const auto& [symbol, frameLits] : symbolFrameLits_) {
    if (frame >= frameLits.size()) {
      throw std::runtime_error("Frame index is out of range");
    }
    leafLits.emplace(symbol, frameLits[frame]);
  }
  return leafLits;
}

FrameFormulaEncoder::FrameFormulaEncoder(
    SATSolverWrapper& solver,
    std::unordered_map<size_t, int> leafLits)
    : solver_(solver), leafLits_(std::move(leafLits)) {}

int FrameFormulaEncoder::getConstLit(bool value) {
  auto& cache = value ? trueLit_ : falseLit_;
  if (cache.has_value()) {
    return *cache;
  }
  int lit = newSolverLiteral(solver_);
  solver_.addClause({value ? lit : -lit});
  cache = lit;
  return lit;
}

int FrameFormulaEncoder::encode(BoolExpr* expr) {
  if (expr == nullptr) {
    throw std::invalid_argument("FrameFormulaEncoder::encode: null expr");
  }

  struct StackFrame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };

  // Encode iteratively so large BoolExpr DAGs do not rely on recursion depth.
  std::stack<StackFrame> stack;
  stack.push({expr, false});

  while (!stack.empty()) {
    auto current = stack.top();
    stack.pop();
    BoolExpr* node = current.expr;

    if (nodeToLit_.find(node) != nodeToLit_.end()) {
      continue;
    }

    if (node->getOp() == Op::VAR) {
      if (node->getId() == 0) {
        nodeToLit_.emplace(node, getConstLit(false));
      } else if (node->getId() == 1) {
        nodeToLit_.emplace(node, getConstLit(true));
      } else {
        auto it = leafLits_.find(node->getId());
        if (it == leafLits_.end()) {
          throw std::runtime_error("Missing leaf literal for symbol " +
                                   std::to_string(node->getId()));
        }
        nodeToLit_.emplace(node, it->second);
      }
      continue;
    }

    if (!current.visited) {
      stack.push({node, true});
      if (node->getRight()) {
        stack.push({node->getRight(), false});
      }
      if (node->getLeft()) {
        stack.push({node->getLeft(), false});
      }
      continue;
    }

    const int leftLit = node->getLeft() ? nodeToLit_.at(node->getLeft()) : 0;
    const int rightLit = node->getRight() ? nodeToLit_.at(node->getRight()) : 0;
    const int lit = newSolverLiteral(solver_);

    // Standard Tseitin clauses for the BoolExpr node at this frame.
    switch (node->getOp()) {
      case Op::NOT:
        solver_.addClause({-lit, -leftLit});
        solver_.addClause({lit, leftLit});
        break;
      case Op::AND:
        solver_.addClause({-lit, leftLit});
        solver_.addClause({-lit, rightLit});
        solver_.addClause({lit, -leftLit, -rightLit});
        break;
      case Op::OR:
        solver_.addClause({-leftLit, lit});
        solver_.addClause({-rightLit, lit});
        solver_.addClause({-lit, leftLit, rightLit});
        break;
      case Op::XOR:
        solver_.addClause({-lit, -leftLit, -rightLit});
        solver_.addClause({-lit, leftLit, rightLit});
        solver_.addClause({lit, -leftLit, rightLit});
        solver_.addClause({lit, leftLit, -rightLit});
        break;
      case Op::VAR:
      case Op::NONE:
      default:
        throw std::runtime_error("Unsupported BoolExpr operator in frame encoder");
    }

    nodeToLit_.emplace(node, lit);
  }

  return nodeToLit_.at(expr);
}

void addLiteralEquivalence(SATSolverWrapper& solver, int lhs, int rhs) {
  solver.addClause({-lhs, rhs});
  solver.addClause({lhs, -rhs});
}

int createXorLiteral(SATSolverWrapper& solver, int lhs, int rhs) {
  const int lit = newSolverLiteral(solver);
  solver.addClause({-lit, -lhs, -rhs});
  solver.addClause({-lit, lhs, rhs});
  solver.addClause({lit, -lhs, rhs});
  solver.addClause({lit, lhs, -rhs});
  return lit;
}

void addSimplePathConstraint(SATSolverWrapper& solver,
                             const FrameVariableStore& variables,
                             const std::vector<size_t>& stateSymbols,
                             size_t numFrames) {
  addSimplePathConstraint(solver, variables, stateSymbols, 0, numFrames);
}

void addSimplePathConstraint(SATSolverWrapper& solver,
                             const FrameVariableStore& variables,
                             const std::vector<size_t>& stateSymbols,
                             size_t firstFrame,
                             size_t numFrames) {
  if (stateSymbols.empty()) {
    return;
  }

  // Slide-48 refinement: every pair of frames must differ in at least one
  // state bit, which rules out cyclic paths in the induction step.
  const size_t lastFrame = firstFrame + numFrames;
  for (size_t i = firstFrame; i < lastFrame; ++i) {
    for (size_t j = i + 1; j < lastFrame; ++j) {
      std::vector<int> diffClause;
      diffClause.reserve(stateSymbols.size());
      for (const auto symbol : stateSymbols) {
        diffClause.push_back(createXorLiteral(
            solver,
            variables.getLiteral(symbol, i),
            variables.getLiteral(symbol, j)));
      }
      solver.addClause(diffClause);
    }
  }
}

}  // namespace KEPLER_FORMAL::SEC
