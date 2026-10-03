// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "proof/InternalRelations.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "common/SecDiag.h"
#include "model/SequentialDesignModel.h"
#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

namespace {

// One bit per simulation pattern.
using Words = std::unordered_map<size_t, uint64_t>;

uint64_t evaluateWords(BoolExpr* root, Words& leaves, std::mt19937_64& random,
                       std::unordered_map<BoolExpr*, uint64_t>& values) {
  std::vector<BoolExpr*> stack{root};
  while (!stack.empty()) {
    BoolExpr* node = stack.back();
    if (values.contains(node)) {
      stack.pop_back();
      continue;
    }
    if (node->getOp() == Op::VAR) {
      // Constants are Var(0) and Var(1). A leaf without a value is an input.
      const size_t id = node->getId();
      values[node] = id == 0 ? 0 : id == 1 ? ~uint64_t{0}
                                           : leaves.try_emplace(id, random()).first->second;
      stack.pop_back();
      continue;
    }
    BoolExpr* left = node->getLeft();
    BoolExpr* right = node->getOp() == Op::NOT ? nullptr : node->getRight();
    const bool leftReady = values.contains(left);
    const bool rightReady = right == nullptr || values.contains(right);
    if (!leftReady || !rightReady) {
      if (!leftReady) stack.push_back(left);
      if (!rightReady) stack.push_back(right);
      continue;
    }
    const uint64_t a = values.at(left);
    const uint64_t b = right == nullptr ? 0 : values.at(right);
    values[node] = node->getOp() == Op::NOT ? ~a : node->getOp() == Op::AND ? a & b
                 : node->getOp() == Op::OR ? a | b : a ^ b;
    stack.pop_back();
  }
  return values.at(root);
}

uint64_t relationWord(const InternalRelationCandidate& candidate, const Words& state,
                      bool allowXEquality) {
  uint64_t holds = ~uint64_t{0};
  for (const auto& [lhs, rhs] : candidate.equalities) {
    holds &= ~(state.at(lhs) ^ state.at(rhs));
  }
  if (!allowXEquality) {
    for (const auto& value : candidate.values) {
      holds &= state.at(value.mayBeOne) ^ state.at(value.mayBeZero);
    }
  }
  return holds;
}

// Rewrites a transition under the hypotheses: a register merged into another
// reads that register's leaf.  The factory is hash-consed, so two transitions
// that are the same formula after the rewrite come back as one node.
//
// A gate with three or more inputs is a chain of two-input nodes, and swapping
// its pins, as physical optimization does, re-brackets that chain.  Short
// AND/OR chains are therefore rebuilt from their sorted operands, which makes
// the result independent of the pin order.  Longer chains keep their shape.
BoolExpr* canonicalUnderMerge(BoolExpr* root,
                              const std::unordered_map<size_t, size_t>& merged,
                              std::unordered_map<BoolExpr*, BoolExpr*>& canonical) {
  constexpr size_t kChainOperandLimit = 16;
  std::vector<BoolExpr*> chain;
  std::vector<BoolExpr*> operands;
  std::vector<BoolExpr*> stack{root};
  while (!stack.empty()) {
    BoolExpr* node = stack.back();
    if (node == nullptr || canonical.contains(node)) {
      stack.pop_back();
      continue;
    }
    if (node->getOp() == Op::VAR) {
      const auto it = merged.find(node->getId());
      canonical[node] = it == merged.end() || it->second == node->getId()
                            ? node : BoolExpr::Var(it->second);
      stack.pop_back();
      continue;
    }
    BoolExpr* left = node->getLeft();
    BoolExpr* right = node->getOp() == Op::NOT ? nullptr : node->getRight();
    const bool leftReady = canonical.contains(left);
    const bool rightReady = right == nullptr || canonical.contains(right);
    if (!leftReady || !rightReady) {
      if (!leftReady) stack.push_back(left);
      if (!rightReady) stack.push_back(right);
      continue;
    }
    BoolExpr* a = canonical.at(left);
    BoolExpr* b = right == nullptr ? nullptr : canonical.at(right);
    const Op op = node->getOp();
    if (op == Op::NOT || op == Op::XOR) {
      canonical[node] = op == Op::NOT ? BoolExpr::Not(a) : BoolExpr::Xor(a, b);
      stack.pop_back();
      continue;
    }
    operands.clear();
    chain.assign({a, b});
    while (!chain.empty() && operands.size() <= kChainOperandLimit) {
      BoolExpr* term = chain.back();
      chain.pop_back();
      if (term->getOp() == op) {
        chain.push_back(term->getLeft());
        chain.push_back(term->getRight());
      } else {
        operands.push_back(term);
      }
    }
    if (operands.size() > kChainOperandLimit) {
      operands.assign({a, b});
    }
    std::sort(operands.begin(), operands.end(),
              [](const BoolExpr* lhs, const BoolExpr* rhs) { return *lhs < *rhs; });
    operands.erase(std::unique(operands.begin(), operands.end()), operands.end());
    BoolExpr* rebuilt = operands.front();
    for (size_t i = 1; i < operands.size(); ++i) {
      rebuilt = op == Op::AND ? BoolExpr::And(rebuilt, operands[i])
                              : BoolExpr::Or(rebuilt, operands[i]);
    }
    canonical[node] = rebuilt;
    stack.pop_back();
  }
  return root == nullptr ? nullptr : canonical.at(root);
}

// Numbers the same rewrite instead of building it: two transitions get one
// number exactly when they are the same formula after the rewrite.  The
// expression factory never frees a node, so building the rewrite of a whole
// design there keeps a second copy of its logic for the rest of the run.  This
// table lives in a few flat arrays that go back to the system with it.
class RewriteTable {
 public:
  explicit RewriteTable(const std::unordered_map<size_t, size_t>& merged)
      : merged_(merged), nodes_(1), slots_(size_t{1} << 16), memo_(size_t{1} << 16) {}

  // Zero means no expression.
  uint32_t numberOf(BoolExpr* root) {
    if (root == nullptr) {
      return 0;
    }
    stack_.assign({root});
    while (!stack_.empty()) {
      BoolExpr* node = stack_.back();
      if (recall(node) != 0) {
        stack_.pop_back();
        continue;
      }
      const Op op = node->getOp();
      if (op == Op::VAR) {
        const auto it = merged_.find(node->getId());
        remember(node, intern({it == merged_.end() ? node->getId() : it->second, 0, op}));
        stack_.pop_back();
        continue;
      }
      BoolExpr* left = node->getLeft();
      BoolExpr* right = op == Op::NOT ? nullptr : node->getRight();
      if (left == nullptr || (op != Op::NOT && right == nullptr)) {
        // Not a formula this table reads: equal only to itself.
        remember(node, intern({reinterpret_cast<uintptr_t>(node), 0, Op::NONE}));
        stack_.pop_back();
        continue;
      }
      const uint32_t a = recall(left);
      const uint32_t b = right == nullptr ? 0 : recall(right);
      if (a == 0 || (right != nullptr && b == 0)) {
        if (a == 0) stack_.push_back(left);
        if (right != nullptr && b == 0) stack_.push_back(right);
        continue;
      }
      remember(node, op == Op::NOT ? negate(a) : op == Op::XOR
                   ? intern({std::min(a, b), std::max(a, b), op}) : chain(op, a, b));
      stack_.pop_back();
    }
    return recall(root);
  }

 private:
  // `a` is a symbol for a leaf and a node number otherwise.
  struct Node {
    uint64_t a = 0;
    uint32_t b = 0;
    Op op = Op::NONE;
    bool operator==(const Node&) const = default;
  };
  struct Memo {
    BoolExpr* expr = nullptr;
    uint32_t number = 0;
  };

  static size_t mix(uint64_t value) {
    value ^= value >> 33;
    value *= 0xFF51AFD7ED558CCDULL;
    return value ^ (value >> 33);
  }

  uint32_t intern(const Node& node) {
    if (nodes_.size() * 2 > slots_.size()) {
      slots_.assign(slots_.size() * 2, 0);
      for (uint32_t number = 1; number < nodes_.size(); ++number) {
        slots_[freeSlot(nodes_[number])] = number;
      }
    }
    const size_t slot = freeSlot(node);
    if (slots_[slot] == 0) {
      slots_[slot] = static_cast<uint32_t>(nodes_.size());
      nodes_.push_back(node);
    }
    return slots_[slot];
  }

  // The slot holding `node`, or the empty one where it belongs.
  size_t freeSlot(const Node& node) const {
    size_t slot = mix(node.a * 0x9E3779B97F4A7C15ULL + node.b * 31 + static_cast<uint64_t>(node.op));
    for (slot &= slots_.size() - 1; slots_[slot] != 0 && !(nodes_[slots_[slot]] == node);
         slot = (slot + 1) & (slots_.size() - 1)) {
    }
    return slot;
  }

  uint32_t negate(uint32_t a) {
    return nodes_[a].op == Op::NOT ? static_cast<uint32_t>(nodes_[a].a) : intern({a, 0, Op::NOT});
  }

  // See canonicalUnderMerge for the chain rule.
  uint32_t chain(Op op, uint32_t a, uint32_t b) {
    constexpr size_t kChainOperandLimit = 16;
    operands_.clear();
    chain_.assign({a, b});
    while (!chain_.empty() && operands_.size() <= kChainOperandLimit) {
      const uint32_t term = chain_.back();
      chain_.pop_back();
      if (nodes_[term].op == op) {
        chain_.push_back(static_cast<uint32_t>(nodes_[term].a));
        chain_.push_back(nodes_[term].b);
      } else {
        operands_.push_back(term);
      }
    }
    if (operands_.size() > kChainOperandLimit) {
      operands_.assign({a, b});
    }
    std::sort(operands_.begin(), operands_.end());
    operands_.erase(std::unique(operands_.begin(), operands_.end()), operands_.end());
    uint32_t rebuilt = operands_.front();
    for (size_t i = 1; i < operands_.size(); ++i) {
      rebuilt = intern({rebuilt, operands_[i], op});
    }
    return rebuilt;
  }

  size_t memoSlot(BoolExpr* expr) const {
    size_t slot = mix(reinterpret_cast<uintptr_t>(expr)) & (memo_.size() - 1);
    while (memo_[slot].expr != nullptr && memo_[slot].expr != expr) {
      slot = (slot + 1) & (memo_.size() - 1);
    }
    return slot;
  }

  uint32_t recall(BoolExpr* expr) const { return memo_[memoSlot(expr)].number; }

  void remember(BoolExpr* expr, uint32_t number) {
    if (++memoSize_ * 4 > memo_.size() * 3) {
      std::vector<Memo> old(memo_.size() * 2);
      old.swap(memo_);
      for (const Memo& entry : old) {
        if (entry.expr != nullptr) {
          memo_[memoSlot(entry.expr)] = entry;
        }
      }
    }
    memo_[memoSlot(expr)] = {expr, number};
  }

  const std::unordered_map<size_t, size_t>& merged_;
  std::vector<Node> nodes_;
  std::vector<uint32_t> slots_;
  std::vector<Memo> memo_;
  size_t memoSize_ = 0;
  std::vector<BoolExpr*> stack_;
  std::vector<uint32_t> chain_;
  std::vector<uint32_t> operands_;
};

// Refines every candidate from one counterexample by replaying it, next to
// random states that also satisfy the hypotheses, for several steps (Mony et
// al., DAC 2005, section 3.2). Pairs of the greatest inductive subset stay
// equal along any such trajectory, so a pair seen differing is outside it.
// This only drops candidates: the certificate remains the SAT proof.
void refuteBySimulation(const KInductionProblem& problem,
                        const TransitionExprResolver& transitions,
                        const std::vector<const InternalRelationCandidate*>& active,
                        const InternalRelationOptions& options,
                        const std::unordered_map<size_t, bool>& counterexample,
                        std::vector<char>& refuted) {
  const bool allowX = options.allowXEqualityInInternalRelations;
  std::mt19937_64 random(1);
  Words state;
  for (size_t symbol : problem.allSymbols) {
    state[symbol] = random();
  }
  for (const auto& value : problem.dualRailStatePairs) {
    // Rails 00 are not a value. Only X/X relations need X states.
    state[value.mayBeZero] = ~state[value.mayBeOne] |
        (allowX ? random() & random() & random() : 0);
  }
  for (const auto& [symbol, bit] : counterexample) {
    state[symbol] = (state[symbol] & ~uint64_t{1}) | (bit ? 1 : 0);
  }
  for (const auto* candidate : active) {
    for (const auto& [lhs, rhs] : candidate->equalities) {
      state[rhs] = state.at(lhs);
    }
  }
  uint64_t valid = ~uint64_t{0};
  for (const auto* candidate : active) {
    valid &= relationWord(*candidate, state, allowX);
  }
  std::vector<size_t> stateSymbols = problem.state0Symbols;
  stateSymbols.insert(stateSymbols.end(), problem.state1Symbols.begin(),
                      problem.state1Symbols.end());
  stateSymbols.insert(stateSymbols.end(), problem.auxiliaryStateSymbols.begin(),
                      problem.auxiliaryStateSymbols.end());
  std::unordered_map<BoolExpr*, uint64_t> values;
  for (size_t step = 0; step < 16 && valid != 0; ++step) {
    values.clear();
    Words next;
    for (size_t symbol : stateSymbols) {
      next[symbol] = state.at(symbol);
    }
    for (const auto* candidate : active) {
      for (const auto& [lhs, rhs] : candidate->equalities) {
        for (size_t symbol : {lhs, rhs}) {
          next[symbol] = evaluateWords(transitions.at(symbol), state, random, values);
        }
      }
    }
    for (const auto& value : problem.dualRailStatePairs) {
      valid &= next.at(value.mayBeOne) | next.at(value.mayBeZero);
    }
    bool progress = false;
    for (size_t i = 0; i < active.size(); ++i) {
      if (!refuted[i] && (~relationWord(*active[i], next, allowX) & valid) != 0) {
        refuted[i] = 1;
        progress = true;
      }
    }
    if (!progress) {
      break;
    }
    state = std::move(next);
  }
}

}  // namespace

std::vector<std::pair<size_t, size_t>> proveInternalRelations(
    const KInductionProblem& problem,
    const std::vector<InternalRelationCandidate>& candidates,
    const InternalRelationOptions& options,
    Config::SolverType solverType) {
  if (!options.learnInternalRelations || candidates.empty()) {
    return {};
  }
  const std::unordered_map<size_t, bool> initial(
      problem.initialStateAssignments.begin(), problem.initialStateAssignments.end());
  TransitionExprResolver transitions(problem);
  std::vector<const InternalRelationCandidate*> active;
  for (const auto& candidate : candidates) {
    bool baseHolds = !candidate.equalities.empty();
    for (const auto& [lhs, rhs] : candidate.equalities) {
      baseHolds &= initial.contains(lhs) && initial.contains(rhs) &&
                   initial.at(lhs) == initial.at(rhs) &&
                   transitions.contains(lhs) && transitions.contains(rhs);
    }
    if (!options.allowXEqualityInInternalRelations) {
      for (const auto& value : candidate.values) {
        baseHolds &= initial.contains(value.mayBeOne) &&
                     initial.contains(value.mayBeZero) &&
                     initial.at(value.mayBeOne) != initial.at(value.mayBeZero);
      }
    }
    if (!baseHolds) {
      continue;
    }
    active.push_back(&candidate);
  }

  // Candidates whose two transitions are one formula once the hypotheses merge
  // their registers hold in the next frame wherever the hypotheses hold in the
  // current one; they are proved by that rewrite alone.  Only the rest is
  // solved, and the rest is given up when its logic is too large for the
  // solver to pay off.  Shared logic counts once, and counting stops at the
  // limit, so a large design is never materialized just to be measured.
  constexpr size_t kLogicNodeLimit = size_t{1} << 23;

  // Hypotheses are applied by literal substitution rather than as solver
  // assumptions: encoding both sides' transitions over the same input
  // literals lets CaDiCaL's congruence closure discharge structurally
  // identical logic without search.  Refinement therefore re-encodes in a
  // fresh solver instead of retracting assumptions.
  //
  // A large design is split rather than skipped (Mishchenko et al., ICCAD
  // 2008, section 3.3): every hypothesis is merged in every partition and each
  // candidate is proved in exactly one, so splitting loses no relation.
  constexpr int kPartitionVariables = 1 << 21;
  const auto dropRefuted = [&active](const std::vector<char>& refuted) {
    size_t kept = 0;
    for (size_t i = 0; i < active.size(); ++i) {
      if (!refuted[i]) {
        active[kept++] = active[i];
      }
    }
    active.resize(kept);
  };
  for (size_t round = 0; round < 64 && !active.empty(); ++round) {
    const auto roundStart = std::chrono::steady_clock::now();
    std::vector<char> refuted(active.size(), 0);
    // Both registers of a hypothesis share one current-frame literal.
    std::unordered_map<size_t, size_t> merged;
    merged.reserve(problem.allSymbols.size());
    for (size_t symbol : problem.allSymbols) {
      merged.emplace(symbol, symbol);
    }
    for (const auto* candidate : active) {
      for (const auto& [lhs, rhs] : candidate->equalities) {
        merged[rhs] = merged.at(lhs);
      }
    }
    std::vector<size_t> pending;
    // The solver reads the rewritten transitions: logic the two registers of a
    // pending candidate have in common is one node there, so it is encoded
    // once and only their real difference is left to search.
    std::unordered_map<size_t, BoolExpr*> rewritten;
    {
      RewriteTable table(merged);
      for (size_t i = 0; i < active.size(); ++i) {
        // Definedness hypotheses are not structural.
        bool structural = options.allowXEqualityInInternalRelations;
        for (const auto& [lhs, rhs] : active[i]->equalities) {
          if (!structural) {
            break;
          }
          const uint32_t left = table.numberOf(transitions.at(lhs));
          structural = left != 0 && left == table.numberOf(transitions.at(rhs));
        }
        if (!structural) {
          pending.push_back(i);
        }
      }
    }
    {
      std::unordered_map<BoolExpr*, BoolExpr*> canonical;
      for (size_t i : pending) {
        for (const auto& [lhs, rhs] : active[i]->equalities) {
          for (size_t symbol : {lhs, rhs}) {
            rewritten[symbol] = canonicalUnderMerge(transitions.at(symbol), merged, canonical);
          }
        }
      }
    }
    const size_t solved = pending.size();
    {
      std::unordered_set<BoolExpr*> counted;
      std::vector<BoolExpr*> stack;
      for (size_t i : pending) {
        for (const auto& [lhs, rhs] : active[i]->equalities) {
          stack.push_back(rewritten.at(lhs));
          stack.push_back(rewritten.at(rhs));
        }
        while (!stack.empty()) {
          BoolExpr* node = stack.back();
          stack.pop_back();
          if (node == nullptr || !counted.insert(node).second) {
            continue;
          }
          stack.push_back(node->getLeft());
          stack.push_back(node->getRight());
        }
        if (counted.size() > kLogicNodeLimit) {
          for (size_t j : pending) {
            refuted[j] = 1;
          }
          pending.clear();
          break;
        }
      }
    }
    for (size_t begin = 0, end = 0; begin < pending.size(); begin = end) {
      // Every limit allows 100 decisions and 1000 ticks for each conflict. A
      // tighter decision or tick limit gives up on queries the conflict limit
      // still allows, and one pair given up unravels every pair that reads it.
      SATSolverWrapper::CadicalWorkBudget budget(100000, 10000000, 100000000);
      SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
      SATSolverWrapper solver(SATSolverWrapper::assumptionSolverTypeFor(solverType));
      // A partition reads a small part of the design, so solver variables
      // are created on first use rather than for every symbol.
      FrameFormulaEncoder current(solver, {}, &merged, true, 0);
      FrameFormulaEncoder next(solver, {}, true);
      std::set<size_t> targets;
      std::vector<int> conclusions;
      std::vector<int> badClause;
      const int firstVariable = solver.newVar();
      for (; end < pending.size() &&
             (end == begin || solver.newVar() - firstVariable < kPartitionVariables);
           ++end) {
        const auto* candidate = active[pending[end]];
        BoolExpr* relation = BoolExpr::createTrue();
        for (const auto& [lhs, rhs] : candidate->equalities) {
          for (size_t symbol : {lhs, rhs}) {
            if (targets.insert(symbol).second) {
              addLiteralEquivalence(solver, next.encode(BoolExpr::Var(symbol)),
                                    current.encode(rewritten.at(symbol)));
            }
          }
          relation = BoolExpr::And(relation, BoolExpr::Not(
              BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs))));
        }
        if (!options.allowXEqualityInInternalRelations) {
          for (const auto& value : candidate->values) {
            relation = BoolExpr::And(relation, BoolExpr::Xor(
                BoolExpr::Var(value.mayBeOne), BoolExpr::Var(value.mayBeZero)));
          }
        }
        // Equalities are tautological after substitution; definedness
        // hypotheses (when X equality is disallowed) remain real constraints.
        solver.addClause({current.encode(relation)});
        conclusions.push_back(next.encode(relation));
        badClause.push_back(-conclusions.back());
        for (const auto& value : candidate->values) {
          solver.addClause({next.encode(BoolExpr::Var(value.mayBeOne)),
                            next.encode(BoolExpr::Var(value.mayBeZero))});
        }
      }
      for (const auto& value : problem.dualRailStatePairs) {
        if (current.leafLits().contains(merged.at(value.mayBeOne)) ||
            current.leafLits().contains(merged.at(value.mayBeZero))) {
          solver.addClause({current.encode(BoolExpr::Var(value.mayBeOne)),
                            current.encode(BoolExpr::Var(value.mayBeZero))});
        }
      }
      const int allTogether = solver.newVar() + 2;
      badClause.push_back(-allTogether);
      solver.addClause(badClause);
      const auto status = solver.solveWithAssumptionsStatus(
          {allTogether}, 10000, 1000000, 10000000);
      if (isSecDiagEnabled()) {
        printf("SEC diag: internal relations round=%zu partition=%zu..%zu vars=%d status=%s "
               "conflicts=%llu decisions=%llu ticks=%llu\n",
               round, begin, end, solver.newVar() - firstVariable,
               status == SATSolverWrapper::SolveStatus::Unsat ? "unsat"
               : status == SATSolverWrapper::SolveStatus::Sat ? "sat" : "unknown",
               static_cast<unsigned long long>(budget.conflictsUsed()),
               static_cast<unsigned long long>(budget.decisionsUsed()),
               static_cast<unsigned long long>(budget.ticksUsed()));
        fflush(stdout);
      }
      if (status == SATSolverWrapper::SolveStatus::Unknown) {
        // Each pair is then its own obligation with its own budget, and only
        // the undecided ones are given up (Mony et al., sections 2 and 4.1).
        for (size_t i = begin; i < end; ++i) {
          refuted[pending[i]] = solver.solveWithAssumptionsStatus(
                                    {-conclusions[i - begin]}, 1000, 100000, 1000000) !=
                                SATSolverWrapper::SolveStatus::Unsat;
        }
      } else if (status == SATSolverWrapper::SolveStatus::Sat) {
        std::unordered_map<size_t, bool> counterexample;
        for (const auto& [symbol, literal] : current.leafLits()) {
          counterexample[symbol] = solver.getLiteralValue(literal);
        }
        refuteBySimulation(problem, transitions, active, options, counterexample, refuted);
        for (size_t i = begin; i < end; ++i) {
          refuted[pending[i]] |= !solver.getLiteralValue(conclusions[i - begin]);
        }
      }
    }
    if (isSecDiagEnabled()) {
      printf("SEC diag: internal relations round=%zu active=%zu structural=%zu refuted=%zu "
             "seconds=%.1f\n", round, active.size(), active.size() - solved,
             static_cast<size_t>(std::count(refuted.begin(), refuted.end(), 1)),
             std::chrono::duration<double>(std::chrono::steady_clock::now() - roundStart).count());
      fflush(stdout);
    }
    // The survivors are a certificate only once a whole round refutes nothing:
    // dropping any hypothesis weakens every other partition's proof.
    if (std::find(refuted.begin(), refuted.end(), 1) == refuted.end()) {
      std::set<std::pair<size_t, size_t>> proved;
      for (const auto* candidate : active) {
        proved.insert(candidate->equalities.begin(), candidate->equalities.end());
      }
      return {proved.begin(), proved.end()};
    }
    dropRefuted(refuted);
  }
  // An unfinished fixed point is not a certificate for any surviving guess.
  return {};
}

void learnInternalStateRelations(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    KInductionProblem& problem,
    const InternalRelationOptions& options,
    Config::SolverType solverType,
    bool diagnostics) {
  if (!options.learnInternalRelations) {
    return;
  }
  const size_t width = problem.usesDualRailStateEncoding ? 2 : 1;
  const std::array<const SequentialDesignModel*, 2> models{&model0, &model1};
  const std::array<const std::vector<size_t>*, 2> symbols{
      &problem.state0Symbols, &problem.state1Symbols};
  using Location = std::pair<size_t, size_t>;
  std::vector<InternalRelationCandidate> candidates;
  auto append = [&](Location lhs, Location rhs) {
    InternalRelationCandidate candidate;
    for (size_t rail = 0; rail < width; ++rail) {
      candidate.equalities.emplace_back(
          symbols[lhs.first]->at(width * lhs.second + rail),
          symbols[rhs.first]->at(width * rhs.second + rail));
    }
    if (width == 2) {
      for (const auto& location : {lhs, rhs}) {
        candidate.values.push_back({
            symbols[location.first]->at(2 * location.second),
            symbols[location.first]->at(2 * location.second + 1)});
      }
    }
    candidates.push_back(std::move(candidate));
  };
  // Names and shared drivers are guesses only. Neither participates in the
  // certificate: the SAT learner checks the original transitions and boot.
  std::unordered_map<std::string, Location> namedStates;
  for (size_t side = 0; side < models.size(); ++side) {
    const auto& model = *models[side];
    std::unordered_map<BoolExpr*, Location> sharedDrivers;
    for (size_t i = 0; i < model.stateBits.size(); ++i) {
      const auto& key = model.stateBits[i];
      const Location location{side, i};
      if (auto name = model.displayNameByKey.find(key);
          name != model.displayNameByKey.end() && !name->second.empty()) {
        auto [existing, inserted] = namedStates.emplace(name->second, location);
        if (!inserted) {
          append(existing->second, location);
        }
      }
      auto [driver, inserted] = sharedDrivers.emplace(
          model.nextStateExprByStateKey.at(key), location);
      if (!inserted) {
        append(driver->second, location);
      }
    }
  }
  const auto proved = proveInternalRelations(problem, candidates, options, solverType);
  problem.sameFrameStateEqualityPairs0.insert(
      problem.sameFrameStateEqualityPairs0.end(), proved.begin(), proved.end());
  // The conjunction costs nodes that are never freed for every relation, and
  // only the exact IMC engine reads it.
  if (!problem.usesLargeDualRailImc()) {
    for (const auto& [lhs, rhs] : proved) {
      BoolExpr* equality = BoolExpr::Not(BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs)));
      problem.learnedInternalRelationInvariant = problem.learnedInternalRelationInvariant == nullptr
          ? equality : BoolExpr::And(problem.learnedInternalRelationInvariant, equality);
    }
  }
  if (diagnostics) {
    printf("SEC summary: internal_relation_candidates=%zu proved_rail_equalities=%zu "
           "allow_x_equality_in_internal_relations=%d\n",
           candidates.size(), proved.size(), options.allowXEqualityInInternalRelations);
    fflush(stdout);
  }
}

}  // namespace KEPLER_FORMAL::SEC
