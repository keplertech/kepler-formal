// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "proof/InternalRelations.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <unordered_map>
#include <unordered_set>

#include "common/BoolExprUtils.h"
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

// A binary obligation is the combinational equivalence of the design under
// the hypotheses, the query a LEC solves with no limit; it is allowed a
// hundred times the budget of a rail check.
constexpr uint64_t kBinaryBudgetScale = 100;

// The binary functions of both designs over the proof symbols: a register
// by its may-be-one rail, an input by its shared symbol. The lazy store
// keeps every next-state function before dual-rail lifting; lifted, the two
// rails of one register are two unrelated copies of its logic to a solver.
class BinaryFormulas {
 public:
  explicit BinaryFormulas(const LazyTransitionStore& store) : store_(store) {
    for (size_t design = 0; design < 2; ++design) {
      symbolMap_[design] = store.localToCombinedByDesign[design];
      for (const auto& [local, rails] :
           store.dualRailStateByLocalSymbolByDesign[design]) {
        symbolMap_[design][local] = rails.mayBeOne;
      }
    }
  }

  bool isRailOne(size_t symbol) const {
    const auto it = store_.sourceByStateSymbol.find(symbol);
    return it != store_.sourceByStateSymbol.end() &&
           it->second.rail == LazyTransitionRail::DualRailOne;
  }

  BoolExpr* remap(size_t design, BoolExpr* local) {
    return remapBoolExprVariables(local, symbolMap_[design], memo_[design]);
  }

  BoolExpr* transition(size_t symbol) {
    auto& expr = transitionBySymbol_[symbol];
    if (expr == nullptr) {
      const auto& source = store_.sourceByStateSymbol.at(symbol);
      expr = remap(source.designIndex, source.localExpr);
    }
    return expr;
  }

 private:
  const LazyTransitionStore& store_;
  std::array<std::unordered_map<size_t, size_t>, 2> symbolMap_;
  std::array<std::unordered_map<BoolExpr*, BoolExpr*>, 2> memo_;
  std::unordered_map<size_t, BoolExpr*> transitionBySymbol_;
};

// Outputs whose two functions agree in binary semantics over the merged
// registers are one Boolean function, so their exact ternary values agree
// as well; the rail obligation of such an output is discharged.
size_t dischargeBinaryEqualOutputs(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const std::vector<std::pair<size_t, size_t>>& proved,
    BinaryFormulas& binary,
    Config::SolverType solverType,
    KInductionProblem& problem) {
  std::unordered_map<size_t, size_t> merged;
  const auto representative = [&merged](size_t symbol) {
    for (auto it = merged.find(symbol); it != merged.end(); it = merged.find(symbol)) {
      symbol = it->second;
    }
    return symbol;
  };
  for (const auto& [lhs, rhs] : proved) {
    const size_t a = representative(lhs);
    const size_t b = representative(rhs);
    if (a != b) {
      merged[std::max(a, b)] = std::min(a, b);
    }
  }
  for (auto& [_, target] : merged) {
    target = representative(target);
  }
  const std::array<const SequentialDesignModel*, 2> models{&model0, &model1};
  std::array<std::unordered_map<std::string, SignalKey>, 2> keyByName;
  for (size_t design = 0; design < 2; ++design) {
    for (const auto& key : models[design]->allObservedOutputs) {
      if (const auto name = models[design]->displayNameByKey.find(key);
          name != models[design]->displayNameByKey.end()) {
        keyByName[design].emplace(name->second, key);
      }
    }
  }
  const size_t count = std::min(problem.observedOutputNames.size(),
                                problem.observedOutputExprs0.size());
  std::vector<char> discharged(count, 0);
  std::vector<size_t> pending;
  std::vector<std::pair<BoolExpr*, BoolExpr*>> functions(count);
  std::unordered_map<BoolExpr*, BoolExpr*> canonical;
  for (size_t i = 0; i < count; ++i) {
    std::array<BoolExpr*, 2> function{nullptr, nullptr};
    for (size_t design = 0; design < 2; ++design) {
      const auto key = keyByName[design].find(problem.observedOutputNames[i]);
      if (key == keyByName[design].end()) {
        break;
      }
      const auto expr = models[design]->observedOutputExprByKey.find(key->second);
      if (expr == models[design]->observedOutputExprByKey.end()) {
        break;
      }
      function[design] =
          canonicalUnderMerge(binary.remap(design, expr->second), merged, canonical);
    }
    if (function[0] == nullptr || function[1] == nullptr) {
      continue;
    }
    if (function[0] == function[1]) {
      discharged[i] = 1;
      continue;
    }
    functions[i] = {function[0], function[1]};
    pending.push_back(i);
  }
  if (!pending.empty()) {
    SATSolverWrapper::CadicalWorkBudget budget(
        100000 * kBinaryBudgetScale, 10000000 * kBinaryBudgetScale,
        100000000 * kBinaryBudgetScale);
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
    SATSolverWrapper solver(SATSolverWrapper::assumptionSolverTypeFor(solverType));
    // The functions already read the representative registers.
    FrameFormulaEncoder current(solver, {}, true, 0);
    std::vector<int> differences;
    std::vector<int> badClause;
    for (size_t i : pending) {
      differences.push_back(current.encode(
          BoolExpr::Xor(functions[i].first, functions[i].second)));
      badClause.push_back(differences.back());
    }
    const int anyDifference = solver.newVar() + 2;
    badClause.push_back(-anyDifference);
    solver.addClause(badClause);
    const auto status = solver.solveWithAssumptionsStatus(
        {anyDifference}, 10000 * kBinaryBudgetScale, 1000000 * kBinaryBudgetScale,
        10000000 * kBinaryBudgetScale);
    for (size_t j = 0; j < pending.size(); ++j) {
      discharged[pending[j]] =
          status == SATSolverWrapper::SolveStatus::Unsat ||
          solver.solveWithAssumptionsStatus({differences[j]}, 10000, 1000000, 10000000) ==
              SATSolverWrapper::SolveStatus::Unsat;
    }
  }
  // Discharged like a structurally identical output, see
  // dischargeStructurallyEquivalentOutputs.
  if (problem.structurallyProvedOutputs.size() != problem.observedOutputExprs0.size()) {
    problem.structurallyProvedOutputs.assign(problem.observedOutputExprs0.size(), false);
  }
  size_t dischargedCount = 0;
  for (size_t i = 0; i < count; ++i) {
    if (!discharged[i]) {
      continue;
    }
    ++dischargedCount;
    problem.observedOutputExprs0[i] = BoolExpr::createTrue();
    problem.observedOutputExprs1[i] = BoolExpr::createTrue();
    if (i < problem.dualRailOutputStrictEqualityExprs.size()) {
      problem.dualRailOutputStrictEqualityExprs[i] = BoolExpr::createTrue();
    }
    problem.structurallyProvedOutputs[i] = true;
  }
  if (dischargedCount > 0) {
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
  }
  return dischargedCount;
}

}  // namespace

std::vector<std::pair<size_t, size_t>> proveInternalRelations(
    const KInductionProblem& problem,
    const std::vector<InternalRelationCandidate>& candidates,
    const InternalRelationOptions& options,
    Config::SolverType solverType,
    bool commonStartBinary) {
  if (!options.learnInternalRelations || candidates.empty()) {
    return {};
  }
  const std::unordered_map<size_t, bool> initial(
      problem.initialStateAssignments.begin(), problem.initialStateAssignments.end());
  TransitionExprResolver transitions(problem);
  // Complement candidates pair unlike rails and are not retried in binary
  // semantics.
  std::optional<BinaryFormulas> binary;
  if (commonStartBinary) {
    if (problem.lazyTransitions == nullptr) {
      return {};
    }
    binary.emplace(*problem.lazyTransitions);
  }
  const auto transitionOf = [&](size_t symbol) {
    return binary ? binary->transition(symbol) : transitions.at(symbol);
  };
  const auto checked = [&](const InternalRelationCandidate& candidate) {
    return std::span(candidate.equalities)
        .first(commonStartBinary ? 1 : candidate.equalities.size());
  };
  std::vector<const InternalRelationCandidate*> active;
  for (const auto& candidate : candidates) {
    bool baseHolds = !candidate.equalities.empty();
    if (binary) {
      baseHolds &= binary->isRailOne(candidate.equalities.front().first) &&
                   binary->isRailOne(candidate.equalities.front().second);
    }
    for (const auto& [lhs, rhs] : candidate.equalities) {
      baseHolds &= initial.contains(lhs) && initial.contains(rhs) &&
                   initial.at(lhs) == initial.at(rhs) &&
                   transitions.contains(lhs) && transitions.contains(rhs);
    }
    // Binary equality from the common initial state needs no defined values.
    if (!options.allowXEqualityInInternalRelations && !commonStartBinary) {
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
  constexpr int kPartitionVariables = 1 << 17;
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
        for (const auto& [lhs, rhs] : checked(*active[i])) {
          if (!structural) {
            break;
          }
          const uint32_t left = table.numberOf(transitionOf(lhs));
          structural = left != 0 && left == table.numberOf(transitionOf(rhs));
        }
        if (!structural) {
          pending.push_back(i);
        }
      }
    }
    {
      std::unordered_map<BoolExpr*, BoolExpr*> canonical;
      for (size_t i : pending) {
        for (const auto& [lhs, rhs] : checked(*active[i])) {
          for (size_t symbol : {lhs, rhs}) {
            rewritten[symbol] = canonicalUnderMerge(transitionOf(symbol), merged, canonical);
          }
        }
      }
    }
    const size_t solved = pending.size();
    {
      std::unordered_set<BoolExpr*> counted;
      std::vector<BoolExpr*> stack;
      for (size_t i : pending) {
        for (const auto& [lhs, rhs] : checked(*active[i])) {
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
      const uint64_t scale = commonStartBinary ? kBinaryBudgetScale : 1;
      SATSolverWrapper::CadicalWorkBudget budget(
          100000 * scale, 10000000 * scale, 100000000 * scale);
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
        for (const auto& [lhs, rhs] : checked(*candidate)) {
          for (size_t symbol : {lhs, rhs}) {
            if (targets.insert(symbol).second) {
              addLiteralEquivalence(solver, next.encode(BoolExpr::Var(symbol)),
                                    current.encode(rewritten.at(symbol)));
            }
          }
          relation = BoolExpr::And(relation, BoolExpr::Not(
              BoolExpr::Xor(BoolExpr::Var(lhs), BoolExpr::Var(rhs))));
        }
        if (!options.allowXEqualityInInternalRelations && !commonStartBinary) {
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
        if (commonStartBinary) {
          continue;
        }
        for (const auto& value : candidate->values) {
          solver.addClause({next.encode(BoolExpr::Var(value.mayBeOne)),
                            next.encode(BoolExpr::Var(value.mayBeZero))});
        }
      }
      for (const auto& value : problem.dualRailStatePairs) {
        if (commonStartBinary) {
          break;
        }
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
          {allTogether}, 10000 * scale, 1000000 * scale, 10000000 * scale);
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
        // The retries share the unscaled partition budget: the scaled one
        // pays for a single large monolithic query, and a partition can hold
        // tens of thousands of pairs (register files kept as flops), so
        // retries charged to the scaled budget would run for hours.
        SATSolverWrapper::CadicalWorkBudget retryBudget(
            100000, 10000000, 100000000);
        SATSolverWrapper::ScopedCadicalWorkBudget retryScope(retryBudget);
        for (size_t i = begin; i < end; ++i) {
          const auto pairStatus = solver.solveWithAssumptionsStatus(
              {-conclusions[i - begin]}, 1000, 100000, 1000000);
          refuted[pending[i]] = pairStatus != SATSolverWrapper::SolveStatus::Unsat;
          if (refuted[pending[i]] && isSecDiagEnabled() && i - begin < 3) {
            printf("SEC diag: internal relation given up (%s): %s\n",
                   pairStatus == SATSolverWrapper::SolveStatus::Sat ? "sat" : "unknown",
                   active[pending[i]]->label.c_str());
          }
        }
      } else if (status == SATSolverWrapper::SolveStatus::Sat) {
        std::unordered_map<size_t, bool> counterexample;
        for (const auto& [symbol, literal] : current.leafLits()) {
          counterexample[symbol] = solver.getLiteralValue(literal);
        }
        if (!commonStartBinary) {
          refuteBySimulation(problem, transitions, active, options, counterexample, refuted);
        }
        for (size_t i = begin; i < end; ++i) {
          const bool refutedHere = !solver.getLiteralValue(conclusions[i - begin]);
          refuted[pending[i]] |= refutedHere;
          if (refutedHere && isSecDiagEnabled() && i - begin < 3) {
            printf("SEC diag: internal relation refuted: %s\n",
                   active[pending[i]]->label.c_str());
          }
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
  // Every relation is proposed once, whichever guess finds it.
  std::set<std::pair<size_t, size_t>> proposed;
  // `swapped` pairs the rails crosswise: in the dual-rail encoding the
  // complement of a value is its rails exchanged (X stays X), so two registers
  // that always hold complementary values are related without any negation.
  // The binary encoding has no such pair; complements are not proposed there.
  auto append = [&](Location lhs, Location rhs, bool swapped) {
    if (swapped && width != 2) {
      return false;
    }
    InternalRelationCandidate candidate;
    for (size_t rail = 0; rail < width; ++rail) {
      candidate.equalities.emplace_back(
          symbols[lhs.first]->at(width * lhs.second + rail),
          symbols[rhs.first]->at(width * rhs.second + (swapped ? 1 - rail : rail)));
    }
    if (!proposed.insert(candidate.equalities.front()).second) {
      return false;
    }
    if (width == 2) {
      for (const auto& location : {lhs, rhs}) {
        candidate.values.push_back({
            symbols[location.first]->at(2 * location.second),
            symbols[location.first]->at(2 * location.second + 1)});
      }
    }
    if (isSecDiagEnabled()) {
      for (const auto& location : {lhs, rhs}) {
        const auto& model = *models[location.first];
        const auto name = model.displayNameByKey.find(model.stateBits.at(location.second));
        candidate.label += (candidate.label.empty() ? "" : " ~ ") +
            (name == model.displayNameByKey.end() ? std::string("?") : name->second);
      }
    }
    candidates.push_back(std::move(candidate));
    return true;
  };
  // Names, net names and shared drivers are guesses only. None participates
  // in the certificate: the SAT learner checks the original transitions and
  // boot.  Net names survive flows that re-cluster or rename instances, such
  // as multibit flop clustering, where pin names do not.
  std::unordered_map<std::string, Location> namedStates;
  std::unordered_map<std::string, Location> netNamedStates;
  size_t nameCandidates = 0;
  size_t netCandidates = 0;
  size_t driverCandidates = 0;
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
          nameCandidates += append(existing->second, location, false);
        }
      }
      if (auto net = model.outputNetNameByKey.find(key);
          net != model.outputNetNameByKey.end() && !net->second.empty()) {
        auto [existing, inserted] = netNamedStates.emplace(net->second, location);
        if (!inserted) {
          netCandidates += append(existing->second, location, false);
        }
      }
      auto [driver, inserted] = sharedDrivers.emplace(
          model.nextStateExprByStateKey.at(key), location);
      if (!inserted) {
        driverCandidates += append(driver->second, location, false);
      }
    }
  }
  auto proved = proveInternalRelations(problem, candidates, options, solverType);

  // Propagation to a fixed point (van Eijk, TCAD 2000): a register one stage
  // behind proved pairs has, once those pairs are merged, the same next-state
  // formula as its counterpart, so it is found by rewriting every transition
  // under the proved relations and numbering the result.  Equal numbers
  // propose an equality; in dual rail, numbers equal across the rails propose
  // a complement (the complemented edges of Kuehlmann et al., TCAD 2002).
  // Each round reaches one register stage deeper and re-proves the whole
  // candidate set, so the result stays a jointly inductive certificate.
  constexpr size_t kPropagationRoundLimit = 32;
  // The first round runs with nothing merged: that is the cross-design
  // structural guess, registers whose transitions already read the same
  // inputs.
  size_t propagatedCandidates = 0;
  size_t propagationRounds = 0;
  {
    TransitionExprResolver transitions(problem);
    for (; propagationRounds < kPropagationRoundLimit; ++propagationRounds) {
      std::unordered_map<size_t, size_t> merged;
      merged.reserve(problem.allSymbols.size());
      for (size_t symbol : problem.allSymbols) {
        merged.emplace(symbol, symbol);
      }
      auto representative = [&merged](size_t symbol) {
        while (merged.at(symbol) != symbol) {
          merged[symbol] = merged.at(merged.at(symbol));
          symbol = merged.at(symbol);
        }
        return symbol;
      };
      for (const auto& [lhs, rhs] : proved) {
        const size_t a = representative(lhs);
        const size_t b = representative(rhs);
        if (a != b) {
          merged[std::max(a, b)] = std::min(a, b);
        }
      }
      for (size_t symbol : problem.allSymbols) {
        merged[symbol] = representative(symbol);
      }
      struct Entry {
        Location location;
        bool swapped;
      };
      std::unordered_map<uint64_t, std::vector<Entry>> buckets;
      {
        RewriteTable table(merged);
        for (size_t side = 0; side < models.size(); ++side) {
          for (size_t i = 0; i < models[side]->stateBits.size(); ++i) {
            std::array<uint32_t, 2> numbers{0, 0};
            bool complete = true;
            for (size_t rail = 0; rail < width && complete; ++rail) {
              const size_t symbol = symbols[side]->at(width * i + rail);
              complete = transitions.contains(symbol) &&
                         (numbers[rail] = table.numberOf(transitions.at(symbol))) != 0;
            }
            if (!complete) {
              continue;
            }
            // Rails ordered by number: two registers whose rails match
            // crosswise land in one bucket with opposite orientation.
            const bool swapped = width == 2 && numbers[0] > numbers[1];
            const uint64_t key = width == 2
                ? (uint64_t{std::min(numbers[0], numbers[1])} << 32) | std::max(numbers[0], numbers[1])
                : numbers[0];
            buckets[key].push_back({{side, i}, swapped});
          }
        }
      }
      size_t added = 0;
      for (const auto& [key, entries] : buckets) {
        const Entry& first = entries.front();
        const size_t firstSymbol = symbols[first.location.first]->at(width * first.location.second);
        for (size_t j = 1; j < entries.size(); ++j) {
          const Entry& other = entries[j];
          const size_t otherSymbol = symbols[other.location.first]->at(width * other.location.second);
          if (representative(firstSymbol) == representative(otherSymbol)) {
            continue;
          }
          const bool swapped = first.swapped != other.swapped;
          // Both rails equal means value and complement coincide (X held):
          // there is nothing to orient, propose the plain equality only.
          if (swapped && (key >> 32) == (key & 0xffffffffu)) {
            continue;
          }
          added += append(first.location, other.location, swapped);
        }
      }
      if (added == 0) {
        break;
      }
      propagatedCandidates += added;
      auto next = proveInternalRelations(problem, candidates, options, solverType);
      if (next.size() <= proved.size()) {
        break;
      }
      proved = std::move(next);
    }
  }

  size_t binaryRoundCandidates = 0;
  size_t binaryRoundProved = 0;
  size_t dischargedOutputs = 0;
  if (options.binaryRelationRound && problem.usesDualRailStateEncoding &&
      problem.lazyTransitions != nullptr) {
    // Shapes of one function propagate X differently, so a true register
    // pair can fail the rail check. Retry those pairs in binary semantics,
    // which proves equality from any common initial state; design 2's
    // register then takes design 1's transition in the proof problem.
    std::unordered_set<size_t> provedSymbols;
    for (const auto& [lhs, rhs] : proved) {
      provedSymbols.insert(lhs);
      provedSymbols.insert(rhs);
    }
    // Every candidate is a hypothesis again: a rail-equal pair is binary
    // equal, and the others need it in the current frame.
    for (const auto& candidate : candidates) {
      for (const auto& [lhs, rhs] : candidate.equalities) {
        if (!provedSymbols.contains(lhs) || !provedSymbols.contains(rhs)) {
          ++binaryRoundCandidates;
          break;
        }
      }
    }
    const auto provedBinary =
        proveInternalRelations(problem, candidates, options, solverType, true);
    for (const auto& [lhs, rhs] : provedBinary) {
      if (provedSymbols.contains(lhs)) {
        continue;
      }
      problem.sharedTransitionPairs.emplace_back(lhs, rhs);
      proved.push_back({lhs, rhs});
      ++binaryRoundProved;
    }
    BinaryFormulas binary(*problem.lazyTransitions);
    dischargedOutputs = dischargeBinaryEqualOutputs(
        model0, model1, proved, binary, solverType, problem);
  }
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
           "allow_x_equality_in_internal_relations=%d name_candidates=%zu "
           "net_name_candidates=%zu driver_candidates=%zu propagated_candidates=%zu "
           "propagation_rounds=%zu\n",
           candidates.size(), proved.size(), options.allowXEqualityInInternalRelations,
           nameCandidates, netCandidates, driverCandidates, propagatedCandidates,
           propagationRounds);
    printf("SEC summary: binary_relation_round candidates=%zu proved_rail_equalities=%zu "
           "discharged_outputs=%zu/%zu\n",
           binaryRoundCandidates, binaryRoundProved, dischargedOutputs,
           problem.observedOutputExprs0.size());
    fflush(stdout);
  }
  dischargeStructurallyEquivalentOutputs(problem, diagnostics);
}

size_t dischargeStructurallyEquivalentOutputs(KInductionProblem& problem,
                                              bool diagnostics) {
  const size_t outputs = problem.observedOutputExprs0.size();
  if (outputs == 0 || problem.observedOutputExprs1.size() != outputs) {
    return 0;
  }
  const bool dualRail = problem.usesDualRailStateEncoding;
  if (dualRail && (problem.dualRailObservedOutputs0.size() != outputs ||
                   problem.dualRailObservedOutputs1.size() != outputs)) {
    return 0;
  }

  // Merge: one representative per equality class of the certified same-frame
  // equalities (the proved register pairs and the definitional complement
  // pairs), so that a register of design 1 is read as its design-0 partner.
  std::unordered_map<size_t, size_t> merged;
  auto representative = [&merged](size_t symbol) {
    auto it = merged.find(symbol);
    while (it != merged.end() && it->second != symbol) {
      symbol = it->second;
      it = merged.find(symbol);
    }
    return symbol;
  };
  for (const auto* pairs :
       {&problem.sameFrameStateEqualityPairs0, &problem.sameFrameStateEqualityPairs1}) {
    for (const auto& [lhs, rhs] : *pairs) {
      const size_t a = representative(lhs);
      const size_t b = representative(rhs);
      if (a != b) {
        merged[std::max(a, b)] = std::min(a, b);
      }
    }
  }
  for (auto& [symbol, root] : merged) {
    root = representative(symbol);
  }

  // Hash: two cones get one number exactly when they are the same formula
  // under the merge.  Numbering never builds the rewritten logic.
  RewriteTable table(merged);
  if (problem.structurallyProvedOutputs.size() != outputs) {
    problem.structurallyProvedOutputs.assign(outputs, false);
  }
  size_t discharged = 0;
  for (size_t i = 0; i < outputs; ++i) {
    if (problem.structurallyProvedOutputs[i]) {
      ++discharged;
      continue;
    }
    bool same = false;
    if (dualRail) {
      const auto& rails0 = problem.dualRailObservedOutputs0[i];
      const auto& rails1 = problem.dualRailObservedOutputs1[i];
      const uint32_t one0 = table.numberOf(rails0.mayBeOne);
      const uint32_t zero0 = table.numberOf(rails0.mayBeZero);
      same = one0 != 0 && zero0 != 0 &&
             one0 == table.numberOf(rails1.mayBeOne) &&
             zero0 == table.numberOf(rails1.mayBeZero);
    } else {
      const uint32_t number0 = table.numberOf(problem.observedOutputExprs0[i]);
      same = number0 != 0 && number0 == table.numberOf(problem.observedOutputExprs1[i]);
    }
    if (!same) {
      continue;
    }
    // The cones are one formula over shared inputs and merged registers, so
    // the output pair is equal in every state the equalities hold in, that is
    // in every reachable state.  Nothing is left for an engine to prove.
    problem.observedOutputExprs0[i] = BoolExpr::createTrue();
    problem.observedOutputExprs1[i] = BoolExpr::createTrue();
    if (problem.dualRailOutputStrictEqualityExprs.size() == outputs) {
      problem.dualRailOutputStrictEqualityExprs[i] = BoolExpr::createTrue();
    }
    problem.structurallyProvedOutputs[i] = true;
    ++discharged;
  }
  if (discharged != 0) {
    // Both encodings build the property as the conjunction of the per-output
    // equalities and share it with the induction property.
    BoolExpr* property = BoolExpr::createTrue();
    for (size_t i = 0; i < outputs; ++i) {
      property = BoolExpr::And(
          property,
          makeEqualityExpr(problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]));
    }
    problem.property = BoolExpr::simplify(property);
    problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
    problem.inductionProperty = problem.property;
    problem.inductionBad = problem.bad;
  }
  if (diagnostics) {
    printf("SEC summary: structurally_proved_outputs=%zu of %zu merged_equalities=%zu\n",
           discharged, outputs, merged.size());
    fflush(stdout);
  }
  return discharged;
}

void mergeProvedStateRelations(KInductionProblem& problem, bool diagnostics) {
  if (problem.sameFrameStateEqualityPairs0.empty() &&
      problem.sameFrameStateEqualityPairs1.empty()) {
    return;
  }
  // The merge changes every formula's support, which lazily remapped problems
  // derive from design-local sources, so such a problem is materialized first.
  // Materializing a very large design costs the memory the lazy store exists
  // to save; those keep their relations as frame constraints only.
  constexpr size_t kMergeStateSymbolLimit = size_t{1} << 17;
  const size_t stateCount = problem.state0Symbols.size() +
                            problem.state1Symbols.size() +
                            problem.auxiliaryStateSymbols.size();
  if (problem.lazyTransitions != nullptr && stateCount > kMergeStateSymbolLimit) {
    if (diagnostics) {
      printf("SEC diag: merged state symbols skipped: %zu state symbols\n", stateCount);
      fflush(stdout);
    }
    return;
  }
  // Union-find over the equality pairs; the smallest symbol represents a class.
  std::unordered_map<size_t, size_t> parent;
  const auto find = [&parent](size_t symbol) {
    size_t root = symbol;
    for (auto it = parent.find(root); it != parent.end() && it->second != root;
         it = parent.find(root)) {
      root = it->second;
    }
    while (symbol != root) {
      const size_t next = parent.at(symbol);
      parent[symbol] = root;
      symbol = next;
    }
    return root;
  };
  for (const auto* pairs : {&problem.sameFrameStateEqualityPairs0,
                            &problem.sameFrameStateEqualityPairs1}) {
    for (const auto& [lhs, rhs] : *pairs) {
      const size_t a = find(lhs);
      const size_t b = find(rhs);
      if (a != b) {
        parent[std::max(a, b)] = std::min(a, b);
      }
    }
  }
  std::unordered_map<size_t, size_t> merged;
  for (const auto& [symbol, _] : parent) {
    if (const size_t representative = find(symbol); representative != symbol) {
      merged.emplace(symbol, representative);
    }
  }
  if (merged.empty()) {
    return;
  }
  if (problem.lazyTransitions != nullptr) {
    TransitionExprResolver transitions(problem);
    auto& store = *problem.lazyTransitions;
    const auto materialize = [&](const std::vector<size_t>& symbols,
                                 std::vector<std::pair<size_t, BoolExpr*>>& eager) {
      for (const size_t symbol : symbols) {
        if (store.sourceByStateSymbol.find(symbol) != store.sourceByStateSymbol.end()) {
          eager.emplace_back(symbol, transitions.at(symbol));
        }
      }
    };
    materialize(problem.state0Symbols, problem.transitions0);
    materialize(problem.state1Symbols, problem.transitions1);
    store.sourceByStateSymbol.clear();
    store.remappedByStateSymbol.clear();
    store.supportByStateSymbol.clear();
    store.nodeCountByStateSymbol.clear();
    for (auto& memo : store.remapMemoByDesign) {
      memo.clear();
    }
    for (auto& memo : store.dualRailRemapMemoByDesign) {
      memo.clear();
    }
  }
  if (!problem.sharedTransitionPairs.empty()) {
    std::unordered_map<size_t, BoolExpr*> transition0BySymbol;
    for (const auto& [symbol, expr] : problem.transitions0) {
      transition0BySymbol.emplace(symbol, expr);
    }
    std::unordered_map<size_t, BoolExpr*> sharedBySymbol;
    for (const auto& [lhs, rhs] : problem.sharedTransitionPairs) {
      if (const auto it = transition0BySymbol.find(lhs); it != transition0BySymbol.end()) {
        sharedBySymbol.emplace(rhs, it->second);
      }
    }
    for (auto& [symbol, expr] : problem.transitions1) {
      if (const auto it = sharedBySymbol.find(symbol); it != sharedBySymbol.end()) {
        expr = it->second;
      }
    }
  }
  std::unordered_map<BoolExpr*, BoolExpr*> canonical;
  const auto rewrite = [&](BoolExpr*& expr) {
    if (expr != nullptr) {
      expr = canonicalUnderMerge(expr, merged, canonical);
    }
  };
  for (auto* transitions : {&problem.transitions0, &problem.transitions1,
                            &problem.auxiliaryTransitions}) {
    for (auto& [_, expr] : *transitions) {
      rewrite(expr);
    }
  }
  for (auto* exprs : {&problem.observedOutputExprs0, &problem.observedOutputExprs1,
                      &problem.dualRailOutputStrictEqualityExprs}) {
    for (auto& expr : *exprs) {
      rewrite(expr);
    }
  }
  rewrite(problem.property);
  rewrite(problem.bad);
  rewrite(problem.inductionProperty);
  rewrite(problem.inductionBad);
  problem.inductionTransitionSupportCache.reset();
  if (diagnostics) {
    size_t sharedOutputs = 0;
    for (size_t i = 0; i < problem.observedOutputExprs0.size() &&
                       i < problem.observedOutputExprs1.size(); ++i) {
      sharedOutputs += problem.observedOutputExprs0[i] == problem.observedOutputExprs1[i] ? 1 : 0;
    }
    printf("SEC diag: merged state symbols=%zu identical_output_formulas=%zu/%zu\n",
           merged.size(), sharedOutputs, problem.observedOutputExprs0.size());
    fflush(stdout);
  }
}
}  // namespace KEPLER_FORMAL::SEC
