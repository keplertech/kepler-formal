// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "search/SimulationSeededSearch.h"

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>

#include "kinduction/SatEncoding.h"
#include "proof/DualRailEncoding.h"
#include "proof/TransitionExprResolver.h"
#include "../../sat/SATSolverWrapper.h"

namespace KEPLER_FORMAL::SEC {
namespace {

constexpr uint32_t kNoNode = std::numeric_limits<uint32_t>::max();
constexpr size_t kNoState = std::numeric_limits<size_t>::max();
constexpr uint8_t kUnknownBits = 3;

enum NodeOp : uint8_t {
  kLeaf = 0,
  kConstZero,
  kConstOne,
  kNot,
  kAnd,
  kOr,
  kXor,
};

// Two-rail gate semantics, identical to buildDualRailBoolExpr.
inline uint8_t ternaryNot(uint8_t a) {
  return static_cast<uint8_t>(((a & 1u) << 1) | ((a >> 1) & 1u));
}
inline uint8_t ternaryAnd(uint8_t a, uint8_t b) {
  return static_cast<uint8_t>(((a & b) & 1u) | ((a | b) & 2u));
}
inline uint8_t ternaryOr(uint8_t a, uint8_t b) {
  return static_cast<uint8_t>(((a | b) & 1u) | ((a & b) & 2u));
}
inline uint8_t ternaryXor(uint8_t a, uint8_t b) {
  const uint8_t a1 = a & 1u;
  const uint8_t a0 = (a >> 1) & 1u;
  const uint8_t b1 = b & 1u;
  const uint8_t b0 = (b >> 1) & 1u;
  return static_cast<uint8_t>(
      ((a1 & b0) | (a0 & b1)) | (((a1 & b1) | (a0 & b0)) << 1));
}

// Renames the variables listed in `varMap` and leaves every other variable
// as it is (the shared remap helper requires a complete map).
BoolExpr* renameVariables(
    BoolExpr* root,
    const std::unordered_map<size_t, size_t>& varMap,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo) {
  if (root == nullptr) {
    return nullptr;  // LCOV_EXCL_LINE
  }
  struct Frame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };
  std::vector<Frame> stack;
  stack.push_back({root, false});
  while (!stack.empty()) {
    const Frame current = stack.back();
    stack.pop_back();
    BoolExpr* expr = current.expr;
    if (expr == nullptr || memo.find(expr) != memo.end()) {
      continue;
    }
    const Op op = expr->getOp();
    if (op == Op::VAR) {
      const auto it = varMap.find(expr->getId());
      memo.emplace(expr, it == varMap.end() || it->second == expr->getId()
                             ? expr
                             : BoolExpr::Var(it->second));
      continue;
    }
    if (!current.visited) {
      stack.push_back({expr, true});
      if (op != Op::NOT && expr->getRight() != nullptr &&
          memo.find(expr->getRight()) == memo.end()) {
        stack.push_back({expr->getRight(), false});
      }
      if (expr->getLeft() != nullptr && memo.find(expr->getLeft()) == memo.end()) {
        stack.push_back({expr->getLeft(), false});
      }
      continue;
    }
    BoolExpr* left = memo.at(expr->getLeft());
    BoolExpr* right = op == Op::NOT ? nullptr : memo.at(expr->getRight());
    BoolExpr* renamed = nullptr;
    switch (op) {
      case Op::NOT:
        renamed = BoolExpr::Not(left);
        break;
      case Op::AND:
        renamed = BoolExpr::And(left, right);
        break;
      case Op::OR:
        renamed = BoolExpr::Or(left, right);
        break;
      case Op::XOR:
        renamed = BoolExpr::Xor(left, right);
        break;
      default:
        throw std::runtime_error(  // LCOV_EXCL_LINE
            "Unsupported BoolExpr operator in variable renaming");  // LCOV_EXCL_LINE
    }
    memo.emplace(expr, renamed);
  }
  return memo.at(root);
}

double elapsedSeconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now() - start)
      .count();
}

}  // namespace

TernaryModelSimulator::TernaryModelSimulator(const SequentialDesignModel& model)
    : model_(model) {
  size_t maxVar = 0;
  for (const auto& [key, var] : model.inputVarByKey) {
    maxVar = std::max(maxVar, var);
  }
  ensureVarCapacity(maxVar);
  state_.assign(model.stateBits.size(), TernaryValue::Unknown);
  toggled_.assign(model.stateBits.size(), false);
  nextRootByState_.assign(model.stateBits.size(), kNoNode);
  for (size_t index = 0; index < model.stateBits.size(); ++index) {
    const auto& key = model.stateBits[index];
    if (const auto varIt = model.inputVarByKey.find(key);
        varIt != model.inputVarByKey.end()) {
      ensureVarCapacity(varIt->second);
      stateIndexByLocalVar_[varIt->second] = index;
    }
    if (const auto nextIt = model.nextStateExprByStateKey.find(key);
        nextIt != model.nextStateExprByStateKey.end() &&
        nextIt->second != nullptr) {
      nextRootByState_[index] = compile(nextIt->second);
    }
  }
  for (const auto& [key, expr] : model.observedOutputExprByKey) {
    if (expr != nullptr) {
      outputRootByKey_.emplace(key, compile(expr));
    }
  }
  values_.assign(nodes_.size(), kUnknownBits);
  nodeIndexByExpr_.clear();
  resetToInitialState();
}

void TernaryModelSimulator::ensureVarCapacity(size_t localVar) {
  if (localVar < leafNodeByLocalVar_.size()) {
    return;
  }
  leafNodeByLocalVar_.resize(localVar + 1, kNoNode);
  stateIndexByLocalVar_.resize(localVar + 1, kNoState);
  inputValueByLocalVar_.resize(localVar + 1, kUnknownBits);
}

uint32_t TernaryModelSimulator::compile(BoolExpr* root) {
  if (const auto it = nodeIndexByExpr_.find(root); it != nodeIndexByExpr_.end()) {
    return it->second;
  }
  struct Frame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };
  std::vector<Frame> stack;
  stack.push_back({root, false});
  while (!stack.empty()) {
    const Frame current = stack.back();
    stack.pop_back();
    BoolExpr* expr = current.expr;
    if (expr == nullptr || nodeIndexByExpr_.find(expr) != nodeIndexByExpr_.end()) {
      continue;
    }
    const Op op = expr->getOp();
    if (op == Op::VAR) {
      const size_t id = expr->getId();
      Node node;
      node.left = kNoNode;
      node.right = kNoNode;
      if (id == 0) {
        node.op = kConstZero;
      } else if (id == 1) {
        node.op = kConstOne;
      } else {
        ensureVarCapacity(id);
        if (leafNodeByLocalVar_[id] != kNoNode) {
          nodeIndexByExpr_.emplace(expr, leafNodeByLocalVar_[id]);
          continue;
        }
        node.op = kLeaf;
        node.left = static_cast<uint32_t>(id);
      }
      const auto index = static_cast<uint32_t>(nodes_.size());
      nodes_.push_back(node);
      nodeIndexByExpr_.emplace(expr, index);
      if (node.op == kLeaf) {
        leafNodeByLocalVar_[id] = index;
      }
      continue;
    }
    if (op == Op::NONE) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          "Unsupported BoolExpr operator in ternary simulation");  // LCOV_EXCL_LINE
    }
    if (!current.visited) {
      stack.push_back({expr, true});
      if (op != Op::NOT && expr->getRight() != nullptr &&
          nodeIndexByExpr_.find(expr->getRight()) == nodeIndexByExpr_.end()) {
        stack.push_back({expr->getRight(), false});
      }
      if (expr->getLeft() != nullptr &&
          nodeIndexByExpr_.find(expr->getLeft()) == nodeIndexByExpr_.end()) {
        stack.push_back({expr->getLeft(), false});
      }
      continue;
    }
    Node node;
    node.left = nodeIndexByExpr_.at(expr->getLeft());
    node.right = op == Op::NOT ? kNoNode : nodeIndexByExpr_.at(expr->getRight());
    switch (op) {
      case Op::NOT:
        node.op = kNot;
        break;
      case Op::AND:
        node.op = kAnd;
        break;
      case Op::OR:
        node.op = kOr;
        break;
      case Op::XOR:
        node.op = kXor;
        break;
      default:
        throw std::runtime_error(  // LCOV_EXCL_LINE
            "Unsupported BoolExpr operator in ternary simulation");  // LCOV_EXCL_LINE
    }
    const auto index = static_cast<uint32_t>(nodes_.size());
    nodes_.push_back(node);
    nodeIndexByExpr_.emplace(expr, index);
  }
  return nodeIndexByExpr_.at(root);
}

void TernaryModelSimulator::resetToInitialState() {
  for (size_t index = 0; index < state_.size(); ++index) {
    const auto valueIt = model_.initialStateValueByKey.find(model_.stateBits[index]);
    state_[index] = valueIt == model_.initialStateValueByKey.end()
                        ? TernaryValue::Unknown
                        : (valueIt->second ? TernaryValue::One : TernaryValue::Zero);
  }
  combinationalValid_ = false;
}

void TernaryModelSimulator::setInput(size_t localVar, TernaryValue value) {
  ensureVarCapacity(localVar);
  inputValueByLocalVar_[localVar] = static_cast<uint8_t>(value);
  combinationalValid_ = false;
}

void TernaryModelSimulator::clearInputs() {
  std::fill(inputValueByLocalVar_.begin(), inputValueByLocalVar_.end(), kUnknownBits);
  combinationalValid_ = false;
}

void TernaryModelSimulator::setState(const std::vector<TernaryValue>& state) {
  if (state.size() != state_.size()) {
    throw std::runtime_error("Ternary simulator state size mismatch");  // LCOV_EXCL_LINE
  }
  state_ = state;
  combinationalValid_ = false;
}

void TernaryModelSimulator::evaluateCombinational() {
  for (size_t index = 0; index < nodes_.size(); ++index) {
    const Node& node = nodes_[index];
    switch (node.op) {
      case kLeaf: {
        const size_t localVar = node.left;
        const size_t stateIndex = stateIndexByLocalVar_[localVar];
        values_[index] = stateIndex != kNoState
                             ? static_cast<uint8_t>(state_[stateIndex])
                             : inputValueByLocalVar_[localVar];
        break;
      }
      case kConstZero:
        values_[index] = static_cast<uint8_t>(TernaryValue::Zero);
        break;
      case kConstOne:
        values_[index] = static_cast<uint8_t>(TernaryValue::One);
        break;
      case kNot:
        values_[index] = ternaryNot(values_[node.left]);
        break;
      case kAnd:
        values_[index] = ternaryAnd(values_[node.left], values_[node.right]);
        break;
      case kOr:
        values_[index] = ternaryOr(values_[node.left], values_[node.right]);
        break;
      case kXor:
        values_[index] = ternaryXor(values_[node.left], values_[node.right]);
        break;
      default:
        break;  // LCOV_EXCL_LINE
    }
  }
  combinationalValid_ = true;
}

bool TernaryModelSimulator::step() {
  evaluateCombinational();
  bool changed = false;
  for (size_t index = 0; index < state_.size(); ++index) {
    const uint32_t root = nextRootByState_[index];
    const TernaryValue next = root == kNoNode
                                  ? TernaryValue::Unknown
                                  : static_cast<TernaryValue>(values_[root]);
    if (next != state_[index]) {
      changed = true;
      // A toggle is a change between two known values; becoming known (reset,
      // first write) or unknown is initialization, not activity.
      if (isKnownTernaryValue(next) && isKnownTernaryValue(state_[index])) {
        toggled_[index] = true;
      }
    }
    state_[index] = next;
  }
  combinationalValid_ = false;
  return changed;
}

void TernaryModelSimulator::clearToggles() {
  std::fill(toggled_.begin(), toggled_.end(), false);
}

size_t TernaryModelSimulator::knownStateBitCount() const {
  size_t known = 0;
  for (const auto value : state_) {
    known += isKnownTernaryValue(value) ? 1 : 0;
  }
  return known;
}

TernaryValue TernaryModelSimulator::evaluateOutput(const SignalKey& outputKey) {
  const auto rootIt = outputRootByKey_.find(outputKey);
  if (rootIt == outputRootByKey_.end()) {
    return TernaryValue::Unknown;
  }
  if (!combinationalValid_) {
    evaluateCombinational();
  }
  return static_cast<TernaryValue>(values_[rootIt->second]);
}

namespace {

struct RailLits {
  int one = 0;
  int zero = 0;
};

// Small Tseitin builder with constant folding for the two-rail word logic of
// the memory ports. Literals follow the SAT wrapper convention (ids >= 2).
class RailCnfBuilder {
 public:
  explicit RailCnfBuilder(SATSolverWrapper& solver)
      : solver_(solver), trueLit_(solver.newVar() + 2) {
    solver_.addClause({trueLit_});
  }

  int trueLit() const { return trueLit_; }
  int falseLit() const { return -trueLit_; }

  int andLit(int a, int b) {
    if (a == falseLit() || b == falseLit() || a == -b) {
      return falseLit();
    }
    if (a == trueLit() || a == b) {
      return b;
    }
    if (b == trueLit()) {
      return a;
    }
    const int out = solver_.newVar() + 2;
    solver_.addClause({-out, a});
    solver_.addClause({-out, b});
    solver_.addClause({out, -a, -b});
    return out;
  }

  int orLit(int a, int b) { return -andLit(-a, -b); }

  RailLits constRail(bool value) const {
    return value ? RailLits{trueLit(), falseLit()} : RailLits{falseLit(), trueLit()};
  }

  RailLits railAnd(const RailLits& a, const RailLits& b) {
    return {andLit(a.one, b.one), orLit(a.zero, b.zero)};
  }

  // May-be-equal / may-differ rails of two ternary bits.
  RailLits railEqual(const RailLits& a, const RailLits& b) {
    return {orLit(andLit(a.one, b.one), andLit(a.zero, b.zero)),
            orLit(andLit(a.one, b.zero), andLit(a.zero, b.one))};
  }

  RailLits wordEqual(const std::vector<RailLits>& a, const std::vector<RailLits>& b) {
    RailLits result{trueLit(), falseLit()};
    const size_t width = std::min(a.size(), b.size());
    for (size_t bit = 0; bit < width; ++bit) {
      const RailLits equal = railEqual(a[bit], b[bit]);
      result.one = andLit(result.one, equal.one);
      result.zero = orLit(result.zero, equal.zero);
    }
    return result;
  }

  // May-be-equal rail of a ternary word against a constant address.
  int wordMayEqualConstant(const std::vector<RailLits>& a, size_t value) {
    int result = trueLit();
    for (size_t bit = 0; bit < a.size(); ++bit) {
      const bool expected = ((value >> bit) & size_t{1}) != 0;
      result = andLit(result, expected ? a[bit].one : a[bit].zero);
    }
    return result;
  }

  RailLits ite(const RailLits& condition, const RailLits& whenTrue, const RailLits& whenFalse) {
    return {orLit(andLit(condition.one, whenTrue.one), andLit(condition.zero, whenFalse.one)),
            orLit(andLit(condition.one, whenTrue.zero), andLit(condition.zero, whenFalse.zero))};
  }

 private:
  SATSolverWrapper& solver_;
  int trueLit_;
};

RailLits encodeRails(FrameFormulaEncoder& encoder, const DualRailBoolExpr& expr) {
  return {encoder.encode(expr.mayBeOne), encoder.encode(expr.mayBeZero)};
}

std::vector<RailLits> encodeRailWord(
    FrameFormulaEncoder& encoder, const std::vector<DualRailBoolExpr>& word) {
  std::vector<RailLits> lits;
  lits.reserve(word.size());
  for (const auto& bit : word) {
    lits.push_back(encodeRails(encoder, bit));
  }
  return lits;
}

struct WriteEvent {
  RailLits enable;
  std::vector<RailLits> address;
  std::vector<RailLits> data;
  std::vector<RailLits> mask;  // empty: every bit is written
};

struct MemoryFrameEvents {
  std::vector<WriteEvent> writes;  // port order
  std::optional<RailLits> reset;
};

struct ReadOutputRef {
  size_t memory = 0;
  size_t port = 0;
  size_t bit = 0;
};

struct FrontierBoundedSearchResult {
  size_t badFrame = 0;
  // inputsByFrame[frame][alignedIndex]
  std::vector<std::vector<bool>> inputsByFrame;
  // The hit was not on the target itself but on an enlargement cube or an
  // intermediate target supplied by the caller.
  bool viaIntermediateTarget = false;
};

// Target enlargement (Baumgartner, Kuehlmann and Abraham, CAV 2002) for the
// bounded search. Structural form: `structuralSteps` extra frames in which
// only the target's backward cone is unrolled, so the check at the last full
// frame covers every state within that many cycles of the target. State-set
// form: up to `cubeBudget` generalized preimage cubes per step, each a
// conjunction of rail literals from which the target is reached under some
// input, enumerated by SAT (McMillan, CAV 2002 style model enumeration with
// assumption-based lifting) and added to the target at every frame.
struct TargetEnlargement {
  size_t structuralSteps = 0;
  size_t cubeSteps = 0;
  size_t cubeBudget = 0;
};

using StateCube = std::vector<std::pair<size_t, bool>>;  // (rail symbol, value)

void logSeeded(bool diagnostics, const char* format, ...)
    __attribute__((format(printf, 2, 3)));

void logSeeded(bool diagnostics, const char* format, ...) {
  if (!diagnostics) {
    return;
  }
  va_list args;
  va_start(args, format);
  fputs("SEC seeded search: ", stdout);
  vfprintf(stdout, format, args);
  fputc('\n', stdout);
  fflush(stdout);
  va_end(args);
}

// Diagnostics: how a set of design-1 registers spreads over the hierarchy.
// Prints the most populated instance prefixes (first `depth` path elements).
void logRegisterHistogram(bool diagnostics, const char* what, const SequentialDesignModel& model,
                          const std::vector<size_t>& indices, size_t depth = 3,
                          size_t maxRows = 12, size_t maxNames = 16) {
  if (!diagnostics || indices.empty()) {
    return;
  }
  std::map<std::string, size_t> counts;
  for (const auto index : indices) {
    if (index >= model.stateBits.size()) {
      continue;
    }
    const auto it = model.displayNameByKey.find(model.stateBits[index]);
    std::string name = it == model.displayNameByKey.end() ? std::string("?") : it->second;
    size_t pos = 0;
    for (size_t element = 0; element < depth; ++element) {
      const auto dot = name.find('.', pos);
      if (dot == std::string::npos) {
        pos = name.size();
        break;
      }
      pos = dot + 1;
    }
    counts[name.substr(0, pos == name.size() ? pos : pos - 1)] += 1;
  }
  std::vector<std::pair<std::string, size_t>> rows(counts.begin(), counts.end());
  std::stable_sort(rows.begin(), rows.end(),
                   [](const auto& a, const auto& b) { return a.second > b.second; });
  for (size_t i = 0; i < rows.size() && i < maxRows; ++i) {
    logSeeded(diagnostics, "%s by instance: %s = %zu", what, rows[i].first.c_str(), rows[i].second);
  }
  if (indices.size() <= maxNames) {
    for (const auto index : indices) {
      const auto it = model.displayNameByKey.find(model.stateBits[index]);
      logSeeded(diagnostics, "%s register: %s", what,
                it == model.displayNameByKey.end() ? "?" : it->second.c_str());
    }
  }
}

bool supportTouches(const std::set<size_t>& support, const std::unordered_set<size_t>& symbols) {
  for (const auto symbol : support) {
    if (symbols.find(symbol) != symbols.end()) {
      return true;
    }
  }
  return false;
}

// Maps one design's local variables into a private rail symbol space.
class LocalRailMapper final : public DualRailVariableMapper {
 public:
  LocalRailMapper(const std::vector<DualRailSymbolPair>& railsByLocalVar,
                  const std::vector<size_t>& inputSymbolByLocalVar)
      : railsByLocalVar_(railsByLocalVar),
        inputSymbolByLocalVar_(inputSymbolByLocalVar) {}
  ~LocalRailMapper() override = default;

  DualRailBoolExpr mapVariable(size_t symbol) override {
    if (symbol < railsByLocalVar_.size() &&
        railsByLocalVar_[symbol].mayBeOne != kNoState) {
      return {BoolExpr::Var(railsByLocalVar_[symbol].mayBeOne),
              BoolExpr::Var(railsByLocalVar_[symbol].mayBeZero)};
    }
    if (symbol < inputSymbolByLocalVar_.size() &&
        inputSymbolByLocalVar_[symbol] != kNoState) {
      return {BoolExpr::Var(inputSymbolByLocalVar_[symbol]),
              BoolExpr::Not(BoolExpr::Var(inputSymbolByLocalVar_[symbol]))};
    }
    throw std::runtime_error(  // LCOV_EXCL_LINE
        "Seeded search formula contains unpublished internal support v" +  // LCOV_EXCL_LINE
        std::to_string(symbol));  // LCOV_EXCL_LINE
  }

 private:
  const std::vector<DualRailSymbolPair>& railsByLocalVar_;
  const std::vector<size_t>& inputSymbolByLocalVar_;
};

// Both values binary and different.
BoolExpr* binaryDifferent(const DualRailBoolExpr& a, const DualRailBoolExpr& b) {
  BoolExpr* aOne = BoolExpr::And(a.mayBeOne, BoolExpr::Not(a.mayBeZero));
  BoolExpr* aZero = BoolExpr::And(a.mayBeZero, BoolExpr::Not(a.mayBeOne));
  BoolExpr* bOne = BoolExpr::And(b.mayBeOne, BoolExpr::Not(b.mayBeZero));
  BoolExpr* bZero = BoolExpr::And(b.mayBeZero, BoolExpr::Not(b.mayBeOne));
  return BoolExpr::Or(BoolExpr::And(aOne, bZero), BoolExpr::And(aZero, bOne));
}

// The two designs run in lockstep from the selected state until some register
// or output takes different next values. On one design's state, design 2's
// next-state formula of a paired register, renamed to design 1's symbols,
// either is design 1's formula or differs; the disjunction of "binary and
// different" over the differing pairs is the divergence predicate. Searching
// it on one copy of the design avoids the two-copy unrolling, whose every
// bounded check is a full combinational equivalence proof.
struct DivergenceProblem {
  KInductionProblem problem;                   // design 1 only, private symbols
  std::unordered_map<size_t, bool> frontier;   // rail -> selected-state value
  size_t pairedRegisters = 0;
  size_t unpairedRegisters = 0;                // design-2 registers without a partner
  size_t differingRegisters = 0;
  size_t incomparableRegisters = 0;            // differ but read unpaired registers
  size_t differingOutputs = 0;
  bool frontierDiffers = false;                // a paired register already differs
  std::vector<std::string> unpairedNames;      // design-2 registers without a partner
  std::vector<std::string> incomparableNames;  // design-1 side of incomparable pairs
  std::vector<size_t> differingIndices;        // design-1 state indices of differing pairs
  std::vector<DualRailSymbolPair> railsByStateIndex;  // design-1 state index -> rails
  // The register part of the divergence predicate (no output terms): its cone
  // stays local to the changed logic, which makes it the guide for waypoints.
  BoolExpr* registerDivergence = nullptr;
  // The origins of the register differences: the deepest subformulas where
  // design 1's and the renamed design 2's next-state functions stop matching
  // (the changed logic itself), as "binary and different" terms. Its support
  // is the changed logic's fan-in only, which makes it the preferred guide.
  BoolExpr* originDivergence = nullptr;
  size_t originPairs = 0;

  std::unordered_map<size_t, bool> frontierFor(const std::vector<TernaryValue>& state0) const {
    std::unordered_map<size_t, bool> result;
    result.reserve(2 * railsByStateIndex.size());
    for (size_t index = 0; index < railsByStateIndex.size() && index < state0.size(); ++index) {
      const auto bits = static_cast<uint8_t>(state0[index]);
      result[railsByStateIndex[index].mayBeOne] = (bits & 1u) != 0;
      result[railsByStateIndex[index].mayBeZero] = (bits & 2u) != 0;
    }
    return result;
  }
};

DivergenceProblem buildDivergenceProblem(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs,
    const std::vector<TernaryValue>& state0,
    const std::vector<TernaryValue>& state1,
    const std::vector<std::optional<bool>>& resetActiveValue) {
  DivergenceProblem out;
  const size_t inputCount = alignedInputs.names.size();

  // Design 2 -> design 1 renaming. Unpaired design-2 symbols move above every
  // design-1 symbol so they can never coincide with one by accident.
  size_t maxVar0 = 0;
  for (const auto& [key, var] : model0.inputVarByKey) {
    maxVar0 = std::max(maxVar0, var);
  }
  const size_t offset = maxVar0 + 2;
  std::unordered_map<size_t, size_t> rename;
  rename.reserve(model1.inputVarByKey.size());
  for (const auto& [key, var] : model1.inputVarByKey) {
    rename[var] = var + offset;
  }
  for (size_t index = 0; index < inputCount; ++index) {
    const auto it0 = model0.inputVarByKey.find(alignedInputs.keys0[index]);
    const auto it1 = model1.inputVarByKey.find(alignedInputs.keys1[index]);
    if (it0 != model0.inputVarByKey.end() && it1 != model1.inputVarByKey.end()) {
      rename[it1->second] = it0->second;
    }
  }

  // Register pairing by output net name first (RTL signal names survive a
  // local edit), then by display name (which carries instance numbering that
  // shifts where a module changed). Names are only hints: the formula
  // comparison below decides what differs.
  std::unordered_map<std::string, size_t> indexByDisplayName1;
  std::unordered_map<std::string, size_t> indexByNetName1;
  const auto addName = [](std::unordered_map<std::string, size_t>& byName,
                          const std::string& name, size_t index) {
    const auto [it, inserted] = byName.emplace(name, index);
    if (!inserted) {
      it->second = kNoState;  // ambiguous name
    }
  };
  for (size_t index = 0; index < model1.stateBits.size(); ++index) {
    const auto& key = model1.stateBits[index];
    if (const auto it = model1.displayNameByKey.find(key); it != model1.displayNameByKey.end()) {
      addName(indexByDisplayName1, it->second, index);
    }
    if (const auto it = model1.outputNetNameByKey.find(key);
        it != model1.outputNetNameByKey.end() && !it->second.empty()) {
      addName(indexByNetName1, it->second, index);
    }
  }
  const auto lookup = [](const std::unordered_map<std::string, size_t>& byName,
                         const std::string& name) {
    const auto it = byName.find(name);
    return it == byName.end() ? kNoState : it->second;
  };
  // Candidate partners of every design-1 register.
  std::vector<size_t> netPartner(model0.stateBits.size(), kNoState);
  std::vector<size_t> displayPartner(model0.stateBits.size(), kNoState);
  for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
    const auto& key0 = model0.stateBits[index0];
    if (const auto it = model0.outputNetNameByKey.find(key0);
        it != model0.outputNetNameByKey.end() && !it->second.empty()) {
      netPartner[index0] = lookup(indexByNetName1, it->second);
    }
    if (const auto it = model0.displayNameByKey.find(key0);
        it != model0.displayNameByKey.end()) {
      displayPartner[index0] = lookup(indexByDisplayName1, it->second);
    }
  }
  std::vector<size_t> partnerOf0(model0.stateBits.size(), kNoState);
  std::vector<size_t> partnerOf1(model1.stateBits.size(), kNoState);
  const auto assign = [&](size_t index0, size_t index1) {
    partnerOf0[index0] = index1;
    partnerOf1[index1] = index0;
  };
  const auto release = [&](size_t index0) {
    if (partnerOf0[index0] != kNoState) {
      partnerOf1[partnerOf0[index0]] = kNoState;
      partnerOf0[index0] = kNoState;
    }
  };
  for (const auto* candidates : {&netPartner, &displayPartner}) {
    for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
      const size_t index1 = (*candidates)[index0];
      if (partnerOf0[index0] == kNoState && index1 != kNoState &&
          partnerOf1[index1] == kNoState) {
        assign(index0, index1);
      }
    }
  }
  const auto localVar1 = [&](size_t index1) {
    return model1.inputVarByKey.at(model1.stateBits[index1]);
  };
  const auto localVar0 = [&](size_t index0) {
    return model0.inputVarByKey.at(model0.stateBits[index0]);
  };
  const auto nextExpr = [](const SequentialDesignModel& model, size_t index) -> BoolExpr* {
    const auto it = model.nextStateExprByStateKey.find(model.stateBits[index]);
    return it == model.nextStateExprByStateKey.end() ? nullptr : it->second;
  };
  const auto renamingFromPairs = [&]() {
    std::unordered_map<size_t, size_t> map = rename;
    for (size_t index0 = 0; index0 < partnerOf0.size(); ++index0) {
      if (partnerOf0[index0] != kNoState) {
        map[localVar1(partnerOf0[index0])] = localVar0(index0);
      }
    }
    return map;
  };
  // A net name follows the wire a register drives, a display name follows the
  // instance; a local edit can change either. Where the net-name partner's
  // formula differs but the display-name partner's is identical, the latter
  // is the register that did not change.
  for (size_t pass = 0; pass < 3; ++pass) {
    const auto currentRename = renamingFromPairs();
    std::unordered_map<BoolExpr*, BoolExpr*> passMemo;
    bool changed = false;
    for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
      const size_t index1 = partnerOf0[index0];
      const size_t alternative = displayPartner[index0];
      if (index1 == kNoState || alternative == kNoState || alternative == index1) {
        continue;
      }
      BoolExpr* next0 = nextExpr(model0, index0);
      BoolExpr* next1 = nextExpr(model1, index1);
      BoolExpr* nextAlternative = nextExpr(model1, alternative);
      if (next0 == nullptr || next1 == nullptr || nextAlternative == nullptr ||
          renameVariables(next1, currentRename, passMemo) == next0) {
        continue;
      }
      const size_t displaced0 = partnerOf1[alternative];
      if (displaced0 != kNoState) {
        BoolExpr* displacedNext0 = nextExpr(model0, displaced0);
        if (displacedNext0 != nullptr &&
            renameVariables(nextAlternative, currentRename, passMemo) == displacedNext0) {
          continue;  // the alternative already has an identical partner
        }
      }
      auto trial = currentRename;
      trial[localVar1(alternative)] = localVar0(index0);
      trial[localVar1(index1)] = localVar1(index1) + offset;
      std::unordered_map<BoolExpr*, BoolExpr*> trialMemo;
      if (renameVariables(nextAlternative, trial, trialMemo) != next0) {
        continue;
      }
      if (displaced0 != kNoState) {
        release(displaced0);
      }
      release(index0);
      assign(index0, alternative);
      changed = true;
    }
    if (!changed) {
      break;
    }
    // Registers displaced by a re-pairing take any free candidate.
    for (const auto* candidates : {&netPartner, &displayPartner}) {
      for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
        const size_t index1 = (*candidates)[index0];
        if (partnerOf0[index0] == kNoState && index1 != kNoState &&
            partnerOf1[index1] == kNoState) {
          assign(index0, index1);
        }
      }
    }
  }
  std::vector<std::pair<size_t, size_t>> pairs;
  std::vector<bool> used1(model1.stateBits.size(), false);
  for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
    if (partnerOf0[index0] != kNoState) {
      pairs.emplace_back(index0, partnerOf0[index0]);
      used1[partnerOf0[index0]] = true;
    }
  }
  for (const auto& [index0, index1] : pairs) {
    const auto& key0 = model0.stateBits[index0];
    rename[model1.inputVarByKey.at(model1.stateBits[index1])] =
        model0.inputVarByKey.at(key0);
    if (isKnownTernaryValue(state0[index0]) && isKnownTernaryValue(state1[index1]) &&
        state0[index0] != state1[index1]) {
      out.frontierDiffers = true;
    }
  }
  out.pairedRegisters = pairs.size();
  out.unpairedRegisters = model1.stateBits.size() - pairs.size();
  for (size_t index = 0; index < model1.stateBits.size(); ++index) {
    if (used1[index]) {
      continue;
    }
    const auto& key = model1.stateBits[index];
    std::string name;
    if (const auto it = model1.displayNameByKey.find(key); it != model1.displayNameByKey.end()) {
      name = it->second;
    }
    if (const auto it = model1.outputNetNameByKey.find(key); it != model1.outputNetNameByKey.end()) {
      name += " net=" + it->second;
    }
    out.unpairedNames.push_back(name);
  }

  // Private rail symbols for design 1.
  std::vector<DualRailSymbolPair> railsByLocalVar(maxVar0 + 1, {kNoState, kNoState});
  std::vector<size_t> inputSymbolByLocalVar(maxVar0 + 1, kNoState);
  auto& problem = out.problem;
  problem.usesDualRailStateEncoding = true;
  problem.initialCondition = BoolExpr::createTrue();
  for (size_t index = 0; index < inputCount; ++index) {
    const size_t symbol = 2 + index;
    problem.inputSymbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    if (const auto it = model0.inputVarByKey.find(alignedInputs.keys0[index]);
        it != model0.inputVarByKey.end()) {
      inputSymbolByLocalVar[it->second] = symbol;
    }
    if (resetActiveValue[index].has_value()) {
      problem.resetBootstrapInputs.emplace_back(symbol, *resetActiveValue[index]);
    }
  }
  std::unordered_map<SignalKey, size_t, SignalKeyHash> indexByKey0;
  indexByKey0.reserve(model0.stateBits.size());
  size_t nextSymbol = 2 + inputCount;
  for (size_t index = 0; index < model0.stateBits.size(); ++index) {
    const auto& key = model0.stateBits[index];
    indexByKey0.emplace(key, index);
    const DualRailSymbolPair rails{nextSymbol, nextSymbol + 1};
    nextSymbol += 2;
    railsByLocalVar[model0.inputVarByKey.at(key)] = rails;
    out.railsByStateIndex.push_back(rails);
    problem.state0Symbols.push_back(rails.mayBeOne);
    problem.state0Symbols.push_back(rails.mayBeZero);
    problem.allSymbols.push_back(rails.mayBeOne);
    problem.allSymbols.push_back(rails.mayBeZero);
    problem.dualRailStatePairs.push_back(rails);
    const auto bits = static_cast<uint8_t>(state0[index]);
    out.frontier[rails.mayBeOne] = (bits & 1u) != 0;
    out.frontier[rails.mayBeZero] = (bits & 2u) != 0;
  }
  problem.totalStateCount = problem.state0Symbols.size();
  const auto railsOfKey0 = [&](const SignalKey& key) -> std::optional<DualRailSymbolPair> {
    const auto it = indexByKey0.find(key);
    if (it == indexByKey0.end()) {
      return std::nullopt;
    }
    return railsByLocalVar[model0.inputVarByKey.at(key)];
  };

  LocalRailMapper mapper(railsByLocalVar, inputSymbolByLocalVar);
  std::unordered_map<BoolExpr*, DualRailBoolExpr> memo;
  const auto lift = [&](BoolExpr* expr) { return buildDualRailBoolExpr(expr, mapper, memo); };
  for (size_t index = 0; index < model0.stateBits.size(); ++index) {
    const auto& key = model0.stateBits[index];
    const auto nextIt = model0.nextStateExprByStateKey.find(key);
    if (nextIt == model0.nextStateExprByStateKey.end() || nextIt->second == nullptr) {
      continue;
    }
    const auto rails = *railsOfKey0(key);
    const auto next = lift(nextIt->second);
    problem.transitions0.emplace_back(rails.mayBeOne, next.mayBeOne);
    problem.transitions0.emplace_back(rails.mayBeZero, next.mayBeZero);
  }

  auto memories = std::make_shared<std::vector<ProofMemoryDescription>>();
  const auto liftWord = [&](const std::vector<BoolExpr*>& exprs) {
    std::vector<DualRailBoolExpr> lifted;
    lifted.reserve(exprs.size());
    for (auto* expr : exprs) {
      lifted.push_back(lift(expr));
    }
    return lifted;
  };
  for (const auto& memory : model0.structuredMemories) {
    ProofMemoryDescription description;
    description.width = memory.width;
    description.depth = memory.depth;
    description.abits = memory.abits;
    bool complete = true;
    for (const auto& cellKeys : memory.cellKeys) {
      std::vector<DualRailSymbolPair> cellRails;
      for (const auto& cellKey : cellKeys) {
        const auto rails = railsOfKey0(cellKey);
        if (!rails.has_value()) {
          complete = false;  // LCOV_EXCL_LINE
          break;  // LCOV_EXCL_LINE
        }
        cellRails.push_back(*rails);
      }
      description.cellRails.push_back(std::move(cellRails));
    }
    if (!complete) {
      continue;  // LCOV_EXCL_LINE
    }
    for (const auto& port : memory.readPorts) {
      ProofMemoryReadPort readPort;
      readPort.address = liftWord(port.addressExprs);
      for (const auto& outputKey : port.outputKeys) {
        readPort.outputRails.push_back(railsOfKey0(outputKey));
      }
      description.readPorts.push_back(std::move(readPort));
    }
    for (const auto& port : memory.writePorts) {
      description.writePorts.push_back(
          {liftWord(port.addressExprs), liftWord(port.dataExprs),
           liftWord(port.maskExprs), liftWord(port.enableExprs)});
    }
    if (memory.resetAssertedExpr != nullptr) {
      description.resetAsserted = lift(memory.resetAssertedExpr);
    }
    memories->push_back(std::move(description));
  }
  if (!memories->empty()) {
    problem.structuredMemories = std::move(memories);
  }

  // Divergence predicate over the differing register pairs and outputs.
  std::unordered_map<BoolExpr*, BoolExpr*> renameMemo;
  const auto readsUnpaired = [&](BoolExpr* renamed) {
    for (const auto var : renamed->getSupportVars()) {
      if (var >= offset) {
        return true;
      }
    }
    return false;
  };
  BoolExpr* diverge = BoolExpr::createFalse();
  // Lockstep descent of two hash-consed formulas: matching nodes are skipped,
  // matching operators recurse, and the first mismatch is an origin pair.
  BoolExpr* origin = BoolExpr::createFalse();
  std::unordered_set<BoolExpr*> originSeen;
  struct PairHash {
    size_t operator()(const std::pair<BoolExpr*, BoolExpr*>& pair) const noexcept {
      return std::hash<BoolExpr*>()(pair.first) * 31u ^ std::hash<BoolExpr*>()(pair.second);
    }
  };
  std::unordered_set<std::pair<BoolExpr*, BoolExpr*>, PairHash> visitedPairs;
  constexpr size_t kMaxOriginPairs = 4096;
  const auto collectOrigins = [&](BoolExpr* left, BoolExpr* right) {
    std::vector<std::pair<BoolExpr*, BoolExpr*>> stack{{left, right}};
    while (!stack.empty()) {
      const auto [u, v] = stack.back();
      stack.pop_back();
      if (u == v || u == nullptr || v == nullptr || !visitedPairs.insert({u, v}).second) {
        continue;
      }
      const bool sameShape = u->getOp() == v->getOp() && u->getOp() != Op::VAR &&
                             (u->getLeft() == nullptr) == (v->getLeft() == nullptr) &&
                             (u->getRight() == nullptr) == (v->getRight() == nullptr);
      if (sameShape) {
        // Descend only where the operand correspondence is unambiguous: one
        // operand identical on both sides (either order, every binary
        // operator here is commutative), the other pair carries the
        // difference. When both operands differ the pairing could be wrong,
        // so this node pair is the origin.
        BoolExpr* const ul = u->getLeft();
        BoolExpr* const ur = u->getRight();
        BoolExpr* const vl = v->getLeft();
        BoolExpr* const vr = v->getRight();
        if (u->getOp() == Op::NOT || ur == nullptr) {
          stack.emplace_back(ul, vl);
          continue;
        }
        if (ul == vl) {
          stack.emplace_back(ur, vr);
          continue;
        }
        if (ur == vr) {
          stack.emplace_back(ul, vl);
          continue;
        }
        if (ul == vr) {
          stack.emplace_back(ur, vl);
          continue;
        }
        if (ur == vl) {
          stack.emplace_back(ul, vr);
          continue;
        }
      }
      if (out.originPairs < kMaxOriginPairs && originSeen.insert(u).second) {
        ++out.originPairs;
        origin = BoolExpr::Or(origin, binaryDifferent(lift(u), lift(v)));
      }
    }
  };
  for (const auto& [index0, index1] : pairs) {
    const auto next0It = model0.nextStateExprByStateKey.find(model0.stateBits[index0]);
    const auto next1It = model1.nextStateExprByStateKey.find(model1.stateBits[index1]);
    if (next0It == model0.nextStateExprByStateKey.end() ||
        next1It == model1.nextStateExprByStateKey.end() ||
        next0It->second == nullptr || next1It->second == nullptr) {
      continue;
    }
    BoolExpr* renamed = renameVariables(next1It->second, rename, renameMemo);
    if (renamed == next0It->second) {
      continue;
    }
    if (readsUnpaired(renamed)) {
      ++out.incomparableRegisters;
      if (out.incomparableNames.size() < 12) {
        const auto it = model0.displayNameByKey.find(model0.stateBits[index0]);
        out.incomparableNames.push_back(
            it == model0.displayNameByKey.end() ? std::string("?") : it->second);
      }
      continue;
    }
    ++out.differingRegisters;
    out.differingIndices.push_back(index0);
    diverge = BoolExpr::Or(diverge, binaryDifferent(lift(next0It->second), lift(renamed)));
    collectOrigins(next0It->second, renamed);
  }
  out.registerDivergence = BoolExpr::simplify(diverge);
  out.originDivergence = BoolExpr::simplify(origin);
  for (size_t index = 0; index < alignedOutputs.names.size(); ++index) {
    const auto out0It = model0.observedOutputExprByKey.find(alignedOutputs.keys0[index]);
    const auto out1It = model1.observedOutputExprByKey.find(alignedOutputs.keys1[index]);
    if (out0It == model0.observedOutputExprByKey.end() ||
        out1It == model1.observedOutputExprByKey.end() ||
        out0It->second == nullptr || out1It->second == nullptr) {
      continue;  // LCOV_EXCL_LINE
    }
    BoolExpr* renamed = renameVariables(out1It->second, rename, renameMemo);
    if (renamed == out0It->second || readsUnpaired(renamed)) {
      continue;
    }
    ++out.differingOutputs;
    diverge = BoolExpr::Or(diverge, binaryDifferent(lift(out0It->second), lift(renamed)));
  }
  problem.bad = BoolExpr::simplify(diverge);
  problem.property = BoolExpr::simplify(BoolExpr::Not(problem.bad));
  problem.description = "seeded search divergence predicate";
  return out;
}

// Per-frame state-domain constraints shared by every unrolling here: valid
// dual-rail values, same-frame equalities and complements, released reset
// controls (the simulated prefix is at least as long as the reset bootstrap).
void addFrameDomainConstraints(
    SATSolverWrapper& solver,
    const FrameVariableStore& variables,
    const KInductionProblem& problem,
    size_t frame) {
  for (const auto& pair : problem.dualRailStatePairs) {
    if (variables.hasSymbol(pair.mayBeOne) && variables.hasSymbol(pair.mayBeZero)) {
      solver.addClause({variables.getLiteral(pair.mayBeOne, frame),
                        variables.getLiteral(pair.mayBeZero, frame)});
    }
  }
  for (const auto* pairs : {&problem.sameFrameStateEqualityPairs0,
                            &problem.sameFrameStateEqualityPairs1}) {
    for (const auto& [lhs, rhs] : *pairs) {
      if (variables.hasSymbol(lhs) && variables.hasSymbol(rhs)) {
        addLiteralEquivalence(
            solver, variables.getLiteral(lhs, frame), variables.getLiteral(rhs, frame));
      }
    }
  }
  for (const auto* pairs : {&problem.complementedStatePairs0,
                            &problem.complementedStatePairs1}) {
    for (const auto& [primary, complemented] : *pairs) {
      if (variables.hasSymbol(primary) && variables.hasSymbol(complemented)) {
        addLiteralEquivalence(
            solver, variables.getLiteral(primary, frame),
            -variables.getLiteral(complemented, frame));
      }
    }
  }
  for (const auto& [symbol, activeValue] : problem.resetBootstrapInputs) {
    if (variables.hasSymbol(symbol)) {
      const int literal = variables.getLiteral(symbol, frame);
      solver.addClause({activeValue ? -literal : literal});
    }
  }
}

int cubeLiteralAtFrame(
    RailCnfBuilder& rails,
    const FrameVariableStore& variables,
    const StateCube& cube,
    size_t frame) {
  int literal = rails.trueLit();
  for (const auto& [symbol, value] : cube) {
    if (!variables.hasSymbol(symbol)) {
      continue;  // LCOV_EXCL_LINE
    }
    const int railLiteral = variables.getLiteral(symbol, frame);
    literal = rails.andLit(literal, value ? railLiteral : -railLiteral);
  }
  return literal;
}

// Backward cone of a target over the problem's state symbols: layer 0 is the
// target's own state support, layer j+1 adds every register a layer-j
// register's next-state function reads; a memory read output pulls in the
// registers its memory's ports read. Memory cells are not state here.
std::vector<std::unordered_set<size_t>> targetConeLayers(
    const KInductionProblem& problem,
    BoolExpr* target,
    size_t maxLayers) {
  const std::vector<ProofMemoryDescription> noMemories;
  const auto& memories = problem.structuredMemories != nullptr
                             ? *problem.structuredMemories
                             : noMemories;
  std::unordered_set<size_t> cellRailSymbols;
  std::unordered_map<size_t, size_t> memoryByReadOutput;
  for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
    for (const auto& cell : memories[memoryIndex].cellRails) {
      for (const auto& rails : cell) {
        cellRailSymbols.insert(rails.mayBeOne);
        cellRailSymbols.insert(rails.mayBeZero);
      }
    }
    for (const auto& port : memories[memoryIndex].readPorts) {
      for (const auto& rails : port.outputRails) {
        if (rails.has_value()) {
          memoryByReadOutput.emplace(rails->mayBeOne, memoryIndex);
          memoryByReadOutput.emplace(rails->mayBeZero, memoryIndex);
        }
      }
    }
  }
  std::unordered_set<size_t> stateSymbols;
  for (const auto symbol : problem.combinedStateSymbols()) {
    if (cellRailSymbols.find(symbol) == cellRailSymbols.end()) {
      stateSymbols.insert(symbol);
    }
  }
  const auto addStateSupport = [&](const std::set<size_t>& support,
                                   std::unordered_set<size_t>& out) {
    for (const auto symbol : support) {
      if (stateSymbols.find(symbol) != stateSymbols.end()) {
        out.insert(symbol);
      }
    }
  };
  std::vector<std::unordered_set<size_t>> memoryPortSupport(memories.size());
  for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
    const auto& memory = memories[memoryIndex];
    const auto addRails = [&](const DualRailBoolExpr& expr) {
      addStateSupport(expr.mayBeOne->getSupportVars(), memoryPortSupport[memoryIndex]);
      addStateSupport(expr.mayBeZero->getSupportVars(), memoryPortSupport[memoryIndex]);
    };
    for (const auto& port : memory.readPorts) {
      for (const auto& bit : port.address) {
        addRails(bit);
      }
    }
    for (const auto& port : memory.writePorts) {
      for (const auto* word : {&port.address, &port.data, &port.mask, &port.enable}) {
        for (const auto& bit : *word) {
          addRails(bit);
        }
      }
    }
    if (memory.resetAsserted.has_value()) {
      addRails(*memory.resetAsserted);
    }
  }
  const TransitionExprResolver transitions(problem);
  std::vector<std::unordered_set<size_t>> layers(1);
  addStateSupport(target->getSupportVars(), layers[0]);
  for (size_t layer = 0; layer < maxLayers; ++layer) {
    std::unordered_set<size_t> next = layers[layer];
    for (const auto symbol : layers[layer]) {
      if (const auto readIt = memoryByReadOutput.find(symbol);
          readIt != memoryByReadOutput.end()) {
        const auto& support = memoryPortSupport[readIt->second];
        next.insert(support.begin(), support.end());
      } else if (transitions.contains(symbol)) {
        addStateSupport(transitions.support(symbol), next);
      }
    }
    const bool grew = next.size() > layers[layer].size();
    layers.push_back(std::move(next));
    if (!grew) {
      break;
    }
  }
  return layers;
}

// Optional controls of the bounded search used by the abstraction-guided
// search: a localized model (only kept states have transitions; the others are
// free after frame 0), pinned rail values as solver assumptions with the
// failed subset returned on UNSAT, and per-frame values of chosen symbols.
struct BoundedSearchExtras {
  const std::unordered_set<size_t>* keptStates = nullptr;
  const std::vector<std::tuple<size_t, size_t, bool>>* pins = nullptr;  // (symbol, frame, value)
  std::vector<std::tuple<size_t, size_t, bool>>* failedPinsOut = nullptr;
  const std::vector<size_t>* traceSymbols = nullptr;
  std::vector<std::vector<std::pair<size_t, bool>>>* traceValuesOut = nullptr;  // [frame][(symbol, value)]
  bool forceAssumptionSolver = false;
  bool* unsatOut = nullptr;
};

// Bounded search from one concrete rail-valued state. Memory arrays are not
// unrolled cell by cell: a read returns the data of the latest earlier write
// to the same address and otherwise the frontier contents (read-over-write,
// Burch and Dill, CAV 1994; Velev and Bryant, TACAS 1998 for the ternary
// form).
std::optional<FrontierBoundedSearchResult> runFrontierBoundedSearch(
    const KInductionProblem& problem,
    const std::unordered_map<size_t, bool>& frontierBySymbol,
    const std::unordered_map<size_t, size_t>& mergedRails,
    size_t depth,
    size_t checkEvery,
    const TargetEnlargement& enlargement,
    BoolExpr* intermediateTarget,
    KEPLER_FORMAL::Config::SolverType solverType,
    bool diagnostics,
    std::vector<std::pair<size_t, bool>>* frameZeroStateOut = nullptr,
    bool targetAtLastFrameOnly = false,
    size_t conflictLimit = 0,
    const std::vector<std::pair<size_t, bool>>& heldInputs = {},
    const BoundedSearchExtras& extras = {}) {
  const auto start = std::chrono::steady_clock::now();
  const bool finalCheckOnly = checkEvery == 0 || targetAtLastFrameOnly;
  if (extras.unsatOut != nullptr) {
    *extras.unsatOut = false;
  }
  std::unordered_map<BoolExpr*, BoolExpr*> mergeMemo;
  const auto mergeExpr = [&](BoolExpr* expr) {
    return mergedRails.empty() ? expr : renameVariables(expr, mergedRails, mergeMemo);
  };
  const auto mergeRails = [&](const DualRailBoolExpr& expr) {
    return DualRailBoolExpr{mergeExpr(expr.mayBeOne), mergeExpr(expr.mayBeZero)};
  };
  const auto isMerged = [&](size_t symbol) {
    return mergedRails.find(symbol) != mergedRails.end();
  };
  const std::vector<ProofMemoryDescription> noMemories;
  const auto& memories = problem.structuredMemories != nullptr
                             ? *problem.structuredMemories
                             : noMemories;

  std::unordered_set<size_t> cellRailSymbols;
  std::unordered_map<size_t, ReadOutputRef> readOutputBySymbol;
  for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
    const auto& memory = memories[memoryIndex];
    for (const auto& cell : memory.cellRails) {
      for (const auto& rails : cell) {
        cellRailSymbols.insert(rails.mayBeOne);
        cellRailSymbols.insert(rails.mayBeZero);
      }
    }
    for (size_t portIndex = 0; portIndex < memory.readPorts.size(); ++portIndex) {
      const auto& port = memory.readPorts[portIndex];
      for (size_t bit = 0; bit < port.outputRails.size(); ++bit) {
        if (!port.outputRails[bit].has_value()) {
          continue;
        }
        readOutputBySymbol.emplace(
            port.outputRails[bit]->mayBeOne, ReadOutputRef{memoryIndex, portIndex, bit});
        readOutputBySymbol.emplace(
            port.outputRails[bit]->mayBeZero, ReadOutputRef{memoryIndex, portIndex, bit});
      }
    }
  }

  std::vector<size_t> stateSymbols;
  for (const auto symbol : problem.combinedStateSymbols()) {
    if (cellRailSymbols.find(symbol) == cellRailSymbols.end() && !isMerged(symbol)) {
      stateSymbols.push_back(symbol);
    }
  }
  const std::unordered_set<size_t> stateSymbolSet(stateSymbols.begin(), stateSymbols.end());
  const TransitionExprResolver transitions(problem);
  std::vector<size_t> transitionStates;
  for (const auto symbol : stateSymbols) {
    if (readOutputBySymbol.find(symbol) != readOutputBySymbol.end() ||
        !transitions.contains(symbol)) {
      continue;
    }
    if (supportTouches(transitions.support(symbol), cellRailSymbols)) {
      logSeeded(diagnostics,
                "state v%zu reads memory cells outside a read port; search skipped",
                symbol);
      return std::nullopt;
    }
    transitionStates.push_back(symbol);
  }
  if (problem.bad == nullptr) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  {
    const auto badSupport = problem.bad->getSupportVars();
    if (supportTouches(badSupport, cellRailSymbols)) {
      logSeeded(diagnostics, "property reads memory cells directly; search skipped");
      return std::nullopt;
    }
  }
  // With design 2's shared registers renamed, outputs with identical cones
  // fold away and the property keeps only the logic that differs.
  BoolExpr* bad = BoolExpr::simplify(mergeExpr(problem.bad));
  if (bad == BoolExpr::createFalse()) {
    logSeeded(diagnostics, "every output cone is shared after pairing; no difference possible");
    return std::nullopt;
  }

  std::vector<size_t> symbols = problem.inputSymbols;
  symbols.insert(symbols.end(), stateSymbols.begin(), stateSymbols.end());

  // State-set target enlargement: generalized preimage cubes of the target.
  // Memory cells are not search variables, so a problem with structured
  // memories keeps the structural form only.
  std::vector<StateCube> cubes;
  if (enlargement.cubeSteps > 0 && enlargement.cubeBudget > 0) {
    if (!memories.empty()) {
      logSeeded(diagnostics,
                "state-set enlargement skipped: structured memories keep the structural form");
    } else {
      const auto enlargeStart = std::chrono::steady_clock::now();
      std::vector<StateCube> levelTargets;  // cubes of the previous level
      for (size_t level = 1; level <= enlargement.cubeSteps; ++level) {
        const auto levelStart = std::chrono::steady_clock::now();
        const auto enumSolverType = SATSolverWrapper::assumptionSolverTypeFor(solverType);
        SATSolverWrapper enumSolver(enumSolverType);
        SATSolverWrapper liftSolver(enumSolverType);
        enumSolver.configureForSecIncrementalBoundedSearch();
        liftSolver.configureForSecIncrementalBoundedSearch();
        FrameVariableStore enumVariables(enumSolver, symbols, 2);
        FrameVariableStore liftVariables(liftSolver, symbols, 2);
        RailCnfBuilder enumRails(enumSolver);
        RailCnfBuilder liftRails(liftSolver);
        // Both instances: one free step; the target at frame 1 is the original
        // target or one of the previous level's cubes.
        const auto build = [&](SATSolverWrapper& instance,
                               const FrameVariableStore& instanceVariables,
                               RailCnfBuilder& instanceRails) {
          addFrameDomainConstraints(instance, instanceVariables, problem, 0);
          addFrameDomainConstraints(instance, instanceVariables, problem, 1);
          FrameFormulaEncoder stepEncoder(instance, instanceVariables.makeLeafLits(0));
          for (const auto symbol : transitionStates) {
            addLiteralEquivalence(
                instance,
                instanceVariables.getLiteral(symbol, 1),
                stepEncoder.encode(mergeExpr(transitions.at(symbol))));
          }
          FrameFormulaEncoder targetEncoder(instance, instanceVariables.makeLeafLits(1));
          int target = targetEncoder.encode(bad);
          for (const auto& cube : levelTargets) {
            target = instanceRails.orLit(
                target, cubeLiteralAtFrame(instanceRails, instanceVariables, cube, 1));
          }
          return target;
        };
        const int enumTarget = build(enumSolver, enumVariables, enumRails);
        const int liftTarget = build(liftSolver, liftVariables, liftRails);
        liftSolver.addClause({-liftTarget});
        for (const auto& cube : cubes) {
          std::vector<int> blocking;
          for (const auto& [symbol, value] : cube) {
            const int literal = enumVariables.getLiteral(symbol, 0);
            blocking.push_back(value ? -literal : literal);
          }
          enumSolver.addClause(blocking);
        }
        std::vector<StateCube> levelCubes;
        bool complete = false;
        size_t literalTotal = 0;
        while (cubes.size() + levelCubes.size() < enlargement.cubeBudget) {
          const auto status = enumSolver.solveWithAssumptionsStatus({enumTarget});
          if (status == SATSolverWrapper::SolveStatus::Unsat) {
            complete = true;
            break;
          }
          if (status != SATSolverWrapper::SolveStatus::Sat) {
            break;  // LCOV_EXCL_LINE
          }
          // The model's state is one preimage point; its inputs are fixed and
          // the state literals not needed to force the target are dropped.
          StateCube cube;
          std::vector<int> assumptions;
          cube.reserve(stateSymbols.size());
          for (const auto symbol : stateSymbols) {
            const bool value = enumSolver.getLiteralValue(enumVariables.getLiteral(symbol, 0));
            cube.emplace_back(symbol, value);
            const int literal = liftVariables.getLiteral(symbol, 0);
            assumptions.push_back(value ? literal : -literal);
          }
          for (size_t frame = 0; frame < 2; ++frame) {
            for (const auto symbol : problem.inputSymbols) {
              const bool value =
                  enumSolver.getLiteralValue(enumVariables.getLiteral(symbol, frame));
              const int literal = liftVariables.getLiteral(symbol, frame);
              assumptions.push_back(value ? literal : -literal);
            }
          }
          StateCube generalized = cube;
          if (liftSolver.solveWithAssumptionsStatus(assumptions) ==
              SATSolverWrapper::SolveStatus::Unsat) {
            const auto failedList = liftSolver.failedAssumptions();
            const std::unordered_set<int> failed(failedList.begin(), failedList.end());
            generalized.clear();
            for (const auto& [symbol, value] : cube) {
              const int literal = liftVariables.getLiteral(symbol, 0);
              if (failed.count(value ? literal : -literal) != 0) {
                generalized.emplace_back(symbol, value);
              }
            }
          }
          literalTotal += generalized.size();
          std::vector<int> blocking;
          for (const auto& [symbol, value] : generalized) {
            const int literal = enumVariables.getLiteral(symbol, 0);
            blocking.push_back(value ? -literal : literal);
          }
          enumSolver.addClause(blocking);
          levelCubes.push_back(std::move(generalized));
        }
        logSeeded(diagnostics,
                  "enlargement level=%zu cubes=%zu complete=%d avg_literals=%.1f seconds=%.1f",
                  level, levelCubes.size(), complete ? 1 : 0,
                  levelCubes.empty() ? 0.0
                                     : static_cast<double>(literalTotal) /
                                           static_cast<double>(levelCubes.size()),
                  elapsedSeconds(levelStart));
        if (levelCubes.empty()) {
          break;
        }
        cubes.insert(cubes.end(), levelCubes.begin(), levelCubes.end());
        levelTargets = std::move(levelCubes);
        if (!complete) {
          break;  // a partial level is still sound, but its preimage would be too
        }
      }
      logSeeded(diagnostics, "enlargement cubes=%zu (%.1f s)", cubes.size(),
                elapsedSeconds(enlargeStart));
    }
  }

  // Structural enlargement: backward cone layers of the (enlarged) target
  // over the state symbols; tail frames unroll only what the target needs.
  const size_t tailFrames = enlargement.structuralSteps;
  const size_t totalFrames = depth + tailFrames;
  std::vector<std::unordered_set<size_t>> layers(tailFrames + 1);
  if (tailFrames > 0) {
    const auto addMergedStateSupport = [&](const std::set<size_t>& support,
                                           std::unordered_set<size_t>& out) {
      for (const auto raw : support) {
        const auto mergedIt = mergedRails.find(raw);
        const size_t symbol = mergedIt == mergedRails.end() ? raw : mergedIt->second;
        if (stateSymbolSet.find(symbol) != stateSymbolSet.end()) {
          out.insert(symbol);
        }
      }
    };
    std::vector<std::unordered_set<size_t>> memoryPortSupport(memories.size());
    for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
      const auto& memory = memories[memoryIndex];
      const auto addRails = [&](const DualRailBoolExpr& expr) {
        addMergedStateSupport(expr.mayBeOne->getSupportVars(), memoryPortSupport[memoryIndex]);
        addMergedStateSupport(expr.mayBeZero->getSupportVars(), memoryPortSupport[memoryIndex]);
      };
      for (const auto& port : memory.readPorts) {
        for (const auto& bit : port.address) {
          addRails(bit);
        }
      }
      for (const auto& port : memory.writePorts) {
        for (const auto* word : {&port.address, &port.data, &port.mask, &port.enable}) {
          for (const auto& bit : *word) {
            addRails(bit);
          }
        }
      }
      if (memory.resetAsserted.has_value()) {
        addRails(*memory.resetAsserted);
      }
    }
    addMergedStateSupport(bad->getSupportVars(), layers[0]);
    if (intermediateTarget != nullptr) {
      addMergedStateSupport(intermediateTarget->getSupportVars(), layers[0]);
    }
    for (const auto& cube : cubes) {
      for (const auto& [symbol, value] : cube) {
        layers[0].insert(symbol);
      }
    }
    for (size_t layer = 0; layer < tailFrames; ++layer) {
      layers[layer + 1] = layers[layer];
      for (const auto symbol : layers[layer]) {
        if (const auto readIt = readOutputBySymbol.find(symbol);
            readIt != readOutputBySymbol.end()) {
          const auto& support = memoryPortSupport[readIt->second.memory];
          layers[layer + 1].insert(support.begin(), support.end());
        } else if (transitions.contains(symbol)) {
          addMergedStateSupport(transitions.support(symbol), layers[layer + 1]);
        }
      }
    }
    logSeeded(diagnostics, "structural enlargement tail_frames=%zu cone_states=%zu of %zu",
              tailFrames, layers[tailFrames].size(), stateSymbols.size());
  }

  // One final check is a single query, so the configured solver serves it;
  // intermediate checks need assumptions and an incremental solver.
  const auto searchSolverType = finalCheckOnly && !extras.forceAssumptionSolver
      ? solverType
      : SATSolverWrapper::assumptionSolverTypeFor(solverType);
  const bool usesAssumptions =
      searchSolverType != KEPLER_FORMAL::Config::SolverType::KISSAT;
  SATSolverWrapper solver(searchSolverType);
  solver.configureForSecIncrementalBoundedSearch();
  FrameVariableStore variables(solver, symbols, totalFrames + 1);
  RailCnfBuilder rails(solver);

  // Frame 0 is the simulated frontier.
  size_t frontierUnits = 0;
  for (const auto symbol : stateSymbols) {
    if (const auto it = frontierBySymbol.find(symbol); it != frontierBySymbol.end()) {
      const int literal = variables.getLiteral(symbol, 0);
      solver.addClause({it->second ? literal : -literal});
      ++frontierUnits;
    }
  }
  const auto holdInputs = [&](size_t frame) {
    for (const auto& [symbol, value] : heldInputs) {
      if (variables.hasSymbol(symbol)) {
        const int literal = variables.getLiteral(symbol, frame);
        solver.addClause({value ? literal : -literal});
      }
    }
  };
  addFrameDomainConstraints(solver, variables, problem, 0);
  holdInputs(0);

  const auto extractInputs = [&](size_t badFrame, bool viaIntermediateTarget) {
    FrontierBoundedSearchResult result;
    result.badFrame = badFrame;
    result.viaIntermediateTarget = viaIntermediateTarget;
    if (frameZeroStateOut != nullptr) {
      frameZeroStateOut->clear();
      frameZeroStateOut->reserve(stateSymbols.size());
      for (const auto symbol : stateSymbols) {
        frameZeroStateOut->emplace_back(
            symbol, solver.getLiteralValue(variables.getLiteral(symbol, 0)));
      }
    }
    if (extras.traceSymbols != nullptr && extras.traceValuesOut != nullptr) {
      extras.traceValuesOut->assign(badFrame + 1, {});
      for (size_t traceFrame = 0; traceFrame <= badFrame; ++traceFrame) {
        auto& values = (*extras.traceValuesOut)[traceFrame];
        values.reserve(extras.traceSymbols->size());
        for (const auto symbol : *extras.traceSymbols) {
          if (variables.hasSymbol(symbol)) {
            values.emplace_back(
                symbol, solver.getLiteralValue(variables.getLiteral(symbol, traceFrame)));
          }
        }
      }
    }
    result.inputsByFrame.resize(badFrame + 1);
    for (size_t traceFrame = 0; traceFrame <= badFrame; ++traceFrame) {
      auto& inputs = result.inputsByFrame[traceFrame];
      inputs.reserve(problem.inputSymbols.size());
      for (const auto symbol : problem.inputSymbols) {
        inputs.push_back(solver.getLiteralValue(variables.getLiteral(symbol, traceFrame)));
      }
    }
    return result;
  };
  const auto encodeTarget = [&](size_t frame, int& badLit, int& targetLit) {
    FrameFormulaEncoder targetEncoder(solver, variables.makeLeafLits(frame));
    badLit = targetEncoder.encode(bad);
    targetLit = badLit;
    if (intermediateTarget != nullptr) {
      targetLit = rails.orLit(targetLit, targetEncoder.encode(mergeExpr(intermediateTarget)));
    }
    for (const auto& cube : cubes) {
      targetLit = rails.orLit(targetLit, cubeLiteralAtFrame(rails, variables, cube, frame));
    }
  };

  std::vector<int> pendingBadLits;
  std::vector<int> pendingTargetLits;
  size_t firstPendingFrame = 0;
  // The selected state itself may already show the difference.
  {
    int badLit = 0;
    int targetLit = 0;
    encodeTarget(0, badLit, targetLit);
    if (targetAtLastFrameOnly) {
      firstPendingFrame = totalFrames;  // only the last frame's target counts
    } else if (finalCheckOnly) {
      pendingBadLits.push_back(badLit);
      pendingTargetLits.push_back(targetLit);
    } else {
      const auto status = solver.solveWithAssumptionsStatus({targetLit});
      if (status == SATSolverWrapper::SolveStatus::Sat) {
        logSeeded(diagnostics, "difference at the selected state itself");
        return extractInputs(0, !solver.getLiteralValue(badLit));
      }
      if (status != SATSolverWrapper::SolveStatus::Unsat) {
        return std::nullopt;  // LCOV_EXCL_LINE
      }
      solver.addClause({-targetLit});
      firstPendingFrame = 1;
    }
  }

  // Frontier memory contents per cell bit, as rail values.
  struct FrontierCell {
    bool one = true;
    bool zero = true;
  };
  const auto frontierRail = [&](size_t symbol) {
    const auto it = frontierBySymbol.find(symbol);
    return it == frontierBySymbol.end() ? true : it->second;
  };

  std::vector<std::vector<MemoryFrameEvents>> eventsByMemory(memories.size());
  for (auto& events : eventsByMemory) {
    events.reserve(totalFrames);
  }
  const auto mergeRailWord = [&](const std::vector<DualRailBoolExpr>& word) {
    std::vector<DualRailBoolExpr> merged;
    merged.reserve(word.size());
    for (const auto& bit : word) {
      merged.push_back(mergeRails(bit));
    }
    return merged;
  };
  // A memory whose every read output is paired with design 1's is design 1's
  // memory; its chain would only duplicate work.
  std::vector<bool> sharedMemory(memories.size(), false);
  for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
    bool allMerged = true;
    bool anyOutput = false;
    for (const auto& port : memories[memoryIndex].readPorts) {
      for (const auto& rails : port.outputRails) {
        if (rails.has_value()) {
          anyOutput = true;
          allMerged = allMerged && isMerged(rails->mayBeOne);
        }
      }
    }
    sharedMemory[memoryIndex] = anyOutput && allMerged;
  }

  // Pinned rail values: assumptions when the solver supports them (so the
  // failed subset can be read back), otherwise unit clauses.
  std::vector<int> pinLiterals;
  std::unordered_map<int, std::tuple<size_t, size_t, bool>> pinByLiteral;
  if (extras.pins != nullptr) {
    for (const auto& pin : *extras.pins) {
      const auto [symbol, frame, value] = pin;
      if (frame > totalFrames || !variables.hasSymbol(symbol)) {
        continue;
      }
      const int literal = variables.getLiteral(symbol, frame);
      const int signedLiteral = value ? literal : -literal;
      if (usesAssumptions) {
        pinLiterals.push_back(signedLiteral);
        pinByLiteral.emplace(signedLiteral, pin);
      } else {
        solver.addClause({signedLiteral});
      }
    }
  }
  const auto checkPendingTarget = [&](size_t frame) -> std::optional<FrontierBoundedSearchResult> {
    int anyTarget = rails.falseLit();
    for (const int lit : pendingTargetLits) {
      anyTarget = rails.orLit(anyTarget, lit);
    }
    const auto solveStart = std::chrono::steady_clock::now();
    SATSolverWrapper::SolveStatus status = SATSolverWrapper::SolveStatus::Unknown;
    if (usesAssumptions) {
      std::vector<int> assumptions = pinLiterals;
      assumptions.push_back(anyTarget);
      status = solver.solveWithAssumptionsStatus(
          assumptions, conflictLimit == 0 ? -1 : static_cast<int64_t>(conflictLimit));
      if (status == SATSolverWrapper::SolveStatus::Unsat && extras.failedPinsOut != nullptr) {
        extras.failedPinsOut->clear();
        for (const int failed : solver.failedAssumptions()) {
          if (const auto it = pinByLiteral.find(failed); it != pinByLiteral.end()) {
            extras.failedPinsOut->push_back(it->second);
          }
        }
      }
    } else {
      solver.addClause({anyTarget});
      status = conflictLimit == 0
          ? solver.solveStatus()
          : solver.solveWithKissatResourceLimits(static_cast<unsigned>(conflictLimit));
    }
    if (status == SATSolverWrapper::SolveStatus::Unsat && extras.unsatOut != nullptr) {
      *extras.unsatOut = true;
    }
    logSeeded(diagnostics,
              "frame=%zu checked_frames=%zu status=%s solve_seconds=%.1f total_seconds=%.1f",
              frame,
              pendingTargetLits.size(),
              status == SATSolverWrapper::SolveStatus::Sat
                  ? "sat"
                  : status == SATSolverWrapper::SolveStatus::Unsat ? "unsat" : "unknown",
              elapsedSeconds(solveStart),
              elapsedSeconds(start));
    if (status == SATSolverWrapper::SolveStatus::Sat) {
      // The earliest checked frame whose target holds.
      for (size_t i = 0; i < pendingTargetLits.size(); ++i) {
        if (solver.getLiteralValue(pendingTargetLits[i])) {
          return extractInputs(firstPendingFrame + i,
                               !solver.getLiteralValue(pendingBadLits[i]));
        }
      }
      return extractInputs(frame, false);  // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Unsat && usesAssumptions) {
      for (const int lit : pendingTargetLits) {
        solver.addClause({-lit});
      }
    }
    pendingBadLits.clear();
    pendingTargetLits.clear();
    firstPendingFrame = frame + 1;
    return std::nullopt;
  };

  for (size_t frame = 0; frame < totalFrames; ++frame) {
    const auto frameStart = std::chrono::steady_clock::now();
    // Tail frames of the structural enlargement unroll only the states the
    // target still depends on at that distance from the last frame.
    const std::unordered_set<size_t>* needed =
        frame >= depth ? &layers[totalFrames - frame - 1] : nullptr;
    const auto isNeeded = [&](size_t symbol) {
      return (needed == nullptr || needed->find(symbol) != needed->end()) &&
             (extras.keptStates == nullptr ||
              extras.keptStates->find(symbol) != extras.keptStates->end());
    };
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    size_t encodedStates = 0;
    for (const auto symbol : transitionStates) {
      if (!isNeeded(symbol)) {
        continue;
      }
      ++encodedStates;
      addLiteralEquivalence(
          solver,
          variables.getLiteral(symbol, frame + 1),
          encoder.encode(mergeExpr(transitions.at(symbol))));
    }

    for (size_t memoryIndex = 0; memoryIndex < memories.size(); ++memoryIndex) {
      const auto& memory = memories[memoryIndex];
      if (sharedMemory[memoryIndex]) {
        continue;
      }
      bool memoryNeeded = needed == nullptr && extras.keptStates == nullptr;
      for (const auto& port : memory.readPorts) {
        for (const auto& outputRails : port.outputRails) {
          memoryNeeded = memoryNeeded ||
                         (outputRails.has_value() && isNeeded(outputRails->mayBeOne));
        }
      }
      if (!memoryNeeded) {
        continue;
      }
      auto& events = eventsByMemory[memoryIndex];
      MemoryFrameEvents frameEvents;
      for (const auto& port : memory.writePorts) {
        WriteEvent event;
        event.enable = rails.constRail(true);
        for (const auto& enable : port.enable) {
          event.enable = rails.railAnd(event.enable, encodeRails(encoder, mergeRails(enable)));
        }
        event.address = encodeRailWord(encoder, mergeRailWord(port.address));
        event.data = encodeRailWord(encoder, mergeRailWord(port.data));
        event.mask = encodeRailWord(encoder, mergeRailWord(port.mask));
        frameEvents.writes.push_back(std::move(event));
      }
      if (memory.resetAsserted.has_value()) {
        frameEvents.reset = encodeRails(encoder, mergeRails(*memory.resetAsserted));
      }
      events.push_back(std::move(frameEvents));

      for (const auto& port : memory.readPorts) {
        bool portNeeded = false;
        for (const auto& outputRails : port.outputRails) {
          portNeeded = portNeeded ||
                       (outputRails.has_value() && isNeeded(outputRails->mayBeOne));
        }
        if (!portNeeded) {
          continue;
        }
        const auto address = encodeRailWord(encoder, mergeRailWord(port.address));
        // Address decode against the frontier contents, built on demand.
        std::vector<int> mayEqualCell;
        const auto decode = [&]() {
          if (mayEqualCell.empty()) {
            mayEqualCell.reserve(memory.depth);
            for (size_t cell = 0; cell < memory.depth; ++cell) {
              mayEqualCell.push_back(rails.wordMayEqualConstant(address, cell));
            }
          }
        };
        // Write conditions of every earlier and same-cycle write event.
        std::vector<std::vector<RailLits>> writeHits(events.size());
        for (size_t eventFrame = 0; eventFrame < events.size(); ++eventFrame) {
          for (const auto& write : events[eventFrame].writes) {
            writeHits[eventFrame].push_back(
                rails.railAnd(write.enable, rails.wordEqual(write.address, address)));
          }
        }
        for (size_t bit = 0; bit < port.outputRails.size(); ++bit) {
          if (!port.outputRails[bit].has_value() ||
              isMerged(port.outputRails[bit]->mayBeOne) ||
              !variables.hasSymbol(port.outputRails[bit]->mayBeOne) ||
              !isNeeded(port.outputRails[bit]->mayBeOne)) {
            continue;
          }
          // Frontier contents at the read address.
          RailLits value;
          {
            bool allOne = true;
            bool allZero = true;
            bool anyOne = false;
            bool anyZero = false;
            std::vector<FrontierCell> cells(memory.depth);
            for (size_t cell = 0; cell < memory.depth; ++cell) {
              const auto& cellRails = memory.cellRails[cell][bit];
              cells[cell].one = frontierRail(cellRails.mayBeOne);
              cells[cell].zero = frontierRail(cellRails.mayBeZero);
              allOne = allOne && cells[cell].one;
              allZero = allZero && cells[cell].zero;
              anyOne = anyOne || cells[cell].one;
              anyZero = anyZero || cells[cell].zero;
            }
            if (allOne || !anyOne) {
              value.one = allOne ? rails.trueLit() : rails.falseLit();
            } else {
              decode();
              value.one = rails.falseLit();
              for (size_t cell = 0; cell < memory.depth; ++cell) {
                if (cells[cell].one) {
                  value.one = rails.orLit(value.one, mayEqualCell[cell]);
                }
              }
            }
            if (allZero || !anyZero) {
              value.zero = allZero ? rails.trueLit() : rails.falseLit();
            } else {
              decode();
              value.zero = rails.falseLit();
              for (size_t cell = 0; cell < memory.depth; ++cell) {
                if (cells[cell].zero) {
                  value.zero = rails.orLit(value.zero, mayEqualCell[cell]);
                }
              }
            }
          }
          for (size_t eventFrame = 0; eventFrame < events.size(); ++eventFrame) {
            const auto& frameEvent = events[eventFrame];
            for (size_t writeIndex = 0; writeIndex < frameEvent.writes.size(); ++writeIndex) {
              const auto& write = frameEvent.writes[writeIndex];
              if (bit >= write.data.size()) {
                continue;  // LCOV_EXCL_LINE
              }
              RailLits hit = writeHits[eventFrame][writeIndex];
              if (bit < write.mask.size()) {
                hit = rails.railAnd(hit, write.mask[bit]);
              }
              value = rails.ite(hit, write.data[bit], value);
            }
            if (frameEvent.reset.has_value()) {
              value = rails.ite(*frameEvent.reset, rails.constRail(false), value);
            }
          }
          addLiteralEquivalence(
              solver, variables.getLiteral(port.outputRails[bit]->mayBeOne, frame + 1), value.one);
          addLiteralEquivalence(
              solver, variables.getLiteral(port.outputRails[bit]->mayBeZero, frame + 1), value.zero);
        }
      }
    }

    addFrameDomainConstraints(solver, variables, problem, frame + 1);
    holdInputs(frame + 1);
    int badLit = 0;
    int targetLit = 0;
    if (!targetAtLastFrameOnly || frame + 1 == totalFrames) {
      encodeTarget(frame + 1, badLit, targetLit);
      pendingBadLits.push_back(badLit);
      pendingTargetLits.push_back(targetLit);
    }
    logSeeded(diagnostics,
              "frame=%zu encoded states=%zu%s encode_seconds=%.1f frontier_units=%zu",
              frame + 1, encodedStates, needed == nullptr ? "" : " (target cone)",
              elapsedSeconds(frameStart), frontierUnits);
    // Every intermediate check is an unsatisfiability proof the search does
    // not need, so by default the target is checked once, at the last frame.
    const size_t checked = frame + 1;
    if (checked == totalFrames ||
        (!finalCheckOnly && (checked <= 2 || checked % checkEvery == 0))) {
      if (auto result = checkPendingTarget(checked)) {
        return result;
      }
    }
  }
  return std::nullopt;
}

struct PrefixTrial {
  size_t index = 0;
  // Registers that took two different known values during the prefix (toggle
  // coverage). A design held in reset toggles nothing, one stuck in an
  // initialization loop toggles only that loop's counters, and a design that
  // completed its start-up and went on to work toggles the most.
  size_t toggledRegisters = 0;
  size_t knownBits = 0;
  // Cycles in which some state bit changed.
  size_t activeCycles = 0;
  std::vector<bool> firstVector;
  std::vector<bool> secondVector;
  std::vector<TernaryValue> state0;
  std::vector<TernaryValue> state1;
};

struct SuffixTrace {
  std::vector<std::vector<bool>> inputsByFrame;  // [frame][alignedIndex]
  size_t mismatchFrame = 0;
  std::vector<KInductionResult::SignalMismatch> mismatches;
};

}  // namespace

std::optional<KInductionResult::CounterexampleWitness>
runSimulationSeededSearch(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& rails0ByKey,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& rails1ByKey,
    const KInductionProblem& problem,
    const SimulationSeededSearchOptions& options,
    KEPLER_FORMAL::Config::SolverType solverType,
    bool diagnostics) {
  if (!options.enabled()) {
    return std::nullopt;
  }
  if (!problem.usesDualRailStateEncoding) {
    logSeeded(diagnostics, "only the dual-rail encoding is supported; search skipped");
    return std::nullopt;
  }
  const size_t inputCount = alignedInputs.names.size();
  if (problem.inputSymbols.size() != inputCount ||
      alignedInputs.keys0.size() != inputCount ||
      alignedInputs.keys1.size() != inputCount) {
    logSeeded(diagnostics, "aligned inputs and proof inputs disagree; search skipped");  // LCOV_EXCL_LINE
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  const auto start = std::chrono::steady_clock::now();

  // Reset controls follow the bootstrap spec: active for its cycles, then
  // released. Everything else is a held random vector.
  std::vector<std::optional<bool>> resetActiveValue(inputCount);
  for (const auto& [symbol, activeValue] : problem.resetBootstrapInputs) {
    for (size_t index = 0; index < inputCount; ++index) {
      if (problem.inputSymbols[index] == symbol) {
        resetActiveValue[index] = activeValue;
      }
    }
  }
  const size_t prefixCycles = std::max(options.prefixCycles, problem.resetBootstrapCycles);
  const size_t holdCycles = std::min(options.holdCycles, prefixCycles);

  std::vector<size_t> localInput0(inputCount, kNoState);
  std::vector<size_t> localInput1(inputCount, kNoState);
  for (size_t index = 0; index < inputCount; ++index) {
    if (const auto it = model0.inputVarByKey.find(alignedInputs.keys0[index]);
        it != model0.inputVarByKey.end()) {
      localInput0[index] = it->second;
    }
    if (const auto it = model1.inputVarByKey.find(alignedInputs.keys1[index]);
        it != model1.inputVarByKey.end()) {
      localInput1[index] = it->second;
    }
  }

  TernaryModelSimulator simulator0(model0);
  TernaryModelSimulator simulator1(model1);
  logSeeded(diagnostics,
            "simulators ready states0=%zu states1=%zu inputs=%zu seconds=%.1f",
            simulator0.stateCount(), simulator1.stateCount(), inputCount,
            elapsedSeconds(start));

  const auto inputVectorAtCycle = [&](const PrefixTrial& trial, size_t cycle) {
    std::vector<bool> vector = cycle < holdCycles ? trial.firstVector : trial.secondVector;
    for (size_t index = 0; index < inputCount; ++index) {
      if (resetActiveValue[index].has_value()) {
        vector[index] = cycle < problem.resetBootstrapCycles ? *resetActiveValue[index]
                                                            : !*resetActiveValue[index];
      }
    }
    return vector;
  };
  const auto applyInputs = [&](const std::vector<bool>& vector) {
    for (size_t index = 0; index < inputCount; ++index) {
      const TernaryValue value = vector[index] ? TernaryValue::One : TernaryValue::Zero;
      if (localInput0[index] != kNoState) {
        simulator0.setInput(localInput0[index], value);
      }
      if (localInput1[index] != kNoState) {
        simulator1.setInput(localInput1[index], value);
      }
    }
  };
  const auto collectMismatches = [&]() {
    std::vector<KInductionResult::SignalMismatch> mismatches;
    for (size_t index = 0; index < alignedOutputs.names.size(); ++index) {
      const auto value0 = simulator0.evaluateOutput(alignedOutputs.keys0[index]);
      const auto value1 = simulator1.evaluateOutput(alignedOutputs.keys1[index]);
      if (isKnownTernaryValue(value0) && isKnownTernaryValue(value1) && value0 != value1) {
        mismatches.push_back({alignedOutputs.names[index],
                              value0 == TernaryValue::One,
                              value1 == TernaryValue::One});
      }
    }
    return mismatches;
  };

  std::mt19937_64 generator(options.seed);
  const auto randomVector = [&]() {
    std::vector<bool> vector(inputCount);
    for (size_t index = 0; index < inputCount; ++index) {
      vector[index] = (generator() & 1u) != 0;
    }
    for (size_t index = 0; index < inputCount; ++index) {
      if (resetActiveValue[index].has_value()) {
        vector[index] = !*resetActiveValue[index];
      }
    }
    return vector;
  };

  std::optional<PrefixTrial> best;
  std::vector<PrefixTrial> trials;
  for (size_t trialIndex = 0; trialIndex < options.trials; ++trialIndex) {
    PrefixTrial trial;
    trial.firstVector.resize(inputCount);
    trial.secondVector.resize(inputCount);
    for (size_t index = 0; index < inputCount; ++index) {
      trial.firstVector[index] = (generator() & 1u) != 0;
      trial.secondVector[index] = (generator() & 1u) != 0;
    }
    simulator0.resetToInitialState();
    simulator1.resetToInitialState();
    std::vector<bool> toggled0(simulator0.stateCount(), false);
    std::vector<bool> toggled1(simulator1.stateCount(), false);
    std::vector<TernaryValue> previous0 = simulator0.state();
    std::vector<TernaryValue> previous1 = simulator1.state();
    const auto noteToggles = [](const std::vector<TernaryValue>& previous,
                                const std::vector<TernaryValue>& current,
                                std::vector<bool>& toggled) {
      for (size_t index = 0; index < current.size(); ++index) {
        if (!toggled[index] && isKnownTernaryValue(previous[index]) &&
            isKnownTernaryValue(current[index]) && previous[index] != current[index]) {
          toggled[index] = true;
        }
      }
    };
    for (size_t cycle = 0; cycle < prefixCycles; ++cycle) {
      applyInputs(inputVectorAtCycle(trial, cycle));
      const bool changed0 = simulator0.step();
      const bool changed1 = simulator1.step();
      trial.activeCycles += (changed0 || changed1) ? 1 : 0;
      if (changed0) {
        noteToggles(previous0, simulator0.state(), toggled0);
        previous0 = simulator0.state();
      }
      if (changed1) {
        noteToggles(previous1, simulator1.state(), toggled1);
        previous1 = simulator1.state();
      }
    }
    trial.toggledRegisters = static_cast<size_t>(std::count(toggled0.begin(), toggled0.end(), true)) +
                             static_cast<size_t>(std::count(toggled1.begin(), toggled1.end(), true));
    trial.knownBits = simulator0.knownStateBitCount() + simulator1.knownStateBitCount();
    logSeeded(diagnostics,
              "trial=%zu toggled_registers=%zu known_state_bits=%zu active_cycles=%zu seconds=%.1f",
              trialIndex, trial.toggledRegisters, trial.knownBits, trial.activeCycles,
              elapsedSeconds(start));
    trial.index = trialIndex;
    trial.state0 = simulator0.state();
    trial.state1 = simulator1.state();
    trials.push_back(std::move(trial));
  }
  if (trials.empty()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  // Toggle coverage first: it separates a design that completed its start-up
  // from one stuck in reset or in an initialization loop, which a count of
  // known bits or of active cycles does not. Known bits then prefer the state
  // with the fewest unknowns, and active cycles break the last ties.
  const auto trialKey = [](const PrefixTrial& t) {
    return std::make_tuple(t.toggledRegisters, t.knownBits, t.activeCycles);
  };
  std::vector<size_t> ranking(trials.size());
  std::iota(ranking.begin(), ranking.end(), 0);
  std::stable_sort(ranking.begin(), ranking.end(), [&](size_t a, size_t b) {
    return trialKey(trials[a]) > trialKey(trials[b]);
  });
  best = trials[ranking.front()];
  logSeeded(diagnostics,
            "prefix_cycles=%zu best_toggled_registers=%zu known_state_bits=%zu of %zu "
            "active_cycles=%zu",
            prefixCycles, best->toggledRegisters, best->knownBits,
            simulator0.stateCount() + simulator1.stateCount(), best->activeCycles);
  if (diagnostics) {
    // The single-bit inputs of the selected prefix, first and second held
    // vector: these are the controls (enables, requests, invalidations) that
    // decide whether the prefix reached a useful state.
    std::unordered_map<std::string, size_t> bitsPerPortName;
    const auto portName = [](const std::string& name) {
      const auto bracket = name.rfind('[');
      return bracket == std::string::npos ? name : name.substr(0, bracket);
    };
    for (const auto& name : alignedInputs.names) {
      ++bitsPerPortName[portName(name)];
    }
    std::string held;
    for (size_t index = 0; index < inputCount; ++index) {
      if (bitsPerPortName[portName(alignedInputs.names[index])] != 1) {
        continue;
      }
      held += (held.empty() ? "" : " ") + portName(alignedInputs.names[index]) + "=" +
              (best->firstVector[index] ? "1" : "0") + "/" +
              (best->secondVector[index] ? "1" : "0");
    }
    logSeeded(diagnostics, "selected prefix single-bit inputs (first/second vector): %s",
              held.c_str());
    // Registers the selected state leaves unknown, by RTL net name: data that
    // was never written, or control that a held input never resolved.
    const auto nameOfState = [&](const SequentialDesignModel& model, size_t index) {
      const auto& key = model.stateBits[index];
      if (const auto it = model.outputNetNameByKey.find(key);
          it != model.outputNetNameByKey.end() && !it->second.empty()) {
        return it->second;
      }
      if (const auto it = model.displayNameByKey.find(key); it != model.displayNameByKey.end()) {
        return it->second;
      }
      return std::string("?");
    };
    std::string unknown;
    size_t unknownCount = 0;
    for (size_t index = 0; index < best->state0.size() && index < model0.stateBits.size();
         ++index) {
      if (isKnownTernaryValue(best->state0[index])) {
        continue;
      }
      ++unknownCount;
      if (unknownCount <= 48) {
        unknown += (unknown.empty() ? "" : " ") + nameOfState(model0, index);
      }
    }
    logSeeded(diagnostics, "selected prefix unknown registers in design 1: %zu%s%s",
              unknownCount, unknown.empty() ? "" : ": ", unknown.c_str());
  }

  // Replays `given` input frames from the two states, then tries random
  // continuations, and returns the first trace with an output mismatch. The
  // states after the given frames are reported for a solver-based follow-up.
  const auto replayForMismatch = [&](const std::vector<TernaryValue>& from0,
                                     const std::vector<TernaryValue>& from1,
                                     const std::vector<std::vector<bool>>& given,
                                     size_t extendFrames,
                                     size_t extendTries,
                                     std::vector<TernaryValue>& after0,
                                     std::vector<TernaryValue>& after1)
      -> std::optional<SuffixTrace> {
    simulator0.setState(from0);
    simulator1.setState(from1);
    SuffixTrace trace;
    for (size_t frame = 0; frame < given.size(); ++frame) {
      applyInputs(given[frame]);
      trace.inputsByFrame.push_back(given[frame]);
      trace.mismatches = collectMismatches();
      if (!trace.mismatches.empty()) {
        trace.mismatchFrame = frame;
        return trace;
      }
      simulator0.step();
      simulator1.step();
    }
    after0 = simulator0.state();
    after1 = simulator1.state();
    for (size_t attempt = 0; attempt < extendTries; ++attempt) {
      simulator0.setState(after0);
      simulator1.setState(after1);
      SuffixTrace extended = trace;
      for (size_t frame = 0; frame < extendFrames; ++frame) {
        const auto vector = randomVector();
        applyInputs(vector);
        extended.inputsByFrame.push_back(vector);
        extended.mismatches = collectMismatches();
        if (!extended.mismatches.empty()) {
          extended.mismatchFrame = extended.inputsByFrame.size() - 1;
          return extended;
        }
        simulator0.step();
        simulator1.step();
      }
    }
    return std::nullopt;
  };
  const auto frontierOf = [&](const std::vector<TernaryValue>& state0,
                              const std::vector<TernaryValue>& state1) {
    std::unordered_map<size_t, bool> frontier;
    frontier.reserve(2 * (state0.size() + state1.size()));
    const auto add = [&](const SequentialDesignModel& model,
                         const std::vector<TernaryValue>& state,
                         const std::unordered_map<SignalKey, DualRailSymbolPair,
                                                  SignalKeyHash>& railsByKey) {
      for (size_t index = 0; index < model.stateBits.size(); ++index) {
        const auto railsIt = railsByKey.find(model.stateBits[index]);
        if (railsIt == railsByKey.end()) {
          continue;  // LCOV_EXCL_LINE
        }
        const auto bits = static_cast<uint8_t>(state[index]);
        frontier[railsIt->second.mayBeOne] = (bits & 1u) != 0;
        frontier[railsIt->second.mayBeZero] = (bits & 2u) != 0;
      }
    };
    add(model0, state0, rails0ByKey);
    add(model1, state1, rails1ByKey);
    return frontier;
  };
  TargetEnlargement enlargement;
  enlargement.structuralSteps = options.enlargeSteps;
  enlargement.cubeSteps = options.enlargeCubes > 0 ? options.enlargeSteps : 0;
  enlargement.cubeBudget = options.enlargeCubes;
  const TargetEnlargement noEnlargement;

  // Single-bit inputs that feed the next state of most registers are global
  // controls (resets, scan and test enables). Letting a bounded search drive
  // them turns every hop into a reset; they keep the value the selected prefix
  // ran with (held inputs as in simulation-seeded formal search). Bits of a
  // bus are data however many registers they feed: at module scope a memory
  // response bus reaches most of the registers, and holding it would freeze
  // the very values the search has to choose.
  std::vector<std::pair<size_t, bool>> heldInputIndices;  // (aligned index, value)
  {
    std::vector<size_t> registersFedByInput(inputCount, 0);
    const auto lastVector = inputVectorAtCycle(*best, prefixCycles - 1);
    std::unordered_map<std::string, size_t> bitsPerPort;
    const auto portOf = [](const std::string& name) {
      const auto bracket = name.rfind('[');
      return bracket == std::string::npos ? name : name.substr(0, bracket);
    };
    for (size_t index = 0; index < inputCount; ++index) {
      ++bitsPerPort[portOf(alignedInputs.names[index])];
    }
    std::unordered_map<size_t, size_t> indexByLocalInput0;
    for (size_t index = 0; index < inputCount; ++index) {
      if (localInput0[index] != kNoState) {
        indexByLocalInput0.emplace(localInput0[index], index);
      }
    }
    // Memory cells and read outputs all depend on the write controls of their
    // memory; they say nothing about global controls and are left out.
    std::unordered_set<SignalKey, SignalKeyHash> memoryKeys;
    for (const auto& memory : model0.structuredMemories) {
      for (const auto& cellKeys : memory.cellKeys) {
        memoryKeys.insert(cellKeys.begin(), cellKeys.end());
      }
      for (const auto& port : memory.readPorts) {
        memoryKeys.insert(port.outputKeys.begin(), port.outputKeys.end());
      }
    }
    size_t registers = 0;
    for (const auto& key : model0.stateBits) {
      const auto nextIt = model0.nextStateExprByStateKey.find(key);
      if (nextIt == model0.nextStateExprByStateKey.end() || nextIt->second == nullptr ||
          memoryKeys.count(key) != 0) {
        continue;
      }
      ++registers;
      for (const auto var : nextIt->second->getSupportVars()) {
        if (const auto it = indexByLocalInput0.find(var); it != indexByLocalInput0.end()) {
          ++registersFedByInput[it->second];
        }
      }
    }
    for (size_t index = 0; index < inputCount; ++index) {
      if (registers > 0 && registersFedByInput[index] * 2 >= registers &&
          bitsPerPort[portOf(alignedInputs.names[index])] == 1 &&
          !resetActiveValue[index].has_value()) {
        heldInputIndices.emplace_back(index, lastVector[index]);
        logSeeded(diagnostics, "held input %s=%d (feeds %zu of %zu registers)",
                  alignedInputs.names[index].c_str(), lastVector[index] ? 1 : 0,
                  registersFedByInput[index], registers);
      }
    }
  }
  const auto heldFor = [&](const KInductionProblem& target) {
    std::vector<std::pair<size_t, bool>> held;
    for (const auto& [index, value] : heldInputIndices) {
      if (index < target.inputSymbols.size()) {
        held.emplace_back(target.inputSymbols[index], value);
      }
    }
    return held;
  };
  const auto heldProblemInputs = heldFor(problem);

  // Replays `frames` input vectors from the two states; returns the states
  // reached.
  const auto statesAfter = [&](const std::vector<TernaryValue>& from0,
                               const std::vector<TernaryValue>& from1,
                               const std::vector<std::vector<bool>>& inputs,
                               size_t frames,
                               std::vector<TernaryValue>& after0,
                               std::vector<TernaryValue>& after1) {
    simulator0.setState(from0);
    simulator1.setState(from1);
    for (size_t frame = 0; frame < frames && frame < inputs.size(); ++frame) {
      applyInputs(inputs[frame]);
      simulator0.step();
      simulator1.step();
    }
    after0 = simulator0.state();
    after1 = simulator1.state();
  };
  // A hit on an enlarged target lies within the enlargement's step count of
  // the real target: finish the trace with an exact search from the state
  // reached, replacing the arbitrary inputs of the hit frame.
  const auto finishEnlarged =
      [&](FrontierBoundedSearchResult result,
          const std::vector<TernaryValue>& from0,
          const std::vector<TernaryValue>& from1,
          const auto& exactSearchFrom) -> std::optional<FrontierBoundedSearchResult> {
    if (!result.viaIntermediateTarget) {
      return result;
    }
    std::vector<TernaryValue> reached0;
    std::vector<TernaryValue> reached1;
    statesAfter(from0, from1, result.inputsByFrame, result.badFrame, reached0, reached1);
    logSeeded(diagnostics,
              "enlarged target hit at suffix frame %zu; exact search of %zu frames from there",
              result.badFrame, enlargement.cubeSteps);
    const auto rest = exactSearchFrom(reached0, reached1, enlargement.cubeSteps);
    if (!rest.has_value()) {
      logSeeded(diagnostics, "enlarged target hit could not be completed; discarded");  // LCOV_EXCL_LINE
      return std::nullopt;  // LCOV_EXCL_LINE
    }
    FrontierBoundedSearchResult combined;
    combined.inputsByFrame.assign(result.inputsByFrame.begin(),
                                  result.inputsByFrame.begin() + result.badFrame);
    combined.inputsByFrame.insert(combined.inputsByFrame.end(),
                                  rest->inputsByFrame.begin(), rest->inputsByFrame.end());
    combined.badFrame = result.badFrame + rest->badFrame;
    return combined;
  };

  // Two-copy bounded search from a pair of states, confirmed by replay.
  const auto twoCopyExact = [&](const std::vector<TernaryValue>& from0,
                                const std::vector<TernaryValue>& from1,
                                size_t depth) {
    return runFrontierBoundedSearch(
        problem, frontierOf(from0, from1), {}, depth, options.checkEvery, noEnlargement,
        nullptr, solverType, diagnostics, nullptr, false, 0, heldProblemInputs);
  };
  const auto twoCopySearch = [&](const std::vector<TernaryValue>& from0,
                                 const std::vector<TernaryValue>& from1,
                                 size_t depth) -> std::optional<SuffixTrace> {
    auto search = runFrontierBoundedSearch(
        problem, frontierOf(from0, from1), {}, depth, options.checkEvery, enlargement,
        nullptr, solverType, diagnostics, nullptr, false, 0, heldProblemInputs);
    if (search.has_value()) {
      search = finishEnlarged(*search, from0, from1, twoCopyExact);
    }
    if (!search.has_value()) {
      return std::nullopt;
    }
    std::vector<TernaryValue> after0;
    std::vector<TernaryValue> after1;
    auto trace = replayForMismatch(from0, from1, search->inputsByFrame, 0, 0, after0, after1);
    if (!trace.has_value()) {
      logSeeded(diagnostics, "candidate trace not confirmed by simulation; discarded");  // LCOV_EXCL_LINE
    }
    return trace;
  };

  // The guide for abstraction and waypoints: the origins of the register
  // differences when known, else the register differences, else the whole
  // predicate.
  const auto guideOf = [](const DivergenceProblem& divergence) {
    if (divergence.originDivergence != nullptr &&
        divergence.originDivergence != BoolExpr::createFalse()) {
      return divergence.originDivergence;
    }
    if (divergence.registerDivergence != nullptr &&
        divergence.registerDivergence != BoolExpr::createFalse()) {
      return divergence.registerDivergence;
    }
    return divergence.problem.bad;
  };

  // Waypoint chaining on the single-copy divergence problem (lighthouses:
  // Ganai, Aziz and Kuehlmann, DAC 1999; Ho et al., ICCAD 2000; distance
  // guidance: Kuehlmann, McMillan and Brayton, ICCAD 1999; Shyam and Bertacco,
  // DATE 2006). A target state is sampled as the start of a short free-start
  // search that ends in the divergence, so its registers are consistent with
  // reaching it; each hop is a short search that brings one more register of
  // the predicate's cone to its target value without losing an aligned one.
  // When the sample gives no reachable direction, a hop may instead change
  // any cone register that has not changed yet. Hops are chained from the
  // concrete state each one reaches; the real target is part of every goal.
  struct WaypointChainOutcome {
    std::optional<FrontierBoundedSearchResult> divergence;  // the real target was reached
    std::vector<std::vector<bool>> inputs;                   // chain inputs so far
    std::vector<TernaryValue> state0;
    std::vector<TernaryValue> state1;
  };
  const auto runWaypointChain = [&](const DivergenceProblem& divergence) -> WaypointChainOutcome {
    WaypointChainOutcome outcome;
    outcome.state0 = best->state0;
    outcome.state1 = best->state1;
    const auto chainStart = std::chrono::steady_clock::now();
    constexpr size_t kMaxConeLayers = 32;
    constexpr size_t kMaxTargetSamples = 4;
    BoolExpr* guide = guideOf(divergence);
    const auto layers = targetConeLayers(divergence.problem, guide, kMaxConeLayers);
    std::unordered_map<size_t, size_t> indexByRail;
    indexByRail.reserve(2 * divergence.railsByStateIndex.size());
    for (size_t index = 0; index < divergence.railsByStateIndex.size(); ++index) {
      indexByRail.emplace(divergence.railsByStateIndex[index].mayBeOne, index);
      indexByRail.emplace(divergence.railsByStateIndex[index].mayBeZero, index);
    }
    // layerIndices[L]: registers within L steps of the predicate (cumulative).
    std::vector<std::vector<size_t>> layerIndices(layers.size());
    for (size_t layer = 0; layer < layers.size(); ++layer) {
      std::unordered_set<size_t> seen;
      for (const auto symbol : layers[layer]) {
        if (const auto it = indexByRail.find(symbol);
            it != indexByRail.end() && seen.insert(it->second).second) {
          layerIndices[layer].push_back(it->second);
        }
      }
    }
    const std::vector<size_t>& coneIndices = layerIndices.back();
    logSeeded(diagnostics, "waypoints: cone layers=%zu registers=%zu",
              layers.size(), coneIndices.size());

    // Registers that changed during the prefix are not flip waypoints.
    simulator0.resetToInitialState();
    simulator1.resetToInitialState();
    simulator0.clearToggles();
    for (size_t cycle = 0; cycle < prefixCycles; ++cycle) {
      applyInputs(inputVectorAtCycle(*best, cycle));
      simulator0.step();
      simulator1.step();
    }
    std::vector<bool> retired = simulator0.toggled();

    // Memory contents of the selected state: a sampled start state keeps them,
    // every register is free.
    std::unordered_map<size_t, bool> cellFrontier;
    if (divergence.problem.structuredMemories != nullptr) {
      for (const auto& memory : *divergence.problem.structuredMemories) {
        for (const auto& cell : memory.cellRails) {
          for (const auto& rails : cell) {
            for (const auto symbol : {rails.mayBeOne, rails.mayBeZero}) {
              if (const auto it = divergence.frontier.find(symbol); it != divergence.frontier.end()) {
                cellFrontier.emplace(symbol, it->second);
              }
            }
          }
        }
      }
    }
    std::vector<TernaryValue> target(outcome.state0.size(), TernaryValue::Unknown);
    size_t samples = 0;
    const auto nameOf = [&](size_t index) {
      const auto& key = model0.stateBits[index];
      const auto it = model0.displayNameByKey.find(key);
      std::string name = it == model0.displayNameByKey.end() ? std::string("?") : it->second;
      if (const auto netIt = model0.outputNetNameByKey.find(key);
          netIt != model0.outputNetNameByKey.end() && !netIt->second.empty()) {
        name += " (" + netIt->second + ")";
      }
      return name;
    };
    const auto sampleTarget = [&]() {
      const size_t depth = options.waypointDepth + 2 * samples;
      std::vector<std::pair<size_t, bool>> startState;
      const auto sample = runFrontierBoundedSearch(
          divergence.problem, cellFrontier, {}, depth, options.checkEvery, noEnlargement,
          nullptr, solverType, false, &startState, /*targetAtLastFrameOnly=*/true, 0,
          heldFor(divergence.problem));
      ++samples;
      if (!sample.has_value()) {
        logSeeded(diagnostics, "waypoints: no state reaches the target within %zu frames", depth);
        return false;
      }
      std::fill(target.begin(), target.end(), TernaryValue::Unknown);
      std::unordered_map<size_t, bool> value(startState.begin(), startState.end());
      size_t known = 0;
      for (const auto index : coneIndices) {
        const auto& rails = divergence.railsByStateIndex[index];
        const auto oneIt = value.find(rails.mayBeOne);
        const auto zeroIt = value.find(rails.mayBeZero);
        if (oneIt == value.end() || zeroIt == value.end()) {
          continue;  // LCOV_EXCL_LINE
        }
        if (oneIt->second && !zeroIt->second) {
          target[index] = TernaryValue::One;
          ++known;
        } else if (zeroIt->second && !oneIt->second) {
          target[index] = TernaryValue::Zero;
          ++known;
        }
      }
      logSeeded(diagnostics,
                "waypoints: target sample %zu (a state %zu frames before the divergence) fixes %zu of %zu cone registers",
                samples, sample->badFrame, known, coneIndices.size());
      return true;
    };
    if (!sampleTarget()) {
      return outcome;
    }
    const auto equalsTarget = [&](size_t index) {
      const auto& rails = divergence.railsByStateIndex[index];
      BoolExpr* one = BoolExpr::Var(rails.mayBeOne);
      BoolExpr* zero = BoolExpr::Var(rails.mayBeZero);
      return target[index] == TernaryValue::One ? BoolExpr::And(one, BoolExpr::Not(zero))
                                                : BoolExpr::And(zero, BoolExpr::Not(one));
    };
    const auto flipExpr = [&](size_t index) {
      const auto& rails = divergence.railsByStateIndex[index];
      BoolExpr* one = BoolExpr::Var(rails.mayBeOne);
      BoolExpr* zero = BoolExpr::Var(rails.mayBeZero);
      BoolExpr* binaryOne = BoolExpr::And(one, BoolExpr::Not(zero));
      BoolExpr* binaryZero = BoolExpr::And(zero, BoolExpr::Not(one));
      return outcome.state0[index] == TernaryValue::One ? binaryZero
             : outcome.state0[index] == TernaryValue::Zero ? binaryOne
                                                           : BoolExpr::Or(binaryOne, binaryZero);
    };
    const auto alignedCount = [&]() {
      size_t count = 0;
      for (const auto index : coneIndices) {
        count += (target[index] != TernaryValue::Unknown &&
                  outcome.state0[index] == target[index]) ? 1 : 0;
      }
      return count;
    };

    enum class HopMode { Directed, DirectedRelaxed, Flip };
    size_t layer = std::min(options.waypointLayers, layers.size() - 1);
    HopMode mode = HopMode::Directed;
    size_t hop = 0;
    // Registers the chain itself brought to their target value: a directed
    // hop must not undo them. Registers aligned from the start are left free,
    // since holding thousands of them makes hops hard to refute.
    std::unordered_set<size_t> gained;
    while (hop < options.waypointHops) {
      std::vector<size_t> candidates;
      std::vector<size_t> aligned;
      for (const auto index : layerIndices[layer]) {
        if (mode == HopMode::Flip) {
          if (!retired[index]) {
            candidates.push_back(index);
          }
        } else if (target[index] != TernaryValue::Unknown) {
          if (outcome.state0[index] == target[index]) {
            if (gained.count(index) != 0) {
              aligned.push_back(index);
            }
          } else {
            candidates.push_back(index);
          }
        }
      }
      const auto nextAttempt = [&]() {
        // Directed, directed without holding aligned registers, flips; then a
        // wider layer; then another sample.
        if (mode == HopMode::Directed) {
          mode = HopMode::DirectedRelaxed;
          return true;
        }
        if (mode == HopMode::DirectedRelaxed) {
          mode = HopMode::Flip;
          return true;
        }
        mode = HopMode::Directed;
        if (layer + 1 < layers.size()) {
          ++layer;
          return true;
        }
        if (samples < kMaxTargetSamples && sampleTarget()) {
          layer = std::min(options.waypointLayers, layers.size() - 1);
          return true;
        }
        return false;
      };
      if (candidates.empty()) {
        if (nextAttempt()) {
          continue;
        }
        logSeeded(diagnostics, "waypoints: no waypoint left; stopping");
        break;
      }
      BoolExpr* waypoint = BoolExpr::createFalse();
      for (const auto index : candidates) {
        waypoint = BoolExpr::Or(waypoint, mode == HopMode::Flip ? flipExpr(index) : equalsTarget(index));
      }
      if (mode == HopMode::Directed) {
        for (const auto index : aligned) {
          waypoint = BoolExpr::And(waypoint, equalsTarget(index));
        }
      }
      if (diagnostics && candidates.size() <= 8) {
        for (const auto index : candidates) {
          logSeeded(diagnostics, "waypoints:   candidate %s now=%d target=%d",
                    nameOf(index).c_str(), static_cast<int>(outcome.state0[index]),
                    static_cast<int>(target[index]));
        }
      }
      logSeeded(diagnostics,
                "waypoints: hop=%zu layer=%zu mode=%s candidates=%zu aligned=%zu depth=%zu e.g. %s (%.1f s)",
                hop, layer,
                mode == HopMode::Directed ? "directed"
                : mode == HopMode::DirectedRelaxed ? "directed-relaxed" : "flip",
                candidates.size(), aligned.size(), options.waypointDepth,
                nameOf(candidates.front()).c_str(), elapsedSeconds(chainStart));
      const auto result = runFrontierBoundedSearch(
          divergence.problem, divergence.frontierFor(outcome.state0), {},
          options.waypointDepth, options.checkEvery, noEnlargement, waypoint, solverType,
          diagnostics, nullptr, false, options.waypointConflicts, heldFor(divergence.problem));
      ++hop;
      if (!result.has_value()) {
        if (nextAttempt()) {
          continue;
        }
        logSeeded(diagnostics, "waypoints: nothing reachable within %zu frames; stopping",
                  options.waypointDepth);
        break;
      }
      if (!result->viaIntermediateTarget) {
        outcome.inputs.insert(outcome.inputs.end(), result->inputsByFrame.begin(),
                              result->inputsByFrame.end());
        FrontierBoundedSearchResult total;
        total.badFrame = outcome.inputs.size() - 1;
        total.inputsByFrame = outcome.inputs;
        outcome.divergence = std::move(total);
        logSeeded(diagnostics,
                  "waypoints: target reached after %zu hops at suffix frame %zu (%.1f s)",
                  hop, outcome.inputs.size() - 1, elapsedSeconds(chainStart));
        return outcome;
      }
      const size_t alignedBefore = alignedCount();
      simulator0.setState(outcome.state0);
      simulator1.setState(outcome.state1);
      std::vector<bool> changed(outcome.state0.size(), false);
      for (size_t frame = 0; frame < result->badFrame; ++frame) {
        applyInputs(result->inputsByFrame[frame]);
        const auto before = simulator0.state();
        simulator0.step();
        simulator1.step();
        const auto& after = simulator0.state();
        for (size_t index = 0; index < changed.size(); ++index) {
          changed[index] = changed[index] || before[index] != after[index];
        }
        outcome.inputs.push_back(result->inputsByFrame[frame]);
      }
      size_t newlyRetired = 0;
      for (size_t index = 0; index < retired.size(); ++index) {
        if (changed[index] && !retired[index]) {
          retired[index] = true;
          ++newlyRetired;
        }
      }
      outcome.state0 = simulator0.state();
      outcome.state1 = simulator1.state();
      for (const auto index : coneIndices) {
        if (changed[index] && target[index] != TernaryValue::Unknown &&
            outcome.state0[index] == target[index]) {
          gained.insert(index);
        }
      }
      const size_t alignedAfter = alignedCount();
      logSeeded(diagnostics,
                "waypoints: hop %zu reached a waypoint at frame %zu; aligned %zu -> %zu of %zu, newly changed=%zu",
                hop, result->badFrame, alignedBefore, alignedAfter, coneIndices.size(), newlyRetired);
      const bool progressed = mode == HopMode::Flip ? newlyRetired > 0 : alignedAfter > alignedBefore;
      if (progressed) {
        mode = HopMode::Directed;
        layer = std::min(options.waypointLayers, layers.size() - 1);
      } else if (!nextAttempt()) {
        logSeeded(diagnostics, "waypoints: no progress; stopping");
        break;
      }
    }
    if (hop >= options.waypointHops) {
      logSeeded(diagnostics, "waypoints: hop budget of %zu exhausted", options.waypointHops);
    }
    return outcome;
  };

  // Abstraction-guided search on the single-copy divergence problem
  // (localization: Kurshan 1994; counterexample-guided refinement: Clarke,
  // Grumberg, Jha, Lu and Veith, CAV 2000; guidance by abstract traces: Nanshi
  // and Somenzi, DAC 2006). The abstract model keeps the registers within a
  // few steps of the divergence predicate and leaves every other register
  // free. Its trace pins the kept part's boundary registers frame by frame in
  // the concrete search, which then only has to find inputs that reproduce
  // them. When the pinned trace is infeasible, the registers whose pins took
  // part in the refutation join the abstract model, and the loop repeats.
  const auto runAbstractionGuided =
      [&](const DivergenceProblem& divergence) -> std::optional<FrontierBoundedSearchResult> {
    const auto abstractionStart = std::chrono::steady_clock::now();
    constexpr size_t kMaxConeLayers = 32;
    BoolExpr* guide = guideOf(divergence);
    const auto layers = targetConeLayers(divergence.problem, guide, kMaxConeLayers);
    const TransitionExprResolver transitions(divergence.problem);
    std::unordered_map<size_t, size_t> indexByRail;
    indexByRail.reserve(2 * divergence.railsByStateIndex.size());
    for (size_t index = 0; index < divergence.railsByStateIndex.size(); ++index) {
      indexByRail.emplace(divergence.railsByStateIndex[index].mayBeOne, index);
      indexByRail.emplace(divergence.railsByStateIndex[index].mayBeZero, index);
    }
    // Memory cells are never model registers here: reads go through the port
    // chains, so cells are neither kept nor boundary symbols.
    std::unordered_set<size_t> cellRails;
    if (divergence.problem.structuredMemories != nullptr) {
      for (const auto& memory : *divergence.problem.structuredMemories) {
        for (const auto& cell : memory.cellRails) {
          for (const auto& rails : cell) {
            cellRails.insert(rails.mayBeOne);
            cellRails.insert(rails.mayBeZero);
          }
        }
      }
    }
    std::unordered_set<size_t> kept;
    size_t keptRegisters = 0;
    const auto keepRegister = [&](size_t index) {
      const auto& rails = divergence.railsByStateIndex[index];
      if (cellRails.count(rails.mayBeOne) != 0) {
        return;
      }
      if (kept.insert(rails.mayBeOne).second) {
        ++keptRegisters;
      }
      kept.insert(rails.mayBeZero);
    };
    // abstractionLayers = 1 keeps the predicate's own support (layer 0).
    size_t layer = std::min(options.abstractionLayers - 1, layers.size() - 1);
    const auto keepLayer = [&](size_t which) {
      for (const auto symbol : layers[which]) {
        if (const auto it = indexByRail.find(symbol); it != indexByRail.end()) {
          keepRegister(it->second);
        }
      }
    };
    // Instance scope: the component that holds the differing registers keeps
    // its logic, its environment is free (compositional reasoning with the
    // environment abstracted: Clarke, Long and McMillan, LICS 1989; Grumberg
    // and Long, TOPLAS 1994). The instance is the deepest common ancestor of
    // the differing registers, widened until it holds more than them.
    bool instanceScoped = false;
    if (options.abstractionByInstance && !divergence.differingIndices.empty()) {
      const auto nameOf = [&](size_t index) -> std::string {
        if (index >= model0.stateBits.size()) {
          return std::string();
        }
        const auto it = model0.displayNameByKey.find(model0.stateBits[index]);
        return it == model0.displayNameByKey.end() ? std::string() : it->second;
      };
      const auto parentOf = [](const std::string& name) {
        const auto dot = name.rfind('.');
        return dot == std::string::npos ? std::string() : name.substr(0, dot);
      };
      const auto commonInstance = [](const std::string& a, const std::string& b) {
        size_t common = 0;  // length of the common leading elements
        size_t pos = 0;
        while (true) {
          const auto dotA = a.find('.', pos);
          const auto dotB = b.find('.', pos);
          const auto endA = dotA == std::string::npos ? a.size() : dotA;
          const auto endB = dotB == std::string::npos ? b.size() : dotB;
          if (endA != endB || a.compare(pos, endA - pos, b, pos, endB - pos) != 0) {
            break;
          }
          common = endA;
          if (dotA == std::string::npos || dotB == std::string::npos) {
            break;
          }
          pos = dotA + 1;
        }
        return a.substr(0, common);
      };
      std::string scope;
      bool first = true;
      for (const auto index : divergence.differingIndices) {
        const auto parent = parentOf(nameOf(index));
        scope = first ? parent : commonInstance(scope, parent);
        first = false;
      }
      const std::unordered_set<size_t> differingSet(divergence.differingIndices.begin(),
                                                    divergence.differingIndices.end());
      std::vector<size_t> under;
      while (true) {
        under.clear();
        bool beyondDiffering = false;
        for (size_t index = 0; index < divergence.railsByStateIndex.size(); ++index) {
          if (cellRails.count(divergence.railsByStateIndex[index].mayBeOne) != 0) {
            continue;
          }
          const auto name = nameOf(index);
          if (!scope.empty() &&
              (name.size() <= scope.size() || name.compare(0, scope.size(), scope) != 0 ||
               name[scope.size()] != '.')) {
            continue;
          }
          under.push_back(index);
          beyondDiffering = beyondDiffering || differingSet.count(index) == 0;
        }
        if (beyondDiffering || scope.empty()) {
          break;
        }
        scope = parentOf(scope);
      }
      for (const auto index : under) {
        keepRegister(index);
      }
      instanceScoped = true;
      logSeeded(diagnostics, "abstraction: instance scope '%s' keeps %zu registers",
                scope.empty() ? "(top)" : scope.c_str(), under.size());
    }
    if (!instanceScoped) {
      keepLayer(layer);
    }
    KInductionProblem abstractProblem = divergence.problem;
    abstractProblem.bad = guide;
    abstractProblem.property = BoolExpr::simplify(BoolExpr::Not(guide));
    const auto cutSymbols = [&]() {
      std::vector<size_t> cut;
      std::unordered_set<size_t> seen;
      const auto addSupport = [&](const std::set<size_t>& support) {
        for (const auto symbol : support) {
          if (indexByRail.count(symbol) != 0 && kept.count(symbol) == 0 &&
              cellRails.count(symbol) == 0 && seen.insert(symbol).second) {
            cut.push_back(symbol);
          }
        }
      };
      for (const auto symbol : kept) {
        if (transitions.contains(symbol)) {
          addSupport(transitions.support(symbol));
        }
      }
      if (divergence.problem.structuredMemories != nullptr) {
        for (const auto& memory : *divergence.problem.structuredMemories) {
          bool keptOutput = false;
          for (const auto& port : memory.readPorts) {
            for (const auto& rails : port.outputRails) {
              keptOutput = keptOutput || (rails.has_value() && kept.count(rails->mayBeOne) != 0);
            }
          }
          if (!keptOutput) {
            continue;
          }
          const auto addRails = [&](const DualRailBoolExpr& expr) {
            addSupport(expr.mayBeOne->getSupportVars());
            addSupport(expr.mayBeZero->getSupportVars());
          };
          for (const auto& port : memory.readPorts) {
            for (const auto& bit : port.address) {
              addRails(bit);
            }
          }
          for (const auto& port : memory.writePorts) {
            for (const auto* word : {&port.address, &port.data, &port.mask, &port.enable}) {
              for (const auto& bit : *word) {
                addRails(bit);
              }
            }
          }
          if (memory.resetAsserted.has_value()) {
            addRails(*memory.resetAsserted);
          }
        }
      }
      return cut;
    };
    const auto widen = [&]() {
      if (layer + 1 >= layers.size()) {
        return false;
      }
      ++layer;
      keepLayer(layer);
      return true;
    };
    for (size_t iteration = 0; iteration < options.abstractionIterations; ++iteration) {
      const auto cut = cutSymbols();
      logSeeded(diagnostics,
                "abstraction: iteration=%zu kept_registers=%zu of %zu boundary_registers=%zu (%.1f s)",
                iteration, keptRegisters,
                divergence.railsByStateIndex.size() - cellRails.size() / 2, cut.size() / 2,
                elapsedSeconds(abstractionStart));
      if (diagnostics) {
        std::vector<size_t> keptIndices;
        std::vector<size_t> cutIndices;
        for (const auto symbol : kept) {
          if (const auto it = indexByRail.find(symbol);
              it != indexByRail.end() &&
              divergence.railsByStateIndex[it->second].mayBeOne == symbol) {
            keptIndices.push_back(it->second);
          }
        }
        for (const auto symbol : cut) {
          if (const auto it = indexByRail.find(symbol);
              it != indexByRail.end() &&
              divergence.railsByStateIndex[it->second].mayBeOne == symbol) {
            cutIndices.push_back(it->second);
          }
        }
        std::sort(keptIndices.begin(), keptIndices.end());
        std::sort(cutIndices.begin(), cutIndices.end());
        logRegisterHistogram(diagnostics, "abstraction: kept", model0, keptIndices);
        logRegisterHistogram(diagnostics, "abstraction: boundary", model0, cutIndices);
      }
      std::vector<std::vector<std::pair<size_t, bool>>> trace;
      BoundedSearchExtras abstractExtras;
      abstractExtras.keptStates = &kept;
      abstractExtras.traceSymbols = &cut;
      abstractExtras.traceValuesOut = &trace;
      // Abstract traces start short and grow with each refinement, so the
      // abstract search checks the target incrementally at every frame.
      abstractExtras.forceAssumptionSolver = true;
      const auto abstractResult = runFrontierBoundedSearch(
          abstractProblem, divergence.frontier, {}, options.searchDepth,
          options.checkEvery == 0 ? 1 : options.checkEvery, noEnlargement, nullptr, solverType,
          diagnostics, nullptr, false, 0, heldFor(divergence.problem), abstractExtras);
      if (!abstractResult.has_value()) {
        logSeeded(diagnostics,
                  "abstraction: the abstract model has no divergence within %zu frames (%.1f s)",
                  options.searchDepth, elapsedSeconds(abstractionStart));
        return std::nullopt;
      }
      if (abstractResult->badFrame == 0) {
        logSeeded(diagnostics,
                  "abstraction: the guide already holds at the frontier; not applicable (%.1f s)",
                  elapsedSeconds(abstractionStart));
        return std::nullopt;
      }
      // Pin the boundary registers' binary values along the abstract trace.
      std::vector<std::tuple<size_t, size_t, bool>> pins;
      for (size_t frame = 1; frame < trace.size(); ++frame) {
        std::unordered_map<size_t, bool> value(trace[frame].begin(), trace[frame].end());
        std::unordered_set<size_t> done;
        for (const auto symbol : cut) {
          const auto indexIt = indexByRail.find(symbol);
          if (indexIt == indexByRail.end() || !done.insert(indexIt->second).second) {
            continue;
          }
          const auto& rails = divergence.railsByStateIndex[indexIt->second];
          const auto oneIt = value.find(rails.mayBeOne);
          const auto zeroIt = value.find(rails.mayBeZero);
          if (oneIt == value.end() || zeroIt == value.end() || oneIt->second == zeroIt->second) {
            continue;  // unknown or absent: no pin
          }
          pins.emplace_back(rails.mayBeOne, frame, oneIt->second);
          pins.emplace_back(rails.mayBeZero, frame, zeroIt->second);
        }
      }
      logSeeded(diagnostics, "abstraction: abstract trace of %zu frames; %zu pinned rail values",
                abstractResult->badFrame, pins.size());
      // Lift the pins to those the abstract model needs: with every pin and
      // every input of the abstract trace assumed, the negated guide at the
      // last frame is unsatisfiable, and the failed assumptions are the
      // relevant ones (minimal assignments by unsatisfiable core: Ravi and
      // Somenzi, TACAS 2004; Een, Mishchenko and Brayton, FMCAD 2011).
      {
        std::vector<std::tuple<size_t, size_t, bool>> fixed = pins;
        for (size_t frame = 0; frame < abstractResult->inputsByFrame.size(); ++frame) {
          const auto& inputs = abstractResult->inputsByFrame[frame];
          for (size_t i = 0; i < inputs.size() && i < abstractProblem.inputSymbols.size(); ++i) {
            fixed.emplace_back(abstractProblem.inputSymbols[i], frame, inputs[i]);
          }
        }
        KInductionProblem liftProblem = abstractProblem;
        liftProblem.bad = abstractProblem.property;
        liftProblem.property = abstractProblem.bad;
        std::vector<std::tuple<size_t, size_t, bool>> needed;
        bool liftUnsat = false;
        BoundedSearchExtras liftExtras;
        liftExtras.keptStates = &kept;
        liftExtras.pins = &fixed;
        liftExtras.failedPinsOut = &needed;
        liftExtras.forceAssumptionSolver = true;
        liftExtras.unsatOut = &liftUnsat;
        runFrontierBoundedSearch(liftProblem, divergence.frontier, {}, abstractResult->badFrame,
                                 0, noEnlargement, nullptr, solverType, false, nullptr,
                                 /*targetAtLastFrameOnly=*/true, options.abstractionConflicts,
                                 heldFor(divergence.problem), liftExtras);
        if (liftUnsat) {
          std::vector<std::tuple<size_t, size_t, bool>> lifted;
          for (const auto& pin : needed) {
            if (indexByRail.count(std::get<0>(pin)) != 0) {
              lifted.push_back(pin);
            }
          }
          logSeeded(diagnostics, "abstraction: lifted pins %zu -> %zu", pins.size(),
                    lifted.size());
          pins = std::move(lifted);
        } else {
          logSeeded(diagnostics, "abstraction: pin lifting inconclusive; keeping every pin");
        }
      }
      std::vector<std::tuple<size_t, size_t, bool>> failed;
      bool unsat = false;
      BoundedSearchExtras concreteExtras;
      concreteExtras.pins = &pins;
      concreteExtras.failedPinsOut = &failed;
      concreteExtras.forceAssumptionSolver = true;
      concreteExtras.unsatOut = &unsat;
      const auto concrete = runFrontierBoundedSearch(
          divergence.problem, divergence.frontier, {}, abstractResult->badFrame, 0,
          noEnlargement, nullptr, solverType, diagnostics, nullptr, false,
          options.abstractionConflicts, heldFor(divergence.problem), concreteExtras);
      if (concrete.has_value()) {
        logSeeded(diagnostics, "abstraction: concrete trace found at suffix frame %zu (%.1f s)",
                  concrete->badFrame, elapsedSeconds(abstractionStart));
        return concrete;
      }
      if (!unsat) {
        logSeeded(diagnostics, "abstraction: concretization inconclusive; widening the model");
        if (!widen()) {
          break;
        }
        continue;
      }
      // The abstract trace relied on every lifted pin, so the registers behind
      // all of them join the model, not only those of the failed subset: a
      // failed core is often a single pin and would refine one register per
      // round.
      size_t added = 0;
      std::unordered_set<size_t> failedFrames;
      for (const auto& [symbol, frame, value] : failed) {
        failedFrames.insert(frame);
      }
      for (const auto& [symbol, frame, value] : pins) {
        if (const auto it = indexByRail.find(symbol); it != indexByRail.end()) {
          if (kept.count(divergence.railsByStateIndex[it->second].mayBeOne) == 0) {
            keepRegister(it->second);
            ++added;
          }
        }
      }
      logSeeded(diagnostics,
                "abstraction: pinned trace infeasible; %zu failed pins over %zu frames; the %zu "
                "lifted pins add %zu registers",
                failed.size(), failedFrames.size(), pins.size(), added);
      if (added == 0 && !widen()) {
        break;
      }
    }
    logSeeded(diagnostics, "abstraction: gave up after the iteration budget (%.1f s)",
              elapsedSeconds(abstractionStart));
    return std::nullopt;
  };

  // Seed screening: a short, budgeted bounded search from every end state in
  // rank order. The ranking is a heuristic; a held control input can leave the
  // best-ranked state in an initialization loop from which the target is
  // deep, while another trial ended a few cycles from it. The first state that
  // diverges within the probe depth becomes the start state and its trace the
  // suffix (short searches from several simulated seeds: Ho et al., ICCAD
  // 2000).
  std::optional<FrontierBoundedSearchResult> probedSearch;
  if (options.probeDepth > 0 && trials.size() > 1) {
    const auto probeStart = std::chrono::steady_clock::now();
    const auto heldForTrial = [&](const KInductionProblem& target, const PrefixTrial& trial) {
      const auto lastVector = inputVectorAtCycle(trial, prefixCycles - 1);
      std::vector<std::pair<size_t, bool>> held;
      for (const auto& [index, value] : heldInputIndices) {
        (void)value;
        if (index < target.inputSymbols.size()) {
          held.emplace_back(target.inputSymbols[index], lastVector[index]);
        }
      }
      return held;
    };
    for (size_t rank = 0; rank < ranking.size(); ++rank) {
      const PrefixTrial& trial = trials[ranking[rank]];
      const auto candidate = buildDivergenceProblem(
          model0, model1, alignedInputs, alignedOutputs, trial.state0, trial.state1,
          resetActiveValue);
      if (candidate.frontierDiffers) {
        logSeeded(diagnostics, "probe: trial %zu (rank %zu) already differs in its end state",
                  trial.index, rank);
        best = trial;
        break;
      }
      if (candidate.problem.bad == BoolExpr::createFalse()) {
        continue;
      }
      bool unsat = false;
      BoundedSearchExtras probeExtras;
      probeExtras.unsatOut = &unsat;
      const auto probe = runFrontierBoundedSearch(
          candidate.problem, candidate.frontierFor(trial.state0), {}, options.probeDepth,
          options.checkEvery, noEnlargement, nullptr, solverType, false, nullptr, false,
          options.probeConflicts, heldForTrial(candidate.problem, trial), probeExtras);
      logSeeded(diagnostics,
                "probe: trial %zu (rank %zu, toggled=%zu known=%zu) depth %zu: %s (%.1f s)",
                trial.index, rank, trial.toggledRegisters, trial.knownBits, options.probeDepth,
                probe.has_value() ? "diverges" : unsat ? "unsat" : "budget exhausted",
                elapsedSeconds(probeStart));
      if (probe.has_value()) {
        best = trial;
        probedSearch = probe;
        break;
      }
    }
    if (!probedSearch.has_value()) {
      logSeeded(diagnostics,
                "probe: no end state diverges within %zu frames; deep search from the best-ranked state",
                options.probeDepth);
    }
    // The held controls keep the values of the chosen trial.
    const auto lastVector = inputVectorAtCycle(*best, prefixCycles - 1);
    for (auto& [index, value] : heldInputIndices) {
      value = lastVector[index];
    }
    logSeeded(diagnostics, "probe: start state is trial %zu (%.1f s)", best->index,
              elapsedSeconds(probeStart));
  }

  std::optional<SuffixTrace> trace;
  const auto divergence = buildDivergenceProblem(
      model0, model1, alignedInputs, alignedOutputs, best->state0, best->state1,
      resetActiveValue);
  logSeeded(diagnostics,
            "paired registers=%zu unpaired=%zu differing=%zu incomparable=%zu "
            "differing_outputs=%zu origin_pairs=%zu frontier_differs=%d (%.1f s)",
            divergence.pairedRegisters, divergence.unpairedRegisters,
            divergence.differingRegisters, divergence.incomparableRegisters,
            divergence.differingOutputs, divergence.originPairs,
            divergence.frontierDiffers ? 1 : 0, elapsedSeconds(start));
  if (diagnostics) {
    for (size_t i = 0; i < divergence.unpairedNames.size() && i < 64; ++i) {
      logSeeded(diagnostics, "unpaired design-2 register: %s", divergence.unpairedNames[i].c_str());
    }
    for (const auto& name : divergence.incomparableNames) {
      logSeeded(diagnostics, "incomparable pair (design 1): %s", name.c_str());
    }
    logRegisterHistogram(diagnostics, "differing", model0, divergence.differingIndices, 3, 12,
                         1024);
  }
  if (divergence.frontierDiffers) {
    // The designs already disagree in the selected state: look for the
    // output difference directly on both copies.
    trace = twoCopySearch(best->state0, best->state1, options.searchDepth);
  } else if (divergence.problem.bad == BoolExpr::createFalse()) {
    logSeeded(diagnostics, "no comparable differing next-state logic; search skipped");
    return std::nullopt;
  } else {
    const auto divergenceExact = [&](const std::vector<TernaryValue>& from0,
                                     const std::vector<TernaryValue>& from1,
                                     size_t depth) {
      (void)from1;
      return runFrontierBoundedSearch(
          divergence.problem, divergence.frontierFor(from0), {}, depth,
          options.checkEvery, noEnlargement, nullptr, solverType, diagnostics, nullptr,
          false, 0, heldFor(divergence.problem));
    };
    std::optional<FrontierBoundedSearchResult> search;
    std::vector<std::vector<bool>> chainInputs;
    std::vector<TernaryValue> from0 = best->state0;
    std::vector<TernaryValue> from1 = best->state1;
    if (probedSearch.has_value()) {
      search = probedSearch;
    }
    if (!search.has_value() && options.abstractionLayers > 0) {
      search = runAbstractionGuided(divergence);
    }
    if (!search.has_value() && options.waypointHops > 0) {
      auto outcome = runWaypointChain(divergence);
      if (outcome.divergence.has_value()) {
        search = std::move(outcome.divergence);
      } else {
        // The chain prepared a state closer to the target; the plain search
        // continues from there.
        chainInputs = std::move(outcome.inputs);
        from0 = std::move(outcome.state0);
        from1 = std::move(outcome.state1);
        logSeeded(diagnostics, "waypoints: direct search of %zu frames from the chain's state after %zu frames",
                  options.searchDepth, chainInputs.size());
      }
    }
    if (!search.has_value()) {
      search = runFrontierBoundedSearch(
          divergence.problem, divergence.frontierFor(from0), {}, options.searchDepth,
          options.checkEvery, enlargement, nullptr, solverType, diagnostics, nullptr, false,
          0, heldFor(divergence.problem));
      if (search.has_value()) {
        search = finishEnlarged(*search, from0, from1, divergenceExact);
      }
      if (search.has_value() && !chainInputs.empty()) {
        search->inputsByFrame.insert(search->inputsByFrame.begin(), chainInputs.begin(),
                                     chainInputs.end());
        search->badFrame += chainInputs.size();
      }
    }
    if (!search.has_value()) {
      logSeeded(diagnostics, "no divergence within depth %zu after the prefix (%.1f s)",
                options.searchDepth, elapsedSeconds(start));
      return std::nullopt;
    }
    logSeeded(diagnostics, "designs diverge after suffix frame %zu (%.1f s)",
              search->badFrame, elapsedSeconds(start));
    // The divergence is internal first; let it reach the outputs.
    std::vector<TernaryValue> diverged0;
    std::vector<TernaryValue> diverged1;
    constexpr size_t kExtendTries = 16;
    trace = replayForMismatch(best->state0, best->state1, search->inputsByFrame,
                              options.searchDepth, kExtendTries, diverged0, diverged1);
    if (!trace.has_value() && !diverged0.empty()) {
      constexpr size_t kFollowUpDepth = 8;
      logSeeded(diagnostics, "no output difference by simulation; bounded follow-up of %zu frames",
                kFollowUpDepth);
      if (auto followUp = twoCopySearch(
              diverged0, diverged1, std::min(options.searchDepth, kFollowUpDepth))) {
        SuffixTrace combined;
        combined.inputsByFrame = search->inputsByFrame;
        combined.inputsByFrame.insert(combined.inputsByFrame.end(),
                                      followUp->inputsByFrame.begin(),
                                      followUp->inputsByFrame.end());
        combined.mismatchFrame = search->inputsByFrame.size() + followUp->mismatchFrame;
        combined.mismatches = std::move(followUp->mismatches);
        trace = std::move(combined);
      }
    }
  }
  if (!trace.has_value()) {
    logSeeded(diagnostics, "no confirmed output difference (%.1f s)", elapsedSeconds(start));
    return std::nullopt;
  }

  KInductionResult::CounterexampleWitness witness;
  witness.badFrame = prefixCycles + trace->mismatchFrame;
  witness.outputMismatches = std::move(trace->mismatches);
  const auto inputName = [&](size_t index) {
    return index < problem.environmentInputNames.size()
               ? problem.environmentInputNames[index]
               : alignedInputs.names[index];
  };
  const auto toAssignments = [&](const std::vector<bool>& vector) {
    std::vector<KInductionResult::FrameInputAssignments::Assignment> assignments;
    assignments.reserve(inputCount);
    for (size_t index = 0; index < inputCount; ++index) {
      assignments.push_back({inputName(index), vector[index]});
    }
    return assignments;
  };
  for (size_t cycle = 0; cycle < prefixCycles; ++cycle) {
    const auto vector = inputVectorAtCycle(*best, cycle);
    if (!witness.heldInputSegments.empty() &&
        witness.heldInputSegments.back().values == vector) {
      witness.heldInputSegments.back().lastFrame = cycle;
      continue;
    }
    KInductionResult::HeldInputSegment segment;
    segment.firstFrame = cycle;
    segment.lastFrame = cycle;
    segment.values = vector;
    segment.assignments = toAssignments(vector);
    witness.heldInputSegments.push_back(std::move(segment));
  }
  for (size_t frame = 0; frame <= trace->mismatchFrame; ++frame) {
    KInductionResult::FrameInputAssignments frameInputs;
    frameInputs.frame = prefixCycles + frame;
    frameInputs.assignments = toAssignments(trace->inputsByFrame[frame]);
    witness.inputTrace.push_back(std::move(frameInputs));
  }
  logSeeded(diagnostics, "confirmed difference at cycle %zu (%.1f s)",
            witness.badFrame, elapsedSeconds(start));
  return witness;
}

}  // namespace KEPLER_FORMAL::SEC
