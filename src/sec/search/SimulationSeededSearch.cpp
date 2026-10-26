// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "search/SimulationSeededSearch.h"

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
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
    changed = changed || next != state_[index];
    state_[index] = next;
  }
  combinationalValid_ = false;
  return changed;
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
};

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
  std::vector<std::pair<size_t, size_t>> pairs;
  std::vector<bool> used1(model1.stateBits.size(), false);
  std::vector<bool> used0(model0.stateBits.size(), false);
  for (const bool byNetName : {true, false}) {
    for (size_t index0 = 0; index0 < model0.stateBits.size(); ++index0) {
      if (used0[index0]) {
        continue;
      }
      const auto& key0 = model0.stateBits[index0];
      size_t index1 = kNoState;
      if (byNetName) {
        if (const auto it = model0.outputNetNameByKey.find(key0);
            it != model0.outputNetNameByKey.end() && !it->second.empty()) {
          index1 = lookup(indexByNetName1, it->second);
        }
      } else if (const auto it = model0.displayNameByKey.find(key0);
                 it != model0.displayNameByKey.end()) {
        index1 = lookup(indexByDisplayName1, it->second);
      }
      if (index1 == kNoState || used1[index1]) {
        continue;
      }
      used1[index1] = true;
      used0[index0] = true;
      pairs.emplace_back(index0, index1);
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
    diverge = BoolExpr::Or(diverge, binaryDifferent(lift(next0It->second), lift(renamed)));
  }
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
    KEPLER_FORMAL::Config::SolverType solverType,
    bool diagnostics) {
  const auto start = std::chrono::steady_clock::now();
  const bool finalCheckOnly = checkEvery == 0;
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

  // One final check is a single query, so the configured solver serves it;
  // intermediate checks need assumptions and an incremental solver.
  const auto searchSolverType = finalCheckOnly
      ? solverType
      : SATSolverWrapper::assumptionSolverTypeFor(solverType);
  const bool usesAssumptions =
      searchSolverType != KEPLER_FORMAL::Config::SolverType::KISSAT;
  SATSolverWrapper solver(searchSolverType);
  solver.configureForSecIncrementalBoundedSearch();
  std::vector<size_t> symbols = problem.inputSymbols;
  symbols.insert(symbols.end(), stateSymbols.begin(), stateSymbols.end());
  FrameVariableStore variables(solver, symbols, depth + 1);
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

  const auto addFrameDomain = [&](size_t frame) {
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
    // The simulated prefix is at least as long as the reset bootstrap, so
    // reset controls stay released, as in every post-bootstrap frame.
    for (const auto& [symbol, activeValue] : problem.resetBootstrapInputs) {
      if (variables.hasSymbol(symbol)) {
        const int literal = variables.getLiteral(symbol, frame);
        solver.addClause({activeValue ? -literal : literal});
      }
    }
  };
  addFrameDomain(0);

  const auto extractInputs = [&](size_t badFrame) {
    FrontierBoundedSearchResult result;
    result.badFrame = badFrame;
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

  std::vector<int> pendingBadLits;
  size_t firstPendingFrame = 0;
  // The selected state itself may already show the difference.
  {
    FrameFormulaEncoder badEncoder(solver, variables.makeLeafLits(0));
    const int badLit = badEncoder.encode(bad);
    if (finalCheckOnly) {
      pendingBadLits.push_back(badLit);
    } else {
      const auto status = solver.solveWithAssumptionsStatus({badLit});
      if (status == SATSolverWrapper::SolveStatus::Sat) {
        logSeeded(diagnostics, "difference at the selected state itself");
        return extractInputs(0);
      }
      if (status != SATSolverWrapper::SolveStatus::Unsat) {
        return std::nullopt;  // LCOV_EXCL_LINE
      }
      solver.addClause({-badLit});
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
    events.reserve(depth);
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

  const auto checkPendingBad = [&](size_t frame) -> std::optional<FrontierBoundedSearchResult> {
    int anyBad = rails.falseLit();
    for (const int lit : pendingBadLits) {
      anyBad = rails.orLit(anyBad, lit);
    }
    const auto solveStart = std::chrono::steady_clock::now();
    SATSolverWrapper::SolveStatus status = SATSolverWrapper::SolveStatus::Unknown;
    if (usesAssumptions) {
      status = solver.solveWithAssumptionsStatus({anyBad});
    } else {
      solver.addClause({anyBad});
      status = solver.solveStatus();
    }
    logSeeded(diagnostics,
              "frame=%zu checked_frames=%zu status=%s solve_seconds=%.1f total_seconds=%.1f",
              frame,
              pendingBadLits.size(),
              status == SATSolverWrapper::SolveStatus::Sat
                  ? "sat"
                  : status == SATSolverWrapper::SolveStatus::Unsat ? "unsat" : "unknown",
              elapsedSeconds(solveStart),
              elapsedSeconds(start));
    if (status == SATSolverWrapper::SolveStatus::Sat) {
      // The earliest checked frame whose property fails.
      for (size_t i = 0; i < pendingBadLits.size(); ++i) {
        if (solver.getLiteralValue(pendingBadLits[i])) {
          return extractInputs(firstPendingFrame + i);
        }
      }
      return extractInputs(frame);  // LCOV_EXCL_LINE
    }
    if (status == SATSolverWrapper::SolveStatus::Unsat && usesAssumptions) {
      for (const int lit : pendingBadLits) {
        solver.addClause({-lit});
      }
    }
    pendingBadLits.clear();
    firstPendingFrame = frame + 1;
    return std::nullopt;
  };

  for (size_t frame = 0; frame < depth; ++frame) {
    const auto frameStart = std::chrono::steady_clock::now();
    FrameFormulaEncoder encoder(solver, variables.makeLeafLits(frame));
    for (const auto symbol : transitionStates) {
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
              !variables.hasSymbol(port.outputRails[bit]->mayBeOne)) {
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

    addFrameDomain(frame + 1);
    FrameFormulaEncoder badEncoder(solver, variables.makeLeafLits(frame + 1));
    pendingBadLits.push_back(badEncoder.encode(bad));
    logSeeded(diagnostics, "frame=%zu encoded encode_seconds=%.1f frontier_units=%zu",
              frame + 1, elapsedSeconds(frameStart), frontierUnits);
    // Every intermediate check is an unsatisfiability proof the search does
    // not need, so the property is checked for the first frames, then every
    // few frames, and at the last frame.
    const size_t checked = frame + 1;
    if (checked == depth ||
        (!finalCheckOnly && (checked <= 2 || checked % checkEvery == 0))) {
      if (auto result = checkPendingBad(checked)) {
        return result;
      }
    }
  }
  return std::nullopt;
}

struct PrefixTrial {
  size_t knownBits = 0;
  // Cycles in which some state bit changed: a tie-breaker that prefers a
  // design that kept running over one held in reset.
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
    for (size_t cycle = 0; cycle < prefixCycles; ++cycle) {
      applyInputs(inputVectorAtCycle(trial, cycle));
      const bool changed0 = simulator0.step();
      const bool changed1 = simulator1.step();
      trial.activeCycles += (changed0 || changed1) ? 1 : 0;
    }
    trial.knownBits = simulator0.knownStateBitCount() + simulator1.knownStateBitCount();
    logSeeded(diagnostics, "trial=%zu known_state_bits=%zu active_cycles=%zu seconds=%.1f",
              trialIndex, trial.knownBits, trial.activeCycles, elapsedSeconds(start));
    const bool better = !best.has_value() || trial.knownBits > best->knownBits ||
                        (trial.knownBits == best->knownBits &&
                         trial.activeCycles > best->activeCycles);
    if (better) {
      trial.state0 = simulator0.state();
      trial.state1 = simulator1.state();
      best = std::move(trial);
    }
  }
  if (!best.has_value()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  logSeeded(diagnostics,
            "prefix_cycles=%zu best_known_state_bits=%zu of %zu active_cycles=%zu",
            prefixCycles, best->knownBits,
            simulator0.stateCount() + simulator1.stateCount(), best->activeCycles);

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
  // Two-copy bounded search from a pair of states, confirmed by replay.
  const auto twoCopySearch = [&](const std::vector<TernaryValue>& from0,
                                 const std::vector<TernaryValue>& from1,
                                 size_t depth) -> std::optional<SuffixTrace> {
    const auto search = runFrontierBoundedSearch(
        problem, frontierOf(from0, from1), {}, depth, options.checkEvery, solverType,
        diagnostics);
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

  std::optional<SuffixTrace> trace;
  const auto divergence = buildDivergenceProblem(
      model0, model1, alignedInputs, alignedOutputs, best->state0, best->state1,
      resetActiveValue);
  logSeeded(diagnostics,
            "paired registers=%zu unpaired=%zu differing=%zu incomparable=%zu "
            "differing_outputs=%zu frontier_differs=%d (%.1f s)",
            divergence.pairedRegisters, divergence.unpairedRegisters,
            divergence.differingRegisters, divergence.incomparableRegisters,
            divergence.differingOutputs, divergence.frontierDiffers ? 1 : 0,
            elapsedSeconds(start));
  if (diagnostics) {
    for (size_t i = 0; i < divergence.unpairedNames.size() && i < 64; ++i) {
      logSeeded(diagnostics, "unpaired design-2 register: %s", divergence.unpairedNames[i].c_str());
    }
    for (const auto& name : divergence.incomparableNames) {
      logSeeded(diagnostics, "incomparable pair (design 1): %s", name.c_str());
    }
  }
  if (divergence.frontierDiffers) {
    // The designs already disagree in the selected state: look for the
    // output difference directly on both copies.
    trace = twoCopySearch(best->state0, best->state1, options.searchDepth);
  } else if (divergence.problem.bad == BoolExpr::createFalse()) {
    logSeeded(diagnostics, "no comparable differing next-state logic; search skipped");
    return std::nullopt;
  } else {
    const auto search = runFrontierBoundedSearch(
        divergence.problem, divergence.frontier, {}, options.searchDepth,
        options.checkEvery, solverType, diagnostics);
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
