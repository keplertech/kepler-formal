// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "phase/NormalizedTransitionSystem.h"

#include <algorithm>
#include <functional>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

#include "BoolExpr.h"
#include "DNL.h"
#include "NLDB.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLPath.h"
#include "BuildPrimaryOutputClauses.h"

namespace KEPLER_FORMAL::SEC::PHASE {

namespace {

using BooleanExpression = naja::NL::SNLDesignModeling::BooleanExpression;
using SequentialModel = naja::NL::SNLDesignModeling::SequentialModel;
using DNLID = naja::DNL::DNLID;

std::mutex& frontendMutex() {
  static std::mutex mutex;
  return mutex;
}

bool isConstant(const ExpressionArena& arena, ExprID id, bool value) {
  if (!arena.isValid(id)) {
    return false;
  }
  const auto& node = arena.node(id);
  return node.op == ExpressionOp::Constant && node.constant == value;
}

bool isNegationPair(const ExpressionArena& arena, ExprID lhs, ExprID rhs) {
  if (!arena.isValid(lhs) || !arena.isValid(rhs)) {
    return false;
  }
  const auto& lhsNode = arena.node(lhs);
  const auto& rhsNode = arena.node(rhs);
  return (lhsNode.op == ExpressionOp::Not && lhsNode.left == rhs) ||
         (rhsNode.op == ExpressionOp::Not && rhsNode.left == lhs);
}

SignalKey terminalKey(const naja::DNL::DNLTerminalFull& terminal) {
  SignalKey key;
  const auto names = terminal.getDNLInstance().getPath().getPathNames();
  key.first.reserve(names.size() + 1);
  for (const auto& name : names) {
    key.first.push_back(name.getID());
  }
  key.first.push_back(terminal.getSnlBitTerm()->getName().getID());
  key.second.push_back(static_cast<naja::NL::NLID::DesignObjectID>(
      terminal.getSnlBitTerm()->getBit()));
  return key;
}

std::string terminalName(const naja::DNL::DNLTerminalFull& terminal) {
  std::ostringstream out;
  const auto names = terminal.getDNLInstance().getPath().getPathNames();
  for (const auto& name : names) {
    out << name.getString() << ".";
  }
  out << terminal.getSnlBitTerm()->getName().getString() << "["
      << terminal.getSnlBitTerm()->getBit() << "]";
  return out.str();
}

std::string instanceName(const naja::DNL::DNLInstanceFull& instance) {
  std::string name = instance.getFullPath();
  while (!name.empty() && (name.back() == '.' || name.back() == '/')) {
    name.pop_back();
  }
  if (name.empty()) {
    name = "<instance#" + std::to_string(instance.getID()) + ">";
  }
  return name;
}

std::optional<std::pair<size_t, bool>> stateOutputReference(
    const BooleanExpression& expression) {
  using Operator = BooleanExpression::Operator;
  if (!expression.isValid()) {
    return std::nullopt;
  }
  const auto* node = &expression.nodes[expression.root];
  bool complemented = false;
  if (node->operation == Operator::Not && node->operands.size() == 1 &&
      node->operands.front() < expression.nodes.size()) {
    complemented = true;
    node = &expression.nodes[node->operands.front()];
  }
  if (node->operation != Operator::State) {
    return std::nullopt;
  }
  return std::pair<size_t, bool>{node->state, complemented};
}

void collectTermReferences(
    const BooleanExpression& expression,
    std::unordered_set<const naja::NL::SNLBitTerm*>& terms) {
  using Operator = BooleanExpression::Operator;
  for (const auto& node : expression.nodes) {
    if (node.operation == Operator::Term && node.term != nullptr) {
      terms.insert(node.term);
    }
  }
}

struct PendingState {
  size_t instanceIndex = 0;
  size_t modelStateIndex = 0;
  bool primaryOutputComplemented = false;
  SignalKey key;
  std::string displayName;
  size_t systemStateIndex = std::numeric_limits<size_t>::max();
};

struct PendingInstance {
  const SequentialModel* model = nullptr;
  std::string path;
  std::string modelName;
  std::unordered_map<const naja::NL::SNLBitTerm*, DNLID> termIDs;
  std::vector<std::optional<std::pair<DNLID, bool>>> primaryOutputs;
  std::vector<size_t> systemStateByModelState;
  std::vector<std::pair<DNLID, std::pair<size_t, bool>>> modeledOutputs;
};

struct FrontendCleanup {
  naja::NL::NLUniverse* universe = nullptr;
  naja::NL::NLDB* previousTopDB = nullptr;
  naja::NL::NLDB* selectedDB = nullptr;
  naja::NL::SNLDesign* previousSelectedDBTop = nullptr;

  ~FrontendCleanup() {
    naja::DNL::destroy();
    if (selectedDB != nullptr) {
      selectedDB->setTopDesign(previousSelectedDBTop);
    }
    if (universe != nullptr) {
      universe->setTopDB(previousTopDB);
    }
  }
};

ClearPresetValue convertClearPresetValue(
    naja::NL::SNLDesignModeling::SequentialState::ClearPresetValue value) {
  using Source = naja::NL::SNLDesignModeling::SequentialState::ClearPresetValue;
  switch (value) {
    case Source::Zero:
      return ClearPresetValue::Zero;
    case Source::One:
      return ClearPresetValue::One;
    case Source::Hold:
      return ClearPresetValue::Hold;
    case Source::Toggle:
      return ClearPresetValue::Toggle;
    case Source::Unknown:
      return ClearPresetValue::Unknown;
  }
  return ClearPresetValue::Unknown;
}

TernaryValue invertTernary(TernaryValue value) {
  switch (value) {
    case TernaryValue::Zero:
      return TernaryValue::One;
    case TernaryValue::One:
      return TernaryValue::Zero;
    case TernaryValue::Unknown:
      return TernaryValue::Unknown;
  }
  return TernaryValue::Unknown;
}

ExprID compileBoolExpr(
    BoolExpr* root,
    const std::unordered_map<size_t, ExprID>& sourceVariableExpressions,
    ExpressionArena& arena, std::unordered_map<BoolExpr*, ExprID>& memo,
    std::string* error) {
  if (root == nullptr || !root->isValid()) {
    if (error != nullptr) {
      *error = "invalid Boolean expression";
    }
    return InvalidExpr;
  }

  struct Frame {
    BoolExpr* expression = nullptr;
    bool visited = false;
  };
  std::vector<Frame> stack{{root, false}};
  while (!stack.empty()) {
    const Frame frame = stack.back();
    stack.pop_back();
    BoolExpr* expression = frame.expression;
    if (memo.find(expression) != memo.end()) {
      continue;
    }
    if (expression->getOp() == Op::VAR) {
      if (expression->getId() == 0) {
        memo.emplace(expression, arena.constant(false));
      } else if (expression->getId() == 1) {
        memo.emplace(expression, arena.constant(true));
      } else {
        const auto variable =
            sourceVariableExpressions.find(expression->getId());
        if (variable == sourceVariableExpressions.end()) {
          if (error != nullptr) {
            *error = "unpublished Boolean variable " +
                     std::to_string(expression->getId());
          }
          return InvalidExpr;
        }
        memo.emplace(expression, variable->second);
      }
      continue;
    }
    if (expression->getOp() == Op::NONE) {
      if (error != nullptr) {
        *error = "unsupported Boolean expression operator";
      }
      return InvalidExpr;
    }
    if (!frame.visited) {
      stack.push_back({expression, true});
      if (expression->getRight() != nullptr &&
          memo.find(expression->getRight()) == memo.end()) {
        stack.push_back({expression->getRight(), false});
      }
      if (expression->getLeft() != nullptr &&
          memo.find(expression->getLeft()) == memo.end()) {
        stack.push_back({expression->getLeft(), false});
      }
      continue;
    }

    const ExprID left = memo.at(expression->getLeft());
    ExprID result = InvalidExpr;
    switch (expression->getOp()) {
      case Op::NOT:
        result = arena.makeNot(left);
        break;
      case Op::AND:
        result = arena.makeAnd(left, memo.at(expression->getRight()));
        break;
      case Op::OR:
        result = arena.makeOr(left, memo.at(expression->getRight()));
        break;
      case Op::XOR:
        result = arena.makeXor(left, memo.at(expression->getRight()));
        break;
      case Op::VAR:
      case Op::NONE:
      default:
        if (error != nullptr) {
          *error = "unsupported Boolean expression operator";
        }
        return InvalidExpr;
    }
    memo.emplace(expression, result);
  }
  return memo.at(root);
}

}  // namespace

ExpressionArena::ExpressionArena() {
  nodes_.push_back(
      {ExpressionOp::Constant, 0, false, InvalidExpr, InvalidExpr});
  nodes_.push_back({ExpressionOp::Constant, 0, true, InvalidExpr, InvalidExpr});
}

ExprID ExpressionArena::variable(size_t variableIndex) {
  if (const auto found = variableNodes_.find(variableIndex);
      found != variableNodes_.end()) {
    return found->second;
  }
  if (nodes_.size() >= InvalidExpr) {
    throw std::overflow_error("phase expression arena exhausted");
  }
  const ExprID id = static_cast<ExprID>(nodes_.size());
  nodes_.push_back(
      {ExpressionOp::Variable, variableIndex, false, InvalidExpr, InvalidExpr});
  variableNodes_.emplace(variableIndex, id);
  return id;
}

ExprID ExpressionArena::makeNot(ExprID operand) {
  if (!isValid(operand)) {
    return InvalidExpr;
  }
  if (isConstant(*this, operand, false)) {
    return constant(true);
  }
  if (isConstant(*this, operand, true)) {
    return constant(false);
  }
  if (node(operand).op == ExpressionOp::Not) {
    return node(operand).left;
  }
  if (nodes_.size() >= InvalidExpr) {
    throw std::overflow_error("phase expression arena exhausted");
  }
  const ExprID id = static_cast<ExprID>(nodes_.size());
  nodes_.push_back({ExpressionOp::Not, 0, false, operand, InvalidExpr});
  return id;
}

ExprID ExpressionArena::appendBinary(ExpressionOp op, ExprID lhs, ExprID rhs) {
  if (!isValid(lhs) || !isValid(rhs)) {
    return InvalidExpr;
  }
  if (rhs < lhs) {
    std::swap(lhs, rhs);
  }
  if (nodes_.size() >= InvalidExpr) {
    throw std::overflow_error("phase expression arena exhausted");
  }
  const ExprID id = static_cast<ExprID>(nodes_.size());
  nodes_.push_back({op, 0, false, lhs, rhs});
  return id;
}

ExprID ExpressionArena::makeAnd(ExprID lhs, ExprID rhs) {
  if (!isValid(lhs) || !isValid(rhs)) {
    return InvalidExpr;
  }
  if (isConstant(*this, lhs, false) || isConstant(*this, rhs, false)) {
    return constant(false);
  }
  if (isConstant(*this, lhs, true)) {
    return rhs;
  }
  if (isConstant(*this, rhs, true)) {
    return lhs;
  }
  if (lhs == rhs) {
    return lhs;
  }
  if (isNegationPair(*this, lhs, rhs)) {
    return constant(false);
  }
  return appendBinary(ExpressionOp::And, lhs, rhs);
}

ExprID ExpressionArena::makeOr(ExprID lhs, ExprID rhs) {
  if (!isValid(lhs) || !isValid(rhs)) {
    return InvalidExpr;
  }
  if (isConstant(*this, lhs, true) || isConstant(*this, rhs, true)) {
    return constant(true);
  }
  if (isConstant(*this, lhs, false)) {
    return rhs;
  }
  if (isConstant(*this, rhs, false)) {
    return lhs;
  }
  if (lhs == rhs) {
    return lhs;
  }
  if (isNegationPair(*this, lhs, rhs)) {
    return constant(true);
  }
  return appendBinary(ExpressionOp::Or, lhs, rhs);
}

ExprID ExpressionArena::makeXor(ExprID lhs, ExprID rhs) {
  if (!isValid(lhs) || !isValid(rhs)) {
    return InvalidExpr;
  }
  if (isConstant(*this, lhs, false)) {
    return rhs;
  }
  if (isConstant(*this, rhs, false)) {
    return lhs;
  }
  if (isConstant(*this, lhs, true)) {
    return makeNot(rhs);
  }
  if (isConstant(*this, rhs, true)) {
    return makeNot(lhs);
  }
  if (lhs == rhs) {
    return constant(false);
  }
  if (isNegationPair(*this, lhs, rhs)) {
    return constant(true);
  }
  return appendBinary(ExpressionOp::Xor, lhs, rhs);
}

ExprID ExpressionArena::makeIte(ExprID condition, ExprID whenTrue,
                                ExprID whenFalse) {
  if (!isValid(condition) || !isValid(whenTrue) || !isValid(whenFalse)) {
    return InvalidExpr;
  }
  if (isConstant(*this, condition, true)) {
    return whenTrue;
  }
  if (isConstant(*this, condition, false)) {
    return whenFalse;
  }
  if (whenTrue == whenFalse) {
    return whenTrue;
  }
  return makeOr(makeAnd(condition, whenTrue),
                makeAnd(makeNot(condition), whenFalse));
}

size_t NormalizedTransitionSystem::addInput(SignalKey key,
                                            std::string displayName) {
  const size_t variable = variables_.size();
  variables_.push_back(
      {std::move(key), std::move(displayName), VariableKind::Input});
  inputVariables_.push_back(variable);
  return variable;
}

size_t NormalizedTransitionSystem::addState(SignalKey key,
                                            std::string displayName,
                                            TernaryValue initialValue) {
  const size_t variable = variables_.size();
  variables_.push_back(
      {std::move(key), std::move(displayName), VariableKind::State});
  const size_t state = states_.size();
  states_.push_back({variable, initialValue, InvalidExpr});
  return state;
}

void NormalizedTransitionSystem::setNextState(size_t stateIndex,
                                              ExprID expression) {
  states_.at(stateIndex).nextState = expression;
}

size_t NormalizedTransitionSystem::addLatch(NormalizedLatch latch) {
  const size_t index = latches_.size();
  latches_.push_back(std::move(latch));
  return index;
}

size_t NormalizedTransitionSystem::addOutput(NormalizedOutput output) {
  const size_t index = outputs_.size();
  outputs_.push_back(std::move(output));
  return index;
}

void NormalizedTransitionSystem::addDiagnostic(
    NormalizationDiagnostic diagnostic) {
  diagnostics_.push_back(std::move(diagnostic));
}

const NormalizedVariable& NormalizedTransitionSystem::stateVariable(
    size_t stateIndex) const {
  return variables_.at(states_.at(stateIndex).variable);
}

bool NormalizedTransitionSystem::hasFatalDiagnostics() const {
  return std::any_of(diagnostics_.begin(), diagnostics_.end(),
                     [](const auto& diagnostic) { return diagnostic.fatal; });
}

bool NormalizedTransitionSystem::validate(std::string* reason) const {
  auto fail = [&](std::string message) {
    if (reason != nullptr) {
      *reason = std::move(message);
    }
    return false;
  };
  if (expressions_.nodes().size() < 2 || !isConstant(expressions_, 0, false) ||
      !isConstant(expressions_, 1, true)) {
    return fail("expression DAG has invalid constant nodes");
  }
  for (size_t nodeIndex = 0; nodeIndex < expressions_.nodes().size();
       ++nodeIndex) {
    const auto& node = expressions_.nodes()[nodeIndex];
    switch (node.op) {
      case ExpressionOp::Constant:
        break;
      case ExpressionOp::Variable:
        if (node.variable >= variables_.size()) {
          return fail("expression references an out-of-range variable");
        }
        break;
      case ExpressionOp::Not:
        if (node.left >= nodeIndex) {
          return fail("expression DAG is not topologically ordered");
        }
        break;
      case ExpressionOp::And:
      case ExpressionOp::Or:
      case ExpressionOp::Xor:
        if (node.left >= nodeIndex || node.right >= nodeIndex) {
          return fail("expression DAG is not topologically ordered");
        }
        break;
    }
  }
  std::vector<uint8_t> variableUse(variables_.size(), 0);
  for (const size_t variable : inputVariables_) {
    if (variable >= variables_.size() ||
        variables_[variable].kind != VariableKind::Input) {
      return fail("input list references an invalid input variable");
    }
    if (variableUse[variable] != 0) {
      return fail("input list references a variable more than once");
    }
    variableUse[variable] = 1;
  }
  for (const auto& state : states_) {
    if (state.variable >= variables_.size() ||
        variables_[state.variable].kind != VariableKind::State) {
      return fail("state references an invalid state variable");
    }
    if (variableUse[state.variable] != 0) {
      return fail("state list references a variable more than once");
    }
    variableUse[state.variable] = 1;
    if (!expressions_.isValid(state.nextState)) {
      return fail("state has no valid next-state expression");
    }
  }
  if (std::find(variableUse.begin(), variableUse.end(), uint8_t{0}) !=
      variableUse.end()) {
    return fail("variable is missing from the input or state index");
  }
  std::vector<uint8_t> latchByState(states_.size(), 0);
  for (const auto& latch : latches_) {
    if (latch.stateIndex >= states_.size() ||
        !expressions_.isValid(latch.data) ||
        !expressions_.isValid(latch.transparency)) {
      return fail("latch has incomplete normalized semantics");
    }
    if ((latch.clear.has_value() && !expressions_.isValid(*latch.clear)) ||
        (latch.preset.has_value() && !expressions_.isValid(*latch.preset))) {
      return fail("latch has an invalid clear or preset expression");
    }
    if (latchByState[latch.stateIndex] != 0) {
      return fail("state is described by more than one latch");
    }
    latchByState[latch.stateIndex] = 1;
  }
  for (const auto& output : outputs_) {
    if (!expressions_.isValid(output.expression)) {
      return fail("observed output has no valid expression");
    }
  }
  return true;
}

const char* toString(TernaryValue value) {
  switch (value) {
    case TernaryValue::Zero:
      return "0";
    case TernaryValue::One:
      return "1";
    case TernaryValue::Unknown:
      return "X";
  }
  return "X";
}

NormalizedTransitionSystem collectNormalizedTransitionSystem(
    naja::NL::SNLDesign* top, const NormalizationOptions& options) {
  if (top == nullptr) {
    throw std::invalid_argument(
        "collectNormalizedTransitionSystem: null top design");
  }
  std::lock_guard<std::mutex> frontendLock(frontendMutex());
  auto* universe = naja::NL::NLUniverse::get();
  if (universe == nullptr) {
    throw std::runtime_error(
        "collectNormalizedTransitionSystem: NLUniverse is not available");
  }
  auto* selectedDB = top->getDB();
  FrontendCleanup cleanup{
      universe, universe->getTopDB(), selectedDB,
      selectedDB == nullptr ? nullptr : selectedDB->getTopDesign()};
  naja::DNL::destroy();
  universe->setTopDesign(top);

  NormalizedTransitionSystem system;
  KEPLER_FORMAL::BuildPrimaryOutputClauses boundary;
  boundary.setRetainDnl(true);
  boundary.collect();
  auto* dnl = naja::DNL::get();

  std::vector<PendingInstance> instances;
  std::vector<PendingState> states;
  std::unordered_set<DNLID> requestedTerms;

  for (const DNLID leaf : dnl->getLeaves()) {
    const auto& instance = dnl->getDNLInstanceFromID(leaf);
    auto* modelDesign = instance.getSNLModel();
    if (modelDesign == nullptr) {
      continue;
    }
    if (naja::NL::SNLDesignModeling::hasMemoryInterface(modelDesign)) {
      system.addDiagnostic(
          {instanceName(instance),
           "structured memories are not yet normalized by phase analysis",
           true});
      continue;
    }
    if (!naja::NL::SNLDesignModeling::hasSequentialModel(modelDesign)) {
      bool hasSequentialArc = false;
      for (DNLID termID = instance.getTermIndexes().first;
           termID != naja::DNL::DNLID_MAX &&
           termID <= instance.getTermIndexes().second;
           ++termID) {
        const auto& term = dnl->getDNLTerminalFromID(termID);
        if (!naja::NL::SNLDesignModeling::getOutputRelatedClocks(
                 term.getSnlBitTerm())
                 .empty()) {
          hasSequentialArc = true;
          break;
        }
      }
      if (hasSequentialArc) {
        system.addDiagnostic(
            {instanceName(instance),
             "sequential primitive has no usable Naja sequential model", true});
      }
      continue;
    }

    const auto& model =
        naja::NL::SNLDesignModeling::getSequentialModel(modelDesign);
    if (!model.isValid()) {
      system.addDiagnostic(
          {instanceName(instance), "invalid Naja sequential model", true});
      continue;
    }

    PendingInstance pending;
    pending.model = &model;
    pending.path = instanceName(instance);
    pending.modelName = modelDesign->getName().getString();
    pending.primaryOutputs.resize(model.states.size());
    pending.systemStateByModelState.assign(model.states.size(),
                                           std::numeric_limits<size_t>::max());
    for (DNLID termID = instance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= instance.getTermIndexes().second;
         ++termID) {
      const auto& term = dnl->getDNLTerminalFromID(termID);
      if (!term.isNull()) {
        pending.termIDs.emplace(term.getSnlBitTerm(), termID);
      }
    }

    for (const auto& output : model.outputs) {
      const auto term = pending.termIDs.find(output.term);
      const auto reference = stateOutputReference(output.function);
      if (term == pending.termIDs.end() || !reference.has_value() ||
          reference->first >= model.states.size()) {
        system.addDiagnostic(
            {pending.path,
             "sequential output is not a direct state or complemented-state "
             "function",
             true});
        continue;
      }
      pending.modeledOutputs.push_back(
          {term->second, {reference->first, reference->second}});
      auto& primary = pending.primaryOutputs[reference->first];
      if (!primary.has_value() || (primary->second && !reference->second)) {
        primary = std::pair<DNLID, bool>{term->second, reference->second};
      }
    }

    const size_t pendingInstanceIndex = instances.size();
    instances.push_back(std::move(pending));
    auto& stored = instances.back();
    for (size_t stateIndex = 0; stateIndex < model.states.size();
         ++stateIndex) {
      if (!stored.primaryOutputs[stateIndex].has_value()) {
        system.addDiagnostic({stored.path,
                              "sequential state " + std::to_string(stateIndex) +
                                  " has no modeled physical output",
                              true});
        continue;
      }
      const auto [termID, complemented] = *stored.primaryOutputs[stateIndex];
      const auto& term = dnl->getDNLTerminalFromID(termID);
      states.push_back({pendingInstanceIndex, stateIndex, complemented,
                        terminalKey(term), terminalName(term)});
    }

    std::unordered_set<const naja::NL::SNLBitTerm*> referencedModelTerms;
    if (model.kind == SequentialModel::Kind::Latch) {
      collectTermReferences(model.clockedOn, referencedModelTerms);
    }
    for (const auto& state : model.states) {
      collectTermReferences(state.nextState, referencedModelTerms);
      if (state.clear.has_value()) {
        collectTermReferences(*state.clear, referencedModelTerms);
      }
      if (state.preset.has_value()) {
        collectTermReferences(*state.preset, referencedModelTerms);
      }
    }
    for (const auto* modelTerm : referencedModelTerms) {
      const auto physical = stored.termIDs.find(modelTerm);
      if (physical == stored.termIDs.end()) {
        system.addDiagnostic(
            {stored.path,
             "sequential expression references an unmapped model terminal",
             true});
        continue;
      }
      requestedTerms.insert(physical->second);
    }
  }

  std::sort(states.begin(), states.end(), [](const auto& lhs, const auto& rhs) {
    if (SignalKeyLess{}(lhs.key, rhs.key)) {
      return true;
    }
    if (SignalKeyLess{}(rhs.key, lhs.key)) {
      return false;
    }
    return lhs.displayName < rhs.displayName;
  });
  for (auto& state : states) {
    TernaryValue initial = TernaryValue::Unknown;
    if (const auto override =
            options.initialValueByDisplayName.find(state.displayName);
        override != options.initialValueByDisplayName.end()) {
      initial = override->second;
      // The transition variable represents the Naja model state, while the
      // only physical output available for naming may be its complement.
      if (state.primaryOutputComplemented) {
        initial = invertTernary(initial);
      }
    }
    state.systemStateIndex =
        system.addState(state.key, state.displayName, initial);
    instances[state.instanceIndex]
        .systemStateByModelState[state.modelStateIndex] =
        state.systemStateIndex;
  }

  std::vector<DNLID> topOutputTerms;
  if (options.collectObservedOutputs) {
    const auto& topInstance = dnl->getTop();
    for (DNLID termID = topInstance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= topInstance.getTermIndexes().second;
         ++termID) {
      const auto& term = dnl->getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() !=
          naja::NL::SNLBitTerm::Direction::Input) {
        topOutputTerms.push_back(termID);
        requestedTerms.insert(termID);
      }
    }
  }

  std::vector<DNLID> materializationInputs = boundary.getInputs();
  for (const auto& pending : instances) {
    for (const auto& output : pending.modeledOutputs) {
      materializationInputs.push_back(output.first);
    }
  }
  std::sort(materializationInputs.begin(), materializationInputs.end());
  materializationInputs.erase(
      std::unique(materializationInputs.begin(), materializationInputs.end()),
      materializationInputs.end());
  std::vector<DNLID> materializationOutputs(requestedTerms.begin(),
                                            requestedTerms.end());
  std::sort(materializationOutputs.begin(), materializationOutputs.end());

  KEPLER_FORMAL::BuildPrimaryOutputClauses materializer;
  materializer.setRetainDnl(true);
  materializer.setStopAtOpaqueInternalOutputs(true);
  materializer.setInputs(materializationInputs);
  materializer.setOutputs(materializationOutputs);
  materializer.build();

  const auto& sourceTermVars = materializer.getTermDNLID2VarID();
  std::unordered_map<size_t, ExprID> sourceVariableExpressions;
  for (const auto& pending : instances) {
    for (const auto& output : pending.modeledOutputs) {
      const DNLID termID = output.first;
      const size_t modelState = output.second.first;
      const bool complemented = output.second.second;
      if (modelState >= pending.systemStateByModelState.size()) {
        continue;
      }
      const size_t systemState = pending.systemStateByModelState[modelState];
      if (systemState == std::numeric_limits<size_t>::max() ||
          termID >= sourceTermVars.size() || sourceTermVars[termID] < 2) {
        continue;
      }
      ExprID expression =
          system.variableExpression(system.states()[systemState].variable);
      if (complemented) {
        expression = system.expressions().makeNot(expression);
      }
      sourceVariableExpressions.insert_or_assign(sourceTermVars[termID],
                                                 expression);
    }
  }

  struct PendingInput {
    DNLID termID = naja::DNL::DNLID_MAX;
    size_t sourceVariable = 0;
    SignalKey key;
    std::string displayName;
  };
  std::vector<PendingInput> pendingInputs;
  std::unordered_set<size_t> seenSourceVariables;
  for (const DNLID termID : materializationInputs) {
    if (termID >= sourceTermVars.size() || sourceTermVars[termID] < 2 ||
        sourceVariableExpressions.find(sourceTermVars[termID]) !=
            sourceVariableExpressions.end() ||
        !seenSourceVariables.insert(sourceTermVars[termID]).second) {
      continue;
    }
    const auto& term = dnl->getDNLTerminalFromID(termID);
    pendingInputs.push_back({termID, sourceTermVars[termID], terminalKey(term),
                             terminalName(term)});
  }
  std::sort(pendingInputs.begin(), pendingInputs.end(),
            [](const auto& lhs, const auto& rhs) {
              if (SignalKeyLess{}(lhs.key, rhs.key)) {
                return true;
              }
              if (SignalKeyLess{}(rhs.key, lhs.key)) {
                return false;
              }
              return lhs.displayName < rhs.displayName;
            });
  for (auto& input : pendingInputs) {
    const size_t variable =
        system.addInput(std::move(input.key), std::move(input.displayName));
    sourceVariableExpressions.emplace(input.sourceVariable,
                                      system.variableExpression(variable));
  }

  // Preserve the sorted root order while copying the temporary BoolExpr DAG.
  // Expression IDs are part of the returned representation, so deriving them
  // from unordered-map iteration would make otherwise identical collections
  // library-implementation dependent.
  std::vector<std::pair<DNLID, BoolExpr*>> materializedBoolExprs;
  materializedBoolExprs.reserve(materializationOutputs.size());
  const auto& materializedPOs = materializer.getPOs();
  for (size_t index = 0; index < materializationOutputs.size(); ++index) {
    BoolExpr* expression =
        index < materializedPOs.size() ? materializedPOs[index] : nullptr;
    if (expression != nullptr && expression->isValid()) {
      materializedBoolExprs.emplace_back(materializationOutputs[index],
                                         expression);
      continue;
    }
    const DNLID termID = materializationOutputs[index];
    const auto skipped = materializer.getSkippedOutputs().find(termID);
    if (skipped == materializer.getSkippedOutputs().end()) {
      system.addDiagnostic(
          {terminalName(dnl->getDNLTerminalFromID(termID)),
           "failed to materialize combinational root: builder returned no "
           "valid expression",
           true});
    }
  }
  for (const auto& [termID, skipped] : materializer.getSkippedOutputs()) {
    system.addDiagnostic(
        {terminalName(dnl->getDNLTerminalFromID(termID)),
         "failed to materialize combinational root: " + skipped.detail, true});
  }

  std::unordered_map<BoolExpr*, ExprID> boolMemo;
  std::unordered_map<DNLID, ExprID> termExpressions;
  for (const auto& [termID, boolExpression] : materializedBoolExprs) {
    std::string error;
    const ExprID expression =
        compileBoolExpr(boolExpression, sourceVariableExpressions,
                        system.expressions(), boolMemo, &error);
    if (expression == InvalidExpr) {
      system.addDiagnostic({terminalName(dnl->getDNLTerminalFromID(termID)),
                            "failed to own combinational expression: " + error,
                            true});
      continue;
    }
    termExpressions.emplace(termID, expression);
  }

  auto compileModelExpression = [&](const BooleanExpression& expression,
                                    const PendingInstance& pending,
                                    std::string* error) -> ExprID {
    using Operator = BooleanExpression::Operator;
    if (!expression.isValid()) {
      if (error != nullptr) {
        *error = "invalid Naja Boolean expression";
      }
      return InvalidExpr;
    }
    std::vector<ExprID> memo(expression.nodes.size(), InvalidExpr);
    std::vector<uint8_t> active(expression.nodes.size(), 0);
    std::function<ExprID(size_t)> compileNode = [&](size_t nodeID) -> ExprID {
      if (nodeID >= expression.nodes.size()) {
        if (error != nullptr) {
          *error = "Naja expression operand is out of range";
        }
        return InvalidExpr;
      }
      if (memo[nodeID] != InvalidExpr) {
        return memo[nodeID];
      }
      if (active[nodeID] != 0) {
        if (error != nullptr) {
          *error = "Naja expression contains a cycle";
        }
        return InvalidExpr;
      }
      active[nodeID] = 1;
      const auto& node = expression.nodes[nodeID];
      ExprID result = InvalidExpr;
      switch (node.operation) {
        case Operator::Constant:
          result = system.expressions().constant(node.constant);
          break;
        case Operator::Term: {
          const auto physical = pending.termIDs.find(node.term);
          if (physical == pending.termIDs.end()) {
            if (error != nullptr) {
              *error = "Naja expression references an unmapped terminal";
            }
            break;
          }
          const auto materialized = termExpressions.find(physical->second);
          if (materialized == termExpressions.end()) {
            if (error != nullptr) {
              *error = "missing materialized terminal expression";
            }
            break;
          }
          result = materialized->second;
          break;
        }
        case Operator::State:
          if (node.state >= pending.systemStateByModelState.size() ||
              pending.systemStateByModelState[node.state] ==
                  std::numeric_limits<size_t>::max()) {
            if (error != nullptr) {
              *error = "Naja expression references an unmapped state";
            }
            break;
          }
          result = system.variableExpression(
              system.states()[pending.systemStateByModelState[node.state]]
                  .variable);
          break;
        case Operator::Not:
          if (node.operands.size() != 1) {
            if (error != nullptr) {
              *error = "Naja NOT expression has invalid arity";
            }
            break;
          }
          result = system.expressions().makeNot(compileNode(node.operands[0]));
          break;
        case Operator::And:
        case Operator::Or:
        case Operator::Xor: {
          if (node.operands.empty()) {
            if (error != nullptr) {
              *error = "Naja Boolean operation has no operands";
            }
            break;
          }
          result = compileNode(node.operands.front());
          for (size_t operand = 1;
               operand < node.operands.size() && result != InvalidExpr;
               ++operand) {
            const ExprID rhs = compileNode(node.operands[operand]);
            if (node.operation == Operator::And) {
              result = system.expressions().makeAnd(result, rhs);
            } else if (node.operation == Operator::Or) {
              result = system.expressions().makeOr(result, rhs);
            } else {
              result = system.expressions().makeXor(result, rhs);
            }
          }
          break;
        }
      }
      active[nodeID] = 0;
      memo[nodeID] = result;
      return result;
    };
    return compileNode(expression.root);
  };

  for (auto& pending : instances) {
    for (size_t stateIndex = 0;
         stateIndex < pending.systemStateByModelState.size(); ++stateIndex) {
      const size_t systemState = pending.systemStateByModelState[stateIndex];
      if (systemState == std::numeric_limits<size_t>::max()) {
        continue;
      }
      const auto& stateModel = pending.model->states[stateIndex];
      std::string error;
      ExprID data =
          compileModelExpression(stateModel.nextState, pending, &error);
      if (data == InvalidExpr) {
        system.addDiagnostic({pending.path,
                              "failed to normalize state " +
                                  std::to_string(stateIndex) + ": " + error,
                              true});
        continue;
      }
      const ExprID current =
          system.variableExpression(system.states()[systemState].variable);
      ExprID next = data;
      ExprID transparency = InvalidExpr;
      if (pending.model->kind == SequentialModel::Kind::Latch) {
        error.clear();
        transparency =
            compileModelExpression(pending.model->clockedOn, pending, &error);
        if (transparency == InvalidExpr) {
          system.addDiagnostic(
              {pending.path, "failed to normalize latch transparency: " + error,
               true});
          continue;
        }
        next = system.expressions().makeIte(transparency, data, current);
      }

      std::optional<ExprID> clear;
      std::optional<ExprID> preset;
      if (stateModel.clear.has_value()) {
        error.clear();
        clear = compileModelExpression(*stateModel.clear, pending, &error);
        if (*clear == InvalidExpr) {
          system.addDiagnostic(
              {pending.path, "failed to normalize clear: " + error, true});
          continue;
        }
      }
      if (stateModel.preset.has_value()) {
        error.clear();
        preset = compileModelExpression(*stateModel.preset, pending, &error);
        if (*preset == InvalidExpr) {
          system.addDiagnostic(
              {pending.path, "failed to normalize preset: " + error, true});
          continue;
        }
      }

      if (clear.has_value() && preset.has_value()) {
        ExprID simultaneous = InvalidExpr;
        switch (convertClearPresetValue(stateModel.clearPresetValue)) {
          case ClearPresetValue::Zero:
            simultaneous = system.expressions().constant(false);
            break;
          case ClearPresetValue::One:
            simultaneous = system.expressions().constant(true);
            break;
          case ClearPresetValue::Hold:
            simultaneous = current;
            break;
          case ClearPresetValue::Toggle:
            simultaneous = system.expressions().makeNot(current);
            break;
          case ClearPresetValue::Unknown:
            system.addDiagnostic(
                {pending.path, "simultaneous clear/preset behavior is unknown",
                 true});
            break;
        }
        if (simultaneous == InvalidExpr) {
          continue;
        }
        next = system.expressions().makeIte(
            *clear,
            system.expressions().makeIte(*preset, simultaneous,
                                         system.expressions().constant(false)),
            system.expressions().makeIte(
                *preset, system.expressions().constant(true), next));
      } else if (clear.has_value()) {
        next = system.expressions().makeIte(
            *clear, system.expressions().constant(false), next);
      } else if (preset.has_value()) {
        next = system.expressions().makeIte(
            *preset, system.expressions().constant(true), next);
      }

      system.setNextState(systemState, next);
      if (pending.model->kind == SequentialModel::Kind::Latch) {
        NormalizedLatch latch{
            systemState,
            pending.path,
            pending.modelName,
            data,
            transparency,
            clear,
            preset,
            convertClearPresetValue(stateModel.clearPresetValue),
            {}};
        for (const auto& [termID, outputState] : pending.modeledOutputs) {
          if (outputState.first != stateIndex) {
            continue;
          }
          const auto& term = dnl->getDNLTerminalFromID(termID);
          latch.outputs.push_back(
              {terminalKey(term), terminalName(term), outputState.second});
        }
        std::sort(latch.outputs.begin(), latch.outputs.end(),
                  [](const auto& lhs, const auto& rhs) {
                    if (SignalKeyLess{}(lhs.key, rhs.key)) {
                      return true;
                    }
                    if (SignalKeyLess{}(rhs.key, lhs.key)) {
                      return false;
                    }
                    return lhs.displayName < rhs.displayName;
                  });
        system.addLatch(std::move(latch));
      }
    }
  }

  for (const DNLID termID : topOutputTerms) {
    const auto expression = termExpressions.find(termID);
    if (expression == termExpressions.end()) {
      continue;
    }
    const auto& term = dnl->getDNLTerminalFromID(termID);
    system.addOutput(
        {terminalKey(term), terminalName(term), expression->second});
  }

  std::string validationReason;
  if (!system.hasFatalDiagnostics() && !system.validate(&validationReason)) {
    system.addDiagnostic(
        {top->getName().getString(),
         "normalized transition system is invalid: " + validationReason, true});
  }
  return system;
}

}  // namespace KEPLER_FORMAL::SEC::PHASE
