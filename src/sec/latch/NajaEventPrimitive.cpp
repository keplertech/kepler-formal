// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/NajaEventPrimitive.h"

#include <algorithm>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include "NLBitDependencies.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using Modeling = naja::NL::SNLDesignModeling;
using Expression = Modeling::BooleanExpression;
using Operator = Expression::Operator;
using Term = naja::NL::SNLBitTerm;

// Copy the expression and resolve pin identities once. No Naja objects are read
// by worker callbacks, so compact extraction can release the source netlist.
struct Formula {
  Expression expression;
  std::vector<size_t> pins;
  bool operator()(const Bits& input, const Bits& storage) const {
    Bits values(expression.nodes.size());
    for (size_t i = 0; i < values.size(); ++i) {
      const auto& node = expression.nodes[i];
      switch (node.operation) {
        case Operator::Constant: values[i] = node.constant; break;
        case Operator::Term: values[i] = input.at(pins[i]); break;
        case Operator::State: values[i] = storage.at(node.state); break;
        case Operator::Not: values[i] = !values.at(node.operands[0]); break;
        case Operator::And: case Operator::Or: case Operator::Xor: {
          bool value = node.operation == Operator::And;
          for (auto operand : node.operands) {
            if (node.operation == Operator::And) value &= values[operand];
            else if (node.operation == Operator::Or) value |= values[operand];
            else value ^= values[operand];
          }
          values[i] = value;
          break;
        }
      }
    }
    return values.at(expression.root);
  }
};

Formula formula(const Expression& source, const std::map<const Term*, size_t>& pins,
                size_t storageBits, bool allowState = true) {
  if (!source.isValid()) throw std::runtime_error("missing sequential expression");
  Formula result{source, std::vector<size_t>(source.nodes.size())};
  for (size_t i = 0; i < source.nodes.size(); ++i) {
    const auto& node = source.nodes[i];
    if (node.operation == Operator::Term) {
      const auto pin = pins.find(node.term);
      if (pin == pins.end()) throw std::runtime_error("expression references a non-input pin");
      result.pins[i] = pin->second;
    } else if (node.operation == Operator::State &&
               (!allowState || node.state >= storageBits)) {
      throw std::runtime_error("invalid state reference in sequential expression");
    }
    if ((node.operation == Operator::Not && node.operands.size() != 1) ||
        ((node.operation == Operator::And || node.operation == Operator::Or ||
          node.operation == Operator::Xor) && node.operands.empty()))
      throw std::runtime_error("invalid sequential expression arity");
    for (auto child : node.operands)
      if (child >= i) throw std::runtime_error("cyclic or unordered sequential expression");
  }
  // The copied expression no longer needs pointers into the source netlist.
  for (auto& node : result.expression.nodes) node.term = nullptr;
  return result;
}

struct StoredRule {
  Formula data;
  std::optional<Formula> clear, preset;
  Modeling::SequentialState::ClearPresetValue conflict;
};

struct SequentialRule {
  bool latch = false;
  Formula control;
  std::vector<StoredRule> states;
  std::vector<Formula> outputs;

  Bits output(const Bits& state, const Bits& pins) const {
    Bits result;
    for (const auto& expression : outputs) result.push_back(expression(pins, state));
    return result;
  }
  Reaction react(const Bits& old, const Bits& before, const Bits& now,
                 std::optional<size_t> changed, bool bootstrap) const {
    Reaction result;
    result.storage = old;
    const bool open = control(now, old);
    const bool capture = latch ? open :
        (!bootstrap && changed.has_value() && !control(before, old) && open);
    for (size_t i = 0; i < states.size(); ++i) {
      const auto& rule = states[i];
      const bool clear = rule.clear && (*rule.clear)(now, old);
      const bool preset = rule.preset && (*rule.preset)(now, old);
      if (clear && preset) {
        using Value = Modeling::SequentialState::ClearPresetValue;
        switch (rule.conflict) {
          case Value::Zero: result.storage[i] = 0; break;
          case Value::One: result.storage[i] = 1; break;
          case Value::Hold: break;
          case Value::Toggle:
            result.error = true;
            result.reason = "state-dependent simultaneous clear/preset toggle is unsupported";
            break;
          case Value::Unknown:
            result.error = true;
            result.reason = "undefined simultaneous clear/preset";
            break;
        }
      } else if (clear) result.storage[i] = 0;
      else if (preset) result.storage[i] = 1;
      else if (capture) result.storage[i] = rule.data(now, old);
    }
    result.outputs = output(result.storage, now);
    return result;
  }
};

struct Table {
  naja::NL::SNLTruthTable truth;
  std::vector<size_t> pins;
  bool operator()(const Bits& input) const {
    using Type = naja::NL::SNLTruthTable::GenericType;
    if (!truth.isGeneric()) {
      size_t index = 0;
      for (size_t i = 0; i < pins.size(); ++i) index |= size_t(input[pins[i]]) << i;
      return truth.bits().bit(index);
    }
    const auto type = truth.getGenericType();
    bool value = type == Type::AND || type == Type::NAND;
    for (size_t pin : pins) {
      if (type == Type::AND || type == Type::NAND) value &= input[pin];
      else if (type == Type::OR || type == Type::NOR) value |= input[pin];
      else value ^= input[pin];
    }
    return (type == Type::NAND || type == Type::NOR || type == Type::XNOR) ? !value : value;
  }
};
}  // namespace

Primitive makeNajaEventPrimitive(naja::NL::SNLInstance* instance, std::string path,
                                const std::map<const Term*, size_t>& nets) {
  if (!instance) throw std::runtime_error("missing leaf instance");
  Primitive primitive;
  primitive.name = std::move(path);
  std::map<const Term*, size_t> pins;
  std::map<size_t, size_t> pinByOrder;
  std::vector<Term*> outputs;
  for (auto* term : instance->getModel()->getBitTerms()) {
    if (term->getDirection() == Term::Direction::Input) {
      pins.emplace(term, primitive.inputs.size());
      pinByOrder.emplace(term->getOrderID(), primitive.inputs.size());
      primitive.inputs.push_back(nets.at(term));
    } else if (term->getDirection() == Term::Direction::Output) {
      primitive.outputs.push_back(nets.at(term));
      outputs.push_back(term);
    } else throw std::runtime_error("bidirectional or undefined primitive pin");
  }
  if (outputs.empty()) throw std::runtime_error("outputless primitive has no supported Boolean event model");
  if (Modeling::hasSequentialModel(instance->getModel())) {
    const auto& model = Modeling::getSequentialModel(instance->getModel());
    if (!model.isValid()) throw std::runtime_error("invalid explicit sequential model");
    if (model.kind != Modeling::SequentialModel::Kind::Latch &&
        model.kind != Modeling::SequentialModel::Kind::FlipFlop)
      throw std::runtime_error("unknown sequential element kind");
    std::set<const Term*> declaredOutputs;
    for (const auto& output : model.outputs)
      if (!output.term || output.term->getDesign() != instance->getModel() ||
          output.term->getDirection() != Term::Direction::Output ||
          !declaredOutputs.insert(output.term).second)
        throw std::runtime_error("invalid or duplicate sequential output association");
    primitive.storageBits = model.states.size();
    auto rule = std::make_shared<SequentialRule>();
    rule->latch = model.kind == Modeling::SequentialModel::Kind::Latch;
    rule->control = formula(model.clockedOn, pins, model.states.size(), false);
    for (const auto& state : model.states) {
      // Internal state references are not pin events. Latch data feedback and
      // state-dependent async controls need additional activation semantics;
      // never certify them from the pin-only graph as if they were quiescent.
      StoredRule item{formula(state.nextState, pins, model.states.size(), !rule->latch), {}, {}, state.clearPresetValue};
      if (state.clear) item.clear = formula(*state.clear, pins, model.states.size(), false);
      if (state.preset) item.preset = formula(*state.preset, pins, model.states.size(), false);
      rule->states.push_back(std::move(item));
    }
    for (auto* term : outputs) {
      const auto it = std::find_if(model.outputs.begin(), model.outputs.end(),
          [term](const auto& output) { return output.term == term; });
      if (it == model.outputs.end()) throw std::runtime_error("missing physical sequential output");
      rule->outputs.push_back(formula(it->function, pins, model.states.size()));
    }
    primitive.react = [rule](const Bits& state, const Bits& before, const Bits& now,
                             std::optional<size_t> changed, bool bootstrap) {
      return rule->react(state, before, now, changed, bootstrap);
    };
    primitive.initialOutputValues = [rule](const Bits& state, const Bits& pins) {
      return rule->output(state, pins);
    };
  } else {
    std::vector<Table> tables;
    for (auto* output : outputs) {
      Table table{Modeling::getTruthTable(instance, output->getOrderID()), {}};
      if (!table.truth.isInitialized()) throw std::runtime_error("no explicit sequential model or truth table");
      using Type = naja::NL::SNLTruthTable::GenericType;
      const auto type = table.truth.getGenericType();
      if (type == Type::TABLE_SELECT || type == Type::DIVMOD)
        throw std::runtime_error("generic arithmetic/table-select primitive not supported by Boolean event adapter");
      for (auto dependency : naja::NL::NLBitDependencies::decodeBits(table.truth.getDependencies())) {
        const auto pin = pinByOrder.find(dependency);
        if (pin == pinByOrder.end()) throw std::runtime_error("truth table references a non-input pin");
        table.pins.push_back(pin->second);
      }
      if (table.pins.size() != table.truth.size()) throw std::runtime_error("truth table dependency arity mismatch");
      if (!table.truth.isGeneric() && table.pins.size() >= sizeof(size_t) * 8)
        throw std::runtime_error("truth table exceeds event index width");
      tables.push_back(std::move(table));
    }
    primitive.react = [tables = std::move(tables)](const Bits&, const Bits&, const Bits& input,
                                                 std::optional<size_t>, bool) {
      Reaction result;
      for (const auto& table : tables) result.outputs.push_back(table(input));
      return result;
    };
  }
  return primitive;
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
