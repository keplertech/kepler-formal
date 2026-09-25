// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "SNLDesignModeling.h"
#include "latch/LatchSymbolicModel.h"

namespace KEPLER_FORMAL::SEC::LATCH::detail {

// The concrete and symbolic callbacks share the same copied, validated Naja
// bytecode and pin mapping. No source-netlist pointers are accessed here.
inline BoolExpr* symbolicFormula(
    const naja::NL::SNLDesignModeling::BooleanExpression& expression,
    const std::vector<size_t>& pins, const SymbolicBits& input,
    const SymbolicBits& storage) {
  using Operator = naja::NL::SNLDesignModeling::BooleanExpression::Operator;
  SymbolicBits values(expression.nodes.size());
  for (size_t i = 0; i < values.size(); ++i) {
    const auto& node = expression.nodes[i];
    switch (node.operation) {
      case Operator::Constant: values[i] = BoolExpr::Var(node.constant ? 1 : 0); break;
      case Operator::Term: values[i] = input.at(pins.at(i)); break;
      case Operator::State: values[i] = storage.at(node.state); break;
      case Operator::Not: values[i] = BoolExpr::Not(values.at(node.operands.at(0))); break;
      case Operator::And: case Operator::Or: case Operator::Xor: {
        auto* value = BoolExpr::Var(node.operation == Operator::And ? 1 : 0);
        for (auto operand : node.operands) {
          if (node.operation == Operator::And) value = BoolExpr::And(value, values.at(operand));
          else if (node.operation == Operator::Or) value = BoolExpr::Or(value, values.at(operand));
          else value = BoolExpr::Xor(value, values.at(operand));
        }
        values[i] = value;
        break;
      }
      default: throw std::runtime_error("unknown sequential expression operator");
    }
  }
  return values.at(expression.root);
}

template<class Rule>
SymbolicBits symbolicSequentialOutputs(const Rule& rule, const SymbolicBits& state,
                                       const SymbolicBits& pins) {
  SymbolicBits result;
  for (const auto& expression : rule.outputs) result.push_back(expression(pins, state));
  return result;
}

template<class Rule>
SymbolicReaction symbolicSequentialReaction(const Rule& cell, const SymbolicBits& old,
    const SymbolicBits& before, const SymbolicBits& now,
    std::optional<size_t> changed, bool bootstrap) {
  using Value = naja::NL::SNLDesignModeling::SequentialState::ClearPresetValue;
  SymbolicReaction result;
  auto* zero = BoolExpr::createFalse();
  auto* one = BoolExpr::createTrue();
  auto* open = cell.control(now, old);
  auto* capture = cell.latch ? open : (!bootstrap && changed.has_value()
      ? BoolExpr::And(BoolExpr::Not(cell.control(before, old)), open) : zero);
  for (size_t i = 0; i < cell.states.size(); ++i) {
    const auto& rule = cell.states[i];
    auto* clear = rule.clear ? (*rule.clear)(now, old) : zero;
    auto* preset = rule.preset ? (*rule.preset)(now, old) : zero;
    auto* both = BoolExpr::And(clear, preset);
    auto* conflict = old.at(i);
    switch (rule.conflict) {
      case Value::Zero: conflict = zero; break;
      case Value::One: conflict = one; break;
      case Value::Hold: break;
      case Value::Toggle: case Value::Unknown:
        result.error = BoolExpr::Or(result.error, both);
        break;
      default: throw std::runtime_error("unknown asynchronous conflict behavior");
    }
    auto* normal = symbolicMux(capture, rule.data(now, old), old.at(i));
    result.storage.push_back(symbolicMux(both, conflict,
        symbolicMux(clear, zero, symbolicMux(preset, one, normal))));
  }
  result.outputs = symbolicSequentialOutputs(cell, result.storage, now);
  return result;
}

inline BoolExpr* symbolicTruthTable(const naja::NL::SNLTruthTable& truth,
    const std::vector<size_t>& pins, const SymbolicBits& input) {
  using Type = naja::NL::SNLTruthTable::GenericType;
  if (!truth.isGeneric()) {
    // This is expansion of the provided local cell table, never a circuit's
    // input/state space. Refuse pathological tables before allocating a DAG.
    if (pins.size() > 16) throw Limit("symbolic local truth table exceeds 16 inputs");
    SymbolicBits layer;
    const size_t width = size_t(1) << pins.size();
    layer.reserve(width);
    for (size_t i = 0; i < width; ++i) layer.push_back(BoolExpr::Var(truth.bits().bit(i) ? 1 : 0));
    for (size_t pin = 0; pin < pins.size(); ++pin) {
      for (size_t i = 0; i < layer.size() / 2; ++i)
        layer[i] = symbolicMux(input.at(pins[pin]), layer[2*i + 1], layer[2*i]);
      layer.resize(layer.size() / 2);
    }
    return layer.front();
  }
  const auto type = truth.getGenericType();
  auto* value = BoolExpr::Var(type == Type::AND || type == Type::NAND ? 1 : 0);
  for (size_t pin : pins) {
    if (type == Type::AND || type == Type::NAND) value = BoolExpr::And(value, input.at(pin));
    else if (type == Type::OR || type == Type::NOR) value = BoolExpr::Or(value, input.at(pin));
    else if (type == Type::XOR || type == Type::XNOR) value = BoolExpr::Xor(value, input.at(pin));
    else throw std::runtime_error("unsupported symbolic generic truth table");
  }
  return type == Type::NAND || type == Type::NOR || type == Type::XNOR ? BoolExpr::Not(value) : value;
}

}  // namespace KEPLER_FORMAL::SEC::LATCH::detail
