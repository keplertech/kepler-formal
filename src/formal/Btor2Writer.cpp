// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "Btor2Writer.h"

#include <limits>
#include <ostream>
#include <stdexcept>

#include "BoolExpr.h"

namespace KEPLER_FORMAL {

Btor2Writer::Btor2Writer(std::ostream& out) : out_(out) {
  emit("sort bitvec 1", NodeKind::Other);
}

std::string Btor2Writer::symbolSuffix(const std::string& symbol) {
  if (symbol.empty()) {
    return {};
  }
  // Keep readable net names, but encode whitespace, comment delimiters,
  // non-ASCII bytes and '%' itself to avoid collisions or injected lines.
  static constexpr char hex[] = "0123456789ABCDEF";
  std::string result = " ";
  for (unsigned char c : symbol) {
    if (c > 32 && c < 127 && c != ';' && c != '%') {
      result += static_cast<char>(c);
    } else {
      result += '%';
      result += hex[c >> 4];
      result += hex[c & 15];
    }
  }
  return result;
}

Btor2Writer::NodeId Btor2Writer::emit(const std::string& body, NodeKind kind) {
  if (nodes_.size() >= static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
    throw std::overflow_error("BTOR2 node ID limit exceeded");
  }
  const NodeId id = nodes_.size();
  const std::string line = std::to_string(id) + " " + body + '\n';
  out_.write(line.data(), static_cast<std::streamsize>(line.size()));
  if (!out_) {
    throw std::runtime_error("Failed to write BTOR2 output");
  }
  nodes_.push_back(kind);
  dependsOnInput_.push_back(kind == NodeKind::Input);
  return id;
}

void Btor2Writer::requireValue(NodeId node) const {
  if (node >= nodes_.size() || nodes_[node] == NodeKind::Other) {
    throw std::invalid_argument("Invalid BTOR2 value node " + std::to_string(node));
  }
}

void Btor2Writer::requireState(NodeId node) const {
  if (node >= nodes_.size() || nodes_[node] != NodeKind::State) {
    throw std::invalid_argument("Invalid BTOR2 state node " + std::to_string(node));
  }
}

Btor2Writer::NodeId Btor2Writer::input(const std::string& symbol) {
  return emit("input " + std::to_string(boolSort_) + symbolSuffix(symbol),
              NodeKind::Input);
}

Btor2Writer::NodeId Btor2Writer::state(const std::string& symbol) {
  return emit("state " + std::to_string(boolSort_) + symbolSuffix(symbol),
              NodeKind::State);
}

Btor2Writer::NodeId Btor2Writer::constant(bool value) {
  NodeId& result = constants_[value ? 1 : 0];
  if (result == 0) {
    result = emit("const " + std::to_string(boolSort_) + (value ? " 1" : " 0"),
                  NodeKind::Value);
  }
  return result;
}

void Btor2Writer::bindVariable(size_t symbol, NodeId node) {
  requireValue(node);
  if (symbol < 2 || symbol == std::numeric_limits<size_t>::max()) {
    throw std::invalid_argument("Cannot bind reserved BoolExpr variable " +
                                std::to_string(symbol));
  }
  const auto result = variables_.emplace(symbol, node);
  if (!result.second && result.first->second != node) {
    throw std::invalid_argument("BoolExpr variable already bound: " +
                                std::to_string(symbol));
  }
}

Btor2Writer::NodeId Btor2Writer::logicalNot(NodeId value) {
  requireValue(value);
  const NodeId result = emit(
      "not " + std::to_string(boolSort_) + " " + std::to_string(value),
      NodeKind::Value);
  dependsOnInput_[result] = dependsOnInput_[value];
  return result;
}

Btor2Writer::NodeId Btor2Writer::binary(const char* op, NodeId left, NodeId right) {
  requireValue(left);
  requireValue(right);
  const NodeId result = emit(
      std::string(op) + " " + std::to_string(boolSort_) + " " +
          std::to_string(left) + " " + std::to_string(right),
      NodeKind::Value);
  dependsOnInput_[result] = dependsOnInput_[left] || dependsOnInput_[right];
  return result;
}

Btor2Writer::NodeId Btor2Writer::logicalAnd(NodeId left, NodeId right) {
  return binary("and", left, right);
}

Btor2Writer::NodeId Btor2Writer::logicalOr(NodeId left, NodeId right) {
  return binary("or", left, right);
}

Btor2Writer::NodeId Btor2Writer::logicalXor(NodeId left, NodeId right) {
  return binary("xor", left, right);
}

Btor2Writer::NodeId Btor2Writer::equal(NodeId left, NodeId right) {
  return binary("eq", left, right);
}

Btor2Writer::NodeId Btor2Writer::expression(BoolExpr* expr) {
  struct Frame {
    const BoolExpr* node;
    bool visited;
  };
  std::vector<Frame> stack{{expr, false}};
  std::unordered_set<const BoolExpr*> active;
  while (!stack.empty()) {
    const Frame frame = stack.back();
    stack.pop_back();
    const BoolExpr* node = frame.node;
    if (node == nullptr || !node->isValid()) {
      throw std::invalid_argument("Cannot export invalid BoolExpr to BTOR2");
    }
    if (expressions_.count(node) != 0) {
      continue;
    }
    const Op op = node->getOp();
    if (op == Op::VAR) {
      const size_t symbol = node->getId();
      if (symbol < 2) {
        expressions_.emplace(node, constant(symbol == 1));
      } else {
        const auto binding = variables_.find(symbol);
        if (binding == variables_.end()) {
          throw std::invalid_argument("Unbound BoolExpr variable in BTOR2 export: " +
                                      std::to_string(symbol));
        }
        expressions_.emplace(node, binding->second);
      }
      continue;
    }
    if (op != Op::NOT && op != Op::AND && op != Op::OR && op != Op::XOR) {
      throw std::invalid_argument("Unsupported BoolExpr operator in BTOR2 export");
    }
    if (node->getLeft() == nullptr ||
        (op != Op::NOT && node->getRight() == nullptr) ||
        (op == Op::NOT && node->getRight() != nullptr)) {
      throw std::invalid_argument("Malformed BoolExpr operands in BTOR2 export");
    }
    if (!frame.visited) {
      if (!active.insert(node).second) {
        throw std::invalid_argument("Cyclic BoolExpr in BTOR2 export");
      }
      stack.push_back({node, true});
      if (op != Op::NOT) {
        stack.push_back({node->getRight(), false});
      }
      stack.push_back({node->getLeft(), false});
      continue;
    }
    const NodeId left = expressions_.at(node->getLeft());
    NodeId result = 0;
    if (op == Op::NOT) {
      result = logicalNot(left);
    } else {
      const NodeId right = expressions_.at(node->getRight());
      switch (op) {
        case Op::AND: result = logicalAnd(left, right); break;
        case Op::OR: result = logicalOr(left, right); break;
        case Op::XOR: result = logicalXor(left, right); break;
        default: throw std::logic_error("Unexpected BTOR2 binary operator");
      }
    }
    expressions_.emplace(node, result);
    active.erase(node);
  }
  return expressions_.at(expr);
}

Btor2Writer::NodeId Btor2Writer::stateRelation(
    const char* op, NodeId state, NodeId value,
    std::unordered_set<NodeId>& assigned) {
  requireState(state);
  requireValue(value);
  if (assigned.count(state) != 0) {
    throw std::invalid_argument(std::string("Duplicate BTOR2 ") + op +
                                " for state " + std::to_string(state));
  }
  const NodeId result = emit(std::string(op) + " " + std::to_string(boolSort_) +
                                 " " + std::to_string(state) + " " +
                                 std::to_string(value),
                             NodeKind::Other);
  assigned.insert(state);
  return result;
}

Btor2Writer::NodeId Btor2Writer::init(NodeId state, NodeId value) {
  requireState(state);
  requireValue(value);
  if (value >= state) {
    throw std::invalid_argument("BTOR2 initialization value must precede its state");
  }
  if (dependsOnInput_[value]) {
    throw std::invalid_argument("BTOR2 initialization cannot depend on inputs");
  }
  return stateRelation("init", state, value, initialized_);
}

Btor2Writer::NodeId Btor2Writer::next(NodeId state, NodeId value) {
  return stateRelation("next", state, value, transitioned_);
}

Btor2Writer::NodeId Btor2Writer::property(
    const char* op, NodeId value, const std::string& symbol) {
  requireValue(value);
  return emit(std::string(op) + " " + std::to_string(value) + symbolSuffix(symbol),
              NodeKind::Other);
}

Btor2Writer::NodeId Btor2Writer::constraint(NodeId value) {
  return property("constraint", value, {});
}

Btor2Writer::NodeId Btor2Writer::bad(NodeId value, const std::string& symbol) {
  return property("bad", value, symbol);
}

Btor2Writer::NodeId Btor2Writer::output(NodeId value, const std::string& symbol) {
  return property("output", value, symbol);
}

void Btor2Writer::comment(const std::string& text) {
  std::string line = "; ";
  for (unsigned char c : text) {
    line += (c >= 32 && c != 127) ? static_cast<char>(c) : ' ';
  }
  line += '\n';
  out_.write(line.data(), static_cast<std::streamsize>(line.size()));
  if (!out_) {
    throw std::runtime_error("Failed to write BTOR2 comment");
  }
}

}  // namespace KEPLER_FORMAL
