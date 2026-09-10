// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace KEPLER_FORMAL {

class BoolExpr;

// Streaming BTOR2 writer for an already bit-blasted Boolean transition system.
// Bind all BoolExpr variables before exporting expressions that reference them.
// Expressions must remain alive (and their cache must not be cleared) while the
// writer is used. Every emitted value has the one-bit bit-vector sort.
class Btor2Writer {
 public:
  using NodeId = uint64_t;

  explicit Btor2Writer(std::ostream& out);
  Btor2Writer(const Btor2Writer&) = delete;
  Btor2Writer& operator=(const Btor2Writer&) = delete;

  NodeId input(const std::string& symbol = {});
  NodeId state(const std::string& symbol = {});
  NodeId constant(bool value);
  void bindVariable(size_t symbol, NodeId node);
  NodeId expression(BoolExpr* expr);

  NodeId logicalNot(NodeId value);
  NodeId logicalAnd(NodeId left, NodeId right);
  NodeId logicalOr(NodeId left, NodeId right);
  NodeId logicalXor(NodeId left, NodeId right);
  NodeId equal(NodeId left, NodeId right);

  // BTOR2 requires the initialization value to precede the state declaration
  // and forbids inputs in its expression cone. Declare constants/states used
  // for initialization before declaring the state they initialize.
  NodeId init(NodeId state, NodeId value);
  NodeId next(NodeId state, NodeId value);
  NodeId constraint(NodeId value);
  NodeId bad(NodeId value, const std::string& symbol = {});
  NodeId output(NodeId value, const std::string& symbol = {});
  void comment(const std::string& text);

 private:
  enum class NodeKind { Other, Value, State, Input };
  static constexpr NodeId boolSort_ = 1;

  static std::string symbolSuffix(const std::string& symbol);
  NodeId emit(const std::string& body, NodeKind kind);
  void requireValue(NodeId node) const;
  void requireState(NodeId node) const;
  NodeId binary(const char* op, NodeId left, NodeId right);
  NodeId property(const char* op, NodeId value, const std::string& symbol);
  NodeId stateRelation(const char* op, NodeId state, NodeId value,
                       std::unordered_set<NodeId>& assigned);

  std::ostream& out_;
  // Index zero is unused; sort and directive IDs are not expression values.
  std::vector<NodeKind> nodes_{NodeKind::Other};
  std::vector<bool> dependsOnInput_{false};
  NodeId constants_[2] = {0, 0};
  std::unordered_map<size_t, NodeId> variables_;
  std::unordered_map<const BoolExpr*, NodeId> expressions_;
  std::unordered_set<NodeId> initialized_;
  std::unordered_set<NodeId> transitioned_;
};

}  // namespace KEPLER_FORMAL
