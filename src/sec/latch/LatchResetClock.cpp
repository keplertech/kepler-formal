// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchResetClock.h"

#include <deque>
#include <stdexcept>
#include <unordered_map>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
struct Route {
  enum class Kind { Unknown, Constant, Literal };
  Kind kind = Kind::Unknown;
  size_t input = 0;
  // A constant value, or the inversion of the external input literal.
  bool bit = false;
};

Route constant(bool bit) { return {Route::Kind::Constant, 0, bit}; }
Route invert(Route value) {
  if (value.kind != Route::Kind::Unknown) value.bit = !value.bit;
  return value;
}
bool isConstant(const Route& value, bool bit) {
  return value.kind == Route::Kind::Constant && value.bit == bit;
}
Route combine(Op op, Route a, Route b) {
  if (op == Op::AND) {
    if (isConstant(a, false) || isConstant(b, false)) return constant(false);
    if (isConstant(a, true)) return b;
    if (isConstant(b, true)) return a;
  } else if (op == Op::OR) {
    if (isConstant(a, true) || isConstant(b, true)) return constant(true);
    if (isConstant(a, false)) return b;
    if (isConstant(b, false)) return a;
  } else if (op == Op::XOR) {
    if (a.kind == Route::Kind::Constant) return a.bit ? invert(b) : b;
    if (b.kind == Route::Kind::Constant) return b.bit ? invert(a) : a;
  } else throw std::invalid_argument("invalid clock routing operator");
  if (a.kind == Route::Kind::Literal && b.kind == Route::Kind::Literal &&
      a.input == b.input) {
    if (op == Op::XOR) return constant(a.bit != b.bit);
    if (a.bit == b.bit) return a;
    return constant(op == Op::OR);
  }
  return {};
}

Route expressionRoute(BoolExpr* root, const Primitive& primitive,
                      const std::vector<Route>& routes) {
  std::unordered_map<BoolExpr*, Route> memo;
  std::vector<std::pair<BoolExpr*, bool>> pending{{root, false}};
  while (!pending.empty()) {
    const auto [node, ready] = pending.back();
    pending.pop_back();
    if (!node || !node->isValid()) throw std::invalid_argument("missing or invalid clock routing expression");
    if (memo.count(node)) continue;
    if (node->getOp() == Op::VAR) {
      if (node->getId() < 2) memo.emplace(node, constant(node->getId() != 0));
      else {
        const size_t pin = node->getId() - 2;
        if (pin >= primitive.inputs.size())
          throw std::invalid_argument("clock routing expression references a non-input symbol");
        memo.emplace(node, routes.at(primitive.inputs[pin]));
      }
      continue;
    }
    const auto op = node->getOp();
    if (op != Op::NOT && op != Op::AND && op != Op::OR && op != Op::XOR)
      throw std::invalid_argument("invalid clock routing operator");
    if (!ready) {
      pending.emplace_back(node, true);
      if (op != Op::NOT) pending.emplace_back(node->getRight(), false);
      pending.emplace_back(node->getLeft(), false);
      continue;
    }
    const auto left = memo.at(node->getLeft());
    memo.emplace(node, op == Op::NOT ? invert(left) :
        combine(op, left, memo.at(node->getRight())));
  }
  return memo.at(root);
}

ResetClockDiscovery unsupported(std::string detail) {
  return {ResetClockDiscovery::Status::Unsupported, {}, {}, std::move(detail)};
}
}  // namespace

ResetClockDiscovery discoverResetClock(
    const Network& network, const std::vector<ResetClockPrimitive>& metadata) {
  try {
    if (metadata.size() != network.primitives.size())
      return unsupported("automatic reset clock discovery requires metadata for every primitive");
    if (!network.constantByNet.empty() && network.constantByNet.size() != network.netCount)
      return unsupported("invalid constant net table during automatic reset clock discovery");
    std::vector<Route> routes(network.netCount);
    std::vector<bool> hasSource(network.netCount, false);
    std::vector<std::vector<size_t>> consumers(network.netCount);
    for (size_t i = 0; i < network.externalInputs.size(); ++i) {
      const size_t net = network.externalInputs[i];
      if (net >= network.netCount || hasSource[net])
        return unsupported("invalid or duplicate external net during automatic reset clock discovery");
      hasSource[net] = true;
      routes[net] = {Route::Kind::Literal, i, false};
    }
    for (size_t net = 0; net < network.constantByNet.size(); ++net) {
      if (!network.constantByNet[net].has_value()) continue;
      if (hasSource[net]) return unsupported("external reset clock input is also a constant net");
      hasSource[net] = true;
      routes[net] = constant(*network.constantByNet[net]);
    }
    size_t edgeCells = 0;
    std::deque<size_t> work;
    std::vector<bool> queued(network.primitives.size(), false);
    for (size_t cell = 0; cell < network.primitives.size(); ++cell) {
      const auto& primitive = network.primitives[cell];
      const auto& meta = metadata[cell];
      if (meta.kind == ResetClockPrimitive::Kind::Unsupported)
        return unsupported("unclassified primitive prevents automatic reset clock discovery: " + primitive.name);
      if (meta.kind == ResetClockPrimitive::Kind::FlipFlop) ++edgeCells;
      if (meta.kind == ResetClockPrimitive::Kind::Combinational &&
          (primitive.storageBits != 0 || meta.outputs.size() != primitive.outputs.size()))
        return unsupported("incomplete combinational clock routing metadata: " + primitive.name);
      for (size_t net : primitive.inputs) {
        if (net >= network.netCount)
          return unsupported("out-of-range input net during automatic reset clock discovery");
        if (meta.kind == ResetClockPrimitive::Kind::Combinational) consumers[net].push_back(cell);
      }
      for (size_t net : primitive.outputs) {
        if (net >= network.netCount || hasSource[net])
          return unsupported("multiple drivers or invalid output net during automatic reset clock discovery");
        hasSource[net] = true;
      }
      if (meta.kind == ResetClockPrimitive::Kind::Combinational) {
        queued[cell] = true;
        work.push_back(cell);
      }
    }
    if (edgeCells == 0)
      return {ResetClockDiscovery::Status::NoEdgeClock, {}, {},
          "reset cycles require an explicit flip-flop clock; latch enables do not define a clock cycle"};

    // Route facts only become more precise (unknown -> constant/literal). Gates
    // revisit only when an input fact is discovered, so routing cycles terminate
    // without choosing an arbitrary clock or recursing through the netlist.
    while (!work.empty()) {
      const size_t cell = work.front();
      work.pop_front();
      queued[cell] = false;
      const auto& primitive = network.primitives[cell];
      for (size_t pin = 0; pin < primitive.outputs.size(); ++pin) {
        const size_t net = primitive.outputs[pin];
        if (routes[net].kind != Route::Kind::Unknown) continue;
        const auto route = expressionRoute(metadata[cell].outputs[pin], primitive, routes);
        if (route.kind == Route::Kind::Unknown) continue;
        routes[net] = route;
        for (size_t consumer : consumers[net]) if (!queued[consumer]) {
          queued[consumer] = true;
          work.push_back(consumer);
        }
      }
    }

    std::optional<size_t> carrier;
    for (size_t cell = 0; cell < network.primitives.size(); ++cell) {
      if (metadata[cell].kind != ResetClockPrimitive::Kind::FlipFlop) continue;
      const auto& primitive = network.primitives[cell];
      const auto clock = expressionRoute(metadata[cell].clock, primitive, routes);
      if (clock.kind == Route::Kind::Constant)
        return unsupported("flip-flop clock is constant, so automatic reset cycles cannot clock it: " + primitive.name);
      if (clock.kind != Route::Kind::Literal)
        return unsupported("flip-flop clock is gated, state-generated, cyclic, or has no unique external carrier: " + primitive.name);
      if (carrier && *carrier != clock.input)
        return unsupported("multiple independent flip-flop clock roots require an explicit reset clock protocol");
      carrier = clock.input;
    }
    return {ResetClockDiscovery::Status::Resolved, network.externalInputs.at(*carrier), carrier,
        "one external reset clock carrier resolved from every explicit flip-flop clock"};
  } catch (const std::exception& error) {
    return unsupported(std::string("automatic reset clock discovery failed: ") + error.what());
  }
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
