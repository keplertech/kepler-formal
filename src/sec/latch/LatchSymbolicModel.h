// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "BoolExpr.h"
#include "latch/LatchEventModel.h"

namespace KEPLER_FORMAL::SEC::LATCH {

using SymbolicBits = std::vector<BoolExpr*>;
using FreshSymbol = std::function<BoolExpr*()>;

struct SymbolicReaction {
  SymbolicBits storage, outputs;
  BoolExpr* error = BoolExpr::createFalse();
};

// These callbacks have the same pin-local contract as Primitive, but construct
// Boolean DAGs. They must be pure and safe for concurrent evaluation.
struct SymbolicPrimitive {
  std::function<SymbolicReaction(const SymbolicBits&, const SymbolicBits&,
      const SymbolicBits&, std::optional<size_t>, bool)> react;
  std::function<SymbolicBits(const SymbolicBits&, const SymbolicBits&)> initialOutputs;
};

struct SymbolicNetwork {
  Network reference;
  std::vector<SymbolicPrimitive> primitives;
};

struct SymbolicState {
  SymbolicBits current, previous;
  std::vector<SymbolicBits> storage;
  SymbolicBits active;
  BoolExpr* bootstrap = BoolExpr::createFalse();
  BoolExpr* error = BoolExpr::createFalse();
};

BoolExpr* symbolicMux(BoolExpr* condition, BoolExpr* yes, BoolExpr* no);
SymbolicBits flattenSymbolicBoundary(const SymbolicState& state);

// Total symbolic counterpart of EventModel. Choice encodings must cover every
// pin ordering, with unused codes selecting a legal ordering (not pruning).
// Fresh symbols are requested serially in stable primitive order, before any
// parallel work. Wave boundaries retain shared producer choices across fanout.
class SymbolicEventModel {
 public:
  explicit SymbolicEventModel(SymbolicNetwork network, Limits limits = {}, size_t workers = 0);
  const SymbolicNetwork& network() const { return network_; }
  SymbolicState boundary(const SymbolicBits& current,
                         const std::vector<SymbolicBits>& storage) const;
  SymbolicState bootstrap(const SymbolicBits& inputs,
      const std::vector<SymbolicBits>& storage, const SymbolicBits& seeds) const;
  SymbolicState admit(const SymbolicState& boundary, const SymbolicBits& inputs) const;
  // For reference/proof construction fresh must allocate globally fresh Boolean
  // variables. After independent choice-independence certification, constants
  // may instead select one legal ordering for a deterministic compiled result.
  SymbolicState wave(const SymbolicState& state, const FreshSymbol& fresh) const;
  BoolExpr* stable(const SymbolicState& state) const;
  // A candidate boundary invariant: all primitive output and level/async rules
  // are consistent. The compiler must prove BOOT establishes it and episodes
  // preserve it; it is not an assumed reachability restriction.
  BoolExpr* boundaryInvariant(const SymbolicState& state) const;
 private:
  SymbolicNetwork network_;
  Limits limits_;
  size_t workers_;
  std::vector<std::vector<size_t>> consumers_;
};

// Exact truth-table lifting for small test/reference primitives, not the
// production frontend for wide gates. The latter supplies native DAG callbacks.
SymbolicPrimitive liftPrimitive(const Primitive& primitive, size_t maxLocalBits = 12);

}  // namespace KEPLER_FORMAL::SEC::LATCH
