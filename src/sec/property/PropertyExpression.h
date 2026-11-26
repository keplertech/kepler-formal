// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "BoolExpr.h"

namespace KEPLER_FORMAL::SEC {

// Parsed form of a property expression: either an atom or a
// parenthesized list, for example `(implies (leaf 0) (leaf 1))`.
struct PropertyTerm {
  bool isAtom = false;
  std::string atom;
  std::vector<PropertyTerm> items;

  // Operator name of a list, empty for an atom or a list that starts with a
  // nested list.
  const std::string& op() const;
};

// Throws std::invalid_argument when the text is not one well-formed term.
PropertyTerm parsePropertyTerm(const std::string& text);

// One state bit of a property monitor. Every monitor bit starts at zero and
// its next value is `next`; the caller adds the reset and disable gating.
struct PropertyMonitorBit {
  size_t var = 0;
  BoolExpr* next = nullptr;
};

struct CompiledProperty {
  // Assertion and assumption view: 1 in a step where a violation is detected.
  // Cover view: 1 in a step where a match completes.
  BoolExpr* signal = nullptr;
  std::vector<PropertyMonitorBit> monitor;
  // Non-empty when the property is outside the supported subset; signal and
  // monitor are then meaningless.
  std::string unsupportedReason;

  bool supported() const { return unsupportedReason.empty(); }
};

// Turns the temporal structure of a property into combinational logic over
// its Boolean leaves plus monitor state.
//
// Scope: safety properties with weak semantics on one clock, where one step
// of the model is one tick of that clock. A new evaluation attempt starts in
// every step. Sequences on the left of an implication are matched with one
// monitor bit per delay stage, so overlapping attempts share the bits. A
// sequence on the right of an implication is supported when its failure is
// decided without tracking alternatives: exact delays, plus a final bounded
// delay window. Everything else is reported through unsupportedReason, never
// approximated.
class PropertyCompiler {
 public:
  using VarAllocator = std::function<size_t()>;

  PropertyCompiler(std::vector<BoolExpr*> leaves, VarAllocator allocateVar);

  // Violation detector for assert, assume and restrict.
  CompiledProperty compileViolation(const PropertyTerm& property);
  // Match detector for cover.
  CompiledProperty compileCover(const PropertyTerm& property);

 private:
  struct Unsupported {
    std::string reason;
  };

  BoolExpr* leaf(const PropertyTerm& term) const;
  bool isBoolean(const PropertyTerm& term) const;
  BoolExpr* boolean(const PropertyTerm& term) const;

  BoolExpr* delayOne(BoolExpr* signal);
  BoolExpr* delayExact(BoolExpr* signal, size_t steps);
  BoolExpr* delayWindow(BoolExpr* signal, size_t low, size_t high, bool unbounded);

  BoolExpr* matchSequence(const PropertyTerm& sequence, BoolExpr* start);
  BoolExpr* failProperty(const PropertyTerm& property, BoolExpr* start);
  BoolExpr* failSequence(const PropertyTerm& sequence, BoolExpr* start);

  CompiledProperty finish(BoolExpr* signal);

  std::vector<BoolExpr*> leaves_;
  VarAllocator allocateVar_;
  std::vector<PropertyMonitorBit> monitor_;
};

}  // namespace KEPLER_FORMAL::SEC
