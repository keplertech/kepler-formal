// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/SignalKey.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL::SEC::PHASE {

using ExprID = uint32_t;
constexpr ExprID InvalidExpr = std::numeric_limits<ExprID>::max();

enum class TernaryValue : uint8_t {
  Zero,
  One,
  Unknown,
};

enum class ExpressionOp : uint8_t {
  Constant,
  Variable,
  Not,
  And,
  Or,
  Xor,
};

struct ExpressionNode {
  ExpressionOp op = ExpressionOp::Constant;
  size_t variable = 0;
  bool constant = false;
  ExprID left = InvalidExpr;
  ExprID right = InvalidExpr;
};

// An owned, immutable-after-construction Boolean DAG.  The SNL frontend copies
// every temporary BoolExpr into this representation before releasing DNL, so
// the phase analysis is independent of process-global frontend caches.
class ExpressionArena {
 public:
  ExpressionArena();

  ExprID constant(bool value) const { return value ? 1u : 0u; }
  ExprID variable(size_t variableIndex);
  ExprID makeNot(ExprID operand);
  ExprID makeAnd(ExprID lhs, ExprID rhs);
  ExprID makeOr(ExprID lhs, ExprID rhs);
  ExprID makeXor(ExprID lhs, ExprID rhs);
  ExprID makeIte(ExprID condition, ExprID whenTrue, ExprID whenFalse);

  const ExpressionNode& node(ExprID id) const { return nodes_.at(id); }
  const std::vector<ExpressionNode>& nodes() const { return nodes_; }
  bool isValid(ExprID id) const { return id < nodes_.size(); }

 private:
  ExprID appendBinary(ExpressionOp op, ExprID lhs, ExprID rhs);

  std::vector<ExpressionNode> nodes_;
  std::unordered_map<size_t, ExprID> variableNodes_;
};

enum class VariableKind : uint8_t {
  Input,
  State,
};

struct NormalizedVariable {
  SignalKey key;
  std::string displayName;
  VariableKind kind = VariableKind::Input;
};

struct NormalizedStateBit {
  size_t variable = 0;
  TernaryValue initialValue = TernaryValue::Unknown;
  ExprID nextState = InvalidExpr;
};

enum class ClearPresetValue : uint8_t {
  Zero,
  One,
  Hold,
  Toggle,
  Unknown,
};

struct NormalizedStateOutput {
  SignalKey key;
  std::string displayName;
  // False is Q / State(i); true is QN / Not(State(i)).
  bool complemented = false;
};

struct NormalizedLatch {
  size_t stateIndex = 0;
  std::string instancePath;
  std::string modelName;
  ExprID data = InvalidExpr;
  // Canonical active-high transparency predicate. Active-low gate polarity is
  // represented by a Not node in this expression rather than a separate clock.
  ExprID transparency = InvalidExpr;
  std::optional<ExprID> clear;
  std::optional<ExprID> preset;
  ClearPresetValue simultaneousClearPreset = ClearPresetValue::Unknown;
  // All modeled physical outputs for this state, including Q and QN aliases.
  std::vector<NormalizedStateOutput> outputs;
};

struct NormalizedOutput {
  SignalKey key;
  std::string displayName;
  ExprID expression = InvalidExpr;
};

struct NormalizationDiagnostic {
  std::string object;
  std::string message;
  bool fatal = false;
};

struct NormalizationOptions {
  // Naja sequential models do not carry reset-state values.  Callers that
  // have an environment/reset protocol may provide the resulting physical
  // output value here, keyed by the collector's display name.  When a model
  // exposes only a complemented state output, collection converts that value
  // back to the logical Naja model state.
  std::unordered_map<std::string, TernaryValue> initialValueByDisplayName;
  bool collectObservedOutputs = true;
};

class NormalizedTransitionSystem {
 public:
  size_t addInput(SignalKey key, std::string displayName);
  size_t addState(SignalKey key, std::string displayName,
                  TernaryValue initialValue = TernaryValue::Unknown);
  void setNextState(size_t stateIndex, ExprID expression);
  size_t addLatch(NormalizedLatch latch);
  size_t addOutput(NormalizedOutput output);
  void addDiagnostic(NormalizationDiagnostic diagnostic);

  ExpressionArena& expressions() { return expressions_; }
  const ExpressionArena& expressions() const { return expressions_; }
  const std::vector<NormalizedVariable>& variables() const {
    return variables_;
  }
  const std::vector<size_t>& inputVariables() const { return inputVariables_; }
  const std::vector<NormalizedStateBit>& states() const { return states_; }
  const std::vector<NormalizedLatch>& latches() const { return latches_; }
  const std::vector<NormalizedOutput>& outputs() const { return outputs_; }
  const std::vector<NormalizationDiagnostic>& diagnostics() const {
    return diagnostics_;
  }

  const NormalizedVariable& stateVariable(size_t stateIndex) const;
  ExprID variableExpression(size_t variableIndex) {
    return expressions_.variable(variableIndex);
  }
  bool hasFatalDiagnostics() const;
  bool validate(std::string* reason = nullptr) const;

 private:
  ExpressionArena expressions_;
  std::vector<NormalizedVariable> variables_;
  std::vector<size_t> inputVariables_;
  std::vector<NormalizedStateBit> states_;
  std::vector<NormalizedLatch> latches_;
  std::vector<NormalizedOutput> outputs_;
  std::vector<NormalizationDiagnostic> diagnostics_;
};

// Collect the complete implicitly-clocked transition representation required
// by Section II of Bjesse/Kukula. The frontend is serialized internally
// because NLUniverse and DNL are global; all returned data is owned and safe to
// analyze concurrently afterwards. Unsupported state or logic is diagnosed as
// fatal instead of being abstracted into an unconstrained input.
NormalizedTransitionSystem collectNormalizedTransitionSystem(
    naja::NL::SNLDesign* top, const NormalizationOptions& options = {});

const char* toString(TernaryValue value);

}  // namespace KEPLER_FORMAL::SEC::PHASE
