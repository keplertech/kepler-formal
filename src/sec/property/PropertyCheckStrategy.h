// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

#include "../../config/Config.h"
#include "BoolExpr.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace naja::NL {
class SNLDesign;
class SNLInstance;
}

namespace KEPLER_FORMAL::SEC {

enum class PropertyKind {
  Assert,
  Assume,  // also `restrict`
  Cover,
};

enum class PropertyStatus {
  Proved,        // assert: holds in every reachable step
  Failed,        // assert: a counterexample was found
  Covered,       // cover: a trace reaches it
  Unreachable,   // cover: proved impossible
  Assumed,       // assume: used as a constraint on every trace
  Inconclusive,  // neither a proof nor a trace within the bound
  Unsupported,   // outside the supported subset; never counted as proved
};

const char* propertyKindName(PropertyKind kind);
const char* propertyStatusName(PropertyStatus status);

struct PropertyResult {
  // Hierarchical path of the property from the top design.
  std::string name;
  // "file:line" when the frontend recorded it.
  std::string source;
  PropertyKind kind = PropertyKind::Assert;
  PropertyStatus status = PropertyStatus::Unsupported;
  // Proof bound, or the cycle of the failure or of the cover hit.
  size_t bound = 0;
  // Reason for an unsupported or inconclusive result, or the input trace.
  std::string detail;
};

struct PropertyCheckResult {
  std::vector<PropertyResult> properties;
  // Non-empty when the design itself could not be modelled; no property was
  // checked in that case.
  std::string error;

  size_t count(PropertyStatus status) const;
};

// What the strategy hands a translator for one captured property.
struct PropertyTranslationRequest {
  // The naja_formal_property instance that holds the property.
  const naja::NL::SNLInstance* instance = nullptr;
  // Value in the current step of every DATA bit of the instance, bit 0 first.
  std::vector<BoolExpr*> data;
  // Value `bit` had in the previous step. Each call adds one history register
  // that runs in every step, whatever the reset or disable condition.
  std::function<BoolExpr*(BoolExpr* bit)> previous;
};

// A captured property reduced to logic over its DATA bits.
struct TranslatedProperty {
  // Signal whose rising edge advances the property.
  BoolExpr* clock = nullptr;
  // 1 disables the property in the current step.
  BoolExpr* disable = nullptr;
  // Temporal structure over the leaves, as PropertyExpression.h describes.
  std::string expression;
  std::vector<BoolExpr*> leaves;
  // Non-empty when the property cannot be translated; the rest is unused.
  std::string unsupportedReason;
};

// A frontend captures a property without giving it a meaning. The translator
// supplies that meaning for the source language of the property.
using PropertyTranslator =
    std::function<TranslatedProperty(const PropertyTranslationRequest&)>;

// Checks the properties a frontend captured in the design as
// naja_formal_property instances. Each one is translated, compiled into a
// 1-bit signal plus monitor state and proved with the selected sequential
// engine.
//
// Model: one step is one rising edge of the single design clock, values are
// two-valued, and properties are checked from the first step after the reset
// window. Assumptions constrain every step of a trace up to and including the
// step under check.
class PropertyCheckStrategy {
 public:
  PropertyCheckStrategy(
      naja::NL::SNLDesign* top,
      PropertyTranslator translator,
      KEPLER_FORMAL::Config::SolverType solverType =
          KEPLER_FORMAL::Config::getSolverType(),
      SecEngine engine = SecEngine::Pdr,
      SecResetSpec resetSpec = {});

  PropertyCheckResult run(size_t maxK) const;

 private:
  naja::NL::SNLDesign* top_;
  PropertyTranslator translator_;
  KEPLER_FORMAL::Config::SolverType solverType_;
  SecEngine engine_;
  SecResetSpec resetSpec_;
};

}  // namespace KEPLER_FORMAL::SEC
