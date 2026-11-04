// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "BoolExpr.h"
#include "SNLTruthTableTree.h"
#include <tbb/concurrent_unordered_map.h>

namespace KEPLER_FORMAL {

/// A primitive output reads a value the Boolean cone cannot hold: a
/// VECTOR_SELECT whose offset can leave its vector reads X at that output.
/// Thrown when no unknown-read variable base is set (see Tree2BoolExpr).
class UnknownValueReadError : public std::runtime_error {
 public:
  UnknownValueReadError(const std::string& what, naja::DNL::DNLID isoID)
      : std::runtime_error(what), isoID_(isoID) {}
  naja::DNL::DNLID getIsoID() const { return isoID_; }

 private:
  naja::DNL::DNLID isoID_ = naja::DNL::DNLID_MAX;
};

/// Convert a truth-table tree directly into a BoolExpr
class Tree2BoolExpr {
 public:
  static BoolExpr* convert(const SNLTruthTableTree& tree,
                                           const std::vector<size_t>& varNames);
  static tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> iso2boolExpr_;

  // Unknown reads.  A VECTOR_SELECT output bit whose offset can leave the
  // vector reads X there, which naja reports as no value.  With a variable
  // base set, that output's iso reads the variable `base + iso`, and the
  // clause builder publishes every (iso, variable) pair so SEC models it as
  // permanently-X state, like an X literal (docs/sec-sequential-models.md).
  // Without a base, converting such a cone throws UnknownValueReadError and
  // the clause builder skips the output.
  static void setUnknownReadVarBase(std::optional<size_t> base);
  static std::optional<size_t> getUnknownReadVarBase();
  // The variable for an unknown read at `isoID`, recording the pair; throws
  // UnknownValueReadError when no base is set or the output drives no net.
  static size_t unknownReadVar(naja::DNL::DNLID isoID, const std::string& what);
  // The (iso, variable) pairs recorded since the base was last set, in iso
  // order; clears them.
  static std::vector<std::pair<naja::DNL::DNLID, size_t>> takeUnknownReadVars();
};

}  // namespace KEPLER_FORMAL
