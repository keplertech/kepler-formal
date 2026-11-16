// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BoolExpr.h"
#include "SNLTruthTableTree.h"
#include <tbb/concurrent_unordered_map.h>

namespace KEPLER_FORMAL {

/// Convert a truth-table tree directly into a BoolExpr
class Tree2BoolExpr {
 public:
  static BoolExpr* convert(const SNLTruthTableTree& tree,
                                           const std::vector<size_t>& varNames);
  static tbb::concurrent_unordered_map<naja::DNL::DNLID, BoolExpr*> iso2boolExpr_;

  // Unknown reads.  A VECTOR_SELECT output bit whose offset can leave the
  // vector reads X there, which naja reports as no value.  That output's iso
  // reads the variable `base + iso`, and the pair is recorded with the reason.
  // A clause builder that models X publishes the pairs, so SEC models each as
  // permanently-X state, like an X literal (docs/sec-sequential-models.md); a
  // builder without X skips the outputs whose expressions read such a
  // variable.  The clause builder sets the base for the duration of a build.
  static void setUnknownReadVarBase(std::optional<size_t> base);
  // The variable for an unknown read at `isoID`, recording the pair and
  // `what` reads X there.  A cone is only converted inside a clause build and
  // a primitive output inside a cone drives a net, so a call without a base or
  // without an iso is an internal error.
  static size_t unknownReadVar(naja::DNL::DNLID isoID, const std::string& what);
  // The (iso, variable) pairs recorded since the base was last set, in iso
  // order; clears them.
  static std::vector<std::pair<naja::DNL::DNLID, size_t>> takeUnknownReadVars();
  // What reads X at each of those isos; clears them.
  static std::unordered_map<naja::DNL::DNLID, std::string>
  takeUnknownReadReasons();
};

}  // namespace KEPLER_FORMAL
