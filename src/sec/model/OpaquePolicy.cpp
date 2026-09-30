// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "model/OpaquePolicy.h"

#include <algorithm>
#include <optional>

#include "../../config/Config.h"
#include "../../utils/DesignBoundary.h"
#include "SNLDesign.h"
#include "SNLPath.h"
#include "common/SignalKey.h"
#include "model/SequentialDesignModel.h"

namespace KEPLER_FORMAL::SEC {

void recordOpaqueOutputlessCells(SequentialDesignModel& model,
                                 const naja::DNL::DNLFull& dnl,
                                 const LeafBoundary* boundary) {
  if (!Config::getErrorOnOpaque()) return;
  for (const auto leafID : dnl.getLeaves()) {
    const auto& instance = dnl.getDNLInstanceFromID(leafID);
    if (instance.isTop() || (boundary && boundary->containsInstance(leafID))) continue;
    const auto* cell = instance.getSNLModel();
    if (!cell || !cell->isLeaf()) continue;
    bool hasOutput = false;
    for (auto* term : cell->getBitTerms()) {
      if (term->getDirection() != naja::NL::SNLTerm::Direction::Input) {
        hasOutput = true;
        break;
      }
    }
    if (hasOutput) continue;
    SignalKey key;
    for (const auto& part : instance.getPath().getPathNames()) {
      key.first.push_back(part.getID());
    }
    key.first.push_back(uint64_t{1} << 60);
    key.second.push_back(0);
    const auto path = instance.getFullPath();
    model.displayNameByKey.insert_or_assign(key, path);
    model.connectivitySkipInfoByKey.insert_or_assign(key,
        ConnectivitySkipInfo{ConnectivitySkipOrigin::OpaqueInternal,
            "opaque outputless cell `" + path + "` (model `" +
            cell->getName().getString() + "`): no usable SEC output model"});
  }
}

void applyOpaquePolicy(SequentialDesignModel& model,
                       const std::string& topName, size_t side) {
  if (!Config::getErrorOnOpaque()) {
    return;
  }
  std::optional<std::string> first;
  for (const auto& [key, info] : model.connectivitySkipInfoByKey) {
    if (info.origin != ConnectivitySkipOrigin::OpaqueInternal) {
      continue;
    }
    const auto name = model.displayNameByKey.find(key);
    const auto signal = name == model.displayNameByKey.end()
        ? signalKeyToString(key) : name->second;
    const std::string diagnostic =
        "SEC error-on-opaque: design " + std::to_string(side + 1) +
        " (`" + topName + "`), signal `" + signal + "`: " + info.detail;
    // Hash-map traversal and extraction worker ordering must not determine
    // which diagnostic the user sees first.
    if (!first || diagnostic < *first) {
      first = diagnostic;
    }
  }
  if (first && std::find(model.unsupportedReasons.begin(),
                         model.unsupportedReasons.end(), *first) ==
                   model.unsupportedReasons.end()) {
    model.unsupportedReasons.push_back(*first);
  }
}

}  // namespace KEPLER_FORMAL::SEC
