// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "DNL.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL {

using BoundaryPair = std::pair<std::string, std::string>;
using BoundaryPairs = std::vector<BoundaryPair>;

// A virtual verification port. No corresponding SNL port is created.
struct BoundaryPort {
  size_t pairIndex = 0;
  std::string pinName;
  int32_t bit = 0;
  bool isInput = false;
  size_t width = 1;
  int32_t msb = 0;
  int32_t lsb = 0;
  std::string topTermName;
};

class BoundarySelection {
 public:
  BoundarySelection(naja::NL::SNLDesign* top,
                    const BoundaryPairs& pairs,
                    size_t side);
  const std::vector<BoundaryPort>& getPorts() const { return ports_; }

 private:
  std::vector<BoundaryPort> ports_;
};

// Read-only metadata for selected leaf pins in the original DNL.
class LeafBoundary {
 public:
  using DNLID = naja::DNL::DNLID;
  LeafBoundary(const naja::DNL::DNLFull& dnl,
               const BoundaryPairs& pairs, size_t side);
  const std::vector<BoundaryPort>& getPorts() const { return ports_; }
  const std::vector<DNLID>& getInputs() const { return inputs_; }
  const std::vector<DNLID>& getOutputs() const { return outputs_; }
  bool containsInstance(DNLID id) const { return instances_.count(id) != 0; }
  const BoundaryPort* getPort(DNLID id) const;
  bool isInput(DNLID id) const {
    const auto* port = getPort(id);
    return port && !port->isInput;
  }
 private:
  std::vector<BoundaryPort> ports_;
  std::unordered_map<DNLID, size_t> portIndices_;
  std::vector<DNLID> inputs_, outputs_;
  std::unordered_set<DNLID> instances_;
};

// Check that the two independently selected interfaces expose the same
// boundary pin bits, directions, and bus shapes.  Throws std::invalid_argument
// on the first mismatch.
void validateBoundaryInterfaces(const std::vector<BoundaryPort>& left,
                                const std::vector<BoundaryPort>& right);

}  // namespace KEPLER_FORMAL
