// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL {

using BoundaryPair = std::pair<std::string, std::string>;
using BoundaryPairs = std::vector<BoundaryPair>;

// Value-owned description of one promoted instance pin bit.  isInput refers
// to the selected instance pin: an input pin is promoted to a top output,
// while an output pin is promoted to a top input.
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

class BoundaryDesign {
 public:
  // The source design, its database, and its universe must outlive this
  // object.  A process-global DNL built from getTop() must be destroyed or
  // exchanged out before this object releases its scratch design.
  BoundaryDesign(naja::NL::SNLDesign* top,
                 const BoundaryPairs& pairs,
                 size_t side);
  ~BoundaryDesign();

  BoundaryDesign(const BoundaryDesign&) = delete;
  BoundaryDesign& operator=(const BoundaryDesign&) = delete;
  BoundaryDesign(BoundaryDesign&&) noexcept;
  BoundaryDesign& operator=(BoundaryDesign&&) noexcept;

  naja::NL::SNLDesign* getTop() const;
  const std::vector<BoundaryPort>& getPorts() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Check that the two independently transformed designs expose the same
// boundary pin bits, directions, and bus shapes.  Throws std::invalid_argument
// on the first mismatch.
void validateBoundaryInterfaces(const std::vector<BoundaryPort>& left,
                                const std::vector<BoundaryPort>& right);

}  // namespace KEPLER_FORMAL
