// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "BorrowedLatchOptions.h"

#include <limits>
#include <stdexcept>

namespace KEPLER_FORMAL {

SEC::LATCH::SupportOptions BorrowedLatchOptions::validated(
    bool isSec, bool hasLeafBoundaries) const {
  SEC::LATCH::SupportOptions result;
  if (!enabled) {
    if (inputChanges || initialInputs || initialStorage || workers || maxWaves ||
        maxStates || maxTransactions || maxSymbolicNodes || maxSatConflicts || maxSatDecisions) {
      throw std::invalid_argument("latch event tuning requires latch_support=true");
    }
    return result;
  }
  if (!isSec) throw std::invalid_argument("latch_support is only supported for SEC");
  if (!inputChanges || !initialInputs || !initialStorage) {
    throw std::invalid_argument(
        "latch_support requires explicit latch_input_changes, latch_initial_inputs and latch_initial_storage");
  }
  if (*inputChanges != LatchInputChanges::Any && *inputChanges != LatchInputChanges::Single) {
    throw std::invalid_argument("latch_input_changes must be any or single");
  }
  if (hasLeafBoundaries) {
    throw std::invalid_argument(
        "latch_support requires the complete top interface, not selected leaf boundaries");
  }
  if (workers && *workers > size_t(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("latch_workers must fit a nonnegative int");
  }
  if ((maxWaves && !*maxWaves) || (maxStates && !*maxStates) ||
      (maxTransactions && !*maxTransactions) || (maxSymbolicNodes && !*maxSymbolicNodes) ||
      (maxSatConflicts && !*maxSatConflicts) || (maxSatDecisions && !*maxSatDecisions)) {
    throw std::invalid_argument("latch resource limits must be positive integers");
  }
  if ((maxSatConflicts && *maxSatConflicts > std::numeric_limits<unsigned>::max()) ||
      (maxSatDecisions && *maxSatDecisions > std::numeric_limits<unsigned>::max())) {
    throw std::invalid_argument("latch SAT limits must fit a positive unsigned int");
  }
  result.enabled = true;
  result.singleInputChange = *inputChanges == LatchInputChanges::Single;
  result.initialInputs = initialInputs;
  result.initialStorage = initialStorage;
  if (workers) result.workers = *workers;
  if (maxWaves) result.limits.maxWaves = *maxWaves;
  if (maxStates) result.limits.maxBoundaryStates = *maxStates;
  if (maxTransactions) result.limits.maxTransactions = *maxTransactions;
  if (maxSymbolicNodes) result.maxSymbolicNodes = *maxSymbolicNodes;
  if (maxSatConflicts) result.maxSatConflicts = static_cast<unsigned>(*maxSatConflicts);
  if (maxSatDecisions) result.maxSatDecisions = static_cast<unsigned>(*maxSatDecisions);
  return result;
}

}  // namespace KEPLER_FORMAL
