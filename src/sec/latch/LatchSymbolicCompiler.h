// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "latch/LatchSymbolicModel.h"
#include "latch/LatchSettlingCompiler.h"

namespace KEPLER_FORMAL::SEC::LATCH {

struct SymbolicCompileOptions {
  Bits initialInputs;
  std::vector<Bits> initialStorage;
  bool singleExternalInputChange = false;
  size_t maxWaves = 256;
  size_t maxNodes = 2000000;
  unsigned maxSatConflicts = 500000;
  unsigned maxSatDecisions = 5000000;
  size_t workers = 0;
};

struct SymbolicMacro {
  // Full boundary layout: current nets, then each primitive's storage. At a
  // boundary previous=current, active=BOOT=error=0, reconstructing full history.
  std::vector<size_t> stateSymbols, inputSymbols;
  Bits initialState;
  SymbolicBits nextState, observedNets;
  size_t bootstrapWaves = 0, transitionWaves = 0;
  // Retain the proved admission contract and its current-net projection. An
  // encoder must not broaden the environment or relabel remembered inputs.
  bool singleExternalInputChange = false;
  std::vector<size_t> externalInputNets;
};

struct SymbolicCompileResult {
  CertificationStatus status = CertificationStatus::Invalid;
  std::string detail;
  std::optional<SymbolicMacro> model;
  bool certified() const { return status == CertificationStatus::Certified && model.has_value(); }
};

// Prove universal progress and complete-boundary uniqueness with independent
// choice/seed copies, plus initialization and inductive boundary closure. Only
// then eliminate choices and retain K unfolded waves as deterministic logic.
// SAT/UNKNOWN over the candidate invariant is unproved, not a reachable defect.
SymbolicCompileResult compileSymbolicNetwork(const SymbolicNetwork& network,
                                            const SymbolicCompileOptions& options);

}  // namespace KEPLER_FORMAL::SEC::LATCH
