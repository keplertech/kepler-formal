// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include "KeplerFormalDriver.h"
#include "Config.h"
#include "strategy/SequentialEquivalenceStrategy.h"

namespace KEPLER_FORMAL {

enum class BorrowedVerificationMode { LEC, SEC };

// File loading, preprocessing, compact mode and scope cleaning deliberately
// have no representation here: none may take ownership of a borrowed design.
struct BorrowedDesignOptions {
  BorrowedVerificationMode mode = BorrowedVerificationMode::LEC;
  Config::SolverType solver = Config::SolverType::KISSAT;
  size_t maxK = 32;
  SEC::SecEngine secEngine = SEC::SecEngine::Pdr;
  SEC::SecEncoding secEncoding = SEC::SecEncoding::DualRailSteady;
  bool allowBoundaryMismatch = false;
  bool reportSkippedOutputs = false;
  std::string logFile;
  std::string logLevel;
};

// Verify two live designs in the current Naja universe, without cloning,
// serializing, deleting or editing their netlists. The binding must validate
// runtime identity and retain the design wrappers; the caller must prevent
// concurrent Naja access (including reads/reset) for the synchronous call. This entry
// point is serialized and non-reentrant, but does not lock arbitrary Naja APIs.
// Temporary DNL, ordering IDs, top selections and expression caches are scoped
// to this call. Engine errors are returned as Error/exitCode=1 with a reason.
int verifyBorrowedDesigns(naja::NL::SNLDesign* design0,
                         naja::NL::SNLDesign* design1,
                         const BorrowedDesignOptions& options,
                         RunResult& result);

}  // namespace KEPLER_FORMAL
