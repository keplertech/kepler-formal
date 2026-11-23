// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "../../config/Config.h"
#include "../../utils/DesignBoundary.h"
#include "export/Btor2ExportOptions.h"
#include "proof/InternalRelations.h"

namespace naja::NL {
class SNLDesign;
}

namespace KEPLER_FORMAL::SEC {

enum class SecEngine {
  KInduction,
  Imc,
  Pdr,
};

enum class SecEncoding {
  Binary,
  DualRailSteady,
};

enum class SequentialEquivalenceStatus {
  Equivalent,
  PartiallyProved,
  Different,
  Inconclusive,
  Unsupported,
  Exported,
};

struct SecResetPortSpec {
  std::string name;
  bool activeValue = true;
};

struct SecResetSpec {
  size_t cycles = 0;
  std::vector<SecResetPortSpec> ports;

  bool enabled() const {
    return cycles != 0 && !ports.empty();
  }
};

struct ExtractedBoundaryReportEntry {  // LCOV_EXCL_LINE
  std::string design;
  std::string signal;
  std::vector<std::string> roles;
  std::string connectivitySkip;
};

struct SequentialEquivalenceUnprovenOutput { // LCOV_EXCL_LINE
  size_t index = 0;
  std::string name;
};

struct SequentialEquivalenceProofProgress { // LCOV_EXCL_LINE
  std::string engineLabel;
  size_t provenOutputs = 0;
  size_t totalOutputs = 0;
  std::vector<SequentialEquivalenceUnprovenOutput> unprovenOutputs;
};

// The state bits of one design that a dual-rail proof may still hold as X when
// the outputs are first compared.
struct SequentialEquivalenceUnknownState {  // LCOV_EXCL_LINE
  size_t unknownBits = 0;
  size_t stateBits = 0;
  // The first few of them, by name.
  std::vector<std::string> examples;
};

// What a dual-rail proof starts from. The proof takes every state bit at its
// boot value, X where the design gives none, through the reset window; a bit
// the window cannot reach with a value stays X. Empty in the binary encoding,
// which skips the outputs reading such state.
struct SequentialEquivalenceResetReport {  // LCOV_EXCL_LINE
  // Zero when no reset window was given.
  size_t resetCycles = 0;
  // The binary passes (binary_relation_round) were allowed to run. They start
  // matched unknown registers of the two designs from one common value.
  bool binaryPasses = false;
  std::array<SequentialEquivalenceUnknownState, 2> unknownState;

  bool leavesUnknownState() const {
    return unknownState[0].unknownBits != 0 || unknownState[1].unknownBits != 0;
  }
};

struct SequentialEquivalenceResult {  // LCOV_EXCL_LINE
  SequentialEquivalenceStatus status = SequentialEquivalenceStatus::Unsupported;
  size_t bound = 0;
  std::string reason;
  size_t coveredOutputs = 0;
  size_t totalOutputs = 0;
  std::optional<SequentialEquivalenceProofProgress> proofProgress;
  std::vector<std::string> skippedObservedOutputs;
  std::vector<std::string> resetUnanchoredSkippedOutputs;
  std::vector<std::string> multiClockDomainSkippedOutputs;
  std::vector<std::string> opaqueCellSkippedOutputs;
  std::vector<ExtractedBoundaryReportEntry> extractedBoundaryReports;
  SequentialEquivalenceResetReport resetReport;

  double outputCoveragePercent() const {
    if (totalOutputs == 0) {
      // LCOV_EXCL_START
      return 0.0;
      // LCOV_EXCL_STOP
    }
    return (100.0 * static_cast<double>(coveredOutputs)) /
           static_cast<double>(totalOutputs);
  }
};

struct SequentialDesignModel;

// Builds a combined SEC problem from two sequential designs and discharges it
// with the selected SEC proof engine over the extracted transition system.
class SequentialEquivalenceStrategy {
 public:
  SequentialEquivalenceStrategy(
      naja::NL::SNLDesign* top0,
      naja::NL::SNLDesign* top1,
      KEPLER_FORMAL::Config::SolverType solverType =
          KEPLER_FORMAL::Config::getSolverType(),
      SecEngine secEngine = SecEngine::Pdr,
      SecEncoding encoding = SecEncoding::DualRailSteady,
      SecResetSpec resetSpec = {},
      Btor2ExportOptions exportOptions = {});

  SequentialEquivalenceResult run(size_t maxK) const;
  void setBoundaryPairs(const BoundaryPairs& pairs) { boundaryPairs_ = pairs; }
  void setInternalRelationOptions(InternalRelationOptions options) {
    internalRelationOptions_ = options;
  }
  SequentialEquivalenceResult runExtractedModels(
      const SequentialDesignModel& model0,
      const SequentialDesignModel& model1,
      size_t maxK) const;

 private:
  naja::NL::SNLDesign* top0_;
  naja::NL::SNLDesign* top1_;
  KEPLER_FORMAL::Config::SolverType solverType_;
  SecEngine secEngine_;
  SecEncoding encoding_;
  SecResetSpec resetSpec_;
  Btor2ExportOptions exportOptions_;
  BoundaryPairs boundaryPairs_;
  InternalRelationOptions internalRelationOptions_;
};

namespace detail {

SequentialEquivalenceProofProgress buildSecEngineProofProgress(
    const std::string& engineLabel,
    const std::vector<std::string>& observedOutputNames,
    size_t totalOutputCount,
    size_t provenOutputCount);

std::vector<std::string> buildSecEngineProofProgressDiagLines(
    const std::string& engineLabel,
    const std::vector<std::string>& observedOutputNames,
    size_t totalOutputCount,
    size_t provenOutputCount);

}  // namespace detail

}  // namespace KEPLER_FORMAL::SEC
