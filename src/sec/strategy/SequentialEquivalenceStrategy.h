// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
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
  // Memories the proof compared at their interface instead of cell by cell.
  // Their interface pins are proved equal and left out of the output counts.
  std::vector<std::string> abstractedMemories;
  // Why the interface comparison could not carry the proof. The status of such
  // a result is Inconclusive, never a counterexample; run() then repeats the
  // proof on bit-level memory models and keeps this text on the final result.
  std::string memoryAbstractionFailure;
  // Set when the counterexample of a binary proof on bit-level memories is
  // not a binary-defined difference: it depends on uninitialized state, such
  // as memory cells that were never written. The status is then Inconclusive
  // and `reason` holds the trace.
  bool counterexampleDependsOnUninitializedState = false;

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
struct SequentialExtractOptions;

enum class SecProgressLevel {
  Info,
  Warning,
};

// A memory with more bits than this is reported with a warning when it is
// modeled bit by bit: one state variable per bit does not scale.
inline constexpr size_t kLargeBitLevelMemoryBits = size_t{1} << 16;

// Builds a combined SEC problem from two sequential designs and discharges it
// with the selected SEC proof engine over the extracted transition system.
class SequentialEquivalenceStrategy {
 public:
  // Receives the progress lines of a run: model sizes, memory handling and
  // the start of each long phase. Nothing is reported without a logger.
  using ProgressLogger =
      std::function<void(SecProgressLevel, const std::string&)>;
  // Extracts one design, `side` 0 or 1, with the given options. Compact SEC
  // loads and releases the design inside this call; identical designs may
  // return the same model for both sides.
  using ModelExtractor =
      std::function<std::shared_ptr<const SequentialDesignModel>(
          size_t side, const SequentialExtractOptions& options)>;

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
  // Memory abstraction: a memory found under the same path and with the same
  // primitive in both designs is compared at its interface first; its cells
  // are modeled only when that comparison does not carry the proof. Starts
  // from Config::getSecMemoryAbstraction(), which is on by default.
  void setMemoryAbstraction(bool enabled) { memoryAbstraction_ = enabled; }
  void setProgressLogger(ProgressLogger logger) {
    progressLogger_ = std::move(logger);
  }
  // Same as run(), with the caller extracting the models.
  SequentialEquivalenceResult runWithModelExtractor(
      const ModelExtractor& extractor,
      size_t maxK) const;
  // Proves the two models as given. Models extracted with abstracted memories
  // are accepted when their memories pair up; a comparison that cannot carry
  // the proof returns Inconclusive with memoryAbstractionFailure set.
  SequentialEquivalenceResult runExtractedModels(
      const SequentialDesignModel& model0,
      const SequentialDesignModel& model1,
      size_t maxK) const;
  // The reset port of `resetSpec` that names no input of `top`, as the error
  // text SEC reports for it. Lets a caller reject a mistyped reset port right
  // after parsing, before any model is built.
  static std::optional<std::string> findUnknownResetPort(
      const naja::NL::SNLDesign* top,
      const SecResetSpec& resetSpec);

 private:
  void logProgress(SecProgressLevel level, const std::string& message) const;

  naja::NL::SNLDesign* top0_;
  naja::NL::SNLDesign* top1_;
  KEPLER_FORMAL::Config::SolverType solverType_;
  SecEngine secEngine_;
  SecEncoding encoding_;
  SecResetSpec resetSpec_;
  Btor2ExportOptions exportOptions_;
  BoundaryPairs boundaryPairs_;
  InternalRelationOptions internalRelationOptions_;
  bool memoryAbstraction_ = true;
  ProgressLogger progressLogger_;
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
