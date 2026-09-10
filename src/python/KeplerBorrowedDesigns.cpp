// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerBorrowedDesigns.h"

#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "BoolExprCache.h"
#include "DNL.h"
#include "KeplerFormalUtils.h"
#include "MiterStrategy.h"
#include "NLDB.h"
#include "NLUniverse.h"
#include "SNLBitTerm.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "Tree2BoolExpr.h"

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace KEPLER_FORMAL {
namespace {

// NLUniverse::setTopDesign also changes the owning DB's top selection. Save
// those selections separately, including absent tops and an absent top DB.
class BorrowedNajaState {
 public:
  BorrowedNajaState(naja::NL::NLUniverse* universe,
                    naja::NL::SNLDesign* first,
                    naja::NL::SNLDesign* second)
      : universe_(universe), topDB_(universe->getTopDB()) {
    for (auto* db : universe->getDBs()) {
      dbTops_.emplace_back(db, db->getTopDesign());
    }
    std::unordered_set<naja::NL::SNLDesign*> visited;
    std::vector<naja::NL::SNLDesign*> pending{first, second};
    while (!pending.empty()) {
      auto* design = pending.back();
      pending.pop_back();
      if (!visited.insert(design).second) {
        continue;
      }
      for (auto* term : design->getBitTerms()) {
        termOrderIDs_.emplace_back(term, term->getOrderID());
      }
      for (auto* instance : design->getInstances()) {
        instanceOrderIDs_.emplace_back(instance, instance->getOrderID());
        pending.push_back(instance->getModel());
      }
    }
    // The exchange does not lazily construct a graph, unlike DNL::get(). Do
    // this last so allocation failures above leave the caller untouched.
    dnl_ = naja::DNL::exchange(nullptr);
  }

  ~BorrowedNajaState() {
    naja::DNL::destroy();
    for (const auto& entry : termOrderIDs_) {
      entry.first->setOrderID(entry.second);
    }
    for (const auto& entry : instanceOrderIDs_) {
      entry.first->setOrderID(entry.second);
    }
    for (const auto& entry : dbTops_) {
      entry.first->setTopDesign(entry.second);
    }
    universe_->setTopDB(topDB_);
    naja::DNL::exchange(dnl_);
  }

  BorrowedNajaState(const BorrowedNajaState&) = delete;
  BorrowedNajaState& operator=(const BorrowedNajaState&) = delete;

 private:
  naja::NL::NLUniverse* universe_;
  naja::NL::NLDB* topDB_;
  naja::DNL::DNLFull* dnl_ = nullptr;
  std::vector<std::pair<naja::NL::NLDB*, naja::NL::SNLDesign*>> dbTops_;
  std::vector<std::pair<naja::NL::SNLBitTerm*, naja::NL::NLID::DesignObjectID>>
      termOrderIDs_;
  std::vector<std::pair<naja::NL::SNLInstance*, naja::NL::NLID::DesignObjectID>>
      instanceOrderIDs_;
};

class BorrowedExpressionState {
 public:
  BorrowedExpressionState() {
    previousIsoExpressions_.swap(Tree2BoolExpr::iso2boolExpr_);
  }
  ~BorrowedExpressionState() {
    Tree2BoolExpr::iso2boolExpr_.clear();
    previousIsoExpressions_.swap(Tree2BoolExpr::iso2boolExpr_);
  }

 private:
  BoolExprCache::ScopedContext expressions_;
  decltype(Tree2BoolExpr::iso2boolExpr_) previousIsoExpressions_;
};

class BorrowedRunState {
 public:
  BorrowedRunState()
      : solver_(Config::getSolverType()),
        reportSkipped_(Config::getReportSkippedPOs()),
        defaultLogger_(spdlog::default_logger()),
        miterLogger_(spdlog::get("miter_logger")),
        fallbackLogger_(spdlog::get("miter_logger_fallback")),
        miterLogFile_(MiterStrategy::logFileName_) {}

  ~BorrowedRunState() {
    MiterStrategy::cleanupProcessState();
    MiterStrategy::logFileName_ = std::move(miterLogFile_);
    Config::setSolverType(solver_);
    Config::setReportSkippedPOs(reportSkipped_);
    restoreLogger("miter_logger", miterLogger_);
    restoreLogger("miter_logger_fallback", fallbackLogger_);
    spdlog::set_default_logger(defaultLogger_);
  }

 private:
  static void restoreLogger(const char* name,
                            const std::shared_ptr<spdlog::logger>& previous) {
    if (spdlog::get(name) != previous) {
      spdlog::drop(name);
      if (previous) {
        spdlog::register_logger(previous);
      }
    }
  }

  Config::SolverType solver_;
  bool reportSkipped_;
  std::shared_ptr<spdlog::logger> defaultLogger_;
  std::shared_ptr<spdlog::logger> miterLogger_;
  std::shared_ptr<spdlog::logger> fallbackLogger_;
  std::string miterLogFile_;
};

void configureLogger(const BorrowedDesignOptions& options, RunResult& result) {
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
  if (options.mode == BorrowedVerificationMode::SEC && !options.logFile.empty()) {
    const auto parent = std::filesystem::path(options.logFile).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
    sinks.push_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(
        options.logFile, true));
    result.logFile = options.logFile;
  }
  auto logger = std::make_shared<spdlog::logger>(
      "kepler_formal_borrowed_logger", sinks.begin(), sinks.end());
  logger->set_level(options.logLevel == "debug" ? spdlog::level::debug
                                               : spdlog::level::info);
  logger->flush_on(spdlog::level::info);
  spdlog::set_default_logger(logger);
}

void assignSecResult(const SEC::SequentialEquivalenceResult& proof,
                     RunResult& result) {
  result.bound = proof.bound;
  result.reason = proof.reason;
  result.coveredOutputs = proof.coveredOutputs;
  result.totalOutputs = proof.totalOutputs;
  result.skippedObservedOutputs = proof.skippedObservedOutputs;
  if (proof.proofProgress) {
    result.provenOutputs = proof.proofProgress->provenOutputs;
    for (const auto& output : proof.proofProgress->unprovenOutputs) {
      result.unprovenOutputs.push_back(output.name);
    }
  } else if (proof.status == SEC::SequentialEquivalenceStatus::Equivalent ||
             proof.status == SEC::SequentialEquivalenceStatus::PartiallyProved) {
    result.provenOutputs = proof.coveredOutputs;
  }
  switch (proof.status) {
    case SEC::SequentialEquivalenceStatus::Equivalent:
      result.status = RunStatus::Equivalent;
      result.exitCode = kSecProvedExitCode;
      break;
    case SEC::SequentialEquivalenceStatus::PartiallyProved:
      result.status = RunStatus::PartiallyProved;
      result.exitCode = kSecPartiallyProvedExitCode;
      break;
    case SEC::SequentialEquivalenceStatus::Different:
      result.status = RunStatus::Different;
      result.exitCode = kSecCounterexampleExitCode;
      break;
    case SEC::SequentialEquivalenceStatus::Inconclusive:
      result.status = RunStatus::Inconclusive;
      result.exitCode = kSecInconclusiveExitCode;
      break;
    case SEC::SequentialEquivalenceStatus::Unsupported:
      result.status = RunStatus::Unsupported;
      result.exitCode = kSecInconclusiveExitCode;
      break;
  }
}

}  // namespace

int verifyBorrowedDesigns(naja::NL::SNLDesign* design0,
                         naja::NL::SNLDesign* design1,
                         const BorrowedDesignOptions& options,
                         RunResult& result) {
  result = RunResult{};
  result.inputFormat = "naja_design";
  result.verification = options.mode == BorrowedVerificationMode::SEC ? "sec" : "lec";
  static std::mutex mutex;
  static thread_local bool inProgress = false;
  if (inProgress) {
    result.reason = "Borrowed-design verification is not reentrant";
    return result.exitCode;
  }
  struct ReentrancyGuard {
    bool& active;
    explicit ReentrancyGuard(bool& flag) : active(flag) { active = true; }
    ~ReentrancyGuard() { active = false; }
  } reentrancyGuard(inProgress);
  std::lock_guard<std::mutex> lock(mutex);

  try {
    auto* universe = naja::NL::NLUniverse::get();
    if (!universe || !design0 || !design1) {
      throw std::invalid_argument("Borrowed designs need a live Naja universe and two live designs");
    }
    for (auto* design : {design0, design1}) {
      if (universe->getSNLDesign(design->getReference()) != design) {
        throw std::invalid_argument("Borrowed design does not belong to the active Naja runtime");
      }
    }
    BorrowedNajaState najaState(universe, design0, design1);
    Config::ScopedVerificationContext verificationContext;
    BorrowedExpressionState expressionState;
    BorrowedRunState runState;
    Config::setSolverType(options.solver);
    Config::setReportSkippedPOs(options.reportSkippedOutputs);
    configureLogger(options, result);
    if (options.mode == BorrowedVerificationMode::SEC) {
      SEC::SequentialEquivalenceStrategy strategy(
          design0, design1, options.solver, options.secEngine, options.secEncoding);
      const auto proof = strategy.run(options.maxK);
      assignSecResult(proof, result);
      if (options.reportSkippedOutputs) {
        writeBoundaryTermsReport("boundary_terms.txt", proof.extractedBoundaryReports);
        writeResetUnanchoredSkippedOutputsReport(
            "skipped_reset_unanchored_pos.txt", proof.resetUnanchoredSkippedOutputs);
        writeMultiClockDomainSkippedOutputsReport(
            "skipped_multi_clock_domain_pos.txt", proof.multiClockDomainSkippedOutputs);
        writeOpaqueCellSkippedOutputsReport(
            "skipped_opaque_cells_pos.txt", proof.opaqueCellSkippedOutputs);
      }
      SPDLOG_INFO("Borrowed SEC result: {} at k = {} ({}/{} outputs covered)",
                  runStatusName(result.status), result.bound,
                  result.coveredOutputs, result.totalOutputs);
      if (!result.reason.empty()) {
        SPDLOG_INFO("{}", result.reason);
      }
    } else {
      MiterStrategy strategy(design0, design1, options.logFile);
      strategy.setAllowBoundaryMismatch(options.allowBoundaryMismatch);
      strategy.init();
      result.logFile = MiterStrategy::getActualLogFileName();
      result.status = strategy.run(false) ? RunStatus::Equivalent : RunStatus::Different;
      // Match the CLI's LEC exit convention: a completed comparison returns 0
      // even for Different; the structured status carries the verdict.
      result.exitCode = 0;
    }
  } catch (const std::exception& error) {
    result.status = RunStatus::Error;
    result.exitCode = 1;
    result.reason = error.what();
  } catch (...) {
    result.status = RunStatus::Error;
    result.exitCode = 1;
    result.reason = "Unknown exception during borrowed-design verification";
  }
  return result.exitCode;
}

}  // namespace KEPLER_FORMAL
