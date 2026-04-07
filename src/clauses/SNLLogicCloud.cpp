// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLLogicCloud.h"
#include <tbb/tbb_allocator.h>
#include <cassert>
#include <sstream>
#include "NajaProperty.h"
#include "SNLDesignModeling.h"
#include "SNLBitNet.h"
#include "NLDB0.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLScalarTerm.h"
#include "tbb/concurrent_vector.h"
#include "tbb/enumerable_thread_specific.h"
#include "SNLPath.h"
#include "Tree2BoolExpr.h"
#include "../config/Config.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <mutex>
#include <ostream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>



//#define DEBUG_CHECKS
//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

namespace {

bool shouldReportSkippedPOs() {
  return KEPLER_FORMAL::Config::getReportSkippedPOs();
}

bool canUseCachedIsoShortcut(const naja::DNL::DNLIso& iso,
                             naja::DNL::DNLID input) {
  if (iso.getIsoID() == naja::DNL::DNLID_MAX || iso.getDrivers().empty() ||
      iso.getDrivers().front() != input) {
    return false;
  }
  return Tree2BoolExpr::iso2boolExpr_.find(iso.getIsoID()) !=
         Tree2BoolExpr::iso2boolExpr_.end();
}

const char* kSkippedMultiDriverPOReport = "skipped_multi_driver_pos.txt";
const char* kSkippedNoDriverPOReport = "skipped_no_driver_pos.txt";
const char* kSkippedLogicalLoopPOReport = "skipped_logical_loop_pos.txt";

	void initializeSkippedPOReportFiles() {
	  static std::once_flag once;
	  if (!shouldReportSkippedPOs()) {
	    return; // LCOV_EXCL_LINE
	  }
  std::call_once(once, []() {
    std::ofstream(kSkippedMultiDriverPOReport, std::ios::trunc);
    std::ofstream(kSkippedNoDriverPOReport, std::ios::trunc);
    std::ofstream(kSkippedLogicalLoopPOReport, std::ios::trunc);
  });
}

DNLID getReportIsoID(const DNLFull* dnl, DNLID currentInput, DNLID mergeTerm) {
	  auto getIsoID = [&](DNLID termID) -> DNLID {
	    if (termID == DNLID_MAX) {
	      return DNLID_MAX; // LCOV_EXCL_LINE
	    }
    return dnl->getDNLTerminalFromID(termID).getIsoID();
  };

  const DNLID currentIsoID = getIsoID(currentInput);
  if (currentIsoID != DNLID_MAX) {
    return currentIsoID;
  }
  return getIsoID(mergeTerm); // LCOV_EXCL_LINE
}

struct SkippedPOReportEvent {
  const DNLFull* dnl = nullptr;
  DNLID rootTerm = DNLID_MAX;
  DNLID currentInput = DNLID_MAX;
  DNLID mergeTerm = DNLID_MAX;
  DNLID reportIsoID = DNLID_MAX;
  const char* reason = nullptr;
  const char* reportFile = nullptr;
  std::vector<DNLID, tbb::tbb_allocator<DNLID>> loopTerms;
};

using SkippedPOReportEventVector =
    std::vector<SkippedPOReportEvent, tbb::tbb_allocator<SkippedPOReportEvent>>;

thread_local SkippedPOReportEventVector skippedPOReportEvents;

SkippedPOReportEventVector& getSkippedPOReportEvents() {
  return skippedPOReportEvents;
}

using SkippedPOReportEventsETS =
    tbb::enumerable_thread_specific<SkippedPOReportEventVector*>;
SkippedPOReportEventsETS skippedPOReportEventsETS;

void recordSkippedPOReportEvent(const SkippedPOReportEvent& event) {
  auto& events = getSkippedPOReportEvents();
  if (events.empty()) {
    skippedPOReportEventsETS.local() = &events;
  }
  events.emplace_back(event);
}

	void appendCloudTermName(std::ostream& out, const DNLFull* dnl, DNLID termID) {
	  if (termID == DNLID_MAX) {
	    out << "<invalid>"; // LCOV_EXCL_LINE
	    return; // LCOV_EXCL_LINE
	  }
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.getDNLInstance().getSNLModel()) {
    out << term.getDNLInstance().getSNLModel()->getName().getString() << ":";
  }
  const auto path = term.getDNLInstance().getPath().getPathNames();
  for (size_t i = 0; i < path.size(); ++i) {
    out << path[i].getString() << ".";
  }
  out << term.getSnlBitTerm()->getName().getString()
      << term.getSnlBitTerm()->getBit()
      << " (term_id=" << termID
      << ", iso=" << term.getIsoID() << ")";
}

template <typename TermIDs>
void appendTermsToReport(std::ostream& out,
                         const DNLFull* dnl,
                         const char* label,
                         const TermIDs& termIDs) {
  out << label << ": [";
  for (size_t i = 0; i < termIDs.size(); ++i) {
    appendCloudTermName(out, dnl, termIDs[i]);
    if (i + 1 != termIDs.size()) {
      out << ", ";
    }
  }
  out << "]";
}

	void appendIsoDetailsToReport(std::ostream& out,
	                              const DNLFull* dnl,
	                              DNLID termID,
	                              const char* label) {
	  if (termID == DNLID_MAX) {
	    return; // LCOV_EXCL_LINE
	  }
	  const auto& term = dnl->getDNLTerminalFromID(termID);
	  if (term.getIsoID() == DNLID_MAX) {
	    return; // LCOV_EXCL_LINE
	  }
  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  out << " " << label << "_iso=" << term.getIsoID() << " ";
  appendTermsToReport(out, dnl, "readers", iso.getReaders());
  out << " ";
  appendTermsToReport(out, dnl, "drivers", iso.getDrivers());
}

void reportCloudSkippedRoot(const DNLFull* dnl,
                            DNLID rootTerm,
                            DNLID currentInput,
                            DNLID mergeTerm,
                            const char* reason,
                            const char* reportFile,
                            const std::vector<DNLID>* loopTerms = nullptr) {
  if (!shouldReportSkippedPOs()) {
    return;
  }
  SkippedPOReportEvent event;
  event.dnl = dnl;
  event.rootTerm = rootTerm;
  event.currentInput = currentInput;
  event.mergeTerm = mergeTerm;
  event.reportIsoID = getReportIsoID(dnl, currentInput, mergeTerm);
  event.reason = reason;
  event.reportFile = reportFile;
  if (loopTerms) {
    event.loopTerms.assign(loopTerms->begin(), loopTerms->end());
  }
  recordSkippedPOReportEvent(event);
}

}  // namespace

typedef std::pair<
    std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>,
    size_t>
    IterationInputsETSPair;

thread_local IterationInputsETSPair currentIterationInputsETS;

IterationInputsETSPair& getCurrentIterationInputsETS() {
  return currentIterationInputsETS;
}

thread_local IterationInputsETSPair newIterationInputsETS;

IterationInputsETSPair& getNewIterationInputsETS() {
  return newIterationInputsETS;
}

void clearCurrentIterationInputsETS() {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  currentIterationInputs.first.clear();
}

void pushBackCurrentIterationInputsETS(naja::DNL::DNLID input) {
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  currentIterationInputs.first.emplace_back(input);
}

size_t sizeOfCurrentIterationInputsETS() {
  return getCurrentIterationInputsETS().first.size();
}

void copyCurrentIterationInputsETS(std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& res) {
  res.clear();
  auto& current = getCurrentIterationInputsETS();
  res = std::move(current.first);
}

void clearNewIterationInputsETS() {
  auto& newIterationInputs = getNewIterationInputsETS();
  newIterationInputs.first.clear();
}

void pushBackNewIterationInputsETS(naja::DNL::DNLID input) {
  getNewIterationInputsETS().first.emplace_back(input);
}

bool emptyNewIterationInputsETS() {
  return getNewIterationInputsETS().first.empty();
}

size_t sizeOfNewIterationInputsETS() {
  return getNewIterationInputsETS().first.size();
}

void copyNewIterationInputsETStoCurrent() {
  auto& newIterationInputs = getNewIterationInputsETS();
  auto& currentIterationInputs = getCurrentIterationInputsETS();
  #ifdef DEBUG_CHECKS
  size_t newSize = newIterationInputs.first.size();
  #endif
  currentIterationInputs = std::move(newIterationInputs);
  #ifdef DEBUG_CHECKS
  assert(currentIterationInputs.first.size() == newSize &&
         "copyNewIterationInputsETStoCurrent: size mismatch after copy");
  #endif
}

thread_local std::pair<
    std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
                tbb::tbb_allocator<std::pair<naja::DNL::DNLID,
                                            naja::DNL::DNLID>>>,
    size_t>
    inputsToMergeETS;

struct PairHash {
  size_t operator()(const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& p) const noexcept {
    uint64_t a = static_cast<uint64_t>(p.first);
    uint64_t b = static_cast<uint64_t>(p.second);
    return (a * 11400714819323198485ull) ^ (b + 0x9e3779b97f4a7c15ull + (a<<6) + (a>>2));
    }
  };
  struct PairEq {
    bool operator()(const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& x, const std::pair<naja::DNL::DNLID,naja::DNL::DNLID>& y) const noexcept {
      return x.first == y.first && x.second == y.second;
    }
  };
using HandledSet = std::unordered_set<
    std::pair<naja::DNL::DNLID,naja::DNL::DNLID>,
    PairHash, PairEq,
    tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>;

std::pair<std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
                      tbb::tbb_allocator<std::pair<naja::DNL::DNLID,
                                                  naja::DNL::DNLID>>>, size_t>&
getInputsToMergeETS() {
  return inputsToMergeETS;
}

void clearInputsToMergeETS() {
  auto& inputsToMerge = getInputsToMergeETS();
  inputsToMerge.first.clear();
}

void pushBackInputsToMergeETS(
    const std::pair<naja::DNL::DNLID, naja::DNL::DNLID>& input) {
  getInputsToMergeETS().first.emplace_back(input);
}

size_t sizeOfInputsToMergeETS() {
  return getInputsToMergeETS().first.size();
}

// 2 level vector visited terms pair - 1st: termID, 2nd: termID
typedef std::vector<
    std::unordered_set<naja::DNL::DNLID, 
                                                 std::hash<naja::DNL::DNLID>,
                                                 std::equal_to<naja::DNL::DNLID>,
                                                 tbb::tbb_allocator<naja::DNL::DNLID>>,
    tbb::tbb_allocator<std::unordered_set<naja::DNL::DNLID, 
                                                 std::hash<naja::DNL::DNLID>,
                                                 std::equal_to<naja::DNL::DNLID>,
                                                 tbb::tbb_allocator<naja::DNL::DNLID>>>>
    VisitedTermsPairsVec;

thread_local VisitedTermsPairsVec visitedTermsPairsETS;
thread_local HandledSet visitedTermsPairsETSSet;

void clearVisitedTermsPairsETS() {
  visitedTermsPairsETSSet.clear();
}

thread_local std::pair<naja::DNL::DNLID, naja::DNL::DNLID> tempPairETS;

bool isPairVisitedETS(naja::DNL::DNLID termA,
                              naja::DNL::DNLID termB) {
  tempPairETS.first = termA;
  tempPairETS.second = termB;
  if (!(visitedTermsPairsETSSet.insert(tempPairETS)).second) {
    return true;
  }
  return false;
}

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return PIs_[termID];
}

bool SNLLogicCloud::isOutput(naja::DNL::DNLID termID) {
  return POs_[termID];
}

void SNLLogicCloud::compute() {
  clearNewIterationInputsETS();
  clearCurrentIterationInputsETS();
  struct TruthTableKey {
    const SNLDesign* design;
    size_t flatTermID;

    bool operator==(const TruthTableKey& other) const {
      return design == other.design && flatTermID == other.flatTermID;
    }
  };
  struct TruthTableKeyHash {
    size_t operator()(const TruthTableKey& key) const {
      return std::hash<const SNLDesign*>{}(key.design) ^
             (std::hash<size_t>{}(key.flatTermID) << 1);
    }
  };
  std::unordered_map<const SNLDesign*, size_t> truthTableCountCache;
  std::unordered_map<TruthTableKey, SNLTruthTable, TruthTableKeyHash>
      truthTableCache;
  auto getTruthTableCountCached = [&](const SNLDesign* model) {
    auto it = truthTableCountCache.find(model);
    if (it != truthTableCountCache.end()) {
      return it->second;
    }
    const size_t count = SNLDesignModeling::getTruthTableCount(model);
    truthTableCountCache.emplace(model, count);
    const auto& entry = truthTableCountCache.find(model);
    return count;
  };
  auto getTruthTableCached = [&](const SNLDesign* model, size_t flatTermID) {
    const TruthTableKey key{model, flatTermID};
    const auto& it = truthTableCache.find(key);
    if (it != truthTableCache.end()) {
      return it->second;
    }
    auto tt = SNLDesignModeling::getTruthTable(model, flatTermID);
    const auto& entry = truthTableCache.emplace(key, tt);
    return entry.first->second;
  };
  auto formatTermName = [&](naja::DNL::DNLID termID) {
    const auto& term = dnl_.getDNLTerminalFromID(termID);
    std::ostringstream fullName;
    fullName << "term_id=" << termID;
    if (term.getSnlBitTerm()) {
      fullName << ", term=" << term.getSnlBitTerm()->getName().getString()
               << ", flat_term_id=" << term.getSnlBitTerm()->getOrderID()
               << ", bit=" << term.getSnlBitTerm()->getBit();
    }
    if (term.getDNLInstance().getSNLModel()) {
      fullName << ", model="
               << term.getDNLInstance().getSNLModel()->getName().getString();
    }
    return fullName.str();
  };
  const bool captureFrontierHistory =
      std::getenv("KEPLER_CAPTURE_FRONTIER_HISTORY") != nullptr;
  std::vector<std::string> frontierHistory;
  // LCOV_EXCL_START
  auto appendTermList = [&](std::ostringstream& out,
                            const auto& termIDs,
                            size_t limit = 24) {
    const size_t capped = std::min(termIDs.size(), limit);
    for (size_t i = 0; i < capped; ++i) {
      out << "#" << i << "{" << formatTermName(termIDs[i]) << "}";
      if (i + 1 != capped) {
        out << ", ";
      }
    }
    if (termIDs.size() > capped) {
      out << ", ... +" << (termIDs.size() - capped) << " more";
    }
  };
  auto appendInstNonOutputs = [&](std::ostringstream& out,
                                  naja::DNL::DNLID instID,
                                  size_t limit = 24) {
    if (instID == naja::DNL::DNLID_MAX) {
      out << "<pi-or-constant>";
      return;
    }
    const auto& inst = dnl_.getDNLInstanceFromID(instID);
    size_t emitted = 0;
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; ++termID) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() ==
          SNLBitTerm::Direction::Output) {
        continue;
      }
      if (emitted != 0) {
        out << ", ";
      }
      out << formatTermName(termID);
      ++emitted;
      if (emitted == limit) {
        out << ", ...";
        return;
      }
    }
    if (emitted == 0) {
      out << "<none>";
    }
  };
  auto appendMergeListDetailed = [&](std::ostringstream& out,
                                     size_t limit = 24) {
    const auto& merges = getInputsToMergeETS().first;
    const size_t capped = std::min(merges.size(), limit);
    for (size_t i = 0; i < capped; ++i) {
      out << "#" << i << "{inst=" << merges[i].first
          << ", term=" << formatTermName(merges[i].second)
          << ", inst_non_outputs=[";
      appendInstNonOutputs(out, merges[i].first);
      out << "]}";
      if (i + 1 != capped) {
        out << ", ";
      }
    }
    if (merges.size() > capped) {
      out << ", ... +" << (merges.size() - capped) << " more";
    }
  };
  auto buildIterationSnapshot = [&](size_t iter) {
    std::ostringstream snapshot;
    snapshot << "iter " << iter << ": current_inputs=[";
    appendTermList(snapshot, getCurrentIterationInputsETS().first);
    snapshot << "] inputs_to_merge=[";
    appendMergeListDetailed(snapshot);
    snapshot << "] next_inputs=[";
    appendTermList(snapshot, getNewIterationInputsETS().first);
    snapshot << "]";
    return snapshot.str();
  };
  auto collectRelevantInstanceInputs = [&](naja::DNL::DNLID driver) {
    const auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    const auto* model = inst.getSNLModel();

    std::vector<naja::DNL::DNLID> allNonOutputTerms;
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; ++termID) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() ==
          SNLBitTerm::Direction::Output) {
        continue;
      }
      allNonOutputTerms.push_back(termID);
    }

    if (!NLDB0::isMux2(model)) {
      return allNonOutputTerms;
    }

    std::unordered_map<size_t, naja::DNL::DNLID> orderIDToTermID;
    orderIDToTermID.reserve(allNonOutputTerms.size());
    for (const auto termID : allNonOutputTerms) {
      orderIDToTermID.emplace(
          dnl_.getDNLTerminalFromID(termID).getSnlBitTerm()->getOrderID(),
          termID);
    }

    const auto* driverTerm = dnl_.getDNLTerminalFromID(driver).getSnlBitTerm();
    const auto driverBit = static_cast<NLID::Bit>(driverTerm->getBit());
    const auto* inputA = NLDB0::getMux2InputA(model);
    const auto* inputB = NLDB0::getMux2InputB(model);
    const auto* select = NLDB0::getMux2Select(model);
    if (!inputA || !inputB || !select) {
      throw std::runtime_error("SNLLogicCloud failed to resolve wide mux inputs");
    }
    std::vector<naja::DNL::DNLID> relevantTerms;
    relevantTerms.reserve(3);
    for (const auto* bitTerm :
         {inputA->getBit(driverBit), inputB->getBit(driverBit)}) {
      if (!bitTerm) {
        std::ostringstream error;
        error << "SNLLogicCloud failed to resolve wide mux bit "
              << static_cast<size_t>(driverBit) << " for driver "
              << formatTermName(driver) << " in model "
              << model->getName().getString();
        throw std::runtime_error(error.str());
      }
      const auto it = orderIDToTermID.find(bitTerm->getOrderID());
      if (it == orderIDToTermID.end()) {
        std::ostringstream error;
        error << "SNLLogicCloud failed to map wide mux dependency "
              << bitTerm->getOrderID() << " for driver "
              << formatTermName(driver)
              << " in model " << model->getName().getString();
        throw std::runtime_error(error.str());
      }
      relevantTerms.push_back(it->second);
    }
    const auto selectIt = orderIDToTermID.find(select->getOrderID());
    if (selectIt == orderIDToTermID.end()) {
      std::ostringstream error;
      error << "SNLLogicCloud failed to map wide mux select dependency "
            << select->getOrderID() << " for driver "
            << formatTermName(driver) << " in model "
            << model->getName().getString();
      throw std::runtime_error(error.str());
    }
    relevantTerms.push_back(selectIt->second);
    return relevantTerms;
  };
  // LCOV_EXCL_STOP
	  auto throwIfTruthTableArityMismatch = [&](naja::DNL::DNLID driver) {
	    const auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
	    const auto* model = inst.getSNLModel();
	    if (getTruthTableCountCached(model) == 0) {
	      return; // LCOV_EXCL_LINE
	    }
	    const auto& tt = getTruthTableCached(
	        model, dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID());
	    if (!tt.isInitialized()) {
	      return; // LCOV_EXCL_LINE
	    }
    const auto relevantInputs = collectRelevantInstanceInputs(driver);
    const size_t expectedInputCount =
        NLDB0::isMux2(model) ? size_t{3} : tt.size();
    if (expectedInputCount == relevantInputs.size()) {
      return;
    }

    size_t modelNonOutputCount = 0;
    std::ostringstream modelNonOutputTerms;
    bool firstModelTerm = true;
    for (const auto* term : model->getBitTerms()) {
      if (term->getDirection() != SNLBitTerm::Direction::Output) {
        ++modelNonOutputCount;
        if (!firstModelTerm) {
          modelNonOutputTerms << ", ";
        }
        firstModelTerm = false;
        modelNonOutputTerms << "flat_term_id=" << term->getOrderID()
                            << ", bit=" << term->getBit();
      }
    }

    std::ostringstream instanceNonOutputTerms;
    bool firstInstanceTerm = true;
    for (const auto termID : relevantInputs) {
      if (!firstInstanceTerm) {
        instanceNonOutputTerms << ", ";
      }
      firstInstanceTerm = false;
      instanceNonOutputTerms << formatTermName(termID);
    }

    std::ostringstream error;
    error << "SNLLogicCloud arity mismatch for model "
          << model->getName().getString() << ": TT arity=" << tt.size()
          << ", model non-output term count=" << modelNonOutputCount
          << ", instance non-output term count=" << relevantInputs.size()
          << ", driver=" << formatTermName(driver)
          << ", model non-output terms=[" << modelNonOutputTerms.str()
          << "], instance non-output terms=["
          << instanceNonOutputTerms.str() << "], tt=" << tt.toString();
    throw std::runtime_error(error.str());
  };
  // LCOV_EXCL_START
  auto throwIfFrontierMismatch = [&](size_t iter) {
    const size_t currentCount = sizeOfCurrentIterationInputsETS();
    const size_t mergeCount = sizeOfInputsToMergeETS();
    const size_t borderCount = table_.getBorderLeavesSize();
    if (currentCount == borderCount && mergeCount == borderCount) {
      return;
    }

    constexpr size_t kMaxEntries = 24;
    auto appendInputList = [&](std::ostringstream& error) {
      error << " current_inputs=[";
      const auto& current = getCurrentIterationInputsETS().first;
      const size_t limit = std::min(current.size(), kMaxEntries);
      for (size_t i = 0; i < limit; ++i) {
        const auto input = current[i];
        const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
            dnl_.getDNLTerminalFromID(input).getIsoID());
        error << "#" << i << "{"
              << formatTermName(input)
              << ", iso=" << iso.getIsoID()
              << ", drivers=" << iso.getDrivers().size()
              << ", readers=" << iso.getReaders().size()
              << "}";
        if (i + 1 != limit) {
          error << ", ";
        }
      }
      if (current.size() > limit) {
        error << ", ... +" << (current.size() - limit) << " more";
      }
      error << "]";
    };
    auto appendMergeList = [&](std::ostringstream& error) {
      error << " inputs_to_merge=[";
      const auto& merges = getInputsToMergeETS().first;
      const size_t limit = std::min(merges.size(), kMaxEntries);
      for (size_t i = 0; i < limit; ++i) {
        error << "#" << i << "{inst=" << merges[i].first
              << ", term=" << formatTermName(merges[i].second) << "}";
        if (i + 1 != limit) {
          error << ", ";
        }
      }
      if (merges.size() > limit) {
        error << ", ... +" << (merges.size() - limit) << " more";
      }
      error << "]";
    };
    auto appendDuplicateMergeTerms = [&](std::ostringstream& error) {
      std::map<naja::DNL::DNLID, size_t> counts;
      for (const auto& merge : getInputsToMergeETS().first) {
        ++counts[merge.second];
      }
      bool emitted = false;
      error << " duplicate_merge_terms=[";
      size_t shown = 0;
      for (const auto& [termID, count] : counts) {
        if (count <= 1) {
          continue;
        }
        if (shown == kMaxEntries) {
          error << "...";
          emitted = true;
          break;
        }
        if (shown != 0) {
          error << ", ";
        }
        error << formatTermName(termID) << " x" << count;
        emitted = true;
        ++shown;
      }
      if (!emitted) {
        error << "<none>";
      }
      error << "]";
    };
    auto appendDuplicateCurrentInputs = [&](std::ostringstream& error) {
      std::map<naja::DNL::DNLID, size_t> counts;
      for (const auto input : getCurrentIterationInputsETS().first) {
        ++counts[input];
      }
      bool emitted = false;
      error << " duplicate_current_inputs=[";
      size_t shown = 0;
      for (const auto& [termID, count] : counts) {
        if (count <= 1) {
          continue;
        }
        if (shown == kMaxEntries) {
          error << "...";
          emitted = true;
          break;
        }
        if (shown != 0) {
          error << ", ";
        }
        error << formatTermName(termID) << " x" << count;
        emitted = true;
        ++shown;
      }
      if (!emitted) {
        error << "<none>";
      }
      error << "]";
    };

    std::ostringstream error;
    error << "SNLLogicCloud frontier mismatch before concat at iter " << iter
          << ": current_inputs=" << currentCount
          << ", inputs_to_merge=" << mergeCount
          << ", border_leaves=" << borderCount;
    appendInputList(error);
    appendMergeList(error);
    appendDuplicateCurrentInputs(error);
    appendDuplicateMergeTerms(error);
    if (!frontierHistory.empty()) {
      error << " history=[";
      for (size_t i = 0; i < frontierHistory.size(); ++i) {
        if (i != 0) {
          error << " || ";
        }
        error << frontierHistory[i];
      }
      error << "]";
    }
    throw std::runtime_error(error.str());
  };
  auto throwIfNextFrontierMismatch = [&](size_t iter) {
    const size_t nextCount = sizeOfNewIterationInputsETS();
    const size_t borderCount = table_.getBorderLeavesSize();
    if (nextCount == borderCount) {
      return;
    }
    constexpr size_t kMaxEntries = 24;
    std::ostringstream error;
    error << "SNLLogicCloud next frontier mismatch after concat at iter "
          << iter << ": next_inputs=" << nextCount
          << ", border_leaves=" << borderCount;

    error << " next_inputs=[";
    const auto& nextInputs = getNewIterationInputsETS().first;
    const size_t nextLimit = std::min(nextInputs.size(), kMaxEntries);
    for (size_t i = 0; i < nextLimit; ++i) {
      const auto input = nextInputs[i];
      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      error << "#" << i << "{"
            << formatTermName(input)
            << ", iso=" << iso.getIsoID()
            << ", drivers=" << iso.getDrivers().size()
            << ", readers=" << iso.getReaders().size()
            << "}";
      if (i + 1 != nextLimit) {
        error << ", ";
      }
    }
    if (nextInputs.size() > nextLimit) {
      error << ", ... +" << (nextInputs.size() - nextLimit) << " more";
    }
    error << "]";

    error << " from_current_inputs=[";
    const auto& current = getCurrentIterationInputsETS().first;
    const size_t currentLimit = std::min(current.size(), kMaxEntries);
    for (size_t i = 0; i < currentLimit; ++i) {
      error << "#" << i << "{" << formatTermName(current[i]) << "}";
      if (i + 1 != currentLimit) {
        error << ", ";
      }
    }
    if (current.size() > currentLimit) {
      error << ", ... +" << (current.size() - currentLimit) << " more";
    }
    error << "]";

    error << " from_inputs_to_merge=[";
    const auto& merges = getInputsToMergeETS().first;
    const size_t mergeLimit = std::min(merges.size(), kMaxEntries);
    for (size_t i = 0; i < mergeLimit; ++i) {
      error << "#" << i << "{inst=" << merges[i].first
            << ", term=" << formatTermName(merges[i].second) << "}";
      if (i + 1 != mergeLimit) {
        error << ", ";
      }
    }
    if (merges.size() > mergeLimit) {
      error << ", ... +" << (merges.size() - mergeLimit) << " more";
    }
    error << "]";
    error << " failing_iteration_detail=[" << buildIterationSnapshot(iter)
          << "]";
    if (!frontierHistory.empty()) {
      error << " history=[";
      for (size_t i = 0; i < frontierHistory.size(); ++i) {
        if (i != 0) {
          error << " || ";
        }
        error << frontierHistory[i];
      }
      error << "]";
    }
    throw std::runtime_error(error.str());
  };
  // LCOV_EXCL_STOP
  DEBUG_LOG("---- Begin!!\n");
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() ||
      isOutput(seedOutputTerm_)) {
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(seedOutputTerm_).getIsoID());
    // LCOV_EXCL_START
    if (iso.getDrivers().size() > 1) {
      #ifdef DEBUG_PRINTS
      for (const auto& driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n", dnl_.getDNLTerminalFromID(driver)
                                      .getSnlBitTerm()
                                      ->getName()
                                      .getString()
                                      .c_str());
      }
      #endif
      throw std::runtime_error("Seed output term is not a single driver");
    } else if (iso.getDrivers().empty()) {
      std::string termName = dnl_.getDNLTerminalFromID(seedOutputTerm_)
                                 .getSnlBitTerm()
                                 ->getName()
                                 .getString();
      std::string error =
          "Seed output term '" + termName + "' has no drivers";
      throw std::runtime_error(error);
    }
    // LCOV_EXCL_STOP
    const auto& driver = iso.getDrivers().front();
    auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    if (isInput(driver)) {
      pushBackCurrentIterationInputsETS(driver);
      table_ = SNLTruthTableTree(inst.getID(), driver,
                                 SNLTruthTableTree::Node::Type::P);
      return;
    }
    throwIfTruthTableArityMismatch(driver);
    DEBUG_LOG("Instance name: %s\n",
              inst.getSNLInstance()->getName().getString().c_str());
    for (const auto termID : collectRelevantInstanceInputs(driver)) {
      pushBackNewIterationInputsETS(termID);
      DEBUG_LOG("Add input with id: %zu\n", termID);
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), driver);
    auto* model = inst.getSNLModel();
    assert(SNLDesignModeling::getTruthTable(model, 
                dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID())
            .isInitialized() &&
        "Truth table is not initialized");
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    const auto& inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);
    for (DNLID termID = inst.getTermIndexes().first;
         termID <= inst.getTermIndexes().second; termID++) {
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        // newIterationInputs.emplace_back(termID);
        pushBackNewIterationInputsETS(termID);
        DEBUG_LOG("Add input with id: %zu\n", termID);
      }
    }
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), seedOutputTerm_);
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  }

  if (emptyNewIterationInputsETS()) {
    DEBUG_LOG("No inputs found for seed output term %zu\n", seedOutputTerm_);
    return;
  }
	  // LCOV_EXCL_START
	  if (captureFrontierHistory) {
	    std::ostringstream seedInfo;
	    seedInfo << "seed_output={" << formatTermName(seedOutputTerm_) << "}"
	             << " initial_inputs=[";
	    appendTermList(seedInfo, getNewIterationInputsETS().first);
	    seedInfo << "]";
	    frontierHistory.emplace_back(seedInfo.str());
	  }
	  // LCOV_EXCL_STOP

  bool reachedPIs = true;
  size_t size = sizeOfNewIterationInputsETS();
  
  for (size_t i = 0; i < size; i++) {
    auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
      dnl_.getDNLTerminalFromID(getNewIterationInputsETS().first[i]).getIsoID());
    if (!isInput(
            getNewIterationInputsETS().first
                [i]) /* && !isOutput(getNewIterationInputsETS().first[i])*/) {
      reachedPIs = false;
      break;
    }
  } // allocator for buckets

  // HandledSet handledTerms;
  // handledTerms.reserve(naja::DNL::get()->getDNLTerms().size() / 4);
  clearVisitedTermsPairsETS();
  size_t iter = 0;

  while (!reachedPIs) {
    // Originally computation of reachedPIs have been handled in the end of the loop,
    // but by adding isConstant on isos as part of the check, it had to be moved to the beginning.
    // Why? because before we cached the inputs of the leaves we meet and then we look at the drivers in the next loop iteration
    // and then we run the checks on the drivers in order to know if we need to iterate again, after the drivers were contacted to the cloud.
    // Now, we also check isos of inputs, which means that if an input has an iso that is constant, we will stop.
    // In this case, we have to make the check in the next loop iteration in order to force concating the constants to the cloud.
    reachedPIs = true;
    size_t sizeOfNewInputs = sizeOfNewIterationInputsETS();
    for (size_t i = 0; i < sizeOfNewInputs; i++) {
      const auto input = getNewIterationInputsETS().first[i];
      auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      const bool cachedAsInput = canUseCachedIsoShortcut(iso, input);
      if (!isInput(input) && !cachedAsInput && !iso.isConstant()) {
        reachedPIs = false;
        break;
      }
    }
    DEBUG_LOG("---iter %lu---\n", iter);
    DEBUG_LOG("Current iteration inputs size: %zu\n",
              sizeOfNewIterationInputsETS());
    copyNewIterationInputsETStoCurrent();

    clearNewIterationInputsETS();
    DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n",
              table_.size(), sizeOfCurrentIterationInputsETS());
    clearInputsToMergeETS();
    size_t sizeOfCurrentInputs = sizeOfCurrentIterationInputsETS();
    for (size_t i = 0; i < sizeOfCurrentInputs; i++) {
      const auto& input = getCurrentIterationInputsETS().first[i];
      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          dnl_.getDNLTerminalFromID(input).getIsoID());
      if (isInput(input) || iso.isConstant()) {
        pushBackNewIterationInputsETS(input);
        DEBUG_LOG("Adding input id: %zu %s\n", input,
                  dnl_.getDNLTerminalFromID(input)
                      .getSnlBitTerm()
                      ->getName()
                      .getString()
                      .c_str());
        pushBackInputsToMergeETS(
            {naja::DNL::DNLID_MAX, input});  // Placeholder for PI/PO
        continue;
      }

      
      DEBUG_LOG("number of drivers: %zu\n", iso.getDrivers().size());

      for (const auto& driver : iso.getDrivers()) {
        DEBUG_LOG("Driver: %s\n", dnl_.getDNLTerminalFromID(driver)
                                      .getSnlBitTerm()
                                      ->getName()
                                      .getString()
                                      .c_str());
      }

      auto appendTerms = [&](std::string& error,
                             const char* label,
                             const auto& termIDs) {
        error += label;
        error += ": [";
        for (size_t i = 0; i < termIDs.size(); ++i) {
          error += formatTermName(termIDs[i]);
          if (i + 1 != termIDs.size()) {
            error += ", ";
          }
        }
        error += "]";
      };
      auto formatNet = [](const naja::NL::SNLBitNet* net) {
        std::string description = "design=";
        description += net->getDesign()->getName().getString();
        description += " name=";
        description += net->getName().getString();
        description += " type=";
        description += net->getType().getString();
        description += " is_assign=";
        description += net->getType().isAssign() ? "true" : "false";
        description += " is_supply=";
        description += net->getType().isSupply() ? "true" : "false";
        description += " is_constant0=";
        description += net->isConstant0() ? "true" : "false";
        description += " is_constant1=";
        description += net->isConstant1() ? "true" : "false";
        description += " model_is_assign=";
        description += net->getDesign()->isAssign() ? "true" : "false";
        description += " properties=[";
        bool first = true;
        for (auto* property : net->getProperties()) {
          if (!first) {
            description += ", ";
          }
          first = false;
          description += property->getName();
          description += "=";
          description += property->getString();
        }
        description += "]";
        return description;
      };
      auto appendNets = [&](std::string& error,
                            const char* label,
                            const std::set<naja::NL::SNLBitNet*>& nets) {
        error += label;
        error += ": [";
        size_t i = 0;
        for (const auto* net : nets) {
          error += formatNet(net);
          if (++i != nets.size()) {
            error += ", ";
          }
        }
        error += "]";
      };
      auto appendComplexIsoDetails = [&](std::string& error) {
        naja::DNL::DNLComplexIso complexIso(iso.getIsoID());
        dnl_.getCustomIso(iso.getIsoID(), complexIso);
        error += " ";
        appendTerms(error, "complex_readers", complexIso.getReaders());
        error += " ";
        appendTerms(error, "complex_drivers", complexIso.getDrivers());
        error += " ";
        appendTerms(error, "complex_hier_terms", complexIso.getHierTerms());
        error += " ";
        appendNets(error, "complex_nets", complexIso.getNets());
      };

      if (iso.getDrivers().size() >= 1) {
        // proper error with names of all the drivers
        // throw an error and separate names by comma
        if (iso.getDrivers().size() > 1) {
          reportCloudSkippedRoot(&dnl_, seedOutputTerm_, input, DNLID_MAX,
                                 "its iso has multiple drivers during cloud expansion",
                                 kSkippedMultiDriverPOReport);
          table_ = SNLTruthTableTree();
          return;
        }
      } else if (iso.getDrivers().empty()) {
        if (!iso.isConstant()) {
          reportCloudSkippedRoot(&dnl_, seedOutputTerm_, input, DNLID_MAX,
                                 "its iso has no drivers during cloud expansion",
                                 kSkippedNoDriverPOReport);
          table_ = SNLTruthTableTree();
          return;
        }
	        // LCOV_EXCL_START
	        pushBackNewIterationInputsETS(input);
	        pushBackInputsToMergeETS(
	            {naja::DNL::DNLID_MAX, input});  // Placeholder for PI/PO
	        continue;
	        // LCOV_EXCL_STOP
	      }
      const auto& driver = iso.getDrivers().front();
      
      if (isInput(driver) || (canUseCachedIsoShortcut(iso, driver) && iter > 0)) {
        pushBackNewIterationInputsETS(driver);
        DEBUG_LOG(
            "- %lu After analyzing input %s(%lu), addings driver %s(%lu) is a "
            "primary input\n",
            iter,
            dnl_.getDNLTerminalFromID(input)
                .getSnlBitTerm()
                ->getName()
                .getString()
                .c_str(),
            input,
            dnl_.getDNLTerminalFromID(driver)
                .getSnlBitTerm()
                ->getName()
                .getString()
                .c_str(),
            driver);
        pushBackInputsToMergeETS(
            {naja::DNL::DNLID_MAX, driver});  // Placeholder for PI/PO
        continue;
      }

      const auto& inst = dnl_.getDNLInstanceFromID(
          dnl_.getDNLTerminalFromID(driver).getDNLInstance().getID());
      throwIfTruthTableArityMismatch(driver);

      DEBUG_LOG("Adding driver id: %zu %s(%s)\n", driver,
                dnl_.getDNLTerminalFromID(driver)
                    .getSnlBitTerm()
                    ->getName()
                    .getString()
                    .c_str(),
                dnl_.getDNLTerminalFromID(driver)
                    .getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      pushBackInputsToMergeETS({inst.getID(), driver});

      for (const auto termID : collectRelevantInstanceInputs(driver)) {
        if (isPairVisitedETS(driver, termID)) {
          DEBUG_LOG(
              "#### iter %lu 1 Term (%zu) %s of inst %s already handled, "
              "skipping\n",
              iter, input,
              naja::DNL::get()
                  ->getDNLTerminalFromID(input)
                  .getSnlBitTerm()
                  ->getName()
                  .getString()
                  .c_str(),
              naja::DNL::get()
                  ->getDNLTerminalFromID(input)
                  .getDNLInstance()
                  .getSNLModel()
                  ->getName()
                  .getString()
                  .c_str());
          continue;
        }
        pushBackNewIterationInputsETS(termID);
      }
    }

    if (sizeOfInputsToMergeETS() == 0) {
      break;
    }

    DEBUG_LOG("--- Merging truth tables with %zu inputs\n",
              sizeOfInputsToMergeETS());
    throwIfFrontierMismatch(iter);
    {
      const auto& merges = getInputsToMergeETS().first;
      const auto& currentInputs = getCurrentIterationInputsETS().first;
      for (size_t i = 0; i < merges.size(); ++i) {
        if (merges[i].first == naja::DNL::DNLID_MAX) {
          continue;
        }
        std::vector<DNLID> loopTerms;
        if (table_.findAncestorLoopForBorderLeaf(i, merges[i].second,
                                                 loopTerms)) {
          reportCloudSkippedRoot(
              &dnl_, seedOutputTerm_, currentInputs[i], merges[i].second,
              "a logical loop was detected during cloud expansion",
              kSkippedLogicalLoopPOReport, &loopTerms);
          table_ = SNLTruthTableTree();
          return;
        }
      }
    }
    table_.concatFull(getInputsToMergeETS().first,
                      sizeOfInputsToMergeETS());
    throwIfNextFrontierMismatch(iter);
	    // LCOV_EXCL_START
	    if (captureFrontierHistory) {
	      frontierHistory.emplace_back(buildIterationSnapshot(iter));
	    }
	    // LCOV_EXCL_STOP
    DEBUG_LOG("--- End of iteration %zu\n", iter);
    iter++;
  }

  copyNewIterationInputsETStoCurrent();
  #ifdef DEBUG_CHECKS
  size_t finalSize = sizeOfCurrentIterationInputsETS();
  #endif
  copyCurrentIterationInputsETS(currentIterationInputs_);
  #ifdef DEBUG_CHECKS
  assert(finalSize == currentIterationInputs_.size() &&
         "compute: size mismatch after final copy");
  //assert(currentIterationInputs_.size() == sizeOfCurrentIterationInputsETS());
  for (const auto& input : currentIterationInputs_) {
    auto iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(input).getIsoID());
    const bool cachedAsInput = canUseCachedIsoShortcut(iso, input);
    assert(isInput(input) || cachedAsInput || iso.isConstant());
  }
  #endif
}

void SNLLogicCloud::flushSkippedPOReports() {
  if (!shouldReportSkippedPOs()) {
    return;
  }

  initializeSkippedPOReportFiles();
  std::ofstream multiDriverOut(kSkippedMultiDriverPOReport, std::ios::app);
  std::ofstream noDriverOut(kSkippedNoDriverPOReport, std::ios::app);
  std::ofstream logicalLoopOut(kSkippedLogicalLoopPOReport, std::ios::app);
  std::unordered_set<uint64_t> seenIsoKeys;

  auto makeIsoKey = [](const DNLFull* dnl, DNLID isoID) -> uint64_t {
    return (static_cast<uint64_t>(reinterpret_cast<uintptr_t>(dnl)) << 32) ^
           static_cast<uint64_t>(isoID);
  };

	  for (auto* eventsPtr : skippedPOReportEventsETS) {
	    if (eventsPtr == nullptr) {
	      continue; // LCOV_EXCL_LINE
	    }
    auto& events = *eventsPtr;
    for (const auto& event : events) {
      std::ostream* out = &logicalLoopOut;
      const std::string_view reportFile(event.reportFile);
	      if (reportFile == std::string_view(kSkippedMultiDriverPOReport)) {
	        out = &multiDriverOut; // LCOV_EXCL_LINE
	      } else if (reportFile == std::string_view(kSkippedNoDriverPOReport)) {
	        out = &noDriverOut;
	      }

	      if (!(*out)) {
	        continue; // LCOV_EXCL_LINE
	      }

      const bool isFirstForIso =
          event.reportIsoID == DNLID_MAX ||
          seenIsoKeys.insert(makeIsoKey(event.dnl, event.reportIsoID)).second;

      *out << "Skipping cloud root ";
      appendCloudTermName(*out, event.dnl, event.rootTerm);
      *out << " because " << event.reason;
      if (event.reportIsoID != DNLID_MAX) {
        *out << ". iso=" << event.reportIsoID;
      }
	      if (!isFirstForIso) {
	        // LCOV_EXCL_START
	        if (event.reportIsoID != DNLID_MAX) {
	          *out << ". See first encounter of iso=" << event.reportIsoID
	               << " for details";
	        }
	        *out << "\n\n";
	        continue;
	        // LCOV_EXCL_STOP
	      }

      *out << ". current_input=";
      appendCloudTermName(*out, event.dnl, event.currentInput);
      if (event.mergeTerm != DNLID_MAX) {
        *out << " merge_term=";
        appendCloudTermName(*out, event.dnl, event.mergeTerm);
      }
      appendIsoDetailsToReport(*out, event.dnl, event.currentInput,
                               "current_input");
      if (event.mergeTerm != DNLID_MAX) {
        appendIsoDetailsToReport(*out, event.dnl, event.mergeTerm,
                                 "merge_term");
      }
      if (!event.loopTerms.empty()) {
        *out << " ";
        appendTermsToReport(*out, event.dnl, "loop_terms", event.loopTerms);
      }
      *out << "\n\n";
    }
    events.clear();
  }
}
