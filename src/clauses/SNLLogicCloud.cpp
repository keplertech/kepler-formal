// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "SNLLogicCloud.h"
#include <tbb/tbb_allocator.h>
#include <cassert>
#include <sstream>
#include "NajaProperty.h"
#include "NLBitDependencies.h"
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
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;

namespace {

constexpr size_t InvalidFlatTermID = std::numeric_limits<size_t>::max();

struct ModelInputLayoutKey {
  const SNLDesign* model = nullptr;

  bool operator==(const ModelInputLayoutKey& other) const {
    return model == other.model;
  }
};

struct ModelInputLayoutKeyHash {
  size_t operator()(const ModelInputLayoutKey& key) const {
    return std::hash<const SNLDesign*>{}(key.model);
  }
};

struct ModelInputLayout {
  bool isMux2 = false;
  bool isTableSelect = false;
  bool isAssign = false;
  size_t bitTermCount = 0;
  std::vector<size_t> nonOutputTermFlatIDs;
  std::unordered_map<NLID::Bit, size_t> muxInputAFlatIDs;
  std::unordered_map<NLID::Bit, size_t> muxInputBFlatIDs;
  size_t muxSelectFlatID = InvalidFlatTermID;
};

struct TruthTableKey {
  const SNLDesign* design = nullptr;
  const naja::NL::SNLInstance* instance = nullptr;
  size_t flatTermID = 0;

  bool operator==(const TruthTableKey& other) const {
    return design == other.design && instance == other.instance &&
           flatTermID == other.flatTermID;
  }
};

struct TruthTableKeyHash {
  size_t operator()(const TruthTableKey& key) const {
    return std::hash<const SNLDesign*>{}(key.design) ^
           (std::hash<const naja::NL::SNLInstance*>{}(key.instance) << 1) ^
           (std::hash<size_t>{}(key.flatTermID) << 2);
  }
};

thread_local std::unordered_map<ModelInputLayoutKey,
                                ModelInputLayout,
                                ModelInputLayoutKeyHash>
    modelInputLayoutCache;
thread_local std::unordered_map<TruthTableKey, size_t, TruthTableKeyHash>
    truthTableCountCache;
thread_local std::unordered_map<TruthTableKey, SNLTruthTable, TruthTableKeyHash>
    truthTableCache;
thread_local std::vector<uint32_t> expandedTableTermEpochs;
thread_local uint32_t expandedTableTermEpoch = 1;
thread_local std::unordered_map<naja::DNL::DNLID, naja::DNL::DNLID>
    transparentLoopTargetCacheTL;
thread_local std::unordered_set<naja::DNL::DNLID>
    transparentLoopVisitedTermsTL;
thread_local std::vector<naja::DNL::DNLID> loopTermsScratchTL;

ModelInputLayoutKey makeModelInputLayoutKey(const SNLDesign* model) {
  return ModelInputLayoutKey{model};
}

uint64_t mixSignature(uint64_t seed, uint64_t value) {
  return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

uint64_t pointerSignature(const void* ptr) {
  return static_cast<uint64_t>(reinterpret_cast<uintptr_t>(ptr));
}

uint64_t modelSignature(const SNLDesign* model) {
  uint64_t signature = pointerSignature(model);
  if (model == nullptr) {
    return signature;
  }
  const auto id = model->getNLID();
  signature = mixSignature(signature, model->getName().getID());
  signature = mixSignature(signature, static_cast<uint64_t>(id.dbID_));
  signature = mixSignature(signature, static_cast<uint64_t>(id.libraryID_));
  signature = mixSignature(signature, static_cast<uint64_t>(id.designID_));
  signature = mixSignature(signature, static_cast<uint64_t>(model->getID()));
  return signature;
}

uint64_t instanceSignature(const DNLInstanceFull& inst) {
  uint64_t signature = mixSignature(inst.getID(), inst.getParentID());
  signature = mixSignature(signature, inst.getTermIndexes().first);
  signature = mixSignature(signature, inst.getTermIndexes().second);
  signature = mixSignature(signature, inst.getChildren().first);
  signature = mixSignature(signature, inst.getChildren().second);
  signature = mixSignature(signature, pointerSignature(inst.getSNLInstance()));
  signature = mixSignature(signature, modelSignature(inst.getSNLModel()));
  return signature;
}

uint64_t termSignature(const DNLTerminalFull& term) {
  uint64_t signature = mixSignature(term.getID(), term.getIsoID());
  signature = mixSignature(signature,
                           pointerSignature(term.getSnlBitTerm()));
  signature = mixSignature(signature, pointerSignature(term.getSnlTerm()));
  signature = mixSignature(signature, instanceSignature(term.getDNLInstance()));
  if (const auto* bitTerm = term.getSnlBitTerm()) {
    signature = mixSignature(signature, bitTerm->getOrderID());
    signature = mixSignature(signature,
        static_cast<uint64_t>(bitTerm->getDirection()));
    signature = mixSignature(signature,
        static_cast<uint64_t>(bitTerm->getBit()));
    signature = mixSignature(signature, bitTerm->getName().getID());
    signature = mixSignature(signature, modelSignature(bitTerm->getDesign()));
  }
  return signature;
}

uint64_t dnlContextSignature(const DNLFull& dnl) {
  const auto& instances = dnl.getDNLInstances();
  const auto& terms = dnl.getDNLTerms();
  uint64_t signature = pointerSignature(&dnl);
  signature = mixSignature(signature, pointerSignature(instances.data()));
  signature = mixSignature(signature, instances.size());
  signature = mixSignature(signature, instances.capacity());
  signature = mixSignature(signature, pointerSignature(terms.data()));
  signature = mixSignature(signature, terms.size());
  signature = mixSignature(signature, terms.capacity());
  signature = mixSignature(signature, dnl.getNBterms());
  if (!instances.empty()) {
    signature = mixSignature(signature, instanceSignature(instances.front()));
    signature = mixSignature(signature, instanceSignature(instances.back()));
    signature =
        mixSignature(signature, instanceSignature(instances[instances.size() / 2]));
  }
  if (dnl.getNBterms() > 0) {
    signature =
        mixSignature(signature, termSignature(dnl.getDNLTerminalFromID(0)));
    const auto lastTerm = dnl.getNBterms() - 1;
    signature = mixSignature(
        signature, termSignature(dnl.getDNLTerminalFromID(lastTerm)));
    signature = mixSignature(
        signature, termSignature(dnl.getDNLTerminalFromID(lastTerm / 2)));
  }
  return signature;
}

void refreshPerDnlCaches(const DNLFull& dnl) {
  thread_local uint64_t lastDnlContextSignature = 0;
  const uint64_t currentSignature = dnlContextSignature(dnl);
  if (currentSignature == lastDnlContextSignature) {
    return;
  }
  modelInputLayoutCache.clear();
  truthTableCountCache.clear();
  truthTableCache.clear();
  expandedTableTermEpochs.clear();
  expandedTableTermEpoch = 1;
  transparentLoopTargetCacheTL.clear();
  transparentLoopVisitedTermsTL.clear();
  loopTermsScratchTL.clear();
  lastDnlContextSignature = currentSignature;
}

std::shared_ptr<const std::vector<DNLID>> getSharedTermIsoIDCache(
    const DNLFull& dnl) {
  static std::mutex mutex;
  static uint64_t cachedSignature = 0;
  static std::shared_ptr<const std::vector<DNLID>> cachedIsoIDs;

  const uint64_t currentSignature = dnlContextSignature(dnl);
  std::lock_guard<std::mutex> lock(mutex);
  if (cachedIsoIDs && cachedSignature == currentSignature) {
    return cachedIsoIDs;
  }

  auto isoIDs = std::make_shared<std::vector<DNLID>>();
  const size_t termCount = dnl.getNBterms();
  isoIDs->reserve(termCount);
  for (size_t i = 0; i < termCount; ++i) {
    isoIDs->push_back(
        dnl.getDNLTerminalFromID(static_cast<DNLID>(i)).getIsoID());
  }
  cachedSignature = currentSignature;
  cachedIsoIDs = std::move(isoIDs);
  return cachedIsoIDs;
}

const ModelInputLayout& getModelInputLayout(const DNLFull& dnl,
                                            const SNLDesign* model) {
  (void)dnl;
  const ModelInputLayoutKey key = makeModelInputLayoutKey(model);
  auto it = modelInputLayoutCache.find(key);
  if (it != modelInputLayoutCache.end()) {
    return it->second;
  }

  ModelInputLayout layout;
  if (model != nullptr) {
    layout.isMux2 = NLDB0::isMux2(model);
    layout.isTableSelect = NLDB0::isTableSelect(model);
    layout.isAssign = NLDB0::isAssign(model);
    layout.nonOutputTermFlatIDs.reserve(model->getBitTerms().size());
    for (const auto* term : model->getBitTerms()) {
      ++layout.bitTermCount;
      if (term->getDirection() != SNLBitTerm::Direction::Output) {
        layout.nonOutputTermFlatIDs.push_back(term->getOrderID());
      }
    }
    if (layout.isMux2) {
      if (const auto* muxInputA = NLDB0::getMux2InputA(model)) {
        for (const auto* bitTerm : muxInputA->getBusBits()) {
          layout.muxInputAFlatIDs.emplace(bitTerm->getBit(),
                                          bitTerm->getOrderID());
        }
      }
      if (const auto* muxInputB = NLDB0::getMux2InputB(model)) {
        for (const auto* bitTerm : muxInputB->getBusBits()) {
          layout.muxInputBFlatIDs.emplace(bitTerm->getBit(),
                                          bitTerm->getOrderID());
        }
      }
      if (const auto* muxSelect = NLDB0::getMux2Select(model)) {
        layout.muxSelectFlatID = muxSelect->getOrderID();
      }
    }
  }

  const auto inserted =
      modelInputLayoutCache.emplace(key, std::move(layout));
  return inserted.first->second;
}

size_t getTruthTableCountCached(const SNLDesign* model) {
  const TruthTableKey key{model, nullptr, InvalidFlatTermID};
  auto it = truthTableCountCache.find(key);
  if (it != truthTableCountCache.end()) {
    return it->second;
  }
  const size_t count = SNLDesignModeling::getTruthTableCount(model);
  truthTableCountCache.emplace(key, count);
  return count;
}

const SNLTruthTable& getTruthTableCached(const SNLDesign* model,
                                         const naja::NL::SNLInstance* instance,
                                         size_t flatTermID) {
  const bool usesInstanceTable =
      SNLDesignModeling::hasTruthTableFromParameter(model, flatTermID);
  const TruthTableKey key{
      model, usesInstanceTable ? instance : nullptr, flatTermID};
  const auto& it = truthTableCache.find(key);
  if (it != truthTableCache.end()) {
    return it->second;
  }
  auto tt = usesInstanceTable
                ? SNLDesignModeling::getTruthTable(instance, flatTermID)
                : SNLDesignModeling::getTruthTable(model, flatTermID);
  const auto& entry = truthTableCache.emplace(key, std::move(tt));
  return entry.first->second;
}

bool canUseCachedIsoShortcut(const naja::DNL::DNLIso& iso,
                             naja::DNL::DNLID driver) {
  if (iso.getIsoID() == naja::DNL::DNLID_MAX || iso.getDrivers().empty() ||
      iso.getDrivers().front() != driver) {
    return false;
  }
  const auto cached = Tree2BoolExpr::iso2boolExpr_.find(iso.getIsoID());
  return cached != Tree2BoolExpr::iso2boolExpr_.end() &&
         cached->second != nullptr && cached->second->isValid();
}

uint32_t nextExpandedTableTermEpoch() {
  ++expandedTableTermEpoch;
  // LCOV_EXCL_START
  if (expandedTableTermEpoch == 0) {
    // LCOV_DISABLED_START
    std::fill(expandedTableTermEpochs.begin(), expandedTableTermEpochs.end(),
              0);
    expandedTableTermEpoch = 1;
  }
  // LCOV_DISABLED_STOP
  // LCOV_EXCL_STOP
  return expandedTableTermEpoch;
}

bool markExpandedTableTermThisIteration(naja::DNL::DNLID termID,
                                        uint32_t epoch,
                                        size_t termCount) {
  const size_t index = static_cast<size_t>(termID);
  const size_t requiredSize = std::max(termCount, index + 1);
  if (expandedTableTermEpochs.size() < requiredSize) {
    expandedTableTermEpochs.resize(requiredSize, 0);
  }
  uint32_t& slot = expandedTableTermEpochs[index];
  if (slot == epoch) {
    return false;
  }
  slot = epoch;
  return true;
}

bool shouldReportSkippedPOs() {
  return KEPLER_FORMAL::Config::getReportSkippedPOs();
}

const char* kSkippedMultiDriverPOReport = "skipped_multi_driver_pos.txt";
const char* kSkippedNoDriverPOReport = "skipped_no_driver_pos.txt";
const char* kSkippedLogicalLoopPOReport = "skipped_logical_loop_pos.txt";

void initializeSkippedPOReportFiles() {
  static std::once_flag once;
  if (!shouldReportSkippedPOs()) {
    // LCOV_EXCL_START
    return; // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
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
      // LCOV_EXCL_START
      return DNLID_MAX; // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    return dnl->getDNLTerminalFromID(termID).getIsoID();
  };

  const DNLID currentIsoID = getIsoID(currentInput);
  if (currentIsoID != DNLID_MAX) {
    return currentIsoID;
  }
  // LCOV_EXCL_START
  return getIsoID(mergeTerm); // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

struct SkippedPOReportEvent {  // LCOV_EXCL_LINE
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

const char* getSnlDirectionName(naja::NL::SNLBitTerm::Direction direction) {
  switch (direction) {
    case naja::NL::SNLBitTerm::Direction::Undefined:
      // LCOV_EXCL_START
      return "Undefined";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    case naja::NL::SNLBitTerm::Direction::Input:
      return "Input";
    case naja::NL::SNLBitTerm::Direction::Output:
      return "Output";
    case naja::NL::SNLBitTerm::Direction::InOut:
      // LCOV_EXCL_START
      return "InOut";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
  }
  // LCOV_EXCL_START
  return "Unknown";  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

std::string getSnlModelName(const DNLTerminalFull& term) {
  const auto* design = term.getSnlBitTerm()->getDesign();
  return design ? design->getName().getString() : "<unknown>";
}

void appendCloudTermName(std::ostream& out, const DNLFull* dnl, DNLID termID) {
  if (termID == DNLID_MAX) {
    // LCOV_EXCL_START
    out << "<invalid>"; // LCOV_EXCL_LINE
    return; // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
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
      << " (model=" << getSnlModelName(term)
      << ", direction="
      << getSnlDirectionName(term.getSnlBitTerm()->getDirection())
      << ")"
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
    // LCOV_EXCL_START
    return; // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.getIsoID() == DNLID_MAX) {
    // LCOV_EXCL_START
    return; // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
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
    IterationInputsTLPair;

thread_local IterationInputsTLPair currentIterationInputsTL;

IterationInputsTLPair& getCurrentIterationInputsTL() {
  return currentIterationInputsTL;
}

thread_local IterationInputsTLPair newIterationInputsTL;

IterationInputsTLPair& getNewIterationInputsTL() {
  return newIterationInputsTL;
}

void clearCurrentIterationInputsTL() {
  auto& currentIterationInputs = getCurrentIterationInputsTL();
  currentIterationInputs.first.clear();
}

size_t sizeOfCurrentIterationInputsTL() {
  return getCurrentIterationInputsTL().first.size();
}

void copyCurrentIterationInputsTL(std::vector<naja::DNL::DNLID, tbb::tbb_allocator<naja::DNL::DNLID>>& res) {
  res.clear();
  auto& current = getCurrentIterationInputsTL();
  res = std::move(current.first);
}

void clearNewIterationInputsTL() {
  auto& newIterationInputs = getNewIterationInputsTL();
  newIterationInputs.first.clear();
}

size_t sizeOfNewIterationInputsTL() {
  return getNewIterationInputsTL().first.size();
}

void copyNewIterationInputsTLToCurrent() {
  auto& newIterationInputs = getNewIterationInputsTL();
  auto& currentIterationInputs = getCurrentIterationInputsTL();
  #ifdef DEBUG_CHECKS
  size_t newSize = newIterationInputs.first.size();
  #endif
  currentIterationInputs = std::move(newIterationInputs);
  #ifdef DEBUG_CHECKS
  assert(currentIterationInputs.first.size() == newSize &&
         "copyNewIterationInputsTLToCurrent: size mismatch after copy");
  #endif
}

thread_local std::pair<
    std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
                tbb::tbb_allocator<std::pair<naja::DNL::DNLID,
                                            naja::DNL::DNLID>>>,
    size_t>
    inputsToMergeTL;

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
getInputsToMergeTL() {
  return inputsToMergeTL;
}

void clearInputsToMergeTL() {
  auto& inputsToMerge = getInputsToMergeTL();
  inputsToMerge.first.clear();
}

size_t sizeOfInputsToMergeTL() {
  return getInputsToMergeTL().first.size();
}

thread_local HandledSet visitedTermsPairsTLSet;

void clearVisitedTermsPairsTL() {
  visitedTermsPairsTLSet.clear();
}

bool SNLLogicCloud::isInput(naja::DNL::DNLID termID) {
  return PIs_[termID];
}

bool SNLLogicCloud::isOutput(naja::DNL::DNLID termID) {
  return POs_[termID];
}

naja::DNL::DNLID SNLLogicCloud::getIsoIDCached(
    naja::DNL::DNLID termID,
    const std::shared_ptr<const std::vector<naja::DNL::DNLID>>&
        termIsoIDs) const {
  if (termID == naja::DNL::DNLID_MAX) {
    // LCOV_EXCL_START
    return naja::DNL::DNLID_MAX;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  const size_t index = static_cast<size_t>(termID);
  if (index < termIsoIDs->size()) {
    return (*termIsoIDs)[index];
  }
  // LCOV_EXCL_START
  return dnl_.getDNLTerminalFromID(termID).getIsoID();  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

std::string SNLLogicCloud::formatTermName(naja::DNL::DNLID termID) const {
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
}

void SNLLogicCloud::appendTermList(std::ostream& out,
                                   const TermIDVector& termIDs,
                                   size_t limit) const {
  const size_t capped = std::min(termIDs.size(), limit);
  for (size_t i = 0; i < capped; ++i) {
    out << "#" << i << "{" << formatTermName(termIDs[i]) << "}";
    if (i + 1 != capped) {
      out << ", ";
    }
  }
  if (termIDs.size() > capped) {
    // LCOV_EXCL_START
    out << ", ... +" << (termIDs.size() - capped) << " more";  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
}

// LCOV_EXCL_START
void SNLLogicCloud::appendInstNonOutputs(std::ostream& out,
                                         naja::DNL::DNLID instID,
                                         size_t limit) const {
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
      // LCOV_DISABLED_START
      out << ", ";
    }
    // LCOV_DISABLED_STOP
    out << formatTermName(termID);
    ++emitted;
    if (emitted == limit) {
      // LCOV_DISABLED_START
      out << ", ...";
      return;
      // LCOV_DISABLED_STOP
    }
  }
  if (emitted == 0) {
    // LCOV_DISABLED_START
    out << "<none>";
  }
  // LCOV_DISABLED_STOP
}

void SNLLogicCloud::appendMergeListDetailed(std::ostream& out,
                                            size_t limit) const {
  const auto& merges = getInputsToMergeTL().first;
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
    // LCOV_DISABLED_START
    out << ", ... +" << (merges.size() - capped) << " more";
  }
  // LCOV_DISABLED_STOP
}
// LCOV_EXCL_STOP

std::string SNLLogicCloud::buildIterationSnapshot(size_t iter) const {
  std::ostringstream snapshot;
  snapshot << "iter " << iter << ": current_inputs=[";
  appendTermList(snapshot, getCurrentIterationInputsTL().first);
  snapshot << "] inputs_to_merge=[";
  appendMergeListDetailed(snapshot);
  snapshot << "] next_inputs=[";
  appendTermList(snapshot, getNewIterationInputsTL().first);
  snapshot << "]";
  return snapshot.str();
}

naja::DNL::DNLID SNLLogicCloud::resolveInstanceInputTerm(
    const naja::DNL::DNLInstanceFull& inst,
    size_t flatTermID,
    naja::DNL::DNLID driver,
    const char* role) const {
  if (flatTermID == InvalidFlatTermID) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    std::ostringstream error;
    error << "SNLLogicCloud failed to resolve " << role << " for driver "
          << formatTermName(driver) << " in model "
          << inst.getSNLModel()->getName().getString();
    throw std::runtime_error(error.str());
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  const auto& termIndexes = inst.getTermIndexes();
  if (termIndexes.first == DNLID_MAX ||
      termIndexes.second < termIndexes.first ||
      flatTermID > termIndexes.second - termIndexes.first) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    std::ostringstream error;
    error << "SNLLogicCloud failed to map " << role
          << " flat_term_id=" << flatTermID << " for driver "
          << formatTermName(driver) << " in model "
          << inst.getSNLModel()->getName().getString();
    throw std::runtime_error(error.str());
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  const auto& inputTerm =
      dnl_.getDNLTerminalFromID(termIndexes.first + flatTermID);
  if (inputTerm.isNull() || inputTerm.getSnlBitTerm() == nullptr ||
      inputTerm.getSnlBitTerm()->getOrderID() != flatTermID) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    std::ostringstream error;
    error << "SNLLogicCloud failed to map " << role
          << " flat_term_id=" << flatTermID << " for driver "
          << formatTermName(driver) << " in model "
          << inst.getSNLModel()->getName().getString();
    throw std::runtime_error(error.str());
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  return inputTerm.getID();  // LCOV_EXCL_LINE
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

size_t SNLLogicCloud::getRelevantInstanceInputCount(
    naja::DNL::DNLID driver) const {
  const auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
  const auto& layout = getModelInputLayout(dnl_, inst.getSNLModel());
  if (layout.isTableSelect) {
    const auto* model = inst.getSNLModel();
    const auto& tt = getTruthTableCached(
        model,
        inst.getSNLInstance(),
        dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID());
    return naja::NL::NLBitDependencies::countBitsForVector(
        tt.getDependencies());
  }
  return layout.isMux2 ? size_t{3} : layout.nonOutputTermFlatIDs.size();
}

void SNLLogicCloud::appendRelevantInstanceInputs(
    naja::DNL::DNLID driver,
    TermIDVector& relevantTerms) const {
  const auto& driverTerm = dnl_.getDNLTerminalFromID(driver);
  const auto& inst = driverTerm.getDNLInstance();
  const auto* model = inst.getSNLModel();
  const auto& layout = getModelInputLayout(dnl_, model);
  if (layout.isTableSelect) {
    const auto& tt = getTruthTableCached(
        model, inst.getSNLInstance(), driverTerm.getSnlBitTerm()->getOrderID());
    const auto deps = naja::NL::NLBitDependencies::decodeBits(
        tt.getDependencies());
    for (size_t flatTermID : deps) {
      relevantTerms.emplace_back(resolveInstanceInputTerm(
          inst, flatTermID, driver, "table select dependency"));
    }
    return;
  }
  if (!layout.isMux2) {
    for (size_t flatTermID : layout.nonOutputTermFlatIDs) {
      relevantTerms.emplace_back(resolveInstanceInputTerm(
          inst, flatTermID, driver, "instance dependency"));
    }
    return;
  }

  if (layout.muxInputAFlatIDs.empty() || layout.muxInputBFlatIDs.empty() ||
      layout.muxSelectFlatID == InvalidFlatTermID) {
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::runtime_error(
    // LCOV_DISABLED_STOP
        "SNLLogicCloud failed to resolve wide mux inputs");
    // LCOV_EXCL_STOP
  }

  const auto driverBit =
      static_cast<NLID::Bit>(driverTerm.getSnlBitTerm()->getBit());
  auto findMuxInput = [&](const std::unordered_map<NLID::Bit, size_t>& bits,
                          const char* role) {
    const auto it = bits.find(driverBit);
    if (it == bits.end()) {
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      std::ostringstream error;
      error << "SNLLogicCloud failed to resolve " << role << " bit "
            << driverBit << " for driver " << formatTermName(driver)
            << " in model " << model->getName().getString();
      throw std::runtime_error(error.str());
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    return it->second;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  };  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP
  relevantTerms.emplace_back(resolveInstanceInputTerm(
      inst, findMuxInput(layout.muxInputAFlatIDs, "wide mux input A"), driver,
      "wide mux input A"));
  relevantTerms.emplace_back(resolveInstanceInputTerm(
      inst, findMuxInput(layout.muxInputBFlatIDs, "wide mux input B"), driver,
      "wide mux input B"));
  relevantTerms.emplace_back(resolveInstanceInputTerm(
      inst, layout.muxSelectFlatID, driver, "wide mux select"));
}

SNLLogicCloud::TermIDVector SNLLogicCloud::collectRelevantInstanceInputs(
    naja::DNL::DNLID driver) const {
  TermIDVector relevantTerms;
  relevantTerms.reserve(getRelevantInstanceInputCount(driver));
  appendRelevantInstanceInputs(driver, relevantTerms);
  return relevantTerms;
}

void SNLLogicCloud::throwIfTruthTableArityMismatch(
    naja::DNL::DNLID driver) const {
  const auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
  const auto* model = inst.getSNLModel();
  if (getTruthTableCountCached(model) == 0) {
    // LCOV_EXCL_START
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  const auto& tt = getTruthTableCached(
      model,
      inst.getSNLInstance(),
      dnl_.getDNLTerminalFromID(driver).getSnlBitTerm()->getOrderID());
  if (!tt.isInitialized()) {
    // LCOV_EXCL_START
    return;  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }
  const auto& layout = getModelInputLayout(dnl_, model);
  const size_t expectedInputCount = layout.isMux2 ? size_t{3} : tt.size();
  const size_t actualInputCount = getRelevantInstanceInputCount(driver);
  if (expectedInputCount == actualInputCount) {
    return;
  }

  std::ostringstream modelNonOutputTerms;
  bool firstModelTerm = true;
  for (size_t flatTermID : layout.nonOutputTermFlatIDs) {
    if (!firstModelTerm) {
      modelNonOutputTerms << ", ";
    }
    firstModelTerm = false;
    modelNonOutputTerms << "flat_term_id=" << flatTermID;
  }

  const auto relevantInputs = collectRelevantInstanceInputs(driver);
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
        << ", model non-output term count="
        << layout.nonOutputTermFlatIDs.size()
        << ", instance non-output term count=" << actualInputCount
        << ", driver=" << formatTermName(driver)
        << ", model non-output terms=[" << modelNonOutputTerms.str()
        << "], instance non-output terms=[" << instanceNonOutputTerms.str()
        << "], tt=" << tt.toString();
  throw std::runtime_error(error.str());
}

void SNLLogicCloud::throwIfFrontierMismatch(
    size_t iter,
    const std::vector<std::string>& frontierHistory) const {
  const size_t currentCount = sizeOfCurrentIterationInputsTL();
  const size_t mergeCount = sizeOfInputsToMergeTL();
  const size_t borderCount = table_.getBorderLeavesSize();
  if (currentCount == borderCount && mergeCount == borderCount) {
    return;
  }

  // LCOV_EXCL_START
  // LCOV_DISABLED_START
  constexpr size_t kMaxEntries = 24;
  std::ostringstream error;
  error << "SNLLogicCloud frontier mismatch before concat at iter " << iter
        << ": current_inputs=" << currentCount
        << ", inputs_to_merge=" << mergeCount
        << ", border_leaves=" << borderCount;
        // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  error << " current_inputs=[";
  const auto& current = getCurrentIterationInputsTL().first;
  const size_t currentLimit = std::min(current.size(), kMaxEntries);
  for (size_t i = 0; i < currentLimit; ++i) {
    const auto input = current[i];
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(input).getIsoID());
    error << "#" << i << "{" << formatTermName(input)
          << ", iso=" << iso.getIsoID()
          << ", drivers=" << iso.getDrivers().size()
          << ", readers=" << iso.getReaders().size() << "}";
    if (i + 1 != currentLimit) {
      error << ", ";
    }
  }
  if (current.size() > currentLimit) {
    error << ", ... +" << (current.size() - currentLimit) << " more";
  }
  error << "]";
  // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  error << " inputs_to_merge=[";
  const auto& merges = getInputsToMergeTL().first;
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
  // LCOV_DISABLED_STOP

  auto appendDuplicateTerms =
      // LCOV_DISABLED_START
      [&](const char* label, const TermIDVector& terms) {
        std::map<naja::DNL::DNLID, size_t> counts;
        for (const auto termID : terms) {
          ++counts[termID];
          // LCOV_DISABLED_STOP
        }
        // LCOV_DISABLED_START
        bool emitted = false;
        error << " " << label << "=[";
        size_t shown = 0;
        for (const auto& [termID, count] : counts) {
          if (count <= 1) {
            continue;
            // LCOV_DISABLED_STOP
          }
          // LCOV_DISABLED_START
          if (shown == kMaxEntries) {
            error << "...";
            emitted = true;
            break;
            // LCOV_DISABLED_STOP
          }
          // LCOV_DISABLED_START
          if (shown != 0) {
            error << ", ";
          }
          error << formatTermName(termID) << " x" << count;
          emitted = true;
          ++shown;
          // LCOV_DISABLED_STOP
        }
        // LCOV_DISABLED_START
        if (!emitted) {
          error << "<none>";
        }
        error << "]";
      };
      // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  TermIDVector mergeTerms;
  mergeTerms.reserve(merges.size());
  for (const auto& merge : merges) {
    mergeTerms.push_back(merge.second);
    // LCOV_DISABLED_STOP
  }
  // LCOV_DISABLED_START
  appendDuplicateTerms("duplicate_current_inputs", current);
  appendDuplicateTerms("duplicate_merge_terms", mergeTerms);
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
  // LCOV_DISABLED_STOP
  // LCOV_EXCL_STOP
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

void SNLLogicCloud::throwIfNextFrontierMismatch(
    size_t iter,
    const std::vector<std::string>& frontierHistory) const {
  const size_t nextCount = sizeOfNewIterationInputsTL();
  const size_t borderCount = table_.getBorderLeavesSize();
  if (nextCount == borderCount) {
    return;
  }
  // LCOV_EXCL_START
  // LCOV_DISABLED_START
  constexpr size_t kMaxEntries = 24;
  std::ostringstream error;
  error << "SNLLogicCloud next frontier mismatch after concat at iter "
        << iter << ": next_inputs=" << nextCount
        << ", border_leaves=" << borderCount;
        // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  error << " next_inputs=[";
  const auto& nextInputs = getNewIterationInputsTL().first;
  const size_t nextLimit = std::min(nextInputs.size(), kMaxEntries);
  for (size_t i = 0; i < nextLimit; ++i) {
    const auto input = nextInputs[i];
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        dnl_.getDNLTerminalFromID(input).getIsoID());
    error << "#" << i << "{" << formatTermName(input)
          << ", iso=" << iso.getIsoID()
          << ", drivers=" << iso.getDrivers().size()
          << ", readers=" << iso.getReaders().size() << "}";
    if (i + 1 != nextLimit) {
      error << ", ";
    }
  }
  if (nextInputs.size() > nextLimit) {
    error << ", ... +" << (nextInputs.size() - nextLimit) << " more";
  }
  error << "]";
  // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  error << " from_current_inputs=[";
  const auto& current = getCurrentIterationInputsTL().first;
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
  // LCOV_DISABLED_STOP

  // LCOV_DISABLED_START
  error << " from_inputs_to_merge=[";
  const auto& merges = getInputsToMergeTL().first;
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
  // LCOV_DISABLED_STOP
  // LCOV_EXCL_STOP
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

naja::DNL::DNLID SNLLogicCloud::resolveTransparentLoopTarget(
    naja::DNL::DNLID termID,
    const std::shared_ptr<const std::vector<naja::DNL::DNLID>>&
        termIsoIDs) const {
  auto& transparentLoopTargetCache = transparentLoopTargetCacheTL;
  const auto cachedIt = transparentLoopTargetCache.find(termID);
  if (cachedIt != transparentLoopTargetCache.end()) {
    return cachedIt->second;
  }
  auto& visitedTerms = transparentLoopVisitedTermsTL;
  visitedTerms.clear();
  naja::DNL::DNLID currentTermID = termID;
  while (currentTermID != naja::DNL::DNLID_MAX &&
         visitedTerms.insert(currentTermID).second) {
    const auto& currentTerm = dnl_.getDNLTerminalFromID(currentTermID);
    if (currentTerm.isNull()) {
      // LCOV_EXCL_START
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (currentTerm.getSnlBitTerm()->getDirection() !=
        SNLBitTerm::Direction::Output) {
      const auto isoID = getIsoIDCached(currentTermID, termIsoIDs);
      if (isoID == naja::DNL::DNLID_MAX) {
        // LCOV_EXCL_START
        break;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant() || iso.getDrivers().size() != 1) {
        break;
      }
      currentTermID = iso.getDrivers().front();
      continue;
    }

    const auto* model = currentTerm.getDNLInstance().getSNLModel();
    if (model == nullptr || !getModelInputLayout(dnl_, model).isAssign) {
      break;
    }

    if (getRelevantInstanceInputCount(currentTermID) != 1) {
      // LCOV_EXCL_START
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    const auto relevantTerms = collectRelevantInstanceInputs(currentTermID);
    currentTermID = relevantTerms.empty() ? naja::DNL::DNLID_MAX
                                          : relevantTerms.front();
  }
  transparentLoopTargetCache.emplace(termID, currentTermID);
  return currentTermID;
}

bool SNLLogicCloud::rejectOpaqueInternalOutput(naja::DNL::DNLID termID) {
  if (!stopAtOpaqueInternalOutputs_ || isInput(termID)) {
    return false;
  }

  const auto& term = dnl_.getDNLTerminalFromID(termID);
  if (term.isNull() || term.isTopPort() ||
      term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
    return false;  // LCOV_EXCL_LINE - callers pass an internal iso driver.
  }
  const auto truthTable = SNLDesignModeling::getTruthTable(
      term.getDNLInstance().getSNLInstance(),
      term.getSnlBitTerm()->getOrderID());

  std::string reason =
      "no initialized combinational truth table or usable sequential model";
  if (truthTable.isInitialized()) {
    const auto& layout =
        getModelInputLayout(dnl_, term.getDNLInstance().getSNLModel());
    const size_t expectedInputCount =
        layout.isMux2 ? size_t{3} : truthTable.size();
    const size_t actualInputCount = getRelevantInstanceInputCount(termID);
    if (expectedInputCount == actualInputCount) {
      return false;
    }

    std::ostringstream arityReason;
    arityReason
        << "combinational truth table arity does not match instance inputs "
        << "(TT arity=" << truthTable.size()
        << ", instance input count=" << actualInputCount << ")";
    reason = arityReason.str();
  } else {
    const auto relatedClocks =
        SNLDesignModeling::getOutputRelatedClocks(term.getSnlBitTerm());
    if (!relatedClocks.empty()) {
      const auto* model = term.getDNLInstance().getSNLModel();
      if (model == nullptr || !SNLDesignModeling::hasSequentialModel(model)) {
        reason = "Missing Naja sequential model";
      } else if (SNLDesignModeling::getSequentialModel(model).kind ==
                 SNLDesignModeling::SequentialModel::Kind::Latch) {
        reason = "Naja latch sequential models are not supported by SEC";
      } else {
        reason = "the sequential output has no usable SEC model";
      }
    } else {
      const auto* model = term.getDNLInstance().getSNLModel();
      if (model != nullptr && SNLDesignModeling::hasSequentialModel(model)) {
        reason = "the sequential output has no usable SEC model";
      }
    }
  }

  std::string instance = term.getDNLInstance().getFullPath();
  while (!instance.empty() &&
         (instance.back() == '/' || instance.back() == '.')) {
    instance.pop_back();
  }
  if (instance.empty()) {
    // DNL uses the instance ID when an SNL instance has no name.
    instance = "<unnamed internal instance>";  // LCOV_EXCL_LINE
  }
  std::ostringstream detail;
  detail << "opaque internal cell `" << instance << "`";
  if (const auto* model = term.getDNLInstance().getSNLModel()) {
    detail << " (model `" << model->getName().getString() << "`)";
  }
  detail << " pin `" << term.getSnlBitTerm()->getName().getString() << "["
         << term.getSnlBitTerm()->getBit() << "]`: " << reason;

  skipReason_ = SkipReason::OpaqueInternal;
  skipReasonText_ = detail.str();
  opaqueInternalTerm_ = termID;
  table_ = SNLTruthTableTree();
  return true;
}

void SNLLogicCloud::compute() {
  refreshPerDnlCaches(dnl_);
  clearNewIterationInputsTL();
  clearCurrentIterationInputsTL();
  auto& currentIterationInputs = getCurrentIterationInputsTL().first;
  auto& newIterationInputs = getNewIterationInputsTL().first;
  auto& inputsToMerge = getInputsToMergeTL().first;
  const auto termIsoIDs = getSharedTermIsoIDCache(dnl_);
  const bool captureFrontierHistory =
      std::getenv("KEPLER_CAPTURE_FRONTIER_HISTORY") != nullptr;
  std::vector<std::string> frontierHistory;
  DEBUG_LOG("---- Begin!!\n");
  if (dnl_.getDNLTerminalFromID(seedOutputTerm_).isTopPort() ||
      isOutput(seedOutputTerm_)) {
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
        getIsoIDCached(seedOutputTerm_, termIsoIDs));
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
      // LCOV_DISABLED_START
      throw std::runtime_error("Seed output term is not a single driver");
      // LCOV_DISABLED_STOP
    } else if (iso.getDrivers().empty()) {
      // LCOV_DISABLED_START
      std::string termName = dnl_.getDNLTerminalFromID(seedOutputTerm_)
                                 .getSnlBitTerm()
                                 ->getName()
                                 .getString();
                                 // LCOV_DISABLED_STOP
      std::string error =
          // LCOV_DISABLED_START
          "Seed output term '" + termName + "' has no drivers";
      throw std::runtime_error(error);
    }
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
    const auto& driver = iso.getDrivers().front();
    auto& inst = dnl_.getDNLTerminalFromID(driver).getDNLInstance();
    if (rejectOpaqueInternalOutput(driver)) {
      return;
    }
    if (isInput(driver)) {
      currentIterationInputs.emplace_back(driver);
      table_ = SNLTruthTableTree(inst.getID(), driver,
                                 SNLTruthTableTree::Node::Type::P);
      return;
    }
    throwIfTruthTableArityMismatch(driver);
    DEBUG_LOG("Instance name: %s\n",
              inst.getSNLInstance()->getName().getString().c_str());
    appendRelevantInstanceInputs(driver, newIterationInputs);
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    table_ = SNLTruthTableTree(inst.getID(), driver);
    auto* model = inst.getSNLModel();
    assert(getTruthTableCached(
               model,
               inst.getSNLInstance(),
               dnl_.getDNLTerminalFromID(driver)
                   .getSnlBitTerm()
                   ->getOrderID())
               .isInitialized() &&
        "Truth table is not initialized");
    assert(table_.isInitialized() &&
           "Truth table for seed output term is not initialized");
  } else {
    // LCOV_EXCL_START
    const auto& inst = dnl_.getDNLInstanceFromID(seedOutputTerm_);  // LCOV_EXCL_LINE
    for (DNLID termID = inst.getTermIndexes().first;  // LCOV_EXCL_LINE
         termID <= inst.getTermIndexes().second; termID++) {  // LCOV_EXCL_LINE
      const DNLTerminalFull& term = dnl_.getDNLTerminalFromID(termID);  // LCOV_EXCL_LINE
      if (term.getSnlBitTerm()->getDirection() !=  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          SNLBitTerm::Direction::Output) {
        // newIterationInputs.emplace_back(termID);
        // LCOV_EXCL_START
        newIterationInputs.emplace_back(termID);  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        DEBUG_LOG("Add input with id: %zu\n", termID);
      // LCOV_EXCL_START
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    DEBUG_LOG("model name: %s\n",
              inst.getSNLModel()->getName().getString().c_str());
    // LCOV_EXCL_START
    table_ = SNLTruthTableTree(inst.getID(), seedOutputTerm_);  // LCOV_EXCL_LINE
    assert(table_.isInitialized() &&  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
           "Truth table for seed output term is not initialized");
  }

  if (newIterationInputs.empty()) {
    DEBUG_LOG("No inputs found for seed output term %zu\n", seedOutputTerm_);
    return;  // LCOV_EXCL_LINE
  }
  // LCOV_EXCL_START
  if (captureFrontierHistory) {
    std::ostringstream seedInfo;
    seedInfo << "seed_output={" << formatTermName(seedOutputTerm_) << "}"
             << " initial_inputs=[";
    appendTermList(seedInfo, getNewIterationInputsTL().first);
    seedInfo << "]";
    frontierHistory.emplace_back(seedInfo.str());
  }
  // LCOV_EXCL_STOP

  bool reachedPIs = true;
  size_t size = newIterationInputs.size();
  
  for (size_t i = 0; i < size; i++) {
    if (!isInput(
            newIterationInputs
                [i]) /* && !isOutput(newIterationInputs[i])*/) {
      reachedPIs = false;
      break;
    }
  // LCOV_EXCL_START
  } // allocator for buckets LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  // HandledSet handledTerms;
  // handledTerms.reserve(naja::DNL::get()->getDNLTerms().size() / 4);
  clearVisitedTermsPairsTL();
  size_t iter = 0;
  transparentLoopTargetCacheTL.clear();

  while (!reachedPIs) {
    // Originally computation of reachedPIs have been handled in the end of the loop,
    // but by adding isConstant on isos as part of the check, it had to be moved to the beginning.
    // Why? because before we cached the inputs of the leaves we meet and then we look at the drivers in the next loop iteration
    // and then we run the checks on the drivers in order to know if we need to iterate again, after the drivers were contacted to the cloud.
    // Now, we also check isos of inputs, which means that if an input has an iso that is constant, we will stop.
    // In this case, we have to make the check in the next loop iteration in order to force concating the constants to the cloud.
    reachedPIs = true;
    size_t sizeOfNewInputs = newIterationInputs.size();
    for (size_t i = 0; i < sizeOfNewInputs; i++) {
      const auto input = newIterationInputs[i];
      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          getIsoIDCached(input, termIsoIDs));
      const bool cachedAsInput = canUseCachedIsoShortcut(iso, input);
      if (!isInput(input) && !cachedAsInput && !iso.isConstant()) {
        reachedPIs = false;
        break;
      }
    }
    DEBUG_LOG("---iter %lu---\n", iter);
    DEBUG_LOG("Current iteration inputs size: %zu\n",
              newIterationInputs.size());
    copyNewIterationInputsTLToCurrent();

    clearNewIterationInputsTL();
    DEBUG_LOG("table size: %zu, currentIterationInputs_ size: %zu\n",
              table_.size(), currentIterationInputs.size());
    clearInputsToMergeTL();
    const uint32_t expandedTableTermEpochThisIteration =
        nextExpandedTableTermEpoch();
    size_t sizeOfCurrentInputs = currentIterationInputs.size();
    inputsToMerge.reserve(sizeOfCurrentInputs);
    newIterationInputs.reserve(sizeOfCurrentInputs);
    transparentLoopTargetCacheTL.reserve(transparentLoopTargetCacheTL.size() +
                                         sizeOfCurrentInputs);
    for (size_t i = 0; i < sizeOfCurrentInputs; i++) {
      const auto& input = currentIterationInputs[i];
      const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
          getIsoIDCached(input, termIsoIDs));
      if (isInput(input) || canUseCachedIsoShortcut(iso, input) ||
          iso.isConstant()) {
        newIterationInputs.emplace_back(input);
        DEBUG_LOG("Adding input id: %zu %s\n", input,
                  dnl_.getDNLTerminalFromID(input)
                      .getSnlBitTerm()
                      ->getName()
                      .getString()
                      .c_str());
        inputsToMerge.emplace_back(naja::DNL::DNLID_MAX,
                                   input);  // Placeholder for PI/PO
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

      if (iso.getDrivers().size() >= 1) {
        // proper error with names of all the drivers
        // throw an error and separate names by comma
        if (iso.getDrivers().size() > 1) {
          skipReason_ = SkipReason::MultiDriver;
          skipReasonText_ = "its iso has multiple drivers during cloud expansion";
          reportCloudSkippedRoot(&dnl_, seedOutputTerm_, input, DNLID_MAX,
                                 "its iso has multiple drivers during cloud expansion",
                                 kSkippedMultiDriverPOReport);
          table_ = SNLTruthTableTree();
          return;
        }
      } else if (iso.getDrivers().empty()) {
        if (!iso.isConstant()) {
          skipReason_ = SkipReason::NoDriver;
          skipReasonText_ = "its iso has no drivers during cloud expansion";
          reportCloudSkippedRoot(&dnl_, seedOutputTerm_, input, DNLID_MAX,
                                 "its iso has no drivers during cloud expansion",
                                 kSkippedNoDriverPOReport);
          table_ = SNLTruthTableTree();
          return;
        }
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        newIterationInputs.emplace_back(input);
        inputsToMerge.emplace_back(naja::DNL::DNLID_MAX,
                                   input);  // Placeholder for PI/PO
        continue;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }
      const auto& driver = iso.getDrivers().front();
      if (rejectOpaqueInternalOutput(driver)) {
        return;
      }
      if (isInput(driver) || canUseCachedIsoShortcut(iso, driver)) {
        newIterationInputs.emplace_back(driver);
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
        inputsToMerge.emplace_back(naja::DNL::DNLID_MAX,
                                   driver);  // Placeholder for PI/PO
        continue;
      }

      const auto& driverTerm = dnl_.getDNLTerminalFromID(driver);
      if (driverTerm.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        // LCOV_EXCL_START
        // LCOV_DISABLED_START
        skipReason_ = SkipReason::NoDriver;
        skipReasonText_ =
        // LCOV_DISABLED_STOP
            "encountered an internal frontier that was not collected as a PI";
        // LCOV_DISABLED_START
        reportCloudSkippedRoot(
            &dnl_, seedOutputTerm_, input, driver,
            // LCOV_DISABLED_STOP
            "encountered an internal frontier that was not collected as a PI",
            // LCOV_DISABLED_START
            kSkippedNoDriverPOReport);
        table_ = SNLTruthTableTree();
        return;
        // LCOV_DISABLED_STOP
        // LCOV_EXCL_STOP
      }

      const auto& inst =
          dnl_.getDNLInstanceFromID(driverTerm.getDNLInstance().getID());
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
      inputsToMerge.emplace_back(inst.getID(), driver);

      const bool alreadyInTree = table_.hasTableTerm(driver);
      const bool firstExpansionInThisBatch =
          markExpandedTableTermThisIteration(
              driver, expandedTableTermEpochThisIteration, dnl_.getNBterms());
      const bool tableWillExpand =
          !alreadyInTree && firstExpansionInThisBatch;
      if (!tableWillExpand) {
        // Reused truth-table nodes keep their existing children, so adding the
        // same instance inputs again would desynchronize the Boolean frontier
        // from the tree border leaves.
        continue;
      }

      appendRelevantInstanceInputs(driver, newIterationInputs);
    }

    if (inputsToMerge.empty()) {
      // LCOV_EXCL_START
      break;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    DEBUG_LOG("--- Merging truth tables with %zu inputs\n",
              inputsToMerge.size());
    throwIfFrontierMismatch(iter, frontierHistory);
    {
      const auto& merges = getInputsToMergeTL().first;
      const auto& currentInputs = getCurrentIterationInputsTL().first;
      for (size_t i = 0; i < merges.size(); ++i) {
        if (merges[i].first == naja::DNL::DNLID_MAX) {
          continue;
        }
        auto loopTargetTerm = merges[i].second;
        if (!table_.hasTableTerm(loopTargetTerm)) {
          const auto& mergeTerm = dnl_.getDNLTerminalFromID(loopTargetTerm);
          const auto* model = mergeTerm.getDNLInstance().getSNLModel();
          if (model == nullptr ||
              !getModelInputLayout(dnl_, model).isAssign) {
            continue;
          }
          loopTargetTerm =
              resolveTransparentLoopTarget(loopTargetTerm, termIsoIDs);
          if (!table_.hasTableTerm(loopTargetTerm)) {
            continue;
          }
        }
        auto& loopTerms = loopTermsScratchTL;
        loopTerms.clear();
        // Normal cloud expansion must stay conservative: if this merge would
        // reconnect to any transparent alias already above the border leaf, the
        // cone contains a combinational loop. SEC's structured-memory path can
        // opt into cutting a specific internal loop frontier before reaching
        // this strict reporting path.
        if (table_.findAncestorLoopForBorderLeaf(
                i, loopTargetTerm, loopTerms)) {
          skipReason_ = SkipReason::LogicalLoop;
          skipReasonText_ =
              "a logical loop was detected during cloud expansion";
          reportCloudSkippedRoot(
              &dnl_, seedOutputTerm_, currentInputs[i], merges[i].second,
              "a logical loop was detected during cloud expansion",
              kSkippedLogicalLoopPOReport, &loopTerms);
          table_ = SNLTruthTableTree();
          return;
        }
      }
    }
    table_.concatFull(inputsToMerge, inputsToMerge.size());
    throwIfNextFrontierMismatch(iter, frontierHistory);
    // LCOV_EXCL_START
    if (captureFrontierHistory) {
      frontierHistory.emplace_back(buildIterationSnapshot(iter));
    }
    // LCOV_EXCL_STOP
    DEBUG_LOG("--- End of iteration %zu\n", iter);
    iter++;
  }

  copyNewIterationInputsTLToCurrent();
  #ifdef DEBUG_CHECKS
  size_t finalSize = currentIterationInputs.size();
  #endif
  copyCurrentIterationInputsTL(currentIterationInputs_);
  #ifdef DEBUG_CHECKS
  assert(finalSize == currentIterationInputs_.size() &&
         "compute: size mismatch after final copy");
  for (const auto& input : currentIterationInputs_) {
    const auto& iso = dnl_.getDNLIsoDB().getIsoFromIsoIDconst(
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
      // LCOV_EXCL_START
      continue; // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    auto& events = *eventsPtr;
    for (const auto& event : events) {
      std::ostream* out = &logicalLoopOut;
      const std::string_view reportFile(event.reportFile);
      if (reportFile == std::string_view(kSkippedMultiDriverPOReport)) {
        // LCOV_EXCL_START
        out = &multiDriverOut; // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      } else if (reportFile == std::string_view(kSkippedNoDriverPOReport)) {
        out = &noDriverOut;
      }

      if (!(*out)) {
        // LCOV_EXCL_START
        continue; // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
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
        // LCOV_DISABLED_START
        if (event.reportIsoID != DNLID_MAX) {
          *out << ". See first encounter of iso=" << event.reportIsoID
               << " for details";
        }
        *out << "\n\n";
        continue;
        // LCOV_DISABLED_STOP
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
