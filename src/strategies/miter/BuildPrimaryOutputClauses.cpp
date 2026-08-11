// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "BuildPrimaryOutputClauses.h"
#include "DNL.h"
#include "NLDB0.h"
#include "SNLDesignModeling.h"
#include "SNLLogicCloud.h"
#include "Tree2BoolExpr.h"
#include "SNLPath.h"
#include "NajaProperty.h"
#include "../../config/Config.h"
#include <algorithm>
#include <fstream>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string_view>
#include <thread>
#include <tbb/global_control.h>

#ifdef DEBUG_PRINTS
#define DEBUG_LOG(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_LOG(fmt, ...)
#endif

using namespace KEPLER_FORMAL;
using namespace naja::DNL;
using namespace naja::NL;

namespace {

constexpr BuildPrimaryOutputClauses::PathComponentID kUnnamedPathComponentTag =
    BuildPrimaryOutputClauses::PathComponentID{1} << 63;

const char* getSnlDirectionName(SNLBitTerm::Direction direction) {
  switch (direction) {
    case SNLBitTerm::Direction::Undefined:
      // LCOV_EXCL_START
      return "Undefined";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    case SNLBitTerm::Direction::Input:
      return "Input";
    case SNLBitTerm::Direction::Output:
      return "Output";
    case SNLBitTerm::Direction::InOut:
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

void appendTerminalName(std::ostream& out, const DNLTerminalFull& term) {
  const auto path = term.getDNLInstance().getPath().getPathNames();
  for (size_t i = 0; i < path.size(); ++i) {
    out << path[i].getString() << ".";
  }
  out << term.getSnlBitTerm()->getName().getString()
      << term.getSnlBitTerm()->getBit()
      << " (model=" << getSnlModelName(term)
      << ", direction="
      << getSnlDirectionName(term.getSnlBitTerm()->getDirection())
      << ")";
}

bool shouldReportSkippedPOs() {
  return KEPLER_FORMAL::Config::getReportSkippedPOs();
}

const char* kSkippedMultiDriverPOReport = "skipped_multi_driver_pos.txt";
const char* kSkippedNoDriverPOReport = "skipped_no_driver_pos.txt";

	struct SparseReportedIDs {
	  bool mark(const DNLFull* dnl, size_t nBits, DNLID id) {
	    if (id == DNLID_MAX) {
	      // LCOV_EXCL_START
	      return true; // LCOV_EXCL_LINE
	      // LCOV_EXCL_STOP
	    }
	    if (owner != dnl) {
	      for (uint32_t wordIndex : touchedWords) {
	        // LCOV_EXCL_START
	        words[wordIndex] = 0; // LCOV_EXCL_LINE
	        // LCOV_EXCL_STOP
	      }
      touchedWords.clear();
      owner = dnl;
    }
    ensureCapacity(std::max(nBits, static_cast<size_t>(id + 1)));
    const size_t wordIndex = id >> 6;
    const uint64_t mask = uint64_t{1} << (id & 63);
    uint64_t& word = words[wordIndex];
    if (word & mask) {
      return false;
    }
    if (word == 0) {
      touchedWords.push_back(static_cast<uint32_t>(wordIndex));
    }
    word |= mask;
    return true;
  }

 private:
  void ensureCapacity(size_t nBits) {
    const size_t wordCount = (nBits + 63) / 64;
    if (words.size() < wordCount) {
      words.resize(wordCount, 0);
    }
  }

  const DNLFull* owner = nullptr;
  std::vector<uint64_t> words;
  std::vector<uint32_t> touchedWords;
};

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
  });
}

bool shouldEmitSkippedPOReport(const DNLFull* dnl,
                               DNLID isoID) {
  static std::mutex mutex;
  static SparseReportedIDs reportedIsos;

	  if (isoID == DNLID_MAX) {
	    // LCOV_EXCL_START
	    return true; // LCOV_EXCL_LINE
	    // LCOV_EXCL_STOP
	  }

  std::lock_guard<std::mutex> lock(mutex);
  return reportedIsos.mark(dnl, dnl->getDNLIsoDB().getNumIsos() + 1, isoID);
}

std::mutex& getReportWriteMutex(const char* reportFile) {
  static std::mutex multiDriverMutex;
  static std::mutex noDriverMutex;
  return std::string_view(reportFile) == kSkippedMultiDriverPOReport
             ? multiDriverMutex
             : noDriverMutex;
}

void appendIsoTermName(std::ostream& out, const DNLFull* dnl, DNLID termID) {
  const auto& term = dnl->getDNLTerminalFromID(termID);
  out << term.getDNLInstance().getSNLModel()->getName().getString() << ":";
  const auto path = term.getDNLInstance().getPath().getPathNames();
  for (size_t i = 0; i < path.size(); ++i) {
    out << path[i].getString() << ".";
  }
  out << term.getSnlBitTerm()->getName().getString()
      << term.getSnlBitTerm()->getBit()
      << " (model=" << getSnlModelName(term)
      << ", direction="
      << getSnlDirectionName(term.getSnlBitTerm()->getDirection())
      << ")";
}

void appendNetReport(std::ostream& out, const SNLBitNet* net) {
  out << "design=" << net->getDesign()->getName().getString()
      << " name=" << net->getName().getString()
      << " type=" << net->getType().getString()
      << " is_assign=" << (net->getType().isAssign() ? "true" : "false")
      << " is_supply=" << (net->getType().isSupply() ? "true" : "false")
      << " is_constant0=" << (net->isConstant0() ? "true" : "false")
      << " is_constant1=" << (net->isConstant1() ? "true" : "false")
      << " model_is_assign=" << (net->getDesign()->isAssign() ? "true" : "false")
      << " properties=[";
  bool first = true;
  for (auto* property : net->getProperties()) {
    if (!first) {
      out << ", ";
    }
    first = false;
    out << property->getName() << "=" << property->getString();
  }
  out << "]";
}

template <typename TermIDs>
void appendTermsToReport(std::ostream& out,
                         const DNLFull* dnl,
                         const char* label,
                         const TermIDs& termIDs) {
  out << label << ": [";
  for (size_t i = 0; i < termIDs.size(); ++i) {
    appendIsoTermName(out, dnl, termIDs[i]);
    if (i + 1 != termIDs.size()) {
      out << ", ";
    }
  }
  out << "]";
}

void appendNetsToReport(std::ostream& out,
                        const char* label,
                        const std::set<SNLBitNet*>& nets) {
  out << label << ": [";
  size_t i = 0;
  for (const auto* net : nets) {
    appendNetReport(out, net);
    // LCOV_EXCL_START
    if (++i != nets.size()) {
      // LCOV_DISABLED_START
      out << ", ";
    }
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }
  out << "]";
}

void mergeDependencies(std::vector<uint64_t>& merged,
                       const std::vector<uint64_t>& deps) {
  if (merged.size() < deps.size()) {
    merged.resize(deps.size(), 0);
  }
  for (size_t i = 0; i < deps.size(); ++i) {
    merged[i] |= deps[i];
  }
}

bool containsDependencyBit(const std::vector<uint64_t>& deps, uint64_t orderID) {
  const size_t chunkIndex = orderID / 64;
  if (chunkIndex >= deps.size()) {
    return false;
  }
  return (deps[chunkIndex] & (uint64_t{1} << (orderID % 64))) != 0;
}

BuildPrimaryOutputClauses::SkippedOutputInfo makeSkippedOutputInfo(
    BuildPrimaryOutputClauses::SkippedOutputReason reason,
    std::string detail) {
  return {reason, std::move(detail)};
}

void reportSkippedPO(const DNLFull* dnl,
                     const DNLTerminalFull& term,
                     const char* reason,
                     const char* reportFile) {
  if (!shouldReportSkippedPOs()) {
    return;
  }
  initializeSkippedPOReportFiles();
  const bool isFirstForIso = shouldEmitSkippedPOReport(dnl, term.getIsoID());
  std::lock_guard<std::mutex> lock(getReportWriteMutex(reportFile));
  std::ofstream out(reportFile, std::ios::app);
  if (out) {
    out << "Skipping PO ";
    appendTerminalName(out, term);
    out << " of model "
        << term.getSnlBitTerm()->getDesign()->getName().getString()
        << " because " << reason;

    if (term.getIsoID() != DNLID_MAX) {
      out << ". iso=" << term.getIsoID();
    }

    if (!isFirstForIso) {
      if (term.getIsoID() != DNLID_MAX) {
        out << ". See first encounter of iso=" << term.getIsoID()
            << " for details";
      }
    } else if (term.getIsoID() != DNLID_MAX) {
      out << " ";
      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
      appendTermsToReport(out, dnl, "readers", iso.getReaders());
      out << " ";
      appendTermsToReport(out, dnl, "drivers", iso.getDrivers());
      if (std::string_view(reportFile) == kSkippedNoDriverPOReport) {
        naja::DNL::DNLComplexIso complexIso(term.getIsoID());
        dnl->getCustomIso(term.getIsoID(), complexIso);
        out << " ";
        appendTermsToReport(out, dnl, "complex_readers", complexIso.getReaders());
        out << " ";
        appendTermsToReport(out, dnl, "complex_drivers", complexIso.getDrivers());
        out << " ";
        appendTermsToReport(out, dnl, "complex_hier_terms", complexIso.getHierTerms());
        out << " ";
        appendNetsToReport(out, "complex_nets", complexIso.getNets());
      }
    }
    out << "\n\n";
  }
}

}  // namespace

BuildPrimaryOutputClauses::PathNameIDs BuildPrimaryOutputClauses::getPathNameIDs(
    const DNLInstanceFull& instance) const {
  const auto instanceID = instance.getID();
  {
    std::lock_guard<std::mutex> lock(pathNameIDsCacheMutex_);
    const auto cached = instancePathNameIDsCache_.find(instanceID);
    if (cached != instancePathNameIDsCache_.end()) {
      return cached->second;
    }
  }

  PathNameIDs ids;
  const auto pathInstances = instance.getPath().getInstances();
  ids.reserve(pathInstances.size());
  for (const auto* pathInstance : pathInstances) {
    const auto& name = pathInstance->getName();
    if (!name.empty()) {
      ids.push_back(static_cast<BuildPrimaryOutputClauses::PathComponentID>(
          name.getID()));
      continue;
    }

    ids.push_back(kUnnamedPathComponentTag |
                  static_cast<BuildPrimaryOutputClauses::PathComponentID>(
                      pathInstance->getID()));
  }

  std::lock_guard<std::mutex> lock(pathNameIDsCacheMutex_);
  const auto [cached, _] =
      instancePathNameIDsCache_.emplace(instanceID, std::move(ids));
  return cached->second;
}

BuildPrimaryOutputClauses::PathKey
BuildPrimaryOutputClauses::getTerminalPathKey(
    const DNLTerminalFull& terminal) const {
  auto pathIDs = getPathNameIDs(terminal.getDNLInstance());
  pathIDs.push_back(static_cast<PathComponentID>(
      terminal.getSnlBitTerm()->getName().getID()));
  PathObjectIDs objectIDs = {
      static_cast<NLID::DesignObjectID>(terminal.getSnlBitTerm()->getBit())};
  return {std::move(pathIDs), std::move(objectIDs)};
}

std::vector<DNLID> BuildPrimaryOutputClauses::collectInputs() {
  std::vector<DNLID> inputs;
  auto dnl = get();
  DNLInstanceFull top = dnl->getTop();

  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX && termId <= top.getTermIndexes().second; termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
      DEBUG_LOG("Collecting input %s\n",
                term.getSnlBitTerm()->getName().getString().c_str());
      assert(termId < naja::DNL::get()->getDNLTerms().size());
      inputs.emplace_back(termId);
    }
  }

  for (DNLID leaf : dnl->getLeaves()) {
    auto iter = modelCache_.find(dnl->getDNLInstanceFromID(leaf).getSNLModel());
    const DNLInstanceFull& instance = dnl->getDNLInstanceFromID(leaf);
    if ((iter != modelCache_.end()) && iter->second.analyzedPIs) {
      const auto& cache = iter->second;
      if (cache.PIs.empty()) {
        continue;
      }
      for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
         const auto& term = dnl->getDNLTerminalFromID(termId);
          if (cache.PIs.find(term.getSnlBitTerm()) != cache.PIs.end()) {
            inputs.emplace_back(termId);
          }
      }
      continue;
    }
    modelCache_[instance.getSNLModel()].analyzedPIs = true;
    size_t numberOfInputs = 0, numberOfOutputs = 0;
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Output) {
        numberOfInputs++;
      }
      if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input) {
        numberOfOutputs++;
      }
    }

    if (numberOfInputs == 0 && numberOfOutputs > 0) {
      // no inputs primtive, assuming generator so adding as PI
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          assert(termId < naja::DNL::get()->getDNLTerms().size());
          inputs.emplace_back(termId);
          modelCache_[instance.getSNLModel()].PIs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
      continue;
    }

    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      auto related =
          SNLDesignModeling::getClockRelatedOutputs(term.getSnlBitTerm());
      if (!related.empty()) {
        isSequential = true;
        for (auto bitTerm : related) {
          seqBitTerms.emplace_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          // LCOV_EXCL_START
          assert(termId < naja::DNL::get()->getDNLTerms().size());  // LCOV_EXCL_LINE
          inputs.emplace_back(termId);  // LCOV_EXCL_LINE
          modelCache_[instance.getSNLModel()].PIs.insert(  // LCOV_EXCL_LINE
              term.getSnlBitTerm());  // LCOV_EXCL_LINE
              // LCOV_EXCL_STOP
          DEBUG_LOG(
              "Collecting seq input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        // LCOV_EXCL_START
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
    }
    if (!isSequential) {
      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Input) {
          const auto tt = SNLDesignModeling::getTruthTable(
              instance.getSNLInstance(), term.getSnlBitTerm()->getOrderID());
          if (!tt.isInitialized()) {
            assert(termId < naja::DNL::get()->getDNLTerms().size());
            inputs.emplace_back(termId);
            modelCache_[instance.getSNLModel()].PIs.insert(
                term.getSnlBitTerm());
            DEBUG_LOG("Collecting input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
          
          if (tt.all0() ||
              tt.all1()) {
            assert(termId < naja::DNL::get()->getDNLTerms().size());
            inputs.emplace_back(termId);
            modelCache_[instance.getSNLModel()].PIs.insert(
                term.getSnlBitTerm());
            DEBUG_LOG("Collecting constant input %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
        }
      }
      continue;
    }
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Input) {
        if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                      term.getSnlBitTerm()) != seqBitTerms.end()) {
          assert(termId < naja::DNL::get()->getDNLTerms().size());
          inputs.emplace_back(termId);
          modelCache_[instance.getSNLModel()].PIs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting seq input %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
    }
  }
  std::set<DNLID> inputSet(inputs.begin(), inputs.end());
  inputs.clear();
  inputs.assign(inputSet.begin(), inputSet.end());
  DEBUG_LOG("Collected %zu inputs\n", inputs.size());
  return inputs;
}

std::vector<DNLID> BuildPrimaryOutputClauses::collectOutputs() {
  std::vector<DNLID> outputs;
  std::set<DNLID> outputsSet;
  skippedOutputs_.clear();
  auto dnl = get();
  DNLInstanceFull top = dnl->getTop();

  for (DNLID termId = top.getTermIndexes().first;
       termId != DNLID_MAX && termId <= top.getTermIndexes().second; termId++) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
    if (term.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input) {
      outputsSet.insert(termId);
      DEBUG_LOG(
          "Collecting top output %s of model %s\n",
          term.getSnlBitTerm()->getName().getString().c_str(),
          term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
    }
  }
  for (DNLID leaf : dnl->getLeaves()) {
    const DNLInstanceFull& instance = dnl->getDNLInstanceFromID(leaf);
    auto iter = modelCache_.find(instance.getSNLModel());
    if ((iter != modelCache_.end()) && iter->second.analyzedPOs) {
      const auto& cache = iter->second;
      if (cache.POs.empty()) {
        continue;
      }
      for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
         const auto& term = dnl->getDNLTerminalFromID(termId);
          if (cache.POs.find(term.getSnlBitTerm()) != cache.POs.end()) {
            outputsSet.insert(termId);
          }
      }
      continue;
    }
    modelCache_[instance.getSNLModel()].analyzedPOs = true;
    bool isSequential = false;
    std::vector<SNLBitTerm*> seqBitTerms;

    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      auto related =
          SNLDesignModeling::getClockRelatedInputs(term.getSnlBitTerm());
      if (!related.empty()) {
        isSequential = true;
        for (auto bitTerm : related) {
          seqBitTerms.emplace_back(bitTerm);
        }
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          outputsSet.insert(termId);
          modelCache_[instance.getSNLModel()].POs.insert(
              term.getSnlBitTerm());
          DEBUG_LOG(
              "Collecting seq output %s of model %s\n",
              term.getSnlBitTerm()->getName().getString().c_str(),
              term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
        }
      }
    }

    if (!isSequential) {
      if (NLDB0::isMux2(instance.getSNLModel())) {
        continue;
      }

      auto& cache = modelCache_[instance.getSNLModel()];
      if (!cache.analyzedPODependencies) {
        for (DNLID tId = instance.getTermIndexes().first;
             tId != DNLID_MAX && tId <= instance.getTermIndexes().second;
             tId++) {
          const DNLTerminalFull& tTerm = dnl->getDNLTerminalFromID(tId);
          if (tTerm.getSnlBitTerm()->getDirection() ==
              SNLBitTerm::Direction::Input) {
            continue;
          }
          const auto tt = SNLDesignModeling::getTruthTable(
              tTerm.getSnlBitTerm()->getDesign(),
              tTerm.getSnlBitTerm()->getOrderID());
          if (tt.isInitialized()) {
            mergeDependencies(cache.poTruthTableDependencies,
                              tt.getDependencies());
          }
        }
        cache.analyzedPODependencies = true;
      }

      for (DNLID termId = instance.getTermIndexes().first;
           termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
           termId++) {
        const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
        if (term.getSnlBitTerm()->getDirection() !=
            SNLBitTerm::Direction::Output) {
          if (!containsDependencyBit(cache.poTruthTableDependencies,
                                     term.getSnlBitTerm()->getOrderID())) {
            outputsSet.insert(termId);
            modelCache_[instance.getSNLModel()].POs.insert(
                term.getSnlBitTerm());
            DEBUG_LOG("Collecting output %s of model %s\n",
                      term.getSnlBitTerm()->getName().getString().c_str(),
                      term.getSnlBitTerm()
                          ->getDesign()
                          ->getName()
                          .getString()
                          .c_str());
          }
        }
      }
      continue;
    }
    for (DNLID termId = instance.getTermIndexes().first;
         termId != DNLID_MAX && termId <= instance.getTermIndexes().second;
         termId++) {
      const DNLTerminalFull& term = dnl->getDNLTerminalFromID(termId);
      if (term.getSnlBitTerm()->getDirection() !=
          SNLBitTerm::Direction::Output) {
        if (std::find(seqBitTerms.begin(), seqBitTerms.end(),
                      term.getSnlBitTerm()) != seqBitTerms.end()) {
          outputsSet.insert(termId);
          modelCache_[instance.getSNLModel()].POs.insert(
              term.getSnlBitTerm());
        }
        DEBUG_LOG(
            "Collecting seq output %s of model %s\n",
            term.getSnlBitTerm()->getName().getString().c_str(),
            term.getSnlBitTerm()->getDesign()->getName().getString().c_str());
      }
    }
  }
  outputs.clear();
  // keep only terminals who are connected to nets
  for (const auto& out : outputsSet) {
    const DNLTerminalFull& term = dnl->getDNLTerminalFromID(out);
    if (term.getIsoID() == DNLID_MAX) {
      DEBUG_LOG("Skipping output %s of model %s as it is not connected to any net\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    if (term.getIsoID() != DNLID_MAX && 
      dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID()).getDrivers().empty()) {
      skippedOutputs_[out] = makeSkippedOutputInfo(
          SkippedOutputReason::NoDriver, "its iso has no drivers");
      reportSkippedPO(
          dnl, term, "its iso has no drivers", kSkippedNoDriverPOReport);
      DEBUG_LOG("Skipping output %s of model %s as it is not connected to any net\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    if (term.getIsoID() != DNLID_MAX && 
      dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID()).getDrivers().size() > 1) {
      skippedOutputs_[out] = makeSkippedOutputInfo(
          SkippedOutputReason::MultiDriver, "its iso has multiple drivers");
      reportSkippedPO(
          dnl,
          term,
          "its iso has multiple drivers",
          kSkippedMultiDriverPOReport);
      DEBUG_LOG("Skipping output %s of model %s as it is driven by multiple drivers\n",
                term.getSnlBitTerm()->getName().getString().c_str(),
                term.getSnlBitTerm()
                    ->getDesign()
                    ->getName()
                    .getString()
                    .c_str());
      continue;
    }
    outputs.emplace_back(out);
  }
  return outputs;
}

void BuildPrimaryOutputClauses::collect() {
  inputs_ = collectInputs();
  for (const auto& input : inputs_) {
    PathKey key = getTerminalPathKey(naja::DNL::get()->getDNLTerminalFromID(input));
    inputsMap_[std::move(key)]  =
            input;
  }
  outputs_ = collectOutputs();
  for (const auto& output : outputs_) {
    PathKey key = getTerminalPathKey(naja::DNL::get()->getDNLTerminalFromID(output));
    outputsMap_[std::move(key)]  =
            output;
    DEBUG_LOG("Output collected: %s\n", naja::DNL::get()
                                         ->getDNLTerminalFromID(output)
                                         .getSnlBitTerm()
                                         ->getName()
                                         .getString()
                                         .c_str());
  }
  POs_.resize(outputs_.size());
}

std::vector<BuildPrimaryOutputClauses::PathKey>
BuildPrimaryOutputClauses::getLecBoundaryInputs() const {
  std::vector<PathKey> boundaries;
  for (const auto& [path, input] : inputsMap_) {
    const auto& term = get()->getDNLTerminalFromID(input);
    if (term.isTopPort() ||
        !SNLDesignModeling::getOutputRelatedClocks(term.getSnlBitTerm()).empty()) {
      boundaries.emplace_back(path);
    }
  }
  return boundaries;
}

void BuildPrimaryOutputClauses::initVarNames() {
  termDNLID2varID_.resize(naja::DNL::get()->getDNLTerms().size(), (size_t)-1);
  for (size_t i = 0; i < inputs_.size(); ++i) {
    // Get Truth Table for terminal
    const DNLTerminalFull& tTerm = naja::DNL::get()->getDNLTerminalFromID(inputs_[i]);
    // If direction is input, skip
    if (!tTerm.isTopPort() &&
        tTerm.getSnlBitTerm()->getDirection() != SNLBitTerm::Direction::Input) {
      const auto tt = SNLDesignModeling::getTruthTable(
          tTerm.getDNLInstance().getSNLInstance(),
          tTerm.getSnlBitTerm()->getOrderID());
      if (tt.isInitialized()) {
        if (tt.all0()) {
          termDNLID2varID_[inputs_[i]] = 0;
          continue;
        } else if (tt.all1()) {
          termDNLID2varID_[inputs_[i]] = 1;
          continue;
        }
      // LCOV_EXCL_START
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    termDNLID2varID_[inputs_[i]] =
        i + 2;  // +2 to avoid 0 and 1 which are reserved for constants
  }
  for (DNLID constIsoID : naja::DNL::get()->getDNLIsoDB().getConstant0Isos()) {
    const auto& constIso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(constIsoID);
    for (auto termID : constIso.getReaders()) {
        termDNLID2varID_[termID] = 0;
    }
    for (auto termID : constIso.getDrivers()) {
        termDNLID2varID_[termID] = 0;
    }
  }
  for (DNLID constIsoID : naja::DNL::get()->getDNLIsoDB().getConstant1Isos()) {
    const auto& constIso = naja::DNL::get()->getDNLIsoDB().getIsoFromIsoIDconst(constIsoID);
    for (auto termID : constIso.getReaders()) {
        termDNLID2varID_[termID] = 1;
    }
    for (auto termID : constIso.getDrivers()) {
        termDNLID2varID_[termID] = 1;
    }
  }
}

void BuildPrimaryOutputClauses::build() {
  naja::DNL::get();
  POs_.clear();
  POs_ = tbb::concurrent_vector<BoolExpr*>(outputs_.size());
  initVarNames();
  size_t processedOutputs = 0;
  tbb::task_arena arena(20);
  IsPIs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
  for (auto pi : inputs_) {
    if (pi >= IsPIs_.size()) {
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      std::string error = "PI " + std::to_string(pi) + " is out of range";
      throw std::runtime_error(error);
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    IsPIs_[pi] = true;
  }
  IsPOs_ = std::vector<bool>(naja::DNL::get()->getNBterms(), false);
  for (auto po : outputs_) {
    if (po >= IsPOs_.size()) {
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      std::string error = "PO " + std::to_string(po) + " is out of range";
      throw std::runtime_error(error);
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    IsPOs_[po] = true;
  }

  std::vector<size_t> representativeForOutput(outputs_.size());
  std::vector<size_t> representativeOutputs;
  representativeOutputs.reserve(outputs_.size());
  std::unordered_map<DNLID, size_t> firstOutputForIso;
  firstOutputForIso.reserve(outputs_.size());
  for (size_t i = 0; i < outputs_.size(); ++i) {
    const DNLID isoID = get()->getDNLTerminalFromID(outputs_[i]).getIsoID();
    if (isoID != DNLID_MAX) {
      auto [it, inserted] = firstOutputForIso.emplace(isoID, i);
      if (!inserted) {
        representativeForOutput[i] = it->second;
        continue;
      }
    }
    representativeForOutput[i] = i;
    representativeOutputs.emplace_back(i);
  }

  auto processOutput = [&](size_t i) {
    DNLID out = outputs_[i];
    #ifdef DEBUG_PRINTS
    printf("Procssing output %zu/%zu: %s\n", ++processedOutputs,
           outputs_.size(),
           get()
               ->getDNLTerminalFromID(out)
               .getSnlBitTerm()
               ->getName()
               .getString()
               .c_str());
    #endif

    DNLID isoID = get()->getDNLTerminalFromID(out).getIsoID();
    DEBUG_LOG("isoID: %zu\n", isoID);
    auto cachedIt = Tree2BoolExpr::iso2boolExpr_.find(isoID);
    if (isoID != DNLID_MAX &&
        cachedIt != Tree2BoolExpr::iso2boolExpr_.end() &&
        cachedIt->second != nullptr && cachedIt->second->isValid()) {
      POs_[i] = cachedIt->second;
      #ifdef DEBUG_CHECKS
      assert(POs_[i] != nullptr);
      #endif
      #ifdef DEBUG_PRINTS
      printf("Reusing iso output %s for output %s\n",
             POs_[i]->toString().c_str(),
             get()
                 ->getDNLTerminalFromID(out)
                 .getSnlBitTerm()
                 ->getName()
                 .getString()
                 .c_str());
      #endif
      return;
    }
    
    SNLLogicCloud cloud(out, IsPIs_, IsPOs_);
    #ifdef DEBUG_CHECKS
    auto startComp = std::chrono::steady_clock::now();
    #endif
    cloud.compute();
    #ifdef DEBUG_CHECKS
    auto endComp = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_comp = endComp - startComp;
    printf("Computation time for %lu: %f seconds\n", i, elapsed_seconds_comp.count());
    #endif
#ifdef DEBUG_CHECKS
    assert(cloud.SNLDesignModeling::getTruthTable().isInitialized());
#endif
    assert(POs_.size() - 1 >= i);
    #ifdef DEBUG_CHECKS
    auto startFin = std::chrono::steady_clock::now();    
    #endif
    if (cloud.getTruthTable().isValid()) {
      cloud.getTruthTable().finalize();
    }
    // LCOV_EXCL_START
    if (cloud.getTruthTable().isValid()) {
      auto hasCachedIsoExpression = [](DNLID termID) {
        if (termID == DNLID_MAX) {
          // LCOV_DISABLED_START
          return false;
          // LCOV_DISABLED_STOP
        }
        const auto& term = get()->getDNLTerminalFromID(termID);
        const DNLID termIsoID = term.getIsoID();
        if (termIsoID == DNLID_MAX) {
          // LCOV_DISABLED_START
          return false;
          // LCOV_DISABLED_STOP
        }
        const auto cached = Tree2BoolExpr::iso2boolExpr_.find(termIsoID);
        return cached != Tree2BoolExpr::iso2boolExpr_.end() &&
               cached->second != nullptr && cached->second->isValid();
      };
      DNLID unmappedInput = DNLID_MAX;
      for (const auto inputTermID : cloud.getInputs()) {
        if (inputTermID == DNLID_MAX) {
          // LCOV_DISABLED_START
          continue;  // LCOV_EXCL_LINE
          // LCOV_DISABLED_STOP
        }
        if (inputTermID >= termDNLID2varID_.size() ||
            termDNLID2varID_[inputTermID] == static_cast<size_t>(-1)) {
          if (hasCachedIsoExpression(inputTermID)) {
            continue;
          }
          // LCOV_DISABLED_START
          unmappedInput = inputTermID;
          break;
          // LCOV_DISABLED_STOP
        }
      }
      if (unmappedInput != DNLID_MAX) {
        // LCOV_DISABLED_START
        POs_[i] = BoolExpr::createInvalid();
        std::ostringstream detail;
        detail << "encountered internal frontier term "
               << unmappedInput
               << " that was not collected as a primary input";
               // LCOV_DISABLED_STOP
        {
          // LCOV_DISABLED_START
          std::lock_guard<std::mutex> lock(skippedOutputsMutex_);
          skippedOutputs_[out] = makeSkippedOutputInfo(
              SkippedOutputReason::NoDriver, detail.str());
        }
        reportSkippedPO(
            get(),
            get()->getDNLTerminalFromID(out),
            detail.str().c_str(),
            kSkippedNoDriverPOReport);
      }
      // LCOV_DISABLED_STOP
    }
    // LCOV_EXCL_STOP
    #ifdef DEBUG_CHECKS
    auto endFin = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_fin = endFin - startFin;
    printf("Finalization time for %lu: %f seconds\n", i, elapsed_seconds_fin.count());
    #endif
    #ifdef DEBUG_CHECKS
    auto startConv = std::chrono::steady_clock::now();
    #endif
    if (cloud.getTruthTable().isValid() && POs_[i] == nullptr) {
      POs_[i] = Tree2BoolExpr::convert(cloud.getTruthTable(), termDNLID2varID_);
    } else if (POs_[i] == nullptr) {
      POs_[i] = BoolExpr::createInvalid();
      SkippedOutputReason skipReason = SkippedOutputReason::None;
      switch (cloud.getSkipReason()) {
        case SNLLogicCloud::SkipReason::NoDriver:
          skipReason = SkippedOutputReason::NoDriver;
          break;
        case SNLLogicCloud::SkipReason::MultiDriver:
          skipReason = SkippedOutputReason::MultiDriver;
          break;
        case SNLLogicCloud::SkipReason::LogicalLoop:
          skipReason = SkippedOutputReason::LogicalLoop;
          break;
        // LCOV_EXCL_START
        case SNLLogicCloud::SkipReason::None:  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        default:
          // LCOV_EXCL_START
          break;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
      }
      if (skipReason != SkippedOutputReason::None) {
        std::lock_guard<std::mutex> lock(skippedOutputsMutex_);
        skippedOutputs_[out] = makeSkippedOutputInfo(
            skipReason, cloud.getSkipReasonText());
      }
    }
    #ifdef DEBUG_CHECKS
    auto endConv = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_conv = endConv - startConv;
    printf("Conversion time for %lu: %f seconds\n", i, elapsed_seconds_conv.count());
    #endif
    if (isoID != DNLID_MAX &&
        POs_[i] != nullptr && POs_[i]->isValid()) {
      // Publish only fully-built expressions; do not expose a placeholder entry.
      auto insertResult = Tree2BoolExpr::iso2boolExpr_.insert({isoID, POs_[i]});
      if (!insertResult.second) {
        POs_[i] = insertResult.first->second;
      }
    }
    cloud.destroy();
  };
  struct ScopedIsoCacheReset {
    ScopedIsoCacheReset() {
      Tree2BoolExpr::iso2boolExpr_.clear();
    }
    ~ScopedIsoCacheReset() {
      Tree2BoolExpr::iso2boolExpr_.clear();
    }
  } scopedIsoCacheReset;
  if (getenv("KEPLER_NO_MT")) {
    for (size_t i : representativeOutputs) {
      processOutput(i);
    }
  } else {
    // compute grain safely
    size_t n = representativeOutputs.size();
    size_t default_grain = 1000;
    size_t computed = (n >= 1000) ? (n / 1000) : 1; // never zero
    size_t grain = std::max<size_t>(computed, default_grain); // or clamp as you prefer

    tbb::parallel_for(
      tbb::blocked_range<DNLID>(0, n, grain),
      [&](const tbb::blocked_range<DNLID>& r) {
        for (DNLID i = r.begin(); i < r.end(); ++i) {
          processOutput(representativeOutputs[i]);
        }
      },
      tbb::static_partitioner()
    );
  }
  for (size_t i = 0; i < outputs_.size(); ++i) {
    const size_t representative = representativeForOutput[i];
    if (representative == i) {
      continue;
    }
    POs_[i] = POs_[representative];
    const auto skippedIt = skippedOutputs_.find(outputs_[representative]);
    if (skippedIt != skippedOutputs_.end()) {
      skippedOutputs_[outputs_[i]] = skippedIt->second;
    }
  }
  SNLLogicCloud::flushSkippedPOReports();
  if (!retainDnl_) {
    destroy();  // Clean up DNL instance unless the caller is batching builds.
  }
}

void BuildPrimaryOutputClauses::setInputs2InputsIDs() {
  inputs2inputsIDs_.clear();
  for (const auto& input : inputs_) {
    if (get()->getDNLTerminalFromID(input).isNull()) {
      // LCOV_EXCL_START
      throw std::runtime_error("Input terminal is null");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    const DNLInstanceFull& currentInstance =
        get()->getDNLTerminalFromID(input).getDNLInstance();
    PathKey& pair = inputs2inputsIDs_[input];
    pair.first = getPathNameIDs(currentInstance);
    pair.first.emplace_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getName().getID());
    pair.second.emplace_back(
        get()->getDNLTerminalFromID(input).getSnlBitTerm()->getBit());
  }
}

void BuildPrimaryOutputClauses::setOutputs2OutputsIDs() {
  outputs2outputsIDs_.clear();
  for (const auto& output : outputs_) {
    const DNLInstanceFull& currentInstance =
        get()->getDNLTerminalFromID(output).getDNLInstance();
    PathKey& pair = outputs2outputsIDs_[output];
    pair.first = getPathNameIDs(currentInstance);
    pair.first.emplace_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getName().getID());
    pair.second.emplace_back(
        get()->getDNLTerminalFromID(output).getSnlBitTerm()->getBit());
  }
}
