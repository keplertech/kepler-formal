// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "strategy/SequentialEquivalenceStrategy.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "DNL.h"
#include "NLUniverse.h"
#include "SNLPath.h"
#include "common/BoolExprUtils.h"
#include "common/AlignedSignals.h"
#include "common/SecDiag.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "imc/IMCEngine.h"
#include "kinduction/KInductionEngine.h"
#include "model/SequentialDesignModel.h"
#include "pdr/PDREngine.h"
#include "strategy/ReachableStateInvariant.h"
#include "strategy/StructuralStateInvariant.h"

namespace KEPLER_FORMAL::SEC {

// Overall SEC strategy pipeline:
// 1. Extract both designs into the normalized sequential model used by SEC.
// 2. Align environment inputs and observed outputs by stable external names.
// 3. Infer internal state correspondences structurally, not by register names.
// 4. Build reset/init reachable-state strengthening for startup anchoring.
// 5. Remap both designs into one shared SAT symbol space.
// 6. Build the checked SEC property and the stronger proof invariant.
// 7. Hand the combined transition system to the selected top-level engine and
//    translate its result back into user-facing SEC diagnostics.

namespace {

std::string joinReasons(const std::vector<std::string>& reasons) {
  std::ostringstream oss;
  for (size_t i = 0; i < reasons.size(); ++i) {
    if (i) {
      oss << " | ";
    }
    oss << reasons[i];
  }
  return oss.str();
}

std::string describeMismatchedNames(const std::vector<std::string>& lhs,
                                    const std::vector<std::string>& rhs,
                                    const char* label) {
  std::ostringstream oss;
  oss << "Mismatched " << label << " sets";
  if (!lhs.empty()) {
    oss << " lhs=[";
    for (size_t i = 0; i < lhs.size(); ++i) {
      if (i) {
        oss << ", ";
      }
      oss << lhs[i];
    }
    oss << "]";
  }
  if (!rhs.empty()) {
    oss << " rhs=[";
    for (size_t i = 0; i < rhs.size(); ++i) {
      if (i) {
        oss << ", ";
      }
      oss << rhs[i];
    }
    oss << "]";
  }
  return oss.str();
}

std::string formatBoolValue(bool value) {
  return value ? "1" : "0";
}

std::string normalizeSignalBaseName(const std::string& name) {
  std::string base = name;
  const auto bracket = base.find('[');
  if (bracket != std::string::npos) {
    base = base.substr(0, bracket);
  }
  std::transform(base.begin(), base.end(), base.begin(), [](unsigned char ch) {
    return static_cast<char>(std::toupper(ch));
  });
  return base;
}

std::optional<bool> getResetAssertionValue(const std::string& displayName) {
  const std::string normalized = normalizeSignalBaseName(displayName);
  if (normalized == "RESET" || normalized == "RST") {
    return true;
  }
  if (normalized == "RESET_N" || normalized == "RESETN" ||
      normalized == "RST_N" || normalized == "RSTN") {
    return false;
  }
  return std::nullopt;
}

SignalKey getTerminalPathKey(const naja::DNL::DNLTerminalFull& terminal) {
  SignalKey key;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  key.first.reserve(pathNames.size() + 1);
  for (const auto& name : pathNames) {
    key.first.push_back(stableSignalKeyNameID(name.getString()));
  }
  key.first.push_back(
      stableSignalKeyNameID(terminal.getSnlBitTerm()->getName().getString()));
  key.second.push_back(
      static_cast<naja::NL::NLID::DesignObjectID>(terminal.getSnlBitTerm()->getBit()));
  return key;
}

std::string getTerminalDisplayName(const naja::DNL::DNLTerminalFull& terminal) {
  std::ostringstream oss;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  for (const auto& name : pathNames) {
    oss << name.getString() << ".";
  }
  oss << terminal.getSnlBitTerm()->getName().getString() << "["
      << terminal.getSnlBitTerm()->getBit() << "]";
  return oss.str();
}

std::string formatStringList(const std::vector<std::string>& values, size_t limit) {
  if (values.empty()) {
    return "<none>";
  }

  std::ostringstream oss;
  const size_t printed = std::min(values.size(), limit);
  for (size_t i = 0; i < printed; ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  if (values.size() > printed) {
    oss << ", ... +" << (values.size() - printed) << " more";
  }
  return oss.str();
}

std::vector<std::string> setDifference(const std::set<std::string>& lhs,
                                       const std::set<std::string>& rhs) {
  std::vector<std::string> diff;
  std::set_difference(
      lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(diff));
  return diff;
}

std::string describeConnectivitySkipOrigin(ConnectivitySkipOrigin origin) {
  switch (origin) {
    case ConnectivitySkipOrigin::NoDriver:
      return "no-driver";
    case ConnectivitySkipOrigin::MultiDriver:
      return "multi-driver";
  }
  return "connectivity";  // LCOV_EXCL_LINE
}

std::string describeConnectivitySkipInfo(const ConnectivitySkipInfo& info) {
  std::ostringstream oss;
  oss << describeConnectivitySkipOrigin(info.origin) << " connectivity: "
      << info.detail;
  return oss.str();
}

struct OutputCoverageSelection {
  AlignedSignals checkedOutputs;
  std::vector<std::string> skippedOutputs;
  size_t totalOutputs = 0;
};

struct ScopedDnlContext {
  explicit ScopedDnlContext(naja::NL::SNLDesign* top)
      : universe_(naja::NL::NLUniverse::get()),
        previousTop_(universe_ ? universe_->getTopDesign() : nullptr) {
    if (universe_ == nullptr) {
      throw std::runtime_error("NLUniverse not created for SEC cone tracing");
    }

    naja::DNL::destroy();
    universe_->setTopDesign(top);
    dnl_ = naja::DNL::get();
  }

  ~ScopedDnlContext() {
    naja::DNL::destroy();
    if (universe_ != nullptr && previousTop_ != nullptr) {
      universe_->setTopDesign(previousTop_);
    }
  }

  naja::DNL::DNLFull* dnl() const {
    return dnl_;
  }

 private:
  naja::NL::NLUniverse* universe_ = nullptr;
  naja::NL::SNLDesign* previousTop_ = nullptr;
  naja::DNL::DNLFull* dnl_ = nullptr;
};

std::optional<naja::DNL::DNLID> findTermByDisplayName(
    naja::DNL::DNLFull* dnl,
    const std::string& signalName) {
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    if (getTerminalDisplayName(term) == signalName) {
      return termID;
    }
  }
  return std::nullopt;
}

std::optional<naja::DNL::DNLID> findTermByKey(naja::DNL::DNLFull* dnl,
                                              const SignalKey& key) {
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    if (getTerminalPathKey(term) == key) {
      return termID;
    }
  }
  return std::nullopt;
}

std::vector<naja::DNL::DNLID> resolveTermsByKey(
    naja::DNL::DNLFull* dnl,
    const std::vector<SignalKey>& keys) {
  std::vector<naja::DNL::DNLID> resolved;
  resolved.reserve(keys.size());
  for (const auto& key : keys) {
    if (auto termID = findTermByKey(dnl, key); termID.has_value()) {
      resolved.push_back(*termID);
    }
  }
  return resolved;
}

std::string formatConeTerm(naja::DNL::DNLFull* dnl, naja::DNL::DNLID termID) {
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.isNull()) {
    return "<null>";
  }
  if (term.getIsoID() == naja::DNL::DNLID_MAX) {
    return getTerminalDisplayName(term);
  }

  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  if (iso.isConstant0()) {
    return "Constant 0";
  }
  if (iso.isConstant1()) {
    return "Constant 1";
  }
  return getTerminalDisplayName(term);
}

struct ConeTrace {
  std::vector<std::vector<std::string>> levels;
  std::set<std::string> allTerms;
};

ConeTrace buildConeTrace(naja::DNL::DNLFull* dnl,
                         naja::DNL::DNLID seedTermID,
                         const std::vector<naja::DNL::DNLID>& environmentInputs) {
  ConeTrace trace;
  std::vector<bool> isEnvironmentInput(dnl->getDNLTerms().size(), false);
  for (const auto termID : environmentInputs) {
    if (termID < isEnvironmentInput.size()) {
      isEnvironmentInput[termID] = true;
    }
  }

  const auto seedIsoID = dnl->getDNLTerminalFromID(seedTermID).getIsoID();
  if (seedIsoID == naja::DNL::DNLID_MAX) {
    return trace;
  }

  std::vector<naja::DNL::DNLID> currentIsos = {seedIsoID};
  std::unordered_set<naja::DNL::DNLID> visitedIsos;
  while (!currentIsos.empty()) {
    std::set<std::string> levelTerms;
    std::vector<naja::DNL::DNLID> nextIsos;

    for (const auto isoID : currentIsos) {
      if (isoID == naja::DNL::DNLID_MAX || !visitedIsos.insert(isoID).second) {
        continue;
      }

      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant0()) {
        levelTerms.insert("Constant 0");
        continue;
      }
      if (iso.isConstant1()) {
        levelTerms.insert("Constant 1");
        continue;
      }

      for (const auto driver : iso.getDrivers()) {
        if (driver == naja::DNL::DNLID_MAX) {
          continue;
        }

        const auto& driverTerm = dnl->getDNLTerminalFromID(driver);
        if (driverTerm.isNull()) {
          continue;
        }

        levelTerms.insert(formatConeTerm(dnl, driver));
        if (driver < isEnvironmentInput.size() && isEnvironmentInput[driver]) {
          continue;
        }

        const auto& inst = driverTerm.getDNLInstance();
        for (naja::DNL::DNLID termID = inst.getTermIndexes().first;
             termID != naja::DNL::DNLID_MAX && termID <= inst.getTermIndexes().second;
             ++termID) {
          const auto& term = dnl->getDNLTerminalFromID(termID);
          if (term.isNull()) {
            continue;
          }
          if (term.getSnlBitTerm()->getDirection() ==
              naja::NL::SNLBitTerm::Direction::Output) {
            continue;
          }
          if (term.getIsoID() != naja::DNL::DNLID_MAX) {
            nextIsos.push_back(term.getIsoID());
          }
        }
      }
    }

    if (!levelTerms.empty()) {
      std::vector<std::string> orderedTerms(levelTerms.begin(), levelTerms.end());
      trace.allTerms.insert(orderedTerms.begin(), orderedTerms.end());
      trace.levels.push_back(std::move(orderedTerms));
    }

    std::sort(nextIsos.begin(), nextIsos.end());
    nextIsos.erase(std::unique(nextIsos.begin(), nextIsos.end()), nextIsos.end());
    currentIsos = std::move(nextIsos);
  }

  return trace;
}

std::string formatConeLevels(const ConeTrace& trace) {
  constexpr size_t kMaxLevels = 12;
  constexpr size_t kMaxTermsPerLevel = 12;

  if (trace.levels.empty()) {
    return "    <no traced cone terms>\n";
  }

  std::ostringstream oss;
  const size_t printedLevels = std::min(trace.levels.size(), kMaxLevels);
  for (size_t level = 0; level < printedLevels; ++level) {
    oss << "    step " << level << ": "
        << formatStringList(trace.levels[level], kMaxTermsPerLevel) << "\n";
  }
  if (trace.levels.size() > printedLevels) {
    oss << "    ... +" << (trace.levels.size() - printedLevels)
        << " more trace steps\n";
  }
  return oss.str();
}

struct ConeDiffReport {
  ConeTrace trace;
  std::string error;
};

ConeDiffReport buildConeDiffReport(naja::NL::SNLDesign* top,
                                   const std::string& differenceSignal,
                                   const std::vector<SignalKey>& environmentInputs) {
  ConeDiffReport report;
  ScopedDnlContext dnlContext(top);
  auto* dnl = dnlContext.dnl();

  const auto seedTermID = findTermByDisplayName(dnl, differenceSignal);
  if (!seedTermID.has_value()) {
    report.error =
        "could not resolve the differing SEC signal back into the DNL";
    return report;
  }

  report.trace = buildConeTrace(
      dnl, *seedTermID, resolveTermsByKey(dnl, environmentInputs));
  return report;
}

std::string formatConeTraceback(const KInductionResult::CounterexampleWitness& witness,
                                const SequentialDesignModel& model0,
                                const SequentialDesignModel& model1,
                                naja::NL::SNLDesign* top0,
                                naja::NL::SNLDesign* top1) {
  if (witness.outputMismatches.empty()) {
    return "";
  }
  const auto& differencePoint = witness.outputMismatches.front();

  std::ostringstream oss;
  oss << "Traceback for first differing point `" << differencePoint.signal
      << "` at cycle " << witness.badFrame << ":\n";

  try {
    const auto report0 = buildConeDiffReport(
        top0, differencePoint.signal, model0.environmentInputs);
    const auto report1 = buildConeDiffReport(
        top1, differencePoint.signal, model1.environmentInputs);

    if (!report0.error.empty() || !report1.error.empty()) {
      oss << "  Cone traceback unavailable: ";
      if (!report0.error.empty()) {
        oss << "design0 " << report0.error;
      }
      if (!report0.error.empty() && !report1.error.empty()) {
        oss << "; ";
      }
      if (!report1.error.empty()) {
        oss << "design1 " << report1.error;
      }
      oss << "\n";
      return oss.str();
    }

    oss << "  design0 cone to environment inputs:\n"
        << formatConeLevels(report0.trace);
    oss << "  design1 cone to environment inputs:\n"
        << formatConeLevels(report1.trace);

    constexpr size_t kMaxDiffTerms = 20;
    const auto onlyInDesign0 =
        setDifference(report0.trace.allTerms, report1.trace.allTerms);
    const auto onlyInDesign1 =
        setDifference(report1.trace.allTerms, report0.trace.allTerms);
    oss << "  cone terms only in design0: "
        << formatStringList(onlyInDesign0, kMaxDiffTerms) << "\n";
    oss << "  cone terms only in design1: "
        << formatStringList(onlyInDesign1, kMaxDiffTerms) << "\n";
  } catch (const std::exception& e) {
    oss << "  Cone traceback unavailable: " << e.what() << "\n";
  }

  return oss.str();
}

std::string formatCounterexampleWitness(const KInductionResult& result,
                                        const SequentialDesignModel& model0,
                                        const SequentialDesignModel& model1,
                                        naja::NL::SNLDesign* top0,
                                        naja::NL::SNLDesign* top1) {
  if (!result.witness.has_value()) {
    return "";
  }

  const auto& witness = *result.witness;
  std::ostringstream oss;
  oss << "Counterexample reaches the first bad frame at cycle "
      << witness.badFrame << ".\n";

  if (witness.inputTrace.empty()) {
    oss << "Input trace: <none>\n";
  } else {
    oss << "Input trace:\n";
    for (const auto& frame : witness.inputTrace) {
      oss << "  cycle " << frame.frame << ": ";
      if (frame.assignments.empty()) {
        oss << "<no environment inputs>";
      } else {
        for (size_t i = 0; i < frame.assignments.size(); ++i) {
          if (i) {
            oss << ", ";
          }
          oss << frame.assignments[i].signal << "="
              << formatBoolValue(frame.assignments[i].value);
        }
      }
      oss << "\n";
    }
  }

  if (!witness.outputMismatches.empty()) {
    oss << "Observed output mismatches at cycle " << witness.badFrame << ":\n";
    for (const auto& mismatch : witness.outputMismatches) {
      oss << "  " << mismatch.signal << ": design0="
          << formatBoolValue(mismatch.design0Value)
          << ", design1=" << formatBoolValue(mismatch.design1Value) << "\n";
    }
  }

  oss << formatConeTraceback(witness, model0, model1, top0, top1);

  return oss.str();
}

std::map<SignalKey, std::string, SignalKeyLess> buildKeyToNameMap(
    const std::vector<SignalKey>& keys,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames,
    const char* label) {
  std::map<SignalKey, std::string, SignalKeyLess> byKey;
  for (const auto& key : keys) {
    const auto nameIt = displayNames.find(key);
    if (nameIt == displayNames.end()) {
      throw std::runtime_error(
          std::string("Missing display name for SEC ") + label);
    }
    const auto [_, inserted] = byKey.emplace(key, nameIt->second);
    if (!inserted) {
      throw std::runtime_error(
          std::string("Duplicate SEC ") + label + " key `" + signalKeyToString(key) + "`");
    }
  }
  return byKey;
}

std::vector<std::string> sortedNames(
    const std::map<SignalKey, std::string, SignalKeyLess>& byKey) {
  std::vector<std::string> names;
  names.reserve(byKey.size());
  for (const auto& [_, name] : byKey) {
    names.push_back(name);
  }
  return names;
}

AlignedSignals alignSignalsByName(
    const std::vector<SignalKey>& keys0,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames0,
    const std::vector<SignalKey>& keys1,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames1,
    const char* label) {
  // Inputs/outputs are part of the user-visible contract, so SEC requires an
  // exact stable-key match before any proof engine is allowed to run.
  const auto byKey0 = buildKeyToNameMap(keys0, displayNames0, label);
  const auto byKey1 = buildKeyToNameMap(keys1, displayNames1, label);
  if (byKey0.size() != byKey1.size()) {
    throw std::runtime_error(
        describeMismatchedNames(sortedNames(byKey0), sortedNames(byKey1), label));
  }

  auto it0 = byKey0.begin();
  auto it1 = byKey1.begin();
  for (; it0 != byKey0.end() && it1 != byKey1.end(); ++it0, ++it1) {
    if (it0->first != it1->first) {
      throw std::runtime_error(
          describeMismatchedNames(sortedNames(byKey0), sortedNames(byKey1), label));
    }
  }

  AlignedSignals aligned;
  aligned.names.reserve(byKey0.size());
  aligned.keys0.reserve(byKey0.size());
  aligned.keys1.reserve(byKey0.size());
  for (const auto& [key, displayName] : byKey0) {
    aligned.names.push_back(displayName);
    aligned.keys0.push_back(key);
    aligned.keys1.push_back(key);
  }
  return aligned;
}

OutputCoverageSelection selectCoveredObservedOutputs(
    const AlignedSignals& allObservedOutputs,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1) {
  // Connectivity skips are the only SEC skips we allow here. Unsupported
  // primitive semantics should already have stopped extraction earlier.
  OutputCoverageSelection selection;
  selection.totalOutputs = allObservedOutputs.names.size();
  selection.checkedOutputs.names.reserve(allObservedOutputs.names.size());
  selection.checkedOutputs.keys0.reserve(allObservedOutputs.names.size());
  selection.checkedOutputs.keys1.reserve(allObservedOutputs.names.size());

  for (size_t i = 0; i < allObservedOutputs.names.size(); ++i) {
    const auto& key0 = allObservedOutputs.keys0[i];
    const auto& key1 = allObservedOutputs.keys1[i];
    const auto& name = allObservedOutputs.names[i];

    const auto skip0 = model0.connectivitySkipInfoByKey.find(key0);
    const auto skip1 = model1.connectivitySkipInfoByKey.find(key1);
    if (skip0 != model0.connectivitySkipInfoByKey.end() ||
        skip1 != model1.connectivitySkipInfoByKey.end()) {
      std::vector<std::string> reasons;
      if (skip0 != model0.connectivitySkipInfoByKey.end()) {
        reasons.push_back(
            "design0 " + describeConnectivitySkipInfo(skip0->second));
      }
      if (skip1 != model1.connectivitySkipInfoByKey.end()) {
        reasons.push_back(
            "design1 " + describeConnectivitySkipInfo(skip1->second));
      }
      selection.skippedOutputs.push_back(
          name + ": " + joinReasons(reasons));
      continue;
    }

    if (model0.observedOutputExprByKey.find(key0) ==
            model0.observedOutputExprByKey.end() ||
        model1.observedOutputExprByKey.find(key1) ==
            model1.observedOutputExprByKey.end()) {
      throw std::runtime_error(
          "Missing observed output expression for aligned SEC output `" +
          name + "`");
    }

    selection.checkedOutputs.names.push_back(name);
    selection.checkedOutputs.keys0.push_back(key0);
    selection.checkedOutputs.keys1.push_back(key1);
  }

  return selection;
}

SequentialEquivalenceResult makeSecResult(
    SequentialEquivalenceStatus status,
    size_t bound,
    std::string reason,
    const OutputCoverageSelection& coverage,
    std::vector<std::string> abstractedSequentialBoundaries = {}) {
  SequentialEquivalenceResult result;
  result.status = status;
  result.bound = bound;
  result.reason = std::move(reason);
  result.coveredOutputs = coverage.checkedOutputs.names.size();
  result.totalOutputs = coverage.totalOutputs;
  result.skippedObservedOutputs = coverage.skippedOutputs;
  result.abstractedSequentialBoundaries =
      std::move(abstractedSequentialBoundaries);
  return result;
}

template <typename MapT>
void assignSymbols(const std::vector<SignalKey>& keys,
                   MapT& output,
                   std::vector<size_t>& allSymbols,
                   size_t& nextSymbol) {
  for (const auto& key : keys) {
    output.emplace(key, nextSymbol);
    allSymbols.push_back(nextSymbol);
    ++nextSymbol;
  }
}

std::unordered_map<size_t, size_t> buildLocalToCombinedMap(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& inputSymbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& stateSymbols) {
  // Each design is extracted independently, so its BoolExpr variables must be
  // rewritten into a single shared symbol space before we can compare them.
  std::unordered_map<size_t, size_t> localToCombined;
  for (const auto& [key, localVar] : model.inputVarByKey) {
    if (auto it = inputSymbols.find(key); it != inputSymbols.end()) {
      localToCombined.emplace(localVar, it->second);
      continue;
    }
    if (auto it = stateSymbols.find(key); it != stateSymbols.end()) {
      localToCombined.emplace(localVar, it->second);
    }
  }
  return localToCombined;
}

}  // namespace

SequentialEquivalenceStrategy::SequentialEquivalenceStrategy(
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    KEPLER_FORMAL::Config::SolverType solverType,
    SecEngine secEngine)
    : top0_(top0), top1_(top1), solverType_(solverType), secEngine_(secEngine) {}

SequentialEquivalenceResult SequentialEquivalenceStrategy::run(size_t maxK) const {
  const bool secDiagEnabled = isSecDiagEnabled();
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: start run");
  }

  // Step 1: extract both tops into the same normalized SEC representation.
  SequentialDesignModel model0 = SequentialDesignModel::extract(top0_);
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: extracted design0");
  }
  std::vector<std::string> abstractedSequentialBoundaries;
  for (const auto& description : model0.abstractedSequentialBoundaries) {
    abstractedSequentialBoundaries.push_back("design0 " + description);
  }
  if (model0.hasUnsupportedFeatures()) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(model0.unsupportedReasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries);
  }
  SequentialDesignModel model1 = SequentialDesignModel::extract(top1_);
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: extracted design1");
  }
  abstractedSequentialBoundaries.reserve(
      model0.abstractedSequentialBoundaries.size() +
      model1.abstractedSequentialBoundaries.size());
  for (const auto& description : model1.abstractedSequentialBoundaries) {
    abstractedSequentialBoundaries.push_back("design1 " + description);
  }
  if (model1.hasUnsupportedFeatures()) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(model1.unsupportedReasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries);
  }

  if (model0.hasUnsupportedFeatures() || model1.hasUnsupportedFeatures()) {
    std::vector<std::string> reasons = model0.unsupportedReasons;
    reasons.insert(
        reasons.end(),
        model1.unsupportedReasons.begin(),
        model1.unsupportedReasons.end());
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(reasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries);
  }

  // Step 2: SEC only makes sense when the observable interfaces align by the
  // user-visible term names, not by parser-local object IDs.
  AlignedSignals alignedInputs;
  AlignedSignals alignedAllOutputs;
  AlignedSignals alignedOutputs;
  AlignedSignals inductiveStateEqualities;
  OutputCoverageSelection outputCoverage;
  try {
    if (secDiagEnabled) {
      emitSecDiag("SEC diag: aligning inputs/outputs");
    }
    alignedInputs = alignSignalsByName(
        model0.environmentInputs,
        model0.displayNameByKey,
        model1.environmentInputs,
        model1.displayNameByKey,
        "environment input");
    alignedAllOutputs = alignSignalsByName(
        model0.allObservedOutputs,
        model0.displayNameByKey,
        model1.allObservedOutputs,
        model1.displayNameByKey,
        "observed output");
    outputCoverage = selectCoveredObservedOutputs(
        alignedAllOutputs, model0, model1);
    alignedOutputs = outputCoverage.checkedOutputs;
    if (alignedOutputs.names.empty()) {
      return makeSecResult(
          SequentialEquivalenceStatus::Unsupported,
          0,
          "No aligned observed outputs remain after skipping cones with no-driver or "
          "multi-driver connectivity.",
          outputCoverage,
          abstractedSequentialBoundaries);
    }
    if (secDiagEnabled) {
      emitSecDiag(
          "SEC diag: checked_outputs=",
          alignedOutputs.names.size(),
          " total_outputs=",
          outputCoverage.totalOutputs,
          " skipped=",
          outputCoverage.skippedOutputs.size());
    }
    alignedOutputs = alignSignalsByName(
        model0.observedOutputs,
        model0.displayNameByKey,
        model1.observedOutputs,
        model1.displayNameByKey,
        "observed output");
    if (alignedOutputs.names.size() != outputCoverage.checkedOutputs.names.size()) {
      throw std::runtime_error(
          "Internal SEC error: checked observed outputs and extractor-visible observed "
          "outputs disagree after connectivity skipping");
    }
    if (secDiagEnabled) {
      emitSecDiag("SEC diag: inferring inductive state equalities");
    }
    // Internal-state correspondence must be inferred from the transition
    // structure itself. Matching register names is not strong enough for SEC:
    // equivalent RTL can rename or retime state freely.
    inductiveStateEqualities = inferStructurallyEquivalentStatePairs(
        model0, model1, alignedInputs);
    if (secDiagEnabled) {
      emitSecDiag("SEC diag: inferred inductive state equalities");
    }
  } catch (const std::exception& e) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        e.what(),
        outputCoverage,
        abstractedSequentialBoundaries);
  }

  if (secDiagEnabled) {
    emitSecDiag(
        "SEC diag: aligned_inputs=",
        alignedInputs.names.size(),
        " aligned_outputs=",
        alignedOutputs.names.size(),
        " inductive_state_equalities=",
        inductiveStateEqualities.names.size(),
        " state_bits0=",
        model0.stateBits.size(),
        " state_bits1=",
        model1.stateBits.size());
  }

  KInductionProblem problem;
  problem.environmentInputNames = alignedInputs.names;
  problem.observedOutputNames = alignedOutputs.names;

  // Step 3: create the shared symbol space used by the combined SAT problem.
  // Inputs are shared, while each design keeps its own private state vector.
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputSymbols0;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputSymbols1;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> state0Symbols;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> state1Symbols;
  size_t nextSymbol = 2;

  assignSymbols(model0.stateBits, state0Symbols, problem.allSymbols, nextSymbol);
  assignSymbols(model1.stateBits, state1Symbols, problem.allSymbols, nextSymbol);

  for (size_t i = 0; i < alignedInputs.names.size(); ++i) {
    const size_t symbol = nextSymbol++;
    inputSymbols0.emplace(alignedInputs.keys0[i], symbol);
    inputSymbols1.emplace(alignedInputs.keys1[i], symbol);
    problem.allSymbols.push_back(symbol);
    problem.inputSymbols.push_back(symbol);
    if (auto assertedValue = getResetAssertionValue(alignedInputs.names[i]);
        assertedValue.has_value()) {
      problem.resetBootstrapInputs.emplace_back(symbol, *assertedValue);
    }
  }
  for (const auto& key : model0.stateBits) {
    problem.state0Symbols.push_back(state0Symbols.at(key));
  }
  for (const auto& key : model1.stateBits) {
    problem.state1Symbols.push_back(state1Symbols.at(key));
  }
  for (const auto& relation : model0.complementedStateRelations) {
    if (state0Symbols.find(relation.primaryKey) != state0Symbols.end() &&
        state0Symbols.find(relation.complementedKey) != state0Symbols.end()) {
      problem.complementedStatePairs0.emplace_back(
          state0Symbols.at(relation.primaryKey),
          state0Symbols.at(relation.complementedKey));
    }
  }
  for (const auto& relation : model1.complementedStateRelations) {
    if (state1Symbols.find(relation.primaryKey) != state1Symbols.end() &&
        state1Symbols.find(relation.complementedKey) != state1Symbols.end()) {
      problem.complementedStatePairs1.emplace_back(
          state1Symbols.at(relation.primaryKey),
          state1Symbols.at(relation.complementedKey));
    }
  }

  const auto localToCombined0 =
      buildLocalToCombinedMap(model0, inputSymbols0, state0Symbols);
  const auto localToCombined1 =
      buildLocalToCombinedMap(model1, inputSymbols1, state1Symbols);

  std::unordered_map<BoolExpr*, BoolExpr*> remapMemo0;
  std::unordered_map<BoolExpr*, BoolExpr*> remapMemo1;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> remappedOutputs0;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> remappedOutputs1;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> remappedNext0;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> remappedNext1;

  // Step 4: rewrite both designs' formulas into that shared symbol space.
  for (size_t i = 0; i < alignedOutputs.names.size(); ++i) {
    const auto& key0 = alignedOutputs.keys0[i];
    const auto& key1 = alignedOutputs.keys1[i];
    remappedOutputs0.emplace(
        key0,
        remapBoolExprVariables(
            model0.observedOutputExprByKey.at(key0), localToCombined0, remapMemo0));
    remappedOutputs1.emplace(
        key1,
        remapBoolExprVariables(
            model1.observedOutputExprByKey.at(key1), localToCombined1, remapMemo1));
    problem.observedOutputExprs0.push_back(remappedOutputs0.at(key0));
    problem.observedOutputExprs1.push_back(remappedOutputs1.at(key1));
  }
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: remapped observed outputs");
  }

  for (const auto& key : model0.stateBits) {
    remappedNext0.emplace(
        key,
        remapBoolExprVariables(
            model0.nextStateExprByStateKey.at(key), localToCombined0, remapMemo0));
  }
  for (const auto& key : model1.stateBits) {
    remappedNext1.emplace(
        key,
        remapBoolExprVariables(
            model1.nextStateExprByStateKey.at(key), localToCombined1, remapMemo1));
  }
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: remapped next-state formulas");
  }

  // Step 5: if reset/init data is available, build the explicit frame-0 state
  // constraint before we hand the problem to k-induction.
  BoolExpr* initialCondition = BoolExpr::createTrue();
  auto addInitialStateAssignments =
      [&](const std::unordered_map<SignalKey, bool, SignalKeyHash>& initialValues,
          const std::unordered_map<SignalKey, size_t, SignalKeyHash>& stateSymbols) {
        for (const auto& [key, value] : initialValues) {
          const auto symbolIt = stateSymbols.find(key);
          if (symbolIt == stateSymbols.end()) {
            continue;
          }
          BoolExpr* literal = BoolExpr::Var(symbolIt->second);
          initialCondition = BoolExpr::And(
              initialCondition, value ? literal : BoolExpr::Not(literal));
          ++problem.initializedStateCount;
        }
      };
  addInitialStateAssignments(model0.initialStateValueByKey, state0Symbols);
  addInitialStateAssignments(model1.initialStateValueByKey, state1Symbols);
  problem.totalStateCount = problem.state0Symbols.size() + problem.state1Symbols.size();
  if (problem.hasExplicitInitialState()) {
    problem.initialCondition = BoolExpr::simplify(initialCondition);
  }
  const ReachableStateInvariant reachableInvariant = buildReachableStateInvariant(
      model0, model1, alignedInputs, inductiveStateEqualities, secDiagEnabled);
  // Reachable-state strengthening tells SEC which state equalities are already
  // justified at startup and which extra state values become known while reset
  // is held asserted for the bootstrap window.

  const AlignedSignals& initialStateCorrespondence =
      reachableInvariant.initialStateCorrespondence;
  for (size_t i = 0; i < initialStateCorrespondence.names.size(); ++i) {
    problem.initialStateEqualityPairs.emplace_back(
        state0Symbols.at(initialStateCorrespondence.keys0[i]),
        state1Symbols.at(initialStateCorrespondence.keys1[i]));
  }

  const AlignedSignals& anchoredStateEqualities =
      reachableInvariant.anchoredStateEqualities;
  for (const auto& [key, value] : reachableInvariant.bootstrapValues0) {
    if (state0Symbols.find(key) != state0Symbols.end()) {
      problem.bootstrapStateAssignments.emplace_back(state0Symbols.at(key), value);
    }
  }
  for (const auto& [key, value] : reachableInvariant.bootstrapValues1) {
    if (state1Symbols.find(key) != state1Symbols.end()) {
      problem.bootstrapStateAssignments.emplace_back(state1Symbols.at(key), value);
    }
  }
  problem.resetBootstrapCycles = reachableInvariant.bootstrapCycles;
  for (size_t i = 0; i < anchoredStateEqualities.names.size(); ++i) {
    problem.inductiveStateEqualityPairs.emplace_back(
        state0Symbols.at(anchoredStateEqualities.keys0[i]),
        state1Symbols.at(anchoredStateEqualities.keys1[i]));
    if (!problem.resetBootstrapInputs.empty()) {
      problem.bootstrapStateEqualityPairs.emplace_back(
          state0Symbols.at(anchoredStateEqualities.keys0[i]),
          state1Symbols.at(anchoredStateEqualities.keys1[i]));
    }
  }

  const auto [abstractOutputMap0, abstractOutputMap1] = buildAbstractTransitionMaps(
      model0, model1, alignedInputs, anchoredStateEqualities);
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: built abstract transition maps");
  }

  // Step 6: build the SEC proof obligations. The checked SEC property remains
  // pure observed-output equality, while induction uses a stronger invariant
  // made of the anchored internal-state equalities plus any remaining output
  // obligations not already implied by that state correspondence.
  BoolExpr* property = BoolExpr::createTrue();
  for (const auto& key : model0.stateBits) {
    problem.transitions0.emplace_back(state0Symbols.at(key), remappedNext0.at(key));
  }
  for (const auto& key : model1.stateBits) {
    problem.transitions1.emplace_back(state1Symbols.at(key), remappedNext1.at(key));
  }

  // Keep the SEC property honest: the base case always checks every observed
  // output, while the induction step uses a stronger invariant whose state
  // equalities are explicit BoolExpr clauses instead of out-of-band SAT
  // assumptions.
  BoolExpr* inductionProperty = BoolExpr::createTrue();
  for (size_t i = 0; i < anchoredStateEqualities.names.size(); ++i) {
    inductionProperty = BoolExpr::And(
        inductionProperty,
        makeEqualityExpr(
            BoolExpr::Var(state0Symbols.at(anchoredStateEqualities.keys0[i])),
            BoolExpr::Var(state1Symbols.at(anchoredStateEqualities.keys1[i]))));
  }
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    const auto outputEquality = makeEqualityExpr(
        problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]);
    property = BoolExpr::And(property, outputEquality);

    const auto& key0 = alignedOutputs.keys0[i];
    const auto& key1 = alignedOutputs.keys1[i];
    if (areEquivalentUnderAbstractMaps(
            model0.observedOutputExprByKey.at(key0),
            model1.observedOutputExprByKey.at(key1),
            abstractOutputMap0,
            abstractOutputMap1)) {
      continue;
    }
    inductionProperty = BoolExpr::And(inductionProperty, outputEquality);
  }

  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = BoolExpr::simplify(inductionProperty);
  problem.inductionBad = BoolExpr::simplify(BoolExpr::Not(problem.inductionProperty));

  problem.description = "SEC property with aligned observed outputs";
  if (secDiagEnabled) {
    emitSecDiag("SEC diag: built SEC and induction properties");
  }

  if (secDiagEnabled) {
    emitSecDiag(
        "SEC diag: property_is_true=",
        problem.property == BoolExpr::createTrue(),
        " induction_property_is_true=",
        problem.inductionProperty == BoolExpr::createTrue(),
        " bad_is_false=",
        problem.bad == BoolExpr::createFalse(),
        " induction_bad_is_false=",
        problem.inductionBad == BoolExpr::createFalse(),
        " reset_bootstrap_inputs=",
        problem.resetBootstrapInputs.size(),
        " bootstrap_cycles=",
        problem.resetBootstrapCycles,
        " bootstrap_equalities=",
        problem.bootstrapStateEqualityPairs.size(),
        " inductive_equalities=",
        problem.inductiveStateEqualityPairs.size());
  }

  // Step 7: hand the combined transition system to the selected proof engine.
  if (secDiagEnabled) {
    emitSecDiag(
        "SEC diag: entering ",
        secEngine_ == SecEngine::Pdr
            ? "pdr engine"
            : (secEngine_ == SecEngine::Imc
                   ? "imc engine"
                   : (secEngine_ == SecEngine::KInduction
                          ? "classic k-induction engine"
                          : "legacy engine")));
  }

  if (secEngine_ == SecEngine::Pdr) {
    // Reuse the existing base-case machinery so the new PDR path keeps the
    // same user-facing counterexample reporting for k=0 and combinational SEC.
    KInductionEngine baseline(problem, solverType_);
    const auto baselineResult = baseline.run(0);
    switch (baselineResult.status) {
      case KInductionStatus::Equivalent:
        return makeSecResult(
            SequentialEquivalenceStatus::Equivalent,
            baselineResult.bound,
            "",
            outputCoverage,
            abstractedSequentialBoundaries);
      case KInductionStatus::Different:
        return makeSecResult(
            SequentialEquivalenceStatus::Different,
            baselineResult.bound,
            formatCounterexampleWitness(baselineResult, model0, model1, top0_, top1_),
            outputCoverage,
            abstractedSequentialBoundaries);
      case KInductionStatus::Inconclusive:
      default:
        break;
    }

    PDREngine pdrEngine(problem, solverType_);
    const auto pdrResult = pdrEngine.run(maxK);
    switch (pdrResult.status) {
      case PDRStatus::Equivalent:
        return makeSecResult(
            SequentialEquivalenceStatus::Equivalent,
            pdrResult.bound,
            "",
            outputCoverage,
            abstractedSequentialBoundaries);
      case PDRStatus::Different: {
        // PDR proves reachability depth, but the legacy base-case solver already
        // knows how to reconstruct the concrete SEC witness and traceback.
        KInductionEngine witnessEngine(problem, solverType_);
        const auto witnessResult = witnessEngine.run(pdrResult.bound);
        const std::string details =
            witnessResult.status == KInductionStatus::Different
                ? formatCounterexampleWitness(witnessResult, model0, model1, top0_, top1_)
                : "PDR found a counterexample at k = " + std::to_string(pdrResult.bound);
        return makeSecResult(
            SequentialEquivalenceStatus::Different,
            pdrResult.bound,
            details,
            outputCoverage,
            abstractedSequentialBoundaries);
      }
      case PDRStatus::Inconclusive:
      default:
        return makeSecResult(
            SequentialEquivalenceStatus::Inconclusive,
            pdrResult.bound,
            "Reached max_k without a proof or counterexample",
            outputCoverage,
            abstractedSequentialBoundaries);
    }
  }

  if (secEngine_ == SecEngine::KInduction) {
    // Classic k-induction keeps the usual "search first, then prove" shape:
    // find bounded counterexamples, then try to close the proof by induction.
    KInductionEngine engine(problem, solverType_);
    const auto result = engine.run(maxK);
    switch (result.status) {
      case KInductionStatus::Equivalent:
        return makeSecResult(
            SequentialEquivalenceStatus::Equivalent,
            result.bound,
            "",
            outputCoverage,
            abstractedSequentialBoundaries);
      case KInductionStatus::Different:
        {
        return makeSecResult(
            SequentialEquivalenceStatus::Different,
            result.bound,
            result.witness.has_value()
                ? formatCounterexampleWitness(
                      result, model0, model1, top0_, top1_)
                : "Classic k-induction found a counterexample at k = " +
                      std::to_string(result.bound),
            outputCoverage,
            abstractedSequentialBoundaries);
        }
      case KInductionStatus::Inconclusive:
      default:
        return makeSecResult(
            SequentialEquivalenceStatus::Inconclusive,
            result.bound,
            "Reached max_k without a proof or counterexample",
            outputCoverage,
            abstractedSequentialBoundaries);
    }
  }

  if (secEngine_ == SecEngine::Imc) {
    // IMC still uses the shared bounded counterexample path, but its proof side
    // grows a reachable/interpolated frontier instead of using an induction
    // step over increasing simple paths.
    IMCEngine engine(problem, solverType_);
    const auto result = engine.run(maxK);
    switch (result.status) {
      case IMCStatus::Equivalent:
        return makeSecResult(
            SequentialEquivalenceStatus::Equivalent,
            result.bound,
            "",
            outputCoverage,
            abstractedSequentialBoundaries);
      case IMCStatus::Different:
        {
          const KInductionResult witnessResult{
              KInductionStatus::Different, result.bound, result.witness};
        return makeSecResult(
            SequentialEquivalenceStatus::Different,
            result.bound,
            result.witness.has_value()
                ? formatCounterexampleWitness(
                      witnessResult, model0, model1, top0_, top1_)
                : "IMC found a counterexample at k = " +
                      std::to_string(result.bound),
            outputCoverage,
            abstractedSequentialBoundaries);
        }
      case IMCStatus::Inconclusive:
      default:
        return makeSecResult(
            SequentialEquivalenceStatus::Inconclusive,
            result.bound,
            "Reached max_k without a proof or counterexample",
            outputCoverage,
            abstractedSequentialBoundaries);
    }
  }

  KInductionProblem legacyProblem = problem;
  ExactInterpolantSynthesizer interpolantSynthesizer(legacyProblem, solverType_);
  if (auto interpolant =
          interpolantSynthesizer.deriveOneStepReachableStateInvariant();
      interpolant.has_value()) {
    // Preserve the historical legacy behavior: it folds an exact small-state
    // interpolant into the induction invariant before the proof engines run.
    legacyProblem.inductionProperty = BoolExpr::simplify(
        BoolExpr::And(legacyProblem.inductionProperty, *interpolant));
    legacyProblem.inductionBad =
        BoolExpr::simplify(BoolExpr::Not(legacyProblem.inductionProperty));
  }

  KInductionEngine engine(legacyProblem, solverType_);
  const auto result = engine.run(maxK);
  switch (result.status) {
    case KInductionStatus::Equivalent:
      return makeSecResult(
          SequentialEquivalenceStatus::Equivalent,
          result.bound,
          "",
          outputCoverage,
          abstractedSequentialBoundaries);
    case KInductionStatus::Different:
      return makeSecResult(
          SequentialEquivalenceStatus::Different,
          result.bound,
          formatCounterexampleWitness(result, model0, model1, top0_, top1_),
          outputCoverage,
          abstractedSequentialBoundaries);
    case KInductionStatus::Inconclusive:
      break;
    default:
      break;
  }

  return makeSecResult(
      SequentialEquivalenceStatus::Inconclusive,
      result.bound,
      "Reached max_k without a proof or counterexample",
      outputCoverage,
      abstractedSequentialBoundaries);
}

namespace detail {

std::string formatStringListForTest(const std::vector<std::string>& values,
                                    size_t limit) {
  return formatStringList(values, limit);
}

std::string formatConeLevelsForTest(
    const std::vector<std::vector<std::string>>& levels) {
  ConeTrace trace;
  trace.levels = levels;
  for (const auto& level : levels) {
    trace.allTerms.insert(level.begin(), level.end());
  }
  return formatConeLevels(trace);
}

std::string formatCounterexampleWitnessForTest(
    const KInductionResult& result,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1) {
  return formatCounterexampleWitness(result, model0, model1, top0, top1);
}

AlignedSignals alignSignalsByNameForTest(
    const std::vector<SignalKey>& keys0,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames0,
    const std::vector<SignalKey>& keys1,
    const std::unordered_map<SignalKey, std::string, SignalKeyHash>& displayNames1,
    const char* label) {
  return alignSignalsByName(keys0, displayNames0, keys1, displayNames1, label);
}

SignalKey getTerminalPathKeyForTest(
    const naja::DNL::DNLTerminalFull& terminal) {
  return getTerminalPathKey(terminal);
}

std::optional<naja::DNL::DNLID> findTermByKeyForTest(
    naja::DNL::DNLFull* dnl,
    const SignalKey& key) {
  return findTermByKey(dnl, key);
}

}  // namespace detail

}  // namespace KEPLER_FORMAL::SEC
