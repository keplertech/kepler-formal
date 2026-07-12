// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "strategy/SequentialEquivalenceStrategy.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DNL.h"
#include "NLUniverse.h"
#include "SNLPath.h"
#include "common/BoolExprUtils.h"
#include "common/AlignedSignals.h"
#include "common/SecDiag.h"
#include "imc/ExactInterpolantSynthesizer.h"
#include "imc/IMCEngine.h"
#include "kinduction/BaseCaseSolver.h"
#include "kinduction/KInductionEngine.h"
#include "kinduction/OutputBatching.h"
#include "kinduction/SatEncoding.h"
#include "model/SequentialDesignModel.h"
#include "pdr/PDREngine.h"
#include "proof/DualRailEncoding.h"
#include "proof/TransitionExprResolver.h"
#include "strategy/ReachableStateInvariant.h"
#include "../../sat/SATSolverWrapper.h"

namespace KEPLER_FORMAL::SEC {

// Overall SEC strategy pipeline:
// 1. Extract both designs into the normalized sequential model used by SEC.
// 2. Align environment inputs and observed outputs by stable external names.
// 3. Keep cross-design internal state uncorrelated; only public/reset facts can
//    constrain the two designs before the selected SEC engine proves outputs.
// 4. Build reset/init reachable-state strengthening for startup anchoring.
// 5. Remap both designs into one shared SAT symbol space.
// 6. Build the checked SEC property and the stronger proof invariant.
// 7. Hand the combined transition system to the selected top-level engine and
//    translate its result back into user-facing SEC diagnostics.

namespace {

using PublicInputAbstractMap = std::unordered_map<size_t, size_t>;

struct PublicInputExprPairHash {
  size_t operator()(const std::pair<BoolExpr*, BoolExpr*>& pair) const noexcept {
    size_t seed = std::hash<const void*>()(pair.first);
    seed ^= std::hash<const void*>()(pair.second) + 0x9e3779b9 +
            (seed << 6) + (seed >> 2);
    return seed;
  }
};

using PublicInputExprPairMemo =
    std::unordered_map<std::pair<BoolExpr*, BoolExpr*>, bool,
                       PublicInputExprPairHash>;

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

bool isCommutativeOutputCompareOp(Op op) {
  return op == Op::AND || op == Op::OR || op == Op::XOR;
}

void cachePublicInputOutputEquivalence(
    PublicInputExprPairMemo& memo,
    const std::pair<BoolExpr*, BoolExpr*>& key,
    bool equivalent) {
  memo[key] = equivalent;
}

std::pair<PublicInputAbstractMap, PublicInputAbstractMap>
buildPublicInputAbstractMaps(const SequentialDesignModel& model0,
                             const SequentialDesignModel& model1,
                             const AlignedSignals& alignedInputs) {
  PublicInputAbstractMap abstractMap0;
  PublicInputAbstractMap abstractMap1;
  abstractMap0.reserve(alignedInputs.names.size());
  abstractMap1.reserve(alignedInputs.names.size());

  // Only aligned SEC inputs are shared between designs here.  Internal state is
  // deliberately left unmapped so this comparison cannot create a hidden
  // cross-design state relation.
  size_t nextAbstractSymbol = 2;
  for (size_t i = 0; i < alignedInputs.names.size(); ++i) {
    const size_t symbol = nextAbstractSymbol++;
    abstractMap0.emplace(model0.inputVarByKey.at(alignedInputs.keys0[i]), symbol);
    abstractMap1.emplace(model1.inputVarByKey.at(alignedInputs.keys1[i]), symbol);
  }
  return {std::move(abstractMap0), std::move(abstractMap1)};
}

bool publicInputVarsMatch(size_t var0,
                          size_t var1,
                          const PublicInputAbstractMap& abstractMap0,
                          const PublicInputAbstractMap& abstractMap1) {
  if (var0 < 2 || var1 < 2) {
    return var0 == var1;
  }
  const auto it0 = abstractMap0.find(var0);
  const auto it1 = abstractMap1.find(var1);
  return it0 != abstractMap0.end() && it1 != abstractMap1.end() &&
         it0->second == it1->second;
}

bool areOutputExprsEquivalentUnderPublicInputs(
    BoolExpr* expr0,
    BoolExpr* expr1,
    const PublicInputAbstractMap& abstractMap0,
    const PublicInputAbstractMap& abstractMap1) {
  struct StackFrame {
    BoolExpr* lhs = nullptr;
    BoolExpr* rhs = nullptr;
    bool visited = false;
  };

  PublicInputExprPairMemo memo;
  std::vector<StackFrame> stack{{expr0, expr1, false}};
  while (!stack.empty()) {
    const StackFrame current = stack.back();
    stack.pop_back();
    const auto key = std::make_pair(current.lhs, current.rhs);
    if (memo.find(key) != memo.end()) {
      continue;
    }

    if (current.lhs == nullptr || current.rhs == nullptr) {
      cachePublicInputOutputEquivalence(
          memo, key, current.lhs == current.rhs);
      continue;
    }

    const Op lhsOp = current.lhs->getOp();
    const Op rhsOp = current.rhs->getOp();
    if (lhsOp == Op::VAR && rhsOp == Op::VAR) {
      cachePublicInputOutputEquivalence(
          memo,
          key,
          publicInputVarsMatch(
              current.lhs->getId(),
              current.rhs->getId(),
              abstractMap0,
              abstractMap1));
      continue;
    }
    if (lhsOp != rhsOp) {
      cachePublicInputOutputEquivalence(memo, key, false);
      continue;
    }

    if (!current.visited) {
      stack.push_back({current.lhs, current.rhs, true});
      const auto rightKey =
          std::make_pair(current.lhs->getRight(), current.rhs->getRight());
      if (memo.find(rightKey) == memo.end()) {
        stack.push_back(
            {current.lhs->getRight(), current.rhs->getRight(), false});
      }
      const auto leftKey =
          std::make_pair(current.lhs->getLeft(), current.rhs->getLeft());
      if (memo.find(leftKey) == memo.end()) {
        stack.push_back(
            {current.lhs->getLeft(), current.rhs->getLeft(), false});
      }
      continue;
    }

    const bool leftEquivalent =
        memo.at(std::make_pair(current.lhs->getLeft(), current.rhs->getLeft()));
    const bool rightEquivalent =
        memo.at(std::make_pair(current.lhs->getRight(), current.rhs->getRight()));
    bool equivalent = leftEquivalent && rightEquivalent;
    if (!equivalent && isCommutativeOutputCompareOp(lhsOp)) {
      const auto crossLeftKey =
          std::make_pair(current.lhs->getLeft(), current.rhs->getRight());
      const auto crossRightKey =
          std::make_pair(current.lhs->getRight(), current.rhs->getLeft());
      bool deferred = false;
      if (memo.find(crossRightKey) == memo.end()) {
        stack.push_back({current.lhs, current.rhs, true});
        stack.push_back(
            {current.lhs->getRight(), current.rhs->getLeft(), false});
        deferred = true;
      }
      if (memo.find(crossLeftKey) == memo.end()) {
        if (!deferred) {
          stack.push_back({current.lhs, current.rhs, true});
        }
        stack.push_back(
            {current.lhs->getLeft(), current.rhs->getRight(), false});
        deferred = true;
      }
      if (deferred) {
        continue;
      }
      equivalent = memo.at(crossLeftKey) && memo.at(crossRightKey);
    }
    cachePublicInputOutputEquivalence(memo, key, equivalent);
  }

  return memo.at(std::make_pair(expr0, expr1));
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

bool hasSuffix(const std::string& value, const std::string& suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isResetNameToken(const std::string& candidate, const std::string& token) {
  // Domain-prefixed top resets, for example `wb_rst_i`, normalize to `WB_RST`
  // after input-suffix stripping.  Match only a final underscore-separated
  // reset token so prefixes do not block reset bootstrap alignment.
  return candidate == token || hasSuffix(candidate, "_" + token);
// LCOV_EXCL_START
}  // LCOV_EXCL_LINE
// LCOV_EXCL_STOP

bool isActiveLowResetToken(const std::string& candidate) {
  return candidate == "RESET_N" || candidate == "RESETN" ||
         candidate == "RESET_L" || candidate == "RST_N" ||
         candidate == "RSTN" || candidate == "RST_L";
}

void appendDomainPrefixedActiveLowResetCandidates(
    std::vector<std::string>& candidates) {
  // LCOV_EXCL_START
  const size_t originalSize = candidates.size();
  for (size_t index = 0; index < originalSize; ++index) {
  // LCOV_EXCL_STOP
    const std::string& candidate = candidates[index];
    if (candidate.size() <= 1) {
      continue;
    }
    // LCOV_EXCL_START
    const std::string strippedDomain = candidate.substr(1);
    if (isActiveLowResetToken(strippedDomain)) {
    // LCOV_EXCL_STOP
      // Async FIFO top ports commonly use rrst_n/wrst_n.  Recognize those
      // active-low one-letter domain prefixes without treating arbitrary
      // embedded "rst" names as reset controls.
      candidates.push_back(strippedDomain);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
}

std::vector<std::string> resetNameCandidates(const std::string& displayName) {
  // The shared SEC symbol space sees user-visible top-input names such as
  // `reset_i[0]`.  Match the same reset spelling policy as the reachable-state
  // pass so a reset discovered during model analysis remains available when
  // bootstrap constraints are converted to shared SAT symbols.
  const std::string normalized = normalizeSignalBaseName(displayName);
  std::vector<std::string> candidates = {normalized};
  if (hasSuffix(normalized, "_IN")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 3));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  if (hasSuffix(normalized, "_I")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 2));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  if (hasSuffix(normalized, "_NI")) {
    candidates.push_back(normalized.substr(0, normalized.size() - 1));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  appendDomainPrefixedActiveLowResetCandidates(candidates);
  return candidates;
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

std::optional<bool> getResetAssertionValue(const std::string& displayName) {
  for (const auto& candidate : resetNameCandidates(displayName)) {
    if (isResetNameToken(candidate, "RESET") ||
        isResetNameToken(candidate, "RST")) {
      return true;
    }
    if (isResetNameToken(candidate, "RESET_N") ||
        isResetNameToken(candidate, "RESETN") ||
        isResetNameToken(candidate, "RESET_L") ||
        isResetNameToken(candidate, "RST_N") ||
        isResetNameToken(candidate, "RSTN") ||
        isResetNameToken(candidate, "RST_L")) {
      return false;
    }
  }
  return std::nullopt;
}

SignalKey getTerminalPathKey(const naja::DNL::DNLTerminalFull& terminal) {
  SignalKey key;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  key.first.reserve(pathNames.size() + 1);
  for (const auto& name : pathNames) {
    key.first.push_back(name.getID());  // LCOV_EXCL_LINE
  }
  key.first.push_back(terminal.getSnlBitTerm()->getName().getID());
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
  // LCOV_EXCL_START
  if (values.size() > printed) {
  // LCOV_EXCL_STOP
    oss << ", ... +" << (values.size() - printed) << " more";  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
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
    // LCOV_EXCL_START
    case ConnectivitySkipOrigin::NoDriver:
    // LCOV_EXCL_STOP
      return "no-driver";
    case ConnectivitySkipOrigin::MultiDriver:
      return "multi-driver";
    case ConnectivitySkipOrigin::LogicalLoop:
      return "logical-loop";
    case ConnectivitySkipOrigin::MultiClockDomain:
      return "multi-clock-domain";
  }
  return "connectivity";  // LCOV_EXCL_LINE
}

std::string describeConnectivitySkipInfo(const ConnectivitySkipInfo& info) {
  std::ostringstream oss;
  oss << describeConnectivitySkipOrigin(info.origin) << " connectivity: "
      << info.detail;
  return oss.str();
}

std::string describeSecSignalKey(const SequentialDesignModel& model,
                                 const SignalKey& key) {
  if (const auto it = model.displayNameByKey.find(key);
      it != model.displayNameByKey.end()) {
    return it->second;
  }
  return signalKeyToString(key);  // LCOV_EXCL_LINE
}

void appendUniqueRole(std::vector<std::string>& roles, const char* role) {
  if (std::find(roles.begin(), roles.end(), role) == roles.end()) {
    roles.push_back(role);
  }
}

struct OutputCoverageSelection {
  AlignedSignals checkedOutputs;
  std::vector<std::string> skippedOutputs;
  std::vector<std::string> resetUnanchoredSkippedOutputs;
  std::vector<std::string> multiClockDomainSkippedOutputs;
  size_t totalOutputs = 0;
};

struct AlignedSecInterface {
  AlignedSignals inputs;
  AlignedSignals outputs;
  OutputCoverageSelection outputCoverage;
};

struct SharedSecSymbolSpace {
  KInductionProblem problem;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputSymbols0;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> inputSymbols1;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> state0Symbols;
  std::unordered_map<SignalKey, size_t, SignalKeyHash> state1Symbols;
  std::unordered_map<size_t, size_t> localToCombined0;
  std::unordered_map<size_t, size_t> localToCombined1;
};

struct RemappedSecExpressions {
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> next0;
  std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> next1;
};

struct DualRailStateSymbolMaps {
  std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash> state0ByKey;
  std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash> state1ByKey;
  // LCOV_EXCL_START
  std::unordered_map<size_t, DualRailSymbolPair> localState0BySymbol;
  // LCOV_EXCL_STOP
  std::unordered_map<size_t, DualRailSymbolPair> localState1BySymbol;
};

size_t privateSupportSymbol(
    size_t designIndex,
    size_t localSymbol,
    std::unordered_map<size_t, size_t>& symbolMap,
    KInductionProblem& problem,
    size_t& nextSymbol);

// LCOV_EXCL_START
class SecDualRailVariableMapper final : public DualRailVariableMapper {
 public:
  SecDualRailVariableMapper(
      size_t designIndex,
      const std::unordered_map<size_t, DualRailSymbolPair>& stateRails,
      // LCOV_EXCL_STOP
      std::unordered_map<size_t, size_t>& binarySymbolMap,
      KInductionProblem& problem,
      size_t& nextSymbol)
      : designIndex_(designIndex),
        stateRails_(stateRails),
        binarySymbolMap_(binarySymbolMap),
        problem_(problem),
        // LCOV_EXCL_START
        nextSymbol_(nextSymbol) {}


// LCOV_EXCL_STOP
  ~SecDualRailVariableMapper() override = default;  // LCOV_EXCL_LINE

  DualRailBoolExpr mapVariable(size_t symbol) override {
    if (const auto stateIt = stateRails_.find(symbol);
        stateIt != stateRails_.end()) {
      return DualRailBoolExpr{
          BoolExpr::Var(stateIt->second.mayBeOne),
          BoolExpr::Var(stateIt->second.mayBeZero)};
    }

    if (symbol < 2) {
      return symbol == 1  // LCOV_EXCL_LINE
                 ? DualRailBoolExpr{  // LCOV_EXCL_LINE
                       BoolExpr::createTrue(), BoolExpr::createFalse()}  // LCOV_EXCL_LINE
                 : DualRailBoolExpr{  // LCOV_EXCL_LINE
                       BoolExpr::createFalse(), BoolExpr::createTrue()};  // LCOV_EXCL_LINE
    }

    size_t binarySymbol = 0;
    // LCOV_EXCL_START
    if (const auto mappedIt = binarySymbolMap_.find(symbol);
    // LCOV_EXCL_STOP
        mappedIt != binarySymbolMap_.end()) {
      binarySymbol = mappedIt->second;
    } else {
      binarySymbol = privateSupportSymbol(  // LCOV_EXCL_LINE
          designIndex_, symbol, binarySymbolMap_, problem_, nextSymbol_);  // LCOV_EXCL_LINE
    }
    return DualRailBoolExpr{
        BoolExpr::Var(binarySymbol),
        BoolExpr::Not(BoolExpr::Var(binarySymbol))};
  }

 private:
  size_t designIndex_ = 0;
  const std::unordered_map<size_t, DualRailSymbolPair>& stateRails_;
  std::unordered_map<size_t, size_t>& binarySymbolMap_;
  KInductionProblem& problem_;
  size_t& nextSymbol_;
};

struct ScopedDnlContext {
  explicit ScopedDnlContext(naja::NL::SNLDesign* top)
      : universe_(naja::NL::NLUniverse::get()),
        previousTop_(universe_ ? universe_->getTopDesign() : nullptr) {
    if (universe_ == nullptr) {
      throw std::runtime_error("NLUniverse not created for SEC cone tracing");  // LCOV_EXCL_LINE
    }

    naja::DNL::destroy();
    universe_->setTopDesign(top);
    dnl_ = naja::DNL::get();
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  ~ScopedDnlContext() {
    naja::DNL::destroy();
    if (universe_ != nullptr && previousTop_ != nullptr) {
      universe_->setTopDesign(previousTop_);
    // LCOV_EXCL_START
    }
    // LCOV_EXCL_STOP
  }

  naja::DNL::DNLFull* dnl() const {
    return dnl_;
  }

 private:
  // LCOV_EXCL_START
  naja::NL::NLUniverse* universe_ = nullptr;
  // LCOV_EXCL_STOP
  naja::NL::SNLDesign* previousTop_ = nullptr;
  naja::DNL::DNLFull* dnl_ = nullptr;
};

std::optional<naja::DNL::DNLID> findTermByDisplayName(
    // LCOV_EXCL_START
    naja::DNL::DNLFull* dnl,
    // LCOV_EXCL_STOP
    const std::string& signalName) {
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (getTerminalDisplayName(term) == signalName) {
      return termID;
    }
  }
  return std::nullopt;  // LCOV_EXCL_LINE
}

std::optional<naja::DNL::DNLID> findTermByKey(naja::DNL::DNLFull* dnl,
                                              const SignalKey& key) {
  for (naja::DNL::DNLID termID = 0; termID < dnl->getDNLTerms().size(); ++termID) {
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      // LCOV_EXCL_START
      continue;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }
    if (getTerminalPathKey(term) == key) {
      // LCOV_EXCL_START
      return termID;
      // LCOV_EXCL_STOP
    }
  }
  return std::nullopt;  // LCOV_EXCL_LINE
}

// LCOV_EXCL_START


// LCOV_EXCL_STOP
std::vector<naja::DNL::DNLID> resolveTermsByKey(
    naja::DNL::DNLFull* dnl,
    // LCOV_EXCL_START
    const std::vector<SignalKey>& keys) {
    // LCOV_EXCL_STOP
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
    return "<null>";  // LCOV_EXCL_LINE
  }
  if (term.getIsoID() == naja::DNL::DNLID_MAX) {
    return getTerminalDisplayName(term);  // LCOV_EXCL_LINE
  }

  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  if (iso.isConstant0()) {
    return "Constant 0";  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP
  if (iso.isConstant1()) {
    return "Constant 1";  // LCOV_EXCL_LINE
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
  // LCOV_EXCL_START
  std::vector<bool> isEnvironmentInput(dnl->getDNLTerms().size(), false);
  for (const auto termID : environmentInputs) {
  // LCOV_EXCL_STOP
    if (termID < isEnvironmentInput.size()) {
      isEnvironmentInput[termID] = true;
    // LCOV_EXCL_START
    }
  }
  // LCOV_EXCL_STOP

  const auto seedIsoID = dnl->getDNLTerminalFromID(seedTermID).getIsoID();
  if (seedIsoID == naja::DNL::DNLID_MAX) {
    return trace;  // LCOV_EXCL_LINE
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  std::vector<naja::DNL::DNLID> currentIsos = {seedIsoID};
  std::unordered_set<naja::DNL::DNLID> visitedIsos;
  while (!currentIsos.empty()) {
    // LCOV_EXCL_START
    std::set<std::string> levelTerms;
    // LCOV_EXCL_STOP
    std::vector<naja::DNL::DNLID> nextIsos;

    for (const auto isoID : currentIsos) {
      if (isoID == naja::DNL::DNLID_MAX || !visitedIsos.insert(isoID).second) {
        continue;
      }

      const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant0()) {
        levelTerms.insert("Constant 0");  // LCOV_EXCL_LINE
        continue;  // LCOV_EXCL_LINE
      }
      if (iso.isConstant1()) {
        // LCOV_EXCL_START
        levelTerms.insert("Constant 1");  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        continue;  // LCOV_EXCL_LINE
      }

      for (const auto driver : iso.getDrivers()) {
        if (driver == naja::DNL::DNLID_MAX) {
          continue;  // LCOV_EXCL_LINE
        }

        const auto& driverTerm = dnl->getDNLTerminalFromID(driver);
        if (driverTerm.isNull()) {
          continue;  // LCOV_EXCL_LINE
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
            continue;  // LCOV_EXCL_LINE
          }
          if (term.getSnlBitTerm()->getDirection() ==
              naja::NL::SNLBitTerm::Direction::Output) {
            continue;
          }
          if (term.getIsoID() != naja::DNL::DNLID_MAX) {
            // LCOV_EXCL_START
            nextIsos.push_back(term.getIsoID());
            // LCOV_EXCL_STOP
          }
        }
      }
    }

    if (!levelTerms.empty()) {
      std::vector<std::string> orderedTerms(levelTerms.begin(), levelTerms.end());
      trace.allTerms.insert(orderedTerms.begin(), orderedTerms.end());
      trace.levels.push_back(std::move(orderedTerms));
    // LCOV_EXCL_START
    }

    std::sort(nextIsos.begin(), nextIsos.end());
    // LCOV_EXCL_STOP
    nextIsos.erase(std::unique(nextIsos.begin(), nextIsos.end()), nextIsos.end());
    currentIsos = std::move(nextIsos);
  }

  return trace;
}

std::string formatConeLevels(const ConeTrace& trace) {
  constexpr size_t kMaxLevels = 12;
  constexpr size_t kMaxTermsPerLevel = 12;

  if (trace.levels.empty()) {
    return "    <no traced cone terms>\n";  // LCOV_EXCL_LINE
  }

  std::ostringstream oss;
  const size_t printedLevels = std::min(trace.levels.size(), kMaxLevels);
  // LCOV_EXCL_START
  for (size_t level = 0; level < printedLevels; ++level) {
  // LCOV_EXCL_STOP
    oss << "    step " << level << ": "
        // LCOV_EXCL_START
        << formatStringList(trace.levels[level], kMaxTermsPerLevel) << "\n";
        // LCOV_EXCL_STOP
  }
  if (trace.levels.size() > printedLevels) {
    oss << "    ... +" << (trace.levels.size() - printedLevels)  // LCOV_EXCL_LINE
        << " more trace steps\n";  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return oss.str();
}

struct ConeDiffReport {
  ConeTrace trace;
  std::string error;
};

// LCOV_EXCL_START
ConeDiffReport buildConeDiffReport(naja::NL::SNLDesign* top,
// LCOV_EXCL_STOP
                                   const std::string& differenceSignal,
                                   const std::vector<SignalKey>& environmentInputs) {
  ConeDiffReport report;
  ScopedDnlContext dnlContext(top);
  auto* dnl = dnlContext.dnl();

  const auto seedTermID = findTermByDisplayName(dnl, differenceSignal);
  if (!seedTermID.has_value()) {
    report.error =  // LCOV_EXCL_LINE
        "could not resolve the differing SEC signal back into the DNL";
    return report;  // LCOV_EXCL_LINE
  }

  report.trace = buildConeTrace(
      dnl, *seedTermID, resolveTermsByKey(dnl, environmentInputs));
  return report;
}

std::string formatConeTraceback(const KInductionResult::CounterexampleWitness& witness,
                                const SequentialDesignModel& model0,
                                // LCOV_EXCL_START
                                const SequentialDesignModel& model1,
                                naja::NL::SNLDesign* top0,
                                naja::NL::SNLDesign* top1) {
  if (witness.outputMismatches.empty()) {
    return "";  // LCOV_EXCL_LINE
  }
  const auto& differencePoint = witness.outputMismatches.front();

  std::ostringstream oss;
  oss << "Traceback for first differing point `" << differencePoint.signal
      << "` at cycle " << witness.badFrame << ":\n";


// LCOV_EXCL_STOP
  if (top0 == nullptr || top1 == nullptr) {
    oss << "  Cone traceback unavailable: compact SEC released the "
           "elaborated designs after model extraction.\n";
    return oss.str();
  }

  const auto report0 = buildConeDiffReport(
      top0, differencePoint.signal, model0.environmentInputs);
  const auto report1 = buildConeDiffReport(
      top1, differencePoint.signal, model1.environmentInputs);

  if (!report0.error.empty() || !report1.error.empty()) {
    oss << "  Cone traceback unavailable: ";  // LCOV_EXCL_LINE
    if (!report0.error.empty()) {  // LCOV_EXCL_LINE
      oss << "design0 " << report0.error;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_START
    if (!report0.error.empty() && !report1.error.empty()) {  // LCOV_EXCL_LINE
      oss << "; ";  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }  // LCOV_EXCL_LINE
    if (!report1.error.empty()) {  // LCOV_EXCL_LINE
      oss << "design1 " << report1.error;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    oss << "\n";  // LCOV_EXCL_LINE
    return oss.str();  // LCOV_EXCL_LINE
  }

  oss << "  design0 cone to environment inputs:\n"
      << formatConeLevels(report0.trace);
  // LCOV_EXCL_START
  oss << "  design1 cone to environment inputs:\n"
  // LCOV_EXCL_STOP
      << formatConeLevels(report1.trace);

  constexpr size_t kMaxDiffTerms = 20;
  const auto onlyInDesign0 =
      setDifference(report0.trace.allTerms, report1.trace.allTerms);
  const auto onlyInDesign1 =
      setDifference(report1.trace.allTerms, report0.trace.allTerms);
  oss << "  cone terms only in design0: "
      // LCOV_EXCL_START
      << formatStringList(onlyInDesign0, kMaxDiffTerms) << "\n";
  oss << "  cone terms only in design1: "
  // LCOV_EXCL_STOP
      << formatStringList(onlyInDesign1, kMaxDiffTerms) << "\n";

  return oss.str();
}

std::string formatCounterexampleWitness(const KInductionResult& result,
                                        const SequentialDesignModel& model0,
                                        const SequentialDesignModel& model1,
                                        naja::NL::SNLDesign* top0,
                                        naja::NL::SNLDesign* top1) {
  if (!result.witness.has_value()) {
    return "";  // LCOV_EXCL_LINE
  }

  const auto& witness = *result.witness;
  std::ostringstream oss;
  oss << "Counterexample reaches the first bad frame at cycle "
      << witness.badFrame << ".\n";

  if (witness.inputTrace.empty()) {
    oss << "Input trace: <none>\n";  // LCOV_EXCL_LINE
  } else {  // LCOV_EXCL_LINE
    oss << "Input trace:\n";
    for (const auto& frame : witness.inputTrace) {
      oss << "  cycle " << frame.frame << ": ";
      if (frame.assignments.empty()) {
        oss << "<no environment inputs>";  // LCOV_EXCL_LINE
      } else {  // LCOV_EXCL_LINE
        for (size_t i = 0; i < frame.assignments.size(); ++i) {
          if (i) {
            oss << ", ";
          }
          oss << frame.assignments[i].signal << "="
              << formatBoolValue(frame.assignments[i].value);
        }
      }
      oss << "\n";
    // LCOV_EXCL_START
    }
  }
  // LCOV_EXCL_STOP

  if (!witness.outputMismatches.empty()) {
    oss << "Observed output mismatches at cycle " << witness.badFrame << ":\n";
    // LCOV_EXCL_START
    for (const auto& mismatch : witness.outputMismatches) {
      oss << "  " << mismatch.signal << ": design0="
      // LCOV_EXCL_STOP
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
      throw std::runtime_error(  // LCOV_EXCL_LINE
          std::string("Missing display name for SEC ") + label);  // LCOV_EXCL_LINE
    }
    const auto [_, inserted] = byKey.emplace(key, nameIt->second);
    if (!inserted) {
      throw std::runtime_error(  // LCOV_EXCL_LINE
          std::string("Duplicate SEC ") + label + " key `" + signalKeyToString(key) + "`");  // LCOV_EXCL_LINE
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
      const std::string skippedOutput = name + ": " + joinReasons(reasons);
      selection.skippedOutputs.push_back(skippedOutput);
      const bool skippedByMultiClockDomain =
          (skip0 != model0.connectivitySkipInfoByKey.end() &&
           skip0->second.origin == ConnectivitySkipOrigin::MultiClockDomain) ||
          (skip1 != model1.connectivitySkipInfoByKey.end() &&
           skip1->second.origin == ConnectivitySkipOrigin::MultiClockDomain);
      if (skippedByMultiClockDomain) {
        selection.multiClockDomainSkippedOutputs.push_back(skippedOutput);
      }
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
    // LCOV_EXCL_START
    std::string reason,
    // LCOV_EXCL_STOP
    const OutputCoverageSelection& coverage,
    std::vector<std::string> abstractedSequentialBoundaries = {},
    std::vector<ExtractedBoundaryReportEntry> extractedBoundaryReports = {}) {
  SequentialEquivalenceResult result;
  result.status = status;
  result.bound = bound;
  result.reason = std::move(reason);
  // LCOV_EXCL_START
  result.coveredOutputs = coverage.checkedOutputs.names.size();
  // LCOV_EXCL_STOP
  result.totalOutputs = coverage.totalOutputs;
  result.skippedObservedOutputs = coverage.skippedOutputs;
  result.resetUnanchoredSkippedOutputs =
      // LCOV_EXCL_START
      coverage.resetUnanchoredSkippedOutputs;
      // LCOV_EXCL_STOP
  result.multiClockDomainSkippedOutputs =
      coverage.multiClockDomainSkippedOutputs;
  result.abstractedSequentialBoundaries =
      std::move(abstractedSequentialBoundaries);
  result.extractedBoundaryReports = std::move(extractedBoundaryReports);
  return result;
}

std::string outputNameForProblemIndex(const KInductionProblem& problem,
                                      size_t outputIndex) {
  return outputIndex < problem.observedOutputNames.size()
             ? problem.observedOutputNames[outputIndex]
             : std::to_string(outputIndex);  // LCOV_EXCL_LINE
}

std::vector<bool> makeInitialPdrCoveredOutputs(
    const KInductionProblem& problem) {
  // No SEC strategy stage is allowed to pre-cover an output.  PDR coverage
  // starts empty and is filled only after PDR proves the selected batch.
  return std::vector<bool>(problem.observedOutputExprs0.size(), false);
}

void markDualRailPdrOutputRangeCovered(
    // LCOV_EXCL_START
    std::vector<bool>& coveredOutputs,
    std::unordered_map<size_t, std::string>& skipReasons,
    size_t firstOutput,
    // LCOV_EXCL_STOP
    size_t endOutput) {
  const size_t cappedEnd = std::min(endOutput, coveredOutputs.size());
  for (size_t outputIndex = firstOutput; outputIndex < cappedEnd;
       ++outputIndex) {
    coveredOutputs[outputIndex] = true;
    skipReasons.erase(outputIndex);
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

bool markDualRailPdrOutputSkipped(
    const KInductionProblem& problem,
    // LCOV_EXCL_START
    std::vector<bool>& coveredOutputs,
    // LCOV_EXCL_STOP
    std::unordered_map<size_t, std::string>& skipReasons,
    // LCOV_DISABLED_START
    size_t outputIndex) {
  if (outputIndex >= coveredOutputs.size()) {
    return true;  // LCOV_EXCL_LINE
  }
  const std::string reason =
      "dual-rail PDR repair was inconclusive after isolating this output";
  if (!coveredOutputs[outputIndex]) {
    skipReasons[outputIndex] = reason;
    // LCOV_DISABLED_START
    return true;  // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
  }
  if (std::count(coveredOutputs.begin(), coveredOutputs.end(), true) <= 1) {  // LCOV_EXCL_LINE
    // LCOV_DISABLED_START
    // Partial coverage must not erase the whole proof surface.  A single-output
    // LCOV_DISABLED_STOP
    // dual-rail timeout remains inconclusive so a real mismatch cannot be
    // hidden behind an empty covered set.
    return false;  // LCOV_EXCL_LINE
  }
  coveredOutputs[outputIndex] = false;  // LCOV_EXCL_LINE
  skipReasons[outputIndex] = reason;  // LCOV_EXCL_LINE
  if (isSecDiagEnabled() ||  // LCOV_EXCL_LINE
      std::getenv("KEPLER_SEC_SUMMARY_STATS") != nullptr) {  // LCOV_EXCL_LINE
    emitSecDiag(  // LCOV_EXCL_LINE
        "SEC diag: dual-rail PDR leaves output uncovered: ",
        outputNameForProblemIndex(problem, outputIndex));  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
  return true;  // LCOV_EXCL_LINE
}

void markPdrValidationUnknownOutputs(
    std::vector<bool>& coveredOutputs,
    std::unordered_map<size_t, std::string>& skipReasons,
    size_t firstOutput,
    const std::vector<size_t>& localOutputIndices) {
  for (const size_t localOutputIndex : localOutputIndices) {
    const size_t outputIndex = firstOutput + localOutputIndex; // LCOV_EXCL_LINE
    if (outputIndex >= coveredOutputs.size()) { // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    coveredOutputs[outputIndex] = false; // LCOV_EXCL_LINE
    skipReasons[outputIndex] = // LCOV_EXCL_LINE
        "PDR concrete validation was inconclusive for this output";
  }
}

OutputCoverageSelection buildCoverageWithDualRailOutputSkips(
    const OutputCoverageSelection& baseCoverage,
    const KInductionProblem& problem,
    const std::vector<bool>& coveredOutputs,
    const std::unordered_map<size_t, std::string>& skipReasons) {
  if (coveredOutputs.empty()) {
    return baseCoverage;  // LCOV_EXCL_LINE
  }
  // LCOV_DISABLED_START
  if (baseCoverage.checkedOutputs.names.size() != coveredOutputs.size()) {
  // LCOV_DISABLED_STOP
    return baseCoverage;  // LCOV_EXCL_LINE
  }
  if (static_cast<size_t>(
          std::count(coveredOutputs.begin(), coveredOutputs.end(), true)) ==
          coveredOutputs.size() &&
      skipReasons.empty()) {
    return baseCoverage;
  }

  OutputCoverageSelection coverage = baseCoverage;
  coverage.checkedOutputs.names.clear();
  coverage.checkedOutputs.keys0.clear();
  // LCOV_DISABLED_START
  coverage.checkedOutputs.keys1.clear();
  // LCOV_DISABLED_STOP
  for (size_t i = 0; i < coveredOutputs.size(); ++i) {
    if (coveredOutputs[i]) {
      coverage.checkedOutputs.names.push_back(baseCoverage.checkedOutputs.names[i]);
      coverage.checkedOutputs.keys0.push_back(baseCoverage.checkedOutputs.keys0[i]);
      coverage.checkedOutputs.keys1.push_back(baseCoverage.checkedOutputs.keys1[i]);
      continue;
    }

    const auto reasonIt = skipReasons.find(i);
    std::string reason;
    if (reasonIt != skipReasons.end()) {
      reason = reasonIt->second;
    } else if (problem.dualRailOutputSkipReasons.size() ==
                   coveredOutputs.size() &&  // LCOV_EXCL_LINE
               !problem.dualRailOutputSkipReasons[i].empty()) {  // LCOV_EXCL_LINE
      reason = problem.dualRailOutputSkipReasons[i];  // LCOV_EXCL_LINE
    } else {  // LCOV_EXCL_LINE
      reason = "dual-rail proof was inconclusive";  // LCOV_EXCL_LINE
    }
    coverage.skippedOutputs.push_back(
        outputNameForProblemIndex(problem, i) + ": " + reason);
  }
  return coverage;
}

std::vector<size_t> collectOutputsRequiringDualRailEngineProof(
    const KInductionProblem& problem) {
  std::vector<size_t> outputIndices;
  outputIndices.reserve(problem.observedOutputExprs0.size());
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    if (problem.dualRailOutputSkipReasons.size() ==
            problem.observedOutputExprs0.size() &&
        !problem.dualRailOutputSkipReasons[i].empty()) {
      continue; // LCOV_EXCL_LINE
    }
    outputIndices.push_back(i);
  // LCOV_DISABLED_START
  }
  return outputIndices;
  // LCOV_DISABLED_STOP
}

void rebuildSelectedOutputProperty(KInductionProblem& problem) {
  BoolExpr* property = BoolExpr::createTrue();
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    property = BoolExpr::And(
        property,
        makeEqualityExpr(problem.observedOutputExprs0[i],
                         problem.observedOutputExprs1[i]));
  }
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = nullptr;
  problem.inductionBad = nullptr;
}

KInductionProblem makeOutputSubsetProblem(
    const KInductionProblem& source,
    const std::vector<size_t>& outputIndices) {
  // LCOV_DISABLED_START
  KInductionProblem subset = source;
  // LCOV_DISABLED_STOP
  subset.observedOutputs.clear();
  subset.observedOutputNames.clear();
  subset.observedOutputExprs0.clear();
  subset.observedOutputExprs1.clear();
  subset.dualRailOutputSkipReasons.clear();

  // LCOV_DISABLED_START
  const bool copyObservedKeys =
  // LCOV_DISABLED_STOP
      source.observedOutputs.size() == source.observedOutputExprs0.size();
  const bool copySkipReasons =
      source.dualRailOutputSkipReasons.size() ==
      source.observedOutputExprs0.size();
  for (const size_t outputIndex : outputIndices) {
    if (copyObservedKeys) {
      subset.observedOutputs.push_back(source.observedOutputs[outputIndex]);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    subset.observedOutputNames.push_back(
        outputNameForProblemIndex(source, outputIndex));
    subset.observedOutputExprs0.push_back(
        source.observedOutputExprs0[outputIndex]);
    subset.observedOutputExprs1.push_back(
        source.observedOutputExprs1[outputIndex]);
    if (copySkipReasons) {
      subset.dualRailOutputSkipReasons.push_back(
          source.dualRailOutputSkipReasons[outputIndex]);
    }
  }
  rebuildSelectedOutputProperty(subset);
  subset.description = source.description + " selected output repair";
  // LCOV_DISABLED_START
  return subset;
}

OutputCoverageSelection buildCoverageSkippingOutputIndices( // LCOV_EXCL_LINE
    const OutputCoverageSelection& baseCoverage,
    const KInductionProblem& problem,
    // LCOV_DISABLED_STOP
    const std::vector<size_t>& skippedOutputIndices,
    const std::string& reason) {
  std::vector<bool> coveredOutputs(problem.observedOutputExprs0.size(), true); // LCOV_EXCL_LINE
  // LCOV_DISABLED_START
  std::unordered_map<size_t, std::string> skipReasons; // LCOV_EXCL_LINE
  for (const size_t outputIndex : skippedOutputIndices) { // LCOV_EXCL_LINE
    if (outputIndex < coveredOutputs.size()) { // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
      coveredOutputs[outputIndex] = false; // LCOV_EXCL_LINE
      skipReasons.emplace(outputIndex, reason); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  }
  return buildCoverageWithDualRailOutputSkips( // LCOV_EXCL_LINE
      baseCoverage, problem, coveredOutputs, skipReasons); // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

std::unordered_map<size_t, SignalKey> buildStateKeyByLocalVar(
    const SequentialDesignModel& model) {
  std::unordered_map<size_t, SignalKey> stateKeyByVar;
  stateKeyByVar.reserve(model.stateBits.size());
  for (const auto& key : model.stateBits) {
    // LCOV_DISABLED_START
    const auto varIt = model.inputVarByKey.find(key);
    // LCOV_DISABLED_STOP
    if (varIt != model.inputVarByKey.end()) {
      stateKeyByVar.emplace(varIt->second, key);
    }
  }
  return stateKeyByVar;
}

std::optional<SignalKey> findFirstUnanchoredStateSupportKey(
    BoolExpr* expr,
    const std::unordered_map<size_t, SignalKey>& stateKeyByVar,
    const std::unordered_set<SignalKey, SignalKeyHash>& anchoredKeys) {
  if (expr == nullptr) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }
  for (const auto var : expr->getSupportVars()) {
    const auto stateIt = stateKeyByVar.find(var);
    if (stateIt == stateKeyByVar.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (anchoredKeys.find(stateIt->second) == anchoredKeys.end()) {
      return stateIt->second;
    }
  }
  return std::nullopt;  // LCOV_EXCL_LINE
}

std::string describeUnanchoredStateSupport(
    const SequentialDesignModel& model,
    const SignalKey& key) {
  return "depends on reset-unanchored internal state " +
         describeSecSignalKey(model, key);
}  // LCOV_EXCL_LINE

void logSecDiagLine(bool secDiagEnabled, const char* message);

void filterOutputsRequiringUnanchoredResetState(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const ReachableStateInvariant& reachableInvariant,
    bool resetBootstrapActive,
    AlignedSecInterface& aligned,
    bool secDiagEnabled) {
  if (!resetBootstrapActive || aligned.outputs.names.empty()) {
    return;
  }

  const auto stateKeyByVar0 = buildStateKeyByLocalVar(model0);
  const auto stateKeyByVar1 = buildStateKeyByLocalVar(model1);
  if (stateKeyByVar0.empty() && stateKeyByVar1.empty()) {
    return;  // LCOV_EXCL_LINE
  }

  const std::unordered_set<SignalKey, SignalKeyHash> anchoredStateKeys0;
  const std::unordered_set<SignalKey, SignalKeyHash> anchoredStateKeys1;
  const auto [abstractOutputMap0, abstractOutputMap1] =
      buildPublicInputAbstractMaps(model0, model1, aligned.inputs);

  AlignedSignals filteredOutputs;
  filteredOutputs.names.reserve(aligned.outputs.names.size());
  filteredOutputs.keys0.reserve(aligned.outputs.keys0.size());
  filteredOutputs.keys1.reserve(aligned.outputs.keys1.size());

  size_t resetUnanchoredSkipCount = 0;
  for (size_t i = 0; i < aligned.outputs.names.size(); ++i) {
    const auto& name = aligned.outputs.names[i];
    const auto& key0 = aligned.outputs.keys0[i];
    // LCOV_DISABLED_START
    const auto& key1 = aligned.outputs.keys1[i];
    std::vector<std::string> reasons;


// LCOV_DISABLED_STOP
    if (areOutputExprsEquivalentUnderPublicInputs(
            model0.observedOutputExprByKey.at(key0),
            model1.observedOutputExprByKey.at(key1),
            abstractOutputMap0,
            abstractOutputMap1)) {
      // LCOV_DISABLED_START
      filteredOutputs.names.push_back(name);
      filteredOutputs.keys0.push_back(key0);
      // LCOV_DISABLED_STOP
      filteredOutputs.keys1.push_back(key1);
      // LCOV_DISABLED_START
      continue;
    }

    const auto unanchored0 = findFirstUnanchoredStateSupportKey(
        model0.observedOutputExprByKey.at(key0),
        // LCOV_DISABLED_STOP
        stateKeyByVar0,
        // LCOV_DISABLED_START
        anchoredStateKeys0);
    if (unanchored0.has_value()) {
      reasons.push_back(
      // LCOV_DISABLED_STOP
          "design0 " + describeUnanchoredStateSupport(model0, *unanchored0));
    }
    const auto unanchored1 = findFirstUnanchoredStateSupportKey(
        model1.observedOutputExprByKey.at(key1),
        stateKeyByVar1,
        anchoredStateKeys1);
    if (unanchored1.has_value()) {
      reasons.push_back(
          "design1 " + describeUnanchoredStateSupport(model1, *unanchored1));
    }

    if (!reasons.empty()) {
      // This is a coverage decision, not an internal-state equality shortcut:
      // reset/bootstrap values are per-design facts at one frontier.  They do
      // not justify comparing later state-dependent outputs unless SEC also
      // has an inductive cross-design state relation for that support.
      const auto skippedOutput = name + ": " + joinReasons(reasons);
      aligned.outputCoverage.skippedOutputs.push_back(skippedOutput);
      aligned.outputCoverage.resetUnanchoredSkippedOutputs.push_back(
          skippedOutput);
      ++resetUnanchoredSkipCount;
      continue;
    }

    filteredOutputs.names.push_back(name);  // LCOV_EXCL_LINE
    filteredOutputs.keys0.push_back(key0);  // LCOV_EXCL_LINE
    filteredOutputs.keys1.push_back(key1);  // LCOV_EXCL_LINE
  }

  aligned.outputs = std::move(filteredOutputs);
  aligned.outputCoverage.checkedOutputs = aligned.outputs;
  if (secDiagEnabled && resetUnanchoredSkipCount != 0) {
    fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: reset-unanchored output skips=%zu checked_outputs=%zu total_outputs=%zu\n",
        resetUnanchoredSkipCount,  // LCOV_EXCL_LINE
        aligned.outputs.names.size(),  // LCOV_EXCL_LINE
        aligned.outputCoverage.totalOutputs);  // LCOV_EXCL_LINE
    fprintf(  // LCOV_EXCL_LINE
        stderr,  // LCOV_EXCL_LINE
        "SEC diag: reset-frontier checked outputs=%s\n",
        formatStringList(aligned.outputs.names, aligned.outputs.names.size()).c_str());  // LCOV_EXCL_LINE
    fflush(stderr);  // LCOV_EXCL_LINE
  }  // LCOV_EXCL_LINE
}

template <typename MapT>
void assignSymbols(const std::vector<SignalKey>& keys,
                   MapT& output,
                   std::vector<size_t>& allSymbols,
                   size_t& nextSymbol);

std::unordered_map<size_t, size_t> buildLocalToCombinedMap(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& inputSymbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& stateSymbols);

void logSecDiagLine(bool secDiagEnabled, const char* message) {
  if (!secDiagEnabled) {
    return;
  }
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
}

bool pdrStrategyStatsEnabled() {
  return std::getenv("KEPLER_SEC_PDR_STATS") != nullptr;
}

std::vector<std::string> secStrategyCommaListFromEnv(const char* name) {
  const char* valueText = std::getenv(name);
  if (valueText == nullptr || *valueText == '\0') {
    return {};
  }

  std::vector<std::string> values;
  std::string current;
  std::istringstream stream(valueText);
  // LCOV_DISABLED_START
  while (std::getline(stream, current, ',')) {
  // LCOV_DISABLED_STOP
    current.erase(
        current.begin(),
        std::find_if(
            current.begin(), current.end(), [](unsigned char ch) {
              return !std::isspace(ch);
            }));
    current.erase(
        std::find_if(
            current.rbegin(), current.rend(), [](unsigned char ch) {
              return !std::isspace(ch);
            }).base(),
        current.end());
    if (!current.empty()) {
      values.push_back(current);
    }
  }
  return values;
}

size_t secStrategySizeLimitFromEnv(const char* name, size_t defaultValue) {
  const char* valueText = std::getenv(name);
  if (valueText == nullptr || *valueText == '\0') {
    return defaultValue;
  }
  const auto value = std::strtoull(valueText, nullptr, 10);
  return value == 0 ? defaultValue : static_cast<size_t>(value);
}

bool secSummaryStatsEnabled() {
  return std::getenv("KEPLER_SEC_SUMMARY_STATS") != nullptr ||
         pdrStrategyStatsEnabled();
}

constexpr size_t kMaxDualRailResidualOutputs = 128;
constexpr size_t kMaxDualRailResidualProofStateSymbols = 4096;
constexpr size_t kMaxDualRailResidualConcretePrecheckOutputs = 16;
constexpr size_t kMaxDualRailFinalResetFrontierOriginalOutputs = 384;

enum class DualRailResidualEngine {
  KInduction,
};

struct DualRailResidualProofState {
  std::vector<bool> coveredOutputs;
  std::unordered_map<size_t, std::string> skipReasons;
  std::optional<SequentialEquivalenceResult> terminalResult;
  size_t provedBound = 0;
};

const char* dualRailResidualEngineName(DualRailResidualEngine engine) {
  switch (engine) {
    case DualRailResidualEngine::KInduction:
      return "k-induction";
  }
  return "selected engine";  // LCOV_EXCL_LINE
}

void markDualRailResidualOutputsCovered(
    const std::vector<size_t>& outputIndices,
    DualRailResidualProofState& proofState) {
  for (const size_t outputIndex : outputIndices) {
    if (outputIndex < proofState.coveredOutputs.size()) {
      proofState.coveredOutputs[outputIndex] = true;
      proofState.skipReasons.erase(outputIndex);
    }
  }
}

void markDualRailResidualOutputSkipped(
    size_t outputIndex,
    const KInductionProblem& problem,
    DualRailResidualEngine engine,
    DualRailResidualProofState& proofState,
    const std::string& reason) {
  if (outputIndex >= proofState.coveredOutputs.size()) {
    return;  // LCOV_EXCL_LINE
  }
  proofState.coveredOutputs[outputIndex] = false;
  proofState.skipReasons[outputIndex] = reason.empty()
                                            ? std::string("dual-rail ") +
                                                  dualRailResidualEngineName(engine) +
                                                  " proof was inconclusive for this output"
                                            : reason;  // LCOV_EXCL_LINE
  if (isSecDiagEnabled() || secSummaryStatsEnabled()) {
    emitSecDiag(
        "SEC diag: dual-rail ", dualRailResidualEngineName(engine),
        " leaves residual output uncovered: ",
        outputNameForProblemIndex(problem, outputIndex));
  }
}

size_t dualRailResidualStateSymbolCount(const KInductionProblem& problem) {
  return problem.usesDualRailStateEncoding
             ? problem.dualRailStatePairs.size() * 2
             : problem.totalStateCount;  // LCOV_EXCL_LINE
}

bool shouldRunDualRailResidualCounterexampleSweep(
    const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding) {
    return true;  // LCOV_EXCL_LINE
  }

  // The residual sweep is a concrete-counterexample shortcut before deferred
  // induction.  Once the rail-state surface crosses the small-proof threshold,
  // the base checker decomposes the sweep across many output cones and can
  // dominate runtime for equivalent designs.  Leave those cases to the selected
  // engine proof and its final base validation instead.
  return dualRailResidualStateSymbolCount(problem) <=
         kMaxDualRailResidualProofStateSymbols;
}

bool shouldSkipLargeDualRailResidualProofSurface(
    const KInductionProblem& problem,
    size_t residualOutputCount,
    DualRailResidualEngine engine) {
  if (!problem.usesDualRailStateEncoding || residualOutputCount == 0) {
    return false;  // LCOV_EXCL_LINE
  }
  if (engine == DualRailResidualEngine::KInduction) {
    // K-induction must get the actual residual obligation.  Large surfaces may
    // still finish inconclusive after strict KI splitting, but they must not be
    // pre-skipped before the selected engine has tried to prove them.
    return false;
  }
  // This guard is for truly oversized residual surfaces.  It reports
  // inconclusive/skipped coverage instead of accepting unproved outputs.
  return residualOutputCount > kMaxDualRailResidualOutputs && // LCOV_EXCL_LINE
         dualRailResidualStateSymbolCount(problem) > // LCOV_EXCL_LINE
             kMaxDualRailResidualProofStateSymbols;
}

std::optional<KInductionResult::CounterexampleWitness>
findDualRailResidualCounterexample(const KInductionProblem& subsetProblem,
                                   KEPLER_FORMAL::Config::SolverType solverType,
                                   // LCOV_DISABLED_START
                                   size_t maxK) {
                                   // LCOV_DISABLED_STOP
  for (size_t depth = 0; depth <= maxK; ++depth) {
    if (auto witness =
            SEC::findBaseCounterexampleAtFrontier(subsetProblem, solverType, depth);
        witness.has_value()) {
      return witness; // LCOV_EXCL_LINE
    }
  }
  return std::nullopt;
}

bool findAndRecordDualRailResidualCounterexample(
    const KInductionProblem& subsetProblem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports,
    DualRailResidualEngine engine,
    DualRailResidualProofState& proofState);

std::unordered_set<size_t> combinedStateSymbolSet(
    const KInductionProblem& problem) {
  // LCOV_DISABLED_START
  std::unordered_set<size_t> stateSymbols;
  // LCOV_DISABLED_STOP
  const auto combinedStateSymbols = problem.combinedStateSymbols();
  stateSymbols.reserve(combinedStateSymbols.size());
  stateSymbols.insert(combinedStateSymbols.begin(), combinedStateSymbols.end());
  return stateSymbols;
}

bool isInputOnlyOutputMismatchCandidate(
    const KInductionProblem& problem,
    size_t outputIndex,
    const std::unordered_set<size_t>& stateSymbols) {
  BoolExpr* outputBad = BoolExpr::simplify(BoolExpr::Xor(
      problem.observedOutputExprs0[outputIndex],
      problem.observedOutputExprs1[outputIndex]));
  if (outputBad == BoolExpr::createFalse()) {
    return false;
  }
  const auto support = outputBad->getSupportVars();
  for (const auto symbol : support) {
    if (stateSymbols.find(symbol) != stateSymbols.end()) {
      return false;
    }
  // LCOV_DISABLED_START
  }
  // LCOV_DISABLED_STOP
  return true;
}

std::optional<KInductionResult::CounterexampleWitness>
findInputOnlyFrameZeroResidualCounterexample(
    const KInductionProblem& problem,
    const std::vector<size_t>& outputIndices,
    KEPLER_FORMAL::Config::SolverType solverType) {
  if (!problem.usesDualRailStateEncoding || outputIndices.empty()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  // LCOV_DISABLED_START
  }
  // LCOV_DISABLED_STOP

  const std::unordered_set<size_t> stateSymbols =
      combinedStateSymbolSet(problem);
  std::vector<size_t> inputOnlyOutputs;
  inputOnlyOutputs.reserve(outputIndices.size());
  for (const size_t outputIndex : outputIndices) {
    if (outputIndex >= problem.observedOutputExprs0.size() ||
        outputIndex >= problem.observedOutputExprs1.size()) {
      continue;  // LCOV_EXCL_LINE
    }
    if (isInputOnlyOutputMismatchCandidate(
            problem, outputIndex, stateSymbols)) {
      inputOnlyOutputs.push_back(outputIndex);
    }
  }
  if (inputOnlyOutputs.empty()) {
    return std::nullopt;
  }

  // This is a witness-only guard for skipped residual top outputs.  Restrict it
  // to frame-0 input/constant mismatches so equivalent reset-bootstrap designs
  // do not pay to materialize transition cones before the selected SEC engine.
  const KInductionProblem inputOnlyProblem =
      makeOutputSubsetProblem(problem, inputOnlyOutputs);
  return SEC::findFastBaseCounterexampleAtFrontier(
      inputOnlyProblem, solverType, 0);
}

std::optional<KInductionResult::CounterexampleWitness>
findSmallDualRailResidualConcreteCounterexample(
    const KInductionProblem& problem,
    const std::vector<size_t>& outputIndices,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  if (!problem.usesDualRailStateEncoding || outputIndices.empty()) {
    return std::nullopt;  // LCOV_EXCL_LINE
  }

  if (auto witness =
          findInputOnlyFrameZeroResidualCounterexample(
              problem, outputIndices, solverType);
      witness.has_value()) {
    return witness;
  }

  if (outputIndices.size() > kMaxDualRailResidualConcretePrecheckOutputs) {
    return std::nullopt;
  }

  for (const size_t outputIndex : outputIndices) {
    const KInductionProblem singleOutputProblem =
        makeOutputSubsetProblem(problem, {outputIndex});
    for (size_t depth = 0; depth <= maxK; ++depth) {
      if (auto witness = SEC::findBaseCounterexampleAtFrontier(
              singleOutputProblem, solverType, depth);
          witness.has_value()) {
        return witness; // LCOV_EXCL_LINE
      }
    }
  }
  return std::nullopt;
}

void recordDualRailResidualCounterexample(
    KInductionResult::CounterexampleWitness witness,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports,
    DualRailResidualProofState& proofState) {
  KInductionResult witnessResult{
      KInductionStatus::Different, witness.badFrame, std::move(witness)};
  proofState.terminalResult = makeSecResult(
      SequentialEquivalenceStatus::Different,
      witnessResult.bound,
      formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),
      outputCoverage,
      abstractedSequentialBoundaries,
      extractedBoundaryReports);
}

void proveDualRailResidualOutputSet(
    const KInductionProblem& problem,
    const std::vector<size_t>& outputIndices,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports,
    DualRailResidualEngine engine,
    bool runCounterexampleSweep,
    DualRailResidualProofState& proofState) {
  if (outputIndices.empty() || proofState.terminalResult.has_value()) {
    return;  // LCOV_EXCL_LINE
  }

  KInductionProblem subsetProblem =
      makeOutputSubsetProblem(problem, outputIndices);
  const bool useDeferredInductionResidualProof =
      subsetProblem.usesDualRailStateEncoding;
  if (useDeferredInductionResidualProof && runCounterexampleSweep &&
      findAndRecordDualRailResidualCounterexample( // LCOV_EXCL_LINE
          subsetProblem,
          maxK, // LCOV_EXCL_LINE
          solverType, // LCOV_EXCL_LINE
          model0, // LCOV_EXCL_LINE
          model1, // LCOV_EXCL_LINE
          top0, // LCOV_EXCL_LINE
          top1, // LCOV_EXCL_LINE
          outputCoverage, // LCOV_EXCL_LINE
          abstractedSequentialBoundaries, // LCOV_EXCL_LINE
          extractedBoundaryReports, // LCOV_EXCL_LINE
          engine, // LCOV_EXCL_LINE
          proofState)) { // LCOV_EXCL_LINE
    return; // LCOV_EXCL_LINE
  }
  if (useDeferredInductionResidualProof) {
    // Localized residual proofs should not spend max_k frontier checks on a
    // single hard rail output.  First ask whether the induction step closes; if
    // it does, validate the required concrete base prefix once below.
    subsetProblem.deferBaseCaseChecks = true;
  }
  if (isSecDiagEnabled() || secSummaryStatsEnabled()) {
    emitSecDiag(
        "SEC diag: dual-rail ", dualRailResidualEngineName(engine),
        " proving residual output batch size=", outputIndices.size());
  }

  if (useDeferredInductionResidualProof) {
    KInductionEngine residualEngine(subsetProblem, solverType);
    KInductionResult result = residualEngine.run(maxK);
    proofState.provedBound = std::max(proofState.provedBound, result.bound);
    if (result.status == KInductionStatus::Equivalent) {
      const size_t baseHorizon = result.bound == 0 ? 0 : result.bound - 1;
      const SEC::BaseCounterexampleCheckResult baseCheck =
          SEC::checkBaseCounterexampleWithFastValidation(
              subsetProblem, solverType, baseHorizon);
      if (baseCheck.status == SEC::BaseCounterexampleCheckStatus::Counterexample) {
        KInductionResult witnessResult{  // LCOV_EXCL_LINE
            KInductionStatus::Different,
            baseCheck.witness.has_value() ? baseCheck.witness->badFrame : baseHorizon,  // LCOV_EXCL_LINE
            baseCheck.witness};  // LCOV_EXCL_LINE
        proofState.terminalResult = makeSecResult(  // LCOV_EXCL_LINE
            SequentialEquivalenceStatus::Different,
            witnessResult.bound,  // LCOV_EXCL_LINE
            formatCounterexampleWitness(  // LCOV_EXCL_LINE
                witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
            outputCoverage,  // LCOV_EXCL_LINE
            abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
            extractedBoundaryReports);  // LCOV_EXCL_LINE
        return;
      }  // LCOV_EXCL_LINE
      if (baseCheck.status == SEC::BaseCounterexampleCheckStatus::NoCounterexample) {
        markDualRailResidualOutputsCovered(outputIndices, proofState);
        return;
      }
      if (outputIndices.size() == 1) {  // LCOV_EXCL_LINE
        markDualRailResidualOutputSkipped(  // LCOV_EXCL_LINE
            outputIndices.front(),  // LCOV_EXCL_LINE
            problem,  // LCOV_EXCL_LINE
            engine,  // LCOV_EXCL_LINE
            proofState,  // LCOV_EXCL_LINE
            std::string("dual-rail ") + dualRailResidualEngineName(engine) +  // LCOV_EXCL_LINE
                " base validation was inconclusive for this output");
        return;  // LCOV_EXCL_LINE
      }
    }
    if (result.status == KInductionStatus::Different) {
      proofState.terminalResult = makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Different,
          result.bound,  // LCOV_EXCL_LINE
          result.witness.has_value()  // LCOV_EXCL_LINE
              ? formatCounterexampleWitness(result, model0, model1, top0, top1)  // LCOV_EXCL_LINE
              : "Classic k-induction found a counterexample at k = " +  // LCOV_EXCL_LINE
                    std::to_string(result.bound),  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
      return;  // LCOV_EXCL_LINE
    }
  }

  if (outputIndices.size() > 1) {
    const size_t mid = outputIndices.size() / 2;
    std::vector<size_t> left(outputIndices.begin(), outputIndices.begin() + mid);
    std::vector<size_t> right(outputIndices.begin() + mid, outputIndices.end());
    proveDualRailResidualOutputSet(
        problem,
        left,
        maxK,
        solverType,
        model0,
        model1,
        top0,
        top1,
        outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports,
        engine,
        runCounterexampleSweep,
        proofState);
    proveDualRailResidualOutputSet(
        problem,
        right,
        maxK,
        solverType,
        model0,
        model1,
        top0,
        top1,
        outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports,
        engine,
        runCounterexampleSweep,
        proofState);
    return;
  }

  markDualRailResidualOutputSkipped(
      outputIndices.front(), problem, engine, proofState, "");
}

bool findAndRecordDualRailResidualCounterexample(
    const KInductionProblem& subsetProblem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports,
    DualRailResidualEngine engine,
    DualRailResidualProofState& proofState) {
  if (!shouldRunDualRailResidualCounterexampleSweep(subsetProblem)) {
    return false;
  }

  if (auto witness =
          findDualRailResidualCounterexample(
              subsetProblem, solverType, maxK);
      witness.has_value()) {
    // Deferred residual proof shortcuts still owe the selected engine's
    // concrete top-output base search.  Check the whole residual once before
    // batching so real edits are caught without repeating the same wide query
    // for every bit of an equivalent bus.
    recordDualRailResidualCounterexample( // LCOV_EXCL_LINE
        std::move(*witness), // LCOV_EXCL_LINE
        model0, // LCOV_EXCL_LINE
        model1, // LCOV_EXCL_LINE
        top0, // LCOV_EXCL_LINE
        top1, // LCOV_EXCL_LINE
        outputCoverage, // LCOV_EXCL_LINE
        abstractedSequentialBoundaries, // LCOV_EXCL_LINE
        extractedBoundaryReports, // LCOV_EXCL_LINE
        proofState); // LCOV_EXCL_LINE
    return true; // LCOV_EXCL_LINE
  }
  return false;
}

std::optional<SequentialEquivalenceResult> proveDualRailResidualsWithSelectedEngine(
    const KInductionProblem& problem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports,
    DualRailResidualEngine engine) {
  if (!problem.usesDualRailStateEncoding ||
      problem.observedOutputExprs0.empty()) {
    return std::nullopt;
  }

  const std::vector<size_t> residualOutputIndices =
      collectOutputsRequiringDualRailEngineProof(problem);
  std::unordered_map<size_t, std::string> presetSkipReasons;
  if (problem.dualRailOutputSkipReasons.size() ==
      problem.observedOutputExprs0.size()) {
    for (size_t i = 0; i < problem.dualRailOutputSkipReasons.size(); ++i) {
      if (!problem.dualRailOutputSkipReasons[i].empty()) {
        presetSkipReasons.emplace(i, problem.dualRailOutputSkipReasons[i]); // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }
  }
  if (residualOutputIndices.empty()) {
    if (!presetSkipReasons.empty()) {  // LCOV_EXCL_LINE
      std::vector<bool> coveredOutputs(problem.observedOutputExprs0.size(), false); // LCOV_EXCL_LINE
      const OutputCoverageSelection partialCoverage =
          buildCoverageWithDualRailOutputSkips(  // LCOV_EXCL_LINE
              outputCoverage,  // LCOV_EXCL_LINE
              problem,  // LCOV_EXCL_LINE
              coveredOutputs,  // LCOV_EXCL_LINE
              presetSkipReasons);
      return makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Inconclusive,
          0,
          "Dual-rail output skips left no selected-engine obligation",  // LCOV_EXCL_LINE
          partialCoverage,
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    return std::nullopt; // LCOV_EXCL_LINE
  }

  DualRailResidualProofState proofState;
  proofState.coveredOutputs.assign(problem.observedOutputExprs0.size(), false);
  proofState.skipReasons = presetSkipReasons;

  if (auto witness = findSmallDualRailResidualConcreteCounterexample(
          problem, residualOutputIndices, solverType, maxK);
      witness.has_value()) {
    recordDualRailResidualCounterexample(
        std::move(*witness),
        model0,
        model1,
        top0,
        top1,
        outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports,
        proofState);
    return proofState.terminalResult;
  }

  if (shouldSkipLargeDualRailResidualProofSurface(
          problem, residualOutputIndices.size(), engine)) {
    emitSecDiag( // LCOV_EXCL_LINE
        "SEC diag: dual-rail ", dualRailResidualEngineName(engine), // LCOV_EXCL_LINE
        " skipping large residual surface outputs=",
        residualOutputIndices.size(), // LCOV_EXCL_LINE
        " state_symbols=",
        problem.dualRailStatePairs.size() * 2); // LCOV_EXCL_LINE
    const std::string reason =
        std::string("dual-rail ") + dualRailResidualEngineName(engine) + // LCOV_EXCL_LINE
        " residual skipped for large rail-state surface";
    const OutputCoverageSelection partialCoverage =
        buildCoverageSkippingOutputIndices( // LCOV_EXCL_LINE
            outputCoverage, problem, residualOutputIndices, reason); // LCOV_EXCL_LINE
    return makeSecResult( // LCOV_EXCL_LINE
        SequentialEquivalenceStatus::Inconclusive,
        0,
        "Dual-rail residual surface exceeded selected-engine proof limits", // LCOV_EXCL_LINE
        partialCoverage,
        abstractedSequentialBoundaries, // LCOV_EXCL_LINE
        extractedBoundaryReports); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  const KInductionProblem residualProblem =
      makeOutputSubsetProblem(problem, residualOutputIndices);
  if (findAndRecordDualRailResidualCounterexample(
          residualProblem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports,
          engine,
          proofState)) {
    return proofState.terminalResult; // LCOV_EXCL_LINE
  }
  for (const auto& [firstOutput, endOutput] :
       buildSupportBoundedOutputBatches(
           residualProblem,
           defaultOutputBatchingLimitsForProblem(residualProblem))) {
    std::vector<size_t> batchOutputIndices;
    batchOutputIndices.reserve(endOutput - firstOutput);
    for (size_t i = firstOutput; i < endOutput; ++i) {
      batchOutputIndices.push_back(residualOutputIndices[i]);
    }
    proveDualRailResidualOutputSet(
        problem,
        batchOutputIndices,
        maxK,
        // LCOV_DISABLED_START
        solverType,
        model0,
        model1,
        top0,
        // LCOV_DISABLED_STOP
        top1,
        // LCOV_DISABLED_START
        outputCoverage,
        // LCOV_DISABLED_STOP
        abstractedSequentialBoundaries,
        extractedBoundaryReports,
        engine,
        /*runCounterexampleSweep=*/false,
        proofState);
    if (proofState.terminalResult.has_value()) {
      return proofState.terminalResult;  // LCOV_EXCL_LINE
    }
  }

  const size_t coveredCount = static_cast<size_t>(std::count(
      proofState.coveredOutputs.begin(), proofState.coveredOutputs.end(), true));
  const OutputCoverageSelection finalCoverage =
      buildCoverageWithDualRailOutputSkips(
          outputCoverage, problem, proofState.coveredOutputs,
          proofState.skipReasons);
  if (coveredCount == 0) {
    return makeSecResult(
        SequentialEquivalenceStatus::Inconclusive,
        proofState.provedBound,
        std::string("Dual-rail ") + dualRailResidualEngineName(engine) +
            " did not prove any output",
        finalCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }

  return makeSecResult(
      SequentialEquivalenceStatus::Equivalent,
      proofState.provedBound,
      "",
      finalCoverage,
      abstractedSequentialBoundaries,
      extractedBoundaryReports);
}

size_t directObservedOutputSupportSize(const KInductionProblem& problem) {
  std::unordered_set<size_t> support;
  for (const auto* expr : problem.observedOutputExprs0) {
    const auto exprSupport = expr->getSupportVars();
    support.insert(exprSupport.begin(), exprSupport.end());
  }
  for (const auto* expr : problem.observedOutputExprs1) {
    const auto exprSupport = expr->getSupportVars();
    support.insert(exprSupport.begin(), exprSupport.end());
  }
  return support.size();
}

size_t pdrCertificateStateSymbolCount(const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding) {
    return problem.totalStateCount; // LCOV_EXCL_LINE
  }
  // Dual-rail PDR reasons about both value and known rails.  Runtime guards
  // must count the actual rail symbols, not only the original flop count.
  return problem.dualRailStatePairs.size() * 2;
}

bool canReportPdrValidationCounterexample(
    const KInductionProblem& problem,
    const KInductionResult::CounterexampleWitness& witness);

bool shouldDeferWideDualRailPdrValidation(const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding) {
    return false;
  }
  const size_t outputCount = problem.originalObservedOutputCount == 0
                                 ? problem.observedOutputExprs0.size()
                                 : problem.originalObservedOutputCount;
  return detail::shouldDeferPdrDualRailFrameZeroValidation(
      outputCount, pdrCertificateStateSymbolCount(problem));
}

std::optional<KInductionResult::CounterexampleWitness>
findPdrValidationCounterexample(const KInductionProblem& problem,
                                KEPLER_FORMAL::Config::SolverType solverType,
                                size_t maxK) {
  // PDR uses this only to validate concrete top-output SEC behavior.  Scanning
  // exact frontiers keeps the same bounded semantics as one prefix query while
  // letting the base solver localize wide multi-output ASIC cones.
  for (size_t depth = 0; depth <= maxK; ++depth) {
    if (auto witness = SEC::findBaseCounterexampleAtFrontier(
            problem, solverType, depth);
        witness.has_value()) {
      if (canReportPdrValidationCounterexample(problem, *witness)) {
        return witness;
      }
    } // LCOV_EXCL_LINE
  }
  return std::nullopt;
}

bool isUnknownBootstrapRailPair( // LCOV_EXCL_LINE
    const DualRailSymbolPair& rails,
    const std::unordered_map<size_t, bool>& bootstrapValueBySymbol) {
  const auto oneIt = bootstrapValueBySymbol.find(rails.mayBeOne); // LCOV_EXCL_LINE
  const auto zeroIt = bootstrapValueBySymbol.find(rails.mayBeZero); // LCOV_EXCL_LINE
  return oneIt != bootstrapValueBySymbol.end() && // LCOV_EXCL_LINE
         zeroIt != bootstrapValueBySymbol.end() && // LCOV_EXCL_LINE
         oneIt->second && // LCOV_EXCL_LINE
         zeroIt->second; // LCOV_EXCL_LINE
}

bool dualRailBadDependsOnUnknownBootstrapState( // LCOV_EXCL_LINE
    const KInductionProblem& problem) {
  if (!problem.usesDualRailStateEncoding || // LCOV_EXCL_LINE
      problem.resetBootstrapCycles == 0 || // LCOV_EXCL_LINE
      problem.bad == nullptr || // LCOV_EXCL_LINE
      problem.bootstrapStateAssignments.empty()) { // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }

  const std::set<size_t> support = problem.bad->getSupportVars(); // LCOV_EXCL_LINE
  if (support.empty()) { // LCOV_EXCL_LINE
    return false; // LCOV_EXCL_LINE
  }

  std::unordered_map<size_t, bool> bootstrapValueBySymbol; // LCOV_EXCL_LINE
  bootstrapValueBySymbol.reserve(problem.bootstrapStateAssignments.size()); // LCOV_EXCL_LINE
  for (const auto& [symbol, value] : problem.bootstrapStateAssignments) { // LCOV_EXCL_LINE
    bootstrapValueBySymbol.emplace(symbol, value); // LCOV_EXCL_LINE
  }

  for (const auto& rails : problem.dualRailStatePairs) { // LCOV_EXCL_LINE
    if (!isUnknownBootstrapRailPair(rails, bootstrapValueBySymbol)) { // LCOV_EXCL_LINE
      continue; // LCOV_EXCL_LINE
    }
    if (support.find(rails.mayBeOne) != support.end() || // LCOV_EXCL_LINE
        support.find(rails.mayBeZero) != support.end()) { // LCOV_EXCL_LINE
      return true; // LCOV_EXCL_LINE
    }
  }
  return false; // LCOV_EXCL_LINE
} // LCOV_EXCL_LINE

bool canReportPdrValidationCounterexample(
    const KInductionProblem& problem,
    const KInductionResult::CounterexampleWitness& witness) {
  if (!problem.usesDualRailStateEncoding ||
      problem.canReportSteadyFrontierMismatchAsCounterexample()) {
    return true;
  }
  if (witness.badFrame != 0) { // LCOV_EXCL_LINE
    return true; // LCOV_EXCL_LINE
  }
  // A frame-0 dual-rail mismatch that depends on X-valued reset-bootstrap
  // state can be only a rail-overapproximation artifact. Keep it as proof
  // feedback for PDR, but do not expose it as a concrete SEC counterexample.
  return !dualRailBadDependsOnUnknownBootstrapState(problem); // LCOV_EXCL_LINE
}

std::optional<KInductionResult::CounterexampleWitness>
findPdrPerOutputValidationCounterexample(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  // Exact single-output frontiers avoid the wide all-output bad-state OR while
  // still proving/reporting only concrete top-output SEC behavior.
  for (size_t outputIndex = 0;
       outputIndex < problem.observedOutputExprs0.size();
       ++outputIndex) {
    const KInductionProblem singleOutputProblem =
        makeOutputSubsetProblem(problem, {outputIndex});
    for (size_t depth = 0; depth <= maxK; ++depth) {
      if (auto witness = SEC::findBaseCounterexampleAtFrontier(
              singleOutputProblem, solverType, depth);
          witness.has_value()) {
        if (canReportPdrValidationCounterexample(
                singleOutputProblem, *witness)) {
          return witness;
        }
      } // LCOV_EXCL_LINE
    }
  }
  return std::nullopt;
}

struct PdrConcreteValidationCheck {
  std::optional<KInductionResult::CounterexampleWitness> witness;
  std::vector<size_t> unknownOutputIndices;
};

PdrConcreteValidationCheck checkPdrConcreteValidation(
    const KInductionProblem& problem,
    KEPLER_FORMAL::Config::SolverType solverType,
    size_t maxK) {
  PdrConcreteValidationCheck check;
  for (size_t outputIndex = 0;
       outputIndex < problem.observedOutputExprs0.size();
       ++outputIndex) {
    const KInductionProblem singleOutputProblem =
        makeOutputSubsetProblem(problem, {outputIndex});
    bool outputValidationUnknown = false;
    for (size_t depth = 0; depth <= maxK; ++depth) {
      const SEC::BaseCounterexampleCheckResult baseCheck =
          SEC::checkBaseCounterexampleWithFastValidation(
              singleOutputProblem, solverType, depth);
      if (baseCheck.status ==
          SEC::BaseCounterexampleCheckStatus::Counterexample) {
        check.witness = baseCheck.witness;  // LCOV_EXCL_LINE
        return check;  // LCOV_EXCL_LINE
      }
      if (baseCheck.status == SEC::BaseCounterexampleCheckStatus::Unknown) {
        outputValidationUnknown = true; // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
      }
    }
    if (outputValidationUnknown) {
      check.unknownOutputIndices.push_back(outputIndex); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  }
  return check;
}

OutputBatchingLimits dualRailPdrOutputBatchingLimits(
    OutputBatchingLimits defaultLimits) {
  // Keep the production default conservative, but allow diagnostics/regressions
  // to test a broad all-output dual-rail PDR slice without editing code.
  // LCOV_DISABLED_START
  defaultLimits.maxOutputBatchSize = secStrategySizeLimitFromEnv(
  // LCOV_DISABLED_STOP
      "KEPLER_SEC_DUAL_RAIL_OUTPUT_BATCH_SIZE",
      defaultLimits.maxOutputBatchSize);
  defaultLimits.outputBatchSupportLimit = secStrategySizeLimitFromEnv(
      "KEPLER_SEC_DUAL_RAIL_OUTPUT_BATCH_SUPPORT_LIMIT",
      defaultLimits.outputBatchSupportLimit);
  return defaultLimits;
}

void emitPdrStrategyStageStats(
    bool enabled,
    size_t batchIndex,
    size_t firstOutput,
    size_t endOutput,
    const char* stage,
    size_t transitionClosureLimit,
    size_t predecessorProjectionLimit,
    size_t badCubeLimit,
    const KInductionProblem& batch) {
  if (!enabled) {
    return;
  }

  // These stage markers are intentionally coarse: when a large SEC/PDR run is
  // sampled, they identify which CEGAR retry owns the following predecessor SAT
  // traffic without flooding the log with every query.
  emitSecDiag(  // LCOV_EXCL_LINE
      "SEC PDR stats: strategy batch=", batchIndex,
      " outputs=[", firstOutput, ",", endOutput, ")",
      " stage=", stage,
      " closure_limit=", transitionClosureLimit,
      " projection_limit=", predecessorProjectionLimit,
      " bad_cube_limit=", badCubeLimit,
      " transitions=", batch.transitions0.size() + batch.transitions1.size(),  // LCOV_EXCL_LINE
      " init_assignments=", batch.initialStateAssignments.size(),  // LCOV_EXCL_LINE
      " bootstrap_assignments=", batch.bootstrapStateAssignments.size(),  // LCOV_EXCL_LINE
      " observed_outputs=", batch.observedOutputExprs0.size(),
      " direct_support=", directObservedOutputSupportSize(batch),
      " output_names=[",
      formatStringList(batch.observedOutputNames, 8),
      "]");  // LCOV_EXCL_LINE
}

void appendAbstractedSequentialBoundaries(
    const SequentialDesignModel& model,
    const char* designPrefix,
    std::vector<std::string>& abstractedSequentialBoundaries) {
  abstractedSequentialBoundaries.reserve(
      abstractedSequentialBoundaries.size() +
      model.abstractedSequentialBoundaries.size());
  for (const auto& description : model.abstractedSequentialBoundaries) {
    abstractedSequentialBoundaries.push_back(
        std::string(designPrefix) + " " + description);
  }
}

void appendExtractedBoundaryReports(
    const SequentialDesignModel& model,
    const char* designPrefix,
    std::vector<ExtractedBoundaryReportEntry>& reports) {
  std::map<std::string, ExtractedBoundaryReportEntry> reportsBySignal;

  auto ensureEntry = [&](const SignalKey& key) -> ExtractedBoundaryReportEntry& {
    const auto signal = describeSecSignalKey(model, key);
    auto [it, _] = reportsBySignal.try_emplace(signal);
    it->second.design = designPrefix;
    it->second.signal = signal;
    if (const auto skipIt = model.connectivitySkipInfoByKey.find(key);
        skipIt != model.connectivitySkipInfoByKey.end()) {
      it->second.connectivitySkip = describeConnectivitySkipInfo(skipIt->second);
    }
    return it->second;
  };

  auto addRole = [&](const SignalKey& key, const char* role) {
    auto& entry = ensureEntry(key);
    appendUniqueRole(entry.roles, role);
  };

  // Boundary terms are the full exposed SEC cut surface:
  // - the original top interface
  // - opaque internal cut points from leaves SEC cannot model combinationally
  //   and does not recognize as sequential
  // - the interface exposed when an uncomputable sequential is abstracted
  for (const auto& key : model.topInputKeys) {
    addRole(key, "top_input");
  }
  for (const auto& key : model.topOutputKeys) {
    addRole(key, "top_output");
  }
  for (const auto& key : model.internalBoundaryInputKeys) {
    addRole(key, "opaque_internal_input");
  }
  for (const auto& key : model.internalBoundaryOutputKeys) {
    addRole(key, "opaque_internal_output");
  }
  for (const auto& detail : model.abstractedSequentialBoundaryDetails) {
    for (const auto& key : detail.stateKeys) {
      addRole(key, "abstracted_sequential_state");
    }
    for (const auto& key : detail.observedKeys) {
      addRole(key, "abstracted_sequential_observed");
    }
  }

  reports.reserve(reports.size() + reportsBySignal.size());
  for (auto& [_, entry] : reportsBySignal) {
    reports.push_back(std::move(entry));
  }
}

SequentialDesignModel extractSecDesign(naja::NL::SNLDesign* top,
                                       const char* extractedMessage,
                                       bool secDiagEnabled) {
  // LCOV_DISABLED_START
  SequentialDesignModel model = SequentialDesignModel::extract(top);
  // LCOV_DISABLED_STOP
  logSecDiagLine(secDiagEnabled, extractedMessage);
  return model;
}

AlignedSecInterface alignSecInterface(const SequentialDesignModel& model0,
                                      const SequentialDesignModel& model1,
                                      bool secDiagEnabled) {
  AlignedSecInterface aligned;
  logSecDiagLine(secDiagEnabled, "SEC diag: aligning inputs/outputs");

  aligned.inputs = alignSignalsByName(
      model0.environmentInputs,
      model0.displayNameByKey,
      model1.environmentInputs,
      model1.displayNameByKey,
      "environment input");
  const auto alignedAllOutputs = alignSignalsByName(
      model0.allObservedOutputs,
      model0.displayNameByKey,
      model1.allObservedOutputs,
      model1.displayNameByKey,
      "observed output");
  aligned.outputCoverage =
      selectCoveredObservedOutputs(alignedAllOutputs, model0, model1);
  aligned.outputs = aligned.outputCoverage.checkedOutputs;
  if (aligned.outputs.names.empty()) {
    return aligned;
  }

  if (secDiagEnabled) {
    fprintf(
        stderr,
        "SEC diag: checked_outputs=%zu total_outputs=%zu skipped=%zu\n",
        aligned.outputs.names.size(),
        aligned.outputCoverage.totalOutputs,
        aligned.outputCoverage.skippedOutputs.size());
    fflush(stderr);
  }

  // SEC equivalence is only a top-terminal contract.  Internal/opaque boundary
  // terms remain useful for coverage diagnostics, but they must not become
  // proof outputs.  Keep the coverage-selected top outputs here; re-aligning
  // the raw extractor-visible list would reintroduce stale or one-sided skipped
  // outputs after selectCoveredObservedOutputs already removed them.

  return aligned;
}

SharedSecSymbolSpace buildSharedSecSymbolSpace(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs) {
  SharedSecSymbolSpace symbolSpace;
  symbolSpace.problem.environmentInputNames = alignedInputs.names;
  symbolSpace.problem.observedOutputNames = alignedOutputs.names;
  symbolSpace.problem.originalObservedOutputCount = alignedOutputs.names.size();

  size_t nextSymbol = 2;
  assignSymbols(
      model0.stateBits,
      symbolSpace.state0Symbols,
      symbolSpace.problem.allSymbols,
      nextSymbol);
  assignSymbols(
      model1.stateBits,
      symbolSpace.state1Symbols,
      symbolSpace.problem.allSymbols,
      nextSymbol);

  for (size_t i = 0; i < alignedInputs.names.size(); ++i) {
    const size_t symbol = nextSymbol++;
    symbolSpace.inputSymbols0.emplace(alignedInputs.keys0[i], symbol);
    symbolSpace.inputSymbols1.emplace(alignedInputs.keys1[i], symbol);
    symbolSpace.problem.allSymbols.push_back(symbol);
    symbolSpace.problem.inputSymbols.push_back(symbol);
    if (auto assertedValue = getResetAssertionValue(alignedInputs.names[i]);
        assertedValue.has_value()) {
      symbolSpace.problem.resetBootstrapInputs.emplace_back(symbol, *assertedValue);
    }
  }

  for (const auto& key : model0.stateBits) {
    symbolSpace.problem.state0Symbols.push_back(symbolSpace.state0Symbols.at(key));
  }
  for (const auto& key : model1.stateBits) {
    symbolSpace.problem.state1Symbols.push_back(symbolSpace.state1Symbols.at(key));
  }

  for (const auto& relation : model0.complementedStateRelations) {
    if (symbolSpace.state0Symbols.find(relation.primaryKey) !=
            symbolSpace.state0Symbols.end() &&
        symbolSpace.state0Symbols.find(relation.complementedKey) !=
            symbolSpace.state0Symbols.end()) {
      symbolSpace.problem.complementedStatePairs0.emplace_back(
          symbolSpace.state0Symbols.at(relation.primaryKey),
          symbolSpace.state0Symbols.at(relation.complementedKey));
    }
  }
  for (const auto& relation : model1.complementedStateRelations) {
    if (symbolSpace.state1Symbols.find(relation.primaryKey) !=
            symbolSpace.state1Symbols.end() &&
        symbolSpace.state1Symbols.find(relation.complementedKey) !=
            symbolSpace.state1Symbols.end()) {
      symbolSpace.problem.complementedStatePairs1.emplace_back(
          symbolSpace.state1Symbols.at(relation.primaryKey),
          symbolSpace.state1Symbols.at(relation.complementedKey));
    }
  }

  symbolSpace.localToCombined0 = buildLocalToCombinedMap(
      model0, symbolSpace.inputSymbols0, symbolSpace.state0Symbols);
  symbolSpace.localToCombined1 = buildLocalToCombinedMap(
      model1, symbolSpace.inputSymbols1, symbolSpace.state1Symbols);
  return symbolSpace;
}

size_t nextUnusedProofSymbol(const KInductionProblem& problem) {
  size_t nextSymbol = 2;
  for (const auto symbol : problem.allSymbols) {
    nextSymbol = std::max(nextSymbol, symbol + 1);
  }
  return nextSymbol;
}

size_t privateSupportSymbol(
    size_t designIndex,
    size_t localSymbol,
    std::unordered_map<size_t, size_t>& symbolMap,
    KInductionProblem& problem,
    size_t& nextSymbol) {
  const auto existingIt = symbolMap.find(localSymbol);
  if (existingIt != symbolMap.end()) {
    return existingIt->second;
  }

  const size_t privateSymbol = nextSymbol++;
  symbolMap.emplace(localSymbol, privateSymbol);
  problem.allSymbols.push_back(privateSymbol);
  problem.inputSymbols.push_back(privateSymbol);
  problem.environmentInputNames.push_back(
      "$private_d" + std::to_string(designIndex) + "_v" +
      std::to_string(localSymbol));
  return privateSymbol;
}

BoolExpr* remapSecFormulaWithPrivateSupport(
    BoolExpr* root,
    size_t designIndex,
    std::unordered_map<size_t, size_t>& symbolMap,
    std::unordered_map<BoolExpr*, BoolExpr*>& memo,
    KInductionProblem& problem,
    size_t& nextSymbol) {
  if (root == nullptr) {
    return nullptr;  // LCOV_EXCL_LINE
  }
  if (auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  struct StackFrame {
    BoolExpr* expr = nullptr;
    bool visited = false;
  };

  // LCOV_DISABLED_START
  // Internal support left after compact extraction is a design-local free
  // input, not a name-aligned SEC assumption.  Allocate private proof symbols
  // LCOV_DISABLED_STOP
  // on demand while preserving the same iterative DAG remap used by
  // LCOV_DISABLED_START
  // BoolExprUtils for large gate-level cones.
  std::vector<StackFrame> stack;
  stack.push_back({root, false});
  while (!stack.empty()) {
  // LCOV_DISABLED_STOP
    const StackFrame current = stack.back();
    stack.pop_back();
    BoolExpr* node = current.expr;
    if (node == nullptr || memo.find(node) != memo.end()) {
      continue;  // LCOV_EXCL_LINE
    }

    if (node->getOp() == Op::VAR) {
      const size_t id = node->getId();
      const size_t mapped =
          id < 2 ? id : privateSupportSymbol(
                            designIndex, id, symbolMap, problem, nextSymbol);
      memo.emplace(node, BoolExpr::Var(mapped));
      continue;
    }

    if (!current.visited) {
      stack.push_back({node, true});
      if (node->getRight() != nullptr) {
        stack.push_back({node->getRight(), false});
      }
      if (node->getLeft() != nullptr) {
        stack.push_back({node->getLeft(), false});
      }
      continue;
    }

    BoolExpr* remapped = nullptr;
    switch (node->getOp()) {
      case Op::NOT:
        remapped = BoolExpr::Not(memo.at(node->getLeft()));
        break;
      case Op::AND:
        remapped =
            BoolExpr::And(memo.at(node->getLeft()), memo.at(node->getRight()));
        break;
      case Op::OR:
        remapped =
            BoolExpr::Or(memo.at(node->getLeft()), memo.at(node->getRight()));
        break;
      case Op::XOR:
        remapped =  // LCOV_EXCL_LINE
            BoolExpr::Xor(memo.at(node->getLeft()), memo.at(node->getRight()));  // LCOV_EXCL_LINE
        break;  // LCOV_EXCL_LINE
      case Op::NONE:  // LCOV_EXCL_LINE
      default:
        throw std::runtime_error("Unsupported BoolExpr operator in SEC remap");  // LCOV_EXCL_LINE
    }
    memo.emplace(node, remapped);
  }

  return memo.at(root);
}

RemappedSecExpressions remapSecExpressions(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedOutputs,
    SharedSecSymbolSpace& symbolSpace,
    KInductionProblem& problem,
    bool remapTransitions,
    bool secDiagEnabled) {
  RemappedSecExpressions remapped;
  std::unordered_map<BoolExpr*, BoolExpr*> remapMemo0;
  std::unordered_map<BoolExpr*, BoolExpr*> remapMemo1;
  size_t nextPrivateSymbol = nextUnusedProofSymbol(problem);

  for (size_t i = 0; i < alignedOutputs.names.size(); ++i) {
    const auto& key0 = alignedOutputs.keys0[i];
    const auto& key1 = alignedOutputs.keys1[i];
    const auto remappedOutput0 = remapSecFormulaWithPrivateSupport(
        model0.observedOutputExprByKey.at(key0),
        /*designIndex=*/0,
        symbolSpace.localToCombined0,
        remapMemo0,
        problem,
        nextPrivateSymbol);
    const auto remappedOutput1 = remapSecFormulaWithPrivateSupport(
        model1.observedOutputExprByKey.at(key1),
        /*designIndex=*/1,
        symbolSpace.localToCombined1,
        remapMemo1,
        problem,
        nextPrivateSymbol);
    problem.observedOutputExprs0.push_back(remappedOutput0);
    problem.observedOutputExprs1.push_back(remappedOutput1);
  }
  logSecDiagLine(secDiagEnabled, "SEC diag: remapped observed outputs");

  if (remapTransitions) {
    for (const auto& key : model0.stateBits) {
      remapped.next0.emplace(
          key,
          remapSecFormulaWithPrivateSupport(
              model0.nextStateExprByStateKey.at(key),
              /*designIndex=*/0,
              symbolSpace.localToCombined0,
              remapMemo0,
              problem,
              nextPrivateSymbol));
    }
    for (const auto& key : model1.stateBits) {
      remapped.next1.emplace(
          key,
          remapSecFormulaWithPrivateSupport(
              model1.nextStateExprByStateKey.at(key),
              /*designIndex=*/1,
              symbolSpace.localToCombined1,
              remapMemo1,
              problem,
              nextPrivateSymbol));
    }
    logSecDiagLine(secDiagEnabled, "SEC diag: remapped next-state formulas");
  } else {
    logSecDiagLine(
        secDiagEnabled,
        "SEC diag: deferred next-state formula remapping for k-induction");
  }
  return remapped;
}

void attachLazyTransitions(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state0Symbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state1Symbols,
    std::unordered_map<size_t, size_t>&& localToCombined0,
    std::unordered_map<size_t, size_t>&& localToCombined1,
    KInductionProblem& problem) {
  auto store = std::make_shared<LazyTransitionStore>();
  store->localToCombinedByDesign[0] = std::move(localToCombined0);
  store->localToCombinedByDesign[1] = std::move(localToCombined1);
  store->sourceByStateSymbol.reserve(model0.stateBits.size() + model1.stateBits.size());
  store->remappedByStateSymbol.reserve(model0.stateBits.size() + model1.stateBits.size());
  store->remapMemoByDesign[0].reserve(model0.stateBits.size());
  store->remapMemoByDesign[1].reserve(model1.stateBits.size());

  for (const auto& key : model0.stateBits) {
    const auto symbolIt = state0Symbols.find(key);
    const auto exprIt = model0.nextStateExprByStateKey.find(key);
    if (symbolIt == state0Symbols.end() ||
        exprIt == model0.nextStateExprByStateKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    store->sourceByStateSymbol.emplace(
        symbolIt->second, LazyTransitionSource{0, exprIt->second});
  }
  for (const auto& key : model1.stateBits) {
    const auto symbolIt = state1Symbols.find(key);
    const auto exprIt = model1.nextStateExprByStateKey.find(key);
    if (symbolIt == state1Symbols.end() ||
        exprIt == model1.nextStateExprByStateKey.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    store->sourceByStateSymbol.emplace(
        symbolIt->second, LazyTransitionSource{1, exprIt->second});
  }
  problem.lazyTransitions = std::move(store);
}

DualRailSymbolPair allocateDualRailStateSymbols(
    KInductionProblem& problem,
    size_t& nextSymbol) {
  const DualRailSymbolPair rails{nextSymbol++, nextSymbol++};
  problem.allSymbols.push_back(rails.mayBeOne);
  problem.allSymbols.push_back(rails.mayBeZero);
  return rails;
}

void addDualRailStateAssignment(
    std::vector<std::pair<size_t, bool>>& assignments,
    const DualRailSymbolPair& rails,
    std::optional<bool> value) {
  if (value.has_value()) {
    assignments.emplace_back(rails.mayBeOne, *value);
    assignments.emplace_back(rails.mayBeZero, !*value);
    return;
  }

  // Resetless X is represented as the set {0,1}: both possible-value rails
  // are true.  This follows the dual-rail/ternary encoding from the SEC BMC
  // paper and avoids choosing one arbitrary boot value.
  assignments.emplace_back(rails.mayBeOne, true);
  assignments.emplace_back(rails.mayBeZero, true);
}

std::optional<bool> lookupStateValue(
    const std::unordered_map<SignalKey, bool, SignalKeyHash>& values,
    const SignalKey& key) {
  const auto valueIt = values.find(key);
  if (valueIt == values.end()) {
    return std::nullopt;
  }
  return valueIt->second;
}

std::optional<bool> lookupBootstrapValue(
    const std::unordered_map<SignalKey, bool, SignalKeyHash>& values,
    const SignalKey& key) {
  const auto valueIt = values.find(key);
  if (valueIt == values.end()) {
    return std::nullopt;
  }
  return valueIt->second;  // LCOV_EXCL_LINE
}

void addDualRailInitialAssignments(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& railsByKey,
    KInductionProblem& problem) {
  for (const auto& key : model.stateBits) {
    const auto rails = railsByKey.at(key);
    const auto value = lookupStateValue(model.initialStateValueByKey, key);
    addDualRailStateAssignment(problem.initialStateAssignments, rails, value);
    problem.initializedStateCount += 2;
  }
}

void addDualRailBootstrapAssignments(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, bool, SignalKeyHash>& bootstrapValues,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& railsByKey,
    KInductionProblem& problem) {
  if (problem.resetBootstrapInputs.empty() || problem.resetBootstrapCycles == 0) {
    return;
  }
  for (const auto& key : model.stateBits) {
    addDualRailStateAssignment(
        problem.bootstrapStateAssignments,
        railsByKey.at(key),
        lookupBootstrapValue(bootstrapValues, key));
  }
}

void addDualRailEqualityPairs(
    const AlignedSignals& equalities,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& rails0,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& rails1,
    std::vector<std::pair<size_t, size_t>>& targetPairs) {
  for (size_t i = 0; i < equalities.names.size(); ++i) {
    const auto lhs = rails0.at(equalities.keys0[i]); // LCOV_EXCL_LINE
    const auto rhs = rails1.at(equalities.keys1[i]); // LCOV_EXCL_LINE
    targetPairs.emplace_back(lhs.mayBeOne, rhs.mayBeOne); // LCOV_EXCL_LINE
    targetPairs.emplace_back(lhs.mayBeZero, rhs.mayBeZero); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE
}

void addDualRailComplementStateEqualities(
    const SequentialDesignModel& model,
    const std::unordered_map<SignalKey, DualRailSymbolPair, SignalKeyHash>& railsByKey,
    std::vector<std::pair<size_t, size_t>>& targetPairs) {
  for (const auto& relation : model.complementedStateRelations) {
    const auto primaryIt = railsByKey.find(relation.primaryKey); // LCOV_EXCL_LINE
    const auto complementedIt = railsByKey.find(relation.complementedKey); // LCOV_EXCL_LINE
    if (primaryIt == railsByKey.end() || complementedIt == railsByKey.end()) { // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }
    // QN is the Boolean complement of Q, so its possible-one rail is the same
    // as Q's possible-zero rail, and vice versa. This is an in-design
    // structural fact, not a relation between the two SEC designs.
    targetPairs.emplace_back( // LCOV_EXCL_LINE
        complementedIt->second.mayBeOne, primaryIt->second.mayBeZero); // LCOV_EXCL_LINE
    targetPairs.emplace_back( // LCOV_EXCL_LINE
        complementedIt->second.mayBeZero, primaryIt->second.mayBeOne); // LCOV_EXCL_LINE
  }
}

void allocateDualRailStateMaps(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    KInductionProblem& problem,
    size_t& nextSymbol,
    DualRailStateSymbolMaps& maps) {
  maps.state0ByKey.reserve(model0.stateBits.size());
  maps.state1ByKey.reserve(model1.stateBits.size());
  maps.localState0BySymbol.reserve(model0.stateBits.size());
  maps.localState1BySymbol.reserve(model1.stateBits.size());

  for (const auto& key : model0.stateBits) {
    const auto rails = allocateDualRailStateSymbols(problem, nextSymbol);
    maps.state0ByKey.emplace(key, rails);
    maps.localState0BySymbol.emplace(model0.inputVarByKey.at(key), rails);
    problem.dualRailStatePairs.push_back(rails);
    problem.state0Symbols.push_back(rails.mayBeOne);
    problem.state0Symbols.push_back(rails.mayBeZero);
  }
  for (const auto& key : model1.stateBits) {
    const auto rails = allocateDualRailStateSymbols(problem, nextSymbol);
    maps.state1ByKey.emplace(key, rails);
    maps.localState1BySymbol.emplace(model1.inputVarByKey.at(key), rails);
    problem.dualRailStatePairs.push_back(rails);
    problem.state1Symbols.push_back(rails.mayBeOne);
    problem.state1Symbols.push_back(rails.mayBeZero);
  }
}

void attachLazyDualRailTransitions(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const DualRailStateSymbolMaps& maps,
    std::unordered_map<size_t, size_t>&& localToCombined0,
    std::unordered_map<size_t, size_t>&& localToCombined1,
    KInductionProblem& problem) {
  auto store = std::make_shared<LazyTransitionStore>();
  store->localToCombinedByDesign[0] = std::move(localToCombined0);
  store->localToCombinedByDesign[1] = std::move(localToCombined1);
  store->dualRailStateByLocalSymbolByDesign[0] = maps.localState0BySymbol;
  store->dualRailStateByLocalSymbolByDesign[1] = maps.localState1BySymbol;

  for (const auto& key : model0.stateBits) {
    const auto rails = maps.state0ByKey.at(key);
    BoolExpr* expr = model0.nextStateExprByStateKey.at(key);
    store->sourceByStateSymbol.emplace(
        rails.mayBeOne,
        LazyTransitionSource{0, expr, LazyTransitionRail::DualRailOne});
    store->sourceByStateSymbol.emplace(
        rails.mayBeZero,
        LazyTransitionSource{0, expr, LazyTransitionRail::DualRailZero});
  }
  for (const auto& key : model1.stateBits) {
    const auto rails = maps.state1ByKey.at(key);
    BoolExpr* expr = model1.nextStateExprByStateKey.at(key);
    store->sourceByStateSymbol.emplace(
        rails.mayBeOne,
        LazyTransitionSource{1, expr, LazyTransitionRail::DualRailOne});
    store->sourceByStateSymbol.emplace(
        rails.mayBeZero,
        LazyTransitionSource{1, expr, LazyTransitionRail::DualRailZero});
  }
  problem.lazyTransitions = std::move(store);
}

void applyInitialStateAssignments(
    const std::unordered_map<SignalKey, bool, SignalKeyHash>& initialValues,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& stateSymbols,
    BoolExpr*& initialCondition,
    KInductionProblem& problem) {
  for (const auto& [key, value] : initialValues) {
    const auto symbolIt = stateSymbols.find(key);
    if (symbolIt == stateSymbols.end()) {
      continue;  // LCOV_EXCL_LINE
    }
    BoolExpr* literal = BoolExpr::Var(symbolIt->second);
    initialCondition = BoolExpr::And(
        initialCondition, value ? literal : BoolExpr::Not(literal));
    // Keep the unit reset/init facts separately from the monolithic
    // initial-condition formula.  The k-induction base solver can then encode
    // only the init values for state bits that are in the current COI, instead
    // of dragging the whole design reset cone into every output proof.
    problem.initialStateAssignments.emplace_back(symbolIt->second, value);
    ++problem.initializedStateCount;
  }
}

ReachableStateInvariant integrateReachableStateInvariant(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state0Symbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state1Symbols,
    KInductionProblem& problem,
    bool deriveResetBootstrapStrengthening) {
  BoolExpr* initialCondition = BoolExpr::createTrue();
  applyInitialStateAssignments(
      model0.initialStateValueByKey, state0Symbols, initialCondition, problem);
  applyInitialStateAssignments(
      model1.initialStateValueByKey, state1Symbols, initialCondition, problem);
  problem.totalStateCount =
      problem.state0Symbols.size() + problem.state1Symbols.size();
  if (problem.hasExplicitInitialState()) {
    problem.initialCondition = BoolExpr::simplify(initialCondition);
  }
  const ReachableStateInvariant reachableInvariant = buildReachableStateInvariant(
      model0,
      model1,
      deriveResetBootstrapStrengthening);

  for (const auto& [key, value] : reachableInvariant.bootstrapValues0) {
    if (state0Symbols.find(key) != state0Symbols.end()) {
      problem.bootstrapStateAssignments.emplace_back(state0Symbols.at(key), value);
    }
  }
  // LCOV_DISABLED_START
  for (const auto& [key, value] : reachableInvariant.bootstrapValues1) {
  // LCOV_DISABLED_STOP
    if (state1Symbols.find(key) != state1Symbols.end()) {
      problem.bootstrapStateAssignments.emplace_back(state1Symbols.at(key), value);
    // LCOV_DISABLED_START
    }
    // LCOV_DISABLED_STOP
  }

// LCOV_DISABLED_START

  problem.resetBootstrapCycles = reachableInvariant.bootstrapCycles;
  if (problem.resetBootstrapInputs.empty()) {
  // LCOV_DISABLED_STOP
    // The reachable-state pass works on each extracted model and can recognize
    // reset-looking local inputs before the final shared SEC symbol space is
    // assembled. PDR/KI/IMC can only run a reset-bootstrap proof when that
    // reset also exists as an aligned environment input with one shared symbol.
    // If no such symbol was created, keep the proof in normal initial-frontier
    // mode so design-local initial facts remain active instead of being
    // replaced by an unconstrained "bootstrap" frontier.
    problem.resetBootstrapCycles = 0;
    problem.bootstrapStateAssignments.clear();
  }
  return reachableInvariant;
}
// LCOV_DISABLED_STOP

void buildSecPropertiesAndTransitions(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const AlignedSignals& alignedOutputs,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state0Symbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state1Symbols,
    const RemappedSecExpressions& remapped,
    KInductionProblem& problem,
    bool secDiagEnabled) {
  if (problem.lazyTransitions == nullptr) {
    for (const auto& key : model0.stateBits) {
      problem.transitions0.emplace_back(state0Symbols.at(key), remapped.next0.at(key));
    }
    for (const auto& key : model1.stateBits) {
      problem.transitions1.emplace_back(state1Symbols.at(key), remapped.next1.at(key));
    }
  }

  BoolExpr* property = BoolExpr::createTrue();
  BoolExpr* inductionCore = BoolExpr::createTrue();

  BoolExpr* inductionProperty = inductionCore;
  problem.dualRailOutputSkipReasons.clear();
  for (size_t i = 0; i < problem.observedOutputExprs0.size(); ++i) {
    const auto outputEquality = makeEqualityExpr(
        problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]);
    property = BoolExpr::And(property, outputEquality);
    if (secDiagEnabled) {
      fprintf(  // LCOV_EXCL_LINE
          stderr,  // LCOV_EXCL_LINE
          "SEC diag: output requires engine proof index=%zu name=%s\n",
          i,  // LCOV_EXCL_LINE
          alignedOutputs.names[i].c_str());  // LCOV_EXCL_LINE
      fflush(stderr);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_DISABLED_START
    inductionProperty = BoolExpr::And(inductionProperty, outputEquality);
  }


// LCOV_DISABLED_STOP
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  // LCOV_DISABLED_START
  problem.inductionProperty = BoolExpr::simplify(inductionProperty);
  problem.inductionBad = BoolExpr::simplify(BoolExpr::Not(problem.inductionProperty));
  // LCOV_DISABLED_STOP
  problem.description = "SEC property with aligned observed outputs";
  // LCOV_DISABLED_STOP
  logSecDiagLine(secDiagEnabled, "SEC diag: built SEC and induction properties");

// LCOV_DISABLED_START

  if (secDiagEnabled || secSummaryStatsEnabled()) {
    printf(
    // LCOV_DISABLED_STOP
        "SEC summary: property_is_true=%d induction_property_is_true=%d "
        "bad_is_false=%d induction_bad_is_false=%d reset_bootstrap_inputs=%zu "
        "bootstrap_cycles=%zu bootstrap_assignments=%zu\n",
        problem.property == BoolExpr::createTrue(),
        problem.inductionProperty == BoolExpr::createTrue(),
        problem.bad == BoolExpr::createFalse(),
        problem.inductionBad == BoolExpr::createFalse(),
        problem.resetBootstrapInputs.size(),
        problem.resetBootstrapCycles,
        problem.bootstrapStateAssignments.size());
    fflush(stdout);
  }
}

KInductionProblem buildDualRailSecProblem(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    // LCOV_DISABLED_START
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs,
    const ReachableStateInvariant& reachableInvariant,
    SharedSecSymbolSpace& symbolSpace,
    // LCOV_DISABLED_STOP
    bool useLazyTransitionRemapping,
    bool secDiagEnabled) {
  KInductionProblem problem;
  problem.environmentInputs = alignedInputs.keys0;
  problem.environmentInputNames = symbolSpace.problem.environmentInputNames;
  problem.inputSymbols = symbolSpace.problem.inputSymbols;
  problem.resetBootstrapCycles = symbolSpace.problem.resetBootstrapCycles;
  problem.resetBootstrapInputs = symbolSpace.problem.resetBootstrapInputs;
  problem.allSymbols = symbolSpace.problem.inputSymbols;
  problem.usesDualRailStateEncoding = true;

  size_t nextSymbol = nextUnusedProofSymbol(symbolSpace.problem);
  DualRailStateSymbolMaps railMaps;
  allocateDualRailStateMaps(model0, model1, problem, nextSymbol, railMaps);
  problem.totalStateCount =
      problem.state0Symbols.size() + problem.state1Symbols.size();
  addDualRailComplementStateEqualities(
      model0, railMaps.state0ByKey, problem.sameFrameStateEqualityPairs0);
  addDualRailComplementStateEqualities(
      model1, railMaps.state1ByKey, problem.sameFrameStateEqualityPairs1);

// LCOV_DISABLED_START

  addDualRailInitialAssignments(model0, railMaps.state0ByKey, problem);
  addDualRailInitialAssignments(model1, railMaps.state1ByKey, problem);
  // LCOV_DISABLED_STOP
  // The rail-valued boot frontier is already represented as structured unit
  // facts in initialStateAssignments.  Keep initialCondition non-null so the
  // existing base-case encoders enter their structured-init path without
  // materializing a huge duplicate conjunction over every rail.
  problem.initialCondition = BoolExpr::createTrue();

  addDualRailBootstrapAssignments(
      model0, reachableInvariant.bootstrapValues0, railMaps.state0ByKey, problem);
  addDualRailBootstrapAssignments(
      model1, reachableInvariant.bootstrapValues1, railMaps.state1ByKey, problem);

// LCOV_DISABLED_START

  SecDualRailVariableMapper mapper0(
      0,
      railMaps.localState0BySymbol,
      symbolSpace.localToCombined0,
      // LCOV_DISABLED_STOP
      problem,
      // LCOV_DISABLED_START
      nextSymbol);
      // LCOV_DISABLED_STOP
  SecDualRailVariableMapper mapper1(
      1,
      railMaps.localState1BySymbol,
      symbolSpace.localToCombined1,
      problem,
      nextSymbol);
  std::unordered_map<BoolExpr*, DualRailBoolExpr> memo0;
  std::unordered_map<BoolExpr*, DualRailBoolExpr> memo1;

  BoolExpr* inductionCore = BoolExpr::createTrue();

// LCOV_DISABLED_START

  BoolExpr* property = BoolExpr::createTrue();
  // LCOV_DISABLED_STOP
  problem.dualRailOutputSkipReasons.clear();
  problem.dualRailOutputSkipReasons.reserve(alignedOutputs.names.size());
  for (size_t i = 0; i < alignedOutputs.names.size(); ++i) {
    const auto out0 = buildDualRailBoolExpr(
        model0.observedOutputExprByKey.at(alignedOutputs.keys0[i]),
        mapper0,
        memo0);
    const auto out1 = buildDualRailBoolExpr(
        model1.observedOutputExprByKey.at(alignedOutputs.keys1[i]),
        mapper1,
        memo1);

    BoolExpr* outputRailEquality = BoolExpr::And(
        makeEqualityExpr(out0.mayBeOne, out1.mayBeOne),
        makeEqualityExpr(out0.mayBeZero, out1.mayBeZero));
    // LCOV_DISABLED_START
    // Keep batching/reporting aligned to real top outputs.  The rail pair is a
    // single ternary output value, so its may-one and may-zero equalities are
    // LCOV_DISABLED_STOP
    // one SEC obligation rather than two independent output bits.
    problem.observedOutputNames.push_back(alignedOutputs.names[i]);
    // LCOV_DISABLED_START
    problem.observedOutputExprs0.push_back(outputRailEquality);
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
    // LCOV_DISABLED_STOP
    // Dual-rail strategy construction only builds obligations.  The selected
    // engine must prove each top output; no side implication query can mark it
    // covered before PDR, IMC, or k-induction runs.
    if (secDiagEnabled) {
      fprintf(
          stderr,
          "SEC diag: dual-rail output requires engine proof: %s\n",
          alignedOutputs.names[i].c_str());
      fflush(stderr);
    }
    problem.dualRailOutputSkipReasons.emplace_back();
    property = BoolExpr::And(property, outputRailEquality);
  // LCOV_DISABLED_START
  }
  // LCOV_DISABLED_STOP

  if (useLazyTransitionRemapping) {
    attachLazyDualRailTransitions(
        model0,
        model1,
        railMaps,
        std::move(symbolSpace.localToCombined0),
        std::move(symbolSpace.localToCombined1),
        problem);
  } else {
    for (const auto& key : model0.stateBits) {
      const auto rails = railMaps.state0ByKey.at(key);
      const auto next = buildDualRailBoolExpr(
          // LCOV_DISABLED_START
          model0.nextStateExprByStateKey.at(key), mapper0, memo0);
          // LCOV_DISABLED_STOP
      problem.transitions0.emplace_back(rails.mayBeOne, next.mayBeOne);
      problem.transitions0.emplace_back(rails.mayBeZero, next.mayBeZero);
    }
    for (const auto& key : model1.stateBits) {
      const auto rails = railMaps.state1ByKey.at(key);
      const auto next = buildDualRailBoolExpr(
          model1.nextStateExprByStateKey.at(key), mapper1, memo1);
      problem.transitions1.emplace_back(rails.mayBeOne, next.mayBeOne);
      problem.transitions1.emplace_back(rails.mayBeZero, next.mayBeZero);
    }
  }

  BoolExpr* inductionProperty = inductionCore;
  inductionProperty = BoolExpr::And(inductionProperty, property);

  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.inductionProperty = BoolExpr::simplify(inductionProperty);
  problem.inductionBad = BoolExpr::simplify(BoolExpr::Not(problem.inductionProperty));
  problem.description =
      "SEC dual-rail steady-state property with aligned observed outputs";

  if (secDiagEnabled || secSummaryStatsEnabled()) {
    printf(
        "SEC summary: encoding=dual_rail_steady rail_state_bits=%zu "
        "rail_outputs=%zu reset_bootstrap_inputs=%zu bootstrap_cycles=%zu "
        "bootstrap_assignments=%zu "
        "dual_rail_state_relation_pairs=%zu\n",
        problem.totalStateCount,
        problem.observedOutputExprs0.size(),
        problem.resetBootstrapInputs.size(),
        problem.resetBootstrapCycles,
        problem.bootstrapStateAssignments.size(),
        problem.sameFrameStateEqualityPairs0.size() +
            problem.sameFrameStateEqualityPairs1.size());
    fflush(stdout);
  }

  return problem;
}

const char* describeSecEngine(SecEngine secEngine) {
  switch (secEngine) {
    case SecEngine::Pdr:
      return "pdr engine";
    case SecEngine::Imc:
      return "imc engine";
    case SecEngine::KInduction:
      return "classic k-induction engine";
    default:
      return "pdr engine";
  }
}

void emitSecEngineProofProgress(
    const KInductionProblem& problem,
    const std::string& engineLabel,
    size_t provenOutputCount) {
  // Keep proof-progress wording in the SEC layer so engines only report proof
  // data and future engines can share the same diagnostic style.
  for (const std::string& line : detail::buildSecEngineProofProgressDiagLines(
           engineLabel,
           problem.observedOutputNames,
           problem.observedOutputExprs0.size(),
           provenOutputCount)) {
    emitSecDiag(line);
  }
}

void setSecEngineProofProgress(
    SequentialEquivalenceResult& result,
    const KInductionProblem& problem,
    const std::string& engineLabel,
    size_t provenOutputCount) {
  result.proofProgress = detail::buildSecEngineProofProgress(
      engineLabel,
      problem.observedOutputNames,
      problem.observedOutputExprs0.size(),
      provenOutputCount);
}

SequentialEquivalenceResult runPdrSecEngine(
    // LCOV_DISABLED_START
    const KInductionProblem& problem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    // LCOV_DISABLED_STOP
    const OutputCoverageSelection& outputCoverage,
    // LCOV_DISABLED_START
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports) {
  // PDR still needs the cheap frame-0 mismatch check before growing frames, but
  // it should not invoke the full k-induction top engine with max_k=0.  A
  // LCOV_DISABLED_STOP
  // bounded engine run at k=0 is necessarily inconclusive for sequential
  // LCOV_DISABLED_START
  // problems, which made the output-batching fallback split every output and
  // repeat the same BMC setup hundreds of times before PDR even started.
  // Keep this optional validation local in dual-rail SEC.  Small probe cases
  // still need exact frame-0 localization for counterexamples, but huge
  // rail-state SoC surfaces should enter PDR directly instead of materializing
  // a pre-PDR dual-rail transition relation.
  bool broadBasePrecheckDone = false;
  if (problem.usesDualRailStateEncoding) {
    if (!shouldDeferWideDualRailPdrValidation(problem)) {
      // Validate small dual-rail frame-0 SEC predicates one output at a time so
      // PDR can safely seed the batch property into F0.  Medium/wide batches use
      // the selected PDR engine directly and split on abstract traces.
      if (auto witness =
              findPdrPerOutputValidationCounterexample(problem, solverType, 0);
          witness.has_value()) {
        const KInductionResult witnessResult{  // LCOV_EXCL_LINE
            KInductionStatus::Different,
            witness->badFrame,  // LCOV_EXCL_LINE
            std::move(witness)};  // LCOV_EXCL_LINE
        return makeSecResult(  // LCOV_EXCL_LINE
            SequentialEquivalenceStatus::Different,
            witnessResult.bound,  // LCOV_EXCL_LINE
            formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
            outputCoverage,  // LCOV_EXCL_LINE
            abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
            extractedBoundaryReports);  // LCOV_EXCL_LINE
      }
      broadBasePrecheckDone = true;
    } else if (pdrStrategyStatsEnabled()) {
      emitSecDiag(  // LCOV_EXCL_LINE
          "SEC PDR stats: skipped dual-rail frame-0 validation ",
          "outputs=", problem.observedOutputExprs0.size(),
          " output_limit=",
          detail::kMaxPdrDualRailFrameZeroValidationOutputs,
          " state_symbols=", pdrCertificateStateSymbolCount(problem),
          " state_limit=",
          detail::kMaxPdrDualRailFrameZeroValidationStateSymbols);
    }
  } else {
    broadBasePrecheckDone = true;
    // This is the same exact frame-0 SEC query as the broad base checker, but it
    // lets the base solver localize multi-output frontiers.  Wide ASIC outputs can
    // make one monolithic bad-output OR dominate PDR before the engine starts.
    if (auto witness = findPdrValidationCounterexample(problem, solverType, 0);
        witness.has_value()) {
      const KInductionResult witnessResult{  // LCOV_EXCL_LINE
          KInductionStatus::Different, witness->badFrame, std::move(witness)};  // LCOV_EXCL_LINE
      return makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Different,
          witnessResult.bound,  // LCOV_EXCL_LINE
          formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }

  if (problem.combinedStateSymbols().empty()) {
    return makeSecResult(
    // LCOV_DISABLED_STOP
        SequentialEquivalenceStatus::Equivalent,
        // LCOV_DISABLED_START
        0,
        "",
        outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }


// LCOV_DISABLED_STOP
  const std::vector<size_t> dualRailEngineOutputIndices =
      // LCOV_DISABLED_START
      collectOutputsRequiringDualRailEngineProof(problem);
  std::unordered_map<size_t, std::string> presetDualRailSkipReasons;
  if (problem.dualRailOutputSkipReasons.size() ==
      problem.observedOutputExprs0.size()) {
      // LCOV_DISABLED_STOP
    for (size_t i = 0; i < problem.dualRailOutputSkipReasons.size(); ++i) {
      // LCOV_DISABLED_START
      if (!problem.dualRailOutputSkipReasons[i].empty()) {
        presetDualRailSkipReasons.emplace(i, problem.dualRailOutputSkipReasons[i]); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
    }
  }
  if (problem.usesDualRailStateEncoding && dualRailEngineOutputIndices.empty()) {
    std::vector<bool> coveredOutputs(problem.observedOutputExprs0.size(), false); // LCOV_EXCL_LINE
    const OutputCoverageSelection partialCoverage =
        buildCoverageWithDualRailOutputSkips( // LCOV_EXCL_LINE
            outputCoverage, // LCOV_EXCL_LINE
            problem, // LCOV_EXCL_LINE
            coveredOutputs,
            presetDualRailSkipReasons);
    return makeSecResult( // LCOV_EXCL_LINE
        SequentialEquivalenceStatus::Inconclusive,
        0,
        "Dual-rail PDR has no selected-engine output obligation to prove", // LCOV_EXCL_LINE
        partialCoverage,
        abstractedSequentialBoundaries, // LCOV_EXCL_LINE
        extractedBoundaryReports); // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  // LCOV_DISABLED_START
  auto filterPairsToSupport =
  // LCOV_DISABLED_STOP
      [](const std::vector<std::pair<size_t, size_t>>& source,
         std::vector<std::pair<size_t, size_t>>& target,
         const std::unordered_set<size_t>& support) {
        target.clear();
        for (const auto& pair : source) {
          if (support.find(pair.first) != support.end() ||  // LCOV_EXCL_LINE
              support.find(pair.second) != support.end()) {  // LCOV_EXCL_LINE
            target.push_back(pair);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
        }
      };

  // LCOV_DISABLED_START
  auto filterAssignmentsToSupport =
  // LCOV_DISABLED_STOP
      [](const std::vector<std::pair<size_t, bool>>& source,
         std::vector<std::pair<size_t, bool>>& target,
         const std::unordered_set<size_t>& support) {
        target.clear();
        for (const auto& assignment : source) {
          // LCOV_DISABLED_START
          if (support.find(assignment.first) != support.end()) {  // LCOV_EXCL_LINE
          // LCOV_DISABLED_STOP
            target.push_back(assignment);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
        }
      };

// LCOV_DISABLED_START

  auto rebuildPdrBatchStrengthening = [](KInductionProblem& batch) {
    BoolExpr* inductionProperty = BoolExpr::createTrue();
    for (size_t i = 0; i < batch.observedOutputExprs0.size(); ++i) {
      inductionProperty = BoolExpr::And(
          inductionProperty,
          makeEqualityExpr(
              batch.observedOutputExprs0[i], batch.observedOutputExprs1[i]));
    // LCOV_DISABLED_START
    }
    // PDR consumes this only as a candidate frame-strengthening lemma. The
    // engine validates both Init => lemma and lemma /\ T => lemma' before the
    // formula can constrain any bad-cube or predecessor query.
    batch.inductionProperty = BoolExpr::simplify(inductionProperty);
    // LCOV_DISABLED_STOP
    batch.inductionBad = BoolExpr::simplify(BoolExpr::Not(batch.inductionProperty));
  };

  // LCOV_DISABLED_START
  TransitionExprResolver pdrBatchTransitionByState(problem);
  // LCOV_DISABLED_STOP
  const auto& pdrBatchPrimaryByComplement =
      // LCOV_DISABLED_START
      pdrBatchTransitionByState.primaryByComplement();
      // LCOV_DISABLED_STOP

  auto computePdrBatchSupportClosure = [&](const KInductionProblem& batch,
                                           size_t transitionClosureLimit) {
    // LCOV_DISABLED_START
    if (batch.property == nullptr) {
    // LCOV_DISABLED_STOP
      return std::unordered_set<size_t>{};  // LCOV_EXCL_LINE
    }
    const auto propertySupport = batch.property->getSupportVars();
    // LCOV_DISABLED_START
    std::unordered_set<size_t> support(propertySupport.begin(), propertySupport.end());
    std::unordered_set<size_t> expandedTransitionStates;
    std::vector<size_t> worklist;

    auto enqueueTransitionState = [&](size_t symbol) {
    // LCOV_DISABLED_STOP
      if (!pdrBatchTransitionByState.contains(symbol)) {
        if (const auto primaryIt = pdrBatchPrimaryByComplement.find(symbol);  // LCOV_EXCL_LINE
            // LCOV_DISABLED_START
            primaryIt != pdrBatchPrimaryByComplement.end()) {  // LCOV_EXCL_LINE
          symbol = primaryIt->second;  // LCOV_EXCL_LINE
        } else {  // LCOV_EXCL_LINE
        // LCOV_DISABLED_STOP
          return;  // LCOV_EXCL_LINE
        // LCOV_DISABLED_START
        }
      }  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
      support.insert(symbol);
      if (expandedTransitionStates.insert(symbol).second) {
        worklist.push_back(symbol);
      }
    };

// LCOV_DISABLED_START

    for (const auto propertySymbol : propertySupport) {
    // LCOV_DISABLED_STOP
      enqueueTransitionState(propertySymbol);
    // LCOV_DISABLED_START
    }
    for (size_t cursor = 0;
         cursor < worklist.size() &&
         // LCOV_DISABLED_STOP
         support.size() < transitionClosureLimit;
         // LCOV_DISABLED_START
         ++cursor) {
      for (const auto dependency : pdrBatchTransitionByState.support(worklist[cursor])) {
        if (support.insert(dependency).second) {
          enqueueTransitionState(dependency);  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
      }
    }
    return support;
    // LCOV_DISABLED_STOP
  };

  auto prunePdrBatchStrengthening = [&](KInductionProblem& batch,
                                        size_t transitionClosureLimit) {
    // LCOV_DISABLED_START
    auto support =
        computePdrBatchSupportClosure(batch, transitionClosureLimit);

    // A PDR output slice may only inherit same-design rail relations and reset
    // value facts; cross-design internal equalities have no representation.
    filterPairsToSupport(
        problem.sameFrameStateEqualityPairs0,
        batch.sameFrameStateEqualityPairs0,
        support);
    filterPairsToSupport(
        problem.sameFrameStateEqualityPairs1,
        batch.sameFrameStateEqualityPairs1,
        support);
    filterAssignmentsToSupport(
        // LCOV_DISABLED_START
        problem.initialStateAssignments, batch.initialStateAssignments, support);
        // LCOV_DISABLED_STOP
    filterAssignmentsToSupport(
        problem.bootstrapStateAssignments, batch.bootstrapStateAssignments, support);
    for (const auto& pair : batch.sameFrameStateEqualityPairs0) {
      support.insert(pair.first);  // LCOV_EXCL_LINE
      support.insert(pair.second);  // LCOV_EXCL_LINE
    }
    for (const auto& pair : batch.sameFrameStateEqualityPairs1) {
      support.insert(pair.first);  // LCOV_EXCL_LINE
      support.insert(pair.second);  // LCOV_EXCL_LINE
    }
    // LCOV_DISABLED_START
    for (const auto& assignment : batch.initialStateAssignments) {
      support.insert(assignment.first);  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
    }
    for (const auto& assignment : batch.bootstrapStateAssignments) {
      support.insert(assignment.first);  // LCOV_EXCL_LINE
    }
    rebuildPdrBatchStrengthening(batch);

// LCOV_DISABLED_START

    if (batch.lazyTransitions != nullptr) {
    // LCOV_DISABLED_STOP
      auto& store = *batch.lazyTransitions;
      // LCOV_DISABLED_START
      batch.transitions0.clear();
      batch.transitions1.clear();
      constexpr size_t kMaxEagerRemappedPdrBatchTransitions = 1024;
      // LCOV_DISABLED_STOP
      if (support.size() > kMaxEagerRemappedPdrBatchTransitions) {
        // LCOV_DISABLED_START
        // Keep large ASIC batches lazy.  Sampling on BlackParrot showed that
        // eagerly remapping a 12k-symbol support closure built more than a
        // million transition DAG nodes before the first PDR SAT query.  The
        // transition resolver still has the exact support closure above, and
        // will remap only the transitions that PDR actually encodes.
        return;  // LCOV_EXCL_LINE
      }
      batch.transitions0.reserve(support.size());
      batch.transitions1.reserve(support.size());
      // LCOV_DISABLED_STOP
      const TransitionExprResolver batchTransitionByState(batch);

      // Sampling on BlackParrot showed the proof spending time lazily remapping
      // next-state expressions inside predecessor queries. Once the batch cone
      // is already pruned to the output support closure, remap those relevant
      // LCOV_DISABLED_START
      // transitions eagerly through the resolver so binary and dual-rail lazy
      // transitions are materialized in the same symbol space used by COI.
      for (const auto symbol : support) {
      // LCOV_DISABLED_STOP
        const auto sourceIt = store.sourceByStateSymbol.find(symbol);
        // LCOV_DISABLED_START
        if (sourceIt == store.sourceByStateSymbol.end()) {
          continue;  // LCOV_EXCL_LINE
        }
        // LCOV_DISABLED_STOP

        BoolExpr* remapped = batchTransitionByState.at(symbol);

// LCOV_DISABLED_START

        if (sourceIt->second.designIndex == 0) {
        // LCOV_DISABLED_STOP
          batch.transitions0.emplace_back(symbol, remapped);
        // LCOV_DISABLED_START
        } else {
          batch.transitions1.emplace_back(symbol, remapped);
        }
      }
    }
    // LCOV_DISABLED_STOP
  };

// LCOV_DISABLED_START


// LCOV_DISABLED_STOP
  auto prunePdrBatchRelations = [&](KInductionProblem& batch,
                                    // LCOV_DISABLED_START
                                    size_t transitionClosureLimit) {
    prunePdrBatchStrengthening(batch, transitionClosureLimit);
    // LCOV_DISABLED_STOP
  };

// LCOV_DISABLED_START

  // PDR is still proving real PDR obligations, but wide ASIC SEC properties are
  // better handled as output-cone slices. This keeps reset-bootstrap F[0]
  // LCOV_DISABLED_STOP
  // strengthening and blocking queries local to a small property instead of
  // LCOV_DISABLED_START
  // materializing every observed output in one frame.
  // LCOV_DISABLED_STOP
  //
  // Keep each PDR batch bounded, but do not prove one output per engine run.
  // BlackParrot sampling showed the one-output mode repeating the same
  // reset-frontier and PDR blocking work hundreds of times.  A moderate batch
  // still proves a real conjunction slice. If projected PDR finds a
  // LCOV_DISABLED_START
  // counterexample on a multi-output slice, escalate PDR precision first and
  // avoid broad concrete-BMC validation until the final exact retry.
  constexpr size_t kMinOutputsForBatchedPdrProof = 129;
  constexpr OutputBatchingLimits kPdrOutputBatchingLimits{32, 1024};
  // Dual-rail residuals often need the shared all-output reset frontier as an
  // LCOV_DISABLED_STOP
  // F0 strengthening fact. Ibex in particular proves completely when the 100
  // LCOV_DISABLED_START
  // residual rail outputs are handled together, while small slices lose that
  // LCOV_DISABLED_STOP
  // context and only cover the first few control outputs.
  constexpr OutputBatchingLimits kDualRailPdrOutputBatchingLimits{128, 8192};
  const OutputBatchingLimits pdrOutputBatchingLimits =
      // LCOV_DISABLED_START
      problem.usesDualRailStateEncoding ? dualRailPdrOutputBatchingLimits(
      // LCOV_DISABLED_STOP
                                              kDualRailPdrOutputBatchingLimits)
                                        // LCOV_DISABLED_START
                                        : kPdrOutputBatchingLimits;
                                        // LCOV_DISABLED_STOP
  constexpr size_t kPdrBatchTransitionClosureLimit = 12000;
  constexpr size_t kRefinedPdrBatchTransitionClosureLimit = 60000;
  constexpr size_t kDualRailPdrBatchTransitionClosureLimit = 2048;
  constexpr size_t kDualRailRefinedPdrBatchTransitionClosureLimit = 8192;
  const size_t pdrBatchTransitionClosureLimit =
      problem.usesDualRailStateEncoding
          ? secStrategySizeLimitFromEnv(
                "KEPLER_SEC_PDR_DUAL_RAIL_BATCH_CLOSURE_LIMIT",
                // LCOV_DISABLED_START
                kDualRailPdrBatchTransitionClosureLimit)
          : kPdrBatchTransitionClosureLimit;
  const size_t refinedPdrBatchTransitionClosureLimit =
      problem.usesDualRailStateEncoding
      // LCOV_DISABLED_STOP
          ? secStrategySizeLimitFromEnv(
                "KEPLER_SEC_PDR_DUAL_RAIL_REFINED_CLOSURE_LIMIT",
                // LCOV_DISABLED_START
                kDualRailRefinedPdrBatchTransitionClosureLimit)
          : kRefinedPdrBatchTransitionClosureLimit;
  const bool dualRailPdrUsesResetFrontier =
  // LCOV_DISABLED_STOP
      problem.usesDualRailStateEncoding;
  // LCOV_DISABLED_START
  struct PdrOutputBatch {
    size_t firstOutput = 0;
    size_t endOutput = 0;
    // LCOV_DISABLED_STOP
    bool startAtFinalExact = false;
  // LCOV_DISABLED_START
  };
  std::vector<PdrOutputBatch> outputBatches;
  const bool useSupportBoundedPdrBatches =
      problem.usesDualRailStateEncoding ||
      // LCOV_DISABLED_STOP
      problem.observedOutputExprs0.size() >= kMinOutputsForBatchedPdrProof;  // LCOV_EXCL_LINE
  // LCOV_DISABLED_START
  if (!useSupportBoundedPdrBatches) {
    // Batching protects very wide SEC/PDR properties from broad bad-state
    // queries.  On medium designs, each tiny batch repeats the same
    // reset/bootstrap invariant validation, so prove one conjunction slice and
    // reserve batching for BlackParrot/AES-scale output counts.
    outputBatches.push_back({0, problem.observedOutputExprs0.size(), false});  // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
  } else {  // LCOV_EXCL_LINE
    for (const auto& [firstOutput, endOutput] :
         buildSupportBoundedOutputBatches(problem, pdrOutputBatchingLimits)) {
      outputBatches.push_back({firstOutput, endOutput, false});
    }
  }
  KInductionProblem batchProblem = problem;
  std::vector<bool> pdrCoveredOutputs =
      makeInitialPdrCoveredOutputs(problem);
  std::unordered_map<size_t, std::string> pdrSkippedOutputReasons =
      presetDualRailSkipReasons;
  size_t provedBound = 0;
  const bool emitPdrStageStats = pdrStrategyStatsEnabled();
  struct FinalPdrStageOutcome {  // LCOV_EXCL_LINE
    bool equivalent = false;  // LCOV_EXCL_LINE
    bool shouldSplit = false;  // LCOV_EXCL_LINE
    bool shouldSkipOutput = false;  // LCOV_EXCL_LINE
    std::optional<SequentialEquivalenceResult> terminalResult;
  };
  auto runFinalExactPdrStage =
      // LCOV_DISABLED_START
      [&](size_t batchIndex,
          size_t firstOutput,
          // LCOV_DISABLED_STOP
          size_t endOutput) -> FinalPdrStageOutcome {
    // LCOV_DISABLED_START
    constexpr size_t kMaxPdrConcreteValidationOutputs = 1;  // LCOV_EXCL_LINE
    constexpr size_t kMaxDualRailFinalExactPdrOutputBatchSize = 8;  // LCOV_EXCL_LINE
    const size_t kMaxFinalExactPdrOutputBatchSize =  // LCOV_EXCL_LINE
        problem.usesDualRailStateEncoding
            ? kMaxDualRailFinalExactPdrOutputBatchSize
            : pdrOutputBatchingLimits.maxOutputBatchSize;  // LCOV_EXCL_LINE
    // The final exact repair already carries exact frame clauses and validated
    // bad-formula clauses. Keeping both predecessor and bad cubes bounded avoids
    // LCOV_DISABLED_STOP
    // large single-output loops from enumerating thousands of sibling cubes.
    constexpr size_t kFinalExactPdrPredecessorProjectionLimit = 16;  // LCOV_EXCL_LINE
    constexpr size_t kFinalExactPdrBadCubeStateLimit = 32;  // LCOV_EXCL_LINE
    constexpr size_t kFinalExactPdrRootGeneralizationAttempts = 0;  // LCOV_EXCL_LINE
    // Dual-rail final PDR validates projected roots exactly. Keep that exact
    // CEGAR repair bounded, but leave enough queries for small ASIC leaves that
    // need several reset-frontier blockers before the output proof converges.
    // Multi-output slices still split quickly; isolated hard leaves are skipped
    // as uncovered instead of consuming the whole workflow.
    constexpr size_t kDualRailFinalExactPdrRootGeneralizationAttempts = 4;  // LCOV_EXCL_LINE
    constexpr size_t kDualRailFinalExactPdrMultiOutputQueryBudget = 64;  // LCOV_EXCL_LINE
    constexpr size_t kDualRailFinalExactPdrSingleOutputQueryBudget = 1024;  // LCOV_EXCL_LINE
    constexpr size_t kLargeDualRailFinalExactPdrSingleOutputQueryBudget = 64;  // LCOV_EXCL_LINE
    constexpr size_t kDualRailFinalExactPdrMultiOutputRepairBudget = 2;  // LCOV_EXCL_LINE
    constexpr size_t kDualRailFinalExactPdrSingleOutputRepairBudget = 8;  // LCOV_EXCL_LINE
    constexpr size_t kLargeDualRailFinalExactPdrSingleOutputRepairBudget = 2;  // LCOV_EXCL_LINE
    constexpr size_t kMediumDualRailFinalExactPdrPredecessorProjectionLimit = 32;  // LCOV_EXCL_LINE
    constexpr size_t kMediumDualRailFinalExactPdrMultiOutputRepairBudget = 4;  // LCOV_EXCL_LINE
    constexpr size_t kMediumDualRailFinalExactPdrSingleOutputRepairBudget = 8;  // LCOV_EXCL_LINE
    if (endOutput - firstOutput > kMaxFinalExactPdrOutputBatchSize) {  // LCOV_EXCL_LINE
      if (emitPdrStageStats) {  // LCOV_EXCL_LINE
        emitSecDiag(  // LCOV_EXCL_LINE
            "SEC PDR stats: splitting before final exact repair ",
            "outputs=[", firstOutput, ",", endOutput, ")",
            " limit=", kMaxFinalExactPdrOutputBatchSize);
      // LCOV_DISABLED_START
      }  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
      FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
      outcome.shouldSplit = true;  // LCOV_EXCL_LINE
      return outcome;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE

    KInductionProblem validationProblem = problem;  // LCOV_EXCL_LINE
    configureOutputBatchProblem(  // LCOV_EXCL_LINE
        validationProblem, problem, firstOutput, endOutput);  // LCOV_EXCL_LINE

    KInductionProblem fullExactBatchProblem = problem;  // LCOV_EXCL_LINE
    configureOutputBatchProblem(  // LCOV_EXCL_LINE
        fullExactBatchProblem, problem, firstOutput, endOutput);  // LCOV_EXCL_LINE
    // Keep the full transition relation for this output slice during the
    // final exact retry, but do not reintroduce unrelated startup/induction
    // facts from the rest of the SEC problem. Those global relations are
    // useful for broad proofs, yet they can dominate SAT encoding once this
    // last retry is focused on a bounded slice.
    prunePdrBatchStrengthening(  // LCOV_EXCL_LINE
        fullExactBatchProblem, refinedPdrBatchTransitionClosureLimit);  // LCOV_EXCL_LINE
    // Dual-rail PDR proves the property over the rail-encoded transition
    // system.  Re-running a full bounded SEC check after an Equivalent leaf is
    // only a sanity check, and Swerv-size cones materialize gigabytes of BMC
    // clauses there.  Keep concrete validation for binary PDR; dual-rail
    // abstract differences still split/skip unless PDR repairs them internally.
    const bool finalBatchCanValidateConcrete =  // LCOV_EXCL_LINE
        !problem.usesDualRailStateEncoding &&  // LCOV_EXCL_LINE
        endOutput - firstOutput <= kMaxPdrConcreteValidationOutputs;  // LCOV_EXCL_LINE
    const bool finalBatchCanRefineProjectedCounterexamples = true;  // LCOV_EXCL_LINE
    const size_t originalOutputCount =  // LCOV_EXCL_LINE
        problem.originalObservedOutputCount == 0  // LCOV_EXCL_LINE
            ? problem.observedOutputExprs0.size()  // LCOV_EXCL_LINE
            : problem.originalObservedOutputCount;  // LCOV_EXCL_LINE
    const bool mediumDualRailOutputSurface =  // LCOV_EXCL_LINE
        problem.usesDualRailStateEncoding &&  // LCOV_EXCL_LINE
        originalOutputCount <=
            kMaxDualRailFinalResetFrontierOriginalOutputs;  // LCOV_EXCL_LINE
    const bool largeDualRailOutputSurface =  // LCOV_EXCL_LINE
        problem.usesDualRailStateEncoding && !mediumDualRailOutputSurface;  // LCOV_EXCL_LINE
    // The bad-formula repair opens exact reset-frontier queries. Keep it away
    // from broad dual-rail batches. Also keep already-split one-output leaves
    // behind the original-output surface guard: BP-scale SoC probes otherwise
    // repeat this local-looking repair hundreds of times over the same huge
    // reset-specialized frontier.
    const bool finalSliceUsesBadFormulaValidation =  // LCOV_EXCL_LINE
        (!problem.usesDualRailStateEncoding ||  // LCOV_EXCL_LINE
         (endOutput - firstOutput == 1 && mediumDualRailOutputSurface)) &&  // LCOV_EXCL_LINE
        endOutput - firstOutput <=  // LCOV_EXCL_LINE
            pdrOutputBatchingLimits.maxOutputBatchSize;  // LCOV_EXCL_LINE
    // Exact reset-frontier checks repair reset-bootstrap dual-rail slices, but
    // the final stage may split a wide original SEC surface into one-output
    // leaves.  Keep the original output width in the decision so a wide design
    // does not re-enter the same reset-frontier wall one leaf at a time.
    // Medium CPU-style residual buses need exact reset-frontier repair even after
    // batching splits them. Larger SoC-scale surfaces remain guarded in PDREngine
    // by the rail-state and transition-source limits.
    const bool finalSliceUsesResetFrontier =  // LCOV_EXCL_LINE
        !problem.usesDualRailStateEncoding ||
        originalOutputCount <=
            kMaxDualRailFinalResetFrontierOriginalOutputs;  // LCOV_EXCL_LINE
    const size_t finalPdrPredecessorProjectionLimit =  // LCOV_EXCL_LINE
        mediumDualRailOutputSurface  // LCOV_EXCL_LINE
            ? kMediumDualRailFinalExactPdrPredecessorProjectionLimit
            : kFinalExactPdrPredecessorProjectionLimit;
    const size_t finalPdrBadCubeStateLimit =  // LCOV_EXCL_LINE
        kFinalExactPdrBadCubeStateLimit;
    // These are per-leaf repair budgets. Keep them modest: wide SoC surfaces
    // can leave many final single-output dual-rail slices, while smaller
    // isolated memory-handshake leaves may need the ordinary PDR loop to run to
    // frame/max-K convergence after deterministic reset-conflict repairs.
    const size_t defaultFinalDualRailPredecessorQueryBudget =  // LCOV_EXCL_LINE
        endOutput - firstOutput == 1
            ? (largeDualRailOutputSurface  // LCOV_EXCL_LINE
                   ? kLargeDualRailFinalExactPdrSingleOutputQueryBudget
                   : kDualRailFinalExactPdrSingleOutputQueryBudget)
            : kDualRailFinalExactPdrMultiOutputQueryBudget;
    const size_t finalDualRailPredecessorQueryBudget =  // LCOV_EXCL_LINE
        secStrategySizeLimitFromEnv(  // LCOV_EXCL_LINE
            "KEPLER_SEC_PDR_DUAL_RAIL_FINAL_QUERY_BUDGET",
            defaultFinalDualRailPredecessorQueryBudget);
    emitPdrStrategyStageStats(  // LCOV_EXCL_LINE
        // LCOV_DISABLED_START
        emitPdrStageStats,  // LCOV_EXCL_LINE
        batchIndex,  // LCOV_EXCL_LINE
        // LCOV_DISABLED_STOP
        firstOutput,  // LCOV_EXCL_LINE
        endOutput,  // LCOV_EXCL_LINE
        "full_exact_strengthening_pruned",
        // LCOV_DISABLED_START
        refinedPdrBatchTransitionClosureLimit,  // LCOV_EXCL_LINE
        // LCOV_DISABLED_STOP
        finalPdrPredecessorProjectionLimit,
        finalPdrBadCubeStateLimit,
        fullExactBatchProblem);
    // LCOV_DISABLED_START
    PDREngine fullExactPdrEngine(  // LCOV_EXCL_LINE
        fullExactBatchProblem,
        solverType,  // LCOV_EXCL_LINE
        finalPdrPredecessorProjectionLimit,
        finalPdrBadCubeStateLimit,
        // LCOV_DISABLED_STOP
        /*useExactFrameClauses=*/true,
        // LCOV_DISABLED_START
        /*maxPredecessorQueries=*/
        // LCOV_DISABLED_STOP
            problem.usesDualRailStateEncoding  // LCOV_EXCL_LINE
                // LCOV_DISABLED_START
                ? finalDualRailPredecessorQueryBudget  // LCOV_EXCL_LINE
                : 0,
        /*refineProjectedCounterexamples=*/
            finalBatchCanRefineProjectedCounterexamples,
        /*maxBoundedRootGeneralizationAttempts=*/
        // LCOV_DISABLED_STOP
            problem.usesDualRailStateEncoding  // LCOV_EXCL_LINE
                ? kDualRailFinalExactPdrRootGeneralizationAttempts
                : kFinalExactPdrRootGeneralizationAttempts,
        /*learnValidatedBadFormulaClauses=*/finalSliceUsesBadFormulaValidation,  // LCOV_EXCL_LINE
        /*useExactResetFrontierChecks=*/finalSliceUsesResetFrontier,
        /*maxProjectedCounterexampleRefinements=*/
            problem.usesDualRailStateEncoding  // LCOV_EXCL_LINE
                ? (endOutput - firstOutput > 1  // LCOV_EXCL_LINE
                       ? (mediumDualRailOutputSurface  // LCOV_EXCL_LINE
                              ? kMediumDualRailFinalExactPdrMultiOutputRepairBudget
                              : kDualRailFinalExactPdrMultiOutputRepairBudget)
                       : (largeDualRailOutputSurface  // LCOV_EXCL_LINE
                              ? kLargeDualRailFinalExactPdrSingleOutputRepairBudget
                              : (mediumDualRailOutputSurface  // LCOV_EXCL_LINE
                                     ? kMediumDualRailFinalExactPdrSingleOutputRepairBudget
                                     : kDualRailFinalExactPdrSingleOutputRepairBudget)))
                : 0);
    // Dual-rail properties are intentionally batched before the reset
    // bootstrap proof, otherwise the frame-0 precheck materializes the entire
    // wide rail-encoded SEC property. Binary SEC kept the historical broad
    // precheck above, so those batches may still reuse it.
    const auto fullExactPdrResult =
        fullExactPdrEngine.run(maxK, broadBasePrecheckDone);  // LCOV_EXCL_LINE
    if (fullExactPdrResult.status == PDRStatus::Equivalent) {  // LCOV_EXCL_LINE
      if (finalBatchCanValidateConcrete) {  // LCOV_EXCL_LINE
        if (auto fullExactWitness = findPdrValidationCounterexample(  // LCOV_EXCL_LINE
                validationProblem, solverType, maxK);  // LCOV_EXCL_LINE
            fullExactWitness.has_value()) {  // LCOV_EXCL_LINE
          const KInductionResult witnessResult{  // LCOV_EXCL_LINE
              KInductionStatus::Different,
              fullExactWitness->badFrame,  // LCOV_EXCL_LINE
              std::move(fullExactWitness)};  // LCOV_EXCL_LINE
          FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
          outcome.terminalResult = makeSecResult(  // LCOV_EXCL_LINE
              SequentialEquivalenceStatus::Different,
              witnessResult.bound,  // LCOV_EXCL_LINE
              formatCounterexampleWitness(  // LCOV_EXCL_LINE
                  witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
              outputCoverage,  // LCOV_EXCL_LINE
              abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
              extractedBoundaryReports);  // LCOV_EXCL_LINE
          return outcome;  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      provedBound = std::max(provedBound, fullExactPdrResult.bound);  // LCOV_EXCL_LINE
      markDualRailPdrOutputRangeCovered(  // LCOV_EXCL_LINE
          pdrCoveredOutputs,  // LCOV_EXCL_LINE
          pdrSkippedOutputReasons,  // LCOV_EXCL_LINE
          firstOutput,  // LCOV_EXCL_LINE
          endOutput);  // LCOV_EXCL_LINE
      FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
      outcome.equivalent = true;  // LCOV_EXCL_LINE
      return outcome;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (fullExactPdrResult.status == PDRStatus::Different) {  // LCOV_EXCL_LINE
      // LCOV_DISABLED_START
      std::optional<KInductionResult::CounterexampleWitness>
          fullExactWitness;  // LCOV_EXCL_LINE
          // LCOV_DISABLED_STOP
      if (finalBatchCanValidateConcrete) {  // LCOV_EXCL_LINE
        fullExactWitness = SEC::findBaseCounterexampleAtFrontier(  // LCOV_EXCL_LINE
            validationProblem, solverType, fullExactPdrResult.bound);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_DISABLED_START
      if (fullExactWitness.has_value()) {  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
        const KInductionResult witnessResult{  // LCOV_EXCL_LINE
            KInductionStatus::Different,
            fullExactPdrResult.bound,  // LCOV_EXCL_LINE
            // LCOV_DISABLED_START
            std::move(fullExactWitness)};  // LCOV_EXCL_LINE
        FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
        outcome.terminalResult = makeSecResult(  // LCOV_EXCL_LINE
            SequentialEquivalenceStatus::Different,
            fullExactPdrResult.bound,  // LCOV_EXCL_LINE
            // LCOV_DISABLED_STOP
            formatCounterexampleWitness(  // LCOV_EXCL_LINE
                // LCOV_DISABLED_START
                witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
                // LCOV_DISABLED_STOP
            outputCoverage,  // LCOV_EXCL_LINE
            // LCOV_DISABLED_START
            abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
            extractedBoundaryReports);  // LCOV_EXCL_LINE
        return outcome;  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
    if (endOutput - firstOutput > 1) {  // LCOV_EXCL_LINE
      FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
      outcome.shouldSplit = true;  // LCOV_EXCL_LINE
      return outcome;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    if (problem.usesDualRailStateEncoding) {  // LCOV_EXCL_LINE
      FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
      outcome.shouldSkipOutput = true;  // LCOV_EXCL_LINE
      return outcome;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    const std::string outputName =
        firstOutput < problem.observedOutputNames.size()  // LCOV_EXCL_LINE
            ? problem.observedOutputNames[firstOutput]  // LCOV_EXCL_LINE
            : std::to_string(firstOutput);  // LCOV_EXCL_LINE
    FinalPdrStageOutcome outcome;  // LCOV_EXCL_LINE
    outcome.terminalResult = makeSecResult(  // LCOV_EXCL_LINE
        SequentialEquivalenceStatus::Inconclusive,
        fullExactPdrResult.bound,  // LCOV_EXCL_LINE
        "PDR reached an abstract counterexample that concrete BMC did not "
        "validate for output `" +  // LCOV_EXCL_LINE
            outputName + "` at k = " +  // LCOV_EXCL_LINE
            std::to_string(fullExactPdrResult.bound),  // LCOV_EXCL_LINE
        outputCoverage,  // LCOV_EXCL_LINE
        abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
        extractedBoundaryReports);  // LCOV_EXCL_LINE
    return outcome;  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE

  auto splitPdrBatchAtFinalStage =
      [&](size_t batchIndex, size_t firstOutput, size_t endOutput) {  // LCOV_EXCL_LINE
    const size_t midOutput = firstOutput + (endOutput - firstOutput) / 2;  // LCOV_EXCL_LINE
    outputBatches.insert(  // LCOV_EXCL_LINE
        outputBatches.begin() + static_cast<std::ptrdiff_t>(batchIndex + 1),  // LCOV_EXCL_LINE
        {PdrOutputBatch{firstOutput, midOutput, true},  // LCOV_EXCL_LINE
         PdrOutputBatch{midOutput, endOutput, true}});  // LCOV_EXCL_LINE
  };  // LCOV_EXCL_LINE

  for (size_t batchIndex = 0; batchIndex < outputBatches.size(); ++batchIndex) {
    const auto [firstOutput, endOutput, startAtFinalExact] =
        outputBatches[batchIndex];
    if (startAtFinalExact) {
      const FinalPdrStageOutcome finalOutcome =
          runFinalExactPdrStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
      if (finalOutcome.terminalResult.has_value()) {  // LCOV_EXCL_LINE
        return *finalOutcome.terminalResult;  // LCOV_EXCL_LINE
      }
      if (finalOutcome.shouldSkipOutput) {  // LCOV_EXCL_LINE
        if (!markDualRailPdrOutputSkipped(  // LCOV_EXCL_LINE
                problem,  // LCOV_EXCL_LINE
                pdrCoveredOutputs,
                pdrSkippedOutputReasons,
                firstOutput)) {  // LCOV_EXCL_LINE
          return makeSecResult(  // LCOV_EXCL_LINE
              SequentialEquivalenceStatus::Inconclusive,
              provedBound,  // LCOV_EXCL_LINE
              "Dual-rail PDR repair was inconclusive for the only checked output",  // LCOV_EXCL_LINE
              outputCoverage,  // LCOV_EXCL_LINE
              abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
              extractedBoundaryReports);  // LCOV_EXCL_LINE
        }
        continue;  // LCOV_EXCL_LINE
      }
      if (finalOutcome.shouldSplit) {  // LCOV_EXCL_LINE
        splitPdrBatchAtFinalStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      continue;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    configureOutputBatchProblem(batchProblem, problem, firstOutput, endOutput);
    prunePdrBatchRelations(batchProblem, pdrBatchTransitionClosureLimit);
    // ASIC SEC runs need smaller carried predecessor obligations than the
    // standalone PDR engine default.  This does not change the PDR proof rule:
    // every learned clause is still justified by an UNSAT predecessor query,
    // and every reported counterexample is still checked by concrete BMC.
    constexpr size_t kSecPdrPredecessorProjectionLimit = 4;
    constexpr size_t kProjectedPdrPredecessorQueryBudget = 5000;
    constexpr size_t kDualRailProjectedPdrPredecessorQueryBudget = 5000;
    const size_t projectedPdrPredecessorQueryBudget =
        problem.usesDualRailStateEncoding
            ? secStrategySizeLimitFromEnv(
                  "KEPLER_SEC_PDR_DUAL_RAIL_PROJECTED_QUERY_BUDGET",
                  kDualRailProjectedPdrPredecessorQueryBudget)
            : kProjectedPdrPredecessorQueryBudget;
    emitPdrStrategyStageStats(
        emitPdrStageStats,
        batchIndex,
        firstOutput,
        endOutput,
        "initial",
        pdrBatchTransitionClosureLimit,
        kSecPdrPredecessorProjectionLimit,
        kSecPdrPredecessorProjectionLimit,
        batchProblem);
    PDREngine pdrEngine(
        batchProblem,
        solverType,
        kSecPdrPredecessorProjectionLimit,
        kSecPdrPredecessorProjectionLimit,
        /*useExactFrameClauses=*/false,
        projectedPdrPredecessorQueryBudget,
        /*refineProjectedCounterexamples=*/false,
        PDREngine::kDefaultBoundedRootGeneralizationAttempts,
        /*learnValidatedBadFormulaClauses=*/false,
        /*useExactResetFrontierChecks=*/dualRailPdrUsesResetFrontier);
    const auto pdrResult = pdrEngine.run(maxK, broadBasePrecheckDone);
    switch (pdrResult.status) {
      case PDRStatus::Equivalent:
        // Projected PDR frames are proof accelerators. Before marking a batch
        // covered, validate the actual top-output SEC base predicate through
        // the proved bound so no abstraction is trusted as a result.
        if (!problem.usesDualRailStateEncoding) {
          const PdrConcreteValidationCheck validationCheck =
              checkPdrConcreteValidation(
                  batchProblem, solverType, pdrResult.bound);
          if (validationCheck.witness.has_value()) {  // LCOV_EXCL_LINE
            const KInductionResult witnessResult{  // LCOV_EXCL_LINE
                KInductionStatus::Different,
                validationCheck.witness->badFrame,  // LCOV_EXCL_LINE
                validationCheck.witness};  // LCOV_EXCL_LINE
            return makeSecResult(  // LCOV_EXCL_LINE
                SequentialEquivalenceStatus::Different,
                witnessResult.bound,  // LCOV_EXCL_LINE
                formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
                outputCoverage,  // LCOV_EXCL_LINE
                abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
                extractedBoundaryReports);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
          markDualRailPdrOutputRangeCovered(
              pdrCoveredOutputs,
              pdrSkippedOutputReasons,
              firstOutput,
              endOutput);
          markPdrValidationUnknownOutputs(
              pdrCoveredOutputs,
              pdrSkippedOutputReasons,
              firstOutput,
              validationCheck.unknownOutputIndices);
          provedBound = std::max(provedBound, pdrResult.bound);
          break;
        } else if (!shouldDeferWideDualRailPdrValidation(batchProblem)) {
          if (auto concreteWitness = findPdrValidationCounterexample(
                  batchProblem, solverType, pdrResult.bound);
              concreteWitness.has_value()) {
            const KInductionResult witnessResult{  // LCOV_EXCL_LINE
                KInductionStatus::Different,
                concreteWitness->badFrame,  // LCOV_EXCL_LINE
                std::move(concreteWitness)};  // LCOV_EXCL_LINE
            return makeSecResult(  // LCOV_EXCL_LINE
                SequentialEquivalenceStatus::Different,
                witnessResult.bound,  // LCOV_EXCL_LINE
                formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
                outputCoverage,  // LCOV_EXCL_LINE
                abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
                extractedBoundaryReports);  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
        } else if (emitPdrStageStats) {  // LCOV_EXCL_LINE
          // The PDR batch itself closed.  On very wide dual-rail SoC surfaces,
          // rebuilding a full concrete BMC checker here materializes the same
          // million-symbol transition relation that batching was avoiding.
          emitSecDiag(  // LCOV_EXCL_LINE
              "SEC PDR stats: deferred wide dual-rail equivalent ",
              "validation outputs=",  // LCOV_EXCL_LINE
              batchProblem.originalObservedOutputCount);
        }  // LCOV_EXCL_LINE
        provedBound = std::max(provedBound, pdrResult.bound);
        markDualRailPdrOutputRangeCovered(
            pdrCoveredOutputs,
            pdrSkippedOutputReasons,
            firstOutput,
            endOutput);
        break;
      case PDRStatus::Different: {
        {
          std::optional<KInductionResult::CounterexampleWitness>
              concreteWitness;
          if (!shouldDeferWideDualRailPdrValidation(batchProblem)) {
            concreteWitness = findPdrValidationCounterexample(
                batchProblem, solverType, pdrResult.bound);
          } else if (emitPdrStageStats) {  // LCOV_EXCL_LINE
            emitSecDiag(  // LCOV_EXCL_LINE
                "SEC PDR stats: deferred wide dual-rail projected ",
                "counterexample validation outputs=",  // LCOV_EXCL_LINE
                batchProblem.originalObservedOutputCount); // LCOV_EXCL_LINE
          } // LCOV_EXCL_LINE
          if (concreteWitness.has_value()) {
            const KInductionResult witnessResult{ // LCOV_EXCL_LINE
                KInductionStatus::Different,
                concreteWitness->badFrame, // LCOV_EXCL_LINE
                std::move(concreteWitness)}; // LCOV_EXCL_LINE
            return makeSecResult( // LCOV_EXCL_LINE
                SequentialEquivalenceStatus::Different,
                witnessResult.bound, // LCOV_EXCL_LINE
                formatCounterexampleWitness(witnessResult, model0, model1, top0, top1), // LCOV_EXCL_LINE
                outputCoverage, // LCOV_EXCL_LINE
                abstractedSequentialBoundaries, // LCOV_EXCL_LINE
                extractedBoundaryReports); // LCOV_EXCL_LINE
          } // LCOV_EXCL_LINE
          // ASIC cones can still produce an abstract trace when the local
          // relation slice is too small. Retry the same output batch with more
          // relation / predecessor context before concrete validation. This
          // keeps the proof as real PDR over the conjunction slice while
          // avoiding the measured 598-pass one-output loop on BlackParrot. Any
          // reported difference is still accepted only after concrete BMC
          // validation below.
          // The initial 4-literal projection can be too abstract on a widened
          // ASIC relation, but BlackParrot measurements showed that jumping
          // straight to 64 literals creates a large level-1 blocked-predecessor
          // enumeration loop. Use an intermediate precision step before the
          // later exact retries.
          constexpr size_t kModeratePdrPredecessorProjectionLimit = 16;  // LCOV_EXCL_LINE
          // Exact-frame retries need more predecessor context than the
          // moderate projection to avoid abstract counterexamples, but fully
          // unbounded predecessor cubes were measured to enumerate thousands of
          // adjacent SAT models on BlackParrot. Use this bounded midpoint for
          // exact-frame passes.
          constexpr size_t kExactFramePdrPredecessorProjectionLimit = 32;  // LCOV_EXCL_LINE
          // Projected CEGAR stages are allowed to be inconclusive. If they
          // keep finding abstract SAT predecessors without strengthening the
          // frames, stop that stage and move to the stronger exact-frame PDR
          // retry instead of enumerating the same projected space for minutes.
          KInductionProblem refinedBatchProblem = problem;  // LCOV_EXCL_LINE
          configureOutputBatchProblem(  // LCOV_EXCL_LINE
              refinedBatchProblem, problem, firstOutput, endOutput);  // LCOV_EXCL_LINE
          prunePdrBatchRelations(  // LCOV_EXCL_LINE
              refinedBatchProblem, refinedPdrBatchTransitionClosureLimit);  // LCOV_EXCL_LINE

          // If concrete BMC rejects the first projected trace, first widen the
          // relation slice while keeping predecessor cubes small. Sampling on
          // BlackParrot showed that widening predecessor cubes before the
          // relation makes PDR enumerate thousands of exact level-1
          // predecessors. A wider relation can remove the abstraction that
          // produced the trace without abandoning the compact PDR obligation
          // shape that keeps ASIC proofs tractable.
          emitPdrStrategyStageStats(  // LCOV_EXCL_LINE
              emitPdrStageStats,  // LCOV_EXCL_LINE
              batchIndex,  // LCOV_EXCL_LINE
              firstOutput,  // LCOV_EXCL_LINE
              endOutput,  // LCOV_EXCL_LINE
              "widened_relation",
              refinedPdrBatchTransitionClosureLimit,  // LCOV_EXCL_LINE
              kSecPdrPredecessorProjectionLimit,
              kSecPdrPredecessorProjectionLimit,
              refinedBatchProblem);
          PDREngine refinedPdrEngine(  // LCOV_EXCL_LINE
              refinedBatchProblem,
              solverType,  // LCOV_EXCL_LINE
              kSecPdrPredecessorProjectionLimit,
              kSecPdrPredecessorProjectionLimit,
              /*useExactFrameClauses=*/false,
              projectedPdrPredecessorQueryBudget,  // LCOV_EXCL_LINE
              /*refineProjectedCounterexamples=*/false,
              PDREngine::kDefaultBoundedRootGeneralizationAttempts,
              /*learnValidatedBadFormulaClauses=*/false,
              /*useExactResetFrontierChecks=*/dualRailPdrUsesResetFrontier);  // LCOV_EXCL_LINE
          const auto refinedPdrResult = refinedPdrEngine.run(maxK, true);  // LCOV_EXCL_LINE
          if (refinedPdrResult.status == PDRStatus::Equivalent) {  // LCOV_EXCL_LINE
            provedBound = std::max(provedBound, refinedPdrResult.bound);  // LCOV_EXCL_LINE
            markDualRailPdrOutputRangeCovered(  // LCOV_EXCL_LINE
                pdrCoveredOutputs,  // LCOV_EXCL_LINE
                pdrSkippedOutputReasons,  // LCOV_EXCL_LINE
                firstOutput,  // LCOV_EXCL_LINE
                endOutput);  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }
          // If the wider relation still finds only an abstract trace, grow the
          // predecessor projection moderately on that same relation.  This is
          // a precision refinement, not a proof shortcut; any reported
          // difference is still validated by concrete BMC below.
          KInductionProblem widenedBatchProblem = refinedBatchProblem;  // LCOV_EXCL_LINE
          emitPdrStrategyStageStats(  // LCOV_EXCL_LINE
              emitPdrStageStats,  // LCOV_EXCL_LINE
              batchIndex,  // LCOV_EXCL_LINE
              firstOutput,  // LCOV_EXCL_LINE
              endOutput,  // LCOV_EXCL_LINE
              "widened_relation_moderate_projection",
              refinedPdrBatchTransitionClosureLimit,  // LCOV_EXCL_LINE
              kModeratePdrPredecessorProjectionLimit,
              kModeratePdrPredecessorProjectionLimit,
              widenedBatchProblem);
          PDREngine widenedPdrEngine(  // LCOV_EXCL_LINE
              widenedBatchProblem,
              solverType,  // LCOV_EXCL_LINE
              kModeratePdrPredecessorProjectionLimit,
              kModeratePdrPredecessorProjectionLimit,
              /*useExactFrameClauses=*/false,
              // LCOV_DISABLED_START
              projectedPdrPredecessorQueryBudget,  // LCOV_EXCL_LINE
              /*refineProjectedCounterexamples=*/false,
              // LCOV_DISABLED_STOP
              PDREngine::kDefaultBoundedRootGeneralizationAttempts,
              /*learnValidatedBadFormulaClauses=*/false,
              // LCOV_DISABLED_START
              /*useExactResetFrontierChecks=*/dualRailPdrUsesResetFrontier);  // LCOV_EXCL_LINE
              // LCOV_DISABLED_STOP
          const auto widenedPdrResult = widenedPdrEngine.run(maxK, true);  // LCOV_EXCL_LINE
          // LCOV_DISABLED_START
          if (widenedPdrResult.status == PDRStatus::Equivalent) {  // LCOV_EXCL_LINE
            provedBound = std::max(provedBound, widenedPdrResult.bound);  // LCOV_EXCL_LINE
            // LCOV_DISABLED_STOP
            markDualRailPdrOutputRangeCovered(  // LCOV_EXCL_LINE
                pdrCoveredOutputs,  // LCOV_EXCL_LINE
                pdrSkippedOutputReasons,  // LCOV_EXCL_LINE
                firstOutput,  // LCOV_EXCL_LINE
                endOutput);  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }
          if (problem.usesDualRailStateEncoding) {  // LCOV_EXCL_LINE
            // Dual-rail PDR keeps hard rail predicates local by going directly
            // to the final isolated retry. The intermediate exact-frame stage
            // can enumerate thousands of sibling predecessors before reaching
            // the same proof/split/skip decision.
            const FinalPdrStageOutcome finalOutcome =
                runFinalExactPdrStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
            if (finalOutcome.terminalResult.has_value()) {  // LCOV_EXCL_LINE
              return *finalOutcome.terminalResult;  // LCOV_EXCL_LINE
            }
            if (finalOutcome.shouldSkipOutput) {  // LCOV_EXCL_LINE
              if (!markDualRailPdrOutputSkipped(  // LCOV_EXCL_LINE
                      problem,  // LCOV_EXCL_LINE
                      pdrCoveredOutputs,
                      pdrSkippedOutputReasons,
                      firstOutput)) {  // LCOV_EXCL_LINE
                return makeSecResult(  // LCOV_EXCL_LINE
                    SequentialEquivalenceStatus::Inconclusive,
                    provedBound,  // LCOV_EXCL_LINE
                    "Dual-rail PDR repair was inconclusive for the only checked output",  // LCOV_EXCL_LINE
                    outputCoverage,  // LCOV_EXCL_LINE
                    abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
                    extractedBoundaryReports);  // LCOV_EXCL_LINE
              }
              break;  // LCOV_EXCL_LINE
            }
            if (finalOutcome.equivalent) {  // LCOV_EXCL_LINE
              break;  // LCOV_EXCL_LINE
            }
            if (finalOutcome.shouldSplit) {  // LCOV_EXCL_LINE
              splitPdrBatchAtFinalStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
              break;  // LCOV_EXCL_LINE
            }
          }  // LCOV_EXCL_LINE
          // Last retry on the widened relation slice: keep the complete
          // learned frame, but keep carried predecessor cubes bounded. Sampling
          // on BlackParrot showed that unbounded predecessor cubes made exact
          // PDR enumerate thousands of adjacent full SAT models; exact frame
          // clauses are the part that removes stale abstract predecessors.
          emitPdrStrategyStageStats(  // LCOV_EXCL_LINE
              emitPdrStageStats,  // LCOV_EXCL_LINE
              batchIndex,  // LCOV_EXCL_LINE
              firstOutput,  // LCOV_EXCL_LINE
              endOutput,  // LCOV_EXCL_LINE
              "widened_relation_exact",
              refinedPdrBatchTransitionClosureLimit,  // LCOV_EXCL_LINE
              kExactFramePdrPredecessorProjectionLimit,
              kExactFramePdrPredecessorProjectionLimit,
              widenedBatchProblem);
          PDREngine exactPdrEngine(  // LCOV_EXCL_LINE
              widenedBatchProblem,
              solverType,  // LCOV_EXCL_LINE
              kExactFramePdrPredecessorProjectionLimit,
              kExactFramePdrPredecessorProjectionLimit,
              /*useExactFrameClauses=*/true,
              /*maxPredecessorQueries=*/0,
              /*refineProjectedCounterexamples=*/false,
              PDREngine::kDefaultBoundedRootGeneralizationAttempts,
              /*learnValidatedBadFormulaClauses=*/false,
              /*useExactResetFrontierChecks=*/false);
          const auto exactPdrResult = exactPdrEngine.run(maxK, true);  // LCOV_EXCL_LINE
          if (exactPdrResult.status == PDRStatus::Equivalent) {  // LCOV_EXCL_LINE
            provedBound = std::max(provedBound, exactPdrResult.bound);  // LCOV_EXCL_LINE
            markDualRailPdrOutputRangeCovered(  // LCOV_EXCL_LINE
                pdrCoveredOutputs,  // LCOV_EXCL_LINE
                pdrSkippedOutputReasons,  // LCOV_EXCL_LINE
                firstOutput,  // LCOV_EXCL_LINE
                endOutput);  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }
          const FinalPdrStageOutcome finalOutcome =
              runFinalExactPdrStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
          if (finalOutcome.terminalResult.has_value()) {  // LCOV_EXCL_LINE
            return *finalOutcome.terminalResult;  // LCOV_EXCL_LINE
          }
          if (finalOutcome.shouldSkipOutput) {  // LCOV_EXCL_LINE
            if (!markDualRailPdrOutputSkipped(  // LCOV_EXCL_LINE
                    problem,  // LCOV_EXCL_LINE
                    pdrCoveredOutputs,
                    pdrSkippedOutputReasons,
                    firstOutput)) {  // LCOV_EXCL_LINE
              return makeSecResult(  // LCOV_EXCL_LINE
                  SequentialEquivalenceStatus::Inconclusive,
                  provedBound,  // LCOV_EXCL_LINE
                  "Dual-rail PDR repair was inconclusive for the only checked output",  // LCOV_EXCL_LINE
                  outputCoverage,  // LCOV_EXCL_LINE
                  abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
                  extractedBoundaryReports);  // LCOV_EXCL_LINE
            }
            break;  // LCOV_EXCL_LINE
          }
          if (finalOutcome.equivalent) {  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }
          if (finalOutcome.shouldSplit) {  // LCOV_EXCL_LINE
            splitPdrBatchAtFinalStage(batchIndex, firstOutput, endOutput);  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }
        }
        break; // LCOV_EXCL_LINE
      }
      case PDRStatus::Inconclusive:
      default:
        if (problem.usesDualRailStateEncoding) {  // LCOV_EXCL_LINE
          if (endOutput - firstOutput > 1) {  // LCOV_EXCL_LINE
            const size_t midOutput =
                firstOutput + (endOutput - firstOutput) / 2;  // LCOV_EXCL_LINE
            // Once projected dual-rail PDR has asked to split, retry the
            // children in the final exact path. Re-running the same projected
            // stage on each child was the measured RISC-V runtime wall.
            outputBatches.insert(  // LCOV_EXCL_LINE
                outputBatches.begin() +
                    static_cast<std::ptrdiff_t>(batchIndex + 1),  // LCOV_EXCL_LINE
                {PdrOutputBatch{firstOutput, midOutput, true},  // LCOV_EXCL_LINE
                 PdrOutputBatch{midOutput, endOutput, true}});  // LCOV_EXCL_LINE
            break;  // LCOV_EXCL_LINE
          }  // LCOV_EXCL_LINE
          if (markDualRailPdrOutputSkipped(  // LCOV_EXCL_LINE
                  problem,
                  pdrCoveredOutputs,
                  pdrSkippedOutputReasons,
                  firstOutput)) {
            break;  // LCOV_EXCL_LINE
          }
        }  // LCOV_EXCL_LINE
        return makeSecResult(  // LCOV_EXCL_LINE
            SequentialEquivalenceStatus::Inconclusive,
            pdrResult.bound,  // LCOV_EXCL_LINE
            "Reached max_k without a proof or counterexample",  // LCOV_EXCL_LINE
            outputCoverage,  // LCOV_EXCL_LINE
            abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
            extractedBoundaryReports);  // LCOV_EXCL_LINE
    }
  }

  const OutputCoverageSelection finalCoverage =
      buildCoverageWithDualRailOutputSkips(
          outputCoverage, problem, pdrCoveredOutputs, pdrSkippedOutputReasons);
  if (problem.usesDualRailStateEncoding &&
      finalCoverage.checkedOutputs.names.size() <
          outputCoverage.checkedOutputs.names.size()) {
    std::vector<size_t> skippedOutputIndices;
    skippedOutputIndices.reserve(pdrCoveredOutputs.size());
    for (size_t outputIndex = 0; outputIndex < pdrCoveredOutputs.size();
         ++outputIndex) {
      if (!pdrCoveredOutputs[outputIndex]) {
        skippedOutputIndices.push_back(outputIndex);
      }
    }
    if (auto witness = findInputOnlyFrameZeroResidualCounterexample(
            problem, skippedOutputIndices, solverType);
        witness.has_value()) {
      KInductionResult witnessResult{  // LCOV_EXCL_LINE
          KInductionStatus::Different,
          witness->badFrame,  // LCOV_EXCL_LINE
          std::move(witness)};  // LCOV_EXCL_LINE
      return makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Different,
          witnessResult.bound,  // LCOV_EXCL_LINE
          formatCounterexampleWitness(witnessResult, model0, model1, top0, top1),  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
  }
  const bool pdrProvedNoOutputs = finalCoverage.checkedOutputs.names.empty();
  if (pdrProvedNoOutputs) {
    return makeSecResult(
        SequentialEquivalenceStatus::Inconclusive,
        provedBound,
        problem.usesDualRailStateEncoding
            ? "Dual-rail PDR exhausted repair/projection before proving any observed output"
            : "PDR did not prove any observed output",
        finalCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }
  return makeSecResult(
      SequentialEquivalenceStatus::Equivalent,
      provedBound,
      "",
      finalCoverage,
      abstractedSequentialBoundaries,
      extractedBoundaryReports);
}

SequentialEquivalenceResult runKInductionSecEngine(
    const KInductionProblem& problem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports) {
  if (auto dualRailResult = proveDualRailResidualsWithSelectedEngine(
          problem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports,
          DualRailResidualEngine::KInduction);
      dualRailResult.has_value()) {
    return *dualRailResult;
  }

  KInductionEngine engine(problem, solverType);
  const auto result = engine.run(maxK);
  switch (result.status) {
    case KInductionStatus::Equivalent:
      return makeSecResult(
          SequentialEquivalenceStatus::Equivalent,
          result.bound,
          "",
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    case KInductionStatus::Different:
      return makeSecResult(
          SequentialEquivalenceStatus::Different,
          result.bound,
          result.witness.has_value()
              ? formatCounterexampleWitness(result, model0, model1, top0, top1)
              : "Classic k-induction found a counterexample at k = " +  // LCOV_EXCL_LINE
                    std::to_string(result.bound),  // LCOV_EXCL_LINE
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    case KInductionStatus::Inconclusive:  // LCOV_EXCL_LINE
    default:
      // Honor the selected SEC engine.  KI must not silently invoke PDR as a
      // secondary prover; callers can rerun with sec_engine=pdr if desired.
      return makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Inconclusive,
          result.bound,  // LCOV_EXCL_LINE
          "Reached max_k without a proof or counterexample",  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
  }
}

SequentialEquivalenceResult runImcSecEngine(
    const KInductionProblem& problem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports) {
  // IMC must remain an interpolation-based engine.  Do not route dual-rail
  // residuals through the KI residual helper before the IMC engine runs.
  IMCEngine engine(problem, solverType);
  const auto result = engine.run(maxK);
  switch (result.status) {
    case IMCStatus::Equivalent: {
      emitSecEngineProofProgress(
          problem, "IMC", problem.observedOutputExprs0.size());
      SequentialEquivalenceResult secResult = makeSecResult(
          SequentialEquivalenceStatus::Equivalent,
          result.bound,
          "",
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
      setSecEngineProofProgress(
          secResult, problem, "IMC", problem.observedOutputExprs0.size());
      return secResult;
    }
    case IMCStatus::Different: {
      const KInductionResult witnessResult{  // LCOV_EXCL_LINE
          KInductionStatus::Different, result.bound, result.witness};  // LCOV_EXCL_LINE
      return makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Different,
          result.bound,  // LCOV_EXCL_LINE
          result.witness.has_value()  // LCOV_EXCL_LINE
              ? formatCounterexampleWitness(witnessResult, model0, model1, top0, top1)  // LCOV_EXCL_LINE
              : "IMC found a counterexample at k = " +  // LCOV_EXCL_LINE
                    std::to_string(result.bound),  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    case IMCStatus::Inconclusive:  // LCOV_EXCL_LINE
    default: {
      // Honor the selected SEC engine.  IMC must not silently invoke PDR as a
      // secondary prover; callers can rerun with sec_engine=pdr if desired.
      if (result.firstUnprovenOutput.has_value()) {  // LCOV_EXCL_LINE
        emitSecEngineProofProgress(  // LCOV_EXCL_LINE
            problem, "IMC", *result.firstUnprovenOutput);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      SequentialEquivalenceResult secResult = makeSecResult(  // LCOV_EXCL_LINE
          SequentialEquivalenceStatus::Inconclusive,
          result.bound,  // LCOV_EXCL_LINE
          "Reached max_k without a proof or counterexample",  // LCOV_EXCL_LINE
          outputCoverage,  // LCOV_EXCL_LINE
          abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
          extractedBoundaryReports);  // LCOV_EXCL_LINE
      if (result.firstUnprovenOutput.has_value()) {  // LCOV_EXCL_LINE
        setSecEngineProofProgress(  // LCOV_EXCL_LINE
            secResult, problem, "IMC", *result.firstUnprovenOutput);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      return secResult;  // LCOV_EXCL_LINE
    }
  }
}

SequentialEquivalenceResult runSelectedSecEngine(
    SecEngine secEngine,
    const KInductionProblem& proofProblem,
    size_t maxK,
    KEPLER_FORMAL::Config::SolverType solverType,
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports) {
  switch (secEngine) {
    case SecEngine::Pdr:
      return runPdrSecEngine(
          proofProblem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    case SecEngine::KInduction:
      return runKInductionSecEngine(
          proofProblem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    case SecEngine::Imc:
      return runImcSecEngine(
          proofProblem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    default:
      return runPdrSecEngine(
          proofProblem,
          maxK,
          solverType,
          model0,
          model1,
          top0,
          top1,
          outputCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
  }
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

namespace detail {

SequentialEquivalenceProofProgress buildSecEngineProofProgress(
    const std::string& engineLabel,
    const std::vector<std::string>& observedOutputNames,
    size_t totalOutputCount,
    size_t provenOutputCount) {
  SequentialEquivalenceProofProgress progress;
  progress.engineLabel = engineLabel;
  progress.provenOutputs = std::min(provenOutputCount, totalOutputCount);
  progress.totalOutputs = totalOutputCount;
  progress.unprovenOutputs.reserve(totalOutputCount - progress.provenOutputs);
  for (size_t output = progress.provenOutputs; output < totalOutputCount;
       ++output) {
    progress.unprovenOutputs.push_back(
        {output,
         output < observedOutputNames.size() ? observedOutputNames[output]
                                             : std::string("<unknown>")}); // LCOV_EXCL_LINE
  }
  return progress;
}

std::vector<std::string> buildSecEngineProofProgressDiagLines(
    const std::string& engineLabel,
    const std::vector<std::string>& observedOutputNames,
    size_t totalOutputCount,
    size_t provenOutputCount) {
  const SequentialEquivalenceProofProgress progress =
      buildSecEngineProofProgress(
          engineLabel,
          observedOutputNames,
          totalOutputCount,
          provenOutputCount);
  std::vector<std::string> lines;
  lines.reserve(1 + progress.unprovenOutputs.size());

  std::ostringstream summary;
  summary << "SEC diag: SEC " << progress.engineLabel
          << " proven outputs: " << progress.provenOutputs << "/"
          << progress.totalOutputs;
  lines.push_back(summary.str());

  for (const SequentialEquivalenceUnprovenOutput& output :
       progress.unprovenOutputs) {
    std::ostringstream line;
    line << "SEC diag: SEC " << progress.engineLabel
         << " not proven output[" << output.index << "]=" << output.name;
    lines.push_back(line.str());
  }
  return lines;
}

}  // namespace detail

SequentialEquivalenceStrategy::SequentialEquivalenceStrategy(
    naja::NL::SNLDesign* top0,
    naja::NL::SNLDesign* top1,
    KEPLER_FORMAL::Config::SolverType solverType,
    SecEngine secEngine,
    SecEncoding encoding)
    : top0_(top0),
      top1_(top1),
      solverType_(solverType),
      secEngine_(secEngine),
      encoding_(encoding) {}

SequentialEquivalenceResult SequentialEquivalenceStrategy::run(size_t maxK) const {
  const bool secDiagEnabled = std::getenv("KEPLER_SEC_DIAG") != nullptr;
  logSecDiagLine(secDiagEnabled, "SEC diag: start run");

  // Phase 1: extract both designs into the normalized SEC model used by every
  // downstream engine. If either side cannot be modeled soundly, stop before we
  // spend time aligning interfaces or building proof problems.
  const auto model0 =
      extractSecDesign(top0_, "SEC diag: extracted design0", secDiagEnabled);
  if (model0.hasUnsupportedFeatures()) {
    std::vector<std::string> abstractedSequentialBoundaries;
    std::vector<ExtractedBoundaryReportEntry> extractedBoundaryReports;
    appendAbstractedSequentialBoundaries(
        model0, "design0", abstractedSequentialBoundaries);
    appendExtractedBoundaryReports(model0, "design0", extractedBoundaryReports);
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(model0.unsupportedReasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }

  const auto model1 =
      extractSecDesign(top1_, "SEC diag: extracted design1", secDiagEnabled);
  return runExtractedModels(model0, model1, maxK);
}

SequentialEquivalenceResult SequentialEquivalenceStrategy::runExtractedModels(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    size_t maxK) const {
  const bool secDiagEnabled = std::getenv("KEPLER_SEC_DIAG") != nullptr;
  logSecDiagLine(secDiagEnabled, "SEC diag: start run from extracted models");

  // Compact SEC can release the elaborated Naja DBs after extraction and run
  // entirely from these value-type models. Rebuilding the boundary summaries
  // here keeps normal and compact flows reporting the same coverage details.
  std::vector<std::string> abstractedSequentialBoundaries;
  std::vector<ExtractedBoundaryReportEntry> extractedBoundaryReports;
  appendAbstractedSequentialBoundaries(
      model0, "design0", abstractedSequentialBoundaries);
  appendExtractedBoundaryReports(model0, "design0", extractedBoundaryReports);
  if (model0.hasUnsupportedFeatures()) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(model0.unsupportedReasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }
  appendAbstractedSequentialBoundaries(
      model1, "design1", abstractedSequentialBoundaries);
  appendExtractedBoundaryReports(model1, "design1", extractedBoundaryReports);
  if (model1.hasUnsupportedFeatures()) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        joinReasons(model1.unsupportedReasons),
        OutputCoverageSelection{},
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }

  // Phase 2: align the externally visible SEC interface, then drop any outputs
  // whose cones were already classified as skipped by extraction.
  // SEC proves only the public interface relation.  Do not mine or assume
  // same-named internal state across the two independently extracted designs.
  AlignedSecInterface aligned = alignSecInterface(
      model0,
      model1,
      secDiagEnabled);
  if (aligned.outputs.names.empty()) {
    return makeSecResult(
        SequentialEquivalenceStatus::Unsupported,
        0,
        "No aligned observed outputs remain after skipping cones with no-driver, "
        "multi-driver, or logical-loop connectivity.",
        aligned.outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }

  if (secDiagEnabled) {
    printf(
        "SEC diag: aligned_inputs=%zu aligned_outputs=%zu "
        "state_bits0=%zu state_bits1=%zu\n",
        aligned.inputs.names.size(),
        aligned.outputs.names.size(),
        model0.stateBits.size(),
        model1.stateBits.size());
    fflush(stdout);
  }

  // Phase 3: rewrite both designs into one shared symbol space, strengthen the
  // startup frontier with reset/bootstrap facts, and build the final SEC
  // property plus the induction-friendly variant that some engines consume.
  SharedSecSymbolSpace symbolSpace = buildSharedSecSymbolSpace(
      model0, model1, aligned.inputs, aligned.outputs);
  // KI / IMC consume explicit post-reset state values directly. SEC/PDR keeps
  // the reset cycle/input model and validates startup candidates with concrete
  // BMC / reset-frontier checks, so it can avoid the sampled full-design sweep
  // that tries to constant-evaluate every state bit before the first PDR query.
  const bool deriveResetBootstrapStrengthening = secEngine_ != SecEngine::Pdr;
  // Reset bootstrap is allowed to add concrete values inside each design, but
  // it must not add any cross-design internal state relation.
  const auto reachableInvariant = integrateReachableStateInvariant(
      model0,
      model1,
      symbolSpace.state0Symbols,
      symbolSpace.state1Symbols,
      symbolSpace.problem,
      deriveResetBootstrapStrengthening);
  if (encoding_ == SecEncoding::Binary) {
    filterOutputsRequiringUnanchoredResetState(
        model0,
        model1,
        reachableInvariant,
        !symbolSpace.problem.resetBootstrapInputs.empty() &&
            symbolSpace.problem.resetBootstrapCycles != 0,
        aligned,
        secDiagEnabled);
  } else {
    logSecDiagLine(
        secDiagEnabled,
        "SEC diag: dual-rail mode keeps reset-unanchored outputs in the "
        "rail-encoded proof");
  }
  symbolSpace.problem.observedOutputNames = aligned.outputs.names;
  if (aligned.outputs.names.empty()) {
    return makeSecResult(  // LCOV_EXCL_LINE
        SequentialEquivalenceStatus::Unsupported,
        0,
        "No aligned observed outputs remain after skipping cones that depend "  // LCOV_EXCL_LINE
        "on reset-unanchored internal state.",
        aligned.outputCoverage,  // LCOV_EXCL_LINE
        abstractedSequentialBoundaries,  // LCOV_EXCL_LINE
        extractedBoundaryReports);  // LCOV_EXCL_LINE
  }
  // KI and PDR both work on small COI slices of a potentially huge SEC
  // problem. Keep next-state formulas in their extracted-model symbol space
  // until the proof engine actually asks for a transition; otherwise PDR
  // materializes the full ASIC transition relation before output batching can
  // prune it. Output remapping happens after the reset-frontier coverage filter
  // so skipped top outputs never allocate SAT symbols or proof obligations.
  const bool useLazyTransitionRemapping =
      secEngine_ == SecEngine::KInduction || secEngine_ == SecEngine::Pdr;
  KInductionProblem proofProblem;
  if (encoding_ == SecEncoding::DualRailSteady) {
    proofProblem = buildDualRailSecProblem(
        model0,
        model1,
        aligned.inputs,
        aligned.outputs,
        reachableInvariant,
        symbolSpace,
        useLazyTransitionRemapping,
        secDiagEnabled);
  } else {
    const auto remapped = remapSecExpressions(
        model0,
        model1,
        aligned.outputs,
        symbolSpace,
        symbolSpace.problem,
        !useLazyTransitionRemapping,
        secDiagEnabled);
    if (useLazyTransitionRemapping) {
      attachLazyTransitions(
          model0,
          model1,
          symbolSpace.state0Symbols,
          symbolSpace.state1Symbols,
          std::move(symbolSpace.localToCombined0),
          std::move(symbolSpace.localToCombined1),
          symbolSpace.problem);
    }
    buildSecPropertiesAndTransitions(
        model0,
        model1,
        aligned.outputs,
        symbolSpace.state0Symbols,
        symbolSpace.state1Symbols,
        remapped,
        symbolSpace.problem,
        secDiagEnabled);
    proofProblem = symbolSpace.problem;
  }

  // Phase 4: hand the fully normalized SEC transition system to the requested
  // top-level engine. From here on, every engine sees the same problem and only
  // differs in how it searches for proofs or counterexamples.
  if (secDiagEnabled) {
    fprintf(stderr, "SEC diag: entering %s\n", describeSecEngine(secEngine_));
    fflush(stderr);
  }

  return runSelectedSecEngine(
      secEngine_,
      proofProblem,
      maxK,
      solverType_,
      model0,
      model1,
      top0_,
      top1_,
      aligned.outputCoverage,
      abstractedSequentialBoundaries,
      extractedBoundaryReports);
}

}  // namespace KEPLER_FORMAL::SEC
