// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "strategy/SequentialEquivalenceStrategy.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iterator>
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

#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__GLIBC__)
#include <malloc.h>
#endif

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
#include "../../sat/SATSolverWrapper.h"

namespace KEPLER_FORMAL::SEC {

// Overall SEC strategy pipeline:
// 1. Extract both designs into the normalized sequential model used by SEC.
// 2. Align environment inputs and observed outputs by stable external names.
// 3. Keep cross-design internal state uncorrelated; only public facts constrain
//    the two designs before the selected SEC engine proves outputs.
// 4. Remap both designs into one shared SAT symbol space.
// 5. Build F[0] from the extracted initial predicate.
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

void releasePdrBatchAllocatorPages() {
  // Recursive output batching destroys property-local PDR solvers between
  // ranges. Return their free allocator pages before the next exact range
  // starts so only the deliberately shared F[0] cache remains resident.
#if defined(__APPLE__)
  malloc_zone_pressure_relief(nullptr, 0);
#elif defined(__GLIBC__)
  malloc_trim(0);
#endif
}

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
  if (result.status == SequentialEquivalenceStatus::Equivalent &&
      result.coveredOutputs > 0 &&
      result.coveredOutputs < result.totalOutputs) {
    result.status = SequentialEquivalenceStatus::PartiallyProved;
  }
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

void markPdrOutputRangeCovered(
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

void markPdrOutputRangeSkipped(
    std::vector<bool>& coveredOutputs,
    std::unordered_map<size_t, std::string>& skipReasons,
    size_t firstOutput,
    size_t endOutput,
    const std::string& reason) {
  const size_t cappedEnd = std::min(endOutput, coveredOutputs.size());
  for (size_t outputIndex = firstOutput; outputIndex < cappedEnd;
       ++outputIndex) {
    coveredOutputs[outputIndex] = false;
    skipReasons[outputIndex] = reason;
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
  subset.dualRailOutputStrictEqualityExprs.clear();
  subset.dualRailOutputSkipReasons.clear();

  // LCOV_DISABLED_START
  const bool copyObservedKeys =
  // LCOV_DISABLED_STOP
      source.observedOutputs.size() == source.observedOutputExprs0.size();
  const bool copySkipReasons =
      source.dualRailOutputSkipReasons.size() ==
      source.observedOutputExprs0.size();
  const bool copyStrictEqualityExprs =
      source.dualRailOutputStrictEqualityExprs.size() ==
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
    if (copyStrictEqualityExprs) {
      subset.dualRailOutputStrictEqualityExprs.push_back(
          source.dualRailOutputStrictEqualityExprs[outputIndex]);
    }
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

void filterOutputsRequiringUninitializedState(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    bool hasIncompleteInitialState,
    AlignedSecInterface& aligned,
    bool secDiagEnabled) {
  if (!hasIncompleteInitialState || aligned.outputs.names.empty()) {
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
      // Binary SEC cannot distinguish an initialization-only mismatch from a
      // concrete design mismatch without an initial value for this state.
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
        "SEC diag: reset-anchored checked outputs=%s\n",
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

enum class DualRailResidualEngine {
  KInduction,
  Imc,
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
    case DualRailResidualEngine::Imc:
      return "IMC";
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

void markDualRailResidualOutputsSkipped(
    const std::vector<size_t>& outputIndices,
    const KInductionProblem& problem,
    DualRailResidualEngine engine,
    DualRailResidualProofState& proofState,
    const std::string& reason) {
  for (const size_t outputIndex : outputIndices) {
    markDualRailResidualOutputSkipped(
        outputIndex, problem, engine, proofState, reason);
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

  // This is a witness-only guard for skipped residual top outputs. Restrict it
  // to frame-0 input/constant mismatches so stateful residuals do not pay to
  // materialize transition cones before the selected SEC engine.
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
        std::string("Dual-rail ") +
            dualRailResidualEngineName(engine) +
            " did not prove any output",
        finalCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }

  if (coveredCount != proofState.coveredOutputs.size()) {
    return makeSecResult(
        SequentialEquivalenceStatus::PartiallyProved,
        proofState.provedBound,
        std::string("Dual-rail ") + dualRailResidualEngineName(engine) +
            " proved " + std::to_string(coveredCount) + " of " +
            std::to_string(proofState.coveredOutputs.size()) +
            " observed outputs; remaining outputs are inconclusive",
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
        "SEC diag: deferred next-state formula remapping");
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

BoolExpr* buildDualRailBinaryDefinedExpr(const DualRailBoolExpr& value) {
  // In the paper's encoding, 01 and 10 are binary values while 11 is X.
  // Legal-state constraints separately exclude the empty value 00.
  return BoolExpr::simplify(
      BoolExpr::Xor(value.mayBeOne, value.mayBeZero));
}

struct DualRailOutputProperties {
  BoolExpr* guardedEquality = nullptr;
  BoolExpr* strictEquality = nullptr;
};

DualRailOutputProperties buildDualRailOutputProperties(
    const DualRailBoolExpr& value0,
    const DualRailBoolExpr& value1) {
  BoolExpr* bothValuesDefined = BoolExpr::simplify(BoolExpr::And(
      buildDualRailBinaryDefinedExpr(value0),
      buildDualRailBinaryDefinedExpr(value1)));
  BoolExpr* strictEquality = BoolExpr::simplify(BoolExpr::And(
      makeEqualityExpr(value0.mayBeOne, value1.mayBeOne),
      makeEqualityExpr(value0.mayBeZero, value1.mayBeZero)));
  // Steady-state dual-rail SEC ignores cycles where either value is X and
  // rejects only opposite binary values. Strict rail equality remains metadata
  // for shared exact query surfaces.
  BoolExpr* binaryMismatch = BoolExpr::And(
      bothValuesDefined,
      BoolExpr::Xor(value0.mayBeOne, value1.mayBeOne));
  return {
      BoolExpr::simplify(BoolExpr::Not(binaryMismatch)),
      strictEquality};
}

const char* pdrStatusName(PDRStatus status) {
  switch (status) {
    case PDRStatus::Equivalent:
      return "equivalent";
    case PDRStatus::Different:
      return "different";
    case PDRStatus::Inconclusive:
    default:
      return "inconclusive";
  }
}

constexpr PDRQueryLimits kDualRailPdrBatchProbeLimits{
    /*predecessorConflictLimit=*/10 * 1000,
    /*predecessorDecisionLimit=*/150 * 1000,
    /*blockingConflictLimit=*/10 * 1000,
    /*blockingDecisionLimit=*/150 * 1000,
    /*predecessorEncodingNodeLimit=*/5 * 1000 * 1000,
    /*predecessorNodeHintTargetLimit=*/512};

// Per-call limits alone still allow one hard output to issue thousands of
// expensive CaDiCaL queries. Four full predecessor allowances give a singleton
// room to block several hard obligations while bounding its cumulative work.
constexpr size_t kDefaultDualRailPdrSingletonConflictBudget =
    1000 * 1000;
constexpr size_t kDefaultDualRailPdrSingletonDecisionBudget =
    40 * 1000 * 1000;
// CaDiCaL ticks count clause-cache-line work, bounding propagation-heavy
// queries that can do little visible decision or conflict work.
constexpr size_t kDefaultDualRailPdrSingletonTickBudget =
    100 * 1000 * 1000;

PDRResult runPdrOutputBatch(const PDREngine& engine,
                            size_t maxFrames,
                            BoolExpr* property,
                            size_t outputCount,
                            bool boundSingletonCadicalWork) {
  if (outputCount != 1) {
    // A broad UNKNOWN only schedules exact child properties. Full per-query
    // limits are reserved for singleton leaves, so failed probes cannot
    // dominate them.
    return engine.run(maxFrames, property, kDualRailPdrBatchProbeLimits);
  }
  if (!boundSingletonCadicalWork) {
    return engine.run(maxFrames, property);
  }

  SATSolverWrapper::CadicalWorkBudget budget(
      secStrategySizeLimitFromEnv(
          "KEPLER_SEC_PDR_DUAL_RAIL_SINGLETON_CONFLICT_BUDGET",
          kDefaultDualRailPdrSingletonConflictBudget),
      secStrategySizeLimitFromEnv(
          "KEPLER_SEC_PDR_DUAL_RAIL_SINGLETON_DECISION_BUDGET",
          kDefaultDualRailPdrSingletonDecisionBudget),
      secStrategySizeLimitFromEnv(
          "KEPLER_SEC_PDR_DUAL_RAIL_SINGLETON_TICK_BUDGET",
          kDefaultDualRailPdrSingletonTickBudget));
  PDRResult result;
  {
    SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
    result = engine.run(maxFrames, property);
  }
  if (pdrStrategyStatsEnabled()) {
    emitSecDiag(
        "SEC PDR stats: singleton CaDiCaL cumulative budget ",
        "conflicts=", budget.conflictsUsed(), "/", budget.conflictLimit(),
        " decisions=", budget.decisionsUsed(), "/", budget.decisionLimit(),
        " ticks=", budget.ticksUsed(), "/", budget.tickLimit(),
        " exhausted=", budget.exhausted() ? 1 : 0);
  }
  return result;
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

void integrateInitialState(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state0Symbols,
    const std::unordered_map<SignalKey, size_t, SignalKeyHash>& state1Symbols,
    KInductionProblem& problem) {
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
}

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
        "bad_is_false=%d induction_bad_is_false=%d\n",
        problem.property == BoolExpr::createTrue(),
        problem.inductionProperty == BoolExpr::createTrue(),
        problem.bad == BoolExpr::createFalse(),
        problem.inductionBad == BoolExpr::createFalse());
    fflush(stdout);
  }
}

KInductionProblem buildDualRailSecProblem(
    const SequentialDesignModel& model0,
    const SequentialDesignModel& model1,
    // LCOV_DISABLED_START
    const AlignedSignals& alignedInputs,
    const AlignedSignals& alignedOutputs,
    SharedSecSymbolSpace& symbolSpace,
    // LCOV_DISABLED_STOP
    bool useLazyTransitionRemapping,
    bool secDiagEnabled) {
  KInductionProblem problem;
  problem.environmentInputs = alignedInputs.keys0;
  problem.environmentInputNames = symbolSpace.problem.environmentInputNames;
  problem.inputSymbols = symbolSpace.problem.inputSymbols;
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
  // The rail-valued initial predicate is represented as structured unit
  // facts in initialStateAssignments.  Keep initialCondition non-null so the
  // existing base-case encoders enter their structured-init path without
  // materializing a huge duplicate conjunction over every rail.
  problem.initialCondition = BoolExpr::createTrue();

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
  problem.dualRailOutputStrictEqualityExprs.clear();
  problem.dualRailOutputStrictEqualityExprs.reserve(
      alignedOutputs.names.size());
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

    const auto outputProperties =
        buildDualRailOutputProperties(out0, out1);
    // LCOV_DISABLED_START
    // Keep batching/reporting aligned to real top outputs.  The rail pair is a
    // single ternary output value, so its may-one and may-zero equalities are
    // LCOV_DISABLED_STOP
    // one SEC obligation rather than two independent output bits.
    problem.observedOutputNames.push_back(alignedOutputs.names[i]);
    // LCOV_DISABLED_START
    problem.observedOutputExprs0.push_back(outputProperties.guardedEquality);
    problem.observedOutputExprs1.push_back(BoolExpr::createTrue());
    problem.dualRailOutputStrictEqualityExprs.push_back(
        outputProperties.strictEquality);
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
    property = BoolExpr::And(property, outputProperties.guardedEquality);
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
        "rail_outputs=%zu "
        "dual_rail_state_relation_pairs=%zu\n",
        problem.totalStateCount,
        problem.observedOutputExprs0.size(),
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
      return "pdr engine";  // LCOV_EXCL_LINE
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

void setExactSecEngineProofProgress(
    SequentialEquivalenceResult& result,
    const KInductionProblem& problem,
    const std::string& engineLabel,
    const std::vector<bool>& coveredOutputs) {
  SequentialEquivalenceProofProgress progress;
  progress.engineLabel = engineLabel;
  progress.totalOutputs = coveredOutputs.size();
  progress.provenOutputs = static_cast<size_t>(
      std::count(coveredOutputs.begin(), coveredOutputs.end(), true));
  for (size_t outputIndex = 0; outputIndex < coveredOutputs.size();
       ++outputIndex) {
    if (!coveredOutputs[outputIndex]) {
      progress.unprovenOutputs.push_back(
          {outputIndex, outputNameForProblemIndex(problem, outputIndex)});
    }
  }
  result.proofProgress = std::move(progress);
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

  // Split wide SEC conjunctions into exact output obligations. Each PDR run
  // still uses the complete transition system and exact F[0].
  // LCOV_DISABLED_STOP
  //
  // Keep PDR batches bounded. A resource-limited dual-rail batch is split
  // recursively, and only singleton leaves receive the full query budget.
  constexpr size_t kMinOutputsForBatchedPdrProof = 129;
  constexpr OutputBatchingLimits kPdrOutputBatchingLimits{32, 1024};
  constexpr OutputBatchingLimits kDualRailPdrOutputBatchingLimits{128, 8192};
  const OutputBatchingLimits pdrOutputBatchingLimits =
      // LCOV_DISABLED_START
      problem.usesDualRailStateEncoding ? dualRailPdrOutputBatchingLimits(
      // LCOV_DISABLED_STOP
                                              kDualRailPdrOutputBatchingLimits)
                                        // LCOV_DISABLED_START
                                        : kPdrOutputBatchingLimits;
                                        // LCOV_DISABLED_STOP
  // LCOV_DISABLED_START
  struct PdrOutputBatch {
    size_t firstOutput = 0;
    size_t endOutput = 0;
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
    // queries. On medium designs, each tiny batch repeats the same frame and
    // transition setup, so prove one conjunction slice and reserve batching
    // for BlackParrot/AES-scale output counts.
    outputBatches.push_back({0, problem.observedOutputExprs0.size()});  // LCOV_EXCL_LINE
    // LCOV_DISABLED_STOP
  } else {  // LCOV_EXCL_LINE
    for (const auto& [firstOutput, endOutput] :
         buildSupportBoundedOutputBatches(problem, pdrOutputBatchingLimits)) {
      outputBatches.push_back({firstOutput, endOutput});
    }
  }
  std::vector<bool> pdrCoveredOutputs =
      makeInitialPdrCoveredOutputs(problem);
  std::unordered_map<size_t, std::string> pdrSkippedOutputReasons =
      presetDualRailSkipReasons;
  size_t provedBound = 0;
  bool stopAfterInconclusiveBatch = false;
  std::shared_ptr<PDRExactInitCache> exactInitCache;
  if (problem.usesDualRailStateEncoding) {
    exactInitCache =
        std::make_shared<PDRExactInitCache>(problem, solverType);
  }
  KInductionProblem exactBatchProblem = problem;
  for (size_t batchIndex = 0; batchIndex < outputBatches.size(); ++batchIndex) {
    const auto [firstOutput, endOutput] = outputBatches[batchIndex];
    configureOutputBatchProblem(
        exactBatchProblem, problem, firstOutput, endOutput);
    if (isSecDiagEnabled()) {
      // Output batches may split after an inconclusive result. Record the live
      // range so one fully diagnostic performance run identifies the slow PDR
      // slice without changing batching or proof order.
      emitSecDiag(
          "SEC diag: PDR steady-state check begin index=",
          batchIndex,
          " pending_batches=",
          outputBatches.size(),
          " output_range=",
          firstOutput,
          "..",
          endOutput);
    }
    PDRResult pdrResult;
    {
      PDREngine pdrEngine(
          exactBatchProblem, solverType, 0, exactInitCache);
      pdrResult = runPdrOutputBatch(
          pdrEngine, maxK, exactBatchProblem.property,
          endOutput - firstOutput,
          problem.usesDualRailStateEncoding);
    }
    releasePdrBatchAllocatorPages();
    if (isSecDiagEnabled()) {
      emitSecDiag(
          "SEC diag: PDR steady-state check end index=",
          batchIndex,
          " output_range=",
          firstOutput,
          "..",
          endOutput,
          " status=",
          pdrStatusName(pdrResult.status),
          " bound=",
          pdrResult.bound);
    }
    switch (pdrResult.status) {
      case PDRStatus::Equivalent:
        provedBound = std::max(provedBound, pdrResult.bound);
        markPdrOutputRangeCovered(
            pdrCoveredOutputs,
            pdrSkippedOutputReasons,
            firstOutput,
            endOutput);
        break;
      case PDRStatus::Different:
        return makeSecResult(
            SequentialEquivalenceStatus::Different,
            pdrResult.bound,
            "Exact PDR found a defined-value counterexample at k = " +
                std::to_string(pdrResult.bound),
            outputCoverage,
            abstractedSequentialBoundaries,
            extractedBoundaryReports);
      case PDRStatus::Inconclusive:
      default:
        provedBound = std::max(provedBound, pdrResult.bound);
        if (problem.usesDualRailStateEncoding) {
          if (endOutput - firstOutput > 1) {
            const size_t midOutput =
                firstOutput + (endOutput - firstOutput) / 2;
            outputBatches.insert(
                outputBatches.begin() +
                    static_cast<std::ptrdiff_t>(batchIndex + 1),
                {PdrOutputBatch{firstOutput, midOutput},
                 PdrOutputBatch{midOutput, endOutput}});
            break;
          }
          markPdrOutputRangeSkipped(
              pdrCoveredOutputs,
              pdrSkippedOutputReasons,
              firstOutput,
              endOutput,
              "dual-rail PDR steady-state proof was inconclusive");
          break;
        }
        for (size_t outputIndex = 0;
             outputIndex < pdrCoveredOutputs.size();
             ++outputIndex) {
          if (!pdrCoveredOutputs[outputIndex]) {
            pdrSkippedOutputReasons.emplace(
                outputIndex, "exact PDR proof was inconclusive");
          }
        }
        stopAfterInconclusiveBatch = true;
        break;
    }
    if (stopAfterInconclusiveBatch) {
      break;
    }
  }

  {
    const OutputCoverageSelection finalCoverage =
        buildCoverageWithDualRailOutputSkips(
            outputCoverage, problem, pdrCoveredOutputs, pdrSkippedOutputReasons);
    const size_t coveredOutputCount = static_cast<size_t>(
        std::count(pdrCoveredOutputs.begin(), pdrCoveredOutputs.end(), true));
    if (finalCoverage.checkedOutputs.names.empty()) {
      const OutputCoverageSelection& noProofCoverage =
          problem.usesDualRailStateEncoding ? finalCoverage : outputCoverage;
      return makeSecResult(
          SequentialEquivalenceStatus::Inconclusive,
          provedBound,
          problem.usesDualRailStateEncoding
              ? "Exact dual-rail PDR did not prove any observed output"
              : "Exact PDR did not prove any observed output",
          noProofCoverage,
          abstractedSequentialBoundaries,
          extractedBoundaryReports);
    }
    if (coveredOutputCount != pdrCoveredOutputs.size()) {
      return makeSecResult(
          SequentialEquivalenceStatus::PartiallyProved,
          provedBound,
          std::string("Exact ") +
              (problem.usesDualRailStateEncoding ? "dual-rail " : "") +
              "PDR proved " +
              std::to_string(coveredOutputCount) + " of " +
              std::to_string(pdrCoveredOutputs.size()) +
              " observed outputs; remaining outputs are inconclusive",
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

SequentialEquivalenceResult finishDualRailImcProof(
    const KInductionProblem& problem,
    const IMCResult& guardedResult,
    const OutputCoverageSelection& outputCoverage,
    const std::vector<std::string>& abstractedSequentialBoundaries,
    const std::vector<ExtractedBoundaryReportEntry>& extractedBoundaryReports) {
  DualRailResidualProofState proofState;
  proofState.coveredOutputs.assign(
      problem.observedOutputExprs0.size(), false);
  proofState.provedBound = guardedResult.bound;

  size_t guardedProvedPrefix = 0;
  if (guardedResult.status == IMCStatus::Equivalent) {
    guardedProvedPrefix = problem.observedOutputExprs0.size();
  } else if (guardedResult.firstUnprovenOutput.has_value()) {
    guardedProvedPrefix = std::min(
        *guardedResult.firstUnprovenOutput,
        problem.observedOutputExprs0.size());
  }

  std::vector<size_t> provedOutputIndices;
  provedOutputIndices.reserve(guardedProvedPrefix);
  for (size_t outputIndex = 0;
       outputIndex < problem.observedOutputExprs0.size();
       ++outputIndex) {
    if (problem.dualRailOutputSkipReasons.size() ==
            problem.observedOutputExprs0.size() &&
        !problem.dualRailOutputSkipReasons[outputIndex].empty()) {
      proofState.skipReasons.emplace(
          outputIndex, problem.dualRailOutputSkipReasons[outputIndex]);
      continue;
    }
    if (outputIndex < guardedProvedPrefix) {
      provedOutputIndices.push_back(outputIndex);
    } else {
      proofState.skipReasons.emplace(
          outputIndex,
          "dual-rail IMC steady-state proof was inconclusive");
    }
  }

  markDualRailResidualOutputsCovered(provedOutputIndices, proofState);

  const size_t coveredCount = static_cast<size_t>(std::count(
      proofState.coveredOutputs.begin(),
      proofState.coveredOutputs.end(),
      true));
  const OutputCoverageSelection finalCoverage =
      buildCoverageWithDualRailOutputSkips(
          outputCoverage,
          problem,
          proofState.coveredOutputs,
          proofState.skipReasons);
  SequentialEquivalenceStatus status = SequentialEquivalenceStatus::Equivalent;
  std::string reason;
  if (coveredCount == 0) {
    status = SequentialEquivalenceStatus::Inconclusive;
    reason = "Dual-rail IMC did not prove any observed output";
  } else if (coveredCount != proofState.coveredOutputs.size()) {
    status = SequentialEquivalenceStatus::PartiallyProved;
    reason =
        "Dual-rail IMC proved " + std::to_string(coveredCount) + " of " +
        std::to_string(proofState.coveredOutputs.size()) +
        " observed outputs; remaining outputs are inconclusive";
  }

  SequentialEquivalenceResult secResult = makeSecResult(
      status,
      proofState.provedBound,
      std::move(reason),
      finalCoverage,
      abstractedSequentialBoundaries,
      extractedBoundaryReports);
  setExactSecEngineProofProgress(
      secResult, problem, "IMC", proofState.coveredOutputs);
  return secResult;
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
  if (problem.usesDualRailStateEncoding &&
      result.status != IMCStatus::Different) {
    return finishDualRailImcProof(
        problem,
        result,
        outputCoverage,
        abstractedSequentialBoundaries,
        extractedBoundaryReports);
  }
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
      const size_t provenOutputCount =  // LCOV_EXCL_LINE
          result.firstUnprovenOutput.value_or(0);  // LCOV_EXCL_LINE
      const SequentialEquivalenceStatus status =  // LCOV_EXCL_LINE
          provenOutputCount > 0  // LCOV_EXCL_LINE
              ? SequentialEquivalenceStatus::PartiallyProved  // LCOV_EXCL_LINE
              : SequentialEquivalenceStatus::Inconclusive;  // LCOV_EXCL_LINE
      SequentialEquivalenceResult secResult = makeSecResult(  // LCOV_EXCL_LINE
          status,
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
      // Defensive fallback for corrupted enum values; public parsing rejects
      // unknown SEC engines before strategy construction.
      // LCOV_EXCL_START
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
      // LCOV_EXCL_STOP
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

  // Phase 3: rewrite both designs into one shared symbol space and preserve the
  // exact extracted initial predicate as IC3/PDR's F[0].
  SharedSecSymbolSpace symbolSpace = buildSharedSecSymbolSpace(
      model0, model1, aligned.inputs, aligned.outputs);
  integrateInitialState(
      model0,
      model1,
      symbolSpace.state0Symbols,
      symbolSpace.state1Symbols,
      symbolSpace.problem);
  if (encoding_ == SecEncoding::Binary) {
    const bool hasIncompleteInitialState =
        model0.initialStateValueByKey.size() != model0.stateBits.size() ||
        model1.initialStateValueByKey.size() != model1.stateBits.size();
    filterOutputsRequiringUninitializedState(
        model0,
        model1,
        hasIncompleteInitialState,
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
  // use it. Output remapping happens after the startup coverage filter so
  // skipped top outputs never allocate SAT symbols or proof obligations.
  const bool useLazyTransitionRemapping =
      secEngine_ == SecEngine::KInduction || secEngine_ == SecEngine::Pdr;
  KInductionProblem proofProblem;
  if (encoding_ == SecEncoding::DualRailSteady) {
    proofProblem = buildDualRailSecProblem(
        model0,
        model1,
        aligned.inputs,
        aligned.outputs,
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
