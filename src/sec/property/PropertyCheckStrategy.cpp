// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "property/PropertyCheckStrategy.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "BoolExpr.h"
#include "NLDB0.h"
#include "SNLBusTerm.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLRTLInfos.h"
#include "clocks/SecClockModel.h"
#include "model/SequentialDesignModel.h"
#include "property/BoolExprOps.h"
#include "property/PropertyExpression.h"

namespace KEPLER_FORMAL::SEC {

namespace {

using namespace PropertyLogic;
using naja::NL::NLDB0;
using naja::NL::NLName;
using naja::NL::SNLDesign;
using naja::NL::SNLInstance;

// Synthetic keys, next to the boundary (bit 62) and X-constant (bit 61) tags
// of the extracted model.
constexpr uint64_t kPropertyStateTag = uint64_t{1} << 60;
constexpr uint64_t kPropertyOutputTag = uint64_t{1} << 59;
constexpr uint64_t kPropertyInitialTag = uint64_t{1} << 58;
// One state bit per reset cycle keeps the reset sequence trivially correct.
constexpr size_t kMaxResetCycles = 4096;
// Display name of the free value a register without initial value starts at.
constexpr const char* kInitialValuePrefix = "__kepler_initial(";

SignalKey propertyStateKey(size_t index) {
  return {{kPropertyStateTag, static_cast<uint64_t>(index)}, {0}};
}

SignalKey propertyOutputKey(size_t index) {
  return {{kPropertyOutputTag, static_cast<uint64_t>(index)}, {0}};
}

SignalKey propertyInitialKey(size_t index) {
  return {{kPropertyInitialTag, static_cast<uint64_t>(index)}, {0}};
}

struct PropertySite {
  // Top-relative path in the boundary path syntax, and its readable form.
  std::string path;
  std::string displayName;
  SNLInstance* instance = nullptr;
  PropertyKind kind = PropertyKind::Assert;
  // "file:line" of the property, empty when the frontend recorded none.
  std::string source;
  // Number of DATA bits.
  size_t width = 0;
};

std::string escapePathComponent(const std::string& name) {
  std::string escaped;
  escaped.reserve(name.size());
  for (const char character : name) {
    if (character == '\\' || character == '/') {
      escaped += '\\';
    }
    escaped += character;
  }
  return escaped;
}

bool containsProperties(const SNLDesign* design,
                        std::unordered_map<const SNLDesign*, bool>& memo) {
  if (const auto found = memo.find(design); found != memo.end()) {
    return found->second;
  }
  memo.emplace(design, false);
  bool result = false;
  for (const auto* instance : design->getInstances()) {
    const auto* model = instance->getModel();
    if (NLDB0::isFormalProperty(model) ||
        (!model->isLeaf() && containsProperties(model, memo))) {
      result = true;
      break;
    }
  }
  memo[design] = result;
  return result;
}

PropertyKind toPropertyKind(NLDB0::FormalProperty::Kind kind) {
  switch (kind) {
    case NLDB0::FormalProperty::Kind::Assert:
      return PropertyKind::Assert;
    case NLDB0::FormalProperty::Kind::Assume:
      return PropertyKind::Assume;
    case NLDB0::FormalProperty::Kind::Cover:
      return PropertyKind::Cover;
  }
  return PropertyKind::Assert;  // LCOV_EXCL_LINE
}

void collectProperties(SNLDesign* design,
                       const std::string& pathPrefix,
                       const std::string& namePrefix,
                       std::unordered_map<const SNLDesign*, bool>& memo,
                       std::vector<PropertySite>& sites,
                       std::vector<std::string>& unaddressable) {
  for (auto* instance : design->getInstances()) {
    auto* model = instance->getModel();
    const bool isProperty = NLDB0::isFormalProperty(model);
    if (!isProperty && (model->isLeaf() || !containsProperties(model, memo))) {
      continue;
    }
    if (!isProperty) {
      if (instance->isUnnamed()) {
        unaddressable.push_back(namePrefix + "<unnamed instance of " +
                                model->getName().getString() + ">");
        continue;
      }
      const auto name = instance->getName().getString();
      collectProperties(model, pathPrefix + escapePathComponent(name) + "/",
                        namePrefix + name + "/", memo, sites, unaddressable);
      continue;
    }

    PropertySite site;
    site.instance = instance;
    site.kind = toPropertyKind(NLDB0::getFormalProperty(instance).kind);
    site.width = NLDB0::getFormalPropertyWidth(instance);
    uint32_t line = 0;
    if (const auto* infos = instance->getRTLInfos();
        infos != nullptr && infos->hasSourceLoc()) {
      const auto& location = *infos->getSourceLoc();
      line = location.line;
      site.source = location.file.getString() + ":" + std::to_string(line);
    }
    if (instance->isUnnamed()) {
      // An assertion without a label has no name to report or to select it
      // by: name it after its kind and line, once for the whole module.
      const std::string base =
          std::string(propertyKindName(site.kind)) + "_line" + std::to_string(line);
      std::string candidate = base;
      for (size_t suffix = 1; design->getInstance(NLName(candidate)) != nullptr;
           ++suffix) {
        candidate = base + "_" + std::to_string(suffix);
      }
      instance->setName(NLName(candidate));
    }
    const auto name = instance->getName().getString();
    site.path = pathPrefix + escapePathComponent(name);
    site.displayName = namePrefix + name;
    sites.push_back(std::move(site));
  }
}

struct MonitorState {
  SignalKey key;
  size_t var = 0;
  BoolExpr* next = nullptr;
  std::string displayName;
};

struct PreparedProperty {
  PropertyResult result;
  bool usable = false;
  // Assert and assume: 1 in a step where a violation is detected.
  // Cover: 1 in a step where the cover is hit.
  BoolExpr* signal = nullptr;
  // Attempts in progress: cleared while the property is not active.
  std::vector<MonitorState> monitor;
  // Past values the property reads: these registers always run.
  std::vector<MonitorState> history;
};

// Every state bit of a checked model has an initial value: the engines then
// start from that one state and assume nothing about the first step.
void addState(SequentialDesignModel& model,
              const MonitorState& state,
              bool initialValue) {
  model.stateBits.push_back(state.key);
  model.inputVarByKey[state.key] = state.var;
  model.displayNameByKey[state.key] = state.displayName;
  model.nextStateExprByStateKey[state.key] = state.next;
  model.initialStateValueByKey[state.key] = initialValue;
}

// The single clock a property may be sampled on: the positive phase of the
// clock that drives the design state. Returns the reason when it is not.
std::optional<std::string> describeClockProblem(const SequentialDesignModel& model,
                                                BoolExpr* clock) {
  bool inverted = false;
  while (clock->getOp() == Op::NOT) {
    clock = clock->getLeft();
    inverted = !inverted;
  }
  if (clock->getOp() != Op::VAR || clock == constTrue() || clock == constFalse()) {
    return "the property clock is not a plain clock signal";
  }

  std::vector<SignalKey> stateDomains;
  for (const auto& [key, event] : model.clockEventByStateKey) {
    if (std::find(stateDomains.begin(), stateDomains.end(), event.domain) ==
        stateDomains.end()) {
      stateDomains.push_back(event.domain);
    }
  }
  const auto carrier = std::find_if(
      model.clockCarrierClasses.begin(), model.clockCarrierClasses.end(),
      [&](const ClockCarrierClass& candidate) {
        return candidate.varID == clock->getId();
      });
  if (carrier == model.clockCarrierClasses.end()) {
    if (stateDomains.empty()) {
      // No clocked state: every step is a tick of whatever clock is named.
      return inverted ? std::optional<std::string>(
                            "properties sampled on the falling clock edge are "
                            "not supported yet")
                      : std::nullopt;
    }
    return "the property clock is not the clock of the design state";
  }
  const bool negative = (carrier->phase == ClockPhase::Neg) != inverted;
  if (negative) {
    return "properties sampled on the falling clock edge are not supported yet";
  }
  for (const auto& domain : stateDomains) {
    if (domain != carrier->domain) {
      return "the design has state in another clock domain than the property "
             "clock";
    }
  }
  return std::nullopt;
}

std::string joinReasons(const std::vector<std::string>& reasons) {
  std::string joined;
  for (const auto& reason : reasons) {
    if (!joined.empty()) {
      joined += "; ";
    }
    joined += reason;
  }
  return joined;
}

// The equivalence engines report a trace as text: an "Input trace:" block with
// one "  cycle <n>: name=value, ..." line per step. Rewrite it for a property:
// reset cycles are labelled as such, later cycles count from the first one
// after reset, clocks are dropped, and the free initial value of a register
// is shown once, in the first step.
std::string formatInputTrace(const std::string& reason,
                             size_t resetCycles,
                             const std::vector<std::string>& hiddenNames) {
  const auto begin = reason.find("Input trace:");
  if (begin == std::string::npos) {
    return "No input trace is available.";
  }
  const auto end = reason.find("Observed output mismatches", begin);
  const std::string block = reason.substr(
      begin, end == std::string::npos ? std::string::npos : end - begin);

  std::string formatted;
  size_t position = 0;
  while (position < block.size()) {
    auto lineEnd = block.find('\n', position);
    if (lineEnd == std::string::npos) {
      lineEnd = block.size();
    }
    const std::string line = block.substr(position, lineEnd - position);
    position = lineEnd + 1;

    const std::string prefix = "  cycle ";
    const auto colon = line.find(':');
    if (line.compare(0, prefix.size(), prefix) != 0 ||
        colon == std::string::npos) {
      if (!line.empty()) {
        formatted += line + "\n";
      }
      continue;
    }
    size_t cycle = 0;
    try {
      cycle = static_cast<size_t>(
          std::stoull(line.substr(prefix.size(), colon - prefix.size())));
    } catch (const std::exception&) {
      formatted += line + "\n";
      continue;
    }
    std::string values;
    size_t tokenStart = colon + 1;
    while (tokenStart < line.size()) {
      auto tokenEnd = line.find(", ", tokenStart);
      if (tokenEnd == std::string::npos) {
        tokenEnd = line.size();
      }
      std::string token = line.substr(tokenStart, tokenEnd - tokenStart);
      tokenStart = tokenEnd + 2;
      while (!token.empty() && token.front() == ' ') {
        token.erase(token.begin());
      }
      const std::string name = token.substr(0, token.rfind('='));
      const bool isInitial = name.rfind(kInitialValuePrefix, 0) == 0;
      if (token.empty() || (isInitial && cycle != 0) ||
          std::find(hiddenNames.begin(), hiddenNames.end(), name) !=
              hiddenNames.end()) {
        continue;
      }
      if (isInitial) {
        // "__kepler_initial(reg)=1" reads better as "initial reg=1".
        const size_t open = std::string(kInitialValuePrefix).size();
        const auto close = token.rfind(')');
        if (close != std::string::npos && close > open) {
          token = "initial " + token.substr(open, close - open) +
                  token.substr(close + 1);
        }
      }
      values += (values.empty() ? " " : ", ") + token;
    }
    formatted += cycle < resetCycles
                     ? "  reset cycle " + std::to_string(cycle) + ":"
                     : "  cycle " + std::to_string(cycle - resetCycles) + ":";
    formatted += values + "\n";
  }
  while (!formatted.empty() && formatted.back() == '\n') {
    formatted.pop_back();
  }
  return formatted;
}

}  // namespace

const char* propertyKindName(PropertyKind kind) {
  switch (kind) {
    case PropertyKind::Assert:
      return "assert";
    case PropertyKind::Assume:
      return "assume";
    case PropertyKind::Cover:
      return "cover";
  }
  return "assert";  // LCOV_EXCL_LINE
}

const char* propertyStatusName(PropertyStatus status) {
  switch (status) {
    case PropertyStatus::Proved:
      return "proved";
    case PropertyStatus::Failed:
      return "failed";
    case PropertyStatus::Covered:
      return "covered";
    case PropertyStatus::Unreachable:
      return "unreachable";
    case PropertyStatus::Assumed:
      return "assumed";
    case PropertyStatus::Inconclusive:
      return "inconclusive";
    case PropertyStatus::Unsupported:
      return "unsupported";
  }
  return "unsupported";  // LCOV_EXCL_LINE
}

size_t PropertyCheckResult::count(PropertyStatus status) const {
  return static_cast<size_t>(std::count_if(
      properties.begin(), properties.end(),
      [status](const PropertyResult& property) {
        return property.status == status;
      }));
}

PropertyCheckStrategy::PropertyCheckStrategy(
    naja::NL::SNLDesign* top,
    PropertyTranslator translator,
    KEPLER_FORMAL::Config::SolverType solverType,
    SecEngine engine,
    SecResetSpec resetSpec)
    : top_(top),
      translator_(std::move(translator)),
      solverType_(solverType),
      engine_(engine),
      resetSpec_(std::move(resetSpec)) {}

PropertyCheckResult PropertyCheckStrategy::run(size_t maxK) const {
  PropertyCheckResult result;
  if (top_ == nullptr) {
    result.error = "no top design";
    return result;
  }
  const bool diag = std::getenv("KEPLER_FPV_DIAG") != nullptr;

  std::vector<PropertySite> sites;
  std::vector<std::string> unaddressable;
  {
    std::unordered_map<const SNLDesign*, bool> memo;
    collectProperties(top_, "", "", memo, sites, unaddressable);
  }
  for (const auto& name : unaddressable) {
    PropertyResult property;
    property.name = name;
    property.detail =
        "properties below an unnamed instance cannot be addressed";
    result.properties.push_back(std::move(property));
  }
  if (sites.empty()) {
    return result;
  }

  // Every captured property is a leaf whose input pins the extraction
  // publishes as observed signals, like a user-selected boundary instance.
  BoundaryPairs pairs;
  pairs.reserve(sites.size());
  for (const auto& site : sites) {
    pairs.emplace_back(site.path, site.path);
  }
  SequentialDesignModel model;
  std::vector<BoundaryPort> ports;
  try {
    ports = BoundarySelection(top_, pairs, 0).getPorts();
    model = SequentialDesignModel::extract(top_, pairs, 0);
  } catch (const std::exception& error) {
    result.error = error.what();
    return result;
  }
  if (model.hasUnsupportedFeatures()) {
    result.error = joinReasons(model.unsupportedReasons);
    return result;
  }

  std::map<std::tuple<size_t, std::string, int32_t>, SignalKey> keyByPin;
  for (const auto& port : ports) {
    keyByPin.emplace(std::make_tuple(port.pairIndex, port.pinName, port.bit),
                     SequentialDesignModel::boundarySignalKey(port));
  }
  const auto pinValue = [&](size_t site, const char* pin, int32_t bit,
                            std::string& problem) -> BoolExpr* {
    const auto key = keyByPin.find(std::make_tuple(site, std::string(pin), bit));
    if (key == keyByPin.end()) {
      problem = std::string("property pin ") + pin + " was not found";
      return nullptr;
    }
    if (const auto expr = model.observedOutputExprByKey.find(key->second);
        expr != model.observedOutputExprByKey.end() && expr->second != nullptr) {
      return expr->second;
    }
    const auto skip = model.connectivitySkipInfoByKey.find(key->second);
    problem = std::string("the logic driving ") + pin +
              (skip == model.connectivitySkipInfoByKey.end()
                   ? " is not modelled"
                   : " is not modelled: " + skip->second.detail);
    return nullptr;
  };

  // The reset sequence and the initial state are part of the checked model, so
  // the engines see an ordinary transition system with one fully defined
  // initial state and no assumption about any step:
  //
  //   step 0 .. N-1   reset ports held active (N = reset cycles)
  //   step N ..       reset ports held inactive; properties are checked
  //
  // A register without an initial value may start anywhere: in step 0 its
  // value is a fresh input.
  size_t nextVar = model.nextSyntheticVarID();
  size_t syntheticStateCount = 0;
  const auto newState = [&](const std::string& displayName, BoolExpr* next) {
    MonitorState state;
    state.key = propertyStateKey(syntheticStateCount++);
    state.var = nextVar++;
    state.next = next;
    state.displayName = displayName;
    return state;
  };

  std::vector<MonitorState> sequenceState;
  sequenceState.push_back(newState("__kepler_started[0]", constTrue()));
  BoolExpr* started = BoolExpr::Var(sequenceState.front().var);
  BoolExpr* inReset = constFalse();
  const size_t resetCycles = resetSpec_.enabled() ? resetSpec_.cycles : 0;
  if (resetCycles > kMaxResetCycles) {
    result.error = "a reset window longer than " +
                   std::to_string(kMaxResetCycles) +
                   " cycles is not supported";
    return result;
  }
  if (resetCycles != 0) {
    // `reached` is 1 from step `cycle` on.
    BoolExpr* reached = started;
    for (size_t cycle = 2; cycle <= resetCycles; ++cycle) {
      sequenceState.push_back(newState(
          "__kepler_reset_stage_" + std::to_string(cycle) + "[0]", reached));
      reached = BoolExpr::Var(sequenceState.back().var);
    }
    inReset = notExpr(reached);
  }

  std::unordered_map<size_t, BoolExpr*> replacements;
  std::unordered_set<SignalKey, SignalKeyHash> resetPortKeys;
  for (const auto& port : resetSpec_.ports) {
    if (resetCycles == 0) {
      break;
    }
    const SignalKey* found = nullptr;
    for (const auto& key : model.environmentInputs) {
      const auto name = model.displayNameByKey.find(key);
      if (name != model.displayNameByKey.end() &&
          (name->second == port.name || name->second == port.name + "[0]") &&
          model.inputVarByKey.find(key) != model.inputVarByKey.end()) {
        found = &key;
        break;
      }
    }
    if (found == nullptr) {
      result.error = "reset port `" + port.name +
                     "` was not found among the top-level inputs";
      return result;
    }
    resetPortKeys.insert(*found);
    replacements[model.inputVarByKey.at(*found)] =
        port.activeValue ? inReset : notExpr(inReset);
  }

  struct InitialInput {
    SignalKey key;
    size_t var = 0;
    std::string displayName;
  };
  std::vector<InitialInput> initialInputs;
  std::unordered_map<SignalKey, bool, SignalKeyHash> initialValues =
      model.initialStateValueByKey;
  {
    // A complemented output starts as the inverse of its primary one.
    std::unordered_map<SignalKey, SignalKey, SignalKeyHash> primaryOf;
    for (const auto& relation : model.complementedStateRelations) {
      primaryOf.emplace(relation.complementedKey, relation.primaryKey);
    }
    std::unordered_map<SignalKey, BoolExpr*, SignalKeyHash> freeValueByKey;
    const auto freeValue = [&](const SignalKey& key) {
      if (const auto found = freeValueByKey.find(key);
          found != freeValueByKey.end()) {
        return found->second;
      }
      InitialInput input;
      input.key = propertyInitialKey(initialInputs.size());
      input.var = nextVar++;
      // Prefer the name of the register's output net: the instance behind
      // an RTL register has no source name.
      std::string registerName = "?";
      if (const auto net = model.outputNetNameByKey.find(key);
          net != model.outputNetNameByKey.end() && !net->second.empty()) {
        registerName = net->second;
      } else if (const auto name = model.displayNameByKey.find(key);
                 name != model.displayNameByKey.end()) {
        registerName = name->second;
      }
      input.displayName =
          std::string(kInitialValuePrefix) + registerName + ")";
      BoolExpr* value = BoolExpr::Var(input.var);
      initialInputs.push_back(std::move(input));
      freeValueByKey.emplace(key, value);
      return value;
    };
    for (const auto& key : model.stateBits) {
      const auto var = model.inputVarByKey.find(key);
      if (var == model.inputVarByKey.end() ||
          initialValues.find(key) != initialValues.end()) {
        continue;
      }
      BoolExpr* initial = nullptr;
      if (const auto primary = primaryOf.find(key); primary != primaryOf.end()) {
        if (const auto value = initialValues.find(primary->second);
            value != initialValues.end()) {
          initialValues[key] = !value->second;
          continue;
        }
        initial = notExpr(freeValue(primary->second));
      } else {
        initial = freeValue(key);
      }
      replacements[var->second] =
          orExpr(andExpr(started, BoolExpr::Var(var->second)),
                 andExpr(notExpr(started), initial));
    }
  }
  std::unordered_map<BoolExpr*, BoolExpr*> substitutionMemo;
  const auto substitute = [&](BoolExpr* expr) {
    return substituteBoolExprVariableExpressions(expr, replacements,
                                                 substitutionMemo);
  };

  std::vector<PreparedProperty> prepared;
  prepared.reserve(sites.size());
  for (size_t siteIndex = 0; siteIndex < sites.size(); ++siteIndex) {
    const auto& site = sites[siteIndex];
    PreparedProperty property;
    property.result.name = site.displayName;
    property.result.source = site.source;
    property.result.kind = site.kind;
    const auto unsupported = [&](const std::string& reason) {
      property.result.status = PropertyStatus::Unsupported;
      property.result.detail = reason;
      prepared.push_back(std::move(property));
    };
    if (!translator_) {
      unsupported("no translator is available for this property");
      continue;
    }

    // The property reads the design through the same substitution as the
    // design itself.
    PropertyTranslationRequest request;
    request.instance = site.instance;
    std::string problem;
    const std::string dataName =
        NLDB0::getFormalPropertyData(site.instance->getModel())
            ->getName()
            .getString();
    for (size_t bit = 0; problem.empty() && bit < site.width; ++bit) {
      BoolExpr* value = pinValue(siteIndex, dataName.c_str(),
                                 static_cast<int32_t>(bit), problem);
      request.data.push_back(value != nullptr ? substitute(value) : nullptr);
    }
    if (!problem.empty()) {
      unsupported(problem);
      continue;
    }
    // A past value starts free, like a register without initial value, and is
    // sampled in every step: neither reset nor disable stops it.
    request.previous = [&](BoolExpr* bit) {
      MonitorState state;
      state.key = propertyStateKey(syntheticStateCount);
      state.var = nextVar++;
      state.next = bit;
      state.displayName = "__kepler_property_state_" +
                          std::to_string(syntheticStateCount) + "[0]";
      ++syntheticStateCount;
      InitialInput input;
      input.key = propertyInitialKey(initialInputs.size());
      input.var = nextVar++;
      input.displayName = std::string(kInitialValuePrefix) + site.displayName +
                          " past value " +
                          std::to_string(property.history.size()) + ")";
      BoolExpr* value =
          orExpr(andExpr(started, BoolExpr::Var(state.var)),
                 andExpr(notExpr(started), BoolExpr::Var(input.var)));
      initialInputs.push_back(std::move(input));
      property.history.push_back(std::move(state));
      return value;
    };
    const TranslatedProperty translated = translator_(request);
    if (!translated.unsupportedReason.empty()) {
      unsupported(translated.unsupportedReason);
      continue;
    }
    if (auto clockProblem = describeClockProblem(model, translated.clock)) {
      unsupported(*clockProblem);
      continue;
    }

    PropertyTerm term;
    try {
      term = parsePropertyTerm(translated.expression);
    } catch (const std::exception& error) {
      unsupported(error.what());
      continue;
    }
    const auto& leaves = translated.leaves;
    PropertyCompiler compiler(leaves, [&nextVar]() { return nextVar++; });
    const CompiledProperty compiled =
        property.result.kind == PropertyKind::Cover
            ? compiler.compileCover(term)
            : compiler.compileViolation(term);
    if (!compiled.supported()) {
      unsupported(compiled.unsupportedReason);
      continue;
    }
    if (diag) {
      fprintf(stderr,
              "FPV diag: property=%s kind=%s leaves=%zu monitor_bits=%zu "
              "expr=%s\n",
              site.displayName.c_str(), propertyKindName(property.result.kind),
              leaves.size(), compiled.monitor.size(),
              translated.expression.c_str());
      fflush(stderr);
    }
    // A property makes no claim during the reset window or while disabled,
    // and drops its pending attempts then.
    BoolExpr* active = andExpr(notExpr(inReset), notExpr(translated.disable));
    property.signal = BoolExpr::simplify(andExpr(compiled.signal, active));
    for (const auto& bit : compiled.monitor) {
      MonitorState state;
      state.key = propertyStateKey(syntheticStateCount);
      state.var = bit.var;
      state.next = BoolExpr::simplify(andExpr(active, bit.next));
      state.displayName = "__kepler_property_state_" +
                          std::to_string(syntheticStateCount) + "[0]";
      ++syntheticStateCount;
      property.monitor.push_back(std::move(state));
    }
    property.usable = true;
    prepared.push_back(std::move(property));
  }

  // Shared part of every check: the design with its reset sequence, without
  // its own outputs, plus the monitors of the assumptions.
  SequentialDesignModel base;
  base.syntheticVarFloor = model.syntheticVarFloor;
  std::vector<std::string> hiddenInputNames;
  for (const auto& key : model.environmentInputs) {
    const auto var = model.inputVarByKey.find(key);
    if (var == model.inputVarByKey.end() ||
        resetPortKeys.find(key) != resetPortKeys.end()) {
      continue;
    }
    base.environmentInputs.push_back(key);
    base.inputVarByKey.emplace(key, var->second);
    if (const auto name = model.displayNameByKey.find(key);
        name != model.displayNameByKey.end()) {
      base.displayNameByKey.emplace(key, name->second);
      if (std::find(model.clockCarrierVarIDs.begin(),
                    model.clockCarrierVarIDs.end(),
                    var->second) != model.clockCarrierVarIDs.end()) {
        hiddenInputNames.push_back(name->second);
      }
    }
  }
  for (const auto& input : initialInputs) {
    base.environmentInputs.push_back(input.key);
    base.inputVarByKey.emplace(input.key, input.var);
    base.displayNameByKey.emplace(input.key, input.displayName);
  }
  base.topInputKeys = base.environmentInputs;
  for (const auto& key : model.stateBits) {
    const auto var = model.inputVarByKey.find(key);
    const auto next = model.nextStateExprByStateKey.find(key);
    if (var == model.inputVarByKey.end() ||
        next == model.nextStateExprByStateKey.end() || next->second == nullptr) {
      continue;  // LCOV_EXCL_LINE
    }
    MonitorState state;
    state.key = key;
    state.var = var->second;
    state.next = substitute(next->second);
    if (const auto name = model.displayNameByKey.find(key);
        name != model.displayNameByKey.end()) {
      state.displayName = name->second;
    }
    const auto initial = initialValues.find(key);
    addState(base, state, initial != initialValues.end() && initial->second);
  }
  for (const auto& state : sequenceState) {
    addState(base, state, false);
  }

  BoolExpr* assumptionViolated = constFalse();
  size_t unsupportedAssumptions = 0;
  for (auto& property : prepared) {
    if (property.result.kind != PropertyKind::Assume) {
      continue;
    }
    if (!property.usable) {
      ++unsupportedAssumptions;
      continue;
    }
    property.result.status = PropertyStatus::Assumed;
    assumptionViolated = orExpr(assumptionViolated, property.signal);
    for (const auto& state : property.monitor) {
      addState(base, state, false);
    }
    for (const auto& state : property.history) {
      addState(base, state, false);
    }
  }
  // An assumption broken in an earlier step rules the trace out for good.
  BoolExpr* constraintBroken = assumptionViolated;
  if (assumptionViolated != constFalse()) {
    MonitorState sticky = newState("__kepler_assumption_broken[0]", nullptr);
    constraintBroken = orExpr(BoolExpr::Var(sticky.var), assumptionViolated);
    sticky.next = BoolExpr::simplify(andExpr(notExpr(inReset), constraintBroken));
    addState(base, sticky, false);
  }

  for (size_t index = 0; index < prepared.size(); ++index) {
    auto& property = prepared[index];
    if (!property.usable || property.result.kind == PropertyKind::Assume) {
      continue;
    }
    SequentialDesignModel candidate = base;
    for (const auto& state : property.monitor) {
      addState(candidate, state, false);
    }
    for (const auto& state : property.history) {
      addState(candidate, state, false);
    }
    // The checked signal is 1 in a legal step: no violation (or no cover hit)
    // on a trace that kept every assumption so far.
    const SignalKey outputKey = propertyOutputKey(index);
    BoolExpr* bad = andExpr(property.signal, notExpr(constraintBroken));
    const std::string outputName = property.result.name + "[0]";
    candidate.observedOutputs.push_back(outputKey);
    candidate.allObservedOutputs.push_back(outputKey);
    candidate.topOutputKeys.push_back(outputKey);
    candidate.displayNameByKey[outputKey] = outputName;
    candidate.observedOutputExprByKey[outputKey] =
        BoolExpr::simplify(notExpr(bad));

    // The reference side holds the checked signal at constant 1 over the same
    // inputs, so the equivalence engines prove "always 1" for the candidate.
    SequentialDesignModel reference;
    reference.environmentInputs = candidate.environmentInputs;
    reference.topInputKeys = candidate.topInputKeys;
    reference.syntheticVarFloor = candidate.syntheticVarFloor;
    for (const auto& key : candidate.environmentInputs) {
      reference.inputVarByKey.emplace(key, candidate.inputVarByKey.at(key));
      if (const auto name = candidate.displayNameByKey.find(key);
          name != candidate.displayNameByKey.end()) {
        reference.displayNameByKey.emplace(key, name->second);
      }
    }
    reference.observedOutputs.push_back(outputKey);
    reference.allObservedOutputs.push_back(outputKey);
    reference.topOutputKeys.push_back(outputKey);
    reference.displayNameByKey[outputKey] = outputName;
    reference.observedOutputExprByKey[outputKey] = constTrue();

    SequentialEquivalenceResult sec;
    try {
      // No reset specification here: the reset sequence is already in the
      // model, and every state bit has an initial value.
      SequentialEquivalenceStrategy strategy(
          nullptr, nullptr, solverType_, engine_, SecEncoding::Binary,
          SecResetSpec{});
      InternalRelationOptions relationOptions;
      relationOptions.learnInternalRelations = false;
      relationOptions.binaryRelationRound = false;
      strategy.setInternalRelationOptions(relationOptions);
      sec = strategy.runExtractedModels(candidate, reference, maxK);
    } catch (const std::exception& error) {
      property.result.status = PropertyStatus::Inconclusive;
      property.result.detail = error.what();
      continue;
    }

    const bool isCover = property.result.kind == PropertyKind::Cover;
    property.result.bound = sec.bound;
    switch (sec.status) {
      case SequentialEquivalenceStatus::Equivalent:
        property.result.status =
            isCover ? PropertyStatus::Unreachable : PropertyStatus::Proved;
        break;
      case SequentialEquivalenceStatus::Different:
        property.result.status =
            isCover ? PropertyStatus::Covered : PropertyStatus::Failed;
        // Steps are counted from the first cycle after the reset window.
        property.result.bound =
            sec.bound >= resetCycles ? sec.bound - resetCycles : 0;
        property.result.detail =
            formatInputTrace(sec.reason, resetCycles, hiddenInputNames);
        if (unsupportedAssumptions != 0) {
          property.result.detail +=
              "\nNote: " + std::to_string(unsupportedAssumptions) +
              " unsupported assumption(s) were not applied, so this trace may "
              "break them.";
        }
        break;
      case SequentialEquivalenceStatus::Unsupported: {
        property.result.status = PropertyStatus::Unsupported;
        std::vector<std::string> reasons = sec.skippedObservedOutputs;
        if (reasons.empty()) {
          reasons.push_back(sec.reason);
        }
        property.result.detail = joinReasons(reasons);
        break;
      }
      default:
        property.result.status = PropertyStatus::Inconclusive;
        property.result.detail = sec.reason;
        break;
    }
  }

  for (auto& property : prepared) {
    result.properties.push_back(std::move(property.result));
  }
  return result;
}

}  // namespace KEPLER_FORMAL::SEC
