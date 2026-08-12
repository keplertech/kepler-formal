// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "formal/C2RtlEquivalenceStrategy.h"

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "BoolExpr.h"
#include "DNL.h"
#include "NLDB0.h"
#include "NLUniverse.h"
#include "SNLBitTerm.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLPath.h"
#include "clocks/SecClockModel.h"
#include "common/SignalKey.h"
#include "kinduction/KInductionProblem.h"
#include "model/SequentialDesignModel.h"
#include "pdr/PDREngine.h"

namespace KEPLER_FORMAL::C2RTL {

namespace {

using KEPLER_FORMAL::SEC::ClockEvent;
using KEPLER_FORMAL::SEC::ClockPhase;
using KEPLER_FORMAL::SEC::KInductionProblem;
using KEPLER_FORMAL::SEC::PDRResult;
using KEPLER_FORMAL::SEC::PDRStatus;
using KEPLER_FORMAL::SEC::SequentialDesignModel;
using KEPLER_FORMAL::SEC::SignalKey;
using KEPLER_FORMAL::SEC::SignalKeyHash;

class UnsupportedC2Rtl final : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

std::string joinStrings(const std::vector<std::string> &values,
                        const char *separator) {
  std::ostringstream oss;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      oss << separator;
    }
    oss << values[i];
  }
  return oss.str();
}

SignalKey getTerminalPathKey(const naja::DNL::DNLTerminalFull &terminal) {
  SignalKey key;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  key.first.reserve(pathNames.size() + 1);
  for (const auto &name : pathNames) {
    key.first.push_back(name.getID());
  }
  key.first.push_back(terminal.getSnlBitTerm()->getName().getID());
  key.second.push_back(static_cast<naja::NL::NLID::DesignObjectID>(
      terminal.getSnlBitTerm()->getBit()));
  return key;
}

std::string terminalDisplayName(const naja::DNL::DNLTerminalFull &terminal) {
  std::ostringstream oss;
  const auto pathNames = terminal.getDNLInstance().getPath().getPathNames();
  for (const auto &name : pathNames) {
    oss << name.getString() << ".";
  }
  oss << terminal.getSnlBitTerm()->getName().getString() << "["
      << terminal.getSnlBitTerm()->getBit() << "]";
  return oss.str();
}

struct ParsedBitName {
  std::string base;
  size_t bit = 0;
};

ParsedBitName parseBitName(const std::string &name, const char *role) {
  const size_t open = name.rfind('[');
  if (open == std::string::npos || name.empty() || name.back() != ']' ||
      open == 0 || open + 2 > name.size()) {
    throw UnsupportedC2Rtl(std::string("C2RTL cannot parse ") + role +
                           " bit name `" + name + "`; expected name[index]");
  }
  const std::string bitText = name.substr(open + 1, name.size() - open - 2);
  if (bitText.empty() ||
      bitText.find_first_not_of("0123456789") != std::string::npos) {
    throw UnsupportedC2Rtl(std::string("C2RTL cannot parse ") + role +
                           " bit name `" + name +
                           "`; expected a non-negative bit index");
  }
  size_t bit = 0;
  try {
    bit = static_cast<size_t>(std::stoull(bitText));
  } catch (const std::exception &) {
    throw UnsupportedC2Rtl(std::string("C2RTL bit index is out of range in ") +
                           role + " `" + name + "`");
  }
  return {name.substr(0, open), bit};
}

struct SignalBus {
  std::map<size_t, SignalKey> bits;
};

using SignalBuses = std::map<std::string, SignalBus>;

SignalBuses collectSignalBuses(const SequentialDesignModel &model,
                               const std::vector<SignalKey> &keys,
                               const char *role) {
  SignalBuses buses;
  for (const auto &key : keys) {
    const auto nameIt = model.displayNameByKey.find(key);
    if (nameIt == model.displayNameByKey.end()) {
      throw UnsupportedC2Rtl(std::string("C2RTL extracted ") + role +
                             " has no display name");
    }
    const ParsedBitName parsed = parseBitName(nameIt->second, role);
    auto &bits = buses[parsed.base].bits;
    if (!bits.emplace(parsed.bit, key).second) {
      throw UnsupportedC2Rtl(std::string("C2RTL found duplicate ") + role +
                             " bit `" + nameIt->second + "`");
    }
  }
  return buses;
}

size_t contiguousBusWidth(const std::string &name, const SignalBus &bus,
                          const char *role) {
  if (bus.bits.empty() || bus.bits.begin()->first != 0 ||
      bus.bits.rbegin()->first + 1 != bus.bits.size()) {
    throw UnsupportedC2Rtl(
        std::string("C2RTL requires contiguous zero-based ") + role +
        " bits for `" + name + "`");
  }
  return bus.bits.size();
}

std::string displayNameForKey(const SequentialDesignModel &model,
                              const SignalKey &key, const char *role) {
  const auto it = model.displayNameByKey.find(key);
  if (it == model.displayNameByKey.end()) {
    throw UnsupportedC2Rtl(std::string("C2RTL extracted ") + role +
                           " has no display name");
  }
  return it->second;
}

bool containsKey(const std::vector<SignalKey> &keys, const SignalKey &key) {
  return std::find(keys.begin(), keys.end(), key) != keys.end();
}

class ScopedDnl {
public:
  explicit ScopedDnl(naja::NL::SNLDesign *top)
      : universe_(naja::NL::NLUniverse::get()),
        previousTop_(universe_ == nullptr ? nullptr
                                          : universe_->getTopDesign()) {
    if (universe_ == nullptr) {
      throw UnsupportedC2Rtl(
          "C2RTL control classification requires an Naja universe");
    }
    naja::DNL::destroy();
    universe_->setTopDesign(top);
    dnl_ = naja::DNL::get();
    if (dnl_ == nullptr) {
      throw UnsupportedC2Rtl(
          "C2RTL could not elaborate the RTL for control classification");
    }
  }

  ~ScopedDnl() {
    naja::DNL::destroy();
    if (universe_ != nullptr && previousTop_ != nullptr) {
      universe_->setTopDesign(previousTop_);
    }
  }

  naja::DNL::DNLFull *get() const { return dnl_; }

private:
  naja::NL::NLUniverse *universe_ = nullptr;
  naja::NL::SNLDesign *previousTop_ = nullptr;
  naja::DNL::DNLFull *dnl_ = nullptr;
};

struct ControlSource {
  SignalKey key;
  std::string name;
  bool inverted = false;
};

ControlSource traceUnaryControlSource(naja::DNL::DNLFull *dnl,
                                      naja::DNL::DNLID startTermID,
                                      const char *controlRole) {
  std::unordered_set<naja::DNL::DNLID> visited;
  naja::DNL::DNLID currentTermID = startTermID;
  bool inverted = false;

  while (currentTermID != naja::DNL::DNLID_MAX &&
         visited.insert(currentTermID).second) {
    const auto &current = dnl->getDNLTerminalFromID(currentTermID);
    if (current.isNull()) {
      break;
    }
    if (current.isTopPort()) {
      if (current.getSnlBitTerm()->getDirection() !=
          naja::NL::SNLBitTerm::Direction::Input) {
        throw UnsupportedC2Rtl(std::string("C2RTL ") + controlRole +
                               " does not originate at a top-level input");
      }
      return {getTerminalPathKey(current), terminalDisplayName(current),
              inverted};
    }

    if (current.getSnlBitTerm()->getDirection() !=
        naja::NL::SNLBitTerm::Direction::Output) {
      const auto isoID = current.getIsoID();
      if (isoID == naja::DNL::DNLID_MAX) {
        break;
      }
      const auto &iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
      if (iso.isConstant()) {
        throw UnsupportedC2Rtl(
            std::string("C2RTL does not support a constant ") + controlRole);
      }
      if (iso.getDrivers().size() != 1) {
        throw UnsupportedC2Rtl(std::string("C2RTL requires one unambiguous ") +
                               controlRole + " driver");
      }
      currentTermID = iso.getDrivers().front();
      continue;
    }

    const auto *cell = current.getDNLInstance().getSNLModel();
    if (cell == nullptr) {
      break;
    }
    const bool isBuffer = naja::NL::SNLDesignModeling::isBuf(cell) ||
                          naja::NL::NLDB0::isAssign(cell);
    const bool isInverter = naja::NL::SNLDesignModeling::isInv(cell);
    if (!isBuffer && !isInverter) {
      throw UnsupportedC2Rtl(
          std::string("C2RTL ") + controlRole +
          " is generated by unsupported combinational logic at `" +
          current.getDNLInstance().getFullPath() + "`");
    }

    std::vector<naja::NL::SNLBitTerm *> inputs;
    for (auto *input : naja::NL::SNLDesignModeling::getCombinatorialInputs(
             current.getSnlBitTerm())) {
      if (input != nullptr &&
          input->getDirection() != naja::NL::SNLBitTerm::Direction::Output) {
        inputs.push_back(input);
      }
    }
    if (inputs.size() != 1) {
      throw UnsupportedC2Rtl(std::string("C2RTL cannot trace ") + controlRole +
                             " through a non-unary cell");
    }
    const auto &inputTerm =
        current.getDNLInstance().getTerminalFromBitTerm(inputs.front());
    if (inputTerm.isNull()) {
      break;
    }
    inverted ^= isInverter;
    currentTermID = inputTerm.getID();
  }

  throw UnsupportedC2Rtl(std::string("C2RTL cannot trace ") + controlRole +
                         " to one top-level input");
}

struct ResetControl {
  SignalKey key;
  std::string name;
  bool activeHigh = true;
};

struct RtlControls {
  SignalKey clockKey;
  std::string clockName;
  ClockPhase clockPhase = ClockPhase::Pos;
  std::optional<ResetControl> reset;
};

RtlControls classifyRtlControls(naja::NL::SNLDesign *top,
                                const SequentialDesignModel &model) {
  if (model.stateBits.empty()) {
    throw UnsupportedC2Rtl(
        "C2RTL temporal comparison requires clocked RTL state");
  }

  std::optional<ClockEvent> commonClock;
  for (const auto &stateKey : model.stateBits) {
    const auto eventIt = model.clockEventByStateKey.find(stateKey);
    if (eventIt == model.clockEventByStateKey.end()) {
      throw UnsupportedC2Rtl(
          "C2RTL could not classify the clock event for state `" +
          displayNameForKey(model, stateKey, "state") + "`");
    }
    const ClockEvent &event = eventIt->second;
    if (!KEPLER_FORMAL::SEC::clockEventIsUngated(event)) {
      throw UnsupportedC2Rtl(
          "C2RTL does not support gated clocks, clock enables, or stalls");
    }
    if (!commonClock.has_value()) {
      commonClock = event;
    } else if (commonClock->domain != event.domain ||
               commonClock->phase != event.phase) {
      throw UnsupportedC2Rtl(
          "C2RTL requires exactly one clock domain and one active edge");
    }
  }
  if (!commonClock.has_value() ||
      !containsKey(model.topInputKeys, commonClock->domain)) {
    throw UnsupportedC2Rtl(
        "C2RTL clock is generated internally or is not one top-level input");
  }

  const ParsedBitName parsedClock = parseBitName(
      displayNameForKey(model, commonClock->domain, "clock"), "clock");
  const SignalBuses topInputs =
      collectSignalBuses(model, model.topInputKeys, "RTL input");
  const auto clockBusIt = topInputs.find(parsedClock.base);
  if (clockBusIt == topInputs.end() || clockBusIt->second.bits.size() != 1) {
    throw UnsupportedC2Rtl("C2RTL clock must be a scalar top-level input");
  }

  ScopedDnl scopedDnl(top);
  auto *dnl = scopedDnl.get();
  std::unordered_set<SignalKey, SignalKeyHash> requiredStateKeys(
      model.stateBits.begin(), model.stateBits.end());
  std::unordered_set<SignalKey, SignalKeyHash> foundStateKeys;
  std::unordered_set<naja::DNL::DNLID> relevantInstances;
  for (naja::DNL::DNLID termID = 0; termID < dnl->getNBterms(); ++termID) {
    const auto &term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull() || term.isTopPort()) {
      continue;
    }
    const SignalKey key = getTerminalPathKey(term);
    if (requiredStateKeys.find(key) != requiredStateKeys.end()) {
      foundStateKeys.insert(key);
      relevantInstances.insert(term.getDNLInstance().getID());
    }
  }
  if (foundStateKeys.size() != requiredStateKeys.size()) {
    throw UnsupportedC2Rtl("C2RTL cannot map every extracted state bit to a "
                           "modeled RTL flip-flop");
  }

  size_t resettableInstances = 0;
  size_t resetlessInstances = 0;
  std::optional<ResetControl> commonReset;
  for (const auto instanceID : relevantInstances) {
    const auto &instance = dnl->getDNLInstanceFromID(instanceID);
    const auto *primitive = instance.getSNLModel();
    if (primitive == nullptr ||
        !naja::NL::SNLDesignModeling::hasSequentialModel(primitive)) {
      throw UnsupportedC2Rtl(
          "C2RTL encountered RTL state without Naja sequential metadata");
    }
    const auto &sequential =
        naja::NL::SNLDesignModeling::getSequentialModel(primitive);
    if (sequential.kind !=
        naja::NL::SNLDesignModeling::SequentialModel::Kind::FlipFlop) {
      throw UnsupportedC2Rtl("C2RTL does not support latches");
    }

    std::vector<naja::DNL::DNLID> clockPins;
    std::vector<naja::DNL::DNLID> resetPins;
    for (naja::DNL::DNLID termID = instance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= instance.getTermIndexes().second;
         ++termID) {
      const auto &term = dnl->getDNLTerminalFromID(termID);
      if (term.isNull()) {
        continue;
      }
      using Role = naja::NL::SNLDesignModeling::SNLTermRole;
      const Role role =
          naja::NL::SNLDesignModeling::getTermRole(term.getSnlTerm());
      if (role == Role::Clock) {
        clockPins.push_back(termID);
      } else if (role == Role::AsyncReset || role == Role::SyncReset) {
        resetPins.push_back(termID);
      } else if (role == Role::AsyncSet || role == Role::SyncSet) {
        throw UnsupportedC2Rtl("C2RTL does not support set/preset controls");
      } else if (role == Role::Enable) {
        throw UnsupportedC2Rtl(
            "C2RTL does not support register enables, stalls, or handshakes");
      }
    }
    if (clockPins.size() != 1) {
      throw UnsupportedC2Rtl(
          "C2RTL requires exactly one structurally modeled clock pin per "
          "relevant flip-flop");
    }
    if (resetPins.size() > 1) {
      throw UnsupportedC2Rtl("C2RTL does not support compound reset controls");
    }

    bool everyStateCleared = true;
    bool anyStateCleared = false;
    for (const auto &state : sequential.states) {
      if (state.preset.has_value()) {
        throw UnsupportedC2Rtl("C2RTL does not support state preset behavior");
      }
      anyStateCleared |= state.clear.has_value();
      everyStateCleared &= state.clear.has_value();
    }
    if (resetPins.empty()) {
      if (anyStateCleared) {
        throw UnsupportedC2Rtl(
            "C2RTL found reset behavior without a recognized reset pin");
      }
      ++resetlessInstances;
      continue;
    }
    if (!everyStateCleared) {
      throw UnsupportedC2Rtl(
          "C2RTL does not support a reset that covers only part of a state "
          "primitive");
    }

    const auto &resetPin = dnl->getDNLTerminalFromID(resetPins.front());
    const auto activeLevel =
        naja::NL::SNLDesignModeling::getResetActiveLevel(resetPin.getSnlTerm());
    if (activeLevel == naja::NL::SNLDesignModeling::SNLActiveLevel::NA) {
      throw UnsupportedC2Rtl(
          "C2RTL reset polarity is not available in Naja metadata");
    }
    const ControlSource source =
        traceUnaryControlSource(dnl, resetPins.front(), "reset");
    const bool pinActiveHigh =
        activeLevel == naja::NL::SNLDesignModeling::SNLActiveLevel::High;
    const ResetControl reset{
        source.key, source.name,
        static_cast<bool>(pinActiveHigh ^ source.inverted)};
    if (!containsKey(model.topInputKeys, reset.key)) {
      throw UnsupportedC2Rtl("C2RTL reset is not one top-level input");
    }
    if (!commonReset.has_value()) {
      commonReset = reset;
    } else if (commonReset->key != reset.key ||
               commonReset->activeHigh != reset.activeHigh) {
      throw UnsupportedC2Rtl(
          "C2RTL requires one reset source with one polarity for all relevant "
          "state");
    }
    ++resettableInstances;
  }

  if (resettableInstances != 0 && resetlessInstances != 0) {
    throw UnsupportedC2Rtl(
        "C2RTL does not support a reset that covers only part of the relevant "
        "RTL state");
  }

  RtlControls controls;
  controls.clockKey = commonClock->domain;
  controls.clockName = displayNameForKey(model, commonClock->domain, "clock");
  controls.clockPhase = commonClock->phase;
  controls.reset = std::move(commonReset);
  return controls;
}

struct InputAlignment {
  std::string name;
  SignalKey referenceKey;
  SignalKey implementationKey;
};

std::vector<InputAlignment>
alignDataInputs(const SequentialDesignModel &reference,
                const SequentialDesignModel &implementation,
                const RtlControls &controls) {
  const SignalBuses referenceInputs =
      collectSignalBuses(reference, reference.topInputKeys, "reference input");
  SignalBuses implementationInputs = collectSignalBuses(
      implementation, implementation.topInputKeys, "RTL input");

  const ParsedBitName clockName = parseBitName(controls.clockName, "clock");
  auto clockIt = implementationInputs.find(clockName.base);
  if (clockIt == implementationInputs.end() ||
      clockIt->second.bits.size() != 1 ||
      clockIt->second.bits.begin()->second != controls.clockKey) {
    throw UnsupportedC2Rtl(
        "C2RTL could not remove the classified clock from RTL data inputs");
  }
  implementationInputs.erase(clockIt);

  if (controls.reset.has_value()) {
    const ParsedBitName resetName = parseBitName(controls.reset->name, "reset");
    auto resetIt = implementationInputs.find(resetName.base);
    if (resetIt == implementationInputs.end() ||
        resetIt->second.bits.size() != 1 ||
        resetIt->second.bits.begin()->second != controls.reset->key) {
      throw UnsupportedC2Rtl(
          "C2RTL could not remove the classified reset from RTL data inputs");
    }
    implementationInputs.erase(resetIt);
  }

  std::vector<std::string> missingInRtl;
  std::vector<std::string> extraInRtl;
  for (const auto &[name, _] : referenceInputs) {
    if (implementationInputs.find(name) == implementationInputs.end()) {
      missingInRtl.push_back(name);
    }
  }
  for (const auto &[name, _] : implementationInputs) {
    if (referenceInputs.find(name) == referenceInputs.end()) {
      extraInRtl.push_back(name);
    }
  }
  if (!missingInRtl.empty() || !extraInRtl.empty()) {
    std::ostringstream reason;
    reason << "C2RTL data input mapping is not exact";
    if (!missingInRtl.empty()) {
      reason << "; missing in RTL: " << joinStrings(missingInRtl, ", ");
    }
    if (!extraInRtl.empty()) {
      reason << "; unmatched RTL inputs: " << joinStrings(extraInRtl, ", ")
             << " (unrecognized controls require explicit future protocol "
                "support)";
    }
    throw UnsupportedC2Rtl(reason.str());
  }

  std::vector<InputAlignment> aligned;
  for (const auto &[name, referenceBus] : referenceInputs) {
    const auto &implementationBus = implementationInputs.at(name);
    contiguousBusWidth(name, referenceBus, "reference input");
    contiguousBusWidth(name, implementationBus, "RTL input");
    if (referenceBus.bits.size() != implementationBus.bits.size()) {
      throw UnsupportedC2Rtl(
          "C2RTL input width mismatch for `" + name +
          "`: reference=" + std::to_string(referenceBus.bits.size()) +
          ", RTL=" + std::to_string(implementationBus.bits.size()));
    }
    for (const auto &[bit, referenceKey] : referenceBus.bits) {
      const auto implementationBit = implementationBus.bits.find(bit);
      if (implementationBit == implementationBus.bits.end()) {
        throw UnsupportedC2Rtl("C2RTL input bit mapping is not exact for `" +
                               name + "`");
      }
      aligned.push_back({name + "[" + std::to_string(bit) + "]", referenceKey,
                         implementationBit->second});
    }
  }
  return aligned;
}

struct OutputAlignment {
  std::string name;
  size_t delay = 0;
  size_t width = 0;
  std::vector<std::optional<SignalKey>> referenceBits;
  std::vector<std::optional<SignalKey>> implementationBits;
};

std::vector<OutputAlignment>
alignOutputs(const SequentialDesignModel &reference,
             const SequentialDesignModel &implementation,
             const C2RtlEquivalenceOptions &options) {
  const SignalBuses referenceOutputs = collectSignalBuses(
      reference, reference.topOutputKeys, "reference output");
  const SignalBuses implementationOutputs = collectSignalBuses(
      implementation, implementation.topOutputKeys, "RTL output");

  std::vector<std::string> missingInRtl;
  std::vector<std::string> extraInRtl;
  for (const auto &[name, _] : referenceOutputs) {
    if (implementationOutputs.find(name) == implementationOutputs.end()) {
      missingInRtl.push_back(name);
    }
  }
  for (const auto &[name, _] : implementationOutputs) {
    if (referenceOutputs.find(name) == referenceOutputs.end()) {
      extraInRtl.push_back(name);
    }
  }
  if (!missingInRtl.empty() || !extraInRtl.empty()) {
    std::ostringstream reason;
    reason << "C2RTL output mapping is not exact";
    if (!missingInRtl.empty()) {
      reason << "; missing in RTL: " << joinStrings(missingInRtl, ", ");
    }
    if (!extraInRtl.empty()) {
      reason << "; extra RTL outputs: " << joinStrings(extraInRtl, ", ");
    }
    throw UnsupportedC2Rtl(reason.str());
  }

  std::vector<std::string> missingDelays;
  std::vector<std::string> extraDelays;
  for (const auto &[name, _] : referenceOutputs) {
    if (options.outputDelays.find(name) == options.outputDelays.end()) {
      missingDelays.push_back(name);
    }
  }
  for (const auto &[name, _] : options.outputDelays) {
    if (referenceOutputs.find(name) == referenceOutputs.end()) {
      extraDelays.push_back(name);
    }
  }
  if (!missingDelays.empty()) {
    throw UnsupportedC2Rtl(
        "C2RTL output delay setting is required for every compared output; "
        "missing: " +
        joinStrings(missingDelays, ", "));
  }
  if (!extraDelays.empty()) {
    std::sort(extraDelays.begin(), extraDelays.end());
    throw UnsupportedC2Rtl(
        "C2RTL output delay mapping contains unknown outputs: " +
        joinStrings(extraDelays, ", "));
  }

  std::vector<OutputAlignment> aligned;
  aligned.reserve(referenceOutputs.size());
  for (const auto &[name, referenceBus] : referenceOutputs) {
    const auto &implementationBus = implementationOutputs.at(name);
    const size_t referenceWidth =
        contiguousBusWidth(name, referenceBus, "reference output");
    const size_t implementationWidth =
        contiguousBusWidth(name, implementationBus, "RTL output");
    const size_t normalizedWidth =
        std::max(referenceWidth, implementationWidth);

    OutputAlignment output;
    output.name = name;
    output.delay = options.outputDelays.at(name);
    output.width = normalizedWidth;
    output.referenceBits.resize(normalizedWidth);
    output.implementationBits.resize(normalizedWidth);
    for (const auto &[bit, key] : referenceBus.bits) {
      output.referenceBits[bit] = key;
    }
    for (const auto &[bit, key] : implementationBus.bits) {
      output.implementationBits[bit] = key;
    }
    aligned.push_back(std::move(output));
  }
  return aligned;
}

void ensureModelIsUsable(const SequentialDesignModel &model,
                         const char *label) {
  if (model.hasUnsupportedFeatures()) {
    throw UnsupportedC2Rtl(std::string("C2RTL ") + label +
                           " extraction is unsupported: " +
                           joinStrings(model.unsupportedReasons, "; "));
  }
  if (!model.abstractedSequentialBoundaries.empty() ||
      !model.internalBoundaryInputKeys.empty() ||
      !model.internalBoundaryOutputKeys.empty()) {
    throw UnsupportedC2Rtl(std::string("C2RTL ") + label +
                           " contains an abstracted or opaque boundary");
  }
  if (!model.skippedObservedOutputs.empty() ||
      model.observedOutputs.size() != model.topOutputKeys.size()) {
    throw UnsupportedC2Rtl(
        std::string("C2RTL requires complete coverage of every ") + label +
        " top-level output");
  }
  for (const auto &key : model.topOutputKeys) {
    if (model.observedOutputExprByKey.find(key) ==
        model.observedOutputExprByKey.end()) {
      throw UnsupportedC2Rtl(
          std::string("C2RTL is missing a Boolean function for ") + label +
          " output `" + displayNameForKey(model, key, "output") + "`");
    }
  }
}

std::set<size_t> collectFormulaVariables(BoolExpr *root,
                                         const char *formulaRole) {
  if (root == nullptr) {
    throw UnsupportedC2Rtl(std::string("C2RTL has a null ") + formulaRole +
                           " formula");
  }

  std::set<size_t> variables;
  std::unordered_set<BoolExpr *> visited;
  std::vector<BoolExpr *> stack{root};
  while (!stack.empty()) {
    BoolExpr *node = stack.back();
    stack.pop_back();
    if (node == nullptr || !visited.insert(node).second) {
      continue;
    }
    switch (node->getOp()) {
    case Op::VAR:
      if (node->getId() >= 2) {
        variables.insert(node->getId());
      }
      break;
    case Op::NOT:
      stack.push_back(node->getLeft());
      break;
    case Op::AND:
    case Op::OR:
    case Op::XOR:
      stack.push_back(node->getLeft());
      stack.push_back(node->getRight());
      break;
    case Op::NONE:
    default:
      throw UnsupportedC2Rtl(std::string("C2RTL ") + formulaRole +
                             " uses an unsupported Boolean operator");
    }
  }
  return variables;
}

void ensureFeedForwardRtlState(const SequentialDesignModel &model) {
  std::unordered_map<size_t, size_t> stateIndexByVar;
  stateIndexByVar.reserve(model.stateBits.size());
  for (size_t index = 0; index < model.stateBits.size(); ++index) {
    const auto varIt = model.inputVarByKey.find(model.stateBits[index]);
    if (varIt == model.inputVarByKey.end() ||
        !stateIndexByVar.emplace(varIt->second, index).second) {
      throw UnsupportedC2Rtl(
          "C2RTL cannot identify every RTL state variable uniquely");
    }
  }

  std::vector<std::vector<size_t>> dependents(model.stateBits.size());
  std::vector<size_t> dependencyCount(model.stateBits.size(), 0);
  for (size_t target = 0; target < model.stateBits.size(); ++target) {
    const auto nextIt =
        model.nextStateExprByStateKey.find(model.stateBits[target]);
    if (nextIt == model.nextStateExprByStateKey.end()) {
      throw UnsupportedC2Rtl(
          "C2RTL is missing an RTL next-state function for `" +
          displayNameForKey(model, model.stateBits[target], "state") + "`");
    }
    for (const size_t variable :
         collectFormulaVariables(nextIt->second, "RTL next-state")) {
      const auto sourceIt = stateIndexByVar.find(variable);
      if (sourceIt == stateIndexByVar.end()) {
        continue;
      }
      dependents[sourceIt->second].push_back(target);
      ++dependencyCount[target];
    }
  }

  std::vector<size_t> ready;
  ready.reserve(model.stateBits.size());
  for (size_t index = 0; index < dependencyCount.size(); ++index) {
    if (dependencyCount[index] == 0) {
      ready.push_back(index);
    }
  }
  size_t visited = 0;
  while (!ready.empty()) {
    const size_t source = ready.back();
    ready.pop_back();
    ++visited;
    for (const size_t target : dependents[source]) {
      if (--dependencyCount[target] == 0) {
        ready.push_back(target);
      }
    }
  }
  if (visited != model.stateBits.size()) {
    throw UnsupportedC2Rtl(
        "C2RTL requires feed-forward fixed-latency RTL state; state feedback "
        "indicates a register hold, enable, stall, handshake, or other "
        "stateful behavior");
  }
}

BoolExpr *remapFormulaStrict(BoolExpr *root,
                             const std::unordered_map<size_t, size_t> &symbols,
                             std::unordered_map<BoolExpr *, BoolExpr *> &memo,
                             const char *formulaRole) {
  if (root == nullptr) {
    throw UnsupportedC2Rtl(std::string("C2RTL has a null ") + formulaRole +
                           " formula");
  }
  if (const auto it = memo.find(root); it != memo.end()) {
    return it->second;
  }

  struct Visit {
    BoolExpr *expr = nullptr;
    bool childrenVisited = false;
  };
  std::vector<Visit> stack{{root, false}};
  while (!stack.empty()) {
    const Visit visit = stack.back();
    stack.pop_back();
    BoolExpr *node = visit.expr;
    if (memo.find(node) != memo.end()) {
      continue;
    }
    if (node->getOp() == Op::VAR) {
      const size_t local = node->getId();
      if (local < 2) {
        memo.emplace(node, BoolExpr::Var(local));
        continue;
      }
      const auto mapped = symbols.find(local);
      if (mapped == symbols.end()) {
        throw UnsupportedC2Rtl(
            std::string("C2RTL ") + formulaRole +
            " depends on an unpublished or unrecognized symbol v" +
            std::to_string(local));
      }
      memo.emplace(node, BoolExpr::Var(mapped->second));
      continue;
    }
    if (!visit.childrenVisited) {
      stack.push_back({node, true});
      if (node->getRight() != nullptr) {
        stack.push_back({node->getRight(), false});
      }
      if (node->getLeft() != nullptr) {
        stack.push_back({node->getLeft(), false});
      }
      continue;
    }

    BoolExpr *remapped = nullptr;
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
      remapped =
          BoolExpr::Xor(memo.at(node->getLeft()), memo.at(node->getRight()));
      break;
    case Op::NONE:
    default:
      throw UnsupportedC2Rtl(std::string("C2RTL ") + formulaRole +
                             " uses an unsupported Boolean operator");
    }
    memo.emplace(node, remapped);
  }
  return memo.at(root);
}

void addInitialAssignment(KInductionProblem &problem, size_t symbol, bool value,
                          BoolExpr *&initialCondition) {
  problem.initialStateAssignments.emplace_back(symbol, value);
  BoolExpr *literal = BoolExpr::Var(symbol);
  initialCondition =
      BoolExpr::And(initialCondition, value ? literal : BoolExpr::Not(literal));
  ++problem.initializedStateCount;
}

BoolExpr *makeEquality(BoolExpr *lhs, BoolExpr *rhs) {
  return BoolExpr::Not(BoolExpr::Xor(lhs, rhs));
}

struct BuiltC2RtlProblem {
  KInductionProblem problem;
  size_t comparedBits = 0;
};

BuiltC2RtlProblem buildProblem(const SequentialDesignModel &reference,
                               const SequentialDesignModel &implementation,
                               const RtlControls &controls,
                               const std::vector<InputAlignment> &inputs,
                               const std::vector<OutputAlignment> &outputs) {
  BuiltC2RtlProblem built;
  auto &problem = built.problem;
  size_t nextSymbol = 2;

  std::unordered_map<SignalKey, size_t, SignalKeyHash> rtlStateSymbols;
  std::unordered_map<size_t, size_t> referenceSymbols;
  std::unordered_map<size_t, size_t> implementationSymbols;

  for (const auto &key : implementation.stateBits) {
    const size_t symbol = nextSymbol++;
    rtlStateSymbols.emplace(key, symbol);
    implementationSymbols.emplace(implementation.inputVarByKey.at(key), symbol);
    problem.state1Symbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
  }

  for (const auto &input : inputs) {
    const auto referenceVar = reference.inputVarByKey.find(input.referenceKey);
    const auto implementationVar =
        implementation.inputVarByKey.find(input.implementationKey);
    if (referenceVar == reference.inputVarByKey.end() ||
        implementationVar == implementation.inputVarByKey.end()) {
      throw UnsupportedC2Rtl("C2RTL data input `" + input.name +
                             "` has no extracted Boolean variable");
    }
    const size_t symbol = nextSymbol++;
    referenceSymbols.emplace(referenceVar->second, symbol);
    implementationSymbols.emplace(implementationVar->second, symbol);
    problem.inputSymbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    problem.environmentInputNames.push_back(input.name);
  }

  std::optional<size_t> resetSymbol;
  if (controls.reset.has_value()) {
    const auto resetVar =
        implementation.inputVarByKey.find(controls.reset->key);
    if (resetVar == implementation.inputVarByKey.end()) {
      throw UnsupportedC2Rtl("C2RTL reset has no extracted Boolean variable");
    }
    resetSymbol = nextSymbol++;
    implementationSymbols.emplace(resetVar->second, *resetSymbol);
    problem.inputSymbols.push_back(*resetSymbol);
    problem.allSymbols.push_back(*resetSymbol);
    problem.environmentInputNames.push_back(controls.reset->name);
  }

  std::unordered_map<BoolExpr *, BoolExpr *> referenceMemo;
  std::unordered_map<BoolExpr *, BoolExpr *> implementationMemo;
  for (const auto &key : implementation.stateBits) {
    const auto nextIt = implementation.nextStateExprByStateKey.find(key);
    if (nextIt == implementation.nextStateExprByStateKey.end()) {
      throw UnsupportedC2Rtl(
          "C2RTL is missing an RTL next-state function for `" +
          displayNameForKey(implementation, key, "state") + "`");
    }
    problem.transitions1.emplace_back(
        rtlStateSymbols.at(key),
        remapFormulaStrict(nextIt->second, implementationSymbols,
                           implementationMemo, "RTL next-state"));
  }

  BoolExpr *initialCondition = BoolExpr::createTrue();
  for (const auto &[key, value] : implementation.initialStateValueByKey) {
    const auto symbolIt = rtlStateSymbols.find(key);
    if (symbolIt != rtlStateSymbols.end()) {
      addInitialAssignment(problem, symbolIt->second, value, initialCondition);
    }
  }

  struct ComparedBit {
    std::string name;
    size_t delay = 0;
    BoolExpr *referenceExpr = nullptr;
    BoolExpr *implementationExpr = nullptr;
    std::vector<size_t> history;
  };
  std::vector<ComparedBit> comparedBits;
  size_t maxDelay = 0;
  for (const auto &output : outputs) {
    maxDelay = std::max(maxDelay, output.delay);
    for (size_t bit = 0; bit < output.width; ++bit) {
      BoolExpr *referenceExpr = BoolExpr::createFalse();
      if (output.referenceBits[bit].has_value()) {
        referenceExpr = remapFormulaStrict(
            reference.observedOutputExprByKey.at(*output.referenceBits[bit]),
            referenceSymbols, referenceMemo, "reference output");
      }
      BoolExpr *implementationExpr = BoolExpr::createFalse();
      if (output.implementationBits[bit].has_value()) {
        implementationExpr = remapFormulaStrict(
            implementation.observedOutputExprByKey.at(
                *output.implementationBits[bit]),
            implementationSymbols, implementationMemo, "RTL output");
      }
      comparedBits.push_back({output.name + "[" + std::to_string(bit) + "]",
                              output.delay,
                              referenceExpr,
                              implementationExpr,
                              {}});
    }
  }

  for (auto &bit : comparedBits) {
    bit.history.reserve(bit.delay);
    for (size_t stage = 0; stage < bit.delay; ++stage) {
      const size_t symbol = nextSymbol++;
      bit.history.push_back(symbol);
      problem.auxiliaryStateSymbols.push_back(symbol);
      problem.allSymbols.push_back(symbol);
      addInitialAssignment(problem, symbol, false, initialCondition);
    }
  }

  std::vector<size_t> validSymbols;
  validSymbols.reserve(maxDelay);
  for (size_t stage = 0; stage < maxDelay; ++stage) {
    const size_t symbol = nextSymbol++;
    validSymbols.push_back(symbol);
    problem.auxiliaryStateSymbols.push_back(symbol);
    problem.allSymbols.push_back(symbol);
    addInitialAssignment(problem, symbol, false, initialCondition);
  }

  BoolExpr *resetActive = nullptr;
  BoolExpr *notResetActive = BoolExpr::createTrue();
  if (resetSymbol.has_value()) {
    resetActive = controls.reset->activeHigh
                      ? BoolExpr::Var(*resetSymbol)
                      : BoolExpr::Not(BoolExpr::Var(*resetSymbol));
    notResetActive = BoolExpr::Not(resetActive);
  }

  for (const auto &bit : comparedBits) {
    for (size_t stage = 0; stage < bit.history.size(); ++stage) {
      BoolExpr *source = stage == 0 ? bit.referenceExpr
                                    : BoolExpr::Var(bit.history[stage - 1]);
      if (resetActive != nullptr) {
        source = BoolExpr::And(notResetActive, source);
      }
      problem.auxiliaryTransitions.emplace_back(bit.history[stage], source);
    }
  }
  for (size_t stage = 0; stage < validSymbols.size(); ++stage) {
    BoolExpr *nextValid = stage == 0 ? BoolExpr::createTrue()
                                     : BoolExpr::Var(validSymbols[stage - 1]);
    if (resetActive != nullptr) {
      nextValid = BoolExpr::And(notResetActive, nextValid);
    }
    problem.auxiliaryTransitions.emplace_back(validSymbols[stage], nextValid);
  }

  BoolExpr *property = BoolExpr::createTrue();
  for (const auto &bit : comparedBits) {
    BoolExpr *delayedReference =
        bit.delay == 0 ? bit.referenceExpr : BoolExpr::Var(bit.history.back());
    BoolExpr *bitProperty =
        makeEquality(delayedReference, bit.implementationExpr);
    if (bit.delay != 0) {
      bitProperty = BoolExpr::Or(
          BoolExpr::Not(BoolExpr::Var(validSymbols[bit.delay - 1])),
          bitProperty);
    }
    if (resetActive != nullptr) {
      bitProperty = BoolExpr::Or(resetActive, bitProperty);
    }
    property = BoolExpr::And(property, bitProperty);
  }

  problem.initialCondition = BoolExpr::simplify(initialCondition);
  problem.totalStateCount =
      problem.state1Symbols.size() + problem.auxiliaryStateSymbols.size();
  problem.property = BoolExpr::simplify(property);
  problem.bad = BoolExpr::simplify(BoolExpr::Not(problem.property));
  problem.originalObservedOutputCount = comparedBits.size();
  problem.description = "C2RTL delayed combinational-reference safety property";
  built.comparedBits = comparedBits.size();
  return built;
}

C2RtlEquivalenceResult unsupportedResult(const std::string &reason) {
  C2RtlEquivalenceResult result;
  result.status = C2RtlEquivalenceStatus::Unsupported;
  result.reason = reason;
  return result;
}

} // namespace

C2RtlEquivalenceStrategy::C2RtlEquivalenceStrategy(
    naja::NL::SNLDesign *combinationalReference,
    naja::NL::SNLDesign *clockedImplementation,
    KEPLER_FORMAL::Config::SolverType solverType,
    C2RtlEquivalenceOptions options)
    : combinationalReference_(combinationalReference),
      clockedImplementation_(clockedImplementation), solverType_(solverType),
      options_(std::move(options)) {}

C2RtlEquivalenceResult C2RtlEquivalenceStrategy::run(size_t maxFrames) const {
  if (combinationalReference_ == nullptr || clockedImplementation_ == nullptr) {
    return unsupportedResult("C2RTL requires both parsed RTL designs");
  }

  try {
    const SequentialDesignModel reference =
        SequentialDesignModel::extract(combinationalReference_);
    const SequentialDesignModel implementation =
        SequentialDesignModel::extract(clockedImplementation_);
    ensureModelIsUsable(reference, "combinational reference");
    ensureModelIsUsable(implementation, "clocked RTL");
    if (!reference.stateBits.empty()) {
      throw UnsupportedC2Rtl(
          "C2RTL requires the XLS-generated reference RTL to be purely "
          "combinational");
    }
    ensureFeedForwardRtlState(implementation);

    const RtlControls controls =
        classifyRtlControls(clockedImplementation_, implementation);
    const auto inputs = alignDataInputs(reference, implementation, controls);
    const auto outputs = alignOutputs(reference, implementation, options_);
    const BuiltC2RtlProblem built =
        buildProblem(reference, implementation, controls, inputs, outputs);

    KEPLER_FORMAL::SEC::PDREngine engine(built.problem, solverType_);
    const PDRResult proof = engine.run(maxFrames);

    C2RtlEquivalenceResult result;
    result.bound = proof.bound;
    result.comparedBits = built.comparedBits;
    result.comparedOutputs = outputs.size();
    result.clock = controls.clockName + " " +
                   KEPLER_FORMAL::SEC::clockPhaseName(controls.clockPhase);
    if (controls.reset.has_value()) {
      result.reset = controls.reset->name;
      result.resetActiveHigh = controls.reset->activeHigh;
    }
    for (const auto &output : outputs) {
      result.outputDelays.emplace_back(output.name, output.delay);
    }
    switch (proof.status) {
    case PDRStatus::Equivalent:
      result.status = C2RtlEquivalenceStatus::Equivalent;
      result.reason = "PDR proved every configured delayed output relation";
      break;
    case PDRStatus::Different:
      result.status = C2RtlEquivalenceStatus::Different;
      result.reason =
          "PDR found a trace that violates a configured delayed output "
          "relation";
      break;
    case PDRStatus::Inconclusive:
    default:
      result.status = C2RtlEquivalenceStatus::Inconclusive;
      result.reason =
          "PDR did not prove or refute the C2RTL property within the "
          "configured frame bound";
      break;
    }
    return result;
  } catch (const UnsupportedC2Rtl &error) {
    return unsupportedResult(error.what());
  } catch (const std::exception &error) {
    C2RtlEquivalenceResult result;
    result.status = C2RtlEquivalenceStatus::Inconclusive;
    result.reason =
        std::string("C2RTL proof construction failed: ") + error.what();
    return result;
  }
}

} // namespace KEPLER_FORMAL::C2RTL
