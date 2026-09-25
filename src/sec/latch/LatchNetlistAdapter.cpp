// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchNetlistAdapter.h"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include "DNL.h"
#include "NLDB.h"
#include "NLName.h"
#include "NLUniverse.h"
#include "SNLDesign.h"
#include "SNLDesignModeling.h"
#include "SNLInstance.h"
#include "SNLPath.h"
#include "latch/LatchBoundaryEncoding.h"
#include "latch/LatchConstantNet.h"
#include "latch/LatchDependencyGraph.h"
#include "latch/LatchInputHistory.h"
#include "latch/LatchResetAdapter.h"
#include "latch/LatchResetClock.h"
#include "latch/LatchSupportOptions.h"
#include "latch/NajaEventPrimitive.h"
#include "latch/LatchSymbolicEncoding.h"
#include "model/OpaquePolicy.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {
using Term = naja::DNL::DNLTerminalFull;
using Direction = naja::NL::SNLBitTerm::Direction;
constexpr size_t absent = size_t(-1);

SignalKey key(const Term& term) {
  SignalKey result;
  for (const auto& name : term.getDNLInstance().getPath().getPathNames())
    result.first.push_back(name.getID());
  result.first.push_back(term.getSnlBitTerm()->getName().getID());
  result.second.push_back(term.getSnlBitTerm()->getBit());
  return result;
}
SignalKey syntheticKey(size_t category, size_t object, size_t bit = 0) {
  return {{uint64_t(1) << 61, category, object},
          {static_cast<naja::NL::NLID::DesignObjectID>(bit)}};
}
std::string name(const Term& term) {
  std::string result;
  for (const auto& part : term.getDNLInstance().getPath().getPathNames())
    result += part.getString() + ".";
  return result + term.getSnlBitTerm()->getName().getString() + "[" +
      std::to_string(term.getSnlBitTerm()->getBit()) + "]";
}
struct Port { SignalKey key; std::string name; size_t net; };
struct CellInfo { SignalKey key; std::string name, error; };

// Keep the caller's flattened graph and selected top intact (including borrowed
// designs). The compiled result owns no pointers into this temporary graph.
struct DnlScope {
  naja::NL::NLUniverse* universe = naja::NL::NLUniverse::get();
  naja::NL::NLDB* previousDb = universe->getTopDB();
  naja::NL::NLDB* selectedDb;
  naja::NL::SNLDesign* previousSelectedTop;
  naja::DNL::DNLFull* previousDnl = nullptr;
  std::vector<std::pair<naja::NL::SNLBitTerm*, naja::NL::NLID::DesignObjectID>> termOrders;
  std::vector<std::pair<naja::NL::SNLInstance*, naja::NL::NLID::DesignObjectID>> instanceOrders;
  explicit DnlScope(naja::NL::SNLDesign* top)
      : selectedDb(top->getDB()), previousSelectedTop(selectedDb->getTopDesign()) {
    // Flattening assigns ordering IDs on shared source models. Restoring only
    // the graph pointer would leave a caller's cached graph with stale IDs.
    std::set<naja::NL::SNLDesign*> visited;
    std::vector<naja::NL::SNLDesign*> pending{top};
    while (!pending.empty()) {
      auto* design = pending.back();
      pending.pop_back();
      if (!visited.insert(design).second) continue;
      for (auto* term : design->getBitTerms()) termOrders.emplace_back(term, term->getOrderID());
      for (auto* instance : design->getInstances()) {
        instanceOrders.emplace_back(instance, instance->getOrderID());
        pending.push_back(instance->getModel());
      }
    }
    previousDnl = naja::DNL::exchange(nullptr);
    universe->setTopDesign(top);
  }
  ~DnlScope() {
    naja::DNL::destroy();
    for (const auto& [term, order] : termOrders) term->setOrderID(order);
    for (const auto& [instance, order] : instanceOrders) instance->setOrderID(order);
    selectedDb->setTopDesign(previousSelectedTop);
    universe->setTopDB(previousDb);
    naja::DNL::exchange(previousDnl);
  }
};

void opaque(SequentialDesignModel& model, const SignalKey& key, const std::string& name,
            const std::string& reason) {
  model.displayNameByKey.insert_or_assign(key, name);
  model.connectivitySkipInfoByKey.insert_or_assign(key,
      ConnectivitySkipInfo{ConnectivitySkipOrigin::OpaqueInternal, reason});
}
BoolExpr* variable(SequentialDesignModel& model, const SignalKey& key,
                   const std::string& name, bool state, size_t& nextVar) {
  model.displayNameByKey.emplace(key, name);
  model.inputVarByKey.emplace(key, nextVar);
  (state ? model.stateBits : model.environmentInputs).push_back(key);
  return BoolExpr::Var(nextVar++);
}

SequentialDesignModel extract(naja::NL::SNLDesign* top, size_t side) {
  const auto options = supportOptions();
  SequentialDesignModel model;
  model.eventContract = std::string("boolean-epochs-v1;") +
      (options.singleInputChange ? "single;" : "any;") +
      "initial_inputs=" + std::to_string(*options.initialInputs) +
      ";initial_storage=" + std::to_string(*options.initialStorage);
  DnlScope scope(top);
  const auto* dnl = naja::DNL::get();
  DriverlessConstantResolver driverlessConstants(*dnl);
  Network network;
  std::map<naja::DNL::DNLID, size_t> byIso;
  std::vector<std::string> netErrors;
  auto net = [&](const Term& term) {
    const auto isoID = term.getIsoID();
    if (isoID != naja::DNL::DNLID_MAX) {
      if (const auto found = byIso.find(isoID); found != byIso.end()) return found->second;
    }
    const size_t index = network.netCount++;
    network.constantByNet.push_back({});
    netErrors.emplace_back();
    if (isoID == naja::DNL::DNLID_MAX) {
      const auto constant = driverlessConstants.resolve(term);
      network.constantByNet.back() = constant.value;
      if (constant.conflictingOrUnknown)
        netErrors.back() = "X/Z or conflicting constant net in Boolean event model: " + name(term);
      else if (!constant.value)
        netErrors.back() = "unconnected signal " + name(term);
      return index;
    }
    byIso.emplace(isoID, index);
    const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(isoID);
    if (iso.isConstant0() || iso.isConstant1()) network.constantByNet.back() = iso.isConstant1();
    else if (iso.isConstantX() || iso.isConstantZ() || iso.getType() == naja::DNL::DNLIso::AMBIGUOUS)
      netErrors.back() = "X/Z or conflicting constant net in Boolean event model: " + name(term);
    else if (iso.getDrivers().size() != 1) netErrors.back() = "missing or multiple signal drivers: " + name(term);
    return index;
  };
  std::vector<Port> inputs, outputs;
  const auto& topInstance = dnl->getTop();
  for (auto* bit : top->getBitTerms()) {
    const auto& term = topInstance.getTerminalFromBitTerm(bit);
    Port port{key(term), name(term), net(term)};
    model.displayNameByKey.emplace(port.key, port.name);
    if (bit->getDirection() == Direction::Input) {
      inputs.push_back(port); model.topInputKeys.push_back(port.key);
    } else if (bit->getDirection() == Direction::Output) {
      outputs.push_back(port); model.topOutputKeys.push_back(port.key);
      model.allObservedOutputs.push_back(port.key);
    } else model.unsupportedReasons.push_back("Boolean event model requires unidirectional top ports");
  }
  std::sort(inputs.begin(), inputs.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
  std::set<size_t> external;
  for (const auto& input : inputs) {
    if (!external.insert(input.net).second || network.constantByNet[input.net])
      model.unsupportedReasons.push_back("aliased/constant external input in event contract: " + input.name);
    network.externalInputs.push_back(input.net);
  }
  std::vector<CellInfo> cells;
  std::vector<SymbolicPrimitive> symbolicPrimitives;
  std::vector<ResetClockPrimitive> resetClocks;
  for (auto leaf : dnl->getLeaves()) {
    const auto& instance = dnl->getDNLInstanceFromID(leaf);
    if (instance.isTop()) continue;
    std::map<const naja::NL::SNLBitTerm*, size_t> pins;
    Primitive fallback;
    fallback.name = instance.getFullPath();
    CellInfo info{syntheticKey(3, cells.size()), fallback.name, {}};
    for (auto* bit : instance.getSNLModel()->getBitTerms()) {
      const auto& term = instance.getTerminalFromBitTerm(bit);
      const size_t index = net(term);
      pins.emplace(bit, index);
      if (bit->getDirection() == Direction::Input) fallback.inputs.push_back(index);
      else if (bit->getDirection() == Direction::Output) {
        fallback.outputs.push_back(index); info.key = key(term); info.name = name(term);
      } else info.error = "bidirectional primitive pin: " + name(term);
    }
    SymbolicPrimitive symbolic;
    ResetClockPrimitive resetClock;
    try {
      auto primitive = makeNajaEventPrimitive(instance.getSNLInstance(), fallback.name, pins, &symbolic, &resetClock);
      // Constant equipotentials already have an authoritative source. A cell
      // tied into one is not silently treated as a second varying writer.
      for (auto output : primitive.outputs)
        if (network.constantByNet[output]) throw std::runtime_error("primitive output connected to a constant net");
      network.primitives.push_back(std::move(primitive));
    } catch (const std::exception& error) {
      info.error = fallback.name + ": " + error.what();
      resetClock = {};
      network.primitives.push_back(std::move(fallback));
    }
    symbolicPrimitives.push_back(std::move(symbolic));
    resetClocks.push_back(std::move(resetClock));
    cells.push_back(std::move(info));
  }
  // A top wire observes a stored external level in selector mode. An observer
  // buffer has no consumers, so this extra observation wave cannot affect any
  // storage or capture event in the source circuit.
  for (auto& output : outputs) {
    if (!external.count(output.net) || network.constantByNet[output.net]) continue;
    const size_t observed = network.netCount++;
    network.constantByNet.push_back({}); netErrors.emplace_back();
    network.primitives.push_back(combinational("<top observation>", {output.net}, {observed},
        [](const Bits& value) { return value; }));
    SymbolicPrimitive observer;
    observer.react = [](const SymbolicBits&, const SymbolicBits&, const SymbolicBits& value,
                        std::optional<size_t>, bool) { return SymbolicReaction{{}, value}; };
    symbolicPrimitives.push_back(std::move(observer));
    ResetClockPrimitive observationClock;
    observationClock.kind = ResetClockPrimitive::Kind::Combinational;
    observationClock.outputs = {BoolExpr::Var(2)};
    resetClocks.push_back(std::move(observationClock));
    cells.push_back({syntheticKey(3, cells.size()), output.name, {}});
    output.net = observed;
  }
  if (!model.unsupportedReasons.empty()) return model;

  // Independent components are closed under ALL primitive input/output arcs,
  // including clock/enable/asynchronous controls. Shared read-only PIs are not
  // edges. This deliberately over-groups rather than dropping cross-island
  // pulses; waves inside a component still execute in parallel.
  const auto graph = analyzeDependencies(network);
  std::vector<size_t> owner(network.netCount, absent);
  for (size_t i = 0; i < network.primitives.size(); ++i) {
    for (auto output : network.primitives[i].outputs) {
      if (owner[output] != absent) {
        netErrors[output] = "multiple primitive writers";
      } else owner[output] = i;
    }
  }

  size_t nextVar = 2;
  auto resetInterface = std::make_shared<EventResetInterface>();
  resetInterface->singleInputChange = options.singleInputChange;
  resetInterface->maxCompositionNodes = options.maxSymbolicNodes;
  for (const auto& input : inputs) {
    resetInterface->inputKeys.push_back(input.key);
    resetInterface->inputNames.push_back(input.name);
  }
  resetInterface->currentInputs.assign(inputs.size(), BoolExpr::Var(*options.initialInputs ? 1 : 0));
  std::vector<BoolExpr*> inputExpressions, selector;
  BoolExpr* eventValue = nullptr;
  if (!options.singleInputChange) {
    for (const auto& input : inputs)
      inputExpressions.push_back(variable(model, input.key, input.name, false, nextVar));
  } else {
    // SEC aligns keys, not labels. Retain the exact original input keys as
    // unused interface sentinels so selector positions cannot silently denote
    // different pins on the two sides. Only selector/value drive transactions.
    for (const auto& input : inputs)
      variable(model, input.key, "$event.interface." + input.name, false, nextVar);
    const size_t bits = boundaryEncodingBits(inputs.size() + 1);
    for (size_t i = 0; i < bits; ++i)
      selector.push_back(variable(model, syntheticKey(1, i + 1), "$event.select[" + std::to_string(i) + "]", false, nextVar));
    eventValue = variable(model, syntheticKey(1, bits + 1), "$event.value", false, nextVar);
    for (auto* bit : selector) resetInterface->selectorSymbols.push_back(bit->getId());
    resetInterface->valueSymbol = eventValue->getId();
    inputExpressions.assign(inputs.size(), BoolExpr::createFalse());
  }

  for (size_t componentID = 0; componentID < graph.components.size(); ++componentID) {
    const auto& members = graph.components[componentID];
    std::set<size_t> used;
    std::string failure;
    for (auto member : members) {
      const auto& primitive = network.primitives[member];
      used.insert(primitive.inputs.begin(), primitive.inputs.end());
      used.insert(primitive.outputs.begin(), primitive.outputs.end());
      if (failure.empty()) failure = cells[member].error;
    }
    for (auto index : used)
      if (failure.empty() && !netErrors[index].empty()) failure = netErrors[index];
    Network local;
    std::map<size_t, size_t> localNet;
    for (auto index : used) {
      localNet.emplace(index, local.netCount++);
      local.constantByNet.push_back(network.constantByNet[index]);
    }
    std::vector<size_t> globalInputIndices;
    std::vector<BoolExpr*> localInputs;
    for (size_t i = 0; i < inputs.size(); ++i) {
      if (!used.count(inputs[i].net)) continue;
      globalInputIndices.push_back(i);
      local.externalInputs.push_back(localNet.at(inputs[i].net));
      localInputs.push_back(inputExpressions[i]);
    }
    std::vector<SymbolicPrimitive> localSymbolic;
    for (auto member : members) {
      auto primitive = network.primitives[member];
      for (auto& index : primitive.inputs) index = localNet.at(index);
      for (auto& index : primitive.outputs) index = localNet.at(index);
      local.primitives.push_back(std::move(primitive));
      localSymbolic.push_back(symbolicPrimitives[member]);
    }
    std::optional<SymbolicMacro> symbolic;
    std::string symbolicFailure;
    if (failure.empty()) {
      SymbolicCompileOptions compile;
      compile.initialInputs.assign(local.externalInputs.size(), *options.initialInputs);
      compile.singleExternalInputChange = options.singleInputChange;
      compile.maxWaves = options.limits.maxWaves;
      compile.maxNodes = options.maxSymbolicNodes;
      compile.maxSatConflicts = options.maxSatConflicts;
      compile.maxSatDecisions = options.maxSatDecisions;
      compile.workers = options.workers;
      for (const auto& primitive : local.primitives)
        compile.initialStorage.emplace_back(primitive.storageBits, uint8_t(*options.initialStorage));
      auto result = compileSymbolicNetwork({local, std::move(localSymbolic)}, compile);
      if (result.certified()) symbolic = std::move(result.model);
      else symbolicFailure = result.detail;
    }
    std::optional<TransitionTable> table;
    if (failure.empty() && !symbolic) {
      try {
        EventModel reference(local, {}, options.workers);
        CompileOptions compile;
        compile.initialInputs.assign(local.externalInputs.size(), *options.initialInputs);
        compile.singleExternalInputChange = options.singleInputChange;
        compile.limits = options.limits;
        for (const auto& primitive : local.primitives)
          compile.initialStorage.emplace_back(primitive.storageBits, uint8_t(*options.initialStorage));
        auto result = compileTransitionTable(reference, compile);
        if (result.certified()) table = std::move(result.table);
        else failure = std::string(certificationStatusName(result.status)) + ": " + result.detail +
            "; symbolic certificate: " + symbolicFailure;
      } catch (const std::exception& error) { failure = error.what(); }
    }
    if (table && (table->initials.size() != 1 || table->initials[0].boundary >= table->boundaries.size()))
      failure = "event adapter requires a single explicitly initialized boundary";
    if (!failure.empty()) {
      const std::string reason = "event component " + cells[members.front()].name + ": " + failure;
      for (auto member : members) opaque(model, cells[member].key, cells[member].name, reason);
      for (const auto& output : outputs)
        if (owner[output.net] != absent && graph.componentOf[owner[output.net]] == componentID) {
          model.skippedObservedOutputs.push_back(output.key);
          opaque(model, output.key, output.name, reason);
        }
      continue;
    }
    std::vector<BoolExpr*> state;
    std::vector<SignalKey> stateKeys;
    const size_t stateWidth = symbolic ? symbolic->stateSymbols.size() : boundaryEncodingBits(table->boundaries.size());
    for (size_t bit = 0; bit < stateWidth; ++bit) {
      const auto stateKey = syntheticKey(2, componentID, bit);
      stateKeys.push_back(stateKey);
      state.push_back(variable(model, stateKey, "$event.component[" + std::to_string(componentID) + "].state[" + std::to_string(bit) + "]", true, nextVar));
      model.initialStateValueByKey.emplace(stateKey, symbolic ? symbolic->initialState.at(bit) : (table->initials[0].boundary >> bit) & 1);
    }
    const auto encoded = symbolic
        ? encodeSymbolicMacro(*symbolic, local, state, localInputs, options.singleInputChange, selector, eventValue, globalInputIndices)
        : encodeBoundaryTable(*table, local, state, localInputs, selector, eventValue, globalInputIndices);
    for (size_t bit = 0; bit < state.size(); ++bit)
      model.nextStateExprByStateKey.emplace(stateKeys[bit], encoded.nextState[bit]);
    for (size_t i = 0; i < globalInputIndices.size(); ++i)
      resetInterface->currentInputs[globalInputIndices[i]] = symbolic ? state.at(local.externalInputs[i])
          : finiteInputHistory(*table, state, local.externalInputs[i]);
    for (const auto& output : outputs) {
      if (owner[output.net] == absent || graph.componentOf[owner[output.net]] != componentID) continue;
      model.observedOutputs.push_back(output.key);
      model.observedOutputExprByKey.emplace(output.key, encoded.observedNets.at(localNet.at(output.net)));
    }
  }
  for (const auto& output : outputs) {
    if (owner[output.net] != absent) continue;
    if (network.constantByNet[output.net]) {
      model.observedOutputs.push_back(output.key);
      model.observedOutputExprByKey.emplace(output.key,
          *network.constantByNet[output.net] ? BoolExpr::createTrue() : BoolExpr::createFalse());
    } else {
      model.skippedObservedOutputs.push_back(output.key);
      opaque(model, output.key, output.name, netErrors[output.net].empty() ? "event output has no driver" : netErrors[output.net]);
    }
  }
  recordOpaqueOutputlessCells(model, *dnl);
  // Discovery failure limits automatic reset-cycle expansion only, not ordinary
  // event semantics. No clock is guessed from a latch enable or a signal name.
  const auto clock = discoverResetClock(network, resetClocks);
  resetInterface->clockInputIndex = clock.externalInputIndex;
  if (!clock.resolved()) resetInterface->clockError = clock.detail;
  model.eventResetInterface = std::move(resetInterface);
  applyOpaquePolicy(model, top->getName().getString(), side);
  return model;
}
}  // namespace

std::optional<SequentialDesignModel> extractEventDesign(
    naja::NL::SNLDesign* top, const BoundaryPairs& pairs, size_t side) {
  const auto& options = supportOptions();
  if (!options.enabled) return {};
  if (!options.initialInputs && !options.initialStorage) return {};
  if (!options.initialInputs || !options.initialStorage || !pairs.empty()) {
    SequentialDesignModel result;
    result.unsupportedReasons.push_back(!pairs.empty()
        ? "event semantics do not yet support selected leaf boundaries"
        : "event semantics require explicit Boolean initial_inputs and initial_storage");
    return result;
  }
  return extract(top, side);
}
}  // namespace KEPLER_FORMAL::SEC::LATCH
