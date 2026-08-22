// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "model/SecNetlistChecks.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "SNLDesignModeling.h"

namespace KEPLER_FORMAL::SEC {

namespace {

using DNLID = naja::DNL::DNLID;
using BooleanExpression = naja::NL::SNLDesignModeling::BooleanExpression;
using DependencyGraph = std::vector<std::vector<DNLID>>;

void addDependency(DependencyGraph& dependencies, DNLID source, DNLID output) {
  if (source == naja::DNL::DNLID_MAX || output == naja::DNL::DNLID_MAX ||
      source >= dependencies.size() || output >= dependencies.size()) {
    return;
  }
  dependencies[source].push_back(output);
}

struct ExpressionDependencies {
  std::vector<const naja::NL::SNLBitTerm*> terms;
  std::vector<size_t> states;
};

ExpressionDependencies
collectExpressionDependencies(const BooleanExpression& expression) {
  ExpressionDependencies dependencies;
  if (!expression.isValid()) {
    return dependencies;
  }

  std::vector<uint8_t> visited(expression.nodes.size(), 0);
  std::vector<BooleanExpression::NodeID> work{expression.root};
  std::unordered_set<const naja::NL::SNLBitTerm*> seenTerms;
  std::unordered_set<size_t> seenStates;
  while (!work.empty()) {
    const auto nodeID = work.back();
    work.pop_back();
    if (nodeID >= expression.nodes.size() || visited[nodeID] != 0) {
      continue;
    }
    visited[nodeID] = 1;
    const auto& node = expression.nodes[nodeID];
    if (node.operation == BooleanExpression::Operator::Term &&
        node.term != nullptr && seenTerms.insert(node.term).second) {
      dependencies.terms.push_back(node.term);
    } else if (node.operation == BooleanExpression::Operator::State &&
               seenStates.insert(node.state).second) {
      dependencies.states.push_back(node.state);
    }
    work.insert(work.end(), node.operands.begin(), node.operands.end());
  }
  return dependencies;
}

using InstanceTermMap = std::unordered_map<const naja::NL::SNLBitTerm*, DNLID>;

void addExpressionDependencies(
    const BooleanExpression& expression, const std::vector<DNLID>& targets,
    const InstanceTermMap& instanceTerms,
    const std::vector<std::vector<DNLID>>& stateSources,
    DependencyGraph& dependencies) {
  if (targets.empty()) {
    return;
  }
  const auto sources = collectExpressionDependencies(expression);
  for (const auto* sourceTerm : sources.terms) {
    const auto sourceIt = instanceTerms.find(sourceTerm);
    if (sourceIt == instanceTerms.end()) {
      continue;
    }
    for (const auto target : targets) {
      addDependency(dependencies, sourceIt->second, target);
    }
  }
  for (const auto sourceState : sources.states) {
    if (sourceState >= stateSources.size()) {
      continue;
    }
    for (const auto source : stateSources[sourceState]) {
      for (const auto target : targets) {
        addDependency(dependencies, source, target);
      }
    }
  }
}

std::optional<DNLID> findInstanceTerm(const InstanceTermMap& instanceTerms,
                                      const naja::NL::SNLBitTerm* term) {
  if (term == nullptr) {
    return std::nullopt;
  }
  const auto it = instanceTerms.find(term);
  return it == instanceTerms.end() ? std::nullopt
                                   : std::optional<DNLID>{it->second};
}

bool supportsStructuredMemoryModel(
    const naja::NL::SNLDesignModeling::MemoryInterface& interface) {
  if (!interface.isValid()) {
    return false;  // LCOV_EXCL_LINE - Naja rejects invalid interfaces on set.
  }
  for (const auto& readPort : interface.readPorts) {
    if (readPort.address.size() != interface.abits ||
        readPort.data.size() != interface.width) {
      return false;
    }
  }
  for (const auto& writePort : interface.writePorts) {
    if (writePort.address.size() != interface.abits ||
        writePort.data.size() != interface.width ||
        (!writePort.mask.empty() && writePort.mask.size() != interface.width) ||
        !writePort.extraWriteInputs.empty()) {
      return false;
    }
  }
  return true;
}

bool isDisabledMemoryWriteEnable(naja::DNL::DNLFull* dnl, DNLID termID) {
  if (termID == naja::DNL::DNLID_MAX) {
    return true;  // LCOV_EXCL_LINE - callers pass only mapped instance terms.
  }
  const auto& term = dnl->getDNLTerminalFromID(termID);
  if (term.isNull() || term.getIsoID() == naja::DNL::DNLID_MAX) {
    return true;  // LCOV_EXCL_LINE - mapped terms receive an iso from DNL.
  }
  const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
  return iso.isConstant0() || (!iso.isConstant() && iso.getDrivers().empty());
}

void addCombinationalDependencies(const naja::DNL::DNLInstanceFull& instance,
                                  const InstanceTermMap& instanceTerms,
                                  const std::unordered_set<DNLID>& opaqueTerms,
                                  DependencyGraph& dependencies) {
  for (const auto& [outputBitTerm, outputTermID] : instanceTerms) {
    if (outputBitTerm->getDirection() ==
            naja::NL::SNLBitTerm::Direction::Input ||
        opaqueTerms.find(outputTermID) != opaqueTerms.end()) {
      continue;
    }
    const auto truthTable = naja::NL::SNLDesignModeling::getTruthTable(
        instance.getSNLInstance(), outputBitTerm->getOrderID());
    if (!truthTable.isInitialized()) {
      continue;
    }
    for (auto* inputBitTerm :
         naja::NL::SNLDesignModeling::getCombinatorialInputs(
             const_cast<naja::NL::SNLBitTerm*>(outputBitTerm))) {
      const auto inputIt = instanceTerms.find(inputBitTerm);
      if (inputIt != instanceTerms.end()) {
        addDependency(dependencies, inputIt->second, outputTermID);
      }
    }
  }
}

void addSequentialDependencies(const naja::DNL::DNLInstanceFull& instance,
                               const InstanceTermMap& instanceTerms,
                               const std::unordered_set<DNLID>& opaqueTerms,
                               DependencyGraph& dependencies) {
  const auto* primitive = instance.getSNLModel();
  if (primitive == nullptr ||
      !naja::NL::SNLDesignModeling::hasSequentialModel(primitive)) {
    return;
  }
  const auto& model =
      naja::NL::SNLDesignModeling::getSequentialModel(primitive);
  if (model.kind !=
      naja::NL::SNLDesignModeling::SequentialModel::Kind::FlipFlop) {
    return;
  }

  std::vector<std::vector<DNLID>> physicalOutputsByState(model.states.size());
  std::vector<std::vector<DNLID>> modeledOutputsByState(model.states.size());
  std::vector<DNLID> allModeledOutputs;
  for (const auto& output : model.outputs) {
    const auto outputIt = instanceTerms.find(output.term);
    if (outputIt == instanceTerms.end()) {
      continue;
    }
    const auto outputDependencies =
        collectExpressionDependencies(output.function);
    for (const auto state : outputDependencies.states) {
      if (state >= model.states.size()) {
        continue;
      }
      physicalOutputsByState[state].push_back(outputIt->second);
      if (opaqueTerms.find(outputIt->second) == opaqueTerms.end()) {
        modeledOutputsByState[state].push_back(outputIt->second);
        allModeledOutputs.push_back(outputIt->second);
      }
    }
  }

  std::vector<std::vector<DNLID>> stateSources = modeledOutputsByState;
  for (size_t state = 0; state < stateSources.size(); ++state) {
    if (stateSources[state].empty()) {
      stateSources[state] = physicalOutputsByState[state];
    }
  }
  std::sort(allModeledOutputs.begin(), allModeledOutputs.end());
  allModeledOutputs.erase(
      std::unique(allModeledOutputs.begin(), allModeledOutputs.end()),
      allModeledOutputs.end());

  addExpressionDependencies(model.clockedOn, allModeledOutputs, instanceTerms,
                            stateSources, dependencies);
  for (size_t state = 0; state < model.states.size(); ++state) {
    auto& targets = modeledOutputsByState[state];
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
    const auto& sequentialState = model.states[state];
    addExpressionDependencies(sequentialState.nextState, targets, instanceTerms,
                              stateSources, dependencies);
    if (sequentialState.clear.has_value()) {
      addExpressionDependencies(*sequentialState.clear, targets, instanceTerms,
                                stateSources, dependencies);
    }
    if (sequentialState.preset.has_value()) {
      addExpressionDependencies(*sequentialState.preset, targets, instanceTerms,
                                stateSources, dependencies);
    }
  }
}

void addStructuredMemoryDependencies(
    naja::DNL::DNLFull* dnl, const naja::DNL::DNLInstanceFull& instance,
    const InstanceTermMap& instanceTerms,
    const std::unordered_set<DNLID>& opaqueTerms,
    DependencyGraph& dependencies) {
  const auto* primitive = instance.getSNLModel();
  if (primitive == nullptr ||
      !naja::NL::SNLDesignModeling::hasMemoryInterface(primitive)) {
    return;
  }
  const auto memory = naja::NL::SNLDesignModeling::getMemoryInterface(
      instance.getSNLInstance());
  if (!supportsStructuredMemoryModel(memory)) {
    return;
  }

  std::vector<DNLID> allReadOutputs;
  std::vector<std::vector<DNLID>> readOutputs(memory.readPorts.size());
  for (size_t port = 0; port < memory.readPorts.size(); ++port) {
    for (const auto* dataTerm : memory.readPorts[port].data) {
      const auto output = findInstanceTerm(instanceTerms, dataTerm);
      if (output.has_value() &&
          opaqueTerms.find(*output) == opaqueTerms.end()) {
        readOutputs[port].push_back(*output);
        allReadOutputs.push_back(*output);
      }
    }
  }
  std::sort(allReadOutputs.begin(), allReadOutputs.end());
  allReadOutputs.erase(
      std::unique(allReadOutputs.begin(), allReadOutputs.end()),
      allReadOutputs.end());

  auto addTermsToTargets = [&](const auto& terms,
                               const std::vector<DNLID>& targets) {
    for (const auto* term : terms) {
      const auto source = findInstanceTerm(instanceTerms, term);
      if (!source.has_value()) {
        continue;  // LCOV_EXCL_LINE - Naja validates memory-term ownership.
      }
      for (const auto target : targets) {
        addDependency(dependencies, *source, target);
      }
    }
  };
  for (size_t port = 0; port < memory.readPorts.size(); ++port) {
    addTermsToTargets(memory.readPorts[port].address, readOutputs[port]);
  }
  for (const auto& writePort : memory.writePorts) {
    const bool disabled =
        std::any_of(writePort.enables.begin(), writePort.enables.end(),
                    [&](const auto* enableTerm) {
                      const auto enable =
                          findInstanceTerm(instanceTerms, enableTerm);
                      return !enable.has_value() ||
                             isDisabledMemoryWriteEnable(dnl, *enable);
                    });
    if (disabled) {
      continue;
    }
    addTermsToTargets(writePort.address, allReadOutputs);
    addTermsToTargets(writePort.data, allReadOutputs);
    addTermsToTargets(writePort.mask, allReadOutputs);
    addTermsToTargets(writePort.enables, allReadOutputs);
  }
  if (memory.resetMode != naja::NL::SNLDesignModeling::MemoryResetMode::None) {
    const auto source = findInstanceTerm(instanceTerms, memory.reset);
    if (!source.has_value()) {
      return;
    }
    for (const auto target : allReadOutputs) {
      addDependency(dependencies, *source, target);
    }
  }
}

DependencyGraph buildLocalDependencyGraph(
    naja::DNL::DNLFull* dnl, const std::unordered_set<DNLID>& opaqueTerms) {
  DependencyGraph dependencies(dnl->getNBterms());
  for (const auto leafID : dnl->getLeaves()) {
    const auto& instance = dnl->getDNLInstanceFromID(leafID);
    InstanceTermMap instanceTerms;
    for (DNLID termID = instance.getTermIndexes().first;
         termID != naja::DNL::DNLID_MAX &&
         termID <= instance.getTermIndexes().second;
         ++termID) {
      const auto& term = dnl->getDNLTerminalFromID(termID);
      if (!term.isNull()) {
        instanceTerms.emplace(term.getSnlBitTerm(), termID);
      }
    }
    addCombinationalDependencies(instance, instanceTerms, opaqueTerms,
                                 dependencies);
    addSequentialDependencies(instance, instanceTerms, opaqueTerms,
                              dependencies);
    addStructuredMemoryDependencies(dnl, instance, instanceTerms, opaqueTerms,
                                    dependencies);
  }
  for (auto& outputs : dependencies) {
    std::sort(outputs.begin(), outputs.end());
    outputs.erase(std::unique(outputs.begin(), outputs.end()), outputs.end());
  }
  return dependencies;
}

std::vector<DNLID> findReachedTopOutputs(naja::DNL::DNLFull* dnl,
                                         const DependencyGraph& dependencies,
                                         DNLID seed,
                                         std::vector<uint8_t>& visited,
                                         std::vector<DNLID>& work) {
  std::vector<DNLID> reached;
  auto enqueue = [&](DNLID termID) {
    if (termID < visited.size() && visited[termID] == 0) {
      visited[termID] = 1;
      work.push_back(termID);
    }
  };
  enqueue(seed);
  for (size_t next = 0; next < work.size(); ++next) {
    const DNLID termID = work[next];
    const auto& term = dnl->getDNLTerminalFromID(termID);
    if (term.isNull()) {
      continue;
    }
    if (term.isTopPort() && term.getSnlBitTerm()->getDirection() !=
                                naja::NL::SNLBitTerm::Direction::Input) {
      reached.push_back(termID);
      continue;
    }
    for (const auto output : dependencies[termID]) {
      enqueue(output);
    }
    if (term.getSnlBitTerm()->getDirection() ==
            naja::NL::SNLBitTerm::Direction::Input ||
        term.getIsoID() == naja::DNL::DNLID_MAX) {
      continue;
    }
    const auto& iso = dnl->getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
    for (const auto reader : iso.getReaders()) {
      enqueue(reader);
    }
  }
  std::sort(reached.begin(), reached.end());
  reached.erase(std::unique(reached.begin(), reached.end()), reached.end());
  return reached;
}

}  // namespace

std::vector<OpaqueReachedTopOutput>
SecNetlistChecks::findTopOutputsReachedByOpaqueTerminals(
    std::vector<OpaqueTerminalSeed> opaqueSeeds) const {
  if (dnl_ == nullptr || opaqueSeeds.empty()) {
    return {};
  }
  std::sort(
      opaqueSeeds.begin(), opaqueSeeds.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.termID < rhs.termID; });
  opaqueSeeds.erase(std::unique(opaqueSeeds.begin(), opaqueSeeds.end(),
                                [](const auto& lhs, const auto& rhs) {
                                  return lhs.termID == rhs.termID;
                                }),
                    opaqueSeeds.end());

  std::unordered_set<DNLID> opaqueTerms;
  opaqueTerms.reserve(opaqueSeeds.size());
  for (const auto& seed : opaqueSeeds) {
    opaqueTerms.insert(seed.termID);
  }
  const auto dependencies = buildLocalDependencyGraph(dnl_, opaqueTerms);

  std::vector<std::vector<DNLID>> reachedBySeed(opaqueSeeds.size());
  tbb::parallel_for(tbb::blocked_range<size_t>(0, opaqueSeeds.size()),
                    [&](const tbb::blocked_range<size_t>& range) {
                      std::vector<uint8_t> visited(dnl_->getNBterms(), 0);
                      std::vector<DNLID> work;
                      for (size_t seedIndex = range.begin();
                           seedIndex != range.end(); ++seedIndex) {
                        reachedBySeed[seedIndex] = findReachedTopOutputs(
                            dnl_, dependencies, opaqueSeeds[seedIndex].termID,
                            visited, work);
                        for (const auto termID : work) {
                          visited[termID] = 0;
                        }
                        work.clear();
                      }
                    });

  std::map<DNLID, OpaqueTerminalSeed> sourceByTopOutput;
  for (size_t seedIndex = 0; seedIndex < opaqueSeeds.size(); ++seedIndex) {
    for (const auto topOutput : reachedBySeed[seedIndex]) {
      sourceByTopOutput.try_emplace(topOutput, opaqueSeeds[seedIndex]);
    }
  }
  std::vector<OpaqueReachedTopOutput> result;
  result.reserve(sourceByTopOutput.size());
  for (auto& [topOutput, source] : sourceByTopOutput) {
    result.push_back({topOutput, std::move(source)});
  }
  return result;
}

}  // namespace KEPLER_FORMAL::SEC
