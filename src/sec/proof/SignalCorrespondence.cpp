// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "proof/SignalCorrespondence.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <iterator>
#include <map>
#include <random>
#include <unordered_map>

#include "kinduction/SatEncoding.h"
#include "proof/TransitionExprResolver.h"

namespace KEPLER_FORMAL::SEC {

namespace {

// Simulation restarts from the start state this many times. A run goes on
// until the partition has stopped refining over two looks in a row, each
// look this many frames after the last, or until a limit is reached.
constexpr size_t kSimulationRuns = 4;
constexpr size_t kSimulationLook = 64;
constexpr size_t kSimulationFrameLimit = size_t{1} << 13;
// Node evaluations over all runs.
constexpr uint64_t kSimulationWorkLimit = uint64_t{1} << 31;
// A query has no tick limit of its own, see SignalCorrespondenceLimits.
constexpr int64_t kNoTickLimit = -1;

constexpr uint8_t kUnknown = 2;
constexpr uint32_t kNone = ~uint32_t{0};
constexpr uint64_t kAllPatterns = ~uint64_t{0};

enum class Gate : uint8_t { False, True, Leaf, Not, And, Or, Xor };

// The logic of both designs as one flat list in topological order, so that a
// time frame is evaluated by a single pass over an array, in any value domain.
class Circuit {
 public:
  explicit Circuit(size_t nodeLimit) : nodeLimit_(nodeLimit) {}

  struct Node {
    Gate gate = Gate::False;
    // Operand nodes; for a leaf, `left` is the slot of its symbol.
    uint32_t left = 0;
    uint32_t right = 0;
  };

  // Adds the cone of `root`. False when the cone holds an operator this pass
  // does not read or the node limit is reached.
  bool add(BoolExpr* root, uint32_t& index) {
    stack_.assign({root});
    while (!stack_.empty()) {
      BoolExpr* expr = stack_.back();
      if (index_.contains(expr)) {
        stack_.pop_back();
        continue;
      }
      Node node;
      const Op op = expr->getOp();
      if (op == Op::VAR) {
        const size_t symbol = expr->getId();
        if (symbol == static_cast<size_t>(-1)) {
          return false;
        }
        node.gate = symbol == 0 ? Gate::False : symbol == 1 ? Gate::True : Gate::Leaf;
        if (node.gate == Gate::Leaf) {
          node.left = slotOf(symbol);
        }
      } else {
        BoolExpr* left = expr->getLeft();
        BoolExpr* right = op == Op::NOT ? nullptr : expr->getRight();
        const bool readable = op == Op::NOT || op == Op::AND || op == Op::OR || op == Op::XOR;
        if (!readable || left == nullptr || (op != Op::NOT && right == nullptr)) {
          return false;
        }
        const auto leftIt = index_.find(left);
        const auto rightIt = right == nullptr ? index_.end() : index_.find(right);
        const bool rightReady = right == nullptr || rightIt != index_.end();
        if (leftIt == index_.end() || !rightReady) {
          if (leftIt == index_.end()) stack_.push_back(left);
          if (!rightReady) stack_.push_back(right);
          continue;
        }
        node.gate = op == Op::NOT ? Gate::Not : op == Op::AND ? Gate::And
                  : op == Op::OR ? Gate::Or : Gate::Xor;
        node.left = leftIt->second;
        node.right = right == nullptr ? 0 : rightIt->second;
      }
      if (nodes_.size() >= nodeLimit_) {
        return false;
      }
      index_.emplace(expr, static_cast<uint32_t>(nodes_.size()));
      nodes_.push_back(node);
      exprs_.push_back(expr);
      stack_.pop_back();
    }
    index = index_.at(root);
    return true;
  }

  size_t size() const { return nodes_.size(); }
  const Node& node(uint32_t index) const { return nodes_[index]; }
  BoolExpr* expr(uint32_t index) const { return exprs_[index]; }
  // Leaf symbols, by slot.
  const std::vector<size_t>& symbols() const { return symbols_; }

  uint32_t slotOf(size_t symbol) {
    const auto [it, inserted] = slots_.emplace(symbol, static_cast<uint32_t>(symbols_.size()));
    if (inserted) {
      symbols_.push_back(symbol);
    }
    return it->second;
  }

 private:
  size_t nodeLimit_;
  std::vector<Node> nodes_;
  std::vector<BoolExpr*> exprs_;
  std::unordered_map<BoolExpr*, uint32_t> index_;
  std::vector<size_t> symbols_;
  std::unordered_map<size_t, uint32_t> slots_;
  std::vector<BoolExpr*> stack_;
};

// One bit per simulation pattern.
void evaluate(const Circuit& circuit, const std::vector<uint64_t>& leaves,
              std::vector<uint64_t>& values) {
  values.resize(circuit.size());
  for (uint32_t i = 0; i < circuit.size(); ++i) {
    const auto& node = circuit.node(i);
    switch (node.gate) {
      case Gate::False: values[i] = 0; break;
      case Gate::True: values[i] = kAllPatterns; break;
      case Gate::Leaf: values[i] = leaves[node.left]; break;
      case Gate::Not: values[i] = ~values[node.left]; break;
      case Gate::And: values[i] = values[node.left] & values[node.right]; break;
      case Gate::Or: values[i] = values[node.left] | values[node.right]; break;
      case Gate::Xor: values[i] = values[node.left] ^ values[node.right]; break;
    }
  }
}

// 0, 1 or unknown: what a signal is known to be when some leaves are not.
void evaluateTernary(const Circuit& circuit, const std::vector<uint8_t>& leaves,
                     std::vector<uint8_t>& values) {
  values.resize(circuit.size());
  for (uint32_t i = 0; i < circuit.size(); ++i) {
    const auto& node = circuit.node(i);
    const uint8_t a = node.gate == Gate::Leaf || node.gate == Gate::False ||
                              node.gate == Gate::True ? 0 : values[node.left];
    const uint8_t b = node.gate == Gate::And || node.gate == Gate::Or ||
                              node.gate == Gate::Xor ? values[node.right] : 0;
    switch (node.gate) {
      case Gate::False: values[i] = 0; break;
      case Gate::True: values[i] = 1; break;
      case Gate::Leaf: values[i] = leaves[node.left]; break;
      case Gate::Not: values[i] = a == kUnknown ? kUnknown : 1 - a; break;
      case Gate::And:
        values[i] = a == 0 || b == 0 ? 0 : a == 1 && b == 1 ? 1 : kUnknown;
        break;
      case Gate::Or:
        values[i] = a == 1 || b == 1 ? 1 : a == 0 && b == 0 ? 0 : kUnknown;
        break;
      case Gate::Xor:
        values[i] = a == kUnknown || b == kUnknown ? kUnknown : a ^ b;
        break;
    }
  }
}

struct Register {
  size_t symbol = 0;
  // Leaf slot of the symbol.
  uint32_t slot = 0;
  // The register read as a signal.
  uint32_t node = 0;
  // Its next-state function.
  uint32_t next = 0;
};

// The product machine of van Eijk, section III.
struct Machine {
  explicit Machine(size_t nodeLimit) : circuit(nodeLimit) {}

  Circuit circuit;
  std::vector<Register> registers;
  std::unordered_map<size_t, uint32_t> registerBySymbol;
  // The register a leaf slot reads, kNone for an input.
  std::vector<uint32_t> registerOfSlot;
  // Per leaf slot: the value a reset port holds after the reset window,
  // kUnknown for every other leaf.
  std::vector<uint8_t> idleReset;
  std::vector<std::pair<size_t, bool>> resetPorts;
  size_t resetCycles = 0;
  uint32_t constantTrue = 0;
  // Per observed output, the node pairs (design 1, design 2) that must be
  // equal; empty for an output that is already proved.
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> outputs;
  // A machine of the binary functions starts at boot: per register its value
  // there (0, 1 or unknown) and the register it shares that value with.
  bool binary = false;
  std::vector<uint8_t> boot;
  std::vector<uint32_t> commonRoot;
};

// The leaf slot of every register, and the value the reset ports rest at.
void indexLeaves(const KInductionProblem& problem, Machine& machine) {
  const Circuit& circuit = machine.circuit;
  machine.registerOfSlot.assign(circuit.symbols().size(), kNone);
  for (uint32_t i = 0; i < machine.registers.size(); ++i) {
    machine.registerOfSlot[machine.registers[i].slot] = i;
  }
  machine.idleReset.assign(circuit.symbols().size(), kUnknown);
  if (problem.hasResetBootstrap() && problem.resetBootstrapCycles != 0) {
    machine.resetCycles = problem.resetBootstrapCycles;
    machine.resetPorts = problem.resetBootstrapInputs;
    for (uint32_t slot = 0; slot < circuit.symbols().size(); ++slot) {
      for (const auto& [symbol, active] : machine.resetPorts) {
        if (circuit.symbols()[slot] == symbol) {
          machine.idleReset[slot] = active ? 0 : 1;
        }
      }
    }
  }
}

bool buildMachine(const KInductionProblem& problem, const SignalCorrespondenceLimits& limits,
                  Machine& machine) {
  const TransitionExprResolver transitions(problem);
  if (problem.state0Symbols.size() + problem.state1Symbols.size() +
          problem.auxiliaryStateSymbols.size() > limits.stateSymbols) {
    return false;
  }
  Circuit& circuit = machine.circuit;
  if (!circuit.add(BoolExpr::createTrue(), machine.constantTrue)) {
    return false;
  }
  for (const auto* symbols : {&problem.state0Symbols, &problem.state1Symbols,
                              &problem.auxiliaryStateSymbols}) {
    for (const size_t symbol : *symbols) {
      if (!transitions.contains(symbol) || machine.registerBySymbol.contains(symbol)) {
        continue;
      }
      Register added;
      added.symbol = symbol;
      if (!circuit.add(BoolExpr::Var(symbol), added.node) ||
          !circuit.add(transitions.at(symbol), added.next)) {
        return false;
      }
      added.slot = circuit.slotOf(symbol);
      machine.registerBySymbol.emplace(symbol, static_cast<uint32_t>(machine.registers.size()));
      machine.registers.push_back(added);
    }
  }
  const size_t outputs = problem.observedOutputExprs0.size();
  machine.outputs.resize(outputs);
  const bool dualRail = problem.usesDualRailStateEncoding;
  const bool hasCones = problem.observedOutputExprs1.size() == outputs &&
      (!dualRail || (problem.dualRailObservedOutputs0.size() == outputs &&
                     problem.dualRailObservedOutputs1.size() == outputs));
  for (size_t i = 0; hasCones && i < outputs; ++i) {
    if (i < problem.structurallyProvedOutputs.size() && problem.structurallyProvedOutputs[i]) {
      continue;
    }
    std::vector<std::pair<BoolExpr*, BoolExpr*>> cones;
    if (dualRail) {
      cones.emplace_back(problem.dualRailObservedOutputs0[i].mayBeOne,
                         problem.dualRailObservedOutputs1[i].mayBeOne);
      cones.emplace_back(problem.dualRailObservedOutputs0[i].mayBeZero,
                         problem.dualRailObservedOutputs1[i].mayBeZero);
    } else {
      cones.emplace_back(problem.observedOutputExprs0[i], problem.observedOutputExprs1[i]);
    }
    for (const auto& [cone0, cone1] : cones) {
      std::pair<uint32_t, uint32_t> nodes;
      if (cone0 == nullptr || cone1 == nullptr || !circuit.add(cone0, nodes.first) ||
          !circuit.add(cone1, nodes.second)) {
        return false;
      }
      machine.outputs[i].push_back(nodes);
    }
  }
  indexLeaves(problem, machine);
  return true;
}

// The product machine of the binary functions. A register the input gives no
// function for is an input of the machine, which only weakens what it proves.
bool buildBinaryMachine(const KInductionProblem& problem,
                        const BinarySignalCorrespondenceInput& input,
                        const SignalCorrespondenceLimits& limits, Machine& machine) {
  machine.binary = true;
  Circuit& circuit = machine.circuit;
  if (input.transitions.size() > limits.stateSymbols ||
      !circuit.add(BoolExpr::createTrue(), machine.constantTrue)) {
    return false;
  }
  for (const auto& [symbol, next] : input.transitions) {
    if (next == nullptr || machine.registerBySymbol.contains(symbol)) {
      continue;
    }
    Register added;
    added.symbol = symbol;
    if (!circuit.add(BoolExpr::Var(symbol), added.node) || !circuit.add(next, added.next)) {
      return false;
    }
    added.slot = circuit.slotOf(symbol);
    machine.registerBySymbol.emplace(symbol, static_cast<uint32_t>(machine.registers.size()));
    machine.registers.push_back(added);
  }
  const size_t outputs = problem.observedOutputExprs0.size();
  machine.outputs.resize(outputs);
  for (size_t i = 0; i < outputs && i < input.outputs.size(); ++i) {
    const auto& [function0, function1] = input.outputs[i];
    if ((i < problem.structurallyProvedOutputs.size() && problem.structurallyProvedOutputs[i]) ||
        function0 == nullptr || function1 == nullptr) {
      continue;
    }
    std::pair<uint32_t, uint32_t> nodes;
    if (!circuit.add(function0, nodes.first) || !circuit.add(function1, nodes.second)) {
      return false;
    }
    machine.outputs[i].push_back(nodes);
  }
  indexLeaves(problem, machine);
  // A register is known at boot when its rails hold a binary value there.
  std::unordered_map<size_t, bool> assigned;
  for (const auto& [symbol, value] : problem.initialStateAssignments) {
    assigned[symbol] = value;
  }
  machine.boot.assign(machine.registers.size(), kUnknown);
  for (const auto& rails : problem.dualRailStatePairs) {
    const auto owner = machine.registerBySymbol.find(rails.mayBeOne);
    const auto one = assigned.find(rails.mayBeOne);
    const auto zero = assigned.find(rails.mayBeZero);
    if (owner != machine.registerBySymbol.end() && one != assigned.end() &&
        zero != assigned.end() && one->second != zero->second) {
      machine.boot[owner->second] = one->second ? 1 : 0;
    }
  }
  machine.commonRoot.resize(machine.registers.size());
  for (uint32_t i = 0; i < machine.commonRoot.size(); ++i) {
    machine.commonRoot[i] = i;
  }
  const auto find = [&machine](uint32_t index) {
    while (machine.commonRoot[index] != index) {
      machine.commonRoot[index] = machine.commonRoot[machine.commonRoot[index]];
      index = machine.commonRoot[index];
    }
    return index;
  };
  for (const auto& [lhs, rhs] : input.commonStart) {
    const auto a = machine.registerBySymbol.find(lhs);
    const auto b = machine.registerBySymbol.find(rhs);
    if (a != machine.registerBySymbol.end() && b != machine.registerBySymbol.end() &&
        machine.boot[a->second] == kUnknown && machine.boot[b->second] == kUnknown) {
      const uint32_t rootA = find(a->second);
      const uint32_t rootB = find(b->second);
      machine.commonRoot[std::max(rootA, rootB)] = std::min(rootA, rootB);
    }
  }
  for (uint32_t i = 0; i < machine.commonRoot.size(); ++i) {
    machine.commonRoot[i] = find(i);
  }
  return true;
}

struct StartStates {
  // Every register at the first compared frame: 0, 1 or unknown.
  std::vector<uint8_t> start;
  // The distinct register vectors of the frames before it.
  std::vector<std::vector<uint8_t>> earlier;
};

// The boot state, carried through the reset window with the reset ports at
// their active value and every other input unknown. An unknown register is
// left free, so the start state covers every state the window can leave.
StartStates propagateResetWindow(const Machine& machine, const KInductionProblem& problem) {
  StartStates states;
  std::vector<uint8_t> current(machine.registers.size(), kUnknown);
  for (const auto& [symbol, value] : problem.initialStateAssignments) {
    if (const auto it = machine.registerBySymbol.find(symbol);
        it != machine.registerBySymbol.end()) {
      current[it->second] = value ? 1 : 0;
    }
  }
  std::vector<uint8_t> leaves(machine.circuit.symbols().size());
  std::vector<uint8_t> values;
  for (size_t cycle = 0; cycle < machine.resetCycles; ++cycle) {
    states.earlier.push_back(current);
    for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
      const uint32_t owner = machine.registerOfSlot[slot];
      leaves[slot] = owner != kNone ? current[owner]
                   : machine.idleReset[slot] == kUnknown ? kUnknown
                   : static_cast<uint8_t>(1 - machine.idleReset[slot]);
    }
    evaluateTernary(machine.circuit, leaves, values);
    std::vector<uint8_t> next(current.size());
    for (size_t i = 0; i < next.size(); ++i) {
      next[i] = values[machine.registers[i].next];
    }
    if (next == current) {
      // Every later frame of the window is this one.
      break;
    }
    current = std::move(next);
  }
  states.start = std::move(current);
  return states;
}

// The classes of the relation being refined. A member's phase tells whether
// it equals its class or the complement, van Eijk's polarity normalisation.
struct Classes {
  std::vector<std::vector<uint32_t>> members;
  std::vector<uint8_t> phase;

  uint64_t normalised(const std::vector<uint64_t>& values, uint32_t node) const {
    return values[node] ^ (phase[node] != 0 ? kAllPatterns : 0);
  }
};

// Splits every class of `classes` by the values its members take. Parts of a
// split class stay in place and adjacent; single signals are dropped.
bool split(std::vector<std::vector<uint32_t>>& classes, const Classes& relation,
           const std::vector<uint64_t>& values) {
  bool changed = false;
  std::vector<std::vector<uint32_t>> result;
  std::vector<std::pair<uint64_t, std::vector<uint32_t>>> parts;
  for (auto& members : classes) {
    parts.clear();
    for (const uint32_t node : members) {
      const uint64_t key = relation.normalised(values, node);
      auto part = std::find_if(parts.begin(), parts.end(),
                               [key](const auto& entry) { return entry.first == key; });
      if (part == parts.end()) {
        parts.emplace_back(key, std::vector<uint32_t>{});
        part = parts.end() - 1;
      }
      part->second.push_back(node);
    }
    changed |= parts.size() > 1;
    for (auto& [_, part] : parts) {
      if (part.size() > 1) {
        result.push_back(std::move(part));
      }
    }
  }
  classes = std::move(result);
  return changed;
}

uint64_t mix(uint64_t value) {
  value ^= value >> 33;
  value *= 0xFF51AFD7ED558CCDULL;
  value ^= value >> 33;
  value *= 0xC4CEB9FE1A85EC53ULL;
  return value ^ (value >> 33);
}

// Registers that certified equalities tie together start with one value.
std::vector<uint32_t> equalityRoots(const Machine& machine, const KInductionProblem& problem) {
  std::vector<uint32_t> root(machine.registers.size());
  for (uint32_t i = 0; i < root.size(); ++i) {
    root[i] = i;
  }
  const auto find = [&root](uint32_t index) {
    while (root[index] != index) {
      root[index] = root[root[index]];
      index = root[index];
    }
    return index;
  };
  for (const auto* pairs : {&problem.sameFrameStateEqualityPairs0,
                            &problem.sameFrameStateEqualityPairs1}) {
    for (const auto& [lhs, rhs] : *pairs) {
      const auto a = machine.registerBySymbol.find(lhs);
      const auto b = machine.registerBySymbol.find(rhs);
      if (a != machine.registerBySymbol.end() && b != machine.registerBySymbol.end()) {
        const uint32_t rootA = find(a->second);
        const uint32_t rootB = find(b->second);
        root[std::max(rootA, rootB)] = std::min(rootA, rootB);
      }
    }
  }
  for (uint32_t i = 0; i < root.size(); ++i) {
    root[i] = find(i);
  }
  return root;
}

// Sequential simulation with random inputs from the start state gives the
// first partition (van Eijk, section V): signals that differ on a reachable
// state are not equivalent. The value in the first pattern fixes the phase.
//
// Every signal the simulation fails to separate costs the solver a
// counterexample later, and a deep structure, a memory behind an address
// counter for example, is only reached after many frames. A run therefore
// goes on while the partition still refines.
//
// A machine of the binary functions starts at boot, with the registers of a
// common pair on one value, and is taken through the reset window first.
Classes simulate(const Machine& machine, const KInductionProblem& problem,
                 const std::vector<uint32_t>& signals, const std::vector<uint8_t>& start) {
  Classes relation;
  relation.phase.assign(machine.circuit.size(), 0);
  const std::vector<uint32_t> root =
      machine.binary ? machine.commonRoot : equalityRoots(machine, problem);
  std::mt19937_64 random(1);
  std::vector<uint64_t> signature(signals.size(), 0);
  std::vector<uint64_t> sorted;
  std::vector<uint64_t> state(machine.registers.size());
  std::vector<uint64_t> leaves(machine.circuit.symbols().size());
  std::vector<uint64_t> values;
  const auto distinctSignatures = [&]() {
    sorted = signature;
    std::sort(sorted.begin(), sorted.end());
    return static_cast<size_t>(std::unique(sorted.begin(), sorted.end()) - sorted.begin());
  };
  uint64_t work = 0;
  for (size_t run = 0; run < kSimulationRuns && work < kSimulationWorkLimit; ++run) {
    for (size_t i = 0; i < state.size(); ++i) {
      state[i] = start[i] == kUnknown ? random() : start[i] != 0 ? kAllPatterns : 0;
    }
    for (size_t i = 0; i < state.size(); ++i) {
      if (start[i] != kUnknown) {
        state[root[i]] = state[i];
      }
    }
    for (size_t i = 0; i < state.size(); ++i) {
      state[i] = state[root[i]];
    }
    for (const auto& rails : problem.dualRailStatePairs) {
      const auto one = machine.registerBySymbol.find(rails.mayBeOne);
      const auto zero = machine.registerBySymbol.find(rails.mayBeZero);
      if (one == machine.registerBySymbol.end() || zero == machine.registerBySymbol.end()) {
        continue;
      }
      // Rails 00 are not a value.
      if (start[zero->second] == kUnknown) {
        state[zero->second] |= ~state[one->second];
      } else if (start[one->second] == kUnknown) {
        state[one->second] |= ~state[zero->second];
      }
    }
    for (size_t cycle = 0; machine.binary && cycle < machine.resetCycles; ++cycle) {
      for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
        const uint32_t owner = machine.registerOfSlot[slot];
        leaves[slot] = owner != kNone ? state[owner]
                     : machine.idleReset[slot] == kUnknown ? random()
                     : machine.idleReset[slot] != 0 ? 0 : kAllPatterns;
      }
      evaluate(machine.circuit, leaves, values);
      work += machine.circuit.size();
      for (size_t i = 0; i < state.size(); ++i) {
        state[i] = values[machine.registers[i].next];
      }
    }
    size_t classesAtLastLook = 0;
    size_t quietLooks = 0;
    for (size_t frame = 0;
         frame < kSimulationFrameLimit && quietLooks < 2 && work < kSimulationWorkLimit;
         ++frame) {
      for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
        const uint32_t owner = machine.registerOfSlot[slot];
        leaves[slot] = owner != kNone ? state[owner]
                     : machine.idleReset[slot] == kUnknown ? random()
                     : machine.idleReset[slot] != 0 ? kAllPatterns : 0;
      }
      evaluate(machine.circuit, leaves, values);
      work += machine.circuit.size();
      if (run == 0 && frame == 0) {
        for (const uint32_t node : signals) {
          relation.phase[node] = values[node] & 1;
        }
      }
      for (size_t i = 0; i < signals.size(); ++i) {
        signature[i] = mix(signature[i] + relation.normalised(values, signals[i]));
      }
      for (size_t i = 0; i < state.size(); ++i) {
        state[i] = values[machine.registers[i].next];
      }
      if ((frame + 1) % kSimulationLook == 0) {
        const size_t classes = distinctSignatures();
        quietLooks = classes == classesAtLastLook ? quietLooks + 1 : 0;
        classesAtLastLook = classes;
      }
    }
  }
  std::unordered_map<uint64_t, uint32_t> classBySignature;
  for (size_t i = 0; i < signals.size(); ++i) {
    const auto [it, inserted] = classBySignature.emplace(
        signature[i], static_cast<uint32_t>(relation.members.size()));
    if (inserted) {
      relation.members.emplace_back();
    }
    relation.members[it->second].push_back(signals[i]);
  }
  std::erase_if(relation.members, [](const auto& members) { return members.size() < 2; });
  return relation;
}

// What every frame from the first compared one on satisfies in the engines'
// model: the reset ports rest at their inactive value, no register holds the
// empty rail value, and the equalities certified earlier hold. A frame of the
// reset window holds the reset ports active; the binary functions have no
// rails, and the certified equalities are rail equalities.
void addFrameConstraints(SATSolverWrapper& solver, FrameFormulaEncoder& frame,
                         const Machine& machine, const KInductionProblem& problem,
                         bool resetActive) {
  for (const auto& [symbol, active] : machine.resetPorts) {
    const int literal = frame.encode(BoolExpr::Var(symbol));
    solver.addClause({active == resetActive ? literal : -literal});
  }
  if (machine.binary) {
    return;
  }
  for (const auto& rails : problem.dualRailStatePairs) {
    solver.addClause({frame.encode(BoolExpr::Var(rails.mayBeOne)),
                      frame.encode(BoolExpr::Var(rails.mayBeZero))});
  }
  for (const auto* pairs : {&problem.sameFrameStateEqualityPairs0,
                            &problem.sameFrameStateEqualityPairs1}) {
    for (const auto& [lhs, rhs] : *pairs) {
      addLiteralEquivalence(solver, frame.encode(BoolExpr::Var(lhs)),
                            frame.encode(BoolExpr::Var(rhs)));
    }
  }
  for (const auto& [first, second] : problem.sameFrameStateImplicationPairs) {
    solver.addClause({-frame.encode(BoolExpr::Var(first)),
                      frame.encode(BoolExpr::Var(second))});
  }
  for (const auto& [first, second] : problem.sameFrameStateDisjunctionPairs) {
    solver.addClause({frame.encode(BoolExpr::Var(first)),
                      frame.encode(BoolExpr::Var(second))});
  }
}

// One solver answers many queries, and every query adds variables. The
// profile of the PDR queries fits: no inprocessing, so the solver introduces
// no variables of its own between two queries.
void configureForManyQueries(SATSolverWrapper& solver) {
  solver.configureForSecPdrQuery();
}

// True exactly when `node` differs from `first`, the first member of its
// class, in the frame `literals` describes.
int differenceLiteral(SATSolverWrapper& solver, const std::vector<int>& literals,
                      const Classes& relation, uint32_t first, uint32_t node) {
  const bool samePhase = relation.phase[node] == relation.phase[first];
  return createXorLiteral(solver, literals[node],
                          samePhase ? literals[first] : -literals[first]);
}

// Can a member of one of these classes differ from the first member of its
// class? One query for several classes, so a round that changes nothing
// costs a query per group and not one per class.
SATSolverWrapper::SolveStatus checkGroup(
    SATSolverWrapper& solver, const std::vector<int>& literals, const Classes& relation,
    const std::vector<std::vector<uint32_t>>& classes, const SignalCorrespondenceLimits& limits,
    size_t& solverCalls) {
  std::vector<int> clause;
  for (const auto& members : classes) {
    for (size_t i = 1; i < members.size(); ++i) {
      clause.push_back(differenceLiteral(solver, literals, relation, members.front(), members[i]));
    }
  }
  const int anyDifference = solver.newVar() + 2;
  clause.push_back(-anyDifference);
  solver.addClause(clause);
  ++solverCalls;
  return solver.solveWithAssumptionsStatus(
      {anyDifference}, limits.groupConflicts, limits.groupDecisions, kNoTickLimit);
}

enum class ClassStatus {
  // No member can differ from the first one.
  Proved,
  // The solver holds a model in which a member differs.
  Refuted,
};

// Asks whether a member of the class can differ from its first member in the
// frame `literals` describes. Members the solver cannot decide are removed
// and reported through `dropped`: an undecided signal is not assumed.
ClassStatus checkClass(SATSolverWrapper& solver, const std::vector<int>& literals,
                       const Classes& relation, const SignalCorrespondenceLimits& limits,
                       std::vector<uint32_t>& members, bool& dropped, size_t& solverCalls) {
  const uint32_t first = members.front();
  std::vector<int> differences;
  std::vector<int> clause;
  for (size_t i = 1; i < members.size(); ++i) {
    differences.push_back(differenceLiteral(solver, literals, relation, first, members[i]));
    clause.push_back(differences.back());
  }
  const int anyDifference = solver.newVar() + 2;
  clause.push_back(-anyDifference);
  solver.addClause(clause);
  ++solverCalls;
  const auto status = solver.solveWithAssumptionsStatus(
      {anyDifference}, limits.classConflicts, limits.classDecisions, kNoTickLimit);
  if (status == SATSolverWrapper::SolveStatus::Sat) {
    return ClassStatus::Refuted;
  }
  if (status == SATSolverWrapper::SolveStatus::Unsat) {
    return ClassStatus::Proved;
  }
  // Undecided together: each member is then its own obligation.
  std::vector<uint32_t> kept{first};
  for (size_t i = 1; i < members.size(); ++i) {
    ++solverCalls;
    const auto memberStatus = solver.solveWithAssumptionsStatus(
        {differences[i - 1]}, limits.classConflicts, limits.classDecisions, kNoTickLimit);
    if (memberStatus == SATSolverWrapper::SolveStatus::Sat) {
      kept.insert(kept.end(), members.begin() + i, members.end());
      members = std::move(kept);
      return ClassStatus::Refuted;
    }
    if (memberStatus == SATSolverWrapper::SolveStatus::Unsat) {
      kept.push_back(members[i]);
    } else {
      dropped = true;
    }
  }
  members = std::move(kept);
  return ClassStatus::Proved;
}

// Keeps the classes whose members are equal in the frame `literals`
// describes, splitting on every counterexample. `replay` turns the model the
// solver holds into the values the signals take in that frame; `proved` is
// told every class that is kept. `changed` is set when a class was split or
// lost a member. False when the solver budget ran out.
//
// Classes are asked in groups. A group the solver cannot decide in one query
// is halved, down to single classes, which are asked member by member.
template <typename Replay, typename Proved>
bool refineClasses(SATSolverWrapper& solver, const std::vector<int>& literals,
                   const SignalCorrespondenceLimits& limits,
                   const SATSolverWrapper::CadicalWorkBudget& budget, Classes& relation,
                   bool& changed, size_t& solverCalls, Replay replay, Proved proved) {
  std::vector<std::vector<uint32_t>> pending = std::move(relation.members);
  std::vector<std::vector<uint32_t>> kept;
  std::vector<std::vector<std::vector<uint32_t>>> groups;
  std::vector<uint64_t> values;
  const auto keep = [&](std::vector<uint32_t>& members) {
    if (members.size() > 1) {
      proved(members);
      kept.push_back(std::move(members));
    }
  };
  // One class, member by member if the solver cannot decide it at once.
  const auto single = [&](std::vector<uint32_t>& members) {
    bool dropped = false;
    const auto status =
        checkClass(solver, literals, relation, limits, members, dropped, solverCalls);
    changed |= dropped;
    if (status == ClassStatus::Proved) {
      keep(members);
      return;
    }
    changed = true;
    replay(values);
    pending.push_back(std::move(members));
    if (!split(pending, relation, values)) {
      // The replay disagrees with the solver; the class is not kept.
      pending.pop_back();  // LCOV_EXCL_LINE
    }
  };
  while (!pending.empty()) {
    const size_t count = std::min(std::max<size_t>(limits.groupClasses, 1), pending.size());
    groups.assign(1, {});
    groups.back().assign(std::make_move_iterator(pending.end() - count),
                         std::make_move_iterator(pending.end()));
    pending.resize(pending.size() - count);
    while (!groups.empty()) {
      std::vector<std::vector<uint32_t>> group = std::move(groups.back());
      groups.pop_back();
      if (group.size() == 1) {
        single(group.front());
      } else {
        const auto status = checkGroup(solver, literals, relation, group, limits, solverCalls);
        if (status == SATSolverWrapper::SolveStatus::Unsat) {
          for (auto& members : group) {
            keep(members);
          }
        } else if (status == SATSolverWrapper::SolveStatus::Sat) {
          // The counterexample splits a class of this group, and whatever
          // else it separates among the classes still to ask.
          replay(values);
          changed = true;
          for (auto& waiting : groups) {
            std::move(waiting.begin(), waiting.end(), std::back_inserter(pending));
          }
          groups.clear();
          const size_t first = pending.size();
          std::move(group.begin(), group.end(), std::back_inserter(pending));
          if (!split(pending, relation, values)) {
            // The replay disagrees with the solver; ask the group class by class.
            // LCOV_EXCL_START
            std::vector<std::vector<uint32_t>> again(
                std::make_move_iterator(pending.begin() + first),
                std::make_move_iterator(pending.end()));
            pending.resize(first);
            for (auto& members : again) {
              single(members);
            }
            // LCOV_EXCL_STOP
          }
        } else {
          const auto middle = group.begin() + group.size() / 2;
          groups.emplace_back(std::make_move_iterator(group.begin()),
                              std::make_move_iterator(middle));
          groups.emplace_back(std::make_move_iterator(middle),
                              std::make_move_iterator(group.end()));
        }
      }
      if (budget.exhausted()) {
        return false;
      }
    }
  }
  relation.members = std::move(kept);
  return true;
}

// The solver's model of one frame as simulation leaves, every pattern alike.
// A leaf the frame never encoded keeps the value it has in `leaves`.
void readModel(const SATSolverWrapper& solver, const FrameFormulaEncoder& frame,
               const Circuit& circuit, std::vector<uint64_t>& leaves) {
  for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
    const auto literal = frame.leafLits().find(circuit.symbols()[slot]);
    if (literal != frame.leafLits().end()) {
      leaves[slot] = solver.getLiteralValue(literal->second) ? kAllPatterns : 0;
    }
  }
}

// Consecutive time frames of the machine in one solver. Frame 0 is free; a
// register of every later frame is its next-state function in the frame
// before it. The first `resetFrames` frames are frames of the reset window.
struct Unrolling {
  std::deque<FrameFormulaEncoder> frames;
  // Per frame, the literal of every signal of the relation.
  std::vector<std::vector<int>> literals;
};

void unroll(SATSolverWrapper& solver, const Machine& machine, const KInductionProblem& problem,
            const Classes& relation, size_t count, Unrolling& unrolling,
            size_t resetFrames = 0) {
  for (size_t f = 0; f < count; ++f) {
    auto& frame = unrolling.frames.emplace_back(solver, std::unordered_map<size_t, int>{}, true);
    auto& literals = unrolling.literals.emplace_back(machine.circuit.size(), 0);
    for (const auto& members : relation.members) {
      for (const uint32_t node : members) {
        literals[node] = frame.encode(machine.circuit.expr(node));
      }
    }
    addFrameConstraints(solver, frame, machine, problem, f < resetFrames);
  }
  // Last frame first: tying a frame's registers to the frame before it adds
  // registers there.
  for (size_t f = count; f-- > 1;) {
    const std::vector<std::pair<size_t, int>> leaves(unrolling.frames[f].leafLits().begin(),
                                                     unrolling.frames[f].leafLits().end());
    for (const auto& [symbol, literal] : leaves) {
      if (const auto it = machine.registerBySymbol.find(symbol);
          it != machine.registerBySymbol.end()) {
        addLiteralEquivalence(solver, literal,
                              unrolling.frames[f - 1].encode(
                                  machine.circuit.expr(machine.registers[it->second].next)));
      }
    }
  }
}

// The model the solver holds, replayed on the logic from frame 0 to frame
// `last`. Frame 0 reads its leaves from the model; a leaf that frame never
// encoded keeps the value it has in `leaves`. A later frame takes its
// registers from the frame before it and its inputs from the model.
void replay(const SATSolverWrapper& solver, const Machine& machine, const Unrolling& unrolling,
            size_t last, std::vector<uint64_t>& leaves, std::vector<uint64_t>& values) {
  readModel(solver, unrolling.frames[0], machine.circuit, leaves);
  evaluate(machine.circuit, leaves, values);
  for (size_t f = 1; f <= last; ++f) {
    for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
      leaves[slot] = machine.idleReset[slot] == 1 ? kAllPatterns : 0;
    }
    readModel(solver, unrolling.frames[f], machine.circuit, leaves);
    for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
      if (const uint32_t owner = machine.registerOfSlot[slot]; owner != kNone) {
        leaves[slot] = values[machine.registers[owner].next];
      }
    }
    evaluate(machine.circuit, leaves, values);
  }
}

// Equation 4, over `frames` time frames: members of a class are equal in the
// start state and in the frames that follow it, for every input. False when
// the solver budget ran out.
//
// A machine of the binary functions is given its boot state: the reset
// window is unrolled in front of the frames the classes are checked in, so
// the start state is exactly what the window leaves.
bool refineAtStart(const Machine& machine, const KInductionProblem& problem,
                   const std::vector<uint8_t>& start, size_t frames,
                   Config::SolverType solverType, const SignalCorrespondenceLimits& limits,
                   const SATSolverWrapper::CadicalWorkBudget& budget, Classes& relation,
                   size_t& solverCalls) {
  SATSolverWrapper solver(SATSolverWrapper::assumptionSolverTypeFor(solverType));
  configureForManyQueries(solver);
  Unrolling unrolling;
  const size_t window = machine.binary ? machine.resetCycles : 0;
  unroll(solver, machine, problem, relation, window + frames, unrolling, window);
  const auto bootLiteral = [&](size_t index) {
    return unrolling.frames[0].encode(BoolExpr::Var(machine.registers[index].symbol));
  };
  for (size_t i = 0; i < machine.registers.size(); ++i) {
    if (start[i] != kUnknown) {
      const int literal = bootLiteral(i);
      solver.addClause({start[i] != 0 ? literal : -literal});
    } else if (machine.binary && machine.commonRoot[i] != i) {
      addLiteralEquivalence(solver, bootLiteral(i), bootLiteral(machine.commonRoot[i]));
    }
  }

  std::vector<uint64_t> leaves(machine.circuit.symbols().size());
  for (size_t f = 0; f < frames; ++f) {
    const auto replayFrame = [&](std::vector<uint64_t>& values) {
      // A leaf outside every encoded cone: a register keeps its start value
      // (unknown reads 1: X when on both rails, and one value for the
      // registers of a common pair), an input reads 0.
      for (uint32_t slot = 0; slot < leaves.size(); ++slot) {
        const uint32_t owner = machine.registerOfSlot[slot];
        leaves[slot] = owner != kNone && start[owner] != 0 ? kAllPatterns : 0;
      }
      replay(solver, machine, unrolling, window + f, leaves, values);
    };
    bool changed = false;
    if (!refineClasses(solver, unrolling.literals[window + f], limits, budget, relation, changed,
                       solverCalls, replayFrame, [](const std::vector<uint32_t>&) {})) {
      return false;
    }
  }
  return true;
}

// Equation 5, once, over `frames` time frames: with every class assumed in
// that many consecutive frames, the members of a class are equal in the
// next. `changed` tells whether a class was split or lost a member, in which
// case the assumption was too strong and the step is repeated. False when the
// solver budget ran out.
bool refineOverSteps(const Machine& machine, const KInductionProblem& problem, size_t frames,
                     Config::SolverType solverType, const SignalCorrespondenceLimits& limits,
                     const SATSolverWrapper::CadicalWorkBudget& budget, Classes& relation,
                     bool& changed, size_t& solverCalls) {
  SATSolverWrapper solver(SATSolverWrapper::assumptionSolverTypeFor(solverType));
  configureForManyQueries(solver);
  Unrolling unrolling;
  unroll(solver, machine, problem, relation, frames + 1, unrolling);
  // The state compatibility condition of the relation: plain equalities of
  // the assumed frames. Nothing is merged.
  const auto assume = [&](const std::vector<int>& literals,
                          const std::vector<uint32_t>& members) {
    const uint32_t first = members.front();
    for (size_t i = 1; i < members.size(); ++i) {
      const bool samePhase = relation.phase[members[i]] == relation.phase[first];
      addLiteralEquivalence(solver, literals[members[i]],
                            samePhase ? literals[first] : -literals[first]);
    }
  };
  for (size_t f = 0; f < frames; ++f) {
    for (const auto& members : relation.members) {
      assume(unrolling.literals[f], members);
    }
  }

  std::vector<uint64_t> leaves(machine.circuit.symbols().size());
  const auto replayStep = [&](std::vector<uint64_t>& values) {
    std::fill(leaves.begin(), leaves.end(), 0);
    replay(solver, machine, unrolling, frames, leaves, values);
  };
  // A class proved under this round's assumption may be used by the queries
  // that follow in the round.
  const auto& proved = unrolling.literals[frames];
  return refineClasses(solver, proved, limits, budget, relation, changed, solverCalls, replayStep,
                       [&](const std::vector<uint32_t>& members) { assume(proved, members); });
}

std::vector<uint32_t> selectSignals(const Machine& machine, SignalCorrespondenceScope scope) {
  std::vector<uint32_t> signals;
  if (scope == SignalCorrespondenceScope::AllSignals) {
    signals.resize(machine.circuit.size());
    for (uint32_t i = 0; i < signals.size(); ++i) {
      signals[i] = i;
    }
    return signals;
  }
  signals.push_back(machine.constantTrue);
  for (const auto& entry : machine.registers) {
    signals.push_back(entry.node);
    signals.push_back(entry.next);
  }
  for (const auto& pairs : machine.outputs) {
    for (const auto& [node0, node1] : pairs) {
      signals.push_back(node0);
      signals.push_back(node1);
    }
  }
  std::sort(signals.begin(), signals.end());
  signals.erase(std::unique(signals.begin(), signals.end()), signals.end());
  return signals;
}

// Registers proved equal that may be handed to the engines: complete values
// (both rails in the dual-rail encoding, plain or crosswise for a complement)
// that are also known and equal in every frame before the start state.
void collectStateEqualities(const Machine& machine, const KInductionProblem& problem,
                            const StartStates& states, const Classes& relation,
                            const std::vector<uint32_t>& classOf,
                            SignalCorrespondenceResult& result) {
  const auto heldEarlier = [&states](uint32_t a, uint32_t b) {
    return std::all_of(states.earlier.begin(), states.earlier.end(),
                       [a, b](const std::vector<uint8_t>& frame) {
                         return frame[a] != kUnknown && frame[a] == frame[b];
                       });
  };
  // Already certified pairs are not reported again.
  const std::vector<uint32_t> root = equalityRoots(machine, problem);
  const auto report = [&](uint32_t a, uint32_t b) {
    if (root[a] != root[b]) {
      result.stateEqualities.emplace_back(machine.registers[a].symbol,
                                          machine.registers[b].symbol);
    }
  };
  using Key = std::array<uint32_t, 4>;
  const auto keyOf = [&](uint32_t one, uint32_t zero) {
    const uint32_t nodeOne = machine.registers[one].node;
    const uint32_t nodeZero = machine.registers[zero].node;
    return Key{classOf[nodeOne], relation.phase[nodeOne], classOf[nodeZero],
               relation.phase[nodeZero]};
  };
  std::map<Key, std::pair<uint32_t, uint32_t>> firstByKey;
  if (problem.usesDualRailStateEncoding) {
    for (const auto& rails : problem.dualRailStatePairs) {
      const auto oneIt = machine.registerBySymbol.find(rails.mayBeOne);
      const auto zeroIt = machine.registerBySymbol.find(rails.mayBeZero);
      if (oneIt == machine.registerBySymbol.end() || zeroIt == machine.registerBySymbol.end()) {
        continue;
      }
      const uint32_t one = oneIt->second;
      const uint32_t zero = zeroIt->second;
      const Key key = keyOf(one, zero);
      if (key[0] == kNone || key[2] == kNone) {
        continue;
      }
      if (const auto same = firstByKey.find(key); same != firstByKey.end()) {
        if (heldEarlier(same->second.first, one) && heldEarlier(same->second.second, zero)) {
          report(same->second.first, one);
          report(same->second.second, zero);
        }
        continue;
      }
      // The complement of a value is its rails exchanged.
      if (const auto crossed = firstByKey.find(Key{key[2], key[3], key[0], key[1]});
          crossed != firstByKey.end()) {
        if (heldEarlier(crossed->second.first, zero) &&
            heldEarlier(crossed->second.second, one)) {
          report(crossed->second.first, zero);
          report(crossed->second.second, one);
        }
        continue;
      }
      firstByKey.emplace(key, std::pair{one, zero});
    }
    return;
  }
  std::map<std::pair<uint32_t, uint8_t>, uint32_t> firstByClass;
  for (uint32_t i = 0; i < machine.registers.size(); ++i) {
    const uint32_t node = machine.registers[i].node;
    if (classOf[node] == kNone) {
      continue;
    }
    const auto [it, inserted] =
        firstByClass.emplace(std::pair{classOf[node], relation.phase[node]}, i);
    if (!inserted && heldEarlier(it->second, i)) {
      report(it->second, i);
    }
  }
}

}  // namespace

SignalCorrespondenceResult proveSignalCorrespondence(
    const KInductionProblem& problem,
    SignalCorrespondenceScope scope,
    Config::SolverType solverType,
    const SignalCorrespondenceLimits& limits,
    size_t inductionFrames,
    const BinarySignalCorrespondenceInput* binary) {
  SignalCorrespondenceResult result;
  result.provedOutputs.assign(problem.observedOutputExprs0.size(), false);
  result.differingOutputs.assign(problem.observedOutputExprs0.size(), false);
  const size_t frames = std::max<size_t>(inductionFrames, 1);
  Machine machine(limits.nodes);
  // The step encodes one frame more than it assumes.
  if (!(binary != nullptr ? buildBinaryMachine(problem, *binary, limits, machine)
                          : buildMachine(problem, limits, machine)) ||
      machine.circuit.size() > limits.nodes / (frames + 1) ||
      (scope == SignalCorrespondenceScope::AllSignals &&
       machine.circuit.size() > limits.allSignalsNodes)) {
    result.stoppedBy = SignalCorrespondenceStop::Model;
    return result;
  }
  const std::vector<uint32_t> signals = selectSignals(machine, scope);
  result.signals = signals.size();
  const StartStates states = machine.binary ? StartStates{machine.boot, {}}
                                            : propagateResetWindow(machine, problem);
  Classes relation = simulate(machine, problem, signals, states.start);

  SATSolverWrapper::CadicalWorkBudget budget(
      limits.budgetConflicts, limits.budgetDecisions, limits.budgetTicks);
  SATSolverWrapper::ScopedCadicalWorkBudget budgetScope(budget);
  // The class of every signal, and the output cones the relation holds
  // equal, as it stands when indexed.
  std::vector<uint32_t> classOf;
  const auto indexClasses = [&]() {
    classOf.assign(machine.circuit.size(), kNone);
    for (uint32_t i = 0; i < relation.members.size(); ++i) {
      for (const uint32_t node : relation.members[i]) {
        classOf[node] = i;
      }
    }
  };
  const auto equalOutputs = [&](size_t output) {
    const auto& pairs = machine.outputs[output];
    return !pairs.empty() &&
           std::all_of(pairs.begin(), pairs.end(), [&](const auto& pair) {
             return pair.first == pair.second ||
                    (classOf[pair.first] != kNone && classOf[pair.first] == classOf[pair.second] &&
                     relation.phase[pair.first] == relation.phase[pair.second]);
           });
  };
  bool withinBudget = refineAtStart(machine, problem, states.start, frames, solverType, limits,
                                    budget, relation, result.solverCalls);
  indexClasses();
  for (size_t i = 0; withinBudget && i < machine.outputs.size(); ++i) {
    // Separated before any induction step: a difference in a reachable state.
    result.differingOutputs[i] = !machine.outputs[i].empty() && !equalOutputs(i);
  }
  bool fixedPoint = false;
  while (withinBudget && !fixedPoint && result.rounds < limits.rounds) {
    bool changed = false;
    withinBudget = refineOverSteps(machine, problem, frames, solverType, limits, budget,
                                   relation, changed, result.solverCalls);
    ++result.rounds;
    fixedPoint = withinBudget && !changed;
  }
  result.conflicts = budget.conflictsUsed();
  result.decisions = budget.decisionsUsed();
  result.ticks = budget.ticksUsed();
  if (!fixedPoint) {
    // An unfinished fixed point is not a certificate for any class.
    result.stoppedBy =
        withinBudget ? SignalCorrespondenceStop::Rounds
        : result.ticks >= limits.budgetTicks ? SignalCorrespondenceStop::BudgetTicks
        : result.conflicts >= limits.budgetConflicts ? SignalCorrespondenceStop::BudgetConflicts
        : SignalCorrespondenceStop::BudgetDecisions;
    return result;
  }

  result.completed = true;
  result.classes = relation.members.size();
  indexClasses();
  for (size_t i = 0; i < machine.outputs.size(); ++i) {
    result.provedOutputs[i] = equalOutputs(i);
  }
  if (!machine.binary) {
    collectStateEqualities(machine, problem, states, relation, classOf, result);
  }
  return result;
}

namespace {

const char* stopName(SignalCorrespondenceStop stop) {
  static constexpr const char* names[] = {"none", "model", "rounds",
                                          "conflicts", "decisions", "ticks"};
  return names[static_cast<size_t>(stop)];
}

}  // namespace

SignalCorrespondenceLearning learnSignalCorrespondence(
    const KInductionProblem& problem,
    Config::SolverType solverType,
    bool diagnostics,
    SignalCorrespondenceLimits limits,
    const BinarySignalCorrespondenceInput* binary) {
  SignalCorrespondenceLearning learned;
  if (binary != nullptr) {
    limits.inductionFrames = limits.binaryInductionFrames;
    limits.budgetTicks = limits.binaryBudgetTicks;
  }
  const size_t outputs = problem.observedOutputExprs0.size();
  learned.provedOutputs.assign(outputs, false);
  // Outputs whose cones differ in a reachable state: no run can prove them.
  std::vector<bool> differing(outputs, false);
  const auto open = [&]() {
    size_t count = 0;
    for (size_t i = 0; i < outputs; ++i) {
      const bool earlier = i < problem.structurallyProvedOutputs.size() &&
                           problem.structurallyProvedOutputs[i];
      count += !earlier && !learned.provedOutputs[i] && !differing[i] ? 1 : 0;
    }
    return count;
  };
  // One budget for all the runs: each run gets what the earlier ones left.
  const auto spend = [](uint64_t& left, uint64_t used) { left = used < left ? left - used : 0; };
  const auto budgetLeft = [&limits]() {
    return limits.budgetConflicts != 0 && limits.budgetDecisions != 0 && limits.budgetTicks != 0;
  };
  for (const auto scope : {SignalCorrespondenceScope::RegistersAndOutputs,
                           SignalCorrespondenceScope::AllSignals}) {
    if (binary != nullptr && scope == SignalCorrespondenceScope::AllSignals) {
      break;
    }
    // Deeper induction proves more and costs more; a depth that gave up is
    // not followed by a deeper one.
    for (size_t frames = 1;
         frames <= limits.inductionFrames && open() != 0 && budgetLeft(); frames *= 2) {
      const size_t openBefore = open();
      const auto begin = std::chrono::steady_clock::now();
      const auto result =
          proveSignalCorrespondence(problem, scope, solverType, limits, frames, binary);
      spend(limits.budgetConflicts, result.conflicts);
      spend(limits.budgetDecisions, result.decisions);
      spend(limits.budgetTicks, result.ticks);
      const size_t provedBefore =
          std::count(learned.provedOutputs.begin(), learned.provedOutputs.end(), true);
      for (size_t i = 0; i < outputs; ++i) {
        learned.provedOutputs[i] = learned.provedOutputs[i] || result.provedOutputs[i];
        differing[i] = differing[i] || result.differingOutputs[i];
      }
      // Every completed run reports invariants; duplicates are harmless.
      learned.stateEqualities.insert(learned.stateEqualities.end(),
                                     result.stateEqualities.begin(),
                                     result.stateEqualities.end());
      if (diagnostics) {
        printf("SEC summary: signal_correspondence scope=%s frames=%zu completed=%d signals=%zu "
               "classes=%zu rounds=%zu solver_calls=%zu proved_outputs=%zu/%zu "
               "state_equalities=%zu seconds=%.1f semantics=%s stopped_by=%s\n",
               scope == SignalCorrespondenceScope::AllSignals ? "all_signals"
                                                              : "registers_and_outputs",
               frames, result.completed, result.signals, result.classes, result.rounds,
               result.solverCalls,
               static_cast<size_t>(std::count(learned.provedOutputs.begin(),
                                              learned.provedOutputs.end(), true)) - provedBefore,
               openBefore,
               result.stateEqualities.size(),
               std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count(),
               binary != nullptr ? "binary" : "rails", stopName(result.stoppedBy));
        fflush(stdout);
      }
      if (!result.completed) {
        break;
      }
    }
  }
  std::sort(learned.stateEqualities.begin(), learned.stateEqualities.end());
  learned.stateEqualities.erase(
      std::unique(learned.stateEqualities.begin(), learned.stateEqualities.end()),
      learned.stateEqualities.end());
  return learned;
}

}  // namespace KEPLER_FORMAL::SEC
