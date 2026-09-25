// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "latch/LatchSymbolicCompiler.h"

#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "kinduction/SatEncoding.h"

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

struct Failure {
  CertificationStatus status;
  std::string detail;
};

[[noreturn]] void unproved(const std::string& detail) {
  throw Failure{CertificationStatus::UnprovedBound, detail};
}

[[noreturn]] void resource(const std::string& detail) {
  throw Failure{CertificationStatus::ResourceLimit, detail};
}

// Count the cumulative DAG footprint, including proof-only copies. Traversal is
// iterative and shared nodes are visited once, so deeply unfolded paths do not
// consume the C++ stack. All leaves must belong to this compilation's allocator.
class DagBudget {
 public:
  DagBudget(size_t limit, const size_t& nextSymbol)
      : limit_(limit), nextSymbol_(nextSymbol) {}

  void inspect(const SymbolicBits& roots) {
    std::vector<BoolExpr*> pending(roots.begin(), roots.end());
    while (!pending.empty()) {
      auto* node = pending.back();
      pending.pop_back();
      if (!node || !node->isValid()) {
        throw std::invalid_argument("Invalid Boolean DAG in symbolic latch reference");
      }
      if (!seen_.insert(node).second) continue;
      if (seen_.size() > limit_) resource("Symbolic latch DAG exceeds maxNodes");
      if (node->getOp() == Op::VAR) {
        if (node->getId() >= nextSymbol_) {
          throw std::invalid_argument("Symbolic primitive introduced an unallocated leaf");
        }
      } else {
        pending.push_back(node->getLeft());
        if (node->getRight()) pending.push_back(node->getRight());
      }
    }
  }

  void inspect(BoolExpr* root) { inspect(SymbolicBits{root}); }

  void inspect(const SymbolicState& state) {
    inspect(state.current);
    inspect(state.previous);
    for (const auto& storage : state.storage) inspect(storage);
    inspect(state.active);
    inspect(state.bootstrap);
    inspect(state.error);
  }

 private:
  size_t limit_;
  const size_t& nextSymbol_;
  std::unordered_set<BoolExpr*> seen_;
};

SymbolicBits constants(const Bits& values) {
  SymbolicBits result;
  result.reserve(values.size());
  for (const auto bit : values) {
    if (bit > 1) throw std::invalid_argument("Non-Boolean symbolic initialization");
    result.push_back(BoolExpr::Var(bit));
  }
  return result;
}

BoolExpr* differs(const SymbolicBits& left, const SymbolicBits& right) {
  if (left.size() != right.size()) {
    throw std::invalid_argument("Symbolic boundary layout changed during propagation");
  }
  auto* difference = BoolExpr::createFalse();
  for (size_t i = 0; i < left.size(); ++i) {
    difference = BoolExpr::Or(difference, BoolExpr::Xor(left[i], right[i]));
  }
  return difference;
}

BoolExpr* allowedInput(const Network& network, const SymbolicState& boundary,
                       const SymbolicBits& inputs, bool singleChange) {
  if (!singleChange) return BoolExpr::createTrue();
  // At-most-one changed original external bit, in linear-size prefix form.
  // The condition concerns admission, never the internal pin-change schedule.
  auto* seenChange = BoolExpr::createFalse();
  auto* multiple = BoolExpr::createFalse();
  for (size_t i = 0; i < inputs.size(); ++i) {
    auto* change = BoolExpr::Xor(boundary.current.at(network.externalInputs[i]), inputs[i]);
    multiple = BoolExpr::Or(multiple, BoolExpr::And(seenChange, change));
    seenChange = BoolExpr::Or(seenChange, change);
  }
  return BoolExpr::Not(multiple);
}

// Evaluate ground BOOT expressions without recursive BoolExpr::evaluate.
Bits groundValues(const SymbolicBits& roots) {
  std::unordered_map<BoolExpr*, bool> values;
  std::vector<std::pair<BoolExpr*, bool>> pending;
  for (auto* root : roots) {
    pending.emplace_back(root, false);
    while (!pending.empty()) {
      const auto [node, visited] = pending.back();
      pending.pop_back();
      if (values.contains(node)) continue;
      if (node->getOp() == Op::VAR) {
        if (node->getId() > 1) {
          throw std::invalid_argument("Canonical symbolic BOOT still contains a free input");
        }
        values.emplace(node, node->getId() == 1);
      } else if (!visited) {
        pending.emplace_back(node, true);
        if (node->getRight()) pending.emplace_back(node->getRight(), false);
        pending.emplace_back(node->getLeft(), false);
      } else {
        const bool left = values.at(node->getLeft());
        switch (node->getOp()) {
          case Op::NOT: values.emplace(node, !left); break;
          case Op::AND: values.emplace(node, left && values.at(node->getRight())); break;
          case Op::OR: values.emplace(node, left || values.at(node->getRight())); break;
          case Op::XOR: values.emplace(node, left != values.at(node->getRight())); break;
          default: throw std::invalid_argument("Unsupported symbolic Boolean operator");
        }
      }
    }
  }
  Bits result;
  for (auto* root : roots) result.push_back(values.at(root));
  return result;
}

void checkFinalLeaves(const SymbolicBits& roots, const SymbolicMacro& macro) {
  std::unordered_set<size_t> retained(macro.stateSymbols.begin(), macro.stateSymbols.end());
  retained.insert(macro.inputSymbols.begin(), macro.inputSymbols.end());
  std::unordered_set<BoolExpr*> seen;
  std::vector<BoolExpr*> pending(roots.begin(), roots.end());
  while (!pending.empty()) {
    auto* node = pending.back();
    pending.pop_back();
    if (!seen.insert(node).second) continue;
    if (node->getOp() == Op::VAR) {
      if (node->getId() > 1 && !retained.contains(node->getId())) {
        throw std::invalid_argument("Uneliminated seed or ordering choice in compiled macrostate");
      }
    } else {
      pending.push_back(node->getLeft());
      if (node->getRight()) pending.push_back(node->getRight());
    }
  }
}

class Compiler {
 public:
  Compiler(const SymbolicNetwork& network, const SymbolicCompileOptions& options)
      : network_(network), options_(options), model_(network, {}, options.workers),
        dag_(options.maxNodes, nextSymbol_),
        satBudget_(options.maxSatConflicts, options.maxSatDecisions,
                   std::numeric_limits<uint64_t>::max()) {}

  SymbolicMacro run() {
    if (!options_.maxNodes || !options_.maxSatConflicts || !options_.maxSatDecisions) {
      resource("Symbolic certification requires a nonzero DAG and SAT work budget");
    }
    const auto& reference = network_.reference;
    if (options_.initialInputs.size() != reference.externalInputs.size() ||
        options_.initialStorage.size() != reference.primitives.size()) {
      throw std::invalid_argument("Explicit symbolic initialization has the wrong shape");
    }
    const auto bootInputs = constants(options_.initialInputs);
    std::vector<SymbolicBits> bootStorage;
    for (size_t i = 0; i < reference.primitives.size(); ++i) {
      if (options_.initialStorage[i].size() != reference.primitives[i].storageBits) {
        throw std::invalid_argument("Explicit symbolic primitive storage has the wrong width");
      }
      bootStorage.push_back(constants(options_.initialStorage[i]));
    }

    SymbolicMacro result;
    result.singleExternalInputChange = options_.singleExternalInputChange;
    result.externalInputNets = reference.externalInputs;
    SymbolicBits current = variables(reference.netCount, &result.stateSymbols);
    std::vector<SymbolicBits> storage;
    for (const auto& primitive : reference.primitives) {
      storage.push_back(variables(primitive.storageBits, &result.stateSymbols));
    }
    const auto input = variables(reference.externalInputs.size(), &result.inputSymbols);
    const auto boundary = model_.boundary(current, storage);
    dag_.inspect(boundary);
    auto* invariant = model_.boundaryInvariant(boundary);
    dag_.inspect(invariant);

    // Intended BOOT inputs/storage are shared, auxiliary net seeds and all wave
    // choices are independent. The second copy is not tied to the first entry.
    auto bootA = model_.bootstrap(bootInputs, bootStorage, variables(reference.netCount));
    auto bootB = model_.bootstrap(bootInputs, bootStorage, variables(reference.netCount));
    result.bootstrapWaves = settle(bootA, BoolExpr::createTrue(), "BOOT progress");
    advance(bootB, result.bootstrapWaves, false);
    requireUnsat(BoolExpr::Not(model_.stable(bootB)), "independent BOOT progress");
    requireUnsat(differs(flattenSymbolicBoundary(bootA), flattenSymbolicBoundary(bootB)),
                 "BOOT seed/order independence", CertificationStatus::OrderDependent);
    requireUnsat(BoolExpr::Not(model_.boundaryInvariant(bootA)),
                 "BOOT does not establish the candidate boundary invariant");

    auto* assumption = BoolExpr::And(invariant, allowedInput(
        reference, boundary, input, options_.singleExternalInputChange));
    dag_.inspect(assumption);
    // Admission is a total functional construction, including invalid/error
    // handling. Both copies share q/u, but never any internal ordering choices.
    auto nextA = model_.admit(boundary, input);
    auto nextB = model_.admit(boundary, input);
    result.transitionWaves = settle(nextA, assumption, "boundary episode progress");
    advance(nextB, result.transitionWaves, false);
    requireUnsat(BoolExpr::And(assumption, BoolExpr::Not(model_.stable(nextB))),
                 "independent boundary episode progress");
    requireUnsat(BoolExpr::And(assumption,
                 differs(flattenSymbolicBoundary(nextA), flattenSymbolicBoundary(nextB))),
                 "boundary outcome uniqueness over the candidate invariant");
    requireUnsat(BoolExpr::And(assumption, BoolExpr::Not(model_.boundaryInvariant(nextA))),
                 "boundary invariant closure");

    // A legal canonical pin ordering is selected only after universal progress,
    // uniqueness and inductive closure succeeded. Unused choice codes are legal
    // in the reference; assigning zero therefore removes no possible behavior.
    auto canonicalBoot = model_.bootstrap(bootInputs, bootStorage,
        SymbolicBits(reference.netCount, BoolExpr::createFalse()));
    advance(canonicalBoot, result.bootstrapWaves, true);
    result.initialState = groundValues(flattenSymbolicBoundary(canonicalBoot));
    auto canonicalNext = model_.admit(boundary, input);
    advance(canonicalNext, result.transitionWaves, true);
    result.nextState = flattenSymbolicBoundary(canonicalNext);
    result.observedNets = canonicalNext.current;
    if (result.initialState.size() != result.stateSymbols.size() ||
        result.nextState.size() != result.stateSymbols.size()) {
      throw std::invalid_argument("Symbolic macro layout does not preserve the complete boundary");
    }
    dag_.inspect(result.nextState);
    checkFinalLeaves(result.nextState, result);
    return result;
  }

 private:
  BoolExpr* fresh() {
    if (nextSymbol_ == std::numeric_limits<size_t>::max()) {
      resource("Symbolic variable identifier space exhausted");
    }
    auto* symbol = BoolExpr::Var(nextSymbol_++);
    dag_.inspect(symbol);
    return symbol;
  }

  SymbolicBits variables(size_t count, std::vector<size_t>* ids = nullptr) {
    if (count > options_.maxNodes) resource("Symbolic interface exceeds the DAG budget");
    SymbolicBits result;
    result.reserve(count);
    for (size_t i = 0; i < count; ++i) {
      auto* variable = fresh();
      result.push_back(variable);
      if (ids) ids->push_back(variable->getId());
    }
    return result;
  }

  SATSolverWrapper::SolveStatus solve(BoolExpr* formula) {
    dag_.inspect(formula);
    if (formula == BoolExpr::createFalse()) return SATSolverWrapper::SolveStatus::Unsat;
    if (formula == BoolExpr::createTrue()) return SATSolverWrapper::SolveStatus::Sat;
    if (satBudget_.exhausted()) resource("Symbolic settling cumulative SAT budget exhausted");
    SATSolverWrapper::ScopedCadicalWorkBudget scope(satBudget_);
    SATSolverWrapper solver(Config::SolverType::CADICAL);
    FrameFormulaEncoder encoder(solver, {}, true);
    solver.addClause({encoder.encode(formula)});
    const auto answer = solver.solveWithResourceLimits(
        options_.maxSatConflicts, options_.maxSatDecisions);
    if (answer == SATSolverWrapper::SolveStatus::Unknown) {
      resource("Symbolic settling SAT obligation returned UNKNOWN under its work budget");
    }
    return answer;
  }

  void requireUnsat(BoolExpr* failure, const std::string& obligation,
                    CertificationStatus status = CertificationStatus::UnprovedBound) {
    if (solve(failure) != SATSolverWrapper::SolveStatus::Unsat) {
      throw Failure{status, obligation +
          "; certificate not established (candidate boundary states need not be reachable)"};
    }
  }

  void advance(SymbolicState& state, size_t waves, bool canonical) {
    dag_.inspect(state);
    for (size_t wave = 0; wave < waves; ++wave) {
      state = model_.wave(state, canonical
          ? FreshSymbol([] { return BoolExpr::createFalse(); })
          : FreshSymbol([this] { return fresh(); }));
      dag_.inspect(state);
    }
  }

  size_t settle(SymbolicState& state, BoolExpr* assumption, const std::string& obligation) {
    size_t depth = 0;
    dag_.inspect(state);
    for (;;) {
      auto* failure = BoolExpr::And(assumption, BoolExpr::Not(model_.stable(state)));
      if (solve(failure) == SATSolverWrapper::SolveStatus::Unsat) return depth;
      if (depth == options_.maxWaves) {
        unproved(obligation + " unproved through maxWaves; a larger bound or stronger "
            "inductive invariant may be needed, or behavior may not settle");
      }
      const size_t next = depth == 0 ? 1 :
          depth > options_.maxWaves / 2 ? options_.maxWaves : depth * 2;
      advance(state, next - depth, false);
      depth = next;
    }
  }

  const SymbolicNetwork& network_;
  const SymbolicCompileOptions& options_;
  SymbolicEventModel model_;
  size_t nextSymbol_ = 2;
  DagBudget dag_;
  SATSolverWrapper::CadicalWorkBudget satBudget_;
};

}  // namespace

SymbolicCompileResult compileSymbolicNetwork(const SymbolicNetwork& network,
                                            const SymbolicCompileOptions& options) {
  SymbolicCompileResult result;
  try {
    result.model = Compiler(network, options).run();
    result.status = CertificationStatus::Certified;
  } catch (const Failure& failure) {
    result.status = failure.status;
    result.detail = failure.detail;
  } catch (const Limit& limit) {
    result.status = CertificationStatus::ResourceLimit;
    result.detail = limit.what();
  } catch (const std::bad_alloc&) {
    result.status = CertificationStatus::ResourceLimit;
    result.detail = "Memory exhausted during symbolic settling certification";
  } catch (const std::exception& error) {
    result.status = CertificationStatus::Invalid;
    result.detail = error.what();
  }
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::LATCH
