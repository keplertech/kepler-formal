// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "../../config/Config.h"
#include "kinduction/KInductionProblem.h"

namespace KEPLER_FORMAL::SEC {

// The set F of van Eijk (TCAD 2000, section III): the signals of both designs
// among which correspondences are looked for. Any set is sound; a larger one
// proves more and costs more.
enum class SignalCorrespondenceScope {
  // Every register, every compared output and every register's next-state
  // function.
  RegistersAndOutputs,
  // Every node of the next-state and output logic.
  AllSignals,
};

// Engineering limits of the pass, not techniques from the papers. It encodes
// two time frames of both designs, so it is meant for the designs the engines
// struggle with, not for the very large ones the register learner already
// handles structurally. The defaults are the limits of a normal run; when one
// is reached the pass gives up and reports nothing.
struct SignalCorrespondenceLimits {
  size_t stateSymbols = size_t{1} << 14;
  size_t nodes = size_t{1} << 21;
  // Above this many nodes the all-signals scope is not tried.
  size_t allSignalsNodes = size_t{1} << 18;
  size_t rounds = 64;
  // The deepest induction the learner tries: 1, 2, 4, ... time frames.
  size_t inductionFrames = 8;
  // The same on the binary functions, whose logic is half the size. Sixteen
  // frames span two machine cycles of eight clock phases, the reach of a
  // value one instruction stores and the next one reads.
  size_t binaryInductionFrames = 16;
  // Classes asked in one query.
  size_t groupClasses = 64;
  // One query about a group of classes, and one about a single class. A query
  // is limited by its conflicts and decisions only: propagation alone costs
  // ticks in proportion to the design, so a fixed tick limit would leave
  // every query on a larger design undecided.
  int64_t groupConflicts = 20000;
  int64_t groupDecisions = 2000000;
  int64_t classConflicts = 10000;
  int64_t classDecisions = 1000000;
  // The pass as a whole. Each query is small; this only stops a design the
  // pass does not suit from holding up the engines for more than a couple of
  // minutes. Ticks are the limit that binds: they follow the time spent.
  uint64_t budgetConflicts = 1000000;
  uint64_t budgetDecisions = 1000000000;
  uint64_t budgetTicks = 6000000000ULL;
  // The ticks of the runs on the binary functions, which go deeper.
  uint64_t binaryBudgetTicks = 12000000000ULL;
};

// What ended a run that did not complete.
enum class SignalCorrespondenceStop {
  None,
  // The model is over a size limit or holds an operator the pass does not
  // read.
  Model,
  Rounds,
  BudgetConflicts,
  BudgetDecisions,
  BudgetTicks,
};

// The binary functions of both designs, for a run in binary semantics.
struct BinarySignalCorrespondenceInput {
  // A register, by the symbol the functions read it by, and its next-state
  // function.
  std::vector<std::pair<size_t, BoolExpr*>> transitions;
  // Per observed output, the function of each design; null when a design has
  // none.
  std::vector<std::pair<BoolExpr*, BoolExpr*>> outputs;
  // Register pairs that hold one value at boot when both are unknown there.
  std::vector<std::pair<size_t, size_t>> commonStart;
};

struct SignalCorrespondenceResult {
  // False when the pass gave up (size or solver budget, or an operator it does
  // not read); nothing is reported then.
  bool completed = false;
  // The limit that ended a run that did not complete.
  SignalCorrespondenceStop stoppedBy = SignalCorrespondenceStop::None;
  // State symbol pairs that are equal in every frame, the boot state and the
  // reset window included, so they can be handed to the engines like the
  // pairs of the register learner.
  std::vector<std::pair<size_t, size_t>> stateEqualities;
  // One flag per observed output: its two cones have the same value in every
  // frame from the first compared one on.
  std::vector<bool> provedOutputs;
  // One flag per observed output: its two cones take different values in a
  // state the start state reaches, in simulation or in the solver's check of
  // the first frames. No larger scope and no deeper induction proves such an
  // output, so a caller need not try again for it.
  std::vector<bool> differingOutputs;
  size_t signals = 0;
  size_t classes = 0;
  size_t rounds = 0;
  size_t solverCalls = 0;
  // Solver work spent, for a caller that shares one budget between runs.
  uint64_t conflicts = 0;
  uint64_t decisions = 0;
  uint64_t ticks = 0;
};

// Signal correspondence (C. A. J. van Eijk, "Sequential equivalence checking
// based on structural similarities", IEEE TCAD 19(7), 2000), decided with a
// SAT solver instead of BDDs (Bjesse and Claessen, "SAT-based verification
// without state space traversal", FMCAD 2000).
//
// Signals are grouped by random sequential simulation from the start state,
// each normalised in polarity so that signals holding opposite values share a
// class. A class is kept when its members are equal in the start state for
// every input (eq. 4) and when, with every class assumed in one time frame,
// its members are equal in the next (eq. 5); a counterexample splits the
// classes and the step is repeated until nothing splits. What is left is a
// signal correspondence relation, so its members are equal in every reachable
// state.
//
// The start state is the boot state, or with a reset bootstrap the state the
// reset window leaves, computed by propagating constants through the window
// with every other input unknown. After the window the reset ports are held
// at their inactive value, as in the engines.
//
// Assumed equalities are plain constraints of the current frame. Nothing is
// merged before it is proved.
//
// With `inductionFrames` above 1 the induction is over that many time frames
// (Sheeran, Singh and Stalmarck, FMCAD 2000; Bjesse and Claessen combine it
// with the method above): the classes must hold in the first frames from the
// start state, and, assumed in that many consecutive frames, in the next. A
// relation that only closes over a whole multi-cycle operation, a value read
// in one phase and written back in a later one, needs it.
//
// In the dual-rail encoding a class may relate signals that are unknown, X on
// both sides. A caller that disallows X equality must not run the pass.
//
// With `binary` the pass reads the binary functions of both designs instead
// of their rail-valued lifting. The lifting is a shape-dependent
// over-approximation: two shapes of one function can give X and a defined
// value for the same unknown inputs, and no class relates the two then. The
// run starts at boot, where the unknown registers of a `commonStart` pair
// hold one value (van Eijk's designated initial state; Savoj et al., FMCAD
// 2010, for the common unknown start), and goes through the reset window
// before the first frame it compares. An output proved this way has no
// binary-defined mismatch in the rail model: a defined rail value is the
// value of every binary completion of the unknown boot state (Seger and
// Bryant, Formal Methods in System Design 6(2), 1995), the run covers at
// least one completion, so two defined outputs agree. State equalities are
// not reported: they hold from the common boot states only.
SignalCorrespondenceResult proveSignalCorrespondence(
    const KInductionProblem& problem,
    SignalCorrespondenceScope scope,
    Config::SolverType solverType,
    const SignalCorrespondenceLimits& limits = {},
    size_t inductionFrames = 1,
    const BinarySignalCorrespondenceInput* binary = nullptr);

struct SignalCorrespondenceLearning {
  // One flag per observed output.
  std::vector<bool> provedOutputs;
  std::vector<std::pair<size_t, size_t>> stateEqualities;
};

// Runs the pass for the outputs no earlier step proved: on registers and
// outputs with induction over 1, 2, 4, ... time frames, then on all signals
// the same way, while outputs remain and the budget lasts. An output whose
// cones differ in a reachable state is not tried again. The problem is not
// changed; the caller discharges the proved outputs and records the state
// equalities. With `binary` the runs read the binary functions, on registers
// and outputs only, up to `binaryInductionFrames` within `binaryBudgetTicks`.
SignalCorrespondenceLearning learnSignalCorrespondence(
    const KInductionProblem& problem,
    Config::SolverType solverType,
    bool diagnostics,
    SignalCorrespondenceLimits limits = {},
    const BinarySignalCorrespondenceInput* binary = nullptr);

}  // namespace KEPLER_FORMAL::SEC
