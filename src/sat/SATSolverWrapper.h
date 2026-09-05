#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>
#include <string>
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <memory>

#include "simp/SimpSolver.h"
#include "cadical.hpp"
#include "craigtracer.hpp"
extern "C" {
  #include "kissat.h"
}

#include "../config/Config.h"

//#define USE_KISSAT

class SATSolverWrapper {
public:
  enum class SolveStatus {
    Sat,
    Unsat,
    Unknown,
  };

  // One PDR output can use several incremental CaDiCaL instances.  Share this
  // counter across them so many individually bounded queries cannot make the
  // output's total SAT work unbounded.
  class CadicalWorkBudget {
   public:
    CadicalWorkBudget(uint64_t conflictLimit,
                      uint64_t decisionLimit,
                      uint64_t tickLimit)
        : conflictLimit_(conflictLimit),
          decisionLimit_(decisionLimit),
          tickLimit_(tickLimit) {}

    uint64_t conflictLimit() const { return conflictLimit_; }
    uint64_t decisionLimit() const { return decisionLimit_; }
    uint64_t tickLimit() const { return tickLimit_; }
    uint64_t conflictsUsed() const { return conflictsUsed_; }
    uint64_t decisionsUsed() const { return decisionsUsed_; }
    uint64_t ticksUsed() const { return ticksUsed_; }

    uint64_t remainingConflicts() const {
      return conflictsUsed_ >= conflictLimit_
                 ? 0
                 : conflictLimit_ - conflictsUsed_;
    }

    uint64_t remainingDecisions() const {
      return decisionsUsed_ >= decisionLimit_
                 ? 0
                 : decisionLimit_ - decisionsUsed_;
    }

    uint64_t remainingTicks() const {
      return ticksUsed_ >= tickLimit_ ? 0 : tickLimit_ - ticksUsed_;
    }

    bool exhausted() const {
      return remainingConflicts() == 0 ||
             remainingDecisions() == 0 ||
             remainingTicks() == 0;
    }

   private:
    void consume(uint64_t conflicts, uint64_t decisions, uint64_t ticks) {
      conflictsUsed_ += std::min(conflicts, remainingConflicts());
      decisionsUsed_ += std::min(decisions, remainingDecisions());
      ticksUsed_ += std::min(ticks, remainingTicks());
    }

    uint64_t conflictLimit_ = 0;
    uint64_t decisionLimit_ = 0;
    uint64_t tickLimit_ = 0;
    uint64_t conflictsUsed_ = 0;
    uint64_t decisionsUsed_ = 0;
    uint64_t ticksUsed_ = 0;

    friend class SATSolverWrapper;
  };

  // PDR is currently serial per output.  A thread-local scope lets every
  // property-local CaDiCaL solver charge the same output budget without
  // coupling the generic SAT API to the SEC problem classes.
  class ScopedCadicalWorkBudget {
   public:
    explicit ScopedCadicalWorkBudget(CadicalWorkBudget& budget)
        : previous_(activeCadicalWorkBudget_) {
      activeCadicalWorkBudget_ = &budget;
    }

    ~ScopedCadicalWorkBudget() {
      activeCadicalWorkBudget_ = previous_;
    }

    ScopedCadicalWorkBudget(const ScopedCadicalWorkBudget&) = delete;
    ScopedCadicalWorkBudget& operator=(const ScopedCadicalWorkBudget&) = delete;

   private:
    CadicalWorkBudget* previous_ = nullptr;
  };

  enum class CraigVariablePartition {
    ALocal,
    BLocal,
    Global,
  };

  enum class CraigClausePartition {
    A,
    B,
  };

  struct CraigInterpolantCnf { // LCOV_EXCL_LINE
    enum class Type {
      None,
      ConstantFalse,
      ConstantTrue,
      Normal,
    };

    Type type = Type::None; // LCOV_EXCL_LINE
    int firstAuxiliaryVariable = 0; // LCOV_EXCL_LINE
    std::vector<std::vector<int>> clauses;
  };

  static KEPLER_FORMAL::Config::SolverType assumptionSolverTypeFor(
      KEPLER_FORMAL::Config::SolverType chosenSolverType) {
    return chosenSolverType == KEPLER_FORMAL::Config::SolverType::GLUCOSE
               ? KEPLER_FORMAL::Config::SolverType::GLUCOSE
               : KEPLER_FORMAL::Config::SolverType::CADICAL;
  }

  explicit SATSolverWrapper(KEPLER_FORMAL::Config::SolverType type = KEPLER_FORMAL::Config::SolverType::CADICAL)
    : solverType_(type) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      glucoseSolver_ = std::make_unique<Glucose::SimpSolver>();
      // Keep embedded solver calls quiet. The command-line Glucose default is
      // already low, but SimpSolver still prints eliminated-clause summaries at
      // verbosity 0.
      glucoseSolver_->verbosity = -1;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      kissatSolver_ = kissat_init();
      // Embedded SEC runs can create thousands of short-lived Kissat queries.
      // Keep those solver instances quiet so regressions are not dominated by
      // progress-report I/O.  Some embedded Kissat builds do not expose
      // application-only options such as "quiet", so this must remain
      // best-effort instead of aborting the proof setup.
      setKissatOptionIfSupported(static_cast<kissat*>(kissatSolver_), "quiet", 1);
      kissatNumVars_ = 0;
      kissatReservedVars_ = 0;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      cadicalSolver_ = std::make_unique<CaDiCaL::Solver>(); // LCOV_EXCL_LINE
      cadicalNumVars_ = 0; // LCOV_EXCL_LINE
      cadicalReservedVars_ = 0; // LCOV_EXCL_LINE
      if (CaDiCaL::Solver::is_valid_option("quiet")) { // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        cadicalSolver_->set("quiet", 1);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    } else { // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      throw std::invalid_argument("Unknown solver type");
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
  }

  ~SATSolverWrapper() {
    if (cadicalCraigTracer_ != nullptr && cadicalSolver_ != nullptr) {
      cadicalSolver_->disconnect_proof_tracer(cadicalCraigTracer_.get()); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT && kissatSolver_) {
      kissat_release(static_cast<kissat*>(kissatSolver_));
    }
  }

  void enableCraigInterpolation() { // LCOV_EXCL_LINE
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::CADICAL) { // LCOV_EXCL_LINE
      throw std::runtime_error( // LCOV_EXCL_LINE
          "Craig interpolation requires the CaDiCaL solver backend");
    }
    if (cadicalCraigTracer_ != nullptr) { // LCOV_EXCL_LINE
      return; // LCOV_EXCL_LINE
    }
    // CaDiCaL's Craig tracer requires proof-preserving CNF construction. BVA
    // introduces variables without a caller-visible A/B/global partition.
    cadicalSolver_->set("factor", 0); // LCOV_EXCL_LINE
    cadicalCraigTracer_ = std::make_unique<CaDiCraig::CraigTracer>(); // LCOV_EXCL_LINE
    cadicalSolver_->connect_proof_tracer(cadicalCraigTracer_.get(), true); // LCOV_EXCL_LINE
    cadicalCraigTracer_->set_craig_construction( // LCOV_EXCL_LINE
        CaDiCraig::CraigConstruction::ASYMMETRIC);
  } // LCOV_EXCL_LINE

  void setCraigVariablePartition(CraigVariablePartition partition) { // LCOV_EXCL_LINE
    craigVariablePartition_ = partition; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  void setCraigClausePartition(CraigClausePartition partition) { // LCOV_EXCL_LINE
    craigClausePartition_ = partition; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  CraigInterpolantCnf createCraigInterpolant() { // LCOV_EXCL_LINE
    if (cadicalCraigTracer_ == nullptr) { // LCOV_EXCL_LINE
      throw std::runtime_error("Craig interpolation was not enabled"); // LCOV_EXCL_LINE
    }
    std::vector<std::vector<int>> cadicalClauses; // LCOV_EXCL_LINE
    int nextVariable = cadicalNumVars_ + 1; // LCOV_EXCL_LINE
    const auto result = cadicalCraigTracer_->create_craig_interpolant( // LCOV_EXCL_LINE
        CaDiCraig::CraigInterpolant::ASYMMETRIC,
        cadicalClauses,
        nextVariable);

    CraigInterpolantCnf interpolant; // LCOV_EXCL_LINE
    interpolant.firstAuxiliaryVariable = cadicalNumVars_ + 2; // LCOV_EXCL_LINE
    switch (result) { // LCOV_EXCL_LINE
      case CaDiCraig::CraigCnfType::NONE:
        interpolant.type = CraigInterpolantCnf::Type::None; // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
      case CaDiCraig::CraigCnfType::CONSTANT0:
        interpolant.type = CraigInterpolantCnf::Type::ConstantFalse; // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
      case CaDiCraig::CraigCnfType::CONSTANT1:
        interpolant.type = CraigInterpolantCnf::Type::ConstantTrue; // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
      case CaDiCraig::CraigCnfType::NORMAL:
        interpolant.type = CraigInterpolantCnf::Type::Normal; // LCOV_EXCL_LINE
        break; // LCOV_EXCL_LINE
    }

    interpolant.clauses.reserve(cadicalClauses.size()); // LCOV_EXCL_LINE
    for (const auto& clause : cadicalClauses) { // LCOV_EXCL_LINE
      std::vector<int> externalClause; // LCOV_EXCL_LINE
      externalClause.reserve(clause.size()); // LCOV_EXCL_LINE
      for (const int literal : clause) { // LCOV_EXCL_LINE
        const int externalVariable = std::abs(literal) + 1; // LCOV_EXCL_LINE
        externalClause.push_back( // LCOV_EXCL_LINE
            literal > 0 ? externalVariable : -externalVariable); // LCOV_EXCL_LINE
      }
      interpolant.clauses.push_back(std::move(externalClause)); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
    return interpolant; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  // Create a new variable (returns 0-based index)
  int newVar() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      return glucoseSolver_->newVar();
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // Kissat does not require explicit variable creation, but SEC encoders
      // create variables before streaming clauses. Reserving at creation time
      // keeps Kissat from growing its internal variable arrays inside
      // kissat_add(), which is the hottest path for large PDR queries.
      reserveVars(static_cast<size_t>(kissatNumVars_) + 1);
      return kissatNumVars_++;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) { // LCOV_EXCL_LINE
      reserveVars(static_cast<size_t>(cadicalNumVars_) + 1); // LCOV_EXCL_LINE
      const int variable = cadicalNumVars_++; // LCOV_EXCL_LINE
      if (cadicalCraigTracer_ != nullptr) { // LCOV_EXCL_LINE
        cadicalCraigTracer_->label_variable( // LCOV_EXCL_LINE
            variable + 1, cadicalCraigVariableType(craigVariablePartition_)); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
      return variable; // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }

  void reserveVars(size_t numVars) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      if (numVars <= static_cast<size_t>(kissatReservedVars_)) {
        return;
      }
      const size_t currentReserved = static_cast<size_t>(kissatReservedVars_);
      const size_t growthSlack =
          std::max(currentReserved / 2, static_cast<size_t>(4096));
      const size_t targetVars =
          currentReserved == 0
              ? numVars
              : std::max(numVars, currentReserved + growthSlack);
      // Kissat allocates variable storage lazily while clauses are added. PDR
      // creates many short-lived solvers with thousands of frame variables, so
      // reserving the known frame-variable prefix avoids repeated internal
      // growth during Tseitin clause emission.
      kissat_reserve(static_cast<kissat*>(kissatSolver_),
                     static_cast<int>(targetVars));
      kissatReservedVars_ = static_cast<int>(targetVars);
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      if (numVars <= static_cast<size_t>(cadicalReservedVars_)) { // LCOV_EXCL_LINE
        return; // LCOV_EXCL_LINE
      }
      const size_t currentReserved = static_cast<size_t>(cadicalReservedVars_); // LCOV_EXCL_LINE
      const size_t growthSlack = // LCOV_EXCL_LINE
          std::max(currentReserved / 2, static_cast<size_t>(4096)); // LCOV_EXCL_LINE
      const size_t targetVars = // LCOV_EXCL_LINE
          currentReserved == 0 // LCOV_EXCL_LINE
              ? numVars // LCOV_EXCL_LINE
              : std::max(numVars, currentReserved + growthSlack); // LCOV_EXCL_LINE
      if (targetVars > static_cast<size_t>(std::numeric_limits<int>::max())) { // LCOV_EXCL_LINE
        // LCOV_EXCL_START
        throw std::runtime_error("CaDiCaL variable reservation exceeds int range");  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      cadicalSolver_->resize(static_cast<int>(targetVars)); // LCOV_EXCL_LINE
      cadicalReservedVars_ = static_cast<int>(targetVars); // LCOV_EXCL_LINE
    } // LCOV_EXCL_LINE
  }

  void reserveAdditionalVars(size_t additionalVars) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // Formula encoders create Tseitin variables after the frame variables
      // have already been allocated. Reserving from the wrapper's current
      // variable count targets Kissat's variable-vector growth directly without
      // changing the external literal numbering convention.
      reserveVars(static_cast<size_t>(kissatNumVars_) + additionalVars);
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      reserveVars(static_cast<size_t>(cadicalNumVars_) + additionalVars);
    }
  }

  // Add a clause, literals are signed ints:
  // external convention: 0=false, 1=true, vars are ±(var_id + 2)
  void addClause(const std::vector<int>& lits) {
    addClauseRange(lits);
  }

  void addClause(std::initializer_list<int> lits) {
    addClauseRange(lits);
  }

 private:
  template <typename ClauseRange>
  void addClauseRange(const ClauseRange& lits) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      Glucose::vec<Glucose::Lit> clause;
      for (int lit : lits) {
        if (lit == 0 || lit == 1) {
          // We should never see raw consts here: they are encoded via forced vars.
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Constant literal (0/1) passed to Glucose clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        int v = std::abs(lit);
        int var = v - 2;  // external ±(var+2) -> internal var index
        if (var < 0) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Invalid literal (<2) passed to Glucose clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        while (var >= glucoseSolver_->nVars())
          // LCOV_EXCL_START
          glucoseSolver_->newVar();  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        clause.push((lit > 0) ? Glucose::mkLit(var) : ~Glucose::mkLit(var));
      }
      glucoseSolver_->addClause(clause);
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // Kissat expects ±(var+1), 0 terminates a clause.
      for (int lit : lits) {
        if (lit == 0 || lit == 1) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Constant literal (0/1) passed to Kissat clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        int v = std::abs(lit);
        int var = v - 2;  // external ±(var+2) -> internal var index
        if (var < 0) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Invalid literal (<2) passed to Kissat clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        if (var >= kissatNumVars_) {
          // LCOV_EXCL_START
          kissatNumVars_ = var + 1;
        }
        // LCOV_EXCL_STOP
        int kissatLit = (lit > 0 ? var + 1 : -(var + 1)); // ±(var+1)
        kissat_add(static_cast<kissat*>(kissatSolver_), kissatLit);
      }
      kissat_add(static_cast<kissat*>(kissatSolver_), 0); // end of clause
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      if (cadicalCraigTracer_ != nullptr) { // LCOV_EXCL_LINE
        cadicalCraigTracer_->label_clause( // LCOV_EXCL_LINE
            ++cadicalCraigClauseId_, // LCOV_EXCL_LINE
            craigClausePartition_ == CraigClausePartition::A // LCOV_EXCL_LINE
                ? CaDiCraig::CraigClauseType::A_CLAUSE
                : CaDiCraig::CraigClauseType::B_CLAUSE); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
      for (int lit : lits) {
        if (lit == 0 || lit == 1) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Constant literal (0/1) passed to CaDiCaL clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        int v = std::abs(lit);
        int var = v - 2;  // external ±(var+2) -> internal var index
        if (var < 0) {
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Invalid literal (<2) passed to CaDiCaL clause");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        if (var >= cadicalNumVars_) {
          // LCOV_EXCL_START
          cadicalNumVars_ = var + 1;  // LCOV_EXCL_LINE
          reserveVars(static_cast<size_t>(cadicalNumVars_));  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        int cadicalLit = (lit > 0 ? var + 1 : -(var + 1)); // ±(var+1)
        cadicalSolver_->add(cadicalLit);
      }
      cadicalSolver_->add(0); // end of clause
    } else {
      // LCOV_EXCL_START
      // LCOV_DISABLED_START
      throw std::runtime_error("Unknown solver type");
      // LCOV_DISABLED_STOP
      // LCOV_EXCL_STOP
    }
  }

 public:
  bool solve() {
    return solveStatus() == SolveStatus::Sat;
  }

  SolveStatus solveStatus() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      return glucoseSolver_->solve() ? SolveStatus::Sat : SolveStatus::Unsat;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      int res = kissat_solve(static_cast<kissat*>(kissatSolver_));
      if (res == 10) { // 10 = SAT
        return SolveStatus::Sat;
      }
      if (res == 20) { // 20 = UNSAT
        return SolveStatus::Unsat;
      }
      // LCOV_EXCL_START
      return SolveStatus::Unknown;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {  // LCOV_EXCL_LINE
      lastAssumptionSolveStatus_ = SolveStatus::Unknown;  // LCOV_EXCL_LINE
      lastAssumptions_.clear();  // LCOV_EXCL_LINE
      return solveCadicalWithResourceLimits(
          /*conflictLimit=*/-1,
          /*decisionLimit=*/-1,
          /*tickLimit=*/-1,
          /*useCumulativeBudget=*/false);
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  }

  SolveStatus solveWithKissatResourceLimits(
      unsigned conflictLimit,
      unsigned decisionLimit =
          std::numeric_limits<unsigned>::max()) {
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // LCOV_EXCL_START
      throw std::runtime_error(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          "Kissat resource limits requested for non-Kissat solver");
    }
    auto* solver = static_cast<kissat*>(kissatSolver_);
    if (conflictLimit != std::numeric_limits<unsigned>::max()) {
      kissat_set_conflict_limit(solver, conflictLimit);
    }
    if (decisionLimit != std::numeric_limits<unsigned>::max()) {
      kissat_set_decision_limit(solver, decisionLimit);
    }
    return solveStatus();
  // LCOV_EXCL_START
  }  // LCOV_EXCL_LINE
  // LCOV_EXCL_STOP

  // LCOV_EXCL_START
  SolveStatus solveWithResourceLimits(
  // LCOV_EXCL_STOP
      unsigned conflictLimit,
      unsigned decisionLimit =
          std::numeric_limits<unsigned>::max()) {
    // LCOV_EXCL_START
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      return solveWithKissatResourceLimits(conflictLimit, decisionLimit);
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      lastAssumptionSolveStatus_ = SolveStatus::Unknown;
      lastAssumptions_.clear();
      return solveCadicalWithResourceLimits(
          conflictLimit == std::numeric_limits<unsigned>::max()
              ? -1
              : static_cast<int64_t>(conflictLimit),
          decisionLimit == std::numeric_limits<unsigned>::max()
              ? -1
              : static_cast<int64_t>(decisionLimit),
          /*tickLimit=*/-1,
          /*useCumulativeBudget=*/true);
      // LCOV_EXCL_STOP
    }
    // LCOV_EXCL_START
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {  // LCOV_EXCL_LINE
      Glucose::vec<Glucose::Lit> noAssumptions;  // LCOV_EXCL_LINE
      if (conflictLimit != std::numeric_limits<unsigned>::max()) {  // LCOV_EXCL_LINE
        glucoseSolver_->setConfBudget(conflictLimit);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      if (decisionLimit != std::numeric_limits<unsigned>::max()) {  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
        // Glucose has no decision budget here; a propagation budget is the
        // closest local limiter and preserves UNKNOWN as non-proof.
        // LCOV_EXCL_START
        glucoseSolver_->setPropBudget(decisionLimit);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      const auto result = glucoseSolver_->solveLimited(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
          noAssumptions,
          /*do_simp=*/false,
          /*turn_off_simp=*/true);
      // LCOV_EXCL_START
      glucoseSolver_->budgetOff();  // LCOV_EXCL_LINE
      if (Glucose::toInt(result) == 0) {  // LCOV_EXCL_LINE
        return SolveStatus::Sat;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (Glucose::toInt(result) == 1) {  // LCOV_EXCL_LINE
        return SolveStatus::Unsat;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      return SolveStatus::Unknown;  // LCOV_EXCL_LINE
    }  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  // LCOV_EXCL_START
  }
  // LCOV_EXCL_STOP

  SolveStatus solveWithAssumptionsStatus( // LCOV_EXCL_LINE
      const std::vector<int>& assumptions,
      int64_t conflictLimit = -1,
      int64_t propagationLimit = -1,
      int64_t tickLimit = -1) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) { // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      Glucose::vec<Glucose::Lit> glucoseAssumptions;  // LCOV_EXCL_LINE
      for (int lit : assumptions) {  // LCOV_EXCL_LINE
        if (lit == 0 || lit == 1) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Constant literal (0/1) passed as Glucose assumption");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        const int var = std::abs(lit) - 2;  // LCOV_EXCL_LINE
        if (var < 0) {  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Invalid literal (<2) passed as Glucose assumption");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        while (var >= glucoseSolver_->nVars()) {  // LCOV_EXCL_LINE
          glucoseSolver_->newVar();  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        glucoseAssumptions.push((lit > 0) ? Glucose::mkLit(var)  // LCOV_EXCL_LINE
                                          : ~Glucose::mkLit(var));  // LCOV_EXCL_LINE
                                          // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      if (conflictLimit >= 0 || propagationLimit >= 0) {  // LCOV_EXCL_LINE
        if (conflictLimit >= 0) {  // LCOV_EXCL_LINE
          glucoseSolver_->setConfBudget(conflictLimit);  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        if (propagationLimit >= 0) {  // LCOV_EXCL_LINE
          glucoseSolver_->setPropBudget(propagationLimit);  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        const auto result = glucoseSolver_->solveLimited(  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
            glucoseAssumptions,
            /*do_simp=*/false,
            /*turn_off_simp=*/true);
        // LCOV_EXCL_START
        glucoseSolver_->budgetOff();  // LCOV_EXCL_LINE
        if (Glucose::toInt(result) == 0) {  // LCOV_EXCL_LINE
          return SolveStatus::Sat;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        if (Glucose::toInt(result) == 1) {  // LCOV_EXCL_LINE
          return SolveStatus::Unsat;  // LCOV_EXCL_LINE
          // LCOV_EXCL_STOP
        }
        // LCOV_EXCL_START
        return SolveStatus::Unknown;  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // Repeated CEGAR reachability checks reuse the same solver and vary only
      // assumptions. Running variable elimination on each assumption solve
      // dominates those small queries, so keep this path in plain CDCL mode.
      // LCOV_EXCL_START
      return glucoseSolver_->solve(  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
                 glucoseAssumptions,
                 /*do_simp=*/false,
                 /*turn_off_simp=*/true)
                 ? SolveStatus::Sat
                 : SolveStatus::Unsat;
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) { // LCOV_EXCL_LINE
      // This vendored Kissat exposes only the partial IPASIR API and has no
      // assumption call. Callers that need repeated assumption solves should use
      // CaDiCaL for that local incremental query.
      // LCOV_EXCL_START
      if (assumptions.empty()) {  // LCOV_EXCL_LINE
        return solveStatus();  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
      }
      // LCOV_EXCL_START
      throw std::runtime_error("Kissat assumptions are not available in this build");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) { // LCOV_EXCL_LINE
      const bool useCumulativeBudget =
          conflictLimit >= 0 || propagationLimit >= 0 || tickLimit >= 0;
      lastAssumptions_ = assumptions; // LCOV_EXCL_LINE
      lastAssumptionSolveStatus_ = SolveStatus::Unknown; // LCOV_EXCL_LINE
      // Do not enqueue assumptions if this output has no SAT work left.
      // CaDiCaL consumes assumptions only on solve(), so an early return after
      // assume() would otherwise leak them into the next output's shared solver.
      if (useCumulativeBudget && activeCadicalWorkBudget_ != nullptr &&
          activeCadicalWorkBudget_->exhausted()) {
        lastAssumptions_.clear();
        return lastAssumptionSolveStatus_;
      }
      for (int lit : assumptions) { // LCOV_EXCL_LINE
        if (lit == 0 || lit == 1) { // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Constant literal (0/1) passed as CaDiCaL assumption");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        const int var = std::abs(lit) - 2; // LCOV_EXCL_LINE
        if (var < 0) { // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          // LCOV_DISABLED_START
          throw std::runtime_error("Invalid literal (<2) passed as CaDiCaL assumption");
          // LCOV_DISABLED_STOP
          // LCOV_EXCL_STOP
        }
        if (var >= cadicalNumVars_) { // LCOV_EXCL_LINE
          // LCOV_EXCL_START
          cadicalNumVars_ = var + 1;  // LCOV_EXCL_LINE
          reserveVars(static_cast<size_t>(cadicalNumVars_));  // LCOV_EXCL_LINE
        }  // LCOV_EXCL_LINE
        // LCOV_EXCL_STOP
        cadicalSolver_->assume(lit > 0 ? var + 1 : -(var + 1)); // LCOV_EXCL_LINE
      }
      lastAssumptionSolveStatus_ = solveCadicalWithResourceLimits(
          conflictLimit,
          // CaDiCaL has no propagation budget. Use its decision budget as the
          // closest deterministic limiter for callers that pass both.
          propagationLimit,
          tickLimit,
          useCumulativeBudget);
      if (lastAssumptionSolveStatus_ != SolveStatus::Unsat) { // LCOV_EXCL_LINE
        lastAssumptions_.clear(); // LCOV_EXCL_LINE
      } // LCOV_EXCL_LINE
      return lastAssumptionSolveStatus_; // LCOV_EXCL_LINE
    }
    // LCOV_EXCL_START
    // LCOV_DISABLED_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_DISABLED_STOP
    // LCOV_EXCL_STOP
  } // LCOV_EXCL_LINE

  bool solveWithAssumptions(const std::vector<int>& assumptions) {
    return solveWithAssumptionsStatus(assumptions) == SolveStatus::Sat;
  }

  std::vector<int> failedAssumptions() const { // LCOV_EXCL_LINE
    std::vector<int> failed; // LCOV_EXCL_LINE
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) { // LCOV_EXCL_LINE
      // LCOV_EXCL_START
      failed.reserve(glucoseSolver_->conflict.size());  // LCOV_EXCL_LINE
      for (int i = 0; i < glucoseSolver_->conflict.size(); ++i) {  // LCOV_EXCL_LINE
        const auto lit = glucoseSolver_->conflict[i];  // LCOV_EXCL_LINE
        const int externalVar = Glucose::var(lit) + 2;  // LCOV_EXCL_LINE
        const int conflictLiteral =  // LCOV_EXCL_LINE
            Glucose::sign(lit) ? -externalVar : externalVar;  // LCOV_EXCL_LINE
            // LCOV_EXCL_STOP
        // Glucose exposes the final conflict clause over the negated failed
        // assumptions. Return the caller's original assumption literals so SEC
        // can map them directly back to PDR cube literals.
        // LCOV_EXCL_START
        failed.push_back(-conflictLiteral);  // LCOV_EXCL_LINE
      }  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL && // LCOV_EXCL_LINE
               lastAssumptionSolveStatus_ == SolveStatus::Unsat) { // LCOV_EXCL_LINE
      failed.reserve(lastAssumptions_.size()); // LCOV_EXCL_LINE
      for (int lit : lastAssumptions_) { // LCOV_EXCL_LINE
        const int var = std::abs(lit) - 2; // LCOV_EXCL_LINE
        const int cadicalLit = lit > 0 ? var + 1 : -(var + 1); // LCOV_EXCL_LINE
        if (cadicalSolver_->failed(cadicalLit)) { // LCOV_EXCL_LINE
          failed.push_back(lit); // LCOV_EXCL_LINE
        } // LCOV_EXCL_LINE
      }
    } // LCOV_EXCL_LINE
    return failed; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  bool getLiteralValue(int lit) const {
    if (lit == 0) {
      // LCOV_EXCL_START
      return false;
      // LCOV_EXCL_STOP
    }
    if (lit == 1) {
      // LCOV_EXCL_START
      return true;
      // LCOV_EXCL_STOP
    }

    const int external = std::abs(lit);
    const int var = external - 2;
    if (var < 0) {
      // LCOV_EXCL_START
      throw std::runtime_error("Invalid literal passed to getLiteralValue");
      // LCOV_EXCL_STOP
    }

    bool positiveValue = false;
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      // LCOV_EXCL_START
      const auto value = glucoseSolver_->modelValue(Glucose::mkLit(var));
      if (Glucose::toInt(value) == 2) {
        positiveValue = false;  // LCOV_EXCL_LINE
      } else {  // LCOV_EXCL_LINE
        positiveValue = Glucose::toInt(value) == 0;
        // LCOV_EXCL_STOP
      }
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      const int value = kissat_value(static_cast<kissat*>(kissatSolver_), var + 1);
      if (value == 0) {
        positiveValue = false;
      } else {
        positiveValue = value > 0;
      }
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      // LCOV_EXCL_START
      const int value = cadicalSolver_->val(var + 1);  // LCOV_EXCL_LINE
      positiveValue = value > 0;  // LCOV_EXCL_LINE
    } else {  // LCOV_EXCL_LINE
      throw std::runtime_error("Unknown solver type");  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    return lit > 0 ? positiveValue : !positiveValue;
  }

  void configureForSecConeProof(size_t coneSymbols = 0) {
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // LCOV_EXCL_START
      return;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    auto* solver = static_cast<kissat*>(kissatSolver_);
    // SEC proof obligations are normally expected to be UNSAT when the designs
    // are equivalent, so bias Kissat toward its UNSAT/stable search profile.
    // Keep the rest of the default solver machinery for small and medium COIs:
    // on those instances preprocessing and clause minimization can be the real
    // speedup, and disabling them globally makes BlackParrot-style proofs worse.
    setKissatOptionOrThrow(solver, "stable", 2);

    // The threshold is intentionally based on the engine's COI size, not on a
    // design name. Medium SEC cones still benefit from Kissat's normal
    // congruence/probing passes because those passes recognize duplicated
    // gate-level structure; asap7_jpeg_lvt regressed badly when this cutoff was
    // lowered into the medium-cone range. Reserve the stripped-down profile for
    // truly huge cones where the prepasses themselves dominate memory and
    // runtime.
    constexpr size_t kLargeSecConeSymbolThreshold = 100000;
    if (coneSymbols < kLargeSecConeSymbolThreshold) {
      return;
    }

    // Large SEC cones are already sliced and Tseitin-encoded by the engine.
    // Profiles on very large ASIC obligations showed Kissat spending most of
    // its time in speculative preprocessing before reaching CDCL. For only
    // those large cones, skip the speculative passes and keep the query focused
    // on SAT search.
    // LCOV_EXCL_START
    setKissatOptionOrThrow(solver, "preprocess", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "simplify", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "preprocesscongruence", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "preprocessprobe", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "congruence", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "probe", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "probeinit", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "eliminate", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "eliminateinit", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "lucky", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "luckyearly", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "luckylate", 0);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
    // Recursive learned-clause shrinking can dominate these very large shallow
    // equivalence cones, so disable it only for this large-cone profile.
    // LCOV_EXCL_START
    setKissatOptionOrThrow(solver, "minimize", 0);  // LCOV_EXCL_LINE
    setKissatOptionOrThrow(solver, "shrink", 0);  // LCOV_EXCL_LINE
    // LCOV_EXCL_STOP
  }

  void configureForSecDualRailConeProof(size_t coneSymbols = 0) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      // LCOV_EXCL_START
      configureForSecPdrQuery(coneSymbols);
      return;
      // LCOV_EXCL_STOP
    }

    configureForSecConeProof(coneSymbols);
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // LCOV_EXCL_START
      return;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // Dual-rail formulas double each state bit into may-one/may-zero rails.
    // Medium cones can therefore look small by symbol count while still
    // triggering expensive Kissat probe/vivify passes.  Keep the generic SEC
    // profile unchanged, but use the direct large-cone path earlier for
    // dual-rail induction leaves.
    constexpr size_t kDualRailMediumConeSymbolThreshold = 4096;
    if (coneSymbols < kDualRailMediumConeSymbolThreshold) {
      return;
    }

    // LCOV_EXCL_START
    auto* solver = static_cast<kissat*>(kissatSolver_);
    setKissatOptionOrThrow(solver, "preprocess", 0);
    setKissatOptionOrThrow(solver, "simplify", 0);
    setKissatOptionOrThrow(solver, "preprocesscongruence", 0);
    setKissatOptionOrThrow(solver, "preprocessprobe", 0);
    setKissatOptionOrThrow(solver, "congruence", 0);
    setKissatOptionOrThrow(solver, "probe", 0);
    setKissatOptionOrThrow(solver, "probeinit", 0);
    setKissatOptionOrThrow(solver, "eliminate", 0);
    setKissatOptionOrThrow(solver, "eliminateinit", 0);
    setKissatOptionOrThrow(solver, "lucky", 0);
    setKissatOptionOrThrow(solver, "luckyearly", 0);
    setKissatOptionOrThrow(solver, "luckylate", 0);
    setKissatOptionOrThrow(solver, "minimize", 0);
    setKissatOptionOrThrow(solver, "shrink", 0);
    // LCOV_EXCL_STOP
  }

  void configureForSecPdrQuery(size_t coneSymbols = 0) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      // LCOV_EXCL_START
      auto* solver = cadicalSolver_.get();
      // LCOV_EXCL_STOP
      // CaDiCaL is the default local solver for assumption-capable validation
      // queries. Short-lived callers only need a quick SAT/UNSAT answer, so
      // avoid expensive inprocessing and recursive clause polishing that
      // samples showed dominating deeper sky130hs_ibex frontier checks.
      // LCOV_EXCL_START
      setCadicalOptionIfSupported(solver, "inprocessing", 0);
      setCadicalOptionIfSupported(solver, "compact", 0);
      setCadicalOptionIfSupported(solver, "arenacompact", 0);
      setCadicalOptionIfSupported(solver, "elim", 0);
      setCadicalOptionIfSupported(solver, "fastelim", 0);
      setCadicalOptionIfSupported(solver, "preprocesslight", 0);
      setCadicalOptionIfSupported(solver, "probe", 0);
      setCadicalOptionIfSupported(solver, "inprobing", 0);
      setCadicalOptionIfSupported(solver, "congruence", 0);
      setCadicalOptionIfSupported(solver, "decompose", 0);
      setCadicalOptionIfSupported(solver, "deduplicate", 0);
      setCadicalOptionIfSupported(solver, "factor", 0);
      setCadicalOptionIfSupported(solver, "subsume", 0);
      setCadicalOptionIfSupported(solver, "sweep", 0);
      setCadicalOptionIfSupported(solver, "lucky", 0);
      setCadicalOptionIfSupported(solver, "luckyearly", 0);
      setCadicalOptionIfSupported(solver, "luckylate", 0);
      setCadicalOptionIfSupported(solver, "minimize", 0);
      setCadicalOptionIfSupported(solver, "shrink", 0);
      setCadicalOptionIfSupported(solver, "rephase", 0);
      setCadicalOptionIfSupported(solver, "walk", 0);
      // LCOV_EXCL_STOP
      (void)coneSymbols;
      // LCOV_EXCL_START
      return;
      // LCOV_EXCL_STOP
    }

    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // LCOV_EXCL_START
      return;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    auto* solver = static_cast<kissat*>(kissatSolver_);
    // PDR issues many short-lived predecessor queries while blocking one cube
    // at a time.  These obligations are frequently SAT because PDR walks
    // backward through reachable predecessors before eventually learning an
    // UNSAT blocking clause.  Use focused/SAT-oriented CDCL here; the stable
    // UNSAT-oriented Kissat mode periodically rephases via local search, which
    // dominated ASIC PDR predecessor samples without helping those SAT walks.
    // The generic SEC cone profile below remains available for k-induction and
    // IMC, whose large monolithic queries can benefit from preprocessing.
    setKissatOptionOrThrow(solver, "stable", 0);
    setKissatOptionOrThrow(solver, "target", 2);
    setKissatOptionOrThrow(solver, "restartint", 10);
    setKissatOptionOrThrow(solver, "restartreusetrail", 0);

    setKissatOptionOrThrow(solver, "rephase", 0);
    setKissatOptionOrThrow(solver, "walkeffort", 0);
    setKissatOptionOrThrow(solver, "lucky", 0);
    setKissatOptionOrThrow(solver, "luckyearly", 0);
    setKissatOptionOrThrow(solver, "luckylate", 0);
    setKissatOptionOrThrow(solver, "minimize", 0);
    setKissatOptionOrThrow(solver, "shrink", 0);
    // PDR predecessor queries are rebuilt from engine-level learned frame
    // clauses, so Kissat's own on-the-fly clause strengthening adds heavy
    // per-conflict bookkeeping without preserving state across queries. Keep
    // this disabled only for PDR; k-induction/IMC still use the generic SEC
    // profile where in-solver strengthening can help monolithic UNSAT proofs.
    setKissatOptionOrThrow(solver, "otfs", 0);

    (void)coneSymbols;
    setKissatOptionOrThrow(solver, "preprocess", 0);
    setKissatOptionOrThrow(solver, "simplify", 0);
    setKissatOptionOrThrow(solver, "preprocesscongruence", 0);
    setKissatOptionOrThrow(solver, "preprocessprobe", 0);
    setKissatOptionOrThrow(solver, "congruence", 0);
    setKissatOptionOrThrow(solver, "probe", 0);
    setKissatOptionOrThrow(solver, "probeinit", 0);
    setKissatOptionOrThrow(solver, "eliminate", 0);
    setKissatOptionOrThrow(solver, "eliminateinit", 0);
  }

  void configureForSecPdrPersistentQuery(size_t coneSymbols = 0) {
    configureForSecPdrQuery(coneSymbols);
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::CADICAL) {
      return;  // LCOV_EXCL_LINE
    }
    // A shared PDR solver receives permanent units that retire old property
    // selectors. Re-enable CaDiCaL's exact clause/variable compaction so those
    // satisfied contexts do not accumulate across output batches. All other
    // PDR query settings and every logical clause remain unchanged.
    auto* solver = cadicalSolver_.get();
    setCadicalOptionIfSupported(solver, "compact", 1);
    setCadicalOptionIfSupported(solver, "arenacompact", 1);
  }

  void configureForSecLocalBooleanCheck(size_t coneSymbols = 0) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::CADICAL) {
      auto* solver = cadicalSolver_.get();
      // Structural SEC validators rebuild many small one-shot Boolean queries.
      // Use CaDiCaL as the assumption-capable default, but keep it on a direct
      // CDCL path so startup inference does not spend minutes in gate
      // extraction and inprocessing before the real proof engine starts.
      setCadicalOptionIfSupported(solver, "inprocessing", 0);
      setCadicalOptionIfSupported(solver, "compact", 0);
      setCadicalOptionIfSupported(solver, "arenacompact", 0);
      setCadicalOptionIfSupported(solver, "elim", 0);
      setCadicalOptionIfSupported(solver, "fastelim", 0);
      setCadicalOptionIfSupported(solver, "preprocesslight", 0);
      setCadicalOptionIfSupported(solver, "probe", 0);
      setCadicalOptionIfSupported(solver, "inprobing", 0);
      setCadicalOptionIfSupported(solver, "congruence", 0);
      setCadicalOptionIfSupported(solver, "decompose", 0);
      setCadicalOptionIfSupported(solver, "deduplicate", 0);
      setCadicalOptionIfSupported(solver, "factor", 0);
      setCadicalOptionIfSupported(solver, "subsume", 0);
      setCadicalOptionIfSupported(solver, "sweep", 0);
      setCadicalOptionIfSupported(solver, "lucky", 0);
      setCadicalOptionIfSupported(solver, "luckyearly", 0);
      setCadicalOptionIfSupported(solver, "luckylate", 0);
      setCadicalOptionIfSupported(solver, "minimize", 0);
      setCadicalOptionIfSupported(solver, "shrink", 0);
      setCadicalOptionIfSupported(solver, "rephase", 0);
      setCadicalOptionIfSupported(solver, "walk", 0);
      (void)coneSymbols;
      return;
    }

    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // LCOV_EXCL_START
      return;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    auto* solver = static_cast<kissat*>(kissatSolver_);
    // Local Boolean checks are short-lived validators used while building SEC
    // invariants and output slices. They are not the main proof engine, and
    // sampled sky130hs_ibex runs spent minutes in Kissat's speculative
    // preprocessing before PDR even started. Keep these checks on a direct CDCL
    // path; if the validator cannot decide quickly, the caller falls back to a
    // weaker but still sound SEC obligation.
    setKissatOptionOrThrow(solver, "stable", 0);
    setKissatOptionOrThrow(solver, "target", 2);
    setKissatOptionOrThrow(solver, "restartint", 10);
    setKissatOptionOrThrow(solver, "restartreusetrail", 0);
    setKissatOptionOrThrow(solver, "rephase", 0);
    setKissatOptionOrThrow(solver, "walkeffort", 0);
    setKissatOptionOrThrow(solver, "lucky", 0);
    setKissatOptionOrThrow(solver, "luckyearly", 0);
    setKissatOptionOrThrow(solver, "luckylate", 0);
    setKissatOptionOrThrow(solver, "minimize", 0);
    setKissatOptionOrThrow(solver, "shrink", 0);
    setKissatOptionOrThrow(solver, "otfs", 0);
    setKissatOptionOrThrow(solver, "preprocess", 0);
    setKissatOptionOrThrow(solver, "simplify", 0);
    setKissatOptionOrThrow(solver, "preprocesscongruence", 0);
    setKissatOptionOrThrow(solver, "preprocessprobe", 0);
    setKissatOptionOrThrow(solver, "congruence", 0);
    setKissatOptionOrThrow(solver, "probe", 0);
    setKissatOptionOrThrow(solver, "probeinit", 0);
    setKissatOptionOrThrow(solver, "eliminate", 0);
    setKissatOptionOrThrow(solver, "eliminateinit", 0);
    (void)coneSymbols;
  }

  // LCOV_EXCL_START
  void configureForSecResetExpressionProof(size_t coneSymbols = 0) {
    if (solverType_ != KEPLER_FORMAL::Config::SolverType::KISSAT) {
      return;  // LCOV_EXCL_LINE
      // LCOV_EXCL_STOP
    }

    // LCOV_EXCL_START
    auto* solver = static_cast<kissat*>(kissatSolver_);
    // LCOV_EXCL_STOP
    // Reset-expression checks are optional local UNSAT proofs over the
    // symbolic reset image.  They are reached only after PDR has produced a
    // candidate startup conflict, so use Kissat's proof-oriented stable mode
    // instead of the SAT-oriented predecessor-walk profile.  AES sampling
    // showed otherwise useful 500-symbol reset-image proofs spending their
    // runtime in focused CDCL propagation.
    // LCOV_EXCL_START
    setKissatOptionOrThrow(solver, "stable", 2);
    setKissatOptionOrThrow(solver, "target", 1);
    setKissatOptionOrThrow(solver, "rephase", 0);
    setKissatOptionOrThrow(solver, "walkeffort", 0);
    setKissatOptionOrThrow(solver, "lucky", 0);
    setKissatOptionOrThrow(solver, "luckyearly", 0);
    setKissatOptionOrThrow(solver, "luckylate", 0);
    // LCOV_EXCL_STOP
    // These local proofs are already bounded by an engine-level conflict
    // limit. Sampling the AES reset-image proof with a larger cap showed the
    // extra time going into recursive learned-clause shrinking/minimization
    // rather than useful propagation, so keep this shortcut on a cheap CDCL
    // path and let UNKNOWN fall through to exact validation.
    // LCOV_EXCL_START
    setKissatOptionOrThrow(solver, "minimize", 0);
    setKissatOptionOrThrow(solver, "shrink", 0);
    // LCOV_EXCL_STOP

    (void)coneSymbols;
    // Reset-expression solvers are short-lived and rebuilt for each PDR
    // candidate.  Sampling showed even moderate AES reset-image proofs spending
    // most wall time in speculative probe/sweep/kitten preprocessing, so keep
    // these local proof solvers on a direct CDCL path for every cone size.
    // LCOV_EXCL_START
    setKissatOptionOrThrow(solver, "preprocess", 0);
    setKissatOptionOrThrow(solver, "simplify", 0);
    setKissatOptionOrThrow(solver, "preprocesscongruence", 0);
    setKissatOptionOrThrow(solver, "preprocessprobe", 0);
    setKissatOptionOrThrow(solver, "congruence", 0);
    setKissatOptionOrThrow(solver, "probe", 0);
    setKissatOptionOrThrow(solver, "probeinit", 0);
    setKissatOptionOrThrow(solver, "eliminate", 0);
    setKissatOptionOrThrow(solver, "eliminateinit", 0);
  }
  // LCOV_EXCL_STOP

private:
  static CaDiCraig::CraigVarType cadicalCraigVariableType( // LCOV_EXCL_LINE
      CraigVariablePartition partition) {
    switch (partition) { // LCOV_EXCL_LINE
      case CraigVariablePartition::ALocal:
        return CaDiCraig::CraigVarType::A_LOCAL; // LCOV_EXCL_LINE
      case CraigVariablePartition::BLocal:
        return CaDiCraig::CraigVarType::B_LOCAL; // LCOV_EXCL_LINE
      case CraigVariablePartition::Global:
        return CaDiCraig::CraigVarType::GLOBAL; // LCOV_EXCL_LINE
    }
    return CaDiCraig::CraigVarType::A_LOCAL; // LCOV_EXCL_LINE
  } // LCOV_EXCL_LINE

  static bool setKissatOptionIfSupported(kissat* solver, const char* name, int value) {
    kissat_set_option(solver, name, value);
    return kissat_get_option(solver, name) == value;
  }

  static void setKissatOptionOrThrow(kissat* solver, const char* name, int value) {
    // Kissat option availability varies across the embedded builds used by
    // CMake, Bazel, and CI. These profiles are performance hints only; falling
    // back to the solver default is better than aborting an otherwise valid SEC
    // proof when an option such as `quiet` is not present in that build.
    (void)setKissatOptionIfSupported(solver, name, value);
  }

  static bool setCadicalOptionIfSupported( // LCOV_EXCL_LINE
      CaDiCaL::Solver* solver, const char* name, int value) {
    return CaDiCaL::Solver::is_valid_option(name) && solver->set(name, value); // LCOV_EXCL_LINE
  }

  static int64_t effectiveCadicalLimit(int64_t queryLimit,
                                       uint64_t remainingBudget) {
    const int64_t boundedRemaining = static_cast<int64_t>(
        std::min<uint64_t>(
            remainingBudget,
            static_cast<uint64_t>(std::numeric_limits<int>::max())));
    return queryLimit < 0
               ? boundedRemaining
               : std::min(queryLimit, boundedRemaining);
  }

  static int cadicalApiLimit(int64_t limit) {
    return limit < 0
               ? -1
               : static_cast<int>(
                     std::min<int64_t>(
                         limit, std::numeric_limits<int>::max()));
  }

  SolveStatus solveCadicalWithResourceLimits(
      int64_t conflictLimit,
      int64_t decisionLimit,
      int64_t tickLimit,
      bool useCumulativeBudget) {
    CadicalWorkBudget* budget =
        useCumulativeBudget ? activeCadicalWorkBudget_ : nullptr;
    if (budget != nullptr) {
      if (budget->exhausted()) {
        return SolveStatus::Unknown;
      }
      conflictLimit =
          effectiveCadicalLimit(conflictLimit, budget->remainingConflicts());
      decisionLimit =
          effectiveCadicalLimit(decisionLimit, budget->remainingDecisions());
      tickLimit =
          effectiveCadicalLimit(tickLimit, budget->remainingTicks());
    }

    // CaDiCaL retains incremental limits. Set all three on every solve so an
    // exact query cannot inherit a previous resource-aware query's limits.
    cadicalSolver_->limit("conflicts", cadicalApiLimit(conflictLimit));
    cadicalSolver_->limit("decisions", cadicalApiLimit(decisionLimit));
    cadicalSolver_->limit("ticks", cadicalApiLimit(tickLimit));

    int64_t conflictsBefore = 0;
    int64_t decisionsBefore = 0;
    int64_t ticksBefore = 0;
    if (budget != nullptr) {
      conflictsBefore = cadicalSolver_->get_statistic_value("conflicts");
      decisionsBefore = cadicalSolver_->get_statistic_value("decisions");
      ticksBefore = cadicalSolver_->get_statistic_value("ticks");
    }
    const int result = cadicalSolver_->solve();
    if (budget != nullptr) {
      const int64_t conflictsAfter =
          cadicalSolver_->get_statistic_value("conflicts");
      const int64_t decisionsAfter =
          cadicalSolver_->get_statistic_value("decisions");
      const int64_t ticksAfter =
          cadicalSolver_->get_statistic_value("ticks");
      budget->consume(
          conflictsAfter > conflictsBefore
              ? static_cast<uint64_t>(conflictsAfter - conflictsBefore)
              : 0,
          decisionsAfter > decisionsBefore
              ? static_cast<uint64_t>(decisionsAfter - decisionsBefore)
              : 0,
          ticksAfter > ticksBefore
              ? static_cast<uint64_t>(ticksAfter - ticksBefore)
              : 0);
    }

    if (result == 10) {
      return SolveStatus::Sat;
    }
    if (result == 20) {
      return SolveStatus::Unsat;
    }
    return SolveStatus::Unknown;
  }

  inline static thread_local CadicalWorkBudget* activeCadicalWorkBudget_ =
      nullptr;
  KEPLER_FORMAL::Config::SolverType solverType_;
  std::unique_ptr<Glucose::SimpSolver> glucoseSolver_;
  std::unique_ptr<CaDiCaL::Solver> cadicalSolver_;
  std::unique_ptr<CaDiCraig::CraigTracer> cadicalCraigTracer_;
  CraigVariablePartition craigVariablePartition_ =
      CraigVariablePartition::ALocal;
  CraigClausePartition craigClausePartition_ = CraigClausePartition::A;
  int cadicalCraigClauseId_ = 0;
  void* kissatSolver_ = nullptr;
  int kissatNumVars_ = 0;
  int kissatReservedVars_ = 0;
  int cadicalNumVars_ = 0;
  int cadicalReservedVars_ = 0;
  std::vector<int> lastAssumptions_;
  SolveStatus lastAssumptionSolveStatus_ = SolveStatus::Unknown;
};
