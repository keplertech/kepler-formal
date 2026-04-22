#pragma once

#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <memory>

#include "simp/SimpSolver.h"
extern "C" {
  #include "kissat.h"
}

#include "../config/Config.h"

//#define USE_KISSAT

class SATSolverWrapper {
public:

  explicit SATSolverWrapper(KEPLER_FORMAL::Config::SolverType type = KEPLER_FORMAL::Config::SolverType::GLUCOSE)
    : solverType_(type) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      glucoseSolver_ = std::make_unique<Glucose::SimpSolver>();
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      kissatSolver_ = kissat_init();
      kissatNumVars_ = 0;
    } else {
      // LCOV_EXCL_START
      throw std::invalid_argument("Unknown solver type");
      // LCOV_EXCL_STOP
    }
  }

  ~SATSolverWrapper() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT && kissatSolver_) {
      kissat_release(static_cast<kissat*>(kissatSolver_));
    }
  }

  // Create a new variable (returns 0-based index)
  int newVar() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      return glucoseSolver_->newVar();
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // Kissat does not require explicit variable creation, but we track max var.
      return kissatNumVars_++;
    }
    // LCOV_EXCL_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_EXCL_STOP
  }

  // Add a clause, literals are signed ints:
  // external convention: 0=false, 1=true, vars are ±(var_id + 2)
  void addClause(const std::vector<int>& lits) {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      Glucose::vec<Glucose::Lit> clause;
      for (int lit : lits) {
        if (lit == 0 || lit == 1) {
          // We should never see raw consts here: they are encoded via forced vars.
          // LCOV_EXCL_START
          throw std::runtime_error("Constant literal (0/1) passed to Glucose clause");
          // LCOV_EXCL_STOP
        }
        int v = std::abs(lit);
        int var = v - 2;  // external ±(var+2) -> internal var index
        if (var < 0) {
          // LCOV_EXCL_START
          throw std::runtime_error("Invalid literal (<2) passed to Glucose clause");
          // LCOV_EXCL_STOP
        }
        while (var >= glucoseSolver_->nVars())
          glucoseSolver_->newVar();
        clause.push((lit > 0) ? Glucose::mkLit(var) : ~Glucose::mkLit(var));
      }
      glucoseSolver_->addClause(clause);
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      // Kissat expects ±(var+1), 0 terminates a clause.
      for (int lit : lits) {
        if (lit == 0 || lit == 1) {
          // LCOV_EXCL_START
          throw std::runtime_error("Constant literal (0/1) passed to Kissat clause");
          // LCOV_EXCL_STOP
        }
        int v = std::abs(lit);
        int var = v - 2;  // external ±(var+2) -> internal var index
        if (var < 0) {
          // LCOV_EXCL_START
          throw std::runtime_error("Invalid literal (<2) passed to Kissat clause");
          // LCOV_EXCL_STOP
        }
        if (var >= kissatNumVars_) {
          kissatNumVars_ = var + 1;
        }
        int kissatLit = (lit > 0 ? var + 1 : -(var + 1)); // ±(var+1)
        kissat_add(static_cast<kissat*>(kissatSolver_), kissatLit);
      }
      kissat_add(static_cast<kissat*>(kissatSolver_), 0); // end of clause
    } else {
      // LCOV_EXCL_START
      throw std::runtime_error("Unknown solver type");
      // LCOV_EXCL_STOP
    }
  }

  bool solve() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      return glucoseSolver_->solve();
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      int res = kissat_solve(static_cast<kissat*>(kissatSolver_));
      return res == 10; // 10 = SAT, 20 = UNSAT
    }
    // LCOV_EXCL_START
    throw std::runtime_error("Unknown solver type");
    // LCOV_EXCL_STOP
  }

  bool getLiteralValue(int lit) const {
    if (lit == 0) {
      return false;
    }
    if (lit == 1) {
      return true;
    }

    const int external = std::abs(lit);
    const int var = external - 2;
    if (var < 0) {
      throw std::runtime_error("Invalid literal passed to getLiteralValue");
    }

    bool positiveValue = false;
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      const auto value = glucoseSolver_->modelValue(Glucose::mkLit(var));
      if (Glucose::toInt(value) == 2) {
        positiveValue = false;
      } else {
        positiveValue = Glucose::toInt(value) == 0;
      }
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      const int value = kissat_value(static_cast<kissat*>(kissatSolver_), var + 1);
      if (value == 0) {
        positiveValue = false;
      } else {
        positiveValue = value > 0;
      }
    } else {
      throw std::runtime_error("Unknown solver type");
    }

    return lit > 0 ? positiveValue : !positiveValue;
  }

  void* getSolver() {
    if (solverType_ == KEPLER_FORMAL::Config::SolverType::GLUCOSE) {
      return glucoseSolver_.get();
    } else if (solverType_ == KEPLER_FORMAL::Config::SolverType::KISSAT) {
      return kissatSolver_;
    }
    return nullptr;
  }

  KEPLER_FORMAL::Config::SolverType getSolverType() const { return solverType_; }

private:
  KEPLER_FORMAL::Config::SolverType solverType_;
  std::unique_ptr<Glucose::SimpSolver> glucoseSolver_;
  void* kissatSolver_ = nullptr;
  int kissatNumVars_ = 0;
};
