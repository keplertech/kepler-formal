

#pragma once

#include <mutex>

namespace KEPLER_FORMAL {

class Config {
public:
  enum SolverType {
    KISSAT,
    GLUCOSE,
    CADICAL
  };

  // Delete copy/move to enforce singleton semantics
  Config(const Config&) = delete;
  Config& operator=(const Config&) = delete;
  Config(Config&&) = delete;
  Config& operator=(Config&&) = delete;

  // Static configuration API
  // LCOV_EXCL_START
  static void setSolverType(SolverType type) {
    solverType_ = type;
  }
  // LCOV_EXCL_STOP

  static SolverType getSolverType() {
    return solverType_;
  }

  static void setReportSkippedPOs(bool enabled) {
    reportSkippedPOs_ = enabled;
  }

  static bool getReportSkippedPOs() {
    return reportSkippedPOs_;
  }

private:
  Config() = default;
  ~Config() = default;

  inline static SolverType solverType_ = KISSAT;
  inline static bool reportSkippedPOs_ = false;
};

} // namespace kepler
