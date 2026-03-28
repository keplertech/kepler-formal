

#pragma once

#include <mutex>

namespace KEPLER_FORMAL {

class Config {
public:
  enum SolverType {
    KISSAT,
    GLUCOSE
  };

  // Delete copy/move to enforce singleton semantics
  Config(const Config&) = delete;
  Config& operator=(const Config&) = delete;
  Config(Config&&) = delete;
  Config& operator=(Config&&) = delete;

  // Static configuration API
  static void setSolverType(SolverType type) {
    getInstance().solverType_ = type;
  }

  static SolverType getSolverType() {
    return getInstance().solverType_;
  }

  static void setReportSkippedPOs(bool enabled) {
    getInstance().reportSkippedPOs_ = enabled;
  }

  static bool getReportSkippedPOs() {
    return getInstance().reportSkippedPOs_;
  }

private:
  Config() = default;
  ~Config() = default;

  static Config& getInstance() {
    static Config instance;
    return instance;
  }

  SolverType solverType_ = KISSAT;
  bool reportSkippedPOs_ = false;
};

} // namespace kepler
