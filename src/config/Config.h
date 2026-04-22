

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
    solverType_ = type;
  }

  static SolverType getSolverType() {
    return solverType_;
  }

  static void setReportSkippedPOs(bool enabled) {
    reportSkippedPOs_ = enabled;
  }

  static bool getReportSkippedPOs() {
    return reportSkippedPOs_;
  }

  static void setSecTreatUncomputableSeqAsBoundary(bool enabled) {
    secTreatUncomputableSeqAsBoundary_ = enabled;
  }

  static bool getSecTreatUncomputableSeqAsBoundary() {
    return secTreatUncomputableSeqAsBoundary_;
  }

private:
  Config() = default;
  ~Config() = default;

  inline static SolverType solverType_ = KISSAT;
  inline static bool reportSkippedPOs_ = false;
  inline static bool secTreatUncomputableSeqAsBoundary_ = true;
};

} // namespace kepler
