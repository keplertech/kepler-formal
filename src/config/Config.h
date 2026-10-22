

#pragma once

#include <mutex>
#include <atomic>
#include <cstdint>
#include <stdexcept>

namespace KEPLER_FORMAL {

class Config {
public:
  enum SolverType {
    KISSAT,
    GLUCOSE,
    CADICAL
  };

  // Embedding calls need a fresh cache identity even when their live designs
  // and allocator-reused scratch objects have the same addresses as before.
  // Ordinary standalone execution remains in generation zero. Calls using
  // this context must be serialized, including all of their worker tasks.
  class ScopedVerificationContext {
   public:
    ScopedVerificationContext() {
      const auto generation = nextVerificationGeneration_.fetch_add(1);
      if (generation == 0) {
        throw std::overflow_error("Verification cache generation exhausted");
      }
      previous_ = verificationGeneration_.exchange(generation);
    }
    ~ScopedVerificationContext() {
      verificationGeneration_.store(previous_);
    }
    ScopedVerificationContext(const ScopedVerificationContext&) = delete;
    ScopedVerificationContext& operator=(const ScopedVerificationContext&) = delete;

   private:
    uint64_t previous_ = 0;
  };

  static uint64_t getVerificationGeneration() {
    return verificationGeneration_.load();
  }

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

  // SEC latch abstraction: clock-driven latches are reduced by phase
  // abstraction instead of staying opaque. On by default; off restores the
  // previous behavior where every latch output is opaque.
  static void setSecLatchAbstraction(bool enabled) {
    secLatchAbstraction_ = enabled;
  }

  static bool getSecLatchAbstraction() {
    return secLatchAbstraction_;
  }

  // SEC data-enabled latch model: a latch whose enable is data logic rather
  // than a clock phase is state with next value `enable ? data : q`, and its
  // readers see that value, so it is transparent within the step. On by
  // default; off leaves such latches opaque.
  static void setSecDataLatchModel(bool enabled) {
    secDataLatchModel_ = enabled;
  }

  static bool getSecDataLatchModel() {
    return secDataLatchModel_;
  }

private:
  Config() = default;
  ~Config() = default;

  inline static SolverType solverType_ = KISSAT;
  inline static bool reportSkippedPOs_ = false;
  inline static bool secLatchAbstraction_ = true;
  inline static bool secDataLatchModel_ = true;
  inline static std::atomic<uint64_t> nextVerificationGeneration_{1};
  inline static std::atomic<uint64_t> verificationGeneration_{0};
};

} // namespace kepler
