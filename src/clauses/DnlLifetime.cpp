#include "DnlLifetime.h"

#include <atomic>

#include "DNL.h"

namespace KEPLER_FORMAL {

namespace {
std::atomic<uint64_t> dnlEpochCounter{1};
}  // namespace

uint64_t dnlEpoch() {
  return dnlEpochCounter.load(std::memory_order_acquire);
}

void destroyDnl() {
  naja::DNL::destroy();
  dnlEpochCounter.fetch_add(1, std::memory_order_acq_rel);
}

}  // namespace KEPLER_FORMAL
