// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "KeplerNajaState.h"

#ifdef KEPLER_USE_PUBLISHED_NAJAEDA
namespace naja::DNL {
// NajaEDA 0.7.24 exports this storage but has no public exchange function.
// This version-specific bridge is enabled only with the published adapter.
#ifdef _WIN32
__declspec(dllimport) extern DNLFull* dnlFull_;
#else
extern DNLFull* dnlFull_;
#endif
}  // namespace naja::DNL
#endif

namespace KEPLER_FORMAL {

naja::DNL::DNLFull* exchangeNajaDNL(naja::DNL::DNLFull* replacement) noexcept {
#ifdef KEPLER_USE_PUBLISHED_NAJAEDA
  auto* previous = naja::DNL::dnlFull_;
  naja::DNL::dnlFull_ = replacement;
  return previous;
#else
  return naja::DNL::exchange(replacement);
#endif
}

}  // namespace KEPLER_FORMAL
