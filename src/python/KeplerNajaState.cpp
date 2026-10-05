// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "KeplerNajaState.h"

namespace naja::DNL {
// NajaEDA exports this storage but has no public exchange function. The
// provider identity checks guarantee the linked build owns this symbol.
#ifdef _WIN32
__declspec(dllimport) extern DNLFull* dnlFull_;
#else
extern DNLFull* dnlFull_;
#endif
}  // namespace naja::DNL

namespace KEPLER_FORMAL {

naja::DNL::DNLFull* exchangeNajaDNL(naja::DNL::DNLFull* replacement) noexcept {
  auto* previous = naja::DNL::dnlFull_;
  naja::DNL::dnlFull_ = replacement;
  return previous;
}

}  // namespace KEPLER_FORMAL
