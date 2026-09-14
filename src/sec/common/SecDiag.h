// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#pragma once

#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>
#ifdef _WIN32
#include <climits>
#include <cstdio>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace KEPLER_FORMAL::SEC {

inline bool isSecDiagEnabled() {
  return std::getenv("KEPLER_SEC_DIAG") != nullptr;
}

// LCOV_EXCL_START
inline bool isSecDiagOutputEnabled() {
// LCOV_EXCL_STOP
  // Keep normal SEC runs stderr-clean: wrappers such as Tcl exec can treat any
  // raw stderr as failure even when Kepler exits successfully.
  // LCOV_EXCL_START
  return isSecDiagEnabled() ||
         std::getenv("KEPLER_SEC_KI_DIAG") != nullptr ||
         std::getenv("KEPLER_SEC_KI_COI_DIAG") != nullptr ||
         std::getenv("KEPLER_SEC_PDR_RESET_SHORTCUT_DIAG") != nullptr ||
         std::getenv("KEPLER_SEC_PDR_STATS") != nullptr ||
         std::getenv("KEPLER_SEC_PDR_TRACE") != nullptr ||
         std::getenv("KEPLER_SEC_SUMMARY_STATS") != nullptr;
         // LCOV_EXCL_STOP
}

// LCOV_EXCL_START
inline void appendSecDiagPart(std::ostringstream& stream, const char* value) {
  stream << (value != nullptr ? value : "<null>");
}
// LCOV_EXCL_STOP

inline void appendSecDiagPart(std::ostringstream& stream, char* value) {
  stream << (value != nullptr ? value : "<null>");
}

// Diagnostic formatting templates instantiate many one-off call shapes from
// optional debug paths; line coverage is tracked at the call sites instead.
// LCOV_EXCL_START
template <typename T>
inline void appendSecDiagPart(std::ostringstream& stream, T&& value) {
  stream << std::forward<T>(value);
}

template <typename... Args>
inline void emitSecDiag(Args&&... args) {
  if (!isSecDiagOutputEnabled()) {
    return;
  }
  std::ostringstream stream;
  (appendSecDiagPart(stream, std::forward<Args>(args)), ...);
  stream << '\n';
  const std::string message = stream.str();
  const char* data = message.data();
  size_t remaining = message.size();
  while (remaining > 0) {
#ifdef _WIN32
    const auto count = static_cast<unsigned int>(
        remaining > INT_MAX ? INT_MAX : remaining);
    const int written = ::_write(::_fileno(stderr), data, count);
#else
    const ssize_t written = ::write(STDERR_FILENO, data, remaining);
#endif
    if (written <= 0) {
      // LCOV_DISABLED_START
      break;  // LCOV_EXCL_LINE
      // LCOV_DISABLED_STOP
    }
    data += written;
    remaining -= static_cast<size_t>(written);
  }
}
// LCOV_EXCL_STOP

}  // namespace KEPLER_FORMAL::SEC
