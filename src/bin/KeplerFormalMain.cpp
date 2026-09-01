// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerFormalDriver.h"

int main(int argc, char **argv) {
  const int rc = KeplerFormalMain(argc, argv);
#if defined(__SANITIZE_ADDRESS__)
  KEPLER_FORMAL::cleanupKeplerFormalState();
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
  KEPLER_FORMAL::cleanupKeplerFormalState();
#endif
#endif
  return rc;
}
