# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# CMake's platform initialization resets CMAKE_SHARED_LIBRARY_PREFIX during
# project().  External Naja and its subprojects include this file immediately
# after each project() call, which keeps every Naja shared library in the
# nested runtime distinct from the libraries used by Kepler Formal itself.
if(KEPLER_NESTED_NAJAEDA_BUILD)
  set(CMAKE_SHARED_LIBRARY_PREFIX "libkepler_najaeda_")
endif()
