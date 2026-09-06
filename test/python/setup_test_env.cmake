# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

file(REMOVE_RECURSE "${DST}")
file(MAKE_DIRECTORY "${DST}")
if(DEFINED NAJAEDA_SRC AND
   NOT "${NAJAEDA_SRC}" STREQUAL "" AND
   EXISTS "${NAJAEDA_SRC}")
  file(MAKE_DIRECTORY "${DST}/najaeda")
  file(COPY "${NAJAEDA_SRC}/" DESTINATION "${DST}/najaeda")
endif()
# Add the ordinary Kepler package files after the independently staged nested
# NajaEDA tree. Their installed paths do not overlap.
file(COPY "${SRC}/" DESTINATION "${DST}")
file(COPY "${NATIVE}" DESTINATION "${DST}")
