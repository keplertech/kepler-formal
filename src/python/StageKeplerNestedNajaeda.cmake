# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

foreach(required_variable IN ITEMS
    SOURCE_PACKAGE_DIR
    DESTINATION_PACKAGE_DIR
    SHARED_LIBRARY_PREFIX
    SHARED_LIBRARY_SUFFIX
    MODULE_SUFFIX
    STAMP_FILE)
  if(NOT DEFINED ${required_variable} OR "${${required_variable}}" STREQUAL "")
    message(FATAL_ERROR "StageKeplerNestedNajaeda.cmake requires ${required_variable}")
  endif()
endforeach()

if(NOT IS_DIRECTORY "${SOURCE_PACKAGE_DIR}")
  message(FATAL_ERROR
    "The isolated NajaEDA install did not produce ${SOURCE_PACKAGE_DIR}")
endif()

get_filename_component(destination_parent "${DESTINATION_PACKAGE_DIR}" DIRECTORY)
get_filename_component(stamp_parent "${STAMP_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${destination_parent}" "${stamp_parent}")
file(REMOVE_RECURSE "${DESTINATION_PACKAGE_DIR}")
file(MAKE_DIRECTORY "${DESTINATION_PACKAGE_DIR}")
file(COPY "${SOURCE_PACKAGE_DIR}/"
  DESTINATION "${DESTINATION_PACKAGE_DIR}"
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
  PATTERN "*~" EXCLUDE
)

# Upstream NajaEDA is normally a top-level package and a few modules therefore
# use absolute imports.  Rewrite only the staged copy so it cannot resolve or
# create a process-global top-level `najaeda` alias.
file(GLOB_RECURSE python_sources
  LIST_DIRECTORIES FALSE
  "${DESTINATION_PACKAGE_DIR}/*.py"
)
foreach(python_source IN LISTS python_sources)
  file(READ "${python_source}" contents)
  set(rewritten "${contents}")
  string(REPLACE
    "from najaeda import"
    "from kepler_formal.najaeda import"
    rewritten
    "${rewritten}"
  )
  string(REPLACE
    "from najaeda."
    "from kepler_formal.najaeda."
    rewritten
    "${rewritten}"
  )
  if(NOT rewritten STREQUAL contents)
    file(WRITE "${python_source}" "${rewritten}")
  endif()
endforeach()

if(NOT EXISTS "${DESTINATION_PACKAGE_DIR}/__init__.py")
  message(FATAL_ERROR "The staged NajaEDA package has no __init__.py")
endif()

file(GLOB naja_modules
  "${DESTINATION_PACKAGE_DIR}/naja*${MODULE_SUFFIX}"
)
if(NOT naja_modules)
  message(FATAL_ERROR
    "The staged NajaEDA package has no native naja module ending in ${MODULE_SUFFIX}")
endif()

file(GLOB isolated_libraries
  "${DESTINATION_PACKAGE_DIR}/${SHARED_LIBRARY_PREFIX}*${SHARED_LIBRARY_SUFFIX}"
)
if(NOT isolated_libraries)
  message(FATAL_ERROR
    "The staged NajaEDA package has no ${SHARED_LIBRARY_PREFIX} runtime libraries")
endif()

file(WRITE "${STAMP_FILE}"
  "isolated NajaEDA package: ${DESTINATION_PACKAGE_DIR}\n"
)
