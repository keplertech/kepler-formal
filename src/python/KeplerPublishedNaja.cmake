# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# Opt-in compatibility for the unmodified, published NajaEDA 0.7.24 wheels.
# Only headers are obtained from its verified source archive; all runtime
# targets below point to the libraries owned by the installed Python package.
set(KEPLER_NAJAEDA_SOURCE_ARCHIVE "" CACHE FILEPATH "Optional cached NajaEDA 0.7.24 source archive")
set(published_naja_args --cmake --output-dir "${CMAKE_BINARY_DIR}/published-najaeda")
if(KEPLER_NAJAEDA_SOURCE_ARCHIVE)
  list(APPEND published_naja_args --source-archive "${KEPLER_NAJAEDA_SOURCE_ARCHIVE}")
endif()
execute_process(
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/published_najaeda.py" ${published_naja_args}
  RESULT_VARIABLE provider_status
  OUTPUT_VARIABLE provider_config
  ERROR_VARIABLE provider_error
)
if(NOT provider_status EQUAL 0)
  message(FATAL_ERROR "Cannot use published NajaEDA: ${provider_error}")
endif()
set(provider_config_file "${CMAKE_BINARY_DIR}/published-najaeda/provider.cmake")
file(WRITE "${provider_config_file}" "${provider_config}")
include("${provider_config_file}")
find_package(TBB REQUIRED)
find_path(NajaEDA_BOOST_INCLUDE_DIR boost/intrusive/set.hpp REQUIRED)

foreach(dependency IN ITEMS tbb tbbmalloc)
  if(DEFINED NajaEDA_${dependency}_LIBRARY)
    get_filename_component(provider_name "${NajaEDA_${dependency}_LIBRARY}" NAME)
    foreach(suffix IN ITEMS "" _RELEASE _DEBUG _RELWITHDEBINFO _MINSIZEREL)
      set_property(TARGET TBB::${dependency} PROPERTY IMPORTED_LOCATION${suffix} "${NajaEDA_${dependency}_LIBRARY}")
      set_property(TARGET TBB::${dependency} PROPERTY IMPORTED_SONAME${suffix} "${provider_name}")
      if(WIN32)
        set_property(TARGET TBB::${dependency} PROPERTY IMPORTED_IMPLIB${suffix} "${NajaEDA_${dependency}_IMPLIB}")
      endif()
    endforeach()
  endif()
endforeach()

foreach(runtime IN ITEMS naja_nl naja_dnl naja_bne naja_opt naja_metrics naja_python)
  if(TARGET ${runtime})
    message(FATAL_ERROR "Published NajaEDA cannot share an existing ${runtime} target")
  endif()
  add_library(${runtime} SHARED IMPORTED GLOBAL)
  set_target_properties(${runtime} PROPERTIES
    IMPORTED_LOCATION "${NajaEDA_${runtime}_LIBRARY}"
    SYSTEM FALSE
    IMPORTED_NO_SYSTEM TRUE
    INTERFACE_INCLUDE_DIRECTORIES "${NajaEDA_INCLUDE_DIRS};${NajaEDA_BOOST_INCLUDE_DIR}"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${NajaEDA_BOOST_INCLUDE_DIR}"
    INTERFACE_COMPILE_FEATURES cxx_std_20)
  if(WIN32)
    set_property(TARGET ${runtime} PROPERTY IMPORTED_IMPLIB "${NajaEDA_${runtime}_IMPLIB}")
  endif()
endforeach()
target_link_libraries(naja_dnl INTERFACE naja_nl TBB::tbb TBB::tbbmalloc)
target_link_libraries(naja_bne INTERFACE naja_nl naja_dnl)
target_link_libraries(naja_opt INTERFACE naja_nl naja_dnl naja_bne)
target_link_libraries(naja_metrics INTERFACE naja_nl naja_dnl)
target_link_libraries(naja_python INTERFACE naja_nl naja_dnl naja_bne naja_opt naja_metrics Python3::Module)
foreach(component IN ITEMS naja_core naja_nl_dump naja_snl_liberty naja_snl_systemverilog naja_snl_verilog naja_snl_visual)
  add_library(${component} INTERFACE IMPORTED GLOBAL)
  if(component STREQUAL "naja_core")
    target_link_libraries(${component} INTERFACE naja_nl)
  else()
    target_link_libraries(${component} INTERFACE naja_python)
  endif()
endforeach()

function(najaeda_fixup_consumer target)
  if(APPLE)
    set(helper "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/published_najaeda.py")
    set_property(TARGET "${target}" APPEND PROPERTY LINK_DEPENDS "${helper}")
    add_custom_command(TARGET "${target}" POST_BUILD
      COMMAND "${Python3_EXECUTABLE}" "${helper}" --fixup-consumer "$<TARGET_FILE:${target}>"
      VERBATIM)
  endif()
endfunction()
