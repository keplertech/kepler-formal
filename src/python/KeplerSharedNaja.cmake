# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

# A Python build links the NajaEDA runtime installed in the build interpreter
# instead of compiling the vendored Naja. The default provider is a wheel built
# from thirdparty/naja, which ships its headers. KEPLER_USE_PUBLISHED_NAJAEDA
# selects the unmodified 0.7.24 release; its headers come from the verified
# source archive. Every runtime target below points to the installed package's
# own libraries; this build never modifies the provider.
set(KEPLER_NAJAEDA_SOURCE_ARCHIVE "" CACHE FILEPATH "Optional cached NajaEDA 0.7.24 source archive")
set(naja_provider_args --cmake --output-dir "${CMAKE_BINARY_DIR}/naja-provider")
if(KEPLER_USE_PUBLISHED_NAJAEDA)
  list(APPEND naja_provider_args --published)
  if(KEPLER_NAJAEDA_SOURCE_ARCHIVE)
    list(APPEND naja_provider_args --source-archive "${KEPLER_NAJAEDA_SOURCE_ARCHIVE}")
  endif()
endif()
# Ask the interpreter itself: a CMake override alone must not select a
# different runtime than Python will import.
execute_process(
  COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/naja_provider.py" ${naja_provider_args}
  RESULT_VARIABLE provider_status
  OUTPUT_VARIABLE provider_config
  ERROR_VARIABLE provider_error
)
if(NOT provider_status EQUAL 0)
  message(FATAL_ERROR
    "Python builds require an installed NajaEDA provider in ${Python3_EXECUTABLE} "
    "(see docs/python-api.md).\n${provider_error}")
endif()
set(provider_config_file "${CMAKE_BINARY_DIR}/naja-provider/provider.cmake")
file(WRITE "${provider_config_file}" "${provider_config}")
include("${provider_config_file}")

# These build-policy targets normally come from the vendored Naja project.
foreach(policy_target IN ITEMS coverage_config sanitizers_config)
  if(NOT TARGET "${policy_target}")
    add_library("${policy_target}" INTERFACE)
  endif()
endforeach()
find_path(Boost_INCLUDE_DIR NAMES boost/version.hpp REQUIRED)
set(Boost_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")
find_package(TBB REQUIRED)
if(WIN32)
  # The provider's repaired TBB import libraries are linked explicitly. Do not
  # let TBB headers also request the original tbb12.lib (or its debug variant).
  # Usage requirements must reach static dependencies as well as the module.
  foreach(dependency IN ITEMS tbb tbbmalloc)
    target_compile_definitions(TBB::${dependency} INTERFACE
      __TBB_NO_IMPLICIT_LINKAGE=1)
  endforeach()
endif()

# Keep the development package's TBB headers, but link the allocator and
# runtime bundled with the provider rather than introducing a second copy.
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
    message(FATAL_ERROR "The NajaEDA provider cannot share an existing ${runtime} target")
  endif()
  add_library(${runtime} SHARED IMPORTED GLOBAL)
  set_target_properties(${runtime} PROPERTIES
    IMPORTED_LOCATION "${NajaEDA_${runtime}_LIBRARY}"
    # The provider's headers must precede any older Naja headers on CPATH, so
    # they keep ordinary -I precedence instead of -isystem.
    SYSTEM FALSE
    IMPORTED_NO_SYSTEM TRUE
    INTERFACE_INCLUDE_DIRECTORIES "${NajaEDA_INCLUDE_DIRS};${Boost_INCLUDE_DIR}"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}"
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
# The vendored build's static components live inside the provider's shared
# libraries: naja_core in naja_nl, the format and dump components in naja_python.
foreach(component IN ITEMS naja_core naja_nl_dump naja_snl_liberty naja_snl_systemverilog naja_snl_verilog naja_snl_visual)
  add_library(${component} INTERFACE IMPORTED GLOBAL)
  if(component STREQUAL "naja_core")
    target_link_libraries(${component} INTERFACE naja_nl)
  else()
    target_link_libraries(${component} INTERFACE naja_python)
  endif()
endforeach()

# Repair-renamed Mach-O install IDs need an explicit post-link relocation in
# build trees. Call this from the directory that owns the consumer target.
function(najaeda_fixup_consumer target)
  if(APPLE)
    set(helper "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/naja_provider.py")
    set_property(TARGET "${target}" APPEND PROPERTY LINK_DEPENDS "${helper}")
    add_custom_command(TARGET "${target}" POST_BUILD
      COMMAND "${Python3_EXECUTABLE}" "${helper}" --fixup-consumer "$<TARGET_FILE:${target}>"
      VERBATIM)
  endif()
endfunction()
