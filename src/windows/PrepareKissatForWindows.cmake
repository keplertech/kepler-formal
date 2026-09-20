# Copyright 2026 keplertech.io
# SPDX-License-Identifier: Apache-2.0

function(kepler_prepare_kissat_for_windows source_dir destination_dir)
  # Quoted includes resolve beside each .c file, so copy sources and headers
  # together. Leave the pinned dependency untouched and track source updates.
  file(GLOB inputs CONFIGURE_DEPENDS "${source_dir}/*.c" "${source_dir}/*.h")
  # A dependency update may remove or rename files; do not keep compiling them.
  file(GLOB previous_copies "${destination_dir}/*.c" "${destination_dir}/*.h")
  foreach(previous IN LISTS previous_copies)
    get_filename_component(name "${previous}" NAME)
    if(NOT EXISTS "${source_dir}/${name}")
      file(REMOVE "${previous}")
    endif()
  endforeach()
  foreach(input IN LISTS inputs)
    get_filename_component(name "${input}" NAME)
    if(name STREQUAL "watch.h" OR name STREQUAL "clause.h")
      # Kissat stores these mixed bool/unsigned bitfields in unsigned words.
      # clang-cl's MSVC layout allocates separate units for different types;
      # -mno-ms-bitfields does not override that target's record layout.
      file(READ "${input}" contents)
      set(pattern "bool ([a-z_]+) : 1;")
      # Exclude semicolons from matches: CMake treats them as list separators.
      string(REGEX MATCHALL "bool [a-z_]+ : 1" fields "${contents}")
      list(LENGTH fields field_count)
      if(name STREQUAL "watch.h")
        set(expected_count 4) # binary tag in both structs and endian branches
      else()
        set(expected_count 8) # clause flags between glue and used
      endif()
      if(NOT field_count EQUAL expected_count)
        message(FATAL_ERROR "Kissat ${name} changed; review its Windows bitfield layout")
      endif()
      string(REGEX REPLACE "${pattern}" "unsigned \\1 : 1;" contents "${contents}")
      file(CONFIGURE OUTPUT "${destination_dir}/${name}" CONTENT "${contents}" @ONLY)
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${input}")
    else()
      configure_file("${input}" "${destination_dir}/${name}" COPYONLY)
    endif()
  endforeach()
endfunction()
