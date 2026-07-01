# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set(METRON_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/metron)
set(METROLIB_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/metrolib)
set(MATCHERONI_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/matcheroni)
set(XLS_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/xls)

option(KEPLER_ENABLE_XLS_FRONTEND "Enable embedded XLS/xlscc C-to-RTL frontend" OFF)
option(KEPLER_BUILD_XLS_FRONTEND
    "Build the embedded XLS/xlscc wrapper from thirdparty/xls with Bazel" OFF)
set(KEPLER_XLS_C2RTL_LIBRARY "" CACHE FILEPATH
    "Path to the prebuilt libkepler_xls_c2rtl wrapper library")

if(APPLE)
    set(KEPLER_XLS_C2RTL_LIBRARY_NAME libkepler_xls_c2rtl.dylib)
else()
    set(KEPLER_XLS_C2RTL_LIBRARY_NAME libkepler_xls_c2rtl.so)
endif()
set(KEPLER_XLS_C2RTL_DEFAULT_LIBRARY
    ${XLS_ROOT}/bazel-bin/xls/contrib/kepler/${KEPLER_XLS_C2RTL_LIBRARY_NAME})

if(KEPLER_BUILD_XLS_FRONTEND)
    set(KEPLER_ENABLE_XLS_FRONTEND ON CACHE BOOL
        "Enable embedded XLS/xlscc C-to-RTL frontend" FORCE)
    find_program(KEPLER_BAZEL_EXECUTABLE bazelisk)
    if(NOT KEPLER_BAZEL_EXECUTABLE)
        find_program(KEPLER_BAZEL_EXECUTABLE bazel)
    endif()
    if(NOT KEPLER_BAZEL_EXECUTABLE)
        message(FATAL_ERROR
            "KEPLER_BUILD_XLS_FRONTEND=ON requires bazelisk or bazel in PATH")
    endif()
    add_custom_command(
        OUTPUT ${KEPLER_XLS_C2RTL_DEFAULT_LIBRARY}
        COMMAND ${KEPLER_BAZEL_EXECUTABLE}
                --output_user_root=${CMAKE_BINARY_DIR}/xls-bazel
                build
                -c
                opt
                //xls/contrib/kepler:${KEPLER_XLS_C2RTL_LIBRARY_NAME}
        WORKING_DIRECTORY ${XLS_ROOT}
        COMMENT "Building embedded XLS/xlscc C-to-RTL wrapper"
        VERBATIM
    )
    add_custom_target(kepler_xls_c2rtl_bazel
        DEPENDS ${KEPLER_XLS_C2RTL_DEFAULT_LIBRARY})
    set(KEPLER_XLS_C2RTL_LIBRARY ${KEPLER_XLS_C2RTL_DEFAULT_LIBRARY}
        CACHE FILEPATH "Path to the libkepler_xls_c2rtl wrapper library" FORCE)
elseif(NOT KEPLER_XLS_C2RTL_LIBRARY AND EXISTS ${KEPLER_XLS_C2RTL_DEFAULT_LIBRARY})
    set(KEPLER_XLS_C2RTL_LIBRARY ${KEPLER_XLS_C2RTL_DEFAULT_LIBRARY}
        CACHE FILEPATH "Path to the libkepler_xls_c2rtl wrapper library" FORCE)
endif()

add_library(kepler_c2rtl
    MetronTranslator.cpp
    XlsTranslator.cpp
    ${METRON_ROOT}/metron/CContext.cpp
    ${METRON_ROOT}/metron/CInstance.cpp
    ${METRON_ROOT}/metron/CLexer.cpp
    ${METRON_ROOT}/metron/CNode.cpp
    ${METRON_ROOT}/metron/CParser.cpp
    ${METRON_ROOT}/metron/CScope.cpp
    ${METRON_ROOT}/metron/CSourceFile.cpp
    ${METRON_ROOT}/metron/CSourceRepo.cpp
    ${METRON_ROOT}/metron/CToken.cpp
    ${METRON_ROOT}/metron/Checker.cpp
    ${METRON_ROOT}/metron/Cursor.cpp
    ${METRON_ROOT}/metron/Dumper.cpp
    ${METRON_ROOT}/metron/Emitter.cpp
    ${METRON_ROOT}/metron/MtUtils.cpp
    ${METRON_ROOT}/metron/Tracer.cpp
    ${METRON_ROOT}/metron/main/main_new.cpp
    ${METRON_ROOT}/metron/metron_tools.cpp
    ${METRON_ROOT}/metron/nodes/CNodeClass.cpp
    ${METRON_ROOT}/metron/nodes/CNodeDeclaration.cpp
    ${METRON_ROOT}/metron/nodes/CNodeField.cpp
    ${METRON_ROOT}/metron/nodes/CNodeFunction.cpp
    ${METRON_ROOT}/metron/nodes/CNodeStruct.cpp
    ${METRON_ROOT}/metron/nodes/CNodeTranslationUnit.cpp
    ${METRON_ROOT}/metron/nodes/CNodeUnion.cpp
    ${METROLIB_ROOT}/metrolib/core/Blobs.cpp
    ${METROLIB_ROOT}/metrolib/core/Check.cpp
    ${METROLIB_ROOT}/metrolib/core/Constants.cpp
    ${METROLIB_ROOT}/metrolib/core/Dumper.cpp
    ${METROLIB_ROOT}/metrolib/core/Err.cpp
    ${METROLIB_ROOT}/metrolib/core/FieldInfo.cpp
    ${METROLIB_ROOT}/metrolib/core/File.cpp
    ${METROLIB_ROOT}/metrolib/core/Log.cpp
    ${METROLIB_ROOT}/metrolib/core/Platform.cpp
    ${METROLIB_ROOT}/metrolib/core/Result.cpp
    ${METROLIB_ROOT}/metrolib/core/StateStack.cpp
    ${METROLIB_ROOT}/metrolib/core/Tests.cpp
    ${METROLIB_ROOT}/metrolib/core/Types.cpp
    ${METROLIB_ROOT}/metrolib/core/Utils.cpp
)

target_include_directories(kepler_c2rtl
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
    SYSTEM PRIVATE
        ${METRON_ROOT}
        ${METRON_ROOT}/symlinks
        ${METROLIB_ROOT}
        ${MATCHERONI_ROOT}
)

target_compile_options(kepler_c2rtl PRIVATE
    $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-everything>
    $<$<CXX_COMPILER_ID:GNU>:-w>
)

target_compile_definitions(kepler_c2rtl
    PRIVATE
        KEPLER_METRON_ROOT="${METRON_ROOT}"
        KEPLER_XLS_ROOT="${XLS_ROOT}"
)

if(KEPLER_ENABLE_XLS_FRONTEND)
    if(NOT KEPLER_XLS_C2RTL_LIBRARY)
        message(FATAL_ERROR
            "KEPLER_ENABLE_XLS_FRONTEND=ON requires KEPLER_XLS_C2RTL_LIBRARY "
            "to point at the XLS wrapper library, or KEPLER_BUILD_XLS_FRONTEND=ON")
    endif()
    add_library(kepler_xls_c2rtl_runtime SHARED IMPORTED GLOBAL)
    set_target_properties(kepler_xls_c2rtl_runtime PROPERTIES
        IMPORTED_LOCATION "${KEPLER_XLS_C2RTL_LIBRARY}")
    if(TARGET kepler_xls_c2rtl_bazel)
        add_dependencies(kepler_xls_c2rtl_runtime kepler_xls_c2rtl_bazel)
        add_dependencies(kepler_c2rtl kepler_xls_c2rtl_bazel)
    endif()
    target_include_directories(kepler_c2rtl
        SYSTEM PRIVATE
            ${XLS_ROOT}
    )
    target_compile_definitions(kepler_c2rtl PRIVATE KEPLER_ENABLE_XLS_FRONTEND=1)
    target_link_libraries(kepler_c2rtl PRIVATE kepler_xls_c2rtl_runtime)
endif()
