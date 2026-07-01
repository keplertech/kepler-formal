# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

set(METRON_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/metron)
set(METROLIB_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/metrolib)
set(MATCHERONI_ROOT ${CMAKE_SOURCE_DIR}/thirdparty/matcheroni)

add_library(kepler_c2rtl
    MetronTranslator.cpp
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
)
