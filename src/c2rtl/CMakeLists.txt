# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_library(kepler_xls_c2rtl STATIC
  XlsC2Rtl.cpp
  XlsPipelineScheduleShim.cpp
  XlsPassProfileShim.cpp
  XlsProcElaborationStub.cpp
)

target_compile_definitions(kepler_xls_c2rtl PRIVATE
  KEPLER_C2RTL_LOCAL_PIPELINE_SIGNAL_NAME=1
)

# Compatibility sources are generated under the build tree before they are
# added to the C2RTL source list.
set(XLS_C2RTL_COMPAT_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/compat)

set(XLS_C2RTL_XLS_SOURCES
  ${KEPLER_XLS_ROOT}/xls/codegen/block_metrics.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/bdd_io_analysis.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/clone_nodes_into_block_handler.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/concurrent_stage_groups.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_checker.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_options.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_util.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_wrapper_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/expression_flattening.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/ffi_instantiation_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/finite_state_machine.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/lint_annotate.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/mark_channel_fifos_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/maybe_materialize_fifos_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/module_builder.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/module_signature.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/mulp_combining_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/name_legalization_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/node_expressions.cc
  ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/codegen/op_override_impls.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/port_legalization_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/priority_select_reduction_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/ram_configuration.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/register_chaining_analysis.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/register_combining_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/register_legalization_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/side_effect_condition_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/signature_generation_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/signature_generator.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/trace_verbosity_pass.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/verilog_conversion.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/vast/fold_vast_constants.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/vast/infer_vast_types.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/vast/vast.cc
  ${KEPLER_XLS_ROOT}/xls/codegen/vast/verilog_keywords.cc
  ${KEPLER_XLS_ROOT}/xls/common/case_converters.cc
  ${KEPLER_XLS_ROOT}/xls/common/attribute_data.cc
  ${KEPLER_XLS_ROOT}/xls/common/file/filesystem.cc
  ${KEPLER_XLS_ROOT}/xls/common/file/temp_file.cc
  ${KEPLER_XLS_ROOT}/xls/common/file/path.cc
  ${KEPLER_XLS_ROOT}/xls/common/indent.cc
  ${KEPLER_XLS_ROOT}/xls/common/logging/log_lines.cc
  ${KEPLER_XLS_ROOT}/xls/common/math_util.cc
  ${KEPLER_XLS_ROOT}/xls/common/status/error_code_to_status.cc
  ${KEPLER_XLS_ROOT}/xls/common/status/ret_check.cc
  ${KEPLER_XLS_ROOT}/xls/common/status/status_builder.cc
  ${KEPLER_XLS_ROOT}/xls/common/status/status_builder_oss.cc
  ${KEPLER_XLS_ROOT}/xls/common/stopwatch.cc
  ${KEPLER_XLS_ROOT}/xls/common/strerror.cc
  ${KEPLER_XLS_ROOT}/xls/common/string_to_int.cc
  ${KEPLER_XLS_ROOT}/xls/common/symbolized_stacktrace.cc
  ${KEPLER_XLS_ROOT}/xls/common/to_string_utils.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/cc_parser.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/continuations.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/expr_clone.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/generate_fsm.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/node_manipulation.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/tracked_bvalue.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translate_block.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translate_io.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translate_loops.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translate_metadata.cc
  ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/contrib/xlscc/translator.cc
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translator_types.cc
  ${KEPLER_XLS_ROOT}/xls/data_structures/binary_decision_diagram.cc
  ${KEPLER_XLS_ROOT}/xls/data_structures/binary_search.cc
  ${KEPLER_XLS_ROOT}/xls/data_structures/inline_bitmap.cc
  ${KEPLER_XLS_ROOT}/xls/data_structures/leaf_type_tree.cc
  ${KEPLER_XLS_ROOT}/xls/data_structures/transitive_closure.cc
  ${KEPLER_XLS_ROOT}/xls/dev_tools/link_to_source.cc
  ${KEPLER_XLS_ROOT}/xls/estimators/delay_model/delay_estimator.cc
  ${KEPLER_XLS_ROOT}/xls/estimators/delay_model/delay_estimators.cc
  ${KEPLER_XLS_ROOT}/xls/estimators/delay_model/ffi_delay_estimator.cc
  ${KEPLER_XLS_ROOT}/xls/interpreter/function_interpreter.cc
  ${KEPLER_XLS_ROOT}/xls/interpreter/ir_interpreter.cc
  ${KEPLER_XLS_ROOT}/xls/ir/big_int.cc
  ${KEPLER_XLS_ROOT}/xls/ir/bit_push_buffer.cc
  ${KEPLER_XLS_ROOT}/xls/ir/bits.cc
  ${KEPLER_XLS_ROOT}/xls/ir/bits_ops.cc
  ${KEPLER_XLS_ROOT}/xls/ir/block.cc
  ${KEPLER_XLS_ROOT}/xls/ir/block_elaboration.cc
  ${KEPLER_XLS_ROOT}/xls/ir/call_graph.cc
  ${KEPLER_XLS_ROOT}/xls/ir/caret.cc
  ${KEPLER_XLS_ROOT}/xls/ir/channel.cc
  ${KEPLER_XLS_ROOT}/xls/ir/channel_ops.cc
  ${KEPLER_XLS_ROOT}/xls/ir/clone_package.cc
  ${KEPLER_XLS_ROOT}/xls/ir/code_template.cc
  ${KEPLER_XLS_ROOT}/xls/ir/dfs_visitor.cc
  ${KEPLER_XLS_ROOT}/xls/ir/elaborated_block_dfs_visitor.cc
  ${KEPLER_XLS_ROOT}/xls/ir/elaboration.cc
  ${KEPLER_XLS_ROOT}/xls/ir/events.cc
  ${KEPLER_XLS_ROOT}/xls/ir/foreign_function.cc
  ${KEPLER_XLS_ROOT}/xls/ir/format_preference.cc
  ${KEPLER_XLS_ROOT}/xls/ir/format_strings.cc
  ${KEPLER_XLS_ROOT}/xls/ir/function.cc
  ${KEPLER_XLS_ROOT}/xls/ir/function_base.cc
  ${KEPLER_XLS_ROOT}/xls/ir/function_builder.cc
  ${KEPLER_XLS_ROOT}/xls/ir/instantiation.cc
  ${KEPLER_XLS_ROOT}/xls/ir/interval.cc
  ${KEPLER_XLS_ROOT}/xls/ir/interval_ops.cc
  ${KEPLER_XLS_ROOT}/xls/ir/interval_set.cc
  ${KEPLER_XLS_ROOT}/xls/ir/ir_annotator.cc
  ${KEPLER_XLS_ROOT}/xls/ir/ir_parser.cc
  ${KEPLER_XLS_ROOT}/xls/ir/ir_scanner.cc
  ${KEPLER_XLS_ROOT}/xls/ir/keyword_args.cc
  ${KEPLER_XLS_ROOT}/xls/ir/name_uniquer.cc
  ${KEPLER_XLS_ROOT}/xls/ir/node.cc
  ${KEPLER_XLS_ROOT}/xls/ir/node_util.cc
  ${KEPLER_XLS_ROOT}/xls/ir/nodes.cc
  ${KEPLER_XLS_ROOT}/xls/ir/number_parser.cc
  ${KEPLER_XLS_ROOT}/xls/ir/op.cc
  ${KEPLER_XLS_ROOT}/xls/ir/package.cc
  ${KEPLER_XLS_ROOT}/xls/ir/partial_information.cc
  ${KEPLER_XLS_ROOT}/xls/ir/partial_ops.cc
  ${KEPLER_XLS_ROOT}/xls/ir/proc.cc
  ${KEPLER_XLS_ROOT}/xls/ir/proc_conversion.cc
  ${KEPLER_XLS_ROOT}/xls/ir/proc_instantiation.cc
  ${KEPLER_XLS_ROOT}/xls/ir/register.cc
  ${KEPLER_XLS_ROOT}/xls/ir/scheduled_builder.cc
  ${KEPLER_XLS_ROOT}/xls/ir/state_element.cc
  ${KEPLER_XLS_ROOT}/xls/ir/ternary.cc
  ${KEPLER_XLS_ROOT}/xls/ir/topo_sort.cc
  ${KEPLER_XLS_ROOT}/xls/ir/type.cc
  ${KEPLER_XLS_ROOT}/xls/ir/type_manager.cc
  ${KEPLER_XLS_ROOT}/xls/ir/value.cc
  ${KEPLER_XLS_ROOT}/xls/ir/value_builder.cc
  ${KEPLER_XLS_ROOT}/xls/ir/value_conversion_utils.cc
  ${KEPLER_XLS_ROOT}/xls/ir/value_flattening.cc
  ${KEPLER_XLS_ROOT}/xls/ir/value_utils.cc
  ${KEPLER_XLS_ROOT}/xls/ir/verifier.cc
  ${KEPLER_XLS_ROOT}/xls/ir/verify_node.cc
  ${KEPLER_XLS_ROOT}/xls/passes/array_simplification_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/basic_simplification_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/bdd_query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/passes/constant_folding_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/cse_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/dce_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/dataflow_simplification_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/identity_removal_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/inlining_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/lazy_ternary_query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/passes/node_dependency_analysis.cc
  ${KEPLER_XLS_ROOT}/xls/passes/optimization_pass.cc
  ${KEPLER_XLS_ROOT}/xls/passes/optimization_pass_registry.cc
  ${KEPLER_XLS_ROOT}/xls/passes/partial_info_query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/passes/pass_base.cc
  ${KEPLER_XLS_ROOT}/xls/passes/query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/passes/stateless_query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/passes/union_query_engine.cc
  ${KEPLER_XLS_ROOT}/xls/solvers/solver.cc
  ${KEPLER_XLS_ROOT}/xls/solvers/z3_ir_translator.cc
  ${KEPLER_XLS_ROOT}/xls/solvers/z3_op_translator.cc
  ${KEPLER_XLS_ROOT}/xls/solvers/z3_utils.cc
)

set(XLS_C2RTL_PROTO_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${XLS_C2RTL_PROTO_OUT_DIR})

file(MAKE_DIRECTORY ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/passes)
file(MAKE_DIRECTORY ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/codegen)
file(MAKE_DIRECTORY ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/contrib/xlscc)
file(MAKE_DIRECTORY ${XLS_C2RTL_COMPAT_OUT_DIR}/ortools/graph)
file(WRITE ${XLS_C2RTL_COMPAT_OUT_DIR}/ortools/graph/graph.h
  "#ifndef ORTOOLS_GRAPH_GRAPH_H_\n"
  "#define ORTOOLS_GRAPH_GRAPH_H_\n"
  "namespace util {\n"
  "template <typename NodeIndexType, typename ArcIndexType>\n"
  "class ReverseArcStaticGraph {};\n"
  "}  // namespace util\n"
  "#endif  // ORTOOLS_GRAPH_GRAPH_H_\n"
)

# Keep thirdparty/xls pristine while smoothing over standalone CMake builds with
# dependency combinations that the XLS Bazel build does not expose directly.
file(READ ${KEPLER_XLS_ROOT}/xls/passes/pass_base.h XLS_C2RTL_PASS_BASE_H)
string(REPLACE
  "const PassBase<OptionsT, ContextT...>::InvariantChecker* checker)"
  "const typename PassBase<OptionsT, ContextT...>::InvariantChecker* checker)"
  XLS_C2RTL_PASS_BASE_H
  "${XLS_C2RTL_PASS_BASE_H}"
)
string(REPLACE
  "using InvariantChecker = PassBase<OptionsT, ContextT...>::InvariantChecker;"
  "using InvariantChecker = typename PassBase<OptionsT, ContextT...>::InvariantChecker;"
  XLS_C2RTL_PASS_BASE_H
  "${XLS_C2RTL_PASS_BASE_H}"
)
file(WRITE ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/passes/pass_base.h
  "${XLS_C2RTL_PASS_BASE_H}"
)

file(READ ${KEPLER_XLS_ROOT}/xls/passes/pipeline_generator.h XLS_C2RTL_PIPELINE_GENERATOR_H)
string(REPLACE
  "if (element.has_pass_name()) {"
  "if (element.type_case() == PassPipelineProto::Element::kPassName) {"
  XLS_C2RTL_PIPELINE_GENERATOR_H
  "${XLS_C2RTL_PIPELINE_GENERATOR_H}"
)
file(WRITE ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/passes/pipeline_generator.h
  "${XLS_C2RTL_PIPELINE_GENERATOR_H}"
)

file(READ ${KEPLER_XLS_ROOT}/xls/codegen/op_override_impls.cc
  XLS_C2RTL_OP_OVERRIDE_IMPLS_CC)
string(REPLACE
  "std::string_view piece(fmt_string);"
  "re2::StringPiece piece(fmt_string.data(), fmt_string.size());"
  XLS_C2RTL_OP_OVERRIDE_IMPLS_CC
  "${XLS_C2RTL_OP_OVERRIDE_IMPLS_CC}"
)
file(WRITE ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/codegen/op_override_impls.cc
  "${XLS_C2RTL_OP_OVERRIDE_IMPLS_CC}"
)

file(READ ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/translator.cc
  XLS_C2RTL_TRANSLATOR_CC)
string(REPLACE
  "std::string_view asm_view(sasm);"
  "re2::StringPiece asm_view(sasm.data(), sasm.size());"
  XLS_C2RTL_TRANSLATOR_CC
  "${XLS_C2RTL_TRANSLATOR_CC}"
)
file(WRITE ${XLS_C2RTL_COMPAT_OUT_DIR}/xls/contrib/xlscc/translator.cc
  "${XLS_C2RTL_TRANSLATOR_CC}"
)

set(XLS_C2RTL_PROTOS
  ${CMAKE_CURRENT_SOURCE_DIR}/proto/google/api/field_behavior.proto
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/hls_block.proto
  ${KEPLER_XLS_ROOT}/xls/contrib/xlscc/metadata_output.proto
  ${KEPLER_XLS_ROOT}/xls/codegen/codegen_residual_data.proto
  ${KEPLER_XLS_ROOT}/xls/codegen/module_signature.proto
  ${KEPLER_XLS_ROOT}/xls/codegen/verilog_line_map.proto
  ${KEPLER_XLS_ROOT}/xls/codegen/xls_metrics.proto
  ${KEPLER_XLS_ROOT}/xls/ir/channel.proto
  ${KEPLER_XLS_ROOT}/xls/ir/evaluator_result.proto
  ${KEPLER_XLS_ROOT}/xls/ir/foreign_function_data.proto
  ${KEPLER_XLS_ROOT}/xls/ir/op.proto
  ${KEPLER_XLS_ROOT}/xls/ir/ram_rewrite.proto
  ${KEPLER_XLS_ROOT}/xls/ir/transform_metrics.proto
  ${KEPLER_XLS_ROOT}/xls/ir/xls_ir_interface.proto
  ${KEPLER_XLS_ROOT}/xls/ir/xls_type.proto
  ${KEPLER_XLS_ROOT}/xls/ir/xls_value.proto
  ${KEPLER_XLS_ROOT}/xls/jit/type_layout.proto
  ${KEPLER_XLS_ROOT}/xls/netlist/netlist.proto
  ${KEPLER_XLS_ROOT}/xls/passes/pass_metrics.proto
  ${KEPLER_XLS_ROOT}/xls/passes/pass_pipeline.proto
  ${KEPLER_XLS_ROOT}/xls/passes/optimization_pass_pipeline.proto
  ${KEPLER_XLS_ROOT}/xls/scheduling/pipeline_schedule.proto
)

target_sources(kepler_xls_c2rtl PRIVATE
  ${XLS_C2RTL_XLS_SOURCES}
  ${XLS_C2RTL_PROTOS}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(kepler_xls_c2rtl PRIVATE
    -fno-lto
    -O0
    -Wno-deprecated-declarations
    -Wno-inconsistent-missing-override
    -Wno-nullability-completeness
    -Wno-unused-result
    -Wno-unused-value
  )
endif()

set(XLS_C2RTL_PROTO_IMPORT_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/proto
  ${KEPLER_XLS_ROOT}
  ${Protobuf_INCLUDE_DIRS}
)

set(XLS_C2RTL_PROTO_IMPORT_ARGS)
foreach(_import_dir IN LISTS XLS_C2RTL_PROTO_IMPORT_DIRS)
  list(APPEND XLS_C2RTL_PROTO_IMPORT_ARGS -I ${_import_dir})
endforeach()

if(TARGET protobuf::protoc)
  set(XLS_C2RTL_PROTOC $<TARGET_FILE:protobuf::protoc>)
else()
  set(XLS_C2RTL_PROTOC ${Protobuf_PROTOC_EXECUTABLE})
endif()
if(NOT XLS_C2RTL_PROTOC)
  message(FATAL_ERROR "protoc is required for XLS C2RTL")
endif()

set(XLS_C2RTL_PROTO_GENERATED_SRCS)
foreach(_proto IN LISTS XLS_C2RTL_PROTOS)
  set(_proto_rel "")
  foreach(_import_dir IN LISTS XLS_C2RTL_PROTO_IMPORT_DIRS)
    file(RELATIVE_PATH _candidate "${_import_dir}" "${_proto}")
    if(NOT _candidate MATCHES "^\\.\\." AND NOT IS_ABSOLUTE "${_candidate}")
      set(_proto_rel "${_candidate}")
      break()
    endif()
  endforeach()
  if(_proto_rel STREQUAL "")
    message(FATAL_ERROR "Could not map proto ${_proto} to an import root")
  endif()

  get_filename_component(_proto_rel_dir "${_proto_rel}" DIRECTORY)
  get_filename_component(_proto_name_we "${_proto_rel}" NAME_WE)
  if(_proto_rel_dir STREQUAL "")
    set(_proto_out_base "${XLS_C2RTL_PROTO_OUT_DIR}/${_proto_name_we}")
  else()
    set(_proto_out_base "${XLS_C2RTL_PROTO_OUT_DIR}/${_proto_rel_dir}/${_proto_name_we}")
  endif()

  add_custom_command(
    OUTPUT
      ${_proto_out_base}.pb.cc
      ${_proto_out_base}.pb.h
    COMMAND ${XLS_C2RTL_PROTOC}
    ARGS
      --experimental_allow_proto3_optional
      --cpp_out ${XLS_C2RTL_PROTO_OUT_DIR}
      ${XLS_C2RTL_PROTO_IMPORT_ARGS}
      ${_proto}
    DEPENDS ${_proto}
    COMMENT "Running cpp protocol buffer compiler on ${_proto}"
    VERBATIM
  )
  list(APPEND XLS_C2RTL_PROTO_GENERATED_SRCS
    ${_proto_out_base}.pb.cc
    ${_proto_out_base}.pb.h
  )
endforeach()

target_sources(kepler_xls_c2rtl PRIVATE ${XLS_C2RTL_PROTO_GENERATED_SRCS})

target_include_directories(kepler_xls_c2rtl
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${KEPLER_XLS_ROOT}
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/compat
    ${XLS_C2RTL_PROTO_OUT_DIR}
    ${LLVM_INCLUDE_DIRS}
    ${CLANG_INCLUDE_DIRS}
    ${CLANG_AST_INCLUDE_DIR}
    ${OPENSSL_INCLUDE_DIR}
)

target_include_directories(kepler_xls_c2rtl BEFORE PRIVATE
  ${XLS_C2RTL_COMPAT_OUT_DIR}
)

target_link_libraries(kepler_xls_c2rtl
  PRIVATE
    protobuf::libprotobuf
    OpenSSL::Crypto
    cppitertools
    clang-cpp
    LLVM
    re2::re2
    z3::z3
    absl::base
    absl::btree
    absl::cleanup
    absl::die_if_null
    absl::flat_hash_map
    absl::flat_hash_set
    absl::inlined_vector
    absl::log
    absl::check
    absl::status
    absl::statusor
    absl::strings
    absl::str_format
    absl::stacktrace
    absl::synchronization
    absl::span
    absl::symbolize
    absl::time
)
