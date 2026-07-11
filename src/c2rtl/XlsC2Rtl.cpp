// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "xls/contrib/kepler/kepler_xls_c2rtl.h"

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "xls/codegen/codegen_options.h"
#include "xls/codegen/combinational_generator.h"
#include "xls/common/file/filesystem.h"
#include "xls/common/status/status_macros.h"
#include "xls/contrib/xlscc/hls_block.pb.h"
#include "xls/contrib/xlscc/translator.h"
#include "xls/contrib/xlscc/translator_types.h"
#include "xls/ir/function_base.h"
#include "xls/ir/op.h"
#include "xls/ir/package.h"
#include "xls/ir/proc.h"
#include "xls/passes/array_simplification_pass.h"
#include "xls/passes/basic_simplification_pass.h"
#include "xls/passes/constant_folding_pass.h"
#include "xls/passes/cse_pass.h"
#include "xls/passes/dce_pass.h"
#include "xls/passes/inlining_pass.h"
#include "xls/passes/optimization_pass.h"
#include "xls/passes/pass_base.h"

namespace {

absl::Status ValidateOptions(const KeplerXlsC2RtlOptions& options) {
  if (options.input_path == nullptr || std::string_view(options.input_path).empty()) {
    return absl::InvalidArgumentError("input_path is required");
  }
  if (options.output_path == nullptr || std::string_view(options.output_path).empty()) {
    return absl::InvalidArgumentError("output_path is required");
  }
  if (options.top == nullptr || std::string_view(options.top).empty()) {
    return absl::InvalidArgumentError("top is required");
  }
  if (options.include_path_count != 0 && options.include_paths == nullptr) {
    return absl::InvalidArgumentError(
        "include_paths is null but include_path_count is non-zero");
  }
  return absl::OkStatus();
}

std::vector<std::string> BuildClangArgs(const KeplerXlsC2RtlOptions& options) {
  std::vector<std::string> clang_args;
  clang_args.reserve(options.include_path_count);
  for (size_t i = 0; i < options.include_path_count; ++i) {
    if (options.include_paths[i] == nullptr ||
        std::string_view(options.include_paths[i]).empty()) {
      continue;
    }
    clang_args.push_back(absl::StrCat("-I", options.include_paths[i]));
  }
  return clang_args;
}

absl::StatusOr<std::string> GeneratePackageIr(
    const KeplerXlsC2RtlOptions& options, xls::Package& package) {
  xlscc::Translator translator(
      /*error_on_init_interval=*/false,
      /*error_on_uninitialized=*/false,
      /*generate_new_fsm=*/false,
      /*merge_states=*/true,
      /*split_states_on_channel_ops=*/true,
      xlscc::DebugIrTraceFlags_None,
      /*max_unroll_iters=*/1000,
      /*warn_unroll_iters=*/100,
      /*z3_rlimit=*/100000,
      xlscc::IOOpOrdering::kNone);

  XLS_RETURN_IF_ERROR(translator.SelectTop(options.top));

  std::vector<std::string> clang_arg_storage = BuildClangArgs(options);
  std::vector<std::string_view> clang_args;
  clang_args.reserve(clang_arg_storage.size());
  for (const auto& arg : clang_arg_storage) {
    clang_args.push_back(arg);
  }

  XLS_RETURN_IF_ERROR(translator.ScanFile(
      options.input_path, clang_args.empty()
                              ? absl::Span<std::string_view>()
                              : absl::MakeSpan(clang_args)));
  XLS_ASSIGN_OR_RETURN(std::string entry_name, translator.GetEntryFunctionName());

  if (options.block_proto_path != nullptr &&
      !std::string_view(options.block_proto_path).empty()) {
    xlscc::HLSBlock block;
    XLS_RETURN_IF_ERROR(xls::ParseTextProtoFile(options.block_proto_path, &block));
    xlscc::ChannelOptions channel_options;
    XLS_ASSIGN_OR_RETURN(
        xls::Proc * proc,
        translator.GenerateIR_Block(&package, block,
                                    /*top_level_init_interval=*/1,
                                    channel_options));
    XLS_RETURN_IF_ERROR(package.SetTop(proc));
    translator.AddSourceInfoToPackage(package);
    return proc->name();
  }

  absl::flat_hash_map<const clang::NamedDecl*, xlscc::ChannelBundle>
      top_channel_injections;
  XLS_ASSIGN_OR_RETURN(
      xlscc::GeneratedFunction * top_function,
      translator.GenerateIR_Top_Function(&package, top_channel_injections));
  (void)top_function;
  XLS_RETURN_IF_ERROR(package.SetTopByName(entry_name));
  translator.AddSourceInfoToPackage(package);
  return entry_name;
}

template <typename PassT>
absl::StatusOr<bool> RunOptimizationPass(
    xls::Package* package, const xls::OptimizationPassOptions& pass_options,
    xls::PassResults* pass_results, xls::OptimizationContext* pass_context) {
  PassT pass;
  return pass.Run(package, pass_options, pass_results, *pass_context);
}

absl::Status RunPreCodegenPasses(xls::Package* package) {
  xls::OptimizationPassOptions pass_options;
  xls::PassResults pass_results;
  xls::OptimizationContext pass_context;

  bool changed = true;
  int64_t iteration = 0;
  while (changed && iteration < 20) {
    changed = false;
    ++iteration;

    XLS_ASSIGN_OR_RETURN(bool pass_changed,
                         RunOptimizationPass<xls::InliningPass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
    XLS_ASSIGN_OR_RETURN(pass_changed,
                         RunOptimizationPass<xls::ArraySimplificationPass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
    XLS_ASSIGN_OR_RETURN(pass_changed,
                         RunOptimizationPass<xls::ConstantFoldingPass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
    XLS_ASSIGN_OR_RETURN(pass_changed,
                         RunOptimizationPass<xls::BasicSimplificationPass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
    XLS_ASSIGN_OR_RETURN(pass_changed,
                         RunOptimizationPass<xls::CsePass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
    XLS_ASSIGN_OR_RETURN(pass_changed,
                         RunOptimizationPass<xls::DeadCodeEliminationPass>(
                             package, pass_options, &pass_results,
                             &pass_context));
    changed |= pass_changed;
  }

  if (changed) {
    return absl::InternalError(
        "XLS C2RTL pre-codegen optimization did not converge");
  }
  return absl::OkStatus();
}

absl::Status TranslateToVerilog(const KeplerXlsC2RtlOptions& options) {
  XLS_RETURN_IF_ERROR(ValidateOptions(options));

  xls::Package package("kepler_xls_c2rtl");
  XLS_ASSIGN_OR_RETURN(std::string top_name, GeneratePackageIr(options, package));

  XLS_RETURN_IF_ERROR(RunPreCodegenPasses(&package));

  std::optional<xls::FunctionBase*> top = package.GetTop();
  if (!top.has_value()) {
    return absl::InternalError("XLS package has no top after C synthesis");
  }

  xls::verilog::CodegenOptions codegen_options;
  codegen_options.entry(top_name);
  codegen_options.use_system_verilog(options.use_system_verilog != 0);
  codegen_options.array_index_bounds_checking(false);
  codegen_options.SetOpOverride(
      xls::Op::kSMul,
      xls::verilog::OpOverrideAssignment(
          "assign {output} = $unsigned($signed({input0}) * $signed({input1}))"));
  codegen_options.SetOpOverride(
      xls::Op::kShra,
      xls::verilog::OpOverrideAssignment(
          "assign {output} = $unsigned($signed({input0}) >>> {input1})"));
  if (options.module_name != nullptr &&
      !std::string_view(options.module_name).empty()) {
    codegen_options.module_name(options.module_name);
  }

  XLS_ASSIGN_OR_RETURN(
      xls::verilog::CodegenResult codegen_result,
      xls::verilog::GenerateCombinationalModule(top.value(), codegen_options));

  return xls::SetFileContents(options.output_path, codegen_result.verilog_text);
}

char* CopyCString(std::string_view text) {
  char* result = static_cast<char*>(std::malloc(text.size() + 1));
  if (result == nullptr) {
    return nullptr;
  }
  std::memcpy(result, text.data(), text.size());
  result[text.size()] = '\0';
  return result;
}

}  // namespace

extern "C" int kepler_xls_c2rtl_translate(
    const KeplerXlsC2RtlOptions* options, char** error_message) {
  if (error_message != nullptr) {
    *error_message = nullptr;
  }
  if (options == nullptr) {
    if (error_message != nullptr) {
      *error_message = CopyCString("options is null");
    }
    return 1;
  }

  absl::Status status = TranslateToVerilog(*options);
  if (status.ok()) {
    return 0;
  }
  if (error_message != nullptr) {
    *error_message = CopyCString(status.ToString());
  }
  return 1;
}

extern "C" void kepler_xls_c2rtl_free(char* text) { std::free(text); }
