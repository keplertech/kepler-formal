// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerXlsC2Rtl.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "xls/codegen/codegen_options.h"
#include "xls/codegen/combinational_generator.h"
#include "xls/codegen/pipeline_generator.h"
#include "xls/codegen/vast/vast.h"
#include "xls/common/file/filesystem.h"
#include "xls/common/status/status_macros.h"
#include "xls/contrib/xlscc/hls_block.pb.h"
#include "xls/contrib/xlscc/translator.h"
#include "xls/contrib/xlscc/translator_types.h"
#include "xls/ir/function.h"
#include "xls/ir/function_base.h"
#include "xls/ir/nodes.h"
#include "xls/ir/op.h"
#include "xls/ir/package.h"
#include "xls/ir/proc.h"
#include "xls/ir/type.h"
#include "xls/ir/value_flattening.h"
#include "xls/passes/basic_simplification_pass.h"
#include "xls/passes/constant_folding_pass.h"
#include "xls/passes/cse_pass.h"
#include "xls/passes/dce_pass.h"
#include "xls/passes/inlining_pass.h"
#include "xls/passes/optimization_pass.h"
#include "xls/passes/pass_base.h"
#include "xls/scheduling/pipeline_schedule.h"

namespace {

struct ReferenceOutputMetadata {
  std::string source_name;
  std::string parameter_name;
  int64_t parameter_index = -1;
};

struct FunctionInterfaceMetadata {
  xls::Function* function = nullptr;
  std::vector<ReferenceOutputMetadata> reference_outputs;
};

struct PackageIrResult {
  std::string top_name;
  std::optional<FunctionInterfaceMetadata> function_interface;
};

struct AdapterOutput {
  std::string name;
  int64_t width = 0;
  int64_t lsb = 0;
};

struct ReferenceOutputAdapter {
  std::string public_module_name;
  std::string implementation_module_name;
  std::vector<int64_t> omitted_parameter_indices;
  std::vector<AdapterOutput> outputs;
  int64_t packed_output_width = 0;
};

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

std::optional<FunctionInterfaceMetadata> CaptureFunctionInterface(
    xlscc::GeneratedFunction* top_function) {
  if (top_function == nullptr || top_function->xls_func == nullptr ||
      top_function->clang_decl == nullptr ||
      !top_function->clang_decl->getReturnType()->isVoidType() ||
      !top_function->io_ops.empty() || top_function->this_lvalue != nullptr) {
    return std::nullopt;
  }

  FunctionInterfaceMetadata interface;
  interface.function = top_function->xls_func;
  for (const xlscc::GeneratedParamInfo& param_info :
       top_function->param_infos) {
    if (!param_info.return_val) {
      continue;
    }
    if (!param_info.has_func_param || param_info.param == nullptr) {
      return std::nullopt;
    }

    const std::string source_name = param_info.param->getNameAsString();
    if (source_name.empty()) {
      return std::nullopt;
    }
    absl::StatusOr<xls::Param*> parameter =
        interface.function->GetParamByName(source_name);
    if (!parameter.ok()) {
      return std::nullopt;
    }
    absl::StatusOr<int64_t> parameter_index =
        interface.function->GetParamIndex(*parameter);
    if (!parameter_index.ok()) {
      return std::nullopt;
    }
    interface.reference_outputs.push_back(ReferenceOutputMetadata{
        source_name,
        std::string((*parameter)->name()),
        *parameter_index,
    });
  }

  if (interface.reference_outputs.empty() ||
      top_function->return_value_count != interface.reference_outputs.size()) {
    return std::nullopt;
  }
  return interface;
}

absl::StatusOr<PackageIrResult> GeneratePackageIr(
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
    return PackageIrResult{proc->name(), std::nullopt};
  }

  absl::flat_hash_map<const clang::NamedDecl*, xlscc::ChannelBundle>
      top_channel_injections;
  XLS_ASSIGN_OR_RETURN(
      xlscc::GeneratedFunction * top_function,
      translator.GenerateIR_Top_Function(&package, top_channel_injections));
  XLS_RETURN_IF_ERROR(package.SetTopByName(entry_name));
  translator.AddSourceInfoToPackage(package);
  return PackageIrResult{entry_name, CaptureFunctionInterface(top_function)};
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

  XLS_ASSIGN_OR_RETURN(bool pass_changed,
                       RunOptimizationPass<xls::InliningPass>(
                           package, pass_options, &pass_results,
                           &pass_context));
  (void)pass_changed;

  bool changed = true;
  for (int64_t iteration = 0; changed && iteration < 4; ++iteration) {
    changed = false;

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

  return absl::OkStatus();
}

absl::StatusOr<std::optional<ReferenceOutputAdapter>>
BuildReferenceOutputAdapter(
    const FunctionInterfaceMetadata& interface,
    const KeplerXlsC2RtlOptions& options,
    std::string_view top_name) {
  if (interface.function == nullptr) {
    return std::nullopt;
  }

  ReferenceOutputAdapter adapter;
  const std::string requested_module_name =
      options.module_name != nullptr &&
              !std::string_view(options.module_name).empty()
          ? std::string(options.module_name)
          : std::string(top_name);
  adapter.public_module_name = xls::verilog::SanitizeVerilogIdentifier(
      requested_module_name, options.use_system_verilog != 0);
  adapter.implementation_module_name =
      xls::verilog::SanitizeVerilogIdentifier(
          absl::StrCat(adapter.public_module_name, "__xls_impl"),
          options.use_system_verilog != 0);

  for (const ReferenceOutputMetadata& output :
       interface.reference_outputs) {
    if (output.parameter_index < 0 ||
        output.parameter_index >= interface.function->params().size()) {
      return absl::InternalError(
          "XLSCC reference output parameter index is out of range");
    }
    xls::Param* parameter = interface.function->param(output.parameter_index);
    if (parameter->name() != output.parameter_name) {
      return absl::InternalError(
          "XLSCC reference output parameter order changed unexpectedly");
    }
    // A non-const C++ reference is an in/out parameter. It is safe to expose it
    // as an output-only RTL port only when the lowered function never reads the
    // incoming value.
    if (!parameter->users().empty() ||
        interface.function->HasImplicitUse(parameter)) {
      return std::nullopt;
    }
    adapter.omitted_parameter_indices.push_back(output.parameter_index);
  }

  xls::Type* return_type = interface.function->return_type();
  adapter.packed_output_width = return_type->GetFlatBitCount();
  if (adapter.packed_output_width <= 0) {
    return absl::InternalError("XLSCC reference outputs have zero width");
  }

  absl::flat_hash_set<std::string> output_names;
  if (interface.reference_outputs.size() == 1) {
    const std::string output_name = xls::verilog::SanitizeVerilogIdentifier(
        interface.reference_outputs.front().source_name,
        options.use_system_verilog != 0);
    adapter.outputs.push_back(
        AdapterOutput{output_name, adapter.packed_output_width, 0});
    output_names.insert(output_name);
  } else {
    if (!return_type->IsTuple()) {
      return absl::InternalError(
          "XLSCC reference outputs did not lower to a tuple");
    }
    const xls::TupleType* tuple_type = return_type->AsTupleOrDie();
    if (tuple_type->size() != interface.reference_outputs.size()) {
      return absl::InternalError(
          "XLSCC reference output tuple has an unexpected arity");
    }
    for (int64_t i = 0; i < tuple_type->size(); ++i) {
      const int64_t width = tuple_type->element_type(i)->GetFlatBitCount();
      if (width <= 0) {
        return absl::InternalError(
            "XLSCC reference output tuple contains a zero-width value");
      }
      const std::string output_name =
          xls::verilog::SanitizeVerilogIdentifier(
              interface.reference_outputs[i].source_name,
              options.use_system_verilog != 0);
      if (!output_names.insert(output_name).second) {
        return absl::InvalidArgumentError(
            "C/C++ reference output names collide after Verilog sanitization");
      }
      adapter.outputs.push_back(AdapterOutput{
          output_name,
          width,
          xls::GetFlatBitIndexOfElement(tuple_type, i),
      });
    }
  }

  return adapter;
}

std::string FormatPortDeclaration(std::string_view direction,
                                  int64_t width,
                                  std::string_view name) {
  if (width == 1) {
    return absl::StrCat("  ", direction, " wire ", name);
  }
  return absl::StrCat(
      "  ", direction, " wire [", width - 1, ":0] ", name);
}

std::string FormatZero(int64_t width) {
  return absl::StrCat(width, "'b0");
}

absl::StatusOr<std::string> EmitReferenceOutputAdapter(
    const ReferenceOutputAdapter& adapter,
    const xls::verilog::CodegenResult& codegen_result) {
  const auto data_outputs = codegen_result.signature.data_outputs();
  if (data_outputs.size() != 1 ||
      data_outputs.front().width() != adapter.packed_output_width) {
    return absl::InternalError(
        "XLS codegen produced an unexpected packed output interface");
  }

  const auto data_inputs = codegen_result.signature.data_inputs();
  std::vector<bool> omit_input(data_inputs.size(), false);
  for (int64_t parameter_index : adapter.omitted_parameter_indices) {
    if (parameter_index < 0 || parameter_index >= data_inputs.size()) {
      return absl::InternalError(
          "XLS codegen input interface no longer matches the C parameters");
    }
    omit_input[parameter_index] = true;
  }

  std::vector<std::string> declarations;
  absl::flat_hash_set<std::string> public_names;
  for (int64_t i = 0; i < data_inputs.size(); ++i) {
    const auto& input = data_inputs[i];
    if (input.width() == 0 || omit_input[i]) {
      continue;
    }
    if (!public_names.insert(input.name()).second) {
      return absl::InvalidArgumentError(
          "C/C++ interface names collide after Verilog sanitization");
    }
    declarations.push_back(
        FormatPortDeclaration("input", input.width(), input.name()));
  }
  for (const AdapterOutput& output : adapter.outputs) {
    if (!public_names.insert(output.name).second) {
      return absl::InvalidArgumentError(
          "C/C++ input and output names collide after Verilog sanitization");
    }
    declarations.push_back(
        FormatPortDeclaration("output", output.width, output.name));
  }

  const std::string packed_output_name = "__xls_packed_out";
  std::vector<std::string> connections;
  for (int64_t i = 0; i < data_inputs.size(); ++i) {
    const auto& input = data_inputs[i];
    if (input.width() == 0) {
      continue;
    }
    const std::string value =
        omit_input[i] ? FormatZero(input.width()) : input.name();
    connections.push_back(
        absl::StrCat("    .", input.name(), "(", value, ")"));
  }
  connections.push_back(absl::StrCat(
      "    .", data_outputs.front().name(), "(", packed_output_name, ")"));

  std::string wrapper = absl::StrCat(
      "\nmodule ", adapter.public_module_name, "(\n",
      absl::StrJoin(declarations, ",\n"), "\n);\n");
  if (adapter.packed_output_width == 1) {
    absl::StrAppend(&wrapper, "  wire ", packed_output_name, ";\n");
  } else {
    absl::StrAppend(&wrapper,
                    "  wire [", adapter.packed_output_width - 1, ":0] ",
                    packed_output_name, ";\n");
  }
  absl::StrAppend(
      &wrapper,
      "  ", codegen_result.signature.module_name(), " __xls_impl (\n",
      absl::StrJoin(connections, ",\n"), "\n  );\n");
  for (const AdapterOutput& output : adapter.outputs) {
    std::string value;
    if (adapter.packed_output_width == 1) {
      value = packed_output_name;
    } else if (output.width == 1) {
      value = absl::StrCat(packed_output_name, "[", output.lsb, "]");
    } else {
      value = absl::StrCat(
          packed_output_name, "[", output.lsb + output.width - 1, ":",
          output.lsb, "]");
    }
    absl::StrAppend(
        &wrapper, "  assign ", output.name, " = ", value, ";\n");
  }
  absl::StrAppend(&wrapper, "endmodule\n");
  return wrapper;
}

absl::Status TranslateToVerilog(const KeplerXlsC2RtlOptions& options) {
  XLS_RETURN_IF_ERROR(ValidateOptions(options));

  xls::Package package("kepler_xls_c2rtl");
  XLS_ASSIGN_OR_RETURN(PackageIrResult package_ir,
                       GeneratePackageIr(options, package));
  const std::string& top_name = package_ir.top_name;

  XLS_RETURN_IF_ERROR(RunPreCodegenPasses(&package));

  std::optional<ReferenceOutputAdapter> reference_output_adapter;
  if (package_ir.function_interface.has_value()) {
    XLS_ASSIGN_OR_RETURN(
        reference_output_adapter,
        BuildReferenceOutputAdapter(
            *package_ir.function_interface, options, top_name));
  }

  std::optional<xls::FunctionBase*> top = package.GetTop();
  if (!top.has_value()) {
    return absl::InternalError("XLS package has no top after C synthesis");
  }

  auto configure_common_codegen_options =
      [&](xls::verilog::CodegenOptions& codegen_options) {
        codegen_options.entry(top_name);
        codegen_options.use_system_verilog(options.use_system_verilog != 0);
        codegen_options.array_index_bounds_checking(false);
        codegen_options.SetOpOverride(
            xls::Op::kSMul,
            xls::verilog::OpOverrideAssignment(
                "assign {output} = $unsigned($signed({input0}) * "
                "$signed({input1}))"));
        codegen_options.SetOpOverride(
            xls::Op::kShra,
            xls::verilog::OpOverrideAssignment(
                "assign {output} = $unsigned($signed({input0}) >>> "
                "{input1})"));
        if (reference_output_adapter.has_value()) {
          codegen_options.module_name(
              reference_output_adapter->implementation_module_name);
        } else if (options.module_name != nullptr &&
                   !std::string_view(options.module_name).empty()) {
          codegen_options.module_name(options.module_name);
        }
      };

  xls::verilog::CodegenResult codegen_result;
  if (top.value()->IsProc()) {
    xls::verilog::CodegenOptions codegen_options =
        xls::verilog::BuildPipelineOptions();
    configure_common_codegen_options(codegen_options);
    codegen_options.use_system_verilog(false);
    codegen_options.clock_name("clk");
    codegen_options.reset("rst", /*asynchronous=*/false,
                          /*active_low=*/false, /*reset_data_path=*/false);
    codegen_options.emit_as_pipeline(false);
    XLS_ASSIGN_OR_RETURN(
        xls::PipelineSchedule schedule,
        xls::PipelineSchedule::SingleStage(top.value()));
    XLS_ASSIGN_OR_RETURN(
        codegen_result,
        xls::verilog::ToPipelineModuleText(schedule, top.value(),
                                           codegen_options));
  } else {
    xls::verilog::CodegenOptions codegen_options;
    configure_common_codegen_options(codegen_options);
    XLS_ASSIGN_OR_RETURN(
        codegen_result,
        xls::verilog::GenerateCombinationalModule(top.value(),
                                                  codegen_options));
  }

  std::string verilog_text = std::move(codegen_result.verilog_text);
  if (reference_output_adapter.has_value()) {
    XLS_ASSIGN_OR_RETURN(
        std::string wrapper,
        EmitReferenceOutputAdapter(*reference_output_adapter, codegen_result));
    absl::StrAppend(&verilog_text, wrapper);
  }
  return xls::SetFileContents(options.output_path, verilog_text);
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
