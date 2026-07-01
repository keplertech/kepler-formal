// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "XlsTranslator.h"

#ifdef KEPLER_ENABLE_XLS_FRONTEND
#include "xls/contrib/kepler/kepler_xls_c2rtl.h"
#endif

#include <fstream>
#include <stdexcept>
#include <vector>

namespace KEPLER_FORMAL::C2RTL {
namespace {

void writeStartLogs(const XlsTranslationOptions& options) {
  std::ofstream stdoutLog(options.stdoutLog, std::ios::out | std::ios::trunc);
  if (stdoutLog) {
    stdoutLog << "kepler-formal embedded XLS C2RTL frontend\n";
    stdoutLog << "input: " << options.inputPath.string() << "\n";
    stdoutLog << "output: " << options.outputPath.string() << "\n";
    stdoutLog << "top: " << options.top << "\n";
    if (options.moduleName) {
      stdoutLog << "module_name: " << *options.moduleName << "\n";
    }
    if (options.blockProto) {
      stdoutLog << "block_proto: " << options.blockProto->string() << "\n";
    }
    stdoutLog << "include_paths:\n";
    for (const auto& includePath : options.includePaths) {
      stdoutLog << "  " << includePath.string() << "\n";
    }
  }
}

void writeDisabledLog(const XlsTranslationOptions& options) {
  std::ofstream stderrLog(options.stderrLog, std::ios::out | std::ios::trunc);
  if (stderrLog) {
    stderrLog
        << "The XLS frontend is selected and thirdparty/xls is vendored, but "
           "the XLS translator libraries are not linked into kepler-formal yet.\n"
        << "The intended implementation is in-process: xlscc translation API, "
           "XLS IR optimization API, and XLS codegen API, with no xlscc/opt_main/"
           "codegen_main subprocesses.\n";
    if (options.blockProto) {
      stderrLog << "XLS block textproto requested: "
                << options.blockProto->string() << "\n";
    }
  }
}

#ifdef KEPLER_ENABLE_XLS_FRONTEND
std::vector<std::string> buildIncludePathStrings(
    const XlsTranslationOptions& options) {
  std::vector<std::string> includeStrings;
  includeStrings.reserve(options.includePaths.size() + 1);
  for (const auto& includePath : options.includePaths) {
    includeStrings.push_back(includePath.string());
  }
#ifdef KEPLER_XLS_ROOT
  includeStrings.push_back(
      (std::filesystem::path(KEPLER_XLS_ROOT) /
       "xls/contrib/xlscc/synth_only")
          .string());
#endif
  return includeStrings;
}
#endif

}  // namespace

void translateWithXls(const XlsTranslationOptions& options) {
  writeStartLogs(options);

#ifdef KEPLER_ENABLE_XLS_FRONTEND
  const auto includeStrings = buildIncludePathStrings(options);
  std::vector<const char*> includePaths;
  includePaths.reserve(includeStrings.size());
  for (const auto& includeString : includeStrings) {
    includePaths.push_back(includeString.c_str());
  }

  const std::string inputPath = options.inputPath.string();
  const std::string outputPath = options.outputPath.string();
  const std::string moduleName = options.moduleName.value_or("");
  const std::string blockProto =
      options.blockProto ? options.blockProto->string() : std::string();

  KeplerXlsC2RtlOptions cOptions;
  cOptions.input_path = inputPath.c_str();
  cOptions.output_path = outputPath.c_str();
  cOptions.top = options.top.c_str();
  cOptions.module_name = moduleName.c_str();
  cOptions.block_proto_path = blockProto.c_str();
  cOptions.include_paths = includePaths.data();
  cOptions.include_path_count = includePaths.size();
  cOptions.use_system_verilog = 1;

  char* errorMessage = nullptr;
  const int rc = kepler_xls_c2rtl_translate(&cOptions, &errorMessage);
  if (rc != 0) {
    std::string message =
        errorMessage != nullptr ? std::string(errorMessage)
                                : std::string("unknown XLS frontend error");
    kepler_xls_c2rtl_free(errorMessage);
    std::ofstream stderrLog(options.stderrLog, std::ios::out | std::ios::trunc);
    if (stderrLog) {
      stderrLog << message << "\n";
    }
    throw std::runtime_error(message);
  }
  kepler_xls_c2rtl_free(errorMessage);
#else
  writeDisabledLog(options);

  throw std::runtime_error(
      "embedded XLS frontend is not linked yet; thirdparty/xls is present, but "
      "kepler_c2rtl still needs the XLS xlscc/optimizer/codegen libraries wired "
      "into the build");
#endif
}

}  // namespace KEPLER_FORMAL::C2RTL
