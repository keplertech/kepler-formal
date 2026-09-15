// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#ifndef KEPLER_C2RTL_KEPLER_XLS_C2RTL_H_
#define KEPLER_C2RTL_KEPLER_XLS_C2RTL_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct KeplerXlsC2RtlOptions {
  const char* input_path;
  const char* output_path;
  const char* top;
  const char* module_name;
  const char* block_proto_path;
  const char* const* include_paths;
  size_t include_path_count;
  int use_system_verilog;
} KeplerXlsC2RtlOptions;

int kepler_xls_c2rtl_translate(const KeplerXlsC2RtlOptions* options,
                               char** error_message);
void kepler_xls_c2rtl_free(char* text);

#ifdef __cplusplus
}
#endif

#endif  // KEPLER_C2RTL_KEPLER_XLS_C2RTL_H_
