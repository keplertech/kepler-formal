// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "xls/passes/tools/passes_profile.h"

namespace xls {

void RecordPassEntry(std::string_view short_name) {}

void RecordPassAnnotation(std::string_view key,
                          std::variant<std::string_view, int64_t> contents) {}

void ExitPass(bool changed) {}

}  // namespace xls
