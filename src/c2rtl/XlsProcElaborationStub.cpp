// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "xls/ir/proc_elaboration.h"

#include "absl/status/status.h"
#include "absl/status/statusor.h"

namespace xls {

absl::StatusOr<ProcElaboration> ProcElaboration::Elaborate(Proc *) {
  return absl::UnimplementedError(
      "Proc elaboration is not linked into kepler-formal C2RTL CMake builds");
}

absl::Span<ProcInstance *const> ProcElaboration::GetInstances(Proc *) const {
  return {};
}

absl::Span<ChannelInstance *const>
ProcElaboration::GetInstances(Channel *) const {
  return {};
}

} // namespace xls
