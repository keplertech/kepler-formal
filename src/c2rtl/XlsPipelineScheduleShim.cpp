// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "xls/scheduling/pipeline_schedule.h"

#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_join.h"
#include "absl/types/span.h"
#include "xls/common/status/status_macros.h"
#include "xls/ir/function.h"
#include "xls/ir/function_base.h"
#include "xls/ir/node.h"
#include "xls/ir/nodes.h"
#include "xls/ir/proc.h"
#include "xls/ir/topo_sort.h"

namespace xls {
namespace {

int64_t MaximumCycle(const ScheduleCycleMap& cycle_map) {
  int64_t max_cycle = 0;
  for (const auto& [_, cycle] : cycle_map) {
    max_cycle = std::max(max_cycle, cycle);
  }
  return max_cycle;
}

}  // namespace

PipelineSchedule::PipelineSchedule(
    FunctionBase* function_base, ScheduleCycleMap cycle_map,
    std::optional<int64_t> min_clock_period_ps,
    std::optional<int64_t> target_clock_period_ps)
    : function_base_(function_base),
      cycle_map_(std::move(cycle_map)),
      min_clock_period_ps_(min_clock_period_ps),
      target_clock_period_ps_(target_clock_period_ps) {}

absl::StatusOr<PipelineSchedule> PipelineSchedule::Create(
    FunctionBase* function_base, ScheduleCycleMap cycle_map, Options options) {
  PipelineSchedule schedule(function_base, std::move(cycle_map),
                            options.min_clock_period_ps,
                            options.target_clock_period_ps);

  int64_t max_cycle = MaximumCycle(schedule.cycle_map_);
  if (options.length.has_value()) {
    max_cycle = *options.length - 1;
  }
  schedule.cycle_to_nodes_.resize(max_cycle + 1);
  for (const auto& [node, cycle] : schedule.cycle_map_) {
    schedule.cycle_to_nodes_[cycle].push_back(node);
  }

  absl::flat_hash_map<Node*, int64_t> topo_index;
  XLS_ASSIGN_OR_RETURN(std::vector<Node*> topo_sort_nodes,
                       TopoSort(function_base));
  for (int64_t i = 0; i < static_cast<int64_t>(topo_sort_nodes.size()); ++i) {
    topo_index[topo_sort_nodes[i]] = i;
  }
  for (std::vector<Node*>& nodes_in_cycle : schedule.cycle_to_nodes_) {
    std::sort(nodes_in_cycle.begin(), nodes_in_cycle.end(),
              [&](Node* a, Node* b) {
                return topo_index.at(a) < topo_index.at(b);
              });
  }
  return schedule;
}

absl::StatusOr<PipelineSchedule> PipelineSchedule::SingleStage(
    FunctionBase* function) {
  ScheduleCycleMap cycle_map;
  for (Node* node : function->nodes()) {
    cycle_map.emplace(node, 0);
  }
  return PipelineSchedule::Create(function, std::move(cycle_map));
}

absl::Span<Node* const> PipelineSchedule::nodes_in_cycle(int64_t cycle) const {
  if (cycle < cycle_to_nodes_.size()) {
    return cycle_to_nodes_[cycle];
  }
  return absl::Span<Node* const>();
}

bool PipelineSchedule::IsLiveOutOfCycle(Node* node, int64_t c) const {
  Function* as_func = dynamic_cast<Function*>(function_base_);
  if (cycle(node) > c || c >= length() - 1) {
    return false;
  }
  if (as_func != nullptr && node == as_func->return_value()) {
    return true;
  }
  for (Node* user : node->users()) {
    if (cycle(user) <= c) {
      continue;
    }
    if (user->Is<Next>()) {
      Next* next = user->As<Next>();
      if (next->predicate() != node && next->value() != node) {
        continue;
      }
    }
    return true;
  }
  return false;
}

std::vector<Node*> PipelineSchedule::GetLiveOutOfCycle(int64_t c) const {
  std::vector<Node*> live_out;
  for (int64_t i = 0; i <= c; ++i) {
    for (Node* node : nodes_in_cycle(i)) {
      if (IsLiveOutOfCycle(node, c)) {
        live_out.push_back(node);
      }
    }
  }
  return live_out;
}

std::string PipelineSchedule::ToString() const {
  std::vector<std::string> cycles;
  cycles.reserve(length());
  for (int64_t cycle = 0; cycle < length(); ++cycle) {
    cycles.push_back(absl::StrFormat("cycle %d: %d nodes", cycle,
                                     nodes_in_cycle(cycle).size()));
  }
  return absl::StrFormat("PipelineSchedule for `%s` (length=%d): %s",
                         function_base()->name(), length(),
                         absl::StrJoin(cycles, ", "));
}

std::string PackageSchedule::ToString() const {
  std::vector<std::string> schedules;
  for (FunctionBase* fb : GetScheduledFunctionBases()) {
    schedules.push_back(GetSchedule(fb).ToString());
  }
  return absl::StrFormat("PackageSchedule for `%s`: %s", package_->name(),
                         absl::StrJoin(schedules, "\n"));
}

}  // namespace xls
