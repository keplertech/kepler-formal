// Copyright 2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0
#include "latch/LatchDependencyGraph.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <utility>

namespace KEPLER_FORMAL::SEC::LATCH {
namespace {

class DisjointSets {
 public:
  explicit DisjointSets(size_t size) : parent_(size), size_(size, 1) {
    std::iota(parent_.begin(), parent_.end(), 0);
  }
  size_t root(size_t member) {
    size_t result = member;
    while (result != parent_[result]) result = parent_[result];
    while (member != result) {
      const auto next = parent_[member];
      parent_[member] = result;
      member = next;
    }
    return result;
  }
  void join(size_t a, size_t b) {
    a = root(a);
    b = root(b);
    if (a == b) return;
    if (size_[a] < size_[b] || (size_[a] == size_[b] && a > b)) std::swap(a, b);
    parent_[b] = a;
    size_[a] += size_[b];
  }
 private:
  std::vector<size_t> parent_, size_;
};

void sortGroups(std::vector<std::vector<size_t>>& groups) {
  for (auto& group : groups) std::sort(group.begin(), group.end());
  std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
    return a.front() < b.front();
  });
}

}  // namespace

DependencyGraph analyzeDependencies(const Network& network) {
  const auto validNet = [&](size_t net) {
    if (net >= network.netCount) throw std::invalid_argument("Dependency graph net index out of range");
  };
  for (auto net : network.externalInputs) validNet(net);
  const size_t count = network.primitives.size();
  std::vector<std::vector<size_t>> writers(network.netCount), outgoing(count), incoming(count);
  for (size_t i = 0; i < count; ++i) {
    for (auto net : network.primitives[i].inputs) validNet(net);
    for (auto net : network.primitives[i].outputs) {
      validNet(net);
      writers[net].push_back(i);
    }
  }
  DisjointSets islands(count);
  for (auto& producers : writers) {
    // The input traversal visits cells in index order; repeated output pins can
    // repeat a writer, but must not create duplicate arcs or artificial loops.
    producers.erase(std::unique(producers.begin(), producers.end()), producers.end());
    for (size_t i = 1; i < producers.size(); ++i) islands.join(producers[0], producers[i]);
  }
  for (size_t consumer = 0; consumer < count; ++consumer) {
    for (auto net : network.primitives[consumer].inputs) {
      const auto& producers = writers[net];
      if (!producers.empty()) islands.join(producers[0], consumer);
      for (auto producer : producers) outgoing[producer].push_back(consumer);
    }
  }
  for (size_t source = 0; source < count; ++source) {
    auto& targets = outgoing[source];
    std::sort(targets.begin(), targets.end());
    targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
    for (auto target : targets) incoming[target].push_back(source);
  }

  DependencyGraph result;
  std::map<size_t, std::vector<size_t>> groups;
  for (size_t i = 0; i < count; ++i) groups[islands.root(i)].push_back(i);
  for (auto& [root, members] : groups) result.components.push_back(std::move(members));
  sortGroups(result.components);
  result.componentOf.resize(count);
  for (size_t i = 0; i < result.components.size(); ++i)
    for (auto member : result.components[i]) result.componentOf[member] = i;

  // Iterative Kosaraju: first collect DFS exit order, then traverse reverse arcs.
  // Explicit child cursors preserve exit order without recursive call frames.
  std::vector<uint8_t> visited(count, false);
  std::vector<size_t> finished;
  finished.reserve(count);
  std::vector<std::pair<size_t, size_t>> stack;
  for (size_t start = 0; start < count; ++start) {
    if (visited[start]) continue;
    visited[start] = true;
    stack.emplace_back(start, 0);
    while (!stack.empty()) {
      auto& [node, child] = stack.back();
      if (child == outgoing[node].size()) {
        finished.push_back(node);
        stack.pop_back();
      } else {
        const auto target = outgoing[node][child++];
        if (!visited[target]) {
          visited[target] = true;
          stack.emplace_back(target, 0);
        }
      }
    }
  }
  std::fill(visited.begin(), visited.end(), false);
  std::vector<size_t> pending;
  for (auto item = finished.rbegin(); item != finished.rend(); ++item) {
    if (visited[*item]) continue;
    std::vector<size_t> members;
    pending.push_back(*item);
    visited[*item] = true;
    while (!pending.empty()) {
      const auto node = pending.back();
      pending.pop_back();
      members.push_back(node);
      for (auto source : incoming[node]) {
        if (!visited[source]) {
          visited[source] = true;
          pending.push_back(source);
        }
      }
    }
    if (members.size() > 1 || std::binary_search(outgoing[members[0]].begin(),
                                                outgoing[members[0]].end(), members[0]))
      result.feedbackComponents.push_back(std::move(members));
  }
  sortGroups(result.feedbackComponents);
  return result;
}

}  // namespace KEPLER_FORMAL::SEC::LATCH
