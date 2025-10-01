#include "BoolExprCache.h"
#include "BoolExpr.h"
#include <tbb/concurrent_unordered_map.h>
#include <memory>
#include <type_traits>

namespace KEPLER_FORMAL {

// ensure Op is hashable
struct OpHash {
  size_t operator()(Op op) const noexcept {
    using UT = std::underlying_type_t<Op>;
    return std::hash<UT>()(static_cast<UT>(op));
  }
};

// define NodeID to match BoolExpr::nodeID type (use size_t here; change if different)
using NodeID = size_t;

// clearer nested aliases: final value is std::shared_ptr<BoolExpr>
using Level4 = tbb::concurrent_unordered_map<NodeID, std::shared_ptr<BoolExpr>, std::hash<NodeID>, std::equal_to<NodeID>>;
using Level3 = tbb::concurrent_unordered_map<NodeID, Level4, std::hash<NodeID>, std::equal_to<NodeID>>;
using Level2 = tbb::concurrent_unordered_map<size_t, Level3, std::hash<size_t>, std::equal_to<size_t>>;
using NotTable = tbb::concurrent_unordered_map<Op, Level2, OpHash, std::equal_to<Op>>;

struct BoolExprCache::Impl {
  NotTable notTable;
};

BoolExprCache::Impl& BoolExprCache::impl() {
  static Impl instance;
  return instance;
}

std::shared_ptr<BoolExpr> BoolExprCache::getExpression(Key const& k) {
  NodeID lid = k.l ? k.l->getId() : NodeID{0};
  NodeID rid = k.r ? k.r->getId() : NodeID{0};

  auto &tbl = impl().notTable;

  // 1) find Op
  auto it1 = tbl.find(k.op);
  if (it1 != tbl.end()) {
    auto &lvl2 = it1->second; // Level2
    auto it2 = lvl2.find(k.varId);
    if (it2 != lvl2.end()) {
      auto &lvl3 = it2->second; // Level3
      auto it3 = lvl3.find(lid);
      if (it3 != lvl3.end()) {
        auto &lvl4 = it3->second; // Level4
        auto it4 = lvl4.find(rid);
        if (it4 != lvl4.end()) {
          // it4->second is std::shared_ptr<BoolExpr>
          return it4->second;
        }
      }
    }
  }
  std::shared_ptr<BoolExpr> L = k.l ? k.l->shared_from_this() : nullptr;
  std::shared_ptr<BoolExpr> R = k.r ? k.r->shared_from_this() : nullptr;
  // not found: create via the declared factory that takes Key
  // createNode is declared as: static std::shared_ptr<BoolExpr> createNode(BoolExprCache::Key const& k);
  // ensure BoolExpr grants friendship to BoolExprCache or createNode is public
  auto ptr = std::shared_ptr<BoolExpr>(
        new BoolExpr(k.op, k.varId, std::move(L), std::move(R))
    );

  // insert into nested maps (operator[] will create missing intermediate maps)
  impl().notTable[k.op][k.varId][lid][rid] = ptr;

  return ptr;
}

} // namespace KEPLER_FORMAL
