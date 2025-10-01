#include "BoolExprCache.h"
#include <tbb/concurrent_vector.h>
#include <memory>
#include <type_traits>
#include "BoolExpr.h"

namespace KEPLER_FORMAL {

// ensure Op is hashable
struct OpHash {
  size_t operator()(Op op) const noexcept {
    using UT = std::underlying_type_t<Op>;
    return std::hash<UT>()(static_cast<UT>(op));
  }
};

// define NodeID to match BoolExpr::nodeID type (use size_t here; change if
// different)
using NodeID = size_t;

// clearer nested aliases: final value is std::shared_ptr<BoolExpr>
using Level4 = tbb::concurrent_vector<std::shared_ptr<BoolExpr>>;
using Level3 = tbb::concurrent_vector<Level4>;
using Level2 = tbb::concurrent_vector<Level3>;
using NotTable = tbb::concurrent_vector<Level2>;

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

  auto& tbl = impl().notTable;

  // 1) find Op
  if ((size_t)k.op < tbl.size()) {
    auto& lvl2 = tbl.at((size_t)k.op);  // Level2
    if ((size_t)k.varId < lvl2.size()) {
      auto& lvl3 = lvl2.at((size_t)k.varId);  // Level3
      if ((size_t)lid < lvl3.size()) {
        auto& lvl4 = lvl3.at((size_t)lid);  // Level4
        if ((size_t)rid < lvl4.size()) {
          // it4->second is std::shared_ptr<BoolExpr>
          if (lvl4.at((size_t)rid).get() == nullptr) {
            // Replace entry with new node
            std::shared_ptr<BoolExpr> L =
                k.l ? k.l->shared_from_this() : nullptr;
            std::shared_ptr<BoolExpr> R =
                k.r ? k.r->shared_from_this() : nullptr;
            lvl4.at((size_t)rid) = std::shared_ptr<BoolExpr>(
                new BoolExpr(k.op, k.varId, std::move(L), std::move(R)));
          }
          return lvl4.at((size_t)rid);
        }
      }
    }
  }
  std::shared_ptr<BoolExpr> L = k.l ? k.l->shared_from_this() : nullptr;
  std::shared_ptr<BoolExpr> R = k.r ? k.r->shared_from_this() : nullptr;
  for (size_t i = tbl.size(); i <= static_cast<size_t>(k.op); ++i) {
    tbl.push_back(Level2());
  }
  auto& lvl2 = tbl.at((size_t) k.op);
  for (size_t i = lvl2.size(); i <= k.varId; ++i) {
    lvl2.push_back(Level3());
  }
  auto& lvl3 = lvl2.at(k.varId);
  for (size_t i = lvl3.size(); i <= lid; ++i) {
    lvl3.push_back(Level4());
  }
  auto& lvl4 = lvl3.at(lid);
  for (size_t i = lvl4.size(); i <= rid; ++i) {
    lvl4.push_back(nullptr);
  }
  // not found: create via the declared factory that takes Key
  // createNode is declared as: static std::shared_ptr<BoolExpr>
  // createNode(BoolExprCache::Key const& k); ensure BoolExpr grants friendship
  // to BoolExprCache or createNode is public
  auto ptr = std::shared_ptr<BoolExpr>(
      new BoolExpr(k.op, k.varId, std::move(L), std::move(R)));

  // insert into nested maps (operator[] will create missing intermediate maps)
  impl().notTable.at((size_t)k.op).at(k.varId).at(lid).at(rid) = ptr;

  return ptr;
}

}  // namespace KEPLER_FORMAL
