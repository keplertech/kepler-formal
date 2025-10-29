#include "BoolExprCache.h"
#include <tbb/concurrent_unordered_map.h>
#include <tbb/tbb_allocator.h>
#include <memory>
#include <cassert>
#include <atomic>
#include <tuple>
#include <cstdint>
#include "BoolExpr.h"

namespace KEPLER_FORMAL {

// atomic id counter
std::atomic<size_t> BoolExprCache::lastID_{1};
size_t BoolExprCache::numMiss_ = 0;
size_t BoolExprCache::numQuaries_ = 0;
size_t BoolExprCache::numHit_ = 0;


// Tuple key: (op,varId,lid,rid) using pointer identity for children
using TupleKey = std::tuple<uint32_t, uint64_t, uint64_t, uint64_t>;

// hasher for TupleKey
struct TupleKeyHasher {
  size_t operator()(TupleKey const& t) const noexcept {
    uint64_t a = static_cast<uint64_t>(std::get<0>(t));
    uint64_t b = std::get<1>(t);
    uint64_t c = std::get<2>(t);
    uint64_t d = std::get<3>(t);
    uint64_t x = a;
    x = (x * 11400714819323198485ull) ^ (b + 0x9e3779b97f4a7c15ull + (x<<6) + (x>>2));
    x ^= (c + 0x9e3779b97f4a7c15ull + (x<<6) + (x>>2));
    x ^= (d + 0x9e3779b97f4a7c15ull + (x<<6) + (x>>2));
    return static_cast<size_t>(x ^ (x >> 32));
  }
};

struct TupleKeyEq {
  bool operator()(TupleKey const& a, TupleKey const& b) const noexcept { return a == b; }
};

using ValueT = std::shared_ptr<BoolExpr>;
using PairT = std::pair<const TupleKey, ValueT>;
using TbbAlloc = tbb::tbb_allocator<PairT>;

// concurrent map using tbb allocator
using SingleMap = tbb::concurrent_unordered_map<TupleKey, ValueT, TupleKeyHasher, TupleKeyEq, TbbAlloc>;

struct BoolExprCache::Impl {
  SingleMap table;
};

BoolExprCache::Impl& BoolExprCache::impl() {
  static Impl instance;
  return instance;
}

static inline TupleKey make_tuple_key(Op op, size_t varId, std::shared_ptr<BoolExpr> lptr, std::shared_ptr<BoolExpr> rptr) noexcept {
  // use pointer identity as integer; nullptr -> 0
  auto lid = reinterpret_cast<uint64_t>(lptr.get());
  auto rid = reinterpret_cast<uint64_t>(rptr.get());
  return TupleKey{
    static_cast<uint32_t>(op),
    static_cast<uint64_t>(varId),
    lid,
    rid
  };
}

std::shared_ptr<BoolExpr> BoolExprCache::getExpression(Key const& k) {
  std::shared_ptr<BoolExpr> lptr = k.l;
  std::shared_ptr<BoolExpr> rptr = k.r;
  TupleKey tk = make_tuple_key(k.op, k.varId, lptr, rptr);
  if (k.l == nullptr || k.r == nullptr) {
    if (k.l == nullptr) {
      tk = make_tuple_key(k.op, k.varId, rptr, nullptr);
    } else if (k.r == nullptr) {
      tk = make_tuple_key(k.op, k.varId, lptr, nullptr);
    } else {
      // both null
      tk = make_tuple_key(k.op, k.varId, nullptr, nullptr);
    }
  } else if (*k.l <= *k.r) {
    // enforce canonical order for commutative ops
    tk = make_tuple_key(k.op, k.varId, rptr, lptr);
  } else {
    tk = make_tuple_key(k.op, k.varId, lptr, rptr);
  }

  auto &tbl = impl().table;
  
  numQuaries_ += 1;
  // quick lookup
  auto it = tbl.find(tk);
  if (it != tbl.end()) {
    size_t id = lastID_.fetch_add(1, std::memory_order_relaxed) + 1;
    //it->second->setIndex(id);
    assert(it->second != nullptr);
    numHit_ += 1;
    //printf("######### numHit: %lu\n", numHit_);
    return it->second;
  }

  // construct new BoolExpr. We need shared_ptr owners for children if they exist.
  std::shared_ptr<BoolExpr> L = lptr ? lptr : nullptr;
  std::shared_ptr<BoolExpr> R = rptr ? rptr : nullptr;

  // use new because constructor may be non-public
  std::shared_ptr<BoolExpr> newptr(new BoolExpr(k.op, k.varId, L, R));

  // assign id atomically
  size_t id = lastID_.fetch_add(1, std::memory_order_relaxed) + 1;
  //newptr->setIndex(id);

  // insert; if another thread inserted concurrently, use that one
  auto pr = tbl.insert({tk, newptr});
  // print size of cache
  
  if (!pr.second) {
    //delete newptr;
    // destroy newptr
    
    return pr.first->second;
  }
  numMiss_ += 1;
  //printf("miss rate: %lf\n", (double) numHit_ / (double) numMiss_);
  //printf("size of cache: %lu\n", tbl.size());

  return newptr;
}

void BoolExprCache::destroy() {
  // delete all stored std::shared_ptr<BoolExpr>
  // for (auto& kv : impl().table) {
  //   delete kv.second;
  // }
  impl().table.clear();
}

} // namespace KEPLER_FORMAL
