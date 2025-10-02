#pragma once

#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <mutex>
#include "tbb/concurrent_unordered_map.h"
#include "BoolExprCache.h"
#include <cassert>

namespace KEPLER_FORMAL {


/// A hash-consed Boolean expression DAG with eager constant-folding,
/// now protected for concurrent calls.
class BoolExpr : public std::enable_shared_from_this<BoolExpr> {
    friend class BoolExprCache;
public:
    // Convenient constants
    static std::shared_ptr<BoolExpr> createFalse() { return Var(0); }
    static std::shared_ptr<BoolExpr> createTrue()  { return Var(1); }

    // Factory methods (canonical, fold constants, share structure)
    static std::shared_ptr<BoolExpr> Var(size_t id);
    static std::shared_ptr<BoolExpr> Not(std::shared_ptr<BoolExpr> a);
    static std::shared_ptr<BoolExpr> And(std::shared_ptr<BoolExpr> a,
                                         std::shared_ptr<BoolExpr> b);
    static std::shared_ptr<BoolExpr> Or(std::shared_ptr<BoolExpr> a,
                                        std::shared_ptr<BoolExpr> b);
    static std::shared_ptr<BoolExpr> Xor(std::shared_ptr<BoolExpr> a,
                                         std::shared_ptr<BoolExpr> b);

    // Print and stringify
    void Print(std::ostream& out) const;
    std::string toString() const;

    // Evaluate under a map from var-ID → bool
    bool evaluate(const std::unordered_map<size_t,bool>& env) const;

    // Accessors
    Op getOp()    const { return op_; }
    size_t getId() const { return varID_; }
    const std::shared_ptr<BoolExpr>& getLeft()  const { return left_; }
    const std::shared_ptr<BoolExpr>& getRight() const { return right_; }
        std::string getName() const {
        if (op_ != Op::VAR)
            throw std::logic_error("getName: not a variable");
        if (varID_ == 0)
            return "FALSE";
        if (varID_ == 1)
            return "TRUE";
        return "x" + std::to_string(varID_);
    }
    // default constructor
    BoolExpr() = default;
    void setIndex(size_t idx) { index_ = idx; }
    size_t getIndex() const { assert(index_ != (size_t) -1); return index_; }
private:
    // Private ctor: use factory methods
    BoolExpr(Op op, size_t id,
             std::shared_ptr<BoolExpr> l,
             std::shared_ptr<BoolExpr> r);

    Op     op_;
    size_t varID_;
    std::shared_ptr<BoolExpr> left_, right_;
    size_t index_ = (size_t) -1;

    static std::string OpToString(Op);

    static constexpr uint64_t HASH_SEED = 0x9e3779b97f4a7c15ULL;

    static inline uint64_t splitmix64(uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
    }
    struct KeyHash {
        size_t operator()(BoolExprCache::Key const& k) const noexcept {
            const uint64_t s = HASH_SEED;
            uint64_t a = splitmix64(uint64_t(std::hash<int>()(int(k.op))) + s);
            uint64_t b = splitmix64(uint64_t(std::hash<size_t>()(k.varId)) ^ (s<<1));
            uint64_t c = splitmix64(uint64_t(std::uintptr_t(k.l)) + (s>>1));
            uint64_t d = splitmix64(uint64_t(std::uintptr_t(k.r)) ^ (s<<2));
            uint64_t acc = a;
            acc = splitmix64(acc ^ (b + 0x9e3779b97f4a7c15ULL + (acc<<6) + (acc>>2)));
            acc ^= (c + 0x9e3779b97f4a7c15ULL + (acc<<6) + (acc>>2));
            acc = splitmix64(acc ^ d);
            return static_cast<size_t>(acc);
        }
    };
    struct KeyEq {
        bool operator()(BoolExprCache::Key const& a, BoolExprCache::Key const& b) const noexcept {
            return a.op    == b.op
                && a.varId == b.varId
                && a.l     == b.l
                && a.r     == b.r;
        }
    };

    // Global weak-map: BoolExprCache::Key → shared instance
    // guarded by tableMutex_ on every access
    static tbb::concurrent_unordered_map<BoolExprCache::Key,
                              std::weak_ptr<BoolExpr>,
                              KeyHash,
                              KeyEq> table_;

    // Interning constructor (caller must lock tableMutex_)
    static std::shared_ptr<BoolExpr> createNode(BoolExprCache::Key const& k);
};

} // namespace KEPLER_FORMAL
