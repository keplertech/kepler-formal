#pragma once

#include <memory>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <functional>
#include <mutex>
#include "tbb/concurrent_unordered_map.h"

namespace KEPLER_FORMAL {

enum class Op { VAR, AND, OR, NOT, XOR };

/// A hash-consed Boolean expression DAG with eager constant-folding,
/// now protected for concurrent calls.
class BoolExpr : public std::enable_shared_from_this<BoolExpr> {
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
    size_t getId() const { return id_; }
    const std::shared_ptr<BoolExpr>& getLeft()  const { return left_; }
    const std::shared_ptr<BoolExpr>& getRight() const { return right_; }
        std::string getName() const {
        if (op_ != Op::VAR)
            throw std::logic_error("getName: not a variable");
        if (id_ == 0)
            return "FALSE";
        if (id_ == 1)
            return "TRUE";
        return "x" + std::to_string(id_);
    }
    // default constructor
    BoolExpr() = default;
private:
    // Private ctor: use factory methods
    BoolExpr(Op op, size_t id,
             std::shared_ptr<BoolExpr> l,
             std::shared_ptr<BoolExpr> r);

    Op     op_;
    size_t id_;
    std::shared_ptr<BoolExpr> left_, right_;

    static std::string OpToString(Op);

    // Interning key
    struct Key {
        Op        op;
        size_t    varId;
        BoolExpr* l;
        BoolExpr* r;
    };

    struct KeyHash {
        size_t operator()(Key const& k) const noexcept {
            auto h0 = std::hash<int>()(int(k.op));
            auto h1 = std::hash<size_t>()(k.varId);
            auto h2 = std::hash<const void*>()((const void*)k.l);
            auto h3 = std::hash<const void*>()((const void*)k.r);
            return h0 ^ (h1<<1) ^ (h2<<2) ^ (h3<<3);
        }
    };

    struct KeyEq {
        bool operator()(Key const& a, Key const& b) const noexcept {
            return a.op    == b.op
                && a.varId == b.varId
                && a.l     == b.l
                && a.r     == b.r;
        }
    };

    // Global weak-map: Key → shared instance
    // guarded by tableMutex_ on every access
    static tbb::concurrent_unordered_map<Key,
                              std::weak_ptr<BoolExpr>,
                              KeyHash,
                              KeyEq> table_;

    // Interning constructor (caller must lock tableMutex_)
    static std::shared_ptr<BoolExpr> createNode(Key const& k);
};

} // namespace KEPLER_FORMAL
