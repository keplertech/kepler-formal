#include "BoolExpr.h"
#include <cassert>

namespace KEPLER_FORMAL {

// static definitions
tbb::concurrent_unordered_map<BoolExprCache::Key,
                   std::weak_ptr<BoolExpr>,
                   BoolExpr::KeyHash,
                   BoolExpr::KeyEq>
    BoolExpr::table_{};

/// Private ctor
BoolExpr::BoolExpr(Op op, size_t id,
                   std::shared_ptr<BoolExpr> l,
                   std::shared_ptr<BoolExpr> r)
  : op_(op), varID_(id), left_(std::move(l)), right_(std::move(r))
{}

/// Intern+construct a new node if needed
std::shared_ptr<BoolExpr>
BoolExpr::createNode(BoolExprCache::Key const& k) {
    // Caller already holds lock on tableMutex_
    // print the size in GB of table_
    //printf("BoolExpr table size: %.2f GB\n",
    //       (table_.size() * (sizeof(Key) + sizeof(std::weak_ptr<BoolExpr>))) / (1024.0*1024.0*1024.0));

    /*auto it = table_.find(k);
    if (it != table_.end()) {
        if (auto existing = it->second.lock())
            return existing;
    }
    // retrieve shared_ptr to children via non-const shared_from_this()
    std::shared_ptr<BoolExpr> L = k.l ? k.l->shared_from_this() : nullptr;
    std::shared_ptr<BoolExpr> R = k.r ? k.r->shared_from_this() : nullptr;

    auto ptr = std::shared_ptr<BoolExpr>(
        new BoolExpr(k.op, k.varId, std::move(L), std::move(R))
    );
    table_.emplace(k, ptr);
    return ptr;*/
    return BoolExprCache::getExpression(k);
}

// Factory methods with eager folding & sharing

std::shared_ptr<BoolExpr> BoolExpr::Var(size_t id) {
    BoolExprCache::Key k{Op::VAR, id, nullptr, nullptr};
    return createNode(k);
}

std::shared_ptr<BoolExpr> BoolExpr::Not(std::shared_ptr<BoolExpr> a) {
    // constant-fold
    if (a->op_ == Op::VAR && a->varID_ < 2)
        return Var(1 - a->varID_);
    // double negation
    if (a->op_ == Op::NOT)
        return a->left_;
    BoolExprCache::Key k{Op::NOT, 0, a.get(), nullptr};
    return createNode(k);
}

std::shared_ptr<BoolExpr> BoolExpr::And(
    std::shared_ptr<BoolExpr> a,
    std::shared_ptr<BoolExpr> b)
{
    // constant-fold
    if ((a->op_ == Op::VAR && a->varID_ == 0) ||
        (b->op_ == Op::VAR && b->varID_ == 0))
        return Var(0);
    if (a->op_ == Op::VAR && a->varID_ == 1) return b;
    if (b->op_ == Op::VAR && b->varID_ == 1) return a;
    if (a.get() == b.get())              return a;
    if (a->op_==Op::NOT && a->left_.get()==b.get()) return Var(0);
    if (b->op_==Op::NOT && b->left_.get()==a.get()) return Var(0);

    // canonical order
    if (b.get() < a.get()) std::swap(a, b);
    BoolExprCache::Key k{Op::AND, 0, a.get(), b.get()};
    return createNode(k);
}

std::shared_ptr<BoolExpr> BoolExpr::Or(
    std::shared_ptr<BoolExpr> a,
    std::shared_ptr<BoolExpr> b)
{
    if ((a->op_ == Op::VAR && a->varID_ == 1) ||
        (b->op_ == Op::VAR && b->varID_ == 1))
        return Var(1);
    if (a->op_ == Op::VAR && a->varID_ == 0) return b;
    if (b->op_ == Op::VAR && b->varID_ == 0) return a;
    if (a.get() == b.get())             return a;
    if (a->op_==Op::NOT && a->left_.get()==b.get()) return Var(1);
    if (b->op_==Op::NOT && b->left_.get()==a.get()) return Var(1);

    if (b.get() < a.get()) std::swap(a, b);
    BoolExprCache::Key k{Op::OR, 0, a.get(), b.get()};
    return createNode(k);
}

std::shared_ptr<BoolExpr> BoolExpr::Xor(
    std::shared_ptr<BoolExpr> a,
    std::shared_ptr<BoolExpr> b)
{
    if (a->op_ == Op::VAR && a->varID_ == 0)     return b;
    if (b->op_ == Op::VAR && b->varID_ == 0)     return a;
    if (a->op_ == Op::VAR && a->varID_ == 1)     return Not(b);
    if (b->op_ == Op::VAR && b->varID_ == 1)     return Not(a);
    if (a.get() == b.get())                  return Var(0);

    if (b.get() < a.get()) std::swap(a, b);
    BoolExprCache::Key k{Op::XOR, 0, a.get(), b.get()};
    return createNode(k);
}

// Print routines unchanged…

void BoolExpr::Print(std::ostream& out) const { 
    switch (op_) {
        case Op::VAR:
            out << varID_;
            break;
        case Op::NOT:
            out << "¬";
            if (left_->op_ != Op::VAR)
                out << "(";
            left_->Print(out);
            if (left_->op_ != Op::VAR)
                out << ")";
            break;
        case Op::AND:
        case Op::OR:
        case Op::XOR:
            if (left_->op_ != Op::VAR)
                out << "(";
            left_->Print(out);
            if (left_->op_ != Op::VAR)
                out << ")";
            out << " " << OpToString(op_) << " ";
            if (right_->op_ != Op::VAR)
                out << "(";
            right_->Print(out);
            if (right_->op_ != Op::VAR)
                out << ")";
            break;
        default:
            assert(false && "unknown BoolExpr op");
    }
}
std::string BoolExpr::toString() const     { 
    // print content to string
    std::ostringstream oss;
    Print(oss);
    return oss.str();
}
//bool BoolExpr::evaluate(const std::unordered_map<size_t,bool>& env) const { /* … */ }
std::string BoolExpr::OpToString(Op op) { 
    switch (op) {
        case Op::VAR: return "VAR";
        case Op::NOT: return "NOT";
        case Op::AND: return "AND";
        case Op::OR:  return "OR";
        case Op::XOR: return "XOR";
        default:      return "UNKNOWN";
    }
}

} // namespace KEPLER_FORMAL
