#include "BoolExpr.h"
#include <cassert>
#include <unordered_map>

namespace KEPLER_FORMAL {

// static definitions
tbb::concurrent_unordered_map<BoolExprCache::Key,
                   std::weak_ptr<BoolExpr>,
                   BoolExpr::KeyHash,
                   BoolExpr::KeyEq>
    BoolExpr::table_{};

/// Private ctor
BoolExpr::BoolExpr(Op op, size_t id,
                   BoolExpr* a,
                   BoolExpr* b)
  : op_(op), varID_(id)/*, left_(l) , right_(r)*/ {
    if (b == nullptr) {
        if (a == nullptr && op != Op::VAR) {
            throw std::invalid_argument("BoolExpr: non-VAR with null children");
        }
        right_ = nullptr;
        left_  = a;
    } else if (*b <= *a) {
        left_  = b;
        right_ = a;
    } else {
        left_  = a;
        right_ = b;
    }
}

/// Intern+construct a new node if needed
BoolExpr*
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
    BoolExpr* L = k.l ? k.l->shared_from_this() : nullptr;
    BoolExpr* R = k.r ? k.r->shared_from_this() : nullptr;

    auto ptr = BoolExpr*(
        new BoolExpr(k.op, k.varId, std::move(L), std::move(R))
    );
    table_.emplace(k, ptr);
    return ptr;*/
    return BoolExprCache::getExpression(k);
}

// Factory methods with eager folding & sharing

BoolExpr* BoolExpr::Var(size_t id) {
    BoolExprCache::Key k{Op::VAR, id, nullptr, nullptr};
    return createNode(k);
}

BoolExpr* BoolExpr::Not(BoolExpr* a) {
    // constant-fold
    if (a->op_ == Op::VAR && a->varID_ < 2)
        return Var(1 - a->varID_);
    // double negation
    if (a->op_ == Op::NOT)
        return a->left_;
    BoolExprCache::Key k{Op::NOT, 0, a, nullptr};
    return createNode(k);
}

BoolExpr* BoolExpr::And(
    BoolExpr* a,
    BoolExpr* b)
{
    // constant-fold
    if ((a->op_ == Op::VAR && a->varID_ == 0) ||
        (b->op_ == Op::VAR && b->varID_ == 0))
        return Var(0);
    if (a->op_ == Op::VAR && a->varID_ == 1) return b;
    if (b->op_ == Op::VAR && b->varID_ == 1) return a;
    if (a == b)              return a;
    if (a->op_==Op::NOT && a->left_==b) return Var(0);
    if (b->op_==Op::NOT && b->left_==a) return Var(0);

    // canonical order
    if (b < a) std::swap(a, b);
    BoolExprCache::Key k{Op::AND, 0, a, b};
    return createNode(k);
}

BoolExpr* BoolExpr::Or(
    BoolExpr* a,
    BoolExpr* b)
{
    if ((a->op_ == Op::VAR && a->varID_ == 1) ||
        (b->op_ == Op::VAR && b->varID_ == 1))
        return Var(1);
    if (a->op_ == Op::VAR && a->varID_ == 0) return b;
    if (b->op_ == Op::VAR && b->varID_ == 0) return a;
    if (a == b)             return a;
    if (a->op_==Op::NOT && a->left_==b) return Var(1);
    if (b->op_==Op::NOT && b->left_==a) return Var(1);

    if (b < a) std::swap(a, b);
    BoolExprCache::Key k{Op::OR, 0, a, b};
    return createNode(k);
}

BoolExpr* BoolExpr::Xor(
    BoolExpr* a,
    BoolExpr* b)
{
    if (a->op_ == Op::VAR && a->varID_ == 0)     return b;
    if (b->op_ == Op::VAR && b->varID_ == 0)     return a;
    if (a->op_ == Op::VAR && a->varID_ == 1)     return Not(b);
    if (b->op_ == Op::VAR && b->varID_ == 1)     return Not(a);
    if (a == b)                  return Var(0);

    if (b < a) std::swap(a, b);
    BoolExprCache::Key k{Op::XOR, 0, a, b};
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

// replace previous isConstFalse/isConstTrue and Simplify implementation with this:

static inline bool isConstFalse(BoolExpr* e) {
    return e->getOp() == Op::VAR && e->getId() == 0;
}
static inline bool isConstTrue(BoolExpr* e) {
    return e->getOp() == Op::VAR && e->getId() == 1;
}

BoolExpr* BoolExpr::simplify(BoolExpr* e) {
    if (!e) return nullptr;
    if (e->getOp() == Op::VAR) return e;

    std::unordered_map<BoolExpr*, BoolExpr*> memo;
    std::vector<BoolExpr*> stack;
    std::unordered_map<BoolExpr*, int> state;
    std::vector<BoolExpr*> order;

    stack.push_back(e);
    while (!stack.empty()) {
        BoolExpr* n = stack.back();
        stack.pop_back();
        auto itst = state.find(n);
        if (itst == state.end()) {
            state[n] = 1;
            stack.push_back(n);
            if (n->getRight()) stack.push_back(n->getRight());
            if (n->getLeft())  stack.push_back(n->getLeft());
        } else {
            order.push_back(n);
        }
    }

    for (BoolExpr* node : order) {
        switch (node->getOp()) {
        case Op::NOT: {
            BoolExpr* a = memo.count(node->getLeft()) ? memo[node->getLeft()] : node->getLeft();
            if (isConstFalse(a)) { memo[node] = Var(1); break; }
            if (isConstTrue(a))  { memo[node] = Var(0); break; }
            if (a->getOp() == Op::NOT) { memo[node] = a->getLeft(); break; }
            memo[node] = Not(a);
            break;
        }
        case Op::AND: {
            BoolExpr* A = memo.count(node->getLeft()) ? memo[node->getLeft()] : node->getLeft();
            BoolExpr* B = memo.count(node->getRight()) ? memo[node->getRight()] : node->getRight();

            if (isConstFalse(A) || isConstFalse(B)) { memo[node] = Var(0); break; }
            if (isConstTrue(A)) { memo[node] = B; break; }
            if (isConstTrue(B)) { memo[node] = A; break; }
            if (A == B) { memo[node] = A; break; }
            if ((A->getOp() == Op::NOT && A->getLeft() == B) ||
                (B->getOp() == Op::NOT && B->getLeft() == A)) {
                memo[node] = Var(0); break;
            }
            // Factory will canonicalize order
            memo[node] = And(A, B);
            break;
        }
        case Op::OR: {
            BoolExpr* A = memo.count(node->getLeft()) ? memo[node->getLeft()] : node->getLeft();
            BoolExpr* B = memo.count(node->getRight()) ? memo[node->getRight()] : node->getRight();

            if (isConstTrue(A) || isConstTrue(B)) { memo[node] = Var(1); break; }
            if (isConstFalse(A)) { memo[node] = B; break; }
            if (isConstFalse(B)) { memo[node] = A; break; }
            if (A == B) { memo[node] = A; break; }
            if ((A->getOp() == Op::NOT && A->getLeft() == B) ||
                (B->getOp() == Op::NOT && B->getLeft() == A)) {
                memo[node] = Var(1); break;
            }
            memo[node] = Or(A, B);
            break;
        }
        case Op::XOR: {
            BoolExpr* A = memo.count(node->getLeft()) ? memo[node->getLeft()] : node->getLeft();
            BoolExpr* B = memo.count(node->getRight()) ? memo[node->getRight()] : node->getRight();

            if (isConstFalse(A)) { memo[node] = B; break; }
            if (isConstFalse(B)) { memo[node] = A; break; }
            if (isConstTrue(A))  { memo[node] = Not(B); break; }
            if (isConstTrue(B))  { memo[node] = Not(A); break; }
            if (A == B) { memo[node] = Var(0); break; }
            memo[node] = Xor(A, B);
            break;
        }
        default:
            memo[node] = node;
            break;
        }
    }

    return memo.count(e) ? memo[e] : e;
}



} // namespace KEPLER_FORMAL
