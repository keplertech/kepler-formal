#include "CompressedBoolExpr.h"
#include <stdexcept>
// for function
#include <functional>

using namespace KEPLER_FORMAL;

// A 4‐tuple signature: (op, varId, leftIdx, rightIdx)
using Sig = std::tuple<Op, size_t, size_t, size_t>;

struct SigHash {
    size_t operator()(Sig const& s) const noexcept {
        auto h0 = std::hash<int>()(static_cast<int>(std::get<0>(s)));
        auto h1 = std::hash<size_t>()(std::get<1>(s));
        auto h2 = std::hash<size_t>()(std::get<2>(s));
        auto h3 = std::hash<size_t>()(std::get<3>(s));
        // combine arbitrarily
        return h0 ^ (h1 << 1) ^ (h2 << 2) ^ (h3 << 3);
    }
};

struct SigEq {
    bool operator()(Sig const& a, Sig const& b) const noexcept {
        return a == b;
    }
};

CompressedBoolExpr
CompressedBoolExpr::compress(BoolExpr* expr) {
    CompressedBoolExpr out;

    // map a BoolExpr* to its node‐index in nodes_
    std::unordered_map<const BoolExpr*, size_t> visited;
    visited.reserve(1024);

    // intern table: signature → unique node index
    std::unordered_map<Sig, size_t, SigHash, SigEq> intern;
    intern.reserve(1024);

    // recursive post‐order
    std::function<size_t(BoolExpr*)> dfs;
    dfs = [&](auto const& bexpr) -> size_t {
        const BoolExpr* raw = bexpr;

        // already processed?
        auto itv = visited.find(raw);
        if (itv != visited.end())
            return itv->second;

        // compute child indices
        size_t leftIdx = 0, rightIdx = 0;
        Op op = bexpr->getOp();

        if (op == Op::NOT) {
            leftIdx = dfs(bexpr->getLeft());
        }
        else if (op == Op::AND || op == Op::OR || op == Op::XOR) {
            leftIdx  = dfs(bexpr->getLeft());
            rightIdx = dfs(bexpr->getRight());
        }
        // VAR: varId = bexpr->getId()

        size_t varId = (op == Op::VAR ? bexpr->getId() : 0);
        Sig sig{op, varId, leftIdx, rightIdx};

        // lookup or create
        auto iti = intern.find(sig);
        size_t myIndex;
        if (iti != intern.end()) {
            myIndex = iti->second;
        } else {
            myIndex = out.nodes_.size();
            intern.emplace(sig, myIndex);
            out.nodes_.push_back({op, varId, leftIdx, rightIdx});
        }

        visited.emplace(raw, myIndex);
        return myIndex;
    };

    out.root_    = dfs(expr);
    return out;
}

void CompressedBoolExpr::save(std::ostream& out) const {
    // header: number of nodes and the root index
    out << nodes_.size() << ' ' << root_ << '\n';
    // each node: op, varId, left, right
    for (auto const& n : nodes_) {
        out << static_cast<int>(n.op) << ' '
            << n.varId  << ' '
            << n.left   << ' '
            << n.right  << '\n';
    }
}

CompressedBoolExpr CompressedBoolExpr::load(std::istream& in) {
    CompressedBoolExpr out;
    size_t N;
    in >> N >> out.root_;
    out.nodes_.resize(N);

    for (size_t i = 0; i < N; ++i) {
        int opi; 
        in >> opi 
           >> out.nodes_[i].varId
           >> out.nodes_[i].left
           >> out.nodes_[i].right;

        if (!in) 
            throw std::runtime_error("CompressedBoolExpr::load: parse error");

        out.nodes_[i].op = static_cast<Op>(opi);
    }
    return out;
}