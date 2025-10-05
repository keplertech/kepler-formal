#pragma once

#include "BoolExpr.h"
#include <vector>
#include <tuple>
#include <unordered_map>
#include <ostream>
#include <istream>

namespace KEPLER_FORMAL {

/// A flat, de‐duplicated representation of a BoolExpr
/// that is ideal for serialization.
///
/// After compression, identical subexpressions share one Node.
/// The tree pointers vanish; only integer indices remain.
struct CompressedBoolExpr {
    /// A single node in the compressed DAG.
    ///  - op == VAR:  varId holds the variable ID, left/right are unused (0).
    ///  - op == NOT:  left holds the child index, right is unused.
    ///  - op ∈ {AND, OR, XOR}: left/right hold child indices.
    struct Node {
        Op     op;
        size_t varId;
        size_t left;
        size_t right;
    };

    /// All unique nodes, indexed [0 .. nodes_.size()-1]
    std::vector<Node> nodes_;

    /// Index of the root node
    size_t root_ = 0;

    /// Build a compressed form from a shared_ptr<BoolExpr> tree.
    static CompressedBoolExpr compress(
        BoolExpr* expr);

    /// Serialize in a simple text format:
    ///   <N> <root>
    ///   <op> <varId> <left> <right>    (repeated N lines)
    void save(std::ostream& out) const;

    /// Load the format written by save().
    static CompressedBoolExpr load(std::istream& in);
};

} // namespace KEPLER_FORMAL
