//===----------------------------------------------------------------------===//
// SNLTruthTableTree.cpp
//===----------------------------------------------------------------------===//

#include "SNLTruthTableTree.h"
#include <stdexcept>    // for std::invalid_argument, std::out_of_range
#include <algorithm>    // for std::sort
#include <utility>      // for std::move

using namespace naja::NL;

namespace KEPLER_FORMAL {

//
// InputNode
//

InputNode::InputNode(size_t idx)
  : inputIndex(idx)
{}

bool InputNode::eval(const std::vector<bool>& extInputs) const
{
    if (inputIndex >= extInputs.size()) {
        throw std::out_of_range("InputNode: external‐input index out of range");
    }
    return extInputs[inputIndex];
}

//
// TableNode
//

TableNode::TableNode(SNLTruthTable t)
  : table(std::move(t))
{}

void TableNode::addChild(std::unique_ptr<ITTNode> child)
{
    children.push_back(std::move(child));
}

bool TableNode::eval(const std::vector<bool>& extInputs) const
{
    uint32_t arity = table.size();
    if (children.size() != arity) {
        throw std::logic_error("TableNode: children count != table.size()");
    }

    std::vector<bool> args;
    args.reserve(arity);
    for (auto const& c : children) {
        args.push_back(c->eval(extInputs));
    }

    // build index into the truth‐table
    uint32_t idx = 0;
    for (uint32_t i = 0; i < arity; ++i) {
        if (args[i]) idx |= (1u << i);
    }
    return table.bits().bit(idx);
}

//
// SNLTruthTableTree
//

SNLTruthTableTree::SNLTruthTableTree()
  : root_(nullptr)
  , numExternalInputs_(0)
{}

SNLTruthTableTree::SNLTruthTableTree(SNLTruthTable table) {
    // determine arity of the table
    uint32_t arity = table.size();

    // create the root TableNode
    auto tblNode = std::make_unique<TableNode>(std::move(table));

    // attach InputNode(0)…InputNode(arity-1)
    for (uint32_t i = 0; i < arity; ++i) {
        tblNode->addChild(std::make_unique<InputNode>(i));
    }

    // install root and record external-input count
    root_              = std::move(tblNode);
    numExternalInputs_ = arity;

    // prepare for future concat()/concatFull()
    updateBorderLeaves();
}

SNLTruthTableTree::SNLTruthTableTree(std::unique_ptr<ITTNode> root,
                                     size_t numExternalInputs)
  : root_(std::move(root))
  , numExternalInputs_(numExternalInputs)
{
    updateBorderLeaves();
}

size_t SNLTruthTableTree::size() const
{
    return numExternalInputs_;
}

bool SNLTruthTableTree::eval(const std::vector<bool>& extInputs) const
{
    if (extInputs.size() != numExternalInputs_) {
        throw std::invalid_argument(
            "SNLTruthTableTree::eval: wrong input vector size");
    }
    return root_->eval(extInputs);
}

void SNLTruthTableTree::updateBorderLeaves()
{
    borderLeaves_.clear();

    struct Frame {
        TableNode* parent;
        size_t      childPos;
        ITTNode*    node;
    };

    std::vector<Frame> stk;
    stk.reserve(64);
    stk.push_back({ nullptr, 0, root_.get() });

    while (!stk.empty()) {
        auto f = stk.back();
        stk.pop_back();

        if (auto inp = dynamic_cast<InputNode*>(f.node)) {
            borderLeaves_.push_back(
              { f.parent, f.childPos, inp->inputIndex }
            );
        }
        else if (auto tbl = dynamic_cast<TableNode*>(f.node)) {
            for (size_t i = tbl->children.size(); i-- > 0; ) {
                stk.push_back({ tbl, i, tbl->children[i].get() });
            }
        }
    }

    std::sort(borderLeaves_.begin(),
              borderLeaves_.end(),
              [](auto const& a, auto const& b) {
                  return a.extIndex < b.extIndex;
              });
}

void SNLTruthTableTree::concatBody(size_t borderIndex, naja::NL::SNLTruthTable table) {
    uint32_t arity = table.size();
    // 2) treat borderIndex as a slot index or ext‐input index
    size_t idx = borderIndex;
    if (idx >= borderLeaves_.size()) {
        bool found = false;
        for (size_t i = 0; i < borderLeaves_.size(); ++i) {
            if (borderLeaves_[i].extIndex == borderIndex) {
                idx = i;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::out_of_range("concat: leafIndex out of range");
        }
    }

    // --- original splice body ---
    auto leaf = borderLeaves_[idx];

    std::unique_ptr<ITTNode> oldLeaf;
    if (leaf.parent) {
        oldLeaf = std::move(leaf.parent->children[leaf.childPos]);
    } else {
        oldLeaf = std::move(root_);
    }

    auto newNode = std::make_unique<TableNode>(std::move(table));
    for (uint32_t i = 0; i < arity; ++i) {
        if (i == 0) {
            newNode->addChild(std::move(oldLeaf));
            continue;
        }
        newNode->addChild(
            std::make_unique<InputNode>(numExternalInputs_ + (i - 1))
        );
    }

    if (leaf.parent) {
        leaf.parent->children[leaf.childPos] = std::move(newNode);
    } else {
        root_ = std::move(newNode);
    }

    
    // --- end splice body ---
}

void SNLTruthTableTree::concat(size_t borderIndex, SNLTruthTable table)
{
    uint32_t arity = table.size();
    // 1) refresh leaf‐slot info
    updateBorderLeaves();

    concatBody(borderIndex, std::move(table));   
    numExternalInputs_ += (arity - 1);
    // 3) rebuild slots for future calls
    updateBorderLeaves();
}

/// Splice a batch of truth‐tables onto the first N border leaves.
/// tables[i] is attached to slot i.  Any leftover leaves remain.
/// Processing is done in reverse order so that earlier splices do not
/// invalidate higher‐index slots.
void SNLTruthTableTree::concatFull(std::vector<SNLTruthTable> tables)
{   
    size_t newNumInptuts = numExternalInputs_;
    //printf("--- concatFull: current numExternalInputs = %zu\n", newNumInptuts);
    updateBorderLeaves();
    size_t n = tables.size();
    if (n > borderLeaves_.size()) {
        throw std::invalid_argument(
            "concatFull: more tables than border leaves"
        );
    }
    for (size_t i = 0; i < n; ++i) {
        concatBody(i, std::move(tables[i]));
        
        //printf("Truth table: %s\n", tables[i].getString().c_str());
        //printf("tables[%zu].size() = %zu\n", i, tables[i].size());
        ////printf("tables[%zu].size() - 1 = %zu\n", i, tables[i].size() - 1);
        newNumInptuts =  newNumInptuts + tables[i].size() - 1;
        //printf("After concatBody %zu: current numExternalInputs = %zu\n", i, newNumInptuts);
    }
    numExternalInputs_ = newNumInptuts;
    //printf("--- concatFull: current numExternalInputs = %zu\n", newNumInptuts);
    updateBorderLeaves();
}

bool SNLTruthTableTree::isInitialized() const
{
    // Iterate through the tree and check for each node that the tt accosiated is initialized
    if (!root_) return false;
    struct Frame {
        ITTNode* node;
    };
    std::vector<Frame> stk;
    stk.push_back({ root_.get() });
    while (!stk.empty()) {
        auto f = stk.back();
        stk.pop_back();

        if (auto tbl = dynamic_cast<TableNode*>(f.node)) {
            if (!tbl->table.isInitialized()) return false;
            for (size_t i = 0; i < tbl->children.size(); ++i) {
                stk.push_back({ tbl->children[i].get() });
            }
        }
    }
    return true;
}

void SNLTruthTableTree::print() const
{
    if (!root_) {
        //printf("<empty tree>\n");
        return;
    }

    struct Frame {
        ITTNode* node;
        size_t   depth;
    };

    std::vector<Frame> stk;
    stk.push_back({ root_.get(), 0 });

    while (!stk.empty()) {
        auto f = stk.back();
        stk.pop_back();

        for (size_t i = 0; i < f.depth; ++i) {
            //printf("  ");
        }

        if (auto inp = dynamic_cast<InputNode*>(f.node)) {
            //printf("InputNode(%zu)\n", inp->inputIndex);
        }
        else if (auto tbl = dynamic_cast<TableNode*>(f.node)) {
            //printf("TableNode(size=%u, bits=%s)\n",
                   //tbl->table.size(),
                   //tbl->table.getString().c_str());
            for (size_t i = tbl->children.size(); i-- > 0; ) {
                stk.push_back({ tbl->children[i].get(), f.depth + 1 });
            }
        }
        else {
            //printf("Unknown ITTNode subtype\n");
        }
    }
}

} // namespace naja::NL
