// Tree2BoolExpr.cpp

#include "Tree2BoolExpr.h"
#include "SNLTruthTableTree.h"     // for SNLTruthTableTree::Node
#include "SNLTruthTable.h"
#include "BoolExpr.h"
#include "DNL.h"

#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <vector>
#include <utility>   // for std::pair
#include <bitset>    // needed for bit-manipulation
#include <tbb/tbb_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_vector.h>

using namespace naja::NL;
using namespace KEPLER_FORMAL;

typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> TermsPair;
tbb::enumerable_thread_specific<TermsPair> termsETS;
tbb::concurrent_vector<TermsPair*> termsETSvector;

void initTermsETS() {
    if (termsETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
        for (size_t i = termsETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
            termsETSvector.push_back(nullptr);
        }
    }
    if (termsETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
        termsETSvector[tbb::this_task_arena::current_thread_index()] = &termsETS.local();
    }
}

TermsPair& getTErmsETS() {
    //initTermsETS();
    return *termsETSvector[tbb::this_task_arena::current_thread_index()];
}

size_t sizeOfTermsETS() {
    return getTErmsETS().second;
}

void clearTermsETS() {
    getTErmsETS().second = 0;
}

void pushBackTermsETS(const std::shared_ptr<BoolExpr>& term) {
    auto& vec = getTErmsETS().first;
    auto& sz = getTErmsETS().second;
    if (vec.size() > sz) {
        vec[sz] = term;
        sz++;
        return;
    }
    vec.push_back(term);
    sz++;
}

void reserveTermsETS(size_t n) {
    if (getTErmsETS().first.size() >= n)
        return;
    getTErmsETS().first.reserve(n);
}

bool emptyTermsETS() {
    return getTErmsETS().second == 0;
}

// same for std::vector<bool, tbb::tbb_allocator<bool>> 
typedef std::pair<std::vector<bool, tbb::tbb_allocator<bool>>, size_t> RelevantPair;
tbb::enumerable_thread_specific<RelevantPair> relevantETS;
tbb::concurrent_vector<RelevantPair*> relevantETSvector;

void initRelevantETS() {
    if (relevantETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
        for (size_t i = relevantETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
            relevantETSvector.push_back(nullptr);
        }
    }
    if (relevantETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
        relevantETSvector[tbb::this_task_arena::current_thread_index()] = &relevantETS.local();
    }
}

RelevantPair& getRelevantETS() {
    //initRelevantETS();
    return *relevantETSvector[tbb::this_task_arena::current_thread_index()];
}

size_t sizeOfRelevantETS() {
    return getRelevantETS().second;
}

void clearRelevantETS() {
    getRelevantETS().second = 0;
}

void pushBackRelevantETS(bool b) {
    auto& vec = getRelevantETS().first;
    auto& sz = getRelevantETS().second;
    if (vec.size() > sz) {
        vec[sz] = b;
        sz++;
        return;
    }
    vec.push_back(b);
    sz++;
}

void setRelevantETS(size_t i, bool b) {
    if (i >= getRelevantETS().second) {
       assert(false && "setRelevantETS: index out of range");
    }
    getRelevantETS().first[i] = b;
}

bool getRelevantETS(size_t i) {
    if (i >= getRelevantETS().second) {
        throw std::out_of_range("getRelevantETS: index out of range");
    }
    return getRelevantETS().first[i];
}

void reserveRelevantETSwithFalse(size_t n) {
    auto& vec = getRelevantETS().first;
    auto& sz = getRelevantETS().second;
    if (vec.size() >= n) {
        vec.assign(n, false);
        sz = n;
        return;
    }
    size_t oldSize = vec.size();
    vec.resize(n, false);
    vec.assign(n, false);
    sz = n;
}

// do same for std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> memo;
typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> MemoPair;
tbb::enumerable_thread_specific<MemoPair> memoETS;
tbb::concurrent_vector<MemoPair*> memoETSvector;

void initMemoETS() {
    if (memoETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
        for (size_t i = memoETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
            memoETSvector.push_back(nullptr);
        }
    }
    if (memoETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
        memoETSvector[tbb::this_task_arena::current_thread_index()] = &memoETS.local();
    }
}

MemoPair& getMemoETS() {
    //initMemoETS();
    return *memoETSvector[tbb::this_task_arena::current_thread_index()];
}

size_t sizeOfMemoETS() {
    return getMemoETS().second;
}

void clearMemoETS() {
    getMemoETS().second = 0;
}

void pushBackMemoETS(const std::shared_ptr<BoolExpr>& expr) {
    auto& vec = getMemoETS().first;
    auto& sz = getMemoETS().second;
    if (vec.size() > sz) {
        vec[sz] = expr;
        sz++;
        return;
    }
    vec.push_back(expr);
    sz++;
}

void reserveMemoETS(size_t n) {
    auto& vec = getMemoETS().first;
    auto& sz = getMemoETS().second;
    if (vec.size() >= n) {
        sz = n;
        vec.assign(n, nullptr);
        return;
    }
    vec.resize(n);
    sz = n;
    vec.assign(n, nullptr);
}

void setMemoETS(size_t i, const std::shared_ptr<BoolExpr>& expr) {
    if (i >= getMemoETS().second) {
       assert(false && "setMemoETS: index out of range");
    }
    getMemoETS().first[i] = expr;
}

const std::shared_ptr<BoolExpr>& getMemoETS(size_t i) {
    if (i >= getMemoETS().second) {
        assert(false && "getMemoETS: index out of range");
    }
    return getMemoETS().first[i];
}

// same for std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> childF;
typedef std::pair<std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>, size_t> ChildFETSPair;
tbb::enumerable_thread_specific<ChildFETSPair> childFETS;
tbb::concurrent_vector<ChildFETSPair*> childFETSvector;

void initChildFETS() {
    if (childFETSvector.size() <= tbb::this_task_arena::current_thread_index()) {
        for (size_t i = childFETSvector.size(); i <= tbb::this_task_arena::current_thread_index(); i++) {
            childFETSvector.push_back(nullptr);
        }
    }
    if (childFETSvector[tbb::this_task_arena::current_thread_index()] == nullptr) {
        childFETSvector[tbb::this_task_arena::current_thread_index()] = &childFETS.local();
    }
}

ChildFETSPair& getChildFETS() {
    //initChildFETS();
    return *childFETSvector[tbb::this_task_arena::current_thread_index()];
}

size_t sizeOfChildFETS() {
    return getChildFETS().second;
}

void clearChildFETS() {
    getChildFETS().second = 0;
}

void pushBackChildFETS(const std::shared_ptr<BoolExpr>& expr) {
    auto& vec = getChildFETS().first;
    auto& sz = getChildFETS().second;
    if (vec.size() > sz) {
        vec[sz] = expr;
        sz++;
        return;
    }
    vec.push_back(expr);
    sz++;
}

void reserveChildFETS(size_t n) {
    auto& vec = getChildFETS().first;
    auto& sz = getChildFETS().second;
    if (vec.size() >= n) {
        sz = n;
        vec.assign(n, nullptr);
        return;
    }
    vec.resize(n);
    sz = n;
    vec.assign(n, nullptr);
}

const std::shared_ptr<BoolExpr>& getChildFETS(size_t i) {
    if (i >= getChildFETS().second) {
        assert(false && "getChildFETS: index out of range");
    }
    return getChildFETS().first[i];
}

void setChildFETS(size_t i, const std::shared_ptr<BoolExpr>& expr) {
    if (i >= getChildFETS().second) {
       assert(false && "setChildFETS: index out of range");
    }
    getChildFETS().first[i] = expr;
}


size_t toSizeT(const std::string& s) {
    if (s.empty()) {
        assert(false && "toSizeT: empty string");
    }

    size_t result = 0;
    const size_t max = std::numeric_limits<size_t>::max();

    for (unsigned char uc : s) {
        if (!std::isdigit(uc)) {
            throw std::invalid_argument(
                "toSizeT: invalid character '" + std::string(1, static_cast<char>(uc)) + "' in input"
            );
        }
        size_t digit = static_cast<size_t>(uc - '0');

        // Check for overflow: result * 10 + digit > max
        if (result > (max - digit) / 10) {
            throw std::out_of_range("toSizeT: value out of range for size_t");
        }

        result = result * 10 + digit;
    }

    return result;
}

// Fold a list of literals into a single AND
static std::shared_ptr<BoolExpr>
mkAnd(const std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>& lits)
{
    if (lits.empty())
        return BoolExpr::createTrue();
    auto cur = lits[0];
    for (size_t i = 1; i < lits.size(); ++i)
        cur = BoolExpr::And(cur, lits[i]);
    return cur;
}

// Fold a list of terms into a single OR
static std::shared_ptr<BoolExpr>
mkOr(const std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>>& terms)
{
    if (terms.empty())
        return BoolExpr::createFalse();
    auto cur = terms[0];
    for (size_t i = 1; i < terms.size(); ++i)
        cur = BoolExpr::Or(cur, terms[i]);
    return cur;
}

std::shared_ptr<BoolExpr>
Tree2BoolExpr::convert(
    const SNLTruthTableTree&                        tree,
    const std::vector<size_t>& varNames)
{
    initChildFETS();
    initMemoETS();
    initRelevantETS();
    initTermsETS();
    const auto* root = tree.getRoot();
    if (!root)
        return nullptr;

    // 1) find maxID
    size_t maxID = 0;
    {
        using NodePtr = const SNLTruthTableTree::Node*;
        std::vector<NodePtr, tbb::tbb_allocator<NodePtr>> dfs;
        dfs.reserve(128);
        dfs.push_back(root);
        while (!dfs.empty()) {
            NodePtr n = dfs.back(); dfs.pop_back();
            maxID = std::max(maxID, n->nodeID);
            if (n->type == SNLTruthTableTree::Node::Type::Table ||
                n->type == SNLTruthTableTree::Node::Type::P)
            {
                // iterate children by const ref to avoid copies
                for (const auto& cptr : n->children)
                    dfs.push_back(cptr.get());
            }
        }
    }

    // 2) memo table
    //std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> memo;
    //memo.resize(maxID + 1);
    clearMemoETS();
    reserveMemoETS(maxID + 1);
    // 3) post-order build
    using Frame = std::pair<const SNLTruthTableTree::Node*, bool>;
    std::vector<Frame, tbb::tbb_allocator<Frame>> stack;
    stack.reserve(maxID + 1);
    stack.emplace_back(root, false);

    while (!stack.empty()) {
        Frame f = stack.back(); stack.pop_back();
        const SNLTruthTableTree::Node* node = f.first;
        bool visited = f.second;
        size_t id = node->nodeID;

        if (!visited) {
            //if (memo[id])
            if (getMemoETS(id).get() != nullptr)
                continue;

            if (node->type == SNLTruthTableTree::Node::Type::Table ||
                node->type == SNLTruthTableTree::Node::Type::P)
            {
                stack.emplace_back(node, true);
                for (const auto& c : node->children)
                    stack.emplace_back(c.get(), false);
            }
            else {
                assert(node->type == SNLTruthTableTree::Node::Type::Input);
                auto parent = node->parent.lock();
                assert(parent && parent->type == SNLTruthTableTree::Node::Type::P);
                if (parent->termid >= varNames.size()) {
                    printf("varNames size: %zu, parent termid: %u\n", varNames.size(), parent->termid);
                    assert(parent->termid < varNames.size());
                }
                //memo[id] = BoolExpr::Var(varNames[parent->termid]);
                setMemoETS(id, BoolExpr::Var(varNames[parent->termid]));
            }
        } else {
            // post-visit for Table / P
            const SNLTruthTable& tbl = node->getTruthTable();
            uint32_t k = tbl.size();
            uint64_t rows = uint64_t{1} << k;

            if (tbl.all0()) {
                //memo[id] = BoolExpr::createFalse();
                setMemoETS(id, BoolExpr::createFalse());
            }
            else if (tbl.all1()) {
                //memo[id] = BoolExpr::createTrue();
                setMemoETS(id, BoolExpr::createTrue());
            }
            else {
                // gather children
                //std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> childF;
                //childF.resize(k);
                clearChildFETS();
                reserveChildFETS(k);
                for (uint32_t i = 0; i < k; ++i) {
                    size_t cid = node->children[i]->nodeID;
                    //childF[i] = memo[cid];
                    //childF[i] = getMemoETS(cid);
                    setChildFETS(i, getMemoETS(cid));
                }

                // find which inputs actually matter
                //std::vector<bool, tbb::tbb_allocator<bool>> relevant;
                clearRelevantETS();
                reserveRelevantETSwithFalse(k);
                //relevant.assign(k, false);
                for (uint32_t j = 0; j < k; ++j) {
                    for (uint64_t m = 0; m < rows; ++m) {
                        bool b0 = tbl.bits().bit(m);
                        bool b1 = tbl.bits().bit(m ^ (uint64_t{1} << j));
                        //if (b0 != b1) { relevant[j] = true; break; }
                        if (b0 != b1) { setRelevantETS(j, true); break; }
                    }
                }

                // collect the indices of relevant vars
                std::vector<uint32_t, tbb::tbb_allocator<uint32_t>> relIdx;
                for (uint32_t j = 0; j < k; ++j) {
                    //if (relevant[j]) relIdx.push_back(j);
                    if (getRelevantETS(j)) relIdx.push_back(j);
                }

                // if nothing matters, fall back to constant-false
                if (relIdx.empty()) {
                    //memo[id] = BoolExpr::createFalse();
                    setMemoETS(id, BoolExpr::createFalse());
                }
                else {
                    // build the DNF terms
                    //std::vector<std::shared_ptr<BoolExpr>, tbb::tbb_allocator<std::shared_ptr<BoolExpr>>> terms;
                    //terms.reserve(static_cast<size_t>(rows));
                    clearTermsETS();
                    reserveTermsETS(static_cast<size_t>(rows));

                    for (uint64_t m = 0; m < rows; ++m) {
                        if (!tbl.bits().bit(m)) continue;

                        std::shared_ptr<BoolExpr> term;
                        bool firstLit = true;
                        std::shared_ptr<BoolExpr> lit;
                        for (uint32_t j : relIdx) {
                            bool bit1 = ((m >> j) & 1) != 0;
                            lit = bit1
                                      ? /*childF[j]*/ getChildFETS(j)
                                      : /*BoolExpr::Not(childF[j])*/ BoolExpr::Not(getChildFETS(j));

                            if (firstLit) {
                                term = lit;
                                firstLit = false;
                            } else {
                                term = BoolExpr::And(term, lit);
                            }
                        }

                        // only push if we actually got a literal
                        if (term) {
                            //terms.push_back(std::move(term));
                            pushBackTermsETS(std::move(term));
                        }
                    }

                    // guard against an empty terms list
                    if (emptyTermsETS()) {
                        //memo[id] = BoolExpr::createFalse();
                        setMemoETS(id, BoolExpr::createFalse());
                    }
                    else {
                        // fold into OR
                        std::shared_ptr<BoolExpr> expr = getTErmsETS().first[0];
                        //for (size_t t = 1; t < terms.size(); ++t)
                        //    expr = BoolExpr::Or(expr, terms[t]);
                        for (size_t t = 1; t < sizeOfTermsETS(); ++t) {
                            expr = BoolExpr::Or(expr, getTErmsETS().first[t]);
                        }
                        //memo[id] = expr;
                        setMemoETS(id, expr);
                    }
                }
            }
        }
    }

    // 4) return root
    //return memo[root->nodeID];
    return getMemoETS(root->nodeID);
}
