#ifndef SNLTRUTHTABLETREE_H
#define SNLTRUTHTABLETREE_H

#include "SNLTruthTable.h"
#include <vector>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <atomic>
#include <tbb/tbb_allocator.h>
#include <unordered_set>
#include "DNL.h"
#include "SNLDesignModeling.h"

namespace KEPLER_FORMAL {

// Compact id-based truth-table tree (no pointer mirrors)
class SNLTruthTableTree {
public:
  struct Node {
    uint32_t inputIndex = std::numeric_limits<uint32_t>::max();
    // debug counters
    uint32_t nodeID = std::numeric_limits<uint32_t>::max();   // incremental debug id for diagnostics (unchanged)
    //uint32_t id_   = std::numeric_limits<uint32_t>::max(); // actual id (index + kIdOffset)
    uint32_t parentId = std::numeric_limits<uint32_t>::max();

    // owner pointer for lookup convenience in Node methods
    SNLTruthTableTree* tree = nullptr;

    naja::DNL::DNLID termid = naja::DNL::DNLID_MAX;

    // canonical id-based topology: children ids and parent id
    
    std::vector<uint32_t, tbb::tbb_allocator<uint32_t>> childrenIds;
    enum class Type : uint8_t { Input = 0, Table = 1, P = 2 } type;

    explicit Node(uint32_t idx, SNLTruthTableTree* t);
    Node(SNLTruthTableTree* t,
         naja::DNL::DNLID instid,
         naja::DNL::DNLID term,
         Type type_ = Type::Table);

    Node(const Node& other) = delete;
    ~Node();

    bool eval(const std::vector<bool>& extInputs) const;
    void addChildId(uint32_t childId);
    const SNLTruthTable& getTruthTable() const;
  };

  static constexpr uint32_t kReservedId0 = 0u;
  static constexpr uint32_t kReservedId1 = 1u;
  static constexpr uint32_t kIdOffset = 2u; // id = index + kIdOffset
  static constexpr uint32_t kInvalidId = std::numeric_limits<uint32_t>::max();

  SNLTruthTableTree();
  SNLTruthTableTree(naja::DNL::DNLID instid, naja::DNL::DNLID termid, Node::Type type = Node::Type::Table);

  size_t size() const;
  bool eval(const std::vector<bool>& extInputs) const;

  void concat(size_t borderIndex,
              naja::DNL::DNLID instid,
              naja::DNL::DNLID termid);

  void concatFull(const std::vector<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>,
            tbb::tbb_allocator<std::pair<naja::DNL::DNLID, naja::DNL::DNLID>>>& tables);

  uint32_t getRootId() const { return rootId_; }
  const std::shared_ptr<Node>& getRootShared() const { return nodeFromId(rootId_); }
  const std::shared_ptr<Node>& getRoot() const { return getRootShared(); }

  const std::shared_ptr<Node>& nodeFromId(uint32_t id) const;
  bool isInitialized() const;
  void print() const;
  void simplify();
  void destroy();

  size_t getNumNodes() const { return nodes_.size(); }

  // allocateNode guarantees id assignment before publishing node in nodes_
  uint32_t allocateNode(const std::shared_ptr<Node>& np);

  // finalize repairs and validates the tree; must be called once after build and before traversal
  // It will remap children/parent ids to canonical ids and throw on unresolved references
  void finalize();

private:
  struct BorderLeaf {
    uint32_t parentId;
    size_t childPos;
    size_t extIndex;
  };

  const Node& concatBody(size_t borderIndex,
                         naja::DNL::DNLID instid,
                         naja::DNL::DNLID termid);

  void updateBorderLeaves();

  std::vector<std::shared_ptr<Node>, tbb::tbb_allocator<std::shared_ptr<Node>>> nodes_;
  uint32_t rootId_ = kInvalidId;
  size_t numExternalInputs_ = 0;
  std::vector<BorderLeaf, tbb::tbb_allocator<BorderLeaf>> borderLeaves_;
  size_t lastID_ = 2;       // debug counter for nodeID assignment
  static const SNLTruthTable PtableHolder_;
};

} // namespace KEPLER_FORMAL

#endif // SNLTRUTHTABLETREE_H
