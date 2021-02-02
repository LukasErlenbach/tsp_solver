#include "branching_node.hpp"

/**
 * @file one_tree.hpp
 * @brief Implementation of class @c OneTree.
 * **/

namespace Core {

/////////////////////////////////////////////
//! \c OneTree definitions
/////////////////////////////////////////////

OneTree::OneTree(const PredVec &preds, const Degrees &degrees, const double cost)
    : _preds(preds), _degrees(degrees), _cost(cost) {}

// this test is sufficient since PredVecs are assumed to be connected (since
// they encode OneTrees)
bool OneTree::is_tour() const {
  for (Degree d : _degrees) {
    if (d != 2) {
      return false;
    }
  }
  return true;
}

// greedily finds a node with degree > 2
NodeId OneTree::hd_node() const {
  for (NodeId node_id = 0; node_id < _degrees.size(); ++node_id) {
    if (_degrees[node_id] > 2) {
      return node_id;
    }
  }
  // this case should never occur, we only call this after checking with
  // is_tour()
  throw std::runtime_error("Error in OneTree::hd_node(): No node with degree > 2.");
  return invalid_node_id;
}

NodeIds OneTree::neighbors(const NodeId node_id) const {
  NodeIds neighbors = {};
  neighbors.reserve(_degrees[node_id]);
  if (node_id == 0) {
    // -> first 2 entrys are neighbors
    neighbors.push_back(_preds[0]);
    neighbors.push_back(_preds[1]);
  }
  if (node_id == _preds[0] or node_id == _preds[1]) {
    // -> adjacent to node 0
    neighbors.push_back(0);
  }
  for (NodeId other_id = 2; other_id < _preds.size(); ++other_id) {
    if (_preds[other_id] == node_id) {
      // -> other_id is successor in the OneTree
      neighbors.push_back(other_id);
      continue;
    }
    if (other_id == node_id) {
      // -> other_id is predecessor in the OneTree
      neighbors.push_back(_preds[node_id]);
    }
  }

  return neighbors;
}
} // namespace Core
