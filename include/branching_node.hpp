#ifndef BRANCHING_NODE_HPP
#define BRANCHING_NODE_HPP

/**
 * @file branching_node.hpp
 * @brief This file implements the class @c BranchingNode and the class @c OneTree.
 *        It also provides an minimal encoding of OneTrees @c PredVec together with a distance
 *        function.
 *        Implementation is split into the files branching_node.cpp and one_tree.cpp.
 * **/

#include "graph.hpp"
#include <numeric>

namespace Core {
// imports from graph class
using Dist = TSP::Dist;
using NodeId = TSP::NodeId;
using NodeIds = std::vector<NodeId>;
NodeId constexpr invalid_node_id = std::numeric_limits<NodeId>::max();

/**
 * PredVec
 * @brief This class provides a minimal encoding of the edges of a OneTree.
 * For a graph with n nodes a PredVec is a vector of NodeIds of size n.
 * The first two entries encode the nodes which are adjacent to the node with id 0.
 * The following n - 2 entries encode a (min) spanning tree rooted at node 1,
 * i.e. for 1 < i < n - 1 the entry pred[i] is the predecessor of node i in that tree.
 * **/
using PredVec = NodeIds;
Dist dist(const PredVec &preds, const TSP::Graph &graph);

using size_type = std::size_t;
using Edge = std::pair<NodeId, NodeId>;
using Edges = std::vector<Edge>;
using Degree = int;
using Degrees = std::vector<Degree>;
// this matrix is used to store the forbidden edges
using BoolMatrix = std::vector<std::vector<bool>>;
// -1 is used as a constant as t_0 (i.e. the initial step length) is always > 0 (if set)
double constexpr T_ZERO_NOT_SET = double(-1.);
using Cost = double;
Cost constexpr infinite_cost = std::numeric_limits<Cost>::infinity();

class BranchingNode;
/**
 * @class OneTree
 * @brief A @c OneTree holds a @c PredVec encoding of the OneTree edges, i.e. a vector of node
 * degrees and the lower bound which it yields.
 * **/
class OneTree {
public:
  /** @brief Create a @c OneTree object. **/
  OneTree(const PredVec &preds, const Degrees &degrees, const double cost);
  /** @return True iff OneTree is a tour (-> checks degrees).**/
  bool is_tour() const;
  /** @return A node with degree > 2.**/
  NodeId hd_node() const;
  /** @return The neighbors of a node within the OneTree.**/
  NodeIds neighbors(const NodeId) const;
  /** @return The PredVec encoding of the OneTree. (-> see @c PredVec)**/
  const PredVec &preds() const;
  /** @return The lower bound times @c factor (rounded up, as Cost is int).**/
  Cost cost(double factor) const;

private:
  const PredVec _preds;
  const Degrees _degrees;
  const Cost _cost;
}; // class OneTree

/**
 * @class BranchingNode
 * @brief A @c BranchingNode holds: - a const reference to the input graph
 *                                  - a vector of lambda values
 *                                  - adjacency lists for the required edges
 *                                  - a @c BoolMatrix for the forbidden edges
 * @warning Calling compute_lb() changes the lambda vector.
 *
 * A adjacency list is used for the required edges as the length of each list is bounded by two and
 * the access time is therefore nearly constant.
 * The forbidden edges are stored in a BoolMatrix with also constant time access.
 * **/
class BranchingNode {
public:
  /** @brief Create an @c BranchingNode with no required nor forbidden edges and lambda == 0.**/
  BranchingNode(const TSP::Graph &graph);
  /** @brief Create an @c BranchingNode with given objects.**/
  BranchingNode(const TSP::Graph &graph, const std::vector<NodeIds> &is_required,
                const BoolMatrix &is_forbidden, const std::vector<double> &lambda);
  /** @brief Create an copy of *this with additional required and/or forbidden edges.**/
  BranchingNode create_child(const Edges &required_edges, const Edges &forbidden_edges) const;
  /** @return Two times MST length.**/
  Cost compute_ub() const;
  /** @return A @c OneTree with approx. maximizes the Held Karp lower bound.
   * @brief Use N iterations and t0 as initial step length.
   * **/
  OneTree compute_lb(const size_type N, const double t0);

  /** @return All required neighbors of @c node.**/
  NodeIds required_neighbors(const NodeId node) const;
  /** @return All forbidden neighbors of @c node.**/
  NodeIds forbidden_neighbors(const NodeId node) const;
  /** @return True iff corresponding edge is required. **/
  bool is_required(const NodeId, const NodeId) const;
  /** @return True iff corresponding edge is forbidden. **/
  bool is_forbidden(const NodeId, const NodeId) const;

  // public const access
  const TSP::Graph &graph() const;
  const std::vector<double> &lambda() const;
  const std::vector<NodeIds> &is_required() const;
  const BoolMatrix &is_forbidden() const;

private:
  /** @return A @c PredVec encoding of minimum OneTree wrt. the (by lambda) modified costs.**/
  PredVec mot() const;
  /** @return A @c PredVec encoding of minimum spanning wrt. the (by lambda) modified costs
   * excluding the first node. **/
  PredVec mst() const;

  /** @return By lambda modified distance. If (with_cases == true) returns inf for forbidden edges
   * and -inf for required edges. **/
  Cost edge_cost(const NodeId node1, const NodeId node2, bool with_cases) const;
  Cost edge_cost(const Edge &edge, bool with_cases) const;
  Cost edge_costs(const PredVec &preds, bool with_cases) const;

  const size_type _n;
  const TSP::Graph &_graph;
  std::vector<double> _lambda;
  const std::vector<NodeIds> _is_required;
  const BoolMatrix _is_forbidden;
}; // class BranchingNode

// BEGIN: Inlince Section
inline Cost OneTree::cost(double factor) const { return std::ceil(_cost * factor); };
inline const PredVec &OneTree::preds() const { return _preds; };

inline const TSP::Graph &BranchingNode::graph() const { return _graph; };
inline const std::vector<double> &BranchingNode::lambda() const { return _lambda; };
inline const std::vector<NodeIds> &BranchingNode::is_required() const { return _is_required; };
inline const BoolMatrix &BranchingNode::is_forbidden() const { return _is_forbidden; };
// assumes that there are no more than 2 required edges
inline bool BranchingNode::is_required(const NodeId node1, const NodeId node2) const {
  if (_is_required[node1].empty()) {
    return false;
  }
  return (_is_required[node1].front() == node2 or _is_required[node1].back() == node2);
};
inline bool BranchingNode::is_forbidden(const NodeId node1, const NodeId node2) const {
  return _is_forbidden[node1][node2];
};
// END:: Inline Section
} // namespace Core

#endif
