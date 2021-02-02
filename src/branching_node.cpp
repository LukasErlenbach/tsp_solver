#include "branching_node.hpp"
/**
 * @file branching_node.hpp
 * @brief Implementation of class @c BranchingNode and global functions
 * **/
namespace Core {

/////////////////////////////////////////////
//! global functions
/////////////////////////////////////////////

Dist dist(const PredVec &preds, const TSP::Graph &graph) {
  Dist sum = graph.dist(0, preds[0]);
  sum += graph.dist(0, preds[1]);
  for (NodeId i = 2; i < preds.size(); ++i) {
    sum += graph.dist(i, preds[i]);
  }
  return sum;
}

/////////////////////////////////////////////
//! \c BranchingNode definitions
/////////////////////////////////////////////

// this constuctor is used for the root node
BranchingNode::BranchingNode(const TSP::Graph &graph)
    : _n(graph.num_nodes()), _graph(graph), _lambda(std::vector<double>(graph.num_nodes(), 0.)),
      _is_required(std::vector<NodeIds>(graph.num_nodes(), NodeIds({}))),
      _is_forbidden(BoolMatrix(graph.num_nodes(), std::vector<bool>(graph.num_nodes(), false))) {}

// this constuctor is used for creating child nodes in the branching tree
BranchingNode::BranchingNode(const TSP::Graph &graph, const std::vector<NodeIds> &is_required,
                             const BoolMatrix &is_forbidden, const std::vector<double> &lambda)
    : _n(graph.num_nodes()), _graph(graph), _lambda(lambda), _is_required(is_required),
      _is_forbidden(is_forbidden) {}

// creates a new BranchingNode which is a copy of *this with additional
// required and/or forbidden edges
BranchingNode BranchingNode::create_child(const Edges &required_edges,
                                          const Edges &forbidden_edges) const {
  // init encodings for the new node
  std::vector<NodeIds> is_required = _is_required;
  BoolMatrix is_forbidden = _is_forbidden;
  // make those edges required
  for (Edge edge : required_edges) {
    is_required[edge.first].push_back(edge.second);
    is_required[edge.second].push_back(edge.first);
    const NodeIds adj_to_neighbor = required_neighbors(edge.second);
    // neighbor is already incident to another required edge?
    // -> forbid all other edges
    if (adj_to_neighbor.size() == 1) {
      for (NodeId to_forbid = 0; to_forbid < _n; ++to_forbid) {
        if (to_forbid != adj_to_neighbor[0] and to_forbid != edge.first) {
          is_forbidden[edge.second][to_forbid] = true;
          is_forbidden[to_forbid][edge.second] = true;
        }
      }
    }
  }
  // require 2 edges at a node? -> forbid all other edges
  if (required_edges.size() == 2) {
    const NodeId i = required_edges[0].first;
    const NodeId r1 = required_edges[0].second;
    const NodeId r2 = required_edges[1].second;
    for (NodeId to_forbid = 0; to_forbid < _n; ++to_forbid) {
      if (to_forbid != r1 and to_forbid != r2) {
        is_forbidden[i][to_forbid] = true;
        is_forbidden[to_forbid][i] = true;
      }
    }
  }
  // make those edges forbidden
  for (Edge edge : forbidden_edges) {
    is_forbidden[edge.first][edge.second] = true;
    is_forbidden[edge.second][edge.first] = true;
  }
  return BranchingNode(_graph, is_required, is_forbidden, _lambda);
}

// returns 2 times MST length
Cost BranchingNode::compute_ub() const {
  PredVec tree = mst();
  // find cheapest edge from node 0 since it is excluded in mst computation
  Cost min_edge_cost = edge_cost(0, 1, false);
  for (NodeId node_id = 2; node_id < _n; ++node_id) {
    const Cost cost = edge_cost(0, node_id, false);
    if (cost < min_edge_cost) {
      min_edge_cost = cost;
    }
  }
  return 2 * (edge_costs(tree, false) + min_edge_cost);
}

// helper function: computed degrees from a PredVec encoding
Degrees to_degrees(const PredVec &preds) {
  Degrees d(preds.size(), 0);
  for (NodeId id = 2; id < preds.size(); ++id) {
    ++d[id];
    ++d[preds[id]];
  }
  d[0] += 2;
  ++d[preds[0]];
  ++d[preds[1]];
  return d;
}

// compute Held Karp lower bound in N steps with initial step length t0
OneTree BranchingNode::compute_lb(const size_type N, const double t0) {
  double t = t0;
  PredVec first_tree = mot();
  // Graph without node 0 is not connected -> discard
  if (first_tree.size() == 0) {
    return OneTree({}, {}, infinite_cost);
  }
  // this case occurs in the first iteration
  // T_ZERO_NOT_SET is defined in branching_node.hpp
  // set (t0), delta, delta_delta according paper VJ[1980]
  if (t == T_ZERO_NOT_SET) {
    t = edge_costs(first_tree, false) / (2 * _n);
  }
  double delta = 3 * t / (2 * N);
  const double delta_delta = t / (N * N - N);

  Cost max_cost = 0.;
  PredVec max_tree = {};
  std::vector<double> best_lambda{};

  // include degrees from previous iteration in the update to damp oscillations
  Degrees last_d = to_degrees(first_tree);

  for (size_type i = 1; i < N; ++i) {
    // compute min OneTree for fixed lambda
    PredVec tree = mot();
    // check for maximum
    Cost tree_cost = edge_costs(tree, false);

    if (tree_cost > max_cost) {
      max_cost = tree_cost;
      max_tree = tree;
      best_lambda = _lambda;
    }
    Degrees d = to_degrees(tree);
    // update lambda
    for (size_type idx = 0; idx < d.size(); ++idx) {
      _lambda[idx] += t * (0.6 * (d[idx] - 2) + 0.4 * (last_d[idx] - 2));
    }
    // update last degrees, step length and delta
    last_d = std::move(d);
    t -= delta;
    delta -= delta_delta;
  }
  _lambda = best_lambda;
  return OneTree(max_tree, to_degrees(max_tree), std::ceil(max_cost));
}

// return adjacency list
NodeIds BranchingNode::required_neighbors(const NodeId node_id) const {
  return _is_required[node_id];
}

// first build adjacency list, then return it
NodeIds BranchingNode::forbidden_neighbors(const NodeId node_id) const {
  std::vector<NodeId> neighbors = {};
  for (size_type idx = 0; idx < _n; ++idx) {
    if (_is_forbidden[node_id][idx]) {
      neighbors.push_back(idx);
    }
  }
  return neighbors;
}

// return minimal length OneTree
PredVec BranchingNode::mot() const {
  // start with a MST (which excludes the node 0)
  // the MST is encoded in entries 2,...,n-1
  // -> use entry 0 and 1 below for the edges from node 0
  PredVec mot = mst();
  if (mot.size() == 0) {
    // -> graph is not connected
    return {};
  }
  // find cheapest 2 edges and connect node 0 to MST
  Cost min_cost = edge_cost(0, 1, true);
  NodeId min_id = 1;
  Cost min_cost_2 = edge_cost(0, 2, true);
  NodeId min_id_2 = 2;
  if (min_cost > min_cost_2) {
    std::swap(min_cost, min_cost_2);
    std::swap(min_id, min_id_2);
  }
  for (NodeId id = 3; id < _n; ++id) {
    const Cost cost = edge_cost(0, id, true);
    if (cost < min_cost) {
      min_cost_2 = min_cost;
      min_id_2 = min_id;
      min_cost = cost;
      min_id = id;
      continue;
    }
    if (cost < min_cost_2) {
      min_cost_2 = cost;
      min_id_2 = id;
    }
  }
  // add edges accordingly
  mot[0] = min_id;
  mot[1] = min_id_2;
  return mot;
}

// excludes node 0
// uses prims algorithm to find a MST
PredVec BranchingNode::mst() const {
  // for every node the cheapest connecting edge to the tree is stored
  // -> init in O(n), update n-2 times in O(n)
  // --> O(n^2) runtime
  using CostNodePair = std::pair<Cost, NodeId>;
  const CostNodePair invalid_pair = std::make_pair(infinite_cost, invalid_node_id);
  // for cheapest connection
  std::vector<CostNodePair> cpst_conn(_n, invalid_pair);

  PredVec mst(_n, 0);
  // edge to be added next
  Edge next_edge;
  Cost next_cost = infinite_cost;

  // gets called on the node which we just added to the tree
  // -> scan all edges to nodes that are not covered yet and update their connection costs
  // -> while doing so find the cheapest connection to an uncovered node for the next iteration
  auto update_min = [&](const NodeId node) {
    cpst_conn[node] = invalid_pair;
    // exclude node 0
    // exclude node 1 as well because this is the root of the mst and therefore always covered
    for (NodeId other_node = 2; other_node < _n; ++other_node) {
      // mst[other_node] != 0 -> predecessor set -> already covered by the tree
      if (other_node == node or mst[other_node] != 0) {
        continue;
      }
      const Cost cost = edge_cost(node, other_node, true);
      // Edge (node, other_node) is a cheaper connection to the tree?
      // -> update cpst_conn[other_node]
      if (cost < cpst_conn[other_node].first) {
        cpst_conn[other_node] = std::make_pair(cost, node);
      }
      // Since all cheapest edges are scanned for all uncovered node this finds the next edge
      // to add to the tree
      if (cpst_conn[other_node].first < next_cost) {
        next_cost = cpst_conn[other_node].first;
        next_edge = std::make_pair(cpst_conn[other_node].second, other_node);
      }
    }
  };

  // Prims algorithm
  // starting from node 1
  // WARNING: the starting node is arbitrary in general, BUT since the PredVec format depends on
  // the root of the MST there are parts in the code that depend on this to start at node 1
  const NodeId cur_node = 1;
  update_min(cur_node);
  for (size_type iter = 0; iter + 2 < _n; ++iter) {
    // graph is not connected (wrt forbidden and required edges)
    // this case can occur for example if we require a circle smaller than n nodes
    if (next_cost == infinite_cost) {
      return {};
    }
    mst[next_edge.second] = next_edge.first;
    next_cost = infinite_cost;
    update_min(next_edge.second);
  }
  return mst;
}

// if with_cases is true return inf for forbidden and -inf for required edges
Cost BranchingNode::edge_cost(const NodeId node1, const NodeId node2, bool with_cases) const {
  if (with_cases and _is_forbidden[node1][node2]) {
    return infinite_cost;
  }
  if (with_cases and is_required(node1, node2)) {
    return -infinite_cost;
  }
  return Cost(_graph.dist(node1, node2)) + _lambda[node1] + _lambda[node2];
}

Cost BranchingNode::edge_cost(const Edge &edge, bool with_cases) const {
  return edge_cost(edge.first, edge.second, with_cases);
}

// respects the PredVec encoding
Cost BranchingNode::edge_costs(const PredVec &preds, bool with_cases) const {
  Cost sum = 0.;
  sum += edge_cost(0, preds[0], with_cases);
  sum += edge_cost(0, preds[1], with_cases);
  for (size_type idx = 2; idx < preds.size(); ++idx) {
    sum += edge_cost(idx, preds[idx], with_cases);
  }
  return sum;
}
} // namespace Core
