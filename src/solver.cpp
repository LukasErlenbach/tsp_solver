#include "solver.hpp"
#include <fstream>
#include <queue>

/**
 * @file solver.cpp
 * @brief Implementation of class @c Solver.
 * **/

namespace Solver {

/////////////////////////////////////////////
//! \c Solver definitions
/////////////////////////////////////////////

Solver::Solver(const TSP::Graph &graph) : _graph(graph) {}

// check all incident edges and return as soon as two not required edges are found
Edges Solver::find_branching_edges(const NodeId node, const BranchingNode &bnode,
                                   const NodeIds &adj_nodes) const {
  Edges branching_edges = {};
  for (NodeId adj_node : adj_nodes) {
    if (not bnode.is_required(node, adj_node)) {
      branching_edges.push_back(std::make_pair(node, adj_node));
      if (branching_edges.size() == 2) {
        break;
      }
    }
  }
  return branching_edges;
}

// return false iff there exists a node in the graph with only 0 or 1 allowed edges
bool Solver::check_valid(const BranchingNode &bnode) const {
  const size_type n = _graph.num_nodes();
  for (NodeId node = 0; node < n; ++node) {
    if (bnode.forbidden_neighbors(node).size() > n - 2) {
      return false;
    }
  }
  return true;
}

// write tour in TSPLIB format to file named "filename"
void Solver::write_tour_to_file(const Tour &tour, const std::string &filename) const {
  const std::size_t n = tour.size();
  std::ofstream file;
  file.open(filename);
  if (file.is_open()) {
    file << "TYPE : TOUR\n";
    file << "DIMENSION : " << n << "\n";
    file << "TOUR_SECTION\n";

    // the first Ids in tour are the nodes adjacent to node 0 in the tour
    // -> write them in some orientation
    file << TSP::to_tsplib_id(tour[0]) << "\n";
    file << "1\n";
    file << TSP::to_tsplib_id(tour[1]) << "\n";
    // the nodes tour[3],...,tour.back() encode the predecessors of the index nodes in the
    // mst rooted at vertex 1 (which is just a path, since tour is a tour)
    // -> to get the correct orientation first write the path from tour[1] to 1
    NodeId runner = tour[1];
    while (tour[runner] != NodeId(1)) {
      file << TSP::to_tsplib_id(tour[runner]) << "\n";
      runner = tour[runner];
    }
    runner = tour[0];
    // if node 1 is adjacent to node 0 we are already done
    if (runner == NodeId(1)) {
      file << "-1\nEOF\n";
      file.close();
    }
    // otherwise first collect to path from tour[0] to 1 and then write it backwards
    file << "2\n";
    NodeIds path_to_root({});
    while (tour[runner] != NodeId(1)) {
      path_to_root.push_back(TSP::to_tsplib_id(tour[runner]));
      runner = tour[runner];
    }
    for (int backwards_idx = path_to_root.size() - 1; backwards_idx >= 0; --backwards_idx) {
      file << path_to_root[backwards_idx] << "\n";
    }
    file << "-1\nEOF\n";
    file.close();
  } else {
    throw std::runtime_error("Unable to open solution file: " + filename);
  }
}

// run the Branch and Bound algorithm to determine a shortest TSP tour on _graph
Tour Solver::run(const double eps) const {
  // first handle corner cases
  const size_type n = _graph.num_nodes();
  if (n < 3) {
    throw std::runtime_error("Can not compute a TSP tour on less than 3 points.");
  }
  if (n == 3) {
    std::cout << _graph.dist(0, 1) + _graph.dist(1, 2) + _graph.dist(0, 2) << std::endl;
    return {1, 2, 1};
  }

  //
  // Setup: define all parameters needed during the Branch and Bound algorithm
  //        - N_root, N_child: number of iterations for the HK bound
  //        - t0: initial step length for the HK bound
  //          -> from paper by Volgenant & Jonker [1980]
  //        - U: the current upper bound
  //        - root: first BranchingNode (no forbidden or required edges)
  //        - root_tree: opt OneTree for on root (i.e. on K_n, with dists from _graph)
  //        - Q: priority which stores pairs of ids and lower bounds which correspond
  //             to the BranchingNode and OneTree in the above vectors with that id
  //          -> this is used to perform the best bound heuristic for B&B
  //        - best_tour: best tours found so far
  //
  const size_type N_root = std::ceil(n * n / 50) + n + 15;
  const size_type N_child = std::ceil(n / 4) + 5;
  BranchingNode root(_graph);
  // 2 times MST (+1 to make sure we run at least 1 iteration even if this is already optimal)
  Cost U = root.compute_ub() + 1.;
  OneTree root_tree = root.compute_lb(N_root, T_ZERO_NOT_SET);
  const double t0 = 1. / (2. * n) *
                    std::accumulate(root.lambda().begin(), root.lambda().end(), 0.,
                                    [](double sum, double l) { return sum + std::fabs(l); });
  if (t0 == 0.0) {
    throw std::runtime_error("Step size (t0) of zero is not valid.");
  }
  std::vector<BranchingNode> branching_nodes({root});
  std::vector<OneTree> trees({root_tree});

  using CostIdPair = std::pair<Cost, size_type>;
  auto CostIdCmp = [](const CostIdPair &left, const CostIdPair &right) {
    return left.first > right.first;
  };
  std::priority_queue<CostIdPair, std::vector<CostIdPair>, decltype(CostIdCmp)> Q(CostIdCmp);

  Q.push(std::make_pair(root_tree.cost(1. - eps), 0));
  Tour best_tour;

  std::size_t num_proc{0};
  double lb{0};

  //
  // Branch and Bound implementation
  //
  while (not Q.empty()) {
    const size_type cur_idx = Q.top().second;
    Q.pop();
    ++num_proc;
    // logging output to cerr
    auto status = [&]() {
      std::cerr << "LB " << lb << ", UB " << U << ", #processed nodes " << num_proc
                << ", #nodes in queue " << Q.size() << std::endl;
    };

    auto handle_tour = [&](OneTree const &tree) {
      // the OneTree is a optimal tour (wrt. R and F) -> check if we can update U
      const Dist tour_length = Core::dist(tree.preds(), _graph);
      if (tour_length < U) {
        U = tour_length;
        best_tour = tree.preds();
        status();
      }
    };

    // lower bound > U? -> discard node
    auto const &cur_tree = trees[cur_idx];
    if (trees[cur_idx].cost(1. - eps) < U) {
      auto const cost = trees[cur_idx].cost(1. - eps);
      if (cost > lb) {
        lb = cost;
        status();
      }

      if (cur_tree.is_tour()) {
        handle_tour(cur_tree);
        continue;
      }

      // i is the node that we want to branch on
      // -> first collect adjacent required and tree edges
      // -> then find branching edges
      const NodeId i = trees[cur_idx].hd_node();
      const NodeIds adj_req = branching_nodes[cur_idx].required_neighbors(i);
      const NodeIds adj_tree = trees[cur_idx].neighbors(i);
      Edges e = find_branching_edges(i, branching_nodes[cur_idx], adj_tree);

      // Creates a child node of the given BranchingNode by extending the required edges
      // by req, and the forbidden edges by fbd
      // -> also check if the child is valid and not to costly
      auto update_Q = [&](const BranchingNode &bnode, const Edges &req, const Edges &fbd) {
        BranchingNode new_bnode = bnode.create_child(req, fbd);
        const OneTree new_tree = new_bnode.compute_lb(N_child, t0);
        if (new_tree.cost(1. - eps) < U and check_valid(new_bnode)) {
          if (new_tree.is_tour()) {
            handle_tour(new_tree);
          } else {
            branching_nodes.push_back(new_bnode);
            trees.push_back(new_tree);
            Q.push(std::make_pair(new_tree.cost(1. - eps), trees.size() - 1));
          }
        }
      };
      // actual branching
      // no incident required edges? -> split into 3 BranchingNodes
      //                             -> else only 2
      if (adj_req.size() == 0) {
        update_Q(branching_nodes[cur_idx], {}, {e[1]});
        update_Q(branching_nodes[cur_idx], {e[1]}, {e[0]});
        update_Q(branching_nodes[cur_idx], {e[1], e[0]}, {});
      } else {
        // the edge r is already required
        // but we give it to the BranchingNode again so that it forbids
        // all but the now 2 required edges
        const Edge r = std::make_pair(i, adj_req[0]);
        update_Q(branching_nodes[cur_idx], {}, {e[1]});
        update_Q(branching_nodes[cur_idx], {e[1], r}, {e[0]});
      }
    }
  }
  // Branching finished -> output length
  std::cout << "minimal TSP tour length: " << Core::dist(best_tour, _graph) << std::endl;
  return best_tour;
}
} // namespace Solver
