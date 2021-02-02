#ifndef SOLVER_HPP
#define SOLVER_HPP

/**
 * @file solver.hpp
 *
 * @brief This file implements a TSP solver class @c Solver.
 *        The Solver uses the class @c BranchingNode which is specified in the file
 *        branching_node.hpp.
 *        The Solver performs the outer Branch and Bound algorithm while using the lower and upper
 *        bounds provided by the @c BranchingNode class.
 * **/

#include "branching_node.hpp"
#include <iostream>

namespace Solver {
using size_type = std::size_t;
/** used for lower and upper bounds
 -> need float type for HeldKarp bound **/
using Cost = double;
/** used only for distances, i.e. sums of edge lengths **/
using Dist = Core::Dist;

// use definition of PredVec in file branching_node.hpp
using PredVec = Core::PredVec;
// encode Tours as PredVec
using Tour = PredVec;
using NodeId = Core::NodeId;
using NodeIds = Core::NodeIds;
using BranchingNode = Core::BranchingNode;
using OneTree = Core::OneTree;

using Edge = Core::Edge;
using Edges = Core::Edges;

// definition in file branching_node.hpp
// -> this is used in the first iteration of the branch and bound algorithm
double constexpr T_ZERO_NOT_SET = Core::T_ZERO_NOT_SET;

/**
 * @class Solver
 * @brief The Solver stores a const reference to the input @c Graph.
 * **/
class Solver {
public:
  /** @brief Create Solver object with access to a @c Graph.**/
  Solver(const TSP::Graph &graph);
  /**
   * @return A shortest TSP Tour in the @c Graph _graph.
   * @brief Performs the whole Branch and Bound business using the best-bound heuristic.
   *        Argument @c eps is used to scale down lower bounds by (1. - eps).
   *        For more details on the bounds see @c BranchingNode.
   * **/
  Tour run(const double eps) const;
  /**
   * @brief Write @c tour in TSPLIB format to the file with name "filename".
   * (throws if file can not be created)
   * **/
  void write_tour_to_file(const Tour &tour, const std::string &filename) const;

private:
  /**
   * @return Two edges from @c node to @adj_nodes, which are not forbidden in @bnode.
   * @brief Selects two arbitrary not yet required edges. (see comment in branching_node.cpp)
   * **/
  Edges find_branching_edges(const NodeId node, const BranchingNode &bnode,
                             const NodeIds &adj_nodes) const;
  /**
   * @return True if no node has less than two non forbidden edges and false otherwise.
   * @brief We use this to check whether a @c BranchingNode @c bnode has to many forbidden edges
   * at a node, if so we can discard the branch represented by the @c bnode.
   * **/
  bool check_valid(const BranchingNode &bnode) const;
  const TSP::Graph &_graph;
}; // class Solver
} // namespace Solver
#endif
