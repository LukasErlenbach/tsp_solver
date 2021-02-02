#ifndef GRAPH_HPP
#define GRAPH_HPP

/**
 * @file graph.hpp
 *
 *  @brief This file implements a simple class @c Graph to model complete graphs induced by points
 * in the plane. It supports floating coordinates and rounded Euclidean length.
 * **/
#include <cmath>
#include <cstddef> // std::size_t
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace TSP {

using size_type = std::size_t;
using NodeId = size_type;
using TsplibId = size_type;
using NodeIds = std::vector<NodeId>;
/** This is used for rounded Euclidean length. **/
using Dist = int;

/** Useful constant different from the id of any actual node: **/
NodeId constexpr invalid_node_id = std::numeric_limits<NodeId>::max();
TsplibId constexpr invalid_tpslib_id = std::numeric_limits<TsplibId>::max();

NodeId from_tsplib_id(TsplibId const tsplib_id); //!< Subtracts 1 (throws if @c tsplib_id is 0)
TsplibId to_tsplib_id(NodeId const node_id);     //!< Adds 1 (throws if overflow would occur)

/**
 * @class Graph
 *
 * @brief A @c Graph holds 2 arrays with float x (resp. y) coordinates as well as an complete
 * distance Matrix with rounded Euclidean length.
 *
 * @note Provides a function to read in from a file in TSPLIB format.
 *
 * **/
class Graph {
public:
  /** @brief Create an empty graph on num_nodes nodes. **/
  Graph(const NodeId num_nodes);

  /** @brief Create an graph from TSPLIB file. Throws if format is incorrect.**/
  void init_from_file(const std::string &filename);

  /** @return The number of nodes in the graph.
   *  @warning Assumes _x.size() == _y.size() (which should always hold).**/
  size_type num_nodes() const;

  /** @return Rounded Euclidean length from the dist matrix _dist.
   * @warning Does not perform bound checks.
   * **/
  Dist dist(NodeId node1, NodeId node2) const;

private:
  /** @brief Set coordinates of node @c id to x and y in _x and _y.
   * @warning Does not perform bound checks.
   * **/
  void set_coords(const NodeId id, const double x, const double y);
  /** @brief Compute complete distance matrix _dist.
   * @warning Does not perform bound checks.
   **/
  void init_dist();
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<std::vector<Dist>> _dist;
}; // class Graph

// BEGIN: Inline section
inline size_type Graph::num_nodes() const { return _x.size(); }
// END: Inline section

} // namespace TSP

#endif /* GRAPH_HPP */
