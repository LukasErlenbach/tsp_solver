/**
 *  @file graph.cpp
 *
 *  @brief Implementation of class @c Graph.
 * **/

#include "graph.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

namespace TSP {

/////////////////////////////////////////////
//! global functions
/////////////////////////////////////////////

NodeId from_tsplib_id(TsplibId const tsplib_id) {
  if (tsplib_id == 0) {
    throw std::runtime_error("Invalid (0) TSPLIB id.");
  }

  return static_cast<NodeId>(tsplib_id - 1);
}

TsplibId to_tsplib_id(NodeId const node_id) {
  if (node_id == std::numeric_limits<NodeId>::max()) {
    throw std::runtime_error("Invalid (inf) node id.");
  }

  return static_cast<TsplibId>(node_id + 1);
}

/////////////////////////////////////////////
//! \c Graph definitions
/////////////////////////////////////////////

Graph::Graph(NodeId const num_nodes)
    : _x(std::vector<double>(num_nodes, 0)), _y(std::vector<double>(num_nodes, 0)),
      _dist(std::vector<std::vector<Dist>>(num_nodes, std::vector<Dist>(num_nodes, 0.))) {}

// removes colons and replace consecutive spaces with a single one
// -> used for parsing the TSPLIB file below
void remove_colon_and_space(std::string &str) {
  str.erase(std::remove(str.begin(), str.end(), ':'), str.end());
  std::string::iterator new_end = std::unique(
      str.begin(), str.end(), [=](char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); });
  str.erase(new_end, str.end());
}

// reads in a graph from a file named "filename" which is expected to be in TSPLIB format
void Graph::init_from_file(const std::string &filename) {
  std::string line;
  std::ifstream file;
  file.open(filename);
  if (file.is_open()) {
    // find problem size and resize member vectors accordingly
    while (getline(file, line)) {
      remove_colon_and_space(line);
      std::istringstream line_buffer(line);
      std::string code_word;
      if (line_buffer >> code_word) {
        if (code_word != "DIMENSION") {
          continue;
        }
        size_type n;
        if (line_buffer >> n) {
          _x.resize(n, 0);
          _y.resize(n, 0);
          _dist.resize(n, std::vector<Dist>(n, 0.));
        } else {
          throw std::runtime_error("Error reading input file. (-> expect TSPLIB format!) ");
        }
        break;
      } else {
        throw std::runtime_error("Error reading input file. (-> expect TSPLIB format!)");
      }
    }
    // find coordinates and read them in
    while (getline(file, line)) {
      remove_colon_and_space(line);
      std::istringstream line_buffer(line);
      std::string code_word;
      if (line_buffer >> code_word) {
        if (code_word != "NODE_COORD_SECTION") {
          continue;
        }
        // read in coords
        for (size_type num_lines = 0; num_lines < num_nodes(); ++num_lines) {
          if (not getline(file, line)) {
            throw std::runtime_error("Error reading input file. (-> expect TSPLIB format!)");
            break;
          }
          NodeId i;
          double x, y;
          remove_colon_and_space(line);
          std::istringstream line_buffer(line);
          if (line_buffer >> i >> x >> y) {
            set_coords(from_tsplib_id(i), x, y);
          } else {
            throw std::runtime_error("Error reading input file. (-> expect TSPLIB format!)");
          }
        }
      }
    }
    // init distance matrix
    init_dist();
  } else {
    throw std::runtime_error("Error reading input file. Cannot open file. ");
  }
}

void Graph::set_coords(const NodeId node_id, const double x, const double y) {
  if (node_id >= _x.size()) {
    throw std::runtime_error("Tried to set coords of non existing node. ");
  } else {
    _x[node_id] = x;
    _y[node_id] = y;
  }
}

Dist Graph::dist(NodeId id1, NodeId id2) const { return _dist[id1][id2]; }

void Graph::init_dist() {
  for (NodeId id1 = 0; id1 < num_nodes(); ++id1) {
    for (NodeId id2 = id1; id2 < num_nodes(); ++id2) {
      // rounded Euclidean length
      // -> checked correctness with the official TSPLIB documentation (which provides length of
      // some canonical tours)
      const Dist dist = std::lround(std::sqrt((_x[id1] - _x[id2]) * (_x[id1] - _x[id2]) +
                                              (_y[id1] - _y[id2]) * (_y[id1] - _y[id2])));
      _dist[id1][id2] = dist;
      _dist[id2][id1] = dist;
    }
  }
}
} // namespace TSP
