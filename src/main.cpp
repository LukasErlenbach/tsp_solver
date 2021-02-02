/**
  @file main.cpp

  @brief Handle program arguments and run solver.
 **/

#include "solver.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using Tour = Solver::Tour;

int main(int argc, char **argv) {
  if (argc != 3 and argc != 5) {
    std::cout << "\nThis program is an exact TSP solver based on Volgenant and Jonker [1980].\n";

    std::cout << "Expected Input:\n\n"
              << "--instance file.tsp [--solution opt_tour.tsp]\n\n"
              << "-> file.tsp is expected to be in TSPLIP format\n"
              << "-> if provided, an optimal tour will be written to opt_tour.tsp\n"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::string argv1 = argv[1];
  if (argv1 != "--instance") {
    std::cout << "Unexpected keyword \"" << argv1
              << "\" at 1. argument. (Run w/o arguments for help.)" << std::endl;
    return EXIT_FAILURE;
  }
  TSP::Graph graph(0);
  // try to init graph from file
  try {
    graph.init_from_file(argv[2]);
  } catch (const std::runtime_error &error) {
    std::cout << error.what() << std::endl;
    return EXIT_FAILURE;
  }
  // try to run solver on input graph
  // throws if graph.num_nodes() < 3
  try {
    Solver::Solver solver(graph);
    // scale lower bounds with double(1. - 0.001)
    Tour opt = solver.run(0.001);
    if (argc == 5) {
      std::string argv3 = argv[3];
      if (argv3 != "--solution") {
        std::cout << "Unexpected keyword \"" << argv3
                  << "\" at 3. argument. (Run w/o arguments for help.)" << std::endl;
        return EXIT_FAILURE;
      }
      std::string argv4 = argv[4];
      solver.write_tour_to_file(opt, argv4);
    }
  } catch (const std::runtime_error &error) {
    std::cout << error.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
