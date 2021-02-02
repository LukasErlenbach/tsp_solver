# Exact Solver for the Traveling Salesman Problem (TSP)

_(last modified 01.02.2021 by Lukas Erlenbach, [LinkedIn profile](www.linkedin.com/in/lukas-erlenbach))_

## Introduction

The project contains a C++ solver for the well-known NP-hard TSP problem. Given an undirected graph and edge costs, it finds a minimal cost tour that visits every node exactly once. As the problem is NP-hard, computing an optimal solution is difficult. The here presented program can solve instances (given in the TSPLIB format) up to about 100 nodes in a few seconds using a Branch-and-Bound approach suggested by [1].

## Solving Approach

The solver uses a Branch-and-Bound approach as suggested by Volgenant and Jonker [1].

It can be shown, that a lower bound to the TSP problem can be computed using modified edge costs and 1-trees (this is known a Held-Karp lower bound). An upper bound is found during the branching, every time a computed 1-tree is 2-regular, i.e. it is a tour.

The branching is implemented by successively making edges forbidden or required. The order is processed branching nodes is determined by _best-choice_, i.e. the branching tree is implemented as priority queue.

For more details on the choices of parameters, modified edge costs, etc. please refer to [1]. Emphasis was put on a fast implementation and a minimal encoding of 1-trees.

## Installation

Make sure you have installed recent versions of:

    -    cmake (tested 3.19.3, >3.0 should work)
    -    g++ compiler (tested 10.2.0)
    - OR clang compiler (tested 11.0.1)
After cloning the repository create a new directory and run cmake inside it:

        mkdir build
        cd build
        cmake ../
        make

You can also tell cmake to use another compiler by calling:

        cmake -D CMAKE_CXX_COMPILER={YOUR_COMPILER_NAME} ../

This should compile the code and create an executable inside the build directory named _build/tsp_solver_.

## Usage

After compiling the program, it can be used to compute an minimal cost TSP tour wrt. rounded Euclidean edge costs for graphs given in the TSPLIB format. A few example instances are provided in the _instances/_ directory. The instances are taken from and the results can be checked at the [webpage of the Zuse Institute Berlin](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/stsp-sol.html). The program can be called via:

        ./build/tsp_solver --instance file.tsp [--solution opt_tour.tsp]

The first argument is mandatory and _file.tsp_ is expected to encode the problem instance in the TSPLIB format. The second argument is optional and if given, an optimal tour will be written to _opt_tour.tsp_.

The program writes the cost of an optimal tour to _std::cout_, additional logging information like current lower/upper bounds and state of the branching tree are written to _std::cerr_.

Example:

        ./build/tsp_solver --instance instances/rd100.tsp
Output:

        minimal TSP tour length: 7910

## References
[1] Volgenant and Jonker, "A branch and bound algorithm for the symmetric traveling salesman problem based on 1-tree relaxation", European Journal of Operational Research, 1982.