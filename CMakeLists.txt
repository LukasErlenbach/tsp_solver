cmake_minimum_required(VERSION 3.0)

# set the project name
project(TSP_Solver)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS="-MMD -Werror -Wall -Wextra -pedantic -o")

# bring the headers into the project
include_directories(include)

# add all sources
file(GLOB SOURCES "src/*.cpp")

# add the executable
add_executable(tsp_solver ${SOURCES})
