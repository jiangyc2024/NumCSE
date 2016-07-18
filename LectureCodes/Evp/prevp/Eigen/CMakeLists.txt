project(prevp)
cmake_minimum_required(VERSION 3.0.2)

include_directories(../../utils)
include_directories(../../prbuildA/Eigen)

add_executable_numcse(main main.cpp)
add_resources(../../resources/pagerank)
