project(prstochstim)
cmake_minimum_required(VERSION 3.0.2)

include_directories(../../utils)

add_executable_numcse(main main.cpp)
add_resources(../../resources/pagerank)
