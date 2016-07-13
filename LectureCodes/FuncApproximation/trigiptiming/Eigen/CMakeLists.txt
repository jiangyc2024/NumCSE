project(trigiptiming)
cmake_minimum_required(VERSION 2.8)

add_definitions(-O3)
include_directories(../../trigipequid/Eigen)
include_directories(../../trigpolycoeff/Eigen)
include_directories(../../../../Utils)

add_executable_numcse(main ./trigiptiming.cpp)
