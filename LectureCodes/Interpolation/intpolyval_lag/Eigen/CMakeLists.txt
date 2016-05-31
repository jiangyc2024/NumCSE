project(intpolyval_lag)
cmake_minimum_required(VERSION 2.8)

include_directories(../../lagrangepoly/Eigen/) 
add_executable_numcse(main main.cpp)
