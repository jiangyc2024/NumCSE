project(chebpolplot)
cmake_minimum_required(VERSION 2.8)

include_directories(../../chebpolmult/Eigen)
add_executable_numcse(main main.cpp)
