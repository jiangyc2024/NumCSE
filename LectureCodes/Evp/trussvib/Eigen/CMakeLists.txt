project(trussvib)
cmake_minimum_required(VERSION 2.8)

include_directories(../../bridge/Eigen)
include_directories(../../utils/)

add_executable_numcse(main main.cpp)
