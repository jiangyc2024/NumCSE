project(trigpolyval)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../Interpolation/barycentricformula/Eigen)

add_executable_numcse(main main.cpp)
