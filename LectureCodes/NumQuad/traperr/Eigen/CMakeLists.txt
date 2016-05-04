project(traperr)
cmake_minimum_required(VERSION 2.8)

include_directories(../../trapezoidal/Eigen/)
add_executable_numcse(main traperr.cpp)
