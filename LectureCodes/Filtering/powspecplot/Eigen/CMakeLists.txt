project(powspecplot)
cmake_minimum_required(VERSION 2.8)

include_directories(../../signalgen/Eigen/)
add_executable_numcse(main powspecplot.cpp)
