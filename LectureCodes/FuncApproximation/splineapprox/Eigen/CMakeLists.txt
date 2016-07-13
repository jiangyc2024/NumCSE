project(splineapprox)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../../Utils/)
add_executable_numcse(main splineapproxdriver.cpp)
