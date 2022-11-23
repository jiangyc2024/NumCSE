project(numquad)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../../Utils/)

add_executable_numcse(main main.cpp)

add_executable_numcse(integrate integrate.cpp)
