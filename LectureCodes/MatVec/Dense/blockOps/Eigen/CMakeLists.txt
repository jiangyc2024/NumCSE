cmake_minimum_required(VERSION 2.8)
project(blockOps)

include(../../../../GetModules.cmake)
get_modules("Eigen")

include_directories(${DIRS})
add_executable(main main.cpp)
