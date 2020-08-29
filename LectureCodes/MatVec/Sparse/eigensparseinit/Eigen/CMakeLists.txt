cmake_minimum_required(VERSION 2.8)
project(sparsetiming)

add_executable_numcse(main sparsetiming.cpp)

include_directories( ${CMAKE_SOURCE_DIR}/Utils )