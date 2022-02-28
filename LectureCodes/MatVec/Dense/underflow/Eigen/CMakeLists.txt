cmake_minimum_required(VERSION 2.8)
project(underflow)

add_executable_numcse(main underflow.cpp)
add_compile_definitions(_USE_MATH_DEFINES)
