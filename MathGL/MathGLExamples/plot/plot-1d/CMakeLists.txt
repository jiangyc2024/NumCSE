cmake_minimum_required(VERSION 2.8)
project(plot)

## correct compiler options for MathGL ##
add_definitions(-std=gnu++11 -lmgl)

## find MathGL2 ##
# this is a bit nasty but I dont know a better way to get to 
# the correct directory (../../cmake) yet
message("Current project dir: ${PROJECT_SOURCE_DIR}")
get_filename_component(PRT ${PROJECT_SOURCE_DIR} PATH)
get_filename_component(PPRT ${PRT} PATH)
get_filename_component(PPPRT ${PPRT} PATH)
# set correct path to FindMathGL2.cmake
set(CMAKE_MODULE_PATH "${PPPRT}/cmake/")
message("Module path: ${CMAKE_MODULE_PATH}")
# find package using ../../cmake/FindMathGL2.cmake
find_package(MathGL2 2.0.0 REQUIRED)
include_directories(${MATHGL2_INCLUDE_DIRS})

## build executable ##
add_executable(main plot.cpp)
target_link_libraries(main ${MATHGL2_LIBRARIES})
