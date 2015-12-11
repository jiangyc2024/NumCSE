cmake_minimum_required(VERSION 2.8)

## add compiler flags needed for MathGL and O3 optimisation ##
add_definitions(-std=gnu++11 -lmgl -O3)

project(Interpolation)

message("Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ../../cmake/ )
message("Current module path: ${CMAKE_MODULE_PATH}")
find_package(MathGL2 2.0.0 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(${MATHGL2_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
add_library(INTP ${CMAKE_CURRENT_SOURCE_DIR}/interpol.cpp)

add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(main INTP)
target_link_libraries(main ${MATHGL2_LIBRARIES})
