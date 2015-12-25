cmake_minimum_required(VERSION 2.8)
project(FigureClass)

add_definitions(-std=gnu++11)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules)
find_package(Eigen3 REQUIRED)
find_package(MathGL2 2.0.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR} ${MATHGL2_INCLUDE_DIRS})

add_library(Figure figure.cpp)
add_executable(main main.cpp)
target_link_libraries(main Figure ${MATHGL2_LIBRARIES})
