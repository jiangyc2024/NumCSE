cmake_minimum_required(VERSION 2.8)

project(FFT_filter)

## definitions needed for MathGL ##
add_definitions(-std=gnu++11 -lmgl)
#add_definitions(-O3) #optimizing

## find MathGL and Eigen ##
# set path to FindEigen3.cmake and FindMathGL2.cmake
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../CMake/Modules)

find_package(Eigen3 REQUIRED)
find_package(MathGL2 2.0.0 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR} ${MATHGL2_INCLUDE_DIRS})

## build executable and link libraries ##
add_executable(main filter.cpp)
target_link_libraries(main ${MATHGL2_LIBRARIES})

