cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_COMPILER /usr/bin/g++)
set(CMAKE_C_COMPILER /usr/bin/gcc)
project(mylib)

add_definitions(-std=gnu++11 -lmgl)

include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(MathGL2 2.0.0 REQUIRED)
include_directories(${MATHGL2_INCLUDE_DIRS})

add_library(mylib ${CMAKE_CURRENT_SOURCE_DIR}/mylibrary.cpp)
add_executable(main ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(main mylib)
target_link_libraries(main ${MATHGL2_LIBRARIES})
