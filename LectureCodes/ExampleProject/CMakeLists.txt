cmake_minimum_required(VERSION 2.8)
project(ExampleProject)

include(../CMakeLists.txt)
get_dirs("Eigen Figure")
include_directories(${DIRS})

add_executable(main main.cpp)
target_link_libraries(main ${LIBS})
