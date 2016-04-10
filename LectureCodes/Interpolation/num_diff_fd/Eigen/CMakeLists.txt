project(num_diff_fd)
cmake_minimum_required(VERSION 2.8)

include(../../../GetModules.cmake)
get_modules("Figure")
include_directories(${DIRS})
add_executable(main numericaldifferentiation.cpp)
target_link_libraries(main ${LIBS})
