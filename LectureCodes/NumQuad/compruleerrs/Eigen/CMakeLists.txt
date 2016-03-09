project(compruleerrs)
cmake_minimum_required(VERSION 2.8)

include(../../../GetModules.cmake)
get_modules("Eigen Figure")

include_directories(${DIRS})
add_executable(main compruleerrs.cpp)
target_link_libraries(main ${LIBS})
