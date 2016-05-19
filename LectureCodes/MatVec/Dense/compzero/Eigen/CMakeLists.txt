cmake_minimum_required(VERSION 2.8)
project(compzeros)

include(../../../../GetModules.cmake)
get_modules("Eigen Figure")

include_directories(${DIRS})
add_executable(main main.cpp)
target_link_libraries(main ${LIBS})
