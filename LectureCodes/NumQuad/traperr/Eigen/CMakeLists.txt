project(traperr)
cmake_minimum_required(VERSION 2.8)

include(../../../GetModules.cmake)

get_modules("Eigen Figure")

include_directories(${DIRS} ../../trapezoidal/Eigen/)
add_executable(main traperr.cpp)
target_link_libraries(main ${LIBS})
