project(polytiming)
cmake_minimum_required(VERSION 2.8)

include_directories(
  ../../AitkenNeville/Eigen/ # ANipoleval
  ../../barycentricformula/Eigen/ # intpolval
  ../../intpolyval_lag/Eigen/ # intpolval_lag
  ../../ipoleval/Eigen/ # ipoleval
  ../../../../Utils/ # timer
  ../../lagrangepoly/Eigen/ # needed by intpolval_lag.hpp
  )
#add_definitions(-O3) # todo add warning if CMAKE_RELEASE_TYPE!=RELEASE instead of enforcing O3
add_executable_numcse(main polytiming.cpp)
