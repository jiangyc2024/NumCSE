project(polytiming)
cmake_minimum_required(VERSION 2.8)

set(LOCDIR ${CMAKE_CURRENT_SOURCE_DIR})
include(../../../GetModules.cmake)
get_modules("Eigen Figure")
include_directories(${DIRS} 
  ${LOCDIR}/../../AitkenNeville/Eigen/ # ANipoleval
  ${LOCDIR}/../../barycentricformula/Eigen/ # intpolval
  ${LOCDIR}/../../intpolyval_lag/Eigen/ # intpolval_lag
  ${LOCDIR}/../../ipoleval/Eigen/ # ipoleval
  ${LOCDIR}/../../../FunctionCollection/ # timer
   )
add_definitions(-O3)
add_executable(main polytiming.cpp)
target_link_libraries(main ${LIBS} Timer)
