cmake_minimum_required(VERSION 2.8)
project(arrowsys)

include(../../../../GetModules.cmake)
get_modules("Eigen Figure")

include_directories(${DIRS}
  ../../../../../Utils/ # timer
  )
add_executable(main main.cpp)
target_link_libraries(main ${LIBS})
