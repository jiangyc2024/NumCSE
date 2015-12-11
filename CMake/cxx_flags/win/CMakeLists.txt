message(STATUS "Compiler settings: WIN")

set( CMAKE_CXX_FLAGS 
  "${CMAKE_CXX_FLAGS} -std=c++11 -g -D NO_INITIALIZER_LIST_SUPPORT" )
