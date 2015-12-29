##############################################################
# this CMakeLists.txt contains the function 'get_dirs':      #
# Usage:                                                     #
#   get_dirs("Eigen and MathGL")                             #
#   include_directories(${DIRS})                             #
#   add_exectutable(main main.cpp)                           #
#   target_link_libraries(main ${LIBS})                      #
##############################################################

cmake_minimum_required(VERSION 2.8)
project(DefaultSettings)

## set path to the Find<LIB>.cmake files ##
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../CMake/Modules)
message("-- Current module path: ${CMAKE_MODULE_PATH}")

## to find Eigen, MathGL and Figure: ##
# PRE: valid variable 
# POST: if variable contains "eigen", "figure"/"plot" or "mathgl" the according
#       directories/libraries will be included in DIRS/LIBS
function(get_dirs includes)

  # convert input string to lower case to make the function case insensitive
  string(TOLOWER "${includes}" includes_lower)

  # check whether the terms "eigen", "figure, "plot" or "mathgl" are in the input
  string(REGEX MATCH "(eigen)" include_eigen "${includes_lower}")
  string(REGEX MATCH "(mathgl)" include_mathgl "${includes_lower}")
  string(REGEX MATCH "(figure)" include_figure "${includes_lower}")
  string(REGEX MATCH "(plot)" include_plot "${includes_lower}")

  # if "eigen" was in the input then find it and add the directory to DIRS
  if (include_eigen)
    find_package(Eigen3 REQUIRED)
    set(DIRS ${DIRS} ${EIGEN3_INCLUDE_DIR})
    message("-- Function GET_DIRS : Included Eigen3 directory in variable DIRS")
  endif()

  # if "figure" or "plot" was in the input then find it and add the directory
  # to DIRS and the library to LIBS
  if (include_figure OR include_plot)
    add_definitions(-lFigure)
    set(include_mathgl "mathgl") # Figure needs MathGL!
    find_package(Figure REQUIRED)
    set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
    message("-- Function GET_DIRS : Included Figure directory in variable DIRS")
    set(LIBS ${LIBS} ${FIGURE_LIBRARY})
    message("-- Function GET_DIRS : Included Figure library in variable LIBS")
  endif()

  # if "mathgl" was in the input then find it and add the directories to DIRS
  # and the libraries to LIBS
  if (include_mathgl)
    add_definitions(-std=gnu++11 -lmgl)  # MathGL needs the GNU compiler 
    find_package(MathGL2 2.0.0 REQUIRED)
    set(DIRS ${DIRS} ${MATHGL2_INCLUDE_DIRS})
    message("-- Function GET_DIRS : Included MathGL2 directory in variable DIRS")
    set(LIBS ${LIBS} ${MATHGL2_LIBRARIES})
    message("-- Function GET_DIRS : Included MathGL2 library in variable LIBS")
  endif()

  # make variables DIRS and LIBS available for files which include this CMakeLists
  set(DIRS ${DIRS} PARENT_SCOPE)
  set(LIBS ${LIBS} PARENT_SCOPE)

endfunction(get_dirs)
