##############################################################
# this CMake-file contains the function 'get_modules':       #
# Usage:                                                     #
#   get_modules("Eigen and MathGL")                          #
#   include_directories(${DIRS})                             #
#   add_exectutable(main main.cpp)                           #
#   target_link_libraries(main ${LIBS})                      #
##############################################################

cmake_minimum_required(VERSION 2.8)
project(GetModules)

## set path to the Find<LIB>.cmake files ##
set(LECTURE_CODES_DIR ${CMAKE_CURRENT_LIST_DIR}) # paths to Figure will be relative to this path
set(CMAKE_MODULE_PATH ${LECTURE_CODES_DIR}/../CMake/Modules)
message(STATUS "Current module path: ${CMAKE_MODULE_PATH}")

## to find Eigen, MathGL and Figure: ##
# PRE: valid variable 
# POST: if variable contains "eigen", "figure" or "mathgl" the according
#       directories/libraries will be included in DIRS/LIBS
function(get_modules includes)

  # convert input string to lower case to make the function case insensitive
  string(TOLOWER "${includes}" includes_lower)

  # check whether the terms "eigen", "figure" or "mathgl" are in the input
  string(REGEX MATCH "(eigen)" include_eigen "${includes_lower}")
  string(REGEX MATCH "(mathgl)" include_mathgl "${includes_lower}")
  string(REGEX MATCH "(figure)" include_figure "${includes_lower}")

  add_definitions(-std=gnu++11) # using gnu and not c++11 because mathgl needs the gnu standard

  # ---------------------------- FIGURE -------------------------------- #
  # if "figure" was in the input then find it and add the directory to DIRS and the library to LIBS
  if (include_figure)
  
    set(include_mathgl true) # Figure needs MathGL
    add_definitions(-lmgl) # compiler definitions

    # try to find Figure with FindFigure.cmake
    find_package(Figure QUIET)
    
    if(FIGURE_FOUND)
      set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
      message(STATUS "Function GET_MODULES: Included Figure directory in variable DIRS")
      set(LIBS ${LIBS} ${FIGURE_LIBRARY})
      message(STATUS "Function GET_MODULES: Included Figure library in variable LIBS")
    # case if Figure is not found by FindFigure.cmake - try to get it from MathGL/FigureClass
    else() 
      set(FIGURE_INCLUDE_DIR ${LECTURE_CODES_DIR}/../MathGL/FigureClass/src) # directory which should contains the source files
      message(STATUS "Trying to get it from ${FIGURE_INCLUDE_DIR} ...")

      # check if necessary files exist in MathGL/FigureClass
      foreach(FIGURE_FILE ${FIGURE_FILE_LIST})
        if (NOT EXISTS ${FIGURE_INCLUDE_DIR}/${FIGURE_FILE})
          message(FATAL_ERROR "Could not find necessary files to build Figure library! Try cloning the git repo again or contact someone.")
        endif()
      endforeach()
      message(STATUS "Found necessary Figure files: ${FIGURE_INCLUDE_DIR}")

      # mark Figure library to be built in 'make' stage: we need Eigen and MathGL so fnd the packages and libraries and include them
      unset(inlude_eigen) # we're already getting those packages here, so we don't need to do it again later
      unset(include_mathgl)
      find_package(Eigen3 REQUIRED)
      find_package(MathGL2 2.0.0 REQUIRED)
      include_directories(${EIGEN3_INCLUDE_DIR} ${MATHGL2_INCLUDE_DIRS})

      add_library(Figure STATIC ${FIGURE_INCLUDE_DIR}/figure.cpp)
      target_link_libraries(Figure ${MATHGL2_LIBRARIES})

      # as libFigure.a was not built yet (this happens when '$ make' is executed) we need to call
      # target_link_libraries(main Figure) and *not* target_link_libraries(main libFigure.a)
      set(DIRS ${DIRS} ${EIGEN3_INCLUDE_DIR} ${MATHGL2_INCLUDE_DIRS} ${FIGURE_INCLUDE_DIR})
      set(LIBS ${LIBS} ${MATHGL2_LIBRARIES} Figure)
      message(STATUS "Function GET_MODULES: Included Eigen3, MathGL2 and Figure directories in variable DIRS")
      message(STATUS "Function GET_MODULES: Included MathGL2 library in variable LIBS")
      message(STATUS "Function GET_MODULES: Figure library marked to be built in build stage")

    endif()
    
  endif()

  # ---------------------------- EIGEN --------------------------------- #
  # if "eigen" was in the input then find it and add the directory to DIRS
  if (include_eigen)
    find_package(Eigen3 REQUIRED)
    set(DIRS ${DIRS} ${EIGEN3_INCLUDE_DIR})
    message(STATUS "Function GET_MODULES: Included Eigen3 directory in variable DIRS")
  endif()

  # ---------------------------- MATHGL -------------------------------- #
  # if "mathgl" was in the input then find it and add the directories to DIRS
  # and the libraries to LIBS
  if (include_mathgl)
    add_definitions(-lmgl)  # MathGL needs the GNU compiler 
    find_package(MathGL2 2.0.0 REQUIRED)
    set(DIRS ${DIRS} ${MATHGL2_INCLUDE_DIRS})
    message(STATUS "Function GET_MODULES: Included MathGL2 directory in variable DIRS")
    set(LIBS ${LIBS} ${MATHGL2_LIBRARIES})
    message(STATUS "Function GET_MODULES: Included MathGL2 library in variable LIBS")
  endif()

  # make variables DIRS and LIBS available for files which include this CMakeLists
  set(DIRS ${DIRS} PARENT_SCOPE)
  set(LIBS ${LIBS} PARENT_SCOPE)

endfunction(get_modules)
