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
set(LIBS_PATH ${LECTURE_CODES_DIR}/../third_party)

# create libs folder
file(MAKE_DIRECTORY ${LIBS_PATH})

# external lib module
include(ExternalProject)

message(STATUS "Current module path: ${CMAKE_MODULE_PATH}")

## include flags
# include(${LECTURE_CODES_DIR}/../CMake/cxx_flags/CMakeLists.txt)

function(link_modules includes)


endfunction(link_modules)

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
    set(include_eigen true) # Figure needs Eigen

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

   
      add_library(Figure STATIC ${FIGURE_INCLUDE_DIR}/figure.cpp)
      
      add_dependencies(Figure Eigen)
      add_dependencies(Figure MathGL) 

      # as libFigure.a was not built yet (this happens when '$ make' is executed) we need to call
      # target_link_libraries(main Figure) and *not* target_link_libraries(main libFigure.a)
      set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
      set(LIBS ${LIBS} Figure)

    endif()
    
  endif()

  # ---------------------------- EIGEN --------------------------------- #
  # if "eigen" was in the input then find it and add the directory to DIRS
  if (include_eigen)
    find_package(Eigen3)

    if (NOT EIGEN3_FOUND)

	#  download if not on the local system
	message("-- Download Eigen3 to ${LIBS_PATH}/Eigen")
	ExternalProject_Add(
	    Eigen
	    URL http://bitbucket.org/eigen/eigen/get/3.2.7.zip
 	    URL_MD5 724e02b1b80c4b1c3b455384e8a96e90 #prevents re-downloading
	    SOURCE_DIR ${LIBS_PATH}/Eigen
	    DOWNLOAD_DIR ${LIBS_PATH}/temp
	    CONFIGURE_COMMAND ""
	    BUILD_COMMAND ""
	    INSTALL_COMMAND "")
	set(EIGEN3_INCLUDE_DIR ${LIBS_PATH}/Eigen)

    endif()


    set(DIRS ${DIRS} ${EIGEN3_INCLUDE_DIR})
    message(STATUS "Function GET_MODULES: Included Eigen3 directory in variable DIRS")
  endif()

  # ---------------------------- MATHGL -------------------------------- #
  # if "mathgl" was in the input then find it and add the directories to DIRS
  # and the libraries to LIBS
  if (include_mathgl)
    add_definitions(-lmgl)  # MathGL needs the GNU compiler 
    find_package(MathGL2 2.0.0)

    if (NOT MATHGL2_FOUND)

	#  download if not on the local system
	message("-- Download MathGl to ${LIBS_PATH}/mathgl_install")
	ExternalProject_Add(
	    MathGL
	    URL http://downloads.sourceforge.net/mathgl/mathgl-2.3.3.tar.gz
	    URL_MD5 c37d6f42d4897675bf89fae635aa6868  #prevents re-downloading
	    SOURCE_DIR ${LIBS_PATH}/mathgl_source
	    BINARY_DIR ${LIBS_PATH}/mathgl_binary
	    CMAKE_ARGS -DCMAKE_CXX_STANDARD=11 -Denable-openmp=OFF -DMGL_HAVE_TYPEOF=0 -DMGL_HAVE_C99_COMPLEX=0 -DMGL_LIB_INSTALL_DIR=${LIBS_PATH}/mathgl_install/lib/ -DMGL_CGI_PATH=${LIBS_PATH}/mathgl_install/share/mathgl -DCMAKE_INSTALL_PREFIX=${LIBS_PATH}/mathgl_install
	    INSTALL_DIR ${LIBS_PATH}/mathgl_install
	    DOWNLOAD_DIR ${LIBS_PATH}/temp
	    )
	set(MATHGL2_INCLUDE_DIRS ${LIBS_PATH}/mathgl_install/include)
	set(MATHGL2_LIBRARIES "mgl")
    endif()

    set(DIRS ${DIRS} ${MATHGL2_INCLUDE_DIRS})
    message(STATUS "Function GET_MODULES: Included MathGL2 directory in variable DIRS")
    set(LIBS ${LIBS} ${MATHGL2_LIBRARIES})
    message(STATUS "Function GET_MODULES: Included MathGL2 library in variable LIBS")

 
  endif()

  # make variables DIRS and LIBS available for files which include this CMakeLists
  set(DIRS ${DIRS} PARENT_SCOPE)
  set(LIBS ${LIBS} PARENT_SCOPE)

endfunction(get_modules)

