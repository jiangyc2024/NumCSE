project (NumCSE)
cmake_minimum_required(VERSION 2.8)

#
# setup cmake
#
# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

# store binaries in ./bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# load macros
include(CMake/Macros/subdirlist.cmake)
include(CMake/Macros/add_all_subdirectories.cmake)
include(CMake/Macros/add_executable_numcse.cmake)

#
# setup compiler
#
# use c++11
if(CMAKE_VERSION VERSION_GREATER 3.1.0 OR CMAKE_VERSION VERSION_EQUAL 3.1.0)
    # only valid for new platforms
    set(CMAKE_CXX_STANDARD 11)
else()
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	# for older cmake versions
	# (note, this CXX flag is only valid for clang and gcc, for MSVC it is not needed)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
endif()

#
# include external libraries
#
# mathgl
# todo: download
find_package(MathGL2)
if (${MATHGL2_FOUND})
	# patch mgl2/config.h
	file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/mathgl_patched_headers/mgl2)
	file(READ ${MATHGL2_INCLUDE_DIRS}/mgl2/config.h MATHGL_PATCHED_CONFIG_H)
	STRING(REGEX REPLACE "#define MGL_HAVE_TYPEOF[ \t]+1" "#define MGL_HAVE_TYPEOF 0 // patched"
	                     MATHGL_PATCHED_CONFIG_H ${MATHGL_PATCHED_CONFIG_H})
	STRING(REGEX REPLACE "#define MGL_HAVE_C99_COMPLEX[ \t]+1" "#define MGL_HAVE_C99_COMPLEX 0 // patched"
	                     MATHGL_PATCHED_CONFIG_H ${MATHGL_PATCHED_CONFIG_H})

	file(WRITE ${PROJECT_BINARY_DIR}/mathgl_patched_headers/mgl2/config.h ${MATHGL_PATCHED_CONFIG_H} )
	include_directories(${PROJECT_BINARY_DIR}/mathgl_patched_headers/) # use patched config.h for mathgl
endif()
#  include mathgl headers
include_directories(${MATHGL2_INCLUDE_DIRS})

# eigen
find_package(Eigen3)
if (${EIGEN3_FOUND})
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# mathgl figure
# try to find Figure with FindFigure.cmake
find_package(Figure QUIET)

if(FIGURE_FOUND)
  set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
  message(STATUS "Function GET_MODULES: Included Figure directory in variable DIRS")
  set(LIBS ${LIBS} ${FIGURE_LIBRARY})
  message(STATUS "Function GET_MODULES: Included Figure library in variable LIBS")
# case if Figure is not found by FindFigure.cmake - try to get it from MathGL/FigureClass
else() 
  set(FIGURE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/MathGL/FigureClass/src) # directory which should contains the source files
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

include_directories(${FIGURE_INCLUDE_DIR})

include_directories(${FIGURE_INCLUDE_DIR}/../)

include_directories(${CMAKE_SOURCE_DIR}/LectureCodes/FunctionCollection) # oh dear

#
# include subdirectories
#
add_subdirectory(LectureCodes)