project(WeightedGaussQuadrature)
cmake_minimum_required(VERSION 2.8)

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

# enable some warnings
add_definitions(-pedantic -Wall)

# disable some warnings
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # (this is only applicable on clang)            
    # ignore some mathgl warnings
    add_definitions(-Wno-return-type-c-linkage -Wno-keyword-macro -Wno-missing-braces)
    #set(EXTERNAL_PROJECT_CXX_FLAGS ${EXTERNAL_PROJECT_CXX_FLAGS} -Wno-return-type-c-linkage)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # disable some warnings on gcc
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6)
        add_definitions(-Wno-ignored-attributes -Wno-misleading-indentation)
        #set(EXTERNAL_PROJECT_CXX_FLAGS ${EXTERNAL_PROJECT_CXX_FLAGS} -Wno-ignored-attributes -Wno-misleading-indentation)
    endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_definitions(-Wno-deprecated-declarations -Wno-unknown-pragmas)
    #set(EXTERNAL_PROJECT_CXX_FLAGS ${EXTERNAL_PROJECT_CXX_FLAGS} -Wno-deprecated-declarations)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake-modules/")

include(ExternalProject)

find_package(Eigen3)

if (${EIGEN3_FOUND})
  include_directories(${EIGEN3_INCLUDE_DIR})

  add_custom_target(Eigen) # dependency dummy
else()
  SET(DOWNLOADING_EIGEN ON)
  #  if not found system wide download
  message("-- Downloading Eigen3")
  ExternalProject_Add(
      Eigen
      URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
      DOWNLOAD_NO_PROGRESS 1
      CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)
endif()

find_package(MathGL2)
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  # If mathgl is installed, we do not trust that it is build using clang (in 99% of the cases,
  # it is not). Therefore, we *always* build mathgl when we use clang
  set(MATHGL2_FOUND OFF)
endif()
if (${MATHGL2_FOUND})
  add_custom_target(MathGL) # dependency dummy

  # patch mgl2/config.h
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/mathgl_patched_headers/mgl2)
  file(READ ${MATHGL2_INCLUDE_DIRS}/mgl2/config.h MATHGL_PATCHED_CONFIG_H)
  STRING(REGEX REPLACE "#define MGL_HAVE_TYPEOF[ \t]+1" "#define MGL_HAVE_TYPEOF 0 // patched"
                       MATHGL_PATCHED_CONFIG_H ${MATHGL_PATCHED_CONFIG_H})
  STRING(REGEX REPLACE "#define MGL_HAVE_C99_COMPLEX[ \t]+1" "#define MGL_HAVE_C99_COMPLEX 0 // patched"
                       MATHGL_PATCHED_CONFIG_H ${MATHGL_PATCHED_CONFIG_H})

  file(WRITE ${PROJECT_BINARY_DIR}/mathgl_patched_headers/mgl2/config.h ${MATHGL_PATCHED_CONFIG_H} )
  include_directories(${PROJECT_BINARY_DIR}/mathgl_patched_headers/) # use patched config.h for mathgl
else()
  set(DOWNLOADING_MGL ON)
  message("-- Downloading MathGl")
  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mathgl_install)
    ExternalProject_Add(
        MathGL
        URL http://downloads.sourceforge.net/mathgl/mathgl-2.3.3.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathgl_source
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathgl_binary
        DOWNLOAD_NO_PROGRESS 1
        CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_CXX_STANDARD=11 -Denable-openmp=OFF -DMGL_HAVE_TYPEOF=0 -DMGL_HAVE_C99_COMPLEX=0 -DMGL_LIB_INSTALL_DIR=${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/lib/ -DMGL_CGI_PATH=${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/share/mathgl -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/mathgl_install -DCMAKE_INSTALL_NAME_DIR=@rpath
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/mathgl_install
        )
  else()
    add_custom_target(MathGL)
    message("Skipped download of mathgl because its already present. Delete `${CMAKE_CURRENT_BINARY_DIR}/mathgl_install` to redownload it.")
  endif()

  set(MATHGL2_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/include)
  set(MATHGL2_LIBRARY_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}mgl${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MATHGL2_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/mathgl_install/lib/${MATHGL2_LIBRARY_NAME}")

endif()
include_directories(${MATHGL2_INCLUDE_DIRS})

find_package(Figure QUIET)

if(FIGURE_FOUND)
  set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
  message(STATUS "Function GET_MODULES: Included Figure directory in variable DIRS")
  set(LIBS ${LIBS} ${FIGURE_LIBRARY})
  message(STATUS "Function GET_MODULES: Included Figure library in variable LIBS")
# case if Figure is not found by FindFigure.cmake - try to get it from MathGL/FigureClass
else() 
  set(FIGURE_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/figure) # directory which should contains the source files
  message(STATUS "Trying to get it from ${FIGURE_INCLUDE_DIR} ...")

  # check if necessary files exist in MathGL/FigureClass
  foreach(FIGURE_FILE ${FIGURE_FILE_LIST})
    if (NOT EXISTS ${FIGURE_INCLUDE_DIR}/${FIGURE_FILE})
      message(FATAL_ERROR "Could not find necessary files to build Figure library! Try cloning the git repo again or contact someone.")
    endif()
  endforeach()
  message(STATUS "Found necessary Figure files: ${FIGURE_INCLUDE_DIR}")

  add_library(Figure STATIC ${FIGURE_INCLUDE_DIR}/figure.cpp)
  #add_library(Figure SHARED ${FIGURE_INCLUDE_DIR}/figure.cpp)
  #target_link_libraries(Figure ${MATHGL2_LIBRARIES})
  
  add_dependencies(Figure Eigen)
  add_dependencies(Figure MathGL) 

  # as libFigure.a was not built yet (this happens when '$ make' is executed) we need to call
  # target_link_libraries(main Figure) and *not* target_link_libraries(main libFigure.a)
  set(DIRS ${DIRS} ${FIGURE_INCLUDE_DIR})
  set(LIBS ${LIBS} Figure)

endif()
include_directories(${FIGURE_INCLUDE_DIR})
include_directories(${FIGURE_INCLUDE_DIR}/..)

### WeightedGaussQuadrature

# make quadU
set(SRCS quadU.cpp)
add_executable(quadU ${SRCS})
target_link_libraries(quadU Figure ${MATHGL2_LIBRARIES})

