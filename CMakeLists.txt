project (NumCSE)
cmake_minimum_required(VERSION 2.8)

#
# setup cmake
#
# enforce policy that target_link_libraries must refer
# to a valid target. This avoids misconception of the
# add_executable_numcse command in combination with a
# manual call to target_link_libraries
if(POLICY CMP0016)
  cmake_policy(SET CMP0016 NEW)
endif()

# use install rpath already on built time
# this is used to make the install name of the mathgl (shared)
# library (already) on built time relative instead absolute
# note: disabled since absolute paths are fine for now
#SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# set cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

# store binaries in ./bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# load macros
include(CMake/Macros/subdirlist.cmake)
include(CMake/Macros/add_all_subdirectories.cmake)
include(CMake/Macros/get_target_name_numcse.cmake)
include(CMake/Macros/add_executable_numcse.cmake)
include(CMake/Macros/add_resources.cmake)
include(CMake/Macros/set_eigen_fft_backend.cmake)

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

#set(EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_CXX_FLAGS="${EXTERNAL_PROJECT_CXX_FLAGS}")

#
# include external libraries
#
include(ExternalProject)

# eigen
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
      URL http://bitbucket.org/eigen/eigen/get/3.2.7.zip
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen
      INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen_install
      DOWNLOAD_NO_PROGRESS 1
      CMAKE_ARGS ${EXTERNAL_PROJECT_CMAKE_ARGS_PREFIX} -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/Eigen_install)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/Eigen_install/include/eigen3)
endif()

include_directories(${CMAKE_SOURCE_DIR}/MathGL)

# Python includes
# TODO Generic include maybe with FindPython.cmake, this one only works for
# my machine
set(PYHTONHOME /usr/local/Cellar/python/3.7.3/Frameworks/Python.framework/Versions/3.7)
set(NUMPYHOME /usr/local/lib/python3.7/site-packages/numpy/core)
include_directories(${PYHTONHOME}/include/python3.7m)
include_directories(${NUMPYHOME}/include)


include_directories(${CMAKE_SOURCE_DIR}/Utils) # oh dear

#
# include subdirectories
#
add_subdirectory(LectureCodes)
add_subdirectory(Assignments)
add_subdirectory(Utils)
add_subdirectory(CppTutorial)
add_subdirectory(FFT/BackendTimings)
