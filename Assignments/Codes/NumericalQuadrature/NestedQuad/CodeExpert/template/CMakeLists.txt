
project(NestedQuad_proj)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -fdiagnostics-color=always -std=c++11 -lpthread -lutil -ldl -Xlinker -export-dynamic")

find_package(Python3 COMPONENTS Interpreter Development NumPy)
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

set(SRCS gaussquad.hpp laserquad.hpp main.cpp matplotlibcpp.h)
add_executable(NestedQuad ${SRCS})
set_target_properties(NestedQuad PROPERTIES OUTPUT_NAME "a.out")
target_link_libraries (NestedQuad Eigen3::Eigen)
target_link_libraries (NestedQuad Python3::Python)
target_link_libraries (NestedQuad Python3::NumPy)
