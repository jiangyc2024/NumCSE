project(ProblemSheet1)
cmake_minimum_required(VERSION 2.6.2)

# Tested with g++ and clang
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    # Use -O3 for strong optimization. We use stdc++11. Compile with Wall to get all warnings.
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -O3 -Wextra -Wno-deprecated-declarations -pedantic" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# If your compiler already finds and includes Eigen you can safely remove the following two lines
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

### Problem Sheet 2
# make gramschmidt: problem 1
set(SRCS solveLyapunov.cpp)
add_executable(solveLyapunov ${SRCS}) 

# make gramschmidt: problem 2
set(SRCS block_lu_decomp.cpp)
add_executable(block_lu_decomp ${SRCS}) 

# make gramschmidt: problem 3
set(SRCS banded_matrix.cpp)
add_executable(banded_matrix ${SRCS}) 

# make gramschmidt: problem 3
set(SRCS solvepermb.cpp)
add_executable(solvepermb ${SRCS}) 

### Problem Sheet 3
# Problem 5
set(SRCS rankoneinvit.cpp)
add_executable(rankoneinvit ${SRCS}) 

# Problem 5
set(SRCS impedancemap.cpp)
add_executable(impedancemap ${SRCS}) 

# Problem 5
set(SRCS impedancemap_stack.cpp)
add_executable(impedancemap_stack ${SRCS}) 
