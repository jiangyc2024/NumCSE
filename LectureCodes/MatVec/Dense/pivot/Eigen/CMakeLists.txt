cmake_minimum_required(VERSION 2.8)
project(pivot)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../gausselimsolve/Eigen
	../../lufak/Eigen
	)
