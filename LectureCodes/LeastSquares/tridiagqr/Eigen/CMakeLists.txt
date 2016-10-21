cmake_minimum_required(VERSION 2.8)
project(tridiagqr)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../givensrot/Eigen
	)
