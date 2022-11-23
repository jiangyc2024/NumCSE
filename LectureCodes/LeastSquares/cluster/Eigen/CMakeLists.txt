cmake_minimum_required(VERSION 2.8)
project(cluster)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../pas/Eigen
	../../lloydmax/Eigen)
	
