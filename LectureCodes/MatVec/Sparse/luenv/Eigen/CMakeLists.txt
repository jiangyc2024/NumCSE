cmake_minimum_required(VERSION 2.8)
project(luenv)

add_executable_numcse(luenv_main main.cpp)


include_directories(${DIRS}
	../../rowbandwidth/Eigen
	../../substenv/Eigen
	)

