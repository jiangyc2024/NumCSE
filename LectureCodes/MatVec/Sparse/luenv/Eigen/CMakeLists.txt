project(luenv)

add_executable_numcse(main main.cpp)


include_directories(${DIRS}
	../../rowbandwidth/Eigen
	../../substenv/Eigen
	)

