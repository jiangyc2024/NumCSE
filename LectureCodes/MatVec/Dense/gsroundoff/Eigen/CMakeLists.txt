project(gsroundoff)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../gramschmidt/Eigen
	)
