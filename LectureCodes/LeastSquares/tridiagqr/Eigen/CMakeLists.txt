project(tridiagqr)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../givensrot/Eigen
	)
