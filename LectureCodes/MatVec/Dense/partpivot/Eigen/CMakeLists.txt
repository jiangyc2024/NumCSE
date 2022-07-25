project(partpivot)

add_executable_numcse(main main.cpp)

include_directories(${DIRS}
	../../gausselimsolve/Eigen
	../../lufak/Eigen
	)
