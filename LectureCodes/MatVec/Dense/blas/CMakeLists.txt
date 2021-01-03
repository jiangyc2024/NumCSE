cmake_minimum_required(VERSION 3.02)

#using the built-in FindBLAS.cmake module here, not something in NumCSE/CMake/Modules
find_package(BLAS)
project(blas)

if(BLAS_FOUND)
	message("Found blas (installed with ${BLAS_LIBRARIES})")
	
	set(Sources BlasTimings.cpp ColumnMajorMatrix.cpp ColumnMajorMatrix_Multiplication_dot.cpp ColumnMajorMatrix_Multiplication_gemm.cpp ColumnMajorMatrix_Multiplication_gemv.cpp ColumnMajorMatrix_Multiplication_noblas.cpp)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -D_MAC_OS")
	
	add_executable_numcse(mult_g++O3 ${Sources})
	get_target_name_numcse(mult_g++O3 target_name)
	target_link_libraries(${target_name} ${BLAS_LIBRARIES})

	add_executable_numcse(test blas1.cpp)
	get_target_name_numcse(test target_name)
	target_link_libraries(${target_name} ${BLAS_LIBRARIES})
else()
	message("Blas package not found. Please install by 'apt-get install libblas-dev' (Ubuntu) or 'dnf install openblas-devel' (Fedora)")
endif()
