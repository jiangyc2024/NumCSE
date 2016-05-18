option(WITH_MKL "Build with MKL" ON)

if(WITH_MKL)
	get_target_name_numcse(ncse_denselinalg target_name)

	# create betl executable
	add_executable_numcse(ncse_denselinalg main.cpp)

	######################################################################
	# Dynamic linker has to be told where to find libraries
	# setenv DYLD_LIBRARY_PATH "/opt/intel/mkl/lib:/opt/intel/lib"
	######################################################################

	# Specify libraries to link with
	target_link_libraries(${target_name} ${MKL_LIBRARIES})

	# Specify defines for MKL
	ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)
endif(WITH_MKL)