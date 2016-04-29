option(WITH_MKL "Build with MKL" ON)

if(WITH_MKL)
	# define the executable name
	set (EXECNAME ncse_denselinalg)

	# create betl executable
	add_executable_numcse(${EXECNAME} main.cpp)

	######################################################################
	# Dynamic linker has to be told where to find libraries
	# setenv DYLD_LIBRARY_PATH "/opt/intel/mkl/lib:/opt/intel/lib"
	######################################################################

	# Specify libraries to link with
	target_link_libraries (${EXECNAME} ${MKL_LIBRARIES})

	# Specify defines for MKL
	ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)
endif(WITH_MKL)