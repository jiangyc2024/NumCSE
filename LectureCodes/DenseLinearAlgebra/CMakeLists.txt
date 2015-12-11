# define the executable name
set (EXECNAME ncse_denselinalg)
message(STATUS "Executable ${EXECNAME}")

# create betl executable
add_executable (${EXECNAME}
  ${NCSE_EIGEN_LECTURECODE_PATH}/DenseLinearAlgebra/main.cpp
  )

######################################################################
# Dynamic linker has to be told where to find libraries
# setenv DYLD_LIBRARY_PATH "/opt/intel/mkl/lib:/opt/intel/lib"
######################################################################

# Specify libraries to link with
target_link_libraries (${EXECNAME} ${ALL_LIBS} ${MKL_LIBRARIES})

# Specify defines for MKL
ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)

# No idea what this is good for
set_target_properties (${EXECNAME}
 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_PATH})
