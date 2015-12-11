# define the executable name
set (EXECNAME ncse_eigentutorial)
message(STATUS "Executable ${EXECNAME}")

# create betl executable
add_executable (${EXECNAME}
  ${NCSE_EIGEN_LECTURECODE_PATH}/EigenTutorial/main.cpp
  )

# Specify libraries to link with
target_link_libraries (${EXECNAME} ${ALL_LIBS})

# No idea what this is good for
set_target_properties (${EXECNAME}
 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_PATH})
