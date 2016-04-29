macro (add_executable_numcse _name)
	message(STATUS "Executable ${_name}")

	# create directory that will contain the binary
	set( TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

	string(REPLACE "${PROJECT_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY})

    # invoke built-in add_executable
    add_executable(${ARGV})
    if (TARGET ${_name})
    	set_target_properties(${_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}")
        target_link_libraries(${_name} ${MATHGL2_LIBRARIES})
    	#add_dependencies(${_name} Eigen)
    	#add_dependencies(${_name} MathGL)
    endif()
endmacro()