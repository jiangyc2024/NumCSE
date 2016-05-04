include(CMake/Macros/list_replace.cmake)

macro (add_executable_numcse _name)
    # cmake targets need to be unique
    # as such we prefix every executable with the current project name
    SET(parsed_name "${PROJECT_NAME}_${_name}") # add prefix
    set(argv_parsed ${ARGV})
    LIST_REPLACE(argv_parsed 0 "${parsed_name}")

	# create directory that will contain the binary
	set( TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

	string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY})

    # invoke built-in add_executable
    add_executable(${argv_parsed})
    if (TARGET ${parsed_name})
    	set_target_properties(${parsed_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}")
        target_link_libraries(${parsed_name} ${MATHGL2_LIBRARIES} Figure)
    	#add_dependencies(${parsed_name} Eigen)
    	#add_dependencies(${parsed_name} MathGL)
    endif()
endmacro()