include(CMake/Macros/list_replace.cmake)
include(CMake/Macros/get_target_name_numcse.cmake)

macro (add_executable_numcse _name)
    # cmake targets need to be unique
    # as such we prepend a prefix from the source path
    # all of this is done by get_executable_name_numcse
    get_target_name_numcse(${_name} target_name)
    #  reassemble argv for add_executable macro with new target name
    set(argv_parsed ${ARGV})
    LIST_REPLACE(argv_parsed 0 "${target_name}")

	# create directory that will contain the binary
	set( TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # retrieve relative path of runtime output directory (e.g the directory to which the target is compiled to)
	string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY})
	#	message(${TARGET_RUNTIME_OUTPUT_DIRECTORY})

    # invoke built-in add_executable
    add_executable(${argv_parsed})
    # link with mathgl and figure class
    if (TARGET ${target_name})
    	set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}"
        OUTPUT_NAME "${_name}")

        target_link_libraries(${target_name} Figure ${MATHGL2_LIBRARIES})
    	#add_dependencies(${target_name} Eigen)
    	#add_dependencies(${target_name} MathGL)
    endif()
endmacro()
