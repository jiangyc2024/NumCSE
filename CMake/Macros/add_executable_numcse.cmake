include(CMake/Macros/list_replace.cmake)

macro (add_executable_numcse _name)
    # cmake targets need to be unique
    # as such we prepend a prefix from the source path
    #  convert to relative path
    string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/+" "" target_prefix "${CMAKE_CURRENT_SOURCE_DIR}")
    #  stip eigen suffix
    string(REGEX REPLACE "/Eigen$" "" target_prefix "${target_prefix}")
    #  convert CamelCase to hyphen-seperated-lower-case
    string(REGEX REPLACE "([a-z])([A-Z])" "\\1-\\2" target_prefix "${target_prefix}")
    string(REGEX REPLACE "/" "_" target_prefix "${target_prefix}")
    string(TOLOWER "${target_prefix}" target_prefix)
    #  build target name from prefix and original target name 
    set(target_name "${target_prefix}_${_name}")
    #  reassemble argv for add_executable macro with new target name
    set(argv_parsed ${ARGV})
    LIST_REPLACE(argv_parsed 0 "${target_name}")

	# create directory that will contain the binary
	set( TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # retrieve relative path of runtime output directory (e.g the directory to which the target is compiled to)
	string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY})

    # invoke built-in add_executable
    add_executable(${argv_parsed})
    # link with mathgl and figure class
    if (TARGET ${target_name})
    	set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}"
        OUTPUT_NAME "${_name}")

        target_link_libraries(${target_name} ${MATHGL2_LIBRARIES} Figure)
    	#add_dependencies(${target_name} Eigen)
    	#add_dependencies(${target_name} MathGL)
    endif()
endmacro()