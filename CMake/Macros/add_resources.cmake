
macro (add_resources _folder)

	add_custom_command(TARGET ${target_name} 
		PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory 
		${CMAKE_CURRENT_SOURCE_DIR}/${_folder} $<TARGET_FILE_DIR:${target_name}>)

endmacro()
