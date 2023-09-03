MACRO(ADD_ALL_SUBDIRECTORIES directory)
  subdirlist(SUBDIRS ${directory})

  foreach(SUBDIR ${SUBDIRS})
    if(EXISTS "${directory}/${SUBDIR}/CMakeLists.txt")
      add_subdirectory(${SUBDIR})
    else()
		# check whether the subdirectory contains cpp files
		execute_process(COMMAND bash -c "find ${directory}/${SUBDIR} | grep .cpp"
						  RESULT_VARIABLE contains_cpp_code
						  OUTPUT_QUIET)

		# if the directory contains cpp files we print a warning
		if(${contains_cpp_code} EQUAL 0)
			message("Skipping ${directory}/${SUBDIR} (containing cpp code) because no CMakeLists.txt file was found")
		else() # otherwise just a debug message
			message("Skipping ${directory}/${SUBDIR} since it did not contain any cpp code (e.g. not ported yet)")
		endif()		
    endif()
  endforeach()
ENDMACRO()
