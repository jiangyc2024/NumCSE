subdirlist(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

foreach(SUBDIR ${SUBDIRS})
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SUBDIR}/CMakeLists.txt")
		add_subdirectory(${SUBDIR})
	else()
		message(WARNING "Skipping ${SUBDIR} because no CMakeLists.txt file was found")
	endif()
endforeach()