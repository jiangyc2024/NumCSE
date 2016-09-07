cmake_minimum_required(VERSION 2.8)
project(expnumdiff)

find_package(MPFR)
find_package(GMP)
find_package(MPFR_CXX)

get_target_name_numcse(main target_name)

if(TRUE  )
	include_directories(${GMP_INCLUDES} ${MPFR_INCLUDES})

	# create executable
	add_executable_numcse(main main.cpp)

	# link with mpfr gm
	target_link_libraries(${target_name} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
else()
	message("Skipped ${target_name} because MPFR was not found.")
endif()
