MACRO(set_eigen_fft_backend numcse_target_name backend)
	if (${backend} STREQUAL "FFTW")
		# find the fftw library if possible
		find_package(FFTW)

		if (FFTW_FOUND)
			# set include paths
			include_directories(${FFTW_INCLUDES})

			# link
			get_target_name_numcse(${numcse_target_name} target_name)
			target_link_libraries(${target_name} ${FFTW_LIBRARIES})

			# tell eigen to use fftw
			add_definitions(-DEIGEN_FFTW_DEFAULT)

			message("Eigen FFT backend set to FFTW for target `${target_name}`")
		else()
			message( WARNING "FFT backend `FFTW` not found. Fallback to default." )
		endif()
	elseif(${backend} STREQUAL "Kiss FFT")
		# kiss fft is already the default
    else()
        message( WARNING "Invalid FFT backend ${backend}. Fallback to default." )
    endif()
ENDMACRO()
