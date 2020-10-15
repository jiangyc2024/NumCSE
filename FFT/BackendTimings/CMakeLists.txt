# add executable
add_executable_numcse(main fftbackendtime.cpp)
set_eigen_fft_backend(main "Kiss FFT")

