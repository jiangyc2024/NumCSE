# add executable
add_executable_numcse(main ffttime.cpp)
set_eigen_fft_backend(main "Kiss FFT")

