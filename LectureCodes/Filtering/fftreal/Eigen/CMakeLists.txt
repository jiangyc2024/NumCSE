project(fftreal)

add_executable_numcse(main main.cpp)
set_eigen_fft_backend(main "Kiss FFT")
