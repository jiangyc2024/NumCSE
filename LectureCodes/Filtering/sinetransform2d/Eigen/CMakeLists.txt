project(sinetransform2d)
cmake_minimum_required(VERSION 3.0.2)

add_executable_numcse(main main.cpp)
set_eigen_fft_backend(main "Kiss FFT")
