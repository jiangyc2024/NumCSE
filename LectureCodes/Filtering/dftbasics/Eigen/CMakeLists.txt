cmake_minimum_required(VERSION 2.8)
project(dftbasics)

add_executable_numcse(main main.cpp)
set_eigen_fft_backend(main "FFTW")
