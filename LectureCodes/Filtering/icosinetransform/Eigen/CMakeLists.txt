project(icosinetransform)
cmake_minimum_required(VERSION 2.8)

include_directories(../../cosinetransform/Eigen/
                    )

add_executable_numcse(main main.cpp)
set_eigen_fft_backend(main "FFTW")
