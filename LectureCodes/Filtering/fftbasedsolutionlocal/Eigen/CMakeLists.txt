project(fftbasedsolutionlocal)
cmake_minimum_required(VERSION 2.8)

include_directories(../../sinetransform2d/Eigen/
                    )

add_executable_numcse(main main.cpp)
