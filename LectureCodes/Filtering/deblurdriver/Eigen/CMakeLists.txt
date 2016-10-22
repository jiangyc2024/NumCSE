project(deblurdriver)
cmake_minimum_required(VERSION 2.8)

include_directories(../../blur/Eigen/
                    ../../deblur/Eigen/
                    ../../psf/Eigen/
                    )

add_executable_numcse(main deblurdriver.cpp)
