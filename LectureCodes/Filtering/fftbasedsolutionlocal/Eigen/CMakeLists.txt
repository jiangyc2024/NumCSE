project(fftbasedsolutionlocal)

include_directories(../../sinetransform2d/Eigen/
                    )

add_executable_numcse(main main.cpp)
