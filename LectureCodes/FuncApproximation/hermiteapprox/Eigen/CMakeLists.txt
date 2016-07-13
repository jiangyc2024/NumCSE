project(hermiteapprox)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../../Utils
                    ../../../Interpolation/hermloceval/Eigen)
add_executable_numcse(main main.cpp)
