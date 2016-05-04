project(numquad)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../FunctionCollection/
                    ../../gaussquad/Eigen/)
add_executable_numcse(main main.cpp)
