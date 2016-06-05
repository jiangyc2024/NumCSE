project(numquaderrs)
cmake_minimum_required(VERSION 2.8)

include_directories(../../../../Utils/
                    ../../gaussquad/Eigen/
                    ../../numquad/Eigen/ )
add_executable_numcse(main numquaderrs.cpp)
