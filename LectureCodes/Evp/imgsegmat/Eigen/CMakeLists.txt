project(imgsegmat)
cmake_minimum_required(VERSION 2.8)

include_directories(../../imread/Eigen)

add_executable_numcse(main main.cpp)
add_resources(../../resources/segmentation)
