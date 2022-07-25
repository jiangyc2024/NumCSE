project(imgsegmat)

include_directories(../../imread/Eigen)

add_executable_numcse(main main.cpp)
add_resources(../../resources/segmentation)
