project(seg)

include_directories(../../imread/Eigen)
include_directories(../../imgsegmat/Eigen)
include_directories(../../../../Utils)

add_executable_numcse(main main.cpp)
add_resources(../../resources/segmentation)
