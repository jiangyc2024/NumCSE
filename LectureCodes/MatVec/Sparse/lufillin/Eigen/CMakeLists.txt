cmake_minimum_required(VERSION 2.8)
project(sparseformat)

include(../../../../GetModules.cmake)
get_modules("Eigen Figure")

include_directories(${DIRS}  
		../../../../../Utils/) # spdiags
add_executable(main main.cpp)
target_link_libraries(main ${LIBS})
