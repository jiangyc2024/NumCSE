project(prpowitsim)
cmake_minimum_required(VERSION 2.8)

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})
link_directories(${MKL_LIBRARIES})
target_link_libraries(<module>
	mkl_intel_lp64
	mkl_sequential
	mkl_core
)

include(../../../GetModules.cmake)
get_modules("Eigen Figure")
include_directories(${DIRS})
add_executable(main main.cpp)
add_custom_command(TARGET main PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/../../resources $<TARGET_FILE_DIR:main>)

target_link_libraries(main ${LIBS})
