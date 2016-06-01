project(eigtiming)
cmake_minimum_required(VERSION 2.8)

# SET(GCC_MKL_LINKER_FLAGS    "-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl")
# SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_MKL_LINKER_FLAGS}" )

include(../../../GetModules.cmake)
get_modules("Eigen Figure")
include_directories(${DIRS})
add_executable(main main.cpp)
target_link_libraries(main ${LIBS})
