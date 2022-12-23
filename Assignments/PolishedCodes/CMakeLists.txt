add_all_subdirectories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable_numcse(polyevalnewton polyevalnewton.cpp)

add_executable_numcse(fftmisc fftmisc.cpp)

add_executable_numcse(quadformulas quadformulas.cpp)
