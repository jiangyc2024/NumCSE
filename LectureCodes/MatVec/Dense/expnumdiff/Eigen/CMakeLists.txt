cmake_minimum_required(VERSION 2.8)
project(expnumdiff)

# create executable
add_executable_numcse(expnumdiff_main main.cpp)

# link with mpfr gmp
get_target_name_numcse(expnumdiff_main target_name)
target_link_libraries(${target_name} ${LIBS} mpfr gmp)
