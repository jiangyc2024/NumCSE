find_package(Boost REQUIRED)

include_directories("${Boost_INCLUDE_DIRS}")

add_executable_numcse(test test.cpp)
add_executable_numcse(timing timing.cpp)
