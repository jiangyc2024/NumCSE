find_package(Boost REQUIRED)

include_directories("${Boost_INCLUDE_DIRS}")

add_definitions(-DSOLUTION=1)
add_definitions(-DINTERNAL=1)

add_executable_numcse(gram_schmidt_solution gramschmidt.cpp)

add_definitions(-DSOLUTION=0)
add_definitions(-DINTERNAL=0)

add_executable_numcse(gram_schmidt_template gramschmidt.cpp)
