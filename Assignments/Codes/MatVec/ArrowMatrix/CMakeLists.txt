find_package(Boost REQUIRED)

include_directories("${Boost_INCLUDE_DIRS}")

add_definitions(-DSOLUTION=1)
add_definitions(-DINTERNAL=1)

add_executable_numcse(arrowmatvec arrowmatvec.cpp)
