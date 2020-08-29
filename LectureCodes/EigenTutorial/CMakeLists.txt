# message(WARNING "NCSE Eigen Tutorial is skipped currently!")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Wno-deprecated-declarations -pedantic" )

add_executable_numcse (accesscolrow accesscolrow.cpp)
add_executable_numcse (blockaccess blockaccess.cpp)
add_executable_numcse (blockinit blockinit.cpp)
add_executable_numcse (eigentypedemo eigentypedemo.cpp)
add_executable_numcse (expressiontemplates expressiontemplates.cpp)
add_executable_numcse (initdense initdense.cpp)
add_executable_numcse (mapping mapping.cpp)
add_executable_numcse (matarray matarray.cpp)
add_executable_numcse (storageorder storageorder.cpp)
add_executable (ncse_eigenreductions eigenreductions.cpp)
add_executable (ncse_eigenoptiming eigenoptiming.cpp)

# This is the old code 
add_executable (ncse_eigentutorial main.cpp)
