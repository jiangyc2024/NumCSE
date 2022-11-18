project(numquad)

include_directories(../../../../Utils/)

add_executable_numcse(main main.cpp)

add_executable_numcse(integrate integrate.cpp)
