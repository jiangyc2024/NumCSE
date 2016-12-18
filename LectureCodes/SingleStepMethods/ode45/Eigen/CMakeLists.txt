project(ss_ode45)
cmake_minimum_required(VERSION 2.8)

add_executable_numcse(main main.cpp)
add_executable_numcse(ode45test ode45test.cpp)
add_executable_numcse(ode45stiff ode45stiff.cpp)
