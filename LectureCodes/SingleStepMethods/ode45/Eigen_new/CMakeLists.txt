project(Ode45)
cmake_minimum_required(VERSION 2.8)

# define executables
add_executable_numcse(ode45test_new ode45test.cc ode45.h)
add_executable_numcse(ode45stiff_new ode45stiff.cc ode45.h)
add_executable_numcse(odeintssctrltest_new odeintssctrltest.cc odeintssctrl.h)
add_executable_numcse(embeddedrkssm_new embeddedrkssm.cc)
