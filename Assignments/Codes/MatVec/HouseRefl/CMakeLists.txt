include_directories(../../../Utils)

add_executable_numcse(houserefl_solution houserefl.cpp)
get_target_name_numcse(houserefl_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(houserefl_template houserefl.cpp)
get_target_name_numcse(houserefl_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
