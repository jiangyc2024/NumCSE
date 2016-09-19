include_directories(../../../Utils)

add_executable_numcse(mat_pow_solution matPow.cpp)
get_target_name_numcse(mat_pow_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(mat_pow_template matPow.cpp)
get_target_name_numcse(mat_pow_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
