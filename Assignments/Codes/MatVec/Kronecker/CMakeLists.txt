add_executable_numcse(kron_solution kron.cpp)
get_target_name_numcse(kron_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(kron_template kron.cpp)
get_target_name_numcse(kron_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
