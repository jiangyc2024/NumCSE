include_directories(../../../Utils)

add_executable_numcse(structured_matrix_solution multAmin.cpp)
get_target_name_numcse(structured_matrix_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(structured_matrix_template multAmin.cpp)
get_target_name_numcse(structured_matrix_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
