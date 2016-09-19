include_directories(../../../Utils)

add_executable_numcse(arrow_mat_vec_solution arrowmatvec.cpp)
get_target_name_numcse(arrow_mat_vec_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(arrow_mat_vec_template arrowmatvec.cpp)
get_target_name_numcse(arrow_mat_vec_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)
