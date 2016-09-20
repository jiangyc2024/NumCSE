
add_executable_numcse(test_solution test.cpp)
get_target_name_numcse(test_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(test_template test.cpp)
get_target_name_numcse(test_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)

add_executable_numcse(timing_solution timing.cpp)
get_target_name_numcse(timing_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(timing_template timing.cpp)
get_target_name_numcse(timing_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)
