include_directories(../../../Utils)

add_executable_numcse(solution solution/main.cpp)
get_target_name_numcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(template template/main.cpp)
get_target_name_numcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)

add_tests_numcse(tests template/tests.cpp)
file(MAKE_DIRECTORY ${TARGET_RUNTIME_OUTPUT_DIRECTORY}/cx_out)

