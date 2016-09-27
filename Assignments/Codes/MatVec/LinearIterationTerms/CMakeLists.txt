include_directories(../../../Utils)

add_executable_numcse(solution getit.cpp)
get_target_name_numcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)
