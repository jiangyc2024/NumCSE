include_directories(../../../Utils)

add_executable_numcse(solution toeplitz.cpp)
get_target_name_numcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(template toeplitz.cpp)
get_target_name_numcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)

add_custom_target(${SOLUTION}-extra_files SOURCES "description.json")
