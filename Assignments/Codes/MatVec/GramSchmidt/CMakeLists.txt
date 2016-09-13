find_package(Boost REQUIRED)

include_directories("${Boost_INCLUDE_DIRS}")

add_executable_numcse(gram_schmidt_solution gramschmidt.cpp)
get_target_name_numcse(gram_schmidt_solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(gram_schmidt_template gramschmidt.cpp)
get_target_name_numcse(gram_schmidt_template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=0)
