include_directories(../../../Utils)

#copy files from source dir to binary dir
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/solution/image.pgm
        DESTINATION ${CMAKE_BINARY_DIR}/bin/Assignments/Codes/Filtering/Autofocus/CodeExpert)
#make output directory
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/Assignments/Codes/Filtering/Autofocus/CodeExpert/cx_out)

add_executable_numcse(solution solution/main.cpp)
get_target_name_numcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)
#use O3 flag exclusively for this task
set_target_properties(${SOLUTION} PROPERTIES COMPILE_FLAGS "-O3")

add_executable_numcse(template template/main.cpp)
get_target_name_numcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)
#use O3 flag exclusively for this task
set_target_properties(${TEMPLATE} PROPERTIES COMPILE_FLAGS "-O3")

add_tests_numcse(tests template/tests.cpp)
get_target_name_numcse(tests TESTS)
#use O3 flag exclusively for this task
set_target_properties(${TESTS} PROPERTIES COMPILE_FLAGS "-O3")
