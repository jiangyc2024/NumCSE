include_directories(../../../Utils)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/solution/decay.txt
        DESTINATION ${CMAKE_BINARY_DIR}/bin/Assignments/PolishedCodes/LeastSquares/Radioactive)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/Assignments/PolishedCodes/LeastSquares/Radioactive/cx_out)

add_executable_numcse(solution solution/main.cpp)
get_target_name_numcse(solution SOLUTION)
target_compile_definitions(${SOLUTION} PRIVATE SOLUTION=1 INTERNAL=1)

add_executable_numcse(template template/main.cpp)
get_target_name_numcse(template TEMPLATE)
target_compile_definitions(${TEMPLATE} PRIVATE SOLUTION=0 INTERNAL=1)

add_tests_numcse(tests template/tests.cpp)
