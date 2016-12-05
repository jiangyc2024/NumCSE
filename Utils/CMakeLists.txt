add_subdirectory(examples)

file(GLOB_RECURSE HEADERS *.hpp)

add_custom_target(utils SOURCES ${HEADERS} README.md)
