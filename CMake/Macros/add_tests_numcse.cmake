include(CMake/Macros/get_target_name_numcse.cmake)

# similiar to add_executable_numcse.cmake

macro (add_tests_numcse _name)
  get_target_name_numcse(copy target_name_copy)
  get_target_name_numcse(${_name} target_name)

  set(argv_parsed ${ARGV})
  LIST_REPLACE(argv_parsed 0 "${target_name}")

  set(TARGET_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  string(REPLACE "${CMAKE_SOURCE_DIR}" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" TARGET_RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}")

  add_custom_target(${target_name_copy} ALL COMMAND python3 ${PROJECT_SOURCE_DIR}/Testing/copy_and_tweak.py ${CMAKE_CURRENT_SOURCE_DIR} ${TARGET_RUNTIME_OUTPUT_DIRECTORY})
  add_executable(${argv_parsed})
  add_dependencies(${target_name} ${target_name_copy})
  set_target_properties(${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TARGET_RUNTIME_OUTPUT_DIRECTORY}" OUTPUT_NAME "${_name}")
  add_test(NAME ${target_name} COMMAND ${target_name})

  # link with python3
  target_link_libraries(${target_name} ${Python3_LIBRARIES})
  
  # link with copy.hpp
  target_include_directories(${target_name} PRIVATE ${TARGET_RUNTIME_OUTPUT_DIRECTORY})

endmacro()
