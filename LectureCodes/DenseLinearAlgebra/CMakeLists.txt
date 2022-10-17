set(USE_MKL ON)
find_package(MKL)

get_target_name_numcse(ncse_denselinalg target_name)

if(MKL_FOUND)
  message("${target_name} Compiling with MKL")
  # create executable
  add_executable_numcse(ncse_denselinalg main.cpp)
  
  # Specify libraries to link with
  target_link_libraries(${target_name} ${MKL_LIBRARIES})
  
  # Specify defines for MKL
  ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)
else()
  message("${target_name} MKL was not found.")
  add_executable_numcse(ncse_denselinalg main.cpp)
endif(MKL_FOUND)

add_executable_numcse(ncse_cost main_cost.cpp)
add_executable_numcse(ncse_cancel main_cancel.cpp)
