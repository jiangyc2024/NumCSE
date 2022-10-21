set(USE_MKL ON)
find_package(MKL)

get_target_name_numcse(ncse_denselinalg target_name)

if(MKL_FOUND)
  message("${target_name} Compiling with MKL")

  add_executable_numcse(ncse_denselinalg main.cpp)

  set_target_properties(${target_name} PROPERTIES 
    
    COMPILE_DEFINITIONS "EIGEN_USE_MKL_ALL"
    COMPILE_OPTIONS "-fopenmp;-m64"
    LINK_OPTIONS "-fopenmp;-m64;-Wl,--no-as-needed"
    LINK_LIBRARIES "${MKL_LP_LIBRARY};${MKL_CORE_LIBRARY};${MKL_GNUTHREAD_LIBRARY};pthread"
  )

  target_include_directories(${target_name} PRIVATE ${MKL_INCLUDE_DIR})

else()
  message("${target_name} MKL was not found.")
  add_executable_numcse(ncse_denselinalg main.cpp)
endif(MKL_FOUND)

add_executable_numcse(ncse_cost main_cost.cpp)
add_executable_numcse(ncse_cancel main_cancel.cpp)
