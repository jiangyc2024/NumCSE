#
# get a unique target (executable) name
#  the target name is build by
#  1. splitting the current source directory into chunks
#  2. where each chunk is converted 
#     from CamelCase to hyphen-seperated-lower-case
#  3. dropping Eigen suffix
#  4. afterwards all chunks are concatinated by underscores.
#
#  e.g. calling get_executable_name_numcse(main) in 
#  LectureCodes/NumQuad/numquaderrs/Eigen/CMakeLists.txt
#  will return lecture-codes_num-quad_traperr_main
macro(get_target_name_numcse base_name output_variable)
    #  convert to relative path
    string(REGEX REPLACE "${CMAKE_SOURCE_DIR}/+" "" target_prefix "${CMAKE_CURRENT_SOURCE_DIR}")
    #  stip eigen suffix
    string(REGEX REPLACE "/Eigen$" "" target_prefix "${target_prefix}")
    #  convert CamelCase to hyphen-seperated-lower-case
    string(REGEX REPLACE "([a-z])([A-Z])" "\\1-\\2" target_prefix "${target_prefix}")
    string(REGEX REPLACE "/" "_" target_prefix "${target_prefix}")
    string(TOLOWER "${target_prefix}" target_prefix)
    #  build target name from prefix and original target name 
    set(target_name "${target_prefix}_${base_name}")
    # return value
    set(${output_variable} "${target_name}")
endmacro()