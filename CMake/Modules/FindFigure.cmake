# ================================= FindFigure.cmake =======================================
# once finshed the following variables will be initialized:
#   FIGURE_INCLUDE_DIR : directory which contains the Figure header ( figure.hpp )
#   FIGURE_LIBRARY     : Figure library ( libFigure.a )
#   FIGURE_LIBRARY_DIR : directory which contains the Figure library ( libFigure.a )
#   FIGURE_NOT_FOUND   : if Figure wasn't found this variable contains further information, if it was found it's not defined
# ==========================================================================================================================
#   Typical usage:    find_package( Figure REQUIRED)
#                     include_directories( ${FIGURE_INCLUDE_DIR} )
#                     add_executable( main my_main_file.cpp )
#                     target_link_libraries( main ${FIGURE_LIBRARY} )
# ==========================================================================================================================

# get path to figure header
find_path(FIGURE_INCLUDE_DIR NAMES figure.hpp DOC "Figure header")
if (FIGURE_INCLUDE_DIR)
  message(STATUS "Found Figure header in: ${FIGURE_INCLUDE_DIR}")
else ()
  message(STATUS "Could not find header of Figure ( figure.hpp ). Did you try installing it with sudo make install?")
endif()

# find libFigure.a
find_library(FIGURE_LIBRARY NAMES libFigure.a DOC "Figure library")
if (FIGURE_LIBRARY)
  get_filename_component(FIGURE_LIBRARY_DIR ${FIGURE_LIBRARY} PATH)
  message(STATUS "Found Figure library in: ${FIGURE_LIBRARY_DIR}")
else()
  message(STATUS "Could not find library of Figure ( libFigure.a ). Did you try installing it with sudo make install?")
endif()

# set FIGURE_NOT_FOUND if necessary
if (NOT FIGURE_INCLUDE_DIR AND FIGURE_LIBRARY)
  set(FIGURE_NOT_FOUND "Found Figure library but couldn't find according header - try reinstalling!")
elseif(FIGURE_INCLUDE_DIR AND NOT FIGURE_LIBRARY)
  set(FIGURE_NOT_FOUND "Found Figure header but couldn't find according library - try reinstalling!")
elseif (NOT FIGURE_INCLUDE_DIR AND NOT FIGURE_LIBRARY)
  set(FIGURE_NOT_FOUND "Couldn't find either Figure library not header - try reinstalling!")
endif()

# print error message if necessary
if (FIGURE_NOT_FOUND)
  if(Figure_FIND_REQUIRED) # if Figure was marked as REQUIRED and wasn't found throw an FATAL_ERROR
    message(FATAL_ERROR "${FIGURE_NOT_FOUND}")
  else() # otherwise just print a STATUS message
    message(STATUS "${FIGURE_NOT_FOUND}")
  endif()
endif()
