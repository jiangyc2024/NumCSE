# ================================= FindFigure.cmake =======================================
# once finshed the following variables will be initialized:
#   FIGURE_INCLUDE_DIR : directory which contains all Figure files
#   FIGURE_LIBRARY : libFigure.a
#   FIGURE_FOUND : true, if all parts of the library have been found. Otherwise false.
#   FIGURE_FILE_LIST : list of all needed files 
# ==========================================================================================
#   Typical usage:    find_package( Figure REQUIRED)
#                     include_directories( ${FIGURE_INCLUDE_DIR} )
#                     add_executable( main my_main_file.cpp )
#                     target_link_libraries( main Figure )
# ==========================================================================================

if ( DEBUG )
  message( STATUS "Running FindFigure.cmake with DEBUG option .." )
endif()

# setting variable FIGURE_FILE_LIST
set( FIGURE_FILE_LIST figure.hpp
                      figure.cpp
                      FigureConfig.hpp
                      MglLabel.hpp
                      MglPlot.hpp
                      MglStyle.hpp
                      )


## ------------ find header files ------------- ##
# get paths to figure files
set( FIGURE_PATH_SUFFIX figure )
find_path( FIGURE_HPP NAMES figure.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "Figure header" )
find_path( FIGURE_CPP NAMES figure.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "Figure source" )
find_path( FIGURE NAMES Figure PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "Figure include" )
find_path( FIGURECONFIG_HPP NAMES FigureConfig.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "Figure config" )

find_path( MGL_LABEL_HPP NAMES MglLabel.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "MglLabel" )
find_path( MGL_PLOT_HPP NAMES MglPlot.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "MglPlot" )
find_path( MGL_STYLE_HPP NAMES MglStyle.hpp PATH_SUFFIXES ${FIGURE_PATH_SUFFIX} DOC "MglStyle" )

set( FIGURE_PATHS ${FIGURE_HPP} 
                  ${FIGURECONFIG_HPP}
                  ${MGL_LABEL_HPP}
                  ${MGL_PLOT_HPP}
                  ${MGL_STYLE_HPP}
                  )

if ( DEBUG )
  message( STATUS "Needed files are: figure.hpp FigureConfig.hpp MglLabel.hpp MglPlot.hpp MglStyle.hpp" )
endif()

# check if the files are all in the correct place
set( PATHS_ARE_VALID true )
foreach( FILE_PATH ${FIGURE_PATHS} )

  if( NOT ${FILE_PATH} STREQUAL ${FIGURE_HPP} ) 
    set( PATHS_ARE_VALID false )

    if ( DEBUG ) # this will print which file exactly wasnt found
      message( STATUS "${FILE_PATH}" )
    endif()
  
  endif()

endforeach()

# if all filepaths are valid set FIGURE_INCLUDE_DIR, otherwise throw error/warning
if( PATHS_ARE_VALID )
  
  message( STATUS "Found Figure in ${FIGURE_HPP}" )
  set( FIGURE_HEADERS_FOUND true )
  set( FIGURE_INCLUDE_DIR ${FIGURE_HPP} )

else()
  
  set( FIGURE_HEADERS_FOUND false )
  if( Figure_REQUIRED )
    message( FATAL_ERROR "Couldn't find all necessary header files for Figure, maybe try (re-)installing with administrator rights?" )
    message( "   Run cmake with the '-DDEBUG=1' option for more information on which files are missing." )
  else()
    message( STATUS "FindFigure.cmake couldn't find all necessary header files for Figure. Not stopping as Figure is not marked REQUIRED" )
  endif()

endif()

## ------------ find libFigure.a -------------- ##
find_library( FIGURE_LIBRARY NAMES libFigure.a DOC "Figure library" )
if ( FIGURE_LIBRARY )
  get_filename_component( FIGURE_LIBRARY_DIR ${FIGURE_LIBRARY} PATH )
  message( STATUS "Found Figure library: ${FIGURE_LIBRARY_DIR}" )
else()
  if ( Figure_REQUIRED )
    message( FATAL_ERROR "Couldn't find libFigure.a, maybe try (re-)installing with administrator rights?" )
  else()
    message( STATUS "Couldn't find libFigure.a, maybe try (re-)installing with administrator rights?" )
  endif()
endif()

# if the library and the header files were found set FIGURE_FOUND to true
if ( FIGURE_LIBRARY AND FIGURE_HEADERS_FOUND )
  set( FIGURE_FOUND true )
else()
  set( FIGURE_FOUND false )
endif()
