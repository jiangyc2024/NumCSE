
# SET LINEAR ALGEBRA SETTINGS

set( LINALG_DEFINE "_MATRIX" )
#set( LINALG_DEFINE "EIGEN3" )


#-------------------------------------------------------------------------------
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D ${LINALG_DEFINE}")
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../../Modules" )

if ( LINALG_DEFINE STREQUAL "EIGEN3" )
  # find linear algebra package
  find_package( Eigen3 )
  include_directories( ${EIGEN3_INCLUDE_DIR} )
  message( STATUS "Eigen3 includes: ${EIGEN3_INCLUDE_DIR}" )
endif ()

message(STATUS "Linear algebra settings: ${LINALG_DEFINE}")



