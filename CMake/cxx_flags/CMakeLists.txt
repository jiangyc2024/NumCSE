###########################################################################
#
# SET COMPILER FLAGS
#
###########################################################################

###########################################################################
#
# OS DEPENDENT SECTION 
#
###########################################################################
if (UNIX AND NOT APPLE) 
  include ("${CMAKE_CURRENT_LIST_DIR}/unix/CMakeLists.txt")
endif (UNIX AND NOT APPLE)

if (APPLE)
  include ("${CMAKE_CURRENT_LIST_DIR}/apple/CMakeLists.txt")
endif (APPLE)

if (WIN32)
  include ("${CMAKE_CURRENT_LIST_DIR}/win/CMakeLists.txt")
endif (WIN32)

message( STATUS "CXX flags set to ${CMAKE_CXX_FLAGS}" )

###########################################################################

