# message(WARNING "NCSE Eigen Tutorial is skipped currently!")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Wno-deprecated-declarations -pedantic" )

add_executable (ncse_eigentutorial main.cpp)
add_executable (ncse_eigenreductions eigenreductions.cpp)

