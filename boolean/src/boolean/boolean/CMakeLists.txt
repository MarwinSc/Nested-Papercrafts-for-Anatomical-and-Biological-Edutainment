# CMakeList.txt : CMake project for anatomicaledutainer, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project(${PROJECT_NAME} CXX)

find_package(CGAL)

add_library(${PROJECT_NAME} STATIC
"boolean.h"
"boolean.cpp")

target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(boolean_interface SHARED
"boolean.h"
"boolean.cpp")

target_link_libraries(boolean_interface CGAL::CGAL)
target_include_directories(boolean_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
