cmake_minimum_required (VERSION 3.8)

project(mim C CXX)

find_package(CGAL)
include_directories(".")

add_library(mim SHARED "mesh_inside.cpp")

target_link_libraries(mim CGAL::CGAL)
target_include_directories(mim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})