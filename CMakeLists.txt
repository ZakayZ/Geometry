cmake_minimum_required(VERSION 3.21)
project(Geometry)

set(CMAKE_CXX_STANDARD 20)

find_package(SFML 2.5.1 COMPONENTS graphics REQUIRED)

add_executable(Geometry main.cpp)

target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR})

target_link_libraries(Geometry sfml-graphics)