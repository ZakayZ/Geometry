cmake_minimum_required(VERSION 3.21)
project(Geomerty)

set(CMAKE_CXX_STANDARD 20)

find_package(SFML 2.5.1 COMPONENTS graphics REQUIRED)

add_executable(Geomerty main.cpp)

target_include_directories(Geomerty PUBLIC ${PROJECT_SOURCE_DIR})

target_link_libraries(Geomerty sfml-graphics)