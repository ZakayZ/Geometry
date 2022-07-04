cmake_minimum_required(VERSION 3.21)
project(Geomerty)

set(CMAKE_CXX_STANDARD 20)

add_executable(Geomerty main.cpp)

target_include_directories(Geomerty PUBLIC ${PROJECT_SOURCE_DIR})