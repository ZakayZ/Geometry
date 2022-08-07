cmake_minimum_required(VERSION 3.21)
project(Geometry)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

set(CMAKE_CXX_STANDARD 17)

add_executable(Geometry app.cpp)

target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(Geometry PUBLIC ${PROJECT_SOURCE_DIR}/Tools)

###################################################### SFML

find_package(SFML 2.5.1 COMPONENTS graphics window REQUIRED)

target_link_libraries(Geometry sfml-graphics sfml-window)

###################################################### GLFW

find_package(glfw3 3.3 REQUIRED)

target_link_libraries(Geometry glfw)

###################################################### ImGUI

add_subdirectory(${PROJECT_SOURCE_DIR}/ImGui)

target_link_libraries(Geometry ImGui)
#target_link_libraries(Geometry ImGuiSFML)