#################################################### ImGui + GLFW

set(SOURCES
        imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp
        Imp/imgui_impl_glfw.cpp Imp/imgui_impl_opengl2.cpp)

add_library(ImGui SHARED ${SOURCES})


### OpenGL
find_package(OpenGL REQUIRED)
target_link_libraries(ImGui OpenGL::GL)

### GLFW
find_package(glfw3 REQUIRED)
target_link_libraries(ImGui glfw)

target_include_directories(ImGui PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Imp")
target_include_directories(ImGui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


##################################################### ImGui + SFML

set(SFMLGUI_SOURCES
        imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_tables.cpp imgui_widgets.cpp
        Imp/imgui-SFML.cpp)

add_library(ImGuiSFML SHARED ${SFMLGUI_SOURCES})

### SFML
find_package(SFML 2.5.1 COMPONENTS graphics window REQUIRED)
target_link_libraries(ImGuiSFML sfml-graphics sfml-window)

### OpenGL
find_package(OpenGL REQUIRED)
target_link_libraries(ImGuiSFML OpenGL::GL)

target_include_directories(ImGuiSFML PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Imp")
target_include_directories(ImGuiSFML PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})