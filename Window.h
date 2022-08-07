//
// Created by Artem Novikov on 07.08.2022.
//

#ifndef GEOMETRY__WINDOW_H_
#define GEOMETRY__WINDOW_H_

#include <cstdio>
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include <array>

#include "GeometryEntities/Vector.h"

static void glfw_error_callback(int error, const char* description) {
  fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

class Window {
 public:
  Window(int window_width, int window_height) {
    glfwSetErrorCallback(glfw_error_callback);
    assert(glfwInit());
    window_ = glfwCreateWindow(window_width, window_height, "Dear ImGui GLFW+OpenGL2 example", nullptr, nullptr);
    assert(window_ != nullptr);
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // Enable vsync
  }

  ~Window() {
    glfwDestroyWindow(window_);
    glfwTerminate();
  }

  auto ShareWindow() { return window_; }

  Vector2i GetPosition() {
    Vector2i position;
    glfwGetWindowPos(window_, &position[0], &position[1]);
    return position;
  }

  Vector2i GetSize() {
    Vector2i size;
    glfwGetWindowPos(window_, &size[0], &size[1]);
    return size;
  }

  bool IsOpen() const {
    return !glfwWindowShouldClose(window_);
  }

  void StartFrame() {
    glfwPollEvents();
  }

  void Clear() {
    int display_w, display_h;
    glfwGetFramebufferSize(window_, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(clear_color[0] * clear_color[3],
                 clear_color[1] * clear_color[3],
                 clear_color[2] * clear_color[3],
                 clear_color[3]);
    glClear(GL_COLOR_BUFFER_BIT);
  }

  void EndFrame() {
    glfwSwapBuffers(window_);
  }

  template <class Func>
  void SetResizeCallback(Func function) {
    glfwSetWindowSizeCallback(window_, function); /// TODO if doesnt work use Frame buffer size
  }

  template <class Func>
  void SetMoveCallback(Func function) {
    glfwSetWindowPosCallback(window_, function);
  }

 private:
  GLFWwindow* window_;
  std::array<float, 4> clear_color = {0.45f, 0.55f, 0.60f, 1.00f};
};

#endif //GEOMETRY__WINDOW_H_
