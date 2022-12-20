//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY__USERINPUTS_H_
#define GEOMETRY__USERINPUTS_H_

#include "imgui.h"
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>
#include "GeometryMaster.h"
#include "Inspect/InspectNavigateTool.h"

class UserInputs {
 public:
  UserInputs(ImGuiIO& io, GLFWwindow* window) : io_(io), window_(window) {}

  bool HasMouseInput() const { return !io_.WantCaptureMouse; }

  bool HasKeyboardInput() const { return !io_.WantCaptureKeyboard; }

  bool HasRefreshInput() const { return ImGui::IsKeyPressed(526); }

  bool HasScaleInput() const { return io_.MouseWheel != 0.f; }

  bool HasClicked() const { return ImGui::IsMouseClicked(0); }

  bool IsDown() const { return ImGui::IsMouseDown(0); }

  bool HasReleased() const { return ImGui::IsMouseReleased(0); }

  bool HasChangeTool() const { return user_type_ != ToolType::None; }

  ToolType GetToolType() const { return user_type_; }

  float GetScale() const { return 1.f + io_.MouseWheel * 0.1f; }

  Vector2f GetMousePos() const {
    Vector2f relative_mouse_position = GetCursorPosition();
    Vector2f window_size = GetSize();
    auto frame_mouse_pos = relative_mouse_position.InvertedScaled(window_size) * 2.f;
    frame_mouse_pos[1] *= -1.f;
    frame_mouse_pos += Vector2f(-1.f, 1.f);
    return frame_mouse_pos;
  }

  void SetToolType(ToolType type) { user_type_ = type; }

 private:
  Vector2d GetCursorPosition() const {
    Vector2d position;
    glfwGetCursorPos(window_, &position[0], &position[1]);
    return position;
  }

  Vector2i GetSize() const {
    Vector2i size;
    glfwGetWindowSize(window_, &size[0], &size[1]);
    return size;
  }

  ImGuiIO& io_;
  GLFWwindow* window_;
  bool change_tool_;
  ToolType user_type_ = ToolType::None;
};

#endif //GEOMETRY__USERINPUTS_H_
