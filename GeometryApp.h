//
// Created by Artem Novikov on 07.08.2022.
//

#ifndef GEOMETRY__GEOMETRYAPP_H_
#define GEOMETRY__GEOMETRYAPP_H_

#include "Window.h"
#include "GeometryMaster.h"
#include "UI.h"

void WindowResizeCallback(GLFWwindow*, int width, int height);

void WindowMoveCallback(GLFWwindow*, int new_x, int new_y);

class GeometryApp {
 public:
  ~GeometryApp() = default;

  void RunApp() {
    while (window_.IsOpen()) {
      StartFrame();
      RenderFrame();
      EndFrame();
    }
  }

 private:
  friend class AppManager;
  friend void WindowResizeCallback(GLFWwindow*, int, int);
  friend void WindowMoveCallback(GLFWwindow*, int, int);

  GeometryApp(int window_width, int widow_height)
      : window_(window_width, widow_height), ui_(window_.ShareWindow()),
        geometry_master_(window_width, widow_height), user_inputs_(ui_.GetIO()) {
    window_.SetMoveCallback(WindowMoveCallback);
    window_.SetResizeCallback(WindowResizeCallback);
  }

  void StartFrame() {
    window_.StartFrame();
    ui_.StartFrame();
    ProcessInputs();
  }

  void ProcessInputs() { /// TODO refactor
    if (!user_inputs_.WantCaptureMouse) {
      /// Mouse position
      Vector2f mouse_position = {user_inputs_.MousePos.x, user_inputs_.MousePos.y};
      Vector2f window_position = window_.GetPosition();
      Vector2f window_size = window_.GetSize();
      auto relative_mouse_pos = mouse_position - window_position;
      auto frame_mouse_pos = relative_mouse_pos.InvertedScaled(window_size) * 2.f - Vector2f(1.f, 1.f);
      frame_mouse_pos = geometry_master_.ProcessHover(frame_mouse_pos);
      if (ImGui::IsMouseClicked(0)) {
        geometry_master_.ProcessPressed(frame_mouse_pos);
      }
      if (ImGui::IsMouseDown(0)) {
        geometry_master_.ProcessDown(frame_mouse_pos);
      }
      if (ImGui::IsMouseReleased(0)) {
        geometry_master_.ProcessReleased(frame_mouse_pos);
      }
    }

    if (!user_inputs_.WantCaptureKeyboard) {
      /// ESC button
      if (ImGui::IsKeyPressed(526)) {
        geometry_master_.ProcessRefresh();
      }
    }
  }

  void RenderFrame() {
    window_.Clear();
    ui_.Render();
    geometry_master_.Render();
  }

  void EndFrame() {
    ui_.EndFrame();
    window_.EndFrame();
  }

  void WindowResized(const Vector2f& new_size) {
    geometry_master_.ProcessWindowResize(new_size);
  }

  void WindowMoved(const Vector2f& new_position) {
    Vector2f old_position = window_.GetPosition();
    geometry_master_.ProcessWindowShift(new_position - old_position);
  }

  Window window_;
  UI ui_;
  GeometryMaster geometry_master_;
  ImGuiIO& user_inputs_;
};

class AppManager {
 public:
  AppManager() = default;

  AppManager(int window_width, int widow_height) {
    if (app_ == nullptr) {
      app_ = new GeometryApp(window_width, widow_height);
    }
  }

  GeometryApp& GetApp() {
    return *app_;
  }

  void DestroyApp() {
    delete app_;
  }

 private:
  static GeometryApp* app_;
};

void WindowResizeCallback(GLFWwindow*, int width, int height) {
  AppManager app_manager;
  auto& app = app_manager.GetApp();
  app.WindowResized(Vector2f(width, height));
}

void WindowMoveCallback(GLFWwindow*, int new_x, int new_y) {
  AppManager app_manager;
  auto& app = app_manager.GetApp();
  app.WindowMoved(Vector2f(new_x, new_y));
}

#endif //GEOMETRY__GEOMETRYAPP_H_
