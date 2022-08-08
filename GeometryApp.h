//
// Created by Artem Novikov on 07.08.2022.
//

#ifndef GEOMETRY__GEOMETRYAPP_H_
#define GEOMETRY__GEOMETRYAPP_H_

#include "Window.h"
#include "GeometryMaster.h"
#include "UI.h"

class GeometryApp {
 public:
  GeometryApp(int window_width, int widow_height)
      : window_(window_width, widow_height), ui_(window_.ShareWindow()),
        geometry_master_(window_width, widow_height), user_inputs_(ui_.GetIO()) {}

  ~GeometryApp() = default;

  void RunApp() {
    while (window_.IsOpen()) {
      StartFrame();
      RenderFrame();
      EndFrame();
    }
  }

 private:
  void StartFrame() {
    window_.StartFrame();
    ui_.StartFrame();
    ProcessInputs();
  }

  void ProcessInputs() {
    if (window_.IsResized()) {
      OnWindowResize();
    }

    if (window_.IsMoved()) {
      OnWindowMove();
    }

    if (window_.IsFocused()) {  /// TODO refactor
      /// Mouse position
      if (!user_inputs_.WantCaptureMouse) {
        Vector2f relative_mouse_position = window_.GetCursorPosition();
        Vector2f window_size = window_.GetSize();
        auto frame_mouse_pos = relative_mouse_position.InvertedScaled(window_size) * 2.f;
        frame_mouse_pos[1] *= -1.f;
        frame_mouse_pos += Vector2f(-1.f, 1.f);
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


        /// Mouse wheel
        if (user_inputs_.MouseWheel != 0) {
          geometry_master_.ProcessWindowScale(user_inputs_.MouseWheel);
        }
      }

      /// ESC button
      if (!user_inputs_.WantCaptureKeyboard) {
        if (ImGui::IsKeyPressed(526)) {
          geometry_master_.ProcessRefresh();
        }
      }
    }
  }

  void RenderFrame() {
    window_.Clear();
    geometry_master_.Render();
    ui_.Render();
  }

  void EndFrame() {
    ui_.EndFrame();
    window_.EndFrame();
  }

  void OnWindowResize() {
    geometry_master_.ProcessWindowResize(window_.GetSize());
    window_.UpdateSize();
  }

  void OnWindowMove() {
    Vector2f shift = window_.GetPosition() - window_.GetOldPosition();
    shift[1] *= -1;
    geometry_master_.ProcessWindowShift(shift);
    window_.UpdatePosition();
  }

  Window window_;
  UI ui_;
  GeometryMaster geometry_master_;
  ImGuiIO& user_inputs_;
};

#endif //GEOMETRY__GEOMETRYAPP_H_
