//
// Created by Artem Novikov on 07.08.2022.
//

#ifndef GEOMETRY__GEOMETRYAPP_H_
#define GEOMETRY__GEOMETRYAPP_H_

#include <unordered_map>

#include "Window.h"
#include "GeometryMaster.h"
#include "UI.h"
#include "UserInputs.h"

class GeometryApp {
 public:
  GeometryApp(int window_width, int widow_height)
      : window_(window_width, widow_height), ui_(window_.ShareWindow()),
        geometry_master_(window_width, widow_height),
        user_inputs_(ui_.GetIO(), window_.ShareWindow()) {
    ui_.SetUserInputs(user_inputs_);
  }

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
      if (user_inputs_.HasMouseInput()) {
        auto mouse_pos = user_inputs_.GetMousePos();
        auto snapped_mouse_pos = geometry_master_.ProcessHover(mouse_pos);
        if (user_inputs_.HasClicked()) {
          geometry_master_.ProcessPressed(snapped_mouse_pos);
        }
        if (user_inputs_.IsDown()) {
          geometry_master_.ProcessDown(snapped_mouse_pos);
        }
        if (user_inputs_.HasReleased()) {
          geometry_master_.ProcessReleased(snapped_mouse_pos);
        }

        /// Mouse wheel
        if (user_inputs_.HasScaleInput()) {
          geometry_master_.ProcessWindowScale(user_inputs_.GetScale(), snapped_mouse_pos);
        }
      }

      /// ESC button
      if (user_inputs_.HasKeyboardInput()) {
        if (user_inputs_.HasRefreshInput()) {
          geometry_master_.ProcessRefresh();
        }
      }
    }

    if (user_inputs_.HasChangeTool()) {
      geometry_master_.ChangeTool(user_inputs_.GetToolType());
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
  UserInputs user_inputs_;
};

#endif //GEOMETRY__GEOMETRYAPP_H_
