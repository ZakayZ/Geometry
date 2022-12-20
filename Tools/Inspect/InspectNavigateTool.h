//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY_TOOLS_INSPECT_INSPECTNAVIGATETOOL_H_
#define GEOMETRY_TOOLS_INSPECT_INSPECTNAVIGATETOOL_H_

#include "Tool.h"
#include "Renderer.h"

class InspectNavigateTool : public Tool {
 public:
  InspectNavigateTool(Geometry2D<float>& output_geometry, Renderer& move_renderer)
      : Tool(output_geometry), move_renderer_(move_renderer) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    initial_position_ = last_position_ = move_renderer_.MapCursorToWindow(clicked_pos);
    tracking = true;
  }

  void ProcessDown(const Point2f& cursor_pos) override {
    if (tracking) {
      auto screen_pos = move_renderer_.MapCursorToWindow(cursor_pos);
      move_renderer_.RegisterCoordinateSystemShift(last_position_ - screen_pos);
      last_position_ = screen_pos;
    }
  }

  void ProcessReleased(const Point2f& cursor_pos) override {
    ProcessDown(cursor_pos);
    tracking = false;
  }

  void Refresh() override {
    move_renderer_.RegisterCoordinateSystemShift(last_position_ - initial_position_);
    tracking = false;
  }

 private:
  Renderer& move_renderer_;
  Point2f last_position_;
  Point2f initial_position_;
  bool tracking = false;
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTNAVIGATETOOL_H_
