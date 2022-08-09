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

  virtual ~InspectNavigateTool() = default;

  virtual void ProcessPressed(const Point2f& clicked_pos) {
    last_position_ = clicked_pos;
    initial_position_ = clicked_pos;
    tracking = true;
  }

  virtual void ProcessDown(const Point2f& cursor_pos) {
    if (tracking) {
      move_renderer_.RegisterCoordinateSystemShift(cursor_pos - last_position_);
      last_position_ = cursor_pos;
    }
  }

  virtual void ProcessReleased(const Point2f& cursor_pos) {
    ProcessDown(cursor_pos);
    tracking = false;
  }

  virtual void Refresh() {
    ProcessReleased(initial_position_);
  }

 private:
  Renderer& move_renderer_;
  Point2f last_position_;
  Point2f initial_position_;
  bool tracking = false;
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTNAVIGATETOOL_H_
