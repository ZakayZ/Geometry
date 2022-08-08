//
// Created by Artem Novikov on 06.08.2022.
//

#ifndef GEOMETRY__GEOMETRYMASTER_H_
#define GEOMETRY__GEOMETRYMASTER_H_

#include <memory>
#include <unordered_map>

#include "Renderer.h"
#include "Geometry2D.h"
#include "Tool.h"
#include "Create/CreatePointTool.h"
#include "Create/CreateLineTool.h"

class GeometryMaster {
 public:
  GeometryMaster(size_t window_width, size_t window_height)
      : handled_geometry_(), renderer_(window_width, window_height),
        chosen_tool_(std::make_unique<CreateLineTool>(handled_geometry_)) {}
  ~GeometryMaster() = default;

  void ProcessWindowResize(const Vector2f& new_size) { /// in pixels
    renderer_.RegisterWindowScale(new_size);
  }

  void ProcessWindowShift(const Vector2f shift) {  /// in pixels
    renderer_.RegisterWindowShift(shift);
  }

  void ProcessWindowScale(float scale_change) {
    renderer_.RegisterCoordinateSystemScale(Vector2f(1.f, 1.f) * (1 + scale_change / 10.f));
  }

  Point2f ProcessHover(const Point2f& cursor_position) { /// relative to window
    auto geometry_position = renderer_.MapCursorToGeometry(cursor_position);
    return renderer_.MapCursorToWindow(
        chosen_tool_->ProcessHover(geometry_position, vicinity_ * renderer_.GetScale()));
  }

  void ProcessPressed(const Point2f& click_position) { /// relative to window
    auto geometry_position = renderer_.MapCursorToGeometry(click_position);
    chosen_tool_->ProcessPressed(geometry_position);
  }

  void ProcessDown(const Point2f& click_position) { /// relative to window
    auto geometry_position = renderer_.MapCursorToGeometry(click_position);
    chosen_tool_->ProcessDown(geometry_position);
  }

  void ProcessReleased(const Point2f& click_position) { /// relative to window
    auto geometry_position = renderer_.MapCursorToGeometry(click_position);
    chosen_tool_->ProcessReleased(geometry_position);
  }

  void Render() {
    renderer_.RenderCoordinateSystem();
    for (const auto& object : handled_geometry_.GetData()) {
      renderer_.Render(object, styles_[object]);
    }
  }

  template <class NewTool, typename = std::enable_if_t<std::is_base_of_v<Tool, NewTool>>>
  void ChangeTool() {  /// Add specialization to enable new tools that require different constructor
    chosen_tool_ = std::make_unique<NewTool>(handled_geometry_);
  }

  void ProcessRefresh() {
    chosen_tool_->Refresh();
  }

 private:
  Geometry2D<float> handled_geometry_;
  Renderer renderer_;
  std::unique_ptr<Tool> chosen_tool_;

  std::unordered_map<std::shared_ptr<Void2f>, Style> styles_;
  float vicinity_ = 0.1f;
};

#endif //GEOMETRY__GEOMETRYMASTER_H_
