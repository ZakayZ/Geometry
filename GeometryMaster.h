//
// Created by Artem Novikov on 06.08.2022.
//

#ifndef GEOMETRY__GEOMETRYMASTER_H_
#define GEOMETRY__GEOMETRYMASTER_H_

#include <memory>
#include <unordered_map>

#include "CoordinateSystem.h"
#include "Renderer.h"
#include "Geometry2D.h"
#include "Tools.h"
#include "Tool.h"
#include "Inspect/InspectMoveTool.h"
#include "Inspect/InspectNavigateTool.h"
#include "Inspect/InspectDeleteTool.h"

#include "Construct/ConstructMidpoint.h"
#include "Construct/ConstructProjection.h"

#include "Create/CreateLineTool.h"
#include "Create/CreatePointTool.h"
#include "Create/CreateSegmentTool.h"

class GeometryMaster {
 public:
  GeometryMaster(size_t window_width, size_t window_height)
      : handled_geometry_(), renderer_(window_width, window_height),
        chosen_tool_(nullptr) { ChangeTool(ToolType::Inspect_Navigate); }

  ~GeometryMaster() = default;

  void ProcessWindowResize(const Vector2f& new_size) { /// in pixels
    renderer_.RegisterWindowScale(new_size);
  }

  void ProcessWindowShift(const Vector2f shift) {  /// in pixels
    renderer_.RegisterWindowShift(shift);
  }

  void ProcessWindowScale(float scale) {
    renderer_.RegisterCoordinateSystemScale(Vector2f(1.f, 1.f) * scale);
  }

  Point2f ProcessHover(const Point2f& cursor_position) { /// relative to window
    auto geometry_position = renderer_.MapCursorToGeometry(cursor_position);
    return renderer_.MapCursorToWindow(
        chosen_tool_->ProcessHover(geometry_position, renderer_.GetPointSize()));
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
    renderer_.RenderCoordinateSystem(coordinate_system_);
    for (const auto& object : handled_geometry_.GetData()) {
      renderer_.Render(object, styles_[object]);
    }
  }

  void ChangeTool(ToolType new_type) {  /// Add specialization to enable new tools that require different constructor
    switch (new_type) {
      case ToolType::None: {
        break;
      }
      case ToolType::Inspect_Navigate: {
        chosen_tool_ = std::make_unique<InspectNavigateTool>(handled_geometry_, renderer_);
        break;
      }
      case ToolType::Inspect_Move: {
        chosen_tool_ = std::make_unique<InspectMoveTool>(handled_geometry_);
        break;
      }
      case ToolType::Inspect_Delete: {
//        chosen_tool_ = std::make_unique<InspectDeleteTool>(handled_geometry_);
        break;
      }
      case ToolType::Construct_Midpoint: {
        chosen_tool_ = std::make_unique<ConstructMidpoint>(handled_geometry_);
        break;
      }
      case ToolType::Construct_Project: {
        chosen_tool_ = std::make_unique<ConstructProjection>(handled_geometry_);
        break;
      }
      case ToolType::Create_Point: {
        chosen_tool_ = std::make_unique<CreatePointTool>(handled_geometry_);
        break;
      }
      case ToolType::Create_Line: {
        chosen_tool_ = std::make_unique<CreateLineTool>(handled_geometry_);
        break;
      }
      case ToolType::Create_Segment: {
        chosen_tool_ = std::make_unique<CreateSegmentTool>(handled_geometry_);
        break;
      }
    }
  }

  void ProcessRefresh() {
    chosen_tool_->Refresh();
  }

 private:
  Geometry2D<float> handled_geometry_;
  Renderer renderer_;
  std::unique_ptr<Tool> chosen_tool_;
  CoordinateSystem coordinate_system_;

  std::unordered_map<std::shared_ptr<Void2f>, Style> styles_;
  float vicinity_ = 0.1f;
};

#endif //GEOMETRY__GEOMETRYMASTER_H_
