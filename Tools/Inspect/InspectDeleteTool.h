//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_
#define GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_

#include "Renderer.h"
#include "Tool.h"

class InspectDeleteTool : public Tool {
 public:
  InspectDeleteTool(Geometry2D<float>& output_geometry, std::unordered_map<std::shared_ptr<Void2f>, Style>& styles)
      : Tool(output_geometry), styles_(styles) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    auto selected_objects = output_geometry_.Selected(clicked_pos);
    auto selected_points = Filter(selected_objects, {Entity::Point});
    if (!selected_points.empty()) {
      selected_objects = std::move(selected_points);
    }

    if (selected_objects.empty()) { return; }

    auto object_to_delete = ClosestObject(clicked_pos, selected_objects);
    styles_.erase(object_to_delete);
    output_geometry_.Erase(object_to_delete);
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_objects = output_geometry_.Selected(cursor_pos, vicinity);
    auto selected_points = Filter(selected_objects, {Entity::Point});
    if (!selected_points.empty()) {
      selected_objects = std::move(selected_points);
    }

    auto new_position = ClosestPoint(cursor_pos, selected_objects);
    return new_position;
  }

 private:
  std::unordered_map<std::shared_ptr<Void2f>, Style> styles_;
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTDELETETOOL_H_
