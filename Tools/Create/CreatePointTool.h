//
// Created by Artem Novikov on 05.08.2022.
//

#ifndef GEOMETRY_TOOLS_CREATE_CREATEPOINTTOOL_H_
#define GEOMETRY_TOOLS_CREATE_CREATEPOINTTOOL_H_

#include "Tool.h"

class CreatePointTool : public Tool {
 public:
  CreatePointTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    auto selected_objects = Tool::output_geometry_.Selected(clicked_pos);
    if (Filter(selected_objects, {Entity::Point}).empty()) { output_geometry_.Push(clicked_pos); }
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_objects = Tool::output_geometry_.Selected(cursor_pos, vicinity);
    Point2f new_position = ClosestPoint(cursor_pos, selected_objects);
    return new_position;
  }

  void Refresh() override {}

 private:
};

#endif //GEOMETRY_TOOLS_CREATE_CREATEPOINTTOOL_H_
