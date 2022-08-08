//
// Created by Artem Novikov on 05.08.2022.
//

#ifndef GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTMIDPOINT_H_
#define GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTMIDPOINT_H_

#include "Tool.h"

class CreateMidpoint : public Tool {
 public:
  CreateMidpoint(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    auto selected_points = Filter(output_geometry_.Selected(clicked_pos), {Entity::Point});
    if (selected_points.empty()) { return; }

    std::shared_ptr<Point2f> new_point;
    float min_distance = INFINITY;
    for (const auto& point_ptr : selected_points) {
      auto distance = clicked_pos.Distance(static_cast<const Point2f&>(*point_ptr));
      if (distance < min_distance) {
        min_distance = distance;
        new_point = std::static_pointer_cast<Point2f>(point_ptr);
      }
    }
    if (first_point_ == nullptr) {
      first_point_ = new_point;
    } else {
      output_geometry_.Push(Point2f((*first_point_ + *new_point) * 0.5f));
      first_point_ = nullptr;
    }
  }

  Point2f ProcessHover(const Point2f& cursor_pos) override {
    auto selected_objects = output_geometry_.Selected(cursor_pos);
    selected_objects = Filter(selected_objects, {Entity::Point});
    return ClosestPoint(cursor_pos, selected_objects);;
  }

  void Refresh() override {
    first_point_ = nullptr;
  }
 private:
  std::shared_ptr<Point2f> first_point_;
};

#endif //GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTMIDPOINT_H_
