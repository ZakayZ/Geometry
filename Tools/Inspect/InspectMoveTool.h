//
// Created by Artem Novikov on 09.08.2022.
//

#ifndef GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_
#define GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_

#include "Tool.h"

class InspectMoveTool : public Tool {
 public:
  InspectMoveTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    start_position_ = clicked_pos;
    last_position_ = clicked_pos;
    auto selected_objects = output_geometry_.Selected(clicked_pos);
    auto selected_points = Filter(selected_objects, {Entity::Point});
    if (selected_objects.empty()) { return; }

    if (!selected_points.empty()) {
      pivot_points.push_back(std::static_pointer_cast<Point2f>(ClosestObject(clicked_pos, selected_points)));
    } else {
      auto closest_object = ClosestObject(clicked_pos, selected_objects);
      auto dependencies = FindDependencies(closest_object);
      for (auto& dependency : dependencies) {
        pivot_points.push_back(dependency);
      }
    }

    for (auto& point : pivot_points) {
      connected_objects_.push_back(output_geometry_.Selected(*point));
      connected_objects_.back().erase(
          std::find(connected_objects_.back().begin(), connected_objects_.back().end(), point));
    }
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_points = Filter(output_geometry_.Selected(cursor_pos, vicinity), {Entity::Point});
    auto closest_point = ClosestObject(cursor_pos, selected_points);
    if (closest_point != nullptr
        && std::find(pivot_points.begin(), pivot_points.end(), closest_point) == pivot_points.end()) {
      return static_cast<Point2f&>(*closest_point);
    }

    return cursor_pos;
  }

  void ProcessDown(const Point2f& cursor_pos) override {
    Vector2f shift = cursor_pos - last_position_;
    last_position_ = cursor_pos;
    for (size_t point_index = 0; point_index < pivot_points.size(); ++point_index) {
      for (auto& object : connected_objects_[point_index]) {
        MoveIfPivot(*pivot_points[point_index], object, shift);
      }
      *pivot_points[point_index] += shift;
    }
  }

  void ProcessReleased(const Point2f& cursor_pos) override {
    ProcessDown(cursor_pos);
    pivot_points.clear();
    connected_objects_.clear();
  }

  void Refresh() override {
    ProcessReleased(start_position_);
  }

 private:
  Point2f start_position_;
  Point2f last_position_;
  std::vector<std::shared_ptr<Point2f>> pivot_points;
  std::vector<std::list<std::shared_ptr<Void2f>>> connected_objects_;
};

#endif //GEOMETRY_TOOLS_INSPECT_INSPECTMOVETOOL_H_
