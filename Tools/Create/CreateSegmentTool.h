//
// Created by Artem Novikov on 05.08.2022.
//

#ifndef GEOMETRY_TOOLS_CREATE_CREATESEGMENTTOOL_H_
#define GEOMETRY_TOOLS_CREATE_CREATESEGMENTTOOL_H_

#include "Tool.h"

class CreateLineTool : public Tool {
 public:
  CreateLineTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessInput(const Point2f& clicked_pos) override {
    if (segment_ == nullptr) {
      output_geometry_.Push(Point2f(clicked_pos));
      first_point_ = std::static_pointer_cast<Point2f>(*--output_geometry_.end());
      output_geometry_.Push(Segment2f(clicked_pos, clicked_pos));
      segment_ = std::static_pointer_cast<Segment2f>(*--output_geometry_.end());
    } else {
      first_point_ = nullptr;
      output_geometry_.Push(segment_->GetRight());
      segment_ = nullptr;
    }
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_objects = Tool::output_geometry_.Selected(cursor_pos, vicinity);
    auto selected_points = Filter(selected_objects, {Entity::Point});
    if (segment_ != nullptr && !selected_points.empty()) {
      selected_objects = std::move(selected_points);
    }
    auto new_position = ClosestPoint(cursor_pos, selected_objects);
    if (segment_ != nullptr) {
      segment_->GetRight() = new_position;
    }
    return new_position;
  }

  void Refresh() override {
    if (segment_ != nullptr) {
      output_geometry_.Erase(first_point_);
      output_geometry_.Erase(segment_);
    }
    first_point_ = nullptr;
    segment_ = nullptr;
  }

 private:
  std::shared_ptr<Point2f> first_point_;
  std::shared_ptr<Segment2f> segment_;
};

#endif //GEOMETRY_TOOLS_CREATE_CREATESEGMENTTOOL_H_
