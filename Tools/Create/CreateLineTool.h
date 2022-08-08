//
// Created by Artem Novikov on 05.08.2022.
//

#ifndef GEOMETRY_TOOLS_CREATE_CREATELINETOOL_H_
#define GEOMETRY_TOOLS_CREATE_CREATELINETOOL_H_

#include "Tool.h"

class CreateLineTool : public Tool {
 public:
  CreateLineTool(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    if (line_ == nullptr) {
      output_geometry_.Push(clicked_pos);
      origin_ = std::static_pointer_cast<Point2f>(*--output_geometry_.end());
      output_geometry_.Push(Line2f(clicked_pos, clicked_pos));
      line_ = std::static_pointer_cast<Line2f>(*--output_geometry_.end());
    } else {
      output_geometry_.Push(clicked_pos);
      origin_ = nullptr;
      line_ = nullptr;
    }
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_objects = output_geometry_.Selected(cursor_pos, vicinity);
    auto itself_iter = std::find(selected_objects.begin(), selected_objects.end(), line_);
    if (line_ != nullptr && itself_iter != selected_objects.end()) {
      selected_objects.erase(itself_iter);
    }

    auto selected_points = Filter(selected_objects, {Entity::Point});
    if (line_ != nullptr && !selected_points.empty()) {
      selected_objects = std::move(selected_points); /// prioritize points
    }

    auto new_position = ClosestPoint(cursor_pos, selected_objects);
    if (line_ != nullptr && new_position != line_->GetOrigin()) {
      line_->GetDirection() = new_position - line_->GetOrigin();
    }
    return new_position;
  }

  void Refresh() override {
    if (line_ != nullptr) {
      output_geometry_.Erase(origin_);
      output_geometry_.Erase(line_);
    }
    origin_ = nullptr;
    line_ = nullptr;
  }
 private:
  std::shared_ptr<Point2f> origin_;
  std::shared_ptr<Line2f> line_;
};

#endif //GEOMETRY_TOOLS_CREATE_CREATELINETOOL_H_
