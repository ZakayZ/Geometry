//
// Created by Artem Novikov on 06.08.2022.
//

#ifndef GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTPROJECTION_H_
#define GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTPROJECTION_H_

#include "Tool.h"

class ConstructProjection : public Tool {
 public:
  ConstructProjection(Geometry2D<float>& output_geometry) : Tool(output_geometry) {}

  void ProcessPressed(const Point2f& clicked_pos) override {
    auto selected_objects = output_geometry_.Selected(clicked_pos);
    selected_objects = object_to_project_ != nullptr
                       ? Filter(selected_objects, {Entity::Line, Entity::Segment})
                       : Filter(selected_objects, {Entity::Point, Entity::Segment});
    if (selected_objects.empty()) { return; }
    auto closest_object = ClosestObject(clicked_pos, selected_objects);
    if (object_to_project_ == nullptr) {
      object_to_project_ = closest_object;
      return;
    }

    switch (object_to_project_->GetType()) {    /// TODO virtual project
      case Entity::Point: {
        if (closest_object->GetType() == Entity::Segment) {
          output_geometry_.Push(static_cast<const Segment2f&>(*closest_object)
                                    .Projection(static_cast<const Point2f&>(*object_to_project_)));
        } else {
          output_geometry_.Push(static_cast<const Line2f&>(*closest_object)
                                    .Projection(static_cast<const Point2f&>(*object_to_project_)));
        }
        break;
      }
      case Entity::Segment: {
        if (closest_object->GetType() == Entity::Segment) {
          output_geometry_.Push(static_cast<const Segment2f&>(*closest_object)
                                    .Projection(static_cast<const Segment2f&>(*object_to_project_)));
        } else {
          output_geometry_.Push(static_cast<const Line2f&>(*closest_object)
                                    .Projection(static_cast<const Point2f&>(*object_to_project_)));
        }
        break;
      }
    }
    object_to_project_ = nullptr;
  }

  Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) override {
    auto selected_objects = output_geometry_.Selected(cursor_pos, vicinity);
    if (object_to_project_ != nullptr) {
      selected_objects = Filter(selected_objects, {Entity::Line, Entity::Segment});
    } else {
      selected_objects = Filter(selected_objects, {Entity::Point, Entity::Segment});
    }
    return ClosestPoint(cursor_pos, selected_objects);;
  }

  void Refresh() override {
    object_to_project_ = nullptr;
  }

 private:
  std::shared_ptr<Void2f> object_to_project_;
};

#endif //GEOMETRY_TOOLS_CONSTRUCT_CONSTRUCTPROJECTION_H_
