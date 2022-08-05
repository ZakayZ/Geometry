//
// Created by Artem Novikov on 05.08.2022.
//

#ifndef GEOMETRY_TOOLS_TOOL_H_
#define GEOMETRY_TOOLS_TOOL_H_

#include "../Geometry2D.h"

class Tool {
 public:
  Tool(Geometry2D<float>& output_geometry) : output_geometry_(output_geometry) {}
  virtual ~Tool() = 0;

  virtual void ProcessInput(const Point2f& clicked_pos) = 0;

  virtual Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) = 0;

  virtual void Refresh() = 0;

 protected:
  static Point2f ClosestPoint(const Point2f& point, const std::list<std::shared_ptr<Void2f>>& objects) {
    float min_distance = INFINITY;
    Point2f closest_point = point;
    for (const auto& object_ptr : objects) {
      Point2f snapped_point;
      switch (object_ptr->GetType()) {
        case Entity::Point: {
          snapped_point = *static_cast<Point2f*>(object_ptr.get());
        }
        case Entity::Line: {
          snapped_point = static_cast<Line2f*>(object_ptr.get())->Projection(point);
        }
        case Entity::Segment: {
          snapped_point = static_cast<Segment2f*>(object_ptr.get())->Projection(point);
        } /// TODO other types
        default: {
          snapped_point = point;
          snapped_point[0] += min_distance;
        }
      }
      auto distance = snapped_point.Distance(point);
      if (distance < min_distance) {
        min_distance = distance;
        closest_point = snapped_point;
      }
    }
    return closest_point;
  }

  static std::list<std::shared_ptr<Void2f>> Filter(
      const std::list<std::shared_ptr<Void2f>>& objects, const std::list<Entity>& types) {
    std::list<std::shared_ptr<Void2f>> filtered_objects;
    for (const auto& object : objects) {
      for (auto type : types) {
        if (object->GetType() == type) {
          filtered_objects.push_back(object);
          break;
        }
      }
    }
    return filtered_objects;
  }

  Geometry2D<float>& output_geometry_;
};

#endif //GEOMETRY_TOOLS_TOOL_H_