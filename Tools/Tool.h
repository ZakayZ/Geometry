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

  virtual void ProcessPressed(const Point2f& clicked_pos) = 0;

  virtual Point2f ProcessHover(const Point2f& cursor_pos, float vicinity) { return cursor_pos; }

  virtual void ProcessDown(const Point2f& cursor_pos) {}

  virtual void ProcessReleased(const Point2f& cursor_pos) {}

  virtual void Refresh() {}

 protected:
  static Point2f ClosestPoint(const Point2f& point, const std::list<std::shared_ptr<Void2f>>& objects) {
    float min_distance = INFINITY;
    Point2f closest_point = point;
    for (const auto& object_ptr : objects) {
      Point2f snapped_point;
      switch (object_ptr->GetType()) {
        case Entity::Point: {
          snapped_point = *static_cast<Point2f*>(object_ptr.get());
          break;
        }
        case Entity::Line: {
          snapped_point = static_cast<Line2f*>(object_ptr.get())->Projection(point);
          break;
        }
        case Entity::Segment: {
          snapped_point = static_cast<Segment2f*>(object_ptr.get())->Projection(point);
          break;
        } /// TODO other types
        default: {
          snapped_point = point;
          snapped_point[0] += 3 * min_distance;
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

  static std::shared_ptr<Void2f> ClosestObject(
      const Point2f& point, const std::list<std::shared_ptr<Void2f>>& objects) {
    float min_distance = INFINITY;
    std::shared_ptr<Void2f> closest_object = nullptr;
    for (const auto& object_ptr : objects) {
      auto distance = object_ptr->Distance(point);
      if (distance < min_distance) {
        min_distance = distance;
        closest_object = object_ptr;
      }
    }
    return closest_object;
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

  static void MoveIfPivot(const Point2f& pivot_point, std::shared_ptr<Void2f>& object, const Vector2f& shift) {
    switch (object->GetType()) {
      case Entity::Point: {
        auto& point = static_cast<Point2f&>(*object);
        if (point == pivot_point) {
          point += shift;
        }
        break;
      }

      case Entity::Line: {
        auto& line = static_cast<Line2f&>(*object);
        if (line.GetOrigin() == pivot_point) {
          line.GetDirection() -= shift;
          line.GetOrigin() += shift;
        }

        if (line.GetPoint(1) == pivot_point) {
          line.GetDirection() += shift;
        }
        break;
      }

      case Entity::Segment: {
        auto& segment = static_cast<Segment2f&>(*object);
        if (segment.GetLeft() == pivot_point) {
          segment.GetLeft() += shift;
        }

        if (segment.GetRight() == pivot_point) {
          segment.GetRight() += shift;
        }
        break;
      }
    }
  }

  std::list<std::shared_ptr<Point2f>> FindDependencies(std::shared_ptr<Void2f>& object) {
    std::list<std::shared_ptr<Point2f>> dependencies;
    auto points = Filter(output_geometry_.GetData(), {Entity::Point});
    switch (object->GetType()) {
      case Entity::Point: {
        break;
      }

      case Entity::Line: {
        auto line = static_cast<Line2f&>(*object);
        auto point1 = ClosestObject(line.GetPoint(1.f), points); /// order matters
        if (point1 != nullptr) { dependencies.push_back(std::static_pointer_cast<Point2f>(point1)); }

        auto point2 = ClosestObject(line.GetOrigin(), points);
        if (point2 != nullptr) { dependencies.push_back(std::static_pointer_cast<Point2f>(point2)); }
        break;
      }

      case Entity::Segment: {
        auto segment = static_cast<Segment2f&>(*object);
        auto point1 = ClosestObject(segment.GetLeft(), points);
        if (point1 != nullptr) { dependencies.push_back(std::static_pointer_cast<Point2f>(point1)); }

        auto point2 = ClosestObject(segment.GetRight(), points);
        if (point2 != nullptr) { dependencies.push_back(std::static_pointer_cast<Point2f>(point2)); }
        break;
      }
    }
    return dependencies;
  }

  Geometry2D<float>& output_geometry_;
};

Tool::~Tool() {}

#endif //GEOMETRY_TOOLS_TOOL_H_
