//
// Created by Artem Novikov on 16.05.2022.
//

#include <array>
#include "Void.h"
#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_POINT_H_
#define GEOMERTY_GEOMETRY_POINT_H_

enum class PointRelationship {
  Identical,
  None,
};

template <typename T, size_t dim>
class Point : public Void<T, dim>, public Vector<T, dim> {
 public:
  /// construction
  Point() : Void<T, dim>(Entity::Point) {}
  template <typename... Args>
  Point(Args&& ... args);
  Point(std::initializer_list<T> list);
  Point(const Point<T, dim>& other) = default;
  Point(Point<T, dim>&& other) noexcept = default;
  explicit Point(const Vector<T, dim>& vec) : Void<T, dim>::type_(Entity::Point), Vector<T, dim>(vec) {}
  explicit Point(Vector<T, dim>&& vec) noexcept: Void<T, dim>::type_(Entity::Point), Vector<T, dim>(std::move(vec)) {}
  Point& operator=(const Point& other) = default;
  Point& operator=(Point&& other) noexcept = default;

  /// TODO has different affine transformation

  T SquaredDistance(const Point<T, dim>& point) const;

  T Distance(const Point<T, dim>& point) const;

  std::unique_ptr<Void<T, dim>> Intersection(const Point<T, dim>& point) const;
};

template <typename T>
using Point1 = Point<T, 1>;

using Point1i = Point1<int>;
using Point1f = Point1<float>;
using Point1d = Point1<double>;

template <typename T>
using Point2 = Point<T, 2>;

using Point2i = Point2<int>;
using Point2f = Point2<float>;
using Point2d = Point2<double>;

template <typename T>
using Point3 = Point<T, 3>;

using Point3i = Point3<int>;
using Point3f = Point3<float>;
using Point3d = Point3<double>;

template <typename T, size_t dim>
PointRelationship FindRelationship(const Point<T, dim>& a, const Point<T, dim>& b);

/////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////////

template <typename T, size_t dim>
template <typename... Args>
Point<T, dim>::Point(Args&& ... args)
    : Void<T, dim>(Entity::Point), Vector<T, dim>(std::forward<Args>(args)...) {}

template <typename T, size_t dim>
Point<T, dim>::Point(std::initializer_list<T> list) : Void<T, dim>(Entity::Point), Vector<T, dim>(list) {}

template <typename T, size_t dim>
T Point<T, dim>::SquaredDistance(const Point<T, dim>& point) const {
  return (point - *this).SquaredLength();
}

template <typename T, size_t dim>
T Point<T, dim>::Distance(const Point<T, dim>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Point<T, dim>::Intersection(const Point<T, dim>& point) const {
  if (point == *this) {
    return std::unique_ptr<Point<T, dim>>(point);
  }
  return std::unique_ptr<Void<T, dim>>();
}

template <typename T, size_t dim>
PointRelationship FindRelationship(const Point<T, dim>& a, const Point<T, dim>& b) {
  return a == b ? PointRelationship::Identical : PointRelationship::None;
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
