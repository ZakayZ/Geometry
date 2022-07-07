//
// Created by Artem Novikov on 16.05.2022.
//

#include <array>
#include "GeometricEntity.h"
#include "Void.h"
#include "Vector.h"
#include "Transform.h"

#ifndef GEOMERTY_GEOMETRY_POINT_H_
#define GEOMERTY_GEOMETRY_POINT_H_

enum class PointRelationship {
  Identical,
  None,
};

template <typename T, size_t Dimension>
class Point : public Vector<T, Dimension> {
 public:
  /// construction
  Point() {}
  template <typename... Args>
  Point(Args&& ... args);
  template <typename U, template <typename, typename...> class Container, typename... Args,
      typename = std::enable_if_t<is_iterable_v<Container<U, Args...>>>>
  explicit Point(const Container<U, Args...>& data);
  Point(std::initializer_list<T> list);
  Point(const Point<T, Dimension>& other) = default;
  Point(Point<T, Dimension>&& other) noexcept = default;
  explicit Point(const Vector<T, Dimension>& vec) : Vector<T, Dimension>(vec) {}
  explicit Point(Vector<T, Dimension>&& vec) noexcept: Vector<T, Dimension>(std::move(vec)) {}
  Point& operator=(const Point& other) = default;
  Point& operator=(Point&& other) noexcept = default;

  /// getters and setters
  Entity GetType() const { return Entity::Point; }
  size_t GetDimension() const { return Dimension; }

  /// calc

  T SquaredDistance(const Point<T, Dimension>& point) const;

  T Distance(const Point<T, Dimension>& point) const;

  GeometryEntity Intersection(const Point<T, Dimension>& point) const;

  template <size_t OutputDimension>
  Point<T, OutputDimension> ApplyTransform(const Transform<T, Dimension, OutputDimension>& transform) const;
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

template <typename T, size_t Dimension>
PointRelationship FindRelationship(const Point<T, Dimension>& a, const Point<T, Dimension>& b);

/////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////////

template <typename T, size_t Dimension>
template <typename... Args>
Point<T, Dimension>::Point(Args&& ... args) : Vector<T, Dimension>(std::forward<Args>(args)...) {}

template <typename T, size_t Dimension>
template <typename U, template <typename, typename...> class Container, typename... Args, typename>
Point<T, Dimension>::Point(const Container<U, Args...>& data) : Vector<T, Dimension>(data) {}

template <typename T, size_t Dimension>
Point<T, Dimension>::Point(std::initializer_list<T> list) : Vector<T, Dimension>(list) {}

template <typename T, size_t Dimension>
T Point<T, Dimension>::SquaredDistance(const Point<T, Dimension>& point) const {
  return (point - *this).SquaredLength();
}

template <typename T, size_t Dimension>
T Point<T, Dimension>::Distance(const Point<T, Dimension>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t Dimension>
GeometryEntity Point<T, Dimension>::Intersection(const Point<T, Dimension>& point) const {
  if (point == *this) {
    return GeometryEntity(point);
  }
  return MakeGeometryEntity<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
template <size_t OutputDimension>
Point<T, OutputDimension> Point<T, Dimension>::ApplyTransform(
    const Transform<T, Dimension, OutputDimension>& transform) const {
  return transform.GetMatrix() * (*this) + transform.GetShift();
}

template <typename T, size_t Dimension>
PointRelationship FindRelationship(const Point<T, Dimension>& a, const Point<T, Dimension>& b) {
  return a == b ? PointRelationship::Identical : PointRelationship::None;
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
