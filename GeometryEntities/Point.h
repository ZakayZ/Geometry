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

template <typename T, size_t Dimension>
class Point : public Void<T, Dimension>, public Vector<T, Dimension> {
 public:
  /// construction
  Point() : Void<T, Dimension>(Entity::Point) {}
  template <typename... Args>
  Point(Args&& ... args);
  template <typename U, template <typename, typename...> class Container, typename... Args,
      typename = std::enable_if_t<is_iterable_v<Container<U, Args...>>>>
  explicit Point(const Container<U, Args...>& data);
  Point(std::initializer_list<T> list);
  Point(const Point<T, Dimension>& other) = default;
  Point(Point<T, Dimension>&& other) noexcept = default;
  explicit Point(const Vector<T, Dimension>& vec) : Void<T, Dimension>::type_(Entity::Point),
                                                    Vector<T, Dimension>(vec) {}
  explicit Point(Vector<T, Dimension>&& vec) noexcept: Void<T, Dimension>::type_(Entity::Point),
                                                       Vector<T, Dimension>(std::move(vec)) {}
  Point& operator=(const Point& other) = default;
  Point& operator=(Point&& other) noexcept = default;

  /// TODO has different affine transformation

  T SquaredDistance(const Point<T, Dimension>& point) const;

  T Distance(const Point<T, Dimension>& point) const;

  std::unique_ptr<Void<T, Dimension>> Intersection(const Point<T, Dimension>& point) const;
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
Point<T, Dimension>::Point(Args&& ... args)
    : Void<T, Dimension>(Entity::Point), Vector<T, Dimension>(std::forward<Args>(args)...) {}

template <typename T, size_t Dimension>
template <typename U, template <typename, typename...> class Container, typename... Args, typename>
Point<T, Dimension>::Point(const Container<U, Args...>& data)
    : Void<T, Dimension>(Entity::Point), Vector<T, Dimension>(data) {}

template <typename T, size_t Dimension>
Point<T, Dimension>::Point(std::initializer_list<T> list) : Void<T, Dimension>(Entity::Point),
                                                            Vector<T, Dimension>(list) {}

template <typename T, size_t Dimension>
T Point<T, Dimension>::SquaredDistance(const Point<T, Dimension>& point) const {
  return (point - *this).SquaredLength();
}

template <typename T, size_t Dimension>
T Point<T, Dimension>::Distance(const Point<T, Dimension>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Point<T, Dimension>::Intersection(const Point<T, Dimension>& point) const {
  if (point == *this) {
    return std::unique_ptr<Point<T, Dimension>>(point);
  }
  return std::unique_ptr<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
PointRelationship FindRelationship(const Point<T, Dimension>& a, const Point<T, Dimension>& b) {
  return a == b ? PointRelationship::Identical : PointRelationship::None;
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
