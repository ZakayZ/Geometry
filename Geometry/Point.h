//
// Created by Artem Novikov on 16.05.2022.
//

#include <array>
#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_POINT_H_
#define GEOMERTY_GEOMETRY_POINT_H_

template <typename T, size_t dim>
class Point : public Vector<T, dim> {
 public:
  enum class Relationship {
    Identical,
    None,
  };

  using Vector<T, dim>::Vector;
  using Vector<T, dim>::operator=;

  explicit Point(const Vector<T, dim>& v);
  /// TODO has different transformation

  T SquaredDistance(const Point<T, dim>& point) const;

  T Distance(const Point<T, dim>& point) const;
};

template <typename T, size_t dim>
Point<T, dim>::Point(const Vector<T, dim>& v) : Point::Vector(v) {}

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
typename Point<T, dim>::Relationship FindRelationShip(const Point<T, dim>& a, const Point<T, dim>& b);

/////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////////

template <typename T, size_t dim>
T Point<T, dim>::SquaredDistance(const Point<T, dim>& point) const  {
  T result = (point - *this).SquaredLength();
  return result;
}

template <typename T, size_t dim>
T Point<T, dim>::Distance(const Point<T, dim>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t dim>
typename Point<T, dim>::Relationship FindRelationShip(const Point<T, dim>& a, const Point<T, dim>& b) {
  return a == b ? Point<T, dim>::Relationship::Identical : Point<T, dim>::Relationship::None;
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
