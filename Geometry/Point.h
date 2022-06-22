//
// Created by Artem Novikov on 16.05.2022.
//

#include <array>
#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_POINT_H_
#define GEOMERTY_GEOMETRY_POINT_H_

template <typename T, size_t dim>
class Point : public Vector<T, dim> {
  using Vector<T, dim>::Vector;
  using Vector<T, dim>::operator=;
  /// TODO has different transformation
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
T SquaredDistance(const Point<T, dim>& a, const Point<T, dim>& b);

template <typename T, size_t dim>
T Distance(const Point<T, dim>& a, const Point<T, dim>& b);

/////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////////

template <typename T, size_t dim>
T SquaredDistance(const Point<T, dim>& a, const Point<T, dim>& b) {
  T result = (a - b).SquaredLength();
  return result;
}

template <typename T, size_t dim>
T Distance(const Point<T, dim>& a, const Point<T, dim>& b) {
  T result = (a - b).Length();
  return result;
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
