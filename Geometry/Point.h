//
// Created by Artem Novikov on 16.05.2022.
//

#include <array>
#include "Comparator.h"

#ifndef GEOMERTY_GEOMETRY_POINT_H_
#define GEOMERTY_GEOMETRY_POINT_H_

template <typename T, size_t dim>
union Point {
  std::array<T, dim> cords;
};

template <typename T>
union Point<T, 1> {
  std::array<T, 1> cords;
  struct {
    T x;
  };
};

template <typename T>
using Point1 = Point<T, 1>;

using Point1i = Point1<int>;
using Point1f = Point1<float>;
using Point1d = Point1<double>;

template <typename T>
union Point<T, 2> {
  std::array<T, 2> cords;
  struct {
    T x;
    T y;
  };
};

template <typename T>
using Point2 = Point<T, 2>;

using Point2i = Point2<int>;
using Point2f = Point2<float>;
using Point2d = Point2<double>;

template <typename T>
union Point<T, 3> {
  std::array<T, 3> cords;
  struct {
    T x;
    T y;
    T z;
  };
};

template <typename T>
using Point3 = Point<T, 3>;

using Point3i = Point3<int>;
using Point3f = Point3<float>;
using Point3d = Point3<double>;

template <typename T, size_t dim>
bool operator==(const Point<T, dim>& a, const Point<T, dim>& b) {
  bool ok = true;
  for(size_t i = 0; i < dim; ++i){
    ok = ok && Comparator<T>::Equal(a.cords[i], b.cords[i]);
  }
  return ok;
}

template <typename T, size_t dim>
bool operator!=(const Point<T, dim>& a, const Point<T, dim>& b) {
  return !(a == b);
}

#endif //GEOMERTY_GEOMETRY_POINT_H_
