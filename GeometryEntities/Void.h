//
// Created by Artem Novikov on 28.06.2022.
//

#include <cstddef>
#include "Entities.h"

#ifndef GEOMERTY_GEOMETRY_VOID_H_
#define GEOMERTY_GEOMETRY_VOID_H_

template <typename T, size_t Dimension>
class Void {
 public:
  Void() = default;
  Entity GetType() const{ return Entity::Void; }
  size_t GetDimension() const{ return Dimension; }
  ~Void() = default;
};

template <typename T>
using Void2 = Void<T, 2>;

using Void2i = Void2<int>;
using Void2f = Void2<float>;
using Void2d = Void2<double>;

template <typename T>
using Void3 = Void<T, 3>;

using Void3i = Void3<int>;
using Void3f = Void3<float>;
using Void3d = Void3<double>;

#endif //GEOMERTY_GEOMETRY_VOID_H_
