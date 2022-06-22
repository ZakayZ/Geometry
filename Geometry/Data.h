//
// Created by Artem Novikov on 22.06.2022.
//

#include <array>
#include "Comparator.h"

#ifndef GEOMERTY_GEOMETRY_DATA_H_
#define GEOMERTY_GEOMETRY_DATA_H_

template <typename T, size_t dim>
union Data {
  std::array<T, dim> cords;
};

template <typename T>
union Data<T, 1> {
  std::array<T, 1> cords;
  struct {
    T x;
  };
};


template <typename T>
union Data<T, 2> {
  std::array<T, 2> cords;
  struct {
    T x;
    T y;
  };
};

template <typename T>
union Data<T, 3> {
  std::array<T, 3> cords;
  struct {
    T x;
    T y;
    T z;
  };
};

template <typename T, size_t dim>
bool operator==(const Data<T, dim>& a, const Data<T, dim>& b) {
  bool ok = true;
  for(size_t i = 0; i < dim; ++i){
    ok = ok && Comparator<T>::Equal(a.cords[i], b.cords[i]);
  }
  return ok;
}

template <typename T, size_t dim>
bool operator!=(const Data<T, dim>& a, const Data<T, dim>& b) {
  return !(a == b);
}

#endif //GEOMERTY_GEOMETRY_DATA_H_
