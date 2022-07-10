//
// Created by Artem Novikov on 22.06.2022.
//

#include <array>
#include "Comparator.h"

#ifndef GEOMETRY_GEOMETRY_DATA_H_
#define GEOMETRY_GEOMETRY_DATA_H_

template <typename T, size_t Dimension>
union Data {
  std::array<T, Dimension> cords;
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

template <typename T, size_t Dimension>
bool operator==(const Data<T, Dimension>& a, const Data<T, Dimension>& b) {
  bool ok = true;
  for(size_t i = 0; i < Dimension; ++i){
    ok = ok && Comparator<T>::Equal(a.cords[i], b.cords[i]);
  }
  return ok;
}

template <typename T, size_t Dimension>
bool operator!=(const Data<T, Dimension>& a, const Data<T, Dimension>& b) {
  return !(a == b);
}

#endif //GEOMETRY_GEOMETRY_DATA_H_
