//
// Created by Artem Novikov on 16.05.2022.
//

#include <cmath>

#ifndef GEOMERTY_GEOMETRY_COMPARATOR_H_
#define GEOMERTY_GEOMETRY_COMPARATOR_H_

template <typename T>
class Comparator {
 public:
  static bool Equal(const T& a, const T& b) {
    return a == b;
  }
 private:
};

template <>
class Comparator<double> {
 public:
  static bool Equal(double a, double b) {
    return std::abs(a) < sigma
           ? std::abs(a - b) < sigma
           : std::abs((a - b) / std::max(std::abs(a), std::abs(b))) < eps;
  }

 private:
  constexpr static const double sigma = 1e-4;
  constexpr static const double eps = 1e-4;
};

template <>
class Comparator<float> {
 public:
  static bool Equal(float a, float b) {
    return std::abs(a) < sigma
           ? std::abs(a - b) < sigma
           : std::abs((a - b) / std::max(std::abs(a), std::abs(b))) < eps;
  }

 private:
  constexpr static const float sigma = 1e-2;
  constexpr static const float eps = 1e-2;
};

#endif //GEOMERTY_GEOMETRY_COMPARATOR_H_