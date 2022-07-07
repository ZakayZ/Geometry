//
// Created by Artem Novikov on 05.07.2022.
//

#include <vector>

#ifndef GEOMERTY_SPLINES_BEZIERHELPER_H_
#define GEOMERTY_SPLINES_BEZIERHELPER_H_

struct Casteljau {};

template <size_t N, size_t K>
struct BinomialCoefficient {
  static_assert(K <= N);
  static const size_t value = BinomialCoefficient<N - 1, K - 1>::value + BinomialCoefficient<N - 1, K>::value;
};

template <size_t N>
struct BinomialCoefficient<N, 0> {
  static const size_t value = 1;
};

template <size_t N>
struct BinomialCoefficient<N, N> {
  static const size_t value = 1;
};

template <size_t N, size_t K>
static const size_t binomial_coefficient = BinomialCoefficient<N, K>::value;

template <typename T>
class QuadraticSolver {
 public:
  inline static std::vector<T> Solve(const T& a, const T& b, const T& c) {
    std::vector<T> roots;
    roots.reserve(2);
    T determinant = b * b - 4 * a * c;
    if (determinant == 0) {
      roots.push_back(b / (2 * a));
    }
    if (determinant > 0) {
      roots.push_back((b - std::sqrt(determinant)) / (2 * a));
      roots.push_back((b + std::sqrt(determinant)) / (2 * a));
    }
    return roots;
  }
};

#endif //GEOMERTY_SPLINES_BEZIERHELPER_H_
