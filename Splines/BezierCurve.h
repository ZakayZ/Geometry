//
// Created by Artem Novikov on 01.07.2022.
//

#include <array>
#include <vector>
#include "GeometryEntities/Vector.h"

#ifndef GEOMERTY_SPLINES_BEZIERCURVE_H_
#define GEOMERTY_SPLINES_BEZIERCURVE_H_

struct Casteljau {};

template <size_t N, size_t K>
struct BinomialCoefficient {
  static_assert(K <= N);
  static const size_t value = BinomialCoefficient<N - 1, K - 1>::value + BinomialCoefficient<N - 1, K>::value;
};

template <size_t K>
struct BinomialCoefficient<1, K> {
  static_assert(K <= 1);
  static const size_t value = 1;
};

template <>
struct BinomialCoefficient<0, 0> {
  static const size_t value = 1;
};

template <size_t N, size_t K>
static const size_t binomial_coefficient = BinomialCoefficient<N, K>::value;

template <typename T, size_t Dimension, size_t Power>
class BezierCurve {
 public:
  /// constructors
  template <typename... Args, typename = std::enable_if_t<sizeof...(Args) == Power>>
  BezierCurve(Args&& ... args);
  BezierCurve(std::initializer_list<Vector<T, Dimension>> list);
  template <typename U, template <typename, typename...> class Container, typename... Args>
  explicit BezierCurve(const Container<Vector<U, Dimension>, Args...>& data);
  BezierCurve(const BezierCurve& other) = default;
  BezierCurve(BezierCurve&& other) noexcept = default;
  template <typename U, template <typename, typename...> class Container, typename... Args>
  BezierCurve& operator=(const Container<Vector<U, Dimension>, Args...>& data);
  ~BezierCurve() = default;
  BezierCurve& operator=(const BezierCurve& other) = default;
  BezierCurve& operator=(BezierCurve&& other) noexcept = default;

  /// getters and setter
  Vector<T, Dimension>& GetPivot(size_t index) { return pivots_[index]; }
  const Vector<T, Dimension>& GetPivot(size_t index) const { return pivots_[index]; }
  Vector<T, Dimension> GetPoint(const T& value) const;
  Vector<T, Dimension> GetPoint(const T& value, Casteljau) const;
  Vector<T, Dimension> GetVelocity(const T& value) const;
  Vector<T, Dimension> GetAcceleration(const T& value) const;
  template <bool Access = Dimension <= 3, typename = std::enable_if_t<Access>>
  T GetCurvature(const T& value) const;
  template <bool Access = Dimension <= 3, typename = std::enable_if_t<Access>>
  Vector<T, Dimension> GetNormal(const T& value) const;

  std::vector<Vector<T, Dimension>> GetDivision(size_t divisions);
  std::vector<Vector<T, Dimension>> GetDivision(size_t divisions, Casteljau);

 private:
  template <typename U>
  void PushVector(U&& value);
  template <typename U, typename... Args>
  void PushVector(U&& value, Args&& ...args);
  static inline std::array<std::pair<T, T>, Power + 1> CalcPowers(const T& value);
  static inline Vector<T, Dimension> LinearInterpolation(
      const Vector<T, Dimension>& a, Vector<T, Dimension>& b, const T& value);

  template <size_t Index>
  inline Vector<T, Dimension> PointSum(const T& value, const std::array<std::pair<T, T>, Power + 1>& powers);

  template <size_t Index>
  inline Vector<T, Dimension> VelocitySum(const T& value, const std::array<std::pair<T, T>, Power + 1>& powers);

  template <size_t Index>
  inline Vector<T, Dimension> AccelerationSum(const T& value, const std::array<std::pair<T, T>, Power + 1>& powers);

  std::array<Vector<T, Dimension>, Power + 1> pivots_;
};

template <typename T, size_t Dimension>
using QuadraticBezier = BezierCurve<T, Dimension, 2>;

using QuadraticBezier2f = QuadraticBezier<float, 2>;
using QuadraticBezier2d = QuadraticBezier<double, 2>;
using QuadraticBezier3f = QuadraticBezier<float, 2>;
using QuadraticBezier3d = QuadraticBezier<double, 2>;

template <typename T, size_t Dimension>
using CubicBezier = BezierCurve<T, Dimension, 3>;

using CubicBezier2f = CubicBezier<float, 2>;
using CubicBezier2d = CubicBezier<double, 2>;
using CubicBezier3f = CubicBezier<float, 2>;
using CubicBezier3d = CubicBezier<double, 2>;

template <size_t Divisions, typename T, size_t Dimension, size_t Power>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Power>& curve);

template <size_t Divisions, typename T, size_t Dimension, size_t Power>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Power>& curve, Casteljau);

///////////////////////////////////////////////////DEFINITION///////////////////////////////////////////////////////////

template <typename T, size_t Dimension, size_t Power>
template <typename... Args, typename>
BezierCurve<T, Dimension, Power>::BezierCurve(Args&& ... args) {
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension, size_t Power>
BezierCurve<T, Dimension, Power>::BezierCurve(std::initializer_list<Vector<T, Dimension>> list) {
  assert(list.size() == Power + 1);
  size_t index = 0;
  for (auto& el : list) {
    pivots_[index++] = el;
  }
}

template <typename T, size_t Dimension, size_t Power>
template <typename U, template <typename, typename...> class Container, typename... Args>
BezierCurve<T, Dimension, Power>::BezierCurve(const Container<Vector<U, Dimension>, Args...>& data) {
  assert(data.size() == Power + 1);
  size_t index = 0;
  for (auto& el : data) {
    pivots_[index++] = el;
  }
}

template <typename T, size_t Dimension, size_t Power>
template <typename U, template <typename, typename...> class Container, typename... Args>
BezierCurve<T, Dimension, Power>& BezierCurve<T, Dimension, Power>::operator=(
    const Container<Vector<U, Dimension>, Args...>& data) {
  BezierCurve<T, Dimension, Power> copy = data;
  std::swap(copy.pivots_, pivots_);
  return *this;
}

template <typename T, size_t Dimension, size_t Power>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::GetPoint(const T& value) const {
  return PointSum<Power>(value);
}

template <typename T, size_t Dimension, size_t Power>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::GetVelocity(const T& value) const {
  return VelocitySum<Power>(value);
}

template <typename T, size_t Dimension, size_t Power>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::GetAcceleration(const T& value) const {
  return AccelerationSum<Power>(value);
}

template <typename T, size_t Dimension, size_t Power>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::GetPoint(const T& value, Casteljau) const {
  std::array<Vector<T, Dimension>, Power + 1> data = pivots_;
  for (size_t i = 0; i < Power; ++i) {
    for (size_t j = 0; j < Power - i; ++j) {
      data[j] = LinearInterpolation(data[j], data[j + 1], value);
    }
  }
  return data[0];
}

template <typename T, size_t Dimension, size_t Power>
std::vector<Vector<T, Dimension>> BezierCurve<T, Dimension, Power>::GetDivision(size_t divisions) {
  std::vector<Vector<T, Dimension>> division(divisions);
  T delta = T(1) / divisions;
  T value = 0;
  for (size_t i = 0; i < divisions; ++i) {
    division[i] = GetPoint(value);
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Power>
std::vector<Vector<T, Dimension>> BezierCurve<T, Dimension, Power>::GetDivision(size_t divisions, Casteljau) {
  std::vector<Vector<T, Dimension>> division(divisions);
  T delta = T(1) / divisions;
  T value = 0;
  for (size_t i = 0; i < divisions; ++i) {
    division[i] = GetPoint(value, Casteljau());
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Power>
template <bool Access, typename>
T BezierCurve<T, Dimension, Power>::GetCurvature(const T& value) const {
  if constexpr (Dimension == 2) {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return std::sqrt(acceleration.SquaredLength() / std::pow(1 + velocity.SquaredLength(), 3));
  } else {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return std::sqrt(CrossProduct(velocity, acceleration).SquaredLength() / std::pow(velocity.SquaredLength(), 3));
  }
}

template <typename T, size_t Dimension, size_t Power>
template <bool Access, typename>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::GetNormal(const T& value) const {
  if constexpr (Dimension == 2) {
    auto velocity = GetVelocity(value);
    return {-velocity.data().y, velocity.data().x};
  } else {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return acceleration - acceleration * (acceleration * velocity / acceleration.SquaredLength());
  }
}

template <typename T, size_t Dimension, size_t Power>
template <typename U>
void BezierCurve<T, Dimension, Power>::PushVector(U&& value) {
  pivots_[Power] = value;
}

template <typename T, size_t Dimension, size_t Power>
template <typename U, typename... Args>
void BezierCurve<T, Dimension, Power>::PushVector(U&& value, Args&& ... args) {
  pivots_[Power - sizeof...(Args)] = value;
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension, size_t Power>
std::array<std::pair<T, T>, Power + 1> BezierCurve<T, Dimension, Power>::CalcPowers(const T& value) {
  std::array<std::pair<T, T>, Power + 1> data;
  data[0] = {1, 1};
  for (size_t i = 1; i <= Power; ++i) {
    data[i].first = data[i - 1].first * (1 - value);
    data[i].second = data[i - 1].second * value;
  }
  return data;
}

template <typename T, size_t Dimension, size_t Power>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::LinearInterpolation(
    const Vector<T, Dimension>& a, Vector<T, Dimension>& b, const T& value) {
  return (1 - value) * a + value * b;
}

template <typename T, size_t Dimension, size_t Power>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::PointSum(
    const T& value, const std::array<std::pair<T, T>, Power + 1>& powers) {
  if constexpr (Index == 0) {
    return pivots_[Index] * powers[Power].first;
  } else { /// TODO find out what is better
//    auto sum = PointSum<Index - 1>(value);
//    sum += pivots_[Index]
//        * (std::pow(1 - value, Power - Index) * std::pow(value, Index) * binomial_coefficient<Power, Index>);
//    return sum;
    return PointSum<Index - 1>(value) + pivots_[Index]
        * (powers[Power - Index].first * powers[Index].second * binomial_coefficient<Power, Index>);
  }
}

template <typename T, size_t Dimension, size_t Power>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::VelocitySum(
    const T& value, const std::array<std::pair<T, T>, Power + 1>& powers) {
  if constexpr (Index == 0) {
    return pivots_[Index] * (-powers[Power - 1].first * Power);
  } else if constexpr (Index == Power) {
    return pivots_[Index] * (powers[Power - 1].second * Power);
  } else {
    return pivots_[Index] * (binomial_coefficient<Power, Index> * (
        Index * powers[Power - Index].first * powers[Index - 1].second -
            (Power - Index) * powers[Power - Index - 1].first * powers[Index].second)) +
        VelocitySum<Index - 1>(value);
  }
}

template <typename T, size_t Dimension, size_t Power>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Power>::AccelerationSum(
    const T& value, const std::array<std::pair<T, T>, Power + 1>& powers) {
  if constexpr (Index == 0) {
    return pivots_[Index] * (powers[Power - 2].first * Power * (Power - 1));
  } else if constexpr (Index == 1) {
    return pivots_[Index] * ((Power - 1) * (Power - 2) * value * powers[Power - 3].first
        - 2 * (Power - 1) * powers[Power - 2].first);
  } else if constexpr (Index == Power) {
    return pivots_[Index] * (Power * (Power - 1) * powers[Power - 2].second);
  } else if constexpr (Index + 1 == Power) {
    return pivots_[Index] * ((Power - 1) * (Power - 2) * (1 - value) * powers[Power - 3].second
        - 2 * (Power - 1) * powers[Power - 2].first);
  } else {
    return pivots_[Index] * (binomial_coefficient<Power, Index> * (
        (Power - Index) * (Power - Index - 1) * powers[Power - Index - 2].first * powers[Index].second -
            2 * (Power - Index) * Index * powers[Power - Index - 1].first * powers[Index - 1].second +
            Index * (Index - 1) * powers[Index - 2].second * powers[Power - Index].first)) +
        VelocitySum<Index - 1>(value);
  }
}

template <size_t Divisions, typename T, size_t Dimension, size_t Power>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Power>& curve) {
  std::array<Vector<T, Dimension>, Divisions> division;
  T delta = T(1) / Divisions;
  T value = 0;
  for (size_t i = 0; i < Divisions; ++i) {
    division[i] = curve.GetPoint(value);
    value += delta;
  }
  return division;
}

template <size_t Divisions, typename T, size_t Dimension, size_t Power>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Power>& curve, Casteljau) {
  std::array<Vector<T, Dimension>, Divisions> division;
  T delta = 1 / Divisions;
  T value = 0;
  for (size_t i = 0; i < Divisions; ++i) {
    division[i] = curve.GetPoint(value, Casteljau());
    value += delta;
  }
  return division;
}

#endif //GEOMERTY_SPLINES_BEZIERCURVE_H_
