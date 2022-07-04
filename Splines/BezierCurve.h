//
// Created by Artem Novikov on 01.07.2022.
//

#include <array>
#include <vector>
#include "GeometryEntities/Vector.h"
#include "GeometryEntities/BoundaryBox.h"

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

template <typename T, size_t Dimension, size_t Degree>
class BezierCurve {
 public:
  /// constructors
  template <typename... Args, typename = std::enable_if_t<sizeof...(Args) == Degree>>
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
  Vector<T, Dimension>& GetControl(size_t index) { return controls_[index]; }
  const Vector<T, Dimension>& GetControl(size_t index) const { return controls_[index]; }
  Vector<T, Dimension> GetPoint(const T& value) const;
  Vector<T, Dimension> GetPoint(const T& value, Casteljau) const;
  Vector<T, Dimension> GetVelocity(const T& value) const;
  Vector<T, Dimension> GetAcceleration(const T& value) const;
  template <bool Access = Dimension <= 3, typename = std::enable_if_t<Access>>
  T GetCurvature(const T& value) const;
  template <bool Access = Dimension <= 3, typename = std::enable_if_t<Access>>
  Vector<T, Dimension> GetNormal(const T& value) const;

  /// calc
  std::vector<Vector<T, Dimension>> GetDivision(size_t divisions);
  std::vector<Vector<T, Dimension>> GetDivision(size_t divisions, Casteljau);
  void Shift(const Vector<T, Dimension>& shift) const;
  template <bool Access = Degree <= 3, typename = std::enable_if_t<Access>>
  BoundaryBox<T, Dimension> GetBoundaryBox() const;

 private:
  template <typename U>
  void PushVector(U&& value);
  template <typename U, typename... Args>
  void PushVector(U&& value, Args&& ...args);
  static inline std::array<std::pair<T, T>, Degree + 1> CalcDegrees(const T& value);
  static inline Vector<T, Dimension> LinearInterpolation(
      const Vector<T, Dimension>& a, Vector<T, Dimension>& b, const T& value);

  template <size_t Index>
  inline Vector<T, Dimension> PointSum(const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers);

  template <size_t Index>
  inline Vector<T, Dimension> VelocitySum(const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers);

  template <size_t Index>
  inline Vector<T, Dimension> AccelerationSum(const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers);

  std::array<Vector<T, Dimension>, Degree + 1> controls_;
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

template <size_t Divisions, typename T, size_t Dimension, size_t Degree>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Degree>& curve);

template <size_t Divisions, typename T, size_t Dimension, size_t Degree>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Degree>& curve, Casteljau);

template <typename T, size_t Dimension, size_t Degree>
BezierCurve<T, Dimension, Degree> Shifted(
    const BezierCurve<T, Dimension, Degree>& curve, const Vector<T, Dimension>& shift);

///////////////////////////////////////////////////DEFINITION///////////////////////////////////////////////////////////

template <typename T, size_t Dimension, size_t Degree>
template <typename... Args, typename>
BezierCurve<T, Dimension, Degree>::BezierCurve(Args&& ... args) {
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension, size_t Degree>
BezierCurve<T, Dimension, Degree>::BezierCurve(std::initializer_list<Vector<T, Dimension>> list) {
  assert(list.size() == Degree + 1);
  size_t index = 0;
  for (auto& el : list) {
    controls_[index++] = el;
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <typename U, template <typename, typename...> class Container, typename... Args>
BezierCurve<T, Dimension, Degree>::BezierCurve(const Container<Vector<U, Dimension>, Args...>& data) {
  assert(data.size() == Degree + 1);
  size_t index = 0;
  for (auto& el : data) {
    controls_[index++] = el;
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <typename U, template <typename, typename...> class Container, typename... Args>
BezierCurve<T, Dimension, Degree>& BezierCurve<T, Dimension, Degree>::operator=(
    const Container<Vector<U, Dimension>, Args...>& data) {
  BezierCurve<T, Dimension, Degree> copy = data;
  std::swap(copy.controls_, controls_);
  return *this;
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::GetPoint(const T& value) const {
  return PointSum<Degree>(value);
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::GetVelocity(const T& value) const {
  return VelocitySum<Degree>(value);
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::GetAcceleration(const T& value) const {
  return AccelerationSum<Degree>(value);
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::GetPoint(const T& value, Casteljau) const {
  std::array<Vector<T, Dimension>, Degree + 1> data = controls_;
  for (size_t i = 0; i < Degree; ++i) {
    for (size_t j = 0; j < Degree - i; ++j) {
      data[j] = LinearInterpolation(data[j], data[j + 1], value);
    }
  }
  return data[0];
}

template <typename T, size_t Dimension, size_t Degree>
template <bool Access, typename>
T BezierCurve<T, Dimension, Degree>::GetCurvature(const T& value) const {
  if constexpr (Dimension == 2) {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return std::sqrt(acceleration.SquaredLength() / std::pow(1 + velocity.SquaredLength(), 3));
  }
  if constexpr (Dimension == 3) {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return std::sqrt(CrossProduct(velocity, acceleration).SquaredLength() / std::pow(velocity.SquaredLength(), 3));
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <bool Access, typename>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::GetNormal(const T& value) const {
  if constexpr (Dimension == 2) {
    auto velocity = GetVelocity(value);
    return {-velocity.data().y, velocity.data().x};
  }
  if constexpr (Dimension == 3) {
    auto velocity = GetVelocity(value);
    auto acceleration = GetAcceleration(value);
    return acceleration - acceleration * (acceleration * velocity / acceleration.SquaredLength());
  }
}

template <typename T, size_t Dimension, size_t Degree>
std::vector<Vector<T, Dimension>> BezierCurve<T, Dimension, Degree>::GetDivision(size_t divisions) {
  std::vector<Vector<T, Dimension>> division(divisions);
  T delta = T(1) / (divisions - 1);
  T value = 0;
  for (size_t i = 0; i < divisions; ++i) {
    division[i] = GetPoint(value);
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
std::vector<Vector<T, Dimension>> BezierCurve<T, Dimension, Degree>::GetDivision(size_t divisions, Casteljau) {
  std::vector<Vector<T, Dimension>> division(divisions);
  T delta = T(1) / (divisions - 1);
  T value = 0;
  for (size_t i = 0; i < divisions; ++i) {
    division[i] = GetPoint(value, Casteljau());
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
void BezierCurve<T, Dimension, Degree>::Shift(const Vector<T, Dimension>& shift) const {
  for (auto& pivot : controls_) {
    pivot += shift;
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <bool Access, typename>
BoundaryBox<T, Dimension> BezierCurve<T, Dimension, Degree>::GetBoundaryBox() const {
  if constexpr (Degree == 2) {
    for (size_t i = 0; i < Dimension; ++i) {
      T t_dim = (controls_[0][i] - controls_[1][i]) / (controls_[0][i] + controls_[2][i] - 3 * controls_[1][i]);


    }
  }
  if constexpr (Degree == 3) {

  }
}

template <typename T, size_t Dimension, size_t Degree>
template <typename U>
void BezierCurve<T, Dimension, Degree>::PushVector(U&& value) {
  controls_[Degree] = value;
}

template <typename T, size_t Dimension, size_t Degree>
template <typename U, typename... Args>
void BezierCurve<T, Dimension, Degree>::PushVector(U&& value, Args&& ... args) {
  controls_[Degree - sizeof...(Args)] = value;
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension, size_t Degree>
std::array<std::pair<T, T>, Degree + 1> BezierCurve<T, Dimension, Degree>::CalcDegrees(const T& value) {
  std::array<std::pair<T, T>, Degree + 1> data;
  data[0] = {1, 1};
  for (size_t i = 1; i <= Degree; ++i) {
    data[i].first = data[i - 1].first * (1 - value);
    data[i].second = data[i - 1].second * value;
  }
  return data;
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::LinearInterpolation(
    const Vector<T, Dimension>& a, Vector<T, Dimension>& b, const T& value) {
  return (1 - value) * a + value * b;
}

template <typename T, size_t Dimension, size_t Degree>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::PointSum(
    const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers) {
  if constexpr (Index == 0) {
    return controls_[Index] * powers[Degree].first;
  } else { /// TODO find out what is better
//    auto sum = PointSum<Index - 1>(value);
//    sum += controls_[Index]
//        * (std::pow(1 - value, Degree - Index) * std::pow(value, Index) * binomial_coefficient<Degree, Index>);
//    return sum;
    return PointSum<Index - 1>(value) + controls_[Index]
        * (powers[Degree - Index].first * powers[Index].second * binomial_coefficient<Degree, Index>);
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::VelocitySum(
    const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers) {
  if constexpr (Index == 0) {
    return controls_[Index] * (-powers[Degree - 1].first * Degree);
  } else if constexpr (Index == Degree) {
    return controls_[Index] * (powers[Degree - 1].second * Degree);
  } else {
    return controls_[Index] * (binomial_coefficient<Degree, Index> * (
        Index * powers[Degree - Index].first * powers[Index - 1].second -
            (Degree - Index) * powers[Degree - Index - 1].first * powers[Index].second)) +
        VelocitySum<Index - 1>(value);
  }
}

template <typename T, size_t Dimension, size_t Degree>
template <size_t Index>
Vector<T, Dimension> BezierCurve<T, Dimension, Degree>::AccelerationSum(
    const T& value, const std::array<std::pair<T, T>, Degree + 1>& powers) {
  if constexpr (Index == 0) {
    return controls_[Index] * (powers[Degree - 2].first * Degree * (Degree - 1));
  } else if constexpr (Index == 1) {
    return controls_[Index] * ((Degree - 1) * (Degree - 2) * value * powers[Degree - 3].first
        - 2 * (Degree - 1) * powers[Degree - 2].first);
  } else if constexpr (Index == Degree) {
    return controls_[Index] * (Degree * (Degree - 1) * powers[Degree - 2].second);
  } else if constexpr (Index + 1 == Degree) {
    return controls_[Index] * ((Degree - 1) * (Degree - 2) * (1 - value) * powers[Degree - 3].second
        - 2 * (Degree - 1) * powers[Degree - 2].first);
  } else {
    return controls_[Index] * (binomial_coefficient<Degree, Index> * (
        (Degree - Index) * (Degree - Index - 1) * powers[Degree - Index - 2].first * powers[Index].second -
            2 * (Degree - Index) * Index * powers[Degree - Index - 1].first * powers[Index - 1].second +
            Index * (Index - 1) * powers[Index - 2].second * powers[Degree - Index].first)) +
        VelocitySum<Index - 1>(value);
  }
}

template <size_t Divisions, typename T, size_t Dimension, size_t Degree>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Degree>& curve) {
  static_assert(Divisions > 1);
  std::array<Vector<T, Dimension>, Divisions> division;
  T delta = T(1) / (Divisions - 1);
  T value = 0;
  for (size_t i = 0; i < Divisions; ++i) {
    division[i] = curve.GetPoint(value);
    value += delta;
  }
  return division;
}

template <size_t Divisions, typename T, size_t Dimension, size_t Degree>
std::array<Vector<T, Dimension>, Divisions> GetDivision(const BezierCurve<T, Dimension, Degree>& curve, Casteljau) {
  static_assert(Divisions > 1);
  std::array<Vector<T, Dimension>, Divisions> division;
  T delta = T(1) / (Divisions - 1);
  T value = 0;
  for (size_t i = 0; i < Divisions; ++i) {
    division[i] = curve.GetPoint(value, Casteljau());
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
BezierCurve<T, Dimension, Degree> Shifted(
    const BezierCurve<T, Dimension, Degree>& curve, const Vector<T, Dimension>& shift) {
  auto copy = curve;
  copy.Shift(shift);
  return copy;
}

#endif //GEOMERTY_SPLINES_BEZIERCURVE_H_
