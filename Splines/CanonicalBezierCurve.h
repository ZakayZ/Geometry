//
// Created by Artem Novikov on 04.07.2022.
//

#include <vector>
#include "BezierCurve.h"

#ifndef GEOMETRY_SPLINES_CANONICALBEZIERCURVE_H_
#define GEOMETRY_SPLINES_CANONICALBEZIERCURVE_H_

template <typename T, size_t Dimension, size_t Degree>
class CanonicalBezierCurve : public BezierCurve<T, Dimension, Degree> {
 public:
  /// constructors
  CanonicalBezierCurve(const CanonicalBezierCurve& other) = default;
  CanonicalBezierCurve(CanonicalBezierCurve&& other) noexcept = default;
  explicit CanonicalBezierCurve(const BezierCurve<T, Dimension, Degree>& curve, size_t cuts = 2);
  explicit CanonicalBezierCurve(BezierCurve<T, Dimension, Degree>&& curve, size_t cuts = 2) noexcept;
  template <typename... Args, typename = std::enable_if_t<sizeof...(Args) == Degree + 1>>
  CanonicalBezierCurve(size_t cuts, Args&& ... args); /// Remake
  CanonicalBezierCurve(std::initializer_list<Vector<T, Dimension>> list, size_t cuts = 2);
  template <typename U, template <typename, typename...> class Container, typename... Args>
  explicit CanonicalBezierCurve(const Container<Vector<U, Dimension>, Args...>& data, size_t cuts = 2);
  ~CanonicalBezierCurve() = default;
  CanonicalBezierCurve& operator=(const CanonicalBezierCurve& other) = default;
  CanonicalBezierCurve& operator=(CanonicalBezierCurve&& other) noexcept = default;

  /// getters and setter
  Point<T, Dimension> GetPoint(const T& length) const;
  Vector<T, Dimension> GetVelocity(const T& length) const;
  Vector<T, Dimension> GetAcceleration(const T& length) const;
  Point<T, Dimension> GetPoint(const T& length, Casteljau) const;
  T GetCurvature(const T& length) const;
  Vector<T, Dimension> GetNormal(const T& length) const;

  /// calc
  void UpdateCuts();
  inline T GetLength(const T& value = 1) const;
  inline T GetArcLength(const T& length) const;
  CanonicalBezierCurve<T, Dimension, Degree - 1> GetDerivative() const;
  template <typename... Args>
  std::vector<Point<T, Dimension>> GetDivision(size_t divisions, Args... args) const;
  void ResizeCuts(size_t cuts_size);

 private:
  std::vector<T> cuts;
};

template <typename T, size_t Dimension>
using CanonicalQuadraticBezier = CanonicalBezierCurve<T, Dimension, 2>;

using CanonicalQuadraticBezier2f = CanonicalQuadraticBezier<float, 2>;
using CanonicalQuadraticBezier2d = CanonicalQuadraticBezier<double, 2>;
using CanonicalQuadraticBezier3f = CanonicalQuadraticBezier<float, 2>;
using CanonicalQuadraticBezier3d = CanonicalQuadraticBezier<double, 2>;

template <typename T, size_t Dimension>
using CanonicalCubicBezier = CanonicalBezierCurve<T, Dimension, 3>;

using CanonicalCubicBezier2f = CanonicalCubicBezier<float, 2>;
using CanonicalCubicBezier2d = CanonicalCubicBezier<double, 2>;
using CanonicalCubicBezier3f = CanonicalCubicBezier<float, 2>;
using CanonicalCubicBezier3d = CanonicalCubicBezier<double, 2>;

template <size_t Divisions, typename T, size_t Dimension, size_t Degree, typename... Args>
std::array<Vector<T, Dimension>, Divisions> GetDivision(
    const CanonicalBezierCurve<T, Dimension, Degree>& curve, Args... args);

//////////////////////////////////////////////////////////DEFINITION////////////////////////////////////////////////////

template <typename T, size_t Dimension, size_t Degree>
CanonicalBezierCurve<T, Dimension, Degree>::CanonicalBezierCurve(
    const BezierCurve<T, Dimension, Degree>& curve, size_t cuts): BezierCurve<T, Dimension, Degree>(curve) {
  UpdateCuts(cuts);
}

template <typename T, size_t Dimension, size_t Degree>
CanonicalBezierCurve<T, Dimension, Degree>::CanonicalBezierCurve(BezierCurve<T, Dimension, Degree>&& curve, size_t cuts)
noexcept: BezierCurve<T, Dimension, Degree>(std::move(curve)) {
  UpdateCuts(cuts);
}

template <typename T, size_t Dimension, size_t Degree>
template <typename... Args, typename>
CanonicalBezierCurve<T, Dimension, Degree>::CanonicalBezierCurve(size_t cuts, Args&& ... args)
    : BezierCurve<T, Dimension, Degree>(std::forward<Args>(args)...) {
  UpdateCuts(cuts);
}

template <typename T, size_t Dimension, size_t Degree>
CanonicalBezierCurve<T, Dimension, Degree>::CanonicalBezierCurve(
    std::initializer_list<Vector<T, Dimension>> list, size_t cuts): BezierCurve<T, Dimension, Degree>(list) {
  UpdateCuts(cuts);
}

template <typename T, size_t Dimension, size_t Degree>
template <typename U, template <typename, typename...> class Container, typename... Args>
CanonicalBezierCurve<T, Dimension, Degree>::CanonicalBezierCurve(
    const Container<Vector<U, Dimension>, Args...>& data, size_t cuts):BezierCurve<T, Dimension, Degree>(data) {
  UpdateCuts(cuts);
}

template <typename T, size_t Dimension, size_t Degree>
Point<T, Dimension> CanonicalBezierCurve<T, Dimension, Degree>::GetPoint(const T& length) const {
  return BezierCurve<T, Dimension, Degree>::GetPoint(GetArcLength(length));
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> CanonicalBezierCurve<T, Dimension, Degree>::GetVelocity(const T& length) const {
  return Normalised(BezierCurve<T, Dimension, Degree>::GetVelocity(GetArcLength(length)));
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> CanonicalBezierCurve<T, Dimension, Degree>::GetAcceleration(const T& length) const {
  T t = GetArcLength(length);
  auto velocity = BezierCurve<T, Dimension, Degree>::GetVelocity(t);
  auto acceleration = BezierCurve<T, Dimension, Degree>::GetAcceleration(t);
  auto squared_length = velocity.SquaredLength();
  return (acceleration - velocity * (DotProduct(acceleration, velocity) / squared_length)) / squared_length;
}

template <typename T, size_t Dimension, size_t Degree>
Point<T, Dimension> CanonicalBezierCurve<T, Dimension, Degree>::GetPoint(const T& length, Casteljau) const {
  return BezierCurve<T, Dimension, Degree>::GetPoint(GetArcLength(length), Casteljau());
}

template <typename T, size_t Dimension, size_t Degree>
T CanonicalBezierCurve<T, Dimension, Degree>::GetCurvature(const T& length) const {
  return BezierCurve<T, Dimension, Degree>::GetCurvature(GetArcLength(length));
}

template <typename T, size_t Dimension, size_t Degree>
Vector<T, Dimension> CanonicalBezierCurve<T, Dimension, Degree>::GetNormal(const T& length) const {
  return BezierCurve<T, Dimension, Degree>::GetNormal(GetArcLength(length));
}

template <typename T, size_t Dimension, size_t Degree>
void CanonicalBezierCurve<T, Dimension, Degree>::UpdateCuts() {
  T delta = T(1) / (cuts.size() - 1);
  T value = 0;
  Vector<T, Dimension> previous_point = BezierCurve<T, Dimension, Degree>::GetPoint(value);
  for (size_t i = 1; i < cuts.size(); ++i) {
    value += delta;
    auto current_point = BezierCurve<T, Dimension, Degree>::GetPoint(value);
    cuts[i] = cuts[i - 1] + (current_point - previous_point).Length();
    previous_point = std::move(current_point);
  }
}

template <typename T, size_t Dimension, size_t Degree>
T CanonicalBezierCurve<T, Dimension, Degree>::GetLength(const T& value) const {
  T delta = T(1) / (cuts.size() - 1);
  auto index = static_cast<size_t>((cuts.size() - 1) * value);
  return cuts[index] + (cuts[index + 1] - cuts[index]) * (value - index * delta) / delta;
}

template <typename T, size_t Dimension, size_t Degree>
T CanonicalBezierCurve<T, Dimension, Degree>::GetArcLength(const T& length) const {
  T delta = T(1) / (cuts.size() - 1);
  auto lower = std::lower_bound(cuts.begin(), cuts.end(), length);
  if (lower == cuts.end()) { lower -= 2; }
  return delta * (static_cast<T>(lower - cuts.begin()) + (length - *lower) / (*(lower + 1) - *lower));
}

template <typename T, size_t Dimension, size_t Degree>
CanonicalBezierCurve<T, Dimension, Degree - 1> CanonicalBezierCurve<T, Dimension, Degree>::GetDerivative() const {
  return {BezierCurve<T, Dimension, Degree>::GetDerivative(), cuts.size()};
}

template <typename T, size_t Dimension, size_t Degree>
template <typename... Args>
std::vector<Point<T, Dimension>> CanonicalBezierCurve<T, Dimension, Degree>::GetDivision(
    size_t divisions, Args... args) const {
  std::vector<Point<T, Dimension>> division(divisions);
  T delta = cuts.back() / (divisions - 1);
  T value = 0;
  for (size_t i = 0; i < divisions; ++i) {
    division[i] = GetPoint(value, args...);
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
void CanonicalBezierCurve<T, Dimension, Degree>::ResizeCuts(size_t cuts_size) {
  cuts.resize(cuts_size);
  UpdateCuts();
}

template <size_t Divisions, typename T, size_t Dimension, size_t Degree, typename... Args>
std::array<Vector<T, Dimension>, Divisions> GetDivision(
    const CanonicalBezierCurve<T, Dimension, Degree>& curve, Args... args) {
  static_assert(Divisions > 1);
  std::array<Vector<T, Dimension>, Divisions> division;
  T delta = curve.cuts.back() / (Divisions - 1);
  T value = 0;
  for (size_t i = 0; i < Divisions; ++i) {
    division[i] = curve.GetPoint(value, args...);
    value += delta;
  }
  return division;
}

#endif //GEOMETRY_SPLINES_CANONICALBEZIERCURVE_H_
