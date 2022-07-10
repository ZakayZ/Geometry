//
// Created by Artem Novikov on 10.07.2022.
//

#include <vector>
#include "CanonicalBezierCurve.h"
#include "GeometryEntities/Transform.h"

#ifndef GEOMERTY_SPLINES_CANONICALSPLINE_H_
#define GEOMERTY_SPLINES_CANONICALSPLINE_H_

template <typename T, size_t Dimension, size_t Degree>
class CanonicalSpline {
 public:
  using iterator = typename std::vector<BezierCurve<T, Dimension, Degree>>::iterator;
  using const_iterator = typename std::vector<BezierCurve<T, Dimension, Degree>>::const_iterator;

  CanonicalSpline() = default;
  CanonicalSpline(const CanonicalSpline& other) = default;
  CanonicalSpline(CanonicalSpline&& other) noexcept = default;
  CanonicalSpline(const std::vector<CanonicalBezierCurve<T, Dimension, Degree>>& curves) : curves_(curves) {}
  CanonicalSpline(std::vector<CanonicalBezierCurve<T, Dimension, Degree>>&& curves) : curves_(std::move(curves)) {}
  CanonicalSpline(std::initializer_list<CanonicalBezierCurve<T, Dimension, Degree>> list) : curves_(list) {}
  ~CanonicalSpline() = default;
  CanonicalSpline& operator=(const CanonicalSpline& other) = default;
  CanonicalSpline& operator=(CanonicalSpline&& other) noexcept = default;

  /// getters and setters
  iterator begin() { return curves_.begin(); }
  const_iterator begin() const { return curves_.cbegin(); }
  const_iterator cbegin() { return curves_.cbegin(); }
  iterator end() { return curves_.end(); }
  const_iterator end() const { return curves_.cend(); }
  const_iterator cend() { return curves_.cend(); }
  BezierCurve<T, Dimension, Degree>& GetCurve(size_t index) { return curves_[index]; }
  const BezierCurve<T, Dimension, Degree>& GetCurve(size_t index) const { return curves_[index]; }
  size_t GetSize() const { return curves_.size(); }
  T GetLength() const;

  /// calc
  std::vector<Point<T, Dimension>> GetDivision(size_t divisions) const;
  template <size_t OutputDimension>
  CanonicalSpline<T, OutputDimension, Degree> ApplyTransform(
      const Transform<T, Dimension, OutputDimension>& transform) const;
 private:
  std::vector<CanonicalBezierCurve<T, Dimension, Degree>> curves_;
};

/////////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////

template <typename T, size_t Dimension, size_t Degree>
T CanonicalSpline<T, Dimension, Degree>::GetLength() const {
  T length = 0;
  for (const auto& curve : curves_) {
    length += curve.GetLength();
  }
  return length;
}

template <typename T, size_t Dimension, size_t Degree>
std::vector<Point<T, Dimension>> CanonicalSpline<T, Dimension, Degree>::GetDivision(size_t divisions) const {
  std::vector<Point<T, Dimension>> division(divisions);
  T delta = GetLength() / (divisions - 1);
  T value = 0;
  size_t curve_index = 0;
  for (size_t i = 0; i < divisions; ++i) {
    if (curves_[curve_index].GetLength() <= value) {
      value -= curves_[curve_index].GetLength();
      ++curve_index;
    }
    division[i] = curves_[curve_index].GetPoint(value);
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
template <size_t OutputDimension>
CanonicalSpline<T, OutputDimension, Degree> CanonicalSpline<T, Dimension, Degree>::ApplyTransform(
    const Transform<T, Dimension, OutputDimension>& transform) const {
  std::vector<CanonicalBezierCurve<T, OutputDimension, Degree>> new_spline;
  new_spline.reserve(curves_.size());
  for (const auto& curve : curves_) {
    new_spline.push_back(transform(curve));
  }
  return Spline<T, OutputDimension, Degree>(std::move(new_spline));
}

#endif //GEOMERTY_SPLINES_CANONICALSPLINE_H_
