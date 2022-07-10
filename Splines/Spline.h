//
// Created by Artem Novikov on 08.07.2022.
//

#include <vector>
#include "BezierCurve.h"
#include "GeometryEntities/Transform.h"

#ifndef GEOMERTY_SPLINES_SPLINE_H_
#define GEOMERTY_SPLINES_SPLINE_H_

template <typename T, size_t Dimension, size_t Degree>
class Spline {
 public:
  using iterator = typename std::vector<BezierCurve<T, Dimension, Degree>>::iterator;
  using const_iterator = typename std::vector<BezierCurve<T, Dimension, Degree>>::const_iterator;

  Spline() = default;
  Spline(const Spline& other) = default;
  Spline(Spline&& other) noexcept = default;
  Spline(const std::vector<BezierCurve<T, Dimension, Degree>>& curves) : curves_(curves) {}
  Spline(std::vector<BezierCurve<T, Dimension, Degree>>&& curves) : curves_(std::move(curves)) {}
  Spline(std::initializer_list<BezierCurve<T, Dimension, Degree>> list) : curves_(list) {}
  ~Spline() = default;
  Spline& operator=(const Spline& other) = default;
  Spline& operator=(Spline&& other) noexcept = default;

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

  /// calc
  std::vector<Point<T, Dimension>> GetDivision(size_t divisions) const;
  template <size_t OutputDimension>
  Spline<T, OutputDimension, Degree> ApplyTransform(
      const Transform<T, Dimension, OutputDimension>& transform) const;

 private:
  std::vector<BezierCurve<T, Dimension, Degree>> curves_;
};

/////////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////

template <typename T, size_t Dimension, size_t Degree>
std::vector<Point<T, Dimension>> Spline<T, Dimension, Degree>::GetDivision(size_t divisions) const {
  std::vector<Point<T, Dimension>> division(divisions);
  T delta = T(curves_.size()) / (divisions - 1);
  T value = 0;
  size_t curve_index = 0;
  for (size_t i = 0; i < divisions; ++i) {
    if (value >= 1) {
      value -= 1;
      ++curve_index;
    }
    division[i] = curves_[curve_index].GetPoint(value);
    value += delta;
  }
  return division;
}

template <typename T, size_t Dimension, size_t Degree>
template <size_t OutputDimension>
Spline<T, OutputDimension, Degree> Spline<T, Dimension, Degree>::ApplyTransform(
    const Transform<T, Dimension, OutputDimension>& transform) const {
  std::vector<BezierCurve<T, OutputDimension, Degree>> new_spline;
  new_spline.reserve(curves_.size());
  for (const auto& curve : curves_) {
    new_spline.push_back(transform(curve));
  }
  return Spline<T, OutputDimension, Degree>(std::move(new_spline));
}

#endif //GEOMERTY_SPLINES_SPLINE_H_
