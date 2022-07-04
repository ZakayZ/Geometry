//
// Created by Artem Novikov on 04.07.2022.
//

#include "Point.h"

#ifndef GEOMERTY_GEOMETRYENTITIES_BOUNDARYBOX_H_
#define GEOMERTY_GEOMETRYENTITIES_BOUNDARYBOX_H_

template <typename T, size_t Dimension>
class BoundaryBox {
 public:
  /// constructors
  BoundaryBox(const Point<T, Dimension>& point_l, const Point<T, Dimension>& point_r);
  BoundaryBox(const BoundaryBox& other) = default;
  BoundaryBox(BoundaryBox&& other) noexcept = default;
  ~BoundaryBox() = default;
  BoundaryBox& operator=(const BoundaryBox& other) = default;
  BoundaryBox& operator=(BoundaryBox&& other) noexcept = default;

  /// getters and setters
  Point<T, Dimension>& GetLeft() { return point_l_; }
  Point<T, Dimension>& GetRight() { return point_r_; }
  const Point<T, Dimension>& GetLeft() const { return point_l_; }
  const Point<T, Dimension>& GetRight() const { return point_r_; }

  /// calc
  bool Inside(const Point<T, Dimension>& point) const;
  bool Inside(const BoundaryBox& box) const;
  bool Intersects(const BoundaryBox& box) const;
 private:
  Point<T, Dimension> point_l_;
  Point<T, Dimension> point_r_;
};

template <typename T>
using BoundaryBox2 = BoundaryBox<T, 2>;

using BoundaryBox2i = BoundaryBox2<int>;
using BoundaryBox2f = BoundaryBox2<float>;
using BoundaryBox2d = BoundaryBox2<double>;

template <typename T>
using BoundaryBox3 = BoundaryBox<T, 3>;

using BoundaryBox3i = BoundaryBox3<int>;
using BoundaryBox3f = BoundaryBox3<float>;
using BoundaryBox3d = BoundaryBox3<double>;

//////////////////////////////////////////////////////DEFINITION////////////////////////////////////////////////////////

template <typename T, size_t Dimension>
BoundaryBox<T, Dimension>::BoundaryBox(const Point<T, Dimension>& point_l, const Point<T, Dimension>& point_r)
    : point_l_(point_l), point_r_(point_r) {}

template <typename T, size_t Dimension>
bool BoundaryBox<T, Dimension>::Inside(const Point<T, Dimension>& point) const {
  bool is_intersecting = true;
  for (size_t i = 0; i < Dimension; ++i) {
    is_intersecting &= GetLeft()[i] <= point[i] && point[i] <= GetRight()[i];
  }
  return is_intersecting;
}

template <typename T, size_t Dimension>
bool BoundaryBox<T, Dimension>::Inside(const BoundaryBox& box) const {
  bool is_inside = true;
  for (size_t i = 0; i < Dimension; ++i) {
    is_inside &= box.GetLeft()[i] <= GetLeft()[i] && GetRight()[i] <= box.GetRight()[i];
  }
  return is_inside;
}

template <typename T, size_t Dimension>
bool BoundaryBox<T, Dimension>::Intersects(const BoundaryBox& box) const {
  bool is_intersecting = true;
  for (size_t i = 0; i < Dimension; ++i) {
    is_intersecting &= std::max(box.GetLeft()[i], GetLeft()[i]) <= std::min(box.GetRight()[i], GetRight()[i]);
  }
  return is_intersecting;
}

#endif //GEOMERTY_GEOMETRYENTITIES_BOUNDARYBOX_H_
