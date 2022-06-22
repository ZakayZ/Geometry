//
// Created by Artem Novikov on 16.05.2022.
//

#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_LINE_H_
#define GEOMERTY_GEOMETRY_LINE_H_

template <typename T, size_t dim>
class Line {
 public:
  enum class Relationship {
    Parallel,
    Intersecting,
    Skew,
    Identical,
  };

  /// construction
  Line() = default;
  Line(const Point<T, dim>& p1, const Point<T, dim>& p2);
  Line(const Vector<T, dim>& p, const Vector<T, dim>& a);
  Line(const Line& other) = default;
  Line(Line&& other) noexcept = default;
  Line& operator=(const Line& other) = default;
  Line& operator=(Line&& other) noexcept = default;

  /// setters and getters
  Vector<T, dim> operator[](T t) const { return origin_ + direction_ * t; }
  const Vector<T, dim>& GetOrigin() const { return origin_; }
  const Vector<T, dim>& GetDirection() const { return direction_; }

 private:
  Vector<T, dim> origin_;
  Vector<T, dim> direction_;
};

/// Relationships

template <typename T, size_t dim>
typename Line<T, dim>::Relationship FindRelationship(const Line<T, dim>& a, const Line<T, dim>& b) {
  auto direction_relationship = FindRelationship(a.GetDirection(), b.GetDirection());
  if (direction_relationship == Vector<T, dim>::Relationship::Parallel) {
    return FindRelationship(a.GetDirection(), b.GetOrigin() - a.GetOrigin()) == Vector<T, dim>::Relationship::Parallel
           ? Line<T, dim>::Relationship::Identical
           : Line<T, dim>::Relationship::Parallel;
  }

  Vector<T, dim> dif = a.GetOrigin() - b.GetOrigin();
  size_t idx1 = dim;
  size_t idx2 = dim;
  for (size_t i = 0; i < dim; ++i) {
    if (a.GetDirection()[0]) {

    }
  }
}

#endif //GEOMERTY_GEOMETRY_LINE_H_
