//
// Created by Artem Novikov on 16.05.2022.
//

#include "Vector.h"
#include "Point.h"

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
  Line(const Point<T, dim>& point, const Vector<T, dim>& direction);
  Line(const Line& other) = default;
  Line(Line&& other) noexcept = default;
  Line& operator=(const Line& other) = default;
  Line& operator=(Line&& other) noexcept = default;

  /// setters and getters
  Point<T, dim> GetPoint(T t) const { return origin_ + direction_ * t; }
  Point<T, dim> operator[](T t) const { return origin_ + direction_ * t; }
  Point<T, dim>& GetOrigin() { return origin_; }
  Vector<T, dim>& GetDirection() { return direction_; }
  const Point<T, dim>& GetOrigin() const { return origin_; }
  const Vector<T, dim>& GetDirection() const { return direction_; }

  /// line calc
  void Normalise();
  bool Contains(const Point<T, dim>& point) const;
  T SquaredDistance(const Point<T, dim>& point) const;
  T Distance(const Point<T, dim>& point) const;
  T SquaredDistance(const Line<T, dim>& line) const;
  T Distance(const Line<T, dim>& line) const;

 private:
  Point<T, dim> origin_;
  Vector<T, dim> direction_;
};

template <typename T>
using Line2 = Line<T, 2>;

using Line2i = Line2<int>;
using Line2f = Line2<float>;
using Line2d = Line2<double>;

template <typename T>
using Line3 = Line<T, 3>;

using Line3i = Line3<int>;
using Line3f = Line3<float>;
using Line3d = Line3<double>;

/// Relationships

template <typename T, size_t dim>
bool operator==(const Line<T, dim>& a, const Line<T, dim>& b);

template <typename T, size_t dim>
bool operator!=(const Line<T, dim>& a, const Line<T, dim>& b);

template <typename T, size_t dim>
typename Line<T, dim>::Relationship FindRelationship(const Line<T, dim>& a, const Line<T, dim>& b);

/////////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////

template <typename T, size_t dim>
Line<T, dim>::Line(const Point<T, dim>& p1, const Point<T, dim>& p2) : origin_(p1), direction_(p2 - p1) {}

template <typename T, size_t dim>
Line<T, dim>::Line(const Point<T, dim>& point, const Vector<T, dim>& direction)
    : origin_(point), direction_(direction) {}

template <typename T, size_t dim>
void Line<T, dim>::Normalise() {
  direction_.Normalise();
}

template <typename T, size_t dim>
bool Line<T, dim>::Contains(const Point<T, dim>& point) const {
  return FindRelationship(point - origin_, direction_) == Vector<T, dim>::Relationship::Parallel;
}

template <typename T, size_t dim>
T Line<T, dim>::SquaredDistance(const Point<T, dim>& point) const {
  T t = (point - origin_) * direction_ / direction_.SquaredLength();
  T squared_distance = point.SquaredDistance(GetPoint(t));
  return squared_distance;
}

template <typename T, size_t dim>
T Line<T, dim>::Distance(const Point<T, dim>& point) const {
  T distance = std::sqrt(SquaredDistance(point));
  return distance;
}

template <typename T, size_t dim>
T Line<T, dim>::SquaredDistance(const Line<T, dim>& line) const {
  if constexpr(dim == 3) {
    T cross = CrossProduct(line.direction_, direction_);
    if (cross == 0) {
      T squared_distance =
          CrossProduct(line.origin_ - origin_, direction_).SquaredLength() / direction_.SquaredLength();
      return squared_distance;
    }
    return std::pow(MixedProduct(line.origin_ - origin_, line.direction_, direction_) / cross, 2);
  } else {
    T sq = direction_ * line.direction_;
    Vector<T, dim> delta = origin_ - line.origin_;
    T t = (sq * delta * direction_ - delta * line.direction_ * direction_.SquaredLength()) /
        (direction_.SquaredLength() * line.direction_.SquaredLength() - sq * sq);
    return line.SquaredDistance(GetPoint(t));
  }
}

template <typename T, size_t dim>
T Line<T, dim>::Distance(const Line<T, dim>& line) const {
  if constexpr(dim == 3) {
    T cross = CrossProduct(line.direction_, direction_);
    if (cross == 0) {
      T distance =
          std::sqrt(CrossProduct(line.origin_ - origin_, direction_).SquaredLength() / direction_.SquaredLength());
      return distance;
    }
    return MixedProduct(line.origin_ - origin_, line.direction_, direction_) / cross;
  } else {
    T distance = std::sqrt(SquaredDistance(line));
    return distance;
  }
}

template <typename T, size_t dim>
bool operator==(const Line<T, dim>& a, const Line<T, dim>& b) {
  return b.Contains(a.GetOrigin())
      && FindRelationship(a.GetDirection(), b.GetDirection()) == Vector<T, dim>::Relationship::Paralel;
}

template <typename T, size_t dim>
bool operator!=(const Line<T, dim>& a, const Line<T, dim>& b) {
  return !(a == b);
}

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
  /// TODO Need matrix
}

#endif //GEOMERTY_GEOMETRY_LINE_H_
