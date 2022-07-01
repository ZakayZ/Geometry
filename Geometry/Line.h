//
// Created by Artem Novikov on 16.05.2022.
//

#include "Void.h"
#include "Vector.h"
#include "Point.h"
#include "Segment.h"

#ifndef GEOMERTY_GEOMETRY_LINE_H_
#define GEOMERTY_GEOMETRY_LINE_H_

enum class LineRelationship {
  Parallel,
  Intersecting,
  Skew,
  Identical,
};

template <typename T, size_t dim>
class Line : public Void<T, dim> {
 public:
  /// construction
  Line() : Void<T, dim>(Entity::Line) {};
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
  T SquaredDistance(const Segment<T, dim>& segment) const;
  T Distance(const Segment<T, dim>& segment) const;
  T SquaredDistance(const Line<T, dim>& line) const;
  T Distance(const Line<T, dim>& line) const;

  std::unique_ptr<Void<T, dim>> Intersection(const Point<T, dim>& point) const;
  std::unique_ptr<Void<T, dim>> Intersection(const Segment<T, dim>& segment) const;
  std::unique_ptr<Void<T, dim>> Intersection(const Line<T, dim>& line) const;

  template <bool Temp = dim == 2, typename = std::enable_if_t<Temp>>
  Point2<T> Projection(const Point2<T>& point, const Vector2<T>& a) const;
  template <bool Temp = dim == 2, typename = std::enable_if_t<Temp>>
  Segment2<T> Projection(const Segment2<T>& segment, const Vector2<T>& a) const;

  Point<T, dim> Projection(const Point<T, dim>& point) const;
  Segment<T, dim> Projection(const Segment<T, dim>& segment) const;

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
LineRelationship FindRelationship(const Line<T, dim>& a, const Line<T, dim>& b);

/////////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////

template <typename T, size_t dim>
Line<T, dim>::Line(const Point<T, dim>& p1, const Point<T, dim>& p2)
    : Void<T, dim>(Entity::Line), origin_(p1), direction_(p2 - p1) {}

template <typename T, size_t dim>
Line<T, dim>::Line(const Point<T, dim>& point, const Vector<T, dim>& direction)
    : Void<T, dim>(Entity::Line), origin_(point), direction_(direction) {}

template <typename T, size_t dim>
void Line<T, dim>::Normalise() {
  direction_.Normalise();
}

template <typename T, size_t dim>
bool Line<T, dim>::Contains(const Point<T, dim>& point) const {
  return FindRelationship(point - origin_, direction_) == VectorRelationship::Parallel;
}

template <typename T, size_t dim>
T Line<T, dim>::SquaredDistance(const Point<T, dim>& point) const {
  T t = (point - origin_) * direction_ / direction_.SquaredLength();
  return point.SquaredDistance(GetPoint(t));
}

template <typename T, size_t dim>
T Line<T, dim>::Distance(const Point<T, dim>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t dim>
T Line<T, dim>::SquaredDistance(const Segment<T, dim>& segment) const {
  Vector<T, dim> r = segment.GetLeft() - origin_;
  Vector<T, dim> a = segment.GetRight() - segment.GetLeft();
  T dot = a * direction_;
  T len1 = direction_.SquaredLength();
  T len2 = a.SquaredLength();
  if (!Comparator<T>::Equal(len1 * len2, dot * dot)) {
    T t = -(dot * direction_ * r - len1 * a * r) / (len1 * len2 - dot * dot);
    if (t >= 0 && t <= 1) {
      return SquaredDistance(segment.GetPoint(t));
    }
  }
  return std::min(SquaredDistance(segment.GetLeft()), SquaredDistance(segment.GetRight()));
}

template <typename T, size_t dim>
T Line<T, dim>::Distance(const Segment<T, dim>& segment) const {
  return std::sqrt(SquaredDistance(segment));
}

template <typename T, size_t dim>
T Line<T, dim>::SquaredDistance(const Line<T, dim>& line) const {
  if constexpr(dim == 3) {
    T cross = CrossProduct(line.direction_, direction_).SquaredLength();
    if (Comparator<T>::IsZero(cross)) {
      return CrossProduct(line.origin_ - origin_, direction_).SquaredLength() / direction_.SquaredLength();
    }
    return std::pow(MixedProduct(line.origin_ - origin_, line.direction_, direction_), 2) / cross;
  } else {
    T dot = direction_ * line.direction_;
    Vector<T, dim> delta = origin_ - line.origin_;
    T t = (dot * delta * direction_ - delta * line.direction_ * direction_.SquaredLength()) /
        (direction_.SquaredLength() * line.direction_.SquaredLength() - dot * dot);
    return line.SquaredDistance(GetPoint(t));
  }
}

template <typename T, size_t dim>
T Line<T, dim>::Distance(const Line<T, dim>& line) const {
  return std::sqrt(SquaredDistance(line));
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Line<T, dim>::Intersection(const Point<T, dim>& point) const {
  if (Contains(point)) {
    return std::make_unique<Point<T, dim>>(point);
  }
  return std::make_unique<Void<T, dim>>();
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Line<T, dim>::Intersection(const Segment<T, dim>& segment) const {
  Vector<T, dim> delta_r = segment.GetLeft() - origin_;
  Vector<T, dim> segment_dir = segment.GetRight() - segment.GetLeft();
  T dot = segment_dir * direction_;
  T len1 = direction_.SquaredLength();
  T len2 = segment_dir.SquaredLength();
  if (!Comparator<T>::Equal(len1 * len2, dot * dot)) {
    T t = (dot * direction_ * delta_r - len1 * segment_dir * delta_r) / (len1 * len2 - dot * dot);
    if (t >= 0 && t <= 1) {
      return std::make_unique<Point<T, dim>>(segment.GetPoint(t));
    }
  }
  if (Contains(segment.GetLeft()) && Contains(segment.GetRight())) {
    return std::make_unique<Segment<T, dim>>(segment);
  }
  return std::make_unique<Void<T, dim>>();
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Line<T, dim>::Intersection(const Line<T, dim>& line) const {
  if (FindRelationship(direction_, line.direction_) == VectorRelationship::Parallel) {
    if (Contains(line.origin_)) {
      return std::make_unique<Line<T, dim>>(line);
    }
    return std::make_unique<Void<T, dim>>();
  }
  T dot = direction_ * line.direction_;
  Vector<T, dim> delta = origin_ - line.origin_;
  T t = (dot * delta * direction_ - delta * line.direction_ * direction_.SquaredLength()) /
      (direction_.SquaredLength() * line.direction_.SquaredLength() - dot * dot);
  if (Comparator<T>::IsZero(line.Distance(GetPoint(t)))) {
    return std::make_unique<Point<T, dim>>(GetPoint(t));
  }
  return std::make_unique<Void<T, dim>>();
}

template <typename T, size_t dim>
template <bool, typename>
Point2<T> Line<T, dim>::Projection(const Point2<T>& point, const Vector2<T>& a) const {
  Vector2<T> delta_r = point - origin_;
  T t = (delta_r[0] * a[1] - delta_r[1] * a[0]) / (direction_[0] * a[1] - direction_[1] * a[0]);
  return GetPoint(t);
}

template <typename T, size_t dim>
template <bool, typename>
Segment2<T> Line<T, dim>::Projection(const Segment2<T>& segment, const Vector2<T>& a) const {
  return {Projection(segment.GetLeft(), a), Projection(segment.GetRight(), a)};
}

template <typename T, size_t dim>
Point<T, dim> Line<T, dim>::Projection(const Point<T, dim>& point) const {
  T t = (point - origin_) * direction_ / direction_.SquaredLength();
  return GetPoint(t);
}

template <typename T, size_t dim>
Segment<T, dim> Line<T, dim>::Projection(const Segment<T, dim>& segment) const {
  return {Projection(segment.GetLeft()), Projection(segment.GetRight())};
}

template <typename T, size_t dim>
bool operator==(const Line<T, dim>& a, const Line<T, dim>& b) {
  return b.Contains(a.GetOrigin())
      && FindRelationship(a.GetDirection(), b.GetDirection()) == VectorRelationship::Parallel;
}

template <typename T, size_t dim>
bool operator!=(const Line<T, dim>& a, const Line<T, dim>& b) {
  return !(a == b);
}

template <typename T, size_t dim>
LineRelationship FindRelationship(const Line<T, dim>& a, const Line<T, dim>& b) {
  if (FindRelationship(a.GetDirection(), b.GetDirection()) == VectorRelationship::Parallel) {
    return FindRelationship(a.GetDirection(), b.GetOrigin() - a.GetOrigin()) == VectorRelationship::Parallel
           ? LineRelationship::Identical
           : LineRelationship::Parallel;
  }

  if (Comparator<T>::IsZero(a.Distance(b))) {
    return LineRelationship::Intersecting;
  }
  return LineRelationship::Skew;
}

#endif //GEOMERTY_GEOMETRY_LINE_H_
