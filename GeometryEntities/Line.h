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

template <typename T, size_t Dimension>
class Line : public Void<T, Dimension> {
 public:
  /// construction
  Line() : Void<T, Dimension>(Entity::Line) {};
  Line(const Point<T, Dimension>& p1, const Point<T, Dimension>& p2);
  Line(const Point<T, Dimension>& point, const Vector<T, Dimension>& direction);
  Line(const Line& other) = default;
  Line(Line&& other) noexcept = default;
  Line& operator=(const Line& other) = default;
  Line& operator=(Line&& other) noexcept = default;

  /// setters and getters
  Point<T, Dimension> GetPoint(T t) const { return origin_ + direction_ * t; }
  Point<T, Dimension> operator[](T t) const { return origin_ + direction_ * t; }
  Point<T, Dimension>& GetOrigin() { return origin_; }
  Vector<T, Dimension>& GetDirection() { return direction_; }
  const Point<T, Dimension>& GetOrigin() const { return origin_; }
  const Vector<T, Dimension>& GetDirection() const { return direction_; }

  /// line calc
  void Normalise();
  bool Contains(const Point<T, Dimension>& point) const;
  T SquaredDistance(const Point<T, Dimension>& point) const;
  T Distance(const Point<T, Dimension>& point) const;
  T SquaredDistance(const Segment<T, Dimension>& segment) const;
  T Distance(const Segment<T, Dimension>& segment) const;
  T SquaredDistance(const Line<T, Dimension>& line) const;
  T Distance(const Line<T, Dimension>& line) const;

  std::unique_ptr<Void<T, Dimension>> Intersection(const Point<T, Dimension>& point) const;
  std::unique_ptr<Void<T, Dimension>> Intersection(const Segment<T, Dimension>& segment) const;
  std::unique_ptr<Void<T, Dimension>> Intersection(const Line<T, Dimension>& line) const;

  template <bool Temp = Dimension == 2, typename = std::enable_if_t<Temp>>
  Point2<T> Projection(const Point2<T>& point, const Vector2<T>& a) const;
  template <bool Temp = Dimension == 2, typename = std::enable_if_t<Temp>>
  Segment2<T> Projection(const Segment2<T>& segment, const Vector2<T>& a) const;

  Point<T, Dimension> Projection(const Point<T, Dimension>& point) const;
  Segment<T, Dimension> Projection(const Segment<T, Dimension>& segment) const;

 private:
  Point<T, Dimension> origin_;
  Vector<T, Dimension> direction_;
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

template <typename T, size_t Dimension>
bool operator==(const Line<T, Dimension>& a, const Line<T, Dimension>& b);

template <typename T, size_t Dimension>
bool operator!=(const Line<T, Dimension>& a, const Line<T, Dimension>& b);

template <typename T, size_t Dimension>
LineRelationship FindRelationship(const Line<T, Dimension>& a, const Line<T, Dimension>& b);

/////////////////////////////////////////////////////DEFINITION/////////////////////////////////////////////////////////

template <typename T, size_t Dimension>
Line<T, Dimension>::Line(const Point<T, Dimension>& p1, const Point<T, Dimension>& p2)
    : Void<T, Dimension>(Entity::Line), origin_(p1), direction_(p2 - p1) {}

template <typename T, size_t Dimension>
Line<T, Dimension>::Line(const Point<T, Dimension>& point, const Vector<T, Dimension>& direction)
    : Void<T, Dimension>(Entity::Line), origin_(point), direction_(direction) {}

template <typename T, size_t Dimension>
void Line<T, Dimension>::Normalise() {
  direction_.Normalise();
}

template <typename T, size_t Dimension>
bool Line<T, Dimension>::Contains(const Point<T, Dimension>& point) const {
  return FindRelationship(point - origin_, direction_) == VectorRelationship::Parallel;
}

template <typename T, size_t Dimension>
T Line<T, Dimension>::SquaredDistance(const Point<T, Dimension>& point) const {
  T t = (point - origin_) * direction_ / direction_.SquaredLength();
  return point.SquaredDistance(GetPoint(t));
}

template <typename T, size_t Dimension>
T Line<T, Dimension>::Distance(const Point<T, Dimension>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t Dimension>
T Line<T, Dimension>::SquaredDistance(const Segment<T, Dimension>& segment) const {
  Vector<T, Dimension> r = segment.GetLeft() - origin_;
  Vector<T, Dimension> a = segment.GetRight() - segment.GetLeft();
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

template <typename T, size_t Dimension>
T Line<T, Dimension>::Distance(const Segment<T, Dimension>& segment) const {
  return std::sqrt(SquaredDistance(segment));
}

template <typename T, size_t Dimension>
T Line<T, Dimension>::SquaredDistance(const Line<T, Dimension>& line) const {
  if constexpr(Dimension == 3) {
    T cross = CrossProduct(line.direction_, direction_).SquaredLength();
    if (Comparator<T>::IsZero(cross)) {
      return CrossProduct(line.origin_ - origin_, direction_).SquaredLength() / direction_.SquaredLength();
    }
    return std::pow(MixedProduct(line.origin_ - origin_, line.direction_, direction_), 2) / cross;
  } else {
    T dot = direction_ * line.direction_;
    Vector<T, Dimension> delta = origin_ - line.origin_;
    T t = (dot * delta * direction_ - delta * line.direction_ * direction_.SquaredLength()) /
        (direction_.SquaredLength() * line.direction_.SquaredLength() - dot * dot);
    return line.SquaredDistance(GetPoint(t));
  }
}

template <typename T, size_t Dimension>
T Line<T, Dimension>::Distance(const Line<T, Dimension>& line) const {
  return std::sqrt(SquaredDistance(line));
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Line<T, Dimension>::Intersection(const Point<T, Dimension>& point) const {
  if (Contains(point)) {
    return std::make_unique<Point<T, Dimension>>(point);
  }
  return std::make_unique<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Line<T, Dimension>::Intersection(const Segment<T, Dimension>& segment) const {
  Vector<T, Dimension> delta_r = segment.GetLeft() - origin_;
  Vector<T, Dimension> segment_dir = segment.GetRight() - segment.GetLeft();
  T dot = segment_dir * direction_;
  T len1 = direction_.SquaredLength();
  T len2 = segment_dir.SquaredLength();
  if (!Comparator<T>::Equal(len1 * len2, dot * dot)) {
    T t = (dot * direction_ * delta_r - len1 * segment_dir * delta_r) / (len1 * len2 - dot * dot);
    if (t >= 0 && t <= 1) {
      return std::make_unique<Point<T, Dimension>>(segment.GetPoint(t));
    }
  }
  if (Contains(segment.GetLeft()) && Contains(segment.GetRight())) {
    return std::make_unique<Segment<T, Dimension>>(segment);
  }
  return std::make_unique<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Line<T, Dimension>::Intersection(const Line<T, Dimension>& line) const {
  if (FindRelationship(direction_, line.direction_) == VectorRelationship::Parallel) {
    if (Contains(line.origin_)) {
      return std::make_unique<Line<T, Dimension>>(line);
    }
    return std::make_unique<Void<T, Dimension>>();
  }
  T dot = direction_ * line.direction_;
  Vector<T, Dimension> delta = origin_ - line.origin_;
  T t = (dot * delta * direction_ - delta * line.direction_ * direction_.SquaredLength()) /
      (direction_.SquaredLength() * line.direction_.SquaredLength() - dot * dot);
  if (Comparator<T>::IsZero(line.Distance(GetPoint(t)))) {
    return std::make_unique<Point<T, Dimension>>(GetPoint(t));
  }
  return std::make_unique<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
template <bool, typename>
Point2<T> Line<T, Dimension>::Projection(const Point2<T>& point, const Vector2<T>& a) const {
  Vector2<T> delta_r = point - origin_;
  T t = (delta_r[0] * a[1] - delta_r[1] * a[0]) / (direction_[0] * a[1] - direction_[1] * a[0]);
  return GetPoint(t);
}

template <typename T, size_t Dimension>
template <bool, typename>
Segment2<T> Line<T, Dimension>::Projection(const Segment2<T>& segment, const Vector2<T>& a) const {
  return {Projection(segment.GetLeft(), a), Projection(segment.GetRight(), a)};
}

template <typename T, size_t Dimension>
Point<T, Dimension> Line<T, Dimension>::Projection(const Point<T, Dimension>& point) const {
  T t = (point - origin_) * direction_ / direction_.SquaredLength();
  return GetPoint(t);
}

template <typename T, size_t Dimension>
Segment<T, Dimension> Line<T, Dimension>::Projection(const Segment<T, Dimension>& segment) const {
  return {Projection(segment.GetLeft()), Projection(segment.GetRight())};
}

template <typename T, size_t Dimension>
bool operator==(const Line<T, Dimension>& a, const Line<T, Dimension>& b) {
  return b.Contains(a.GetOrigin())
      && FindRelationship(a.GetDirection(), b.GetDirection()) == VectorRelationship::Parallel;
}

template <typename T, size_t Dimension>
bool operator!=(const Line<T, Dimension>& a, const Line<T, Dimension>& b) {
  return !(a == b);
}

template <typename T, size_t Dimension>
LineRelationship FindRelationship(const Line<T, Dimension>& a, const Line<T, Dimension>& b) {
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
