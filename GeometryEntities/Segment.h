//
// Created by Artem Novikov on 28.06.2022.
//

#include "Void.h"
#include "Vector.h"
#include "BoundaryBox.h"
#include "Point.h"

#ifndef GEOMERTY_GEOMETRY_SEGMENT_H_
#define GEOMERTY_GEOMETRY_SEGMENT_H_

enum class SegmentRelationship {
  Parallel,
  Intersecting,
  Skew,
  Identical,
};

template <typename T, size_t Dimension>
class Segment : public Void<T, Dimension> {
 public:
  /// construction
  Segment() : Void<T, Dimension>(Entity::Segment) {};
  Segment(const Point<T, Dimension>& p1, const Point<T, Dimension>& p2);
  Segment(const Segment& other) = default;
  Segment(Segment&& other) noexcept = default;
  Segment& operator=(const Segment& other) = default;
  Segment& operator=(Segment&& other) noexcept = default;
  ~Segment() = default;

  /// setters and getters
  Point<T, Dimension> GetPoint(T t) const { return point_l_ * (1 - t) + point_r_ * t; }
  Point<T, Dimension> operator[](T t) const { return GetPoint(t); }
  Point<T, Dimension>& GetLeft() { return point_l_; }
  const Point<T, Dimension>& GetLeft() const { return point_l_; }
  Point<T, Dimension>& GetRight() { return point_r_; }
  const Point<T, Dimension>& GetRight() const { return point_r_; }
  Vector<T, Dimension> GetDirection() const { return point_r_ - point_l_; }
  Point<T, Dimension> GetMidpoint() const { return (point_l_ + point_r_) / 2; }

  /// segment calc
  bool Contains(const Point<T, Dimension>& point) const;

  T SquaredDistance(const Point<T, Dimension>& point) const;
  T SquaredDistance(const Segment<T, Dimension>& segment) const;

  T Distance(const Point<T, Dimension>& point) const;
  T Distance(const Segment<T, Dimension>& segment) const;

  std::unique_ptr<Void<T, Dimension>> Intersection(const Point<T, Dimension>& point) const;
  std::unique_ptr<Void<T, Dimension>> Intersection(const Segment<T, Dimension>& segment) const;

  Point<T, Dimension> Projection(const Point<T, Dimension>& point) const;
  Segment<T, Dimension> Projection(const Segment<T, Dimension>& segment) const;

  bool Intersects(const BoundaryBox<T, Dimension>& box) const;
  bool Inside(const BoundaryBox<T, Dimension>& box) const;

 private:
  inline std::pair<T, T> FindBoxIntersection(const BoundaryBox<T, Dimension>& box) const;

  Point<T, Dimension> point_l_;
  Point<T, Dimension> point_r_;
};

template <typename T>
using Segment1 = Segment<T, 1>;

using Segment1i = Segment1<int>;
using Segment1f = Segment1<float>;
using Segment1d = Segment1<double>;

template <typename T>
using Segment2 = Segment<T, 2>;

using Segment2i = Segment2<int>;
using Segment2f = Segment2<float>;
using Segment2d = Segment2<double>;

template <typename T>
using Segment3 = Segment<T, 3>;

using Segment3i = Segment3<int>;
using Segment3f = Segment3<float>;
using Segment3d = Segment3<double>;

/// Relationships

template <typename T, size_t Dimension>
bool operator==(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b);

template <typename T, size_t Dimension>
bool operator!=(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b);

template <typename T, size_t Dimension>
SegmentRelationship FindRelationship(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b);

////////////////////////////////////////////////////DEFINITION//////////////////////////////////////////////////////////

template <typename T, size_t Dimension>
Segment<T, Dimension>::Segment(const Point<T, Dimension>& p1, const Point<T, Dimension>& p2)
    : Void<T, Dimension>(Entity::Segment), point_l_(p1), point_r_(p2) {}

template <typename T, size_t Dimension>
bool Segment<T, Dimension>::Contains(const Point<T, Dimension>& point) const {
  return (point - point_l_) * GetDirection() >= 0 &&
      (point - point_r_) * GetDirection() <= 0 &&
      FindRelationship(point - point_l_, GetDirection()) == VectorRelationship::Parallel;
}

template <typename T, size_t Dimension>
T Segment<T, Dimension>::SquaredDistance(const Point<T, Dimension>& point) const {
  Vector<T, Dimension> dir = GetDirection();
  if ((point - point_l_) * dir >= 0 &&
      (point - point_r_) * dir <= 0) {
    T t = (point - point_l_) * dir / dir.SquaredLength();
    return point.SquaredDistance(GetPoint(t));
  }
  return std::min(point.SquaredDistance(point_l_), point.SquaredDistance(point_r_));
}

template <typename T, size_t Dimension>
T Segment<T, Dimension>::Distance(const Point<T, Dimension>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t Dimension>
T Segment<T, Dimension>::SquaredDistance(const Segment<T, Dimension>& segment) const {
  Vector<T, Dimension> r = segment.point_l_ - point_l_;
  Vector<T, Dimension> a1 = GetDirection();
  Vector<T, Dimension> a2 = segment.GetDirection();
  T dot = a1 * a2;
  T len1 = a1.SquaredLength();
  T len2 = a2.SquaredLength();
  if (!Comparator<T>::Equal(len1 * len2, dot * dot)) {
    T t1 = -(dot * a2 * r - len2 * a1 * r) / (len1 * len2 - dot * dot);
    T t2 = (dot * a1 * r - len1 * a2 * r) / (len1 * len2 - dot * dot);
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
      return GetPoint(t1).SquaredDistance(segment.GetPoint(t2));
    }
  }
  return std::min(std::min(SquaredDistance(segment.point_l_), SquaredDistance(segment.point_r_)),
                  std::min(segment.SquaredDistance(point_l_), segment.SquaredDistance(point_r_)));
}

template <typename T, size_t Dimension>
T Segment<T, Dimension>::Distance(const Segment<T, Dimension>& segment) const {
  return std::sqrt(SquaredDistance(segment));
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Segment<T, Dimension>::Intersection(const Point<T, Dimension>& point) const {
  if (Contains(point)) {
    return std::make_unique<Point<T, Dimension>>(point);
  }
  return std::make_unique<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
std::unique_ptr<Void<T, Dimension>> Segment<T, Dimension>::Intersection(const Segment<T, Dimension>& segment) const {
  Vector<T, Dimension> r = segment.point_l_ - point_l_;
  Vector<T, Dimension> a1 = GetDirection();
  Vector<T, Dimension> a2 = segment.GetDirection();
  T dot = a1 * a2;
  T len1 = a1.SquaredLength();
  T len2 = a2.SquaredLength();
  if (Comparator<T>::Equal(len1 * len2, dot * dot)) { /// REFACTOR:))))))
    if (Contains(segment.GetLeft())) {
      if (Contains(segment.GetRight())) {
        return std::make_unique<Segment<T, Dimension>>(segment);
      }
      if (segment.Contains(GetLeft())) {
        return std::make_unique<Segment<T, Dimension>>(segment.GetLeft(), GetLeft());
      }
      if (segment.Contains(GetRight())) {
        return std::make_unique<Segment<T, Dimension>>(segment.GetLeft(), GetRight());
      }
    }
    if (Contains(segment.GetRight())) {
      if (segment.Contains(GetLeft())) {
        return std::make_unique<Segment<T, Dimension>>(segment.GetRight(), GetLeft());
      }
      if (segment.Contains(GetRight())) {
        return std::make_unique<Segment<T, Dimension>>(segment.GetRight(), GetRight());
      }
    }
    if (segment.Contains(GetLeft()) && segment.Contains(GetRight())) {
      return std::make_unique<Segment<T, Dimension>>(GetLeft(), GetRight());
    }
  } else {
    T t1 = -(dot * a2 * r - len2 * a1 * r) / (len1 * len2 - dot * dot);
    T t2 = (dot * a1 * r - len1 * a2 * r) / (len1 * len2 - dot * dot);
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
        && Comparator<T>::IsZero(GetPoint(t1).Distance(segment.GetPoint(t2)))) {
      return std::make_unique<Point<T, Dimension>>(GetPoint(t1));
    }
  }
  return std::make_unique<Void<T, Dimension>>();
}

template <typename T, size_t Dimension>
Point<T, Dimension> Segment<T, Dimension>::Projection(const Point<T, Dimension>& point) const {
  auto dir = GetDirection();
  T t = DotProduct(point - point_l_, dir) / dir.SquaredLength();
  if (t <= 0) { return point_l_; }
  if (t >= 1) { return point_r_; }
  return GetPoint(t);
}

template <typename T, size_t Dimension>
Segment<T, Dimension> Segment<T, Dimension>::Projection(const Segment<T, Dimension>& segment) const {
  return {Projection(segment.point_l_), Projection(segment.point_r_)};
}

template <typename T, size_t Dimension>
bool operator==(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b) {
  return a.GetLeft() == b.GetLeft() && a.GetRight() == b.GetRight() ||
      a.GetRight() == b.GetLeft() && a.GetLeft() == b.GetRight();
}

template <typename T, size_t Dimension>
bool Segment<T, Dimension>::Intersects(const BoundaryBox<T, Dimension>& box) const {
  auto[t_min, t_max] = FindBoxIntersection(box);
  return t_min <= t_max && std::max(t_min, 0) <= std::min(t_max, 1);
}

template <typename T, size_t Dimension>
bool Segment<T, Dimension>::Inside(const BoundaryBox<T, Dimension>& box) const {
  auto[t_min, t_max] = FindBoxIntersection(box);
  return t_min <= t_max && 0 <= t_min && 1 <= t_max;
}

template <typename T, size_t Dimension>
std::pair<T, T> Segment<T, Dimension>::FindBoxIntersection(const BoundaryBox<T, Dimension>& box) const {
  T t_min;
  T t_max;
  auto direction = GetDirection();
  Vector<T, Dimension>& origin = point_l_;
  if ((box.GetRight()[0] - box.GetLeft()[0]) * direction[0] > 0) {
    t_min = (box.GetLeft()[0] - origin[0]) / direction[0];
    t_max = (box.GetRight()[0] - origin[0]) / direction[0];
  } else {
    t_min = (box.GetRight()[0] - origin[0]) / direction[0];
    t_max = (box.GetLeft()[0] - origin[0]) / direction[0];
  }
  for (size_t i = 1; i < Dimension; ++i) {
    if ((box.GetRight()[i] - box.GetLeft()[i]) * direction[i] > 0) {
      t_min = std::max(t_min, (box.GetLeft()[i] - origin[i]) / direction[i]);
      t_max = std::min(t_max, (box.GetRight()[i] - origin[i]) / direction[i]);
    } else {
      t_min = std::max(t_min, (box.GetRight()[i] - origin[i]) / direction[i]);
      t_max = std::min(t_max, (box.GetLeft()[i] - origin[i]) / direction[i]);
    }
  }
  return {t_min, t_max};
}

template <typename T, size_t Dimension>
bool operator!=(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b) {
  return !(a == b);
}

template <typename T, size_t Dimension>
SegmentRelationship FindRelationship(const Segment<T, Dimension>& a, const Segment<T, Dimension>& b) {
  if (FindRelationship(a.GetDirection(), b.GetDirection()) == VectorRelationship::Parallel) {
    if (FindRelationship(a.GetDirection(), b.GetLeft() - a.GetLeft()) == VectorRelationship::Parallel) {
      if (a == b) {
        return SegmentRelationship::Identical;
      }
      if (a.Contains(b.GetLeft()) || a.Contains(b.GetRight()) || b.Contains(a.GetLeft()) || b.Contains(a.GetRight())) {
        return SegmentRelationship::Intersecting;
      }
    }
    return SegmentRelationship::Parallel;;
  }
  if (Comparator<T>::IsZero(a.Distance(b))) {
    return SegmentRelationship::Intersecting;
  }
  return SegmentRelationship::Skew;
}

#endif //GEOMERTY_GEOMETRY_SEGMENT_H_
