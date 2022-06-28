//
// Created by Artem Novikov on 28.06.2022.
//

#include "Vector.h"
#include "Point.h"

#ifndef GEOMERTY_GEOMETRY_SEGMENT_H_
#define GEOMERTY_GEOMETRY_SEGMENT_H_

template <typename T, size_t dim>
class Segment {
 public:
  enum class Relationship {
    Parallel,
    Intersecting,
    None,
    Identical,
  };

  /// construction
  Segment() = default;
  Segment(const Point<T, dim>& p1, const Point<T, dim>& p2);
  Segment(const Segment& other) = default;
  Segment(Segment&& other) noexcept = default;
  Segment& operator=(const Segment& other) = default;
  Segment& operator=(Segment&& other) noexcept = default;

  /// setters and getters
  Point<T, dim> GetPoint(T t) const { return point_l_ * t + point_r_ * (1 - t); }
  Point<T, dim> operator[](T t) const { return GetPoint(t); }
  Point<T, dim>& GetLeft() { return point_l_; }
  const Point<T, dim>& GetLeft() const { return point_l_; }
  Point<T, dim>& GetRight() { return point_r_; }
  const Point<T, dim>& GetRight() const { return point_r_; }

  /// segment calc
  bool Contains(const Point<T, dim>& point) const;
  T SquaredDistance(const Point<T, dim>& point) const;
  T Distance(const Point<T, dim>& point) const;
  T SquaredDistance(const Segment<T, dim>& segment) const;
  T Distance(const Segment<T, dim>& segment) const;

 private:
  Point<T, dim> point_l_;
  Point<T, dim> point_r_;
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

template <typename T, size_t dim>
bool operator==(const Segment<T, dim>& a, const Segment<T, dim>& b);

template <typename T, size_t dim>
bool operator!=(const Segment<T, dim>& a, const Segment<T, dim>& b);

template <typename T, size_t dim>
typename Segment<T, dim>::Relationship FindRelationship(const Segment<T, dim>& a, const Segment<T, dim>& b);

////////////////////////////////////////////////////DEFINITION//////////////////////////////////////////////////////////

template <typename T, size_t dim>
Segment<T, dim>::Segment(const Point<T, dim>& p1, const Point<T, dim>& p2) : point_l_(p1), point_r_(p2) {}

template <typename T, size_t dim>
bool Segment<T, dim>::Contains(const Point<T, dim>& point) const {
  return (point - point_l_) * (point_r_ - point_l_) >= 0 &&
      (point - point_r_) * (point_l_ - point_r_) >= 0 &&
      FindRelationship(point - point_l_, point_r_ - point_l_) == Vector<T, dim>::Relationship::Parallel;
}

template <typename T, size_t dim>
T Segment<T, dim>::SquaredDistance(const Point<T, dim>& point) const {
  if ((point - point_l_) * (point_r_ - point_l_) >= 0 &&
      (point - point_r_) * (point_l_ - point_r_) >= 0) {
    Vector<T, dim> dir = point_r_ - point_l_;
    T t = (point - point_l_) * dir / dir.SquaredLength();
    return point.SquaredDistance(GetPoint(t));
  }
  return std::min(point.SquaredDistance(point_l_), point.SquaredDistance(point_r_));
}

template <typename T, size_t dim>
T Segment<T, dim>::Distance(const Point<T, dim>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T, size_t dim>
T Segment<T, dim>::SquaredDistance(const Segment<T, dim>& segment) const {
  Vector<T, dim> r = segment.point_l_ - point_l_;
  Vector<T, dim> a1 = point_r_ - point_l_;
  Vector<T, dim> a2 = segment.point_r_ - segment.point_l_;
  T dot = a1 * a2;
  T len1 = a1.SquaredLength();
  T len2 = a2.SquaredLength();
  if (!Comparator<T>::Equal(len1 * len2, dot * dot)) {
    T t1 = (dot * a2 * r - len2 * a1 * r) / (len1 * len2 - dot * dot);
    T t2 = -(dot * a1 * r - len1 * a2 * r) / (len1 * len2 - dot * dot);
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
      return GetPoint(t1).SquaredDistance(segment.GetPoint(t2));
    }
  }
  return std::min(std::min(SquaredDistance(segment.point_l_), SquaredDistance(segment.point_r_)),
                  std::min(segment.SquaredDistance(point_l_), segment.SquaredDistance(point_r_)));
}

template <typename T, size_t dim>
T Segment<T, dim>::Distance(const Segment<T, dim>& segment) const {
  return std::sqrt(SquaredDistance(segment));
}

template <typename T, size_t dim>
bool operator==(const Segment<T, dim>& a, const Segment<T, dim>& b) {
  return a.GetLeft() == b.GetLeft() && a.GetRight() == b.GetRight() ||
      a.GetRight() == b.GetLeft() && a.GetLeft() == b.GetRight();
}

template <typename T, size_t dim>
bool operator!=(const Segment<T, dim>& a, const Segment<T, dim>& b) {
  return !(a == b);
}

template <typename T, size_t dim>
typename Segment<T, dim>::Relationship FindRelationship(const Segment<T, dim>& a, const Segment<T, dim>& b) {
  /// TODO Need matrix
}

#endif //GEOMERTY_GEOMETRY_SEGMENT_H_
