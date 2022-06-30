//
// Created by Artem Novikov on 28.06.2022.
//

#include "Void.h"
#include "Vector.h"
#include "Point.h"

#ifndef GEOMERTY_GEOMETRY_SEGMENT_H_
#define GEOMERTY_GEOMETRY_SEGMENT_H_

template <typename T, size_t dim>
class Segment : public Void<T, dim> {
 public:
  enum class Relationship {
    Parallel,
    Intersecting,
    Skew,
    Identical,
  };

  /// construction
  Segment() : Void<T, dim>(Entity::Segment) {};
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
  Vector<T, dim> GetDirection() const { return point_r_ - point_l_; }
  Point<T, dim> GetMidpoint() const { return (point_l_ + point_r_) / 2; }

  /// segment calc
  bool Contains(const Point<T, dim>& point) const;

  T SquaredDistance(const Point<T, dim>& point) const;
  T SquaredDistance(const Segment<T, dim>& segment) const;

  T Distance(const Point<T, dim>& point) const;
  T Distance(const Segment<T, dim>& segment) const;

  std::unique_ptr<Void<T, dim>> Intersection(const Point<T, dim>& point) const;
  std::unique_ptr<Void<T, dim>> Intersection(const Segment<T, dim>& segment) const;

  std::unique_ptr<Void<T, dim>> Projection(const Point3<T>& point) const;

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
Segment<T, dim>::Segment(const Point<T, dim>& p1, const Point<T, dim>& p2)
    : Void<T, dim>(Entity::Segment), point_l_(p1), point_r_(p2) {}

template <typename T, size_t dim>
bool Segment<T, dim>::Contains(const Point<T, dim>& point) const {
  return (point - point_l_) * GetDirection() >= 0 &&
      (point - point_r_) * GetDirection() <= 0 &&
      FindRelationship(point - point_l_, GetDirection()) == Vector<T, dim>::Relationship::Parallel;
}

template <typename T, size_t dim>
T Segment<T, dim>::SquaredDistance(const Point<T, dim>& point) const {
  Vector<T, dim> dir = GetDirection();
  if ((point - point_l_) * dir >= 0 &&
      (point - point_r_) * dir <= 0) {
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
  Vector<T, dim> a1 = GetDirection();
  Vector<T, dim> a2 = segment.GetDirection();
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
std::unique_ptr<Void<T, dim>> Segment<T, dim>::Intersection(const Point<T, dim>& point) const {
  if (Contains(point)) {
    return std::unique_ptr<Point<T, dim>>(point);
  }
  return std::unique_ptr<Void<T, dim>>();
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Segment<T, dim>::Intersection(const Segment<T, dim>& segment) const {
  Vector<T, dim> r = segment.point_l_ - point_l_;
  Vector<T, dim> a1 = GetDirection();
  Vector<T, dim> a2 = segment.GetDirection();
  T dot = a1 * a2;
  T len1 = a1.SquaredLength();
  T len2 = a2.SquaredLength();
  if (Comparator<T>::Equal(len1 * len2, dot * dot)) { /// REFACTOR:))))))
    if (Contains(segment.GetLeft())) {
      if (Contains(segment.GetRight())) {
        return std::unique_ptr<Segment<T, dim>>(segment);
      }
      if (segment.Contains(GetLeft())) {
        return std::unique_ptr<Segment<T, dim>>(segment.GetLeft(), GetLeft());
      }
      if (segment.Contains(GetRight())) {
        return std::unique_ptr<Segment<T, dim>>(segment.GetLeft(), GetRight());
      }
    }
    if (Contains(segment.GetRight())) {
      if (segment.Contains(GetLeft())) {
        return std::unique_ptr<Segment<T, dim>>(segment.GetRight(), GetLeft());
      }
      if (segment.Contains(GetRight())) {
        return std::unique_ptr<Segment<T, dim>>(segment.GetRight(), GetRight());
      }
    }
    if (segment.Contains(GetLeft()) && segment.Contains(GetRight())) {
      return std::unique_ptr<Segment<T, dim>>(GetLeft(), GetRight());
    }
  } else {
    T t1 = (dot * a2 * r - len2 * a1 * r) / (len1 * len2 - dot * dot);
    T t2 = -(dot * a1 * r - len1 * a2 * r) / (len1 * len2 - dot * dot);
    if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1
        && Comparator<T>::IsZero(GetPoint(t1).Distance(segment.GetPoint(t2)))) {
      return std::unique_ptr<Point<T, dim>>(GetPoint(t1));
    }
  }
  return std::unique_ptr<Void<T, dim>>();
}

template <typename T, size_t dim>
std::unique_ptr<Void<T, dim>> Segment<T, dim>::Projection(const Point3<T>& point) const {
  auto dir = GetDirection();
  T t = DotProduct(point - point_l_, dir) / dir.SquaredLength();
  if (t >= 0 && t <= 1) {
    return std::unique_ptr<Point<T, dim>>(GetPoint(t));
  }
  return std::unique_ptr<Void<T, dim>>();
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
  if (FindRelationship(a.GetDirection(), b.GetDirection()) == Vector<T, dim>::Relationship::Parallel) {
    if (FindRelationship(a.GetDirection(), b.GetOrigin() - a.GetOrigin()) == Vector<T, dim>::Relationship::Parallel) {
      if (a == b) {
        return Segment<T, dim>::Relationship::Identical;
      }
      if (a.Contains(b.GetLeft()) || a.Contains(b.GetRight()) || b.Contains(a.GetLeft()) || b.Contains(a.GetRight())) {
        return Segment<T, dim>::Relationship::Intersecting;
      }
    }
    return Segment<T, dim>::Relationship::Parallel;;
  }
  if (Comparator<T>::IsZero(a.Distance(b))) {
    return Segment<T, dim>::Relationship::Intersecting;
  }
  return Segment<T, dim>::Relationship::Skew;
}

#endif //GEOMERTY_GEOMETRY_SEGMENT_H_
