//
// Created by Artem Novikov on 24.06.2022.
//

#include "Void.h"
#include "Vector.h"
#include "Point.h"
#include "Line.h"
#include "Segment.h"

#ifndef GEOMERTY_GEOMETRY_PLANE_H_
#define GEOMERTY_GEOMETRY_PLANE_H_

template <typename T>
class Plane : public Void<T, 3> {
 public:
  enum class Relationship {
    Parallel,
    Intersecting,
    Identical,
  };

  /// construction
  Plane() : Void<T, 3>(Entity::Plane) {}
  Plane(const Point3<T>& origin, const Vector3<T>& abscissa, const Vector3<T>& ordinates);
  Plane(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3);
  Plane(const Plane& other) = default;
  Plane(Plane&& other) noexcept = default;
  Plane& operator=(const Plane& other) = default;
  Plane& operator=(Plane&& other) noexcept = default;

  /// setters and getters
  Point3<T> GetPoint(T x, T y) const { return origin_ + abscissa_ * x + ordinate_ * y; }
  Point3<T>& GetOrigin() { return origin_; }
  const Point3<T>& GetOrigin() const { return origin_; }
  Point3<T>& GetAbscissa() { return abscissa_; }
  const Point3<T>& GetAbscissa() const { return abscissa_; }
  Point3<T>& GetOrdinate() { return ordinate_; }
  const Point3<T>& GetOrdinate() const { return ordinate_; }
  Vector3<T> GetNormal() const { return CrossProduct(abscissa_, ordinate_); }

  /// plane calc
  void Normalise();
  bool Contains(const Point3<T>& point) const;
  bool Contains(const Segment3<T>& segment) const;
  bool Contains(const Line3<T>& line) const;

  T SquaredDistance(const Point3<T>& point) const;
  T SquaredDistance(const Segment3<T>& segment) const;
  T SquaredDistance(const Line3<T>& line) const;
  T SquaredDistance(const Plane<T>& plane) const;

  T Distance(const Point3<T>& point) const;
  T Distance(const Segment3<T>& segment) const;
  T Distance(const Line3<T>& line) const;
  T Distance(const Plane<T>& plane) const;

  std::unique_ptr<Void3<T>> Intersection(const Point3<T>& point) const;
  std::unique_ptr<Void3<T>> Intersection(const Segment3<T>& segment) const;
  std::unique_ptr<Void3<T>> Intersection(const Line3<T>& line) const;
  std::unique_ptr<Void3<T>> Intersection(const Plane<T>& plane) const;

  Point3<T> Projection(const Point3<T>& point, const Vector3<T>& a) const;
  Segment3<T> Projection(const Segment3<T>& segment, const Vector3<T>& a) const;
  Line3<T> Projection(const Line3<T>& line, const Vector3<T>& a) const;

  Point3<T> Projection(const Point3<T>& point) const;
  Segment3<T> Projection(const Segment3<T>& segment) const;
  Line3<T> Projection(const Line3<T>& line) const;

 private:
  Point3<T> origin_;
  Vector3<T> abscissa_;
  Vector3<T> ordinate_;
};

/// Relationships

template <typename T>
bool operator==(const Plane<T>& a, const Plane<T>& b);

template <typename T>
bool operator!=(const Plane<T>& a, const Plane<T>& b);

template <typename T>
typename Plane<T>::Relationship FindRelationship(const Plane<T>& a, const Plane<T>& b);

//////////////////////////////////////////////////////DEFINITION////////////////////////////////////////////////////////

template <typename T>
Plane<T>::Plane(const Point3<T>& origin, const Vector3<T>& abscissa, const Vector3<T>& ordinates)
    : Void<T, 3>(Entity::Plane), origin_(origin), abscissa_(abscissa), ordinate_(ordinates) {}

template <typename T>
Plane<T>::Plane(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3)
    : Void<T, 3>(Entity::Plane), origin_(p1), abscissa_(p2 - p1), ordinate_(p3 - p1) {}

template <typename T>
void Plane<T>::Normalise() {
  abscissa_.Normalise();
  ordinate_.Normalise();
}

template <typename T>
bool Plane<T>::Contains(const Point3<T>& point) const {
  return Comparator<T>::IsZero(MixedProduct(point - origin_, abscissa_, ordinate_));
}

template <typename T>
bool Plane<T>::Contains(const Segment3<T>& segment) const {
  return Contains(segment.GetLeft()) && Contains(segment.GetRight());
}

template <typename T>
bool Plane<T>::Contains(const Line3<T>& line) const {
  return Comparator<T>::IsZero(MixedProduct(line.GetOrigin() - origin_, abscissa_, ordinate_)) &&
      Comparator<T>::IsZero(MixedProduct(line.GetDirection() - origin_, abscissa_, ordinate_));
}
template <typename T>
T Plane<T>::SquaredDistance(const Point3<T>& point) const {
  return std::pow(MixedProduct(point - origin_, abscissa_, ordinate_), 2)
      / CrossProduct(abscissa_, ordinate_).SquaredLength();
}
template <typename T>
T Plane<T>::Distance(const Point3<T>& point) const {
  return std::sqrt(SquaredDistance(point));
}

template <typename T>
T Plane<T>::SquaredDistance(const Segment3<T>& segment) const {
  if (Comparator<T>::IsZero(MixedProduct(segment.GetDirection(), abscissa_, ordinate_))) {
    return SquaredDistance(segment.GetLeft());
  }
  T m1 = MixedProduct(segment.GetLeft(), abscissa_, ordinate_);
  T m2 = MixedProduct(segment.GetRight(), abscissa_, ordinate_);
  if (m1 * m2 < 0) {
    return 0;
  }
  return std::min(m1 * m1, m2 * m2) / CrossProduct(abscissa_, ordinate_).SquaredLength();
}

template <typename T>
T Plane<T>::Distance(const Segment3<T>& segment) const {
  return std::sqrt(SquaredDistance(segment));
}

template <typename T>
T Plane<T>::SquaredDistance(const Line3<T>& line) const {
  return std::pow(Distance(line), 2);
}

template <typename T>
T Plane<T>::Distance(const Line3<T>& line) const {
  return Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_)) ?
         Distance(line.GetOrigin()) : 0;
}

template <typename T>
T Plane<T>::SquaredDistance(const Plane<T>& plane) const {
  return std::pow(Distance(plane), 2);
}

template <typename T>
T Plane<T>::Distance(const Plane<T>& plane) const {
  return FindRelationship(CrossProduct(abscissa_, ordinate_), CrossProduct(plane.GetAbscissa(), plane.GetOrdinate()))
             == Vector3<T>::Relationship::Parallel ? Distance(plane.GetOrigin()) : 0;
}

template <typename T>
std::unique_ptr<Void3<T>> Plane<T>::Intersection(const Point3<T>& point) const {
  if (Contains(point)) {
    return std::make_unique<Point3<T>>(point);
  }
  return std::make_unique<Void3<T>>();
}

template <typename T>
std::unique_ptr<Void3<T>> Plane<T>::Intersection(const Segment3<T>& segment) const {
  T m1 = MixedProduct(segment.GetRight() - origin_, abscissa_, ordinate_);
  T m2 = MixedProduct(segment.GetLeft() - origin_, abscissa_, ordinate_);
  if (Comparator<T>::IsZero(m1) && Comparator<T>::IsZero(m2)) {
    return std::unique_ptr<Segment3<T>>(segment);
  }
  if (m1 * m2 < 0) {
    T t = MixedProduct(origin_ - segment.GetLeft(), abscissa_, ordinate_)
        / MixedProduct(segment.GetDirection(), abscissa_, ordinate_);
    return std::unique_ptr<Point3<T>>(segment.GetPoint(t));
  }
  return std::make_unique<Void3<T>>();
}

template <typename T>
std::unique_ptr<Void3<T>> Plane<T>::Intersection(const Line3<T>& line) const {
  if (Contains(line.GetOrigin())) {
    if (Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_))) {
      return std::make_unique<Line3<T>>(line);
    }
    return std::make_unique<Point3<T>>(line.GetOrigin());
  }
  if (Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_))) {
    return std::make_unique<Void3<T>>();
  }
  T t = MixedProduct(origin_ - line.GetOrigin(), abscissa_, ordinate_)
      / MixedProduct(line.GetDirection(), abscissa_, ordinate_);
  return std::make_unique<Point3<T>>(line.GetPoint(t));
}

template <typename T>
std::unique_ptr<Void3<T>> Plane<T>::Intersection(const Plane<T>& plane) const {
  switch (FindRelationship(*this, plane)) {
    case Plane<T>::Relationship::Intersecting: {
      Vector3<T> dir = CrossProduct(CrossProduct(abscissa_, ordinate_), CrossProduct(plane.abscissa_, plane.ordinate_));
      Line3<T> l(plane.origin_, plane.abscissa_);
      if (Comparator<T>::IsZero(MixedProduct(plane.abscissa_, abscissa_, ordinate_))) {
        l.GetDirection() = plane.ordinate_;
      }
      T t = MixedProduct(origin_ - l.GetOrigin(), abscissa_, ordinate_)
          / MixedProduct(l.GetDirection(), abscissa_, ordinate_);
      return std::make_unique<Line3<T>>(l.GetPoint(t), dir);
    }
    case Plane<T>::Relationship::Identical: { return std::make_unique<Plane<T>>(plane); }
    case Plane<T>::Relationship::Parallel: { return std::make_unique<Void3<T>>(); }
  }
}

template <typename T>
Point3<T> Plane<T>::Projection(const Point3<T>& point, const Vector3<T>& a) const {
  T t = MixedProduct(origin_ - point, abscissa_, ordinate_)
      / MixedProduct(a, abscissa_, ordinate_);
  return point + a * t;
}

template <typename T>
Segment3<T> Plane<T>::Projection(const Segment3<T>& segment, const Vector3<T>& a) const {
  return {Projection(segment.GetLeft(), a), Projection(segment.GetRight(), a)};
}

template <typename T>
Line3<T> Plane<T>::Projection(const Line3<T>& line, const Vector3<T>& a) const {
  return {Projection(line.GetOrigin(), a), Projection(line.GetPoint(1), a)};
}

template <typename T>
Point3<T> Plane<T>::Projection(const Point3<T>& point) const {
  auto normal = CrossProduct(abscissa_, ordinate_);
  return point - normal * (DotProduct(point - origin_, normal) / normal.SquaredLength());
}

template <typename T>
Segment3<T> Plane<T>::Projection(const Segment3<T>& segment) const {
  return {Projection(segment.GetLeft()), Projection(segment.GetRight())};
}

template <typename T>
Line3<T> Plane<T>::Projection(const Line3<T>& line) const {
  return {Projection(line.GetOrigin()), Projection(line.GetPoint(1))};
}

template <typename T>
bool operator==(const Plane<T>& a, const Plane<T>& b) {
  return FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                          CrossProduct(b.GetAbscissa(), b.GetOrdinate())) ==
      Vector3<T>::Relationship::Parallel && a.Contains(b.GetOrigin());
}

template <typename T>
bool operator!=(const Plane<T>& a, const Plane<T>& b) {
  return !a.Contains(b.GetOrigin()) || FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                                                        CrossProduct(b.GetAbscissa(), b.GetOrdinate())) !=
      Vector3<T>::Relationship::Parallel;
}

template <typename T>
typename Plane<T>::Relationship FindRelationship(const Plane<T>& a, const Plane<T>& b) {
  if (FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                       CrossProduct(b.GetAbscissa(), b.GetOrdinate())) == Vector3<T>::Relationship::Parallel) {
    if (a.Contains(b.GetOrigin())) {
      return Plane<T>::Relationship::Identical;
    }
    return Plane<T>::Relationship::Parallel;
  }
  return Plane<T>::Relationship::Intersecting;
}

#endif //GEOMERTY_GEOMETRY_PLANE_H_
