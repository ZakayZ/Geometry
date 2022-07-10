//
// Created by Artem Novikov on 24.06.2022.
//

#include "GeometricEntity.h"
#include "Void.h"
#include "Vector.h"
#include "Point.h"
#include "Line.h"
#include "Segment.h"
#include "BoundaryBox.h"
#include "Transform.h"

#ifndef GEOMETRY_GEOMETRY_PLANE_H_
#define GEOMETRY_GEOMETRY_PLANE_H_

enum class PlaneRelationship {
  Parallel,
  Intersecting,
  Identical,
};

template <typename T>
class Plane {
 public:
  /// construction
  Plane() = default;
  Plane(const Point3<T>& origin, const Vector3<T>& abscissa, const Vector3<T>& ordinates);
  Plane(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3);
  Plane(const Plane& other) = default;
  Plane(Plane&& other) noexcept = default;
  Plane& operator=(const Plane& other) = default;
  Plane& operator=(Plane&& other) noexcept = default;

  /// setters and getters
  Entity GetType() const { return Entity::Plane; }
  size_t GetDimension() const { return 3; }
  Point3<T> GetPoint(T x, T y) const { return origin_ + abscissa_ * x + ordinate_ * y; }
  Point3<T>& GetOrigin() { return origin_; }
  const Point3<T>& GetOrigin() const { return origin_; }
  Vector3<T>& GetAbscissa() { return abscissa_; }
  const Vector3<T>& GetAbscissa() const { return abscissa_; }
  Vector3<T>& GetOrdinate() { return ordinate_; }
  const Vector3<T>& GetOrdinate() const { return ordinate_; }
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

  GeometryEntity Intersection(const Point3<T>& point) const;
  GeometryEntity Intersection(const Segment3<T>& segment) const;
  GeometryEntity Intersection(const Line3<T>& line) const;
  GeometryEntity Intersection(const Plane<T>& plane) const;

  Point3<T> Projection(const Point3<T>& point, const Vector3<T>& a) const;
  Segment3<T> Projection(const Segment3<T>& segment, const Vector3<T>& a) const;
  Line3<T> Projection(const Line3<T>& line, const Vector3<T>& a) const;

  Point3<T> Projection(const Point3<T>& point) const;
  Segment3<T> Projection(const Segment3<T>& segment) const;
  Line3<T> Projection(const Line3<T>& line) const;

  bool Intersects(const BoundaryBox3<T>& box) const;

  template <size_t OutputDimension>
  Segment<T, OutputDimension> ApplyTransform(const Transform<T, 3, OutputDimension>& transform) const;

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
PlaneRelationship FindRelationship(const Plane<T>& a, const Plane<T>& b);

//////////////////////////////////////////////////////DEFINITION////////////////////////////////////////////////////////

template <typename T>
Plane<T>::Plane(const Point3<T>& origin, const Vector3<T>& abscissa, const Vector3<T>& ordinates)
    : origin_(origin), abscissa_(abscissa), ordinate_(ordinates) {}

template <typename T>
Plane<T>::Plane(const Point3<T>& p1, const Point3<T>& p2, const Point3<T>& p3)
    : origin_(p1), abscissa_(p2 - p1), ordinate_(p3 - p1) {}

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
  return Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_)) ?
         SquaredDistance(line.GetOrigin()) : 0;
}

template <typename T>
T Plane<T>::Distance(const Line3<T>& line) const {
  return std::sqrt(SquaredDistance(line));
}

template <typename T>
T Plane<T>::SquaredDistance(const Plane<T>& plane) const {
  return FindRelationship(CrossProduct(abscissa_, ordinate_), CrossProduct(plane.GetAbscissa(), plane.GetOrdinate()))
             == VectorRelationship::Parallel ? SquaredDistance(plane.GetOrigin()) : 0;
}

template <typename T>
T Plane<T>::Distance(const Plane<T>& plane) const {
  return std::sqrt(SquaredDistance(plane));
}

template <typename T>
GeometryEntity Plane<T>::Intersection(const Point3<T>& point) const {
  if (Contains(point)) {
    return GeometryEntity(point);
  }
  return MakeGeometryEntity<Void3<T>>();
}

template <typename T>
GeometryEntity Plane<T>::Intersection(const Segment3<T>& segment) const {
  T m1 = MixedProduct(segment.GetRight() - origin_, abscissa_, ordinate_);
  T m2 = MixedProduct(segment.GetLeft() - origin_, abscissa_, ordinate_);
  if (Comparator<T>::IsZero(m1) && Comparator<T>::IsZero(m2)) {
    return GeometryEntity(segment);
  }
  if (m1 * m2 <= 0) {
    T t = MixedProduct(origin_ - segment.GetLeft(), abscissa_, ordinate_)
        / MixedProduct(segment.GetDirection(), abscissa_, ordinate_);
    return GeometryEntity(segment.GetPoint(t));
  }
  return MakeGeometryEntity<Void3<T>>();
}

template <typename T>
GeometryEntity Plane<T>::Intersection(const Line3<T>& line) const {
  if (Contains(line.GetOrigin())) {
    if (Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_))) {
      return GeometryEntity(line);
    }
    return GeometryEntity(line.GetOrigin());
  }
  if (Comparator<T>::IsZero(MixedProduct(line.GetDirection(), abscissa_, ordinate_))) {
    return MakeGeometryEntity<Void3<T>>();
  }
  T t = MixedProduct(origin_ - line.GetOrigin(), abscissa_, ordinate_)
      / MixedProduct(line.GetDirection(), abscissa_, ordinate_);
  return GeometryEntity(line.GetPoint(t));
}

template <typename T>
GeometryEntity Plane<T>::Intersection(const Plane<T>& plane) const {
  switch (FindRelationship(*this, plane)) {
    case PlaneRelationship::Intersecting: {
      Vector3<T> dir = CrossProduct(CrossProduct(abscissa_, ordinate_), CrossProduct(plane.abscissa_, plane.ordinate_));
      Line3<T> l(plane.origin_, plane.abscissa_);
      if (Comparator<T>::IsZero(MixedProduct(plane.abscissa_, abscissa_, ordinate_))) {
        l.GetDirection() = plane.ordinate_;
      }
      T t = MixedProduct(origin_ - l.GetOrigin(), abscissa_, ordinate_)
          / MixedProduct(l.GetDirection(), abscissa_, ordinate_);
      return MakeGeometryEntity<Line3<T>>(l.GetPoint(t), std::move(dir));
    }
    case PlaneRelationship::Identical: { return GeometryEntity(plane); }
    case PlaneRelationship::Parallel: { return MakeGeometryEntity<Void3<T>>(); }
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
bool Plane<T>::Intersects(const BoundaryBox3<T>& box) const {
  auto p1 = box.GetLeft();
  auto p2 = box.GetRight();
  if (MixedProduct(p1 - origin_, abscissa_, ordinate_) * MixedProduct(p1 - origin_, abscissa_, ordinate_) <= 0) {
    return true;
  }
  p1[0] = box.GetRight()[0];
  p2[0] = box.GetLeft()[0];
  if (MixedProduct(p1 - origin_, abscissa_, ordinate_) * MixedProduct(p1 - origin_, abscissa_, ordinate_) <= 0) {
    return true;
  }
  p1[1] = box.GetRight()[1];
  p2[1] = box.GetLeft()[1];
  if (MixedProduct(p1 - origin_, abscissa_, ordinate_) * MixedProduct(p1 - origin_, abscissa_, ordinate_) <= 0) {
    return true;
  }
  p1[0] = box.GetLeft()[0];
  p2[0] = box.GetRight()[0];
  if (MixedProduct(p1 - origin_, abscissa_, ordinate_) * MixedProduct(p1 - origin_, abscissa_, ordinate_) <= 0) {
    return true;
  }
  return false;
}

template <typename T>
template <size_t OutputDimension>
Segment<T, OutputDimension> Plane<T>::ApplyTransform(const Transform<T, 3, OutputDimension>& transform) const {
  return {transform(origin_), transform(abscissa_), transform(ordinate_)};
}

template <typename T>
bool operator==(const Plane<T>& a, const Plane<T>& b) {
  return FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                          CrossProduct(b.GetAbscissa(), b.GetOrdinate())) ==
      VectorRelationship::Parallel && a.Contains(b.GetOrigin());
}

template <typename T>
bool operator!=(const Plane<T>& a, const Plane<T>& b) {
  return !a.Contains(b.GetOrigin()) || FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                                                        CrossProduct(b.GetAbscissa(), b.GetOrdinate())) !=
      VectorRelationship::Parallel;
}

template <typename T>
PlaneRelationship FindRelationship(const Plane<T>& a, const Plane<T>& b) {
  if (FindRelationship(CrossProduct(a.GetAbscissa(), a.GetOrdinate()),
                       CrossProduct(b.GetAbscissa(), b.GetOrdinate())) == VectorRelationship::Parallel) {
    if (a.Contains(b.GetOrigin())) {
      return PlaneRelationship::Identical;
    }
    return PlaneRelationship::Parallel;
  }
  return PlaneRelationship::Intersecting;
}

#endif //GEOMETRY_GEOMETRY_PLANE_H_
