//
// Created by Artem Novikov on 24.06.2022.
//

#include "Vector.h"
#include "Point.h"

#ifndef GEOMERTY_GEOMETRY_PLANE_H_
#define GEOMERTY_GEOMETRY_PLANE_H_

template <typename T>
class Plane {
 public:
  enum class Relationship {
    Parallel,
    Intersecting,
    Identical,
  };

  /// construction
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
  bool Contains(const Line3<T>& line) const;
  T SquaredDistance(const Point3<T>& point) const;
  T Distance(const Point3<T>& point) const;
  T SquaredDistance(const Line3<T>& line) const;
  T Distance(const Line3<T>& line) const;
  T SquaredDistance(const Plane<T>& plane) const;
  T Distance(const Plane<T>& plane) const;

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
  return Comparator<T>::Equal(MixedProduct(point - origin_, abscissa_, ordinate_), 0);
}

template <typename T>
bool Plane<T>::Contains(const Line3<T>& line) const {
  return Comparator<T>::Equal(MixedProduct(line.GetOrigin() - origin_, abscissa_, ordinate_), 0) &&
      Comparator<T>::Equal(MixedProduct(line.GetDirection() - origin_, abscissa_, ordinate_), 0);
}
template <typename T>
T Plane<T>::SquaredDistance(const Point3<T>& point) const {
  return std::pow(Distance(point), 2);
}
template <typename T>
T Plane<T>::Distance(const Point3<T>& point) const {
  T distance = MixedProduct(point - origin_, abscissa_, ordinate_) / CrossProduct(abscissa_, ordinate_);
  return distance;
}

template <typename T>
T Plane<T>::SquaredDistance(const Line3<T>& line) const {
  T squared_distance = std::pow(Distance(line), 2);
  return squared_distance;
}

template <typename T>
T Plane<T>::Distance(const Line3<T>& line) const {
  T distance = Comparator<T>::Equal(CrossProduct(line.GetDirection(), abscissa_, ordinate_), 0) ?
               Distance(line.GetOrigin()) : 0;
  return distance;
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
