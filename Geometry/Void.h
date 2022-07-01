//
// Created by Artem Novikov on 28.06.2022.
//

#ifndef GEOMERTY_GEOMETRY_VOID_H_
#define GEOMERTY_GEOMETRY_VOID_H_

enum class Entity {
  Void,
  Point,
  Segment,
  Line,
  Plane,
};

template <typename T, size_t dim>
class Void {
 public:
  Void() : type_(Entity::Void) {}
  Void(Entity type) : type_(type) {}
  [[nodiscard]] Entity GetType() const { return type_; }
  virtual  ~Void() = default;
 protected:
  Entity type_ = Entity::Void;
};

template <typename T>
using Void2 = Void<T, 2>;

using Void2i = Void2<int>;
using Void2f = Void2<float>;
using Void2d = Void2<double>;

template <typename T>
using Void3 = Void<T, 3>;

using Void3i = Void3<int>;
using Void3f = Void3<float>;
using Void3d = Void3<double>;

#endif //GEOMERTY_GEOMETRY_VOID_H_
