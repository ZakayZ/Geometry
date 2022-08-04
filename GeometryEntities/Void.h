//
// Created by Artem Novikov on 04.08.2022.
//

#ifndef GEOMETRY_GEOMETRYENTITIES_ENTITY_H_
#define GEOMETRY_GEOMETRYENTITIES_ENTITY_H_

#include <cstddef>
#include <stdexcept>
#include <memory>
#include "Transform.h"
#include "Entities.h"

template <typename T, size_t Dimension>
class Void {
 public:
  [[nodiscard]] virtual Entity GetType() const { return Entity::Void; }

  virtual bool Contains(const Void& object) const { return false; }

  virtual T SquaredDistance(const Void& object) const {
    throw std::runtime_error("One object is Void");
    return 0;
  }
  virtual T Distance(const Void& object) const {
    throw std::runtime_error("One object is Void");
    return 0;
  }
  virtual std::unique_ptr<Void> Intersection(const Void& object) const { return std::make_unique<Void>(); }

  virtual void ApplyTransform(const Transform<T, Dimension>& transform) {}

  virtual ~Void() = default;
};

template <typename T>
using Void1 = Void<T, 1>;

using Void1i = Void1<int>;
using Void1f = Void1<float>;
using Void1d = Void1<double>;

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

#endif //GEOMETRY_GEOMETRYENTITIES_ENTITY_H_
