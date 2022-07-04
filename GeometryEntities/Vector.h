//
// Created by Artem Novikov on 16.05.2022.
//

#include <iostream>
#include <vector>
#include <cmath>
#include "is_iterable.h"
#include "Data.h"

#ifndef GEOMERTY_GEOMETRY_VECTOR_H_
#define GEOMERTY_GEOMETRY_VECTOR_H_

enum class VectorRelationship {
  Parallel,
  Orthogonal,
  Identical,
  None,
};

template <typename T, size_t Dimension>
class Vector {
 public:
  /// construction
  Vector() = default;

  Vector(const Data<T, Dimension>& data);
  template <typename U>
  Vector(const Data<U, Dimension>& data);

  template <typename U, template <typename, typename...> class Container, typename... Args,
      typename = std::enable_if_t<is_iterable_v<Container<U, Args...>>>>
  explicit Vector(const Container<U, Args...>& data);

  Vector(std::initializer_list<T> list);

  template <typename... Args, typename = std::enable_if_t<sizeof...(Args) == Dimension>>
  Vector(Args&& ... args);

  Vector(const Vector& other) noexcept = default;

  Vector(Vector&& other) noexcept = default;

  template <typename U>
  Vector(const Vector<U, Dimension>& other);

  ~Vector() = default;

  Vector& operator=(const Vector& other) = default;

  Vector& operator=(Vector&& other) noexcept = default;

  /// calc
  T SquaredLength() const noexcept;

  T Length() const noexcept;

  Vector& Normalise();

  Vector& ClampMagnitude(const T& new_magnitude);

  /// getters and setters
  Data<T, Dimension>& data() noexcept { return point_; }

  const Data<T, Dimension>& data() const noexcept { return point_; }

  T& operator[](size_t index) { return point_.cords[index]; }

  const T& operator[](size_t index) const { return point_.cords[index]; }

  T& at(size_t index);

  const T& at(size_t index) const;

  /// arithmetic
  Vector& operator+=(const Vector& other);
  Vector& operator-=(const Vector& other);
  Vector& operator*=(const Vector& other);
  Vector& operator*=(const T& scalar);
  Vector& operator/=(const T& scalar);

 private:
  template <typename U>
  void PushCord(U&& value);
  template <typename U, typename... Args>
  void PushCord(U&& value, Args&& ...args);

  Data<T, Dimension> point_;
};

template <typename T>
using Vector1 = Vector<T, 1>;

using Vector1i = Vector1<int>;
using Vector1f = Vector1<float>;
using Vector1d = Vector1<double>;

template <typename T>
using Vector2 = Vector<T, 2>;

using Vector2i = Vector2<int>;
using Vector2f = Vector2<float>;
using Vector2d = Vector2<double>;

template <typename T>
using Vector3 = Vector<T, 3>;

using Vector3i = Vector3<int>;
using Vector3f = Vector3<float>;
using Vector3d = Vector3<double>;

/// arithmetic

template <typename T, size_t Dimension>
Vector<T, Dimension> operator+(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
Vector<T, Dimension> operator-(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
T operator*(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, Dimension> operator*(const Vector<T, Dimension>& v, const U& scalar);

template <typename T, size_t Dimension, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, Dimension> operator*(const U& scalar, const Vector<T, Dimension>& v);

template <typename T, size_t Dimension, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, Dimension> operator/(const Vector<T, Dimension>& v, const U& scalar);

/// vector calc

template <typename T, size_t Dimension>
T DotProduct(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T>
Vector<T, 3> CrossProduct(const Vector<T, 3>& a, const Vector<T, 3>& b);

template <typename T>
T CrossProduct(const Vector<T, 2>& a, const Vector<T, 2>& b);

template <typename T>
T MixedProduct(const Vector<T, 3>& a, const Vector<T, 3>& b, const Vector<T, 3>& c);

template <typename T, size_t Dimension>
Vector<T, Dimension> Normalised(const Vector<T, Dimension>& v);

template <typename T, size_t Dimension>
double Cos(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
double Sin(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
double Angle(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

/// relation

template <typename T, size_t Dimension>
bool operator==(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
bool operator!=(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

template <typename T, size_t Dimension>
VectorRelationship FindRelationship(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b);

/// stream

template <typename T, size_t Dimension>
std::ostream& operator<<(std::ostream& out, const Vector<T, Dimension>& v);

template <typename T, size_t Dimension>
std::istream& operator>>(std::istream& in, Vector<T, Dimension>& v);

////////////////////////////////////////////////////DEFINITION//////////////////////////////////////////////////////////

template <typename T, size_t Dimension>
Vector<T, Dimension>::Vector(const Data<T, Dimension>& data) : point_(data) {}

template <typename T, size_t Dimension>
template <typename U>
Vector<T, Dimension>::Vector(const Data<U, Dimension>& data) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) = data.cords[i];
  }
}

template <typename T, size_t Dimension>
template <typename U, template <typename, typename...> class Container, typename... Args, typename>
Vector<T, Dimension>::Vector(const Container<U, Args...>& data) {
  assert(data.size() <= Dimension);
  size_t index = 0;
  for(auto& el: data){
    operator[](index++) = el;
  }
}

template <typename T, size_t Dimension>
Vector<T, Dimension>::Vector(std::initializer_list<T> list) {
  assert(list.size() <= Dimension);
  size_t index = 0;
  for (const auto& el : list) {
    operator[](index++) = el;
  }
}

template <typename T, size_t Dimension>
template <typename... Args, typename P>
Vector<T, Dimension>::Vector(Args&& ... args) {
  PushCord(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension>
template <typename U>
Vector<T, Dimension>::Vector(const Vector<U, Dimension>& other) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) = other[i];
  }
}

template <typename T, size_t Dimension>
T Vector<T, Dimension>::SquaredLength() const noexcept {
  T squared_length = DotProduct(*this, *this);
  return squared_length;
}

template <typename T, size_t Dimension>
T Vector<T, Dimension>::Length() const noexcept {
  return std::sqrt(SquaredLength());
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::Normalise() {
  *this /= Length();
  return *this;
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::ClampMagnitude(const T& new_magnitude) {
  *this *= new_magnitude / Length();
  return *this;
}

template <typename T, size_t Dimension>
T& Vector<T, Dimension>::at(size_t index) {
  if (index >= Dimension) {
    throw std::out_of_range("bad index");
  }
  return operator[](index);
}

template <typename T, size_t Dimension>
const T& Vector<T, Dimension>::at(size_t index) const {
  if (index >= Dimension) {
    throw std::out_of_range("bad index");
  }
  return operator[](index);
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::operator+=(const Vector& other) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) += other[i];
  }
  return *this;
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::operator-=(const Vector& other) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) -= other[i];
  }
  return *this;
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::operator*=(const Vector& other) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) *= other[i];
  }
  return *this;
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::operator*=(const T& scalar) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) *= scalar;
  }
  return *this;
}

template <typename T, size_t Dimension>
Vector<T, Dimension>& Vector<T, Dimension>::operator/=(const T& scalar) {
  for (size_t i = 0; i < Dimension; ++i) {
    operator[](i) /= scalar;
  }
  return *this;
}

template <typename T, size_t Dimension>
template <typename U>
void Vector<T, Dimension>::PushCord(U&& value) {
  operator[](Dimension - 1) = value;
}

template <typename T, size_t Dimension>
template <typename U, typename... Args>
void Vector<T, Dimension>::PushCord(U&& value, Args&& ... args) {
  operator[](Dimension - sizeof...(Args) - 1) = value;
  PushCord(std::forward<Args>(args)...);
}

template <typename T, size_t Dimension>
Vector<T, Dimension> operator+(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  auto copy = a;
  copy += b;
  return copy;
}

template <typename T, size_t Dimension>
Vector<T, Dimension> operator-(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  auto copy = a;
  copy -= b;
  return copy;
}

template <typename T, size_t Dimension>
T operator*(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  T product = a[0] * b[0];
  for (size_t i = 1; i < Dimension; ++i) {
    product += a[i] * b[i];
  }
  return product;
}

template <typename T, size_t Dimension, typename U, typename>
Vector<T, Dimension> operator*(const Vector<T, Dimension>& v, const U& scalar) {
  auto copy = v;
  copy *= scalar;
  return copy;
}

template <typename T, size_t Dimension, typename U, typename>
Vector<T, Dimension> operator*(const U& scalar, const Vector<T, Dimension>& v) {
  auto copy = v;
  copy *= scalar;
  return copy;
}

template <typename T, size_t Dimension, typename U, typename>
Vector<T, Dimension> operator/(const Vector<T, Dimension>& v, const U& scalar) {
  auto copy = v;
  copy /= scalar;
  return copy;
}

template <typename T, size_t Dimension>
T DotProduct(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  auto dot = a * b;
  return dot;
}

template <typename T>
Vector<T, 3> CrossProduct(const Vector<T, 3>& a, const Vector<T, 3>& b) {
  Vector<T, 3> c(a.data().y * b.data().z - a.data().z * b.data().y,
                 a.data().z * b.data().x - a.data().x * b.data().z,
                 a.data().x * b.data().y - a.data().y * b.data().x);
  return c;
}

template <typename T>
T CrossProduct(const Vector<T, 2>& a, const Vector<T, 2>& b) {
  T result = a.data().x * b.data().y - a.data().y * b.data().x;
  return result;
}

template <typename T>
T MixedProduct(const Vector<T, 3>& a, const Vector<T, 3>& b, const Vector<T, 3>& c) {
  T result = a.data().x * (b.data().y * c.data().z - b.data().z * c.data().y) +
      a.data().y * (b.data().z * c.data().x - b.data().x * c.data().z) +
      a.data().z * (b.data().x * c.data().y - b.data().y * c.data().x);
  return result;
}

template <typename T, size_t Dimension>
Vector<T, Dimension> Normalised(const Vector<T, Dimension>& v) {
  Vector<T, Dimension> copy = v;
  copy.Normalise();
  return copy;
}

template <typename T, size_t Dimension>
double Cos(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  double dot = DotProduct(a, b);
  double cos = std::sqrt(dot * dot / a.SquaredLength() / b.SquaredLength());
  return cos;
}

template <typename T, size_t Dimension>
double Sin(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  double cos = Cos(a, b);
  return std::sqrt(1 - cos * cos);
}

template <typename T>
double Sin(const Vector<T, 3>& a, const Vector<T, 3>& b) {
  double cross = CrossProduct(a, b).SquaredLength();
  double sin = std::sqrt(cross / a.SquaredLength() / b.SquaredLength());
  return sin;
}

template <typename T, size_t Dimension>
double Angle(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  return std::acos(Cos(a, b));
}

template <typename T, size_t Dimension>
bool operator==(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  return a.data() == b.data();
}

template <typename T, size_t Dimension>
bool operator!=(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  return a.data() != b.data();
}

template <typename T, size_t Dimension>
VectorRelationship FindRelationship(const Vector<T, Dimension>& a, const Vector<T, Dimension>& b) {
  T dot = a * b;
  if (Comparator<T>::Equal(dot * dot, a.SquaredLength() * b.SquaredLength())) {
    return VectorRelationship::Parallel;
  }

  if (Comparator<T>::IsZero(dot)) {
    return VectorRelationship::Orthogonal;
  }

  return VectorRelationship::None;
}

template <typename T, size_t Dimension>
std::ostream& operator<<(std::ostream& out, const Vector<T, Dimension>& v) {
  out << '(' << v.data().cords[0];
  for (size_t i = 1; i < Dimension; ++i) {
    out << ", " << v[i];
  }
  out << ')';
  return out;
}

template <typename T, size_t Dimension>
std::istream& operator>>(std::istream& in, Vector<T, Dimension>& v) {
  for (size_t i = 0; i < Dimension; ++i) {
    in >> v[i];
  }
  return in;
}

#endif //GEOMERTY_GEOMETRY_VECTOR_H_