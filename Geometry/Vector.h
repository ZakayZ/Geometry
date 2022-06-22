//
// Created by Artem Novikov on 16.05.2022.
//

#include <iostream>
#include <vector>
#include <cmath>
#include "Data.h"

#ifndef GEOMERTY_GEOMETRY_VECTOR_H_
#define GEOMERTY_GEOMETRY_VECTOR_H_

template <typename... Args>
static const size_t size_v = sizeof...(Args);

template <typename T, size_t dim>
class Vector {
 public:
  enum class Relationship {
    Parallel,
    Orthogonal,
    Identical,
    None,
  };

  /// construction
  Vector() = default;

  Vector(const Data<T, dim>& data);
  template <typename U>
  Vector(const Data<U, dim>& data);

  template <typename U>
  Vector(const std::vector<U>& data);

  Vector(const std::initializer_list<T>& data);

  template <typename... Args, typename = std::enable_if_t<size_v<Args...> == dim>>
  Vector(Args&& ... args);

  Vector(const Vector& other) noexcept = default;

  Vector(Vector&& other) noexcept = default;

  template <typename U>
  Vector(const Vector<U, dim>& other);

  ~Vector() = default;

  Vector& operator=(const Vector& other) = default;

  Vector& operator=(Vector&& other) noexcept = default;

  /// calc
  T SquaredLength() const noexcept;

  T Length() const noexcept;

  Vector& Normalise();

  Vector& ClampMagnitude(const T& new_magnitude);

  /// getters and setters
  Data<T, dim>& data() noexcept { return point_; }

  const Data<T, dim>& data() const noexcept { return point_; }

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

  Data<T, dim> point_;
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

template <typename T, size_t dim>
Vector<T, dim> operator+(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
Vector<T, dim> operator-(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
T operator*(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, dim> operator*(const Vector<T, dim>& v, const U& scalar);

template <typename T, size_t dim, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, dim> operator*(const U& scalar, const Vector<T, dim>& v);

template <typename T, size_t dim, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Vector<T, dim> operator/(const Vector<T, dim>& v, const U& scalar);

/// vector calc

template <typename T, size_t dim>
T DotProduct(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T>
Vector<T, 3> CrossProduct(const Vector<T, 3>& a, const Vector<T, 3>& b);

template <typename T>
T CrossProduct(const Vector<T, 2>& a, const Vector<T, 2>& b);

template <typename T>
T MixedProduct(const Vector<T, 3>& a, const Vector<T, 3>& b, const Vector<T, 3>& c);

template <typename T, size_t dim>
Vector<T, dim> Normalised(const Vector<T, dim>& v);

template <typename T, size_t dim>
double Cos(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
double Sin(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
double Angle(const Vector<T, dim>& a, const Vector<T, dim>& b);

/// relation

template <typename T, size_t dim>
bool operator==(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
bool operator!=(const Vector<T, dim>& a, const Vector<T, dim>& b);

template <typename T, size_t dim>
typename Vector<T, dim>::Relationship FindRelationShip(const Vector<T, dim>& a, const Vector<T, dim>& b);

/// stream

template <typename T, size_t dim>
std::ostream& operator<<(std::ostream& out, const Vector<T, dim>& v);

template <typename T, size_t dim>
std::istream& operator>>(std::istream& in, Vector<T, dim>& v);

////////////////////////////////////////////////////DEFINITION//////////////////////////////////////////////////////////

template <typename T, size_t dim>
Vector<T, dim>::Vector(const Data<T, dim>& data) : point_(data) {}

template <typename T, size_t dim>
template <typename U>
Vector<T, dim>::Vector(const Data<U, dim>& data) {
  for (size_t i = 0; i < dim; ++i) {
    point_.cords[i] = data.cords[i];
  }
}

template <typename T, size_t dim>
template <typename U>
Vector<T, dim>::Vector(const std::vector<U>& data) {
  for (size_t i = 0; i < dim; ++i) {
    point_.cords[i] = data[i];
  }
}

template <typename T, size_t dim>
Vector<T, dim>::Vector(const std::initializer_list<T>& data) {
  size_t index = 0;
  for (const auto& el : data) {
    point_.cords[index] = el;
    ++index;
  }
}

template <typename T, size_t dim>
template <typename... Args, typename P>
Vector<T, dim>::Vector(Args&& ... args) {
  PushCord(std::forward<Args>(args)...);
}

template <typename T, size_t dim>
template <typename U>
Vector<T, dim>::Vector(const Vector<U, dim>& other) {
  for (size_t i = 0; i < dim; ++i) {
    point_.cords[i] = other[i];
  }
}

template <typename T, size_t dim>
T Vector<T, dim>::SquaredLength() const noexcept {
  T squared_length = DotProduct(*this, *this);
  return squared_length;
}

template <typename T, size_t dim>
T Vector<T, dim>::Length() const noexcept {
  return std::sqrt(SquaredLength());
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::Normalise() {
  *this /= Length();
  return *this;
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::ClampMagnitude(const T& new_magnitude) {
  *this *= new_magnitude / Length();
  return *this;
}

template <typename T, size_t dim>
T& Vector<T, dim>::at(size_t index) {
  if (index >= dim) {
    throw std::out_of_range("bad index");
  }
  return operator[](index);
}

template <typename T, size_t dim>
const T& Vector<T, dim>::at(size_t index) const {
  if (index >= dim) {
    throw std::out_of_range("bad index");
  }
  return operator[](index);
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::operator+=(const Vector& other) {
  for (size_t i = 0; i < dim; ++i) {
    operator[](i) += other[i];
  }
  return *this;
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::operator-=(const Vector& other) {
  for (size_t i = 0; i < dim; ++i) {
    operator[](i) -= other[i];
  }
  return *this;
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::operator*=(const Vector& other) {
  for (size_t i = 0; i < dim; ++i) {
    operator[](i) *= other[i];
  }
  return *this;
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::operator*=(const T& scalar) {
  for (size_t i = 0; i < dim; ++i) {
    operator[](i) *= scalar;
  }
  return *this;
}

template <typename T, size_t dim>
Vector<T, dim>& Vector<T, dim>::operator/=(const T& scalar) {
  for (size_t i = 0; i < dim; ++i) {
    operator[](i) /= scalar;
  }
  return *this;
}

template <typename T, size_t dim>
template <typename U>
void Vector<T, dim>::PushCord(U&& value) {
  point_.cords[dim - 1] = value;
}

template <typename T, size_t dim>
template <typename U, typename... Args>
void Vector<T, dim>::PushCord(U&& value, Args&& ... args) {
  point_.cords[dim - sizeof...(Args) - 1] = value;
  PushCord(std::forward<Args>(args)...);
}

template <typename T, size_t dim>
Vector<T, dim> operator+(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  auto copy = a;
  copy += b;
  return copy;
}

template <typename T, size_t dim>
Vector<T, dim> operator-(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  auto copy = a;
  copy -= b;
  return copy;
}

template <typename T, size_t dim>
T operator*(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  T product = a[0] * b[0];
  for (size_t i = 1; i < dim; ++i) {
    product += a[i] * b[i];
  }
  return product;
}

template <typename T, size_t dim, typename U, typename>
Vector<T, dim> operator*(const Vector<T, dim>& v, const U& scalar) {
  auto copy = v;
  copy *= scalar;
  return copy;
}

template <typename T, size_t dim, typename U, typename>
Vector<T, dim> operator*(const U& scalar, const Vector<T, dim>& v) {
  auto copy = v;
  copy *= scalar;
  return copy;
}

template <typename T, size_t dim, typename U, typename>
Vector<T, dim> operator/(const Vector<T, dim>& v, const U& scalar) {
  auto copy = v;
  copy /= scalar;
  return copy;
}

template <typename T, size_t dim>
T DotProduct(const Vector<T, dim>& a, const Vector<T, dim>& b) {
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

template <typename T, size_t dim>
Vector<T, dim> Normalised(const Vector<T, dim>& v) {
  Vector<T, dim> copy = v;
  copy.Normalise();
  return copy;
}

template <typename T, size_t dim>
double Cos(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  double dot = DotProduct(a, b);
  double cos = std::sqrt(dot * dot / a.SquaredLength() / b.SquaredLength());
  return cos;
}

template <typename T, size_t dim>
double Sin(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  double cos = Cos(a, b);
  return std::sqrt(1 - cos * cos);
}

template <typename T>
double Sin(const Vector<T, 3>& a, const Vector<T, 3>& b) {
  double cross = CrossProduct(a, b).SquaredLength();
  double sin = std::sqrt(cross / a.SquaredLength() / b.SquaredLength());
  return sin;
}

template <typename T, size_t dim>
double Angle(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  return std::acos(Cos(a, b));
}

template <typename T, size_t dim>
bool operator==(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  return a.data() == b.data();
}

template <typename T, size_t dim>
bool operator!=(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  return a.data() != b.data();
}

template <typename T, size_t dim>
typename Vector<T, dim>::Relationship FindRelationShip(const Vector<T, dim>& a, const Vector<T, dim>& b) {
  T dot = a * b;
  if (Comparator<T>::Equal(dot * dot, a.SquaredLength() * b.SquaredLength())) {
    return Vector<T, dim>::Relationship::Parallel;
  }

  if (Comparator<T>::Equal(dot, 0)) {
    return Vector<T, dim>::Relationship::Orthogonal;
  }

  return Vector<T, dim>::Relationship::None;
}

template <typename T, size_t dim>
std::ostream& operator<<(std::ostream& out, const Vector<T, dim>& v) {
  out << '(' << v.data().cords[0];
  for (size_t i = 1; i < dim; ++i) {
    out << ", " << v[i];
  }
  out << ')';
  return out;
}

template <typename T, size_t dim>
std::istream& operator>>(std::istream& in, Vector<T, dim>& v) {
  for (size_t i = 0; i < dim; ++i) {
    in >> v[i];
  }
  return in;
}

#endif //GEOMERTY_GEOMETRY_VECTOR_H_