//
// Created by Artem Novikov on 07.07.2022.
//

#include "Vector.h"
#include "Matrix.h"

#ifndef GEOMETRY_GEOMETRY_TRANSFORM_H_
#define GEOMETRY_GEOMETRY_TRANSFORM_H_

template <typename T, size_t InputDimension, size_t OutputDimension = InputDimension>
class Transform {
 public:
  /// construction
  Transform(const Transform& other) = default;
  Transform(Transform&& other) noexcept = default;
  Transform(const Matrix<T, InputDimension, OutputDimension>& matrix,
            const Vector<T, OutputDimension>& shift = Vector<T, OutputDimension>());
  ~Transform() = default;
  Transform& operator=(const Transform& other) = default;
  Transform& operator=(Transform&& other) noexcept = default;

  /// calls
  template <typename Type>
  decltype(auto) operator()(const Type& value) const;
  Vector<T, OutputDimension> operator()(const Vector<T, InputDimension>& value) const;

  /// arithmetic
  Transform& operator+=(const Transform<T, InputDimension, OutputDimension>& other);
  Transform& operator-=(const Transform<T, InputDimension, OutputDimension>& other);
  Transform& operator*=(const Transform<T, InputDimension, InputDimension>& other);

  /// getters and setters
  const Matrix<T, InputDimension, OutputDimension>& GetMatrix() const { return matrix_; }
  Matrix<T, InputDimension, OutputDimension>& GetMatrix() { return matrix_; }
  const Vector<T, OutputDimension>& GetShift() const { return shift_; }
  Vector<T, OutputDimension>& GetShift() { return shift_; }

  /// calc
  bool Affine() const;
  bool Orthogonal() const;
  template <bool Access = InputDimension == OutputDimension, typename = std::enable_if_t<Access>>
  Transform& Invert();
 private:
  Matrix<T, InputDimension, OutputDimension> matrix_;
  Vector<T, OutputDimension> shift_;
};

template <typename T, size_t InputDimension, size_t TransitDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator*(const Transform<T, TransitDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, TransitDimension>& b);

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator+(const Transform<T, InputDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, OutputDimension>& b);

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator-(const Transform<T, InputDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, OutputDimension>& b);

template <typename T, size_t Dimension>
Transform<T, Dimension> Inverted(const Transform<T, Dimension>& transform);

template <typename T, size_t Dimension>
Transform<T, Dimension> Shift(const Vector<T, Dimension>& shift);

template <typename T, size_t Dimension>
Transform<T, Dimension> Flip(size_t direction);

template <typename T, size_t Dimension>
Transform<T, Dimension> Rotate(size_t abscissa, size_t ordinate, const T& angle);

//////////////////////////////////////////////////////////DEFINITION////////////////////////////////////////////////////

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension>::Transform(const Matrix<T, InputDimension, OutputDimension>& matrix,
                                                         const Vector<T, OutputDimension>& shift)
    : matrix_(matrix), shift_(shift) {}

template <typename T, size_t InputDimension, size_t OutputDimension>
Vector<T, OutputDimension> Transform<T, InputDimension, OutputDimension>::operator()(
    const Vector<T, InputDimension>& value) const {
  return matrix_ * value;
}

template <typename T, size_t InputDimension, size_t OutputDimension>
template <typename Type>
decltype(auto) Transform<T, InputDimension, OutputDimension>::operator()(const Type& value) const {
  return value.ApplyTransform(*this);
}

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension>& Transform<T, InputDimension, OutputDimension>::operator+=(
    const Transform<T, InputDimension, OutputDimension>& other) {
  matrix_ += other.matrix_;
  shift_ += other.shift_;
  return *this;
}

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension>& Transform<T, InputDimension, OutputDimension>::operator-=(
    const Transform<T, InputDimension, OutputDimension>& other) {
  matrix_ -= other.matrix_;
  shift_ -= other.shift_;
  return *this;
}

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension>& Transform<T, InputDimension, OutputDimension>::operator*=(
    const Transform<T, InputDimension, InputDimension>& other) {
  shift_ += matrix_ * other.shift_;
  matrix_ *= other.matrix_;
  return *this;
}

template <typename T, size_t InputDimension, size_t OutputDimension>
bool Transform<T, InputDimension, OutputDimension>::Affine() const {
  if constexpr (InputDimension == OutputDimension) {
    return !Comparator<T>::IsZero(matrix_.Determinant());
  } else {
    return false;
  }
}

template <typename T, size_t InputDimension, size_t OutputDimension>
bool Transform<T, InputDimension, OutputDimension>::Orthogonal() const {
  return false; //// TODO
}

template <typename T, size_t InputDimension, size_t OutputDimension>
template <bool, typename>
Transform<T, InputDimension, OutputDimension>& Transform<T, InputDimension, OutputDimension>::Invert() {
  matrix_.Invert();
  shift_ = -matrix_ * shift_;
  return *this;
}

template <typename T, size_t InputDimension, size_t TransitDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator*(const Transform<T, TransitDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, TransitDimension>& b) {
  return {a.GetMatrix() * b.GetMatrix(), a.GetShift() + a.GetMatrix() * b.GetShift()};
}

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator+(const Transform<T, InputDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, OutputDimension>& b) {
  auto copy = a;
  copy += b;
  return copy;
}

template <typename T, size_t InputDimension, size_t OutputDimension>
Transform<T, InputDimension, OutputDimension> operator-(const Transform<T, InputDimension, OutputDimension>& a,
                                                        const Transform<T, InputDimension, OutputDimension>& b) {
  auto copy = a;
  copy -= b;
  return copy;
}

template <typename T, size_t Dimension>
Transform<T, Dimension> Inverted(const Transform<T, Dimension>& transform) {
  auto copy = transform;
  copy.Invert();
  return copy;
}

template <typename T, size_t Dimension>
Transform<T, Dimension> Shift(const Vector<T, Dimension>& shift) {
  return Transform(Identity<Dimension>(), shift);
}

template <typename T, size_t Dimension>
Transform<T, Dimension> Flip(size_t direction) {
  auto matrix = Identity<Dimension>();
  matrix[direction][direction] *= -1;
  return Transform(matrix, Vector<T, Dimension>());
}

template <typename T, size_t Dimension>
Transform<T, Dimension> Rotate(size_t abscissa, size_t ordinate, const T& angle) {
  auto matrix = Identity<Dimension>();
  matrix[abscissa][abscissa] = std::cos(angle);
  matrix[ordinate][abscissa] = -std::sin(angle);
  matrix[ordinate][ordinate] = matrix[abscissa][abscissa];
  matrix[abscissa][ordinate] = -matrix[ordinate][abscissa];
  return Transform(matrix, Vector<T, Dimension>());
}

#endif //GEOMETRY_GEOMETRY_TRANSFORM_H_
