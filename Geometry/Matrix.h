//
// Created by Artem Novikov on 30.06.2022.
//

#include <iomanip>
#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_MATRIX_H_
#define GEOMERTY_GEOMETRY_MATRIX_H_

template <typename T, size_t N, size_t M = N>
class Matrix {
 public:
  /// construction
  Matrix();
  Matrix(const Matrix& other) = default;
  Matrix(Matrix&& other) noexcept = default;
  Matrix(std::initializer_list<Vector<T, N>> list);
  template <typename... Args, typename = std::enable_if_t<size_v<Args...> == M>>
  Matrix(Args&& ... args);
  ~Matrix() = default;
  Matrix& operator=(const Matrix& other) = default;
  Matrix& operator=(Matrix&& other) noexcept = default;

  /// getters and setters
  Vector<T, N>& operator[](size_t column) { return data_[column]; }
  const Vector<T, N>& operator[](size_t column) const { return data_[column]; }

  /// calc
  void Triangulate();
  void Diagonalize();
  template <typename = std::enable_if_t<N == M>>
  void Transpose();
  size_t Rank() const;
  template <typename = std::enable_if_t<N == M>>
  T Determinant() const;
  template <typename = std::enable_if_t<N == M>>
  T Trace() const;
  template <typename = std::enable_if_t<N == M>>
  void Invert();

  /// arithmetic
  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  template <typename = std::enable_if_t<N == M>>
  Matrix& operator*=(const Matrix<T, N, M>& other);
  Matrix& operator*=(const T& scalar);
  template <typename = std::enable_if_t<N == M>>
  Matrix& operator+=(const T& scalar);
  template <typename = std::enable_if_t<N == M>>
  Matrix& operator-=(const T& scalar);
 private:
  template <typename U>
  void PushVector(U&& value);
  template <typename U, typename... Args>
  void PushVector(U&& value, Args&& ...args);

  bool IsZeroRow(size_t row) const;
  void ScaleRow(size_t row, T scale);
  void AddRow(size_t row, size_t add, T scale);
  void SwapRows(size_t row1, size_t row2);

  std::array<Vector<T, N>, M> data_;
};

/// calc

template <typename T, size_t N, size_t M>
Matrix<T, M, N> Transposed(const Matrix<T, N, M>& matrix);

template <typename T, size_t N>
Matrix<T, N> Triangulated(const Matrix<T, N>& matrix);

template <typename T, size_t N>
Matrix<T, N> Diagonalized(const Matrix<T, N>& matrix);

template <typename T, size_t N>
Matrix<T, N> Inverted(const Matrix<T, N>& matrix);

/// arithmetic

template <typename T, size_t N, size_t M>
Matrix<T, N, M> operator+(const Matrix<T, N, M>& a, const Matrix<T, N, M>& b);

template <typename T, size_t N, size_t M>
Matrix<T, N, M> operator-(const Matrix<T, N, M>& a, const Matrix<T, N, M>& b);

template <typename T, size_t N, size_t M, size_t K>
Matrix<T, N, K> operator*(const Matrix<T, N, M>& a, const Matrix<T, M, K>& b);

template <typename T, size_t N, size_t M, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Matrix<T, N, M> operator*(const Matrix<T, N, M>& matrix, const U& scalar);

template <typename T, size_t N, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Matrix<T, N> operator+(const Matrix<T, N>& matrix, const U& scalar);

template <typename T, size_t N, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Matrix<T, N> operator-(const Matrix<T, N>& matrix, const U& scalar);

/// stream

template <typename T, size_t N, size_t M>
std::ostream& operator<<(std::ostream& out, const Matrix<T, N, M>& m);

template <typename T, size_t N, size_t M>
std::istream& operator>>(std::istream& in, Matrix<T, N, M>& m);

///////////////////////////////////////////////////////DEFINITION///////////////////////////////////////////////////////

template <typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix() {
  if constexpr(N == M) {
    for (size_t i = 0; i < N; ++i) {
      data_[i][i] = 1;
    }
  }
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>::Matrix(std::initializer_list<Vector<T, N>> list) {
  size_t index = 0;
  for (auto& el : list) {
    data_[index] = el;
    ++index;
  }
}

template <typename T, size_t N, size_t M>
template <typename... Args, typename>
Matrix<T, N, M>::Matrix(Args&& ... args) {
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t N, size_t M>
void Matrix<T, N, M>::Triangulate() {
  for (size_t row = 0; row < std::min(N, M); ++row) {
    if (Comparator<T>::IsZero(data_[row][row])) {
      data_[row][row] = 0; /// To make it zero for sure
      for (size_t non_zero_row = row + 1; non_zero_row < N; ++non_zero_row) {
        if (Comparator<T>::IsZero(data_[row][non_zero_row])) {
          SwapRows(row, non_zero_row);
          break;
        }
      }
    }
    if (!Comparator<T>::IsZero(data_[row][row])) {
      for (size_t i = row + 1; i < N; ++i) {
        AddRow(i, row, -data_[row][i] / data_[row][row]);
        data_[row][i] = 0;
      }
    }
  }
}

template <typename T, size_t N, size_t M>
void Matrix<T, N, M>::Diagonalize() {
  Triangulate();
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i + 1; j < M; ++j) {
      data_[i][j] = 0;
    }
  }
}

template <typename T, size_t N, size_t M>
template <typename>
void Matrix<T, N, M>::Transpose() {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(data_[i][j], data_[j][i]);
    }
  }
}

template <typename T, size_t N, size_t M>
size_t Matrix<T, N, M>::Rank() const {
  auto tri = Triangulated(*this);
  size_t non_trivial_rows = 0;
  while (tri.IsZeroRow(non_trivial_rows)) {
    ++non_trivial_rows;
  }
  return non_trivial_rows;
}

template <typename T, size_t N, size_t M>
template <typename>
T Matrix<T, N, M>::Determinant() const {
  auto tri = Triangulated(*this);
  T det = 1;
  for (size_t i = 0; i < N; ++i) {
    det *= tri[i][i];
  }
  return det;
}

template <typename T, size_t N, size_t M>
template <typename>
T Matrix<T, N, M>::Trace() const {
  T sum = 0;
  for (size_t i = 0; i < N; ++i) {
    sum += data_[i][i];
  }
  return sum;
}

template <typename T, size_t N, size_t M>
template <typename>
void Matrix<T, N, M>::Invert() {
  Matrix<T, N> identity_matrix;
  for (size_t row = 0; row < std::min(N, M); ++row) {
    if (Comparator<T>::IsZero(data_[row][row])) {
      data_[row][row] = 0; /// To make it zero for sure
      for (size_t non_zero_row = row + 1; non_zero_row < N; ++non_zero_row) {
        if (Comparator<T>::IsZero(data_[row][non_zero_row])) {
          SwapRows(row, non_zero_row);
          identity_matrix.SwapRows(row, non_zero_row);
          break;
        }
      }
    }
    if (!Comparator<T>::IsZero(data_[row][row])) {
      for (size_t i = 0; i < N; ++i) {
        if (i != row) {
          T scale = -data_[row][i] / data_[row][row];
          identity_matrix.AddRow(i, row, scale);
          AddRow(i, row, scale);
          data_[row][i] = 0;
        }
      }
      identity_matrix.ScaleRow(row, data_[row][row]);
      ScaleRow(row, 1 / data_[row][row]);
    }
  }
  *this = std::move(identity_matrix);
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator+=(const Matrix& other) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      data_[i][j] += other[i][j];
    }
  }
  return *this;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator-=(const Matrix& other) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      data_[i][j] -= other[i][j];
    }
  }
  return *this;
}

template <typename T, size_t N, size_t M>
template <typename>
Matrix<T, N, M>& Matrix<T, N, M>::operator*=(const Matrix<T, N, M>& other) {
  auto prev = std::move(*this);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      data_[i][j] = 0;
      for (size_t k = 0; k < N; ++k) {
        data_[i][j] += prev[i][k] * other[k][j];
      }
    }
  }
  return *this;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator*=(const T& scalar) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      data_[i][j] *= scalar;
    }
  }
  return *this;
}

template <typename T, size_t N, size_t M>
template <typename>
Matrix<T, N, M>& Matrix<T, N, M>::operator+=(const T& scalar) {
  for (size_t i = 0; i < N; ++i) {
    data_[i][i] += scalar;
  }
  return *this;
}

template <typename T, size_t N, size_t M>
template <typename>
Matrix<T, N, M>& Matrix<T, N, M>::operator-=(const T& scalar) {
  for (size_t i = 0; i < N; ++i) {
    data_[i][i] -= scalar;
  }
  return *this;
}

template <typename T, size_t N, size_t M>
Matrix<T, M, N> Transposed(const Matrix<T, N, M>& matrix) {
  Matrix<T, M, N> matrix_t;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      matrix_t[j][i] = matrix[i][j];
    }
  }
  return matrix_t;
}

template <typename T, size_t N>
Matrix<T, N> Triangulated(const Matrix<T, N>& matrix) {
  auto tri = matrix;
  tri.Triangulate();
  return tri;
}

template <typename T, size_t N>
Matrix<T, N> Diagonalized(const Matrix<T, N>& matrix) {
  auto diag = matrix;
  diag.Diagonalize();
  return diag;
}

template <typename T, size_t N>
Matrix<T, N> Inverted(const Matrix<T, N>& matrix) {
  auto inv = matrix;
  inv.Invert();
  return inv;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M> operator+(const Matrix<T, N, M>& a, const Matrix<T, N, M>& b) {
  auto sum = a;
  sum += b;
  return sum;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M> operator-(const Matrix<T, N, M>& a, const Matrix<T, N, M>& b) {
  auto sum = a;
  sum -= b;
  return sum;
}

template <typename T, size_t N, size_t M, size_t K>
Matrix<T, N, K> operator*(const Matrix<T, N, M>& a, const Matrix<T, M, K>& b) {
  Matrix<T, N, K> mult;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < K; ++j) {
      mult[i][j] = 0;
      for (size_t k = 0; k < M; ++k) {
        mult[i][j] += a[i][k] + b[k][j];
      }
    }
  }
  return mult;
}

template <typename T, size_t N, size_t M, typename U, typename>
Matrix<T, N, M> operator*(const Matrix<T, N, M>& matrix, const U& scalar) {
  auto result = matrix;
  result *= scalar;
  return result;
}

template <typename T, size_t N, typename U, typename>
Matrix<T, N> operator+(const Matrix<T, N>& matrix, const U& scalar) {
  auto result = matrix;
  result += scalar;
  return result;
}

template <typename T, size_t N, typename U, typename>
Matrix<T, N> operator-(const Matrix<T, N>& matrix, const U& scalar) {
  auto result = matrix;
  result -= scalar;
  return result;
}

template <typename T, size_t N, size_t M>
template <typename U>
void Matrix<T, N, M>::PushVector(U&& value) {
  data_[M - 1] = value;
}

template <typename T, size_t N, size_t M>
template <typename U, typename... Args>
void Matrix<T, N, M>::PushVector(U&& value, Args&& ... args) {
  data_[M - sizeof...(Args) - 1] = value;
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t N, size_t M>
bool Matrix<T, N, M>::IsZeroRow(size_t row) const {
  bool ans = true;
  for (size_t j = 0; j < M; ++j) {
    ans &= Comparator<T>::IsZero(data_[j][row]);
  }
  return ans;
}

template <typename T, size_t N, size_t M>
void Matrix<T, N, M>::ScaleRow(size_t row, T scale) {
  for (size_t j = 0; j < M; ++j) {
    data_[j][row] *= scale;
  }
}

template <typename T, size_t N, size_t M>
void Matrix<T, N, M>::AddRow(size_t row, size_t add, T scale) {
  for (size_t j = 0; j < M; ++j) {
    data_[j][row] += data_[j][add] * scale;
  }
}

template <typename T, size_t N, size_t M>
void Matrix<T, N, M>::SwapRows(size_t row1, size_t row2) {
  for (size_t j = 0; j < M; ++j) {
    std::swap(data_[j][row1], data_[j][row2]);
  }
}

template <typename T, size_t N, size_t M>
std::ostream& operator<<(std::ostream& out, const Matrix<T, N, M>& m) {
  for (size_t i = 0; i < N; ++i) {
    out << "| ";
    for (size_t j = 0; j < M; ++j) {
      out << std::setprecision(2) << m[i][j] << ' ';
    }
    out << "|\n";
  }
  return out;
}

template <typename T, size_t N, size_t M>
std::istream& operator>>(std::istream& in, Matrix<T, N, M>& m){
  for(size_t i = 0; i < N; ++i){
    for(size_t j = 0; j < M; ++j){
      in >> m[i][j];
    }
  }
  return in;
}

#endif //GEOMERTY_GEOMETRY_MATRIX_H_
