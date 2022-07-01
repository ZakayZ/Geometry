//
// Created by Artem Novikov on 30.06.2022.
//

#include <iomanip>
#include "Vector.h"

#ifndef GEOMERTY_GEOMETRY_MATRIX_H_
#define GEOMERTY_GEOMETRY_MATRIX_H_

template <typename T, size_t Rows, size_t Columns = Rows>
class Matrix {
 public:
  /// construction
  Matrix();
  Matrix(const Matrix& other) = default;
  Matrix(Matrix&& other) noexcept = default;
  Matrix(std::initializer_list<Vector<T, Rows>> list);
  template <typename... Args, bool Temp = sizeof...(Args) == Columns, typename = std::enable_if_t<Temp>>
  Matrix(Args&& ... args);
  ~Matrix() = default;
  Matrix& operator=(const Matrix& other) = default;
  Matrix& operator=(Matrix&& other) noexcept = default;

  /// getters and setters
  Vector<T, Rows>& operator[](size_t column) { return data_[column]; }
  const Vector<T, Rows>& operator[](size_t column) const { return data_[column]; }

  /// calc
  void Triangulate();
  void Diagonalize();
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  void Transpose();
  size_t Rank() const;
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  T Determinant() const;
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  T Trace() const;
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  void Invert();

  /// arithmetic
  Matrix& operator+=(const Matrix& other);
  Matrix& operator-=(const Matrix& other);
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  Matrix& operator*=(const Matrix<T, Columns, Rows>& other);
  Matrix& operator*=(const T& scalar);
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  Matrix& operator+=(const T& scalar);
  template <bool Temp = Rows == Columns, typename = std::enable_if_t<Temp>>
  Matrix& operator-=(const T& scalar);
 private:
  template <typename U>
  void PushVector(U&& value);
  template <typename U, typename... Args>
  void PushVector(U&& value, Args&& ...args);

  [[nodiscard]] bool IsZeroRow(size_t row) const;
  void ScaleRow(size_t row, T scale);
  void AddRow(size_t row, size_t add, T scale);
  void SwapRows(size_t row1, size_t row2);

  std::array<Vector<T, Rows>, Columns> data_;
};

/// calc

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Columns, Rows> Transposed(const Matrix<T, Rows, Columns>& matrix);

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Triangulated(const Matrix<T, Rows, Columns>& matrix);

template <typename T, size_t N>
Matrix<T, N> Diagonalized(const Matrix<T, N>& matrix);

template <typename T, size_t N>
Matrix<T, N> Inverted(const Matrix<T, N>& matrix);

/// arithmetic

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator+(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b);

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator-(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b);

template <typename T, size_t Rows, size_t Columns, size_t OtherColumns>
Matrix<T, Rows, OtherColumns> operator*(const Matrix<T, Rows, Columns>& a, const Matrix<T, Columns, OtherColumns>& b);

template <typename T, size_t Rows, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T,
                                                                                                                   U>>>
Matrix<T, Rows, Columns> operator*(const Matrix<T, Rows, Columns>& matrix, const U& scalar);

template <typename T, size_t Rows, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T,
                                                                                                                   U>>>
Matrix<T, Rows, Columns> operator*(const U& scalar, const Matrix<T, Rows, Columns>& matrix);

template <typename T, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Matrix<T, Columns> operator+(const Matrix<T, Columns>& matrix, const U& scalar);

template <typename T, size_t Rows, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T,
                                                                                                                   U>>>
Matrix<T, Rows, Columns> operator+(const U& scalar, const Matrix<T, Rows, Columns>& matrix);

template <typename T, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T, U>>>
Matrix<T, Columns> operator-(const Matrix<T, Columns>& matrix, const U& scalar);

template <typename T, size_t Rows, size_t Columns, typename U, typename = std::enable_if_t<std::is_constructible_v<T,
                                                                                                                   U>>>
Matrix<T, Rows, Columns> operator-(const U& scalar, const Matrix<T, Rows, Columns>& matrix);

/// relationship

template <typename T, size_t Rows, size_t Columns>
bool operator==(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b);

template <typename T, size_t Rows, size_t Columns>
bool operator!=(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b);

/// stream

template <typename T, size_t Rows, size_t Columns>
std::ostream& operator<<(std::ostream& out, const Matrix<T, Rows, Columns>& m);

template <typename T, size_t Rows, size_t Columns>
std::istream& operator>>(std::istream& in, Matrix<T, Rows, Columns>& m);

///////////////////////////////////////////////////////DEFINITION///////////////////////////////////////////////////////

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns>::Matrix() {
  if constexpr(Columns == Rows) {
    for (size_t i = 0; i < Columns; ++i) {
      for (size_t j = 0; j < Columns; ++j) {
        data_[i][j] = 0;
      }
      data_[i][i] = 1;
    }
  }
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns>::Matrix(std::initializer_list<Vector<T, Rows>> list) {
  size_t index = 0;
  for (auto& el : list) {
    data_[index] = el;
    ++index;
  }
}

template <typename T, size_t Rows, size_t Columns>
template <typename... Args, bool, typename>
Matrix<T, Rows, Columns>::Matrix(Args&& ... args) {
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Rows, size_t Columns>
void Matrix<T, Rows, Columns>::Triangulate() {
  for (size_t row = 0; row < std::min(Rows, Columns); ++row) {
    if (Comparator<T>::IsZero(data_[row][row])) {
      data_[row][row] = 0;
      for (size_t non_zero_row = row + 1; non_zero_row < Rows; ++non_zero_row) {
        if (Comparator<T>::IsZero(data_[row][non_zero_row])) {
          SwapRows(row, non_zero_row);
          break;
        }
      }
    }
    if (!Comparator<T>::IsZero(data_[row][row])) {
      for (size_t other_row = row + 1; other_row < Rows; ++other_row) {
        AddRow(other_row, row, -data_[row][other_row] / data_[row][row]);
        data_[row][other_row] = 0;
      }
    }
  }
}

template <typename T, size_t Rows, size_t Columns>
void Matrix<T, Rows, Columns>::Diagonalize() {
  Triangulate();
  for (size_t row = 0; row < Rows; ++row) {
    if (!Comparator<T>::IsZero(data_[row][row])) {
      for (size_t other_row = 0; other_row < row; ++other_row) {
        AddRow(other_row, row, -data_[row][other_row] / data_[row][row]);
        data_[row][other_row] = 0;
      }
    }
  }
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
void Matrix<T, Rows, Columns>::Transpose() {
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < i; ++j) {
      std::swap(data_[i][j], data_[j][i]);
    }
  }
}

template <typename T, size_t Rows, size_t Columns>
size_t Matrix<T, Rows, Columns>::Rank() const {
  auto tri = Triangulated(*this);
  size_t non_trivial_rows = 0;
  while (!Comparator<T>::IsZero(tri[non_trivial_rows][non_trivial_rows])
      && non_trivial_rows < std::min(Rows, Columns)) {
    ++non_trivial_rows;
  }
  return non_trivial_rows;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
T Matrix<T, Rows, Columns>::Determinant() const {
  auto tri = Triangulated(*this);
  T det = 1;
  for (size_t i = 0; i < Columns; ++i) {
    det *= tri[i][i];
  }
  return det;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
T Matrix<T, Rows, Columns>::Trace() const {
  T sum = 0;
  for (size_t i = 0; i < Columns; ++i) {
    sum += data_[i][i];
  }
  return sum;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
void Matrix<T, Rows, Columns>::Invert() {
  Matrix<T, Columns> identity_matrix;
  for (size_t row = 0; row < std::min(Rows, Columns); ++row) {
    if (Comparator<T>::IsZero(data_[row][row])) {
      data_[row][row] = 0;
      for (size_t non_zero_row = row + 1; non_zero_row < Rows; ++non_zero_row) {
        if (Comparator<T>::IsZero(data_[row][non_zero_row])) {
          SwapRows(row, non_zero_row);
          identity_matrix.SwapRows(row, non_zero_row);
          break;
        }
      }
    }
    if (!Comparator<T>::IsZero(data_[row][row])) {
      for (size_t other_row = 0; other_row < Rows; ++other_row) {
        if (other_row != row) {
          T scale = -data_[row][other_row] / data_[row][row];
          identity_matrix.AddRow(other_row, row, scale);
          AddRow(other_row, row, scale);
          data_[row][other_row] = 0;
        }
      }
      T scale = 1 / data_[row][row];
      identity_matrix.ScaleRow(row, scale);
      ScaleRow(row, scale);
    }
  }
  *this = std::move(identity_matrix);
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator+=(const Matrix& other) {
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      data_[i][j] += other[i][j];
    }
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator-=(const Matrix& other) {
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      data_[i][j] -= other[i][j];
    }
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator*=(const Matrix<T, Columns, Rows>& other) {
  const static size_t N = Rows;
  auto prev = std::move(*this);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      data_[i][j] = 0;
      for (size_t k = 0; k < N; ++k) {
        data_[i][j] += prev[k][i] * other[j][k];
      }
    }
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator*=(const T& scalar) {
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      data_[i][j] *= scalar;
    }
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator+=(const T& scalar) {
  for (size_t i = 0; i < Columns; ++i) {
    data_[i][i] += scalar;
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
template <bool, typename>
Matrix<T, Rows, Columns>& Matrix<T, Rows, Columns>::operator-=(const T& scalar) {
  for (size_t i = 0; i < Columns; ++i) {
    data_[i][i] -= scalar;
  }
  return *this;
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Columns, Rows> Transposed(const Matrix<T, Rows, Columns>& matrix) {
  Matrix<T, Columns, Rows> matrix_t;
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      matrix_t[j][i] = matrix[i][j];
    }
  }
  return matrix_t;
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Triangulated(const Matrix<T, Rows, Columns>& matrix) {
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

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator+(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b) {
  auto sum = a;
  sum += b;
  return sum;
}

template <typename T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator-(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b) {
  auto sum = a;
  sum -= b;
  return sum;
}

template <typename T, size_t Rows, size_t Columns, size_t OtherColumns>
Matrix<T, Rows, OtherColumns> operator*(const Matrix<T, Rows, Columns>& a, const Matrix<T, Columns, OtherColumns>& b) {
  Matrix<T, Rows, OtherColumns> mult;
  for (size_t i = 0; i < OtherColumns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      mult[i][j] = 0;
      for (size_t k = 0; k < Columns; ++k) {
        mult[i][j] += b[j][k] * a[k][i];
      }
    }
  }
  return mult;
}

template <typename T, size_t Rows, size_t Columns, typename U, typename>
Matrix<T, Rows, Columns> operator*(const Matrix<T, Rows, Columns>& matrix, const U& scalar) {
  auto result = matrix;
  result *= scalar;
  return result;
}

template <typename T, size_t Rows, size_t Columns, typename U, typename>
Matrix<T, Rows, Columns> operator*(const U& scalar, const Matrix<T, Rows, Columns>& matrix) {
  return matrix * scalar;
}

template <typename T, size_t Columns, typename U, typename>
Matrix<T, Columns> operator+(const Matrix<T, Columns>& matrix, const U& scalar) {
  auto result = matrix;
  result += scalar;
  return result;
}

template <typename T, size_t Rows, size_t Columns, typename U, typename>
Matrix<T, Rows, Columns> operator+(const U& scalar, const Matrix<T, Rows, Columns>& matrix) {
  return matrix + scalar;
}

template <typename T, size_t Columns, typename U, typename>
Matrix<T, Columns> operator-(const Matrix<T, Columns>& matrix, const U& scalar) {
  auto result = matrix;
  result -= scalar;
  return result;
}

template <typename T, size_t Rows, size_t Columns, typename U, typename>
Matrix<T, Rows, Columns> operator-(const U& scalar, const Matrix<T, Rows, Columns>& matrix) {
  return matrix - scalar;
}

template <typename T, size_t Rows, size_t Columns>
template <typename U>
void Matrix<T, Rows, Columns>::PushVector(U&& value) {
  data_[Rows - 1] = value;
}

template <typename T, size_t Rows, size_t Columns>
template <typename U, typename... Args>
void Matrix<T, Rows, Columns>::PushVector(U&& value, Args&& ... args) {
  data_[Rows - sizeof...(Args) - 1] = value;
  PushVector(std::forward<Args>(args)...);
}

template <typename T, size_t Rows, size_t Columns>
bool Matrix<T, Rows, Columns>::IsZeroRow(size_t row) const {
  bool ans = true;
  for (size_t j = 0; j < Columns; ++j) {
    ans &= Comparator<T>::IsZero(data_[j][row]);
  }
  return ans;
}

template <typename T, size_t Rows, size_t Columns>
void Matrix<T, Rows, Columns>::ScaleRow(size_t row, T scale) {
  for (size_t j = 0; j < Columns; ++j) {
    data_[j][row] *= scale;
  }
}

template <typename T, size_t Rows, size_t Columns>
void Matrix<T, Rows, Columns>::AddRow(size_t row, size_t add, T scale) {
  for (size_t j = 0; j < Columns; ++j) {
    data_[j][row] += data_[j][add] * scale;
  }
}

template <typename T, size_t Rows, size_t Columns>
void Matrix<T, Rows, Columns>::SwapRows(size_t row1, size_t row2) {
  for (size_t j = 0; j < Columns; ++j) {
    std::swap(data_[j][row1], data_[j][row2]);
  }
}

template <typename T, size_t Rows, size_t Columns>
bool operator==(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b) {
  for (size_t i = 0; i < Columns; ++i) {
    for (size_t j = 0; j < Rows; ++j) {
      if (!Comparator<T>::Equal(a[i][j], b[i][j])) {
        return false;
      }
    }
  }
  return true;
}

template <typename T, size_t Rows, size_t Columns>
bool operator!=(const Matrix<T, Rows, Columns>& a, const Matrix<T, Rows, Columns>& b) {
  return !(a == b);
}

template <typename T, size_t Rows, size_t Columns>
std::ostream& operator<<(std::ostream& out, const Matrix<T, Rows, Columns>& m) {
  for (size_t i = 0; i < Rows; ++i) {
    out << "| ";
    for (size_t j = 0; j < Columns; ++j) {
      if (Comparator<T>::IsZero(m[j][i])) {
        out << 0 << ' ';
      } else {
        out << m[j][i] << ' ';
      }

    }
    out << "|\n";
  }
  return out;
}

template <typename T, size_t Rows, size_t Columns>
std::istream& operator>>(std::istream& in, Matrix<T, Rows, Columns>& m) {
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      in >> m[j][i];
    }
  }
  return in;
}

#endif //GEOMERTY_GEOMETRY_MATRIX_H_
