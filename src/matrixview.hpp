#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <cmath>
#include <iostream>
#include <matrixexpr.hpp>
#include <vector.hpp>
#include "../myASC-HPC/src/simd.hpp"


namespace ASC_bla {
  enum ORDERING { ColMajor, RowMajor };

  template<typename T, ORDERING ORD>
  class MatrixView : public MatrixExpr<MatrixView<T, ORD>> {
  protected:
    size_t m_rows, m_cols, m_dist;
    T *m_data;

  public:
    MatrixView(const MatrixView &) = default;

    MatrixView(size_t rows, size_t cols, size_t dist, T *data)
        : m_data(data), m_rows(rows), m_cols(cols), m_dist(dist) {}

    MatrixView(size_t rows, size_t cols, T *data) : MatrixView(rows, cols, ORD == RowMajor ? cols : rows, data) {}

    auto row(size_t i) const {
      assert(i < m_rows);
      if constexpr (ORD == RowMajor) {
        return VectorView<T>(m_cols, &m_data[i * m_dist]);
      }
      if constexpr (ORD == ColMajor) {
        return VectorView<T, size_t>(m_cols, m_dist, &m_data[i]);
      }
    }

    auto col(size_t i) const {
      assert(i < m_cols);
      if constexpr (ORD == RowMajor) {
        return VectorView<T, size_t>(m_rows, m_dist, &m_data[i]);
      }
      if constexpr (ORD == ColMajor) {
        return VectorView<T>(m_cols, &m_data[i * m_dist]);
      }
    }

    auto rows(size_t first, size_t next) const {
      assert(first <= next && next <= m_rows);
      if constexpr (ORD == RowMajor) {
        return MatrixView<T, ORD>(next - first, m_cols, &m_data[first * m_dist]);
      }
      if constexpr (ORD == ColMajor) {
        return MatrixView<T, ORD>(next - first, m_cols, m_dist, &m_data[first]);
      }
    }

    auto cols(size_t first, size_t next) const {
      assert(first <= next && next <= m_cols);
      if constexpr (ORD == RowMajor) {
        return MatrixView<T, ORD>(m_rows, next - first, m_dist, &m_data[first]);
      }
      if constexpr (ORD == ColMajor) {
        return MatrixView<T, ORD>(m_rows, next - first, &m_data[first * m_dist]);
      }
    }

    auto transpose() const {
      if constexpr (ORD == RowMajor) {
        return MatrixView<T, ColMajor>(m_cols, m_rows, m_data);
      }
      if constexpr (ORD == ColMajor) {
        return MatrixView<T, RowMajor>(m_cols, m_rows, m_data);
      }
    }


    MatrixView &operator=(const MatrixView<T, ORD> &m2) {
      assert(m_rows == m2.rows());
      assert(m_cols == m2.cols());
      for (int i = 0; i < rows()/2; i+=2) {
        for (int j = 0; j < cols()/2; j+=2) {
          (*this)(i, j) = m2(i, j);
          (*this)(i, j+1) = m2(i, j+1);
          (*this)(i+1, j) = m2(i+1, j);
          (*this)(i+1, j+1) = m2(i+1, j+1);
        }
      }
      // residual indices 
      for (int i = (cols()/2)*2; i < rows(); i++) {
        for (int j = (cols()/2)*2; j < cols(); j++) {
          (*this)(i, j) = m2(i, j);
        }
      }
      return *this;
    }

    template<typename TB>
    MatrixView &operator=(const MatrixExpr<TB> &m2) {
      assert(m_rows == m2.rows());
      assert(m_cols == m2.cols());
      for (int i = 0; i < rows()/2; i+=2) {
        for (int j = 0; j < cols()/2; j+=2) {
          (*this)(i, j) = m2(i, j);
          (*this)(i, j+1) = m2(i, j+1);
          (*this)(i+1, j) = m2(i+1, j);
          (*this)(i+1, j+1) = m2(i+1, j+1);
        }
      }
      // residual indices 
      for (int i = (cols()/2)*2; i < rows(); i++) {
        for (int j = (cols()/2)*2; j < cols(); j++) {
          (*this)(i, j) = m2(i, j);
        }
      }


      return *this;
    }

    MatrixView &operator=(T scal) {
      for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < cols(); j++) {
          *this(i, j) = scal;
        }
      }
      return *this;
    }

    size_t rows() const { return m_rows; }

    size_t cols() const { return m_cols; }

    size_t dist() const { return m_dist;}

    T *data() const { return m_data; }


    T &operator()(size_t i, size_t j) {
      if constexpr (ORD == RowMajor) {
        return m_data[i * m_dist + j];
      }
      if constexpr (ORD == ColMajor) {
        return m_data[i + m_dist * j];
      }
    }

    const T &operator()(size_t i, size_t j) const {
      if constexpr (ORD == RowMajor) {
        return m_data[i * m_dist + j];
      }
      if constexpr (ORD == ColMajor) {
        return m_data[i + m_dist * j];
      }
    }
  };

  template<typename T, ORDERING ORD>
  class Matrix : public MatrixView<T, ORD> {
    typedef MatrixView<T, ORD> BASE;
    using BASE::m_rows;
    using BASE::m_cols;
    using BASE::m_dist;
    using BASE::m_data;
  public:
    Matrix(size_t rows, size_t cols)
        : MatrixView<T, ORD>(rows, cols, new T[rows * cols]) {}

    Matrix(const Matrix &other)
        : Matrix(other.rows(), other.cols()) {
      *this = other;
    }

    template<typename TB>
    Matrix(const MatrixExpr<TB> &m2) : Matrix(m2.rows(), m2.cols()) {
      *this = m2;
    }

    Matrix(Matrix &&other) noexcept
        : MatrixView<T, ORD>(0, 0, 0, nullptr) {
      std::swap(m_data, other.m_data);
      std::swap(m_rows, other.m_rows);
      std::swap(m_cols, other.m_cols);
      std::swap(m_dist, other.m_dist);
    }

    // assign from any MatrixExpr (deep-copy element-wise)
    template<typename TB>
    Matrix &operator=(const MatrixExpr<TB> &m2) {
      MatrixView<T, ORD>::operator=(m2);
      return *this;
    }

    // copy assign (deep-copy)
    Matrix &operator=(const Matrix &other) {
      assert(m_rows == other.m_rows);
      assert(m_cols == other.m_cols);
      *this = other;
      return *this;
    }

    // move assign (steal ownership)
    Matrix &operator=(Matrix &&other) noexcept {
      m_rows = 0;
      m_cols = 0;
      m_dist = 0;
      m_data = nullptr;
      std::swap(m_data, other.m_data);
      std::swap(m_rows, other.m_rows);
      std::swap(m_cols, other.m_cols);
      std::swap(m_dist, other.m_dist);
      return *this;
    }

    ~Matrix() {
      if (m_data) {
        delete[] this->m_data;
        m_data = nullptr;
      }
    }

    T *getRawDataDanger() {
      return m_data;
    }

    Matrix<T, ORD> inverse() const {
      if (m_rows != m_cols) {
        throw std::invalid_argument("Error: Matrix must be square for inversion");
      }
      size_t n = m_rows;
      Matrix<T, ORD> aug(n, 2 * n);

      // Copy A into left half and identity into right half
      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
          aug(i, j) = (*this)(i, j);
          aug(i, j + n) = (i == j) ? 1.0 : 0.0;
        }
      }

      // Gauss-Jordan elimination
      for (size_t pivot = 0; pivot < n; pivot++) {
        // Find the row with largest pivot element (partial pivoting)
        size_t max_row = pivot;
        T max_val = std::abs(aug(pivot, pivot));
        for (size_t i = pivot + 1; i < n; i++) {
          T val = std::abs(aug(i, pivot));
          if (val > max_val) {
            max_val = val;
            max_row = i;
          }
        }

        if (max_val < 1e-14) {
          throw std::runtime_error("Error: Matrix is singular, cannot invert");
        }

        // Swap rows if needed
        if (max_row != pivot) {
          for (size_t j = 0; j < 2 * n; j++) {
            std::swap(aug(pivot, j), aug(max_row, j));
          }
        }

        // Scale pivot row to make diagonal element = 1
        T pivot_val = aug(pivot, pivot);
        auto pivot_row = aug.row(pivot);
        pivot_row = (1.0 / pivot_val) * pivot_row;

        // Eliminate column in all other rows
        for (size_t i = 0; i < n; i++) {
          if (i != pivot) {
            T factor = aug(i, pivot);
            auto current_row = aug.row(i);
            current_row = current_row - factor * pivot_row;
          }
        }
      }

      // Extract the inverse from the right half
      Matrix<T, ORD> inv(n, n);
      inv = aug.cols(n, 2 * n);

      return inv;
    }
  };

  template<typename T>
  std::ostream &operator<<(std::ostream &ost, const MatrixExpr<T> &M) {
    for (size_t i = 0; i < M.rows(); i++) {
      if (M.cols() > 0) {
        ost << M(i, 0);
      }
      for (size_t j = 1; j < M.cols(); j++) {
        ost << ", " << M(i, j);
      }
      ost << "\n";
    }

    return ost;
  }
};
#endif