#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <cmath>
#include <iostream>
#include <matrixexpr.hpp>
#include <vector.hpp>

namespace ASC_bla {
  enum ORDERING {
    ColMajor, RowMajor
  };

  template<typename T, ORDERING ORD = RowMajor>
  class MatrixView : public MatrixExpr<MatrixView<T, ORD>> {
  protected:
    size_t m_rows, m_cols, m_dist;
    T *m_data;

  public:
    MatrixView(const MatrixView &) = default;

    MatrixView(size_t rows, size_t cols, size_t dist, T *data)
        : m_data(data), m_rows(rows), m_cols(cols), m_dist(dist) {
    }

    MatrixView(size_t rows, size_t cols, T *data) : MatrixView(rows, cols, ORD == RowMajor ? cols : rows, data) {
    }

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

    template<typename TB>
    MatrixView &operator=(const MatrixExpr<TB> &m2) {
      assert(m_rows == m2.rows());
      assert(m_cols == m2.cols());
      for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < cols(); j++) {
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

    using BASE::operator=;

    Matrix(const Matrix &m) : Matrix(m.rows(), m.cols()) { *this = m; }

    Matrix &operator=(Matrix &&m2) {
      std::swap(this->m_rows, m2.m_rows);
      std::swap(this->m_cols, m2.m_cols);
      std::swap(this->m_data, m2.m_data);
      return *this;
    }

    template<typename TB>
    Matrix(const MatrixExpr<TB> &m)
        : Matrix(m.rows(), m.cols()) {
      *this = m;
    }

    ~Matrix() { delete[] this->m_data; }

    T *getRawDataDanger() {
      return m_data;
    }
  };

// template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
// auto operator+(const Matrix<T1, ORD1> &M1, const Matrix<T2, ORD2> &M2)
//     -> Matrix<decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
// {
//   if (M1.Height() != M2.Height() || M1.Width() != M2.Width())
//   {
//     throw std::invalid_argument("Error: Matrix dimensions do not match!\n");
//   }

//   typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES; // common
//   type casting Matrix<TRES, ORD1> sum(M1.Height(), M1.Width()); // use
//   ORDERING of first matrix for (size_t i = 0; i < M1.Height(); i++)
//   {
//     for (size_t j = 0; j < M1.Width(); j++)
//     {
//       sum(i, j) = static_cast<TRES>(M1(i, j)) + static_cast<TRES>(M2(i, j));
//     }
//   }
//   return sum;
// }

// template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
// auto operator-(const Matrix<T1, ORD1> &M1, const Matrix<T2, ORD2> &M2)
//     -> Matrix<decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
// {
//   if (M1.Height() != M2.Height() || M1.Width() != M2.Width())
//   {
//     throw std::invalid_argument("Error: Matrix dimensions do not match!\n");
//   }

//   typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES; // common
//   type casting Matrix<TRES, ORD1> sum(M1.Height(), M1.Width()); // use
//   ORDERING of first matrix for (size_t i = 0; i < M1.Height(); i++)
//   {
//     for (size_t j = 0; j < M1.Width(); j++)
//     {
//       sum(i, j) = static_cast<TRES>(M1(i, j)) - static_cast<TRES>(M2(i, j));
//     }
//   }
//   return sum;
// }

// template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
// auto operator*(const Matrix<T1, ORD1> &M1, const Matrix<T2, ORD2> &M2)
//     -> Matrix<decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
// {
//   if (M1.Width() != M2.Height())
//   {
//     throw std::invalid_argument("Error: Matrix dimensions do not match!\n");
//   }

//   typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;
//   Matrix<TRES, ORD1> prod(M1.Height(), M2.Width());
//   for (size_t i = 0; i < M1.Height(); ++i)
//   {
//     for (size_t j = 0; j < M2.Width(); ++j)
//     {
//       TRES sum = 0; // T gives type of matrix entries
//       for (size_t k = 0; k < M1.Width(); ++k)
//       {
//         sum += M1(i, k) * M2(k, j);
//       }
//       prod(i, j) = sum;
//     }
//   }

//   return prod;
// }

// template <typename T1, typename T2, ORDERING ORD1>
// auto operator*(const Matrix<T1, ORD1> &M, const Vector<T2> &v)
//     -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
// {
//   if (M.Width() != v.Size())
//   {
//     throw std::invalid_argument("Error: Matrix Width must match vector size!
//     \n");
//   }

//   typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;
//   Vector<TRES> prod(M.Height());
//   for (size_t i = 0; i < M.Height(); i++)
//   {
//     TRES sum = 0;
//     for (size_t j = 0; j < M.Width(); j++)
//     {
//       sum += M(i, j) * v(j);
//     }
//     prod(i) = sum;
//   }
//   return prod;
// }

// template <typename T, ORDERING ORD>
// size_t Matrix<T, ORD>::max_in_Col_index(size_t j) const
// {
//   if (j >= width)
//   {
//     throw std::invalid_argument("Error: Column Index must not be greater than
//     Matrix width");
//   }
//   size_t max_index = 0;
//   for (size_t i = 1; i < height; i++)
//   {
//     if ((*this)(i, j) > (*this)(max_index, j))
//     {
//       max_index = i;
//     }
//   }
//   return max_index;
// }

// template <typename T, ORDERING ORD>
// std::pair<Matrix<T, ORD>, Matrix<T, ORD>> Matrix<T, ORD>::LU()
// {
//   if (height != width)
//   {
//     throw std::invalid_argument("Error: Matrix must by quadratic! \n");
//   }

//   size_t n = height;
//   Matrix<T, ORD> A(*this); // get Matrix itself
//   Matrix<T, ORD> U(n, n);
//   Matrix<T, ORD> L(n, n);

//   for (size_t k = 0; k < height; k++)
//   {
//     T temp = 0;
//     for (size_t j = k; j < height; k++)
//     {
//       for (size_t l = 0; k - 1; l++)
//       {
//         temp += L(k, l) * U(l, j);
//       }
//       U(k, j) = A(k, j) - temp;
//     }
//     if (U(k, k) == 0)
//     {
//       throw std::invalid_argument("Error: Matrix is singular");
//     }
//     else
//       for (size_t j = k; k < height; j++)
//       {
//         T temp = 0;
//         for (size_t l = 0; l < k - 1; l++)
//         {
//           temp += L(j, l) * U(l, k);
//         }
//         U(j, k) = 1 / U(k, k) * (A(j, k) - temp);
//       }
//   }

//   return std::pair(L, U);
// }

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

