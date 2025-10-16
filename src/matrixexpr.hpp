#ifndef FILE_MATRIX_EXPRESSION
#define FILE_MATRIX_EXPRESSION

#include <cassert>

namespace ASC_bla {

  template<typename T>
  class MatrixExpr {
  public:
    auto derived() const { return static_cast<const T &> (*this); }

    size_t rows() const { return derived().rows(); }

    size_t cols() const { return derived().cols(); }

    auto operator()(size_t i, size_t j) const { return derived()(i, j); }
  };

  // ***************** Sum of two vectors *****************

  template<typename TA, typename TB>
  class SumMatrixExpr : public MatrixExpr<SumMatrixExpr<TA, TB>> {
    TA a;
    TB b;
  public:
    SumMatrixExpr(TA _a, TB _b) : a(_a), b(_b) {}

    auto operator()(size_t i, size_t j) const { return a(i, j) + b(i, j); }

    size_t rows() const { return a.rows(); }

    size_t cols() const { return a.cols(); }
  };

  template<typename TA, typename TB>
  auto operator+(const MatrixExpr<TA> &a, const MatrixExpr<TB> &b) {
    assert (a.rows() == b.rows());
    assert (a.cols() == b.cols());
    return SumMatrixExpr(a.derived(), b.derived());
  }

  // ***************** Scalar multiplication of a matrix *****************
  template<typename TA, typename TS>
  class ScalMultMatrixExpr : public MatrixExpr<ScalMultMatrixExpr<TA, TS>> {
    TA a;
    TS s;
  public:
    ScalMultMatrixExpr(TA _a, TS _s) : a(_a), s(_s) {}
    auto operator()(size_t i, size_t j) const { return s * a(i, j); }
    size_t rows() const { return a.rows(); }
    size_t cols() const { return a.cols(); }
  };

  template<typename TA, typename TS>
  auto operator*(const MatrixExpr<TA> &a, TS s) {
    return ScalMultMatrixExpr(a.derived(), s);
  }

  template<typename TA, typename TS>
  auto operator*(TS s, const MatrixExpr<TA> &a) {
    return a * s;
  }

  // ***************** Matrix multiplication *****************
  template<typename TA, typename TB>
  class MatMultMatrixExpr : public MatrixExpr<MatMultMatrixExpr<TA, TB>> {
    TA a;
    TB b;
  public:
    MatMultMatrixExpr(TA _a, TB _b) : a(_a), b(_b) {}

    auto operator()(size_t i, size_t j) const {
      decltype(a(0, 0) * b(0, 0)) sum = 0;
      for (size_t k = 0; k < a.cols(); k++) {
        sum += a(i, k) * b(k, j);
      }
      return sum;
    }

    size_t rows() const { return a.rows(); }
    size_t cols() const { return b.cols(); }
  };

  template<typename TA, typename TB>
  auto operator*(const MatrixExpr<TA> &a, const MatrixExpr<TB> &b) {
    assert (a.cols() == b.rows());
    return MatMultMatrixExpr(a.derived(), b.derived());
  }

}

#endif