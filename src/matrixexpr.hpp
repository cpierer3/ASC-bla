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

  // ***************** Sum of two matrices *****************
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

  // ***************** Difference of two matrices *****************
  template<typename TA, typename TB>
  class DifMatrixExpr : public MatrixExpr<DifMatrixExpr<TA, TB>> {
    TA a;
    TB b;
  public:
    DifMatrixExpr(TA _a, TB _b) : a(_a), b(_b) {}
    auto operator()(size_t i, size_t j) const { return a(i, j) - b(i, j); }
    size_t rows() const { return a.rows(); }
    size_t cols() const { return a.cols(); }
  };

  template<typename TA, typename TB>
  auto operator-(const MatrixExpr<TA> &a, const MatrixExpr<TB> &b) {
      assert (a.rows() == b.rows());
      assert (a.cols() == b.cols());
      return DifMatrixExpr(a.derived(), b.derived());
  }


  // ***************** Matrix times Matrix *****************
  template<typename TA, typename TB>
  class MatMultMatrixExpr : public MatrixExpr<MatMultMatrixExpr<TA, TB>> {
    TA a;
    TB b;
  public:
    MatMultMatrixExpr(TA _a, TB _b) : a(_a), b(_b) { }

    auto operator()(size_t i, size_t j) const {
        return dot(a.row(i), b.col(j));
    }

    size_t rows() const { return a.rows(); }
    size_t cols() const { return b.cols(); }
  };

  template<typename TA, typename TB>
  auto operator*(const MatrixExpr<TA> &a, const MatrixExpr<TB> &b) {
      assert (a.cols() == b.rows());
      return MatMultMatrixExpr(a.derived(), b.derived());
  }

  // ***************** Scaling a matrix *****************
  template <typename TSCAL, typename TM>
  class ScaleMatrixExpr : public MatrixExpr<ScaleMatrixExpr<TSCAL,TM>>
  {
    TSCAL scal;
    TM mat;
  public:
    ScaleMatrixExpr(TSCAL _scal, TM _mat) : scal(_scal), mat(_mat) { }
    auto operator() (size_t i, size_t j) const { return scal*mat(i,j); }
    size_t rows() const { return mat.rows(); }
    size_t cols() const { return mat.cols(); }
  };

  template <typename T>
  auto operator* (double scal, const MatrixExpr<T> & v)
  {
      return ScaleMatrixExpr(scal, v.derived());
  }

//  // ***************** Scalar multiplication of a matrix *****************
//  template<typename TA, typename TS>
//  class ScalMultMatrixExpr : public MatrixExpr<ScalMultMatrixExpr<TA, TS>> {
//    TA s;
//    MatrixExpr<TS> m;
//  public:
//    ScalMultMatrixExpr(TA _s, MatrixExpr<TS> _m) : m(_m), s(_s) {}
//    auto operator()(size_t i, size_t j) const { return m(i, j) * s; }
//    size_t rows() const { return m.rows(); }
//    size_t cols() const { return m.cols(); }
//  };
//
//  template<typename TA, typename TS>
//  auto operator*(TS s, const MatrixExpr<TA> &m) {
//    return ScalMultMatrixExpr(m.derived(), s);
//  }
}

#endif
