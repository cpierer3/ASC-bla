#ifndef FILE_EXPRESSION
#define FILE_EXPRESSION

#include <cassert>
#include "matrixexpr.hpp"
namespace ASC_bla
{

  template <typename T>
  class VecExpr
  {
  public:
    auto derived() const { return static_cast<const T&> (*this); }//casts the VecExpr into the class that we are initiated with
    size_t size() const { return derived().size(); }
    auto operator() (size_t i) const { return derived()(i); }//supposedly casts parts of the vector
  };

  // ***************** Sum of two vectors *****************

  template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    SumVecExpr (TA _a, TB _b) : a(_a), b(_b) { }
    auto operator() (size_t i) const { return a(i)+b(i); }
    size_t size() const { return a.size(); }
  };

  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    assert (a.size() == b.size());
    return SumVecExpr(a.derived(), b.derived());
  }

  // ***************** Difference of two vectors *****************

  template <typename TA, typename TB>
  class DifVecExpr : public VecExpr<DifVecExpr<TA,TB>>
  {
    TA a;
    TB b;
  public:
    DifVecExpr (TA _a, TB _b) : a(_a), b(_b) { }
    auto operator() (size_t i) const { return a(i)-b(i); }
    size_t size() const { return a.size(); }      
  };
  
  template <typename TA, typename TB>
  auto operator- (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    assert (a.size() == b.size());
    return DifVecExpr(a.derived(), b.derived());
  }


// ***************** Scaling a vector *****************
  
  template <typename TSCAL, typename TV>
  class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL,TV>>
  {
    TSCAL scal;
    TV vec;
  public:
    ScaleVecExpr (TSCAL _scal, TV _vec) : scal(_scal), vec(_vec) { }
    auto operator() (size_t i) const { return scal*vec(i); }
    size_t size() const { return vec.size(); }      
  };
  
  template <typename T>
  auto operator* (double scal, const VecExpr<T> & v)
  {
    return ScaleVecExpr(scal, v.derived());
  }

 
  template <typename TA, typename TB>
  auto dot (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    assert (a.size() == b.size());

    using elemtypeA = typename std::invoke_result<TA,size_t>::type; //elementypeA is return type of TA being called with argument size t, so TA(t)
    using elemtypeB = typename std::invoke_result<TB,size_t>::type;
    using TSUM = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());

    TSUM sum = 0;
    for (size_t i = 0; i < a.size(); i++)
      sum += a(i)*b(i);
    return sum;
  }

  // ***************** Output operator *****************

  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const VecExpr<T> & v)
  {
    if (v.size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.size(); i++)
      ost << ", " << v(i);
    return ost;
  }


  // ***************** Matrix times vector *****************
  template<typename TA, typename TB>
  class MatTimesVecExpr : public VecExpr<MatTimesVecExpr <TA, TB>> {
    TA a;
    TB b;
  public:
    MatTimesVecExpr (TA _a, TB _b) : a(_a), b(_b) {}

    auto operator()(size_t i) const {
      return dot(a.row(i), b);
    }

    // provide size so VecExpr::size() works
    size_t size() const { return a.rows(); }

  };

  // accept const matrix expressions
  template<typename TA, typename TB>
  auto operator*(const MatrixExpr<TA> &a, const VecExpr<TB> &b) {
    assert (a.cols() == b.size());
    return MatTimesVecExpr(a.derived(), b.derived());
  }
}
 
#endif
