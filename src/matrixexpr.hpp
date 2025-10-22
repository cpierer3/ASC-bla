#ifndef FILE_EXPRESSION
#define FILE_EXPRESSION

#include <cassert>

namespace ASC_bla{

template <typename T>
class MatrixExpr 
{
public:
    auto derived() const { return static_cast<const T&> (*this);}
    size_t height() const { return derived().height();}
    size_t width() const { return derived().width();}
    auto operator () (size_t i, size_t j) const { return derived()(i, j)}   

};

// ************sum of two matrices
template <typename TA, typename TB>
class SumMatrixExpr: public MatrixExpr<SumMatrixExpr<TA, TB>>
{
    TA A;
    TB B;
public:
    SumMatExpr (TA _A, TB _B) : A(_A), B(_B) {}
    auto operator() (size_t i, size_t j) const { return A(i,j)+B(i,j); } //??
    size_t height() const { return A.height();}
    size_t width() const { return A.width();}  
};

template <typename TA, typename TB>
auto operator+ (const MatrixExpr<TA> & a, const MatrixExpr <TB> &b)
{
    assert (A.height() == B.height() && A.width() == B.width());
    return  SumMatrixExpr(A.derived(), B.derived());
}

// *************scaling a matrix

template <typename TSCAL, typename TM>
class ScaleMatrixExpr : public MatrixExpr<ScaleMatrixExpr<TSCAL, TM>>
{
    TSCAL scal;
    TM mat;

public:
    ScaleMatrixExpr (TSCAL _scal, TM _mat) : scal(_scal), mat(_mat) {}
    auto operator() (size_t i, size_t j) const {return scal*mat(i,j);}
    size_t height() const { return mat.height();}
    size_t width() const { return mat.width();}
    };

template <typename TA>
auto operator* (double scal, const MatrixExpr<TA> &A)
{
    return ScaleMatrixExpr(scal, A.derived());
}



template <typename TA, typename TB>
class ProdMatrixExpr : public MatrixExpr<ProdMatrixExpr<TA, TB>>
{
TA A;
TB B;
// macht es sinn da schöberls vector dot verwenden= uch glaube iwie nicht
public:
    ProdMatrixExpr (TA _A, TB _B) : A(_A), B(_B) {}
    auto operator() (site_t i, size_t j) const
    {
        using elemtypeA = typename std::invoke_result<TA,size_t,size_t>::type;
        using elemtypeB = typename std::invoke_result<TB,size_t, size_t>::type;
        using TSUM = decltype(std::declval<elemtypeA>()*std::declval<elemtypeB>());//glaub nicht 
        TSUM sum = 0;
        for (size_t k = 0; k< A.height; k++){
            sum += A(i,k)*B(k,j);
        }
    }
    size_t height() const { return A.height();}
    size_t width() const { return A.width();}
};
template <typename TA, typename TB>
auto operator* (const MatrixExpr<TA> &A, const MatrixExpr<TB> &B)
{
    assert (A.width()== B.height());
    return ProdMatrixExpr(A.derived(), B.derived());
}
}

#endif