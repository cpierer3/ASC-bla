#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include "matrixexpr.hpp"


namespace ASC_bla
{

enum ORDERING { ColMajor, RowMajor};
template <typename T, ORDERING ORD>
class MatrixView {
protected:
  size_t m_rows, m_cols, m_dist;
  T * m_data;

public:
    MatrixView() = default;
    MatrixView(const MatrixView &) = default;

    template <ORDERING ORD2>
    MatrixView (const MatrixView <T,ORD2> & M2)
    : m_rows(M2.Rows()), m_cols(M2.Cols()), m_dist(M2.dist()) {}

    MatrixView (size_t rows, size_t cols, T* data)
    : m_rows(rows), m_cols(cols), m_data(data) {
        if constexpr (ORD == RowMajor){
            m_dist = m_cols;
        }
        else { m_dist = m_rows; }
    }
    MatrixView (size_t rows, size_t cols, size_t dist, T* data)
    : m_data(data), m_rows(rows), m_cols(cols), m_dist(dist) {}

    template <typename TB>
    MatrixView & operator= (const MatrixExpr <TB> & M2)
    {
        assert (m_rows == M2.rows() && m_cols == M2.cols());
        for (size_t i = 0; i < m_rows; i++){
            for (size_t j = 0; j < m_cols; j++){
                if constexpr (ORD == RowMajor){
                    m_data[i*m_dist+j] = M2(i,j);
                }
                else
                    m_data[i+j*m_dist] = M2(i,j);
            }
        }
        return *this;
    }

    MatrixView & operator= (T scal)
    {
        for (size_t i = 0; i < m_rows; i++){
            for (size_t j = 0; j < m_cols; j++){
                if constexpr (ORD == RowMajor)
                m_data[i*m_dist+j] = scal;
                else
                m_data[i+m_dist*j]= scal;
            }
        }
        return *this;
    }

    T * data() const{ return m_data;}
    size_t cols() const { return m_cols;}
    size_t rows() const { return m_rows;}
    size_t dist() const { return m_dist;}

    T & operator() (size_t i, size_t j) {
        if constexpr (ORD == RowMajor){
            return m_data[i*m_dist+j];
        }
        else
            return m_data[i+m_dist*j];
    }
};
 template <typename T, ORDERING ORD>
 class Matrix : public MatrixView<T, ORD>
 {
    typedef MatrixView<T, ORD> BASE;
    using BASE::m_rows;
    using BASE::m_cols;
    using BASE::m_data;
public:
    Matrix (size_t rows, size_t cols)
        : MatrixView <T, ORD> (rows, cols, new T[rows*cols]) { }

    Matrix (const Matrix & M)                   //delegating consr
        : Matrix(M.rows(), M.cols())
    {
        *this = M;
    }

    Matrix (Matrix && M)
        : MatrixView<T, ORD> (0, nullptr)
    {
        std::swap(m_rows, M.m_rows);
        std::swap(m_cols, M.m_cols);
        std::swap(m_data, M.m_data);
    }

    template <typename TB>
    Matrix (const MatrixExpr <TB> & M)
        : Matrix(M.rows(), M.cols())
    {
        *this = M;
    }

    ~Matrix () { delete [] m_data; }

    using BASE::operator=;
    Matrix & operator= (const Matrix & M2)
    {
        assert (m_rows == M2.m_rows && m_cols == M2.m_cols);
        for (size_t i = 0; i < m_rows; i++){
            for (size_t j = 0; j < m_cols; j++){
                if constexpr (ORD == RowMajor)
                    m_data[i*m_dist+j] = M2(i,j);
            }
        }
    }

    Matrix & operator= (Matrix && M2)
    {
        std::swap(m_rows, M2.m_rows);
        std::swap(m_cols, M2.m_cols);
        std::swap(m_data, M2.m_data);
        return *this
    }

    size_t Rows() const { return m_rows;}
    size_t Cols() const { return m_cols;}

    T & operator() (size_t i, size_t j)
    {
        if constexpr (ORD == RowMajor)
            return m_data[i * m_cols + j];
        else
            return m_data[j * m_rows + i];
    }
    const T& operator() (site_t i, size_tj) const 
    {
        if constexpr (ORD == RowMajor)
            return m_data[i * m_cols + j];
        else
            return m_data[j * m_rows + i];
    }
 };

  template <typename T, ORDERING ORD>
  std::ostream & operator<< (std::ostream & ost, const Matrix<T, ORD> & M)
  {
    for (size_t i = 0; i < M.Height(); i++){
        for (size_t j = 0; j < M.Width(); j++){
            ost << ", " << M(i,j);
        }
        ost << "\n";
    }
      
    return ost;
  }
}



#endif