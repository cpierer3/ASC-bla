#ifndef FILE_MATRIX
#define FILE_MATRIX

#include <iostream>
#include <vector.hpp>
#include <cmath>

namespace ASC_bla {

enum ORDERING { ColMajor, RowMajor };
template <typename T, ORDERING ORD>
class Matrix {
    size_t height;
    size_t width;
    T * data;

public:
    Matrix (size_t _height, size_t _width) 
        : height(_height), width(_width), data(new T[_height*_width]) { ; }
    

    Matrix (const Matrix & m)
      : Matrix(m.Height(), m.Width())
    {
      *this = m;
    }
    
    std::pair<Matrix<T, ORD>, Matrix<T, ORD>> LU();
    size_t max_in_Col_index(size_t j) const;

    Matrix & operator=(const Matrix & M2)
    {
      for (size_t i = 0; i < height; i++){
        for (size_t j = 0; j< width; j++){
        if constexpr (ORD == RowMajor){
          data[i*width+j] = M2(i,j);
        }
        if constexpr (ORD == ColMajor){
          data[i+height*j] = M2(i,j);
        }
        }
    }
    return *this;
    }

    Matrix & operator= (Matrix && m2)
    {
      std::swap(height, m2.height);
      std::swap(width, m2.width);
      std::swap(data, m2.data);
      return *this;
    }

~Matrix () { delete [] data; }


size_t Height() const { return height; }
size_t Width() const { return width;}
T & operator()(size_t i, size_t j) {
  if constexpr (ORD == RowMajor){
    return data[i*width+j];
  }
  if constexpr (ORD == ColMajor){
    return data[i+height*j];
  }
}
    
const T & operator() (size_t i, size_t j) const {
    if constexpr (ORD == RowMajor){
    return data[i*width+j];
  }
  if constexpr (ORD == ColMajor){
    return data[i+height*j];
  }
}
};

  template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
  auto operator+ (const Matrix<T1, ORD1> & M1, const Matrix<T2, ORD2> & M2)
    -> Matrix <decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
  {
    if (M1.Height() != M2.Height() || M1.Width() != M2.Width()){
      throw std::invalid_argument("Error: Matrix dimensions do not match!\n"); 
  }

    typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;        //common type casting
    Matrix<TRES, ORD1> sum(M1.Height(), M1.Width());             // use ORDERING of first matrix
    for (size_t i = 0; i < M1.Height(); i++){
        for (size_t j =0; j < M1.Width(); j++){
            sum(i, j) = static_cast<TRES>(M1(i,j)) + static_cast<TRES>(M2(i,j));
            }
    }
    return sum;
  }


  template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
  auto operator- (const Matrix<T1, ORD1> & M1, const Matrix<T2, ORD2> & M2)
    -> Matrix <decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
  {
    if (M1.Height() != M2.Height() || M1.Width() != M2.Width()){
      throw std::invalid_argument("Error: Matrix dimensions do not match!\n"); 
  }

    typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;        //common type casting
    Matrix<TRES, ORD1> sum(M1.Height(), M1.Width());             // use ORDERING of first matrix
    for (size_t i = 0; i < M1.Height(); i++){
        for (size_t j =0; j < M1.Width(); j++){
            sum(i, j) = static_cast<TRES>(M1(i,j)) - static_cast<TRES>(M2(i,j));
            }
    }
    return sum;
  }


template <typename T1, typename T2, ORDERING ORD1, ORDERING ORD2>
auto operator*(const Matrix<T1, ORD1>& M1, const Matrix<T2, ORD2> & M2)
  -> Matrix <decltype(std::declval<T1>() + std::declval<T2>()), ORD1>
{
    if (M1.Width() != M2.Height()){ 
      throw std::invalid_argument("Error: Matrix dimensions do not match!\n"); 
  }
    
    typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;
    Matrix <TRES, ORD1> prod(M1.Height(), M2.Width());
    for (size_t i = 0; i < M1.Height(); ++i) {
        for (size_t j = 0; j < M2.Width(); ++j) {
            TRES sum = 0;         // T gives type of matrix entries
            for (size_t k = 0; k < M1.Width(); ++k) {
                sum += M1(i, k) * M2(k, j);
            }
            prod(i, j) = sum;
        }
    }

    return prod;  
}


template <typename T1, typename T2, ORDERING ORD1>
auto operator*(const Matrix<T1, ORD1> & M, const Vector <T2> & v)
  -> Vector<decltype(std::declval<T1>() + std::declval<T2>())>
{
  if (M.Width() != v.Size()){
    throw std::invalid_argument("Error: Matrix Width must match vector size! \n");
  }

  typedef decltype(std::declval<T1>() + std::declval<T2>()) TRES;
  Vector<TRES> prod(M.Height());
  for (size_t i = 0; i < M.Height(); i++){
      TRES sum = 0;
    for (size_t j = 0; j< M.Width(); j++){
      sum += M(i,j)* v(j);
    }
    prod(i) = sum;
  }
  return prod;
}

template <typename T, ORDERING ORD>
size_t Matrix <T, ORD>::max_in_Col_index(size_t j) const{
 if (j >= width){
  throw std::invalid_argument("Error: Column Index must not be greater than Matrix width");
 }
 size_t max_index = 0;
 for (size_t i = 1; i< height; i++){
  if((*this)(i,j)> (*this)(max_index,j)){
    max_index = i;
  }
 }
 return max_index;
}

template <typename T, ORDERING ORD>
std::pair<Matrix<T, ORD>, Matrix<T, ORD>> Matrix <T, ORD>::LU() {
  if (height != width) {
    throw std::invalid_argument("Error: Matrix must by quadratic! \n");
  }

  size_t n = height;
  Matrix <T, ORD> A(*this);   //get Matrix itself
  Matrix <T, ORD> U(n,n);
  Matrix <T, ORD> L(n,n);

  for (size_t k = 0; k< height; k++){
    T temp = 0;
    for (size_t j = k; j< height; k++){
      for(size_t l = 0; k-1; l++ ){
      temp += L(k,l)*U(l,j);
      }
      U(k,j) = A(k,j) - temp;
    }
    if (U(k,k) ==0){
      throw std::invalid_argument("Error: Matrix is singular");
    }
    else
      for (size_t j = k; k< height; j ++){
        T temp = 0;
        for (size_t l = 0; l< k-1; l++){
          temp += L(j,l)*U(l,k);
        }
        U(j,k) = 1/U(k,k)*(A(j,k)-temp);
      } 
  }

    return std::pair(L,U);
  }


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
};
#endif