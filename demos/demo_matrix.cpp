#include <iostream>
#include <matrix.hpp>

namespace bla = ASC_bla;

int main()
{
  size_t n = 3;
  size_t m = 3;
  bla::Matrix<double, ASC_bla::RowMajor> A(n, m);
  bla::Matrix<double, ASC_bla::ColMajor> B(n, m);

  for (size_t i = 0; i < A.Width(); i++){
    for (size_t j=0; j< A.Height(); j++ )
        {
        A(i, j) = j+1;
        B(i, j) = (i+1)*(j+1);
        }
    }
  auto C = B-A;

  
  std::cout << "x+y = " << "\n" << C << std::endl;
}

