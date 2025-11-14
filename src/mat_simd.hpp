#ifndef FILE_MATSIMD
#define FILE_MATSIMD

#include <matrixview.hpp>
#include "../myASC-HPC/src/simd.hpp"

namespace ASC_bla {
  template<typename T>
  Matrix<T, ORDERING::RowMajor> Multi(const MatrixView<T, ORDERING::RowMajor> A, const MatrixView<T, ORDERING::RowMajor> B) {
    Matrix<T, ORDERING::RowMajor> C(A.rows(), B.cols());

    for (size_t i = 0; i < A.rows() / 2; i += 2) {
      auto r1 = A.row(i);
      auto r2 = A.row(i + 1);
      for (size_t j = 0; j < B.cols() / 2; j += 2) {
        auto c1 = B.col(j);
        auto c2 = B.col(j+1);

        auto a00 = 0.0, a01=0.0, a10=0.0, a11 = 0.0;
        for (size_t k=0; k<A.cols(); k++) {
          a00 += r1(k) * c1(k);
          a01 += r1(k) * c2(k);
          a10 += r2(k) * c1(k);
          a11 += r2(k) * c2(k);
        }

        C(i, j) = a00;
        C(i, j + 1) = a01;
        C(i + 1, j) = a10;
        C(i + 1, j + 1) = a11;
      }
    }
  }
}

#endif